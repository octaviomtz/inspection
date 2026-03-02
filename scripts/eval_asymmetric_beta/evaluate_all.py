#!/usr/bin/env python3
"""
Evaluation script for O+ (Asymmetric beta-Scaling v2).

Compares the new feature against three baselines using:
  1. DockQ v2 (interface quality: DockQ, iRMSD, LRMSD, fnat, F1)
  2. CDR-specific metrics (H3/L3 RMSD, contact ratio)
  3. Boltz confidence metrics (confidence_score, iptm, complex_plddt)

Usage:
    conda activate dockq2
    python evaluate_all.py \
        --predictions_dir /path/to/predictions_examples \
        --ground_truth_dir /path/to/pdb_minimized \
        --cdrs_csv /path/to/cdrs.csv \
        --output_dir /path/to/evaluation_results

The predictions_dir should contain subdirectories:
    antigen_cut/                                           (baseline)
    antigen_cut_vhvl_msa_pocket_ab_boltz_post_2023_omm/   (pocket restraints)
    antigen_cut_contact_restraints/                        (contact restraints)
    new_feature_cdr3_beta/                                 (O+ asymmetric beta)
"""

import argparse
import ast
import csv
import json
import os
import re
import subprocess
import sys
import tempfile
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.Align import PairwiseAligner
from Bio.Data.IUPACData import protein_letters_3to1
from Bio.PDB import PDBParser, MMCIFParser, Superimposer


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

METHOD_BASELINE = "baseline"
METHOD_POCKET = "pocket_restraints"
METHOD_CONTACT = "contact_restraints"
METHOD_ASYM_BETA = "asymmetric_beta"

CONTACT_DISTANCE = 5.0  # Angstroms for contact definition

CAPRI_THRESHOLDS = {
    "incorrect": (0.0, 0.23),
    "acceptable": (0.23, 0.49),
    "medium": (0.49, 0.80),
    "high": (0.80, 1.01),
}


# ---------------------------------------------------------------------------
# Utility: residue name -> one-letter code
# ---------------------------------------------------------------------------

def resname_to_one(resname: str) -> str:
    return protein_letters_3to1.get(resname.capitalize(), "X")


# ---------------------------------------------------------------------------
# Step 1: File discovery
# ---------------------------------------------------------------------------

def extract_complex_name_baseline(folder_name: str) -> str:
    """boltz_results_7TRH_HBG -> 7TRH_HBG"""
    return folder_name.replace("boltz_results_", "")


def extract_complex_name_pocket(folder_name: str) -> str:
    """boltz_results_restraint_to_A_7TRH_HBG_B_W_109 -> 7TRH_HBG"""
    m = re.match(r"boltz_results_restraint_to_\w_(\w+_\w+)_\w_\w+_\d+", folder_name)
    if m:
        return m.group(1)
    return None


def extract_complex_name_contact(folder_name: str) -> str:
    """boltz_results_restraint_7TRH_HBG_hbond_23 -> 7TRH_HBG"""
    m = re.match(r"boltz_results_restraint_(\w+_\w+)_(?:hbond|hydrophobic|salt_bridge)_\d+", folder_name)
    if m:
        return m.group(1)
    return None


def extract_complex_name_asym_beta(folder_name: str) -> str:
    """boltz_results_7TRH_HBG_cdr3_beta -> 7TRH_HBG
    Also handles: boltz_results_7TRH_HBG_asymmetric_beta, boltz_results_7TRH_HBG"""
    name = folder_name.replace("boltz_results_", "")
    # Remove known suffixes
    for suffix in ["_cdr3_beta", "_asymmetric_beta"]:
        if name.endswith(suffix):
            name = name[: -len(suffix)]
    return name


def find_prediction_files(result_dir: Path):
    """Find all structure prediction files (.pdb or .cif) in a result directory.

    Returns list of dicts: {model_idx, pred_path, confidence_path}
    """
    predictions = []
    pred_dir = result_dir / "predictions"
    if not pred_dir.exists():
        return predictions

    for subdir in pred_dir.iterdir():
        if not subdir.is_dir():
            continue
        for f in sorted(subdir.iterdir()):
            m = re.match(r".*_model_(\d+)\.(pdb|cif)$", f.name)
            if m:
                model_idx = int(m.group(1))
                conf_name = f"confidence_{subdir.name}_model_{model_idx}.json"
                conf_path = subdir / conf_name
                predictions.append({
                    "model_idx": model_idx,
                    "pred_path": str(f),
                    "confidence_path": str(conf_path) if conf_path.exists() else None,
                })
    return predictions


def discover_files(predictions_dir: Path, ground_truth_dir: Path):
    """Scan all prediction folders and build a manifest.

    Returns a DataFrame with columns:
        complex, method, variant, model_idx, pred_path, confidence_path, gt_path
    """
    rows = []
    gt_files = {f.stem: str(f) for f in ground_truth_dir.glob("*.pdb")}

    method_configs = [
        (METHOD_BASELINE, "antigen_cut", extract_complex_name_baseline),
        (METHOD_POCKET, "antigen_cut_vhvl_msa_pocket_ab_boltz_post_2023_omm", extract_complex_name_pocket),
        (METHOD_CONTACT, "antigen_cut_contact_restraints", extract_complex_name_contact),
        (METHOD_ASYM_BETA, "new_feature_cdr3_beta", extract_complex_name_asym_beta),
    ]

    for method_name, subdir_name, name_extractor in method_configs:
        method_dir = predictions_dir / subdir_name
        if not method_dir.exists():
            print(f"  Warning: {method_dir} not found, skipping {method_name}")
            continue

        for result_folder in sorted(method_dir.iterdir()):
            if not result_folder.is_dir():
                continue

            complex_name = name_extractor(result_folder.name)
            if complex_name is None:
                continue

            gt_path = gt_files.get(complex_name)
            if gt_path is None:
                continue

            preds = find_prediction_files(result_folder)
            for pred in preds:
                rows.append({
                    "complex": complex_name,
                    "method": method_name,
                    "variant": result_folder.name,
                    "model_idx": pred["model_idx"],
                    "pred_path": pred["pred_path"],
                    "confidence_path": pred["confidence_path"],
                    "gt_path": gt_path,
                })

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Step 2: DockQ evaluation
# ---------------------------------------------------------------------------

def run_dockq(pred_path: str, gt_path: str, json_out: str) -> dict:
    """Run DockQ and return parsed results."""
    cmd = ["DockQ", pred_path, gt_path, "--json", json_out]
    result = subprocess.run(cmd, capture_output=True, text=True)

    if not os.path.exists(json_out):
        print(f"  DockQ failed for {pred_path}: {result.stderr[:200]}")
        return None

    with open(json_out) as f:
        return json.load(f)


def parse_dockq_results(dockq_data: dict, complex_name: str, method: str,
                        variant: str, model_idx: int) -> list:
    """Parse DockQ JSON output into flat rows."""
    rows = []
    if dockq_data is None:
        return rows

    # DockQ v2 JSON has interfaces nested under 'best_result'
    best_result = dockq_data.get("best_result", {})

    for interface_key, interface_data in best_result.items():
        if not isinstance(interface_data, dict):
            continue

        chain1 = interface_data.get("chain1", "")
        chain2 = interface_data.get("chain2", "")
        interface_label = f"{chain1}:{chain2}"

        row = {
            "complex": complex_name,
            "method": method,
            "variant": variant,
            "model_idx": model_idx,
            "interface": interface_label,
            "DockQ": interface_data.get("DockQ"),
            "iRMSD": interface_data.get("iRMSD"),
            "LRMSD": interface_data.get("LRMSD"),
            "fnat": interface_data.get("fnat"),
            "fnonnat": interface_data.get("fnonnat"),
            "F1": interface_data.get("F1"),
            "clashes": interface_data.get("clashes"),
        }
        rows.append(row)

    return rows


def evaluate_dockq(manifest: pd.DataFrame, output_dir: Path) -> pd.DataFrame:
    """Run DockQ for all predictions in the manifest."""
    print("\n=== Step 2: DockQ Evaluation ===")
    dockq_dir = output_dir / "dockq_json"
    dockq_dir.mkdir(parents=True, exist_ok=True)

    all_rows = []
    total = len(manifest)
    for i, (_, row) in enumerate(manifest.iterrows()):
        json_name = f"{row['variant']}_model_{row['model_idx']}.json"
        json_path = str(dockq_dir / json_name)

        if os.path.exists(json_path):
            with open(json_path) as f:
                dockq_data = json.load(f)
        else:
            print(f"  [{i+1}/{total}] DockQ: {row['method']} / {row['complex']} model_{row['model_idx']}")
            dockq_data = run_dockq(row["pred_path"], row["gt_path"], json_path)

        parsed = parse_dockq_results(
            dockq_data, row["complex"], row["method"],
            row["variant"], row["model_idx"],
        )
        all_rows.extend(parsed)

    df = pd.DataFrame(all_rows)
    return df


# ---------------------------------------------------------------------------
# Step 3: CDR-specific metrics
# ---------------------------------------------------------------------------

def load_cdrs(cdrs_csv: str) -> dict:
    """Load CDR definitions from cdrs.csv.

    Returns dict: complex_name -> {
        'cdr3_h': [list of 0-indexed positions],
        'cdr3_l': [list of 0-indexed positions],
        'cdr1_h', 'cdr2_h', 'cdr1_l', 'cdr2_l': similar lists
    }
    """
    cdrs = {}
    with open(cdrs_csv) as f:
        reader = csv.DictReader(f)
        for row in reader:
            name = row["complex"]
            cdrs[name] = {
                "cdr3_h": ast.literal_eval(row["cdr3_h"]),
                "cdr3_l": ast.literal_eval(row["cdr3_l"]),
                "cdr1_h": ast.literal_eval(row["cdr1_h"]),
                "cdr2_h": ast.literal_eval(row["cdr2_h"]),
                "cdr1_l": ast.literal_eval(row["cdr1_l"]),
                "cdr2_l": ast.literal_eval(row["cdr2_l"]),
            }
    return cdrs


def parse_structure(path: str):
    """Parse PDB or CIF file, return BioPython structure."""
    if path.endswith(".cif"):
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)
    return parser.get_structure("s", path)


def get_chain_residues(structure, chain_id: str):
    """Get ordered list of standard residues for a chain."""
    model = structure[0]
    if chain_id not in model:
        return []
    chain = model[chain_id]
    return [r for r in chain if r.id[0] == " "]


def get_chain_sequence(residues):
    """Extract one-letter sequence from residue list."""
    return "".join(resname_to_one(r.resname) for r in residues)


def align_residues(pred_residues, gt_residues):
    """Align two residue lists by sequence, return matched (pred_res, gt_res) pairs.

    Returns a dict mapping 0-indexed pred position -> gt residue object.
    """
    pred_seq = get_chain_sequence(pred_residues)
    gt_seq = get_chain_sequence(gt_residues)

    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -0.5

    alignments = aligner.align(pred_seq, gt_seq)
    if len(alignments) == 0:
        return {}

    aln = alignments[0]

    # Build mapping from pred_pos -> gt_pos using aligned pairs
    mapping = {}
    for pred_idx, gt_idx in aln.indices.T if hasattr(aln, 'indices') else []:
        pass

    # Use the aligned attribute to build position mapping
    pred_aligned, gt_aligned = aln[0], aln[1]

    pred_pos = 0
    gt_pos = 0
    mapping = {}
    for p_char, g_char in zip(str(aln).split("\n")[0] if False else "x", "y"):
        pass

    # More robust: use aln.indices
    try:
        indices = aln.indices  # shape (2, alignment_length)
        for i in range(indices.shape[1]):
            p_idx = indices[0, i]
            g_idx = indices[1, i]
            if p_idx >= 0 and g_idx >= 0:
                mapping[int(p_idx)] = gt_residues[int(g_idx)]
    except AttributeError:
        # Fallback: use aligned property
        aligned = aln.aligned
        for (ps, pe), (gs, ge) in zip(aligned[0], aligned[1]):
            for offset in range(pe - ps):
                mapping[ps + offset] = gt_residues[gs + offset]

    return mapping


def get_ca_coord(residue):
    """Get CA atom coordinates from a residue, or None."""
    if "CA" in residue:
        return residue["CA"].get_vector().get_array()
    return None


def compute_rmsd(coords1, coords2):
    """Compute RMSD between two sets of coordinates (numpy arrays)."""
    if len(coords1) == 0:
        return float("nan")
    diff = np.array(coords1) - np.array(coords2)
    return np.sqrt(np.mean(np.sum(diff ** 2, axis=1)))


def compute_cdr_rmsd(pred_path: str, gt_path: str, cdr_positions: list,
                     chain_id: str) -> float:
    """Compute CDR RMSD after aligning on antibody framework.

    1. Align prediction to GT on non-CDR CA atoms of the antibody chain
    2. Compute RMSD on CDR CA atoms

    Parameters
    ----------
    pred_path : path to prediction structure
    gt_path : path to ground truth structure
    cdr_positions : list of 0-indexed CDR positions in the chain sequence
    chain_id : 'B' for heavy, 'C' for light

    Returns
    -------
    float : CDR RMSD in Angstroms
    """
    pred_struct = parse_structure(pred_path)
    gt_struct = parse_structure(gt_path)

    pred_residues = get_chain_residues(pred_struct, chain_id)
    gt_residues = get_chain_residues(gt_struct, chain_id)

    if not pred_residues or not gt_residues:
        return float("nan")

    # Build alignment mapping: pred_position -> gt_residue
    mapping = align_residues(pred_residues, gt_residues)

    cdr_set = set(cdr_positions)

    # Collect framework (non-CDR) CA atoms for alignment
    framework_pred_atoms = []
    framework_gt_atoms = []
    # Collect CDR CA atoms for RMSD
    cdr_pred_coords = []
    cdr_gt_coords = []

    # Also collect all CDR positions from all CDR regions for framework exclusion
    # We'll get all CDR positions for this chain
    all_cdr_set = cdr_set  # For now, just the target CDR

    for pred_pos in range(len(pred_residues)):
        if pred_pos not in mapping:
            continue

        pred_res = pred_residues[pred_pos]
        gt_res = mapping[pred_pos]

        pred_ca = get_ca_coord(pred_res)
        gt_ca = get_ca_coord(gt_res)
        if pred_ca is None or gt_ca is None:
            continue

        if pred_pos in cdr_set:
            cdr_pred_coords.append(pred_ca)
            cdr_gt_coords.append(gt_ca)
        else:
            framework_pred_atoms.append(pred_res["CA"])
            framework_gt_atoms.append(gt_res["CA"])

    if len(framework_pred_atoms) < 3 or len(cdr_pred_coords) == 0:
        return float("nan")

    # Superimpose on framework
    sup = Superimposer()
    sup.set_atoms(framework_gt_atoms, framework_pred_atoms)
    sup.apply(pred_struct[0].get_atoms())

    # Re-extract CDR CA coords after superposition
    cdr_pred_coords_aligned = []
    cdr_gt_coords_final = []
    for pred_pos in sorted(cdr_set):
        if pred_pos >= len(pred_residues) or pred_pos not in mapping:
            continue
        pred_res = pred_residues[pred_pos]
        gt_res = mapping[pred_pos]
        pred_ca = get_ca_coord(pred_res)
        gt_ca = get_ca_coord(gt_res)
        if pred_ca is not None and gt_ca is not None:
            cdr_pred_coords_aligned.append(pred_ca)
            cdr_gt_coords_final.append(gt_ca)

    return compute_rmsd(cdr_pred_coords_aligned, cdr_gt_coords_final)


def compute_contacts(structure, chain_id_a: str, chain_id_b: str,
                     residue_positions_b: list = None, distance: float = 5.0):
    """Count inter-chain contacts.

    A contact exists when any heavy atom of a residue in chain_a is within
    `distance` of any heavy atom of a residue in chain_b.

    Parameters
    ----------
    structure : BioPython structure
    chain_id_a : antigen chain
    chain_id_b : antibody chain
    residue_positions_b : if given, only count contacts for these 0-indexed
                          positions in chain_b. If None, count all.
    distance : contact distance threshold in Angstroms

    Returns
    -------
    set of (res_id_a, res_pos_b) contact pairs
    """
    model = structure[0]
    if chain_id_a not in model or chain_id_b not in model:
        return set()

    res_a = get_chain_residues(structure, chain_id_a)
    res_b = get_chain_residues(structure, chain_id_b)

    if residue_positions_b is not None:
        pos_set = set(residue_positions_b)
        res_b_filtered = [(i, r) for i, r in enumerate(res_b) if i in pos_set]
    else:
        res_b_filtered = list(enumerate(res_b))

    contacts = set()
    dist_sq = distance ** 2

    for ra in res_a:
        atoms_a = [a for a in ra if a.element != "H"]
        for pos_b, rb in res_b_filtered:
            atoms_b = [a for a in rb if a.element != "H"]
            found = False
            for aa in atoms_a:
                if found:
                    break
                coord_a = aa.get_vector().get_array()
                for ab in atoms_b:
                    diff = coord_a - ab.get_vector().get_array()
                    if np.dot(diff, diff) < dist_sq:
                        contacts.add((ra.id[1], pos_b))
                        found = True
                        break

    return contacts


def compute_contact_metrics(pred_path: str, gt_path: str,
                            cdr3_h_positions: list, cdr3_l_positions: list):
    """Compute H3/L3 contact metrics.

    Returns dict with:
        h3_contacts_pred, h3_contacts_native, h3_native_recovered,
        l3_contacts_pred, l3_contacts_native, l3_native_recovered,
        h3_l3_ratio_pred, h3_l3_ratio_native
    """
    pred_struct = parse_structure(pred_path)
    gt_struct = parse_structure(gt_path)

    # We need to map pred CDR positions to GT CDR positions for native contact recovery
    pred_res_b = get_chain_residues(pred_struct, "B")
    gt_res_b = get_chain_residues(gt_struct, "B")
    pred_res_c = get_chain_residues(pred_struct, "C")
    gt_res_c = get_chain_residues(gt_struct, "C")

    mapping_b = align_residues(pred_res_b, gt_res_b)
    mapping_c = align_residues(pred_res_c, gt_res_c)

    # Map CDR positions from pred numbering to GT numbering
    gt_cdr3_h = [mapping_b[p].id[1] for p in cdr3_h_positions
                 if p in mapping_b]
    gt_cdr3_l = [mapping_c[p].id[1] for p in cdr3_l_positions
                 if p in mapping_c]

    # Compute native contacts (GT)
    native_h3 = compute_contacts(gt_struct, "A", "B", cdr3_h_positions)
    native_l3 = compute_contacts(gt_struct, "A", "C", cdr3_l_positions)

    # Compute predicted contacts
    pred_h3 = compute_contacts(pred_struct, "A", "B", cdr3_h_positions)
    pred_l3 = compute_contacts(pred_struct, "A", "C", cdr3_l_positions)

    h3_pred_count = len(pred_h3)
    l3_pred_count = len(pred_l3)
    h3_native_count = len(native_h3)
    l3_native_count = len(native_l3)

    h3_l3_ratio_pred = h3_pred_count / l3_pred_count if l3_pred_count > 0 else float("nan")
    h3_l3_ratio_native = h3_native_count / l3_native_count if l3_native_count > 0 else float("nan")

    return {
        "h3_contacts_pred": h3_pred_count,
        "l3_contacts_pred": l3_pred_count,
        "h3_contacts_native": h3_native_count,
        "l3_contacts_native": l3_native_count,
        "h3_l3_ratio_pred": h3_l3_ratio_pred,
        "h3_l3_ratio_native": h3_l3_ratio_native,
    }


def evaluate_cdr_metrics(manifest: pd.DataFrame, cdrs: dict) -> pd.DataFrame:
    """Compute CDR-specific metrics for all predictions."""
    print("\n=== Step 3: CDR Metrics ===")
    all_rows = []
    total = len(manifest)

    for i, (_, row) in enumerate(manifest.iterrows()):
        complex_name = row["complex"]
        if complex_name not in cdrs:
            continue

        cdr = cdrs[complex_name]
        print(f"  [{i+1}/{total}] CDR: {row['method']} / {complex_name} model_{row['model_idx']}")

        try:
            h3_rmsd = compute_cdr_rmsd(
                row["pred_path"], row["gt_path"],
                cdr["cdr3_h"], "B",
            )
        except Exception as e:
            print(f"    H3 RMSD failed: {e}")
            h3_rmsd = float("nan")

        try:
            l3_rmsd = compute_cdr_rmsd(
                row["pred_path"], row["gt_path"],
                cdr["cdr3_l"], "C",
            )
        except Exception as e:
            print(f"    L3 RMSD failed: {e}")
            l3_rmsd = float("nan")

        try:
            contact_metrics = compute_contact_metrics(
                row["pred_path"], row["gt_path"],
                cdr["cdr3_h"], cdr["cdr3_l"],
            )
        except Exception as e:
            print(f"    Contact metrics failed: {e}")
            contact_metrics = {
                "h3_contacts_pred": float("nan"),
                "l3_contacts_pred": float("nan"),
                "h3_contacts_native": float("nan"),
                "l3_contacts_native": float("nan"),
                "h3_l3_ratio_pred": float("nan"),
                "h3_l3_ratio_native": float("nan"),
            }

        all_rows.append({
            "complex": complex_name,
            "method": row["method"],
            "variant": row["variant"],
            "model_idx": row["model_idx"],
            "h3_rmsd": h3_rmsd,
            "l3_rmsd": l3_rmsd,
            **contact_metrics,
        })

    return pd.DataFrame(all_rows)


# ---------------------------------------------------------------------------
# Step 4: Confidence metrics
# ---------------------------------------------------------------------------

def evaluate_confidence(manifest: pd.DataFrame) -> pd.DataFrame:
    """Extract Boltz confidence metrics from JSON files."""
    print("\n=== Step 4: Confidence Metrics ===")
    all_rows = []

    for _, row in manifest.iterrows():
        conf_path = row.get("confidence_path")
        if conf_path is None or not os.path.exists(conf_path):
            continue

        with open(conf_path) as f:
            conf = json.load(f)

        all_rows.append({
            "complex": row["complex"],
            "method": row["method"],
            "variant": row["variant"],
            "model_idx": row["model_idx"],
            "confidence_score": conf.get("confidence_score"),
            "iptm": conf.get("iptm"),
            "ptm": conf.get("ptm"),
            "complex_plddt": conf.get("complex_plddt"),
            "complex_iplddt": conf.get("complex_iplddt"),
            "complex_ipde": conf.get("complex_ipde"),
        })

    return pd.DataFrame(all_rows)


# ---------------------------------------------------------------------------
# Step 5: Ensemble diversity
# ---------------------------------------------------------------------------

def compute_ensemble_diversity(manifest: pd.DataFrame, cdrs: dict) -> pd.DataFrame:
    """Compute pairwise CDR-H3 RMSD within each (complex, method, variant) ensemble."""
    print("\n=== Step 5: Ensemble Diversity ===")
    all_rows = []

    groups = manifest.groupby(["complex", "method", "variant"])
    for (complex_name, method, variant), group in groups:
        if complex_name not in cdrs:
            continue
        if len(group) < 2:
            continue

        cdr3_h = cdrs[complex_name]["cdr3_h"]

        # Load all structures
        structures = {}
        for _, row in group.iterrows():
            try:
                s = parse_structure(row["pred_path"])
                res = get_chain_residues(s, "B")
                ca_coords = []
                for pos in sorted(cdr3_h):
                    if pos < len(res):
                        ca = get_ca_coord(res[pos])
                        if ca is not None:
                            ca_coords.append(ca)
                if ca_coords:
                    structures[row["model_idx"]] = np.array(ca_coords)
            except Exception:
                pass

        if len(structures) < 2:
            continue

        # Pairwise RMSD
        pairwise_rmsds = []
        for (i, coords_i), (j, coords_j) in combinations(structures.items(), 2):
            if coords_i.shape == coords_j.shape:
                rmsd = compute_rmsd(coords_i, coords_j)
                pairwise_rmsds.append(rmsd)

        if pairwise_rmsds:
            all_rows.append({
                "complex": complex_name,
                "method": method,
                "variant": variant,
                "mean_pairwise_h3_rmsd": np.mean(pairwise_rmsds),
                "n_models": len(structures),
            })

    return pd.DataFrame(all_rows)


# ---------------------------------------------------------------------------
# Step 6: Aggregation
# ---------------------------------------------------------------------------

def bootstrap_ci(values, n_boot=1000, alpha=0.05):
    """Compute mean and 95% bootstrap CI."""
    values = np.array(values)
    values = values[~np.isnan(values)]
    if len(values) == 0:
        return float("nan"), float("nan"), float("nan")

    boot_means = []
    for _ in range(n_boot):
        sample = np.random.choice(values, size=len(values), replace=True)
        boot_means.append(np.mean(sample))
    boot_means = np.array(boot_means)

    return (
        np.mean(values),
        np.percentile(boot_means, 100 * alpha / 2),
        np.percentile(boot_means, 100 * (1 - alpha / 2)),
    )


def select_best_variant(df: pd.DataFrame, metric: str, higher_is_better=True) -> pd.DataFrame:
    """For methods with multiple variants per complex, select the best variant.

    For each (complex, method), pick the variant whose best model has the
    highest (or lowest) metric value.
    """
    if df.empty:
        return df

    results = []
    for (cmplx, method), group in df.groupby(["complex", "method"]):
        valid = group[metric].dropna()
        if valid.empty:
            results.append(group)
            continue
        if higher_is_better:
            best_idx = valid.idxmax()
        else:
            best_idx = valid.idxmin()
        best_var = group.loc[best_idx, "variant"]
        results.append(group[group["variant"] == best_var])

    if not results:
        return df.iloc[:0]
    return pd.concat(results, ignore_index=True)


def select_top1_model(df: pd.DataFrame, confidence_df: pd.DataFrame) -> pd.DataFrame:
    """Select the model with highest confidence_score per (complex, method, variant)."""
    if confidence_df.empty or df.empty:
        return df.iloc[:0]

    # Find top model per (complex, method, variant)
    top_models = (
        confidence_df
        .sort_values("confidence_score", ascending=False)
        .groupby(["complex", "method", "variant"])
        .first()
        .reset_index()[["complex", "method", "variant", "model_idx"]]
    )

    merged = df.merge(top_models, on=["complex", "method", "variant", "model_idx"])
    return merged


def aggregate_dockq(dockq_df: pd.DataFrame, confidence_df: pd.DataFrame,
                    output_dir: Path):
    """Aggregate DockQ results and produce summary tables."""
    print("\n=== Step 6: Aggregation ===")

    if dockq_df.empty:
        print("  No DockQ results to aggregate.")
        return pd.DataFrame()

    # For restraint methods, select best variant per complex
    dockq_best = select_best_variant(dockq_df, "DockQ", higher_is_better=True)

    summary_rows = []

    for method in dockq_best["method"].unique():
        method_df = dockq_best[dockq_best["method"] == method]

        # --- Overall DockQ (average across interfaces per model, then aggregate) ---
        model_avg = (
            method_df.groupby(["complex", "method", "variant", "model_idx"])["DockQ"]
            .mean()
            .reset_index()
        )

        # Oracle: best model per complex
        oracle_df = model_avg.loc[model_avg.groupby("complex")["DockQ"].idxmax()]
        mean_o, lo_o, hi_o = bootstrap_ci(oracle_df["DockQ"].values)

        # Average: mean across all models per complex, then across complexes
        avg_per_complex = model_avg.groupby("complex")["DockQ"].mean()
        mean_a, lo_a, hi_a = bootstrap_ci(avg_per_complex.values)

        # Top-1
        top1_df = select_top1_model(model_avg, confidence_df)
        if not top1_df.empty:
            mean_t, lo_t, hi_t = bootstrap_ci(top1_df["DockQ"].values)
        else:
            mean_t, lo_t, hi_t = float("nan"), float("nan"), float("nan")

        # CAPRI categories (on oracle)
        n_total = len(oracle_df)
        for cat, (low, high) in CAPRI_THRESHOLDS.items():
            frac = ((oracle_df["DockQ"] >= low) & (oracle_df["DockQ"] < high)).mean()
            summary_rows.append({
                "method": method, "metric": f"CAPRI_{cat}",
                "selection": "oracle", "mean": frac, "ci_low": None, "ci_high": None,
                "n": n_total,
            })

        # DockQ >= 0.23, >= 0.49
        frac_023 = (oracle_df["DockQ"] >= 0.23).mean()
        frac_049 = (oracle_df["DockQ"] >= 0.49).mean()

        for sel, (m, lo, hi) in [("oracle", (mean_o, lo_o, hi_o)),
                                  ("average", (mean_a, lo_a, hi_a)),
                                  ("top1", (mean_t, lo_t, hi_t))]:
            summary_rows.append({
                "method": method, "metric": "DockQ_overall",
                "selection": sel, "mean": m, "ci_low": lo, "ci_high": hi,
                "n": n_total,
            })

        summary_rows.append({
            "method": method, "metric": "DockQ>=0.23",
            "selection": "oracle", "mean": frac_023, "ci_low": None, "ci_high": None,
            "n": n_total,
        })
        summary_rows.append({
            "method": method, "metric": "DockQ>=0.49",
            "selection": "oracle", "mean": frac_049, "ci_low": None, "ci_high": None,
            "n": n_total,
        })

        # --- Per-interface breakdown ---
        for iface in ["A:B", "A:C", "B:C"]:
            iface_df = method_df[method_df["interface"] == iface]
            if iface_df.empty:
                continue
            oracle_iface = iface_df.loc[
                iface_df.groupby("complex")["DockQ"].idxmax()
            ]
            mean_i, lo_i, hi_i = bootstrap_ci(oracle_iface["DockQ"].values)
            summary_rows.append({
                "method": method, "metric": f"DockQ_{iface}",
                "selection": "oracle", "mean": mean_i, "ci_low": lo_i, "ci_high": hi_i,
                "n": len(oracle_iface),
            })

            # iRMSD, LRMSD, fnat for this interface (oracle)
            for submetric in ["iRMSD", "LRMSD", "fnat", "F1"]:
                vals = oracle_iface[submetric].dropna().values
                m_s, lo_s, hi_s = bootstrap_ci(vals)
                summary_rows.append({
                    "method": method, "metric": f"{submetric}_{iface}",
                    "selection": "oracle", "mean": m_s, "ci_low": lo_s, "ci_high": hi_s,
                    "n": len(vals),
                })

    return pd.DataFrame(summary_rows)


def aggregate_cdr_metrics(cdr_df: pd.DataFrame, output_dir: Path) -> pd.DataFrame:
    """Aggregate CDR metrics."""
    if cdr_df.empty:
        return pd.DataFrame()

    cdr_best = select_best_variant(cdr_df, "h3_rmsd", higher_is_better=False)

    summary_rows = []
    for method in cdr_best["method"].unique():
        method_df = cdr_best[cdr_best["method"] == method]

        # Oracle: best H3 RMSD per complex
        oracle = method_df.loc[method_df.groupby("complex")["h3_rmsd"].idxmin()]

        for metric in ["h3_rmsd", "l3_rmsd", "h3_l3_ratio_pred",
                        "h3_contacts_pred", "l3_contacts_pred"]:
            vals = oracle[metric].dropna().values
            m, lo, hi = bootstrap_ci(vals)
            summary_rows.append({
                "method": method, "metric": metric,
                "selection": "oracle", "mean": m, "ci_low": lo, "ci_high": hi,
                "n": len(vals),
            })

    return pd.DataFrame(summary_rows)


# ---------------------------------------------------------------------------
# Step 7: Plotting
# ---------------------------------------------------------------------------

def generate_plots(dockq_df: pd.DataFrame, cdr_df: pd.DataFrame,
                   confidence_df: pd.DataFrame, summary_df: pd.DataFrame,
                   output_dir: Path):
    """Generate comparison plots."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("  matplotlib not available, skipping plots.")
        return

    plots_dir = output_dir / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)

    if dockq_df.empty:
        print("  No data for plots.")
        return

    # --- Plot 1: DockQ comparison bar plot ---
    dockq_summary = summary_df[
        (summary_df["metric"] == "DockQ_overall") &
        (summary_df["selection"] == "oracle")
    ]
    if not dockq_summary.empty:
        fig, ax = plt.subplots(figsize=(8, 5))
        methods = dockq_summary["method"].values
        means = dockq_summary["mean"].values
        ci_lo = dockq_summary["ci_low"].values.astype(float)
        ci_hi = dockq_summary["ci_high"].values.astype(float)
        yerr = np.array([means - ci_lo, ci_hi - means])
        yerr = np.clip(yerr, 0, None)

        colors = ["#55C2FF", "#86E935", "#FC8AD9", "#FFB55A"]
        bars = ax.bar(range(len(methods)), means, yerr=yerr, capsize=5,
                      color=colors[:len(methods)], edgecolor="black", linewidth=0.5)
        ax.set_xticks(range(len(methods)))
        ax.set_xticklabels(methods, rotation=15, ha="right")
        ax.set_ylabel("DockQ (oracle)")
        ax.set_title("Overall DockQ Comparison")
        ax.set_ylim(0, 1)
        plt.tight_layout()
        plt.savefig(plots_dir / "dockq_comparison.pdf")
        plt.close()
        print(f"  Saved dockq_comparison.pdf")

    # --- Plot 2: Per-interface DockQ ---
    iface_summary = summary_df[
        summary_df["metric"].str.startswith("DockQ_A:") |
        summary_df["metric"].str.startswith("DockQ_B:")
    ]
    if not iface_summary.empty:
        fig, axes = plt.subplots(1, 3, figsize=(14, 5), sharey=True)
        for idx, iface in enumerate(["A:B", "A:C", "B:C"]):
            ax = axes[idx]
            idf = iface_summary[iface_summary["metric"] == f"DockQ_{iface}"]
            if idf.empty:
                continue
            methods = idf["method"].values
            means = idf["mean"].values
            ci_lo = idf["ci_low"].values.astype(float)
            ci_hi = idf["ci_high"].values.astype(float)
            yerr = np.clip(np.array([means - ci_lo, ci_hi - means]), 0, None)
            ax.bar(range(len(methods)), means, yerr=yerr, capsize=4,
                   color=colors[:len(methods)], edgecolor="black", linewidth=0.5)
            ax.set_xticks(range(len(methods)))
            ax.set_xticklabels(methods, rotation=20, ha="right", fontsize=8)
            ax.set_title(f"Interface {iface}")
            ax.set_ylim(0, 1)
        axes[0].set_ylabel("DockQ (oracle)")
        plt.suptitle("Per-Interface DockQ Comparison")
        plt.tight_layout()
        plt.savefig(plots_dir / "interface_breakdown.pdf")
        plt.close()
        print(f"  Saved interface_breakdown.pdf")

    # --- Plot 3: CDR RMSD comparison ---
    if not cdr_df.empty:
        cdr_best = select_best_variant(cdr_df, "h3_rmsd", higher_is_better=False)
        oracle_cdr = cdr_best.loc[cdr_best.groupby(["complex", "method"])["h3_rmsd"].idxmin()]

        fig, axes = plt.subplots(1, 2, figsize=(10, 5))
        for idx, (metric, label) in enumerate([("h3_rmsd", "CDR-H3 RMSD"),
                                                 ("l3_rmsd", "CDR-L3 RMSD")]):
            ax = axes[idx]
            methods_list = sorted(oracle_cdr["method"].unique())
            positions = []
            labels = []
            for j, meth in enumerate(methods_list):
                vals = oracle_cdr[oracle_cdr["method"] == meth][metric].dropna().values
                if len(vals) > 0:
                    positions.append(vals)
                    labels.append(meth)
            if positions:
                bp = ax.boxplot(positions, tick_labels=labels, patch_artist=True)
                for patch, color in zip(bp["boxes"], colors):
                    patch.set_facecolor(color)
                ax.set_ylabel(f"{label} (A)")
                ax.set_title(label)
                plt.setp(ax.get_xticklabels(), rotation=20, ha="right", fontsize=8)
        plt.suptitle("CDR RMSD Comparison (framework-aligned)")
        plt.tight_layout()
        plt.savefig(plots_dir / "cdr_rmsd_comparison.pdf")
        plt.close()
        print(f"  Saved cdr_rmsd_comparison.pdf")

    # --- Plot 4: H3/L3 contact ratio ---
    if not cdr_df.empty:
        fig, ax = plt.subplots(figsize=(8, 5))
        methods_list = sorted(oracle_cdr["method"].unique())
        positions_ratio = []
        labels_ratio = []
        for meth in methods_list:
            vals = oracle_cdr[oracle_cdr["method"] == meth]["h3_l3_ratio_pred"].dropna().values
            if len(vals) > 0:
                positions_ratio.append(vals)
                labels_ratio.append(meth)
        if positions_ratio:
            bp = ax.boxplot(positions_ratio, tick_labels=labels_ratio, patch_artist=True)
            for patch, color in zip(bp["boxes"], colors):
                patch.set_facecolor(color)
            # Add reference line for biological target
            ax.axhline(y=1.5, color="red", linestyle="--", alpha=0.5, label="Target ratio (1.5)")
            ax.axhline(y=2.0, color="red", linestyle=":", alpha=0.5, label="Target ratio (2.0)")
            ax.set_ylabel("H3/L3 Contact Ratio")
            ax.set_title("H3/L3 Antigen Contact Ratio")
            ax.legend()
            plt.setp(ax.get_xticklabels(), rotation=20, ha="right", fontsize=8)
        plt.tight_layout()
        plt.savefig(plots_dir / "contact_ratio.pdf")
        plt.close()
        print(f"  Saved contact_ratio.pdf")

    # --- Plot 5: Per-complex DockQ scatter (O+ vs baseline) ---
    dockq_best = select_best_variant(dockq_df, "DockQ", higher_is_better=True)
    model_avg = (
        dockq_best.groupby(["complex", "method", "variant", "model_idx"])["DockQ"]
        .mean().reset_index()
    )
    oracle_per_complex = model_avg.loc[
        model_avg.groupby(["complex", "method"])["DockQ"].idxmax()
    ]

    baseline_vals = oracle_per_complex[oracle_per_complex["method"] == METHOD_BASELINE]
    asym_vals = oracle_per_complex[oracle_per_complex["method"] == METHOD_ASYM_BETA]

    if not baseline_vals.empty and not asym_vals.empty:
        merged = baseline_vals[["complex", "DockQ"]].merge(
            asym_vals[["complex", "DockQ"]],
            on="complex", suffixes=("_baseline", "_asym_beta"),
        )
        if not merged.empty:
            fig, ax = plt.subplots(figsize=(6, 6))
            ax.scatter(merged["DockQ_baseline"], merged["DockQ_asym_beta"],
                      alpha=0.7, s=40, edgecolor="black", linewidth=0.5)
            lims = [0, 1]
            ax.plot(lims, lims, "k--", alpha=0.3, label="y=x")
            ax.set_xlabel("Baseline DockQ (oracle)")
            ax.set_ylabel("Asymmetric Beta DockQ (oracle)")
            ax.set_title("Per-Complex DockQ: O+ vs Baseline")
            ax.set_xlim(lims)
            ax.set_ylim(lims)
            ax.legend()
            ax.set_aspect("equal")
            plt.tight_layout()
            plt.savefig(plots_dir / "per_complex_dockq.pdf")
            plt.close()
            print(f"  Saved per_complex_dockq.pdf")


# ---------------------------------------------------------------------------
# Step 8: Statistical tests
# ---------------------------------------------------------------------------

def compute_statistical_tests(dockq_df: pd.DataFrame, output_dir: Path) -> pd.DataFrame:
    """Compute paired Wilcoxon signed-rank tests between O+ and each baseline."""
    from scipy.stats import wilcoxon

    if dockq_df.empty:
        return pd.DataFrame()

    # Get oracle DockQ per complex per method
    dockq_best = select_best_variant(dockq_df, "DockQ", higher_is_better=True)
    model_avg = (
        dockq_best.groupby(["complex", "method", "variant", "model_idx"])["DockQ"]
        .mean().reset_index()
    )
    oracle = model_avg.loc[
        model_avg.groupby(["complex", "method"])["DockQ"].idxmax()
    ][["complex", "method", "DockQ"]]

    asym_data = oracle[oracle["method"] == METHOD_ASYM_BETA].set_index("complex")["DockQ"]

    test_rows = []
    for baseline_method in [METHOD_BASELINE, METHOD_POCKET, METHOD_CONTACT]:
        baseline_data = oracle[oracle["method"] == baseline_method].set_index("complex")["DockQ"]

        common = asym_data.index.intersection(baseline_data.index)
        if len(common) < 5:
            continue

        a = asym_data.loc[common].values
        b = baseline_data.loc[common].values

        stat, p_value = wilcoxon(a, b, alternative="two-sided")
        delta = np.mean(a - b)
        _, ci_lo, ci_hi = bootstrap_ci(a - b)

        test_rows.append({
            "comparison": f"{METHOD_ASYM_BETA} vs {baseline_method}",
            "n_complexes": len(common),
            "mean_delta_DockQ": delta,
            "delta_ci_low": ci_lo,
            "delta_ci_high": ci_hi,
            "wilcoxon_stat": stat,
            "p_value": p_value,
            "p_value_bonferroni": min(p_value * 3, 1.0),
        })

    return pd.DataFrame(test_rows)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Evaluate O+ (Asymmetric beta-Scaling) against baselines",
    )
    parser.add_argument(
        "--predictions_dir", type=Path, required=True,
        help="Directory containing prediction subfolders (antigen_cut, new_feature_cdr3_beta, etc.)",
    )
    parser.add_argument(
        "--ground_truth_dir", type=Path, required=True,
        help="Directory containing ground truth PDB files (pdb_minimized/)",
    )
    parser.add_argument(
        "--cdrs_csv", type=Path, required=True,
        help="Path to cdrs.csv with CDR region definitions",
    )
    parser.add_argument(
        "--output_dir", type=Path, default=Path("evaluation_results"),
        help="Output directory for results and plots",
    )
    parser.add_argument(
        "--skip_dockq", action="store_true",
        help="Skip DockQ computation (use existing results)",
    )
    parser.add_argument(
        "--skip_cdr", action="store_true",
        help="Skip CDR metric computation",
    )
    parser.add_argument(
        "--skip_plots", action="store_true",
        help="Skip plot generation",
    )
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    # ---- Step 1: File discovery ----
    print("=== Step 1: File Discovery ===")
    manifest = discover_files(args.predictions_dir, args.ground_truth_dir)
    manifest.to_csv(args.output_dir / "file_manifest.csv", index=False)
    print(f"  Found {len(manifest)} prediction files across {manifest['method'].nunique()} methods")
    for method in manifest["method"].unique():
        method_df = manifest[manifest["method"] == method]
        n_complexes = method_df["complex"].nunique()
        n_models = len(method_df)
        print(f"    {method}: {n_complexes} complexes, {n_models} models")

    if manifest.empty:
        print("\n  No predictions found. Exiting.")
        return

    # ---- Step 2: DockQ ----
    dockq_path = args.output_dir / "dockq_results.csv"
    if args.skip_dockq and dockq_path.exists():
        print("\n  Loading existing DockQ results...")
        dockq_df = pd.read_csv(dockq_path)
    else:
        dockq_df = evaluate_dockq(manifest, args.output_dir)
        if not dockq_df.empty:
            dockq_df.to_csv(dockq_path, index=False)
            print(f"  Saved {len(dockq_df)} DockQ rows to {dockq_path}")

    # ---- Step 3: CDR metrics ----
    cdr_path = args.output_dir / "cdr_metrics.csv"
    if args.skip_cdr and cdr_path.exists():
        print("\n  Loading existing CDR results...")
        cdr_df = pd.read_csv(cdr_path)
    else:
        cdrs = load_cdrs(str(args.cdrs_csv))
        cdr_df = evaluate_cdr_metrics(manifest, cdrs)
        if not cdr_df.empty:
            cdr_df.to_csv(cdr_path, index=False)
            print(f"  Saved {len(cdr_df)} CDR rows to {cdr_path}")

    # ---- Step 4: Confidence metrics ----
    confidence_df = evaluate_confidence(manifest)
    conf_path = args.output_dir / "confidence_metrics.csv"
    if not confidence_df.empty:
        confidence_df.to_csv(conf_path, index=False)
        print(f"  Saved {len(confidence_df)} confidence rows to {conf_path}")

    # ---- Step 5: Ensemble diversity ----
    cdrs = load_cdrs(str(args.cdrs_csv))
    diversity_df = compute_ensemble_diversity(manifest, cdrs)
    div_path = args.output_dir / "ensemble_diversity.csv"
    if not diversity_df.empty:
        diversity_df.to_csv(div_path, index=False)
        print(f"  Saved {len(diversity_df)} diversity rows to {div_path}")

    # ---- Step 6: Aggregation ----
    summary_dockq = aggregate_dockq(dockq_df, confidence_df, args.output_dir)
    summary_cdr = aggregate_cdr_metrics(cdr_df, args.output_dir)
    summary_df = pd.concat([summary_dockq, summary_cdr], ignore_index=True)
    summary_path = args.output_dir / "summary_table.csv"
    if not summary_df.empty:
        summary_df.to_csv(summary_path, index=False)
        print(f"\n  Saved summary to {summary_path}")

    # Print summary table
    print("\n" + "=" * 70)
    print("SUMMARY RESULTS")
    print("=" * 70)
    if not summary_df.empty:
        key_metrics = summary_df[
            summary_df["metric"].isin([
                "DockQ_overall", "DockQ_A:B", "DockQ_A:C",
                "DockQ>=0.23", "DockQ>=0.49",
                "h3_rmsd", "l3_rmsd", "h3_l3_ratio_pred",
            ])
        ].copy()
        if not key_metrics.empty:
            for _, row in key_metrics.iterrows():
                ci = ""
                if pd.notna(row.get("ci_low")) and pd.notna(row.get("ci_high")):
                    ci = f" [{row['ci_low']:.3f}, {row['ci_high']:.3f}]"
                print(f"  {row['method']:25s} | {row['metric']:20s} | "
                      f"{row['selection']:8s} | {row['mean']:.4f}{ci} (n={row['n']})")

    # ---- Step 7: Statistical tests ----
    stat_tests = compute_statistical_tests(dockq_df, args.output_dir)
    if not stat_tests.empty:
        stat_path = args.output_dir / "statistical_tests.csv"
        stat_tests.to_csv(stat_path, index=False)
        print(f"\n  Statistical tests saved to {stat_path}")
        print("\n  Paired Wilcoxon Tests (DockQ oracle):")
        for _, row in stat_tests.iterrows():
            print(f"    {row['comparison']}: delta={row['mean_delta_DockQ']:.4f} "
                  f"p={row['p_value']:.4f} (Bonf: {row['p_value_bonferroni']:.4f}) "
                  f"n={row['n_complexes']}")

    # ---- Step 8: Plots ----
    if not args.skip_plots:
        print("\n=== Step 7: Generating Plots ===")
        generate_plots(dockq_df, cdr_df, confidence_df, summary_df, args.output_dir)

    print("\n=== Evaluation Complete ===")
    print(f"  All results saved to: {args.output_dir}")


if __name__ == "__main__":
    main()
