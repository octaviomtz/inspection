#!/usr/bin/env python3
"""
Evaluation script for Strategy L+ — CDR3 Beta-Scaling v2.

Compares the L+ v2 feature against 3 baselines across ~48 antibody-antigen
complexes using DockQ v2, CDR3 backbone RMSD (H3 and L3 separately),
CDR3 pLDDT, ensemble diversity, and Boltz confidence metrics.

Usage:
    conda activate dockq2
    python evaluate_cdr3beta_v2.py \
        --predictions_dir predictions_examples \
        --gt_dir pdb_minimized \
        --cdrs_csv examples/cdrs.csv \
        --output_dir evaluation_v2

Optional overrides:
    --feature_folder new_feature_cdr3beta_v2   (subfolder name inside predictions_dir)
    --feature_label  cdr3_beta_scaling_v2       (label used in output tables/plots)

Requirements (in dockq2 environment):
    pip install matplotlib scipy biopython
"""

import argparse
import ast
import csv
import json
import os
import subprocess
import sys
import warnings
from collections import defaultdict
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

try:
    from Bio.PDB import PDBParser, MMCIFParser, Superimposer, PDBIO, Select
except ImportError:
    sys.exit("BioPython is required. Install with: pip install biopython")

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

# Method definitions: (label, folder_name, multi_setting)
# The new feature entry uses defaults that can be overridden via CLI.
BASELINE_METHODS = [
    ("baseline_antigen_cut", "antigen_cut", False),
    ("baseline_contact_restraints", "antigen_cut_contact_restraints", True),
    ("baseline_pocket", "antigen_cut_vhvl_msa_pocket_ab_boltz_post_2023_omm", True),
]

DEFAULT_FEATURE_FOLDER = "new_feature_cdr3beta_v2"
DEFAULT_FEATURE_LABEL = "cdr3_beta_scaling_v2"

BACKBONE_ATOMS = {"N", "CA", "C", "O"}
MAX_MODELS = 5  # Boltz produces up to 5 models per run


# ---------------------------------------------------------------------------
# Utility: CIF → PDB conversion
# ---------------------------------------------------------------------------

def ensure_pdb(filepath: str, tmp_dir: Path) -> str:
    """
    If filepath is a .cif file, convert it to PDB format and return the
    path to the .pdb file. Otherwise return filepath unchanged.
    """
    if not filepath.endswith(".cif"):
        return filepath

    pdb_path = tmp_dir / (Path(filepath).stem + ".pdb")
    if pdb_path.exists():
        return str(pdb_path)

    try:
        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure("s", filepath)
        io = PDBIO()
        io.set_structure(structure)
        io.save(str(pdb_path))
        return str(pdb_path)
    except Exception as e:
        print(f"  [CIF→PDB ERROR] {filepath}: {e}")
        return filepath  # Fall back to original; DockQ may still handle CIF


# ---------------------------------------------------------------------------
# Utility: CDR data loading
# ---------------------------------------------------------------------------

def load_cdr_data(cdrs_csv: str) -> dict:
    """Load CDR indices from CSV. Returns dict[complex_name] -> cdr_info."""
    cdr_data = {}
    with open(cdrs_csv) as f:
        reader = csv.DictReader(f)
        for row in reader:
            name = row["complex"].strip()
            cdr_data[name] = {
                "heavy_seq": row.get("heavy", ""),
                "light_seq": row.get("light", ""),
                "cdr3_h": ast.literal_eval(row["cdr3_h"]),  # 0-indexed
                "cdr3_l": ast.literal_eval(row["cdr3_l"]),  # 0-indexed
                "cdr1_h": ast.literal_eval(row["cdr1_h"]),
                "cdr2_h": ast.literal_eval(row["cdr2_h"]),
                "cdr1_l": ast.literal_eval(row["cdr1_l"]),
                "cdr2_l": ast.literal_eval(row["cdr2_l"]),
            }
    return cdr_data


# ---------------------------------------------------------------------------
# Step 1: Discover all prediction files
# ---------------------------------------------------------------------------

def discover_predictions(
    predictions_dir: str,
    gt_dir: str,
    feature_folder: str,
    feature_label: str,
) -> list[dict]:
    """
    Walk the predictions directory and build a list of evaluation tasks.

    Each task is a dict with:
        method, complex_name, setting, model_idx, pred_file, gt_file,
        confidence_file, plddt_file, pae_file
    """
    methods = BASELINE_METHODS + [(feature_label, feature_folder, False)]

    gt_files = {p.stem: str(p) for p in Path(gt_dir).glob("*.pdb")}
    tasks = []

    for method_label, folder_name, multi_setting in methods:
        method_dir = Path(predictions_dir) / folder_name
        if not method_dir.exists():
            print(f"  [WARN] Method dir not found: {method_dir}")
            continue

        for result_dir in sorted(method_dir.iterdir()):
            if not result_dir.is_dir():
                continue

            # Extract complex name from folder name
            complex_name = _extract_complex_name(result_dir.name, gt_files)
            if complex_name is None:
                print(f"  [WARN] Cannot match complex for: {result_dir.name}")
                continue

            gt_file = gt_files.get(complex_name)
            if gt_file is None:
                print(f"  [WARN] No ground truth for: {complex_name}")
                continue

            # Determine setting name for multi-setting methods
            setting = result_dir.name if multi_setting else "default"

            # Find the predictions subfolder
            preds_dir = result_dir / "predictions"
            if not preds_dir.exists():
                continue

            # The predictions subfolder may contain one directory per run
            for run_subdir in preds_dir.iterdir():
                if not run_subdir.is_dir():
                    continue

                # Find all model files (pdb or cif)
                model_files = sorted(
                    list(run_subdir.glob("*_model_*.pdb"))
                    + list(run_subdir.glob("*_model_*.cif"))
                )

                for mf in model_files:
                    model_idx = _extract_model_idx(mf.name)
                    if model_idx is None:
                        continue

                    stem_base = mf.stem  # e.g. "7TRH_HBG_model_0"
                    conf_file = run_subdir / f"confidence_{stem_base}.json"
                    plddt_file = run_subdir / f"plddt_{stem_base}.npz"
                    pae_file = run_subdir / f"pae_{stem_base}.npz"

                    tasks.append({
                        "method": method_label,
                        "complex_name": complex_name,
                        "setting": setting,
                        "model_idx": model_idx,
                        "pred_file": str(mf),
                        "gt_file": gt_file,
                        "confidence_file": str(conf_file) if conf_file.exists() else None,
                        "plddt_file": str(plddt_file) if plddt_file.exists() else None,
                        "pae_file": str(pae_file) if pae_file.exists() else None,
                    })

    return tasks


def _extract_complex_name(folder_name: str, gt_files: dict) -> str | None:
    """Extract the complex name (e.g. '7TRH_HBG') from a results folder name."""
    for gt_name in gt_files:
        if gt_name in folder_name:
            return gt_name
    return None


def _extract_model_idx(filename: str) -> int | None:
    """Extract model index from filename like '7TRH_HBG_model_0.pdb'."""
    parts = filename.replace(".pdb", "").replace(".cif", "").split("_")
    for i, p in enumerate(parts):
        if p == "model" and i + 1 < len(parts):
            try:
                return int(parts[i + 1])
            except ValueError:
                pass
    return None


# ---------------------------------------------------------------------------
# Step 2: DockQ evaluation
# ---------------------------------------------------------------------------

def run_dockq(pred_file: str, gt_file: str, json_out: str, tmp_dir: Path) -> dict | None:
    """Run DockQ and return parsed JSON results. Converts CIF→PDB if needed."""
    pred_for_dockq = ensure_pdb(pred_file, tmp_dir)

    try:
        cmd = ["DockQ", pred_for_dockq, gt_file, "--json", json_out]
        result = subprocess.run(
            cmd, capture_output=True, text=True, timeout=120,
        )
        if result.returncode != 0:
            print(f"  [DockQ ERROR] {pred_file}: {result.stderr[:200]}")
            return None
        with open(json_out) as f:
            return json.load(f)
    except Exception as e:
        print(f"  [DockQ EXCEPTION] {pred_file}: {e}")
        return None


def parse_dockq_results(dockq_json: dict) -> dict:
    """Extract key metrics from DockQ JSON output."""
    result = {
        "dockq_total": dockq_json.get("best_dockq", np.nan),
    }

    best = dockq_json.get("best_result", {})

    # Extract per-interface metrics for all chain pairs
    for pair_key, pair_label in [("AB", "AB"), ("AC", "AC"), ("BC", "BC")]:
        if pair_key in best:
            iface = best[pair_key]
            result[f"dockq_{pair_label}"] = iface.get("DockQ", np.nan)
            result[f"irmsd_{pair_label}"] = iface.get("iRMSD", np.nan)
            result[f"lrmsd_{pair_label}"] = iface.get("LRMSD", np.nan)
            result[f"fnat_{pair_label}"] = iface.get("fnat", np.nan)
            result[f"f1_{pair_label}"] = iface.get("F1", np.nan)
            result[f"fnonnat_{pair_label}"] = iface.get("fnonnat", np.nan)
        else:
            for metric in ["dockq", "irmsd", "lrmsd", "fnat", "f1", "fnonnat"]:
                result[f"{metric}_{pair_label}"] = np.nan

    # Average DockQ over antigen-antibody interfaces (AB + AC)
    ab = result.get("dockq_AB", np.nan)
    ac = result.get("dockq_AC", np.nan)
    vals = [v for v in [ab, ac] if not np.isnan(v)]
    result["dockq_antigen_ab"] = np.mean(vals) if vals else np.nan

    return result


def capri_class(dockq: float) -> str:
    """Classify DockQ score into CAPRI quality categories."""
    if np.isnan(dockq):
        return "N/A"
    if dockq >= 0.80:
        return "High"
    if dockq >= 0.49:
        return "Medium"
    if dockq >= 0.23:
        return "Acceptable"
    return "Incorrect"


# ---------------------------------------------------------------------------
# Step 3: CDR3 backbone RMSD
# ---------------------------------------------------------------------------

def parse_structure(filepath: str):
    """Parse a PDB or CIF file, return Bio.PDB Structure."""
    filepath = str(filepath)
    if filepath.endswith(".cif"):
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)
    return parser.get_structure("s", filepath)


def get_chain_residues(structure, chain_id: str) -> list:
    """Get list of standard residues for a given chain, ordered by index."""
    model = structure[0]
    if chain_id not in model:
        return []
    chain = model[chain_id]
    return [r for r in chain.get_residues() if r.id[0] == " "]


def compute_cdr3_rmsd(pred_file: str, gt_file: str, cdr_info: dict) -> dict:
    """
    Compute CDR3 backbone RMSD after aligning on antibody framework.

    Strategy:
      1. Identify framework residues (all antibody residues MINUS all CDRs)
      2. Superimpose prediction onto ground truth using framework Cα atoms
      3. Compute CDR3-H and CDR3-L backbone (N/CA/C/O) RMSD after superposition

    Returns separate cdr3_h_rmsd, cdr3_l_rmsd, and combined.
    """
    result = {
        "cdr3_h_rmsd": np.nan,
        "cdr3_l_rmsd": np.nan,
        "cdr3_combined_rmsd": np.nan,
    }

    try:
        pred_struct = parse_structure(pred_file)
        gt_struct = parse_structure(gt_file)
    except Exception as e:
        print(f"  [PARSE ERROR] {pred_file}: {e}")
        return result

    all_cdr_h = set(cdr_info["cdr1_h"] + cdr_info["cdr2_h"] + cdr_info["cdr3_h"])
    all_cdr_l = set(cdr_info["cdr1_l"] + cdr_info["cdr2_l"] + cdr_info["cdr3_l"])
    cdr3_h_set = set(cdr_info["cdr3_h"])
    cdr3_l_set = set(cdr_info["cdr3_l"])

    pred_B = get_chain_residues(pred_struct, "B")
    pred_C = get_chain_residues(pred_struct, "C")
    gt_B = get_chain_residues(gt_struct, "B")
    gt_C = get_chain_residues(gt_struct, "C")

    if not pred_B or not pred_C or not gt_B or not gt_C:
        print(f"  [CHAIN ERROR] Missing chains in {pred_file} or {gt_file}")
        return result

    # Collect framework Cα atoms for superposition
    fw_atoms_pred = []
    fw_atoms_gt = []

    n_heavy = min(len(pred_B), len(gt_B))
    n_light = min(len(pred_C), len(gt_C))

    for i in range(n_heavy):
        if i not in all_cdr_h and "CA" in pred_B[i] and "CA" in gt_B[i]:
            fw_atoms_pred.append(pred_B[i]["CA"])
            fw_atoms_gt.append(gt_B[i]["CA"])

    for i in range(n_light):
        if i not in all_cdr_l and "CA" in pred_C[i] and "CA" in gt_C[i]:
            fw_atoms_pred.append(pred_C[i]["CA"])
            fw_atoms_gt.append(gt_C[i]["CA"])

    if len(fw_atoms_pred) < 10:
        print(f"  [ALIGN ERROR] Too few framework atoms ({len(fw_atoms_pred)}) for {pred_file}")
        return result

    # Superimpose prediction onto ground truth using framework Cα
    sup = Superimposer()
    sup.set_atoms(fw_atoms_gt, fw_atoms_pred)
    sup.apply(list(pred_struct[0].get_atoms()))

    # After superposition: collect CDR3 backbone coords
    cdr3_h_pred_coords = []
    cdr3_h_gt_coords = []
    cdr3_l_pred_coords = []
    cdr3_l_gt_coords = []

    for i in range(n_heavy):
        if i in cdr3_h_set:
            for atom_name in BACKBONE_ATOMS:
                if atom_name in pred_B[i] and atom_name in gt_B[i]:
                    cdr3_h_pred_coords.append(pred_B[i][atom_name].get_vector().get_array())
                    cdr3_h_gt_coords.append(gt_B[i][atom_name].get_vector().get_array())

    for i in range(n_light):
        if i in cdr3_l_set:
            for atom_name in BACKBONE_ATOMS:
                if atom_name in pred_C[i] and atom_name in gt_C[i]:
                    cdr3_l_pred_coords.append(pred_C[i][atom_name].get_vector().get_array())
                    cdr3_l_gt_coords.append(gt_C[i][atom_name].get_vector().get_array())

    if cdr3_h_pred_coords:
        diff = np.array(cdr3_h_pred_coords) - np.array(cdr3_h_gt_coords)
        result["cdr3_h_rmsd"] = float(np.sqrt(np.mean(np.sum(diff ** 2, axis=1))))

    if cdr3_l_pred_coords:
        diff = np.array(cdr3_l_pred_coords) - np.array(cdr3_l_gt_coords)
        result["cdr3_l_rmsd"] = float(np.sqrt(np.mean(np.sum(diff ** 2, axis=1))))

    all_pred = cdr3_h_pred_coords + cdr3_l_pred_coords
    all_gt = cdr3_h_gt_coords + cdr3_l_gt_coords
    if all_pred:
        diff = np.array(all_pred) - np.array(all_gt)
        result["cdr3_combined_rmsd"] = float(np.sqrt(np.mean(np.sum(diff ** 2, axis=1))))

    return result


# ---------------------------------------------------------------------------
# Step 4: CDR3 pLDDT
# ---------------------------------------------------------------------------

def compute_cdr3_plddt(
    pred_file: str, plddt_file: str | None, cdr_info: dict
) -> dict:
    """Extract CDR3-specific pLDDT values from Boltz npz file."""
    result = {
        "cdr3_h_plddt": np.nan,
        "cdr3_l_plddt": np.nan,
        "cdr3_combined_plddt": np.nan,
    }

    if plddt_file is None or not os.path.exists(plddt_file):
        return result

    try:
        plddt_data = np.load(plddt_file)["plddt"]
    except Exception as e:
        print(f"  [PLDDT ERROR] {plddt_file}: {e}")
        return result

    try:
        struct = parse_structure(pred_file)
        chain_a = get_chain_residues(struct, "A")
        chain_b = get_chain_residues(struct, "B")
    except Exception:
        return result

    n_a = len(chain_a)
    n_b = len(chain_b)

    # CDR3-H: residues in heavy chain (chain B), offset by chain A length
    cdr3_h_indices = [n_a + idx for idx in cdr_info["cdr3_h"]]
    # CDR3-L: residues in light chain (chain C), offset by chain A + chain B
    cdr3_l_indices = [n_a + n_b + idx for idx in cdr_info["cdr3_l"]]

    total_tokens = len(plddt_data)
    h_vals = [plddt_data[i] for i in cdr3_h_indices if i < total_tokens]
    l_vals = [plddt_data[i] for i in cdr3_l_indices if i < total_tokens]

    if h_vals:
        result["cdr3_h_plddt"] = float(np.mean(h_vals))
    if l_vals:
        result["cdr3_l_plddt"] = float(np.mean(l_vals))
    if h_vals or l_vals:
        result["cdr3_combined_plddt"] = float(np.mean(h_vals + l_vals))

    return result


# ---------------------------------------------------------------------------
# Step 5: Boltz confidence metrics
# ---------------------------------------------------------------------------

def extract_confidence(confidence_file: str | None) -> dict:
    """Extract Boltz confidence metrics from JSON file."""
    defaults = {
        "confidence_score": np.nan,
        "complex_plddt": np.nan,
        "ptm": np.nan,
        "iptm": np.nan,
        "protein_iptm": np.nan,
    }
    if confidence_file is None or not os.path.exists(confidence_file):
        return defaults

    try:
        with open(confidence_file) as f:
            data = json.load(f)
        return {
            "confidence_score": data.get("confidence_score", np.nan),
            "complex_plddt": data.get("complex_plddt", np.nan),
            "ptm": data.get("ptm", np.nan),
            "iptm": data.get("iptm", np.nan),
            "protein_iptm": data.get("protein_iptm", np.nan),
        }
    except Exception:
        return defaults


# ---------------------------------------------------------------------------
# Step 6: Ensemble diversity
# ---------------------------------------------------------------------------

def compute_ensemble_diversity(pred_files: list[str], cdr_info: dict) -> dict:
    """
    Compute pairwise CDR3 Cα RMSD among models of the same complex.

    Aligns each pair on antibody framework Cα, then measures CDR3 Cα RMSD.
    Reports separate H3 and L3 diversity in addition to combined.
    """
    result = {
        "ensemble_mean_cdr3_rmsd": np.nan,
        "ensemble_max_cdr3_rmsd": np.nan,
        "ensemble_mean_cdr3h_rmsd": np.nan,
        "ensemble_mean_cdr3l_rmsd": np.nan,
        "ensemble_frac_above_1A": np.nan,
        "ensemble_n_models": len(pred_files),
    }

    if len(pred_files) < 2:
        return result

    structures = []
    for pf in pred_files:
        try:
            structures.append((pf, parse_structure(pf)))
        except Exception:
            continue

    if len(structures) < 2:
        return result

    all_cdr_h = set(cdr_info["cdr1_h"] + cdr_info["cdr2_h"] + cdr_info["cdr3_h"])
    all_cdr_l = set(cdr_info["cdr1_l"] + cdr_info["cdr2_l"] + cdr_info["cdr3_l"])
    cdr3_h_set = set(cdr_info["cdr3_h"])
    cdr3_l_set = set(cdr_info["cdr3_l"])

    pairwise_combined = []
    pairwise_h = []
    pairwise_l = []

    for (pf1, s1), (pf2, s2) in combinations(structures, 2):
        combined, h_rmsd, l_rmsd = _pairwise_cdr3_rmsd(
            s1, s2, all_cdr_h, all_cdr_l, cdr3_h_set, cdr3_l_set
        )
        if combined is not None:
            pairwise_combined.append(combined)
        if h_rmsd is not None:
            pairwise_h.append(h_rmsd)
        if l_rmsd is not None:
            pairwise_l.append(l_rmsd)

    if pairwise_combined:
        result["ensemble_mean_cdr3_rmsd"] = float(np.mean(pairwise_combined))
        result["ensemble_max_cdr3_rmsd"] = float(np.max(pairwise_combined))
        result["ensemble_frac_above_1A"] = float(
            np.mean([r >= 1.0 for r in pairwise_combined])
        )
    if pairwise_h:
        result["ensemble_mean_cdr3h_rmsd"] = float(np.mean(pairwise_h))
    if pairwise_l:
        result["ensemble_mean_cdr3l_rmsd"] = float(np.mean(pairwise_l))

    return result


def _pairwise_cdr3_rmsd(
    s1, s2, all_cdr_h, all_cdr_l, cdr3_h_set, cdr3_l_set
) -> tuple:
    """
    Compute CDR3 Cα RMSD between two structures after framework alignment.
    Returns (combined_rmsd, h_rmsd, l_rmsd) — any can be None if atoms missing.
    """
    b1 = get_chain_residues(s1, "B")
    c1 = get_chain_residues(s1, "C")
    b2 = get_chain_residues(s2, "B")
    c2 = get_chain_residues(s2, "C")

    if not b1 or not c1 or not b2 or not c2:
        return None, None, None

    fw_atoms_1 = []
    fw_atoms_2 = []

    n_heavy = min(len(b1), len(b2))
    n_light = min(len(c1), len(c2))

    for i in range(n_heavy):
        if i not in all_cdr_h and "CA" in b1[i] and "CA" in b2[i]:
            fw_atoms_1.append(b1[i]["CA"])
            fw_atoms_2.append(b2[i]["CA"])
    for i in range(n_light):
        if i not in all_cdr_l and "CA" in c1[i] and "CA" in c2[i]:
            fw_atoms_1.append(c1[i]["CA"])
            fw_atoms_2.append(c2[i]["CA"])

    if len(fw_atoms_1) < 10:
        return None, None, None

    sup = Superimposer()
    sup.set_atoms(fw_atoms_1, fw_atoms_2)
    sup.apply(list(s2[0].get_atoms()))

    # Collect CDR3 Cα coords after alignment
    h_coords_1, h_coords_2 = [], []
    l_coords_1, l_coords_2 = [], []

    for i in range(n_heavy):
        if i in cdr3_h_set and "CA" in b1[i] and "CA" in b2[i]:
            h_coords_1.append(b1[i]["CA"].get_vector().get_array())
            h_coords_2.append(b2[i]["CA"].get_vector().get_array())
    for i in range(n_light):
        if i in cdr3_l_set and "CA" in c1[i] and "CA" in c2[i]:
            l_coords_1.append(c1[i]["CA"].get_vector().get_array())
            l_coords_2.append(c2[i]["CA"].get_vector().get_array())

    h_rmsd = None
    l_rmsd = None
    combined_rmsd = None

    if h_coords_1:
        diff = np.array(h_coords_1) - np.array(h_coords_2)
        h_rmsd = float(np.sqrt(np.mean(np.sum(diff ** 2, axis=1))))
    if l_coords_1:
        diff = np.array(l_coords_1) - np.array(l_coords_2)
        l_rmsd = float(np.sqrt(np.mean(np.sum(diff ** 2, axis=1))))

    all_c1 = h_coords_1 + l_coords_1
    all_c2 = h_coords_2 + l_coords_2
    if all_c1:
        diff = np.array(all_c1) - np.array(all_c2)
        combined_rmsd = float(np.sqrt(np.mean(np.sum(diff ** 2, axis=1))))

    return combined_rmsd, h_rmsd, l_rmsd


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

def run_evaluation(args):
    """Run the full evaluation pipeline."""

    predictions_dir = args.predictions_dir
    gt_dir = args.gt_dir
    cdrs_csv = args.cdrs_csv
    output_dir = Path(args.output_dir)
    feature_folder = args.feature_folder
    feature_label = args.feature_label

    results_dir = output_dir / "results"
    plots_dir = output_dir / "plots"
    dockq_json_dir = output_dir / "dockq_json"
    tmp_dir = output_dir / "tmp_pdb"
    results_dir.mkdir(parents=True, exist_ok=True)
    plots_dir.mkdir(parents=True, exist_ok=True)
    dockq_json_dir.mkdir(parents=True, exist_ok=True)
    tmp_dir.mkdir(parents=True, exist_ok=True)

    print(f"Feature folder : {feature_folder}")
    print(f"Feature label  : {feature_label}")

    # Load CDR data
    print("\nLoading CDR data...")
    cdr_data = load_cdr_data(cdrs_csv)
    print(f"  Loaded CDR info for {len(cdr_data)} complexes")

    # Discover predictions
    print("\nDiscovering predictions...")
    tasks = discover_predictions(
        predictions_dir, gt_dir, feature_folder, feature_label
    )
    print(f"  Found {len(tasks)} prediction files")

    if not tasks:
        print("No predictions found. Check your --predictions_dir path.")
        return

    method_counts = defaultdict(int)
    for t in tasks:
        method_counts[t["method"]] += 1
    for m, c in method_counts.items():
        print(f"    {m}: {c} files")

    # -----------------------------------------------------------------------
    # Run evaluations
    # -----------------------------------------------------------------------
    all_results = []

    for i, task in enumerate(tasks):
        complex_name = task["complex_name"]
        method = task["method"]
        setting = task["setting"]
        model_idx = task["model_idx"]

        print(
            f"\n[{i+1}/{len(tasks)}] {method} / {complex_name} / "
            f"model_{model_idx} ({setting})"
        )

        row = {
            "method": method,
            "complex_name": complex_name,
            "setting": setting,
            "model_idx": model_idx,
            "pred_file": task["pred_file"],
        }

        # --- DockQ ---
        safe_setting = setting.replace("/", "_").replace("\\", "_")
        dockq_out = dockq_json_dir / f"{method}_{complex_name}_{safe_setting}_model{model_idx}.json"
        dockq_json_data = run_dockq(
            task["pred_file"], task["gt_file"], str(dockq_out), tmp_dir
        )
        if dockq_json_data:
            row.update(parse_dockq_results(dockq_json_data))
        else:
            row.update({
                "dockq_total": np.nan, "dockq_antigen_ab": np.nan,
                "dockq_AB": np.nan, "dockq_AC": np.nan, "dockq_BC": np.nan,
            })

        # --- CDR3 RMSD ---
        if complex_name in cdr_data:
            row.update(compute_cdr3_rmsd(
                task["pred_file"], task["gt_file"], cdr_data[complex_name]
            ))
        else:
            row.update({
                "cdr3_h_rmsd": np.nan, "cdr3_l_rmsd": np.nan,
                "cdr3_combined_rmsd": np.nan,
            })

        # --- CDR3 pLDDT ---
        if complex_name in cdr_data:
            row.update(compute_cdr3_plddt(
                task["pred_file"], task["plddt_file"], cdr_data[complex_name]
            ))
        else:
            row.update({
                "cdr3_h_plddt": np.nan, "cdr3_l_plddt": np.nan,
                "cdr3_combined_plddt": np.nan,
            })

        # --- Confidence metrics ---
        row.update(extract_confidence(task["confidence_file"]))

        all_results.append(row)

    # Convert to DataFrame
    df = pd.DataFrame(all_results)
    df.to_csv(results_dir / "all_per_model_results.csv", index=False)
    print(f"\nSaved per-model results: {results_dir / 'all_per_model_results.csv'}")

    # -----------------------------------------------------------------------
    # Ensemble diversity (per method × complex × setting)
    # -----------------------------------------------------------------------
    print("\nComputing ensemble diversity...")
    diversity_rows = []
    groups = df.groupby(["method", "complex_name", "setting"])
    for (method, cname, setting), group in groups:
        if cname not in cdr_data:
            continue
        pred_files = sorted(group["pred_file"].tolist())
        div = compute_ensemble_diversity(pred_files, cdr_data[cname])
        div.update({"method": method, "complex_name": cname, "setting": setting})
        diversity_rows.append(div)

    df_div = pd.DataFrame(diversity_rows)
    if not df_div.empty:
        df_div.to_csv(results_dir / "ensemble_diversity.csv", index=False)
        print(f"  Saved: {results_dir / 'ensemble_diversity.csv'}")

    # -----------------------------------------------------------------------
    # Best-of-N selection (by confidence_score)
    # -----------------------------------------------------------------------
    print("\nSelecting best model per complex (by confidence_score)...")
    df_best = (
        df.sort_values("confidence_score", ascending=False)
        .groupby(["method", "complex_name", "setting"])
        .first()
        .reset_index()
    )
    df_best.to_csv(results_dir / "best_model_results.csv", index=False)
    print(f"  Saved: {results_dir / 'best_model_results.csv'}")

    # For multi-setting methods, also pick best setting per complex
    df_best_setting = (
        df_best.sort_values("confidence_score", ascending=False)
        .groupby(["method", "complex_name"])
        .first()
        .reset_index()
    )
    df_best_setting.to_csv(results_dir / "best_setting_results.csv", index=False)

    # -----------------------------------------------------------------------
    # Summary table
    # -----------------------------------------------------------------------
    print("\nGenerating summary table...")
    summary = generate_summary(df_best_setting)
    summary.to_csv(results_dir / "summary_table.csv")
    print(f"  Saved: {results_dir / 'summary_table.csv'}")
    print("\n" + summary.to_string())

    # -----------------------------------------------------------------------
    # Statistical tests
    # -----------------------------------------------------------------------
    print("\nStatistical comparisons (Wilcoxon signed-rank test)...")
    stats_results = run_statistical_tests(df_best_setting, feature_label)
    if stats_results:
        df_stats = pd.DataFrame(stats_results)
        df_stats.to_csv(results_dir / "statistical_tests.csv", index=False)
        print(f"  Saved: {results_dir / 'statistical_tests.csv'}")
        print(df_stats.to_string(index=False))

    # -----------------------------------------------------------------------
    # Plots
    # -----------------------------------------------------------------------
    print("\nGenerating plots...")
    try:
        generate_plots(df_best_setting, df_div, plots_dir, feature_label)
        print(f"  Saved plots to: {plots_dir}")
    except ImportError:
        print("  [WARN] matplotlib not available. Install with: pip install matplotlib")

    # -----------------------------------------------------------------------
    # Final report
    # -----------------------------------------------------------------------
    generate_report(
        df_best_setting, df_div, summary, stats_results, output_dir, feature_label
    )

    print(f"\nEvaluation complete. Results in: {output_dir}")


# ---------------------------------------------------------------------------
# Summary and statistics
# ---------------------------------------------------------------------------

def generate_summary(df: pd.DataFrame) -> pd.DataFrame:
    """Generate aggregate summary table per method."""
    metrics = [
        "dockq_antigen_ab", "dockq_AB", "dockq_AC", "dockq_BC",
        "irmsd_AB", "irmsd_AC", "lrmsd_AB", "lrmsd_AC",
        "fnat_AB", "fnat_AC",
        "cdr3_h_rmsd", "cdr3_l_rmsd", "cdr3_combined_rmsd",
        "cdr3_h_plddt", "cdr3_l_plddt", "cdr3_combined_plddt",
        "confidence_score", "iptm", "ptm", "complex_plddt",
    ]

    rows = []
    for method, group in df.groupby("method"):
        row = {"method": method, "n_complexes": len(group)}
        for m in metrics:
            if m in group.columns:
                vals = group[m].dropna()
                row[f"{m}_mean"] = vals.mean() if len(vals) > 0 else np.nan
                row[f"{m}_std"] = vals.std() if len(vals) > 0 else np.nan
                row[f"{m}_median"] = vals.median() if len(vals) > 0 else np.nan
        rows.append(row)

    return pd.DataFrame(rows).set_index("method")


def run_statistical_tests(df: pd.DataFrame, feature_label: str) -> list[dict]:
    """Run Wilcoxon signed-rank tests comparing new feature vs each baseline."""
    try:
        from scipy.stats import wilcoxon
    except ImportError:
        print("  [WARN] scipy not available. Install with: pip install scipy")
        return []

    metrics_to_test = [
        ("dockq_antigen_ab", "higher"),   # higher is better
        ("dockq_AB", "higher"),
        ("dockq_AC", "higher"),
        ("cdr3_combined_rmsd", "lower"),  # lower is better
        ("cdr3_h_rmsd", "lower"),
        ("cdr3_l_rmsd", "lower"),
        ("iptm", "neutral"),              # monitor for inflation
    ]

    results = []
    df_new = df[df["method"] == feature_label]

    if df_new.empty:
        print(f"  [WARN] No results found for method '{feature_label}'")
        return results

    for baseline_method in df["method"].unique():
        if baseline_method == feature_label:
            continue

        df_base = df[df["method"] == baseline_method]

        # Find common complexes
        common = set(df_new["complex_name"]) & set(df_base["complex_name"])
        if len(common) < 5:
            continue

        for metric, direction in metrics_to_test:
            if metric not in df.columns:
                continue

            new_vals = (
                df_new[df_new["complex_name"].isin(common)]
                .set_index("complex_name")[metric]
                .dropna()
            )
            base_vals = (
                df_base[df_base["complex_name"].isin(common)]
                .set_index("complex_name")[metric]
                .dropna()
            )

            common_idx = new_vals.index.intersection(base_vals.index)
            if len(common_idx) < 5:
                continue

            x_new = new_vals.loc[common_idx].values
            x_base = base_vals.loc[common_idx].values

            try:
                stat, pval = wilcoxon(x_new, x_base)
            except ValueError:
                continue

            if direction == "higher":
                wins = int(np.sum(x_new > x_base))
                losses = int(np.sum(x_new < x_base))
            elif direction == "lower":
                wins = int(np.sum(x_new < x_base))
                losses = int(np.sum(x_new > x_base))
            else:  # neutral: just report
                wins = int(np.sum(x_new > x_base))
                losses = int(np.sum(x_new < x_base))

            results.append({
                "new_feature": feature_label,
                "baseline": baseline_method,
                "metric": metric,
                "direction": direction,
                "n_common": len(common_idx),
                "new_mean": float(np.mean(x_new)),
                "base_mean": float(np.mean(x_base)),
                "delta_mean": float(np.mean(x_new) - np.mean(x_base)),
                "wins": wins,
                "losses": losses,
                "ties": len(common_idx) - wins - losses,
                "wilcoxon_stat": float(stat),
                "p_value": float(pval),
            })

    return results


# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------

def generate_plots(
    df: pd.DataFrame,
    df_div: pd.DataFrame,
    plots_dir: Path,
    feature_label: str,
):
    """Generate evaluation plots."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # All unique methods in data, preserving a consistent order
    method_order = (
        ["baseline_antigen_cut", "baseline_contact_restraints", "baseline_pocket"]
        + [feature_label]
    )
    colors_map = {
        "baseline_antigen_cut": "#1f77b4",
        "baseline_contact_restraints": "#ff7f0e",
        "baseline_pocket": "#2ca02c",
        feature_label: "#d62728",
    }
    methods_present = [m for m in method_order if m in df["method"].values]

    def _make_boxplot(ax, metric, title, ylabel):
        data = [df[df["method"] == m][metric].dropna().values for m in methods_present]
        labels = [m.replace("baseline_", "").replace("_", "\n") for m in methods_present]
        bp = ax.boxplot(data, tick_labels=labels, patch_artist=True)
        for patch, m in zip(bp["boxes"], methods_present):
            patch.set_facecolor(colors_map.get(m, "gray"))
            patch.set_alpha(0.7)
        ax.set_title(title)
        ax.set_ylabel(ylabel)

    # --- DockQ box plots ---
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    for ax, metric, title in zip(
        axes,
        ["dockq_antigen_ab", "dockq_AB", "dockq_AC"],
        ["DockQ (Ab-Ag avg)", "DockQ (A-B)", "DockQ (A-C)"],
    ):
        _make_boxplot(ax, metric, title, "DockQ score")
        ax.axhline(y=0.23, color="gray", linestyle="--", alpha=0.5)
        ax.axhline(y=0.49, color="gray", linestyle="-.", alpha=0.5)
        ax.axhline(y=0.80, color="gray", linestyle=":", alpha=0.5)
    plt.tight_layout()
    plt.savefig(plots_dir / "dockq_comparison_boxplot.png", dpi=150, bbox_inches="tight")
    plt.close()

    # --- CDR3 RMSD box plots (H3, L3, combined) ---
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    for ax, metric, title in zip(
        axes,
        ["cdr3_h_rmsd", "cdr3_l_rmsd", "cdr3_combined_rmsd"],
        ["CDR3-H RMSD (Å)", "CDR3-L RMSD (Å)", "CDR3 Combined RMSD (Å)"],
    ):
        _make_boxplot(ax, metric, title, "RMSD (Å)")
    plt.tight_layout()
    plt.savefig(plots_dir / "cdr3_rmsd_comparison.png", dpi=150, bbox_inches="tight")
    plt.close()

    # --- CDR3 pLDDT comparison ---
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    for ax, metric, title in zip(
        axes,
        ["cdr3_h_plddt", "cdr3_l_plddt", "cdr3_combined_plddt"],
        ["CDR3-H pLDDT", "CDR3-L pLDDT", "CDR3 Combined pLDDT"],
    ):
        _make_boxplot(ax, metric, title, "pLDDT")
    plt.tight_layout()
    plt.savefig(plots_dir / "cdr3_plddt_comparison.png", dpi=150, bbox_inches="tight")
    plt.close()

    # --- Scatter: new feature vs baseline_antigen_cut ---
    baseline = "baseline_antigen_cut"
    if feature_label in df["method"].values and baseline in df["method"].values:
        df_new = df[df["method"] == feature_label].set_index("complex_name")
        df_base = df[df["method"] == baseline].set_index("complex_name")
        common = df_new.index.intersection(df_base.index)

        if len(common) > 0:
            fig, axes = plt.subplots(1, 3, figsize=(15, 5))
            for ax, metric, title in zip(
                axes,
                ["dockq_antigen_ab", "cdr3_combined_rmsd", "iptm"],
                ["DockQ (Ab-Ag)", "CDR3 RMSD", "iPTM"],
            ):
                if metric not in df_new.columns:
                    continue
                x = df_base.loc[common, metric].values
                y = df_new.loc[common, metric].values
                mask = ~(np.isnan(x) | np.isnan(y))
                ax.scatter(x[mask], y[mask], alpha=0.7, edgecolors="black", linewidths=0.5)
                if mask.any():
                    lims = [
                        min(np.nanmin(x[mask]), np.nanmin(y[mask])) * 0.95,
                        max(np.nanmax(x[mask]), np.nanmax(y[mask])) * 1.05,
                    ]
                    ax.plot(lims, lims, "k--", alpha=0.5)
                ax.set_xlabel(f"Baseline (antigen_cut)")
                ax.set_ylabel(f"L+ v2")
                ax.set_title(title)
            plt.tight_layout()
            plt.savefig(
                plots_dir / "scatter_vs_baseline.png", dpi=150, bbox_inches="tight",
            )
            plt.close()

    # --- Ensemble diversity bar chart ---
    if not df_div.empty and "ensemble_mean_cdr3_rmsd" in df_div.columns:
        div_methods = [m for m in method_order if m in df_div["method"].values]
        if div_methods:
            fig, axes = plt.subplots(1, 3, figsize=(18, 5))
            for ax, metric, title in zip(
                axes,
                ["ensemble_mean_cdr3_rmsd", "ensemble_mean_cdr3h_rmsd", "ensemble_mean_cdr3l_rmsd"],
                ["CDR3 Combined Diversity", "CDR3-H Diversity", "CDR3-L Diversity"],
            ):
                if metric not in df_div.columns:
                    continue
                means = []
                stds = []
                for m in div_methods:
                    vals = df_div[df_div["method"] == m][metric].dropna()
                    means.append(vals.mean() if len(vals) > 0 else 0)
                    stds.append(vals.std() if len(vals) > 0 else 0)
                x = range(len(div_methods))
                labels = [m.replace("baseline_", "").replace("_", "\n") for m in div_methods]
                ax.bar(
                    x, means, yerr=stds, capsize=5,
                    color=[colors_map.get(m, "gray") for m in div_methods], alpha=0.7,
                )
                ax.axhline(y=1.0, color="red", linestyle="--", alpha=0.5, label="1.0 Å")
                ax.set_xticks(x)
                ax.set_xticklabels(labels)
                ax.set_ylabel("Mean pairwise CDR3 Cα RMSD (Å)")
                ax.set_title(title)
                ax.legend()
            plt.tight_layout()
            plt.savefig(plots_dir / "ensemble_diversity.png", dpi=150, bbox_inches="tight")
            plt.close()

    # --- Per-complex DockQ delta heatmap ---
    _plot_per_complex_heatmap(df, plots_dir, feature_label)


def _plot_per_complex_heatmap(df: pd.DataFrame, plots_dir: Path, feature_label: str):
    """Plot per-complex ΔDockQ heatmap (new feature minus each baseline)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    metric = "dockq_antigen_ab"

    if feature_label not in df["method"].values or metric not in df.columns:
        return

    df_new = df[df["method"] == feature_label].set_index("complex_name")[metric]
    baselines = [m for m in df["method"].unique() if m != feature_label]

    if not baselines:
        return

    delta_data = {}
    for bl in baselines:
        df_bl = df[df["method"] == bl].set_index("complex_name")[metric]
        common = df_new.index.intersection(df_bl.index)
        delta_data[bl] = {
            c: df_new[c] - df_bl[c]
            for c in common
            if not (np.isnan(df_new.get(c, np.nan)) or np.isnan(df_bl.get(c, np.nan)))
        }

    complexes = sorted(set().union(*(d.keys() for d in delta_data.values())))
    if not complexes:
        return

    matrix = np.full((len(complexes), len(baselines)), np.nan)
    for j, bl in enumerate(baselines):
        for i, c in enumerate(complexes):
            matrix[i, j] = delta_data[bl].get(c, np.nan)

    fig_h = max(4, len(complexes) * 0.3)
    fig, ax = plt.subplots(figsize=(max(6, len(baselines) * 3), fig_h))
    im = ax.imshow(matrix, cmap="RdYlGn", aspect="auto", vmin=-0.3, vmax=0.3)
    ax.set_xticks(range(len(baselines)))
    ax.set_xticklabels(
        [b.replace("baseline_", "") for b in baselines], rotation=45, ha="right"
    )
    ax.set_yticks(range(len(complexes)))
    ax.set_yticklabels(complexes, fontsize=7)
    ax.set_title(f"ΔDockQ ({feature_label} − baseline)")
    plt.colorbar(im, ax=ax, label="ΔDockQ")
    plt.tight_layout()
    plt.savefig(plots_dir / "per_complex_heatmap.png", dpi=150, bbox_inches="tight")
    plt.close()


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------

def generate_report(
    df: pd.DataFrame,
    df_div: pd.DataFrame,
    summary: pd.DataFrame,
    stats_results: list[dict],
    output_dir: Path,
    feature_label: str,
):
    """Generate a markdown evaluation report."""
    lines = [
        "# Evaluation Report: Strategy L+ — CDR3 Beta-Scaling v2",
        "",
        "## 1. Overview",
        "",
        f"- **Total complexes evaluated**: {df['complex_name'].nunique()}",
        f"- **Methods compared**: {df['method'].nunique()}",
        f"- **New feature label**: {feature_label}",
        "",
    ]

    lines += [
        "## 2. Methods",
        "",
        "| Method | # Complexes |",
        "|--------|-------------|",
    ]
    for method, group in df.groupby("method"):
        lines.append(f"| {method} | {group['complex_name'].nunique()} |")
    lines.append("")

    # DockQ summary
    lines += ["## 3. DockQ Results (best model per complex, selected by confidence_score)", ""]
    dockq_cols = [c for c in summary.columns if "dockq" in c and "_mean" in c]
    if dockq_cols:
        header = "| Method | " + " | ".join(c.replace("_mean", "") for c in dockq_cols) + " |"
        sep = "|--------|" + "|".join(["--------"] * len(dockq_cols)) + "|"
        lines += [header, sep]
        for method in summary.index:
            vals = " | ".join(
                f"{summary.loc[method, c]:.3f}" if not np.isnan(summary.loc[method, c]) else "N/A"
                for c in dockq_cols
            )
            lines.append(f"| {method} | {vals} |")
    lines.append("")

    # CDR3 RMSD summary (H3 and L3 separately — key for L+.2 validation)
    lines += ["## 4. CDR3 Backbone RMSD (Å) — framework-aligned", ""]
    rmsd_cols = [c for c in summary.columns if "rmsd" in c and "_mean" in c and "ensemble" not in c]
    if rmsd_cols:
        header = "| Method | " + " | ".join(c.replace("_mean", "") for c in rmsd_cols) + " |"
        sep = "|--------|" + "|".join(["--------"] * len(rmsd_cols)) + "|"
        lines += [header, sep]
        for method in summary.index:
            vals = " | ".join(
                f"{summary.loc[method, c]:.2f}" if not np.isnan(summary.loc[method, c]) else "N/A"
                for c in rmsd_cols
            )
            lines.append(f"| {method} | {vals} |")
    lines.append("")

    # CDR3 pLDDT
    lines += ["## 5. CDR3 pLDDT", ""]
    plddt_cols = [c for c in summary.columns if "plddt" in c and "_mean" in c]
    if plddt_cols:
        header = "| Method | " + " | ".join(c.replace("_mean", "") for c in plddt_cols) + " |"
        sep = "|--------|" + "|".join(["--------"] * len(plddt_cols)) + "|"
        lines += [header, sep]
        for method in summary.index:
            vals = " | ".join(
                f"{summary.loc[method, c]:.3f}" if not np.isnan(summary.loc[method, c]) else "N/A"
                for c in plddt_cols
            )
            lines.append(f"| {method} | {vals} |")
    lines.append("")

    # iptm inflation check (L+.1 motivation)
    lines += ["## 6. Confidence Metrics (iPTM inflation check)", ""]
    conf_cols = [c for c in summary.columns if c in ["iptm_mean", "ptm_mean", "confidence_score_mean"]]
    if conf_cols:
        header = "| Method | " + " | ".join(c.replace("_mean", "") for c in conf_cols) + " |"
        sep = "|--------|" + "|".join(["--------"] * len(conf_cols)) + "|"
        lines += [header, sep]
        for method in summary.index:
            vals = " | ".join(
                f"{summary.loc[method, c]:.3f}" if not np.isnan(summary.loc[method, c]) else "N/A"
                for c in conf_cols
            )
            lines.append(f"| {method} | {vals} |")
    lines.append("")

    # Ensemble diversity
    if not df_div.empty:
        lines += ["## 7. Ensemble CDR3 Diversity", ""]
        lines += [
            "| Method | CDR3 Mean (Å) | CDR3-H Mean (Å) | CDR3-L Mean (Å) | Frac ≥ 1.0 Å |",
            "|--------|--------------|-----------------|-----------------|---------------|",
        ]
        for method, group in df_div.groupby("method"):
            def _fmt(col):
                v = group[col].mean() if col in group.columns else np.nan
                return f"{v:.2f}" if not np.isnan(v) else "N/A"
            lines.append(
                f"| {method} | {_fmt('ensemble_mean_cdr3_rmsd')} | "
                f"{_fmt('ensemble_mean_cdr3h_rmsd')} | "
                f"{_fmt('ensemble_mean_cdr3l_rmsd')} | "
                f"{_fmt('ensemble_frac_above_1A')} |"
            )
        lines.append("")

    # Statistical tests
    if stats_results:
        lines += ["## 8. Statistical Tests (Wilcoxon signed-rank, paired)", ""]
        lines += ["| Baseline | Metric | Direction | n | Δ mean | Wins | Losses | p-value |"]
        lines += ["|----------|--------|-----------|---|--------|------|--------|---------|"]
        for r in stats_results:
            sig = "**" if r["p_value"] < 0.05 else ""
            lines.append(
                f"| {r['baseline']} | {r['metric']} | {r['direction']} | "
                f"{r['n_common']} | {r['delta_mean']:+.3f} | "
                f"{r['wins']} | {r['losses']} | "
                f"{sig}{r['p_value']:.4f}{sig} |"
            )
        lines.append("")

    # CAPRI classification
    lines += ["## 9. CAPRI Classification (dockq_antigen_ab, best model)", ""]
    lines += ["| Method | Incorrect | Acceptable | Medium | High |"]
    lines += ["|--------|-----------|------------|--------|------|"]
    for method, group in df.groupby("method"):
        classes = group["dockq_antigen_ab"].apply(capri_class)
        counts = classes.value_counts()
        lines.append(
            f"| {method} | {counts.get('Incorrect', 0)} | "
            f"{counts.get('Acceptable', 0)} | {counts.get('Medium', 0)} | "
            f"{counts.get('High', 0)} |"
        )

    report_path = output_dir / "evaluation_report.md"
    with open(report_path, "w") as f:
        f.write("\n".join(lines))
    print(f"  Saved report: {report_path}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Evaluate Strategy L+ (CDR3 Beta-Scaling v2) against baselines"
    )
    parser.add_argument(
        "--predictions_dir", default="predictions_examples",
        help="Root directory containing method subdirectories",
    )
    parser.add_argument(
        "--gt_dir", default="pdb_minimized",
        help="Directory containing ground truth PDB files",
    )
    parser.add_argument(
        "--cdrs_csv", default="examples/cdrs.csv",
        help="Path to CDR indices CSV file",
    )
    parser.add_argument(
        "--output_dir", default="evaluation_v2",
        help="Output directory for results, plots, and reports",
    )
    parser.add_argument(
        "--feature_folder", default=DEFAULT_FEATURE_FOLDER,
        help=(
            f"Subfolder name inside predictions_dir for the L+ v2 results "
            f"(default: {DEFAULT_FEATURE_FOLDER})"
        ),
    )
    parser.add_argument(
        "--feature_label", default=DEFAULT_FEATURE_LABEL,
        help=(
            f"Method label used in output tables and plots "
            f"(default: {DEFAULT_FEATURE_LABEL})"
        ),
    )
    args = parser.parse_args()

    warnings.filterwarnings("ignore", category=UserWarning)
    run_evaluation(args)


if __name__ == "__main__":
    main()
