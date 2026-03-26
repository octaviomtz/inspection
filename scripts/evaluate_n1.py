#!/usr/bin/env python3
"""
Evaluation script for N1: A+Y Hybrid (FK Particles + Hierarchical Timing) vs baselines.

Computes DockQ metrics and CDR-H3 RMSD for all methods and generates comparison
tables + plots.

Must be run with: conda activate dockq2 && python scripts/evaluate_n1.py

Methods:
  B0 - Unconstrained (antigen_cut)
  B1 - Contact restraints (antigen_cut_contact_restraints)
  B2 - Pocket restraints (antigen_cut_vhvl_msa_pocket_ab_boltz_post_2023_omm)
  N1 - Hybrid FK+Hierarchical (new_feature_hybrid_fk)
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
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon, spearmanr

from Bio.PDB import PDBParser, MMCIFParser, Superimposer
from Bio.PDB.Polypeptide import PPBuilder
from Bio.Align import PairwiseAligner

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

PROJECT_ROOT = Path(__file__).resolve().parent.parent

DEFAULT_GROUND_TRUTH_DIR = PROJECT_ROOT / "pdb_minimized"
DEFAULT_PREDICTIONS_DIR = PROJECT_ROOT / "predictions_examples"
DEFAULT_CDR_CSV = PROJECT_ROOT / "examples" / "cdrs.csv"

METHOD_SUBDIRS = {
    "B0": "antigen_cut",
    "B1": "antigen_cut_contact_restraints",
    "B2": "antigen_cut_vhvl_msa_pocket_ab_boltz_post_2023_omm",
    "N1": "new_feature_hybrid_fk",
}

METHOD_DIRS = {}
GROUND_TRUTH_DIR = DEFAULT_GROUND_TRUTH_DIR


def _build_method_dirs(predictions_dir: Path) -> dict:
    return {m: predictions_dir / subdir for m, subdir in METHOD_SUBDIRS.items()}


METHOD_LABELS = {
    "B0": "Unconstrained",
    "B1": "Contact restraints",
    "B2": "Pocket restraints",
    "N1": "Hybrid FK+Hier.",
}

CAPRI_THRESHOLDS = [
    ("Incorrect",  0.00, 0.23),
    ("Acceptable", 0.23, 0.49),
    ("Medium",     0.49, 0.80),
    ("High",       0.80, 1.01),
]

INTERFACES_AB_AC = ["AB", "AC"]
INTERFACE_BC = "BC"

ALL_METHODS = ["B0", "B1", "B2", "N1"]

METHOD_COLORS = {
    "B0": "#7f8c8d",   # gray
    "B1": "#2980b9",   # blue
    "B2": "#8e44ad",   # purple
    "N1": "#e74c3c",   # red
}

METHOD_ORDER = ["B0", "B1", "B2", "N1"]


# ---------------------------------------------------------------------------
# CDR loading
# ---------------------------------------------------------------------------

def load_cdr_definitions(cdr_csv_path):
    """
    Load CDR definitions from cdrs.csv.

    Returns dict: complex_id -> {
        "cdr3_h": [list of 1-indexed residue positions in heavy chain],
        "cdr3_l": [list of 1-indexed residue positions in light chain],
    }
    """
    cdrs = {}
    with open(cdr_csv_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            cid = row["complex"]
            cdr3_h = ast.literal_eval(row["cdr3_h"])
            cdr3_l = ast.literal_eval(row["cdr3_l"])
            cdrs[cid] = {
                "cdr3_h": cdr3_h,
                "cdr3_l": cdr3_l,
            }
    return cdrs


# ---------------------------------------------------------------------------
# Discovery: find prediction files for each method
# ---------------------------------------------------------------------------

def get_complex_ids():
    """Return sorted list of complex IDs from ground truth folder."""
    ids = []
    for f in sorted(GROUND_TRUTH_DIR.glob("*.pdb")):
        ids.append(f.stem)
    return ids


def _find_prediction_files(pred_dir):
    """Find all structure files (.pdb or .cif) under a predictions directory."""
    files = list(pred_dir.rglob("*_model_*.pdb")) + list(pred_dir.rglob("*_model_*.cif"))
    return sorted(files)


def _find_confidence_file(model_file):
    """Given a model file, find its corresponding confidence JSON."""
    d = model_file.parent
    name = model_file.stem
    conf_path = d / f"confidence_{name}.json"
    if conf_path.exists():
        return conf_path
    return None


def discover_b0(complex_ids):
    """
    B0: antigen_cut/boltz_results_{CID}/predictions/{CID}/
    One setting per complex, up to 5 models.
    """
    results = {}
    base = METHOD_DIRS["B0"]
    for cid in complex_ids:
        result_dir = base / f"boltz_results_{cid}" / "predictions" / cid
        if not result_dir.exists():
            continue
        models = _find_prediction_files(result_dir)
        if models:
            results[cid] = [{"setting": "unconstrained", "models": models}]
    return results


def discover_b1(complex_ids):
    """
    B1: antigen_cut_contact_restraints/
        boltz_results_restraint_{CID}_{type}_{N}/predictions/restraint_{CID}_{type}_{N}/
    Multiple settings per complex.
    """
    results = {}
    base = METHOD_DIRS["B1"]
    if not base.exists():
        return results

    for entry in sorted(base.iterdir()):
        if not entry.is_dir() or not entry.name.startswith("boltz_results_restraint_"):
            continue
        setting_id = entry.name.replace("boltz_results_", "")
        match = re.match(r"restraint_(\w{4}_\w{3})_(.*)", setting_id)
        if not match:
            continue
        cid = match.group(1)
        if cid not in complex_ids:
            continue
        pred_dir = entry / "predictions" / setting_id
        models = _find_prediction_files(pred_dir)
        if not models:
            continue
        if cid not in results:
            results[cid] = []
        results[cid].append({"setting": setting_id, "models": models})

    return results


def discover_b2(complex_ids):
    """
    B2: antigen_cut_vhvl_msa_pocket_ab_boltz_post_2023_omm/
        boltz_results_restraint_to_A_{CID}_{chain}_{res}_{num}/
        predictions/restraint_to_A_{CID}_{chain}_{res}_{num}/
    Multiple settings per complex.
    """
    results = {}
    base = METHOD_DIRS["B2"]
    if not base.exists():
        return results

    for entry in sorted(base.iterdir()):
        if not entry.is_dir() or not entry.name.startswith("boltz_results_restraint_to_A_"):
            continue
        setting_id = entry.name.replace("boltz_results_", "")
        match = re.match(r"restraint_to_A_(\w{4}_\w{3})_(.*)", setting_id)
        if not match:
            continue
        cid = match.group(1)
        if cid not in complex_ids:
            continue
        pred_dir = entry / "predictions" / setting_id
        models = _find_prediction_files(pred_dir)
        if not models:
            continue
        if cid not in results:
            results[cid] = []
        results[cid].append({"setting": setting_id, "models": models})

    return results


def discover_n1(complex_ids):
    """
    N1: new_feature_hybrid_fk/
    Supports layouts:
      - boltz_results_{CID}/predictions/{CID}/
      - boltz_results_{CID}_hybrid_fk/predictions/{CID}_hybrid_fk/
    One setting per complex.
    """
    results = {}
    base = METHOD_DIRS["N1"]
    if not base.exists():
        return results

    for cid in complex_ids:
        # Try suffix first
        for suffix in ["_hybrid_fk", ""]:
            result_dir = base / f"boltz_results_{cid}{suffix}" / "predictions" / f"{cid}{suffix}"
            if result_dir.exists():
                models = _find_prediction_files(result_dir)
                if models:
                    results[cid] = [{"setting": "hybrid_fk", "models": models}]
                    break
    return results


DISCOVER_FUNCS = {
    "B0": discover_b0,
    "B1": discover_b1,
    "B2": discover_b2,
    "N1": discover_n1,
}


def discover_all(complex_ids, methods=None):
    """Return dict: method -> {complex_id -> [{setting, models}]}."""
    if methods is None:
        methods = ALL_METHODS
    return {m: DISCOVER_FUNCS[m](complex_ids) for m in methods}


# ---------------------------------------------------------------------------
# DockQ execution
# ---------------------------------------------------------------------------

def run_dockq(model_path, native_path, output_json_path):
    """Run DockQ and write JSON output. Returns parsed JSON or None on failure."""
    cmd = [
        "DockQ", str(model_path), str(native_path),
        "--json", str(output_json_path),
    ]
    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True, timeout=120,
        )
        if result.returncode != 0:
            print(f"  WARNING: DockQ failed for {model_path.name}: "
                  f"{result.stderr.strip()[:200]}", file=sys.stderr)
            return None
        with open(output_json_path) as f:
            return json.load(f)
    except (subprocess.TimeoutExpired, FileNotFoundError) as e:
        print(f"  WARNING: DockQ error for {model_path.name}: {e}", file=sys.stderr)
        return None


def extract_interface_metrics(dockq_json):
    """
    Extract per-interface metrics from DockQ JSON output.
    Returns dict: {interface_key: {metric: value}} or None.
    """
    if dockq_json is None or "best_result" not in dockq_json:
        return None

    metrics = {}
    for iface_key, iface_data in dockq_json["best_result"].items():
        metrics[iface_key] = {
            "DockQ": iface_data.get("DockQ", np.nan),
            "iRMSD": iface_data.get("iRMSD", np.nan),
            "LRMSD": iface_data.get("LRMSD", np.nan),
            "fnat": iface_data.get("fnat", np.nan),
            "fnonnat": iface_data.get("fnonnat", np.nan),
            "F1": iface_data.get("F1", np.nan),
        }
    return metrics


def load_confidence(model_path):
    """Load confidence JSON for a model file. Returns dict or None."""
    conf_path = _find_confidence_file(model_path)
    if conf_path is None:
        return None
    try:
        with open(conf_path) as f:
            return json.load(f)
    except (json.JSONDecodeError, OSError):
        return None


# ---------------------------------------------------------------------------
# CDR-H3 RMSD computation
# ---------------------------------------------------------------------------

_pdb_parser = PDBParser(QUIET=True)
_cif_parser = MMCIFParser(QUIET=True)
_ppb = PPBuilder()


def _load_structure(path):
    """Load a PDB or CIF file into a BioPython structure."""
    path = Path(path)
    if path.suffix.lower() == ".cif":
        return _cif_parser.get_structure(path.stem, str(path))
    return _pdb_parser.get_structure(path.stem, str(path))


def _get_chain_sequence_and_ca(structure, chain_id):
    """
    Extract 1-letter sequence and CA atom coordinates for a chain.
    Returns:
        residues: list of Bio.PDB.Residue objects (protein only, ordered)
        seq: str, 1-letter amino acid sequence
        ca_atoms: list of Bio.PDB.Atom objects (CA)
    """
    from Bio.PDB.Polypeptide import protein_letters_3to1

    chain = structure[0][chain_id]
    residues = []
    seq_chars = []
    ca_atoms = []

    for res in chain.get_residues():
        # Skip non-protein residues
        if res.id[0] != " ":
            continue
        resname = res.get_resname().strip()
        if resname in protein_letters_3to1:
            one_letter = protein_letters_3to1[resname]
        else:
            one_letter = "X"
        if "CA" not in res:
            continue
        residues.append(res)
        seq_chars.append(one_letter)
        ca_atoms.append(res["CA"])

    return residues, "".join(seq_chars), ca_atoms


def _align_sequences(seq1, seq2):
    """
    Align two sequences and return index mapping: seq1_idx -> seq2_idx.
    Returns dict mapping 0-based indices from seq1 to seq2.
    """
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -0.5

    alignments = aligner.align(seq1, seq2)
    if not alignments:
        return {}

    aln = alignments[0]
    # aln.aligned returns (2, N_blocks, 2): blocks for seq1 and seq2
    blocks_s1 = aln.aligned[0]
    blocks_s2 = aln.aligned[1]

    mapping = {}
    for (s1_start, s1_end), (s2_start, s2_end) in zip(blocks_s1, blocks_s2):
        for i, j in zip(range(s1_start, s1_end), range(s2_start, s2_end)):
            mapping[i] = j

    return mapping


def compute_cdr_h3_rmsd(pred_path, gt_path, cdr3_h_positions):
    """
    Compute CDR-H3 RMSD after superimposing on heavy-chain framework.

    Args:
        pred_path: Path to predicted structure (PDB or CIF).
        gt_path: Path to ground truth PDB.
        cdr3_h_positions: list of 1-indexed residue positions for CDR-H3
                          in the heavy chain (from cdrs.csv).

    Returns:
        float: CDR-H3 CA RMSD in Angstroms, or np.nan on failure.
    """
    try:
        pred_struct = _load_structure(pred_path)
        gt_struct = _load_structure(gt_path)
    except Exception as e:
        print(f"  WARNING: Could not load structures for CDR-H3 RMSD: {e}",
              file=sys.stderr)
        return np.nan

    # Extract chain B (heavy chain) from both
    try:
        pred_res, pred_seq, pred_ca = _get_chain_sequence_and_ca(pred_struct, "B")
        gt_res, gt_seq, gt_ca = _get_chain_sequence_and_ca(gt_struct, "B")
    except (KeyError, Exception) as e:
        print(f"  WARNING: Could not extract chain B: {e}", file=sys.stderr)
        return np.nan

    if len(pred_ca) < 10 or len(gt_ca) < 10:
        return np.nan

    # Align prediction and GT sequences to handle different numbering/lengths
    mapping = _align_sequences(pred_seq, gt_seq)
    if not mapping:
        return np.nan

    # CDR3 positions are 1-indexed in the prediction sequence.
    # Convert to 0-indexed in the prediction residue list.
    cdr3_set = set(p - 1 for p in cdr3_h_positions)  # 0-indexed

    # Split into framework and CDR-H3 based on prediction indexing
    framework_pred_atoms = []
    framework_gt_atoms = []
    cdr3_pred_atoms = []
    cdr3_gt_atoms = []

    for pred_idx, gt_idx in mapping.items():
        if pred_idx >= len(pred_ca) or gt_idx >= len(gt_ca):
            continue
        p_atom = pred_ca[pred_idx]
        g_atom = gt_ca[gt_idx]

        if pred_idx in cdr3_set:
            cdr3_pred_atoms.append(p_atom)
            cdr3_gt_atoms.append(g_atom)
        else:
            framework_pred_atoms.append(p_atom)
            framework_gt_atoms.append(g_atom)

    if len(framework_pred_atoms) < 10 or len(cdr3_pred_atoms) < 3:
        return np.nan

    # Superimpose using framework CA atoms
    sup = Superimposer()
    sup.set_atoms(framework_gt_atoms, framework_pred_atoms)
    # Apply rotation/translation to ALL prediction atoms in the chain
    all_pred_atoms = [a for res in pred_struct[0]["B"].get_atoms() for a in [res] if False]
    # Actually just apply to the CDR-H3 prediction atoms
    sup.apply(list(pred_struct[0]["B"].get_atoms()))

    # Compute RMSD on CDR-H3 CA atoms after superposition
    deltas = []
    for p_atom, g_atom in zip(cdr3_pred_atoms, cdr3_gt_atoms):
        diff = p_atom.get_vector() - g_atom.get_vector()
        deltas.append(diff.norm() ** 2)

    if not deltas:
        return np.nan

    rmsd = np.sqrt(np.mean(deltas))
    return rmsd


# ---------------------------------------------------------------------------
# Core evaluation logic
# ---------------------------------------------------------------------------

def evaluate_models(model_paths, native_path, dockq_cache_dir):
    """
    Evaluate a list of model files against a native structure.
    Returns list of dicts, one per model, with DockQ metrics and confidence.
    """
    results = []
    for model_path in model_paths:
        parent_parts = model_path.parts
        try:
            idx = list(parent_parts).index("predictions_examples")
            unique_parts = parent_parts[idx + 1:]
        except ValueError:
            unique_parts = parent_parts[-4:]
        cache_key = "_".join(unique_parts).replace(".pdb", "").replace(".cif", "")
        json_out = dockq_cache_dir / f"dockq_{cache_key}.json"

        if json_out.exists():
            with open(json_out) as f:
                dockq_data = json.load(f)
        else:
            dockq_data = run_dockq(model_path, native_path, json_out)

        iface_metrics = extract_interface_metrics(dockq_data)
        conf_data = load_confidence(model_path)

        entry = {
            "model_path": str(model_path),
            "model_name": model_path.name,
            "interfaces": iface_metrics,
            "confidence": conf_data,
        }

        # Compute composite antibody-antigen DockQ
        if iface_metrics:
            ab_dockq = iface_metrics.get("AB", {}).get("DockQ", np.nan)
            ac_dockq = iface_metrics.get("AC", {}).get("DockQ", np.nan)
            vals = [v for v in [ab_dockq, ac_dockq] if not np.isnan(v)]
            entry["dockq_ab_ac"] = np.mean(vals) if vals else np.nan
            entry["dockq_bc"] = iface_metrics.get("BC", {}).get("DockQ", np.nan)
        else:
            entry["dockq_ab_ac"] = np.nan
            entry["dockq_bc"] = np.nan

        # Extract confidence score
        if conf_data:
            entry["confidence_score"] = conf_data.get("confidence_score", np.nan)
            entry["iptm"] = conf_data.get("iptm", np.nan)
            pair_iptm = conf_data.get("pair_chains_iptm", {})
            entry["iptm_ab"] = pair_iptm.get("0", {}).get("1", np.nan)
            entry["iptm_ac"] = pair_iptm.get("0", {}).get("2", np.nan)
            entry["iptm_bc"] = pair_iptm.get("1", {}).get("2", np.nan)
            entry["fk_energy"] = conf_data.get("fk_energy", np.nan)
            entry["composite_score"] = conf_data.get("composite_score", np.nan)
        else:
            entry["confidence_score"] = np.nan
            entry["iptm"] = np.nan
            entry["iptm_ab"] = np.nan
            entry["iptm_ac"] = np.nan
            entry["iptm_bc"] = np.nan
            entry["fk_energy"] = np.nan
            entry["composite_score"] = np.nan

        results.append(entry)

    return results


def select_best_model(model_results, criterion="confidence"):
    """
    Select the best model from a list of evaluated models.
    criterion: "confidence" (practical), "oracle" (best DockQ),
               or "composite" (FK energy + confidence).
    """
    if not model_results:
        return None

    if criterion == "oracle":
        valid = [m for m in model_results if not np.isnan(m["dockq_ab_ac"])]
        if not valid:
            return model_results[0]
        return max(valid, key=lambda m: m["dockq_ab_ac"])
    elif criterion == "composite":
        valid = [m for m in model_results
                 if not np.isnan(m.get("composite_score", np.nan))]
        if valid:
            return max(valid, key=lambda m: m["composite_score"])
        # Fall back to confidence
        valid = [m for m in model_results if not np.isnan(m["confidence_score"])]
        if not valid:
            return model_results[0]
        return max(valid, key=lambda m: m["confidence_score"])
    else:  # confidence
        valid = [m for m in model_results if not np.isnan(m["confidence_score"])]
        if not valid:
            return model_results[0]
        return max(valid, key=lambda m: m["confidence_score"])


def evaluate_method(method_id, method_data, complex_ids, dockq_cache_dir,
                    cdr_defs=None):
    """
    Evaluate all complexes for a single method.
    Returns dict: complex_id -> {
        "confidence_selected": best_model_entry,
        "oracle_selected": best_model_entry,
        "composite_selected": best_model_entry,
        "all_models": [...],
        "n_settings": int,
        "n_models": int,
        "cdr_h3_rmsd": float (on confidence-selected model),
    }
    """
    results = {}

    for cid in complex_ids:
        if cid not in method_data:
            continue

        native = GROUND_TRUTH_DIR / f"{cid}.pdb"
        if not native.exists():
            continue

        settings = method_data[cid]
        all_models = []

        for setting_info in settings:
            model_results = evaluate_models(
                setting_info["models"], native, dockq_cache_dir,
            )
            all_models.extend(model_results)

        if not all_models:
            continue

        conf_selected = select_best_model(all_models, "confidence")
        oracle_selected = select_best_model(all_models, "oracle")
        composite_selected = select_best_model(all_models, "composite")

        # Compute CDR-H3 RMSD on confidence-selected and oracle-selected models
        cdr_h3_rmsd_conf = np.nan
        cdr_h3_rmsd_oracle = np.nan
        if cdr_defs and cid in cdr_defs:
            cdr3_h = cdr_defs[cid]["cdr3_h"]
            if conf_selected and cdr3_h:
                cdr_h3_rmsd_conf = compute_cdr_h3_rmsd(
                    conf_selected["model_path"], str(native), cdr3_h,
                )
            if oracle_selected and cdr3_h:
                cdr_h3_rmsd_oracle = compute_cdr_h3_rmsd(
                    oracle_selected["model_path"], str(native), cdr3_h,
                )

        results[cid] = {
            "confidence_selected": conf_selected,
            "oracle_selected": oracle_selected,
            "composite_selected": composite_selected,
            "all_models": all_models,
            "n_settings": len(settings),
            "n_models": len(all_models),
            "cdr_h3_rmsd_conf": cdr_h3_rmsd_conf,
            "cdr_h3_rmsd_oracle": cdr_h3_rmsd_oracle,
        }

    return results


# ---------------------------------------------------------------------------
# Aggregation and reporting
# ---------------------------------------------------------------------------

def capri_class(dockq):
    """Classify DockQ score into CAPRI quality category."""
    if np.isnan(dockq):
        return "N/A"
    for label, lo, hi in CAPRI_THRESHOLDS:
        if lo <= dockq < hi:
            return label
    return "N/A"


def build_summary_table(all_results, selection="confidence", methods=None):
    """
    Build per-complex summary DataFrame.
    """
    if methods is None:
        methods = ALL_METHODS
    sel_key = f"{selection}_selected"
    rmsd_key = f"cdr_h3_rmsd_{selection}" if selection in ("conf", "oracle") else "cdr_h3_rmsd_conf"
    if selection == "oracle":
        rmsd_key = "cdr_h3_rmsd_oracle"
    elif selection in ("confidence", "composite"):
        rmsd_key = "cdr_h3_rmsd_conf"

    rows = []
    all_cids = set()
    for method_data in all_results.values():
        all_cids.update(method_data.keys())

    for cid in sorted(all_cids):
        row = {"complex_id": cid}
        for method_id in methods:
            prefix = method_id
            if method_id in all_results and cid in all_results[method_id]:
                entry = all_results[method_id][cid]
                best = entry[sel_key]
                if best and best["interfaces"]:
                    ab = best["interfaces"].get("AB", {})
                    ac = best["interfaces"].get("AC", {})
                    bc = best["interfaces"].get("BC", {})
                    row[f"{prefix}_DockQ_AB"] = ab.get("DockQ", np.nan)
                    row[f"{prefix}_DockQ_AC"] = ac.get("DockQ", np.nan)
                    row[f"{prefix}_DockQ_BC"] = bc.get("DockQ", np.nan)
                    row[f"{prefix}_DockQ_AbAg"] = best["dockq_ab_ac"]
                    row[f"{prefix}_iRMSD_AB"] = ab.get("iRMSD", np.nan)
                    row[f"{prefix}_LRMSD_AB"] = ab.get("LRMSD", np.nan)
                    row[f"{prefix}_fnat_AB"] = ab.get("fnat", np.nan)
                    row[f"{prefix}_F1_AB"] = ab.get("F1", np.nan)
                    row[f"{prefix}_iRMSD_AC"] = ac.get("iRMSD", np.nan)
                    row[f"{prefix}_fnat_AC"] = ac.get("fnat", np.nan)
                    row[f"{prefix}_confidence"] = best["confidence_score"]
                    row[f"{prefix}_iptm"] = best["iptm"]
                    row[f"{prefix}_fk_energy"] = best.get("fk_energy", np.nan)
                    row[f"{prefix}_composite_score"] = best.get("composite_score", np.nan)
                    row[f"{prefix}_n_models"] = entry["n_models"]
                    row[f"{prefix}_n_settings"] = entry["n_settings"]
                else:
                    for col in ["DockQ_AB", "DockQ_AC", "DockQ_BC", "DockQ_AbAg",
                                "iRMSD_AB", "LRMSD_AB", "fnat_AB", "F1_AB",
                                "iRMSD_AC", "fnat_AC", "confidence", "iptm",
                                "fk_energy", "composite_score"]:
                        row[f"{prefix}_{col}"] = np.nan
                    row[f"{prefix}_n_models"] = entry["n_models"]
                    row[f"{prefix}_n_settings"] = entry["n_settings"]

                # CDR-H3 RMSD
                row[f"{prefix}_CDR_H3_RMSD"] = entry.get(rmsd_key, np.nan)

        rows.append(row)

    return pd.DataFrame(rows)


def build_aggregate_table(df, methods=None):
    """Build method-level aggregate statistics from per-complex DataFrame."""
    if methods is None:
        methods = ALL_METHODS
    agg_rows = []
    for m in methods:
        row = {"Method": f"{m}: {METHOD_LABELS[m]}"}
        for col_suffix, label in [
            ("DockQ_AbAg", "DockQ_AbAg"), ("DockQ_AB", "DockQ_AB"),
            ("DockQ_AC", "DockQ_AC"), ("DockQ_BC", "DockQ_BC"),
            ("iRMSD_AB", "iRMSD_AB"), ("fnat_AB", "fnat_AB"),
            ("iptm", "iptm"), ("CDR_H3_RMSD", "CDR_H3_RMSD"),
        ]:
            col = f"{m}_{col_suffix}"
            if col in df.columns:
                vals = df[col].dropna()
                row[f"{label}_mean"] = vals.mean() if len(vals) > 0 else np.nan
                row[f"{label}_median"] = vals.median() if len(vals) > 0 else np.nan
                row[f"{label}_std"] = vals.std() if len(vals) > 0 else np.nan
                row[f"{label}_n"] = len(vals)
            else:
                row[f"{label}_mean"] = np.nan
                row[f"{label}_median"] = np.nan
                row[f"{label}_std"] = np.nan
                row[f"{label}_n"] = 0
        agg_rows.append(row)

    return pd.DataFrame(agg_rows)


def build_capri_table(df, methods=None):
    """Build CAPRI quality distribution table based on DockQ_AbAg."""
    if methods is None:
        methods = ALL_METHODS
    rows = []
    for m in methods:
        col = f"{m}_DockQ_AbAg"
        row = {"Method": f"{m}: {METHOD_LABELS[m]}"}
        if col not in df.columns:
            for label, _, _ in CAPRI_THRESHOLDS:
                row[label] = "N/A"
            row["N"] = 0
            rows.append(row)
            continue

        vals = df[col].dropna()
        n = len(vals)
        row["N"] = n
        for label, lo, hi in CAPRI_THRESHOLDS:
            count = ((vals >= lo) & (vals < hi)).sum()
            row[label] = f"{count}/{n} ({100*count/n:.0f}%)" if n > 0 else "N/A"
        rows.append(row)

    return pd.DataFrame(rows)


def run_statistical_tests(df, primary_nf="N1", methods=None):
    """
    Run Wilcoxon signed-rank tests: primary NF vs each baseline.
    Tests on both DockQ_AbAg and CDR_H3_RMSD.
    """
    if methods is None:
        methods = ALL_METHODS
    results = []
    comparisons = [m for m in methods if m != primary_nf]

    for metric, higher_better in [("DockQ_AbAg", True), ("CDR_H3_RMSD", False)]:
        nf_col = f"{primary_nf}_{metric}"
        if nf_col not in df.columns:
            continue

        for baseline in comparisons:
            bl_col = f"{baseline}_{metric}"
            if bl_col not in df.columns:
                continue

            paired = df[[nf_col, bl_col]].dropna()
            n = len(paired)
            if n < 5:
                results.append({
                    "Metric": metric,
                    "Comparison": f"{primary_nf} vs {baseline}",
                    "N_paired": n,
                    f"{primary_nf}_mean": paired[nf_col].mean() if n > 0 else np.nan,
                    "Baseline_mean": paired[bl_col].mean() if n > 0 else np.nan,
                    "Delta_mean": np.nan,
                    "Delta_median": np.nan,
                    "p_value": np.nan,
                    "significant": "N/A (too few pairs)",
                })
                continue

            delta = paired[nf_col] - paired[bl_col]
            try:
                stat, p = wilcoxon(delta, alternative="two-sided")
            except ValueError:
                stat, p = np.nan, np.nan

            results.append({
                "Metric": metric,
                "Comparison": f"{primary_nf} vs {baseline}",
                "N_paired": n,
                f"{primary_nf}_mean": paired[nf_col].mean(),
                "Baseline_mean": paired[bl_col].mean(),
                "Delta_mean": delta.mean(),
                "Delta_median": delta.median(),
                "p_value": p,
                "significant": "Yes" if p < 0.05 else "No",
            })

    return pd.DataFrame(results)


def compute_confidence_correlation(all_results, methods=None):
    """Compute Spearman correlation between confidence score and DockQ_AbAg."""
    if methods is None:
        methods = ALL_METHODS
    rows = []
    for method_id in methods:
        if method_id not in all_results:
            continue
        confs = []
        dockqs = []
        for cid, data in all_results[method_id].items():
            for model in data["all_models"]:
                c = model["confidence_score"]
                d = model["dockq_ab_ac"]
                if not np.isnan(c) and not np.isnan(d):
                    confs.append(c)
                    dockqs.append(d)
        if len(confs) >= 5:
            rho, p = spearmanr(confs, dockqs)
        else:
            rho, p = np.nan, np.nan
        rows.append({
            "Method": f"{method_id}: {METHOD_LABELS[method_id]}",
            "N_models": len(confs),
            "Spearman_rho": rho,
            "p_value": p,
        })

    # FK energy correlation
    for method_id in methods:
        if method_id not in all_results:
            continue
        energies = []
        dockqs = []
        for cid, data in all_results[method_id].items():
            for model in data["all_models"]:
                e = model.get("fk_energy", np.nan)
                d = model["dockq_ab_ac"]
                if not np.isnan(e) and not np.isnan(d):
                    energies.append(e)
                    dockqs.append(d)
        if len(energies) >= 5:
            rho, p = spearmanr(energies, dockqs)
            rows.append({
                "Method": f"{method_id} (FK energy): {METHOD_LABELS[method_id]}",
                "N_models": len(energies),
                "Spearman_rho": rho,
                "p_value": p,
            })

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def _available_methods(df, metric_suffix="DockQ_AbAg"):
    return [m for m in METHOD_ORDER if f"{m}_{metric_suffix}" in df.columns
            and df[f"{m}_{metric_suffix}"].notna().any()]


def _savefig(fig, path):
    fig.savefig(path, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  Plot saved: {path.name}")


def plot_aggregate_bars(df_conf, df_oracle, output_dir):
    """Bar chart comparing mean DockQ_AbAg across methods, confidence vs oracle."""
    methods = _available_methods(df_conf)
    if not methods:
        return

    conf_means = [df_conf[f"{m}_DockQ_AbAg"].dropna().mean() for m in methods]
    orac_means = [df_oracle[f"{m}_DockQ_AbAg"].dropna().mean() for m in methods]
    labels = [f"{m}\n{METHOD_LABELS[m]}" for m in methods]
    colors = [METHOD_COLORS[m] for m in methods]

    x = np.arange(len(methods))
    w = 0.35

    fig, ax = plt.subplots(figsize=(max(8, len(methods) * 1.8), 5))
    bars1 = ax.bar(x - w / 2, conf_means, w, label="Confidence selection",
                   color=colors, edgecolor="white", linewidth=0.8)
    bars2 = ax.bar(x + w / 2, orac_means, w, label="Oracle selection",
                   color=colors, edgecolor="white", linewidth=0.8, alpha=0.55,
                   hatch="//")

    for bar in list(bars1) + list(bars2):
        h = bar.get_height()
        ax.text(bar.get_x() + bar.get_width() / 2, h + 0.005, f"{h:.3f}",
                ha="center", va="bottom", fontsize=8)

    ax.set_ylabel("Mean DockQ (Ab-Ag)")
    ax.set_title("Aggregate DockQ by Method")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=9)
    ax.legend(frameon=False)
    ax.set_ylim(0, max(max(conf_means), max(orac_means)) * 1.25)
    ax.spines[["top", "right"]].set_visible(False)

    _savefig(fig, output_dir / "plot_aggregate_dockq.png")


def plot_capri_stacked(df_conf, df_oracle, output_dir):
    """Stacked bar chart of CAPRI quality distribution."""
    tier_colors = {
        "Incorrect":  "#e74c3c",
        "Acceptable": "#f39c12",
        "Medium":     "#2ecc71",
        "High":       "#2980b9",
    }

    for selection, df in [("confidence", df_conf), ("oracle", df_oracle)]:
        methods = _available_methods(df)
        if not methods:
            continue

        fig, ax = plt.subplots(figsize=(max(7, len(methods) * 1.6), 5))
        labels = [f"{m}\n{METHOD_LABELS[m]}" for m in methods]
        x = np.arange(len(methods))
        bottoms = np.zeros(len(methods))

        for tier_label, lo, hi in CAPRI_THRESHOLDS:
            fracs = []
            for m in methods:
                vals = df[f"{m}_DockQ_AbAg"].dropna()
                n = len(vals)
                count = ((vals >= lo) & (vals < hi)).sum()
                fracs.append(count / n * 100 if n > 0 else 0)
            fracs = np.array(fracs)
            bars = ax.bar(x, fracs, 0.6, bottom=bottoms, label=tier_label,
                          color=tier_colors[tier_label], edgecolor="white",
                          linewidth=0.5)
            for i, (frac, bot) in enumerate(zip(fracs, bottoms)):
                if frac >= 8:
                    ax.text(x[i], bot + frac / 2, f"{frac:.0f}%",
                            ha="center", va="center", fontsize=8,
                            color="white", fontweight="bold")
            bottoms += fracs

        ax.set_ylabel("Percentage of complexes")
        ax.set_title(f"CAPRI Quality Distribution ({selection} selection)")
        ax.set_xticks(x)
        ax.set_xticklabels(labels, fontsize=9)
        ax.set_ylim(0, 105)
        ax.legend(loc="upper right", frameon=False, fontsize=8)
        ax.spines[["top", "right"]].set_visible(False)

        _savefig(fig, output_dir / f"plot_capri_{selection}.png")


def plot_boxplot_dockq(df_conf, df_oracle, output_dir):
    """Box plots of per-complex DockQ_AbAg distribution per method."""
    for selection, df in [("confidence", df_conf), ("oracle", df_oracle)]:
        methods = _available_methods(df)
        if not methods:
            continue

        data = []
        labels = []
        colors = []
        for m in methods:
            vals = df[f"{m}_DockQ_AbAg"].dropna().values
            data.append(vals)
            labels.append(f"{m}\n{METHOD_LABELS[m]}")
            colors.append(METHOD_COLORS[m])

        fig, ax = plt.subplots(figsize=(max(7, len(methods) * 1.6), 5))
        bp = ax.boxplot(data, patch_artist=True, widths=0.5,
                        medianprops=dict(color="black", linewidth=1.5))
        for patch, color in zip(bp["boxes"], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)

        # Overlay individual points
        for i, vals in enumerate(data):
            jitter = np.random.default_rng(42).uniform(-0.12, 0.12, len(vals))
            ax.scatter(np.full(len(vals), i + 1) + jitter, vals,
                       s=12, alpha=0.4, color=colors[i], edgecolors="none", zorder=3)

        ax.set_xticklabels(labels, fontsize=9)
        ax.set_ylabel("DockQ (Ab-Ag)")
        ax.set_title(f"DockQ Distribution ({selection} selection)")
        ax.spines[["top", "right"]].set_visible(False)

        _savefig(fig, output_dir / f"plot_boxplot_{selection}.png")


def plot_paired_scatter(df, output_dir, nf_method="N1"):
    """Scatter NF vs B0 per-complex DockQ_AbAg (confidence selection)."""
    nf_col = f"{nf_method}_DockQ_AbAg"
    b0_col = "B0_DockQ_AbAg"
    if nf_col not in df.columns or b0_col not in df.columns:
        return
    paired = df[[b0_col, nf_col]].dropna()
    if len(paired) < 3:
        return

    fig, ax = plt.subplots(figsize=(6, 6))
    lim = max(paired[b0_col].max(), paired[nf_col].max()) * 1.1
    lim = max(lim, 0.1)
    ax.plot([0, lim], [0, lim], "k--", alpha=0.3, linewidth=1, zorder=0)
    ax.scatter(paired[b0_col], paired[nf_col],
               s=30, alpha=0.7, color=METHOD_COLORS[nf_method], edgecolors="white",
               linewidth=0.4, zorder=2)

    n_above = (paired[nf_col] > paired[b0_col]).sum()
    n_below = (paired[nf_col] < paired[b0_col]).sum()
    n_equal = (paired[nf_col] == paired[b0_col]).sum()
    ax.text(0.05, 0.92, f"{nf_method} wins: {n_above}\nB0 wins: {n_below}\nTied: {n_equal}",
            transform=ax.transAxes, fontsize=9, va="top",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="wheat", alpha=0.5))

    ax.set_xlabel(f"B0 ({METHOD_LABELS['B0']}) DockQ Ab-Ag")
    ax.set_ylabel(f"{nf_method} ({METHOD_LABELS[nf_method]}) DockQ Ab-Ag")
    ax.set_title(f"Per-Complex: {nf_method} vs B0 (confidence selection)")
    ax.set_xlim(-0.02, lim)
    ax.set_ylim(-0.02, lim)
    ax.set_aspect("equal")
    ax.spines[["top", "right"]].set_visible(False)

    _savefig(fig, output_dir / f"plot_scatter_{nf_method.lower()}_vs_b0.png")


def plot_delta_waterfall(df, output_dir, nf_method="N1"):
    """Waterfall chart of per-complex DockQ delta (NF - B0), sorted."""
    nf_col = f"{nf_method}_DockQ_AbAg"
    b0_col = "B0_DockQ_AbAg"
    if nf_col not in df.columns or b0_col not in df.columns:
        return
    paired = df[["complex_id", b0_col, nf_col]].dropna()
    if len(paired) < 3:
        return

    paired = paired.copy()
    paired["delta"] = paired[nf_col] - paired[b0_col]
    paired = paired.sort_values("delta")

    fig, ax = plt.subplots(figsize=(max(8, len(paired) * 0.22), 5))
    colors = ["#2ecc71" if d > 0 else "#e74c3c" if d < 0 else "#95a5a6"
              for d in paired["delta"]]
    ax.bar(range(len(paired)), paired["delta"], color=colors, edgecolor="white",
           linewidth=0.3)
    ax.axhline(0, color="black", linewidth=0.8)

    ax.set_xticks(range(len(paired)))
    ax.set_xticklabels(paired["complex_id"], rotation=90, fontsize=6)
    ax.set_ylabel(f"Delta DockQ Ab-Ag ({nf_method} - B0)")
    ax.set_title(f"Per-Complex Improvement: {nf_method} vs B0 (confidence selection)")
    ax.spines[["top", "right"]].set_visible(False)

    mean_delta = paired["delta"].mean()
    ax.axhline(mean_delta, color=METHOD_COLORS[nf_method], linestyle="--",
               linewidth=1, alpha=0.7)
    ax.text(len(paired) - 1, mean_delta, f" mean={mean_delta:+.3f}",
            va="bottom", fontsize=8, color=METHOD_COLORS[nf_method])

    _savefig(fig, output_dir / f"plot_delta_{nf_method.lower()}_vs_b0.png")


def plot_confidence_vs_dockq(all_results, output_dir):
    """Scatter plot of confidence score vs DockQ_AbAg for each method."""
    methods = [m for m in METHOD_ORDER if m in all_results]
    if not methods:
        return

    n = len(methods)
    fig, axes = plt.subplots(1, n, figsize=(4.5 * n, 4.5), squeeze=False)
    axes = axes[0]

    for ax, method_id in zip(axes, methods):
        confs, dockqs = [], []
        for cid, data in all_results[method_id].items():
            for model in data["all_models"]:
                c = model["confidence_score"]
                d = model["dockq_ab_ac"]
                if not np.isnan(c) and not np.isnan(d):
                    confs.append(c)
                    dockqs.append(d)

        if len(confs) < 5:
            ax.set_title(f"{method_id}: {METHOD_LABELS[method_id]}\n(too few points)")
            continue

        ax.scatter(confs, dockqs, s=15, alpha=0.45, color=METHOD_COLORS[method_id],
                   edgecolors="none")
        rho, p = spearmanr(confs, dockqs)
        ax.set_xlabel("Confidence score")
        ax.set_ylabel("DockQ Ab-Ag")
        ax.set_title(f"{method_id}: {METHOD_LABELS[method_id]}\n"
                     f"rho={rho:.3f}, p={p:.1e}, n={len(confs)}")
        ax.spines[["top", "right"]].set_visible(False)

    fig.suptitle("Confidence Calibration", fontsize=13, y=1.02)
    fig.tight_layout()
    _savefig(fig, output_dir / "plot_confidence_calibration.png")


def plot_cdr_h3_rmsd_bars(df_conf, df_oracle, output_dir):
    """Bar chart comparing mean CDR-H3 RMSD across methods."""
    for selection, df in [("confidence", df_conf), ("oracle", df_oracle)]:
        methods = _available_methods(df, "CDR_H3_RMSD")
        if not methods:
            continue

        means = [df[f"{m}_CDR_H3_RMSD"].dropna().mean() for m in methods]
        stds = [df[f"{m}_CDR_H3_RMSD"].dropna().std() for m in methods]
        labels = [f"{m}\n{METHOD_LABELS[m]}" for m in methods]
        colors = [METHOD_COLORS[m] for m in methods]

        fig, ax = plt.subplots(figsize=(max(7, len(methods) * 1.6), 5))
        bars = ax.bar(range(len(methods)), means, 0.55, yerr=stds,
                      color=colors, edgecolor="white", linewidth=0.8,
                      capsize=4, error_kw=dict(lw=1.2))
        for bar, val in zip(bars, means):
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.15,
                    f"{val:.2f}", ha="center", va="bottom", fontsize=9)

        ax.set_xticks(range(len(methods)))
        ax.set_xticklabels(labels, fontsize=9)
        ax.set_ylabel("Mean CDR-H3 RMSD (A, lower is better)")
        ax.set_title(f"CDR-H3 RMSD by Method ({selection} selection)")
        ax.spines[["top", "right"]].set_visible(False)

        _savefig(fig, output_dir / f"plot_cdr_h3_rmsd_{selection}.png")


def plot_cdr_h3_rmsd_scatter(df, output_dir, nf_method="N1"):
    """Scatter N1 vs B0 per-complex CDR-H3 RMSD (confidence selection)."""
    nf_col = f"{nf_method}_CDR_H3_RMSD"
    b0_col = "B0_CDR_H3_RMSD"
    if nf_col not in df.columns or b0_col not in df.columns:
        return
    paired = df[[b0_col, nf_col]].dropna()
    if len(paired) < 3:
        return

    fig, ax = plt.subplots(figsize=(6, 6))
    lim = max(paired[b0_col].max(), paired[nf_col].max()) * 1.1
    lim = max(lim, 0.5)
    ax.plot([0, lim], [0, lim], "k--", alpha=0.3, linewidth=1, zorder=0)
    ax.scatter(paired[b0_col], paired[nf_col],
               s=30, alpha=0.7, color=METHOD_COLORS[nf_method], edgecolors="white",
               linewidth=0.4, zorder=2)

    # For RMSD, below diagonal = improvement (lower is better)
    n_below = (paired[nf_col] < paired[b0_col]).sum()
    n_above = (paired[nf_col] > paired[b0_col]).sum()
    n_equal = (paired[nf_col] == paired[b0_col]).sum()
    ax.text(0.05, 0.92,
            f"{nf_method} better: {n_below}\nB0 better: {n_above}\nTied: {n_equal}",
            transform=ax.transAxes, fontsize=9, va="top",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="wheat", alpha=0.5))

    ax.set_xlabel(f"B0 ({METHOD_LABELS['B0']}) CDR-H3 RMSD (A)")
    ax.set_ylabel(f"{nf_method} ({METHOD_LABELS[nf_method]}) CDR-H3 RMSD (A)")
    ax.set_title(f"CDR-H3 RMSD: {nf_method} vs B0 (confidence selection)")
    ax.set_xlim(-0.1, lim)
    ax.set_ylim(-0.1, lim)
    ax.set_aspect("equal")
    ax.spines[["top", "right"]].set_visible(False)

    _savefig(fig, output_dir / f"plot_cdr_h3_rmsd_scatter_{nf_method.lower()}_vs_b0.png")


def plot_irmsd_bars(df_conf, output_dir):
    """Bar chart comparing mean iRMSD_AB across methods (confidence selection)."""
    methods = [m for m in METHOD_ORDER if f"{m}_iRMSD_AB" in df_conf.columns
               and df_conf[f"{m}_iRMSD_AB"].notna().any()]
    if not methods:
        return

    means = [df_conf[f"{m}_iRMSD_AB"].dropna().mean() for m in methods]
    labels = [f"{m}\n{METHOD_LABELS[m]}" for m in methods]
    colors = [METHOD_COLORS[m] for m in methods]

    fig, ax = plt.subplots(figsize=(max(7, len(methods) * 1.6), 5))
    bars = ax.bar(range(len(methods)), means, 0.55, color=colors,
                  edgecolor="white", linewidth=0.8)
    for bar, val in zip(bars, means):
        ax.text(bar.get_x() + bar.get_width() / 2, val + 0.15, f"{val:.1f}",
                ha="center", va="bottom", fontsize=9)

    ax.set_xticks(range(len(methods)))
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel("Mean iRMSD (Ab chain, lower is better)")
    ax.set_title("Interface RMSD by Method (confidence selection)")
    ax.spines[["top", "right"]].set_visible(False)

    _savefig(fig, output_dir / "plot_irmsd.png")


def plot_fnat_bars(df_conf, output_dir):
    """Bar chart comparing mean fnat_AB across methods (confidence selection)."""
    methods = [m for m in METHOD_ORDER if f"{m}_fnat_AB" in df_conf.columns
               and df_conf[f"{m}_fnat_AB"].notna().any()]
    if not methods:
        return

    means = [df_conf[f"{m}_fnat_AB"].dropna().mean() for m in methods]
    labels = [f"{m}\n{METHOD_LABELS[m]}" for m in methods]
    colors = [METHOD_COLORS[m] for m in methods]

    fig, ax = plt.subplots(figsize=(max(7, len(methods) * 1.6), 5))
    bars = ax.bar(range(len(methods)), means, 0.55, color=colors,
                  edgecolor="white", linewidth=0.8)
    for bar, val in zip(bars, means):
        ax.text(bar.get_x() + bar.get_width() / 2, val + 0.005, f"{val:.3f}",
                ha="center", va="bottom", fontsize=9)

    ax.set_xticks(range(len(methods)))
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel("Mean fnat (fraction native contacts)")
    ax.set_title("Native Contact Recovery by Method (confidence selection)")
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_ylim(0, max(means) * 1.3)

    _savefig(fig, output_dir / "plot_fnat.png")


def plot_per_complex_heatmap(df, output_dir):
    """Heatmap of DockQ_AbAg per complex per method (confidence selection)."""
    methods = _available_methods(df)
    if len(methods) < 2:
        return

    cols = [f"{m}_DockQ_AbAg" for m in methods]
    sub = df[["complex_id"] + cols].dropna(subset=cols, how="all").copy()
    if len(sub) < 3:
        return

    sort_col = cols[0]
    sub = sub.sort_values(sort_col, ascending=False)

    matrix = sub[cols].values
    labels_y = sub["complex_id"].values
    labels_x = [f"{m}: {METHOD_LABELS[m]}" for m in methods]

    fig_height = max(5, len(labels_y) * 0.25)
    fig, ax = plt.subplots(figsize=(3 + len(methods) * 1.5, fig_height))
    im = ax.imshow(matrix, aspect="auto", cmap="RdYlGn", vmin=0, vmax=1)

    ax.set_xticks(range(len(labels_x)))
    ax.set_xticklabels(labels_x, fontsize=9, rotation=30, ha="right")
    ax.set_yticks(range(len(labels_y)))
    ax.set_yticklabels(labels_y, fontsize=6)
    ax.set_title("DockQ Ab-Ag per Complex (confidence selection)")

    cbar = fig.colorbar(im, ax=ax, shrink=0.6, pad=0.02)
    cbar.set_label("DockQ Ab-Ag", fontsize=9)

    _savefig(fig, output_dir / "plot_heatmap_dockq.png")


def plot_fk_energy_vs_dockq(all_results, output_dir):
    """Scatter plot of FK energy vs DockQ for methods that have FK energy."""
    methods_with_energy = []
    for m in METHOD_ORDER:
        if m not in all_results:
            continue
        has_energy = False
        for cid, data in all_results[m].items():
            for model in data["all_models"]:
                if not np.isnan(model.get("fk_energy", np.nan)):
                    has_energy = True
                    break
            if has_energy:
                break
        if has_energy:
            methods_with_energy.append(m)

    if not methods_with_energy:
        return

    n = len(methods_with_energy)
    fig, axes = plt.subplots(1, n, figsize=(5 * n, 4.5), squeeze=False)
    axes = axes[0]

    for ax, method_id in zip(axes, methods_with_energy):
        energies, dockqs = [], []
        for cid, data in all_results[method_id].items():
            for model in data["all_models"]:
                e = model.get("fk_energy", np.nan)
                d = model["dockq_ab_ac"]
                if not np.isnan(e) and not np.isnan(d):
                    energies.append(e)
                    dockqs.append(d)

        if len(energies) < 5:
            ax.set_title(f"{method_id}: FK energy vs DockQ\n(too few points)")
            continue

        ax.scatter(energies, dockqs, s=15, alpha=0.45, color=METHOD_COLORS[method_id],
                   edgecolors="none")
        rho, p = spearmanr(energies, dockqs)
        ax.set_xlabel("FK Energy (lower = better)")
        ax.set_ylabel("DockQ Ab-Ag")
        ax.set_title(f"{method_id}: {METHOD_LABELS[method_id]}\n"
                     f"rho={rho:.3f}, p={p:.1e}, n={len(energies)}")
        ax.spines[["top", "right"]].set_visible(False)

    fig.suptitle("FK Energy vs DockQ", fontsize=13, y=1.02)
    fig.tight_layout()
    _savefig(fig, output_dir / "plot_fk_energy_vs_dockq.png")


def generate_all_plots(all_results, output_dir):
    """Generate all evaluation plots and save to output_dir."""
    print(f"\n{'='*70}")
    print("  GENERATING PLOTS")
    print(f"{'='*70}")

    df_conf = build_summary_table(all_results, "confidence")
    df_oracle = build_summary_table(all_results, "oracle")

    plot_aggregate_bars(df_conf, df_oracle, output_dir)
    plot_capri_stacked(df_conf, df_oracle, output_dir)
    plot_boxplot_dockq(df_conf, df_oracle, output_dir)
    plot_paired_scatter(df_conf, output_dir, "N1")
    plot_delta_waterfall(df_conf, output_dir, "N1")
    plot_confidence_vs_dockq(all_results, output_dir)
    plot_fk_energy_vs_dockq(all_results, output_dir)
    plot_cdr_h3_rmsd_bars(df_conf, df_oracle, output_dir)
    plot_cdr_h3_rmsd_scatter(df_conf, output_dir, "N1")
    plot_irmsd_bars(df_conf, output_dir)
    plot_fnat_bars(df_conf, output_dir)
    plot_per_complex_heatmap(df_conf, output_dir)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Evaluate N1: Hybrid FK+Hierarchical vs baselines using DockQ v2.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--predictions-dir", "-p", type=Path, default=DEFAULT_PREDICTIONS_DIR,
        help="Parent directory containing per-method prediction subdirectories "
             "(default: predictions_examples/)",
    )
    parser.add_argument(
        "--ground-truth-dir", "-g", type=Path, default=DEFAULT_GROUND_TRUTH_DIR,
        help="Directory containing ground truth PDB files (default: pdb_minimized/)",
    )
    parser.add_argument(
        "--output-dir", "-o", type=Path,
        default=PROJECT_ROOT / "evaluation_results_n1",
        help="Directory to write evaluation outputs (default: evaluation_results_n1/)",
    )
    parser.add_argument(
        "--cdr-csv", type=Path, default=DEFAULT_CDR_CSV,
        help="Path to CDR definitions CSV (default: examples/cdrs.csv)",
    )
    parser.add_argument(
        "--methods", nargs="+", default=ALL_METHODS,
        choices=ALL_METHODS,
        help="Methods to evaluate (default: all)",
    )
    parser.add_argument(
        "--complexes", nargs="+", default=None,
        help="Specific complex IDs to evaluate (default: all available)",
    )
    parser.add_argument(
        "--skip-dockq", action="store_true",
        help="Skip DockQ computation, use cached results only",
    )
    parser.add_argument(
        "--skip-cdr-rmsd", action="store_true",
        help="Skip CDR-H3 RMSD computation",
    )
    parser.add_argument(
        "--primary-nf", default="N1", choices=ALL_METHODS,
        help="Primary new feature method for statistical tests (default: N1)",
    )
    args = parser.parse_args()

    # Set global dirs from CLI args
    global METHOD_DIRS, GROUND_TRUTH_DIR
    METHOD_DIRS = _build_method_dirs(args.predictions_dir.resolve())
    GROUND_TRUTH_DIR = args.ground_truth_dir.resolve()

    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    dockq_cache_dir = output_dir / "dockq_cache"
    dockq_cache_dir.mkdir(parents=True, exist_ok=True)

    methods = args.methods

    print(f"Predictions dir: {args.predictions_dir.resolve()}")
    print(f"Ground truth dir: {GROUND_TRUTH_DIR}")

    # Load CDR definitions
    cdr_defs = None
    if not args.skip_cdr_rmsd and args.cdr_csv.exists():
        cdr_defs = load_cdr_definitions(args.cdr_csv)
        print(f"CDR definitions loaded: {len(cdr_defs)} complexes from {args.cdr_csv}")
    elif not args.skip_cdr_rmsd:
        print(f"WARNING: CDR CSV not found at {args.cdr_csv}, skipping CDR-H3 RMSD",
              file=sys.stderr)

    # 1. Discover complexes and predictions
    all_complex_ids = get_complex_ids()
    if args.complexes:
        all_complex_ids = [c for c in all_complex_ids if c in args.complexes]

    print(f"Ground truth complexes: {len(all_complex_ids)}")

    discoveries = discover_all(all_complex_ids, methods)

    for method_id in methods:
        method_data = discoveries[method_id]
        n_complexes = len(method_data)
        n_models = sum(
            sum(len(s["models"]) for s in settings)
            for settings in method_data.values()
        )
        print(f"  {method_id} ({METHOD_LABELS[method_id]}): "
              f"{n_complexes} complexes, {n_models} models")

    # 2. Run DockQ evaluation (+ CDR-H3 RMSD)
    all_results = {}
    for method_id in methods:
        method_data = discoveries[method_id]
        if not method_data:
            print(f"\n  Skipping {method_id}: no predictions found")
            continue

        print(f"\nEvaluating {method_id} ({METHOD_LABELS[method_id]})...")

        if args.skip_dockq:
            print("  (using cached DockQ results only)")

        all_results[method_id] = evaluate_method(
            method_id, method_data, all_complex_ids, dockq_cache_dir,
            cdr_defs=cdr_defs,
        )
        print(f"  Evaluated {len(all_results[method_id])} complexes")

    if not all_results:
        print("\nNo predictions found for any method. Nothing to evaluate.")
        sys.exit(0)

    # 3. Build tables for each selection criterion
    for selection in ["confidence", "oracle"]:
        print(f"\n{'='*70}")
        print(f"  MODEL SELECTION: {selection.upper()}")
        print(f"{'='*70}")

        df = build_summary_table(all_results, selection, methods)
        csv_path = output_dir / f"per_complex_{selection}.csv"
        df.to_csv(csv_path, index=False, float_format="%.4f")
        print(f"\nPer-complex results saved to: {csv_path}")

        # Aggregate table
        agg_df = build_aggregate_table(df, methods)
        agg_csv = output_dir / f"aggregate_{selection}.csv"
        agg_df.to_csv(agg_csv, index=False, float_format="%.4f")

        print(f"\n--- Aggregate Results ({selection} selection) ---")
        display_cols = ["Method", "DockQ_AbAg_mean", "DockQ_AbAg_median",
                        "DockQ_AB_mean", "DockQ_AC_mean", "DockQ_BC_mean",
                        "iRMSD_AB_mean", "fnat_AB_mean", "iptm_mean",
                        "CDR_H3_RMSD_mean", "DockQ_AbAg_n"]
        display_cols = [c for c in display_cols if c in agg_df.columns]
        print(agg_df[display_cols].to_string(index=False,
              float_format=lambda x: f"{x:.4f}"))

        # CAPRI table
        capri_df = build_capri_table(df, methods)
        capri_csv = output_dir / f"capri_{selection}.csv"
        capri_df.to_csv(capri_csv, index=False)
        print(f"\n--- CAPRI Quality Distribution ({selection} selection) ---")
        print(capri_df.to_string(index=False))

        # Statistical tests
        stat_df = run_statistical_tests(df, primary_nf=args.primary_nf,
                                        methods=methods)
        if not stat_df.empty:
            stat_csv = output_dir / f"statistical_tests_{selection}.csv"
            stat_df.to_csv(stat_csv, index=False, float_format="%.4f")
            print(f"\n--- Statistical Tests ({selection} selection) ---")
            print(stat_df.to_string(index=False,
                  float_format=lambda x: f"{x:.4f}"))

    # 4. Confidence calibration
    print(f"\n{'='*70}")
    print("  CONFIDENCE CALIBRATION")
    print(f"{'='*70}")
    corr_df = compute_confidence_correlation(all_results, methods)
    if not corr_df.empty:
        corr_csv = output_dir / "confidence_correlation.csv"
        corr_df.to_csv(corr_csv, index=False, float_format="%.4f")
        print(corr_df.to_string(index=False, float_format=lambda x: f"{x:.4f}"))

    # 5. Generate plots
    generate_all_plots(all_results, output_dir)

    # 6. Summary
    print(f"\n{'='*70}")
    print(f"  All results saved to: {output_dir}/")
    print(f"{'='*70}")
    print("Files:")
    for f in sorted(output_dir.glob("*.csv")) + sorted(output_dir.glob("*.png")):
        print(f"  {f.name}")


if __name__ == "__main__":
    main()
