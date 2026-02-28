#!/usr/bin/env python
"""
evaluate_k_plus.py — Evaluation script for K+ region-specific beta-scaling.

Compares K+ against baselines (B1/B2/B3) on two axes:
  Axis A: Epitope prediction accuracy (Precision, Recall, F1, MCC)
  Axis B: Docking / interface quality (DockQ, iRMSD, LRMSD, fnat, ab-aligned ag RMSD)

Usage:
  conda run -n boltz python evaluate_k_plus.py --predictions_dir predictions_examples --gt_dir pdb_minimized
"""

import argparse
import json
import logging
import os
import re
import subprocess
import sys
import math
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

from Bio.PDB import PDBParser, MMCIFParser

try:
    from scipy.spatial import cKDTree
except ImportError:
    cKDTree = None

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

METHOD_DIRS = {
    "B1": "antigen_cut",
    "B2": "antigen_cut_contact_restraints",
    "B3": "antigen_cut_vhvl_msa_pocket_ab_boltz_post_2023_omm",
    "K+": "new_feature_cdr3_beta",
}

DOCKQ_QUALITY = {
    "Incorrect": (0.0, 0.23),
    "Acceptable": (0.23, 0.49),
    "Medium": (0.49, 0.80),
    "High": (0.80, 1.01),
}

ANTIGEN_CHAIN = "A"
HEAVY_CHAIN = "B"
LIGHT_CHAIN = "C"
ANTIBODY_CHAINS = (HEAVY_CHAIN, LIGHT_CHAIN)

# ---------------------------------------------------------------------------
# 1. Data Discovery
# ---------------------------------------------------------------------------


def _extract_complex_name_b1(folder_name):
    """boltz_results_{COMPLEX} -> COMPLEX"""
    m = re.match(r"^boltz_results_(.+)$", folder_name)
    return m.group(1) if m else None


def _extract_complex_name_b2(folder_name):
    """boltz_results_restraint_{COMPLEX}_{TYPE}_{POS} -> COMPLEX, config_label"""
    m = re.match(r"^boltz_results_restraint_(.+?)_(hbond|hydrophobic|salt_bridge)_(\d+)$", folder_name)
    if m:
        return m.group(1), f"{m.group(2)}_{m.group(3)}"
    return None, None


def _extract_complex_name_b3(folder_name):
    """boltz_results_restraint_to_A_{COMPLEX}_{CHAIN}_{RES}_{POS} -> COMPLEX, config_label"""
    m = re.match(r"^boltz_results_restraint_to_A_(.+?)_([A-Z])_([A-Z]+)_(\d+)$", folder_name)
    if m:
        return m.group(1), f"{m.group(2)}_{m.group(3)}_{m.group(4)}"
    return None, None


def _extract_complex_name_kplus(folder_name):
    """boltz_results_{COMPLEX}_cdr3_beta -> COMPLEX"""
    m = re.match(r"^boltz_results_(.+?)_cdr3_beta$", folder_name)
    return m.group(1) if m else None


def _find_prediction_files(boltz_dir):
    """Find model and confidence files inside a boltz_results_* directory."""
    predictions_dir = boltz_dir / "predictions"
    if not predictions_dir.exists():
        return [], []

    model_files = []
    confidence_files = []
    for subdir in sorted(predictions_dir.iterdir()):
        if not subdir.is_dir():
            continue
        for f in sorted(subdir.iterdir()):
            if f.name.startswith("confidence_") and f.suffix == ".json":
                confidence_files.append(f)
            elif f.suffix in (".pdb", ".cif") and "_model_" in f.name:
                model_files.append(f)
    return model_files, confidence_files


def discover_predictions(predictions_dir):
    """Scan predictions_dir for all methods and their prediction folders.

    Returns: {method: {complex_name: [prediction_info_dict, ...]}}
    """
    predictions_dir = Path(predictions_dir)
    results = {}

    for method, dirname in METHOD_DIRS.items():
        method_dir = predictions_dir / dirname
        if not method_dir.exists():
            log.warning("Method directory not found: %s", method_dir)
            results[method] = {}
            continue

        complexes = defaultdict(list)
        for boltz_dir in sorted(method_dir.iterdir()):
            if not boltz_dir.is_dir() or not boltz_dir.name.startswith("boltz_results"):
                continue

            model_files, confidence_files = _find_prediction_files(boltz_dir)
            if not model_files:
                log.debug("Skipping %s (no prediction files)", boltz_dir.name)
                continue

            config_label = "default"
            if method == "B1":
                complex_name = _extract_complex_name_b1(boltz_dir.name)
            elif method == "B2":
                complex_name, config_label = _extract_complex_name_b2(boltz_dir.name)
            elif method == "B3":
                complex_name, config_label = _extract_complex_name_b3(boltz_dir.name)
            elif method == "K+":
                complex_name = _extract_complex_name_kplus(boltz_dir.name)
            else:
                continue

            if complex_name is None:
                log.warning("Could not parse complex name from %s", boltz_dir.name)
                continue

            file_format = model_files[0].suffix  # .pdb or .cif
            complexes[complex_name].append({
                "folder_path": boltz_dir,
                "complex_name": complex_name,
                "config_label": config_label,
                "model_files": model_files,
                "confidence_files": confidence_files,
                "file_format": file_format,
            })

        results[method] = dict(complexes)
        n_complexes = len(complexes)
        n_configs = sum(len(v) for v in complexes.values())
        if n_complexes > 0:
            log.info("  %s: %d complex(es), %d config(s)", method, n_complexes, n_configs)
        else:
            log.info("  %s: no predictions found", method)

    return results


# ---------------------------------------------------------------------------
# 2. Model Selection
# ---------------------------------------------------------------------------


def _load_confidence(json_path):
    """Load a confidence JSON and return its dict."""
    with open(json_path) as f:
        return json.load(f)


def select_best_model(prediction_info):
    """Select the model with the highest confidence_score.

    Returns: (model_file_path, confidence_dict) or (None, None) if no files.
    """
    model_files = prediction_info["model_files"]
    confidence_files = prediction_info["confidence_files"]

    if not model_files:
        return None, None

    # Build a map from model index to files
    model_map = {}
    for mf in model_files:
        m = re.search(r"_model_(\d+)", mf.name)
        if m:
            model_map.setdefault(int(m.group(1)), {})["model"] = mf

    for cf in confidence_files:
        m = re.search(r"_model_(\d+)", cf.name)
        if m:
            idx = int(m.group(1))
            if idx in model_map:
                model_map[idx]["confidence_file"] = cf

    best_score = -1
    best_model = None
    best_conf = None

    for idx in sorted(model_map.keys()):
        entry = model_map[idx]
        if "confidence_file" not in entry:
            continue
        conf = _load_confidence(entry["confidence_file"])
        score = conf.get("confidence_score", 0)
        if score > best_score:
            best_score = score
            best_model = entry["model"]
            best_conf = conf

    if best_model is None and model_files:
        # Fallback: pick first model without confidence
        best_model = model_files[0]
        best_conf = {}

    return best_model, best_conf


# ---------------------------------------------------------------------------
# 3. DockQ Scoring
# ---------------------------------------------------------------------------


def _dockq_cache_path(output_dir, method, complex_name, config_label, model_name):
    cache_dir = Path(output_dir) / "dockq_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    return cache_dir / f"{method}_{complex_name}_{config_label}_{model_name}.json"


def run_dockq(model_path, native_path, output_dir, method, complex_name,
              config_label, skip_dockq=False):
    """Run DockQ and return parsed results dict.

    Uses caching to avoid re-running.
    """
    model_name = Path(model_path).stem
    cache_file = _dockq_cache_path(output_dir, method, complex_name, config_label, model_name)

    if cache_file.exists():
        log.debug("Using cached DockQ: %s", cache_file)
        with open(cache_file) as f:
            return json.load(f)

    if skip_dockq:
        log.warning("DockQ skipped (no cache) for %s %s %s", method, complex_name, config_label)
        return None

    dockq_json = cache_file.with_suffix(".dockq_out.json")
    cmd = [
        "conda", "run", "-n", "dockq2",
        "DockQ", str(model_path), str(native_path),
        "--json", str(dockq_json),
    ]

    log.info("Running DockQ: %s vs %s", Path(model_path).name, Path(native_path).name)
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
    except subprocess.TimeoutExpired:
        log.error("DockQ timed out for %s", model_path)
        return None
    except FileNotFoundError:
        log.error("DockQ command not found. Is dockq2 conda env available?")
        return None

    if proc.returncode != 0:
        log.error("DockQ failed (rc=%d): %s\n%s", proc.returncode, proc.stdout, proc.stderr)
        return None

    if not dockq_json.exists():
        log.error("DockQ JSON output not found at %s", dockq_json)
        return None

    with open(dockq_json) as f:
        raw = json.load(f)

    # Parse into a cleaner structure
    parsed = _parse_dockq_output(raw)

    # Cache the parsed result
    with open(cache_file, "w") as f:
        json.dump(parsed, f, indent=2)

    # Clean up the raw DockQ output
    dockq_json.unlink(missing_ok=True)

    return parsed


def _parse_dockq_output(raw):
    """Parse DockQ JSON output into a structured dict.

    DockQ v2 JSON format:
      { "model": ..., "native": ..., "best_dockq": float,
        "best_result": { "AB": {...}, "AC": {...}, "BC": {...} } }
    """
    result = {"interfaces": {}, "global_dockq": None}

    # Global DockQ (sum of per-interface DockQ values)
    best_dockq = raw.get("best_dockq")
    best_result = raw.get("best_result", {})

    # Per-interface results
    for interface_key, iface_data in best_result.items():
        if not isinstance(iface_data, dict) or "DockQ" not in iface_data:
            continue
        result["interfaces"][interface_key] = {
            "DockQ": iface_data.get("DockQ"),
            "iRMSD": iface_data.get("iRMSD"),
            "LRMSD": iface_data.get("LRMSD"),
            "fnat": iface_data.get("fnat"),
            "F1": iface_data.get("F1"),
            "fnonnat": iface_data.get("fnonnat"),
        }

    # Compute global DockQ as average across interfaces
    dockq_vals = [v["DockQ"] for v in result["interfaces"].values() if v.get("DockQ") is not None]
    result["global_dockq"] = sum(dockq_vals) / len(dockq_vals) if dockq_vals else None

    # Compute AB+AC average DockQ
    ab = result["interfaces"].get("AB", {}).get("DockQ")
    ac = result["interfaces"].get("AC", {}).get("DockQ")
    vals = [v for v in (ab, ac) if v is not None]
    result["ab_ac_avg_dockq"] = sum(vals) / len(vals) if vals else None

    return result


# ---------------------------------------------------------------------------
# 4. Epitope Extraction
# ---------------------------------------------------------------------------


def _load_structure(path):
    """Load a structure using BioPython PDBParser or MMCIFParser."""
    path = Path(path)
    if path.suffix == ".cif":
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)
    return parser.get_structure(path.stem, str(path))


def _is_heavy_atom(atom):
    """Check if atom is a heavy atom (not hydrogen/deuterium)."""
    elem = atom.element.strip().upper()
    if elem in ("H", "D"):
        return False
    # Fallback: some PDBs don't populate element field properly
    name = atom.get_name().strip()
    if not elem and name and name[0] in ("H", "D") and name[0:2] not in ("HE", "HG", "HO"):
        return False
    return True


def _get_chain_atoms(structure, chain_id, heavy_only=True):
    """Get atoms from a specific chain. Returns list of atoms."""
    atoms = []
    for model in structure:
        for chain in model:
            if chain.id == chain_id:
                for residue in chain:
                    # Skip hetero residues (water, ligands)
                    if residue.id[0] != " ":
                        continue
                    for atom in residue:
                        if heavy_only and not _is_heavy_atom(atom):
                            continue
                        atoms.append(atom)
        break  # Only first model
    return atoms


def _get_residues(structure, chain_id):
    """Get standard residues from a chain. Returns list of residues."""
    residues = []
    for model in structure:
        for chain in model:
            if chain.id == chain_id:
                for residue in chain:
                    if residue.id[0] != " ":
                        continue
                    residues.append(residue)
        break
    return residues


def extract_epitope(structure, cutoff=5.0):
    """Extract epitope residues: antigen residues within cutoff of any antibody heavy atom.

    Returns: set of (chain_id, resseq) tuples for epitope residues,
             and total number of antigen residues.
    """
    # Collect antibody heavy atom coordinates
    ab_atoms = []
    for chain_id in ANTIBODY_CHAINS:
        ab_atoms.extend(_get_chain_atoms(structure, chain_id, heavy_only=True))

    if not ab_atoms:
        log.warning("No antibody atoms found in structure")
        return set(), 0

    ab_coords = np.array([a.get_vector().get_array() for a in ab_atoms])

    # Build KDTree for antibody atoms
    if cKDTree is not None:
        tree = cKDTree(ab_coords)
    else:
        tree = None

    # Check each antigen residue
    antigen_residues = _get_residues(structure, ANTIGEN_CHAIN)
    epitope = set()

    for res in antigen_residues:
        res_id = (ANTIGEN_CHAIN, res.id[1])
        heavy_atoms = [a for a in res if _is_heavy_atom(a)]
        if not heavy_atoms:
            continue

        res_coords = np.array([a.get_vector().get_array() for a in heavy_atoms])

        if tree is not None:
            # Use KDTree for efficiency
            dists, _ = tree.query(res_coords)
            min_dist = dists.min()
        else:
            # Brute-force fallback
            diff = res_coords[:, np.newaxis, :] - ab_coords[np.newaxis, :, :]
            dists = np.sqrt((diff ** 2).sum(axis=2))
            min_dist = dists.min()

        if min_dist < cutoff:
            epitope.add(res_id)

    return epitope, len(antigen_residues)


# ---------------------------------------------------------------------------
# 5. Epitope Metrics
# ---------------------------------------------------------------------------


def compute_epitope_metrics(gt_epitope, pred_epitope, n_antigen_residues):
    """Compute epitope prediction metrics.

    Args:
        gt_epitope: set of (chain, resseq) ground truth epitope residues
        pred_epitope: set of (chain, resseq) predicted epitope residues
        n_antigen_residues: total number of antigen residues

    Returns: dict with Precision, Recall, F1, MCC
    """
    tp = len(gt_epitope & pred_epitope)
    fp = len(pred_epitope - gt_epitope)
    fn = len(gt_epitope - pred_epitope)
    tn = n_antigen_residues - tp - fp - fn
    tn = max(tn, 0)  # Safety

    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0.0

    denom = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    mcc = (tp * tn - fp * fn) / denom if denom > 0 else 0.0

    return {
        "TP": tp,
        "FP": fp,
        "FN": fn,
        "TN": tn,
        "Precision": precision,
        "Recall": recall,
        "F1": f1,
        "MCC": mcc,
        "n_gt_epitope": len(gt_epitope),
        "n_pred_epitope": len(pred_epitope),
        "n_antigen_residues": n_antigen_residues,
    }


# ---------------------------------------------------------------------------
# 6. Antibody-Aligned Antigen RMSD
# ---------------------------------------------------------------------------


def _get_ca_coords_ordered(structure, chain_id):
    """Get CA coordinates for a chain in sequential order.

    Returns: list of numpy arrays (one per residue with a CA atom)
    """
    coords = []
    for model in structure:
        for chain in model:
            if chain.id != chain_id:
                continue
            for residue in chain:
                if residue.id[0] != " ":
                    continue
                if "CA" in residue:
                    coords.append(residue["CA"].get_vector().get_array())
        break
    return coords


def compute_ab_aligned_ag_rmsd(pred_structure, gt_structure):
    """Compute antigen RMSD after superimposing on antibody CA atoms.

    Matches residues by sequential position within each chain (not by residue
    number), since GT and predictions may use different numbering schemes.

    1. Extract CA atoms from antibody chains (B, C) in both structures
    2. Superimpose prediction onto GT using antibody CAs (matched by position)
    3. Apply transformation to antigen CAs
    4. Compute RMSD of antigen CAs

    Returns: RMSD value (float) or None on failure
    """
    # Collect paired antibody CA coords by matching sequential positions
    pred_ab_coords_list = []
    gt_ab_coords_list = []
    for chain_id in ANTIBODY_CHAINS:
        pred_cas = _get_ca_coords_ordered(pred_structure, chain_id)
        gt_cas = _get_ca_coords_ordered(gt_structure, chain_id)
        n = min(len(pred_cas), len(gt_cas))
        for i in range(n):
            pred_ab_coords_list.append(pred_cas[i])
            gt_ab_coords_list.append(gt_cas[i])

    if len(pred_ab_coords_list) < 3:
        log.warning("Too few common antibody CA atoms (%d) for superposition",
                     len(pred_ab_coords_list))
        return None

    pred_ab_coords = np.array(pred_ab_coords_list)
    gt_ab_coords = np.array(gt_ab_coords_list)

    # Collect paired antigen CA coords
    pred_ag_cas = _get_ca_coords_ordered(pred_structure, ANTIGEN_CHAIN)
    gt_ag_cas = _get_ca_coords_ordered(gt_structure, ANTIGEN_CHAIN)
    n_ag = min(len(pred_ag_cas), len(gt_ag_cas))
    if n_ag < 1:
        log.warning("No common antigen CA atoms for RMSD")
        return None

    pred_ag_coords = np.array(pred_ag_cas[:n_ag])
    gt_ag_coords = np.array(gt_ag_cas[:n_ag])

    # Compute optimal rotation using Kabsch algorithm (SVD-based)
    # Center the coordinates
    pred_ab_center = pred_ab_coords.mean(axis=0)
    gt_ab_center = gt_ab_coords.mean(axis=0)
    pred_ab_centered = pred_ab_coords - pred_ab_center
    gt_ab_centered = gt_ab_coords - gt_ab_center

    # SVD of cross-covariance matrix
    H = pred_ab_centered.T @ gt_ab_centered
    U, S, Vt = np.linalg.svd(H)

    # Handle reflection case
    d = np.linalg.det(Vt.T @ U.T)
    sign_matrix = np.diag([1, 1, np.sign(d)])
    rot = Vt.T @ sign_matrix @ U.T

    # Apply transformation to antigen coords
    transformed_ag = (pred_ag_coords - pred_ab_center) @ rot.T + gt_ab_center

    # Compute RMSD
    diff = transformed_ag - gt_ag_coords
    rmsd = np.sqrt((diff ** 2).sum(axis=1).mean())
    return float(rmsd)


# ---------------------------------------------------------------------------
# 7. Confidence Extraction
# ---------------------------------------------------------------------------


def extract_confidence(confidence_dict):
    """Extract key confidence metrics from a confidence JSON dict."""
    if not confidence_dict:
        return {}
    return {
        "confidence_score": confidence_dict.get("confidence_score"),
        "iptm": confidence_dict.get("iptm"),
        "ptm": confidence_dict.get("ptm"),
        "complex_plddt": confidence_dict.get("complex_plddt"),
        "complex_iplddt": confidence_dict.get("complex_iplddt"),
    }


# ---------------------------------------------------------------------------
# 8. Result Aggregation
# ---------------------------------------------------------------------------


def aggregate_results(all_results):
    """Aggregate per-complex results into a DataFrame.

    Args:
        all_results: list of dicts, each with method/complex/metrics

    Returns: pd.DataFrame
    """
    if not all_results:
        return pd.DataFrame()

    df = pd.DataFrame(all_results)
    return df


def compute_summary(df):
    """Compute summary statistics grouped by method."""
    if df.empty:
        return pd.DataFrame()

    numeric_cols = [
        "DockQ_global", "DockQ_AB_AC", "iRMSD_AB_AC",
        "fnat_AB_AC", "ab_aligned_ag_rmsd",
        "Precision", "Recall", "F1", "MCC",
        "confidence_score", "iptm", "complex_plddt", "complex_iplddt",
    ]
    existing = [c for c in numeric_cols if c in df.columns]

    summary = df.groupby("method")[existing].agg(["mean", "median", "std", "count"])
    return summary


# ---------------------------------------------------------------------------
# 9. Output Generation
# ---------------------------------------------------------------------------


def _dockq_quality_label(dockq):
    """Classify DockQ score into quality category."""
    if dockq is None:
        return "N/A"
    for label, (lo, hi) in DOCKQ_QUALITY.items():
        if lo <= dockq < hi:
            return label
    return "N/A"


def print_summary_tables(df, summary):
    """Print formatted markdown tables to stdout."""
    if df.empty:
        print("\nNo results to display.")
        return

    methods_present = sorted(df["method"].unique())

    # Table 1: Epitope Prediction (Axis A)
    print("\n## Table 1: Epitope Prediction (Axis A)")
    print("| Method | Precision | Recall | F1 | MCC |")
    print("|--------|-----------|--------|----|-----|")
    for method in methods_present:
        mdf = df[df["method"] == method]
        p = mdf["Precision"].mean()
        r = mdf["Recall"].mean()
        f = mdf["F1"].mean()
        m = mdf["MCC"].mean()
        print(f"| {method:20s} | {p:.3f}     | {r:.3f}  | {f:.3f} | {m:.3f} |")

    # Table 2: Docking Quality (Axis B)
    print("\n## Table 2: Docking Quality (Axis B)")
    dockq_cols = ["DockQ_AB_AC", "DockQ_global", "iRMSD_AB_AC", "fnat_AB_AC", "ab_aligned_ag_rmsd"]
    avail_cols = [c for c in dockq_cols if c in df.columns]
    header = "| Method |" + "|".join(f" {c} " for c in avail_cols) + "|"
    sep = "|--------|" + "|".join("---" for _ in avail_cols) + "|"
    print(header)
    print(sep)
    for method in methods_present:
        mdf = df[df["method"] == method]
        vals = []
        for c in avail_cols:
            v = mdf[c].mean()
            if pd.notna(v):
                vals.append(f" {v:.3f} ")
            else:
                vals.append(" N/A ")
        print(f"| {method:20s} |" + "|".join(vals) + "|")

    # Table 3: Confidence Metrics
    print("\n## Table 3: Confidence Metrics")
    conf_cols = ["confidence_score", "iptm", "complex_plddt", "complex_iplddt"]
    avail_conf = [c for c in conf_cols if c in df.columns]
    header = "| Method |" + "|".join(f" {c} " for c in avail_conf) + "|"
    sep = "|--------|" + "|".join("---" for _ in avail_conf) + "|"
    print(header)
    print(sep)
    for method in methods_present:
        mdf = df[df["method"] == method]
        vals = []
        for c in avail_conf:
            v = mdf[c].mean()
            if pd.notna(v):
                vals.append(f" {v:.4f} ")
            else:
                vals.append(" N/A ")
        print(f"| {method:20s} |" + "|".join(vals) + "|")

    # Table 4: Per-Complex Breakdown
    print("\n## Table 4: Per-Complex Breakdown")
    print("| Complex | Method | DockQ (AB+AC) | DockQ (global) | F1 | Quality |")
    print("|---------|--------|---------------|----------------|-----|---------|")
    for _, row in df.sort_values(["complex", "method"]).iterrows():
        dq_abac = f"{row['DockQ_AB_AC']:.3f}" if pd.notna(row.get("DockQ_AB_AC")) else "N/A"
        dq_glob = f"{row['DockQ_global']:.3f}" if pd.notna(row.get("DockQ_global")) else "N/A"
        f1 = f"{row['F1']:.3f}" if pd.notna(row.get("F1")) else "N/A"
        qual = _dockq_quality_label(row.get("DockQ_global"))
        print(f"| {row['complex']:15s} | {row['method']:6s} | {dq_abac:>13s} | {dq_glob:>14s} | {f1:>5s} | {qual:9s} |")


def save_csv(df, summary, output_dir):
    """Save results and summary to CSV files."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    results_csv = output_dir / "results_per_complex.csv"
    df.to_csv(results_csv, index=False)
    log.info("Saved per-complex results to %s", results_csv)

    if not summary.empty:
        summary_csv = output_dir / "summary_by_method.csv"
        summary.to_csv(summary_csv)
        log.info("Saved summary to %s", summary_csv)


def generate_plots(df, output_dir):
    """Generate evaluation plots. Skips silently if matplotlib is unavailable."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        log.info("matplotlib not available, skipping plots")
        return

    output_dir = Path(output_dir) / "plots"
    output_dir.mkdir(parents=True, exist_ok=True)
    methods = sorted(df["method"].unique())

    # Plot 1: Mean DockQ per method
    if "DockQ_AB_AC" in df.columns:
        fig, ax = plt.subplots(figsize=(8, 5))
        means = [df[df["method"] == m]["DockQ_AB_AC"].mean() for m in methods]
        stds = [df[df["method"] == m]["DockQ_AB_AC"].std() for m in methods]
        stds = [s if not np.isnan(s) else 0 for s in stds]
        bars = ax.bar(methods, means, yerr=stds, capsize=5, color=["#4C72B0", "#DD8452", "#55A868", "#C44E52"][:len(methods)])
        ax.set_ylabel("DockQ (AB+AC avg)")
        ax.set_title("Docking Quality: DockQ (AB+AC interface)")
        ax.set_ylim(0, 1)
        for bar, val in zip(bars, means):
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.02,
                    f"{val:.3f}", ha="center", va="bottom", fontsize=10)
        fig.tight_layout()
        fig.savefig(output_dir / "dockq_bar.png", dpi=150)
        plt.close(fig)
        log.info("Saved dockq_bar.png")

    # Plot 2: Mean epitope F1 per method
    if "F1" in df.columns:
        fig, ax = plt.subplots(figsize=(8, 5))
        means = [df[df["method"] == m]["F1"].mean() for m in methods]
        stds = [df[df["method"] == m]["F1"].std() for m in methods]
        stds = [s if not np.isnan(s) else 0 for s in stds]
        bars = ax.bar(methods, means, yerr=stds, capsize=5, color=["#4C72B0", "#DD8452", "#55A868", "#C44E52"][:len(methods)])
        ax.set_ylabel("Epitope F1")
        ax.set_title("Epitope Prediction: F1 Score")
        ax.set_ylim(0, 1)
        for bar, val in zip(bars, means):
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.02,
                    f"{val:.3f}", ha="center", va="bottom", fontsize=10)
        fig.tight_layout()
        fig.savefig(output_dir / "epitope_f1_bar.png", dpi=150)
        plt.close(fig)
        log.info("Saved epitope_f1_bar.png")

    # Plot 3: K+ vs B1 DockQ scatter
    if "K+" in methods and "B1" in methods and "DockQ_AB_AC" in df.columns:
        b1_df = df[df["method"] == "B1"].set_index("complex")
        kp_df = df[df["method"] == "K+"].set_index("complex")
        common = sorted(set(b1_df.index) & set(kp_df.index))
        if common:
            fig, ax = plt.subplots(figsize=(6, 6))
            b1_vals = [b1_df.loc[c, "DockQ_AB_AC"] for c in common]
            kp_vals = [kp_df.loc[c, "DockQ_AB_AC"] for c in common]
            ax.scatter(b1_vals, kp_vals, s=80, zorder=5)
            for c, bv, kv in zip(common, b1_vals, kp_vals):
                ax.annotate(c, (bv, kv), textcoords="offset points",
                            xytext=(5, 5), fontsize=8)
            lims = [0, 1]
            ax.plot(lims, lims, "k--", alpha=0.3, label="y=x")
            ax.set_xlabel("B1 DockQ (AB+AC)")
            ax.set_ylabel("K+ DockQ (AB+AC)")
            ax.set_title("K+ vs B1: DockQ")
            ax.set_xlim(lims)
            ax.set_ylim(lims)
            ax.legend()
            fig.tight_layout()
            fig.savefig(output_dir / "kplus_vs_b1_dockq.png", dpi=150)
            plt.close(fig)
            log.info("Saved kplus_vs_b1_dockq.png")

    # Plot 4: DockQ quality distribution
    if "DockQ_global" in df.columns:
        fig, ax = plt.subplots(figsize=(8, 5))
        categories = list(DOCKQ_QUALITY.keys())
        x = np.arange(len(categories))
        width = 0.8 / max(len(methods), 1)
        colors = ["#4C72B0", "#DD8452", "#55A868", "#C44E52"]

        for i, method in enumerate(methods):
            mdf = df[df["method"] == method]
            counts = []
            for cat in categories:
                lo, hi = DOCKQ_QUALITY[cat]
                n = ((mdf["DockQ_global"] >= lo) & (mdf["DockQ_global"] < hi)).sum()
                counts.append(n)
            ax.bar(x + i * width, counts, width, label=method,
                   color=colors[i % len(colors)])

        ax.set_xticks(x + width * (len(methods) - 1) / 2)
        ax.set_xticklabels(categories)
        ax.set_ylabel("Count")
        ax.set_title("DockQ Quality Distribution")
        ax.legend()
        fig.tight_layout()
        fig.savefig(output_dir / "dockq_distribution.png", dpi=150)
        plt.close(fig)
        log.info("Saved dockq_distribution.png")


# ---------------------------------------------------------------------------
# 10. Main CLI
# ---------------------------------------------------------------------------


def evaluate_complex(method, complex_name, prediction_info, gt_path, output_dir,
                     cutoff, skip_dockq):
    """Evaluate a single method+complex+config. Returns a result dict."""
    # Select best model
    model_path, confidence_dict = select_best_model(prediction_info)
    if model_path is None:
        log.warning("No model file for %s / %s / %s", method, complex_name,
                     prediction_info["config_label"])
        return None

    log.info("Evaluating %s / %s [%s] -> %s",
             method, complex_name, prediction_info["config_label"], model_path.name)

    result = {
        "method": method,
        "complex": complex_name,
        "config_label": prediction_info["config_label"],
        "model_file": str(model_path.name),
    }

    # DockQ
    dockq = run_dockq(model_path, gt_path, output_dir, method, complex_name,
                      prediction_info["config_label"], skip_dockq=skip_dockq)
    if dockq:
        result["DockQ_global"] = dockq.get("global_dockq")
        result["DockQ_AB_AC"] = dockq.get("ab_ac_avg_dockq")
        # Per-interface
        for iface_id in ("AB", "AC", "BC"):
            iface = dockq.get("interfaces", {}).get(iface_id, {})
            if iface:
                result[f"DockQ_{iface_id}"] = iface.get("DockQ")
                result[f"iRMSD_{iface_id}"] = iface.get("iRMSD")
                result[f"fnat_{iface_id}"] = iface.get("fnat")
        # Average iRMSD for AB+AC
        irmsd_ab = result.get("iRMSD_AB")
        irmsd_ac = result.get("iRMSD_AC")
        vals = [v for v in (irmsd_ab, irmsd_ac) if v is not None]
        result["iRMSD_AB_AC"] = sum(vals) / len(vals) if vals else None
        # Average fnat for AB+AC
        fnat_ab = result.get("fnat_AB")
        fnat_ac = result.get("fnat_AC")
        vals = [v for v in (fnat_ab, fnat_ac) if v is not None]
        result["fnat_AB_AC"] = sum(vals) / len(vals) if vals else None
    else:
        result["DockQ_global"] = None
        result["DockQ_AB_AC"] = None

    # Load structures for epitope + RMSD
    try:
        pred_structure = _load_structure(model_path)
        gt_structure = _load_structure(gt_path)
    except Exception as e:
        log.error("Failed to load structures: %s", e)
        return result

    # Epitope extraction
    gt_epitope, n_ag = extract_epitope(gt_structure, cutoff=cutoff)
    pred_epitope, _ = extract_epitope(pred_structure, cutoff=cutoff)

    epitope_metrics = compute_epitope_metrics(gt_epitope, pred_epitope, n_ag)
    result.update(epitope_metrics)

    # Antibody-aligned antigen RMSD
    ag_rmsd = compute_ab_aligned_ag_rmsd(pred_structure, gt_structure)
    result["ab_aligned_ag_rmsd"] = ag_rmsd

    # Confidence metrics
    conf = extract_confidence(confidence_dict)
    result.update(conf)

    return result


def main():
    parser = argparse.ArgumentParser(
        description="Evaluate K+ region-specific beta-scaling against baselines")
    parser.add_argument("--predictions_dir", default="predictions_examples",
                        help="Directory containing method subdirectories (default: predictions_examples)")
    parser.add_argument("--gt_dir", default="pdb_minimized",
                        help="Directory containing ground-truth PDB files (default: pdb_minimized)")
    parser.add_argument("--output_dir", default="evaluation_output",
                        help="Output directory for results (default: evaluation_output)")
    parser.add_argument("--cutoff", type=float, default=5.0,
                        help="Distance cutoff in Angstroms for epitope definition (default: 5.0)")
    parser.add_argument("--skip_dockq", action="store_true",
                        help="Skip DockQ computation (use cached results only)")
    parser.add_argument("--skip_plots", action="store_true",
                        help="Skip plot generation")
    parser.add_argument("--verbose", "-v", action="store_true",
                        help="Enable debug logging")
    args = parser.parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    predictions_dir = Path(args.predictions_dir)
    gt_dir = Path(args.gt_dir)
    output_dir = Path(args.output_dir)

    if not predictions_dir.exists():
        log.error("Predictions directory not found: %s", predictions_dir)
        sys.exit(1)
    if not gt_dir.exists():
        log.error("Ground truth directory not found: %s", gt_dir)
        sys.exit(1)

    output_dir.mkdir(parents=True, exist_ok=True)

    # Build GT lookup: complex_name -> pdb_path
    gt_files = {}
    for f in sorted(gt_dir.glob("*.pdb")):
        gt_files[f.stem] = f
    log.info("Found %d ground-truth structures in %s", len(gt_files), gt_dir)

    # Discover predictions
    log.info("Discovering predictions in %s ...", predictions_dir)
    predictions = discover_predictions(predictions_dir)

    # Evaluate
    all_results = []
    for method, complexes in predictions.items():
        for complex_name, configs in complexes.items():
            if complex_name not in gt_files:
                log.warning("No ground truth for %s (method %s), skipping", complex_name, method)
                continue

            gt_path = gt_files[complex_name]

            for pred_info in configs:
                result = evaluate_complex(
                    method, complex_name, pred_info, gt_path, output_dir,
                    args.cutoff, args.skip_dockq,
                )
                if result:
                    all_results.append(result)

    # Aggregate
    df = aggregate_results(all_results)
    if df.empty:
        log.warning("No results collected. Check your predictions directory.")
        sys.exit(0)

    summary = compute_summary(df)

    # Output
    print_summary_tables(df, summary)
    save_csv(df, summary, output_dir)

    if not args.skip_plots:
        generate_plots(df, output_dir)

    log.info("Evaluation complete. Results in %s", output_dir)


if __name__ == "__main__":
    main()
