#!/usr/bin/env python
"""
Evaluation script for Hierarchical Steering (Idea Y).

Compares hierarchical steering against 3 baselines across 48 antibody-antigen complexes.

Metrics computed:
  1. DockQ (total, Ab-Ag only, per-interface) via DockQ v2 Python API
  2. CDR loop RMSD (H1, H2, H3, L1, L2, L3) after framework alignment
  3. Epitope prediction (precision, recall, F1, MCC at 5A and 8A)
  4. Boltz confidence metrics (confidence_score, iptm, ptm, iplddt)

Usage:
  /home/oc/anaconda3/envs/dockq2/bin/python scripts/eval/evaluate_hierarchical_steering.py \
      --ground_truth_dir /mnt/c/Users/octav/Documents/claude/proteinEBM/boltz/pdb_minimized \
      --predictions_dir predictions_examples \
      --cdrs_csv examples/cdrs.csv \
      --out_dir evaluation_results

Requires: conda environment dockq2 (DockQ 2.1.3, biopython, numpy, pandas, scipy, matplotlib)
"""

import argparse
import ast
import csv
import json
import math
import os
import sys
import warnings
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")

# ---------------------------------------------------------------------------
# DockQ import
# ---------------------------------------------------------------------------
from DockQ.DockQ import load_PDB, run_on_all_native_interfaces

# ---------------------------------------------------------------------------
# BioPython imports
# ---------------------------------------------------------------------------
from Bio.PDB import PDBParser, MMCIFParser, Superimposer
from Bio.SeqUtils import seq1
from Bio.Align import PairwiseAligner


# ============================================================================
# CONFIGURATION
# ============================================================================

# Method definitions: label -> (folder_pattern, naming_convention)
# Each method specifies how to discover prediction folders and map to complexes
METHOD_CONFIGS = {
    "baseline": {
        "subfolder": "antigen_cut",
        "prefix": "boltz_results_",
        "suffix": "",
        "multi_setting": False,
    },
    "contact_restraints": {
        "subfolder": "antigen_cut_contact_restraints",
        "prefix": "boltz_results_restraint_",
        "suffix": "",
        "multi_setting": True,
    },
    "pocket_guided": {
        "subfolder": "antigen_cut_vhvl_msa_pocket_ab_boltz_post_2023_omm",
        "prefix": "boltz_results_restraint_to_A_",
        "suffix": "",
        "multi_setting": True,
    },
    "hierarchical": {
        "subfolder": "new_feature_hierarchical",
        "prefix": "boltz_results_",
        "suffix": "_hierarchical",
        "multi_setting": False,
    },
}

GROUND_TRUTH_DIR = None  # Set by CLI
CHAIN_MAP = {"A": "A", "B": "B", "C": "C"}

# Antigen chain, heavy chain, light chain
AG_CHAIN = "A"
HEAVY_CHAIN = "B"
LIGHT_CHAIN = "C"


# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

def load_cdrs(cdrs_csv):
    """Load CDR definitions from CSV. Returns dict: complex_name -> cdr_dict."""
    cdrs = {}
    with open(cdrs_csv) as f:
        reader = csv.DictReader(f)
        for row in reader:
            name = row["complex"]
            cdr_dict = {}
            for key in ["cdr1_h", "cdr2_h", "cdr3_h", "cdr1_l", "cdr2_l", "cdr3_l"]:
                cdr_dict[key] = ast.literal_eval(row[key])
            cdrs[name] = cdr_dict
    return cdrs


def parse_structure(filepath):
    """Parse PDB or CIF structure file with BioPython."""
    filepath = str(filepath)
    if filepath.endswith(".cif"):
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)
    return parser.get_structure("s", filepath)


def get_chain_residues(structure, chain_id):
    """Get standard residues (excluding heteroatoms/water) for a chain."""
    try:
        chain = structure[0][chain_id]
    except KeyError:
        return []
    return [r for r in chain.get_residues() if r.get_id()[0] == " "]


def get_chain_sequence(residues):
    """Get one-letter amino acid sequence from residue list."""
    return "".join(seq1(r.get_resname()) for r in residues)


def build_residue_mapping(gt_residues, pred_residues):
    """Build a mapping from pred residue list index to GT residue list index
    using pairwise sequence alignment.

    Returns dict: pred_list_idx (0-based) -> gt_list_idx (0-based)
    """
    gt_seq = get_chain_sequence(gt_residues)
    pred_seq = get_chain_sequence(pred_residues)

    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -0.5

    alignments = list(aligner.align(gt_seq, pred_seq))
    if not alignments:
        return {}

    a = alignments[0]
    mapping = {}
    for (gt_start, gt_end), (pred_start, pred_end) in zip(a.aligned[0], a.aligned[1]):
        for gi, pi in zip(range(gt_start, gt_end), range(pred_start, pred_end)):
            mapping[pi] = gi

    return mapping


# ============================================================================
# DISCOVERY: Find prediction folders for each method × complex
# ============================================================================

def discover_predictions(predictions_dir, complex_names):
    """Discover available prediction folders for each method and complex.

    Returns: dict[method][complex_name] -> list of prediction folder paths
    """
    results = {}
    for method, config in METHOD_CONFIGS.items():
        method_dir = Path(predictions_dir) / config["subfolder"]
        if not method_dir.exists():
            print(f"  Warning: method directory not found: {method_dir}")
            results[method] = {}
            continue

        method_results = defaultdict(list)
        for folder in sorted(method_dir.iterdir()):
            if not folder.is_dir():
                continue
            folder_name = folder.name

            # Try to match this folder to a complex
            for cname in complex_names:
                if config["multi_setting"]:
                    # Multi-setting: folder contains complex name somewhere
                    if cname in folder_name:
                        method_results[cname].append(folder)
                        break
                else:
                    # Single setting: exact match with prefix/suffix
                    expected = f"{config['prefix']}{cname}{config['suffix']}"
                    if folder_name == expected:
                        method_results[cname].append(folder)
                        break

        results[method] = dict(method_results)
        n_complexes = len(method_results)
        n_total = sum(len(v) for v in method_results.values())
        print(f"  {method}: {n_complexes} complexes, {n_total} total folders")

    return results


def find_prediction_files(run_folder):
    """Given a run folder (boltz_results_XXX/), find the predictions subfolder
    and return list of (model_idx, structure_path, confidence_path) tuples.
    """
    pred_dir = run_folder / "predictions"
    if not pred_dir.exists():
        return []

    # The predictions folder contains one subfolder named after the run
    subfolders = [d for d in pred_dir.iterdir() if d.is_dir()]
    if not subfolders:
        return []

    inner = subfolders[0]
    models = []
    for model_idx in range(10):  # check up to 10 models
        # Try both .pdb and .cif
        struct_path = None
        for ext in [".pdb", ".cif"]:
            candidate = inner / f"{inner.name}_model_{model_idx}{ext}"
            if candidate.exists():
                struct_path = candidate
                break
        if struct_path is None:
            continue

        conf_path = inner / f"confidence_{inner.name}_model_{model_idx}.json"
        if not conf_path.exists():
            conf_path = None

        models.append((model_idx, struct_path, conf_path))

    return models


# ============================================================================
# METRIC 1: DockQ
# ============================================================================

def compute_dockq(model_path, native_path, chain_map=None):
    """Compute DockQ metrics using DockQ v2 Python API.

    Returns dict with:
      - total_dockq: average DockQ across all interfaces
      - ab_ag_dockq: average DockQ for Ab-Ag interfaces only (AB + AC)
      - per_interface: dict of interface -> {DockQ, iRMSD, LRMSD, fnat, F1, fnonnat}
    """
    if chain_map is None:
        chain_map = CHAIN_MAP

    try:
        model = load_PDB(str(model_path))
        native = load_PDB(str(native_path))
        result_dict, global_dockq = run_on_all_native_interfaces(
            model, native, chain_map=chain_map
        )
    except Exception as e:
        print(f"    DockQ error for {model_path}: {e}")
        return None

    if not result_dict:
        return None

    per_interface = {}
    dockq_values = []
    ab_ag_dockqs = []

    for iface, data in result_dict.items():
        per_interface[iface] = {
            "DockQ": data["DockQ"],
            "iRMSD": data["iRMSD"],
            "LRMSD": data["LRMSD"],
            "fnat": data["fnat"],
            "F1": data["F1"],
            "fnonnat": data.get("fnonnat", 0),
            "clashes": data.get("clashes", 0),
        }
        dockq_values.append(data["DockQ"])
        # Ab-Ag interfaces: those involving antigen chain A
        if AG_CHAIN in iface:
            ab_ag_dockqs.append(data["DockQ"])

    total_dockq = np.mean(dockq_values) if dockq_values else 0.0
    ab_ag_dockq = np.mean(ab_ag_dockqs) if ab_ag_dockqs else 0.0

    return {
        "total_dockq": total_dockq,
        "ab_ag_dockq": ab_ag_dockq,
        "per_interface": per_interface,
    }


# ============================================================================
# METRIC 2: CDR RMSD
# ============================================================================

def compute_cdr_rmsd(model_path, native_path, cdr_indices):
    """Compute per-CDR Cα RMSD after framework alignment.

    cdr_indices: dict with keys cdr1_h,...,cdr3_l; values are 1-indexed positions
    Returns dict: cdr_name -> RMSD in Angstroms, plus 'framework_rmsd'
    """
    try:
        model_struct = parse_structure(model_path)
        native_struct = parse_structure(native_path)
    except Exception as e:
        print(f"    Parse error: {e}")
        return None

    # Build sequence alignment mappings for heavy and light chains
    mappings = {}
    for chain_id in [HEAVY_CHAIN, LIGHT_CHAIN]:
        gt_res = get_chain_residues(native_struct, chain_id)
        pred_res = get_chain_residues(model_struct, chain_id)
        if not gt_res or not pred_res:
            return None
        mappings[chain_id] = (gt_res, pred_res, build_residue_mapping(gt_res, pred_res))

    # Collect all CDR positions (0-indexed into pred residue list)
    all_cdr_positions = {HEAVY_CHAIN: set(), LIGHT_CHAIN: set()}
    cdr_chain_map = {
        "cdr1_h": HEAVY_CHAIN, "cdr2_h": HEAVY_CHAIN, "cdr3_h": HEAVY_CHAIN,
        "cdr1_l": LIGHT_CHAIN, "cdr2_l": LIGHT_CHAIN, "cdr3_l": LIGHT_CHAIN,
    }
    for cdr_name, chain_id in cdr_chain_map.items():
        for pos in cdr_indices[cdr_name]:
            all_cdr_positions[chain_id].add(pos - 1)  # Convert 1-indexed to 0-indexed

    # Collect framework Cα atoms (non-CDR antibody residues)
    framework_model_atoms = []
    framework_native_atoms = []
    for chain_id in [HEAVY_CHAIN, LIGHT_CHAIN]:
        gt_res, pred_res, mapping = mappings[chain_id]
        cdr_pos_set = all_cdr_positions[chain_id]
        for pred_idx, gt_idx in mapping.items():
            if pred_idx in cdr_pos_set:
                continue
            pred_r = pred_res[pred_idx]
            gt_r = gt_res[gt_idx]
            if "CA" in pred_r and "CA" in gt_r:
                framework_model_atoms.append(pred_r["CA"])
                framework_native_atoms.append(gt_r["CA"])

    if len(framework_model_atoms) < 10:
        print(f"    Too few framework atoms: {len(framework_model_atoms)}")
        return None

    # Superimpose on framework
    sup = Superimposer()
    sup.set_atoms(framework_native_atoms, framework_model_atoms)
    sup.apply(list(model_struct.get_atoms()))

    framework_rmsd = sup.rms

    # Compute per-CDR RMSD
    results = {"framework_rmsd": framework_rmsd}
    for cdr_name, chain_id in cdr_chain_map.items():
        gt_res, pred_res, mapping = mappings[chain_id]
        positions = [p - 1 for p in cdr_indices[cdr_name]]  # 0-indexed
        diffs_sq = []
        for pred_idx in positions:
            gt_idx = mapping.get(pred_idx)
            if gt_idx is None:
                continue
            pred_r = pred_res[pred_idx]
            gt_r = gt_res[gt_idx]
            if "CA" in pred_r and "CA" in gt_r:
                m_coord = pred_r["CA"].get_vector().get_array()
                n_coord = gt_r["CA"].get_vector().get_array()
                diffs_sq.append(np.sum((m_coord - n_coord) ** 2))
        if diffs_sq:
            results[cdr_name] = float(np.sqrt(np.mean(diffs_sq)))

    return results


# ============================================================================
# METRIC 3: Epitope prediction
# ============================================================================

def get_contact_residues(structure, chain1_id, chain2_id, threshold=5.0):
    """Get residues on chain1 that contact chain2 (heavy atom distance < threshold).

    Returns set of residue list indices (0-indexed) on chain1 that are in contact.
    """
    chain1_res = get_chain_residues(structure, chain1_id)
    chain2_res = get_chain_residues(structure, chain2_id)
    if not chain1_res or not chain2_res:
        return set(), len(chain1_res)

    # Build KD-tree approach: collect all heavy atoms from chain2
    chain2_coords = []
    for res in chain2_res:
        for atom in res.get_atoms():
            if atom.element != "H":
                chain2_coords.append(atom.get_vector().get_array())
    if not chain2_coords:
        return set(), len(chain1_res)
    chain2_coords = np.array(chain2_coords)

    contact_indices = set()
    for idx, res in enumerate(chain1_res):
        for atom in res.get_atoms():
            if atom.element == "H":
                continue
            coord = atom.get_vector().get_array()
            dists = np.sqrt(np.sum((chain2_coords - coord) ** 2, axis=1))
            if np.min(dists) < threshold:
                contact_indices.add(idx)
                break

    return contact_indices, len(chain1_res)


def get_ca_contact_residues(structure, chain1_id, chain2_id, threshold=8.0):
    """Get residues on chain1 that contact chain2 (CA-CA distance < threshold).

    Returns set of residue list indices (0-indexed) on chain1.
    """
    chain1_res = get_chain_residues(structure, chain1_id)
    chain2_res = get_chain_residues(structure, chain2_id)
    if not chain1_res or not chain2_res:
        return set(), len(chain1_res)

    # Collect CA coords from chain2
    chain2_ca = []
    for res in chain2_res:
        if "CA" in res:
            chain2_ca.append(res["CA"].get_vector().get_array())
    if not chain2_ca:
        return set(), len(chain1_res)
    chain2_ca = np.array(chain2_ca)

    contact_indices = set()
    for idx, res in enumerate(chain1_res):
        if "CA" not in res:
            continue
        coord = res["CA"].get_vector().get_array()
        dists = np.sqrt(np.sum((chain2_ca - coord) ** 2, axis=1))
        if np.min(dists) < threshold:
            contact_indices.add(idx)

    return contact_indices, len(chain1_res)


def compute_epitope_metrics(model_path, native_path):
    """Compute epitope prediction metrics.

    Epitope = antigen residues contacting the antibody.
    Compares model vs native epitope at 5A (heavy atom) and 8A (CA-CA).
    """
    try:
        model_struct = parse_structure(model_path)
        native_struct = parse_structure(native_path)
    except Exception as e:
        print(f"    Parse error: {e}")
        return None

    results = {}
    for label, contact_fn, threshold in [
        ("5A", get_contact_residues, 5.0),
        ("8A", get_ca_contact_residues, 8.0),
    ]:
        # Get antigen epitope residues: antigen residues contacting either Ab chain
        native_epitope_h, n_ag_native = contact_fn(native_struct, AG_CHAIN, HEAVY_CHAIN, threshold)
        native_epitope_l, _ = contact_fn(native_struct, AG_CHAIN, LIGHT_CHAIN, threshold)
        native_epitope = native_epitope_h | native_epitope_l

        model_epitope_h, n_ag_model = contact_fn(model_struct, AG_CHAIN, HEAVY_CHAIN, threshold)
        model_epitope_l, _ = contact_fn(model_struct, AG_CHAIN, LIGHT_CHAIN, threshold)
        model_epitope = model_epitope_h | model_epitope_l

        # Note: native and model antigen may differ in length due to structure differences
        # Use the set of all possible residue indices as the universe
        n_ag = max(n_ag_native, n_ag_model) if n_ag_native > 0 else n_ag_model

        tp = len(native_epitope & model_epitope)
        fp = len(model_epitope - native_epitope)
        fn = len(native_epitope - model_epitope)
        tn = n_ag - tp - fp - fn
        tn = max(tn, 0)

        precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
        recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
        f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0.0

        # MCC
        denom = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
        mcc = (tp * tn - fp * fn) / denom if denom > 0 else 0.0

        results[f"epitope_precision_{label}"] = precision
        results[f"epitope_recall_{label}"] = recall
        results[f"epitope_f1_{label}"] = f1
        results[f"epitope_mcc_{label}"] = mcc
        results[f"n_true_epitope_{label}"] = len(native_epitope)
        results[f"n_pred_epitope_{label}"] = len(model_epitope)

    return results


# ============================================================================
# METRIC 4: Confidence metrics
# ============================================================================

def load_confidence(conf_path):
    """Load confidence metrics from JSON file."""
    if conf_path is None or not Path(conf_path).exists():
        return None
    with open(conf_path) as f:
        data = json.load(f)
    return {
        "confidence_score": data.get("confidence_score", 0),
        "iptm": data.get("iptm", 0),
        "ptm": data.get("ptm", 0),
        "complex_plddt": data.get("complex_plddt", 0),
        "complex_iplddt": data.get("complex_iplddt", 0),
    }


# ============================================================================
# MODEL SELECTION
# ============================================================================

def select_models(models):
    """Select top-1 model by confidence_score.

    Args:
        models: list of (model_idx, struct_path, conf_path)

    Returns: (top1_idx, confidences_list) where top1_idx indexes into models list
    """
    confidences = []
    for model_idx, struct_path, conf_path in models:
        conf = load_confidence(conf_path)
        confidences.append(conf)

    # Select top-1 by confidence_score
    scores = [c["confidence_score"] if c else -1 for c in confidences]
    top1_idx = int(np.argmax(scores))

    return top1_idx, confidences


# ============================================================================
# MAIN EVALUATION LOOP
# ============================================================================

def evaluate_complex_method(complex_name, run_folders, native_path, cdr_indices):
    """Evaluate all models for a complex under one method (possibly multiple settings).

    Returns list of dicts, one per setting, each containing:
      - All metrics for top-1 and oracle model
      - Per-model DockQ for oracle/average computation
    """
    settings_results = []

    for run_folder in run_folders:
        models = find_prediction_files(run_folder)
        if not models:
            continue

        # Model selection
        top1_idx, confidences = select_models(models)

        # Evaluate all models for DockQ (needed for oracle and average)
        all_dockq = []
        for i, (model_idx, struct_path, conf_path) in enumerate(models):
            dq = compute_dockq(struct_path, native_path)
            all_dockq.append(dq)

        # Top-1 metrics
        top1_model_idx, top1_struct, top1_conf = models[top1_idx]
        top1_dockq = all_dockq[top1_idx]

        # Oracle: best total_dockq
        valid_dockq = [(i, dq) for i, dq in enumerate(all_dockq) if dq is not None]
        if valid_dockq:
            oracle_idx = max(valid_dockq, key=lambda x: x[1]["total_dockq"])[0]
            oracle_dockq = all_dockq[oracle_idx]
        else:
            oracle_idx = top1_idx
            oracle_dockq = top1_dockq

        # Average DockQ
        valid_totals = [dq["total_dockq"] for dq in all_dockq if dq is not None]
        valid_ab_ag = [dq["ab_ag_dockq"] for dq in all_dockq if dq is not None]
        avg_dockq = np.mean(valid_totals) if valid_totals else None
        avg_ab_ag_dockq = np.mean(valid_ab_ag) if valid_ab_ag else None

        # CDR RMSD (top-1 model only)
        cdr_rmsd = compute_cdr_rmsd(top1_struct, native_path, cdr_indices)

        # Epitope metrics (top-1 model only)
        epitope = compute_epitope_metrics(top1_struct, native_path)

        # Compile results
        result = {
            "complex": complex_name,
            "run_folder": str(run_folder.name),
            "n_models": len(models),
            "top1_model_idx": top1_model_idx,
            "oracle_model_idx": models[oracle_idx][0] if valid_dockq else None,
        }

        # DockQ metrics
        if top1_dockq:
            result["top1_total_dockq"] = top1_dockq["total_dockq"]
            result["top1_ab_ag_dockq"] = top1_dockq["ab_ag_dockq"]
            for iface, idata in top1_dockq["per_interface"].items():
                for metric, val in idata.items():
                    result[f"top1_{iface}_{metric}"] = val
        if oracle_dockq:
            result["oracle_total_dockq"] = oracle_dockq["total_dockq"]
            result["oracle_ab_ag_dockq"] = oracle_dockq["ab_ag_dockq"]
        if avg_dockq is not None:
            result["avg_total_dockq"] = avg_dockq
            result["avg_ab_ag_dockq"] = avg_ab_ag_dockq

        # CDR RMSD metrics
        if cdr_rmsd:
            for key, val in cdr_rmsd.items():
                result[f"cdr_{key}"] = val

        # Epitope metrics
        if epitope:
            result.update(epitope)

        # Confidence metrics (top-1)
        if confidences[top1_idx]:
            for key, val in confidences[top1_idx].items():
                result[f"conf_{key}"] = val

        settings_results.append(result)

    return settings_results


def aggregate_multi_setting(settings_results, aggregation="best"):
    """For methods with multiple settings per complex, aggregate to one result.

    aggregation: 'best' (oracle over settings) or 'average'
    """
    if not settings_results:
        return None
    if len(settings_results) == 1:
        return settings_results[0]

    if aggregation == "best":
        # Pick setting with best top1_total_dockq
        valid = [r for r in settings_results if "top1_total_dockq" in r]
        if not valid:
            return settings_results[0]
        return max(valid, key=lambda r: r["top1_total_dockq"])
    else:
        # Average across settings
        result = dict(settings_results[0])
        numeric_keys = [k for k, v in result.items() if isinstance(v, (int, float))]
        for key in numeric_keys:
            vals = [r[key] for r in settings_results if key in r and r[key] is not None]
            if vals:
                result[key] = np.mean(vals)
        result["n_settings"] = len(settings_results)
        return result


# ============================================================================
# STATISTICAL ANALYSIS
# ============================================================================

def compute_statistics(df, metric, method1="hierarchical", method2="baseline"):
    """Compute paired Wilcoxon signed-rank test between two methods on a metric."""
    df1 = df[df["method"] == method1][["complex", metric]].dropna()
    df2 = df[df["method"] == method2][["complex", metric]].dropna()
    merged = pd.merge(df1, df2, on="complex", suffixes=("_1", "_2"))
    if len(merged) < 5:
        return {"n": len(merged), "p_value": None, "statistic": None}

    vals1 = merged[f"{metric}_1"].values
    vals2 = merged[f"{metric}_2"].values
    diff = vals1 - vals2

    # Count wins/ties/losses
    wins = int(np.sum(diff > 0.001))
    losses = int(np.sum(diff < -0.001))
    ties = len(diff) - wins - losses

    try:
        stat, p_val = stats.wilcoxon(vals1, vals2, alternative="two-sided")
    except Exception:
        stat, p_val = None, None

    return {
        "n": len(merged),
        "mean_diff": float(np.mean(diff)),
        "median_diff": float(np.median(diff)),
        "wins": wins,
        "losses": losses,
        "ties": ties,
        "statistic": stat,
        "p_value": p_val,
    }


def bootstrap_ci(values, n_boot=1000, alpha=0.05):
    """Compute mean and 95% bootstrap CI."""
    values = np.array(values)
    values = values[~np.isnan(values)]
    if len(values) == 0:
        return np.nan, np.nan, np.nan
    boot_means = [np.random.choice(values, len(values), replace=True).mean()
                  for _ in range(n_boot)]
    return (
        float(np.mean(values)),
        float(np.percentile(boot_means, 100 * alpha / 2)),
        float(np.percentile(boot_means, 100 * (1 - alpha / 2))),
    )


# ============================================================================
# PLOTTING
# ============================================================================

def generate_plots(df, out_dir):
    """Generate evaluation plots."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("  matplotlib not available, skipping plots")
        return

    figures_dir = Path(out_dir) / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)

    methods = [m for m in ["baseline", "contact_restraints", "pocket_guided", "hierarchical"]
               if m in df["method"].unique()]
    colors = {"baseline": "#888888", "contact_restraints": "#E88B23",
              "pocket_guided": "#931652", "hierarchical": "#188F52"}

    # --- Plot 1: DockQ bar plot with bootstrap CIs ---
    fig, axes = plt.subplots(1, 3, figsize=(14, 5))
    for ax, metric, title in zip(
        axes,
        ["top1_total_dockq", "top1_ab_ag_dockq", "oracle_total_dockq"],
        ["Top-1 Total DockQ", "Top-1 Ab-Ag DockQ", "Oracle Total DockQ"],
    ):
        means, lowers, uppers = [], [], []
        for method in methods:
            vals = df[df["method"] == method][metric].dropna().values
            m, lo, hi = bootstrap_ci(vals)
            means.append(m)
            lowers.append(m - lo)
            uppers.append(hi - m)
        x = np.arange(len(methods))
        bars = ax.bar(x, means, yerr=[lowers, uppers], capsize=4,
                      color=[colors.get(m, "#555") for m in methods])
        ax.set_xticks(x)
        ax.set_xticklabels([m.replace("_", "\n") for m in methods], fontsize=8)
        ax.set_title(title)
        ax.set_ylim(0, 1)
    plt.tight_layout()
    plt.savefig(figures_dir / "dockq_barplot.png", dpi=150)
    plt.close()

    # --- Plot 2: DockQ boxplot ---
    fig, ax = plt.subplots(figsize=(8, 5))
    data_to_plot = [df[df["method"] == m]["top1_ab_ag_dockq"].dropna().values for m in methods]
    bp = ax.boxplot(data_to_plot, labels=[m.replace("_", "\n") for m in methods],
                    patch_artist=True)
    for patch, method in zip(bp["boxes"], methods):
        patch.set_facecolor(colors.get(method, "#555"))
        patch.set_alpha(0.7)
    ax.set_ylabel("Ab-Ag DockQ (top-1)")
    ax.set_title("DockQ Distribution Across Complexes")
    plt.tight_layout()
    plt.savefig(figures_dir / "dockq_boxplot.png", dpi=150)
    plt.close()

    # --- Plot 3: Scatter - hierarchical vs baseline ---
    if "hierarchical" in df["method"].unique() and "baseline" in df["method"].unique():
        df_h = df[df["method"] == "hierarchical"][["complex", "top1_ab_ag_dockq"]].dropna()
        df_b = df[df["method"] == "baseline"][["complex", "top1_ab_ag_dockq"]].dropna()
        merged = pd.merge(df_h, df_b, on="complex", suffixes=("_hier", "_base"))
        if len(merged) > 0:
            fig, ax = plt.subplots(figsize=(6, 6))
            ax.scatter(merged["top1_ab_ag_dockq_base"], merged["top1_ab_ag_dockq_hier"],
                       c=colors["hierarchical"], alpha=0.7, edgecolors="black", linewidths=0.5)
            ax.plot([0, 1], [0, 1], "k--", alpha=0.5)
            ax.set_xlabel("Baseline Ab-Ag DockQ")
            ax.set_ylabel("Hierarchical Ab-Ag DockQ")
            ax.set_title("Hierarchical vs Baseline (per complex)")
            ax.set_xlim(0, 1)
            ax.set_ylim(0, 1)
            ax.set_aspect("equal")
            plt.tight_layout()
            plt.savefig(figures_dir / "dockq_scatter_hierarchical_vs_baseline.png", dpi=150)
            plt.close()

    # --- Plot 4: CDR RMSD bar plot ---
    cdr_cols = ["cdr_cdr1_h", "cdr_cdr2_h", "cdr_cdr3_h",
                "cdr_cdr1_l", "cdr_cdr2_l", "cdr_cdr3_l"]
    cdr_labels = ["CDR-H1", "CDR-H2", "CDR-H3", "CDR-L1", "CDR-L2", "CDR-L3"]
    available_cdr_cols = [c for c in cdr_cols if c in df.columns]
    if available_cdr_cols:
        fig, ax = plt.subplots(figsize=(10, 5))
        x = np.arange(len(available_cdr_cols))
        width = 0.8 / len(methods)
        for i, method in enumerate(methods):
            mdf = df[df["method"] == method]
            means = [mdf[col].dropna().mean() for col in available_cdr_cols]
            stds = [mdf[col].dropna().std() for col in available_cdr_cols]
            offset = x - 0.4 + width * (i + 0.5)
            ax.bar(offset, means, width, yerr=stds, label=method.replace("_", " "),
                   color=colors.get(method, "#555"), capsize=3)
        ax.set_xticks(x)
        ax.set_xticklabels([cdr_labels[cdr_cols.index(c)] for c in available_cdr_cols])
        ax.set_ylabel("Cα RMSD (Å)")
        ax.set_title("CDR Loop RMSD (framework-aligned)")
        ax.legend(fontsize=8)
        plt.tight_layout()
        plt.savefig(figures_dir / "cdr_rmsd_barplot.png", dpi=150)
        plt.close()

    # --- Plot 5: Epitope F1 bar plot ---
    epitope_cols = ["epitope_f1_5A", "epitope_f1_8A"]
    epitope_labels = ["Epitope F1 (5Å)", "Epitope F1 (8Å)"]
    available_ep = [c for c in epitope_cols if c in df.columns]
    if available_ep:
        fig, ax = plt.subplots(figsize=(8, 5))
        x = np.arange(len(available_ep))
        width = 0.8 / len(methods)
        for i, method in enumerate(methods):
            mdf = df[df["method"] == method]
            means = [mdf[col].dropna().mean() for col in available_ep]
            stds = [mdf[col].dropna().std() for col in available_ep]
            offset = x - 0.4 + width * (i + 0.5)
            ax.bar(offset, means, width, yerr=stds, label=method.replace("_", " "),
                   color=colors.get(method, "#555"), capsize=3)
        ax.set_xticks(x)
        ax.set_xticklabels([epitope_labels[epitope_cols.index(c)] for c in available_ep])
        ax.set_ylabel("F1 Score")
        ax.set_title("Epitope Prediction Quality")
        ax.legend(fontsize=8)
        plt.tight_layout()
        plt.savefig(figures_dir / "epitope_f1_barplot.png", dpi=150)
        plt.close()

    # --- Plot 6: Confidence vs DockQ scatter ---
    if "conf_iptm" in df.columns and "top1_ab_ag_dockq" in df.columns:
        fig, ax = plt.subplots(figsize=(7, 5))
        for method in methods:
            mdf = df[df["method"] == method].dropna(subset=["conf_iptm", "top1_ab_ag_dockq"])
            if len(mdf) == 0:
                continue
            ax.scatter(mdf["conf_iptm"], mdf["top1_ab_ag_dockq"],
                       c=colors.get(method, "#555"), label=method.replace("_", " "),
                       alpha=0.6, edgecolors="black", linewidths=0.3)
        ax.set_xlabel("iptm (predicted)")
        ax.set_ylabel("Ab-Ag DockQ (actual)")
        ax.set_title("Confidence Calibration")
        ax.legend(fontsize=8)
        plt.tight_layout()
        plt.savefig(figures_dir / "confidence_vs_dockq.png", dpi=150)
        plt.close()

    print(f"  Plots saved to {figures_dir}")


# ============================================================================
# MAIN
# ============================================================================

def main():
    parser = argparse.ArgumentParser(description="Evaluate Hierarchical Steering")
    parser.add_argument("--ground_truth_dir", type=str, required=True,
                        help="Path to pdb_minimized/ folder with ground truth PDBs")
    parser.add_argument("--predictions_dir", type=str, required=True,
                        help="Path to predictions_examples/ folder")
    parser.add_argument("--cdrs_csv", type=str, required=True,
                        help="Path to cdrs.csv with CDR definitions")
    parser.add_argument("--out_dir", type=str, default="evaluation_results",
                        help="Output directory for results")
    parser.add_argument("--methods", type=str, nargs="+",
                        default=list(METHOD_CONFIGS.keys()),
                        help="Methods to evaluate")
    parser.add_argument("--complexes", type=str, nargs="*", default=None,
                        help="Specific complexes to evaluate (default: all)")
    parser.add_argument("--skip_cdr_rmsd", action="store_true",
                        help="Skip CDR RMSD computation")
    parser.add_argument("--skip_epitope", action="store_true",
                        help="Skip epitope evaluation")
    parser.add_argument("--skip_plots", action="store_true",
                        help="Skip plot generation")
    args = parser.parse_args()

    global GROUND_TRUTH_DIR
    GROUND_TRUTH_DIR = Path(args.ground_truth_dir)
    predictions_dir = Path(args.predictions_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Load CDR definitions
    print("Loading CDR definitions...")
    cdrs = load_cdrs(args.cdrs_csv)
    complex_names = sorted(cdrs.keys())
    if args.complexes:
        complex_names = [c for c in args.complexes if c in cdrs]
    print(f"  {len(complex_names)} complexes to evaluate")

    # Verify ground truths exist
    gt_available = []
    for cname in complex_names:
        gt_path = GROUND_TRUTH_DIR / f"{cname}.pdb"
        if gt_path.exists():
            gt_available.append(cname)
    print(f"  {len(gt_available)} complexes with ground truth PDBs")
    complex_names = gt_available

    # Discover predictions
    print("\nDiscovering prediction folders...")
    predictions = discover_predictions(predictions_dir, complex_names)

    # Filter methods to those requested
    active_methods = [m for m in args.methods if m in predictions and predictions[m]]
    print(f"\nActive methods: {active_methods}")

    # Main evaluation loop
    all_results = []
    for cname in complex_names:
        native_path = str(GROUND_TRUTH_DIR / f"{cname}.pdb")
        cdr_indices = cdrs[cname]

        for method in active_methods:
            run_folders = predictions[method].get(cname, [])
            if not run_folders:
                continue

            print(f"\n  Evaluating {cname} / {method} ({len(run_folders)} settings)...")
            settings_results = evaluate_complex_method(
                cname, run_folders, native_path, cdr_indices
            )

            if not settings_results:
                print(f"    No results for {cname}/{method}")
                continue

            if METHOD_CONFIGS[method]["multi_setting"]:
                # For multi-setting: save best-of-settings result
                best = aggregate_multi_setting(settings_results, "best")
                if best:
                    best["method"] = method
                    best["setting_aggregation"] = "best"
                    all_results.append(best)
            else:
                # Single setting
                for res in settings_results:
                    res["method"] = method
                    res["setting_aggregation"] = "single"
                    all_results.append(res)

    if not all_results:
        print("\nNo results collected. Check prediction folders.")
        return

    # Build DataFrame
    df = pd.DataFrame(all_results)

    # Save detailed per-complex results
    df.to_csv(out_dir / "per_complex_all_metrics.csv", index=False)
    print(f"\nDetailed results saved to {out_dir / 'per_complex_all_metrics.csv'}")

    # ---- Summary table ----
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)

    summary_metrics = [
        ("top1_total_dockq", "Top-1 Total DockQ", True),
        ("top1_ab_ag_dockq", "Top-1 Ab-Ag DockQ", True),
        ("oracle_total_dockq", "Oracle Total DockQ", True),
        ("cdr_cdr3_h", "CDR-H3 RMSD (Å)", False),
        ("epitope_f1_5A", "Epitope F1 (5Å)", True),
        ("epitope_f1_8A", "Epitope F1 (8Å)", True),
        ("conf_iptm", "iptm", True),
    ]

    summary_rows = []
    for method in active_methods:
        mdf = df[df["method"] == method]
        row = {"method": method, "n_complexes": len(mdf)}

        for col, label, higher_better in summary_metrics:
            if col not in mdf.columns:
                continue
            vals = mdf[col].dropna().values
            if len(vals) == 0:
                continue
            mean, ci_lo, ci_hi = bootstrap_ci(vals)
            row[f"{label}_mean"] = round(mean, 4)
            row[f"{label}_ci"] = f"[{ci_lo:.4f}, {ci_hi:.4f}]"
            row[f"{label}_median"] = round(float(np.median(vals)), 4)

        # Success rates
        if "top1_total_dockq" in mdf.columns:
            vals = mdf["top1_total_dockq"].dropna().values
            row["dockq>0.23 (%)"] = round(100 * np.mean(vals >= 0.23), 1)
            row["dockq>0.49 (%)"] = round(100 * np.mean(vals >= 0.49), 1)

        summary_rows.append(row)

    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(out_dir / "summary_table.csv", index=False)

    # Print summary
    print(f"\n{'Method':<22} {'N':>3} {'Top1 DockQ':>11} {'Ab-Ag DQ':>10} "
          f"{'CDR-H3':>8} {'Epi F1':>8} {'iptm':>7} {'DQ>0.23':>8} {'DQ>0.49':>8}")
    print("-" * 100)
    for _, row in summary_df.iterrows():
        print(f"{row['method']:<22} {row.get('n_complexes', 0):>3} "
              f"{row.get('Top-1 Total DockQ_mean', '-'):>11} "
              f"{row.get('Top-1 Ab-Ag DockQ_mean', '-'):>10} "
              f"{row.get('CDR-H3 RMSD (Å)_mean', '-'):>8} "
              f"{row.get('Epitope F1 (5Å)_mean', '-'):>8} "
              f"{row.get('iptm_mean', '-'):>7} "
              f"{row.get('dockq>0.23 (%)', '-'):>8} "
              f"{row.get('dockq>0.49 (%)', '-'):>8}")

    # ---- Statistical tests ----
    print("\n" + "=" * 80)
    print("STATISTICAL TESTS (Hierarchical vs each baseline)")
    print("=" * 80)

    test_metrics = ["top1_total_dockq", "top1_ab_ag_dockq", "cdr_cdr3_h",
                    "epitope_f1_5A", "conf_iptm"]
    stat_rows = []
    for method in active_methods:
        if method == "hierarchical":
            continue
        for metric in test_metrics:
            if metric not in df.columns:
                continue
            result = compute_statistics(df, metric, "hierarchical", method)
            result["comparison"] = f"hierarchical vs {method}"
            result["metric"] = metric
            stat_rows.append(result)
            if result["p_value"] is not None:
                sig = "***" if result["p_value"] < 0.001 else "**" if result["p_value"] < 0.01 else "*" if result["p_value"] < 0.05 else "ns"
                print(f"  {metric:<25} vs {method:<22} "
                      f"diff={result['mean_diff']:+.4f}  "
                      f"W/T/L={result['wins']}/{result['ties']}/{result['losses']}  "
                      f"p={result['p_value']:.4f} {sig}")

    if stat_rows:
        stat_df = pd.DataFrame(stat_rows)
        stat_df.to_csv(out_dir / "statistical_tests.csv", index=False)

    # ---- Plots ----
    if not args.skip_plots:
        print("\nGenerating plots...")
        generate_plots(df, out_dir)

    print(f"\nAll results saved to {out_dir}/")


if __name__ == "__main__":
    main()
