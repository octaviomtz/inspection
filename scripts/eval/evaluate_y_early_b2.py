#!/usr/bin/env python
"""
Evaluation script for Y Early-Phase + B2 Contact Restraints (Experiment 1).

Compares the combination feature against three baselines across ~46 antibody-antigen
complexes.

Methods evaluated:
  baseline         — vanilla Boltz-2, no steering
  contact_restraints — B2: oracle contact restraints (multi-config, best-of)
  pocket_guided    — B3: pocket restraints (multi-config, best-of)
  y_early_b2       — NEW: Y early-phase beta-scaling + B2 contact restraints
                      (multi-config, same restraint configs as B2, best-of)

Metrics:
  1. DockQ (total, Ab-Ag, per-interface) via DockQ v2 Python API
  2. CDR loop RMSD (H1-H3, L1-L3) after antibody-framework alignment
  3. Epitope prediction quality (precision, recall, F1, MCC at 5A and 8A)
  4. Boltz confidence metrics (confidence_score, iptm, ptm, iplddt)

Model selection: top-1 by confidence_score (primary); oracle (best by DockQ)
and mean-over-5 also reported.

Statistical tests: paired Wilcoxon signed-rank (y_early_b2 vs each baseline),
bootstrap 95% CIs on means.

Usage (dockq2 conda env required):
  /home/oc/anaconda3/envs/dockq2/bin/python scripts/eval/evaluate_y_early_b2.py \\
      --ground_truth_dir pdb_minimized \\
      --predictions_dir predictions_examples \\
      --cdrs_csv examples/cdrs.csv \\
      --out_dir evaluation_results_y_early_b2

  # Override prediction subfolder names:
      --y_early_b2_subfolder my_y_early_b2_results
      --hierarchical_v2_subfolder my_hierarchical_results  (optional Y+ alone comparison)

Requires: conda activate dockq2
  (DockQ 2.x, biopython, numpy, pandas, scipy, matplotlib)
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
# DockQ import (requires dockq2 conda env)
# ---------------------------------------------------------------------------
try:
    from DockQ.DockQ import load_PDB, run_on_all_native_interfaces
except ImportError:
    print("ERROR: DockQ not found. Run with: /home/oc/anaconda3/envs/dockq2/bin/python")
    sys.exit(1)

# ---------------------------------------------------------------------------
# BioPython imports
# ---------------------------------------------------------------------------
from Bio.PDB import PDBParser, MMCIFParser, Superimposer
from Bio.SeqUtils import seq1
from Bio.Align import PairwiseAligner


# ============================================================================
# CONFIGURATION
# ============================================================================

# Chains: A = antigen, B = heavy, C = light (standard for this dataset)
AG_CHAIN = "A"
HEAVY_CHAIN = "B"
LIGHT_CHAIN = "C"
CHAIN_MAP = {"A": "A", "B": "B", "C": "C"}

# Method configs: label -> discovery parameters
# subfolder: folder inside predictions_dir
# prefix/suffix: naming convention for the run folder
# multi_setting: True if a complex may have multiple run folders (different restraint types)
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
    # Y+ alone (optional comparison — subfolder overridable via CLI)
    # multi_setting=True so substring matching finds e.g. boltz_results_7TRH_HBG_hierarchical
    "hierarchical_v2": {
        "subfolder": "new_feature_hierarchical",
        "prefix": "boltz_results_",
        "suffix": "",
        "multi_setting": True,
    },
    # NEW: Y early-phase + B2 contact restraints (Experiment 1)
    # multi_setting=True: same restraint configs as B2; best-of selection per complex
    "y_early_b2": {
        "subfolder": "y_early_b2_contact_restraints",
        "prefix": "boltz_results_restraint_",
        "suffix": "",
        "multi_setting": True,
    },
}

# Display colours for plots
METHOD_COLORS = {
    "baseline": "#888888",
    "contact_restraints": "#E88B23",
    "pocket_guided": "#931652",
    "hierarchical_v2": "#188F52",
    "y_early_b2": "#1A6B8A",     # teal/blue — distinct from all others
}
METHOD_DISPLAY = {
    "baseline": "Baseline",
    "contact_restraints": "Contact\nRestraints",
    "pocket_guided": "Pocket\nGuided",
    "hierarchical_v2": "Y+\n(hier. v2)",
    "y_early_b2": "Y-early\n+ B2",
}


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
    """Parse PDB or CIF file with BioPython."""
    filepath = str(filepath)
    if filepath.endswith(".cif"):
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)
    return parser.get_structure("s", filepath)


def get_chain_residues(structure, chain_id):
    """Return standard residues (no HETATM/water) for a chain."""
    try:
        chain = structure[0][chain_id]
    except KeyError:
        return []
    return [r for r in chain.get_residues() if r.get_id()[0] == " "]


def get_chain_sequence(residues):
    """One-letter amino acid sequence from residue list."""
    return "".join(seq1(r.get_resname()) for r in residues)


def build_residue_mapping(gt_residues, pred_residues):
    """Pairwise local sequence alignment to map pred indices to GT indices.

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
# PREDICTION DISCOVERY
# ============================================================================

def discover_predictions(predictions_dir, complex_names):
    """Find prediction run folders for each method × complex.

    Returns: dict[method][complex_name] -> list[Path]
    """
    results = {}
    for method, config in METHOD_CONFIGS.items():
        method_dir = Path(predictions_dir) / config["subfolder"]
        if not method_dir.exists():
            print(f"  Warning: folder not found for method '{method}': {method_dir}")
            results[method] = {}
            continue

        method_results = defaultdict(list)
        for folder in sorted(method_dir.iterdir()):
            if not folder.is_dir():
                continue
            name = folder.name
            for cname in complex_names:
                if config["multi_setting"]:
                    if cname in name:
                        method_results[cname].append(folder)
                        break
                else:
                    expected = f"{config['prefix']}{cname}{config['suffix']}"
                    if name == expected:
                        method_results[cname].append(folder)
                        break

        results[method] = dict(method_results)
        n_c = len(method_results)
        n_t = sum(len(v) for v in method_results.values())
        print(f"  {method:<22}: {n_c} complexes, {n_t} run folders")

    return results


def find_prediction_files(run_folder):
    """Given a boltz_results_XXX/ folder, return list of
    (model_idx, structure_path, confidence_path) for all available models.
    """
    pred_dir = run_folder / "predictions"
    if not pred_dir.exists():
        return []

    subfolders = [d for d in pred_dir.iterdir() if d.is_dir()]
    if not subfolders:
        return []

    inner = subfolders[0]  # e.g. predictions/restraint_7TRH_HBG_hbond_19/
    models = []
    for model_idx in range(10):
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

def compute_dockq(model_path, native_path):
    """Compute DockQ using DockQ v2 Python API.

    Returns dict:
      total_dockq  — mean DockQ across all interfaces
      ab_ag_dockq  — mean DockQ for interfaces involving antigen chain A
      per_interface — {interface_label: {DockQ, iRMSD, LRMSD, fnat, F1, ...}}
    """
    try:
        model = load_PDB(str(model_path))
        native = load_PDB(str(native_path))
        result_dict, _ = run_on_all_native_interfaces(
            model, native, chain_map=CHAIN_MAP
        )
    except Exception as e:
        print(f"    DockQ error ({Path(model_path).name}): {e}")
        return None

    if not result_dict:
        return None

    per_interface = {}
    all_dockq = []
    ab_ag_dockq = []

    for iface, data in result_dict.items():
        per_interface[iface] = {
            "DockQ": data["DockQ"],
            "iRMSD": data["iRMSD"],
            "LRMSD": data["LRMSD"],
            "fnat": data["fnat"],
            "F1": data["F1"],
            "fnonnat": data.get("fnonnat", 0.0),
            "clashes": data.get("clashes", 0),
        }
        all_dockq.append(data["DockQ"])
        if AG_CHAIN in iface:
            ab_ag_dockq.append(data["DockQ"])

    return {
        "total_dockq": float(np.mean(all_dockq)) if all_dockq else 0.0,
        "ab_ag_dockq": float(np.mean(ab_ag_dockq)) if ab_ag_dockq else 0.0,
        "per_interface": per_interface,
    }


# ============================================================================
# METRIC 2: CDR RMSD (after antibody-framework alignment)
# ============================================================================

def compute_cdr_rmsd(model_path, native_path, cdr_indices):
    """Compute per-CDR Cα RMSD after superimposing on the antibody framework.

    Strategy: Superimpose on non-CDR Cα atoms of heavy + light chains, then
    measure Cα RMSD for each CDR loop in that aligned frame.

    Returns dict: 'framework_rmsd', 'cdr1_h', 'cdr2_h', 'cdr3_h',
                  'cdr1_l', 'cdr2_l', 'cdr3_l'
    """
    try:
        model_struct = parse_structure(model_path)
        native_struct = parse_structure(native_path)
    except Exception as e:
        print(f"    Parse error: {e}")
        return None

    # Sequence-align each antibody chain and build index mappings
    mappings = {}
    for chain_id in [HEAVY_CHAIN, LIGHT_CHAIN]:
        gt_res = get_chain_residues(native_struct, chain_id)
        pred_res = get_chain_residues(model_struct, chain_id)
        if not gt_res or not pred_res:
            return None
        mapping = build_residue_mapping(gt_res, pred_res)
        mappings[chain_id] = (gt_res, pred_res, mapping)

    # All CDR positions (0-indexed into the pred residue list)
    cdr_chain_map = {
        "cdr1_h": HEAVY_CHAIN, "cdr2_h": HEAVY_CHAIN, "cdr3_h": HEAVY_CHAIN,
        "cdr1_l": LIGHT_CHAIN, "cdr2_l": LIGHT_CHAIN, "cdr3_l": LIGHT_CHAIN,
    }
    all_cdr_positions = {HEAVY_CHAIN: set(), LIGHT_CHAIN: set()}
    for cdr_name, chain_id in cdr_chain_map.items():
        for pos in cdr_indices[cdr_name]:
            all_cdr_positions[chain_id].add(pos - 1)  # 1-indexed → 0-indexed

    # Framework atoms: non-CDR Cα atoms matched across both structures
    fw_pred_atoms, fw_native_atoms = [], []
    for chain_id in [HEAVY_CHAIN, LIGHT_CHAIN]:
        gt_res, pred_res, mapping = mappings[chain_id]
        cdr_pos = all_cdr_positions[chain_id]
        for pred_idx, gt_idx in mapping.items():
            if pred_idx in cdr_pos:
                continue
            pr = pred_res[pred_idx]
            gr = gt_res[gt_idx]
            if "CA" in pr and "CA" in gr:
                fw_pred_atoms.append(pr["CA"])
                fw_native_atoms.append(gr["CA"])

    if len(fw_pred_atoms) < 10:
        print(f"    Too few framework atoms ({len(fw_pred_atoms)}), skipping RMSD")
        return None

    # Superimpose on framework
    sup = Superimposer()
    sup.set_atoms(fw_native_atoms, fw_pred_atoms)
    sup.apply(list(model_struct.get_atoms()))

    results = {"framework_rmsd": float(sup.rms)}

    # Per-CDR RMSD in the now-aligned frame
    for cdr_name, chain_id in cdr_chain_map.items():
        gt_res, pred_res, mapping = mappings[chain_id]
        positions = [p - 1 for p in cdr_indices[cdr_name]]
        sq_diffs = []
        for pred_idx in positions:
            gt_idx = mapping.get(pred_idx)
            if gt_idx is None:
                continue
            pr = pred_res[pred_idx]
            gr = gt_res[gt_idx]
            if "CA" in pr and "CA" in gr:
                diff = pr["CA"].get_vector().get_array() - gr["CA"].get_vector().get_array()
                sq_diffs.append(float(np.dot(diff, diff)))
        if sq_diffs:
            results[cdr_name] = float(np.sqrt(np.mean(sq_diffs)))

    return results


# ============================================================================
# METRIC 3: Epitope prediction
# ============================================================================

def _heavy_atom_contacts(structure, chain1_id, chain2_id, threshold):
    """Residues on chain1 within `threshold` Å (heavy-atom) of chain2."""
    c1 = get_chain_residues(structure, chain1_id)
    c2 = get_chain_residues(structure, chain2_id)
    if not c1 or not c2:
        return set(), len(c1)

    c2_coords = np.array([
        atom.get_vector().get_array()
        for res in c2 for atom in res.get_atoms()
        if atom.element != "H"
    ])
    if len(c2_coords) == 0:
        return set(), len(c1)

    contacts = set()
    for idx, res in enumerate(c1):
        for atom in res.get_atoms():
            if atom.element == "H":
                continue
            coord = atom.get_vector().get_array()
            if np.sqrt(np.sum((c2_coords - coord) ** 2, axis=1)).min() < threshold:
                contacts.add(idx)
                break
    return contacts, len(c1)


def _ca_contacts(structure, chain1_id, chain2_id, threshold):
    """Residues on chain1 with a Cα within `threshold` Å of any Cα on chain2."""
    c1 = get_chain_residues(structure, chain1_id)
    c2 = get_chain_residues(structure, chain2_id)
    if not c1 or not c2:
        return set(), len(c1)

    c2_ca = np.array([
        res["CA"].get_vector().get_array()
        for res in c2 if "CA" in res
    ])
    if len(c2_ca) == 0:
        return set(), len(c1)

    contacts = set()
    for idx, res in enumerate(c1):
        if "CA" not in res:
            continue
        coord = res["CA"].get_vector().get_array()
        if np.sqrt(np.sum((c2_ca - coord) ** 2, axis=1)).min() < threshold:
            contacts.add(idx)
    return contacts, len(c1)


def compute_epitope_metrics(model_path, native_path):
    """Epitope = antigen residues contacting the antibody.

    Evaluates at two thresholds:
      5Å  — heavy-atom distance (more stringent)
      8Å  — Cα-Cα distance (coarser, standard for epitope mapping)

    Returns dict of precision/recall/F1/MCC at each threshold.
    """
    try:
        model_struct = parse_structure(model_path)
        native_struct = parse_structure(native_path)
    except Exception as e:
        print(f"    Parse error: {e}")
        return None

    results = {}
    for label, contact_fn, threshold in [
        ("5A", _heavy_atom_contacts, 5.0),
        ("8A", _ca_contacts, 8.0),
    ]:
        nat_h, n_ag_nat = contact_fn(native_struct, AG_CHAIN, HEAVY_CHAIN, threshold)
        nat_l, _ = contact_fn(native_struct, AG_CHAIN, LIGHT_CHAIN, threshold)
        native_epitope = nat_h | nat_l

        pred_h, n_ag_pred = contact_fn(model_struct, AG_CHAIN, HEAVY_CHAIN, threshold)
        pred_l, _ = contact_fn(model_struct, AG_CHAIN, LIGHT_CHAIN, threshold)
        model_epitope = pred_h | pred_l

        n_ag = max(n_ag_nat, n_ag_pred)
        tp = len(native_epitope & model_epitope)
        fp = len(model_epitope - native_epitope)
        fn = len(native_epitope - model_epitope)
        tn = max(n_ag - tp - fp - fn, 0)

        precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
        recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
        f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0.0
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
# METRIC 4: Confidence
# ============================================================================

def load_confidence(conf_path):
    """Load Boltz confidence metrics from JSON."""
    if conf_path is None or not Path(conf_path).exists():
        return None
    with open(conf_path) as f:
        data = json.load(f)
    return {
        "confidence_score": data.get("confidence_score", 0.0),
        "iptm": data.get("iptm", 0.0),
        "ptm": data.get("ptm", 0.0),
        "complex_plddt": data.get("complex_plddt", 0.0),
        "complex_iplddt": data.get("complex_iplddt", 0.0),
    }


# ============================================================================
# MODEL SELECTION
# ============================================================================

def select_models(models):
    """Select top-1 model by confidence_score.

    Returns (top1_list_idx, list_of_confidence_dicts).
    """
    confidences = [load_confidence(conf_path) for _, _, conf_path in models]
    scores = [c["confidence_score"] if c else -1.0 for c in confidences]
    return int(np.argmax(scores)), confidences


# ============================================================================
# PER-COMPLEX EVALUATION
# ============================================================================

def evaluate_complex_method(
    complex_name, run_folders, native_path, cdr_indices,
    skip_cdr_rmsd=False, skip_epitope=False,
):
    """Evaluate all models for one complex under one method.

    For multi-setting methods (contact_restraints, pocket_guided, y_early_b2),
    each run_folder is one restraint configuration; we later pick the best.

    Returns list of result dicts (one per run_folder).
    """
    settings_results = []

    for run_folder in run_folders:
        models = find_prediction_files(run_folder)
        if not models:
            continue

        top1_idx, confidences = select_models(models)

        # DockQ for all models (needed for oracle / average)
        all_dockq = [compute_dockq(sp, native_path) for _, sp, _ in models]

        top1_model_idx, top1_struct, _ = models[top1_idx]
        top1_dockq = all_dockq[top1_idx]

        # Oracle: model with highest total DockQ
        valid_pairs = [(i, dq) for i, dq in enumerate(all_dockq) if dq is not None]
        if valid_pairs:
            oracle_idx = max(valid_pairs, key=lambda x: x[1]["total_dockq"])[0]
        else:
            oracle_idx = top1_idx
        oracle_dockq = all_dockq[oracle_idx]

        valid_totals = [dq["total_dockq"] for dq in all_dockq if dq is not None]
        valid_ab_ag = [dq["ab_ag_dockq"] for dq in all_dockq if dq is not None]

        cdr_rmsd = None if skip_cdr_rmsd else compute_cdr_rmsd(top1_struct, native_path, cdr_indices)
        epitope = None if skip_epitope else compute_epitope_metrics(top1_struct, native_path)

        result = {
            "complex": complex_name,
            "run_folder": run_folder.name,
            "n_models": len(models),
            "top1_model_idx": top1_model_idx,
            "oracle_model_idx": models[oracle_idx][0],
        }

        # DockQ metrics
        if top1_dockq:
            result["top1_total_dockq"] = top1_dockq["total_dockq"]
            result["top1_ab_ag_dockq"] = top1_dockq["ab_ag_dockq"]
            for iface, idata in top1_dockq["per_interface"].items():
                for k, v in idata.items():
                    result[f"top1_{iface}_{k}"] = v
        if oracle_dockq:
            result["oracle_total_dockq"] = oracle_dockq["total_dockq"]
            result["oracle_ab_ag_dockq"] = oracle_dockq["ab_ag_dockq"]
        if valid_totals:
            result["avg_total_dockq"] = float(np.mean(valid_totals))
            result["avg_ab_ag_dockq"] = float(np.mean(valid_ab_ag))

        # CDR RMSD
        if cdr_rmsd:
            for k, v in cdr_rmsd.items():
                result[f"cdr_{k}"] = v

        # Epitope
        if epitope:
            result.update(epitope)

        # Confidence
        if confidences[top1_idx]:
            for k, v in confidences[top1_idx].items():
                result[f"conf_{k}"] = v

        settings_results.append(result)

    return settings_results


def best_of_settings(settings_results):
    """For multi-setting methods, pick the setting with highest top1_total_dockq."""
    if not settings_results:
        return None
    if len(settings_results) == 1:
        return settings_results[0]
    valid = [r for r in settings_results if "top1_total_dockq" in r]
    return max(valid, key=lambda r: r["top1_total_dockq"]) if valid else settings_results[0]


# ============================================================================
# STATISTICS
# ============================================================================

def wilcoxon_test(df, metric, method1, method2):
    """Paired Wilcoxon signed-rank test of method1 vs method2 on metric."""
    d1 = df[df["method"] == method1][["complex", metric]].dropna()
    d2 = df[df["method"] == method2][["complex", metric]].dropna()
    merged = pd.merge(d1, d2, on="complex", suffixes=("_1", "_2"))
    if len(merged) < 5:
        return {"n": len(merged), "p_value": None, "statistic": None}

    v1 = merged[f"{metric}_1"].values
    v2 = merged[f"{metric}_2"].values
    diff = v1 - v2

    wins = int(np.sum(diff > 1e-4))
    losses = int(np.sum(diff < -1e-4))
    ties = len(diff) - wins - losses

    try:
        stat, p_val = stats.wilcoxon(v1, v2, alternative="two-sided")
    except Exception:
        stat, p_val = None, None

    return {
        "n": len(merged),
        "mean_diff": float(np.mean(diff)),
        "median_diff": float(np.median(diff)),
        "wins": wins, "losses": losses, "ties": ties,
        "statistic": stat,
        "p_value": p_val,
    }


def bootstrap_ci(values, n_boot=1000):
    """Mean and 95% bootstrap CI."""
    arr = np.asarray(values, dtype=float)
    arr = arr[~np.isnan(arr)]
    if len(arr) == 0:
        return float("nan"), float("nan"), float("nan")
    boot = [np.random.choice(arr, len(arr), replace=True).mean() for _ in range(n_boot)]
    return float(arr.mean()), float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))


def capri_rates(values):
    """Fraction of complexes in each CAPRI quality category."""
    arr = np.asarray(values, dtype=float)
    arr = arr[~np.isnan(arr)]
    if len(arr) == 0:
        return {}
    return {
        "incorrect_pct": 100 * float(np.mean(arr < 0.23)),
        "acceptable_pct": 100 * float(np.mean((arr >= 0.23) & (arr < 0.49))),
        "medium_pct": 100 * float(np.mean((arr >= 0.49) & (arr < 0.80))),
        "high_pct": 100 * float(np.mean(arr >= 0.80)),
        "medium_or_high_pct": 100 * float(np.mean(arr >= 0.49)),
    }


# ============================================================================
# PLOTTING
# ============================================================================

def generate_plots(df, out_dir, active_methods):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("  matplotlib not available, skipping plots")
        return

    fig_dir = Path(out_dir) / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    colors = [METHOD_COLORS.get(m, "#555555") for m in active_methods]
    labels = [METHOD_DISPLAY.get(m, m) for m in active_methods]

    # ---- 1. DockQ bar chart (3 panels) ----
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    for ax, col, title in zip(
        axes,
        ["top1_total_dockq", "top1_ab_ag_dockq", "oracle_total_dockq"],
        ["Top-1 Total DockQ", "Top-1 Ab-Ag DockQ", "Oracle Total DockQ"],
    ):
        means, lowers, uppers = [], [], []
        for m in active_methods:
            vals = df[df["method"] == m][col].dropna().values
            mn, lo, hi = bootstrap_ci(vals)
            means.append(mn); lowers.append(mn - lo); uppers.append(hi - mn)
        x = np.arange(len(active_methods))
        ax.bar(x, means, yerr=[lowers, uppers], capsize=4, color=colors)
        ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=8)
        ax.set_ylim(0, min(1, max(means) * 1.4 + 0.1))
        ax.set_ylabel("DockQ"); ax.set_title(title)
        for xi, mn in zip(x, means):
            ax.text(xi, mn + 0.01, f"{mn:.3f}", ha="center", va="bottom", fontsize=7)
    fig.suptitle("DockQ Metrics — Y Early-Phase + B2 vs Baselines", fontweight="bold")
    plt.tight_layout()
    plt.savefig(fig_dir / "dockq_barplot.png", dpi=150)
    plt.close()

    # ---- 2. DockQ boxplot ----
    fig, ax = plt.subplots(figsize=(8, 5))
    data_bp = [df[df["method"] == m]["top1_ab_ag_dockq"].dropna().values for m in active_methods]
    bp = ax.boxplot(data_bp, labels=labels, patch_artist=True)
    for patch, c in zip(bp["boxes"], colors):
        patch.set_facecolor(c); patch.set_alpha(0.7)
    ax.set_ylabel("Ab-Ag DockQ (top-1 by confidence)")
    ax.set_title("DockQ Distribution Across Complexes")
    plt.tight_layout()
    plt.savefig(fig_dir / "dockq_boxplot.png", dpi=150)
    plt.close()

    # ---- 3a. Scatter: Y_early_B2 vs B2 (primary comparison) ----
    if "y_early_b2" in active_methods and "contact_restraints" in active_methods:
        d_new = df[df["method"] == "y_early_b2"][["complex", "top1_ab_ag_dockq"]].dropna()
        d_b2 = df[df["method"] == "contact_restraints"][["complex", "top1_ab_ag_dockq"]].dropna()
        merged = pd.merge(d_new, d_b2, on="complex", suffixes=("_new", "_b2"))
        if len(merged) > 0:
            fig, ax = plt.subplots(figsize=(6, 6))
            ax.scatter(
                merged["top1_ab_ag_dockq_b2"], merged["top1_ab_ag_dockq_new"],
                c=METHOD_COLORS["y_early_b2"], alpha=0.7, edgecolors="black", linewidths=0.5,
            )
            for _, row in merged.iterrows():
                diff = row["top1_ab_ag_dockq_new"] - row["top1_ab_ag_dockq_b2"]
                if abs(diff) > 0.10:
                    ax.annotate(row["complex"], (row["top1_ab_ag_dockq_b2"], row["top1_ab_ag_dockq_new"]),
                                fontsize=5, alpha=0.7)
            ax.plot([0, 1], [0, 1], "k--", alpha=0.4, label="y=x (no change)")
            n_above = int((merged["top1_ab_ag_dockq_new"] > merged["top1_ab_ag_dockq_b2"]).sum())
            n_below = int((merged["top1_ab_ag_dockq_new"] < merged["top1_ab_ag_dockq_b2"]).sum())
            ax.set_xlabel("B2 Contact Restraints Ab-Ag DockQ")
            ax.set_ylabel("Y-early + B2 Ab-Ag DockQ")
            ax.set_title(f"Y-early+B2 vs B2 (per complex)\n"
                         f"Above diagonal = improvement  ({n_above}↑ / {n_below}↓)")
            ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_aspect("equal")
            ax.legend(fontsize=8)
            plt.tight_layout()
            plt.savefig(fig_dir / "scatter_y_early_b2_vs_b2.png", dpi=150)
            plt.close()

    # ---- 3b. Scatter: Y_early_B2 vs baseline ----
    if "y_early_b2" in active_methods and "baseline" in active_methods:
        d_new = df[df["method"] == "y_early_b2"][["complex", "top1_ab_ag_dockq"]].dropna()
        d_base = df[df["method"] == "baseline"][["complex", "top1_ab_ag_dockq"]].dropna()
        merged = pd.merge(d_new, d_base, on="complex", suffixes=("_new", "_base"))
        if len(merged) > 0:
            fig, ax = plt.subplots(figsize=(6, 6))
            ax.scatter(
                merged["top1_ab_ag_dockq_base"], merged["top1_ab_ag_dockq_new"],
                c=METHOD_COLORS["y_early_b2"], alpha=0.7, edgecolors="black", linewidths=0.5,
            )
            for _, row in merged.iterrows():
                diff = row["top1_ab_ag_dockq_new"] - row["top1_ab_ag_dockq_base"]
                if abs(diff) > 0.15:
                    ax.annotate(row["complex"], (row["top1_ab_ag_dockq_base"], row["top1_ab_ag_dockq_new"]),
                                fontsize=5, alpha=0.7)
            ax.plot([0, 1], [0, 1], "k--", alpha=0.4, label="y=x")
            ax.set_xlabel("Baseline Ab-Ag DockQ")
            ax.set_ylabel("Y-early + B2 Ab-Ag DockQ")
            ax.set_title("Y-early+B2 vs Baseline (per complex)\nAbove diagonal = improvement")
            ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_aspect("equal")
            ax.legend(fontsize=8)
            plt.tight_layout()
            plt.savefig(fig_dir / "scatter_y_early_b2_vs_baseline.png", dpi=150)
            plt.close()

    # ---- 4. CDR RMSD bar chart ----
    cdr_cols = ["cdr_cdr1_h", "cdr_cdr2_h", "cdr_cdr3_h", "cdr_cdr1_l", "cdr_cdr2_l", "cdr_cdr3_l"]
    cdr_labels = ["CDR-H1", "CDR-H2", "CDR-H3", "CDR-L1", "CDR-L2", "CDR-L3"]
    avail_cdr = [c for c in cdr_cols if c in df.columns]
    if avail_cdr:
        fig, ax = plt.subplots(figsize=(11, 5))
        x = np.arange(len(avail_cdr))
        w = 0.8 / len(active_methods)
        for i, (m, c) in enumerate(zip(active_methods, colors)):
            mdf = df[df["method"] == m]
            means = [mdf[col].dropna().mean() for col in avail_cdr]
            stds = [mdf[col].dropna().std() for col in avail_cdr]
            offset = x - 0.4 + w * (i + 0.5)
            ax.bar(offset, means, w, yerr=stds, label=METHOD_DISPLAY.get(m, m),
                   color=c, capsize=3)
        ax.set_xticks(x)
        ax.set_xticklabels([cdr_labels[cdr_cols.index(c)] for c in avail_cdr])
        ax.set_ylabel("Cα RMSD (Å)")
        ax.set_title("CDR Loop RMSD after Framework Alignment\n(lower is better)")
        ax.legend(fontsize=8)
        plt.tight_layout()
        plt.savefig(fig_dir / "cdr_rmsd_barplot.png", dpi=150)
        plt.close()

    # ---- 4b. CDR-H3 RMSD scatter: Y_early_B2 vs B2 (primary hypothesis) ----
    h3_col = "cdr_cdr3_h"
    if (h3_col in df.columns
            and "y_early_b2" in active_methods
            and "contact_restraints" in active_methods):
        d_new = df[df["method"] == "y_early_b2"][["complex", h3_col]].dropna()
        d_b2 = df[df["method"] == "contact_restraints"][["complex", h3_col]].dropna()
        merged = pd.merge(d_new, d_b2, on="complex", suffixes=("_new", "_b2"))
        if len(merged) > 0:
            fig, ax = plt.subplots(figsize=(6, 6))
            ax.scatter(
                merged[f"{h3_col}_b2"], merged[f"{h3_col}_new"],
                c=METHOD_COLORS["y_early_b2"], alpha=0.7, edgecolors="black", linewidths=0.5,
            )
            for _, row in merged.iterrows():
                diff = row[f"{h3_col}_new"] - row[f"{h3_col}_b2"]
                if abs(diff) > 0.5:
                    ax.annotate(row["complex"], (row[f"{h3_col}_b2"], row[f"{h3_col}_new"]),
                                fontsize=5, alpha=0.7)
            max_val = max(merged[[f"{h3_col}_new", f"{h3_col}_b2"]].max()) * 1.05
            ax.plot([0, max_val], [0, max_val], "k--", alpha=0.4, label="y=x (no change)")
            n_better = int((merged[f"{h3_col}_new"] < merged[f"{h3_col}_b2"]).sum())
            n_worse = int((merged[f"{h3_col}_new"] > merged[f"{h3_col}_b2"]).sum())
            ax.set_xlabel("B2 CDR-H3 RMSD (Å)")
            ax.set_ylabel("Y-early + B2 CDR-H3 RMSD (Å)")
            ax.set_title(f"CDR-H3 RMSD: Y-early+B2 vs B2\n"
                         f"Below diagonal = improvement  ({n_better}↓ / {n_worse}↑)")
            ax.set_xlim(0, max_val); ax.set_ylim(0, max_val); ax.set_aspect("equal")
            ax.legend(fontsize=8)
            plt.tight_layout()
            plt.savefig(fig_dir / "scatter_cdrh3_rmsd_y_early_b2_vs_b2.png", dpi=150)
            plt.close()

    # ---- 5. Epitope F1 ----
    ep_cols = ["epitope_f1_5A", "epitope_f1_8A"]
    ep_labels = ["Epitope F1 (5Å)", "Epitope F1 (8Å)"]
    avail_ep = [c for c in ep_cols if c in df.columns]
    if avail_ep:
        fig, ax = plt.subplots(figsize=(8, 5))
        x = np.arange(len(avail_ep))
        w = 0.8 / len(active_methods)
        for i, (m, c) in enumerate(zip(active_methods, colors)):
            mdf = df[df["method"] == m]
            means = [mdf[col].dropna().mean() for col in avail_ep]
            stds = [mdf[col].dropna().std() for col in avail_ep]
            offset = x - 0.4 + w * (i + 0.5)
            ax.bar(offset, means, w, yerr=stds, label=METHOD_DISPLAY.get(m, m),
                   color=c, capsize=3)
        ax.set_xticks(x); ax.set_xticklabels([ep_labels[ep_cols.index(c)] for c in avail_ep])
        ax.set_ylabel("F1 Score"); ax.set_title("Epitope Prediction Quality (higher is better)")
        ax.legend(fontsize=8)
        plt.tight_layout()
        plt.savefig(fig_dir / "epitope_f1_barplot.png", dpi=150)
        plt.close()

    # ---- 6. Confidence calibration (iptm vs DockQ) ----
    if "conf_iptm" in df.columns and "top1_ab_ag_dockq" in df.columns:
        fig, ax = plt.subplots(figsize=(7, 5))
        for m, c in zip(active_methods, colors):
            sub = df[df["method"] == m].dropna(subset=["conf_iptm", "top1_ab_ag_dockq"])
            if len(sub) == 0:
                continue
            corr = sub["conf_iptm"].corr(sub["top1_ab_ag_dockq"])
            ax.scatter(sub["conf_iptm"], sub["top1_ab_ag_dockq"],
                       c=c, label=f"{METHOD_DISPLAY.get(m, m)} (r={corr:.2f})",
                       alpha=0.6, edgecolors="black", linewidths=0.3)
        ax.set_xlabel("iptm (predicted confidence)")
        ax.set_ylabel("Ab-Ag DockQ (actual)")
        ax.set_title("Confidence Calibration")
        ax.legend(fontsize=7)
        plt.tight_layout()
        plt.savefig(fig_dir / "confidence_calibration.png", dpi=150)
        plt.close()

    # ---- 7. CAPRI stacked bar ----
    capri_data = {}
    for m in active_methods:
        vals = df[df["method"] == m]["top1_ab_ag_dockq"].dropna().values
        capri_data[m] = capri_rates(vals)

    if capri_data:
        cats = ["incorrect_pct", "acceptable_pct", "medium_pct", "high_pct"]
        cat_labels = ["Incorrect\n(<0.23)", "Acceptable\n(0.23–0.49)",
                      "Medium\n(0.49–0.80)", "High\n(>0.80)"]
        cat_colors = ["#d73027", "#fc8d59", "#fee090", "#4dac26"]
        x = np.arange(len(active_methods))
        bottoms = np.zeros(len(active_methods))
        fig, ax = plt.subplots(figsize=(9, 5))
        for cat, cat_col, cat_lbl in zip(cats, cat_colors, cat_labels):
            vals = [capri_data[m].get(cat, 0) for m in active_methods]
            ax.bar(x, vals, bottom=bottoms, label=cat_lbl, color=cat_col)
            bottoms += np.array(vals)
        ax.set_xticks(x); ax.set_xticklabels(labels)
        ax.set_ylabel("% of complexes"); ax.set_ylim(0, 105)
        ax.set_title("CAPRI Quality Classification (Ab-Ag DockQ)")
        ax.legend(loc="upper right", fontsize=8)
        plt.tight_layout()
        plt.savefig(fig_dir / "capri_stacked_bar.png", dpi=150)
        plt.close()

    print(f"  Plots saved to {fig_dir}/")


# ============================================================================
# MAIN
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Evaluate Y Early-Phase + B2 Contact Restraints vs baselines",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--ground_truth_dir", required=True,
                        help="Path to pdb_minimized/ with ground truth PDBs")
    parser.add_argument("--predictions_dir", required=True,
                        help="Path to predictions_examples/ parent folder")
    parser.add_argument("--cdrs_csv", required=True,
                        help="Path to cdrs.csv with CDR definitions")
    parser.add_argument("--out_dir", default="evaluation_results_y_early_b2",
                        help="Output directory (default: evaluation_results_y_early_b2)")
    parser.add_argument("--y_early_b2_subfolder", default="y_early_b2_contact_restraints",
                        help="Subfolder name inside predictions_dir for Y-early+B2 results "
                             "(default: y_early_b2_contact_restraints)")
    parser.add_argument("--hierarchical_v2_subfolder", default="new_feature_hierarchical",
                        help="Subfolder name for Y+ alone results, if present "
                             "(default: new_feature_hierarchical)")
    parser.add_argument("--methods", nargs="+",
                        default=["baseline", "contact_restraints", "pocket_guided", "y_early_b2"],
                        help="Methods to include. Default excludes hierarchical_v2 (Y+ alone). "
                             "Add 'hierarchical_v2' to include it as a comparison point.")
    parser.add_argument("--complexes", nargs="*", default=None,
                        help="Restrict to specific complex names (default: all in cdrs.csv)")
    parser.add_argument("--skip_cdr_rmsd", action="store_true",
                        help="Skip CDR RMSD computation (faster)")
    parser.add_argument("--skip_epitope", action="store_true",
                        help="Skip epitope metric computation (faster)")
    parser.add_argument("--skip_plots", action="store_true",
                        help="Skip plot generation")
    args = parser.parse_args()

    # Apply subfolder overrides
    METHOD_CONFIGS["y_early_b2"]["subfolder"] = args.y_early_b2_subfolder
    METHOD_CONFIGS["hierarchical_v2"]["subfolder"] = args.hierarchical_v2_subfolder

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Load CDR definitions
    print("Loading CDR definitions...")
    cdrs = load_cdrs(args.cdrs_csv)
    complex_names = sorted(cdrs.keys())
    if args.complexes:
        complex_names = [c for c in args.complexes if c in cdrs]
    print(f"  {len(complex_names)} complexes loaded")

    # Verify ground truths
    gt_dir = Path(args.ground_truth_dir)
    complex_names = [c for c in complex_names if (gt_dir / f"{c}.pdb").exists()]
    print(f"  {len(complex_names)} complexes with ground truth PDB")

    # Discover predictions
    print("\nDiscovering prediction folders...")
    predictions = discover_predictions(Path(args.predictions_dir), complex_names)

    active_methods = [m for m in args.methods if m in predictions and predictions[m]]
    print(f"\nActive methods: {active_methods}")

    # Evaluation loop
    all_results = []
    n_total = len(complex_names) * len(active_methods)
    done = 0
    for cname in complex_names:
        native_path = str(gt_dir / f"{cname}.pdb")
        cdr_indices = cdrs[cname]

        for method in active_methods:
            run_folders = predictions[method].get(cname, [])
            if not run_folders:
                done += 1
                continue

            print(f"  [{done+1}/{n_total}] {cname} / {method} "
                  f"({len(run_folders)} setting{'s' if len(run_folders)>1 else ''})...")

            settings = evaluate_complex_method(
                cname, run_folders, native_path, cdr_indices,
                skip_cdr_rmsd=args.skip_cdr_rmsd,
                skip_epitope=args.skip_epitope,
            )
            done += 1

            if not settings:
                print(f"    No results for {cname}/{method}")
                continue

            if METHOD_CONFIGS[method]["multi_setting"]:
                best = best_of_settings(settings)
                if best:
                    best["method"] = method
                    all_results.append(best)
            else:
                for r in settings:
                    r["method"] = method
                    all_results.append(r)

    if not all_results:
        print("\nNo results collected. Check folder names and paths.")
        return

    df = pd.DataFrame(all_results)
    df.to_csv(out_dir / "per_complex_all_metrics.csv", index=False)
    print(f"\nDetailed results → {out_dir / 'per_complex_all_metrics.csv'}")

    # ---- Summary table ----
    print("\n" + "=" * 90)
    print("SUMMARY TABLE")
    print("=" * 90)

    summary_metrics = [
        ("top1_total_dockq",      "Top-1 Total DockQ",   True),
        ("top1_ab_ag_dockq",      "Top-1 Ab-Ag DockQ",   True),
        ("oracle_total_dockq",    "Oracle Total DockQ",   True),
        ("avg_total_dockq",       "Mean Total DockQ",     True),
        ("cdr_cdr3_h",            "CDR-H3 RMSD (Å)",     False),
        ("cdr_cdr3_l",            "CDR-L3 RMSD (Å)",     False),
        ("epitope_f1_5A",         "Epitope F1 (5Å)",      True),
        ("epitope_f1_8A",         "Epitope F1 (8Å)",      True),
        ("conf_iptm",             "iptm",                 True),
        ("conf_confidence_score", "confidence_score",     True),
    ]

    summary_rows = []
    for method in active_methods:
        mdf = df[df["method"] == method]
        row = {"method": method, "n_complexes": len(mdf)}
        for col, label, _ in summary_metrics:
            if col not in mdf.columns:
                continue
            vals = mdf[col].dropna().values
            if len(vals) == 0:
                continue
            mean, ci_lo, ci_hi = bootstrap_ci(vals)
            row[f"{col}_mean"] = round(mean, 4)
            row[f"{col}_ci95"] = f"[{ci_lo:.4f}, {ci_hi:.4f}]"
            row[f"{col}_median"] = round(float(np.median(vals)), 4)
        # CAPRI rates
        if "top1_ab_ag_dockq" in mdf.columns:
            for k, v in capri_rates(mdf["top1_ab_ag_dockq"].dropna().values).items():
                row[k] = round(v, 1)
        summary_rows.append(row)

    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(out_dir / "summary_table.csv", index=False)

    # Console print
    hdr = (f"{'Method':<24} {'N':>3}  {'Top1 DQ':>8}  {'Ab-Ag DQ':>9}  "
           f"{'H3 RMSD':>8}  {'EpiF1@8':>8}  {'iptm':>6}  "
           f"{'DQ>=0.23':>8}  {'Med+Hi%':>8}")
    print(hdr)
    print("-" * len(hdr))
    for _, row in summary_df.iterrows():
        print(
            f"{row['method']:<24} {row.get('n_complexes', 0):>3}  "
            f"{row.get('top1_total_dockq_mean', float('nan')):>8.4f}  "
            f"{row.get('top1_ab_ag_dockq_mean', float('nan')):>9.4f}  "
            f"{row.get('cdr_cdr3_h_mean', float('nan')):>8.4f}  "
            f"{row.get('epitope_f1_8A_mean', float('nan')):>8.4f}  "
            f"{row.get('conf_iptm_mean', float('nan')):>6.4f}  "
            f"{row.get('acceptable_pct', 0) + row.get('medium_pct', 0) + row.get('high_pct', 0):>8.1f}  "
            f"{row.get('medium_or_high_pct', float('nan')):>8.1f}"
        )

    # ---- Statistical tests ----
    # Primary: y_early_b2 vs each other method
    # Also include pairwise y_early_b2 vs contact_restraints (the key comparison)
    print("\n" + "=" * 90)
    print("STATISTICAL TESTS — y_early_b2 vs each other method (paired Wilcoxon)")
    print("=" * 90)
    print(f"  {'Metric':<28}  {'vs':<24}  {'mean diff':>10}  "
          f"{'W/T/L':>10}  {'p-value':>9}  {'sig':>4}")
    print("  " + "-" * 90)

    test_metrics = [
        "top1_total_dockq", "top1_ab_ag_dockq",
        "cdr_cdr3_h", "cdr_cdr3_l",
        "epitope_f1_5A", "epitope_f1_8A", "conf_iptm",
    ]
    stat_rows = []

    # y_early_b2 vs every other active method
    if "y_early_b2" in active_methods:
        for other_method in [m for m in active_methods if m != "y_early_b2"]:
            for metric in test_metrics:
                if metric not in df.columns:
                    continue
                res = wilcoxon_test(df, metric, "y_early_b2", other_method)
                res["comparison"] = f"y_early_b2 vs {other_method}"
                res["metric"] = metric
                stat_rows.append(res)
                if res["p_value"] is not None:
                    p = res["p_value"]
                    sig = "***" if p < 0.001 else ("**" if p < 0.01 else ("*" if p < 0.05 else "ns"))
                    w, t, l = res["wins"], res["ties"], res["losses"]
                    print(f"  {metric:<28}  vs {other_method:<24}  "
                          f"{res['mean_diff']:>+10.4f}  "
                          f"{w:>3}/{t:>3}/{l:>3}  "
                          f"{p:>9.4f}  {sig:>4}")

    if stat_rows:
        pd.DataFrame(stat_rows).to_csv(out_dir / "statistical_tests.csv", index=False)
        print(f"\nStatistical tests → {out_dir / 'statistical_tests.csv'}")

    # ---- Plots ----
    if not args.skip_plots:
        print("\nGenerating plots...")
        generate_plots(df, out_dir, active_methods)

    print(f"\nAll outputs in: {out_dir}/")


if __name__ == "__main__":
    main()
