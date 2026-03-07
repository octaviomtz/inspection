#!/usr/bin/env python
"""
Evaluation script for Strategy W (Embedding-Based Interface Steering).

Compares predictions from Strategy W against three baselines:
  B1: antigen_cut (no steering)
  B2: antigen_cut_contact_restraints (coordinate-space restraints)
  B3: antigen_cut_vhvl_msa_pocket (pocket constraints)
  W:  new_feature_embed_interface (embedding-based interface steering)

Runs in conda environment: dockq2
  conda activate dockq2
  python evaluate_embed_interface.py

Evaluation tiers:
  Tier 1: DockQ (CAPRI-standard interface docking quality)
  Tier 2: Epitope recovery (residue-level contact prediction)
  Tier 3: Buried surface area (SASA-based interface properties)
  Tier 4: Confidence correlation analysis

Usage:
  python evaluate_embed_interface.py                    # Run all steps
  python evaluate_embed_interface.py --step index       # Step 1: build index
  python evaluate_embed_interface.py --step dockq       # Step 2: run DockQ
  python evaluate_embed_interface.py --step epitope     # Step 3: epitope metrics
  python evaluate_embed_interface.py --step sasa        # Step 4: SASA/BSA metrics
  python evaluate_embed_interface.py --step summary     # Step 5: aggregate + plots
"""

import argparse
import copy
import json
import math
import os
import subprocess
import sys
import tempfile
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from Bio.PDB import MMCIFParser, NeighborSearch, PDBParser
from Bio.PDB.SASA import ShrakeRupley

warnings.filterwarnings("ignore")

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

BASE_DIR = Path(__file__).resolve().parent
GROUND_TRUTH_DIR = BASE_DIR / "pdb_minimized"
PREDICTIONS_DIR = BASE_DIR / "predictions_examples"
OUTPUT_DIR = BASE_DIR / "evaluation_results"

# Condition folder configuration
CONDITIONS = {
    "B1_no_steering": {
        "dir": PREDICTIONS_DIR / "antigen_cut",
        "prefix": "boltz_results_",
        "suffix": "",
        "multi_variant": False,
    },
    "B2_contact_restraints": {
        "dir": PREDICTIONS_DIR / "antigen_cut_contact_restraints",
        "prefix": "boltz_results_restraint_",
        "suffix": "",
        "multi_variant": True,
    },
    "B3_pocket_msa": {
        "dir": PREDICTIONS_DIR / "antigen_cut_vhvl_msa_pocket_ab_boltz_post_2023_omm",
        "prefix": "boltz_results_restraint_to_A_",
        "suffix": "",
        "multi_variant": True,
    },
    "W_embed_interface": {
        "dir": PREDICTIONS_DIR / "new_feature_embed_interface",
        "prefix": "boltz_results_",
        "suffix": "_embed_interface",
        "multi_variant": False,
    },
}

# Contact distance threshold for epitope definition (Angstroms)
EPITOPE_CONTACT_THRESHOLD = 4.5

# Chain mapping: A=antigen, B=heavy, C=light
ANTIGEN_CHAIN = "A"
HEAVY_CHAIN = "B"
LIGHT_CHAIN = "C"
ALL_CHAINS = [ANTIGEN_CHAIN, HEAVY_CHAIN, LIGHT_CHAIN]

# All 48 complexes
COMPLEXES = sorted([p.stem for p in GROUND_TRUTH_DIR.glob("*.pdb")])


# ---------------------------------------------------------------------------
# Step 1: Discover and index all predictions
# ---------------------------------------------------------------------------

def extract_complex_name_from_folder(folder_name, condition_key):
    """Extract the complex name (e.g. 7TRH_HBG) from a result folder name.

    Returns (complex_name, variant_label) where variant_label is None for
    single-variant conditions, or a string like 'hbond_23' for multi-variant.
    """
    cfg = CONDITIONS[condition_key]

    name = folder_name
    if cfg["prefix"] and name.startswith(cfg["prefix"]):
        name = name[len(cfg["prefix"]):]
    if cfg["suffix"] and name.endswith(cfg["suffix"]):
        name = name[: -len(cfg["suffix"])]

    if not cfg["multi_variant"]:
        return name, None

    # For multi-variant conditions, the complex name is embedded in the folder
    # name. Try to match against known complexes.
    for cplx in COMPLEXES:
        if name.startswith(cplx):
            variant = name[len(cplx):].lstrip("_")
            return cplx, variant if variant else "default"
        # B3 pattern: 7TRH_HBG_B_W_109 -> complex=7TRH_HBG, variant=B_W_109
        if cplx in name:
            idx = name.index(cplx)
            variant = name[idx + len(cplx):].lstrip("_")
            return cplx, variant if variant else "default"

    return name, None


def find_model_files(predictions_subdir):
    """Find all structure model files (.pdb or .cif) in a predictions subdirectory.

    Returns list of (model_idx, model_path, confidence_path).
    """
    results = []
    if not predictions_subdir.exists():
        return results

    for f in sorted(predictions_subdir.iterdir()):
        if f.suffix not in (".pdb", ".cif"):
            continue
        if f.name.startswith("confidence_") or f.name.startswith("pae_") or \
           f.name.startswith("plddt_") or f.name.startswith("pde_"):
            continue

        # Extract model index from name like 7TRH_HBG_model_0.pdb
        stem = f.stem
        model_idx = None
        if "_model_" in stem:
            try:
                model_idx = int(stem.split("_model_")[-1])
            except ValueError:
                model_idx = 0
        else:
            model_idx = 0

        # Find corresponding confidence file
        conf_name = f"confidence_{stem}.json"
        conf_path = f.parent / conf_name
        if not conf_path.exists():
            conf_path = None

        results.append((model_idx, f, conf_path))

    return results


def load_confidence(conf_path):
    """Load confidence metrics from a JSON file."""
    if conf_path is None or not conf_path.exists():
        return {}
    with open(conf_path) as f:
        return json.load(f)


def build_index():
    """Step 1: Discover and index all predictions across all conditions."""
    print("=" * 70)
    print("STEP 1: Building prediction index")
    print("=" * 70)

    rows = []
    for condition_key, cfg in CONDITIONS.items():
        cond_dir = cfg["dir"]
        if not cond_dir.exists():
            print(f"  WARNING: {condition_key} directory not found: {cond_dir}")
            continue

        for result_folder in sorted(cond_dir.iterdir()):
            if not result_folder.is_dir():
                continue

            complex_name, variant = extract_complex_name_from_folder(
                result_folder.name, condition_key
            )

            if complex_name not in COMPLEXES:
                print(f"  WARNING: complex '{complex_name}' from folder "
                      f"'{result_folder.name}' not in ground truth set, skipping")
                continue

            # Find the inner predictions directory
            pred_dir = result_folder / "predictions"
            if not pred_dir.exists():
                continue

            # The predictions folder contains a subdirectory named after the YAML
            inner_dirs = [d for d in pred_dir.iterdir() if d.is_dir()]
            if not inner_dirs:
                continue

            for inner_dir in inner_dirs:
                models = find_model_files(inner_dir)
                for model_idx, model_path, conf_path in models:
                    conf = load_confidence(conf_path)
                    rows.append({
                        "complex": complex_name,
                        "condition": condition_key,
                        "variant": variant or "default",
                        "model_idx": model_idx,
                        "model_path": str(model_path),
                        "confidence_path": str(conf_path) if conf_path else "",
                        "confidence_score": conf.get("confidence_score", float("nan")),
                        "complex_iplddt": conf.get("complex_iplddt", float("nan")),
                        "complex_plddt": conf.get("complex_plddt", float("nan")),
                        "iptm": conf.get("iptm", float("nan")),
                        "ptm": conf.get("ptm", float("nan")),
                    })

    df = pd.DataFrame(rows)
    if df.empty:
        print("  No predictions found!")
        return df

    # Summary
    for condition in df["condition"].unique():
        sub = df[df["condition"] == condition]
        n_complexes = sub["complex"].nunique()
        n_models = len(sub)
        n_variants = sub["variant"].nunique()
        print(f"  {condition}: {n_complexes} complexes, {n_models} models, "
              f"{n_variants} variants")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    out_path = OUTPUT_DIR / "evaluation_index.csv"
    df.to_csv(out_path, index=False)
    print(f"\n  Saved index to {out_path}")
    return df


# ---------------------------------------------------------------------------
# Step 2: Run DockQ on all predictions
# ---------------------------------------------------------------------------

def run_dockq_single(model_path, native_path, json_out_path):
    """Run DockQ on a single prediction vs native, return parsed results."""
    cmd = [
        "DockQ", str(model_path), str(native_path),
        "--mapping", "ABC:ABC",
        "--json", str(json_out_path),
    ]
    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True, timeout=120
        )
    except subprocess.TimeoutExpired:
        print(f"    DockQ timeout for {model_path}")
        return None
    except FileNotFoundError:
        print("    ERROR: DockQ command not found. Activate dockq2 environment.")
        return None

    if result.returncode != 0:
        print(f"    DockQ failed for {model_path}: {result.stderr[:200]}")
        return None

    if not Path(json_out_path).exists():
        return None

    with open(json_out_path) as f:
        dockq_data = json.load(f)

    return dockq_data


def parse_dockq_result(dockq_data, complex_name, condition, variant, model_idx):
    """Parse DockQ JSON output into a flat row dict."""
    if dockq_data is None:
        return None

    best = dockq_data.get("best_result", {})

    row = {
        "complex": complex_name,
        "condition": condition,
        "variant": variant,
        "model_idx": model_idx,
        "total_dockq": dockq_data.get("best_dockq", float("nan")),
    }

    # Extract per-interface metrics
    for interface_key, iface_label in [("AB", "AB"), ("AC", "AC"), ("BC", "BC")]:
        idata = best.get(interface_key, {})
        row[f"DockQ_{iface_label}"] = idata.get("DockQ", float("nan"))
        row[f"iRMSD_{iface_label}"] = idata.get("iRMSD", float("nan"))
        row[f"LRMSD_{iface_label}"] = idata.get("LRMSD", float("nan"))
        row[f"fnat_{iface_label}"] = idata.get("fnat", float("nan"))
        row[f"fnonnat_{iface_label}"] = idata.get("fnonnat", float("nan"))
        row[f"F1_{iface_label}"] = idata.get("F1", float("nan"))

    # Compute antigen-antibody average DockQ (AB + AC only, excluding BC)
    dq_ab = row.get("DockQ_AB", float("nan"))
    dq_ac = row.get("DockQ_AC", float("nan"))
    if not math.isnan(dq_ab) and not math.isnan(dq_ac):
        row["DockQ_antigen_avg"] = (dq_ab + dq_ac) / 2
    else:
        row["DockQ_antigen_avg"] = float("nan")

    return row


def run_dockq_all(index_df):
    """Step 2: Run DockQ for all indexed predictions."""
    print("\n" + "=" * 70)
    print("STEP 2: Running DockQ")
    print("=" * 70)

    dockq_dir = OUTPUT_DIR / "dockq_json"
    dockq_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    total = len(index_df)
    for i, (_, entry) in enumerate(index_df.iterrows()):
        model_path = Path(entry["model_path"])
        native_path = GROUND_TRUTH_DIR / f"{entry['complex']}.pdb"

        if not model_path.exists():
            print(f"  [{i+1}/{total}] SKIP (model missing): {model_path}")
            continue
        if not native_path.exists():
            print(f"  [{i+1}/{total}] SKIP (native missing): {native_path}")
            continue

        json_name = (
            f"{entry['complex']}_{entry['condition']}_{entry['variant']}"
            f"_model{entry['model_idx']}.json"
        )
        json_out = dockq_dir / json_name

        print(f"  [{i+1}/{total}] {entry['condition']} / {entry['complex']} "
              f"/ model_{entry['model_idx']}",
              end=" ... ", flush=True)

        dockq_data = run_dockq_single(model_path, native_path, json_out)
        row = parse_dockq_result(
            dockq_data,
            entry["complex"], entry["condition"],
            entry["variant"], entry["model_idx"],
        )
        if row is not None:
            rows.append(row)
            print(f"DockQ_AB={row['DockQ_AB']:.3f}, "
                  f"DockQ_AC={row['DockQ_AC']:.3f}")
        else:
            print("FAILED")

    df = pd.DataFrame(rows)
    if not df.empty:
        out_path = OUTPUT_DIR / "dockq_results.csv"
        df.to_csv(out_path, index=False)
        print(f"\n  Saved DockQ results to {out_path}")
    return df


# ---------------------------------------------------------------------------
# Step 3: Epitope recovery metrics
# ---------------------------------------------------------------------------

def load_structure(path):
    """Load a structure from PDB or CIF file."""
    path = Path(path)
    if path.suffix == ".cif":
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)
    return parser.get_structure(path.stem, str(path))


def get_heavy_atoms(chain):
    """Get heavy (non-hydrogen) atoms from a chain."""
    return [a for a in chain.get_atoms() if a.element != "H"]


def compute_epitope_residues(structure, antigen_chain_id, ab_chain_ids,
                             threshold=EPITOPE_CONTACT_THRESHOLD):
    """Find antigen residues in contact with antibody chains.

    Returns set of (chain_id, resseq) tuples for contacted antigen residues.
    """
    model = structure[0]

    try:
        ag_chain = model[antigen_chain_id]
    except KeyError:
        return set()

    ab_atoms = []
    for chain_id in ab_chain_ids:
        try:
            ab_atoms.extend(get_heavy_atoms(model[chain_id]))
        except KeyError:
            pass

    if not ab_atoms:
        return set()

    ns = NeighborSearch(ab_atoms)
    epitope = set()
    for atom in get_heavy_atoms(ag_chain):
        neighbors = ns.search(atom.get_vector().get_array(), threshold, "A")
        if neighbors:
            res = atom.get_parent()
            epitope.add(res.get_id()[1])

    return epitope


def compute_epitope_metrics(true_epitope, pred_epitope, all_antigen_residues):
    """Compute epitope classification metrics."""
    true_set = set(true_epitope)
    pred_set = set(pred_epitope)
    all_set = set(all_antigen_residues)

    tp = len(true_set & pred_set)
    fp = len(pred_set - true_set)
    fn = len(true_set - pred_set)
    tn = len(all_set - true_set - pred_set)

    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0.0
    jaccard = tp / (tp + fp + fn) if (tp + fp + fn) > 0 else 0.0

    # MCC
    denom = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    mcc = (tp * tn - fp * fn) / denom if denom > 0 else 0.0

    return {
        "epitope_precision": precision,
        "epitope_recall": recall,
        "epitope_f1": f1,
        "epitope_mcc": mcc,
        "epitope_jaccard": jaccard,
        "n_true_epitope": len(true_set),
        "n_pred_epitope": len(pred_set),
        "n_tp": tp,
        "n_fp": fp,
        "n_fn": fn,
    }


def run_epitope_all(index_df):
    """Step 3: Compute epitope recovery metrics for best model per condition/complex."""
    print("\n" + "=" * 70)
    print("STEP 3: Computing epitope recovery metrics")
    print("=" * 70)

    # Select best model per (complex, condition, variant) by complex_iplddt
    best_df = select_best_models(index_df)

    # Cache ground truth epitopes
    gt_cache = {}

    rows = []
    total = len(best_df)
    for i, (_, entry) in enumerate(best_df.iterrows()):
        complex_name = entry["complex"]
        model_path = Path(entry["model_path"])
        native_path = GROUND_TRUTH_DIR / f"{complex_name}.pdb"

        if not model_path.exists() or not native_path.exists():
            continue

        print(f"  [{i+1}/{total}] {entry['condition']} / {complex_name}", end=" ... ", flush=True)

        # Ground truth epitope (cached)
        if complex_name not in gt_cache:
            native_struct = load_structure(native_path)
            gt_epitope = compute_epitope_residues(
                native_struct, ANTIGEN_CHAIN, [HEAVY_CHAIN, LIGHT_CHAIN]
            )
            # Also get all antigen residue numbers
            try:
                all_ag_res = [r.get_id()[1] for r in native_struct[0][ANTIGEN_CHAIN].get_residues()]
            except KeyError:
                all_ag_res = []
            gt_cache[complex_name] = (gt_epitope, all_ag_res)

        gt_epitope, all_ag_res = gt_cache[complex_name]

        # Predicted epitope
        pred_struct = load_structure(model_path)
        pred_epitope = compute_epitope_residues(
            pred_struct, ANTIGEN_CHAIN, [HEAVY_CHAIN, LIGHT_CHAIN]
        )

        # Metrics
        metrics = compute_epitope_metrics(gt_epitope, pred_epitope, all_ag_res)
        metrics["complex"] = complex_name
        metrics["condition"] = entry["condition"]
        metrics["variant"] = entry["variant"]
        metrics["model_idx"] = entry["model_idx"]
        rows.append(metrics)

        print(f"F1={metrics['epitope_f1']:.3f}, "
              f"Prec={metrics['epitope_precision']:.3f}, "
              f"Rec={metrics['epitope_recall']:.3f}")

    df = pd.DataFrame(rows)
    if not df.empty:
        out_path = OUTPUT_DIR / "epitope_results.csv"
        df.to_csv(out_path, index=False)
        print(f"\n  Saved epitope results to {out_path}")
    return df


# ---------------------------------------------------------------------------
# Step 4: SASA / BSA metrics
# ---------------------------------------------------------------------------

def compute_chain_sasa(structure, chain_ids):
    """Compute total SASA for specified chains in context of the full complex."""
    model = structure[0]
    sr = ShrakeRupley()
    sr.compute(model, level="R")

    total_sasa = 0.0
    for chain_id in chain_ids:
        try:
            for res in model[chain_id].get_residues():
                total_sasa += res.sasa
        except KeyError:
            pass
    return total_sasa


def compute_isolated_chain_sasa(structure, chain_id):
    """Compute SASA of a single chain in isolation (no other chains present).

    Creates a temporary structure with only the specified chain.
    """
    model = structure[0]
    try:
        chain = model[chain_id]
    except KeyError:
        return 0.0

    # Build a new structure with only this chain
    from Bio.PDB import StructureBuilder
    sb = StructureBuilder.StructureBuilder()
    sb.init_structure("isolated")
    sb.init_model(0)
    sb.init_seg(" ")

    isolated_struct = sb.get_structure()
    isolated_model = isolated_struct[0]

    chain_copy = chain.copy()
    isolated_model.add(chain_copy)

    sr = ShrakeRupley()
    sr.compute(isolated_model, level="R")

    total_sasa = 0.0
    for res in isolated_model[chain_id].get_residues():
        total_sasa += res.sasa
    return total_sasa


def compute_bsa_metrics(structure):
    """Compute buried surface area for antigen-antibody interface.

    BSA = SASA_unbound(antigen) + SASA_unbound(antibody) - SASA_bound(complex)
    """
    model = structure[0]

    # Compute SASA for full complex
    sr = ShrakeRupley()
    sr.compute(model, level="R")

    # Collect per-residue SASA in complex
    ag_sasa_bound = {}
    try:
        for res in model[ANTIGEN_CHAIN].get_residues():
            ag_sasa_bound[res.get_id()[1]] = res.sasa
    except KeyError:
        return None

    # Compute antigen SASA in isolation
    ag_sasa_unbound_total = compute_isolated_chain_sasa(structure, ANTIGEN_CHAIN)

    # Compute antibody SASA in isolation (B + C together)
    from Bio.PDB import StructureBuilder
    sb = StructureBuilder.StructureBuilder()
    sb.init_structure("ab_only")
    sb.init_model(0)
    sb.init_seg(" ")
    ab_struct = sb.get_structure()
    ab_model = ab_struct[0]
    for cid in [HEAVY_CHAIN, LIGHT_CHAIN]:
        try:
            ab_model.add(model[cid].copy())
        except KeyError:
            pass

    sr2 = ShrakeRupley()
    sr2.compute(ab_model, level="R")
    ab_sasa_unbound_total = 0.0
    for cid in [HEAVY_CHAIN, LIGHT_CHAIN]:
        try:
            for res in ab_model[cid].get_residues():
                ab_sasa_unbound_total += res.sasa
        except KeyError:
            pass

    # SASA in complex
    ag_sasa_bound_total = sum(ag_sasa_bound.values())
    ab_sasa_bound_total = 0.0
    for cid in [HEAVY_CHAIN, LIGHT_CHAIN]:
        try:
            for res in model[cid].get_residues():
                ab_sasa_bound_total += res.sasa
        except KeyError:
            pass

    # BSA
    bsa = (ag_sasa_unbound_total + ab_sasa_unbound_total) - (
        ag_sasa_bound_total + ab_sasa_bound_total
    )

    # Compute per-residue unbound SASA for antigen to count buried residues
    ag_iso_struct = load_structure.__wrapped__(structure) if hasattr(load_structure, '__wrapped__') else None
    # Simpler: recompute isolated antigen per-residue
    sb3 = StructureBuilder.StructureBuilder()
    sb3.init_structure("ag_only")
    sb3.init_model(0)
    sb3.init_seg(" ")
    ag_only = sb3.get_structure()
    ag_only_model = ag_only[0]
    try:
        ag_only_model.add(model[ANTIGEN_CHAIN].copy())
    except KeyError:
        return None

    sr3 = ShrakeRupley()
    sr3.compute(ag_only_model, level="R")
    ag_sasa_unbound = {}
    for res in ag_only_model[ANTIGEN_CHAIN].get_residues():
        ag_sasa_unbound[res.get_id()[1]] = res.sasa

    # Count antigen residues that become buried (delta SASA >= 1 A^2)
    buried_count = 0
    for resid in ag_sasa_unbound:
        if resid in ag_sasa_bound:
            delta = ag_sasa_unbound[resid] - ag_sasa_bound[resid]
            if delta >= 1.0:
                buried_count += 1

    return {
        "bsa": bsa,
        "ag_sasa_unbound": ag_sasa_unbound_total,
        "ag_sasa_bound": ag_sasa_bound_total,
        "ab_sasa_unbound": ab_sasa_unbound_total,
        "ab_sasa_bound": ab_sasa_bound_total,
        "buried_residue_count": buried_count,
    }


def run_sasa_all(index_df):
    """Step 4: Compute SASA/BSA metrics for best model per condition/complex."""
    print("\n" + "=" * 70)
    print("STEP 4: Computing SASA/BSA metrics")
    print("=" * 70)

    best_df = select_best_models(index_df)

    # Also compute for ground truth (cached)
    gt_bsa_cache = {}

    rows = []
    total = len(best_df)
    for i, (_, entry) in enumerate(best_df.iterrows()):
        complex_name = entry["complex"]
        model_path = Path(entry["model_path"])
        native_path = GROUND_TRUTH_DIR / f"{complex_name}.pdb"

        if not model_path.exists() or not native_path.exists():
            continue

        print(f"  [{i+1}/{total}] {entry['condition']} / {complex_name}", end=" ... ", flush=True)

        # Ground truth BSA (cached)
        if complex_name not in gt_bsa_cache:
            native_struct = load_structure(native_path)
            gt_bsa = compute_bsa_metrics(native_struct)
            gt_bsa_cache[complex_name] = gt_bsa

        gt_bsa = gt_bsa_cache[complex_name]

        # Prediction BSA
        pred_struct = load_structure(model_path)
        pred_bsa = compute_bsa_metrics(pred_struct)

        if pred_bsa is None:
            print("FAILED")
            continue

        row = {
            "complex": complex_name,
            "condition": entry["condition"],
            "variant": entry["variant"],
            "model_idx": entry["model_idx"],
            "pred_bsa": pred_bsa["bsa"],
            "pred_buried_count": pred_bsa["buried_residue_count"],
        }

        if gt_bsa is not None:
            row["native_bsa"] = gt_bsa["bsa"]
            row["native_buried_count"] = gt_bsa["buried_residue_count"]
            row["bsa_ratio"] = pred_bsa["bsa"] / gt_bsa["bsa"] if gt_bsa["bsa"] > 0 else float("nan")
            row["buried_count_ratio"] = (
                pred_bsa["buried_residue_count"] / gt_bsa["buried_residue_count"]
                if gt_bsa["buried_residue_count"] > 0 else float("nan")
            )

        rows.append(row)
        print(f"BSA={pred_bsa['bsa']:.0f} Å², buried={pred_bsa['buried_residue_count']}")

    df = pd.DataFrame(rows)
    if not df.empty:
        out_path = OUTPUT_DIR / "sasa_results.csv"
        df.to_csv(out_path, index=False)
        print(f"\n  Saved SASA results to {out_path}")
    return df


# ---------------------------------------------------------------------------
# Helper: model selection
# ---------------------------------------------------------------------------

def select_best_models(index_df):
    """Select the best model per (complex, condition, variant) by complex_iplddt.

    Falls back to first model if no confidence data is available.
    """
    if index_df.empty:
        return index_df

    best_indices = []
    for _, group in index_df.groupby(["complex", "condition", "variant"]):
        valid = group.dropna(subset=["complex_iplddt"])
        if not valid.empty:
            best_indices.append(valid["complex_iplddt"].idxmax())
        else:
            best_indices.append(group.index[0])

    return index_df.loc[best_indices].reset_index(drop=True)


def aggregate_multi_variant(df, metric_cols, mode="best"):
    """For multi-variant conditions (B2, B3), aggregate to one row per complex.

    mode='best': keep variant with highest DockQ_antigen_avg
    mode='mean': average across variants
    """
    result_rows = []
    for condition in df["condition"].unique():
        cond_df = df[df["condition"] == condition]
        cfg = CONDITIONS.get(condition, {})

        if cfg.get("multi_variant", False) and mode == "best":
            for complex_name in cond_df["complex"].unique():
                cplx_df = cond_df[cond_df["complex"] == complex_name]
                # Pick variant with best primary metric
                primary = "DockQ_antigen_avg" if "DockQ_antigen_avg" in cplx_df.columns else metric_cols[0]
                valid = cplx_df.dropna(subset=[primary])
                if not valid.empty:
                    result_rows.append(valid.loc[valid[primary].idxmax()].to_dict())
                else:
                    result_rows.append(cplx_df.iloc[0].to_dict())
        elif cfg.get("multi_variant", False) and mode == "mean":
            for complex_name in cond_df["complex"].unique():
                cplx_df = cond_df[cond_df["complex"] == complex_name]
                row = {"complex": complex_name, "condition": condition, "variant": "mean"}
                for col in metric_cols:
                    if col in cplx_df.columns:
                        row[col] = cplx_df[col].mean()
                result_rows.append(row)
        else:
            for _, row in cond_df.iterrows():
                result_rows.append(row.to_dict())

    return pd.DataFrame(result_rows)


# ---------------------------------------------------------------------------
# Step 5: Aggregate and compare
# ---------------------------------------------------------------------------

def capri_category(dockq):
    """Classify DockQ score into CAPRI categories."""
    if math.isnan(dockq):
        return "N/A"
    if dockq >= 0.80:
        return "High"
    if dockq >= 0.49:
        return "Medium"
    if dockq >= 0.23:
        return "Acceptable"
    return "Incorrect"


def run_summary(dockq_df, epitope_df, sasa_df, index_df):
    """Step 5: Aggregate results and generate comparison plots."""
    print("\n" + "=" * 70)
    print("STEP 5: Aggregating results and generating plots")
    print("=" * 70)

    plots_dir = OUTPUT_DIR / "evaluation_plots"
    plots_dir.mkdir(parents=True, exist_ok=True)

    # --- Merge DockQ with confidence data ---
    if not dockq_df.empty:
        # Select best model per (complex, condition, variant) from DockQ results
        dockq_best = select_best_dockq(dockq_df)

        # For multi-variant conditions, pick best variant per complex
        dockq_agg_best = aggregate_multi_variant(
            dockq_best,
            ["DockQ_AB", "DockQ_AC", "DockQ_antigen_avg", "total_dockq"],
            mode="best",
        )
        dockq_agg_mean = aggregate_multi_variant(
            dockq_best,
            ["DockQ_AB", "DockQ_AC", "DockQ_antigen_avg", "total_dockq"],
            mode="mean",
        )

        print("\n--- DockQ Summary (best model, best variant) ---")
        print_dockq_summary(dockq_agg_best)

        print("\n--- DockQ Summary (best model, mean variant for B2/B3) ---")
        print_dockq_summary(dockq_agg_mean)

        # Pairwise comparisons: W vs each baseline
        w_label = "W_embed_interface"
        for baseline in ["B1_no_steering", "B2_contact_restraints", "B3_pocket_msa"]:
            print_pairwise_comparison(dockq_agg_best, w_label, baseline, "DockQ_antigen_avg")

        # Generate plots
        plot_dockq_comparison(dockq_agg_best, plots_dir)
        plot_per_complex_delta(dockq_agg_best, w_label, plots_dir)
        plot_capri_categories(dockq_agg_best, plots_dir)

    # --- Epitope summary ---
    if epitope_df is not None and not epitope_df.empty:
        epi_agg = aggregate_multi_variant(
            epitope_df,
            ["epitope_f1", "epitope_precision", "epitope_recall", "epitope_mcc", "epitope_jaccard"],
            mode="best",
        )
        print("\n--- Epitope Recovery Summary ---")
        print_epitope_summary(epi_agg)
        plot_epitope_comparison(epi_agg, plots_dir)

    # --- SASA summary ---
    if sasa_df is not None and not sasa_df.empty:
        sasa_agg = aggregate_multi_variant(
            sasa_df,
            ["pred_bsa", "bsa_ratio", "pred_buried_count", "buried_count_ratio"],
            mode="best",
        )
        print("\n--- SASA / BSA Summary ---")
        print_sasa_summary(sasa_agg)
        plot_sasa_comparison(sasa_agg, plots_dir)

    # --- Confidence correlation ---
    if not dockq_df.empty and not index_df.empty:
        print("\n--- Confidence Correlation ---")
        compute_confidence_correlation(dockq_df, index_df)

    # Save combined summary
    save_combined_summary(dockq_df, epitope_df, sasa_df)

    print(f"\n  Plots saved to {plots_dir}/")


def select_best_dockq(dockq_df):
    """Select best model per (complex, condition, variant) by DockQ_antigen_avg."""
    if dockq_df.empty:
        return dockq_df

    best_indices = []
    for _, group in dockq_df.groupby(["complex", "condition", "variant"]):
        valid = group.dropna(subset=["DockQ_antigen_avg"])
        if not valid.empty:
            best_indices.append(valid["DockQ_antigen_avg"].idxmax())
        else:
            best_indices.append(group.index[0])

    return dockq_df.loc[best_indices].reset_index(drop=True)


def print_dockq_summary(df):
    """Print DockQ summary table per condition."""
    metrics = ["DockQ_AB", "DockQ_AC", "DockQ_antigen_avg", "DockQ_BC",
               "total_dockq", "fnat_AB", "fnat_AC", "iRMSD_AB", "iRMSD_AC"]
    for condition in sorted(df["condition"].unique()):
        sub = df[df["condition"] == condition]
        print(f"\n  {condition} ({len(sub)} complexes):")
        for m in metrics:
            if m in sub.columns:
                vals = sub[m].dropna()
                if len(vals) > 0:
                    print(f"    {m:25s}: mean={vals.mean():.3f}  "
                          f"median={vals.median():.3f}  std={vals.std():.3f}")

        # CAPRI categories for antigen interfaces
        for iface in ["DockQ_AB", "DockQ_AC"]:
            if iface in sub.columns:
                cats = sub[iface].dropna().apply(capri_category).value_counts()
                cats_str = ", ".join(f"{k}:{v}" for k, v in sorted(cats.items()))
                print(f"    CAPRI {iface}: {cats_str}")


def print_pairwise_comparison(df, condition_a, condition_b, metric):
    """Print pairwise comparison between two conditions."""
    a = df[df["condition"] == condition_a][["complex", metric]].rename(
        columns={metric: "a"}
    )
    b = df[df["condition"] == condition_b][["complex", metric]].rename(
        columns={metric: "b"}
    )
    merged = a.merge(b, on="complex", how="inner")

    if merged.empty:
        print(f"\n  {condition_a} vs {condition_b}: no overlapping complexes")
        return

    merged["delta"] = merged["a"] - merged["b"]
    n = len(merged)
    wins = (merged["delta"] > 0.02).sum()
    losses = (merged["delta"] < -0.02).sum()
    ties = n - wins - losses

    print(f"\n  {condition_a} vs {condition_b} ({metric}, n={n}):")
    print(f"    Mean Δ: {merged['delta'].mean():+.4f}")
    print(f"    Median Δ: {merged['delta'].median():+.4f}")
    print(f"    Wins/Ties/Losses: {wins}/{ties}/{losses}")

    if n >= 5:
        try:
            stat_result = stats.wilcoxon(
                merged["a"].values, merged["b"].values, alternative="two-sided"
            )
            print(f"    Wilcoxon p-value: {stat_result.pvalue:.4f}")
        except Exception:
            print("    Wilcoxon test: could not compute (insufficient data)")


def print_epitope_summary(df):
    """Print epitope recovery summary."""
    metrics = ["epitope_f1", "epitope_precision", "epitope_recall",
               "epitope_mcc", "epitope_jaccard"]
    for condition in sorted(df["condition"].unique()):
        sub = df[df["condition"] == condition]
        print(f"\n  {condition} ({len(sub)} complexes):")
        for m in metrics:
            if m in sub.columns:
                vals = sub[m].dropna()
                if len(vals) > 0:
                    print(f"    {m:25s}: mean={vals.mean():.3f}  "
                          f"median={vals.median():.3f}")


def print_sasa_summary(df):
    """Print SASA/BSA summary."""
    metrics = ["pred_bsa", "bsa_ratio", "pred_buried_count", "buried_count_ratio"]
    for condition in sorted(df["condition"].unique()):
        sub = df[df["condition"] == condition]
        print(f"\n  {condition} ({len(sub)} complexes):")
        for m in metrics:
            if m in sub.columns:
                vals = sub[m].dropna()
                if len(vals) > 0:
                    print(f"    {m:25s}: mean={vals.mean():.3f}  "
                          f"median={vals.median():.3f}")


def compute_confidence_correlation(dockq_df, index_df):
    """Compute correlation between model confidence and DockQ quality."""
    # Merge confidence from index with DockQ scores
    merged = dockq_df.merge(
        index_df[["complex", "condition", "variant", "model_idx",
                   "confidence_score", "complex_iplddt"]],
        on=["complex", "condition", "variant", "model_idx"],
        how="left",
    )

    for condition in sorted(merged["condition"].unique()):
        sub = merged[merged["condition"] == condition]

        for conf_metric, dockq_metric in [
            ("complex_iplddt", "DockQ_antigen_avg"),
            ("confidence_score", "total_dockq"),
        ]:
            valid = sub.dropna(subset=[conf_metric, dockq_metric])
            if len(valid) < 5:
                continue
            try:
                r_pearson, p_pearson = stats.pearsonr(valid[conf_metric], valid[dockq_metric])
                r_spearman, p_spearman = stats.spearmanr(valid[conf_metric], valid[dockq_metric])
                print(f"  {condition}: {conf_metric} vs {dockq_metric}")
                print(f"    Pearson r={r_pearson:.3f} (p={p_pearson:.4f}), "
                      f"Spearman ρ={r_spearman:.3f} (p={p_spearman:.4f})")
            except Exception:
                pass


def save_combined_summary(dockq_df, epitope_df, sasa_df):
    """Save a combined summary CSV with one row per (complex, condition)."""
    dfs = []
    if dockq_df is not None and not dockq_df.empty:
        dfs.append(dockq_df)
    if epitope_df is not None and not epitope_df.empty:
        join_cols = ["complex", "condition", "variant", "model_idx"]
        if dfs:
            dfs[0] = dfs[0].merge(epitope_df, on=join_cols, how="outer",
                                   suffixes=("", "_epi"))
        else:
            dfs.append(epitope_df)
    if sasa_df is not None and not sasa_df.empty:
        join_cols = ["complex", "condition", "variant", "model_idx"]
        if dfs:
            dfs[0] = dfs[0].merge(sasa_df, on=join_cols, how="outer",
                                   suffixes=("", "_sasa"))
        else:
            dfs.append(sasa_df)

    if dfs:
        combined = dfs[0]
        out_path = OUTPUT_DIR / "evaluation_summary.csv"
        combined.to_csv(out_path, index=False)
        print(f"\n  Saved combined summary to {out_path}")


# ---------------------------------------------------------------------------
# Plotting functions
# ---------------------------------------------------------------------------

CONDITION_COLORS = {
    "B1_no_steering": "#4C72B0",
    "B2_contact_restraints": "#55A868",
    "B3_pocket_msa": "#C44E52",
    "W_embed_interface": "#8172B2",
}

CONDITION_SHORT = {
    "B1_no_steering": "B1: No steering",
    "B2_contact_restraints": "B2: Contact restr.",
    "B3_pocket_msa": "B3: Pocket+MSA",
    "W_embed_interface": "W: Embed. interface",
}


def plot_dockq_comparison(df, plots_dir):
    """Box plot comparing DockQ across conditions."""
    conditions = [c for c in CONDITION_SHORT if c in df["condition"].unique()]
    if len(conditions) < 2:
        return

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    for ax, metric, title in zip(
        axes,
        ["DockQ_AB", "DockQ_AC", "DockQ_antigen_avg"],
        ["DockQ: Heavy–Antigen (AB)", "DockQ: Light–Antigen (AC)",
         "DockQ: Antigen avg (AB+AC)/2"],
    ):
        data = []
        labels = []
        colors = []
        for cond in conditions:
            vals = df[df["condition"] == cond][metric].dropna().values
            if len(vals) > 0:
                data.append(vals)
                labels.append(CONDITION_SHORT[cond])
                colors.append(CONDITION_COLORS[cond])

        if not data:
            continue

        bp = ax.boxplot(data, labels=labels, patch_artist=True, widths=0.6)
        for patch, color in zip(bp["boxes"], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)
        ax.set_title(title, fontsize=11)
        ax.set_ylabel("DockQ")
        ax.tick_params(axis="x", rotation=30)
        ax.set_ylim(-0.05, 1.05)
        ax.axhline(y=0.23, color="gray", linestyle="--", alpha=0.5, linewidth=0.8)
        ax.axhline(y=0.49, color="gray", linestyle="--", alpha=0.5, linewidth=0.8)
        ax.axhline(y=0.80, color="gray", linestyle="--", alpha=0.5, linewidth=0.8)

    plt.tight_layout()
    plt.savefig(plots_dir / "dockq_boxplot_comparison.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("  Saved dockq_boxplot_comparison.png")


def plot_per_complex_delta(df, w_condition, plots_dir):
    """Bar plot showing per-complex ΔDockQ (W vs B1)."""
    b1_condition = "B1_no_steering"
    metric = "DockQ_antigen_avg"

    w_df = df[df["condition"] == w_condition][["complex", metric]].rename(
        columns={metric: "w"}
    )
    b1_df = df[df["condition"] == b1_condition][["complex", metric]].rename(
        columns={metric: "b1"}
    )
    merged = w_df.merge(b1_df, on="complex", how="inner")
    if merged.empty:
        return

    merged["delta"] = merged["w"] - merged["b1"]
    merged = merged.sort_values("delta", ascending=True)

    fig, ax = plt.subplots(figsize=(14, max(6, len(merged) * 0.3)))
    colors = ["#55A868" if d > 0 else "#C44E52" for d in merged["delta"]]
    ax.barh(range(len(merged)), merged["delta"].values, color=colors, alpha=0.8)
    ax.set_yticks(range(len(merged)))
    ax.set_yticklabels(merged["complex"].values, fontsize=7)
    ax.set_xlabel("ΔDockQ (W − B1)")
    ax.set_title("Per-complex DockQ improvement: W vs B1 (no steering)")
    ax.axvline(x=0, color="black", linewidth=0.8)

    plt.tight_layout()
    plt.savefig(plots_dir / "delta_dockq_W_vs_B1.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("  Saved delta_dockq_W_vs_B1.png")


def plot_capri_categories(df, plots_dir):
    """Stacked bar chart of CAPRI categories per condition."""
    conditions = [c for c in CONDITION_SHORT if c in df["condition"].unique()]
    if len(conditions) < 2:
        return

    metric = "DockQ_antigen_avg"
    categories = ["Incorrect", "Acceptable", "Medium", "High"]
    cat_colors = {"Incorrect": "#d9534f", "Acceptable": "#f0ad4e",
                  "Medium": "#5bc0de", "High": "#5cb85c"}

    fig, ax = plt.subplots(figsize=(10, 5))
    x = range(len(conditions))
    width = 0.6
    bottom = np.zeros(len(conditions))

    for cat in categories:
        counts = []
        for cond in conditions:
            vals = df[df["condition"] == cond][metric].dropna()
            n_total = len(vals)
            n_cat = (vals.apply(capri_category) == cat).sum()
            counts.append(n_cat / n_total * 100 if n_total > 0 else 0)
        ax.bar(x, counts, width, bottom=bottom, label=cat, color=cat_colors[cat], alpha=0.85)
        bottom += counts

    ax.set_xticks(x)
    ax.set_xticklabels([CONDITION_SHORT[c] for c in conditions], rotation=20)
    ax.set_ylabel("% of complexes")
    ax.set_title("CAPRI Quality Categories (antigen interface)")
    ax.legend(loc="upper right")
    ax.set_ylim(0, 105)

    plt.tight_layout()
    plt.savefig(plots_dir / "capri_categories.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("  Saved capri_categories.png")


def plot_epitope_comparison(df, plots_dir):
    """Box plot comparing epitope metrics across conditions."""
    conditions = [c for c in CONDITION_SHORT if c in df["condition"].unique()]
    if len(conditions) < 2:
        return

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    for ax, metric, title in zip(
        axes,
        ["epitope_f1", "epitope_precision", "epitope_recall"],
        ["Epitope F1", "Epitope Precision", "Epitope Recall"],
    ):
        data = []
        labels = []
        colors = []
        for cond in conditions:
            vals = df[df["condition"] == cond][metric].dropna().values
            if len(vals) > 0:
                data.append(vals)
                labels.append(CONDITION_SHORT[cond])
                colors.append(CONDITION_COLORS[cond])

        if not data:
            continue

        bp = ax.boxplot(data, labels=labels, patch_artist=True, widths=0.6)
        for patch, color in zip(bp["boxes"], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)
        ax.set_title(title, fontsize=11)
        ax.set_ylabel(title)
        ax.tick_params(axis="x", rotation=30)
        ax.set_ylim(-0.05, 1.05)

    plt.tight_layout()
    plt.savefig(plots_dir / "epitope_boxplot_comparison.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("  Saved epitope_boxplot_comparison.png")


def plot_sasa_comparison(df, plots_dir):
    """Box plot comparing BSA and buried residue count across conditions."""
    conditions = [c for c in CONDITION_SHORT if c in df["condition"].unique()]
    if len(conditions) < 2:
        return

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for ax, metric, title in zip(
        axes,
        ["bsa_ratio", "buried_count_ratio"],
        ["BSA ratio (pred/native)", "Buried residue ratio (pred/native)"],
    ):
        data = []
        labels = []
        colors = []
        for cond in conditions:
            vals = df[df["condition"] == cond][metric].dropna().values
            if len(vals) > 0:
                data.append(vals)
                labels.append(CONDITION_SHORT[cond])
                colors.append(CONDITION_COLORS[cond])

        if not data:
            continue

        bp = ax.boxplot(data, labels=labels, patch_artist=True, widths=0.6)
        for patch, color in zip(bp["boxes"], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)
        ax.set_title(title, fontsize=11)
        ax.set_ylabel("Ratio")
        ax.tick_params(axis="x", rotation=30)
        ax.axhline(y=1.0, color="gray", linestyle="--", alpha=0.5)

    plt.tight_layout()
    plt.savefig(plots_dir / "sasa_boxplot_comparison.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("  Saved sasa_boxplot_comparison.png")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Evaluate Strategy W (Embedding-Based Interface Steering)"
    )
    parser.add_argument(
        "--step",
        choices=["index", "dockq", "epitope", "sasa", "summary", "all"],
        default="all",
        help="Which evaluation step to run (default: all)",
    )
    args = parser.parse_args()

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Step 1: Build index
    index_path = OUTPUT_DIR / "evaluation_index.csv"
    if args.step in ("index", "all"):
        index_df = build_index()
    elif index_path.exists():
        index_df = pd.read_csv(index_path)
        print(f"Loaded existing index from {index_path} ({len(index_df)} entries)")
    else:
        print("ERROR: No index found. Run with --step index first.")
        sys.exit(1)

    if index_df.empty:
        print("\nNo predictions found. Nothing to evaluate.")
        sys.exit(0)

    # Step 2: DockQ
    dockq_path = OUTPUT_DIR / "dockq_results.csv"
    dockq_df = pd.DataFrame()
    if args.step in ("dockq", "all"):
        dockq_df = run_dockq_all(index_df)
    elif dockq_path.exists():
        dockq_df = pd.read_csv(dockq_path)
        print(f"Loaded existing DockQ results from {dockq_path} ({len(dockq_df)} entries)")

    # Step 3: Epitope
    epitope_path = OUTPUT_DIR / "epitope_results.csv"
    epitope_df = pd.DataFrame()
    if args.step in ("epitope", "all"):
        epitope_df = run_epitope_all(index_df)
    elif epitope_path.exists():
        epitope_df = pd.read_csv(epitope_path)
        print(f"Loaded existing epitope results from {epitope_path}")

    # Step 4: SASA
    sasa_path = OUTPUT_DIR / "sasa_results.csv"
    sasa_df = pd.DataFrame()
    if args.step in ("sasa", "all"):
        sasa_df = run_sasa_all(index_df)
    elif sasa_path.exists():
        sasa_df = pd.read_csv(sasa_path)
        print(f"Loaded existing SASA results from {sasa_path}")

    # Step 5: Summary
    if args.step in ("summary", "all"):
        run_summary(dockq_df, epitope_df, sasa_df, index_df)

    print("\n" + "=" * 70)
    print("EVALUATION COMPLETE")
    print(f"Results in: {OUTPUT_DIR}/")
    print("=" * 70)


if __name__ == "__main__":
    main()
