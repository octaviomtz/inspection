#!/usr/bin/env python3
"""Evaluation script for Strategy D+ (Enhanced Canonical Ensemble).

Compares D+ against three baselines using DockQ (interface quality),
CDR-specific RMSD (loop accuracy), and confidence correlation.

Usage:
    conda run -n dockq2 python scripts/evaluate_canonical_ensemble.py \
        --predictions_dir predictions_examples \
        --gt_dir pdb_minimized \
        --cdrs_csv examples/cdrs.csv \
        --output_dir evaluation_results

See EVALUATION_PLAN.md for full methodology.
"""

import argparse
import ast
import csv
import itertools
import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.Align import PairwiseAligner
from Bio.Data.IUPACData import protein_letters_3to1
from Bio.PDB import MMCIFParser, PDBParser, Superimposer
from DockQ.DockQ import load_PDB, run_on_all_native_interfaces
from scipy.stats import spearmanr, wilcoxon

warnings.filterwarnings("ignore")

# ============================================================================
# Constants
# ============================================================================

THREE_TO_ONE = {k.upper(): v.upper() for k, v in protein_letters_3to1.items()}
STANDARD_AA = set(THREE_TO_ONE.keys())

CAPRI_THRESHOLDS = [
    (0.80, "High"),
    (0.49, "Medium"),
    (0.23, "Acceptable"),
    (0.00, "Incorrect"),
]

METHOD_LABELS = {
    "baseline1": "Boltz2 vanilla",
    "baseline2_best": "Contact restraints (best)",
    "baseline2_median": "Contact restraints (median)",
    "baseline3_best": "Pocket+MSA (best)",
    "baseline3_median": "Pocket+MSA (median)",
    "dplus": "D+ (ours)",
}

# ============================================================================
# File discovery
# ============================================================================


def get_complex_names(gt_dir):
    """Get list of complex names from ground truth directory."""
    return sorted(p.stem for p in Path(gt_dir).glob("*.pdb"))


def discover_predictions(predictions_dir, complex_name):
    """Discover all prediction files for a complex across all methods.

    Returns dict: method_key -> list of (structure_path, confidence_path, variant_name)
    """
    pdir = Path(predictions_dir)
    results = {
        "baseline1": [],
        "baseline2": [],
        "baseline3": [],
        "dplus": [],
    }

    # Baseline 1: antigen_cut/boltz_results_{complex}/predictions/{complex}/
    b1_dir = pdir / "antigen_cut"
    if b1_dir.exists():
        for run_dir in b1_dir.glob(f"boltz_results_{complex_name}"):
            pred_dir = run_dir / "predictions" / complex_name
            if pred_dir.exists():
                _collect_models(pred_dir, complex_name, "baseline1", results)

    # Baseline 2: antigen_cut_contact_restraints/boltz_results_restraint_{complex}_*/
    b2_dir = pdir / "antigen_cut_contact_restraints"
    if b2_dir.exists():
        for run_dir in b2_dir.glob(f"boltz_results_restraint_{complex_name}_*"):
            variant = run_dir.name.replace("boltz_results_", "")
            # predictions subdir name matches the variant
            pred_dir = run_dir / "predictions" / variant
            if pred_dir.exists():
                _collect_models(pred_dir, variant, "baseline2", results)

    # Baseline 3: antigen_cut_vhvl_msa_pocket_.../boltz_results_restraint_to_A_{complex}_*/
    b3_dir = pdir / "antigen_cut_vhvl_msa_pocket_ab_boltz_post_2023_omm"
    if b3_dir.exists():
        for run_dir in b3_dir.glob(f"boltz_results_restraint_to_A_{complex_name}_*"):
            variant = run_dir.name.replace("boltz_results_", "")
            pred_dir = run_dir / "predictions" / variant
            if pred_dir.exists():
                _collect_models(pred_dir, variant, "baseline3", results)

    # D+ (ours): new_feature_cdr3_beta/boltz_results_{complex}*/
    # Also check for canonical_ensemble folder pattern
    for dplus_dirname in [
        "new_feature_cdr3_beta",
        "canonical_ensemble",
        "test_canonical_ensemble",
    ]:
        dplus_dir = pdir / dplus_dirname
        if not dplus_dir.exists():
            continue
        for run_dir in dplus_dir.glob(f"boltz_results_{complex_name}*"):
            # Find the predictions subfolder (name varies)
            pred_parent = run_dir / "predictions"
            if not pred_parent.exists():
                continue
            for pred_dir in pred_parent.iterdir():
                if pred_dir.is_dir():
                    variant = pred_dir.name
                    _collect_models(pred_dir, variant, "dplus", results)

    return results


def _collect_models(pred_dir, variant, method_key, results):
    """Find all model files (PDB or CIF) and their confidence JSONs in a predictions dir."""
    for struct_file in sorted(pred_dir.glob("*_model_*.pdb")) + sorted(
        pred_dir.glob("*_model_*.cif")
    ):
        # Extract model index from filename like {name}_model_{idx}.pdb
        stem = struct_file.stem
        model_idx = int(stem.split("_model_")[-1])
        conf_file = pred_dir / f"confidence_{variant}_model_{model_idx}.json"
        plddt_file = pred_dir / f"plddt_{variant}_model_{model_idx}.npz"
        results[method_key].append(
            {
                "struct_path": str(struct_file),
                "conf_path": str(conf_file) if conf_file.exists() else None,
                "plddt_path": str(plddt_file) if plddt_file.exists() else None,
                "model_idx": model_idx,
                "variant": variant,
            }
        )


# ============================================================================
# DockQ evaluation
# ============================================================================


def run_dockq(model_path, native_path):
    """Run DockQ on a model-native pair.

    Returns dict with per-interface and global metrics, or None on failure.
    """
    try:
        model_struct = load_PDB(model_path)
        native_struct = load_PDB(native_path)

        # Determine chain IDs present in both
        model_chain_ids = {c.id for c in model_struct}
        native_chain_ids = {c.id for c in native_struct}
        common = sorted(model_chain_ids & native_chain_ids)

        if len(common) < 2:
            return None

        chain_map = {c: c for c in common}
        result_dict, total_dockq = run_on_all_native_interfaces(
            model_struct, native_struct, chain_map=chain_map
        )

        if not result_dict:
            return None

        global_dockq = total_dockq / len(result_dict)

        out = {"GlobalDockQ": global_dockq}
        for iface, metrics in result_dict.items():
            for key in ["DockQ", "fnat", "fnonnat", "iRMSD", "LRMSD", "F1"]:
                out[f"{key}_{iface}"] = metrics[key]

        return out
    except Exception as e:
        print(f"  DockQ error for {model_path}: {e}", file=sys.stderr)
        return None


# ============================================================================
# CDR RMSD evaluation
# ============================================================================


def parse_cdrs(cdrs_csv):
    """Parse CDR definitions from cdrs.csv.

    Returns dict: complex_name -> {cdr_name: list_of_1indexed_positions}
    """
    cdrs = {}
    with open(cdrs_csv) as f:
        reader = csv.DictReader(f)
        for row in reader:
            name = row["complex"]
            cdrs[name] = {
                "H1": ast.literal_eval(row["cdr1_h"]),
                "H2": ast.literal_eval(row["cdr2_h"]),
                "H3": ast.literal_eval(row["cdr3_h"]),
                "L1": ast.literal_eval(row["cdr1_l"]),
                "L2": ast.literal_eval(row["cdr2_l"]),
                "L3": ast.literal_eval(row["cdr3_l"]),
            }
    return cdrs


def _get_sequence_and_ca(structure, chain_id):
    """Extract sequence, CA atoms, and residue IDs for a chain."""
    seq, ca_atoms, res_ids = "", [], []
    for chain in structure[0]:
        if chain.id == chain_id:
            for res in chain:
                if res.id[0] == " " and res.get_resname() in STANDARD_AA:
                    seq += THREE_TO_ONE[res.get_resname()]
                    res_ids.append(res.id[1])
                    if "CA" in res:
                        ca_atoms.append(res["CA"])
                    else:
                        ca_atoms.append(None)
    return seq, ca_atoms, res_ids


def _align_sequences(seq_native, seq_model):
    """Pairwise sequence alignment. Returns list of (native_idx, model_idx) pairs."""
    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -0.5

    alignments = aligner.align(seq_native, seq_model)
    if not alignments:
        return []

    aln = alignments[0]
    pairs = []
    for (ns, ne), (ms, me) in zip(aln.aligned[0], aln.aligned[1]):
        for ni, mi in zip(range(ns, ne), range(ms, me)):
            pairs.append((ni, mi))
    return pairs


def compute_cdr_rmsd(model_path, native_path, cdr_defs):
    """Compute framework-aligned CDR RMSD for a model-native pair.

    cdr_defs: {cdr_name: list_of_1indexed_positions}
        H1/H2/H3 are positions in chain B (model numbering, 1-indexed)
        L1/L2/L3 are positions in chain C (model numbering, 1-indexed)

    Returns dict: {cdr_name: rmsd_value, "all_cdr": combined_rmsd} or None on failure.
    """
    pdb_parser = PDBParser(QUIET=True)
    cif_parser = MMCIFParser(QUIET=True)

    try:
        if str(model_path).endswith(".cif"):
            model_struct = cif_parser.get_structure("m", str(model_path))
        else:
            model_struct = pdb_parser.get_structure("m", str(model_path))
        native_struct = pdb_parser.get_structure("n", str(native_path))
    except Exception:
        return None

    cdr_chain_map = {"H1": "B", "H2": "B", "H3": "B", "L1": "C", "L2": "C", "L3": "C"}

    # For each antibody chain, align sequences and classify residues
    framework_atoms_native = []
    framework_atoms_model = []
    cdr_atoms = {name: {"native": [], "model": []} for name in cdr_defs}

    for chain_id in ["B", "C"]:
        n_seq, n_ca, n_resids = _get_sequence_and_ca(native_struct, chain_id)
        m_seq, m_ca, m_resids = _get_sequence_and_ca(model_struct, chain_id)

        if not n_seq or not m_seq:
            return None

        pairs = _align_sequences(n_seq, m_seq)
        if not pairs:
            return None

        # Determine which model residue indices (0-based in sequence) are CDR
        cdr_set_for_chain = set()
        for cdr_name, positions in cdr_defs.items():
            if cdr_chain_map[cdr_name] != chain_id:
                continue
            # positions are 1-indexed model residue numbers
            for pos in positions:
                # Convert to 0-based model sequence index
                if pos - 1 < len(m_resids):
                    # Find the sequence index where model residue ID == pos
                    try:
                        seq_idx = m_resids.index(pos)
                        cdr_set_for_chain.add((seq_idx, cdr_name))
                    except ValueError:
                        pass

        for ni, mi in pairs:
            if n_ca[ni] is None or m_ca[mi] is None:
                continue

            # Check if this model residue is in any CDR
            cdr_hit = None
            for seq_idx, cdr_name in cdr_set_for_chain:
                if mi == seq_idx:
                    cdr_hit = cdr_name
                    break

            if cdr_hit:
                cdr_atoms[cdr_hit]["native"].append(n_ca[ni])
                cdr_atoms[cdr_hit]["model"].append(m_ca[mi])
            else:
                framework_atoms_native.append(n_ca[ni])
                framework_atoms_model.append(m_ca[mi])

    if len(framework_atoms_native) < 10:
        return None

    # Superimpose on framework
    sup = Superimposer()
    try:
        sup.set_atoms(framework_atoms_native, framework_atoms_model)
    except Exception:
        return None

    # Apply transformation to ALL model atoms
    sup.apply(model_struct[0].get_atoms())

    # Compute per-CDR RMSD
    results = {}
    all_diffs = []
    for cdr_name, atoms in cdr_atoms.items():
        if not atoms["native"] or not atoms["model"]:
            results[cdr_name] = np.nan
            continue
        # After superimposition, get updated model coordinates
        coords_n = np.array([a.get_vector().get_array() for a in atoms["native"]])
        coords_m = np.array([a.get_vector().get_array() for a in atoms["model"]])
        diffs = np.sum((coords_n - coords_m) ** 2, axis=1)
        results[cdr_name] = np.sqrt(np.mean(diffs))
        all_diffs.extend(diffs.tolist())

    if all_diffs:
        results["all_cdr"] = np.sqrt(np.mean(all_diffs))
    else:
        results["all_cdr"] = np.nan

    results["framework_rmsd"] = sup.rms
    return results


# ============================================================================
# Confidence reading
# ============================================================================


def read_confidence(conf_path):
    """Read confidence JSON and return key metrics."""
    if conf_path is None or not Path(conf_path).exists():
        return None
    try:
        with open(conf_path) as f:
            data = json.load(f)
        return {
            "confidence_score": data.get("confidence_score"),
            "complex_plddt": data.get("complex_plddt"),
            "iptm": data.get("iptm"),
            "ptm": data.get("ptm"),
        }
    except Exception:
        return None


# ============================================================================
# Ensemble diversity
# ============================================================================


def compute_ensemble_diversity(model_paths, native_path, cdr_defs):
    """Compute pairwise CDR-H3 RMSD across an ensemble of models.

    Returns dict with mean, std, min, max pairwise RMSD, or None.
    """
    if len(model_paths) < 2:
        return None

    pdb_parser = PDBParser(QUIET=True)
    cif_parser = MMCIFParser(QUIET=True)
    native_struct = pdb_parser.get_structure("n", str(native_path))

    # For each model: framework-align to native, then extract CDR-H3 CA coords
    h3_coords_list = []

    h3_positions = set(cdr_defs.get("H3", []))
    if not h3_positions:
        return None

    for model_path in model_paths:
        try:
            if str(model_path).endswith(".cif"):
                model_struct = cif_parser.get_structure("m", str(model_path))
            else:
                model_struct = pdb_parser.get_structure("m", str(model_path))

            # Get sequences for chain B
            n_seq, n_ca, n_resids = _get_sequence_and_ca(native_struct, "B")
            m_seq, m_ca, m_resids = _get_sequence_and_ca(model_struct, "B")

            if not n_seq or not m_seq:
                continue

            pairs = _align_sequences(n_seq, m_seq)
            if not pairs:
                continue

            # Also need chain C for framework alignment
            n_seq_c, n_ca_c, _ = _get_sequence_and_ca(native_struct, "C")
            m_seq_c, m_ca_c, m_resids_c = _get_sequence_and_ca(model_struct, "C")
            pairs_c = _align_sequences(n_seq_c, m_seq_c) if n_seq_c and m_seq_c else []

            # Classify framework vs CDR for chain B
            h3_model_seq_indices = set()
            for pos in h3_positions:
                if pos in m_resids:
                    h3_model_seq_indices.add(m_resids.index(pos))

            fw_native, fw_model = [], []
            h3_native, h3_model = [], []
            for ni, mi in pairs:
                if n_ca[ni] is None or m_ca[mi] is None:
                    continue
                if mi in h3_model_seq_indices:
                    h3_native.append(n_ca[ni])
                    h3_model.append(m_ca[mi])
                else:
                    fw_native.append(n_ca[ni])
                    fw_model.append(m_ca[mi])

            # Add chain C framework
            for ni, mi in pairs_c:
                if n_ca_c[ni] is None or m_ca_c[mi] is None:
                    continue
                fw_native.append(n_ca_c[ni])
                fw_model.append(m_ca_c[mi])

            if len(fw_native) < 10 or not h3_model:
                continue

            sup = Superimposer()
            sup.set_atoms(fw_native, fw_model)
            sup.apply(model_struct[0].get_atoms())

            coords = np.array([a.get_vector().get_array() for a in h3_model])
            h3_coords_list.append(coords)
        except Exception:
            continue

    if len(h3_coords_list) < 2:
        return None

    # Pairwise RMSD between all models
    # Coords must have same length (same CDR-H3 residues aligned)
    min_len = min(len(c) for c in h3_coords_list)
    h3_coords_list = [c[:min_len] for c in h3_coords_list]

    pairwise_rmsds = []
    for i, j in itertools.combinations(range(len(h3_coords_list)), 2):
        diff = h3_coords_list[i] - h3_coords_list[j]
        rmsd = np.sqrt(np.mean(np.sum(diff**2, axis=1)))
        pairwise_rmsds.append(rmsd)

    pairwise_rmsds = np.array(pairwise_rmsds)
    return {
        "mean": float(np.mean(pairwise_rmsds)),
        "std": float(np.std(pairwise_rmsds)),
        "min": float(np.min(pairwise_rmsds)),
        "max": float(np.max(pairwise_rmsds)),
        "n_models": len(h3_coords_list),
    }


# ============================================================================
# Aggregation helpers
# ============================================================================


def classify_capri(dockq):
    """Classify DockQ into CAPRI quality category."""
    for threshold, label in CAPRI_THRESHOLDS:
        if dockq >= threshold:
            return label
    return "Incorrect"


def aggregate_baselines_multi_variant(rows):
    """For baselines with multiple variants per complex, aggregate.

    Returns (best_model0_row, oracle_row, median_model0_value).
    """
    if not rows:
        return None, None, None

    # best model_0: among model_idx==0, pick highest GlobalDockQ
    model0_rows = [r for r in rows if r["model_idx"] == 0]
    if not model0_rows:
        model0_rows = rows  # fallback: use all

    best_model0 = max(model0_rows, key=lambda r: r.get("GlobalDockQ", -1))
    median_val = np.median([r.get("GlobalDockQ", 0) for r in model0_rows])

    # oracle: best across all variants and models
    oracle = max(rows, key=lambda r: r.get("GlobalDockQ", -1))

    return best_model0, oracle, float(median_val)


# ============================================================================
# Main evaluation
# ============================================================================


def evaluate(predictions_dir, gt_dir, cdrs_csv, output_dir):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    complex_names = get_complex_names(gt_dir)
    cdr_definitions = parse_cdrs(cdrs_csv)

    print(f"Found {len(complex_names)} complexes in {gt_dir}")
    print(f"CDR definitions for {len(cdr_definitions)} complexes")

    # ---- Collect all per-model results ----
    all_rows = []

    for ci, cname in enumerate(complex_names):
        print(f"\n[{ci+1}/{len(complex_names)}] Processing {cname}...")
        native_path = Path(gt_dir) / f"{cname}.pdb"
        if not native_path.exists():
            print(f"  Ground truth not found: {native_path}")
            continue

        predictions = discover_predictions(predictions_dir, cname)
        cdr_def = cdr_definitions.get(cname)

        for method_key in ["baseline1", "baseline2", "baseline3", "dplus"]:
            models = predictions[method_key]
            if not models:
                continue

            print(f"  {method_key}: {len(models)} model files")

            for model_info in models:
                struct_path = model_info["struct_path"]
                row = {
                    "complex": cname,
                    "method": method_key,
                    "model_idx": model_info["model_idx"],
                    "variant": model_info["variant"],
                    "struct_path": struct_path,
                }

                # DockQ
                dockq_result = run_dockq(struct_path, str(native_path))
                if dockq_result:
                    row.update(dockq_result)

                # CDR RMSD
                if cdr_def:
                    cdr_result = compute_cdr_rmsd(
                        struct_path, str(native_path), cdr_def
                    )
                    if cdr_result:
                        for k, v in cdr_result.items():
                            row[f"cdr_{k}"] = v

                # Confidence
                conf = read_confidence(model_info["conf_path"])
                if conf:
                    row.update({f"pred_{k}": v for k, v in conf.items()})

                all_rows.append(row)

    if not all_rows:
        print("\nNo results collected. Check that prediction directories have model files.")
        return

    df = pd.DataFrame(all_rows)

    # Save raw results
    raw_path = output_dir / "raw_results.csv"
    df.to_csv(raw_path, index=False)
    print(f"\nRaw results saved to {raw_path} ({len(df)} rows)")

    # ---- Table 1: DockQ per-interface (best model_0) ----
    print("\n" + "=" * 80)
    print("TABLE 1: DockQ Per-Interface (Best Model)")
    print("=" * 80)
    table1_rows = []
    for cname in complex_names:
        for method_key in ["baseline1", "dplus", "baseline2", "baseline3"]:
            subset = df[(df["complex"] == cname) & (df["method"] == method_key)]
            if subset.empty:
                continue

            if method_key in ("baseline2", "baseline3"):
                # Multi-variant: pick best model_0
                m0 = subset[subset["model_idx"] == 0]
                if m0.empty:
                    m0 = subset
                best = m0.loc[m0["GlobalDockQ"].idxmax()] if "GlobalDockQ" in m0 else m0.iloc[0]
            else:
                # Single variant: pick model_0
                m0 = subset[subset["model_idx"] == 0]
                best = m0.iloc[0] if not m0.empty else subset.iloc[0]

            t1row = {
                "complex": cname,
                "method": method_key,
            }
            for col in [
                "GlobalDockQ",
                "DockQ_AB",
                "DockQ_AC",
                "DockQ_BC",
                "fnat_AB",
                "iRMSD_AB",
                "LRMSD_AB",
            ]:
                t1row[col] = best.get(col, np.nan)
            table1_rows.append(t1row)

    df_t1 = pd.DataFrame(table1_rows)
    df_t1.to_csv(output_dir / "table1_dockq_per_interface.csv", index=False)
    print(df_t1.to_string(index=False, float_format="%.4f"))

    # ---- Table 2: CDR RMSD (best model_0) ----
    print("\n" + "=" * 80)
    print("TABLE 2: CDR RMSD (Best Model, Framework-Aligned)")
    print("=" * 80)
    table2_rows = []
    for cname in complex_names:
        for method_key in ["baseline1", "dplus", "baseline2", "baseline3"]:
            subset = df[(df["complex"] == cname) & (df["method"] == method_key)]
            if subset.empty:
                continue

            if method_key in ("baseline2", "baseline3"):
                m0 = subset[subset["model_idx"] == 0]
                if m0.empty:
                    m0 = subset
                best = m0.loc[m0["GlobalDockQ"].idxmax()] if "GlobalDockQ" in m0 else m0.iloc[0]
            else:
                m0 = subset[subset["model_idx"] == 0]
                best = m0.iloc[0] if not m0.empty else subset.iloc[0]

            t2row = {"complex": cname, "method": method_key}
            for cdr in ["H1", "H2", "H3", "L1", "L2", "L3", "all_cdr"]:
                t2row[cdr] = best.get(f"cdr_{cdr}", np.nan)
            t2row["framework"] = best.get("cdr_framework_rmsd", np.nan)
            table2_rows.append(t2row)

    df_t2 = pd.DataFrame(table2_rows)
    df_t2.to_csv(output_dir / "table2_cdr_rmsd.csv", index=False)
    print(df_t2.to_string(index=False, float_format="%.3f"))

    # ---- Table 3: Aggregated results (mean ± std) ----
    print("\n" + "=" * 80)
    print("TABLE 3: Aggregated Results (Mean ± Std)")
    print("=" * 80)
    table3_rows = []

    for method_key in ["baseline1", "dplus", "baseline2", "baseline3"]:
        # For each complex, get the best model_0 row
        per_complex_best = []
        for cname in complex_names:
            subset = df[(df["complex"] == cname) & (df["method"] == method_key)]
            if subset.empty:
                continue
            if method_key in ("baseline2", "baseline3"):
                m0 = subset[subset["model_idx"] == 0]
                if m0.empty:
                    m0 = subset
                if "GlobalDockQ" in m0.columns and m0["GlobalDockQ"].notna().any():
                    best = m0.loc[m0["GlobalDockQ"].idxmax()]
                else:
                    best = m0.iloc[0]
            else:
                m0 = subset[subset["model_idx"] == 0]
                best = m0.iloc[0] if not m0.empty else subset.iloc[0]
            per_complex_best.append(best)

        if not per_complex_best:
            continue

        df_best = pd.DataFrame(per_complex_best)
        label = METHOD_LABELS.get(
            method_key + ("_best" if method_key.startswith("baseline") and method_key != "baseline1" else ""),
            method_key,
        )

        t3row = {"method": label, "n_complexes": len(df_best)}
        for col, name in [
            ("GlobalDockQ", "DockQ_Global"),
            ("DockQ_AB", "DockQ_AB"),
            ("DockQ_AC", "DockQ_AC"),
            ("DockQ_BC", "DockQ_BC"),
            ("cdr_H3", "CDR-H3"),
            ("cdr_all_cdr", "All-CDR"),
        ]:
            if col in df_best.columns:
                vals = df_best[col].dropna()
                if len(vals) > 0:
                    t3row[name] = f"{vals.mean():.4f} ± {vals.std():.4f}"
                else:
                    t3row[name] = "N/A"
            else:
                t3row[name] = "N/A"

        table3_rows.append(t3row)

    df_t3 = pd.DataFrame(table3_rows)
    df_t3.to_csv(output_dir / "table3_aggregated.csv", index=False)
    print(df_t3.to_string(index=False))

    # ---- Table 3b: Oracle results (best across all models) ----
    print("\n" + "=" * 80)
    print("TABLE 3b: Oracle Results (Best Across All Models)")
    print("=" * 80)
    table3b_rows = []
    for method_key in ["baseline1", "dplus", "baseline2", "baseline3"]:
        per_complex_oracle = []
        for cname in complex_names:
            subset = df[(df["complex"] == cname) & (df["method"] == method_key)]
            if subset.empty:
                continue
            if "GlobalDockQ" in subset.columns and subset["GlobalDockQ"].notna().any():
                oracle = subset.loc[subset["GlobalDockQ"].idxmax()]
            else:
                oracle = subset.iloc[0]
            per_complex_oracle.append(oracle)

        if not per_complex_oracle:
            continue

        df_oracle = pd.DataFrame(per_complex_oracle)
        t3b_row = {"method": method_key, "n_complexes": len(df_oracle)}
        for col, name in [
            ("GlobalDockQ", "DockQ_Global"),
            ("DockQ_AB", "DockQ_AB"),
            ("cdr_H3", "CDR-H3"),
        ]:
            if col in df_oracle.columns:
                vals = df_oracle[col].dropna()
                if len(vals) > 0:
                    t3b_row[name] = f"{vals.mean():.4f} ± {vals.std():.4f}"
                else:
                    t3b_row[name] = "N/A"
            else:
                t3b_row[name] = "N/A"
        table3b_rows.append(t3b_row)

    df_t3b = pd.DataFrame(table3b_rows)
    df_t3b.to_csv(output_dir / "table3b_oracle.csv", index=False)
    print(df_t3b.to_string(index=False))

    # ---- Table 4: CAPRI quality classification ----
    print("\n" + "=" * 80)
    print("TABLE 4: CAPRI Quality Classification (% of complexes, best model_0)")
    print("=" * 80)
    table4_rows = []
    for method_key in ["baseline1", "dplus", "baseline2", "baseline3"]:
        per_complex = []
        for cname in complex_names:
            subset = df[(df["complex"] == cname) & (df["method"] == method_key)]
            if subset.empty:
                continue
            if method_key in ("baseline2", "baseline3"):
                m0 = subset[subset["model_idx"] == 0]
                if m0.empty:
                    m0 = subset
                if "GlobalDockQ" in m0.columns and m0["GlobalDockQ"].notna().any():
                    best = m0.loc[m0["GlobalDockQ"].idxmax()]
                else:
                    continue
            else:
                m0 = subset[subset["model_idx"] == 0]
                if m0.empty:
                    continue
                best = m0.iloc[0]
            if pd.notna(best.get("GlobalDockQ")):
                per_complex.append(best["GlobalDockQ"])

        if not per_complex:
            continue

        n = len(per_complex)
        categories = [classify_capri(d) for d in per_complex]
        t4row = {"method": method_key, "n_complexes": n}
        for label in ["High", "Medium", "Acceptable", "Incorrect"]:
            count = categories.count(label)
            t4row[label] = f"{count} ({100*count/n:.1f}%)"
        table4_rows.append(t4row)

    df_t4 = pd.DataFrame(table4_rows)
    df_t4.to_csv(output_dir / "table4_capri.csv", index=False)
    print(df_t4.to_string(index=False))

    # ---- Table 5: Ensemble diversity (D+ and Baseline 1) ----
    print("\n" + "=" * 80)
    print("TABLE 5: Ensemble Diversity (CDR-H3 Pairwise RMSD)")
    print("=" * 80)
    table5_rows = []
    for cname in complex_names:
        native_path = Path(gt_dir) / f"{cname}.pdb"
        cdr_def = cdr_definitions.get(cname)
        if not cdr_def:
            continue

        predictions = discover_predictions(predictions_dir, cname)

        for method_key, label in [("dplus", "D+"), ("baseline1", "Baseline1")]:
            models = predictions[method_key]
            if len(models) < 2:
                continue

            model_paths = [m["struct_path"] for m in sorted(models, key=lambda x: x["model_idx"])]
            diversity = compute_ensemble_diversity(
                model_paths, str(native_path), cdr_def
            )
            if diversity:
                table5_rows.append(
                    {
                        "complex": cname,
                        "method": label,
                        "n_models": diversity["n_models"],
                        "mean_pairwise_H3_rmsd": diversity["mean"],
                        "std": diversity["std"],
                        "min": diversity["min"],
                        "max": diversity["max"],
                    }
                )

    if table5_rows:
        df_t5 = pd.DataFrame(table5_rows)
        df_t5.to_csv(output_dir / "table5_diversity.csv", index=False)
        print(df_t5.to_string(index=False, float_format="%.3f"))
    else:
        print("  (No ensembles with >=2 models found)")

    # ---- Table 6: Confidence vs DockQ correlation ----
    print("\n" + "=" * 80)
    print("TABLE 6: Confidence vs Actual DockQ Correlation")
    print("=" * 80)
    table6_rows = []
    for method_key in ["baseline1", "dplus", "baseline2", "baseline3"]:
        subset = df[df["method"] == method_key].dropna(
            subset=["GlobalDockQ", "pred_confidence_score"]
        )
        if len(subset) < 5:
            continue

        rho, pval = spearmanr(subset["pred_confidence_score"], subset["GlobalDockQ"])
        table6_rows.append(
            {
                "method": method_key,
                "n_samples": len(subset),
                "spearman_rho": f"{rho:.4f}",
                "p_value": f"{pval:.4e}",
            }
        )

        # Also test iptm and complex_plddt
        for pred_col, name in [
            ("pred_iptm", "iptm"),
            ("pred_complex_plddt", "complex_plddt"),
        ]:
            sub2 = df[df["method"] == method_key].dropna(
                subset=["GlobalDockQ", pred_col]
            )
            if len(sub2) >= 5:
                rho2, pval2 = spearmanr(sub2[pred_col], sub2["GlobalDockQ"])
                table6_rows.append(
                    {
                        "method": f"{method_key} ({name})",
                        "n_samples": len(sub2),
                        "spearman_rho": f"{rho2:.4f}",
                        "p_value": f"{pval2:.4e}",
                    }
                )

    if table6_rows:
        df_t6 = pd.DataFrame(table6_rows)
        df_t6.to_csv(output_dir / "table6_confidence_correlation.csv", index=False)
        print(df_t6.to_string(index=False))
    else:
        print("  (Not enough data for correlation)")

    # ---- Statistical tests: D+ vs baselines ----
    print("\n" + "=" * 80)
    print("STATISTICAL TESTS: Paired Wilcoxon Signed-Rank (D+ vs Baselines)")
    print("=" * 80)
    stat_rows = []

    # Build per-complex best-model_0 tables for each method
    method_dockq = {}
    for method_key in ["baseline1", "dplus", "baseline2", "baseline3"]:
        per_complex = {}
        for cname in complex_names:
            subset = df[(df["complex"] == cname) & (df["method"] == method_key)]
            if subset.empty:
                continue
            if method_key in ("baseline2", "baseline3"):
                m0 = subset[subset["model_idx"] == 0]
                if m0.empty:
                    m0 = subset
                if "GlobalDockQ" in m0.columns and m0["GlobalDockQ"].notna().any():
                    per_complex[cname] = m0["GlobalDockQ"].max()
            else:
                m0 = subset[subset["model_idx"] == 0]
                if not m0.empty and pd.notna(m0.iloc[0].get("GlobalDockQ")):
                    per_complex[cname] = m0.iloc[0]["GlobalDockQ"]
        method_dockq[method_key] = per_complex

    dplus_scores = method_dockq.get("dplus", {})
    for baseline_key in ["baseline1", "baseline2", "baseline3"]:
        baseline_scores = method_dockq.get(baseline_key, {})
        common = sorted(set(dplus_scores.keys()) & set(baseline_scores.keys()))
        if len(common) < 5:
            print(f"  D+ vs {baseline_key}: <5 common complexes ({len(common)}), skipping")
            continue

        d_vals = [dplus_scores[c] for c in common]
        b_vals = [baseline_scores[c] for c in common]
        diff = np.array(d_vals) - np.array(b_vals)

        try:
            stat, pval = wilcoxon(diff, alternative="two-sided")
        except ValueError:
            stat, pval = np.nan, np.nan

        stat_row = {
            "comparison": f"D+ vs {baseline_key}",
            "n_common": len(common),
            "D+_mean": f"{np.mean(d_vals):.4f}",
            "baseline_mean": f"{np.mean(b_vals):.4f}",
            "mean_diff": f"{np.mean(diff):.4f}",
            "wilcoxon_stat": f"{stat:.1f}",
            "p_value": f"{pval:.4e}",
            "D+_wins": int(np.sum(diff > 0)),
            "ties": int(np.sum(diff == 0)),
            "D+_loses": int(np.sum(diff < 0)),
        }
        stat_rows.append(stat_row)
        print(
            f"  D+ vs {baseline_key}: mean_diff={np.mean(diff):.4f}, "
            f"p={pval:.4e}, wins/ties/losses={stat_row['D+_wins']}/{stat_row['ties']}/{stat_row['D+_loses']}"
        )

    if stat_rows:
        df_stat = pd.DataFrame(stat_rows)
        df_stat.to_csv(output_dir / "statistical_tests.csv", index=False)

    # ---- Generate plots ----
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        generate_plots(df, complex_names, method_dockq, output_dir)
    except ImportError:
        print("\nmatplotlib not available, skipping plots")

    print(f"\nAll results saved to {output_dir}/")


# ============================================================================
# Plotting
# ============================================================================


def generate_plots(df, complex_names, method_dockq, output_dir):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # Plot 1: DockQ comparison bar chart
    methods_present = [m for m in ["baseline1", "dplus", "baseline2", "baseline3"] if method_dockq.get(m)]
    if len(methods_present) >= 2:
        fig, ax = plt.subplots(figsize=(8, 5))
        means = []
        stds = []
        labels = []
        for m in methods_present:
            vals = list(method_dockq[m].values())
            means.append(np.mean(vals))
            stds.append(np.std(vals))
            labels.append(METHOD_LABELS.get(m + ("_best" if m.startswith("baseline") and m != "baseline1" else ""), m))

        x = np.arange(len(labels))
        bars = ax.bar(x, means, yerr=stds, capsize=5, color=["#4C72B0", "#DD8452", "#55A868", "#C44E52"][: len(labels)])
        ax.set_ylabel("Global DockQ")
        ax.set_title("Global DockQ by Method (best model_0)")
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=15, ha="right")
        ax.set_ylim(0, 1)
        for bar, m in zip(bars, means):
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.02, f"{m:.3f}", ha="center", fontsize=9)
        plt.tight_layout()
        plt.savefig(output_dir / "plot1_dockq_comparison.png", dpi=150)
        plt.close()
        print(f"  Saved plot1_dockq_comparison.png")

    # Plot 2: Per-complex DockQ comparison (D+ vs Baseline 1)
    dplus_scores = method_dockq.get("dplus", {})
    b1_scores = method_dockq.get("baseline1", {})
    common = sorted(set(dplus_scores.keys()) & set(b1_scores.keys()))
    if len(common) >= 3:
        fig, ax = plt.subplots(figsize=(10, 5))
        x = np.arange(len(common))
        width = 0.35
        d_vals = [dplus_scores[c] for c in common]
        b_vals = [b1_scores[c] for c in common]
        ax.bar(x - width / 2, b_vals, width, label="Boltz2 vanilla", color="#4C72B0")
        ax.bar(x + width / 2, d_vals, width, label="D+ (ours)", color="#DD8452")
        ax.set_ylabel("Global DockQ")
        ax.set_title("Per-Complex DockQ: D+ vs Boltz2 Vanilla")
        ax.set_xticks(x)
        ax.set_xticklabels(common, rotation=90, fontsize=7)
        ax.legend()
        ax.set_ylim(0, 1)
        plt.tight_layout()
        plt.savefig(output_dir / "plot2_per_complex_dockq.png", dpi=150)
        plt.close()
        print(f"  Saved plot2_per_complex_dockq.png")

    # Plot 3: CDR-H3 RMSD comparison
    cdr_cols = ["cdr_H3"]
    methods_with_cdr = []
    for m in ["baseline1", "dplus"]:
        subset = df[df["method"] == m]
        if not subset.empty and "cdr_H3" in subset.columns and subset["cdr_H3"].notna().any():
            methods_with_cdr.append(m)

    if len(methods_with_cdr) >= 2:
        fig, ax = plt.subplots(figsize=(8, 5))
        data_to_plot = []
        plot_labels = []
        for m in methods_with_cdr:
            # Get best model_0 per complex
            vals = []
            for cname in complex_names:
                subset = df[(df["complex"] == cname) & (df["method"] == m) & (df["model_idx"] == 0)]
                if not subset.empty and pd.notna(subset.iloc[0].get("cdr_H3")):
                    vals.append(subset.iloc[0]["cdr_H3"])
            if vals:
                data_to_plot.append(vals)
                plot_labels.append(METHOD_LABELS.get(m, m))

        if data_to_plot:
            ax.boxplot(data_to_plot, labels=plot_labels)
            ax.set_ylabel("CDR-H3 RMSD (Å)")
            ax.set_title("CDR-H3 RMSD Distribution (framework-aligned)")
            plt.tight_layout()
            plt.savefig(output_dir / "plot3_cdr_h3_rmsd.png", dpi=150)
            plt.close()
            print(f"  Saved plot3_cdr_h3_rmsd.png")

    # Plot 4: Confidence vs DockQ scatter
    for method_key, color in [("baseline1", "#4C72B0"), ("dplus", "#DD8452")]:
        subset = df[df["method"] == method_key].dropna(subset=["GlobalDockQ", "pred_confidence_score"])
        if len(subset) < 3:
            continue
        fig, ax = plt.subplots(figsize=(6, 5))
        ax.scatter(subset["pred_confidence_score"], subset["GlobalDockQ"], alpha=0.7, color=color)
        ax.set_xlabel("Predicted Confidence Score")
        ax.set_ylabel("Global DockQ")
        ax.set_title(f"Confidence vs DockQ ({method_key})")
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        # Add diagonal
        ax.plot([0, 1], [0, 1], "k--", alpha=0.3)
        rho, pval = spearmanr(subset["pred_confidence_score"], subset["GlobalDockQ"])
        ax.text(0.05, 0.95, f"ρ={rho:.3f}, p={pval:.3e}", transform=ax.transAxes, fontsize=9, va="top")
        plt.tight_layout()
        plt.savefig(output_dir / f"plot4_confidence_vs_dockq_{method_key}.png", dpi=150)
        plt.close()
        print(f"  Saved plot4_confidence_vs_dockq_{method_key}.png")


# ============================================================================
# Entry point
# ============================================================================


def main():
    parser = argparse.ArgumentParser(
        description="Evaluate Strategy D+ (Enhanced Canonical Ensemble)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # With sample data already present
  conda run -n dockq2 python scripts/evaluate_canonical_ensemble.py

  # Custom paths
  conda run -n dockq2 python scripts/evaluate_canonical_ensemble.py \\
      --predictions_dir /path/to/predictions \\
      --gt_dir /path/to/pdb_minimized \\
      --cdrs_csv /path/to/cdrs.csv \\
      --output_dir /path/to/output
        """,
    )
    parser.add_argument(
        "--predictions_dir",
        default="predictions_examples",
        help="Directory containing method subdirectories (default: predictions_examples)",
    )
    parser.add_argument(
        "--gt_dir",
        default="pdb_minimized",
        help="Directory containing ground truth PDB files (default: pdb_minimized)",
    )
    parser.add_argument(
        "--cdrs_csv",
        default="examples/cdrs.csv",
        help="Path to CDR definitions CSV (default: examples/cdrs.csv)",
    )
    parser.add_argument(
        "--output_dir",
        default="evaluation_results",
        help="Output directory for results (default: evaluation_results)",
    )
    args = parser.parse_args()
    evaluate(args.predictions_dir, args.gt_dir, args.cdrs_csv, args.output_dir)


if __name__ == "__main__":
    main()
