#!/usr/bin/env python
"""
Evaluation script for Strategy Q+ v2 (Iterative Epitope Refinement v2) vs baselines.

Computes:
  - Axis A (Primary): Epitope prediction accuracy (Precision, Recall, F1, MCC)
  - Axis B (Secondary): Docking/interface quality (DockQ, iRMSD, LRMSD, fnat)
  - Axis C (Supplementary): Antibody-aligned antigen RMSD
  - Axis D (Supplementary): Per-residue confidence (CDR3 pLDDT, PAE)
  - Statistical comparisons (Wilcoxon signed-rank) vs B1 and Q v1
  - Summary tables (CSV) and plots (PNG)

Usage:
    conda activate dockq2
    python evaluation/evaluate_q_plus.py \
        --pred_dir predictions_examples \
        --native_dir pdb_minimized \
        --cdr_csv examples/cdrs.csv \
        --out_dir evaluation/results_q_plus

    # With Q v1 comparison:
    python evaluation/evaluate_q_plus.py \
        --pred_dir predictions_examples \
        --native_dir pdb_minimized \
        --cdr_csv examples/cdrs.csv \
        --qv1_selected ../../worktree_steering_cdr3/q_epitope_refinement/evaluation/results_q/selected_models_aggregated.csv \
        --out_dir evaluation/results_q_plus

See EVALUATION_STRATEGY.md for full details.
"""

import argparse
import ast
import csv
import json
import logging
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.PDB import PDBParser, NeighborSearch
from DockQ.DockQ import load_PDB, run_on_all_native_interfaces
from scipy import stats

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

CHAIN_MAP = {"A": "A", "B": "B", "C": "C"}

# Distance threshold for epitope extraction (Angstroms)
EPITOPE_THRESHOLD = 5.0
EPITOPE_THRESHOLD_PERMISSIVE = 8.0

# Method definitions: folder names and discovery patterns
METHOD_DEFS = {
    "B1": {
        "folder": "antigen_cut",
        "prefix": "boltz_results_",
        "multi_config": False,
    },
    "B2": {
        "folder": "antigen_cut_contact_restraints",
        "prefix": "boltz_results_restraint_",
        "multi_config": True,
    },
    "B3": {
        "folder": "antigen_cut_vhvl_msa_pocket_ab_boltz_post_2023_omm",
        "prefix": "boltz_results_restraint_to_A_",
        "multi_config": True,
    },
    "Q+": {
        "folder": "new_feature_epitope_refinement",
        "prefix": "boltz_results_",
        "multi_config": False,
    },
}


# ---------------------------------------------------------------------------
# CDR loading
# ---------------------------------------------------------------------------


def load_cdrs(cdrs_csv: Path) -> dict:
    """Load CDR residue indices from cdrs.csv.

    Returns dict: complex_name -> {
        'cdr1_h': [0-indexed positions], 'cdr2_h': ..., 'cdr3_h': ...,
        'cdr1_l': ..., 'cdr2_l': ..., 'cdr3_l': ...,
        'heavy_len': int, 'light_len': int,
    }
    """
    cdrs = {}
    with open(cdrs_csv) as f:
        reader = csv.DictReader(f)
        for row in reader:
            cname = row["complex"]
            cdrs[cname] = {
                "cdr1_h": ast.literal_eval(row["cdr1_h"]),
                "cdr2_h": ast.literal_eval(row["cdr2_h"]),
                "cdr3_h": ast.literal_eval(row["cdr3_h"]),
                "cdr1_l": ast.literal_eval(row["cdr1_l"]),
                "cdr2_l": ast.literal_eval(row["cdr2_l"]),
                "cdr3_l": ast.literal_eval(row["cdr3_l"]),
                "heavy_len": len(row["heavy"]),
                "light_len": len(row["light"]),
            }
    return cdrs


# ---------------------------------------------------------------------------
# Discovery: find prediction files and map to ground truths
# ---------------------------------------------------------------------------


def get_known_complexes(native_dir: Path) -> list[str]:
    """Get list of complex names from ground truth PDB files."""
    return [pdb.stem for pdb in sorted(native_dir.glob("*.pdb"))]


def extract_complex_name(folder_name: str, known_complexes: list[str]) -> str | None:
    """Extract complex name from a boltz_results folder name by matching
    against the known complex list (longest match first)."""
    for cname in sorted(known_complexes, key=len, reverse=True):
        if cname in folder_name:
            return cname
    return None


def extract_sub_experiment(folder_name: str, complex_name: str, method: str) -> str:
    """Extract sub-experiment label for B2/B3 methods."""
    if method in ("B2", "B3"):
        after = folder_name.split(complex_name + "_", 1)
        if len(after) > 1:
            return after[1]
    return ""


def discover_predictions(pred_dir: Path, known_complexes: list[str]) -> pd.DataFrame:
    """Discover all prediction files across all methods.

    Returns DataFrame with columns:
        method, complex, sub_experiment, model_idx, model_path,
        confidence_path, plddt_path, pae_path, pde_path
    """
    records = []

    for method_label, mdef in METHOD_DEFS.items():
        method_dir = pred_dir / mdef["folder"]
        if not method_dir.exists():
            log.warning("Method dir not found: %s", method_dir)
            continue

        for result_dir in sorted(method_dir.iterdir()):
            if not result_dir.is_dir():
                continue

            complex_name = extract_complex_name(result_dir.name, known_complexes)
            if complex_name is None:
                log.warning(
                    "Could not extract complex name from %s", result_dir.name
                )
                continue

            sub_exp = extract_sub_experiment(
                result_dir.name, complex_name, method_label
            )

            pred_subdir = result_dir / "predictions"
            if not pred_subdir.exists():
                continue

            inner_dirs = [d for d in pred_subdir.iterdir() if d.is_dir()]
            if not inner_dirs:
                continue

            for inner_dir in inner_dirs:
                model_files = sorted(inner_dir.glob("*_model_*.pdb")) + sorted(
                    inner_dir.glob("*_model_*.cif")
                )
                if not model_files:
                    continue

                for mf in model_files:
                    m = re.search(r"model_(\d+)", mf.name)
                    if not m:
                        continue
                    model_idx = int(m.group(1))

                    stem = mf.stem
                    parent = mf.parent
                    conf_path = parent / f"confidence_{stem}.json"
                    plddt_path = parent / f"plddt_{stem}.npz"
                    pae_path = parent / f"pae_{stem}.npz"
                    pde_path = parent / f"pde_{stem}.npz"

                    records.append(
                        {
                            "method": method_label,
                            "complex": complex_name,
                            "sub_experiment": sub_exp,
                            "model_idx": model_idx,
                            "model_path": str(mf),
                            "confidence_path": str(conf_path)
                            if conf_path.exists()
                            else "",
                            "plddt_path": str(plddt_path)
                            if plddt_path.exists()
                            else "",
                            "pae_path": str(pae_path) if pae_path.exists() else "",
                            "pde_path": str(pde_path) if pde_path.exists() else "",
                        }
                    )

    df = pd.DataFrame(records)
    log.info(
        "Discovered %d prediction files across %d methods",
        len(df),
        df["method"].nunique() if len(df) > 0 else 0,
    )
    for method, group in df.groupby("method"):
        complexes = group["complex"].nunique()
        models = len(group)
        log.info("  %s: %d complexes, %d models", method, complexes, models)
    return df


# ---------------------------------------------------------------------------
# Confidence metrics extraction
# ---------------------------------------------------------------------------


def extract_confidence(confidence_path: str) -> dict:
    """Extract confidence metrics from a Boltz2 confidence JSON."""
    if not confidence_path or not Path(confidence_path).exists():
        return {}
    try:
        with open(confidence_path) as f:
            data = json.load(f)
    except Exception as e:
        log.error("Failed to read confidence JSON %s: %s", confidence_path, e)
        return {}

    result = {}
    for key in [
        "confidence_score",
        "ptm",
        "iptm",
        "ligand_iptm",
        "protein_iptm",
        "complex_plddt",
        "complex_iplddt",
        "complex_pde",
        "complex_ipde",
    ]:
        result[key] = data.get(key, np.nan)

    pair_iptm = data.get("pair_chains_iptm", {})
    result["iptm_AB"] = pair_iptm.get("0", {}).get("1", np.nan)  # antigen-heavy
    result["iptm_AC"] = pair_iptm.get("0", {}).get("2", np.nan)  # antigen-light
    result["iptm_BC"] = pair_iptm.get("1", {}).get("2", np.nan)  # heavy-light

    return result


def extract_all_confidence(predictions_df: pd.DataFrame) -> pd.DataFrame:
    """Extract confidence metrics for all predictions."""
    records = []
    for _, row in predictions_df.iterrows():
        conf = extract_confidence(row["confidence_path"])
        record = {
            "method": row["method"],
            "complex": row["complex"],
            "sub_experiment": row["sub_experiment"],
            "model_idx": row["model_idx"],
        }
        record.update(conf)
        records.append(record)
    return pd.DataFrame(records)


# ---------------------------------------------------------------------------
# DockQ evaluation (Axis B)
# ---------------------------------------------------------------------------


def run_dockq(model_path: str, native_path: str) -> dict:
    """Run DockQ on a single model-native pair using the Python API.

    Returns dict with per-interface DockQ/fnat/iRMSD/LRMSD/F1 and GlobalDockQ.
    """
    try:
        model = load_PDB(model_path)
        native = load_PDB(native_path)
        result_dict, best_dockq_sum = run_on_all_native_interfaces(
            model, native, chain_map=CHAIN_MAP
        )
    except Exception as e:
        log.error("DockQ failed for %s: %s", model_path, e)
        return {"GlobalDockQ": np.nan, "error": str(e)}

    n_interfaces = len(result_dict)
    global_dockq = best_dockq_sum / n_interfaces if n_interfaces > 0 else np.nan

    out = {"GlobalDockQ": global_dockq, "n_interfaces": n_interfaces}

    for iface_key, metrics in result_dict.items():
        for metric_name in ["DockQ", "fnat", "iRMSD", "LRMSD", "fnonnat", "F1"]:
            out[f"{iface_key}_{metric_name}"] = metrics.get(metric_name, np.nan)

    # Extract specific interfaces for antibody-antigen evaluation
    ba = result_dict.get("BA", result_dict.get("AB", {}))
    ca = result_dict.get("CA", result_dict.get("AC", {}))
    bc = result_dict.get("CB", result_dict.get("BC", {}))

    for prefix, iface in [("BA", ba), ("CA", ca), ("BC", bc)]:
        out[f"{prefix}_DockQ"] = iface.get("DockQ", np.nan)
        out[f"{prefix}_fnat"] = iface.get("fnat", np.nan)
        out[f"{prefix}_iRMSD"] = iface.get("iRMSD", np.nan)
        out[f"{prefix}_LRMSD"] = iface.get("LRMSD", np.nan)
        out[f"{prefix}_F1"] = iface.get("F1", np.nan)

    # Compute AB+AC average DockQ (antigen-antibody interfaces)
    ba_dq = out.get("BA_DockQ", np.nan)
    ca_dq = out.get("CA_DockQ", np.nan)
    if not np.isnan(ba_dq) and not np.isnan(ca_dq):
        out["AB_AC_avg_DockQ"] = (ba_dq + ca_dq) / 2
    elif not np.isnan(ba_dq):
        out["AB_AC_avg_DockQ"] = ba_dq
    elif not np.isnan(ca_dq):
        out["AB_AC_avg_DockQ"] = ca_dq
    else:
        out["AB_AC_avg_DockQ"] = np.nan

    return out


def classify_capri(dockq: float) -> str:
    """Classify DockQ score into CAPRI quality category."""
    if np.isnan(dockq):
        return "Error"
    if dockq >= 0.80:
        return "High"
    if dockq >= 0.49:
        return "Medium"
    if dockq >= 0.23:
        return "Acceptable"
    return "Incorrect"


def evaluate_dockq(predictions_df: pd.DataFrame, native_dir: Path) -> pd.DataFrame:
    """Run DockQ on all prediction-native pairs."""
    results = []
    total = len(predictions_df)
    for i, (_, row) in enumerate(predictions_df.iterrows()):
        native_path = native_dir / f"{row['complex']}.pdb"
        if not native_path.exists():
            log.warning("Native not found: %s", native_path)
            continue

        log.info(
            "[%d/%d] DockQ: %s %s model_%d",
            i + 1,
            total,
            row["method"],
            row["complex"],
            row["model_idx"],
        )
        dockq_result = run_dockq(row["model_path"], str(native_path))
        dockq_result["CAPRI"] = classify_capri(
            dockq_result.get("GlobalDockQ", np.nan)
        )
        dockq_result["CAPRI_AB_AC"] = classify_capri(
            dockq_result.get("AB_AC_avg_DockQ", np.nan)
        )

        result = {
            "method": row["method"],
            "complex": row["complex"],
            "sub_experiment": row["sub_experiment"],
            "model_idx": row["model_idx"],
        }
        result.update(dockq_result)
        results.append(result)

    return pd.DataFrame(results)


# ---------------------------------------------------------------------------
# Epitope prediction evaluation (Axis A)
# ---------------------------------------------------------------------------


def extract_epitope_from_structure(
    pdb_path: str, threshold: float = 5.0
) -> set[int]:
    """Extract epitope residue indices from a PDB/CIF structure.

    Epitope = antigen (chain A) residues with any heavy atom within
    `threshold` Angstroms of any antibody (chains B+C) heavy atom.

    Returns set of 1-indexed residue positions within chain A.
    """
    parser = PDBParser(QUIET=True)
    path = Path(pdb_path)

    try:
        if path.suffix == ".cif":
            from Bio.PDB.MMCIFParser import MMCIFParser
            parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure("s", pdb_path)
    except Exception as e:
        log.error("Failed to parse %s: %s", pdb_path, e)
        return set()

    model = structure[0]

    # Collect antibody heavy atoms (chains B + C)
    ab_atoms = []
    for chain_id in ["B", "C"]:
        if chain_id not in model:
            continue
        for residue in model[chain_id].get_residues():
            if residue.id[0] != " ":
                continue
            for atom in residue.get_atoms():
                ab_atoms.append(atom)

    if not ab_atoms:
        log.warning("No antibody atoms found in %s", pdb_path)
        return set()

    ns = NeighborSearch(ab_atoms)

    epitope = set()
    if "A" not in model:
        log.warning("No chain A (antigen) in %s", pdb_path)
        return set()

    for residue in model["A"].get_residues():
        if residue.id[0] != " ":
            continue
        for atom in residue.get_atoms():
            if ns.search(atom.coord, threshold):
                epitope.add(residue.id[1])
                break

    return epitope


def compute_epitope_metrics(
    gt_epitope: set[int],
    pred_epitope: set[int],
    n_antigen_residues: int,
) -> dict:
    """Compute epitope prediction accuracy metrics.

    Args:
        gt_epitope: Ground truth epitope residue indices
        pred_epitope: Predicted epitope residue indices
        n_antigen_residues: Total number of antigen residues (for MCC)

    Returns dict with precision, recall, F1, MCC, jaccard.
    """
    if not gt_epitope and not pred_epitope:
        return {"precision": 1.0, "recall": 1.0, "f1": 1.0, "mcc": 1.0, "jaccard": 1.0}
    if not gt_epitope:
        return {"precision": 0.0, "recall": 0.0, "f1": 0.0, "mcc": 0.0, "jaccard": 0.0}
    if not pred_epitope:
        return {"precision": 0.0, "recall": 0.0, "f1": 0.0, "mcc": 0.0, "jaccard": 0.0}

    tp = len(gt_epitope & pred_epitope)
    fp = len(pred_epitope - gt_epitope)
    fn = len(gt_epitope - pred_epitope)
    tn = n_antigen_residues - tp - fp - fn

    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0.0

    # Matthews Correlation Coefficient
    denom = np.sqrt(float((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)))
    mcc = (tp * tn - fp * fn) / denom if denom > 0 else 0.0

    # Jaccard (IoU)
    union = tp + fp + fn
    jaccard = tp / union if union > 0 else 0.0

    return {"precision": precision, "recall": recall, "f1": f1, "mcc": mcc, "jaccard": jaccard}


def count_antigen_residues(pdb_path: str) -> int:
    """Count the number of standard residues in chain A."""
    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure("s", pdb_path)
    except Exception:
        return 0
    model = structure[0]
    if "A" not in model:
        return 0
    return sum(1 for r in model["A"].get_residues() if r.id[0] == " ")


def evaluate_epitope(
    predictions_df: pd.DataFrame, native_dir: Path, threshold: float = 5.0
) -> pd.DataFrame:
    """Evaluate epitope prediction for all models."""
    gt_cache = {}
    ag_size_cache = {}

    results = []
    total = len(predictions_df)
    for i, (_, row) in enumerate(predictions_df.iterrows()):
        complex_name = row["complex"]
        native_path = native_dir / f"{complex_name}.pdb"

        if not native_path.exists():
            continue

        if complex_name not in gt_cache:
            gt_cache[complex_name] = extract_epitope_from_structure(
                str(native_path), threshold
            )
            ag_size_cache[complex_name] = count_antigen_residues(str(native_path))
            log.info(
                "  GT epitope for %s: %d residues (of %d antigen residues)",
                complex_name,
                len(gt_cache[complex_name]),
                ag_size_cache[complex_name],
            )

        log.info(
            "[%d/%d] Epitope: %s %s model_%d",
            i + 1,
            total,
            row["method"],
            complex_name,
            row["model_idx"],
        )

        pred_epitope = extract_epitope_from_structure(
            row["model_path"], threshold
        )

        metrics = compute_epitope_metrics(
            gt_cache[complex_name],
            pred_epitope,
            ag_size_cache[complex_name],
        )

        result = {
            "method": row["method"],
            "complex": complex_name,
            "sub_experiment": row["sub_experiment"],
            "model_idx": row["model_idx"],
            "gt_epitope_size": len(gt_cache[complex_name]),
            "pred_epitope_size": len(pred_epitope),
            "n_antigen_residues": ag_size_cache[complex_name],
        }
        result.update(metrics)
        results.append(result)

    return pd.DataFrame(results)


# ---------------------------------------------------------------------------
# Antibody-aligned antigen RMSD (Axis C)
# ---------------------------------------------------------------------------


def compute_ab_aligned_ag_rmsd(model_path: str, native_path: str) -> float:
    """Compute antibody-aligned antigen RMSD.

    1. Superimpose predicted antibody (chains B+C) CA atoms onto native antibody CAs.
    2. Apply same transform to predicted antigen (chain A) CAs.
    3. Return RMSD of antigen CAs.
    """
    pdb_parser = PDBParser(QUIET=True)

    try:
        if Path(model_path).suffix == ".cif":
            from Bio.PDB.MMCIFParser import MMCIFParser
            model_parser = MMCIFParser(QUIET=True)
            pred_struct = model_parser.get_structure("pred", model_path)
        else:
            pred_struct = pdb_parser.get_structure("pred", model_path)
        native_struct = pdb_parser.get_structure("native", native_path)
    except Exception as e:
        log.error("Failed to parse structures for RMSD: %s", e)
        return np.nan

    pred_model = pred_struct[0]
    native_model = native_struct[0]

    def get_ca_coords(model, chain_ids):
        coords = []
        for cid in chain_ids:
            if cid not in model:
                continue
            for res in model[cid].get_residues():
                if res.id[0] != " " or "CA" not in res:
                    continue
                coords.append(res["CA"].get_vector().get_array())
        return np.array(coords)

    pred_ab_ca = get_ca_coords(pred_model, ["B", "C"])
    native_ab_ca = get_ca_coords(native_model, ["B", "C"])
    pred_ag_ca = get_ca_coords(pred_model, ["A"])
    native_ag_ca = get_ca_coords(native_model, ["A"])

    if len(pred_ab_ca) == 0 or len(native_ab_ca) == 0:
        return np.nan
    if len(pred_ag_ca) == 0 or len(native_ag_ca) == 0:
        return np.nan

    min_ab = min(len(pred_ab_ca), len(native_ab_ca))
    pred_ab_ca = pred_ab_ca[:min_ab]
    native_ab_ca = native_ab_ca[:min_ab]

    min_ag = min(len(pred_ag_ca), len(native_ag_ca))
    pred_ag_ca = pred_ag_ca[:min_ag]
    native_ag_ca = native_ag_ca[:min_ag]

    # SVD-based alignment on antibody
    pred_ab_mean = pred_ab_ca.mean(axis=0)
    native_ab_mean = native_ab_ca.mean(axis=0)
    pred_ab_centered = pred_ab_ca - pred_ab_mean
    native_ab_centered = native_ab_ca - native_ab_mean

    H = pred_ab_centered.T @ native_ab_centered
    U, S, Vt = np.linalg.svd(H)
    d = np.linalg.det(Vt.T @ U.T)
    sign_matrix = np.diag([1, 1, d])
    R = Vt.T @ sign_matrix @ U.T

    # Apply rotation to antigen
    pred_ag_transformed = (pred_ag_ca - pred_ab_mean) @ R.T + native_ab_mean

    diff = pred_ag_transformed - native_ag_ca
    rmsd = np.sqrt(np.mean(np.sum(diff**2, axis=1)))
    return float(rmsd)


# ---------------------------------------------------------------------------
# Per-residue metrics: CDR3 pLDDT and PAE (Axis D)
# ---------------------------------------------------------------------------


def extract_perres_metrics(
    predictions_df: pd.DataFrame, cdrs: dict
) -> pd.DataFrame:
    """Extract CDR3-specific pLDDT and CDR3-antigen PAE from NPZ files."""
    records = []
    for _, row in predictions_df.iterrows():
        record = {
            "method": row["method"],
            "complex": row["complex"],
            "sub_experiment": row["sub_experiment"],
            "model_idx": row["model_idx"],
        }

        cdr_info = cdrs.get(row["complex"])

        # pLDDT
        if row["plddt_path"] and Path(row["plddt_path"]).exists():
            try:
                plddt_data = np.load(row["plddt_path"])
                plddt = plddt_data["plddt"]
                record["mean_plddt"] = float(np.mean(plddt))

                if cdr_info:
                    n_total = len(plddt)
                    heavy_len = cdr_info["heavy_len"]
                    light_len = cdr_info["light_len"]
                    antigen_len = n_total - heavy_len - light_len

                    if antigen_len > 0:
                        h_start = antigen_len
                        l_start = antigen_len + heavy_len

                        record["plddt_antigen"] = float(
                            np.mean(plddt[:antigen_len])
                        )
                        record["plddt_heavy"] = float(
                            np.mean(plddt[h_start : h_start + heavy_len])
                        )
                        record["plddt_light"] = float(
                            np.mean(plddt[l_start : l_start + light_len])
                        )

                        # CDR3 pLDDT (indices are 0-based within each chain)
                        cdr3_h_global = [h_start + i for i in cdr_info["cdr3_h"]]
                        cdr3_l_global = [l_start + i for i in cdr_info["cdr3_l"]]

                        cdr3_h_global = [
                            i for i in cdr3_h_global if i < n_total
                        ]
                        cdr3_l_global = [
                            i for i in cdr3_l_global if i < n_total
                        ]

                        if cdr3_h_global:
                            record["plddt_cdr3h"] = float(
                                np.mean(plddt[cdr3_h_global])
                            )
                        if cdr3_l_global:
                            record["plddt_cdr3l"] = float(
                                np.mean(plddt[cdr3_l_global])
                            )
                        if cdr3_h_global or cdr3_l_global:
                            all_cdr3 = cdr3_h_global + cdr3_l_global
                            record["plddt_cdr3_all"] = float(
                                np.mean(plddt[all_cdr3])
                            )
            except Exception as e:
                log.error("Failed to load pLDDT %s: %s", row["plddt_path"], e)

        # PAE - CDR3-antigen PAE
        if row["pae_path"] and Path(row["pae_path"]).exists() and cdr_info:
            try:
                pae_data = np.load(row["pae_path"])
                pae = pae_data["pae"]
                n_total = pae.shape[0]
                heavy_len = cdr_info["heavy_len"]
                light_len = cdr_info["light_len"]
                antigen_len = n_total - heavy_len - light_len

                if antigen_len > 0:
                    h_start = antigen_len
                    l_start = antigen_len + heavy_len

                    # Interface PAE: antigen <-> antibody
                    pae_ag_ab = pae[:antigen_len, h_start:]
                    pae_ab_ag = pae[h_start:, :antigen_len]
                    interface_pae = np.concatenate(
                        [pae_ag_ab.flatten(), pae_ab_ag.flatten()]
                    )
                    record["interface_pae"] = float(np.mean(interface_pae))

                    # CDR3-antigen PAE
                    cdr3_h_global = [h_start + i for i in cdr_info["cdr3_h"]]
                    cdr3_l_global = [l_start + i for i in cdr_info["cdr3_l"]]
                    cdr3_all = [
                        i for i in cdr3_h_global + cdr3_l_global if i < n_total
                    ]
                    if cdr3_all:
                        pae_cdr3_ag = pae[
                            np.ix_(cdr3_all, range(antigen_len))
                        ]
                        pae_ag_cdr3 = pae[
                            np.ix_(range(antigen_len), cdr3_all)
                        ]
                        cdr3_ag_pae = np.concatenate(
                            [pae_cdr3_ag.flatten(), pae_ag_cdr3.flatten()]
                        )
                        record["cdr3_antigen_pae"] = float(
                            np.mean(cdr3_ag_pae)
                        )
            except Exception as e:
                log.error("Failed to load PAE %s: %s", row["pae_path"], e)

        records.append(record)

    return pd.DataFrame(records)


# ---------------------------------------------------------------------------
# Model selection & aggregation
# ---------------------------------------------------------------------------


def select_best_models(
    dockq_df: pd.DataFrame,
    confidence_df: pd.DataFrame,
    epitope_df: pd.DataFrame,
) -> pd.DataFrame:
    """Select best models per method-complex using confidence and oracle strategies.

    Returns one row per method-complex-sub_experiment with metrics for both
    selection strategies.
    """
    merge_keys = ["method", "complex", "sub_experiment", "model_idx"]

    for df in [dockq_df, confidence_df, epitope_df]:
        if "sub_experiment" in df.columns:
            df["sub_experiment"] = df["sub_experiment"].fillna("").astype(str)

    merged = pd.merge(dockq_df, confidence_df, on=merge_keys, how="left")
    if len(epitope_df) > 0:
        merged = pd.merge(merged, epitope_df, on=merge_keys, how="left")

    records = []
    for (method, cplx, sub_exp), group in merged.groupby(
        ["method", "complex", "sub_experiment"]
    ):
        record = {
            "method": method,
            "complex": cplx,
            "sub_experiment": sub_exp,
            "n_models": len(group),
        }

        # ------- Oracle: best by GlobalDockQ -------
        best_oracle_idx = group["GlobalDockQ"].idxmax()
        if not pd.isna(best_oracle_idx):
            orow = group.loc[best_oracle_idx]
            record["oracle_model_idx"] = int(orow["model_idx"])
            record["oracle_GlobalDockQ"] = orow["GlobalDockQ"]
            record["oracle_AB_AC_avg_DockQ"] = orow.get("AB_AC_avg_DockQ", np.nan)
            record["oracle_CAPRI"] = orow.get("CAPRI", "")
            record["oracle_BA_fnat"] = orow.get("BA_fnat", np.nan)
            record["oracle_BA_iRMSD"] = orow.get("BA_iRMSD", np.nan)
            if "f1" in orow:
                record["oracle_epitope_f1"] = orow.get("f1", np.nan)
                record["oracle_epitope_precision"] = orow.get("precision", np.nan)
                record["oracle_epitope_recall"] = orow.get("recall", np.nan)
                record["oracle_epitope_mcc"] = orow.get("mcc", np.nan)

        # ------- Confidence: best by confidence_score -------
        conf_vals = group["confidence_score"].dropna()
        if len(conf_vals) > 0:
            best_conf_idx = group.loc[conf_vals.index, "confidence_score"].idxmax()
            crow = group.loc[best_conf_idx]
            record["conf_model_idx"] = int(crow["model_idx"])
            record["conf_GlobalDockQ"] = crow["GlobalDockQ"]
            record["conf_AB_AC_avg_DockQ"] = crow.get("AB_AC_avg_DockQ", np.nan)
            record["conf_CAPRI"] = crow.get("CAPRI", "")
            record["conf_CAPRI_AB_AC"] = crow.get("CAPRI_AB_AC", "")
            record["conf_confidence_score"] = crow["confidence_score"]
            record["conf_iptm"] = crow.get("iptm", np.nan)
            record["conf_complex_plddt"] = crow.get("complex_plddt", np.nan)
            record["conf_complex_iplddt"] = crow.get("complex_iplddt", np.nan)
            record["conf_iptm_AB"] = crow.get("iptm_AB", np.nan)
            record["conf_iptm_AC"] = crow.get("iptm_AC", np.nan)
            record["conf_iptm_BC"] = crow.get("iptm_BC", np.nan)
            record["conf_BA_DockQ"] = crow.get("BA_DockQ", np.nan)
            record["conf_CA_DockQ"] = crow.get("CA_DockQ", np.nan)
            record["conf_BC_DockQ"] = crow.get("BC_DockQ", np.nan)
            record["conf_BA_fnat"] = crow.get("BA_fnat", np.nan)
            record["conf_BA_iRMSD"] = crow.get("BA_iRMSD", np.nan)
            record["conf_BA_LRMSD"] = crow.get("BA_LRMSD", np.nan)
            record["conf_model_path"] = crow.get("model_path", "")

            if "f1" in crow:
                record["conf_epitope_f1"] = crow.get("f1", np.nan)
                record["conf_epitope_precision"] = crow.get("precision", np.nan)
                record["conf_epitope_recall"] = crow.get("recall", np.nan)
                record["conf_epitope_mcc"] = crow.get("mcc", np.nan)
                record["conf_epitope_jaccard"] = crow.get("jaccard", np.nan)
                record["conf_gt_epitope_size"] = crow.get("gt_epitope_size", np.nan)
                record["conf_pred_epitope_size"] = crow.get("pred_epitope_size", np.nan)

        records.append(record)

    return pd.DataFrame(records)


def aggregate_for_multi_config(selected_df: pd.DataFrame) -> pd.DataFrame:
    """For methods with multiple configs (B2, B3), compute best-oracle and mean.

    Single-config methods (B1, Q+) pass through unchanged.
    """
    result_parts = []

    for method, method_group in selected_df.groupby("method"):
        mdef = METHOD_DEFS.get(method, {})

        if not mdef.get("multi_config", False):
            result_parts.append(method_group)
            continue

        best_rows = []
        mean_rows = []
        for cplx, cplx_group in method_group.groupby("complex"):
            best_idx = cplx_group["conf_GlobalDockQ"].idxmax()
            if not pd.isna(best_idx):
                best_row = cplx_group.loc[best_idx].copy()
                best_row["method"] = f"{method}_best"
                best_rows.append(best_row)

            numeric_cols = cplx_group.select_dtypes(include=[np.number]).columns
            mean_vals = cplx_group[numeric_cols].mean()
            mean_row = cplx_group.iloc[0].copy()
            for col in numeric_cols:
                mean_row[col] = mean_vals[col]
            mean_row["method"] = f"{method}_mean"
            mean_row["sub_experiment"] = "mean"
            mean_rows.append(mean_row)

        if best_rows:
            result_parts.append(pd.DataFrame(best_rows))
        if mean_rows:
            result_parts.append(pd.DataFrame(mean_rows))

    if result_parts:
        return pd.concat(result_parts, ignore_index=True)
    return pd.DataFrame()


# ---------------------------------------------------------------------------
# Summary aggregation
# ---------------------------------------------------------------------------


def aggregate_summary(selected_df: pd.DataFrame) -> pd.DataFrame:
    """Compute aggregate summary statistics per method."""
    summary_records = []

    for method, group in selected_df.groupby("method"):
        record = {"method": method, "n_complexes": len(group)}

        # DockQ (confidence-selected)
        for prefix, dockq_col in [
            ("conf_global", "conf_GlobalDockQ"),
            ("conf_ab_ac", "conf_AB_AC_avg_DockQ"),
        ]:
            vals = group[dockq_col].dropna() if dockq_col in group.columns else pd.Series()
            if len(vals) > 0:
                record[f"{prefix}_mean_DockQ"] = float(vals.mean())
                record[f"{prefix}_median_DockQ"] = float(vals.median())
                record[f"{prefix}_std_DockQ"] = float(vals.std())
                record[f"{prefix}_pct_acceptable"] = float(
                    (vals >= 0.23).mean() * 100
                )
                record[f"{prefix}_pct_medium"] = float(
                    (vals >= 0.49).mean() * 100
                )
                record[f"{prefix}_pct_high"] = float(
                    (vals >= 0.80).mean() * 100
                )

        # Oracle DockQ
        oracle_dockq = group["oracle_GlobalDockQ"].dropna() if "oracle_GlobalDockQ" in group.columns else pd.Series()
        if len(oracle_dockq) > 0:
            record["oracle_mean_DockQ"] = float(oracle_dockq.mean())
            record["oracle_median_DockQ"] = float(oracle_dockq.median())

        # Epitope metrics (confidence-selected)
        for metric in ["epitope_f1", "epitope_precision", "epitope_recall", "epitope_mcc", "epitope_jaccard"]:
            col = f"conf_{metric}"
            if col in group.columns:
                vals = group[col].dropna()
                if len(vals) > 0:
                    record[f"mean_{metric}"] = float(vals.mean())
                    record[f"median_{metric}"] = float(vals.median())
                    record[f"std_{metric}"] = float(vals.std())

        # Confidence metrics
        for metric in ["conf_confidence_score", "conf_iptm", "conf_complex_plddt", "conf_complex_iplddt"]:
            if metric in group.columns:
                vals = group[metric].dropna()
                if len(vals) > 0:
                    record[f"mean_{metric}"] = float(vals.mean())

        summary_records.append(record)

    return pd.DataFrame(summary_records)


# ---------------------------------------------------------------------------
# Statistical tests
# ---------------------------------------------------------------------------


def pairwise_wilcoxon(
    selected_df: pd.DataFrame,
    reference: str = "B1",
    metrics: list[str] | None = None,
) -> pd.DataFrame:
    """Run paired Wilcoxon signed-rank tests comparing each method to reference."""
    if metrics is None:
        metrics = [
            "conf_GlobalDockQ",
            "conf_AB_AC_avg_DockQ",
            "conf_epitope_f1",
            "conf_epitope_mcc",
        ]

    results = []
    ref_data = selected_df[selected_df["method"] == reference]

    for method in selected_df["method"].unique():
        if method == reference:
            continue
        method_data = selected_df[selected_df["method"] == method]

        for metric in metrics:
            if metric not in ref_data.columns or metric not in method_data.columns:
                continue

            ref_vals = ref_data[["complex", metric]].rename(
                columns={metric: "ref_val"}
            )
            method_vals = method_data[["complex", metric]].rename(
                columns={metric: "method_val"}
            )
            merged = pd.merge(ref_vals, method_vals, on="complex", how="inner")
            merged = merged.dropna(subset=["ref_val", "method_val"])

            if len(merged) < 3:
                results.append(
                    {
                        "method": method,
                        "vs": reference,
                        "metric": metric,
                        "n_paired": len(merged),
                        "statistic": np.nan,
                        "p_value": np.nan,
                        "median_delta": np.nan,
                        "mean_delta": np.nan,
                    }
                )
                continue

            delta = merged["method_val"].values - merged["ref_val"].values

            try:
                stat, p_val = stats.wilcoxon(delta, alternative="two-sided")
            except Exception:
                stat, p_val = np.nan, np.nan

            results.append(
                {
                    "method": method,
                    "vs": reference,
                    "metric": metric,
                    "n_paired": len(merged),
                    "statistic": stat,
                    "p_value": p_val,
                    "median_delta": float(np.median(delta)),
                    "mean_delta": float(np.mean(delta)),
                    "wins": int((delta > 0).sum()),
                    "losses": int((delta < 0).sum()),
                    "ties": int((delta == 0).sum()),
                }
            )

    return pd.DataFrame(results)


# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------


def _make_scatter(ax, df, x_method, y_method, metric, x_label, y_label, title, color):
    """Helper: paired scatter plot of metric for two methods."""
    x_data = df[df["method"] == x_method][["complex", metric]].rename(
        columns={metric: "x_val"}
    )
    y_data = df[df["method"] == y_method][["complex", metric]].rename(
        columns={metric: "y_val"}
    )
    scatter_df = pd.merge(x_data, y_data, on="complex").dropna()

    if len(scatter_df) == 0:
        return False

    ax.scatter(
        scatter_df["x_val"],
        scatter_df["y_val"],
        s=50,
        c=color,
        edgecolors="black",
        linewidths=0.5,
        zorder=3,
    )
    lim = [0, 1.05]
    ax.plot(lim, lim, "k--", alpha=0.5, label="y=x (no change)")
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_title(title)
    ax.legend(fontsize=8)
    ax.set_xlim(lim)
    ax.set_ylim(lim)

    for _, r in scatter_df.iterrows():
        ax.annotate(
            r["complex"],
            (r["x_val"], r["y_val"]),
            fontsize=5,
            alpha=0.6,
        )
    return True


def _make_delta_waterfall(ax, df, ref_method, test_method, metric, ylabel, title):
    """Helper: per-complex delta waterfall bar chart."""
    ref_data = df[df["method"] == ref_method][["complex", metric]].rename(
        columns={metric: "ref_val"}
    )
    test_data = df[df["method"] == test_method][["complex", metric]].rename(
        columns={metric: "test_val"}
    )
    delta_df = pd.merge(ref_data, test_data, on="complex").dropna()
    delta_df["delta"] = delta_df["test_val"] - delta_df["ref_val"]
    delta_df = delta_df.sort_values("delta")

    if len(delta_df) == 0:
        return False

    bar_colors = ["#2ca02c" if d >= 0 else "#d62728" for d in delta_df["delta"]]
    ax.bar(range(len(delta_df)), delta_df["delta"], color=bar_colors, edgecolor="black", linewidth=0.3)
    ax.set_xticks(range(len(delta_df)))
    ax.set_xticklabels(delta_df["complex"], rotation=90, fontsize=6)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.axhline(0, color="black", linewidth=0.5)

    wins = (delta_df["delta"] > 0).sum()
    losses = (delta_df["delta"] < 0).sum()
    ties = (delta_df["delta"] == 0).sum()
    ax.text(
        0.02, 0.98,
        f"Wins: {wins}  Losses: {losses}  Ties: {ties}",
        transform=ax.transAxes,
        fontsize=8,
        verticalalignment="top",
        bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5),
    )
    return True


def generate_plots(selected_df: pd.DataFrame, output_dir: Path):
    """Generate all evaluation plots."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        log.warning("matplotlib not available, skipping plots")
        return

    fig_dir = output_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    methods_present = sorted(selected_df["method"].unique())
    method_order = [
        m for m in [
            "B1", "B2_best", "B2_mean", "B3_best", "B3_mean", "Q_v1", "Q+"
        ] if m in methods_present
    ]
    if not method_order:
        method_order = methods_present

    colors = {
        "B1": "#1f77b4",
        "B2_best": "#ff7f0e",
        "B2_mean": "#ffbb78",
        "B3_best": "#2ca02c",
        "B3_mean": "#98df8a",
        "Q_v1": "#9467bd",
        "Q+": "#d62728",
    }

    def _boxplot(metric_col, ylabel, title, filename):
        if metric_col not in selected_df.columns:
            return
        data_to_plot = []
        labels = []
        for m in method_order:
            vals = selected_df[selected_df["method"] == m][metric_col].dropna().values
            if len(vals) > 0:
                data_to_plot.append(vals)
                labels.append(m)
        if not data_to_plot:
            return

        fig, ax = plt.subplots(figsize=(10, 6))
        bp = ax.boxplot(data_to_plot, tick_labels=labels, patch_artist=True)
        for patch, method in zip(bp["boxes"], labels):
            patch.set_facecolor(colors.get(method, "gray"))
            patch.set_alpha(0.7)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        if "DockQ" in metric_col:
            ax.axhline(0.23, color="orange", linestyle="--", alpha=0.5, label="Acceptable")
            ax.axhline(0.49, color="green", linestyle="--", alpha=0.5, label="Medium")
            ax.axhline(0.80, color="blue", linestyle="--", alpha=0.5, label="High")
            ax.legend(loc="upper right", fontsize=8)
        if "f1" in metric_col or "mcc" in metric_col or "jaccard" in metric_col:
            ax.set_ylim(-0.05, 1.05)

        # Add median annotation
        for i, (d, label) in enumerate(zip(data_to_plot, labels)):
            median_val = np.median(d)
            ax.text(i + 1, median_val, f" {median_val:.3f}", fontsize=7, va="bottom", ha="left")

        plt.tight_layout()
        plt.savefig(fig_dir / filename, dpi=150)
        plt.close()
        log.info("Saved %s", filename)

    # ------------------------------------------------------------------
    # 1. DockQ box plot (AB+AC interface)
    # ------------------------------------------------------------------
    _boxplot(
        "conf_AB_AC_avg_DockQ",
        "DockQ (AB+AC avg, confidence-selected)",
        "Docking Quality: Antigen-Antibody Interfaces",
        "dockq_boxplot_ab_ac.png",
    )

    # ------------------------------------------------------------------
    # 2. Epitope F1 box plot
    # ------------------------------------------------------------------
    _boxplot(
        "conf_epitope_f1",
        "Epitope F1 (confidence-selected)",
        "Epitope Prediction Accuracy (F1)",
        "epitope_f1_boxplot.png",
    )

    # ------------------------------------------------------------------
    # 3. Epitope MCC box plot
    # ------------------------------------------------------------------
    _boxplot(
        "conf_epitope_mcc",
        "Epitope MCC (confidence-selected)",
        "Epitope Prediction Accuracy (MCC)",
        "epitope_mcc_boxplot.png",
    )

    # ------------------------------------------------------------------
    # 4. DockQ scatter: Q+ vs B1
    # ------------------------------------------------------------------
    if "Q+" in methods_present and "B1" in methods_present:
        fig, ax = plt.subplots(figsize=(7, 7))
        ok = _make_scatter(
            ax, selected_df, "B1", "Q+",
            "conf_AB_AC_avg_DockQ",
            "B1 (vanilla Boltz2) DockQ",
            "Q+ (Epitope Refinement v2) DockQ",
            "Q+ vs B1: DockQ (AB+AC avg)",
            colors["Q+"],
        )
        if ok:
            plt.tight_layout()
            plt.savefig(fig_dir / "dockq_scatter_qplus_vs_b1.png", dpi=150)
            log.info("Saved dockq_scatter_qplus_vs_b1.png")
        plt.close()

    # ------------------------------------------------------------------
    # 5. Epitope F1 scatter: Q+ vs B1
    # ------------------------------------------------------------------
    if "Q+" in methods_present and "B1" in methods_present:
        fig, ax = plt.subplots(figsize=(7, 7))
        ok = _make_scatter(
            ax, selected_df, "B1", "Q+",
            "conf_epitope_f1",
            "B1 (vanilla) Epitope F1",
            "Q+ Epitope F1",
            "Q+ vs B1: Epitope Prediction F1",
            colors["Q+"],
        )
        if ok:
            plt.tight_layout()
            plt.savefig(fig_dir / "epitope_f1_scatter_qplus_vs_b1.png", dpi=150)
            log.info("Saved epitope_f1_scatter_qplus_vs_b1.png")
        plt.close()

    # ------------------------------------------------------------------
    # 6. DockQ scatter: Q+ vs Q v1
    # ------------------------------------------------------------------
    if "Q+" in methods_present and "Q_v1" in methods_present:
        fig, ax = plt.subplots(figsize=(7, 7))
        ok = _make_scatter(
            ax, selected_df, "Q_v1", "Q+",
            "conf_AB_AC_avg_DockQ",
            "Q v1 DockQ",
            "Q+ v2 DockQ",
            "Q+ v2 vs Q v1: DockQ (AB+AC avg)",
            colors["Q+"],
        )
        if ok:
            plt.tight_layout()
            plt.savefig(fig_dir / "dockq_scatter_qplus_vs_qv1.png", dpi=150)
            log.info("Saved dockq_scatter_qplus_vs_qv1.png")
        plt.close()

    # ------------------------------------------------------------------
    # 7. Epitope F1 scatter: Q+ vs Q v1
    # ------------------------------------------------------------------
    if "Q+" in methods_present and "Q_v1" in methods_present:
        fig, ax = plt.subplots(figsize=(7, 7))
        ok = _make_scatter(
            ax, selected_df, "Q_v1", "Q+",
            "conf_epitope_f1",
            "Q v1 Epitope F1",
            "Q+ v2 Epitope F1",
            "Q+ v2 vs Q v1: Epitope Prediction F1",
            colors["Q+"],
        )
        if ok:
            plt.tight_layout()
            plt.savefig(fig_dir / "epitope_f1_scatter_qplus_vs_qv1.png", dpi=150)
            log.info("Saved epitope_f1_scatter_qplus_vs_qv1.png")
        plt.close()

    # ------------------------------------------------------------------
    # 8. CAPRI stacked bar
    # ------------------------------------------------------------------
    capri_order = ["High", "Medium", "Acceptable", "Incorrect"]
    capri_colors = {
        "High": "#2ca02c",
        "Medium": "#98df8a",
        "Acceptable": "#ffbb78",
        "Incorrect": "#d62728",
    }

    if "conf_CAPRI_AB_AC" in selected_df.columns:
        capri_data = selected_df.groupby("method")["conf_CAPRI_AB_AC"].value_counts().unstack(fill_value=0)
        for cat in capri_order:
            if cat not in capri_data.columns:
                capri_data[cat] = 0
        capri_data = capri_data[capri_order]
        capri_pct = capri_data.div(capri_data.sum(axis=1), axis=0) * 100
        capri_pct = capri_pct.reindex(
            [m for m in method_order if m in capri_pct.index]
        )

        if len(capri_pct) > 0:
            fig, ax = plt.subplots(figsize=(10, 5))
            bottom = np.zeros(len(capri_pct))
            for cat in capri_order:
                vals = capri_pct[cat].values
                ax.bar(
                    capri_pct.index,
                    vals,
                    bottom=bottom,
                    label=cat,
                    color=capri_colors[cat],
                    edgecolor="white",
                    linewidth=0.5,
                )
                bottom += vals
            ax.set_ylabel("Percentage (%)")
            ax.set_title("CAPRI Quality Distribution (AB+AC avg, confidence-selected)")
            ax.legend(loc="upper right")
            ax.set_ylim(0, 100)
            plt.tight_layout()
            plt.savefig(fig_dir / "capri_stacked_bar.png", dpi=150)
            plt.close()
            log.info("Saved capri_stacked_bar.png")

    # ------------------------------------------------------------------
    # 9. Delta waterfall: DockQ Q+ - B1
    # ------------------------------------------------------------------
    if "Q+" in methods_present and "B1" in methods_present:
        fig, ax = plt.subplots(figsize=(max(10, 48 * 0.4), 6))
        ok = _make_delta_waterfall(
            ax, selected_df, "B1", "Q+",
            "conf_AB_AC_avg_DockQ",
            "DockQ Delta (Q+ - B1)",
            "Per-Complex DockQ Improvement: Q+ vs B1",
        )
        if ok:
            plt.tight_layout()
            plt.savefig(fig_dir / "delta_waterfall_dockq.png", dpi=150)
            log.info("Saved delta_waterfall_dockq.png")
        plt.close()

    # ------------------------------------------------------------------
    # 10. Delta waterfall: Epitope F1 Q+ - B1
    # ------------------------------------------------------------------
    if "Q+" in methods_present and "B1" in methods_present:
        fig, ax = plt.subplots(figsize=(max(10, 48 * 0.4), 6))
        ok = _make_delta_waterfall(
            ax, selected_df, "B1", "Q+",
            "conf_epitope_f1",
            "Epitope F1 Delta (Q+ - B1)",
            "Per-Complex Epitope F1 Improvement: Q+ vs B1",
        )
        if ok:
            plt.tight_layout()
            plt.savefig(fig_dir / "delta_waterfall_epitope_f1.png", dpi=150)
            log.info("Saved delta_waterfall_epitope_f1.png")
        plt.close()

    # ------------------------------------------------------------------
    # 11. Confidence vs DockQ scatter
    # ------------------------------------------------------------------
    if "conf_confidence_score" in selected_df.columns and "conf_GlobalDockQ" in selected_df.columns:
        fig, ax = plt.subplots(figsize=(8, 6))
        for method in method_order:
            mdata = selected_df[selected_df["method"] == method]
            valid = mdata.dropna(subset=["conf_confidence_score", "conf_GlobalDockQ"])
            if len(valid) > 0:
                ax.scatter(
                    valid["conf_confidence_score"],
                    valid["conf_GlobalDockQ"],
                    label=method,
                    color=colors.get(method, "gray"),
                    alpha=0.7,
                    s=40,
                    edgecolors="black",
                    linewidths=0.3,
                )
        ax.set_xlabel("Confidence Score")
        ax.set_ylabel("DockQ (Global)")
        ax.set_title("Confidence Score vs DockQ (best-by-confidence models)")
        ax.legend()
        plt.tight_layout()
        plt.savefig(fig_dir / "confidence_vs_dockq.png", dpi=150)
        plt.close()
        log.info("Saved confidence_vs_dockq.png")


# ---------------------------------------------------------------------------
# Q v1 results loading
# ---------------------------------------------------------------------------


def load_qv1_results(qv1_csv_path: Path) -> pd.DataFrame:
    """Load Q v1 selected_models_aggregated.csv and relabel method as Q_v1."""
    df = pd.read_csv(qv1_csv_path)
    # Keep only the Q rows (not B1/B2/B3 which we already have)
    df = df[df["method"] == "Q"].copy()
    df["method"] = "Q_v1"
    return df


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(
        description="Evaluate Strategy Q+ v2 (Iterative Epitope Refinement v2) vs baselines"
    )
    parser.add_argument(
        "--pred_dir",
        type=Path,
        required=True,
        help="Root predictions directory containing method subfolders",
    )
    parser.add_argument(
        "--native_dir",
        type=Path,
        required=True,
        help="Directory with ground truth PDB files (pdb_minimized/)",
    )
    parser.add_argument(
        "--cdr_csv",
        type=Path,
        required=True,
        help="Path to cdrs.csv with CDR definitions",
    )
    parser.add_argument(
        "--out_dir",
        type=Path,
        default=Path("evaluation/results_q_plus"),
        help="Output directory for results",
    )
    parser.add_argument(
        "--qv1_selected",
        type=Path,
        default=None,
        help="Path to Q v1 selected_models_aggregated.csv for v1 vs v2 comparison",
    )
    parser.add_argument(
        "--epitope_threshold",
        type=float,
        default=5.0,
        help="Distance threshold for epitope extraction (Angstroms)",
    )
    parser.add_argument(
        "--skip_dockq",
        action="store_true",
        help="Skip DockQ/epitope computation (use cached results)",
    )
    parser.add_argument(
        "--skip_plots",
        action="store_true",
        help="Skip plot generation",
    )
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    # ---- Step 0: Load CDR definitions ----
    log.info("Loading CDR definitions from %s", args.cdr_csv)
    cdrs = load_cdrs(args.cdr_csv)
    log.info("Loaded CDR info for %d complexes", len(cdrs))

    # ---- Step 1: Discover predictions ----
    log.info("Discovering prediction files in %s", args.pred_dir)
    known_complexes = get_known_complexes(args.native_dir)
    log.info("Found %d ground truth complexes", len(known_complexes))

    predictions_df = discover_predictions(args.pred_dir, known_complexes)
    if len(predictions_df) == 0:
        log.error("No predictions found! Check --pred_dir path.")
        sys.exit(1)
    predictions_df.to_csv(args.out_dir / "predictions_discovered.csv", index=False)

    # ---- Step 2: Extract confidence metrics ----
    log.info("Extracting confidence metrics...")
    confidence_df = extract_all_confidence(predictions_df)
    confidence_df.to_csv(args.out_dir / "confidence_all.csv", index=False)

    # ---- Step 3: Run DockQ (Axis B) ----
    dockq_cache = args.out_dir / "dockq_all.csv"
    if args.skip_dockq and dockq_cache.exists():
        log.info("Loading cached DockQ results from %s", dockq_cache)
        dockq_df = pd.read_csv(dockq_cache)
    else:
        log.info("Running DockQ on all predictions...")
        dockq_df = evaluate_dockq(predictions_df, args.native_dir)
        dockq_df.to_csv(dockq_cache, index=False)
    log.info("DockQ computed for %d models", len(dockq_df))

    # ---- Step 4: Evaluate epitope prediction (Axis A) ----
    epitope_cache = args.out_dir / "epitope_all.csv"
    if args.skip_dockq and epitope_cache.exists():
        log.info("Loading cached epitope results from %s", epitope_cache)
        epitope_df = pd.read_csv(epitope_cache)
    else:
        log.info("Evaluating epitope prediction (threshold=%.1f A)...", args.epitope_threshold)
        epitope_df = evaluate_epitope(predictions_df, args.native_dir, args.epitope_threshold)
        epitope_df.to_csv(epitope_cache, index=False)
    log.info("Epitope evaluated for %d models", len(epitope_df))

    # ---- Step 5: Extract per-residue metrics (Axis D) ----
    log.info("Extracting per-residue metrics (CDR3 pLDDT, PAE)...")
    perres_df = extract_perres_metrics(predictions_df, cdrs)
    perres_df.to_csv(args.out_dir / "perres_all.csv", index=False)

    # ---- Step 6: Model selection ----
    log.info("Selecting best models (confidence + oracle)...")
    selected_df = select_best_models(dockq_df, confidence_df, epitope_df)
    selected_df.to_csv(args.out_dir / "selected_models.csv", index=False)

    # ---- Step 7: Handle multi-config methods (B2, B3) ----
    log.info("Aggregating multi-config methods...")
    selected_agg = aggregate_for_multi_config(selected_df)
    selected_agg.to_csv(args.out_dir / "selected_models_aggregated.csv", index=False)

    # ---- Step 8: Compute antibody-aligned antigen RMSD (Axis C) ----
    log.info("Computing antibody-aligned antigen RMSD for confidence-selected models...")
    rmsd_records = []
    for _, row in selected_agg.iterrows():
        model_path = row.get("conf_model_path", "")
        if not model_path or pd.isna(model_path) or not Path(model_path).exists():
            continue
        native_path = args.native_dir / f"{row['complex']}.pdb"
        if not native_path.exists():
            continue
        rmsd = compute_ab_aligned_ag_rmsd(model_path, str(native_path))
        rmsd_records.append({
            "method": row["method"],
            "complex": row["complex"],
            "ab_aligned_ag_rmsd": rmsd,
        })
    if rmsd_records:
        rmsd_df = pd.DataFrame(rmsd_records)
        rmsd_df.to_csv(args.out_dir / "ab_aligned_ag_rmsd.csv", index=False)
        selected_agg = pd.merge(
            selected_agg,
            rmsd_df,
            on=["method", "complex"],
            how="left",
        )

    # ---- Step 9: Load Q v1 results (if provided) ----
    if args.qv1_selected and args.qv1_selected.exists():
        log.info("Loading Q v1 results from %s", args.qv1_selected)
        qv1_df = load_qv1_results(args.qv1_selected)
        log.info("Loaded Q v1 results for %d complexes", len(qv1_df))
        selected_agg = pd.concat([selected_agg, qv1_df], ignore_index=True)
        # Ensure sub_experiment is string after concat
        if "sub_experiment" in selected_agg.columns:
            selected_agg["sub_experiment"] = selected_agg["sub_experiment"].fillna("").astype(str)

    # Save final aggregated data (including Q v1 if present)
    selected_agg.to_csv(args.out_dir / "selected_models_aggregated.csv", index=False)

    # ---- Step 10: Summary tables ----
    log.info("Computing summary statistics...")
    summary_df = aggregate_summary(selected_agg)
    summary_df.to_csv(args.out_dir / "summary.csv", index=False)

    # Print summary
    print("\n" + "=" * 80)
    print("EVALUATION SUMMARY — Strategy Q+ v2 (Iterative Epitope Refinement v2)")
    print("=" * 80)

    print("\n--- Table 1: Docking Quality (confidence-selected, AB+AC avg) ---")
    docking_cols = [
        "method", "n_complexes", "conf_ab_ac_mean_DockQ", "conf_ab_ac_median_DockQ",
        "conf_ab_ac_std_DockQ", "conf_ab_ac_pct_acceptable", "conf_ab_ac_pct_medium",
        "conf_ab_ac_pct_high",
    ]
    available_cols = [c for c in docking_cols if c in summary_df.columns]
    if available_cols:
        print(summary_df[available_cols].to_string(index=False, float_format="%.4f"))

    print("\n--- Table 2: Epitope Prediction (confidence-selected) ---")
    epitope_cols = [
        "method", "mean_epitope_precision", "mean_epitope_recall",
        "mean_epitope_f1", "median_epitope_f1", "mean_epitope_mcc",
        "mean_epitope_jaccard",
    ]
    available_cols = [c for c in epitope_cols if c in summary_df.columns]
    if available_cols:
        print(summary_df[available_cols].to_string(index=False, float_format="%.4f"))

    print("\n--- Table 3: Confidence Metrics ---")
    conf_cols = [
        "method", "mean_conf_confidence_score", "mean_conf_iptm",
        "mean_conf_complex_plddt", "mean_conf_complex_iplddt",
    ]
    available_cols = [c for c in conf_cols if c in summary_df.columns]
    if available_cols:
        print(summary_df[available_cols].to_string(index=False, float_format="%.4f"))

    print("\n--- Table 4: Oracle DockQ (upper bound) ---")
    oracle_cols = ["method", "oracle_mean_DockQ", "oracle_median_DockQ"]
    available_cols = [c for c in oracle_cols if c in summary_df.columns]
    if available_cols:
        print(summary_df[available_cols].to_string(index=False, float_format="%.4f"))

    # Per-complex breakdown
    per_complex_cols = [
        "method", "complex", "conf_AB_AC_avg_DockQ", "conf_epitope_f1",
        "conf_epitope_mcc", "conf_confidence_score",
    ]
    available_cols = [c for c in per_complex_cols if c in selected_agg.columns]
    if available_cols:
        pivot_data = selected_agg[available_cols].copy()
        pivot_data.to_csv(args.out_dir / "per_complex_breakdown.csv", index=False)

    # ---- Step 11: Statistical tests ----
    log.info("Running statistical tests (Wilcoxon)...")
    wilcoxon_results = []

    # Q+ vs B1
    w_vs_b1 = pairwise_wilcoxon(selected_agg, reference="B1")
    if len(w_vs_b1) > 0:
        wilcoxon_results.append(w_vs_b1)

    # Q+ vs Q_v1 (if available)
    if "Q_v1" in selected_agg["method"].values:
        w_vs_qv1 = pairwise_wilcoxon(selected_agg, reference="Q_v1")
        # Keep only the Q+ row
        w_vs_qv1 = w_vs_qv1[w_vs_qv1["method"] == "Q+"]
        if len(w_vs_qv1) > 0:
            wilcoxon_results.append(w_vs_qv1)

    if wilcoxon_results:
        wilcoxon_df = pd.concat(wilcoxon_results, ignore_index=True)
    else:
        wilcoxon_df = pd.DataFrame()
    wilcoxon_df.to_csv(args.out_dir / "wilcoxon_tests.csv", index=False)

    print("\n--- Statistical Tests (Wilcoxon) ---")
    if len(wilcoxon_df) > 0:
        print(wilcoxon_df.to_string(index=False, float_format="%.4f"))

    # ---- Step 12: Hard-zero analysis (Q+ helps where B1 fails) ----
    if "Q+" in selected_agg["method"].values and "B1" in selected_agg["method"].values:
        b1_f1 = selected_agg[selected_agg["method"] == "B1"][
            ["complex", "conf_epitope_f1"]
        ].rename(columns={"conf_epitope_f1": "b1_f1"})
        qp_f1 = selected_agg[selected_agg["method"] == "Q+"][
            ["complex", "conf_epitope_f1"]
        ].rename(columns={"conf_epitope_f1": "qp_f1"})
        hard_zero_df = pd.merge(b1_f1, qp_f1, on="complex").dropna()

        b1_zeros = hard_zero_df[hard_zero_df["b1_f1"] == 0]
        rescued = b1_zeros[b1_zeros["qp_f1"] > 0]

        print(f"\n--- Hard-Zero Analysis ---")
        print(f"Complexes where B1 has F1=0: {len(b1_zeros)}")
        print(f"Of those, Q+ rescued (F1>0):  {len(rescued)}")
        if len(rescued) > 0:
            print(f"Rescued complexes: {', '.join(rescued['complex'].tolist())}")
            print(f"Mean F1 of rescued: {rescued['qp_f1'].mean():.4f}")

    # ---- Step 13: Plots ----
    if not args.skip_plots:
        log.info("Generating plots...")
        generate_plots(selected_agg, args.out_dir)

    print(f"\nAll results saved to {args.out_dir}/")
    print("  - predictions_discovered.csv")
    print("  - confidence_all.csv")
    print("  - dockq_all.csv")
    print("  - epitope_all.csv")
    print("  - perres_all.csv")
    print("  - selected_models.csv")
    print("  - selected_models_aggregated.csv")
    print("  - ab_aligned_ag_rmsd.csv")
    print("  - summary.csv")
    print("  - wilcoxon_tests.csv")
    print("  - per_complex_breakdown.csv")
    if not args.skip_plots:
        print("  - figures/*.png")


if __name__ == "__main__":
    main()
