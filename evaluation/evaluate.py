#!/usr/bin/env python
"""
Evaluation script for Strategy G+ vs baselines.

Computes DockQ, confidence metrics, per-residue metrics, ensemble diversity,
model selection, and generates summary tables + plots.

Usage:
    conda activate dockq2
    python evaluation/evaluate.py --pred_dir /path/to/predictions_examples --native_dir /path/to/pdb_minimized

See EVALUATION_STRATEGY.md for full details.
"""

import argparse
import ast
import csv
import json
import logging
import re
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.PDB import MMCIFParser, PDBParser, Superimposer
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

# Method definitions: (label, folder_name, naming_pattern)
# naming_pattern describes how the boltz_results folder maps to complex names
METHOD_DEFS = {
    "B1": {
        "folder": "antigen_cut",
        "prefix": "boltz_results_",
        "suffix": "",
    },
    "B2": {
        "folder": "antigen_cut_vhvl_msa_pocket_ab_boltz_post_2023_omm",
        "prefix": "boltz_results_restraint_to_A_",
        "suffix": "",  # has extra _CHAIN_RES_NUM suffix
    },
    "B3": {
        "folder": "antigen_cut_contact_restraints",
        "prefix": "boltz_results_restraint_",
        "suffix": "",  # has extra _TYPE_NUM suffix
    },
    "Gplus": {
        "folder": "new_feature_cdr3_beta",
        "prefix": "boltz_results_",
        "suffix": "_cdr3_beta",
    },
}


# ---------------------------------------------------------------------------
# Discovery: find prediction files and map to ground truths
# ---------------------------------------------------------------------------


def get_known_complexes(native_dir: Path) -> list[str]:
    """Get list of complex names from ground truth PDB files."""
    complexes = []
    for pdb in sorted(native_dir.glob("*.pdb")):
        complexes.append(pdb.stem)
    return complexes


def extract_complex_name(folder_name: str, known_complexes: list[str]) -> str | None:
    """Extract complex name from a boltz_results folder name by matching
    against the known complex list."""
    for cname in sorted(known_complexes, key=len, reverse=True):
        if cname in folder_name:
            return cname
    return None


def extract_sub_experiment(folder_name: str, complex_name: str, method: str) -> str:
    """Extract sub-experiment label for B2/B3 methods."""
    if method == "B2":
        # boltz_results_restraint_to_A_7TRH_HBG_B_W_109 -> B_W_109
        after = folder_name.split(complex_name + "_", 1)
        if len(after) > 1:
            return after[1]
    elif method == "B3":
        # boltz_results_restraint_7TRH_HBG_hbond_23 -> hbond_23
        after = folder_name.split(complex_name + "_", 1)
        if len(after) > 1:
            return after[1]
    return ""


def discover_predictions(pred_dir: Path, known_complexes: list[str]) -> pd.DataFrame:
    """Discover all prediction files across all methods.

    Returns DataFrame with columns:
        method, complex, sub_experiment, model_idx, model_path, confidence_path,
        plddt_path, pae_path, pde_path
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

            # Find the predictions subfolder - it's the only dir inside predictions/
            pred_subdir = result_dir / "predictions"
            if not pred_subdir.exists():
                continue

            inner_dirs = [d for d in pred_subdir.iterdir() if d.is_dir()]
            if not inner_dirs:
                continue

            for inner_dir in inner_dirs:
                # Find model files (PDB or CIF)
                model_files = sorted(inner_dir.glob("*_model_*.pdb")) + sorted(
                    inner_dir.glob("*_model_*.cif")
                )
                if not model_files:
                    continue

                for mf in model_files:
                    # Extract model index from filename like 7TRH_HBG_model_0.pdb
                    m = re.search(r"model_(\d+)", mf.name)
                    if not m:
                        continue
                    model_idx = int(m.group(1))

                    # Derive confidence/plddt/pae/pde paths
                    stem = mf.stem  # e.g. 7TRH_HBG_model_0
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
        df["method"].nunique(),
    )
    for method, group in df.groupby("method"):
        complexes = group["complex"].nunique()
        models = len(group)
        log.info("  %s: %d complexes, %d models", method, complexes, models)
    return df


# ---------------------------------------------------------------------------
# DockQ evaluation
# ---------------------------------------------------------------------------


def run_dockq(model_path: str, native_path: str) -> dict:
    """Run DockQ on a single model-native pair using the Python API.

    Returns dict with: GlobalDockQ, per-interface DockQ/fnat/iRMSD/LRMSD,
    and the best interface metrics (BA = antigen-heavy).
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

    # Extract per-interface metrics
    for iface_key, metrics in result_dict.items():
        for metric_name in ["DockQ", "fnat", "iRMSD", "LRMSD", "fnonnat", "F1"]:
            out[f"{iface_key}_{metric_name}"] = metrics.get(metric_name, np.nan)

    # Extract the antibody-antigen interface specifically (BA and CA)
    # BA = heavy-antigen, CA = light-antigen
    ba = result_dict.get("BA", result_dict.get("AB", {}))
    ca = result_dict.get("CA", result_dict.get("AC", {}))

    out["BA_DockQ"] = ba.get("DockQ", np.nan)
    out["BA_fnat"] = ba.get("fnat", np.nan)
    out["BA_iRMSD"] = ba.get("iRMSD", np.nan)
    out["BA_LRMSD"] = ba.get("LRMSD", np.nan)
    out["CA_DockQ"] = ca.get("DockQ", np.nan)
    out["CA_fnat"] = ca.get("fnat", np.nan)
    out["CA_iRMSD"] = ca.get("iRMSD", np.nan)
    out["CA_LRMSD"] = ca.get("LRMSD", np.nan)

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
        dockq_result["CAPRI"] = classify_capri(dockq_result.get("GlobalDockQ", np.nan))

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

    # Pairwise chain iPTM
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
# Per-residue metrics (pLDDT, PAE)
# ---------------------------------------------------------------------------


def load_cdrs(cdrs_csv: Path) -> dict:
    """Load CDR residue indices from cdrs.csv.

    Returns dict: complex_name -> {
        'cdr3_h': [list of 0-indexed residue positions in heavy chain],
        'cdr3_l': [list of 0-indexed residue positions in light chain],
        'heavy_len': int,
        'light_len': int,
    }
    """
    cdrs = {}
    with open(cdrs_csv) as f:
        reader = csv.DictReader(f)
        for row in reader:
            cname = row["complex"]
            heavy_seq = row["heavy"]
            light_seq = row["light"]

            # Parse CDR3 indices (they are 0-indexed positions within each chain)
            cdr3_h = ast.literal_eval(row["cdr3_h"])
            cdr3_l = ast.literal_eval(row["cdr3_l"])

            cdrs[cname] = {
                "cdr3_h": cdr3_h,
                "cdr3_l": cdr3_l,
                "heavy_len": len(heavy_seq),
                "light_len": len(light_seq),
            }
    return cdrs


def extract_perres_metrics(
    predictions_df: pd.DataFrame, cdrs: dict
) -> pd.DataFrame:
    """Extract per-residue metrics from NPZ files.

    Computes CDR3-specific pLDDT and interface PAE.
    """
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
                    # Chain A (antigen) comes first, then B (heavy), then C (light)
                    # We need to figure out antigen length:
                    # total = antigen_len + heavy_len + light_len
                    heavy_len = cdr_info["heavy_len"]
                    light_len = cdr_info["light_len"]
                    antigen_len = n_total - heavy_len - light_len

                    if antigen_len > 0:
                        # pLDDT indices for each chain
                        ag_start, ag_end = 0, antigen_len
                        h_start, h_end = antigen_len, antigen_len + heavy_len
                        l_start, l_end = (
                            antigen_len + heavy_len,
                            n_total,
                        )

                        record["plddt_antigen"] = float(
                            np.mean(plddt[ag_start:ag_end])
                        )
                        record["plddt_heavy"] = float(np.mean(plddt[h_start:h_end]))
                        record["plddt_light"] = float(np.mean(plddt[l_start:l_end]))

                        # CDR3 pLDDT (indices are 0-based within each chain)
                        cdr3_h_global = [h_start + i for i in cdr_info["cdr3_h"]]
                        cdr3_l_global = [l_start + i for i in cdr_info["cdr3_l"]]

                        # Filter valid indices
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

        # PAE - compute interface PAE (antibody-antigen block)
        if row["pae_path"] and Path(row["pae_path"]).exists() and cdr_info:
            try:
                pae_data = np.load(row["pae_path"])
                pae = pae_data["pae"]
                n_total = pae.shape[0]
                heavy_len = cdr_info["heavy_len"]
                light_len = cdr_info["light_len"]
                antigen_len = n_total - heavy_len - light_len

                if antigen_len > 0:
                    ag_start, ag_end = 0, antigen_len
                    h_start, h_end = antigen_len, antigen_len + heavy_len
                    l_start, l_end = antigen_len + heavy_len, n_total

                    # Interface PAE: antigen rows, antibody cols (and vice versa)
                    pae_ag_ab = pae[ag_start:ag_end, h_start:l_end]
                    pae_ab_ag = pae[h_start:l_end, ag_start:ag_end]
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
                        pae_cdr3_ag = pae[np.ix_(cdr3_all, range(ag_start, ag_end))]
                        pae_ag_cdr3 = pae[np.ix_(range(ag_start, ag_end), cdr3_all)]
                        cdr3_ag_pae = np.concatenate(
                            [pae_cdr3_ag.flatten(), pae_ag_cdr3.flatten()]
                        )
                        record["cdr3_antigen_pae"] = float(np.mean(cdr3_ag_pae))
            except Exception as e:
                log.error("Failed to load PAE %s: %s", row["pae_path"], e)

        records.append(record)

    return pd.DataFrame(records)


# ---------------------------------------------------------------------------
# Ensemble diversity
# ---------------------------------------------------------------------------


def compute_pairwise_ca_rmsd(
    model_paths: list[str], chains: str = "ABC"
) -> tuple[float, float, int]:
    """Compute pairwise CA-RMSD between all model files.

    Returns (mean_rmsd, std_rmsd, n_pairs).
    """
    parsers = {"pdb": PDBParser(QUIET=True), "cif": MMCIFParser(QUIET=True)}

    def load_ca_atoms(path: str) -> list:
        p = Path(path)
        ext = p.suffix.lstrip(".")
        parser = parsers.get(ext, parsers["pdb"])
        structure = parser.get_structure("s", str(p))
        model = structure[0]
        atoms = []
        for chain_id in chains:
            if chain_id in model:
                for residue in model[chain_id]:
                    if "CA" in residue:
                        atoms.append(residue["CA"])
        return atoms

    structures = []
    for mp in model_paths:
        try:
            atoms = load_ca_atoms(mp)
            structures.append(atoms)
        except Exception as e:
            log.error("Failed to load %s for RMSD: %s", mp, e)

    if len(structures) < 2:
        return np.nan, np.nan, 0

    rmsds = []
    for (i, atoms1), (j, atoms2) in combinations(enumerate(structures), 2):
        if len(atoms1) != len(atoms2):
            log.warning(
                "Atom count mismatch: model %d has %d, model %d has %d",
                i,
                len(atoms1),
                j,
                len(atoms2),
            )
            continue
        try:
            sup = Superimposer()
            sup.set_atoms(atoms1, atoms2)
            rmsds.append(sup.rms)
        except Exception as e:
            log.error("Superimposer failed for pair %d-%d: %s", i, j, e)

    if not rmsds:
        return np.nan, np.nan, 0
    return float(np.mean(rmsds)), float(np.std(rmsds)), len(rmsds)


def compute_ensemble_diversity(predictions_df: pd.DataFrame) -> pd.DataFrame:
    """Compute ensemble diversity for each method-complex combination."""
    records = []

    for (method, cplx, sub_exp), group in predictions_df.groupby(
        ["method", "complex", "sub_experiment"]
    ):
        if len(group) < 2:
            records.append(
                {
                    "method": method,
                    "complex": cplx,
                    "sub_experiment": sub_exp,
                    "n_models": len(group),
                    "mean_ca_rmsd": np.nan,
                    "std_ca_rmsd": np.nan,
                    "n_pairs": 0,
                }
            )
            continue

        model_paths = group["model_path"].tolist()
        mean_rmsd, std_rmsd, n_pairs = compute_pairwise_ca_rmsd(model_paths)

        records.append(
            {
                "method": method,
                "complex": cplx,
                "sub_experiment": sub_exp,
                "n_models": len(group),
                "mean_ca_rmsd": mean_rmsd,
                "std_ca_rmsd": std_rmsd,
                "n_pairs": n_pairs,
            }
        )

    return pd.DataFrame(records)


# ---------------------------------------------------------------------------
# Model selection & aggregation
# ---------------------------------------------------------------------------


def select_best_models(
    dockq_df: pd.DataFrame, confidence_df: pd.DataFrame
) -> pd.DataFrame:
    """Apply model selection strategies: confidence-based and oracle.

    Returns DataFrame with one row per method-complex-sub_experiment,
    containing the selected model index and metrics for both strategies.
    """
    # Merge DockQ and confidence data
    merge_keys = ["method", "complex", "sub_experiment", "model_idx"]
    merged = pd.merge(dockq_df, confidence_df, on=merge_keys, how="left")

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

        # All-models stats
        dockq_vals = group["GlobalDockQ"].dropna()
        if len(dockq_vals) > 0:
            record["mean_DockQ"] = float(dockq_vals.mean())
            record["std_DockQ"] = float(dockq_vals.std())
            record["n_acceptable"] = int((dockq_vals >= 0.23).sum())

        # Oracle selection: best DockQ
        best_oracle_idx = group["GlobalDockQ"].idxmax()
        if not pd.isna(best_oracle_idx):
            oracle_row = group.loc[best_oracle_idx]
            record["oracle_model_idx"] = int(oracle_row["model_idx"])
            record["oracle_DockQ"] = oracle_row["GlobalDockQ"]
            record["oracle_CAPRI"] = oracle_row.get("CAPRI", "")
            for col in ["BA_fnat", "BA_iRMSD", "BA_LRMSD", "CA_fnat", "CA_iRMSD", "CA_LRMSD"]:
                if col in oracle_row:
                    record[f"oracle_{col}"] = oracle_row[col]

        # Confidence-based selection: best confidence_score
        conf_vals = group["confidence_score"].dropna()
        if len(conf_vals) > 0:
            best_conf_idx = group.loc[conf_vals.index, "confidence_score"].idxmax()
            conf_row = group.loc[best_conf_idx]
            record["conf_model_idx"] = int(conf_row["model_idx"])
            record["conf_DockQ"] = conf_row["GlobalDockQ"]
            record["conf_CAPRI"] = conf_row.get("CAPRI", "")
            record["conf_confidence_score"] = conf_row["confidence_score"]
            record["conf_iptm"] = conf_row.get("iptm", np.nan)
            record["conf_complex_iplddt"] = conf_row.get("complex_iplddt", np.nan)
            for col in ["BA_fnat", "BA_iRMSD", "BA_LRMSD", "CA_fnat", "CA_iRMSD", "CA_LRMSD"]:
                if col in conf_row:
                    record[f"conf_{col}"] = conf_row[col]

        # iPTM-based selection
        iptm_vals = group["iptm"].dropna()
        if len(iptm_vals) > 0:
            best_iptm_idx = group.loc[iptm_vals.index, "iptm"].idxmax()
            iptm_row = group.loc[best_iptm_idx]
            record["iptm_model_idx"] = int(iptm_row["model_idx"])
            record["iptm_DockQ"] = iptm_row["GlobalDockQ"]

        records.append(record)

    return pd.DataFrame(records)


def aggregate_summary(selected_df: pd.DataFrame) -> pd.DataFrame:
    """Compute aggregate summary statistics per method."""
    summary_records = []

    for method, group in selected_df.groupby("method"):
        record = {"method": method, "n_complexes": len(group)}

        # Best-of-N (confidence selection)
        conf_dockq = group["conf_DockQ"].dropna()
        if len(conf_dockq) > 0:
            record["conf_mean_DockQ"] = float(conf_dockq.mean())
            record["conf_median_DockQ"] = float(conf_dockq.median())
            record["conf_pct_acceptable"] = float((conf_dockq >= 0.23).mean() * 100)
            record["conf_pct_medium"] = float((conf_dockq >= 0.49).mean() * 100)
            record["conf_pct_high"] = float((conf_dockq >= 0.80).mean() * 100)

        # Oracle selection
        oracle_dockq = group["oracle_DockQ"].dropna()
        if len(oracle_dockq) > 0:
            record["oracle_mean_DockQ"] = float(oracle_dockq.mean())
            record["oracle_median_DockQ"] = float(oracle_dockq.median())
            record["oracle_pct_acceptable"] = float(
                (oracle_dockq >= 0.23).mean() * 100
            )
            record["oracle_pct_medium"] = float((oracle_dockq >= 0.49).mean() * 100)
            record["oracle_pct_high"] = float((oracle_dockq >= 0.80).mean() * 100)

        # Mean across all models
        mean_dockq = group["mean_DockQ"].dropna()
        if len(mean_dockq) > 0:
            record["all_models_mean_DockQ"] = float(mean_dockq.mean())

        summary_records.append(record)

    return pd.DataFrame(summary_records)


# ---------------------------------------------------------------------------
# Statistical tests
# ---------------------------------------------------------------------------


def pairwise_wilcoxon(
    selected_df: pd.DataFrame, reference: str = "B1", metric: str = "conf_DockQ"
) -> pd.DataFrame:
    """Run paired Wilcoxon signed-rank tests comparing each method to reference."""
    ref_data = selected_df[selected_df["method"] == reference][
        ["complex", metric]
    ].rename(columns={metric: f"{reference}_{metric}"})

    results = []
    for method in selected_df["method"].unique():
        if method == reference:
            continue
        method_data = selected_df[selected_df["method"] == method][
            ["complex", metric]
        ].rename(columns={metric: f"{method}_{metric}"})

        merged = pd.merge(ref_data, method_data, on="complex", how="inner")
        if len(merged) < 3:
            results.append(
                {
                    "method": method,
                    "vs": reference,
                    "n_paired": len(merged),
                    "statistic": np.nan,
                    "p_value": np.nan,
                    "median_delta": np.nan,
                }
            )
            continue

        ref_vals = merged[f"{reference}_{metric}"].values
        method_vals = merged[f"{method}_{metric}"].values
        delta = method_vals - ref_vals

        try:
            stat, p_val = stats.wilcoxon(delta, alternative="two-sided")
        except Exception:
            stat, p_val = np.nan, np.nan

        results.append(
            {
                "method": method,
                "vs": reference,
                "n_paired": len(merged),
                "statistic": stat,
                "p_value": p_val,
                "median_delta": float(np.median(delta)),
                "mean_delta": float(np.mean(delta)),
            }
        )

    return pd.DataFrame(results)


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------


def generate_plots(
    selected_df: pd.DataFrame,
    dockq_df: pd.DataFrame,
    confidence_df: pd.DataFrame,
    diversity_df: pd.DataFrame,
    output_dir: Path,
):
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
    colors = {"B1": "#1f77b4", "B2": "#ff7f0e", "B3": "#2ca02c", "Gplus": "#d62728"}

    # ------------------------------------------------------------------
    # 1. DockQ bar chart per method (confidence-selected best-of-N)
    # ------------------------------------------------------------------
    summary = selected_df.groupby("method")["conf_DockQ"].agg(["mean", "std", "count"])
    summary = summary.reindex([m for m in ["B1", "B2", "B3", "Gplus"] if m in summary.index])
    if len(summary) > 0:
        fig, ax = plt.subplots(figsize=(8, 5))
        bars = ax.bar(
            summary.index,
            summary["mean"],
            yerr=summary["std"],
            capsize=5,
            color=[colors.get(m, "gray") for m in summary.index],
            edgecolor="black",
            linewidth=0.5,
        )
        ax.set_ylabel("DockQ (confidence-selected best-of-N)")
        ax.set_title("DockQ by Method")
        ax.set_ylim(0, 1)
        for bar, val in zip(bars, summary["mean"]):
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height() + 0.02,
                f"{val:.3f}",
                ha="center",
                va="bottom",
                fontsize=10,
            )
        plt.tight_layout()
        plt.savefig(fig_dir / "dockq_bar.png", dpi=150)
        plt.close()
        log.info("Saved dockq_bar.png")

    # ------------------------------------------------------------------
    # 2. DockQ scatter: G+ vs B1
    # ------------------------------------------------------------------
    if "Gplus" in methods_present and "B1" in methods_present:
        gplus = selected_df[selected_df["method"] == "Gplus"][
            ["complex", "conf_DockQ"]
        ].rename(columns={"conf_DockQ": "Gplus_DockQ"})
        b1 = selected_df[selected_df["method"] == "B1"][
            ["complex", "conf_DockQ"]
        ].rename(columns={"conf_DockQ": "B1_DockQ"})
        scatter_df = pd.merge(b1, gplus, on="complex")

        if len(scatter_df) > 0:
            fig, ax = plt.subplots(figsize=(6, 6))
            ax.scatter(
                scatter_df["B1_DockQ"],
                scatter_df["Gplus_DockQ"],
                s=50,
                c=colors["Gplus"],
                edgecolors="black",
                linewidths=0.5,
                zorder=3,
            )
            lim = [0, max(1, scatter_df[["B1_DockQ", "Gplus_DockQ"]].max().max() + 0.05)]
            ax.plot(lim, lim, "k--", alpha=0.5, label="y=x")
            ax.set_xlabel("B1 (no steering) DockQ")
            ax.set_ylabel("G+ DockQ")
            ax.set_title("G+ vs Baseline 1")
            ax.legend()
            ax.set_xlim(lim)
            ax.set_ylim(lim)

            for _, r in scatter_df.iterrows():
                ax.annotate(
                    r["complex"],
                    (r["B1_DockQ"], r["Gplus_DockQ"]),
                    fontsize=6,
                    alpha=0.7,
                )
            plt.tight_layout()
            plt.savefig(fig_dir / "dockq_scatter_gplus_vs_b1.png", dpi=150)
            plt.close()
            log.info("Saved dockq_scatter_gplus_vs_b1.png")

    # ------------------------------------------------------------------
    # 3. CAPRI classification stacked bar
    # ------------------------------------------------------------------
    capri_order = ["High", "Medium", "Acceptable", "Incorrect"]
    capri_colors = {
        "High": "#2ca02c",
        "Medium": "#98df8a",
        "Acceptable": "#ffbb78",
        "Incorrect": "#d62728",
    }

    capri_data = selected_df.groupby("method")["conf_CAPRI"].value_counts().unstack(fill_value=0)
    # Reorder columns
    for cat in capri_order:
        if cat not in capri_data.columns:
            capri_data[cat] = 0
    capri_data = capri_data[capri_order]
    capri_pct = capri_data.div(capri_data.sum(axis=1), axis=0) * 100
    capri_pct = capri_pct.reindex([m for m in ["B1", "B2", "B3", "Gplus"] if m in capri_pct.index])

    if len(capri_pct) > 0:
        fig, ax = plt.subplots(figsize=(8, 5))
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
        ax.set_title("CAPRI Quality Distribution (confidence-selected)")
        ax.legend(loc="upper right")
        ax.set_ylim(0, 100)
        plt.tight_layout()
        plt.savefig(fig_dir / "capri_stacked_bar.png", dpi=150)
        plt.close()
        log.info("Saved capri_stacked_bar.png")

    # ------------------------------------------------------------------
    # 4. Confidence vs DockQ correlation
    # ------------------------------------------------------------------
    merged_conf = pd.merge(
        dockq_df[["method", "complex", "sub_experiment", "model_idx", "GlobalDockQ"]],
        confidence_df[
            ["method", "complex", "sub_experiment", "model_idx", "confidence_score"]
        ],
        on=["method", "complex", "sub_experiment", "model_idx"],
        how="inner",
    )
    valid = merged_conf.dropna(subset=["GlobalDockQ", "confidence_score"])
    if len(valid) > 0:
        fig, ax = plt.subplots(figsize=(7, 5))
        for method in methods_present:
            mdata = valid[valid["method"] == method]
            ax.scatter(
                mdata["confidence_score"],
                mdata["GlobalDockQ"],
                label=method,
                color=colors.get(method, "gray"),
                alpha=0.7,
                s=30,
                edgecolors="black",
                linewidths=0.3,
            )
        ax.set_xlabel("Confidence Score")
        ax.set_ylabel("DockQ")
        ax.set_title("Confidence Score vs DockQ (all models)")
        ax.legend()
        plt.tight_layout()
        plt.savefig(fig_dir / "confidence_vs_dockq.png", dpi=150)
        plt.close()
        log.info("Saved confidence_vs_dockq.png")

    # ------------------------------------------------------------------
    # 5. Ensemble diversity box plot
    # ------------------------------------------------------------------
    div_valid = diversity_df.dropna(subset=["mean_ca_rmsd"])
    if len(div_valid) > 0:
        fig, ax = plt.subplots(figsize=(8, 5))
        methods_in_div = [m for m in ["B1", "B2", "B3", "Gplus"] if m in div_valid["method"].values]
        data_to_plot = [
            div_valid[div_valid["method"] == m]["mean_ca_rmsd"].values
            for m in methods_in_div
        ]
        bp = ax.boxplot(data_to_plot, tick_labels=methods_in_div, patch_artist=True)
        for patch, method in zip(bp["boxes"], methods_in_div):
            patch.set_facecolor(colors.get(method, "gray"))
            patch.set_alpha(0.7)
        ax.set_ylabel("Mean Pairwise CA-RMSD (A)")
        ax.set_title("Ensemble Diversity")
        plt.tight_layout()
        plt.savefig(fig_dir / "ensemble_diversity_boxplot.png", dpi=150)
        plt.close()
        log.info("Saved ensemble_diversity_boxplot.png")

    # ------------------------------------------------------------------
    # 6. Per-complex DockQ comparison (sorted bar chart)
    # ------------------------------------------------------------------
    if len(methods_present) >= 2:
        pivot = selected_df.pivot_table(
            index="complex", columns="method", values="conf_DockQ"
        )
        # Sort by first available baseline
        sort_col = "B1" if "B1" in pivot.columns else pivot.columns[0]
        pivot = pivot.sort_values(sort_col, ascending=True)

        if len(pivot) > 0:
            fig, ax = plt.subplots(figsize=(max(10, len(pivot) * 0.5), 6))
            x = np.arange(len(pivot))
            width = 0.8 / len(methods_present)
            for i, method in enumerate(
                [m for m in ["B1", "B2", "B3", "Gplus"] if m in pivot.columns]
            ):
                offset = (i - len(methods_present) / 2 + 0.5) * width
                vals = pivot[method].values
                ax.bar(
                    x + offset,
                    vals,
                    width,
                    label=method,
                    color=colors.get(method, "gray"),
                    edgecolor="black",
                    linewidth=0.3,
                )
            ax.set_xticks(x)
            ax.set_xticklabels(pivot.index, rotation=45, ha="right", fontsize=7)
            ax.set_ylabel("DockQ")
            ax.set_title("Per-Complex DockQ (confidence-selected)")
            ax.legend()
            ax.set_ylim(0, 1)
            plt.tight_layout()
            plt.savefig(fig_dir / "dockq_per_complex.png", dpi=150)
            plt.close()
            log.info("Saved dockq_per_complex.png")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def parse_args():
    parser = argparse.ArgumentParser(
        description="Evaluate Strategy G+ predictions vs baselines using DockQ"
    )
    parser.add_argument(
        "--pred_dir",
        type=Path,
        required=True,
        help="Path to predictions_examples/ directory containing method subfolders",
    )
    parser.add_argument(
        "--native_dir",
        type=Path,
        required=True,
        help="Path to pdb_minimized/ directory with ground truth PDB files",
    )
    parser.add_argument(
        "--cdrs_csv",
        type=Path,
        default=None,
        help="Path to cdrs.csv (default: auto-detect from repo)",
    )
    parser.add_argument(
        "--output_dir",
        type=Path,
        default=Path("evaluation/results"),
        help="Output directory for results (default: evaluation/results)",
    )
    parser.add_argument(
        "--methods",
        nargs="+",
        default=None,
        help="Methods to evaluate (default: all found). Choices: B1 B2 B3 Gplus",
    )
    parser.add_argument(
        "--skip_dockq",
        action="store_true",
        help="Skip DockQ computation (load from existing CSV)",
    )
    parser.add_argument(
        "--skip_plots",
        action="store_true",
        help="Skip plot generation",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    # Find CDR CSV
    cdrs_csv = args.cdrs_csv
    if cdrs_csv is None:
        # Try common locations relative to script
        candidates = [
            Path(__file__).resolve().parent.parent / "examples" / "cdrs.csv",
            Path("examples/cdrs.csv"),
        ]
        for c in candidates:
            if c.exists():
                cdrs_csv = c
                break
    if cdrs_csv and cdrs_csv.exists():
        log.info("Loading CDRs from %s", cdrs_csv)
        cdrs = load_cdrs(cdrs_csv)
        log.info("Loaded CDR info for %d complexes", len(cdrs))
    else:
        log.warning("CDR CSV not found, per-residue CDR metrics will be skipped")
        cdrs = {}

    # ---- Step 1: Discover predictions ----
    log.info("=" * 60)
    log.info("Step 1: Discovering predictions")
    log.info("=" * 60)
    known_complexes = get_known_complexes(args.native_dir)
    log.info("Found %d ground truth complexes", len(known_complexes))

    predictions_df = discover_predictions(args.pred_dir, known_complexes)
    if predictions_df.empty:
        log.error("No predictions found! Check --pred_dir path.")
        sys.exit(1)

    # Filter methods if specified
    if args.methods:
        predictions_df = predictions_df[predictions_df["method"].isin(args.methods)]
        log.info("Filtered to methods: %s", args.methods)

    predictions_df.to_csv(output_dir / "discovered_predictions.csv", index=False)

    # ---- Step 2: DockQ evaluation ----
    log.info("=" * 60)
    log.info("Step 2: DockQ evaluation")
    log.info("=" * 60)
    dockq_csv = output_dir / "dockq_results.csv"
    if args.skip_dockq and dockq_csv.exists():
        log.info("Loading existing DockQ results from %s", dockq_csv)
        dockq_df = pd.read_csv(dockq_csv)
    else:
        dockq_df = evaluate_dockq(predictions_df, args.native_dir)
        dockq_df.to_csv(dockq_csv, index=False)
        log.info("Saved DockQ results to %s", dockq_csv)

    # ---- Step 3: Confidence metrics ----
    log.info("=" * 60)
    log.info("Step 3: Extracting confidence metrics")
    log.info("=" * 60)
    confidence_df = extract_all_confidence(predictions_df)
    confidence_df.to_csv(output_dir / "confidence_metrics.csv", index=False)
    log.info("Saved confidence metrics (%d rows)", len(confidence_df))

    # ---- Step 4: Per-residue metrics ----
    log.info("=" * 60)
    log.info("Step 4: Extracting per-residue metrics")
    log.info("=" * 60)
    perres_df = extract_perres_metrics(predictions_df, cdrs)
    perres_df.to_csv(output_dir / "perres_metrics.csv", index=False)
    log.info("Saved per-residue metrics (%d rows)", len(perres_df))

    # ---- Step 5: Ensemble diversity ----
    log.info("=" * 60)
    log.info("Step 5: Computing ensemble diversity")
    log.info("=" * 60)
    diversity_df = compute_ensemble_diversity(predictions_df)
    diversity_df.to_csv(output_dir / "diversity_metrics.csv", index=False)
    log.info("Saved diversity metrics (%d rows)", len(diversity_df))

    # ---- Step 6: Model selection & aggregation ----
    log.info("=" * 60)
    log.info("Step 6: Model selection & aggregation")
    log.info("=" * 60)
    selected_df = select_best_models(dockq_df, confidence_df)
    selected_df.to_csv(output_dir / "selected_models.csv", index=False)
    log.info("Saved selected models (%d rows)", len(selected_df))

    summary_df = aggregate_summary(selected_df)
    summary_df.to_csv(output_dir / "summary_table.csv", index=False)

    # Print summary
    log.info("")
    log.info("=" * 60)
    log.info("SUMMARY")
    log.info("=" * 60)
    print("\n" + summary_df.to_string(index=False))
    print()

    # ---- Step 7: Statistical tests ----
    log.info("=" * 60)
    log.info("Step 7: Statistical tests")
    log.info("=" * 60)
    if "B1" in selected_df["method"].values:
        wilcoxon_df = pairwise_wilcoxon(selected_df, reference="B1")
        wilcoxon_df.to_csv(output_dir / "wilcoxon_tests.csv", index=False)
        print("\nWilcoxon signed-rank tests (vs B1):")
        print(wilcoxon_df.to_string(index=False))
        print()
    else:
        log.info("Skipping Wilcoxon tests (B1 not present)")

    # ---- Step 8: Plots ----
    if not args.skip_plots:
        log.info("=" * 60)
        log.info("Step 8: Generating plots")
        log.info("=" * 60)
        generate_plots(selected_df, dockq_df, confidence_df, diversity_df, output_dir)

    # ---- Per-complex comparison table ----
    log.info("=" * 60)
    log.info("Per-complex comparison")
    log.info("=" * 60)
    comparison = selected_df.pivot_table(
        index="complex",
        columns="method",
        values=["conf_DockQ", "oracle_DockQ", "conf_CAPRI"],
        aggfunc="first",
    )
    comparison.to_csv(output_dir / "per_complex_comparison.csv")
    print("\nPer-complex DockQ (confidence-selected):")
    dockq_pivot = selected_df.pivot_table(
        index="complex", columns="method", values="conf_DockQ", aggfunc="first"
    )
    print(dockq_pivot.to_string())
    print()

    log.info("=" * 60)
    log.info("Evaluation complete! Results saved to %s", output_dir)
    log.info("=" * 60)


if __name__ == "__main__":
    main()
