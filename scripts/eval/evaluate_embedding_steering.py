"""
Evaluation script for Strategy V: Embedding-Space CDR3 Steering.

Compares Strategy V against baselines using:
  - Tier 1: DockQ (interface quality, CAPRI standard)
  - Tier 2: CDR3-specific RMSD (after framework alignment)
  - Tier 3: Boltz2 confidence metrics (ipTM, pLDDT, etc.)
  - Tier 4: Epitope contact analysis (precision/recall/F1)

Usage:
    conda activate boltz
    python scripts/eval/evaluate_embedding_steering.py \
        --gt_dir pdb_minimized \
        --pred_dir predictions_examples \
        --cdrs_csv examples/cdrs.csv \
        --out_dir eval_results_embedding_steering

DockQ is called via `conda run -n dockq2 DockQ ...` subprocess.
Plots require matplotlib (available in dockq2 env, not boltz). Use --skip_plots or
run from dockq2 env for plots. All other metrics work in the boltz env.
"""

import argparse
import ast
import csv
import json
import os
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

import numpy as np
from Bio.Data.IUPACData import protein_letters_3to1
from Bio.PDB import PDBParser, Superimposer
from Bio.PDB.MMCIFParser import MMCIFParser
from scipy import stats


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

CHAINS_ANTIGEN = "A"
CHAINS_HEAVY = "B"
CHAINS_LIGHT = "C"

CAPRI_THRESHOLDS = {
    "incorrect": (0.0, 0.23),
    "acceptable": (0.23, 0.49),
    "medium": (0.49, 0.80),
    "high": (0.80, 1.01),
}

DOCKQ_CONDA_ENV = "dockq2"

# Methods and their folder patterns
METHODS = {
    "B1": {
        "folder": "antigen_cut",
        "pattern": "boltz_results_{complex}",
        "pred_subdir": "predictions/{complex}",
        "file_pattern": "{complex}_model_{i}.pdb",
        "multi_setting": False,
    },
    "V": {
        "folder": "new_feature_cdr3_beta",
        "pattern": "boltz_results_{complex}_cdr3_beta",
        "pred_subdir": "predictions/{complex}_cdr3_beta",
        "file_pattern": "{complex}_cdr3_beta_model_{i}",  # .pdb or .cif
        "multi_setting": False,
    },
    "B2": {
        "folder": "antigen_cut_contact_restraints",
        "pattern": "boltz_results_restraint_{complex}_*",
        "multi_setting": True,
    },
    "B3": {
        "folder": "antigen_cut_vhvl_msa_pocket_ab_boltz_post_2023_omm",
        "pattern": "boltz_results_restraint_to_A_{complex}_*",
        "multi_setting": True,
    },
}


# ---------------------------------------------------------------------------
# Helpers: Structure I/O
# ---------------------------------------------------------------------------

def parse_structure(filepath):
    """Parse a PDB or CIF file and return a Bio.PDB Structure."""
    filepath = str(filepath)
    if filepath.endswith(".cif"):
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)
    return parser.get_structure("s", filepath)


def residue_to_aa(residue):
    """Convert a Bio.PDB Residue to single-letter amino acid code."""
    try:
        return protein_letters_3to1[residue.resname.capitalize()]
    except KeyError:
        return "X"


def get_chain_sequence(chain):
    """Extract amino acid sequence and list of residues from a chain."""
    residues = [r for r in chain.get_residues() if r.id[0] == " "]
    seq = "".join(residue_to_aa(r) for r in residues)
    return seq, residues


def get_ca_atom(residue):
    """Get CA atom from a residue, or None if missing."""
    if "CA" in residue:
        return residue["CA"]
    return None


# ---------------------------------------------------------------------------
# Helpers: Sequence alignment for residue mapping
# ---------------------------------------------------------------------------

def find_subsequence_offset(short_seq, long_seq):
    """Find where short_seq starts in long_seq. Returns offset or -1."""
    idx = long_seq.find(short_seq)
    if idx >= 0:
        return idx
    # Try with minor mismatches (first/last residue differences from clipping)
    for start in range(min(5, len(short_seq))):
        sub = short_seq[start:]
        idx = long_seq.find(sub)
        if idx >= 0:
            return idx - start
    return -1


def align_chains_by_sequence(pred_chain, gt_chain):
    """Align prediction chain residues to ground truth chain residues by sequence.

    Returns list of (pred_residue, gt_residue) pairs for matched positions.
    Handles different numbering, different chain lengths (GT may include
    constant regions), and sequence variants (insertions/deletions/mutations).

    Uses pairwise sequence alignment to handle non-identical sequences.
    """
    pred_seq, pred_residues = get_chain_sequence(pred_chain)
    gt_seq, gt_residues = get_chain_sequence(gt_chain)

    # Try simple substring match first (fast path for identical sequences)
    if len(pred_seq) <= len(gt_seq):
        offset = find_subsequence_offset(pred_seq, gt_seq)
        if offset >= 0:
            pairs = []
            for i, pred_res in enumerate(pred_residues):
                gt_idx = offset + i
                if 0 <= gt_idx < len(gt_residues):
                    pairs.append((pred_res, gt_residues[gt_idx]))
            return pairs
    else:
        offset = find_subsequence_offset(gt_seq, pred_seq)
        if offset >= 0:
            pairs = []
            for i, gt_res in enumerate(gt_residues):
                pred_idx = offset + i
                if 0 <= pred_idx < len(pred_residues):
                    pairs.append((pred_residues[pred_idx], gt_res))
            return pairs

    # Fallback: pairwise alignment for non-identical sequences
    from Bio.Align import PairwiseAligner
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -0.5
    aligner.extend_gap_score = -0.1
    alignment = aligner.align(pred_seq, gt_seq)[0]
    aligned_pred, aligned_gt = alignment[0], alignment[1]

    pairs = []
    pred_idx = 0
    gt_idx = 0
    for p_char, g_char in zip(aligned_pred, aligned_gt):
        if p_char != "-" and g_char != "-":
            pairs.append((pred_residues[pred_idx], gt_residues[gt_idx]))
        if p_char != "-":
            pred_idx += 1
        if g_char != "-":
            gt_idx += 1

    return pairs


# ---------------------------------------------------------------------------
# CDR3 RMSD computation
# ---------------------------------------------------------------------------

def parse_cdr_indices(cdr_str):
    """Parse CDR index string like '[26, 27, 28]' to list of 1-indexed ints."""
    return ast.literal_eval(cdr_str)


def load_cdrs(cdrs_csv):
    """Load CDR definitions from csv. Returns dict keyed by complex name."""
    cdrs = {}
    with open(cdrs_csv) as f:
        reader = csv.DictReader(f)
        for row in reader:
            name = row["complex"]
            cdrs[name] = {
                "cdr1_h": parse_cdr_indices(row["cdr1_h"]),
                "cdr2_h": parse_cdr_indices(row["cdr2_h"]),
                "cdr3_h": parse_cdr_indices(row["cdr3_h"]),
                "cdr1_l": parse_cdr_indices(row["cdr1_l"]),
                "cdr2_l": parse_cdr_indices(row["cdr2_l"]),
                "cdr3_l": parse_cdr_indices(row["cdr3_l"]),
            }
    return cdrs


def compute_cdr3_rmsd(pred_path, gt_path, cdr_info):
    """Compute CDR3 RMSD after framework superposition.

    Args:
        pred_path: Path to prediction structure (PDB or CIF)
        gt_path: Path to ground truth PDB
        cdr_info: Dict with cdr1_h, cdr2_h, cdr3_h, cdr1_l, cdr2_l, cdr3_l
                  (1-indexed residue positions within each chain)

    Returns:
        Dict with cdr3_h_rmsd, cdr3_l_rmsd, cdr3_combined_rmsd, framework_rmsd
    """
    pred_struct = parse_structure(pred_path)
    gt_struct = parse_structure(gt_path)

    pred_model = pred_struct[0]
    gt_model = gt_struct[0]

    # Get chains
    pred_heavy = pred_model[CHAINS_HEAVY]
    pred_light = pred_model[CHAINS_LIGHT]
    gt_heavy = gt_model[CHAINS_HEAVY]
    gt_light = gt_model[CHAINS_LIGHT]

    # Align residues by sequence
    heavy_pairs = align_chains_by_sequence(pred_heavy, gt_heavy)
    light_pairs = align_chains_by_sequence(pred_light, gt_light)

    if not heavy_pairs or not light_pairs:
        return {"cdr3_h_rmsd": np.nan, "cdr3_l_rmsd": np.nan,
                "cdr3_combined_rmsd": np.nan, "framework_rmsd": np.nan}

    # Build index mapping: 1-indexed position in prediction chain -> pair index
    # Prediction residues are numbered 1..N (per-chain numbering)
    pred_heavy_seq, pred_heavy_residues = get_chain_sequence(pred_heavy)
    pred_light_seq, pred_light_residues = get_chain_sequence(pred_light)

    # Map: 1-indexed position -> index in heavy_pairs / light_pairs
    # For each pair, determine what 1-indexed position the pred residue has
    heavy_pos_to_pair = {}
    for pair_idx, (pred_res, gt_res) in enumerate(heavy_pairs):
        # Find 1-indexed position of pred_res in pred_heavy_residues
        try:
            pos_in_chain = pred_heavy_residues.index(pred_res) + 1  # 1-indexed
            heavy_pos_to_pair[pos_in_chain] = pair_idx
        except ValueError:
            continue

    light_pos_to_pair = {}
    for pair_idx, (pred_res, gt_res) in enumerate(light_pairs):
        try:
            pos_in_chain = pred_light_residues.index(pred_res) + 1
            light_pos_to_pair[pos_in_chain] = pair_idx
        except ValueError:
            continue

    # Collect all CDR indices (for identifying framework = non-CDR)
    all_cdr_h = set(cdr_info["cdr1_h"] + cdr_info["cdr2_h"] + cdr_info["cdr3_h"])
    all_cdr_l = set(cdr_info["cdr1_l"] + cdr_info["cdr2_l"] + cdr_info["cdr3_l"])
    cdr3_h_set = set(cdr_info["cdr3_h"])
    cdr3_l_set = set(cdr_info["cdr3_l"])

    # Collect framework and CDR3 CA atoms
    framework_pred_atoms = []
    framework_gt_atoms = []
    cdr3_h_pred_atoms = []
    cdr3_h_gt_atoms = []
    cdr3_l_pred_atoms = []
    cdr3_l_gt_atoms = []

    # Heavy chain
    for pos, pair_idx in heavy_pos_to_pair.items():
        pred_res, gt_res = heavy_pairs[pair_idx]
        ca_pred = get_ca_atom(pred_res)
        ca_gt = get_ca_atom(gt_res)
        if ca_pred is None or ca_gt is None:
            continue

        if pos in cdr3_h_set:
            cdr3_h_pred_atoms.append(ca_pred)
            cdr3_h_gt_atoms.append(ca_gt)
        elif pos not in all_cdr_h:
            framework_pred_atoms.append(ca_pred)
            framework_gt_atoms.append(ca_gt)

    # Light chain
    for pos, pair_idx in light_pos_to_pair.items():
        pred_res, gt_res = light_pairs[pair_idx]
        ca_pred = get_ca_atom(pred_res)
        ca_gt = get_ca_atom(gt_res)
        if ca_pred is None or ca_gt is None:
            continue

        if pos in cdr3_l_set:
            cdr3_l_pred_atoms.append(ca_pred)
            cdr3_l_gt_atoms.append(ca_gt)
        elif pos not in all_cdr_l:
            framework_pred_atoms.append(ca_pred)
            framework_gt_atoms.append(ca_gt)

    if len(framework_pred_atoms) < 10:
        return {"cdr3_h_rmsd": np.nan, "cdr3_l_rmsd": np.nan,
                "cdr3_combined_rmsd": np.nan, "framework_rmsd": np.nan}

    # Kabsch alignment on framework
    sup = Superimposer()
    pred_fw_coords = np.array([a.get_vector().get_array() for a in framework_pred_atoms])
    gt_fw_coords = np.array([a.get_vector().get_array() for a in framework_gt_atoms])
    sup.set_atoms(framework_gt_atoms, framework_pred_atoms)
    framework_rmsd = sup.rms

    # Apply rotation/translation to CDR3 atoms
    rot, tran = sup.rotran

    def transform_and_rmsd(pred_atoms, gt_atoms):
        if not pred_atoms:
            return np.nan
        pred_coords = np.array([a.get_vector().get_array() for a in pred_atoms])
        gt_coords = np.array([a.get_vector().get_array() for a in gt_atoms])
        transformed = np.dot(pred_coords, rot) + tran
        diff = transformed - gt_coords
        return np.sqrt(np.mean(np.sum(diff ** 2, axis=1)))

    cdr3_h_rmsd = transform_and_rmsd(cdr3_h_pred_atoms, cdr3_h_gt_atoms)
    cdr3_l_rmsd = transform_and_rmsd(cdr3_l_pred_atoms, cdr3_l_gt_atoms)

    # Combined CDR3 RMSD
    all_cdr3_pred = cdr3_h_pred_atoms + cdr3_l_pred_atoms
    all_cdr3_gt = cdr3_h_gt_atoms + cdr3_l_gt_atoms
    cdr3_combined_rmsd = transform_and_rmsd(all_cdr3_pred, all_cdr3_gt)

    return {
        "cdr3_h_rmsd": cdr3_h_rmsd,
        "cdr3_l_rmsd": cdr3_l_rmsd,
        "cdr3_combined_rmsd": cdr3_combined_rmsd,
        "framework_rmsd": framework_rmsd,
    }


# ---------------------------------------------------------------------------
# DockQ via subprocess
# ---------------------------------------------------------------------------

def run_dockq(pred_path, gt_path):
    """Run DockQ and return parsed results dict.

    Returns dict with keys like 'AB', 'AC', 'BC' each containing DockQ metrics,
    plus 'global_dockq' (mean DockQ across interfaces).
    """
    with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as f:
        json_path = f.name

    try:
        cmd = [
            "conda", "run", "-n", DOCKQ_CONDA_ENV,
            "DockQ", str(pred_path), str(gt_path), "--json", json_path,
        ]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)

        if not os.path.exists(json_path) or os.path.getsize(json_path) == 0:
            print(f"  DockQ failed for {pred_path}: {result.stderr[:200]}")
            return None

        with open(json_path) as f:
            data = json.load(f)

        # Extract per-interface results
        best = data.get("best_result", {})
        output = {}
        dockq_values = []
        for interface_key, metrics in best.items():
            output[interface_key] = {
                "DockQ": metrics.get("DockQ", np.nan),
                "iRMSD": metrics.get("iRMSD", np.nan),
                "LRMSD": metrics.get("LRMSD", np.nan),
                "fnat": metrics.get("fnat", np.nan),
                "F1": metrics.get("F1", np.nan),
                "fnonnat": metrics.get("fnonnat", np.nan),
            }
            dockq_values.append(metrics.get("DockQ", 0.0))

        output["global_dockq"] = np.mean(dockq_values) if dockq_values else np.nan
        return output

    except (subprocess.TimeoutExpired, Exception) as e:
        print(f"  DockQ error for {pred_path}: {e}")
        return None
    finally:
        if os.path.exists(json_path):
            os.unlink(json_path)


# ---------------------------------------------------------------------------
# Confidence metrics
# ---------------------------------------------------------------------------

def extract_confidence(pred_dir, prefix):
    """Extract confidence metrics from JSON and pLDDT from npz.

    Args:
        pred_dir: Directory containing confidence_*.json and plddt_*.npz
        prefix: File prefix (e.g., '7TRH_HBG_model_0' or '7TRH_HBG_cdr3_beta_model_0')

    Returns:
        Dict with confidence metrics.
    """
    result = {}

    # Confidence JSON
    json_path = pred_dir / f"confidence_{prefix}.json"
    if json_path.exists():
        with open(json_path) as f:
            data = json.load(f)
        result["confidence_score"] = data.get("confidence_score", np.nan)
        result["iptm"] = data.get("iptm", np.nan)
        result["ptm"] = data.get("ptm", np.nan)
        result["complex_plddt"] = data.get("complex_plddt", np.nan)
        result["complex_iplddt"] = data.get("complex_iplddt", np.nan)
        result["protein_iptm"] = data.get("protein_iptm", np.nan)

    # pLDDT npz
    plddt_path = pred_dir / f"plddt_{prefix}.npz"
    if plddt_path.exists():
        plddt_data = np.load(plddt_path)
        plddt = plddt_data["plddt"]
        result["mean_plddt"] = float(np.mean(plddt))

    return result


def extract_cdr3_plddt(pred_dir, prefix, cdr_info, pred_path):
    """Extract per-CDR3 pLDDT values.

    The pLDDT array has one value per residue across all chains in order.
    We need to know the chain lengths to index into it.
    """
    plddt_path = pred_dir / f"plddt_{prefix}.npz"
    if not plddt_path.exists():
        return {}

    plddt_data = np.load(plddt_path)
    plddt = plddt_data["plddt"]

    # Get chain lengths from the structure
    struct = parse_structure(pred_path)
    model = struct[0]
    chain_lengths = {}
    for chain in model:
        _, residues = get_chain_sequence(chain)
        chain_lengths[chain.id] = len(residues)

    # Compute offsets
    chain_order = sorted(chain_lengths.keys())
    offsets = {}
    offset = 0
    for cid in chain_order:
        offsets[cid] = offset
        offset += chain_lengths[cid]

    result = {}

    # CDR3-H pLDDT (chain B)
    if CHAINS_HEAVY in offsets:
        h_offset = offsets[CHAINS_HEAVY]
        cdr3_h_indices = [h_offset + (i - 1) for i in cdr_info["cdr3_h"]]  # 1-indexed to 0-indexed
        valid_h = [i for i in cdr3_h_indices if i < len(plddt)]
        if valid_h:
            result["cdr3_h_plddt"] = float(np.mean(plddt[valid_h]))

    # CDR3-L pLDDT (chain C)
    if CHAINS_LIGHT in offsets:
        l_offset = offsets[CHAINS_LIGHT]
        cdr3_l_indices = [l_offset + (i - 1) for i in cdr_info["cdr3_l"]]
        valid_l = [i for i in cdr3_l_indices if i < len(plddt)]
        if valid_l:
            result["cdr3_l_plddt"] = float(np.mean(plddt[valid_l]))

    # Antibody pLDDT (chains B+C)
    ab_indices = []
    for cid in [CHAINS_HEAVY, CHAINS_LIGHT]:
        if cid in offsets:
            ab_indices.extend(range(offsets[cid], offsets[cid] + chain_lengths[cid]))
    valid_ab = [i for i in ab_indices if i < len(plddt)]
    if valid_ab:
        result["antibody_plddt"] = float(np.mean(plddt[valid_ab]))

    # Antigen pLDDT (chain A)
    if CHAINS_ANTIGEN in offsets:
        ag_indices = list(range(offsets[CHAINS_ANTIGEN],
                                offsets[CHAINS_ANTIGEN] + chain_lengths[CHAINS_ANTIGEN]))
        valid_ag = [i for i in ag_indices if i < len(plddt)]
        if valid_ag:
            result["antigen_plddt"] = float(np.mean(plddt[valid_ag]))

    return result


# ---------------------------------------------------------------------------
# Epitope contact analysis (Tier 4)
# ---------------------------------------------------------------------------

def compute_contacts(structure, chain_ab_ids, chain_ag_id, cutoff=5.0):
    """Compute antibody-antigen contacts (heavy atom distance <= cutoff).

    Returns set of (ag_resid, ab_chain, ab_resid) tuples.
    """
    model = structure[0]
    ag_chain = model[chain_ag_id]
    contacts = set()

    ag_residues = [r for r in ag_chain.get_residues() if r.id[0] == " "]

    for ab_chain_id in chain_ab_ids:
        ab_chain = model[ab_chain_id]
        ab_residues = [r for r in ab_chain.get_residues() if r.id[0] == " "]

        for ag_res in ag_residues:
            ag_heavy = [a for a in ag_res.get_atoms() if a.element != "H"]
            for ab_res in ab_residues:
                ab_heavy = [a for a in ab_res.get_atoms() if a.element != "H"]
                min_dist = float("inf")
                for a1 in ag_heavy:
                    for a2 in ab_heavy:
                        d = a1 - a2
                        if d < min_dist:
                            min_dist = d
                            if d <= cutoff:
                                break
                    if min_dist <= cutoff:
                        break
                if min_dist <= cutoff:
                    # Use sequence position as identifier (chain-agnostic)
                    contacts.add((chain_ag_id, ag_res.id[1], ab_chain_id, ab_res.id[1]))

    return contacts


def contact_metrics(pred_contacts, gt_contacts):
    """Compute precision, recall, F1 for contact prediction.

    We compare at the (ag_resid, ab_chain, ab_resid) level.
    Since numbering differs, we compare by position-within-chain index.
    """
    if not gt_contacts:
        return {"contact_precision": np.nan, "contact_recall": np.nan, "contact_f1": np.nan,
                "epitope_precision": np.nan, "epitope_recall": np.nan}

    # Contact-level
    tp = len(pred_contacts & gt_contacts)
    fp = len(pred_contacts - gt_contacts)
    fn = len(gt_contacts - pred_contacts)

    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0.0

    # Epitope residue-level (antigen side only)
    pred_epitope = {(c[0], c[1]) for c in pred_contacts}
    gt_epitope = {(c[0], c[1]) for c in gt_contacts}
    ep_tp = len(pred_epitope & gt_epitope)
    ep_fp = len(pred_epitope - gt_epitope)
    ep_fn = len(gt_epitope - pred_epitope)

    ep_precision = ep_tp / (ep_tp + ep_fp) if (ep_tp + ep_fp) > 0 else 0.0
    ep_recall = ep_tp / (ep_tp + ep_fn) if (ep_tp + ep_fn) > 0 else 0.0

    return {
        "contact_precision": precision,
        "contact_recall": recall,
        "contact_f1": f1,
        "epitope_precision": ep_precision,
        "epitope_recall": ep_recall,
    }


def compute_contact_analysis(pred_path, gt_path):
    """Compute epitope contact analysis between prediction and ground truth.

    GT and prediction chains may have different lengths (GT includes constant
    regions). We align chains by sequence to build a mapping from GT position
    indices to prediction position indices, then compare contacts using the
    prediction's position space.
    """
    pred_struct = parse_structure(pred_path)
    gt_struct = parse_structure(gt_path)

    pred_model = pred_struct[0]
    gt_model = gt_struct[0]

    # Build GT position -> pred position mapping per chain via sequence alignment
    gt_to_pred_pos = {}  # (chain_id, gt_pos) -> pred_pos
    for chain_id in [CHAINS_ANTIGEN, CHAINS_HEAVY, CHAINS_LIGHT]:
        if chain_id not in pred_model or chain_id not in gt_model:
            continue
        pairs = align_chains_by_sequence(pred_model[chain_id], gt_model[chain_id])
        pred_seq, pred_residues = get_chain_sequence(pred_model[chain_id])
        gt_seq, gt_residues = get_chain_sequence(gt_model[chain_id])

        for pred_res, gt_res in pairs:
            try:
                pred_pos = pred_residues.index(pred_res)
                gt_pos = gt_residues.index(gt_res)
                gt_to_pred_pos[(chain_id, gt_pos)] = pred_pos
            except ValueError:
                continue

    # Get contacts using raw residue IDs
    pred_raw = compute_contacts(pred_struct, [CHAINS_HEAVY, CHAINS_LIGHT], CHAINS_ANTIGEN)
    gt_raw = compute_contacts(gt_struct, [CHAINS_HEAVY, CHAINS_LIGHT], CHAINS_ANTIGEN)

    # Normalize prediction contacts to 0-indexed position-within-chain
    def normalize_pred_contacts(raw_contacts, model):
        chain_maps = {}
        for chain in model:
            residues = [r for r in chain.get_residues() if r.id[0] == " "]
            chain_maps[chain.id] = {r.id[1]: i for i, r in enumerate(residues)}
        normalized = set()
        for ag_cid, ag_rid, ab_cid, ab_rid in raw_contacts:
            ag_pos = chain_maps.get(ag_cid, {}).get(ag_rid)
            ab_pos = chain_maps.get(ab_cid, {}).get(ab_rid)
            if ag_pos is not None and ab_pos is not None:
                normalized.add((ag_cid, ag_pos, ab_cid, ab_pos))
        return normalized

    pred_norm = normalize_pred_contacts(pred_raw, pred_model)

    # Normalize GT contacts by mapping GT positions -> pred positions
    gt_chain_maps = {}
    for chain in gt_model:
        residues = [r for r in chain.get_residues() if r.id[0] == " "]
        gt_chain_maps[chain.id] = {r.id[1]: i for i, r in enumerate(residues)}

    gt_norm = set()
    for ag_cid, ag_rid, ab_cid, ab_rid in gt_raw:
        ag_gt_pos = gt_chain_maps.get(ag_cid, {}).get(ag_rid)
        ab_gt_pos = gt_chain_maps.get(ab_cid, {}).get(ab_rid)
        if ag_gt_pos is None or ab_gt_pos is None:
            continue
        ag_pred_pos = gt_to_pred_pos.get((ag_cid, ag_gt_pos))
        ab_pred_pos = gt_to_pred_pos.get((ab_cid, ab_gt_pos))
        if ag_pred_pos is not None and ab_pred_pos is not None:
            gt_norm.add((ag_cid, ag_pred_pos, ab_cid, ab_pred_pos))

    return contact_metrics(pred_norm, gt_norm)


# ---------------------------------------------------------------------------
# Discovery: find prediction files
# ---------------------------------------------------------------------------

def discover_predictions_b1(pred_base, complexes):
    """Discover B1 (vanilla Boltz2) predictions."""
    folder = pred_base / METHODS["B1"]["folder"]
    results = {}

    for cname in complexes:
        result_dir = folder / f"boltz_results_{cname}"
        pred_dir = result_dir / "predictions" / cname
        if not pred_dir.exists():
            continue

        models = []
        for i in range(5):
            pdb_path = pred_dir / f"{cname}_model_{i}.pdb"
            if pdb_path.exists():
                models.append({
                    "path": pdb_path,
                    "prefix": f"{cname}_model_{i}",
                    "pred_dir": pred_dir,
                    "model_idx": i,
                })
        if models:
            results[cname] = models

    return results


def discover_predictions_v(pred_base, complexes):
    """Discover V (embedding steering) predictions."""
    folder = pred_base / METHODS["V"]["folder"]
    results = {}

    for cname in complexes:
        result_dir = folder / f"boltz_results_{cname}_cdr3_beta"
        pred_dir = result_dir / "predictions" / f"{cname}_cdr3_beta"
        if not pred_dir.exists():
            continue

        models = []
        for i in range(5):
            # Try both .cif and .pdb
            for ext in [".cif", ".pdb"]:
                path = pred_dir / f"{cname}_cdr3_beta_model_{i}{ext}"
                if path.exists():
                    models.append({
                        "path": path,
                        "prefix": f"{cname}_cdr3_beta_model_{i}",
                        "pred_dir": pred_dir,
                        "model_idx": i,
                    })
                    break
        if models:
            results[cname] = models

    return results


def discover_predictions_multi(pred_base, method_key, complexes):
    """Discover predictions for multi-setting baselines (B2, B3)."""
    folder = pred_base / METHODS[method_key]["folder"]
    if not folder.exists():
        return {}

    results = {}  # complex -> list of {setting, models: [...]}

    for cname in complexes:
        settings = []
        # Find all result folders for this complex
        if method_key == "B2":
            prefix = f"boltz_results_restraint_{cname}_"
        else:  # B3
            prefix = f"boltz_results_restraint_to_A_{cname}_"

        for entry in sorted(folder.iterdir()):
            if not entry.is_dir() or not entry.name.startswith(prefix.rstrip("_")):
                continue

            # Check if this specific entry matches the complex
            # For B2: boltz_results_restraint_{complex}_{type}_{num}
            # For B3: boltz_results_restraint_to_A_{complex}_{chain}_{res}_{num}
            setting_name = entry.name

            # Find predictions subfolder
            pred_subdir = entry / "predictions"
            if not pred_subdir.exists():
                continue

            # The inner folder name varies
            inner_dirs = [d for d in pred_subdir.iterdir() if d.is_dir()]
            if not inner_dirs:
                continue

            for inner_dir in inner_dirs:
                models = []
                for i in range(5):
                    for ext in [".pdb", ".cif"]:
                        path = inner_dir / f"{inner_dir.name}_model_{i}{ext}"
                        if path.exists():
                            models.append({
                                "path": path,
                                "prefix": f"{inner_dir.name}_model_{i}",
                                "pred_dir": inner_dir,
                                "model_idx": i,
                                "setting": setting_name,
                            })
                            break
                if models:
                    settings.append({"setting": setting_name, "models": models})

        if settings:
            results[cname] = settings

    return results


# ---------------------------------------------------------------------------
# Single prediction evaluation
# ---------------------------------------------------------------------------

def evaluate_single_prediction(pred_info, gt_path, cdr_info, run_dockq_flag=True,
                                run_cdr3_flag=True, run_contacts_flag=True):
    """Evaluate a single prediction against ground truth.

    Returns dict with all metrics.
    """
    pred_path = pred_info["path"]
    pred_dir = pred_info["pred_dir"]
    prefix = pred_info["prefix"]

    result = {
        "pred_path": str(pred_path),
        "model_idx": pred_info["model_idx"],
    }

    # Tier 1: DockQ
    if run_dockq_flag:
        dockq = run_dockq(pred_path, gt_path)
        if dockq:
            result["global_dockq"] = dockq.get("global_dockq", np.nan)
            for iface in ["AB", "AC", "BC"]:
                if iface in dockq:
                    for metric in ["DockQ", "iRMSD", "LRMSD", "fnat", "F1", "fnonnat"]:
                        result[f"{iface}_{metric}"] = dockq[iface].get(metric, np.nan)

    # Tier 2: CDR3 RMSD
    if run_cdr3_flag and cdr_info:
        try:
            cdr3_results = compute_cdr3_rmsd(pred_path, gt_path, cdr_info)
            result.update(cdr3_results)
        except Exception as e:
            print(f"  CDR3 RMSD failed for {pred_path}: {e}")
            result.update({"cdr3_h_rmsd": np.nan, "cdr3_l_rmsd": np.nan,
                           "cdr3_combined_rmsd": np.nan, "framework_rmsd": np.nan})

    # Tier 3: Confidence
    conf = extract_confidence(pred_dir, prefix)
    result.update(conf)

    if cdr_info:
        try:
            cdr3_plddt = extract_cdr3_plddt(pred_dir, prefix, cdr_info, pred_path)
            result.update(cdr3_plddt)
        except Exception as e:
            print(f"  CDR3 pLDDT extraction failed for {pred_path}: {e}")

    # Tier 4: Contact analysis
    if run_contacts_flag:
        try:
            contact_results = compute_contact_analysis(pred_path, gt_path)
            result.update(contact_results)
        except Exception as e:
            print(f"  Contact analysis failed for {pred_path}: {e}")

    return result


# ---------------------------------------------------------------------------
# Model selection
# ---------------------------------------------------------------------------

def select_best_model(model_results):
    """Select best model by confidence_score."""
    best = None
    best_score = -1
    for r in model_results:
        score = r.get("confidence_score", -1)
        if not np.isnan(score) and score > best_score:
            best_score = score
            best = r
    return best if best is not None else (model_results[0] if model_results else None)


# ---------------------------------------------------------------------------
# CAPRI classification
# ---------------------------------------------------------------------------

def capri_category(dockq):
    """Classify DockQ score into CAPRI category."""
    if np.isnan(dockq):
        return "N/A"
    for cat, (lo, hi) in CAPRI_THRESHOLDS.items():
        if lo <= dockq < hi:
            return cat
    return "incorrect"


# ---------------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------------

def compute_statistics(method_a_values, method_b_values, metric_name):
    """Compute comparison statistics between two methods."""
    # Filter to complexes with both values
    pairs = [(a, b) for a, b in zip(method_a_values, method_b_values)
             if not np.isnan(a) and not np.isnan(b)]
    if not pairs:
        return {}

    a_vals, b_vals = zip(*pairs)
    a_vals, b_vals = np.array(a_vals), np.array(b_vals)
    deltas = a_vals - b_vals

    result = {
        "n_pairs": len(pairs),
        "mean_a": float(np.mean(a_vals)),
        "mean_b": float(np.mean(b_vals)),
        "mean_delta": float(np.mean(deltas)),
        "median_delta": float(np.median(deltas)),
        "std_delta": float(np.std(deltas)),
    }

    # For RMSD metrics, lower is better (wins = a < b)
    # For DockQ/confidence, higher is better (wins = a > b)
    higher_is_better = "rmsd" not in metric_name.lower()
    if higher_is_better:
        result["wins"] = int(np.sum(deltas > 0.001))
        result["losses"] = int(np.sum(deltas < -0.001))
    else:
        result["wins"] = int(np.sum(deltas < -0.001))
        result["losses"] = int(np.sum(deltas > 0.001))
    result["ties"] = len(pairs) - result["wins"] - result["losses"]

    # Wilcoxon signed-rank test (paired, non-parametric)
    if len(pairs) >= 5:
        try:
            stat, pval = stats.wilcoxon(a_vals, b_vals, alternative="two-sided")
            result["wilcoxon_stat"] = float(stat)
            result["wilcoxon_pval"] = float(pval)
        except ValueError:
            pass

    return result


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------

def format_val(v, fmt=".3f"):
    """Format a numeric value for display."""
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return "—"
    return f"{v:{fmt}}"


def write_summary_table(results_by_method, complexes, out_dir):
    """Write primary results summary table."""
    # Metrics to report
    metrics = [
        ("global_dockq", "DockQ (Global)", ".3f"),
        ("AB_DockQ", "DockQ (A-B)", ".3f"),
        ("AC_DockQ", "DockQ (A-C)", ".3f"),
        ("BC_DockQ", "DockQ (B-C)", ".3f"),
        ("cdr3_h_rmsd", "CDR3-H RMSD", ".2f"),
        ("cdr3_l_rmsd", "CDR3-L RMSD", ".2f"),
        ("cdr3_combined_rmsd", "CDR3 Combined RMSD", ".2f"),
        ("confidence_score", "Confidence", ".3f"),
        ("iptm", "ipTM", ".3f"),
        ("ptm", "pTM", ".3f"),
        ("complex_plddt", "Complex pLDDT", ".3f"),
        ("cdr3_h_plddt", "CDR3-H pLDDT", ".3f"),
        ("cdr3_l_plddt", "CDR3-L pLDDT", ".3f"),
        ("contact_f1", "Contact F1", ".3f"),
        ("epitope_recall", "Epitope Recall", ".3f"),
    ]

    lines = []
    lines.append("=" * 100)
    lines.append("EVALUATION RESULTS: Strategy V (Embedding-Space CDR3 Steering)")
    lines.append("=" * 100)
    lines.append("")

    # Aggregate table
    header = f"{'Metric':<25}"
    for method in results_by_method:
        header += f" | {method:>20}"
    lines.append(header)
    lines.append("-" * len(header))

    for key, label, fmt in metrics:
        row = f"{label:<25}"
        for method, best_results in results_by_method.items():
            values = [r.get(key, np.nan) for r in best_results.values() if r]
            values = [v for v in values if v is not None and not np.isnan(v)]
            if values:
                mean = np.mean(values)
                sem = np.std(values) / np.sqrt(len(values)) if len(values) > 1 else 0
                row += f" | {mean:{fmt}} +/- {sem:{fmt}}"
            else:
                row += f" | {'—':>20}"
        lines.append(row)

    lines.append("")

    # CAPRI distribution
    lines.append("CAPRI Quality Distribution (DockQ Global, best-by-confidence):")
    lines.append("-" * 60)
    for method, best_results in results_by_method.items():
        cats = defaultdict(int)
        for r in best_results.values():
            if r:
                cats[capri_category(r.get("global_dockq", np.nan))] += 1
        cat_str = ", ".join(f"{k}: {v}" for k, v in sorted(cats.items()))
        lines.append(f"  {method}: {cat_str}")

    lines.append("")

    # Statistical comparison (V vs each baseline)
    if "V" in results_by_method:
        v_results = results_by_method["V"]
        for baseline in results_by_method:
            if baseline == "V":
                continue
            b_results = results_by_method[baseline]

            lines.append(f"Statistical Comparison: V vs {baseline}")
            lines.append("-" * 60)

            common = sorted(set(v_results.keys()) & set(b_results.keys()))
            if not common:
                lines.append("  No common complexes")
                continue

            for key, label, fmt in metrics:
                v_vals = [v_results[c].get(key, np.nan) if v_results[c] else np.nan for c in common]
                b_vals = [b_results[c].get(key, np.nan) if b_results[c] else np.nan for c in common]
                stat = compute_statistics(v_vals, b_vals, key)
                if stat:
                    pval_str = f"p={stat.get('wilcoxon_pval', np.nan):.4f}" if "wilcoxon_pval" in stat else "N/A"
                    lines.append(
                        f"  {label:<25}: delta={stat['mean_delta']:+{fmt}}  "
                        f"W/T/L={stat['wins']}/{stat['ties']}/{stat['losses']}  "
                        f"{pval_str}"
                    )
            lines.append("")

    summary_text = "\n".join(lines)
    print(summary_text)

    summary_path = out_dir / "summary.txt"
    with open(summary_path, "w") as f:
        f.write(summary_text)
    print(f"\nSummary written to {summary_path}")


def write_per_complex_csv(results_by_method, complexes, out_dir):
    """Write detailed per-complex CSV results."""
    # Collect all metric keys
    all_keys = set()
    for method_results in results_by_method.values():
        for r in method_results.values():
            if r:
                all_keys.update(r.keys())
    all_keys -= {"pred_path", "model_idx"}
    metric_keys = sorted(all_keys)

    csv_path = out_dir / "per_complex_results.csv"
    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)

        # Header
        header = ["complex", "method"] + metric_keys
        writer.writerow(header)

        for cname in sorted(complexes):
            for method, best_results in results_by_method.items():
                r = best_results.get(cname)
                if r is None:
                    continue
                row = [cname, method]
                for key in metric_keys:
                    val = r.get(key, "")
                    if isinstance(val, float) and np.isnan(val):
                        val = ""
                    row.append(val)
                writer.writerow(row)

    print(f"Per-complex CSV written to {csv_path}")


def write_all_models_csv(all_results, out_dir):
    """Write CSV with all models (not just best) for secondary analysis."""
    all_keys = set()
    for method_results in all_results.values():
        for cname, models in method_results.items():
            for r in models:
                all_keys.update(r.keys())
    all_keys -= {"pred_path"}
    metric_keys = sorted(all_keys)

    csv_path = out_dir / "all_models_results.csv"
    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        header = ["complex", "method"] + metric_keys
        writer.writerow(header)

        for method, method_results in all_results.items():
            for cname in sorted(method_results.keys()):
                for r in method_results[cname]:
                    row = [cname, method]
                    for key in metric_keys:
                        val = r.get(key, "")
                        if isinstance(val, float) and np.isnan(val):
                            val = ""
                        row.append(val)
                    writer.writerow(row)

    print(f"All-models CSV written to {csv_path}")


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def generate_plots(results_by_method, complexes, out_dir):
    """Generate evaluation plots."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available, skipping plots")
        return

    plot_dir = out_dir / "plots"
    plot_dir.mkdir(exist_ok=True)

    methods = list(results_by_method.keys())
    colors = {"B1": "#1f77b4", "B2": "#ff7f0e", "B3": "#2ca02c", "V": "#d62728"}

    # --- Plot 1: DockQ box plots ---
    fig, axes = plt.subplots(1, 4, figsize=(16, 4), sharey=True)
    for ax, (interface, title) in zip(axes, [
        ("global_dockq", "Global"), ("AB_DockQ", "A-B (Ag-Heavy)"),
        ("AC_DockQ", "A-C (Ag-Light)"), ("BC_DockQ", "B-C (Heavy-Light)")
    ]):
        data = []
        labels = []
        for m in methods:
            vals = [results_by_method[m].get(c, {}).get(interface, np.nan)
                    if results_by_method[m].get(c) else np.nan
                    for c in complexes]
            vals = [v for v in vals if not np.isnan(v)]
            if vals:
                data.append(vals)
                labels.append(m)
        if data:
            bp = ax.boxplot(data, labels=labels, patch_artist=True)
            for patch, label in zip(bp["boxes"], labels):
                patch.set_facecolor(colors.get(label, "#999999"))
                patch.set_alpha(0.7)
        ax.set_title(title)
        ax.set_ylabel("DockQ" if ax == axes[0] else "")
        ax.axhline(y=0.23, color="gray", linestyle="--", alpha=0.5, linewidth=0.8)
        ax.axhline(y=0.49, color="gray", linestyle="--", alpha=0.5, linewidth=0.8)
        ax.axhline(y=0.80, color="gray", linestyle="--", alpha=0.5, linewidth=0.8)
        ax.set_ylim(0, 1.05)
    plt.tight_layout()
    plt.savefig(plot_dir / "dockq_boxplots.png", dpi=150)
    plt.close()

    # --- Plot 2: Scatter V vs B1 ---
    if "V" in results_by_method and "B1" in results_by_method:
        fig, ax = plt.subplots(figsize=(6, 6))
        common = sorted(set(results_by_method["V"].keys()) & set(results_by_method["B1"].keys()))
        x_vals, y_vals = [], []
        for c in common:
            v_r = results_by_method["V"].get(c)
            b_r = results_by_method["B1"].get(c)
            if v_r and b_r:
                xv = b_r.get("global_dockq", np.nan)
                yv = v_r.get("global_dockq", np.nan)
                if not np.isnan(xv) and not np.isnan(yv):
                    x_vals.append(xv)
                    y_vals.append(yv)

        if x_vals:
            ax.scatter(x_vals, y_vals, alpha=0.7, s=30, color=colors["V"])
            ax.plot([0, 1], [0, 1], "k--", alpha=0.3)
            ax.set_xlabel("B1 (Vanilla) DockQ")
            ax.set_ylabel("V (Embed Steer) DockQ")
            ax.set_title("DockQ: V vs B1")
            ax.set_xlim(0, 1.05)
            ax.set_ylim(0, 1.05)
            ax.set_aspect("equal")
            n_above = sum(1 for x, y in zip(x_vals, y_vals) if y > x + 0.001)
            n_below = sum(1 for x, y in zip(x_vals, y_vals) if y < x - 0.001)
            ax.text(0.05, 0.95, f"V wins: {n_above}\nB1 wins: {n_below}",
                    transform=ax.transAxes, va="top", fontsize=9)
        plt.tight_layout()
        plt.savefig(plot_dir / "scatter_v_vs_b1.png", dpi=150)
        plt.close()

    # --- Plot 3: CAPRI category bar chart ---
    fig, ax = plt.subplots(figsize=(8, 4))
    cat_order = ["high", "medium", "acceptable", "incorrect"]
    cat_colors = {"high": "#2ca02c", "medium": "#98df8a", "acceptable": "#ffbb78", "incorrect": "#ff7f0e"}
    x = np.arange(len(methods))
    width = 0.6
    bottoms = np.zeros(len(methods))

    for cat in cat_order:
        heights = []
        for m in methods:
            total = sum(1 for r in results_by_method[m].values() if r)
            if total == 0:
                heights.append(0)
                continue
            count = sum(1 for r in results_by_method[m].values()
                        if r and capri_category(r.get("global_dockq", np.nan)) == cat)
            heights.append(count / total * 100)
        ax.bar(x, heights, width, bottom=bottoms, label=cat, color=cat_colors[cat])
        bottoms += heights

    ax.set_xticks(x)
    ax.set_xticklabels(methods)
    ax.set_ylabel("Percentage (%)")
    ax.set_title("CAPRI Quality Distribution")
    ax.legend(loc="upper right")
    plt.tight_layout()
    plt.savefig(plot_dir / "capri_distribution.png", dpi=150)
    plt.close()

    # --- Plot 4: CDR3 RMSD violin/box plot ---
    fig, axes = plt.subplots(1, 3, figsize=(12, 4), sharey=True)
    for ax, (metric, title) in zip(axes, [
        ("cdr3_h_rmsd", "CDR3-H"), ("cdr3_l_rmsd", "CDR3-L"),
        ("cdr3_combined_rmsd", "CDR3 Combined")
    ]):
        data = []
        labels = []
        for m in methods:
            vals = [results_by_method[m].get(c, {}).get(metric, np.nan)
                    if results_by_method[m].get(c) else np.nan
                    for c in complexes]
            vals = [v for v in vals if not np.isnan(v)]
            if vals:
                data.append(vals)
                labels.append(m)
        if data:
            bp = ax.boxplot(data, labels=labels, patch_artist=True)
            for patch, label in zip(bp["boxes"], labels):
                patch.set_facecolor(colors.get(label, "#999999"))
                patch.set_alpha(0.7)
        ax.set_title(title)
        ax.set_ylabel("RMSD (Å)" if ax == axes[0] else "")
    plt.tight_layout()
    plt.savefig(plot_dir / "cdr3_rmsd_boxplots.png", dpi=150)
    plt.close()

    # --- Plot 5: Stratified analysis by difficulty ---
    if "B1" in results_by_method and len(results_by_method) > 1:
        fig, ax = plt.subplots(figsize=(8, 4))
        difficulty_bins = {"easy": [], "medium": [], "hard": []}
        for c in complexes:
            b1_r = results_by_method["B1"].get(c)
            if not b1_r:
                continue
            dq = b1_r.get("global_dockq", np.nan)
            if np.isnan(dq):
                continue
            if dq >= 0.80:
                difficulty_bins["easy"].append(c)
            elif dq >= 0.49:
                difficulty_bins["medium"].append(c)
            else:
                difficulty_bins["hard"].append(c)

        diff_labels = ["easy", "medium", "hard"]
        x = np.arange(len(diff_labels))
        bar_width = 0.8 / len(methods)

        for i, m in enumerate(methods):
            means = []
            for diff in diff_labels:
                vals = [results_by_method[m].get(c, {}).get("global_dockq", np.nan)
                        if results_by_method[m].get(c) else np.nan
                        for c in difficulty_bins[diff]]
                vals = [v for v in vals if not np.isnan(v)]
                means.append(np.mean(vals) if vals else 0)
            ax.bar(x + i * bar_width, means, bar_width, label=m,
                   color=colors.get(m, "#999999"), alpha=0.8)

        ax.set_xticks(x + bar_width * (len(methods) - 1) / 2)
        ax.set_xticklabels([f"{d}\n(n={len(difficulty_bins[d])})" for d in diff_labels])
        ax.set_ylabel("Mean DockQ")
        ax.set_title("DockQ by Difficulty (based on B1 performance)")
        ax.legend()
        plt.tight_layout()
        plt.savefig(plot_dir / "stratified_difficulty.png", dpi=150)
        plt.close()

    print(f"Plots saved to {plot_dir}/")


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Evaluate embedding-space CDR3 steering")
    parser.add_argument("--gt_dir", type=Path, default=Path("pdb_minimized"),
                        help="Directory with ground truth PDB files")
    parser.add_argument("--pred_dir", type=Path, default=Path("predictions_examples"),
                        help="Base directory with prediction folders")
    parser.add_argument("--cdrs_csv", type=Path, default=Path("examples/cdrs.csv"),
                        help="CSV with CDR definitions per complex")
    parser.add_argument("--out_dir", type=Path, default=Path("eval_results_embedding_steering"),
                        help="Output directory for results")
    parser.add_argument("--methods", nargs="+", default=["B1", "V"],
                        help="Methods to evaluate (default: B1 V)")
    parser.add_argument("--skip_dockq", action="store_true",
                        help="Skip DockQ computation (use cached results)")
    parser.add_argument("--skip_contacts", action="store_true",
                        help="Skip contact analysis (Tier 4)")
    parser.add_argument("--skip_plots", action="store_true",
                        help="Skip plot generation")
    parser.add_argument("--complexes", nargs="+", default=None,
                        help="Specific complexes to evaluate (default: all available)")
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    # Load CDR definitions
    print("Loading CDR definitions...")
    cdrs = load_cdrs(args.cdrs_csv)

    # Discover ground truth
    gt_files = {p.stem: p for p in args.gt_dir.glob("*.pdb")}
    all_complexes = sorted(gt_files.keys())
    if args.complexes:
        all_complexes = [c for c in args.complexes if c in gt_files]
    print(f"Found {len(all_complexes)} ground truth structures")

    # Discover predictions per method
    print("\nDiscovering predictions...")
    discovered = {}
    for method in args.methods:
        if method == "B1":
            discovered["B1"] = discover_predictions_b1(args.pred_dir, all_complexes)
        elif method == "V":
            discovered["V"] = discover_predictions_v(args.pred_dir, all_complexes)
        elif method in ("B2", "B3"):
            multi = discover_predictions_multi(args.pred_dir, method, all_complexes)
            if multi:
                discovered[method] = multi
        print(f"  {method}: {len(discovered.get(method, {}))} complexes with predictions")

    # Evaluate each method
    all_results = {}  # method -> {complex -> [model_results]}
    best_results = {}  # method -> {complex -> best_model_result}

    for method in args.methods:
        if method not in discovered:
            continue

        print(f"\n{'='*60}")
        print(f"Evaluating {method}...")
        print(f"{'='*60}")

        method_all = {}
        method_best = {}

        method_data = discovered[method]

        if method in ("B2", "B3"):
            # Multi-setting: flatten all settings' models
            for cname, settings in method_data.items():
                gt_path = gt_files.get(cname)
                if not gt_path:
                    continue
                cdr_info = cdrs.get(cname)
                print(f"\n  {cname} ({len(settings)} settings)...")

                all_models = []
                for setting_info in settings:
                    for model_info in setting_info["models"]:
                        print(f"    Evaluating {model_info['path'].name}...")
                        r = evaluate_single_prediction(
                            model_info, gt_path, cdr_info,
                            run_dockq_flag=not args.skip_dockq,
                            run_contacts_flag=not args.skip_contacts,
                        )
                        r["setting"] = setting_info["setting"]
                        all_models.append(r)

                method_all[cname] = all_models
                method_best[cname] = select_best_model(all_models)
        else:
            # Single-setting: list of models per complex
            for cname, models in method_data.items():
                gt_path = gt_files.get(cname)
                if not gt_path:
                    continue
                cdr_info = cdrs.get(cname)
                print(f"\n  {cname} ({len(models)} models)...")

                model_results = []
                for model_info in models:
                    print(f"    Evaluating {model_info['path'].name}...")
                    r = evaluate_single_prediction(
                        model_info, gt_path, cdr_info,
                        run_dockq_flag=not args.skip_dockq,
                        run_contacts_flag=not args.skip_contacts,
                    )
                    model_results.append(r)

                method_all[cname] = model_results
                method_best[cname] = select_best_model(model_results)

        all_results[method] = method_all
        best_results[method] = method_best

    # Reporting
    print(f"\n{'='*60}")
    print("RESULTS")
    print(f"{'='*60}\n")

    complexes_with_results = set()
    for m_results in best_results.values():
        complexes_with_results.update(m_results.keys())
    complexes_with_results = sorted(complexes_with_results)

    write_summary_table(best_results, complexes_with_results, args.out_dir)
    write_per_complex_csv(best_results, complexes_with_results, args.out_dir)
    write_all_models_csv(all_results, args.out_dir)

    # Cache raw results as JSON for reuse
    cache_path = args.out_dir / "raw_results.json"
    serializable = {}
    for method, method_results in all_results.items():
        serializable[method] = {}
        for cname, models in method_results.items():
            serializable[method][cname] = []
            for r in models:
                sr = {}
                for k, v in r.items():
                    if isinstance(v, (np.floating, np.integer)):
                        sr[k] = float(v)
                    elif isinstance(v, Path):
                        sr[k] = str(v)
                    else:
                        sr[k] = v
                serializable[method][cname].append(sr)
    with open(cache_path, "w") as f:
        json.dump(serializable, f, indent=2, default=str)
    print(f"Raw results cached to {cache_path}")

    # Plots
    if not args.skip_plots:
        generate_plots(best_results, complexes_with_results, args.out_dir)

    print("\nEvaluation complete!")


if __name__ == "__main__":
    main()
