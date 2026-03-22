"""Iterative Epitope Refinement (Strategy Q+ v2).

Multi-round epitope discovery pipeline:
- Round 1: Broad exploration with many samples → contact heatmap
- Round 2 (optional): Focused exploration on hotspot residues
- Round 3: Final refinement with tight pocket constraints

Implements Q+ improvements:
- Q+.1: More Round 1 diversity (configurable samples + noise scale)
- Q+.2: Confidence-weighted contact maps
- Q+.3: Ensemble-based hotspot selection (consensus voting)
- Q+.4: Adaptive rounds (skip Round 2 if entropy is low)
"""

import logging
from dataclasses import dataclass
from typing import Optional

import numpy as np
import torch
from torch import Tensor

logger = logging.getLogger(__name__)


@dataclass
class EpitopeRefinementConfig:
    """Configuration for epitope refinement pipeline."""

    contact_threshold: float = 8.0
    round1_samples: int = 10
    round1_noise_scale: float = 1.1
    entropy_threshold: float = 2.0
    hotspot_thresholds: list[float] = None
    min_consensus: int = 2

    def __post_init__(self):
        if self.hotspot_thresholds is None:
            self.hotspot_thresholds = [0.1, 0.2, 0.3]

    @classmethod
    def from_feats(cls, feats: dict[str, Tensor]) -> Optional["EpitopeRefinementConfig"]:
        """Extract config from feature dict produced by featurizer."""
        if "epitope_contact_threshold" not in feats:
            return None

        hotspot_thresholds = feats["epitope_hotspot_thresholds"]
        # Remove batch dimension if present
        if hotspot_thresholds.dim() > 1:
            hotspot_thresholds = hotspot_thresholds[0]
        if hotspot_thresholds.dim() > 0:
            hotspot_thresholds = hotspot_thresholds.tolist()
        else:
            hotspot_thresholds = [hotspot_thresholds.item()]

        # Handle batch dimension - squeeze any leading batch dims for scalar tensors
        def _scalar(t):
            return t.flatten()[0].item()

        return cls(
            contact_threshold=_scalar(feats["epitope_contact_threshold"]),
            round1_samples=int(_scalar(feats["epitope_round1_samples"])),
            round1_noise_scale=_scalar(feats["epitope_round1_noise_scale"]),
            entropy_threshold=_scalar(feats["epitope_entropy_threshold"]),
            hotspot_thresholds=hotspot_thresholds,
            min_consensus=int(_scalar(feats["epitope_min_consensus"])),
        )


def extract_contacts_from_coords(
    coords: Tensor,
    antigen_atom_idx: Tensor,
    cdr_atom_idx: Tensor,
    threshold: float = 8.0,
) -> Tensor:
    """Extract CDR-antigen contacts from predicted coordinates.

    Parameters
    ----------
    coords : Tensor
        Predicted atom coordinates [num_samples, num_atoms, 3]
    antigen_atom_idx : Tensor
        Indices of antigen CA atoms [N_antigen]
    cdr_atom_idx : Tensor
        Indices of CDR CA atoms [N_cdr]
    threshold : float
        Contact distance threshold in Angstroms

    Returns
    -------
    Tensor
        Per-antigen-residue contact counts [num_samples, N_antigen], normalized to [0,1]
    """
    if coords.dim() == 2:
        coords = coords.unsqueeze(0)

    num_samples = coords.shape[0]
    n_antigen = len(antigen_atom_idx)

    # Extract antigen and CDR CA coordinates
    ag_coords = coords[:, antigen_atom_idx, :]  # [S, N_ag, 3]
    cdr_coords = coords[:, cdr_atom_idx, :]  # [S, N_cdr, 3]

    # Compute pairwise distances
    # [S, N_ag, 1, 3] - [S, 1, N_cdr, 3] → [S, N_ag, N_cdr]
    dists = torch.cdist(ag_coords, cdr_coords)  # [S, N_ag, N_cdr]

    # Count contacts per antigen residue (any CDR atom within threshold)
    contacts = (dists < threshold).any(dim=-1).float()  # [S, N_ag]

    return contacts


def build_confidence_weighted_heatmap(
    contacts: Tensor,
    confidence_scores: Tensor,
) -> Tensor:
    """Build confidence-weighted contact heatmap (Q+.2).

    Parameters
    ----------
    contacts : Tensor
        Per-sample contact map [num_samples, N_antigen]
    confidence_scores : Tensor
        Confidence score per sample [num_samples]

    Returns
    -------
    Tensor
        Confidence-weighted heatmap [N_antigen]
    """
    # Weight each sample's contacts by its confidence score
    weights = confidence_scores.float()

    # Avoid division by zero
    weight_sum = weights.sum()
    if weight_sum < 1e-8:
        return contacts.mean(dim=0)

    # Weighted average
    weighted_contacts = (contacts * weights.unsqueeze(-1)).sum(dim=0) / weight_sum
    return weighted_contacts


def ensemble_hotspot_selection(
    heatmap: Tensor,
    thresholds: list[float],
    min_consensus: int = 2,
) -> Tensor:
    """Select hotspot residues using ensemble voting (Q+.3).

    A residue is a hotspot if it exceeds threshold in at least
    min_consensus threshold levels.

    Parameters
    ----------
    heatmap : Tensor
        Normalized contact heatmap [N_antigen]
    thresholds : list[float]
        Threshold values for voting
    min_consensus : int
        Minimum number of thresholds a residue must exceed

    Returns
    -------
    Tensor
        Boolean mask of hotspot residues [N_antigen]
    """
    votes = torch.zeros_like(heatmap)
    for thresh in thresholds:
        votes += (heatmap > thresh).float()

    hotspots = votes >= min_consensus
    return hotspots


def compute_contact_entropy(heatmap: Tensor) -> float:
    """Compute Shannon entropy of contact heatmap (Q+.4).

    Low entropy = concentrated signal (clear hotspot).
    High entropy = dispersed signal (unclear).

    Parameters
    ----------
    heatmap : Tensor
        Normalized contact heatmap [N_antigen]

    Returns
    -------
    float
        Shannon entropy value
    """
    # Normalize to probability distribution
    total = heatmap.sum()
    if total < 1e-8:
        return float('inf')  # No contacts → maximum uncertainty

    p = heatmap / total
    # Filter zero entries
    p = p[p > 0]
    entropy = -(p * torch.log(p)).sum().item()
    return entropy


def generate_pocket_constraints_from_hotspots(
    hotspot_mask: Tensor,
    antigen_res_idx: Tensor,
    antigen_chain_id: int,
    cdr_chain_ids: list[int],
    max_distance: float = 8.0,
) -> list[tuple[int, list[tuple[int, int]], float, bool]]:
    """Generate Boltz pocket constraints from hotspot residues.

    Parameters
    ----------
    hotspot_mask : Tensor
        Boolean mask of hotspot antigen residues [N_antigen]
    antigen_res_idx : Tensor
        0-indexed residue indices for antigen residues [N_antigen]
    antigen_chain_id : int
        Chain ID of the antigen
    cdr_chain_ids : list[int]
        Chain IDs of the CDR-containing chains (binders)
    max_distance : float
        Maximum distance for pocket constraint

    Returns
    -------
    list
        List of pocket constraints in Boltz format:
        [(binder_chain_id, contacts, max_distance, force), ...]
    """
    hotspot_indices = torch.where(hotspot_mask)[0]
    if len(hotspot_indices) == 0:
        return []

    # Get the residue indices of hotspot residues (0-indexed)
    hotspot_res_indices = antigen_res_idx[hotspot_indices].tolist()

    # Create contact list: [(chain_id, res_idx), ...]
    contacts = [(antigen_chain_id, int(res_idx)) for res_idx in hotspot_res_indices]

    # Create pocket constraint for each binder chain
    pocket_constraints = []
    for binder_chain_id in cdr_chain_ids:
        pocket_constraints.append((binder_chain_id, contacts, max_distance, True))

    return pocket_constraints


def analyze_round_results(
    all_coords: Tensor,
    all_confidence: Tensor,
    feats: dict[str, Tensor],
    config: EpitopeRefinementConfig,
) -> tuple[Tensor, Tensor, float, bool]:
    """Analyze results from a prediction round.

    Returns
    -------
    tuple of:
        - heatmap: confidence-weighted contact heatmap [N_antigen]
        - hotspot_mask: boolean mask of hotspot residues [N_antigen]
        - entropy: contact entropy value
        - should_skip_round2: whether to skip Round 2 (Q+.4)
    """
    antigen_atom_idx = feats["epitope_antigen_atom_index"]
    cdr_atom_idx = feats["epitope_cdr_atom_index"]

    # Remove batch dimension if present
    if antigen_atom_idx.dim() > 1:
        antigen_atom_idx = antigen_atom_idx[0]
    if cdr_atom_idx.dim() > 1:
        cdr_atom_idx = cdr_atom_idx[0]

    # Move indices to same device as coords
    device = all_coords.device
    antigen_atom_idx = antigen_atom_idx.to(device)
    cdr_atom_idx = cdr_atom_idx.to(device)

    # Extract contacts
    contacts = extract_contacts_from_coords(
        all_coords, antigen_atom_idx, cdr_atom_idx,
        threshold=config.contact_threshold,
    )

    # Build confidence-weighted heatmap (Q+.2)
    heatmap = build_confidence_weighted_heatmap(contacts, all_confidence)

    # Ensemble-based hotspot selection (Q+.3)
    hotspot_mask = ensemble_hotspot_selection(
        heatmap, config.hotspot_thresholds, config.min_consensus,
    )

    # Compute entropy for adaptive round decision (Q+.4)
    entropy = compute_contact_entropy(heatmap)
    should_skip = entropy < config.entropy_threshold

    n_hotspots = hotspot_mask.sum().item()
    logger.info(
        f"Epitope analysis: {n_hotspots} hotspot residues, "
        f"entropy={entropy:.3f}, skip_round2={should_skip}"
    )

    return heatmap, hotspot_mask, entropy, should_skip
