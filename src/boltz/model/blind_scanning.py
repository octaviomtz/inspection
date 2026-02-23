"""Blind scanning utilities for epitope mapping.

Provides functions to partition an antigen surface into overlapping regions,
build pair masks for region-specific β-scaling, compute contacts from
predicted coordinates, and aggregate contacts into an epitope propensity
heatmap.
"""

import torch
from torch import Tensor


def partition_antigen_sequence(
    antigen_token_indices: Tensor,
    num_regions: int,
    overlap_fraction: float = 0.5,
) -> list[Tensor]:
    """Partition antigen token indices into overlapping sequential windows.

    Parameters
    ----------
    antigen_token_indices : Tensor
        [N_antigen] long tensor of token indices for antigen residues.
    num_regions : int
        Number of regions to partition into.
    overlap_fraction : float
        Fraction of window size to overlap between consecutive windows.

    Returns
    -------
    list[Tensor]
        List of index tensors, one per region.
    """
    n = len(antigen_token_indices)
    if n == 0 or num_regions <= 0:
        return []

    # Clamp num_regions to at most n
    num_regions = min(num_regions, n)

    if num_regions == 1:
        return [antigen_token_indices]

    # Compute window size and stride
    # stride * (num_regions - 1) + window_size = n  (approximately)
    # window_size = stride / (1 - overlap_fraction)
    stride = max(1, n // num_regions)
    window_size = max(1, int(stride / (1.0 - overlap_fraction)))
    window_size = min(window_size, n)  # Can't exceed total

    regions = []
    for i in range(num_regions):
        start = i * stride
        end = min(start + window_size, n)
        if start >= n:
            start = n - 1
            end = n
        regions.append(antigen_token_indices[start:end])

    return regions


def build_region_pair_masks(
    cdr_mask: Tensor,
    region_indices: Tensor,
    antigen_mask: Tensor,
    n_tokens: int,
) -> tuple[Tensor, Tensor]:
    """Build emphasis and de-emphasis pair masks for a scan region.

    The emphasis mask covers CDR-region pairs (interactions to emphasize).
    The de-emphasis mask covers CDR-antigen pairs outside the region.

    Parameters
    ----------
    cdr_mask : Tensor
        [N_tokens] boolean mask for CDR residues.
    region_indices : Tensor
        [N_region] long tensor of token indices for the current region.
    antigen_mask : Tensor
        [N_tokens] boolean mask for all antigen residues.
    n_tokens : int
        Total number of tokens.

    Returns
    -------
    tuple[Tensor, Tensor]
        (emph_mask, deemph_mask) as [N_tokens, N_tokens] bool tensors.
    """
    # Build region mask
    region_mask = torch.zeros(n_tokens, dtype=torch.bool, device=cdr_mask.device)
    region_mask[region_indices] = True

    # Emphasis: CDR x Region + Region x CDR (symmetric)
    emph_mask = (cdr_mask.unsqueeze(1) & region_mask.unsqueeze(0)) | \
                (region_mask.unsqueeze(1) & cdr_mask.unsqueeze(0))

    # De-emphasis: CDR x (Antigen \ Region) + (Antigen \ Region) x CDR
    non_region_antigen = antigen_mask & ~region_mask
    deemph_mask = (cdr_mask.unsqueeze(1) & non_region_antigen.unsqueeze(0)) | \
                  (non_region_antigen.unsqueeze(1) & cdr_mask.unsqueeze(0))

    return emph_mask, deemph_mask


def compute_contacts_from_coords(
    coords: Tensor,
    cdr_mask: Tensor,
    antigen_mask: Tensor,
    threshold: float,
) -> Tensor:
    """Compute binary contacts between CDR and antigen from predicted coords.

    Uses C-alpha (token center) coordinates. A contact exists when any
    CDR residue is within threshold distance of an antigen residue.

    Parameters
    ----------
    coords : Tensor
        [N_atoms, 3] predicted atom coordinates (use token centers).
    cdr_mask : Tensor
        [N_tokens] boolean mask for CDR residues.
    antigen_mask : Tensor
        [N_tokens] boolean mask for antigen residues.
    threshold : float
        Distance threshold in Angstroms.

    Returns
    -------
    Tensor
        [N_antigen] binary contact tensor (1 if any CDR residue within threshold).
    """
    cdr_indices = torch.where(cdr_mask)[0]
    antigen_indices = torch.where(antigen_mask)[0]

    if len(cdr_indices) == 0 or len(antigen_indices) == 0:
        return torch.zeros(antigen_mask.sum().item(), dtype=torch.float32, device=coords.device)

    # Get coordinates for CDR and antigen tokens
    # coords may be [N_atoms, 3]; we need token-level coordinates
    # Use the first atom per token (token center / C-alpha)
    cdr_coords = coords[cdr_indices]  # [N_cdr, 3]
    antigen_coords = coords[antigen_indices]  # [N_antigen, 3]

    # Compute pairwise distances: [N_cdr, N_antigen]
    dists = torch.cdist(cdr_coords.unsqueeze(0).float(), antigen_coords.unsqueeze(0).float()).squeeze(0)

    # Contact if minimum distance to any CDR residue < threshold
    min_dists = dists.min(dim=0).values  # [N_antigen]
    contacts = (min_dists < threshold).float()

    return contacts


def aggregate_contacts(
    all_contacts: list[Tensor],
    n_antigen: int,
) -> Tensor:
    """Aggregate contacts from all region scans into epitope propensity.

    Parameters
    ----------
    all_contacts : list[Tensor]
        List of [N_antigen] contact tensors from each region scan.
    n_antigen : int
        Number of antigen residues.

    Returns
    -------
    Tensor
        [N_antigen] float propensity in [0, 1].
    """
    if not all_contacts:
        return torch.zeros(n_antigen, dtype=torch.float32)

    # Stack and compute fraction of scans with contact
    stacked = torch.stack(all_contacts, dim=0)  # [N_regions, N_antigen]
    propensity = stacked.mean(dim=0)  # [N_antigen]

    return propensity
