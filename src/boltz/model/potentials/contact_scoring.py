"""
Contact scoring utilities for steering during diffusion.

Computes CDR-antigen contact scores to monitor steering progress.
"""

from typing import Optional, Tuple

import torch
import torch.nn.functional as F


def compute_cdr_antigen_contacts(
    atom_coords: torch.Tensor,
    feats: dict,
    cdr_indices: Optional[torch.Tensor] = None,
    antigen_indices: Optional[torch.Tensor] = None,
    contact_distance_threshold: float = 4.5,
) -> Tuple[float, int]:
    """
    Compute contact score between CDR and antigen regions.

    Parameters
    ----------
    atom_coords : torch.Tensor
        Atom coordinates [batch, num_atoms, 3].
    feats : dict
        Feature dictionary containing region masks.
    cdr_indices : torch.Tensor, optional
        Indices of CDR atoms. If None, extracted from feats.
    antigen_indices : torch.Tensor, optional
        Indices of antigen atoms. If None, extracted from feats.
    contact_distance_threshold : float
        Distance threshold for counting a contact (Å).

    Returns
    -------
    contact_score : float
        Fraction of CDR atoms within contact distance of antigen.
    num_contacts : int
        Total number of atom pairs in contact.
    """
    if atom_coords.dim() == 3:
        # Multiple samples: [batch, num_atoms, 3]
        atom_coords = atom_coords[0]  # Use first sample

    # Extract CDR and antigen indices from features
    if cdr_indices is None and "cdr_indices" in feats:
        cdr_indices = feats["cdr_indices"]

    if antigen_indices is None and "antigen_indices" in feats:
        antigen_indices = feats["antigen_indices"]

    if cdr_indices is None or antigen_indices is None:
        # Fallback: try to identify from chain information
        # This is a simplified version; actual implementation should extract from feats
        return 0.0, 0

    # Compute pairwise distances between CDR and antigen atoms
    cdr_coords = atom_coords[cdr_indices]  # [n_cdr, 3]
    antigen_coords = atom_coords[antigen_indices]  # [n_antigen, 3]

    # Compute pairwise distances [n_cdr, n_antigen]
    distances = torch.cdist(cdr_coords, antigen_coords, p=2)

    # Count contacts (distances below threshold)
    contacts = (distances < contact_distance_threshold).float()
    num_contacts = contacts.sum().item()
    num_cdr_atoms = cdr_coords.shape[0]

    # Contact score: fraction of CDR atoms in contact
    contact_score = num_contacts / (num_cdr_atoms * antigen_coords.shape[0]) if num_cdr_atoms > 0 else 0.0

    return float(contact_score), int(num_contacts)


def extract_cdr_and_antigen_indices(feats: dict) -> Tuple[Optional[torch.Tensor], Optional[torch.Tensor]]:
    """
    Extract CDR and antigen atom indices from features.

    Parameters
    ----------
    feats : dict
        Feature dictionary.

    Returns
    -------
    cdr_indices : torch.Tensor or None
        Indices of CDR atoms.
    antigen_indices : torch.Tensor or None
        Indices of antigen atoms.
    """
    cdr_indices = None
    antigen_indices = None

    # Try to find CDR3 regions from constraint information
    if "inference_cdr3_constraints" in feats:
        constraints = feats["inference_cdr3_constraints"]
        if constraints is not None and len(constraints) > 0:
            # Extract atom indices from CDR3 constraints
            cdr_indices = torch.cat([c["atom_indices"] for c in constraints])

    # Try to find antigen indices from chain information
    if "mol_type" in feats:
        # mol_type: 0 = protein, 1 = RNA, 2 = ligand, etc.
        # Antigen is typically a different chain
        mol_type = feats["mol_type"]
        if mol_type is not None:
            # Find atoms that belong to antigen chain
            # This is a simplified version; actual implementation depends on chain info
            antigen_mask = mol_type > 0  # Non-protein molecules
            antigen_indices = torch.where(antigen_mask)[0]

    return cdr_indices, antigen_indices


def compute_contact_score_from_pae(
    pae: torch.Tensor,
    cdr_token_indices: Optional[torch.Tensor] = None,
    antigen_token_indices: Optional[torch.Tensor] = None,
    pae_threshold: float = 5.0,
) -> float:
    """
    Compute contact score from predicted aligned error (PAE).

    More robust than distance-based contacts when atomic coordinates are uncertain.

    Parameters
    ----------
    pae : torch.Tensor
        PAE matrix [num_tokens, num_tokens].
    cdr_token_indices : torch.Tensor, optional
        Indices of CDR tokens.
    antigen_token_indices : torch.Tensor, optional
        Indices of antigen tokens.
    pae_threshold : float
        PAE threshold for considering atoms as in contact.

    Returns
    -------
    float
        Contact score based on PAE values.
    """
    if pae is None or pae.shape[0] == 0:
        return 0.0

    if cdr_token_indices is None or antigen_token_indices is None:
        return 0.0

    # Extract PAE submatrix for CDR-antigen interface
    pae_interface = pae[
        cdr_token_indices[:, None],
        antigen_token_indices[None, :],
    ]

    # Count high-confidence contacts (low PAE)
    contacts = (pae_interface < pae_threshold).float()
    contact_score = contacts.mean().item()

    return contact_score
