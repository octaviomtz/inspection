"""
Contact tracking and epitope heatmap accumulation.

Computes and accumulates contact information across multiple predictions
to generate epitope propensity maps.
"""

from typing import Dict, List, Optional, Tuple
import numpy as np
import torch


def compute_contacts(
    cdr_coords: np.ndarray,
    antigen_coords: np.ndarray,
    contact_threshold: float = 8.0,
    atoms_per_residue: int = 10,  # Average atoms per residue
) -> np.ndarray:
    """
    Compute which antigen residues are in contact with CDR residues.

    Args:
        cdr_coords: CDR atom coordinates, shape (num_cdr_atoms, 3)
        antigen_coords: Antigen atom coordinates, shape (num_antigen_atoms, 3)
        contact_threshold: Distance threshold for contact in Angstroms
        atoms_per_residue: Average atoms per residue (for atom->residue aggregation)

    Returns:
        Binary contact array, shape (num_antigen_residues,)
    """
    # Compute pairwise distances
    # Shape: (num_cdr_atoms, num_antigen_atoms)
    diff = cdr_coords[:, np.newaxis, :] - antigen_coords[np.newaxis, :, :]  # (N_cdr, N_ag, 3)
    distances = np.linalg.norm(diff, axis=2)  # (N_cdr, N_ag)

    # Find minimum distance for each antigen atom
    min_distances = distances.min(axis=0)  # (N_ag,)

    # Create binary contact array at atom level
    atom_contacts = (min_distances <= contact_threshold).astype(np.float32)

    # Aggregate to residue level
    num_atoms = len(atom_contacts)
    num_residues = (num_atoms + atoms_per_residue - 1) // atoms_per_residue  # Round up

    residue_contacts = np.zeros(num_residues, dtype=np.float32)
    for res_idx in range(num_residues):
        atom_start = res_idx * atoms_per_residue
        atom_end = min((res_idx + 1) * atoms_per_residue, num_atoms)
        # Residue is in contact if any of its atoms are in contact
        residue_contacts[res_idx] = np.any(atom_contacts[atom_start:atom_end])

    return residue_contacts


def weight_contacts_by_confidence(
    contacts: np.ndarray,
    confidence_scores: np.ndarray,
) -> np.ndarray:
    """
    Weight contact contributions by confidence scores.

    Args:
        contacts: Binary contact array, shape (num_residues,)
        confidence_scores: Confidence scores (pLDDT/pAE), shape (num_residues,)

    Returns:
        Weighted contact array, shape (num_residues,)
    """
    # Normalize confidence to [0, 1]
    conf_normalized = np.clip(confidence_scores / 100.0, 0.0, 1.0)

    # Weight contacts by confidence
    weighted_contacts = contacts * conf_normalized

    return weighted_contacts


class ContactHeatmapAccumulator:
    """Accumulates contacts across multiple predictions into an epitope heatmap."""

    def __init__(self, num_antigen_residues: int):
        """
        Initialize accumulator.

        Args:
            num_antigen_residues: Number of residues in antigen
        """
        self.num_residues = num_antigen_residues
        self.heatmap = np.zeros(num_antigen_residues, dtype=np.float32)
        self.contact_counts = np.zeros(num_antigen_residues, dtype=np.int32)
        self.num_predictions = 0

    def add_prediction(
        self,
        contacts: np.ndarray,
        confidence_scores: Optional[np.ndarray] = None,
    ):
        """
        Add contacts from a single prediction to the heatmap.

        Args:
            contacts: Binary contact array, shape (num_residues,)
            confidence_scores: Optional confidence scores, shape (num_residues,)
        """
        if contacts.shape[0] != self.num_residues:
            raise ValueError(
                f"Contact array shape {contacts.shape[0]} does not match "
                f"expected {self.num_residues}"
            )

        # Weight by confidence if provided
        if confidence_scores is not None:
            weighted = weight_contacts_by_confidence(contacts, confidence_scores)
        else:
            weighted = contacts

        self.heatmap += weighted
        self.contact_counts += (contacts > 0.5).astype(np.int32)
        self.num_predictions += 1

    def get_heatmap(self, normalize: bool = True) -> np.ndarray:
        """
        Get the accumulated epitope heatmap.

        Args:
            normalize: Whether to normalize by number of predictions

        Returns:
            Epitope propensity heatmap, shape (num_residues,)
        """
        if self.num_predictions == 0:
            return self.heatmap.copy()

        if normalize:
            return self.heatmap / max(1, self.num_predictions)
        else:
            return self.heatmap.copy()

    def get_contact_frequency(self) -> np.ndarray:
        """
        Get frequency of contacts (0-1 range).

        Returns:
            Contact frequency per residue, shape (num_residues,)
        """
        if self.num_predictions == 0:
            return self.contact_counts.astype(np.float32)

        return self.contact_counts.astype(np.float32) / self.num_predictions

    def reset(self):
        """Reset the accumulator."""
        self.heatmap.fill(0.0)
        self.contact_counts.fill(0)
        self.num_predictions = 0

    def get_epitope_hotspots(
        self,
        percentile: float = 80.0,
    ) -> List[int]:
        """
        Get indices of high-confidence epitope hotspots.

        Args:
            percentile: Percentile threshold for hotspot definition (default: 80)

        Returns:
            List of residue indices in hotspots
        """
        heatmap = self.get_heatmap(normalize=True)
        threshold = np.percentile(heatmap, percentile)
        hotspot_indices = np.where(heatmap >= threshold)[0]
        return hotspot_indices.tolist()

    def get_statistics(self) -> dict:
        """
        Get scanning statistics.

        Returns:
            Dictionary of statistics
        """
        heatmap = self.get_heatmap(normalize=True)
        hotspots = self.get_epitope_hotspots(percentile=80)
        freq = self.get_contact_frequency()

        return {
            "num_predictions": self.num_predictions,
            "num_residues": self.num_residues,
            "num_hotspots": len(hotspots),
            "heatmap_statistics": {
                "min": float(np.min(heatmap)) if len(heatmap) > 0 else 0.0,
                "max": float(np.max(heatmap)) if len(heatmap) > 0 else 0.0,
                "mean": float(np.mean(heatmap)),
                "median": float(np.median(heatmap)),
                "std": float(np.std(heatmap)),
            },
            "contact_frequency_statistics": {
                "min": float(np.min(freq)) if len(freq) > 0 else 0.0,
                "max": float(np.max(freq)) if len(freq) > 0 else 0.0,
                "mean": float(np.mean(freq)),
                "median": float(np.median(freq)),
            },
        }
