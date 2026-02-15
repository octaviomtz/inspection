"""
Epitope region scanning orchestration for K+ method.

Orchestrates running multiple predictions with region-specific β-scaling
to discover antibody epitopes on unknown antigens.
"""

from typing import Dict, List, Optional, Tuple
import numpy as np
import torch
from pathlib import Path

from boltz.steering.antigen_regions import (
    partition_antigen_spherical,
    create_all_region_masks,
)
from boltz.steering.contact_tracking import ContactHeatmapAccumulator


class EpitopeScanningParams:
    """Parameters for epitope region scanning."""

    def __init__(
        self,
        antigen_chain_idx: int,
        antigen_num_residues: int,
        num_regions: int = 10,
        beta_emphasis: float = 0.5,
        beta_deemphasis: float = -0.3,
        contact_threshold: float = 8.0,
        confidence_weighting: bool = True,
    ):
        """
        Initialize epitope scanning parameters.

        Args:
            antigen_chain_idx: Index of antigen chain in structure
            antigen_num_residues: Number of residues in antigen
            num_regions: Number of regions to scan (default: 10)
            beta_emphasis: β value for emphasized region (default: 0.5)
            beta_deemphasis: β value for other regions (default: -0.3)
            contact_threshold: Contact distance threshold in Å (default: 8.0)
            confidence_weighting: Weight by pLDDT/pAE (default: true)
        """
        self.antigen_chain_idx = antigen_chain_idx
        self.antigen_num_residues = antigen_num_residues
        self.num_regions = num_regions
        self.beta_emphasis = beta_emphasis
        self.beta_deemphasis = beta_deemphasis
        self.contact_threshold = contact_threshold
        self.confidence_weighting = confidence_weighting

        # Partition antigen surface into regions
        self.regions = partition_antigen_spherical(
            num_residues=antigen_num_residues,
            num_regions=num_regions,
        )

        # Create region masks
        self.region_masks = create_all_region_masks(
            self.regions,
            antigen_num_residues,
        )

    def get_beta_config_for_region(self, region_id: int) -> np.ndarray:
        """
        Get beta configuration for emphasizing a specific region.

        Args:
            region_id: ID of region to emphasize

        Returns:
            Beta values for all regions, shape (num_regions,)
        """
        beta_config = np.full(self.num_regions, self.beta_deemphasis, dtype=np.float32)
        beta_config[region_id] = self.beta_emphasis
        return beta_config

    def get_region_indices(self, region_id: int) -> List[int]:
        """Get residue indices for a specific region."""
        return self.regions[region_id].residue_indices


class EpitopeScanningCoordinator:
    """Coordinates epitope region scanning across multiple predictions."""

    def __init__(self, params: EpitopeScanningParams):
        """
        Initialize coordinator.

        Args:
            params: EpitopeScanningParams object
        """
        self.params = params
        self.accumulator = ContactHeatmapAccumulator(
            num_antigen_residues=params.antigen_num_residues
        )
        self.region_predictions = {}  # Store predictions for each region
        self.region_contact_maps = {}  # Store contact maps for each region

    def add_prediction(
        self,
        region_id: int,
        cdr_coords: np.ndarray,
        antigen_coords: np.ndarray,
        plddt_scores: Optional[np.ndarray] = None,
    ):
        """
        Add prediction from a region scan.

        Args:
            region_id: Which region was emphasized
            cdr_coords: CDR atom coordinates, shape (num_atoms, 3)
            antigen_coords: Antigen atom coordinates, shape (num_atoms, 3)
            plddt_scores: pLDDT scores for antigen residues (optional)
        """
        from boltz.steering.contact_tracking import compute_contacts, weight_contacts_by_confidence

        # Compute contacts for this prediction
        contacts = compute_contacts(
            cdr_coords,
            antigen_coords,
            contact_threshold=self.params.contact_threshold,
        )

        # Weight by confidence if available
        if self.params.confidence_weighting and plddt_scores is not None:
            weighted_contacts = weight_contacts_by_confidence(contacts, plddt_scores)
        else:
            weighted_contacts = contacts

        # Store for later analysis
        self.region_contact_maps[region_id] = weighted_contacts

        # Add to accumulator
        self.accumulator.add_prediction(
            contacts=weighted_contacts,
            confidence_scores=plddt_scores if self.params.confidence_weighting else None,
        )

    def get_epitope_heatmap(self) -> np.ndarray:
        """Get the accumulated epitope propensity heatmap."""
        return self.accumulator.get_heatmap(normalize=True)

    def get_epitope_hotspots(self, percentile: float = 80.0) -> List[int]:
        """Get epitope hotspot residues by percentile."""
        return self.accumulator.get_epitope_hotspots(percentile=percentile)

    def get_contact_frequency(self) -> np.ndarray:
        """Get per-residue contact frequency."""
        return self.accumulator.get_contact_frequency()

    def get_statistics(self) -> Dict:
        """Get scanning statistics."""
        heatmap = self.get_epitope_heatmap()
        hotspots = self.get_epitope_hotspots(percentile=80)
        freq = self.get_contact_frequency()

        return {
            "num_regions_scanned": self.params.num_regions,
            "total_predictions": self.accumulator.num_predictions,
            "antigen_residues": self.params.antigen_num_residues,
            "num_hotspots": len(hotspots),
            "mean_contact_propensity": float(np.mean(heatmap)),
            "max_contact_propensity": float(np.max(heatmap)),
            "min_contact_propensity": float(np.min(heatmap)),
            "mean_contact_frequency": float(np.mean(freq)),
        }

    def reset(self):
        """Reset the coordinator for reuse."""
        self.accumulator.reset()
        self.region_predictions.clear()
        self.region_contact_maps.clear()


def extract_antigen_coordinates(
    structure_dict: dict,
    antigen_chain_idx: int,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Extract antigen coordinates from structure dictionary.

    Args:
        structure_dict: Structure data dictionary from prediction
        antigen_chain_idx: Index of antigen chain

    Returns:
        Tuple of (antigen_coords, antigen_residue_plddt)
    """
    # Assuming structure_dict contains coordinates and pLDDT scores
    # This will be adapted based on actual Boltz2 output structure
    coords = structure_dict.get("coordinates", None)
    plddt = structure_dict.get("plddt", None)

    if coords is None:
        raise ValueError("No coordinates found in structure dictionary")

    # Extract antigen chain coordinates
    # (structure varies by Boltz2 version, may need adaptation)
    if isinstance(coords, dict):
        antigen_coords = coords.get(antigen_chain_idx, None)
    else:
        # Assumes coords is array-like with shape (num_chains, num_atoms, 3)
        antigen_coords = coords[antigen_chain_idx]

    # Extract pLDDT if available
    antigen_plddt = None
    if plddt is not None:
        if isinstance(plddt, dict):
            antigen_plddt = plddt.get(antigen_chain_idx, None)
        else:
            antigen_plddt = plddt[antigen_chain_idx]

    return np.array(antigen_coords, dtype=np.float32), antigen_plddt


def save_epitope_heatmap_json(
    heatmap: np.ndarray,
    hotspots: List[int],
    statistics: Dict,
    output_path: Path,
):
    """
    Save epitope heatmap and hotspots to JSON.

    Args:
        heatmap: Epitope propensity heatmap
        hotspots: Hotspot residue indices
        statistics: Scanning statistics
        output_path: Path to save JSON file
    """
    import json

    data = {
        "epitope_heatmap": heatmap.tolist(),
        "epitope_hotspots": hotspots,
        "statistics": statistics,
    }

    with open(output_path, "w") as f:
        json.dump(data, f, indent=2)
