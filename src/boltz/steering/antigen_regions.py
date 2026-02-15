"""
Antigen region partitioning for epitope discovery.

Provides methods to partition antigen surface into regions for
region-specific beta-scaling during epitope scanning.
"""

from typing import List, Tuple
import numpy as np


class AntigenRegion:
    """Represents a region on the antigen surface."""

    def __init__(self, region_id: int, residue_indices: List[int], center: np.ndarray = None, radius: float = None):
        """
        Initialize an antigen region.

        Args:
            region_id: Unique identifier for this region
            residue_indices: List of residue indices (0-indexed) in this region
            center: 3D coordinates of region center (optional)
            radius: Radius of region in Angstroms (optional)
        """
        self.region_id = region_id
        self.residue_indices = residue_indices
        self.center = center
        self.radius = radius

    def __repr__(self):
        return f"AntigenRegion(id={self.region_id}, residues={len(self.residue_indices)}, center={self.center})"


def partition_antigen_spherical(
    num_residues: int,
    num_regions: int = 10,
    overlap_ratio: float = 0.2,
) -> List[AntigenRegion]:
    """
    Partition antigen into overlapping spherical regions.

    Args:
        num_residues: Total number of residues in antigen
        num_regions: Number of regions to create (default: 10)
        overlap_ratio: Fraction of residues to overlap between regions (default: 0.2)

    Returns:
        List of AntigenRegion objects
    """
    regions = []
    region_size = num_residues / num_regions
    overlap_size = int(region_size * overlap_ratio)

    for i in range(num_regions):
        start = max(0, int(i * region_size) - overlap_size // 2)
        end = min(num_residues, int((i + 1) * region_size) + overlap_size // 2)
        residue_indices = list(range(start, end))

        region = AntigenRegion(
            region_id=i,
            residue_indices=residue_indices,
        )
        regions.append(region)

    return regions


def create_region_mask(
    region: AntigenRegion,
    num_residues: int,
) -> np.ndarray:
    """
    Create a binary mask for a region.

    Args:
        region: AntigenRegion object
        num_residues: Total number of residues in antigen

    Returns:
        Binary mask array (1.0 for residues in region, 0.0 otherwise)
    """
    mask = np.zeros(num_residues, dtype=np.float32)
    mask[region.residue_indices] = 1.0
    return mask


def create_all_region_masks(
    regions: List[AntigenRegion],
    num_residues: int,
) -> np.ndarray:
    """
    Create masks for all regions.

    Args:
        regions: List of AntigenRegion objects
        num_residues: Total number of residues in antigen

    Returns:
        Array of masks with shape (num_regions, num_residues)
    """
    masks = np.zeros((len(regions), num_residues), dtype=np.float32)
    for i, region in enumerate(regions):
        masks[i] = create_region_mask(region, num_residues)
    return masks
