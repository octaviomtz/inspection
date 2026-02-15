"""
Configuration for region-specific beta-scaling during epitope scanning.

Handles the configuration of beta values for different antigen regions
to guide the diffusion process toward epitope discovery.
"""

from typing import Dict, List, Optional
from dataclasses import dataclass


@dataclass
class RegionBetaConfig:
    """Configuration for beta-scaling on a single antigen region."""

    region_id: int
    beta_emphasis: float = 0.5  # Beta value to emphasize this region
    beta_deemphasis: float = -0.3  # Beta value for other regions


@dataclass
class EpitopeScanningConfig:
    """Configuration for epitope region scanning with beta-scaling."""

    # Region definition
    antigen_chain: str  # e.g., 'A'
    num_regions: int = 10  # Number of regions to scan
    region_overlap_ratio: float = 0.2  # Overlap between regions

    # Beta-scaling parameters
    beta_emphasis: float = 0.5  # Beta value when emphasizing a region
    beta_deemphasis: float = -0.3  # Beta value for de-emphasized regions

    # Contact scoring and accumulation
    contact_threshold: float = 8.0  # Distance threshold for contacts (Angstroms)
    confidence_weighting: bool = True  # Weight contacts by pLDDT/pAE

    # Output
    accumulate_heatmap: bool = True  # Accumulate contact heatmap across predictions
    save_per_region: bool = False  # Save predictions for each region separately

    def get_region_beta_config(self, emphasize_region_id: int) -> Dict[int, float]:
        """
        Get beta configuration for a specific region emphasis.

        Args:
            emphasize_region_id: ID of region to emphasize

        Returns:
            Dict mapping region_id -> beta_value
        """
        config = {}
        for i in range(self.num_regions):
            if i == emphasize_region_id:
                config[i] = self.beta_emphasis
            else:
                config[i] = self.beta_deemphasis
        return config


def parse_epitope_scanning_constraint(constraint_dict: dict) -> Optional[EpitopeScanningConfig]:
    """
    Parse epitope_region_scanning constraint from YAML.

    Args:
        constraint_dict: Constraint dictionary from YAML

    Returns:
        EpitopeScanningConfig object or None if not a valid epitope scanning constraint
    """
    if "epitope_region_scanning" not in constraint_dict:
        return None

    scanning_config = constraint_dict["epitope_region_scanning"]

    config = EpitopeScanningConfig(
        antigen_chain=scanning_config.get("antigen_chain"),
        num_regions=scanning_config.get("num_regions", 10),
        region_overlap_ratio=scanning_config.get("region_overlap_ratio", 0.2),
        beta_emphasis=scanning_config.get("beta_emphasis", 0.5),
        beta_deemphasis=scanning_config.get("beta_deemphasis", -0.3),
        contact_threshold=scanning_config.get("contact_threshold", 8.0),
        confidence_weighting=scanning_config.get("confidence_weighting", True),
        accumulate_heatmap=scanning_config.get("accumulate_heatmap", True),
        save_per_region=scanning_config.get("save_per_region", False),
    )

    return config
