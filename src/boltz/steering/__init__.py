"""
Steering mechanisms for Boltz2 inference.

This module provides various steering mechanisms for guiding Boltz2 predictions:
- Region-specific β-scaling for epitope discovery (K+)
- CDR3 conformation steering
- Antigen orientation optimization
"""

from .antigen_regions import (
    AntigenRegion,
    partition_antigen_spherical,
    create_region_mask,
    create_all_region_masks,
)
from .beta_scaling_config import (
    RegionBetaConfig,
    EpitopeScanningConfig,
    parse_epitope_scanning_constraint,
)
from .contact_tracking import (
    compute_contacts,
    weight_contacts_by_confidence,
    ContactHeatmapAccumulator,
)
from .epitope_scanning import (
    EpitopeScanningParams,
    EpitopeScanningCoordinator,
    extract_antigen_coordinates,
    save_epitope_heatmap_json,
)
from .output_processing import (
    save_epitope_results,
    load_epitope_heatmap,
    load_epitope_hotspots,
)

__all__ = [
    "AntigenRegion",
    "partition_antigen_spherical",
    "create_region_mask",
    "create_all_region_masks",
    "RegionBetaConfig",
    "EpitopeScanningConfig",
    "parse_epitope_scanning_constraint",
    "compute_contacts",
    "weight_contacts_by_confidence",
    "ContactHeatmapAccumulator",
    "EpitopeScanningParams",
    "EpitopeScanningCoordinator",
    "extract_antigen_coordinates",
    "save_epitope_heatmap_json",
    "save_epitope_results",
    "load_epitope_heatmap",
    "load_epitope_hotspots",
]
