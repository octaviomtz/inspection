"""
Output processing for epitope scanning results.

Handles saving epitope heatmaps, hotspots, and statistics to JSON files.
"""

from pathlib import Path
import json
import numpy as np
from typing import Dict, List


def save_epitope_results(
    output_dir: Path,
    complex_name: str,
    heatmap: np.ndarray,
    hotspots: List[int],
    statistics: Dict,
    region_contact_maps: Dict = None,
):
    """
    Save epitope scanning results to JSON files.

    Args:
        output_dir: Directory to save results
        complex_name: Name of the complex
        heatmap: Epitope propensity heatmap (shape: num_residues)
        hotspots: List of hotspot residue indices
        statistics: Dictionary of scanning statistics
        region_contact_maps: Optional contact maps per region
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    # Save epitope heatmap
    heatmap_data = {
        "epitope_heatmap": heatmap.tolist() if isinstance(heatmap, np.ndarray) else list(heatmap),
        "epitope_hotspots": hotspots,
        "heatmap_statistics": {
            "min": float(np.min(heatmap)),
            "max": float(np.max(heatmap)),
            "mean": float(np.mean(heatmap)),
            "median": float(np.median(heatmap)),
            "std": float(np.std(heatmap)),
        },
    }

    heatmap_path = output_dir / f"{complex_name}_epitope_heatmap.json"
    with open(heatmap_path, "w") as f:
        json.dump(heatmap_data, f, indent=2)

    # Save hotspots
    hotspots_data = {
        "hotspot_residues": hotspots,
        "num_hotspots": len(hotspots),
        "hotspot_percentile": 80,  # Default percentile
    }

    hotspots_path = output_dir / f"{complex_name}_epitope_hotspots.json"
    with open(hotspots_path, "w") as f:
        json.dump(hotspots_data, f, indent=2)

    # Save statistics
    stats_path = output_dir / f"{complex_name}_epitope_statistics.json"
    with open(stats_path, "w") as f:
        json.dump(statistics, f, indent=2)

    # Optionally save per-region contact maps
    if region_contact_maps:
        contact_maps_data = {
            f"region_{region_id}": contact_map.tolist()
            if isinstance(contact_map, np.ndarray)
            else list(contact_map)
            for region_id, contact_map in region_contact_maps.items()
        }

        contact_maps_path = output_dir / f"{complex_name}_epitope_region_contacts.json"
        with open(contact_maps_path, "w") as f:
            json.dump(contact_maps_data, f, indent=2)

    return {
        "heatmap_file": str(heatmap_path),
        "hotspots_file": str(hotspots_path),
        "statistics_file": str(stats_path),
    }


def load_epitope_heatmap(heatmap_file: Path) -> np.ndarray:
    """Load epitope heatmap from JSON file."""
    with open(heatmap_file, "r") as f:
        data = json.load(f)
    return np.array(data["epitope_heatmap"], dtype=np.float32)


def load_epitope_hotspots(hotspots_file: Path) -> List[int]:
    """Load epitope hotspots from JSON file."""
    with open(hotspots_file, "r") as f:
        data = json.load(f)
    return data["hotspot_residues"]
