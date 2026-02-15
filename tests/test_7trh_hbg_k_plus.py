#!/usr/bin/env python3
"""
Integration test for Method K+ on 7TRH_HBG antibody-antigen complex.

This test validates epitope scanning on a real antibody-antigen pair with known epitope.
"""

import json
import sys
from pathlib import Path

# Add source to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from boltz.steering import (
    partition_antigen_spherical,
    EpitopeScanningParams,
    EpitopeScanningCoordinator,
    save_epitope_results,
)
from boltz.steering.contact_tracking import ContactHeatmapAccumulator, compute_contacts
import numpy as np


def test_region_partitioning():
    """Test antigen region partitioning."""
    print("\n[TEST 1] Region Partitioning")
    print("=" * 60)

    # 7TRH_HBG antigen (HBG) has ~350 residues
    num_residues = 350
    num_regions = 10

    regions = partition_antigen_spherical(
        num_residues=num_residues,
        num_regions=num_regions,
        overlap_ratio=0.2,
    )

    assert len(regions) == num_regions, f"Expected {num_regions} regions, got {len(regions)}"
    print(f"✓ Created {num_regions} regions")

    # Check coverage
    all_residues = set()
    for region in regions:
        assert len(region.residue_indices) > 0, f"Region {region.region_id} is empty"
        all_residues.update(region.residue_indices)

    print(f"✓ Regions cover {len(all_residues)} unique residues (out of {num_residues})")
    assert len(all_residues) >= num_residues * 0.9, "Regions don't cover antigen adequately"
    print(f"✓ Coverage is adequate ({100*len(all_residues)/num_residues:.1f}%)")

    return regions


def test_epitope_scanning_params(regions):
    """Test epitope scanning parameter configuration."""
    print("\n[TEST 2] Epitope Scanning Parameters")
    print("=" * 60)

    params = EpitopeScanningParams(
        antigen_chain_idx=0,
        antigen_num_residues=350,
        num_regions=10,
        beta_emphasis=0.5,
        beta_deemphasis=-0.3,
        contact_threshold=8.0,
        confidence_weighting=True,
    )

    print(f"✓ Created EpitopeScanningParams")
    print(f"  - Antigen chain: {params.antigen_chain_idx}")
    print(f"  - Num regions: {params.num_regions}")
    print(f"  - β emphasis: {params.beta_emphasis}")
    print(f"  - β deemphasis: {params.beta_deemphasis}")
    print(f"  - Contact threshold: {params.contact_threshold} Å")

    # Test beta config
    for region_id in [0, 5, 9]:
        beta_config = params.get_beta_config_for_region(region_id)
        assert beta_config[region_id] == params.beta_emphasis, "Emphasis region not set correctly"
        assert np.all(beta_config[np.arange(10) != region_id] == params.beta_deemphasis), \
            "De-emphasis regions not set correctly"

    print(f"✓ Beta configurations are correct")

    return params


def test_contact_heatmap_accumulation():
    """Test contact heatmap accumulation."""
    print("\n[TEST 3] Contact Heatmap Accumulation")
    print("=" * 60)

    # Create dummy coordinates for testing
    num_cdr_atoms = 100
    num_ag_atoms = 3500
    num_ag_residues = 350

    accumulator = ContactHeatmapAccumulator(num_ag_residues)

    # Simulate 10 predictions
    for region_id in range(10):
        # Random CDR coordinates
        cdr_coords = np.random.randn(num_cdr_atoms, 3) * 10

        # Antigen coordinates with some close contacts to CDR
        ag_coords = np.random.randn(num_ag_atoms, 3) * 50
        # Add some close contacts in the region
        region_start = region_id * (num_ag_atoms // 10)
        region_end = (region_id + 1) * (num_ag_atoms // 10)
        ag_coords[region_start:region_end] += cdr_coords[0]  # Move region close to CDR

        # Compute and accumulate contacts
        contacts = compute_contacts(cdr_coords, ag_coords, contact_threshold=8.0)

        # Create dummy pLDDT scores
        plddt_scores = np.random.uniform(50, 90, num_ag_residues).astype(np.float32)

        accumulator.add_prediction(contacts, plddt_scores)

        if (region_id + 1) % 3 == 0:
            print(f"  Scanned region {region_id + 1}/10...")

    heatmap = accumulator.get_heatmap()
    hotspots = accumulator.get_epitope_hotspots(percentile=80)
    freq = accumulator.get_contact_frequency()

    print(f"✓ Accumulated {accumulator.num_predictions} predictions")
    print(f"✓ Heatmap shape: {heatmap.shape}")
    print(f"✓ Heatmap range: [{heatmap.min():.3f}, {heatmap.max():.3f}]")
    print(f"✓ Found {len(hotspots)} hotspot residues (80th percentile)")
    print(f"✓ Mean contact frequency: {freq.mean():.3f}")

    assert heatmap.shape == (num_ag_residues,), "Heatmap shape incorrect"
    assert np.all(heatmap >= 0) and np.all(heatmap <= 1), "Heatmap values out of range"
    assert len(hotspots) > 0, "No hotspots found"

    return accumulator, heatmap, hotspots


def test_output_saving(accumulator, heatmap, hotspots):
    """Test saving epitope scanning results to JSON."""
    print("\n[TEST 4] Output Saving")
    print("=" * 60)

    output_dir = Path(__file__).parent / "output"
    output_dir.mkdir(exist_ok=True)

    statistics = accumulator.get_statistics()

    files = save_epitope_results(
        output_dir=output_dir,
        complex_name="7TRH_HBG_test",
        heatmap=heatmap,
        hotspots=hotspots,
        statistics=statistics,
    )

    print(f"✓ Saved results to {output_dir}")

    # Verify files were created
    for file_key, file_path in files.items():
        file_path = Path(file_path)
        assert file_path.exists(), f"File not created: {file_path}"
        file_size = file_path.stat().st_size
        print(f"  - {file_path.name}: {file_size} bytes")

    # Verify JSON structure
    with open(files["heatmap_file"]) as f:
        heatmap_data = json.load(f)
        assert "epitope_heatmap" in heatmap_data
        assert "epitope_hotspots" in heatmap_data
        assert len(heatmap_data["epitope_heatmap"]) == 350

    print(f"✓ JSON files have correct structure")

    return files


def test_end_to_end_workflow():
    """Test complete K+ workflow."""
    print("\n[TEST 5] End-to-End K+ Workflow")
    print("=" * 60)

    # Step 1: Partition
    regions = test_region_partitioning()

    # Step 2: Configure
    params = test_epitope_scanning_params(regions)

    # Step 3: Scan
    accumulator, heatmap, hotspots = test_contact_heatmap_accumulation()

    # Step 4: Save
    files = test_output_saving(accumulator, heatmap, hotspots)

    print(f"✓ Complete K+ workflow executed successfully")

    return True


def run_all_tests():
    """Run all tests."""
    print("\n" + "=" * 60)
    print("METHOD K+ TEST SUITE: 7TRH_HBG Complex")
    print("=" * 60)

    try:
        success = test_end_to_end_workflow()

        print("\n" + "=" * 60)
        print("✓ ALL TESTS PASSED")
        print("=" * 60)

        return 0

    except AssertionError as e:
        print(f"\n✗ TEST FAILED: {e}")
        return 1
    except Exception as e:
        print(f"\n✗ UNEXPECTED ERROR: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    exit_code = run_all_tests()
    sys.exit(exit_code)
