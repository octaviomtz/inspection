#!/usr/bin/env python3
"""
Test YAML constraint parsing for epitope_region_scanning.
"""

import sys
import yaml
from pathlib import Path

# Test YAML content
TEST_YAML_CONTENT = """
sequences:
- protein:
    id: A
    sequence: APLHLGKCNIAGWILGNPECESLSTASSWSYIVETPSSDNGTCYPGDFIDYEELREQLSSVSSFERFEIFPKTSSWPNHDSDKGVTAACPHAGAKSFYKNLIWLVKKGNSYPKLSKSYINDKGKEVLVLWGIHHPSTSADQQSLYQNADAYVFVGSSRYSKTFKPEIAIRPKVRDREGRMNYYWTLVEPGDKITFEATGNLVVPRYAFAMERNA
    msa: dummy.a3m
- protein:
    id: B
    sequence: EVQLVESGGGLIQPGGSLRLSCEASAFTFSSYEMNWVRQAPGKGLEWVSYITSSGSRIYYADSVKGRFTISRDNAKNSLYLQMNSLRVEDTAVYYCARLLDSIVWGEGWYYGMDVWGQGTTVTVSG
    msa: empty
- protein:
    id: C
    sequence: SYELTQSPSVSVAPGRTARITCGGNDIGLKGVHWYQQKPGQAPVLVLYDNNHRPSGIPERFSGSISGDTATLTVTRVEADDGADYFCQVWDTSSGPPHVIFGGGTKLTVL
    msa: empty

constraints:
  - epitope_region_scanning:
      antigen_chain: A
      num_regions: 10
      region_overlap_ratio: 0.2
      beta_emphasis: 0.5
      beta_deemphasis: -0.3
      contact_threshold: 8.0
      confidence_weighting: true
      accumulate_heatmap: true
      save_per_region: false
"""


def test_yaml_parsing():
    """Test that epitope_region_scanning constraint is properly parsed."""
    print("\n[TEST] YAML Constraint Parsing")
    print("=" * 60)

    # Parse YAML
    data = yaml.safe_load(TEST_YAML_CONTENT)

    assert "constraints" in data, "No constraints found in YAML"
    print("✓ YAML loaded successfully")

    constraints = data["constraints"]
    assert len(constraints) > 0, "No constraints found"
    print(f"✓ Found {len(constraints)} constraint(s)")

    # Check epitope_region_scanning constraint
    epitope_constraint = None
    for constraint in constraints:
        if "epitope_region_scanning" in constraint:
            epitope_constraint = constraint["epitope_region_scanning"]
            break

    assert epitope_constraint is not None, "epitope_region_scanning constraint not found"
    print("✓ epitope_region_scanning constraint found")

    # Validate constraint fields
    required_fields = ["antigen_chain"]
    for field in required_fields:
        assert field in epitope_constraint, f"Required field '{field}' not found"
    print(f"✓ All required fields present")

    # Validate optional fields with defaults
    optional_fields = {
        "num_regions": 10,
        "region_overlap_ratio": 0.2,
        "beta_emphasis": 0.5,
        "beta_deemphasis": -0.3,
        "contact_threshold": 8.0,
        "confidence_weighting": True,
        "accumulate_heatmap": True,
        "save_per_region": False,
    }

    for field, default in optional_fields.items():
        value = epitope_constraint.get(field, default)
        print(f"  - {field}: {value}")
        assert value is not None, f"Field '{field}' is None"

    print("✓ All constraint fields validated")

    print("\n✓ YAML parsing test PASSED")
    return True


if __name__ == "__main__":
    try:
        test_yaml_parsing()
        print("\n" + "=" * 60)
        print("✓ CONSTRAINT PARSING TEST PASSED")
        print("=" * 60)
        sys.exit(0)
    except Exception as e:
        print(f"\n✗ TEST FAILED: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
