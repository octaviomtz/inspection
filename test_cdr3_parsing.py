#!/usr/bin/env python3
"""Test CDR3 beta scaling YAML parsing."""

import sys
from pathlib import Path

# Add src to path
repo_root = Path(__file__).parent
sys.path.insert(0, str(repo_root / "src"))

try:
    from boltz.data.parse.yaml import parse_yaml

    yaml_file = repo_root / "yaml_cdr3_beta_scaling/7TRH_HBG_cdr3_beta.yml"

    print(f"Testing YAML parsing: {yaml_file}")
    print(f"File exists: {yaml_file.exists()}")

    if yaml_file.exists():
        # Try parsing the YAML
        try:
            result = parse_yaml(str(yaml_file))
            print("\n✓ YAML parsing successful!")
            print(f"  - Number of targets: {len(result)}")

            if result:
                target = result[0]
                print(f"  - Target ID: {target.record.id}")
                print(f"  - Number of chains: {len(target.record.chains)}")

                # Check inference options
                if target.record.inference_options:
                    inf_opts = target.record.inference_options
                    print(f"\n  Inference Options:")
                    if inf_opts.cdr3_beta_constraints:
                        print(f"    ✓ CDR3 beta constraints found!")
                        for i, (beta, regions) in enumerate(inf_opts.cdr3_beta_constraints):
                            print(f"      Constraint {i}: beta={beta}")
                            for chain_id, start, end in regions:
                                print(f"        - Chain {chain_id}: residues {start}-{end}")
                    else:
                        print(f"    ✗ No CDR3 beta constraints found")
                else:
                    print(f"  ✗ No inference options")
        except Exception as e:
            print(f"\n✗ Error parsing YAML: {e}")
            import traceback
            traceback.print_exc()
    else:
        print(f"✗ YAML file not found at {yaml_file}")

except ImportError as e:
    print(f"✗ Could not import parse_yaml: {e}")
    print("\nTrying alternative import...")

    try:
        from boltz.data.parse import schema
        yaml_file = repo_root / "yaml_cdr3_beta_scaling/7TRH_HBG_cdr3_beta.yml"

        print(f"Testing with schema module: {yaml_file}")
        result = schema.parse_yaml(str(yaml_file))
        print("✓ Schema parsing successful!")
    except Exception as e2:
        print(f"✗ Schema parsing also failed: {e2}")
