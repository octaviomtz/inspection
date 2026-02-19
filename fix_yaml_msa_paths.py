#!/usr/bin/env python3
"""Fix MSA paths in generated YAML files."""

import yaml
from pathlib import Path

yaml_dir = Path("yaml_cdr3_beta_scaling")

for yaml_file in sorted(yaml_dir.glob("*.yml")):
    with open(yaml_file, 'r') as f:
        data = yaml.safe_load(f)

    # Fix MSA path for antigen (first sequence)
    if data.get("sequences"):
        for seq in data["sequences"]:
            for key, val in seq.items():
                if key == "protein" and val.get("id") == "A":
                    # Get complex name from filename
                    complex_name = yaml_file.stem.replace("_cdr3_beta", "")
                    val["msa"] = f"examples/msa_antigen_cut_new/antigen_{complex_name}.a3m"

    # Write back
    with open(yaml_file, 'w') as f:
        yaml.dump(data, f, default_flow_style=False, sort_keys=False)

    print(f"✓ Fixed: {yaml_file.name}")

print(f"\nUpdated all YAML files with correct MSA paths")
