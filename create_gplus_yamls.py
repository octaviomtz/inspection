#!/usr/bin/env python3
"""
Generate Strategy G+ YAML files for all complexes.

This script:
1. Reads the original YAML test set files
2. Extracts CDR information from cdrs.csv
3. Creates new YAML files with antigen_orientation constraints
4. Updates MSA paths to use examples/msa_antigen_cut_new/
"""

import csv
import yaml
from pathlib import Path
from typing import Dict, List, Tuple

# Define paths
YAML_TEST_SET_DIR = Path("yaml_test_set_boltz_vhvl_antigen_cut")
OUTPUT_DIR = Path("yaml_test_set_gplus")
CDR_CSV_FILE = Path("examples/cdrs.csv")
MSA_PATH_TEMPLATE = "examples/msa_antigen_cut_new/antigen_{}.a3m"


def load_cdr_data() -> Dict[str, Dict]:
    """Load CDR information from CSV file."""
    cdr_data = {}
    with open(CDR_CSV_FILE, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            complex_name = row['complex']
            # Parse CDR ranges from string representation like "[26, 27, 28, ...]"
            cdr3_h = eval(row['cdr3_h'])  # noqa: S307
            cdr3_l = eval(row['cdr3_l'])  # noqa: S307

            cdr_data[complex_name] = {
                'cdr3_h_start': cdr3_h[0],
                'cdr3_h_end': cdr3_h[-1],
                'cdr3_l_start': cdr3_l[0],
                'cdr3_l_end': cdr3_l[-1],
            }
    return cdr_data


def load_original_yaml(complex_name: str) -> Dict:
    """Load the original YAML file for a complex."""
    yaml_file = YAML_TEST_SET_DIR / f"{complex_name}.yml"
    if not yaml_file.exists():
        print(f"Warning: YAML file not found for {complex_name}: {yaml_file}")
        return None

    with open(yaml_file, 'r') as f:
        return yaml.safe_load(f)


def create_gplus_yaml(
    original_yaml: Dict,
    complex_name: str,
    cdr_data: Dict,
) -> Dict:
    """Create a new YAML with Strategy G+ constraints."""
    gplus_yaml = yaml.safe_load(yaml.dump(original_yaml))  # Deep copy

    # Update MSA path for antigen chain (chain A)
    for seq in gplus_yaml['sequences']:
        if 'protein' in seq:
            protein = seq['protein']
            if protein.get('id') == 'A':
                protein['msa'] = MSA_PATH_TEMPLATE.format(complex_name)

    # Create antigen_orientation constraint with CDR3 information
    cdr_info = cdr_data[complex_name]

    constraint = {
        'antigen_orientation': {
            'antigen_chain': 'A',
            'contact_threshold': 8.0,
            'cdr3_regions': [
                {
                    'chain': 'B',  # Heavy chain CDR3
                    'start_res': cdr_info['cdr3_h_start'],
                    'end_res': cdr_info['cdr3_h_end'],
                },
                {
                    'chain': 'C',  # Light chain CDR3
                    'start_res': cdr_info['cdr3_l_start'],
                    'end_res': cdr_info['cdr3_l_end'],
                },
            ],
            'force': True,
        }
    }

    # Add constraints section if it doesn't exist
    if 'constraints' not in gplus_yaml:
        gplus_yaml['constraints'] = []
    else:
        gplus_yaml['constraints'] = list(gplus_yaml['constraints'])  # Ensure it's a list

    gplus_yaml['constraints'].append(constraint)

    return gplus_yaml


def main():
    """Main function to generate all YAML files."""
    # Create output directory
    OUTPUT_DIR.mkdir(exist_ok=True)
    print(f"Creating output directory: {OUTPUT_DIR}")

    # Load CDR data
    cdr_data = load_cdr_data()
    print(f"Loaded CDR data for {len(cdr_data)} complexes")

    # Get all YAML files from test set
    yaml_files = sorted(YAML_TEST_SET_DIR.glob("*.yml"))
    print(f"Found {len(yaml_files)} YAML files in test set")

    # Process each complex
    created = 0
    skipped = 0

    for yaml_file in yaml_files:
        complex_name = yaml_file.stem  # Remove .yml extension

        # Check if we have CDR data for this complex
        if complex_name not in cdr_data:
            print(f"⚠ Skipping {complex_name}: No CDR data found")
            skipped += 1
            continue

        # Load original YAML
        original_yaml = load_original_yaml(complex_name)
        if original_yaml is None:
            skipped += 1
            continue

        # Create G+ YAML
        gplus_yaml = create_gplus_yaml(original_yaml, complex_name, cdr_data)

        # Write output file
        output_file = OUTPUT_DIR / f"{complex_name}_strategy_gplus.yml"
        with open(output_file, 'w') as f:
            yaml.dump(gplus_yaml, f, default_flow_style=False, sort_keys=False)

        print(f"✓ Created {output_file.name}")
        created += 1

    # Summary
    print("\n" + "=" * 60)
    print(f"Summary:")
    print(f"  Created: {created} YAML files")
    print(f"  Skipped: {skipped} files")
    print(f"  Output directory: {OUTPUT_DIR}")
    print("=" * 60)


if __name__ == "__main__":
    main()
