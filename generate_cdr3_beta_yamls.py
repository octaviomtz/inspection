#!/usr/bin/env python3
"""Generate YAML files with CDR3 beta scaling for all complexes."""

import csv
import os
from pathlib import Path
from typing import Dict, List, Tuple

# Read CDR information from CSV
def read_cdrs_csv(csv_path: str) -> Dict[str, Dict]:
    """Read CDR information from CSV file."""
    cdrs_data = {}
    with open(csv_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            complex_name = row['complex']
            # Parse CDR indices from string representation of lists
            # e.g., "[97, 98, 99]" -> [97, 98, 99]
            cdr3_h = eval(row['cdr3_h'])  # Heavy chain CDR3
            cdr3_l = eval(row['cdr3_l'])  # Light chain CDR3

            cdrs_data[complex_name] = {
                'heavy': row['heavy'],
                'light': row['light'],
                'cdr3_h': cdr3_h,
                'cdr3_l': cdr3_l,
            }
    return cdrs_data


def read_original_yaml(yaml_path: str) -> Dict:
    """Read original YAML file to get sequences."""
    import yaml
    with open(yaml_path, 'r') as f:
        data = yaml.safe_load(f)
    return data


def get_cdr3_ranges(cdr_indices: List[int]) -> Tuple[int, int]:
    """Convert list of indices to start and end residue numbers."""
    # CDR indices are 0-indexed, but YAML uses 1-indexed
    start = min(cdr_indices) + 1
    end = max(cdr_indices) + 1
    return start, end


def generate_cdr3_beta_yaml(
    complex_name: str,
    original_yaml: Dict,
    cdr_info: Dict,
    output_path: str,
    beta_value: float = 0.3,
) -> None:
    """Generate YAML with CDR3 beta scaling constraint."""

    # Get CDR3 ranges
    h3_start, h3_end = get_cdr3_ranges(cdr_info['cdr3_h'])
    l3_start, l3_end = get_cdr3_ranges(cdr_info['cdr3_l'])

    # Build YAML content
    yaml_lines = ["sequences:"]

    # Add sequences from original YAML
    if 'sequences' in original_yaml:
        for i, seq in enumerate(original_yaml['sequences']):
            yaml_lines.append(f"- {list(seq.keys())[0]}:")
            for key, value in seq[list(seq.keys())[0]].items():
                if key == 'msa':
                    # Update MSA path to use antigen_cut_new for antigen
                    if i == 0:  # Assume first chain is antigen
                        msa_path = f"../msa_antigen_cut_new/antigen_{complex_name}.a3m"
                    else:
                        msa_path = value
                    yaml_lines.append(f"    {key}: {msa_path}")
                else:
                    yaml_lines.append(f"    {key}: {value}")

    # Add CDR3 beta scaling constraint
    yaml_lines.append("")
    yaml_lines.append("constraints:")
    yaml_lines.append("  - cdr3_beta_scaling:")
    yaml_lines.append(f"      beta: {beta_value}")
    yaml_lines.append("      cdr3_regions:")
    yaml_lines.append(f"        - chain: B          # Heavy chain CDR3: residues {h3_start}-{h3_end}")
    yaml_lines.append(f"          start_res: {h3_start}")
    yaml_lines.append(f"          end_res: {h3_end}")
    yaml_lines.append(f"        - chain: C          # Light chain CDR3: residues {l3_start}-{l3_end}")
    yaml_lines.append(f"          start_res: {l3_start}")
    yaml_lines.append(f"          end_res: {l3_end}")

    # Write file
    with open(output_path, 'w') as f:
        f.write('\n'.join(yaml_lines) + '\n')


def main():
    """Generate YAML files for all complexes."""
    # Paths
    repo_root = Path(__file__).parent
    yaml_test_dir = repo_root / "yaml_test_set_boltz_vhvl_antigen_cut"
    examples_dir = repo_root / "examples"
    cdrs_csv = examples_dir / "cdrs.csv"

    # Output directory
    output_dir = repo_root / "yaml_cdr3_beta_scaling"
    output_dir.mkdir(exist_ok=True)

    # Read CDR information
    print("Reading CDR information from CSV...")
    cdrs_data = read_cdrs_csv(str(cdrs_csv))
    print(f"Found {len(cdrs_data)} complexes in CSV")

    # Get list of YAML files to process
    yaml_files = sorted(yaml_test_dir.glob("*.yml"))
    print(f"Found {len(yaml_files)} YAML files in test set")

    # Process each YAML file
    success_count = 0
    for yaml_file in yaml_files:
        complex_name = yaml_file.stem  # Remove .yml extension

        if complex_name not in cdrs_data:
            print(f"⚠ Skipping {complex_name}: not found in CDR CSV")
            continue

        try:
            # Read original YAML
            original_yaml = read_original_yaml(str(yaml_file))

            # Generate new YAML with CDR3 beta scaling
            output_path = output_dir / f"{complex_name}_cdr3_beta.yml"
            generate_cdr3_beta_yaml(
                complex_name,
                original_yaml,
                cdrs_data[complex_name],
                str(output_path),
                beta_value=0.3,
            )
            print(f"✓ Created {output_path.name}")
            success_count += 1

        except Exception as e:
            print(f"✗ Error processing {complex_name}: {e}")

    print(f"\nSuccessfully created {success_count}/{len(yaml_files)} YAML files")
    print(f"Output directory: {output_dir}")


if __name__ == "__main__":
    main()
