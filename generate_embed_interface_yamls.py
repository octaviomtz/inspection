#!/usr/bin/env python3
"""
Generate YAML files for the embedding interface steering feature.

Reads CDR data from cdrs.csv and existing YAML files, then creates new YAMLs
with embedding_interface constraints using CDR3 regions.
"""

import csv
import os
import ast
import yaml

BASE_DIR = "/mnt/c/Users/octav/Documents/claude/proteinEBM/worktree_steering_cdr3/w_embed_interface"
CDR_CSV = os.path.join(BASE_DIR, "examples", "cdrs.csv")
INPUT_YAML_DIR = os.path.join(BASE_DIR, "yaml_test_set_boltz_vhvl_antigen_cut")
OUTPUT_YAML_DIR = os.path.join(BASE_DIR, "yaml_embed_interface")


def parse_cdr_list(cdr_string):
    """Parse a CDR column string like '[97, 98, 99, 100]' into a list of ints."""
    return ast.literal_eval(cdr_string.strip())


def read_cdr_data(csv_path):
    """Read CDR data from CSV, returning a dict keyed by complex name."""
    cdr_data = {}
    with open(csv_path, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            complex_name = row["complex"]
            cdr3_h = parse_cdr_list(row["cdr3_h"])
            cdr3_l = parse_cdr_list(row["cdr3_l"])
            cdr_data[complex_name] = {
                "cdr3_h_start": min(cdr3_h),
                "cdr3_h_end": max(cdr3_h),
                "cdr3_l_start": min(cdr3_l),
                "cdr3_l_end": max(cdr3_l),
            }
    return cdr_data


def read_original_yaml(yaml_path):
    """Read an original YAML file and extract sequences by chain ID."""
    with open(yaml_path, "r") as f:
        data = yaml.safe_load(f)

    sequences = {}
    for entry in data["sequences"]:
        protein = entry["protein"]
        chain_id = protein["id"]
        sequences[chain_id] = protein["sequence"]
    return sequences


def generate_yaml(complex_name, sequences, cdr3_info):
    """Generate the output YAML structure with embedding_interface constraints."""
    output = {
        "sequences": [
            {
                "protein": {
                    "id": "A",
                    "sequence": sequences["A"],
                    "msa": f"examples/msa_antigen_cut_new/antigen_{complex_name}.a3m",
                }
            },
            {
                "protein": {
                    "id": "B",
                    "sequence": sequences["B"],
                    "msa": "empty",
                }
            },
            {
                "protein": {
                    "id": "C",
                    "sequence": sequences["C"],
                    "msa": "empty",
                }
            },
        ],
        "constraints": [
            {
                "embedding_interface": {
                    "antigen_chain": "A",
                    "contact_threshold": 8.0,
                    "cdr3_regions": [
                        {
                            "chain": "B",
                            "start_res": cdr3_info["cdr3_h_start"],
                            "end_res": cdr3_info["cdr3_h_end"],
                        },
                        {
                            "chain": "C",
                            "start_res": cdr3_info["cdr3_l_start"],
                            "end_res": cdr3_info["cdr3_l_end"],
                        },
                    ],
                    "force": True,
                }
            }
        ],
    }
    return output


def write_yaml(output_path, data):
    """Write YAML data to file with clean formatting."""
    with open(output_path, "w") as f:
        yaml.dump(data, f, default_flow_style=False, sort_keys=False, width=200)


def main():
    # Create output directory
    os.makedirs(OUTPUT_YAML_DIR, exist_ok=True)

    # Read CDR data
    cdr_data = read_cdr_data(CDR_CSV)
    print(f"Read CDR data for {len(cdr_data)} complexes from cdrs.csv")

    # Get all existing YAML files
    yaml_files = sorted([f for f in os.listdir(INPUT_YAML_DIR) if f.endswith(".yml")])
    print(f"Found {len(yaml_files)} YAML files in input directory")

    generated = 0
    missing_from_csv = []

    for yaml_file in yaml_files:
        complex_name = yaml_file.replace(".yml", "")
        input_path = os.path.join(INPUT_YAML_DIR, yaml_file)

        if complex_name not in cdr_data:
            missing_from_csv.append(complex_name)
            print(f"  WARNING: {complex_name} not found in cdrs.csv, skipping")
            continue

        # Read original sequences
        sequences = read_original_yaml(input_path)

        # Generate new YAML with embedding_interface constraint
        output_data = generate_yaml(complex_name, sequences, cdr_data[complex_name])

        # Write output
        output_path = os.path.join(OUTPUT_YAML_DIR, yaml_file)
        write_yaml(output_path, output_data)
        generated += 1

    print(f"\n--- Summary ---")
    print(f"Total YAML files in input directory: {len(yaml_files)}")
    print(f"Files generated: {generated}")
    print(f"Complexes missing from cdrs.csv: {len(missing_from_csv)}")
    if missing_from_csv:
        for name in missing_from_csv:
            print(f"  - {name}")


if __name__ == "__main__":
    main()
