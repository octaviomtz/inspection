#!/usr/bin/env python3
"""Generate combined YAML files for Experiment 1: Y early-phase only + B2 Contact Restraints.

Each output YAML adds a hierarchical_steering constraint (CDR3 beta-scaling, early-phase only)
to the existing B2 contact restraint YAML files.

Usage:
    python scripts/generate_y_early_b2_yamls.py
"""

import ast
import csv
import os
import re
import sys
from pathlib import Path


# Paths
REPO_ROOT = Path(__file__).parent.parent
CDR_CSV = Path("/mnt/c/Users/octav/Documents/claude/proteinEBM/boltz/examples/cdrs.csv")
B2_YAML_DIR = REPO_ROOT / "yaml_test_set_boltz_contact_restraints_msa_vhvl_antigen_cut"
OUTPUT_DIR = REPO_ROOT / "yaml_test_set_y_early_b2_contact_restraints_msa_vhvl_antigen_cut"

# Y early-phase parameters (Experiment 1: beta_max in recommended 0.3-0.5 range)
BETA_MAX = 0.4
BETA_SCHEDULE = "linear"
CONTACT_THRESHOLD = 8.0
GUIDANCE_WEIGHT_SCALE = 1.0
ANTIGEN_CHAIN = "A"
HEAVY_CHAIN = "B"
LIGHT_CHAIN = "C"


def load_cdr_data(csv_path: Path) -> dict:
    """Load CDR3 regions from cdrs.csv, keyed by complex name."""
    cdr_data = {}
    with open(csv_path) as f:
        for row in csv.DictReader(f):
            complex_name = row["complex"]
            cdr3_h = ast.literal_eval(row["cdr3_h"])
            cdr3_l = ast.literal_eval(row["cdr3_l"])
            cdr_data[complex_name] = {
                "cdr3_h_start": min(cdr3_h),
                "cdr3_h_end": max(cdr3_h),
                "cdr3_l_start": min(cdr3_l),
                "cdr3_l_end": max(cdr3_l),
            }
    return cdr_data


def make_hierarchical_block(cdr: dict) -> str:
    """Return the YAML text for the hierarchical_steering constraint."""
    return (
        f"- hierarchical_steering:\n"
        f"    antigen_chain: {ANTIGEN_CHAIN}\n"
        f"    contact_threshold: {CONTACT_THRESHOLD}\n"
        f"    cdr_regions:\n"
        f"      - chain: {HEAVY_CHAIN}\n"
        f"        start_res: {cdr['cdr3_h_start']}\n"
        f"        end_res: {cdr['cdr3_h_end']}\n"
        f"      - chain: {LIGHT_CHAIN}\n"
        f"        start_res: {cdr['cdr3_l_start']}\n"
        f"        end_res: {cdr['cdr3_l_end']}\n"
        f"    beta_max: {BETA_MAX}\n"
        f"    beta_schedule: {BETA_SCHEDULE}\n"
        f"    guidance_weight_scale: {GUIDANCE_WEIGHT_SCALE}\n"
        f"    force: true\n"
    )


def process_yaml(input_path: Path, output_path: Path, cdr: dict) -> None:
    """Read a B2 YAML, append hierarchical_steering constraint, write to output."""
    content = input_path.read_text()

    # The hierarchical_steering block goes after the existing constraints block
    hierarchical_block = make_hierarchical_block(cdr)

    # Append to the constraints section
    # The existing YAML ends with a contact constraint; we add our block after it
    if "constraints:" in content:
        content = content.rstrip() + "\n" + hierarchical_block
    else:
        # No constraints section yet — add one
        content = content.rstrip() + "\nconstraints:\n" + hierarchical_block

    output_path.write_text(content)


def main():
    cdr_data = load_cdr_data(CDR_CSV)
    print(f"Loaded CDR3 data for {len(cdr_data)} complexes")

    OUTPUT_DIR.mkdir(exist_ok=True)

    yaml_files = sorted(B2_YAML_DIR.glob("*.yml"))
    print(f"Found {len(yaml_files)} input YAML files")

    skipped = []
    processed = 0

    for yaml_file in yaml_files:
        # Extract complex code from filename: restraint_7TRH_HBG_hbond_19.yml → 7TRH_HBG
        m = re.match(r"restraint_([^_]+_[^_]+)_", yaml_file.name)
        if not m:
            print(f"  WARNING: could not parse complex from {yaml_file.name}, skipping")
            skipped.append(yaml_file.name)
            continue

        complex_code = m.group(1)
        if complex_code not in cdr_data:
            print(f"  WARNING: {complex_code} not found in cdrs.csv, skipping")
            skipped.append(yaml_file.name)
            continue

        # Output filename: same name, different folder
        output_file = OUTPUT_DIR / yaml_file.name
        process_yaml(yaml_file, output_file, cdr_data[complex_code])
        processed += 1

    print(f"\nDone: {processed} files written to {OUTPUT_DIR}")
    if skipped:
        print(f"Skipped {len(skipped)} files: {skipped}")


if __name__ == "__main__":
    main()
