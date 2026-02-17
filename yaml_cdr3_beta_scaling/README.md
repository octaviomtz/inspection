# CDR3 Beta Scaling YAML Files

This directory contains YAML configuration files for all 48 complexes in the test set, with **CDR3-Specific β-Scaling** constraints enabled.

## Overview

Each YAML file includes:
- **Antigen sequence** (chain A) with MSA path: `../msa_antigen_cut_new/antigen_<complex_name>.a3m`
- **Heavy chain antibody** (chain B) with empty MSA
- **Light chain antibody** (chain C) with empty MSA
- **CDR3 beta scaling constraint** with β = 0.3 for both heavy (H3) and light (L3) CDR3 regions

## File Naming Convention

All files follow the pattern: `<COMPLEX_NAME>_cdr3_beta.yml`

Example: `7TRH_HBG_cdr3_beta.yml`, `8CDD_EDB_cdr3_beta.yml`

## CDR3 Regions

CDR3 region boundaries vary by complex based on sequence-specific characteristics:

### Examples:

| Complex | H3 Range | L3 Range |
|---------|----------|----------|
| 7TRH_HBG | 97-115 | 88-100 |
| 8CDD_EDB | 98-116 | 86-94 |
| 8EZ3_HLA | 99-113 | 91-102 |

All CDR3 boundaries are automatically extracted from `examples/cdrs.csv` and converted to 1-indexed residue numbers for YAML specification.

## YAML Structure

Each file has the following structure:

```yaml
sequences:
- protein:
    id: A
    sequence: <antigen_sequence>
    msa: ../msa_antigen_cut_new/antigen_<complex_name>.a3m
- protein:
    id: B
    sequence: <heavy_chain_sequence>
    msa: empty
- protein:
    id: C
    sequence: <light_chain_sequence>
    msa: empty

constraints:
  - cdr3_beta_scaling:
      beta: 0.3
      cdr3_regions:
        - chain: B          # Heavy chain CDR3
          start_res: <H3_START>
          end_res: <H3_END>
        - chain: C          # Light chain CDR3
          start_res: <L3_START>
          end_res: <L3_END>
```

## Usage

Run predictions on any complex using:

```bash
boltz predict yaml_cdr3_beta_scaling/<COMPLEX_NAME>_cdr3_beta.yml
```

For example:

```bash
boltz predict yaml_cdr3_beta_scaling/7TRH_HBG_cdr3_beta.yml
boltz predict yaml_cdr3_beta_scaling/8CDD_EDB_cdr3_beta.yml
```

## Generation

These files were automatically generated using the `generate_cdr3_beta_yamls.py` script:

```bash
python generate_cdr3_beta_yamls.py
```

The script:
1. Reads CDR3 region information from `examples/cdrs.csv`
2. Reads baseline YAML files from `yaml_test_set_boltz_vhvl_antigen_cut/`
3. Generates new YAML files with CDR3 beta scaling constraints
4. Uses β = 0.3 as the default scaling factor (adjustable in script)

## Customization

To modify the beta value for all complexes:

1. Edit the `generate_cdr3_beta_yamls.py` script
2. Change the `beta_value` parameter in the `generate_cdr3_beta_yaml()` call
3. Re-run the script to regenerate all files

## Related Documentation

- **CDR3 Beta Scaling Feature**: See `STEERING_METHODS.md` (Idea L+, Phase 1)
- **Implementation Details**: See implementation guide in `notes_steering/`
- **CDR Information**: See `examples/cdrs.csv` for all CDR region indices

## Files Generated

Total: **48 YAML files** (one for each complex in the test set)

- 7TRH_HBG_cdr3_beta.yml
- 7TRI_ZYB_cdr3_beta.yml
- 7WT9_HAE_cdr3_beta.yml
- 7Y0O_HLA_cdr3_beta.yml
- 7ZOZ_HLA_cdr3_beta.yml
- 8BLQ_ABD_cdr3_beta.yml
- 8BLQ_ECD_cdr3_beta.yml
- 8BYU_HLA_cdr3_beta.yml
- ... and 40 more

See directory listing for complete list.
