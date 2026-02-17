# K+ Region-Specific β-Scaling YAML Test Set

This folder contains YAML configuration files for testing Method K+ (epitope region scanning) on 48 antibody-antigen complexes.

## Configuration Details

Each YAML file includes:

- **Protein sequences** for antigen (A), heavy chain (B), and light chain (C)
- **MSA path** for the antigen: `examples/msa_antigen_cut_new/antigen_<complex_name>.a3m`
- **Epitope region scanning constraint** with the following parameters:
  - `antigen_chain`: A (target antigen)
  - `num_regions`: 10 (partition antigen into 10 overlapping regions)
  - `region_overlap_ratio`: 0.2 (20% overlap between adjacent regions)
  - `beta_emphasis`: 0.5 (boost pair representations for emphasized region)
  - `beta_deemphasis`: -0.3 (suppress pair representations for other regions)
  - `contact_threshold`: 8.0 Å (distance threshold for contact detection)
  - `confidence_weighting`: true (weight contacts by pLDDT/pAE)
  - `accumulate_heatmap`: true (generate epitope propensity heatmap)
  - `save_per_region`: false (don't save individual region predictions)

## Usage

Run Boltz2 with K+ constraint:
```bash
python src/boltz/main.py --config yaml_test_set_boltz_vhvl_antigen_cut_k_plus/<complex_name>.yml
```

Or use CLI flags:
```bash
python src/boltz/main.py \
  --config yaml_test_set_boltz_vhvl_antigen_cut_k_plus/<complex_name>.yml \
  --epitope_scanning \
  --scan_regions 10
```

## Complexes Included

48 PDB complexes covering diverse epitopes:

- 7TRH_HBG, 7TRI_ZYB, 7WT9_HAE, 7Y0O_HLA, 7ZOZ_HLA
- 8BLQ_ABD, 8BLQ_ECD, 8BYU_HLA, 8CDD_EDB, 8CDE_DCB
- 8CXC_HLM, 8CYH_HLM, 8DE3_BCA, 8DTK_CBA, 8E2U_HLA
- 8EAY_HLA, 8EQ6_HLA, 8EZ3_HLA, 8EZ7_HLA, 8EZ8_HLA
- 8FAH_HLA, 8FGX_BAC, 8FXB_HLE, 8GH4_HLE, 8GHP_HLA
- 8GP5_EFX, 8GQ1_HLC, 8HES_HLC, 8HGM_CDB, 8HLB_BCA
- 8HLB_DEA, 8IDN_HLA, 8IV4_ABG, 8IV4_HLG, 8IV5_ABG
- 8IVA_CEG, 8IVX_HLA, 8IX3_HLG, 8J1T_HKF, 8OL9_BAH
- 8OXW_BCA, 8OXX_BCA, 8PE9_HLA, 8SLB_HLA, 8T9Z_HLA
- 8TFR_ABC, 8TRS_AGD, 8X0T_HLA

## Expected Outputs

For each complex, Boltz2 with K+ will generate:
- `<complex_name>_epitope_heatmap.json` - Propensity values and statistics
- `<complex_name>_epitope_hotspots.json` - High-confidence hotspot residues
- `<complex_name>_epitope_statistics.json` - Scanning metrics and coverage statistics
- Standard PDB files from the final prediction

## CDR Information

CDR sequences for each antibody are documented in `examples/cdrs.csv`:
- CDR-H1, CDR-H2, CDR-H3 (heavy chain)
- CDR-L1, CDR-L2, CDR-L3 (light chain)

CDRs are automatically identified from the sequences during processing.
