# Evaluation Plan: Antigen Steering Feature

**Date**: 2026-03-06
**Status**: Planning
**Branch**: `steering_cdr3`

---

## 1. Goal

Evaluate whether the antigen steering feature (`--use_potentials --antigen_steering`) improves antibody-antigen interface quality compared to three baselines. The feature uses Feynman-Kac particle resampling + gradient-based physical guidance to steer antigen orientation toward CDR3 loops during diffusion.

The primary hypothesis: **steering the antigen toward CDR3 regions during diffusion produces better antibody-antigen interfaces than unconstrained prediction or contact-restraint-based approaches.**

---

## 2. Dataset

### Ground Truth
- **Location**: `pdb_minimized/`
- **Count**: 48 complexes
- **Format**: PDB
- **Chain convention**: A = antigen, B = heavy chain (VH), C = light chain (VL)
- **Naming**: `{PDB_ID}_{CHAINS}.pdb` (e.g., `7TRH_HBG.pdb`)

### Methods to Compare

| ID | Method | Location | Status | Format | Models/complex | Settings/complex |
|----|--------|----------|--------|--------|---------------|-----------------|
| B0 | Baseline 1: Unconstrained | `predictions_examples/antigen_cut/` | 2 of 48 done | PDB | 5 | 1 |
| B1 | Baseline 2: Contact restraints | `predictions_examples/antigen_cut_contact_restraints/` | 1 of 48 done (processed only, no predictions) | PDB | 5 | ~4 (hbond_N, hydrophobic_N, salt_bridge_N) |
| B2 | Baseline 3: Pocket restraints | `predictions_examples/antigen_cut_vhvl_msa_pocket_ab_boltz_post_2023_omm/` | 1 of 48 done (processed only, no predictions) | PDB | 5 | ~3 (restraint_to_A_{ID}_{chain}_{res}_{num}) |
| NF | New Feature: Antigen steering | `predictions_examples/new_feature_cdr3_beta/` | 1 of 48 done | CIF | 1 | 1 |

### Current Data Gaps

**Critical**: Most predictions have not been generated yet. Before evaluation can proceed:

1. **B0**: Run remaining 46 complexes unconstrained (no flags)
2. **B1**: Run predictions for all 48 complexes across all settings (processed data exists)
3. **B2**: Run predictions for all 48 complexes across all settings (processed data exists)
4. **NF**: Run remaining 47 complexes with `--use_potentials --antigen_steering --diffusion_samples 5 --num_particles 1 --no_kernels`

**NF output format**: Currently outputs CIF (1 model). Should be re-run with `--diffusion_samples 5 --output_format pdb` for consistency and model selection.

---

## 3. Prediction Generation Commands

### B0: Unconstrained Baseline
```bash
conda activate boltz
boltz predict examples/{PDB_ID}/{PDB_ID}_{CHAINS}.yml \
  --output_format pdb --diffusion_samples 5
```

### B1: Contact Restraints
```bash
conda activate boltz
boltz predict examples/{PDB_ID}/restraint_{PDB_ID}_{CHAINS}_{type}_{N}.yml \
  --output_format pdb --diffusion_samples 5
```

### B2: Pocket Restraints
```bash
conda activate boltz
boltz predict examples/{PDB_ID}/restraint_to_A_{PDB_ID}_{CHAINS}_{chain}_{res}_{num}.yml \
  --output_format pdb --diffusion_samples 5
```

### NF: Antigen Steering
```bash
conda activate boltz
boltz predict examples/{PDB_ID}/{PDB_ID}_{CHAINS}_antigen_steer.yml \
  --output_format pdb --use_potentials --antigen_steering \
  --diffusion_samples 5 --num_particles 1 --no_kernels
```

**Note**: `--num_particles 1` avoids OOM (otherwise total = diffusion_samples x num_particles). With `--num_particles 1`, FK resampling is disabled but gradient-based physical guidance still operates.

---

## 4. Metrics

### 4.1 Primary Metric: DockQ v2

**Tool**: DockQ v2 (`conda activate dockq2`)

**Command**:
```bash
conda run -n dockq2 DockQ <model> <native> --json <output.json>
```

DockQ produces per-interface scores. For a 3-chain complex (A, B, C), it reports 3 interfaces:
- **A-B** (antigen - heavy chain) -- **primary interface of interest**
- **A-C** (antigen - light chain) -- **secondary interface of interest**
- **B-C** (VH - VL) -- internal antibody interface, less relevant for steering evaluation

**DockQ Quality Thresholds** (CAPRI classification):
| Quality | DockQ range |
|---------|-------------|
| Incorrect | 0.00 - 0.23 |
| Acceptable | 0.23 - 0.49 |
| Medium | 0.49 - 0.80 |
| High | >= 0.80 |

**Per-interface sub-metrics**:
- `DockQ`: Combined quality score (0-1)
- `iRMSD`: Interface RMSD (lower is better)
- `LRMSD`: Ligand RMSD (lower is better)
- `fnat`: Fraction of native contacts recovered (higher is better)
- `fnonnat`: Fraction of non-native contacts (lower is better)
- `F1`: Harmonic mean of precision and recall of contacts

**CIF compatibility**: DockQ v2 handles CIF input natively (verified).

### 4.2 Interface-Specific Analysis

Since the steering feature targets the antibody-antigen interface specifically:

1. **Antibody-Antigen DockQ**: Average of A-B and A-C interface DockQ scores
2. **VH-VL DockQ**: B-C interface DockQ (should remain stable; degradation indicates steering damages antibody fold)

Report both separately. The feature should improve (1) without degrading (2).

### 4.3 Confidence Metrics (from Boltz output)

Available in `confidence_{ID}_model_{N}.json`:
- `confidence_score`: Overall confidence (0.8 * iptm + 0.2 * ptm)
- `iptm`: Interface pTM (inter-chain predicted TM-score)
- `protein_iptm`: Same as iptm for protein-only complexes
- `complex_plddt`: Average pLDDT across all residues
- `complex_iplddt`: Interface pLDDT (pLDDT at interface residues)
- `pair_chains_iptm`: Per-chain-pair iptm matrix (chains indexed 0, 1, 2)

**Chain index mapping**: 0 = A (antigen), 1 = B (heavy), 2 = C (light)

Key confidence metrics to extract:
- `pair_chains_iptm["0"]["1"]` = antigen-heavy iptm
- `pair_chains_iptm["0"]["2"]` = antigen-light iptm
- `pair_chains_iptm["1"]["2"]` = VH-VL iptm

### 4.4 Epitope Prediction Accuracy (Optional)

Extract predicted interface residues from each model and compare to ground truth interface:
1. **Ground truth epitope**: Antigen residues within 8A of any antibody atom in native structure
2. **Predicted epitope**: Antigen residues within 8A of any antibody atom in predicted model
3. **Metrics**: Precision, recall, F1 of epitope residue prediction

---

## 5. Model Selection Strategy

When multiple models are generated per complex (e.g., 5 diffusion samples), we need a strategy to select the "best" model for evaluation. Two approaches:

### 5.1 Best-of-N by Confidence (Practical)
Select the model with the highest `confidence_score` (= 0.8 * iptm + 0.2 * ptm).

**Rationale**: This is what a user would do in practice without access to ground truth.

### 5.2 Best-of-N by DockQ (Oracle)
Select the model with the highest antibody-antigen DockQ.

**Rationale**: Upper bound on performance -- shows the best the method CAN produce.

### 5.3 Recommendation
Report **both**. The gap between confidence-selected and oracle-selected reveals how well confidence correlates with actual quality. A smaller gap means the method is more reliable in practice.

---

## 6. Handling Multiple Settings (Baselines 2 and 3)

Baselines 2 and 3 have multiple restraint settings per complex:

**B1 settings per complex**: ~4 (e.g., hbond_23, hbond_7, hydrophobic_5, salt_bridge_3)
**B2 settings per complex**: ~3 (e.g., restraint_to_A_{chain}_{res}_{num})

### Aggregation Strategy

For each complex in B1/B2:
1. Run all settings with 5 diffusion samples each
2. For each setting, select best model by confidence (Section 5.1)
3. **Best-across-settings**: Report the setting that produced the best DockQ (oracle selection across settings)
4. **Confidence-across-settings**: Report the setting whose best model had highest confidence

This gives B1 and B2 a fair advantage: they get to pick the best restraint type. This makes the comparison conservative for NF -- if NF still wins, the result is strong.

---

## 7. Evaluation Pipeline

### Step 1: Generate All Predictions

Run predictions for all 48 complexes across all methods (see Section 3).

### Step 2: Run DockQ on All Predictions

```bash
#!/bin/bash
# evaluate_dockq.sh
# Run from project root

GROUND_TRUTH_DIR="pdb_minimized"
OUTPUT_DIR="evaluation_results"
mkdir -p "$OUTPUT_DIR"

# For each method and each complex:
for MODEL_PDB in <prediction_path>/*.pdb; do
    COMPLEX_ID=$(basename "$MODEL_PDB" | sed 's/_model_[0-9]*.pdb//')
    NATIVE="$GROUND_TRUTH_DIR/${COMPLEX_ID}.pdb"

    if [ -f "$NATIVE" ]; then
        conda run -n dockq2 DockQ "$MODEL_PDB" "$NATIVE" \
          --json "$OUTPUT_DIR/dockq_$(basename $MODEL_PDB .pdb).json"
    fi
done
```

### Step 3: Extract and Aggregate Metrics

For each complex and method:
1. Parse DockQ JSON outputs
2. Extract per-interface metrics (A-B, A-C, B-C)
3. Apply model selection (confidence-based and oracle)
4. Compute antibody-antigen DockQ = mean(DockQ_AB, DockQ_AC)

### Step 4: Statistical Comparison

For each metric, compare NF vs each baseline:
1. **Per-complex comparison**: For each complex, compute delta = NF_metric - baseline_metric
2. **Aggregate**: Mean, median, std of deltas across all 48 complexes
3. **Significance**: Wilcoxon signed-rank test (paired, non-parametric)
4. **CAPRI classification**: Percentage of complexes in each quality tier per method

---

## 8. Output Tables

### Table 1: Per-Method Aggregate Results
| Method | DockQ (AB+AC) | DockQ_AB | DockQ_AC | DockQ_BC | iRMSD_AB | fnat_AB | iptm |
|--------|--------------|----------|----------|----------|----------|---------|------|
| B0: Unconstrained | | | | | | | |
| B1: Contact restraints (best) | | | | | | | |
| B2: Pocket restraints (best) | | | | | | | |
| NF: Antigen steering | | | | | | | |

### Table 2: CAPRI Quality Distribution
| Method | Incorrect (%) | Acceptable (%) | Medium (%) | High (%) |
|--------|--------------|----------------|------------|----------|
| B0 | | | | |
| B1 | | | | |
| B2 | | | | |
| NF | | | | |

### Table 3: Per-Complex Detailed Results
One row per complex, columns for each method's DockQ_AB and DockQ_AC.

---

## 9. Environments

| Environment | Purpose | Activation |
|-------------|---------|------------|
| `boltz` | Run predictions | `conda activate boltz` |
| `dockq2` | Run DockQ evaluation | `conda activate dockq2` |

---

## 10. Key Considerations

### Format Consistency
- B0, B1, B2 output PDB; NF currently outputs CIF
- DockQ handles both, so no conversion needed
- For consistency, re-run NF with `--output_format pdb`

### Fair Comparison
- All methods should use the same `--diffusion_samples 5` for equal sampling budget
- B1 and B2 get advantage of multiple settings (best-across-settings selection)
- NF uses `--num_particles 1` (no FK resampling) to avoid OOM; if GPU memory allows, test with `--num_particles 3` as well

### Antibody Internal Quality
- Monitor B-C (VH-VL) DockQ across methods. If steering degrades antibody fold quality, this is a significant concern even if interface quality improves.

### Confidence Calibration
- Compare confidence scores to actual DockQ. A well-calibrated method should have high correlation between confidence and quality.
- Report Spearman rank correlation of confidence_score vs DockQ for each method.

---

## 11. Stretch Goals

1. **Epitope accuracy analysis**: Which method best predicts native epitope residues?
2. **CDR3 RMSD**: Align on framework, compute CDR3 backbone RMSD to native
3. **Contact map visualization**: Plot predicted vs native contact maps per method
4. **Effect of num_particles**: Compare NF with num_particles=1 vs 3 (if feasible)
5. **Per-complex scatter plots**: DockQ_NF vs DockQ_B0 for each complex to identify winners/losers
