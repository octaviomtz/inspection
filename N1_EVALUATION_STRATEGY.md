# N1: A+Y Hybrid — Evaluation Strategy

**Date**: 2026-03-26
**Feature**: N1 (FK Particles + Hierarchical Timing)
**Branch**: `a_v2_fk_particles`

---

## 1. Evaluation Goals

N1 combines Strategy A (FK particle resampling, best DockQ: +0.046) with Strategy Y (hierarchical timing, only statistically significant CDR-H3 RMSD: -0.30A, p=0.025). The evaluation must answer:

1. **Does the combination improve DockQ over both individual strategies and baselines?** (target: +0.05-0.08 vs B0)
2. **Does CDR-H3 RMSD improve?** (target: -0.3-0.5A vs B0)
3. **Is the improvement statistically significant?** (paired Wilcoxon, p < 0.05)
4. **Does it avoid regressions on easy cases?** (per-complex delta analysis)
5. **How does confidence calibration behave?** (Spearman rho between confidence and DockQ)

---

## 2. Methods to Compare

| ID | Name | Prediction Folder | Description |
|----|------|-------------------|-------------|
| B0 | Unconstrained | `predictions_examples/antigen_cut/` | Standard Boltz2 — 5 models per complex |
| B1 | Contact restraints | `predictions_examples/antigen_cut_contact_restraints/` | Known contact pairs as restraints — multiple settings per complex |
| B2 | Pocket restraints | `predictions_examples/antigen_cut_vhvl_msa_pocket_ab_boltz_post_2023_omm/` | Pocket-based restraints — multiple settings per complex |
| N1 | Hybrid FK+Hierarchical | `predictions_examples/new_feature_hybrid_fk/` | This feature — FK particles + hierarchical timing |

**Note on N1 predictions**: The folder `new_feature_hybrid_fk/` does not yet exist. When running inference, use:
```bash
conda run -n boltz boltz predict \
  yaml_hybrid_fk_hierarchical/<CID>.yml \
  --use_potentials \
  --hybrid_fk_hierarchical \
  --hybrid_late_particles 5 \
  --num_samples 5 \
  --out_dir predictions_examples/new_feature_hybrid_fk/
```
This will produce 5 models per complex, enabling proper confidence-based model selection.

---

## 3. Prediction Folder Structure

### B0 (Unconstrained)
```
antigen_cut/
  boltz_results_{CID}/
    predictions/{CID}/
      {CID}_model_0.pdb ... {CID}_model_4.pdb
      confidence_{CID}_model_0.json ...
```
- 1 setting per complex, 5 PDB models, 5 confidence JSONs.

### B1 (Contact Restraints)
```
antigen_cut_contact_restraints/
  boltz_results_restraint_{CID}_{type}_{N}/
    predictions/restraint_{CID}_{type}_{N}/
      restraint_{CID}_{type}_{N}_model_0.pdb ...
      confidence_restraint_{CID}_{type}_{N}_model_0.json ...
```
- Multiple settings per complex (hbond_23, hbond_7, hydrophobic_5, salt_bridge_3, etc.).
- Each setting has 5 PDB models.
- For evaluation: treat each setting as a separate "model" — oracle selection picks best across all settings and models.

### B2 (Pocket Restraints)
```
antigen_cut_vhvl_msa_pocket_ab_boltz_post_2023_omm/
  boltz_results_restraint_to_A_{CID}_{chain}_{res}_{num}/
    predictions/restraint_to_A_{CID}_{chain}_{res}_{num}/
      restraint_to_A_{CID}_{chain}_{res}_{num}_model_0.pdb ...
```
- Multiple settings per complex (different pocket residues).
- Each setting has 5 PDB models.

### N1 (Hybrid FK+Hierarchical) — Expected
```
new_feature_hybrid_fk/
  boltz_results_{CID}/
    predictions/{CID}/
      {CID}_model_0.cif ... {CID}_model_4.cif
      confidence_{CID}_model_0.json ...
```
- 1 setting per complex, 5 CIF models (Boltz2 outputs CIF by default with steering flags), 5 confidence JSONs.
- **Important**: The reference `new_feature_hierarchical/` folder (Strategy Y) uses CIF format and `_hierarchical` suffix. N1 should NOT use a suffix to keep naming clean — but the discovery function must handle both PDB and CIF formats.

---

## 4. Can We Reuse the FK v2 evaluate.py?

**Short answer**: Yes, with modifications. The reference `evaluate.py` (from `../a_v2_fk_particles_true/scripts/`) provides ~80% of the needed infrastructure.

### What can be reused as-is
- DockQ execution via subprocess with JSON caching
- Model evaluation (confidence parsing, DockQ metric extraction)
- Aggregation functions (summary tables, CAPRI classification, Wilcoxon tests)
- Plotting functions (aggregate bars, CAPRI stacked, boxplot, scatter, waterfall, heatmap)
- Confidence calibration analysis

### What needs to change

| Component | Current (FK v2) | Needed for N1 |
|-----------|-----------------|---------------|
| Method IDs | B0, B1, B2, NF_V2 | B0, B1, B2, N1 |
| Method subfolder | `new_feature_fk_v2` | `new_feature_hybrid_fk` |
| Discovery pattern (N1) | FK-specific (`fk_v2` suffix) | Hybrid-specific (no suffix or `_hybrid_fk` suffix) |
| Model selection | confidence, oracle, composite (FK energy) | confidence, oracle (no composite — FK energy may not be saved in the same way) |
| CDR-H3 RMSD | Not computed | **Must add** — core metric for validating Y's contribution |
| File format handling | PDB + CIF | Same (already handles both via `_find_prediction_files`) |

### Recommendation

**Copy and adapt** the FK v2 evaluate.py rather than modifying it in place. This keeps the FK v2 evaluation intact and avoids branch conflicts. Create `scripts/evaluate_n1.py`.

---

## 5. Metrics

### 5.1 Primary Metrics (from DockQ)

DockQ v2 computes per-interface metrics. For antibody-antigen complexes with chains A (antigen), B (heavy), C (light):

| Metric | Interfaces | Description |
|--------|-----------|-------------|
| **DockQ_AbAg** | AB + AC averaged | Primary docking quality metric (0-1). Geometric mean of fnat, iRMSD, LRMSD components. |
| DockQ_AB | AB only | Heavy chain - antigen interface quality |
| DockQ_AC | AC only | Light chain - antigen interface quality |
| DockQ_BC | BC only | VH-VL pairing quality (should remain high) |
| iRMSD_AB | AB | Interface RMSD — measures local interface geometry |
| LRMSD_AB | AB | Ligand RMSD — measures global antibody placement |
| fnat_AB | AB | Fraction of native contacts recovered |

### 5.2 CAPRI Classification

Based on DockQ_AbAg:

| Class | DockQ Range | Interpretation |
|-------|-------------|----------------|
| Incorrect | < 0.23 | Failed docking |
| Acceptable | 0.23 - 0.49 | Usable model |
| Medium | 0.49 - 0.80 | Good model |
| High | > 0.80 | Near-native |

### 5.3 CDR-H3 RMSD (New for N1)

**This is the key new metric** that validates Strategy Y's contribution.

- **Definition**: RMSD of CDR-H3 backbone atoms (CA) after superimposing the antibody framework region.
- **Alignment**: Superimpose on heavy chain framework residues (not CDR3), then compute RMSD on CDR-H3 residues only.
- **CDR-H3 boundaries**: From `examples/cdrs.csv`, column `cdr3_h` — complex-specific 1-indexed residue lists.
- **Implementation**:
  1. Parse CDR-H3 residue indices from `cdrs.csv` for each complex.
  2. Load prediction and ground truth structures.
  3. Extract heavy chain (B) CA atoms.
  4. Define framework = all heavy chain residues NOT in CDR-H3.
  5. Superimpose prediction on ground truth using framework CA atoms (Kabsch alignment).
  6. Compute RMSD over CDR-H3 CA atoms only.
- **Conda**: Requires BioPython for structure alignment — use `conda run -n boltz` (has BioPython).

### 5.4 Confidence Metrics

| Metric | Source | Description |
|--------|--------|-------------|
| confidence_score | `confidence_*.json` → `confidence_score` | Overall model confidence |
| iptm | `confidence_*.json` → `iptm` | Interface pTM — predicted quality of chain-chain interactions |
| ptm | `confidence_*.json` → `ptm` | Predicted TM-score |

### 5.5 Composite Score (if FK energy available)

If N1 predictions save FK energy (from the late-phase particle resampling), compute:
```
composite = alpha * confidence_score + (1 - alpha) * normalized_fk_energy
```
Check for `fk_energy` in the output. If unavailable, skip composite model selection.

---

## 6. Preprocessing and Alignment

### 6.1 Ground Truth

- Location: `pdb_minimized/{CID}.pdb` — 48 minimized crystal structures.
- Format: PDB, chains A (antigen), B (heavy), C (light).
- No preprocessing needed — these are the reference structures.

### 6.2 Predictions

DockQ v2 handles alignment internally when given model + native as input. No external alignment or preprocessing is needed for DockQ metrics.

For **CDR-H3 RMSD**, the alignment must be done explicitly:
1. Convert CIF to PDB if needed (DockQ handles both, but BioPython's PDB parser is more reliable for PDB format).
2. Superimpose using heavy chain framework CA atoms.
3. Compute CDR-H3 RMSD.

### 6.3 DockQ Execution

```bash
conda run -n dockq2 DockQ <model> <native> --mapping AB:AB AC:AC BC:BC --json
```

- DockQ v2 requires explicit chain mapping.
- JSON output provides per-interface metrics.
- Cache results in `evaluation_results/dockq_cache/{method}_{complex}_{model}.json` to avoid re-computation.

---

## 7. Model Selection Strategies

For methods with multiple models per complex:

| Strategy | How | Purpose |
|----------|-----|---------|
| **Confidence** | Pick model with highest `confidence_score` | Practical — no ground truth needed |
| **Oracle** | Pick model with highest DockQ_AbAg | Upper bound — requires ground truth |
| **Composite** | Pick by composite score (if FK energy available) | FK-aware re-ranking |

For B1 and B2 with multiple settings per complex:
- Pool all models from all settings.
- Confidence/oracle selection picks from the entire pool.
- This is the same approach as in the FK v2 evaluate.py.

---

## 8. Statistical Tests

### 8.1 Paired Wilcoxon Signed-Rank Test

For each metric (DockQ_AbAg, CDR-H3 RMSD), test N1 vs each baseline:
- H0: No difference between N1 and baseline.
- Only include complexes where both methods have predictions.
- Report: W statistic, p-value, effect size (mean difference).
- Threshold: p < 0.05 for significance.

### 8.2 Spearman Rank Correlation

- Confidence vs DockQ: measures calibration quality.
- FK energy vs DockQ (if available): validates composite scoring.

### 8.3 Per-Complex Delta Analysis

- Compute delta = N1_DockQ - B0_DockQ for each complex.
- Count wins/losses/ties.
- Identify complexes where N1 regresses vs B0 — investigate patterns.

---

## 9. Plots

| Plot | Description | Key Insight |
|------|-------------|-------------|
| Aggregate DockQ bars | Mean DockQ_AbAg per method (confidence + oracle) | Overall comparison |
| CAPRI stacked bars | Distribution of quality classes per method | Shift in quality distribution |
| DockQ boxplot | Distribution shape per method | Spread and outliers |
| N1 vs B0 scatter | Per-complex DockQ scatter (diagonal = no difference) | Which complexes improve/regress |
| Delta waterfall | Per-complex DockQ delta sorted | Clear win/loss visualization |
| CDR-H3 RMSD bars | Mean CDR-H3 RMSD per method | Validates Y's contribution |
| CDR-H3 RMSD scatter | N1 vs B0 per-complex CDR-H3 RMSD | Which complexes improve structurally |
| Confidence calibration | Confidence vs DockQ scatter per method | Calibration quality |
| Per-complex heatmap | DockQ matrix (complexes x methods) | Comprehensive per-complex view |
| FK energy vs DockQ | If available — scatter for N1 | Validates composite scoring |

---

## 10. Evaluation Script Architecture

### 10.1 File: `scripts/evaluate_n1.py`

Adapted from the FK v2 evaluate.py with these changes:

```
1. Constants
   - METHOD_SUBDIRS: B0, B1, B2, N1
   - METHOD_LABELS: updated for N1
   - CDR_CSV_PATH: path to examples/cdrs.csv

2. Discovery Functions
   - discover_b0(): unchanged
   - discover_b1(): unchanged
   - discover_b2(): unchanged
   - discover_n1(): new — pattern: boltz_results_{CID}/ with CIF or PDB models

3. DockQ Evaluation (reuse)
   - run_dockq(): unchanged
   - evaluate_single_model(): unchanged
   - select_best_model(): unchanged (confidence/oracle/composite)

4. CDR-H3 RMSD (new)
   - load_cdrs(): parse cdrs.csv into dict[complex_id -> cdr3_h_indices]
   - compute_cdr_h3_rmsd(): BioPython-based alignment + RMSD
   - Must run under conda boltz (BioPython), NOT dockq2

5. Method Evaluation
   - evaluate_method(): add CDR-H3 RMSD to per-model metrics

6. Aggregation (reuse + extend)
   - build_summary_table(): add CDR-H3 RMSD columns
   - build_aggregate_table(): add CDR-H3 RMSD mean/median/std
   - build_capri_table(): unchanged
   - run_statistical_tests(): add CDR-H3 RMSD Wilcoxon test

7. Plots (reuse + extend)
   - All existing plots: unchanged
   - plot_cdr_h3_rmsd_bars(): new
   - plot_cdr_h3_rmsd_scatter(): new (N1 vs B0)

8. Main
   - CLI args: predictions-dir, ground-truth-dir, output-dir, methods, complexes
   - Two-phase execution for conda environments (see Section 11)
```

### 10.2 CDR-H3 RMSD Implementation Detail

```python
def compute_cdr_h3_rmsd(pred_path, gt_path, cdr3_h_indices):
    """
    Compute CDR-H3 RMSD after framework superposition.

    Args:
        pred_path: Path to predicted structure (PDB or CIF)
        gt_path: Path to ground truth PDB
        cdr3_h_indices: List of 0-indexed residue indices for CDR-H3 in chain B

    Returns:
        float: CDR-H3 CA RMSD in Angstroms
    """
    # 1. Load structures with BioPython
    # 2. Extract chain B CA atoms
    # 3. Split into framework (non-CDR3) and CDR3 sets
    # 4. Superimpose using framework CAs (Kabsch)
    # 5. Apply same transform to CDR3 CAs
    # 6. Compute RMSD on CDR3 CAs
```

---

## 11. Conda Environment Strategy

Two environments are needed:

| Environment | Used For | Key Packages |
|-------------|----------|-------------|
| `dockq2` | DockQ v2 computation | DockQ, BioPython |
| `boltz` | CDR-H3 RMSD computation, CIF parsing | BioPython, numpy |

### Option A: Two-Phase Script (Recommended)

Run the script twice with different conda environments:

```bash
# Phase 1: DockQ computation (creates cache)
conda run -n dockq2 python scripts/evaluate_n1.py --phase dockq

# Phase 2: CDR-H3 RMSD + aggregation + plots
conda run -n boltz python scripts/evaluate_n1.py --phase analysis
```

Phase 1 computes DockQ for all models and caches results as JSON files.
Phase 2 reads cached DockQ, computes CDR-H3 RMSD, and generates tables + plots.

### Option B: Single Script with Subprocess

Run everything under `dockq2`, but shell out to `boltz` for CDR-H3 RMSD:
```python
result = subprocess.run(
    ["conda", "run", "-n", "boltz", "python", "-c", cdr_rmsd_script],
    capture_output=True
)
```

**Recommendation**: Option A is cleaner and easier to debug.

---

## 12. Output Structure

```
evaluation_results_n1/
  # CSV Tables
  per_complex_confidence.csv       # Per-complex metrics, confidence selection
  per_complex_oracle.csv           # Per-complex metrics, oracle selection
  aggregate_confidence.csv         # Aggregate summary, confidence selection
  aggregate_oracle.csv             # Aggregate summary, oracle selection
  capri_confidence.csv             # CAPRI distribution, confidence selection
  capri_oracle.csv                 # CAPRI distribution, oracle selection
  statistical_tests_confidence.csv # Wilcoxon test results
  statistical_tests_oracle.csv
  confidence_correlation.csv       # Spearman rho per method
  cdr_h3_rmsd_per_complex.csv     # CDR-H3 RMSD per complex per method

  # Plots
  plot_aggregate_dockq.png
  plot_capri_confidence.png
  plot_capri_oracle.png
  plot_boxplot_confidence.png
  plot_boxplot_oracle.png
  plot_scatter_n1_vs_b0.png
  plot_delta_n1_vs_b0.png
  plot_confidence_calibration.png
  plot_heatmap_dockq.png
  plot_cdr_h3_rmsd_bars.png
  plot_cdr_h3_rmsd_scatter_n1_vs_b0.png
  plot_irmsd.png
  plot_fnat.png

  # Cache
  dockq_cache/
    {method}_{complex}_{model}.json
```

---

## 13. Execution Plan

### Step 1: Generate Predictions (GPU required)
```bash
# Run N1 inference for all 48 complexes (on GPU machine)
for yml in yaml_hybrid_fk_hierarchical/*.yml; do
  conda run -n boltz boltz predict "$yml" \
    --use_potentials \
    --hybrid_fk_hierarchical \
    --hybrid_late_particles 5 \
    --num_samples 5 \
    --out_dir predictions_examples/new_feature_hybrid_fk/
done
```

### Step 2: Create Evaluation Script
Copy and adapt `../a_v2_fk_particles_true/scripts/evaluate.py` → `scripts/evaluate_n1.py`.

### Step 3: Run DockQ Phase
```bash
conda run -n dockq2 python scripts/evaluate_n1.py \
  --predictions-dir predictions_examples/ \
  --ground-truth-dir pdb_minimized/ \
  --output-dir evaluation_results_n1/ \
  --phase dockq
```

### Step 4: Run Analysis Phase
```bash
conda run -n boltz python scripts/evaluate_n1.py \
  --predictions-dir predictions_examples/ \
  --ground-truth-dir pdb_minimized/ \
  --output-dir evaluation_results_n1/ \
  --phase analysis
```

### Step 5: Review Results
- Check `aggregate_confidence.csv` for primary comparison.
- Check `statistical_tests_confidence.csv` for significance.
- Check `cdr_h3_rmsd_per_complex.csv` for Y's contribution.
- Review plots for patterns and regressions.

---

## 14. Success Criteria

| Criterion | Target | Measure |
|-----------|--------|---------|
| DockQ improvement over B0 | +0.05-0.08 | Mean DockQ_AbAg (confidence selection) |
| CDR-H3 RMSD improvement | -0.3-0.5A | Mean CDR-H3 RMSD (confidence selection) |
| Statistical significance | p < 0.05 | Paired Wilcoxon on DockQ_AbAg |
| No regression on easy cases | < 5 complexes regress | Per-complex delta waterfall |
| VH-VL pairing preserved | DockQ_BC >= B0 | Mean DockQ_BC |
| Competitive with B1 | Within 0.05 DockQ | Mean DockQ_AbAg gap to B1 |

---

## 15. Risk Assessment

| Risk | Impact | Mitigation |
|------|--------|------------|
| N1 only has 1 model per complex (test run) | No model selection possible | Run with `--num_samples 5` for evaluation |
| CIF vs PDB format mismatch | DockQ fails or chain mapping errors | Discovery function handles both; DockQ v2 supports CIF |
| CDR-H3 residue indexing off-by-one | Wrong RMSD computation | Validate against cdrs.csv; compare manual vs computed on 1-2 known structures |
| FK energy not saved to disk | Cannot do composite model selection | Fall back to confidence-only selection |
| OOM on large complexes | Missing predictions | Document which complexes fail; report N on all metrics |
| B1/B2 have many more models per complex | Unfair oracle comparison | Report both confidence and oracle; note model count per method |
