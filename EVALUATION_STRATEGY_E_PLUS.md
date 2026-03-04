# Evaluation Strategy: E+ (Enhanced Blind Scanning)

## 1. Feature Goals Recap

E+ (Enhanced Blind Scanning) has **two distinct goals**, and the evaluation must address both:

| Goal | Description |
|------|-------------|
| **Primary: Epitope Discovery** | Identify which antigen residues form the binding interface (epitope). The feature produces a per-residue epitope propensity heatmap via the `epitope_heatmap_*.npz` file. |
| **Secondary: Structure Quality** | The feature also produces standard 3D structure predictions (PDB/CIF). We evaluate whether region-specific β-scaling degrades or improves docking quality compared to baselines. |

This is fundamentally different from other steering methods (which primarily aim to improve docking). E+ is an **epitope prediction method** that happens to use structure prediction as its engine.

---

## 2. Baselines

| Label | Directory | Description |
|-------|-----------|-------------|
| B1 | `antigen_cut` | Vanilla Boltz2 prediction (no constraints) |
| B2 | `antigen_cut_contact_restraints` | Boltz2 with contact restraint guidance (multiple configs per complex) |
| B3 | `antigen_cut_vhvl_msa_pocket_ab_boltz_post_2023_omm` | Boltz2 with pocket+antibody constraints (multiple configs per complex) |
| **E+** | `new_feature_blind_scanning` (name TBD) | Enhanced Blind Scanning |

**Note on B2/B3:** These baselines have multiple restraint configurations per complex (e.g., `hbond_23`, `hydrophobic_5`, etc.). For fair comparison, we evaluate each configuration separately and also report the **best-per-complex** result (the config with highest confidence_score), since B2/B3 require prior knowledge of contact residues while E+ does not.

---

## 3. Model Selection Strategy

Each prediction run produces 5 models (`_model_0` through `_model_4`) with associated confidence files. We use two strategies:

1. **Best model (primary):** Select the model with the highest `confidence_score` from its JSON file. This is the standard approach for Boltz2 and what downstream users would do.
2. **All models (secondary):** For epitope heatmap evaluation, the heatmap is produced once per run (not per model), so model selection applies only to Axis B (structure quality).

---

## 4. Evaluation Axes

### Axis A: Epitope Prediction Accuracy (Primary)

This is the core evaluation for E+. We compare the predicted epitope (from the heatmap) against the ground truth epitope (from crystal structures).

#### 4A.1 Ground Truth Epitope Definition

From each crystal structure in `pdb_minimized/`:
- Extract antigen residues (chain A) whose **heavy atoms** are within a distance cutoff of any antibody **heavy atom** (chains B, C).
- **Default cutoff: 5.0 Å** (standard epitope definition). We also report results at 4.5 Å and 6.0 Å for sensitivity analysis.
- This produces a binary label per antigen residue: epitope (1) or non-epitope (0).

#### 4A.2 E+ Heatmap Evaluation (Continuous)

The E+ output `epitope_heatmap_*.npz` contains:
- `epitope_propensity`: float array in [0, 1] per antigen residue (fraction of region scans that made contact)
- `antigen_token_indices`: mapping from propensity index to token index
- `per_region_contacts`: binary contacts from each of the N region scans

**Metrics (continuous predictions):**

| Metric | Description | Why |
|--------|-------------|-----|
| **AUC-ROC** | Area under ROC curve | Threshold-independent measure of discrimination ability |
| **AUC-PR** | Area under Precision-Recall curve | Better than AUC-ROC for imbalanced data (epitopes are typically 15-30% of antigen) |
| **Spearman correlation** | Rank correlation between propensity and min-distance-to-antibody | Measures if higher propensity correlates with closer proximity |

**Metrics (binarized predictions at optimal threshold):**

| Metric | Description |
|--------|-------------|
| **Precision** | TP / (TP + FP) — fraction of predicted epitope residues that are true |
| **Recall** | TP / (TP + FN) — fraction of true epitope residues recovered |
| **F1** | Harmonic mean of precision and recall |
| **MCC** | Matthews Correlation Coefficient (accounts for class imbalance) |

The optimal threshold is determined by maximizing F1 on the validation set. We also report metrics at fixed thresholds (0.3, 0.5, 0.7).

#### 4A.3 Baseline Epitope Evaluation (Structure-Derived)

For baselines B1/B2/B3 (which don't produce heatmaps), the "predicted epitope" is extracted from the predicted structure using the same distance-based method as the ground truth:
- Extract antigen residues within cutoff of any antibody atom in the **predicted structure**.
- This is a binary prediction, so we only compute Precision, Recall, F1, MCC.

This comparison is fair because it answers: "Does the predicted structure correctly capture which antigen residues are at the interface?"

#### 4A.4 Heatmap-Specific Visualizations

For each complex:
- Plot propensity heatmap alongside ground truth epitope labels (bar chart with colored GT overlay)
- Plot ROC and PR curves

Aggregate:
- Mean AUC-ROC and AUC-PR across all complexes (with confidence intervals)

---

### Axis B: Docking / Interface Quality (Secondary)

Even though E+ is primarily an epitope discovery tool, it also produces 3D structures. We evaluate whether the β-scaling distorts the structure or incidentally improves docking.

#### 4B.1 DockQ Metrics

Computed using the `DockQ` tool (`conda activate dockq2`):

| Metric | Description |
|--------|-------------|
| **DockQ** | Combined docking quality score (0-1). Thresholds: <0.23 Incorrect, 0.23-0.49 Acceptable, 0.49-0.80 Medium, ≥0.80 High |
| **iRMSD** | Interface RMSD (Å) — RMSD of interfacial residues |
| **LRMSD** | Ligand RMSD (Å) — RMSD of the ligand after receptor alignment |
| **fnat** | Fraction of native contacts recovered |
| **fnonnat** | Fraction of predicted contacts that are non-native |
| **F1** | DockQ's contact F1 (harmonic mean of precision and recall of contacts) |

**Interface of interest:** We focus on the **antibody-antigen interfaces** (AB and AC), not the intra-antibody interface (BC). Report:
- `DockQ_AB_AC`: Average DockQ across AB and AC interfaces
- `DockQ_global`: Average across all interfaces
- Per-interface breakdown

#### 4B.2 Antibody-Aligned Antigen RMSD

A more interpretable metric for antibody-antigen prediction:
1. Superimpose predicted structure onto ground truth using antibody Cα atoms (chains B+C) via Kabsch algorithm
2. Compute RMSD of antigen Cα atoms (chain A) after alignment

This directly measures: "Given the antibody is in the right place, how far off is the predicted antigen position?"

#### 4B.3 Confidence Metrics

From Boltz2 confidence JSONs:

| Metric | Description |
|--------|-------------|
| `confidence_score` | Overall confidence (4×pLDDT + iPTM or PTM) / 5 |
| `iptm` | Interface predicted TM-score |
| `ptm` | Predicted TM-score |
| `complex_plddt` | Complex-level pLDDT |
| `complex_iplddt` | Interface-level pLDDT |

---

### Axis C: Efficiency (Tertiary)

The key claim of E+ is a 10-20x speedup over running N separate full predictions. We measure:

| Metric | How |
|--------|-----|
| **Wall-clock time per complex** | Parse from logs or time the run. E+ trunk is computed once + N lightweight region scans. |
| **Speedup ratio** | `time(N separate B1 runs) / time(E+ with N regions)` |
| **GPU memory** | Peak GPU memory during E+ vs. N separate runs |

This is measured empirically on the evaluation server. For the document, we note the expected speedup and validate it.

---

## 5. Postprocessing Pipeline

### 5.1 Structure Alignment (for Axis B)

Before computing RMSD metrics, structures need alignment. The approach:

1. **Parse structures:** Use BioPython `PDBParser`/`MMCIFParser` (both PDB and CIF formats must be supported since different methods may output different formats).
2. **Extract Cα atoms** per chain, matched by sequential position (not residue number, since GT and predictions may differ in numbering).
3. **Kabsch alignment:** Align on antibody Cα atoms (chains B+C), then measure antigen Cα RMSD.

No additional preprocessing (e.g., OpenMM minimization) is needed since we compare raw predictions against already-minimized ground truths.

### 5.2 Epitope Extraction (for Axis A)

1. **From crystal structure (GT):** Heavy-atom distance cutoff between chains A and B+C.
2. **From predicted structure (baselines):** Same method applied to predicted coordinates.
3. **From E+ heatmap:** Direct from `epitope_propensity` array in `.npz` file.
4. **Residue alignment:** Map antigen residues by sequential position (not PDB residue number) to ensure GT and prediction indices match.

### 5.3 Handling Multiple Configs (B2/B3)

B2 and B3 have multiple restraint configurations per complex. For aggregation:
- **Per-config:** Evaluate each config independently (e.g., `hbond_23`, `salt_bridge_3`).
- **Best-per-complex:** Select the config with the highest `confidence_score` and report that as the method's result for the complex.
- **Oracle:** Select the config with the best DockQ (upper bound on method performance).

---

## 6. Environment and Tools

| Environment | Usage |
|-------------|-------|
| `conda activate boltz` | Running E+ predictions, loading Boltz2 outputs, BioPython parsing, NumPy/SciPy analysis |
| `conda activate dockq2` | Running DockQ (called via `subprocess` from the evaluation script with `conda run -n dockq2 DockQ ...`) |

**Required packages (all available in `boltz` env):**
- `biopython` — Structure parsing, Cα extraction
- `numpy` — Array operations, RMSD computation
- `scipy` — `cKDTree` for distance queries, `scipy.stats.spearmanr` for correlation
- `scikit-learn` — `roc_auc_score`, `average_precision_score`, `precision_recall_curve` for AUC metrics
- `pandas` — Results aggregation
- `matplotlib` — Plots

If `scikit-learn` is not in the `boltz` env, we can install it (`pip install scikit-learn`) or implement AUC-ROC/AUC-PR manually using NumPy (straightforward for sorted arrays). No new conda environment is needed.

---

## 7. Output Structure

```
evaluation_output_eplus/
├── results_per_complex.csv          # All metrics per method × complex × config
├── summary_by_method.csv            # Aggregated means/medians/stds per method
├── dockq_cache/                     # Cached DockQ JSON results (avoid re-running)
│   └── {method}_{complex}_{config}_{model}.json
├── plots/
│   ├── auc_roc_bar.png              # Mean AUC-ROC per method
│   ├── auc_pr_bar.png               # Mean AUC-PR per method
│   ├── epitope_f1_bar.png           # Mean epitope F1 per method
│   ├── dockq_bar.png                # Mean DockQ per method
│   ├── dockq_distribution.png       # DockQ quality category distribution
│   ├── eplus_vs_b1_dockq.png        # E+ vs B1 scatter
│   ├── heatmap_{complex}.png        # Per-complex propensity heatmaps (E+ only)
│   └── roc_pr_{complex}.png         # Per-complex ROC and PR curves (E+ only)
└── tables/
    ├── table1_epitope_prediction.md  # Axis A summary
    ├── table2_docking_quality.md     # Axis B summary
    └── table3_confidence.md          # Confidence metrics
```

---

## 8. Summary Tables (Expected Output Format)

### Table 1: Epitope Prediction (Axis A)

| Method | AUC-ROC | AUC-PR | F1 (optimal) | MCC | Precision | Recall |
|--------|---------|--------|--------------|-----|-----------|--------|
| B1     | —       | —      | x.xxx        | x.xxx | x.xxx   | x.xxx  |
| B2 (best) | —    | —      | x.xxx        | x.xxx | x.xxx   | x.xxx  |
| B3 (best) | —    | —      | x.xxx        | x.xxx | x.xxx   | x.xxx  |
| **E+** | **x.xxx** | **x.xxx** | **x.xxx** | **x.xxx** | **x.xxx** | **x.xxx** |

(AUC-ROC/AUC-PR only available for E+ since baselines produce binary predictions)

### Table 2: Docking Quality (Axis B)

| Method | DockQ (AB+AC) | DockQ (global) | iRMSD (AB+AC) | fnat (AB+AC) | Ab-aligned Ag RMSD |
|--------|---------------|----------------|---------------|--------------|---------------------|
| B1     | x.xxx         | x.xxx          | x.xxx         | x.xxx        | x.xxx               |
| B2 (best) | x.xxx      | x.xxx          | x.xxx         | x.xxx        | x.xxx               |
| B3 (best) | x.xxx      | x.xxx          | x.xxx         | x.xxx        | x.xxx               |
| **E+** | x.xxx         | x.xxx          | x.xxx         | x.xxx        | x.xxx               |

### Table 3: Confidence Metrics

| Method | confidence_score | iPTM | complex_pLDDT | complex_ipLDDT |
|--------|-----------------|------|---------------|----------------|
| B1     | x.xxxx          | x.xxxx | x.xxxx      | x.xxxx         |
| ...    | ...             | ...    | ...          | ...            |

---

## 9. Implementation Plan

The evaluation script will be `evaluate_e_plus.py`, adapted from the existing `evaluate_k_plus.py` (from the K+ feature). Key additions:

1. **E+ folder discovery:** Add `_extract_complex_name_eplus()` parser for E+ result folders.
2. **Heatmap loading:** Load `epitope_heatmap_*.npz` from E+ predictions folder.
3. **AUC metrics:** Add `compute_heatmap_auc_metrics()` using scikit-learn or manual implementation.
4. **Correlation:** Add Spearman correlation between propensity and min-distance-to-antibody.
5. **Heatmap visualization:** Per-complex propensity bar charts with GT overlay.
6. **Best-per-complex aggregation:** For B2/B3 multi-config methods.

### CLI Usage

```bash
# Full evaluation (requires dockq2 env for DockQ, boltz env for everything else)
conda run -n boltz python evaluate_e_plus.py \
    --predictions_dir predictions_examples \
    --gt_dir pdb_minimized \
    --output_dir evaluation_output_eplus

# Skip DockQ (use cached or skip Axis B)
conda run -n boltz python evaluate_e_plus.py \
    --predictions_dir predictions_examples \
    --gt_dir pdb_minimized \
    --skip_dockq

# Verbose with custom epitope cutoff
conda run -n boltz python evaluate_e_plus.py \
    --predictions_dir predictions_examples \
    --gt_dir pdb_minimized \
    --cutoff 4.5 \
    --verbose
```

---

## 10. Key Evaluation Questions

The evaluation should answer these questions:

1. **Does E+ correctly identify epitope regions?** (AUC-ROC > 0.7, AUC-PR > 0.5 target)
2. **Does E+ outperform structure-derived epitope extraction from baselines?** (F1 comparison)
3. **Does the region-specific β-scaling degrade structure quality?** (DockQ comparison with B1)
4. **How does E+ compare to constraint-based methods (B2/B3)?** (B2/B3 use oracle contact info; E+ discovers contacts blindly — a lower DockQ is expected, but better epitope discovery from fewer assumptions)
5. **Is the heatmap confidence-calibrated?** (Do higher propensity values correlate with actual epitope residues?)
6. **Is E+ efficient?** (Wall-clock time vs. N×B1)
