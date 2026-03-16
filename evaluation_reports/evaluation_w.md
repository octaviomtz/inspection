# Evaluation Report: Strategy W — Embedding-Based Interface Steering

**Date**: 2026-03-13
**Feature**: `EmbeddingInterfacePotential`
**Evaluation script**: `evaluate_embed_interface.py`
**Dataset**: 47 antibody-antigen complexes with ground truth crystal structures

---

## 1. Original Goal

Strategy W hypothesizes that Boltz2's pair representation `z` (shape `[batch, N_tokens, N_tokens, token_z]`) encodes which CDR-antigen residue pairs are likely to form productive binding interactions. By extracting CDR-antigen pairwise embeddings and using their magnitudes to weight distance-based steering, the method aims to capture **binding specificity** rather than treating all CDR-antigen pairs equally.

The expected improvements were:
1. Better antibody-antigen interface quality (DockQ)
2. More accurate contact/epitope prediction
3. Comparable or better overall structure quality

---

## 2. Experimental Setup

### 2.1 Conditions

| ID | Condition | Description |
|----|-----------|-------------|
| **B1** | No steering | Boltz2 default, no potentials |
| **B2** | Contact restraints | Coordinate-space contact restraints (hbond, hydrophobic, salt_bridge variants) |
| **B3** | Pocket + MSA | Pocket constraints targeting specific residues |
| **W** | Embed. interface | Embedding-based interface steering (this feature) |

Two evaluation runs were performed:
- **w_plus_2/**: W vs B1 only (2 conditions)
- **w_plus_4/**: W vs B1, B2, and B3 (4 conditions)

Both runs produced identical W and B1 results (same predictions evaluated).

### 2.2 Model and Variant Selection

- **Model selection**: For DockQ plots, the best model out of 5 was selected per complex by **best DockQ_antigen_avg** (oracle, i.e. selecting the model that happened to score best). For epitope and SASA metrics, the best model was selected by **best complex_iplddt** (model confidence).
- **Variant aggregation** (B2/B3 only): These baselines have multiple sub-folders per complex (different restraint types or target residues). For plotting, the **best variant** per complex was used (highest DockQ_antigen_avg for DockQ; highest epitope_f1 for epitope).

---

## 3. Results

### 3.1 DockQ Interface Quality

**Aggregate statistics (best model, best variant for B2/B3):**

| Condition | DockQ_AB mean | DockQ_AB median | DockQ_AC mean | DockQ_AC median | DockQ_avg mean | DockQ_avg median |
|-----------|--------------|----------------|--------------|----------------|---------------|-----------------|
| B1 (no steering) | 0.207 | 0.046 | 0.200 | 0.029 | 0.204 | 0.042 |
| B2 (contact restr.) | 0.410 | 0.451 | 0.402 | 0.484 | 0.406 | 0.482 |
| B3 (pocket+MSA) | 0.205 | 0.048 | 0.191 | 0.024 | 0.198 | 0.037 |
| **W (embed. interface)** | **0.213** | **0.045** | **0.202** | **0.032** | **0.207** | **0.043** |

**Paired comparison W vs baselines (DockQ_antigen_avg):**

| Comparison | Mean delta | Median delta | Wins | Ties | Losses | Wilcoxon p |
|------------|-----------|-------------|------|------|--------|-----------|
| W vs B1 | +0.004 | -0.000 | 5 | 31 | 11 | 0.224 |
| W vs B2 (best) | -0.194 | -0.071 | 0 | 9 | 37 | <0.001 |
| W vs B3 (best) | +0.009 | -0.002 | 4 | 24 | 8 | 0.074 |

---

### 3.2 Plot-by-Plot Analysis

#### Figure 1: DockQ Boxplot Comparison

**File**: `evaluation_plots/dockq_boxplot_comparison.png`
**Method**: Best model per complex selected by highest DockQ_antigen_avg (oracle). For B2/B3, best variant per complex. Boxplots of DockQ scores across all 47 complexes.

**Subplot 1 — DockQ: Heavy-Antigen (AB)**: Distributions for B1 and W are virtually identical (medians both near 0.10, IQR 0.01-0.57). B2 contact restraints show a clearly elevated distribution (median ~0.45, upper quartile ~0.70). B3 is comparable to B1/W.

**Subplot 2 — DockQ: Light-Antigen (AC)**: Same pattern. B1 and W overlap almost entirely (medians ~0.09). B2 is substantially better (median ~0.53). B3 is comparable to B1/W.

**Subplot 3 — DockQ: Antigen avg (AB+AC)/2**: Confirms the overall picture. W provides no meaningful separation from B1. B2 is the clear winner. The horizontal dashed lines at 0.23, 0.49, 0.80 indicate CAPRI thresholds; both B1 and W have medians well below the "Acceptable" threshold (0.23).

**Conclusion**: Strategy W does not improve DockQ distributions relative to no steering. Contact restraints (B2) substantially outperform all other methods.

---

#### Figure 2: Per-Complex DockQ Improvement (W vs B1)

**File**: `evaluation_plots/delta_dockq_W_vs_B1.png`
**Method**: For each complex, delta = DockQ_antigen_avg(W) - DockQ_antigen_avg(B1), using best model by DockQ. Sorted horizontal bar chart; green = improvement, red = degradation.

**Analysis**: The plot is nearly symmetric around zero. A few complexes show large improvements (8OL9_BAH: +0.49, 8IDN_HLA: +0.38, 8EZ7_HLA: +0.30, 7TRI_ZYB: +0.19, 7TRH_HBG: +0.14), but these are counterbalanced by comparable degradations (7Y0O_HLA: -0.22, 8IX3_HLG: -0.16, 8HES_HLC: -0.12, 8EZ8_HLA: -0.11). The majority of complexes (~31/47) show negligible change (|delta| <= 0.02).

**Conclusion**: W produces sporadic per-complex improvements that are offset by similar degradations elsewhere. There is no systematic benefit — the wins appear stochastic rather than driven by the embedding signal.

---

#### Figure 3: CAPRI Quality Categories

**File**: `evaluation_plots/capri_categories.png`
**Method**: Each complex classified by DockQ_antigen_avg into Incorrect (<0.23), Acceptable (0.23-0.49), Medium (0.49-0.80), High (>=0.80). Stacked bar chart as percentage of complexes.

**w_plus_2 (B1 vs W):**

| Category | B1 | W |
|----------|-----|-----|
| Incorrect | 68% | 57% |
| Acceptable | 2% | 6% |
| Medium | 28% | 34% |
| High | 2% | 2% |

**w_plus_4 (all 4 conditions):**

| Category | B1 | B2 | B3 | W |
|----------|-----|-----|-----|-----|
| Incorrect | 68% | 22% | 64% | 57% |
| Acceptable | 2% | 21% | 2% | 6% |
| Medium | 28% | 50% | 32% | 34% |
| High | 2% | 7% | 2% | 2% |

**Analysis**: W modestly reduces the "Incorrect" fraction from 68% to 57% and increases "Acceptable" and "Medium" categories. However, this shift is minor — the dominant category remains "Incorrect" for both B1 and W. In contrast, B2 dramatically reduces "Incorrect" to 22% and achieves 7% "High" quality predictions.

**Conclusion**: W produces a small shift toward better CAPRI categories, but the effect is weak. B2 is far more effective at moving complexes into correct docking quality categories.

---

#### Figure 4: Epitope Recovery

**File**: `evaluation_plots/epitope_boxplot_comparison.png`
**Method**: For each complex, ground truth and predicted epitopes defined as antigen residues with any heavy atom within 4.5 A of any antibody heavy atom. Best model selected by complex_iplddt. Boxplots of F1, Precision, Recall.

**Subplot 1 — Epitope F1**: B1 and W are nearly identical (median ~0.31 for B1, ~0.29 for W; mean ~0.37 for both). B2 contact restraints achieve substantially higher F1 (median ~0.72). B3 is similar to B1/W.

**Subplot 2 — Epitope Precision**: Same pattern. B1 and W overlap (median ~0.31 vs 0.30). B2 is clearly superior (median ~0.65).

**Subplot 3 — Epitope Recall**: B1 and W nearly identical (median ~0.31 vs 0.29). B2 again dominates (median ~0.66).

| Condition | F1 mean | F1 median | Precision mean | Recall mean |
|-----------|---------|-----------|----------------|-------------|
| B1 | 0.370 | 0.320 | 0.387 | 0.365 |
| B2 (best) | 0.605 | 0.722 | 0.630 | 0.601 |
| B3 (best) | 0.354 | 0.203 | 0.379 | 0.340 |
| **W** | **0.365** | **0.286** | **0.396** | **0.351** |

**Conclusion**: W does not improve epitope recovery over B1. The embedding weights do not produce more accurate contact predictions at the residue level.

---

#### Figure 5: SASA / Buried Surface Area

**File**: `evaluation_plots/sasa_boxplot_comparison.png`
**Method**: SASA computed using BioPython ShrakeRupley for each chain in isolation (unbound) and full complex (bound). BSA = SASA_unbound - SASA_bound. Ratios are predicted/native. Best model by complex_iplddt.

**Subplot 1 — BSA ratio (pred/native)**: All conditions have median BSA ratios below 1.0 (B1: 0.86, W: 0.89, B2: 0.97, B3: 0.91), indicating predictions generally bury less surface area than ground truth. W is slightly closer to 1.0 than B1 but the difference is minor. B2 is closest to the ideal ratio of 1.0.

**Subplot 2 — Buried residue ratio (pred/native)**: B1 median = 1.00, W median = 0.96. Both are close to ideal. B2 and B3 are also near 1.0. W shows slightly more outliers with high ratios (>2.0).

| Condition | BSA ratio mean | BSA ratio median | Buried count ratio mean | Buried count ratio median |
|-----------|---------------|-----------------|----------------------|------------------------|
| B1 | 0.895 | 0.856 | 1.056 | 1.000 |
| B2 | 0.910 | 0.848 | 1.048 | 1.000 |
| B3 | 0.881 | 0.846 | 1.029 | 1.000 |
| **W** | **0.873** | **0.893** | **1.041** | **0.962** |

**Conclusion**: All conditions produce similarly realistic interface sizes. W does not degrade or improve interface structural properties compared to B1.

---

## 4. Overall Conclusions

### 4.1 Did Strategy W Achieve Its Goals?

**No.** Strategy W failed to deliver meaningful improvements on any of the three evaluation axes:

1. **Interface quality (DockQ)**: No statistically significant improvement over B1 (p=0.22). Mean improvement of +0.004 DockQ is negligible. The method won on only 5/47 complexes while losing on 11.

2. **Epitope recovery**: Essentially identical to B1 (F1: 0.365 vs 0.370). The embedding weights did not improve which antigen residues were contacted.

3. **Structure quality**: Neither improved nor degraded — SASA ratios and buried residue counts are comparable across all conditions.

### 4.2 Context: How Do Other Methods Compare?

B2 (contact restraints) dramatically outperforms all other conditions, achieving:
- 2x higher mean DockQ_avg (0.406 vs 0.207)
- 78% of complexes at Acceptable or better (vs 32% for W)
- 1.6x higher epitope F1 (0.605 vs 0.365)

This confirms that steering antibody-antigen docking is achievable with explicit coordinate-space restraints. B3 (pocket+MSA) performs comparably to B1/W, suggesting that pocket-level constraints alone are also insufficient.

### 4.3 Why Did Strategy W Fail?

Several factors likely contributed:

1. **Pair embedding norms may not encode binding specificity**: The approach assumed that `||z[i,j]||` for CDR-antigen pairs correlates with binding propensity. In practice, pair embedding norms may primarily reflect sequence distance, chain identity, or other structural features rather than interface-specific interaction quality. The weights may be near-uniform across CDR-antigen pairs, providing little differentiation.

2. **Weak gradient signal**: Even if the embeddings contain useful information, the weighting may be too subtle to meaningfully redirect the diffusion trajectory. The flat-bottom distance potential multiplied by normalized [0,1] weights may produce gradients too similar to uniform weighting.

3. **Single-pass embedding extraction**: The embedding weights are extracted once from the trunk output `z` and held fixed across diffusion steps. As the structure evolves during sampling, the initially computed weights become stale.

### 4.4 Suggested Modifications

If pursuing this direction further, the following modifications could be explored:

1. **Validate the embedding signal**: Before modifying the steering pipeline, directly verify that pair embedding norms correlate with true contacts. Extract `z[i,j]` for all CDR-antigen pairs from the ground truth complexes and measure whether high-norm pairs correspond to real contacts (e.g., AUC-ROC). If this correlation is weak, the entire approach is fundamentally limited.

2. **Use a learned scoring head instead of raw norms**: Rather than using `||z[i,j]||` as weights, train a small MLP head on top of pair embeddings to predict contact probability. This could be trained separately on existing antibody-antigen complex data.

3. **Re-extract embeddings during diffusion**: Instead of computing weights once from the initial trunk pass, re-extract embeddings at each (or every N-th) diffusion step so the weights adapt as the structure evolves.

4. **Increase force magnitude**: The current embedding-weighted potential may be too weak relative to other forces. Experiment with larger force multipliers (e.g., 2x-5x) specifically for high-weight pairs to create a stronger directional signal.

5. **Combine with contact restraints**: Use embedding weights to prioritize which contact restraints to apply, rather than replacing them entirely. This hybrid approach could leverage the proven effectiveness of B2-style restraints while using embeddings for intelligent residue selection.

6. **Focus on CDR-H3 only**: CDR-H3 is typically the most critical loop for antigen binding. Restricting the embedding-weighted steering to CDR-H3 pairs only may produce a stronger, more focused signal.
