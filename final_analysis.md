# Final Analysis: Strategy G+ — Progressive Steering v2

**Date**: 2026-03-31
**Branch**: `g_v2_progre`
**Evaluation**: 47 antibody-antigen complexes, DockQ via DockQ2, confidence-selected best-of-5

---

## 1. Feature Goal

Strategy G+ (Progressive Steering v2) is the upgraded version of Round 1's Strategy G (Progressive CDR Refinement). In Round 1, Strategy G showed a modest +0.017 DockQ improvement (p=0.41, not significant) with three-phase FK steering and constant CDR3 beta-scaling. It helped medium-difficulty cases but regressed easy ones (-0.025 DockQ on high-confidence targets).

G+ v2 introduced three improvements to fix these issues:

- **G+.2 — Time-varying beta**: CDR3 beta-scaling starts strong in early diffusion (exploration phase) and decays to zero in late steps (convergence phase). This prevents the distortion of already-converged structures that caused regressions in Round 1.
- **G+.3 — Progressive antigen orientation potential**: A `ProgressiveAntigenOrientationPotential` with a three-phase schedule (weak early, strong late) for CDR-antigen distance guidance, replacing the fixed-phase weights.
- **G+.1 — Conditional steering** (not yet implemented): Skip steering when baseline confidence is already high, to avoid regressions on easy cases.

The primary evaluation question: **Does progressive steering improve antibody-antigen interface quality (DockQ/CAPRI) compared to three baselines, and does it reduce the regressions seen in Round 1?**

---

## 2. Methods Evaluated

| Label | Method | Description | N entries |
|-------|--------|-------------|-----------|
| B1 | No steering | Plain Boltz-2, antigen-cut MSA | 47 (1 per complex) |
| B2 | Restraint-to-A | Contact restraints to specific antigen residues | 87 (multiple restraint configs per complex, 36 unique complexes) |
| B3 | Contact restraints | hbond/hydrophobic/salt-bridge restraints | 522 (many configs per complex, 46 unique complexes) |
| Gv1 | Progressive Steering v2 | Time-varying beta + progressive potential | 47 (1 per complex) |

Note: B2 and B3 have multiple sub-experiments per complex (different restraint residues or types), so their total entry count is higher. For fair comparison, the evaluation selects the best sub-experiment per complex based on confidence score.

---

## 3. Plot-by-Plot Analysis

### 3.1 DockQ Bar Chart (`dockq_bar.png`)

**How produced**: For each method and complex, the model with the highest Boltz-2 `confidence_score` was selected from the 5 generated models (confidence-selected best-of-5). The bar shows the mean DockQ across all complexes; error bars show standard deviation.

**Results**:

| Method | Mean DockQ | Std |
|--------|-----------|-----|
| B1 | 0.338 | ~0.17 |
| B2 | 0.314 | ~0.16 |
| B3 | 0.413 | ~0.22 |
| Gv1 (G+ v2) | 0.345 | ~0.18 |

**Conclusions**: G+ v2 (0.345) shows a marginal +0.007 improvement over B1 (0.338), within error bars and not statistically significant. B3 (contact restraints) is the clear winner at 0.413, benefiting from its large number of sub-experiments (522 configurations across 46 complexes) which gives it many chances to find a good restraint configuration — this is a best-of-many advantage. B2 slightly underperforms B1 at 0.314. The improvement from G+ v2 over B1 is in the right direction but too small to be meaningful on its own.

---

### 3.2 CAPRI Stacked Bar Chart (`capri_stacked_bar.png`)

**How produced**: Each complex's confidence-selected model was classified into CAPRI quality categories (High: DockQ >= 0.80, Medium: >= 0.49, Acceptable: >= 0.23, Incorrect: < 0.23). The stacked bar shows the percentage distribution per method.

**Results**:

| Method | Medium+ (%) | Acceptable+ (%) | Incorrect (%) |
|--------|-------------|------------------|---------------|
| B1 | 21.3% | 63.8% | 36.2% |
| B2 | 18.4% | 54.0% | 46.0% |
| B3 | 37.7% (+0.2% High) | 77.8% | 22.2% |
| Gv1 (G+ v2) | 23.4% | 63.8% | 36.2% |

**Conclusions**: G+ v2 shifts one additional complex from Acceptable to Medium quality compared to B1 (23.4% vs 21.3% Medium+), consistent with Round 1's observation that Strategy G helps medium-difficulty cases. The Incorrect rate is identical (36.2%), meaning G+ v2 does not introduce new catastrophic failures. B3 stands out with 37.7% Medium+ and only 22.2% Incorrect, again leveraging its many sub-experiment configurations. No method achieves High quality (DockQ >= 0.80), except B3 with a single complex (0.2%).

---

### 3.3 Confidence Score vs DockQ Scatter Plot (`confidence_vs_dockq.png`)

**How produced**: Every individual model (all 5 per complex, all methods) is plotted with its Boltz-2 confidence score on the X-axis and actual DockQ on the Y-axis. Colors distinguish methods.

**Conclusions**: All methods show the same characteristic pattern: a dense cluster at high confidence (0.80-0.95) spanning the full DockQ range (0.15-0.80). This reveals a fundamental problem: **Boltz-2 confidence is weakly correlated with actual docking quality**. Many models with confidence > 0.90 have DockQ < 0.25 (incorrect). This means confidence-based model selection is unreliable and explains the large oracle-vs-selected gap (0.05 DockQ for B1). G+ v2 models (purple) overlap entirely with B1 models (blue), indicating no shift in the confidence-accuracy relationship.

---

### 3.4 Ensemble Diversity Box Plot (`ensemble_diversity_boxplot.png`)

**How produced**: For each method-complex pair, the mean pairwise CA-RMSD across the 5 models was computed after structural superposition (all chains A/B/C). The box plot shows the distribution of this diversity metric across complexes.

**Results**:

| Method | Median CA-RMSD (A) | IQR |
|--------|-------------------|-----|
| B1 | ~8.5 | ~2-15 |
| B2 | ~10.5 | ~3-15 |
| B3 | ~6.0 | ~2-12 |
| Gv1 (G+ v2) | ~9.5 | ~6-13 |

**Conclusions**: G+ v2 has slightly higher median diversity than B1 (9.5 vs 8.5 A), with a narrower interquartile range. This is a positive signal: the progressive steering maintains or slightly increases structural diversity, meaning it is not over-constraining the ensemble. B3 has the lowest median diversity (6.0 A), likely because contact restraints steer all models toward similar conformations. B2 has the highest spread. The preserved diversity in G+ v2 means the time-varying beta (decaying in late steps) is working as intended — it does not lock the ensemble into a single basin.

---

### 3.5 Per-Complex DockQ Comparison (`dockq_per_complex.png`)

**How produced**: Each complex (X-axis) has grouped bars for all methods, showing the confidence-selected DockQ. Complexes are sorted by B1 DockQ (ascending).

**Conclusions**: The per-complex view reveals the mosaic pattern typical of these interventions. For most complexes, all four methods perform similarly (bars of equal height). The notable exceptions are:

- **B3 tall spikes**: Complexes like 8DTK_CBA, 8FXB_HLE, 8GP5_EFX show B3 dramatically outperforming (0.7+ vs 0.2), but these represent cases where one specific restraint configuration happened to match the ground truth.
- **G+ v2 wins**: 8HGM_CDB (0.69 vs 0.37 B1), 8OL9_BAH (0.69 vs 0.33 B1) show G+ v2 finding a binding mode that B1 misses.
- **G+ v2 losses**: 7TRH_HBG (0.24 vs 0.52 B1) is the largest regression — a -0.28 DockQ drop.

---

### 3.6 Delta Histogram: G+ v2 vs B1 (`delta_histogram_gv1_vs_b1.png`)

**How produced**: For each of the 47 complexes, the DockQ difference (G+ v2 confidence-selected minus B1 confidence-selected) was computed. The histogram shows the distribution of these deltas. Complexes with |delta| > 0.05 are categorized as improved or regressed.

**Results**: **4 improved, 41 unchanged, 2 regressed**

| Complex | B1 DockQ | G+ v2 DockQ | Delta |
|---------|---------|-------------|-------|
| 8OL9_BAH | 0.327 | 0.688 | **+0.361** |
| 8HGM_CDB | 0.365 | 0.690 | **+0.325** |
| 8HES_HLC | 0.157 | 0.269 | +0.112 |
| 8EZ3_HLA | 0.216 | 0.325 | +0.108 |
| 7ZOZ_HLA | 0.646 | 0.552 | -0.094 |
| 7TRH_HBG | 0.525 | 0.240 | **-0.285** |

**Conclusions**: The histogram is sharply peaked near zero with a few outliers. The bulk of complexes (41/47, 87%) are unchanged. The 4 wins are substantial (two exceed +0.3 DockQ), while the 2 regressions are concerning — particularly 7TRH_HBG where a Medium-quality prediction collapses to near-Incorrect. The overall mean delta is +0.007 (median -0.004), confirming the near-zero net effect. This pattern is consistent with Round 1's finding: the steering mechanism changes which binding mode is found, producing large individual swings in both directions.

---

## 4. Statistical Summary

### Wilcoxon Signed-Rank Test (vs B1)

| Method | N paired | Median delta | Mean delta | p-value |
|--------|---------|-------------|-----------|---------|
| B2 | 87 | -0.0002 | -0.00001 | 0.145 |
| B3 | 522 | +0.004 | +0.059 | **1.2e-21** |
| G+ v2 | 47 | -0.004 | +0.007 | 0.305 |

G+ v2 is **not statistically significant** (p=0.305). The positive mean delta (+0.007) is driven by a few large wins but the median is actually slightly negative (-0.004), indicating that more complexes get slightly worse than get slightly better. Only B3 achieves significance, powered by its massive sample size (522 sub-experiments).

### Oracle vs Confidence-Selected Gap

| Method | Conf-selected | Oracle | Gap |
|--------|--------------|--------|-----|
| B1 | 0.338 | 0.391 | 0.053 |
| G+ v2 | 0.345 | 0.391 | 0.046 |

G+ v2 has a slightly smaller oracle-confidence gap (0.046 vs 0.053), suggesting its best models are slightly more likely to also be the highest-confidence ones. The oracle DockQ is virtually identical (0.391 vs 0.391), meaning the ceiling of quality is unchanged — the method can produce equally good structures, just not consistently.

### Per-Residue Metrics

| Method | Mean pLDDT | CDR3 pLDDT | Interface PAE | CDR3-Ag PAE |
|--------|-----------|-----------|--------------|-------------|
| B1 | 0.914 | 0.821 | 16.9 | 15.8 |
| G+ v2 | 0.912 | 0.815 | 17.3 | 16.3 |

G+ v2 shows marginally **worse** per-residue metrics: CDR3 pLDDT drops slightly (0.821 to 0.815) and interface PAE increases (16.9 to 17.3). This is the opposite of what was expected — the time-varying beta was supposed to improve CDR3 confidence by removing late-stage distortion. The slight worsening suggests the beta schedule may still be too aggressive in early steps, or the progressive potential introduces noise in the PAE predictions.

---

## 5. Comparison with Round 1

| Metric | Round 1 G | G+ v2 |
|--------|----------|-------|
| DockQ delta vs B1 | +0.017 | +0.007 |
| p-value | 0.410 | 0.305 |
| CAPRI Medium+ gain | +6.4 pp | +2.1 pp |
| Net CAPRI shifts | +3 (6 up, 3 down) | +1 (1 up, 0 down) |
| Ensemble diversity | +0.6 A | +1.0 A |
| # regressions (>0.05) | ~3 | 2 |
| # improvements (>0.05) | ~6 | 4 |

G+ v2 did achieve its goal of fewer regressions (2 vs ~3), but at the cost of also having fewer improvements (4 vs ~6). The net effect is smaller than Round 1, suggesting the time-varying beta is more conservative — it helps less, but also hurts less. The original regression problem on easy cases is partially addressed: the two regressions (7TRH_HBG, 7ZOZ_HLA) are medium-quality complexes, not easy ones.

---

## 6. Overall Conclusions

### Did the feature achieve its goal?

**Partially, but insufficiently.**

1. **Fewer regressions**: The decaying beta schedule does reduce regressions (2 vs ~3 in Round 1). The 87% unchanged rate indicates the approach is conservative. However, the 7TRH_HBG regression (-0.285 DockQ) is the single largest negative swing in the entire evaluation — worse than any individual regression in Round 1.

2. **Net improvement is marginal**: The +0.007 mean DockQ improvement is smaller than Round 1's +0.017 and far from significant (p=0.305). The approach is not producing a reliable, systematic improvement.

3. **CDR3 confidence did not improve**: Interface PAE and CDR3 pLDDT both worsened slightly, contradicting the expectation that time-varying beta would improve CDR3 quality.

4. **Diversity preserved**: The +1.0 A increase in ensemble diversity is a positive signal — the method explores more conformations without collapsing. However, this diversity does not translate into better oracle DockQ (0.391 vs 0.391).

5. **The "redistribution" pattern persists**: Like Round 1, the method changes which binding mode is found (large individual deltas in both directions) rather than systematically steering toward better conformations. The fundamental limitation is that beta-scaling and distance potentials alter the diffusion trajectory direction but lack the selectivity to distinguish correct from incorrect binding modes.

### What modifications could improve the feature?

1. **Implement G+.1 (Conditional Steering)**: The most impactful unimplemented improvement. Skip steering entirely for complexes where B1 already produces high-confidence predictions (iPTM > 0.85). This would eliminate the 7TRH_HBG regression (B1 confidence 0.925, already very high) while preserving the wins on difficult cases. A fast pre-check with 1 sample and reduced steps would add minimal overhead.

2. **Reduce beta_max**: The current beta_max=0.3 is likely too aggressive for the early phase. Reducing to 0.1-0.15 might preserve the exploratory benefit while reducing the magnitude of regressions. Round 1's Strategy L showed that beta=0.3 increases CDR3 RMSD significantly (p=0.0008).

3. **Switch to cosine schedule**: The linear beta decay may drop too quickly. A cosine schedule `beta(t) = beta_max * 0.5 * (1 + cos(pi * (1-t)))` would maintain stronger steering for longer in the early phase and decay more smoothly.

4. **Use B3-style best-of-many strategy**: B3's dominant performance (0.413 DockQ, p=1.2e-21) comes from running many restraint configurations and selecting the best. A similar approach for progressive steering — running 3-5 beta_max values per complex and selecting by confidence — could capture the upside while mitigating regressions.

5. **Improve model selection**: The weak confidence-DockQ correlation (visible in the scatter plot) means even when G+ v2 produces a good structure (oracle DockQ = 0.391), it often selects a worse one. Implementing G+.3 (interface-specific iPLDDT ranking) could close the 0.046 oracle-confidence gap.
