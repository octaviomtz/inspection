# Final Analysis: Strategy A+ — FK Particles v2

**Date**: 2026-03-31
**Evaluation**: 4 methods across 48 antibody-antigen complexes using DockQ v2
**Script**: `scripts/evaluate.py` with `conda run -n dockq2`

---

## 1. Goal of Strategy A+ (FK Particles v2)

In Round 1, Strategy A (Feynman-Kac particle steering with 20 particles) was ranked **#1 among all 11 strategies**, achieving the best DockQ improvement over baselines (+0.046 vs B1). However, it had three critical limitations:

1. **25% of complexes were missing** due to GPU out-of-memory errors at 20 particles
2. **Poor confidence calibration** (Spearman rho=0.265) — the model's confidence score could not reliably select the best prediction
3. **No statistical significance** (p=0.610) — improvements were concentrated in ~5 complexes

Strategy A+ was designed with four sub-strategies to address these issues:

| Sub-strategy | Goal |
|---|---|
| **A+.1** Adaptive Particle Reduction | Auto-reduce particle count when GPU memory is insufficient, recovering the missing 25% of complexes |
| **A+.2** Custom Re-Ranking by Steering Energy | Return FK energies from the diffusion loop and use a composite score (confidence + energy + interface pLDDT) to select the best model |
| **A+.3** Separate FK Seeds from Particles | Independent restart groups that resample internally for better diversity |
| **A+.4** Embedding Steering Fallback | Fall back to beta-scaling when particles are too expensive even at minimum count |

---

## 2. Methods Evaluated

| Method | Description | N complexes | N models |
|---|---|---|---|
| **B0** | Unconstrained Boltz-2 (no steering) | 47 | 235 (5 per complex) |
| **B1** | Contact restraints (oracle residue-pair restraints) | 46 | 2,610 (multiple settings x 5 models) |
| **B2** | Pocket restraints (oracle pocket residue restraints) | 36 | 435 (multiple settings x 5 models) |
| **NF_V2** | FK Particles v2 (this feature) | 33 | 165 (5 per complex) |

B1 uses ground-truth contact information and thus represents an **oracle upper bound** rather than a fair baseline. B2 also uses oracle pocket residues. NF_V2 uses only CDR region indices (derivable from sequence) and antigen MSA — no oracle structural information.

**Coverage**: NF_V2 covers 33/48 complexes (69%). This is an improvement over Round 1's 35/47 (74% for paired comparisons), but 15 complexes remain missing despite A+.1's adaptive particle reduction and A+.4's beta-scaling fallback. This suggests the OOM problem was only partially resolved.

---

## 3. Plot-by-Plot Analysis

### 3.1 Aggregate DockQ Bar Chart (`plot_aggregate_dockq.png`)

**How produced**: Mean DockQ Ab-Ag across all complexes available for each method. Two bars per method: solid = confidence selection (model with highest confidence score), hatched = oracle selection (model with highest actual DockQ). Each complex contributes one value — the DockQ of the selected model.

**Results**:

| Method | Confidence selection | Oracle selection |
|---|---|---|
| B0 | 0.179 | 0.248 |
| B1 | 0.291 | 0.479 |
| B2 | 0.190 | 0.269 |
| NF_V2 | 0.172 | 0.289 |

**Analysis**: Under confidence selection (the practical scenario), NF_V2 (0.172) performs **slightly below** B0 (0.179) and well below B1 (0.291). However, under oracle selection, NF_V2 (0.289) **outperforms** both B0 (0.248) and B2 (0.269), and approaches B1's confidence-selected performance. This reveals a critical gap: NF_V2's sample pool contains good predictions, but the confidence score fails to identify them. The oracle-vs-confidence gap for NF_V2 (0.289 - 0.172 = 0.117) is the **largest** among all methods, indicating severe model selection failure.

**Conclusion**: The FK particle mechanism generates structurally improved candidates, but the confidence-based selection pipeline destroys this advantage.

---

### 3.2 CAPRI Quality Distribution — Confidence Selection (`plot_capri_confidence.png`)

**How produced**: For each method, the confidence-selected model per complex is classified into CAPRI quality categories based on DockQ Ab-Ag thresholds: Incorrect (<0.23), Acceptable (0.23-0.49), Medium (0.49-0.80), High (>=0.80).

**Results**:

| Method | Incorrect | Acceptable | Medium | High |
|---|---|---|---|---|
| B0 | 77% | 4% | 17% | 2% |
| B1 | 54% | 13% | 30% | 2% |
| B2 | 72% | 3% | 22% | 3% |
| NF_V2 | 73% | 9% | 15% | 3% |

**Analysis**: NF_V2 has 73% Incorrect, nearly identical to B0 (77%) and B2 (72%). Its Medium-quality rate (15%) is actually the lowest of all methods. B1 is the clear leader with only 54% Incorrect and 30% Medium. NF_V2's slight edge is in Acceptable (9% vs 4% for B0), suggesting it nudges a few complexes past the Incorrect threshold without achieving Medium quality.

**Conclusion**: Under confidence selection, NF_V2 does not meaningfully improve the CAPRI quality distribution over the unconstrained baseline.

---

### 3.3 CAPRI Quality Distribution — Oracle Selection (`plot_capri_oracle.png`)

**How produced**: Same as 3.2 but selecting the model with the best actual DockQ per complex.

**Results**:

| Method | Incorrect | Acceptable | Medium | High |
|---|---|---|---|---|
| B0 | 68% | 2% | 28% | 2% |
| B1 | 22% | 22% | 50% | 7% |
| B2 | 64% | 3% | 31% | 3% |
| NF_V2 | 55% | 18% | 24% | 3% |

**Analysis**: Under oracle selection, NF_V2 reduces the Incorrect rate to **55%**, a substantial improvement over B0 (68%) and B2 (64%). NF_V2 achieves 18% Acceptable — the second-highest after B1 (22%). This confirms that the FK particle ensemble produces better candidate structures (more complexes with at least one acceptable prediction), but the practical pipeline cannot identify them.

**Conclusion**: The FK mechanism succeeds at generating better samples. The failure is in selection, not generation.

---

### 3.4 DockQ Box Plots — Confidence Selection (`plot_boxplot_confidence.png`)

**How produced**: Box plot of DockQ Ab-Ag values across complexes, one value per complex (the confidence-selected model). Individual data points overlaid with jitter.

**Analysis**: B0, B2, and NF_V2 show nearly identical distributions: median near 0.04, IQR from ~0.01 to ~0.25, similar upper whiskers. B1 has a visibly higher median (~0.10) and wider spread into higher DockQ values. NF_V2's distribution is indistinguishable from B0.

**Conclusion**: Under confidence selection, NF_V2 provides no distributional advantage over the unconstrained baseline.

---

### 3.5 DockQ Box Plots — Oracle Selection (`plot_boxplot_oracle.png`)

**How produced**: Same as 3.4 but using oracle-selected models.

**Analysis**: Under oracle selection, NF_V2 shows a noticeably higher median (~0.22) compared to B0 (~0.08) and B2 (~0.08). The upper quartile extends further. B1 remains dominant with the highest median (~0.54). The NF_V2 distribution shifts upward substantially compared to its confidence-selection counterpart, confirming the selection gap.

**Conclusion**: The particle ensemble genuinely improves the best-available predictions. Oracle NF_V2 outperforms oracle B0 and oracle B2.

---

### 3.6 Paired Scatter: NF_V2 vs B0 (`plot_scatter_nf_v2_vs_b0.png`)

**How produced**: Each point is one complex, plotting the confidence-selected DockQ Ab-Ag for B0 (x-axis) vs NF_V2 (y-axis). Points above the diagonal indicate NF_V2 wins. Only the 33 complexes where both methods have predictions are shown.

**Results**: NF_V2 wins on 14 complexes, B0 wins on 19, 0 tied.

**Analysis**: The majority of points cluster near the origin (both methods fail), with B0 winning slightly more often. A few notable outliers: NF_V2 achieves a large improvement on 7TRH_HBG (~0.49 → 0.69), and on several complexes around 0.6 DockQ, NF_V2 slightly underperforms B0. The cluster pattern shows most complexes are too hard for either method under confidence selection.

**Conclusion**: NF_V2 does not consistently beat B0 under confidence selection. The 14-19 win/loss ratio is slightly negative.

---

### 3.7 Delta Waterfall: NF_V2 vs B0 (`plot_delta_nf_v2_vs_b0.png`)

**How produced**: Per-complex DockQ difference (NF_V2 - B0) under confidence selection, sorted from most negative to most positive. Green bars = NF_V2 improvement, red bars = NF_V2 regression.

**Results**: Mean delta = +0.003. Roughly balanced between improvements and regressions, with the largest improvement (~+0.20 on 7TRH_HBG) and the largest regression (~-0.07).

**Analysis**: The mean improvement is negligible (+0.003). The waterfall is nearly symmetric, with no systematic shift toward improvement. A few large positive deltas (right side) are offset by many small negative deltas across the majority of complexes. The pattern suggests NF_V2 helps a few specific cases substantially but introduces small regressions across many others.

**Conclusion**: NF_V2 provides no systematic improvement over B0 under confidence selection. The mean delta of +0.003 is effectively zero.

---

### 3.8 Confidence Calibration (`plot_confidence_calibration.png`)

**How produced**: Four subplots (one per method). Each point is a single model (not per-complex — all 5 models per complex are shown). X-axis = model confidence score, y-axis = actual DockQ Ab-Ag. Spearman correlation (rho) and p-value shown.

**Results**:

| Method | rho | p-value | N models |
|---|---|---|---|
| B0 | 0.380 | 1.7e-09 | 235 |
| B1 | 0.443 | 4.2e-126 | 2,610 |
| B2 | 0.359 | 1.1e-14 | 435 |
| NF_V2 | 0.216 | 5.4e-03 | 165 |

**Analysis**: NF_V2 has the **weakest confidence calibration** (rho=0.216) among all methods — substantially lower than B0 (0.380), B1 (0.443), and B2 (0.359). While still statistically significant (p=0.005), the correlation is too weak for reliable model selection. The scatter plot for NF_V2 shows many high-confidence models with low DockQ, indicating that the FK steering process disrupts the confidence score's ability to discriminate good from bad predictions.

**Conclusion**: FK particle steering **degrades confidence calibration**. This is the primary reason the confidence-selected results are poor despite oracle-selected results being strong. The steering process likely inflates confidence scores for steered trajectories regardless of their actual structural quality.

---

### 3.9 FK Energy vs DockQ (`plot_fk_energy_vs_dockq.png`)

**How produced**: Scatter plot of FK steering energy (lower = better docking potential) vs actual DockQ Ab-Ag for all 165 NF_V2 models. Spearman correlation computed.

**Results**: rho = 0.022, p = 7.8e-01 (not significant).

**Analysis**: There is **zero correlation** between FK energy and DockQ quality. The scatter shows no trend whatsoever — models with very low FK energy (supposedly good steering) have DockQ values spanning the entire range, and high-DockQ models appear at all energy levels. This is a critical finding: the FK potential energies (ContactPotential + SymmetricChainCOMPotential) do not predict actual docking quality.

**Conclusion**: The FK energy is **not a useful signal for re-ranking** (A+.2's core hypothesis). The composite re-ranking score that incorporates FK energy cannot improve model selection because the energy term is uninformative. This directly undermines sub-strategy A+.2.

---

### 3.10 Composite vs Confidence Selection (`plot_composite_vs_confidence_selection.png`)

**How produced**: Each point is one complex. X-axis = DockQ of the confidence-selected model, y-axis = DockQ of the composite-score-selected model (alpha * confidence - beta * normalized_energy + gamma * interface_pLDDT). Only NF_V2 complexes shown.

**Results**: Composite wins 13, Confidence wins 14, Tied 6. Confidence mean: 0.172, Composite mean: 0.169.

**Analysis**: The composite re-ranking performs **identically** to pure confidence selection. Points scatter evenly around the diagonal. The composite score is marginally worse (0.169 vs 0.172 mean), though the difference is negligible. This is expected given the FK energy's zero correlation with DockQ (Section 3.9) — adding a noise term to the selection criterion cannot improve it.

**Conclusion**: Sub-strategy A+.2 (composite re-ranking by steering energy) **does not work**. The FK energy component adds noise rather than signal to the selection process.

---

### 3.11 Interface RMSD (`plot_irmsd.png`)

**How produced**: Mean iRMSD (interface RMSD for the antibody heavy chain interface, AB) across complexes, using confidence-selected models. Lower is better.

**Results**: B0 = 11.5, B1 = 7.6, B2 = 11.3, NF_V2 = 11.1.

**Analysis**: NF_V2 (11.1) is marginally better than B0 (11.5) and B2 (11.3) but far behind B1 (7.6). The improvement of 0.4 angstroms over B0 is negligible in the context of values >11 angstroms (which indicate incorrect docking poses). B1's advantage (7.6) reflects the benefit of oracle contact information.

**Conclusion**: NF_V2 does not meaningfully reduce interface RMSD under confidence selection.

---

### 3.12 Native Contact Recovery — fnat (`plot_fnat.png`)

**How produced**: Mean fraction of native contacts (fnat) for the AB interface across complexes, using confidence-selected models. Higher is better.

**Results**: B0 = 0.181, B1 = 0.289, B2 = 0.184, NF_V2 = 0.162.

**Analysis**: NF_V2 (0.162) is the **worst** of all methods on fnat, below even the unconstrained baseline B0 (0.181). This is surprising — in Round 1, Strategy A achieved 87% of B2's fnat. The regression suggests that the confidence selector is actively choosing models with poor native contacts, or that the steering modifications introduced in A+ (seeds, adaptive reduction, fallback) are less effective than the original Round 1 configuration.

**Conclusion**: NF_V2 slightly worsens native contact recovery under confidence selection, reversing Round 1's advantage.

---

### 3.13 Per-Complex Heatmap (`plot_heatmap_dockq.png`)

**How produced**: Heatmap of DockQ Ab-Ag for every complex (rows) and method (columns), using confidence-selected models. Color scale: red = 0 (incorrect), green = 1 (perfect). White cells = missing data (no prediction for that complex/method).

**Analysis**: The heatmap reveals several patterns:
1. **NF_V2 has the most white cells** (15 missing complexes), confirming coverage limitations
2. **Top performers** (8GQ1, 8CDE, 8CDD, 8J1T, 7WT9, 8BLQ) show green across all methods — these are easy complexes where all methods succeed, and NF_V2 provides no additional benefit
3. **Hard complexes** (bottom half, uniformly red) remain unsolved by all methods including NF_V2
4. **The few complexes where NF_V2 shows a green shift** (e.g., 7TRH_HBG) are the same ones identified in the waterfall plot
5. B1 shows the most green/yellow cells, especially for mid-difficulty complexes

**Conclusion**: NF_V2's improvements are limited to a small subset of complexes. It does not solve any complex that all other methods fail on.

---

## 4. Statistical Tests

### 4.1 Confidence Selection

| Comparison | N paired | NF_V2 mean | Baseline mean | Delta | p-value | Significant? |
|---|---|---|---|---|---|---|
| NF_V2 vs B0 | 33 | 0.172 | 0.170 | +0.003 | 0.447 | No |
| NF_V2 vs B1 | 32 | 0.178 | 0.281 | -0.103 | 0.014 | **Yes (worse)** |
| NF_V2 vs B2 | 24 | 0.206 | 0.178 | +0.028 | 0.768 | No |

### 4.2 Oracle Selection

| Comparison | N paired | NF_V2 mean | Baseline mean | Delta | p-value | Significant? |
|---|---|---|---|---|---|---|
| NF_V2 vs B0 | 33 | 0.289 | 0.242 | +0.047 | 0.126 | No |
| NF_V2 vs B1 | 32 | 0.298 | 0.501 | -0.203 | 0.000 | **Yes (worse)** |
| NF_V2 vs B2 | 24 | 0.299 | 0.280 | +0.019 | 0.768 | No |

**Analysis**: Under confidence selection, NF_V2 is statistically indistinguishable from B0 (p=0.447) and B2 (p=0.768). It is significantly **worse** than B1 (p=0.014), which is expected since B1 uses oracle contact information.

Under oracle selection, NF_V2 shows a promising +0.047 improvement over B0, but this is not statistically significant (p=0.126). Notably, this +0.047 delta matches the Round 1 result (+0.046), suggesting the FK mechanism's underlying effectiveness is unchanged — but the improvements remain concentrated in too few complexes to achieve significance.

---

## 5. Comparison with Round 1

| Metric | Round 1 (Strategy A) | Round 2 (NF_V2) | Change |
|---|---|---|---|
| DockQ improvement vs B0 (confidence) | +0.046* | +0.003 | Worsened |
| DockQ improvement vs B0 (oracle) | N/A | +0.047 | Similar |
| Coverage | 35/47 (74%) | 33/48 (69%) | Slightly worse |
| Confidence calibration (rho) | 0.265 | 0.216 | Worsened |
| CAPRI downgrades | 0 | ~19 complexes where B0 wins | Worsened (no longer a "safety net") |

*Round 1 compared against B1, not B0; and used 20 particles.

The Round 1 evaluation noted that Strategy A had "zero CAPRI downgrades" (it never made things worse), which was its key safety property. NF_V2 loses this property — B0 wins on 19/33 paired complexes under confidence selection.

---

## 6. Sub-Strategy Assessment

| Sub-strategy | Intended goal | Outcome |
|---|---|---|
| **A+.1** Adaptive Particle Reduction | Recover missing 25% of complexes | **Partially failed**: Still missing 31% of complexes (15/48). The auto-reduce + fallback mechanism did not fully solve OOM. |
| **A+.2** Composite Re-Ranking | Better model selection using FK energy | **Failed**: FK energy has zero correlation with DockQ (rho=0.022, p=0.78). Composite score performs identically to confidence (13 vs 14 wins). |
| **A+.3** Separate FK Seeds | Better diversity through independent restarts | **Inconclusive**: Cannot isolate the effect of seeds from the other changes. Overall oracle improvement (+0.047 vs B0) suggests some diversity benefit but it's not statistically significant. |
| **A+.4** Embedding Steering Fallback | Beta-scaling when particles too expensive | **Inconclusive**: The fallback was triggered for some complexes, but we cannot determine whether those complexes benefited or were harmed. |

---

## 7. Overall Conclusions

### What worked
1. **The FK particle mechanism still generates better candidate structures** — the oracle DockQ of NF_V2 (0.289) exceeds B0 (0.248) and B2 (0.269), confirming Round 1's finding that particle-based importance sampling improves structural sampling.
2. **CAPRI oracle analysis** shows NF_V2 reduces the Incorrect rate from 68% to 55%, meaning more complexes have at least one good prediction in their ensemble.

### What failed
1. **Confidence calibration is critically broken** — NF_V2 has the worst calibration of all methods (rho=0.216), which means the practical (confidence-selected) results are far worse than the oracle results. This is the single biggest failure mode.
2. **FK energy is uninformative** — The core hypothesis of A+.2 (that steering energy predicts docking quality) is definitively falsified (rho=0.022). The potentials (ContactPotential + SymmetricChainCOMPotential) steer trajectories but their energy does not correlate with the ground-truth metric.
3. **Coverage was not recovered** — Despite adaptive reduction and fallback, 31% of complexes remain missing, worse than Round 1's 25%.
4. **No statistical significance** — Neither confidence-selected (p=0.447) nor oracle-selected (p=0.126) improvements over B0 reach significance.

### Why NF_V2 underperforms Round 1's Strategy A
The most likely explanation is that the combination of reduced particle counts (from 20 to auto-reduced lower values), the seed mechanism, and the beta-scaling fallback diluted the effectiveness of the original FK approach. Round 1 used a fixed 20 particles; NF_V2's adaptive reduction likely lowered many complexes to 5-10 particles, reducing the ensemble's ability to find good orientations.

---

## 8. Recommended Modifications

Based on this analysis, the following modifications would address the identified failure modes:

1. **Fix confidence calibration (highest priority)**: The FK steering process must not corrupt the confidence score. Options:
   - Run confidence estimation on a clean (unsteered) forward pass after selecting the steered coordinates
   - Train a separate confidence head that accounts for steered trajectories
   - Use a completely different selection criterion (e.g., cluster analysis of the ensemble, consensus scoring)

2. **Replace FK energy with a geometry-based re-ranking signal**: Since FK energy is uninformative, replace it with post-hoc structural metrics that can be computed without ground truth:
   - Interface buried surface area
   - CDR-antigen contact count
   - Rosetta interface energy (computed after prediction)
   - Clash score at the interface

3. **Increase particle count instead of reducing it**: The adaptive reduction likely hurt more than it helped. Instead:
   - Use gradient checkpointing or mixed precision to fit more particles in memory
   - Run complexes in batches with different particle counts based on their size
   - Accept that some very large complexes cannot use FK particles and exclude them rather than degrading them

4. **Investigate the Round 1 configuration**: Run the exact Round 1 A configuration (20 particles, no seeds, no adaptive reduction, no fallback) on the full 48-complex test set to establish whether the regression is due to A+ modifications or simply different test set composition.

5. **Combine with B1-style restraints**: Since B1 dramatically outperforms all methods, a hybrid approach using both FK particles and contact restraints could be promising — using the restraints for global orientation and particles for refinement.

---

**Document Status**: Final evaluation of Strategy A+ (FK Particles v2)
**Last Updated**: 2026-03-31
