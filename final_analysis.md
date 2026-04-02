# N1: A+Y Hybrid (FK Particles + Hierarchical Timing) — Final Analysis

**Date**: 2026-04-02
**Branch**: `a_v2_fk_particles`
**Evaluation script**: `scripts/evaluate_n1.py`

---

## 1. Original Goal

N1 was designed to combine the two most promising Round 1 strategies:

- **Strategy A** (FK particle resampling): Achieved the best DockQ improvement in Round 1 (+0.046 DockQ Ab-Ag, confidence selection), though not statistically significant (p=0.610).
- **Strategy Y** (Hierarchical timing): Achieved the only statistically significant structural improvement in Round 1 (CDR-H3 RMSD -0.30A, p=0.025).

The hypothesis was that combining both would yield additive or synergistic improvements: the early phase (CDR3 beta-scaling on pair representations) would guide global antibody orientation, while the late phase (FK particles with AntigenOrientationPotential) would refine local docking geometry.

**Quantitative targets** (from N1_EVALUATION_STRATEGY.md):
- DockQ Ab-Ag: +0.05 to +0.08 vs B0 (unconstrained baseline)
- CDR-H3 RMSD: -0.3 to -0.5A vs B0
- Statistical significance: p < 0.05 (Wilcoxon signed-rank)
- Regression control: fewer than 5 easy-case regressions

---

## 2. Baselines

| ID | Method | Models/complex | Description |
|----|--------|---------------|-------------|
| B0 | Unconstrained | 5 | Standard Boltz-2, no steering |
| B1 | Contact restraints | ~57 (across multiple settings) | Known contact pairs as restraints |
| B2 | Pocket restraints | ~12 (across multiple settings) | Pocket-based restraints |
| N1 | Hybrid FK+Hier. | 5 | This feature (FK particles + hierarchical timing) |

**Important**: B1 and B2 have many more models per complex than N1 and B0 (B1 averages ~57 models across ~16 settings, B2 averages ~12 across ~3 settings). This gives B1 and B2 a significant advantage under oracle (best-of-N) selection, since more samples increase the chance of a lucky hit. Confidence-based selection (picking the model with highest confidence score) is the fairer comparison across methods.

---

## 3. Plot-by-Plot Analysis

### 3.1 Aggregate DockQ by Method (`plot_aggregate_dockq.png`)

**How produced**: For each complex, either the confidence-selected best model (highest confidence score) or the oracle-selected best model (highest DockQ Ab-Ag) was chosen. The DockQ Ab-Ag was then averaged (mean) across all complexes per method. Two bars per method show confidence vs oracle selection.

**Results**:

| Method | Confidence | Oracle |
|--------|-----------|--------|
| B0 | 0.179 | 0.248 |
| B1 | 0.291 | 0.479 |
| B2 | 0.190 | 0.269 |
| N1 | 0.189 | 0.253 |

**Conclusion**: N1 achieves only +0.0094 DockQ over B0 with confidence selection (target was +0.05 to +0.08). With oracle selection the gap narrows further to +0.0048. B1 dominates both selection methods, largely benefiting from its ~10x higher model count under oracle selection. N1 shows no meaningful improvement over B0 or B2. The combination of strategies A and Y did not produce the expected additive gains.

### 3.2 CAPRI Quality Distribution — Confidence Selection (`plot_capri_confidence.png`)

**How produced**: Each complex's confidence-selected model was classified by CAPRI quality thresholds (Incorrect: DockQ < 0.23, Acceptable: 0.23-0.49, Medium: 0.49-0.80, High: > 0.80). Percentages shown as stacked bars.

**Results**:

| Method | Incorrect | Acceptable | Medium | High |
|--------|-----------|------------|--------|------|
| B0 | 77% | 4% | 17% | 2% |
| B1 | 54% | 13% | 30% | 2% |
| B2 | 72% | 3% | 22% | 3% |
| N1 | 72% | 6% | 19% | 2% |

**Conclusion**: N1 reduces the Incorrect rate by 5 percentage points (77% to 72%) vs B0, with small gains distributed across Acceptable (+2%) and Medium (+2%). This is a marginal improvement. B1 achieves a much larger shift (23 percentage points fewer Incorrect). N1 is essentially on par with B2.

### 3.3 CAPRI Quality Distribution — Oracle Selection (`plot_capri_oracle.png`)

**How produced**: Same as above but using the oracle-selected (best DockQ) model per complex.

**Results**:

| Method | Incorrect | Acceptable | Medium | High |
|--------|-----------|------------|--------|------|
| B0 | 68% | 2% | 28% | 2% |
| B1 | 22% | 22% | 50% | 7% |
| B2 | 64% | 3% | 31% | 3% |
| N1 | 62% | 11% | 23% | 4% |

**Conclusion**: Under oracle selection, N1 drops to 62% Incorrect (vs B0's 68%), a 6-point improvement. The Acceptable category gains the most (+9%), though Medium actually drops slightly compared to B0 (23% vs 28%), suggesting N1's oracle improvements come from pushing a few complexes from Incorrect into Acceptable rather than broadly improving quality. B1's dramatic 22% Incorrect rate (benefiting from ~57 models per complex) confirms that model count is the strongest lever for oracle performance.

### 3.4 DockQ Distribution — Confidence Selection (`plot_boxplot_confidence.png`)

**How produced**: Box plot of DockQ Ab-Ag values across all complexes, one box per method, using confidence-selected models. Overlaid scatter points show individual complexes.

**Conclusion**: N1's distribution is virtually identical to B0 — same median near 0.04, same interquartile range, same whisker extent. B1 shows a noticeably higher median and wider upper quartile. The overlaid scatter reveals that most complexes cluster near zero for all methods, with a few high-DockQ outliers. N1 does not shift the distribution meaningfully.

### 3.5 DockQ Distribution — Oracle Selection (`plot_boxplot_oracle.png`)

**How produced**: Same as above but with oracle-selected models.

**Conclusion**: Under oracle selection, all methods shift upward, but the relative ordering is preserved. B1 clearly separates from the pack with a median around 0.54. N1, B0, and B2 remain overlapping, with N1 showing a marginally higher upper whisker but the same median as B0. The oracle gap between N1 and B0 is not visually distinguishable.

### 3.6 Per-Complex Scatter: N1 vs B0 (`plot_scatter_n1_vs_b0.png`)

**How produced**: Each point is one complex, with B0's confidence-selected DockQ Ab-Ag on the x-axis and N1's on the y-axis. Points above the diagonal indicate N1 improvement; below indicate regression. Uses confidence selection.

**Results**: N1 wins on 23 complexes, B0 wins on 24, with 0 ties (using the plot's threshold).

**Conclusion**: The win/loss ratio is essentially 50/50. Points cluster along the diagonal, indicating N1 and B0 produce very similar results complex-by-complex. A few dramatic outliers exist in both directions — some complexes see large N1 improvement (points far above diagonal) while others see large regressions (points far below). The lack of systematic above-diagonal displacement confirms the non-significant p-value (0.8054).

### 3.7 Per-Complex Improvement Waterfall: N1 vs B0 (`plot_delta_n1_vs_b0.png`)

**How produced**: For each complex, delta = N1_DockQ - B0_DockQ (confidence selection). Complexes sorted by delta from most negative (left, red) to most positive (right, green). Horizontal dashed line at the mean delta (+0.009).

**Conclusion**: The waterfall reveals the core problem: improvements and regressions roughly cancel out. The mean delta of +0.009 is barely positive. The most severe regression is 7ZOZ_HLA (-0.591, from 0.644 to 0.054 — a catastrophic loss of an already well-docked complex). The largest improvement is approximately +0.4 on a hard case. The asymmetry — large regressions on easy cases, moderate improvements on hard cases — suggests the steering mechanism disrupts already-good predictions while occasionally rescuing bad ones.

**Easy-case regression analysis**: Of 11 complexes where B0 achieves DockQ > 0.3, N1 regresses on 8:
- 7ZOZ_HLA: 0.644 to 0.054 (delta -0.591, catastrophic)
- 8CDD_EDB: 0.728 to 0.550 (delta -0.178)
- 8CDE_DCB: 0.744 to 0.680 (delta -0.065)
- 8IV4_ABG: 0.364 to 0.319 (delta -0.045)
- 4 others with deltas -0.017 to -0.021

This exceeds the target of fewer than 5 easy-case regressions.

**Hard-case improvement analysis**: Of 33 complexes where B0 achieves DockQ < 0.1, N1 improves only 4:
- 8EZ8_HLA: 0.009 to 0.416 (delta +0.407)
- 8TRS_AGD: 0.037 to 0.293 (delta +0.257)
- 8EZ3_HLA: 0.020 to 0.083 (delta +0.063)
- 8BYU_HLA: 0.046 to 0.078 (delta +0.031)

The remaining 29 hard cases show negligible change (ties).

### 3.8 Confidence Calibration (`plot_confidence_calibration.png`)

**How produced**: Four subplots, one per method. Each point is a single model (not per-complex — all 5 models plotted). X-axis: confidence score, Y-axis: DockQ Ab-Ag. Spearman rho and p-value annotated. This uses all models, not just confidence-selected ones.

**Results**:

| Method | N models | Spearman rho | p-value |
|--------|----------|-------------|---------|
| B0 | 235 | 0.380 | 1.7e-9 |
| B1 | 2610 | 0.443 | 4.2e-126 |
| B2 | 435 | 0.359 | 1.1e-14 |
| N1 | 235 | 0.416 | 3.1e-11 |

**Subplot analysis**:
- **B0 subplot**: Moderate positive correlation (rho=0.380). Points spread across confidence 0.70-0.95, with DockQ mostly below 0.4. A few high-confidence, high-DockQ outliers.
- **B1 subplot**: Highest correlation (rho=0.443) with 2610 points forming a dense cloud. The high model count produces a smooth correlation curve. Higher-confidence models tend to be better.
- **B2 subplot**: Weakest correlation (rho=0.359). Points are more scattered, suggesting pocket restraints introduce noise in the confidence-DockQ relationship.
- **N1 subplot**: Second-best correlation (rho=0.416), improved over B0 (+0.036). This is a positive finding — the hybrid steering makes confidence a slightly better predictor of quality, even though overall quality didn't improve. The scatter pattern is similar to B0 but with a marginally tighter trend.

**Conclusion**: N1 improves confidence calibration over B0 (rho 0.416 vs 0.380). This is a modest positive result, suggesting the hybrid mechanism makes the model's internal confidence more reflective of actual quality. However, since the overall DockQ didn't improve, better calibration does not translate to better model selection outcomes.

### 3.9 Interface RMSD (`plot_irmsd.png`)

**How produced**: Mean iRMSD of the antibody chain (AB interface) across all complexes using confidence-selected models. Lower is better. Single bar per method.

**Results**: B0=11.5A, B1=7.6A, B2=11.3A, N1=11.4A.

**Conclusion**: N1's interface RMSD (11.4A) is essentially identical to B0 (11.5A) and B2 (11.3A). Only B1 achieves a meaningful reduction (7.6A). The hybrid steering produces no improvement in interface structural accuracy.

### 3.10 Native Contact Recovery (`plot_fnat.png`)

**How produced**: Mean fraction of native contacts (fnat) at the AB interface across all complexes using confidence-selected models. Higher is better. Single bar per method.

**Results**: B0=0.181, B1=0.289, B2=0.184, N1=0.194.

**Conclusion**: N1 shows a small improvement over B0 (+0.013 fnat), recovering about 1.3% more native contacts. B1 remains the best (0.289). The improvement is consistent with the marginal DockQ gains and suggests the FK particles provide a weak but non-zero signal for contact recovery.

### 3.11 Per-Complex Heatmap (`plot_heatmap_dockq.png`)

**How produced**: 47 rows (complexes) x 4 columns (methods). Color scale: dark red = 0.0 (poor), green = 1.0 (excellent). Rows sorted by B0 DockQ descending. Uses confidence-selected models. White cells indicate missing data (B2 was not run on all complexes).

**Conclusion**: The heatmap confirms the overall pattern:
- The top rows (easy complexes) show green/yellow across all methods — N1 generally preserves these but occasionally introduces darker cells (regressions).
- The bottom rows (hard complexes) are uniformly dark red across all methods — N1 fails to rescue most of them.
- B1 shows notably lighter colors across more rows, especially in the mid-range complexes.
- N1 and B0 columns are visually nearly indistinguishable, reinforcing that the hybrid approach produces negligible aggregate change.
- A few individual rows show dramatic color differences (e.g., the 7ZOZ_HLA regression visible as a green-to-red transition between B0 and N1).

---

## 4. Statistical Tests

### Confidence Selection

| Comparison | N | N1 mean | Baseline mean | Delta | p-value | Significant? |
|------------|---|---------|---------------|-------|---------|-------------|
| N1 vs B0 | 47 | 0.189 | 0.179 | +0.009 | 0.805 | No |
| N1 vs B1 | 46 | 0.193 | 0.291 | -0.099 | 0.002 | Yes (N1 loses) |
| N1 vs B2 | 36 | 0.190 | 0.190 | +0.001 | 0.871 | No |

### Oracle Selection

| Comparison | N | N1 mean | Baseline mean | Delta | p-value | Significant? |
|------------|---|---------|---------------|-------|---------|-------------|
| N1 vs B0 | 47 | 0.253 | 0.248 | +0.005 | 0.618 | No |
| N1 vs B1 | 46 | 0.258 | 0.479 | -0.221 | <0.001 | Yes (N1 loses) |
| N1 vs B2 | 36 | 0.253 | 0.269 | -0.016 | 0.085 | No |

**Conclusion**: N1 does not achieve statistically significant improvement over any baseline. The only significant result is that N1 is significantly *worse* than B1. The p-value of 0.805 against B0 indicates the improvements and regressions are effectively random noise.

---

## 5. CDR-H3 RMSD Analysis

CDR-H3 RMSD was not computed for any method (n=0 in all CDR_H3_RMSD columns). This was a secondary goal from Strategy Y (target: -0.3 to -0.5A improvement). Without this metric, we cannot assess whether the hierarchical timing component (early-phase beta-scaling) contributed any structural benefit to CDR-H3 loop accuracy.

This is a critical gap: Strategy Y's only significant result in Round 1 was on CDR-H3 RMSD, and validating whether the hybrid preserved this benefit was a key evaluation goal.

---

## 6. Overall Conclusions

### Goal Achievement Summary

| Target | Achieved | Result |
|--------|----------|--------|
| DockQ +0.05 to +0.08 vs B0 | No | +0.009 (18% of minimum target) |
| CDR-H3 RMSD -0.3 to -0.5A | Unknown | Not computed (n=0) |
| Statistical significance (p<0.05) | No | p=0.805 |
| <5 easy-case regressions | No | 8 regressions out of 11 easy cases |
| Confidence calibration preserved | Yes | rho improved 0.380 to 0.416 |

**The N1 hybrid strategy failed to achieve its primary goals.** The combination of Strategy A (FK particles) and Strategy Y (hierarchical timing) did not produce additive improvements. Instead, the effect is indistinguishable from the unconstrained baseline.

### Why Did the Hybrid Fail?

1. **Optimization conflict**: The early phase (beta-scaling for CDR3 orientation) and late phase (FK particles for docking refinement) address different objectives. Early-phase modifications to pair representations may push the trajectory in a direction that is suboptimal for the late-phase FK potential, effectively creating an optimization tug-of-war.

2. **Easy-case disruption**: The early-phase beta-scaling appears to destabilize complexes that Boltz-2 already handles well. Of 11 easy cases (B0 DockQ > 0.3), 8 regressed with N1, including one catastrophic failure (7ZOZ_HLA: 0.644 to 0.054). The steering mechanism overcorrects on cases that don't need correction.

3. **Insufficient late-phase strength**: With only 5 FK particles and a 50% transition point, the late phase may not have enough sampling capacity to recover from early-phase perturbations or to meaningfully explore the energy landscape.

4. **Diminished individual strategy effects**: In Round 1, Strategy A alone achieved +0.046 DockQ (though not significant), and Strategy Y alone achieved significant CDR-H3 improvement. Combining them appears to dilute both effects rather than amplify them — possibly because each strategy's hyperparameters were tuned independently and their interaction was not optimized.

---

## 7. Proposed Modifications

If the hybrid approach is to be revisited:

1. **Case-adaptive steering**: Only apply N1 steering on complexes where B0 performs poorly (DockQ < 0.1). Leave easy cases unconstrained. This could preserve the 4 hard-case improvements while avoiding the 8 easy-case regressions.

2. **Decouple and compare**: Run Strategy A and Strategy Y independently on the full test set (as in Round 1) to confirm their individual effects before combining. If individual effects are not reproducible, the hybrid has no foundation.

3. **Tune transition timing**: The current 50% transition fraction was set without optimization. Earlier transition (e.g., 30%) gives more steps to FK particles for refinement; later transition (e.g., 70%) gives more time for beta-scaling to establish orientation. A grid search over transition_fraction could reveal the optimal balance.

4. **Increase late-phase capacity**: Use more FK particles (10-20 instead of 5) to improve sampling in the refinement phase. Strategy A's Round 1 results used 20 particles.

5. **Abandon the hybrid**: Given the negligible improvement and the strong performance of B1 (contact restraints), resources may be better spent on strategies that incorporate known contact information rather than attempting blind steering through embedding modifications. B1's consistent superiority (0.291 vs 0.189 confidence DockQ) suggests that explicit structural knowledge outperforms implicit gradient guidance.

6. **Ensemble approach**: Rather than combining A and Y within a single run, generate predictions from each independently and select the best using a re-ranker (as proposed in the Round 1 N15 strategy). This avoids the optimization conflict while capturing complementary benefits.
