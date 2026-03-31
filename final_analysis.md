# Final Analysis — Strategy Y+ (Hierarchical Steering v2)

## Feature Goal

Strategy Y+ aims to improve antibody-antigen docking accuracy in Boltz-2 structure predictions **without using any oracle information** (no experimental contacts, no known epitopes). It uses a two-phase hierarchical approach:

- **Early phase**: Amplify attention between CDR residue pairs in embedding space (time-varying beta-scaling with beta_max=0.8), shaping loop conformations while coordinates are still noisy.
- **Late phase**: Apply coordinate-space potentials that pull the antibody's CDR loops toward the antigen surface and orient the antigen binding face toward the antibody.

The predecessor (Strategy Y, Round 1) was the only strategy out of 11 tested with a statistically significant structural improvement: CDR-H3 RMSD -0.30 A (p=0.025). Y+ strengthens both phases and adds two new capabilities: (Y+.3) epitope-focused guidance from predicted PDB structures, and (Y+.4) per-complex guidance weight scaling.

**Specific targets for Y+** (from STRATEGY_UPGRADES_GUIDE.md):
- CDR-H3 RMSD improvement beyond -0.30 A
- DockQ (Ab-Ag) improvement beyond +0.034
- CAPRI Medium+ rate increase beyond +6.4 percentage points
- Improved epitope prediction accuracy

## Evaluation Setup

- **Complexes**: 47 antibody-antigen complexes (ground truths from pdb_minimized/)
- **Baselines**: (1) Baseline — vanilla Boltz-2, (2) Contact Restraints — oracle contacts, (3) Pocket Guided — pocket-based restraints
- **Model selection**: For each complex, 5 models were predicted. The **top-1 model by confidence_score** was selected for all reported metrics (realistic scenario with no oracle knowledge). Oracle (best by DockQ) and mean-over-5 are also reported.
- **Tool**: DockQ v2 for interface quality; BioPython for CDR RMSD after antibody-framework alignment; epitope contacts at 5A/8A thresholds.
- **Statistical test**: Paired Wilcoxon signed-rank (non-parametric, appropriate for n=47).

---

## Plot-by-Plot Analysis

### Figure 1: DockQ Bar Plot (3 subplots)

**How produced**: Mean DockQ across all complexes, with 95% bootstrap confidence intervals. Top-1 model selected by confidence_score.

**Subplot 1 — Top-1 Total DockQ** (averaged across all chain-chain interfaces):

| Method | Mean | 95% CI |
|--------|------|--------|
| Baseline | 0.338 | [0.292, 0.391] |
| Contact Restraints | 0.499 | [0.441, 0.555] |
| Pocket Guided | 0.364 | [0.301, 0.426] |
| Hierarchical v2 | 0.361 | [0.304, 0.422] |

Y+ (0.361) shows a marginal improvement over baseline (0.338), a delta of +0.023. This is smaller than the +0.026 delta observed in Round 1 for Strategy Y. The confidence intervals overlap substantially. The improvement is **not statistically significant** (p=0.582).

**Subplot 2 — Top-1 Ab-Ag DockQ** (antibody-antigen interfaces only, the most relevant metric):

| Method | Mean | 95% CI |
|--------|------|--------|
| Baseline | 0.179 | [0.111, 0.255] |
| Contact Restraints | 0.408 | [0.323, 0.492] |
| Pocket Guided | 0.217 | [0.129, 0.312] |
| Hierarchical v2 | 0.213 | [0.137, 0.298] |

Y+ (0.213) shows +0.033 improvement over baseline (0.179), nearly identical to Round 1's +0.034 delta. However, this is again **not statistically significant** (p=0.589). Contact restraints remain far ahead at 0.408.

**Subplot 3 — Oracle Total DockQ** (best model of 5 selected by actual DockQ, upper bound):

| Method | Mean |
|--------|------|
| Baseline | 0.391 |
| Contact Restraints | 0.515 |
| Pocket Guided | 0.399 |
| Hierarchical v2 | 0.406 |

Oracle DockQ for Y+ (0.406) is slightly higher than baseline (0.391), indicating that among the 5 samples, Y+ does produce better structures — but the confidence-based selector doesn't always find them.

**Conclusion**: Y+ produces a small positive trend in DockQ but fails to achieve statistical significance. The gap to contact restraints remains very large (~2x on Ab-Ag DockQ).

---

### Figure 2: DockQ Boxplot

**How produced**: Distribution of Ab-Ag DockQ (top-1 by confidence) across all complexes, shown as box-and-whisker plots.

All four methods share a bimodal distribution: a large cluster near DockQ~0 (failed docking) and a spread of successes above 0.2. Key observations:

- **Baseline**: Median at 0.038, with most complexes failing (Q3 ~ 0.21). A few outliers reach 0.6-0.85.
- **Contact Restraints**: Median at 0.414, clearly shifted upward. The interquartile range spans ~0.07-0.68, meaning it rescues many of the failure cases.
- **Pocket Guided**: Median at 0.041, similar to baseline, but with a higher Q3 (~0.50).
- **Hierarchical v2**: Median at 0.033, essentially identical to baseline. However, the upper quartile (~0.51) is slightly higher than baseline's (~0.21), and the whiskers extend to ~0.89.

**Conclusion**: Y+ does not shift the median — the majority of complexes still fail at docking. However, Y+ has a slightly higher upper tail, suggesting it can produce dramatic wins on specific complexes while leaving the hard-failure cluster untouched.

---

### Figure 3: Y+ vs Baseline Scatter Plot

**How produced**: Per-complex Ab-Ag DockQ (top-1 by confidence), Y+ on y-axis vs baseline on x-axis. Points above the diagonal are improvements.

Key observations:
- A dense cluster at the origin (~20 complexes with DockQ near 0 for both methods) — these are global orientation failures that neither method resolves.
- Most points along the diagonal — Y+ and baseline give similar results on most complexes.
- **Dramatic wins** (labeled, far above diagonal):
  - **8FXB_HLE**: 0.009 -> 0.890 (+0.881) — near-complete failure rescued to high quality
  - **8OL9_BAH**: 0.193 -> 0.738 (+0.545)
  - **8HGM_CDB**: 0.228 -> 0.517 (+0.288)
  - **8EZ8_HLA**: 0.009 -> 0.297 (+0.288)
- **Notable losses** (below diagonal):
  - **8EZ7_HLA**: 0.150 -> 0.018 (-0.132)
  - **7ZOZ_HLA**: 0.644 -> 0.533 (-0.111)
- Overall: 16 wins, 20 losses, 11 ties.

**Conclusion**: Y+ produces spectacular improvements on a small number of complexes (4-5 cases with >0.1 gain) but slightly more losses (20) than wins (16) overall. The wins are larger in magnitude than the losses: the mean delta is +0.033 despite more losses, because the wins are dramatic (up to +0.881) while losses are modest (max -0.132). This high-variance, high-upside pattern is characteristic of steering methods that sometimes find the correct binding mode but occasionally push the structure away from a borderline-correct pose.

---

### Figure 4: CDR Loop RMSD Bar Plot

**How produced**: Mean Ca RMSD per CDR loop after superimposing on antibody framework (non-CDR) Ca atoms (top-1 model by confidence). Error bars are standard deviations.

| CDR Loop | Baseline | Contact Restr. | Pocket Guided | Hierarchical v2 | Delta (v2 - base) |
|----------|----------|---------------|---------------|-----------------|-------------------|
| CDR-H1 | 1.36 | 1.27 | 1.34 | 1.33 | -0.03 |
| CDR-H2 | 1.13 | 1.10 | 1.12 | 1.16 | +0.03 |
| **CDR-H3** | **3.00** | **2.35** | **2.64** | **2.84** | **-0.16** |
| CDR-L1 | 1.25 | 1.15 | 1.21 | 1.29 | +0.03 |
| CDR-L2 | 0.95 | 0.96 | 0.90 | 0.95 | +0.00 |
| CDR-L3 | 1.34 | 1.27 | 1.35 | 1.45 | +0.11 |

CDR-H3 (the most critical loop) shows a -0.16 A improvement (3.00 -> 2.84), which is a positive direction but **smaller than Round 1's -0.30 A** and **not statistically significant** (p=0.846). CDR-L3 shows a slight degradation (+0.11 A). Other loops are essentially unchanged.

**Conclusion**: Y+ shows a modest CDR-H3 improvement trend but does not replicate the statistically significant -0.30 A improvement from Round 1. The CDR-L3 degradation is a concern that was not present in Round 1. Contact restraints achieve the best CDR-H3 RMSD (2.35 A), showing what is achievable with oracle information.

---

### Figure 5: CAPRI Quality Classification (Stacked Bar)

**How produced**: Each complex classified by Ab-Ag DockQ (top-1 by confidence) into CAPRI categories: Incorrect (<0.23), Acceptable (0.23-0.49), Medium (0.49-0.80), High (>0.80).

| Method | Incorrect | Acceptable | Medium | High | Medium+High |
|--------|-----------|-----------|--------|------|------------|
| Baseline | 76.6% | 4.3% | 17.0% | 2.1% | 19.1% |
| Contact Restraints | 37.0% | 19.6% | 37.0% | 6.5% | 43.5% |
| Pocket Guided | 69.4% | 2.8% | 25.0% | 2.8% | 27.8% |
| Hierarchical v2 | 68.1% | 4.3% | 23.4% | 4.3% | 27.7% |

Y+ reduces the Incorrect rate from 76.6% to 68.1% (-8.5 pp) and increases Medium+High from 19.1% to 27.7% (+8.6 pp). This is a meaningful shift — roughly 4 additional complexes moved from failure to medium/high quality. The High category doubles from 2.1% to 4.3% (1 -> 2 complexes).

Y+ matches pocket-guided performance (27.8% Medium+High) despite using no oracle information, and exceeds baseline by +8.6 pp — slightly better than Round 1's +6.4 pp.

**Conclusion**: CAPRI classification is the metric where Y+ shows its clearest benefit. However, ~68% of complexes remain in the Incorrect category, indicating that the fundamental problem of global orientation failure is not solved.

---

### Figure 6: Epitope Prediction Quality

**How produced**: Mean F1 score comparing predicted vs native epitope residues (antigen residues contacting the antibody) at 5A (heavy-atom) and 8A (Ca-Ca) thresholds. Top-1 model by confidence.

| Threshold | Baseline | Contact Restr. | Pocket Guided | Hierarchical v2 |
|-----------|----------|---------------|---------------|-----------------|
| 5A | 0.345 | 0.587 | 0.368 | 0.336 |
| 8A | 0.311 | 0.552 | 0.328 | 0.310 |

Y+ (0.336 at 5A, 0.310 at 8A) is slightly below baseline (0.345, 0.311). The differences are small but consistently negative. The statistical test vs baseline shows p=0.116 at 5A — borderline concerning but not significant.

**Conclusion**: Y+ does **not** improve epitope prediction. Despite Y+.3 being designed to incorporate predicted epitope information, the epitope F1 scores are marginally worse than baseline. This was a key sub-goal that was not achieved. Note: Y+.3 (epitope from PDB predictions) was only available for 2 of 47 complexes (7TRH_HBG and 7TRI_ZYB), so the feature was largely untested in this evaluation.

---

### Figure 7: Confidence Calibration

**How produced**: Scatter plot of iptm (predicted confidence) vs Ab-Ag DockQ (actual quality) for each complex, colored by method. Pearson r shown in legend.

| Method | iptm mean | r(iptm, DockQ) |
|--------|-----------|----------------|
| Baseline | 0.738 | 0.25 |
| Contact Restraints | 0.742 | 0.43 |
| Pocket Guided | 0.773 | 0.34 |
| Hierarchical v2 | 0.756 | 0.14 |

Y+ has a slightly higher mean iptm (0.756 vs 0.738, p=0.068) but **worse confidence calibration** (r=0.14 vs 0.25). This means Y+ inflates confidence scores while the actual DockQ does not improve proportionally. Many complexes have high iptm (0.7-0.9) but near-zero DockQ for all methods, but Y+ exacerbates this pattern.

**Conclusion**: Y+ degrades the confidence-DockQ correlation. This is a practical concern: if the model is more confident but not more accurate, confidence-based model selection (which is the realistic deployment scenario) becomes less reliable. This may partially explain why the oracle DockQ (0.406) is noticeably better than the top-1 DockQ (0.361) — the selector is picking worse models with higher confidence.

---

## Summary Table

| Metric | Baseline | Y+ | Delta | p-value | Significant? |
|--------|----------|-----|-------|---------|-------------|
| Top-1 Total DockQ | 0.338 | 0.361 | +0.023 | 0.582 | No |
| Top-1 Ab-Ag DockQ | 0.179 | 0.213 | +0.033 | 0.589 | No |
| CDR-H3 RMSD (A) | 3.00 | 2.84 | -0.16 | 0.846 | No |
| Epitope F1 (5A) | 0.345 | 0.336 | -0.009 | 0.116 | No |
| Epitope F1 (8A) | 0.311 | 0.310 | -0.001 | 0.610 | No |
| iptm | 0.738 | 0.756 | +0.019 | 0.068 | No |
| CAPRI Medium+High | 19.1% | 27.7% | +8.6 pp | — | — |

---

## Overall Conclusions

### Did Y+ achieve its goals?

**Partially, but below expectations.**

1. **CDR-H3 RMSD**: Y+ shows a -0.16 A mean improvement, which is in the right direction but weaker than Round 1's -0.30 A and not statistically significant (p=0.846). The target of exceeding -0.30 A was not met.

2. **DockQ improvement**: +0.033 on Ab-Ag DockQ, nearly identical to Round 1's +0.034. The upgrades in Y+ (stronger beta-scaling, added CDR proximity potential, epitope guidance) did not amplify the effect beyond Round 1 levels. The target of exceeding +0.034 was not met.

3. **CAPRI rates**: +8.6 pp in Medium+High rate, exceeding Round 1's +6.4 pp. This is the one metric where Y+ shows clear improvement over its predecessor.

4. **Epitope prediction**: Slightly worse than baseline. Y+.3 (epitope from predicted PDBs) was only testable on 2 of 47 complexes, making it essentially unevaluated.

5. **Confidence calibration**: Degraded (r=0.14 vs 0.25). Y+ makes the model more confident without proportionally improving accuracy, making model selection harder.

### Key pattern

Y+ produces **high-variance, high-upside** results: spectacular improvements on a few complexes (8FXB_HLE: +0.881, 8OL9_BAH: +0.545) but slightly more regressions than improvements overall (20 losses vs 16 wins). The wins are large enough to shift the mean upward, but the signal is drowned by the noise of the ~68% of complexes that fail at global orientation.

### What could be modified

1. **Address the global orientation failure**: ~68% of complexes remain in the Incorrect category, suggesting the fundamental problem is not loop conformation but global docking orientation. A coarse rigid-body pre-docking step or much stronger early-phase orientation forcing could help. The current beta-scaling operates on CDR-CDR pair attention, which may not be sufficient to reorient the entire antibody relative to the antigen.

2. **Run epitope-guided predictions on all 48 complexes**: Y+.3 was only available for 2 complexes. Running a first-pass prediction and feeding those PDBs back through the pipeline for all complexes would properly test the epitope-focusing mechanism.

3. **Fix the confidence calibration degradation**: The higher iptm without higher DockQ means model selection is impaired. Consider a custom re-ranking function that incorporates interface-specific features rather than relying on iptm alone. Alternatively, combine Y+ with the FK resampling from Strategy A+ to improve particle diversity.

4. **Increase late-phase potential strength selectively**: The CDR-antigen proximity potential (activated at steering_t <= 0.4) may be too weak. The current weight of 0.8 could be increased for complexes where the CDR-antigen distance is large (Y+.4's guidance_weight_scale), but this was set to 1.0 for all complexes in this evaluation. A diagnostic pass to identify hard cases and scale up guidance for them could help.

5. **Combine with other strategies**: Round 2 proposals include N1 (A+Y Hybrid) combining Y's hierarchical timing with A's FK resampling. Given that Y+ alone fails to achieve significance, combining it with complementary mechanisms may be more productive than further tuning Y+ in isolation.

## Summary

Strategy Y+ aims to improve antibody-antigen docking accuracy in Boltz-2 structure predictions **without
using any oracle information** (no experimental contacts, no known epitopes). It uses a two-phase hiera                rchical approach:   
**Early phase**: Amplify attention between CDR residue pairs in embedding space (time-varying beta-sca
ling with beta_max=0.8), shaping loop conformations while coordinates are still noisy.   
 **Late phase**: Apply coordinate-space potentials that pull the antibody's CDR loops toward the antige
n surface and orient the antigen binding face toward the antibody.

The predecessor (Strategy Y, Round 1) was the only strategy out of 11 tested with a statistically signif
icant structural improvement: CDR-H3 RMSD -0.30 A (p=0.025). Y+ strengthens both phases and adds two new capabilities: (Y+.3) epitope-focused guidance from predicted PDB structures, and (Y+.4) per-complex guidance weight scaling.

- DockQ: +0.033 Ab-Ag DockQ improvement over baseline, but not statistically significant (p=0.589). Nearly identical    to Round 1's +0.034 — the Y+ upgrades did not amplify the effect.
- CDR-H3 RMSD: -0.16 A improvement, weaker than Round 1's -0.30 A and not significant (p=0.846).
- CAPRI rates: The one clear win — +8.6 pp in Medium+High, exceeding Round 1's +6.4 pp. 2 additional complexes rescued   from Incorrect to Medium/High.                                                                                         - Dramatic per-complex wins: 8FXB_HLE (+0.881), 8OL9_BAH (+0.545) — but offset by 20 losses vs 16 wins overall.
- Epitope prediction: Slightly worse than baseline. Y+.3 was only testable on 2/47 complexes.
- Confidence calibration: Degraded (r=0.14 vs 0.25) — model becomes overconfident without matching accuracy gains.
The report concludes with 5 specific modification proposals, the most impactful being: addressing the ~68% global       orientation failures, properly testing epitope guidance on all complexes, and combining with FK resampling (Strategy
A+).        