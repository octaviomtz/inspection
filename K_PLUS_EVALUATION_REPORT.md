# K+ Evaluation Report: Region-Specific Beta-Scaling for Antibody-Antigen Docking

## 1. Experimental Setup

**Goal**: K+ (region-specific beta-scaling) aims to discover epitope locations on unknown antigens by emphasizing different antigen surface regions during diffusion. The evaluation tests whether K+ improves (a) epitope prediction accuracy and (b) docking/interface quality compared to baselines.

**Methods evaluated**:
- **B1** (47 complexes): Vanilla Boltz2, no constraints
- **B2** (522 configs across 46 complexes): Contact restraints (hbond, hydrophobic, salt-bridge), ~11 configs/complex
- **B3** (87 configs across 36 complexes): Pocket restraints to specific residues, ~2.4 configs/complex
- **K+** (47 complexes): CDR3 beta-scaling (our feature)

**Model selection**: For each method/complex/config, the model with the highest `confidence_score` was selected from up to 5 candidates. All bar charts show the **mean across complexes** with standard deviation error bars. The scatter plot shows **per-complex best-by-confidence** values. DockQ quality distribution uses **raw counts** (k_plus_2) or counts per config (k_plus_4).

**Ground truth**: 48 crystal structures from `pdb_minimized/`, with epitope defined as antigen residues with any heavy atom within 5 Angstrom of any antibody heavy atom.

---

## 2. k_plus_2: K+ vs B1 (Head-to-Head, 47 Complexes)

### 2.1 DockQ Bar Chart (dockq_bar.png)

**What it shows**: Mean DockQ (AB+AC interface average) for B1 and K+, with standard deviation error bars. Each bar represents the mean of the best-by-confidence model across all 47 complexes.

**Results**: B1 mean = 0.179, K+ mean = 0.182. The difference is negligible (+0.003 in favor of K+). Both methods have large standard deviations (~0.27), indicating enormous per-complex variability. The error bars nearly span the full range from 0 to 0.45.

**Conclusion**: K+ and B1 produce statistically indistinguishable docking quality at the antibody-antigen interface. A paired Wilcoxon signed-rank test confirms this (p = 0.805). The beta-scaling did not degrade docking quality, but it also did not improve it on average.

### 2.2 Epitope F1 Bar Chart (epitope_f1_bar.png)

**What it shows**: Mean epitope F1 score for B1 and K+, with standard deviation error bars. Epitope is extracted from each method's best-by-confidence predicted structure using a 5 Angstrom heavy-atom distance cutoff.

**Results**: B1 mean F1 = 0.345, K+ mean F1 = 0.363. K+ shows a small advantage (+0.018), but with very large standard deviations (~0.34). The Wilcoxon test gives p = 0.464 (not significant). K+ beats B1 on 14/47 complexes, B1 beats K+ on 23/47, and 10 are tied.

**Conclusion**: K+ does not significantly improve epitope prediction over vanilla Boltz2. Both methods fall far short of the design targets (Precision > 0.7, Recall > 0.6, F1 > 0.65). The mean F1 of ~0.35 for both methods is modest, reflecting the inherent difficulty of the task. Many complexes (especially the HLA-type targets) have F1 = 0.0 for both methods.

### 2.3 K+ vs B1 Scatter Plot (kplus_vs_b1_dockq.png)

**What it shows**: Per-complex DockQ (AB+AC) comparison. Each point is one of 47 complexes, with B1 on the x-axis and K+ on the y-axis. The dashed diagonal is the y=x line: points above it indicate K+ superiority, points below indicate B1 superiority.

**Results**: The scatter reveals a striking pattern:

- **Diagonal cluster (majority)**: Most complexes where both methods produce good or poor results fall close to the y=x line (e.g., 8GQ1_HLC at 0.86/0.86, 8CDD_EDB at 0.73/0.71, 7TRH_HBG at 0.49/0.46). On these cases, K+ and B1 are interchangeable.

- **Major K+ wins (above diagonal)**:
  - 8FAH_HLA: B1=0.020, K+=0.610 (delta +0.590)
  - 8BLQ_ECD: B1=0.046, K+=0.598 (delta +0.552)
  - 8E2U_HLA: B1=0.008, K+=0.495 (delta +0.487)
  - 8HGM_CDB: B1=0.228, K+=0.695 (delta +0.466)

- **Major K+ losses (below diagonal)**:
  - 8CDE_DCB: B1=0.744, K+=0.056 (delta -0.688)
  - 8J1T_HKF: B1=0.715, K+=0.027 (delta -0.688)
  - 7ZOZ_HLA: B1=0.644, K+=0.083 (delta -0.561)

- **Bottom-left cluster**: A dense cluster of ~25 complexes where both methods have DockQ < 0.1 (effectively incorrect). These are the "hard" cases that neither method can solve.

**Conclusion**: This is the most informative plot. K+ does not uniformly improve or degrade docking. Instead, it **redistributes successes and failures**: there are a few complexes where K+ dramatically outperforms B1 (converting incorrect predictions to medium-quality ones), and a roughly equal number where K+ catastrophically fails on complexes B1 handles well. The net effect cancels out. The beta-scaling appears to be altering which binding mode the diffusion process converges to, sometimes for better, sometimes for worse.

### 2.4 DockQ Quality Distribution (dockq_distribution.png)

**What it shows**: Number of complexes falling into each DockQ quality category (Incorrect < 0.23, Acceptable 0.23-0.49, Medium 0.49-0.80, High >= 0.80) for B1 and K+. Based on the global DockQ score (average of all 3 interfaces: AB, AC, BC).

**Results**:
| Category | B1 | K+ |
|----------|----|----|
| Incorrect | 17 | 16 |
| Acceptable | 20 | 20 |
| Medium | 10 | 11 |
| High | 0 | 0 |

K+ shifts one complex from Incorrect to Medium compared to B1. Neither method achieves any High-quality predictions. Both methods have ~36% Incorrect and ~43% Acceptable.

As percentages of 47 complexes:
| Category | B1 | K+ |
|----------|----|----|
| Acceptable or better | 63.8% | 66.0% |
| Medium or better | 21.3% | 23.4% |

**Conclusion**: The quality distribution is nearly identical. K+ gains one additional Medium-quality prediction but the difference is minimal.

---

## 3. k_plus_4: K+ vs All Three Baselines (B1, B2, B3)

### 3.1 DockQ Bar Chart (dockq_bar.png)

**What it shows**: Mean DockQ (AB+AC) for all four methods. For B2 and B3, the mean is computed across **all configs** (not oracle), meaning each restraint configuration is treated as a separate data point (522 for B2, 87 for B3).

**Results**: B2 mean = 0.290, B1 mean = 0.179, K+ mean = 0.182, B3 mean = 0.148. B2 leads substantially (+0.11 over B1/K+), while B3 performs worst. All methods show very large standard deviations (0.25-0.30).

**Conclusion**: When averaging across all restraint configurations, B2 (contact restraints) outperforms all other methods. This makes sense: B2 has ~11 configs per complex, and some configurations happen to use restraints near the true epitope, producing good results that raise the mean. B3 performs worst, likely because its pocket restraints are often directed at incorrect residues. K+ and B1 remain essentially tied.

**Important caveat**: B2's advantage is inflated by the multi-config design. It has 11x more "attempts" per complex, and the mean includes the lucky ones. The fair comparison is either (a) per-complex oracle (best config) or (b) per-complex average of B2 configs.

### 3.2 Epitope F1 Bar Chart (epitope_f1_bar.png)

**What it shows**: Mean epitope F1 score for all four methods, with the same averaging as above.

**Results**: B2 = 0.452, K+ = 0.363, B1 = 0.345, B3 = 0.283.

**Conclusion**: B2 again leads on epitope F1, for the same reason: contact restraints that happen to be near the true epitope produce structures with correct binding poses, which in turn yield correct epitope predictions. B3 is worst. K+ has a small edge over B1 (+0.018) but it is not significant.

### 3.3 K+ vs B1 Scatter Plot (kplus_vs_b1_dockq.png)

This plot is identical to the k_plus_2 version since K+ and B1 have the same 47 complexes in both evaluations. See analysis in Section 2.3.

### 3.4 DockQ Quality Distribution (dockq_distribution.png)

**What it shows**: Count of predictions in each DockQ quality category for all four methods. Note that B2 has 522 data points and B3 has 87, vs 47 each for B1 and K+.

**Results** (as percentages, for fair comparison):
| Category | B1 (n=47) | K+ (n=47) | B2 (n=522) | B3 (n=87) |
|----------|-----------|-----------|------------|-----------|
| Incorrect | 36.2% | 34.0% | 22.4% | 46.0% |
| Acceptable | 42.5% | 42.6% | 40.2% | 35.6% |
| Medium | 21.3% | 23.4% | 37.5% | 18.4% |
| High | 0.0% | 0.0% | 0.2% (1 config) | 0.0% |
| **Acceptable+** | **63.8%** | **66.0%** | **77.9%** | **54.0%** |

**Conclusion**: B2 has the highest fraction of Acceptable+ predictions (77.9%), benefiting from its many configurations. B3 is the weakest (54% Acceptable+). K+ and B1 are again nearly identical, with K+ having a marginally higher Acceptable+ rate.

---

## 4. Oracle Analysis (B2 and B3 best-per-complex)

Since B2 and B3 have multiple configurations per complex, a fair comparison involves taking the **best configuration per complex** (oracle) and comparing to K+.

### B2-oracle vs K+
- B2-oracle mean DockQ(AB+AC) = **0.410** vs K+ mean = **0.186**
- B2-oracle wins on **42/46 complexes** (91%)
- Wilcoxon p < 0.000001 (highly significant)
- B2 with oracle selection decisively outperforms K+. However, this requires knowing which restraint configuration is correct -- in practice, the user does not know this.

### B3-oracle vs K+
- B3-oracle mean DockQ(AB+AC) = **0.217** vs K+ mean = **0.199** (on 36 shared complexes)
- B3-oracle wins on **19/36 complexes** (53%)
- Wilcoxon p = 0.592 (not significant)
- Even with oracle selection, B3 barely edges out K+, and the difference is not significant.

---

## 5. Summary Statistics

| Metric | B1 | K+ | B2 (mean) | B2 (oracle) | B3 (mean) | B3 (oracle) |
|--------|----|----|-----------|-------------|-----------|-------------|
| DockQ AB+AC (mean) | 0.179 | 0.182 | 0.290 | 0.410 | 0.148 | 0.217 |
| DockQ global (mean) | 0.338 | 0.341 | 0.413 | -- | 0.314 | -- |
| iRMSD AB+AC (mean) | 11.46 | 11.10 | 8.11 | -- | 12.96 | -- |
| fnat AB+AC (mean) | 0.169 | 0.182 | 0.282 | -- | 0.132 | -- |
| Ab-aligned Ag RMSD (mean) | 31.0 | 29.0 | 24.9 | -- | 32.3 | -- |
| Epitope Precision (mean) | 0.357 | 0.378 | 0.477 | -- | 0.305 | -- |
| Epitope Recall (mean) | 0.345 | 0.363 | 0.441 | -- | 0.272 | -- |
| Epitope F1 (mean) | 0.345 | 0.363 | 0.452 | -- | 0.283 | -- |
| Epitope MCC (mean) | 0.269 | 0.290 | 0.391 | -- | 0.196 | -- |
| confidence_score (mean) | 0.883 | 0.884 | 0.885 | -- | 0.884 | -- |
| iptm (mean) | 0.737 | 0.750 | 0.748 | -- | 0.762 | -- |

---

## 6. Overall Conclusions

### Did K+ achieve its goal?

**No, not in its current form.** K+ was designed to discover epitope locations and potentially improve docking quality through region-specific beta-scaling. The evaluation across 47 antibody-antigen complexes shows:

1. **Epitope prediction (Axis A)**: K+ does not significantly outperform vanilla Boltz2 (B1). Mean F1 improves from 0.345 to 0.363 (p = 0.464, not significant). Both methods fall far short of the design targets (Precision > 0.7, Recall > 0.6).

2. **Docking quality (Axis B)**: K+ is statistically indistinguishable from B1 on DockQ (p = 0.805). The scatter plot reveals that K+ redistributes successes rather than adding new ones -- it dramatically improves some complexes while catastrophically degrading others.

3. **Compared to restraint methods**: B2 (contact restraints) with oracle selection strongly outperforms K+ (DockQ 0.410 vs 0.182), but B2 requires prior knowledge of the correct restraint type/location. B3 (pocket restraints) does not significantly beat K+ even with oracle selection.

### Why does K+ fail?

The scatter plot (Section 2.3) provides the key insight: K+ does not fail uniformly. It shows a **high-variance redistribution pattern** where:
- Some complexes that B1 cannot solve are solved by K+ (e.g., 8FAH_HLA, 8BLQ_ECD)
- Some complexes that B1 solves well are broken by K+ (e.g., 8CDE_DCB, 8J1T_HKF)

This suggests the beta-scaling is strong enough to steer diffusion into different binding modes, but it lacks the selectivity to consistently steer toward the correct one.

### Recommended Modifications

1. **Ensemble over beta-scaling runs**: Instead of producing a single prediction with CDR3 beta-scaling, run K+ with multiple region-emphasis patterns and select the best by confidence score. This would capture K+'s wins while mitigating its losses (similar to B2's multi-config approach). The scatter plot shows that when K+ succeeds, it can reach DockQ > 0.6, suggesting an oracle-over-K+-regions strategy could be powerful.

2. **Softer beta values**: The current beta-scaling may be too aggressive, causing the diffusion to collapse to incorrect modes on some complexes. Reducing `beta_emphasis` from 0.5 to 0.2-0.3 and `beta_deemphasis` from -0.3 to -0.1 could reduce the variance of outcomes.

3. **Confidence-based filtering**: Use the Boltz confidence metrics (iptm, iplddt) to detect when K+ has produced a poor structure and fall back to the B1 prediction. The confidence scores are similar between methods (both ~0.88), so a more discriminative confidence metric may be needed.

4. **Contact heatmap aggregation**: K+ was designed to produce epitope heatmaps by accumulating contacts across region scans. The current evaluation only tests the single best-by-confidence structure. Implementing and evaluating the full heatmap accumulation pipeline may reveal K+'s intended strength.

5. **Combine with restraints**: Use K+'s epitope heatmap to identify candidate epitope regions, then run B2-style contact restraints on those regions. This two-stage approach would combine K+'s exploration capability with B2's proven ability to improve docking when given correct restraints.

6. **Target-specific analysis**: The evaluation reveals that K+ and B1 both fail on the majority of HLA-type targets (long-chain antigens with distant epitopes). Focus K+ development on the target classes where it shows potential (shorter antigens, single-domain complexes).
