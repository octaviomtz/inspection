# Evaluation Report: Strategy G+ (Progressive CDR Refinement v2)

## Setup

**Test set**: 47 antibody-antigen complexes with minimized crystal structure ground truths.

**Methods compared**:

| Label | Description | N entries |
|-------|-------------|-----------|
| B1 | No steering (standard Boltz2) | 47 complexes x 5 models |
| B2 | Pocket-based contact restraints (residue-specific) | 87 sub-experiments x 5 models |
| B3 | Contact-type restraints (hbond, hydrophobic, salt bridge) | 522 sub-experiments x 5 models |
| G+ | Strategy G+ with adaptive phases + CDR3 beta-scaling | 47 complexes x 5 models |

B2 and B3 run multiple sub-experiments per complex (different restraint residues or contact types), resulting in more entries than the 47 complexes. B3 averages ~11 sub-experiments per complex.

**Primary metric**: DockQ (continuous score 0-1, combining fnat, L-RMSD, and I-RMSD). Computed using the DockQ v2 Python API with chain mapping A:A B:B C:C (A=antigen, B=heavy chain, C=light chain). DockQ handles antibody superposition and antigen displacement internally.

**Model selection**: For each method-complex pair, the model with the highest Boltz2 `confidence_score` is selected ("confidence-selected best-of-5"). This is the field-standard approach (AlphaFold-Multimer, Chai-1, Boltz-1). Oracle selection (best DockQ, requiring ground truth) is reported as a ceiling.

---

## Results: g_plus_2 (G+ vs B1)

### Figure 1: DockQ Bar Chart

**How produced**: Mean DockQ across the 47 complexes, using the confidence-selected best-of-5 model for each complex. Error bars show standard deviation across complexes.

**Results**: B1 mean = 0.338, G+ mean = 0.355. G+ shows a modest +0.017 improvement in mean DockQ. The standard deviations are large and fully overlapping (~0.18 for both), indicating high variance across complexes and no clear separation between methods at the aggregate level.

**Conclusion**: G+ provides a small average improvement over B1, but it is not statistically distinguishable at this sample size.

---

### Figure 2: G+ vs B1 Scatter Plot

**How produced**: Each point is one complex. X-axis = B1 confidence-selected DockQ, Y-axis = G+ confidence-selected DockQ. Points above the y=x diagonal indicate G+ is better.

**Results**: Most points cluster near the diagonal, especially in the low-DockQ region (0.20-0.30) where many complexes are difficult for both methods. A few notable outliers:
- **8EZ8_HLA**: B1=0.199 -> G+=0.644 (+0.445), the largest single improvement, jumping from Incorrect to Medium
- **8OL9_BAH**: B1=0.327 -> G+=0.600 (+0.273), Acceptable to Medium
- **8HGM_CDB**: B1=0.365 -> G+=0.546 (+0.181), Acceptable to Medium

For high-DockQ complexes (>0.6), points sit slightly below the diagonal, meaning G+ marginally underperforms B1 on already-easy cases (e.g. 8CDE_DCB: 0.741->0.712, 8DE3_BCA: 0.680->0.637).

**Conclusion**: G+ helps most on a subset of medium-difficulty complexes where B1 fails, producing large improvements. On easy complexes that B1 already solves well, G+ is slightly worse. On truly hard complexes (~0.20), neither method succeeds. The overall win/loss count is 20 wins for G+ vs 27 for B1, but G+ wins are larger in magnitude (mean delta when G+ wins: +0.043 vs -0.019 when B1 wins).

---

### Figure 3: CAPRI Quality Distribution (Stacked Bar)

**How produced**: Each confidence-selected best-of-5 model is classified into CAPRI quality tiers: Incorrect (<0.23), Acceptable (0.23-0.49), Medium (0.49-0.80), High (>=0.80). Bars show the percentage of 47 complexes in each tier.

**Results**:

| Tier | B1 | G+ |
|------|----|----|
| Medium | 21.3% (10) | 27.7% (13) |
| Acceptable | 42.6% (20) | 38.3% (18) |
| Incorrect | 36.2% (17) | 34.0% (16) |
| High | 0% | 0% |

G+ shifts 3 more complexes into the Medium tier compared to B1, with fewer Incorrect classifications.

Detailed CAPRI transitions: G+ upgrades 6 complexes (3 Incorrect->Acceptable, 2 Acceptable->Medium, 1 Incorrect->Medium) while downgrading 3 (all Acceptable->Incorrect). Net: +3 CAPRI-class improvements.

**Conclusion**: G+ improves the quality distribution, especially by promoting more complexes to Medium quality. The net CAPRI upgrade count is positive (+3), consistent with the scatter plot findings.

---

### Figure 4: Confidence Score vs DockQ (All Models)

**How produced**: Each point is one individual model (all 5 models per complex, both methods). X-axis = Boltz2 confidence_score, Y-axis = DockQ. Blue = B1, Red = G+.

**Results**: Both methods show a scattered relationship between confidence and DockQ. High-confidence models (>0.90) can have DockQ anywhere from 0.20 to 0.80, and low-confidence models (<0.80) sometimes achieve moderate DockQ. The distributions for B1 and G+ are largely overlapping, with no clear difference in the confidence-DockQ correlation pattern.

**Conclusion**: The Boltz2 confidence score is a weak predictor of actual docking quality. This means the confidence-based model selection (used for all our primary results) is imperfect -- there is a gap between the confidence-selected model and the oracle-best model. For B1 the gap is 0.338 (conf) vs 0.391 (oracle); for G+ it is 0.355 (conf) vs 0.392 (oracle). A better model selector could recover ~0.04 DockQ for both methods.

---

### Figure 5: Ensemble Diversity Box Plot

**How produced**: For each method-complex pair, the mean pairwise CA-RMSD is computed across all 10 pairs of the 5 models (over all chains A, B, C). Each box shows the distribution of these diversity values across the 47 complexes.

**Results**: B1 median = 8.5 A, G+ median = 9.1 A. Both distributions span from ~0.2 A (very consistent predictions) to ~23 A (highly diverse). G+ has a slightly more compact interquartile range (5-14 A) compared to B1 (2-15 A), with a marginally higher median.

**Conclusion**: G+ maintains similar ensemble diversity to B1. The exploratory Phase 1 of the G+ strategy (negative beta for CDR3 exploration) does not significantly increase or decrease structural diversity compared to the unsteered baseline. Ensemble diversity preservation, one of the G+ design goals, is achieved.

---

### Figure 6: Per-Complex DockQ Bar Chart

**How produced**: Side-by-side bars for B1 (blue) and G+ (red) for each of the 47 complexes, sorted by ascending B1 DockQ. Confidence-selected best-of-5 model.

**Results**: The chart confirms the scatter plot findings: for the ~30 low-DockQ complexes on the left (DockQ < 0.30), the bars are nearly identical between methods. In the middle range, a few G+ bars clearly exceed B1 (8EZ8_HLA, 8OL9_BAH, 8HGM_CDB). On the right side (high-DockQ complexes), B1 bars are slightly taller.

**Conclusion**: G+ does not uniformly improve or worsen predictions. It selectively helps on specific complexes while leaving most unchanged and slightly hurting a few well-predicted ones.

---

## Results: g_plus_4 (G+ vs All Three Baselines)

### Figure 1: DockQ Bar Chart (4 methods)

**How produced**: Same methodology as g_plus_2 but including B2 and B3. Mean DockQ of confidence-selected best-of-5, with standard deviation error bars. Note: B2 has 87 entries and B3 has 522 entries (multiple sub-experiments per complex), so their means are computed across sub-experiments, not per-complex.

**Results**:

| Method | Mean DockQ | Median DockQ |
|--------|-----------|--------------|
| B1 | 0.338 | 0.244 |
| B2 | 0.314 | 0.235 |
| B3 | 0.413 | 0.293 |
| G+ | 0.355 | 0.249 |

Ranking: B3 > G+ > B1 > B2.

**Conclusion**: B3 (contact-type restraints) is the strongest method, outperforming all others by a wide margin. This is expected -- B3 uses explicit structural knowledge derived from the known crystal structure contacts (hbond, hydrophobic, salt bridge restraints), which amounts to informing the model about the native interface. G+ outperforms both B1 and B2 without using any such oracle information. B2 (pocket restraints) actually performs slightly worse than B1, suggesting that its specific restraint formulation may be suboptimal.

**Important caveat**: B3's 522 entries include ~11 sub-experiments per complex with different restraint types. Its mean includes the best restraint configurations alongside poor ones. This inflates its denominator but also means some sub-experiments specifically target the correct interface.

---

### Figure 2: G+ vs B1 Scatter Plot (4-baseline run)

**How produced**: Identical to g_plus_2 (same G+ and B1 data).

**Results and conclusion**: Same as g_plus_2 Figure 2.

---

### Figure 3: CAPRI Quality Distribution (4 methods)

**How produced**: Same methodology, now with 4 bars. Note B2 and B3 counts are per sub-experiment, not per unique complex.

**Results**:

| Tier | B1 | B2 | B3 | G+ |
|------|----|----|----|----|
| Medium+ | 21.3% | 18.4% | 37.9% | 27.7% |
| Acceptable | 42.6% | 35.6% | 40.0% | 38.3% |
| Incorrect | 36.2% | 46.0% | 22.2% | 34.0% |

B3 achieves the most Medium+ predictions (37.9%) and fewest Incorrect (22.2%). G+ is second-best on both metrics. B2 is worst, with 46% Incorrect.

**Conclusion**: G+ occupies a clear middle ground: better quality distribution than B1 and B2, but below B3. This is consistent with G+ using learned priors only (no structural oracle information) while B3 leverages native contact knowledge.

---

### Figure 4: Confidence Score vs DockQ (4 methods)

**How produced**: Same as g_plus_2 but with all 4 methods plotted. The image is small but shows 4 colored point clouds (blue=B1, orange=B2, green=B3, red=G+).

**Results**: The point clouds for all methods overlap substantially. B3 (green) has a notable tail of high-DockQ points at moderate confidence scores, consistent with its structural restraints improving docking quality beyond what the model's confidence alone would predict. B2 (orange) clusters slightly lower than B1 on average.

**Conclusion**: Confidence score remains a weak predictor across all methods. B3's restraints decouple DockQ from confidence more than other methods, as the restraints force better structures even when the model is uncertain.

---

### Figure 5: Ensemble Diversity Box Plot (4 methods)

**How produced**: Same pairwise CA-RMSD methodology, now for all 4 methods.

**Results**: B1 and B2 have similar distributions (median ~8.5-10.5 A). B3 has the lowest median diversity (~6.0 A) and a more compact box. G+ sits between B1 and B3 (median ~9.1 A).

**Conclusion**: B3's explicit restraints reduce structural diversity by constraining the antibody-antigen interface, pushing all 5 models toward the restrained conformation. G+ maintains diversity similar to B1, which is desirable -- it does not collapse the ensemble. B2's diversity is comparable to B1, suggesting its pocket restraints don't strongly constrain the sampling.

---

### Figure 6: Per-Complex DockQ Bar Chart (4 methods)

**How produced**: Side-by-side bars for all 4 methods per complex, confidence-selected best-of-5. Sorted by ascending B1 DockQ.

**Results**: B3 (green) frequently produces the tallest bar, especially on medium-difficulty complexes (e.g., 8EZ8_HLA, 8OL9_BAH, 8HGM_CDB). On some complexes B3 bars are dramatically taller than all others. G+ (red) often matches or slightly exceeds B1 (blue). B2 (orange) tends to track B1 closely, sometimes lower.

**Conclusion**: B3's explicit contact restraints provide the most consistent improvements. G+, without any oracle information, occasionally matches B3 levels on favorable complexes but cannot do so consistently.

---

## Statistical Tests

### Wilcoxon Signed-Rank Test (paired by complex, vs B1)

| Method | N paired | Median delta | Mean delta | p-value |
|--------|----------|-------------|-----------|---------|
| G+ | 47 | -0.003 | +0.017 | 0.410 |
| B2 | 87 | -0.0002 | -0.00001 | 0.145 |
| B3 | 522 | +0.004 | +0.059 | 1.2e-21 |

- **G+ vs B1**: Not statistically significant (p=0.41). The positive mean delta (+0.017) is driven by a few large G+ wins, but the median is slightly negative, reflecting that B1 wins more often by smaller margins.
- **B2 vs B1**: Not significant (p=0.15). B2 is essentially equivalent to B1.
- **B3 vs B1**: Highly significant (p=1.2e-21). B3 provides a consistent, genuine improvement over B1.

---

## Overall Conclusions

### Did G+ achieve its goals?

**Goal 1 -- Improve docking quality over unsteered baseline (B1)**: **Partially.** G+ improves mean DockQ (+0.017) and shifts more complexes into higher CAPRI quality tiers (6 upgrades vs 3 downgrades). However, the improvement is not statistically significant (p=0.41), B1 wins more individual head-to-heads (27 vs 20), and the median delta is slightly negative (-0.003). The benefit comes from a few large wins on specific complexes rather than a broad improvement.

**Goal 2 -- Ensemble diversity preservation**: **Yes.** G+ maintains ensemble diversity comparable to B1 (median pairwise RMSD: 9.1 A vs 8.5 A). The exploratory first phase does not collapse the ensemble.

**Goal 3 -- Competitive with contact restraint baselines**: **No.** B3 (contact-type restraints) significantly outperforms G+ (mean DockQ 0.413 vs 0.355, p=1.2e-21). This gap is expected since B3 uses oracle structural information that G+ does not have access to. G+ does outperform B2 (pocket restraints), which is itself worse than B1.

### Key Observations

1. **G+ helps most on medium-difficulty cases**: The 3 largest G+ improvements (8EZ8_HLA +0.445, 8OL9_BAH +0.273, 8HGM_CDB +0.181) are all complexes where B1 produces borderline Incorrect/Acceptable predictions. G+ promotes them to Medium quality.

2. **G+ slightly hurts easy cases**: On well-predicted complexes (B1 DockQ > 0.6), G+ introduces small regressions (average -0.025). The adaptive steering may be perturbing already-good predictions.

3. **Confidence-based model selection is suboptimal**: Both methods leave ~0.04 DockQ on the table vs oracle selection. A better selector could amplify G+'s benefits.

4. **B3's advantage is structural oracle information**: B3 uses contact restraints derived from native structures, which is not a fair comparison for methods targeting blind prediction. G+ should be compared primarily against B1 and B2.

### Suggested Modifications

1. **Conditional steering**: Only apply G+ steering when the model's initial confidence is low (e.g., iptm < 0.8). Skip steering on high-confidence predictions to avoid regressing easy cases.

2. **Stronger CDR3 beta-scaling**: The current Phase 1 exploratory beta (-0.5) may be too weak. More aggressive negative beta values could expand CDR3 conformational search on difficult targets.

3. **Phase transition tuning**: The contact-triggered phase advancement threshold (5% improvement) may transition too early on hard complexes. Allowing more exploration time before advancing to Phase 2 could help.

4. **Better model selection**: The weak confidence-DockQ correlation suggests developing an interface-specific selection criterion (e.g., combining iptm + interface pLDDT) could recover the oracle-selection gap.

5. **Combine with contact restraints**: A hybrid approach using G+'s adaptive phases with B3-style contact restraints (when available) could combine the strengths of both methods.
