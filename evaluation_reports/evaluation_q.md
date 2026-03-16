# Method Q (Iterative Epitope Refinement) -- Evaluation Report

## Overview

Method Q implements a 3-round iterative prediction approach for antibody-antigen docking:
- **Round 1:** Weak guidance (scale 0.1) over the full antigen surface with exploratory beta-noise
- **Round 2:** Medium guidance (scale 0.5) focused on hotspot residues identified from Round 1
- **Round 3:** Full guidance (scale 1.0) on refined hotspots

The evaluation was performed on 47 antibody-antigen complexes across two axes:
- **Axis A (Primary):** Epitope prediction accuracy (Precision, Recall, F1, MCC)
- **Axis B (Secondary):** Docking/interface quality (DockQ, CAPRI classification)

Two evaluation runs were conducted:
- **q_plus_2:** Q vs B1 (vanilla Boltz2) only
- **q_plus_4:** Q vs all three baselines (B1, B2 contact restraints, B3 pocket restraints)

**Model selection:** For all plots and metrics, the **best-by-confidence** model was used (the model with the highest `confidence_score` out of 5 predictions per complex). This mimics real-world usage where ground truth is unavailable. For B2 and B3 (which have multiple restraint configurations per complex), both the **best** configuration (oracle over configs) and **mean** across configurations are reported.

---

## Plot-by-Plot Analysis

### 1. DockQ Box Plot -- Antigen-Antibody Interfaces

**How produced:** Box plot of the AB+AC average DockQ score (mean of antigen-heavy and antigen-light interface DockQ) for the confidence-selected best model per complex. CAPRI quality thresholds are shown as horizontal dashed lines (Acceptable=0.23, Medium=0.49, High=0.80).

**q_plus_2 (Q vs B1):**
Both B1 and Q show similar distributions with medians near 0.04--0.05 (below the Acceptable threshold). The interquartile ranges are comparable though Q's box extends slightly higher (Q75 ~0.37 vs B1's Q75 ~0.21). Both have a long tail of outliers reaching up to ~0.85. The distributions are dominated by low-DockQ complexes, reflecting the inherent difficulty of blind antibody-antigen docking.

**q_plus_4 (all methods):**
B2_best is the clear winner with a median of ~0.41 and upper quartile reaching ~0.69, substantially above the Acceptable threshold. B2_mean, B3_best, and B3_mean all show intermediate performance. Q and B1 sit at the bottom with nearly identical median values (~0.04--0.05). Q's upper quartile (~0.37) is slightly higher than B1's (~0.21) but the improvement is modest and driven by a few complexes.

**Conclusion:** Q does not meaningfully improve docking quality over the vanilla baseline B1 on the AB+AC interface metric. Both are far behind B2_best, which has the advantage of ground-truth-informed contact restraints.

---

### 2. Epitope F1 Box Plot

**How produced:** Box plot of the epitope prediction F1 score for the confidence-selected best model per complex. The predicted epitope is defined as the set of antigen residues within 5A of any antibody atom in the predicted structure, compared against the same definition applied to the crystal structure.

**q_plus_2 (Q vs B1):**
Q shows a higher median F1 (~0.37) compared to B1 (~0.23). Q's interquartile range (IQR ~0.08--0.65) is narrower than B1's (IQR ~0.00--0.67), and Q's lower whisker sits higher, indicating fewer zero-F1 failures. Both reach similar maximum values (~0.91).

**q_plus_4 (all methods):**
B2_best leads with a median F1 of ~0.67 and a tight IQR. B2_mean and B3_best follow. Q (median ~0.37) sits above B1 (median ~0.23), B3_mean (~0.22), and roughly comparable to B3_best (~0.28). The improvement of Q over B1 in epitope F1 is visible but the gap is not dramatic.

**Conclusion:** Q shows a modest improvement in epitope prediction over B1 (mean F1: 0.381 vs 0.352, median F1: 0.373 vs 0.233), but this is not statistically significant (p=0.75, Wilcoxon). Q falls well short of B2_best (mean F1: 0.583).

---

### 3. DockQ Scatter -- Q vs B1

**How produced:** Scatter plot with one point per complex. X-axis = B1 DockQ (AB+AC avg), Y-axis = Q DockQ (AB+AC avg) for the confidence-selected model. The dashed diagonal line represents y=x (no change). Points above the line indicate Q improved over B1; points below indicate degradation.

**Analysis:**
The majority of points cluster in the bottom-left corner (DockQ < 0.1 for both), where both methods essentially fail -- these are the hard cases where neither method finds the correct binding site. For complexes where B1 already performs well (DockQ > 0.4), Q tends to track closely along the diagonal, indicating no improvement. A few notable exceptions:
- **8EZ8_HLA:** Large improvement (B1 ~0.01, Q ~0.41) -- Q found the binding site when B1 completely missed it
- **8OL9_BAH:** Strong improvement (B1 ~0.19, Q ~0.56)
- **7TRH_HBG:** Notable improvement (B1 ~0.49, Q ~0.67)
- **7ZOZ_HLA:** Degradation (B1 ~0.64, Q ~0.45)
- **8BLQ_ABD:** Slight degradation (B1 ~0.66, Q ~0.62)

**Conclusion:** Q produces mixed results at the per-complex level. It achieves large improvements on ~3--4 complexes but slight degradations on others. Most complexes are unchanged. The overall effect averages to near-zero.

---

### 4. Epitope F1 Scatter -- Q vs B1

**How produced:** Scatter plot with one point per complex. X-axis = B1 epitope F1, Y-axis = Q epitope F1 for the confidence-selected model. Diagonal = no change.

**Analysis:**
The scatter shows a general trend along the diagonal with moderate dispersion. Key observations:
- **Left edge (B1 F1 ~0):** Several points where Q achieves non-zero F1 (0.05--0.85) while B1 has F1=0. These are cases where Q identified at least some epitope residues in complexes where B1 found nothing. This is the most promising pattern -- Q's iterative refinement helps discover binding sites on otherwise intractable targets.
- **Right side (B1 F1 > 0.7):** Points generally cluster near or slightly below the diagonal, suggesting Q doesn't harm epitope prediction where B1 is already successful.
- **Notable outlier:** One point at ~(0.65, 0.08) shows a case where Q badly degraded an otherwise good epitope prediction (likely 8CDE_DCB, where DockQ is high but both methods have low epitope F1 around 0.17--0.18).

**Conclusion:** Q shows a pattern of helping on "hard zero" cases where B1 completely fails, while roughly matching B1 on easier targets. The overall population effect is positive but small (18 wins, 18 losses, 11 ties).

---

### 5. Per-Complex DockQ Delta -- Q vs B1

**How produced:** Bar chart showing DockQ(Q) - DockQ(B1) for each of the 47 complexes, sorted from most negative to most positive. Green bars = Q improved, red bars = Q degraded.

**Analysis:**
The distribution is asymmetric in an interesting way:
- **Negative deltas (left):** ~10 complexes show degradation, with the worst being 7ZOZ_HLA (delta ~ -0.19) and 8BLQ_ABD (delta ~ -0.04). Most negative deltas are small (< 0.05).
- **Near-zero (middle):** ~30 complexes have near-zero delta, meaning Q neither helped nor hurt. This is the dominant pattern.
- **Positive deltas (right):** ~4--5 complexes show significant improvement: 8EZ8_HLA (+0.40), 8OL9_BAH (+0.37), 7TRH_HBG (+0.18), 8TRS_AGD (+0.08), 8DTK_CBA (+0.06).

The positive tail is larger in magnitude than the negative tail -- the biggest wins (+0.40) exceed the biggest losses (-0.19).

**Conclusion:** Q produces a few dramatic improvements on select complexes while causing smaller degradations on others. The net effect is slightly positive (mean delta +0.018) but not statistically significant (Wilcoxon p=0.79). The feature shows potential but is not consistently beneficial.

---

### 6. CAPRI Quality Distribution -- Stacked Bar

**How produced:** Stacked bar chart showing the percentage of complexes in each CAPRI quality category (High >= 0.80, Medium >= 0.49, Acceptable >= 0.23, Incorrect < 0.23) based on the AB+AC avg DockQ of the confidence-selected model.

**q_plus_2 (Q vs B1):**
Both methods are dominated by "Incorrect" (~70--75%). Q shows a slight improvement: ~30% at Acceptable or better vs ~25% for B1. The Medium+High fractions are nearly identical (~21% for Q vs ~19% for B1). Both have ~2% High quality.

**q_plus_4 (all methods):**
B2_best stands out with ~37% Incorrect, ~20% Acceptable, ~37% Medium, and ~7% High. B1 and Q are nearly indistinguishable (~75% Incorrect). B3_best and B3_mean are also similar to B1. This highlights that Q's iterative approach without ground-truth restraints cannot match the performance of methods that incorporate known contact information.

**Conclusion:** Q marginally shifts a few complexes from Incorrect to Acceptable compared to B1, but the CAPRI distributions are largely indistinguishable. The dominant outcome for both Q and B1 remains Incorrect (~70+%).

---

### 7. Confidence Score vs DockQ Scatter

**How produced:** Scatter plot with confidence score (x-axis) vs Global DockQ (y-axis) for the best-by-confidence model per complex, colored by method.

**q_plus_2 (Q vs B1):**
Both methods show a weak positive correlation between confidence and DockQ, with substantial scatter. B1 (blue) and Q (red) points largely overlap. Confidence scores range from ~0.78--0.98 while DockQ varies from ~0.15 to ~0.79. Neither method shows a clear advantage in confidence calibration -- high confidence doesn't reliably predict high DockQ for either.

**q_plus_4 (all methods):**
The cloud becomes denser with all 6 method variants. B2_best (dark orange) points tend to appear at higher DockQ values across the confidence range, but there's no method-specific pattern in confidence. All methods cluster in a similar confidence range (0.78--0.98). The lack of clear clustering by color in the high-DockQ region suggests that confidence is an imperfect proxy for actual docking quality across all methods.

**Conclusion:** Q does not improve the confidence-DockQ correlation compared to B1. The Boltz confidence score remains a noisy predictor of actual docking quality regardless of method.

---

## Numerical Summary

### Docking Quality (AB+AC avg DockQ, confidence-selected)

| Method | Mean | Median | % Acceptable+ | % Medium+ | % High |
|--------|------|--------|---------------|-----------|--------|
| B1 | 0.179 | 0.038 | 23.4% | 19.1% | 2.1% |
| B2_best | 0.408 | 0.414 | 63.0% | 43.5% | 6.5% |
| B2_mean | 0.269 | 0.122 | 45.7% | 23.9% | 2.2% |
| B3_best | 0.217 | 0.041 | 30.6% | 27.8% | 2.8% |
| B3_mean | 0.202 | 0.038 | 27.8% | 27.8% | 2.8% |
| **Q** | **0.198** | **0.049** | **29.8%** | **21.3%** | **2.1%** |

### Epitope Prediction (confidence-selected)

| Method | Mean F1 | Median F1 | Mean Precision | Mean Recall | Mean MCC |
|--------|---------|-----------|---------------|-------------|----------|
| B1 | 0.352 | 0.233 | 0.397 | 0.325 | 0.269 |
| B2_best | 0.583 | 0.669 | 0.660 | 0.535 | 0.540 |
| B2_mean | 0.450 | 0.444 | 0.509 | 0.416 | 0.385 |
| B3_best | 0.369 | 0.281 | 0.428 | 0.330 | 0.288 |
| B3_mean | 0.348 | 0.219 | 0.406 | 0.311 | 0.266 |
| **Q** | **0.381** | **0.373** | **0.444** | **0.347** | **0.305** |

### Statistical Tests (Wilcoxon signed-rank vs B1)

| Method | Metric | p-value | Mean Delta | Wins/Losses/Ties |
|--------|--------|---------|-----------|-------------------|
| Q | Global DockQ | 0.416 | +0.016 | 22/25/0 |
| Q | AB+AC DockQ | 0.789 | +0.018 | 22/25/0 |
| Q | Epitope F1 | 0.753 | +0.029 | 18/18/11 |
| Q | Epitope MCC | 0.782 | +0.036 | 23/22/2 |
| B2_best | Global DockQ | **2.5e-12** | +0.159 | 45/1/0 |
| B2_best | AB+AC DockQ | **4.2e-11** | +0.225 | 43/3/0 |
| B2_best | Epitope F1 | **5.4e-06** | +0.243 | 33/6/7 |

---

## Overall Conclusions

### Did Method Q achieve its goal?

**Partially, but not convincingly.** Method Q aimed to discover epitope locations on unknown antigens through iterative refinement. The results show:

1. **Marginal improvements in epitope prediction:** Q raises mean epitope F1 from 0.352 to 0.381 (+8.2%) and median F1 from 0.233 to 0.373 (+60% relative), suggesting that the iterative hotspot refinement does shift the distribution in the right direction. However, none of these improvements are statistically significant (all p > 0.4).

2. **No meaningful improvement in docking quality:** DockQ scores are essentially unchanged (mean +0.018, median +0.011). The CAPRI distributions are nearly identical. Q wins and losses approximately balance out across complexes.

3. **Sporadic large improvements on individual complexes:** Q produces dramatic improvements on 3--4 complexes (e.g., 8EZ8_HLA: DockQ 0.01 -> 0.41, 8OL9_BAH: 0.19 -> 0.56) while slightly degrading a few others. This suggests the iterative mechanism can work in some cases but is unreliable.

4. **Far behind informed baselines:** B2_best (contact restraints with ground-truth knowledge) vastly outperforms Q on every metric, with highly significant p-values. This is expected since B2 uses oracle information, but it shows how much room for improvement remains.

### Why the limited improvement?

The iterative refinement approach has a fundamental bootstrapping problem: Round 1 uses weak guidance (0.1x) over the full antigen, but on most complexes the vanilla Boltz2 prediction already places the antibody far from the true binding site. The hotspot identification in Round 2 then focuses on contacts from an incorrect pose, potentially reinforcing the wrong epitope. Without a reliable initial signal, the iterative narrowing converges to a local minimum rather than the true epitope.

### Suggested Modifications

1. **Increase Round 1 diversity:** Use more diffusion samples (e.g., 10--20 instead of 5) with higher beta-noise and weaker guidance to genuinely explore the antigen surface. The current 5-sample setup likely doesn't provide enough diversity for the hotspot identification to be meaningful.

2. **Ensemble-based hotspot selection:** Instead of selecting hotspots from a single round, run multiple independent Round 1s and aggregate the contact maps. This would reduce the noise in hotspot identification.

3. **Contact-frequency thresholding instead of quantile:** Rather than taking the top 20% of antigen residues by contact frequency, use an absolute contact threshold. If no residues meet the threshold, fall back to full-antigen guidance rather than forcing narrowing.

4. **Incorporate confidence-aware filtering:** Weight Round 1 contacts by the per-residue pLDDT or iPTM of the generating model. Low-confidence models should contribute less to hotspot identification.

5. **Adaptive round scheduling:** If Round 1 contacts are dispersed uniformly (no clear hotspot), skip Round 2 refinement and go directly to full guidance. The iterative narrowing only helps when there's a genuine signal to refine.

6. **Hybrid with contact restraints:** Combine Q's iterative approach with lightweight contact restraints (e.g., from predicted contacts or evolutionary coupling) to provide a better starting signal for Round 1.
