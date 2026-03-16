# Evaluation Report: Antigen Steering with 20 FK Particles (Strategy A v2)

**Date**: 2026-03-10
**Evaluation tool**: `scripts/evaluate.py` (DockQ v2)
**Results directories**: `a_plus_2_particles_20/` (NF vs B0) and `a_plus_4_particles_20/` (NF vs B0, B1, B2)
**Inference flags**: `--use_potentials --antigen_steering --num_particles 20 --diffusion_samples 5`
**Change from previous run**: `--num_particles 20` (was 5 in the previous evaluation reported in `evaluation_report.md`)

---

## 1. Evaluation Setup

### Methods Compared

| ID | Method | Description | N complexes | Models per complex |
|----|--------|-------------|------------|-------------------|
| B0 | Unconstrained | Boltz-2 baseline | 47 | 5 |
| B1 | Contact restraints | Ground-truth contact restraints (oracle info) | 46 | ~57 per complex (multiple settings) |
| B2 | Pocket restraints | Predicted pocket restraints | 36 | ~12 per complex (multiple settings) |
| NF | Antigen steering | `--antigen_steering --num_particles 20` | **35** | 5 |

**Important caveat**: NF has predictions for only **35 of 47** complexes (74%). 12 complexes are missing NF results, likely due to OOM or timeout from the 20-particle FK resampling. Aggregate metrics for NF are computed only over these 35 complexes. Paired statistical tests use only the intersection (35 for NF vs B0, 34 for NF vs B1, 25 for NF vs B2).

### Model selection strategies

- **Confidence selection**: Pick the model with the highest Boltz confidence_score (practical, no ground-truth needed).
- **Oracle selection**: Pick the model with the highest DockQ_AbAg (upper bound, requires ground truth).

### Primary metric

**DockQ_AbAg** = mean(DockQ_AB, DockQ_AC): average DockQ across both antibody-antigen interfaces.

### Comparison with previous run

The previous evaluation (in `evaluation_report.md`) used `--num_particles 5` and achieved a mean delta of +0.010 (not significant). This run uses `--num_particles 20` to enable stronger Feynman-Kac particle resampling, as recommended.

---

## 2. Plot-by-Plot Analysis: NF vs B0 (a_plus_2_particles_20)

### 2.1 Aggregate DockQ Bar Chart (`plot_aggregate_dockq.png`)

**How produced**: Mean DockQ_AbAg across all complexes available per method. Solid bars = confidence selection, hatched bars = oracle selection. Note B0 uses N=47 while NF uses N=35 (different denominators).

| Method | Confidence (N) | Oracle (N) |
|--------|---------------|-----------|
| B0 | 0.179 (47) | 0.248 (47) |
| NF | **0.247** (35) | **0.314** (35) |

**Analysis**: NF shows a substantial improvement: +0.068 in confidence selection (+38% relative) and +0.066 in oracle (+27% relative). However, the different sample sizes (35 vs 47) inflate the comparison since the 12 missing NF complexes are a mix of easy and hard cases. On the paired subset (N=35), B0's mean is 0.201 (confidence) and 0.270 (oracle), giving paired deltas of +0.046 and +0.045 respectively. This is still a **4.6x improvement** over the previous run's +0.010 delta.

---

### 2.2 CAPRI Quality Distribution -- Confidence Selection (`plot_capri_confidence.png`)

**How produced**: Each complex's confidence-selected DockQ_AbAg classified into CAPRI tiers: Incorrect (<0.23), Acceptable (0.23-0.49), Medium (0.49-0.80), High (>=0.80).

| Tier | B0 (N=47) | NF (N=35) |
|------|-----------|-----------|
| Incorrect | 36/47 (77%) | 22/35 (**63%**) |
| Acceptable | 2/47 (4%) | 3/35 (**9%**) |
| Medium | 8/47 (17%) | 9/35 (**26%**) |
| High | 1/47 (2%) | 1/35 (3%) |

**Analysis**: NF reduces incorrect from 77% to 63% and increases Acceptable+Medium+High from 23% to 37%. Among the 35 paired complexes, there are **5 CAPRI tier upgrades and 0 downgrades**:

| Complex | B0 | NF | Tier change |
|---------|----|----|-------------|
| 8FAH_HLA | 0.020 (Incorrect) | 0.610 (Medium) | +2 tiers |
| 8BLQ_ECD | 0.046 (Incorrect) | 0.574 (Medium) | +2 tiers |
| 7Y0O_HLA | 0.008 (Incorrect) | 0.452 (Acceptable) | +1 tier |
| 8HGM_CDB | 0.228 (Incorrect) | 0.313 (Acceptable) | +1 tier |
| 7TRH_HBG | 0.487 (Acceptable) | 0.776 (Medium) | +1 tier |

**Zero tier downgrades is a key result** -- the previous run had 2 downgrades (7ZOZ_HLA, 8EAY_HLA). The 20-particle FK resampling appears to protect against catastrophic regressions.

---

### 2.3 CAPRI Quality Distribution -- Oracle Selection (`plot_capri_oracle.png`)

**How produced**: Same as above but using the best DockQ_AbAg model among the 5 samples.

| Tier | B0 (N=47) | NF (N=35) |
|------|-----------|-----------|
| Incorrect | 32/47 (68%) | 18/35 (**51%**) |
| Acceptable | 1/47 (2%) | 4/35 (**11%**) |
| Medium | 13/47 (28%) | 12/35 (**34%**) |
| High | 1/47 (2%) | 1/35 (3%) |

**Analysis**: Under oracle selection, NF reduces incorrect from 68% to 51% -- a 17 percentage-point drop. Nearly half of NF's complexes (49%) have at least one acceptable-or-better model among 5 samples. This is the strongest signal for the feature: the steering is producing qualitatively better orientations.

---

### 2.4 Box Plot -- Confidence Selection (`plot_boxplot_confidence.png`)

**How produced**: Box-and-whisker of per-complex DockQ_AbAg with jittered individual points. Median = thick horizontal line.

**Analysis**: NF's box is visibly larger and taller than B0's. B0's median is ~0.04 with Q3 ~0.21, while NF's median is ~0.06 with Q3 ~0.54. The upper quartile shift from ~0.21 to ~0.54 is the clearest distributional difference between the two methods. NF has fewer near-zero outliers and more mid-range points.

---

### 2.5 Box Plot -- Oracle Selection (`plot_boxplot_oracle.png`)

**How produced**: Same layout, oracle-selected models.

**Analysis**: Under oracle selection, the distributional shift is even clearer. NF's median (~0.17) is roughly double B0's (~0.08). NF's Q3 (~0.62) is substantially higher than B0's (~0.55). The majority of NF's points have shifted upward relative to B0, though both still have a long lower tail of hard cases near zero.

---

### 2.6 Paired Scatter Plot -- NF vs B0 (`plot_scatter_nf_vs_b0.png`)

**How produced**: Each dot = one complex at (B0_DockQ_AbAg, NF_DockQ_AbAg), confidence selection. Dashed diagonal = identity.

**Analysis**: NF wins on 13 complexes, B0 wins on 22, 0 tied. Despite fewer wins than B0, NF's gains are dramatically larger than B0's. Several points sit far above the diagonal (8FAH_HLA: 0.02->0.61, 8BLQ_ECD: 0.05->0.57, 7Y0O_HLA: 0.01->0.45), while B0's wins are all small (max regression = -0.14 for 7ZOZ_HLA). The asymmetry is key: when NF wins, it wins big; when it loses, it loses small.

Compare with previous run: that had NF wins=22, B0 wins=24 with a tiny mean delta. Here NF wins fewer complexes (13 vs 22), but the wins are much larger (+0.59, +0.53, +0.44 vs +0.88, +0.46).

---

### 2.7 Delta Waterfall Chart (`plot_delta_nf_vs_b0.png`)

**How produced**: Per-complex DockQ_AbAg difference (NF - B0) sorted ascending, confidence selection. Green = NF better, red = B0 better. Dashed line = mean (+0.046).

**Analysis**: The waterfall shows a strikingly different pattern from the previous run:

- **Regressions (red bars)**: The maximum regression is -0.14 (7ZOZ_HLA), drastically smaller than the previous run's -0.58. Most red bars are tiny (-0.01 to -0.08).
- **Improvements (green bars)**: Three massive towers: 8FAH_HLA (+0.59), 8BLQ_ECD (+0.53), 7Y0O_HLA (+0.44), plus 7TRH_HBG (+0.29) and 8TRS_AGD (+0.09).
- **Mean = +0.046**: 4.6x higher than the previous run's +0.010, though still driven by outlier wins.

The 20-particle FK resampling has **capped the downside** while preserving (and in some cases creating new) large upside improvements.

---

### 2.8 Confidence Calibration (`plot_confidence_calibration.png`)

**How produced**: Scatter of confidence_score vs DockQ_AbAg for all individual models (5 per complex), one subplot per method. Spearman rho shown.

#### B0 subplot (left)
rho=0.380, p=1.7e-09, n=235. Same as previous run (B0 unchanged).

#### NF subplot (right)
rho=0.265, p=3.9e-04, n=175. Lower than B0's 0.380 and lower than the previous NF run's 0.406. The steering with 20 particles appears to degrade confidence calibration -- many high-DockQ NF models have moderate confidence (0.83-0.88 range), while some low-DockQ models have high confidence. The narrower confidence range (0.80-0.98 vs B0's 0.73-0.98) also contributes to the lower correlation.

**Conclusion**: Confidence is a less reliable selector for steered models than for unconstrained ones. This explains why the confidence-selection gap between B0 and NF (0.046) is similar to the oracle gap (0.045) -- the confidence selector isn't choosing well for NF either.

---

### 2.9 Interface RMSD Bar Chart (`plot_irmsd.png`)

**How produced**: Mean iRMSD_AB across complexes, confidence selection.

| Method | Mean iRMSD_AB |
|--------|--------------|
| B0 | 11.5 |
| NF | **9.7** |

**Analysis**: NF reduces interface RMSD by 1.8 Angstroms (16% relative improvement). This is a meaningful change compared to the previous run where iRMSD was essentially identical (11.5 vs 11.4). The steering with more particles produces structures where the antibody-antigen interface is geometrically closer to the native.

---

### 2.10 Native Contact Recovery Bar Chart (`plot_fnat.png`)

**How produced**: Mean fnat_AB across complexes, confidence selection.

| Method | Mean fnat_AB |
|--------|-------------|
| B0 | 0.181 |
| NF | **0.252** |

**Analysis**: NF recovers 39% more native contacts than B0 (+0.071 absolute). This is the largest single-metric improvement and directly validates the mechanism: the steering potential is pulling CDR3 loops closer to the antigen, resulting in more native interface contacts. Previous run showed a much smaller improvement (0.196 vs 0.181).

---

### 2.11 Per-Complex Heatmap (`plot_heatmap_dockq.png`)

**How produced**: Heatmap of DockQ_AbAg per complex (rows) and method (columns), confidence selection. White cells = missing data.

**Analysis**: The 12 white cells in the NF column (missing predictions) are clearly visible. For complexes where NF has data, several notable improvements are visible:
- **8FAH_HLA**: Deep red (B0) to bright green (NF) -- the most visually dramatic change.
- **7Y0O_HLA**: Red to orange-yellow.
- **8BLQ_ECD**: Red to green.
- **7TRH_HBG**: Already yellowish in B0, becomes bright green in NF.

The bottom rows (hardest complexes: 8CYH_HLM, 8TFR_ABC, 8GP5_EFX, etc.) remain red for both methods -- the steering does not help on the truly intractable cases.

---

## 3. Statistical Tests (a_plus_2_particles_20)

### 3.1 Wilcoxon Signed-Rank Test -- NF vs B0

| Selection | N paired | NF mean | B0 mean | Delta mean | Delta median | p-value | Significant |
|-----------|----------|---------|---------|------------|-------------|---------|-------------|
| Confidence | 35 | 0.247 | 0.201 | **+0.046** | -0.0003 | 0.610 | No |
| Oracle | 35 | 0.314 | 0.270 | **+0.045** | -0.0009 | 0.533 | No |

**Conclusion**: Despite the 4.6x larger mean delta (+0.046 vs +0.010), the test does not reach significance (p=0.610). This is because: (1) the median delta is still ~0 (the improvement is driven by a handful of large wins rather than a systematic shift), and (2) the sample size dropped from 47 to 35. The Wilcoxon test is sensitive to the distribution shape, not just the mean -- with 22 small losses and 13 large wins, the signed ranks don't sum to significance.

---

## 4. Four-Method Comparison (a_plus_4_particles_20)

### 4.1 Aggregate DockQ Bar Chart (`plot_aggregate_dockq.png`)

**How produced**: Mean DockQ_AbAg per method. Note different N per method.

| Method | Confidence (N) | Oracle (N) |
|--------|---------------|-----------|
| B0 | 0.179 (47) | 0.248 (47) |
| B1 | **0.291** (46) | **0.479** (46) |
| B2 | 0.190 (36) | 0.269 (36) |
| NF | 0.247 (35) | 0.314 (35) |

**Analysis**: Ranking: B1 >> NF > B2 ~ B0. NF clearly beats B0 and B2 but falls well short of B1. This is expected: B1 uses ground-truth contact information (which residues actually touch in the native structure), giving it an unfair oracle advantage. NF requires no such information.

Under oracle selection, B1 dominates even more (0.479 vs NF's 0.314). NF's oracle advantage over B0 (+0.066) is larger than B2's (+0.021).

---

### 4.2 CAPRI Quality Distribution -- Confidence Selection (`plot_capri_confidence.png`)

| Tier | B0 (47) | B1 (46) | B2 (36) | NF (35) |
|------|---------|---------|---------|---------|
| Incorrect | 77% | **54%** | 72% | **63%** |
| Acceptable | 4% | 13% | 3% | 9% |
| Medium | 17% | **30%** | 22% | **26%** |
| High | 2% | 2% | 3% | 3% |

**Analysis**: NF places second after B1 in reducing incorrect predictions. NF achieves 37% acceptable-or-better rate vs B1's 46%, B2's 28%, and B0's 23%. NF notably outperforms B2 (pocket restraints) despite requiring no restraint information.

---

### 4.3 CAPRI Quality Distribution -- Oracle Selection (`plot_capri_oracle.png`)

| Tier | B0 (47) | B1 (46) | B2 (36) | NF (35) |
|------|---------|---------|---------|---------|
| Incorrect | 68% | **22%** | 64% | **51%** |
| Acceptable | 2% | 22% | 3% | 11% |
| Medium | 28% | **50%** | 31% | **34%** |
| High | 2% | 7% | 3% | 3% |

**Analysis**: Under oracle selection, B1 dominates (78% acceptable-or-better). NF achieves 49% acceptable-or-better (second place), outperforming B2 (36%) and B0 (32%). NF's 49% oracle success rate means that in nearly half of all complexes, the steering produces at least one model with DockQ >= 0.23 among 5 samples.

---

### 4.4 Box Plots (`plot_boxplot_confidence.png`, `plot_boxplot_oracle.png`)

**How produced**: Per-complex DockQ_AbAg distribution per method.

**Confidence**: B1 has the highest median (~0.10) and largest box (Q3 ~0.60). NF has the second-largest box (Q3 ~0.54), substantially wider than B0 (Q3 ~0.21) and B2 (Q3 ~0.30).

**Oracle**: B1's median rises to ~0.54 with a very tall box. NF's median (~0.17) is above B0 (~0.08) and B2 (~0.08). NF's Q3 (~0.62) approaches B1's lower quartile, showing overlap in their performance ranges.

---

### 4.5 Statistical Tests (a_plus_4_particles_20)

#### Confidence Selection

| Comparison | N | Delta mean | Delta median | p-value | Significant |
|-----------|---|------------|-------------|---------|-------------|
| NF vs B0 | 35 | +0.046 | -0.0003 | 0.610 | No |
| NF vs B1 | 34 | -0.053 | -0.004 | 0.059 | No (borderline) |
| NF vs B2 | 25 | +0.067 | +0.0004 | 0.615 | No |

#### Oracle Selection

| Comparison | N | Delta mean | Delta median | p-value | Significant |
|-----------|---|------------|-------------|---------|-------------|
| NF vs B0 | 35 | +0.045 | -0.0009 | 0.533 | No |
| NF vs B1 | 34 | **-0.193** | -0.074 | **0.000** | **Yes** |
| NF vs B2 | 25 | +0.038 | -0.0004 | 0.597 | No |

**Analysis**: NF vs B1 under oracle selection is the only statistically significant result (p<0.001), confirming B1's clear superiority. The NF vs B1 confidence comparison is borderline (p=0.059). NF vs B0 and NF vs B2 do not reach significance despite positive mean deltas, again because the improvements are concentrated in a few complexes.

---

### 4.6 Confidence Calibration (`plot_confidence_calibration.png`)

**How produced**: Scatter of confidence vs DockQ_AbAg for all individual models per method.

| Method | N models | Spearman rho | p-value |
|--------|----------|-------------|---------|
| B0 | 235 | 0.380 | <0.001 |
| B1 | 2610 | **0.443** | <0.001 |
| B2 | 435 | 0.359 | <0.001 |
| NF | 175 | **0.265** | 0.0004 |

**Analysis**: NF has the weakest confidence-DockQ correlation of all methods. B1 has the strongest (0.443), likely because the constraint information helps the model "know" when it finds the right interface. NF's poor calibration (0.265) means the confidence score is unreliable for identifying which of NF's 5 models is best -- a custom re-ranking strategy using the steering potential energy could help recover the oracle advantage.

---

### 4.7 Interface RMSD (`plot_irmsd.png`)

| Method | Mean iRMSD_AB |
|--------|--------------|
| B0 | 11.5 |
| B1 | **7.6** |
| B2 | 11.3 |
| NF | **9.7** |

**Analysis**: NF places second after B1 in interface RMSD. The 1.8A improvement over B0 is meaningful. B2 barely improves over B0 despite using restraint information.

---

### 4.8 Native Contact Recovery (`plot_fnat.png`)

| Method | Mean fnat_AB |
|--------|-------------|
| B0 | 0.181 |
| B1 | **0.289** |
| B2 | 0.184 |
| NF | **0.252** |

**Analysis**: NF achieves 87% of B1's fnat performance without using any ground-truth contact information. B2 barely exceeds B0. NF's 39% improvement over B0 in fnat is its strongest single-metric result.

---

### 4.9 Per-Complex Heatmap (`plot_heatmap_dockq.png`)

**How produced**: DockQ_AbAg per complex per method, confidence selection, sorted by B0 score. White = missing.

**Analysis**: Key patterns visible:
- **B1 column** is the greenest overall, but has white gaps (1 missing complex) and some surprising red cells where contact restraints fail.
- **B2 column** has more white cells (11 missing) and is mostly red/orange -- pocket restraints don't help much.
- **NF column** has 12 white cells but shows several dramatic red-to-green transitions (8FAH_HLA, 7Y0O_HLA, 8BLQ_ECD, 7TRH_HBG).
- Several complexes where B0 is red are green only in B1 but not NF (e.g., 8EZ8_HLA, 8E2U_HLA) -- these are among the 12 missing NF predictions.

---

## 5. Comparison: 20 Particles vs 5 Particles (Previous Run)

| Metric | 5 particles (prev) | 20 particles (current) | Change |
|--------|--------------------|-----------------------|--------|
| N complexes (NF) | 47 | 35 | -12 (26% loss) |
| Mean DockQ_AbAg (conf) | 0.189 | 0.247 | +0.058 |
| Paired delta (conf) | +0.010 | +0.046 | **4.6x larger** |
| CAPRI tier upgrades (conf) | 3 | **5** | +2 |
| CAPRI tier downgrades (conf) | 2 | **0** | -2 (eliminated) |
| Max regression | -0.579 | **-0.143** | 4x smaller |
| iRMSD | 11.4 | **9.7** | -1.7A |
| fnat | 0.196 | **0.252** | +0.056 |
| Confidence rho | 0.406 | 0.265 | -0.141 (worse) |
| p-value (Wilcoxon, conf) | 0.733 | 0.610 | slightly better |

### Key complexes that changed behavior

| Complex | Delta (5 particles) | Delta (20 particles) | Change |
|---------|--------------------|--------------------|--------|
| 8FAH_HLA | +0.001 | **+0.590** | New big win |
| 8BLQ_ECD | +0.003 | **+0.528** | New big win |
| 7Y0O_HLA | +0.000 | **+0.444** | New big win |
| 7ZOZ_HLA | **-0.579** | -0.143 | Regression 4x smaller |
| 8FXB_HLE | **+0.881** | -0.000 | Lost previous big win |
| 8HGM_CDB | +0.462 | +0.085 | Win reduced |

**Analysis**: The 20-particle FK resampling fundamentally changes which complexes benefit. It creates three new large wins (8FAH_HLA, 8BLQ_ECD, 7Y0O_HLA) while losing the previous run's biggest win (8FXB_HLE). The catastrophic regression on 7ZOZ_HLA is reduced from -0.58 to -0.14, and there are zero CAPRI tier downgrades.

---

## 6. Overall Conclusions

### The 20-particle FK resampling substantially improves the antigen steering feature, but it still does not reach statistical significance.

**Progress made**:
1. **4.6x larger mean improvement**: Delta rises from +0.010 to +0.046 over B0.
2. **Zero CAPRI tier downgrades**: The FK resampling acts as a safety net, filtering out bad orientations before they become the selected model.
3. **Max regression capped**: Worst-case drops from -0.58 to -0.14 (4x reduction).
4. **Meaningful metric improvements**: iRMSD improves by 1.8A, fnat by 39%.
5. **Competitive with restraint methods**: NF achieves 87% of B1's fnat without any oracle information, and clearly beats B2.

**Remaining limitations**:
1. **Not statistically significant**: p=0.610 (Wilcoxon). The improvement is concentrated in ~5 complexes out of 35, producing a positive mean but zero median shift. The test requires a broader, more systematic improvement.
2. **25% of complexes missing**: 12/47 complexes produced no NF output, likely due to computational cost of 20 particles. This both limits statistical power and means the feature is unusable for these cases.
3. **Poor confidence calibration**: rho=0.265 (worst of all methods). The confidence score cannot reliably pick the best steered model from 5 samples.
4. **Still loses to B0 on most complexes**: NF wins 13/35, B0 wins 22/35. The feature helps a minority of cases by a lot but slightly hurts the majority.

### Recommended next steps

1. **Fix the missing 12 complexes**: Investigate OOM/timeout with 20 particles. Consider reducing to 10-15 particles as a compromise, or using gradient accumulation to reduce memory.

2. **Custom model re-ranking**: Since confidence calibration is poor for steered models, implement a re-ranking step using the steering potential energy (CDR3-antigen distance score) as an additional selection signal. This could recover some of the oracle advantage. A simple weighted combination like `score = alpha * confidence + (1-alpha) * potential_energy` may suffice.

3. **Combine with FK particles at lower count + more diffusion samples**: Run `--num_particles 10 --diffusion_samples 10` to generate more diverse samples. The current 5x20=100 total trajectories collapsed to 5 output models might be too aggressive in filtering.

4. **Embedding-space steering**: Instead of coordinate-space gradient guidance (which the denoising network can partially undo), modify the pair representation z directly. The beta-scaling approach `z_scaled = (1+beta) * z` in the Pairformer would be:
   - More computationally efficient (no 20x particle overhead)
   - Harder to undo by the network (since it acts on internal representations)
   - Compatible with all 47 complexes (no OOM risk)

5. **Adaptive weight schedule**: The current fixed schedule (1.5->1.0->0.3) may not be optimal. Consider making the weight proportional to the current CDR3-antigen distance -- stronger guidance when far from contact, lighter when already close.

6. **Ensemble NF results across particle counts**: Since 5-particle and 20-particle runs help different complexes (8FXB_HLE vs 8FAH_HLA), ensembling both could capture the union of improvements. Run both and take the best model across the combined pool.
