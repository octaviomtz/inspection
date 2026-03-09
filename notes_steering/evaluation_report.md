# Evaluation Report: Antigen Steering Feature (Strategy A)

**Date**: 2026-03-09
**Evaluation tool**: `scripts/evaluate.py` (DockQ v2 via `conda run -n dockq2`)
**Results directory**: `a_plus_2/` (NF vs B0 baseline)
**Note on previous evaluation**: An earlier evaluation (`a_plus_4/`) incorrectly showed NF = B0 due to a DockQ cache key collision bug in `evaluate.py`. That bug has been fixed and the results below are from the corrected evaluation.

---

## 1. Evaluation Setup

### Methods Compared

| ID | Method | Description | n_models per complex |
|----|--------|-------------|---------------------|
| B0 | Unconstrained | Boltz-2 with no steering (`boltz predict`) | 5 |
| NF | Antigen steering | Boltz-2 with `--use_potentials --antigen_steering` | 5 |

### How the evaluation works

The script runs DockQ v2 on every predicted model (5 per complex per method) against the ground-truth PDB from `pdb_minimized/`. DockQ decomposes each 3-chain antibody-antigen complex into 3 interfaces: AB (antigen-heavy), AC (antigen-light), BC (VH-VL). Two model selection strategies are applied:

- **Confidence selection**: Pick the model with the highest Boltz confidence_score (= 0.8 * iptm + 0.2 * ptm). This is what a user would do without access to ground truth.
- **Oracle selection**: Pick the model with the highest DockQ_AbAg. This is an upper bound showing the best the method can produce among its 5 samples.

### Primary metric

**DockQ_AbAg** = mean(DockQ_AB, DockQ_AC), the average DockQ score across both antibody-antigen interfaces.

### Dataset

47 antibody-antigen complexes. Both B0 and NF have predictions for all 47.

---

## 2. Plot-by-Plot Analysis

### 2.1 Aggregate DockQ Bar Chart (`plot_aggregate_dockq.png`)

**How produced**: For each method, compute the mean DockQ_AbAg across all 47 complexes. Solid bars show confidence-selected models; hatched bars show oracle-selected models.

| Method | Confidence | Oracle |
|--------|-----------|--------|
| B0 | 0.179 | 0.248 |
| NF | 0.189 | 0.258 |
| Delta (NF - B0) | +0.010 | +0.010 |

**Analysis**: NF shows a small +0.010 improvement in mean DockQ_AbAg over B0 under both selection strategies. The improvement is consistent but modest. Oracle scores are ~0.07 higher than confidence scores for both methods, indicating that among 5 samples the best model is often not the most confident one.

---

### 2.2 CAPRI Quality Distribution -- Confidence Selection (`plot_capri_confidence.png`)

**How produced**: Each complex's confidence-selected DockQ_AbAg is classified into CAPRI tiers: Incorrect (<0.23), Acceptable (0.23-0.49), Medium (0.49-0.80), High (>=0.80). Shown as stacked percentage bars.

| Tier | B0 | NF | Change |
|------|----|----|--------|
| Incorrect | 36/47 (77%) | 35/47 (74%) | -3pp |
| Acceptable | 2/47 (4%) | 2/47 (4%) | 0 |
| Medium | 8/47 (17%) | 8/47 (17%) | 0 |
| High | 1/47 (2%) | 2/47 (4%) | +2pp |

**Analysis**: NF reduces the incorrect rate by 3 percentage points and gains one additional High-quality prediction (8FXB_HLE, from 0.009 to 0.890). The Acceptable and Medium tiers remain unchanged. The distributions are very similar overall.

**CAPRI tier changes (3 upgrades, 2 downgrades)**:
- Upgraded: 8FXB_HLE (Incorrect->High, +0.88), 8HGM_CDB (Incorrect->Medium, +0.46), 7TRH_HBG (Acceptable->Medium, +0.006)
- Downgraded: 7ZOZ_HLA (Medium->Incorrect, -0.58), 8EAY_HLA (Medium->Acceptable, -0.04)

---

### 2.3 CAPRI Quality Distribution -- Oracle Selection (`plot_capri_oracle.png`)

**How produced**: Same as above but using the model with the highest DockQ_AbAg (best possible among 5 samples).

| Tier | B0 | NF | Change |
|------|----|----|--------|
| Incorrect | 32/47 (68%) | 29/47 (62%) | -6pp |
| Acceptable | 1/47 (2%) | 5/47 (11%) | +9pp |
| Medium | 13/47 (28%) | 11/47 (23%) | -5pp |
| High | 1/47 (2%) | 2/47 (4%) | +2pp |

**Analysis**: Under oracle selection NF shows its strongest signal. The incorrect rate drops 6 percentage points (68% to 62%), and Acceptable jumps from 2% to 11%. This means NF produces at least one acceptable-or-better model in 38% of complexes vs 32% for B0. The 5 newly upgraded complexes under oracle selection (7TRI_ZYB, 8EZ7_HLA, 8FXB_HLE, 8OL9_BAH, 8TRS_AGD) suggest the steering does explore alternative orientations that sometimes land closer to the native interface.

However, 3 complexes are downgraded (7Y0O_HLA, 8EZ8_HLA, 8HES_HLC: Medium->Incorrect or Acceptable), showing the steering can also destroy good models.

---

### 2.4 Box Plot -- Confidence Selection (`plot_boxplot_confidence.png`)

**How produced**: Box-and-whisker plot of per-complex DockQ_AbAg (47 values per method), with individual points jittered alongside. Median shown as the horizontal line in the box, box spans Q1-Q3.

**Analysis**: Both distributions are extremely right-skewed. The median is near zero for both methods (~0.04), reflecting that most complexes produce incorrect docking. NF's box extends slightly higher (larger Q3), and its outliers reach further (0.89 for 8FXB_HLE). The B0 distribution has a slightly higher maximum non-outlier whisker. Overall the two distributions are visually indistinguishable -- the signal is dominated by a few outlier complexes rather than a systematic shift.

---

### 2.5 Box Plot -- Oracle Selection (`plot_boxplot_oracle.png`)

**How produced**: Same as 2.4 but using oracle-selected models.

**Analysis**: Oracle selection lifts the medians slightly (B0: 0.079, NF: 0.059). Both distributions broaden, showing that each method's 5-sample pool contains some better models. NF's Q3 is higher but its median is slightly lower. The overall shape remains very similar between methods. The takeaway is that neither method consistently produces good docks for the majority of complexes.

---

### 2.6 Paired Scatter Plot -- NF vs B0 (`plot_scatter_nf_vs_b0.png`)

**How produced**: Each dot is one complex, plotted at (B0_DockQ_AbAg, NF_DockQ_AbAg) using confidence-selected models. The dashed diagonal is the identity line (NF = B0); points above it mean NF is better.

**Analysis**: NF wins on 22 complexes, B0 wins on 25, 0 tied. The split is nearly even (47% vs 53%), indicating no systematic advantage in either direction. Most points cluster near (0, 0) -- both methods fail on the same hard complexes. The scatter around the diagonal for complexes with moderate-to-high DockQ (>0.3) is substantial, showing the methods produce meaningfully different structures. Key outliers:

- **Far above diagonal**: 8FXB_HLE (B0=0.009, NF=0.890) -- NF finds the correct orientation where B0 completely fails. 8HGM_CDB (B0=0.228, NF=0.691) -- major improvement.
- **Far below diagonal**: 7ZOZ_HLA (B0=0.644, NF=0.065) -- NF destroys a good B0 prediction. 8OL9_BAH (B0=0.193, NF=0.080) -- moderate regression.

The mean improvement (+0.010) is driven almost entirely by the 8FXB_HLE and 8HGM_CDB outliers.

---

### 2.7 Delta Waterfall Chart (`plot_delta_nf_vs_b0.png`)

**How produced**: Per-complex DockQ_AbAg difference (NF - B0), sorted from most negative to most positive (confidence selection). Green bars = NF improvement, red bars = NF regression. Dashed horizontal line marks the mean delta (+0.010).

**Analysis**: The waterfall reveals a strongly asymmetric distribution:

- **Left tail (regressions)**: One large regression dominates (7ZOZ_HLA at -0.58), followed by several moderate regressions (-0.02 to -0.11). Most red bars are small.
- **Right tail (improvements)**: Two massive gains (8FXB_HLE at +0.88, 8HGM_CDB at +0.46) tower over all other bars. Without these two complexes, the mean delta would be approximately -0.02 (net negative).
- **Middle**: Most complexes (~35 of 47) have deltas between -0.05 and +0.05, indicating negligible change.

The mean delta of +0.010 is not robust -- it depends on two outlier complexes. The median delta is -0.0001 (essentially zero).

---

### 2.8 Confidence Calibration (`plot_confidence_calibration.png`)

**How produced**: Scatter of confidence_score vs DockQ_AbAg for all 235 individual models (47 complexes x 5 models), one subplot per method. Spearman rho and p-value shown in the title.

#### B0 subplot (left)
rho=0.380, p=1.7e-09, n=235. Confidence scores range 0.73-0.98. A moderate positive correlation exists: the highest-confidence models (>0.92) tend to be the highest-DockQ ones, but many high-confidence models (0.85-0.92) have near-zero DockQ. The large cluster at DockQ~0 across all confidence levels shows that confidence alone cannot reliably distinguish good from bad docks.

#### NF subplot (right)
rho=0.406, p=9.3e-11, n=235. Slightly better calibration than B0. The confidence range is wider (0.70-0.98), suggesting the steering introduces more variation in model confidence. The correlation structure is similar to B0 but with a few more high-DockQ points at moderate confidence levels (0.85-0.90), suggesting some steered models achieve good docking despite not being the highest-confidence.

**Conclusion**: Both methods show moderate confidence-DockQ correlation. NF's slightly higher rho (0.406 vs 0.380) indicates the steering does not degrade confidence calibration.

---

### 2.9 Interface RMSD Bar Chart (`plot_irmsd.png`)

**How produced**: Mean iRMSD_AB (antigen-heavy chain interface RMSD) across all 47 complexes, using confidence-selected models. Lower is better.

| Method | Mean iRMSD_AB |
|--------|--------------|
| B0 | 11.5 |
| NF | 11.4 |

**Analysis**: Virtually identical. Both values are very high (>10A), reflecting that most confidence-selected models have poor interface geometry. The 0.1A difference is negligible.

---

### 2.10 Native Contact Recovery Bar Chart (`plot_fnat.png`)

**How produced**: Mean fnat_AB (fraction of native contacts recovered at the antigen-heavy chain interface) across all 47 complexes, using confidence-selected models. Higher is better.

| Method | Mean fnat_AB |
|--------|-------------|
| B0 | 0.181 |
| NF | 0.196 |

**Analysis**: NF shows a small improvement in native contact recovery (+0.015). This is consistent with the DockQ improvement and suggests the steering does bring more antigen residues into proximity with CDR loops, recovering slightly more native interface contacts. However, both values are low (<0.2), meaning the average prediction recovers fewer than 20% of native contacts.

---

### 2.11 Per-Complex Heatmap (`plot_heatmap_dockq.png`)

**How produced**: Heatmap of DockQ_AbAg per complex (rows) per method (columns), confidence selection. Complexes are sorted by B0 score (highest at top). Color scale: dark red = 0 (incorrect), dark green = 1 (perfect).

**Analysis**: The heatmap visually confirms that most complexes are deep red (low DockQ) for both methods. Key observations:

- **Top rows** (high B0): Most retain similar quality under NF (similar green shading). Exception: 7ZOZ_HLA drops from green (0.64) to red (0.07).
- **8FXB_HLE**: Stands out dramatically -- pure red in B0, bright green in NF (0.89). This is the single most striking result.
- **8HGM_CDB**: Also shows a clear color shift from orange/red (B0=0.23) to green (NF=0.69).
- **Bottom half**: Uniformly red for both methods -- the hard complexes remain unsolved.

---

## 3. Statistical Tests

### 3.1 Wilcoxon Signed-Rank Test -- Confidence Selection

| Comparison | N | NF mean | B0 mean | Delta mean | Delta median | p-value | Significant |
|-----------|---|---------|---------|------------|-------------|---------|-------------|
| NF vs B0 | 47 | 0.189 | 0.179 | +0.010 | -0.000 | 0.733 | No |

### 3.2 Wilcoxon Signed-Rank Test -- Oracle Selection

| Comparison | N | NF mean | B0 mean | Delta mean | Delta median | p-value | Significant |
|-----------|---|---------|---------|------------|-------------|---------|-------------|
| NF vs B0 | 47 | 0.258 | 0.248 | +0.010 | -0.001 | 0.375 | No |

**Conclusion**: Neither test reaches significance (p=0.733 and p=0.375). The mean improvement of +0.010 is not statistically distinguishable from zero. The median delta is essentially zero in both cases, confirming that the improvement is driven by a few outlier complexes rather than a systematic shift.

---

## 4. VH-VL Interface Stability

| Metric | B0 | NF |
|--------|----|----|
| Mean DockQ_BC (confidence) | 0.655 | 0.653 |

The antibody fold quality is virtually identical between methods. The steering potential, which acts only on antigen-CDR3 distances, does not damage the VH-VL interface.

---

## 5. Overall Conclusions

### The antigen steering produces different structures but does not achieve statistically significant improvement.

Unlike the previous (buggy) evaluation which showed B0 = NF due to a cache collision, the corrected evaluation confirms the steering **does** alter the predicted structures. However:

1. **Marginal mean improvement**: NF improves mean DockQ_AbAg by +0.010 over B0 (both confidence and oracle), but this is not statistically significant (p=0.733 confidence, p=0.375 oracle).

2. **Nearly even win rate**: NF wins on 22/47 complexes, B0 wins on 24/47, 1 tied (confidence selection). The steering helps about as often as it hurts.

3. **Driven by outliers**: The +0.010 mean delta is almost entirely due to two complexes: 8FXB_HLE (+0.88) and 8HGM_CDB (+0.46). Without these, the mean delta would be negative (-0.02). This suggests the steering occasionally finds dramatically better orientations but more often introduces noise.

4. **High variance**: The steering increases outcome variance -- it can produce spectacular improvements (8FXB_HLE: 0.009 -> 0.890) but also devastating regressions (7ZOZ_HLA: 0.644 -> 0.065).

5. **Oracle CAPRI shows promise**: Under oracle selection, NF reduces incorrect predictions from 68% to 62% and increases Acceptable from 2% to 11%. This means the steering does produce some good alternative orientations, but they are not reliably selected by confidence.

6. **No collateral damage**: VH-VL interface quality (DockQ_BC ~0.65) and confidence calibration (rho ~0.40) are preserved.

### What the results suggest about the mechanism

The steering potential is operating -- structures differ between B0 and NF. But the effect is weak and inconsistent:

- When it works (8FXB_HLE, 8HGM_CDB), it appears to have successfully reoriented the antigen toward the CDR3 loops, finding a qualitatively different and correct binding mode.
- When it fails (7ZOZ_HLA, 8OL9_BAH), it may be pulling the antigen away from a good orientation that the denoising network had found on its own.
- For most complexes (~35/47), the perturbation is too small to change the outcome.

### Recommended modifications

1. **Increase guidance strength**: The current weight schedule (1.5 -> 1.0 -> 0.3) and 20 GD steps may be too weak. Increasing by 5-10x could push the steering from "noise" to "signal" territory. The risk is amplifying regressions too, so this should be paired with FK resampling.

2. **Enable FK resampling**: Running with `--num_particles >= 3` enables the Feynman-Kac particle mechanism, which scores and filters orientations rather than pushing coordinates directly. This is fundamentally more robust than gradient guidance alone because it selects good orientations rather than forcing them.

3. **Move to embedding space**: Optimizing the pair representation z rather than atom coordinates would prevent the denoising network from "erasing" the guidance signal at each step.

4. **Beta-scaling approach**: Applying z_scaled = (1 + beta) * z to pair representations in the Pairformer is computationally cheaper and more direct than coordinate-space guidance.

5. **Reduce potential scope**: Focus the potential on the top-k closest antigen residues (k=10-20) rather than all antigen CA atoms, to concentrate the gradient signal.

6. **Improve model selection**: The oracle CAPRI results show NF produces more acceptable models than B0, but confidence doesn't find them. A re-ranking strategy using the steering potential's energy as an additional selection signal could recover some of the oracle advantage.
