# Analysis Report: L+ (CDR3 Beta-Scaling) Evaluation Results

## Overview

**Feature**: CDR3 Beta-Scaling (Idea L+) — scales attention biases for CDR3 pair interactions in the token transformer by a factor `(1 + beta)` after the LayerNorm+Linear projection. Applied with `beta = 0.3` to all CDR3-H and CDR3-L regions.

**Dataset**: 47 antibody-antigen complexes, 5 diffusion samples per complex, seed 42.

**Evaluations**:
- **l_plus_2**: l_cdr3_beta_scaling vs. 1 baseline (baseline_antigen_cut, n=47)
- **l_plus_4**: l_cdr3_beta_scaling vs. 3 baselines (baseline_antigen_cut n=47, baseline_contact_restraints n=46, baseline_pocket n=36)

---

## Plot-by-Plot Analysis

### 1. DockQ Comparison Boxplot (`dockq_comparison_boxplot.png`)

**How produced**: Boxplot over **all 5 models per complex** (not best-model selection). Three subplots show DockQ(Ab-Ag avg), DockQ(A-B), and DockQ(A-C). Horizontal reference lines mark CAPRI quality thresholds: Acceptable (0.23), Medium (0.49), High (0.80).

#### Subplot 1 — DockQ (Ab-Ag avg):
- **l_cdr3_beta_scaling** (red) shows a slightly higher median and a taller upper quartile compared to **baseline_antigen_cut** (blue), but the difference is small. Both methods have medians near 0.04, well below the Acceptable threshold.
- **baseline_contact_restraints** (orange) is clearly the best performer, with median ~0.10 and a much larger interquartile range reaching into the Acceptable/Medium range.
- **baseline_pocket** (green) is comparable to antigen_cut.

#### Subplot 2 — DockQ (A-B) (antigen vs. heavy chain):
- Same pattern: contact_restraints dominates. l_cdr3_beta_scaling has a slightly higher 75th percentile than antigen_cut but nearly identical median.

#### Subplot 3 — DockQ (A-C) (antigen vs. light chain):
- Similar to A-B. l_cdr3_beta_scaling shows marginal improvement over antigen_cut but remains far below contact_restraints.

**Conclusion**: CDR3 beta-scaling provides a marginal improvement in DockQ distribution over the simplest baseline (antigen_cut) but does not approach the performance of contact_restraints. The improvement is not statistically significant (Wilcoxon p=0.67 for dockq_antigen_ab).

---

### 2. DockQ Scatter Per Baseline (`dockq_scatter_per_baseline.png`)

**How produced**: Each point is one complex using the **best model** (highest DockQ) per complex per method. Scatter compares l_cdr3_beta_scaling (y-axis) vs. baseline_antigen_cut (x-axis). Dashed line = y=x (parity). Three subplots.

#### Subplot 1 — DockQ (Ab-Ag):
- Points scatter roughly evenly around the diagonal. Several complexes show large improvements (points well above the line, e.g., one case jumping from ~0.2 to ~0.75), but others show regressions (points below the line, e.g., one dropping from ~0.35 to ~0.33).
- Win/loss record: 23 wins, 24 losses — essentially a coin flip.
- A few outliers with large positive deltas (e.g., 8HGM_CDB, 8SLB_HLA) are counterbalanced by comparable negative deltas (e.g., 8EZ8_HLA).

#### Subplot 2 — CDR3 RMSD:
- Most points cluster near the diagonal in the 0-8 A range. A few high-RMSD outliers remain high for both methods.
- Slight tendency for points to fall above the diagonal (beta-scaling produces slightly worse RMSD), consistent with the statistical test (10 wins, 37 losses, p=0.0008).

#### Subplot 3 — iPTM:
- Striking pattern: **nearly all points lie above the diagonal**, indicating beta-scaling systematically increases iPTM scores across almost every complex.
- This is the strongest and most statistically significant effect (40 wins, 7 losses, p<0.0001).
- However, iPTM is a model-confidence metric, not a ground-truth quality metric. Higher iPTM without corresponding DockQ improvement suggests the model is becoming **overconfident** rather than more accurate.

**Conclusion**: DockQ shows no systematic improvement (evenly split wins/losses). CDR3 RMSD is slightly degraded. iPTM is strongly inflated, suggesting the beta-scaling makes the model more confident in its CDR3 predictions without improving their accuracy.

---

### 3. CDR3 RMSD Comparison (`cdr3_rmsd_comparison.png`)

**How produced**: Boxplot over **all models** (all 5 samples per complex). Three subplots: CDR3-H RMSD, CDR3-L RMSD, CDR3 Combined RMSD. RMSD is computed by superimposing framework residues and measuring CDR3 backbone Ca displacement.

#### Subplot 1 — CDR3-H RMSD:
- All four methods show similar distributions with median ~3 A and upper whiskers reaching ~8 A. l_cdr3_beta_scaling (red) has a slightly higher median (~3.1 A) than antigen_cut (~2.9 A) and contact_restraints (~2.8 A).
- Some extreme outliers (>25 A) appear for all methods.

#### Subplot 2 — CDR3-L RMSD:
- Distributions are tighter (median ~1.3-1.7 A). l_cdr3_beta_scaling shows a slightly higher median (~1.7 A) compared to all baselines (~1.3-1.5 A).

#### Subplot 3 — CDR3 Combined RMSD:
- Same pattern: l_cdr3_beta_scaling (median ~3.0 A) is marginally worse than all baselines.

**Conclusion**: CDR3 beta-scaling produces **statistically significantly worse** CDR3 RMSD compared to all three baselines (p=0.0008 vs antigen_cut, p=0.0002 vs contact_restraints, p=0.004 vs pocket). The degradation is modest in absolute terms (~0.1-0.4 A increase in mean) but highly consistent across complexes (37/47 losses vs antigen_cut).

---

### 4. CDR3 pLDDT Comparison (`cdr3_plddt_comparison.png`)

**How produced**: Boxplot over **all models**. Three subplots: CDR3-H pLDDT, CDR3-L pLDDT, CDR3 Combined pLDDT. pLDDT is the predicted local distance difference test (model confidence in per-residue accuracy).

#### Subplot 1 — CDR3-H pLDDT:
- l_cdr3_beta_scaling (red) has the lowest median (~0.84) and the widest spread. contact_restraints (orange) has the highest median (~0.86). antigen_cut and pocket are intermediate.

#### Subplot 2 — CDR3-L pLDDT:
- l_cdr3_beta_scaling shows the lowest median (~0.85) while contact_restraints leads (~0.90). The light-chain CDR3 generally has higher confidence across all methods.

#### Subplot 3 — CDR3 Combined pLDDT:
- l_cdr3_beta_scaling has the lowest combined pLDDT (mean 0.820 vs baseline 0.838, contact_restraints 0.865).

**Conclusion**: CDR3 beta-scaling **reduces local confidence** (pLDDT) at CDR3 positions. This is an interesting contrast with the iPTM increase: the model is globally more confident (iPTM) but locally less confident at the CDR3 loops themselves. This suggests beta-scaling disrupts local CDR3 structure prediction while inflating global interface confidence.

---

### 5. Ensemble CDR3 Diversity (`ensemble_diversity.png`)

**How produced**: Bar chart showing **mean pairwise CDR3 Ca RMSD** across the 5 diffusion samples per complex, averaged over all complexes. Error bars show standard deviation. Red dashed line marks 1.0 A structural diversity threshold.

- l_cdr3_beta_scaling achieves the **highest ensemble diversity**: 1.66 A mean (vs 1.51 A antigen_cut, 1.41 A contact_restraints, 1.44 A pocket).
- Fraction of complexes with diversity >= 1.0 A: **0.60** for beta-scaling vs 0.45 for antigen_cut (a 33% relative improvement).
- Error bars are large for all methods, indicating high variance across complexes.

**Conclusion**: This is the **strongest positive result** for CDR3 beta-scaling. The feature achieves its stated design goal of increasing CDR3 conformational diversity. The 10% absolute increase in mean diversity (1.51 -> 1.66 A) and the large jump in fraction of diverse ensembles (45% -> 60%) confirm that modulating attention biases does produce more structurally varied CDR3 conformations.

---

### 6. Per-Complex Heatmap (`per_complex_heatmap.png`)

**How produced**: Shows DeltaDockQ = DockQ(l_cdr3_beta_scaling) - DockQ(baseline) for the **best model** per complex. Rows are complexes, columns are baselines. RdYlGn colormap: green = improvement, red = degradation, yellow = no change. Scale: -0.3 to +0.3.

#### Column 1 — vs. antigen_cut:
- Predominantly yellow-green (slight improvements or neutral). Notable strong improvements (dark green): 8HGM_CDB (+~0.3), 8SLB_HLA (+~0.15), 8EZ7_HLA (+~0.1).
- A few regressions (red): 8EZ8_HLA (~-0.2), 8IVA_CEG (~-0.05).
- Most complexes show very small deltas (pale yellow), consistent with non-significant Wilcoxon test.

#### Column 2 — vs. contact_restraints:
- More red than green. Several strong regressions: 7TRI_ZYB, 7Y0O_HLA, 8DTK_CBA, 8GP5_EFX, 8IV5_ABG (~-0.2 to -0.3).
- A few strong improvements: 8CDE_DCB (+0.3), 8FAH_HLA (~+0.1).
- Overall pattern is clearly negative, consistent with the significant Wilcoxon result (p=0.005).

#### Column 3 — vs. pocket:
- Mixed pattern but more green than red. Several improvements: 8DE3_BCA, 8SLB_HLA, 8EZ8_HLA (+0.1 to +0.2).
- A few regressions: 8GP5_EFX, 8DTK_CBA (~-0.15).
- Overall roughly balanced but slightly positive, consistent with non-significant Wilcoxon (p=0.33).

**Conclusion**: The per-complex heatmap reveals that CDR3 beta-scaling shows **case-specific behavior**: a few complexes benefit substantially while most show negligible change. The feature is clearly inferior to contact_restraints across the board. Against antigen_cut and pocket, the improvements are concentrated in specific complexes rather than being a general trend.

---

## Summary Statistics

| Metric | antigen_cut | contact_restraints | pocket | l_cdr3_beta_scaling |
|--------|-------------|-------------------|--------|---------------------|
| DockQ (Ab-Ag) mean | 0.179 | **0.291** | 0.190 | 0.201 |
| CDR3 Combined RMSD (A) | 3.64 | **3.37** | 3.73 | 3.76 |
| CDR3 Combined pLDDT | 0.838 | **0.865** | 0.848 | 0.820 |
| Ensemble diversity (A) | 1.51 | 1.41 | 1.44 | **1.66** |
| Frac diverse (>=1A) | 0.45 | 0.47 | 0.53 | **0.60** |
| iPTM | 0.737 | 0.834 | 0.794 | **0.828** |
| CAPRI Medium+High | 9 | **15** | 9 | 11 |

**Bold** = best value for that metric.

---

## CAPRI Classification

| Method | Incorrect | Acceptable | Medium | High |
|--------|-----------|------------|--------|------|
| baseline_antigen_cut | 36 | 2 | 8 | 1 |
| baseline_contact_restraints | **25** | **6** | **14** | 1 |
| baseline_pocket | 26 | 1 | 8 | 1 |
| l_cdr3_beta_scaling | 34 | 2 | 10 | 1 |

CDR3 beta-scaling moves 2 complexes from Incorrect to Medium compared to antigen_cut, but remains far behind contact_restraints which achieves 15 Medium+High predictions vs 11.

---

## Overall Conclusions

### Did the feature achieve its goal?

**Partially.** The CDR3 beta-scaling feature was designed to modulate CDR3 structural diversity by scaling pair representation attention biases. The evaluation reveals:

1. **Ensemble diversity: SUCCESS.** The feature achieves the highest CDR3 ensemble diversity (1.66 A vs 1.51 A baseline), with 60% of complexes producing structurally diverse ensembles (vs 45% baseline). This confirms the mechanism works as intended — amplifying CDR3 pair attention biases does produce more varied CDR3 conformations.

2. **Docking quality: MARGINAL.** DockQ improves slightly (0.179 -> 0.201, +12% relative) but the improvement is not statistically significant (p=0.67). The feature does push 2 additional complexes into CAPRI Medium quality. However, it remains substantially inferior to contact_restraints (0.291).

3. **CDR3 structural accuracy: DEGRADED.** CDR3 backbone RMSD worsens significantly (3.64 -> 3.76 A, p=0.0008) with the feature losing on 37/47 complexes. The increased diversity comes at the cost of accuracy — the model explores more conformations but does not converge on better ones.

4. **Model confidence: INFLATED.** iPTM increases dramatically (0.737 -> 0.828, p<0.0001) while pLDDT decreases (0.838 -> 0.820). This divergence between global confidence (iPTM up) and local accuracy (pLDDT down, RMSD up) indicates the scaling creates an artificial inflation of interface confidence without corresponding structural improvement.

### What the feature does well
- Increases CDR3 conformational exploration (higher diversity)
- Provides modest DockQ improvements on specific complexes (case-dependent)
- The mechanism (attention bias scaling) is confirmed to have a real effect on the model

### What the feature does poorly
- No systematic improvement in docking quality
- Degrades CDR3 loop accuracy across most complexes
- Creates a confidence-accuracy mismatch (overconfident iPTM)
- Significantly inferior to contact_restraints on all structural metrics

---

## Recommended Modifications

1. **Reduce beta value.** The current `beta=0.3` may be too aggressive. Try `beta=0.1` or `beta=0.15` to find a regime where diversity increases without degrading accuracy as much. Run a beta sweep: {-0.1, 0.0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5}.

2. **Asymmetric CDR3 scaling (Idea O+).** CDR3-H is structurally more variable and harder to predict than CDR3-L. Apply different beta values: stronger beta on CDR3-H for exploration, weaker or zero on CDR3-L to preserve its typically more accurate conformation.

3. **Expand the scaling mask beyond CDR3-CDR3 pairs.** Currently only CDR3_i x CDR3_j pairs are scaled. Consider also scaling CDR3 x antigen pairs (the interface pairs most relevant to docking). This could steer the model to explore different binding modes rather than just different loop conformations.

4. **Combine with contact restraints.** Contact_restraints is the strongest baseline. Layer CDR3 beta-scaling on top of contact_restraints to see if the diversity benefit adds value when the model already has better docking guidance.

5. **Time-dependent beta scheduling.** Apply high beta (exploration) in early diffusion steps and low beta (refinement) in later steps. This could capture the diversity benefit while allowing the model to converge on accurate structures.

6. **Best-of-N selection using confidence.** Since beta-scaling increases diversity, use it to generate larger ensembles (e.g., 10-20 samples) and then select the best model by confidence ranking. The increased diversity may yield better best-model picks even if the mean quality is slightly lower.

7. **Additive bias injection instead of multiplicative scaling.** The current approach scales existing biases: `bias + beta * mask * bias`. For positions where the existing bias is near zero, this has minimal effect. Consider adding a learned or fixed additive term: `bias + beta * mask * constant` to ensure CDR3 pairs always receive a meaningful perturbation regardless of the base bias magnitude.
