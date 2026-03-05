# Evaluation Report: E+ (Enhanced Blind Scanning)

## Overview

This report analyzes the evaluation results of the E+ (Enhanced Blind Scanning) feature against up to three baselines across 46-47 antibody-antigen complexes. Two evaluation runs were performed:

- **e_plus_2/**: E+ vs B1 (vanilla Boltz2)
- **e_plus_4/**: E+ vs B1, B2 (contact restraints), and B3 (pocket+antibody constraints)

Both were generated using `evaluate_e_plus.py`. All plots use the **best model per complex** (selected by highest `confidence_score` from the 5 Boltz2 models). Baselines B2 and B3 additionally use a **best-per-complex config** strategy (selecting the restraint configuration with the highest confidence score).

---

## Plot-by-Plot Analysis

### 1. AUC-ROC Bar Plot (`auc_roc_bar.png`)

**Methodology:** Mean AUC-ROC across all complexes for each method, with standard deviation error bars. AUC-ROC is computed from the continuous epitope propensity heatmap (E+ only) against the ground truth binary epitope labels (5.0 A heavy-atom distance cutoff). Baselines B1/B2/B3 produce binary epitope predictions from their predicted structures, so AUC-ROC is not applicable to them (bars are absent).

**Results:**
- E+ mean AUC-ROC = **0.645** (std = 0.216)
- B1/B2/B3: not applicable (no continuous predictions)

**Conclusion:** The mean AUC-ROC of 0.645 falls **below the target of 0.7** set in the evaluation strategy. The large standard deviation (0.216) indicates highly variable performance: 46% of complexes achieve AUC-ROC > 0.7, while 65% exceed 0.5 (random baseline). The heatmap provides **better-than-random discrimination** for the majority of complexes, but is not reliably accurate across the full test set.

---

### 2. AUC-PR Bar Plot (`auc_pr_bar.png`)

**Methodology:** Mean Area Under the Precision-Recall curve, computed the same way as AUC-ROC but using PR curves instead. AUC-PR is more informative for imbalanced data (epitopes typically comprise 15-30% of antigen residues).

**Results:**
- E+ mean AUC-PR = **0.304** (std = 0.238)

**Conclusion:** The AUC-PR of 0.304 falls **below the target of 0.5**. This is the more concerning metric: while AUC-ROC can appear reasonable due to the large number of true negatives, the PR curve reveals that when E+ predicts a residue as part of the epitope, it is correct only ~30% of the time on average. The large variance (std = 0.238) means some complexes get good epitope recovery while many do not.

---

### 3. Epitope F1 Bar Plot (`epitope_f1_bar.png`)

**Methodology:** Mean F1 score per method across all complexes. For E+, the F1 is computed at the **optimal threshold** (threshold that maximizes F1 per complex). For baselines, epitope residues are extracted from the predicted structure using the same 5.0 A distance cutoff as ground truth. Error bars show standard deviation.

**Results (e_plus_4, 4 methods):**
| Method | Mean F1 |
|--------|---------|
| B1     | 0.345   |
| B2     | **0.492** |
| B3     | 0.347   |
| E+     | 0.340   |

**Conclusion:** E+ achieves **essentially the same epitope F1 as B1** (0.340 vs 0.345) and B3 (0.347), while B2 (contact restraints, best config) leads with 0.492. This is a critical finding: **E+ does not improve epitope identification over vanilla Boltz2 (B1)** even when using an optimized threshold. B2 benefits from oracle contact information, making it a strong but unfair comparison (B2 requires prior knowledge of contact residues). All methods show high variance (std ~0.25-0.34), reflecting the difficulty of the task.

---

### 4. DockQ Bar Plot (`dockq_bar.png`)

**Methodology:** Mean DockQ score averaged across the antibody-antigen interfaces (AB and AC chains) per method. DockQ combines iRMSD, LRMSD, and fnat into a single score: < 0.23 = Incorrect, 0.23-0.49 = Acceptable, 0.49-0.80 = Medium, >= 0.80 = High.

**Results (e_plus_4):**
| Method | Mean DockQ (AB+AC) |
|--------|-------------------|
| B1     | 0.179             |
| B2     | **0.291**         |
| B3     | 0.190             |
| E+     | 0.088             |

**Conclusion:** E+ has the **lowest docking quality** across all methods. Its mean DockQ of 0.088 is firmly in the "Incorrect" category and represents a **51% degradation** relative to vanilla B1 (0.179). The region-specific beta-scaling significantly distorts the predicted structure. B2 leads with 0.291, benefiting from targeted contact restraints.

---

### 5. DockQ Quality Distribution (`dockq_distribution.png`)

**Methodology:** Count of complexes falling into each DockQ quality category (Incorrect / Acceptable / Medium / High) using the global DockQ score. Grouped bar chart per method.

**Results (e_plus_4):**
| Category    | B1 | B2 | B3 | E+ |
|-------------|----|----|----|----|
| Incorrect   | 17 |  9 | 14 | **37** |
| Acceptable  | 20 | 18 | 13 |  9 |
| Medium      | 10 | 19 |  9 |  0 |
| High        |  0 |  0 |  0 |  0 |

**Conclusion:** 37 out of 46 E+ predictions (80%) are classified as "Incorrect" by DockQ, with only 9 reaching "Acceptable" and none reaching "Medium". In contrast, B1 achieves 30 Acceptable-or-better out of 47, and B2 achieves 37 out of 46. This confirms that E+'s beta-scaling **severely degrades structural quality**.

---

### 6. E+ vs B1: Docking Quality Scatter (`eplus_vs_b1_dockq.png`)

**Methodology:** Per-complex scatter plot of DockQ (AB+AC) with B1 on x-axis and E+ on y-axis. Points above the y=x diagonal indicate E+ outperforms B1. Each point is labeled with the complex name.

**Results:**
- Nearly all points fall **below the diagonal**, confirming systematic degradation by E+.
- Head-to-head: E+ wins on 14 complexes, B1 wins on 16, others are close.
- The few complexes where E+ does well (e.g., 8HGM_CDB, 8HES_HLC) have low B1 DockQ to begin with.
- For complexes where B1 achieves good docking (DockQ > 0.5), E+ consistently achieves lower scores (e.g., 8CDD_EDB: B1=0.75, E+=0.01; 8CDE_DCB: B1=0.78, E+=0.08).

**Conclusion:** E+ does not improve docking for any well-docked complex. Its occasional "wins" occur only when B1 also fails. The beta-scaling perturbation is destructive to structure quality.

---

### 7. E+ vs B1: Epitope F1 Scatter (`eplus_vs_b1_f1.png`)

**Methodology:** Per-complex scatter of Epitope F1 (B1 on x-axis, E+ on y-axis). E+'s F1 uses the optimal threshold from the heatmap; B1's F1 uses structure-derived epitope contacts.

**Results:**
- Points are distributed roughly symmetrically around the diagonal.
- Head-to-head: E+ wins on 20 complexes, B1 wins on 21, 5 ties (within 0.01).
- Notable E+ advantages: 8HGM_CDB (E+=0.80, B1=0.52), 8IVA_CEG (E+=0.71, B1=0.53), 8IV5_ABG (E+=0.60, B1=0.18).
- Notable B1 advantages: 8BYU_HLA (B1=0.65, E+=0.07), 8CDD_EDB (B1=0.95, E+=0.30).
- Several complexes cluster in the bottom-left corner (both methods fail).

**Conclusion:** E+ and B1 are **roughly equivalent** on epitope F1 overall. E+ achieves better F1 on some complexes where B1 mispositions the antibody, but fails on others. This 50/50 split is not a meaningful advantage for E+.

---

## Summary Tables

### Axis A: Epitope Prediction

| Method | AUC-ROC | AUC-PR | F1 (optimal) | MCC   | Precision | Recall |
|--------|---------|--------|--------------|-------|-----------|--------|
| B1     | --      | --     | 0.345        | 0.269 | 0.357     | 0.345  |
| B2 (best) | --   | --     | **0.492**    | **0.436** | **0.508** | **0.492** |
| B3 (best) | --   | --     | 0.347        | 0.269 | 0.366     | 0.338  |
| **E+** | 0.645   | 0.304  | 0.340        | 0.238 | 0.321     | 0.423  |

### Axis B: Docking Quality

| Method | DockQ (AB+AC) | DockQ (global) | iRMSD (AB+AC) | fnat (AB+AC) | Ab-aligned Ag RMSD |
|--------|---------------|----------------|---------------|--------------|---------------------|
| B1     | 0.179         | 0.338          | 11.46         | 0.169        | 30.98               |
| B2 (best) | **0.291**  | **0.416**      | **7.72**      | **0.286**    | **24.86**           |
| B3 (best) | 0.190      | 0.343          | 11.22         | 0.171        | 29.53               |
| **E+** | 0.088         | 0.154          | 14.33         | 0.124        | 32.19               |

### Axis C: Confidence Metrics

| Method | confidence_score | iPTM  | complex_pLDDT | complex_ipLDDT |
|--------|-----------------|-------|---------------|----------------|
| B1     | 0.883           | 0.737 | 0.919         | 0.895          |
| B2     | 0.908           | 0.834 | 0.927         | 0.905          |
| B3     | 0.894           | 0.794 | 0.919         | 0.897          |
| **E+** | **0.520**       | **0.661** | **0.485** | **0.478**      |

E+ has significantly lower confidence scores across all metrics, reflecting the distortion introduced by region-specific beta-scaling during diffusion.

---

## Overall Conclusions

### Did E+ achieve its primary goal (epitope discovery)?

**No.** The E+ heatmap achieves a mean AUC-ROC of 0.645 (target was > 0.7) and AUC-PR of 0.304 (target was > 0.5). When binarized at the optimal threshold, E+ achieves an epitope F1 of 0.340 -- statistically equivalent to simply extracting contacts from a vanilla Boltz2 prediction (B1 F1 = 0.345). The heatmap provides some signal above random, but does not provide a meaningful advantage over using the standard structure prediction output.

### Did E+ maintain structure quality?

**No.** E+ introduces substantial structural degradation. Its mean DockQ (0.088) is less than half of vanilla Boltz2 (0.179) and 80% of predictions are classified as "Incorrect." The confidence scores (0.520 vs 0.883) also confirm the model is less certain about E+ structures.

### Why does E+ underperform?

Several likely factors:

1. **Beta-scaling magnitude may be insufficient.** The emphasis (0.5) and de-emphasis (-0.3) values may not create enough differential signal between regions to produce discriminative contact patterns.

2. **Region partitioning is too coarse.** With only 10 regions and 50% overlap for long antigens, each "region" covers a large fraction of the antigen surface, reducing spatial resolution.

3. **Diffusion model resilience.** The structure module's sampling may partially "undo" the pair representation biases, averaging out the region-specific steering signal.

4. **Contact threshold too permissive.** Using 8.0 A for contact detection during scanning may capture too many spurious contacts, diluting the true signal.

5. **Single-step steering.** Beta-scaling is applied only to the diffusion conditioning step, while the trunk (which contributes most of the structural knowledge) is shared across all regions.

### Suggested Modifications

1. **Increase beta magnitudes.** Try emphasis values of 1.0-2.0 and de-emphasis of -0.5 to -1.0 to create stronger differential signal. Run a sweep of beta values.

2. **Increase region count.** Use 20-40 regions instead of 10 to improve spatial resolution, especially for large antigens. The single-trunk speedup makes this practically free.

3. **Tighten contact threshold.** Reduce from 8.0 A to 5.0-6.0 A during scanning to match the ground truth epitope definition and reduce false positives.

4. **Multi-step steering.** Apply beta-scaling not just in diffusion conditioning but also during recycling/pairformer steps to create a stronger structural bias.

5. **Weighted aggregation.** Instead of binary contact counting, weight contacts by confidence score or by CDR-antigen distance to produce a more calibrated heatmap.

6. **Per-region structure assessment.** Instead of extracting contacts from each region scan, compute a per-region confidence metric and use it to weight the heatmap contribution.

7. **Consensus with baseline.** Use E+ heatmap as a prior to re-weight contacts from B1's vanilla prediction, combining the structural accuracy of B1 with the scanning diversity of E+.
