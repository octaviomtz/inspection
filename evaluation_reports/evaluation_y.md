# Hierarchical Steering (Idea Y) — Evaluation Report

## 1. Feature Goal

Hierarchical steering aims to improve antibody-antigen docking in Boltz-2 by combining two complementary strategies across the diffusion timeline:

1. **Early diffusion (embedding space):** β-scaling of pair representations steers the model toward the correct binding interface before atomic coordinates have converged.
2. **Late diffusion (coordinate space):** Gradient-based potentials (contact, pocket, clash) refine the docked pose once coordinates become meaningful.

The hypothesis is that embedding-space guidance in early steps provides a better "starting point" for the coordinate-space potentials that activate later, yielding better overall docking quality than either strategy alone or baseline Boltz-2.

## 2. Evaluation Setup

- **Benchmark:** 48 antibody-antigen complexes from `examples/cdrs.csv` (47 successfully evaluated).
- **Ground truth:** Energy-minimized PDB structures in `pdb_minimized/`.
- **Chain convention:** A = antigen, B = heavy chain, C = light chain.
- **Model selection:** Top-1 by `confidence_score` (≈ 0.2×ptm + 0.8×iptm); oracle = best DockQ among 5 samples.
- **Multi-setting baselines (contact_restraints, pocket_guided):** Best-of-settings aggregation — each complex was run with multiple restraint/pocket configurations and the setting producing the highest top-1 total DockQ was kept.
- **Statistical tests:** Wilcoxon signed-rank (paired, non-parametric) on per-complex metrics; 95% bootstrap confidence intervals (10,000 resamples).
- **Script:** `scripts/eval/evaluate_hierarchical_steering.py`

### Comparisons

| Evaluation | Methods | n |
|---|---|---|
| `y_plus_2/` | baseline, hierarchical | 47 |
| `y_plus_4/` | baseline, contact_restraints, pocket_guided, hierarchical | 47 (36 for pocket_guided) |

---

## 3. Summary of Results

### 3.1. Head-to-head: Hierarchical vs Baseline (`y_plus_2/`)

| Metric | Baseline | Hierarchical | Δ | p-value |
|---|---|---|---|---|
| Top-1 Total DockQ (mean) | 0.338 | 0.364 | +0.026 | 0.224 |
| Top-1 Ab-Ag DockQ (mean) | 0.179 | 0.213 | +0.034 | 0.459 |
| Oracle Total DockQ (mean) | 0.391 | 0.403 | +0.012 | — |
| CDR-H3 RMSD (Å, mean) | 3.00 | 2.70 | −0.30 | **0.025** |
| Epitope F1 5Å (mean) | 0.345 | 0.348 | +0.003 | 0.544 |
| iptm (mean) | 0.738 | 0.726 | −0.011 | 0.446 |
| DockQ > 0.23 (%) | 63.8 | 63.8 | 0.0 | — |
| DockQ > 0.49 (%) | 21.3 | 27.7 | +6.4 | — |

### 3.2. Four-method comparison (`y_plus_4/`)

| Metric | Baseline | Contact Restraints | Pocket Guided | Hierarchical |
|---|---|---|---|---|
| Top-1 Total DockQ | 0.338 | **0.499** | 0.364 | 0.364 |
| Top-1 Ab-Ag DockQ | 0.179 | **0.408** | 0.217 | 0.213 |
| Oracle Total DockQ | 0.391 | **0.515** | 0.399 | 0.403 |
| CDR-H3 RMSD (Å) | 3.00 | **2.35** | 2.64 | 2.70 |
| Epitope F1 5Å | 0.345 | **0.587** | 0.368 | 0.348 |
| iptm | 0.738 | 0.742 | **0.773** | 0.726 |
| DockQ > 0.23 (%) | 63.8 | **95.7** | 66.7 | 63.8 |
| DockQ > 0.49 (%) | 21.3 | **54.3** | 27.8 | 27.7 |

---

## 4. Per-Plot Analysis

### 4.1. DockQ Bar Plot (`dockq_barplot.png`)

**How produced:** Three side-by-side bar charts. Each bar shows the mean across all complexes for a given method, with 95% bootstrap confidence intervals as error bars. Metrics: Top-1 Total DockQ (confidence-selected best model, all-interface average), Top-1 Ab-Ag DockQ (confidence-selected best model, antibody-antigen interface only), and Oracle Total DockQ (best DockQ model among 5 samples).

#### `y_plus_2/` (2 methods)

- **Top-1 Total DockQ:** Hierarchical (0.364) slightly exceeds baseline (0.338), but the 95% CIs overlap substantially ([0.289, 0.391] vs [0.309, 0.419]).
- **Top-1 Ab-Ag DockQ:** Hierarchical (0.213) vs baseline (0.179). The improvement is more visible here but CIs are very wide due to high variance, especially for hierarchical.
- **Oracle Total DockQ:** Nearly identical (0.391 vs 0.403), indicating that the pool of 5 samples contains similar best-case structures for both methods.

**Conclusion:** Hierarchical steering shows a consistent but small and statistically non-significant improvement across all three DockQ variants. The large confidence intervals reflect high inter-complex variance.

#### `y_plus_4/` (4 methods)

- **Top-1 Total DockQ:** Contact restraints (0.499) clearly separates from the pack. Baseline (0.338), pocket_guided (0.364), and hierarchical (0.364) are clustered together, with the latter two virtually identical.
- **Top-1 Ab-Ag DockQ:** Contact restraints (0.408) is roughly double baseline (0.179). Pocket_guided (0.217) and hierarchical (0.213) show marginal gains over baseline.
- **Oracle Total DockQ:** Contact restraints (0.515) again leads. The other three methods are tightly grouped around 0.39–0.40.

**Conclusion:** Contact restraints is the dominant method by a large margin. Hierarchical and pocket_guided show nearly identical performance and only marginal gains over baseline.

---

### 4.2. DockQ Box Plot (`dockq_boxplot.png`)

**How produced:** Box-and-whisker plot of per-complex top-1 Ab-Ag DockQ values. The box spans Q1–Q3 (interquartile range), the orange line is the median, whiskers extend to 1.5×IQR, and outliers are shown as individual points.

#### `y_plus_2/` (2 methods)

- **Baseline:** Median near 0.04, Q3 ≈ 0.21. A cluster of outliers above 0.60 represents the few complexes where baseline achieves good docking.
- **Hierarchical:** Median still near 0.04, but Q3 expands to ≈ 0.53. The upper whisker extends to ≈ 0.80 (no outliers beyond it), showing that hierarchical pushes more complexes into the medium-quality range.

**Conclusion:** Both methods have the same failure mode — about half the complexes result in near-zero DockQ (incorrect docking). However, hierarchical steering lifts the upper quartile: when it works, it works better. The wider IQR reflects more variable but occasionally much-improved results.

#### `y_plus_4/` (4 methods)

- **Baseline:** Lowest IQR, median ≈ 0.04, most mass near zero.
- **Contact restraints:** Median ≈ 0.41, by far the highest. Q1 ≈ 0.08, Q3 ≈ 0.67. The majority of complexes achieve acceptable or better docking.
- **Pocket guided:** Median ≈ 0.04, but Q3 ≈ 0.50, similar spread to hierarchical.
- **Hierarchical:** Median ≈ 0.04, Q3 ≈ 0.53, closely mirroring pocket_guided.

**Conclusion:** Contact restraints is the only method that moves the median above near-zero. Hierarchical and pocket_guided extend the upper quartile similarly but cannot shift the median — too many complexes remain unsolved.

---

### 4.3. Per-Complex Scatter Plot (`dockq_scatter_hierarchical_vs_baseline.png`)

**How produced:** Each point represents one complex, with baseline Ab-Ag DockQ on the x-axis and hierarchical Ab-Ag DockQ on the y-axis (both top-1, confidence-selected). The dashed diagonal is y = x; points above it indicate hierarchical improvement.

#### `y_plus_2/` and `y_plus_4/` (identical scatter — same paired data)

- **Dense cluster at origin (0, 0):** ~20 complexes where both methods fail completely (DockQ < 0.05). These are inherently difficult targets that neither method can solve.
- **Points above diagonal:** Several complexes show large hierarchical gains, notably:
  - One complex jumps from baseline ≈ 0.22 to hierarchical ≈ 0.78 (the most dramatic improvement).
  - Another from ≈ 0.20 to ≈ 0.66.
  - Several moderate gains in the 0.50–0.75 range.
- **Points below diagonal:** One notable regression from ≈ 0.65 to ≈ 0.07 (large loss). A few smaller regressions exist.
- **High-quality region (> 0.60):** Points cluster near the diagonal, suggesting hierarchical preserves quality when baseline already performs well.

**Conclusion:** Hierarchical steering produces a few dramatic wins (large above-diagonal offsets) and one significant regression. The win count (25) exceeds the loss count (18), but the cluster at zero dominates. The feature helps most on "medium-difficulty" complexes where baseline partially succeeds; it cannot rescue the hardest cases.

---

### 4.4. CDR RMSD Bar Plot (`cdr_rmsd_barplot.png`)

**How produced:** Grouped bar chart showing mean Cα RMSD per CDR loop (H1, H2, H3, L1, L2, L3) after framework alignment (superposition on non-CDR antibody residues). Error bars show standard deviation across complexes. Lower is better.

#### `y_plus_2/` (2 methods)

- **CDR-H3:** Largest improvement — baseline 3.00 Å vs hierarchical 2.70 Å (−0.30 Å). This is the only CDR with a statistically significant difference (p = 0.025 by Wilcoxon signed-rank). CDR-H3 is also the most variable loop (std ≈ 1.5 Å), consistent with it being the longest and most flexible CDR.
- **Other CDRs (H1, H2, L1, L2, L3):** All show small improvements (0.02–0.06 Å) for hierarchical, but none are statistically significant. These loops have lower baseline RMSD (0.95–1.35 Å) and less room for improvement.

**Conclusion:** Hierarchical steering significantly improves CDR-H3 loop modeling (p = 0.025), the most clinically relevant CDR loop for binding specificity. This is arguably the most important positive finding — CDR-H3 conformation is critical for antigen recognition and is notoriously difficult to predict.

#### `y_plus_4/` (4 methods)

- **CDR-H3:** Contact restraints achieves the best RMSD (2.35 Å), followed by pocket_guided (2.64 Å) and hierarchical (2.70 Å), all improving over baseline (3.00 Å).
- **Other CDRs:** All four methods perform similarly (within 0.05–0.10 Å), with contact restraints having a very slight edge across the board.

**Conclusion:** For CDR-H3, contact restraints leads, but hierarchical still delivers meaningful improvement over baseline. The non-H3 loops are well-modeled by all methods and show minimal differentiation.

---

### 4.5. Epitope F1 Bar Plot (`epitope_f1_barplot.png`)

**How produced:** Grouped bar chart showing mean F1 score for epitope residue prediction at two distance thresholds: 5Å (heavy atom contacts) and 8Å (Cα-Cα contacts). Error bars show standard deviation. A residue is predicted as an epitope residue if it contacts an antibody residue within the threshold in the predicted structure. Higher is better.

#### `y_plus_2/` (2 methods)

- **Epitope F1 (5Å):** Baseline 0.345 vs hierarchical 0.348 — virtually identical.
- **Epitope F1 (8Å):** Baseline 0.311 vs hierarchical 0.301 — also virtually identical, with hierarchical marginally lower.
- Both methods show very large standard deviations (error bars extend to near 0), reflecting high inter-complex variance.

**Conclusion:** Hierarchical steering does not improve epitope prediction over baseline. The 5Å Wilcoxon test (p = 0.544) confirms no significant difference. This suggests the steering does not systematically improve the predicted binding site location.

#### `y_plus_4/` (4 methods)

- **Contact restraints** dominates: 0.587 at 5Å, 0.552 at 8Å — nearly double the other methods. This is expected since contact restraints explicitly provide interface information.
- Baseline (0.345), pocket_guided (0.368), and hierarchical (0.348) are tightly grouped at both thresholds.

**Conclusion:** Explicit contact restraints dramatically improve epitope prediction. The information-free methods (baseline, hierarchical) and weakly-informed method (pocket_guided) cannot match this. Hierarchical steering adds no benefit for epitope identification.

---

### 4.6. Confidence Calibration Scatter (`confidence_vs_dockq.png`)

**How produced:** Scatter plot with predicted iptm on the x-axis and actual Ab-Ag DockQ on the y-axis. Each point is one complex, colored by method. A well-calibrated model would show strong positive correlation (high iptm → high DockQ).

#### `y_plus_2/` (2 methods)

- **General pattern:** Poor calibration. Many complexes with high iptm (0.80–0.95) have near-zero DockQ (false confidence). Conversely, some complexes with moderate iptm (0.55–0.65) achieve high DockQ (0.70+).
- **Baseline (gray):** Spread across iptm range 0.35–0.97 with no clear trend.
- **Hierarchical (green):** Similar spread and similar lack of correlation. iptm is slightly lower on average (0.726 vs 0.738) but this is not significant (p = 0.446).

**Conclusion:** Boltz-2's confidence metric (iptm) is poorly calibrated for docking quality in both methods. The model is often confident about incorrect poses and sometimes less confident about correct ones. Hierarchical steering does not improve or degrade confidence calibration.

#### `y_plus_4/` (4 methods)

- **Contact restraints (orange):** Shows the most points in the upper-right region (high iptm and high DockQ), but also has many false-confident cases. iptm mean (0.742) is similar to baseline.
- **Pocket guided (pink):** Highest mean iptm (0.773), but this doesn't translate to better DockQ — many high-iptm points sit at DockQ ≈ 0.
- **Hierarchical (green) and baseline (gray):** Intermixed, confirming similar calibration behavior.

**Conclusion:** iptm is unreliable as a docking quality predictor across all methods. Pocket_guided inflates confidence without improving DockQ. The scatter pattern suggests iptm primarily reflects chain-internal folding quality rather than interface accuracy.

---

## 5. Overall Conclusions

### What hierarchical steering achieves

1. **Statistically significant CDR-H3 improvement** (p = 0.025): The most important CDR loop for antigen recognition is modeled 0.30 Å better on average. This is a meaningful gain for the loop that drives binding specificity.

2. **Trend toward better DockQ:** +0.026 top-1 total DockQ, +0.034 Ab-Ag DockQ, and a 6.4 percentage-point increase in medium-quality predictions (DockQ > 0.49: 21.3% → 27.7%). While not statistically significant at α = 0.05, the direction is consistent.

3. **No degradation:** The method does not significantly hurt any metric compared to baseline. Confidence scores remain in the same range.

### What hierarchical steering does not achieve

1. **No statistically significant DockQ improvement** (p = 0.224): The overall docking quality gain is too small relative to inter-complex variance.

2. **No epitope prediction improvement:** The binding interface is not more accurately identified.

3. **Far below contact restraints:** Contact restraints improve DockQ by +0.161 over baseline (p = 4.9×10⁻¹¹), while hierarchical achieves only +0.026. This is expected — contact restraints supply ground-truth interface information, while hierarchical steering uses no oracle information.

4. **Cannot rescue the hardest cases:** ~20 complexes (≈40%) remain at near-zero DockQ regardless of method. These likely fail at a fundamental level (incorrect global orientation) that local potentials cannot fix.

5. **Performance equivalent to pocket_guided:** Across all metrics, hierarchical and pocket_guided are statistically indistinguishable (all p > 0.37), despite using different mechanisms.

### Is the goal achieved?

**Partially.** The hierarchical approach achieves its primary mechanistic goal — the embedding-space early guidance followed by coordinate-space refinement produces measurably better CDR-H3 conformations and a consistent (though non-significant) trend toward better docking. However, the magnitude of improvement is modest, and the feature falls short of delivering a transformative gain in overall docking quality.

---

## 6. Suggested Modifications

### 6.1. Stronger embedding-space signal

The current β-scaling blending uses a linear schedule that smoothly transitions from steered to neutral conditioning. The effect may be too diluted. Options:
- Increase the β-scaling factor from the current value.
- Extend the embedding-space phase to overlap more with the coordinate-space phase.
- Use a sharper (step-like) transition rather than smooth blending.

### 6.2. Better potentials for the late phase

The coordinate-space potentials currently activate after the embedding phase ends. The potentials may not be strong enough or well-tuned:
- Increase `guidance_weight` and `resampling_weight` in the late phase.
- Add an attractive potential that explicitly pulls the antibody toward the antigen's predicted epitope region.
- Consider using DockQ-like interface loss as a differentiable potential.

### 6.3. Incorporate epitope information

The hierarchical method currently uses no external interface information. Adding even partial epitope hints (e.g., predicted from sequence) could provide the early-phase conditioning with a concrete binding site target, bridging the gap with contact restraints without requiring ground-truth contacts.

### 6.4. Adaptive scheduling per complex

The fixed diffusion timestep thresholds (t_activate, transition window) may not be optimal for all complexes. Easy complexes may benefit from less steering (to avoid disrupting already-good predictions), while hard complexes may need more aggressive intervention.

### 6.5. Focus on the failure cluster

The ~40% of complexes at DockQ ≈ 0 represent the biggest opportunity. These likely have incorrect global binding orientation. A coarse-grained rigid-body search or docking pre-step before diffusion could help place the antibody in the correct hemisphere around the antigen.

### 6.6. Combine with contact restraints

Given that contact restraints are vastly superior, a hybrid approach that uses hierarchical steering's embedding-space guidance together with contact restraint potentials in the late phase could potentially outperform either alone. The embedding phase could pre-orient the complex, allowing the contact restraints to refine from a better starting point.
