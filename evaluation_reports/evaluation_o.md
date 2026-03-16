# Evaluation Report: O+ (Asymmetric Beta-Scaling v2)

## 1. Feature Goal

The **O+ (Asymmetric Beta-Scaling v2)** feature encodes a well-established biological asymmetry in antibody-antigen recognition: CDR-H3 (heavy chain CDR3) dominates antigen binding, contributing approximately 60% of contacts, while CDR-L3 (light chain CDR3) plays a supporting role. The feature applies differential denoising step scaling during Boltz2 diffusion sampling — `beta_H = 0.4` for CDR-H3 and `beta_L = 0.1` for CDR-L3 — so that heavy chain loop atoms receive a larger update per diffusion step (`step * 1.4`) compared to light chain atoms (`step * 1.1`).

**Expected outcomes:**
1. H3/L3 antigen contact ratio closer to the biological target of 1.5–2.0
2. Improved docking quality (DockQ) at the antigen-antibody interface
3. Better CDR-H3 loop accuracy (lower RMSD) without degrading CDR-L3

## 2. Evaluation Setup

### Methods compared
- **asymmetric_beta (O+):** The new feature (47 complexes, 5 models each)
- **baseline:** Unmodified Boltz2 (47 complexes, 5 models each)
- **contact_restraints:** Boltz2 with explicit contact restraints from known interface residues (46 complexes, multiple restraint variants, 5 models each)
- **pocket_restraints:** Boltz2 with pocket-based distance restraints (36 complexes, multiple variants, 5 models each)

### Evaluation runs
- **o_plus_2:** O+ vs. baseline only (2 methods)
- **o_plus_4:** O+ vs. all three baselines (4 methods)

Both runs used the corrected `scripts/eval_asymmetric_beta/evaluate_all.py` script. Interface quality was assessed with DockQ v2 (producing DockQ, iRMSD, LRMSD, fnat, F1). CDR-specific metrics were computed with BioPython-based framework-aligned RMSD and heavy-atom contact counting (<5 A cutoff).

### Model selection strategies
The evaluation reports three selection strategies:
- **Oracle:** Best model per complex (upper bound on performance)
- **Average:** Mean across all 5 models per complex
- **Top-1:** Best model selected by Boltz2 confidence score (iptm)

All plots use oracle selection unless otherwise noted. For methods with multiple restraint variants (contact_restraints, pocket_restraints), the best variant per complex is selected first.

---

## 3. Plot-by-Plot Analysis: o_plus_2 (O+ vs. Baseline)

### Figure 1 — Overall DockQ Comparison (Bar Plot)

**How produced:** For each method, DockQ scores across all interfaces within each model are averaged to produce a single per-model score. The **oracle** (best model per complex) is selected. The bar height is the mean across all 47 complexes; error bars show 95% bootstrap confidence intervals.

**Result:** O+ achieves DockQ = 0.392 [0.341, 0.445] and baseline achieves DockQ = 0.391 [0.337, 0.448]. The bars are virtually identical in height and the confidence intervals overlap almost completely.

**Conclusion:** O+ produces no meaningful improvement over baseline in overall docking quality. The +0.001 DockQ difference is well within statistical noise.

---

### Figure 2 — Per-Interface DockQ Comparison (3 Subplots)

**How produced:** Same oracle selection methodology as Figure 1, but DockQ is reported separately for each chain-pair interface (A:B, A:C, B:C). Oracle is selected per interface independently.

**Subplot A:B (Antigen–Heavy chain):** Both methods show DockQ ~0.259 [0.18, 0.34]. This is the most challenging interface involving direct antigen docking.

**Subplot A:C (Antigen–Light chain):** Both methods show DockQ ~0.247 [0.17, 0.33]. Slightly lower than A:B, consistent with typically weaker antigen-light chain contacts.

**Subplot B:C (Heavy–Light chain):** Both methods show DockQ ~0.685 [0.66, 0.71]. Much higher, reflecting that VH-VL pairing is well predicted by Boltz2's baseline capability.

**Conclusion:** O+ does not improve any specific interface. The antigen-antibody interfaces (A:B, A:C) are poorly predicted (DockQ ~0.25) while the intra-antibody interface (B:C) is well modeled (DockQ ~0.69), indicating that antigen docking — not antibody folding — is the primary challenge. O+ does not address this bottleneck.

---

### Figure 3 — CDR RMSD Comparison (2 Subplots, Boxplots)

**How produced:** Predictions are aligned to ground truth on **framework (non-CDR) CA atoms** using BioPython Superimposer. RMSD is then computed on CDR CA atoms only. The oracle model (lowest H3 RMSD per complex) is selected. Boxplots show the distribution across all 47 complexes.

**Left subplot — CDR-H3 RMSD:**
- O+ median: ~2.05 A, IQR: ~1.2–2.95 A
- Baseline median: ~2.11 A, IQR: ~1.2–3.1 A
- Mean oracle: O+ = 2.24 A, baseline = 2.27 A (delta = -0.023 A)

**Right subplot — CDR-L3 RMSD:**
- O+ median: ~1.02 A, IQR: ~0.7–1.5 A
- Baseline median: ~0.98 A, IQR: ~0.7–1.5 A
- Mean oracle: O+ = 1.22 A, baseline = 1.21 A (delta = +0.014 A)
- Both show outliers around 3.0–3.7 A

**Conclusion:** There is no meaningful difference in CDR loop accuracy. The delta for H3 RMSD (0.023 A in O+'s favor) and L3 RMSD (0.014 A in baseline's favor) are both negligible relative to the ~2.2 A mean and the high variance. The feature does not improve CDR-H3 loop modeling.

---

### Figure 4 — H3/L3 Antigen Contact Ratio (Boxplot)

**How produced:** For the oracle model (lowest H3 RMSD per complex), heavy-atom contacts within 5 A between each CDR and the antigen chain are counted, and the ratio H3_contacts / L3_contacts is computed. Complexes where L3_contacts = 0 are excluded. Reference lines at 1.5 and 2.0 mark the biologically expected H3 dominance ratio.

**Result:**
- O+ median ratio: ~3.0, IQR: ~1–5, mean: 5.39, outliers up to 42
- Baseline median ratio: ~3.2, IQR: ~1–7, mean: 6.65, outliers up to 44

**Conclusion:** Both methods already produce H3-dominant predictions (ratio >> 1), but with extreme variance driven by cases where L3 makes very few contacts. The biological target of 1.5–2.0 is far below both distributions. O+ shows a slightly lower mean ratio (5.39 vs 6.65) and a tighter IQR, suggesting marginally less extreme H3 over-representation, but the wide spread makes this difference not statistically meaningful. The feature does not shift the contact ratio toward the target range.

---

### Figure 5 — Per-Complex DockQ Scatter (O+ vs. Baseline)

**How produced:** Each dot represents one complex. X-axis = baseline oracle DockQ (mean across interfaces, best model per complex). Y-axis = O+ oracle DockQ. The dashed line marks y=x (parity).

**Result:** The 47 points cluster tightly along the y=x diagonal, spanning from ~0.21 to ~0.80. Most points lie very close to parity, with small deviations in both directions. A few complexes in the mid-range (DockQ ~0.3–0.5) show slightly larger deviations above or below the line, but no systematic trend is visible.

**Conclusion:** O+ and baseline produce nearly equivalent docking quality on a per-complex basis. There is no subset of complexes where O+ consistently outperforms baseline. The scatter confirms that the feature's effect is negligible across the entire difficulty spectrum — it neither helps easy cases nor difficult ones.

---

## 4. Plot-by-Plot Analysis: o_plus_4 (O+ vs. Three Baselines)

### Figure 6 — Overall DockQ Comparison (Bar Plot)

**How produced:** Same oracle methodology as Figure 1. For `contact_restraints` and `pocket_restraints`, which have multiple restraint variants per complex, the script first selects the **best variant** (the one whose best model achieves the highest DockQ) before computing the oracle. Bar heights are means across complexes with 95% bootstrap CIs.

**Result:**

| Method | DockQ (oracle) | 95% CI | n |
|--------|---------------|--------|---|
| asymmetric_beta | 0.392 | [0.337, 0.449] | 47 |
| baseline | 0.391 | [0.341, 0.448] | 47 |
| contact_restraints | **0.504** | [0.443, 0.567] | 46 |
| pocket_restraints | 0.396 | [0.331, 0.468] | 36 |

**Conclusion:** Contact restraints are the clear winner (DockQ 0.504), outperforming all other methods by a large margin (+0.112 over O+). O+ and baseline are virtually indistinguishable (delta = +0.001). Pocket restraints (0.396) show a marginal, non-significant improvement over baseline. This demonstrates that explicit structural guidance about the interface is far more effective than implicit step-size modulation.

---

### Figure 7 — Per-Interface DockQ Comparison (3 Subplots)

**How produced:** Same as Figure 2 but with all 4 methods displayed side-by-side.

**Subplot A:B (Antigen–Heavy chain):** Contact restraints (0.419) substantially outperform baseline/O+ (~0.259) and pocket restraints (0.264). This is the interface where contact restraints provide direct structural guidance.

**Subplot A:C (Antigen–Light chain):** Similar pattern — contact restraints (0.412) dominate, while the other three methods cluster near 0.245–0.259.

**Subplot B:C (Heavy–Light chain):** All four methods perform similarly (0.68–0.70), confirming that VH-VL pairing quality is an inherent Boltz2 capability not affected by any steering method.

**Conclusion:** Contact restraints improve both antigen-antibody interfaces (A:B and A:C) by ~0.16 DockQ units over O+/baseline. O+, baseline, and pocket restraints are indistinguishable on all three interfaces. The heavy-light chain interface is well modeled regardless of method.

---

### Figure 8 — CDR RMSD Comparison (2 Subplots, Boxplots)

**How produced:** Same framework-aligned methodology as Figure 3, now showing all 4 methods.

**Left subplot — CDR-H3 RMSD:**

| Method | Mean (oracle) | Median |
|--------|--------------|--------|
| asymmetric_beta | 2.24 A | ~2.05 A |
| baseline | 2.27 A | ~2.10 A |
| contact_restraints | **1.81 A** | ~1.55 A |
| pocket_restraints | 1.96 A | ~2.00 A |

Contact restraints achieve the lowest H3 RMSD (1.81 A mean, ~20% improvement over baseline). O+ and baseline are virtually identical and worse than both restraint-based methods.

**Right subplot — CDR-L3 RMSD:**
All four methods perform similarly (mean 1.2–1.25 A), with no meaningful separation. This is expected since no method specifically targets L3 structural accuracy.

**Conclusion:** O+ provides no improvement in CDR-H3 loop accuracy. Contact restraints, which provide explicit interface guidance, achieve a meaningful H3 RMSD reduction (~0.46 A or ~20% improvement). The fact that L3 RMSD is similar across all methods indicates that CDR-L3 modeling is not a bottleneck.

---

### Figure 9 — H3/L3 Antigen Contact Ratio (Boxplot)

**How produced:** Same as Figure 4 but with all 4 methods. Oracle selection per complex, H3/L3 contact count ratio with reference lines at 1.5 and 2.0.

**Result:**

| Method | Mean H3/L3 Ratio | Median |
|--------|------------------|--------|
| asymmetric_beta | 5.39 | ~3.0 |
| baseline | 6.65 | ~3.2 |
| contact_restraints | 4.69 | ~3.0 |
| pocket_restraints | 4.45 | ~3.0 |

All methods show median ratios between ~2.5–4.5, well above the biological target of 1.5–2.0. Baseline has the highest variance with outliers extending to ~48. O+ and contact restraints show similar distributions. Pocket restraints have the tightest spread.

**Conclusion:** All methods inherently produce H3-dominant predictions far exceeding the biological target. The feature was designed to encode H3 dominance, but the model already exhibits this bias naturally. O+ does not meaningfully modulate the contact ratio toward the desired 1.5–2.0 range.

---

### Figure 10 — Per-Complex DockQ Scatter (O+ vs. Baseline)

**How produced:** Same as Figure 5 — each dot is one complex, plotting O+ oracle DockQ (y-axis) against baseline oracle DockQ (x-axis).

**Result:** Points cluster tightly along the y=x diagonal, with small deviations in both directions. No systematic bias above or below the line is visible. A few mid-range complexes show the largest deviations, but these go both ways (some favor O+, some favor baseline).

**Conclusion:** O+ and baseline are equivalent on a per-complex basis. The feature does not selectively help any subset of complexes (e.g., difficult vs. easy, or particular antigen types).

---

## 5. Statistical Tests

### o_plus_2 (O+ vs. Baseline)

| Comparison | n | Mean Delta DockQ | 95% CI | Wilcoxon p | Bonferroni p |
|-----------|---|-----------------|--------|-----------|-------------|
| O+ vs. baseline | 47 | +0.001 | [-0.009, 0.017] | 0.434 | 1.0 |

Not significant. The delta of +0.001 is negligible, and the CI spans zero. The Wilcoxon signed-rank test (p = 0.43) confirms no paired difference.

### o_plus_4 (O+ vs. Three Baselines)

| Comparison | n | Mean Delta DockQ | 95% CI | Wilcoxon p | Bonferroni p |
|-----------|---|-----------------|--------|-----------|-------------|
| O+ vs. baseline | 47 | +0.001 | [-0.010, 0.016] | 0.434 | 1.0 |
| O+ vs. pocket_restraints | 36 | -0.007 | [-0.015, -0.0003] | 0.203 | 0.609 |
| **O+ vs. contact_restraints** | **46** | **-0.108** | **[-0.156, -0.061]** | **5.6e-9** | **1.7e-8** |

O+ is **significantly worse** than contact restraints (delta = -0.108 DockQ, p < 2e-8 after Bonferroni correction). O+ is not significantly different from pocket restraints (p = 0.61) or baseline (p = 1.0).

---

## 6. Overall Conclusions

### Did the feature achieve its goal?

**No.** The Asymmetric Beta-Scaling v2 (O+) feature did not achieve measurable improvement over the baseline across any evaluated metric:

1. **DockQ is unchanged:** O+ oracle = 0.392 vs. baseline = 0.391 (delta = +0.001, p = 0.43). The per-complex scatter confirms equivalence across the full difficulty range.
2. **CDR-H3 RMSD is unchanged:** O+ = 2.24 A vs. baseline 2.27 A (delta = 0.023 A, negligible).
3. **CDR-L3 RMSD is unchanged:** O+ = 1.22 A vs. baseline 1.21 A (delta = 0.014 A, negligible).
4. **H3/L3 contact ratio is not shifted** toward the biological target. Both methods already over-represent H3 contacts (mean ratio ~5–7x vs. target 1.5–2.0x).
5. **Contact restraints dramatically outperform all methods** (DockQ 0.504 vs. ~0.392, H3 RMSD 1.81 A vs. ~2.25 A), demonstrating that explicit interface guidance is far more effective.

### Why did O+ not produce improvement?

The feature scales the denoising step by `(1 + beta)` uniformly across CDR atoms — `1.4x` for H3, `1.1x` for L3. This is a subtle perturbation that:

1. **Modifies step magnitude, not direction:** Larger steps make CDR-H3 atoms move more per diffusion step, but provide no information about *where* they should move (i.e., toward the antigen). The model's learned denoising direction dominates.
2. **Is too weak:** A 40% step increase is mild relative to the many diffusion steps. The cumulative effect is absorbed by the model's self-correction in subsequent steps.
3. **Has no antigen awareness:** The scaling is applied uniformly to CDR atoms regardless of their spatial relationship to the antigen, unlike contact restraints which encode explicit interface knowledge.
4. **Encodes an already-present bias:** The model naturally produces H3-dominant contact patterns (ratio ~3–7x). Adding a mild H3 emphasis does not address the actual bottleneck, which is accurate antigen docking orientation.

### Recommended Modifications

1. **Increase beta values significantly:** Test beta_H = 2.0–5.0 (3x–6x step size) to determine if there is any dose-response relationship. The current 1.4x may be below a meaningful threshold.

2. **Add directional bias:** Instead of uniformly scaling the denoising step, project the CDR-H3 update onto the antigen-CDR vector and amplify only that directional component. This combines *where* (toward antigen) with *how much* (amplified step).

3. **Apply noise-level-dependent scaling:** Use stronger scaling in early diffusion steps (when global structure is set) and weaker scaling later (when fine details are resolved), since CDR-H3 positioning is a global structural decision.

4. **Combine with contact restraints:** Use asymmetric beta scaling as a complement to contact restraints rather than a standalone method — e.g., apply step scaling only within the context of an active contact potential.

5. **Target the update direction, not magnitude:** Instead of scaling `step * (1 + beta)`, project the denoising update onto the antigen-CDR vector and amplify only that component.

6. **Explore embedding-space steering (Idea V):** Per the EmbedOpt-inspired proposals, optimizing CDR3 representations in Boltz2's embedding/pair representation space may be more robust than coordinate-space step scaling, especially for novel sequences outside the training distribution.
