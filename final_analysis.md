# Final Analysis: Strategy L+ — CDR3 Beta-Scaling v2

**Date**: 2026-03-30
**Complexes evaluated**: 47 antibody-antigen complexes
**Ground truth**: minimized crystal structures (`pdb_minimized/`)
**Script**: `evaluate_cdr3beta_v2.py`

---

## 1. Feature Goal and Motivation

Strategy L (Round 1) applied a uniform beta=0.3 scaling to CDR3-CDR3 token pair representations in the DiffusionConditioning module. It successfully increased CDR3 conformational diversity (ensemble mean RMSD: 1.51→1.66 Å; fraction of diverse ensembles: 45%→60%), demonstrating that attention-bias scaling can actively steer CDR3 loop exploration during diffusion. However, it came with two significant costs: **CDR3 RMSD significantly worsened** (p=0.0008), meaning more diverse but less accurate loops, and **iPTM inflated sharply** (+0.091, p<0.0001) without corresponding structural improvement — a dangerous confidence-accuracy mismatch.

Strategy L+ v2 was designed to fix these problems through four targeted improvements:

| Improvement | Goal |
|---|---|
| **L+.1** Reduce beta: 0.3 → 0.12 | Eliminate CDR3 RMSD regression and iPTM inflation |
| **L+.2** Asymmetric betas: H3=0.15, L3=0.05 | Focus diversity on CDR-H3 (dominant binding loop), protect CDR-L3 geometry |
| **L+.3** CDR3-antigen interface scaling: `cdr3_antigen_beta=0.15` | Attract CDR3 loops toward the antigen, improving docking interface quality |
| **L+.4** Time-dependent beta: `beta(t) = beta_max × t` | Preserve exploration in early diffusion, allow convergence in late steps |

The feature was implemented in full and tested on all 47 complexes, compared against three baselines: `baseline_antigen_cut` (vanilla Boltz-2, B1), `baseline_contact_restraints` (oracle contact restraints, B2), and `baseline_pocket` (pocket-based restraints, B3).

---

## 2. How Results Were Produced

All plots and tables use the **best model per complex**, selected by `confidence_score` from the Boltz-2 JSON output files. Each method generates up to 5 structure predictions per complex; the model with the highest confidence score is used for DockQ, CDR3 RMSD, and pLDDT comparisons. Multi-setting baselines (contact_restraints, pocket) additionally select the best setting per complex before comparison. Ensemble diversity uses all available models per complex.

Structural quality (DockQ) was computed using DockQ v2 (`conda activate dockq2`). CDR3 RMSD was computed by aligning predictions onto ground truth using antibody framework Cα atoms (non-CDR residues of chains B and C), then measuring backbone (N/CA/C/O) RMSD of CDR3 residues. CDR indices were taken from `examples/cdrs.csv`.

---

## 3. Figure-by-Figure Analysis

---

### Figure 1: `dockq_comparison_boxplot.png` — Interface Quality

**How produced**: Best model per complex (by confidence_score), DockQ computed against minimized crystal structure.

**Three subplots**: DockQ averaged over antigen-antibody interfaces (A-B + A-C), DockQ for the heavy chain–antigen pair (A-B), and DockQ for the light chain–antigen pair (A-C).

**Subplot 1 — DockQ (Ab-Ag avg)**: The contact_restraints baseline (orange) is clearly the best-performing method, with a box sitting substantially above the others and median near 0.2. The three remaining methods — antigen_cut (blue), pocket (green), and L+ v2 (red) — are closely clustered, all with medians below 0.23 (CAPRI Incorrect threshold, dotted gray line). L+ v2 is marginally above antigen_cut but the difference is small. All distributions are right-skewed with long tails and outliers near 0.8–0.9 corresponding to the few complexes where any method produces a correct pose.

**Subplot 2 — DockQ (A-B, heavy chain)**: The same ordering holds. Contact_restraints has a noticeably higher median and box. L+ v2 is marginally above antigen_cut. Pocket (green) is essentially identical to antigen_cut.

**Subplot 3 — DockQ (A-C, light chain)**: Identical pattern. Contact_restraints leads; the other three are indistinguishable visually. Note that A-C values are systematically very similar to A-B, reflecting the typical antibody-antigen geometry.

**Conclusion**: L+ v2 achieves a modest but non-significant DockQ improvement over vanilla Boltz-2 (+0.020, p=0.47). Contact restraints remain the dominant baseline (statistically significant advantage over v2, p=0.003). The docking interface quality improvement from L+ v2 is real but too small and too noisy to be conclusive with 47 complexes.

---

### Figure 2: `cdr3_rmsd_comparison.png` — CDR3 Loop Accuracy

**How produced**: Best model per complex (by confidence_score), backbone (N/CA/C/O) RMSD after framework alignment.

**Three subplots**: CDR3-H RMSD (Å), CDR3-L RMSD (Å), CDR3 Combined RMSD (Å).

**Subplot 1 — CDR3-H RMSD**: All four methods are visually nearly identical, with interquartile ranges spanning roughly 1.5–6 Å and medians around 2.5–3 Å. The extreme outliers at ~30 Å are present in all methods and represent structurally failed predictions. L+ v2 (red) is slightly lower than antigen_cut (blue), matching the numerical result (3.81 vs 4.00 Å mean), but the overlapping boxes confirm the non-significance.

**Subplot 2 — CDR3-L RMSD**: All methods nearly indistinguishable. L+ v2 and antigen_cut are essentially identical (2.42 vs 2.41 Å mean). Contact_restraints is marginally lower. The L3 distribution is tighter than H3 (reflecting that L3 is typically easier to predict and shorter).

**Subplot 3 — CDR3 Combined RMSD**: Matches the mean pattern: contact_restraints is best (3.37 Å), followed by v2 (3.52 Å), then antigen_cut (3.64 Å), then pocket (3.73 Å). Again, the differences are small relative to the variance.

**Conclusion**: L+ v2 does not worsen CDR3 RMSD — this is the most important finding relative to Round 1 where L v1 significantly degraded RMSD (p=0.0008). The small improvement is not statistically significant (combined: −0.12 Å, p=0.84), but neutralizing the Round 1 regression is a genuine success. The asymmetric beta strategy (L+.2, H3=0.15 > L3=0.05) does not appear to differentially improve H3 over L3 beyond noise.

---

### Figure 3: `cdr3_plddt_comparison.png` — CDR3 Confidence

**How produced**: Best model per complex (by confidence_score), per-residue pLDDT averaged over CDR3-H and CDR3-L residues from Boltz-2 plddt npz files.

**Three subplots**: CDR3-H pLDDT, CDR3-L pLDDT, CDR3 Combined pLDDT.

**Subplot 1 — CDR3-H pLDDT**: All four methods cluster tightly between 0.75 and 0.95. Contact_restraints (orange) is slightly higher (mean 0.850) than the others (0.817–0.828). L+ v2 (red) is essentially identical to antigen_cut (blue): 0.818 vs 0.817.

**Subplot 2 — CDR3-L pLDDT**: Very tight distributions, all between 0.85 and 0.95. L+ v2 (0.875) is nearly identical to antigen_cut (0.877). Pocket and contact_restraints are slightly higher.

**Subplot 3 — CDR3 Combined pLDDT**: Same pattern. L+ v2 (0.838) matches antigen_cut (0.838) exactly at the mean level. Contact_restraints leads (0.865).

**Conclusion**: The CDR3 pLDDT of L+ v2 is indistinguishable from the vanilla baseline. This is a key positive result: in Round 1, the iPTM score was inflated by +0.091 (false confidence), suggesting the model was becoming overconfident in distorted conformations. Here the confidence is unchanged, confirming that L+.1 (beta reduction 0.3→0.12) and L+.4 (time-decay schedule) successfully eliminated the artificial confidence inflation. The model is neither more nor less confident — and since accuracy is slightly improved (non-significantly), this represents better calibration.

---

### Figure 4: `ensemble_diversity.png` — Conformational Diversity

**How produced**: All 5 models per complex per method; pairwise CDR3 Cα RMSD among models after framework alignment. Bar heights show mean over complexes, error bars show standard deviation.

**Three subplots**: CDR3 Combined diversity, CDR3-H diversity, CDR3-L diversity.

**Subplot 1 — CDR3 Combined Diversity**: All four methods produce mean pairwise diversity above the 1.0 Å threshold (red dashed line). Antigen_cut (1.51 Å, blue) has the highest combined diversity. L+ v2 (1.45 Å, red) is slightly below antigen_cut and similar to contact_restraints (1.41 Å) and pocket (1.44 Å). All error bars are large, indicating complex-to-complex variance dominates.

**Subplot 2 — CDR3-H Diversity**: All methods are above 1.0 Å. Antigen_cut (1.76 Å) leads, followed by L+ v2 (1.68 Å), pocket (1.60 Å), contact_restraints (1.62 Å). L+ v2 retains good H3 diversity but is no longer the highest.

**Subplot 3 — CDR3-L Diversity**: All methods are below 1.0 Å (below the threshold line), indicating low L3 diversity across the board. L+ v2 (0.56 Å) is essentially identical to contact_restraints (0.56 Å) and slightly above antigen_cut (0.53 Å). Pocket has the highest L3 diversity (0.68 Å).

**Conclusion**: L+ v2 maintains adequate CDR3 diversity — above the 1.0 Å threshold for combined and H3 — but does not improve on the vanilla baseline. In Round 1, L v1 increased diversity from 1.51→1.66 Å and fraction above 1.0 Å from 45%→60%, which was its primary achievement. L+ v2 with reduced beta (0.12) and time-decay schedule shows 1.45 Å and ~48% fraction, slightly lower than even the baseline. The diversity benefit of Round 1 is not preserved. This is likely because the weaker and time-decaying beta provides insufficient energy to push CDR3 out of its typical distribution by end of diffusion.

---

### Figure 5: `scatter_vs_baseline.png` — Per-Complex Comparison vs Vanilla Baseline

**How produced**: Best model per complex (by confidence_score); x-axis = antigen_cut value, y-axis = L+ v2 value. Points above the diagonal (y=x) indicate improvement with L+ v2.

**Three subplots**: DockQ (Ab-Ag), CDR3 RMSD, iPTM.

**Subplot 1 — DockQ**: The vast majority of points cluster near the origin (0,0), reflecting that most complexes have near-zero DockQ in both methods. A few high-DockQ outliers appear near the diagonal, suggesting that on complexes where Boltz-2 gets the pose right, L+ v2 neither improves nor harms. A small number of points are clearly above the diagonal (L+ v2 wins) and a few are below (L+ v2 loses). There is no systematic improvement — the pattern mirrors the redistribution effect seen in Round 1 with K and L strategies.

**Subplot 2 — CDR3 RMSD**: Points are concentrated along the diagonal with no visible systematic offset. High-RMSD outliers (failed predictions ~20-25 Å) fall on the diagonal, indicating both methods fail on the same complexes. There is no visible improvement for L+ v2.

**Subplot 3 — iPTM**: Points fall tightly along the diagonal y=x line across the full range (0.4–1.0). This is the strongest visual confirmation that iPTM is no longer inflated: L+ v2 and antigen_cut produce essentially identical confidence scores per complex.

**Conclusion**: L+ v2 produces per-complex DockQ improvements on roughly half the complexes and regressions on the other half — consistent with the redistribution pattern observed with beta-scaling strategies. The iPTM scatter confirms that the inflation problem is fully resolved. The CDR3 RMSD scatter reveals that the feature does not reliably guide individual predictions to better loop geometries.

---

### Figure 6: `per_complex_heatmap.png` — Per-Complex ΔDockQ

**How produced**: ΔDockQ = DockQ(L+ v2) − DockQ(baseline) for each of the 47 complexes. Red = baseline better, green = L+ v2 better, white/cream = tie. Best model per complex.

**Against antigen_cut (B1, left column)**: The heatmap shows a mixed pattern — roughly equal areas of green and red with moderate amplitude (|ΔDockQ| mostly ≤ 0.1). A few complexes show large positive changes (dark green, ΔDockQ ≈ +0.2 to +0.3), mostly in the upper section of the heatmap. Several complexes show moderate regression (ΔDockQ ≈ −0.05 to −0.15, orange-red). This confirms the redistribution nature of the effect: the feature helps some complexes substantially while modestly hurting others.

**Against baseline_contact_restraints (middle column)**: Nearly uniformly red — L+ v2 is worse than contact restraints on the vast majority of complexes, consistent with the statistically significant DockQ deficit (p=0.003). Only a handful of complexes show near-zero or slightly positive delta.

**Against baseline_pocket (right column)**: Mostly cream/beige (negligible difference), with a slight tilt toward green in the upper complexes. L+ v2 is marginally better than the pocket baseline on more complexes than not (+0.012, non-significant), particularly on complexes where the pocket heuristic performs poorly.

**Conclusion**: The heatmap reveals a pattern common to CDR3 beta-scaling strategies: the feature changes the diffusion trajectory strongly enough to find a correct pose on some complexes (green) but occasionally steers away from the correct pose on others (red). It does not systematically improve all complexes uniformly. The feature currently lacks a per-complex selection mechanism that could filter out the regressions and keep only the improvements.

---

## 4. Numerical Summary

| Metric | B1 (antigen_cut) | B2 (contact_restraints) | B3 (pocket) | L+ v2 | Δ vs B1 | p-value |
|--------|-----------------|------------------------|-------------|--------|---------|---------|
| DockQ antigen_ab | 0.179 | **0.291** | 0.190 | 0.199 | +0.020 | 0.466 |
| DockQ A-B (heavy) | 0.183 | **0.295** | 0.196 | 0.208 | +0.025 | 0.331 |
| CDR3-H RMSD (Å) | 4.00 | 3.73 | 4.12 | **3.81** | −0.19 | 0.921 |
| CDR3-L RMSD (Å) | 2.41 | **2.28** | 2.52 | 2.42 | +0.01 | 0.871 |
| CDR3 combined RMSD (Å) | 3.64 | **3.37** | 3.73 | 3.52 | −0.12 | 0.838 |
| iPTM | 0.737 | 0.834 | 0.794 | **0.735** | −0.003 | 0.988 |
| CDR3-H pLDDT | 0.817 | **0.850** | 0.828 | 0.818 | +0.001 | — |
| Ensemble diversity (Å) | **1.51** | 1.41 | 1.44 | 1.45 | −0.06 | — |
| CAPRI Incorrect | 36 | 25 | 26 | **34** | −2 | — |
| CAPRI Medium+ | 9 | **15** | 9 | **12** | +3 | — |

Bold = best among the four methods for that metric.

---

## 5. Did the Feature Achieve Its Goals?

### L+.1 — Reduce beta (0.3 → 0.12): SUCCESS
The iPTM inflation is fully resolved. Round 1 L v1 inflated iPTM by +0.091 (p<0.0001); L+ v2 changes iPTM by −0.003 (p=0.99 vs B1), essentially zero. CDR3 RMSD no longer significantly worsens (p=0.84 vs p=0.0008 in Round 1). The beta reduction successfully eliminated the harmful side effects.

### L+.2 — Asymmetric betas (H3=0.15, L3=0.05): PARTIAL
The numerical data shows CDR3-H RMSD decreases by 0.19 Å (from 4.00 to 3.81) while CDR3-L RMSD stays flat (+0.01 Å). This is the expected directional pattern — H3 improves more than L3 — but neither change is statistically significant. The asymmetric betas appear to be a sensible design but the effect is too small to be conclusive.

### L+.3 — CDR3-antigen interface scaling (cdr3_antigen_beta=0.15): INCONCLUSIVE
Interface quality (DockQ antigen_ab) improved by +0.020 non-significantly. The CAPRI distribution gained 3 Medium-quality predictions (34→34 Incorrect, 9→12 Medium+), which represents a practical improvement in the tail of successes. However, iRMSD and lRMSD from DockQ show no significant changes. Whether L+.3 specifically drove the DockQ trend cannot be isolated from the combined effect of all four improvements.

### L+.4 — Time-dependent beta schedule: AMBIGUOUS
The time-decay schedule was designed to preserve exploration early while allowing convergence late. The absence of iPTM inflation and the maintenance of adequate diversity (1.45 Å vs 1.51 Å for B1) suggest the schedule is not harmful. However, diversity did not increase over B1 (was 1.66 Å with L v1 beta=0.3), suggesting the beta decay may reduce the steering effect too aggressively by end of diffusion.

---

## 6. Overall Conclusions

**L+ v2 is a safer but weaker version of L v1.** The feature successfully eliminates the two main failure modes of Round 1 — iPTM inflation and CDR3 RMSD regression — but at the cost of also losing the primary benefit (increased conformational diversity). The result is a feature that is statistically indistinguishable from the vanilla baseline on all key metrics, with a marginal positive trend on DockQ (+0.020) and a small practical gain in CAPRI Medium-quality predictions (+3 complexes).

**The feature is not yet transformative.** It does not significantly outperform the simplest baseline (antigen_cut) on any metric, and it is significantly inferior to oracle contact restraints on interface quality.

**The redistribution pattern persists.** Like L v1 and K in Round 1, L+ v2 shows complex-specific redistribution: it substantially helps some complexes while modestly hurting others. This pattern suggests that beta-scaling changes the diffusion trajectory without reliably steering it toward the correct binding mode.

---

## 7. Suggested Modifications

Given the analysis above, the following modifications are most likely to improve performance:

1. **Stronger antigen interface scaling (L+.3)**: The `cdr3_antigen_beta=0.15` may still be too weak. Try values in the 0.3–0.5 range for the CDR3-antigen term specifically. Since the time-decay schedule reduces beta to near zero at the end, the effective average antigen scaling over the trajectory is lower than 0.15 × the fraction of the diffusion interval where it's active.

2. **Adjust the time-decay schedule (L+.4)**: The linear decay `beta(t) = beta_max × t` may be too aggressive. Consider a plateau-then-decay schedule (e.g., beta stays at beta_max for t ∈ [0.5, 1.0] and decays for t ∈ [0.0, 0.5]), preserving the full steering effect for longer before the convergence phase.

3. **Per-complex selection / confidence re-calibration**: The scatter plot shows that wins exist but are mixed with losses. Combining L+ v2 with a confidence-aware selection mechanism (e.g., only accept L+ v2 predictions when their confidence exceeds the vanilla baseline confidence) could recover the wins while avoiding the regressions.

4. **Combine with strategy A (FK particles)**: Strategy A (Round 1 rank 1) achieved the strongest DockQ improvement (+0.046) using Feynman-Kac importance sampling for global orientation. Combining A's global pose optimization with L+ v2's CDR3 loop steering could produce complementary gains — A handles the bulk docking geometry, L+ v2 refines the loop-antigen contact.

5. **Separate evaluation of L+.1 alone vs combined**: The current evaluation tests all four improvements together. It is unknown whether L+.3 and L+.4 add value or partially cancel the gains from L+.1/L+.2. An ablation study (L+.1 only, L+.1+.2, L+.1+.2+.3, all four) would identify which combination is optimal.

## summary

The feature in one sentence: L+ v2 is a safer but weaker version of L v1 — it fixes the problems but loses the
  benefits.

  What worked:
  - L+.1 fully resolved iPTM inflation (was +0.091 in Round 1, now −0.003, p=0.99). The confidence problem is gone.
  - CDR3 RMSD no longer significantly worsens (p=0.84 vs p=0.0008 in Round 1). The accuracy regression is neutralized.
  - +3 CAPRI Medium-quality predictions (9→12 Medium+ over B1) — a practical gain.

  What didn't work:
  - Diversity dropped back to baseline (1.45 vs 1.66 Å for L v1; B1 = 1.51 Å). The time-decay schedule appears to reduce
   the beta effect too aggressively by end of diffusion.
  - No significant DockQ improvement (+0.020, p=0.47). The antigen interface push (L+.3) trend is positive but too
  small.
  - Redistribution pattern persists — the heatmap shows the feature helps some complexes and hurts others, with no
  reliable per-complex predictor.

  Top suggested fixes: (1) stronger cdr3_antigen_beta (try 0.3–0.5), (2) plateau-then-decay schedule instead of pure
  linear decay, (3) combine with Strategy A (FK particles) for complementary global+local steering, (4) run ablation to
  isolate which of the four sub-improvements adds value.