# Round 2 Final Evaluation: Upgraded Steering Strategies for Antibody-Antigen Structure Prediction

**Date**: 2026-04-02
**Round 2 strategies evaluated**: N1 (A+Y Hybrid), A+ (FK Particles v2), G+ (Progressive v2), L+ (CDR3 Beta v2), Q+ (Epitope Refinement v2), Y+ (Hierarchical v2)
**Round 1 strategies referenced**: A, D, E, G, K, L, O, Q, V, W, Y (11 total)
**Benchmark**: 47 antibody-antigen complexes, ground truth from minimized crystal structures
**Primary metrics**: DockQ Ab-Ag (docking quality), Epitope F1 (binding site prediction), CDR3 RMSD (loop accuracy)

---

## 1. Baselines

All Round 2 strategies are compared against the same baselines used in Round 1:

| Baseline | Description | DockQ (Ab-Ag) | Epitope F1 | Key Feature |
|----------|-------------|---------------|------------|-------------|
| **B0/B1** | Vanilla Boltz-2, no steering | 0.179 | 0.345 | No oracle info, universal |
| **B2 (Contacts)** | Ground-truth contact restraints | 0.291–0.408 | 0.492–0.583 | Oracle info; many configs/complex |
| **B3 (Pocket)** | Pocket-based restraints | 0.148–0.217 | 0.283–0.369 | Moderate oracle info |

**Reminder on fairness**: B2 runs ~11 restraint configurations per complex using ground-truth interface knowledge. Strategies requiring no oracle information should be compared primarily against B1 (the simple baseline). B2 serves as an informative upper bound showing how far oracle structural knowledge can take performance.

---

## 2. Round 2 Strategy-by-Strategy Analysis

### 2.1 N1: A+Y Hybrid (FK Particles + Hierarchical Timing)

**Parent strategies**: A (Round 1 rank #1) + Y (Round 1 rank #2)
**Goal**: Combine A's FK resampling (best DockQ) with Y's hierarchical timing (only significant CDR-H3 improvement). Early phase: beta-scaling for CDR3 orientation. Late phase: FK particles with coordinate potentials.
**Quantitative targets**: DockQ +0.05 to +0.08 vs B1, CDR-H3 RMSD -0.3 to -0.5A, p<0.05

| Metric | B1 | N1 | Delta | p-value |
|--------|-----|-----|-------|---------|
| DockQ (Ab-Ag, conf) | 0.179 | 0.189 | +0.009 | 0.805 |
| DockQ (Ab-Ag, oracle) | 0.248 | 0.253 | +0.005 | 0.618 |
| CAPRI Incorrect (conf) | 77% | 72% | -5 pp | — |
| Confidence rho | 0.380 | 0.416 | +0.036 | — |
| CDR-H3 RMSD | — | — | Not computed | — |
| Easy-case regressions | — | 8/11 | Target was <5 | — |

**Verdict**: **Failed**. The combination of A and Y did not produce additive gains. DockQ improvement (+0.009) is 18% of the minimum target. The hybrid created an optimization tug-of-war: early-phase beta-scaling pushed trajectories in directions suboptimal for late-phase FK refinement. The catastrophic 7ZOZ_HLA regression (0.644→0.054) and 8 out of 11 easy-case regressions reveal that the steering destabilizes already-good predictions. CDR-H3 RMSD — Y's signature metric — was not even computed, leaving the primary Y contribution unevaluated. The one positive finding is improved confidence calibration (rho 0.380→0.416), but this does not translate to better model selection outcomes.

**Did the upgrade improve over its parents?** No. Round 1 A alone achieved +0.046 DockQ (vs N1's +0.009). Round 1 Y alone achieved significant CDR-H3 RMSD -0.30A (not measured in N1). The combination diluted both effects rather than amplifying them.

---

### 2.2 A+ (FK Particles v2)

**Parent strategy**: A (Round 1 rank #1)
**Goal**: Fix A's three limitations: (1) 25% missing complexes via adaptive particle reduction, (2) poor confidence calibration via composite re-ranking using FK energy, (3) non-significance via better coverage and selection.
**Sub-strategies**: A+.1 (adaptive particle reduction), A+.2 (composite re-ranking), A+.3 (separate FK seeds), A+.4 (embedding fallback)

| Metric | B1 | A+ (NF_V2) | Delta | p-value |
|--------|-----|-----------|-------|---------|
| DockQ (Ab-Ag, conf) | 0.170* | 0.172 | +0.003 | 0.447 |
| DockQ (Ab-Ag, oracle) | 0.242* | 0.289 | +0.047 | 0.126 |
| Coverage | 47/47 | 33/48 (69%) | — | — |
| Confidence rho | 0.380 | 0.216 | -0.164 | — |
| FK energy vs DockQ rho | — | 0.022 | — | p=0.78 |
| fnat (conf) | 0.181 | 0.162 | -0.019 | — |

*Paired B1 means on N=33.

**Verdict**: **Mixed — generation improved, selection failed critically**. The most important finding is the stark confidence-vs-oracle gap: under oracle selection, A+ achieves 0.289 DockQ, outperforming B1 (0.242) and B2 (0.280). Under confidence selection, it drops to 0.172 — worse than B1. This means the FK mechanism successfully generates better candidate structures, but confidence-based selection actively destroys this advantage.

The individual sub-strategy outcomes are informative:
- **A+.1 (adaptive particles)**: Partially failed — still missing 31% of complexes (worse than Round 1's 25%).
- **A+.2 (composite re-ranking)**: **Definitively failed** — FK energy has zero correlation with DockQ (rho=0.022, p=0.78). The core hypothesis that steering energy predicts quality is falsified.
- **A+.3 (FK seeds)**: Inconclusive.
- **A+.4 (embedding fallback)**: Inconclusive.

**Did the upgrade improve over its parent?** No. Round 1 A achieved +0.046 DockQ (conf) with zero CAPRI downgrades. A+ achieves +0.003 with 19/33 regressions. The adaptive particle reduction (lowering from 20 to ~5-10 particles) likely reduced ensemble quality. The "safety net" property of Round 1 A (zero downgrades) is lost.

**Key insight**: FK particle resampling works for generation but not for selection. The oracle DockQ of +0.047 vs B1 matches Round 1's +0.046, confirming the mechanism is reproducibly effective. Solving the selection problem — not the generation problem — is the bottleneck.

---

### 2.3 G+ (Progressive Steering v2)

**Parent strategy**: G (Round 1 rank #3)
**Goal**: Fix G's regressions on easy cases via (G+.2) time-varying beta that decays to zero in late steps and (G+.3) progressive antigen orientation potential. G+.1 (conditional steering) was not implemented.

| Metric | B1 | G+ | Delta | p-value |
|--------|-----|-----|-------|---------|
| DockQ (Total, conf) | 0.338 | 0.345 | +0.007 | 0.305 |
| CAPRI Medium+ (%) | 21.3% | 23.4% | +2.1 pp | — |
| Improved/Unchanged/Regressed | — | 4/41/2 | — | — |
| Ensemble diversity | 8.5A | 9.5A | +1.0A | — |

**Notable individual results**: 8OL9_BAH +0.361, 8HGM_CDB +0.325 (large wins); 7TRH_HBG -0.285 (large regression).

**Verdict**: **Marginal — safer but weaker than Round 1**. G+ achieved its primary goal of fewer regressions (2 vs ~3 in Round 1), but at the cost of fewer improvements (4 vs ~6). The time-varying beta is more conservative: 87% of complexes are unchanged. The mean delta +0.007 is smaller than Round 1's +0.017. The 7TRH_HBG regression (-0.285) is actually worse than any single regression in Round 1.

**Did the upgrade improve over its parent?** Partially. Fewer total regressions (the stated goal), but the net effect is smaller. The decaying beta prevents late-step distortion but also prevents late-step improvement. The unimplemented G+.1 (conditional steering) — which would skip steering on high-confidence cases — was likely the most impactful proposed improvement and its absence limits the evaluation.

---

### 2.4 L+ (CDR3 Beta-Scaling v2)

**Parent strategy**: L (Round 1 rank #4)
**Goal**: Fix L's two failures: (1) significant CDR3 RMSD worsening (p=0.0008), (2) dangerous iPTM inflation (+0.091, p<0.0001). Implemented via reduced beta (0.3→0.12), asymmetric H3/L3 betas, CDR3-antigen interface scaling, and time-dependent beta decay.

| Metric | B1 | L+ | Delta | p-value |
|--------|-----|-----|-------|---------|
| DockQ (Ab-Ag, conf) | 0.179 | 0.199 | +0.020 | 0.466 |
| CDR3-H RMSD | 4.00 | 3.81 | -0.19 | 0.921 |
| CDR3-L RMSD | 2.41 | 2.42 | +0.01 | 0.871 |
| iPTM | 0.737 | 0.735 | -0.003 | 0.988 |
| Ensemble diversity | 1.51A | 1.45A | -0.06A | — |
| CAPRI Medium+ | 9 | 12 | +3 | — |

**Verdict**: **Successful safety fix, but lost the primary benefit**. L+ is the clearest example of the safety-vs-efficacy tradeoff in Round 2:

*What was fixed (success)*:
- iPTM inflation **fully resolved**: from +0.091 (p<0.0001) in Round 1 to -0.003 (p=0.988). The confidence problem is gone.
- CDR3 RMSD regression **neutralized**: from significant worsening (p=0.0008) to no change (p=0.921).
- +3 CAPRI Medium-quality predictions — a practical gain.

*What was lost (cost)*:
- Conformational diversity dropped back to baseline (1.45A vs 1.66A for L v1; B1=1.51A). The primary achievement of Round 1 — demonstrating that beta-scaling can increase CDR3 exploration — is no longer present.
- No significant DockQ improvement (+0.020, p=0.47).

**Did the upgrade improve over its parent?** It depends on the framing. L+ is a better tool in production (no harmful side effects), but L v1 was a more informative scientific result (proved beta-scaling works for diversity). L+ is safer but does not advance the science.

---

### 2.5 Q+ (Iterative Epitope Refinement v2)

**Parent strategy**: Q (Round 1 rank #5)
**Goal**: Fix Q's bootstrapping problem (incorrect Round 1 poses misleading later rounds) via more diversity, confidence-weighted contact maps, ensemble-based hotspot selection, and adaptive rounds.

| Metric | B1 | Q+ | Delta | p-value |
|--------|-----|-----|-------|---------|
| DockQ (Ab-Ag, conf) | 0.179 | 0.076 | **-0.103** | **0.035** |
| Global DockQ | — | — | -0.197 | **1.4e-14** |
| Epitope F1 (mean) | 0.352 | 0.281 | -0.071 | 0.132 |
| Epitope MCC | 0.269 | 0.121 | **-0.149** | **0.002** |
| Confidence score | 0.883 | 0.517 | -0.366 | — |
| CAPRI Medium+ | 19.1% | 0% | -19.1 pp | — |

**Verdict**: **Harmful — significantly worse than baseline on multiple metrics**. Q+ is the worst-performing strategy in either round. It is statistically significantly worse than B1 on Global DockQ (p=1.4e-14, **zero wins out of 47**), Ab-Ag DockQ (p=0.035), and epitope MCC (p=0.002). It achieves 0% Medium-quality or better CAPRI predictions. The confidence score collapses to 0.517 (vs B1's 0.883), and only 1 model was produced per complex (eliminating selection opportunities).

The one positive signal — rescuing 10 of 13 hard-zero epitope cases — is real but overwhelmed by severe degradation on everything else.

**Did the upgrade improve over its parent?** No — it is dramatically worse. Round 1 Q achieved +0.029 epitope F1 (mean) and +0.140 median F1 improvement. Q+ drops to -0.071 mean F1 and -0.149 MCC. The v2 improvements exacerbated the bootstrapping problem rather than fixing it: the multi-round pipeline with more aggressive steering disrupted the model's internal representations, producing universally degraded structures.

**Root cause**: The multi-round pipeline fundamentally disrupts the diffusion model. Feeding intermediate predictions back as constraints for subsequent rounds introduces systematic errors that compound rather than self-correct. The model itself recognizes this — confidence collapses to 0.517, the lowest of any method ever tested.

---

### 2.6 Y+ (Hierarchical Steering v2)

**Parent strategy**: Y (Round 1 rank #2)
**Goal**: Strengthen Y's two phases (stronger beta_max=0.8 in early phase, added CDR proximity potential in late phase), add epitope-focused guidance (Y+.3) and per-complex weight scaling (Y+.4).

| Metric | B1 | Y+ | Delta | p-value |
|--------|-----|-----|-------|---------|
| DockQ (Total, conf) | 0.338 | 0.361 | +0.023 | 0.582 |
| DockQ (Ab-Ag, conf) | 0.179 | 0.213 | +0.033 | 0.589 |
| CDR-H3 RMSD | 3.00 | 2.84 | -0.16 | 0.846 |
| CAPRI Medium+High | 19.1% | 27.7% | +8.6 pp | — |
| Epitope F1 (5A) | 0.345 | 0.336 | -0.009 | 0.116 |
| iptm | 0.738 | 0.756 | +0.019 | 0.068 |
| Confidence calibration r | 0.25 | 0.14 | -0.11 | — |

**Notable individual results**: 8FXB_HLE +0.881 (near-complete failure rescued to high quality), 8OL9_BAH +0.545; but 20 losses vs 16 wins overall.

**Verdict**: **Best Round 2 strategy overall, but below expectations**. Y+ produces the strongest positive trends of any Round 2 strategy: +0.033 Ab-Ag DockQ, +8.6 pp CAPRI Medium+High (its best metric, exceeding Round 1's +6.4 pp), and -0.16A CDR-H3 RMSD. The CAPRI improvement means ~4 additional complexes moved from failure to medium/high quality — including the spectacular 8FXB_HLE rescue (+0.881). Y+ matches pocket-guided performance (27.7% vs 27.8% Medium+) despite using no oracle information.

However, none of the improvements reach statistical significance, and the CDR-H3 RMSD improvement (-0.16A) is weaker than Round 1 Y's -0.30A (p=0.025). Confidence calibration degraded (r=0.14 vs 0.25), meaning the model becomes overconfident without matching accuracy, impairing model selection in practice. Y+.3 (epitope guidance) was only testable on 2/47 complexes, leaving a key sub-goal unevaluated.

**Did the upgrade improve over its parent?** Mixed. Y+ improved CAPRI rates (+8.6 pp vs +6.4 pp) but lost statistical significance on CDR-H3 RMSD (-0.16A, p=0.846 vs -0.30A, p=0.025). The stronger beta (0.8 vs original) and added potentials did not amplify DockQ beyond Round 1 levels (+0.033 vs +0.034). The upgrades neither helped nor hurt the core mechanism — the effect is essentially identical to the original Y.

---

## 3. Global Ranking: All Strategies (Round 1 + Round 2)

Strategies ranked by closeness to achieving their individual goals, quality of improvement, statistical evidence, and absence of harm. Round 2 strategies are marked with (R2).

| Rank | Strategy | Round | Category | Key Achievement | DockQ vs B1 | Significant? | Harm? |
|------|----------|-------|----------|----------------|-------------|-------------|-------|
| 1 | **A (20 FK particles)** | R1 | Docking | Best DockQ improvement (+0.046), 87% of B2's fnat, 0 CAPRI downgrades | +0.046 | No (p=0.61) | 25% missing |
| 2 | **Y (Hierarchical)** | R1 | Docking+CDR3 | **Only significant structural improvement** (CDR-H3 -0.30A, p=0.025) | +0.034 | **Yes** | None |
| 3 | **Y+ (Hierarchical v2)** | R2 | Docking+CDR3 | Best CAPRI shift (+8.6 pp Medium+High), dramatic per-complex wins | +0.033 | No | Calibration ↓ |
| 4 | **G (Progressive)** | R1 | Docking | +3 net CAPRI upgrades, helps medium-difficulty cases | +0.017 | No (p=0.41) | Slight regressions |
| 5 | **L+ (CDR3 Beta v2)** | R2 | CDR3 safety | Fixed iPTM inflation and CDR3 RMSD regression; +3 CAPRI Medium | +0.020 | No | None |
| 6 | **L (CDR3 Beta)** | R1 | CDR3 diversity | Proved beta-scaling increases diversity (45%→60%) | +0.022 | No | CDR3 RMSD ↑ (p=0.0008) |
| 7 | **Q (Iterative)** | R1 | Epitope | Best epitope F1 improvement (+0.029 mean, +0.140 median) | +0.019 | No | Minor |
| 8 | **G+ (Progressive v2)** | R2 | Docking | Fewer regressions than G, but weaker effect | +0.007 | No | 7TRH_HBG -0.285 |
| 9 | **A+ (FK Particles v2)** | R2 | Docking | Oracle DockQ confirms generation works (+0.047); FK energy falsified | +0.003 | No | Calibration ↓↓ |
| 10 | **K (Region Beta)** | R1 | Epitope | Marginal epitope gain, strong redistribution pattern | +0.003 | No | Catastrophic on some |
| 11 | **D (Canonical)** | R1 | CDR diversity | Computational efficiency achieved, quality unchanged | +0.009 | No | Calibration ↓ |
| 12 | **N1 (A+Y Hybrid)** | R2 | Docking | Improved calibration (+0.036 rho); but easy-case regressions | +0.009 | No | 8/11 easy regressions |
| 13 | **W (Embed Interface)** | R1 | Docking | Negligible effect, concept unvalidated | +0.004 | No | No |
| 14 | **O (Asymmetric Beta)** | R1 | CDR balance | Zero effect — model already has H3 bias | +0.001 | No | No |
| 15 | **E (Blind Scanning)** | R1 | Epitope | Severe structural degradation (DockQ -51%) | -0.091 | **Yes (worse)** | Severe |
| 16 | **V (Embed CDR3)** | R1 | CDR3 | Every metric degraded significantly | -0.034 | **Yes (worse)** | Severe |
| 17 | **Q+ (Epitope Refine v2)** | R2 | Epitope | **Worst strategy tested**: 0/47 Global DockQ wins, 0% Medium+ | -0.103 | **Yes (worse)** | Severe |

---

## 4. Did Round 2 Upgrades Improve Over Round 1?

This is the central question. The answer is sobering.

| Upgrade | Parent (R1) | R1 DockQ Delta | R2 DockQ Delta | R1→R2 Change | Improved? |
|---------|-------------|---------------|---------------|--------------|-----------|
| N1 (A+Y) | A + Y | +0.046 (A), +0.034 (Y) | +0.009 | Much worse | **No** |
| A+ | A | +0.046 | +0.003 (conf) | Worse | **No** (but oracle matched) |
| G+ | G | +0.017 | +0.007 | Worse | **No** |
| L+ | L | +0.022 | +0.020 | Similar | **Neutral** (fixed harms) |
| Q+ | Q | +0.019 | -0.103 | Much worse | **No** (harmful) |
| Y+ | Y | +0.034 | +0.033 | Similar | **Neutral** (CAPRI better) |

**Summary**: No Round 2 strategy clearly outperformed its Round 1 parent on the primary DockQ metric. Two strategies (L+, Y+) maintained similar performance while fixing identified problems. Two (N1, G+) underperformed their parents. One (A+) revealed a fundamental selection problem. One (Q+) became harmful.

### Why did the upgrades not improve?

1. **Combination is harder than addition**: N1 assumed A's and Y's benefits would stack. Instead, their mechanisms conflicted — early-phase beta-scaling and late-phase FK particles optimize different objectives, creating a tug-of-war.

2. **Conservative fixes reduce both harm and benefit**: L+ and G+ successfully reduced regressions by weakening the steering signal (lower beta, decaying schedules). But weaker steering also means weaker improvement. The safety-efficacy tradeoff was not overcome.

3. **Confidence calibration is a systemic bottleneck**: A+ demonstrated that FK particles generate better structures (oracle +0.047) but confidence-based selection destroys this advantage. Y+ showed confidence inflation (iptm↑ without DockQ↑). Every steered method degrades the confidence-DockQ correlation, and no Round 2 strategy addressed this.

4. **Multi-round pipelines are destructive**: Q+ proved that feeding intermediate results back through the diffusion model compounds errors rather than correcting them. The bootstrapping problem cannot be solved by iterating within the same model.

5. **Hyperparameter interaction was not optimized**: Each parent strategy's parameters were tuned independently. When combined (N1) or modified (A+), the parameter interactions were not re-optimized, leading to suboptimal configurations.

---

## 5. Cross-Strategy Analysis

### 5.1 The Generation-vs-Selection Gap

The most important finding across all Round 2 evaluations is the **generation-selection dissociation**:

| Strategy | Oracle DockQ vs B1 | Confidence DockQ vs B1 | Gap |
|----------|--------------------|------------------------|-----|
| A+ | +0.047 | +0.003 | 0.044 |
| Y+ | +0.015 (total) | +0.023 (total) | -0.008 |
| N1 | +0.005 | +0.009 | -0.004 |
| G+ | 0.000 (=0.391) | +0.007 | -0.007 |

A+ has the largest generation-selection gap (0.044 DockQ lost to bad selection). This means A+ produces better structures than any non-oracle method, but cannot identify them. The FK energy — the proposed selection signal — has zero predictive power (rho=0.022). This is the single most actionable finding: **solving model selection would immediately unlock DockQ gains that already exist in the steered ensembles**.

### 5.2 The Redistribution Pattern

Every beta-scaling strategy (L, L+, K, G, G+, Y, Y+) shows the same qualitative pattern: dramatic wins on a small subset of complexes offset by modest losses on others. This is not noise — it is a systematic consequence of steering the diffusion trajectory:

- Beta-scaling changes which binding mode the model finds
- Some mode changes are beneficial (pushing toward the correct epitope)
- Some are harmful (pushing away from an already-correct orientation)
- There is no per-complex mechanism to determine the direction of the change a priori

Until a per-complex steering selector is developed (e.g., N3's conditional steering, or running steered + unsteered and picking the better one post-hoc), this redistribution pattern will persist.

### 5.3 The ~68% Failure Floor

Across all strategies and both rounds, approximately 60-77% of complexes remain in the CAPRI "Incorrect" category:

| Method | CAPRI Incorrect |
|--------|----------------|
| B1 (vanilla) | 77% |
| Best R1 strategy (A) | ~72% |
| Best R2 strategy (Y+) | 68% |
| B2 (contacts, oracle) | 37% |

No steering strategy has breached the ~65% floor without oracle information. B2 achieves 37% only because it uses ground-truth contact knowledge and runs many configurations. The ~68% failure cluster represents complexes where the antibody-antigen global orientation is fundamentally wrong — the antibody is pointed at the wrong hemisphere of the antigen. Local steering (CDR3 beta-scaling, proximity potentials) cannot fix a global orientation error. This suggests that a coarse rigid-body pre-docking step or much stronger early-phase orientation forcing is needed to make progress on the hardest cases.

### 5.4 Epitope Discovery Remains Unsolved

Round 1's Q was the most promising epitope-discovery approach (median F1 +0.140). Round 2's Q+ regressed catastrophically. The lesson is that the iterative approach has potential but the execution — feeding intermediate structures back through the diffusion model — is destructive.

A more viable path may be to decouple epitope discovery from structure generation: use the multi-round approach only to identify candidate epitope residues, then feed those as soft restraints into a single clean forward pass. This was proposed as modification #5 in Q+'s analysis and echoes the N2 (Q→B2 pipeline) strategy from the original Round 2 proposals.

Comparing all epitope-focused strategies against B1:

| Strategy | Epitope F1 vs B1 | Round | Assessment |
|----------|-----------------|-------|------------|
| Q | +0.029 (mean), +0.140 (median) | R1 | Most promising, not significant |
| K | +0.018 | R1 | Marginal |
| Y+ | -0.009 | R2 | Slightly worse |
| E | -0.005 | R1 | No improvement + structural harm |
| Q+ | -0.071 | R2 | Significantly worse |

Q v1 remains the best epitope-discovery strategy tested.

### 5.5 Confidence Calibration Degradation is Universal

Every steered method tested in both rounds degrades confidence calibration:

| Method | Confidence-DockQ correlation | Change vs B1 |
|--------|------------------------------|---------------|
| B1 (vanilla) | rho=0.380 / r=0.25 | reference |
| A+ | rho=0.216 | -0.164 |
| Y+ | r=0.14 | -0.11 |
| N1 | rho=0.416 | +0.036 (exception) |
| Round 1 D | rho=0.177 | -0.110 |
| Round 1 A | rho=0.265 | -0.115 |

N1 is the one exception (improved calibration), but its overall DockQ was too poor for this to matter. The pattern is clear: steering the diffusion process disrupts the model's ability to assess its own prediction quality. This is likely because the confidence head was trained on unsteered trajectories; steered structures produce activation patterns that the confidence head has not learned to evaluate.

---

## 6. Comparison Against B1 vs B2: Where Each Strategy Stands

For strategies that do not use oracle information, the primary comparison should be against B1. The comparison against B2 shows the remaining gap.

### Against B1 (Simple Baseline — Fair Comparison)

| Strategy | DockQ vs B1 | Key Non-DockQ Win | Worth Pursuing? |
|----------|-------------|-------------------|-----------------|
| Y+ (R2) | +0.033 | +8.6 pp CAPRI Medium+High | Yes — best non-oracle CAPRI improvement |
| Y (R1) | +0.034 | CDR-H3 -0.30A (p=0.025) | Yes — only significant structural gain |
| L+ (R2) | +0.020 | Fixed iPTM inflation; +3 CAPRI Medium | Yes — safe, no harm |
| A (R1) | +0.046 | 0 CAPRI downgrades | Yes — best DockQ, needs selection fix |
| Q (R1) | +0.019 | Median epitope F1 +0.140 | Yes — best epitope discovery |
| G+ (R2) | +0.007 | Diversity preserved | Weak — needs conditional steering |
| A+ (R2) | +0.003 | Oracle +0.047 (generation works) | Yes — selection fix unlocks gains |

### Against B2 (Contact Restraints — Gap Analysis)

| Metric | Best non-oracle strategy | B2 (contacts) | Gap |
|--------|------------------------|---------------|-----|
| DockQ (conf) | A: 0.247 / Y+: 0.213 | 0.291–0.408 | 0.04–0.20 |
| CAPRI Medium+High | Y+: 27.7% | 43.5% (best) | 16 pp |
| Epitope F1 | Q: 0.381 | 0.583 (best) | 0.20 |
| fnat | A: 0.252 | 0.289 | 0.04 |

The gap to B2 remains substantial. The most promising path to closing it is not through stronger steering of a single method, but through: (1) solving model selection to recover A+'s oracle performance, (2) combining complementary mechanisms (A's orientation + Y's CDR3 + Q's epitope discovery), and (3) running multiple configurations per complex (mimicking B2's multi-config advantage).

---

## 7. Lessons Learned from Round 2

### 7.1 What Round 2 Confirmed

1. **FK particle resampling is the most effective generation mechanism**. Both A (Round 1) and A+ (Round 2) produce better oracle DockQ than any other non-oracle method. The oracle improvement (+0.046-0.047) is reproducible across rounds.

2. **Hierarchical timing (Y) produces the most consistent positive trends**. Y and Y+ both show +0.033-0.034 Ab-Ag DockQ, the most consistent non-oracle improvements tested. The CAPRI rate improves in both rounds.

3. **Beta-scaling is mechanistically effective but lacks selectivity**. L (R1) proved it increases diversity. L+ (R2) proved it can be made safe. But no beta-scaling strategy reliably steers toward the correct binding mode.

4. **The ~68% failure floor is a global orientation problem**, not a local refinement problem. No combination of CDR3 scaling, proximity potentials, or FK particles breaks through this floor.

### 7.2 What Round 2 Revealed (New Insights)

1. **Combining strategies is non-trivial**. N1 proved that naively stacking two successful mechanisms can produce worse results than either alone. Mechanism interactions must be studied and parameters re-optimized for the combined system.

2. **FK steering energy is not a selection signal**. A+.2 definitively falsified the hypothesis that FK potential energy predicts docking quality (rho=0.022). Future selection mechanisms must use different signals.

3. **Multi-round pipelines damage the diffusion model**. Q+ proved that feeding intermediate structures back through the model compounds errors. Iterative refinement must be decoupled from the diffusion process.

4. **Conservative parameter fixes preserve safety but lose efficacy**. L+ and G+ demonstrate that weakening steering to eliminate harm also eliminates benefit. The solution is not weaker universal steering but conditional steering — strong where needed, absent where not.

5. **Confidence calibration is the universal bottleneck**. If one problem could be solved to immediately improve all strategies, it is confidence-based model selection. A+ loses 0.044 DockQ to bad selection. Y+ loses an unknown amount. A custom re-ranker trained on steered trajectory features would benefit every strategy.

---

## 8. Strategies That Can Be Presented as Wins

Based on all evidence from both rounds, the following strategies represent genuine contributions:

### Tier 1: Clear Wins

**1. Strategy Y (Round 1) — Hierarchical Steering**
- **Highlight**: The only strategy across 17 tested with a statistically significant structural improvement (CDR-H3 RMSD -0.30A, p=0.025).
- **Why it matters**: Demonstrates that the two-phase approach (embedding early, coordinate late) can improve the most clinically relevant loop without degrading anything else. The architectural insight — that different diffusion phases need different steering mechanisms — is fundamental.
- **Best supporting evidence**: No degradation on any metric. Positive trends on DockQ (+0.034) and CAPRI (+6.4 pp).

**2. Strategy A (Round 1) — FK Particle Resampling**
- **Highlight**: Best DockQ improvement among all strategies (+0.046 Ab-Ag), confirmed reproducible in Round 2 (A+ oracle: +0.047). Zero CAPRI downgrades (safety net property in Round 1).
- **Why it matters**: Proves that importance sampling with multiple particles is the most effective mechanism for improving antibody-antigen docking without oracle information. The generation capability is real and consistent.
- **Caveat to acknowledge**: Not significant (p=0.61), 25% coverage gap, confidence calibration degrades.

**3. Strategy L+ (Round 2) — Safe CDR3 Beta-Scaling**
- **Highlight**: Fully resolved L v1's harmful side effects (iPTM inflation eliminated: from +0.091 to -0.003; CDR3 RMSD regression neutralized: from p=0.0008 to p=0.921) while maintaining a positive DockQ trend (+0.020) and gaining +3 CAPRI Medium predictions.
- **Why it matters**: Demonstrates that beta-scaling can be made production-safe with proper parameter tuning. Provides a template for how to iterate on a mechanism: fix the harms first, then strengthen the signal.

### Tier 2: Promising Directions

**4. Strategy Y+ (Round 2) — Enhanced Hierarchical**
- **Highlight**: Best CAPRI improvement of any strategy (+8.6 pp Medium+High) including the spectacular 8FXB_HLE rescue (+0.881 DockQ). Matches pocket-guided restraint performance without oracle information.
- **Why it matters**: Shows the hierarchical approach scales — stronger parameters produce larger CAPRI gains. The dramatic per-complex wins prove the mechanism can find correct binding modes on complexes all other methods fail on.
- **Caveat**: Lost CDR-H3 significance; confidence degraded.

**5. Strategy Q (Round 1) — Iterative Epitope Refinement**
- **Highlight**: Best blind epitope discovery method tested (median F1 +0.140, +60% relative improvement). Rescued multiple hard-zero cases.
- **Why it matters**: The only strategy that meaningfully improves epitope prediction on the hardest cases where vanilla Boltz-2 finds nothing. The concept of iterative narrowing is sound — the execution needs decoupling from the diffusion process (as Q+ showed).

**6. A+ Oracle Performance (Round 2) — FK Generation Capability**
- **Highlight**: FK particle ensembles consistently outperform all non-oracle methods under oracle selection (DockQ 0.289 vs B1's 0.242, B2's 0.280).
- **Why it matters**: The generation problem is solved. Better-than-baseline structures exist in the ensemble. Solving selection unlocks +0.047 DockQ immediately. This is the highest-ROI research direction.

### Tier 3: Informative Scientific Results

**7. L v1 (Round 1)**: Proved beta-scaling mechanistically controls CDR3 diversity. This is confirmed science even though the production version (L+) trades the diversity for safety.

**8. A+.2 Falsification**: Definitively proved FK steering energy does not predict DockQ (rho=0.022). This eliminates a hypothesis and redirects selection research toward structural features (buried surface area, interface contacts, Rosetta energy).

---

## 9. Recommended Next Steps

### Immediate (Highest ROI)

1. **Build a custom re-ranking module** (from N6/N12 proposals). Train a lightweight model on Round 1+2 data to predict DockQ from: interface pLDDT, CDR3-antigen contact count, buried surface area, SASA ratio, H3/L3 contact ratio. Cross-validate on 47 complexes. This addresses the universal confidence calibration bottleneck and would immediately improve A+, Y+, and all future strategies.

2. **Implement conditional steering** (N3/G+.1). Run a fast baseline prediction first; skip steering when confidence is high (iPTM > 0.85). This eliminates the easy-case regression problem seen in N1 (8/11 regressions), G+ (7TRH_HBG -0.285), and others. Low implementation cost, directly addresses the redistribution pattern.

3. **Decouple Q's epitope discovery from structure generation**. Use Q v1's iterative approach only to identify candidate epitope residues, then feed those as soft pocket restraints into a clean single-round prediction (the N2 Q→B2 pipeline). This preserves Q's epitope discovery benefit while avoiding Q+'s structural degradation.

### Medium-Term (Combining Complementary Mechanisms)

4. **Ensemble selection across strategies** (N15 concept). Run Y+, A (with fixed particle count), and B1 independently for each complex. Pool all predictions and select the best using the re-ranker from step 1. Different strategies help different complexes — the union of wins should exceed any individual strategy.

5. **Strengthen Y+'s late-phase potentials selectively**. Y+ shows the most consistent trends. The late-phase CDR-antigen proximity potential (weight 0.8) can likely be increased for hard cases (where CDR-antigen distance is large). Combined with conditional steering (skip easy cases), this could amplify wins without regressions.

6. **Increase L+'s CDR3-antigen beta**. L+.3's interface scaling (beta=0.15) was likely too weak. Testing 0.3-0.5 with a plateau-then-decay schedule (maintain full beta for 50% of diffusion, then decay) could recover diversity without the iPTM inflation.

### Research Directions

7. **Address the ~68% failure floor**. The fundamental limit is global orientation. Options include: (a) coarse rigid-body pre-docking (ZDOCK-style FFT scan) to initialize Boltz-2 from 3-5 candidate orientations, (b) much stronger early-phase orientation potentials that force the antibody binding face toward the antigen, (c) antigen surface region elimination strategies that progressively exclude wrong hemispheres.

8. **Investigate why combining A+Y failed**. Run A and Y independently on the exact same test set with the exact same parameters as Round 1 to confirm individual effects are reproducible. If they are, test combinations with a grid search over transition fraction (0.2, 0.3, 0.5, 0.7) and particle count (5, 10, 15).

---

## 10. Final Summary

**Round 2 was a testing round, not a breakthrough round.** No strategy achieved statistical significance on DockQ. No upgrade clearly outperformed its Round 1 parent. But the round produced critical insights:

- **The generation problem is solved** (A/A+ oracle performance reproducibly exceeds baseline)
- **The selection problem is the bottleneck** (confidence calibration degrades universally; FK energy is not a signal)
- **Conservative parameter fixes work** (L+ eliminated harm) **but don't advance efficacy**
- **Naive combination fails** (N1); **iterative pipelines are destructive** (Q+)
- **Y/Y+ is the most reliable mechanism** (consistent positive trends across both rounds)

The path forward is not to develop more steering mechanisms in isolation, but to:
1. **Fix selection** (re-ranker trained on structural features)
2. **Apply steering conditionally** (only when the baseline fails)
3. **Combine strategies at the ensemble level** (pool predictions, select best) rather than within the diffusion loop
4. **Attack the global orientation problem** (the ~68% failure floor that no local steering can breach)

The strategies that should be highlighted as wins are **Y (significant CDR-H3)**, **A (best generation)**, **L+ (safe CDR3 scaling)**, and **Q v1 (best epitope discovery)**. These four address different aspects of the problem and their lessons can be combined — not by running them simultaneously in a single diffusion loop, but by running them independently and selecting the best result per complex.

---

## 11. Combining Strategies with Constraint Baselines

Throughout Rounds 1 and 2, the steering strategies and the constraint baselines (B2 contacts, B3 pocket) were treated as alternatives — each was compared against the vanilla B1 baseline, but they were never tested *together*. This is a significant gap: the constraint baselines and the steering strategies operate through different mechanisms and address different bottlenecks. The right question is not just "does strategy X beat B1?" but also "can strategy X improve B2 or B3 beyond what they achieve alone?"

### 11.1 Why Combination Is Promising

The constraint baselines and the steering strategies modify the diffusion process at different levels:

| Component | Constraint Baselines (B2/B3) | Steering Strategies |
|-----------|------------------------------|---------------------|
| **What they provide** | Explicit residue-pair distance targets (coordinate space) | Attention modulation (embedding space) and/or trajectory-level particle resampling |
| **What they fix** | Global antibody-antigen orientation | CDR3 loop conformations, conformational diversity, trajectory exploration |
| **Where they operate** | Potential energy terms in the diffusion loop (gradient guidance) | Pair representation scaling in DiffusionConditioning (beta-scaling) and/or FK resampling weights |
| **Limitation** | Does not specifically optimize CDR3 loop geometry; fixed restraint configurations cannot adapt | Cannot reliably fix global orientation (~68% failure floor) |

The key insight is that **B2's contact restraints solve the global orientation problem** (reducing Incorrect from 77% to 37%) **but do not specifically target CDR3 loop accuracy**. Conversely, the best steering strategies specifically improve CDR3 conformations (Y: CDR-H3 -0.30A) or explore more trajectory space (A: FK particles), but cannot fix global orientation. These are complementary bottlenecks.

### 11.2 Strategy-by-Strategy Compatibility Assessment

#### Y / Y+ (Hierarchical Steering) + B2 Contacts — **HIGH PROMISE**

**Complementarity**: Y operates in two phases — (1) early-phase beta-scaling on CDR3 pair representations in embedding space, and (2) late-phase coordinate potentials pulling CDR3 toward the antigen. B2 operates entirely in coordinate space via residue-pair distance restraints. The early phase of Y is fully orthogonal to B2: it modifies attention patterns, not coordinates. The late phase partially overlaps with B2 (both pull CDR toward antigen), but Y's CDR-proximity potential targets the CDR3 loops specifically, while B2's restraints target whichever residue pairs were specified — which may not include CDR3-specific contacts.

**Evidence for complementarity**: In the Y+ evaluation, contact restraints achieve CDR-H3 RMSD of 2.35A while Y+ achieves 2.84A (vs B1's 3.00A). These improvements address partially different aspects: B2 improves CDR-H3 primarily by getting the global orientation right (correct orientation = CDR3 near its target), while Y improves CDR-H3 by directly modulating CDR3 attention patterns. Combining them could push CDR-H3 below 2.35A by adding Y's embedding-space CDR3 refinement on top of B2's correctly oriented structure.

**Risk of conflict**: Y's late-phase coordinate potentials could interfere with B2's restraint potentials if both are active on the same residue pairs. Mitigation: disable Y's late-phase potentials and keep only the early-phase beta-scaling when running with B2. This preserves the orthogonal embedding-space benefit while avoiding coordinate-space conflicts.

**Expected benefit**: CDR-H3 RMSD improvement beyond B2 alone; potentially small DockQ gain from better loop conformations in an already correctly-oriented complex.

**Implementation complexity**: Low — simply enable `cdr3_beta_scaling` on the same YAML files that specify contact restraints. Only the early phase needs to be active.

---

#### L+ (CDR3 Beta-Scaling v2) + B2 Contacts — **MEDIUM-HIGH PROMISE**

**Complementarity**: L+ operates purely in embedding space (beta-scaling on CDR3-CDR3 and CDR3-antigen pair representations). B2 operates purely in coordinate space. There is zero mechanism overlap. L+ was shown to be safe (no iPTM inflation, no CDR3 RMSD regression, no harmful side effects). When combined with B2, L+'s CDR3-antigen interface scaling (beta=0.15) could strengthen the attention between CDR3 residues and the antigen residues that B2's restraints are pulling together. This is a "same direction, different lever" combination — B2 pulls the atoms together, L+ amplifies the model's attention to the same interface.

**Evidence**: L+ alone achieved +3 CAPRI Medium predictions over B1 and a positive (non-significant) DockQ trend. B2 alone achieves 43.5% Medium+High. The question is whether L+ can push the B2 predictions that are borderline Acceptable (DockQ 0.23-0.49) into Medium territory by refining the CDR3 loops within B2's correctly-oriented framework.

**Risk**: Minimal. L+ was specifically designed (and verified) to be harmless. The worst case is no additional benefit.

**Expected benefit**: Modest improvement on B2's already-strong performance through better CDR3 loop conformations. Most impactful on complexes where B2 gets the orientation right but the CDR3 loops are imprecise.

**Implementation complexity**: Low — add `cdr3_beta_scaling` and `cdr3_antigen_beta` to existing B2 YAML configurations.

---

#### A / A+ (FK Particles) + B2 Contacts — **MEDIUM PROMISE, HIGH COST**

**Complementarity**: FK particles provide trajectory-level diversity through importance sampling — maintaining multiple structural hypotheses and discarding poor ones. B2 provides structural guidance via contact restraints. Combining them means running multiple particles, each guided by the same contact restraints, and then selecting the best particle. This could help B2 in cases where the restraints are correct but the model finds a local minimum instead of the global one. Multiple particles would explore different ways to satisfy the same restraints.

**Evidence**: A+'s oracle DockQ (0.289) already exceeds B2's confidence-selected DockQ (0.291) on a smaller set. If FK particles ran on top of B2's restraints, the oracle improvement could stack: B2 guides orientation correctly, FK particles explore multiple ways to satisfy that orientation, and the best particle is selected.

**Risk**: Severe memory constraints. B2 already runs ~11 restraint configurations per complex. Adding even 5 FK particles per configuration means 55 forward passes per complex. The OOM problem (already affecting 31% of complexes with A+ alone) would be much worse. Additionally, A+.2 showed that FK energy is uninformative for selection (rho=0.022) — so even if particles improve the oracle, selecting the best particle remains unsolved.

**Expected benefit**: Potentially large oracle DockQ improvement. Practical improvement depends entirely on solving the selection problem first.

**Implementation complexity**: Medium — FK particles already coexist with potentials in the diffusion loop. The main barrier is memory/compute, not code.

**Recommendation**: Only viable after (1) a custom re-ranker is built (Section 9.1) and (2) memory-efficient FK is implemented (gradient checkpointing, reduced particles). Consider running with 3-5 particles on only the B2 configurations that produce borderline results (DockQ 0.1-0.3), where additional exploration is most valuable.

---

#### G+ (Progressive Steering v2) + B2 Contacts — **LOW-MEDIUM PROMISE**

**Complementarity**: G+'s progressive antigen orientation potential (weak early → strong late) substantially overlaps with what B2's restraints already provide. Both guide the antibody toward the antigen. The only non-redundant component is G+'s time-varying CDR3 beta-scaling (exploration early, convergence late), which operates in embedding space.

**Evidence**: G+ showed only +0.007 DockQ over B1 alone. Its progressive potentials are weaker versions of what B2 achieves directly. The time-varying beta component has the same character as Y's early-phase — but G+ bundles it with redundant coordinate potentials.

**Risk**: The coordinate-space potentials in G+ could conflict with B2's restraints (two sets of forces pulling atoms in potentially different directions). Disabling G+'s potentials and keeping only the time-varying beta would make this equivalent to a variant of L+ with a decay schedule.

**Expected benefit**: Marginal at best. G+'s unique contribution (time-varying beta) is already captured by Y and L+ in a cleaner form.

**Recommendation**: Not worth pursuing as a distinct combination. The time-varying beta concept is better implemented through Y or L+ combined with B2.

---

#### N1 (A+Y Hybrid) + B2 Contacts — **LOW PROMISE**

**Complementarity**: N1 already combines two mechanisms (FK particles + hierarchical timing) that proved to conflict with each other. Adding B2's restraints would introduce a third force in coordinate space, further complicating the optimization landscape.

**Evidence**: N1 alone failed to achieve its goals. Its easy-case regressions (8/11 easy cases regressed) suggest the mechanism is unstable. Adding more forces would likely increase instability.

**Risk**: High. Three competing optimization objectives (FK potential energy, CDR3 beta-scaling, contact restraints) in a single diffusion loop.

**Recommendation**: Do not pursue. The individual components (Y + B2, or A + B2) are better tested in isolation.

---

#### Q v1 (Iterative Epitope Refinement) → B2 Pipeline — **HIGH PROMISE (as pipeline, not co-inference)**

**Complementarity**: This is not a case of running Q and B2 simultaneously during inference. Instead, Q's output (predicted epitope residues) would be used to *generate* B2-style contact restraints for complexes where no oracle contact information is available. This is the N2 (Q→B2 pipeline) concept from the Round 2 proposals.

**Evidence**: Q v1 improved median epitope F1 by +0.140 (+60% relative). If those predicted epitope residues are accurate enough, they could serve as contact restraints for a second clean B2-style run. This would bridge the gap between blind prediction and oracle-guided prediction. The question is precision: B2's strength comes from correct contacts. If Q's predictions have >50% precision (i.e., more than half the predicted epitope residues are real), the B2-style restraints would be net beneficial.

**Risk**: Q's epitope predictions may not be precise enough. False-positive epitope residues would create misleading restraints. However, this risk is mitigated by running the B2-style restraints as additional configurations alongside B1 (vanilla), so the worst case is adding unhelpful configurations that get filtered out by model selection.

**Expected benefit**: Potentially large — this could bring B2-level performance (DockQ 0.29-0.41) to complexes without oracle information. Even partial epitope accuracy (F1 ~0.38) could provide a useful subset of correct contacts.

**Implementation complexity**: Medium — requires a pipeline script (run Q v1 → extract top predicted residues → generate restraint YAML → run B2-style inference). No changes to the model code.

**Critical note**: Q+ (the v2 upgrade) should NOT be used for this pipeline. Q+ was harmful and degraded epitope prediction. Use Q v1's approach (which was the best epitope discovery method tested).

---

#### Q+ (Iterative Epitope Refinement v2) + B2 — **NOT RECOMMENDED**

Q+ was statistically significantly worse than B1 on multiple metrics (Global DockQ: p=1.4e-14, epitope MCC: p=0.002). It is harmful. Combining a harmful strategy with B2 would only risk degrading B2's performance. The multi-round pipeline disrupts the model's representations — this disruption would occur even with B2's restraints active.

---

#### K (Region Beta-Scaling, R1) + B3 Pocket — **MEDIUM PROMISE**

**Complementarity**: K emphasizes specific antigen surface regions through beta-scaling, producing a redistribution of binding modes. B3 provides pocket-based restraints that guide the antibody toward a known binding region but with less specificity than B2. K's region emphasis could sharpen B3's broader pocket guidance by amplifying attention to specific sub-regions within the pocket.

**Evidence**: K showed strong per-complex redistribution in Round 1 (some complexes improved by +0.5 DockQ, others regressed by -0.6 DockQ). When combined with B3, the pocket restraints would anchor the global orientation, and K's region emphasis would modulate which specific surface patch within the pocket the model explores. This mitigates K's main failure mode (steering toward the wrong hemisphere) because B3 already constrains the hemisphere.

**Risk**: K's redistribution pattern could still cause regressions even within B3's correctly-oriented framework. But the magnitude of regressions would likely be smaller (B3 prevents catastrophic orientation errors).

**Expected benefit**: Most promising for B3 (pocket) rather than B2 (contacts). Could help B3 achieve more of B2's specificity by narrowing the binding region within the pocket.

---

### 11.3 Ranking: Strategies Most Promising for Combination with Constraint Baselines

| Rank | Strategy | Combine with | Mechanism Orthogonality | Risk | Expected Gain | Priority |
|------|----------|-------------|------------------------|------|---------------|----------|
| 1 | **Y / Y+ (early phase only)** | B2 | **High** — embedding-space CDR3 steering vs coordinate-space restraints | Low (disable late-phase potentials) | CDR-H3 RMSD below B2's 2.35A | **Immediate** |
| 2 | **Q v1 → B2 pipeline** | B2 (as pipeline) | **Full** — Q discovers epitopes, B2 uses them as restraints | Medium (depends on Q's precision) | B2-level DockQ without oracle info | **Immediate** |
| 3 | **L+** | B2 | **High** — embedding-space interface scaling vs coordinate-space restraints | Very low (proven safe) | Modest CDR3 refinement on top of B2 | **Immediate** |
| 4 | **K** | B3 | **Medium** — region emphasis complements pocket guidance | Medium (redistribution) | Sharpen B3's broad pocket to B2-level specificity | **Medium-term** |
| 5 | **A / A+ (FK particles)** | B2 | **Medium** — trajectory diversity vs structural guidance | High (memory, selection) | Large oracle gain; needs re-ranker first | **After re-ranker** |
| 6 | **G+** | B2 | **Low** — progressive potentials redundant with restraints | Low | Marginal (beta component = L+) | **Not recommended** |
| 7 | **N1** | B2 | **Low** — too many conflicting mechanisms | High | Unlikely positive | **Not recommended** |
| 8 | **Q+** | B2 | **None** — harmful strategy | Very high | Negative | **Do not use** |

### 11.4 Recommended Combination Experiments

Based on the analysis above, three combination experiments should be prioritized:

**Experiment 1: Y (early-phase only) + B2 Contact Restraints**
- Run B2's existing YAML configurations with `cdr3_beta_scaling` enabled (beta_max=0.3-0.5 for CDR3, time-decaying schedule).
- Disable Y's late-phase coordinate potentials (they would conflict with B2's restraints).
- Primary metric: CDR-H3 RMSD (target: below B2's 2.35A).
- Secondary metric: DockQ (target: above B2's 0.291 confidence-selected).
- This directly tests whether embedding-space CDR3 refinement stacks on top of coordinate-space contact guidance.

**Experiment 2: L+ + B2 Contact Restraints**
- Run B2's existing YAML configurations with L+'s parameters added (`cdr3_beta=0.12`, `cdr3_antigen_beta=0.15`, `h3_beta=0.15`, `l3_beta=0.05`, time-decay schedule).
- Primary metric: DockQ, CAPRI distribution.
- This is the safest experiment (L+ has no harmful side effects) and tests whether the CDR3-antigen attention amplification strengthens B2's already-good restraints.

**Experiment 3: Q v1 epitope predictions → B2-style restraints (pipeline)**
- Run Q v1 on all 47 complexes to predict epitope residues.
- Convert top predicted residues (F1-thresholded) into B2-format contact restraint YAMLs.
- Run standard Boltz-2 with these generated restraints.
- Primary metric: DockQ vs B1 (target: approach B2's 0.291 without oracle info).
- This tests whether Q's blind epitope discovery is accurate enough to bootstrap B2-level performance — the most impactful potential outcome of any combination.

### 11.5 What About Combining with B3 (Pocket)?

B3 (pocket restraints) is a weaker baseline than B2, achieving DockQ 0.148-0.217 and sometimes even underperforming B1. Strategies that could improve B3 are those that sharpen its broad binding region guidance:

- **K + B3**: K's region emphasis could narrow the binding site within B3's pocket. This addresses B3's main weakness (too broad, not specific enough).
- **L+ + B3**: Same logic as L+ + B2 but with a lower starting point. The CDR3-antigen interface scaling could help B3's pocket-guided structures form better CDR3 contacts.
- **Y early-phase + B3**: CDR3 beta-scaling to improve loop conformations within B3's oriented framework.

B3 combinations are lower priority than B2 combinations because B3's baseline performance is weaker and the potential gains are smaller in absolute terms.

### 11.6 Key Principle: Embedding-Space Strategies Are the Best Candidates

The clearest pattern from this analysis is that **strategies operating in embedding space (beta-scaling: Y, L+, K) are the best candidates for combination with constraint baselines**, which operate in coordinate space. The two mechanism types modify different parts of the diffusion pipeline — beta-scaling modifies pair representation attention biases in `DiffusionConditioning`, while restraint potentials modify atom coordinates via gradient guidance in the diffusion loop. They do not directly interfere with each other.

Conversely, strategies with their own coordinate-space potentials (G+'s progressive potential, Y's late-phase potential, N1's FK potentials + hierarchy) are more likely to conflict with B2/B3's restraint potentials. The recommendation is to strip any coordinate-space components from strategies before combining them with B2/B3, and rely solely on the embedding-space contributions.
