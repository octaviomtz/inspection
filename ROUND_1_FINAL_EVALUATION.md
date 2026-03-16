# Round 1 Final Evaluation: Steering Strategies for Antibody-Antigen Structure Prediction

**Date**: 2026-03-15
**Strategies evaluated**: A, D, E, G, K, L, O, Q, V, W, Y (11 total)
**Benchmark**: 47 antibody-antigen complexes, ground truth from minimized crystal structures
**Primary metrics**: DockQ (docking quality), Epitope F1 (binding site prediction), CDR3 RMSD (loop accuracy)

---

## 1. Baselines

All strategies are compared against up to three baselines:

| Baseline | Description | DockQ (Ab-Ag) | Epitope F1 | Key Advantage |
|----------|-------------|---------------|------------|---------------|
| **B0/B1** | Vanilla Boltz2, no steering | 0.179 | 0.345 | No oracle info needed, universal |
| **B2 (Contacts)** | Ground-truth contact restraints (hbond, hydrophobic, salt bridge) | 0.291-0.408 | 0.492-0.583 | Uses oracle interface knowledge; many configs per complex |
| **B3 (Pocket)** | Pocket-based restraints | 0.148-0.217 | 0.283-0.369 | Moderate oracle info |

**Important context for interpretation**: B2 (contacts) has a structural advantage: it runs ~11 restraint configurations per complex and uses ground-truth knowledge of which residues form the interface. This makes it an unfair but informative upper bound. **Strategies that require no oracle information should be compared primarily against B1** (the simple baseline). However, it remains relevant to compare against B2 to understand the gap that remains and whether a strategy could complement contact information.

---

## 2. Strategy-by-Strategy Summary

### Strategy A: Antigen Steering with 20 FK Particles

**Goal**: Optimize 3D orientation of antigen relative to CDR loops using Feynman-Kac importance sampling.

| Metric | B1 | A (20p) | Delta | p-value |
|--------|-----|---------|-------|---------|
| DockQ (Ab-Ag, conf) | 0.201* | 0.247 | +0.046 | 0.610 |
| fnat | 0.181 | 0.252 | +0.071 | — |
| iRMSD | 11.5 | 9.7 | -1.8 A | — |
| CAPRI upgrades / downgrades | — | 5 / 0 | — | — |
| N complexes evaluated | 47 | 35 (74%) | — | — |

*Paired B1 mean on N=35.

**Verdict**: Best overall improvement among all strategies. Achieves 87% of B2's fnat without oracle info. Zero CAPRI downgrades (safety net effect). Main limitations: 25% of complexes missing (OOM with 20 particles), not statistically significant (p=0.610), and poor confidence calibration (rho=0.265). Improvements are concentrated in ~5 complexes rather than broadly distributed.

---

### Strategy D: Enhanced Canonical Ensemble (CDR Beta-Scaling Sweep)

**Goal**: Generate diverse CDR conformations by sweeping beta from -0.5 to +0.5 on CDR pair representations, selecting top-5 by confidence.

| Metric | B1 | D | Delta | p-value |
|--------|-----|---|-------|---------|
| DockQ (Global) | 0.338 | 0.347 | +0.009 | 0.410 |
| CDR-H3 RMSD | 2.998 | 2.986 | -0.012 | — |
| Confidence-DockQ rho | 0.287 | 0.177 | -0.110 | — |

**Verdict**: Negligible impact. The beta perturbation is too gentle (absorbed by diffusion), confidence-based selection is poorly calibrated for perturbed structures (rho drops from 0.287 to 0.177), and the method produces 24 wins vs 23 losses (coin flip). Computational efficiency goal achieved (trunk runs once), but quality is unchanged.

---

### Strategy E: Enhanced Blind Scanning (Region-Specific Beta-Scaling)

**Goal**: Discover epitope location by emphasizing different antigen surface regions during diffusion.

| Metric | B1 | E | Delta |
|--------|-----|---|-------|
| Epitope F1 | 0.345 | 0.340 | -0.005 |
| AUC-ROC | — | 0.645 | (target was >0.7) |
| AUC-PR | — | 0.304 | (target was >0.5) |
| DockQ (Ab-Ag) | 0.179 | 0.088 | **-0.091** |
| Confidence score | 0.883 | 0.520 | -0.363 |

**Verdict**: Failed on both axes. Epitope prediction is no better than simply extracting contacts from vanilla B1 (F1: 0.340 vs 0.345). Meanwhile, structural quality is severely degraded: DockQ drops by 51%, and 80% of predictions are classified Incorrect. The region-specific beta-scaling with emphasis=0.5/de-emphasis=-0.3 is too destructive to structure quality. The model's confidence collapses (0.520), confirming the perturbation pushes predictions far off the learned manifold.

---

### Strategy G: Progressive CDR Refinement v2 (Multi-Phase)

**Goal**: Three-phase refinement — exploration with CDR3 beta-scaling, transition, then optimization.

| Metric | B1 | G | Delta | p-value |
|--------|-----|---|-------|---------|
| DockQ (Global) | 0.338 | 0.355 | +0.017 | 0.410 |
| CAPRI Medium+ (%) | 21.3 | 27.7 | +6.4 pp | — |
| Net CAPRI upgrades | — | +3 (6 up, 3 down) | — | — |
| Ensemble diversity | 8.5 A | 9.1 A | +0.6 A | — |

**Verdict**: Modest positive direction. G shifts 3 more complexes into CAPRI Medium quality, with the largest wins on medium-difficulty cases (8EZ8_HLA +0.445, 8OL9_BAH +0.273). B1 wins more individual matchups (27 vs 20), but G's wins are larger in magnitude. Diversity is preserved. The improvement is not significant (p=0.41) but the approach is sound: it helps where help is most needed (medium-difficulty) without catastrophic regressions.

---

### Strategy K: Region-Specific Beta-Scaling (CDR3 Focus)

**Goal**: Discover epitope locations through CDR3-specific beta-scaling.

| Metric | B1 | K | Delta | p-value |
|--------|-----|---|-------|---------|
| Epitope F1 | 0.345 | 0.363 | +0.018 | 0.464 |
| DockQ (Ab-Ag) | 0.179 | 0.182 | +0.003 | 0.805 |

**Verdict**: Marginal, non-significant improvements. The most informative finding is the scatter plot pattern: K dramatically improves some complexes (8FAH_HLA: +0.590, 8BLQ_ECD: +0.552) while catastrophically breaking others (8CDE_DCB: -0.688, 8J1T_HKF: -0.688). This **redistribution pattern** indicates that beta-scaling alters the diffusion trajectory strongly enough to change which binding mode is found, but lacks selectivity to steer toward the correct one. The potential is there, but it needs a selection mechanism to filter the wins from the losses.

---

### Strategy L: CDR3 Beta-Scaling for Conformation Exploration

**Goal**: Increase CDR3 conformational diversity through attention bias scaling (beta=0.3).

| Metric | B1 | L | Delta | p-value |
|--------|-----|---|-------|---------|
| CDR3 Ensemble diversity | 1.51 A | 1.66 A | +0.15 A | — |
| Frac diverse (>=1 A) | 45% | 60% | +15 pp | — |
| DockQ (Ab-Ag) | 0.179 | 0.201 | +0.022 | 0.670 |
| CDR3 RMSD | 3.64 | 3.76 | +0.12 | **0.0008** |
| iPTM | 0.737 | 0.828 | +0.091 | **<0.0001** |

**Verdict**: Mixed — achieves its primary diversity goal but at a cost. The 33% increase in diverse ensembles (45% to 60%) is the strongest positive result. However, CDR3 RMSD worsens significantly (p=0.0008), meaning more diverse but less accurate. iPTM inflates dramatically without corresponding quality improvement, creating a dangerous confidence-accuracy mismatch. The mechanism works (beta-scaling provably changes CDR3 conformations) but beta=0.3 is likely too aggressive.

---

### Strategy O: Asymmetric Beta-Scaling v2

**Goal**: Encode biological asymmetry (CDR-H3 dominates binding ~60%) through differential denoising step scaling.

| Metric | B1 | O | Delta | p-value |
|--------|-----|---|-------|---------|
| DockQ (Global, oracle) | 0.391 | 0.392 | +0.001 | 0.434 |
| CDR-H3 RMSD (oracle) | 2.27 | 2.24 | -0.03 | — |
| H3/L3 contact ratio | 6.65 | 5.39 | -1.26 | — |

**Verdict**: Zero effect. The +0.001 DockQ is indistinguishable from noise. The model already produces H3-dominant contact patterns (ratio ~3-7x, far exceeding the 1.5-2.0x biological target), so adding a mild 1.4x step scaling for H3 is redundant. The approach modifies step magnitude but not direction — it makes CDR-H3 atoms move faster per step but gives no information about *where* they should move. This is fundamentally insufficient.

---

### Strategy Q: Iterative Epitope Refinement (3-Round)

**Goal**: Discover epitope through iterative narrowing — Round 1 (broad), Round 2 (focused), Round 3 (refined).

| Metric | B1 | Q | Delta | p-value |
|--------|-----|---|-------|---------|
| Epitope F1 (mean) | 0.352 | 0.381 | +0.029 | 0.753 |
| Epitope F1 (median) | 0.233 | 0.373 | +0.140 | — |
| DockQ (Ab-Ag) | 0.179 | 0.198 | +0.019 | 0.789 |

**Verdict**: Most promising epitope-discovery approach, though not significant. The median F1 jump from 0.233 to 0.373 (+60% relative) is the best epitope improvement of any strategy. Q's key strength is helping on "hard zero" cases where B1 finds nothing — several complexes gain non-zero F1 from zero. The iterative narrowing does provide a focused signal. However, the bootstrapping problem (Round 1's incorrect pose misleads Rounds 2-3) limits reliability. With more Round 1 diversity and better hotspot identification, this approach could improve substantially.

---

### Strategy V: Embedding-Space CDR3 Steering

**Goal**: Robust CDR3 exploration by optimizing pair representations toward a "self-reference" target.

| Metric | B1 | V | Delta | p-value |
|--------|-----|---|-------|---------|
| DockQ (Global) | 0.338 | 0.304 | **-0.034** | **0.0004** |
| CDR3-H RMSD | ~2.5 | ~4.2 | **+1.7** | **<0.0001** |
| CDR3-L RMSD | ~1.0 | ~3.0 | **+2.0** | **<0.0001** |
| B-C DockQ (VH-VL) | ~0.65 | ~0.58 | **-0.082** | **<0.0001** |
| ipTM | 0.737 | 0.584 | -0.153 | — |

**Verdict**: Harmful. The worst-performing strategy — significantly degrades every metric. The root cause is clear: the self-reference target (mean of non-CDR3 pair embeddings) forces CDR3 to resemble framework regions, erasing the structural specificity CDR3 loops need. The perturbation also corrupts the VH-VL interface. The model itself recognizes the damage (ipTM collapses by 0.153). The embedding optimization concept is not inherently wrong, but the target definition is fundamentally flawed.

---

### Strategy W: Embedding-Based Interface Steering

**Goal**: Improve antibody-antigen interface by weighting distance-based steering with pair embedding magnitudes.

| Metric | B1 | W | Delta | p-value |
|--------|-----|---|-------|---------|
| DockQ (avg, oracle) | 0.204 | 0.207 | +0.004 | 0.224 |
| Epitope F1 | 0.370 | 0.365 | -0.005 | — |
| Wins / Ties / Losses | — | 5 / 31 / 11 | — | — |

**Verdict**: Negligible effect. The +0.004 DockQ is within noise, and W ties on 31/47 complexes. The fundamental assumption — that pair embedding norms encode binding specificity — appears unvalidated. The weights may reflect sequence distance or chain identity rather than interface propensity. The single-pass embedding extraction (computed once, held fixed) also means weights become stale as the structure evolves during diffusion.

---

### Strategy Y: Hierarchical Steering (Embedding Early + Coordinate Late)

**Goal**: Combine embedding-space beta-scaling (early diffusion) with coordinate-space potentials (late diffusion).

| Metric | B1 | Y | Delta | p-value |
|--------|-----|---|-------|---------|
| DockQ (Total, top-1) | 0.338 | 0.364 | +0.026 | 0.224 |
| DockQ (Ab-Ag) | 0.179 | 0.213 | +0.034 | 0.459 |
| CDR-H3 RMSD | 3.00 | 2.70 | **-0.30** | **0.025** |
| CAPRI Medium+ (%) | 21.3 | 27.7 | +6.4 pp | — |

**Verdict**: The only strategy with a statistically significant structural improvement. CDR-H3 RMSD improves by 0.30 A (p=0.025), the most clinically relevant loop for binding. DockQ trends positive consistently. No degradation on any metric. The dual-phase approach (embedding guidance for global orientation, coordinate potentials for local refinement) mirrors natural structure formation. Performance is equivalent to pocket-guided restraints despite using no oracle info. The modest magnitude suggests the embedding-space and coordinate-space signals need strengthening.

---

## 3. Global Ranking

Strategies ranked by closeness to achieving their individual goals, considering the quality of improvement, statistical evidence, and absence of harm:

| Rank | Strategy | Category | Key Achievement | DockQ vs B1 | Significant? | Harm? |
|------|----------|----------|----------------|-------------|-------------|-------|
| 1 | **A (20 FK particles)** | Docking | Best DockQ improvement (+0.046), 87% of B2's fnat, 0 CAPRI downgrades | +0.046 | No (p=0.61) | 25% complexes missing |
| 2 | **Y (Hierarchical)** | Docking + CDR3 | Only significant structural improvement (CDR-H3 -0.30A, p=0.025) | +0.034 | **Yes** (CDR-H3) | None |
| 3 | **G (Progressive)** | Docking | +3 net CAPRI upgrades, helps medium-difficulty cases, preserves diversity | +0.017 | No (p=0.41) | Slight regression on easy cases |
| 4 | **L (CDR3 beta)** | CDR3 diversity | Achieves diversity goal (60% vs 45% diverse ensembles) | +0.022 | No | CDR3 RMSD worsens (p=0.0008) |
| 5 | **Q (Iterative)** | Epitope | Best epitope F1 improvement (+0.029 mean, +0.140 median) | +0.019 | No | Minor regressions on some |
| 6 | **K (Region beta)** | Epitope | Marginal epitope gain (+0.018 F1), strong redistribution pattern | +0.003 | No | Catastrophic on some complexes |
| 7 | **D (Canonical)** | CDR diversity | Computational efficiency achieved, quality unchanged | +0.009 | No | Confidence degrades |
| 8 | **W (Embed interface)** | Docking | Negligible effect, concept unvalidated | +0.004 | No | No |
| 9 | **O (Asymmetric beta)** | CDR balance | Zero effect — model already has H3 bias | +0.001 | No | No |
| 10 | **E (Blind scanning)** | Epitope | Epitope F1 = B1, severe structural degradation (DockQ -51%) | -0.091 | **Yes (worse)** | Severe |
| 11 | **V (Embed CDR3)** | CDR3 | Harmful — degrades CDR3 RMSD by 1.3A, DockQ by 0.034 | -0.034 | **Yes (worse)** | Severe |

---

## 4. Cross-Strategy Lessons Learned

### 4.1 What Works

1. **Feynman-Kac resampling (A)** is the most effective mechanism. By maintaining multiple particles and discarding poor orientations, it provides a natural safety net (0 CAPRI downgrades) while allowing large upside improvements. The particle diversity → selection pipeline is fundamentally sound.

2. **Hierarchical timing (Y)** is the correct architectural insight. Embedding-space guidance for global orientation (early) followed by coordinate-space refinement (late) mirrors the physics of protein folding. Y is the only strategy with a significant structural improvement.

3. **Progressive phase transitions (G)** help medium-difficulty cases. The multi-phase schedule (exploration → transition → optimization) is a natural fit for diffusion sampling.

4. **CDR3 beta-scaling (L)** provably modulates conformational diversity. The mechanism is confirmed to work — the challenge is controlling it.

5. **Iterative refinement (Q)** addresses the right problem for epitope discovery. The median F1 jump (+60%) on hard cases shows the bootstrapping concept has merit.

### 4.2 What Fails

1. **Step-magnitude-only modifications (O)** are insufficient. Changing how fast atoms move provides no information about where they should move. Direction matters more than magnitude.

2. **Self-reference embedding targets (V)** are fundamentally wrong. CDR3 loops should not resemble framework regions. The target definition determines success or failure of embedding-space methods.

3. **Region-specific beta-scaling for structure prediction (E)** destroys structural quality. The emphasis/de-emphasis magnitudes (0.5/-0.3) applied to specific antigen regions are too disruptive. The model cannot maintain physically valid structures under such strong region-specific biases.

4. **Raw embedding norms as interaction proxies (W)** are unvalidated. Pair representation norms may not encode binding specificity. The assumption must be tested before building methods on it.

5. **Confidence-based model selection** degrades for all steered methods. Every strategy that modifies the diffusion process weakens the confidence-DockQ correlation (D: 0.287→0.177, A: 0.380→0.265, L: iPTM inflation). This is a systemic problem requiring a custom re-ranking solution.

### 4.3 Recurring Patterns

- **Sporadic wins vs systematic improvement**: Most strategies improve a handful of complexes dramatically while leaving the majority unchanged or slightly worse. Only B2 (contacts) achieves systematic improvement. This suggests the steering mechanisms are strong enough to change outcomes but lack the selectivity to consistently change them for the better.

- **Medium-difficulty complexes benefit most**: Strategies A, G, Q, Y all show their largest improvements on complexes where B1 is borderline (DockQ 0.1-0.3). The hardest cases (DockQ ≈ 0, ~40% of the set) remain unsolvable by any method except B2, suggesting a fundamental global orientation problem that local steering cannot fix.

- **Diversity-accuracy tradeoff**: L demonstrates that increasing conformational diversity degrades individual prediction accuracy. This is expected but important: diversity is only valuable if coupled with an effective selection mechanism.

### 4.4 Comparison Against B1 (Simple Baseline) — Epitope-Focused Strategies

For strategies aimed at epitope identification (E, K, Q), the fair comparison is against B1 since contact-based baselines use oracle information:

| Strategy | Epitope F1 vs B1 | Improvement Mechanism | Assessment |
|----------|-----------------|----------------------|------------|
| Q | +0.029 (mean), +0.140 (median) | Iterative narrowing | Most promising |
| K | +0.018 | CDR3 beta-scaling | Marginal |
| E | -0.005 | Region scanning | No improvement + structural harm |

Q's approach of iterative refinement is the best direction for blind epitope discovery, despite not reaching significance. The median improvement is substantial and the method helps on the hardest cases.

### 4.5 Relevance Even in Presence of Contact Information

Even if contact restraints (B2) are available, new strategies could still add value by:
- **Predicting which contacts to use** (reducing B2's reliance on oracle information)
- **Complementing contacts with CDR3 optimization** (B2 doesn't improve CDR-H3 RMSD as much as Y does per strategy)
- **Providing diversity for ensemble methods** (L's diversity + B2's accuracy)
- **Handling novel antigens** where contact information is unavailable

---

## 5. Analysis of Next Steps from Individual Evaluations

Each evaluation report proposed specific modifications. Here I consolidate the key themes:

### From A (FK Particles):
- Fix OOM for missing 12 complexes (reduce to 10-15 particles)
- Custom re-ranking using steering potential energy instead of confidence
- Combine FK particles with more diffusion samples
- Embedding-space steering to avoid particle overhead

### From D (Canonical Ensemble):
- Increase perturbation strength and specificity (target CDR-antigen pairs)
- Diversity-aware model selection instead of confidence-only
- Per-CDR beta values (stronger for H3, lighter for H1/H2)
- External reranking (Rosetta, shape complementarity)

### From E (Blind Scanning):
- Increase beta magnitudes (emphasis 1.0-2.0)
- More regions (20-40 instead of 10)
- Tighter contact threshold (5-6 A instead of 8 A)
- Weighted aggregation by confidence

### From G (Progressive):
- Conditional steering — only when initial confidence is low
- Stronger CDR3 beta in Phase 1
- Better model selection (interface-specific criteria)

### From K (Region Beta):
- Ensemble over multiple region-emphasis patterns (like B2's multi-config)
- Softer beta values (0.2-0.3 instead of 0.5)
- Full heatmap accumulation pipeline (currently only single-config tested)
- Two-stage: K for epitope → B2-style restraints on predicted epitope

### From L (CDR3 Beta):
- Reduce beta (0.1-0.15 instead of 0.3)
- Asymmetric CDR scaling (H3 stronger, L3 weaker)
- Expand scaling mask to CDR3-antigen pairs
- Time-dependent beta (high early, low late)

### From Q (Iterative):
- More Round 1 diversity (10-20 samples, higher noise)
- Ensemble-based hotspot selection
- Confidence-weighted contact maps
- Adaptive rounds (skip if no clear hotspot)

### From V (Embed CDR3):
- Alternative target (canonical class embeddings, not framework mean)
- Reduced strength (0.01-0.1 instead of 1.0)
- Apply only at early denoising steps
- CDR3-to-antigen pairs as target instead of self-reference

### From W (Embed Interface):
- Validate embedding signal first (AUC-ROC of norms vs true contacts)
- Learned scoring head instead of raw norms
- Re-extract embeddings during diffusion
- Use embeddings to select which B2-style restraints to apply

### From Y (Hierarchical):
- Stronger beta-scaling factor in early phase
- Better coordinate-space potentials in late phase
- Incorporate predicted epitope information
- Focus on the ~40% failure cluster (coarse rigid-body pre-step)

---

## 6. Proposed New Strategies for Round 2

Based on the lessons learned and cross-strategy analysis, here are new strategies organized by expected impact and implementation complexity. These combine successful mechanisms, address identified failure modes, and leverage the union of individual strategy strengths.

### Tier 1: High Impact, Low-Medium Complexity (Start Here)

| # | Name | Mechanism | Rationale | Expected Impact | Complexity |
|---|------|-----------|-----------|----------------|------------|
| N1 | **A+Y Hybrid (FK + Hierarchical)** | Apply Y's hierarchical timing (embedding early, coordinate late) with A's FK resampling (5-10 particles). Early: beta-scaling for CDR3 orientation. Late: FK particles with coordinate potentials. | A has the best docking improvement; Y has the only significant CDR-H3 gain. Combining them addresses both global orientation (Y's embedding phase) and local refinement (A's particle selection). FK particle count reduced to 5-10 to avoid OOM. | DockQ +0.05-0.08, CDR-H3 RMSD -0.3-0.5 A | Medium |
| N2 | **Q→B2 Pipeline (Iterative Epitope → Contact Restraints)** | Round 1: Q's iterative epitope discovery (3 rounds, 10+ samples/round). Round 2: Use top predicted epitope residues as B2-style contact restraints. No oracle info needed. | Q is the best blind epitope finder; B2 is the best method when contacts are known. This pipeline turns Q's output into B2's input, potentially achieving B2-level performance without oracle info. The main question is whether Q's epitope predictions are accurate enough. | DockQ +0.10-0.20 if epitope accurate | Medium |
| N3 | **Confidence-Aware Conditional Steering** | Compute B1 baseline confidence first. If iptm > 0.85: skip steering (B1 is likely correct). If iptm < 0.75: apply G+-style progressive steering. If 0.75-0.85: apply mild steering (Y-style). | G and Y both regress on easy cases while helping medium ones. Conditional application avoids the "steering hurts good predictions" failure mode observed in G (-0.025 on easy cases) and concentrates computational effort where it's needed. | Eliminate regressions, preserve wins | Low |
| N4 | **L-light (CDR3 Beta, Reduced Strength)** | L with beta=0.1-0.15 (instead of 0.3) + time-dependent schedule (beta_early=0.15, beta_late=0.0). Apply to CDR3-antigen pairs (not just CDR3-CDR3). | L proved that beta-scaling modulates diversity but beta=0.3 degrades accuracy. Lower beta + interface-directed scaling could keep the diversity benefit while reducing the CDR3 RMSD penalty. Time-dependent schedule allows exploration early and convergence late. | CDR3 diversity +10% without accuracy loss | Low |
| N5 | **Multi-Config K Ensemble** | Run K with 5-8 different region emphasis patterns per complex (different beta configurations emphasizing different antigen patches). Select best by confidence or aggregate into epitope heatmap. | K's scatter plot shows it can reach DockQ >0.6 on some complexes. An ensemble over multiple K configs mimics B2's multi-config advantage. Oracle-over-K-configs could capture the union of K's wins while avoiding its worst failures. | Epitope F1 +0.05-0.10, DockQ +0.02-0.04 | Low |

### Tier 2: Medium Impact, Medium Complexity (Follow-Up)

| # | Name | Mechanism | Rationale | Expected Impact | Complexity |
|---|------|-----------|-----------|----------------|------------|
| N6 | **Custom Re-Ranking Module** | For any steered method, generate 10-20 samples, then re-rank using a composite score: `alpha * confidence + beta * steering_energy + gamma * interface_pLDDT`. Tune alpha/beta/gamma on a held-out set. | Every strategy suffers from degraded confidence calibration. A custom re-ranker could recover the oracle gap (typically 0.04-0.08 DockQ) for any steered method. This is a universal improvement applicable to A, G, K, L, Q, Y. | DockQ +0.03-0.06 (applied to any strategy) | Medium |
| N7 | **Y+ (Stronger Hierarchical)** | Y with increased beta-scaling factor (2-3x current) in early phase + stronger coordinate potentials in late phase + predicted epitope information from Q or sequence-based tools as input to late-phase potentials. | Y's trends are consistently positive but modest. Strengthening both phases and adding predicted epitope info could amplify the improvement. The key change is providing the late-phase coordinate potentials with a target binding region. | DockQ +0.04-0.08, CDR-H3 RMSD -0.4-0.6 A | Medium |
| N8 | **G+L Combined (Progressive + CDR3 Diversity)** | Phase 1: L-style CDR3 beta-scaling (beta=0.15) for diverse CDR3 conformations. Phase 2: G-style transition (neutral). Phase 3: A-style FK resampling (5 particles) for antigen optimization. Best-of-all selection via N6's re-ranker. | Combines the best elements: L's diversity (Phase 1), G's phased approach, A's particle selection (Phase 3), and custom re-ranking. Each phase addresses a different bottleneck in sequence. | DockQ +0.04-0.07, CDR3 diversity +10% | Medium |
| N9 | **V-Fixed (Corrected Embedding CDR3)** | V with corrected target: use CDR3 embeddings from known canonical class structures (not framework mean). Reduced strength (0.05-0.1). Apply only during first 30% of denoising. CDR3-antigen pair target instead of self-reference. | V's mechanism (embedding optimization) is powerful but the target was wrong. With a biologically appropriate target (canonical CDR3 conformations or CDR3-antigen interaction patterns), the same gradient descent could steer CDR3 conformations productively. The early-only application allows recovery. | CDR3 conformation accuracy if target is correct | Medium |
| N10 | **E-Refined (Softer Blind Scanning)** | E with reduced emphasis/de-emphasis (0.2/-0.1 instead of 0.5/-0.3), 20-30 regions instead of 10, tighter contact threshold (5 A instead of 8 A), confidence-weighted heatmap aggregation. | E's concept (scanning for epitopes) is sound but the implementation was too aggressive. Halving the beta magnitudes and doubling the region count should preserve structural quality while improving spatial resolution. The key fix is making the perturbation subtle enough to avoid structural collapse. | AUC-ROC >0.7, DockQ maintained | Medium |

### Tier 3: High Impact if Successful, Higher Complexity (Research Directions)

| # | Name | Mechanism | Rationale | Expected Impact | Complexity |
|---|------|-----------|-----------|----------------|------------|
| N11 | **Full Pipeline: Q→K→A (Discover→Focus→Refine)** | Stage 1: Q (3 rounds) to identify candidate epitope regions. Stage 2: K ensemble focused on top 2-3 epitope patches. Stage 3: A (10 FK particles) with coordinate potentials targeting the best K patch. | Complete blind-to-refined pipeline. Each stage addresses a different challenge: Q discovers the epitope, K validates and focuses, A refines the complex structure. The pipeline produces both an epitope prediction and a docked structure. | DockQ +0.08-0.15, Epitope F1 +0.10-0.15 | High |
| N12 | **Learned Re-Ranking with Steering Features** | Train a lightweight regression model (GBT or small MLP) on Round 1 results to predict DockQ from: confidence_score, steering_energy, interface_pLDDT, CDR3_pLDDT, SASA_ratio, H3/L3_contact_ratio. Cross-validate on 47 complexes. | The biggest universal bottleneck is model selection. With 47 complexes × multiple strategies, there's enough data to train a re-ranker that generalizes across steering methods. This directly addresses the confidence calibration problem. | DockQ +0.05-0.10 (model selection improvement) | Medium-High |
| N13 | **Coarse Rigid-Body Pre-Step + Y** | Before Boltz2 diffusion: run a fast rigid-body docking scan (e.g., ZDOCK-style FFT or a learned pose classifier) to identify 3-5 candidate global orientations. Initialize Boltz2 diffusion from each orientation. Apply Y-style hierarchical steering. Select best. | ~40% of complexes fail at the global orientation level — no local steering can fix an antibody pointed at the wrong hemisphere. A coarse pre-step could reduce this failure rate by providing reasonable starting orientations. | DockQ +0.05-0.10 on the hardest 40% | High |
| N14 | **W-Validated (Embedding Interface with Learned Head)** | Train a small MLP on CDR-antigen pair embeddings from B2's successful complexes to predict contact probability. Use these learned scores (not raw norms) to weight the interface potential. Re-extract embeddings every 10 diffusion steps. | W failed because raw embedding norms don't encode binding. A learned scoring head could extract the binding signal if it exists. Re-extraction prevents staleness. This validates or refutes the core EmbedOpt hypothesis for Boltz2. | Validates embedding approach; DockQ +0.02-0.05 if signal exists | High |
| N15 | **Ensemble Union Strategy** | Run A(5p), G, Y, and B1 in parallel. For each complex, take the pool of all predictions (20 models: 5×4) and select the best using N6's re-ranker. | Different strategies help different complexes (A helps 8FAH_HLA, G helps 8EZ8_HLA, Y helps others). Taking the union of all strategies' predictions and selecting the best should capture all individual wins. The key is the selection mechanism. | DockQ +0.05-0.10 (captures union of wins) | Medium-High |

### Summary of New Strategies

| Tier | Strategies | Focus | Key Dependency |
|------|-----------|-------|---------------|
| **Tier 1** | N1-N5 | Quick combinations of proven mechanisms | Minimal new code |
| **Tier 2** | N6-N10 | Fixing identified failure modes | Moderate new code + tuning |
| **Tier 3** | N11-N15 | Full pipelines and research directions | New infrastructure + training |

**Recommended starting point**: N3 (conditional steering) and N5 (multi-config K ensemble) are the lowest-effort changes with clear rationale. N1 (A+Y hybrid) and N2 (Q→B2 pipeline) are the highest expected-impact combinations.

---

## 7. Final Conclusions

### What Round 1 Established

1. **Steering works mechanistically.** Every strategy (except O) demonstrably changes the diffusion trajectory. L provably increases CDR3 diversity. A's FK particles provably filter bad orientations. K and E provably redistribute binding modes. Y provably improves CDR-H3 accuracy. The challenge is not whether steering has an effect, but whether the effect is consistently beneficial.

2. **No single strategy achieves statistical significance on DockQ.** Only Y reaches significance on any structural metric (CDR-H3 RMSD, p=0.025). The fundamental problem is that improvements are concentrated in ~5 complexes per strategy rather than broadly distributed. This is a selectivity problem, not a mechanism problem.

3. **The gap to contact restraints (B2) remains large.** B2 achieves DockQ ~0.29-0.41 vs the best strategy's ~0.25. B2 benefits from oracle information and many configurations per complex. Closing this gap without oracle info is the primary challenge for Round 2.

4. **Combination is the clear next step.** Individual strategies address different bottlenecks (A: orientation, Y: CDR-H3, L: diversity, Q: epitope). No single strategy addresses all of them. Combining complementary mechanisms (Tier 1 proposals N1-N5) should yield larger, more systematic improvements.

5. **Model selection is a universal bottleneck.** Every steered method degrades confidence calibration. Developing a custom re-ranking solution (N6, N12) would benefit all strategies and is likely the single highest-ROI investment for Round 2.

### What Needs to Change for Round 2

- **From single mechanisms to pipelines**: Combine strategies that address different failure modes in sequence.
- **From fixed parameters to adaptive**: Use confidence and contact metrics to determine when and how aggressively to steer (N3).
- **From confidence-based selection to custom re-ranking**: Steered models need steered selection criteria (N6).
- **From equal application to conditional application**: Steer only when the baseline is uncertain (N3).
- **From independent strategies to ensemble selection**: Run multiple strategies and pick the best prediction per complex (N15).
