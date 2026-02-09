# Boltz-2 Steering Methods

## Overview

Steering methods in Boltz-2 allow you to guide structure prediction toward desired conformations during the diffusion inference process. Instead of predicting structures based solely on the trained model, steering applies **physical and biological constraints** that nudge the predictions to satisfy specific criteria.

### General Goals

- **Controllability**: Direct the model to generate structures with specific properties
- **Biological Accuracy**: Enforce known biological constraints (CDR-antigen contacts, loop conformations)
- **Physical Validity**: Maintain chemically realistic structures while satisfying constraints
- **Exploration**: Generate multiple samples that satisfy constraints in different ways

### When to Use Steering

- You have **prior knowledge** about the desired structure (e.g., known CDR3 loops from experiments)
- You want to **optimize for specific properties** (e.g., better antigen binding)
- You're working with **constrained design** problems
- You want to **explore alternative conformations** that meet certain criteria

---

## Available Steering Methods

### 1. [[CDR3 Steering|notes_steering/cdr3_steering]]

**Flag**: `--cdr3_steering`

**Goal**: Control the backbone conformation of complementarity-determining region 3 (CDR3) loops in antibodies.

**What it does**: Guides the model to generate CDR3 loops with specific dihedral angles (psi angles), allowing you to enforce known or desired loop geometries.

→ [Learn more about CDR3 Steering](./notes_steering/cdr3_steering.md)

---

### 2. [[Antigen Steering|notes_steering/antigen_steering]]

**Flags**: `--antigen_steering` + `--num_particles` (optional)

**Goal**: Optimize the 3D orientation and position of antigens relative to antibody CDR loops to maximize productive contacts.

**What it does**: Applies spatial constraints to drive the antigen toward geometries that maximize CDR-antigen binding interactions. Uses importance sampling with multiple particles to explore the space.

→ [Learn more about Antigen Steering](./notes_steering/antigen_steering.md)

---

## Common Usage Patterns

### Basic Steering

All steering methods require `--use_potentials` to be enabled:

```bash
boltz predict config.yaml --use_potentials --cdr3_steering
boltz predict config.yaml --use_potentials --antigen_steering
```

### Multiple Independent Samples

Generate N different structures that all satisfy the steering constraints:

```bash
boltz predict config.yaml --use_potentials --cdr3_steering --diffusion_samples 5
```

Each sample explores different valid conformations while respecting the constraints.

### Particle-Based Exploration (Antigen Steering)

Antigen steering can use `--num_particles` for more thorough exploration:

```bash
boltz predict config.yaml --use_potentials --antigen_steering --num_particles 5
```

Higher particle counts = better exploration but higher computational cost.

---

## How Steering Works (Technical Overview)

All steering methods operate during the **iterative denoising phase** of diffusion inference:

```
Input: Protein + Antigen sequences
  ↓
Network Trunk: Generate learned representations
  ↓
Diffusion Loop: Iteratively denoise structure (t=T → t=0)
  ├─ Step 1: Network predicts cleaner coordinates
  ├─ Step 2: [STEERING APPLIED] Constraints modify coordinates
  ├─ Step 3: Noise added back
  └─ Repeat until t=0
  ↓
Output: Final structure satisfying constraints
```

### Two Core Mechanisms

**Gradient Guidance**: At each denoising step, constraint energy gradients push atomic coordinates toward satisfying the constraints.

**Importance Sampling** (optional): Multiple particle copies explored in parallel, weighted by how well they satisfy constraints, poorly performing copies discarded.

---

## Steering Methods Comparison

| Method | Purpose | Particle Support | Key Parameter |
|--------|---------|------------------|----------------|
| [[CDR3 Steering\|notes_steering/cdr3_steering]] | Loop backbone control | No | Target psi angles |
| [[Antigen Steering\|notes_steering/antigen_steering]] | Binding orientation | Yes | `--num_particles` |

---

## Next Steps

- **Start with**: [[CDR3 Steering|notes_steering/cdr3_steering]] for loop design
- **Then try**: [[Antigen Steering|notes_steering/antigen_steering]] for binding optimization
- **See also**: Check `--help` for all steering-related flags

```bash
boltz predict --help | grep -i steer
```

---

## Future Extensions

This framework supports adding more steering methods in the future:

- Secondary structure steering (enforce alpha/beta regions)
- Domain orientation steering (control multi-domain arrangements)
- Loop conformation steering (beyond CDR3)
- Custom energy function steering (user-defined constraints)

---

## Steering Proposals

New steering strategies are under development to improve antibody-antigen binding prediction, particularly for cases with unknown epitope locations. These strategies aim to:

- Explore different CDR conformations (especially CDR-H3 and CDR-L3)
- Simultaneously optimize antigen orientation without knowing the epitope location
- Discover which regions of the antigen are preferred for binding

**[→ See detailed steering proposals and novel strategies](./notes_steering/steering_proposals.md)**

### Summary of Proposed Strategies

- **Strategy A-C**: Existing approaches (reference)
- **Strategy D (Canonical Ensemble)**: Sample from known CDR canonical structures
- **Strategy E (Blind Scanning)**: Systematically scan antigen surface to find epitope
- **Strategy F (Ensemble Modes)**: Generate multiple binding conformations with different CDR geometries
- **Strategy G (Multi-Phase)**: Progressive refinement from exploration to exploitation
- **Strategy H (Entropy-Maximizing)**: Generate diverse conformations with equal binding quality
- **Strategy I (Epitope Scoring)**: Integrate modern ML epitope predictions
- **Strategy J (Coupled CDR3)**: Joint sampling of heavy and light chain CDR3

Each strategy includes implementation complexity analysis, expected benefits, and computational cost estimates.

### Latent Space Scaling-Inspired Ideas

Inspired by recent work on pair representation scaling (Suzuki & Amagasa, 2026), we are exploring novel steering methods that modulate the latent space directly during inference:

**[→ See latent space scaling-inspired steering ideas](./notes_steering/latent_space_scaling_inspired_ideas.md)**

These approaches include:
- **Idea K**: Region-Specific β-Scaling for Epitope Discovery
- **Idea L**: CDR3-Specific β-Scaling for Conformation Exploration
- **Idea M**: Inter-Chain β-Scaling for Interface Optimization
- **Idea N**: Adaptive β-Scheduling During Diffusion
- **Idea O**: Differential β-Scaling Between Heavy and Light Chain CDRs
- **Idea P**: Contact-Guided β-Scaling Feedback Loop

All latent space scaling ideas share key advantages: **computational efficiency**, **no retraining required**, **transparent mechanism**, and **MSA-agnostic** (work with or without evolutionary information).

**Status**: Design Phase - These are proposals for future implementation.

---

## Improvements to Steering Proposals and Novel Ideas

Based on analysis of the FK-Diffusion Steering framework (Horvitz et al., 2501.06848) and Boltz-sample β-scaling (Suzuki & Amagasa, 2026), we have identified significant improvements to the existing strategies (A-J) and enhancements to the β-scaling ideas (K-P), plus several entirely novel hybrid approaches.

### Quick Start
**New Analysis (Feb 2025)**: Based on FK-Diffusion Steering and Boltz-sample papers, we've identified 22 significant improvements and 5 novel ideas. Start here:

1. **[→ Analysis Summary](./notes_steering/ANALYSIS_SUMMARY.md)** - Overview of what was analyzed and key findings
2. **[→ Quick Reference](./notes_steering/improvements_quick_reference.md)** - At-a-glance summary, implementation roadmap, decision guide
3. **[→ Full Technical Details](./notes_steering/improvements_and_new_ideas.md)** - Comprehensive analysis with implementation notes for each improvement
4. **[→ Implementation Guide](./notes_steering/implementation_guide.md)** - Concrete code changes needed for each improvement (without making changes yet)

### Summary of Improvements to Existing Strategies

- **Strategy D+ (Enhanced Canonical Ensemble)**: Use β-scaling instead of explicit potentials for better efficiency
- **Strategy E+ (Enhanced Blind Scanning)**: Use region-specific β-scaling (Idea K) for 10-20x speedup
- **Strategy G+ (Enhanced Progressive Refinement)**: Implement contact-triggered phase transitions and smooth β-scheduling
- **Strategy J+ (Enhanced Coupled CDR3)**: Replace 2D grid with asymmetric β-scaling and beam search

### Summary of Enhancements to β-Scaling Ideas

- **Idea K+ (Region-Specific β-Scaling v2)**: Adaptive region definition using surface properties and confidence weighting
- **Idea L+ (CDR3-Specific β-Scaling v2)**: Two-stage approach combining CDR3 exploration with contact optimization
- **Idea M+ (Interface β-Scaling v2)**: Integration with predicted binding affinity and physics-informed β scheduling
- **Idea N+ (Adaptive β-Scheduling v2)**: Contact-aware dynamic phase transitions with plateau detection
- **Idea O+ (Asymmetric β-Scaling v2)**: Task-specific β values with experimental validation targets
- **Idea P+ (Contact-Guided β-Feedback v2)**: Multi-pass refinement with convergence criteria

### Summary of Novel Hybrid Approaches

- **Idea Q (Iterative Epitope Refinement)**: Three-round FK + β-scaling combination for rapid epitope discovery and refinement
- **Idea R (Multi-Objective Steering)**: Simultaneous optimization of contact, conformation, and clash avoidance
- **Idea S (Template-Guided β-Scheduling)**: Incorporate known antibody-antigen structure knowledge via soft β-biasing
- **Idea T (MSA-Free Steering)**: Strategies for synthetic/designed antibodies with limited evolutionary information
- **Idea U (Conformational Energy Landscape Mapping)**: Map entire accessible structure space via β-sweep for complete characterization

### Implementation Priority

**Phase 1 (Immediate)**: L+, K+, G+, N+ - highest impact, moderate complexity
**Phase 2 (Near-term)**: Q, T, O+ - novel contributions, good scientific value
**Phase 3 (Medium-term)**: R, S, U - advanced approaches, publication-worthy

All improvements are designed to leverage the efficiency and transparency of β-scaling mechanisms and the flexibility of FK-Diffusion resampling frameworks.

---

## EmbedOpt-Inspired Improvements (Robust Embedding Space Steering)

Recent work on **"Robust Inference-Time Steering of Protein Diffusion Models via Embedding Optimization"** (Li et al., 2602.05285) provides critical insights for improving our steering strategies:

### Key Insight: Embedding Space vs. Coordinate Space

**Finding**: Optimizing in embedding space is significantly more robust than coordinate-space optimization, especially for constraints outside the training distribution.

**Benefits**:
- Stable across 2 orders of magnitude of hyperparameter variation
- Better handles out-of-distribution constraints (novel antibodies, antigens)
- Fewer aggressive parameter tuning requirements
- Captures sequence and coevolutionary signals naturally
- Integrates with experimental constraints more naturally

### New EmbedOpt-Inspired Strategies

**[→ See detailed embedding-space ideas V-Z in improvements_and_new_ideas.md](./notes_steering/improvements_and_new_ideas.md)**

These include:
- **Idea V**: Embedding-Space CDR3 Steering (more robust than dihedral potentials)
- **Idea W**: Embedding-Based Interface Steering (more stable than contact distances)
- **Idea X**: Experimental Constraint Integration (cryo-EM, cross-linking, HDX-MS data)
- **Idea Y**: Hierarchical Steering (embedding robustness + coordinate precision)
- **Idea Z**: Soft Experimental Constraints (probabilistic, uncertainty-aware)

### Implementation Priority (EmbedOpt Ideas)

**Phase 2.5 (Week 2-3)**: V, W - Medium complexity, good robustness improvement
**Phase 3 (Week 5+)**: X, Y, Z - High complexity, publication-worthy, experimental data integration

### Key Advantage Over Current Approaches

Current steering methods use coordinate-space potentials (Strategies A-C) or latent-space β-scaling (Ideas K-P). EmbedOpt-inspired methods provide:
1. **Better generalization** to novel/synthetic antibodies
2. **Natural integration** of experimental data (not just contacts)
3. **More stable hyperparameter behavior** (less tuning needed)
4. **Principled handling** of measurement uncertainty
5. **Hierarchical refinement** (coarse-grained embedding control → fine-grained coordinate control)

---

## Comprehensive Steering Methods Comparison Table

**Organized by Implementation Roadmap** with detailed evaluation metrics for each method.

### Phase 1 (Week 1-2): Quick Wins - Highest Impact, Moderate Complexity

| Name | Goal | Mechanism | Use Case | Key Insight | Complexity | Evaluation Metrics |
|------|------|-----------|----------|-------------|-----------|-------------------|
| **L+ (CDR3 β-Scaling v2)** | Explore CDR3 conformations robustly for novel sequences | CDR3-specific pair representation scaling: z[CDR3] *= (1+β); two-stage exploration then optimization | Synthetic/designed antibodies; novel CDR sequences without canonical conformations; generating conformational ensembles | β-scaling in latent space is simpler and more robust than explicit dihedral potentials; works without MSA | **LOW** | RMSD of CDR3 backbone vs. known structures (if available); pLDDT of CDR3 region; ensemble diversity (pairwise CDR3 RMSD ≥ 1.0 Å); contact score with antigen |
| **K+ (Region-Specific β-Scaling v2)** | Discover epitope location on unknown antigens without prior knowledge | Adaptive antigen surface partitioning (curvature/accessibility-based); region-specific β values: emphasize one region (+0.5), de-emphasize others (-0.3); single diffusion pass | Unknown epitope problems; antibody-antigen pairs where binding site unknown; discovering cryptic or secondary epitopes; high-throughput epitope mapping | Region-specific β-scaling in single pass achieves 10-20x speedup vs. running separate passes; confidence weighting improves robustness | **MEDIUM** | Epitope prediction accuracy: precision/recall vs. known crystal structure epitopes; contact heatmap correlation with experimental data (e.g., PISA interface); precision = true_contacts/predicted_contacts; recall = true_contacts/all_true_contacts; secondary epitope detection rate |
| **N+ (Adaptive β-Scheduling v2)** | Smooth transition from exploration to exploitation during diffusion based on actual contact progress | Contact-triggered phase transitions (plateau detection); smooth β(t) scheduling instead of step functions; continuous weighting progression | Any steering task; especially important when phase transition timing is system-dependent | Contact-aware scheduling responds to actual prediction quality (not fixed timesteps); smoother transitions reduce artifacts | **MEDIUM** | Contact score improvement rate per phase; phase transition timing (early/late detection); final contact score vs. baseline; structure quality (pLDDT) at phase transitions |
| **G+ (Progressive Refinement v2)** | Three-phase refinement combining CDR3 exploration with antigen optimization in natural progression | Phase 1 (loose): exploratory β for CDR3 + weak antigen guidance; Phase 2 (medium): neutral CDR3 + moderate contact guidance; Phase 3 (tight): optimizing β + strong contact guidance | General purpose antibody-antigen predictions; cases needing both CDR3 and antigen orientation optimization; balancing diversity and accuracy | Multi-phase approach mimics natural structure formation and allows different optimization strategies at different stages | **MEDIUM** | Contact score per phase (should increase monotonically); structural RMSD if native available; ensemble diversity preservation; pLDDT trajectory throughout diffusion |

### Phase 2 (Week 3-4): Novel Contributions - Medium Complexity, Scientific Value

| Name | Goal | Mechanism | Use Case | Key Insight | Complexity | Evaluation Metrics |
|------|------|-----------|----------|-------------|-----------|-------------------|
| **Q (Iterative Epitope Refinement)** | Rapidly discover epitope and generate high-quality binding structures in 3-pass iterative approach | Round 1 (broad exploration): FK particles with weak guidance + exploratory β; Round 2 (targeted): focus on hotspots from Round 1 with stronger guidance; Round 3 (validation): final optimization with validated parameters | Unknown epitope systems; computational efficiency important; publication-oriented epitope discovery | Three-stage FK+β approach combines exploration, targeting, and refinement; validates results progressively | **MEDIUM** | Round 1: contact distribution breadth (should cover multiple antigen regions); Round 2: convergence speed (should narrow to 2-3 hotspots); Round 3: final contact score; overall speed (should be < baseline E); epitope accuracy vs. known structures |
| **T (MSA-Free Steering)** | Handle synthetic/designed antibodies with no sequence homologs using intrinsic model priors | β-scaling on latent space features that work without MSA; pAE-based contact scoring instead of distance; confidence weighting based on pLDDT only | De novo antibody design; synthetic antibody libraries with no natural MSA; optimization of CDR libraries; novel epitope binders | Latent space β-scaling activates model's internalized priors; pAE is more robust than explicit distances when MSA absent | **MEDIUM** | pLDDT of designed structures (should maintain ≥ 50); predicted binding affinity if gold standard available; ensemble quality when MSA_size varies from 1 to 1000; robustness to MSA quality degradation |
| **O+ (Asymmetric β-Scaling v2)** | Encode biological prior that CDR-H3 dominates binding (≈60% contacts) while CDR-L3 provides support (≈40%) | Differential β values: β_H = +0.4 to +0.5 (emphasize heavy chain), β_L = -0.1 to +0.2 (lighter emphasis on light chain); applied selectively to H3 and L3 regions | General antibody-antigen predictions; especially important for maximizing binding affinity; library design where H3 is primary target | Asymmetric scaling respects known antibody biology; reduces wasted particles on incompatible H3/L3 pairs vs. symmetric approaches | **LOW** | H3 contact score vs. L3 contact score (ratio should be ≈1.5-2.0); ensemble quality ranking by contact score; comparison of β_H=β_L (symmetric) vs. β_H>β_L (asymmetric) |
| **D+ (Enhanced Canonical Ensemble)** | Systematically explore known CDR canonical conformations efficiently using β-scaling instead of explicit potentials | β-scaling applied to CDR1/2/L3 for canonical exploration; CDR-H3 free to explore sub-canonicals; soft β (±0.3) preserves flexibility | Antibodies with canonical CDR1/2/L3; when canonical structure database is available; cases needing systematic conformation search | β-scaling more efficient than explicit potentials; exploits biological structure knowledge (canonicals) while allowing flexibility (especially H3) | **MEDIUM** | Number of distinct canonical combinations explored; contact score for best canonical vs. baseline; speed improvement (should be 20-30% faster than Strategy D); H3 flexibility within canonical framework |
| **E+ (Enhanced Blind Scanning)** | Discover epitope location systematically and efficiently in single diffusion pass with region emphasis | Region-specific β-scaling: emphasize different antigen patches sequentially; accumulate contact heatmap across all patches; multi-pass with adaptive refinement | Unknown epitope problems; high-throughput epitope discovery; screening multiple antigen candidates; when computation time critical (10-20x faster needed) | Single-pass with region-specific β achieves 10-20x speedup over separate runs; combines efficiency (β-scaling) with systematic scanning | **HIGH** | Epitope hotspot identification precision/recall; heatmap correlation with known epitopes; number of novel epitope regions discovered; speed comparison to Strategy E (should be 10-20x faster) |
| **V (Embedding-Space CDR3 Steering)** | More robust CDR3 exploration for novel sequences using embedding optimization instead of coordinate potentials | Extract CDR3 embedding from Pairformer; optimize embedding to match target canonical embedding; gradient descent in embedding space; final denoising refinement | Novel CDR sequences not in training set; synthetic antibody libraries; robust steering needed across wide parameter ranges | Embedding space captures biological signals more stably; robust across 2 orders of magnitude of hyperparameter variation vs. coordinate methods | **MEDIUM** | CDR3 conformation accuracy for novel sequences (RMSD to expected structure if known); stability across β parameter ranges (10-100x variation should give similar results); robustness metric (prediction quality vs. hyperparameter) |
| **W (Embedding-Based Interface Steering)** | Optimize antibody-antigen interface robustly by steering in embedding space (captures binding specificity) vs. coordinate distance | Extract CDR-antigen pairwise embeddings; optimize toward ideal "binding" embedding configuration; gradient descent in interface embedding space; back-propagation through network | Binding affinity improvement; cases with novel antigen types; when robustness more important than precision; integrating learned binding preferences | Embedding space captures binding specificity (not just proximity); natural integration with what model learned; more robust to hyperparameter variation | **MEDIUM-HIGH** | Interface binding energy prediction (if benchmark available); buried surface area in predicted complex; comparison to coordinate-based steering (should be more stable); binding affinity correlation with pLDDT |

### Phase 3 (Week 5+): Advanced Research - High Complexity, Publication-Worthy

| Name | Goal | Mechanism | Use Case | Key Insight | Complexity | Evaluation Metrics |
|------|------|-----------|----------|-------------|-----------|-------------------|
| **R (Multi-Objective Steering)** | Simultaneously optimize multiple objectives: CDR-antigen contacts, physical realism, internal antibody compatibility | Composite β-scaling: β_contact, β_cdr3, β_clash applied simultaneously; pareto-optimal ensemble selection showing trade-offs | Complex systems needing balance (don't want only contacts, must maintain structure quality); design optimization; understanding trade-offs | Multi-objective approach reveals pareto frontier of achievable solutions; more realistic than single-objective optimization | **HIGH** | Pareto frontier characterization (3D: contact vs. pLDDT vs. clash score); hypervolume of pareto front (coverage of objective space); solution diversity along frontier |
| **S (Template-Guided β-Scheduling)** | Incorporate known antibody-antigen structure knowledge as soft bias (not hard constraint) via adaptive β | Use PDB templates (human/synthetic antibodies) as reference; compute deviation from template; β(t) = f(deviation) adaptive schedule (loose when far from template, tight when close) | Humanization tasks; antibody library design; species-specific design; when template structures available; soft biasing preferred over hard constraints | Soft β-biasing more flexible than hard templates; template knowledge acts as prior not requirement; works for novel designs divergent from templates | **MEDIUM** | Structural similarity to template (RMSD); quality of predictions divergent from template (pLDDT for out-of-template conformations); contact score with/without template bias |
| **X (Experimental Constraint Integration)** | Integrate experimental data (cryo-EM maps, cross-linking, HDX-MS, SPR) into structure prediction robustly | Convert experimental constraints to embedding-space objectives; multi-modal constraint fusion with uncertainty weighting; EmbedOpt-style gradient descent in embeddings | Structure determination from electron microscopy; integrating cross-linking mass spec data; combining multiple measurement modalities; cryo-EM antibody-antigen complexes | Embedding space naturally accommodates multiple constraint types; EmbedOpt robustness handles measurement uncertainty; 3-5x accuracy improvement when experimental data available | **HIGH** | RMSD to experimental structure when available; contact agreement with experimental interface definitions; multi-constraint satisfaction score (weighted sum); speed improvement over traditional refinement |
| **U (Energy Landscape Mapping)** | Characterize complete accessible conformational space by sweeping β values and analyzing landscape | Generate 9 predictions (β ∈ [-1.0, 1.0], step 0.25); cluster by structural similarity; identify basins, plateaus, discontinuities; characterize robustness | Understanding achievable structure diversity; identifying robust vs. fragile binding modes; fundamental understanding of design space; research/exploration | Complete landscape map shows which β regions lead to similar structures; identifies robust solutions (large plateaus) vs. fragile (single points) | **LOW** | Contact score vs. β (smooth curve indicates robust; jagged indicates brittle); cluster count and size distribution; plateau detection (large same-score regions); basin depth and width |
| **Y (Hierarchical Steering)** | Combine robustness of embedding steering (early) with precision of coordinate steering (late) in smooth multi-stage approach | Stage 1 (t=T to T/2): EmbedOpt-style embedding optimization, coarse-grained steering; Stage 2 (t=T/2 to 0): coordinate-space potentials, fine-grained optimization; smooth weight transition | Any prediction needing both robustness and precision; novel/out-of-distribution systems; maximizing prediction quality across diverse inputs | Hierarchical approach leverages advantages of both methods; mimics natural structure formation; fewer total hyperparameters than two independent methods | **HIGH** | Quality metrics for each stage (stage 1 robustness, stage 2 precision); transition smoothness (no discontinuities in quality metrics); final structure quality vs. single-stage methods; computational cost (should be < 2×single method) |
| **Z (Soft Experimental Constraints)** | Handle noisy/uncertain experimental measurements probabilistically instead of as hard constraints | Model constraints as probability distributions with uncertainty bands; optimize toward likely regions (not exact values); weight objectives by measurement precision (1/σ²) | Real experimental data with noise; multi-modal measurements with different precision; Bayesian structure determination | Probabilistic approach principled and robust; avoids over-fitting to noisy measurements; natural weighting by measurement confidence | **MEDIUM** | Prediction quality vs. measurement noise level (should degrade gracefully, not cliff); constraint satisfaction spread (should be wide distribution, not point); uncertainty calibration (does model uncertainty match actual error?) |
| **P+ (Contact-Guided β-Feedback)** | Iteratively refine predictions by stratifying results based on contact quality and re-optimizing each stratum | Pass 1 (screening): fast β=-0.5 with 20% diffusion steps; stratify by contact into top/middle/bottom thirds; Pass 2 (refinement): different β values per stratum, full steps; Pass 3 (selection): choose best by contact+confidence | High-accuracy applications; when computational budget allows multiple passes; optimization that needs adaptive strategies | Stratified approach adapts to system-specific contact landscape; automatic convergence detection; avoids wasted computation on poor starting configurations | **MEDIUM** | Convergence speed (passes to plateau); contact score improvement across passes; ensemble diversity between passes; total compute cost vs. single-pass equivalents |
| **M+ (Interface β-Scaling v2)** | Optimize binding interface strength by modulating pair representations between CDR and antigen residues | Interface-specific β: z[i,j]*=(1+β) where i∈CDR, j∈antigen; two-phase: negative β for exploration, positive β for optimization; weight with predicted binding affinity | Binding affinity optimization; interface-specific design; exploring alternative binding modes; when binding strength critical | Direct control of interface properties in latent space; transparent mechanism; integration with predicted binding energy | **LOW** | Interface binding energy prediction; buried surface area (SASA) in complex; buried residue count; comparison to baseline antigen steering (should improve affinity metrics) |
| **J+ (Coupled CDR3 v2)** | Efficiently explore compatible H3/L3 combinations using beam search + asymmetric β instead of exhaustive grid | Beam search Phase 1: sample H3 with β sweep, keep top-3; Phase 2: for each H3, sample L3, score by contact+interface quality, keep top-2; joint scoring penalizes H3/L3 clashes | High-quality predictions where H3/L3 compatibility critical; reducing particle count (50-70% reduction possible); discovering cooperative binding modes | Beam search dramatically more efficient than grid; asymmetric β encodes biological dominance; avoids incompatible H3/L3 pairs | **MEDIUM** | Number of valid H3/L3 combinations (incompatible pairs detected); H3 dominance in contact score (should be ≥1.5× L3); interface quality metrics between H3 and L3 |

### Summary Statistics

- **Total Strategies**: 25 improvements/novel ideas (D through Z, excluding A-C which are baseline)
- **Phase 1 Strategies** (5): L+, K+, N+, G+ (quick wins)
- **Phase 2 Strategies** (8): Q, T, O+, D+, E+, V, W (novel contributions + research value)
- **Phase 3 Strategies** (12): R, S, X, U, Y, Z, P+, M+, J+ (advanced + experimental integration)

**Key Selection Criteria**:
- **Implementation Complexity**: LOW = <1 week, MEDIUM = 1-2 weeks, HIGH = 2-4 weeks
- **Phase Assignment**: Based on complexity, scientific novelty, and estimated impact
- **Evaluation**: Each method has specific quantitative metrics (not subjective assessments)
