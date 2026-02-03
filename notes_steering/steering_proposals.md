# Steering Strategy Proposals

**Status**: Brainstorming / Design Phase
**Focus**: Better antibody-antigen binding prediction, especially with unknown epitope location
**Key Challenge**: Exploring CDR conformations while optimizing antigen orientation

---

## Executive Summary

This document outlines potential steering strategies to improve antibody-antigen binding prediction when the epitope location is unknown. The strategies range from simple extensions of existing ideas (A, B, C) to novel ensemble-based approaches that explore the conformational space of CDR loops while simultaneously optimizing antigen orientation.

The common theme: **enable exploration of CDR conformations while guiding antigen orientation through contact scoring**.

---

## Existing Strategies (Reference)

### Strategy A: Multi-Particle Antigen Rotation Potential
- Generate multiple rotated/translated copies of the antigen as particles
- Score each by proximity of antigen residues to CDR residues
- Use Feynman-Kac resampling to select best orientations
- **Pros**: Leverages existing particle framework; direct exploration
- **Cons**: Limited CDR flexibility; purely rotation-based

### Strategy B: Gradient-Based Antigen Orientation Steering
- Define energy function based on CDR-antigen distances
- Compute gradients that push antigen toward optimal orientation
- Apply rotational/translational updates via guidance mechanism
- **Pros**: Direct, continuous optimization; efficient
- **Cons**: May get stuck in local minima; no CDR exploration

### Strategy C: Coarse-Grained Search + Fine-Grained Steering
- Early diffusion: Sample multiple antigen orientations via center-of-mass distances
- Mid-to-late diffusion: Switch to contact-based fine-grained steering
- **Pros**: Balances exploration and exploitation; two-phase approach
- **Cons**: Requires careful tuning of phase transition

---

## Novel Strategies

### Strategy D: CDR Canonical Ensemble Steering

**Motivation**: Research shows that CDR1, CDR2, and CDR-L3 adopt canonical (fixed) backbone conformations, while only CDR-H3 is highly variable. Exploit this structure.

**Approach**:
1. **Canonical Potential Library**: Pre-compute energy potentials for known canonical CDR conformations
   - For CDR1 and CDR2: 5-10 dominant canonical structures per residue
   - For CDR-L3: ~3-5 canonical cluster centers
   - For CDR-H3: Allow full flexibility (no canonical constraint)

2. **Multi-Canonical Particles**: Run diffusion with particles, each initialized to different canonical conformations
   - Particle 1: Canonical CDR1a + CDR2a + CDR-L3a
   - Particle 2: Canonical CDR1a + CDR2a + CDR-L3b
   - Particle 3: Canonical CDR1b + CDR2b + CDR-L3a
   - etc. (combinatorial sweep)

3. **Joint Antigen-CDR Optimization**:
   - At each step: CDR3 backbone can move freely (for conformation search)
   - Antigen orientation optimized to contact whichever canonical CDR conformation is active
   - Heavy chain CDR3 and light chain CDR3 can explore independently

4. **Resampling**: Particles scored by:
   - CDR-antigen contact quality
   - Preservation of canonical backbone angles (for CDR1/2/L3)
   - Weak particles discarded, good ones replicated

**Expected Benefits**:
- Systematic exploration of CDR conformation space
- Constrain mostly-rigid CDRs, allow flexible CDR-H3
- Biologically plausible (exploits known structure knowledge)
- Can discover which canonical + CDR-H3 combination works best

**Implementation Complexity**: Medium (requires canonical structure database)

---

### Strategy E: Unknown Epitope Explorer (Blind Antigen Scanning)

**Motivation**: When epitope location is unknown, systematically scan the antigen surface to find which regions CDRs naturally contact.

**Approach**:
1. **Antigen Surface Partitioning**: Divide antigen into N overlapping spherical regions
   - Center at different points on antigen surface
   - Radius ~20-30 Å (covers ~2-3 CDRs)

2. **Parallel Emphasis Steering**: For each region:
   - Run diffusion with strong guidance to position that region toward CDR loops
   - Other regions allowed to position freely
   - Use particles to explore multiple modes per region

3. **Contact Heat Map Generation**:
   - Collect which antigen residues contact CDR residues across all runs
   - Build 3D heat map showing "epitope likelihood" per residue
   - Identify clusters of high-contact residues

4. **Epitope Refinement Phase** (optional second pass):
   - Run guided prediction with identified epitope region emphasized
   - Generate final high-quality complex structures

**Expected Benefits**:
- Discovers epitope without prior knowledge
- Generates multiple binding modes per candidate epitope
- Produces contact maps showing interaction patterns
- Can identify cryptic or multiple epitopes

**Implementation Complexity**: Medium-High (requires systematic scanning + contact tracking)

**Computational Cost**: High (N × K particles × diffusion steps)

---

### Strategy F: Conformational Ensemble Binding Mode Discovery

**Motivation**: Instead of finding a single "best" structure, generate an ensemble of equally-good binding conformations that differ in CDR backbone geometry.

**Approach**:
1. **Ensemble Particles**: Run diffusion with multiple particles
   - Each initialized with different random seed (different CDR3 backbone exploration)
   - NOT resampled based on contact (preserve diversity)
   - All particles allowed to explore full conformational space

2. **Weak Guidance**: Apply loose CDR-antigen contact guidance
   - Not forcing optimal contacts (like Strategy B)
   - Just "pulling" toward reasonable binding geometry
   - Soft penalty rather than hard constraint

3. **Diversity Reward** (optional):
   - Add bonus energy term for particles whose CDR conformations differ significantly
   - Encourages exploration of different binding modes
   - Prevents convergence to single conformation

4. **Output**: Ensemble of 5-10 structures showing:
   - Different CDR3 backbone angles
   - Different antigen binding modes
   - All satisfying "reasonable binding" criterion

**Expected Benefits**:
- Captures binding plasticity (how flexible antibody is)
- Shows multiple functional binding modes
- Useful for understanding antibody promiscuity
- More biologically realistic (antibodies aren't rigid)

**Implementation Complexity**: Low-Medium (mostly parameter tuning)

---

### Strategy G: Progressive CDR Refinement (Multi-Phase Steering)

**Motivation**: CDR formation happens in stages; use this to progressively refine structure.

**Approach**:
1. **Phase 1 - Global Search** (t=T to t=T/2, first half of diffusion):
   - Loose CDR3 guidance (allow exploration)
   - Antigen position entirely flexible (find any reasonable orientation)
   - Particles explore broadly without heavy penalties
   - Goal: Sample diverse CDR3 conformations with plausible antigen orientations

2. **Phase 2 - Contact Identification** (t=T/2 to t=T/4):
   - Tighten contact-based guidance
   - Switch on "epitope finding" - identify which antigen regions contacted by each particle
   - Begin steering antigen toward high-contact regions

3. **Phase 3 - Fine-Tuning** (t=T/4 to t=0):
   - Strong CDR-antigen contact optimization
   - Lock in antigen orientation (prevent drift)
   - Refine CDR side-chain packing and backbone angles
   - Optimize hydrogen bonding and van der Waals interactions

**Schedule Implementation**:
```
Phase 1 guidance_weight: 0.1 (loose)    → Phase 2: 0.5 (medium)    → Phase 3: 1.0 (tight)
Phase 1 antigen_penalty: 0.0 (free)    → Phase 2: 0.3 (soft)      → Phase 3: 1.0 (locked)
Phase 1 cdr3_temperature: 100°K (hot)  → Phase 2: 50°K (cool)     → Phase 3: 10°K (frozen)
```

**Expected Benefits**:
- Balances exploration (early) vs. exploitation (late)
- Natural progression mimics structure formation
- Systematic refinement at appropriate timescales
- Parameter-tuned per phase for optimal behavior

**Implementation Complexity**: Medium (requires phase scheduling)

---

### Strategy H: Entropy-Maximizing Ensemble Generation

**Motivation**: Instead of optimizing for maximum contact, generate ensemble of conformations with maximum contact diversity while maintaining binding quality.

**Approach**:
1. **Contact Quality Threshold**: Define minimum acceptable CDR-antigen contact score (e.g., top 20% range)
   - Not trying for absolute best (avoiding overfitting)
   - Just "good enough" binding

2. **Diversity Reward Potential**:
   - Compute pairwise distances between particles in CDR3 backbone angle space
   - Add energy bonus for maximum diversity: `E_diversity = -sum(d_ij)` where d_ij = distance between particles i and j
   - Balance with contact energy: `E_total = E_contact + λ × E_diversity`

3. **Ensemble Output**: Collection of conformations with:
   - Similar binding quality (contact scores within threshold)
   - Maximally different CDR3 geometries
   - Explores "accessible" conformational space

**Expected Benefits**:
- Shows conformational plasticity under binding
- Captures "ensemble binding" phenomenon
- Useful for epitope mapping and antibody engineering
- More robust to mutations than single-structure design

**Implementation Complexity**: Medium (requires diversity metric definition)

---

### Strategy I: Epitope Surface Scoring + CDR3 Targeting

**Motivation**: Explicitly score antigen surface patches and steer CDR3 to high-scoring patches.

**Approach**:
1. **Pre-compute Epitope Propensity**:
   - For each antigen residue: compute propensity to be in antibody-antigen interface
   - Use structure geometry (surface accessibility, cavity depth)
   - Or use existing epitope prediction tools (ESM-based methods mentioned in recent literature)

2. **Heavy-Chain CDR3 Focusing**:
   - CDR-H3 typically makes most contacts
   - Explicitly steer antigen such that high-propensity regions face CDR-H3
   - Use orientation potential: `E = -sum(propensity_i × contact_i)`

3. **Light-Chain CDR3 Positioning**:
   - Secondary steering for CDR-L3
   - Position to support CDR-H3 contacts (cooperative binding)
   - Lower weight than CDR-H3 steering

4. **Validation Phase** (optional):
   - Run without epitope propensity information to avoid bias
   - Compare discovered contacts to propensity predictions
   - Validate that model found the predicted epitope

**Expected Benefits**:
- Uses modern ML predictions (epitope propensity models)
- Focuses on biologically likely regions
- Heavy chain (CDR-H3) gets priority (biologically correct)
- Can validate epitope predictions experimentally

**Implementation Complexity**: Medium-High (requires epitope scoring function)

**Note**: Could integrate with external epitope prediction tools (e.g., ParaSurf, EPP from recent literature)

---

### Strategy J: Coupled CDR3 Heavy-Light Chain Sampling

**Motivation**: Heavy chain CDR3 and light chain CDR3 work cooperatively; should explore jointly, not independently.

**Approach**:
1. **Coupled Sampling Space**:
   - Define 2D subspace: (CDR-H3 backbone angles, CDR-L3 backbone angles)
   - Typical: 30-50 canonical/representative points in this space
   - Create particles at grid points: (H3_canonical_a, L3_canonical_b), (H3_canonical_a, L3_canonical_c), etc.

2. **Joint Scoring**:
   - Score not just CDR-antigen contact individually
   - Score quality of CDR-H3 + CDR-L3 interface (complementarity)
   - Penalize cases where H3 and L3 clash with each other

3. **Particles Explore Together**:
   - Antigen position optimized for both CDR3s simultaneously
   - Good particles: good CDR3-antigen + good CDR3-CDR3 interface
   - Bad particles discarded (resampling)

4. **Output**: Discovers cooperative binding modes where H3 and L3 work together

**Expected Benefits**:
- Captures cooperative binding (H3 + L3 better together)
- Prevents unphysical configurations (H3/L3 clashes)
- More realistic antibody binding geometry
- Fewer "wasted" particles on incompatible H3/L3 pairs

**Implementation Complexity**: Medium (requires joint CDR3 sampling space)

---

## Comparison Matrix

| Strategy | CDR Exploration | Antigen Orientation | Epitope Discovery | Complexity | Compute Cost |
|----------|-----------------|---------------------|-------------------|-----------|--------------|
| A (Multi-particle rotation) | No | Yes (rotation) | Limited | Low | Medium |
| B (Gradient steering) | No | Yes (direct) | No | Low | Low |
| C (Coarse→Fine) | No | Yes (two-phase) | Partial | Medium | Medium |
| **D (Canonical ensemble)** | Yes (limited) | Yes | Partial | Medium | High |
| **E (Blind scanning)** | Yes | Yes | **Yes** | High | Very High |
| **F (Ensemble modes)** | Yes | Yes | Partial | Low | Medium |
| **G (Multi-phase)** | Yes | Yes | Partial | Medium | Medium |
| **H (Entropy-max)** | Yes | Yes | No | Medium | Medium |
| **I (Epitope scoring)** | Yes | Yes | **Yes** (with input) | Medium | Medium |
| **J (Coupled CDR3)** | Yes | Yes | Partial | Medium | High |

---

## Recommended Implementation Roadmap

### Phase 1 (Immediate): Low-Risk Extensions
1. **Strategy D (Canonical Ensemble)**: Builds on existing CDR3 steering
   - Adds combinatorial canonical CDR exploration
   - Medium complexity, good payoff
   - Validates canonical structure hypothesis

2. **Strategy G (Multi-Phase)**: Simple schedule-based approach
   - Reuse existing guidance mechanisms
   - Just add scheduling to guidance weights
   - No new physics needed

### Phase 2 (Near-term): Novel Contributions
3. **Strategy E (Blind Epitope Scanning)**: Original contribution
   - Directly addresses "unknown epitope" problem
   - High visibility for publication
   - High compute cost but parallelizable

4. **Strategy I (Epitope Scoring)**: Integration with modern ML
   - Leverage recent epitope prediction advances (literature references below)
   - Bridges structure prediction and immunoinformatics
   - Medium complexity

### Phase 3 (Advanced): Ensemble-Based Methods
5. **Strategy F (Ensemble Modes)**: Simplest ensemble approach
   - Low complexity, good scientific value
   - Shows conformational plasticity
   - Good for understanding binding mechanisms

6. **Strategy J (Coupled CDR3)**: Sophisticated sampling
   - High-quality predictions
   - Captures antibody biology
   - More compute-intensive

---

## Recent Literature References

The following recent advances (2025) should inform implementation:

- **Epitope Prediction**: [Improved Graph-based Antibody-aware Epitope Prediction with Protein Language Model-based Embeddings](https://www.biorxiv.org/content/10.1101/2025.02.12.637989v1) and [ParaDeep](https://www.frontiersin.org/journals/bioinformatics/articles/10.3389/fbinf.2025.1684042/full) show state-of-the-art ML epitope scoring

- **CDR Conformations**: [Canonical structure research](https://www.tandfonline.com/doi/full/10.1080/19420862.2020.1744328) shows CDR loops as ensembles, with CDR-H3 being most variable

- **Guided Diffusion**: [Monte Carlo Tree Diffusion with Multiple Experts](https://arxiv.org/html/2509.15796) demonstrates exploration-exploitation tradeoffs in diffusion sampling

- **De Novo Design**: [De novo design of epitope-specific antibodies via structure-driven workflow](https://www.nature.com/articles/s41467-025-67361-9) shows integration of structure prediction with design

---

## Key Open Questions

1. **CDR3 Flexibility**: How flexible should CDR-H3 be during guided sampling? (Hard vs. soft constraints)

2. **Particle Count**: Optimal number of particles for each strategy? (3 vs. 5 vs. 10 vs. 20?)

3. **Phase Transitions**: When should diffusion phases transition? (Based on timestep or energy landscape?)

4. **Epitope Definition**: What counts as "found" epitope? (Residue contact threshold? Contact surface area?)

5. **Ensemble Size**: How many conformations needed to characterize binding ensemble? (5 vs. 10 vs. 20?)

---

## Next Steps

1. **Choose 1-2 strategies** from Phase 1 for implementation
2. **Design computational experiments** to validate approach
3. **Benchmark** against existing methods (B, C)
4. **Test** on known antibody-antigen complexes with available structures
5. **Publish** findings (novel steering strategies)
6. **Iterate** to Phase 2 strategies based on results

---

## Document Metadata

- **Created**: 2025-02-03
- **Status**: Design Phase / Brainstorming
- **Contributors**: Claude (AI) + Octavio (User)
- **Related Files**:
  - [[STEERING_METHODS]] (Overview)
  - [[cdr3_steering]] (Existing)
  - [[antigen_steering]] (Existing)
