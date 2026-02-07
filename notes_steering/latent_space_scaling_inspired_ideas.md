# Latent Space Scaling-Inspired Steering Ideas

**Reference**: Suzuki & Amagasa (2026) - "Steering Conformational Sampling in Boltz-2 via Pair Representation Scaling" (Boltz-sample)

**Date**: 2025-02-06

**Status**: Ideation Phase - Novel steering strategies inspired by pair representation β-scaling

---

## Summary of β-Scaling Paper

The Suzuki & Amagasa paper introduces a simple yet powerful steering mechanism:

**Core Idea**: Rescale the latent pair representation `z` by a scalar β at the Pairformer input:
```
z_scaled = (1 + β) * z
```

**Key Findings**:
- **Efficiency**: Negligible computational overhead (unlike MSA clustering)
- **Model Priors**: Works even without MSA - activates internal structural priors
- **Directional Steering**: Sign of β defines distinct search directions (β < 0 vs β > 0)
- **Mechanistic**: Acts like diffusion guidance on signal-to-noise ratio of pairwise couplings
- **Interpretability**: Transparent control parameter, not black-box input perturbation

---

## Novel Ideas for Antibody-Antigen Binding

### Idea K: Region-Specific β-Scaling for Epitope Discovery

**Goal**: Identify epitope without prior knowledge

**Motivation**: The paper shows that β defines search directions in latent space. For antibody-antigen binding, different antigen surface patches should have different "preference" for CDR contact based on the model's learned priors.

**Approach**:
1. **Antigen Partitioning**: Divide antigen surface into overlapping spherical regions (similar to Strategy E from proposals)
2. **Per-Region β Scaling**: For each region R:
   - Apply strong negative β to residues outside R (relax constraints elsewhere)
   - Apply strong positive β to residues in R (enforce constraints)
   - This creates region-specific "emphasis"
3. **Latent Space Implementation**: Instead of explicit spatial constraints, modulate the pair representation for (CDR, antigen_region) pairs
4. **Epitope Heatmap**: Across multiple β-scans, collect which antigen regions naturally contact CDRs when emphasized

**Expected Benefits**:
- Discovers epitope using model's internalized priors (no external predictor needed)
- Systematic scanning via latent space rather than explicit spatial constraints
- Works with or without MSA
- Generates multiple binding modes per candidate epitope region

**Implementation Complexity**: Medium (requires masking scheme for pair-wise scaling)

**Target Goal**: **Epitope Identification** (secondary benefit: binding affinity)

---

### Idea L: CDR3-Specific β-Scaling for Conformation Exploration

**Goal**: Explore CDR3 conformations without explicit dihedral potentials

**Motivation**: The paper's sequence-only results show that β-scaling activates structural priors even without evolutionary signals. CDR3 loops have strong structural priors in the Boltz2 training distribution.

**Approach**:
1. **CDR3-Specific Scaling**: Apply β only to pair features involving CDR3 residues:
   - z_scaled[i,j] = (1 + β) * z[i,j] if (i or j) in CDR3 region
   - z_scaled[i,j] = z[i,j] otherwise
2. **Ensemble Generation**: Sweep β from -0.75 to +0.75 to generate CDR3 conformational ensemble
3. **No Explicit Potentials**: Unlike current CDR3 steering, this requires no dihedral angle constraints
4. **Contact-Aware Selection**: Use simple confidence-based selection (like sign-based selection in paper) to identify conformations that contact antigen well

**Expected Benefits**:
- Simpler than explicit potential-based CDR3 steering (no need to define target angles)
- Leverages model's internalized loop priors
- Works better with limited MSA (since β works without MSA)
- Generates diverse CDR3 conformations exploring latent space

**Implementation Complexity**: Low-Medium (local modification to pair scaling)

**Target Goal**: **Structure Prediction Improvement** + **CDR3 Exploration**

---

### Idea M: Inter-Chain β-Scaling for Antibody-Antigen Interface Optimization

**Goal**: Improve antigen-antibody binding affinity by modulating interface strength

**Motivation**: The paper shows β controls the "signal-to-noise ratio" of pairwise couplings. At the antibody-antigen interface, we could use positive β to strengthen the learned couplings and negative β to explore alternative binding modes.

**Approach**:
1. **Interface-Specific Scaling**: Apply β selectively to pair features between CDR and antigen residues:
   - z_scaled[i,j] = (1 + β) * z[i,j] if (i in CDR and j in antigen) or vice versa
   - z_scaled[i,j] = z[i,j] for intra-chain pairs
2. **Two-Phase Sampling**:
   - **Phase A** (β < 0): Relaxed interface allows exploration of multiple binding geometries
   - **Phase B** (β > 0): Tightened interface optimizes for strong contacts
3. **Ensemble Selection**: Use contact score + pLDDT to select models from both phases

**Expected Benefits**:
- Direct modulation of interface strength in latent space
- Discovers both loose-binding (exploratory) and tight-binding (optimized) modes
- More efficient than particle-based resampling (no importance weighting needed)
- Transparent control of interface properties

**Implementation Complexity**: Low (pair-wise masking based on chain IDs)

**Target Goal**: **Binding Affinity Improvement**

---

### Idea N: Adaptive β-Scheduling During Diffusion

**Goal**: Progressively guide structure toward better binding throughout the diffusion process

**Motivation**: The paper applies β uniformly at Pairformer input of each recycle step. We could adaptively change β as diffusion progresses, combining benefits of early exploration with late-stage optimization.

**Approach**:
1. **Time-Dependent β Schedule**: β(t) = f(t) where t is diffusion timestep
   - Early diffusion (t large): β(t) = -0.5 (exploratory, relax constraints)
   - Mid diffusion (t medium): β(t) = 0 (neutral, model's default behavior)
   - Late diffusion (t small): β(t) = +0.5 (optimizing, enforce interactions)

2. **Contact-Aware Scheduling**: Make schedule adaptive to current contact score
   - If contacts are poor: use negative β to explore
   - If contacts are good: use positive β to refine

3. **Multi-Phase Steering**: Similar to Strategy G (Progressive CDR Refinement) but using latent space scaling instead of explicit potentials

**Expected Benefits**:
- Natural progression from exploration to exploitation
- No explicit energy function needed (uses learned priors)
- Computationally efficient (just modify one scalar per step)
- Balances diversity and accuracy

**Implementation Complexity**: Low-Medium (requires diffusion step tracking in pair scaling)

**Target Goal**: **Binding Affinity Improvement** + **Structure Prediction**

---

### Idea O: Differential β-Scaling Between Heavy and Light Chain CDRs

**Goal**: Improve predictions by recognizing asymmetric roles of CDR-H3 vs CDR-L3

**Motivation**: In antibody biology, heavy chain CDR3 typically dominates binding (makes ~60% of contacts), while light chain CDR3 provides supporting contacts. The paper shows β can separately control different coupling strengths.

**Approach**:
1. **Asymmetric Scaling**:
   - Apply β_H to pair features involving CDR-H3
   - Apply β_L (different value) to pair features involving CDR-L3
   - Use β_H > β_L to reflect CDR-H3 dominance

2. **Biologically-Informed Defaults**:
   - β_H = +0.3 to +0.5 (emphasize heavy chain contacts)
   - β_L = -0.1 to +0.1 (lighter emphasis on light chain)

3. **Ensemble Generation**: Sweep both β_H and β_L to explore the 2D landscape of heavy-light chain combinations

**Expected Benefits**:
- Respects known antibody biology (H3 > L3 in binding)
- Reduces wasted particles on incompatible H3/L3 pairs (compared to joint sampling)
- More efficient than Strategy J (Coupled CDR3) since no explicit sampling grid needed
- Natural way to encode domain knowledge in latent space

**Implementation Complexity**: Medium (requires chain and CDR region tracking)

**Target Goal**: **Binding Affinity Improvement** + **Structure Prediction**

---

### Idea P: Contact-Guided β-Scaling Feedback Loop

**Goal**: Iteratively improve antibody-antigen predictions using contact information

**Motivation**: The paper shows sign of β defines search directions. We could use contact information to dynamically select which search direction (β sign) is beneficial at each iteration.

**Approach**:
1. **Multi-Pass Prediction**:
   - Pass 1: Generate predictions with β = -0.5 (exploratory)
   - Compute CDR-antigen contacts for each prediction
   - Rank predictions by contact quality

2. **Adaptive Second Pass**:
   - For top 30% (good contacts): use β = +0.5 (refine)
   - For bottom 30% (poor contacts): use β = -0.75 (radical exploration)
   - For middle 40%: use β = 0 (neutral refinement)

3. **Final Ensemble**: Combine results, selecting by contact score + pLDDT

**Expected Benefits**:
- Data-driven selection of search direction
- Automatic adaptation to each system (no manual β tuning)
- Leverages model's ability to predict contacts
- Reduces wasted sampling on poor starting conformations

**Implementation Complexity**: Medium (requires contact computation between predictions)

**Target Goal**: **Binding Affinity Improvement**

---

## Comparison to Existing Strategies

| ID | Name | β-Inspired | Mechanism | Latent Space | Simplicity | Epitope | Affinity | Structure |
|----|----|-----------|-----------|--------------|-----------|---------|----------|-----------|
| K | Region-Specific β-Scaling | ✓ | Per-region emphasis | Yes | Medium | ✓✓✓ | ✓✓ | ✓ |
| L | CDR3-Specific β-Scaling | ✓ | CDR3 conformation space | Yes | Low | ✓ | ✓ | ✓✓✓ |
| M | Inter-Chain Interface β | ✓ | Interface modulation | Yes | Low | ✓✓ | ✓✓✓ | ✓✓ |
| N | Adaptive β-Scheduling | ✓ | Time-dependent scaling | Yes | Low | ✓ | ✓✓ | ✓✓✓ |
| O | Asymmetric CDR β-Scaling | ✓ | Heavy>Light weighting | Yes | Low | ✓ | ✓✓✓ | ✓✓ |
| P | Contact-Guided β-Feedback | ✓ | Iterative adaptation | Yes | Medium | ✓✓ | ✓✓✓ | ✓✓ |

---

## Key Advantages of β-Scaling Approach

1. **Efficiency**: Negligible computational overhead (unlike potentials-based steering)
2. **No Retraining**: Works with existing Boltz2 weights
3. **Transparency**: Clear mechanistic interpretation (signal-to-noise ratio of pairwise couplings)
4. **MSA-Agnostic**: Works with or without evolutionary information
5. **Latent Space Native**: Direct control of model's internal representations
6. **Ensemble-Friendly**: Simple confidence-based selection without clustering

---

## Implementation Priority

### Phase 1 (Immediate): Simplest to Implement
1. **Idea L** (CDR3-Specific β-Scaling): Minimal code changes, leverages existing pair masking
2. **Idea M** (Inter-Chain Interface β): Simple interface identification, direct affinity improvement

### Phase 2 (Near-term): Medium Complexity
3. **Idea N** (Adaptive β-Scheduling): Requires diffusion step tracking but straightforward logic
4. **Idea O** (Asymmetric CDR β-Scaling): Adds complexity but biologically motivated

### Phase 3 (Future): More Involved
5. **Idea K** (Region-Specific Epitope Discovery): Requires antigen partitioning logic
6. **Idea P** (Contact-Guided Feedback): Requires inter-prediction communication

---

## Critical Questions

1. **Computational Efficiency**: How much overhead does masking (for region/chain/CDR-specific scaling) add?
2. **β-Range Optimality**: For each idea, what is the optimal β range? (-0.75 to +0.75 works for general case but may differ for specific applications)
3. **Combination Effects**: Can multiple β-scalings be applied simultaneously (e.g., CDR3-specific + inter-chain)? Do they interfere?
4. **MSA Dependency**: How does performance degrade when MSA is limited or absent? (Paper shows improvement but antibodies are well-represented in databases)
5. **Generalization**: Do β values learned on training set transfer to unseen antibody-antigen pairs?

---

## Research Contributions

These ideas represent novel applications of latent space steering to antibody-antigen problems:

- **K**: First systematic latent-space epitope discovery method
- **L**: Alternative to explicit dihedral potentials for CDR3 sampling
- **M**: Direct latent-space control of binding interface
- **N**: Novel diffusion scheduling strategy using latent modulation
- **O**: Biologically-informed latent space asymmetry
- **P**: Feedback-driven iterative improvement

All maintain the key advantages of the β-scaling approach: simplicity, efficiency, and transparency.

---

## References

- Suzuki, S. & Amagasa, T. (2026). Steering Conformational Sampling in Boltz-2 via Pair Representation Scaling. *bioRxiv* 2026.01.23.701250
- Boltz Sample Repository: https://github.com/suzuki-2001/boltz-sample
- Related Work: Classifier-free diffusion guidance (Ho & Salimans, 2022)
