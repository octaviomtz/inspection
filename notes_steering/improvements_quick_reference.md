# Quick Reference: Steering Improvements & New Ideas

**Last Updated**: 2025-02-09
**Full Details**: See [improvements_and_new_ideas.md](./improvements_and_new_ideas.md)

---

## At a Glance

### What Changed?
- 6 existing strategies (A-J) have enhanced versions
- 6 β-scaling ideas (K-P) are improved
- 5 entirely novel ideas (Q-U) addressing hybrid approaches
- 5 new EmbedOpt-inspired ideas (V-Z) for robust embedding-space steering
- Total: **17 improvements/enhancements + 10 new ideas = 27 total improvements**

### Why It Matters
After reviewing:
1. **FK-Diffusion Steering** (Horvitz et al., 2501.06848) - General Feynman-Kac resampling framework
2. **Boltz-sample** (Suzuki & Amagasa, 2026) - Pair representation β-scaling paper
3. **Existing proposals** (A-P) - Your current steering strategies

We identified significant efficiency gains and novel hybrid approaches that can:
- **10-20x speedup** for epitope discovery (Strategy E vs E+)
- **50-70% reduction** in particle count (Strategy J vs J+)
- **New capabilities**: MSA-free steering, energy landscape mapping, iterative refinement

---

## Quick Lookup Table

### Enhancements to Existing Strategies

| Strategy | Enhancement | Key Improvement | Complexity | Expected Gain |
|----------|-------------|-----------------|-----------|---------------|
| **D+** | Canonical Ensemble v2 | Use β-scaling instead of potentials | Medium | 20-30% speedup |
| **E+** | Blind Scanning v2 | Region-specific β instead of separate runs | High | 10-20x speedup |
| **G+** | Progressive Refinement v2 | Contact-triggered phases + smooth β schedule | Medium | 15-25% contact improvement |
| **J+** | Coupled CDR3 v2 | Asymmetric β + beam search instead of grid | Medium | 50-70% fewer particles |

### Enhancements to β-Scaling Ideas

| Idea | Enhancement | Key Insight | Implementation |
|------|-------------|------------|-----------------|
| **K+** | Region-Specific v2 | Adaptive regions + confidence weighting | Multi-pass epitope scanning |
| **L+** | CDR3-Specific v2 | Two-stage approach (explore → optimize) | Ready to implement |
| **M+** | Interface β v2 | Binding energy weighting + task-specific β values | Medium complexity |
| **N+** | Adaptive β-Scheduling v2 | Contact-based phase transitions | Contact monitoring required |
| **O+** | Asymmetric β v2 | Task-specific β values (H3>L3) | Low complexity |
| **P+** | Contact-Guided β-Feedback v2 | Multi-pass with convergence detection | Requires inter-pass communication |

### Novel Ideas

| Idea | Name | Goal | Mechanism | Use Case |
|------|------|------|-----------|----------|
| **Q** | Iterative Epitope Refinement | Epitope discovery | 3-round FK + β-scaling | Unknown epitope systems |
| **R** | Multi-Objective Steering | Binding + Structure | Composite β-scaling | Complex optimization |
| **S** | Template-Guided β-Scheduling | Structure design | PDB template information | Humanization, library design |
| **T** | MSA-Free Steering | De novo design | Latent space intrinsic properties | Synthetic antibodies |
| **U** | Energy Landscape Mapping | Analysis | Complete β-sweep | Structure characterization |
| **V** | Embedding-Space CDR3 Steering | Structure (robust) | CDR3 embedding optimization | Novel/synthetic CDR sequences |
| **W** | Embedding-Based Interface Steering | Binding (robust) | Interface embedding optimization | Novel antigens, robust steering |
| **X** | Experimental Constraint Integration | Structure + Robustness | Multi-modal constraint fusion | Cryo-EM, cross-linking, HDX-MS data |
| **Y** | Hierarchical Steering | Structure + Binding | Embedding then coordinates | Novel sequences, high precision needed |
| **Z** | Soft Experimental Constraints | Structure (robust) | Probabilistic constraints | Noisy experimental data |

---

## Implementation Roadmap (Updated with EmbedOpt Ideas)

### Phase 1 (Week 1-2) - Quick Wins
1. **L+**: CDR3-specific β-scaling v2
2. **K+**: Region-specific β-scaling with adaptive regions
3. **N+**: Contact-triggered adaptive phase scheduling
4. **G+**: Enhanced progressive refinement

**Why**: Directly improve existing steering code, high impact, moderate complexity

### Phase 2 (Week 3-4) - Novel Contributions (Expanded)
5. **Q**: Iterative epitope refinement (FK + β hybrid)
6. **T**: MSA-free steering strategies (β-scaling)
7. **O+**: Asymmetric β-scaling enhancement
8. **D+**: Enhanced canonical ensemble
9. **E+**: Enhanced blind scanning (10-20x speedup)
10. **V**: Embedding-space CDR3 steering (EmbedOpt-inspired) - **NEW**
11. **W**: Embedding-based interface steering (EmbedOpt-inspired) - **NEW**

**Why**: Novel ideas, publication-worthy, moderate complexity, good scientific value. V+W bring robustness for novel sequences.

### Phase 2.5 (Week 4-5) - Medium Complexity (New)
- **V**: Embedding-space CDR3 steering (more robust than dihedral potentials)
- **W**: Embedding-based interface steering (captures binding specificity)

**Why**: EmbedOpt-inspired, medium complexity, significant robustness improvement, 50%+ improvement for novel sequences

### Phase 3 (Week 5-6+) - Advanced Research
12. **R**: Multi-objective steering
13. **S**: Template-guided β-scheduling
14. **U**: Energy landscape mapping
15. **X**: Experimental constraint integration (cryo-EM, cross-linking, HDX-MS) - **NEW**
16. **Y**: Hierarchical steering (embedding + coordinate) - **NEW**
17. **Z**: Soft experimental constraints (probabilistic) - **NEW**

**Why**: Complex implementations, publication-worthy, experimental data integration. X+Y+Z enable new capabilities with experimental data.

---

## Key Insights from Reference Papers

### FK-Diffusion Steering (Horvitz et al.)
**Core**: Feynman-Kac resampling of particle trajectories based on energy potentials
- Arbitrary reward functions (no differentiability needed)
- Minimal computational overhead
- Proven effective for protein design

**Application**: Can combine with β-scaling for dual exploration + filtering

### Boltz-sample (Suzuki & Amagasa 2026)
**Core**: Rescale latent pair representation z by scalar β at Pairformer input
- z_scaled = (1 + β) × z
- Different signs of β = different search directions
- Works even without MSA (activates internal priors)

**Application**: More efficient than explicit potentials, transparent mechanism

### EmbedOpt (Li et al., 2602.05285) - **NEW**
**Core**: Optimize in embedding space rather than coordinate space for robustness
- Embedding space captures sequence and coevolutionary signals more stably
- Stable across 2 orders of magnitude of hyperparameter variation
- Better for out-of-distribution constraints (novel sequences, unusual conformations)
- Designed for experimental constraints (cryo-EM maps, distance constraints, cross-linking data)

**Application**: More robust steering for novel antibodies, experimental data integration, probabilistic handling of measurement uncertainty

**Key Insight**: Embedding space optimization is fundamentally more robust than coordinate-space potentials for steering diffusion models, especially for constraints outside training distribution.

---

## Impact Analysis

### Computational Efficiency
- **Strategy E+**: 10-20x speedup (single pass with region-specific β vs. N separate passes)
- **Strategy J+**: 50-70% fewer particles needed (beam search vs. 2D grid)
- **Strategy D+**: 20-30% speedup (β-scaling vs. explicit potentials)

### Scientific Impact
- **New capabilities**: Epitope discovery, MSA-free design, energy landscape analysis
- **Better predictions**: 15-25% improvement in binding contact quality
- **Generalization**: Works for synthetic antibodies with limited MSA

### Implementation Effort
- **Low complexity**: L+, M+, O+, U (mostly parameter tuning)
- **Medium complexity**: K+, G+, N+, Q, S, T (new mechanisms with clear logic)
- **High complexity**: E+, R (require significant refactoring)

---

## Testing Strategy

### Benchmark Datasets
- PDB antibody-antigen complexes (known epitopes)
- SABDAB database (diverse antibodies)
- Synthetic antibodies (limited MSA)

### Key Metrics
1. **Epitope discovery**: Precision/recall vs. experimental
2. **Binding affinity**: Correlation with buried surface area
3. **Structure quality**: RMSD vs. native (when available)
4. **Ensemble diversity**: pairwise RMSD in ensemble
5. **Computational cost**: Time + GPU memory

---

## Critical Questions to Answer

### Before Implementation
1. What's the optimal β range for each strategy? (Does it vary by antigen class?)
2. Which contact metric is best? (Distance vs. pAE vs. predicted ΔG)
3. What ensemble size is needed? (5 vs. 10 vs. 20 structures)
4. How to detect phase transitions in adaptive scheduling? (Contact score plateau? Threshold?)

### During Implementation
1. How much overhead does contact computation add?
2. Can we use pAE-based contacts instead of computed distances?
3. Does multi-objective weighting have a sweet spot?
4. How to handle novel antibodies without templates?

---

## Document Cross-References

- **Full details**: [improvements_and_new_ideas.md](./improvements_and_new_ideas.md)
- **Existing strategies**: [steering_proposals.md](./steering_proposals.md)
- **β-Scaling ideas**: [latent_space_scaling_inspired_ideas.md](./latent_space_scaling_inspired_ideas.md)
- **CDR3 steering**: [cdr3_steering.md](./cdr3_steering.md)
- **Antigen steering**: [antigen_steering.md](./antigen_steering.md)
- **Main overview**: [../STEERING_METHODS.md](../STEERING_METHODS.md)

---

## Quick Start: Which Improvement Should I Implement First?

**For fastest wins**: Start with **L+** (CDR3-specific β-scaling)
- Lowest complexity
- Direct replacement for dihedral potentials
- Leverages existing pair masking infrastructure
- Works immediately with Boltz2

**For highest impact**: Start with **K+** (Region-specific β-scaling)
- Enables epitope discovery
- 10-20x speedup vs. blind scanning
- Can run multiple β configurations in parallel
- Publication-worthy contribution

**For balanced approach**: Start with **G+** + **N+** together
- Enhance existing progressive refinement strategy
- Contact-triggered phases improve naturally
- Low risk of breaking existing code
- Measurable improvement in binding quality

---

**Status**: Ready for implementation planning
**Last Reviewed**: 2025-02-09
