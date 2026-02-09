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
