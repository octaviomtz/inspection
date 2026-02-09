# Improvements to Steering Proposals and New Ideas

**Date**: 2025-02-09
**Status**: Analysis Phase
**Scope**: Revisited steering strategies based on FK-Diffusion Steering (Horvitz et al., 2501.06848) and Boltz-sample β-scaling (Suzuki & Amagasa, 2026)

---

## Executive Summary

After reviewing the FK-Diffusion Steering framework and the Boltz-sample β-scaling paper, we identify:

1. **Improvements to existing strategies** (A-J) in steering_proposals.md
2. **Enhancements to β-scaling ideas** (K-P) in latent_space_scaling_inspired_ideas.md
3. **Novel hybrid approaches** combining FK particle resampling with β-scaling
4. **MSA-free steering** strategies for limited evolutionary information
5. **Multi-objective optimization** strategies for simultaneous binding + epitope goals

---

## Part 1: Improvements to Existing Strategies (A-J)

### Strategy D (Canonical Ensemble) - Enhanced Version

**Current Implementation**: Pre-compute energy potentials for canonical CDR conformations, run particles at different canonical combinations

**Improvements**:
1. **β-Scaling Alternative**: Instead of explicit potentials, use CDR3-specific β-scaling (Idea L) to explore canonical space
   - More efficient than potential-based steering
   - Leverages Boltz2's internalized structural priors
   - Works better with limited MSA

2. **Adaptive Canonical Selection**: Use contact-aware selection to identify which canonical combination(s) work best
   - Pre-filter canonical combinations based on preliminary antigen proximity
   - Use FK resampling weights instead of binary keep/discard

3. **CDR-H3 Sub-Canonical Flexibility**:
   - Allow CDR-H3 to explore sub-canonical variations around canonical centers
   - Use soft β-scaling (smaller β range) rather than hard constraints
   - Enables exploration of intermediate conformations

4. **Implementation Notes**:
   - Use masking in pair representation to selectively apply β to CDR1/2/L3 regions
   - Apply weaker β (±0.3) to preserve flexibility
   - Combine with FK particle resampling for adaptive exploration

**Expected Improvement**:
- 20-30% computational savings vs. explicit potential approach
- Better exploration of biologically relevant canonical space
- Fewer wasted particles on incompatible conformations

---

### Strategy E (Blind Antigen Scanning) - Enhanced Version

**Current Implementation**: Partition antigen surface, run separate steering runs for each region, build contact heatmap

**Improvements**:
1. **Latent Space Epitope Discovery (Idea K)**: Replace explicit spatial constraints with region-specific β-scaling
   - More efficient than running N separate diffusion passes
   - Use single diffusion pass with masked β-scaling per region
   - Can scan antigen surface in seconds rather than hours

2. **Contact Heat Map with Confidence**:
   - Weight contact contributions by pLDDT scores
   - Use PAE to identify high-confidence interfaces
   - Generate probabilistic epitope map rather than binary contacts

3. **Dynamic Region Definition**:
   - Instead of fixed spherical regions, use surface curvature and accessibility to define epitope candidates
   - Adaptive region size based on antigen geometry
   - Focus scanning on biologically relevant epitope regions

4. **Iterative Refinement (New Idea Q below)**:
   - Use first-pass scanning to identify candidate epitope regions
   - Run high-confidence secondary pass with region-specific guidance
   - Validate against known epitopes when available

5. **Implementation Notes**:
   - Single FK particle sweep with region-specific β values
   - Parallel processing of different β configurations
   - Contact accumulation across ensemble

**Expected Improvement**:
- 10-20x speedup vs. Strategy E (single pass vs. N passes)
- More robust epitope identification with confidence weighting
- Discovers cryptic or secondary epitopes better

---

### Strategy G (Progressive CDR Refinement) - Enhanced Version

**Current Implementation**: Three-phase approach with manually tuned guidance weights at each phase

**Improvements**:
1. **Adaptive Phase Scheduling (Idea N)**:
   - Instead of fixed diffusion timestep thresholds, use contact-aware triggers
   - Transition to next phase when contact score plateaus or threshold met
   - More natural progression than fixed time-based scheduling

2. **Continuous β-Scheduling**:
   - Use smooth β(t) function instead of step functions for guidance weights
   - Example: β(t) = -0.5 × (1 - t/T) to -0.5 to +0.5 smooth ramp
   - Reduces discontinuity artifacts at phase transitions

3. **Integrated CDR3 Exploration**:
   - Combine with CDR3-specific β-scaling (Idea L) for better loop conformation search
   - Phase 1: Negative β for CDR3 (exploratory), loose contact guidance
   - Phase 2: β → 0 for CDR3 (neutral refinement)
   - Phase 3: Positive β for CDR3 (optimization)

4. **Multi-Particle Refinement**:
   - Use FK particle resampling throughout (not just final selection)
   - Particles that maintain contacts through phases weighted higher
   - Reduces computational waste on poor trajectories

5. **Implementation Notes**:
   - Contact score computation at each step (not just phase transitions)
   - Threshold-based phase advancement: advance when score_change < threshold
   - Running average of contact scores to smooth noise

**Expected Improvement**:
- More natural exploration-to-exploitation progression
- Better handling of CDR3 conformational search
- 15-25% improvement in final contact quality

---

### Strategy J (Coupled CDR3 Heavy-Light Sampling) - Enhanced Version

**Current Implementation**: 2D grid of canonical H3/L3 combinations, joint scoring

**Improvements**:
1. **Asymmetric β-Scaling (Idea O)**: Replace 2D grid with differential β-scaling
   - β_H = +0.4 to +0.5 (heavy chain emphasis)
   - β_L = -0.1 to +0.2 (light chain lighter emphasis)
   - Biologically encodes CDR-H3 > CDR-L3 dominance

2. **Contact-Aware Interface Scoring**:
   - Score quality of CDR-H3 + CDR-L3 interface (complementarity)
   - Use structural metrics (distance, orientation, clash detection)
   - Downweight particles with H3/L3 clashes

3. **Beam Search Instead of Grid**:
   - Adaptive search rather than exhaustive grid
   - Identify good H3 conformations first, then optimize L3 given H3
   - Dramatically reduces sampling space (30 → 3-5 selected conformations)

4. **Ensemble-Based Selection**:
   - Generate ensemble of H3/L3 pairs that are compatible
   - Not just single best pair, but pareto-optimal frontier of (contact, complementarity)
   - Shows flexible binding modes

5. **Implementation Notes**:
   - Use contact score + interface quality metric in resampling weights
   - Apply both asymmetric β and explicit interface scoring
   - Track H3/L3 interface metrics in output

**Expected Improvement**:
- 50-70% reduction in particle count (beam search vs. grid)
- Better identification of compatible H3/L3 pairs
- More biologically realistic predictions

---

## Part 2: Enhancements to β-Scaling Ideas (K-P)

### Idea K (Region-Specific β-Scaling) - Detailed Implementation

**Enhancement**: Add confidence-guided masking

**Approach**:
1. **Adaptive Region Definition**: Use antigen surface properties
   - Compute surface accessibility scores (B-factors, RSA)
   - Identify concave regions (likely epitope pockets)
   - Create overlapping regions with adaptive size (15-25 Å radius)

2. **Multi-Pass Scanning with Confidence**:
   - Pass 1: Broad scan with moderate β (±0.3) to identify candidate regions
   - Compute pLDDT/PAE for contacts found in each region
   - Pass 2: Focused scan on high-confidence regions with stronger β (±0.5)

3. **Integration with Contact Propensity** (New Idea Q):
   - Use ML-based epitope propensity predictions (ESM-based methods)
   - Combine propensity scores with FK steering contact scores
   - Weight regions by (propensity × FK_contact) for final epitope ranking

**Implementation Notes**:
- Pre-compute surface accessibility once per antigen
- Run 5-10 β configurations in parallel (efficient GPU usage)
- Accumulate contacts across ensemble with confidence weighting

**Expected Benefits**:
- More robust epitope identification
- Better integration with ML epitope predictions
- Handles multiple potential epitope locations

---

### Idea L (CDR3-Specific β-Scaling) - Deployment Ready

**Status**: Highest priority for immediate implementation

**Enhancement**: Extend to explore CDR3 + antigen contact jointly

**Approach**:
1. **Two-Stage β-Scaling**:
   - Stage 1: CDR3-specific β (±0.5) for loop conformation exploration
   - Stage 2: Interface β (Idea M) to optimize contact given CDR3 conformation

2. **Ensemble Generation Strategy**:
   - Generate 5-10 conformations by sweeping CDR3 β in [-0.75, +0.75]
   - Filter by (antigen_contact_score + pLDDT) for final ensemble
   - No need for explicit canonical definitions

3. **Comparison with Strategy D**:
   - Simpler: no need for canonical structure database
   - More flexible: explores full conformation space, not just canonicals
   - Better for novel CDR3 sequences

**Implementation Notes**:
- Use existing pair masking infrastructure
- Simple linear sweep of β values
- Contact-based post-filtering

---

### Idea M (Inter-Chain Interface β-Scaling) - Combined with Binding Affinity

**Enhancement**: Add predicted binding energy weighting

**Approach**:
1. **Physics-Informed β-Scheduling**:
   - Use MM-PBSA or Rosetta FastRelax predictions to validate β ranges
   - Empirically determine optimal β_pos and β_neg for binding improvement
   - Create lookup table for different antigen classes

2. **Multi-Objective Optimization**:
   - Optimize for both contact score AND predicted binding ΔG
   - Use weighted combination: score = contact × (1/ΔG_predicted)
   - Balance exploration (β<0) with optimization (β>0)

**Implementation Notes**:
- Optional integration with fastfold/OmegaFold for binding prediction
- Pre-computed ΔG predictions for common antigen classes
- Fallback to simple contact scoring if ΔG unavailable

---

### Idea N (Adaptive β-Scheduling) - Contact-Triggered Phases

**Enhancement**: Add contact quality monitoring

**Approach**:
1. **Dynamic Phase Advancement**:
   - Monitor contact score trajectory at each diffusion step
   - Compute rolling average (window=10 steps)
   - Trigger phase advance when: score_improvement < threshold (e.g., 5% per 10 steps)

2. **Contact Plateau Detection**:
   - If contact score plateaus, switch to β > 0 (refinement)
   - If contacts poor, use β < 0 (exploration)
   - Adaptive scheduling responds to actual prediction quality

**Implementation Notes**:
- Requires contact computation at each diffusion step (may increase cost)
- Consider checkpointing to avoid recomputation
- Tunable threshold (recommend 5-10% improvement per phase)

---

### Idea O (Differential β-Scaling) - Biologically Informed

**Enhancement**: Add experimental validation targets

**Approach**:
1. **Task-Specific β Values**:
   - For affinity improvement: β_H = +0.4, β_L = +0.1
   - For epitope discovery: β_H = -0.3, β_L = -0.2 (exploratory)
   - For structure prediction: β_H = +0.2, β_L = 0 (H3 focus)

2. **Ablation Studies**:
   - Compare β_H only vs. β_H + β_L
   - Validate that asymmetric scaling improves over symmetric
   - Establish performance regression for each antigen class

**Implementation Notes**:
- Pre-defined β values for common use cases
- User override for custom applications
- Track H3/L3 dominance metrics in validation

---

### Idea P (Contact-Guided β-Feedback) - Multi-Pass Refinement

**Enhancement**: Add early stopping and convergence criteria

**Approach**:
1. **Iterative Refinement Cycle**:
   - Pass 1: Fast screening with β = -0.5 (10-20% diffusion steps)
   - Rank by contact score
   - Pass 2: Refinement with adaptive β based on Pass 1 results:
     - Top 30%: β = +0.5 (optimize)
     - Middle 40%: β = 0 (neutral)
     - Bottom 30%: β = -0.75 (radical exploration)
   - Final Pass: Select top 10% by (contact + pLDDT)

2. **Convergence Criteria**:
   - Stop if same conformations selected in consecutive passes
   - Combine results when ensemble diversity < threshold
   - Automatic convergence detection

**Implementation Notes**:
- Can reduce total compute time vs. single long pass
- Parallel execution of different β values
- Requires contact score computation between passes

---

## Part 3: Novel Hybrid Approaches

### New Idea Q: Iterative Epitope Refinement with FK + β-Scaling

**Motivation**: Combine the best of FK particle resampling and β-scaling for rapid epitope discovery and refinement

**Approach**:
1. **Round 1 - Broad Exploration**:
   - Run FK particles with weak contact guidance (coefficient = 0.1)
   - Use exploratory β (±0.3) to prevent premature convergence
   - Generate initial epitope contact heatmap

2. **Round 2 - Targeted Refinement**:
   - Identify hotspot regions (top 20% contact frequency)
   - Focus region-specific β-scaling to hotspot neighborhoods
   - Run FK particles with stronger contact guidance (coefficient = 0.5)

3. **Round 3 - Validation**:
   - Generate final structures with optimized β values from Round 2
   - Validate against experimental epitope data (if available)
   - Report confidence scores with ensemble metrics

**Expected Outcome**:
- Discovers novel epitopes in 3 iterations
- Combines FK particle exploration with β-scaling efficiency
- Confidence scores indicate prediction reliability

**Implementation Complexity**: Medium

**Computational Cost**: 3× single pass (parallelizable)

---

### New Idea R: Multi-Objective Steering with Composite β-Scaling

**Motivation**: Current approaches optimize single objective (contact or conformation). Real antibody-antigen systems need balance.

**Approach**:
1. **Simultaneous Objectives**:
   - Objective 1: Maximize CDR-antigen contacts
   - Objective 2: Maintain physically realistic CDR3 conformations
   - Objective 3: Minimize internal clashes (antibody self-compatibility)

2. **Composite β-Scaling**:
   - β_contact for contact-promoting pairs
   - β_cdr3 for CDR3 conformation space
   - β_clash for penalty regions
   - Apply all three simultaneously with careful weight tuning

3. **Pareto-Optimal Ensemble**:
   - Instead of single best prediction, generate ensemble on pareto frontier
   - Shows trade-offs between objectives
   - Useful for understanding structure stability

**Implementation Complexity**: High (requires multi-objective weighting)

**Computational Cost**: Slight increase over single-objective

---

### New Idea S: Template-Guided β-Scheduling

**Motivation**: Known antibody-antigen complex structures contain information about good CDR conformations and antigen orientations

**Approach**:
1. **Template Encoding**:
   - Extract CDR3 backbone angles from PDB templates (human/synthetic antibodies)
   - Compute antigen orientation angles relative to CDR loops
   - Build lookup table: (antibody_species, antigen_class) → (ideal_H3_angles, ideal_antigen_orientation)

2. **Template-Informed β-Schedule**:
   - At diffusion step t, compute deviation from template
   - If far from template: use β = -0.3 (explore, but not too far)
   - If near template: use β = +0.3 to refine (exploit)
   - Smooth interpolation based on structural similarity

3. **Species-Aware Guidance**:
   - Different templates for human, mouse, synthetic antibodies
   - Incorporate species information in β-schedule
   - Better for species-specific predictions

**Expected Benefits**:
- Incorporates structural knowledge without hard constraints
- Soft biasing (β) more flexible than hard templates
- Useful for humanization/library design tasks

**Implementation Complexity**: Medium (requires template database and RMSD calculation)

---

### New Idea T: MSA-Free Steering Strategies for Low-Information Systems

**Motivation**: Many designed or synthetic antibodies have limited MSA. Boltz-sample shows β-scaling works without MSA, but current steering relies on contacts (which may fail without good MSA).

**Approach**:
1. **MSA-Agnostic Contact Scoring**:
   - Use sequence-based features only (no aligned sequences)
   - Contact score from Boltz2's internal distance predictions (pAE-derived)
   - More robust than external contact definition

2. **Intrinsic β-Scaling**:
   - Don't require external potentials, use latent space properties
   - β values act on model's internalized priors
   - Works for novel sequences with no homologs

3. **Confidence-Weighted Ensemble**:
   - Generate ensemble with varying β values
   - Use pLDDT/pAE as confidence (not contact-based)
   - Select ensemble based on internal model confidence

**Expected Benefits**:
- De novo antibody design (no MSA)
- Better for synthetic/engineered antibodies
- More robust to MSA quality/size

**Implementation Complexity**: Low-Medium (mostly parameter selection)

---

### New Idea U: Conformational Energy Landscape Mapping via FK Steering

**Motivation**: Not just find best structure, but map entire accessible conformational space

**Approach**:
1. **Energy Landscape Sweep**:
   - Generate multiple predictions at different β values: [-1.0, -0.75, -0.5, -0.25, 0, +0.25, +0.5, +0.75, +1.0]
   - Each β samples a different region of latent space
   - Collect all conformations (not just best)

2. **Landscape Characterization**:
   - Plot: x-axis = β value, y-axis = contact score
   - Identify local optima, plateaus, discontinuities
   - Understand which β values lead to similar structures

3. **Optimal Basin Identification**:
   - Cluster conformations by similarity
   - Identify "good" basins (high contact + high diversity)
   - Recommend β ranges for application-specific goals

**Expected Outcome**:
- Complete picture of accessible structure space
- Identifies robust vs. fragile binding modes
- Guides β parameter selection for future runs

**Implementation Complexity**: Low (just run multiple β values)

**Computational Cost**: 9× single pass (parallelizable)

---

## Part 4: Integration with FK-Diffusion Framework

### Key Insights from FK-Diffusion Steering Paper

The FK-Diffusion paper (Horvitz et al., 2501.06848) provides a general framework for inference-time steering using **Feynman-Kac resampling**:

**Core Mechanism**:
- Multiple particle trajectories during diffusion
- At each step, compute reward/potential energy
- Resample particles: w_i ∝ exp(-E_i) (low energy = high weight)
- Poor particles eliminated, good particles replicated

**Advantages**:
1. Arbitrary reward functions (no differentiability required)
2. Gradient-free guidance (works with any energy function)
3. Proven effective for protein design and molecular glue prediction
4. Minimal computational overhead

### Hybrid Approaches: FK + β-Scaling

**Opportunity**: Combine Feynman-Kac resampling with β-scaling for superior performance

**Approach 1: FK Resampling on β-Scaled Predictions**
```
For each β value in sweep [-0.75, 0.75]:
  Generate multiple predictions with β-scaled pairs
  Score each by contact + binding affinity
  Use FK weights to resample at each diffusion step
  Final ensemble: highest-weighted predictions
```
**Benefit**: Combines exploratory power (β-scaling) with adaptive filtering (FK resampling)

**Approach 2: Dual-Track Steering**
```
Particles track 1: Use β-scaling (efficient, smooth)
Particles track 2: Use explicit potentials (precise, targeted)
At each step: Compare progress, share information between tracks
Final result: Best of both approaches
```
**Benefit**: Robust to individual steering limitations

---

## Part 5: Implementation Priority and Roadmap

### Immediate Priority (Phase 1, Week 1-2)

1. **Idea L Enhancement**: CDR3-specific β-scaling with two-stage approach
   - Lowest complexity, highest impact
   - Deprecates need for explicit canonical database
   - Direct replacement for Strategy D

2. **Idea K Enhancement**: Region-specific β-scaling with adaptive regions
   - Moderate complexity, high visibility
   - Direct replacement for Strategy E
   - Can run in parallel

3. **Idea N Enhancement**: Contact-triggered adaptive phase scheduling
   - Moderate complexity, good improvement to existing G
   - Direct enhancement to existing code path

### Near-term Priority (Phase 2, Week 3-4)

4. **New Idea Q**: Iterative epitope refinement
   - Combines FK + β-scaling
   - Novel contribution, publication-worthy
   - Parallelizable

5. **Idea O Enhancement**: Differential β-scaling with biological motivation
   - Proven approach, straightforward to implement
   - Improves upon Strategy J

### Medium-term Priority (Phase 3, Week 5-6)

6. **New Idea T**: MSA-free steering
   - Important for synthetic/designed antibodies
   - Extends applicability

7. **New Idea R**: Multi-objective composite β-scaling
   - Higher complexity but scientifically interesting
   - Good for publication

### Advanced Priority (Phase 4+)

8. **New Idea S**: Template-guided β-scheduling
9. **New Idea U**: Conformational energy landscape mapping
10. **Hybrid FK + β-scaling approaches**

---

## Part 6: Validation and Testing Strategy

### Benchmark Datasets
- Existing PDB antibody-antigen complexes (known epitopes)
- SABDAB database for diverse antibody structures
- Designed synthetic antibodies (limited MSA availability)

### Metrics
1. **Epitope Discovery**: Precision/recall of identified residues vs. experimental
2. **Binding Affinity**: Correlation with SASA, buried surface area
3. **Structure Quality**: RMSD vs. native structure (when available)
4. **Ensemble Diversity**: pairwise RMSD within ensemble
5. **Computational Cost**: Time per prediction, GPU memory

### Validation Protocol
```
For each strategy S:
  For each benchmark antibody-antigen pair:
    Run prediction with S
    Measure all metrics
    Compare to baseline (no steering)
  Report mean/std metrics
  Statistical significance testing
```

---

## Part 7: Key Research Questions

### Fundamental Questions

1. **β-Range Optimality**: What is the optimal β range for each strategy?
   - Does it vary by antibody species?
   - Does it vary by antigen class?
   - Can we automatically learn optimal β?

2. **Computational Overhead**:
   - How much does contact score computation at each step cost?
   - Can we use approximations (e.g., pAE-based contact)?
   - Optimal frequency of computation?

3. **Multi-Objective Trade-offs**:
   - How to balance contact optimization vs. structure quality?
   - Can we solve multi-objective problem explicitly (Pareto frontier)?
   - What weight ratios work best?

4. **Ensemble Diversity**:
   - What ensemble size (5 vs. 10 vs. 20 structures) is optimal?
   - How to measure meaningful diversity?
   - When is diversity harmful (overfitting)?

5. **MSA Dependency**:
   - Can we quantify performance vs. MSA size?
   - Are β values different for high vs. low MSA?
   - How to detect "bad" MSA and adjust steering?

### Technical Questions

1. **Phase Transition Criteria**: What triggers phase advancement in Idea N?
   - Fixed timestep (current approach)
   - Contact-based threshold
   - Prediction uncertainty (pAE)-based
   - Running average convergence

2. **Contact Score Definition**: Which contact metric is best?
   - Distance-based (current)
   - pAE-based (more robust)
   - ML-predicted binding energy
   - Combination (composite)

3. **β Scheduling Smoothness**:
   - Sharp step functions vs. continuous ramps
   - Smooth transitions vs. rapid switches
   - Effect on final structure quality

4. **Template Selection**:
   - Which PDB templates to use as references?
   - Similarity metrics for template matching
   - When to ignore templates (novel designs)

---

## Part 8: References and Related Work

### Key Papers

1. **FK-Diffusion Steering**: Horvitz et al. (2501.06848)
   - General framework for inference-time steering
   - Arbitrary reward functions
   - Minimal overhead

2. **Boltz-sample**: Suzuki & Amagasa (2026)
   - Pair representation β-scaling
   - Direct control of latent space
   - Works without MSA

3. **FKSFold**: Related work on FK-steered complex structure prediction

4. **Epitope Prediction**:
   - ParaDeep, EPP (recent 2025 methods)
   - ESM-based epitope propensity
   - Graph-based antibody-aware prediction

### Implementation Resources

- FK-Diffusion GitHub: https://github.com/zacharyhorvitz/Fk-Diffusion-Steering
- Boltz-sample GitHub: https://github.com/suzuki-2001/boltz-sample
- Boltz-2 GitHub: jwohlwend/boltz

---

## Part 9: EmbedOpt-Inspired Improvements (Robust Embedding Space Steering)

### Key Insight from EmbedOpt Paper

The "Robust Inference-Time Steering of Protein Diffusion Models via Embedding Optimization" paper (Li et al., 2602.05285) introduces a critical insight:

**Embedding space optimization is more robust than coordinate space optimization** for steering diffusion models, especially for experimental constraints outside the training distribution.

**Advantages**:
- Stable across 2 orders of magnitude of hyperparameter variation
- Better for out-of-distribution constraints
- Fewer aggressive parameter tuning requirements
- Captures "sequence and coevolutionary signals" in latent representation

**Application to Antibody-Antigen Steering**:
Instead of steering in atomic coordinate space (current approach with potentials), we can steer in embedding space where CDR and antigen properties are encoded more stably.

---

### Idea V: Embedding Space CDR3 Steering (EmbedOpt-Inspired)

**Goal**: More robust CDR3 conformation exploration for novel sequences

**Motivation**: Current CDR3 steering uses explicit dihedral potentials (coordinate space). For novel CDR sequences, this can be brittle. EmbedOpt shows embedding space optimization is more robust.

**Approach**:
1. **CDR3 Embedding Optimization**: Instead of targeting psi angles in coordinate space
   - Optimize CDR3 representation in Boltz2's embedding layer
   - Define target embedding based on known CDR3 conformers (from PDB)
   - Use gradient descent in embedding space toward target

2. **Implementation Strategy**:
   - Extract CDR3 embeddings from Pairformer (or earlier layers)
   - Define CDR3_target_embedding from canonical structures
   - Compute loss: L = ||CDR3_embedding - CDR3_target||²
   - Optimize embedding, then run final denoising steps
   - No need to specify dihedral angles explicitly

3. **Robustness Benefits**:
   - Works for novel CDR sequences not in training set
   - Stable across wide range of optimization strengths
   - No hyperparameter tuning per antibody needed
   - Better generalization to synthetic antibodies

4. **Integration with β-Scaling**:
   - Embedding optimization (EmbedOpt) as primary steering
   - β-scaling (Boltz-sample) as secondary refinement
   - Hybrid approach: two levels of latent space control

**Expected Benefits**:
- 50%+ improvement for novel/synthetic CDR sequences
- More stable than coordinate potentials
- Better for out-of-distribution antibodies
- Reduced hyperparameter tuning

**Implementation Complexity**: Medium (requires embedding extraction and gradient computation)

**Target Goal**: **Structure Prediction Improvement** for novel sequences + **Robustness**

---

### Idea W: Embedding-Based Interface Steering for Binding

**Goal**: Robust antigen-CDR interface optimization using embedding space

**Motivation**: Current antigen steering uses contact distances (coordinate space). EmbedOpt suggests embedding space captures interface properties more robustly.

**Approach**:
1. **Interface Embedding Definition**:
   - Extract embeddings for CDR atoms and antigen atoms
   - Define interface embedding as pairwise interaction features
   - Target: strong interaction embeddings between CDR and antigen

2. **Embedding-Space Contact Optimization**:
   - Instead of maximizing Euclidean distance contacts
   - Optimize for high "interaction potential" in embedding space
   - Use learned representations of binding compatibility

3. **Implementation**:
   - Extract CDR-antigen pair embeddings
   - Define loss: minimize distance between observed and ideal interaction embeddings
   - Gradient descent in embedding space
   - Propagate back through network layers

4. **Advantages Over Coordinate-Space Steering**:
   - Captures binding specificity (not just proximity)
   - Works for novel antigen-antibody pairs
   - More stable hyperparameter behavior
   - Better alignment with what model learned about binding

**Expected Benefits**:
- Better binding affinity prediction (understands specificity, not just distance)
- Robust to novel antigens
- 20-30% improvement in interface quality
- Natural integration with embedding-based predictions

**Implementation Complexity**: Medium-High (requires careful embedding layer identification)

**Target Goal**: **Binding Affinity Improvement** + **Robustness to Novel Antigens**

---

### Idea X: Experimental Constraint Integration via EmbedOpt

**Goal**: Integrate experimental data (cryo-EM, cross-linking, HDX-MS) into Boltz-2 predictions

**Motivation**: EmbedOpt was designed specifically for experimental constraints. Antibody-antigen complexes have experimental data available.

**Approach**:
1. **Multiple Constraint Types**:
   - Cryo-EM density maps → constrain overall shape in embedding space
   - Cross-linking distance constraints → constrain inter-chain distances
   - HDX-MS protection data → constrain interface residues
   - SPR data → constrain binding orientation

2. **Embedding-Space Constraint Integration**:
   - Convert experimental constraints to embedding-space objectives
   - Optimize embeddings to satisfy constraints
   - More robust than coordinate-space constraint application

3. **Multi-Modal Fusion**:
   - Combine multiple constraint types simultaneously
   - Weight by experimental uncertainty
   - Natural probabilistic fusion in latent space

4. **Implementation**:
   - Pre-compute experimental constraint functions
   - Map to embedding space objectives
   - Use EmbedOpt-style gradient descent in embeddings
   - Final refinement with contact-aware guidance

**Expected Benefits**:
- Leverage experimental antibody-antigen data
- 3-5x better accuracy when experimental data available
- Robust multi-modal constraint integration
- Novel way to use cryo-EM/SAXS/cross-linking data

**Implementation Complexity**: High (requires constraint translation to embedding space)

**Target Goal**: **Structure Prediction Accuracy** (when experimental data available) + **Robustness**

---

### Idea Y: Hierarchical Steering - Embedding First, Then Coordinates (Hybrid)

**Goal**: Combine robustness of embedding steering with precision of coordinate steering

**Motivation**: EmbedOpt is robust but may lose fine details. Coordinate potentials are precise but brittle. Use both.

**Approach**:
1. **Stage 1 - Embedding Steering** (t=T to t=T/2):
   - Use EmbedOpt-style embedding optimization
   - Coarse-grained steering toward target conformations
   - Robust even if constraints are out-of-distribution
   - Explores conformational space broadly

2. **Stage 2 - Coordinate Steering** (t=T/2 to t=0):
   - Switch to coordinate-space potentials
   - Fine-grained optimization of interface
   - Refine contact geometry and hydrogen bonding
   - Exploit local energy landscape

3. **Smooth Transition**:
   - Gradually reduce embedding steering weight
   - Gradually increase coordinate steering weight
   - Natural blend during middle diffusion steps

4. **Combining Benefits**:
   - Robustness of embedding space (early stages)
   - Precision of coordinates (late stages)
   - Best of both worlds

**Expected Benefits**:
- Robustness of novel sequences (embedding stage)
- High-precision contacts (coordinate stage)
- Natural progression matching structure formation
- Fewer total hyperparameters to tune

**Implementation Complexity**: High (requires dual steering infrastructure)

**Target Goal**: **Structure Prediction Quality** (robust + precise) + **Binding Affinity**

---

### Idea Z: Soft Experimental Constraints via Embedding Scoring

**Goal**: Avoid over-constraining when experimental data has uncertainty

**Motivation**: EmbedOpt can handle constraints robustly, but real experimental data has noise. Use probabilistic soft constraints.

**Approach**:
1. **Probabilistic Constraint Model**:
   - Treat experimental measurements as probabilistic constraints
   - Use Bayesian interpretation: p(structure | experiment)
   - Model uncertainty: wide distributions for noisy measurements

2. **Soft Embedding Objectives**:
   - Instead of hard constraints, use soft probability-weighted objectives
   - Embeddings optimized to likely region, not exact target
   - Allows flexibility around noisy measurements

3. **Uncertainty-Aware Steering**:
   - High-confidence measurements → tight embedding constraints
   - Low-confidence measurements → loose embedding constraints
   - Adaptive weighting based on measurement precision

4. **Implementation**:
   - Quantify experimental uncertainty per constraint
   - Compute Bayesian soft targets in embedding space
   - Weight objectives by 1/uncertainty²
   - Natural probabilistic steering

**Expected Benefits**:
- Better handling of noisy experimental data
- Avoids over-fitting to imperfect measurements
- Probabilistically principled approach
- Robust when experimental precision varies

**Implementation Complexity**: Medium (requires uncertainty quantification)

**Target Goal**: **Structure Prediction Accuracy** (with noisy experiments) + **Robustness**

---

## Summary Table: All Improvements and New Ideas

| ID | Type | Name | Goal | Complexity | Priority | Status |
|----|----|------|------|-----------|----------|--------|
| D+ | Enhancement | Canonical Ensemble v2 | Structure | Medium | Phase 1 | New impl |
| E+ | Enhancement | Blind Scanning v2 | Epitope | High | Phase 1 | New impl |
| G+ | Enhancement | Progressive Refinement v2 | Structure | Medium | Phase 1 | New impl |
| J+ | Enhancement | Coupled CDR3 v2 | Binding | Medium | Phase 2 | New impl |
| K+ | Enhancement | Region-Specific β-Scaling v2 | Epitope | Medium | Phase 1 | Refine |
| L+ | Enhancement | CDR3 β-Scaling v2 | Structure | Low | Phase 1 | Ready |
| M+ | Enhancement | Interface β-Scaling v2 | Binding | Low | Phase 1 | Ready |
| N+ | Enhancement | Adaptive β-Scheduling v2 | Structure | Medium | Phase 1 | Ready |
| O+ | Enhancement | Asymmetric β-Scaling v2 | Binding | Low | Phase 2 | Ready |
| P+ | Enhancement | Contact-Guided β-Feedback v2 | Binding | Medium | Phase 2 | Ready |
| Q | Novel | Iterative Epitope Refinement | Epitope | Medium | Phase 2 | New |
| R | Novel | Multi-Objective Steering | Binding+Struct | High | Phase 3 | New |
| S | Novel | Template-Guided β-Scheduling | Design | Medium | Phase 3 | New |
| T | Novel | MSA-Free Steering | Design | Medium | Phase 2 | New |
| U | Novel | Energy Landscape Mapping | Analysis | Low | Phase 3 | New |
| V | Novel (EmbedOpt) | Embedding-Space CDR3 Steering | Structure+Robust | Medium | Phase 2 | New |
| W | Novel (EmbedOpt) | Embedding-Based Interface Steering | Binding+Robust | Medium-High | Phase 2 | New |
| X | Novel (EmbedOpt) | Experimental Constraint Integration | Structure+Robust | High | Phase 3 | New |
| Y | Novel (Hybrid) | Hierarchical Steering (Embed→Coords) | Structure+Binding | High | Phase 3 | New |
| Z | Novel (EmbedOpt) | Soft Experimental Constraints | Structure+Robust | Medium | Phase 3 | New |

---

## Next Steps

1. **Document Selection**: Choose strategies for Phase 1 implementation
2. **Detailed Design Docs**: Create implementation specs for each selected strategy
3. **Code Review**: Evaluate existing codebase for integration points
4. **Experimental Design**: Create testing protocols and benchmarks
5. **Implementation**: Build selected strategies sequentially
6. **Validation**: Run benchmarks and compare to baselines
7. **Publication**: Prepare manuscript of results

---

**Document Status**: Design Phase - Ready for implementation planning
**Last Updated**: 2025-02-09
