# Implementation Guide: How to Build the Improvements

**Date**: 2025-02-09
**Purpose**: Concrete implementation suggestions for each improvement (without making code changes yet)
**Scope**: Code structure, module changes, integration points

---

## Part 1: Enhancement Implementations (D+ through P+)

### Strategy D+ - Canonical Ensemble v2 with β-Scaling

**Current Approach**: Explicit dihedral potentials for canonical CDR conformations

**Proposed Change**: Replace potentials with CDR3-specific β-scaling

**Code Changes Needed**:

1. **In pair_scaling_module.py** (or equivalent):
   - Add function `apply_cdr_region_beta_scaling(z, beta, cdr_h3_mask, cdr_l3_mask)`
   - Parameters: β value, masks for CDR regions
   - Apply: `z[cdr_mask] *= (1 + beta)`

2. **In potential_manager.py** (or inference config):
   - Deprecate `CDR3ConformationPotential` class for this use case
   - Add conditional: if use_beta_scaling: skip explicit potentials
   - Track β values used (for logging)

3. **In ensemble_generation.py**:
   - Sweep β from -0.5 to +0.5 in steps of 0.1 (11 values)
   - Generate ensemble of 11 predictions
   - Score by contact + confidence (pLDDT)
   - Select top-K (suggest K=5) for final output

4. **Configuration**:
   - Add CLI flag: `--use_cdr3_beta_scaling` (boolean)
   - Add parameter: `--cdr3_beta_range` (default: [-0.5, +0.5])
   - Add parameter: `--cdr3_beta_steps` (default: 11)

**Testing Approach**:
```python
# Pseudocode for validation
def test_cdr3_beta_scaling():
    for beta in [-0.5, 0, 0.5]:
        pred = boltz_predict(config, cdr3_beta=beta)
        # Check: CDR3 angles should vary across beta values
        assert psi_angle_diversity > threshold
        # Check: Contact scores should improve from negative to positive beta
        contact_scores_should_increase()
```

---

### Strategy E+ - Blind Scanning v2 with Region-Specific β-Scaling

**Current Approach**: Run N separate diffusion passes with different antigen spatial constraints

**Proposed Change**: Single pass with region-specific β-scaling via masking

**Code Changes Needed**:

1. **In antigen_partitioning.py** (new module):
   - Function: `partition_antigen_surface(antigen_coords, method='curvature')`
   - Methods: spherical regions (current), surface curvature, accessibility-based
   - Output: List of region masks (shape: [num_residues])
   - Suggest 8-12 overlapping regions per antigen

2. **In pair_scaling_module.py**:
   - Function: `apply_region_specific_beta_scaling(z, beta_config, antigen_masks)`
   - For each region: apply β to (CDR, antigen_region) pairs
   - Handle overlapping regions (additive or max)
   - Profile: should be O(n_pairs) = minimal overhead

3. **In inference_loop.py**:
   - New path: `inference_with_region_sweep()`
   - Generate N predictions with different region-emphasis β values
   - Each β value emphasizes one region and de-emphasizes others
   - Example: β_region_i = +0.5, β_other = -0.3

4. **In contact_tracking.py**:
   - Accumulate contact heatmap across all N predictions
   - Weight contacts by pLDDT + pAE confidence
   - Generate epitope propensity map per antigen residue
   - Output: coordinates of epitope hotspots

5. **Configuration**:
   - CLI flag: `--epitope_scanning` (boolean)
   - Parameter: `--scan_regions` (default: 10)
   - Parameter: `--region_beta_emphasis` (default: 0.5)
   - Parameter: `--region_beta_deemphasis` (default: -0.3)

**Algorithm Pseudocode**:
```python
def epitope_scanning(ab_seq, ag_seq, num_regions=10):
    ag_masks = partition_antigen(ag_seq, regions=num_regions)
    contact_heatmap = np.zeros(len(ag_seq))

    for region_id in range(num_regions):
        # Emphasize this region
        beta_config = {region_id: 0.5}  # Emphasize
        for other_id in range(num_regions):
            if other_id != region_id:
                beta_config[other_id] = -0.3  # De-emphasize

        # Single diffusion pass with this beta configuration
        pred = boltz_predict_with_beta_scaling(
            ab_seq, ag_seq,
            beta_config=beta_config,
            ag_masks=ag_masks
        )

        # Accumulate contacts
        cdr_coords, ag_coords = extract_structures(pred)
        contacts = compute_contacts(cdr_coords, ag_coords)
        contact_heatmap += weight_by_confidence(contacts, pred.pLDDT, pred.pAE)

    return contact_heatmap / num_regions
```

**Testing**:
```python
def test_epitope_scanning():
    # Known antibody-antigen pair
    contact_heatmap = epitope_scanning(ab_seq, ag_seq)

    # Check: Top contacted regions should match known epitope
    precision, recall = evaluate_against_known_epitope(contact_heatmap)
    assert precision > 0.7 and recall > 0.6
```

---

### Strategy G+ - Progressive Refinement v2 with Contact-Triggered Phases

**Current Approach**: Three fixed phases at predetermined diffusion timesteps

**Proposed Change**: Contact-triggered phase transitions with smooth β-scheduling

**Code Changes Needed**:

1. **In diffusion_loop.py**:
   - Modify loop structure to compute contact scores at each step
   - Track running average of contact improvement
   - Condition on improvement threshold (e.g., 5% improvement per 10 steps)

2. **In phase_scheduling.py** (new module):
   - Class: `AdaptivePhaseScheduler`
   - Method: `get_current_phase(current_contact_score, running_avg_improvement, timestep)`
   - Returns: (phase_id, guidance_weight, beta_value)
   - Smooth interpolation between phases instead of step function

3. **In guidance_weights.py**:
   - Replace step function with smooth ramp
   - Phase 1: β(t) = -0.5, guidance = 0.1
   - Phase 2: β(t) = -0.5 + t_progress × 0.5, guidance = 0.1 + t_progress × 0.4
   - Phase 3: β(t) = 0, guidance = 0.5
   - Where t_progress = fraction of phase completed (0 to 1)

4. **Contact Score Integration**:
   - Compute contact distance at each step: `d_cdr_ag = min_distance(CDR_atoms, Ag_atoms)`
   - Rolling average: `avg = (9 × prev_avg + new_distance) / 10` (window=10)
   - Check improvement: `improvement = (prev_avg - current) / prev_avg`
   - Trigger phase advance when: `improvement < threshold` for N_consecutive steps

5. **Configuration**:
   - CLI flag: `--adaptive_phases` (boolean)
   - Parameter: `--phase_improvement_threshold` (default: 0.05 = 5%)
   - Parameter: `--phase_improvement_window` (default: 10 steps)
   - Parameter: `--contact_computation_stride` (default: 1 = every step, or 2/5 to reduce cost)

**Pseudocode**:
```python
def adaptive_diffusion_loop(model, noise_schedule, config):
    scheduler = AdaptivePhaseScheduler()
    contact_history = deque(maxlen=10)

    for t in range(T, 0, -1):
        # Standard diffusion step
        x_pred = model(x_t, t)

        # Contact-based steering decision (optional: compute every N steps)
        if t % config.contact_stride == 0:
            contacts = compute_cdr_ag_contacts(x_pred)
            contact_history.append(contacts)
            avg_improvement = compute_improvement(contact_history)

            # Adaptive phase scheduling
            phase, guidance_weight, beta_value = scheduler.get_parameters(
                contact_score=contacts,
                improvement=avg_improvement,
                timestep=t
            )

        # Apply steering with current parameters
        if config.use_beta_scaling:
            x_pred = apply_beta_scaling(x_pred, beta_value)

        if config.use_potentials:
            energy = contact_potential(x_pred)
            gradient = compute_gradient(energy, x_pred)
            x_pred = x_pred - guidance_weight * gradient

        # Continue diffusion
        noise = sample_noise()
        x_t_minus_1 = x_pred + noise
```

---

### Strategy J+ - Coupled CDR3 v2 with Asymmetric β-Scaling and Beam Search

**Current Approach**: 2D grid of canonical H3/L3 combinations, exhaustive evaluation

**Proposed Change**: Asymmetric β-scaling with beam search

**Code Changes Needed**:

1. **In asymmetric_beta_scaling.py** (new module):
   - Function: `apply_asymmetric_cdr_beta(z, beta_h, beta_l, h3_mask, l3_mask)`
   - Apply different β values: z[h3_mask] *= (1 + beta_h), z[l3_mask] *= (1 + beta_l)
   - Default values: beta_h=+0.4, beta_l=+0.1

2. **In beam_search_module.py** (new module):
   - Phase 1: Sample H3 conformations with β sweep (-0.3 to +0.5, 5 values)
   - Score each by contact with antigen
   - Keep top-3 H3 conformations (beam width)

   - Phase 2: For each top H3, sample L3 with β sweep
   - Score by (CDR-Ag contact + H3-L3 interface quality)
   - Keep top-2 per H3 (total: 3×2=6 combinations)

3. **In interface_scoring.py**:
   - Compute H3-L3 backbone interface metrics:
     - Minimum distance between H3 and L3 atoms
     - Number of atoms closer than 4.0 Å (clash count)
     - Dihedral angle compatibility
   - Return: interface_quality_score in [0, 1]

4. **In ensemble_generation.py**:
   - Replace grid expansion with beam search
   - Significantly fewer particles needed (6 vs. potentially 30+)
   - Cost reduction proportional to grid size

5. **Configuration**:
   - CLI flag: `--coupled_cdr3_steering` (boolean)
   - Parameter: `--cdr3_beam_width` (default: 3 for H3, 2 for L3)
   - Parameter: `--beta_h_range` (default: [-0.3, 0.5])
   - Parameter: `--beta_l_range` (default: [-0.2, 0.2])
   - Parameter: `--h3_l3_interface_weight` (default: 0.3)

**Pseudocode**:
```python
def coupled_cdr3_steering(ab_seq, ag_seq, config):
    # Phase 1: H3 conformation sampling
    h3_candidates = []
    for beta_h in np.linspace(-0.3, 0.5, 5):
        pred = boltz_predict_with_beta(
            ab_seq, ag_seq,
            beta_h3=beta_h,
            beta_l3=0  # Neutral L3
        )
        contact_score = compute_cdr_ag_contact(pred)
        h3_candidates.append((pred, contact_score, beta_h))

    # Beam search: keep top-3
    h3_beam = sorted(h3_candidates, key=lambda x: x[1], reverse=True)[:3]

    # Phase 2: For each top H3, sample L3
    final_pairs = []
    for h3_pred, h3_contact, h3_beta in h3_beam:
        l3_candidates = []
        for beta_l in np.linspace(-0.2, 0.2, 3):
            # Fix H3 conformation, vary L3
            pred = boltz_predict_with_beta(
                ab_seq, ag_seq,
                beta_h3=h3_beta,  # Keep this H3's beta
                beta_l3=beta_l    # Vary L3
            )

            ag_contact = compute_cdr_ag_contact(pred)
            h3_l3_interface = compute_h3_l3_interface_quality(pred)
            score = ag_contact + config.interface_weight * h3_l3_interface

            l3_candidates.append((pred, score))

        # Keep top-2 L3 for this H3
        l3_beam = sorted(l3_candidates, key=lambda x: x[1], reverse=True)[:2]
        final_pairs.extend(l3_beam)

    return final_pairs
```

---

## Part 2: β-Scaling Enhancement Implementations (K+ through P+)

### Idea K+ - Region-Specific β-Scaling v2

**Same as Strategy E+** (they use the same underlying mechanism)

See Strategy E+ implementation above.

---

### Idea L+ - CDR3-Specific β-Scaling v2 (Ready to Implement)

**Key Feature**: Simplest enhancement, can integrate immediately

**Code Changes Needed**:

1. **In pair_scaling_module.py**:
   - Add function: `apply_cdr3_beta_scaling(z, beta, heavy_chain_mask, light_chain_mask, cdr3_region_bounds)`
   - Bounds for CDR3: H3 residues 96-111, L3 residues 90-98
   - Apply: `z[cdr3_region] *= (1 + beta)` where cdr3_region = union of H3 and L3

2. **New CLI flags**:
   - `--cdr3_beta_scaling` (boolean, simpler than explicit potentials)
   - `--cdr3_beta_value` (float, default: 0.0)
   - `--cdr3_beta_sweep` (boolean) - if True, generate ensemble with β sweep

3. **Integration with existing code**:
   - If both `--cdr3_steering` (explicit psi angles) and `--cdr3_beta_scaling` enabled → use β-scaling only (β is simpler)
   - Log which method used

4. **Ensemble generation**:
   - Option 1: Single prediction with `--cdr3_beta_value 0.3`
   - Option 2: Ensemble with `--cdr3_beta_sweep --cdr3_beta_steps 5` (generates 5 predictions with β ∈ [-0.5, +0.5])

**Testing**:
```python
def test_cdr3_beta_scaling():
    # Single prediction
    pred_pos = boltz_predict(config, cdr3_beta=+0.3)  # Optimize CDR3
    pred_neg = boltz_predict(config, cdr3_beta=-0.3)  # Explore CDR3

    # Check: Conformations should differ significantly
    h3_rmsd = compute_rmsd(pred_pos.cdr3_h, pred_neg.cdr3_h)
    assert h3_rmsd > 1.0  # At least 1 Angstrom difference

    # Check: Both should maintain reasonable structure
    assert pred_pos.pLDDT > 50 and pred_neg.pLDDT > 50
```

---

### Idea M+ - Inter-Chain Interface β-Scaling v2

**Enhancement**: Add predicted binding energy weighting

**Code Changes Needed**:

1. **In interface_beta_scaling.py**:
   - Function: `apply_interface_beta_scaling(z, beta, cdr_mask, antigen_mask)`
   - For pairs (i, j) where i in CDR and j in antigen:
     - `z[i, j] *= (1 + beta)`

2. **Optional: Binding Energy Weighting**:
   - Function: `compute_binding_energy_weighted_beta(pred, antigen_seq)`
   - Use simple SASA-based ΔG estimate (or integrate fastfold if available)
   - Adjust β based on predicted ΔG: if ΔG_pred > threshold, use β = +0.5; else β = 0

3. **Configuration**:
   - CLI flag: `--interface_beta_scaling` (boolean)
   - Parameter: `--interface_beta` (default: 0.3)
   - Parameter: `--use_binding_energy_weighting` (boolean)
   - Parameter: `--binding_energy_threshold` (default: -5.0 kcal/mol estimate)

---

### Idea N+ - Adaptive β-Scheduling v2

**Same as Strategy G+** implementation above (uses contact-triggered phases with smooth β scheduling)

---

### Idea O+ - Asymmetric β-Scaling v2

**Same as Strategy J+** implementation above (β_h vs β_l for heavy vs light chain)

---

### Idea P+ - Contact-Guided β-Feedback v2

**Code Changes Needed**:

1. **In multi_pass_steering.py** (new module):
   - Function: `multi_pass_contact_feedback(ab_seq, ag_seq, num_passes=3)`
   - Implements iterative refinement loop

2. **Pass 1 - Screening**:
   - Fast prediction with β = -0.5 (exploratory)
   - Run for reduced diffusion steps (20% of normal)
   - Compute contact scores

3. **Pass 2 - Stratified Refinement**:
   - Sort Pass 1 results by contact score
   - Top 30%: Run with β = +0.5 (optimize)
   - Middle 40%: Run with β = 0 (neutral refinement)
   - Bottom 30%: Run with β = -0.75 (radical exploration)
   - Full diffusion steps

4. **Pass 3 - Selection**:
   - Combine all Pass 2 results
   - Score by: contact + confidence (pLDDT)
   - Select top-K for final output

5. **Convergence Check**:
   - Compare ensemble diversity across passes
   - If diversity decreases significantly, stop early
   - Convergence metric: pairwise RMSD of final structures

6. **Configuration**:
   - CLI flag: `--contact_guided_feedback` (boolean)
   - Parameter: `--feedback_passes` (default: 3)
   - Parameter: `--pass1_diffusion_fraction` (default: 0.2)
   - Parameter: `--stratification_thresholds` (default: [0.3, 0.7] for [top, bottom])

**Pseudocode**:
```python
def contact_guided_feedback(ab_seq, ag_seq, config):
    # Pass 1: Screening with exploratory β
    pass1_results = []
    for seed in range(5):  # 5 diverse samples
        pred = boltz_predict_with_beta(
            ab_seq, ag_seq,
            beta=-0.5,
            diffusion_steps=int(0.2 * T)  # Only 20% of steps
        )
        contact = compute_contact_score(pred)
        pass1_results.append((pred, contact))

    # Stratify by contact quality
    sorted_results = sorted(pass1_results, key=lambda x: x[1])
    split_point_low = int(0.3 * len(sorted_results))
    split_point_high = int(0.7 * len(sorted_results))

    pass2_results = []

    # Top 30%: Optimize
    for pred, contact in sorted_results[split_point_high:]:
        pred2 = refine_with_beta(pred, ab_seq, ag_seq, beta=+0.5, full_steps=True)
        pass2_results.append(pred2)

    # Middle 40%: Neutral refinement
    for pred, contact in sorted_results[split_point_low:split_point_high]:
        pred2 = refine_with_beta(pred, ab_seq, ag_seq, beta=0.0, full_steps=True)
        pass2_results.append(pred2)

    # Bottom 30%: Radical exploration
    for pred, contact in sorted_results[:split_point_low]:
        pred2 = refine_with_beta(pred, ab_seq, ag_seq, beta=-0.75, full_steps=True)
        pass2_results.append(pred2)

    # Pass 3: Select best
    final_ensemble = sorted(
        pass2_results,
        key=lambda p: compute_contact_score(p) + 0.5 * p.mean_pLDDT
    )[:5]

    return final_ensemble
```

---

## Part 3: Novel Ideas (Q through U)

### Idea Q - Iterative Epitope Refinement

**Combines**: FK particle resampling + β-scaling + contact accumulation

**Code Changes Needed**:

1. **In iterative_epitope_refinement.py** (new module):
   - Implements 3-round algorithm

2. **Round 1 - Broad Exploration**:
   - Generate predictions with weak contact guidance (weight=0.1)
   - Use exploratory β (±0.3) across CDR3
   - Use FK resampling with low pressure (soft weighting)
   - Output: Initial epitope heatmap

3. **Round 2 - Targeted Refinement**:
   - Identify hotspot regions from Round 1 (top 20% contact)
   - Apply region-specific β-scaling focused on hotspots
   - Use FK resampling with stronger contact guidance (weight=0.5)
   - Output: Refined epitope estimate

4. **Round 3 - Validation & High-Resolution**:
   - Run final predictions optimized for epitope regions
   - Generate final structures with validated contacts
   - Output: High-confidence complex structures + epitope map

5. **Configuration**:
   - CLI flag: `--iterative_epitope_refinement` (boolean)
   - Parameter: `--epitope_rounds` (default: 3)
   - Parameter: `--epitope_hotspot_percentile` (default: 0.2 = top 20%)
   - Parameter: `--round1_guidance_weight` (default: 0.1)
   - Parameter: `--round2_guidance_weight` (default: 0.5)

---

### Idea T - MSA-Free Steering Strategies

**Target**: Synthetic/designed antibodies with limited evolutionary information

**Code Changes Needed**:

1. **In msa_free_steering.py** (new module):
   - Detect MSA quality: num_sequences, sequence_diversity
   - If MSA_quality < threshold: activate MSA-free mode

2. **MSA-Free Contact Scoring**:
   - Use pAE-derived distance predictions (already computed by Boltz2)
   - Don't rely on explicit distance calculations
   - Contact defined as: pAE < threshold (e.g., 5 Å)

3. **Intrinsic β-Scaling Priors**:
   - Use larger β ranges (±0.75 instead of ±0.5) to activate internalized priors
   - More aggressive exploration to compensate for lack of MSA

4. **Confidence Weighting**:
   - Use pLDDT directly (not derived from MSA quality)
   - Ensemble selection based on internal confidence only
   - No external epitope predictions

5. **Configuration**:
   - CLI flag: `--msa_free_mode` (boolean, auto-detect or manual override)
   - Parameter: `--msa_quality_threshold` (default: 10 sequences for "good" MSA)
   - Parameter: `--msa_free_beta_scale` (default: 1.5 = 50% larger β)
   - Parameter: `--pae_contact_threshold` (default: 5.0 Å)

---

### Idea U - Energy Landscape Mapping

**Code Changes Needed**:

1. **In landscape_mapper.py** (new module):
   - Sweep β from -1.0 to +1.0 in steps of 0.25 (9 values)
   - Generate 9 predictions, one per β value
   - Collect all conformations (not just best)

2. **Landscape Characterization**:
   - Plot contact_score vs. β value
   - Identify local maxima (optimal β regions)
   - Identify plateaus (insensitive regions)
   - Identify discontinuities (phase transitions)

3. **Conformational Clustering**:
   - Cluster predictions by structural similarity (pairwise RMSD)
   - For each cluster: report average contact score, size, β range

4. **Output**: Visualization of energy landscape
   - Graph: β value on x-axis, contact score on y-axis
   - Annotations: cluster IDs, structure similarity
   - Recommended β ranges for different objectives

5. **Configuration**:
   - CLI flag: `--map_energy_landscape` (boolean)
   - Parameter: `--landscape_beta_range` (default: [-1.0, +1.0])
   - Parameter: `--landscape_beta_steps` (default: 9)
   - Parameter: `--landscape_clustering_threshold` (default: 2.0 Å RMSD)

---

## Part 4: Testing & Validation Infrastructure

### Required Testing Framework

```python
# test_steering_improvements.py

class TestSteeringImprovements:
    """Test suite for all improvements"""

    @pytest.fixture
    def known_complex(self):
        """Load known antibody-antigen complex for testing"""
        # PDB ID: 1BUN, 1A2Y, etc.
        pass

    def test_efficiency_improvement(self, improvement_name):
        """Check that improved method is faster than original"""
        time_original = run_original_strategy()
        time_improved = run_improved_strategy()
        assert time_improved < time_original * 0.9  # 10% improvement minimum

    def test_quality_improvement(self, improvement_name):
        """Check that improved method produces better predictions"""
        pred_original = run_original_strategy()
        pred_improved = run_improved_strategy()

        # Contact score should improve
        assert pred_improved.contact_score >= pred_original.contact_score

        # Structure quality should not degrade
        assert pred_improved.mean_pLDDT >= pred_original.mean_pLDDT * 0.95

    def test_epitope_discovery(self):
        """Test epitope identification accuracy"""
        contact_heatmap = epitope_scanning(ab_seq, ag_seq)
        known_epitope = get_known_epitope()

        precision, recall = evaluate_epitope_prediction(contact_heatmap, known_epitope)
        assert precision > 0.7 and recall > 0.6

    def test_ensemble_diversity(self):
        """Test ensemble generation produces diverse structures"""
        ensemble = generate_ensemble_with_beta_sweep()

        # Compute pairwise RMSD
        for i in range(len(ensemble)):
            for j in range(i+1, len(ensemble)):
                rmsd = compute_rmsd(ensemble[i], ensemble[j])
                assert rmsd > 1.0  # Minimum diversity
```

---

## Part 5: Code Organization & Module Structure

### Recommended Directory Structure

```
boltz/
├── steering/
│   ├── __init__.py
│   ├── core/
│   │   ├── beta_scaling.py          # Core β-scaling functions
│   │   ├── contact_scoring.py       # Contact computation
│   │   ├── potentials.py            # Existing potential functions (CDR3, antigen)
│   │   └── fk_resampling.py         # FK particle resampling
│   │
│   ├── strategies/
│   │   ├── cdr3_beta_scaling.py    # Idea L+ implementation
│   │   ├── region_specific_beta.py # Strategies E+, K+
│   │   ├── adaptive_phases.py      # Strategies G+, N+
│   │   ├── asymmetric_cdr.py       # Strategies J+, O+
│   │   ├── iterative_epitope.py    # Idea Q
│   │   ├── multi_objective.py      # Idea R
│   │   ├── template_guided.py      # Idea S
│   │   ├── msa_free_steering.py    # Idea T
│   │   └── landscape_mapper.py     # Idea U
│   │
│   ├── utils/
│   │   ├── antigen_partitioning.py
│   │   ├── interface_metrics.py
│   │   ├── beam_search.py
│   │   └── visualization.py
│   │
│   └── tests/
│       ├── test_beta_scaling.py
│       ├── test_strategies.py
│       └── test_integration.py
```

---

## Part 6: Integration Points with Existing Code

### Inference Loop Integration

**Location**: `boltz/inference.py` or `boltz/inference_loop.py`

**Changes**:
- Add conditional: if `config.use_beta_scaling`: apply β scaling in pair representation
- Hook into Pairformer input (after MSA processing, before attention computation)
- Minimal code changes: 5-10 lines per integration

### Configuration System

**Location**: `boltz/config.py` or argument parser

**New Configuration Sections**:
```python
@dataclass
class BetaScalingConfig:
    enabled: bool = False
    cdr3_beta: float = 0.0
    cdr3_beta_sweep: bool = False
    interface_beta: float = 0.0
    region_specific_config: Dict = None

@dataclass
class SteeringConfig:
    use_potentials: bool = False
    use_beta_scaling: bool = False
    beta_scaling: BetaScalingConfig = None
    adaptive_phases: bool = False
    iterative_refinement: bool = False
```

### CLI Flag Integration

**Location**: Argument parser (typically `boltz/cli.py` or `bin/boltz`)

**New flags** (with suggested names):
- `--use_beta_scaling` (boolean)
- `--cdr3_beta_value` (float)
- `--cdr3_beta_sweep` (boolean)
- `--interface_beta` (float)
- `--adaptive_phases` (boolean)
- `--epitope_scanning` (boolean)
- `--iterative_epitope_refinement` (boolean)
- `--coupled_cdr3_steering` (boolean)
- `--energy_landscape_map` (boolean)
- `--msa_free_mode` (boolean)

---

## Part 7: Performance Considerations

### Computational Overhead by Improvement

| Improvement | Operation | Overhead | Mitigation |
|-------------|-----------|----------|-----------|
| β-scaling | Pair masking + multiplication | ~1-2% | Vectorized ops, GPU kernels |
| Region-specific β | Additional masks | ~2-5% | Precompute masks once |
| Contact computation at each step | Distance calculation | ~10-20% | Use pAE-based contacts (faster), or compute every N steps |
| FK resampling | Particle weighting | ~5-10% | Already done in particle methods |
| Epitope scanning | N separate passes | ~N× | Mitigated by single pass with β-sweep (10-20x speedup) |

### Recommended Optimizations

1. **GPU-Accelerated β-scaling**:
   - Use torch/TensorFlow operations for masking
   - Fuse with existing pair representation computation

2. **Cached Contact Scores**:
   - Precompute distance matrices (CDR atoms × antigen atoms)
   - Update incrementally rather than full recomputation

3. **Parallel β-Sweeps**:
   - For ensemble generation, run different β values in parallel
   - Separate GPU instances per β value (if resources allow)

---

## Summary: Implementation Order Recommendation

**Week 1**: Implement L+ (CDR3 β-scaling) → Simplest, immediate integration
**Week 2**: Implement K+ (Region-specific) + N+ (Adaptive phases) → Build on L+
**Week 3**: Implement Q (Iterative epitope) → Novel contribution, uses K+
**Week 4**: Implement T (MSA-free) + O+ (Asymmetric β) → Extend to more cases
**Week 5+**: Implement R, S, U → Advanced research features

---

---

## Part 8: Strategy V - Embedding-Space CDR3 Steering (EmbedOpt-Inspired)

### Overview

Strategy V steers CDR3 loop conformations by optimizing the Pairformer's pair representation `z` in embedding space, rather than applying coordinate-space potentials (like the existing `CDR3ConformationPotential`). Inspired by the EmbedOpt paper (Li et al., 2602.05285).

**Key Advantages**:
- More robust for novel CDR sequences not in training data
- Stable across 2 orders of magnitude of hyperparameter variation (vs. brittle coordinate methods)
- Captures biological signals more naturally in learned representation space
- Does NOT require `--use_potentials` (operates in embedding space, not coordinate space)
- Compatible with existing `--cdr3_steering` and `--antigen_steering`

### Integration Point

**Location**: `src/boltz/model/modules/diffusion_conditioning.py`, in the `DiffusionConditioning.forward()` method, after CDR3 beta scaling (line ~118) and before `atom_encoder` call (line ~120).

This is where the pair representation `z` is available and modifiable before it's projected into biases for the atom encoder/decoder and token transformer.

### Mechanism

1. Extract CDR3 pair mask from `feats["embedding_steering_mask"]`
2. Compute target embedding:
   - `self_reference` mode: mean of non-CDR3 valid pairs (pushes CDR3 toward well-ordered patterns learned by the model)
3. Run N gradient descent steps with analytical gradient:
   - Loss: `||z_cdr3 - z_target||^2`
   - Gradient: `2 * (z_cdr3 - z_target)` (closed-form, no autograd needed)
   - Learning rate schedule: `strength * 0.1 * (1 - step/N)` (decreasing)
4. Modified `z` flows into `atom_encoder` and `token_trans_bias`, affecting all 200 denoising steps

### Code Changes Needed

#### 1. `src/boltz/data/types.py` (~line 545)
Add `embedding_steering_constraints` field to `InferenceOptions`:
```python
embedding_steering_constraints: Optional[
    list[tuple[str, float, int, list[tuple[int, int, int]]]]
] = None
```
Format: `(mode, strength, num_opt_steps, [(chain_id, start_res, end_res), ...])`

#### 2. `src/boltz/data/parse/schema.py` (~line 1670)
Add `embedding_steering` constraint parsing alongside existing `cdr3_beta_scaling`:
- Initialize `embedding_steering_constraints = []` at ~line 1517
- Add `elif "embedding_steering" in constraint:` block
- Parse: mode, strength, num_opt_steps, cdr3_regions
- Pass to `InferenceOptions` constructor at ~line 1920

#### 3. `src/boltz/data/feature/featurizerv2.py` (~line 2403)
Add `process_embedding_steering_constraints()` function:
- Creates `embedding_steering_mask` (boolean token mask for CDR3 residues)
- Creates `embedding_steering_strength` (scalar)
- Creates `embedding_steering_num_opt_steps` (scalar)
- Creates `embedding_steering_mode` (integer: 0=self_reference)

Update `Boltz2Featurizer.process()`:
- Add `inference_embedding_steering_constraints` parameter
- Call `process_embedding_steering_constraints()` and merge into return dict

#### 4. `src/boltz/data/module/inferencev2.py` (~line 264)
Unpack `embedding_steering_constraints` from `options` and pass to featurizer.

#### 5. `src/boltz/main.py`
- Add `embedding_steering: bool = False` to `BoltzSteeringParams` (~line 175)
- Add `--embedding_steering` CLI flag (~line 1003)
- Add to predict function signature and steering_args setup

#### 6. `src/boltz/model/modules/diffusion_conditioning.py` (THE CORE CHANGE)
Insert embedding-space optimization after existing CDR3 beta scaling (line 118) and before `atom_encoder` call (line 120).

### YAML Constraint Format

```yaml
constraints:
  - embedding_steering:
      mode: self_reference     # self_reference | canonical
      strength: 1.0            # gradient step size multiplier
      num_opt_steps: 10        # number of optimization steps
      cdr3_regions:
        - chain: B
          start_res: 97
          end_res: 115
        - chain: C
          start_res: 88
          end_res: 100
```

### Key Design Decisions

1. **Standalone flag**: `--embedding_steering` does NOT require `--use_potentials`
2. **Analytical gradients**: Closed-form MSE gradient avoids torch.autograd overhead during inference
3. **Decreasing learning rate**: `lr = strength * 0.1 * (1 - step/N)` ensures convergence
4. **Self-reference target**: Uses mean of non-CDR3 well-ordered pairs as target, leveraging the model's learned representations

### Testing

```bash
# Embedding steering on 7TRH_HBG
conda activate boltz
boltz predict examples/7TRH/7TRH_HBG_embedding_steer.yml \
    --embedding_steering --out_dir ./output_embed_steer --recycling_steps 3 --sampling_steps 50

# Baseline comparison
boltz predict examples/7TRH/7TRH_HBG.yml --out_dir ./output_baseline --recycling_steps 3 --sampling_steps 50
```

---

**Status**: Ready for implementation
**Last Updated**: 2025-02-23
