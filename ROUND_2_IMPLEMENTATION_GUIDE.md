# Round 2 Implementation Guide: New Strategies N1-N15

**Date**: 2026-03-16
**Purpose**: Concrete implementation suggestions for each new strategy from `ROUND_1_FINAL_EVALUATION.md` (without making code changes yet)
**Scope**: Code structure, module changes, integration points
**Codebase root**: `/mnt/c/Users/octav/Documents/claude/proteinEBM/boltz/src/boltz/`

---

## Architecture Overview

All Round 1 strategies follow the same pattern of modifying a small set of core files:

| File | Role | What Gets Modified |
|------|------|--------------------|
| `main.py` (lines 165-178) | CLI flags + `BoltzSteeringParams` dataclass | New flags and steering params |
| `model/potentials/potentials.py` | Potential classes + `get_potentials()` factory | New potential classes, factory logic |
| `model/potentials/schedules.py` | Parameter scheduling | New schedule types |
| `model/modules/diffusionv2.py` (`AtomDiffusion.sample()`, lines 295-505) | Diffusion loop with FK resampling + gradient guidance | Loop logic, phase transitions, re-ranking |
| `model/modules/diffusion_conditioning.py` (`DiffusionConditioning.forward()`, lines 83-139) | Beta-scaling on attention biases | New mask types, scaling logic |
| `data/parse/schema.py` (lines 1514-1707) | YAML constraint parsing | New constraint blocks |
| `data/feature/featurizerv2.py` (lines 2170-2403) | Feature tensor generation from constraints | New `process_*()` methods |
| `data/types.py` | `InferenceOptions` fields | New constraint fields |
| `data/module/inferencev2.py` | Constraint wiring to featurizer | New data paths |
| `model/models/boltz2.py` | `predict_step()`, passes `z_trunk` and steering_args | New steering orchestration |

Each new strategy should be developed in its own **git worktree** (e.g., `worktree_steering_cdr3/n1_fk_hierarchical/`) branching from the common base.

---

## Tier 1: High Impact, Low-Medium Complexity

---

### N1: A+Y Hybrid (FK Particles + Hierarchical Timing)

**Goal**: Combine A's FK resampling (best DockQ) with Y's hierarchical timing (best CDR-H3 RMSD).

**Architecture**: Two-phase diffusion with FK resampling in the late phase only. Early phase uses beta-scaling for global orientation; late phase activates FK particles with coordinate potentials for refinement.

#### Files to Modify

**1. `main.py`** — Add CLI flags and steering params

```python
# In BoltzSteeringParams (line 165):
@dataclass
class BoltzSteeringParams:
    # ... existing fields ...
    hybrid_fk_hierarchical: bool = False
    hybrid_transition_fraction: float = 0.5    # Switch phase at 50% of steps
    hybrid_early_beta: float = 0.3             # Beta-scaling in early phase
    hybrid_late_particles: int = 5             # FK particles in late phase

# In click options (~line 1004):
@click.option("--hybrid_fk_hierarchical", is_flag=True,
    help="Enable A+Y hybrid: beta-scaling early + FK resampling late.")
@click.option("--hybrid_transition_fraction", type=float, default=0.5)
@click.option("--hybrid_early_beta", type=float, default=0.3)
@click.option("--hybrid_late_particles", type=int, default=5)
```

**2. `model/modules/diffusionv2.py`** — Two-phase logic in `AtomDiffusion.sample()`

The key change is in the main denoising loop (line 350). Currently it applies FK resampling + gradient guidance uniformly across all steps. For N1, we split the loop into two phases:

```python
# Inside sample() method, after sigmas_and_gammas are computed (line 337):

if steering_args.get("hybrid_fk_hierarchical", False):
    transition_step = int(num_sampling_steps * steering_args["hybrid_transition_fraction"])
else:
    transition_step = 0  # No phase transition, all steps are "late phase"

# In the main loop (line 350):
for step_idx, (sigma_tm, sigma_t, gamma) in enumerate(sigmas_and_gammas):
    is_early_phase = steering_args.get("hybrid_fk_hierarchical", False) and step_idx < transition_step

    # ... existing random augmentation code (lines 351-377) ...

    # EARLY PHASE: Beta-scaling only (no FK, no gradient guidance)
    if is_early_phase:
        # Beta-scaling is handled in DiffusionConditioning.forward() via feats
        # Just do standard denoising — no FK energy, no gradient update
        # Skip FK resampling block (lines 398-441)
        # Skip gradient guidance block (lines 443-473)
        # Skip particle resampling block (lines 475-505)
        pass  # Standard denoising step only

    # LATE PHASE: FK resampling + gradient guidance (existing code)
    else:
        # Activate FK particles at transition point
        if step_idx == transition_step and steering_args.get("hybrid_fk_hierarchical"):
            # Replicate current atom_coords into num_particles copies
            multiplicity = multiplicity * steering_args["hybrid_late_particles"]
            atom_coords = atom_coords.repeat_interleave(steering_args["hybrid_late_particles"], 0)
            atom_mask = atom_mask.repeat_interleave(steering_args["hybrid_late_particles"], 0)
            # Initialize FK tracking tensors
            energy_traj = torch.empty((multiplicity, 0), device=self.device)
            resample_weights = torch.ones(multiplicity, device=self.device).reshape(
                -1, steering_args["hybrid_late_particles"])

        # ... existing FK + gradient guidance code (lines 398-505) ...
```

**3. `model/modules/diffusion_conditioning.py`** — Phase-aware beta-scaling

The existing CDR3 beta-scaling (line 134-137) applies uniformly. For N1, the beta value should be non-zero only during the early phase:

```python
# In DiffusionConditioning.forward() (line 83):
# The feats dict should carry a "current_phase" or "diffusion_step_fraction" signal
# Option A: Pass step_fraction via feats from the diffusion loop
# Option B: Use a time-dependent beta schedule

# Modify line 137:
if cdr3_pair_mask_4d is not None:
    # Apply beta only during early phase (step_fraction > transition)
    step_fraction = feats.get("diffusion_step_fraction", torch.tensor(1.0))
    transition = feats.get("hybrid_transition_fraction", torch.tensor(0.5))
    # Smooth ramp: full beta early, zero late
    phase_weight = torch.clamp((step_fraction - transition) / 0.1, 0.0, 1.0)
    effective_beta = beta_val * (1.0 - phase_weight)
    token_trans_bias = token_trans_bias + effective_beta * cdr3_pair_mask_4d * token_trans_bias
```

To pass `diffusion_step_fraction` into feats, modify the diffusion loop in `diffusionv2.py` to inject it before the network forward call (line 388):

```python
network_condition_kwargs["feats"]["diffusion_step_fraction"] = torch.tensor(
    1.0 - step_idx / num_sampling_steps, device=self.device
)
```

**4. YAML / Featurizer** — Reuse existing `cdr3_beta_scaling` constraint format with the hybrid flag. No new YAML parsing needed; the hybrid behavior is controlled entirely by CLI flags.

#### Data Flow

```
CLI: --hybrid_fk_hierarchical --hybrid_early_beta 0.3 --hybrid_late_particles 5
  ↓
BoltzSteeringParams → sample() in diffusionv2.py
  ↓
Step 0 to transition_step:
  DiffusionConditioning applies beta=0.3 to CDR3 attention biases
  No FK resampling, no gradient guidance
  ↓
Step transition_step:
  Replicate atom_coords into 5 particles
  Initialize FK tracking tensors
  ↓
Step transition_step to end:
  FK resampling + gradient guidance (AntigenOrientationPotential, etc.)
  Standard potentials applied with full weights
  Final resampling selects best particle
```

#### Testing

```python
def test_n1_hybrid():
    # Phase 1: should run without FK overhead
    # Phase 2: should produce multiple particles
    # Final output: single best structure from FK selection
    # Check: CDR-H3 RMSD should improve (Y's contribution)
    # Check: DockQ should improve (A's contribution)
    # Check: No OOM for 5 particles (vs A's 20 that caused OOM)
```

#### Estimated Complexity: Medium (main logic is in diffusionv2.py loop refactoring)

---

### N2: Q→B2 Pipeline (Iterative Epitope → Contact Restraints)

**Goal**: Use Q's iterative epitope discovery output as input to B2-style contact restraints, achieving B2-level performance without oracle info.

**Architecture**: Two-stage external pipeline. Stage 1 runs Q to predict epitope residues. Stage 2 generates new YAML files with contact restraints targeting those residues, then runs standard Boltz2 with B2-style restraints.

#### Implementation Approach

This is **not a model code change** — it's a **pipeline script** that orchestrates two Boltz2 runs. This makes it the simplest to implement.

**1. New script: `scripts/pipeline_q_to_b2.py`**

```python
"""
Pipeline: Q (iterative epitope refinement) → B2 (contact restraints)

Stage 1: Run Method Q to identify candidate epitope residues
Stage 2: Convert epitope predictions to contact restraints YAML
Stage 3: Run Boltz2 with generated contact restraints
"""

def stage1_run_method_q(input_yaml_dir, output_dir, boltz_cmd):
    """Run Method Q predictions on all complexes."""
    # For each YAML in input_yaml_dir:
    #   boltz predict <yaml> --use_potentials --iterative_epitope_refinement \
    #     --diffusion_samples 5 --output <output_dir>/q_predictions/
    pass

def stage2_extract_epitope(q_predictions_dir, gt_threshold=5.0):
    """Extract predicted epitope residues from Q's output structures.

    For each complex:
    1. Load the best-by-confidence predicted structure
    2. Compute CDR-antigen contacts (heavy-atom distance < threshold)
    3. Identify antigen residues with contacts → predicted epitope
    4. Return dict: {complex_name: [epitope_residue_ids]}
    """
    pass

def stage3_generate_contact_yaml(original_yaml_dir, epitope_predictions, output_yaml_dir):
    """Generate B2-style contact restraint YAML files.

    For each complex:
    1. Read the original YAML (sequences, chain definitions)
    2. Add contact constraints between CDR residues and predicted epitope residues
    3. Generate multiple constraint variants:
       - hbond restraints on predicted epitope
       - hydrophobic restraints on predicted epitope
       - salt_bridge restraints on predicted epitope
    4. Write to output_yaml_dir/
    """
    # YAML format for contact constraints (matching existing B2 format):
    # constraints:
    #   - contact:
    #       atoms: [[chain_A, res_id, atom_name], [chain_B, res_id, atom_name]]
    #       max_distance: 5.0
    pass

def stage4_run_with_restraints(restraint_yaml_dir, output_dir, boltz_cmd):
    """Run Boltz2 with generated contact restraints."""
    # boltz predict <yaml> --use_potentials --diffusion_samples 5 --output <output_dir>/
    pass

def main():
    stage1_run_method_q(...)
    epitope_preds = stage2_extract_epitope(...)
    stage3_generate_contact_yaml(..., epitope_preds, ...)
    stage4_run_with_restraints(...)
```

**2. Contact restraint YAML generation**

The key function is `stage3_generate_contact_yaml`. It needs to produce YAML in the exact format that `schema.py` (lines 1565-1597) expects for `contact_constraints`:

```yaml
# Generated YAML for complex 8FAH_HLA
version: 2
sequences:
  # ... original sequence definitions ...
constraints:
  - contact:
      connections:
        - [A, 42, CA]   # Predicted epitope residue
        - [B, 103, CA]  # CDR-H3 residue
      max_distance: 8.0
  - contact:
      connections:
        - [A, 45, CA]
        - [B, 105, CA]
      max_distance: 8.0
  # ... more predicted contacts ...
```

**3. No model code changes needed** — this leverages existing B2 infrastructure entirely.

#### Data Flow

```
Input: 47 YAML configs (no oracle info)
  ↓
Stage 1: Run Q (existing --iterative_epitope_refinement flag)
  → 47 predicted structures with epitope contacts
  ↓
Stage 2: Extract epitope residues from predicted structures
  → Dict: {complex: [residue_ids]}
  ↓
Stage 3: Generate B2-style contact restraint YAMLs
  → 47 × 3 new YAML files (hbond, hydrophobic, salt_bridge variants)
  ↓
Stage 4: Run Boltz2 with restraints (existing --use_potentials)
  → Final predictions with contact-guided docking
```

#### Risk Analysis

The pipeline's success depends entirely on Q's epitope prediction accuracy. From Round 1:
- Q's epitope F1 = 0.381 (mean), 0.373 (median)
- If Q predicts the wrong epitope, Stage 4 will dock to the wrong region
- **Mitigation**: Generate restraints at multiple confidence thresholds (top 10%, 20%, 30% of contacts), run all three, and select the best by Boltz2 confidence

#### Estimated Complexity: Low-Medium (pipeline script, no model changes)

---

### N3: Confidence-Aware Conditional Steering

**Goal**: Only apply steering when the baseline confidence is low, avoiding regressions on easy cases.

**Architecture**: Two-pass approach. Pass 1 runs vanilla Boltz2 and checks confidence. Pass 2 applies steering only if confidence is below a threshold.

#### Implementation Approach A: External Script (Simplest)

**New script: `scripts/conditional_steering.py`**

```python
"""
Conditional Steering Pipeline:
Pass 1: Run vanilla Boltz2, check confidence (iptm)
Pass 2: If iptm < threshold, run with steering; else keep Pass 1 result
"""

def run_conditional_steering(yaml_dir, output_dir, steering_method, threshold_high=0.85, threshold_low=0.75):
    for yaml_file in yaml_dir.glob("*.yml"):
        # Pass 1: Vanilla
        result_vanilla = boltz_predict(yaml_file, use_potentials=False, diffusion_samples=5)
        iptm = result_vanilla.confidence_scores["iptm"]

        if iptm >= threshold_high:
            # High confidence → keep vanilla result, skip steering
            copy_result(result_vanilla, output_dir)
        elif iptm >= threshold_low:
            # Medium confidence → apply mild steering (Y-style hierarchical)
            result_steered = boltz_predict(yaml_file,
                use_potentials=True, hierarchical_steering=True)
            # Select better result by confidence
            select_better(result_vanilla, result_steered, output_dir)
        else:
            # Low confidence → apply aggressive steering (G+-style progressive)
            result_steered = boltz_predict(yaml_file,
                use_potentials=True, progressive_cdr_refinement=True)
            select_better(result_vanilla, result_steered, output_dir)
```

#### Implementation Approach B: Single-Pass with Internal Conditioning (More Efficient)

Modify the diffusion loop to check confidence after an initial "assessment" period and then decide whether to activate steering.

**1. `model/modules/diffusionv2.py`** — Add confidence assessment at a checkpoint

```python
# Inside sample() method, after the first ~20% of denoising steps:
ASSESSMENT_FRACTION = 0.2
assessment_step = int(num_sampling_steps * ASSESSMENT_FRACTION)

for step_idx, (sigma_tm, sigma_t, gamma) in enumerate(sigmas_and_gammas):
    # ... standard denoising ...

    if step_idx == assessment_step and steering_args.get("conditional_steering"):
        # Quick confidence estimate from current denoised coords
        # Use iptm proxy: compute CDR-antigen contact count
        contact_count = compute_quick_contact_count(
            atom_coords_denoised, network_condition_kwargs["feats"])

        if contact_count > steering_args["conditional_contact_threshold"]:
            # Good trajectory → disable steering for remaining steps
            steering_active = False
        else:
            # Poor trajectory → activate steering
            steering_active = True

    # Apply steering only if active
    if steering_active:
        # ... FK resampling + gradient guidance (existing code) ...
```

**2. `main.py`** — Add CLI flags

```python
# In BoltzSteeringParams:
conditional_steering: bool = False
conditional_contact_threshold: int = 5  # Min CDR-antigen contacts to skip steering

# CLI options:
@click.option("--conditional_steering", is_flag=True,
    help="Only apply steering when initial prediction has low contact quality.")
```

**3. Helper function for quick contact assessment** — Add to `potentials.py` or a new utility:

```python
def compute_quick_contact_count(atom_coords, feats, threshold=8.0):
    """Count CDR-antigen atom pairs within threshold distance.
    Uses existing antigen_atom_index and cdr_atom_index from feats.
    Returns: int count of contacts.
    """
    if "antigen_atom_index" not in feats or "cdr_atom_index" not in feats:
        return 0
    ag_idx = feats["antigen_atom_index"][0]
    cdr_idx = feats["cdr_atom_index"][0]
    ag_coords = atom_coords[0, ag_idx]   # [N_ag, 3]
    cdr_coords = atom_coords[0, cdr_idx] # [N_cdr, 3]
    dists = torch.cdist(ag_coords, cdr_coords)  # [N_ag, N_cdr]
    return (dists < threshold).sum().item()
```

#### Estimated Complexity: Low (Approach A) / Medium (Approach B)

---

### N4: L-light (CDR3 Beta, Reduced Strength + Time-Dependent)

**Goal**: L's diversity benefit with reduced accuracy penalty. Lower beta (0.1-0.15) + time-dependent schedule + CDR3-antigen pair scaling.

**Architecture**: Modifies the existing CDR3 beta-scaling in `DiffusionConditioning` to use a time-dependent schedule and expanded mask.

#### Files to Modify

**1. `data/parse/schema.py`** — Extend `cdr3_beta_scaling` constraint

Add new fields to the existing CDR3 beta-scaling constraint block (lines 1670-1703):

```python
# In parse_cdr3_beta_scaling() or equivalent:
# New optional fields:
#   beta_schedule: "constant" | "time_dependent" (default: "constant")
#   beta_early: float (beta value for early diffusion steps)
#   beta_late: float (beta value for late diffusion steps, default: 0.0)
#   include_antigen_pairs: bool (scale CDR3-antigen pairs too, default: false)
```

**2. `data/feature/featurizerv2.py`** — Expand mask generation in `process_cdr3_beta_constraints()` (line 2357)

```python
def process_cdr3_beta_constraints(self, ...):
    # Existing: creates cdr3_token_mask for CDR3-CDR3 pairs
    # New: if include_antigen_pairs, also create cdr3_antigen_token_mask
    #   cdr3_antigen_mask[i,j] = True if (i in CDR3 and j in antigen) or vice versa
    # This mask targets the CDR3-antigen interface pairs for scaling

    feats["cdr3_token_mask"] = cdr3_mask           # [n_tokens]
    feats["cdr3_beta_value"] = beta_value           # scalar
    feats["cdr3_beta_early"] = beta_early           # scalar (new)
    feats["cdr3_beta_late"] = beta_late             # scalar (new)
    feats["cdr3_beta_schedule"] = schedule_type     # str (new)
    feats["cdr3_antigen_token_mask"] = ag_mask      # [n_tokens] (new, optional)
```

**3. `model/modules/diffusion_conditioning.py`** — Time-dependent beta + expanded mask

Replace the fixed beta application (line 134-137) with a schedule-aware version:

```python
# After line 98 (existing beta extraction):
if "cdr3_token_mask" in feats and "cdr3_beta_value" in feats:
    cdr3_mask = feats["cdr3_token_mask"].to(z.device).to(torch.bool)
    beta_val = float(feats["cdr3_beta_value"].item())

    # NEW: Time-dependent beta schedule
    schedule = feats.get("cdr3_beta_schedule", "constant")
    if schedule == "time_dependent":
        beta_early = float(feats.get("cdr3_beta_early", torch.tensor(beta_val)).item())
        beta_late = float(feats.get("cdr3_beta_late", torch.tensor(0.0)).item())
        # step_fraction comes from diffusion loop (1.0=start, 0.0=end)
        t = float(feats.get("diffusion_step_fraction", torch.tensor(0.5)).item())
        # Linear interpolation: beta_early at t=1.0, beta_late at t=0.0
        beta_val = beta_early * t + beta_late * (1.0 - t)

    if abs(beta_val) > 1e-6:
        # CDR3-CDR3 pair mask (existing)
        cdr3_mask_i = cdr3_mask.unsqueeze(-1)
        cdr3_mask_j = cdr3_mask.unsqueeze(-2)
        cdr3_pair_mask = (cdr3_mask_i & cdr3_mask_j)

        # NEW: CDR3-Antigen pair mask (optional expansion)
        if "cdr3_antigen_token_mask" in feats:
            ag_mask = feats["cdr3_antigen_token_mask"].to(z.device).to(torch.bool)
            ag_mask_i = ag_mask.unsqueeze(-1)
            ag_mask_j = ag_mask.unsqueeze(-2)
            # CDR3-antigen pairs: CDR3_i × antigen_j OR antigen_i × CDR3_j
            cross_mask = (cdr3_mask_i & ag_mask_j) | (ag_mask_i & cdr3_mask_j)
            combined_mask = cdr3_pair_mask | cross_mask
        else:
            combined_mask = cdr3_pair_mask

        cdr3_pair_mask_4d = combined_mask.unsqueeze(0).unsqueeze(-1).float()

# Apply at line 137:
if cdr3_pair_mask_4d is not None:
    token_trans_bias = token_trans_bias + beta_val * cdr3_pair_mask_4d * token_trans_bias
```

**4. `model/modules/diffusionv2.py`** — Inject step fraction into feats

```python
# Before the network forward call (line 388), inject time info:
network_condition_kwargs["feats"]["diffusion_step_fraction"] = torch.tensor(
    1.0 - step_idx / num_sampling_steps, device=self.device
)
```

#### YAML Format

```yaml
constraints:
  - cdr3_beta_scaling:
      chains:
        - chain: B
          cdr3_start: 97
          cdr3_end: 115
        - chain: C
          cdr3_start: 88
          cdr3_end: 100
      beta: 0.15                    # Reduced from 0.3
      beta_schedule: time_dependent  # NEW
      beta_early: 0.15              # NEW: exploration in early steps
      beta_late: 0.0                # NEW: convergence in late steps
      include_antigen_pairs: true    # NEW: scale CDR3-antigen pairs too
      antigen_chain: A               # NEW: required if include_antigen_pairs
```

#### Estimated Complexity: Low (extends existing CDR3 beta infrastructure)

---

### N5: Multi-Config K Ensemble

**Goal**: Run K with 5-8 different region emphasis patterns per complex, mimicking B2's multi-config advantage.

**Architecture**: External pipeline script that generates multiple YAML configs with different CDR3 beta-scaling parameters and aggregates results.

#### Implementation

**1. New script: `scripts/k_ensemble.py`**

```python
"""
Multi-Config K Ensemble Pipeline:
Generate 5-8 YAML variants per complex with different beta configurations.
Run all variants, aggregate into epitope heatmap, select best predictions.
"""

def generate_k_configs(base_yaml, output_dir, n_configs=6):
    """Generate N variant YAML configs with different beta values.

    Configs:
    1. beta = +0.1 (mild emphasis)
    2. beta = +0.2 (moderate emphasis)
    3. beta = +0.3 (strong emphasis — original L value)
    4. beta = -0.1 (mild de-emphasis / exploration)
    5. beta = -0.2 (moderate exploration)
    6. beta = +0.15 with include_antigen_pairs=true (interface-focused)
    """
    beta_configs = [
        {"beta": 0.1, "include_antigen_pairs": False},
        {"beta": 0.2, "include_antigen_pairs": False},
        {"beta": 0.3, "include_antigen_pairs": False},
        {"beta": -0.1, "include_antigen_pairs": False},
        {"beta": -0.2, "include_antigen_pairs": False},
        {"beta": 0.15, "include_antigen_pairs": True},
    ]
    # For each config, copy base YAML and modify cdr3_beta_scaling block
    # Write to output_dir/complex_name_configN.yml

def aggregate_results(predictions_dir, n_configs):
    """For each complex:
    1. Load all N config predictions
    2. Pool all 5×N models (5 samples per config × N configs)
    3. Select best by confidence → primary output
    4. Accumulate contact heatmap across all models → epitope prediction
    """
    pass

def main():
    generate_k_configs(...)
    run_all_configs(...)  # boltz predict for each YAML variant
    aggregate_results(...)
```

**2. No model code changes needed** — leverages existing `cdr3_beta_scaling` feature with different parameter values.

**3. Evaluation enhancement** — Modify `scripts/evaluate.py` to support multi-config aggregation (similar to how B2/B3 are already handled with best-per-complex config selection).

#### Key Design Decision

The critical question is **how to select the best config per complex** without oracle information. Options:
- **Max confidence**: Pick the config whose best model has the highest `confidence_score`
- **Contact count**: Pick the config with the most CDR-antigen contacts in the predicted structure
- **Consensus**: Accumulate contacts across all configs, pick the config most consistent with consensus epitope

#### Estimated Complexity: Low (pipeline script only)

---

## Tier 2: Medium Impact, Medium Complexity

---

### N6: Custom Re-Ranking Module

**Goal**: Replace Boltz2's confidence score with a composite re-ranking metric for steered predictions.

**Architecture**: Post-prediction re-ranking that combines multiple signals. Can be applied to any steering method's output.

#### Files to Modify

**1. New module: `model/reranking.py`**

```python
"""
Custom re-ranking for steered predictions.
Replaces Boltz2's confidence_score with a composite metric.
"""

import torch
import numpy as np

class SteeringReranker:
    """Re-rank steered models using composite score."""

    def __init__(self, alpha=0.3, beta=0.4, gamma=0.3):
        """
        Weights for composite score:
        alpha: confidence_score weight
        beta: interface_quality weight (from steering energy or contacts)
        gamma: interface_plddt weight
        """
        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma

    def compute_interface_quality(self, atom_coords, feats):
        """Compute interface quality from predicted structure.

        Metrics:
        1. CDR-antigen contact count (< 8A threshold)
        2. Minimum CDR-antigen distance (lower = better)
        3. Buried surface area proxy (number of close contacts < 5A)
        """
        if "antigen_atom_index" not in feats or "cdr_atom_index" not in feats:
            return 0.0
        ag_idx = feats["antigen_atom_index"][0]
        cdr_idx = feats["cdr_atom_index"][0]
        ag_coords = atom_coords[0, ag_idx]
        cdr_coords = atom_coords[0, cdr_idx]
        dists = torch.cdist(ag_coords, cdr_coords)

        n_contacts_8A = (dists < 8.0).sum().float()
        n_contacts_5A = (dists < 5.0).sum().float()
        min_dist = dists.min()

        # Normalize to [0, 1]
        contact_score = torch.sigmoid(n_contacts_8A / 50.0 - 1.0)
        proximity_score = torch.sigmoid(5.0 - min_dist)

        return 0.6 * contact_score + 0.4 * proximity_score

    def compute_interface_plddt(self, plddt, feats):
        """Extract mean pLDDT at CDR and nearby antigen residues."""
        if "cdr3_token_mask" not in feats:
            return plddt.mean()
        cdr_mask = feats["cdr3_token_mask"][0].bool()
        return plddt[0, cdr_mask].mean()

    def rerank(self, models, feats_list):
        """Re-rank list of (atom_coords, confidence_score, plddt) tuples.

        Returns: indices sorted by composite score (best first).
        """
        scores = []
        for atom_coords, confidence, plddt, feats in zip(*models, feats_list):
            interface_quality = self.compute_interface_quality(atom_coords, feats)
            interface_plddt = self.compute_interface_plddt(plddt, feats)

            composite = (
                self.alpha * confidence +
                self.beta * interface_quality +
                self.gamma * interface_plddt
            )
            scores.append(composite)

        return torch.argsort(torch.tensor(scores), descending=True)
```

**2. Integration in `model/models/boltz2.py`** — Apply re-ranking in `predict_step()`

After generating all diffusion samples (line varies), instead of selecting by `confidence_score` alone:

```python
# After all samples are generated:
if steering_args.get("use_custom_reranking", False):
    reranker = SteeringReranker(
        alpha=steering_args.get("rerank_confidence_weight", 0.3),
        beta=steering_args.get("rerank_interface_weight", 0.4),
        gamma=steering_args.get("rerank_plddt_weight", 0.3),
    )
    best_indices = reranker.rerank(all_samples, all_feats)
    # Reorder samples by composite score
    all_samples = [all_samples[i] for i in best_indices]
```

**3. `main.py`** — CLI flags

```python
# In BoltzSteeringParams:
use_custom_reranking: bool = False
rerank_confidence_weight: float = 0.3
rerank_interface_weight: float = 0.4
rerank_plddt_weight: float = 0.3
```

#### Estimated Complexity: Medium (new module + integration in predict_step)

---

### N7: Y+ (Stronger Hierarchical)

**Goal**: Amplify Y's consistently positive trends by strengthening both phases and adding predicted epitope info.

**Architecture**: Extension of Y's existing hierarchical steering with three changes: (1) stronger beta-scaling, (2) stronger coordinate potentials, (3) epitope-guided late-phase potentials.

#### Files to Modify

The Y implementation already exists in the `y_hierarchical_steering` worktree. N7 extends it:

**1. `model/modules/diffusion_conditioning.py`** — Increase beta factor

In Y's existing implementation, the early-phase beta-scaling uses the value from the YAML constraint. For N7, increase the effective beta:

```yaml
# Y+ YAML: increase beta from Y's value (e.g., 0.3) to 0.5-0.7
constraints:
  - cdr3_beta_scaling:
      beta: 0.5   # Increased from Y's ~0.3
      beta_schedule: time_dependent
      beta_early: 0.5
      beta_late: 0.0
      include_antigen_pairs: true  # From N4
```

**2. `model/potentials/potentials.py`** — Stronger late-phase potentials

Modify the `get_potentials()` function (line 802) to accept amplification parameters:

```python
# In get_potentials(), for AntigenOrientationPotential (line 938):
if boltz2 and steering_args.get("antigen_steering", False):
    guidance_scale = steering_args.get("late_phase_guidance_scale", 1.0)
    potentials.append(
        AntigenOrientationPotential(
            parameters={
                "guidance_interval": 2,
                "guidance_weight": PiecewiseStepFunction(
                    thresholds=[0.3, 0.7],
                    values=[1.5 * guidance_scale, 1.0 * guidance_scale, 0.3 * guidance_scale]
                ),
                "resampling_weight": PiecewiseStepFunction(
                    thresholds=[0.5],
                    values=[1.0 * guidance_scale, 0.5 * guidance_scale]
                ),
                "union_lambda": ExponentialInterpolation(start=8.0, end=0.0, alpha=-2.0),
            }
        )
    )
```

**3. Epitope-guided late-phase** — New `EpitopeGuidedPotential`

This is a new potential class that uses a pre-computed set of "likely epitope" residues to provide directional guidance in the late phase:

```python
class EpitopeGuidedPotential(FlatBottomPotential, DistancePotential):
    """Guides CDR atoms toward a predicted set of epitope residues.

    Unlike AntigenOrientationPotential which targets ALL antigen residues,
    this focuses force on a subset of antigen residues predicted as epitope.

    Features expected:
    - predicted_epitope_index: [N_epitope] indices of predicted epitope atoms
    - cdr_atom_index: [N_cdr] indices of CDR atoms
    """

    def compute_args(self, feats, parameters):
        if "predicted_epitope_index" not in feats:
            # Fall back to full antigen if no epitope prediction
            return super().compute_args(feats, parameters)

        epitope_idx = feats["predicted_epitope_index"][0]
        cdr_idx = feats["cdr_atom_index"][0]
        # Create pairs between CDR and epitope residues only
        # ... (same pattern as AntigenOrientationPotential lines 754-799)
```

**4. `main.py`** — CLI flags

```python
# In BoltzSteeringParams:
late_phase_guidance_scale: float = 1.0    # 1.0 = Y original, 2.0 = doubled
predicted_epitope_residues: str = None     # Comma-separated residue IDs or file path
```

#### Estimated Complexity: Medium (extends existing Y worktree)

---

### N8: G+L Combined (Progressive + CDR3 Diversity)

**Goal**: Phase 1 uses L's CDR3 beta-scaling for diversity, Phase 3 uses A's FK resampling for optimization.

**Architecture**: Extends G's existing progressive refinement with L's beta-scaling in Phase 1 and A's FK particles in Phase 3. Essentially a three-phase version of N1.

#### Files to Modify

**1. `model/modules/diffusionv2.py`** — Three-phase loop

```python
# Phase definitions (steering_args parameters):
# Phase 1: CDR3 beta-scaling (L-style) — steps 0 to phase1_end
# Phase 2: Neutral transition — steps phase1_end to phase2_end
# Phase 3: FK resampling + gradient guidance (A-style) — steps phase2_end to end

phase1_end = int(num_sampling_steps * steering_args.get("gl_phase1_fraction", 0.33))
phase2_end = int(num_sampling_steps * steering_args.get("gl_phase2_fraction", 0.66))

for step_idx, (sigma_tm, sigma_t, gamma) in enumerate(sigmas_and_gammas):
    # Determine current phase
    if step_idx < phase1_end:
        current_phase = "exploration"
        # CDR3 beta active (via feats), no FK, no gradient guidance
    elif step_idx < phase2_end:
        current_phase = "transition"
        # CDR3 beta ramps to 0, no FK, mild gradient guidance
    else:
        current_phase = "optimization"
        # CDR3 beta = 0, FK particles active, full gradient guidance

    # Inject phase info into feats for DiffusionConditioning
    network_condition_kwargs["feats"]["current_phase"] = current_phase
    network_condition_kwargs["feats"]["diffusion_step_fraction"] = (
        1.0 - step_idx / num_sampling_steps
    )

    # Phase 3: Activate FK particles at transition
    if step_idx == phase2_end and steering_args.get("gl_combined"):
        # Replicate coords into particles (same logic as N1)
        ...

    # Phase 3: Apply FK resampling + gradient guidance
    if current_phase == "optimization":
        # ... existing FK + gradient code ...
```

**2. `model/modules/diffusion_conditioning.py`** — Phase-aware beta

The time-dependent beta from N4 handles this automatically: `beta_early=0.15, beta_late=0.0` means beta is active in early steps and fades out.

**3. `main.py`** — CLI flags

```python
# In BoltzSteeringParams:
gl_combined: bool = False
gl_phase1_fraction: float = 0.33
gl_phase2_fraction: float = 0.66
gl_phase3_particles: int = 5
gl_phase1_beta: float = 0.15
```

#### Estimated Complexity: Medium (similar to N1 but with three phases)

---

### N9: V-Fixed (Corrected Embedding CDR3)

**Goal**: Fix V's failed self-reference target by using biologically appropriate CDR3 targets.

**Architecture**: Modifies V's existing embedding optimization in `DiffusionConditioning` to use a different target computation.

#### Root Cause of V's Failure

V failed because it pushed CDR3 pair embeddings toward the mean of **non-CDR3 framework** pair embeddings. CDR3 loops are inherently different from framework — forcing them to resemble framework erases their specificity.

#### Fix: Alternative Target Definitions

**Option A: CDR3-Antigen Interface Target**
- Target = mean of CDR3-to-antigen pair embeddings (z[cdr3, antigen])
- Intuition: "make CDR3 pairs look like they're interacting with the antigen"
- This steers CDR3 toward binding rather than toward generic structure

**Option B: Self-Consistency Target**
- Target = CDR3 pair embeddings from the current trunk output z (no modification)
- Apply a small perturbation: z_target = z_cdr3 + epsilon * direction
- Where direction is toward antigen residues in embedding space
- Much milder than V's approach

**Option C: Canonical CDR3 Template Embeddings (requires data)**
- Pre-compute z_cdr3 from known canonical CDR3 structures
- Store as reference embeddings per canonical class
- Target = matching canonical class embedding

#### Files to Modify

**1. `model/modules/diffusion_conditioning.py`** — Replace target computation

V's existing code (in the y_hierarchical/v_embed worktrees) computes:
```python
z_target = z[~cdr3_mask].mean(dim=...)  # Framework mean — THIS IS WRONG
```

Replace with Option A:

```python
# NEW: CDR3-Antigen interface target
if "cdr3_antigen_token_mask" in feats:
    ag_mask = feats["cdr3_antigen_token_mask"].to(z.device).to(torch.bool)
    # Target: mean of CDR3-row, antigen-column pair embeddings
    cdr3_ag_pairs = z[:, cdr3_mask][:, :, ag_mask]  # [batch, n_cdr3, n_ag, z_dim]
    z_target = cdr3_ag_pairs.mean(dim=(1, 2))  # [batch, z_dim]
    # This target represents "what binding CDR3-antigen pairs look like"
```

**2. Key parameter changes**:
- Reduce strength from 1.0 to 0.05-0.1
- Reduce optimization steps from 10 to 1-3
- Apply only during first 30% of denoising steps
- Only modify CDR3-CDR3 pairs (not CDR3-framework)

```python
# Modified optimization:
if step_fraction > 0.7:  # Only first 30% of steps (step_fraction goes from 1.0 to 0.0)
    strength = steering_args.get("embed_cdr3_strength", 0.05)  # Reduced from 1.0
    num_opt_steps = steering_args.get("embed_cdr3_opt_steps", 2)  # Reduced from 10
    lr = strength * 0.01 * (1 - opt_step / num_opt_steps)  # Slower learning rate

    # Only modify CDR3-CDR3 pair embeddings (not CDR3-framework!)
    z_cdr3 = z[:, cdr3_mask][:, :, cdr3_mask]
    gradient = 2.0 * (z_cdr3 - z_target)  # MSE gradient
    z[:, cdr3_mask][:, :, cdr3_mask] -= lr * gradient
```

**3. `main.py`** — CLI flags

```python
# In BoltzSteeringParams:
embed_cdr3_steering: bool = False
embed_cdr3_target: str = "interface"    # "interface" | "self_consistency" | "canonical"
embed_cdr3_strength: float = 0.05
embed_cdr3_opt_steps: int = 2
embed_cdr3_early_only_fraction: float = 0.3
```

#### Estimated Complexity: Medium (modifying V's existing code with new targets)

---

### N10: E-Refined (Softer Blind Scanning)

**Goal**: Fix E's structural degradation by reducing beta magnitudes and increasing spatial resolution.

**Architecture**: Uses existing E's region-specific beta-scanning but with softer parameters and more regions.

#### Key Parameter Changes from E

| Parameter | E (failed) | N10 (proposed) |
|-----------|-----------|----------------|
| Beta emphasis | 0.5 | **0.15-0.2** |
| Beta de-emphasis | -0.3 | **-0.05 to -0.1** |
| Number of regions | 10 | **20-30** |
| Contact threshold | 8.0 A | **5.0-6.0 A** |
| Heatmap aggregation | Unweighted contact count | **Confidence-weighted** |

#### Files to Modify

**1. In E's existing implementation** — Change default parameters

The E worktree (`e_blind_scanning/`) modified `schema.py` and `featurizerv2.py` to support region-specific beta-scaling. N10 reuses this infrastructure with different parameter values.

**2. `data/parse/schema.py`** — Update defaults in region scanning constraint

```python
# In the epitope_scanning constraint parser:
default_emphasis = 0.2      # Was 0.5
default_deemphasis = -0.1   # Was -0.3
default_num_regions = 25    # Was 10
default_contact_threshold = 5.5  # Was 8.0
```

**3. Contact heatmap with confidence weighting**

In the epitope scanning aggregation logic (likely in a post-processing script):

```python
def aggregate_epitope_heatmap(predictions, feats_list):
    """Confidence-weighted contact heatmap aggregation.

    Instead of: heatmap[residue] += 1 for each contact
    Use: heatmap[residue] += confidence_score * (1.0 / distance)

    This weights contacts from high-confidence predictions more heavily
    and closer contacts more heavily.
    """
    heatmap = np.zeros(n_antigen_residues)
    for pred, feats, confidence in zip(predictions, feats_list, confidences):
        contacts = extract_contacts(pred, feats, threshold=5.5)
        for residue_id, distance in contacts:
            heatmap[residue_id] += confidence * (1.0 / max(distance, 1.0))
    return heatmap / heatmap.max()  # Normalize to [0, 1]
```

**4. Multi-step steering** — Apply beta-scaling during recycling/pairformer steps too

E only applies beta-scaling during diffusion conditioning. For stronger effect without increasing beta magnitude, apply the scaling earlier in the pipeline:

In `model/modules/trunkv2.py` (or wherever the Pairformer module is invoked):

```python
# After each Pairformer layer, optionally apply region-specific scaling
# This requires passing the region mask through the trunk
# Caution: This is a more invasive change
```

**Decision**: Start with parameter changes only (non-invasive). Only add multi-step steering if parameter tuning is insufficient.

#### Estimated Complexity: Low (parameter changes) to Medium (multi-step steering)

---

## Tier 3: High Impact if Successful, Higher Complexity

---

### N11: Full Pipeline Q→K→A (Discover→Focus→Refine)

**Goal**: Complete blind-to-refined pipeline combining three strategies.

**Architecture**: Three-stage external pipeline script.

#### Implementation

**New script: `scripts/pipeline_full_qka.py`**

```python
"""
Full Pipeline: Q (discover) → K (focus) → A (refine)

Stage 1: Q identifies candidate epitope regions (3 rounds)
Stage 2: K ensemble validates and narrows to best epitope patches (5-8 configs)
Stage 3: A with FK particles refines the complex targeting best patch (5-10 particles)
"""

def stage1_epitope_discovery(input_yaml_dir, output_dir):
    """Run Q's iterative epitope refinement.
    Output: Per-complex epitope heatmaps with candidate residue sets.
    """
    # boltz predict --use_potentials --iterative_epitope_refinement
    # Extract top 2-3 epitope patches per complex
    pass

def stage2_epitope_validation(input_yaml_dir, epitope_patches, output_dir):
    """Run K ensemble focused on Q's candidate patches.

    For each complex, for each candidate patch:
      Generate YAML with CDR3 beta-scaling emphasizing that patch
      Run Boltz2 with beta-scaling
    Select the patch with best confidence / contact quality.
    """
    # For each patch:
    #   boltz predict --use_potentials (with cdr3_beta_scaling + antigen pairs)
    pass

def stage3_refinement(input_yaml_dir, best_patches, output_dir):
    """Run A-style FK resampling targeting the validated epitope.

    Generate YAML with:
    - antigen_orientation constraint focused on best patch
    - CDR3 beta-scaling (mild, from N4)
    - FK particles = 5-10

    boltz predict --use_potentials --antigen_steering --num_particles 10
    """
    pass

def main():
    patches = stage1_epitope_discovery(...)
    best = stage2_epitope_validation(..., patches, ...)
    stage3_refinement(..., best, ...)
```

#### Key Integration Points

- Stage 1 output → Stage 2 input: Epitope residue IDs converted to region-specific YAML constraints
- Stage 2 output → Stage 3 input: Best epitope patch converted to `antigen_orientation` constraint YAML
- All stages use existing Boltz2 features (no model code changes)

#### Estimated Complexity: High (pipeline orchestration, error handling between stages)

---

### N12: Learned Re-Ranking with Steering Features

**Goal**: Train a lightweight model to predict DockQ from steering features.

**Architecture**: Gradient-boosted tree or small MLP trained on Round 1 data.

#### Implementation

**1. New script: `scripts/train_reranker.py`**

```python
"""
Train a lightweight re-ranker on Round 1 evaluation data.

Features (per model):
- confidence_score, iptm, ptm (Boltz2 confidence metrics)
- complex_plddt, complex_iplddt (structural confidence)
- cdr3_h_plddt, cdr3_l_plddt (CDR-specific confidence)
- n_contacts_5A, n_contacts_8A (contact counts)
- min_cdr_ag_distance (closest CDR-antigen distance)
- bsa_ratio (buried surface area ratio)
- h3_l3_contact_ratio (H3/L3 dominance)
- steering_method (one-hot: A, G, K, L, Q, Y, baseline)

Target: DockQ_AbAg (from ground truth)

Training: Leave-one-complex-out cross-validation on 47 complexes.
"""

import lightgbm as lgb
from sklearn.model_selection import LeaveOneGroupOut

def extract_features(predictions_dir, eval_results_csv):
    """Extract features from all Round 1 predictions.

    For each strategy × complex × model:
      1. Load predicted structure
      2. Load Boltz2 confidence JSON
      3. Compute contact features
      4. Load DockQ from eval_results_csv (target)
    """
    pass

def train_reranker(features_df, target_col="dockq_abag"):
    """Train LightGBM with LOOCV."""
    logo = LeaveOneGroupOut()
    model = lgb.LGBMRegressor(n_estimators=100, max_depth=4, learning_rate=0.1)
    # Cross-validation...
    return model

def save_reranker(model, output_path):
    """Save model for use in N6's SteeringReranker."""
    pass
```

**2. Data sources**: All evaluation results from Round 1 (47 complexes × 11 strategies × 5 models = ~2585 data points).

**3. Integration**: The trained model replaces the simple linear combination in N6's `SteeringReranker.rerank()` method.

#### Estimated Complexity: Medium-High (data collection, training, cross-validation)

---

### N13: Coarse Rigid-Body Pre-Step + Y

**Goal**: Provide reasonable starting orientations for the ~40% of complexes that fail at global orientation.

**Architecture**: Run a fast rigid-body docking scan before Boltz2 diffusion, then initialize diffusion from the best poses.

#### Implementation Approach

**1. New module: `model/modules/rigid_body_prescan.py`**

```python
"""
Coarse rigid-body orientation pre-scan.

Generates 3-5 candidate global orientations by:
1. Sampling random rotations of antigen relative to antibody
2. Scoring each rotation by CDR-antigen surface complementarity
3. Selecting top-K orientations as diffusion starting points
"""

import torch
from scipy.spatial.transform import Rotation

def generate_candidate_orientations(
    antibody_coords,   # [N_ab, 3] CA coordinates
    antigen_coords,    # [N_ag, 3] CA coordinates
    cdr_mask,          # [N_ab] boolean mask for CDR residues
    n_candidates=100,  # Number of random orientations to try
    n_keep=3,          # Number of best orientations to keep
):
    """Score orientations by CDR-antigen proximity.

    For each random rotation of the antigen:
    1. Rotate antigen around its center of mass
    2. Translate to place antigen COM at CDR COM + offset
    3. Score = number of CDR-antigen CA pairs within 15A
    4. Keep top-K scoring orientations
    """
    cdr_com = antibody_coords[cdr_mask].mean(dim=0)
    ag_com = antigen_coords.mean(dim=0)

    best_orientations = []
    for _ in range(n_candidates):
        R = torch.tensor(Rotation.random().as_matrix(), dtype=torch.float32)
        ag_rotated = (antigen_coords - ag_com) @ R.T + cdr_com + random_offset()
        score = ((torch.cdist(antibody_coords[cdr_mask], ag_rotated) < 15.0)
                 .sum().item())
        best_orientations.append((R, score))

    best_orientations.sort(key=lambda x: x[1], reverse=True)
    return [o[0] for o in best_orientations[:n_keep]]
```

**2. Integration in `model/modules/diffusionv2.py`** — Initialize noise from oriented coordinates

```python
# In sample() method, before the main loop (line 343-345):
if steering_args.get("rigid_body_prescan"):
    orientations = generate_candidate_orientations(
        antibody_coords=get_antibody_coords(network_condition_kwargs["feats"]),
        antigen_coords=get_antigen_coords(network_condition_kwargs["feats"]),
        cdr_mask=get_cdr_mask(network_condition_kwargs["feats"]),
        n_keep=steering_args.get("prescan_n_orientations", 3),
    )
    # Initialize different diffusion samples with different orientations
    # Instead of pure noise, use: oriented_coords + init_sigma * noise
    for i, R in enumerate(orientations):
        ag_atoms = get_antigen_atom_indices(feats)
        atom_coords[i, ag_atoms] = R @ atom_coords[i, ag_atoms].T  # Apply rotation
```

**3. Challenge**: This requires knowing which atoms belong to the antigen before diffusion starts. The information is available in `feats` (antigen chain indices), but the coordinate initialization happens before the network has processed any sequence information. The pre-scan would operate on a rough initial geometry, not on learned representations.

**Alternative**: Use the trunk output `z` to score orientations in embedding space (more aligned with the model's learned priors).

#### Estimated Complexity: High (requires careful coordinate initialization + testing)

---

### N14: W-Validated (Embedding Interface with Learned Head)

**Goal**: Test whether pair embeddings contain binding signal, then learn to extract it.

**Architecture**: Two-phase approach. Phase A validates the embedding signal offline. Phase B (if validated) trains a contact prediction head.

#### Implementation

**Phase A: Validate embedding signal** — New analysis script

```python
# scripts/validate_embedding_signal.py
"""
Test whether ||z[i,j]|| correlates with true contacts for CDR-antigen pairs.

For each of the 47 test complexes:
1. Run vanilla Boltz2, extract z_trunk from the model
2. For all CDR-antigen token pairs (i,j):
   - Compute ||z[i,j]|| (L2 norm of pair embedding)
   - Label: 1 if true contact (heavy atom < 5A in crystal), 0 otherwise
3. Compute AUC-ROC of ||z[i,j]|| as a binary classifier of contacts

If AUC-ROC < 0.6: Embedding norms do NOT encode binding → STOP
If AUC-ROC > 0.7: Signal exists → proceed to Phase B
"""
```

This requires extracting `z_trunk` from the model. The hook is already available: `boltz2.py` passes `z_trunk=z.float()` to the structure module (see W's implementation notes in the implementation guide).

**Phase B: Train contact prediction MLP** — New module

```python
# model/modules/contact_head.py
"""Small MLP trained on pair embeddings to predict contact probability."""

class ContactPredictionHead(nn.Module):
    """Maps pair embedding z[i,j] to contact probability p(contact|z[i,j])."""

    def __init__(self, z_dim=128, hidden_dim=64):
        super().__init__()
        self.mlp = nn.Sequential(
            nn.Linear(z_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, 1),
            nn.Sigmoid(),
        )

    def forward(self, z_pairs):
        """z_pairs: [N_pairs, z_dim] → [N_pairs] contact probabilities."""
        return self.mlp(z_pairs).squeeze(-1)
```

Training data: Extract (z[i,j], contact_label) pairs from all 47 complexes. Train with cross-entropy loss, leave-one-out CV.

**Phase C: Integration with W-style steering**

Replace raw embedding norms in `EmbeddingInterfacePotential` with learned contact probabilities:

```python
# In model/potentials/potentials.py (W's EmbeddingInterfacePotential):
def set_embedding_weights(self, z_trunk, feats):
    # OLD: weights = z_trunk[cdr_idx, ag_idx].norm(dim=-1)
    # NEW: weights = contact_head(z_trunk[cdr_idx, ag_idx])
    contact_head = load_pretrained_contact_head()
    z_pairs = z_trunk[0, cdr_idx][:, ag_idx]  # [N_cdr, N_ag, z_dim]
    weights = contact_head(z_pairs.reshape(-1, z_dim)).reshape(N_cdr, N_ag)
    self.embedding_weights = weights / weights.max()  # Normalize to [0, 1]
```

#### Estimated Complexity: High (validation analysis + MLP training + integration)

---

### N15: Ensemble Union Strategy

**Goal**: Run multiple strategies in parallel and select the best prediction per complex.

**Architecture**: External orchestration script that runs A, G, Y, and baseline in parallel, pools all models, and selects the best using N6's re-ranker.

#### Implementation

**New script: `scripts/ensemble_union.py`**

```python
"""
Ensemble Union: Run multiple strategies, pool all models, re-rank.

For each complex:
1. Run vanilla Boltz2 (5 models)
2. Run A with FK 5 particles (5 models)
3. Run G+ progressive (5 models)
4. Run Y hierarchical (5 models)
Total: 20 models per complex

Re-rank all 20 using N6's SteeringReranker (or N12's learned model).
Output: Top-5 models from the pool.
"""

import subprocess
from pathlib import Path

def run_strategy(yaml_file, strategy, output_dir):
    """Run a single strategy on a single complex."""
    cmd = ["boltz", "predict", str(yaml_file), "--output", str(output_dir)]

    if strategy == "baseline":
        pass  # No extra flags
    elif strategy == "A":
        cmd += ["--use_potentials", "--antigen_steering", "--num_particles", "5"]
    elif strategy == "G":
        cmd += ["--use_potentials", "--progressive_cdr_refinement"]
    elif strategy == "Y":
        cmd += ["--use_potentials", "--hierarchical_steering"]

    subprocess.run(cmd, check=True)

def pool_and_rerank(predictions_dirs, feats, reranker):
    """Pool all models from all strategies, re-rank by composite score.

    1. Load all predicted structures + confidence JSONs
    2. Compute re-ranking features for each model
    3. Apply reranker.rerank() to select top-5
    4. Output top-5 with provenance labels (which strategy produced each)
    """
    all_models = []
    for strategy_dir in predictions_dirs:
        for model_file in strategy_dir.glob("*_model_*.cif"):
            model = load_model(model_file)
            model.strategy = strategy_dir.name
            all_models.append(model)

    rankings = reranker.rerank(all_models, feats)
    return [all_models[i] for i in rankings[:5]]

def main():
    strategies = ["baseline", "A", "G", "Y"]

    for yaml_file in input_dir.glob("*.yml"):
        for strategy in strategies:
            run_strategy(yaml_file, strategy,
                output_dir / f"{yaml_file.stem}_{strategy}")

        # Pool and re-rank
        pred_dirs = [output_dir / f"{yaml_file.stem}_{s}" for s in strategies]
        best_models = pool_and_rerank(pred_dirs, feats, reranker)
```

#### Parallelization

The four strategy runs per complex are independent and can run on separate GPUs or in sequence. For 47 complexes × 4 strategies = 188 Boltz2 runs, this is computationally expensive but embarrassingly parallel.

#### Estimated Complexity: Medium-High (orchestration script + N6 re-ranker dependency)

---

## Cross-Cutting Concerns

### Passing Diffusion Step Information to DiffusionConditioning

Multiple strategies (N1, N4, N7, N8, N9) require the diffusion step fraction to be available in `DiffusionConditioning.forward()`. The cleanest approach:

**In `model/modules/diffusionv2.py`**, inject step info before the network forward call:

```python
# Before line 388 (network forward call in the main loop):
# Inject diffusion step fraction into feats for conditioning modules
network_condition_kwargs["feats"]["diffusion_step_fraction"] = torch.tensor(
    1.0 - step_idx / num_sampling_steps,  # 1.0 at start, 0.0 at end
    device=self.device
)
network_condition_kwargs["feats"]["diffusion_step_idx"] = torch.tensor(
    step_idx, device=self.device
)
```

This is a **one-time change** that enables all time-dependent conditioning strategies.

### FK Particle Activation Mid-Loop

Strategies N1, N8 need to activate FK particles partway through the loop (not from the beginning). This requires:

1. Expanding `atom_coords` dimensionality mid-loop
2. Initializing FK tracking tensors (`energy_traj`, `resample_weights`)
3. Adjusting `atom_mask` and `network_condition_kwargs` for the new multiplicity

The cleanest approach is a helper function:

```python
def activate_fk_particles(atom_coords, atom_mask, n_particles, device):
    """Replicate current state into n_particles copies for FK resampling."""
    new_coords = atom_coords.repeat_interleave(n_particles, dim=0)
    new_mask = atom_mask.repeat_interleave(n_particles, dim=0)
    new_multiplicity = atom_coords.shape[0] * n_particles
    energy_traj = torch.empty((new_multiplicity, 0), device=device)
    resample_weights = torch.ones(new_multiplicity, device=device).reshape(
        -1, n_particles)
    return new_coords, new_mask, energy_traj, resample_weights, new_multiplicity
```

### Worktree Strategy for N1-N15

Each new strategy should be developed in its own git worktree:

```
worktree_steering_cdr3/
├── n1_fk_hierarchical/     # Branch from common base
├── n2_q_to_b2_pipeline/    # Pipeline scripts only
├── n3_conditional_steering/
├── n4_l_light/
├── n5_k_ensemble/           # Pipeline scripts only
├── n6_custom_reranking/
├── n7_y_plus/               # Branch from y_hierarchical_steering
├── n8_gl_combined/
├── n9_v_fixed/              # Branch from v_embed_cdr3_steering
├── n10_e_refined/           # Branch from e_blind_scanning
├── n11_full_pipeline/       # Pipeline scripts only
├── n12_learned_reranking/
├── n13_rigid_body_prescan/
├── n14_w_validated/         # Branch from w_embed_interface
└── n15_ensemble_union/      # Pipeline scripts only
```

Strategies N2, N5, N11, N15 are **pipeline-only** (external scripts, no model changes). They can share the same boltz codebase and don't need separate worktrees for model code.

---

## Implementation Order Recommendation

```
Week 1: N3 (conditional steering) + N5 (K ensemble) + N4 (L-light)
         → Low complexity, addresses regressions and leverages existing features
         → N4's step-fraction injection enables all time-dependent strategies

Week 2: N1 (A+Y hybrid) + N6 (custom re-ranking)
         → Core architectural changes; N6 benefits all strategies

Week 3: N2 (Q→B2 pipeline) + N7 (Y+)
         → Pipeline integration; stronger hierarchical steering

Week 4: N8 (G+L combined) + N10 (E-refined)
         → Three-phase loop; softer blind scanning

Week 5: N12 (learned re-ranking) + N15 (ensemble union)
         → Data-driven selection; aggregate all wins

Week 6+: N9 (V-fixed) + N11 (full pipeline) + N13 (rigid-body) + N14 (W-validated)
         → Research directions; higher risk/reward
```

---

## Summary: Files Modified Per Strategy

| Strategy | main.py | diffusionv2.py | diffusion_conditioning.py | potentials.py | schema.py | featurizerv2.py | New Scripts |
|----------|---------|---------------|--------------------------|---------------|-----------|----------------|-------------|
| N1 | Yes | **Major** | Yes | - | - | - | - |
| N2 | - | - | - | - | - | - | **pipeline_q_to_b2.py** |
| N3 | Yes | Yes | - | - | - | - | conditional_steering.py (opt) |
| N4 | Yes | Minor | **Major** | - | Yes | Yes | - |
| N5 | - | - | - | - | - | - | **k_ensemble.py** |
| N6 | Yes | - | - | - | - | - | **reranking.py** (new module) |
| N7 | Yes | - | - | Yes | Yes | - | - |
| N8 | Yes | **Major** | Minor | - | - | - | - |
| N9 | Yes | - | **Major** | - | - | - | - |
| N10 | - | - | - | - | Minor | - | - |
| N11 | - | - | - | - | - | - | **pipeline_full_qka.py** |
| N12 | - | - | - | - | - | - | **train_reranker.py** |
| N13 | Yes | **Major** | - | - | - | - | **rigid_body_prescan.py** |
| N14 | - | - | - | Yes | - | - | **validate_embedding.py**, **contact_head.py** |
| N15 | - | - | - | - | - | - | **ensemble_union.py** |

**Bold** = primary change location. Minor = small additions (flags, injection of step info).
