# N1: A+Y Hybrid (FK Particles + Hierarchical Timing) — Implementation Report

**Date**: 2026-03-20
**Strategy**: N1 from ROUND_1_FINAL_EVALUATION.md Tier 1
**Branch**: `a_v2_fk_particles`

---

## 1. Strategy Overview

**Goal**: Combine Strategy A's FK particle resampling (best DockQ improvement, +0.046) with Strategy Y's hierarchical timing (only statistically significant CDR-H3 RMSD improvement, -0.30Å, p=0.025).

**Key Insight**: A provides the best docking improvement through particle-based importance sampling, while Y provides the only significant structural improvement through two-phase timing (embedding-space early, coordinate-space late). Combining them should address both global orientation (Y's embedding phase) and local refinement (A's particle selection).

**Architecture**: Two-phase diffusion:
- **Early phase** (steps 0 to transition_step): CDR3 beta-scaling on pair representations for global orientation guidance. No FK resampling, no gradient guidance — just standard denoising with modified attention biases.
- **Late phase** (transition_step to end): FK particles activated with coordinate-space potentials (AntigenOrientationPotential) for fine-grained refinement. Particles replicated at transition point, resampled by energy.

---

## 2. Files to Modify

| # | File | Change Description |
|---|------|--------------------|
| 1 | `src/boltz/main.py` | Add `--hybrid_fk_hierarchical`, `--hybrid_transition_fraction`, `--hybrid_early_beta`, `--hybrid_late_particles` CLI flags + `BoltzSteeringParams` fields |
| 2 | `src/boltz/model/modules/diffusionv2.py` | Two-phase logic in `AtomDiffusion.sample()`: early=denoising only, late=FK+guidance. Replicate particles at transition point. Pass `diffusion_step_fraction` to feats. |
| 3 | `src/boltz/model/modules/diffusion_conditioning.py` | Phase-aware beta-scaling: apply CDR3 beta only during early phase using `diffusion_step_fraction` from feats |
| 4 | `src/boltz/data/parse/schema.py` | Add `hybrid_fk_hierarchical` YAML constraint block (combines antigen_orientation + cdr3_beta_scaling in one constraint) |
| 5 | `src/boltz/data/feature/featurizerv2.py` | Add `process_hybrid_fk_hierarchical_constraints()` that generates both beta-scaling AND antigen orientation features |
| 6 | `src/boltz/data/types.py` | Add `hybrid_fk_hierarchical_constraints` field to `InferenceOptions` |
| 7 | `src/boltz/data/module/inferencev2.py` | Wire hybrid constraint from record to featurizer |
| 8 | `src/boltz/model/potentials/potentials.py` | Update `get_potentials()` to handle hybrid mode (antigen potential active only in late phase) |

---

## 3. Design Decisions

### 3.1 YAML Constraint Format

The hybrid strategy needs BOTH antigen orientation info AND CDR3 beta-scaling info from the YAML. Rather than requiring users to specify two separate constraint blocks, we create a single `hybrid_fk_hierarchical` constraint that encompasses both:

```yaml
constraints:
  - hybrid_fk_hierarchical:
      antigen_chain: A
      contact_threshold: 8.0
      transition_fraction: 0.5     # When to switch phases (0.5 = halfway)
      early_beta: 0.3              # Beta-scaling strength in early phase
      cdr_regions:
        - chain: B                 # Heavy chain CDR3
          start_res: 97
          end_res: 115
        - chain: C                 # Light chain CDR3
          start_res: 88
          end_res: 100
      force: true
```

This reuses the same CDR region specification format as `antigen_orientation` and `cdr3_beta_scaling`. The YAML provides all CDR indices (since they differ per sequence).

### 3.2 CLI Flags vs YAML Config

- `--hybrid_fk_hierarchical` (CLI flag): Enables the hybrid mode
- `--hybrid_late_particles` (CLI, default 5): Number of FK particles in late phase
- `transition_fraction` and `early_beta`: Specified in YAML per-complex (since CDR regions differ)
- `--num_particles`: Ignored when hybrid mode is active (uses `--hybrid_late_particles` instead)

### 3.3 Phase Transition Mechanism

At the transition point:
1. Current single trajectory `atom_coords` is replicated into N particles
2. FK tracking tensors (`energy_traj`, `resample_weights`) are initialized
3. `scaled_guidance_update` is initialized for gradient guidance
4. From this point forward, the existing FK + guidance code runs unchanged

This minimizes code changes — the late phase reuses 100% of the existing FK machinery.

### 3.4 Beta-Scaling Phase Awareness

The `DiffusionConditioning.forward()` currently applies beta-scaling uniformly. For N1, we make it phase-aware:

- Pass `diffusion_step_fraction` (1.0 → 0.0 as diffusion progresses) via `feats` dict
- In conditioning, compute a smooth ramp: full beta early, zero beta late
- The ramp uses `transition_fraction` from feats to know where to cut off
- Smooth transition over 10% of steps to avoid discontinuity

### 3.5 Potentials in Hybrid Mode

In hybrid mode, `get_potentials()` creates the same potentials as antigen_steering mode, but the diffusion loop only invokes them during the late phase (after particles are activated). This means:
- Early phase: No energy computation, no gradient guidance, no resampling
- Late phase: Full FK resampling + gradient guidance with AntigenOrientationPotential

---

## 4. Data Flow

```
User YAML: hybrid_fk_hierarchical constraint with CDR regions + antigen chain
  ↓
schema.py: Parse → (antigen_chain_id, contact_threshold, cdr_regions, transition_fraction, early_beta, force)
  ↓
types.py: InferenceOptions.hybrid_fk_hierarchical_constraints
  ↓
inferencev2.py: Wire to featurizer
  ↓
featurizerv2.py: process_hybrid_fk_hierarchical_constraints() generates:
  - cdr3_token_mask (for beta-scaling)
  - cdr3_beta_value (for beta-scaling)
  - antigen_atom_index (for AntigenOrientationPotential)
  - cdr_atom_index (for AntigenOrientationPotential)
  - antigen_orientation_threshold (for AntigenOrientationPotential)
  - hybrid_transition_fraction (for phase switching)
  ↓
batch dict → GPU → boltz2.py predict_step → diffusionv2.py sample()
  ↓
CLI: --hybrid_fk_hierarchical --hybrid_late_particles 5 --use_potentials
  ↓
main.py: steering_args.hybrid_fk_hierarchical = True, steering_args.hybrid_late_particles = 5
  ↓
diffusionv2.py sample():
  Step 0 to transition_step:
    - DiffusionConditioning applies beta=0.3 to CDR3 pair representations
    - Standard denoising (no FK, no gradient guidance)
  Step transition_step:
    - Replicate atom_coords into 5 particles
    - Initialize FK tensors
  Step transition_step to end:
    - FK resampling + gradient guidance with AntigenOrientationPotential
    - Final resampling selects best particle
  ↓
Output: Single best structure from FK selection
```

---

## 5. Implementation Details

### 5.1 `main.py` Changes

Add to `BoltzSteeringParams`:
```python
hybrid_fk_hierarchical: bool = False
hybrid_late_particles: int = 5
```

Add CLI options:
```python
@click.option("--hybrid_fk_hierarchical", is_flag=True, ...)
@click.option("--hybrid_late_particles", type=int, default=5, ...)
```

Wire in predict():
```python
steering_args.hybrid_fk_hierarchical = hybrid_fk_hierarchical
steering_args.hybrid_late_particles = hybrid_late_particles
```

### 5.2 `diffusionv2.py` Changes

In `sample()`, the key structural change:

1. **Before loop**: If hybrid mode, do NOT multiply multiplicity by particles yet. Only set up potentials.
2. **In loop**: Check `is_early_phase` flag per step. If early → skip FK/guidance blocks. If late → run them.
3. **At transition step**: Replicate coordinates into particles, initialize FK tensors.
4. **After loop**: Return best particle (existing code handles this via final resampling step).

### 5.3 `diffusion_conditioning.py` Changes

Modify beta-scaling block to be time-dependent:
- Read `diffusion_step_fraction` and `hybrid_transition_fraction` from feats
- Compute smooth ramp weight: 1.0 in early phase → 0.0 in late phase
- Multiply beta by ramp weight

### 5.4 `schema.py` Changes

Add parsing for `hybrid_fk_hierarchical` constraint block. Extracts:
- `antigen_chain`, `contact_threshold` (same as antigen_orientation)
- `cdr_regions` with chain/start_res/end_res (same format)
- `transition_fraction` (new, default 0.5)
- `early_beta` (new, default 0.3)
- `force` (same as antigen_orientation)

### 5.5 `featurizerv2.py` Changes

New `process_hybrid_fk_hierarchical_constraints()` function that:
1. Calls existing antigen atom index logic (from process_antigen_orientation_constraints)
2. Calls existing CDR3 mask logic (from process_cdr3_beta_constraints)
3. Adds `hybrid_transition_fraction` tensor to features
4. Returns merged feature dict

### 5.6 `types.py` Changes

Add to InferenceOptions:
```python
hybrid_fk_hierarchical_constraints: Optional[
    list[tuple[int, float, list[tuple[int, int, int]], float, float, bool]]
] = None
```

Tuple: (antigen_chain_id, contact_threshold, cdr_regions, transition_fraction, early_beta, force)

### 5.7 `inferencev2.py` Changes

Extract `hybrid_fk_hierarchical_constraints` from `record.inference_options` and pass to featurizer.

### 5.8 `potentials.py` Changes

In `get_potentials()`, add hybrid-specific antigen potential when `hybrid_fk_hierarchical` is True:
- Same AntigenOrientationPotential as antigen_steering
- But with schedule tuned for late-phase only (stronger since it has fewer steps)

---

## 6. Test Plan

### 6.1 Test Complex: 7TRH_HBG

From cdrs.csv:
- Heavy chain (B): CDR-H3 residues 97-115 (1-indexed)
- Light chain (C): CDR-L3 residues 88-100 (1-indexed)
- Antigen chain: A

YAML config:
```yaml
sequences:
- protein:
    id: A
    sequence: APLHLGKCNIAGWILGNPECESLSTASSWSYIVETPSSDNGTCYPGDFIDYEELREQLSSVSSFERFEIFPKTSSWPNHDSDKGVTAACPHAGAKSFYKNLIWLVKKGNSYPKLSKSYINDKGKEVLVLWGIHHPSTSADQQSLYQNADAYVFVGSSRYSKTFKPEIAIRPKVRDREGRMNYYWTLVEPGDKITFEATGNLVVPRYAFAMERNA
    msa: ../msa_antigen_cut_new/antigen_7TRH_HBG.a3m
- protein:
    id: B
    sequence: EVQLVESGGGLIQPGGSLRLSCEASAFTFSSYEMNWVRQAPGKGLEWVSYITSSGSRIYYADSVKGRFTISRDNAKNSLYLQMNSLRVEDTAVYYCARLLDSIVWGEGWYYGMDVWGQGTTVTVSG
    msa: empty
- protein:
    id: C
    sequence: SYELTQSPSVSVAPGRTARITCGGNDIGLKGVHWYQQKPGQAPVLVLYDNNHRPSGIPERFSGSISGDTATLTVTRVEADDGADYFCQVWDTSSGPPHVIFGGGTKLTVL
    msa: empty

constraints:
  - hybrid_fk_hierarchical:
      antigen_chain: A
      contact_threshold: 8.0
      transition_fraction: 0.5
      early_beta: 0.3
      cdr_regions:
        - chain: B
          start_res: 97
          end_res: 115
        - chain: C
          start_res: 88
          end_res: 100
      force: true
```

### 6.2 Test Command

```bash
conda activate boltz
boltz predict examples/7TRH/7TRH_HBG_hybrid.yml \
    --use_potentials \
    --hybrid_fk_hierarchical \
    --hybrid_late_particles 5 \
    --diffusion_samples 1 \
    --sampling_steps 50 \
    --override \
    --out_dir ./test_hybrid_output
```

Using `--sampling_steps 50` for faster testing (vs default 200).

### 6.3 Success Criteria

1. **Completion**: Prediction completes without errors or OOM
2. **Output**: Valid mmCIF file produced in output directory
3. **Phase behavior**: Early phase runs without FK overhead, late phase uses FK particles
4. **Memory**: 5 particles should stay well within GPU memory (vs A's 20 that caused OOM)

---

## 7. Risk Assessment

| Risk | Likelihood | Impact | Mitigation |
|------|-----------|--------|------------|
| OOM with 5 particles | Low | Medium | Default to 5 particles (vs A's 20); can reduce to 3 |
| Phase transition discontinuity | Medium | Low | Smooth beta ramp over 10% of steps |
| Late phase too short for FK convergence | Medium | Medium | Default transition at 50% leaves 100 steps for FK |
| Beta-scaling in early phase too weak | Low | Low | Tunable via YAML `early_beta` parameter |

---

## 8. Comparison to Existing Strategies

| Aspect | Strategy A | Strategy Y | N1 Hybrid |
|--------|-----------|-----------|-----------|
| FK particles | 20 (all steps) | 0 | 5 (late phase only) |
| Beta-scaling | None | CDR3 beta (all steps) | CDR3 beta (early phase only) |
| Gradient guidance | All steps | All steps | Late phase only |
| Memory usage | High (OOM 25%) | Low | Medium |
| Expected DockQ | +0.046 | +0.034 | +0.05-0.08 |
| Expected CDR-H3 | No change | -0.30Å | -0.3-0.5Å |
