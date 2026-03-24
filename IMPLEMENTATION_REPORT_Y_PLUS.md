# Implementation Report: Strategy Y+ — Hierarchical Steering v2

**Date**: 2026-03-23
**Branch**: `y_v2_hierarch`
**Base**: Round 1 evaluation codebase (with strategies A, D, E, G, K, L, O, Q, V, W, Y)

---

## 1. Executive Summary

Strategy Y was the most promising strategy from Round 1 — the **only one with a statistically significant structural improvement** (CDR-H3 RMSD −0.30 Å, p=0.025). Y+ amplifies both phases of the hierarchical approach:

- **Y+.1**: Stronger time-varying beta-scaling in the early (embedding) phase
- **Y+.2**: Enhanced late-phase coordinate potentials with CDR3-antigen proximity
- **Y+.3**: Incorporate predicted epitope information from Strategy Q's output PDB structures
- **Y+.4**: Guidance weight scaling for hard cases (configurable per-complex)

Y+ reuses infrastructure patterns established by Strategy G+ (progressive steering), specifically the time-varying beta mechanism in the diffusion loop.

---

## 2. Architecture Overview

### Data Flow

```
YAML Input (hierarchical_steering constraint)
    ↓
schema.py: parse_boltz_schema()
    ├─ Parse antigen_chain, cdr_regions, beta_max, beta_schedule, guidance_weight_scale
    ├─ If predicted_epitope_pdb: parse PDB → extract epitope residue indices
    └─ Store as hierarchical_steering_constraints tuple in InferenceOptions
    ↓
inferencev2.py: PredictionDataset.__getitem__()
    └─ Extract from record.inference_options, pass to featurizer
    ↓
featurizerv2.py: process_hierarchical_steering_constraints()
    ├─ Create hierarchical_cdr_token_mask [N_tokens] — CDR residue positions
    ├─ Create hierarchical_antigen_atom_index [N_antigen] — antigen CA indices
    ├─ Create hierarchical_cdr_atom_index [N_cdr] — CDR CA indices
    ├─ Create hierarchical_epitope_token_mask [N_tokens] — epitope residue positions
    ├─ Encode beta_max, beta_schedule, contact_threshold, guidance_weight_scale
    └─ Return feature dict
    ↓
diffusion_conditioning.py: DiffusionConditioning.forward()
    └─ Skip static CDR3 beta scaling when hierarchical active (beta handled per-step)
    ↓
diffusionv2.py: AtomDiffusion.sample() — per diffusion step
    ├─ Compute beta_t from schedule (linear/cosine/step) using steering_t
    ├─ Build CDR pair mask, scale token_trans_bias *= (1 + beta_t * mask)
    ├─ Run forward pass
    └─ Restore original token_trans_bias
    ↓
potentials.py: HierarchicalAntigenOrientationPotential + CDR3AntigenProximityPotential
    ├─ HierarchicalAntigenOrientation: reads from hierarchical_* features,
    │   supports epitope focusing, late-phase schedule, guidance_weight_scale
    └─ CDR3AntigenProximity: late-phase-only, penalizes CDR-antigen distance
```

### Two-Phase Design

| Phase | Diffusion Time | Mechanism | Purpose |
|-------|---------------|-----------|---------|
| **Early** (t=1.0→0.5) | Steps 0–100 | Time-varying CDR beta-scaling in embedding space | Global CDR3 orientation, coarse-grained steering |
| **Late** (t=0.5→0.0) | Steps 100–200 | Coordinate-space potentials (antigen orientation + CDR proximity) | Fine-grained interface refinement |

---

## 3. YAML Format

```yaml
constraints:
  - hierarchical_steering:
      antigen_chain: A                    # Chain ID of the antigen
      contact_threshold: 8.0             # Distance threshold (Å) for contact potential
      cdr_regions:                       # CDR regions for beta-scaling and proximity
        - chain: B                       # Heavy chain
          start_res: 97                  # 1-indexed start
          end_res: 115                   # 1-indexed end
        - chain: C                       # Light chain
          start_res: 88
          end_res: 100
      beta_max: 0.8                      # Max beta scaling factor (default 0.8)
      beta_schedule: linear              # "linear", "cosine", or "step" (default "linear")
      predicted_epitope_pdb: path/to.pdb # Optional: PDB from Q's output for epitope extraction
      epitope_contact_threshold: 10.0    # Optional: Å threshold for epitope extraction from PDB
      guidance_weight_scale: 1.0         # Optional: Scale guidance weights (Y+.4, default 1.0)
      force: true
```

---

## 4. File-by-File Changes

### 4.1 `src/boltz/main.py`

**Changes:**
1. Add `hierarchical_steering: bool = False` to `BoltzSteeringParams` dataclass
2. Add `--hierarchical_steering` CLI flag
3. Add `hierarchical_steering` parameter to `predict()` function signature
4. Wire `steering_args.hierarchical_steering = hierarchical_steering`
5. Add validation: requires `--use_potentials`

### 4.2 `src/boltz/data/types.py`

**Changes:**
Add to `InferenceOptions` dataclass:
```python
# Hierarchical steering constraints:
# (antigen_chain_id, contact_threshold, cdr_regions, beta_max, beta_schedule,
#  epitope_residues, guidance_weight_scale, force)
hierarchical_steering_constraints: Optional[
    list[tuple[int, float, list[tuple[int, int, int]], float, str,
               Optional[list[int]], float, bool]]
] = None
```

### 4.3 `src/boltz/data/parse/schema.py`

**Changes:**
1. Add `extract_epitope_from_pdb()` utility function — parses PDB file, extracts antigen CA atoms near CDR CA atoms within threshold
2. Add `hierarchical_steering` constraint parsing block in the constraint loop
3. Store parsed constraints in `hierarchical_steering_constraints` list
4. Pass to `InferenceOptions`

**PDB Epitope Extraction Logic:**
- Read PDB ATOM records, filter CA atoms
- Group by chain letter
- Compute min distance from each antigen CA to any CDR CA
- Return antigen residue indices where min_dist < epitope_contact_threshold
- Convert to 0-indexed for internal storage

### 4.4 `src/boltz/data/feature/featurizerv2.py`

**Changes:**
1. Add `process_hierarchical_steering_constraints()` function
2. Add `inference_hierarchical_steering_constraints` parameter to `process()`
3. Call it and merge features into return dict

**Feature Tensors Produced:**
| Feature | Shape | Type | Description |
|---------|-------|------|-------------|
| `hierarchical_cdr_token_mask` | [N_tokens] | bool | CDR residue positions |
| `hierarchical_beta_max` | [1] | float32 | Maximum beta value |
| `hierarchical_beta_schedule` | [1] | long | 0=linear, 1=cosine, 2=step |
| `hierarchical_antigen_atom_index` | [N_antigen] | long | Antigen CA atom indices |
| `hierarchical_cdr_atom_index` | [N_cdr] | long | CDR CA atom indices |
| `hierarchical_contact_threshold` | [1] | float32 | Contact distance threshold |
| `hierarchical_epitope_token_mask` | [N_tokens] | bool | Epitope residue positions (from PDB) |
| `hierarchical_guidance_weight_scale` | [1] | float32 | Guidance weight multiplier |

### 4.5 `src/boltz/model/modules/diffusion_conditioning.py`

**Changes:**
Skip static CDR3 beta scaling when hierarchical steering is active:
```python
hierarchical_active = (
    "hierarchical_cdr_token_mask" in feats
    and feats["hierarchical_cdr_token_mask"].any()
)
if not hierarchical_active and "cdr3_token_mask" in feats and ...:
    # existing static beta scaling
```

### 4.6 `src/boltz/model/modules/diffusionv2.py`

**Changes:**
Inside the denoising loop (after `steering_t` computation), add per-step hierarchical beta scaling:

1. Check if `hierarchical_steering` is active in `steering_args`
2. Read CDR mask, beta_max, schedule_type from features
3. Compute `beta_t` from schedule:
   - Linear: `beta_t = beta_max * steering_t`
   - Cosine: `beta_t = beta_max * 0.5 * (1 + cos(π * (1 - steering_t)))`
   - Step: full `beta_max` for t>0.7, 30% for 0.3<t≤0.7, 0 for t≤0.3
4. Build CDR pair mask `[n_tokens, n_tokens]`
5. Scale `token_trans_bias *= (1 + beta_t * mask_4d)`
6. After forward pass, restore original `token_trans_bias`

### 4.7 `src/boltz/model/potentials/potentials.py`

**New Classes:**

**`HierarchicalAntigenOrientationPotential`** (extends `AntigenOrientationPotential`):
- Reads from `hierarchical_antigen_atom_index`, `hierarchical_cdr_atom_index`, `hierarchical_contact_threshold`
- If `hierarchical_epitope_token_mask` is available, restricts antigen atoms to epitope residues only
- Overrides `compute()` to apply `hierarchical_guidance_weight_scale`

**`CDR3AntigenProximityPotential`** (extends `AntigenOrientationPotential`):
- Same atom indices but different purpose: late-phase-only direct proximity enforcement
- Uses a tighter threshold than the orientation potential

**`get_potentials()` update:**
```python
if boltz2 and steering_args.get("hierarchical_steering", False):
    potentials.append(
        HierarchicalAntigenOrientationPotential(
            parameters={
                "guidance_interval": 2,
                "guidance_weight": PiecewiseStepFunction(
                    thresholds=[0.5],
                    values=[0.0, 1.5]  # Zero early, strong late
                ),
                "resampling_weight": PiecewiseStepFunction(
                    thresholds=[0.5],
                    values=[0.0, 1.0]  # Zero early, full late
                ),
                "union_lambda": ExponentialInterpolation(
                    start=8.0, end=0.0, alpha=-2.0
                ),
            }
        )
    )
    potentials.append(
        CDR3AntigenProximityPotential(
            parameters={
                "guidance_interval": 2,
                "guidance_weight": PiecewiseStepFunction(
                    thresholds=[0.4],
                    values=[0.0, 0.8]  # Only active in late phase
                ),
                "resampling_weight": PiecewiseStepFunction(
                    thresholds=[0.4],
                    values=[0.0, 0.5]
                ),
                "union_lambda": ExponentialInterpolation(
                    start=8.0, end=0.0, alpha=-2.0
                ),
            }
        )
    )
```

### 4.8 `src/boltz/data/module/inferencev2.py`

**Changes:**
1. Extract `hierarchical_steering_constraints` from `options`
2. Pass to featurizer via `inference_hierarchical_steering_constraints` kwarg
3. Merge returned features

---

## 5. Improvement Mapping

| Y+ Sub-strategy | Implementation Location | Key Mechanism |
|-----------------|------------------------|---------------|
| **Y+.1** Stronger early beta | `diffusionv2.py` (per-step beta), `diffusion_conditioning.py` (skip static) | `beta_t = beta_max * f(steering_t)`, three schedules |
| **Y+.2** Enhanced late potentials | `potentials.py` (two new potential classes) | Late-phase-only CDR-antigen distance penalty |
| **Y+.3** Predicted epitope input | `schema.py` (PDB parsing), `featurizerv2.py` (epitope mask) | Focus late-phase potentials on predicted epitope residues |
| **Y+.4** Hard case scaling | `potentials.py` (guidance_weight_scale), YAML config | Per-complex guidance weight multiplier |

---

## 6. Test Plan

1. Create YAML for 7TRH_HBG with hierarchical_steering constraint
2. Include CDR3 regions from `examples/cdrs.csv` (H: 97-115, L: 88-100)
3. Point `predicted_epitope_pdb` to `predictions_examples/antigen_cut/boltz_results_7TRH_HBG/predictions/7TRH_HBG/7TRH_HBG_model_0.pdb`
4. Run: `boltz predict test.yaml --use_potentials --hierarchical_steering --output_format pdb --diffusion_samples 1 --override`
5. Verify: prediction completes without errors, output PDB is generated

---

## 7. Dependencies

- **From G+**: Time-varying beta pattern (implemented fresh here following G+'s pattern)
- **From Q**: Predicted PDB structures for epitope extraction (provided as file path in YAML)
- **Existing infrastructure**: `AntigenOrientationPotential`, `FlatBottomPotential`, `DistancePotential`, FK resampling loop
