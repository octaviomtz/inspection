# Idea Y: Hierarchical Steering - Detailed Implementation Guide

**Date**: 2026-03-11
**Status**: Implementation Planning
**Complexity**: HIGH (2-4 weeks)
**Strategy**: Embedding-first (robust), then Coordinates (precise)

---

## 1. Goal

Combine the **robustness** of embedding-space steering (early diffusion) with the **precision** of coordinate-space steering (late diffusion) in a smooth multi-stage approach, yielding higher quality antibody-antigen structure predictions across diverse inputs, including novel/out-of-distribution systems.

---

## 2. Core Insight

Structure prediction via diffusion proceeds from coarse to fine. Early denoising steps establish global topology (which chains are near each other, overall fold); late steps refine precise geometry (exact contacts, hydrogen bonds, side-chain packing). Steering should match this natural progression:

- **Early steps (t=T → ~T/2)**: Embedding-space modifications are robust and capture high-level biological signals (binding compatibility, fold-level interactions). They tolerate out-of-distribution inputs well.
- **Late steps (t=~T/2 → 0)**: Coordinate-space potentials provide precise geometric control (exact contact distances, dihedral angles, steric compatibility). They are brittle for global topology but excellent for local refinement.

By layering both, we get robust global guidance AND precise local optimization.

---

## 3. Architecture Analysis

### 3.1 Existing Data Flow

```
Trunk (runs once per prediction)
  ├── MSA Module → z (pair repr) [batch, n_tokens, n_tokens, token_z]
  ├── Pairformer → z updated, s (single repr) [batch, n_tokens, token_s]
  └── DiffusionConditioning.forward(s, z, feats)
        ├── PairwiseConditioning(z) → z_conditioned
        ├── [CDR3 BETA SCALING]: z_conditioned *= (1 + β) on CDR3 pairs  ← EMBEDDING SPACE
        ├── AtomEncoder(z_conditioned) → q, c, p, to_keys
        ├── token_trans_bias = proj(z_conditioned)  ← pair info enters diffusion
        ├── atom_enc_bias = proj(p)
        └── atom_dec_bias = proj(p)
            → Returns: {q, c, to_keys, atom_enc_bias, atom_dec_bias, token_trans_bias}

Diffusion Loop (runs N denoising steps)
  For each step (sigma_tm, sigma_t, gamma):
    ├── DiffusionModule.forward(coords, sigma, conditioning)
    │     ├── SingleConditioning(sigma, s_trunk, s_inputs) → s_conditioned
    │     ├── AtomAttentionEncoder(coords, q, c, atom_enc_bias) → a
    │     ├── TokenTransformer(a, s_conditioned, token_trans_bias) → a  ← PAIR BIAS USED HERE
    │     └── AtomAttentionDecoder(a, q, c, atom_dec_bias) → r_update
    │         → Returns: denoised coordinates
    ├── [COORDINATE STEERING]: Gradient descent on energy potentials  ← COORDINATE SPACE
    │     ├── CDR3ConformationPotential (dihedral angles)
    │     ├── AntigenOrientationPotential (CDR-antigen distances)
    │     └── Physical potentials (VDW, bonds, chirality, etc.)
    └── [FK RESAMPLING]: Importance-weight particles by energy
```

### 3.2 Critical Observations

1. **Conditioning is computed ONCE** (before the diffusion loop) and reused at every step.
2. **β-scaling currently applies a static β** throughout all denoising steps—it cannot vary per-step.
3. **Coordinate potentials are already per-step** with time-dependent schedules (`PiecewiseStepFunction`, `ExponentialInterpolation`).
4. **`steering_t`** is already computed per step: `steering_t = 1.0 - (step_idx / num_sampling_steps)` (1.0 at start → 0.0 at end).

### 3.3 What Needs to Change

For hierarchical steering, we need **time-varying embedding influence**. Since conditioning is computed once, we use a **dual-conditioning blending** approach:

- Compute conditioning twice: once with β-scaling (embedding-steered), once without (neutral).
- At each denoising step, blend between the two based on a time-dependent weight α(t).
- α(t) = 1.0 early (full embedding steering) → 0.0 late (no embedding steering, coordinate potentials take over).

---

## 4. Implementation Design

### 4.1 New Constraint Type: `hierarchical_steering`

#### YAML Configuration

```yaml
version: 2
sequences:
  - protein:
      id: [H]
      sequence: "EVQLVESGG..."  # Heavy chain
  - protein:
      id: [L]
      sequence: "DIQMTQSPS..."  # Light chain
  - protein:
      id: [A]
      sequence: "MFVFLVLLP..."  # Antigen
constraints:
  - hierarchical_steering:
      # --- CDR regions (user-provided, sequence-dependent) ---
      cdr_regions:
        - chain: H
          start_res: 97       # CDR-H3 start (1-indexed)
          end_res: 115        # CDR-H3 end (1-indexed)
          label: CDR-H3       # Optional label for logging
        - chain: L
          start_res: 88       # CDR-L3 start (1-indexed)
          end_res: 100        # CDR-L3 end (1-indexed)
          label: CDR-L3
        - chain: H
          start_res: 31       # CDR-H1 (optional, for broader interface steering)
          end_res: 35
          label: CDR-H1
        - chain: H
          start_res: 50       # CDR-H2
          end_res: 65
          label: CDR-H2

      # --- Antigen chain ---
      antigen_chain: A

      # --- Embedding stage parameters ---
      embedding_beta: 1.5           # β value for pair representation scaling (default: 1.5)
      embedding_schedule: "cosine"  # How embedding weight decays: "cosine", "linear", "step"
      embedding_end_fraction: 0.5   # Fraction of diffusion where embedding steering ends (default: 0.5 = halfway)

      # --- Coordinate stage parameters ---
      coordinate_start_fraction: 0.3  # Fraction where coordinate steering begins (default: 0.3)
      contact_threshold: 8.0          # Angstrom threshold for CDR-antigen contacts
      coordinate_guidance_weight: 1.0 # Strength multiplier for coordinate potentials

      # --- Transition smoothness ---
      overlap_fraction: 0.2          # Overlap region where both methods active (default: 0.2)
                                     # overlap = embedding_end - coordinate_start
                                     # In this config: 0.5 - 0.3 = 0.2

      # --- Optional: CDR3 conformation targets (for coordinate stage) ---
      cdr3_conformation:
        chain: H
        start_res: 97
        end_res: 115
        conformation: extended       # "extended", "compact", "kinked", or "custom"

      force: true
```

**Key design decisions:**
- All residue indices are **user-provided** (not hardcoded), since they vary across sequences.
- The `cdr_regions` list supports arbitrary CDR definitions (H1/H2/H3/L1/L2/L3 or any subset).
- `embedding_end_fraction` and `coordinate_start_fraction` create an overlap zone where both methods blend.

### 4.2 Time-Dependent Blending Schedules

The steering weight at diffusion progress `t` (where `t` goes from 1.0 at start to 0.0 at end):

```
α_embed(t) = embedding weight    (1.0 early → 0.0 late)
α_coord(t) = coordinate weight   (0.0 early → 1.0 late)
```

**Schedule options:**

```
Linear:
  α_embed(t) = clamp((t - embedding_end) / (1.0 - embedding_end), 0, 1)

Cosine:
  α_embed(t) = 0.5 * (1 + cos(π * clamp((1-t) / embedding_end_fraction, 0, 1)))

Step:
  α_embed(t) = 1.0 if t > embedding_end_fraction else 0.0
```

Visual diagram (t goes from 1.0→0.0 left to right, matching `steering_t` in diffusion loop):

```
steering_t:  1.0                    0.5                    0.0
             |======== EARLY ========|======== LATE ========|

α_embed:     ████████████████████░░░░░░░░░░░░░░░░░░░░░░░░░░  (strong early → zero)
α_coord:     ░░░░░░░░░░░░░░░░░░░░░░░████████████████████████  (zero → strong late)

             |--- embedding only ---|-- overlap --|-- coordinate only --|
```

### 4.3 Dual-Conditioning Blending (Core Mechanism)

**Concept**: Compute diffusion conditioning twice (cheap compared to full trunk), then blend per-step.

```python
# In boltz2.py forward(), before diffusion sampling:

# Conditioning WITH embedding steering (β-scaled pair representations)
feats_with_beta = feats.copy()
feats_with_beta["cdr3_token_mask"] = hierarchical_mask
feats_with_beta["cdr3_beta_value"] = torch.tensor(embedding_beta)

conditioning_steered = self.diffusion_conditioning(
    s_trunk=s, z_trunk=z,
    relative_position_encoding=relative_position_encoding,
    feats=feats_with_beta,
)

# Conditioning WITHOUT embedding steering (baseline)
feats_neutral = feats.copy()
# Remove beta features or set beta=0
if "cdr3_token_mask" in feats_neutral:
    del feats_neutral["cdr3_token_mask"]
if "cdr3_beta_value" in feats_neutral:
    del feats_neutral["cdr3_beta_value"]

conditioning_neutral = self.diffusion_conditioning(
    s_trunk=s, z_trunk=z,
    relative_position_encoding=relative_position_encoding,
    feats=feats_neutral,
)
```

**Per-step blending inside `AtomDiffusion.sample()`:**

```python
# At each denoising step, blend conditioning based on steering_t
alpha = compute_embedding_weight(steering_t, schedule, embedding_end_fraction)

blended_conditioning = {}
for key in conditioning_steered:
    blended_conditioning[key] = (
        alpha * conditioning_steered[key] +
        (1.0 - alpha) * conditioning_neutral[key]
    )

# Use blended_conditioning for this step's network forward pass
atom_coords_denoised_chunk = self.preconditioned_network_forward(
    atom_coords_noisy[sample_ids_chunk],
    t_hat,
    network_condition_kwargs=dict(
        multiplicity=sample_ids_chunk.numel(),
        diffusion_conditioning=blended_conditioning,  # <-- blended
        **other_kwargs,
    ),
)
```

### 4.4 Coordinate-Stage Potentials with Hierarchical Schedule

The existing coordinate potentials already support `PiecewiseStepFunction` schedules. For hierarchical steering, we modify the schedule so potentials ramp up as embedding steering ramps down:

```python
# In get_potentials(), when hierarchical_steering is enabled:

# AntigenOrientationPotential with hierarchical schedule
AntigenOrientationPotential(
    parameters={
        "guidance_interval": 2,
        "guidance_weight": PiecewiseStepFunction(
            # steering_t goes 1.0 → 0.0
            # coordinate_start_fraction=0.3 means potentials start at steering_t=0.7
            thresholds=[0.7, 0.5, 0.3],
            values=[0.0, 0.3, 0.8, 1.5]  # ramp up: off → weak → medium → strong
        ),
        "resampling_weight": PiecewiseStepFunction(
            thresholds=[0.7, 0.5],
            values=[0.0, 0.5, 1.0]  # resampling also ramps up
        ),
        "union_lambda": ExponentialInterpolation(
            start=8.0, end=0.0, alpha=-2.0
        ),
    }
)

# CDR3ConformationPotential (if cdr3_conformation target specified)
CDR3ConformationPotential(
    parameters={
        "guidance_interval": 2,
        "guidance_weight": PiecewiseStepFunction(
            thresholds=[0.7, 0.5],
            values=[0.0, 0.5, 1.0]  # only active in coordinate stage
        ),
        "resampling_weight": PiecewiseStepFunction(
            thresholds=[0.5],
            values=[0.0, 1.0]
        ),
        "buffer": 0.35,
    }
)
```

---

## 5. File-by-File Implementation Plan

### 5.1 `src/boltz/data/parse/schema.py` — Parse YAML Constraint

**Location**: After the `cdr3_beta_scaling` constraint parsing block (~line 1670-1703).

**Add new constraint type** `hierarchical_steering`:

```python
elif "hierarchical_steering" in constraint:
    hs = constraint["hierarchical_steering"]

    # Parse CDR regions (required)
    cdr_regions = []
    for region in hs["cdr_regions"]:
        cdr_chain_id = region["chain"]
        cdr_start = region["start_res"] - 1  # Convert to 0-indexed
        cdr_end = region["end_res"] - 1      # Convert to 0-indexed
        cdr_label = region.get("label", "")
        cdr_regions.append((cdr_chain_id, cdr_start, cdr_end, cdr_label))

    # Parse antigen chain (required)
    antigen_chain = hs["antigen_chain"]

    # Embedding stage parameters
    embedding_beta = hs.get("embedding_beta", 1.5)
    embedding_schedule = hs.get("embedding_schedule", "cosine")
    embedding_end_fraction = hs.get("embedding_end_fraction", 0.5)

    # Coordinate stage parameters
    coordinate_start_fraction = hs.get("coordinate_start_fraction", 0.3)
    contact_threshold = hs.get("contact_threshold", 8.0)
    coordinate_guidance_weight = hs.get("coordinate_guidance_weight", 1.0)

    # Transition
    overlap_fraction = hs.get("overlap_fraction", 0.2)

    # Optional CDR3 conformation target
    cdr3_conformation_target = None
    if "cdr3_conformation" in hs:
        cc = hs["cdr3_conformation"]
        cdr3_conformation_target = {
            "chain": cc["chain"],
            "start_res": cc["start_res"] - 1,
            "end_res": cc["end_res"] - 1,
            "conformation": cc.get("conformation", "extended"),
        }

    force = hs.get("force", True)

    hierarchical_steering_constraints.append({
        "cdr_regions": cdr_regions,
        "antigen_chain": antigen_chain,
        "embedding_beta": embedding_beta,
        "embedding_schedule": embedding_schedule,
        "embedding_end_fraction": embedding_end_fraction,
        "coordinate_start_fraction": coordinate_start_fraction,
        "contact_threshold": contact_threshold,
        "coordinate_guidance_weight": coordinate_guidance_weight,
        "overlap_fraction": overlap_fraction,
        "cdr3_conformation_target": cdr3_conformation_target,
        "force": force,
    })
```

**Add to `InferenceOptions`** dataclass:

```python
@dataclass
class InferenceOptions:
    pocket_constraints: ... = None
    contact_constraints: ... = None
    cdr3_constraints: ... = None
    antigen_orientation_constraints: ... = None
    cdr3_beta_constraints: ... = None
    hierarchical_steering_constraints: ... = None  # NEW
```

### 5.2 `src/boltz/data/parse/featurizerv2.py` — Generate Tensor Features

**Add function** `process_hierarchical_steering_constraints()`:

```python
def process_hierarchical_steering_constraints(
    record,
    tokenized,
    token_array,
    chain_id_to_idx,
):
    """
    Convert hierarchical steering YAML constraints into tensor features.

    Produces:
      - hierarchical_cdr_token_mask: bool [N_tokens] — True for all CDR residues
      - hierarchical_antigen_token_mask: bool [N_tokens] — True for antigen residues
      - hierarchical_cdr_atom_indices: [M] — flat atom indices for CDR atoms
      - hierarchical_antigen_atom_indices: [K] — flat atom indices for antigen atoms
      - hierarchical_embedding_beta: float scalar — β value for pair scaling
      - hierarchical_embedding_schedule: str — "cosine", "linear", or "step"
      - hierarchical_embedding_end_fraction: float scalar
      - hierarchical_coordinate_start_fraction: float scalar
      - hierarchical_contact_threshold: float scalar
      - hierarchical_coordinate_guidance_weight: float scalar
    """
    constraints = record.inference_options.hierarchical_steering_constraints
    if not constraints:
        return {}

    hs = constraints[0]  # Take the first (primary) hierarchical constraint

    n_tokens = len(token_array)
    cdr_token_mask = torch.zeros(n_tokens, dtype=torch.bool)
    antigen_token_mask = torch.zeros(n_tokens, dtype=torch.bool)

    # Build CDR token mask from all specified CDR regions
    for (chain_id, start_res, end_res, label) in hs["cdr_regions"]:
        chain_idx = chain_id_to_idx[chain_id]
        for tok_idx in range(n_tokens):
            tok = token_array[tok_idx]
            if tok.chain_idx == chain_idx and start_res <= tok.res_idx <= end_res:
                cdr_token_mask[tok_idx] = True

    # Build antigen token mask
    antigen_chain_idx = chain_id_to_idx[hs["antigen_chain"]]
    for tok_idx in range(n_tokens):
        tok = token_array[tok_idx]
        if tok.chain_idx == antigen_chain_idx:
            antigen_token_mask[tok_idx] = True

    # Build CDR and antigen atom index lists (for coordinate-stage potentials)
    cdr_atom_indices = []
    antigen_atom_indices = []
    for tok_idx in range(n_tokens):
        atom_start = token_array[tok_idx].atom_start
        atom_end = token_array[tok_idx].atom_end
        if cdr_token_mask[tok_idx]:
            cdr_atom_indices.extend(range(atom_start, atom_end))
        if antigen_token_mask[tok_idx]:
            antigen_atom_indices.extend(range(atom_start, atom_end))

    feats = {
        "hierarchical_cdr_token_mask": cdr_token_mask,
        "hierarchical_antigen_token_mask": antigen_token_mask,
        "hierarchical_cdr_atom_indices": torch.tensor(cdr_atom_indices, dtype=torch.long),
        "hierarchical_antigen_atom_indices": torch.tensor(antigen_atom_indices, dtype=torch.long),
        "hierarchical_embedding_beta": torch.tensor(hs["embedding_beta"], dtype=torch.float32),
        "hierarchical_embedding_schedule": hs["embedding_schedule"],
        "hierarchical_embedding_end_fraction": torch.tensor(hs["embedding_end_fraction"], dtype=torch.float32),
        "hierarchical_coordinate_start_fraction": torch.tensor(hs["coordinate_start_fraction"], dtype=torch.float32),
        "hierarchical_contact_threshold": torch.tensor(hs["contact_threshold"], dtype=torch.float32),
        "hierarchical_coordinate_guidance_weight": torch.tensor(hs["coordinate_guidance_weight"], dtype=torch.float32),
    }

    # Optional: CDR3 conformation target
    if hs.get("cdr3_conformation_target"):
        ct = hs["cdr3_conformation_target"]
        # Reuse existing CDR3 conformation feature generation logic
        feats["hierarchical_cdr3_conformation"] = ct

    return feats
```

**Call this function** in the main featurization pipeline (where `process_cdr3_beta_constraints` is called), and merge the returned dict into `feats`.

### 5.3 `src/boltz/model/modules/diffusion_conditioning.py` — Hierarchical β-Scaling

**Current code** (lines 95-118) applies static β. For hierarchical steering, we modify this to also produce an **unscaled** version of the conditioning, OR we handle the blending at a higher level.

**Recommended approach**: No changes needed to `DiffusionConditioning` itself. Instead, the dual-conditioning logic lives in `boltz2.py` (section 5.5 below), where we call `diffusion_conditioning()` twice: once with β features, once without.

However, we must extend the mask to cover **all CDR regions** (not just CDR3-CDR3 pairs). For hierarchical steering, the embedding mask should include CDR-antigen cross-pairs:

```python
# In DiffusionConditioning.forward(), extend the existing beta-scaling block:

if "hierarchical_cdr_token_mask" in feats and "hierarchical_embedding_beta" in feats:
    cdr_mask = feats["hierarchical_cdr_token_mask"].to(z.device).to(torch.bool)
    ag_mask = feats["hierarchical_antigen_token_mask"].to(z.device).to(torch.bool)
    beta_val = float(feats["hierarchical_embedding_beta"].item())

    if cdr_mask.dim() == 2:
        cdr_mask = cdr_mask.squeeze(0)
    if ag_mask.dim() == 2:
        ag_mask = ag_mask.squeeze(0)

    if abs(beta_val) > 1e-6:
        # Scale CDR-CDR pairs (intra-CDR interactions)
        cdr_pair_mask = cdr_mask.unsqueeze(-1) & cdr_mask.unsqueeze(-2)

        # Scale CDR-Antigen cross-pairs (binding interface)
        cross_pair_mask = (
            (cdr_mask.unsqueeze(-1) & ag_mask.unsqueeze(-2)) |
            (ag_mask.unsqueeze(-1) & cdr_mask.unsqueeze(-2))
        )

        # Combined mask: both intra-CDR and CDR-antigen interface
        combined_mask = cdr_pair_mask | cross_pair_mask

        scaling_tensor = combined_mask.unsqueeze(0).unsqueeze(-1).float()
        scaling_factor = 1.0 + beta_val * scaling_tensor
        z = z * scaling_factor
```

### 5.4 `src/boltz/model/potentials/schedules.py` — Add Cosine Schedule

Add a new schedule class for smooth cosine transitions:

```python
class CosineDecaySchedule(ParameterSchedule):
    """Cosine decay from `start_val` to `end_val` between t_start and t_end.

    t is steering_t (1.0 at beginning of diffusion → 0.0 at end).
    """
    def __init__(self, start_val, end_val, t_start, t_end):
        self.start_val = start_val
        self.end_val = end_val
        self.t_start = t_start   # steering_t where decay begins (e.g., 0.7)
        self.t_end = t_end       # steering_t where decay ends (e.g., 0.3)

    def compute(self, t):
        if t >= self.t_start:
            return self.start_val
        elif t <= self.t_end:
            return self.end_val
        else:
            # Cosine interpolation in the transition zone
            progress = (self.t_start - t) / (self.t_start - self.t_end)
            return self.end_val + 0.5 * (self.start_val - self.end_val) * (1 + math.cos(math.pi * progress))


class LinearRampSchedule(ParameterSchedule):
    """Linear ramp from start_val to end_val between t_start and t_end."""
    def __init__(self, start_val, end_val, t_start, t_end):
        self.start_val = start_val
        self.end_val = end_val
        self.t_start = t_start
        self.t_end = t_end

    def compute(self, t):
        if t >= self.t_start:
            return self.start_val
        elif t <= self.t_end:
            return self.end_val
        else:
            progress = (self.t_start - t) / (self.t_start - self.t_end)
            return self.start_val + (self.end_val - self.start_val) * progress
```

### 5.5 `src/boltz/model/models/boltz2.py` — Dual-Conditioning Logic

**Location**: In the `forward()` method, after trunk processing and before `structure_module.sample()` (~line 534-562).

```python
# After computing trunk representations s, z:

hierarchical_enabled = (
    self.steering_args.get("hierarchical_steering", False) and
    "hierarchical_cdr_token_mask" in feats
)

if hierarchical_enabled:
    # --- Compute STEERED conditioning (with embedding β) ---
    # feats already contains hierarchical_cdr_token_mask and hierarchical_embedding_beta
    # The DiffusionConditioning.forward() will apply β-scaling
    conditioning_steered = self.diffusion_conditioning(
        s_trunk=s,
        z_trunk=z,
        relative_position_encoding=relative_position_encoding,
        feats=feats,
    )

    # --- Compute NEUTRAL conditioning (no β) ---
    feats_neutral = dict(feats)
    # Temporarily remove hierarchical beta features
    saved_beta = feats_neutral.pop("hierarchical_embedding_beta", None)
    saved_cdr_mask = feats_neutral.pop("hierarchical_cdr_token_mask", None)
    saved_ag_mask = feats_neutral.pop("hierarchical_antigen_token_mask", None)

    conditioning_neutral = self.diffusion_conditioning(
        s_trunk=s,
        z_trunk=z,
        relative_position_encoding=relative_position_encoding,
        feats=feats_neutral,
    )

    # Restore feats for downstream use
    feats["hierarchical_embedding_beta"] = saved_beta
    feats["hierarchical_cdr_token_mask"] = saved_cdr_mask
    feats["hierarchical_antigen_token_mask"] = saved_ag_mask

    # Pass both conditionings to the diffusion sampler
    diffusion_conditioning = {
        "steered": conditioning_steered,
        "neutral": conditioning_neutral,
    }
else:
    # Standard single conditioning (existing behavior)
    q, c, to_keys, atom_enc_bias, atom_dec_bias, token_trans_bias = (
        self.diffusion_conditioning(
            s_trunk=s, z_trunk=z,
            relative_position_encoding=relative_position_encoding,
            feats=feats,
        )
    )
    diffusion_conditioning = {
        "q": q, "c": c, "to_keys": to_keys,
        "atom_enc_bias": atom_enc_bias,
        "atom_dec_bias": atom_dec_bias,
        "token_trans_bias": token_trans_bias,
    }

# Pass to sample() as before
struct_out = self.structure_module.sample(
    s_trunk=s.float(),
    s_inputs=s_inputs.float(),
    feats=feats,
    num_sampling_steps=num_sampling_steps,
    atom_mask=feats["atom_pad_mask"].float(),
    multiplicity=diffusion_samples,
    max_parallel_samples=max_parallel_samples,
    steering_args=self.steering_args,
    diffusion_conditioning=diffusion_conditioning,
)
```

### 5.6 `src/boltz/model/modules/diffusionv2.py` — Per-Step Blending in Diffusion Loop

**Location**: Inside `AtomDiffusion.sample()`, modify the denoising loop (~line 350-530).

**Add blending logic** before the network forward pass:

```python
def sample(self, atom_mask, num_sampling_steps=None, multiplicity=1,
           max_parallel_samples=None, steering_args=None,
           **network_condition_kwargs):

    # Detect hierarchical steering mode
    hierarchical_mode = (
        steering_args is not None and
        steering_args.get("hierarchical_steering", False) and
        "steered" in network_condition_kwargs.get("diffusion_conditioning", {})
    )

    if hierarchical_mode:
        conditioning_steered = network_condition_kwargs["diffusion_conditioning"]["steered"]
        conditioning_neutral = network_condition_kwargs["diffusion_conditioning"]["neutral"]
        feats = network_condition_kwargs["feats"]

        # Read schedule parameters from feats
        embedding_schedule = feats.get("hierarchical_embedding_schedule", "cosine")
        embedding_end_fraction = float(feats["hierarchical_embedding_end_fraction"].item())

        # Build the embedding weight schedule
        if embedding_schedule == "cosine":
            embed_schedule = CosineDecaySchedule(
                start_val=1.0, end_val=0.0,
                t_start=1.0, t_end=1.0 - embedding_end_fraction
            )
        elif embedding_schedule == "linear":
            embed_schedule = LinearRampSchedule(
                start_val=1.0, end_val=0.0,
                t_start=1.0, t_end=1.0 - embedding_end_fraction
            )
        else:  # "step"
            embed_schedule = PiecewiseStepFunction(
                thresholds=[1.0 - embedding_end_fraction],
                values=[1.0, 0.0]
            )

    # ... existing initialization code ...

    # Inside the denoising loop:
    for step_idx, (sigma_tm, sigma_t, gamma) in enumerate(sigmas_and_gammas):
        # ... existing augmentation code ...

        steering_t = 1.0 - (step_idx / num_sampling_steps)

        # --- HIERARCHICAL BLENDING ---
        if hierarchical_mode:
            alpha = embed_schedule.compute(steering_t)

            # Blend conditioning tensors
            blended_conditioning = {}
            for key in conditioning_steered:
                blended_conditioning[key] = (
                    alpha * conditioning_steered[key] +
                    (1.0 - alpha) * conditioning_neutral[key]
                )

            # Replace diffusion_conditioning in kwargs for this step
            step_kwargs = dict(network_condition_kwargs)
            step_kwargs["diffusion_conditioning"] = blended_conditioning
        else:
            step_kwargs = network_condition_kwargs

        # ... network forward pass uses step_kwargs ...
        with torch.no_grad():
            for sample_ids_chunk in sample_ids_chunks:
                atom_coords_denoised_chunk = self.preconditioned_network_forward(
                    atom_coords_noisy[sample_ids_chunk],
                    t_hat,
                    network_condition_kwargs=dict(
                        multiplicity=sample_ids_chunk.numel(),
                        **{k: v for k, v in step_kwargs.items()
                           if k != "diffusion_conditioning"},
                        diffusion_conditioning=blended_conditioning if hierarchical_mode
                            else step_kwargs["diffusion_conditioning"],
                    ),
                )
                atom_coords_denoised[sample_ids_chunk] = atom_coords_denoised_chunk

        # ... existing FK resampling and coordinate guidance code ...
        # The coordinate potentials already have their own time-dependent schedules
        # that ramp UP as embedding steering ramps DOWN
```

### 5.7 `src/boltz/model/potentials/potentials.py` — Hierarchical Potential Setup

**Add to `get_potentials()`** (~after line 956):

```python
# Add hierarchical steering potentials
if boltz2 and steering_args.get("hierarchical_steering", False):
    feats = steering_args.get("feats", {})
    coord_start = float(feats.get("hierarchical_coordinate_start_fraction", 0.3))
    coord_weight = float(feats.get("hierarchical_coordinate_guidance_weight", 1.0))

    # steering_t threshold: coordinate_start_fraction=0.3 means
    # coordinate potentials activate at steering_t = 0.7 (= 1.0 - 0.3)
    t_activate = 1.0 - coord_start

    # Antigen orientation potential (CDR-antigen contact optimization)
    potentials.append(
        AntigenOrientationPotential(
            parameters={
                "guidance_interval": 2,
                "guidance_weight": PiecewiseStepFunction(
                    thresholds=[t_activate, t_activate - 0.2, t_activate - 0.4],
                    values=[0.0, 0.3 * coord_weight, 0.8 * coord_weight, 1.5 * coord_weight]
                ),
                "resampling_weight": PiecewiseStepFunction(
                    thresholds=[t_activate, t_activate - 0.2],
                    values=[0.0, 0.5, 1.0]
                ),
                "union_lambda": ExponentialInterpolation(
                    start=8.0, end=0.0, alpha=-2.0
                ),
            }
        )
    )

    # Optional: CDR3 conformation potential if target specified
    if steering_args.get("hierarchical_cdr3_conformation"):
        potentials.append(
            CDR3ConformationPotential(
                parameters={
                    "guidance_interval": 2,
                    "guidance_weight": PiecewiseStepFunction(
                        thresholds=[t_activate],
                        values=[0.0, 0.8 * coord_weight]
                    ),
                    "resampling_weight": PiecewiseStepFunction(
                        thresholds=[t_activate - 0.1],
                        values=[0.0, 1.0]
                    ),
                    "buffer": 0.35,
                }
            )
        )
```

### 5.8 `src/boltz/main.py` — CLI Integration

**Add CLI flag** (~line 989-1010):

```python
@click.option(
    "--hierarchical_steering",
    is_flag=True,
    default=False,
    help="Enable hierarchical steering: embedding-space (early) + coordinate-space (late). "
         "Requires hierarchical_steering constraint in YAML config.",
)
```

**Add to `BoltzSteeringParams`** (~line 165-176):

```python
@dataclass
class BoltzSteeringParams:
    fk_steering: bool = False
    num_particles: int = 3
    fk_lambda: float = 4.0
    fk_resampling_interval: int = 3
    physical_guidance_update: bool = False
    contact_guidance_update: bool = True
    cdr3_steering: bool = False
    antigen_steering: bool = False
    hierarchical_steering: bool = False  # NEW
    num_gd_steps: int = 20
```

**Wire up the flag** (~line 1357-1372):

```python
steering_args.hierarchical_steering = hierarchical_steering

# Hierarchical steering implies both FK and coordinate potentials
if hierarchical_steering:
    steering_args.fk_steering = True
    steering_args.physical_guidance_update = True
    steering_args.contact_guidance_update = True
```

**Validation**:

```python
if hierarchical_steering and not use_potentials:
    raise ValueError("--hierarchical_steering requires --use_potentials")
```

---

## 6. Detailed Algorithm: What Happens at Each Diffusion Step

```
Step 0 (steering_t ≈ 1.0): EARLY - Full embedding steering
  ├─ α_embed = 1.0 → conditioning = 100% steered (β-scaled CDR+antigen pairs)
  ├─ α_coord = 0.0 → no coordinate potentials active
  ├─ Effect: Global topology guided. CDR regions emphasized in pair attention.
  │   Antigen positioned broadly toward CDR neighborhood.
  └─ Robustness: High. Works even for novel sequences.

Steps 1..N/4 (steering_t ≈ 0.75): EARLY - Mostly embedding
  ├─ α_embed ≈ 0.9
  ├─ α_coord ≈ 0.0
  └─ Global fold solidifying under embedding guidance

Steps N/4..N/2 (steering_t ≈ 0.5): TRANSITION - Both active
  ├─ α_embed ≈ 0.5 (cosine decay) → conditioning is 50/50 blend
  ├─ α_coord ≈ 0.3-0.8 → coordinate potentials ramping up
  ├─ Effect: Embedding provides stability, coordinates begin refining contacts
  └─ Smooth handoff prevents discontinuities

Steps N/2..3N/4 (steering_t ≈ 0.25): LATE - Mostly coordinate
  ├─ α_embed ≈ 0.1 → minimal embedding influence
  ├─ α_coord ≈ 1.0 → full coordinate potential strength
  ├─ AntigenOrientationPotential: optimizing CDR-antigen distances
  ├─ CDR3ConformationPotential: enforcing dihedral targets (if specified)
  └─ Physical potentials: VDW, chirality, bonds

Steps 3N/4..N (steering_t ≈ 0.0): FINAL - Coordinate only
  ├─ α_embed = 0.0 → conditioning = 100% neutral (no β)
  ├─ α_coord = 1.0 → full coordinate refinement
  └─ Fine-grained geometry: hydrogen bonds, exact contacts, side-chain packing
```

---

## 7. Cost Analysis

### 7.1 Additional Computation

| Component | Cost | When |
|-----------|------|------|
| Second `DiffusionConditioning.forward()` | ~5-8% of trunk | Once, before diffusion loop |
| Per-step tensor blending (6 tensors, linear interp) | <0.1% per step | Every denoising step |
| Coordinate potentials (existing) | ~10-20% per step | Active from mid-diffusion onward |
| FK resampling (existing) | ~5-10% per step | Every `fk_resampling_interval` steps |

**Total overhead**: ~15-25% over baseline (no steering), which is comparable to existing steering methods. The dual conditioning adds negligible cost since it's just one extra forward pass through the conditioning module (not the full trunk or diffusion model).

### 7.2 Memory

- Storing two sets of conditioning tensors doubles conditioning memory: ~2x for `{q, c, to_keys, atom_enc_bias, atom_dec_bias, token_trans_bias}`.
- These are small relative to model activations (~5-10% of total memory).
- The blended conditioning can be computed in-place to avoid storing a third copy.

---

## 8. Configuration Examples

### 8.1 Default Hierarchical Steering (Recommended Starting Point)

```yaml
constraints:
  - hierarchical_steering:
      cdr_regions:
        - chain: H
          start_res: 97
          end_res: 115
          label: CDR-H3
        - chain: L
          start_res: 88
          end_res: 100
          label: CDR-L3
      antigen_chain: A
      embedding_beta: 1.5
      embedding_schedule: cosine
      embedding_end_fraction: 0.5
      coordinate_start_fraction: 0.3
      contact_threshold: 8.0
      force: true
```

```bash
boltz predict config.yaml --use_potentials --hierarchical_steering
```

### 8.2 Aggressive Embedding (Novel/OOD Antibodies)

For synthetic or unusual antibodies where robustness matters most:

```yaml
constraints:
  - hierarchical_steering:
      cdr_regions:
        - chain: H
          start_res: 97
          end_res: 115
        - chain: L
          start_res: 88
          end_res: 100
      antigen_chain: A
      embedding_beta: 2.5              # Stronger embedding influence
      embedding_schedule: cosine
      embedding_end_fraction: 0.7      # Embedding active for 70% of diffusion
      coordinate_start_fraction: 0.5   # Coordinate starts later
      contact_threshold: 10.0          # Wider contact radius
      force: true
```

### 8.3 Precision-Focused (Known Structure with Refinement Target)

When you have a good initial prediction and want precise contacts:

```yaml
constraints:
  - hierarchical_steering:
      cdr_regions:
        - chain: H
          start_res: 97
          end_res: 115
        - chain: H
          start_res: 31
          end_res: 35
        - chain: H
          start_res: 50
          end_res: 65
        - chain: L
          start_res: 88
          end_res: 100
      antigen_chain: A
      embedding_beta: 0.8             # Lighter embedding
      embedding_schedule: linear
      embedding_end_fraction: 0.3     # Short embedding phase
      coordinate_start_fraction: 0.2  # Coordinate starts early
      coordinate_guidance_weight: 1.5 # Stronger coordinate potentials
      contact_threshold: 6.0          # Tighter contacts
      cdr3_conformation:
        chain: H
        start_res: 97
        end_res: 115
        conformation: extended
      force: true
```

### 8.4 With FK Particles for Diverse Sampling

```bash
boltz predict config.yaml \
  --use_potentials \
  --hierarchical_steering \
  --num_particles 5 \
  --diffusion_samples 3
```

---

## 9. Testing Strategy

### 9.1 Unit Tests

```python
class TestHierarchicalSteering:

    def test_yaml_parsing(self):
        """Verify hierarchical_steering constraint is parsed correctly."""
        config = parse_yaml("test_hierarchical.yaml")
        hs = config.inference_options.hierarchical_steering_constraints[0]
        assert len(hs["cdr_regions"]) >= 2
        assert hs["embedding_beta"] == 1.5
        assert hs["antigen_chain"] == "A"

    def test_feature_generation(self):
        """Verify tensor features are generated correctly."""
        feats = process_hierarchical_steering_constraints(record, tokenized, ...)
        assert feats["hierarchical_cdr_token_mask"].dtype == torch.bool
        assert feats["hierarchical_cdr_token_mask"].sum() > 0
        assert feats["hierarchical_antigen_token_mask"].sum() > 0

    def test_cosine_schedule(self):
        """Verify cosine schedule produces correct values."""
        schedule = CosineDecaySchedule(1.0, 0.0, t_start=1.0, t_end=0.5)
        assert schedule.compute(1.0) == 1.0      # Full strength at start
        assert schedule.compute(0.75) == pytest.approx(0.5, abs=0.05)  # ~50% at midpoint
        assert schedule.compute(0.5) == 0.0       # Zero at end
        assert schedule.compute(0.0) == 0.0       # Zero after end

    def test_dual_conditioning_shapes(self):
        """Verify both conditioning dicts have matching shapes."""
        cond_steered = diffusion_conditioning(s, z, feats_with_beta)
        cond_neutral = diffusion_conditioning(s, z, feats_no_beta)
        for key in cond_steered:
            assert cond_steered[key].shape == cond_neutral[key].shape

    def test_blending_boundary_conditions(self):
        """Verify blending at alpha=0 and alpha=1 gives expected results."""
        for key in cond_steered:
            # alpha=1: fully steered
            blended = 1.0 * cond_steered[key] + 0.0 * cond_neutral[key]
            assert torch.allclose(blended, cond_steered[key])
            # alpha=0: fully neutral
            blended = 0.0 * cond_steered[key] + 1.0 * cond_neutral[key]
            assert torch.allclose(blended, cond_neutral[key])
```

### 9.2 Integration Tests

```python
def test_hierarchical_vs_baseline():
    """Hierarchical steering should produce different structures than baseline."""
    pred_baseline = boltz_predict(config, use_potentials=False)
    pred_hierarchical = boltz_predict(config, use_potentials=True, hierarchical_steering=True)

    cdr_rmsd = compute_cdr_rmsd(pred_baseline, pred_hierarchical)
    assert cdr_rmsd > 0.5  # Structures should differ

def test_hierarchical_vs_coordinate_only():
    """Hierarchical should be at least as good as coordinate-only steering."""
    pred_coord = boltz_predict(config, use_potentials=True, antigen_steering=True)
    pred_hier = boltz_predict(config, use_potentials=True, hierarchical_steering=True)

    contact_coord = compute_contact_score(pred_coord)
    contact_hier = compute_contact_score(pred_hier)
    # Hierarchical should match or beat coordinate-only
    assert contact_hier >= contact_coord * 0.95

def test_hierarchical_robustness():
    """Hierarchical should be more robust to hyperparameter variation."""
    results = []
    for beta in [0.5, 1.0, 1.5, 2.0, 2.5]:
        pred = boltz_predict(config, hierarchical_steering=True, embedding_beta=beta)
        results.append(compute_contact_score(pred))

    # Variance across beta values should be low (robust)
    assert np.std(results) / np.mean(results) < 0.15  # CV < 15%
```

### 9.3 Evaluation Metrics

| Metric | Description | Target |
|--------|-------------|--------|
| CDR-Antigen Contact Score | Number of CDR-antigen atom pairs < 8Å | Higher is better |
| CDR RMSD (if native available) | Backbone RMSD of CDR regions vs. crystal | < 2.0 Å |
| pLDDT | Predicted local distance difference test | > 60 (CDR), > 70 (framework) |
| Transition Smoothness | Max change in contact score between consecutive steps | No discontinuity at transition |
| Robustness (CV) | Coefficient of variation across β values | < 15% |
| Computational Cost | Wall-clock time vs. baseline | < 1.3x baseline |

---

## 10. Key Files Modified (Summary)

| File | Change Type | Lines Changed (approx) |
|------|-------------|----------------------|
| `src/boltz/data/parse/schema.py` | Add `hierarchical_steering` constraint parser | +50 |
| `src/boltz/data/parse/featurizerv2.py` | Add `process_hierarchical_steering_constraints()` | +80 |
| `src/boltz/model/modules/diffusion_conditioning.py` | Extend β-scaling to CDR-antigen cross-pairs | +20 |
| `src/boltz/model/potentials/schedules.py` | Add `CosineDecaySchedule`, `LinearRampSchedule` | +30 |
| `src/boltz/model/potentials/potentials.py` | Add hierarchical potentials to `get_potentials()` | +40 |
| `src/boltz/model/models/boltz2.py` | Dual-conditioning computation | +40 |
| `src/boltz/model/modules/diffusionv2.py` | Per-step blending in `sample()` | +50 |
| `src/boltz/main.py` | CLI flag + `BoltzSteeringParams` | +15 |
| **Total** | | **~325 lines** |

---

## 11. Potential Risks and Mitigations

| Risk | Impact | Mitigation |
|------|--------|------------|
| Blending introduces discontinuity | Structure artifacts at transition | Cosine schedule ensures smooth transition; test with step count sweep |
| Double conditioning OOM on large systems | Memory pressure | In-place blending; profile memory; fallback to single conditioning with scheduled β via feats update |
| β-scaling on cross-pairs too aggressive | Distorted antigen structure | Cap cross-pair β at 0.7 * embedding_beta; monitor antigen pLDDT |
| Coordinate potentials fight embedding | Conflicting guidance directions | Overlap zone is short (20%); both push toward same goal (CDR-antigen binding) |
| Compiled score model breaks with blended conditioning | torch.compile graph breaks | Ensure blended tensors have same shape/dtype; test compile mode |

---

## 12. Future Extensions

1. **EmbedOpt-style gradient descent in embedding space**: Instead of static β-scaling, compute `∂E/∂z` and optimize pair representations directly. Much more expensive but theoretically superior for out-of-distribution.

2. **Adaptive transition point**: Instead of fixed `embedding_end_fraction`, use contact score monitoring (like N+) to detect when embedding steering has plateaued and automatically switch to coordinate mode.

3. **Multi-objective embedding stage**: During the embedding stage, apply different β values to different CDR-antigen region pairs simultaneously (combining with K+ region-specific scaling).

4. **Integration with experimental data (Idea X)**: The embedding stage naturally accommodates experimental constraints (cryo-EM density, cross-linking) by converting them to embedding-space objectives.

---

**Document Status**: Ready for implementation
**Last Updated**: 2026-03-11
