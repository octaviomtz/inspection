# Implementation Guide: Upgraded Strategies from Round 1

**Date**: 2026-03-18
**Purpose**: Concrete implementation suggestions for the upgraded versions of all Round 1 evaluated strategies, based on the "Next Steps" section of `ROUND_1_FINAL_EVALUATION.md`
**Scope**: Code structure, module changes, integration points — no code changes yet

---

## Architecture Reference

The standard modification pattern across all Round 1 strategies involves these files:

| File | Role |
|---|---|
| `src/boltz/main.py` | CLI flags → `BoltzSteeringParams` dataclass |
| `src/boltz/data/parse/schema.py` | YAML constraint parsing |
| `src/boltz/data/feature/featurizerv2.py` | Produces feature tensors passed in `feats` dict |
| `src/boltz/model/modules/diffusion_conditioning.py` | Applies beta-scaling to token-pair attention biases |
| `src/boltz/model/modules/diffusionv2.py` | Core diffusion loop: FK resampling, gradient guidance |
| `src/boltz/model/potentials/potentials.py` | Potential classes + `get_potentials()` factory |
| `src/boltz/model/models/boltz2.py` | Passes `z_trunk` and other kwargs into diffusion |
| `src/boltz/data/module/inferencev2.py` | Wires constraints/feats through to featurizer |

**Current `BoltzSteeringParams`** (in `main.py` lines 165–178):
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
    num_gd_steps: int = 20
```

**Current beta-scaling** (in `diffusion_conditioning.py` lines 95–137):
```python
# Reads feats["cdr3_token_mask"] and feats["cdr3_beta_value"]
# Applies: token_trans_bias += beta_val * (CDR3_i & CDR3_j) * token_trans_bias
```

**Current FK loop** (in `diffusionv2.py`):
- `steering_t = 1.0 - (step_idx / num_sampling_steps)` — fraction remaining (1.0 at start, 0.0 at end)
- Potentials use `compute_parameters(steering_t)` to get time-varying `guidance_weight` and `resampling_weight`
- `PiecewiseStepFunction` and `ExponentialInterpolation` are the schedule primitives

---

## Strategy A+ — FK Particles v2

**Current**: 20 particles (down from proposed higher counts), FK resampling every 3 steps with `ContactPotential` + `SymmetricChainCOMPotential`. 25% of complexes OOM at 20 particles.

---

### A+.1 — Fix OOM: Adaptive Particle Reduction

**Goal**: Recover the 12/47 missing complexes without running out of memory.

**Approach**: Estimate GPU memory before launch and auto-reduce `num_particles` to fit.

**Files to modify**: `src/boltz/main.py`

**Changes**:
1. Add `auto_reduce_particles: bool = False` to `BoltzSteeringParams`
2. In the predict function (around line 1357), before calling `model.sample()`:
   - Query available GPU memory: `torch.cuda.mem_get_info()`
   - Estimate memory cost: `num_atoms × num_particles × 3 × 4 bytes` (rough)
   - If estimated cost > 80% of available memory: halve `num_particles` (min=5)
   - Log the reduction

```python
# Pseudocode in predict()
if steering_args.auto_reduce_particles and steering_args.fk_steering:
    avail_mem, total_mem = torch.cuda.mem_get_info()
    estimated_cost = num_atoms * steering_args.num_particles * 3 * 4 * safety_factor
    while estimated_cost > 0.8 * avail_mem and steering_args.num_particles > 5:
        steering_args.num_particles //= 2
        estimated_cost //= 2
    log(f"Auto-reduced to {steering_args.num_particles} particles")
```

**New CLI flags**:
```
--auto_reduce_particles   (bool) Auto-reduce particle count if OOM risk detected
--fk_min_particles INT    (default: 5) Floor for auto-reduction
```

**Complexity**: Low — no model changes, pure bookkeeping.

---

### A+.2 — Custom Re-Ranking by Steering Energy

**Goal**: Replace iPTM-only model selection with a composite score that incorporates the FK potential energy, which correlates with docking quality better than raw confidence.

**Problem**: The FK loop computes `energy_traj` but this is currently discarded after resampling. It is never returned to the caller or saved to disk.

**Files to modify**:
- `src/boltz/model/modules/diffusionv2.py`: Return final energy per surviving sample alongside coordinates
- `src/boltz/main.py`: Accept and pass through energy scores
- New: `src/boltz/analysis/reranking.py`: Composite scorer

**Changes**:

1. In `diffusionv2.py sample()`, at the end of the loop, collect final energies of the surviving particles and return them alongside `atom_coords_denoised`.

2. In `reranking.py` (new file):
```python
def composite_rerank(predictions, alpha=0.5, beta=0.3, gamma=0.2):
    """
    predictions: list of dicts with keys:
        'iptm': float
        'steering_energy': float (lower = better)
        'interface_plddt': float
    Returns: sorted list, best first
    """
    for p in predictions:
        p['score'] = (
            alpha * p['iptm']
            - beta * p['steering_energy']   # negate: lower energy is better
            + gamma * p['interface_plddt']
        )
    return sorted(predictions, key=lambda p: p['score'], reverse=True)
```

3. In `main.py`: add `--rerank_alpha`, `--rerank_beta`, `--rerank_gamma` flags.

**YAML format**: No YAML change needed — re-ranking is post-processing.

**Complexity**: Low–Medium. Requires threading energy values through the return path.

---

### A+.3 — Separate FK Seeds from FK Particles

**Goal**: Allow independent restarts (seeds) each with their own particle ensemble, increasing diversity without a quadratic memory cost.

**Current structure**: `multiplicity` (total trajectories) = `num_samples × num_particles`. All particles from all seeds are treated identically by the resampling code (reshaped into `[multiplicity / num_particles, num_particles]`).

**Proposed**: Add explicit `num_fk_seeds` controlling how many independent starting points exist. Each seed has its own particle group that resamples internally. Seeds never cross-resample.

**Files to modify**: `src/boltz/main.py`, `src/boltz/model/modules/diffusionv2.py`

**Changes**:

1. Add `num_fk_seeds: int = 1` to `BoltzSteeringParams`.

2. In `diffusionv2.py sample()`:
   - `multiplicity = num_fk_seeds * num_particles`
   - Reshape `resample_weights` as `[num_fk_seeds, num_particles]` (already correct for 1 seed)
   - Each seed resamples only within its own `num_particles` group
   - Return one best candidate per seed (select by lowest energy within each group)

```python
# Current reshape (works for 1 seed):
resample_weights = F.softmax(...).reshape(-1, steering_args["num_particles"])
# With multiple seeds this is already correct since shape is [num_seeds, num_particles]
```

3. Add `--num_fk_seeds INT` to CLI.

**Complexity**: Low. The data structure already supports this — it's mostly a parameter rename and a final selection step.

---

### A+.4 — Embedding Steering as Particle-Free Alternative

**Goal**: Use beta-scaling (Strategy Y's mechanism) instead of FK particles for complexes that would OOM with particles.

**Approach**: If `auto_reduce_particles` reduces particles to minimum and still OOM, fall back to Y-style beta-scaling with no FK particles (`fk_steering=False`, `cdr3_beta_value=set`).

**Files to modify**: `src/boltz/main.py` (fallback logic), no model changes.

**Complexity**: Low. Pure orchestration change.

---

## Strategy D+ — Canonical Ensemble v2

**Current**: `CDR3ConformationPotential` targets dihedral angles of known CDR canonical conformations. Minimal effect (+0.009 DockQ), confidence degrades.

---

### D+.1 — CDR-Antigen Pair Beta-Scaling

**Goal**: Extend the current CDR3-only beta-scaling to include CDR3–antigen token pairs (interface pairs), which are more directly relevant to docking.

**Current mask** in `diffusion_conditioning.py`:
```python
cdr3_pair_mask_4d = (cdr3_mask_i & cdr3_mask_j)  # CDR3-CDR3 pairs only
```

**Proposed**: Add a second mask for CDR3-antigen pairs with a separate beta value.

**Files to modify**:
- `src/boltz/data/feature/featurizerv2.py`: Produce `antigen_token_mask` tensor alongside `cdr3_token_mask`
- `src/boltz/model/modules/diffusion_conditioning.py`: Compute and apply `cdr3_antigen_pair_mask`
- `src/boltz/main.py`: Add `--cdr3_antigen_beta_value FLOAT` flag

**Changes in `diffusion_conditioning.py`**:
```python
# After existing CDR3-CDR3 mask computation:
if "antigen_token_mask" in feats and "cdr3_antigen_beta_value" in feats:
    ag_mask = feats["antigen_token_mask"].to(z.device).to(torch.bool)
    cdr3_ag_beta = float(feats["cdr3_antigen_beta_value"].item())
    if abs(cdr3_ag_beta) > 1e-6:
        # CDR3_i paired with antigen_j, and antigen_i paired with CDR3_j
        interface_mask = (cdr3_mask_i & ag_mask.unsqueeze(-2)) | (ag_mask.unsqueeze(-1) & cdr3_mask_j)
        interface_pair_mask_4d = interface_mask.unsqueeze(0).unsqueeze(-1).float()
        token_trans_bias = token_trans_bias + cdr3_ag_beta * interface_pair_mask_4d * token_trans_bias
```

**Complexity**: Low–Medium. Featurizer needs antigen chain detection.

---

### D+.2 — Per-CDR Beta Values (H3 Stronger, H1/H2 Lighter)

**Goal**: Reflect biological reality where CDR-H3 dominates binding — apply stronger beta-scaling to H3 and lighter to other CDRs.

**Current**: Single `cdr3_token_mask` covers all CDR loops equally.

**Files to modify**:
- `src/boltz/data/feature/featurizerv2.py`: Produce separate `h3_token_mask`, `l3_token_mask`, `other_cdr_token_mask`
- `src/boltz/model/modules/diffusion_conditioning.py`: Apply three separate beta values with three masks
- `src/boltz/main.py`: Add `--h3_beta_value`, `--h1h2_beta_value`, `--l3_beta_value` flags

**YAML format** (new):
```yaml
steering:
  cdr_beta_scaling:
    h3_beta: 0.5
    l3_beta: 0.2
    h1h2_beta: 0.1
```

**Changes in `diffusion_conditioning.py`**:
```python
# Apply each CDR mask with its own beta, accumulate into token_trans_bias
for mask_key, beta_key in [("h3_token_mask", "h3_beta"), ("l3_token_mask", "l3_beta"), ...]:
    if mask_key in feats:
        mask = feats[mask_key].to(torch.bool)
        beta = float(feats[beta_key].item())
        pair_mask = (mask.unsqueeze(-1) & mask.unsqueeze(-2)).unsqueeze(0).unsqueeze(-1).float()
        token_trans_bias = token_trans_bias + beta * pair_mask * token_trans_bias
```

**Complexity**: Medium. Featurizer needs CDR loop boundary detection per chain.

---

### D+.3 — Diversity-Aware Model Selection

**Goal**: Instead of selecting the highest-confidence prediction, select a diverse set that covers different binding modes.

**Files**: New `src/boltz/analysis/diversity_selection.py`, called from evaluate.py.

```python
def diversity_aware_selection(predictions, quality_floor=0.5, k=5):
    """
    Select k predictions that are: (1) above quality_floor by iPTM,
    (2) maximally diverse by pairwise CDR RMSD.
    """
    qualified = [p for p in predictions if p.iptm >= quality_floor]
    if len(qualified) <= k:
        return qualified
    # Greedy max-min diversity selection
    selected = [qualified[0]]  # Start with highest confidence
    while len(selected) < k:
        remaining = [p for p in qualified if p not in selected]
        # Pick candidate maximally distant from all selected
        best = max(remaining, key=lambda p: min(rmsd(p, s) for s in selected))
        selected.append(best)
    return selected
```

**Complexity**: Low. Post-processing only.

---

## Strategy E+ — Blind Scanning v2

**Current**: Region-specific beta emphasis/de-emphasis on antigen surface patches. FAILED — DockQ drops 51%, AUC-ROC=0.645. Root cause: perturbation too aggressive (emphasis=0.5, de-emphasis=-0.3) causing structural collapse.

The N10 (E-Refined) in the new strategies addresses this more comprehensively. This section covers the incremental E+ upgrades from evaluation next steps.

---

### E+.1 — Softer Beta Magnitudes

**Goal**: Prevent the structural collapse that made E harmful. Use emphasis=0.2, de-emphasis=-0.1.

**Files to modify**: YAML config / featurizer parameter only.

**Note**: Section 5 mentions "Increase beta magnitudes (emphasis 1.0-2.0)" but E's failure was caused by aggressive perturbation causing structural collapse. The N10 proposal (softer 0.2/-0.1) is more consistent with the observed failure mode. The correct upgrade is **reduction**, not increase — the higher magnitudes would make sense only if the scanning signal is near zero and needs amplification while keeping DockQ intact, which contradicts the evidence. Recommend testing 0.2 (not 1.0).

**Complexity**: None — parameter only.

---

### E+.2 — More Scanning Regions (20–40 instead of 10)

**Goal**: Higher spatial resolution in the contact heatmap.

**Files to modify**: `src/boltz/data/feature/featurizerv2.py` (region generation) and the run script.

**New CLI flag**: `--scan_num_regions INT` (default: 10, try 20–30)

**Impact**: Linear increase in compute (each region is one forward pass). At 30 regions this is 3× more expensive than current E.

**Complexity**: Low — parameter only if region grid is already parameterized.

---

### E+.3 — Tighter Contact Threshold (5–6 Å instead of 8 Å)

**Goal**: More precise epitope localization, reducing false positives at the periphery.

**Files to modify**: `src/boltz/model/potentials/potentials.py` or wherever the 8 Å threshold is defined in the contact accumulation logic.

**New CLI flag**: `--contact_threshold_angstrom FLOAT` (default: 8.0, try 5.0–6.0)

**Complexity**: Low.

---

### E+.4 — Confidence-Weighted Contact Aggregation

**Goal**: Give more weight to high-confidence predictions when building the contact heatmap.

**Files to modify**: evaluate.py contact accumulation logic.

```python
# Current (unweighted):
contact_heatmap += contacts_from_prediction

# Proposed (confidence-weighted):
contact_heatmap += iptm_score * contacts_from_prediction
contact_weight_sum += iptm_score
contact_heatmap /= contact_weight_sum  # normalize
```

**Complexity**: Low.

---

## Strategy G+ — Progressive Steering v2

**Current**: Three fixed phases using `PiecewiseStepFunction` for guidance weights. CDR3 beta-scaling (constant beta throughout). Helps medium-difficulty cases but regresses on easy cases (−0.025 DockQ).

---

### G+.1 — Conditional Steering (Skip When Confidence High)

**Goal**: Eliminate the regression on easy cases by skipping steering when the model is already confident.

**Mechanism**: Run a fast B1 prediction (no steering, 1 sample) first. Measure iPTM. Gate steering based on result.

**Files to modify**: `src/boltz/main.py` (pre-check logic before main prediction)

```python
# Pseudocode in predict():
if steering_args.conditional_steering:
    # Quick baseline run: 1 sample, no steering, reduced steps
    baseline_iptm = run_baseline_prediction(feats, num_samples=1, reduced_steps=True)

    if baseline_iptm > steering_args.conditional_threshold_high:  # e.g., 0.85
        # Model already confident — skip steering
        steering_args.fk_steering = False
        steering_args.contact_guidance_update = False
    elif baseline_iptm > steering_args.conditional_threshold_low:  # e.g., 0.75
        # Mild case: use Y-style beta-only (no FK)
        steering_args.fk_steering = False
    # else: full steering as planned
```

**New CLI flags**:
```
--conditional_steering         (bool)
--conditional_threshold_high   FLOAT  (default: 0.85)
--conditional_threshold_low    FLOAT  (default: 0.75)
--conditional_baseline_steps   INT    (default: 50, reduced step count for baseline)
```

**Complexity**: Low–Medium. The pre-check adds one forward pass per complex.

---

### G+.2 — Time-Varying Beta (Stronger in Phase 1)

**Goal**: Apply stronger CDR3 beta-scaling in early diffusion steps (Phase 1 = exploration) and decay to zero in late steps (Phase 3 = convergence).

**Problem**: The current `DiffusionConditioning` receives a fixed `cdr3_beta_value` from `feats`, which is constant throughout the run. It has no access to the current diffusion timestep (`steering_t`).

**Key change**: Pass `steering_t` from the diffusion loop into `DiffusionConditioning` so beta can vary.

**Files to modify**:
- `src/boltz/model/modules/diffusionv2.py`: Inject `steering_t` into `network_condition_kwargs` at each step
- `src/boltz/model/modules/diffusion_conditioning.py`: Read `steering_t` from feats and compute `beta(t) = beta_max × steering_t`

**Changes in `diffusionv2.py`** (inside the denoising loop at step_idx):
```python
# After computing steering_t (line 375):
steering_t = 1.0 - (step_idx / num_sampling_steps)

# Add to feats for this step's forward pass:
network_condition_kwargs["feats"]["current_steering_t"] = torch.tensor(steering_t)
```

**Changes in `diffusion_conditioning.py`**:
```python
# Instead of fixed beta_val:
if "current_steering_t" in feats and "cdr3_beta_max" in feats:
    t = float(feats["current_steering_t"].item())
    beta_max = float(feats["cdr3_beta_max"].item())
    beta_val = beta_max * t  # Decays linearly: high at start, 0 at end
```

**New CLI flags**:
```
--cdr3_beta_max FLOAT         (max beta at start of diffusion, e.g. 0.6)
--cdr3_beta_schedule STR      (choices: "linear", "cosine", "step"; default: "linear")
```

**Complexity**: Medium. Requires threading `steering_t` through the conditioning forward pass, which touches the inner diffusion loop.

---

### G+.3 — Interface-Specific Model Selection

**Goal**: Select predictions based on CDR-antigen interface quality rather than global iPTM.

**Files**: evaluate.py, new `src/boltz/analysis/interface_metrics.py`

```python
def compute_interface_plddt(prediction, cdr_residues, antigen_residues, distance_cutoff=8.0):
    """
    Average pLDDT of residues within distance_cutoff of the interface.
    """
    coords = prediction.atom_coords
    plddt = prediction.plddt
    # Find CDR atoms within cutoff of any antigen atom
    dists = pairwise_distances(coords[cdr_residues], coords[antigen_residues])
    interface_cdr_mask = (dists.min(dim=1).values < distance_cutoff)
    return plddt[cdr_residues][interface_cdr_mask].mean()
```

**New CLI flag**: `--rerank_by_interface_plddt` (bool)

**Complexity**: Low.

---

## Strategy K+ — Region Beta-Scaling v2

**Current**: Single region emphasis configuration per run. Results: redistribution pattern (big wins AND big losses), no net significance.

---

### K+.1 — Ensemble Over Multiple Region-Emphasis Patterns

**Goal**: Run 5–8 different region configurations, aggregate contact maps, select best structure from the union.

**Files to modify**: Run script / orchestration logic in `src/boltz/main.py`

**Approach**: Add `--k_num_configs INT` flag. The main prediction loop iterates over `k_num_configs` different beta patterns (each emphasizing a different antigen patch), collects all predictions, and passes them to the re-ranker.

```python
# Pseudocode in predict() or a wrapper:
all_predictions = []
for config_id in range(steering_args.k_num_configs):
    # Generate region emphasis for this config (round-robin over antigen patches)
    feats["region_emphasis_id"] = torch.tensor(config_id)
    preds = model.sample(..., feats=feats)
    all_predictions.extend(preds)
# Select best from union
best = select_by_composite_score(all_predictions)
```

**New CLI flags**:
```
--k_num_configs INT      (default: 1, try 5-8)
--k_emphasis_beta FLOAT  (default: 0.5, try 0.2-0.3 per D+.2 lesson)
```

**Complexity**: Low. Outer loop only; no model changes.

---

### K+.2 — Softer Beta Values

**Goal**: Reduce the big losses seen in K by using less aggressive emphasis.

**Evidence**: K at beta=0.5 caused large swings. Beta=0.2–0.3 should reduce variance.

**Files**: YAML config or `--k_emphasis_beta` CLI flag (no model changes).

**Complexity**: None — parameter only.

---

### K+.3 — Contact Heatmap Accumulation Pipeline

**Goal**: Build a proper per-antigen-residue contact frequency map across all K configs.

**Files**: New `src/boltz/analysis/contact_heatmap.py`

```python
class ContactHeatmap:
    def __init__(self, num_antigen_residues):
        self.heatmap = np.zeros(num_antigen_residues)
        self.weight_sum = 0.0

    def add_prediction(self, pred, cdr_residues, antigen_residues, weight=1.0,
                       contact_threshold_ang=6.0):
        dists = compute_min_distances(pred.coords[cdr_residues], pred.coords[antigen_residues])
        contacts = (dists < contact_threshold_ang).float()  # [num_antigen_residues]
        self.heatmap += weight * contacts.numpy()
        self.weight_sum += weight

    def get_top_residues(self, top_k=10):
        normalized = self.heatmap / self.weight_sum
        return np.argsort(normalized)[-top_k:][::-1]
```

**Complexity**: Low. Analysis-side only.

---

### K+.4 — Two-Stage Pipeline: K Epitope → B2-Style Restraints

**Goal**: Use K's contact heatmap to automatically generate contact restraints (like B2's oracle contacts), eliminating the oracle requirement.

**Files**: New `src/boltz/analysis/epitope_to_constraints.py`

```python
def heatmap_to_yaml_constraints(heatmap, cdr_residues, antigen_chain, cdr_chain,
                                 top_k=5, distance_cutoff=8.0):
    """
    Convert top-k heatmap residues into Boltz2 YAML contact constraints.
    """
    top_antigen_residues = heatmap.get_top_residues(top_k)
    constraints = []
    for ag_res in top_antigen_residues:
        for cdr_res in cdr_residues:
            constraints.append({
                "contact": {
                    "token1": {"chain": cdr_chain, "res": int(cdr_res)},
                    "token2": {"chain": antigen_chain, "res": int(ag_res)},
                    "distance_threshold": distance_cutoff,
                    "force": False,  # Soft constraint
                }
            })
    return constraints
```

**Stage 2 run**: Take the generated YAML constraints → pass to B2-style `ContactPotential`.

**Complexity**: Low–Medium. Requires the contact heatmap pipeline (K+.3) and a YAML generation step between Stage 1 and Stage 2 predictions.

---

## Strategy L+ — CDR3 Beta-Scaling v2

**Current**: beta=0.3 applied uniformly to CDR3-CDR3 token pairs. Diversity INCREASED successfully (1.51→1.66 Å, 45%→60%), but CDR3 RMSD worsened (p=0.0008) and iPTM was inflated.

---

### L+.1 — Reduce Beta (0.1–0.15 instead of 0.3)

**Goal**: Keep the diversity benefit while reducing the accuracy penalty.

**Hypothesis**: beta=0.3 is too strong and distorts CDR3 geometry. At 0.1–0.15 the diversity effect may persist at smaller amplitude.

**Files**: No code change needed — just use `--cdr3_beta_value 0.1` or `0.12` in the run script.

**Complexity**: None — parameter only. Should be the first experiment.

---

### L+.2 — Asymmetric CDR Scaling (H3 Stronger, L3 Weaker)

**Goal**: Focus the diversity push on CDR-H3 (the dominant loop) and reduce impact on CDR-L3 (where diversity worsened accuracy most).

**This upgrade shares implementation with D+.2** (per-CDR beta values).

**Files to modify**: Same as D+.2:
- `src/boltz/data/feature/featurizerv2.py`: Separate H3/L3 masks
- `src/boltz/model/modules/diffusion_conditioning.py`: Apply per-CDR beta

**Suggested values**: H3 beta=0.15, L3 beta=0.05 (or 0.0 initially)

**Complexity**: Medium (shared with D+.2).

---

### L+.3 — Expand Mask to CDR3-Antigen Interface Pairs

**Goal**: Instead of only scaling CDR3-CDR3 attention, also scale CDR3-antigen attention to push the loop toward the antigen rather than just exploring internally.

**This upgrade shares implementation with D+.1** (CDR-antigen pair beta).

**Key difference from D+.1**: In L+ context, the antigen mask direction should **attract** (positive beta on CDR3-antigen pairs), while CDR3-CDR3 scaling handles internal loop geometry.

**Suggested config**: `cdr3_cdr3_beta=0.1` (internal), `cdr3_antigen_beta=0.15` (interface)

**Complexity**: Medium (shared with D+.1).

---

### L+.4 — Time-Dependent Beta Schedule

**Goal**: Higher beta in early diffusion (exploration) decaying to zero late (convergence). Prevents L from distorting final-step geometry.

**This upgrade shares implementation with G+.2** (time-varying beta via `current_steering_t`).

**Suggested schedule**: `beta(t) = beta_max × t` (linear decay from `beta_max` at start to 0 at end)

**Files**: Same as G+.2:
- `src/boltz/model/modules/diffusionv2.py`: Inject `current_steering_t`
- `src/boltz/model/modules/diffusion_conditioning.py`: Compute `beta(t)` from `beta_max` and `t`

**Complexity**: Medium (shared with G+.2 — implement once for both).

---

## Strategy Q+ — Iterative Epitope Refinement v2

**Current**: 3-round iterative loop where contacts from round N inform steering in round N+1. Results: best epitope improvement (F1 0.233→0.373), sporadic DockQ wins, p=0.753 (not significant).

---

### Q+.1 — More Round 1 Diversity (10–20 samples, slightly higher noise)

**Goal**: Get a wider initial contact heatmap so that the Round 2 hotspot target is more robust.

**Current**: Round 1 uses a fixed number of samples (likely ~5) at standard noise scale.

**Files to modify**: Q iteration run script, `src/boltz/main.py`

**New CLI flags**:
```
--q_round1_samples INT           (default: 5, try 10-20)
--q_round1_noise_scale FLOAT     (default: 1.0, try 1.1-1.2 for more diversity)
```

**Implementation**: The `noise_scale` parameter already exists in the diffusion schedule. Passing a higher value for Round 1 increases `eps` magnitude in the noising step (line 377 of `diffusionv2.py`: `eps = sqrt(noise_var) * torch.randn(...)` — this is controlled by `self.noise_scale`).

**Complexity**: Low.

---

### Q+.2 — Confidence-Weighted Contact Maps

**Goal**: Downweight low-confidence predictions when building the contact heatmap so that bad-pose predictions don't pollute the epitope signal.

**Files**: Q iteration analysis code / evaluate.py

```python
# In contact heatmap accumulation:
for pred in round1_predictions:
    weight = pred.iptm  # or pred.interface_plddt
    heatmap.add_prediction(pred, cdr_residues, antigen_residues, weight=weight)
```

**Complexity**: Low.

---

### Q+.3 — Ensemble-Based Hotspot Selection

**Goal**: Select hotspot residues more robustly using consensus across multiple confidence thresholds, not just a single cutoff.

**Files**: New `src/boltz/analysis/hotspot_selection.py`

```python
def robust_hotspot_selection(heatmap, thresholds=[0.1, 0.2, 0.3], min_consensus=2):
    """
    A residue is a hotspot if it exceeds threshold in at least min_consensus
    threshold levels. Reduces sensitivity to single-threshold choice.
    """
    votes = np.zeros(len(heatmap))
    for thresh in thresholds:
        votes += (heatmap > thresh).astype(int)
    return np.where(votes >= min_consensus)[0]
```

**Complexity**: Low.

---

### Q+.4 — Adaptive Rounds (Skip if Hotspot Already Clear)

**Goal**: Avoid unnecessary round 2/3 computation when round 1 already gives a strong epitope signal.

**Logic**: After Round 1, compute contact entropy of the heatmap. If entropy is below a threshold (concentrated signal), skip Round 2 and go directly to Round 3 refinement.

**Files**: Q iteration run script

```python
def contact_entropy(heatmap):
    """Shannon entropy of normalized contact map. Low = concentrated signal."""
    p = heatmap / heatmap.sum()
    p = p[p > 0]
    return -np.sum(p * np.log(p))

# In Q main loop:
if contact_entropy(round1_heatmap) < steering_args.q_entropy_threshold:
    log("Strong hotspot detected in Round 1, skipping Round 2")
    skip_round2 = True
```

**New CLI flag**: `--q_entropy_threshold FLOAT` (default: 2.0, lower = more permissive skip)

**Complexity**: Low.

---

## Strategy V+ — Embedding CDR3 Steering v2

**Current**: `EmbeddingCDR3Potential` (or equivalent) performs gradient descent in embedding space targeting the framework-region mean of CDR3 embeddings. HARMFUL: DockQ −0.034 (p=0.0004), CDR3 RMSD +1.28 Å. Wrong target.

---

### V+.1 — Alternative Target: Canonical Class Embeddings

**Goal**: Replace the framework mean (meaningless average) with embeddings extracted from PDB structures of known CDR3 canonical conformations.

**Implementation strategy**:

1. **Offline (pre-computation)**:
   - Collect PDB structures for each CDR3 canonical class (H3 torso, kinked, extended, etc.)
   - Run Boltz2 trunk on each → extract token embeddings at CDR3 positions
   - Average per canonical class → save as lookup table `canonical_embeddings.pt`

2. **Online (during prediction)**:
   - Detect predicted CDR3 length + sequence features → predict likely canonical class
   - Load corresponding canonical embedding as target
   - Apply gradient: push CDR3 embeddings toward canonical target

**Files to modify**:
- New `src/boltz/data/canonical_embeddings/` directory with pre-computed embeddings
- `src/boltz/model/potentials/potentials.py`: `EmbeddingCDR3Potential.compute()` loads canonical target instead of computing framework mean
- `src/boltz/main.py`: `--v_canonical_class STR` or `--v_canonical_db_path PATH`

**Data flow**:
```
CDR3 sequence → length + hydrophobicity → predicted canonical class
→ load canonical_embeddings[class] → gradient target
```

**Complexity**: Medium–High. Requires offline pre-computation pass and a CDR3 class predictor.

---

### V+.2 — Reduced Gradient Strength (0.01–0.05 instead of 1.0)

**Goal**: Prevent the catastrophic distortion (CDR3 RMSD +1.28 Å) caused by overly strong embedding gradient.

**Files**: `src/boltz/model/potentials/potentials.py` — change `guidance_weight` in `get_potentials()` for the embedding potential.

```python
# Current (problematic):
"guidance_weight": 1.0

# Proposed:
"guidance_weight": PiecewiseStepFunction(
    thresholds=[0.7], values=[0.02, 0.0]  # Only in first 30% of diffusion
)
```

**Complexity**: Low — parameter change.

---

### V+.3 — Apply Only in First 30% of Diffusion Steps

**Goal**: Allow recovery from any distortion by disabling embedding gradient after the early exploration phase.

**Implementation**: Use `PiecewiseStepFunction` for `guidance_weight` and `resampling_weight` that returns 0 when `steering_t < 0.7`.

```python
# In get_potentials():
"guidance_weight": PiecewiseStepFunction(
    thresholds=[0.7],
    values=[v_strength, 0.0]  # Active only when steering_t > 0.7 (first 30%)
)
```

**Note**: `steering_t = 1.0` at step 0, decreasing to 0. So `steering_t > 0.7` corresponds to the first 30% of steps.

**Complexity**: Low — uses existing scheduling primitives.

---

### V+.4 — CDR3-to-Antigen Pair Target

**Goal**: Instead of targeting CDR3 self-embeddings, target CDR3-antigen pair embeddings from known complexes, pushing CDR3 into an "interface-ready" conformation.

**Approach**:
1. Offline: from B2's successful predictions, extract `z[cdr3_i, antigen_j]` pair embeddings
2. Average across successful complexes → target interface embedding profile
3. Online: compute loss = distance of current CDR3-antigen pair embeddings from this target

**Files**:
- New offline analysis script to extract pair embeddings from B2 runs
- `src/boltz/model/potentials/potentials.py`: New `CDR3AntigenPairEmbeddingPotential`
- `src/boltz/model/modules/diffusionv2.py`: Thread `z_trunk` into this potential (already done for W)

**Complexity**: High. Requires offline computation from B2's runs + new potential class.

---

## Strategy W+ — Embedding Interface Steering v2

**Current**: `EmbeddingInterfacePotential` weights distance forces by pair embedding norms `‖z[i,j]‖`. Results: negligible +0.004 DockQ. Embedding norms don't encode binding specificity — they encode token identity and structural context, not interface propensity.

---

### W+.1 — Validate Embedding Signal First

**Goal**: Before investing in a learned head, confirm that any signal exists in pair embeddings for binding prediction.

**This is an analysis step, not a code change.**

**New script**: `src/boltz/analysis/validate_embedding_signal.py`

```python
def validate_embedding_auc(b2_predictions, true_contact_matrix):
    """
    For each CDR-antigen pair (i,j):
      - Feature: norm(z[i,j]) from the pair embedding
      - Label: 1 if true contact (< 8A), 0 otherwise
    Compute AUC-ROC across all pairs across all complexes.
    Target: AUC-ROC > 0.65 to justify learned head.
    """
    features, labels = [], []
    for pred, contacts in zip(b2_predictions, true_contact_matrix):
        for i in cdr_indices:
            for j in antigen_indices:
                features.append(norm(pred.z_trunk[i, j]))
                labels.append(contacts[i, j])
    return roc_auc_score(labels, features)
```

**Decision gate**: If AUC-ROC < 0.60, the embedding norms have no useful signal and W is fundamentally broken. If 0.60–0.70, a learned head may help. If > 0.70, proceed with W+.2.

**Complexity**: Low — analysis only.

---

### W+.2 — Learned Scoring Head Instead of Raw Norms

**Goal**: Replace `weight = norm(z[i,j])` with `weight = MLP(z[i,j])` trained on B2's true contact labels.

**Files to modify**:
- New `src/boltz/model/modules/interface_scoring_head.py`
- `src/boltz/model/potentials/potentials.py`: `EmbeddingInterfacePotential` loads and calls the scoring head

**Architecture**:
```python
class InterfaceScoringHead(nn.Module):
    def __init__(self, token_z=128):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(token_z, 64),
            nn.ReLU(),
            nn.Linear(64, 16),
            nn.ReLU(),
            nn.Linear(16, 1),
            nn.Sigmoid()  # Output: contact probability in [0,1]
        )
    def forward(self, z_ij):  # z_ij: [batch, n_cdr, n_antigen, token_z]
        return self.net(z_ij).squeeze(-1)  # [batch, n_cdr, n_antigen]
```

**Training**: ~47 complexes × multiple seeds × all CDR-antigen pairs ≈ tens of thousands of samples. May be too few for a reliable head. Consider pre-training on all B2 runs across the full dataset.

**New CLI flag**: `--interface_scoring_head_path PATH`

**Complexity**: High. Requires training data from B2 runs.

---

### W+.3 — Re-Extract Embeddings During Diffusion

**Goal**: Prevent the embedding weights from becoming stale as the diffusion trajectory evolves.

**Current**: `z_trunk` is extracted once before the diffusion loop and reused for all 200 steps.

**Problem**: CDR3 embeddings in `z_trunk` are fixed at `t=T` (maximum noise), while coordinates evolve toward structure. By step 100, the pair embeddings no longer reflect the current trajectory.

**Files to modify**: `src/boltz/model/modules/diffusionv2.py`, `src/boltz/model/models/boltz2.py`

**Approach**:
1. In `boltz2.py`, expose the trunk computation as a separate method: `def compute_trunk(feats) -> (s, z)`
2. In `diffusionv2.py`, at every `w_reextract_interval` steps (e.g., every 10 steps), call `compute_trunk()` with current `atom_coords_denoised` injected as pseudo-coordinates, re-extract `z`, call `potential.set_embedding_weights(z)`

**New CLI flag**: `--w_reextract_interval INT` (default: ∞ = never, try 10–20)

**Complexity**: High. Requires a separate trunk forward pass per reextraction interval — adds significant compute.

---

### W+.4 — Use Embeddings to Generate B2-Style Restraints

**Goal**: Instead of using embedding norms as continuous potential weights, use them to identify and generate discrete contact restraints.

**Approach**:
1. After B1 prediction, extract `z_trunk` from the clean structure
2. Score all CDR-antigen pairs with `InterfaceScoringHead` (or norm threshold)
3. Top-k pairs → convert to YAML contact constraints (same format as B2)
4. Re-run with those constraints active

**This is a pipeline approach** (not a steering approach) that converts W's embedding signal into B2-style restraints.

**Files**: New `src/boltz/analysis/embeddings_to_constraints.py`

**Complexity**: Medium. Requires W+.1 validation and W+.2 scoring head.

---

## Strategy Y+ — Hierarchical Steering v2

**Current**: Early phase = CDR3 beta-scaling in `DiffusionConditioning`. Late phase = `AntigenOrientationPotential` in coordinate space. Results: Only significant result in Round 1 — CDR-H3 RMSD −0.30 Å (p=0.025), +0.026 DockQ trend.

Y is the most promising strategy. The upgrades focus on amplifying both phases and adding epitope information.

---

### Y+.1 — Stronger Beta-Scaling in Early Phase

**Goal**: Amplify the embedding-space steering that produced the significant CDR-H3 improvement.

**Current beta value**: From the evaluation context, the beta used in Y was likely 0.3–0.5 (constant). Results suggest it was effective but modest.

**Proposed**: Use time-varying beta (from G+.2 infrastructure) with `beta_max=0.8–1.0` decaying linearly. The early phase should strongly orient CDR3 before coordinate-space potentials take over.

**Implementation**: Identical to G+.2 — this is the same code change:
- `diffusionv2.py`: inject `current_steering_t` into kwargs
- `diffusion_conditioning.py`: compute `beta(t) = beta_max × t`

**Suggested values**: `beta_max=0.8`, linear decay

**Complexity**: Medium (shared with G+.2 — implement once).

---

### Y+.2 — Enhanced Late-Phase Coordinate Potentials

**Goal**: Improve the quality of the AntigenOrientationPotential in the late phase by also adding CDR3 conformation steering.

**Current late phase**: Only `AntigenOrientationPotential` (global antigen orientation).

**Proposed additions**:
1. `CDR3ConformationPotential` in late phase (already implemented, just needs scheduling)
2. A new `CDR3AntigenProximityPotential` that directly penalizes large CDR-antigen distances

**Files**: `src/boltz/model/potentials/potentials.py`, `get_potentials()` factory

```python
# In get_potentials(), add for Y+ mode:
CDR3AntigenProximityPotential(
    parameters={
        "guidance_interval": 2,
        "guidance_weight": PiecewiseStepFunction(
            thresholds=[0.5], values=[0.0, 0.8]  # Only active in late phase (steering_t < 0.5)
        ),
        "resampling_weight": PiecewiseStepFunction(
            thresholds=[0.5], values=[0.0, 0.5]
        ),
        "target_distance_ang": 10.0,  # CDR3 Cα to nearest antigen Cα
    }
)
```

**New potential class** `CDR3AntigenProximityPotential`:
- Energy: soft hinge loss on minimum distance between CDR3 atoms and antigen atoms
- Similar to `ContactPotential` but uses CDR3 mask + antigen mask instead of constraint pairs

**Complexity**: Medium. New potential class following existing pattern.

---

### Y+.3 — Incorporate Predicted Epitope Information

**Goal**: Give the late-phase potentials a specific target region on the antigen, rather than applying isotropic orientation guidance.

**Source of epitope**: Either Q's output from a prior run, or a sequence-based epitope predictor.

**New YAML field**:
```yaml
steering:
  hierarchical:
    early_beta_max: 0.8
    predicted_epitope_residues:
      chain: A
      residues: [45, 46, 47, 48, 67, 68, 102]  # From Q's output
```

**Files to modify**:
- `src/boltz/data/parse/schema.py`: Parse `predicted_epitope_residues`
- `src/boltz/data/feature/featurizerv2.py`: Produce `epitope_token_mask` from parsed residues
- `src/boltz/model/potentials/potentials.py`: `AntigenOrientationPotential` (or new `EpitopeFocusPotential`) reads `epitope_token_mask` from feats and focuses guidance on those residues
- `src/boltz/main.py`: `--predicted_epitope_residues PATH` (path to JSON/YAML with residue list)

**Potential logic change**: Instead of penalizing all CDR-antigen distances equally, focus the potential energy computation on CDR vs. epitope residue distances.

**Complexity**: Medium. Requires schema/featurizer extension (same pattern as existing constraint parsing).

---

### Y+.4 — Orientation Pre-Check for the Hard 40%

**Goal**: Identify complexes where the model places the antibody on the wrong hemisphere (~40% of failures) and apply stronger forcing.

**Diagnostic**: After running B1 prediction, compute CDR-antigen centroid distance. If > 30 Å, the global orientation is likely wrong.

**Files**: `src/boltz/main.py` (pre-check), `src/boltz/model/potentials/potentials.py` (stronger orientation forcing option)

```python
# In predict() after B1 pre-check (extends G+.1 infrastructure):
cdr_centroid = coords[cdr_mask].mean(dim=0)
antigen_centroid = coords[antigen_mask].mean(dim=0)
cdr_ag_distance = (cdr_centroid - antigen_centroid).norm()

if cdr_ag_distance > steering_args.hard_case_distance_threshold:  # e.g., 30.0 A
    # Hard case: increase orientation forcing strength
    steering_args.y_antigen_guidance_weight_scale = 2.0  # 2× stronger
    log(f"Hard case detected (CDR-Ag dist={cdr_ag_distance:.1f}A), increasing Y guidance")
```

**In `get_potentials()`**: Scale `AntigenOrientationPotential` guidance weights by `antigen_guidance_weight_scale` from `steering_args`.

**New CLI flags**:
```
--hard_case_distance_threshold FLOAT  (default: 30.0 Angstrom)
--y_antigen_guidance_weight_scale FLOAT  (default: 1.0, auto-set if hard case)
```

**Complexity**: Low. Reuses the G+.1 pre-check infrastructure.

---

## Cross-Cutting Changes

### Shared Infrastructure (Implement Once)

Several upgrades above share the same underlying code change. Implement each once:

| Infrastructure | Strategies That Need It | Key Change |
|---|---|---|
| Time-varying beta via `current_steering_t` | G+.2, L+.4, Y+.1 | Thread `steering_t` from loop → feats → `DiffusionConditioning` |
| Per-CDR masks (H3/L3 separate) | D+.2, L+.2 | Featurizer: split `cdr3_token_mask` into per-CDR masks |
| CDR-antigen pair mask | D+.1, L+.3 | Featurizer: add `antigen_token_mask`; conditioning: apply interface beta |
| Confidence-based pre-check | G+.1, Y+.4 | `main.py`: run B1 first, measure iPTM + geometry |
| Composite re-ranker | A+.2, K+.1 | `analysis/reranking.py`: weighted score by energy + iptm + interface_plddt |
| Contact heatmap accumulator | E+.4, K+.3, Q+.2 | `analysis/contact_heatmap.py`: confidence-weighted contact map |

---

### Implementation Priority (Strategy-Level)

Each strategy is implemented as a **single round** — all sub-improvements (e.g. A+.1 through A+.4) are built together in one worktree before moving to the next strategy. The table below ranks complete strategy bundles by expected impact, value, and combined complexity, considering dependencies between strategies.

#### Dependency Graph

```
Q+ ──────────────────────────────────────────┐
A+ (self-contained)                          │
G+ ──┬── creates shared infra ──────────┐    │
     │   (time-varying beta,            │    │
     │    confidence pre-check)         │    │
     │                                  ▼    ▼
     └──────────────────────────────── Y+ (uses G+ infra + Q+ output)
                                        │
L+ (uses G+.2 time-varying beta,        │
    implements per-CDR & CDR-Ag masks) ──┘
K+ (self-contained, shares heatmap logic with Q+)
D+ (reuses L+'s mask infrastructure)
E+ (self-contained, parameter-level)
V+ (research, gated by canonical embedding precomputation)
W+ (research, gated by W+.1 signal validation)
```

#### Priority Table

| Priority | Strategy | Upgrades Included | Round 1 Rank | Expected Impact | Bundle Complexity | Rationale |
|---|---|---|---|---|---|---|
| 1 | **Q+** | Q+.1 more R1 diversity, Q+.2 confidence-weighted contacts, Q+.3 ensemble hotspot selection, Q+.4 adaptive rounds | #5 (best epitope) | High | Low | All upgrades are script/analysis-level changes. Strengthens the best blind epitope finder. **Output feeds Y+.3** later, so must come first. Self-contained — no model code changes needed. |
| 2 | **A+** | A+.1 OOM auto-reduction, A+.2 energy-based re-ranking, A+.3 seeds vs particles, A+.4 embedding fallback | #1 (best DockQ) | Very High | Low–Medium | Recovers the 12/47 missing complexes (25% of data). A was already the best strategy; upgrades are orchestration and one medium change (energy return threading). Creates the composite re-ranker (`reranking.py`) reusable by K+ and Y+. Self-contained. |
| 3 | **G+** | G+.1 confidence-based pre-check, G+.2 time-varying beta, G+.3 interface model selection | #3 (+0.017 DockQ) | High | Medium | **Creates critical shared infrastructure**: (1) time-varying beta via `steering_t` injection — needed by L+.4 and Y+.1, (2) confidence pre-check — reused by Y+.4. Fixes the easy-case regression (-0.025 DockQ) that limits G and Y. Must come before L+ and Y+. |
| 4 | **L+** | L+.1 reduce beta to 0.1–0.15, L+.2 asymmetric H3/L3 scaling, L+.3 CDR3-antigen pair mask, L+.4 time-dependent beta schedule | #4 (diversity success) | High | Low–Medium | L+.1 is a zero-code parameter change (run first as sanity check). L+.4 reuses G+.2's time-varying beta infra. L+.2 and L+.3 **implement per-CDR masks and CDR-antigen masks** — shared infrastructure also used by D+ later. Fixes the only strategy that achieved its diversity goal while eliminating its accuracy penalty. |
| 5 | **Y+** | Y+.1 stronger early beta, Y+.2 CDR3-antigen proximity potential, Y+.3 predicted epitope input, Y+.4 orientation pre-check | #2 (only significant) | Very High | Medium–High | Most promising strategy (only significant structural result in Round 1). Depends on G+.2 (time-varying beta) and G+.1 (pre-check), both already built. Y+.3 uses Q+'s epitope output for targeted late-phase guidance. Y+.2 requires a new potential class (`CDR3AntigenProximityPotential`). Placed here because dependencies must exist first. |
| 6 | **K+** | K+.1 multi-config ensemble, K+.2 softer beta values, K+.3 heatmap accumulation pipeline, K+.4 two-stage K→B2 pipeline | #6 (redistribution) | Medium | Low–Medium | Ensemble approach captures the union of K's sporadic wins. K+.3 heatmap pipeline shares logic with Q+.2 (already built). K+.4 converts heatmap into B2-style contact restraints — potentially closes the gap to B2 without oracle info. Self-contained. |
| 7 | **E+** | E+.1 softer beta (0.2), E+.2 more regions (20–30), E+.3 tighter contact threshold (5–6 Å), E+.4 confidence-weighted aggregation | #10 (failed) | Low–Medium | Low | All changes are parameter/config-level or analysis-side. Low investment to test whether E can be rescued with gentler parameters. E+.4 heatmap logic already exists from K+.3/Q+.2. Risk: E may be fundamentally broken regardless of parameters. |
| 8 | **D+** | D+.1 CDR-antigen pair beta, D+.2 per-CDR masks, D+.3 diversity-aware selection | #7 (negligible) | Low–Medium | Low | Most infrastructure already built by L+ (per-CDR masks, CDR-antigen masks). D+.3 diversity selection is a small post-processing addition. Low marginal effort, but D's base mechanism (canonical ensemble) showed minimal effect. |
| 9 | **V+** | V+.1 canonical class embeddings, V+.2 reduced gradient (0.01–0.05), V+.3 early-only (first 30%), V+.4 CDR3-antigen pair targets | #11 (harmful) | Low | High | V was the most harmful strategy (DockQ −0.034, p=0.0004). V+.1 requires offline precomputation of canonical CDR3 embeddings via Boltz2 trunk. V+.2/V+.3 are parameter fixes (low effort). V+.4 is a research direction. Worth pursuing only after higher-priority strategies are evaluated — if Y+.1 (embedding-space via beta-scaling) proves effective, V+ may be unnecessary. |
| 10 | **W+** | W+.1 validate embedding signal (AUC-ROC), W+.2 learned scoring head, W+.3 re-extract embeddings, W+.4 embedding→restraints pipeline | #8 (negligible) | Low | High | **Gated**: W+.1 (analysis-only) must confirm AUC-ROC > 0.60 before any further investment. If signal absent, abandon W entirely. W+.2 requires training data from B2 runs. W+.3 requires trunk re-computation during diffusion (expensive). Lowest priority — embedding norms may simply not encode binding specificity. |

#### Notes on Dependencies and Shared Infrastructure

- **Strategies 1–2 (Q+, A+)** are fully self-contained and can be implemented in parallel if desired.
- **Strategy 3 (G+)** is the infrastructure gateway — it creates the time-varying beta and confidence pre-check mechanisms that strategies 4 (L+) and 5 (Y+) depend on. Must be completed before them.
- **Strategy 4 (L+)** implements the per-CDR mask and CDR-antigen mask infrastructure. These are needed by D+ (#8), so L+ must precede D+.
- **Strategy 5 (Y+)** is the highest-impact strategy but sits at priority 5 because it depends on G+ (infra) and benefits from Q+ (epitope output for Y+.3).
- **Strategies 9–10 (V+, W+)** are research-gated — they should only proceed after confirming that the core embedding hypothesis holds (V+ needs canonical targets, W+ needs signal validation).

---

## Summary Table

| Strategy | Upgrades | Files Changed | Bundle Complexity | Key Risk |
|---|---|---|---|---|
| Q+ | More R1 diversity, confidence-weighted contacts, ensemble hotspots, adaptive rounds | Run script + `analysis/` | Low | Entropy threshold calibration |
| A+ | OOM fix, energy re-ranking, seeds vs particles, embedding fallback | `main.py`, `diffusionv2.py`, new `reranking.py` | Low–Medium | Energy return path threading |
| G+ | Conditional steering, time-varying beta, interface selection | `main.py`, `diffusionv2.py`, `diffusion_conditioning.py` | Medium | Baseline pre-check adds latency |
| L+ | Reduce beta, asymmetric CDR, interface pairs, time schedule | `featurizerv2.py`, `diffusion_conditioning.py`, `main.py` | Low–Medium | Per-CDR mask detection accuracy |
| Y+ | Stronger beta, proximity potential, epitope input, hard-case detect | `main.py`, `potentials.py`, `schema.py`, `featurizerv2.py` | Medium–High | New potential class + schema extension |
| K+ | Multi-config, softer beta, heatmap, two-stage pipeline | `main.py`, new `analysis/` scripts | Low–Medium | Stage 2 YAML generation |
| E+ | Softer beta, more regions, tighter threshold, confidence weighting | Config + `evaluate.py` | Low | E may be fundamentally broken |
| D+ | CDR-antigen beta, per-CDR masks, diversity selection | Reuses L+ infra + new `analysis/diversity_selection.py` | Low | Low independent impact |
| V+ | Canonical targets, reduced strength, early-only, pair targets | `potentials.py`, new canonical data | High | Offline embedding precomputation |
| W+ | Signal validation, learned head, re-extraction, restraint generation | New `analysis/` + `potentials.py` | High | Gated by W+.1 AUC-ROC result |

---

**Document Status**: Planning phase — no code changes made
**Last Updated**: 2026-03-21
