# Hierarchical Steering (Strategy Y+) — Intuitive Explanation

## Goal

### Main goal

Improve the quality of predicted antibody-antigen docked structures — specifically, get the antibody to bind the antigen in the correct orientation and with accurate CDR loop conformations, **without requiring any experimental binding data** (no oracle contacts, no known epitopes).

Vanilla Boltz-2 treats all residues equally during diffusion. It has no built-in understanding that CDR loops are the business end of the antibody, or that they need to physically reach the antigen surface. Hierarchical steering injects this domain knowledge into the diffusion process at the right moments.

### Why this matters

In antibody design and drug discovery, knowing the 3D structure of how an antibody binds its target antigen is critical. Even small improvements in docking accuracy (e.g., CDR-H3 RMSD reduced by 0.3 Å) can be clinically meaningful — it can mean the difference between correctly identifying a binding pose and missing it entirely.

### Sub-goals

1. **Improve CDR loop conformation accuracy** — The CDR loops (especially CDR-H3) are the most variable and hardest-to-predict parts of the antibody. By amplifying attention between CDR residue pairs during early diffusion steps, we guide the network to explore better loop conformations when the structure is still malleable.

2. **Improve global docking orientation** — Round 1 analysis showed that ~40% of complexes fail at the global level: the antibody is pointed at the wrong hemisphere of the antigen, and no amount of local refinement can fix that. Late-phase coordinate potentials explicitly pull the antibody's binding surface toward the antigen.

3. **Focus docking on the correct binding site** — When a predicted structure is available from a previous run, we can extract which antigen residues were near the CDR loops (the predicted epitope) and restrict the docking guidance to that patch. This prevents the antibody from being guided toward the wrong face of the antigen.

4. **Adapt guidance strength per complex** — Not all complexes are equally difficult. Easy cases need gentle guidance (or none), while hard cases (long flexible loops, unusual binding geometries) benefit from stronger forces. The guidance weight scale parameter allows per-complex tuning.

### Evidence from Round 1

Strategy Y (the predecessor of Y+) was the **only strategy out of 11 tested with a statistically significant structural improvement**:

- CDR-H3 RMSD: **−0.30 Å** (p=0.025) — the most clinically relevant metric
- DockQ (antibody-antigen): **+0.034** — consistent positive trend
- CAPRI Medium+ classification: **+6.4 percentage points**
- Zero degradation on any metric — the safest strategy tested
- Achieved results comparable to pocket-guided restraints (Strategy B2), which require oracle contact information that wouldn't be available in a real prediction scenario

The key insight from Round 1: the two-phase approach (embedding-space early, coordinate-space late) mirrors the physics of protein structure formation and is fundamentally sound. Y+ strengthens both phases and adds epitope targeting and per-complex adaptation.

---

## What does vanilla Boltz-2 do?

Boltz-2 predicts protein structures using a **diffusion process**: it starts from random noise and progressively denoises it into a 3D structure over ~200 steps. At each step, a neural network predicts the "clean" structure and nudges the noisy coordinates closer to it.

In vanilla Boltz-2, the model treats all residues equally — it has no special knowledge about which parts of an antibody are most important for binding, or where on the antigen surface the antibody should dock. It just predicts the most likely structure given the sequences.

## What does hierarchical steering add?

Hierarchical steering guides the diffusion process in **two phases**, each targeting a different aspect of the prediction:

```
Step 0                          Step 100                        Step 200
|--- EARLY PHASE ---------------|--- LATE PHASE -----------------|
  CDR beta-scaling (embedding)     Coordinate-space potentials
  "Shape the internal loops"       "Pull antibody toward antigen"
```

The intuition is: **first get the antibody's binding loops into good conformations, then guide the whole complex into a physically reasonable docked pose.**

---

## Phase 1: CDR Beta-Scaling (Early Phase)

### What are CDRs?

CDR (Complementarity-Determining Region) loops are the parts of the antibody that physically contact the antigen. There are six CDR loops (three on the heavy chain, three on the light chain), but CDR-H3 and CDR-L3 are the most variable and most important for binding specificity.

### What does beta-scaling do?

Inside Boltz-2's neural network, there is a matrix called `token_trans_bias` that influences how much attention each residue pair receives during the structure prediction. Beta-scaling **amplifies** the attention between CDR residue pairs:

```
token_trans_bias *= (1 + beta_t * cdr_pair_mask)
```

Where `cdr_pair_mask` is 1 for pairs where both residues are in CDR regions, and 0 otherwise.

**Intuition:** Think of it as telling the network "pay extra attention to how these CDR loops interact with each other." This encourages the network to refine the CDR loop conformations more carefully during the early denoising steps, when the overall structure is still being formed.

### Why time-varying?

The scaling factor `beta_t` decays over time rather than staying constant:

- **Linear schedule:** `beta_t = 0.8 * steering_t` — starts at 0.8, linearly drops to 0
- **Cosine schedule:** smooth S-curve decay
- **Step schedule:** full strength (0.8) for the first 30% of steps, 30% strength for the middle, then off

The reason: in the early steps, the structure is still very noisy and malleable — this is when CDR guidance is most useful. In the late steps, the backbone is mostly settled, and excessive bias would interfere with fine-grained refinement.

### Why embedding space and not coordinate space?

In the early steps, coordinates are extremely noisy (almost random). Applying forces in coordinate space at this stage would be pushing on noise. Instead, beta-scaling operates on the **internal representations** (embeddings) of the network, which carry meaningful structural information even when coordinates are noisy.

---

## Phase 2: Coordinate-Space Potentials (Late Phase)

Once the structure has taken rough shape (~halfway through denoising), hierarchical steering activates two coordinate-space potentials that apply physical forces to guide the docking:

### Potential 1: Antigen Orientation

**What it does:** Encourages the antigen surface to face the antibody's CDR loops.

**How it works:** For each antigen atom, the potential measures its distance to the nearest CDR atom. If the distance exceeds a threshold (default 8 Å), it applies a gradient that pulls the antigen closer. This is aggregated across all antigen-CDR pairs using a soft-union operation (so the potential is satisfied if the antigen is close to *any* CDR atom, not necessarily all of them).

**Phase schedule:**
- Steps 0–100 (steering_t > 0.5): **OFF** — weight = 0.0
- Steps 100–200 (steering_t ≤ 0.5): **ON** — weight = 1.5

### Potential 2: CDR-Antigen Proximity

**What it does:** Ensures CDR residues stay close to the antigen surface.

**How it works:** This is the reverse perspective — for each CDR atom, it checks whether there is a nearby antigen atom. If a CDR residue is floating far from the antigen, a gradient pulls it closer. Uses a tighter threshold (capped at 10 Å).

**Phase schedule:**
- Steps 0–120 (steering_t > 0.4): **OFF**
- Steps 120–200 (steering_t ≤ 0.4): **ON** — weight = 0.8

**Why two potentials?** The first ensures the antigen presents its binding surface. The second ensures the CDR loops actually reach toward that surface. Together, they create a docked pose from both sides.

---

## Epitope Information from Predicted PDBs (Y+.3)

### The problem

The antigen orientation potential, by default, considers **all** antigen surface atoms when computing distances to CDR loops. But real antibodies only bind a specific patch of the antigen called the **epitope** — typically 15–30 residues on the antigen surface.

Using the full antigen surface means the potential might be satisfied if the CDR loops are near *any* part of the antigen, even the wrong side. This dilutes the guidance signal.

### The solution

If you have a predicted structure of the complex (e.g., from a previous Boltz-2 run without steering), you can extract epitope information from it:

1. **Parse the predicted PDB**: find all antigen CA atoms and all CDR CA atoms
2. **Compute distances**: for each antigen residue, check if any CDR atom is within a contact threshold (default 10 Å)
3. **Extract epitope**: antigen residues that have at least one CDR contact are labeled as epitope residues

### How the epitope is used

When epitope residues are available, the **Antigen Orientation Potential** restricts its attention:

- **Without epitope:** considers distances from ALL antigen atoms to CDR atoms
- **With epitope:** considers distances from ONLY epitope antigen atoms to CDR atoms

This focuses the docking guidance on the correct binding site, preventing the antibody from docking to the wrong face of the antigen.

### How to provide epitope information

Pass the `--predicted_epitope_pdbs` flag pointing to a folder of previous prediction results:

```bash
boltz predict examples/yaml_hierarchical_steering/7TRH_HBG.yml \
  --use_potentials \
  --hierarchical_steering \
  --predicted_epitope_pdbs predictions_examples/antigen_cut
```

The PDB path is automatically constructed as:
```
<predicted_epitope_pdbs>/boltz_results_<code>/predictions/<code>/<code>_model_0.pdb
```

Where `<code>` is inferred from the YAML filename (e.g., `7TRH_HBG` from `7TRH_HBG.yml`), or can be overridden with the `--code` flag.

If the PDB file is not found, a warning is printed and the system falls back to using the full antigen surface — it does not fail.

---

## Guidance Weight Scaling (Y+.4)

### The problem

Not all antibody-antigen complexes are equally easy to predict. Some complexes have long, flexible CDR-H3 loops or unusual binding geometries. A one-size-fits-all guidance strength may be too weak for hard cases or too strong for easy ones.

### The solution

The `guidance_weight_scale` parameter (set per-complex in the YAML file) multiplies the energy and gradient of both coordinate-space potentials:

```
energy = base_energy * guidance_weight_scale
gradient = base_gradient * guidance_weight_scale
```

- `guidance_weight_scale = 1.0` — default, no change
- `guidance_weight_scale = 1.5` — 50% stronger guidance for hard cases
- `guidance_weight_scale = 0.5` — weaker guidance if default is too aggressive

This allows tuning per complex without changing the global potential parameters.

---

## Summary: What happens during a prediction

```
Step   steering_t   What's happening
─────  ──────────   ─────────────────────────────────────────────
  0      1.00       CDR beta-scaling at full strength (beta=0.8)
                    Network pays extra attention to CDR loops
                    Coordinate potentials are OFF

 50      0.75       Beta-scaling decaying (beta≈0.6 with linear)
                    Structure emerging from noise
                    Potentials still OFF

100      0.50       Beta-scaling weak (beta≈0.4)
                    Antigen Orientation potential turns ON (weight=1.5)
                    Pulls antigen toward CDR binding site

120      0.40       Beta-scaling very weak (beta≈0.3)
                    CDR Proximity potential turns ON (weight=0.8)
                    Pulls CDR loops toward antigen surface

150      0.25       Beta-scaling minimal (beta≈0.2)
                    Both potentials active — refining the docked pose

200      0.00       Beta-scaling OFF
                    Potentials driving final adjustments
                    Structure finalized
```

## CLI Usage

```bash
# Basic hierarchical steering (no epitope info)
boltz predict examples/yaml_hierarchical_steering/7TRH_HBG.yml \
  --use_potentials \
  --hierarchical_steering \
  --num_workers 0

# With epitope information from previous predictions
boltz predict examples/yaml_hierarchical_steering/7TRH_HBG.yml \
  --use_potentials \
  --hierarchical_steering \
  --predicted_epitope_pdbs predictions_examples/antigen_cut \
  --num_workers 0

# With multiple samples for better coverage
boltz predict examples/yaml_hierarchical_steering/7TRH_HBG.yml \
  --use_potentials \
  --hierarchical_steering \
  --predicted_epitope_pdbs predictions_examples/antigen_cut \
  --diffusion_samples 5 \
  --num_workers 0
```
