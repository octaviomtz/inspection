# Antigen Steering

**Flags**: `--antigen_steering` + `--num_particles` (optional)

**Requirement**: Must be used with `--use_potentials`

---

## Goal

Optimize the 3D orientation and position of antigens relative to antibody CDR loops to maximize productive binding contacts. This is useful when you want to predict antibody-antigen complexes with geometries that emphasize CDR-antigen interactions, or when you want to explore how different antigen orientations affect binding.

---

## What It Does (Intuitive Explanation)

Imagine an antibody CDR region and an antigen as two puzzle pieces. Without steering, the model positions them in a reasonable configuration based on its training. With antigen steering, you **keep rotating and translating the antigen** to find orientations where the CDR loops make the most contact with the antigen surface.

Specifically, antigen steering:

1. **Scores** how well CDR residues contact the antigen at each step
   - Good contact = low energy
   - Poor contact = high energy

2. **Applies** energy gradients that push the antigen toward better positions
   - Antigen moves closer to CDR loops
   - CDR loops orient to face the antigen
   - Geometry optimizes for binding

3. **Samples** multiple particle trajectories (optional) to explore different binding modes
   - Each particle starts with different random noise
   - Particles with poor contacts are downweighted
   - Final structure uses highest-scoring conformation

The result is a predicted complex where the antigen is positioned to maximize CDR-antigen interface contacts.

---

## How It Works (Technical Details)

### Mechanism

Antigen steering uses the `AntigenOrientationPotential` energy function, which:

1. **Identifies** CDR residues:
   - Heavy chain: residues 31-35, 50-65, 95-102 (CDR1, CDR2, CDR3)
   - Light chain: residues 24-34, 50-56, 89-97 (CDR1, CDR2, CDR3)

2. **Computes** contact score between CDR atoms and antigen atoms
   - Distance-based: nearby atoms = better contact
   - Soft threshold: smooth distance function (not hard cutoff)
   - Aggregated over all CDR-antigen pairs

3. **Calculates** energy: `E = -contact_score` (negative because we minimize)
   - Low energy = strong CDR-antigen contacts
   - High energy = poor contacts

4. **Applies** gradient guidance and/or importance sampling:
   - **Gradient mode**: Coordinates modified to improve contact
   - **Particle mode**: Multiple samples scored, best ones kept

### During Diffusion

```
Denoising Step t:
  1. Network predicts structure (x_0_pred)
  2. Compute contact score between CDR and antigen atoms
  3. Calculate energy E = -contact_score
  4. Compute gradient ∇E with respect to all coordinates
  5. Update: x_0_pred = x_0_pred - α·∇E
     (This pushes antigen closer and CDR loops toward it)
  6. Continue to next denoising step
```

OR (with particles):

```
Denoising Step t (with K particles):
  1. K network predictions with different random noise
  2. For each particle: compute contact score
  3. Compute resampling weights: w_i ∝ exp(-E_i)
  4. Resample: keep high-weight particles, discard low-weight
  5. Continue to next denoising step
```

---

## Parameters

### Required

| Parameter | Description |
|-----------|-------------|
| `--antigen_steering` | Enable antigen steering |
| `--use_potentials` | Required; enables potential energy functions |

### Optional

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--num_particles` | 1 | Number of parallel particle trajectories |
| `--diffusion_samples` | 1 | Number of independent final structures |

### Interaction

- **No `--num_particles`**: Gradient guidance only (faster, simpler)
- **`--num_particles N` (N > 1)**: Feynman-Kac resampling (more thorough exploration, more compute)

---

## Usage

### Basic Gradient Guidance

```bash
boltz predict antibody_antigen.yaml --use_potentials --antigen_steering
```

Applies gradient-based steering to optimize antigen-CDR contacts.

### With Particle Resampling

```bash
boltz predict antibody_antigen.yaml \
  --use_potentials \
  --antigen_steering \
  --num_particles 5
```

Explores 5 parallel trajectories, resampling at each step based on contact quality. More thorough but slower.

### Multiple Final Samples

```bash
boltz predict antibody_antigen.yaml \
  --use_potentials \
  --antigen_steering \
  --num_particles 3 \
  --diffusion_samples 5
```

For each of 5 independent samples, uses 3 particles to explore different binding modes.

### Combined with CDR3 Steering

```bash
boltz predict antibody_antigen.yaml \
  --use_potentials \
  --cdr3_steering \
  --antigen_steering \
  --num_particles 3
```

Simultaneously enforce CDR3 backbone angles AND optimize antigen orientation.

---

## How Particles Work

With `--num_particles N`:

1. **Initialization**: N independent copies of random noise initialization
2. **Each denoising step**:
   - All N particles denoised in parallel
   - Each scored by antigen contact potential
   - Resampled: particles with better contacts get higher weight
   - Worse particles may be eliminated or replicated
3. **Final output**: Single best particle selected

This is an importance sampling technique that explores more of the conformational space while filtering for good binding geometries.

### Performance Considerations

| Particles | Speed | Exploration | Recommended For |
|-----------|-------|-------------|-----------------|
| 1 | Fast | Limited | Single binding mode, quick results |
| 3 | Medium | Good | Standard usage, good balance |
| 5 | Slow | Thorough | Critical applications, multiple modes |

---

## When to Use

- **Complex prediction**: Antibody-antigen complexes where binding geometry is important
- **Binding optimization**: Improve CDR-antigen interface contacts
- **Orientation exploration**: Sample different binding modes
- **Interface design**: Optimize for strong, specific interactions

## When NOT to Use

- If antigen position is known (use it as a fixed template instead)
- For monomeric antibody structure prediction (no antigen)
- If you need to explore non-binding orientations
- For very large antigens (computational cost)

---

## Output

When antigen steering is enabled, the output includes:

- **Complex structure**: Antibody-antigen complex with optimized orientation
- **CDR-antigen interface**: Maximized contact surface
- **Confidence scores**: pLDDT, PAE, and pDE metrics
- **Contact statistics** (if available): Number/quality of CDR-antigen contacts

---

## Examples

### Example 1: Optimize Binding Geometry

```bash
boltz predict ab_antigen_complex.yaml \
  --output_format pdb \
  --use_potentials \
  --antigen_steering \
  --num_particles 3 \
  --diffusion_samples 2
```

Generates 2 structures with optimized CDR-antigen contacts (using 3 particles per sample).

### Example 2: Explore Multiple Binding Modes

```bash
boltz predict ab_antigen_complex.yaml \
  --use_potentials \
  --antigen_steering \
  --diffusion_samples 10
```

Generates 10 independent samples, each with optimized but potentially different binding geometries.

### Example 3: Combined CDR3 and Antigen Steering

```bash
boltz predict ab_antigen_complex.yaml \
  --use_potentials \
  --cdr3_steering \
  --antigen_steering \
  --num_particles 5
```

Enforce specific CDR3 backbone angles while simultaneously optimizing antigen orientation.

---

## Troubleshooting

**Q: Antigen isn't getting closer to CDR loops**
- Ensure `--use_potentials` is enabled
- Check that sequences include both antibody and antigen
- Try increasing `--num_particles` for better exploration
- Verify CDR residue numbering is correct

**Q: Results are different each time**
- This is expected! Random diffusion noise initialization causes variation. Use `--diffusion_samples N` to generate multiple samples.

**Q: Computation is very slow**
- High `--num_particles` increases cost. Try reducing to 1-3.
- Check if antigen is very large; consider using a smaller antigen region.

**Q: CDR3 steering conflicts with antigen steering**
- This shouldn't happen; they operate on different residues. If getting unexpected results, check your CDR3 target angles aren't extremely restrictive.

---

## See Also

- [[../STEERING_METHODS|Back to Steering Methods Overview]]
- [[CDR3 Steering|cdr3_steering]] for CDR3 backbone control
- `--use_potentials` documentation for general potential-guided inference
- [Boltz-2 Inference Pipeline](../boltz2_inference_pipeline.md)
