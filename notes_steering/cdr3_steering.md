# CDR3 Steering

**Flag**: `--cdr3_steering`

**Requirement**: Must be used with `--use_potentials`

---

## Goal

Control the backbone conformation of CDR3 (complementarity-determining region 3) loops in antibodies by enforcing specific dihedral angles. This is useful when you have experimental knowledge of desired loop geometries or want to test how well the model can maintain specific loop conformations while predicting the rest of the structure.

---

## What It Does (Intuitive Explanation)

Imagine the CDR3 loop is like a flexible rope. Without steering, the model predicts the most likely rope configuration based on its training. With CDR3 steering, you **pin certain parts of the rope** to specific angles, and the model must work around those constraints.

Specifically, CDR3 steering:

1. **Identifies** CDR3 residues in the antibody:
   - Heavy chain CDR3: residues 96-111
   - Light chain CDR3: residues 90-98

2. **Targets** the backbone dihedral angle (psi/ψ angle) at specific residues within CDR3

3. **Enforces** these angles during structure prediction by penalizing deviations through a potential energy function

The result is a predicted antibody structure where the CDR3 loops maintain the specified backbone angles while the rest of the structure adapts naturally.

---

## How It Works (Technical Details)

### Mechanism

CDR3 steering uses the `CDR3ConformationPotential` energy function, which:

1. **Computes** the current psi dihedral angle at each target residue from predicted coordinates

2. **Calculates** energy penalty: `E = k × (current_angle - target_angle)²`
   - Low energy when angles match target
   - High energy when angles deviate

3. **Applies** gradient guidance: Energy gradients modify atomic coordinates to minimize the penalty

4. **Iterates**: At each denoising step, coordinates are nudged closer to satisfying the target angles

### During Diffusion

```
Denoising Step t:
  1. Network predicts structure (x_0_pred)
  2. Compute psi angle from coordinates
  3. Calculate energy E = (psi - target)²
  4. Compute gradient ∇E with respect to coordinates
  5. Update coordinates: x_0_pred = x_0_pred - α·∇E
  6. Continue to next denoising step
```

The gradient acts like an invisible force pulling atoms into the desired backbone angles.

---

## Usage

### Basic Command

```bash
boltz predict config.yaml --use_potentials --cdr3_steering
```

### With Multiple Samples

```bash
boltz predict config.yaml --use_potentials --cdr3_steering --diffusion_samples 5
```

This generates 5 independent structures, all satisfying the CDR3 steering constraints but with different overall backbones.

### Configuration File

CDR3 target angles are typically specified in the input YAML config file. Example:

```yaml
# config.yaml
protein:
  sequence: "MKTAIL..."
  cdr3_targets:
    heavy_chain:
      - residue: 102
        psi_angle: -60.0  # degrees
      - residue: 105
        psi_angle: -45.0
    light_chain:
      - residue: 94
        psi_angle: -65.0
```

---

## When to Use

- **Antibody engineering**: Enforce known canonical CDR3 conformations
- **Loop design**: Test how well the model maintains specific geometries
- **Structure refinement**: Guide prediction toward experimentally validated loops
- **Conformational analysis**: Generate multiple samples with same CDR3 angles

## When NOT to Use

- If CDR3 conformation is unknown (let the model predict freely)
- If you need very high flexibility (steering restricts exploration)
- For proteins without CDR regions (not applicable)

---

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--cdr3_steering` | False | Enable CDR3 steering |
| `--diffusion_samples N` | 1 | Generate N independent samples |
| CDR3 targets (in config) | None | Target dihedral angles (required when steering enabled) |

---

## Output

When CDR3 steering is enabled, the output includes:

- **Structure files**: PDB format with predicted coordinates
- **CDR3 coordinates**: Heavy and light chain CDR3 loops matching target angles
- **Rest of structure**: Naturally predicted antibody framework and variable regions
- **Confidence scores**: pLDDT and PAE showing prediction confidence

---

## Examples

### Example 1: Enforce Specific CDR3 Geometry

```bash
boltz predict antibody_config.yaml \
  --output_format pdb \
  --use_potentials \
  --cdr3_steering \
  --diffusion_samples 3
```

Generates 3 structures where CDR3 has the specified geometry.

### Example 2: Test Known Loop Conformation

Test whether your experimental CDR3 geometry is compatible with predicted antigen binding:

```bash
# Use experimental CDR3 angles as targets
boltz predict ab_ag_complex.yaml \
  --use_potentials \
  --cdr3_steering
```

---

## Troubleshooting

**Q: My structure has clashes in the CDR3 region**
- The steering angles may be too restrictive. Try relaxing the target angles slightly or use a weaker steering weight.

**Q: CDR3 angles aren't being enforced**
- Ensure `--use_potentials` is enabled (CDR3 steering requires it)
- Check that CDR3 target angles are specified in the config file
- Verify residue numbering matches your input sequence

**Q: Results vary across diffusion samples**
- This is expected! Each sample uses different random noise initialization, so loop angles are enforced but loop side-chain packing may vary.

---

## See Also

- [[../STEERING_METHODS|Back to Steering Methods Overview]]
- [[Antigen Steering|antigen_steering]] for simultaneous antigen orientation control
- [[Steering Proposals|steering_proposals]] - Future strategies including Strategy D (Canonical Ensemble) for enhanced CDR3 conformation exploration
- `--use_potentials` documentation for general potential-guided inference
