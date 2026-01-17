# CDR3 Conformation Steering in Boltz2

## Goal

CDR3 (Complementarity Determining Region 3) is the most variable and structurally diverse loop in antibodies, directly responsible for antigen recognition. This steering potential allows users to bias CDR3 loop conformations toward specific structural states during diffusion sampling, enabling:

- Generation of antibodies with desired binding geometries
- Exploration of alternative CDR3 conformations for the same sequence
- Structural optimization for therapeutic antibody design

## How It Works

The steering uses **backbone psi dihedral angles** (N-CA-C-N') to control CDR3 loop geometry. During diffusion sampling:

1. **Featurization**: Psi dihedral atom indices are extracted for specified CDR3 residues
2. **Energy Computation**: A flat-bottom potential penalizes dihedrals outside target ranges
3. **Guidance**: Feynman-Kac resampling and gradient updates steer the structure toward target conformations

The potential integrates with Boltz2's existing steering infrastructure, applying constraints during the denoising process without retraining the model.

## When Is This Useful?

### Structural Scenarios

| Conformation | Psi Range | Structural Character | Use Case |
|-------------|-----------|---------------------|----------|
| **Extended** | 120° to 160° | Elongated, reaching loop | Large/convex epitopes, enzyme active sites |
| **Compact** | -60° to -20° | Tight turn, retracted | Small/concave epitopes, peptide binding |
| **Kinked** | -20° to 20° | Sharp bend at apex | Cavity penetration, cryptic epitopes |

### Applications

- **Epitope targeting**: Match CDR3 shape to known epitope geometry
- **Conformational sampling**: Generate diverse structures for docking/MD
- **Loop engineering**: Test how sequence changes affect achievable conformations
- **Bispecific design**: Ensure compatible CDR3 geometries across binding arms

## Usage

### Command Line

```bash
# Basic usage with extended CDR3 conformation
boltz predict input.yaml --use_potentials --cdr3_steering --out_dir ./output

# With step-by-step mode for intermediate outputs
boltz predict input.yaml --use_potentials --cdr3_steering --step_by_step --out_dir ./output_stepwise
```

### YAML Constraint Format

```yaml
constraints:
  # Predefined conformation
  - cdr3_conformation:
      chain_id: H        # Chain ID
      start_res: 95      # Start residue (1-indexed)
      end_res: 106       # End residue (1-indexed)
      conformation: extended  # extended, compact, kinked, or custom
      force: true

  # Custom angle ranges (in degrees)
  - cdr3_conformation:
      chain_id: H
      start_res: 95
      end_res: 106
      conformation: custom
      psi_lower: [120, 125, 130, 120, 125, 130, 120, 125, 130, 120]
      psi_upper: [160, 155, 160, 160, 155, 160, 160, 155, 160, 160]
      force: true
```

### Predefined Conformations

| Conformation | Psi Angle Range | Description |
|-------------|-----------------|-------------|
| `extended` | 120° to 160° | Elongated CDR3 loop |
| `compact` | -60° to -20° | Tight turn |
| `kinked` | -20° to 20° | Sharp bend at apex |
| `custom` | User-defined | Specific angle targets per residue |

## Files Modified

| File | Changes |
|------|---------|
| `src/boltz/model/potentials/potentials.py` | `CDR3ConformationPotential` class |
| `src/boltz/data/types.py` | `cdr3_constraints` in `InferenceOptions` |
| `src/boltz/data/parse/schema.py` | YAML parsing for `cdr3_conformation` |
| `src/boltz/data/feature/featurizerv2.py` | `process_cdr3_feature_constraints()` |
| `src/boltz/data/module/inferencev2.py` | Constraint wiring |
| `src/boltz/main.py` | `--cdr3_steering` CLI flag |

## Example

See `examples/cdr3_steering_example.yaml` for a complete working example.
