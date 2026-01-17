# Boltz2 antibody CDR3 steering potential
_Exported on 1/17/2026 at 22:13:41 GMT+1 from Cursor (2.3.34)_

---

**User**

Here is the repository of boltz2 (https://www.biorxiv.org/content/10.1101/2025.06.14.659707v1.full.pdf ) boltz2 uses --use_potentials as a Inference‑time steering potentials. Inference‑time potentials are scoring functions (Feynman–Kac ) that tell the sampler which partial trajectories are “on track” to yield high‑reward samples, allowing the algorithm to selectively amplify promising diffusion paths without retraining the model. Check how to define and integrate a new steering potential to (for instance) force the CDR3 of an antibody to adopt different comformations

---

**Cursor**

I'll help you understand how to define and integrate a new steering potential in Boltz2. Let me explore the codebase to see how the existing potentials are structured and integrated.

Now let me check what features are available that could help identify CDR3 regions:

Let me look at the parse/yaml.py to understand how constraints are defined:

Based on my exploration of the Boltz2 codebase, I can provide you with a comprehensive guide on how to define and integrate a new steering potential to control CDR3 conformations in antibodies.

## Overview of Boltz2 Steering Potentials

Boltz2 implements **Feynman-Kac potentials** as inference-time steering mechanisms. These are scoring functions that guide the diffusion sampling process without retraining the model. The system supports two types of steering:

1. **FK Steering (Feynman-Kac resampling)**: Evaluates energies and resamples particles based on potential scores
2. **Guidance Updates**: Applies gradient-based updates to steer coordinates during diffusion

## Architecture of Potentials

### Key Components

All potentials inherit from the base `Potential` class in `src/boltz/model/potentials/potentials.py`:

```python
class Potential(ABC):
    def __init__(self, parameters: Optional[Dict[str, Union[ParameterSchedule, float, int, bool]]] = None):
        self.parameters = parameters
    
    @abstractmethod
    def compute_function(self, value, *args, negation_mask=None, compute_derivative=False):
        """Compute the energy and optionally its derivative"""
        
    @abstractmethod
    def compute_variable(self, coords, index, compute_gradient=False):
        """Compute the geometric variable (distance, dihedral, etc.)"""
        
    @abstractmethod
    def compute_args(self, feats, parameters):
        """Extract relevant indices and parameters from features"""
```

### Existing Potential Types

1. **DistancePotential**: Controls pairwise distances
2. **DihedralPotential**: Controls torsion angles
3. **ReferencePotential**: Enforces similarity to reference structure
4. **FlatBottomPotential**: Penalizes values outside bounds

## Creating a CDR3 Conformation Potential

Here's how to implement a custom potential to control CDR3 conformations:

### Step 1: Define the Potential Class

```python
# Add to src/boltz/model/potentials/potentials.py

class CDR3ConformationPotential(FlatBottomPotential, DihedralPotential):
    """
    Potential to control CDR3 loop conformations via backbone dihedral angles.
    
    This potential allows you to:
    - Enforce specific phi/psi angles for CDR3 residues
    - Sample different CDR3 conformations by varying angle targets
    - Guide the CDR3 loop into extended, compact, or specific conformations
    """
    
    def compute_args(self, feats, parameters):
        """
        Extract CDR3 dihedral angle indices and target bounds.
        
        Parameters expected in feats:
        - cdr3_dihedral_index: [4, N] tensor of atom indices for phi/psi dihedrals
        - cdr3_target_lower: [N] tensor of lower angle bounds (radians)
        - cdr3_target_upper: [N] tensor of upper angle bounds (radians)
        
        These would be precomputed during input processing based on:
        1. CDR3 region identification (e.g., residues 95-102 in IMGT numbering)
        2. Desired conformation (extended: phi=-120°, psi=140°; 
                              compact: phi=-80°, psi=-40°; etc.)
        """
        
        # Check if CDR3 features are present
        if "cdr3_dihedral_index" not in feats:
            # No CDR3 steering requested
            return torch.empty([4, 0]), None, None, None, None
        
        dihedral_index = feats["cdr3_dihedral_index"][0]  # [4, N_dihedrals]
        
        # Get target angle bounds (with flexibility buffer)
        lower_bounds = feats["cdr3_target_lower"][0].clone()  # [N_dihedrals]
        upper_bounds = feats["cdr3_target_upper"][0].clone()  # [N_dihedrals]
        
        # Apply buffer from parameters (allows some flexibility)
        angular_buffer = parameters.get("angular_buffer", 0.35)  # ~20 degrees default
        lower_bounds -= angular_buffer
        upper_bounds += angular_buffer
        
        # Spring constant for energy penalty
        k = torch.ones_like(lower_bounds) * parameters.get("force_constant", 1.0)
        
        return dihedral_index, (k, lower_bounds, upper_bounds), None, None, None
```

### Step 2: Add CDR3 Feature Extraction

You need to add CDR3 identification and target angle specification during input processing. Add to `src/boltz/data/parse/yaml.py`:

```python
@dataclass
class CDR3Constraint:
    """CDR3 conformation constraint."""
    chain_id: str
    start_res: int  # CDR3 start residue index
    end_res: int    # CDR3 end residue index
    conformation: str  # "extended", "compact", "kinked", or "custom"
    phi_targets: Optional[list[float]] = None  # Custom phi angles
    psi_targets: Optional[list[float]] = None  # Custom psi angles
    force: bool = True  # Enable potential-based enforcement
```

### Step 3: Precompute CDR3 Features

In the featurization stage (`src/boltz/data/feature/featurizerv2.py`), add:

```python
def compute_cdr3_features(self, record, cdr3_constraint):
    """
    Compute CDR3 dihedral indices and target angles.
    
    Returns features to be added to the batch:
    - cdr3_dihedral_index: [4, N] - indices for N-CA-C-N and CA-C-N-CA dihedrals
    - cdr3_target_lower: [N] - lower bounds for each dihedral
    - cdr3_target_upper: [N] - upper bounds for each dihedral
    """
    
    # Find CDR3 residues in the specified chain
    cdr3_residues = []
    for chain in record.chains:
        if chain.entity_id == cdr3_constraint.chain_id:
            for res in chain.residues:
                if cdr3_constraint.start_res <= res.idx <= cdr3_constraint.end_res:
                    cdr3_residues.append(res)
    
    # Build dihedral indices (N-CA-C-N for phi, CA-C-N-CA for psi)
    dihedral_indices = []
    target_lower = []
    target_upper = []
    
    for i, res in enumerate(cdr3_residues[:-1]):
        # Get backbone atom indices
        n_idx = self.get_atom_index(res, "N")
        ca_idx = self.get_atom_index(res, "CA")
        c_idx = self.get_atom_index(res, "C")
        next_n_idx = self.get_atom_index(cdr3_residues[i+1], "N")
        
        # Phi dihedral: C(i-1) - N(i) - CA(i) - C(i)
        # Psi dihedral: N(i) - CA(i) - C(i) - N(i+1)
        dihedral_indices.append([n_idx, ca_idx, c_idx, next_n_idx])
        
        # Set target angles based on conformation type
        if cdr3_constraint.conformation == "extended":
            target_lower.append(-2.27)  # -130°
            target_upper.append(-1.92)  # -110°
        elif cdr3_constraint.conformation == "compact":
            target_lower.append(-1.57)  # -90°
            target_upper.append(-1.22)  # -70°
        # ... etc
    
    return {
        "cdr3_dihedral_index": torch.tensor(dihedral_indices).T,
        "cdr3_target_lower": torch.tensor(target_lower),
        "cdr3_target_upper": torch.tensor(target_upper),
    }
```

### Step 4: Register the Potential

Modify the `get_potentials` function in `src/boltz/model/potentials/potentials.py`:

```python
def get_potentials(steering_args, boltz2=False):
    potentials = []
    
    # Existing potentials...
    if steering_args["fk_steering"] or steering_args["physical_guidance_update"]:
        potentials.extend([
            # ... existing potentials ...
        ])
    
    # Add CDR3 potential if requested
    if steering_args.get("cdr3_steering", False):
        potentials.append(
            CDR3ConformationPotential(
                parameters={
                    "guidance_interval": 2,  # Apply every 2 steps
                    "guidance_weight": PiecewiseStepFunction(
                        thresholds=[0.5],
                        values=[1.0, 0.5]  # Stronger early, relax later
                    ),
                    "resampling_weight": 1.0,
                    "angular_buffer": 0.35,  # ~20° flexibility
                    "force_constant": 1.0,
                }
            )
        )
    
    return potentials
```

### Step 5: YAML Input Format

Users would specify CDR3 constraints in their YAML input:

```yaml
version: 1
sequences:
  - protein:
      id: A
      sequence: EVQLV...  # Heavy chain sequence
      msa: heavy_chain.a3m
  - protein:
      id: B
      sequence: DIQMT...  # Light chain sequence
      msa: light_chain.a3m

constraints:
  - cdr3_conformation:
      chain_id: A
      start_res: 95  # IMGT numbering
      end_res: 102
      conformation: extended  # or "compact", "kinked", "custom"
      force: true
  
  # For custom angles:
  - cdr3_conformation:
      chain_id: A
      start_res: 95
      end_res: 102
      conformation: custom
      phi_targets: [-120, -125, -110, -115, -120, -118, -122]
      psi_targets: [140, 135, 145, 138, 142, 140, 136]
      force: true
```

### Step 6: Command Line Usage

```bash
# Enable potentials with CDR3 steering
boltz predict antibody.yaml \
    --use_potentials \
    --out_dir ./output \
    --diffusion_samples 5 \
    --sampling_steps 200
```

## How It Works

1. **During Preprocessing**: CDR3 regions are identified, and target dihedral angles are computed based on desired conformation
2. **During Diffusion**: At each sampling step:
   - The potential computes current phi/psi angles for CDR3 residues
   - Energy is computed based on deviation from targets
   - **FK Resampling**: Samples with lower energy are preferentially kept
   - **Gradient Updates**: Gradients push coordinates toward target angles
3. **Parameter Schedules**: Steering strength can vary during diffusion (e.g., stronger early, relax later)

## Advantages of This Approach

- **No retraining needed**: Works with existing checkpoints
- **Flexible**: Can sample different CDR3 conformations by changing targets
- **Physically grounded**: Respects backbone geometry constraints
- **Combinable**: Works alongside other potentials (clashes, contacts, etc.)

## Alternative: Distance-Based CDR3 Potential

If you want to control CDR3 shape via distances instead of dihedrals:

```python
class CDR3ShapePotential(FlatBottomPotential, DistancePotential):
    """Control CDR3 loop compactness via Cα-Cα distances"""
    
    def compute_args(self, feats, parameters):
        # Get CDR3 Cα-Cα distance pairs
        pair_index = feats["cdr3_ca_pairs"][0]  # [2, N_pairs]
        
        # Target distances (Å)
        target_distances = feats["cdr3_target_distances"][0]
        
        # Set bounds with tolerance
        tolerance = parameters.get("distance_tolerance", 2.0)
        lower_bounds = target_distances - tolerance
        upper_bounds = target_distances + tolerance
        
        k = torch.ones_like(lower_bounds)
        
        return pair_index, (k, lower_bounds, upper_bounds), None, None, None
```

This approach would let you specify things like "CDR3 apex should be 15Å from base" to control loop extension.

## Key Parameters to Tune

- **guidance_weight**: How strongly to apply the potential (0.1-2.0 typical range)
- **resampling_weight**: Importance for FK particle selection
- **guidance_interval**: How often to apply (1-5 steps)
- **angular_buffer/distance_tolerance**: Flexibility around targets
- **Schedule**: Use `ExponentialInterpolation` or `PiecewiseStepFunction` to vary strength during diffusion

This framework gives you fine-grained control over antibody CDR3 conformations during inference, enabling systematic exploration of the conformational landscape!

