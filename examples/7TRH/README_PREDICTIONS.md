# 7TRH Antibody-Antigen Complex Predictions

Automated prediction suite for the 7TRH antibody-antigen complex with different steering strategies and constraints.

## Quick Start

```bash
cd examples/7TRH

# Full run (all cases, all samples)
bash run_all_predictions.sh

# Quick test (one case per setting, fast validation)
bash run_all_predictions.sh --quick-test

# Custom output directory
bash run_all_predictions.sh /path/to/output
bash run_all_predictions.sh --quick-test /path/to/output
```

Results saved to `./boltz_results/` (or custom directory).

---

## Prediction Cases

| Case | Setting | YAML | Flags | Samples |
|------|---------|------|-------|---------|
| 1 | Baseline | `7TRH_HBG.yml` | None | 5 |
| 2 | Antigen Steering (p=1) | `7TRH_HBG_antigen_steer.yml` | `--use_potentials --antigen_steering --num_particles 1` | 5 |
| 3 | Antigen Steering (p=5) | `7TRH_HBG_antigen_steer.yml` | `--use_potentials --antigen_steering --num_particles 5` | 1 |
| 4+ | Contact Constraints | `restraint_7TRH_HBG_*.yml` | `--use_potentials` | 5 |
| 5+ | Pocket Constraints | `restraint_to_A_*.yml` | `--use_potentials` | 5 |

---

## Case Details

### Case 1: Baseline (No Constraints)
Reference prediction using the model without steering.

### Cases 2-3: Antigen Steering
Optimize antigen position to maximize CDR-antigen contacts.

**Case 2**: Single particle gradient guidance (5 diverse samples)
**Case 3**: Multi-particle importance sampling (1 refined sample)

### Case 4: Contact Constraints
Enforce specific inter-residue distances:
- H-bond constraints
- Hydrophobic interactions
- Salt bridges

Example constraint:
```yaml
constraints:
- contact:
    token1: [B, 102]      # Heavy chain residue 102
    token2: [A, 143]      # Antigen residue 143
    max_distance: 5.0
    force: true
```

### Case 5: Pocket Constraints
Define binding pockets where a binder interacts with pocket residues.

Example constraint:
```yaml
constraints:
- pocket:
    binder: B             # Heavy chain is binder
    contacts:             # Pocket residues on antigen
      - [A, 140]
      - [A, 143]
      - [A, 145]
    max_distance: 8.0
    force: true
```

---

## Quick Test Mode

Use `--quick-test` to validate the pipeline:
```bash
bash run_all_predictions.sh --quick-test
```

Runs 1 representative case per setting:
- Case 1: Baseline (5 samples)
- Case 2: Antigen Steering p=1 (5 samples)
- Case 3: Antigen Steering p=5 (1 sample)
- Case 4: First contact constraint (5 samples)
- Case 5: First pocket constraint (5 samples)

**Expected runtime**: 20-30 minutes (vs. several hours for full run)

---

## Output Structure

```
boltz_results/
├── case1_baseline/
├── case2_antigen_steer_p1/
├── case3_antigen_steer_p5/
├── case4_restraint_7TRH_HBG_hbond_110/
├── case5_restraint_to_A_7TRH_HBG_B_W_109/
├── case6_restraint_to_A_7TRH_HBG_C_I_65/
├── ... (more cases)
├── case*.log
└── ... (log files)
```

Each case directory contains PDB prediction files.

---

## Customization

### Change number of samples
Edit script line 32:
```bash
NUM_SAMPLES=10  # instead of 5
```

### Run only specific cases
Comment out case sections in the script.

### Add new constraints
Simply add YAML files to the directory:
- `restraint_*.yml` - Automatically discovered as contact constraints
- `restraint_to_A_*.yml` - Automatically discovered as pocket constraints

Script auto-discovers new files on re-run.

---

## Reference

- [[../../STEERING_METHODS|Steering Methods]] - Full documentation
- [[../../antigen_steering|Antigen Steering Details]]
- [Prediction Documentation](https://github.com/octaviomtz/inspection/blob/steering_cdr3/docs/prediction.md)
