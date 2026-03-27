# Evaluation Strategy: G+ v2 — Progressive Steering

**Date**: 2026-03-26
**Feature**: Strategy G+ v2 (Progressive Steering) — time-varying CDR3 beta-scaling with three-phase antigen orientation potential
**Branch**: `g_v2_progre`

---

## 1. Goals of the Feature

Strategy G+ v2 addresses three specific problems from Round 1's Strategy G:

| Sub-feature | Goal | What it changes |
|---|---|---|
| G+.2 — Time-varying beta | Stronger CDR3 beta-scaling early (exploration), decaying to zero late (convergence). Prevents late-stage beta from distorting already-converged structures | `token_trans_bias` scaled by `beta(t)` per diffusion step in `diffusionv2.py` |
| G+.3 — Interface-specific potential | Three-phase `ProgressiveAntigenOrientationPotential` with weak-early/strong-late guidance and FK resampling weights | `potentials.py`: CDR-antigen distance potential with progressive schedule |
| G+.1 — Conditional steering (skip when confident) | Eliminate regression on easy cases by skipping steering when baseline iPTM is already high | Not yet implemented — evaluation-time analysis only |

**Primary evaluation question**: Does progressive steering improve antibody-antigen interface quality (DockQ/CAPRI) compared to:
1. **B1** — No steering (plain Boltz-2 with antigen-cut MSA)
2. **B2** — Restraint-to-A (contact restraints to specific antigen residues)
3. **B3** — Contact restraints (hbond/hydrophobic/salt_bridge)
4. **G v1** — Original CDR3 beta-scaling (constant beta, fixed-phase potential) from `g_phase_refinement`

**Secondary questions**:
- Does G+ v2 avoid the regression on easy cases that G v1 showed (−0.025 DockQ)?
- Does the time-varying beta improve CDR3 loop confidence (pLDDT) and interface PAE?
- Is ensemble diversity maintained or reduced?
- How well does confidence-based model selection correlate with actual DockQ quality?

---

## 2. Data Layout

### 2.1 Ground Truth
- **Location**: `pdb_minimized/` (48 complexes)
- **Format**: PDB files named `{complex}.pdb` (e.g., `7TRH_HBG.pdb`)
- **Chain convention**: A = antigen, B = heavy chain, C = light chain

### 2.2 Predictions

All predictions live under `predictions_examples/` in method-specific subfolders:

| Method | Label | Folder | Folder naming pattern | Example |
|---|---|---|---|---|
| Baseline 1 (no steering) | B1 | `antigen_cut/` | `boltz_results_{complex}` | `boltz_results_7TRH_HBG/` |
| Baseline 2 (restraint-to-A) | B2 | `antigen_cut_vhvl_msa_pocket_ab_boltz_post_2023_omm/` | `boltz_results_restraint_to_A_{complex}_{chain}_{res}_{num}` | `boltz_results_restraint_to_A_7TRH_HBG_B_W_109/` |
| Baseline 3 (contact restraints) | B3 | `antigen_cut_contact_restraints/` | `boltz_results_restraint_{complex}_{type}_{num}` | `boltz_results_restraint_7TRH_HBG_hbond_23/` |
| G v1 (constant beta) | Gv1 | `new_feature_cdr3_beta/` | `boltz_results_{complex}_cdr3_beta` | `boltz_results_7TRH_HBG_cdr3_beta/` |
| **G+ v2 (progressive)** | **Gv2** | `new_feature_progressive_steering/` | `boltz_results_{complex}` | `boltz_results_7TRH_HBG/` |

**Note on B2/B3**: These have multiple sub-experiments per complex (different restraint residues/types). The evaluation script handles this by extracting a `sub_experiment` label and evaluating each separately, then selecting the best per complex.

### 2.3 Prediction Files (per complex per method)

Inside each `boltz_results_*/predictions/{name}/`:
- `{name}_model_{0..4}.pdb` or `.cif` — 5 structure predictions
- `confidence_{name}_model_{0..4}.json` — confidence metrics (iPTM, pTM, pLDDT, etc.)
- `plddt_{name}_model_{0..4}.npz` — per-residue pLDDT
- `pae_{name}_model_{0..4}.npz` — predicted aligned error matrix
- `pde_{name}_model_{0..4}.npz` — predicted distance error

**Important**: G v1 (cdr3_beta) currently has only 1 model per complex (CIF format). G+ v2 should be run with `--diffusion_samples 5` for proper ensemble evaluation.

---

## 3. Evaluation Metrics

### 3.1 Primary Metric: DockQ (via DockQ2 package)

**Tool**: [DockQ](https://github.com/wallnerlab/DockQ) — the standard tool for assessing protein interface quality.

**What it computes**:
- **DockQ** (0–1): Composite score combining Fnat, iRMSD, and LRMSD
- **Fnat**: Fraction of native contacts preserved
- **iRMSD**: Interface RMSD (Angstroms)
- **LRMSD**: Ligand RMSD after receptor alignment (Angstroms)
- **F1**: F1 score of native contacts
- **fnonnat**: Fraction of non-native contacts

**Per-interface breakdown**: DockQ reports metrics for each chain pair interface:
- **BA** (heavy–antigen): Most relevant for antibody docking
- **CA** (light–antigen): Secondary interface
- **BC** (heavy–light): Should be consistently good; regression here indicates structural damage

**CAPRI quality classification** (from DockQ score):
| Category | DockQ range |
|---|---|
| High | >= 0.80 |
| Medium | 0.49 – 0.80 |
| Acceptable | 0.23 – 0.49 |
| Incorrect | < 0.23 |

**Why DockQ**: DockQ internally handles structural alignment. It aligns the receptor (antigen) and computes ligand RMSD, which is exactly the right approach for antibody-antigen docking evaluation. No separate alignment preprocessing is needed.

### 3.2 Confidence Metrics

Extracted from Boltz-2 confidence JSON files:
- `confidence_score`: Overall model confidence
- `iptm`: Interface predicted TM-score (key for ranking)
- `ptm`: Predicted TM-score
- `complex_plddt`: Overall predicted LDDT
- `complex_iplddt`: Interface pLDDT
- `pair_chains_iptm`: Pairwise chain iPTM (AB, AC, BC)

### 3.3 Per-Residue Metrics

From NPZ files, using CDR indices from `examples/cdrs.csv`:
- **pLDDT per region**: antigen, heavy chain, light chain, CDR3-H, CDR3-L
- **Interface PAE**: Mean PAE between antibody and antigen residue blocks
- **CDR3-antigen PAE**: Mean PAE between CDR3 residues and antigen

These are particularly relevant for G+ v2: if time-varying beta works correctly, CDR3 pLDDT should improve (less distortion from late-stage beta) while interface PAE should decrease (better CDR-antigen spatial prediction).

### 3.4 Ensemble Diversity

Pairwise CA-RMSD across the 5 models per complex:
- Mean and standard deviation of pairwise RMSD
- Computed over all chains (ABC) after superposition

**Why**: Progressive steering might reduce diversity if the potential over-constrains the ensemble. This metric detects that.

---

## 4. Model Selection Strategies

For each method-complex combination (5 models), select the "best" model using three strategies:

| Strategy | Selection criterion | Purpose |
|---|---|---|
| **Oracle** | Best actual DockQ | Upper bound — what's the best the method can produce? |
| **Confidence** | Highest `confidence_score` from Boltz-2 | Practical selection — what users would do without ground truth |
| **iPTM** | Highest `iptm` | Alternative ranking — interface-specific confidence |

**For B2/B3** (multiple sub-experiments per complex): First select the best model within each sub-experiment, then select the best sub-experiment per complex. Report both the best-sub-experiment result and the all-sub-experiments aggregate.

---

## 5. Statistical Tests

### 5.1 Paired Wilcoxon Signed-Rank Test
- Compare each method against B1 (reference) on confidence-selected DockQ
- Paired by complex (same 48 complexes evaluated across all methods)
- Reports: test statistic, p-value, median delta, mean delta

### 5.2 Per-Complex Delta Analysis
- For each complex: `delta_DockQ = method_DockQ - B1_DockQ`
- Categorize complexes: improved (delta > 0.05), unchanged (|delta| <= 0.05), regressed (delta < -0.05)
- **Key for G+ v2**: Count regressions — the main goal is fewer regressions than G v1

---

## 6. Plots

1. **DockQ bar chart**: Mean DockQ per method (confidence-selected) with error bars
2. **DockQ scatter**: G+ v2 vs B1, G+ v2 vs G v1 — points above diagonal = improvement
3. **CAPRI stacked bar**: Distribution of quality categories per method
4. **Confidence vs DockQ**: Scatter plot to assess if confidence-based selection works
5. **Ensemble diversity box plot**: Pairwise CA-RMSD distribution per method
6. **Per-complex comparison**: Sorted bar chart showing all methods side by side
7. **Delta histogram**: Distribution of DockQ changes (G+ v2 minus B1) — should be right-shifted

---

## 7. Script Reuse Plan

The evaluation script from `g_phase_refinement` branch (`../../worktree_steering_cdr3/g_phase_refinement/evaluation/evaluate.py`) can be reused with the following modifications:

### 7.1 Changes Required

1. **Update `METHOD_DEFS`**: Add G+ v2 entry, rename G v1 label
   ```python
   METHOD_DEFS = {
       "B1": {
           "folder": "antigen_cut",
           "prefix": "boltz_results_",
           "suffix": "",
       },
       "B2": {
           "folder": "antigen_cut_vhvl_msa_pocket_ab_boltz_post_2023_omm",
           "prefix": "boltz_results_restraint_to_A_",
           "suffix": "",
       },
       "B3": {
           "folder": "antigen_cut_contact_restraints",
           "prefix": "boltz_results_restraint_",
           "suffix": "",
       },
       "Gv1": {   # Renamed from "Gplus"
           "folder": "new_feature_cdr3_beta",
           "prefix": "boltz_results_",
           "suffix": "_cdr3_beta",
       },
       "Gv2": {   # NEW — progressive steering
           "folder": "new_feature_progressive_steering",
           "prefix": "boltz_results_",
           "suffix": "",
       },
   }
   ```

2. **Update color map**: Add Gv2 color, differentiate from Gv1
   ```python
   colors = {
       "B1": "#1f77b4",   # blue
       "B2": "#ff7f0e",   # orange
       "B3": "#2ca02c",   # green
       "Gv1": "#9467bd",  # purple
       "Gv2": "#d62728",  # red
   }
   ```

3. **Update plot ordering**: Change `["B1", "B2", "B3", "Gplus"]` to `["B1", "B2", "B3", "Gv1", "Gv2"]` throughout

4. **Add G+ v2 vs G v1 scatter plot**: New comparison plot showing improvement from v1 to v2

### 7.2 What Stays the Same

- Discovery logic (`extract_complex_name` matching against known complexes) — works as-is
- DockQ computation via DockQ Python API — no changes needed
- Confidence extraction — same JSON format
- Per-residue metrics (pLDDT, PAE) — same NPZ format
- Ensemble diversity (pairwise CA-RMSD) — same logic
- Model selection (oracle, confidence, iPTM) — same logic
- Statistical tests (Wilcoxon) — same logic
- CDR loading from `examples/cdrs.csv` — same format

---

## 8. Environment

### 8.1 Primary Environment: `dockq2`
```bash
conda activate dockq2
```
Contains: DockQ package, BioPython, NumPy, Pandas, SciPy, matplotlib

### 8.2 Required Packages
- `DockQ` (from https://github.com/wallnerlab/DockQ) — interface quality scoring
- `biopython` — PDB/CIF parsing, structural superposition
- `numpy`, `pandas` — data manipulation
- `scipy` — statistical tests (Wilcoxon)
- `matplotlib` — plots

No new environment is needed. The `dockq2` conda environment has all required packages.

---

## 9. Running the Evaluation

### 9.1 Generate Predictions (on compute machine)

```bash
# Activate boltz environment
conda activate boltz

# Run progressive steering on all 48 complexes
for yml in yaml_test_set_progressive_steering/*.yml; do
    python -m boltz.main predict "$yml" \
        --out_dir predictions_examples/new_feature_progressive_steering \
        --use_potentials \
        --progressive_steering \
        --diffusion_samples 5 \
        --devices 1
done
```

### 9.2 Run Evaluation

```bash
# Activate dockq2 environment
conda activate dockq2

# Run full evaluation (all methods)
python evaluation/evaluate.py \
    --pred_dir predictions_examples \
    --native_dir pdb_minimized \
    --output_dir evaluation/results \
    --cdrs_csv examples/cdrs.csv

# Or evaluate specific methods only
python evaluation/evaluate.py \
    --pred_dir predictions_examples \
    --native_dir pdb_minimized \
    --output_dir evaluation/results \
    --methods B1 Gv1 Gv2
```

### 9.3 Quick evaluation (skip DockQ if already computed)
```bash
python evaluation/evaluate.py \
    --pred_dir predictions_examples \
    --native_dir pdb_minimized \
    --output_dir evaluation/results \
    --skip_dockq
```

---

## 10. Success Criteria

G+ v2 is considered successful if:

| Criterion | Threshold | Rationale |
|---|---|---|
| Mean DockQ >= B1 | delta >= 0.0 | Must not regress overall |
| Mean DockQ > Gv1 | delta > 0.0 | Must improve on v1 |
| Fewer regressions than Gv1 | < 25% complexes with delta < -0.05 vs B1 | Key G+.1 motivation |
| CDR3 pLDDT >= Gv1 | delta >= 0.0 | Time-varying beta should not hurt CDR3 confidence |
| Interface PAE <= B1 | delta <= 0.0 (lower is better) | Better spatial prediction at interface |
| Medium+ CAPRI % >= B1 | >= 0% delta | No loss of docking quality classification |

---

## 11. Output Files

The evaluation produces these files in `evaluation/results/`:

| File | Description |
|---|---|
| `discovered_predictions.csv` | All found prediction files with paths |
| `dockq_results.csv` | DockQ metrics for every model |
| `confidence_metrics.csv` | Confidence scores for every model |
| `perres_metrics.csv` | Per-residue pLDDT and PAE metrics |
| `diversity_metrics.csv` | Ensemble diversity (pairwise RMSD) |
| `selected_models.csv` | Best model per method-complex (oracle + confidence) |
| `summary_table.csv` | Aggregate summary per method |
| `wilcoxon_tests.csv` | Statistical significance tests |
| `per_complex_comparison.csv` | Side-by-side DockQ per complex per method |
| `figures/dockq_bar.png` | Mean DockQ bar chart |
| `figures/dockq_scatter_gv2_vs_b1.png` | Scatter plot G+ v2 vs B1 |
| `figures/capri_stacked_bar.png` | CAPRI quality distribution |
| `figures/confidence_vs_dockq.png` | Confidence correlation |
| `figures/ensemble_diversity_boxplot.png` | Diversity comparison |
| `figures/dockq_per_complex.png` | Per-complex side-by-side bars |
