# Evaluation Strategy — Strategy Y+ (Hierarchical Steering v2)

## Can we reuse the existing evaluation script?

**Yes, with a small modification.**

The script at `../../worktree_steering_cdr3/y_hierarchical_steering/scripts/eval/evaluate_hierarchical_steering.py` was written to evaluate exactly this class of experiments. It already handles all four methods, computes all relevant metrics, and produces statistical comparisons. The only change needed is to update the `METHOD_CONFIGS` dictionary to point to the correct folder name for the new feature.

**Current entry for hierarchical:**
```python
"hierarchical": {
    "subfolder": "new_feature_hierarchical",       # <-- update this
    "prefix": "boltz_results_",
    "suffix": "_hierarchical",                     # <-- update or remove if not used
    "multi_setting": False,
},
```

**Updated entry for Y+ (exact folder name TBD after running predictions):**
```python
"hierarchical_v2": {
    "subfolder": "hierarchical_steering_v2",       # match actual predictions folder name
    "prefix": "boltz_results_",
    "suffix": "",
    "multi_setting": False,
},
```

No other changes are needed. The prediction folder structure (`boltz_results_{complex}/predictions/{complex}/{complex}_model_*.pdb`) is identical.

---

## Feature Goals (What We Are Evaluating)

Strategy Y+ has a primary goal and three sub-goals that translate directly into metrics:

| Goal | What to measure |
|------|----------------|
| **Primary**: Improve antibody-antigen docking accuracy | DockQ (Ab-Ag interfaces), CAPRI classification |
| **Sub-goal 1**: Improve CDR-H3 loop conformation | CDR RMSD, especially H3 |
| **Sub-goal 2**: Correct global docking orientation | iRMSD, DockQ, fnat |
| **Sub-goal 3**: Focus docking on the correct epitope | Epitope precision/recall/F1 at 5Å and 8Å |

---

## Baselines

Three baselines are compared against Y+:

| Label | Folder | Description |
|-------|--------|-------------|
| **baseline** | `antigen_cut/` | Vanilla Boltz-2, no steering |
| **contact_restraints** | `antigen_cut_contact_restraints/` | Oracle contact restraints (upper bound) |
| **pocket_guided** | `antigen_cut_vhvl_msa_pocket_ab_boltz_post_2023_omm/` | Pocket-guided restraints |

The most important comparison is **Y+ vs baseline** (does it help?) and **Y+ vs contact_restraints** (how far from the oracle upper bound?). Y+ should approach contact_restraints without requiring any experimental data.

---

## Ground Truths

- **Location**: `examples/pdb_minimized/`
- **Format**: `{COMPLEX_NAME}.pdb` (e.g., `7TRH_HBG.pdb`)
- **Coverage**: ~48 complexes (crystal structures)
- **Chains**: A = antigen, B = heavy chain, C = light chain

---

## Prediction Folder Structure

```
predictions_examples/
├── antigen_cut/                                      [baseline]
│   └── boltz_results_{complex}/
│       └── predictions/{complex}/{complex}_model_{0-4}.pdb
│                                    confidence_{complex}_model_{0-4}.json
│                                    pae_{complex}_model_{0-4}.npz
│                                    plddt_{complex}_model_{0-4}.npz
├── antigen_cut_contact_restraints/                   [contact_restraints]
│   └── boltz_results_restraint_{complex}_{setting}/
│       └── predictions/{complex}/...
├── antigen_cut_vhvl_msa_pocket_ab_boltz_post_2023_omm/  [pocket_guided]
│   └── boltz_results_restraint_to_A_{complex}/
│       └── predictions/{complex}/...
└── hierarchical_steering_v2/                         [Y+ new feature]
    └── boltz_results_{complex}/
        └── predictions/{complex}/{complex}_model_{0-4}.pdb
```

---

## Model Selection: Which of the 5 Predictions to Use?

Each run produces 5 independent structure predictions (`model_0` to `model_4`). Three reporting strategies are used:

| Strategy | How | When to use |
|----------|-----|------------|
| **Top-1 by confidence** | Select the model with the highest `confidence_score` from the JSON file | Primary reported number — simulates real usage (no oracle knowledge) |
| **Oracle best** | Select the model with the highest DockQ against the ground truth | Upper bound — shows the best possible result per method |
| **Mean over all 5** | Average metrics across all models | Assesses variance and consistency |

**Recommendation**: Report **top-1 by confidence** as the primary metric (realistic scenario), and include oracle and mean as supplementary. This matches how Round 1 was evaluated and makes comparisons fair.

**Why not pLDDT or PAE for selection?** The `confidence_score` in Boltz-2's JSON already combines pLDDT and PAE into a single calibrated score. It is the canonical selection criterion used throughout this project.

---

## Structural Alignment (Preprocessing)

Raw Boltz-2 predictions are in arbitrary coordinate frames. Before computing RMSD-based metrics, we need to align them to the ground truth.

**Standard approach for antibody-antigen complexes:**

1. **Align on antibody framework** — Superimpose the predicted antibody (heavy + light chain, excluding CDR residues) onto the native antibody framework using Cα atoms.
2. **Evaluate CDR RMSD in this frame** — After aligning frameworks, the CDR loop positions are in the same reference frame as the native, so RMSD measures loop displacement directly.
3. **Evaluate interface metrics (DockQ) after global alignment** — DockQ internally handles alignment (it uses `lDDT-PLI` logic and interface RMSD, which involve superimposing receptor then measuring ligand displacement). The `dockq` command does this automatically.

**Why align on framework, not full chain?** CDR loops are what we are trying to improve. Aligning on the full antibody would mix framework quality with CDR quality. Aligning on the framework isolates whether the CDRs moved to the right positions.

**Implementation**: The existing `evaluate_hierarchical_steering.py` already does this correctly using BioPython's `Superimposer` on framework Cα atoms.

---

## Metrics

### Primary Metrics (reported in summary table)

| Metric | Tool | Description |
|--------|------|-------------|
| **DockQ (Ab-Ag)** | `dockq2` conda env | Standard protein docking quality score for the antibody-antigen interface. Range 0–1; >0.23 = Acceptable, >0.49 = Medium, >0.8 = High (CAPRI thresholds). Captures fnat, iRMSD, LRMSD simultaneously. |
| **CDR-H3 RMSD** | BioPython | Cα RMSD of heavy chain CDR3 loop after framework alignment. The most clinically relevant single metric: CDR-H3 mediates most binding specificity. |
| **CAPRI classification** | Derived from DockQ | Fraction of complexes in each CAPRI category: Incorrect (DockQ<0.23), Acceptable (0.23–0.49), Medium (0.49–0.80), High (>0.80). |

### Secondary Metrics

| Metric | Tool | Description |
|--------|------|-------------|
| **DockQ (total)** | `dockq2` | DockQ averaged across all chain-chain interfaces |
| **iRMSD** | DockQ output | Interface RMSD — deviation of interface residues after receptor superposition |
| **LRMSD** | DockQ output | Ligand RMSD — displacement of the antibody after superimposing antigen |
| **fnat** | DockQ output | Fraction of native contacts recovered |
| **CDR RMSD per loop** | BioPython | Individual RMSD for H1, H2, H3, L1, L2, L3 |
| **Epitope F1 @ 5Å** | BioPython | F1 score comparing predicted vs native epitope residues (contacts within 5Å) |
| **Epitope F1 @ 8Å** | BioPython | Same at 8Å threshold (Cα-Cα distance, appropriate for coarser epitope mapping) |
| **Confidence (iptm)** | From JSON | Inter-chain pLDDT — measures predicted quality of interfaces |

### Statistical Analysis

- **Paired Wilcoxon signed-rank test** (non-parametric): Y+ vs each baseline, per metric
- **Bootstrap 95% confidence intervals** on means
- Report wins/ties/losses counts per complex

The Wilcoxon test is preferred over t-test because DockQ scores are not normally distributed and the sample size (~48) is moderate.

---

## Conda Environments

| Task | Environment | Reason |
|------|-------------|--------|
| Structure evaluation (DockQ) | `conda activate dockq2` | Contains DockQ v2 Python API — the standard tool for protein interface quality assessment |
| CDR RMSD, epitope metrics, plotting | `conda activate boltz` | Contains BioPython, numpy, matplotlib, scipy — all needed for parsing PDBs, alignment, and stats |

The existing evaluation script imports from both: DockQ is called via subprocess or Python API in `dockq2`, while the rest of the analysis runs in `boltz`.

**Recommendation**: Run the full script under `conda activate dockq2`, since DockQ is the bottleneck dependency and BioPython/numpy/matplotlib are typically also installed there.

---

## Running the Evaluation

### Quick check (single complex, to verify script works)
```bash
conda activate dockq2
cd ../../worktree_steering_cdr3/y_hierarchical_steering

python scripts/eval/evaluate_hierarchical_steering.py \
    --ground_truth_dir /path/to/examples/pdb_minimized \
    --predictions_dir /path/to/examples/predictions_examples \
    --cdrs_csv /path/to/examples/cdrs.csv \
    --out_dir evaluation_results_yplus \
    --complexes 7TRH_HBG \
    --methods baseline hierarchical_v2
```

### Full evaluation (all 48 complexes, all 4 methods)
```bash
conda activate dockq2
cd ../../worktree_steering_cdr3/y_hierarchical_steering

python scripts/eval/evaluate_hierarchical_steering.py \
    --ground_truth_dir /path/to/examples/pdb_minimized \
    --predictions_dir /path/to/examples/predictions_examples \
    --cdrs_csv /path/to/examples/cdrs.csv \
    --out_dir evaluation_results_yplus
```

### Required modification before running
Update `METHOD_CONFIGS` in the script to add/update the `hierarchical_v2` entry:
```python
"hierarchical_v2": {
    "subfolder": "hierarchical_steering_v2",   # match the actual folder name used for predictions
    "prefix": "boltz_results_",
    "suffix": "",
    "multi_setting": False,
},
```

---

## What Success Looks Like

Based on Round 1 results (Strategy Y baseline: CDR-H3 RMSD −0.30 Å, +0.034 DockQ) and the Y+ upgrades, we expect:

| Metric | Round 1 Y | Target for Y+ | Stretch goal |
|--------|----------|--------------|--------------|
| CDR-H3 RMSD improvement | −0.30 Å | −0.40 Å | −0.60 Å |
| DockQ (Ab-Ag) improvement | +0.034 | +0.05 | +0.10 |
| CAPRI Medium+ rate | +6.4 pp | +10 pp | +15 pp |
| Epitope F1 @ 8Å | not measured in Y | > baseline | approaching oracle |

**Minimum bar for Y+ to be considered an improvement over Y**: At least two of the three primary metrics (CDR-H3 RMSD, DockQ, CAPRI rate) must show statistically significant improvement (Wilcoxon p < 0.05) over the baseline, with effect size larger than Y.

---

## Output Files Expected

The evaluation script produces:

```
evaluation_results_yplus/
├── per_complex_all_metrics.csv     # Full results, one row per complex per method
├── summary_table.csv               # Mean ± CI, median, CAPRI rates per method
├── statistical_tests.csv           # Wilcoxon tests: Y+ vs each baseline
└── figures/
    ├── dockq_barplot.png           # Mean DockQ per method with error bars
    ├── dockq_boxplot.png           # Distribution of DockQ scores
    ├── dockq_scatter_baseline_vs_hierarchical.png
    ├── cdr_rmsd_barplot.png        # CDR RMSD per loop per method
    ├── epitope_f1_barplot.png      # Epitope F1 @ 5Å and 8Å
    └── confidence_calibration.png  # Confidence vs DockQ scatter
```
