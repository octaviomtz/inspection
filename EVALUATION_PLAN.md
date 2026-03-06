# Evaluation Plan: Strategy D+ (Enhanced Canonical Ensemble)

## 1. Goal of the Feature

Strategy D+ generates diverse CDR conformational predictions by sweeping β-scaling values on pair representations (`z[CDR] *= (1 + β)`) during diffusion conditioning. By sweeping β from -0.5 to +0.5, it produces an ensemble of structures and selects the top-K by confidence score.

**What it aims to improve:**
- Diversity of CDR conformations sampled (especially CDR-H3)
- Antibody-antigen interface quality (better contacts between CDR loops and antigen)
- Computational efficiency (trunk runs once, only diffusion/confidence per β)

**What it does NOT change:**
- The antigen fold itself
- The antibody framework regions
- The VH-VL pairing geometry (only CDR loops are affected)

## 2. Experimental Setup

### 2.1 Methods Being Compared

| Label | Method | Folder Pattern | Description |
|-------|--------|----------------|-------------|
| **Baseline 1** | Boltz2 vanilla | `antigen_cut/boltz_results_{complex}/` | Standard Boltz2 prediction, no constraints |
| **Baseline 2** | Contact restraints | `antigen_cut_contact_restraints/boltz_results_restraint_{complex}_{type}_{id}/` | Boltz2 with specific contact restraints (hbond, hydrophobic, salt_bridge). Multiple runs per complex. |
| **Baseline 3** | Pocket + MSA | `antigen_cut_vhvl_msa_pocket_ab_boltz_post_2023_omm/boltz_results_restraint_to_A_{complex}_{chain}_{res}_{id}/` | Boltz2 with pocket constraints and VH/VL MSA. Multiple runs per complex. |
| **D+ (Ours)** | Canonical ensemble | `new_feature_cdr3_beta/boltz_results_{complex}_cdr3_beta/` | β-sweep canonical ensemble, top-5 by confidence |

### 2.2 Data

- **Test set**: 48 antibody-antigen complexes from `examples/cdrs.csv`
- **Ground truth**: Crystal structures in `pdb_minimized/` (48 PDB files, OpenMM-minimized)
- **Structure format**: Ground truth = PDB (all-atom, hydrogens), Predictions = PDB (Baseline 1) or CIF (D+), heavy atoms only
- **Chain convention**: A = antigen, B = heavy chain, C = light chain (all methods)

### 2.3 Model Selection (Which of the 5 predictions to evaluate?)

Each method produces up to 5 models (model_0 through model_4), already ranked by confidence score (model_0 = best).

**Approach:** Evaluate using **two protocols**:
1. **Best-model (model_0):** Report metrics for the top-ranked prediction only. This measures whether the method's confidence ranking is effective.
2. **Oracle-best:** Report the best metric across all 5 models. This measures the method's ability to generate at least one good structure, regardless of ranking.

For baselines 2 and 3 (which have multiple restraint runs per complex), apply the same logic: pick the best model_0 across all restraint variants (best-of-best), and separately the oracle-best across all variants and models.

## 3. Evaluation Metrics

### 3.1 Primary Metric: DockQ (Interface Quality)

**Tool:** DockQ v2 (`conda activate dockq2`)

DockQ is the standard CAPRI-style metric for protein-protein interface quality. It combines:
- **fnat**: Fraction of native contacts recovered (contact < 5Å)
- **iRMSD**: Interface RMSD (residues within 10Å of interface, heavy atoms)
- **LRMSD**: Ligand RMSD (after receptor alignment)

**Quality thresholds (CAPRI):**
| DockQ Range | Classification |
|------------|---------------|
| < 0.23     | Incorrect     |
| 0.23–0.49  | Acceptable    |
| 0.49–0.80  | Medium        |
| ≥ 0.80     | High          |

**Which interfaces to report:**
- **AB (antigen–heavy chain):** The most relevant for CDR-H steering
- **AC (antigen–light chain):** Relevant for CDR-L steering
- **BC (heavy–light VH-VL):** Should remain stable (D+ shouldn't hurt this)
- **Global DockQ:** Average across all interfaces

**Command:**
```bash
conda run -n dockq2 DockQ <model.pdb|model.cif> <native.pdb> --json <output.json>
```

DockQ v2 handles:
- PDB and CIF input formats
- Automatic chain mapping (via sequence alignment)
- Mismatched chain lengths (GT has full Fab constant domains; predictions have Fv only)

### 3.2 Secondary Metric: CDR-Specific RMSD

**Tool:** BioPython (`conda activate dockq2`, has biopython 1.86)

This measures whether the CDR loop conformations specifically improved. The approach:

1. **Parse** both model and native structures (handling the residue numbering mismatch — GT uses continuous numbering, predictions use per-chain)
2. **Sequence-align** model vs native chains to identify corresponding residues (needed because GT chains are longer — full Fab vs Fv)
3. **Superimpose on framework residues** of chains B and C (all non-CDR residues as defined in `cdrs.csv`)
4. **Compute RMSD on CDR residues only** (Cα atoms), separately for:
   - CDR-H1, CDR-H2, CDR-H3 (chain B)
   - CDR-L1, CDR-L2, CDR-L3 (chain C)
   - All-CDR combined

CDR-H3 RMSD is the most informative single metric because H3 is the most variable loop and the primary target of β-scaling.

### 3.3 Tertiary Metric: Predicted Confidence vs Actual Quality

**Tool:** Python (numpy, pandas) in `conda activate dockq2`

Assess whether the model's self-reported confidence correlates with actual quality:

- **confidence_score** (from confidence JSON) vs **DockQ** (computed)
- **complex_plddt** vs DockQ
- **iptm** (inter-chain predicted TM-score) vs DockQ

Report: Spearman rank correlation between predicted confidence and actual DockQ, across all models of all complexes. This tells us whether the confidence-based ranking in the ensemble is meaningful.

### 3.4 Ensemble Diversity Metric (D+ Specific)

**Tool:** BioPython (`conda activate dockq2`)

Measures whether β-sweeping actually produces diverse CDR conformations:

- **Pairwise CDR RMSD** across the top-5 models (after framework alignment)
- Report: mean, std, and range of pairwise CDR-H3 RMSD within the ensemble
- Compare against Baseline 1's 5-model diversity (which comes from diffusion stochasticity alone)

## 4. Postprocessing Pipeline

### 4.1 Structure Preparation

No explicit alignment needed before DockQ — DockQ v2 performs internal sequence alignment and handles chain mapping automatically. Verified working with both PDB and CIF formats against the GT PDBs.

For CDR-specific RMSD, the pipeline is:

```
For each (complex, method):
  1. Load model structure (PDB or CIF) — BioPython PDBParser or MMCIFParser
  2. Load native structure (PDB) — BioPython PDBParser
  3. Sequence-align model chain B to native chain B (Fv→Fab matching)
     → Identify corresponding residue pairs
  4. Classify residues as framework or CDR using cdrs.csv indices
  5. Extract framework Cα atoms (matched pairs only)
  6. Superimpose model onto native using framework Cα (Superimposer)
  7. Compute RMSD on CDR Cα residues
```

### 4.2 Handling Multiple Runs (Baselines 2 & 3)

Baselines 2 and 3 have multiple restraint variants per complex:
- Baseline 2: `restraint_{complex}_hbond_{id}`, `restraint_{complex}_hydrophobic_{id}`, etc.
- Baseline 3: `restraint_to_A_{complex}_{chain}_{res}_{id}`

**Aggregation strategy:**
- For each complex, collect all variant results
- Report **best model_0** (highest DockQ across variants) — simulates choosing the best restraint
- Report **oracle** (highest DockQ across all variants × all models)
- Report **median** across all variant model_0s — measures average restraint quality

### 4.3 File Discovery Logic

```python
# Baseline 1: one folder per complex
f"antigen_cut/boltz_results_{complex}/predictions/{complex}/{complex}_model_{i}.pdb"

# Baseline 2: multiple folders per complex
f"antigen_cut_contact_restraints/boltz_results_restraint_{complex}_*/predictions/restraint_{complex}_*/*_model_{i}.pdb"

# Baseline 3: multiple folders per complex
f"antigen_cut_vhvl_msa_pocket_.../boltz_results_restraint_to_A_{complex}_*/predictions/restraint_to_A_{complex}_*/*_model_{i}.pdb"

# D+ (ours): one folder per complex
f"new_feature_cdr3_beta/boltz_results_{complex}_cdr3_beta/predictions/{complex}_cdr3_beta/{complex}_cdr3_beta_model_{i}.cif"
```

Note: The exact folder name for D+ results may vary on the evaluation machine. The script should use glob patterns to discover prediction folders.

## 5. Summary Tables to Produce

### Table 1: DockQ Per-Interface (Best Model)

| Complex | Method | DockQ_AB | DockQ_AC | DockQ_BC | DockQ_Global | fnat_AB | iRMSD_AB | LRMSD_AB |
|---------|--------|----------|----------|----------|--------------|---------|----------|----------|
| 7TRH_HBG | Baseline 1 | ... | ... | ... | ... | ... | ... | ... |
| 7TRH_HBG | D+ (ours)  | ... | ... | ... | ... | ... | ... | ... |
| ...     | ...    | ...      | ...      | ...      | ...          | ...     | ...      | ...      |

### Table 2: CDR RMSD (Best Model, Framework-Aligned)

| Complex | Method | H1 | H2 | H3 | L1 | L2 | L3 | All-CDR |
|---------|--------|----|----|----|----|----|----| --------|
| ...     | ...    | ...| ...| ...| ...| ...| ...|  ...    |

### Table 3: Aggregated Results (Mean ± Std Across 48 Complexes)

| Method | DockQ_AB | DockQ_AC | DockQ_Global | CDR-H3 RMSD | All-CDR RMSD | Confidence Correlation |
|--------|----------|----------|--------------|-------------|--------------|----------------------|
| Baseline 1 | ... | ... | ... | ... | ... | ... |
| Baseline 2 (best) | ... | ... | ... | ... | ... | ... |
| Baseline 3 (best) | ... | ... | ... | ... | ... | ... |
| D+ (ours) | ... | ... | ... | ... | ... | ... |

### Table 4: CAPRI Quality Classification (% of Complexes)

| Method | Incorrect | Acceptable | Medium | High |
|--------|-----------|------------|--------|------|
| Baseline 1 | ... | ... | ... | ... |
| D+ (ours)  | ... | ... | ... | ... |

### Table 5: Ensemble Diversity (D+ Specific)

| Complex | Mean Pairwise CDR-H3 RMSD | Range | Baseline 1 Diversity |
|---------|--------------------------|-------|---------------------|
| ...     | ...                      | ...   | ...                 |

## 6. Environment and Dependencies

| Task | Conda Environment | Key Packages |
|------|-------------------|--------------|
| DockQ computation | `dockq2` | DockQ 2.1.3, biopython 1.86, numpy, pandas, scipy |
| CDR RMSD computation | `dockq2` | biopython 1.86, numpy |
| Confidence analysis | `dockq2` | numpy, pandas, scipy (for Spearman correlation) |
| Plotting / reporting | `dockq2` | pandas, numpy (matplotlib if needed — check availability) |

All evaluation can be done in the `dockq2` environment. No new environment needed.

Note: If matplotlib is not installed in `dockq2`, install it with `conda run -n dockq2 pip install matplotlib` for generating plots.

## 7. Evaluation Script Structure

```
scripts/
  evaluate_canonical_ensemble.py   # Main evaluation script
```

**High-level workflow:**
```python
# 1. Discover all prediction files across all methods
# 2. For each (complex, method, model_idx):
#      a. Run DockQ → JSON output
#      b. Compute CDR-specific RMSD (framework-aligned)
#      c. Read confidence JSON
# 3. Aggregate into DataFrames
# 4. Produce summary tables (Tables 1–5)
# 5. Statistical tests: paired Wilcoxon signed-rank (D+ vs each baseline)
```

## 8. Key Considerations

1. **Chain length mismatch**: GT PDBs contain full Fab (constant + variable domains), while predictions contain only Fv (variable domains). DockQ handles this via internal sequence alignment. For custom RMSD calculations, we must align sequences first to find matching residues.

2. **File format**: Baseline 1 outputs PDB; D+ outputs CIF (mmCIF). Both are supported by DockQ and BioPython.

3. **Residue numbering**: GT uses continuous numbering across chains (B starts at 213, C at 435). Predictions use per-chain numbering starting at 1. DockQ handles this automatically; custom code must handle it explicitly.

4. **Statistical significance**: With 48 complexes, use non-parametric paired tests (Wilcoxon signed-rank) to compare D+ against each baseline. Report p-values alongside mean differences.

5. **Multiple restraint variants**: For baselines 2 and 3, the "best restraint" selection gives these methods an oracle advantage. To be fair, compare D+ against both the best-restraint and median-restraint aggregations.
