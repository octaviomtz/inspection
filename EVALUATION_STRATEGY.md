# Q: Iterative Epitope Refinement -- Evaluation Strategy

## 1. Feature Goals Recap

Method Q aims to **discover epitope locations on unknown antigens and generate high-quality antibody-antigen complexes** through a 3-round iterative approach:

- **Round 1 (Broad exploration):** Weak guidance (scale 0.1) over the full antigen surface with exploratory beta-noise. Goal: survey the entire antigen for potential binding sites.
- **Round 2 (Targeted refinement):** Medium guidance (scale 0.5) focused on hotspot residues (top 20% from Round 1). Goal: narrow down to true epitope regions.
- **Round 3 (Validation):** Full guidance (scale 1.0) on refined hotspots. Goal: produce high-confidence final structures.

The method has **two output types**:
1. **Predicted structures** (PDB files) -- assessed for docking quality
2. **Implicit epitope prediction** -- the antigen residues in contact with CDR atoms in the final structure(s) represent the predicted epitope

We therefore evaluate Q on **two axes**:
- **Axis A (Primary): Epitope prediction accuracy** -- can Q correctly identify which antigen residues are in the binding interface?
- **Axis B (Secondary): Docking / interface quality** -- does Q maintain or improve the structural prediction quality compared to baselines?

---

## 2. Baselines

| Label | Folder | Description |
|-------|--------|-------------|
| **B1** | `antigen_cut` | Vanilla Boltz2, no constraints or steering |
| **B2** | `antigen_cut_contact_restraints` | Contact restraints (hbond, hydrophobic, salt-bridge) -- multiple settings per complex |
| **B3** | `antigen_cut_vhvl_msa_pocket_ab_boltz_post_2023_omm` | Pocket restraints to specific residues -- multiple settings per complex |
| **Q** | `new_feature_q_epitope_refinement` (or equivalent) | Our method: iterative epitope refinement |

For B2 and B3, which have multiple restraint configurations per complex:
- Report the **best** result across configurations (oracle upper bound)
- Report the **mean** across configurations (practical expected performance)

---

## 3. Ground Truth

- **Crystal structures:** `pdb_minimized/` (48 complexes, PDB format)
- **Chain convention:** A = antigen, B = heavy chain (VH), C = light chain (VL)
- **CDR indices:** `examples/cdrs.csv` (per-complex CDR1/2/3 residue indices for heavy and light chains)
- **Ground-truth epitope definition:** Antigen residues with any heavy atom within a distance cutoff of any antibody heavy atom. We use **5 A as primary** (strict) and **8 A as secondary** (permissive).

```python
gt_epitope = {res_idx for res in antigen_chain
              if min_heavy_atom_dist(res, antibody_atoms) < 5.0}
```

---

## 4. Model Selection: Which of the 5 Predictions to Evaluate?

Each boltz run produces 5 models (`model_0` to `model_4`), each with a confidence JSON, pLDDT, PAE, and PDE files.

| Strategy | Description | Use |
|----------|-------------|-----|
| **Best-by-confidence** | Model with the highest `confidence_score` from its JSON file | **Primary**. This mimics the user workflow (boltz ranks by confidence). |
| **Best-by-DockQ** | Model with the highest DockQ against ground truth | Oracle upper bound. Shows the best the method *could* achieve if confidence ranking were perfect. |

For each complex, report both. Analysis will compare whether Q improves confidence-DockQ correlation (i.e., does the best-by-confidence model increasingly match the best-by-DockQ model?).

---

## 5. Axis A: Epitope Prediction Accuracy (Primary)

### 5.1 Rationale

Method Q's core goal is epitope discovery. The 3-round iterative process should progressively narrow contacts toward the true epitope. This axis measures whether the final predicted structures place CDR atoms near the correct antigen residues.

### 5.2 Ground-Truth Epitope Extraction

For each complex, extract the set of antigen residue indices that form the interface in the crystal structure:

```python
from Bio.PDB import PDBParser, NeighborSearch

parser = PDBParser(QUIET=True)
structure = parser.get_structure('ref', native_pdb)
model = structure[0]

# Collect all antibody heavy atoms (chains B + C)
ab_atoms = [atom for chain_id in ['B', 'C']
            for atom in model[chain_id].get_atoms()]
ns = NeighborSearch(ab_atoms)

# Epitope = antigen residues with any atom within 5 A of any antibody atom
gt_epitope = set()
for residue in model['A'].get_residues():
    if residue.id[0] != ' ':
        continue
    for atom in residue.get_atoms():
        if ns.search(atom.coord, 5.0):
            gt_epitope.add(residue.id[1])
            break
```

### 5.3 Predicted Epitope Extraction

For Q and all baselines, extract the predicted epitope from the best-by-confidence predicted structure using the same contact criterion:

```python
# Same logic applied to predicted PDB
pred_epitope = set()
for residue in pred_model['A'].get_residues():
    if residue.id[0] != ' ':
        continue
    for atom in residue.get_atoms():
        if ns_pred.search(atom.coord, 5.0):
            pred_epitope.add(residue.id[1])
            break
```

**Note on residue indexing:** The ground truth PDBs have continuous numbering across chains (A: 1-212, B: 213-434, C: 435-648) while predicted PDBs number each chain independently (A: 1-214, B: 1-126, C: 1-110). The evaluation script must extract epitope residues using within-chain iteration (by chain object, not by raw residue number), so indexing is handled implicitly by BioPython's chain/residue objects.

### 5.4 Metrics

| Metric | Formula | Description |
|--------|---------|-------------|
| **Precision** | TP / (TP + FP) | What fraction of predicted epitope residues are true contacts? |
| **Recall** | TP / (TP + FN) | What fraction of true epitope residues are predicted? |
| **F1** | 2 * P * R / (P + R) | Harmonic mean of precision and recall |
| **MCC** | Matthews Correlation Coefficient | Balanced metric robust to class imbalance |

Where TP/FP/FN are defined by comparing the predicted epitope residue set against the ground-truth epitope residue set. The matching is done by residue position within the antigen chain (1-indexed).

### 5.5 Reporting

For each complex x method, report: precision, recall, F1, MCC.
Aggregate across all 48 complexes: **median** and **mean +/- std** for each metric.

---

## 6. Axis B: Docking / Interface Quality (Secondary)

### 6.1 DockQ v2 (Primary Interface Metric)

DockQ is the standard metric for protein-protein docking quality. Run via the `dockq2` conda environment.

```bash
conda run -n dockq2 DockQ <model.pdb> <native.pdb> --json <out.json>
```

DockQ combines:
- **fnat**: Fraction of native contacts recovered
- **iRMSD**: Interface RMSD (backbone atoms at the interface, after superposition)
- **LRMSD**: Ligand RMSD (after superposition on the receptor)
- **F1**: Harmonic mean of fnat and (1 - fnonnat)

DockQ quality thresholds (CAPRI classification):
- **Incorrect:** DockQ < 0.23
- **Acceptable:** 0.23 <= DockQ < 0.49
- **Medium:** 0.49 <= DockQ < 0.80
- **High:** DockQ >= 0.80

### 6.2 Which Interfaces to Report

Each complex has 3 chain-pair interfaces. For Q evaluation:

| Interface | Chains | Purpose |
|-----------|--------|---------|
| **AB** (antigen-heavy) | A-B | Primary: this is the main binding interface Q aims to improve |
| **AC** (antigen-light) | A-C | Primary: second antibody-antigen interface |
| **BC** (heavy-light) | B-C | Sanity check: should not degrade (internal antibody interface) |
| **AB+AC average** | A-(B,C) | Summary of antigen-antibody interface quality |
| **Global DockQ** | all 3 | Overall quality |

### 6.3 Additional Structural Metrics from Boltz Confidence

From the per-model confidence JSON files, extract:

| Metric | Key in JSON | Description |
|--------|-------------|-------------|
| **confidence_score** | `confidence_score` | Boltz composite score (4*plddt + iptm)/5 |
| **iptm** | `iptm` | Interface pTM -- measures inter-chain prediction quality |
| **complex_plddt** | `complex_plddt` | Average pLDDT across all residues |
| **complex_iplddt** | `complex_iplddt` | pLDDT restricted to interface residues |
| **pair_iptm AB** | `pair_chains_iptm["0"]["1"]` | Chain A-B iptm (antigen-heavy) |
| **pair_iptm AC** | `pair_chains_iptm["0"]["2"]` | Chain A-C iptm (antigen-light) |

### 6.4 CDR3-Specific Metrics

Since Q uses CDR3 regions for steering, evaluate the quality of CDR3 predictions specifically:

- **pLDDT_cdr3h**: Average pLDDT over heavy chain CDR3 residues
- **pLDDT_cdr3l**: Average pLDDT over light chain CDR3 residues
- **CDR3-antigen PAE**: Average PAE between CDR3 residues and antigen residues (from PAE .npz files)

CDR3 residue indices come from `examples/cdrs.csv`.

---

## 7. Preprocessing & Alignment

### 7.1 Chain Mapping

DockQ handles alignment internally. However, the ground truth and predicted PDBs may have different residue numbering. DockQ solves this by sequence alignment, so no manual renumbering is needed. Use the automatic chain mapping:

```bash
DockQ model.pdb native.pdb  # auto-maps ABC:ABC
```

If chain mapping fails (rare), specify explicitly:
```bash
DockQ model.pdb native.pdb --mapping ABC:ABC
```

### 7.2 LRMSD Computation (Antibody-Aligned RMSD)

For additional insight beyond DockQ, compute **antibody-aligned antigen RMSD**. This answers: "if we perfectly align the antibody, how far off is the antigen placement?"

```python
from Bio.PDB import Superimposer

# 1. Extract CA atoms from antibody (chains B+C) in both structures
# 2. Superimpose predicted antibody onto native antibody
# 3. Apply same transform to predicted antigen
# 4. Compute RMSD of antigen CA atoms
```

This is conceptually similar to DockQ's LRMSD but gives a single number for the antigen as a whole.

### 7.3 No Additional Preprocessing Required

- DockQ handles structural superposition internally
- BioPython handles PDB parsing for both predicted and native
- No energy minimization or relaxation of predictions is needed (ground truths are already minimized in `pdb_minimized/`)

---

## 8. Evaluation Pipeline

### 8.1 Per-Complex Workflow

For each complex (e.g., `7TRH_HBG`) and each method (B1, B2, B3, Q):

```
1. Locate all prediction folders for this complex+method
2. For each prediction folder:
   a. Read confidence JSONs for models 0-4
   b. Select best-by-confidence model
   c. Run DockQ(best_model.pdb, native.pdb) -> DockQ scores per interface
   d. Extract predicted epitope from best model (5 A threshold)
   e. Extract ground-truth epitope from native (5 A threshold)
   f. Compute epitope precision, recall, F1, MCC
   g. Extract CDR3-specific pLDDT and PAE metrics
3. For B2/B3: aggregate across restraint configurations (best, mean)
```

### 8.2 Aggregation

```
For each method:
  - Compute median, mean, std across all 48 complexes
  - Count CAPRI classifications (Incorrect / Acceptable / Medium / High)
  - Compute statistical tests (Wilcoxon signed-rank, paired) vs B1 baseline
```

### 8.3 Summary Tables

**Table 1: Docking Quality (DockQ)**

| Method | DockQ (AB+AC) | DockQ (global) | iRMSD (AB) | fnat (AB) | CAPRI High% |
|--------|--------------|----------------|------------|-----------|-------------|
| B1 | ... | ... | ... | ... | ... |
| B2 (best) | ... | ... | ... | ... | ... |
| B3 (best) | ... | ... | ... | ... | ... |
| Q | ... | ... | ... | ... | ... |

**Table 2: Epitope Prediction Accuracy**

| Method | Precision | Recall | F1 | MCC |
|--------|-----------|--------|----|-----|
| B1 | ... | ... | ... | ... |
| B2 (best) | ... | ... | ... | ... |
| B3 (best) | ... | ... | ... | ... |
| Q | ... | ... | ... | ... |

**Table 3: Confidence & Per-Residue Metrics**

| Method | confidence_score | iptm | complex_plddt | pLDDT_cdr3h | CDR3-ag PAE |
|--------|-----------------|------|---------------|-------------|-------------|
| B1 | ... | ... | ... | ... | ... |
| Q | ... | ... | ... | ... | ... |

### 8.4 Plots

1. **DockQ box plot**: Per-method DockQ distribution across 48 complexes (AB+AC interface)
2. **Epitope F1 box plot**: Per-method F1 distribution
3. **Scatter: DockQ vs confidence_score**: One point per complex, colored by method. Shows confidence calibration.
4. **Per-complex delta plot**: DockQ(Q) - DockQ(B1) for each complex, sorted. Shows which complexes benefit most.
5. **CAPRI classification bar chart**: Stacked bars showing fraction of Incorrect / Acceptable / Medium / High per method.

---

## 9. Statistical Tests

- **Wilcoxon signed-rank test** (paired, non-parametric): Compare Q vs B1 on matched complexes for each metric. Reports p-value and effect size (rank-biserial correlation).
- **Bootstrap 95% confidence intervals**: For median DockQ and F1 differences.
- **Significance threshold**: p < 0.05 (with Bonferroni correction if testing multiple metrics).

---

## 10. Environment & Dependencies

### 10.1 Primary Environment: `dockq2`

```bash
conda activate dockq2
```

Contains:
- **DockQ** (v2): `DockQ model.pdb native.pdb --json out.json`
- **BioPython**: PDB parsing, NeighborSearch, Superimposer
- **NumPy, SciPy, pandas**: Numerical computation
- Standard Python 3.11

This environment has everything needed for the evaluation. No additional environment is required.

### 10.2 Packages Used

| Package | Use | Available in dockq2? |
|---------|-----|---------------------|
| `DockQ` | Interface quality scoring | Yes |
| `Bio.PDB` (BioPython) | PDB parsing, contact extraction, superposition | Yes |
| `numpy` | Array operations, statistics | Yes |
| `scipy.stats` | Wilcoxon test, bootstrap | Yes |
| `pandas` | DataFrames, CSV I/O, aggregation | Yes |
| `matplotlib` / `seaborn` | Plots | Need to verify, install if missing |

### 10.3 Verification

```bash
conda activate dockq2
python -c "from DockQ.DockQ import load_PDB, run_on_all_native_interfaces; print('DockQ OK')"
python -c "from Bio.PDB import PDBParser, NeighborSearch; print('BioPython OK')"
python -c "import pandas, scipy, numpy; print('Core deps OK')"
```

---

## 11. Folder Structure Reference

### 11.1 Ground Truth
```
pdb_minimized/
  7TRH_HBG.pdb          # Chain A=antigen, B=heavy, C=light
  7TRI_ZYB.pdb           # Residue numbering: continuous across chains
  ...                    # (48 PDB files)
```

### 11.2 Predictions (per method)
```
predictions_examples/
  antigen_cut/                                    # B1
    boltz_results_7TRH_HBG/
      predictions/7TRH_HBG/
        7TRH_HBG_model_0.pdb                     # 5 structure predictions
        ...
        7TRH_HBG_model_4.pdb
        confidence_7TRH_HBG_model_0.json          # 5 confidence files
        ...
        pae_7TRH_HBG_model_0.npz                  # 5 PAE matrices
        pde_7TRH_HBG_model_0.npz                  # 5 PDE matrices
  antigen_cut_contact_restraints/                 # B2
    boltz_results_restraint_7TRH_HBG_hbond_23/    # Multiple configs per complex
    boltz_results_restraint_7TRH_HBG_hbond_7/
    ...
  antigen_cut_vhvl_msa_pocket_.../                # B3
    boltz_results_restraint_to_A_7TRH_HBG_.../    # Multiple configs per complex
    ...
  new_feature_q_epitope_refinement/               # Q
    boltz_results_7TRH_HBG_.../
      predictions/7TRH_HBG/
        ...
```

### 11.3 CDR Definitions
```
examples/cdrs.csv     # complex, heavy_seq, light_seq, cdr1_h, cdr2_h, cdr3_h, cdr1_l, cdr2_l, cdr3_l
                      # CDR indices are 1-based residue positions within each chain
```

---

## 12. Implementation Plan

The evaluation script should be a single Python file (`evaluation/evaluate_q.py`) that:

1. **Discovers** all prediction folders for each method and maps them to complexes
2. **Selects** the best model per complex (by confidence and by DockQ oracle)
3. **Computes** DockQ via the DockQ library (not subprocess) for all interfaces
4. **Extracts** epitope residue sets from predicted and native structures
5. **Computes** epitope accuracy metrics (precision, recall, F1, MCC)
6. **Extracts** confidence metrics from JSON files
7. **Extracts** CDR3-specific pLDDT and CDR3-antigen PAE from npz files
8. **Aggregates** results across complexes with statistics and significance tests
9. **Outputs** summary CSV tables and plots

The script should be modular, reusing patterns from the existing `g_phase_refinement/evaluation/evaluate.py` (METHOD_DEFS dict, extract_complex_name, confidence extraction, DockQ via Python API).

```bash
# Usage
conda activate dockq2
python evaluation/evaluate_q.py \
  --pred_dir /path/to/predictions_examples \
  --native_dir /path/to/pdb_minimized \
  --cdr_csv examples/cdrs.csv \
  --out_dir evaluation/results_q
```
