# K+ Evaluation Strategy: Region-Specific Beta-Scaling

## 1. Feature Goals Recap

K+ aims to **discover epitope locations on unknown antigens** by emphasizing different antigen surface regions during diffusion and accumulating contact heatmaps. It is NOT primarily a structure-quality improvement method; its main output is an **epitope propensity map** over antigen residues. However, the beta-scaling also modifies the diffusion trajectory, so it may also affect antibody-antigen docking quality.

We therefore evaluate K+ on **two axes**:
- **Axis A (primary):** Epitope prediction accuracy (can K+ correctly identify which antigen residues are in the interface?)
- **Axis B (secondary):** Docking / interface quality (does beta-scaling maintain or improve the structural prediction compared to baselines?)

---

## 2. Baselines

| Label | Folder | Description |
|-------|--------|-------------|
| **B1** | `antigen_cut` | Vanilla Boltz2, no constraints |
| **B2** | `antigen_cut_contact_restraints` | Contact restraints (hbond, hydrophobic, salt-bridge) - multiple settings per complex |
| **B3** | `antigen_cut_vhvl_msa_pocket_ab_boltz_post_2023_omm` | Pocket restraints to specific residues - multiple settings per complex |
| **K+** | `new_feature_cdr3_beta` (or the K+ output folder) | Our new feature |

For B2 and B3 there are multiple restraint configurations per complex. For the evaluation we should:
- Report the **best** result across configurations (oracle upper bound)
- Report the **mean** across configurations (practical expected performance)

---

## 3. Ground Truth

- **Crystal structures:** `pdb_minimized/` (48 complexes, PDB format)
- **Chain convention:** A = antigen, B = heavy chain (VH), C = light chain (VL)
- **Ground-truth epitope:** Defined as antigen residues with any heavy atom within a distance cutoff of any antibody heavy atom. Standard cutoffs: **5 A** (strict) and **8 A** (permissive). We use **5 A** as primary and report 8 A as secondary.

---

## 4. Model Selection: Which of the 5 Predictions to Evaluate?

Each run produces up to 5 models (`model_0` to `model_4`). We define two selection strategies:

| Strategy | Description | When to use |
|----------|-------------|-------------|
| **Best-by-confidence** | Select the model with the highest `confidence_score` from the JSON file | Primary. This mimics the user's workflow (boltz already ranks by confidence). |
| **Best-by-DockQ** | Select the model with the highest DockQ against ground truth | Oracle upper bound. Shows the best the method *could* achieve. |

For each complex we report:
- The **best-by-confidence** result (this is what the user would see)
- Optionally the **best-by-DockQ** result (as an upper bound)

---

## 5. Axis A: Epitope Prediction Accuracy

### 5.1 Ground-Truth Epitope Extraction

For each complex, extract the set of antigen residue indices that are within 5 A of any antibody atom in the crystal structure:

```
gt_epitope = {res_idx for res in antigen_chain
              if min_dist(res, antibody_atoms) < 5.0}
```

### 5.2 Predicted Epitope from K+

K+ produces per-region predictions. For each region, we extract antigen-antibody contacts from the predicted structure. Contacts are accumulated across all regions into a **contact heatmap** (vector of length = number of antigen residues). This heatmap can be thresholded to produce a binary predicted epitope set.

For baselines (B1, B2, B3) which do NOT produce epitope heatmaps, we extract the predicted epitope directly from each predicted structure using the same contact criterion as the ground truth (antigen residues within 5 A of antibody atoms). This gives a binary epitope set per prediction.

### 5.3 Metrics

| Metric | Formula | Description |
|--------|---------|-------------|
| **Precision** | TP / (TP + FP) | What fraction of predicted epitope residues are true contacts? |
| **Recall** | TP / (TP + FN) | What fraction of true epitope residues are predicted? |
| **F1** | 2 * P * R / (P + R) | Harmonic mean of precision and recall |
| **MCC** | Matthews correlation coefficient | Balanced metric for imbalanced binary classification |
| **AUC-ROC** | Area under ROC curve | Only for K+ (which produces continuous heatmap scores) |
| **AUC-PR** | Area under precision-recall curve | Better than AUC-ROC for imbalanced data (epitope is ~15% of surface) |

**Targets** (from design docs): Precision > 0.7, Recall > 0.6.

### 5.4 Implementation

```python
# For each complex:
# 1. Extract ground-truth epitope from crystal structure
# 2. For K+: get heatmap, threshold at multiple cutoffs -> PR curve
# 3. For baselines: extract contacts from best predicted structure
# 4. Compute precision, recall, F1, MCC
```

**Contact definition for extracting epitope from predicted structures:**
- Use CA-CA distance < 8 A between antigen and antibody residues (coarser but robust to side-chain conformational noise)
- OR use heavy-atom distance < 5 A (finer but more sensitive to packing errors)
- We use **heavy-atom < 5 A** as primary, report CA < 8 A as secondary.

---

## 6. Axis B: Docking / Interface Quality

### 6.1 DockQ v2 (Primary Interface Metric)

DockQ is the standard metric for protein-protein docking quality. It combines:
- **iRMSD**: Interface RMSD (backbone atoms at the interface)
- **LRMSD**: Ligand RMSD (after superposition on the receptor)
- **fnat**: Fraction of native contacts recovered
- **F1**: Harmonic mean of fnat and (1 - fnonnat)

DockQ thresholds:
- Incorrect: DockQ < 0.23
- Acceptable: 0.23 <= DockQ < 0.49
- Medium: 0.49 <= DockQ < 0.80
- High: DockQ >= 0.80

**Which interfaces to report:**
The complexes have 3 interfaces: AB (antigen-heavy), AC (antigen-light), BC (heavy-light). For K+ evaluation:
- **AB + AC average** (antigen-antibody interfaces): Primary. This is what K+ aims to improve.
- **BC** (heavy-light interface): Sanity check. Should not degrade.
- **Global DockQ** (average of all 3): Overall quality summary.

### 6.2 Additional Structural Metrics

| Metric | Description | Tool |
|--------|-------------|------|
| **iRMSD** | Interface RMSD (from DockQ) | DockQ |
| **LRMSD** | Ligand RMSD after receptor alignment (from DockQ) | DockQ |
| **fnat** | Fraction of native contacts (from DockQ) | DockQ |
| **Complex pLDDT** | Boltz's own confidence | confidence JSON |
| **Complex ipLDDT** | Interface pLDDT | confidence JSON |
| **iptm** | Interface predicted TM-score | confidence JSON |
| **Antibody-aligned antigen RMSD** | Align on antibody, measure antigen RMSD | BioPython Superimposer |

### 6.3 Antibody-Aligned Antigen RMSD (Custom Metric)

A common approach for antibody-antigen docking evaluation:
1. Superimpose predicted structure onto ground truth using **antibody backbone atoms only** (chains B+C CA atoms)
2. Compute **RMSD of antigen CA atoms** after this alignment

This directly measures how well the antigen is placed relative to the antibody.

```python
# Pseudocode:
# 1. Extract CA atoms for chains B, C from prediction and ground truth
# 2. Superimpose prediction onto GT using antibody CAs
# 3. Apply same transform to antigen CAs
# 4. Compute RMSD between transformed predicted antigen CAs and GT antigen CAs
```

---

## 7. Postprocessing Requirements

### 7.1 Chain Matching

- Ground truth and predictions both use chain IDs: A (antigen), B (heavy), C (light)
- The ground truth PDB may contain full antibody; predictions only contain variable domains (VH/VL). DockQ handles sequence alignment internally, so length mismatches are tolerated.
- Verify chain mapping is consistent across all complexes.

### 7.2 File Format Handling

- Baselines produce `.pdb` files
- New feature (K+) may produce `.cif` files
- DockQ accepts both PDB and CIF formats (verified)
- BioPython's `MMCIFParser` handles CIF files

### 7.3 Structure Alignment

- **DockQ performs its own alignment** internally - no pre-alignment needed for DockQ.
- For antibody-aligned antigen RMSD: use BioPython's `Superimposer` on antibody CA atoms.
- For epitope extraction from predictions: no alignment needed (contacts are internal to the structure).

---

## 8. Evaluation Pipeline Design

### Step 1: Data Discovery
```
For each complex_name in pdb_minimized/*.pdb:
    For each method in [B1, B2, B3, K+]:
        Find all result folders matching complex_name
        For each result folder:
            Find prediction files (*.pdb or *.cif)
            Find confidence files
            Find pLDDT files
```

### Step 2: Model Selection
```
For each result folder:
    Load confidence scores for all 5 models
    Select best model by confidence_score
    Record selected model path
```

### Step 3: Epitope Evaluation (Axis A)
```
For each complex:
    Extract GT epitope from crystal structure (5A cutoff)
    For K+:
        Extract predicted contacts from best model
        Compute epitope heatmap if available
        Threshold heatmap -> predicted epitope
    For baselines:
        Extract predicted contacts from best model
    Compute precision, recall, F1, MCC
    For K+ only: compute AUC-PR, AUC-ROC from heatmap
```

### Step 4: Docking Evaluation (Axis B)
```
For each complex, for each method:
    Run DockQ: selected_model vs ground_truth
    Extract: DockQ, iRMSD, LRMSD, fnat, F1 for each interface
    Compute: AB+AC average DockQ (antigen-antibody interface)
    Compute: antibody-aligned antigen RMSD
    Extract: confidence metrics from JSON
```

### Step 5: Aggregation & Comparison
```
For each metric:
    Compute mean, median, std across all complexes
    Compute per-complex delta (K+ - baseline)
    Statistical test: paired Wilcoxon signed-rank test
    Generate comparison table and plots
```

---

## 9. Summary Tables to Produce

### Table 1: Epitope Prediction (Axis A)
| Method | Precision | Recall | F1 | MCC | AUC-PR |
|--------|-----------|--------|----|-----|--------|
| B1 (vanilla) | | | | | N/A |
| B2 (contact restraints, best) | | | | | N/A |
| B3 (pocket restraints, best) | | | | | N/A |
| **K+ (ours)** | | | | | |

### Table 2: Docking Quality (Axis B)
| Method | DockQ (AB+AC) | DockQ (global) | iRMSD (AB+AC) | LRMSD | fnat (AB+AC) | Ab-aligned Ag RMSD |
|--------|---------------|----------------|---------------|-------|--------------|-------------------|
| B1 | | | | | | |
| B2 best | | | | | | |
| B3 best | | | | | | |
| **K+** | | | | | | |

### Table 3: Confidence Metrics
| Method | confidence_score | iptm | complex_plddt | complex_iplddt |
|--------|-----------------|------|---------------|----------------|
| B1 | | | | |
| **K+** | | | | |

### Table 4: Per-Complex Breakdown
One row per complex, showing DockQ and F1 for each method side by side.

---

## 10. Plots to Generate

1. **Bar chart**: Mean DockQ (AB+AC interface) per method with error bars
2. **Bar chart**: Mean epitope F1 per method with error bars
3. **Scatter plot**: K+ DockQ vs B1 DockQ (one point per complex, diagonal = no change)
4. **Scatter plot**: K+ epitope F1 vs B1 epitope F1
5. **Distribution**: DockQ quality categories (incorrect/acceptable/medium/high) per method
6. **PR curve**: Epitope precision-recall curve for K+ (using heatmap thresholds), with baseline F1 points overlaid
7. **Heatmap examples**: For 2-3 representative complexes, show the K+ epitope heatmap on the antigen sequence alongside the ground-truth epitope

---

## 11. Environments & Tools

| Task | Environment | Key Packages |
|------|-------------|-------------|
| DockQ scoring | `conda activate dockq2` | DockQ v2 |
| Epitope extraction, alignment, RMSD | `conda activate boltz` | BioPython 1.84, numpy, scipy |
| Evaluation script | `conda activate boltz` | BioPython, numpy, pandas, matplotlib, scikit-learn (for AUC) |

No new environment is needed. The `boltz` environment has BioPython and all scientific packages. The `dockq2` environment has the DockQ binary. The evaluation script will:
1. Run in the `boltz` environment for all analyses
2. Shell out to `conda run -n dockq2 DockQ ...` for DockQ scoring (or run DockQ as a subprocess)

---

## 12. Handling Baselines with Multiple Configurations

For B2 (contact restraints) and B3 (pocket restraints), each complex may have multiple prediction folders (e.g., different restraint types or residues). We handle this by:

1. **Best-oracle**: For each complex, take the configuration with the highest DockQ. This shows the ceiling when the right restraint is chosen.
2. **Mean**: Average across all configurations. This shows expected performance when the restraint is chosen randomly.
3. **Per-restraint-type** (B2 only): Group by restraint type (hbond, hydrophobic, salt_bridge) and report separately.

For the epitope evaluation, we combine contacts from all configurations to create a union epitope prediction (analogous to K+'s multi-region heatmap).

---

## 13. Statistical Testing

- **Paired Wilcoxon signed-rank test**: Compare K+ vs each baseline on the same complexes. Non-parametric, handles non-normal distributions.
- **Bootstrap confidence intervals**: 95% CI for mean differences.
- **Effect size**: Cohen's d or rank-biserial correlation.
- Report p-values for each metric comparison.

---

## 14. Execution Order

1. Write the evaluation script (`evaluate_k_plus.py`)
2. Run DockQ on all available predictions (subprocess calls)
3. Extract epitope from all predicted and ground-truth structures
4. Compute all metrics
5. Generate summary tables (CSV/markdown)
6. Generate plots (matplotlib)
7. Write summary report

The full evaluation will run on the cluster with all ~48 complexes. The script should be tested on the available subset first.
