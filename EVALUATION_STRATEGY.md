# Evaluation Strategy: L+ (CDR3 β-Scaling v2)

## 1. Feature Goals Recap

The CDR3 β-Scaling feature modulates pair representations in the CDR3 region during diffusion inference, with the goals of:
- **Improving CDR3 loop conformation quality** (structural accuracy vs crystal structures)
- **Improving antibody-antigen interface prediction** (docking quality)
- **Generating meaningful conformational diversity** in the CDR3 region

We compare against **3 baselines** across ~48 antibody-antigen complexes.

---

## 2. Methods Compared

| Label | Folder | Description |
|-------|--------|-------------|
| **Baseline 1** | `antigen_cut` | Standard Boltz2 prediction (no steering) |
| **Baseline 2** | `antigen_cut_contact_restraints` | Contact restraint-guided prediction (multiple restraint types per complex) |
| **Baseline 3** | `antigen_cut_vhvl_msa_pocket_ab_boltz_post_2023_omm` | Pocket-aware antibody prediction with VH/VL MSA |
| **New Feature** | `new_feature_cdr3_beta` | CDR3 β-Scaling (β=0.3) |

**Note on baselines 2 and 3**: These contain multiple result folders per complex (e.g., different restraint types, different pocket residues). Each sub-setting should be evaluated independently AND we should report the **best-of-settings** result per complex to represent the best achievable quality with that method.

---

## 3. Model Selection Strategy

Each method produces **5 model predictions** per complex (model_0 through model_4). For fair evaluation:

### 3.1 Primary: Best-of-5 Selection
Select the best model per complex based on **Boltz confidence score** (`confidence_score` from the JSON file, which combines pTM, ipTM, and pLDDT). This reflects what a practitioner would do: pick the highest-confidence prediction.

### 3.2 Secondary: All-5 Reporting
Also report statistics across all 5 models (mean ± std) to assess prediction consistency.

### 3.3 Rationale
- The confidence score is the model's own ranking metric and is standard practice for AlphaFold/Boltz model selection
- Reporting all-5 shows whether the method is consistently good or just gets lucky once

---

## 4. Evaluation Metrics

### 4.1 DockQ v2 — Antibody-Antigen Interface Quality (PRIMARY)

**Tool**: `DockQ` (conda environment: `dockq2`)
**Why**: DockQ v2 is the standard CAPRI-compliant metric for protein-protein docking quality. It directly measures how well the predicted interface matches the crystal structure.

**Metrics extracted per complex**:

| Metric | Description | Relevance |
|--------|-------------|-----------|
| **DockQ** | Combined quality score (0-1) | Overall docking quality |
| **fnat** | Fraction of native contacts recovered | Contact prediction accuracy |
| **F1** | Harmonic mean of precision/recall of contacts | Balanced contact score |
| **iRMSD** | Interface RMSD (Å) | Local interface geometry |
| **LRMSD** | Ligand RMSD after receptor alignment (Å) | Global positioning accuracy |
| **fnonnat** | Fraction of non-native contacts | False positive contacts |

**Interface decomposition**: DockQ evaluates all pairwise chain interfaces. For our antibody-antigen complexes (A=antigen, B=heavy chain, C=light chain):
- **A-B interface** (antigen–heavy chain): Most relevant for CDR3-H impact
- **A-C interface** (antigen–light chain): Most relevant for CDR3-L impact
- **B-C interface** (VH–VL): Intra-antibody packing (should remain stable)
- **Total DockQ**: Average across all interfaces

We report **A-B and A-C interfaces separately** and **total DockQ** as the headline number.

**CAPRI quality classification** (from DockQ):
- Incorrect: DockQ < 0.23
- Acceptable: 0.23 ≤ DockQ < 0.49
- Medium: 0.49 ≤ DockQ < 0.80
- High: DockQ ≥ 0.80

**Command**:
```bash
conda run -n dockq2 DockQ <prediction> <ground_truth> --json <output.json>
```

DockQ natively handles both `.pdb` and `.cif` formats, so no conversion is needed.

---

### 4.2 CDR3 Backbone RMSD — Loop Conformation Accuracy

**Tool**: BioPython (available in both `boltz` and `dockq2` environments)
**Why**: CDR3 loops are the most variable and functionally critical regions. Their structural accuracy is the direct target of β-scaling.

**Procedure**:
1. Parse prediction and ground truth structures (BioPython `PDBParser`/`MMCIFParser`)
2. Extract **antibody framework residues** (all antibody residues EXCLUDING CDR loops) from chains B and C
3. **Superimpose** prediction onto ground truth using framework Cα atoms (BioPython `Superimposer`)
4. Compute **backbone RMSD** (N, Cα, C, O atoms) of **CDR3-H** and **CDR3-L** residues separately
5. CDR3 residue indices are loaded from `examples/cdrs.csv` (0-indexed in the CSV; residue numbering in structures starts at 1 for each chain within the PDB but may use sequential numbering across chains — handle both cases)

**Metrics**:
- **CDR3-H backbone RMSD** (Å)
- **CDR3-L backbone RMSD** (Å)
- **CDR3-combined backbone RMSD** (Å) — all CDR3 residues together

**Note on alignment**: Aligning on the antibody framework (not the full complex) isolates CDR3 accuracy from global docking errors. This is the standard approach in antibody modeling benchmarks (e.g., ABodyBuilder2, ImmuneBuilder).

---

### 4.3 CDR3 Region pLDDT — Confidence in CDR3 Prediction

**Tool**: NumPy (load `.npz` files)
**Why**: pLDDT reflects the model's confidence in each residue's position. Higher CDR3 pLDDT with β-scaling indicates the model is more certain about CDR3 conformations.

**Procedure**:
1. Load `plddt_<complex>_model_<i>.npz` — shape `[n_tokens]`
2. Map CDR3 residue indices (from `examples/cdrs.csv`) to token positions
   - Chain B (heavy): token offset = len(chain_A_tokens)
   - Chain C (light): token offset = len(chain_A_tokens) + len(chain_B_tokens)
3. Extract pLDDT values for CDR3-H and CDR3-L tokens

**Metrics**:
- **Mean CDR3-H pLDDT**
- **Mean CDR3-L pLDDT**
- **Mean CDR3-combined pLDDT**
- **Complex-wide pLDDT** (from confidence JSON, for reference)

---

### 4.4 Ensemble Diversity — CDR3 Conformational Sampling

**Tool**: BioPython
**Why**: The feature should generate diverse but valid CDR3 conformations. Too little diversity means β-scaling has no effect; too much may indicate instability.

**Procedure**:
1. For each complex, take all 5 models
2. Align each pair of models using antibody framework Cα atoms
3. Compute pairwise CDR3 Cα RMSD between all 5 models → 10 pairs per complex

**Metrics**:
- **Mean pairwise CDR3 RMSD** (Å) — higher = more diversity
- **Max pairwise CDR3 RMSD** (Å)
- **Fraction of pairs with CDR3 RMSD ≥ 1.0 Å** — the diversity threshold from the specification

---

### 4.5 Boltz Confidence Metrics — Model Self-Assessment

**Tool**: JSON parsing
**Why**: These are the model's own quality predictions. They should not degrade with β-scaling if the feature is working correctly.

**Metrics** (directly from confidence JSON):
- **confidence_score** (overall)
- **complex_plddt**
- **ptm** (predicted TM-score)
- **iptm** (interface pTM, most relevant for docking)
- **protein_iptm**
- **pair_chains_iptm** (per-interface ipTM for A-B, A-C, B-C)

---

## 5. Preprocessing Requirements

### 5.1 File Format Handling
- Baselines produce `.pdb` files; new feature produces `.cif` files
- DockQ handles both formats natively — no conversion needed
- BioPython uses `PDBParser` for `.pdb` and `MMCIFParser` for `.cif`

### 5.2 Chain Mapping
- Ground truth and predictions both use chains A (antigen), B (heavy), C (light)
- DockQ auto-detects mapping via sequence alignment (verified to work correctly: `ABC:ABC`)
- For BioPython alignment: match chains by ID

### 5.3 Residue Numbering
- Ground truth PDBs use sequential numbering across chains (A:1-212, B:213-434, C:435-648)
- Predicted PDBs use per-chain numbering starting from 1 (confirmed from file inspection)
- CDR indices in `cdrs.csv` are 0-indexed positions within each chain's sequence
- **Must convert**: CDR CSV index → residue number in the structure (add 1 for predicted structures, add chain offset for ground truth)

### 5.4 Handling Multiple Settings per Baseline
For baselines 2 and 3 (which have multiple restraint types per complex):
- Run DockQ for each setting independently
- Report the **best-of-settings** result per complex (picks the restraint type that works best for each complex)
- Also report **per-setting averages** if useful for understanding which restraint types work best

---

## 6. Statistical Analysis

### 6.1 Per-Complex Comparison
For each complex, compute Δ(metric) = new_feature - baseline for all metrics.

### 6.2 Aggregate Statistics
- **Mean ± std** across all complexes for each metric
- **Median** (robust to outliers)
- **Win/loss/tie counts**: Number of complexes where new feature beats/loses to/ties each baseline

### 6.3 Statistical Tests
- **Wilcoxon signed-rank test** (paired, non-parametric) comparing new feature vs each baseline
- Report p-values for the primary metric (DockQ total, DockQ A-B, DockQ A-C)

### 6.4 Visualization
- Box plots of DockQ scores per method
- Scatter plots: new feature DockQ vs baseline DockQ (per complex)
- CDR3 RMSD comparison bar charts
- Heatmap of per-complex Δ(DockQ) across methods

---

## 7. Summary Table of Metrics

| # | Metric | Source | Level | Primary? |
|---|--------|--------|-------|----------|
| 1 | DockQ (total, A-B, A-C) | DockQ v2 | Per-complex | **Yes** |
| 2 | iRMSD, LRMSD | DockQ v2 | Per-interface | Yes |
| 3 | fnat, F1, fnonnat | DockQ v2 | Per-interface | Yes |
| 4 | CAPRI classification | DockQ v2 | Per-interface | Yes |
| 5 | CDR3-H backbone RMSD | BioPython | Per-complex | **Yes** |
| 6 | CDR3-L backbone RMSD | BioPython | Per-complex | **Yes** |
| 7 | CDR3 pLDDT | npz files | Per-complex | Yes |
| 8 | Ensemble CDR3 diversity | BioPython | Per-complex | Yes |
| 9 | Boltz confidence metrics | JSON files | Per-complex | Secondary |
| 10 | Wilcoxon p-values | scipy | Aggregate | Yes |

---

## 8. Environment and Dependencies

### Primary environment: `dockq2`
Use for: DockQ evaluation, BioPython-based RMSD computation, all analysis scripts

Available packages:
- DockQ v2 (command-line tool)
- BioPython 1.86
- NumPy 1.26.4
- Pandas 3.0.0
- Python 3.12

**Needs installation**: `matplotlib` (for plots) and `scipy` (for statistical tests)
```bash
conda activate dockq2
pip install matplotlib scipy
```

### Secondary environment: `boltz`
Use for: Running Boltz2 predictions only (not for evaluation)

### No new environment needed
All evaluation tasks can be done within `dockq2` with the two additional packages above.

---

## 9. Execution Plan

### Step 1: DockQ Evaluation
For each method × complex × model (0-4):
```bash
conda run -n dockq2 DockQ <prediction_file> pdb_minimized/<complex>.pdb --json <output.json>
```
Collect all JSON outputs into a single DataFrame.

### Step 2: Model Selection
For each method × complex: select best model by `confidence_score` from JSON.

### Step 3: CDR3 RMSD Computation
For each method × complex × model:
- Parse structures with BioPython
- Align on antibody framework Cα
- Compute CDR3-H and CDR3-L backbone RMSD

### Step 4: CDR3 pLDDT Extraction
For each method × complex × model:
- Load pLDDT npz
- Extract CDR3 region values using indices from `cdrs.csv`

### Step 5: Ensemble Diversity
For each method × complex:
- Compute all pairwise CDR3 Cα RMSD among the 5 models

### Step 6: Aggregate and Compare
- Build master results table (rows=complexes, columns=metrics×methods)
- Compute summary statistics and statistical tests
- Generate plots

### Step 7: Report
- Summary table of aggregate metrics per method
- Per-complex comparison heatmap
- Statistical significance results
- Key findings and recommendations

---

## 10. Expected Output Files

```
evaluation/
├── results/
│   ├── dockq_all.csv              # All DockQ results
│   ├── cdr3_rmsd_all.csv          # All CDR3 RMSD results
│   ├── cdr3_plddt_all.csv         # All CDR3 pLDDT results
│   ├── ensemble_diversity.csv     # Ensemble diversity metrics
│   ├── confidence_all.csv         # Boltz confidence metrics
│   └── summary_table.csv          # Aggregate comparison
├── plots/
│   ├── dockq_comparison_boxplot.png
│   ├── dockq_scatter_per_baseline.png
│   ├── cdr3_rmsd_comparison.png
│   ├── cdr3_plddt_comparison.png
│   ├── ensemble_diversity.png
│   └── per_complex_heatmap.png
└── evaluation_report.md           # Final written report
```
