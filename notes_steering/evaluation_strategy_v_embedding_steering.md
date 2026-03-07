# Evaluation Strategy: V (Embedding-Space CDR3 Steering)

**Date**: 2026-03-07
**Feature under evaluation**: Strategy V — Embedding-Space CDR3 Steering
**Baselines**: 3 (see below)
**Ground truth**: 48 crystal structures in `pdb_minimized/`

---

## 1. What Are We Evaluating and Why

### 1.1 Feature Goal Recap

Strategy V optimizes the Pairformer's pair representation `z` in embedding space to steer CDR3 loop conformations. It does **not** apply coordinate-space potentials. From the project documentation, the stated goals are:

- **Primary**: More robust CDR3 conformation exploration for novel sequences
- **Secondary**: Improved antibody-antigen structure prediction quality
- **Claimed advantage**: Stable across 2 orders of magnitude of hyperparameter variation (vs. brittle coordinate methods)

### 1.2 What This Means for Evaluation

The feature modifies the pair representation `z` at every denoising step, which affects the **entire complex** prediction — not just the CDR3 loops. Therefore, we should evaluate:

1. **Overall complex quality** — Did we hurt or help the global prediction?
2. **Antibody-antigen interface quality** — Did we improve the binding interface?
3. **CDR3-specific accuracy** — Did we improve the CDR3 loop conformations?
4. **Robustness** — Is performance consistent across different complexes?

### 1.3 Baselines

| ID | Folder | Description | Models per complex |
|----|--------|-------------|-------------------|
| **B1** | `antigen_cut` | Vanilla Boltz2 (no steering) | 5 PDB |
| **B2** | `antigen_cut_contact_restraints` | Boltz2 with contact restraint steering | 5 PDB (multiple settings per complex) |
| **B3** | `antigen_cut_vhvl_msa_pocket_ab_boltz_post_2023_omm` | Boltz2 with pocket/antibody-specific constraints | 5 PDB (multiple settings per complex) |
| **V** | `new_feature_cdr3_beta` (name may vary for embedding steering) | Strategy V: Embedding-Space CDR3 Steering | 5 PDB or CIF |

**Note**: Baselines B2 and B3 have multiple settings per complex (e.g., different restraint types). For a fair comparison, we should report both the **per-setting** results and the **best-per-complex** results (i.e., the best prediction across all settings for a given complex). This is standard practice — it answers "given all the tricks this baseline can try, what is its best outcome?"

---

## 2. Model Selection Strategy

Each method produces up to 5 structure predictions per complex per setting. We must decide how to pick the representative prediction for each complex.

### Recommended approach: **Best-by-confidence**

For each complex and each method:
1. Read the `confidence_*.json` files for all 5 models
2. Select the model with the **highest `confidence_score`** (which is a composite of pTM and ipTM)
3. Use that single model for all downstream evaluations

**Rationale**:
- `confidence_score` is Boltz2's own ranking metric — it is what a user would actually use to pick the best model
- Using the top-1 by confidence is standard practice in structure prediction benchmarking (e.g., CASP, AlphaFold evaluations)
- Alternatively, we will also report **mean across all 5 models** as a secondary analysis to assess consistency

### Secondary analysis: **All 5 models**

Report mean ± std of all metrics across the 5 models. This captures:
- How consistent is the method? (low std = robust)
- Does the best model come with high-quality alternatives? (useful for ensemble-based workflows)

---

## 3. Evaluation Metrics

### 3.1 Tier 1: DockQ (Primary Metric — Interface Quality)

**Tool**: DockQ v2.1.3 (`conda activate dockq2`)
**What it measures**: Protein-protein docking quality against crystal structure
**Why it's primary**: DockQ is the standard metric for assessing protein complex quality (CAPRI community standard). It directly answers "how good is the predicted antibody-antigen interface?"

DockQ provides per-interface scores. For our 3-chain complexes (A=antigen, B=heavy, C=light):

| Interface | Biological meaning | Relevance |
|-----------|--------------------|-----------|
| **A-B** (antigen–heavy chain) | Primary binding interface | **High** — CDR-H3 dominates antigen contact |
| **A-C** (antigen–light chain) | Secondary binding interface | **High** — CDR-L3 contributes to antigen contact |
| **B-C** (heavy–light chain) | VH-VL pairing | **Medium** — internal antibody quality |
| **Global** | Average of all interfaces | **High** — overall complex quality |

**Metrics extracted per interface:**

| Metric | Description | Range | Key threshold |
|--------|-------------|-------|---------------|
| **DockQ** | Combined docking quality | 0–1 | <0.23 incorrect, 0.23–0.49 acceptable, 0.49–0.80 medium, ≥0.80 high |
| **iRMSD** | Interface RMSD (Å) | 0–∞ | <1Å excellent, <2Å good, <4Å acceptable |
| **LRMSD** | Ligand RMSD after receptor alignment (Å) | 0–∞ | <1Å excellent, <5Å good, <10Å acceptable |
| **fnat** | Fraction of native contacts recovered | 0–1 | Higher is better |
| **F1** | Harmonic mean of precision and recall of contacts | 0–1 | Higher is better |
| **fnonnat** | Fraction of non-native contacts (false positives) | 0–1 | Lower is better |

**How to run**:
```bash
conda activate dockq2
DockQ <prediction.pdb|cif> <ground_truth.pdb> --json output.json
```

DockQ automatically handles chain mapping and different numbering schemes (verified on our data).

### 3.2 Tier 2: CDR3-Specific RMSD

**Tool**: BioPython (available in both `boltz` and `dockq2` envs)
**What it measures**: How accurately the CDR3 loops are predicted
**Why it matters**: Strategy V specifically targets CDR3 conformations. We need to know if CDR3 got better, even if the global interface metrics stay the same.

**Procedure**:
1. **Superpose** the antibody framework (non-CDR residues of chains B+C) between prediction and ground truth using Kabsch alignment on CA atoms
2. **Compute RMSD** of CDR3 CA atoms (CDR3-H and CDR3-L separately) **without further alignment** — this captures both conformation and positioning errors
3. CDR3 residue ranges come from `examples/cdrs.csv` (per-complex, not hardcoded)

**Important note on residue numbering**: Ground truth PDBs use continuous numbering across chains (chain B starts at residue 213, not 1). Predictions use per-chain numbering (chain B starts at 1). The evaluation script must map between these by matching sequences, not residue numbers. BioPython's `Superimposer` on matched CA atoms handles this naturally.

**Metrics**:

| Metric | Description |
|--------|-------------|
| **CDR3-H RMSD (Å)** | CA RMSD of heavy chain CDR3 after framework alignment |
| **CDR3-L RMSD (Å)** | CA RMSD of light chain CDR3 after framework alignment |
| **CDR3-combined RMSD (Å)** | CA RMSD of both CDR3 loops together |

### 3.3 Tier 3: Boltz2 Confidence Metrics

**Tool**: Python (numpy, json) — no special environment needed
**What it measures**: The model's own assessment of prediction quality
**Why it matters**: If Strategy V produces high-confidence predictions, users can trust them; if confidence drops while accuracy improves, that's a calibration concern

**Metrics extracted from `confidence_*.json`**:

| Metric | Description | Expected range |
|--------|-------------|----------------|
| **confidence_score** | Composite (0.8×ipTM + 0.2×pTM) | 0–1, higher is better |
| **iptm** | Interface pTM (inter-chain alignment quality) | 0–1 |
| **ptm** | Predicted TM-score (overall fold quality) | 0–1 |
| **complex_plddt** | Mean pLDDT across complex | 0–1 |
| **complex_iplddt** | Interface pLDDT | 0–1 |

**Metrics extracted from `plddt_*.npz`**:

| Metric | Description |
|--------|-------------|
| **CDR3-H pLDDT** | Mean pLDDT over CDR3-H residues |
| **CDR3-L pLDDT** | Mean pLDDT over CDR3-L residues |
| **Antibody pLDDT** | Mean pLDDT over chains B+C |
| **Antigen pLDDT** | Mean pLDDT over chain A |

### 3.4 Tier 4: Epitope Recovery (Contact Analysis)

**Tool**: BioPython (distance computation)
**What it measures**: Whether the predicted contacts between antibody and antigen match the crystal structure
**Why it matters**: Even if RMSD is high (global misalignment), the predicted interface residues might still be correct — this is crucial for epitope discovery applications.

**Procedure**:
1. Define a contact as any pair of residues (one from antibody, one from antigen) with minimum heavy-atom distance ≤ 5Å (standard cutoff)
2. Compute contacts in both crystal structure and prediction
3. Calculate precision, recall, and F1 for contact prediction

**Metrics**:

| Metric | Description |
|--------|-------------|
| **Contact precision** | Fraction of predicted contacts that are correct |
| **Contact recall** | Fraction of true contacts that are recovered |
| **Contact F1** | Harmonic mean of precision and recall |
| **Epitope residue precision** | Fraction of predicted antigen epitope residues that are correct |
| **Epitope residue recall** | Fraction of true antigen epitope residues that are recovered |

**Note**: fnat and F1 from DockQ already capture this at the interface level. This Tier 4 analysis provides a finer-grained, residue-level view focused on the antigen side (epitope).

---

## 4. Preprocessing Requirements

### 4.1 No Pre-alignment Needed for DockQ

DockQ performs its own internal alignment. No preprocessing required. Simply provide model and native PDB/CIF files.

### 4.2 Framework Alignment for CDR3 RMSD

For CDR3-specific RMSD, we need to superpose on the antibody **framework** (non-CDR residues):

1. Parse both prediction and ground truth structures
2. Match residues by sequence alignment (to handle different numbering)
3. Identify framework residues: all antibody residues that are **not** in any CDR (H1, H2, H3, L1, L2, L3 — all from `cdrs.csv`)
4. Extract CA atoms for framework residues in both structures
5. Kabsch alignment using framework CAs
6. Apply the rotation/translation to the **entire prediction**
7. Compute CDR3 RMSD on the transformed prediction vs. native

**Why framework alignment**: Aligning on the antibody framework isolates CDR3 conformational accuracy from any global positioning error of the antibody relative to the antigen. This is the standard approach in antibody structure prediction benchmarks (e.g., ABlooper, ImmuneBuilder, ABodyBuilder2).

### 4.3 Residue Mapping

Ground truth PDBs use continuous residue numbering across chains. Predictions use per-chain numbering. The mapping must be done by:
1. Extracting the sequence from each chain in both structures
2. Matching chains by sequence identity (not by chain ID, though in our case chain IDs are consistent: A, B, C)
3. Aligning residues within each chain by their position in the sequence

### 4.4 File Format Handling

- Baselines produce `.pdb` files
- Strategy V may produce `.cif` files
- Both DockQ and BioPython handle both formats natively (verified)

---

## 5. Statistical Analysis

### 5.1 Per-Complex Comparison

For each of the 48 complexes, compare Strategy V vs. each baseline:
- Compute Δ(metric) = metric_V - metric_baseline for each metric
- Report the distribution of deltas (mean, median, std)
- Count wins/ties/losses

### 5.2 Aggregate Statistics

| Statistic | Description |
|-----------|-------------|
| **Mean ± SEM** | Average performance across all complexes |
| **Median** | Robust central tendency (less sensitive to outliers) |
| **Wilcoxon signed-rank test** | Paired non-parametric test (appropriate for N=48, no normality assumption) |
| **Win/Loss ratio** | Fraction of complexes where V outperforms each baseline |

### 5.3 CAPRI Quality Categories

Classify each prediction into CAPRI categories based on DockQ:
- **Incorrect**: DockQ < 0.23
- **Acceptable**: 0.23 ≤ DockQ < 0.49
- **Medium**: 0.49 ≤ DockQ < 0.80
- **High**: DockQ ≥ 0.80

Report the **distribution** of categories for each method. An improvement that shifts predictions from "Acceptable" to "Medium" is more valuable than one that improves already-"High" predictions.

### 5.4 Stratified Analysis

Analyze performance stratified by:
1. **Baseline difficulty**: easy (B1 DockQ ≥ 0.80) vs. medium (0.49–0.80) vs. hard (<0.49)
2. **CDR3-H length**: short (≤10 res) vs. long (>10 res) — from `cdrs.csv`
3. **Antigen size**: small (<200 res) vs. large (≥200 res)

This reveals whether Strategy V helps more on hard cases (desired) or easy cases (less useful).

---

## 6. Handling Baselines B2 and B3 (Multiple Settings Per Complex)

Baselines B2 (`antigen_cut_contact_restraints`) and B3 (`antigen_cut_vhvl_msa_pocket_ab_boltz_post_2023_omm`) have multiple result folders per complex (e.g., different restraint types for B2, different pocket residues for B3).

### Approach: Two comparison modes

1. **Best-per-complex**: For each complex, take the setting+model that yields the highest DockQ. This answers: "Given all settings the baseline tried, what was its best outcome?"
2. **All-settings**: Treat each setting as an independent prediction. This gives a fuller picture but inflates sample count.

**Primary comparison** should use best-per-complex, as it is the fairest representation of what the baseline can achieve.

---

## 7. Reporting

### 7.1 Primary Results Table

| Method | DockQ (A-B) | DockQ (A-C) | DockQ (Global) | CDR3-H RMSD | CDR3-L RMSD | ipTM | pLDDT |
|--------|-------------|-------------|----------------|-------------|-------------|------|-------|
| B1 (vanilla) | — | — | — | — | — | — | — |
| B2 (contact) | — | — | — | — | — | — | — |
| B3 (pocket) | — | — | — | — | — | — | — |
| **V (embed steer)** | — | — | — | — | — | — | — |

### 7.2 Visualization

1. **Box plots**: DockQ distribution per method (one panel per interface)
2. **Scatter plot**: DockQ_V vs. DockQ_B1 per complex (diagonal = equal; above = V wins)
3. **Bar chart**: CAPRI category distribution per method
4. **CDR3 RMSD distribution**: Histogram or violin plot per method
5. **Stratified bar chart**: Mean DockQ by difficulty tier per method

### 7.3 Per-Complex Detail Table

Full table with all 48 complexes × all metrics for reproducibility and case-by-case analysis.

---

## 8. Computational Environment

| Task | Environment | Key packages |
|------|------------|--------------|
| **DockQ computation** | `conda activate dockq2` | DockQ 2.1.3, BioPython 1.86 |
| **CDR3 RMSD, pLDDT analysis, contact analysis** | `conda activate boltz` | BioPython 1.84, numpy, scipy |
| **Statistical tests, plotting** | `conda activate boltz` | scipy, matplotlib (or seaborn if available) |

No new conda environment is needed. Both `boltz` and `dockq2` have BioPython and numpy. DockQ is only in `dockq2`. Use `boltz` for everything except DockQ calls.

---

## 9. Evaluation Script Design

A single evaluation pipeline script should:

1. **Discover results**: Walk each baseline/feature folder, find all prediction files, group by complex name
2. **Run DockQ**: For each prediction × ground truth pair, run DockQ and save JSON output
3. **Compute CDR3 RMSD**: Load CDR indices from `cdrs.csv`, superpose on framework, compute CDR3 RMSD
4. **Extract confidence**: Parse confidence JSON and pLDDT npz files
5. **Compute contacts**: Extract antigen-antibody contacts, compute precision/recall/F1
6. **Select best models**: Apply model selection (best-by-confidence, best-by-DockQ)
7. **Aggregate and compare**: Compute summary statistics, run statistical tests
8. **Output**: CSV tables + plots

---

## 10. Summary of Priorities

| Priority | Metric | Why |
|----------|--------|-----|
| **1 (must have)** | DockQ (all interfaces) | Standard benchmark metric for complex quality |
| **2 (must have)** | CDR3-H/L RMSD (after framework alignment) | Directly measures what Strategy V targets |
| **3 (important)** | Confidence metrics (ipTM, pLDDT, CDR3 pLDDT) | Model self-assessment — relevant for user trust |
| **4 (nice to have)** | Epitope contact F1 | Residue-level interface accuracy |
| **5 (nice to have)** | Robustness analysis (stratification) | Characterizes when the method helps vs. hurts |
