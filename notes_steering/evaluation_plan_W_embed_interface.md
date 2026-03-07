# Evaluation Plan: Strategy W — Embedding-Based Interface Steering

**Date**: 2026-03-07
**Feature**: `EmbeddingInterfacePotential` (embedding-weighted CDR-antigen distance steering)
**Goal**: Assess whether using pair embeddings to weight interface steering improves antibody-antigen complex structure prediction compared to baselines.

---

## 1. Feature Goals and What We Are Evaluating

Strategy W addresses a specific hypothesis: **the model's pair representation `z` encodes which CDR-antigen residue pairs are likely to form productive binding interactions**. By weighting distance-based steering with embedding magnitudes, we expect:

1. **Better antibody-antigen interface quality** — the antigen should dock closer to its correct orientation relative to the CDR loops
2. **More accurate contact prediction** — the predicted contacts should better match the crystal structure contacts (epitope recovery)
3. **Comparable or better overall structure quality** — the antibody/antigen folds themselves should not degrade

These three axes define our evaluation tiers below.

---

## 2. Experimental Setup

### 2.1 Conditions (4 total)

| ID | Condition | Folder Pattern | Description |
|----|-----------|----------------|-------------|
| **B1** | Baseline: no steering | `predictions_examples/antigen_cut/boltz_results_<COMPLEX>/` | Boltz2 default, no potentials |
| **B2** | Baseline: contact restraints | `predictions_examples/antigen_cut_contact_restraints/boltz_results_restraint_<COMPLEX>_*/` | Coordinate-space contact restraints (multiple sub-variants per complex: hbond, hydrophobic, salt_bridge) |
| **B3** | Baseline: pocket + MSA | `predictions_examples/antigen_cut_vhvl_msa_pocket_ab_boltz_post_2023_omm/boltz_results_restraint_to_A_<COMPLEX>_*/` | Pocket constraints targeting specific residues |
| **W** | Strategy W | `predictions_examples/new_feature_embed_interface/boltz_results_<COMPLEX>_embed_interface/` | Embedding-based interface steering (our new feature) |

> **Note on B2 and B3**: These baselines have multiple sub-folders per complex (e.g., different restraint types or target residues). For fair comparison, we will report:
> - **Best-of-variants**: the variant with the highest DockQ for the antigen-antibody interface
> - **Mean-across-variants**: average across all variants for robustness

### 2.2 Model Selection Strategy

Each condition produces up to **5 model predictions** per complex (`model_0` through `model_4`). Some conditions may produce fewer (e.g., only 1). Output format may be `.pdb` or `.cif`.

**Selection approaches** (report both):

| Approach | Metric | Rationale |
|----------|--------|-----------|
| **Best-by-confidence** | Highest `confidence_score` from JSON | The model's own ranking; tests whether confidence correlates with quality |
| **Best-by-interface-plddt** | Highest `complex_iplddt` from JSON | Interface-specific confidence; more relevant for docking quality |
| **All-5 mean** | Average metric across all 5 models | Measures consistency/robustness |

**Primary reporting**: Use **best-by-interface-plddt** for main results (since we care about interface quality). Report all-5-mean in supplementary.

### 2.3 Complexes

48 antibody-antigen complexes with ground truth crystal structures in `pdb_minimized/`. Complex naming: `<PDB_ID>_<CHAIN_SUFFIX>` (e.g., `7TRH_HBG`).

**Chain mapping** (consistent across all conditions):
- Chain A = antigen
- Chain B = heavy chain (VH only)
- Chain C = light chain (VL only)

**Important**: The ground truth PDB files contain full-length antibody chains with hydrogens, while predictions contain only variable domains (VH/VL) with heavy atoms only. DockQ handles this mismatch via internal sequence alignment (`--no_align` should NOT be used).

---

## 3. Evaluation Metrics (3 Tiers)

### Tier 1: Interface Docking Quality (Primary — DockQ)

**Tool**: DockQ v2 (`conda activate dockq2`)
**Rationale**: DockQ is the CAPRI-standard metric for protein docking assessment. It combines interface RMSD (iRMSD), ligand RMSD (LRMSD), and fraction of native contacts (fnat) into a single score.

**What to compute per prediction**:

```bash
DockQ <model>.pdb <native>.pdb --mapping ABC:ABC --json <output>.json
```

**Metrics extracted from DockQ JSON**:

| Metric | Interface | Description | Why it matters |
|--------|-----------|-------------|----------------|
| `DockQ` (AB) | Heavy–Antigen | Overall interface quality for VH–antigen | Primary metric for heavy chain binding |
| `DockQ` (AC) | Light–Antigen | Overall interface quality for VL–antigen | Primary metric for light chain binding |
| `DockQ` (BC) | VH–VL | Intra-antibody packing | Sanity check: should not degrade |
| `iRMSD` (AB, AC) | Antigen interfaces | Interface backbone RMSD | Spatial accuracy of binding site |
| `LRMSD` (AB, AC) | Antigen interfaces | Ligand RMSD (after receptor alignment) | Overall docking pose accuracy |
| `fnat` (AB, AC) | Antigen interfaces | Fraction of native contacts recovered | Contact prediction accuracy |
| `fnonnat` (AB, AC) | Antigen interfaces | Fraction of non-native contacts | False positive rate |
| `F1` (AB, AC) | Antigen interfaces | Harmonic mean of precision and recall | Balanced contact metric |
| **Total DockQ** | All 3 interfaces | Average DockQ across all interfaces | Global quality summary |

**CAPRI quality categories** (per interface):

| Category | DockQ range |
|----------|-------------|
| Incorrect | < 0.23 |
| Acceptable | 0.23 – 0.49 |
| Medium | 0.49 – 0.80 |
| High | ≥ 0.80 |

**Aggregate statistics to report**:
- Mean and median DockQ (AB), DockQ (AC), Total DockQ across 48 complexes
- Fraction of complexes in each CAPRI category
- Paired Wilcoxon signed-rank test (W vs each baseline)
- Per-complex improvement: ΔDockQ = DockQ_W − DockQ_baseline

### Tier 2: Epitope Recovery (Contact Prediction)

**Tool**: BioPython (`conda activate dockq2` or `conda activate boltz`)
**Rationale**: Strategy W aims to improve *which* antigen residues are contacted. This tests whether the predicted interface residues match the true epitope.

**Method**:
1. Define **ground truth epitope**: antigen residues with any heavy atom within 4.5 Å of any antibody heavy atom in the crystal structure
2. Define **predicted epitope**: same criterion applied to the predicted structure
3. Compute classification metrics treating each antigen residue as a binary prediction (epitope or not)

**Metrics**:

| Metric | Formula | Description |
|--------|---------|-------------|
| **Epitope Precision** | TP / (TP + FP) | Of predicted contacts, how many are real |
| **Epitope Recall** | TP / (TP + FN) | Of real contacts, how many were predicted |
| **Epitope F1** | 2 × Prec × Rec / (Prec + Rec) | Balanced epitope accuracy |
| **Epitope MCC** | Matthews Correlation Coefficient | Accounts for class imbalance |
| **Epitope overlap (Jaccard)** | \|Pred ∩ True\| / \|Pred ∪ True\| | Set-level similarity |

> Note: fnat/fnonnat from DockQ already capture contact accuracy at the *contact-pair* level. Epitope metrics here operate at the *residue* level and are complementary.

### Tier 3: Interface Structural Properties (SASA-based)

**Tool**: BioPython ShrakeRupley SASA (`conda activate dockq2` or `conda activate boltz`)
**Rationale**: The evaluation metrics from the plan document specify Buried Surface Area (SASA) and buried residue count. These measure whether the predicted interface has a physically realistic size.

**Method**:
1. Compute SASA for each chain in isolation (unbound)
2. Compute SASA for the full complex (bound)
3. **Buried Surface Area (BSA)** = SASA_unbound − SASA_bound (for the antigen-antibody interface)
4. **Buried residue count** = number of antigen residues where SASA decreases by ≥ 1 Å² upon binding

**Metrics**:

| Metric | Description |
|--------|-------------|
| **BSA (Å²)** | Total buried surface area at antigen-antibody interface |
| **BSA ratio** | BSA_predicted / BSA_native |
| **Buried residue count** | Number of antigen residues buried upon binding |
| **Buried residue ratio** | Count_predicted / Count_native |

---

## 4. Confidence Correlation Analysis

**Rationale**: Strategy W should not only improve structure quality but also produce predictions where the model's confidence is aligned with actual quality.

**Metrics**:
- Pearson/Spearman correlation between `complex_iplddt` and DockQ (AB+AC)
- Pearson/Spearman correlation between `confidence_score` and Total DockQ
- Compare these correlations across conditions

---

## 5. Postprocessing Notes

### 5.1 No Manual Alignment Needed for DockQ

DockQ v2 internally performs optimal superposition. It:
- Uses sequence alignment to match model chains to native chains (handles VH/VL vs full-length)
- Computes iRMSD, LRMSD, and fnat with proper alignment
- **Do NOT pre-align structures** before running DockQ

### 5.2 For Custom RMSD Calculations

If computing additional RMSD values beyond DockQ (e.g., CDR3-specific RMSD):
- **Align on the antibody (B+C) CA atoms**, then measure antigen displacement → this gives L-RMSD
- **Align on the antigen (A) CA atoms**, then measure CDR displacement → CDR positioning accuracy
- Use BioPython `Superimposer` with CA atoms only
- Match residues by sequence alignment (not residue numbers) due to length mismatch

### 5.3 File Format Handling

- Predictions may be `.pdb` or `.cif` — check both extensions when scanning folders
- DockQ accepts both formats
- For BioPython: use `PDBParser` for `.pdb`, `MMCIFParser` for `.cif`

---

## 6. Execution Plan

### Step 1: Discover and index all predictions

Scan all four condition folders. For each complex, collect:
- Path to each model file (`.pdb` or `.cif`)
- Path to corresponding confidence JSON
- Path to plddt/pae `.npz` files
- Select best model by `complex_iplddt`

Output: `evaluation_index.csv` with columns: `complex, condition, variant, model_idx, model_path, confidence_path, confidence_score, complex_iplddt`

### Step 2: Run DockQ on all predictions

```bash
conda activate dockq2
# For each row in evaluation_index.csv:
DockQ <model_path> pdb_minimized/<complex>.pdb --mapping ABC:ABC --json <output_json>
```

Output: `dockq_results.csv` with all DockQ metrics per prediction

### Step 3: Compute epitope metrics

```python
# conda activate dockq2  (BioPython available)
# For each complex:
#   1. Extract ground truth epitope from pdb_minimized/<complex>.pdb
#   2. Extract predicted epitope from best model
#   3. Compute precision, recall, F1, MCC, Jaccard
```

Output: `epitope_results.csv`

### Step 4: Compute SASA / BSA metrics

```python
# conda activate dockq2  (BioPython SASA available)
# For each complex:
#   1. Compute SASA unbound (chains isolated) and bound (complex)
#   2. Compute BSA and buried residue count
#   3. Compare predicted vs native
```

Output: `sasa_results.csv`

### Step 5: Aggregate and compare

```python
# Merge all result CSVs
# Compute aggregate statistics per condition
# Statistical tests: paired Wilcoxon for W vs each baseline
# Generate plots: boxplots, scatter plots, bar charts
```

Output: `evaluation_summary.csv`, `evaluation_plots/`

---

## 7. Environment Summary

| Task | Environment | Key Packages |
|------|-------------|-------------|
| DockQ scoring | `conda activate dockq2` | DockQ, BioPython |
| Epitope analysis | `conda activate dockq2` | BioPython (PDBParser, NeighborSearch) |
| SASA computation | `conda activate dockq2` | BioPython (ShrakeRupley SASA) |
| Confidence parsing | `conda activate dockq2` or `boltz` | json, numpy (standard) |
| Aggregation & plots | `conda activate dockq2` or `boltz` | pandas, matplotlib, scipy |

No new conda environment is needed. The `dockq2` environment has everything: DockQ, BioPython (with SASA), numpy, and standard libraries. Use `boltz` only if specific boltz imports are needed.

---

## 8. Expected Outputs Summary

| File | Contents |
|------|----------|
| `evaluation_index.csv` | Index of all predictions with confidence scores |
| `dockq_results.csv` | DockQ metrics for all predictions (per interface) |
| `epitope_results.csv` | Epitope classification metrics per complex |
| `sasa_results.csv` | BSA and buried residue metrics per complex |
| `evaluation_summary.csv` | Aggregate statistics per condition |
| `evaluation_plots/` | Boxplots, scatter plots, comparison figures |

---

## 9. Key Comparisons

The central question for each metric: **Does Strategy W outperform the baselines?**

| Comparison | What it tests |
|------------|---------------|
| **W vs B1** (no steering) | Does embedding steering help at all? |
| **W vs B2** (contact restraints) | Does embedding steering beat explicit coordinate-space contacts? |
| **W vs B3** (pocket + MSA) | Does embedding steering beat pocket-guided docking? |

For B2 and B3 (which have multiple variants per complex), compare W against:
- **Best variant**: generous comparison — does W beat the best possible restraint choice?
- **Mean variant**: practical comparison — does W beat the average restraint performance?

### Statistical Reporting

- Per-complex ΔDockQ (paired differences)
- Wilcoxon signed-rank test p-values
- Win/loss/tie counts (threshold: ΔDockQ > 0.02 for meaningful difference)
- Effect size (median ΔDockQ)
