# Evaluation Strategy: O+ (Asymmetric β-Scaling v2)

## 1. Feature Goals Recap

Asymmetric β-Scaling encodes the biological prior that **CDR-H3 dominates antigen binding (~60% of contacts)** while **CDR-L3 plays a supporting role (~40%)**. It applies differential denoising step scaling during diffusion sampling: `step * (1 + β_H)` for H3 atoms and `step * (1 + β_L)` for L3 atoms, with defaults β_H=0.4, β_L=0.1.

**What this should improve:**
- Antibody-antigen interface quality (H3 converges faster toward predicted structure)
- CDR-H3 structural accuracy (larger steps = stronger denoising signal)
- Overall docking quality by respecting the H3-dominant binding biology

**What this should NOT hurt:**
- Global fold quality (non-CDR regions are unscaled at β=0)
- CDR-L3 accuracy (moderate positive scaling, not aggressive)

## 2. Experimental Setup

### 2.1 Methods to Compare (4 total)

| Label | Description | Folder Pattern |
|-------|-------------|----------------|
| **Baseline (no steering)** | Boltz2 with antigen-cut MSA, no constraints | `antigen_cut/boltz_results_{COMPLEX}/` |
| **Pocket restraints** | Boltz2 with pocket binder constraints | `antigen_cut_vhvl_msa_pocket_ab_boltz_post_2023_omm/boltz_results_restraint_to_A_{COMPLEX}_*` |
| **Contact restraints** | Boltz2 with specific contact restraints | `antigen_cut_contact_restraints/boltz_results_restraint_{COMPLEX}_*` |
| **O+ (Asymmetric β)** | Our new feature | `new_feature_cdr3_beta/boltz_results_{COMPLEX}_*/` (or equivalent) |

### 2.2 Dataset

- **48 antibody-antigen complexes** from the test set
- **Ground truths**: energy-minimized crystal structures in `pdb_minimized/` (PDB format)
- All complexes have chains: A (antigen), B (heavy chain), C (light chain)

### 2.3 Model Selection Strategy

Each method produces **5 structure predictions** (model_0 through model_4). We report three selection strategies:

1. **Top-1 (confidence-ranked)**: Select the model with the highest `confidence_score` from the confidence JSON. This is the most realistic scenario (what a user would pick).
2. **Oracle (best)**: Select the model with the best score for each metric independently. This measures the method's ceiling.
3. **Average (all-5)**: Mean across all 5 models. This measures consistency.

**Implementation**: Read `confidence_{name}_model_{i}.json` → extract `confidence_score` → pick argmax for top-1.

## 3. Evaluation Metrics

### 3.1 Primary Metrics: Interface Quality (DockQ v2)

**Tool**: DockQ v2 (available in `conda activate dockq2`)

**Why**: DockQ is the gold standard for assessing protein-protein docking quality. It combines interface RMSD (iRMSD), ligand RMSD (LRMSD), and fraction of native contacts (fnat) into a single score. It maps directly to CAPRI quality categories.

**Command**:
```bash
DockQ <model.pdb> <native.pdb> --json <output.json>
```
DockQ auto-detects interfaces and computes per-interface scores with optimal chain mapping.

**Metrics extracted per interface**:

| Metric | Description | Why It Matters |
|--------|-------------|----------------|
| **DockQ** | Combined docking quality (0-1) | Overall interface accuracy; CAPRI-equivalent |
| **iRMSD** | Interface RMSD (Å) | How well interface residue positions are predicted |
| **LRMSD** | Ligand RMSD (Å) | How well the "ligand" chain is positioned relative to the "receptor" after receptor alignment |
| **fnat** | Fraction of native contacts recovered | Contact prediction accuracy |
| **fnonnat** | Fraction of non-native contacts | False positive rate |
| **F1** | Harmonic mean of precision and recall of contacts | Balanced contact accuracy |

**Interfaces to evaluate**:
- **A:B (antigen–heavy chain)**: The primary interface our feature targets
- **A:C (antigen–light chain)**: Secondary interface
- **B:C (heavy–light chain)**: Intra-antibody interface (should be preserved, not degraded)
- **Overall (average of all interfaces)**: Global docking quality

**CAPRI quality categories** (from DockQ score):
- Incorrect: DockQ < 0.23
- Acceptable: 0.23 ≤ DockQ < 0.49
- Medium: 0.49 ≤ DockQ < 0.80
- High: DockQ ≥ 0.80

**Aggregate statistics to report**:
- Mean DockQ per method (with 95% CI via bootstrap)
- Fraction of complexes with DockQ ≥ 0.23 (acceptable+)
- Fraction of complexes with DockQ ≥ 0.49 (medium+)
- Per-interface breakdown (A:B vs A:C vs B:C)

### 3.2 Feature-Specific Metric: H3/L3 Contact Score Ratio

**Why**: The biological hypothesis is that H3 dominates binding. If asymmetric β-scaling works, the H3 region should form more native contacts with the antigen than L3, at a ratio of ~1.5–2.0.

**Implementation** (BioPython, `conda activate dockq2`):
1. Parse the predicted PDB and the ground truth PDB
2. Identify CDR-H3 and CDR-L3 residues using indices from `cdrs.csv`
3. Identify antigen residues (chain A)
4. Count contacts: residue pairs where any heavy atom is within 5 Å across the interface
5. Compute:
   - `H3_contacts` = number of H3-antigen native contacts recovered in prediction
   - `L3_contacts` = number of L3-antigen native contacts recovered in prediction
   - `H3/L3 ratio` = H3_contacts / L3_contacts (target: 1.5–2.0)

**What to compare**: The H3/L3 contact ratio across methods. Our feature should produce ratios closer to the native structure's ratio.

### 3.3 CDR-H3 RMSD (Region-Specific Accuracy)

**Why**: Since we apply stronger scaling to H3, we expect its backbone structure to improve.

**Implementation** (BioPython, `conda activate dockq2`):
1. Align prediction to ground truth on the **antibody framework** (heavy chain non-CDR Cα atoms, approximately residues 1–25 and 35–95 for heavy chain)
2. After alignment, compute RMSD of **CDR-H3 Cα atoms only**
3. Similarly compute CDR-L3 RMSD
4. Compare H3 RMSD reduction vs baseline

**Why framework alignment**: Aligning on the full complex or the antigen would confound CDR accuracy with docking accuracy. Framework alignment isolates the CDR loop prediction quality.

### 3.4 Boltz Confidence Metrics (Internal Quality)

**Why**: These are free (already computed) and measure the model's self-assessment.

**Extracted from** `confidence_{name}_model_{i}.json`:

| Metric | Description |
|--------|-------------|
| `confidence_score` | Overall confidence (0-1) |
| `iptm` | Interface predicted TM-score |
| `complex_plddt` | Complex-wide predicted LDDT |
| `complex_iplddt` | Interface predicted LDDT |
| `pair_chains_iptm` | Per-chain-pair interface pTM |

**Key comparison**: Does asymmetric β-scaling produce structures with comparable or better self-confidence? Lower confidence with better DockQ would suggest the model's confidence calibration doesn't capture the improvement (still interesting to report).

### 3.5 Ensemble Diversity

**Why**: The method documentation suggests it should produce meaningful conformational diversity.

**Implementation**:
1. For each complex and method, compute pairwise Cα RMSD across the 5 models (focusing on CDR-H3 region)
2. Report mean pairwise RMSD as a diversity measure
3. Higher diversity with maintained quality = better exploration of conformational space

## 4. Postprocessing Pipeline

### 4.1 File Format Handling

- Ground truths: PDB format (✓ compatible with DockQ)
- Baseline predictions: PDB format (✓)
- New feature predictions: May be CIF format (observed `_model_0.cif` in example)
- **Action**: If CIF, convert to PDB using BioPython before DockQ evaluation, or use DockQ directly (it supports both PDB and CIF)

### 4.2 Chain Mapping

All structures (ground truths and predictions) use consistent chain IDs:
- A = antigen, B = heavy chain, C = light chain
- DockQ mapping: `--mapping ABC:ABC` (or let auto-detect handle it)

### 4.3 Alignment Strategy

DockQ internally handles alignment for its metrics (iRMSD, LRMSD, fnat). No pre-alignment needed for DockQ.

For CDR-H3 RMSD (Section 3.3), we do our own alignment:
1. Extract antibody framework Cα atoms (non-CDR residues of chains B and C)
2. Superimpose prediction onto ground truth using these atoms (BioPython Superimposer)
3. Compute RMSD on CDR-H3 Cα atoms in the aligned coordinates

### 4.4 Handling Multiple Restraint Variants

For **pocket restraints** and **contact restraints**, multiple variants exist per complex (e.g., `restraint_7TRH_HBG_hbond_23`, `restraint_7TRH_HBG_hydrophobic_5`). Strategy:
- Evaluate each variant independently
- For aggregate comparison, report the **best-performing variant per complex** (oracle across restraint types). This gives the baselines their best chance.

## 5. Execution Plan

### 5.1 Environment

**Primary**: `conda activate dockq2`
- DockQ v2 (interface quality)
- BioPython 1.86 (structure parsing, alignment, contact analysis)
- NumPy 1.26, Pandas 3.0, SciPy 1.17 (data analysis)

No additional conda environment needed. The `dockq2` environment has all required packages.

### 5.2 Evaluation Script Structure

```
scripts/eval_asymmetric_beta/
├── evaluate_all.py          # Main orchestrator
├── run_dockq.py             # Run DockQ on all predictions vs ground truths
├── compute_cdr_metrics.py   # CDR-H3/L3 RMSD and contact analysis
├── extract_confidence.py    # Extract Boltz confidence metrics
├── aggregate_results.py     # Aggregate, compute stats, generate tables/plots
└── utils.py                 # Shared utilities (file discovery, chain mapping)
```

### 5.3 Step-by-Step Execution

**Step 1: File discovery and validation**
- Scan all prediction folders, map each to its complex name
- Verify ground truth exists for each complex
- Handle CIF→PDB conversion if needed
- Output: `file_manifest.csv` with columns: `complex, method, model_idx, pred_path, gt_path`

**Step 2: DockQ evaluation**
```bash
# For each (prediction, ground_truth) pair:
conda activate dockq2
DockQ <pred.pdb> <gt.pdb> --json <output.json>
```
- Parse JSON outputs into a single DataFrame
- Output: `dockq_results.csv`

**Step 3: CDR-specific metrics**
- For each prediction: compute H3 RMSD, L3 RMSD, H3/L3 contact ratio
- Uses BioPython + CDR indices from `examples/cdrs.csv`
- Output: `cdr_metrics.csv`

**Step 4: Confidence extraction**
- Read all `confidence_*.json` files
- Extract: confidence_score, iptm, complex_plddt, complex_iplddt
- Output: `confidence_metrics.csv`

**Step 5: Aggregation and visualization**
- Merge all CSVs
- Apply model selection (top-1, oracle, average)
- Compute per-method statistics with 95% bootstrap CIs
- Generate comparison tables and bar plots

### 5.4 Output Artifacts

```
evaluation_results/
├── file_manifest.csv            # All file paths
├── dockq_results.csv            # DockQ per prediction
├── cdr_metrics.csv              # CDR-H3/L3 metrics per prediction
├── confidence_metrics.csv       # Boltz confidence per prediction
├── summary_table.csv            # Aggregated results per method
├── plots/
│   ├── dockq_comparison.pdf     # Bar plot: DockQ across methods
│   ├── cdr_rmsd_comparison.pdf  # Bar plot: H3/L3 RMSD
│   ├── contact_ratio.pdf        # H3/L3 contact ratio distribution
│   ├── per_complex_dockq.pdf    # Per-complex scatter/heatmap
│   └── interface_breakdown.pdf  # A:B vs A:C vs B:C DockQ
└── summary_report.md            # Narrative summary with key findings
```

## 6. Summary Table of Metrics

| Metric | Tool | Level | Primary Question Answered |
|--------|------|-------|---------------------------|
| DockQ (overall) | DockQ v2 | Per-interface | Is the Ab-Ag docking better? |
| DockQ (A:B only) | DockQ v2 | Per-interface | Is the heavy chain–antigen interface better? |
| iRMSD | DockQ v2 | Per-interface | Are interface residues more accurately positioned? |
| LRMSD | DockQ v2 | Per-interface | Is the antibody correctly positioned relative to antigen? |
| fnat | DockQ v2 | Per-interface | Are native contacts recovered? |
| F1 | DockQ v2 | Per-interface | Is contact prediction balanced (precision + recall)? |
| CAPRI category | DockQ v2 | Per-complex | What fraction achieve acceptable/medium/high quality? |
| CDR-H3 RMSD | BioPython | Per-complex | Is the H3 loop structure more accurate? |
| CDR-L3 RMSD | BioPython | Per-complex | Is L3 maintained (not degraded)? |
| H3/L3 contact ratio | BioPython | Per-complex | Does the ratio match biology (~1.5–2.0)? |
| H3 native contact count | BioPython | Per-complex | Does H3 form more correct contacts? |
| confidence_score | Boltz JSON | Per-model | Does the model's self-assessment agree? |
| iptm | Boltz JSON | Per-model | Is interface pTM improved? |
| Pairwise CDR-H3 RMSD | BioPython | Per-ensemble | Is conformational diversity maintained? |

## 7. Statistical Analysis

- **Primary comparison**: Paired test (Wilcoxon signed-rank) between O+ and each baseline, per complex, for DockQ
- **Multiple comparisons**: Report raw p-values + Bonferroni correction (3 baselines)
- **Effect size**: Report mean improvement ΔDockQ and its 95% bootstrap CI
- **Per-complex breakdown**: Show scatter plot of O+ DockQ vs baseline DockQ with identity line to identify which complexes benefit
