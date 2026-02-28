# Evaluation Strategy: Strategy G+ vs Baselines

## 1. Overview

We evaluate **Strategy G+ (Progressive CDR Refinement v2)** against three baselines to quantify the effect of adaptive phase scheduling with CDR3 beta-scaling on antibody-antigen docking quality.

### Methods Under Comparison

| Label | Folder | Description |
|-------|--------|-------------|
| **Baseline 1** (B1) | `antigen_cut/` | No steering. Standard Boltz2 prediction with antigen MSA only |
| **Baseline 2** (B2) | `antigen_cut_vhvl_msa_pocket_ab_boltz_post_2023_omm/` | Pocket-based contact restraints (residue-specific) |
| **Baseline 3** (B3) | `antigen_cut_contact_restraints/` | Contact-type restraints (hbond, hydrophobic, salt bridge) |
| **G+** | `new_feature_cdr3_beta/` | Strategy G+ with adaptive phases + CDR3 beta-scaling |

### Data Locations

```
Ground truth:
  /path/to/boltz_evals/pdb_minimized/               # 48 PDB files (minimized crystal structures)

Predictions:
  /path/to/predictions_examples/
  ├── antigen_cut/                                    # B1: 5 models/complex, PDB format
  ├── antigen_cut_vhvl_msa_pocket_ab_boltz_post_2023_omm/  # B2: 5 models/complex, PDB format
  ├── antigen_cut_contact_restraints/                 # B3: 5 models/complex, PDB format
  └── new_feature_cdr3_beta/                          # G+: 5 models/complex, CIF format
```

### Chain Convention

All structures (ground truth + predictions) use:
- **Chain A**: Antigen
- **Chain B**: Heavy chain (VH)
- **Chain C**: Light chain (VL)

---

## 2. Evaluation Metrics

### 2.1 Primary Metric: DockQ Score

**Tool**: DockQ v2.1.3 (installed in `dockq2` conda environment)

DockQ is a continuous score [0, 1] that combines three CAPRI metrics into a single quality measure:

```
DockQ = (fnat + 1/(1+(LRMSD/8.5)^2) + 1/(1+(iRMSD/1.5)^2)) / 3
```

**Component metrics** (all computed by DockQ internally):

| Metric | Definition |
|--------|------------|
| **fnat** | Fraction of native contacts preserved (contact = any heavy atom pair < 5 A) |
| **L-RMSD** | Ligand RMSD after superimposing receptor (antibody aligned, antigen RMSD measured) |
| **I-RMSD** | Interface RMSD (superimpose only interface residues within 10 A of partner) |

**Quality classification thresholds** (CAPRI standard):

| Quality | DockQ | fnat | L-RMSD (A) | I-RMSD (A) |
|---------|-------|------|------------|------------|
| Incorrect | < 0.23 | < 0.1 | > 10.0 | > 4.0 |
| Acceptable | 0.23-0.49 | >= 0.1 | <= 10.0 | <= 4.0 |
| Medium | 0.49-0.80 | >= 0.3 | <= 5.0 | <= 2.0 |
| High | >= 0.80 | >= 0.5 | <= 1.0 | <= 1.0 |

**Why DockQ**: It handles alignment internally (superimposes antibody BC, measures antigen A displacement), works with both PDB and CIF formats, supports multi-chain complexes, and is the community standard for protein docking evaluation (used in AlphaFold-Multimer, Chai-1, Boltz-1 papers).

**DockQ command**:
```bash
conda activate dockq2
DockQ prediction.pdb native.pdb --mapping BC:BC A:A --json output.json
```

### 2.2 Boltz2 Confidence Metrics (from prediction outputs)

These are extracted directly from the confidence JSON files without additional computation:

| Metric | Source | Description |
|--------|--------|-------------|
| **confidence_score** | confidence JSON | Overall Boltz2 confidence (composite) |
| **iptm** | confidence JSON | Interface predicted TM-score (inter-chain quality) |
| **protein_iptm** | confidence JSON | Protein-specific interface pTM |
| **complex_plddt** | confidence JSON | Mean pLDDT across all residues |
| **complex_iplddt** | confidence JSON | Mean pLDDT at interface residues only |
| **complex_ipde** | confidence JSON | Interface PDE score |
| **pair_chains_iptm** | confidence JSON | Pairwise chain iPTM (A-B, A-C, B-C) |

### 2.3 Per-Residue Quality Metrics (from NPZ files)

| Metric | Source | Shape | Description |
|--------|--------|-------|-------------|
| **pLDDT** | plddt NPZ | (N_residues,) | Per-residue confidence |
| **PAE** | pae NPZ | (N_res, N_res) | Predicted aligned error (pairwise) |
| **PDE** | pde NPZ | (N_res, N_res) | Predicted distance error (pairwise) |

**Derived metrics from per-residue data**:
- **CDR3 pLDDT**: Mean pLDDT over CDR3-H and CDR3-L residues (from cdrs.csv)
- **Interface PAE**: Mean PAE between antibody and antigen residue pairs
- **CDR3-antigen PAE**: Mean PAE between CDR3 residues and antigen residues

### 2.4 G+-Specific Metrics (from STEERING_METHODS.md)

These are the metrics specifically defined for evaluating Strategy G+:

1. **Contact score per phase**: Should increase monotonically across phases 1->2->3
   - Requires logging contact scores during diffusion (not available in post-hoc evaluation)
   - Can be approximated by analyzing final contact quality

2. **Structural RMSD vs native**: L-RMSD and I-RMSD (computed by DockQ)

3. **Ensemble diversity**: Pairwise structural variation across the 5 predictions per complex
   - Measured as pairwise RMSD between all 5 models (10 pairs)
   - Report mean and standard deviation of pairwise RMSD

4. **pLDDT trajectory**: pLDDT progression during diffusion
   - Requires logging during inference (not available post-hoc from saved outputs)
   - Approximate via comparing per-residue pLDDT distributions across models

### 2.5 Ensemble Diversity Metrics

For each complex, given 5 predicted models:

| Metric | Definition |
|--------|------------|
| **Pairwise CA-RMSD** | All-atom or CA RMSD between each pair of 5 models (10 pairs). Report mean +/- std |
| **Interface diversity** | Pairwise I-RMSD between model pairs (using interface residues only) |
| **CDR3 diversity** | Pairwise RMSD of CDR3 backbone atoms across models |
| **Contact diversity** | Variance in fnat across the 5 models (high = inconsistent, low = converged) |

---

## 3. Model Selection Strategy

Each method produces 5 models per complex. We need a selection strategy:

### 3.1 Best-of-5 (Primary)

**Selection criterion**: Select the model with the highest `confidence_score` from the confidence JSON.

**Rationale**: This is the standard in the field (AlphaFold-Multimer, Chai-1, Boltz-1 all report best-of-N). It reflects practical usage where users pick the most confident prediction.

**Alternative selection criteria** (to compare):
- Best by `iptm` (interface-focused selection)
- Best by `complex_iplddt` (interface pLDDT-focused)
- Best by DockQ (oracle/upper bound - requires ground truth, not a real predictor)

### 3.2 All-5 Analysis (Secondary)

Report statistics across all 5 models to assess consistency:
- Mean and standard deviation of DockQ across 5 models
- Fraction of models achieving "Acceptable" or better
- Ensemble diversity (see Section 2.5)

### 3.3 Oracle Selection (Diagnostic)

Select the model with the best DockQ score (using ground truth). This provides the **upper bound** of what the method can achieve if we had a perfect model selector. The gap between confidence-based selection and oracle selection reveals how well the confidence score predicts actual quality.

---

## 4. Evaluation Pipeline

### 4.1 Environment Setup

```bash
conda activate dockq2
# Packages available: DockQ 2.1.3, BioPython 1.86, NumPy 1.26.4, Pandas 3.0.0
```

### 4.2 Pipeline Steps

```
Step 1: Collect metadata
  - Parse all confidence JSONs into a single DataFrame
  - Parse CDR ranges from examples/cdrs.csv
  - Map complex names to ground truth PDB files

Step 2: Run DockQ on all prediction-native pairs
  - For each method x complex x model:
      DockQ prediction native --mapping BC:BC A:A --json result.json
  - Collect: DockQ, fnat, L-RMSD, I-RMSD, CAPRI class

Step 3: Extract per-residue metrics
  - Load pLDDT, PAE, PDE from NPZ files
  - Compute CDR3-specific and interface-specific averages

Step 4: Compute ensemble diversity
  - For each method x complex: pairwise RMSD across 5 models

Step 5: Model selection
  - Apply confidence-based selection (best confidence_score)
  - Apply oracle selection (best DockQ)
  - Record selected model index

Step 6: Aggregate and compare
  - Per-complex comparison table
  - Method-level summary statistics
  - Statistical significance tests
```

### 4.3 File Format Handling

- **Baselines** (B1, B2, B3): PDB format -> DockQ handles directly
- **G+**: CIF format -> DockQ v2 handles CIF natively (BioPython MMCIFParser)
- **Ground truth**: PDB format

If CIF parsing causes issues, fallback conversion:
```python
from Bio.PDB import MMCIFParser, PDBIO
parser = MMCIFParser(QUIET=True)
structure = parser.get_structure('s', 'model.cif')
io = PDBIO()
io.set_structure(structure)
io.save('model.pdb')
```

### 4.4 Baseline 3 Sub-experiments

Baseline 3 (contact restraints) has multiple restraint types per complex:
- `hbond_23`, `hbond_7`, `hydrophobic_5`, `salt_bridge_3`

**Strategy**: Evaluate each sub-experiment separately, then report:
- Best across all restraint types (shows ceiling of contact restraint approach)
- Per-restraint-type averages (shows which restraint type works best)

---

## 5. Statistical Analysis

### 5.1 Per-Complex Comparison

For each complex, report a row with:

| Complex | B1_DockQ | B2_DockQ | B3_DockQ | G+_DockQ | B1_fnat | ... | B1_LRMSD | ... |
|---------|----------|----------|----------|----------|---------|-----|----------|-----|

### 5.2 Aggregate Statistics

| Metric | B1 | B2 | B3 | G+ |
|--------|----|----|----|----|
| Mean DockQ (best-of-5) | | | | |
| Median DockQ (best-of-5) | | | | |
| % Acceptable+ (DockQ >= 0.23) | | | | |
| % Medium+ (DockQ >= 0.49) | | | | |
| % High (DockQ >= 0.80) | | | | |
| Mean fnat | | | | |
| Mean L-RMSD (A) | | | | |
| Mean I-RMSD (A) | | | | |
| Mean iptm | | | | |
| Mean complex_iplddt | | | | |
| Mean ensemble diversity (CA-RMSD) | | | | |

### 5.3 Statistical Tests

- **Paired Wilcoxon signed-rank test**: Compare G+ vs each baseline (paired by complex)
- **Effect size**: Report median DockQ improvement and 95% confidence interval
- **Per-complex delta**: DockQ(G+) - DockQ(Baseline) for each complex

### 5.4 Visualizations

1. **DockQ bar plot**: One bar per method, grouped by complex (or sorted by G+ improvement)
2. **DockQ scatter**: G+ DockQ vs B1 DockQ (each point = one complex), with y=x diagonal
3. **CAPRI classification stacked bar**: % Incorrect/Acceptable/Medium/High per method
4. **Confidence vs DockQ correlation**: Scatter of confidence_score vs DockQ to assess reliability of model selection
5. **Ensemble diversity comparison**: Box plot of pairwise RMSD per method
6. **Per-residue pLDDT heatmap**: pLDDT profiles aligned by sequence for representative complexes
7. **CDR3-specific analysis**: CDR3 pLDDT and CDR3-antigen PAE distributions per method

---

## 6. Implementation Notes

### 6.1 Script Structure

```
evaluation/
├── run_dockq.py              # Run DockQ on all prediction-native pairs
├── extract_confidence.py     # Parse confidence JSONs into DataFrame
├── extract_perres.py         # Extract pLDDT/PAE/PDE per-residue metrics
├── compute_diversity.py      # Ensemble pairwise RMSD
├── aggregate_results.py      # Combine all metrics, model selection, statistics
├── plot_results.py           # Generate all figures
└── results/                  # Output directory
    ├── dockq_results.csv
    ├── confidence_metrics.csv
    ├── perres_metrics.csv
    ├── diversity_metrics.csv
    ├── summary_table.csv
    └── figures/
```

### 6.2 DockQ Batch Script (Pseudocode)

```python
import subprocess, json, pandas as pd
from pathlib import Path

METHODS = {
    'B1': 'antigen_cut',
    'B2': 'antigen_cut_vhvl_msa_pocket_ab_boltz_post_2023_omm',
    'B3': 'antigen_cut_contact_restraints',
    'Gplus': 'new_feature_cdr3_beta',
}

NATIVE_DIR = Path('boltz_evals/pdb_minimized/')
PRED_DIR = Path('predictions_examples/')

results = []
for method_name, method_dir in METHODS.items():
    method_path = PRED_DIR / method_dir
    for complex_dir in method_path.glob('boltz_results_*'):
        complex_name = extract_complex_name(complex_dir.name)
        native_pdb = NATIVE_DIR / f'{complex_name}.pdb'

        pred_dir = complex_dir / 'predictions' / get_pred_subfolder(complex_dir)
        for model_file in sorted(pred_dir.glob('*_model_*.pdb')) + sorted(pred_dir.glob('*_model_*.cif')):
            model_idx = extract_model_index(model_file.name)

            # Run DockQ
            result = subprocess.run(
                ['DockQ', str(model_file), str(native_pdb),
                 '--mapping', 'BC:BC', 'A:A', '--json', '/tmp/dockq_out.json'],
                capture_output=True, text=True
            )
            dockq_data = json.load(open('/tmp/dockq_out.json'))

            results.append({
                'method': method_name,
                'complex': complex_name,
                'model': model_idx,
                'DockQ': dockq_data['DockQ'],
                'fnat': dockq_data['fnat'],
                'LRMSD': dockq_data['LRMS'],
                'iRMSD': dockq_data['iRMS'],
                'CAPRI': dockq_data['CAPRI'],
            })

df = pd.DataFrame(results)
df.to_csv('results/dockq_results.csv', index=False)
```

### 6.3 Ensemble Diversity (Pseudocode)

```python
from Bio.PDB import PDBParser, Superimposer
import numpy as np
from itertools import combinations

def compute_pairwise_rmsd(model_files: list, chain_ids='BC'):
    """Compute pairwise CA-RMSD between all model pairs."""
    parser = PDBParser(QUIET=True)
    structures = [parser.get_structure(f'm{i}', f) for i, f in enumerate(model_files)]

    rmsds = []
    for (i, s1), (j, s2) in combinations(enumerate(structures), 2):
        atoms1 = [a for c in chain_ids for a in s1[0][c].get_atoms() if a.name == 'CA']
        atoms2 = [a for c in chain_ids for a in s2[0][c].get_atoms() if a.name == 'CA']

        sup = Superimposer()
        sup.set_atoms(atoms1, atoms2)
        rmsds.append(sup.rms)

    return np.mean(rmsds), np.std(rmsds)
```

### 6.4 Confidence Extraction (Pseudocode)

```python
import json
from pathlib import Path

def extract_confidence_metrics(json_path):
    with open(json_path) as f:
        data = json.load(f)
    return {
        'confidence_score': data['confidence_score'],
        'ptm': data['ptm'],
        'iptm': data['iptm'],
        'protein_iptm': data['protein_iptm'],
        'complex_plddt': data['complex_plddt'],
        'complex_iplddt': data['complex_iplddt'],
        'complex_ipde': data['complex_ipde'],
        'iptm_AB': data['pair_chains_iptm']['0']['1'],  # antigen-heavy
        'iptm_AC': data['pair_chains_iptm']['0']['2'],  # antigen-light
        'iptm_BC': data['pair_chains_iptm']['1']['2'],  # heavy-light
    }
```

---

## 7. Current Data Availability

Based on inspection of the local prediction folders:

| Method | Complexes Available | Models/Complex | Format | Prediction Files Present |
|--------|-------------------|----------------|--------|-------------------------|
| B1 (antigen_cut) | 2 (7TRH_HBG, 7TRI_ZYB) | 5 | PDB | Yes |
| B2 (pocket) | 3 sub-experiments (7TRH only) | 5 | PDB | No (empty) |
| B3 (contact restraints) | 4 sub-experiments (7TRH only) | 5 | PDB | No (empty) |
| G+ (cdr3_beta) | 1 (7TRH_HBG) | 1 | CIF | Yes |
| Ground truth | 48 | - | PDB | Yes |

**Note**: The local folders contain partial/example data. Full predictions for all 48 complexes are expected to be generated on a compute server. The evaluation pipeline should be designed to handle 48 complexes x 5 models x 4 methods = 960 DockQ evaluations (plus sub-experiments for B2/B3).

---

## 8. Execution Plan

### Phase 1: Prototype on Available Data

1. Run DockQ on B1 vs ground truth for 7TRH_HBG (5 models) and 7TRI_ZYB (5 models)
2. Run DockQ on G+ vs ground truth for 7TRH_HBG (1 model)
3. Extract confidence metrics for all available predictions
4. Validate the pipeline works end-to-end
5. Debug any chain mapping or format issues

### Phase 2: Full Evaluation (after predictions are generated)

1. Generate all predictions (48 complexes x 5 models x 4 methods)
2. Run DockQ batch evaluation
3. Extract all confidence + per-residue metrics
4. Compute ensemble diversity
5. Apply model selection (confidence-based and oracle)
6. Aggregate results and compute statistics

### Phase 3: Analysis and Reporting

1. Generate all plots and tables
2. Run statistical tests (Wilcoxon paired)
3. Identify best/worst cases for G+ vs baselines
4. Analyze CDR3-specific improvements
5. Write results summary

---

## 9. Expected Outcomes

Based on the design of Strategy G+, we hypothesize:

1. **G+ should improve DockQ over B1** (no steering): The adaptive phase scheduling with contact-triggered transitions and CDR3 beta-scaling should produce better antibody-antigen interfaces

2. **G+ should be competitive with B2/B3**: Contact restraints (B2, B3) use explicit structural knowledge; G+ uses learned priors via beta-scaling. G+ may win on diversity while B2/B3 may win on specific contacts

3. **G+ should maintain or improve ensemble diversity**: The exploratory Phase 1 (negative beta) encourages diverse CDR3 conformations before refinement

4. **Confidence scores should correlate with DockQ**: If model selection works, confidence-selected best-of-5 should approach oracle best-of-5

---

## 10. Questions to Resolve Before Implementation

1. **DockQ chain mapping**: Verify that `--mapping BC:BC A:A` works correctly (treat antibody H+L as receptor, antigen as ligand). Test on one example first.

2. **CIF compatibility**: Confirm DockQ v2 handles the Boltz2 CIF output format without errors. The CIF uses ModelCIF dictionary (mmcif_ma.dic) which may differ from standard PDB mmCIF.

3. **B3 sub-experiments**: How to fairly compare multiple restraint types? Report best-of-restraint-types or average?

4. **Missing predictions**: B2 and B3 prediction folders are currently empty. Confirm these will be populated before full evaluation.

5. **Sequence matching**: DockQ uses sequence alignment to match residues between model and native. Verify that the sequences in predictions match the ground truth (same truncation, same residue numbering).

6. **G+ model count**: Currently only 1 model available for G+. The full evaluation needs 5 models per complex for fair comparison. Confirm `--diffusion_samples 5` will be used.
