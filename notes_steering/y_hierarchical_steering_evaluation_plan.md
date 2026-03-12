# Evaluation Plan: Idea Y — Hierarchical Steering

## 1. Feature Goals and What We Are Evaluating

Hierarchical Steering combines **embedding-space steering** (early diffusion, robust for coarse topology) with **coordinate-space potentials** (late diffusion, precise for interface geometry). The feature targets:

1. **Antibody-antigen docking quality** — correct relative orientation of Ab/Ag
2. **Interface accuracy** — native contacts recovered, interface RMSD
3. **CDR loop conformation** — especially CDR-H3, the primary binding determinant
4. **Epitope identification** — does the predicted binding site match the true epitope?
5. **Robustness** — consistent improvement across diverse complexes, not just cherry-picked cases

The evaluation must demonstrate that the hierarchical two-stage approach outperforms each single-stage method alone.

**From STEERING_METHODS.md evaluation criteria for Idea Y**:
- Quality metrics for each stage (stage 1 robustness, stage 2 precision)
- Transition smoothness (no discontinuities in quality metrics)
- Final structure quality vs. single-stage methods
- Computational cost (should be < 2× single method)

---

## 2. Dataset

- **48 antibody-antigen complexes** from `examples/cdrs.csv` (48 data rows + 1 header = 49 lines)
- **Ground truth structures**: 48 PDB files at `/mnt/c/Users/octav/Documents/claude/proteinEBM/boltz/pdb_minimized/` (minimized crystal structures)
- **Chain convention**: A = antigen, B = heavy chain, C = light chain. Both predictions and ground truths use chains A, B, C.
- **CDR definitions**: per-complex from `examples/cdrs.csv` (1-indexed residue positions for H1, H2, H3, L1, L2, L3)

### Ground truth files
All 48 PDB files:
```
7TRH_HBG.pdb  7TRI_ZYB.pdb  7WT9_HAE.pdb  7Y0O_HLA.pdb  7ZOZ_HLA.pdb
8BLQ_ABD.pdb  8BLQ_ECD.pdb  8BYU_HLA.pdb  8CDD_EDB.pdb  8CDE_DCB.pdb
8CXC_HLM.pdb  8CYH_HLM.pdb  8DE3_BCA.pdb  8DTK_CBA.pdb  8E2U_HLA.pdb
8EAY_HLA.pdb  8EQ6_HLA.pdb  8EZ3_HLA.pdb  8EZ7_HLA.pdb  8EZ8_HLA.pdb
8FAH_HLA.pdb  8FGX_BAC.pdb  8FXB_HLE.pdb  8GH4_HLE.pdb  8GHP_HLA.pdb
8GP5_EFX.pdb  8GQ1_HLC.pdb  8HES_HLC.pdb  8HGM_CDB.pdb  8HLB_BCA.pdb
8HLB_DEA.pdb  8IDN_HLA.pdb  8IV4_ABG.pdb  8IV4_HLG.pdb  8IV5_ABG.pdb
8IVA_CEG.pdb  8IVX_HLA.pdb  8IX3_HLG.pdb  8J1T_HKF.pdb  8OL9_BAH.pdb
8OXW_BCA.pdb  8OXX_BCA.pdb  8PE9_HLA.pdb  8SLB_HLA.pdb  8T9Z_HLA.pdb
8TFR_ABC.pdb  8TRS_AGD.pdb  8X0T_HLA.pdb
```

---

## 3. Methods to Compare (4 conditions)

| Label | Folder pattern | Description | Settings/complex |
|-------|---------------|-------------|-----------------|
| **Baseline (antigen_cut)** | `predictions_examples/antigen_cut/boltz_results_<NAME>/` | Boltz2 vanilla — no steering, antigen-cut MSA | 1 |
| **Contact Restraints** | `predictions_examples/antigen_cut_contact_restraints/boltz_results_restraint_<NAME>_<type>_<N>/` | Boltz2 + contact restraint potentials (hbond, hydrophobic, salt bridge) | Multiple |
| **Pocket-guided** | `predictions_examples/antigen_cut_vhvl_msa_pocket_ab_boltz_post_2023_omm/boltz_results_restraint_to_A_<NAME>_<chain>_<res>_<idx>/` | Boltz2 + VH/VL MSA + pocket-based Ab restraints | Multiple |
| **Hierarchical Steering (Y)** | `predictions_examples/new_feature_hierarchical/boltz_results_<NAME>_hierarchical/` | Boltz2 + hierarchical embedding+coordinate steering | 1 |

### Prediction folder structure (per complex per method)
```
boltz_results_<NAME>/
├── lightning_logs/
├── msa/
├── predictions/
│   └── <NAME>/
│       ├── <NAME>_model_0.pdb      # (or .cif for new features)
│       ├── <NAME>_model_1.pdb
│       ├── <NAME>_model_2.pdb
│       ├── <NAME>_model_3.pdb
│       ├── <NAME>_model_4.pdb
│       ├── confidence_<NAME>_model_0.json
│       ├── confidence_<NAME>_model_1.json
│       ├── ... (5 pae, 5 pde, 5 plddt .npz files)
└── processed/
```

### Handling multiple settings per complex
For contact restraints and pocket-guided baselines (which have multiple runs per complex, e.g. different restraint types):
- **Best-per-complex (oracle)**: select the setting that yields the best DockQ — upper bound
- **Average-per-complex**: average DockQ across all settings — expected performance without oracle
- This ensures fair comparison against single-setting methods.

---

## 4. Prediction Ensemble: Model Selection Strategy

Each run produces **5 independent structure predictions** (model_0 through model_4), each with:
- Structure file (.pdb or .cif)
- `confidence_<NAME>_model_<N>.json` with confidence_score, ptm, iptm, plddt, etc.
- `plddt_<NAME>_model_<N>.npz` — per-residue pLDDT array, shape (N_atoms,)
- `pae_<NAME>_model_<N>.npz` — predicted aligned error matrix, shape (N_tokens, N_tokens)
- `pde_<NAME>_model_<N>.npz` — predicted distance error

### Model selection: Report three variants
1. **Top-1 (confidence-ranked)**: Select model with highest `confidence_score` from JSON. This is the standard approach used by Boltz2 (`aggregate_evals.py:compute_boltz_metrics`), AF2/AF3, and the field. **Primary reported metric.**
2. **Oracle (best DockQ)**: Select model with best actual DockQ. Upper bound on method performance.
3. **Average**: Mean ± std across all 5 models. Measures consistency/robustness.

### Confidence score details
The `confidence_score` in the JSON is a weighted combination:
```
confidence_score ≈ 0.2 × ptm + 0.8 × iptm
```
This weighting (from AlphaFold-Multimer) emphasizes interface quality, making it ideal for ranking docking models.

### Existing evaluation infrastructure
The codebase already has:
- `scripts/eval/run_evals.py` — runs OpenStructure Docker (computes lddt, bb-lddt, qs-score, dockq, ics, ips, tm-score). **Requires Docker with OpenStructure image** — may not be available on all machines.
- `scripts/eval/aggregate_evals.py` — aggregates results: top-1/oracle/avg selection, bootstrap CIs, plotting. Computes `dockq_>0.23` and `dockq_>0.49` success rates.

**For our evaluation**: We use DockQ v2 directly (no Docker needed) for interface metrics, plus custom BioPython scripts for CDR RMSD and epitope analysis. The aggregation logic follows the same top-1/oracle pattern as `aggregate_evals.py`.

---

## 5. Evaluation Metrics

### 5.1 DockQ v2 — Primary Metric (Interface Quality)

**Tool**: DockQ v2.1.3 (`/home/oc/anaconda3/envs/dockq2/bin/DockQ`)
**Environment**: `/home/oc/anaconda3/envs/dockq2/`

DockQ is the standard metric for protein-protein docking quality. It combines three CAPRI sub-metrics into a single score [0, 1]:

```
DockQ = (fnat + 1/(1+(iRMSD/1.5)²) + 1/(1+(LRMSD/8.5)²)) / 3
```

| Sub-metric | What it measures |
|------------|-----------------|
| **fnat** | Fraction of native contacts recovered (heavy atoms ≤ 5Å across interface) |
| **iRMSD** | Interface RMSD — backbone RMSD of interface residues (within 10Å of partner in native) after interface superposition |
| **LRMSD** | Ligand RMSD — backbone RMSD of the "ligand" (smaller chain) after superimposing the "receptor" (larger chain). For Ab-Ag: antigen is receptor, Fv is ligand |
| **F1** | Harmonic mean of contact precision and recall (addresses false positives via fnonnat) |
| **fnonnat** | Fraction of predicted contacts that are non-native (false positive rate) |

**Quality classification (CAPRI)**:

| DockQ range | Quality | fnat | iRMSD | LRMSD |
|-------------|---------|------|-------|-------|
| 0.00–0.23 | Incorrect | < 0.1 | > 4 Å | > 10 Å |
| 0.23–0.49 | Acceptable | ≥ 0.1 | ≤ 4 Å | ≤ 10 Å |
| 0.49–0.80 | Medium | ≥ 0.3 | ≤ 2 Å | ≤ 2 Å |
| ≥ 0.80 | High | ≥ 0.5 | ≤ 1 Å | ≤ 1 Å |

**What to compute per complex**:
- **Total DockQ**: average over all 3 interfaces (A-B, A-C, B-C)
- **Ab-Ag DockQ**: average of only antibody-antigen interfaces (A-B and A-C), excluding the VH-VL interface (B-C). This is the most relevant metric for our steering feature since it targets the binding interface.
- **Per-interface breakdown**: DockQ, iRMSD, LRMSD, fnat, F1 for each interface

**Success rates** (following `aggregate_evals.py` convention):
- `dockq_>0.23`: fraction of complexes with acceptable or better docking
- `dockq_>0.49`: fraction of complexes with medium or better docking

**Command**:
```bash
/home/oc/anaconda3/envs/dockq2/bin/DockQ <model> <native.pdb> --json <output.json>
```
DockQ v2 handles both PDB and CIF inputs. It auto-detects chain mapping. Use `--mapping ABC:ABC` when chain IDs are identical in model and native.

### 5.2 CDR Loop RMSD — CDR Conformation Quality

Evaluate CDR loop accuracy separately — CDR-H3 is the most variable and the primary binding determinant.

**Procedure**:
1. Parse model and native structures (BioPython `PDBParser`/`MMCIFParser`)
2. Identify antibody chains (B = heavy, C = light) and CDR residues from `examples/cdrs.csv`
3. **Align on antibody framework Cα atoms** (all VH+VL residues minus CDR residues). This isolates CDR loop errors from global docking errors.
4. After framework superposition, compute **Cα-RMSD for each CDR loop** (H1, H2, H3, L1, L2, L3)

**Why framework alignment**: Standard in the field (ABodyBuilder, ImmuneBuilder, IgFold benchmarks). Aligning on the conserved framework ensures we measure CDR conformation quality independently of docking pose.

**Reference values** (from literature):
- Good prediction: CDR-H3 RMSD < 1.5 Å
- Acceptable: CDR-H3 RMSD < 2.0 Å
- CDR-H3 is always the worst; H1/H2/L1/L2/L3 are typically < 1.5 Å

**Implementation**: BioPython `Superimposer` in `/home/oc/anaconda3/envs/boltz/`

```python
from Bio.PDB import PDBParser, MMCIFParser, Superimposer
import numpy as np

def compute_cdr_rmsd(model_path, native_path, cdr_indices, heavy_chain='B', light_chain='C'):
    """
    cdr_indices: dict with keys 'cdr1_h', 'cdr2_h', 'cdr3_h', 'cdr1_l', 'cdr2_l', 'cdr3_l'
                 values are lists of 1-indexed residue positions
    """
    parser = PDBParser(QUIET=True) if model_path.endswith('.pdb') else MMCIFParser(QUIET=True)
    model_struct = parser.get_structure('model', model_path)
    native_struct = PDBParser(QUIET=True).get_structure('native', native_path)

    # Collect all CDR residue indices (1-indexed)
    all_cdr_h = set(cdr_indices['cdr1_h'] + cdr_indices['cdr2_h'] + cdr_indices['cdr3_h'])
    all_cdr_l = set(cdr_indices['cdr1_l'] + cdr_indices['cdr2_l'] + cdr_indices['cdr3_l'])

    # Extract framework CA atoms (non-CDR residues in Ab chains)
    framework_model_atoms = []
    framework_native_atoms = []
    for chain_id, cdr_set in [(heavy_chain, all_cdr_h), (light_chain, all_cdr_l)]:
        model_chain = model_struct[0][chain_id]
        native_chain = native_struct[0][chain_id]
        for res_m, res_n in zip(model_chain.get_residues(), native_chain.get_residues()):
            res_idx = res_m.get_id()[1]
            if res_idx not in cdr_set and 'CA' in res_m and 'CA' in res_n:
                framework_model_atoms.append(res_m['CA'])
                framework_native_atoms.append(res_n['CA'])

    # Superimpose on framework
    sup = Superimposer()
    sup.set_atoms(framework_native_atoms, framework_model_atoms)
    sup.apply(model_struct.get_atoms())

    # Compute per-CDR RMSD
    results = {}
    for cdr_name, chain_id in [('cdr1_h', heavy_chain), ('cdr2_h', heavy_chain),
                                ('cdr3_h', heavy_chain), ('cdr1_l', light_chain),
                                ('cdr2_l', light_chain), ('cdr3_l', light_chain)]:
        indices = cdr_indices[cdr_name]
        model_chain = model_struct[0][chain_id]
        native_chain = native_struct[0][chain_id]
        diffs = []
        for idx in indices:
            try:
                m_ca = model_chain[(' ', idx, ' ')]['CA'].get_vector().get_array()
                n_ca = native_chain[(' ', idx, ' ')]['CA'].get_vector().get_array()
                diffs.append(np.sum((m_ca - n_ca) ** 2))
            except KeyError:
                continue
        if diffs:
            results[cdr_name] = np.sqrt(np.mean(diffs))
    return results
```

### 5.3 Epitope Prediction — Contact-Based

Evaluate whether the model identifies the correct binding site (epitope) on the antigen.

**Definition**: An antigen residue is an **epitope residue** if any of its heavy atoms are within a distance threshold of any antibody heavy atom.

**Two thresholds**:
- **5.0 Å** (heavy atom): Standard CAPRI/DockQ contact definition — strict
- **8.0 Å** (Cα-Cα): Coarser epitope definition, consistent with `contact_threshold` in our feature config

**Procedure**:
1. Extract epitope residues from native structure (antigen residues contacting Ab)
2. Extract epitope residues from predicted structure
3. Compute per-antigen-residue binary classification metrics

**Metrics**:
| Metric | What it measures |
|--------|-----------------|
| **Precision** | Of predicted epitope residues, fraction that are true epitope |
| **Recall** | Of true epitope residues, fraction that are predicted |
| **F1** | Harmonic mean of precision and recall |
| **MCC** | Matthews Correlation Coefficient — handles class imbalance (epitope is small fraction of antigen surface) |

**Note**: DockQ already computes fnat (recall) and F1 for contacts at 5Å. Our epitope analysis adds the residue-level view (which antigen residues, not which residue pairs) and the 8Å threshold.

**Environment**: `/home/oc/anaconda3/envs/boltz/` (biopython + numpy + scipy)

### 5.4 Boltz Confidence Metrics

Extract from the confidence JSON files — no external tools needed:

| Metric | Description | Relevance |
|--------|-------------|-----------|
| `confidence_score` | Combined Boltz2 confidence (≈0.2×ptm + 0.8×iptm) | Model ranking |
| `iptm` | Interface predicted TM-score | Docking quality self-assessment |
| `ptm` | Predicted TM-score (overall fold quality) | Fold quality |
| `complex_plddt` | Average pLDDT across all residues | Per-residue confidence |
| `complex_iplddt` | Average pLDDT at interface residues | Interface confidence |
| `pair_chains_iptm` | Pairwise iptm between chain pairs | Per-interface confidence |

**Purpose**:
- Verify confidence correlates with actual quality (DockQ). Plot iptm vs DockQ scatter.
- Check that steering doesn't inflate confidence without improving quality (calibration check).
- The `pair_chains_iptm["0"]["1"]` and `pair_chains_iptm["0"]["2"]` values give Ag-VH and Ag-VL interface confidence specifically.

### 5.5 Computational Cost

**Criterion from STEERING_METHODS.md**: Hierarchical steering should cost **< 2× single-method cost**.

**What adds cost**:
- Conditioning is computed twice (steered + neutral) — one extra forward pass through DiffusionConditioning
- Blending is just tensor interpolation (negligible cost)
- Coordinate potentials in late stages (same cost as single-stage coordinate steering)

**Measurement**: Time the `boltz predict` command for baseline vs hierarchical on the same complex. Report wall-clock time ratio.

---

## 6. Postprocessing Requirements

### 6.1 DockQ — No manual alignment needed
DockQ v2 internally handles:
- Sequence alignment between model and native
- Receptor/ligand designation (larger chain = receptor)
- Superposition for iRMSD and LRMSD computation
- Contact identification at 5Å threshold

Just provide model and native files. Both PDB and CIF are supported.

### 6.2 CDR RMSD — Framework alignment required
1. Parse structures with BioPython (`PDBParser` for .pdb, `MMCIFParser` for .cif)
2. Extract Cα atoms for framework residues (Ab chains minus CDR residues)
3. Superimpose model framework onto native framework
4. Compute per-CDR Cα RMSD after superposition (do NOT re-align per CDR)

### 6.3 Epitope prediction — No alignment needed
Contacts are defined by inter-atomic distances within the complex (translation/rotation invariant). Just extract contacts from model and native independently.

### 6.4 Format handling
- Baselines: `.pdb` files — BioPython `PDBParser`, DockQ native support
- Hierarchical steering: `.cif` (mmCIF) files — BioPython `MMCIFParser`, DockQ v2 supports natively
- Ground truths: `.pdb` files
- The evaluation script should auto-detect format by file extension.

### 6.5 Chain mapping
Model and native both use chains A, B, C. DockQ auto-mapping works. For safety, always use `--mapping ABC:ABC`.

---

## 7. Statistical Analysis

### 7.1 Paired comparison
All 48 complexes are evaluated under each method → paired comparison.
- **Wilcoxon signed-rank test** (non-parametric, paired): hierarchical vs. each baseline
- Report p-values for: total DockQ, Ab-Ag DockQ, CDR-H3 RMSD, epitope F1

### 7.2 Bootstrap confidence intervals
Following `aggregate_evals.py`: 1000 bootstrap resamples for 95% CI on the mean.
```python
def bootstrap_ci(values, n_boot=1000, alpha=0.05):
    boot_means = [np.random.choice(values, len(values), replace=True).mean() for _ in range(n_boot)]
    return np.mean(values), np.percentile(boot_means, 100*alpha/2), np.percentile(boot_means, 100*(1-alpha/2))
```

### 7.3 Summary statistics per method
- Mean ± std across complexes
- Median (robust to outliers)
- Win/tie/loss count (hierarchical vs. each baseline per complex)
- Success rates: `dockq_>0.23`, `dockq_>0.49`

### 7.4 Visualization
1. **Bar plot with bootstrap CIs**: Mean metric per method (DockQ, fnat, success rates) — following `aggregate_evals.py` plotting style
2. **Box plot**: Distribution of DockQ across 48 complexes per method
3. **Paired scatter plot**: Hierarchical DockQ (y) vs. Baseline DockQ (x) per complex, with y=x diagonal
4. **Per-CDR bar plot**: Mean RMSD per CDR loop (H1, H2, H3, L1, L2, L3) across methods
5. **Confidence calibration scatter**: iptm vs. DockQ, colored by method

---

## 8. Execution Plan

### Step 1: Generate predictions
Run hierarchical steering for all 48 complexes:
```bash
/home/oc/anaconda3/envs/boltz/bin/python -m boltz.main predict \
    examples/hierarchical_steering/<NAME>_hierarchical.yml \
    --use_potentials --hierarchical_steering \
    --devices 1 --diffusion_samples 5 \
    --out_dir predictions_examples/new_feature_hierarchical/
```
48 YAML configs are already prepared in `examples/hierarchical_steering/`.
Use `--diffusion_samples 5` for fair comparison with baselines (which have 5 models each).

### Step 2: Model selection
For each complex × method:
1. Locate prediction folder and list model files
2. Load all 5 `confidence_*.json` files
3. Select `top1_idx = argmax(confidence_score)`
4. Record confidence metrics for all models

### Step 3: DockQ evaluation
For each complex × method × model:
```bash
/home/oc/anaconda3/envs/dockq2/bin/DockQ \
    <prediction_model.pdb_or_cif> \
    /mnt/c/Users/octav/Documents/claude/proteinEBM/boltz/pdb_minimized/<NAME>.pdb \
    --json evaluation_results/dockq/<METHOD>/<NAME>_model_<N>.json
```
Parse all JSON results → compute top-1, oracle, and average per complex.

### Step 4: CDR RMSD evaluation
For each complex × method (top-1 model only):
- Parse model and native with BioPython
- Load CDR definitions from `examples/cdrs.csv`
- Framework-align and compute per-CDR Cα RMSD
- Environment: `/home/oc/anaconda3/envs/boltz/`

### Step 5: Epitope evaluation
For each complex × method (top-1 model only):
- Extract antigen contact residues at 5Å and 8Å
- Compare model vs native contacts
- Compute precision, recall, F1, MCC

### Step 6: Aggregate and report
- Compute summary statistics per method
- Run Wilcoxon tests and bootstrap CIs
- Generate plots
- Write summary CSV and report

---

## 9. Environments and Dependencies

| Task | Python | Key packages |
|------|--------|-------------|
| DockQ evaluation | `/home/oc/anaconda3/envs/dockq2/bin/python` | DockQ 2.1.3, biopython 1.86, numpy 1.26, pandas 3.0, scipy 1.17 |
| CDR RMSD + epitope analysis | `/home/oc/anaconda3/envs/boltz/bin/python` | biopython 1.84, numpy 1.26, pandas 2.3, scipy 1.13, scikit-learn 1.6 |
| Running predictions | `/home/oc/anaconda3/envs/boltz/bin/python` | boltz 2.2.1 (editable install) |
| Plotting | Either env | matplotlib (install if missing: `pip install matplotlib`) |

**No new environment needed.** Both `dockq2` and `boltz` environments have all required packages.

---

## 10. Output Files

```
evaluation_results/
├── dockq/                              # Raw DockQ JSON outputs
│   ├── baseline/<NAME>_model_<N>.json
│   ├── contact_restraints/<NAME>_<setting>_model_<N>.json
│   ├── pocket_guided/<NAME>_<setting>_model_<N>.json
│   └── hierarchical/<NAME>_model_<N>.json
├── per_complex_dockq.csv               # DockQ per complex × method (top1, oracle, avg)
├── per_complex_dockq_ab_ag.csv         # Ab-Ag interface DockQ only
├── per_complex_cdr_rmsd.csv            # CDR RMSD per complex × method × CDR loop
├── per_complex_epitope.csv             # Epitope P/R/F1/MCC per complex × method
├── per_complex_confidence.csv          # Confidence metrics per complex × method
├── summary_table.csv                   # Aggregated metrics per method (mean, median, CI)
├── statistical_tests.csv               # Wilcoxon p-values
├── figures/
│   ├── dockq_barplot.png               # Mean DockQ with bootstrap CI per method
│   ├── dockq_boxplot.png               # DockQ distribution per method
│   ├── dockq_scatter_hierarchical_vs_baseline.png
│   ├── cdr_rmsd_barplot.png            # Per-CDR RMSD per method
│   ├── epitope_f1_barplot.png          # Epitope F1 per method
│   └── confidence_vs_dockq.png         # Calibration check
└── evaluation_report.md                # Human-readable summary
```

---

## 11. Key Considerations

1. **Chain mapping**: Both model and native use chains A (antigen), B (heavy), C (light). DockQ auto-mapping handles this. Verify with `--mapping ABC:ABC`.

2. **CIF vs PDB format**: Hierarchical steering outputs CIF; baselines output PDB. DockQ v2 and BioPython handle both. Auto-detect by extension.

3. **Multiple settings per baseline**: Contact restraints have 4 settings per complex (hbond_N, hydrophobic_N, salt_bridge_N). Pocket-guided has 3 per complex (different chain/residue targets). Report both best-of-settings (oracle) and average-of-settings.

4. **5 models per complex**: Run hierarchical steering with `--diffusion_samples 5` for fair comparison. Baselines already have 5 models each.

5. **Computational cost**: Hierarchical computes DiffusionConditioning twice. Expected overhead ~1.5×. Measure wall-clock time for ≥3 complexes and report mean ratio.

6. **CDR definitions are per-complex**: Each antibody has different CDR boundaries. The indices in `examples/cdrs.csv` are 1-indexed. The CSV columns are: `complex, heavy, light, cdr1_h, cdr2_h, cdr3_h, cdr1_l, cdr2_l, cdr3_l`. CDR index lists are Python literal format (e.g., `"[26, 27, 28, 29]"`).

7. **Existing scripts**: `scripts/eval/run_evals.py` (Docker/OpenStructure, not needed here) and `scripts/eval/aggregate_evals.py` (aggregation logic, plotting, bootstrap CIs — reuse patterns from this script).

8. **What constitutes a "win"**: For DockQ, iRMSD, LRMSD, CDR RMSD: lower is better (except DockQ where higher is better). For fnat, F1, epitope metrics: higher is better. Ensure sign conventions are correct in win/loss counts.
