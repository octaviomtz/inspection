# Evaluation Plan: Strategy L+ — CDR3 Beta-Scaling v2

## 1. Feature Goals and Evaluation Hypotheses

Strategy L+ v2 builds on the original CDR3 beta-scaling (L) with four targeted improvements. Before defining metrics we need to understand what each improvement is *expected* to change:

| Sub-improvement | Expected effect | How to measure |
|---|---|---|
| **L+.1** Reduce beta 0.3→0.12 | Reduce CDR3 RMSD inflation caused by over-steering; preserve diversity gain | CDR3 backbone RMSD (lower = better geometry) |
| **L+.2** Asymmetric H3/L3 betas (H3=0.15, L3=0.05) | H3 gets pushed more, L3 is less distorted | CDR-H3 RMSD and CDR-L3 RMSD separately |
| **L+.3** CDR3-antigen interface scaling | Attract CDR3 loops toward the correct epitope; improve interface quality | DockQ, interface RMSD (iRMSD), ligand RMSD (Lrms) |
| **L+.4** Time-dependent beta schedule | Better convergence: explore early, refine late; less structural noise | CDR3 RMSD, DockQ, pLDDT confidence |

**Primary hypothesis**: L+ v2 improves interface quality (DockQ ↑) and CDR3 geometry (RMSD ↓) simultaneously, without the pTM/ipTM inflation observed in v1 (beta=0.3).

**Secondary hypothesis**: L+ v2 maintains or improves ensemble diversity (useful for downstream selection), with better diversity-quality tradeoff than baselines.

---

## 2. Baselines

Three baselines are compared (matching the existing `predictions_examples/` folder structure):

| Method | Folder | Description |
|---|---|---|
| `baseline_antigen_cut` | `antigen_cut/` | Standard Boltz-2 with antigen MSA only (no steering) |
| `baseline_contact_restraints` | `antigen_cut_contact_restraints/` | Boltz-2 + predicted contact restraints |
| `baseline_pocket` | `antigen_cut_vhvl_msa_pocket_ab_boltz_post_2023_omm/` | Boltz-2 + pocket restraints, VH/VL MSA |
| **L+ v2** | `new_feature_cdr3beta_v2/` | CDR3 beta-scaling v2 (all 4 improvements) |

---

## 3. Metrics

### 3.1 Primary Metric: DockQ (via `conda activate dockq2`)

**DockQ** ([Wallner 2024](https://github.com/wallnerlab/DockQ)) is the standard for protein-protein interface quality. It combines:
- **Fnat**: fraction of native contacts recovered
- **iRMSD**: RMSD of interface residues after interface alignment
- **Lrms**: ligand RMSD after receptor alignment

DockQ score thresholds (CAPRI classification):
- ≥ 0.8: High quality
- 0.6–0.8: Medium
- 0.23–0.6: Acceptable
- < 0.23: Incorrect

**Run on all antibody-antigen interfaces**: evaluate chain pairs (antibody heavy + antigen, antibody light + antigen, heavy+light). The most informative pair is the one between antibody chains (B, C) and antigen (A).

### 3.2 CDR3 Backbone RMSD

After aligning the predicted structure onto the ground truth using the **antibody framework** (Fv region excluding CDR loops), compute Cα RMSD of:
- CDR-H3 residues only (validates L+.2 H3 focus)
- CDR-L3 residues only (validates L+.2 L3 preservation)
- CDR-H3 + CDR-L3 combined (overall CDR3 quality)

Alignment strategy: align on antibody framework Cα atoms (structurally conserved β-sheet regions), then measure CDR displacement. This isolates CDR positioning error from rigid-body pose error.

CDR indices are available in `examples/cdrs.csv`.

### 3.3 CDR3 pLDDT

Extract per-residue pLDDT from the `plddt_*_model_N.npz` files and average over CDR-H3 and CDR-L3 residues. Higher pLDDT = model is more confident about the CDR3 conformation.

### 3.4 Ensemble Diversity

Pairwise CDR3 Cα RMSD among the 5 models of the same complex (framework-aligned). Reports:
- Mean pairwise diversity
- Max pairwise diversity

High diversity = the steering is exploring conformational space effectively.

### 3.5 Boltz Confidence Metrics

From `confidence_*_model_N.json`:
- `confidence_score` (primary model-ranking metric)
- `complex_plddt`
- `ptm`, `iptm`, `protein_iptm`

These are monitored to detect the pTM/ipTM inflation issue that was observed in v1 (beta=0.3). We expect L+ v2 to have lower/more calibrated iptm than v1.

---

## 4. Best Prediction Selection Strategy

Each method produces 5 structure predictions per complex. Two selection strategies should be evaluated:

1. **Best by confidence_score**: Select the model with the highest `confidence_score` per complex. This reflects the model's self-assessment and is the practical selection strategy users would apply.

2. **Oracle (best by DockQ)**: Select the model with the highest DockQ against the ground truth. This is the upper bound and shows how much is left on the table by confidence-based selection.

Report both. The gap between oracle and confidence-based selection reveals how well Boltz-2 self-calibrates under CDR3 steering.

---

## 5. Alignment and Postprocessing

### 5.1 Alignment Protocol

1. Load predicted structure and ground truth structure with **BioPython** (`PDBParser` / `MMCIFParser`)
2. Map chain IDs from prediction to ground truth using sequence alignment (handle potential chain reordering)
3. Extract Cα atoms from antibody framework residues (VH: residues ~1–96 excluding CDR-H1/H2/H3, VL: ~1–87 excluding CDR-L1/L2/L3 per IMGT numbering; or use the flanking residues not in CDRs CSV)
4. Superpose prediction onto ground truth using framework Cα only (BioPython `Superimposer`)
5. Compute CDR3 Cα RMSD after superposition (do NOT re-align on CDRs)

### 5.2 DockQ Postprocessing

DockQ handles its own alignment internally (aligns on the receptor chain before computing Lrms). Pass:
- Predicted structure (converted to PDB if needed)
- Native/ground truth structure from `pdb_minimized/`
- Chain mapping flags so DockQ aligns the correct chain pairs

### 5.3 File Format: CIF → PDB Conversion

**Critical issue**: The baselines produce `.pdb` files but the new feature (run with current Boltz-2) produces `.cif` files. The evaluation script and DockQ both need consistent input.

Two options:
- **Option A (preferred)**: Convert `.cif` → `.pdb` using BioPython's `MMCIF2Dict` + `PDBIO` as a preprocessing step in the evaluation script
- **Option B**: Update DockQ call to pass `.cif` directly (DockQ ≥2.0 supports CIF via `--pdb1` / `--pdb2` flags)

Option A is preferred since BioPython handles it natively and ensures consistent chain labeling.

---

## 6. Reuse of `evaluate_cdr3_beta.py`

**Short answer**: Yes, the same script can be reused with targeted modifications.

### What the existing script already covers (no changes needed):
- DockQ evaluation on all chain pairs
- CDR3 backbone RMSD (framework-aligned)
- CDR3 pLDDT extraction from npz
- Boltz confidence metric extraction
- Ensemble diversity (pairwise CDR3 RMSD)
- Best-model and best-setting selection
- Wilcoxon signed-rank statistical tests
- Box plots, scatter plots, heatmaps
- Markdown report generation

### Required modifications for L+ v2:

| Issue | Change needed |
|---|---|
| **Feature folder name** | Update `"cdr3_beta_scaling"` folder config entry from `new_feature_cdr3_beta` → `new_feature_cdr3beta_v2` (or whatever the full-evaluation folder is named) |
| **CIF file support** | Add CIF → PDB conversion step when loading prediction files (BioPython `MMCIFParser` → `PDBIO`) |
| **Separate H3 vs L3 RMSD** | Add per-CDR breakdown (H3 RMSD and L3 RMSD separately) to validate L+.2 asymmetric scaling |
| **CDR3-antigen iRMSD** | Optionally add iRMSD specifically for CDR3-antigen interface residues (to validate L+.3) — can be extracted from DockQ output |
| **pTM/ipTM comparison** | Already in script; add explicit comparison against v1 (beta=0.3) if v1 results are available |

---

## 7. Conda Environment

| Task | Environment | Package |
|---|---|---|
| DockQ evaluation | `conda activate dockq2` | `DockQ` (Wallner) |
| CIF parsing, structure alignment, RMSD | `conda activate boltz` | `biopython` |
| Data analysis, statistics | `conda activate boltz` | `numpy`, `scipy`, `pandas`, `matplotlib` |

No new environment needed — `boltz` environment covers all Python needs; `dockq2` is called as a subprocess.

---

## 8. Statistical Analysis

- **Wilcoxon signed-rank test** (paired, non-parametric) comparing L+ v2 vs each baseline on:
  - DockQ (per complex, best-by-confidence model)
  - CDR-H3 RMSD
  - CDR-L3 RMSD
- Report: n_wins, n_losses, n_ties, p-value, effect size (median difference)
- Significance threshold: p < 0.05

---

## 9. Output Files

The evaluation script should produce:

```
eval_output_l_plus_v2/
├── results_per_model.csv         # one row per (complex, method, model)
├── results_best_model.csv        # one row per (complex, method) — best by confidence
├── results_best_setting.csv      # one row per (complex, method) — best setting for multi-run methods
├── diversity.csv                 # ensemble diversity per (complex, method)
├── plots/
│   ├── dockq_boxplot.png
│   ├── cdr3_rmsd_boxplot.png
│   ├── cdrh3_rmsd_boxplot.png    # new for v2
│   ├── cdrl3_rmsd_boxplot.png    # new for v2
│   ├── plddt_comparison.png
│   ├── diversity_barplot.png
│   ├── dockq_delta_heatmap.png
│   └── scatter_*.png
└── evaluation_report.md          # summary with tables and statistical tests
```

---

## 10. Summary Recommendation

1. **Reuse `evaluate_cdr3_beta.py`** with the three modifications listed in Section 6
2. Run `conda activate dockq2` for DockQ subprocess calls, `conda activate boltz` for all Python logic
3. Add a CIF→PDB conversion wrapper so both old (`.pdb`) and new (`.cif`) predictions are handled uniformly
4. Add separate CDR-H3 and CDR-L3 RMSD columns to validate L+.2 asymmetric scaling
5. The primary comparison is L+ v2 vs `baseline_antigen_cut` (the cleanest baseline); secondary comparisons against the other two steered baselines show whether L+ v2 adds value on top of existing steering methods
