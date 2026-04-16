# Evaluation Strategy: Y Early-Phase + B2 Contact Restraints (Experiment 1)

**Date**: 2026-04-15
**Feature**: `hierarchical_steering_early_only` — CDR3 beta-scaling (embedding space, time-decaying) combined with B2 contact restraints (coordinate space)
**Branch**: `y_v2_on_constrains`

---

## 1. Goals of the Feature and What to Measure

### Why this combination was designed
B2 contact restraints fix the global orientation problem: they reduce the fraction of Incorrect predictions from 77% to 37% by physically pulling the known contact residue pairs together. However, B2 does not specifically target CDR loop accuracy — it improves CDR-H3 RMSD only indirectly (by getting the orientation right, CDR3 naturally lands closer to its target).

Y's early-phase beta-scaling improves CDR-H3 RMSD directly (−0.30 Å in Round 1, p=0.025) by amplifying attention between CDR3 residue pairs in the embedding space during the early diffusion steps. This mechanism is fully orthogonal to coordinate-space restraints.

The hypothesis: **B2 gets the antibody pointing at the right place; Y's beta-scaling then ensures CDR3 loops are in a better conformation within that correctly-oriented structure.**

### Primary hypotheses to test
1. **CDR-H3 RMSD**: `Y_early_B2 < B2` (target: below B2's ~2.35 Å). This is the main claim.
2. **DockQ**: `Y_early_B2 ≥ B2` (maintain or improve; do not degrade B2's ~0.291 confidence-selected DockQ).
3. **No regression on safe metrics**: CAPRI tiers should not worsen vs B2; confidence calibration (iptm vs DockQ correlation) should stay comparable.

---

## 2. Baselines

Three baselines, consistent with prior round evaluations:

| Method | Label | Predictions folder | Setting |
|--------|-------|-------------------|---------|
| B1 — Vanilla Boltz-2 | `baseline` | `predictions_examples/antigen_cut/` | Single run per complex |
| B2 — Contact Restraints | `contact_restraints` | `predictions_examples/antigen_cut_contact_restraints/` | Multi-config (several restraint types per complex); pick best |
| B3 — Pocket Guided | `pocket_guided` | `predictions_examples/antigen_cut_vhvl_msa_pocket_ab_boltz_post_2023_omm/` | Multi-config; pick best |
| **New feature** | `y_early_b2` | `predictions_examples/<y_early_b2_subfolder>/` | Multi-config (same restraint configs as B2); pick best |

The new feature is compared primarily against **B2** (the direct parent baseline it extends). B1 is the sanity-check floor and B3 provides context.

An optional 4th comparison point is Y+ alone (`new_feature_hierarchical/`) if predictions exist, to isolate the B2 contribution.

---

## 3. Metrics

### 3.1 Primary: DockQ (via DockQ v2)
DockQ is the community standard for protein–protein docking quality (Basu & Wallner, 2016 + v2 update). It combines iRMSD, LRMSD, and fnat into a single score in [0,1].

- **`top1_ab_ag_dockq`** — DockQ at the antibody-antigen interface, top-1 model by confidence score. This is the primary number. Reported as mean ± 95% bootstrap CI.
- **`oracle_ab_ag_dockq`** — best DockQ among 5 samples (oracle selection). Shows the method's ceiling.
- **`avg_ab_ag_dockq`** — mean DockQ over all 5 samples (shows average trajectory quality).
- **CAPRI classification** — fraction of complexes in each tier (Incorrect <0.23, Acceptable 0.23–0.49, Medium 0.49–0.80, High ≥0.80). Medium+High% is the headline CAPRI stat.

Why DockQ as primary: it captures both interface contact recovery (fnat) and structural accuracy (iRMSD/LRMSD), is standard in the field, and was used consistently across all prior rounds.

### 3.2 Primary: CDR Loop RMSD (after antibody framework alignment)
CDR RMSD after antibody-framework superimposition is the targeted metric for this feature. The alignment strategy:

1. **Superimpose on antibody framework Cα atoms** (non-CDR residues on chains B and C). This removes rigid-body docking error from the RMSD measurement and isolates loop conformation accuracy.
2. In the aligned frame, compute per-CDR Cα RMSD by sequence-aligning prediction to ground truth for each chain.

Reported for all six CDR loops: H1, H2, H3, L1, L2, L3. Key focus:
- **CDR-H3 RMSD** — primary target. B2 baseline is ~2.35 Å; goal is to go below that.
- **CDR-L3 RMSD** — secondary. Y's mechanism acts on all CDR3 loops.

CDR residue positions come from `examples/cdrs.csv` (all six loops per complex).

### 3.3 Secondary: Epitope Prediction Quality
Epitope = antigen residues contacting the antibody in the predicted structure. Compared against the ground truth crystal epitope.

Reported at two thresholds:
- **5 Å (heavy-atom)**: stringent; captures only close contacts
- **8 Å (Cα-Cα)**: coarser; standard for epitope mapping

Metrics: Precision, Recall, F1, MCC per complex, then mean over dataset.

This metric is not the primary target of this feature (B2 already handles orientation well) but serves as a regression check: if Y's beta-scaling somehow disrupts the restraint-driven orientation, epitope F1 will drop.

### 3.4 Secondary: Confidence Metrics
Boltz-2 confidence outputs from the per-model JSON files:
- `confidence_score` = 0.8 × iptm + 0.2 × ptm (used for model selection)
- `iptm` — inter-chain predicted TM-score; main calibration check
- `complex_plddt`, `complex_iplddt`

Check: Spearman correlation between iptm and ab_ag_dockq. If this drops significantly vs B2, the model selection strategy breaks down (as happened with A+).

---

## 4. Model Selection Strategy

With 5 diffusion samples per run:

1. **Top-1 by confidence** (`argmax(confidence_score)`) — primary reported result. Mirrors real-world usage.
2. **Oracle** (`argmax(ab_ag_dockq)`) — upper bound; shows what the method can achieve if selection were perfect.
3. **Mean over 5** — average trajectory quality, independent of selection.

For multi-config methods (B2, B3, Y_early_B2): each complex has multiple run folders (one per restraint configuration). After applying model selection within each run folder, pick the run folder with the highest top-1 DockQ. This mirrors B2's evaluation protocol exactly, ensuring a fair comparison.

**Rationale for using confidence_score (not pLDDT or PAE) for selection**: confidence_score is the standard Boltz-2 output and was used consistently in all prior evaluations. pLDDT and PAE are per-residue/per-pair outputs that don't directly rank full complex quality. Using the same selection criterion across all methods ensures comparability.

---

## 5. Alignment and Post-processing

### For DockQ
DockQ v2 performs its own internal alignment. No pre-alignment needed. The chain map `{A: A, B: B, C: C}` is consistent across all methods and ground truths.

### For CDR RMSD
The evaluation script uses BioPython's `Superimposer` to align on antibody framework Cα atoms (non-CDR residues of chains B and C), then measures per-CDR Cα RMSD in the resulting aligned frame. This is the standard approach for evaluating CDR accuracy in comparative docking experiments.

Sequence alignment (BioPython `PairwiseAligner`) is applied first to map prediction residue indices to ground truth indices, handling any length differences from truncation or numbering offsets.

### Structure format
Both `.pdb` and `.cif` output files are supported. Boltz-2 currently outputs `.cif` files by default; the evaluation script detects extension automatically. Ground truth PDBs in `pdb_minimized/` are `.pdb` format.

---

## 6. Statistical Tests

Paired Wilcoxon signed-rank test (non-parametric, appropriate for non-normal DockQ distributions) on per-complex values.

- Primary: Y_early_B2 vs B2 on `top1_ab_ag_dockq` and `cdr_cdr3_h`
- Secondary: vs B1 and B3 on same metrics
- Also: vs each baseline on `epitope_f1_8A` and `conf_iptm`

Significance thresholds: `*` p<0.05, `**` p<0.01, `***` p<0.001, `ns` otherwise.

Bootstrap 95% confidence intervals (10,000 resamples) on all means.

Report wins/ties/losses alongside p-values.

---

## 7. Prediction Folder Naming Convention

When running boltz on the combined YAML files, the output follows boltz's standard naming:
```
boltz predict yaml_test_set_y_early_b2_contact_restraints_msa_vhvl_antigen_cut/restraint_7TRH_HBG_hbond_19.yml \
  --use_potentials --hierarchical_steering_early_only \
  --diffusion_samples 5 --output_format pdb --outdir <predictions_parent>/<subfolder_name>
```

Output structure per run:
```
<subfolder>/
  boltz_results_restraint_7TRH_HBG_hbond_19/
    predictions/
      restraint_7TRH_HBG_hbond_19/
        restraint_7TRH_HBG_hbond_19_model_0.cif   (or .pdb)
        restraint_7TRH_HBG_hbond_19_model_1.cif
        ...  (5 models total)
        confidence_restraint_7TRH_HBG_hbond_19_model_0.json
        ...
```

The evaluation script discovers run folders by looking for folders containing the complex name inside the subfolder. Multiple restraint configurations (hbond_19, hbond_20, hydrophobic_5, etc.) all contain "7TRH_HBG", so all are found and the best-DockQ setting is selected.

---

## 8. Conda Environment

Use `conda activate dockq2` for the full evaluation.

This environment contains:
- `DockQ` 2.x (Python API: `from DockQ.DockQ import load_PDB, run_on_all_native_interfaces`)
- `biopython` (structure parsing, superimposition, sequence alignment)
- `numpy`, `pandas`, `scipy` (statistics)
- `matplotlib` (plots)

**No new environment is needed.** All required packages are already in `dockq2`.

---

## 9. Can We Reuse `evaluate_hierarchical_v2.py`?

**Yes, with minor additions.** The script was designed for exactly this evaluation pattern:
- It handles multi-setting methods (picks best-DockQ setting per complex) — same logic needed for Y_early_B2
- It handles `.cif` files (auto-detected extension)
- It computes all four metric groups: DockQ, CDR RMSD, epitope, confidence
- It runs paired Wilcoxon tests against all baselines automatically

The only changes needed are adding a 5th method entry for `y_early_b2`:

### Required additions to `evaluate_hierarchical_v2.py`

**1. Add method config** (in `METHOD_CONFIGS`):
```python
"y_early_b2": {
    "subfolder": "y_early_b2_contact_restraints",   # overridable via CLI
    "prefix": "boltz_results_restraint_",
    "suffix": "",
    "multi_setting": True,                           # same as B2
},
```

**2. Add display labels and colors**:
```python
METHOD_COLORS["y_early_b2"] = "#1A6B8A"   # teal/blue — distinct from others
METHOD_DISPLAY["y_early_b2"] = "Y-early\n+ B2"
```

**3. Add CLI argument**:
```python
parser.add_argument("--y_early_b2_subfolder", default="y_early_b2_contact_restraints",
                    help="Subfolder name for Y early-phase + B2 predictions")
```
And wire it: `METHOD_CONFIGS["y_early_b2"]["subfolder"] = args.y_early_b2_subfolder`

**4. Update scatter plot** to also show Y_early_B2 vs B2 (not just vs baseline), since B2 is the primary comparison.

No other changes needed — the DockQ, CDR RMSD, epitope, statistics, and plotting code all work unchanged for the new method.

---

## 10. Run Command (on the evaluation machine)

```bash
conda activate dockq2

python scripts/eval/evaluate_y_early_b2.py \
    --ground_truth_dir pdb_minimized/ \
    --predictions_dir predictions_examples/ \
    --cdrs_csv examples/cdrs.csv \
    --y_early_b2_subfolder y_early_b2_contact_restraints \
    --out_dir evaluation_results_y_early_b2
```

Or using the adapted script with all 4 methods:
```bash
python scripts/eval/evaluate_y_early_b2.py \
    --ground_truth_dir pdb_minimized/ \
    --predictions_dir predictions_examples/ \
    --cdrs_csv examples/cdrs.csv \
    --y_early_b2_subfolder y_early_b2_contact_restraints \
    --out_dir evaluation_results_y_early_b2 \
    --methods baseline contact_restraints pocket_guided y_early_b2
```

To run faster (skip CDR RMSD and epitope for a quick DockQ-only check):
```bash
    --skip_cdr_rmsd --skip_epitope
```

---

## 11. Outputs

| File | Contents |
|------|----------|
| `per_complex_all_metrics.csv` | All metrics for every method × complex combination |
| `summary_table.csv` | Mean, median, 95% CI for each method × metric |
| `statistical_tests.csv` | Wilcoxon results for all method-pairs × metrics |
| `figures/dockq_barplot.png` | Bar chart: Top-1 DockQ, Ab-Ag DockQ, Oracle DockQ |
| `figures/dockq_boxplot.png` | Boxplot: Ab-Ag DockQ distribution across complexes |
| `figures/scatter_y_early_b2_vs_b2.png` | Per-complex scatter: Y_early_B2 vs B2 (primary comparison) |
| `figures/cdr_rmsd_barplot.png` | CDR-H1 through L3 RMSD for all methods |
| `figures/epitope_f1_barplot.png` | Epitope F1 at 5Å and 8Å |
| `figures/confidence_calibration.png` | iptm vs DockQ scatter (calibration check) |
| `figures/capri_stacked_bar.png` | CAPRI tier distribution for all methods |

---

## 12. Interpretation Guide

| Outcome | Interpretation |
|---------|----------------|
| CDR-H3 RMSD Y_early_B2 < B2 (p<0.05) | **Hypothesis confirmed**: beta-scaling adds CDR3 refinement on top of B2's orientation |
| CDR-H3 RMSD Y_early_B2 < B2 (ns) | Trend without significance: probably real but dataset too small (n≈46) |
| CDR-H3 RMSD Y_early_B2 > B2 | Beta-scaling interferes with B2's restraint-guided loops; consider reducing beta_max |
| DockQ Y_early_B2 > B2 | Bonus: combination improves overall docking, not just loops |
| DockQ Y_early_B2 ≈ B2 | Expected; the main gain should be in CDR-H3, not overall DockQ |
| DockQ Y_early_B2 < B2 (large) | Conflict detected: beta-scaling disrupts B2's coordinate-space guidance |
| iptm-DockQ correlation drops | Model selection degrades; the combination confuses Boltz confidence |
| Epitope F1 unchanged or improved | Good; beta-scaling doesn't disrupt the correctly-oriented interface |

**Confidence–oracle gap analysis**: If the gap between oracle DockQ and top-1 DockQ widens vs B2, it means the method generates better structures but the confidence score fails to rank them. This would motivate a composite re-ranker (A+.2 style: `score = α·iptm + (1−α)·restraint_energy`).

---

## 13. Dataset Note

- **48 ground truth PDBs** in `pdb_minimized/` (`.pdb` format, one per complex)
- **48 entries in `cdrs.csv`** with CDR definitions
- **523 YAML files** in `yaml_test_set_y_early_b2_contact_restraints_msa_vhvl_antigen_cut/` (multiple restraint configs per complex, ~11 per complex on average)
- The evaluation script covers only the 46 complexes that are in both the YAML folder and `cdrs.csv`; the 2 extras in cdrs.csv have no YAML counterparts

---

## 14. Summary: What This Evaluation Tests

The single key question: **does CDR3 beta-scaling, operating entirely in embedding space during the early diffusion phase, improve CDR3 loop geometry when applied on top of contact restraints that have already solved the global orientation problem?**

The design ensures clean interpretability:
- Any CDR-H3 RMSD improvement over B2 is attributable entirely to the beta-scaling (since the contact restraints are identical to B2's)
- Any DockQ degradation vs B2 would reveal a conflict between beta-scaling and coordinate-space guidance
- The multi-config best-of selection matches B2's evaluation protocol exactly, ensuring fair comparison
