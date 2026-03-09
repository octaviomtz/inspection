# Evaluation Report: Strategy D+ (Enhanced Canonical Ensemble)

## 1. Original Goal

Strategy D+ aims to improve antibody-antigen structure prediction by generating diverse CDR (Complementarity-Determining Region) conformations through a canonical ensemble approach. The core idea is to apply a scaling factor `z[CDR] *= (1 + beta)` to pair representations during diffusion conditioning, sweeping beta from -0.5 to +0.5 across 11 steps. This produces 11 candidate structures from a single trunk computation, from which the top-5 are selected by confidence score.

**Expected improvements:**
- Greater diversity in CDR loop conformations (especially CDR-H3, the most variable loop)
- Better antibody-antigen interface quality through sampling a broader conformational space
- Computational efficiency: the trunk (MSA + pairformer, ~70% of cost) runs once; only diffusion conditioning, sampling, and confidence scoring run per beta value

**What should remain unchanged:** antigen fold, antibody framework regions, VH-VL pairing geometry.

## 2. Evaluation Setup

- **Test set:** 48 antibody-antigen complexes from `examples/cdrs.csv` (47 successfully evaluated)
- **Ground truth:** OpenMM-minimized crystal structures from `pdb_minimized/`
- **Methods compared:**
  - **Baseline 1 (Boltz2 vanilla):** Standard Boltz2 prediction, 5 models per complex
  - **Baseline 2 (Contact restraints):** Boltz2 with hydrogen bond, hydrophobic, and salt bridge restraints; multiple restraint variants per complex (best variant reported)
  - **Baseline 3 (Pocket+MSA):** Boltz2 with pocket constraints and VH/VL MSA; multiple variants per complex (best variant reported)
  - **D+ (ours):** Canonical ensemble with beta-sweep, top-5 by confidence

- **Evaluation metrics:** DockQ v2 (CAPRI-standard interface quality), CDR-specific RMSD (framework-aligned Ca), ensemble diversity (pairwise CDR-H3 RMSD), and confidence-quality correlation (Spearman rho).

Two evaluation rounds were performed: `d_plus_2/` (D+ vs Baseline 1 only) and `d_plus_4/` (D+ vs all 3 baselines). Both produced identical D+ vs Baseline 1 results. The report focuses primarily on the complete 4-method comparison (`d_plus_4/`).

## 3. Plot-by-Plot Analysis

### 3.1 Plot 1: Global DockQ by Method (Bar Chart)

**How it was produced:** For each complex and method, the top-ranked prediction (model_0) was evaluated with DockQ v2 against the ground truth crystal structure. For Baselines 2 and 3 (which have multiple restraint variants), the best model_0 across all variants was selected. The bar chart shows mean Global DockQ +/- standard deviation across 47 complexes (47 for Baselines 1 and D+, 46 for Baseline 2, 36 for Baseline 3).

**Results:**

| Method | Mean Global DockQ |
|--------|-------------------|
| Boltz2 vanilla | 0.338 +/- 0.184 |
| D+ (ours) | 0.347 +/- 0.176 |
| Contact restraints (best) | 0.499 +/- 0.205 |
| Pocket+MSA (best) | 0.364 +/- 0.197 |

**Conclusion:** D+ shows a negligible improvement of +0.009 over Boltz2 vanilla. This difference is not statistically significant (Wilcoxon p=0.41). Contact restraints (Baseline 2) strongly outperform all other methods with a mean DockQ of 0.499, a +0.161 improvement over D+. Pocket+MSA (Baseline 3) also slightly outperforms D+ at 0.364. The beta-scaling approach does not meaningfully improve global interface quality.

### 3.2 Plot 2: Per-Complex DockQ (Paired Bar Chart)

**How it was produced:** For each of the 47 complexes, the Global DockQ of Boltz2 vanilla (blue) and D+ (orange) best model_0 predictions are plotted side by side, sorted by complex name.

**Results:** The two methods produce nearly identical DockQ values on most complexes. Notable exceptions where D+ improves over vanilla:
- 7Y0O_HLA: D+ 0.478 vs vanilla 0.242 (large gain)
- 8FAH_HLA: D+ 0.634 vs vanilla 0.197 (large gain)
- 8HGM_CDB: D+ 0.402 vs vanilla 0.155 (moderate gain)
- 8OL9_BAH: D+ 0.345 vs vanilla 0.209 (moderate gain)

Notable cases where D+ performs worse:
- 8CDD_EDB: D+ 0.432 vs vanilla 0.714 (large loss)
- 8J1T_HKF: D+ 0.211 vs vanilla 0.714 (large loss)
- 7Y0O_HLA through several complexes show mixed results

**Conclusion:** D+ provides scattered improvements on individual complexes but equally scattered degradations. The wins and losses nearly cancel out (24 wins, 23 losses), confirming the absence of a systematic improvement. The approach helps in some cases and hurts in others, suggesting the beta-scaling effect is case-dependent and unpredictable.

### 3.3 Plot 3: CDR-H3 RMSD Distribution (Box Plot)

**How it was produced:** For each complex and method, the best model_0 prediction was superimposed onto the ground truth using framework residues (non-CDR) of chains B and C. CDR-H3 Ca RMSD was then computed on the framework-aligned structures. The box plot shows the distribution across 47 complexes.

**Results:**

| Method | CDR-H3 RMSD (A) |
|--------|-----------------|
| Boltz2 vanilla | 2.998 +/- 1.993 |
| D+ (ours) | 2.986 +/- 2.058 |

The medians are nearly identical (~2.4 A for vanilla, ~2.6 A for D+). D+ has slightly more extreme outliers (max ~11 A vs ~8.5 A for vanilla).

**Conclusion:** D+ does not improve CDR-H3 loop prediction accuracy. The distributions are virtually indistinguishable. The beta-scaling perturbation does not systematically steer CDR-H3 toward the native conformation. The higher outlier for D+ (11 A) suggests that in some cases the perturbation pushes the loop further from the native state.

### 3.4 Plot 4a: Confidence vs DockQ (Baseline 1)

**How it was produced:** For all 235 predictions (47 complexes x 5 models) of Boltz2 vanilla, the predicted confidence score is plotted against the actual Global DockQ. Spearman rank correlation is computed.

**Results:** rho = 0.287, p = 7.8e-06. Most predictions cluster in the high-confidence range (0.75-1.0) but with widely varying actual quality (DockQ 0.15-0.80). The model is generally overconfident.

**Conclusion:** Boltz2 vanilla's confidence score has a weak but statistically significant positive correlation with actual quality. The model tends to assign high confidence scores regardless of prediction quality, but within that range, relative rankings have some predictive value.

### 3.5 Plot 4b: Confidence vs DockQ (D+)

**How it was produced:** Same methodology as Plot 4a, applied to the 235 D+ predictions (47 complexes x 5 top-K models selected from the beta sweep).

**Results:** rho = 0.177, p = 6.5e-03. The scatter pattern is similar to Baseline 1 but with a weaker correlation. The iptm sub-metric correlation drops to rho = 0.012 (not significant, p = 0.86).

**Conclusion:** D+ degrades the confidence-quality correlation compared to vanilla Boltz2 (rho 0.177 vs 0.287). This is a concerning finding because the entire D+ pipeline relies on confidence scores to select the top-K models from the beta sweep. If confidence is less calibrated for D+ predictions, the top-K selection may not be picking the best structures. The iptm metric becomes essentially uninformative for D+ (rho near zero), suggesting that the beta perturbation disrupts the inter-chain predicted TM-score.

## 4. Additional Quantitative Findings

### 4.1 CAPRI Quality Classification

| Method | High | Medium | Acceptable | Incorrect |
|--------|------|--------|------------|-----------|
| Baseline 1 | 0% | 21.3% | 42.6% | 36.2% |
| D+ (ours) | 0% | 21.3% | 51.1% | 27.7% |
| Contact restraints (best) | 2.2% | 52.2% | 41.3% | 4.3% |
| Pocket+MSA (best) | 0% | 27.8% | 38.9% | 33.3% |

D+ shifts some complexes from "Incorrect" to "Acceptable" (36.2% -> 27.7% incorrect), but none reach "High" quality. Contact restraints dramatically reduce incorrect predictions to only 4.3%.

### 4.2 Oracle Performance (Best Across All 5 Models)

| Method | Oracle DockQ_Global | Oracle DockQ_AB |
|--------|--------------------|-----------------|
| Baseline 1 | 0.391 +/- 0.193 | 0.256 +/- 0.296 |
| D+ (ours) | 0.403 +/- 0.206 | 0.272 +/- 0.300 |
| Contact restraints (best) | 0.545 +/- 0.186 | 0.483 +/- 0.269 |
| Pocket+MSA (best) | 0.404 +/- 0.200 | 0.272 +/- 0.309 |

Oracle performance shows D+ is marginally better than vanilla (+0.012) and on par with Pocket+MSA. Contact restraints remain far superior.

### 4.3 Interface-Specific DockQ (Antigen-Heavy Chain)

The antigen-heavy chain interface (DockQ_AB) is the most relevant for CDR-H3 steering:

| Method | DockQ_AB |
|--------|----------|
| Baseline 1 | 0.183 +/- 0.268 |
| D+ (ours) | 0.199 +/- 0.265 |
| Contact restraints (best) | 0.412 +/- 0.298 |

D+ shows a marginal +0.016 improvement in the antigen-heavy chain interface, but this remains far below contact restraints (+0.213).

### 4.4 Ensemble Diversity

The pairwise CDR-H3 RMSD across 5 models measures structural diversity within each method's ensemble:

| Complex | D+ Mean Pairwise H3 RMSD (A) | Baseline 1 Mean Pairwise H3 RMSD (A) |
|---------|------|------------|
| 7TRH_HBG | 1.63 | 3.77 |
| 8HES_HLC | 1.69 | 5.02 |
| 8SLB_HLA | 3.51 | 6.31 |
| 8DTK_CBA | 7.94 | 4.71 |
| 8OXW_BCA | 1.31 | 0.62 |

Results are mixed. In many cases, D+ actually produces *less* diverse ensembles than vanilla Boltz2 (e.g., 7TRH_HBG: 1.63 vs 3.77 A). This is likely because the confidence-based top-K selection preferentially selects similar high-confidence conformations, while vanilla Boltz2's stochastic sampling naturally produces diverse structures. In some cases D+ produces more diversity (8DTK_CBA: 7.94 vs 4.71).

### 4.5 Statistical Significance (Wilcoxon Signed-Rank Tests)

| Comparison | Mean Diff | p-value | D+ Wins | D+ Losses |
|------------|-----------|---------|---------|-----------|
| D+ vs Baseline 1 | +0.009 | 0.410 | 24 | 23 |
| D+ vs Baseline 2 | -0.149 | 7.2e-12 | 2 | 44 |
| D+ vs Baseline 3 | -0.017 | 0.969 | 16 | 20 |

D+ is not significantly different from Baseline 1 or Baseline 3. It is significantly worse than Baseline 2 (contact restraints), which wins on 44 out of 46 common complexes.

## 5. Overall Conclusions

### 5.1 Did D+ Achieve Its Goals?

**Diversity of CDR conformations:** Partially failed. D+ does not consistently increase ensemble diversity over vanilla Boltz2. The confidence-based top-K selection tends to collapse diversity by preferring structurally similar high-confidence solutions.

**Antibody-antigen interface quality:** Marginally improved but not significantly. The +0.009 DockQ improvement is within noise (p=0.41). Individual complexes show scattered wins and losses with no systematic pattern.

**Computational efficiency:** Achieved. The trunk-once-sweep-many architecture successfully reduces cost to ~(1 + 0.3N) instead of N full predictions. However, the computational savings are moot given the lack of quality improvement.

### 5.2 Why D+ Underperforms

1. **The perturbation is too gentle.** Scaling pair representations by (1 + beta) with beta in [-0.5, +0.5] may not produce sufficient structural variation in the CDR loops. The diffusion process may absorb or smooth out these perturbations.

2. **Confidence-based selection is poorly calibrated.** D+ has a weaker confidence-quality correlation (rho=0.177 vs 0.287), meaning the top-K selection is less effective at identifying truly good structures. The iptm metric becomes completely uninformative (rho=0.012).

3. **The approach lacks physical grounding.** Contact restraints (Baseline 2) embed explicit physical knowledge about the antibody-antigen interface (hydrogen bonds, hydrophobic contacts, salt bridges), which directly guides the prediction toward chemically sensible interactions. D+'s beta-scaling is a generic perturbation of pair representations with no physical interpretation.

4. **Diversity and quality are conflicting objectives.** The confidence module penalizes unusual conformations, so the top-K selection systematically filters out the most diverse (and potentially novel) structures that the beta sweep generates.

### 5.3 Suggested Modifications

1. **Increase perturbation strength and specificity.** Instead of scaling all CDR pair representations uniformly, apply targeted perturbations to CDR-antigen pair representations specifically. Consider larger beta ranges or non-linear scaling functions.

2. **Fix the selection criterion.** Replace or augment the confidence-based top-K selection with a diversity-aware selection (e.g., maximize DockQ-predicted quality subject to a minimum pairwise RMSD constraint). Alternatively, use a composite score that balances confidence and diversity.

3. **Incorporate contact information.** The strong performance of contact restraints suggests that explicit structural knowledge about the binding interface is crucial. D+ could be combined with lightweight contact predictions or evolutionary covariance signals to bias the diffusion toward interface-forming conformations.

4. **Investigate per-CDR beta values.** Rather than a single global beta, sweep independent beta values for each CDR loop (H1, H2, H3, L1, L2, L3) to explore a higher-dimensional conformational space. CDR-H3 may benefit from stronger perturbations than H1/H2.

5. **Calibrate the confidence module.** Fine-tune or re-calibrate the confidence predictor on beta-perturbed structures so it can better distinguish good from bad conformations in the perturbed regime.

6. **Use structure-based reranking.** After generating the ensemble, rerank using external scoring functions (e.g., Rosetta interface energy, shape complementarity) rather than relying solely on the model's internal confidence.
