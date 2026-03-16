# Evaluation Report: Strategy V -- Embedding-Space CDR3 Steering

## Feature Goal

Strategy V steers CDR3 loop conformations by optimizing the Pairformer's pair representation `z` in embedding space during inference. Inspired by the EmbedOpt approach (Li et al.), it computes a target embedding from the mean of non-CDR3 pair representations ("self-reference" mode) and iteratively pushes CDR3 pair embeddings toward that target using analytical gradients (no autograd). The optimization runs for 10 gradient descent steps with a decreasing learning rate (`lr = strength * 0.1 * (1 - step/N)`, strength=1.0) at every denoising step.

The stated goals were:
- **Primary**: More robust CDR3 conformation exploration for novel antibody sequences
- **Secondary**: Improved antibody-antigen structure prediction quality
- **Claimed advantage**: Stability across wide hyperparameter ranges compared to coordinate-space methods

The feature modifies `z` in `DiffusionConditioning.forward()` (after CDR3 beta scaling, before the atom encoder), affecting the entire complex prediction -- not just the CDR3 loops.

## Experimental Setup

- **Baseline (B1)**: Vanilla Boltz2, no steering
- **Feature (V)**: Boltz2 + embedding-space CDR3 steering (self-reference mode, strength=1.0, 10 optimization steps)
- **Test set**: 47 antibody-antigen complexes (out of 48 ground truth structures; 8EQ6_HLA was excluded due to missing predictions)
- **Models per complex**: 5 for both B1 and V
- **Model selection**: Best model per complex selected by highest `confidence_score` (composite of pTM and ipTM)
- **Evaluation script**: `scripts/eval/evaluate_embedding_steering.py`

---

## Figure 1: DockQ Box Plots

**How produced**: For each complex, the best model (by confidence score) was selected. DockQ was computed by running `DockQ` (v2.1.3) comparing each selected prediction against the crystal structure ground truth. The box plots show the distribution of DockQ scores across the 47 complexes. Dashed horizontal lines mark the CAPRI quality thresholds: 0.23 (acceptable), 0.49 (medium), 0.80 (high).

### Subplot 1 -- Global DockQ
The Global DockQ is the mean of all pairwise interface DockQ scores (A-B, A-C, B-C). B1 (median ~0.31) outperforms V (median ~0.24). The B1 interquartile range (IQR) extends higher, with more complexes reaching medium-quality predictions. V's distribution is shifted downward overall.

**Conclusion**: V degrades the overall docking quality. The mean drops from 0.338 (B1) to 0.304 (V), a statistically significant decrease (Wilcoxon p=0.0004, W/T/L=12/0/35 -- B1 wins on 35 of 47 complexes).

### Subplot 2 -- A-B (Antigen-Heavy Chain)
Both methods perform poorly at this interface (medians near 0.05), indicating that the antigen-heavy chain interface is difficult for both. V shows more outliers at higher DockQ values but also a lower box overall.

**Conclusion**: No significant difference (p=0.098). The feature neither helps nor clearly harms this already-poor interface.

### Subplot 3 -- A-C (Antigen-Light Chain)
Similar to A-B, both methods have most predictions below 0.23 (incorrect). The distributions are nearly indistinguishable.

**Conclusion**: No significant difference (p=0.358, delta=+0.003). The feature is neutral on this interface.

### Subplot 4 -- B-C (Heavy-Light Chain)
This is the internal antibody interface. B1 has a tight distribution centered around 0.65, while V drops significantly to ~0.58 with a wider spread. This is the most dramatic degradation across all interfaces.

**Conclusion**: V significantly degrades the heavy-light chain interface (p<0.0001, delta=-0.082, B1 wins 41 of 47 complexes). This is a critical finding: the embedding-space modification disrupts the internal antibody fold.

---

## Figure 2: DockQ Scatter -- V vs B1

**How produced**: Each point represents one complex. The x-coordinate is the B1 (best-by-confidence) Global DockQ; the y-coordinate is the V (best-by-confidence) Global DockQ. Points below the diagonal (dashed line) indicate B1 performs better.

The scatter shows that most points fall below the diagonal: B1 wins on 35 complexes, V wins on 12. The degradation is most pronounced for medium-difficulty complexes (B1 DockQ 0.5--0.8), where V scores consistently drop by 0.05--0.15. For hard complexes (B1 DockQ <0.3), the two methods produce similar results (points cluster near the diagonal). Two notable outliers where V scores ~0.75 (vs B1 ~0.8) exist at the top, and a few hard cases where V slightly improves.

**Conclusion**: B1 wins nearly 3:1 overall. The feature is most damaging on cases where Boltz2 already produces reasonable predictions, while being roughly neutral on the hardest cases.

---

## Figure 3: CAPRI Quality Distribution

**How produced**: Each complex's Global DockQ (best-by-confidence model) was classified into CAPRI quality categories: incorrect (<0.23), acceptable (0.23--0.49), medium (0.49--0.80), high (>=0.80). The stacked bar chart shows the percentage distribution for each method.

- **B1**: 36% incorrect (17), 43% acceptable (20), 21% medium (10), 0% high
- **V**: 51% incorrect (24), 30% acceptable (14), 19% medium (9), 0% high

**Conclusion**: V shifts 7 complexes from acceptable/medium into the incorrect category. The proportion of incorrect predictions grows from 36% to 51%. Neither method achieves high-quality predictions on this test set.

---

## Figure 4: CDR3 RMSD Box Plots

**How produced**: For each complex, the CDR3 RMSD was computed after Kabsch superposition on the antibody framework (non-CDR CA atoms). CDR indices were taken from `examples/cdrs.csv` (per-complex, 1-indexed positions). Sequence alignment (BioPython PairwiseAligner) was used to map prediction residues to ground truth residues, accounting for different chain lengths and numbering schemes. Best-by-confidence model per complex.

### Subplot 1 -- CDR3-H (Heavy Chain CDR3)
B1 median ~2.5 Angstrom with IQR 1.7--4.1. V median ~4.2 with a tighter but higher IQR 3.3--4.8. The entire V distribution is shifted upward.

**Conclusion**: V worsens CDR3-H RMSD by +1.21 Angstrom on average (p<0.0001, B1 wins 38 of 47). This is the opposite of the feature's primary goal.

### Subplot 2 -- CDR3-L (Light Chain CDR3)
B1 median ~1.0 Angstrom (IQR 0.7--1.6), V median ~3.0 (IQR 2.4--3.5). This is the most severe degradation: V nearly triples the light-chain CDR3 error.

**Conclusion**: V catastrophically worsens CDR3-L RMSD (+1.66 Angstrom on average, p<0.0001, B1 wins 46 of 47). Only 1 complex out of 47 saw improvement.

### Subplot 3 -- CDR3 Combined
B1 median ~2.0 Angstrom, V median ~3.7 Angstrom. The combined metric confirms a consistent degradation.

**Conclusion**: V worsens overall CDR3 accuracy by +1.28 Angstrom on average (p<0.0001, B1 wins 42 of 47). The feature fails its primary goal of improving CDR3 loop conformations.

---

## Figure 5: DockQ by Difficulty (Stratified)

**How produced**: Complexes were stratified into difficulty bins based on B1 Global DockQ: easy (>=0.80), medium (0.49--0.80), hard (<0.49). The bar chart shows mean Global DockQ for each method within each difficulty bin.

- **Easy** (n=0): No complexes reached high-quality in B1, so this bin is empty.
- **Medium** (n=10): B1 mean 0.67, V mean 0.53. V drops these by ~0.14 DockQ points.
- **Hard** (n=37): B1 mean 0.25, V mean 0.24. Nearly identical performance.

**Conclusion**: V is most harmful for medium-difficulty complexes where B1 already produces reasonable predictions. For hard complexes (the vast majority of this test set), the feature has negligible impact. The feature does not help on any difficulty tier.

---

## Confidence Metrics

The summary statistics reveal that V produces systematically lower confidence scores across all metrics:

| Metric | B1 | V | Delta |
|--------|-----|-----|-------|
| Confidence | 0.883 | 0.835 | -0.048 |
| ipTM | 0.737 | 0.584 | -0.153 |
| pTM | 0.824 | 0.737 | -0.087 |
| Complex pLDDT | 0.919 | 0.897 | -0.022 |
| CDR3-H pLDDT | 0.816 | 0.689 | -0.127 |
| CDR3-L pLDDT | 0.876 | 0.760 | -0.116 |

All differences are highly significant (p<0.0001). The ipTM drop (-0.153) is particularly large, indicating that the model's own assessment of interface quality degrades substantially. CDR3-specific pLDDT drops by 0.12--0.13, suggesting the model recognizes the CDR3 regions are less well-folded after steering.

**Conclusion**: The model itself "knows" the predictions are worse. Lower confidence scores across the board indicate the embedding-space perturbation pushes `z` away from the learned manifold rather than toward better conformations.

## Contact Analysis

Contact F1 (0.181 vs 0.173, p=0.755) and epitope recall (0.364 vs 0.405, p=0.925) show no statistically significant differences. The feature does not improve or degrade epitope contact prediction.

---

## Overall Conclusions

**Strategy V fails to achieve its goals.** The embedding-space CDR3 steering feature degrades prediction quality on nearly every metric:

1. **CDR3 accuracy worsens** (primary goal failure): CDR3-H RMSD increases by +1.21 Angstrom, CDR3-L by +1.66 Angstrom. The feature was supposed to improve CDR3 conformations but does the opposite.

2. **Overall docking quality degrades**: Global DockQ drops by -0.034 (p=0.0004), with B1 winning on 35 of 47 complexes. The B-C (antibody internal) interface is particularly harmed (-0.082).

3. **Confidence metrics collapse**: ipTM drops by -0.153, CDR3-specific pLDDT drops by ~0.12. The model itself recognizes the perturbation is harmful.

4. **No benefit on any subset**: The feature does not help hard, medium, or easy cases. It is most damaging on medium-difficulty complexes.

### Root Cause Analysis

The self-reference mode computes `z_target` as the mean of all non-CDR3 pair embeddings and pushes CDR3 embeddings toward this average. This is problematic because:

- **CDR3 loops are inherently different from framework regions.** Pushing CDR3 pair representations toward the mean of non-CDR3 (framework) pairs erases the specificity that CDR3 loops need. CDR3 loops are the most variable and structurally diverse parts of the antibody; forcing them to resemble framework regions removes the signal the model needs to predict their unique conformations.

- **The perturbation is applied at every denoising step.** This means the modification compounds across all 200 denoising steps, amplifying the damage.

- **The B-C (heavy-light) interface degradation** suggests the perturbation corrupts pair representations not only within CDR3 but also in how CDR3 interacts with the rest of the antibody.

### Suggested Modifications

If this approach is to be revisited, the following modifications could be explored:

1. **Alternative target**: Instead of self-reference (mean of framework), use CDR3 embeddings from known structures of the same canonical class. The current target is fundamentally wrong for CDR3 biology.

2. **Reduced strength / fewer steps**: The current strength=1.0 with 10 steps may be too aggressive. Sweeping strength in [0.01, 0.1, 0.5] and num_opt_steps in [1, 3, 5] could reveal a regime where the perturbation is mild enough to explore without destroying.

3. **Apply only at early denoising steps**: Rather than modifying `z` at every step, restrict the steering to the first N% of denoising steps, allowing the model to recover in later steps.

4. **Cross-attention target instead of mean**: Use the CDR3-to-antigen pair embeddings as the target, rather than framework self-pairs. This would steer CDR3 toward better binding rather than toward generic well-ordered structure.

5. **Selective pair mask**: Only modify CDR3-CDR3 pairs (not CDR3-framework or CDR3-antigen pairs), to avoid corrupting the interface information that drives docking.
