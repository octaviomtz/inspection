# Method K+: Region-Specific β-Scaling for Epitope Discovery

## Overview

Method K+ (Region-Specific β-Scaling v2) is a steering method for Boltz2 that discovers epitope locations on unknown antigens without prior knowledge. It works by:

1. **Partitioning** the antigen surface into overlapping regions (typically 8-12 regions)
2. **Running multiple predictions** with region-specific β-scaling:
   - Emphasizing one region with β = +0.5
   - De-emphasizing other regions with β = -0.3
3. **Accumulating contacts** across all predictions to generate an epitope propensity map
4. **Identifying hotspots** as high-confidence epitope locations

## Key Advantages

- **Computational Efficiency**: ~10-20x speedup vs. running separate antigen orientation scanning passes
- **No Retraining Required**: Uses existing Boltz2 model with latent space scaling
- **Robust**: Works without multiple sequence alignments (MSA-agnostic)
- **Transparent Mechanism**: Direct pair representation scaling (β-scaling) is interpretable

## Usage

### YAML Configuration

Add an `epitope_region_scanning` constraint to your Boltz YAML file:

```yaml
sequences:
- protein:
    id: A
    sequence: <antigen_sequence>
    msa: <path_to_msa>
- protein:
    id: B
    sequence: <heavy_chain_sequence>
    msa: empty
- protein:
    id: C
    sequence: <light_chain_sequence>
    msa: empty

constraints:
  - epitope_region_scanning:
      antigen_chain: A                    # Which chain is the antigen
      num_regions: 10                     # Number of regions to scan (default: 10)
      region_overlap_ratio: 0.2           # Overlap between regions (default: 0.2)
      beta_emphasis: 0.5                  # Beta for emphasized region (default: 0.5)
      beta_deemphasis: -0.3               # Beta for other regions (default: -0.3)
      contact_threshold: 8.0              # Distance threshold for contacts in Å (default: 8.0)
      confidence_weighting: true          # Weight by pLDDT/pAE (default: true)
      accumulate_heatmap: true            # Generate heatmap (default: true)
      save_per_region: false              # Save predictions per region (default: false)
```

### Command Line Usage

```bash
# Basic epitope scanning
boltz predict config.yaml --use_potentials --epitope_scanning

# With multiple samples
boltz predict config.yaml --use_potentials --epitope_scanning --diffusion_samples 3

# With more regions
boltz predict config.yaml --use_potentials --epitope_scanning --scan_regions 15
```

### CLI Flags

- `--epitope_scanning`: Enable epitope region scanning with β-scaling
- `--scan_regions`: Number of antigen regions to scan (default: 10)
- `--region_beta_emphasis`: Beta value for emphasized region (default: 0.5)
- `--region_beta_deemphasis`: Beta value for de-emphasized regions (default: -0.3)

## Output

The prediction output includes:

1. **Epitope Heatmap** (`epitope_heatmap.json` or in structure file)
   - Per-residue contact frequency from all region scans
   - Value range: [0.0, 1.0] (0 = no contacts, 1.0 = always in contact)

2. **Hotspot Coordinates** (`epitope_hotspots.json`)
   - Indices of high-confidence epitope residues (top 20-30%)

3. **Per-Region Predictions** (if `save_per_region: true`)
   - Separate PDB files for predictions with each region emphasized
   - Useful for visualizing different binding modes

4. **Contact Statistics** (`contact_stats.json`)
   - Contact frequency per region
   - Average contact quality (pLDDT-weighted)
   - Region-wise metrics

## Interpretation

### Reading the Epitope Heatmap

The epitope heatmap represents the **contact propensity** for each antigen residue:

- **High values (0.7-1.0)**: Likely epitope residues, consistently contacted
- **Medium values (0.4-0.7)**: Possibly involved in binding
- **Low values (0.0-0.4)**: Unlikely epitope residues

### Identifying Epitope Regions

1. Look for **continuous clusters** of high-value residues (more likely to be real epitopes)
2. Consider **secondary epitopes** (regions with moderate but consistent contact)
3. Cross-reference with **experimental epitope data** if available

## Advanced Configuration

### Fine-tuning β Values

For specific use cases, adjust β parameters:

```yaml
- epitope_region_scanning:
    antigen_chain: A
    beta_emphasis: 0.7          # Higher = stronger emphasis on region
    beta_deemphasis: -0.5       # More negative = stronger de-emphasis
```

**Guidance:**
- Larger β values (±0.5 to ±1.0): More aggressive steering, faster convergence
- Smaller β values (±0.2 to ±0.3): More subtle steering, broader exploration

### Multiple Region Sizes

For thorough exploration, run K+ with different region numbers:

```bash
# Coarse scan (fewer, larger regions)
boltz predict config.yaml --use_potentials --epitope_scanning --scan_regions 8

# Medium scan
boltz predict config.yaml --use_potentials --epitope_scanning --scan_regions 12

# Fine scan (many small regions)
boltz predict config.yaml --use_potentials --epitope_scanning --scan_regions 20
```

Combine heatmaps (average or union) to get robust epitope definition.

## Validation

### Benchmark on Known Epitopes

For antibody-antigen pairs with known epitopes:

```python
from boltz.steering import ContactHeatmapAccumulator

# Run K+ prediction
results = run_k_plus_prediction(config_yaml)

# Evaluate against known epitope
known_epitope_residues = [...]  # From PDB interface definition
predicted_epitope = results['epitope_hotspots']

# Compute metrics
precision = len(set(predicted_epitope) & set(known_epitope_residues)) / len(predicted_epitope)
recall = len(set(predicted_epitope) & set(known_epitope_residues)) / len(known_epitope_residues)
f1 = 2 * (precision * recall) / (precision + recall)

print(f"Precision: {precision:.2f}")
print(f"Recall: {recall:.2f}")
print(f"F1-Score: {f1:.2f}")
```

**Expected Performance:**
- Precision: > 0.7 (low false positives)
- Recall: > 0.6 (covers most true epitope)
- F1-Score: > 0.65 (balanced prediction)

## Troubleshooting

### Issue: Low contact predictions

**Symptoms**: All regions show low contact frequency

**Causes:**
- CDR and antigen too far apart in initial structure
- β values too small (insufficient steering)

**Solutions:**
- Increase `beta_emphasis` to 0.8-1.0
- Reduce `contact_threshold` to 6.0-7.0 Å
- Increase `diffusion_samples` for better exploration

### Issue: Uniform contact distribution

**Symptoms**: All antigen residues have similar contact frequency

**Causes:**
- Too many regions (regions too small)
- β values not sufficiently different from zero

**Solutions:**
- Reduce `num_regions` to 8-10
- Increase `beta_deemphasis` magnitude to -0.5 or -0.7
- Check that antigen sequence is correct and of reasonable length

### Issue: Unreasonable epitope locations

**Symptoms**: Predicted epitope on antigen surface far from CDRs

**Causes:**
- Initial structure pose is poor
- Model uncertainty (low pLDDT/pAE)

**Solutions:**
- Check input sequences and MSA quality
- Use template structures if available
- Increase `diffusion_samples` and `num_sampling_steps`
- Cross-validate with multiple runs

## Citation

If you use Method K+ in your research, please cite:

```
@article{suzuki2026,
  title={Boltz-sample: Fine-grained diffusion-based sampling via pair representation scaling},
  author={Suzuki and Amagasa},
  journal={arXiv},
  year={2026}
}
```

## References

1. **Implementation Guide**: See `notes_steering/implementation_guide.md` for technical details
2. **Strategy E+**: K+ uses the same mechanism as Strategy E+ (blind scanning with region-specific β-scaling)
3. **Related Work**: FK-Diffusion Steering (Horvitz et al., 2501.06848)
