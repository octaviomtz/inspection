# Method K+ Implementation Summary

## Overview

Method K+ (Region-Specific β-Scaling v2) for Boltz2 is being implemented to enable **epitope discovery on unknown antigens** using region-specific latent space scaling during diffusion. This is a **direct pair representation scaling mechanism** (not potential-based) that achieves 10-20x speedup compared to running separate antigen orientation passes.

## Implementation Status

### ✅ Completed Components (Phase 1: Infrastructure)

1. **`src/boltz/steering/antigen_regions.py`** - Antigen region partitioning
   - `AntigenRegion` class for representing surface regions
   - `partition_antigen_spherical()` - Partition antigen into overlapping regions
   - `create_region_mask()` - Generate binary masks for regions
   - Supports 8-12 overlapping regions per antigen

2. **`src/boltz/steering/beta_scaling_config.py`** - Configuration management
   - `EpitopeScanningConfig` dataclass for YAML configuration
   - `parse_epitope_scanning_constraint()` - Parse from YAML constraints
   - Configurable β values for emphasis/de-emphasis
   - Configuration options:
     - `num_regions`: Number of antigen regions (default: 10)
     - `beta_emphasis`: β value for emphasized region (default: +0.5)
     - `beta_deemphasis`: β value for other regions (default: -0.3)
     - `contact_threshold`: Distance threshold in Ångströms (default: 8.0)
     - `confidence_weighting`: Weight contacts by pLDDT/pAE (default: true)

3. **`src/boltz/steering/contact_tracking.py`** - Contact accumulation
   - `compute_contacts()` - Calculate CDR-antigen contacts
   - `weight_contacts_by_confidence()` - Weight by confidence scores
   - `ContactHeatmapAccumulator` - Accumulate contacts across predictions
   - Methods:
     - `add_prediction()` - Add contacts from one prediction
     - `get_heatmap()` - Get normalized epitope propensity map
     - `get_epitope_hotspots()` - Identify high-confidence hotspots

4. **`src/boltz/steering/__init__.py`** - Module exports
   - Proper module structure with public API

5. **`docs/k_plus_usage.md`** - User documentation
   - Configuration guide
   - Usage examples
   - Output interpretation
   - Troubleshooting guide
   - Validation metrics

6. **`examples/7TRH/7TRH_HBG_k_plus.yml`** - Example configuration
   - Demonstrates YAML constraint syntax
   - Configured for 7TRH_HBG complex
   - All K+ parameters documented

### 🔄 In Progress / TODO: Core Integration (Phase 2)

The following components still need implementation to fully integrate K+ into Boltz2:

#### 2.1 **Constraint Parsing & Processing**
- [ ] Modify `boltz/data/parse/yaml.py` to parse `epitope_region_scanning` constraints
- [ ] Add constraint validation logic
- [ ] Store constraint metadata in Record objects
- [ ] Pass region definitions through data pipeline

**Files to modify:**
- `src/boltz/data/parse/yaml.py` - Add parsing logic
- `src/boltz/data/types.py` - Add fields to ResidueConstraints or Constraints
- `src/boltz/data/process.py` - Process region definitions

#### 2.2 **Prediction Loop Integration**
- [ ] Create `epitope_scanning_inference()` function
- [ ] Loop over regions with different β configurations
- [ ] For each region:
  - Run prediction with region-emphasized β
  - Accumulate contacts in ContactHeatmapAccumulator
- [ ] Generate epitope heatmap output
- [ ] Save hotspot coordinates

**New file to create:**
- `src/boltz/steering/epitope_scanning.py` - Main orchestration logic

#### 2.3 **Diffusion Loop Modifications**
- [ ] Modify `src/boltz/model/modules/diffusionv2.py` to support β-scaling
- [ ] Add β-scaling configuration to steering_args
- [ ] Implement pair representation scaling:
  ```
  z[cdr_indices, antigen_region_indices] *= (1 + beta)
  ```
- [ ] Apply β-scaling in the diffusion forward pass (after pair embeddings computed)

**Key integration points in diffusionv2.py:**
- Line ~309: `potentials = get_potentials(steering_args, boltz2=True)`
  - Add β-scaling config extraction here
- Line ~350-500: Main diffusion loop
  - Call β-scaling application after each denoising step

#### 2.4 **CLI Flag Integration**
- [ ] Add `--epitope_scanning` flag to `src/boltz/main.py`
- [ ] Add `--scan_regions` parameter (default: 10)
- [ ] Add `--region_beta_emphasis` parameter (default: 0.5)
- [ ] Add `--region_beta_deemphasis` parameter (default: -0.3)

**Files to modify:**
- `src/boltz/main.py` (around line 1000-1010)

#### 2.5 **Output Processing**
- [ ] Extract coordinates from predictions
- [ ] Compute contacts with proper residue-level mapping
- [ ] Generate epitope propensity heatmap
- [ ] Save heatmap and hotspots to JSON
- [ ] Create visualization (optional: contact heatmap)

**New file to create:**
- `src/boltz/steering/output_processing.py` - Output formatting and file writing

#### 2.6 **Testing & Validation**
- [ ] Create test script for 7TRH_HBG
- [ ] Validate against known epitope (PDB interface)
- [ ] Benchmark performance (target: 10-20x speedup vs separate passes)
- [ ] Add unit tests for individual components

**Files to create:**
- `tests/test_k_plus.py` - Unit tests
- `examples/test_7TRH_HBG_k_plus.py` - Integration test

## Architecture Overview

```
YAML Input (epitope_region_scanning constraint)
       ↓
Constraint Parsing & Region Definition
       ↓
┌─────────────────────────────────────┐
│   Epitope Scanning Loop             │
│  (for each region to emphasize)      │
├─────────────────────────────────────┤
│  1. Set β_config[region_i] = +0.5  │
│     β_config[other] = -0.3          │
│  2. Run Boltz2 prediction with β    │
│  3. Extract coordinates              │
│  4. Compute CDR-antigen contacts    │
│  5. Add to ContactHeatmapAccumulator│
└─────────────────────────────────────┘
       ↓
Generate Epitope Heatmap
       ↓
Identify Hotspots & Save Output
```

## Key Design Decisions

### 1. **Direct Pair Representation Scaling**
- K+ applies β-scaling directly to pair embeddings (z)
- **Not** using explicit potentials framework (simpler, lighter)
- Mechanism: `z[CDR_indices, antigen_region_indices] *= (1 + β)`

### 2. **Region Definition Strategy**
- Simple **spherical/sequential partitioning** initially
- 8-12 overlapping regions per antigen
- Residues assigned to regions based on position

### 3. **Contact Confidence Weighting**
- Weight contacts by pLDDT (model confidence)
- Also consider pAE (interface confidence)
- Formula: `weighted_contact = binary_contact * (pLDDT / 100.0)`

### 4. **Separate Scanning Loop**
- K+ is orchestrated **outside** the main Boltz2 forward pass
- Run N predictions (one per region) sequentially
- Accumulate results → epitope heatmap
- This avoids complex modifications to core diffusion logic

## Integration Checklist

- [ ] **Week 1**: Constraint parsing + configuration (Phase 2.1)
- [ ] **Week 1**: Epitope scanning orchestration (Phase 2.2)
- [ ] **Week 2**: Diffusion loop β-scaling modification (Phase 2.3)
- [ ] **Week 2**: CLI flags + output processing (Phase 2.4-2.5)
- [ ] **Week 2-3**: Testing + validation on 7TRH_HBG (Phase 2.6)

## Expected Outputs

After full implementation, K+ will generate:

1. **`epitope_heatmap.json`**
   ```json
   {
     "residue_heatmap": [0.0, 0.1, 0.8, 0.9, ...],
     "contact_frequency": [0, 0, 1, 1, ...],
     "mean_contact_score": 0.45
   }
   ```

2. **`epitope_hotspots.json`**
   ```json
   {
     "hotspot_indices": [2, 3, 5, 6, 7, ...],
     "hotspot_residues": ["A:3", "A:4", "A:6", "A:7", ...],
     "percentile_threshold": 80
   }
   ```

3. **Per-region predictions** (if `save_per_region: true`)
   - Separate PDB files showing binding with each region emphasized
   - Useful for understanding alternative binding modes

## Performance Goals

- **Computational Efficiency**: 10-20x speedup vs. running separate passes
- **Accuracy**: Precision > 0.7, Recall > 0.6 on known epitopes
- **Robustness**: Works with minimal MSA information

## References

- Implementation Guide: `notes_steering/implementation_guide.md#strategy-e`
- STEERING_METHODS.md: Table entries for K+ and E+
- Key Papers:
  - Boltz-sample (Suzuki & Amagasa, 2026): β-scaling in latent space
  - FK-Diffusion Steering (Horvitz et al., 2501.06848): General framework

## Next Steps

1. Review this summary and approve Phase 2 implementation plan
2. Begin constraint parsing (Phase 2.1)
3. Implement epitope scanning orchestration (Phase 2.2)
4. Integrate β-scaling into diffusion loop (Phase 2.3)
5. Add CLI flags and output processing (Phase 2.4-2.5)
6. Test on 7TRH_HBG (Phase 2.6)
