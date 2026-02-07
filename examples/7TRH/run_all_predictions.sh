#!/bin/bash

################################################################################
# Boltz Prediction Suite for 7TRH Antibody-Antigen Complex
#
# Usage: bash run_all_predictions.sh [--quick-test] [out_dir]
#
# Examples:
#   bash run_all_predictions.sh                    # Full run, default output
#   bash run_all_predictions.sh --quick-test       # Test one case per setting
#   bash run_all_predictions.sh --quick-test /tmp  # Quick test, custom output
#   bash run_all_predictions.sh /tmp               # Full run, custom output
################################################################################

set -e

# Parse arguments
QUICK_TEST=false
OUTPUT_DIR="./boltz_results"

for arg in "$@"; do
    if [ "$arg" = "--quick-test" ]; then
        QUICK_TEST=true
    else
        OUTPUT_DIR="$arg"
    fi
done

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
NUM_SAMPLES=1
NUM_PARTICLES_SINGLE=1
NUM_PARTICLES_MULTI=5

# Colors
GREEN='\033[0;32m'
BLUE='\033[0;34m'
NC='\033[0m'

mkdir -p "$OUTPUT_DIR"

log_step() { echo -e "${BLUE}[$(date '+%H:%M:%S')]${NC} $1"; }
log_success() { echo -e "${GREEN}[✓]${NC} $1"; }

################################################################################
# Case 1: Baseline
################################################################################
log_step "Case 1: Baseline (no constraints)"
boltz predict "$SCRIPT_DIR/7TRH_HBG.yml" \
  --output_format pdb \
  --diffusion_samples $NUM_SAMPLES \
  --out_dir "$OUTPUT_DIR/case1_baseline" \
  --no_kernels \
  > "$OUTPUT_DIR/case1_baseline.log" 2>&1
log_success "case1_baseline/ ($NUM_SAMPLES samples)"

################################################################################
# Case 2: Antigen Steering - Single Particle
################################################################################
log_step "Case 2: Antigen Steering - Single Particle"
boltz predict "$SCRIPT_DIR/7TRH_HBG_antigen_steer.yml" \
  --output_format pdb \
  --use_potentials \
  --antigen_steering \
  --num_particles $NUM_PARTICLES_SINGLE \
  --diffusion_samples $NUM_SAMPLES \
  --no_kernels \
  --out_dir "$OUTPUT_DIR/case2_antigen_steer_p1" \
  > "$OUTPUT_DIR/case2_antigen_steer_p1.log" 2>&1
log_success "case2_antigen_steer_p1/ ($NUM_SAMPLES samples, p=$NUM_PARTICLES_SINGLE)"

################################################################################
# Case 3: Antigen Steering - Multi-Particle
################################################################################
log_step "Case 3: Antigen Steering - Multi-Particle"
boltz predict "$SCRIPT_DIR/7TRH_HBG_antigen_steer.yml" \
  --output_format pdb \
  --use_potentials \
  --antigen_steering \
  --num_particles $NUM_PARTICLES_MULTI \
  --diffusion_samples 1 \
  --no_kernels \
  --out_dir "$OUTPUT_DIR/case3_antigen_steer_p5" \
  > "$OUTPUT_DIR/case3_antigen_steer_p5.log" 2>&1
log_success "case3_antigen_steer_p5/ (1 sample, p=$NUM_PARTICLES_MULTI)"

################################################################################
# Case 4: Contact Constraints
################################################################################
CONTACT_FILES=("$SCRIPT_DIR"/restraint_7TRH_HBG_*.yml)
if [ ${#CONTACT_FILES[@]} -gt 0 ] && [ -f "${CONTACT_FILES[0]}" ]; then
    log_step "Case 4: Contact Constraints (H-bonds, hydrophobic, salt bridges)"

    if [ "$QUICK_TEST" = true ]; then
        # Test only first contact constraint
        CONTACT_FILES=("${CONTACT_FILES[0]}")
    fi

    CASE_NUM=4
    for constraint_file in "${CONTACT_FILES[@]}"; do
        filename=$(basename "$constraint_file" .yml)
        output_subdir="case${CASE_NUM}_${filename}"

        boltz predict "$constraint_file" \
          --output_format pdb \
          --use_potentials \
          --diffusion_samples $NUM_SAMPLES \
          --no_kernels \
          --out_dir "$OUTPUT_DIR/$output_subdir" \
          > "$OUTPUT_DIR/${output_subdir}.log" 2>&1

        log_success "$output_subdir/ ($NUM_SAMPLES samples)"
        CASE_NUM=$((CASE_NUM + 1))
    done
else
    echo "No contact constraint files found (restraint_7TRH_HBG_*.yml)"
fi

################################################################################
# Case 5+: Pocket Constraints
################################################################################
POCKET_FILES=("$SCRIPT_DIR"/restraint_to_A_*.yml)
if [ ${#POCKET_FILES[@]} -gt 0 ] && [ -f "${POCKET_FILES[0]}" ]; then
    log_step "Case 5+: Pocket Constraints"

    if [ "$QUICK_TEST" = true ]; then
        # Test only first pocket constraint
        POCKET_FILES=("${POCKET_FILES[0]}")
    fi

    CASE_NUM=5
    for constraint_file in "${POCKET_FILES[@]}"; do
        filename=$(basename "$constraint_file" .yml)
        output_subdir="case${CASE_NUM}_${filename}"

        boltz predict "$constraint_file" \
          --output_format pdb \
          --use_potentials \
          --diffusion_samples $NUM_SAMPLES \
          --no_kernels \
          --out_dir "$OUTPUT_DIR/$output_subdir" \
          > "$OUTPUT_DIR/${output_subdir}.log" 2>&1

        log_success "$output_subdir/ ($NUM_SAMPLES samples)"
        CASE_NUM=$((CASE_NUM + 1))
    done
else
    echo "No pocket constraint files found (restraint_to_A_*.yml)"
fi

################################################################################
# Summary
################################################################################
echo ""
log_step "Summary"
TOTAL_PDB=$(find "$OUTPUT_DIR" -name "*.pdb" 2>/dev/null | wc -l)
TOTAL_DIRS=$(find "$OUTPUT_DIR" -mindepth 1 -maxdepth 1 -type d 2>/dev/null | wc -l)

echo "  Results: $OUTPUT_DIR"
echo "  Directories: $TOTAL_DIRS"
echo "  PDB files: $TOTAL_PDB"
echo ""

if [ "$QUICK_TEST" = true ]; then
    echo "Quick test mode completed (1 case per setting)"
else
    echo "Full prediction suite completed"
fi

log_success "Done at $(date '+%H:%M:%S')"
