#!/bin/bash

# ==============================================================================
# MASTER EXECUTION SCRIPT: Synthetic Data Regression Research Pipeline
# ==============================================================================
# Purpose: Execute the complete 4-step regression simulation pipeline
# Usage:   cd ofarres_regression && ./run_all.sh
# ==============================================================================

set -e  # Exit on any error

# ==============================================================================
# NEW: ZOMBIE PROCESS CLEANUP TRAP
# Automatically hunts down and kills orphaned SD workers if you hit Ctrl+C
# ==============================================================================
cleanup() {
    pkill -f "sd_worker.R" 2>/dev/null || true
    pkill -f "od_worker.R" 2>/dev/null || true
}
trap cleanup EXIT INT TERM

# Check if conda is available
if ! command -v conda &> /dev/null; then
    echo "❌ ERROR: conda is not installed or not in PATH"
    exit 1
fi

# Activate conda environment
eval "$(conda shell.bash hook)"

if conda env list | grep -q "synthetic_data"; then
    conda activate synthetic_data
    echo "✓ Activated conda environment: synthetic_data"
else
    echo "❌ ERROR: conda environment 'synthetic_data' not found"
    echo "Please create it first: conda create -n synthetic_data python=3.13"
    exit 1
fi

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

log_info()    { echo -e "${BLUE}[INFO]${NC} $1"; }
log_success() { echo -e "${GREEN}[SUCCESS]${NC} $1"; }
log_warning() { echo -e "${YELLOW}[WARNING]${NC} $1"; }
log_error()   { echo -e "${RED}[ERROR]${NC} $1"; }

# Banner
echo "=============================================================================="
echo "        SYNTHETIC DATA REGRESSION RESEARCH PIPELINE                          "
echo "=============================================================================="
echo ""

PIPELINE_START_TIME=$(date +%s)
log_info "Starting full regression pipeline..."
echo ""

# ==============================================================================
# STEP 0: Validate Environment & NUKE OLD CACHES
# ==============================================================================
log_info "STEP 0: Validating environment and cleaning previous runs..."

if [ ! -f "config/config.json" ]; then
    log_error "config/config.json not found!"
    exit 1
fi
log_success "config/config.json found"

# Aggressive cache nuking
log_warning "Nuking all previously generated data and results..."
rm -f data/original/*.parquet 2>/dev/null || true
rm -f data/synthetic/*.parquet 2>/dev/null || true
rm -f data/*_packed.parquet 2>/dev/null || true
rm -f results/*.parquet 2>/dev/null || true
rm -f results/*.csv 2>/dev/null || true
log_success "Environment cleaned successfully."

# Check Python
if command -v python3 &> /dev/null; then
    PYTHON_CMD="python3"
elif command -v python &> /dev/null; then
    PYTHON_CMD="python"
else
    log_error "Python is not installed or not in PATH"
    exit 1
fi

echo ""

# ==============================================================================
# STEP 1: Generate Original Data (R)
# ==============================================================================
log_info "STEP 1: Generating Original Data (OD) ..."
echo ""

cd 01_generate_OD
START_TIME=$(date +%s)

# --- 1A: Continuous Variables ---
log_info "-> Phase 1A: Generating Continuous Data..."
Rscript 01_generate_original_data.R continuous

if [ $? -ne 0 ]; then
    log_error "Phase 1A (Continuous) OD generation failed!"
    exit 1
fi

# --- 1B: Binary Variables ---
echo ""
log_info "-> Phase 1B: Generating Binary Data (This will take longer)..."
Rscript 01_generate_original_data.R binary

if [ $? -ne 0 ]; then
    log_error "Phase 1B (Binary) OD generation failed!"
    exit 1
fi

END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))
OD_COUNT=$(ls -1 ../data/original/OD_*.parquet 2>/dev/null | wc -l)
log_success "OD generation complete in ${DURATION}s (${OD_COUNT} file(s))"
cd ..
echo ""

# ==============================================================================
# STEP 2: Generate Synthetic Data (R)
# ==============================================================================
log_info "STEP 2: Generating Synthetic Data (SD) via CART, NORM, PMM ..."
echo ""

cd 02_generate_SD
START_TIME=$(date +%s)
Rscript 02_generate_synthetic_data.R

if [ $? -ne 0 ]; then
    log_error "SD generation failed!"
    exit 1
fi

END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))
SD_COUNT=$(ls -1 ../data/synthetic/SD_*.parquet 2>/dev/null | wc -l)
log_success "SD generation complete in ${DURATION}s (${SD_COUNT} file(s))"
cd ..
echo ""

# ==============================================================================
# STEP 3: Parallel OLS Regression Evaluation (Python)
# ==============================================================================
log_info "STEP 3: Running Parallel OLS Regressions ..."
echo ""

cd 03_regression_analysis
START_TIME=$(date +%s)
$PYTHON_CMD 03_regression.py

if [ $? -ne 0 ]; then
    log_error "Regression evaluation failed!"
    exit 1
fi

END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))
log_success "Regression evaluation complete in ${DURATION}s"
cd ..
echo ""

# ==============================================================================
# FINAL SUMMARY
# ==============================================================================
echo "=============================================================================="
echo "                    PIPELINE EXECUTION COMPLETE                              "
echo "=============================================================================="
echo ""
log_success "All steps completed successfully!"
echo ""
echo "📁 Output Locations:"
echo "   ├─ Config:          config/config.json"
echo "   ├─ Original Data:   data/original/ (${OD_COUNT} files)"
echo "   ├─ Synthetic Data:  data/synthetic/ (${SD_COUNT} files)"
echo "   └─ Results Data:    results/aggregated_model_metrics.parquet"
echo ""

PIPELINE_END_TIME=$(date +%s)
TOTAL_RUNTIME=$((PIPELINE_END_TIME - PIPELINE_START_TIME))
MINUTES=$(((TOTAL_RUNTIME % 3600) / 60))
SECONDS=$((TOTAL_RUNTIME % 60))

if [ $MINUTES -gt 0 ]; then
    log_success "Total pipeline runtime: ${MINUTES}m ${SECONDS}s (${TOTAL_RUNTIME} seconds)"
else
    log_success "Total pipeline runtime: ${SECONDS}s"
fi
echo "=============================================================================="