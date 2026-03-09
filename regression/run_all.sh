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
# STEP 0: Validate Environment & SMART CACHE MANAGEMENT
# ==============================================================================
log_info "STEP 0: Validating environment and checking config consistency..."

if [ ! -f "config/config.json" ]; then
    log_error "config/config.json not found!"
    exit 1
fi
log_success "config/config.json found"

# Compute hash of current config.json
CONFIG_HASH=$(md5sum config/config.json | awk '{print $1}')
HASH_FILE="data/.config_hash"
CONFIG_CHANGED=false

if [ -f "$HASH_FILE" ]; then
    PREV_HASH=$(cat "$HASH_FILE")
    if [ "$CONFIG_HASH" != "$PREV_HASH" ]; then
        CONFIG_CHANGED=true
        log_warning "config.json has CHANGED since last run. Nuking all previous data..."
        rm -f data/original/*.parquet 2>/dev/null || true
        rm -f data/synthetic/*.parquet 2>/dev/null || true
        rm -f data/*_packed.parquet 2>/dev/null || true
        rm -f results/*.parquet 2>/dev/null || true
        rm -f results/*.csv 2>/dev/null || true
        log_success "Environment cleaned (config changed)."
    else
        log_success "config.json unchanged — will resume from where we left off."
    fi
else
    # No hash file yet — this is either a true first run or a pre-existing run
    # that never had hashing. NEVER nuke existing data; just save the hash.
    log_info "No previous config hash found. Saving current config fingerprint."
    log_success "Existing data (if any) will be kept — resuming."
fi

# Save the current config hash for future runs
mkdir -p data
echo "$CONFIG_HASH" > "$HASH_FILE"

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

# Count how many OD files we expect from config:
# var_types x N_vals x p_vals x rho_vals x sigma_2_vals (continuous + binary grids)
EXPECTED_CONTINUOUS=$(python3 -c "
import json; c=json.load(open('../config/config.json'))
if 'continuous' in c['simulation']['var_type']:
    print(len(c['simulation']['N']) * len(c['simulation']['p']) * len(c['parameters']['rho']) * len(c['parameters']['sigma_2']))
else:
    print(0)
")
EXPECTED_BINARY=$(python3 -c "
import json; c=json.load(open('../config/config.json'))
if 'binary' in c['simulation']['var_type']:
    print(len(c['simulation']['N']) * len(c['simulation']['p']) * len(c['parameters']['rho']) * len(c['parameters']['sigma_2']) * len(c['parameters']['p1']))
else:
    print(0)
")
EXPECTED_OD=$((EXPECTED_CONTINUOUS + EXPECTED_BINARY))
EXISTING_OD=$(ls -1 ../data/original/OD_*.parquet 2>/dev/null | wc -l)

if [ "$EXISTING_OD" -ge "$EXPECTED_OD" ]; then
    log_success "All $EXPECTED_OD OD file(s) already exist — skipping Step 1."
    OD_COUNT=$EXISTING_OD
else
    if [ "$EXISTING_OD" -gt 0 ]; then
        log_info "Resuming OD generation: $EXISTING_OD / $EXPECTED_OD already complete."
    fi
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
fi

cd ..
echo ""

# ==============================================================================
# STEP 2: Generate Synthetic Data (R)
# ==============================================================================
NUM_METHODS=$(python3 -c "import json; c=json.load(open('config/config.json')); print(len(c['synthesis']['methods']))")
METHOD_NAMES=$(python3 -c "import json; c=json.load(open('config/config.json')); print(', '.join(m.upper() for m in c['synthesis']['methods']))")
EXPECTED_SD=$((OD_COUNT * NUM_METHODS))
EXISTING_SD=$(ls -1 data/synthetic/SD_*.parquet 2>/dev/null | wc -l)

log_info "STEP 2: Generating Synthetic Data (SD) via ${METHOD_NAMES} ..."
echo ""

if [ "$EXISTING_SD" -ge "$EXPECTED_SD" ]; then
    log_success "All $EXPECTED_SD SD file(s) already exist — skipping Step 2."
    SD_COUNT=$EXISTING_SD
else
    if [ "$EXISTING_SD" -gt 0 ]; then
        log_info "Resuming SD generation: $EXISTING_SD / $EXPECTED_SD already complete."
    fi

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
fi
echo ""

# ==============================================================================
# STEP 3: Parallel OLS Regression Evaluation (Python)
# ==============================================================================
log_info "STEP 3: Running Parallel OLS Regressions ..."
echo ""

if [ -f "results/aggregated_model_metrics.parquet" ]; then
    log_success "Regression results already exist — skipping Step 3."
else
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
fi
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