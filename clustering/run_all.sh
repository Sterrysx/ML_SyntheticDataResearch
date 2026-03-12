#!/bin/bash

# ==============================================================================
# MASTER EXECUTION SCRIPT: Synthetic Data Clustering Research Pipeline
# ==============================================================================
# Purpose: Execute the complete 3-step clustering simulation pipeline
# Usage:   cd clustering && ./run_all.sh
# ==============================================================================

set -e  # Exit on any error

# ==============================================================================
# ZOMBIE PROCESS CLEANUP TRAP
# ==============================================================================
cleanup() {
    pkill -f "sd_clust_worker.R" 2>/dev/null || true
    pkill -f "od_clust_worker.R" 2>/dev/null || true
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
echo "         SYNTHETIC DATA CLUSTERING RESEARCH PIPELINE                         "
echo "=============================================================================="
echo ""

PIPELINE_START_TIME=$(date +%s)
log_info "Starting full clustering pipeline..."
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
        find data/original -mindepth 1 -delete 2>/dev/null || true
        find data/synthetic -mindepth 1 -delete 2>/dev/null || true
        find results -mindepth 1 -delete 2>/dev/null || true
        log_success "Environment cleaned (config changed)."
    else
        log_success "config.json unchanged — will resume from where we left off."
    fi
else
    log_info "No previous config hash found. Saving current config fingerprint."
    log_success "Existing data (if any) will be kept — resuming."
fi

# Save the current config hash for future runs
mkdir -p data/original data/synthetic results
echo "$CONFIG_HASH" > "$HASH_FILE"

# Check R
if ! command -v Rscript &> /dev/null; then
    log_error "R is not installed or not in PATH"
    exit 1
fi
log_success "R installation verified"

# Check Python
if command -v python3 &> /dev/null; then
    PYTHON_CMD="python3"
elif command -v python &> /dev/null; then
    PYTHON_CMD="python"
else
    log_error "Python is not installed or not in PATH"
    exit 1
fi
log_success "Python installation verified ($PYTHON_CMD)"

echo ""

# ==============================================================================
# STEP 1: Generate Original Data (R)
# ==============================================================================
log_info "STEP 1: Generating Original Data (OD) ..."
echo ""

cd 01_generate_OD

# Count expected OD files from config:
# p_vals x k_vals x separation_vals x rho_vals x n_distributions
EXPECTED_OD=$($PYTHON_CMD -c "
import json; c=json.load(open('../config/config.json'))
dists = c['parameters']['distribution']
n_dist = len(dists) if isinstance(dists, list) else 1
print(len(c['parameters']['p']) * len(c['parameters']['k']) *
      len(c['parameters']['separation']) * len(c['parameters']['rho']) * n_dist)
")
EXISTING_OD=$(ls -1 ../data/original/OD_*.parquet 2>/dev/null | wc -l)

if [ "$EXISTING_OD" -ge "$EXPECTED_OD" ]; then
    log_success "All $EXPECTED_OD OD file(s) already exist — skipping Step 1."
    OD_COUNT=$EXISTING_OD
else
    if [ "$EXISTING_OD" -gt 0 ]; then
        log_info "Resuming OD generation: $EXISTING_OD / $EXPECTED_OD already complete."
    fi
    START_TIME=$(date +%s)

    Rscript 01_generate_original_data.R

    if [ $? -ne 0 ]; then
        log_error "OD generation failed!"
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
M_SYN=$($PYTHON_CMD -c "import json; c=json.load(open('config/config.json')); print(c['simulation']['m'])")
EXPECTED_SD=$((OD_COUNT * M_SYN))
EXISTING_SD=$(ls -1 data/synthetic/SD_*.parquet 2>/dev/null | wc -l)

log_info "STEP 2: Generating Synthetic Data (SD) via CART — $M_SYN reps per OD file ..."
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
# STEP 3: Clustering Evaluation (Python)
# ==============================================================================
log_info "STEP 3: Running Clustering Evaluation ..."
echo ""

if [ -f "results/clustering_results.parquet" ]; then
    log_success "Clustering results already exist — skipping Step 3."
else
    cd 03_clustering_analysis
    START_TIME=$(date +%s)
    $PYTHON_CMD 03_run_clustering.py

    if [ $? -ne 0 ]; then
        log_error "Clustering evaluation failed!"
        exit 1
    fi

    END_TIME=$(date +%s)
    DURATION=$((END_TIME - START_TIME))
    log_success "Clustering evaluation complete in ${DURATION}s"
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
echo "   ├─ Results:         results/clustering_results.parquet"
echo "   └─ Notebooks:       04_evaluation/"
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
