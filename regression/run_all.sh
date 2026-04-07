#!/bin/bash

# ==============================================================================
# MASTER EXECUTION SCRIPT: Synthetic Data Regression Research Pipeline
# ==============================================================================
# Purpose: Execute the complete regression simulation pipeline (linear + logistic)
# Usage:   cd regression && ./run_all.sh
# ==============================================================================

set -e  # Exit on any error

# ==============================================================================
# ZOMBIE PROCESS CLEANUP TRAP
# Automatically hunts down and kills orphaned workers if you hit Ctrl+C
# ==============================================================================
cleanup() {
    pkill -f "sd_worker.R" 2>/dev/null || true
    pkill -f "od_linear_worker.R" 2>/dev/null || true
    pkill -f "od_logistic_worker.R" 2>/dev/null || true
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
echo "        SYNTHETIC DATA REGRESSION RESEARCH PIPELINE (LINEAR + LOGISTIC)      "
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

# Pre-compute expected counts from config
EXPECTED_LINEAR_CONTINUOUS=$(python3 -c "
import json; c=json.load(open('config/config.json'))
if 'continuous' in c['simulation']['var_type']:
    print(len(c['simulation']['N']) * len(c['simulation']['p']) * len(c['parameters']['rho']) * len(c['parameters']['sigma_2']))
else:
    print(0)
")
EXPECTED_LINEAR_BINARY=$(python3 -c "
import json; c=json.load(open('config/config.json'))
if 'binary' in c['simulation']['var_type']:
    print(len(c['simulation']['N']) * len(c['simulation']['p']) * len(c['parameters']['rho']) * len(c['parameters']['sigma_2']) * len(c['parameters']['p1']))
else:
    print(0)
")
EXPECTED_OD_LINEAR=$((EXPECTED_LINEAR_CONTINUOUS + EXPECTED_LINEAR_BINARY))
EXPECTED_LOGISTIC_CONTINUOUS=$(python3 -c "
import json; c=json.load(open('config/config.json'))
if 'continuous' in c['simulation']['var_type']:
    print(len(c['simulation']['N']) * len(c['simulation']['p']) * len(c['parameters']['rho']) * len(c['parameters']['sigma_2']))
else:
    print(0)
")
EXPECTED_LOGISTIC_BINARY=$(python3 -c "
import json; c=json.load(open('config/config.json'))
if 'binary' in c['simulation']['var_type']:
    print(len(c['simulation']['N']) * len(c['simulation']['p']) * len(c['parameters']['rho']) * len(c['parameters']['sigma_2']) * len(c['parameters']['p1']))
else:
    print(0)
")
EXPECTED_OD_LOGISTIC=$((EXPECTED_LOGISTIC_CONTINUOUS + EXPECTED_LOGISTIC_BINARY))
CONT_METHODS=$(python3 -c "import json; c=json.load(open('config/config.json')); print(', '.join(m.upper() for m in c['synthesis']['continuous_methods']))")
BIN_METHODS=$(python3 -c "import json; c=json.load(open('config/config.json')); print(', '.join(m.upper() for m in c['synthesis']['binary_methods']))")

echo ""

# ==============================================================================
# STEP 1a: Generate Linear Original Data (R)
# ==============================================================================
log_info "STEP 1a: Generating Linear Original Data (OD_linear) ..."
echo ""

cd 01_generate_OD

EXISTING_OD_LINEAR=$(ls -1 ../data/original/OD_linear_*.parquet 2>/dev/null | wc -l)

if [ "$EXISTING_OD_LINEAR" -ge "$EXPECTED_OD_LINEAR" ]; then
    log_success "All $EXPECTED_OD_LINEAR linear OD file(s) already exist — skipping Step 1a."
    OD_LINEAR_COUNT=$EXISTING_OD_LINEAR
else
    if [ "$EXISTING_OD_LINEAR" -gt 0 ]; then
        log_info "Resuming linear OD generation: $EXISTING_OD_LINEAR / $EXPECTED_OD_LINEAR already complete."
    fi
    START_TIME=$(date +%s)

    log_info "-> Phase 1a-i: Generating Continuous Data (linear)..."
    Rscript 01_generate_original_data_linear.R continuous

    if [ $? -ne 0 ]; then
        log_error "Phase 1a-i (Continuous linear) OD generation failed!"
        exit 1
    fi

    echo ""
    log_info "-> Phase 1a-ii: Generating Binary Data (linear)..."
    Rscript 01_generate_original_data_linear.R binary

    if [ $? -ne 0 ]; then
        log_error "Phase 1a-ii (Binary linear) OD generation failed!"
        exit 1
    fi

    END_TIME=$(date +%s)
    DURATION=$((END_TIME - START_TIME))
    OD_LINEAR_COUNT=$(ls -1 ../data/original/OD_linear_*.parquet 2>/dev/null | wc -l)
    log_success "Linear OD generation complete in ${DURATION}s (${OD_LINEAR_COUNT} file(s))"
fi

echo ""

# ==============================================================================
# STEP 1b: Generate Logistic Original Data (R)
# ==============================================================================
log_info "STEP 1b: Generating Logistic Original Data (OD_logistic) ..."
echo ""

EXISTING_OD_LOGISTIC=$(ls -1 ../data/original/OD_logistic_*.parquet 2>/dev/null | wc -l)

if [ "$EXISTING_OD_LOGISTIC" -ge "$EXPECTED_OD_LOGISTIC" ]; then
    log_success "All $EXPECTED_OD_LOGISTIC logistic OD file(s) already exist — skipping Step 1b."
    OD_LOGISTIC_COUNT=$EXISTING_OD_LOGISTIC
else
    if [ "$EXISTING_OD_LOGISTIC" -gt 0 ]; then
        log_info "Resuming logistic OD generation: $EXISTING_OD_LOGISTIC / $EXPECTED_OD_LOGISTIC already complete."
    fi
    START_TIME=$(date +%s)

    log_info "-> Phase 1b-i: Generating Continuous Data (logistic)..."
    Rscript 01_generate_original_data_logistic.R continuous

    if [ $? -ne 0 ]; then
        log_error "Phase 1b-i (Continuous logistic) OD generation failed!"
        exit 1
    fi

    echo ""
    log_info "-> Phase 1b-ii: Generating Binary Data (logistic)..."
    Rscript 01_generate_original_data_logistic.R binary

    if [ $? -ne 0 ]; then
        log_error "Phase 1b-ii (Binary logistic) OD generation failed!"
        exit 1
    fi

    END_TIME=$(date +%s)
    DURATION=$((END_TIME - START_TIME))
    OD_LOGISTIC_COUNT=$(ls -1 ../data/original/OD_logistic_*.parquet 2>/dev/null | wc -l)
    log_success "Logistic OD generation complete in ${DURATION}s (${OD_LOGISTIC_COUNT} file(s))"
fi

OD_COUNT=$((OD_LINEAR_COUNT + OD_LOGISTIC_COUNT))

cd ..
echo ""

# ==============================================================================
# STEP 2: Generate Synthetic Data for ALL OD files (R)
# ==============================================================================
# Case A (linear+cont): 2 arms; Case C (linear+bin): 4 arms
# Case D (logistic+cont): 4 arms; Case B (logistic+bin): 2 arms
EXPECTED_SD_LINEAR=$((EXPECTED_LINEAR_CONTINUOUS * 2 + EXPECTED_LINEAR_BINARY * 4))
EXPECTED_SD_LOGISTIC=$((EXPECTED_LOGISTIC_CONTINUOUS * 4 + EXPECTED_LOGISTIC_BINARY * 2))
EXPECTED_SD=$((EXPECTED_SD_LINEAR + EXPECTED_SD_LOGISTIC))
EXISTING_SD=$(ls -1 data/synthetic/SD_*.parquet 2>/dev/null | wc -l)

log_info "STEP 2: Generating Synthetic Data (SD) [cont: ${CONT_METHODS}, bin: ${BIN_METHODS}] ..."
log_info "  OD total: ${OD_COUNT} (${OD_LINEAR_COUNT} linear + ${OD_LOGISTIC_COUNT} logistic) -> ${EXPECTED_SD} SD files expected"
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
# STEP 3a: Parallel Linear (OLS) Regression (Python)
# ==============================================================================
log_info "STEP 3a: Running Parallel OLS Regressions (linear branch) ..."
echo ""

if [ -f "results/aggregated_model_metrics.parquet" ]; then
    log_success "Linear regression results already exist — skipping Step 3a."
else
    cd 03_regression_analysis
    START_TIME=$(date +%s)
    $PYTHON_CMD 03_regression_linear.py

    if [ $? -ne 0 ]; then
        log_error "Linear regression evaluation failed!"
        exit 1
    fi

    END_TIME=$(date +%s)
    DURATION=$((END_TIME - START_TIME))
    log_success "Linear regression evaluation complete in ${DURATION}s"
    cd ..
fi
echo ""

# ==============================================================================
# STEP 3b: Parallel Logistic Regression (Python)
# ==============================================================================
log_info "STEP 3b: Running Parallel Logistic Regressions (logistic branch) ..."
echo ""

if [ -f "results/aggregated_logistic_metrics.parquet" ]; then
    log_success "Logistic regression results already exist — skipping Step 3b."
else
    cd 03_regression_analysis
    START_TIME=$(date +%s)
    $PYTHON_CMD 03_regression_logistic.py

    if [ $? -ne 0 ]; then
        log_error "Logistic regression evaluation failed!"
        exit 1
    fi

    END_TIME=$(date +%s)
    DURATION=$((END_TIME - START_TIME))
    log_success "Logistic regression evaluation complete in ${DURATION}s"
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
echo "   ├─ Config:            config/config.json"
echo "   ├─ Original Data:     data/original/ (${OD_COUNT} files: ${OD_LINEAR_COUNT} linear + ${OD_LOGISTIC_COUNT} logistic)"
echo "   ├─ Synthetic Data:    data/synthetic/ (${SD_COUNT} files)"
echo "   ├─ Linear Results:    results/aggregated_model_metrics.parquet"
echo "   └─ Logistic Results:  results/aggregated_logistic_metrics.parquet"
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