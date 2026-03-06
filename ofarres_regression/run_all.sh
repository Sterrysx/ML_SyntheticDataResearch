#!/bin/bash

# ==============================================================================
# MASTER EXECUTION SCRIPT: Synthetic Data Regression Research Pipeline
# ==============================================================================
# Purpose: Execute the complete 4-step regression simulation pipeline
# Usage:   cd ofarres_regression && ./run_all.sh
# ==============================================================================

set -e  # Exit on any error

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
# STEP 0: Validate Environment
# ==============================================================================
log_info "STEP 0: Validating environment..."

if [ ! -f "config/config.json" ]; then
    log_error "config/config.json not found!"
    exit 1
fi
log_success "config/config.json found"

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

# Check Python packages
log_info "Checking Python dependencies..."
$PYTHON_CMD -c "import pandas, numpy, sklearn, matplotlib, seaborn, scipy" 2>/dev/null
if [ $? -ne 0 ]; then
    log_warning "Some Python packages are missing. Installing from requirements.txt..."
    pip install -r requirements.txt
fi
log_success "Python dependencies satisfied"

# Check R packages
log_info "Checking R dependencies..."
Rscript -e "pkgs <- c('jsonlite', 'mvtnorm', 'synthpop'); if (!all(pkgs %in% installed.packages()[,'Package'])) { quit(status=1) }" 2>/dev/null
if [ $? -ne 0 ]; then
    log_warning "Some R packages are missing. Installing..."
    Rscript -e "install.packages(c('jsonlite', 'mvtnorm', 'synthpop'), repos='http://cran.us.r-project.org')"
fi
log_success "R dependencies satisfied"

echo ""

# ==============================================================================
# STEP 1: Generate Original Data (R)
# ==============================================================================
log_info "STEP 1: Generating Original Data (OD) ..."
echo ""

cd 01_generate_OD

# Clear old data (Fixed for Parquet)
if [ -d "../data/original" ]; then
    log_warning "Removing old original data..."
    rm -f ../data/original/*.parquet 2>/dev/null || true
fi

START_TIME=$(date +%s)
Rscript 01_generate_original_data.R

if [ $? -ne 0 ]; then
    log_error "OD generation failed!"
    exit 1
fi

END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))

# Fixed counting for Parquet files
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

# Clear old data (Fixed for Parquet)
if [ -d "../data/synthetic" ]; then
    log_warning "Removing old synthetic data..."
    rm -f ../data/synthetic/*.parquet 2>/dev/null || true
fi

START_TIME=$(date +%s)
Rscript 02_generate_synthetic_data.R

if [ $? -ne 0 ]; then
    log_error "SD generation failed!"
    exit 1
fi

END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))

# Fixed counting for Parquet files
SD_COUNT=$(ls -1 ../data/synthetic/SD_*.parquet 2>/dev/null | wc -l)
log_success "SD generation complete in ${DURATION}s (${SD_COUNT} file(s))"
cd ..
echo ""


END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))

log_success "Evaluation complete in ${DURATION}s"
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
echo "   ├─ Models:          models/ (${MODEL_COUNT} pkl files)"
echo "   ├─ Results CSV:     results/evaluation_metrics.csv"
echo "   └─ Plots:           results/plot_*.png"
echo ""
echo "📊 Next Steps:"
echo "   1. Review metrics:  cat results/evaluation_metrics.csv"
echo "   2. View plots:      open results/plot_*.png"
echo "   3. Explore models:  jupyter notebook 03_regression_analysis/"
echo "   4. Explore eval:    jupyter notebook 04_evaluation/"
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