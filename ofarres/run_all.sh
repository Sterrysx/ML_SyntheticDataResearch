#!/bin/bash

# ==============================================================================
# MASTER EXECUTION SCRIPT: Synthetic Data Clustering Research Pipeline
# ==============================================================================
# Purpose: Execute the complete 3-step simulation pipeline
# Author: ML Research Team
# Date: December 19, 2025
# ==============================================================================

set -e  # Exit on any error

# Check if conda is available
if ! command -v conda &> /dev/null; then
    echo "❌ ERROR: conda is not installed or not in PATH"
    exit 1
fi

# Activate conda environment
eval "$(conda shell.bash hook)"

# Check if synthetic_data environment exists and activate it
if conda env list | grep -q "synthetic_data"; then
    conda activate synthetic_data
    echo "✓ Activated conda environment: synthetic_data"
else
    echo "❌ ERROR: conda environment 'synthetic_data' not found"
    echo "Please create it first or update the environment name in this script"
    exit 1
fi

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Logging functions
log_info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

log_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

log_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Banner
echo "=============================================================================="
echo "           SYNTHETIC DATA CLUSTERING RESEARCH PIPELINE                       "
echo "=============================================================================="
echo ""

# Start timer
PIPELINE_START_TIME=$(date +%s)

log_info "Starting full simulation pipeline..."
echo ""

# ==============================================================================
# STEP 0: Validate Environment
# ==============================================================================
log_info "STEP 0: Validating environment..."

# Check if config.json exists
if [ ! -f "config/config.json" ]; then
    log_error "config/config.json not found! Please ensure it exists in the config directory."
    exit 1
fi
log_success "config/config.json found"

# Check R installation
if ! command -v Rscript &> /dev/null; then
    log_error "R is not installed or not in PATH"
    exit 1
fi
log_success "R installation verified"

# Check Python installation
if ! command -v python &> /dev/null && ! command -v python3 &> /dev/null; then
    log_error "Python is not installed or not in PATH"
    exit 1
fi

# Set Python command
if command -v python3 &> /dev/null; then
    PYTHON_CMD="python3"
else
    PYTHON_CMD="python"
fi
log_success "Python installation verified ($PYTHON_CMD)"

# Check for required Python packages
log_info "Checking Python dependencies..."
$PYTHON_CMD -c "import pandas, numpy, sklearn, joblib, tqdm" 2>/dev/null
if [ $? -ne 0 ]; then
    log_warning "Some Python packages are missing. Installing from requirements.txt..."
    pip install -r requirements.txt
fi
log_success "Python dependencies satisfied"

# Check for required R packages
log_info "Checking R dependencies..."
Rscript -e "pkgs <- c('jsonlite', 'mvtnorm', 'synthpop', 'dplyr'); if (!all(pkgs %in% installed.packages()[,'Package'])) { quit(status=1) }" 2>/dev/null
if [ $? -ne 0 ]; then
    log_warning "Some R packages are missing. Installing..."
    Rscript installer_scripts/install_requirements.R
fi
log_success "R dependencies satisfied"

echo ""

# ==============================================================================
# STEP 1: Data Generation (R)
# ==============================================================================
log_info "STEP 1: Generating Real and Synthetic datasets (R)..."
log_info "Expected output: 360 real + 3600 synthetic CSV files (with updated config)"
log_info "Estimated time: 5-10 minutes"
echo ""

cd 01_data_generation

# Clear old data if exists
if [ -d "../data/real" ]; then
    log_warning "Removing old real data..."
    rm -rf ../data/real/*.csv 2>/dev/null || true
fi

if [ -d "../data/synthetic" ]; then
    log_warning "Removing old synthetic data..."
    rm -rf ../data/synthetic/*.csv 2>/dev/null || true
fi

# Run R script
START_TIME=$(date +%s)
Rscript 01_generate_synthetic.R

if [ $? -ne 0 ]; then
    log_error "Data generation failed! Check R script errors above."
    exit 1
fi

END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))

# Verify output
REAL_COUNT=$(ls -1 ../data/real/*.csv 2>/dev/null | wc -l)
SYN_COUNT=$(ls -1 ../data/synthetic/*.csv 2>/dev/null | wc -l)

log_success "Data generation complete in ${DURATION}s"
log_info "Generated files: ${REAL_COUNT} real, ${SYN_COUNT} synthetic"

if [ "$REAL_COUNT" -ne 720 ] || [ "$SYN_COUNT" -ne 7200 ]; then
    log_warning "Expected 720 real and 7200 synthetic files, but got ${REAL_COUNT} and ${SYN_COUNT}"
fi

cd ..
echo ""

# ==============================================================================
# STEP 2: Clustering Simulation (Python)
# ==============================================================================
log_info "STEP 2: Running clustering simulation (Python)..."
log_info "Processing all 3600 synthetic datasets"
log_info "Estimated time: 20-30 minutes"
echo ""

cd 02_clustering

# Remove old results if exists
if [ -f "clustering_results_final.csv" ]; then
    log_warning "Removing old results file..."
    rm clustering_results_final.csv
fi

# Run Python simulation
START_TIME=$(date +%s)
$PYTHON_CMD 02_run_simulation.py

if [ $? -ne 0 ]; then
    log_error "Clustering simulation failed! Check Python script errors above."
    exit 1
fi

END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))

# Verify output
if [ ! -f "clustering_results_final.csv" ]; then
    log_error "Output file clustering_results_final.csv not found!"
    exit 1
fi

RESULT_ROWS=$(wc -l < clustering_results_final.csv)
log_success "Clustering simulation complete in ${DURATION}s"
log_info "Results: ${RESULT_ROWS} rows (expected 3601 including header)"

cd ..
echo ""

# ==============================================================================
# STEP 3: Analysis & Visualization (Jupyter)
# ==============================================================================
log_info "STEP 3: Generating analysis and visualizations..."
log_info "Output: 3 publication-ready plots"
echo ""

cd 03_analysis

# Convert notebook to script and run
log_info "Executing analysis notebook..."
START_TIME=$(date +%s)

$PYTHON_CMD -m jupyter nbconvert --to notebook --execute 03_plot_results.ipynb --output 03_plot_results_executed.ipynb --ExecutePreprocessor.timeout=300

if [ $? -ne 0 ]; then
    log_warning "Notebook execution via nbconvert failed. Trying alternative method..."
    
    # Alternative: Use papermill if available
    if command -v papermill &> /dev/null; then
        papermill 03_plot_results.ipynb 03_plot_results_executed.ipynb
    else
        log_error "Could not execute notebook automatically."
        log_info "Please run the notebook manually: jupyter notebook 03_plot_results.ipynb"
    fi
fi

END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))

# Check for generated plots
PLOT_COUNT=$(ls -1 plot*.png 2>/dev/null | wc -l)
log_success "Analysis complete in ${DURATION}s"
log_info "Generated plots: ${PLOT_COUNT} (expected 3)"

if [ "$PLOT_COUNT" -eq 3 ]; then
    log_success "All visualizations generated successfully!"
    ls -1 plot*.png
else
    log_warning "Expected 3 plots, found ${PLOT_COUNT}"
fi

cd ..
echo ""

# ==============================================================================
# FINAL SUMMARY
# ==============================================================================
echo "=============================================================================="
echo "                          PIPELINE EXECUTION COMPLETE                         "
echo "=============================================================================="
echo ""
log_success "All steps completed successfully!"
echo ""
echo "📁 Output Locations:"
echo "   ├─ Real Data:       data/real/ (${REAL_COUNT} files)"
echo "   ├─ Synthetic Data:  data/synthetic/ (${SYN_COUNT} files)"
echo "   ├─ Results CSV:     02_clustering/clustering_results_final.csv"
echo "   └─ Visualizations:  03_analysis/plot*.png (${PLOT_COUNT} files)"
echo ""
echo "📊 Next Steps:"
echo "   1. Review results: cat 02_clustering/clustering_results_final.csv | head"
echo "   2. View plots: open 03_analysis/plot*.png"
echo "   3. Inspect notebook: jupyter notebook 03_analysis/03_plot_results.ipynb"
echo ""

# Calculate total runtime
PIPELINE_END_TIME=$(date +%s)
TOTAL_RUNTIME=$((PIPELINE_END_TIME - PIPELINE_START_TIME))
HOURS=$((TOTAL_RUNTIME / 3600))
MINUTES=$(((TOTAL_RUNTIME % 3600) / 60))
SECONDS=$((TOTAL_RUNTIME % 60))

if [ $HOURS -gt 0 ]; then
    log_success "Total pipeline runtime: ${HOURS}h ${MINUTES}m ${SECONDS}s (${TOTAL_RUNTIME} seconds)"
elif [ $MINUTES -gt 0 ]; then
    log_success "Total pipeline runtime: ${MINUTES}m ${SECONDS}s (${TOTAL_RUNTIME} seconds)"
else
    log_success "Total pipeline runtime: ${SECONDS}s"
fi

log_info "Pipeline execution log saved to: pipeline.log"
echo "=============================================================================="
