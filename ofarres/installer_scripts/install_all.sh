#!/bin/bash
# ==============================================================================
# Complete Installation Script
# ==============================================================================
# Purpose: Install all R and Python dependencies for the pipeline
# Prerequisites: Conda environment must be created and activated first
# Usage: ./install_all.sh
# ==============================================================================

set -e  # Exit on any error

echo "=========================================="
echo "Installing Pipeline Dependencies"
echo "=========================================="
echo ""

# Check if we're in a conda environment
if [ -z "$CONDA_DEFAULT_ENV" ]; then
    echo "❌ ERROR: No conda environment is active!"
    echo ""
    echo "Please create and activate a conda environment first:"
    echo "  conda create -n synthetic_data python=3.9"
    echo "  conda activate synthetic_data"
    echo ""
    exit 1
fi

echo "✓ Conda environment active: $CONDA_DEFAULT_ENV"
echo ""

# ==============================================================================
# Step 1: Install R packages
# ==============================================================================
echo "=========================================="
echo "Step 1: Installing R packages"
echo "=========================================="
echo ""

if command -v Rscript &> /dev/null; then
    Rscript installer_scripts/install_requirements.R
    if [ $? -ne 0 ]; then
        echo ""
        echo "❌ R package installation failed!"
        exit 1
    fi
else
    echo "❌ ERROR: R is not installed!"
    echo "Please install R first: https://www.r-project.org/"
    exit 1
fi

echo ""

# ==============================================================================
# Step 2: Install Python packages
# ==============================================================================
echo "=========================================="
echo "Step 2: Installing Python packages"
echo "=========================================="
echo ""

if [ -f "../requirements.txt" ]; then
    pip install -r ../requirements.txt
    if [ $? -ne 0 ]; then
        echo ""
        echo "❌ Python package installation failed!"
        exit 1
    fi
else
    echo "❌ ERROR: requirements.txt not found in repository root!"
    exit 1
fi

echo ""

# ==============================================================================
# Step 3: Verify installation
# ==============================================================================
echo "=========================================="
echo "Step 3: Verifying installation"
echo "=========================================="
echo ""

./installer_scripts/verify_setup.sh

if [ $? -eq 0 ]; then
    echo ""
    echo "=========================================="
    echo "✅ Installation complete!"
    echo "=========================================="
    echo ""
    echo "You can now run the pipeline:"
    echo "  ./run_all.sh"
    echo ""
else
    echo ""
    echo "❌ Setup verification failed!"
    echo "Please check the error messages above."
    exit 1
fi
