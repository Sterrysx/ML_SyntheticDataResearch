#!/bin/bash
# ==============================================================================
# Complete Installation Script
# ==============================================================================
# Purpose: Install all R and Python dependencies for the pipeline
# Prerequisites: Conda environment must be created and activated first
# Usage: ./install_all.sh
# ==============================================================================

set -e  # Exit on any error

# Resolve script directory and repository root so the script works
# no matter where it's executed from (keeps path logic "intelligent").
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"


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
    Rscript "$REPO_ROOT/installer_scripts/install_requirements.R"
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

if [ -f "$REPO_ROOT/requirements.txt" ]; then
    pip install -r "$REPO_ROOT/requirements.txt"
    if [ $? -ne 0 ]; then
        echo ""
        echo "❌ Python package installation failed!"
        exit 1
    fi
else
    echo "❌ ERROR: requirements.txt not found in repository root ($REPO_ROOT)!"
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

# Make sure the main runner is executable (harmless if already executable)
if [ -f "$REPO_ROOT/run_all.sh" ]; then
    chmod +x "$REPO_ROOT/run_all.sh" || true
fi

# Run the verifier from the installer_scripts directory so its relative
# path checks resolve consistently (the verifier expects to be run there).
pushd "$REPO_ROOT/installer_scripts" > /dev/null
./verify_setup.sh
VERIFY_EXIT=$?
popd > /dev/null

if [ $VERIFY_EXIT -eq 0 ]; then
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
    exit $VERIFY_EXIT
fi
