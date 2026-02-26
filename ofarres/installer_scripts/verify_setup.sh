#!/bin/bash
# ==============================================================================
# Environment Verification Script
# ==============================================================================
# Purpose: Verify all dependencies are installed before running the pipeline
# Usage:   ./verify_setup.sh
# ==============================================================================

echo "=========================================="
echo "Environment Verification"
echo "=========================================="
echo ""

ALL_OK=true

# Check R installation
echo "Checking R..."
if command -v Rscript &> /dev/null; then
    R_VERSION=$(Rscript --version 2>&1 | head -n1)
    echo "  ✓ R found: $R_VERSION"
else
    echo "  ✗ R not found"
    ALL_OK=false
fi

# Check Python installation
echo "Checking Python..."
# Prefer 'python' (conda env) over 'python3' (system)
if command -v python &> /dev/null; then
    PYTHON_VERSION=$(python --version)
    echo "  ✓ Python found: $PYTHON_VERSION"
    PYTHON_CMD="python"
elif command -v python3 &> /dev/null; then
    PYTHON_VERSION=$(python3 --version)
    echo "  ✓ Python found: $PYTHON_VERSION"
    PYTHON_CMD="python3"
else
    echo "  ✗ Python not found"
    ALL_OK=false
fi

echo ""

# Check R packages
if command -v Rscript &> /dev/null; then
    echo "Checking R packages..."
    Rscript -e "
    pkgs <- c('jsonlite', 'mvtnorm', 'synthpop', 'dplyr')
    installed <- pkgs %in% installed.packages()[,'Package']
    for (i in seq_along(pkgs)) {
        if (installed[i]) {
            cat('  ✓', pkgs[i], '\n')
        } else {
            cat('  ✗', pkgs[i], '- MISSING\n')
        }
    }
    if (!all(installed)) quit(status=1)
    " 2>/dev/null
    
    if [ $? -ne 0 ]; then
        echo "  ⚠ Some R packages missing. Run: Rscript install_requirements.R"
        ALL_OK=false
    fi
    echo ""
fi

# Check Python packages
if [ -n "$PYTHON_CMD" ]; then
    echo "Checking Python packages..."
    
    for pkg in pandas numpy sklearn matplotlib seaborn scipy joblib tqdm jupyter; do
        $PYTHON_CMD -c "import $pkg" 2>/dev/null
        if [ $? -eq 0 ]; then
            echo "  ✓ $pkg"
        else
            echo "  ✗ $pkg - MISSING"
            ALL_OK=false
        fi
    done
    
    if [ "$ALL_OK" = false ]; then
        echo "  ⚠ Some Python packages missing. Run: pip install -r ../requirements.txt"
    fi
    echo ""
fi

# Check file permissions
echo "Checking file permissions..."
if [ -x "../run_all.sh" ]; then
    echo "  ✓ run_all.sh is executable"
else
    echo "  ✗ run_all.sh is not executable"
    echo "    Run: chmod +x run_all.sh"
    ALL_OK=false
fi
echo ""

# Check config file
echo "Checking configuration..."
if [ -f "../config/config.json" ]; then
    echo "  ✓ config/config.json exists"
    
    # Validate JSON
    if [ -n "$PYTHON_CMD" ]; then
        $PYTHON_CMD -c "import json; json.load(open('../config/config.json'))" 2>/dev/null
        if [ $? -eq 0 ]; then
            echo "  ✓ config.json is valid JSON"
        else
            echo "  ✗ config.json has syntax errors"
            ALL_OK=false
        fi
    fi
else
    echo "  ✗ config/config.json not found"
    ALL_OK=false
fi
echo ""

# Check directory structure
echo "Checking directory structure..."
for dir in "../01_data_generation" "../02_clustering" "../03_analysis" "../data/real" "../data/synthetic"; do
    if [ -d "$dir" ]; then
        DIRNAME=$(basename "$dir")
        PARENTDIR=$(basename "$(dirname "$dir")")
        if [ "$PARENTDIR" = "ofarres" ]; then
            echo "  ✓ ${DIRNAME}/"
        else
            echo "  ✓ ${PARENTDIR}/${DIRNAME}/"
        fi
    else
        DIRNAME=$(basename "$dir")
        PARENTDIR=$(basename "$(dirname "$dir")")
        if [ "$PARENTDIR" = "ofarres" ]; then
            echo "  ✗ ${DIRNAME}/ missing"
        else
            echo "  ✗ ${PARENTDIR}/${DIRNAME}/ missing"
        fi
        ALL_OK=false
    fi
done
echo ""

# Final result
echo "=========================================="
if [ "$ALL_OK" = true ]; then
    echo "✅ All checks passed!"
    echo "You can now run: cd .. && ./run_all.sh"
else
    echo "❌ Some checks failed"
    echo "Please fix the issues above before running the pipeline"
    echo "See SETUP.md for installation instructions"
    exit 1
fi
echo "=========================================="
