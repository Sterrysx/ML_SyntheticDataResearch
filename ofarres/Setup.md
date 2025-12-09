# Synthetic Data Clustering Research Pipeline - Setup Guide

## Overview

This directory contains a complete, standalone pipeline for comparing K-Means and Hierarchical Clustering algorithms on synthetic data. This guide will help you set up the environment and run the full simulation.

---

## 📋 Prerequisites

### Required Software

1. **R** (version 4.0 or higher)
   - Download: https://www.r-project.org/
   - Check installation: `R --version`

2. **Python** (version 3.8 or higher)
   - Download: https://www.python.org/
   - Check installation: `python --version` or `python3 --version`

3. **Git** (for cloning/version control)
   - Download: https://git-scm.com/
   - Check installation: `git --version`

### System Requirements

- **RAM:** At least 4 GB (8 GB recommended)
- **Disk Space:** ~500 MB free
- **OS:** Linux, macOS, or Windows (with WSL recommended for Windows)

---

## 🚀 Quick Setup

### Step 1: Install R Packages

Run the provided R installation script:

```bash
Rscript install_requirements.R
```

**Manual installation (if needed):**
```r
install.packages(c("jsonlite", "mvtnorm", "synthpop", "dplyr"), 
                 repos = "http://cran.us.r-project.org")
```

**Required R packages:**
- `jsonlite` - JSON parsing
- `mvtnorm` - Multivariate normal distribution
- `synthpop` - Synthetic data generation (CART method)
- `dplyr` - Data manipulation

---

### Step 2: Install Python Packages

**Option A: Using pip (recommended)**
```bash
pip install -r requirements.txt
```

**Option B: Using conda**
```bash
conda create -n synthetic_data python=3.9
conda activate synthetic_data
pip install -r requirements.txt
```

**Option C: Manual installation**
```bash
pip install pandas numpy scikit-learn matplotlib seaborn scipy joblib tqdm jupyter
```

**Core Python packages:**
- `pandas` - Data manipulation
- `numpy` - Numerical computing
- `scikit-learn` - Clustering algorithms
- `matplotlib` + `seaborn` - Visualization
- `scipy` - Statistical tests
- `joblib` - Parallel processing
- `tqdm` - Progress bars
- `jupyter` - Notebook execution

---

### Step 3: Verify Installation

Run the verification script:

```bash
# Verify R packages
Rscript -e "lapply(c('jsonlite', 'mvtnorm', 'synthpop', 'dplyr'), library, character.only=TRUE); cat('✅ All R packages installed\n')"

# Verify Python packages
python -c "import pandas, numpy, sklearn, matplotlib, seaborn, scipy, joblib, tqdm, jupyter; print('✅ All Python packages installed')"
```

---

## 🎯 Running the Pipeline

### Option 1: Run Everything (Recommended)

Execute the master script that runs all 3 steps automatically:

```bash
./run_all.sh
```

**Expected runtime:** 15-20 minutes

**This will:**
1. Generate 180 real + 1,800 synthetic datasets (~3 min)
2. Run clustering analysis on all 1,800 datasets (~12 min)
3. Generate 3 publication-ready plots (~1 min)

---

### Option 2: Run Steps Individually

#### Step 1: Generate Data (R)
```bash
cd 01_data_generation
Rscript 01_generate_synthetic.R
cd ..
```

**Output:** `data/real/` (180 files), `data/synthetic/` (1,800 files)

---

#### Step 2: Run Clustering Simulation (Python)
```bash
cd 02_clustering
python 02_run_simulation.py
cd ..
```

**Output:** `02_clustering/clustering_results_final.csv` (1,800 rows)

---

#### Step 3: Generate Visualizations (Jupyter)
```bash
cd 03_analysis
jupyter notebook 03_plot_results.ipynb
```

**Or execute programmatically:**
```bash
cd 03_analysis
python -m jupyter nbconvert --to notebook --execute 03_plot_results.ipynb --output 03_plot_results_executed.ipynb
cd ..
```

**Output:** 3 PNG plots in `03_analysis/`

---

## 🔧 Configuration

### Modifying Simulation Parameters

Edit `config.json` to change simulation parameters:

```json
{
  "simulation": {
    "N": 250,              // Sample size per dataset
    "n": 5,                // Number of real datasets per combination
    "m": 10,               // Synthetic datasets per real dataset
    "random_seed_base": 12345
  },
  "parameters": {
    "p": [2, 5, 10],              // Number of features
    "k": [2, 3, 4],               // True number of clusters
    "separation": [0.1, 2, 6, 10], // Cluster separation distance
    "rho": [0.0]                   // Within-cluster correlation
  }
}
```

**Important:** After changing parameters, delete old data:
```bash
rm -rf data/real/* data/synthetic/* 02_clustering/clustering_results_final.csv
```

---

## 📊 Expected Outputs

### Data Files
```
data/
├── real/                    # 180 real datasets
│   └── N250_p10_k2_rho0_sep0.1_rep1.csv
└── synthetic/               # 1,800 synthetic datasets
    └── N250_p10_k2_rho0_sep0.1_rep1_syn5.csv
```

### Results
```
02_clustering/
└── clustering_results_final.csv    # 1,800 rows with metrics
```

### Visualizations
```
03_analysis/
├── plot1_success_rate_comparison.png      # Heatmap: Detection success
├── plot2_performance_delta.png            # Heatmap: HC vs K-Means
└── plot3_quality_distortion_boxplots.png  # Boxplots: Quality degradation
```

---

## 🐛 Troubleshooting

### R Package Installation Fails

**Problem:** Package compilation errors

**Solution:**
```bash
# On Ubuntu/Debian
sudo apt-get install r-base-dev libcurl4-openssl-dev libssl-dev libxml2-dev

# On macOS
brew install gcc
```

---

### Python Import Errors

**Problem:** `ModuleNotFoundError: No module named 'X'`

**Solution:**
```bash
# Ensure you're in the correct environment
which python
pip list

# Reinstall missing package
pip install X
```

---

### Permission Denied on run_all.sh

**Problem:** `bash: ./run_all.sh: Permission denied`

**Solution:**
```bash
chmod +x run_all.sh
./run_all.sh
```

---

### Out of Memory Error

**Problem:** Python script crashes with memory error

**Solution:**
- Reduce `m` in config.json (fewer synthetic datasets)
- Close other applications
- Run on a machine with more RAM

---

### Data Generation Takes Too Long

**Problem:** R script runs for hours

**Solution:**
- Check parameters in config.json
- Reduce grid size (fewer p, k, separation values)
- Ensure synthpop is using CART (not ctree which is slower)

---

### Plots Not Generated

**Problem:** Jupyter notebook doesn't create PNG files

**Solution:**
```bash
# Install missing dependencies
pip install nbconvert

# Or run notebook interactively
cd 03_analysis
jupyter notebook 03_plot_results.ipynb
# Then: Cell → Run All
```

---

## 🔍 Validation Tests

### Verify Data Generation
```bash
# Check file counts
ls data/real/*.csv | wc -l      # Should be 180
ls data/synthetic/*.csv | wc -l  # Should be 1800

# Check row count in a file (should be N+1 for header)
head -n 251 data/real/N250_p2_k2_rho0_sep0.1_rep1.csv | wc -l  # Should be 251
```

### Verify Simulation Results
```python
import pandas as pd
df = pd.read_csv('02_clustering/clustering_results_final.csv')

# Should have 1800 rows
assert len(df) == 1800, f"Expected 1800 rows, got {len(df)}"

# Should have 10 synthetic per real dataset
counts = df.groupby(['N','p','k','sep','rho','rep']).size()
assert counts.min() == 10 and counts.max() == 10

print("✅ Validation passed!")
```

---

## 🌐 Environment-Specific Notes

### Windows (WSL)

If using Windows Subsystem for Linux:

```bash
# Install R
sudo apt update
sudo apt install r-base

# Install Python
sudo apt install python3 python3-pip

# Then follow normal setup
```

---

### macOS

Using Homebrew:

```bash
# Install R
brew install r

# Install Python
brew install python

# Then follow normal setup
```

---

### Linux (Ubuntu/Debian)

```bash
# Install R
sudo apt update
sudo apt install r-base

# Install Python
sudo apt install python3 python3-pip

# Then follow normal setup
```

---

## 📚 Additional Resources

- **Full Documentation:** See `README.md`
- **Quick Reference:** See `QUICKSTART.md`
- **Architecture Details:** See `PROJECT_STRUCTURE.md`
- **Bug Fixes:** See `DIAGNOSIS_REPORT.md`

---

## ✅ Setup Checklist

- [ ] R installed and accessible (`R --version`)
- [ ] Python installed and accessible (`python --version`)
- [ ] R packages installed (`jsonlite`, `mvtnorm`, `synthpop`, `dplyr`)
- [ ] Python packages installed (see `requirements.txt`)
- [ ] `run_all.sh` is executable (`chmod +x run_all.sh`)
- [ ] `config.json` exists and is valid JSON
- [ ] All verification tests pass

---

## 🎓 Next Steps

1. **Run the pipeline:** `./run_all.sh`
2. **Review results:** Check `clustering_results_final.csv`
3. **View plots:** Open PNG files in `03_analysis/`
4. **Modify parameters:** Edit `config.json` and re-run
5. **Read analysis:** Open `03_analysis/03_plot_results.ipynb` in Jupyter

---

## 📧 Support

For issues:
1. Check this SETUP.md first
2. Review `DIAGNOSIS_REPORT.md` for known bugs
3. Check `QUICKSTART.md` for common commands
4. Inspect terminal output for specific errors

---

**Last Updated:** December 9, 2025  
**Pipeline Version:** 2.0  
**Status:** Production Ready
