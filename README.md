# Synthetic Data Clustering Research

**Comparing K-Means and Hierarchical Clustering Performance on Synthetic Data Generated via CART Method**

---

## 📋 Project Overview

This repository contains research extending the work of **Chen Xinnuo** (2025) on evaluating how well synthetic data preserves clustering properties of real data. The extension is developed by **Oriol Farrés** under the supervision of **Jordi Cortés**. The project systematically compares K-Means and Hierarchical Clustering algorithms across 1,800+ datasets to answer critical questions about synthetic data quality.

### Research Questions
1. **Detection**: Can clustering algorithms correctly identify the number of clusters (k) in synthetic data?
2. **Robustness**: Which algorithm (K-Means vs Hierarchical) is more resilient to synthetic data distortion?
3. **Quality**: How much does the synthpop CART method degrade cluster geometry and separability?

---

## 🗂️ Repository Structure

```
ML_SyntheticDataResearch/
│
├── ofarres/                    ← 🎯 MAIN RESEARCH PIPELINE (Current Work)
│   ├── SETUP.md               ← Installation & setup instructions
│   ├── README.md              ← Complete pipeline documentation
│   ├── QUICKSTART.md          ← Quick reference guide
│   ├── run_all.sh             ← Execute complete pipeline
│   ├── requirements.txt       ← Python dependencies
│   │
│   ├── config/                ← Configuration files
│   │   └── config.json        ← Simulation parameters
│   │
│   ├── installer_scripts/     ← Installation & verification scripts
│   │   ├── install_all.sh     ← Master installer (R + Python packages)
│   │   ├── install_requirements.R
│   │   └── verify_setup.sh    ← Dependency checker
│   │
│   ├── 01_data_generation/    ← Step 1: Generate real + synthetic data
│   ├── 02_clustering/         ← Step 2: Run clustering simulations
│   └── 03_analysis/           ← Step 3: Visualization & results
│
├── xinnuo_files/              ← 📚 ORIGINAL RESEARCH (Chen Xinnuo, 2025)
│
├── README.md                  ← This file
└── .gitignore
```

---

## 🚀 Quick Start

### 1. Prerequisites

- **Conda** (Miniconda or Anaconda)
- **R** (version 4.0+)
- **Python** (3.9+)

### 2. Environment Setup

```bash
# Create and activate conda environment
conda create -n synthetic_data python=3.9
conda activate synthetic_data

# Navigate to the pipeline directory
cd ofarres/

# Install all dependencies (R + Python)
./installer_scripts/install_all.sh
```

### 3. Run the Pipeline

```bash
# Execute the complete pipeline (~15-20 minutes)
./run_all.sh
```

**Outputs**: 1,800 datasets analyzed → 3 publication-ready visualizations in `03_analysis/`

---

## 📖 Documentation

Comprehensive documentation is available in the **[ofarres/](ofarres/)** directory:

- **[ofarres/SETUP.md](ofarres/SETUP.md)** - Complete installation and environment setup guide
- **[ofarres/README.md](ofarres/README.md)** - Full pipeline documentation and methodology

---

## 🎓 Citation & Attribution

### Original Research
This work extends the methodology developed by:

> **Chen Xinnuo** (2025). *TFG-EST: Synthetic Data Generation and Clustering Analysis*.  
> Universidad de Barcelona.

Original research materials are preserved in the **[xinnuo_files/](xinnuo_files/)** directory for reference and citation purposes.

### Current Extension
Pipeline implementation and extended analysis by **Oriol Farrés** (2025-2026), under the supervision of **Jordi Cortés**.

**Extended Work:**
- Automated 3-step pipeline implementation (R + Python)
- Systematic parameter sweep across 1,800+ dataset configurations
- Publication-ready visualization framework
- Reproducible research infrastructure with version control

---

## 🧪 What This Pipeline Does

### Step 1: Data Generation
- Generates 180 real datasets with controlled cluster properties (varying N, k, separation, correlation)
- Creates 10 synthetic versions of each real dataset using synthpop CART method
- **Output**: 1,980 total datasets (180 real + 1,800 synthetic)

### Step 2: Clustering Simulation
- Tests both K-Means and Hierarchical Clustering on all datasets
- Evaluates cluster detection accuracy, quality metrics (silhouette, distortion)
- Compares real vs synthetic performance degradation
- **Output**: `clustering_results_final.csv` (1.4 MB)

### Step 3: Visualization & Analysis
- Generates 3 key publication plots:
  1. **Success Rate**: Correct cluster detection frequency
  2. **Performance Delta**: Real vs Synthetic degradation
  3. **Quality Metrics**: Silhouette score & distortion boxplots
- **Output**: 24 publication-ready PNG visualizations

---

## 📊 Key Findings (Preview)

The pipeline reveals that:
- Synthetic data significantly degrades cluster separability
- Hierarchical clustering is more robust to synthetic distortion than K-Means
- Success rates drop 40-60% when using synthetic vs real data

*(For complete results, see [ofarres/03_analysis/](ofarres/03_analysis/))*

---

## 🛠️ Requirements

### R Packages
- `jsonlite` - Config file parsing
- `mvtnorm` - Multivariate normal generation
- `synthpop` - Synthetic data generation (CART method)
- `dplyr` - Data manipulation

### Python Packages
- `pandas`, `numpy` - Data processing
- `scikit-learn` - Clustering algorithms
- `matplotlib`, `seaborn` - Visualization
- `scipy` - Statistical analysis
- `jupyter`, `notebook` - Interactive analysis

*(Install automatically via `installer_scripts/install_all.sh`)*

---

##  License

Research code for academic purposes. Original work by Chen Xinnuo (2025), extended implementation by Oriol Farrés under supervision of Jordi Cortés (2025-2026).

---

**Last Updated**: January 3, 2026  
**Repository**: ML_SyntheticDataResearch  
**Status**: Production-ready pipeline  
**Extended by**: Oriol Farrés (supervised by Jordi Cortés)  
**Original work**: Chen Xinnuo (2025)
