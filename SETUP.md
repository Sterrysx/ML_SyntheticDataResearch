# Environment Setup Guide

This project utilizes both **R** and **Python** to generate synthetic data and evaluate machine learning models. 

## 1. Prerequisites
- **Conda** (Miniconda or Anaconda)
- **R** (version 4.0+)
- **Python** (3.9+)

## 2. Installation

Create and activate the conda environment:
```bash
conda create -n synthetic_data python=3.9
conda activate synthetic_data
```

Install required Python packages:
```bash
pip install pandas numpy scikit-learn matplotlib seaborn scipy pyarrow fastparquet
```

Install required R packages (run within an R console):
```R
install.packages(c("synthpop", "jsonlite", "mvtnorm", "dplyr", "parallel", "arrow"))
```