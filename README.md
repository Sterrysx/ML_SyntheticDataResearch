# Downstream Utility of Synthetic Data

**Structural and Inferential Evaluation in Clustering and Regression**

---

## 📋 Project Overview

This repository contains the comprehensive research codebase evaluating the downstream utility of synthetic data generated via the CART method. The project systematically evaluates synthetic data across two main analytical domains: **Clustering** (structural utility) and **Regression** (inferential utility). 

This work extends the preliminary research of **Chen Xinnuo** (2025) and was developed by **Oriol Farrés** under the supervision of **Jordi Cortés** (2025-2026).

## 🗂️ Repository Structure

The project is organized into two independent simulation pipelines:

```
ML_SyntheticDataResearch/
│
├── SETUP.md                   ← Master environment setup instructions
├── README.md                  ← This file
│
├── clustering/                ← 🎯 STUDY I: CLUSTERING PIPELINE
│   ├── README.md              ← Clustering pipeline documentation
│   ├── run_all.sh             ← Execute clustering pipeline
│   └── ...                    
│
├── regression/                ← 🎯 STUDY II: REGRESSION PIPELINE
│   ├── README.md              ← Regression pipeline documentation
│   ├── run_all.sh             ← Execute regression pipeline
│   └── ...                    
│
└── xinnuo_files/              ← 📚 ORIGINAL EXPLORATORY RESEARCH (Chen Xinnuo, 2025)
```

## 🚀 Quick Start

1. Follow the installation guide in **[SETUP.md](SETUP.md)** to configure the Conda environment (R + Python).
2. Navigate to the pipeline you wish to run:
   - For clustering: `cd clustering/ && ./run_all.sh`
   - For regression: `cd regression/ && ./run_all.sh`

## 🎓 Citation & Attribution

**Primary Authors**: Oriol Farrés (2025-2026), Supervised by Jordi Cortés.  
**Foundational Exploratory Work**: Chen Xinnuo (2025). *TFG-EST: Synthetic Data Generation and Clustering Analysis*.