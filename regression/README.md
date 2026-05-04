# Study II: Regression Utility Pipeline

This directory contains the large-scale simulation pipeline evaluating the inferential utility of synthetic data via **Linear** and **Logistic Regression**.

## Pipeline Execution
```bash
./run_all.sh
```

## Stages
1. **`01_generate_OD/`**: Generates original continuous and binary datasets.
2. **`02_generate_SD/`**: Generates synthetic datasets comparing CART against Parametric generators (Normal/LogReg).
3. **`03_regression_analysis/`**: Fits regression models and extracts coefficients/CIs using a parallel work-stealing queue.
4. **`04_evaluation/`**: Calculates Confidence Interval Overlap (CIO), Bias Ratio, and CI Width Ratio.