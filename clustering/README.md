# Study I: Clustering Utility Pipeline

This directory contains the simulation pipeline for evaluating the structural utility of synthetic data using **K-Means** and **Hierarchical Clustering**.

## Pipeline Execution
```bash
./run_all.sh
```

## Stages
1. **`01_generate_OD/`**: Generates original datasets based on parameters in `config/`.
2. **`02_generate_SD/`**: Generates synthetic datasets using the CART method.
3. **`03_clustering_analysis/`**: Applies K-Means and Hierarchical Clustering to both OD and SD.
4. **`04_evaluation/`**: Evaluates clustering performance (Success Rate, Gini, Centroid Distance, Variance Difference) and generates visualizations.