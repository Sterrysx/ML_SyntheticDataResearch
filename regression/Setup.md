# Setup Guide - Regression Pipeline

This guide covers setup for the regression pipeline in regression/.

## Prerequisites

- uv (recommended) or conda
- R 4.0 or newer
- Git

## 1. Create and activate a uv environment (recommended)

```bash
uv venv --python 3.13
source .venv/bin/activate
uv pip install -r ../requirements.txt
```

## 2. Alternative: conda environment

```bash
conda create -n synthetic_data python=3.13
conda activate synthetic_data
pip install --upgrade pip
pip install -r ../requirements.txt
```

## 3. Install R packages

The R scripts auto-install missing packages on first run. To preinstall:

```bash
Rscript -e "install.packages(c('jsonlite','mvtnorm','MultiDiscreteRNG','synthpop','arrow'), repos='https://cloud.r-project.org')"
```

## 4. Run the pipeline

```bash
cd regression
./run_all.sh
```

## 5. Outputs to expect

- data/original/OD_*.parquet
- data/synthetic/SD_*.parquet
- results/aggregated_model_metrics.parquet
- results/aggregated_logistic_metrics.parquet
- results/summary_table.csv
- results/latest_pipeline_timing.json
- results/latest_sd_timing_breakdown.json

## Troubleshooting

- no active venv: run source .venv/bin/activate before ./run_all.sh.
- conda env not found: recreate synthetic_data and activate it.
- Rscript not found: install R and ensure it is on PATH.
- Permission denied: run chmod +x run_all.sh.
- Long runtimes: reduce simulation.M and/or grid size in config.json.
