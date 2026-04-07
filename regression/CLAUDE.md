# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Running the Pipeline

```bash
cd regression
./run_all.sh          # full run (all 5 steps)
```

For a quick smoke-test, temporarily set `"M": 10` in `config/config.json` before running. The script detects config changes via MD5 hash and automatically nukes `data/` and `results/` before regenerating — so always restore M when done.

Individual steps can be run manually:
```bash
# OD generation (accepts "continuous" or "binary" to run one var_type only)
cd 01_generate_OD
Rscript 01_generate_original_data_linear.R
Rscript 01_generate_original_data_logistic.R continuous

# SD generation
cd 02_generate_SD
Rscript 02_generate_synthetic_data.R

# Regression analysis
cd 03_regression_analysis
python3 03_regression_linear.py
python3 03_regression_logistic.py
```

R packages are auto-installed from CRAN on first run. Python dependencies are in `requirements.txt`; the expected conda env is `synthetic_data`.

## Architecture Overview

The pipeline is a five-step simulation study comparing synthpop synthesis methods across linear and logistic regression. Steps are orchestrated by `run_all.sh` with smart resume: each step checks whether its output files already exist and skips if complete.

```
config.json → run_all.sh → Step 1a (linear OD) → Step 1b (logistic OD)
                                                           ↓
                                              Step 2 (synthetic data)
                                                           ↓
                                  Step 3a (OLS) + Step 3b (Logit)
                                                           ↓
                              results/aggregated_model_metrics.parquet
                              results/aggregated_logistic_metrics.parquet
```

### config.json is the single source of truth

All scenario dimensions live here. `run_all.sh` reads it directly via inline Python to compute expected file counts. The R scripts read it via `jsonlite::fromJSON`. The Python regression scripts load it for paths but don't use simulation parameters directly.

Key parameters: `N`, `p`, `rho`, `sigma_2` (shared by both linear and logistic — flat under `parameters`, not nested), `beta` (intercept + 10 coefficients, sliced to `p+1` per scenario), `continuous_methods`, `binary_methods`.

### Parallelism model (R scripts)

All three R scripts use an identical **work-stealing queue** pattern:

1. Master builds a directory of task RDS files (`tempdir()/od_linear_tasks_todo/task_00001.rds`, …)
2. Master writes a worker script to `tempdir()` and launches `NUM_CORES` independent `Rscript` processes via `system2(wait=FALSE)`
3. Each worker loops: `file.rename(todo/task.rds, doing/task.rds)` — atomic at the OS level, so no locks needed
4. Worker processes task, writes output parquet, removes its file from `doing/`
5. Master polls progress files every 0.5s and redraws an ASCII progress bar

Workers are fully independent processes — they communicate only via the filesystem. The worker code itself is generated as a character vector and written to a temp file; this is intentional (no shared memory, clean subprocess isolation).

`NUM_CORES` is `max(1, detectCores() - 6)` in the OD scripts and hardcoded to 18 in the SD script.

### Parallelism model (Python scripts)

`multiprocessing.Pool(N_WORKERS=18)` with `pool.imap_unordered`. Each worker receives one SD task, reads its own OD+SD parquet pair, and returns a list of result dicts (one per iteration). Low per-process memory because parquets are read inside the worker, not in the master.

### DGP (Data Generating Process)

**Linear:** `ε ~ N(0, σ²)`, `y = X_design · β + ε`

**Logistic:** `ε ~ N(0, σ²)`, `lp = X_design · β + ε`, `prob = 1/(1 + exp(-lp))`, `y ~ Bernoulli(prob)` — this is the **latent variable formulation**: epsilon enters the linear predictor before the logistic link, not as additive noise on the outcome. Both branches use the same sigma_2 values from config.

**Binary X generation:** Correlated binary X is pre-computed via `MultiDiscreteRNG::simBinaryCorr.B()` for each unique (p, rho, p1) combination and cached to disk before workers launch. Workers load this cache and call `genB()` per iteration. The cache is always built in the master process, not in workers.

**Note:** The linear DGP in the current codebase normalises beta (divides by `sqrt(lp_var)`) before generating y. This is intentional — the non-orthogonality of SNR across (p, rho, sigma_2) is a known design choice disclosed in the thesis. Do not remove this normalisation.

### File naming conventions

**OD files:**
```
OD_{linear|logistic}_{continuous|binary}_N{N}_p{p}_rho{rho}_sig{sigma}_p1{p1|NA}.parquet
```

**SD files:**
```
SD_X{XMETHOD}_Y{YMETHOD}_{linear|logistic}_{continuous|binary}_N{N}_p{p}_rho{rho}_sig{sigma}_p1{p1|NA}.parquet
```

X and Y methods are always uppercase in the filename. The scenario portion of the SD filename is the OD filename with `OD_` stripped.

### Arm enumeration (synthesis)

For each OD file, `get_arms()` in `02_generate_synthetic_data.R` determines which method combinations to synthesise. It reads the OD filename to infer types:

| Case | Regression | X type | Y type | Arms | Example |
|------|-----------|--------|--------|------|---------|
| A | Linear | Continuous | Continuous | 2 | XCART_YCART, XNORM_YNORM |
| B | Logistic | Binary | Binary | 2 | XCART_YCART, XLOGREG_YLOGREG |
| C | Linear | Binary | Continuous | 4 | full 2×2 factorial |
| D | Logistic | Continuous | Binary | 4 | full 2×2 factorial |

Cases A and B pair methods (same method for X and Y); Cases C and D take the full cross-product of `bin_methods × cont_methods`.

### OD/SD matching in regression scripts

The Python regression scripts index OD files into `od_map` keyed by `(var_type, N, p, rho, sig_str)`. SD files are matched to their OD via the same 5-tuple extracted by regex from the SD filename. The regex captures `sig` (not `gam`) for both linear and logistic files — the logistic script was updated to match this convention.

### Output schema

**aggregated_model_metrics.parquet (linear):**
`x_method, y_method, var_type, N, p, rho, sigma_2, iter, fit_status_od, fit_status_sd, adj_r2_od, adj_r2_sd, beta_od_0..10, beta_sd_0..10, se_od_0..10, se_sd_0..10, pval_od_0..10, pval_sd_0..10`

**aggregated_logistic_metrics.parquet (logistic):**
Same but with `pseudo_r2` instead of `adj_r2`, plus `accuracy_od, accuracy_sd`. Logistic fits with separation, Hessian inversion failure, or |β| > 20 are recorded as NaN rows with a status string.

## What not to change

- The work-stealing queue logic, parallelism architecture, progress display, parquet I/O, and file discovery logic in any script
- The `file.rename` atomicity pattern — this is what prevents race conditions
- The SD filename convention — both the R synthesis script and Python regression scripts parse it via regex; they must stay in sync
- The `sigma_2` key path in config (`parameters.sigma_2`, flat, not nested under `linear` or `logistic`)
