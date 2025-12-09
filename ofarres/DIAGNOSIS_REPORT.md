# Simulation Pipeline Failure - Diagnosis Report

**Date:** December 9, 2025  
**Investigator:** GitHub Copilot  
**Issue:** Pipeline completing in 8 seconds instead of 10+ minutes

---

## Executive Summary

The pipeline was failing due to **two critical bugs**:

1. **R Script (Data Generation):** Cluster size calculation error causing row mismatch (ALREADY FIXED in current code)
2. **Python Script (Clustering Analysis):** Processing only 10% of synthetic datasets due to incorrect file pairing logic

**Impact:** The Python script was processing only 180 out of 1800 synthetic datasets, completing in seconds rather than minutes and producing incomplete/incorrect results.

---

## Issue #1: R Script - Data Generation Row Mismatch ✅ ALREADY FIXED

### The Problem
When `N` (total sample size) is not evenly divisible by `k` (number of clusters), the original code would create:
- Clusters with `floor(N/k)` samples each
- Leftover samples creating a mismatch

**Example:**
- N=250, k=3
- `floor(250/3) = 83` samples per cluster
- Total: 83 × 3 = **249 samples** (missing 1!)
- Error: `replacement has 250 rows, data has 249`

### The Fix (Lines 70-79)
```r
# EXACT SIZE CALCULATION
base_n <- floor(N / k)
remainder <- N %% k
cluster_sizes <- rep(base_n, k)

if (remainder > 0) {
  cluster_sizes[1:remainder] <- cluster_sizes[1:remainder] + 1
}
```

**How it works:**
- Calculates base cluster size and remainder
- Distributes remainder evenly to first few clusters
- For N=250, k=3: clusters get [84, 83, 83] = exactly 250 samples

**Status:** ✅ This fix is already implemented in the current code. The error logs you saw were from an older version before the fix.

---

## Issue #2: Python Script - Incorrect File Pairing Logic ❌ CRITICAL BUG

### The Problem

**Original Code (Lines 79-84):**
```python
real_files = sorted(glob.glob("../data/real/*.csv"))
syn_files = sorted(glob.glob("../data/synthetic/*.csv"))

print(f"Found {len(real_files)} datasets. Starting processing...")

results = Parallel(n_jobs=-1)(
    delayed(process_dataset_pair)(r_path, s_path) 
    for r_path, s_path in tqdm(zip(real_files, syn_files), total=len(real_files))
)
```

**Why This Was Wrong:**

1. **File Count Mismatch:**
   - Real files: 180
   - Synthetic files: 1,800 (10 per real file, per config `m=10`)
   - `zip()` stops at shortest list → only 180 pairs processed
   - **1,620 synthetic datasets ignored!**

2. **Incorrect Pairing:**
   ```
   Real:     N250_p10_k2_rho0_sep0.1_rep1.csv
   Synthetic files (10 total):
     - N250_p10_k2_rho0_sep0.1_rep1_syn1.csv
     - N250_p10_k2_rho0_sep0.1_rep1_syn2.csv
     - ...
     - N250_p10_k2_rho0_sep0.1_rep1_syn10.csv
   ```
   
   The `zip()` paired real file #1 with synthetic file #1, real file #2 with synthetic file #2, etc., completely ignoring the naming convention and m=10 relationship.

3. **Incorrect Metrics:**
   - Calculated Silhouette on **Real** data: `sil_real, _ = calculate_quality_metrics(X_real, ...)`
   - But we need to evaluate **Synthetic** data quality!
   - Detection was correct (ran on synthetic), but quality metrics were backwards

### The Fix

**New Logic:**
- Process **each synthetic file independently**
- Parse synthetic filename to find its corresponding real file
- Calculate metrics on **synthetic data** and compare to real baseline

**Key Changes:**

1. **Process All Synthetic Files:**
```python
syn_files = sorted(glob.glob("../data/synthetic/*.csv"))
results = Parallel(n_jobs=-1)(
    delayed(process_synthetic_dataset)(syn_path) 
    for syn_path in tqdm(syn_files, desc="Processing")
)
```

2. **Proper Filename Parsing:**
```python
# Example: "N250_p10_k2_rho0_sep0.1_rep1_syn5.csv"
match = re.match(r'N(\d+)_p(\d+)_k(\d+)_rho([\d.]+)_sep([\d.]+)_rep(\d+)_syn(\d+)\.csv', syn_filename)
N, p, k, rho, sep, rep, syn_idx = match.groups()

# Build corresponding real filename
real_filename = f"N{N}_p{p}_k{k}_rho{rho}_sep{sep}_rep{rep}.csv"
```

3. **Correct Metric Calculation:**
```python
# Baseline: Quality of REAL data
sil_real_km, dist_real_km = calculate_quality_metrics(X_real, k_truth, 'kmeans')

# Test: Quality of SYNTHETIC data
sil_syn_km, dist_syn_km = calculate_quality_metrics(X_syn, k_truth, 'kmeans')

# Compare: How much worse is synthetic?
'diff_sil_km': sil_syn_km - sil_real_km
```

---

## Expected Performance Impact

### Before Fix:
- Processes: 180 datasets
- Time: ~8 seconds
- Coverage: 10% of data
- Results: Incorrect (mixed real/synthetic metrics)

### After Fix:
- Processes: 1,800 datasets
- Time: ~10-15 minutes (10x increase)
- Coverage: 100% of data
- Results: Correct (all synthetic evaluated)

---

## Testing Recommendations

1. **Verify File Counts:**
   ```bash
   ls ../data/real/*.csv | wc -l      # Should be 180
   ls ../data/synthetic/*.csv | wc -l  # Should be 1800
   ```

2. **Check Output:**
   - `clustering_results_final.csv` should have 1,800 rows (one per synthetic dataset)
   - Check that `syn_idx` ranges from 1-10 for each real dataset

3. **Validate Metrics:**
   - Success rates should vary based on separation/k
   - Quality differences should be positive (synthetic usually worse than real)

---

## Files Modified

1. ✅ **01_generate_synthetic.R** - No changes needed (already fixed)
2. ✅ **02_run_simulation.py** - Complete rewrite of file pairing and metric logic

---

## Conclusion

The 8-second runtime was caused by the Python script processing only 10% of the data due to incorrect `zip()` pairing of files. The fix processes all 1,800 synthetic datasets correctly, matching each to its real baseline, and properly calculates detection and quality metrics on the synthetic data.

The R script error mentioned in the logs was from an older version and has already been resolved in the current codebase.
