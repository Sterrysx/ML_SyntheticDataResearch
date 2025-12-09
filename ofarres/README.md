# Synthetic Data Clustering Research Pipeline

**Complete 3-Step ML Research Pipeline for Comparing K-Means and Hierarchical Clustering on Synthetic Data**

---

## 📋 Overview

This pipeline systematically evaluates how well synthetic data (generated via `synthpop` CART method) preserves clustering properties of real data. It addresses the critical research question: **Can we trust synthetic data for clustering analysis?**

### Research Questions
1. **Detection:** Can algorithms correctly identify the number of clusters (k) in synthetic data?
2. **Comparison:** Which algorithm (K-Means vs Hierarchical) is more robust to synthetic data distortion?
3. **Quality:** How much does synthpop degrade cluster geometry and separability?

---

## 🏗️ Pipeline Architecture

```
ofarres/
├── config.json                           # ⚙️ Central configuration (single source of truth)
├── run_all.sh                            # 🚀 Master execution script
├── requirements.txt                      # 📦 Python dependencies
│
├── 01_data_generation/
│   └── 01_generate_synthetic.R          # Step 1: Generate real + synthetic data
│
├── 02_clustering/
│   ├── 02_run_simulation.py             # Step 2: Run clustering experiments
│   └── clustering_utils.py              # Helper functions (detection + quality)
│
├── 03_analysis/
│   └── 03_plot_results.ipynb            # Step 3: Generate 3 visualizations
│
└── data/
    ├── real/                             # 180 real datasets
    └── synthetic/                        # 1800 synthetic datasets (10 per real)
```

---

## ⚙️ Configuration (`config.json`)

**Single source of truth for all simulation parameters**

```json
{
  "simulation": {
    "N": 250,                  // Sample size per dataset
    "n": 5,                    // Number of real datasets per combination
    "m": 10,                   // Number of synthetic datasets per real dataset
    "random_seed_base": 12345  // Base seed for reproducibility
  },
  "parameters": {
    "p": [2, 5, 10],                    // Dimensions (features)
    "k": [2, 3, 4],                     // True number of clusters
    "separation": [0.1, 2, 6, 10],      // Cluster separation distance
    "rho": [0.0]                        // Within-cluster correlation
  }
}
```

**Total combinations:** 3 (p) × 3 (k) × 4 (sep) × 1 (rho) = **36 combinations**
- **Real datasets:** 36 × 5 = **180 files**
- **Synthetic datasets:** 180 × 10 = **1,800 files**

---

## 🔧 Step 1: Data Generation (R)

**Script:** `01_data_generation/01_generate_synthetic.R`

### Purpose
Generate multivariate normal clusters with controlled properties, then synthesize them using CART.

### Critical Fix Implemented
**Problem:** Original code failed when N ÷ k had a remainder (e.g., N=250, k=3)
```r
# ❌ OLD (broken): floor(250/3) × 3 = 249 rows
cluster_sizes <- rep(floor(N/k), k)

# ✅ NEW (fixed): Distribute remainder evenly
base_n <- floor(N / k)
remainder <- N %% k
cluster_sizes <- rep(base_n, k)
if (remainder > 0) {
  cluster_sizes[1:remainder] <- cluster_sizes[1:remainder] + 1
}
# For N=250, k=3: [84, 83, 83] = exactly 250 rows ✓
```

### Outputs
- **Real data:** `data/real/N250_p10_k2_rho0_sep0.1_rep1.csv`
- **Synthetic data:** `data/synthetic/N250_p10_k2_rho0_sep0.1_rep1_syn5.csv`

### Execution
```bash
cd 01_data_generation
Rscript 01_generate_synthetic.R
```

**Time:** ~2-5 minutes

---

## 🧪 Step 2: Clustering Simulation (Python)

**Script:** `02_clustering/02_run_simulation.py`  
**Utilities:** `02_clustering/clustering_utils.py`

### Purpose
Evaluate every synthetic dataset against its real counterpart using two tests.

### Two-Part Logic

#### 🎯 Test 1: Detection (Success Rate)
**Question:** Can the algorithm find the correct k?

```python
# Let algorithm guess k (2 to 10) using Silhouette Score
k_found = detect_optimal_k(X_synthetic, method='kmeans')
success = 1 if k_found == k_truth else 0
```

#### 📏 Test 2: Quality (Metric Distortion)
**Question:** How distorted is the cluster geometry?

```python
# Force algorithm to use true k
sil_real, dist_real = calculate_quality_metrics(X_real, k_truth, 'kmeans')
sil_syn, dist_syn = calculate_quality_metrics(X_synthetic, k_truth, 'kmeans')

# Calculate degradation
diff_sil = sil_syn - sil_real  # Negative = synthetic worse
diff_dist = dist_syn - dist_real
```

### Critical Fix Implemented
**Problem:** Original code only processed 180 datasets (10% of data)

```python
# ❌ OLD (broken): zip stops at shortest list (180)
results = [process(r, s) for r, s in zip(real_files, syn_files)]

# ✅ NEW (fixed): Process all 1800 synthetic files
results = [process_synthetic_dataset(syn) for syn in syn_files]

def process_synthetic_dataset(syn_path):
    # Parse: N250_p10_k2_rho0_sep0.1_rep1_syn5.csv
    # Find:  N250_p10_k2_rho0_sep0.1_rep1.csv
    real_path = derive_real_filename(syn_path)
    # ... calculate metrics on SYNTHETIC data
```

### Output
`02_clustering/clustering_results_final.csv` (1,800 rows)

**Columns:**
- Parameters: `N, p, k, rho, sep, rep, syn_idx`
- Detection: `success_kmeans, success_hc`
- Quality: `diff_sil_km, diff_sil_hc, diff_dist_km, diff_dist_hc`
- Raw metrics: `sil_real_km, sil_syn_km, dist_real_km, dist_syn_km`

### Execution
```bash
cd 02_clustering
python 02_run_simulation.py
```

**Time:** ~10-15 minutes (1800 datasets, parallel processing)

---

## 📊 Step 3: Analysis & Visualization (Jupyter)

**Notebook:** `03_analysis/03_plot_results.ipynb`

### Purpose
Generate publication-ready visualizations answering the three research questions.

### Deliverable 1: Success Rate Heatmap
**Question:** How often does each algorithm correctly identify k?

**Method:** Calculate success rate grouped by (separation, k)

**Output:** `plot1_success_rate_comparison.png`
- Side-by-side heatmaps (K-Means | Hierarchical)
- Rows = true k, Columns = separation
- Color = success rate (0-100%)

**Expected Pattern:**
- Low separation (0.1) → low success (clusters overlap)
- High separation (10) → high success (clusters distinct)

---

### Deliverable 2: Performance Delta Heatmap
**Question:** Where does Hierarchical outperform K-Means?

**Method:** Calculate difference: Success_HC - Success_KM

**Output:** `plot2_performance_delta.png`
- Single heatmap with diverging colormap
- Red = HC better, Blue = K-Means better, White = equal
- Identifies algorithm-specific strengths

---

### Deliverable 3: Quality Distortion Boxplots
**Question:** Does synthpop smooth/distort cluster geometry?

**Method:** Plot distribution of Silhouette Score delta (Syn - Real)

**Output:** `plot3_quality_distortion_boxplots.png`
- 4 panels (one per separation level)
- X-axis = true k, Y-axis = Silhouette delta
- Compare K-Means vs Hierarchical

**Expected Pattern:**
- Negative values = synthetic worse than real (expected)
- Larger negative = more distortion
- Should decrease with higher separation

---

## 🚀 Quick Start

### Option 1: Run Full Pipeline (Recommended)
```bash
cd ofarres
./run_all.sh
```

This will:
1. Validate environment (R, Python, packages)
2. Generate 180 real + 1800 synthetic datasets (~3 min)
3. Run clustering simulation on all 1800 (~12 min)
4. Generate 3 visualizations (~1 min)
5. Display summary statistics

**Total time:** ~15-20 minutes

---

### Option 2: Run Steps Individually

```bash
# Step 1: Generate data
cd 01_data_generation
Rscript 01_generate_synthetic.R

# Step 2: Run simulation
cd ../02_clustering
python 02_run_simulation.py

# Step 3: Analyze results
cd ../03_analysis
jupyter notebook 03_plot_results.ipynb
```

---

## 📦 Requirements

### R Packages
```r
install.packages(c("jsonlite", "mvtnorm", "synthpop", "dplyr"))
```

### Python Packages
```bash
pip install -r requirements.txt
```

Contents:
```
pandas>=1.5.0
numpy>=1.23.0
scikit-learn>=1.2.0
matplotlib>=3.6.0
seaborn>=0.12.0
joblib>=1.2.0
tqdm>=4.64.0
jupyter>=1.0.0
scipy>=1.10.0
```

---

## 🔍 Validation & Testing

### Data Generation Validation
```bash
# Check file counts
ls data/real/*.csv | wc -l      # Should be 180
ls data/synthetic/*.csv | wc -l  # Should be 1800

# Verify exact N rows (no rounding errors)
head -n 251 data/real/N250_p10_k2_rho0_sep0.1_rep1.csv | wc -l  # 251 (250 + header)
```

### Simulation Validation
```python
import pandas as pd
df = pd.read_csv('02_clustering/clustering_results_final.csv')

# Check completeness
assert len(df) == 1800, "Should have 1800 rows"

# Check synthetic replication count
assert df.groupby(['N','p','k','sep','rho','rep']).size().min() == 10

# Check success rates are valid
assert df['success_kmeans'].between(0, 1).all()
assert df['success_hc'].between(0, 1).all()
```

---

## 📈 Expected Results

### Success Rates
- **High separation (10):** >90% success for both algorithms
- **Medium separation (2-6):** 40-70% success, HC slightly better
- **Low separation (0.1):** <20% success (clusters indistinguishable)

### Quality Distortion
- **Silhouette delta:** Typically -0.05 to -0.15 (synthetic worse)
- **More distortion when:** Higher k, lower separation
- **HC vs K-Means:** Minimal difference in quality metrics

---

## 🐛 Troubleshooting

### Error: "replacement has 250 rows, data has 249"
✅ **Fixed in current code.** This was the cluster size rounding bug.

### Error: "Real file not found"
Check filename parsing regex in `02_run_simulation.py` line 23.

### Pipeline runs in 8 seconds
Old version bug - ensure you're using the fixed `02_run_simulation.py` that processes all 1800 files.

### Jupyter notebook won't execute
```bash
# Install nbconvert
pip install nbconvert

# Or run interactively
jupyter notebook 03_analysis/03_plot_results.ipynb
```

---

## 📚 Research Context

This pipeline replicates and operationalizes the exploratory work in `/notebooks/02_hierarchicalclustering.ipynb`, with critical stability fixes:

1. **Data generation:** Fixed N÷k rounding error
2. **File pairing:** Process all m synthetic datasets per real dataset
3. **Metric calculation:** Correctly evaluate synthetic data quality
4. **Reproducibility:** Single config file, automated execution

---

## 👥 Authors & Acknowledgments

**Original Research:** ML Synthetic Data Research Team  
**Pipeline Engineering:** Refactored December 2025  
**Reference:** See `DIAGNOSIS_REPORT.md` for detailed debugging history

---

## 📄 License

[Your License Here]

---

## 📧 Contact

For questions or issues, please [contact information].
