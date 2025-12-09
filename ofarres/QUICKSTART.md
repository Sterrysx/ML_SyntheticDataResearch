# Quick Reference Card

## 🚀 Run Complete Pipeline
```bash
cd ofarres
./run_all.sh
```
**Time:** 15-20 minutes  
**Output:** 1800 datasets analyzed, 3 publication plots

---

## 📁 File Structure Quick Map
```
ofarres/
├── config.json              ← Edit parameters here
├── run_all.sh               ← Run this to execute everything
├── README.md                ← Full documentation
│
├── 01_data_generation/
│   └── 01_generate_synthetic.R
│
├── 02_clustering/
│   ├── 02_run_simulation.py
│   ├── clustering_utils.py
│   └── clustering_results_final.csv  ← Results appear here
│
└── 03_analysis/
    ├── 03_plot_results.ipynb
    ├── plot1_success_rate_comparison.png  ← Outputs
    ├── plot2_performance_delta.png
    └── plot3_quality_distortion_boxplots.png
```

---

## ⚙️ Configuration Parameters

**Edit `config.json` to change simulation:**

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `N` | 250 | Sample size per dataset |
| `n` | 5 | Real datasets per combo |
| `m` | 10 | Synthetic per real |
| `p` | [2,5,10] | Number of features |
| `k` | [2,3,4] | True clusters |
| `separation` | [0.1,2,6,10] | Cluster distance |
| `rho` | [0.0] | Within-cluster correlation |

**Total datasets:** `p × k × sep × rho × n × (1 + m)`
- Default: 3×3×4×1×5×11 = 1,980 files (180 real + 1,800 synthetic)

---

## 🔧 Individual Commands

### Step 1: Generate Data
```bash
cd 01_data_generation
Rscript 01_generate_synthetic.R
```
**Output:** `../data/real/` (180 files), `../data/synthetic/` (1,800 files)

### Step 2: Run Clustering
```bash
cd 02_clustering
python 02_run_simulation.py
```
**Output:** `clustering_results_final.csv` (1,800 rows)

### Step 3: Generate Plots
```bash
cd 03_analysis
jupyter notebook 03_plot_results.ipynb
```
**Output:** 3 PNG files (plot1, plot2, plot3)

---

## 📊 Understanding the Results

### CSV Columns (`clustering_results_final.csv`)
```
N, p, k, rho, sep, rep, syn_idx          ← Parameters
success_kmeans, success_hc               ← Detection (0 or 1)
diff_sil_km, diff_sil_hc                 ← Quality delta (Syn - Real)
diff_dist_km, diff_dist_hc               ← Distance delta
sil_real_km, sil_syn_km                  ← Raw Silhouette scores
dist_real_km, dist_syn_km                ← Raw inter-cluster distances
```

### Plot Interpretations

**Plot 1: Success Rate Heatmap**
- Shows detection accuracy for K-Means and HC
- Green = high success, Red = low success
- Expect high success at high separation

**Plot 2: Performance Delta**
- Red = HC outperforms K-Means
- Blue = K-Means outperforms HC
- Identifies algorithm-specific strengths

**Plot 3: Quality Distortion**
- Negative values = synthetic worse than real (expected)
- Larger magnitude = more distortion
- Should decrease with higher separation

---

## 🐛 Common Issues

### "replacement has 250 rows, data has 249"
✅ Fixed in current code (cluster size distribution)

### Pipeline finishes in 8 seconds
❌ Old version - update `02_run_simulation.py` to process all synthetic files

### Missing Python packages
```bash
pip install -r requirements.txt
```

### Missing R packages
```r
install.packages(c("jsonlite", "mvtnorm", "synthpop", "dplyr"))
```

---

## 📈 Quick Results Check

```bash
# Verify data generation
ls data/real/*.csv | wc -l       # Should be 180
ls data/synthetic/*.csv | wc -l  # Should be 1800

# Check results
wc -l 02_clustering/clustering_results_final.csv  # Should be 1801 (including header)

# View sample results
head -5 02_clustering/clustering_results_final.csv

# Check plots
ls -lh 03_analysis/plot*.png  # Should see 3 files
```

---

## 💡 Research Questions Answered

1. **Detection:** Can algorithms find k in synthetic data?
   → See Plot 1 (success rates by separation)

2. **Comparison:** Which algorithm is more robust?
   → See Plot 2 (HC vs K-Means delta)

3. **Quality:** How much distortion does synthpop cause?
   → See Plot 3 (Silhouette score degradation)

---

## 📞 Getting Help

1. Check `README.md` for full documentation
2. Review `DIAGNOSIS_REPORT.md` for debugging history
3. Inspect logs during execution for specific errors

---

**Last Updated:** December 9, 2025
