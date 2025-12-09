# Project Structure Overview

## Complete File Tree

```
ML_SyntheticDataResearch/
└── ofarres/                                  # Main pipeline directory
    │
    ├── config.json                           # ⚙️ Configuration (edit parameters here)
    ├── run_all.sh                            # 🚀 Master execution script (run this!)
    ├── requirements.txt                      # 📦 Python dependencies
    │
    ├── README.md                             # 📖 Complete documentation
    ├── QUICKSTART.md                         # ⚡ Quick reference guide
    ├── DIAGNOSIS_REPORT.md                   # 🔍 Debugging history & fixes
    ├── Setup.md                              # 📋 Original setup notes
    │
    ├── 01_data_generation/                   # STEP 1: Generate datasets
    │   └── 01_generate_synthetic.R           # R script (fixed cluster rounding)
    │
    ├── 02_clustering/                        # STEP 2: Run experiments
    │   ├── 02_run_simulation.py              # Main simulation (fixed file pairing)
    │   ├── clustering_utils.py               # Helper functions (detection + quality)
    │   ├── clustering_results_final.csv      # 📊 Output (1800 rows)
    │   └── __pycache__/                      # Python cache (ignore)
    │
    ├── 03_analysis/                          # STEP 3: Visualize results
    │   ├── 03_plot_results.ipynb             # Jupyter notebook (3 plots)
    │   ├── plot1_success_rate_comparison.png # 📈 Heatmap: Success rates
    │   ├── plot2_performance_delta.png       # 📈 Heatmap: HC vs K-Means
    │   └── plot3_quality_distortion_boxplots.png # 📈 Boxplots: Quality delta
    │
    └── data/                                 # Generated datasets
        ├── real/                             # 180 real CSV files
        │   ├── N250_p2_k2_rho0_sep0.1_rep1.csv
        │   ├── N250_p2_k2_rho0_sep0.1_rep2.csv
        │   └── ... (178 more)
        │
        └── synthetic/                        # 1800 synthetic CSV files
            ├── N250_p2_k2_rho0_sep0.1_rep1_syn1.csv
            ├── N250_p2_k2_rho0_sep0.1_rep1_syn2.csv
            ├── ... (10 per real dataset)
            └── N250_p10_k4_rho0_sep10_rep5_syn10.csv
```

---

## File Descriptions

### Configuration & Execution

| File | Purpose | Who Edits |
|------|---------|-----------|
| `config.json` | Central parameters (N, p, k, separation, etc.) | Researcher |
| `run_all.sh` | Automated pipeline execution | Run only |
| `requirements.txt` | Python package dependencies | Admin |

### Documentation

| File | Purpose | Audience |
|------|---------|----------|
| `README.md` | Complete pipeline documentation | All users |
| `QUICKSTART.md` | Quick reference commands | New users |
| `DIAGNOSIS_REPORT.md` | Bug fixes & debugging history | Developers |
| `Setup.md` | Original setup notes | Legacy |

### Step 1: Data Generation (R)

| File | Purpose | Input | Output |
|------|---------|-------|--------|
| `01_generate_synthetic.R` | Generate multivariate clusters + synthpop | `config.json` | 180 real + 1800 synthetic CSV |

**Key Fix:** Cluster size distribution handles N÷k remainders correctly

### Step 2: Clustering Simulation (Python)

| File | Purpose | Input | Output |
|------|---------|-------|--------|
| `02_run_simulation.py` | Main simulation loop | CSV files | `clustering_results_final.csv` |
| `clustering_utils.py` | Detection & quality functions | - | Called by main |

**Key Fix:** Processes all 1800 synthetic files (not just 180)

### Step 3: Analysis & Visualization (Jupyter)

| File | Purpose | Input | Output |
|------|---------|-------|--------|
| `03_plot_results.ipynb` | Generate 3 plots | `clustering_results_final.csv` | 3 PNG files |

**Deliverables:**
1. `plot1_success_rate_comparison.png` - Detection success heatmaps
2. `plot2_performance_delta.png` - HC vs K-Means comparison
3. `plot3_quality_distortion_boxplots.png` - Quality degradation analysis

### Data Directory

| Directory | Contents | Count |
|-----------|----------|-------|
| `data/real/` | Real datasets (multivariate normal clusters) | 180 |
| `data/synthetic/` | Synthetic datasets (generated via CART) | 1,800 |

**Naming Convention:**
- Real: `N{N}_p{p}_k{k}_rho{rho}_sep{sep}_rep{rep}.csv`
- Synthetic: `N{N}_p{p}_k{k}_rho{rho}_sep{sep}_rep{rep}_syn{idx}.csv`

---

## Data Flow Diagram

```
┌─────────────────────────────────────────────────────────────────┐
│                         config.json                              │
│  { N:250, p:[2,5,10], k:[2,3,4], sep:[0.1,2,6,10], ... }       │
└──────────────────┬──────────────────────────────────────────────┘
                   │
                   ▼
         ┌─────────────────────┐
         │   Step 1: R Script  │
         │ 01_generate_synthetic.R │
         └──────────┬──────────┘
                    │
      ┌─────────────┴─────────────┐
      ▼                           ▼
  data/real/                 data/synthetic/
  (180 files)                (1800 files)
      │                           │
      └─────────────┬─────────────┘
                    │
                    ▼
         ┌─────────────────────┐
         │  Step 2: Python     │
         │ 02_run_simulation.py│
         └──────────┬──────────┘
                    │
                    ▼
          clustering_results_final.csv
                (1800 rows)
                    │
                    ▼
         ┌─────────────────────┐
         │  Step 3: Jupyter    │
         │ 03_plot_results.ipynb│
         └──────────┬──────────┘
                    │
      ┌─────────────┼─────────────┐
      ▼             ▼             ▼
   plot1.png    plot2.png     plot3.png
```

---

## Key Metrics & Variables

### Configuration Parameters
- **N**: Sample size per dataset (default: 250)
- **p**: Number of features/dimensions (default: [2, 5, 10])
- **k**: True number of clusters (default: [2, 3, 4])
- **separation**: Distance between cluster centers (default: [0.1, 2, 6, 10])
- **rho**: Within-cluster correlation (default: [0.0])
- **n**: Real datasets per combination (default: 5)
- **m**: Synthetic datasets per real (default: 10)

### Output Metrics (CSV columns)
- **success_kmeans**: Binary (1 = detected correct k, 0 = failed)
- **success_hc**: Binary (1 = detected correct k, 0 = failed)
- **diff_sil_km**: Silhouette score delta (Synthetic - Real) for K-Means
- **diff_sil_hc**: Silhouette score delta (Synthetic - Real) for HC
- **diff_dist_km**: Inter-cluster distance delta for K-Means
- **diff_dist_hc**: Inter-cluster distance delta for HC
- **sil_real_km**: Raw Silhouette score on real data (K-Means)
- **sil_syn_km**: Raw Silhouette score on synthetic data (K-Means)

---

## Computational Requirements

### Time Complexity
- **Step 1 (R):** O(n × combinations) ≈ 2-5 minutes
- **Step 2 (Python):** O(m × n × combinations) ≈ 10-15 minutes
- **Step 3 (Jupyter):** O(1) ≈ 1-2 minutes
- **Total:** ~15-20 minutes for default config

### Disk Space
- **Real data:** 180 files × ~50 KB = ~9 MB
- **Synthetic data:** 1800 files × ~50 KB = ~90 MB
- **Results CSV:** ~500 KB
- **Plots:** ~3 MB
- **Total:** ~100 MB

### Memory
- **R script:** ~500 MB RAM
- **Python script:** ~2 GB RAM (parallel processing)
- **Jupyter:** ~1 GB RAM

---

## Modification Guide

### To Change Simulation Parameters
Edit `config.json`:
```json
{
  "simulation": {
    "N": 500,              // Increase sample size
    "n": 10,               // More real datasets
    "m": 20                // More synthetic per real
  },
  "parameters": {
    "p": [5, 10, 20],      // Different dimensions
    "k": [3, 4, 5],        // Different k values
    "separation": [1, 5, 10], // Fewer separation levels
    "rho": [0.0, 0.5]      // Add correlation
  }
}
```

### To Add New Clustering Algorithm
1. Edit `clustering_utils.py`: Add `fit_newmethod()` function
2. Edit `02_run_simulation.py`: Add new columns to output
3. Edit `03_plot_results.ipynb`: Add new visualization

### To Change Synthesis Method
Edit `01_generate_synthetic.R` line 142:
```r
syn_res <- syn(real_data, method = "cart", ...)  # Try "ctree", "sample", etc.
```

---

## Version History

| Version | Date | Changes |
|---------|------|---------|
| 1.0 | Dec 2025 | Initial pipeline implementation |
| 1.1 | Dec 2025 | Fixed cluster rounding bug (R) |
| 1.2 | Dec 2025 | Fixed file pairing bug (Python) |
| 2.0 | Dec 2025 | Complete refactoring with docs |

---

**Last Updated:** December 9, 2025
