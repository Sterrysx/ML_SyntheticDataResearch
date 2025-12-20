import pandas as pd
import numpy as np
import json
import glob
import os
import re
from joblib import Parallel, delayed
from tqdm import tqdm
from clustering_utils import detect_optimal_k, calculate_quality_metrics

# 1. LOAD CONFIG
with open("../config.json", "r") as f:
    config = json.load(f)

# 2. WORKER FUNCTION
def process_synthetic_dataset(syn_path):
    """
    Processes ONE synthetic CSV file against its corresponding real file.
    Returns a dictionary with BOTH Heatmap and Boxplot metrics.
    """
    try:
        # Parse Synthetic Filename
        # Example: "N250_p10_k2_rho0_sep0.1_rep1_syn5_normal.csv" or "..._gamma.csv"
        syn_filename = os.path.basename(syn_path)
        
        # Extract parameters using regex (includes distribution suffix)
        match = re.match(r'N(\d+)_p(\d+)_k(\d+)_rho([\d.]+)_sep([\d.]+)_rep(\d+)_syn(\d+)_(normal|gamma)\.csv', syn_filename)
        if not match:
            print(f"⚠️ Could not parse filename: {syn_filename}")
            return None
        
        N, p, k, rho, sep, rep, syn_idx, distribution = match.groups()
        k_truth = int(k)
        
        # Build corresponding real filename (without _synX suffix but WITH distribution)
        real_filename = f"N{N}_p{p}_k{k}_rho{rho}_sep{sep}_rep{rep}_{distribution}.csv"
        real_path = os.path.join("../data/real", real_filename)
        
        if not os.path.exists(real_path):
            print(f"⚠️ Real file not found: {real_filename}")
            return None
        
        # Load Data
        df_real = pd.read_csv(real_path)
        df_syn = pd.read_csv(syn_path)
        
        # Drop non-numeric columns (like 'group' or 'id')
        X_real = df_real.select_dtypes(include=[np.number]).drop(columns=['group'], errors='ignore')
        X_syn = df_syn.select_dtypes(include=[np.number]).drop(columns=['group'], errors='ignore')

        # --- LOGIC PART 1: DETECTION (Heatmap) ---
        # Did the algorithm find the right k when run on SYNTHETIC data?
        k_found_km = detect_optimal_k(X_syn, method='kmeans')
        k_found_hc = detect_optimal_k(X_syn, method='hierarchical')
        
        success_km = 1 if k_found_km == k_truth else 0
        success_hc = 1 if k_found_hc == k_truth else 0
        
        # --- LOGIC PART 2: QUALITY (Boxplot) ---
        # How distorted is the shape when we force the true k on SYNTHETIC data?
        
        # Real Baseline - calculate quality on REAL data
        sil_real_km, dist_real_km = calculate_quality_metrics(X_real, k_truth, 'kmeans')
        sil_real_hc, dist_real_hc = calculate_quality_metrics(X_real, k_truth, 'hierarchical')
        
        # Synthetic Quality - calculate quality on SYNTHETIC data
        sil_syn_km, dist_syn_km = calculate_quality_metrics(X_syn, k_truth, 'kmeans')
        sil_syn_hc, dist_syn_hc = calculate_quality_metrics(X_syn, k_truth, 'hierarchical')
        
        return {
            'N': int(N),
            'p': int(p),
            'k': k_truth,
            'rho': float(rho),
            'sep': float(sep),
            'rep': int(rep),
            'syn_idx': int(syn_idx),
            'distribution': distribution,
            # Heatmap Data (Detection on Synthetic)
            'success_kmeans': success_km,
            'success_hc': success_hc,
            # Boxplot Data (Quality difference: Synthetic - Real)
            'diff_sil_km': sil_syn_km - sil_real_km,
            'diff_sil_hc': sil_syn_hc - sil_real_hc,
            'diff_dist_km': dist_syn_km - dist_real_km,
            'diff_dist_hc': dist_syn_hc - dist_real_hc,
            # Raw metrics for debugging
            'sil_real_km': sil_real_km,
            'sil_syn_km': sil_syn_km,
            'dist_real_km': dist_real_km,
            'dist_syn_km': dist_syn_km
        }
        
    except Exception as e:
        print(f"❌ Error processing {syn_path}: {e}")
        import traceback
        traceback.print_exc()
        return None

# 3. MAIN EXECUTION
if __name__ == "__main__":
    # Get ALL synthetic files (each needs to be evaluated)
    syn_files = sorted(glob.glob("../data/synthetic/*.csv"))
    
    print(f"🚀 Found {len(syn_files)} synthetic datasets to process...")
    
    # Run in Parallel - process EACH synthetic file
    results = Parallel(n_jobs=-1)(
        delayed(process_synthetic_dataset)(syn_path) 
        for syn_path in tqdm(syn_files, desc="Processing")
    )
    
    # Filter out None results and save
    valid_results = [r for r in results if r is not None]
    print(f"✅ Successfully processed {len(valid_results)} / {len(syn_files)} datasets")
    
    df_final = pd.DataFrame(valid_results)
    df_final.to_csv("clustering_results_final.csv", index=False)
    print(f"💾 Results saved to clustering_results_final.csv ({len(df_final)} rows)")