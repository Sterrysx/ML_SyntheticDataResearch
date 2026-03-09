import numpy as np
import pandas as pd
from sklearn.cluster import KMeans, AgglomerativeClustering
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import StandardScaler
from scipy.spatial.distance import cdist

# --- A. CLUSTERING ALGORITHMS ---

def fit_kmeans(data, k, seed=123):
    """Fits K-Means and returns labels + centroids."""
    model = KMeans(n_clusters=k, n_init=25, random_state=seed)
    labels = model.fit_predict(data)
    return labels, model.cluster_centers_

def fit_hc(data, k):
    """Fits Hierarchical (Ward) and returns labels + computed centroids."""
    model = AgglomerativeClustering(n_clusters=k, linkage='ward')
    labels = model.fit_predict(data)
    # HC doesn't give centroids, so we calculate them manually
    centroids = np.array([data[labels == i].mean(axis=0) for i in range(k)])
    return labels, centroids

# --- B. LOGIC TEST 1: DETECTION (For Heatmap) ---

def detect_optimal_k(data, method='kmeans', max_k=10):
    """
    Asks the algorithm to GUESS k (from 2 to max_k).
    Returns the k that maximizes Silhouette Score.
    """
    best_score = -1
    best_k = 2
    
    # We create the Scaler once here
    scaler = StandardScaler()
    data_scaled = scaler.fit_transform(data)
    
    for k in range(2, max_k + 1):
        try:
            if method == 'kmeans':
                labels, _ = fit_kmeans(data_scaled, k)
            else:
                labels, _ = fit_hc(data_scaled, k)
            
            if len(np.unique(labels)) < 2: continue
            
            score = silhouette_score(data_scaled, labels)
            if score > best_score:
                best_score = score
                best_k = k
        except:
            continue
            
    return best_k

# --- C. LOGIC TEST 2: QUALITY (For Boxplots) ---

def calculate_quality_metrics(data, k_fixed, method='kmeans'):
    """
    Forces the algorithm to use 'k_fixed'.
    Returns Silhouette Score and Mean Inter-Cluster Distance.
    """
    scaler = StandardScaler()
    data_scaled = scaler.fit_transform(data)
    
    try:
        if method == 'kmeans':
            labels, centroids = fit_kmeans(data_scaled, k_fixed)
        else:
            labels, centroids = fit_hc(data_scaled, k_fixed)
            
        # 1. Silhouette
        sil = silhouette_score(data_scaled, labels)
        
        # 2. Distance (Geometry check)
        if len(centroids) < 2: dist_mean = 0
        else:
            dists = cdist(centroids, centroids, metric='euclidean')
            dist_mean = dists[np.triu_indices(len(centroids), k=1)].mean()
            
        return sil, dist_mean
    except:
        return np.nan, np.nan