# ==============================================================================
# SCRIPT: 01_generate_synthetic.R
# PURPOSE: Generate Real and Synthetic data using parameters from config.json
# ==============================================================================

# 1. Setup & Imports
# ------------------------------------------------------------------------------
ensure_package <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, repos = "http://cran.us.r-project.org")
    library(pkg, character.only = TRUE)
  }
}

ensure_package("jsonlite")
ensure_package("mvtnorm")
ensure_package("synthpop")
ensure_package("dplyr")

# 2. Load Configuration
# ------------------------------------------------------------------------------
config_path <- "../config.json"
if (!file.exists(config_path)) stop("❌ config.json not found! Please ensure it is in the parent directory.")

# Read JSON
config <- fromJSON(config_path)

N_val <- config$simulation$N
n_real <- config$simulation$n
m_syn  <- config$simulation$m
base_seed <- config$simulation$random_seed_base

# Create Output Directories
dir.create("../data/real", recursive = TRUE, showWarnings = FALSE)
dir.create("../data/synthetic", recursive = TRUE, showWarnings = FALSE)

# 3. Helper Functions
# ------------------------------------------------------------------------------
centroid.generation <- function(k, p, separation) {
  centroids <- matrix(NA, nrow = k, ncol = p)
  centroids[1, ] <- rep(0, p)
  
  if (k > 1) {
    for (i in 2:k) {
      tries <- 0
      repeat {
        vector <- rnorm(p)
        # Normalize and scale
        vector <- vector / sqrt(sum(vector^2)) * separation
        
        # Calculate distances to all existing centroids
        dists <- sqrt(rowSums((centroids[1:(i - 1), , drop = FALSE] - 
                                 matrix(vector, nrow = i - 1, ncol = p, byrow = TRUE))^2))
        
        if (all(dists >= separation)) {
          centroids[i, ] <- vector
          break
        }
        tries <- tries + 1
        if (tries > 1000) stop(paste("Could not generate separated centroids for k=", k, "p=", p))
      }
    }
  }
  return(centroids)
}

data.generation <- function(N, p, rho, separation, k) {
  centroids <- centroid.generation(k, p, separation)
  R <- matrix(rho, nrow = p, ncol = p)
  diag(R) <- 1
  
  # --- FIX: EXACT SIZE CALCULATION ---
  # Distribute N evenly across k clusters, adding remainder to first few
  base_n <- floor(N / k)
  remainder <- N %% k
  cluster_sizes <- rep(base_n, k)
  
  if (remainder > 0) {
    cluster_sizes[1:remainder] <- cluster_sizes[1:remainder] + 1
  }
  # ------------------------------------
  
  clusters_list <- lapply(1:k, function(i) {
    rmvnorm(cluster_sizes[i], mean = centroids[i, ], sigma = R)
  })
  
  data <- do.call(rbind, clusters_list)
  
  df <- as.data.frame(data)
  colnames(df) <- paste0("X", 1:p)
  
  # FIX: Assign groups matching the exact cluster sizes used above
  df$group <- factor(rep(1:k, times = cluster_sizes))
  
  return(df)
}

# 4. Main Simulation Loop
# ------------------------------------------------------------------------------
param_grid <- expand.grid(
  p = config$parameters$p,
  k = config$parameters$k,
  sep = config$parameters$separation,
  rho = config$parameters$rho,
  stringsAsFactors = FALSE
)

total_combos <- nrow(param_grid)
cat(sprintf("🚀 Starting Generation: %d Combinations\n", total_combos))
cat(sprintf("   - Sample Size (N): %d\n", N_val))
cat(sprintf("   - Real Datasets (n): %d\n", n_real))
cat(sprintf("   - Synthetic Reps (m): %d\n", m_syn))

for (i in 1:total_combos) {
  params <- param_grid[i, ]
  
  for (r in 1:n_real) {
    current_seed <- base_seed + (i * 10000) + r
    set.seed(current_seed)
    
    # A. Generate Real Data
    real_data <- tryCatch({
      data.generation(N = N_val, p = params$p, rho = params$rho, 
                      separation = params$sep, k = params$k)
    }, error = function(e) {
      cat(sprintf("⚠️ Error generating real data (p=%d, k=%d): %s\n", params$p, params$k, e$message))
      return(NULL)
    })
    
    if (!is.null(real_data)) {
      real_filename <- sprintf("N%d_p%d_k%d_rho%s_sep%s_rep%d.csv", 
                               N_val, params$p, params$k, params$rho, params$sep, r)
      
      write.csv(real_data, file.path("../data/real", real_filename), row.names = FALSE)
      
      # B. Generate Synthetic Datasets
      for (syn_idx in 1:m_syn) {
        syn_seed <- current_seed + (syn_idx * 500)
        
        # Use tryCatch for synthesis as well to be safe
        tryCatch({
          syn_res <- syn(real_data, method = "cart", m = 1, seed = syn_seed, 
                         print.flag = FALSE, proper = TRUE, cart.minbucket = 10)
          
          syn_data <- syn_res$syn
          
          syn_filename <- sprintf("N%d_p%d_k%d_rho%s_sep%s_rep%d_syn%d.csv", 
                                  N_val, params$p, params$k, params$rho, params$sep, r, syn_idx)
          
          write.csv(syn_data, file.path("../data/synthetic", syn_filename), row.names = FALSE)
        }, error = function(e) {
           cat(sprintf("   ⚠️ Synthesis failed for rep %d: %s\n", syn_idx, e$message))
        })
      }
    }
  }
  
  if (i %% 5 == 0) cat(sprintf("... Processed %d / %d combinations\n", i, total_combos))
}

cat("✅ Data generation complete. Check '../data/real' and '../data/synthetic'.\n")