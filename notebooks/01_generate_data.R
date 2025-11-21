# ==============================================================================
# SCRIPT 1: DATA GENERATION (Run in RStudio)
# ==============================================================================

# 1. Load Libraries
if (!require("mvtnorm")) install.packages("mvtnorm")
if (!require("synthpop")) install.packages("synthpop")
if (!require("dplyr")) install.packages("dplyr")

library(mvtnorm)
library(synthpop)
library(dplyr)

# ==========================================
# 2. Define Core Functions
# ==========================================

centroid.generation <- function(k, p, separation) {
  centroids <- matrix(NA, nrow = k, ncol = p)
  centroids[1, ] <- rep(0, p)
  for (i in 2:k) {
    tries <- 0
    repeat {
      vector <- rnorm(p)
      vector <- vector / sqrt(sum(vector^2)) * separation
      dists <- sqrt(rowSums((centroids[1:(i - 1), , drop = FALSE] - 
                               matrix(vector, nrow = i - 1, ncol = p, byrow = TRUE))^2))
      if (all(dists >= separation)) {
        centroids[i, ] <- vector
        break
      }
      tries <- tries + 1
      if (tries > 1000) stop("Could not generate separated centroids.")
    }
  }
  return(centroids)
}

data.generation <- function(N, p, rho, separation, k){
  centroids <- centroid.generation(k = k, p = p, separation = separation)
  R <- matrix(rho, nrow = p, ncol = p)
  diag(R) <- 1
  sd <- rep(1, p)
  D <- diag(sd)
  Sigma <- D %*% R %*% D
  n <- round(N/k, digits = 0)
  clusters <- lapply(1:k, function(i) {
    rmvnorm(n, mean = centroids[i, ], sigma = Sigma)
  })
  data <- do.call(rbind, clusters)
  data <- data.frame(data, group = factor(rep(1:k, each = n)))
  return(data)
}

# ==========================================
# 3. Main Simulation Loop
# ==========================================

# Output Directory
out_dir <- "simulation_data"
if (!dir.exists(out_dir)) dir.create(out_dir)

# Parameters
N_list <- c(250)
p_list <- c(2, 5, 10)
k_list <- c(2, 3, 4)
rho <- 0.0
sep_list <- c(0.1, 2, 6, 10)
n_simulations <- 50 

# Create Grid
param_grid <- expand.grid(N = N_list, p = p_list, k = k_list, sep = sep_list, 
                          rho = rho, stringsAsFactors = FALSE)

cat(sprintf("Starting generation of %d combinations x %d reps = %d datasets...\n", 
            nrow(param_grid), n_simulations, nrow(param_grid)*n_simulations))

# Loop
counter <- 0
total_ops <- nrow(param_grid) * n_simulations

for(i in 1:nrow(param_grid)) {
  params <- param_grid[i,]
  
  for(seed in 1:n_simulations) {
    set.seed(seed + (i * 1000)) # Unique seed per combo
    
    # 1. Generate Real
    # Using tryCatch to handle impossible geometries (e.g. k=4, p=2, sep=10 might fail to find centroids)
    tryCatch({
      rdata_list <- data.generation(N = params$N, p = params$p, rho = params$rho, 
                                    separation = params$sep, k = params$k)
      real_data <- rdata_list
      
      # 2. Generate Synthetic (CART)
      syn_obj <- syn(real_data, method = "cart", m = 1, print.flag = FALSE)
      syn_data <- syn_obj$syn
      
      # 3. Construct Filename
      # Format: N_p_k_sep_seed_type.csv
      file_base <- sprintf("%s/sim_N%d_p%d_k%d_sep%g_rho%g_seed%d", 
                           out_dir, params$N, params$p, params$k, params$sep, params$rho, seed)
      
      write.csv(real_data, file = paste0(file_base, "_real.csv"), row.names = FALSE)
      write.csv(syn_data, file = paste0(file_base, "_syn.csv"), row.names = FALSE)
      
    }, error = function(e) {
      cat(sprintf("\nSkipping invalid combo: %s\n", e$message))
    })
    
    counter <- counter + 1
    if(counter %% 100 == 0) cat(sprintf("Progress: %d / %d\n", counter, total_ops))
  }
}

cat("Data generation complete.\n")