# ==============================================================================
# SCRIPT: 01_generate_synthetic.R
# PURPOSE: Generate Real and Synthetic data using parameters from config.json
# ==============================================================================

# 1. Setup & Imports
# ------------------------------------------------------------------------------
ensure_package <- function(pkg) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    message(paste("Installing package:", pkg))
    install.packages(pkg, repos = "http://cran.us.r-project.org", dependencies = TRUE, type = "source")
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      stop(paste("Failed to install package:", pkg))
    }
  }
}

# Try to load essential packages
ensure_package("jsonlite")
ensure_package("mvtnorm")
ensure_package("dplyr")

# For synthpop, try binary first, then source
if (!require("synthpop", quietly = TRUE)) {
  message("Installing synthpop package...")
  # Try binary first
  install.packages("synthpop", repos = "http://cran.us.r-project.org", dependencies = TRUE)
  if (!require("synthpop", quietly = TRUE)) {
    stop("Failed to install synthpop. Please install it manually: install.packages('synthpop')")
  }
}

# 2. Load Configuration
# ------------------------------------------------------------------------------
config_path <- "../config/config.json"
if (!file.exists(config_path)) stop("❌ config.json not found! Please ensure it is in ../config/config.json")

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

data.generation <- function(N, p, rho, separation, k, distribution_config) {
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
  
  # Generate data based on distribution type
  # Handle both list and data.frame formats from jsonlite
  if (is.list(distribution_config) && "name" %in% names(distribution_config)) {
    dist_name <- distribution_config$name
  } else if (is.data.frame(distribution_config) && "name" %in% names(distribution_config)) {
    dist_name <- distribution_config$name[1]
  } else {
    stop("Invalid distribution_config format")
  }
  
  if (dist_name == "normal") {
    # Original multivariate normal generation
    clusters_list <- lapply(1:k, function(i) {
      rmvnorm(cluster_sizes[i], mean = centroids[i, ], sigma = R)
    })
    data <- do.call(rbind, clusters_list)
    
  } else if (dist_name == "gamma") {
    # Gamma distribution generation
    # Extract parameters safely
    if (is.list(distribution_config)) {
      shape_param <- distribution_config$shape
      rate_param <- distribution_config$rate
    } else if (is.data.frame(distribution_config)) {
      shape_param <- distribution_config$shape[1]
      rate_param <- distribution_config$rate[1]
    } else {
      stop("Cannot extract gamma parameters")
    }
    
    clusters_list <- lapply(1:k, function(i) {
      # Generate from gamma distribution
      gamma_data <- matrix(rgamma(cluster_sizes[i] * p, shape = shape_param, rate = rate_param),
                          nrow = cluster_sizes[i], ncol = p)
      
      # Apply z-score normalization to each dimension
      gamma_data_normalized <- scale(gamma_data)
      
      # Apply correlation structure
      chol_R <- chol(R)
      correlated_data <- gamma_data_normalized %*% chol_R
      
      # Shift to cluster centroid
      shifted_data <- sweep(correlated_data, 2, centroids[i, ], "+")
      
      return(shifted_data)
    })
    data <- do.call(rbind, clusters_list)
    
  } else {
    stop(paste("Unknown distribution:", dist_name))
  }
  
  df <- as.data.frame(data)
  colnames(df) <- paste0("X", 1:p)
  
  # FIX: Assign groups matching the exact cluster sizes used above
  df$group <- factor(rep(1:k, times = cluster_sizes))
  
  return(df)
}

# 4. Main Simulation Loop
# ------------------------------------------------------------------------------
# Parse distribution configurations
# Note: fromJSON may return a data frame if objects are uniform
dist_configs_raw <- config$parameters$distribution

# Convert to list of lists for consistent access
if (is.data.frame(dist_configs_raw)) {
  # If it's a data frame, convert each row to a list
  dist_configs <- lapply(1:nrow(dist_configs_raw), function(i) {
    as.list(dist_configs_raw[i, , drop = FALSE])
  })
} else if (is.list(dist_configs_raw)) {
  dist_configs <- dist_configs_raw
} else {
  # Fallback to normal if not specified
  dist_configs <- list(list(name = "normal"))
}

# Extract distribution names for display
dist_names <- sapply(dist_configs, function(x) {
  if (is.list(x) && "name" %in% names(x)) {
    return(x$name)
  } else if (is.data.frame(x) && "name" %in% names(x)) {
    return(x$name[1])
  } else {
    return("unknown")
  }
})

param_grid <- expand.grid(
  p = config$parameters$p,
  k = config$parameters$k,
  sep = config$parameters$separation,
  rho = config$parameters$rho,
  dist_idx = 1:length(dist_configs),
  stringsAsFactors = FALSE
)

total_combos <- nrow(param_grid)
cat(sprintf("🚀 Starting Generation: %d Combinations\n", total_combos))
cat(sprintf("   - Sample Size (N): %d\n", N_val))
cat(sprintf("   - Real Datasets (n): %d\n", n_real))
cat(sprintf("   - Synthetic Reps (m): %d\n", m_syn))
cat(sprintf("   - Distributions: %s\n", paste(dist_names, collapse = ", ")))

for (i in 1:total_combos) {
  params <- param_grid[i, ]
  dist_config <- dist_configs[[params$dist_idx]]
  
  # Extract distribution name safely
  if (is.list(dist_config) && "name" %in% names(dist_config)) {
    dist_name <- dist_config$name
  } else if (is.data.frame(dist_config) && "name" %in% names(dist_config)) {
    dist_name <- dist_config$name[1]
  } else {
    dist_name <- "unknown"
  }
  
  for (r in 1:n_real) {
    current_seed <- base_seed + (i * 10000) + r
    set.seed(current_seed)
    
    # A. Generate Real Data
    real_data <- tryCatch({
      data.generation(N = N_val, p = params$p, rho = params$rho, 
                      separation = params$sep, k = params$k,
                      distribution_config = dist_config)
    }, error = function(e) {
      cat(sprintf("⚠️ Error generating real data (p=%d, k=%d, dist=%s): %s\n", 
                  params$p, params$k, dist_name, e$message))
      return(NULL)
    })
    
    if (!is.null(real_data)) {
      real_filename <- sprintf("N%d_p%d_k%d_rho%s_sep%s_rep%d_%s.csv", 
                               N_val, params$p, params$k, params$rho, params$sep, r, dist_name)
      
      write.csv(real_data, file.path("../data/real", real_filename), row.names = FALSE)
      
      # B. Generate Synthetic Datasets
      for (syn_idx in 1:m_syn) {
        syn_seed <- current_seed + (syn_idx * 500)
        
        # Use tryCatch for synthesis as well to be safe
        tryCatch({
          syn_res <- syn(real_data, method = "cart", m = 1, seed = syn_seed, 
                         print.flag = FALSE, proper = TRUE, cart.minbucket = 10)
          
          syn_data <- syn_res$syn
          
          syn_filename <- sprintf("N%d_p%d_k%d_rho%s_sep%s_rep%d_syn%d_%s.csv", 
                                  N_val, params$p, params$k, params$rho, params$sep, r, syn_idx, dist_name)
          
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