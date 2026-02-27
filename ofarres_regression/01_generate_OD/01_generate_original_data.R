# ==============================================================================
# SCRIPT: 01_generate_original_data.R
# PURPOSE: Generate Original Data (OD) for a multiple linear regression
#          simulation using parameters from config.json
# ==============================================================================
#
# Predictors X ~ MVN(0, Sigma) where Sigma has 1s on the diagonal and rho on
# the off-diagonals.  Response y = X %*% beta + epsilon, epsilon ~ N(0, sigma_2).
#
# Now loops over BOTH rho AND sigma_2 from config.json.
#
# Outputs:  ../data/original/OD_N<N>_p<p>_rho<rho>_sigma<sigma>_rep<m>.csv
# ==============================================================================

# 1. Setup & Imports
# ------------------------------------------------------------------------------
ensure_package <- function(pkg) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    message(paste("Installing package:", pkg))
    install.packages(pkg, repos = "http://cran.us.r-project.org",
                     dependencies = TRUE)
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      stop(paste("Failed to install package:", pkg))
    }
  }
}

ensure_package("jsonlite")
ensure_package("mvtnorm")

# 2. Load Configuration
# ------------------------------------------------------------------------------
config_path <- "../config/config.json"
if (!file.exists(config_path)) {
  stop("config.json not found! Please ensure it is at ../config/config.json")
}

config <- fromJSON(config_path)

N            <- config$simulation$N
p            <- config$simulation$p
base_seed    <- config$simulation$random_seed_base
M            <- config$simulation$M
rho_vals     <- config$parameters$rho
sigma_2_vals <- config$parameters$sigma_2
beta         <- config$parameters$beta   # length p + 1 (intercept + p slopes)

cat(sprintf("Configuration loaded: N=%d, p=%d, seed=%d, M=%d\n", N, p, base_seed, M))
cat(sprintf("Rho values:     %s\n", paste(rho_vals, collapse = ", ")))
cat(sprintf("Sigma_2 values: %s\n", paste(sigma_2_vals, collapse = ", ")))
cat(sprintf("Beta vector:    [%s]\n\n", paste(beta, collapse = ", ")))

# 3. Helper Functions
# ------------------------------------------------------------------------------

# Build an equi-correlation covariance matrix (1 on diag, rho off-diag)
make_sigma <- function(p, rho) {
  Sigma <- matrix(rho, nrow = p, ncol = p)
  diag(Sigma) <- 1
  return(Sigma)
}

# Generate one Original Data (OD) dataset
generate_od <- function(N, p, rho, beta, sigma_2, seed) {
  set.seed(seed)

  # Covariance matrix
  Sigma <- make_sigma(p, rho)

  # Predictors X ~ MVN(0, Sigma)
  X <- rmvnorm(n = N, mean = rep(0, p), sigma = Sigma)

  # Design matrix with intercept column
  X_design <- cbind(1, X)   # N x (p+1)

  # Error term
  epsilon <- rnorm(N, mean = 0, sd = sqrt(sigma_2))

  # Response variable
  y <- X_design %*% beta + epsilon

  # Assemble data frame (columns: X1, X2, ..., Xp, y)
  df <- as.data.frame(X)
  colnames(df) <- paste0("X", seq_len(p))
  df$y <- as.numeric(y)

  return(df)
}

# 4. Main Generation Loop (rho x sigma_2 x M repetitions)
# ------------------------------------------------------------------------------
output_dir <- file.path("..", "data", "original")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

total_datasets <- length(rho_vals) * length(sigma_2_vals) * M
cat(sprintf("Generating %d OD dataset(s): %d rho x %d sigma_2 x %d reps...\n\n",
            total_datasets, length(rho_vals), length(sigma_2_vals), M))

count <- 0
for (ri in seq_along(rho_vals)) {
  rho <- rho_vals[ri]

  for (si in seq_along(sigma_2_vals)) {
    sigma_2 <- sigma_2_vals[si]

    for (m in seq_len(M)) {
      # Unique seed for each (rho, sigma_2, rep) combination
      seed <- base_seed + (ri * 10000) + (si * 1000) + m

      od <- generate_od(N = N, p = p, rho = rho, beta = beta,
                        sigma_2 = sigma_2, seed = seed)

      fname <- sprintf("OD_N%d_p%d_rho%.1f_sigma%.1f_rep%d.csv", N, p, rho, sigma_2, m)
      fpath <- file.path(output_dir, fname)
      write.csv(od, fpath, row.names = FALSE)

      count <- count + 1
      if (m == 1 || m == M || m %% 25 == 0) {
        cat(sprintf("[INFO] %s  (rho=%.1f, sigma=%.1f, rep=%d)  [%d/%d]\n",
                    fname, rho, sigma_2, m, count, total_datasets))
      }
    }
  }
}

cat(sprintf("\n[DONE] Original data generation complete. %d files written.\n", count))
