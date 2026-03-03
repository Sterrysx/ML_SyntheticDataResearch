# ==============================================================================
# SCRIPT: 01_generate_original_data.R
# PURPOSE: Generate Original Data (OD) for a multiple linear regression
#          simulation using parameters from config.json.
# ==============================================================================
#
# Supports TWO variable types (controlled by config$simulation$var_type):
#
#   "continuous" — X ~ MVN(0, Sigma) via mvtnorm::rmvnorm.
#   "binary"     — X ~ correlated Bernoulli via MultiDiscreteRNG.
#                  Marginal P(X_j = 1) = p1 (from config); correlation = rho.
#
# Response (always continuous):
#   y = cbind(1, X) %*% beta[1:(p+1)] + N(0, sigma_2)
#
# Parameter grid: var_type x N x p x rho x sigma_2 x (p1 if binary) x 1:M
#
# Parallelisation: parallel::mclapply (fork-based) on all-but-2 cores.
#                  Tasks are processed in batches so live progress is shown.
#                  All parent-env variables (bin_cache, beta_full, etc.) are
#                  inherited automatically by forked workers — no clusterExport.
#
# Output filenames:
#   OD_{var_type}_N{N}_p{p}_rho{rho}_sig{sigma_2}_p1{p1}_iter{padded_m}.csv
#   (continuous files use p1=NA for completeness)
# ==============================================================================

# 1. Setup & Imports
# ------------------------------------------------------------------------------
ensure_package <- function(pkg) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    message(paste("Installing package:", pkg))
    # Use only hard dependencies (Depends/Imports/LinkingTo), NOT Suggests.
    # This avoids pulling in system-level extras like MPI, Redis, etc.
    install.packages(pkg, repos = "http://cran.us.r-project.org",
                     dependencies = c("Depends", "Imports", "LinkingTo"))
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      stop(paste("Failed to install package:", pkg))
    }
  }
}

ensure_package("jsonlite")
ensure_package("mvtnorm")
ensure_package("MultiDiscreteRNG")

library(jsonlite)
library(mvtnorm)
library(MultiDiscreteRNG)
library(parallel)

# 2. Load Configuration
# ------------------------------------------------------------------------------
config_path <- "../config/config.json"
if (!file.exists(config_path)) {
  stop("config.json not found! Please ensure it is at ../config/config.json")
}

config <- fromJSON(config_path)

N_vals       <- config$simulation$N
p_vals       <- config$simulation$p
var_types    <- config$simulation$var_type
base_seed    <- config$simulation$random_seed_base
M            <- config$simulation$M

rho_vals     <- config$parameters$rho
sigma_2_vals <- config$parameters$sigma_2
p1_vals      <- config$parameters$p1
beta_full    <- config$parameters$beta   # length >= max(p) + 1

cat("[INFO] Configuration loaded.\n")
cat(sprintf("  N:         %s\n", paste(N_vals, collapse = ", ")))
cat(sprintf("  p:         %s\n", paste(p_vals, collapse = ", ")))
cat(sprintf("  var_type:  %s\n", paste(var_types, collapse = ", ")))
cat(sprintf("  rho:       %s\n", paste(rho_vals, collapse = ", ")))
cat(sprintf("  sigma_2:   %s\n", paste(sigma_2_vals, collapse = ", ")))
cat(sprintf("  p1:        %s\n", paste(p1_vals, collapse = ", ")))
cat(sprintf("  M:         %d\n", M))
cat(sprintf("  beta:      [%s]  (will be subset to p+1)\n\n",
            paste(beta_full, collapse = ", ")))

# 3. Build full parameter grid
# ------------------------------------------------------------------------------
# For continuous: p1 is irrelevant, so we use a sentinel value of NA.
# For binary:     we expand over all p1 values.

grid_parts <- list()

if ("continuous" %in% var_types) {
  grid_parts[["cont"]] <- expand.grid(
    var_type = "continuous",
    N        = N_vals,
    p        = p_vals,
    rho      = rho_vals,
    sigma_2  = sigma_2_vals,
    p1       = NA_real_,
    m        = seq_len(M),
    stringsAsFactors = FALSE
  )
}

if ("binary" %in% var_types) {
  grid_parts[["bin"]] <- expand.grid(
    var_type = "binary",
    N        = N_vals,
    p        = p_vals,
    rho      = rho_vals,
    sigma_2  = sigma_2_vals,
    p1       = p1_vals,
    m        = seq_len(M),
    stringsAsFactors = FALSE
  )
}

grid <- do.call(rbind, grid_parts)
rownames(grid) <- NULL
n_grid <- nrow(grid)

cat(sprintf("[INFO] Parameter grid: %d rows  (%d scenarios x %d iterations).\n",
            n_grid, n_grid %/% M, M))

# 4. Prepare output directory
# ------------------------------------------------------------------------------
output_dir <- file.path("..", "data", "original")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Padding width for iteration numbers (e.g., 1000 -> 4 digits)
pad_width <- nchar(as.character(M))

# 5. Parallel setup
# ------------------------------------------------------------------------------
n_cores <- max(1L, detectCores() - 2L)
cat(sprintf("[INFO] Using %d cores (mclapply, fork-based).\n\n", n_cores))

# 6. Helper: build equi-correlation matrix
# ------------------------------------------------------------------------------
make_sigma_mat <- function(p, rho) {
  S <- matrix(rho, nrow = p, ncol = p)
  diag(S) <- 1
  S
}

# 7. Pre-compute binary correlation objects (expensive; do once per combo)
# ------------------------------------------------------------------------------
# simBinaryCorr.B is O(100k rows) and only depends on (p, rho, p1).
# We cache one object per unique (p, rho, p1) triple.

bin_combos <- unique(grid[grid$var_type == "binary", c("p", "rho", "p1")])
bin_cache  <- list()

if (nrow(bin_combos) > 0L) {
  cat(sprintf("[INFO] Pre-computing %d binary correlation objects...\n",
              nrow(bin_combos)))
  for (i in seq_len(nrow(bin_combos))) {
    pp   <- bin_combos$p[i]
    rrho <- bin_combos$rho[i]
    pp1  <- bin_combos$p1[i]
    key  <- paste(pp, rrho, pp1, sep = "_")

    Sigma <- make_sigma_mat(pp, rrho)

    # simBinaryCorr.B: calibrate internal Gaussian cutpoint matrix.
    # The package prints to BOTH stdout (cat) and stderr (message).
    # Suppress all of it:
    invisible(capture.output(
      suppressMessages(
        bin_obj <- simBinaryCorr.B(
          B.n.vec    = rep(1, pp),
          B.prob.vec = rep(pp1, pp),
          CorrMat    = Sigma,
          no.rows    = 100000
        )
      ),
      type = "output"
    ))
    bin_cache[[key]] <- bin_obj
    cat(sprintf("  [cache] p=%d  rho=%.1f  p1=%.2f  -> OK\n", pp, rrho, pp1))
  }
  cat("\n")
}

# (no clusterExport needed — mclapply forks the parent process, so bin_cache
#  and make_sigma_mat are available in workers automatically)

# 8. Pre-generate reproducible per-task seeds from base_seed
# ------------------------------------------------------------------------------
set.seed(base_seed)
grid$task_seed <- sample.int(.Machine$integer.max, n_grid)

# 9. Single-task worker (called inside mclapply)
# ------------------------------------------------------------------------------
run_task <- function(row) {
  tryCatch({
    set.seed(row$task_seed)

    vtype    <- row$var_type
    N_cur    <- as.integer(row$N)
    p_cur    <- as.integer(row$p)
    rho_cur  <- row$rho
    sig2_cur <- row$sigma_2
    p1_cur   <- row$p1

    current_beta <- beta_full[1:(p_cur + 1L)]
    Sigma        <- make_sigma_mat(p_cur, rho_cur)

    if (vtype == "continuous") {
      X <- rmvnorm(n = N_cur, mean = rep(0, p_cur), sigma = Sigma)
    } else {
      key     <- paste(p_cur, rho_cur, p1_cur, sep = "_")
      bin_obj <- bin_cache[[key]]
    # genB() may return a data.frame with factor columns, a matrix, or mixed.
    # Flatten → character (handles factors) → numeric → reshape to N × p.
    X_raw   <- genB(no.rows = N_cur, binObj = bin_obj)
    X       <- matrix(as.numeric(as.character(unlist(X_raw))),
                      nrow = N_cur, ncol = p_cur)
    }

    X_design <- cbind(1, X)
    epsilon  <- rnorm(N_cur, mean = 0, sd = sqrt(sig2_cur))
    y        <- as.numeric(X_design %*% current_beta + epsilon)

    df           <- as.data.frame(X)
    colnames(df) <- paste0("X", seq_len(p_cur))
    df$y         <- y

    p1_str   <- if (is.na(p1_cur)) "NA" else sprintf("%.2f", p1_cur)
    m_padded <- formatC(as.integer(row$m), width = pad_width, flag = "0")

    fname <- sprintf("OD_%s_N%d_p%d_rho%.1f_sig%.1f_p1%s_iter%s.csv",
                     vtype, N_cur, p_cur, rho_cur, sig2_cur, p1_str, m_padded)
    write.csv(df, file.path(output_dir, fname), row.names = FALSE)
    fname
  }, error = function(e) {
    # Return the error message as a tagged string so the caller can detect it.
    paste0("ERROR:", conditionMessage(e))
  })
}

# 10. Batched parallel loop with live progress
# ------------------------------------------------------------------------------
# Batch size: enough tasks to keep all cores busy between progress updates.
# n_cores * 50 = ~1100 tasks per print line at 22 cores.
BATCH_SIZE <- n_cores * 50L
n_batches  <- ceiling(n_grid / BATCH_SIZE)

cat(sprintf("[INFO] Processing %d tasks in %d batches (~%d tasks/batch, %d cores)...\n",
            n_grid, n_batches, BATCH_SIZE, n_cores))

start_ts  <- proc.time()[[3]]
done      <- 0L
all_names <- character(n_grid)

for (b in seq_len(n_batches)) {
  i_start <- (b - 1L) * BATCH_SIZE + 1L
  i_end   <- min(b * BATCH_SIZE, n_grid)

  batch_rows <- lapply(seq(i_start, i_end), function(i) grid[i, ])

  batch_res <- mclapply(batch_rows, run_task,
                        mc.cores      = n_cores,
                        mc.preschedule = TRUE)

  n_errors <- 0L
  for (j in seq_along(batch_res)) {
    r <- batch_res[[j]]
    if (inherits(r, "try-error")) {
      # mclapply-level failure (fork crash, signal, etc.)
      n_errors <- n_errors + 1L
      if (n_errors <= 5L) {
        cat(sprintf("\n[WARN] Task %d (mclapply): %s",
                    i_start + j - 1L, as.character(r)))
      }
    } else if (is.character(r) && startsWith(r, "ERROR:")) {
      # Our tryCatch-level failure (R error inside the worker).
      n_errors <- n_errors + 1L
      if (n_errors <= 5L) {
        cat(sprintf("\n[WARN] Task %d: %s",
                    i_start + j - 1L, sub("^ERROR:", "", r)))
      }
    } else if (is.character(r) && length(r) == 1L) {
      all_names[i_start + j - 1L] <- r
    } else {
      n_errors <- n_errors + 1L
      if (n_errors <= 5L) {
        cat(sprintf("\n[WARN] Task %d: unexpected result type %s",
                    i_start + j - 1L, class(r)[1]))
      }
    }
  }
  if (n_errors > 5L) {
    cat(sprintf("\n[WARN] ... and %d more errors in this batch.", n_errors - 5L))
  }

  done    <- i_end
  elapsed <- as.integer(proc.time()[[3]] - start_ts)
  pct     <- round(done / n_grid * 100)
  rate    <- if (elapsed > 0L) round(done / elapsed) else 0L
  eta     <- if (rate > 0L) as.integer((n_grid - done) / rate) else NA_integer_
  eta_str <- if (!is.na(eta)) sprintf(" | ETA: %ds", eta) else ""

  cat(sprintf("\r\033[2K[INFO] %d/%d files (%d%%) | %ds elapsed | ~%d files/s%s",
              done, n_grid, pct, elapsed, rate, eta_str))
}
cat("\n")

elapsed_total <- proc.time()[[3]] - start_ts
n_written     <- sum(nchar(all_names) > 0L)
cat(sprintf("\n[DONE] %d OD files written in %.0f s.\n", n_written, elapsed_total))
