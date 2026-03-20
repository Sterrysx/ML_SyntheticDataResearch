# ==============================================================================
# SCRIPT: 01_generate_original_data_linear.R
# PURPOSE: Generate Original Data (OD) for a multiple linear regression
#          simulation using parameters from config.json.
#          ** FULL DYNAMIC QUEUE (WORK STEALING) APPLIED **
# ==============================================================================

# 1. Setup & Imports
# ------------------------------------------------------------------------------
suppressWarnings(suppressPackageStartupMessages({
  ensure_package <- function(pkg) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      message(paste("Installing package:", pkg))
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
  ensure_package("arrow")

  library(jsonlite)
  library(mvtnorm)
  library(MultiDiscreteRNG)
  library(parallel)
  library(arrow)
}))

# 2. Load Configuration & Command Line Overrides
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
sigma_2_vals <- config$parameters$linear$sigma_2
p1_vals      <- config$parameters$p1
beta_full    <- config$parameters$beta   

args <- commandArgs(trailingOnly = TRUE)
run_mode <- if (length(args) > 0) args[1] else "all"

if (run_mode == "continuous") {
  var_types <- "continuous"
  cat("[INFO] Bash override: Running ONLY continuous variables.\n")
} else if (run_mode == "binary") {
  var_types <- "binary"
  cat("[INFO] Bash override: Running ONLY binary variables.\n")
}

# 3. Build unique scenario grid (No iteration explosion!)
# ------------------------------------------------------------------------------
grid_parts <- list()

if ("continuous" %in% var_types) {
  grid_parts[["cont"]] <- expand.grid(
    var_type = "continuous",
    N        = N_vals,
    p        = p_vals,
    rho      = rho_vals,
    sigma_2  = sigma_2_vals,
    p1       = NA_real_,
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
    stringsAsFactors = FALSE
  )
}

grid <- do.call(rbind, grid_parts)

# SHUFFLE THE GRID: This ensures heavy binary tasks are mixed evenly
set.seed(base_seed)
grid <- grid[sample(nrow(grid)), ]
rownames(grid) <- NULL
n_scenarios <- nrow(grid)

cat(sprintf("[INFO] Parameter grid: %d unique scenarios (x %d iterations per core).\n",
            n_scenarios, M))

# 4. Prepare directories
# ------------------------------------------------------------------------------
output_dir <- file.path("..", "data", "original")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 5. Pre-compute binary correlation objects and save to disk for workers
# ------------------------------------------------------------------------------
NUM_CORES <- max(1L, detectCores() - 6L)
bin_combos <- unique(grid[grid$var_type == "binary", c("p", "rho", "p1")])
bin_cache  <- list()

if (nrow(bin_combos) > 0L) {
  cat(sprintf("[INFO] Pre-computing %d binary correlation objects...\n", nrow(bin_combos)))
  
  make_sigma_mat <- function(p, rho) {
    S <- matrix(rho, nrow = p, ncol = p)
    diag(S) <- 1; S
  }

  for (i in seq_len(nrow(bin_combos))) {
    pp   <- bin_combos$p[i]
    rrho <- bin_combos$rho[i]
    pp1  <- bin_combos$p1[i]
    key  <- paste(pp, rrho, pp1, sep = "_")
    Sigma <- make_sigma_mat(pp, rrho)

    invisible(capture.output(
      suppressMessages(
        bin_obj <- simBinaryCorr.B(
          B.n.vec = rep(1, pp), B.prob.vec = rep(pp1, pp),
          CorrMat = Sigma, no.rows = 100000
        )
      ), type = "output"
    ))
    bin_cache[[key]] <- bin_obj
    cat(sprintf("  [cache] %s -> OK\n", key))
  }
}

# Save cache so background workers can load it natively
cache_file <- file.path(tempdir(), "od_linear_bin_cache.rds")
saveRDS(bin_cache, cache_file)
cat("\n")

# 6. Build Atomic Task Queue
# ------------------------------------------------------------------------------
todo_dir  <- file.path(tempdir(), "od_linear_tasks_todo")
doing_dir <- file.path(tempdir(), "od_linear_tasks_doing")
unlink(todo_dir, recursive = TRUE); unlink(doing_dir, recursive = TRUE)
dir.create(todo_dir, showWarnings = FALSE)
dir.create(doing_dir, showWarnings = FALSE)

# Seed generation per scenario guarantees perfect reproducibility regardless of queue order
set.seed(base_seed)
grid$scenario_seed <- sample.int(.Machine$integer.max, n_scenarios)

for (i in seq_len(n_scenarios)) {
  row <- grid[i, ]
  p1_str <- if (is.na(row$p1)) "NA" else sprintf("%.2f", row$p1)
  scenario_key <- sprintf("OD_linear_%s_N%s_p%s_rho%.1f_sig%.1f_p1%s",
                          row$var_type, row$N, row$p, row$rho, row$sigma_2, p1_str)
  
  task_data <- list(
    scenario_key = scenario_key,
    vtype        = row$var_type,
    N_cur        = as.integer(row$N),
    p_cur        = as.integer(row$p),
    rho_cur      = row$rho,
    sig2_cur     = row$sigma_2,
    p1_cur       = row$p1,
    seed         = row$scenario_seed,
    M            = M,
    beta_full    = beta_full
  )
  saveRDS(task_data, file.path(todo_dir, sprintf("task_%05d.rds", i)))
}

cat(sprintf("[INFO] Dynamic Queue built. %d scenario tasks ready for %d cores.\n\n", n_scenarios, NUM_CORES))

# 7. Generate Independent Worker Script
# ------------------------------------------------------------------------------
worker_script <- file.path(tempdir(), "od_linear_worker.R")
worker_code <- c(
  "args <- commandArgs(trailingOnly = TRUE)",
  "progress_file <- args[1]",
  "result_file <- args[2]",
  "todo_dir <- args[3]",
  "doing_dir <- args[4]",
  "cache_file <- args[5]",
  "",
  "Sys.setenv(OMP_NUM_THREADS = '1')",
  "Sys.setenv(OPENBLAS_NUM_THREADS = '1')",
  "Sys.setenv(MKL_NUM_THREADS = '1')",
  "",
  "suppressWarnings(suppressPackageStartupMessages({",
  "  library(mvtnorm)",
  "  library(MultiDiscreteRNG)",
  "  library(arrow)",
  "}))",
  "arrow::set_cpu_count(1)",
  "",
  "bin_cache <- readRDS(cache_file)",
  "output_dir <- file.path('..', 'data', 'original')",
  "results <- list()",
  "",
  "writeLines('READY', progress_file)",
  "",
  "make_sigma_mat <- function(p, rho) { S <- matrix(rho, nrow = p, ncol = p); diag(S) <- 1; S }",
  "",
  "# WORK STEALING LOOP",
  "repeat {",
  "  available_tasks <- list.files(todo_dir, full.names = TRUE)",
  "  if (length(available_tasks) == 0L) break",
  "",
  "  target_task <- available_tasks[1]",
  "  claimed_task <- file.path(doing_dir, basename(target_task))",
  "",
  "  if (file.rename(target_task, claimed_task)) {",
  "    task <- readRDS(claimed_task)",
  "",
  "    # SKIP if output already exists (resume support)",
  "    out_path <- file.path(output_dir, paste0(task$scenario_key, '.parquet'))",
  "    if (file.exists(out_path)) {",
  "      results[[length(results) + 1L]] <- list(status = 'ok', file = task$scenario_key)",
  "      file.remove(claimed_task)",
  "      next",
  "    }",
  "",
  "    writeLines(paste('OD |', sub('^OD_linear_', '', task$scenario_key)), progress_file)",
  "",
  "    tryCatch({",
  "      set.seed(task$seed)",
  "      current_beta <- task$beta_full[1:(task$p_cur + 1L)]",
  "      Sigma <- make_sigma_mat(task$p_cur, task$rho_cur)",
  "      ",
  "      iter_list <- vector('list', task$M)",
  "      ",
  "      # Run all M iterations for this scenario locally",
  "      for (m in seq_len(task$M)) {",
  "        if (task$vtype == 'continuous') {",
  "          X <- rmvnorm(n = task$N_cur, mean = rep(0, task$p_cur), sigma = Sigma)",
  "        } else {",
  "          bin_obj <- bin_cache[[paste(task$p_cur, task$rho_cur, task$p1_cur, sep = '_')]]",
  "          X_raw <- genB(no.rows = task$N_cur, binObj = bin_obj)",
  "          ",
  "          # FIX: Extract just the data matrix, ignore the correlation matrix at X_raw[[2]]",
  "          X_data <- if (is.list(X_raw) && !is.data.frame(X_raw)) X_raw[[1]] else X_raw",
  "          ",
  "          if (is.list(X_data) || is.data.frame(X_data)) {",
  "            X <- matrix(as.integer(unlist(X_data, use.names = FALSE)) - 1L, nrow = task$N_cur, ncol = task$p_cur)",
  "          } else {",
  "            X <- matrix(as.numeric(X_data), nrow = task$N_cur, ncol = task$p_cur)",
  "          }",
  "        }",
  "        ",
  "        X_design <- cbind(1, X)",
  "        beta_raw  <- current_beta[-1]",
  "        lp_var    <- as.numeric(t(beta_raw) %*% Sigma %*% beta_raw)",
  "        beta_norm <- c(current_beta[1], beta_raw / sqrt(lp_var))",
  "        ",
  "        epsilon <- rnorm(task$N_cur, mean = 0, sd = sqrt(task$sig2_cur))",
  "        y       <- as.numeric(X_design %*% beta_norm + epsilon)",
  "        ",
  "        df           <- as.data.frame(X)",
  "        colnames(df) <- paste0('X', seq_len(task$p_cur))",
  "        df$y         <- y",
  "        df$iter      <- as.integer(m)",
  "        ",
  "        iter_list[[m]] <- df",
  "      }",
  "      ",
  "      combined <- do.call(rbind, iter_list)",
  "      ",
  "      write_parquet(",
  "          combined, ",
  "          file.path(output_dir, paste0(task$scenario_key, '.parquet')),",
  "          compression = 'zstd', ",
  "          chunk_size = task$N_cur",
  "      )",
  "      ",
  "      results[[length(results) + 1L]] <- list(status = 'ok', file = task$scenario_key)",
  "      file.remove(claimed_task)",
  "    }, error = function(e) {",
  "      results[[length(results) + 1L]] <<- list(status = 'error', message = e$message, file = task$scenario_key)",
  "    })",
  "  }",
  "}",
  "",
  "writeLines('DONE', progress_file)",
  "saveRDS(results, result_file)"
)
writeLines(worker_code, worker_script)

# 8. Launch completely independent background jobs
# ------------------------------------------------------------------------------
progress_files <- character(NUM_CORES)
result_files   <- character(NUM_CORES)

for (i in seq_len(NUM_CORES)) {
  progress_files[i] <- file.path(tempdir(), sprintf("od_linear_progress_%02d.txt", i))
  result_files[i]   <- file.path(tempdir(), sprintf("od_linear_result_%02d.rds", i))
  
  writeLines("STARTING", progress_files[i])
  if (file.exists(result_files[i])) file.remove(result_files[i])
  
  system2("Rscript", args = c(worker_script, progress_files[i], result_files[i], todo_dir, doing_dir, cache_file), wait = FALSE)
}

# 9. Live Unified Progress Display
# ------------------------------------------------------------------------------
BAR_WIDTH <- 40L

make_bar <- function(done, total) {
  if (total == 0L) return(strrep("-", BAR_WIDTH))
  filled <- as.integer(round(done / total * BAR_WIDTH))
  paste0(strrep("=", filled), strrep("-", BAR_WIDTH - filled))
}

start_ts <- proc.time()[[3]]

cat(sprintf("Global Progress: [ Waiting... ]\n"))
cat(strrep("-", 80), "\n")
for (i in seq_len(NUM_CORES)) {
  cat(sprintf("   Core %-3d : STARTING\033[K\n", i))
}

repeat {
  Sys.sleep(0.5)
  
  cat(sprintf("\033[%dA", NUM_CORES + 2L))
  
  todo_count  <- length(list.files(todo_dir))
  doing_count <- length(list.files(doing_dir))
  done_count  <- n_scenarios - todo_count - doing_count
  
  pct     <- if (n_scenarios > 0) as.integer(round(done_count / n_scenarios * 100)) else 0
  bar     <- make_bar(done_count, n_scenarios)
  elapsed <- as.integer(proc.time()[[3]] - start_ts)
  
  cat(sprintf("\r\033[2K[INFO] [%s] %3d%% (%d/%d) | %ds elapsed\n", 
              bar, pct, done_count, n_scenarios, elapsed))
  cat(strrep("-", 80), "\n")
  
  n_finish <- 0L
  
  for (i in seq_len(NUM_CORES)) {
    status <- "WAITING"
    
    if (file.exists(progress_files[i])) {
      lines <- readLines(progress_files[i], warn = FALSE)
      if (length(lines) > 0) status <- lines[1]
    }
    
    if (file.exists(result_files[i])) {
      status <- "DONE"
      n_finish <- n_finish + 1L
    }
    
    if (nchar(status) > 65) status <- paste0(substr(status, 1, 62), "...")
    
    cat(sprintf("\r\033[2K   Core %-3d : %s\n", i, status))
  }
  
  if (n_finish == NUM_CORES) break
}

# 10. Collect results 
# ------------------------------------------------------------------------------
count  <- 0L
errors <- 0L
warns  <- character(0)

for (i in seq_len(NUM_CORES)) {
  res_file <- result_files[i]
  
  if (!file.exists(res_file)) {
    errors <- errors + 1L
    warns <- c(warns, sprintf("  [FATAL] Core %d crashed entirely without saving output.", i))
    next
  }
  
  task_results <- readRDS(res_file)
  for (r in task_results) {
    if (r$status == "ok") {
      count <- count + 1L
    } else {
      errors <- errors + 1L
      warns  <- c(warns, sprintf("  [WARN] OD | %s: %s", r$file, r$message))
    }
  }
}

if (length(warns) > 0L) cat("\n", paste(warns, collapse = "\n"), "\n", sep = "")

elapsed_total <- proc.time()[[3]] - start_ts
cat(sprintf("\n[DONE] %d scenario file(s) written, %d errors. Wall time: %.0fs\n",
            count, errors, elapsed_total))