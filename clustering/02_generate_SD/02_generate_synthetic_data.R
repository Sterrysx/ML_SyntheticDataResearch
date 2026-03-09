# ==============================================================================
# SCRIPT: 02_generate_synthetic_data.R
# PURPOSE: Generate Synthetic Data (SD) from each Original Data (OD) parquet
#          file using synthpop (cart method).
#          ** FULL DYNAMIC QUEUE (WORK STEALING) APPLIED **
# ==============================================================================

# 1. Setup & Imports
# ------------------------------------------------------------------------------
suppressWarnings(suppressPackageStartupMessages({
  ensure_package <- function(pkg) {
    if (!(pkg %in% installed.packages()[,"Package"])) {
      message(paste("Installing package:", pkg))
      install.packages(pkg, repos = "http://cran.us.r-project.org", dependencies = TRUE)
      if (!(pkg %in% installed.packages()[,"Package"])) stop(paste("Failed to install:", pkg))
    }
  }

  ensure_package("jsonlite")
  ensure_package("synthpop")
  ensure_package("arrow")

  library(jsonlite)
}))

# 2. Load Configuration
# ------------------------------------------------------------------------------
config_path <- "../config/config.json"
if (!file.exists(config_path)) stop("config.json not found!")

config     <- fromJSON(config_path)
m_syn      <- config$simulation$m             # number of synthetic reps per OD rep
base_seed  <- config$simulation$random_seed_base

cat(sprintf("[INFO] Synthetic reps per OD rep: %d (method = cart)\n", m_syn))

# 3. Discover OD files
# ------------------------------------------------------------------------------
input_dir  <- file.path("..", "data", "original")
output_dir <- file.path("..", "data", "synthetic")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

od_files <- sort(list.files(input_dir, pattern = "^OD_.*\\.parquet$", full.names = TRUE))
if (length(od_files) == 0) {
  stop("[ERROR] No OD parquet files found in: ", input_dir,
       "\n        Run 01_generate_original_data.R first.")
}

# Each OD file generates m_syn SD files
total_tasks <- length(od_files) * m_syn
cat(sprintf("[INFO] Found %d OD file(s) x %d syn reps = %d SD file(s) to generate.\n",
            length(od_files), m_syn, total_tasks))

# 4. Build Atomic Task Queue
# ------------------------------------------------------------------------------
NUM_CORES <- max(1L, parallel::detectCores() - 6L)

todo_dir  <- file.path(tempdir(), "sd_clust_tasks_todo")
doing_dir <- file.path(tempdir(), "sd_clust_tasks_doing")
unlink(todo_dir, recursive = TRUE); unlink(doing_dir, recursive = TRUE)
dir.create(todo_dir, showWarnings = FALSE)
dir.create(doing_dir, showWarnings = FALSE)

task_i <- 1L
set.seed(999)
shuffled_od <- sample(od_files)

for (od_path in shuffled_od) {
  for (syn_idx in seq_len(m_syn)) {
    task_data <- list(file = od_path, syn_idx = syn_idx)
    saveRDS(task_data, file.path(todo_dir, sprintf("task_%05d.rds", task_i)))
    task_i <- task_i + 1L
  }
}

cat(sprintf("[INFO] Dynamic Queue built. %d tasks ready for %d cores.\n\n",
            total_tasks, NUM_CORES))

# 5. Generate Independent Worker Script
# ------------------------------------------------------------------------------
worker_script <- file.path(tempdir(), "sd_clust_worker.R")
worker_code <- c(
  "args <- commandArgs(trailingOnly = TRUE)",
  "progress_file <- args[1]",
  "result_file <- args[2]",
  "todo_dir <- args[3]",
  "doing_dir <- args[4]",
  "",
  "Sys.setenv(OMP_NUM_THREADS = '1')",
  "Sys.setenv(OPENBLAS_NUM_THREADS = '1')",
  "Sys.setenv(MKL_NUM_THREADS = '1')",
  "",
  "suppressWarnings(suppressPackageStartupMessages({",
  "  library(jsonlite)",
  "  library(synthpop)",
  "  library(arrow)",
  "}))",
  "arrow::set_cpu_count(1)",
  "",
  "config <- fromJSON('../config/config.json')",
  "base_seed <- config$simulation$random_seed_base",
  "output_dir <- file.path('..', 'data', 'synthetic')",
  "",
  "full_od_files <- sort(list.files(file.path('..', 'data', 'original'),",
  "                                  pattern = '^OD_.*[.]parquet$', full.names = TRUE))",
  "",
  "results <- list()",
  "writeLines('READY', progress_file)",
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
  "    od_path  <- task$file",
  "    syn_idx  <- task$syn_idx",
  "    od_name  <- basename(od_path)",
  "",
  "    # Build SD output name: SD_cart_<scenario>_syn<idx>.parquet",
  "    scenario_tag <- sub('^OD_', '', sub('[.]parquet$', '', od_name))",
  "    sd_name <- sprintf('SD_cart_%s_syn%d', scenario_tag, syn_idx)",
  "",
  "    out_path <- file.path(output_dir, paste0(sd_name, '.parquet'))",
  "    if (file.exists(out_path)) {",
  "      results[[length(results) + 1L]] <- list(status = 'ok', file = sd_name)",
  "      file.remove(claimed_task)",
  "      next",
  "    }",
  "",
  "    writeLines(paste('CART |', scenario_tag, 'syn', syn_idx), progress_file)",
  "",
  "    tryCatch({",
  "      od_full <- as.data.frame(read_parquet(od_path))",
  "      global_idx <- match(od_path, full_od_files)",
  "",
  "      x_cols <- grep('^X[0-9]+$', names(od_full), value = TRUE)",
  "      reps <- sort(unique(od_full$rep))",
  "",
  "      sd_list <- list()",
  "      for (r in reps) {",
  "        od_rep <- od_full[od_full$rep == r, ]",
  "        od_rep$rep <- NULL",
  "",
  "        syn_seed <- base_seed + (global_idx * 10000L) + (r * 100L) + syn_idx",
  "",
  "        invisible(capture.output({",
  "          syn_obj <- suppressMessages(suppressWarnings(",
  "            syn(od_rep, method = 'cart', m = 1, seed = syn_seed,",
  "                print.flag = FALSE, proper = TRUE, cart.minbucket = 10)",
  "          ))",
  "        }))",
  "        sd <- syn_obj$syn",
  "        sd$rep <- as.integer(r)",
  "        sd_list[[length(sd_list) + 1L]] <- sd",
  "      }",
  "",
  "      combined <- do.call(rbind, sd_list)",
  "      write_parquet(combined, out_path, compression = 'zstd')",
  "",
  "      results[[length(results) + 1L]] <- list(status = 'ok', file = sd_name)",
  "    }, error = function(e) {",
  "      results[[length(results) + 1L]] <<- list(status = 'error',",
  "                                                message = e$message, file = od_name)",
  "    })",
  "",
  "    file.remove(claimed_task)",
  "  }",
  "}",
  "",
  "writeLines('DONE', progress_file)",
  "saveRDS(results, result_file)"
)
writeLines(worker_code, worker_script)

# 6. Launch background workers
# ------------------------------------------------------------------------------
progress_files <- character(NUM_CORES)
result_files   <- character(NUM_CORES)

for (i in seq_len(NUM_CORES)) {
  progress_files[i] <- file.path(tempdir(), sprintf("sd_clust_prog_%02d.txt", i))
  result_files[i]   <- file.path(tempdir(), sprintf("sd_clust_res_%02d.rds", i))

  writeLines("STARTING", progress_files[i])
  if (file.exists(result_files[i])) file.remove(result_files[i])

  system2("Rscript", args = c(worker_script, progress_files[i], result_files[i],
                               todo_dir, doing_dir), wait = FALSE)
}

# 7. Live Unified Progress Display
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
  done_count  <- total_tasks - todo_count - doing_count

  pct     <- if (total_tasks > 0) as.integer(round(done_count / total_tasks * 100)) else 0
  bar     <- make_bar(done_count, total_tasks)
  elapsed <- as.integer(proc.time()[[3]] - start_ts)

  cat(sprintf("\r\033[2K[INFO] [%s] %3d%% (%d/%d) | %ds elapsed\n",
              bar, pct, done_count, total_tasks, elapsed))
  cat(strrep("-", 80), "\n")

  n_finish <- 0L

  for (i in seq_len(NUM_CORES)) {
    status <- "WAITING"
    if (file.exists(progress_files[i])) {
      lines <- readLines(progress_files[i], warn = FALSE)
      if (length(lines) > 0) status <- lines[1]
    }
    if (file.exists(result_files[i])) {
      status <- "DONE"; n_finish <- n_finish + 1L
    }
    if (nchar(status) > 65) status <- paste0(substr(status, 1, 62), "...")
    cat(sprintf("\r\033[2K   Core %-3d : %s\n", i, status))
  }

  if (n_finish == NUM_CORES) break
}

# 8. Collect results
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
    if (r$status == "ok") count <- count + 1L
    else {
      errors <- errors + 1L
      warns <- c(warns, sprintf("  [WARN] %s: %s", r$file, r$message))
    }
  }
}

if (length(warns) > 0L) cat("\n", paste(warns, collapse = "\n"), "\n", sep = "")

elapsed_total <- proc.time()[[3]] - start_ts
cat(sprintf("\n[DONE] %d file(s) written, %d errors.  Wall time: %.0fs\n",
            count, errors, elapsed_total))
