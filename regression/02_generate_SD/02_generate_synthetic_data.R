# ==============================================================================
# SCRIPT: 02_generate_synthetic_data.R
# PURPOSE: Generate Synthetic Data (SD) from each Original Data (OD) file
#          using MULTIPLE synthesis methods from the `synthpop` package.
#          ** FULL DYNAMIC QUEUE (WORK STEALING) APPLIED **
# ==============================================================================

# 1. Setup & Imports
# ------------------------------------------------------------------------------
suppressWarnings(suppressPackageStartupMessages({
  ensure_package <- function(pkg) {
    if (!(pkg %in% installed.packages()[,"Package"])) {
      message(paste("Installing package:", pkg))
      install.packages(pkg, repos = "http://cran.us.r-project.org",
                       dependencies = TRUE)
      if (!(pkg %in% installed.packages()[,"Package"])) {
        stop(paste("Failed to install package:", pkg))
      }
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
if (!file.exists(config_path)) {
  stop("config.json not found! Please ensure it is at ../config/config.json")
}

config       <- fromJSON(config_path)
cont_methods <- config$synthesis$continuous_methods
bin_methods  <- config$synthesis$binary_methods
MAX_ACCEPT_ATTEMPTS <- 1000L
VALIDATE_BATCH_SIZE <- 4L
timing_summary_path <- normalizePath(file.path("..", "results", "latest_sd_timing_breakdown.json"),
                                     winslash = "/", mustWork = FALSE)
python_bin <- Sys.which("python3")
if (python_bin == "") {
  stop("python3 not found on PATH")
}
validator_py <- normalizePath(file.path("..", "03_regression_analysis", "logit_acceptance.py"),
                              winslash = "/", mustWork = TRUE)

cat(sprintf("[INFO] Continuous methods: %s\n", paste(cont_methods, collapse = ", ")))
cat(sprintf("[INFO] Binary methods: %s\n", paste(bin_methods, collapse = ", ")))

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

cat(sprintf("[INFO] Found %d OD file(s) to synthesise.\n", length(od_files)))

# 4. Arm enumeration function
# ------------------------------------------------------------------------------
get_arms <- function(od_name, cont_methods, bin_methods) {
  x_is_binary <- grepl("_binary_",   od_name)
  y_is_binary <- grepl("logistic",   od_name)

  x_methods <- if (x_is_binary) bin_methods  else cont_methods
  y_methods <- if (y_is_binary) bin_methods  else cont_methods

  if (x_is_binary == y_is_binary) {
    # Case A or B: same type, enforce same method for all columns
    arms <- data.frame(x_method = x_methods, y_method = x_methods,
                       stringsAsFactors = FALSE)
  } else {
    # Case C or D: mixed type, full 2x2 factorial
    arms <- expand.grid(x_method = x_methods, y_method = y_methods,
                        stringsAsFactors = FALSE)
  }
  arms
}

# 5. Build Atomic Task Queue Directory structure
# ------------------------------------------------------------------------------
NUM_CORES <- 18L

todo_dir  <- file.path(tempdir(), "tasks_todo")
doing_dir <- file.path(tempdir(), "tasks_doing")
unlink(todo_dir, recursive = TRUE); unlink(doing_dir, recursive = TRUE)
dir.create(todo_dir, showWarnings = FALSE)
dir.create(doing_dir, showWarnings = FALSE)

set.seed(999)
shuffled_od_files <- sample(od_files) # Keep shuffle so binary/continuous mix!

task_i <- 1L
for (od_path in shuffled_od_files) {
  od_name <- basename(od_path)
  arms    <- get_arms(od_name, cont_methods, bin_methods)
  for (arm_idx in seq_len(nrow(arms))) {
    scenario_type <- if (grepl("logistic", od_name)) {
      if (grepl("_binary_", od_name)) "logistic_bin" else "logistic_cont"
    } else {
      if (grepl("_binary_", od_name)) "linear_bin" else "linear_cont"
    }
    task_data <- list(
      file     = od_path,
      x_method = arms$x_method[arm_idx],
      y_method = arms$y_method[arm_idx],
      seed     = as.integer(config$simulation$random_seed_base + task_i * 100000L),
      max_accept_attempts = MAX_ACCEPT_ATTEMPTS,
      validate_batch_size = VALIDATE_BATCH_SIZE,
      scenario_type = scenario_type
    )
    saveRDS(task_data, file.path(todo_dir, sprintf("task_%05d.rds", task_i)))
    task_i <- task_i + 1L
  }
}

total_tasks <- task_i - 1L
cat(sprintf("[INFO] Dynamic Queue built. %d tasks ready for %d cores.\n\n", total_tasks, NUM_CORES))

# 5. Generate Independent Worker Script
# ------------------------------------------------------------------------------
worker_script <- file.path(tempdir(), "sd_worker.R")
worker_code <- c(
  "args <- commandArgs(trailingOnly = TRUE)",
  "progress_file <- args[1]",
  "result_file <- args[2]",
  "todo_dir <- args[3]",
  "doing_dir <- args[4]",
  "python_bin <- args[5]",
  "validator_py <- args[6]",
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
  "full_od_files <- sort(list.files(file.path('..', 'data', 'original'), pattern = '^OD_.*[.]parquet$', full.names = TRUE))",
  "results <- list()",
  "writeLines('READY', progress_file)",
  "",
  "validate_candidates <- function(frame) {",
  "  parquet_path <- tempfile(pattern = 'sd_logit_validate_', fileext = '.parquet')",
  "  json_path    <- tempfile(pattern = 'sd_logit_validate_', fileext = '.json')",
  "  on.exit(unlink(c(parquet_path, json_path), force = TRUE), add = TRUE)",
  "  write_parquet(frame, parquet_path)",
  "  status <- system2(python_bin, args = c(validator_py, '--input', parquet_path, '--output', json_path), stdout = FALSE, stderr = FALSE)",
  "  if (!identical(status, 0L) || !file.exists(json_path)) stop('Logit validator failed')",
  "  out <- fromJSON(json_path)",
  "  as.data.frame(out, stringsAsFactors = FALSE)",
  "}",
  "",
  "# WORK STEALING LOOP: Atomically claim files from todo_dir",
  "repeat {",
  "  available_tasks <- list.files(todo_dir, full.names = TRUE)",
  "  if (length(available_tasks) == 0L) break # No more work!",
  "",
  "  target_task <- available_tasks[1]",
  "  claimed_task <- file.path(doing_dir, basename(target_task))",
  "",
  "  # file.rename is an OS-level ATOMIC operation.",
  "  # If it returns TRUE, this core successfully claimed the file.",
  "  if (suppressWarnings(file.rename(target_task, claimed_task))) {",
  "    task <- readRDS(claimed_task)",
  "    task_start <- proc.time()[[3]]",
  "    od_path <- task$file",
  "    od_name <- basename(od_path)",
  "",
  "    # Extract clean tag (remove 'OD_' and '.parquet')",
  "    scenario_tag <- sub('^OD_', '', sub('\\\\.parquet$', '', od_name))",
  "    arm_tag  <- paste0('X', toupper(task$x_method), '_Y', toupper(task$y_method))",
  "    sd_name  <- paste0('SD_', arm_tag, '_', scenario_tag)",
  "",
  "    # SKIP if output already exists (resume support)",
  "    out_path <- file.path(output_dir, paste0(sd_name, '.parquet'))",
  "    if (file.exists(out_path)) {",
  "      task_elapsed <- proc.time()[[3]] - task_start",
  "      results[[length(results) + 1L]] <- list(status = 'ok', file = sd_name, x_method = task$x_method, y_method = task$y_method, scenario_type = task$scenario_type, elapsed_s = task_elapsed)",
  "      file.remove(claimed_task)",
  "      next",
  "    }",
  "",
  "    writeLines(paste(arm_tag, '|', scenario_tag), progress_file)",
  "",
  "    tryCatch({",
  "      od_full <- as.data.frame(read_parquet(od_path))",
  "      global_idx <- match(od_path, full_od_files)",
  "",
  "      x_is_binary <- grepl('_binary_', od_name)",
  "      y_is_binary <- grepl('logistic',  od_name)",
  "      x_cols <- grep('^X[0-9]+$', names(od_full)[names(od_full) != 'iter'], value = TRUE)",
  "",
  "      iters <- sort(unique(od_full$iter))",
  "      sd_list <- list()",
  "      total_attempts <- 0L",
  "",
  "      for (it in iters) {",
  "        od <- od_full[od_full$iter == it, ]",
  "        od$iter <- NULL",
  "",
  "        # Build per-column method vector",
  "        method_vec <- setNames(rep(task$x_method, ncol(od)), names(od))",
  "        method_vec['y'] <- task$y_method",
  "",
  "        # Type conversions for synthpop",
  "        if (x_is_binary && task$x_method %in% c('cart', 'logreg')) {",
  "          for (col in x_cols) od[[col]] <- as.factor(od[[col]])",
  "        }",
  "        if (y_is_binary) {",
  "          if (task$y_method == 'logreg') {",
  "            od$y <- as.factor(od$y)",
  "          }",
  "        }",
  "",
  "        if (y_is_binary) {",
  "          accepted_sd <- NULL",
  "          attempt_n   <- 0L",
  "          generated_n <- 0L",
  "          while (is.null(accepted_sd)) {",
  "            if (attempt_n >= task$max_accept_attempts) {",
  "              stop(sprintf('Reached max attempts (%d) for SD iter %s in %s', task$max_accept_attempts, it, sd_name))",
  "            }",
  "            batch_n <- min(task$validate_batch_size, task$max_accept_attempts - attempt_n)",
  "            batch_list <- vector('list', batch_n)",
  "            for (b in seq_len(batch_n)) {",
  "              generated_n <- generated_n + 1L",
  "              syn_seed <- ((task$seed + (global_idx * 1000L) + (as.integer(it) * 100L) + generated_n - 2L) %% (.Machine$integer.max - 1L)) + 1L",
  "              invisible(capture.output({",
  "                syn_obj <- suppressMessages(suppressWarnings(",
  "                  syn(od, method = method_vec, seed = syn_seed, print.flag = FALSE, proper = TRUE)",
  "                ))",
  "              }))",
  "              sd <- syn_obj$syn",
  "              if (task$y_method == 'logreg') {",
  "                sd$y <- as.numeric(as.character(sd$y))",
  "              }",
  "              if (x_is_binary && task$x_method %in% c('cart', 'logreg')) {",
  "                for (col in x_cols) sd[[col]] <- as.numeric(as.character(sd[[col]]))",
  "              }",
  "              sd$iter <- as.integer(b)",
  "              batch_list[[b]] <- sd",
  "            }",
  "            batch_frame <- do.call(rbind, batch_list)",
  "            statuses <- validate_candidates(batch_frame)",
  "            accepted_iters <- statuses$iter[statuses$accepted]",
  "            if (length(accepted_iters) > 0L) {",
  "              first_accept <- accepted_iters[1]",
  "              attempt_n <- attempt_n + first_accept",
  "              accepted_sd <- batch_list[[first_accept]]",
  "              accepted_sd$iter <- it",
  "            } else {",
  "              attempt_n <- attempt_n + batch_n",
  "            }",
  "          }",
  "          total_attempts <- total_attempts + attempt_n",
  "          sd_list[[length(sd_list) + 1L]] <- accepted_sd",
  "        } else {",
  "          syn_seed <- base_seed + (global_idx * 1000L) + it",
  "          invisible(capture.output({",
  "            syn_obj <- suppressMessages(suppressWarnings(",
  "              syn(od, method = method_vec, seed = syn_seed, print.flag = FALSE, proper = TRUE)",
  "            ))",
  "          }))",
  "          sd <- syn_obj$syn",
  "          if (x_is_binary && task$x_method %in% c('cart', 'logreg')) {",
  "            for (col in x_cols) sd[[col]] <- as.numeric(as.character(sd[[col]]))",
  "          }",
  "          sd$iter <- it",
  "          sd_list[[length(sd_list) + 1L]] <- sd",
  "        }",
  "      }",
  "",
  "      combined_sd <- do.call(rbind, sd_list)",
  "      if (y_is_binary) {",
  "        combined_sd$attempt_total_sd <- as.integer(total_attempts)",
  "        combined_sd$target_success_sd <- as.integer(length(iters))",
  "      }",
  "      write_parquet(combined_sd, file.path(output_dir, paste0(sd_name, '.parquet')))",
  "",
  "      task_elapsed <- proc.time()[[3]] - task_start",
  "      results[[length(results) + 1L]] <- list(status = 'ok', file = sd_name, x_method = task$x_method, y_method = task$y_method, scenario_type = task$scenario_type, attempts = total_attempts, elapsed_s = task_elapsed)",
  "    }, error = function(e) {",
  "      task_elapsed <- proc.time()[[3]] - task_start",
  "      results[[length(results) + 1L]] <<- list(status = 'error', message = e$message, file = od_name, x_method = task$x_method, y_method = task$y_method, scenario_type = task$scenario_type, elapsed_s = task_elapsed)",
  "    })",
  "",
  "    # Delete the task file now that we are done with it",
  "    file.remove(claimed_task)",
  "  }",
  "}",
  "",
  "writeLines('DONE', progress_file)",
  "saveRDS(results, result_file)"
)
writeLines(worker_code, worker_script)

# 6. Launch completely independent background jobs
# ------------------------------------------------------------------------------
progress_files <- character(NUM_CORES)
result_files   <- character(NUM_CORES)

for (i in seq_len(NUM_CORES)) {
  progress_files[i] <- file.path(tempdir(), sprintf("sd_progress_%02d.txt", i))
  result_files[i]   <- file.path(tempdir(), sprintf("sd_result_%02d.rds", i))
  
  writeLines("STARTING", progress_files[i])
  if (file.exists(result_files[i])) file.remove(result_files[i])
  
  system2("Rscript", args = c(worker_script, progress_files[i], result_files[i], todo_dir, doing_dir, python_bin, validator_py), wait = FALSE)
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

# Initial UI Draw
cat(sprintf("Global Progress: [ Waiting... ]\n"))
cat(strrep("-", 80), "\n")
for (i in seq_len(NUM_CORES)) {
  cat(sprintf("   Core %-3d : STARTING\033[K\n", i))
}

repeat {
  Sys.sleep(0.5)
  
  # Move cursor up by (NUM_CORES + 2) lines to redraw
  cat(sprintf("\033[%dA", NUM_CORES + 2L))
  
  todo_count  <- length(list.files(todo_dir))
  doing_count <- length(list.files(doing_dir))
  done_count  <- total_tasks - todo_count - doing_count
  
  pct     <- if (total_tasks > 0) as.integer(round(done_count / total_tasks * 100)) else 0
  bar     <- make_bar(done_count, total_tasks)
  elapsed <- as.integer(proc.time()[[3]] - start_ts)
  
  # Print Unified Master Bar
  cat(sprintf("\r\033[2K[INFO] [%s] %3d%% (%d/%d) | %ds elapsed\n", 
              bar, pct, done_count, total_tasks, elapsed))
  cat(strrep("-", 80), "\n")
  
  n_finish <- 0L
  
  # Print Individual Core Activity
  for (i in seq_len(NUM_CORES)) {
    status <- "WAITING"
    
    if (file.exists(progress_files[i])) {
      lines <- readLines(progress_files[i], warn = FALSE)
      if (length(lines) > 0) status <- lines[1]
    }
    
    # If the result file exists, the background process definitely exited
    if (file.exists(result_files[i])) {
      status <- "DONE"
      n_finish <- n_finish + 1L
    }
    
    # Truncate very long scenario names so the terminal doesn't line-wrap and break the UI
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
timing_breakdown <- c(linear_cont = 0, linear_bin = 0, logistic_cont = 0, logistic_bin = 0)

for (i in seq_len(NUM_CORES)) {
  res_file <- result_files[i]
  
  if (!file.exists(res_file)) {
    errors <- errors + 1L
    warns <- c(warns, sprintf("  [FATAL] Core %d crashed entirely without saving output.", i))
    next
  }
  
  task_results <- readRDS(res_file)
  for (r in task_results) {
    if (!is.null(r$scenario_type) && !is.null(r$elapsed_s) && r$scenario_type %in% names(timing_breakdown)) {
      timing_breakdown[[r$scenario_type]] <- timing_breakdown[[r$scenario_type]] + as.numeric(r$elapsed_s)
    }
    if (r$status == "ok") {
      count <- count + 1L
    } else {
      errors <- errors + 1L
      warns  <- c(warns, sprintf("  [WARN] X%s_Y%s / %s: %s", r$x_method, r$y_method, r$file, r$message))
    }
  }
}

if (length(warns) > 0L) cat("\n", paste(warns, collapse = "\n"), "\n", sep = "")

elapsed_total <- proc.time()[[3]] - start_ts
total_worker_elapsed <- sum(timing_breakdown)
if (total_worker_elapsed > 0) {
  timing_breakdown_scaled <- timing_breakdown / total_worker_elapsed * elapsed_total
} else {
  timing_breakdown_scaled <- timing_breakdown
}
write(
  "",
  file = timing_summary_path
)
writeLines(
  toJSON(
    list(
      total_wall_s = unname(elapsed_total),
      total_worker_s = unname(total_worker_elapsed),
      scaled_wall_s = list(
        linear_cont = unname(timing_breakdown_scaled[["linear_cont"]]),
        linear_bin = unname(timing_breakdown_scaled[["linear_bin"]]),
        logistic_cont = unname(timing_breakdown_scaled[["logistic_cont"]]),
        logistic_bin = unname(timing_breakdown_scaled[["logistic_bin"]])
      ),
      raw_worker_s = list(
        linear_cont = unname(timing_breakdown[["linear_cont"]]),
        linear_bin = unname(timing_breakdown[["linear_bin"]]),
        logistic_cont = unname(timing_breakdown[["logistic_cont"]]),
        logistic_bin = unname(timing_breakdown[["logistic_bin"]])
      )
    ),
    auto_unbox = TRUE,
    pretty = TRUE
  ),
  con = timing_summary_path
)
cat(sprintf("\n[DONE] %d file(s) written, %d errors.  Wall time: %.0fs\n",
            count, errors, elapsed_total))
