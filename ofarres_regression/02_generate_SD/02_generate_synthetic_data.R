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

config      <- fromJSON(config_path)
syn_methods <- config$synthesis$methods    

cat(sprintf("[INFO] Synthesis methods: %s\n", paste(syn_methods, collapse = ", ")))

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

total <- length(od_files) * length(syn_methods)
cat(sprintf("[INFO] Found %d OD file(s) x %d methods = %d SD file(s) to generate.\n",
            length(od_files), length(syn_methods), total))

# 4. Build Atomic Task Queue Directory structure
# ------------------------------------------------------------------------------
NUM_CORES <- 18L 

todo_dir  <- file.path(tempdir(), "tasks_todo")
doing_dir <- file.path(tempdir(), "tasks_doing")
unlink(todo_dir, recursive = TRUE); unlink(doing_dir, recursive = TRUE)
dir.create(todo_dir, showWarnings = FALSE)
dir.create(doing_dir, showWarnings = FALSE)

task_i <- 1L
set.seed(999) 
shuffled_od_files <- sample(od_files) # Keep shuffle so binary/continuous mix!

for (method in syn_methods) {
  for (od_path in shuffled_od_files) {
    task_data <- list(method = method, file = od_path)
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
  "syn_methods <- config$synthesis$methods",
  "output_dir <- file.path('..', 'data', 'synthetic')",
  "",
  "full_od_files <- sort(list.files(file.path('..', 'data', 'original'), pattern = '^OD_.*[.]parquet$', full.names = TRUE))",
  "results <- list()",
  "writeLines('READY', progress_file)",
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
  "  if (file.rename(target_task, claimed_task)) {",
  "    task <- readRDS(claimed_task)",
  "    od_path <- task$file",
  "    method <- task$method",
  "    od_name <- basename(od_path)",
  "",
  "    # Extract clean tag (remove 'OD_' and '.parquet')",
  "    scenario_tag <- sub('^OD_', '', sub('\\\\.parquet$', '', od_name))",
  "    sd_name <- paste0('SD_', method, '_', scenario_tag)",
  "",
  "    # SKIP if output already exists (resume support)",
  "    out_path <- file.path(output_dir, paste0(sd_name, '.parquet'))",
  "    if (file.exists(out_path)) {",
  "      results[[length(results) + 1L]] <- list(status = 'ok', file = sd_name, method = method)",
  "      file.remove(claimed_task)",
  "      next",
  "    }",
  "",
  "    writeLines(paste(toupper(method), '|', scenario_tag), progress_file)",
  "",
  "    tryCatch({",
  "      od_full <- as.data.frame(read_parquet(od_path))",
  "      global_idx <- match(od_path, full_od_files)",
  "",
  "      is_binary <- grepl('_binary_', od_name)",
  "      x_cols <- grep('^X[0-9]+$', names(od_full)[names(od_full) != 'iter'], value = TRUE)",
  "",
  "      iters <- sort(unique(od_full$iter))",
  "      sd_list <- list()",
  "",
  "      for (it in iters) {",
  "        od <- od_full[od_full$iter == it, ]",
  "        od$iter <- NULL",
  "",
  "        syn_seed <- base_seed + (global_idx * 1000L) + match(method, syn_methods) * 100L + it",
  "",
  "        if (is_binary && method == 'cart') {",
  "          for (col in x_cols) od[[col]] <- as.factor(od[[col]])",
  "        }",
  "",
  "        invisible(capture.output({",
  "          syn_obj <- suppressMessages(suppressWarnings(",
  "            syn(od, method = method, seed = syn_seed, print.flag = FALSE, proper = TRUE)",
  "          ))",
  "        }))",
  "        sd <- syn_obj$syn",
  "",
  "        if (is_binary) {",
  "          for (col in x_cols) {",
  "            if (method == 'cart') {",
  "              sd[[col]] <- as.numeric(as.character(sd[[col]]))",
  "            } else {",
  "              sd[[col]] <- as.numeric(sd[[col]] >= 0.5)",
  "            }",
  "          }",
  "        }",
  "",
  "        sd$iter <- it",
  "        sd_list[[length(sd_list) + 1L]] <- sd",
  "      }",
  "",
  "      combined_sd <- do.call(rbind, sd_list)",
  "      write_parquet(combined_sd, file.path(output_dir, paste0(sd_name, '.parquet')))",
  "",
  "      results[[length(results) + 1L]] <- list(status = 'ok', file = sd_name, method = method)",
  "    }, error = function(e) {",
  "      results[[length(results) + 1L]] <<- list(status = 'error', message = e$message, file = od_name, method = method)",
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
  
  system2("Rscript", args = c(worker_script, progress_files[i], result_files[i], todo_dir, doing_dir), wait = FALSE)
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
      warns  <- c(warns, sprintf("  [WARN] %s / %s: %s", r$method, r$file, r$message))
    }
  }
}

if (length(warns) > 0L) cat("\n", paste(warns, collapse = "\n"), "\n", sep = "")

elapsed_total <- proc.time()[[3]] - start_ts
cat(sprintf("\n[DONE] %d file(s) written, %d errors.  Wall time: %.0fs\n",
            count, errors, elapsed_total))