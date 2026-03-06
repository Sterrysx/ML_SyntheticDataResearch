# ==============================================================================
# SCRIPT: 02_generate_synthetic_data.R
# PURPOSE: Generate Synthetic Data (SD) from each Original Data (OD) file
#          using MULTIPLE synthesis methods from the `synthpop` package.
#          ** FORK-FREE ARCHITECTURE: COMPLETELY IMMUNE TO DEADLOCKS **
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

# 4. Build task list
# ------------------------------------------------------------------------------
NUM_CORES <- 18L 

method_chunks <- list(
  "norm" = 2L,
  "cart" = 3L,
  "pmm"  = 13L
)

split_files <- function(files, n_chunks) {
  if (n_chunks <= 1) return(list(files))
  idx <- seq_along(files)
  split(files, cut(idx, breaks = n_chunks, labels = FALSE))
}

tasks  <- list()
task_i <- 1L

set.seed(999) 
shuffled_od_files <- sample(od_files)

for (m_idx in seq_along(syn_methods)) {
  method <- syn_methods[m_idx]
  n_chunks <- ifelse(!is.null(method_chunks[[method]]), method_chunks[[method]], 2L)
  chunks <- split_files(shuffled_od_files, n_chunks)
  
  for (chunk_idx in seq_len(n_chunks)) {
    tasks[[task_i]] <- list(
      method    = method,
      m_idx     = m_idx,
      chunk_idx = chunk_idx,
      files     = chunks[[chunk_idx]]
    )
    task_i <- task_i + 1L
  }
}

cat(sprintf("[INFO] Custom Dispatch: %d tasks across %d cores.\n", length(tasks), NUM_CORES))
cat(sprintf("[INFO] Distribution: Norm(2), Cart(3), PMM(13)\n\n"))

# 5. Generate Independent Worker Script
# ------------------------------------------------------------------------------
worker_script <- file.path(tempdir(), "sd_worker.R")
worker_code <- c(
  "args <- commandArgs(trailingOnly = TRUE)",
  "task_file <- args[1]",
  "task <- readRDS(task_file)",
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
  "method <- task$method",
  "files <- task$files",
  "progress_file <- task$progress_file",
  "result_file <- task$result_file",
  "n_files <- length(files)",
  "",
  "writeLines(sprintf('0 %d', n_files), progress_file)",
  "",
  "config <- fromJSON('../config/config.json')",
  "base_seed <- config$simulation$random_seed_base",
  "syn_methods <- config$synthesis$methods",
  "output_dir <- file.path('..', 'data', 'synthetic')",
  "",
  "full_od_files <- sort(list.files(file.path('..', 'data', 'original'), pattern = '^OD_.*[.]parquet$', full.names = TRUE))",
  "results <- list()",
  "done <- 0L",
  "",
  "for (od_path in files) {",
  "  od_name <- basename(od_path)",
  "  tryCatch({",
  "    od <- as.data.frame(read_parquet(od_path))",
  "    global_idx <- match(od_path, full_od_files)",
  "    syn_seed <- base_seed + (global_idx * 1000L) + match(method, syn_methods) * 100L",
  "    scenario_tag <- sub('^OD_', '', od_name)",
  "",
  "    is_binary <- grepl('_binary_', od_name)",
  "    x_cols <- grep('^X[0-9]+$', names(od), value = TRUE)",
  "",
  "    # Only CART requires binary variables to be factors. NORM/PMM crash on factors.",
  "    if (is_binary && method == 'cart') {",
  "      for (col in x_cols) od[[col]] <- as.factor(od[[col]])",
  "    }",
  "",
  "    invisible(capture.output({",
  "      syn_obj <- suppressMessages(suppressWarnings(",
  "        syn(od, method = method, seed = syn_seed, print.flag = FALSE, proper = TRUE)",
  "      ))",
  "    }))",
  "    sd <- syn_obj$syn",
  "",
  "    if (is_binary) {",
  "      for (col in x_cols) {",
  "        if (method == 'cart') {",
  "          sd[[col]] <- as.numeric(as.character(sd[[col]]))",
  "        } else {",
  "          # NORM/PMM output continuous predictions for binary inputs, snap them back to 0/1",
  "          sd[[col]] <- as.numeric(sd[[col]] >= 0.5)",
  "        }",
  "      }",
  "    }",
  "",
  "    sd_name <- paste0('SD_', method, '_', scenario_tag)",
  "    write_parquet(sd, file.path(output_dir, sd_name))",
  "",
  "    done <- done + 1L",
  "    results[[length(results) + 1L]] <- list(status = 'ok', file = sd_name, method = method)",
  "    writeLines(sprintf('%d %d', done, n_files), progress_file)",
  "  }, error = function(e) {",
  "    done <<- done + 1L",
  "    writeLines(sprintf('%d %d', done, n_files), progress_file)",
  "    results[[length(results) + 1L]] <<- list(status = 'error', message = e$message, file = od_name, method = method)",
  "  })",
  "}",
  "",
  "saveRDS(results, result_file)"
)
writeLines(worker_code, worker_script)

# 6. Launch completely independent background jobs
# ------------------------------------------------------------------------------
progress_files <- character(length(tasks))

for (i in seq_along(tasks)) {
  tasks[[i]]$progress_file <- file.path(tempdir(), sprintf("sd_progress_%02d.txt", i))
  tasks[[i]]$result_file   <- file.path(tempdir(), sprintf("sd_result_%02d.rds", i))
  progress_files[i]        <- tasks[[i]]$progress_file
  
  # Ensure clean slate
  writeLines("0 0", progress_files[i])
  if (file.exists(tasks[[i]]$result_file)) file.remove(tasks[[i]]$result_file)
  
  # Package up the task and launch the background process
  task_file <- file.path(tempdir(), sprintf("task_%02d.rds", i))
  saveRDS(tasks[[i]], task_file)
  
  system2("Rscript", args = c(worker_script, task_file), wait = FALSE)
}

# 7. Live per-core progress display
# ------------------------------------------------------------------------------
BAR_WIDTH <- 24L

make_bar <- function(done, total) {
  if (total == 0L) return(strrep("-", BAR_WIDTH))
  filled <- as.integer(round(done / total * BAR_WIDTH))
  paste0(strrep("=", filled), strrep("-", BAR_WIDTH - filled))
}

read_progress <- function(pf) {
  tryCatch({
    parts <- as.integer(strsplit(trimws(readLines(pf, n = 1L, warn = FALSE)), " ")[[1]])
    if (length(parts) == 2L) parts else c(0L, 0L)
  }, error = function(e) c(0L, 0L))
}

n_tasks  <- length(tasks)
start_ts <- proc.time()[[3]]

cat(sprintf("%-8s  %-6s  %-7s  %-*s  %4s  %-10s\n",
            "Core", "Method", "Chunk", BAR_WIDTH, "Progress", "Done", "Files"))
cat(strrep("-", 65L), "\n")
for (i in seq_len(n_tasks)) {
  cat(sprintf("Core %-3d  %-6s  chunk-%1d  [%s]   0%%  (0/??)\n",
              i, tasks[[i]]$method, tasks[[i]]$chunk_idx, strrep("-", BAR_WIDTH)))
}

repeat {
  Sys.sleep(0.5)

  elapsed  <- proc.time()[[3]] - start_ts
  n_finish <- 0L

  cat(sprintf("\033[%dA", n_tasks))

  for (i in seq_len(n_tasks)) {
    p <- read_progress(progress_files[i])
    done <- p[1L]; tot <- p[2L]
    
    # The existence of the result file is the absolute truth that the process exited
    fin <- file.exists(tasks[[i]]$result_file)
    
    if (fin) {
      n_finish <- n_finish + 1L
      if (tot == 0L) tot <- 1L 
      done <- tot              
    }
    
    pct   <- if (tot > 0L) as.integer(round(done / tot * 100L)) else 0L
    bar   <- make_bar(done, tot)
    tot_s <- if (tot > 0L) as.character(tot) else "??"
    tag   <- if (fin) " DONE " else sprintf(" %3ds ", as.integer(elapsed))
    
    cat(sprintf("\033[2KCore %-3d  %-6s  chunk-%1d  [%s]  %3d%%  (%d/%s)%s\n",
                i, tasks[[i]]$method, tasks[[i]]$chunk_idx,
                bar, pct, done, tot_s, tag))
  }

  if (n_finish == n_tasks) break
}

# 8. Collect results 
# ------------------------------------------------------------------------------
count  <- 0L
errors <- 0L
warns  <- character(0)

for (i in seq_along(tasks)) {
  res_file <- tasks[[i]]$result_file
  
  if (!file.exists(res_file)) {
    errors <- errors + 1L
    warns <- c(warns, sprintf("  [FATAL] Core %d (%s) crashed entirely without saving output.", i, tasks[[i]]$method))
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