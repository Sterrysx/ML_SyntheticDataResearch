# ==============================================================================
# SCRIPT: 02_generate_synthetic_data.R
# PURPOSE: Generate Synthetic Data (SD) from each Original Data (OD) file
#          using MULTIPLE synthesis methods from the `synthpop` package.
# ==============================================================================
#
# Methods: cart (CART), norm (normal linear regression), pmm (predictive
#          mean matching).  Each OD file produces one SD file per method.
#
# Parallelisation: 6 cores (Ryzen 9 9900X).  Each of the 3 methods is split
#                  into 2 equal file-chunks processed concurrently, so all 6
#                  workers carry the same workload.
#
# WHY CART IS SLOWER:
#   'cart' builds a full decision tree (via rpart) for every predictor column in
#   every dataset -- an O(n log n) recursive-partitioning step per variable.
#   'norm' (Bayesian linear regression) and 'pmm' (predictive mean matching) are
#   algebraically closed-form and are typically 4-10x faster.  The live display
#   below makes this difference visible in real time.
#
# Inputs:   ../data/original/OD_*_rep*.csv    (from Step 1)
# Outputs:  ../data/synthetic/SD_{method}_*_rep*.csv
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
ensure_package("synthpop")
ensure_package("parallel")

library(parallel)

# 2. Load Configuration
# ------------------------------------------------------------------------------
config_path <- "../config/config.json"
if (!file.exists(config_path)) {
  stop("config.json not found! Please ensure it is at ../config/config.json")
}

config      <- fromJSON(config_path)
base_seed   <- config$simulation$random_seed_base
syn_methods <- config$synthesis$methods   # e.g. ["cart", "norm", "pmm"]

cat(sprintf("[INFO] Synthesis methods: %s\n", paste(syn_methods, collapse = ", ")))

# 3. Discover OD files
# ------------------------------------------------------------------------------
input_dir  <- file.path("..", "data", "original")
output_dir <- file.path("..", "data", "synthetic")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

od_files <- sort(list.files(input_dir, pattern = "^OD_.*\\.csv$",
                            full.names = TRUE))

if (length(od_files) == 0) {
  stop("[ERROR] No OD files found in: ", input_dir,
       "\n        Run 01_generate_original_data.R first.")
}

total <- length(od_files) * length(syn_methods)
cat(sprintf("[INFO] Found %d OD file(s) x %d methods = %d SD file(s) to generate.\n",
            length(od_files), length(syn_methods), total))

# 4. Build parallel task list (Customized for Ryzen 9 9900X - 12 Cores)
# ------------------------------------------------------------------------------
# Manual distribution: 2 (norm) + 3 (cart) + 7 (pmm) = 12 total tasks
NUM_CORES <- 12L 

# Define your custom chunk requirements per method
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

# SHUFFLE the files here! This distributes "heavy" (N=1000) and "light" (N=100)
# files evenly across all chunks for ALL methods, perfectly balancing CPU load.
set.seed(999) # Keep chunk assignment deterministic
shuffled_od_files <- sample(od_files)

# Loop through the methods defined in your config
for (m_idx in seq_along(syn_methods)) {
  method <- syn_methods[m_idx]
  
  # Get the chunk count for this method (default to 2 if not in our list)
  n_chunks <- ifelse(!is.null(method_chunks[[method]]), method_chunks[[method]], 2L)
  
  # Use the shuffled list to create chunks
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

# 5. Worker function (runs inside each forked process)
# ------------------------------------------------------------------------------
worker <- function(task, progress_file) {
  suppressPackageStartupMessages({
    library(jsonlite)
    library(synthpop)
  })

  method    <- task$method
  m_idx     <- task$m_idx
  chunk_idx <- task$chunk_idx
  files     <- task$files
  n_files   <- length(files)

  # Signal to the monitor: initialised, total known.
  writeLines(sprintf("0 %d", n_files), progress_file)

  config_path <- "../config/config.json"
  config      <- fromJSON(config_path)
  base_seed   <- config$simulation$random_seed_base
  syn_methods <- config$synthesis$methods
  output_dir  <- file.path("..", "data", "synthetic")

  full_od_files <- sort(list.files(
    file.path("..", "data", "original"),
    pattern    = "^OD_.*\\.csv$",
    full.names = TRUE
  ))

  results <- list()
  done    <- 0L

  for (od_path in files) {
    od_name      <- basename(od_path)
    od           <- read.csv(od_path)
    global_idx   <- match(od_path, full_od_files)

    syn_seed     <- base_seed + (global_idx * 1000L) + match(method, syn_methods) * 100L
    scenario_tag <- sub("^OD_", "", od_name)

    tryCatch({
      # suppressMessages + suppressWarnings silence synthpop's per-variable
      # "Method changed to logreg/polyreg" and "numeric turned to factor"
      # messages that flood the terminal when processing binary OD files.
      # The synthesis itself is unaffected; synthpop still picks the correct
      # method (logreg for 0/1 columns) automatically.
      # capture.output() catches any residual cat() calls from synthpop that
      # would otherwise corrupt the ANSI progress bar display.
      # minnumlevels = 1 suppresses the warning about binary 0/1 columns
      # being numeric with 5 or fewer unique levels.
      invisible(capture.output({
        syn_obj <- suppressMessages(suppressWarnings(
          syn(od, method = method, seed = syn_seed,
              print.flag = FALSE, proper = TRUE, minnumlevels = 1)
        ))
      }))
      sd <- syn_obj$syn

      # When method == "norm", synthpop treats binary predictors as continuous
      # and generates decimal values.  Threshold X1..Xp back to strict 0/1
      # while leaving the continuous y column untouched.
      if (method == "norm") {
        x_cols <- grep("^X[0-9]+$", names(sd), value = TRUE)
        for (col in x_cols) {
          sd[[col]] <- ifelse(sd[[col]] >= 0.5, 1, 0)
        }
      }

      sd_name <- paste0("SD_", method, "_", scenario_tag)
      write.csv(sd, file.path(output_dir, sd_name), row.names = FALSE)

      done <- done + 1L
      results[[length(results) + 1L]] <- list(
        status = "ok",
        file   = sd_name,
        method = method,
        chunk  = chunk_idx
      )
    }, error = function(e) {
      results[[length(results) + 1L]] <<- list(
        status  = "warn",
        file    = od_name,
        method  = method,
        chunk   = chunk_idx,
        message = e$message
      )
    })

    # Update progress file after every dataset (polled by the monitor loop).
    writeLines(sprintf("%d %d", done, n_files), progress_file)
  }

  results
}

# 6. Create per-task progress files and launch workers (non-blocking)
# ------------------------------------------------------------------------------
progress_files <- vapply(
  seq_along(tasks),
  function(i) file.path(tempdir(), sprintf("sd_progress_%02d.txt", i)),
  character(1)
)
for (pf in progress_files) writeLines("0 0", pf)

# mcparallel forks immediately and returns a handle -- does NOT block.
jobs <- lapply(seq_along(tasks), function(i) {
  mcparallel(worker(tasks[[i]], progress_files[i]))
})

# 7. Live per-core progress display (polls every 0.5 s)
# ------------------------------------------------------------------------------
BAR_WIDTH <- 24L

make_bar <- function(done, total) {
  if (total == 0L) return(strrep("-", BAR_WIDTH))
  filled <- as.integer(round(done / total * BAR_WIDTH))
  paste0(strrep("=", filled), strrep("-", BAR_WIDTH - filled))
}

read_progress <- function(pf) {
  tryCatch({
    parts <- as.integer(
      strsplit(trimws(readLines(pf, n = 1L, warn = FALSE)), " ")[[1]]
    )
    if (length(parts) == 2L) parts else c(0L, 0L)
  }, error = function(e) c(0L, 0L))
}

n_tasks  <- length(tasks)
start_ts <- proc.time()[[3]]

# Print header then one placeholder line per core (overwritten in-place).
cat(sprintf("%-8s  %-6s  %-7s  %-*s  %4s  %-10s\n",
            "Core", "Method", "Chunk", BAR_WIDTH, "Progress", "Done", "Files"))
cat(strrep("-", 65L), "\n")
for (i in seq_len(n_tasks)) {
  cat(sprintf("Core %-3d  %-6s  chunk-%1d  [%s]   0%%  (0/??)\n",
              i, tasks[[i]]$method, tasks[[i]]$chunk_idx,
              strrep("-", BAR_WIDTH)))
}

repeat {
  Sys.sleep(0.5)

  elapsed  <- proc.time()[[3]] - start_ts
  n_finish <- 0L

  # Jump cursor up over the n_tasks status lines to overwrite them.
  cat(sprintf("\033[%dA", n_tasks))

  for (i in seq_len(n_tasks)) {
    p     <- read_progress(progress_files[i])
    done  <- p[1L]; tot <- p[2L]
    pct   <- if (tot > 0L) as.integer(round(done / tot * 100L)) else 0L
    bar   <- make_bar(done, tot)
    tot_s <- if (tot > 0L) as.character(tot) else "??"
    fin   <- tot > 0L && done >= tot
    if (fin) n_finish <- n_finish + 1L
    tag   <- if (fin) " DONE " else sprintf(" %3ds ", as.integer(elapsed))
    # \033[2K clears the current terminal line before writing the new content.
    cat(sprintf("\033[2KCore %-3d  %-6s  chunk-%1d  [%s]  %3d%%  (%d/%s)%s\n",
                i, tasks[[i]]$method, tasks[[i]]$chunk_idx,
                bar, pct, done, tot_s, tag))
  }

  if (n_finish == n_tasks) break
}

# 8. Collect results (blocking -- all workers are finished by now)
# ------------------------------------------------------------------------------
all_results <- mccollect(jobs, wait = TRUE)

count  <- 0L
errors <- 0L
warns  <- character(0)

for (task_results in all_results) {
  if (inherits(task_results, "try-error")) {
    cat(sprintf("[ERROR] A worker process failed: %s\n", as.character(task_results)))
    next
  }
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
cat(sprintf("\n[DONE] %d file(s) written, %d warning(s).  Wall time: %.0fs\n",
            count, errors, elapsed_total))
