# ==============================================================================
# SCRIPT: 02_generate_synthetic_data.R
# PURPOSE: Generate Synthetic Data (SD) from each Original Data (OD) file
#          using the CART method from the `synthpop` package.
#          Now supports M repetitions per rho value.
# ==============================================================================
#
# Inputs:   ../data/original/OD_*_rep*.csv    (from Step 1)
# Outputs:  ../data/synthetic/SD_*_rep*.csv
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

# 2. Load Configuration
# ------------------------------------------------------------------------------
config_path <- "../config/config.json"
if (!file.exists(config_path)) {
  stop("config.json not found! Please ensure it is at ../config/config.json")
}

config    <- fromJSON(config_path)
base_seed <- config$simulation$random_seed_base
syn_method <- config$synthesis$method   # "cart"

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

cat(sprintf("[INFO] Found %d OD file(s) to synthesize (method = '%s').\n\n",
            length(od_files), syn_method))

# 4. Synthesis Loop
# ------------------------------------------------------------------------------
total <- length(od_files)
cat(sprintf("[INFO] Processing %d OD file(s)...\n\n", total))

for (idx in seq_along(od_files)) {

  od_path <- od_files[idx]
  od_name <- basename(od_path)
  od      <- read.csv(od_path)

  # Print progress periodically
  if (idx == 1 || idx == total || idx %% 10 == 0) {
    cat(sprintf("[INFO] Synthesizing (%d/%d): %s  (dim: %d x %d)\n",
                idx, total, od_name, nrow(od), ncol(od)))
  }

  syn_seed <- base_seed + (idx * 1000)

  # Generate synthetic data using CART for all variables
  syn_obj <- syn(od, method = syn_method, seed = syn_seed, print.flag = FALSE)

  # Extract the synthetic data frame
  sd <- syn_obj$syn

  # Build output filename: OD_... -> SD_...
  sd_name <- sub("^OD_", "SD_", od_name)
  sd_path <- file.path(output_dir, sd_name)

  write.csv(sd, sd_path, row.names = FALSE)
}

cat(sprintf("\n[DONE] Synthetic data generation complete. %d file(s) written.\n", total))
