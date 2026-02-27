# ==============================================================================
# SCRIPT: 02_generate_synthetic_data.R
# PURPOSE: Generate Synthetic Data (SD) from each Original Data (OD) file
#          using MULTIPLE synthesis methods from the `synthpop` package.
# ==============================================================================
#
# Methods: cart (CART), norm (normal linear regression), pmm (predictive
#          mean matching).  Each OD file produces one SD file per method.
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
cat(sprintf("[INFO] Found %d OD file(s) x %d methods = %d SD file(s) to generate.\n\n",
            length(od_files), length(syn_methods), total))

# 4. Synthesis Loop (each OD × each method)
# ------------------------------------------------------------------------------
count <- 0

for (idx in seq_along(od_files)) {

  od_path <- od_files[idx]
  od_name <- basename(od_path)
  od      <- read.csv(od_path)

  # Strip leading "OD_" to get the scenario tag
  scenario_tag <- sub("^OD_", "", od_name)   # e.g. N1000_p10_rho0.0_sigma1.0_rep1.csv

  for (method in syn_methods) {
    syn_seed <- base_seed + (idx * 1000) + match(method, syn_methods) * 100

    tryCatch({
      syn_obj <- syn(od, method = method, seed = syn_seed,
                     print.flag = FALSE, proper = TRUE)
      sd <- syn_obj$syn

      sd_name <- paste0("SD_", method, "_", scenario_tag)
      sd_path <- file.path(output_dir, sd_name)
      write.csv(sd, sd_path, row.names = FALSE)

      count <- count + 1
      if (count == 1 || count == total || count %% 50 == 0) {
        cat(sprintf("[INFO] (%d/%d) %s\n", count, total, sd_name))
      }
    }, error = function(e) {
      cat(sprintf("[WARN] Synthesis failed (%s, %s): %s\n", method, od_name, e$message))
    })
  }
}

cat(sprintf("\n[DONE] Synthetic data generation complete. %d file(s) written.\n", count))
