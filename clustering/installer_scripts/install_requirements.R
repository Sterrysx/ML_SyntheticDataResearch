# ==============================================================================
# R Package Installation Script for Synthetic Data Clustering Pipeline
# ==============================================================================
# Purpose: Install all required R packages for the data generation step
# Usage:   Rscript install_requirements.R
# ==============================================================================

cat("========================================\n")
cat("R Package Installation\n")
cat("========================================\n\n")

# Define required packages
pkg_list <- c(
  "jsonlite",   # JSON parsing for config.json
  "mvtnorm",    # Multivariate normal distribution generation
  "synthpop",   # Synthetic data generation (CART method)
  "dplyr"       # Data manipulation
)

cat("Required packages:", paste(pkg_list, collapse=", "), "\n\n")

# Check which packages are missing
installed_pkgs <- installed.packages()[, "Package"]
new_pkgs <- pkg_list[!(pkg_list %in% installed_pkgs)]

if (length(new_pkgs) > 0) {
  cat("Installing missing packages:", paste(new_pkgs, collapse=", "), "\n\n")
  
  # Install missing packages
  install.packages(
    new_pkgs, 
    repos = "http://cran.us.r-project.org",
    dependencies = TRUE,
    quiet = FALSE
  )
  
  cat("\n✅ Installation complete!\n")
} else {
  cat("✅ All required packages are already installed!\n")
}

# Verify all packages can be loaded
cat("\nVerifying package installation...\n")
success <- TRUE

for (pkg in pkg_list) {
  if (require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat("  ✓", pkg, "\n")
  } else {
    cat("  ✗", pkg, "- FAILED\n")
    success <- FALSE
  }
}

if (success) {
  cat("\n✅ All packages verified successfully!\n")
  cat("========================================\n")
  cat("Setup complete. Ready to run pipeline.\n")
  cat("========================================\n")
} else {
  cat("\n❌ Some packages failed to load.\n")
  cat("Please install them manually:\n")
  cat("install.packages(c(", paste0("'", pkg_list, "'", collapse=", "), "))\n")
  quit(status = 1)
}