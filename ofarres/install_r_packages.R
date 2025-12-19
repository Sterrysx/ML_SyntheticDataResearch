#!/usr/bin/env Rscript
# ==============================================================================
# SCRIPT: install_r_packages.R
# PURPOSE: Install required R packages for the pipeline
# ==============================================================================

cat("Installing required R packages...\n\n")

# List of required packages
packages <- c("jsonlite", "mvtnorm", "dplyr", "synthpop")

# Function to install if missing
install_if_missing <- function(pkg) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(paste("Installing", pkg, "...\n"))
    install.packages(pkg, repos = "http://cran.us.r-project.org", dependencies = TRUE)
    
    # Verify installation
    if (require(pkg, character.only = TRUE, quietly = TRUE)) {
      cat(paste("✓", pkg, "installed successfully\n\n"))
    } else {
      cat(paste("✗ Failed to install", pkg, "\n\n"))
    }
  } else {
    cat(paste("✓", pkg, "already installed\n\n"))
  }
}

# Install all packages
for (pkg in packages) {
  install_if_missing(pkg)
}

cat("\n===========================================\n")
cat("Package installation complete!\n")
cat("===========================================\n")
