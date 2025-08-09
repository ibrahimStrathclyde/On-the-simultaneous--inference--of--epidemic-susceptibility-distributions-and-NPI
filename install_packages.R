# install_packages.R
# Script to install all required packages for epidemic modeling analysis
# 
# This script checks for and installs all R packages required for the
# simultaneous inference of susceptibility distributions and intervention
# effects from epidemic curves analysis.
#
# Author: Ibrahim Mohammed, Chris Robertson, M. Gabriela M. Gomes
# Date: August 2025
# Repository: epidemic-susceptibility-inference

# Function to install packages if not already installed
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(paste("Installing package:", pkg, "...\n"))
    install.packages(pkg, dependencies = TRUE)
    
    # Verify installation
    if (requireNamespace(pkg, quietly = TRUE)) {
      cat(paste("✓ Package", pkg, "successfully installed.\n"))
    } else {
      cat(paste("✗ Failed to install package:", pkg, "\n"))
      cat("Please try installing manually with: install.packages('", pkg, "')\n", sep = "")
    }
  } else {
    cat(paste("✓ Package", pkg, "already installed.\n"))
  }
}

# Display welcome message
cat("============================================================\n")
cat("EPIDEMIC MODELING ANALYSIS - PACKAGE INSTALLATION\n")
cat("============================================================\n")
cat("Installing required packages for SEIR model analysis...\n\n")

# List of required packages with descriptions
required_packages <- list(
  # Core mathematical and statistical packages
  "deSolve" = "Solving differential equations (SEIR models)",
  "maxLik" = "Maximum likelihood estimation",
  "MASS" = "Statistical functions and distributions",
  
  # Data manipulation and tidying (tidyverse includes dplyr, tidyr, ggplot2, etc.)
  "tidyverse" = "Collection of data science packages (includes dplyr, ggplot2, tidyr)",
  "dplyr" = "Data manipulation and transformation",
  "tidyr" = "Data tidying and reshaping",
  "ggplot2" = "Advanced graphics and plotting",
  
  # Specialized visualization and plotting
  "gridExtra" = "Arranging multiple plots",
  "grid" = "Low-level grid graphics functions",
  "RColorBrewer" = "Color palettes for plots",
  "viridis" = "Perceptually uniform color scales",
  "GGally" = "Extension of ggplot2 for correlation matrices and pair plots",
  
  # Documentation and reporting
  "knitr" = "Dynamic report generation and documentation",
  
  # Additional utilities  
  "reshape2" = "Data reshaping (melt/cast functions)"
)

# Install packages
cat("Checking and installing required packages:\n")
cat("------------------------------------------\n")

for (pkg_name in names(required_packages)) {
  cat(sprintf("%-15s - %s\n", pkg_name, required_packages[[pkg_name]]))
  install_if_missing(pkg_name)
  cat("\n")
}

# Verify all packages can be loaded
cat("============================================================\n")
cat("VERIFICATION: Testing package loading...\n")
cat("============================================================\n")

all_loaded <- TRUE
failed_packages <- character(0)

for (pkg_name in names(required_packages)) {
  tryCatch({
    library(pkg_name, character.only = TRUE, quietly = TRUE)
    cat(paste("✓", pkg_name, "loaded successfully\n"))
  }, error = function(e) {
    cat(paste("✗", pkg_name, "failed to load:", e$message, "\n"))
    all_loaded <<- FALSE
    failed_packages <<- c(failed_packages, pkg_name)
  })
}

# Final status report
cat("\n============================================================\n")
if (all_loaded) {
  cat("SUCCESS: All required packages are installed and working!\n")
  cat("============================================================\n")
  cat("\nYou can now run the analysis scripts:\n")
  cat("1. source('scripts/1_baseline_cases.R')\n")
  cat("2. source('scripts/2_mle_single_epidemic.R')\n")
  cat("3. source('scripts/3_mle_two_epidemics.R')\n")
  cat("4. source('scripts/4_single_epidemic_correlation.R')\n")
  cat("5. source('scripts/5_two_epidemics_correlation.R')\n")
} else {
  cat("WARNING: Some packages failed to install or load properly!\n")
  cat("============================================================\n")
  cat("\nFailed packages:", paste(failed_packages, collapse = ", "), "\n")
  cat("\nTroubleshooting steps:\n")
  cat("1. Check your internet connection\n")
  cat("2. Update R to the latest version\n")
  cat("3. Try installing failed packages manually:\n")
  for (pkg in failed_packages) {
    cat(sprintf("   install.packages('%s')\n", pkg))
  }
  cat("4. Check for system dependencies (especially for mathematical packages)\n")
  cat("5. Consider using install.packages() with repos parameter specified\n")
}

cat("\nFor additional help, see the repository README.md file.\n")

# Clean up
rm(install_if_missing, required_packages, all_loaded, failed_packages)

# Optional: Display session information for debugging
cat("\n============================================================\n")
cat("SESSION INFORMATION (for troubleshooting):\n")
cat("============================================================\n")
print(sessionInfo())