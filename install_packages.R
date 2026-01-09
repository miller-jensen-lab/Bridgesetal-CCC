# Install required R packages for Bridges et al. Fig 6 DBiT-seq analysis
# Run this script once before running Fig6-DBiT-run.R

# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

cat("=== Installing packages for Bridges et al. DBiT-seq analysis ===\n\n")

# Function to install if missing
install_if_missing <- function(pkg, ...) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("Installing %s...\n", pkg))
    install.packages(pkg, quiet = TRUE, ...)
  } else {
    cat(sprintf("%s already installed\n", pkg))
  }
}

# CRAN packages
cran_packages <- c(
  # Core tidyverse
  "tidyverse",
  "magrittr",

  # Data I/O
  "rio",
  "janitor",
  "openxlsx",
  "readxl",
  "data.table",

 # Plotting
  "ggplot2",
  "ggpubr",
  "gridExtra",
  "wesanderson",
  "RColorBrewer",

  # Spatial/image
  "raster",
  "OpenImageR",

  # Other
  "plyr",
  "Matrix",
  "remotes"
)

cat("Installing CRAN packages...\n")
for (pkg in cran_packages) {
  install_if_missing(pkg)
}

# Bioconductor packages (if needed)
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  cat("Installing BiocManager...\n")
  install.packages("BiocManager", quiet = TRUE)
}

# Seurat (from GitHub for latest version)
if (!requireNamespace("Seurat", quietly = TRUE)) {
  cat("\nInstalling Seurat from GitHub...\n")
  remotes::install_github("satijalab/seurat", "seurat5", quiet = TRUE)
} else {
  cat("Seurat already installed\n")
}

# SeuratData (required for SlideSeq class)
if (!requireNamespace("SeuratData", quietly = TRUE)) {
  cat("Installing SeuratData from GitHub...\n
")
  remotes::install_github("satijalab/seurat-data", quiet = TRUE)
} else {
  cat("SeuratData already installed\n")
}

# Verify installation
cat("\n=== Verifying installation ===\n")
required <- c("Seurat", "SeuratData", "tidyverse", "rio", "janitor",
              "openxlsx", "data.table", "ggpubr", "gridExtra", "raster")

all_ok <- TRUE
for (pkg in required) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  ✓ %s\n", pkg))
  } else {
    cat(sprintf("  ✗ %s - FAILED\n", pkg))
    all_ok <- FALSE
  }
}

if (all_ok) {
  cat("\n=== All packages installed successfully! ===\n")
  cat("You can now run: Rscript Fig6-DBiT-run.R\n")
} else {
  cat("\n=== Some packages failed to install ===\n")
  cat("Please install missing packages manually.\n")
}
