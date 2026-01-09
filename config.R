# Configuration file for Bridges et al. DBiT-seq analysis
# Modify these paths to match your local data directory structure

# Base directory for data (relative to project root)
DATA_DIR <- "tmp"

# DBiT-seq data directories
DIR_THIRD_BATCH <- file.path(DATA_DIR, "YUMMER_ThirdBatch_Treatments")
DIR_SECOND_BATCH <- file.path(DATA_DIR, "YUMMER_SecondBatch ")  # Note: trailing space in original dir name

# Input files for each condition
# Control (from second batch)
CTRL_FILE <- "YUM_gex_secondbatch.updated.tsv"  # Original expected: YUM_gex_secondbatch_updated.tsv

# ICB only (from third batch)
ICB_SUBDIR <- "CPI_0610"
ICB_FILE <- "CPI_0610_repaired.csv"

# ICB + CD40ag combination (from third batch)
COMB_SUBDIR <- "CPICD40_0610"
COMB_FILE <- "CPICD40_0610.repaired.csv"  # Original expected: CPICD40_0610_repaired.csv

# Cell type signature genes
SIGNATURE_GENES_FILE <- file.path(DATA_DIR, "celltype-signature-genes.xlsx")

# Output directory
OUTPUT_DIR <- file.path(DATA_DIR, "output")

# Create output directory if it doesn't exist
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}
