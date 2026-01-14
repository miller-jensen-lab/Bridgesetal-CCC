# Configuration for YUMMER 2025 DBiT-seq analysis
# New dataset with 4 conditions: Untreated, CPIlow, CPIcd40, IL12IL18CD45

# Base directory for data (relative to project root)
DATA_DIR <- "tmp"
DIR_2025 <- file.path(DATA_DIR, "YUMMER_2025")

# Untreated (control)
UNTR_SUBDIR <- "Untreated"
UNTR_FILE <- "YUM_untr_0616.updated.tsv"
UNTR_POS <- "positionuntr.txt"

# CPI low dose
CPILOW_SUBDIR <- "CPIlow"
CPILOW_FILE <- "YUM_cpilow_0616.updated.tsv"
CPILOW_POS <- "positioncpilow.txt"

# CPI + CD40
CPICD40_SUBDIR <- "CPIcd40"
CPICD40_FILE <- "YUM_CPICD40_0610_updated.tsv"
CPICD40_POS <- "positionCPICD40.txt"

# IL12/IL18/CD45
IL12_SUBDIR <- "IL12IL18CD45"
IL12_FILE <- "YUM_A_0603.updated.tsv"
IL12_POS <- "positionil12il18cd45.txt"

# Cell type signature genes (shared with original analysis)
SIGNATURE_GENES_FILE <- file.path(DATA_DIR, "celltype-signature-genes.xlsx")

# Output directory
OUTPUT_DIR <- file.path(DIR_2025, "output")

# Create output directory if it doesn't exist
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}
