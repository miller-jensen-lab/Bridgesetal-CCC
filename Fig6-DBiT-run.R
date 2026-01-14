# Fig 6: DBiT-seq spatial transcriptomics analysis
# Modified to use config.R for paths
# Original author: Kate Bridges

# Load configuration
source('config.R')

# Set seed for reproducibility (UMAP, bootstrap sampling, etc.)
set.seed(42)

# Load required libraries
library(ggplot2)
library(plyr)
library(gridExtra)
library(magrittr)
library(tidyr)
# library(raster)  # Not used in this analysis, requires terra/GDAL
# library(OpenImageR)  # Not used in this analysis
library(ggpubr)
library(grid)
library(wesanderson)
library(RColorBrewer)
library(dplyr)
library(tidyverse)
library(Seurat)
library(SeuratData)
library(Matrix)
library(janitor)
library(rio)
library(readxl)
library(data.table)  # Required for fread() in DBiT-func.R
# library(MERINGUE)  # Uncomment if needed for spatial analysis

# Load helper functions
source('DBiT-func.R')

cat("Loading DBiT-seq data...\n")

# Generate Seurat objects for each condition
# Control (from second batch)
cat("Processing Control sample...\n")
ctrl.obj <- generate_Seurat_obj(
  DIR_SECOND_BATCH,
  'YRctrl',
  CTRL_FILE,
  0.9,
  'YRctrl'
)

# ICB only (from third batch)
cat("Processing ICB sample...\n")
icb.obj <- generate_Seurat_obj(
  file.path(DIR_THIRD_BATCH, ICB_SUBDIR),
  'YRicb',
  ICB_FILE,
  0.9,
  'YRicb',
  sept = ',',
  transp = FALSE
)

# ICB + CD40ag combination (from third batch)
cat("Processing ICB+CD40ag sample...\n")
comb.obj <- generate_Seurat_obj(
  file.path(DIR_THIRD_BATCH, COMB_SUBDIR),
  'YRcomb',
  COMB_FILE,
  0.9,
  'YRcomb',
  sept = ',',
  transp = FALSE
)

# ============================================================================
# Figure 6A: Visualize data by cluster
# ============================================================================
cat("\nGenerating Figure 6A: Spatial cluster plots...\n")

pdf(file.path(OUTPUT_DIR, "Fig6A_spatial_clusters.pdf"), width = 12, height = 4)
par(mfrow = c(1, 3))

p1 <- SpatialDimPlot(ctrl.obj, group.by = "seurat_clusters", pt.size.factor = 5) +
  scale_x_continuous(name = "X", expand = expansion(mult = c(0.008, 0.008))) +
  scale_y_continuous(name = "Y", expand = expansion(mult = c(0.008, 0.008))) +
  theme(legend.position = "right") + NoAxes() + ggtitle("Control")

p2 <- SpatialDimPlot(icb.obj, group.by = "seurat_clusters", pt.size.factor = 5) +
  scale_x_continuous(name = "X", expand = expansion(mult = c(0.008, 0.008))) +
  scale_y_continuous(name = "Y", expand = expansion(mult = c(0.008, 0.008))) +
  theme(legend.position = "right") + NoAxes() + ggtitle("ICB")

p3 <- SpatialDimPlot(comb.obj, group.by = "seurat_clusters", pt.size.factor = 5) +
  scale_x_continuous(name = "X", expand = expansion(mult = c(0.008, 0.008))) +
  scale_y_continuous(name = "Y", expand = expansion(mult = c(0.008, 0.008))) +
  theme(legend.position = "right") + NoAxes() + ggtitle("ICB + CD40ag")

print(gridExtra::grid.arrange(p1, p2, p3, ncol = 3))
dev.off()

# ============================================================================
# Load cell type signature genes from scRNA-seq DEGs
# ============================================================================
cat("\nLoading cell type signature genes...\n")
sigs <- import_list(SIGNATURE_GENES_FILE, header = FALSE)

# Limit signature genes to ones consistently detected across all datasets
sigs.detected <- c()

for (g in names(sigs)) {
  top20 <- c()
  i <- 1
  while (length(top20) < 20 && i <= nrow(sigs[[g]])) {
    j <- sigs[[g]]$...1[i]
    if (length(which(rownames(ctrl.obj@assays$SCT) == j)) > 0 &
        length(which(rownames(icb.obj@assays$SCT) == j)) > 0 &
        length(which(rownames(comb.obj@assays$SCT) == j)) > 0) {
      top20 <- c(top20, j)
    }
    i <- i + 1
  }
  sigs.detected[[g]] <- top20
  cat(sprintf("  %s: %d genes\n", g, length(top20)))
}

# Save detected signatures
openxlsx::write.xlsx(sigs.detected, file.path(OUTPUT_DIR, "sigsdetected-DBiT.xlsx"))

# ============================================================================
# Score pixels by cell type signatures
# ============================================================================
cat("\nScoring pixels by cell type signatures...\n")

ctrl.obj <- AddModuleScore(object = ctrl.obj, features = sigs.detected, ctrl = 100, assay = 'SCT', nbin = 10)
icb.obj <- AddModuleScore(object = icb.obj, features = sigs.detected, ctrl = 100, assay = 'SCT', nbin = 10)
comb.obj <- AddModuleScore(object = comb.obj, features = sigs.detected, ctrl = 100, assay = 'SCT', nbin = 10)

# ============================================================================
# Visualize cell type signature scores
# Cluster mapping: Macrophage (Cluster1), Treg (Cluster4), CD8 (Cluster3)
# ============================================================================
cat("\nGenerating cell type signature score plots...\n")

pdf(file.path(OUTPUT_DIR, "Fig6_celltype_signatures.pdf"), width = 12, height = 12)

for (k in c('Cluster1', 'Cluster4', 'Cluster3')) {
  celltype <- switch(k,
    'Cluster1' = 'Macrophage',
    'Cluster3' = 'CD8+ T cell',
    'Cluster4' = 'Treg'
  )

  p1 <- SpatialFeaturePlot(ctrl.obj, features = k, pt.size.factor = 5, min.cutoff = 0, max.cutoff = 0.2) +
    scale_x_continuous(name = "X", expand = expansion(mult = c(0.008, 0.008))) +
    scale_y_continuous(name = "Y", expand = expansion(mult = c(0.008, 0.008))) +
    theme(legend.position = "right") + NoAxes() + ggtitle(paste(celltype, "- Control"))

  p2 <- SpatialFeaturePlot(icb.obj, features = k, pt.size.factor = 5, min.cutoff = 0, max.cutoff = 0.2) +
    scale_x_continuous(name = "X", expand = expansion(mult = c(0.008, 0.008))) +
    scale_y_continuous(name = "Y", expand = expansion(mult = c(0.008, 0.008))) +
    theme(legend.position = "right") + NoAxes() + ggtitle(paste(celltype, "- ICB"))

  p3 <- SpatialFeaturePlot(comb.obj, features = k, pt.size.factor = 5, min.cutoff = 0, max.cutoff = 0.2) +
    scale_x_continuous(name = "X", expand = expansion(mult = c(0.008, 0.008))) +
    scale_y_continuous(name = "Y", expand = expansion(mult = c(0.008, 0.008))) +
    theme(legend.position = "right") + NoAxes() + ggtitle(paste(celltype, "- ICB+CD40ag"))

  print(gridExtra::grid.arrange(p1, p2, p3, ncol = 3))
}

dev.off()

# ============================================================================
# Correlation analysis: Cell type spatial co-localization
# ============================================================================
cat("\nRunning correlation analysis...\n")

# Macrophage-Treg correlation (20x20 grid, 100 bootstrap reps)
cat("  Computing Mac-Treg correlations...\n")
ctrl.mac.treg <- corr_CV(ctrl.obj, 'Cluster1', 'Cluster4', 20, rep = 100)
icb.mac.treg <- corr_CV(icb.obj, 'Cluster1', 'Cluster4', 20, rep = 100)
comb.mac.treg <- corr_CV(comb.obj, 'Cluster1', 'Cluster4', 20, rep = 100)
mac.treg <- tibble('ctrl' = ctrl.mac.treg$corr,
                   'icb' = icb.mac.treg$corr,
                   'comb' = comb.mac.treg$corr)
openxlsx::write.xlsx(mac.treg, file.path(OUTPUT_DIR, "mac-Treg_corr.xlsx"))

# Macrophage-CD8 T cell correlation (20x20 grid, 1000 bootstrap reps)
cat("  Computing Mac-CD8 correlations...\n")
ctrl.mac.cd8 <- corr_CV(ctrl.obj, 'Cluster1', 'Cluster3', 20, rep = 1000)
icb.mac.cd8 <- corr_CV(icb.obj, 'Cluster1', 'Cluster3', 20, rep = 1000)
comb.mac.cd8 <- corr_CV(comb.obj, 'Cluster1', 'Cluster3', 20, rep = 1000)
mac.cd8 <- tibble('ctrl' = ctrl.mac.cd8$corr,
                  'icb' = icb.mac.cd8$corr,
                  'comb' = comb.mac.cd8$corr)
openxlsx::write.xlsx(mac.cd8, file.path(OUTPUT_DIR, "mac-CD8_corr.xlsx"))

# CD8-Treg correlation (20x20 grid, 100 bootstrap reps)
cat("  Computing CD8-Treg correlations...\n")
ctrl.cd8.treg <- corr_CV(ctrl.obj, 'Cluster3', 'Cluster4', 20, rep = 100)
icb.cd8.treg <- corr_CV(icb.obj, 'Cluster3', 'Cluster4', 20, rep = 100)
comb.cd8.treg <- corr_CV(comb.obj, 'Cluster3', 'Cluster4', 20, rep = 100)
cd8.treg <- tibble('ctrl' = ctrl.cd8.treg$corr,
                   'icb' = icb.cd8.treg$corr,
                   'comb' = comb.cd8.treg$corr)
openxlsx::write.xlsx(cd8.treg, file.path(OUTPUT_DIR, "CD8-Treg_corr.xlsx"))

# ============================================================================
# Three-feature correlation: Mac-CD8 vs Mac-Treg (all conditions)
# ============================================================================
cat("  Computing 3-feature correlations (ctrl)...\n")
ctrl.mac.cd8.treg <- corr_CV(ctrl.obj, 'Cluster1', 'Cluster3', 20, rep = 1000, feature3 = 'Cluster4')
cat(sprintf("    ctrl: r = %.3f, p = %.2e\n",
            cor.test(ctrl.mac.cd8.treg$corr, ctrl.mac.cd8.treg$corr.feature3)$estimate,
            cor.test(ctrl.mac.cd8.treg$corr, ctrl.mac.cd8.treg$corr.feature3)$p.value))
openxlsx::write.xlsx(ctrl.mac.cd8.treg, file.path(OUTPUT_DIR, "3-feature_corr_ctrl.xlsx"))

cat("  Computing 3-feature correlations (icb)...\n")
icb.mac.cd8.treg <- corr_CV(icb.obj, 'Cluster1', 'Cluster3', 20, rep = 1000, feature3 = 'Cluster4')
cat(sprintf("    icb:  r = %.3f, p = %.2e\n",
            cor.test(icb.mac.cd8.treg$corr, icb.mac.cd8.treg$corr.feature3)$estimate,
            cor.test(icb.mac.cd8.treg$corr, icb.mac.cd8.treg$corr.feature3)$p.value))
openxlsx::write.xlsx(icb.mac.cd8.treg, file.path(OUTPUT_DIR, "3-feature_corr_icb.xlsx"))

cat("  Computing 3-feature correlations (comb)...\n")
comb.mac.cd8.treg <- corr_CV(comb.obj, 'Cluster1', 'Cluster3', 20, rep = 1000, feature3 = 'Cluster4')
cat(sprintf("    comb: r = %.3f, p = %.2e\n",
            cor.test(comb.mac.cd8.treg$corr, comb.mac.cd8.treg$corr.feature3)$estimate,
            cor.test(comb.mac.cd8.treg$corr, comb.mac.cd8.treg$corr.feature3)$p.value))
openxlsx::write.xlsx(comb.mac.cd8.treg, file.path(OUTPUT_DIR, "3-feature_corr_comb.xlsx"))

# Plot Mac-CD8 vs Mac-Treg correlation (comb only, for backward compatibility)
pdf(file.path(OUTPUT_DIR, "Fig6_mac_cd8_vs_treg_corr.pdf"), width = 6, height = 6)
print(
  comb.mac.cd8.treg %>%
    ggplot(aes(x = corr, y = corr.feature3)) +
    geom_point(alpha = 0.5) +
    theme_bw() +
    theme(aspect.ratio = 1) +
    labs(x = 'Mac-CD8+ correlation', y = 'Mac-Treg correlation',
         title = 'ICB + CD40ag: Spatial correlation of cell type signatures') +
    geom_smooth(method = "lm", color = "red")
)
dev.off()

cat("\n=== Analysis complete! ===\n")
cat(sprintf("Output files saved to: %s\n", OUTPUT_DIR))
