# Fig 6: DBiT-seq spatial transcriptomics analysis - YUMMER 2025 dataset
# 4 conditions: Untreated, CPIlow, CPIcd40, IL12IL18CD45
# Based on Fig6-DBiT-run.R

# Load configuration
source('config-2025.R')

# Set seed for reproducibility (UMAP, bootstrap sampling, etc.)
set.seed(42)

# Load required libraries
library(ggplot2)
library(plyr)
library(gridExtra)
library(magrittr)
library(tidyr)
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
library(data.table)

# Load helper functions
source('DBiT-func.R')

cat("Loading YUMMER 2025 DBiT-seq data...\n")

# Generate Seurat objects for each condition
# Untreated (control)
cat("Processing Untreated sample...\n")
untr.obj <- generate_Seurat_obj(
  file.path(DIR_2025, UNTR_SUBDIR),
  'YRuntr',
  UNTR_FILE,
  0.9,
  'YRuntr',
  pos_file = UNTR_POS
)

# CPI low dose
cat("Processing CPIlow sample...\n")
cpilow.obj <- generate_Seurat_obj(
  file.path(DIR_2025, CPILOW_SUBDIR),
  'YRcpilow',
  CPILOW_FILE,
  0.9,
  'YRcpilow',
  pos_file = CPILOW_POS
)

# CPI + CD40
cat("Processing CPIcd40 sample...\n")
cpicd40.obj <- generate_Seurat_obj(
  file.path(DIR_2025, CPICD40_SUBDIR),
  'YRcpicd40',
  CPICD40_FILE,
  0.9,
  'YRcpicd40',
  pos_file = CPICD40_POS
)

# IL12/IL18/CD45
cat("Processing IL12IL18CD45 sample...\n")
il12.obj <- generate_Seurat_obj(
  file.path(DIR_2025, IL12_SUBDIR),
  'YRil12',
  IL12_FILE,
  0.9,
  'YRil12',
  pos_file = IL12_POS
)

# ============================================================================
# Figure 6A: Visualize data by cluster
# ============================================================================
cat("\nGenerating Figure 6A: Spatial cluster plots...\n")

pdf(file.path(OUTPUT_DIR, "Fig6A_spatial_clusters.pdf"), width = 16, height = 4)
par(mfrow = c(1, 4))

p1 <- SpatialDimPlot(untr.obj, group.by = "seurat_clusters", pt.size.factor = 5) +
  scale_x_continuous(name = "X", expand = expansion(mult = c(0.008, 0.008))) +
  scale_y_continuous(name = "Y", expand = expansion(mult = c(0.008, 0.008))) +
  theme(legend.position = "right") + NoAxes() + ggtitle("Untreated")

p2 <- SpatialDimPlot(cpilow.obj, group.by = "seurat_clusters", pt.size.factor = 5) +
  scale_x_continuous(name = "X", expand = expansion(mult = c(0.008, 0.008))) +
  scale_y_continuous(name = "Y", expand = expansion(mult = c(0.008, 0.008))) +
  theme(legend.position = "right") + NoAxes() + ggtitle("CPI Low")

p3 <- SpatialDimPlot(cpicd40.obj, group.by = "seurat_clusters", pt.size.factor = 5) +
  scale_x_continuous(name = "X", expand = expansion(mult = c(0.008, 0.008))) +
  scale_y_continuous(name = "Y", expand = expansion(mult = c(0.008, 0.008))) +
  theme(legend.position = "right") + NoAxes() + ggtitle("CPI + CD40")

p4 <- SpatialDimPlot(il12.obj, group.by = "seurat_clusters", pt.size.factor = 5) +
  scale_x_continuous(name = "X", expand = expansion(mult = c(0.008, 0.008))) +
  scale_y_continuous(name = "Y", expand = expansion(mult = c(0.008, 0.008))) +
  theme(legend.position = "right") + NoAxes() + ggtitle("IL12/IL18/CD45")

print(gridExtra::grid.arrange(p1, p2, p3, p4, ncol = 4))
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
    if (length(which(rownames(untr.obj@assays$SCT) == j)) > 0 &
        length(which(rownames(cpilow.obj@assays$SCT) == j)) > 0 &
        length(which(rownames(cpicd40.obj@assays$SCT) == j)) > 0 &
        length(which(rownames(il12.obj@assays$SCT) == j)) > 0) {
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

untr.obj <- AddModuleScore(object = untr.obj, features = sigs.detected, ctrl = 100, assay = 'SCT', nbin = 10)
cpilow.obj <- AddModuleScore(object = cpilow.obj, features = sigs.detected, ctrl = 100, assay = 'SCT', nbin = 10)
cpicd40.obj <- AddModuleScore(object = cpicd40.obj, features = sigs.detected, ctrl = 100, assay = 'SCT', nbin = 10)
il12.obj <- AddModuleScore(object = il12.obj, features = sigs.detected, ctrl = 100, assay = 'SCT', nbin = 10)

# ============================================================================
# Visualize cell type signature scores
# Cluster mapping: Macrophage (Cluster1), Treg (Cluster4), CD8 (Cluster3)
# ============================================================================
cat("\nGenerating cell type signature score plots...\n")

pdf(file.path(OUTPUT_DIR, "Fig6_celltype_signatures.pdf"), width = 16, height = 12)

for (k in c('Cluster1', 'Cluster2', 'Cluster3', 'Cluster4')) {
  celltype <- switch(k,
    'Cluster1' = 'Macrophage',
    'Cluster2' = 'mregDC',
    'Cluster3' = 'CD8+ T cell',
    'Cluster4' = 'Treg'
  )

  p1 <- SpatialFeaturePlot(untr.obj, features = k, pt.size.factor = 5, min.cutoff = 0, max.cutoff = 0.2) +
    scale_x_continuous(name = "X", expand = expansion(mult = c(0.008, 0.008))) +
    scale_y_continuous(name = "Y", expand = expansion(mult = c(0.008, 0.008))) +
    theme(legend.position = "right") + NoAxes() + ggtitle(paste(celltype, "- Untreated"))

  p2 <- SpatialFeaturePlot(cpilow.obj, features = k, pt.size.factor = 5, min.cutoff = 0, max.cutoff = 0.2) +
    scale_x_continuous(name = "X", expand = expansion(mult = c(0.008, 0.008))) +
    scale_y_continuous(name = "Y", expand = expansion(mult = c(0.008, 0.008))) +
    theme(legend.position = "right") + NoAxes() + ggtitle(paste(celltype, "- CPI Low"))

  p3 <- SpatialFeaturePlot(cpicd40.obj, features = k, pt.size.factor = 5, min.cutoff = 0, max.cutoff = 0.2) +
    scale_x_continuous(name = "X", expand = expansion(mult = c(0.008, 0.008))) +
    scale_y_continuous(name = "Y", expand = expansion(mult = c(0.008, 0.008))) +
    theme(legend.position = "right") + NoAxes() + ggtitle(paste(celltype, "- CPI+CD40"))

  p4 <- SpatialFeaturePlot(il12.obj, features = k, pt.size.factor = 5, min.cutoff = 0, max.cutoff = 0.2) +
    scale_x_continuous(name = "X", expand = expansion(mult = c(0.008, 0.008))) +
    scale_y_continuous(name = "Y", expand = expansion(mult = c(0.008, 0.008))) +
    theme(legend.position = "right") + NoAxes() + ggtitle(paste(celltype, "- IL12/IL18/CD45"))

  print(gridExtra::grid.arrange(p1, p2, p3, p4, ncol = 4))
}

dev.off()

# ============================================================================
# Correlation analysis: Cell type spatial co-localization
# ============================================================================
cat("\nRunning correlation analysis...\n")

# Macrophage-Treg correlation (20x20 grid, 100 bootstrap reps)
cat("  Computing Mac-Treg correlations...\n")
untr.mac.treg <- corr_CV(untr.obj, 'Cluster1', 'Cluster4', 20, rep = 100)
cpilow.mac.treg <- corr_CV(cpilow.obj, 'Cluster1', 'Cluster4', 20, rep = 100)
cpicd40.mac.treg <- corr_CV(cpicd40.obj, 'Cluster1', 'Cluster4', 20, rep = 100)
il12.mac.treg <- corr_CV(il12.obj, 'Cluster1', 'Cluster4', 20, rep = 100)
mac.treg <- tibble('untr' = untr.mac.treg$corr,
                   'cpilow' = cpilow.mac.treg$corr,
                   'cpicd40' = cpicd40.mac.treg$corr,
                   'il12' = il12.mac.treg$corr)
openxlsx::write.xlsx(mac.treg, file.path(OUTPUT_DIR, "mac-Treg_corr.xlsx"))

# Macrophage-CD8 T cell correlation (20x20 grid, 1000 bootstrap reps)
cat("  Computing Mac-CD8 correlations...\n")
untr.mac.cd8 <- corr_CV(untr.obj, 'Cluster1', 'Cluster3', 20, rep = 1000)
cpilow.mac.cd8 <- corr_CV(cpilow.obj, 'Cluster1', 'Cluster3', 20, rep = 1000)
cpicd40.mac.cd8 <- corr_CV(cpicd40.obj, 'Cluster1', 'Cluster3', 20, rep = 1000)
il12.mac.cd8 <- corr_CV(il12.obj, 'Cluster1', 'Cluster3', 20, rep = 1000)
mac.cd8 <- tibble('untr' = untr.mac.cd8$corr,
                  'cpilow' = cpilow.mac.cd8$corr,
                  'cpicd40' = cpicd40.mac.cd8$corr,
                  'il12' = il12.mac.cd8$corr)
openxlsx::write.xlsx(mac.cd8, file.path(OUTPUT_DIR, "mac-CD8_corr.xlsx"))

# CD8-Treg correlation (20x20 grid, 100 bootstrap reps)
cat("  Computing CD8-Treg correlations...\n")
untr.cd8.treg <- corr_CV(untr.obj, 'Cluster3', 'Cluster4', 20, rep = 100)
cpilow.cd8.treg <- corr_CV(cpilow.obj, 'Cluster3', 'Cluster4', 20, rep = 100)
cpicd40.cd8.treg <- corr_CV(cpicd40.obj, 'Cluster3', 'Cluster4', 20, rep = 100)
il12.cd8.treg <- corr_CV(il12.obj, 'Cluster3', 'Cluster4', 20, rep = 100)
cd8.treg <- tibble('untr' = untr.cd8.treg$corr,
                   'cpilow' = cpilow.cd8.treg$corr,
                   'cpicd40' = cpicd40.cd8.treg$corr,
                   'il12' = il12.cd8.treg$corr)
openxlsx::write.xlsx(cd8.treg, file.path(OUTPUT_DIR, "CD8-Treg_corr.xlsx"))

# mregDC-Macrophage correlation (20x20 grid, 100 bootstrap reps)
cat("  Computing mregDC-Mac correlations...\n")
untr.mregdc.mac <- corr_CV(untr.obj, 'Cluster2', 'Cluster1', 20, rep = 100)
cpilow.mregdc.mac <- corr_CV(cpilow.obj, 'Cluster2', 'Cluster1', 20, rep = 100)
cpicd40.mregdc.mac <- corr_CV(cpicd40.obj, 'Cluster2', 'Cluster1', 20, rep = 100)
il12.mregdc.mac <- corr_CV(il12.obj, 'Cluster2', 'Cluster1', 20, rep = 100)
mregdc.mac <- tibble('untr' = untr.mregdc.mac$corr,
                     'cpilow' = cpilow.mregdc.mac$corr,
                     'cpicd40' = cpicd40.mregdc.mac$corr,
                     'il12' = il12.mregdc.mac$corr)
openxlsx::write.xlsx(mregdc.mac, file.path(OUTPUT_DIR, "mregDC-Mac_corr.xlsx"))

# mregDC-CD8 correlation (20x20 grid, 100 bootstrap reps)
cat("  Computing mregDC-CD8 correlations...\n")
untr.mregdc.cd8 <- corr_CV(untr.obj, 'Cluster2', 'Cluster3', 20, rep = 100)
cpilow.mregdc.cd8 <- corr_CV(cpilow.obj, 'Cluster2', 'Cluster3', 20, rep = 100)
cpicd40.mregdc.cd8 <- corr_CV(cpicd40.obj, 'Cluster2', 'Cluster3', 20, rep = 100)
il12.mregdc.cd8 <- corr_CV(il12.obj, 'Cluster2', 'Cluster3', 20, rep = 100)
mregdc.cd8 <- tibble('untr' = untr.mregdc.cd8$corr,
                     'cpilow' = cpilow.mregdc.cd8$corr,
                     'cpicd40' = cpicd40.mregdc.cd8$corr,
                     'il12' = il12.mregdc.cd8$corr)
openxlsx::write.xlsx(mregdc.cd8, file.path(OUTPUT_DIR, "mregDC-CD8_corr.xlsx"))

# mregDC-Treg correlation (20x20 grid, 100 bootstrap reps)
cat("  Computing mregDC-Treg correlations...\n")
untr.mregdc.treg <- corr_CV(untr.obj, 'Cluster2', 'Cluster4', 20, rep = 100)
cpilow.mregdc.treg <- corr_CV(cpilow.obj, 'Cluster2', 'Cluster4', 20, rep = 100)
cpicd40.mregdc.treg <- corr_CV(cpicd40.obj, 'Cluster2', 'Cluster4', 20, rep = 100)
il12.mregdc.treg <- corr_CV(il12.obj, 'Cluster2', 'Cluster4', 20, rep = 100)
mregdc.treg <- tibble('untr' = untr.mregdc.treg$corr,
                      'cpilow' = cpilow.mregdc.treg$corr,
                      'cpicd40' = cpicd40.mregdc.treg$corr,
                      'il12' = il12.mregdc.treg$corr)
openxlsx::write.xlsx(mregdc.treg, file.path(OUTPUT_DIR, "mregDC-Treg_corr.xlsx"))

# ============================================================================
# Three-feature correlation: Mac-CD8 vs Mac-Treg (all conditions)
# ============================================================================
cat("  Computing 3-feature correlations (untr)...\n")
untr.mac.cd8.treg <- corr_CV(untr.obj, 'Cluster1', 'Cluster3', 20, rep = 1000, feature3 = 'Cluster4')
cat(sprintf("    untr:    r = %.3f, p = %.2e\n",
            cor.test(untr.mac.cd8.treg$corr, untr.mac.cd8.treg$corr.feature3)$estimate,
            cor.test(untr.mac.cd8.treg$corr, untr.mac.cd8.treg$corr.feature3)$p.value))
openxlsx::write.xlsx(untr.mac.cd8.treg, file.path(OUTPUT_DIR, "3-feature_corr_untr.xlsx"))

cat("  Computing 3-feature correlations (cpilow)...\n")
cpilow.mac.cd8.treg <- corr_CV(cpilow.obj, 'Cluster1', 'Cluster3', 20, rep = 1000, feature3 = 'Cluster4')
cat(sprintf("    cpilow:  r = %.3f, p = %.2e\n",
            cor.test(cpilow.mac.cd8.treg$corr, cpilow.mac.cd8.treg$corr.feature3)$estimate,
            cor.test(cpilow.mac.cd8.treg$corr, cpilow.mac.cd8.treg$corr.feature3)$p.value))
openxlsx::write.xlsx(cpilow.mac.cd8.treg, file.path(OUTPUT_DIR, "3-feature_corr_cpilow.xlsx"))

cat("  Computing 3-feature correlations (cpicd40)...\n")
cpicd40.mac.cd8.treg <- corr_CV(cpicd40.obj, 'Cluster1', 'Cluster3', 20, rep = 1000, feature3 = 'Cluster4')
cat(sprintf("    cpicd40: r = %.3f, p = %.2e\n",
            cor.test(cpicd40.mac.cd8.treg$corr, cpicd40.mac.cd8.treg$corr.feature3)$estimate,
            cor.test(cpicd40.mac.cd8.treg$corr, cpicd40.mac.cd8.treg$corr.feature3)$p.value))
openxlsx::write.xlsx(cpicd40.mac.cd8.treg, file.path(OUTPUT_DIR, "3-feature_corr_cpicd40.xlsx"))

cat("  Computing 3-feature correlations (il12)...\n")
il12.mac.cd8.treg <- corr_CV(il12.obj, 'Cluster1', 'Cluster3', 20, rep = 1000, feature3 = 'Cluster4')
cat(sprintf("    il12:    r = %.3f, p = %.2e\n",
            cor.test(il12.mac.cd8.treg$corr, il12.mac.cd8.treg$corr.feature3)$estimate,
            cor.test(il12.mac.cd8.treg$corr, il12.mac.cd8.treg$corr.feature3)$p.value))
openxlsx::write.xlsx(il12.mac.cd8.treg, file.path(OUTPUT_DIR, "3-feature_corr_il12.xlsx"))

# Plot Mac-CD8 vs Mac-Treg correlation (all conditions)
pdf(file.path(OUTPUT_DIR, "Fig6_mac_cd8_vs_treg_corr.pdf"), width = 12, height = 6)

p1 <- untr.mac.cd8.treg %>%
  ggplot(aes(x = corr, y = corr.feature3)) +
  geom_point(alpha = 0.5) +
  theme_bw() +
  theme(aspect.ratio = 1) +
  labs(x = 'Mac-CD8+ correlation', y = 'Mac-Treg correlation',
       title = 'Untreated') +
  geom_smooth(method = "lm", color = "red")

p2 <- cpilow.mac.cd8.treg %>%
  ggplot(aes(x = corr, y = corr.feature3)) +
  geom_point(alpha = 0.5) +
  theme_bw() +
  theme(aspect.ratio = 1) +
  labs(x = 'Mac-CD8+ correlation', y = 'Mac-Treg correlation',
       title = 'CPI Low') +
  geom_smooth(method = "lm", color = "red")

p3 <- cpicd40.mac.cd8.treg %>%
  ggplot(aes(x = corr, y = corr.feature3)) +
  geom_point(alpha = 0.5) +
  theme_bw() +
  theme(aspect.ratio = 1) +
  labs(x = 'Mac-CD8+ correlation', y = 'Mac-Treg correlation',
       title = 'CPI + CD40') +
  geom_smooth(method = "lm", color = "red")

p4 <- il12.mac.cd8.treg %>%
  ggplot(aes(x = corr, y = corr.feature3)) +
  geom_point(alpha = 0.5) +
  theme_bw() +
  theme(aspect.ratio = 1) +
  labs(x = 'Mac-CD8+ correlation', y = 'Mac-Treg correlation',
       title = 'IL12/IL18/CD45') +
  geom_smooth(method = "lm", color = "red")

print(gridExtra::grid.arrange(p1, p2, p3, p4, ncol = 4))
dev.off()

cat("\n=== YUMMER 2025 Analysis complete! ===\n")
cat(sprintf("Output files saved to: %s\n", OUTPUT_DIR))
