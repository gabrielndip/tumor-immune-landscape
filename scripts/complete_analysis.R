# Portfolio 1: Complete Immune Cell Landscape Analysis
# Using 10x Genomics Public PBMC Dataset
# Focus: Demonstrating analytical methodology for tumor immune profiling

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

cat("=== Portfolio 1: Immune Cell Landscape Analysis ===\n")
cat("Dataset: 10x Genomics Public PBMC (Methodology Demonstration)\n\n")

# Paths
results_dir <- "/Volumes/MUKUTA/Projects/asgard_interview/portfolios/portfolio1_tumor_immune_landscape/results"
figures_dir <- "/Volumes/MUKUTA/Projects/asgard_interview/portfolios/portfolio1_tumor_immune_landscape/figures"

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# LOAD DATA - 10x Genomics 2.7k PBMCs
# ============================================================================
cat("STEP 1: Loading data...\n")

# Load pbmc dataset
pbmc.data <- Read10X(data.dir = "/Volumes/MUKUTA/Projects/asgard_interview/portfolios/portfolio1_tumor_immune_landscape/data/filtered_gene_bc_matrices/hg19/")

# Create Seurat object
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "ImmuneProfile", min.cells = 3, min.features = 200)

cat(sprintf("  ✓ Loaded %d cells, %d genes\n\n", ncol(pbmc), nrow(pbmc)))

# ============================================================================
# QC
# ============================================================================
cat("STEP 2: Quality Control...\n")

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

cat(sprintf("  Median genes/cell: %.0f\n", median(pbmc$nFeature_RNA)))
cat(sprintf("  Median UMIs/cell: %.0f\n", median(pbmc$nCount_RNA)))

# QC plots
p1 <- VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
ggsave(file.path(figures_dir, "01_qc_metrics.png"), plot = p1, width = 12, height = 4, dpi = 300)

p2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
p3 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
ggsave(file.path(figures_dir, "02_qc_scatter.png"), plot = p2 + p3, width = 12, height = 5, dpi = 300)

cat("  ✓ QC plots saved\n\n")

# Filter
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
cat(sprintf("  ✓ Filtered to %d cells\n\n", ncol(pbmc)))

# ============================================================================
# NORMALIZATION & SCALING
# ============================================================================
cat("STEP 3: Normalization and scaling...\n")

pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

top10 <- head(VariableFeatures(pbmc), 10)
cat(sprintf("  Top genes: %s\n", paste(top10, collapse = ", ")))

p_var <- VariableFeaturePlot(pbmc)
p_var <- LabelPoints(plot = p_var, points = top10, repel = TRUE)
ggsave(file.path(figures_dir, "03_variable_features.png"), plot = p_var, width = 10, height = 6, dpi = 300)

pbmc <- ScaleData(pbmc, features = rownames(pbmc))
cat("  ✓ Scaling complete\n\n")

# ============================================================================
# PCA & CLUSTERING
# ============================================================================
cat("STEP 4: PCA and clustering...\n")

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

p_elbow <- ElbowPlot(pbmc)
ggsave(file.path(figures_dir, "04_elbow_plot.png"), plot = p_elbow, width = 8, height = 6, dpi = 300)

pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:10)

cat(sprintf("  ✓ Identified %d clusters\n\n", length(unique(pbmc$seurat_clusters))))

# UMAP visualization
p_umap <- DimPlot(pbmc, reduction = "umap", label = TRUE) + NoLegend() +
          ggtitle("Immune Cell Clusters")
ggsave(file.path(figures_dir, "05_umap_clusters.png"), plot = p_umap, width = 10, height = 8, dpi = 300)
cat("  ✓ Saved KEY FIGURE 1: UMAP clusters\n\n")

# ============================================================================
# FIND MARKERS
# ============================================================================
cat("STEP 5: Finding marker genes...\n")

pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(pbmc.markers, file.path(results_dir, "cluster_markers.csv"), row.names = FALSE)

top5 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
cat(sprintf("  ✓ Found %d marker genes\n\n", nrow(pbmc.markers)))

# ============================================================================
# CELL TYPE ANNOTATION
# ============================================================================
cat("STEP 6: Cell type annotation...\n")

# Canonical markers
markers <- c("IL7R", "CCR7", "CD14", "LYZ", "MS4A1", "CD8A", "FCGR3A", "MS4A7", "GNLY", "NKG7", "FCER1A", "CST3", "PPBP")

p_markers <- FeaturePlot(pbmc, features = markers[1:9], ncol = 3)
ggsave(file.path(figures_dir, "06_immune_markers.png"), plot = p_markers, width = 15, height = 12, dpi = 300)
cat("  ✓ Saved KEY FIGURE 2: Immune markers\n\n")

# Annotate clusters
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono", "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
pbmc$celltype <- Idents(pbmc)

p_annotated <- DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(file.path(figures_dir, "07_umap_annotated.png"), plot = p_annotated, width = 10, height = 8, dpi = 300)
cat("  ✓ Saved KEY FIGURE 3: Annotated cell types\n\n")

# ============================================================================
# T CELL ANALYSIS
# ============================================================================
cat("STEP 7: T cell subset analysis...\n")

t_cells <- subset(pbmc, idents = c("Naive CD4 T", "Memory CD4 T", "CD8 T"))
cat(sprintf("  Extracted %d T cells\n", ncol(t_cells)))

# T cell markers
t_markers <- c("IL7R", "CCR7", "CD8A", "CD8B", "CD4", "SELL", "TCF7", "CD69", "CD44")
t_markers <- t_markers[t_markers %in% rownames(t_cells)]

p_tcell <- VlnPlot(t_cells, features = t_markers[1:6], ncol = 3, pt.size = 0.1)
ggsave(file.path(figures_dir, "08_tcell_markers.png"), plot = p_tcell, width = 14, height = 8, dpi = 300)
cat("  ✓ Saved KEY FIGURE 4: T cell analysis\n\n")

# ============================================================================
# SAVE RESULTS
# ============================================================================
cat("STEP 8: Saving results...\n")

saveRDS(pbmc, file.path(results_dir, "pbmc_analyzed.rds"))
cat("  ✓ Saved Seurat object\n\n")

# Summary stats
summary_stats <- data.frame(
  Metric = c("Total cells", "Genes", "Clusters", "Cell types", "T cells", "Myeloid cells", "B cells", "NK cells"),
  Value = c(
    ncol(pbmc),
    nrow(pbmc),
    length(unique(pbmc$seurat_clusters)),
    length(unique(pbmc$celltype)),
    sum(pbmc$celltype %in% c("Naive CD4 T", "Memory CD4 T", "CD8 T")),
    sum(pbmc$celltype %in% c("CD14+ Mono", "FCGR3A+ Mono", "DC")),
    sum(pbmc$celltype == "B"),
    sum(pbmc$celltype == "NK")
  )
)

write.csv(summary_stats, file.path(results_dir, "analysis_summary.csv"), row.names = FALSE)

cat("=== ANALYSIS COMPLETE ===\n\n")
cat("Summary:\n")
print(summary_stats)
cat("\nKey Figures Generated:\n")
cat("  1. 05_umap_clusters.png - Main clustering\n")
cat("  2. 06_immune_markers.png - Marker expression\n")
cat("  3. 07_umap_annotated.png - Annotated cell types\n")
cat("  4. 08_tcell_markers.png - T cell analysis\n")
cat("\nReady for presentation!\n")
