# Portfolio 1: Clustering and Cell Type Identification
# Dimensionality reduction, clustering, and immune cell annotation

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

cat("=== Portfolio 1: Clustering and Cell Type Identification ===\n\n")

# Set paths
results_dir <- "/Volumes/MUKUTA/Projects/asgard_interview/portfolios/portfolio1_tumor_immune_landscape/results"
figures_dir <- "/Volumes/MUKUTA/Projects/asgard_interview/portfolios/portfolio1_tumor_immune_landscape/figures"

# Load QC'd Seurat object
cat("Loading Seurat object...\n")
seurat_obj <- readRDS(file.path(results_dir, "seurat_qc.rds"))
cat(sprintf("  ✓ Loaded object with %d cells\n", ncol(seurat_obj)))

# Filter cells based on QC metrics
cat("\nFiltering cells based on QC thresholds...\n")
seurat_obj <- subset(seurat_obj,
                     subset = nFeature_RNA > 200 &
                             nFeature_RNA < 5000 &
                             percent.mt < 20)

cat(sprintf("  ✓ Retained %d cells after filtering\n", ncol(seurat_obj)))

# Normalization
cat("\nNormalizing data...\n")
seurat_obj <- NormalizeData(seurat_obj,
                           normalization.method = "LogNormalize",
                           scale.factor = 10000)
cat("  ✓ Normalization complete\n")

# Find variable features
cat("\nIdentifying highly variable features...\n")
seurat_obj <- FindVariableFeatures(seurat_obj,
                                  selection.method = "vst",
                                  nfeatures = 2000)

cat(sprintf("  ✓ Identified %d variable features\n",
            length(VariableFeatures(seurat_obj))))

# Top 10 variable genes
top10 <- head(VariableFeatures(seurat_obj), 10)
cat(sprintf("  Top 10: %s\n", paste(top10, collapse = ", ")))

# Plot variable features
p1 <- VariableFeaturePlot(seurat_obj)
p2 <- LabelPoints(plot = p1, points = top10, repel = TRUE)
ggsave(file.path(figures_dir, "03_variable_features.png"),
       plot = p2, width = 10, height = 6, dpi = 300)
cat("  ✓ Saved: 03_variable_features.png\n")

# Scale data
cat("\nScaling data...\n")
all_genes <- rownames(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, features = all_genes)
cat("  ✓ Data scaling complete\n")

# PCA
cat("\nRunning PCA...\n")
seurat_obj <- RunPCA(seurat_obj,
                    features = VariableFeatures(object = seurat_obj),
                    npcs = 50,
                    verbose = FALSE)
cat("  ✓ PCA complete\n")

# Visualize PCA
p_pca <- DimPlot(seurat_obj, reduction = "pca") +
         ggtitle("PCA of scRNA-seq Data")
ggsave(file.path(figures_dir, "04_pca.png"),
       plot = p_pca, width = 8, height = 6, dpi = 300)
cat("  ✓ Saved: 04_pca.png\n")

# Elbow plot
p_elbow <- ElbowPlot(seurat_obj, ndims = 50) +
           ggtitle("PCA Elbow Plot")
ggsave(file.path(figures_dir, "05_pca_elbow.png"),
       plot = p_elbow, width = 8, height = 6, dpi = 300)
cat("  ✓ Saved: 05_pca_elbow.png\n")

# Clustering
cat("\nPerforming clustering...\n")
n_dims <- 30  # Based on elbow plot
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:n_dims)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5, verbose = FALSE)

cat(sprintf("  ✓ Identified %d clusters\n",
            length(unique(seurat_obj$seurat_clusters))))

# UMAP
cat("\nComputing UMAP...\n")
seurat_obj <- RunUMAP(seurat_obj, dims = 1:n_dims, verbose = FALSE)
cat("  ✓ UMAP complete\n")

# Visualize UMAP
p_umap <- DimPlot(seurat_obj,
                  reduction = "umap",
                  label = TRUE,
                  label.size = 5) +
          ggtitle("UMAP: Cell Clusters") +
          NoLegend()

ggsave(file.path(figures_dir, "06_umap_clusters.png"),
       plot = p_umap, width = 10, height = 8, dpi = 300)
cat("  ✓ Saved: 06_umap_clusters.png (KEY FIGURE 1)\n")

# If annotations exist in metadata, plot them
if ("Cell_type" %in% colnames(seurat_obj@meta.data)) {
    p_umap_annot <- DimPlot(seurat_obj,
                           reduction = "umap",
                           group.by = "Cell_type",
                           label = TRUE,
                           repel = TRUE) +
                    ggtitle("UMAP: Annotated Cell Types")

    ggsave(file.path(figures_dir, "07_umap_cell_types.png"),
           plot = p_umap_annot, width = 12, height = 8, dpi = 300)
    cat("  ✓ Saved: 07_umap_cell_types.png\n")
}

# Find cluster markers
cat("\nFinding cluster markers (top genes per cluster)...\n")
cluster_markers <- FindAllMarkers(seurat_obj,
                                 only.pos = TRUE,
                                 min.pct = 0.25,
                                 logfc.threshold = 0.25,
                                 verbose = FALSE)

cat(sprintf("  ✓ Found markers for %d clusters\n",
            length(unique(cluster_markers$cluster))))

# Save marker genes
write.csv(cluster_markers,
          file.path(results_dir, "cluster_markers.csv"),
          row.names = FALSE)
cat("  ✓ Saved: cluster_markers.csv\n")

# Top markers per cluster
top5 <- cluster_markers %>%
        group_by(cluster) %>%
        top_n(n = 5, wt = avg_log2FC)

cat("\nTop 5 markers per cluster:\n")
print(top5 %>% select(cluster, gene, avg_log2FC, pct.1, pct.2))

# Visualize key immune cell markers
cat("\nVisualizing immune cell markers...\n")

immune_markers <- c(
    "Cd3e", "Cd3d",      # T cells
    "Cd8a", "Cd8b1",     # CD8+ T cells
    "Cd4",               # CD4+ T cells
    "Foxp3",             # Regulatory T cells
    "Ncr1", "Klrb1c",    # NK cells
    "Cd19", "Ms4a1",     # B cells
    "Lyz2", "Adgre1",    # Macrophages/Myeloid
    "H2-Aa", "H2-Ab1"    # DCs (MHC-II)
)

# Filter to markers present in dataset
immune_markers_present <- immune_markers[immune_markers %in% rownames(seurat_obj)]

if (length(immune_markers_present) > 0) {
    p_markers <- FeaturePlot(seurat_obj,
                            features = head(immune_markers_present, 9),
                            ncol = 3,
                            pt.size = 0.5)

    ggsave(file.path(figures_dir, "08_immune_markers_umap.png"),
           plot = p_markers, width = 15, height = 12, dpi = 300)
    cat("  ✓ Saved: 08_immune_markers_umap.png (KEY FIGURE 2)\n")
}

# Save processed object
cat("\nSaving clustered Seurat object...\n")
saveRDS(seurat_obj, file.path(results_dir, "seurat_clustered.rds"))
cat("  ✓ Saved: seurat_clustered.rds\n")

cat("\n=== Clustering and Cell Type Identification Complete ===\n")
cat(sprintf("\nSummary:\n"))
cat(sprintf("  - Filtered cells: %d\n", ncol(seurat_obj)))
cat(sprintf("  - Variable features: %d\n", length(VariableFeatures(seurat_obj))))
cat(sprintf("  - PCA dimensions: %d\n", n_dims))
cat(sprintf("  - Clusters identified: %d\n", length(unique(seurat_obj$seurat_clusters))))
cat(sprintf("  - Marker genes found: %d\n", nrow(cluster_markers)))
cat("\nKey figures generated:\n")
cat("  1. 06_umap_clusters.png - Main clustering visualization\n")
cat("  2. 08_immune_markers_umap.png - Immune cell marker expression\n")
