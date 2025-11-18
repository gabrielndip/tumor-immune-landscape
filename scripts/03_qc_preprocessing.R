# Portfolio 1: scRNA-seq QC and Preprocessing
# Analysis of tumor immune landscape

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

cat("=== Portfolio 1: QC and Preprocessing ===\n\n")

# Set paths
data_dir <- "/Volumes/MUKUTA/Projects/asgard_interview/portfolios/portfolio1_tumor_immune_landscape/data"
results_dir <- "/Volumes/MUKUTA/Projects/asgard_interview/portfolios/portfolio1_tumor_immune_landscape/results"
figures_dir <- "/Volumes/MUKUTA/Projects/asgard_interview/portfolios/portfolio1_tumor_immune_landscape/figures"

setwd(data_dir)

# Create output directories if they don't exist
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

# Load data
cat("Loading scRNA-seq data...\n")

# Try to load the normalized matrix (RDS format)
tryCatch({
    # Load normalized expression matrix
    cat("  Loading normalized expression matrix...\n")
    expr_file <- "GSE131907/GSE131907_Lung_Cancer_normalized_log2TPM_matrix.rds.gz"

    if (file.exists(expr_file)) {
        expr_matrix <- readRDS(gzfile(expr_file))
        cat(sprintf("  ✓ Loaded expression matrix: %d genes x %d cells\n",
                    nrow(expr_matrix), ncol(expr_matrix)))
    } else {
        stop("Expression matrix file not found")
    }

    # Load cell annotations
    cat("  Loading cell annotations...\n")
    annot_file <- "GSE131907/GSE131907_Lung_Cancer_cell_annotation.txt.gz"

    if (file.exists(annot_file)) {
        cell_annot <- read.table(gzfile(annot_file), header = TRUE, sep = "\t")
        cat(sprintf("  ✓ Loaded annotations for %d cells\n", nrow(cell_annot)))
        cat(sprintf("    Available annotation columns: %s\n",
                    paste(colnames(cell_annot), collapse = ", ")))
    } else {
        stop("Annotation file not found")
    }

    # Create Seurat object
    cat("\nCreating Seurat object...\n")
    seurat_obj <- CreateSeuratObject(counts = expr_matrix,
                                     project = "LungCancer",
                                     min.cells = 3,
                                     min.features = 200)

    # Add metadata
    cat("Adding cell annotations to metadata...\n")
    # Match cell barcodes
    common_cells <- intersect(colnames(seurat_obj), rownames(cell_annot))
    cat(sprintf("  Matched %d cells between expression and annotation\n", length(common_cells)))

    # Add annotations for matched cells
    for (col in colnames(cell_annot)) {
        seurat_obj@meta.data[[col]] <- NA
        seurat_obj@meta.data[common_cells, col] <- cell_annot[common_cells, col]
    }

    # Calculate QC metrics
    cat("\nCalculating QC metrics...\n")
    # Calculate mitochondrial percentage (mouse: mt-, human: MT-)
    seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^[Mm][Tt]-")

    # Calculate ribosomal percentage
    seurat_obj[["percent.ribo"]] <- PercentageFeatureSet(seurat_obj, pattern = "^Rp[sl]")

    cat(sprintf("  ✓ QC metrics calculated\n"))
    cat(sprintf("    Median genes per cell: %.0f\n", median(seurat_obj$nFeature_RNA)))
    cat(sprintf("    Median UMIs per cell: %.0f\n", median(seurat_obj$nCount_RNA)))
    cat(sprintf("    Median %% MT: %.2f\n", median(seurat_obj$percent.mt)))

    # Create QC plots
    cat("\nGenerating QC visualizations...\n")

    # Violin plots for QC metrics
    p1 <- VlnPlot(seurat_obj,
                  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                  ncol = 3,
                  pt.size = 0.1) +
          plot_annotation(title = "Quality Control Metrics")

    ggsave(file.path(figures_dir, "01_qc_violin_plots.png"),
           plot = p1, width = 14, height = 5, dpi = 300)
    cat("  ✓ Saved: 01_qc_violin_plots.png\n")

    # Feature scatter plots
    p2 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt") +
          ggtitle("UMI Count vs Mitochondrial %")
    p3 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
          ggtitle("UMI Count vs Gene Count")

    p_combined <- p2 + p3
    ggsave(file.path(figures_dir, "02_qc_scatter_plots.png"),
           plot = p_combined, width = 12, height = 5, dpi = 300)
    cat("  ✓ Saved: 02_qc_scatter_plots.png\n")

    # Save QC object
    cat("\nSaving Seurat object with QC metrics...\n")
    saveRDS(seurat_obj, file.path(results_dir, "seurat_qc.rds"))
    cat("  ✓ Saved: seurat_qc.rds\n")

    # Summary statistics
    cat("\n=== QC Summary ===\n")
    cat(sprintf("Total cells: %d\n", ncol(seurat_obj)))
    cat(sprintf("Total genes: %d\n", nrow(seurat_obj)))
    cat(sprintf("\nQC Thresholds Recommendation:\n"))
    cat(sprintf("  nFeature_RNA: 200 - 5000\n"))
    cat(sprintf("  percent.mt: < 20%%\n"))
    cat(sprintf("\nCells passing filters: %d\n",
                sum(seurat_obj$nFeature_RNA > 200 &
                    seurat_obj$nFeature_RNA < 5000 &
                    seurat_obj$percent.mt < 20)))

}, error = function(e) {
    cat("\n✗ Error in data loading/processing:\n")
    cat(sprintf("  %s\n", e$message))
    cat("\nWill use backup dataset strategy\n")
})

cat("\n=== QC and Preprocessing Complete ===\n")
