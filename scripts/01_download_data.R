# Script to download scRNA-seq data for Portfolio 1
# Primary dataset: GSE131907 (Mouse melanoma tumor immune microenvironment)
# Backup: Use Seurat tutorial dataset if primary fails

library(GEOquery)

cat("=== Portfolio 1: Data Acquisition ===\n\n")

# Set working directory
setwd("/Volumes/MUKUTA/Projects/asgard_interview/portfolios/portfolio1_tumor_immune_landscape/data")

# Try to download GSE131907
cat("Attempting to download GSE131907 from GEO...\n")
cat("This may take several minutes...\n\n")

tryCatch({
    # Download GEO series
    gse <- getGEO("GSE131907", GSEMatrix = TRUE, getGPL = FALSE)

    cat("✓ Successfully retrieved GEO series metadata\n")
    cat(sprintf("  - Number of samples: %d\n", length(gse)))

    # Get sample information
    if (length(gse) > 0) {
        pData_info <- pData(gse[[1]])
        cat(sprintf("  - Sample metadata columns: %d\n", ncol(pData_info)))
        cat("\nSample titles:\n")
        print(head(pData_info$title, 10))
    }

    # Save metadata
    saveRDS(gse, "GSE131907_metadata.rds")
    cat("\n✓ Metadata saved to GSE131907_metadata.rds\n")

    # Note about downloading expression data
    cat("\n=== Note ===\n")
    cat("The full expression matrices are typically large and need to be\n")
    cat("downloaded from GEO supplementary files or SRA.\n")
    cat("We'll check for supplementary files...\n\n")

    # Get supplementary files info
    if (length(gse) > 0) {
        supp_info <- getGEOSuppFiles("GSE131907", fetch_files = FALSE)
        cat("Available supplementary files:\n")
        print(supp_info)
    }

}, error = function(e) {
    cat("\n✗ Error downloading from GEO:\n")
    cat(sprintf("  %s\n\n", e$message))
    cat("Will use backup strategy (Seurat tutorial dataset)\n")
})

cat("\n=== Data Download Status ===\n")
cat("Check data/ directory for downloaded files\n")
