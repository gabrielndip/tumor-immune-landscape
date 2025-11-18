# Download the actual expression data and annotations
library(GEOquery)

setwd("/Volumes/MUKUTA/Projects/asgard_interview/portfolios/portfolio1_tumor_immune_landscape/data")

cat("=== Downloading Expression Data ===\n\n")

# Download supplementary files
cat("Downloading GSE131907 supplementary files...\n")
cat("This may take 5-10 minutes depending on file sizes...\n\n")

tryCatch({
    # Download all supplementary files (or fetch_files = TRUE for specific files)
    supp_files <- getGEOSuppFiles("GSE131907",
                                  baseDir = ".",
                                  fetch_files = TRUE,
                                  makeDirectory = TRUE)

    cat("\n✓ Successfully downloaded supplementary files\n")
    cat("Downloaded files:\n")
    print(supp_files)

}, error = function(e) {
    cat("\n✗ Error downloading supplementary files:\n")
    cat(sprintf("  %s\n", e$message))
})

cat("\n=== Download Complete ===\n")
cat("Files saved in: data/GSE131907/\n")
list.files("GSE131907", full.names = FALSE)
