# Download PBMC 3k dataset from 10x Genomics

cat("Downloading PBMC 3k dataset from 10x Genomics...\n\n")

# Download URL
url <- "https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz"
dest_file <- "/Volumes/MUKUTA/Projects/asgard_interview/portfolios/portfolio1_tumor_immune_landscape/data/pbmc3k.tar.gz"
data_dir <- "/Volumes/MUKUTA/Projects/asgard_interview/portfolios/portfolio1_tumor_immune_landscape/data"

cat(sprintf("Downloading from: %s\n", url))
cat(sprintf("Saving to: %s\n\n", dest_file))

download.file(url, destfile = dest_file, method = "curl")

cat("✓ Download complete\n\n")

cat("Extracting files...\n")
setwd(data_dir)
untar(dest_file)

cat("✓ Extraction complete\n\n")

cat("Files ready in: filtered_gene_bc_matrices/hg19/\n")
list.files("filtered_gene_bc_matrices/hg19/")
