# -----------------------------------------------------------------------------
# Script: 01_data_retrieval.R (FIXED)
# Purpose: Download raw/processed data for GSE64810 (RNA-seq) and GSE9660 (Microarray).
#          GSE9660 replaces GSE33000 for a stable, balanced dataset.
# -----------------------------------------------------------------------------

dir.create("data/raw", recursive = TRUE, showWarnings = FALSE)
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("results", recursive = TRUE, showWarnings = FALSE)

if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require("GEOquery")) BiocManager::install("GEOquery")
if (!require("tidyverse")) install.packages("tidyverse")

library(GEOquery)
library(tidyverse)

message("[A] Starting Data Retrieval...")

# --- 1. GSE64810 (RNA-seq) ---
message("Downloading GSE64810 (Metadata)...")
# We still use this to get the metadata (pData)
gse64810_list <- getGEO("GSE64810", GSEMatrix = TRUE, destdir = "data/raw")
gse64810_metadata_obj <- gse64810_list[[1]]

message("Downloading GSE64810 (Supplementary Count Matrix)...")
# The actual count data is in a supplementary file
supp_files_path <- getGEOSuppFiles("GSE64810", baseDir = "data/raw")

# *** FIX: Point to the *correct* file provided by the authors ***
# This is the DESeq2 normalized data, which is what we want.
count_file_path <- file.path("data/raw", "GSE64810", "GSE64810_mlhd_DESeq2_norm_counts_adjust.txt.gz")

if (!file.exists(count_file_path)) {
  # Fallback in case directory structure is different
  message("File not found at expected path. Searching alternate path...")
  alt_path <- list.files("data/raw", pattern = "DESeq2_norm_counts_adjust.txt.gz", 
                         recursive = TRUE, full.names = TRUE)
  if(length(alt_path) > 0) {
    count_file_path <- alt_path[1]
  } else {
    stop("Could not find 'GSE64810_mlhd_DESeq2_norm_counts_adjust.txt.gz'. Download manually to 'data/raw/GSE64810/'")
  }
}

# Read the gzipped count file
message(paste("Reading count matrix from:", count_file_path))
# The first column is GeneID, which read_tsv correctly parses.
gse64810_counts <- read_tsv(count_file_path, comment = "#")

# The file has GeneID as col 1, and samples as others.
# Convert to a matrix with GeneIDs as rownames
gse64810_counts_matrix <- gse64810_counts %>% 
  column_to_rownames(var = "...1") %>% # *** FIX: Use the "...1" column name from the console output ***
  as.matrix()

# --- 2. GSE9660 (Microarray) ---
# *** UPDATED: Using GSE9660. It is balanced and uses the simple hgu133a chip. ***
message("Downloading GSE9660 (Microarray)...")
gse9660_list <- getGEO("GSE9660", GSEMatrix = TRUE, destdir = "data/raw")
gse9660 <- gse9660_list[[1]]

save(gse64810_counts_matrix, gse64810_metadata_obj, gse9660, 
     file = "data/raw/raw_gse_objects.RData")

message("[A] Data Retrieval Complete. Saved to data/raw/raw_gse_objects.RData")