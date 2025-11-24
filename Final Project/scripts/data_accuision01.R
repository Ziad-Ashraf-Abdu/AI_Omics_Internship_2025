# ==============================================================================
# SCRIPT 1: DATA RETRIEVAL (FIXED FOR GSE129473)
# DATASET: GSE129473 (Human Caudate Nucleus, Normalized Counts)
# OUTPUTS: data/raw_counts.csv, data/metadata_clean.csv
# ==============================================================================

library(GEOquery)
library(dplyr)
library(readr)

# Create a data directory
if(!dir.exists("data")) dir.create("data")

message(">>> [1/3] Fetching Clinical Metadata from GEO...")

# Fetch the Series Matrix
gse <- getGEO("GSE129473", GSEMatrix = TRUE)
if (length(gse) > 1) idx <- grep("GPL", attr(gse, "names")) else idx <- 1
gse_data <- gse[[idx]]

# Extract Metadata
metadata <- pData(gse_data)

# ------------------------------------------------------------------------------
# ROBUST CONDITION DETECTION
# ------------------------------------------------------------------------------
# Find the column describing the disease state
target_col_idx <- which(apply(metadata, 2, function(x) any(grepl("Huntington", x, ignore.case = TRUE))))

if(length(target_col_idx) == 0) {
  print(colnames(metadata))
  stop("CRITICAL ERROR: Could not find any column containing 'Huntington' in the metadata.")
}

target_col <- names(target_col_idx)[1]
message(paste(">>> Found disease status info in column:", target_col))

# ------------------------------------------------------------------------------
# CLEANING & SAVING
# ------------------------------------------------------------------------------
clean_metadata <- metadata %>%
  dplyr::mutate(condition = as.factor(ifelse(grepl("Huntington", .[[target_col]], ignore.case = TRUE), "HD", "Control")))

write_csv(clean_metadata, "data/metadata_clean.csv")
message(">>> Metadata saved to 'data/metadata_clean.csv'")

message(">>> [2/3] Fetching Count Matrix...")

# Download supplementary files if not already present
# (The code checks if the file exists locally first to avoid re-downloading)
files_in_dir <- list.files("GSE129473", full.names = TRUE)
if(length(files_in_dir) == 0) {
  getGEOSuppFiles("GSE129473")
}

# ------------------------------------------------------------------------------
# FIXED FILE DETECTION LOGIC
# ------------------------------------------------------------------------------
# Look for "norm" or "csv" instead of just "counts"
count_file_path <- list.files("GSE129473", pattern = "norm", full.names = TRUE)

# Fallback: look for any csv.gz if "norm" isn't found
if(length(count_file_path) == 0) {
  count_file_path <- list.files("GSE129473", pattern = "csv.gz", full.names = TRUE)
}

if(length(count_file_path) == 0) stop("Error: No count file found in GSE129473 folder.")

count_file <- count_file_path[1]
message(paste(">>> Reading count file:", count_file))

# ------------------------------------------------------------------------------
# READING THE DATA (Using read.csv for CSV format)
# ------------------------------------------------------------------------------
# The file is comma-separated, so we use read.csv
raw_counts <- read.csv(gzfile(count_file), row.names = 1)

message(paste(">>> Data loaded. Dimensions:", nrow(raw_counts), "genes x", ncol(raw_counts), "samples"))

# Note: These are likely ALREADY normalized.
# We save them to data/raw_counts.csv to maintain pipeline consistency,
# but treat them as your input expression matrix.
write.csv(raw_counts, "data/raw_counts.csv", row.names = TRUE)

message(">>> [3/3] Data Retrieval Complete. Files saved in 'data/' folder.")