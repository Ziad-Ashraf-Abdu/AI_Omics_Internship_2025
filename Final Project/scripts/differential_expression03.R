# ==============================================================================
# SCRIPT 2: DIFFERENTIAL EXPRESSION ANALYSIS
# FEATURES: 
#   1. Force Numeric Fix (for data loading)
#   2. Robust Labeling (Using Title 'C_' vs 'H_' prefixes)
# INPUTS: data/raw_counts.csv, data/metadata_clean.csv
# OUTPUTS: results/dds_object.rds, results/GSE129473_Full_DEA_Results.csv
# ==============================================================================

# 1. Install/Load Packages
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require("DESeq2", quietly = TRUE)) BiocManager::install("DESeq2", update = FALSE, ask = FALSE)

library(DESeq2)
library(dplyr)
library(readr)
library(tibble)

# Create output directory
if(!dir.exists("results")) dir.create("results")

message(">>> [1/5] Loading Data...")

# Define paths
counts_path <- "data/raw_counts.csv"
meta_path <- "data/metadata_clean.csv"

if(!file.exists(counts_path)) stop("Counts file not found! Run Script 1 first.")
if(!file.exists(meta_path)) stop("Metadata file not found! Run Script 1 first.")

# Read Data
# check.names=FALSE keeps headers exactly as they are (e.g., "C_0002")
raw_counts <- read.csv(counts_path, row.names = 1, check.names = FALSE)
metadata <- read.csv(meta_path) 

message(">>> [2/5] Cleaning & Validating Data...")

# ------------------------------------------------------------------------------
# 1. FORCE NUMERIC CONVERSION (The Numeric Fix)
# ------------------------------------------------------------------------------
# Convert all columns to numeric, suppressing warnings about NAs
raw_counts[] <- lapply(raw_counts, function(x) suppressWarnings(as.numeric(as.character(x))))

# Remove columns that became ALL NAs (text columns)
valid_cols <- colSums(is.na(raw_counts)) < nrow(raw_counts)
if (!all(valid_cols)) {
  bad_cols <- names(raw_counts)[!valid_cols]
  message(paste("WARNING: Removing non-numeric columns:", paste(head(bad_cols), collapse=", ")))
  raw_counts <- raw_counts[, valid_cols]
}

# Rounding for DESeq2 (Integers only)
raw_counts[is.na(raw_counts)] <- 0
raw_counts <- round(as.matrix(raw_counts))

message(paste(">>> Data Dimensions: ", nrow(raw_counts), "genes x", ncol(raw_counts), "samples"))
if(ncol(raw_counts) == 0) stop("CRITICAL ERROR: All columns were removed.")

# ------------------------------------------------------------------------------
# 2. ROBUST ID MATCHING (Metadata vs Counts)
# ------------------------------------------------------------------------------
count_cols <- colnames(raw_counts)
match_col <- NULL
max_overlap <- 0

for (col in colnames(metadata)) {
  common <- intersect(as.character(metadata[[col]]), count_cols)
  if (length(common) > max_overlap) {
    max_overlap <- length(common)
    match_col <- col
  }
}

if (max_overlap == 0) {
  stop("CRITICAL ERROR: No metadata column matches the count matrix headers.")
}

message(paste(">>> Matched ID column:", match_col))

# Filter metadata & Align
common_samples <- intersect(metadata[[match_col]], count_cols)
metadata <- metadata %>% 
  filter(.data[[match_col]] %in% common_samples) %>%
  arrange(match(.data[[match_col]], common_samples))

raw_counts <- raw_counts[, metadata[[match_col]]]

if(!all(colnames(raw_counts) == metadata[[match_col]])) stop("Alignment Error!")

# ------------------------------------------------------------------------------
# 3. ROBUST LABEL CORRECTION (The "C" vs "H" Fix)
# ------------------------------------------------------------------------------
message(">>> [3/5] Verifying Labels based on Titles...")

# Ensure 'title' column exists
if("title" %in% colnames(metadata)) {
  # Logic: 
  # If title starts with 'C' (case insensitive) -> Control
  # If title starts with 'H' (case insensitive) -> HD
  # Keep original if neither matches
  
  metadata$condition <- case_when(
    grepl("^[Cc]", metadata$title) ~ "Control",
    grepl("^[Hh]", metadata$title) ~ "HD",
    TRUE ~ as.character(metadata$condition)
  )
  
  message(">>> Label Distribution after Correction:")
  print(table(metadata$condition))
  
} else {
  message("WARNING: 'title' column not found. Skipping label correction.")
}

# ------------------------------------------------------------------------------
# 4. RUNNING DESEQ2
# ------------------------------------------------------------------------------
message(">>> [4/5] Initializing DESeq2...")

metadata$condition <- as.factor(metadata$condition)
if(length(levels(metadata$condition)) < 2) stop("Error: Condition must have 2 levels (HD vs Control).")

dds <- DESeqDataSetFromMatrix(countData = raw_counts,
                              colData = metadata,
                              design = ~ condition)

# Filter low-count genes
keep <- rowSums(counts(dds) >= 10) >= 10
dds <- dds[keep, ]

# Run Analysis
dds <- DESeq(dds)

# ------------------------------------------------------------------------------
# 5. SAVING OUTPUTS
# ------------------------------------------------------------------------------
message(">>> [5/5] Saving Results...")

res <- results(dds, contrast = c("condition", "HD", "Control"))
res_df <- as.data.frame(res) %>%
  rownames_to_column("GeneID") %>%
  arrange(padj)

saveRDS(dds, "results/dds_object.rds")
write_csv(res_df, "results/GSE129473_Full_DEA_Results.csv")

message(">>> SUCCESS: DEA Complete.")
message(">>> Ground Truth Labels have been corrected based on Sample Titles.")