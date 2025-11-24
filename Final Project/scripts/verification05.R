# ==============================================================================
# SCRIPT 3: PREPARE VALIDATION DATA (VERSION-AGNOSTIC MATCHING)
# PURPOSE: Align Validation Data (GSE129473) to Training Features (GSE64810)
#          ignoring version suffixes (e.g., .10 vs .9)
# INPUTS: 
#   1. results/dds_object.rds (Validation Data)
#   2. HD_ML_Ready_Data.csv   (Training Data)
# OUTPUTS: 
#   1. GSE129473_ML_Ready_Validation.csv
# ==============================================================================

library(DESeq2)
library(dplyr)
library(readr)

# Function to strip version numbers (ENSG000001.10 -> ENSG000001)
strip_version <- function(ids) {
  sub("\\.[0-9]+$", "", ids)
}

message(">>> [1/5] Loading Data...")

if(!file.exists("results/dds_object.rds")) stop("Error: results/dds_object.rds not found.")
dds_val <- readRDS("results/dds_object.rds")

if(!file.exists("HD_ML_Ready_Data.csv")) stop("Error: HD_ML_Ready_Data.csv not found.")
train_data <- read.csv("HD_ML_Ready_Data.csv", nrows = 1, check.names = FALSE)

message(">>> [2/5] Extracting & Cleaning Feature IDs...")

# 1. Get Training Genes and Strip Versions
train_cols_original <- colnames(train_data)
required_genes_original <- train_cols_original[train_cols_original != "Target_Class"]
required_genes_base <- strip_version(required_genes_original)

message(paste(">>> Model requires", length(required_genes_original), "genes."))

# 2. Get Validation Genes and Strip Versions
val_genes_original <- rownames(dds_val)
val_genes_base <- strip_version(val_genes_original)

# Check Overlap
overlap <- intersect(required_genes_base, val_genes_base)
message(paste(">>> Found", length(overlap), "matching genes (ignoring versions)."))

if(length(overlap) < 100) {
  stop("CRITICAL ERROR: Very low gene overlap (<100) even after stripping versions. Check datasets.")
}

message(">>> [3/5] Processing Validation Data...")

# Variance Stabilizing Transformation
vst_val <- vst(dds_val, blind = FALSE)
val_matrix <- assay(vst_val)

# Map Validation Rows to Base IDs
rownames(val_matrix) <- val_genes_base

# ------------------------------------------------------------------------------
# CONSTRUCT FINAL MATRIX
# ------------------------------------------------------------------------------
# We need to reconstruct the matrix using the ORIGINAL Training Column Names
# but pulling data using the BASE Gene IDs.

final_matrix <- matrix(0, nrow = length(required_genes_original), ncol = ncol(val_matrix))
rownames(final_matrix) <- required_genes_original
colnames(final_matrix) <- colnames(val_matrix)

# Fill the matrix
found_count <- 0
for(i in seq_along(required_genes_original)) {
  original_id <- required_genes_original[i]
  base_id <- required_genes_base[i]
  
  if(base_id %in% rownames(val_matrix)) {
    # If multiple validation genes map to same base ID (rare), take the first one
    final_matrix[i, ] <- val_matrix[base_id, , drop=FALSE][1, ]
    found_count <- found_count + 1
  }
}

message(paste(">>> Successfully filled", found_count, "out of", length(required_genes_original), "genes."))

message(">>> [4/5] Final Formatting...")

# Transpose: Rows = Samples, Columns = Genes
ml_input_val <- t(final_matrix)
ml_input_val <- as.data.frame(ml_input_val)

# Add Target Label
ml_input_val$Target_Class <- colData(dds_val)$condition

# Reorder columns to match training data exactly (Target + Original Gene IDs)
ml_input_val <- ml_input_val %>% select(all_of(train_cols_original))

message(">>> [5/5] Verifying Column Order...")

if (identical(colnames(ml_input_val), train_cols_original)) {
  message("SUCCESS: Validation columns match Training columns EXACTLY.")
} else {
  stop("CRITICAL ERROR: Column mismatch!")
}

# Save
outfile <- "GSE129473_ML_Ready_Validation.csv"
write.csv(ml_input_val, outfile, row.names = FALSE)

message("========================================================")
message("VALIDATION DATA PREPARED SUCCESSFULLY.")
message(paste("File saved:", outfile))
message("========================================================")