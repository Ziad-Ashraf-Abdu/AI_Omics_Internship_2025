# -----------------------------------------------------------------------------
# Script: 05_verification.R
# Purpose: Load artifacts from Steps A-D and perform sanity checks.
# *** UPDATED to check for GSE9660 ***
# -----------------------------------------------------------------------------

library(Biobase)
library(limma)

message("=======================================================")
message("   STARTING PIPELINE VERIFICATION")
message("=======================================================")

# --- Check Step A: Data Retrieval ---
if(file.exists("data/raw/raw_gse_objects.RData")){
  load("data/raw/raw_gse_objects.RData")
  message("\n[Step A] Check: Raw Objects Loaded")
  
  # Check dimensions of the count matrix
  d1 <- dim(gse64810_counts_matrix)
  # *** UPDATED: Check gse9660 ***
  d2 <- dim(gse9660)
  
  # Manually get sample count for gse64810 metadata
  gse64810_samples <- ncol(gse64810_metadata_obj)
  
  message(paste("  - GSE64810 (RNA-seq) Dimensions:", d1[1], "Genes x", d1[2], "Samples"))
  message(paste("  - GSE9660 (Array) Dimensions:   ", d2["Features"], "Probes x", d2["Samples"], "Samples"))
  
  if(d1[2] == gse64810_samples) message(paste("    -> PASS: GSE64810 sample count (", d1[2], ") matches methodology."))
  else message(paste("    -> WARNING: Sample count mismatch for GSE64810. Data has", d1[2], "samples, metadata has", gse64810_samples))
  
  if(d2["Samples"] == 20) message(paste("    -> PASS: GSE9660 sample count is 20 (10 HD + 10 Ctrl)."))
  else message(paste("    -> WARNING: GSE9660 sample count is not 20. Is", d2["Samples"]))
  
} else {
  message("\n[Step A] FAIL: 'data/raw/raw_gse_objects.RData' not found.")
}

# --- Check Step B: Preprocessing ---
if(file.exists("data/processed/normalized_data.RData")){
  load("data/processed/normalized_data.RData")
  message("\n[Step B] Check: Preprocessing & Batch Correction")
  
  # Check Class Balance
  print(table(merged_meta$Batch, merged_meta$Class))
  message("  - Look at the table above. You should see samples in both '0' and '1' columns for BOTH datasets.")
  
  # Check Data Range (Log2 check)
  rng <- range(corrected_expr)
  message(paste("  - Data Range: Min =", round(rng[1],2), "| Max =", round(rng[2],2)))
  
  if(rng[2] < 100) message("    -> PASS: Data appears to be Log2 transformed (Max < 100).")
  else message("    -> WARNING: Max value is high. Data might be raw counts (undesirable for ML).")
  
  # Check NAs
  na_count <- sum(is.na(corrected_expr))
  if(na_count == 0) message("    -> PASS: No missing values (NAs) in the matrix.")
  else message(paste("    -> WARNING: Found", na_count, "NAs in the matrix."))
  
} else {
  message("\n[Step B] FAIL: 'data/processed/normalized_data.RData' not found.")
}

# --- Check Step C: Differential Expression ---
if(file.exists("data/processed/dea_results.RData")){
  load("data/processed/dea_results.RData")
  message("\n[Step C] Check: Differential Expression")
  
  # Check number of DEGs
  n_degs <- nrow(deg_list)
  message(paste("  - Significant DEGs found (p < 0.05):", n_degs))
  
  if(n_degs > 100) message("    -> PASS: DEA identified a healthy number of significant genes.")
  else if (n_degs == 0) message("    -> WARNING: No significant genes found. Check normalization or model.")
  else message("    -> NOTE: Very few DEGs found. This might be correct or a sign of weak signal.")
  
  # Check top genes
  message("  - Top 3 Genes by LogFC:")
  
  # *** FIX: Add gene symbols as a column before trying to print them ***
  if(nrow(deg_list) > 0) {
    deg_list$Gene.symbol <- rownames(deg_list)
    print(head(deg_list[, c("logFC", "adj.P.Val", "Gene.symbol")], 3))
  } else {
    message("    -> No DEGs to display.")
  }
  
} else {
  message("\n[Step C] FAIL: 'data/processed/dea_results.RData' not found.")
}

# --- Check Step D: Feature Selection ---
if(file.exists("results/final_biomarker_dataset.csv")){
  df <- read.csv("results/final_biomarker_dataset.csv", row.names = 1)
  message("\n[Step D] Check: Final ML Dataset")
  
  # Check Dimensions
  # Rows = Total Samples, Cols = Features + Class
  dims <- dim(df)
  message(paste("  - Final Dataset:", dims[1], "Samples x", dims[2], "Columns"))
  
  # Check Target Column
  if("Class" %in% colnames(df)) message("    -> PASS: 'Class' column exists.")
  else message("    -> FAIL: 'Class' column missing.")
  
  # Check Feature Count
  # (Columns - 1 Class Column) should be close to 100 or your RFE limit
  feat_count <- dims[2] - 1
  message(paste("  - Feature Count:", feat_count))
  
  if(feat_count <= 100) message("    -> PASS: Feature selection reduced gene count to <= 100.")
  else message("    -> NOTE: Feature count > 100. RFE might have selected a larger optimal set.")
  
} else {
  message("\n[Step D] FAIL: 'results/final_biomarker_dataset.csv' not found.")
}

message("\n=======================================================")
message("   VERIFICATION COMPLETE")
message("================================S=======================")