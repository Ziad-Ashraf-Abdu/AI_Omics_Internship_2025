# ===================================================================
#               AI and Biotechnology / Bioinformatics
# ===================================================================
#                        class_2b.R
#           Differential Gene Expression Analysis
# ===================================================================

# Setting working directory
setwd("C:/AI_Omics_Internship_2025/Module_I")

# Defining directory paths
raw_data_dir <- "raw_data"
results_dir <- "results"

# -------------------------------------------------------------------
#                    Downloading Data from GitHub
# -------------------------------------------------------------------

cat("Downloading data files from GitHub repository...\n")

# Setting URLs for the CSV files
url_1 <- "https://raw.githubusercontent.com/AI-Biotechnology-Bioinformatics/AI_and_Omics_Research_Internship_2025/main/DEGs_Data_1.csv"
url_2 <- "https://raw.githubusercontent.com/AI-Biotechnology-Bioinformatics/AI_and_Omics_Research_Internship_2025/main/DEGs_Data_2.csv"

# Setting File paths for downloaded data
file_path_1 <- file.path(raw_data_dir, "DEGs_Data_1.csv")
file_path_2 <- file.path(raw_data_dir, "DEGs_Data_2.csv")

# Downloading files
tryCatch({
  download.file(url_1, file_path_1, mode = "wb")
  cat("DEGs_Data_1.csv downloaded successfully.\n")
}, error = function(e) {
  cat("Error downloading DEGs_Data_1.csv:", e$message, "\n")
})

tryCatch({
  download.file(url_2, file_path_2, mode = "wb")
  cat("DEGs_Data_2.csv downloaded successfully.\n")
}, error = function(e) {
  cat("Error downloading DEGs_Data_2.csv:", e$message, "\n")
})

# -------------------------------------------------------------------
#                    Function Definition
# -------------------------------------------------------------------

classify_gene <- function(logFC, padj) {
  # Checking if logFC > 1 and padj < 0.05 for upregulation
  if (logFC > 1 & padj < 0.05) {
    return("Upregulated")
  }
  # Checking if logFC < -1 and padj < 0.05 for downregulation
  else if (logFC < -1 & padj < 0.05) {
    return("Downregulated")
  }
  else {
    return("Not_Significant")
  }
}

# Testing the function
# cat("\nTesting classify_gene function:\n")
# cat("Example 1 - logFC = 2.5, padj = 0.01:", classify_gene(2.5, 0.01), "\n")
# cat("Example 2 - logFC = -2.0, padj = 0.03:", classify_gene(-2.0, 0.03), "\n")
# cat("Example 3 - logFC = 0.8, padj = 0.001:", classify_gene(0.8, 0.001), "\n")

# -------------------------------------------------------------------
#                    Processing Loop
# -------------------------------------------------------------------

# List of files to process
files_to_process <- c("DEGs_Data_1.csv", "DEGs_Data_2.csv")

# Preparing empty list to store results
result_list <- list()

# Processing each file in a loop
for (file_name in files_to_process) {
  cat("\n", rep("=", 60), "\n")
  cat("Processing:", file_name, "\n")
  cat(rep("=", 60), "\n")
  
  # Creating input file path
  input_file_path <- file.path(raw_data_dir, file_name)
  
  # Checking if file exists
  if (!file.exists(input_file_path)) {
    cat("Error: File", file_name, "not found in", raw_data_dir, "\n")
    next
  }
  
  # Importing dataset
  data <- read.csv(input_file_path, header = TRUE)
  cat("File imported successfully.\n")
  cat("Dataset dimensions:", nrow(data), "rows,", ncol(data), "columns\n")
  
  # Display column names
  cat("Column names:", paste(names(data), collapse = ", "), "\n")
  
  # Checking structure of the data
  cat("\nDataset structure:\n")
  str(data)
  
  # Checking for missing values in padj column
  missing_padj <- sum(is.na(data$padj))
  cat("\nMissing values in 'padj' column:", missing_padj, "\n")
  
  # Replacing missing padj values with 1
  if (missing_padj > 0) {
    data$padj[is.na(data$padj)] <- 1
    cat("Missing padj values replaced with 1.\n")
  }
  
  # Checking for missing values in logFC column
  missing_logfc <- sum(is.na(data$logFC))
  cat("Missing values in 'logFC' column:", missing_logfc, "\n")
  
  # Initializing status column
  data$status <- character(nrow(data))
  
  # Applying classify_gene function to each row using for-loop
  cat("Classifying genes using for-loop...\n")
  for (i in 1:nrow(data)) {
    data$status[i] <- classify_gene(data$logFC[i], data$padj[i])
  }
  cat("Gene classification completed.\n")
  
  # Generating summary statistics
  cat("\n--- SUMMARY STATISTICS ---\n")
  gene_summary <- table(data$status)
  print(gene_summary)
  
  # Calculating and displaying percentages
  total_genes <- nrow(data)
  cat("\nDetailed Summary for", file_name, ":\n")
  cat("Total genes processed:", total_genes, "\n")
  
  if ("Upregulated" %in% names(gene_summary)) {
    cat("Upregulated genes:", gene_summary["Upregulated"], 
        "(", round(gene_summary["Upregulated"]/total_genes * 100, 1), "%)\n")
  } else {
    cat("Upregulated genes: 0 (0.0%)\n")
  }
  
  if ("Downregulated" %in% names(gene_summary)) {
    cat("Downregulated genes:", gene_summary["Downregulated"], 
        "(", round(gene_summary["Downregulated"]/total_genes * 100, 1), "%)\n")
  } else {
    cat("Downregulated genes: 0 (0.0%)\n")
  }
  
  if ("Not_Significant" %in% names(gene_summary)) {
    cat("Not significant genes:", gene_summary["Not_Significant"], 
        "(", round(gene_summary["Not_Significant"]/total_genes * 100, 1), "%)\n")
  } else {
    cat("Not significant genes: 0 (0.0%)\n")
  }
  
  # Store results in R list
  result_list[[file_name]] <- data
  
  # Save processed file to results folder
  output_file_name <- paste0("Processed_", file_name)
  output_file_path <- file.path(results_dir, output_file_name)
  write.csv(data, output_file_path, row.names = FALSE)
  cat("Results saved to:", output_file_path, "\n")
}

# -------------------------------------------------------------------
#                    Final Combined Analysis
# -------------------------------------------------------------------

cat("\n", rep("=", 70), "\n")
cat("                    FINAL COMBINED ANALYSIS")
cat("\n", rep("=", 70), "\n")

# Accessing individual results
if (length(result_list) >= 2) {
  results_1 <- result_list[["DEGs_Data_1.csv"]]
  results_2 <- result_list[["DEGs_Data_2.csv"]]
  
  # Display top 5 rows of each processed dataset
  cat("\nFirst 5 rows of processed DEGs_Data_1:\n")
  print(head(results_1, 5))
  
  cat("\nFirst 5 rows of processed DEGs_Data_2:\n")
  print(head(results_2, 5))
  
  # Combined analysis
  all_data <- rbind(results_1, results_2)
  combined_summary <- table(all_data$status)
  
  cat("\nCombined Summary Across Both Datasets:\n")
  print(combined_summary)
  
  total_combined <- nrow(all_data)
  cat("\nOverall Statistics:\n")
  cat("Total genes processed:", total_combined, "\n")
  
  if ("Upregulated" %in% names(combined_summary)) {
    cat("Total upregulated genes:", combined_summary["Upregulated"], 
        "(", round(combined_summary["Upregulated"]/total_combined * 100, 1), "%)\n")
  }
  
  if ("Downregulated" %in% names(combined_summary)) {
    cat("Total downregulated genes:", combined_summary["Downregulated"], 
        "(", round(combined_summary["Downregulated"]/total_combined * 100, 1), "%)\n")
  }
  
  if ("Not_Significant" %in% names(combined_summary)) {
    cat("Total not significant genes:", combined_summary["Not_Significant"], 
        "(", round(combined_summary["Not_Significant"]/total_combined * 100, 1), "%)\n")
  }
  
  # Save combined results
  combined_output_path <- file.path(results_dir, "Combined_DEGs_Results.csv")
  write.csv(all_data, combined_output_path, row.names = FALSE)
  cat("\nCombined results saved to:", combined_output_path, "\n")
}