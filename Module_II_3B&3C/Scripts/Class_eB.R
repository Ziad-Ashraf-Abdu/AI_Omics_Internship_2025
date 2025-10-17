# Microarray preprocessing and QC script for GSE79973
# AI & Omics Internship 2025 — Module II assignment
# NOTE TO INSTRUCTORS: Comments are intentionally human and addressed to you.

# --------------------------
# 0. Install & load packages
# --------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# Bioconductor packages
bioc_pkgs <- c("GEOquery", "affy", "arrayQualityMetrics", "genefilter", "AnnotationDbi")
for (p in bioc_pkgs) {
  if (!suppressWarnings(requireNamespace(p, quietly = TRUE))) {
    message("Installing Bioconductor package: ", p)
    BiocManager::install(p, ask = FALSE, update = FALSE)
  }
}

# CRAN packages
cran_pkgs <- c("dplyr", "matrixStats")
for (p in cran_pkgs) {
  if (!suppressWarnings(requireNamespace(p, quietly = TRUE))) {
    message("Installing CRAN package: ", p)
    install.packages(p)
  }
}

# Load libraries 
library(GEOquery)
library(affy)
library(arrayQualityMetrics)
library(dplyr)         #   in assignment
library(matrixStats)   # for rowMedians
# genefilter might be optional on some systems; load if available
if (suppressWarnings(requireNamespace("genefilter", quietly = TRUE))) library(genefilter)
library(AnnotationDbi)

options(stringsAsFactors = FALSE)

# --------------------------
# Parameters / file paths
# --------------------------
gse_acc <- "GSE79973"
results_dir <- "Results"
raw_dir <- "Raw_Data"

if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)
if (!dir.exists(raw_dir)) dir.create(raw_dir, recursive = TRUE)

# default threshold (kept as in teaching script)
default_threshold <- 3.5

# --------------------------
# 1) Download series matrix
# --------------------------
message("Downloading series matrix for ", gse_acc, " — may take a short while...")
gse_list <- getGEO(gse_acc, GSEMatrix = TRUE)
if (length(gse_list) == 0) stop("No GSE matrix found for ", gse_acc)
if (length(gse_list) > 1) message("Multiple ExpressionSets found; using the first one returned by getGEO()")
eset <- gse_list[[1]]

# variables 
expression_data <- exprs(eset)
feature_data    <- fData(eset)
phenotype_data  <- pData(eset)

# check missing values in sample annotation 
na_count_source <- sum(is.na(phenotype_data$source_name_ch1))
message("Missing values in phenotype_data$source_name_ch1: ", na_count_source)
message(sprintf("Series-matrix loaded: %d probes × %d samples", nrow(expression_data), ncol(expression_data)))

# --------------------------
# 2) Attempt to retrieve raw CEL files
# --------------------------
# Find any CEL files under Raw_Data (supports .CEL and .CEL.gz).
# IMPORTANT: correct escaping for R string below (two backslashes)
cel_candidates <- list.files(
  path = raw_dir,
  pattern = "\\.CEL(\\.gz)?$",   # <-- correct pattern (use double backslash)
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)
raw_found <- length(cel_candidates) > 0

if (raw_found) {
  message("Found CEL files locally (showing up to 10):")
  print(head(cel_candidates, 10))
  cel_dir <- normalizePath(dirname(cel_candidates[1]))
  raw_data <- ReadAffy(celfile.path = cel_dir)   # ReadAffy reads all CEL files in the directory
  message("Raw AffyBatch loaded — basic info:")
  print(raw_data)
} else {
  message("No CEL files found under '", raw_dir, "' — continuing with series matrix ExpressionSet")
}

# --------------------------
# 3) QC on raw (if available)
# --------------------------
if (exists("raw_data") && inherits(raw_data, "AffyBatch")) {
  qc_raw_dir <- file.path(results_dir, "QC_Raw_Data")
  if (!dir.exists(qc_raw_dir)) dir.create(qc_raw_dir, recursive = TRUE)
  message("Running arrayQualityMetrics on raw AffyBatch (writes HTML to: ", qc_raw_dir, ")")
  arrayQualityMetrics(expressionset = raw_data,
                      outdir = qc_raw_dir,
                      force = TRUE,
                      do.logtransform = TRUE)
} else {
  message("Skipping raw QC: no raw AffyBatch available.")
}

# --------------------------
# 4) RMA Normalization
# --------------------------
if (exists("raw_data") && inherits(raw_data, "AffyBatch")) {
  message("Performing RMA normalization on raw data...")
  normalized_eset <- rma(raw_data)
} else {
  message("No raw CELs — using series matrix ExpressionSet as normalized_eset (may already be preprocessed)")
  normalized_eset <- eset
}

normalized_data <- exprs(normalized_eset)    
write.csv(normalized_data, file = file.path(results_dir, "normalized_expression_matrix.csv"))
message("Normalized expression matrix saved to Results/normalized_expression_matrix.csv")

# QC after normalization
qc_norm_dir <- file.path(results_dir, "QC_Normalized_Data")
if (!dir.exists(qc_norm_dir)) dir.create(qc_norm_dir, recursive = TRUE)
arrayQualityMetrics(expressionset = normalized_eset, outdir = qc_norm_dir, force = TRUE)

# --------------------------
# 5) Filter low-intensity probes (soft intensity-based filtering)
# --------------------------
probe_median <- matrixStats::rowMedians(as.matrix(normalized_data), na.rm = TRUE)

# Save histogram
png(file = file.path(results_dir, "probe_median_hist.png"), width = 900, height = 650)
hist(probe_median, breaks = 100, freq = FALSE, main = "Median Intensity Distribution", xlab = "Median intensity")
abline(v = default_threshold, col = "black", lwd = 2)
dev.off()
message("Saved probe median histogram to Results/probe_median_hist.png (vertical line at threshold: ", default_threshold, ")")

# Apply intensity threshold
threshold <- default_threshold
keep_idx <- probe_median > threshold
filtered_data <- normalized_data[keep_idx, , drop = FALSE]
message(sprintf("Filtered by intensity: before = %d probes, after = %d probes", nrow(normalized_data), nrow(filtered_data)))

# Rename filtered expression data with sample metadata as  
if (!is.null(rownames(phenotype_data))) colnames(filtered_data) <- rownames(phenotype_data)

# --------------------------
# 6) Additional variance-based filtering
# --------------------------
# If genefilter available use varFilter, otherwise skip this step gracefully
if (suppressWarnings(requireNamespace("genefilter", quietly = TRUE))) {
  filtered_data_var <- genefilter::varFilter(filtered_data, var.cutoff = 0.5, filterByQuantile = TRUE)
  message(sprintf("After variance filtering (varFilter var.cutoff=0.5): %d probes remain", nrow(filtered_data_var)))
} else {
  message("genefilter not available — skipping variance filtering and using intensity-filtered data only")
  filtered_data_var <- filtered_data
}
write.csv(filtered_data_var, file = file.path(results_dir, "filtered_expression_matrix.csv"))
message("Saved filtered expression matrix to Results/filtered_expression_matrix.csv")

# --------------------------
# 7) Phenotype data preparation & groups
# --------------------------
if (!is.null(phenotype_data$source_name_ch1)) {
  message("Preparing experimental groups (assignment expects 'gastric mucosa' and 'gastric adenocarcinoma')")
  sn <- trimws(as.character(phenotype_data$source_name_ch1))
  unique_sn <- unique(sn)
  message("Unique values in source_name_ch1: ", paste(unique_sn, collapse = ", "))
  groups <- factor(sn,
                   levels = c("gastric mucosa", "gastric adenocarcinoma"),
                   labels = c("normal", "cancer"))
  message("Groups factor created. Levels: ", paste(levels(groups), collapse = ", "))
} else {
  message("phenotype_data$source_name_ch1 is not available — cannot create groups as  ")
  groups <- NULL
}

# --------------------------
# 8) Simple outlier detection (distance-based)
# --------------------------
sample_dists <- as.matrix(dist(t(filtered_data_var)))
mean_dists <- apply(sample_dists, 1, mean)
flag_thresh <- mean(mean_dists) + 2 * sd(mean_dists)
outliers <- names(mean_dists)[mean_dists > flag_thresh]
message(sprintf("Outlier detection threshold (mean + 2*sd): %.3f", flag_thresh))
if (length(outliers) > 0) {
  message("Potential outlier samples detected: ", paste(outliers, collapse = ", "))
} else {
  message("No potential outliers detected by this simple distance rule")
}

# --------------------------
# 9) Probe -> gene mapping (annotation)
# --------------------------
platform <- annotation(normalized_eset)
message("Detected platform: ", platform)
ann_pkg <- paste0(platform, ".db")
if (suppressWarnings(requireNamespace(ann_pkg, quietly = TRUE))) {
  library(ann_pkg, character.only = TRUE)
  probes <- rownames(filtered_data_var)
  mapped <- tryCatch({
    AnnotationDbi::select(get(ann_pkg), keys = probes, columns = c("SYMBOL", "GENENAME"), keytype = "PROBEID")
  }, error = function(e) {
    message("Annotation mapping failed: ", conditionMessage(e))
    NULL
  })
  if (!is.null(mapped)) {
    message("Annotation mapping head:")
    print(head(mapped))
  }
} else {
  message("Annotation package ", ann_pkg, " not installed. Install with BiocManager::install(\"", ann_pkg, "\") to map probes to genes.")
}

# --------------------------
# 10) Save processing summary
# --------------------------
summary_file <- file.path(results_dir, "processing_summary.txt")
cat("GEO accession:", gse_acc, "\n", file = summary_file)
cat("Series-matrix probes × samples:", nrow(expression_data), "×", ncol(expression_data), "\n", file = summary_file, append = TRUE)
if (exists("raw_data")) cat("Raw arrays:", ncol(raw_data), "\n", file = summary_file, append = TRUE)
cat("Normalized probes × samples:", nrow(normalized_data), "×", ncol(normalized_data), "\n", file = summary_file, append = TRUE)
cat("After intensity filtering:", nrow(filtered_data), "probes\n", file = summary_file, append = TRUE)
cat("After variance filtering:", nrow(filtered_data_var), "probes\n", file = summary_file, append = TRUE)
cat("Missing phenotype source_name_ch1 values:", na_count_source, "\n", file = summary_file, append = TRUE)
cat("Potential outlier samples (distance rule):", paste(outliers, collapse = ", "), "\n", file = summary_file, append = TRUE)

message("Processing complete.")

