# -----------------------------------------------------------------------------
# Script: 02_preprocessing.R (FIXED)
# Purpose: Log2, Normalize, Merge, ComBat Correction, and BEFORE/AFTER PCA.
#          Uses GSE1750 instead of GSE9660.
# -----------------------------------------------------------------------------

# --- Automatic Package Installation ---
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require("Biobase")) BiocManager::install("Biobase")
if (!require("tidyverse")) install.packages("tidyverse")
if (!require("sva")) BiocManager::install("sva")
# *** UPDATED: We are back to hgu133a.db, which is simple and works. ***
if (!require("hgu133a.db")) BiocManager::install("hgu133a.db")
if (!require("AnnotationDbi")) BiocManager::install("AnnotationDbi")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("org.Hs.eg.db")) BiocManager::install("org.Hs.eg.db")
# --------------------------------------

library(Biobase)
library(tidyverse)
library(sva)          # For ComBat
# *** UPDATED: Load the correct, simple library for GSE1750 (GPL96) ***
library(hgu133a.db) 
library(AnnotationDbi)
library(ggplot2)
library(org.Hs.eg.db) # For mapping GSE64810

load("data/raw/raw_gse_objects.RData")

message("[B] Starting Preprocessing...")

# --- 1. Process GSE64810 (RNA-seq) ---
# ... (This section is correct, no changes needed) ...
# Use the count matrix we loaded in Step 1
exp_64810 <- gse64810_counts_matrix

if(max(exp_64810, na.rm=T) > 100) {
  message("Log2 transforming GSE64810...")
  exp_64810 <- log2(exp_64810 + 1)
} else {
  message("GSE64810 data appears pre-normalized/log-transformed. Skipping log2 transform.")
}

# The rownames are Ensembl IDs (e.g., ENSG00000000003)
# We MUST map them to Gene Symbols
ensembl_ids <- rownames(exp_64810)
# We remove the ".number" from Ensembl IDs if they exist
ensembl_ids_clean <- gsub("\\..*$", "", ensembl_ids)

symbols_64810 <- mapIds(org.Hs.eg.db, 
                        keys=ensembl_ids_clean, 
                        column="SYMBOL", 
                        keytype="ENSEMBL", 
                        multiVals="first")

mapping_64810 <- data.frame(ensembl = names(symbols_64810), 
                            symbol = unname(symbols_64810), 
                            stringsAsFactors = F) %>% na.omit()

# Aggregate duplicates (multiple ensembl IDs mapping to one symbol)
exp_64810_df <- as.data.frame(exp_64810)
# Match using the original (non-clean) ensembl_ids
exp_64810_df$ensembl_clean <- ensembl_ids_clean
exp_64810_df$symbol <- mapping_64810$symbol[match(exp_64810_df$ensembl_clean, mapping_64810$ensembl)]

exp_64810_df <- exp_64810_df %>% 
  dplyr::select(-ensembl_clean) %>% # remove helper column
  filter(!is.na(symbol)) %>% 
  group_by(symbol) %>% 
  summarize_all(mean) # Average duplicates

exp_64810_mat <- as.matrix(exp_64810_df[,-1])
rownames(exp_64810_mat) <- exp_64810_df$symbol

# Get metadata from the object we downloaded
meta_64810 <- pData(gse64810_metadata_obj)
meta_64810$Class <- ifelse(grepl("Huntington", meta_64810$characteristics_ch1.1, ignore.case=T), 1, 0)
meta_64810$Batch <- "GSE64810"


# --- 2. Process GSE1750 (Microarray) ---
# *** UPDATED to use gse1750 ***
exp_1750 <- exprs(gse1750)
probes <- rownames(exp_1750)

# *** MAJOR FIX: Use the hgu133a.db package and PROBEID keytype ***
# This is the standard, simple, correct mapping that will work.
symbols_1750 <- mapIds(hgu133a.db, 
                       keys=probes, 
                       column="SYMBOL", 
                       keytype="PROBEID", # This is the correct, valid keytype for this package
                       multiVals="first")

mapping_1750 <- data.frame(probe=names(symbols_1750), symbol=unname(symbols_1750), stringsAsFactors=F) %>% na.omit()

exp_1750_df <- as.data.frame(exp_1750)
exp_1750_df$symbol <- mapping_1750$symbol[match(rownames(exp_1750_df), mapping_1750$probe)]
exp_1750_df <- exp_1750_df %>% filter(!is.na(symbol)) %>% group_by(symbol) %>% summarize_all(mean)
exp_1750_mat <- as.matrix(exp_1750_df[,-1])
rownames(exp_1750_mat) <- exp_1750_df$symbol

meta_1750 <- pData(gse1750)

# This dataset (GSE1750) reliably uses 'characteristics_ch1'
# for "disease state: Huntington's disease" or "disease state: control"
meta_1750$Class <- ifelse(grepl("Huntington", meta_1750$characteristics_ch1, ignore.case=T), 1, 0)
meta_1750$Batch <- "GSE1750"

# --- 3. Merge ---
# Find common genes (now that both have SYMBOLS)
common_genes <- intersect(rownames(exp_64810_mat), rownames(exp_1750_mat))
message(paste("Merging datasets on", length(common_genes), "common genes..."))

# Subset both matrices to *only* the common genes
merged_expr <- cbind(exp_64810_mat[common_genes,], exp_1750_mat[common_genes,])

# FIX for 'select' error: Use dplyr::select to be explicit
merged_meta <- rbind(
  meta_64810 %>% dplyr::select(Class, Batch), 
  meta_1750 %>% dplyr::select(Class, Batch)
)

# --- 4. PCA: BEFORE Correction ---
message("Generating 'Before' PCA plot...")
pca_before <- prcomp(t(merged_expr), scale. = TRUE) 
df_before <- as.data.frame(pca_before$x)
df_before$Class <- as.factor(merged_meta$Class)
df_before$Batch <- as.factor(merged_meta$Batch)

plot_before <- ggplot(df_before, aes(x=PC1, y=PC2, color=Batch, shape=Class)) +
  geom_point(size=3, alpha=0.7) +
  theme_minimal() +
  labs(title="BEFORE Correction", subtitle="Samples cluster by Batch (Bad)",
       color="Dataset", shape="Class (0=Ctrl, 1=HD)")

ggsave("results/pca_before.png", plot_before)

# --- 5. ComBat Correction ---
message("Running ComBat...")
# Ensure Class is a factor for the model matrix
merged_meta$Class <- as.factor(merged_meta$Class)
mod <- model.matrix(~Class, data=merged_meta)
corrected_expr <- ComBat(dat=merged_expr, batch=merged_meta$Batch, mod=mod)

# --- 6. PCA: AFTER Correction ---
message("Generating 'After' PCA plot...")
pca_after <- prcomp(t(corrected_expr), scale. = TRUE)
df_after <- as.data.frame(pca_after$x)
df_after$Class <- as.factor(merged_meta$Class)
df_after$Batch <- as.factor(merged_meta$Batch)

plot_after <- ggplot(df_after, aes(x=PC1, y=PC2, color=Class, shape=Batch)) +
  geom_point(size=3, alpha=0.7) +
  theme_minimal() +
  labs(title="AFTER Correction", subtitle="Samples cluster by Biology (Good)",
       color="Class (0=Ctrl, 1=HD)", shape="Dataset")

ggsave("results/pca_after.png", plot_after)

save(corrected_expr, merged_meta, file="data/processed/normalized_data.RData")
message("[B] Preprocessing Complete. Plots saved: results/pca_before.png & results/pca_after.png")