# =====================================================================
#         AI and Biotechnology / Bioinformatics
# =====================================================================
#             Module II: Microarray Data Analysis
#                    Assignment: GSE79973
# =====================================================================

# Topics:
# 1. Probe IDs to gene mapping
# 2. Differential Gene Expression Analysis
# 3. Data Visualization

# Set your working directory 
setwd("C:/AI_Omics_Internship_2025/Module_II_3B") 

gc() # Clear memory

# =====================================================================
#### Install and Load Required Packages (Efficiently) ####
# =====================================================================

# List of required CRAN packages
cran_packages <- c("dplyr", "tibble", "ggplot2", "pheatmap")
for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}

# --- Bioconductor Packages ---
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
bioc_packages <- c("limma", "AnnotationDbi", "hgu133plus2.db", "Biobase")
for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg)
}

# Load all the libraries
library(limma)
library(AnnotationDbi)
library(hgu133plus2.db)
library(dplyr)
library(tibble)
library(ggplot2)
library(pheatmap)
library(Biobase)

cat("All required packages are installed and loaded.\n")

# Create directories for results if they don't exist
if (!dir.exists("Results")) dir.create("Results")
if (!dir.exists("Result_Plots")) dir.create("Result_Plots")


# =====================================================================
#### Data Loading and Preparation ####
# =====================================================================

# Load your preprocessed expression and phenotype data for GSE79973
load(".RData") 

# Check the main data object
# From our previous step, we know the main object is 'eset'
# Let's inspect it to understand the experiment design
print(eset)
head(pData(eset)) # Display phenotype data to find the grouping variable


# =====================================================================
#### 1. Probe-to-Gene Mapping ####
# =====================================================================

# Extract expression data from the eset object
expression_matrix <- exprs(eset)

# Get probe IDs from the rownames
probe_ids <- rownames(expression_matrix)

# Map probe IDs to gene symbols using the annotation database
# We'll retrieve GENESYMBOL for each PROBEID
annotation <- AnnotationDbi::select(hgu133plus2.db, 
                                    keys = probe_ids,
                                    columns = "SYMBOL",
                                    keytype = "PROBEID")

# The annotation can have probes with no gene symbol (NA) or one probe mapping to multiple symbols.
# We'll also have multiple probes mapping to the same gene. Let's clean this up.

# Remove probes with no gene symbol
annotation <- annotation %>% filter(!is.na(SYMBOL))

# For genes mapped by multiple probes, we'll select the probe with the highest average expression.
# This is a common strategy to ensure one probe-per-gene.
avg_expression <- data.frame(PROBEID = rownames(expression_matrix), AVG_EXPR = rowMeans(expression_matrix))
annotation <- annotation %>%
  inner_join(avg_expression, by = "PROBEID") %>%
  group_by(SYMBOL) %>%
  summarise(PROBEID = PROBEID[which.max(AVG_EXPR)], .groups = 'drop')

# Filter the main expression matrix to keep only the unique, annotated probes
annotated_expression_matrix <- expression_matrix[annotation$PROBEID, ]

# Assign gene symbols as the new rownames
rownames(annotated_expression_matrix) <- annotation$SYMBOL

cat("\nProbe-to-gene mapping complete. Expression matrix now uses gene symbols.\n")
head(annotated_expression_matrix)

# =====================================================================
#### 2. Differential Gene Expression Analysis with limma ####
# =====================================================================

# Create a factor for the groups from the phenotype data
groups <- factor(pData(eset)$source_name_ch1)
print(levels(groups)) 

# Create the design matrix
design <- model.matrix(~0 + groups)

# --- FIX: Make column names syntactically valid ---
# This changes "gastric adenocarcinoma" to "gastric.adenocarcinoma"
colnames(design) <- make.names(levels(groups))

# Fit the linear model
fit <- lmFit(annotated_expression_matrix, design)

# --- FIX: Define the contrast using your actual, corrected group names ---
# This comparison will find genes up-regulated in adenocarcinoma compared to normal mucosa.
contrast_matrix <- makeContrasts(gastric.adenocarcinoma - gastric.mucosa, levels = design)

# Fit the contrast
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# Get the table of top differentially expressed genes
top_genes <- topTable(fit2, number = Inf, sort.by = "P")

# Save the full results table
write.csv(top_genes, file = "Results/DGE_Results_GSE79973.csv")

cat("\nDifferential expression analysis complete. Results saved to Results/ folder.\n")
head(top_genes)

# =====================================================================
#### 3. Data Visualization ####
# =====================================================================

# --- Volcano Plot ---
# Add a column to classify genes as Up-regulated, Down-regulated, or Not Significant
top_genes <- top_genes %>%
  mutate(Significance = case_when(
    logFC > 1 & adj.P.Val < 0.05  ~ "Up-regulated",
    logFC < -1 & adj.P.Val < 0.05 ~ "Down-regulated",
    TRUE                         ~ "Not Significant"
  ))

# Create the plot
volcano_plot <- ggplot(top_genes, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = Significance), alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Up-regulated" = "red", "Down-regulated" = "blue", "Not Significant" = "grey")) +
  theme_minimal() +
  labs(
    title = "Volcano Plot for GSE79973",
    x = "Log2 Fold Change",
    y = "-log10(Adjusted P-value)"
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed")

# Save the plot
ggsave("Result_Plots/Volcano_Plot_GSE79973.png", plot = volcano_plot, width = 8, height = 6)
cat("\nVolcano plot saved to Result_Plots/ folder.\n")

# --- Heatmap ---
# Select the top 50 most significant genes for the heatmap
top_50_genes <- head(top_genes, 50)
heatmap_matrix <- annotated_expression_matrix[rownames(top_50_genes), ]

# Create an annotation dataframe for the columns (samples)
sample_annotation <- data.frame(Group = groups)
rownames(sample_annotation) <- colnames(heatmap_matrix)

# Generate the heatmap
# We scale by row to see the relative expression changes for each gene
pheatmap(
  heatmap_matrix,
  main = "Heatmap of Top 50 Differentially Expressed Genes",
  scale = "row",
  annotation_col = sample_annotation,
  show_colnames = FALSE, # Hide sample names if there are too many
  show_rownames = TRUE,
  fontsize_row = 8,
  filename = "Result_Plots/Heatmap_Top50_GSE79973.png"
)

cat("\nHeatmap saved to Result_Plots/ folder.\n")
cat("\n\nAnalysis finished successfully! 🎉\n")