# -----------------------------------------------------------------------------
# Script: 03_dea.R (Corresponds to Methodology Section C)
# Purpose: Compare HD vs Control to compute fold-changes and adjusted p-values.
#          Select candidate gene list (~100-500 DEGs).
# -----------------------------------------------------------------------------

# --- Automatic Package Installation ---
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require("limma")) BiocManager::install("limma")
# --------------------------------------

library(limma)
load("data/processed/normalized_data.RData")

message("[C] Starting Differential Expression Analysis (DEA)...")

# Methodology C says "compare HD vs control".
# Since we have merged continuous data (Corrected Array + RNA-seq), we use Limma.
# Limma is the mathematical equivalent to edgeR/DESeq2 for continuous/batch-corrected data.

# 1. Create Design Matrix
design <- model.matrix(~0 + as.factor(merged_meta$Class))
colnames(design) <- c("Control", "HD")

# 2. Fit Linear Model
fit <- lmFit(corrected_expr, design)

# 3. Compute Contrasts (HD vs Control)
contrast <- makeContrasts(HD - Control, levels=design)
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)

# 4. Extract Results (Fold Changes & P-values)
# Methodology C: "select candidate gene list (e.g. top 100-500 DEGs by |log2 FC| and FDR)"
all_res <- topTable(fit2, number=Inf)

# Filter: Adjusted P-value (FDR) < 0.05
# Note: You can adjust this threshold if too few genes pass
deg_list <- all_res[all_res$adj.P.Val < 0.05, ]

# Sort by absolute LogFC to get the most biologically significant ones
deg_list <- deg_list[order(abs(deg_list$logFC), decreasing=TRUE), ]

# Select top 500 candidates as per Methodology C
top_candidates <- head(rownames(deg_list), 500)

message(paste("Number of significant DEGs found:", nrow(deg_list)))
message(paste("Selected top", length(top_candidates), "candidates for Feature Selection."))

save(top_candidates, all_res, deg_list, file="data/processed/dea_results.RData")
write.csv(deg_list, "results/full_dea_results.csv")
message("[C] DEA Complete. Results saved.")