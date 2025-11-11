# -----------------------------------------------------------------------------
# Script: 04_feature_selection.R (Corresponds to Methodology Section D)
# Purpose: Refine features using Recursive Feature Elimination (RFE).
# -----------------------------------------------------------------------------

# --- Automatic Package Installation ---
if (!require("caret")) install.packages("caret")
if (!require("randomForest")) install.packages("randomForest")
if (!require("doParallel")) install.packages("doParallel")
# --------------------------------------

library(caret)
library(randomForest)
library(doParallel) # For parallel processing speed

load("data/processed/normalized_data.RData")
load("data/processed/dea_results.RData")

message("[D] Starting Feature Selection...")

# 1. Prepare Data
# Filter the main expression matrix to keep only the Top 500 candidates from Step C
# Methodology D: "refine the top 100 genes... from DEA results"
candidate_expr <- corrected_expr[rownames(corrected_expr) %in% top_candidates, ]

# Transpose for ML (Rows=Samples, Cols=Genes)
ml_data <- as.data.frame(t(candidate_expr))
ml_data$Class <- as.factor(ifelse(merged_meta$Class == 1, "HD", "Control"))

# 2. Recursive Feature Elimination (RFE)
# Methodology D: "use recursive feature elimination (RFE) or Random Forest importance"

# Set up RFE control
# We will use Random Forest as the RFE function (rfFuncs)
control <- rfeControl(functions=rfFuncs, method="cv", number=5)

# Parallel processing (optional but recommended for RFE)
cl <- makePSOCKcluster(detectCores() - 1)
registerDoParallel(cl)

message("Running RFE (this may take a few minutes)...")
# We check performance for subsets of sizes: 10, 20, 50, 100 genes
results_rfe <- rfe(x = ml_data[, -ncol(ml_data)], 
                   y = ml_data$Class, 
                   sizes = c(10, 20, 50, 100), 
                   rfeControl = control)

stopCluster(cl)

# 3. Extract Top Biomarkers
# Methodology D: "refine the top 100 genes"
top_genes <- predictors(results_rfe)
final_biomarkers <- head(top_genes, 100) # Ensure we cap at top 100 if RFE keeps more

message(paste("RFE Selected", length(final_biomarkers), "optimal features."))
print(head(final_biomarkers))

# 4. Save Final Training Data
# Subset original data to just these final biomarkers + Class
final_training_data <- ml_data[, c(final_biomarkers, "Class")]

write.csv(final_training_data, "results/final_biomarker_dataset.csv", row.names=TRUE)
save(results_rfe, final_biomarkers, file="results/rfe_results.RData")

message("[D] Feature Selection Complete. Final dataset ready: results/final_biomarker_dataset.csv")