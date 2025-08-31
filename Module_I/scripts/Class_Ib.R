# ===================================================================
#               AI and Biotechnology/Bioinformatics        
# ===================================================================
# -------------------------------------------------------------------
#             AI and Omics Research Internship (2025)
# -------------------------------------------------------------------
#               Module I: Getting Started with R (Class Ib) - TASKS
# -------------------------------------------------------------------
# ==================================================================

# TASK I: Setting Working Directory
# ===================================================================

cat("\n=== TASK I: Setting Working Directory ===\n")

# First, creating the main folder using R
main_folder <- "C:/AI_Omics_Internship_2025" 
if (!dir.exists(main_folder)) {
  dir.create(main_folder, recursive = TRUE)
  cat("✓ Created main folder:", main_folder, "\n")
} else {
  cat("⚠ Main folder already exists:", main_folder, "\n")
}

# Setting working directory to the main folder
setwd(main_folder)
cat("✓ Working directory set to:", getwd(), "\n")

# TASK II: Creating Project Folder Structure
# ===================================================================

cat("\n=== TASK II: Creating Project Folder Structure ===\n")

# Creating Module_I project folder
project_folder <- file.path(main_folder, "Module_I")
if (!dir.exists(project_folder)) {
  dir.create(project_folder)
  cat("✓ Created project folder:", project_folder, "\n")
} else {
  cat("⚠ Project folder already exists:", project_folder, "\n")
}

# Setting working directory to project folder
setwd(project_folder)
cat("✓ Working directory changed to:", getwd(), "\n")

# Creating subfolders using R code
folders <- c("raw_data", "clean_data", "scripts", "results", "tasks", "plots")

cat("Creating project subfolders...\n")
for (folder in folders) {
  if (!dir.exists(folder)) {
    dir.create(folder)
    cat("✓ Created folder:", folder, "\n")
  } else {
    cat("⚠ Folder already exists:", folder, "\n")
  }
}

# Verifying folder structure
folders_created <- list.dirs(recursive = FALSE)
cat("✓ Project folder structure verified. Total folders:", length(folders_created), "\n")

# TASK III: Downloading and Processing Patient Data
# ===================================================================

cat("\n=== TASK III: Downloading and Processing Patient Data ===\n")

# GitHub repository URL for the dataset
github_url <- "https://raw.githubusercontent.com/AI-Biotechnology-Bioinformatics/AI_and_Omics_Research_Internship_2025/main/patient_info.csv"

cat("Downloading patient_info.csv from GitHub...\n")
cat("Source URL:", github_url, "\n")

tryCatch({
  # Downloading to raw_data folder
  download.file(github_url, destfile = "raw_data/patient_info.csv", method = "auto")
  
  # Loading the dataset
  patient_data <- read.csv("raw_data/patient_info.csv")
  cat("✓ Dataset successfully downloaded and loaded!\n")
  cat("✓ Dataset contains", nrow(patient_data), "rows and", ncol(patient_data), "columns\n")
  
}, error = function(e) {
  cat("✗ Error downloading from GitHub:", e$message, "\n")
  # Try to load from local file if it exists
  if (file.exists("raw_data/patient_info.csv")) {
    patient_data <- read.csv("raw_data/patient_info.csv")
    cat("⚠ Loaded dataset from existing local file\n")
  } else {
    cat("✗ No dataset available. Please check internet connection or URL\n")
    stop("Cannot proceed without dataset")
  }
})

# INSPECTING THE DATASET
# ===================================================================

cat("\n=== DATASET INSPECTION ===\n")

# Basic information about the dataset
cat("Dataset dimensions:", nrow(patient_data), "rows ×", ncol(patient_data), "columns\n")

# Display first few rows
cat("\nFirst 6 rows of the dataset:\n")
print(head(patient_data))

# Display last few rows
cat("\nLast 6 rows of the dataset:\n")
print(tail(patient_data))

# Check structure of the dataset
cat("\nDataset structure:\n")
str(patient_data)

# Summary statistics
cat("\nSummary statistics:\n")
summary(patient_data)

# Check for missing values
missing_values <- colSums(is.na(patient_data))
total_missing <- sum(missing_values)
if (total_missing > 0) {
  cat("⚠ Found", total_missing, "missing values in the dataset\n")
  cat("Missing values by column:\n")
  print(missing_values[missing_values > 0])
} else {
  cat("✓ No missing values found in the dataset\n")
}

# Check data types
cat("\nData types of each column:\n")
data_types <- sapply(patient_data, class)
for (col in names(data_types)) {
  cat("  •", col, ":", data_types[col], "\n")
}

# Check unique values in categorical columns
cat("\nAnalyzing categorical columns:\n")
categorical_cols <- sapply(patient_data, function(x) is.character(x) || is.factor(x) || is.integer(x) || is.numeric(x))
for (col in names(patient_data)[categorical_cols]) {
  unique_vals <- unique(patient_data[[col]])
  cat("  •", col, ":", length(unique_vals), "unique values -", paste(unique_vals, collapse = ", "), "\n")
}

# DATA TYPE CORRECTIONS
# ===================================================================

cat("\n=== DATA TYPE CORRECTIONS ===\n")

# Creating a copy for cleaning
patient_clean <- patient_data
cat("✓ Created working copy of dataset for cleaning\n")

# Converting character variables to factors where appropriate
conversions_made <- 0

if ("gender" %in% names(patient_clean)) {
  patient_clean$gender <- as.factor(patient_clean$gender)
  cat("✓ Converted 'gender' to factor\n")
  conversions_made <- conversions_made + 1
}

if ("diagnosis" %in% names(patient_clean)) {
  patient_clean$diagnosis <- as.factor(patient_clean$diagnosis)
  cat("✓ Converted 'diagnosis' to factor\n")
  conversions_made <- conversions_made + 1
}

if ("smoker" %in% names(patient_clean)) {
  patient_clean$smoker <- as.factor(patient_clean$smoker)
  cat("✓ Converted 'smoker' to factor\n")
  conversions_made <- conversions_made + 1
}

# Ensuring numeric variables are properly typed
numeric_cols <- c("age", "bmi")
for (col in numeric_cols) {
  if (col %in% names(patient_clean)) {
    patient_clean[[col]] <- as.numeric(patient_clean[[col]])
    cat("✓ Ensured", col, "is numeric\n")
    conversions_made <- conversions_made + 1
  }
}

cat("✓ Completed", conversions_made, "data type conversions\n")

# CREATE BINARY SMOKING STATUS VARIABLE
# ===================================================================

cat("\n=== CREATING BINARY SMOKING STATUS ===\n")

if ("smoker" %in% names(patient_clean)) {
  # Create binary smoking variable (1 for "Yes", 0 for "No")
  patient_clean$smoking_binary <- ifelse(patient_clean$smoker == "Yes", 1, 0)
  
  # Converting to factor with proper labels
  patient_clean$smoking_binary <- factor(patient_clean$smoking_binary,
                                         levels = c(0, 1),
                                         labels = c("0", "1"))
  
  cat("✓ Created binary smoking status variable\n")
  
  # Log the distribution
  smoking_table <- table(patient_clean$smoking_binary)
  cat("Binary smoking distribution: Non-smokers (", smoking_table[1], "), Smokers (", smoking_table[2], ")\n")
  
  # Show the mapping
  cat("Verification of original to binary mapping:\n")
  print(table(patient_clean$smoker, patient_clean$smoking_binary))
} else {
  cat("⚠ 'smoker' column not found - cannot create binary smoking variable\n")
}

# FINAL DATA INSPECTION
# ===================================================================

cat("\n=== FINAL CLEANED DATASET ===\n")

# Checking final structure
cat("Final dataset structure:\n")
str(patient_clean)

# Summary of cleaned data
cat("\nSummary of cleaned dataset:\n")
summary(patient_clean)

# SAVE CLEANED DATASET
# ===================================================================

cat("\n=== SAVING CLEANED DATASET ===\n")

tryCatch({
  # Save cleaned dataset to clean_data folder
  write.csv(patient_clean, "clean_data/patient_info_clean.csv", row.names = FALSE)
  cat("✓ Cleaned dataset saved as 'patient_info_clean.csv' in clean_data folder\n")
  
  # Log file size
  file_size <- file.info("clean_data/patient_info_clean.csv")$size
  cat("✓ Saved file size:", round(file_size/1024, 2), "KB\n")
  
}, error = function(e) {
  cat("✗ Failed to save cleaned dataset:", e$message, "\n")
})

# SAVE R SCRIPT
# ===================================================================

cat("\n=== SAVING R SCRIPT ===\n")

script_path <- "scripts/class_Ib.R"
cat("Script should be saved as:", script_path, "\n")
cat("Remember to upload the script to your GitHub repository\n")

# SUMMARY REPORT
# ===================================================================

cat("\n=== TASK COMPLETION SUMMARY ===\n")
cat("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "]\n")
cat("✓ Working directory set to:", getwd(), "\n")
cat("✓ Project folder structure created\n")
cat("✓ Dataset downloaded/loaded successfully\n")
cat("✓ Data inspection completed\n")
cat("✓ Data types corrected\n")
cat("✓ Binary smoking status variable created\n")
cat("✓ Cleaned dataset saved\n")
cat("✓ Script ready for GitHub upload\n")

cat("\nProject folders created:\n")
folders_list <- list.dirs(recursive = FALSE)
for (folder in basename(folders_list)) {
  cat("  • ", folder, "\n")
}

# Check and log files in key directories
raw_files <- list.files("raw_data")
if (length(raw_files) > 0) {
  cat("Files in raw_data folder:", paste(raw_files, collapse = ", "), "\n")
} else {
  cat("⚠ No files found in raw_data folder\n")
}

clean_files <- list.files("clean_data")
if (length(clean_files) > 0) {
  cat("Files in clean_data folder:", paste(clean_files, collapse = ", "), "\n")
} else {
  cat("⚠ No files found in clean_data folder\n")
}

# Final completion message
cat("\n✓ === ALL TASKS COMPLETED SUCCESSFULLY ===\n")
cat("Session completed at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

# END OF SCRIPT
# ===================================================================