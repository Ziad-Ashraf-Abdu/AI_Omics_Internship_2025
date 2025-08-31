# AI\_Omics\_Internship\_2025

This project contains the code and data for different **Modules** of the AI & Omics Research Internship (2025), where we load, clean, preprocess, and analyze biological datasets in R.

---

## 📂 Folder Structure

```
Module_I/
├── raw_data/              # Unmodified source data (patient_info.csv)
├── clean_data/            # Cleaned datasets
├── scripts/               # R scripts
│   └── class_Ib.R         # Main data-cleaning script
├── plots/                 # Figures generated
├── results/               # Analysis outputs, tables, etc.
└── README.md              # This file

Module_II/
├── raw_data/              # DEGs data (downloaded automatically)
├── results/               # Processed DEGs tables
├── scripts/               # R scripts
│   └── classify_DEGs.R    # DEG classification + summary
└── README.md              # This file
```

---

## 🚀 Getting Started

### Module I: Clinical Data Cleaning

1. **Clone your repo** (or download this folder) into your `AI_Omics_Internship_2025/Module_I` directory.

2. **Open in RStudio**

   * Double-click the `Module_I.Rproj` file to launch your RStudio project.
   * This ensures all file-paths are relative to the project root.

3. **Install prerequisites**

   ```r
   install.packages(c("tidyverse"))
   ```

4. **Run the data-cleaning script**

   * Open **scripts/class\_Ib.R** in RStudio and click **Source**.
   * It will download `patient_info.csv`, inspect and clean it, then write out `clean_data/patient_info_clean.csv`.

---

### Module II: DEG Classification

1. **Open in RStudio**

   * Navigate to the `Module_II` directory and open your project.

2. **Run the DEG classification script**

   * Open **scripts/classify\_DEGs.R** and click **Source**.
   * The script will:

     * Download `DEGs_Data_1.csv` and `DEGs_Data_2.csv` from GitHub.
     * Replace missing `padj` values with 1.
     * Classify each gene as **Upregulated**, **Downregulated**, or **Not\_Significant**.
     * Save processed results in the `results/` folder.
     * Print summary counts using `table()` for each dataset and combined results.

---

## 📝 Script Overview

* **Module I (class\_Ib.R)**

  * Directory setup & data download
  * Inspect structure, types, missingness
  * Convert variables to correct formats
  * Add binary `smoking_binary` flag
  * Save `patient_info_clean.csv`

* **Module II (classify\_DEGs.R)**

  * Define `classify_gene()` function
  * Handle missing values (`padj → 1`)
  * Apply classification with a for-loop
  * Add `status` column to each dataset
  * Save processed DEG tables + combined results
  * Print summaries with `table()`

---

## 📋 Prerequisites

* **R** ≥ 4.0
* **RStudio** (to manage the project files)
* **Internet connection** (to download CSVs)
* R packages:

  * `tidyverse` (for data manipulation; Module I)
  * `utils` (base R, for CSV download; Module II)
