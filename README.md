# AI_Omics_Internship_2025

This project contains the code and data for different **Modules** of the AI & Omics Research Internship (2025), where we load, clean, and preprocess a clinical dataset in R.

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
```

---

## 🚀 Getting Started

1. **Clone your repo** (or download this folder) into your `AI_Omics_Internship_2025/Module_I` directory.

2. **Open in RStudio**
   - Double-click the `Module_I.Rproj` file to launch your RStudio project.
   - This ensures all file-paths are relative to the project root.

3. **Install prerequisites**
   
   In the R console, run:
   ```r
   install.packages(c("tidyverse"))
   ```

4. **Run the data-cleaning script**
   
   In RStudio, open **scripts/class_Ib.R** and click **Source**.
   - It will create the `raw_data/`, `clean_data/`, `scripts/`, etc. folders (if they don't exist).
   - It will download `patient_info.csv`, inspect and clean it, then write out `clean_data/patient_info_clean.csv`.

---

## 📝 Script Overview

- **Directory setup**: Creates and sets your working directories.
- **Download & load**: Grabs the CSV from GitHub into `raw_data/`.
- **Inspect**: Uses `str()`, `summary()`, `head()/tail()`, and `colSums(is.na())`.
- **Type conversions**: Factors for categorical vars, numerics for continuous.
- **Binary smoker flag**: New factor `smoking_binary` (0 = No, 1 = Yes).
- **Save cleaned data**: Outputs `patient_info_clean.csv` in `clean_data/`.

---

## 📋 Prerequisites

- **R** ≥ 4.0
- **RStudio** (to manage the project file)
- **Internet connection** (to download the CSV)
- R packages:
  - `tidyverse` (for data manipulation; you can replace base code if preferred)

---

## 👩‍💻 Next Steps

- If you build additional modules, consider refactoring reusable functions into a utility script (e.g. `scripts/utils.R`).
- Add your own analysis or visualization scripts in the `scripts/` folder and save outputs under `plots/` or `results/`.