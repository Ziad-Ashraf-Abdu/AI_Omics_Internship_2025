# AI_Omics_Internship_2025

This repository contains code and data for different **Modules** of the AI & Omics Research Internship (2025). Each module implements a small reproducible workflow in R: downloading data, cleaning / preprocessing, quality control and producing results for reporting.

---

## 📂 Top-level Folder Structure (updated)

```
Module_I/
Module_II/
Module_II_EB/
README.md                          # this file (updated)
```

---

## 📁 Module_I — Clinical Data Cleaning 

Detailed structure: see `Module_I/README.md`.

Quick start:

* Open `Module_I/Module_I.Rproj` in RStudio
* Install prerequisites: `install.packages(c("tidyverse"))`
* Source `Module_I/scripts/class_Ib.R`

Outputs: cleaned CSV in `Module_I/clean_data/` and plots in `Module_I/plots/`.

---

## 📁 Module_II — DEG Classification

Detailed structure: see `Module_II/README.md`.

Quick start:

* Open `Module_II` project in RStudio
* Source `Module_II/scripts/classify_DEGs.R`
* Outputs written to `Module_II/results/`.

---

## 📁 Module_II_3B

This module implements the full preprocessing workflow for GEO series **GSE79973** (Affymetrix GPL570) as required by the AI & Omics Internship — Module II / Assignment #4.

### Folder layout

```
Module_II_3B/
├── raw_data/                   # optional: place downloaded CELs here (not required)
├── scripts/
│   └── Class_3B.R   # main R script (runs end-to-end)
├── Results/                    # output folder created by the script
│   ├── normalized_expression_matrix.csv
│   ├── filtered_expression_matrix.csv
│   ├── probe_median_hist.png
│   ├── boxplot_normalized.png
│   ├── PCA_normalized.png
│   ├── QC_Normalized_Data/     # arrayQualityMetrics HTML report
│   └── processing_summary.txt
└── README.md
```

### Required R packages

The script will attempt to install missing packages. It assumes R ≥ 4.0. Main packages used:
* tidyverse (for data manipulation; Module I)
* utils (base R, for CSV download; Module II)
* Bioconductor: `GEOquery`, `affy`, `arrayQualityMetrics`, `genefilter`, `AnnotationDbi`, `hgu133plus2.db` (annotation)
* CRAN: `dplyr`, `matrixStats`
**Fix annotation mapping** (important): the platform returned by GEO is `GPL570`; the correct Bioconductor annotation package is `hgu133plus2.db`. The script includes an automatic mapping lookup and will install/load `hgu133plus2.db` if needed.


