# AI_Omics_Internship_2025

<div align="center">

![Cairo University](https://img.shields.io/badge/Cairo%20University-Systems%20%26%20Biomedical%20Engineering-darkred?style=for-the-badge)
![License](https://img.shields.io/badge/License-MIT-blue.svg?style=for-the-badge)
![R](https://img.shields.io/badge/R-276DC3?style=for-the-badge&logo=r&logoColor=white)
![Python](https://img.shields.io/badge/Python-3776AB?style=for-the-badge&logo=python&logoColor=white)

</div>

This repository contains code and data for different modules of the AI & Omics Research Internship (2025). Each module implements a small reproducible workflow in R or Python: downloading data, cleaning/preprocessing, quality control, and producing results for reporting.

---

## 📂 Top-level Folder Structure

```
AI_Omics_Internship_2025/
├── Module_I/           # Clinical Data Cleaning
├── Module_II/          # DEG Classification
├── Module_II_3B&3C/    # Preprocessing & Analysis Workflow (GSE79973)
├── Final Project/      # HD Biomarker Discovery (End-to-End Pipeline)
└── README.md           # This file
```

---

## 📁 Module_I — Clinical Data Cleaning

This module focuses on cleaning and standardizing raw clinical patient data.

### Folder Layout

```
Module_I/
├── clean_data/
│   └── patient_info_clean.csv      # Output: Cleaned dataset
├── logs/
│   └── class_Ib.log                # Execution logs
├── raw_data/
│   └── patient_info.csv            # Input: Raw clinical data
├── scripts/
│   └── Class_Ib.R                  # Main cleaning script
└── Module_I.Rproj
```

### Quick Start

1. Open `Module_I/Module_I.Rproj` in RStudio
2. Install prerequisites:
   ```r
   install.packages(c("tidyverse"))
   ```
3. Source `Module_I/scripts/Class_Ib.R`
4. **Outputs:** The cleaned CSV will be saved to `Module_I/clean_data/`

---

## 📁 Module_II — DEG Classification

This module handles the processing and classification of Differentially Expressed Genes (DEGs) from multiple datasets.

### Folder Layout

```
Module_II/
├── raw_data/
│   ├── DEGs_Data_1.csv             # Input dataset 1
│   └── DEGs_Data_2.csv             # Input dataset 2
├── results/
│   ├── Combined_DEGs_Results.csv   # Merged results
│   ├── Processed_DEGs_Data_1.csv
│   └── Processed_DEGs_Data_2.csv
├── scripts/
│   └── Class_2b.R                  # Main processing script
└── Module_II.Rproj
```

### Quick Start

1. Open `Module_II/Module_II.Rproj` in RStudio
2. Source `Module_II/scripts/Class_2b.R`
3. **Outputs:** Processed and combined results written to `Module_II/results/`

---

## 📁 Module_II_3B&3C

This module implements the full preprocessing and differential expression workflow for GEO series GSE79973 (Affymetrix GPL570).

### Folder Layout

```
Module_II_3B&3C/
├── Scripts/
│   ├── Class_3C.R                  # Main analysis script
│   └── Class_eB.R                  # Supplemental analysis
├── Results/
│   ├── DGE_Results_GSE79973.csv    # Differential Expression results
│   ├── normalized_expression_matrix.csv
│   ├── filtered_expression_matrix.csv
│   ├── processing_summary.txt
│   └── QC_Normalized_Data/         # ArrayQualityMetrics HTML reports
├── Result_Plots/
│   ├── Heatmap_Top50_GSE79973.png
│   └── Volcano_Plot_GSE79973.png
└── Module_II_3B.Rproj
```

### Quick Start

1. Source `Module_II_3B&3C/Scripts/Class_3C.R`
2. **Outputs:** Normalized matrices in `Results/` and visualizations in `Result_Plots/`

---

## 📁 Final Project — HD Biomarker Discovery

This module contains the complete pipeline for biomarker discovery for Huntington's Disease (HD), utilizing both R (for bioinformatics analysis) and Python (for Machine Learning inference).

### Folder Layout

```
Final Project/
├── data/                       # Raw and clean metadata/counts
├── figures/                    # Volcano plots, Heatmaps
├── hd_models/                  # Pre-trained Python models (RF, SVM, Scalers)
├── results/                    # DEA results and predictions
├── scripts/                    # R Bioinformatics Pipeline
│   ├── data_accuision01.R
│   ├── preprocessing02.R
│   ├── differential_expression03.R
│   ├── feature_selection04.R
│   └── verification05.R
├── inference.py                # Python Inference Script
├── requirements.txt            # Python dependencies
└── Final Project.Rproj
```

---

## 🔬 R Bioinformatics Pipeline (Custom Data Instructions)

The default pipeline uses **GSE129473** (Human Caudate Nucleus) for validation, and **GSE64810** for training.

### How to Use Your Own Dataset

> **Note:** Your dataset must be RNA-seq data (not protein).

To run this pipeline on a new dataset, you only need to modify (i.e., replace) the **GSE129473** with your target dataset ID. The subsequent scripts are modular and will adapt to the new data flow.

#### Pipeline Steps

1. **Data Acquisition** (`scripts/data_accuision01.R`):
   - Ensure the script saves the raw counts and metadata to the standard paths expected by the pipeline

2. **Differential Expression** (`scripts/differential_expression03.R`):
   - Run this script to perform DE analysis on the normalized data

3. **Feature Extraction** (`scripts/feature_selection04.R`):
   - Run this script to extract the significant features (biomarkers) for the models

4. **Verification** (`scripts/verification05.R`):
   - Run this validation script last. It does not perform analysis but checks if the previous processes ran smoothly and generated the correct outputs

---

## 🤖 Python Inference Pipeline

The project includes a weighted ensemble classifier (SVM + Random Forest) to predict HD status from processed data.

### 1. Environment Setup

You need a Python environment to run the inference script. You can create one using conda or venv and install the dependencies from the included `requirements.txt`.

#### Option A: Using Conda (Recommended)

```bash
# Create a new environment
conda create -n ai_omics python=3.9 -y
conda activate ai_omics

# Install dependencies
pip install -r "Final Project/requirements.txt"
```

#### Option B: Using Python venv

```bash
# Create and activate virtual environment
python -m venv ai_omics_env

# Windows
ai_omics_env\Scripts\activate

# Mac/Linux
source ai_omics_env/bin/activate

# Install dependencies
pip install -r "Final Project/requirements.txt"
```

### 2. Running Inference

Once your environment is active:

```bash
# Navigate to the Final Project folder
cd "Final Project"

# Run inference on a CSV file
python inference.py "path/to/your_input_data.csv" --output "my_predictions.csv"
```

---

## 🎓 Course Information

This project was developed as part of the **AI & Omics Research Internship 2025** program.

<div align="center">

### 📚 Organized By

**AI, Biotechnology & Bioinformatics Learning Hub**

[![GitHub](https://img.shields.io/badge/GitHub-100000?style=for-the-badge&logo=github&logoColor=white)](https://github.com/AI-Biotechnology-Bioinformatics)
[![Course Materials](https://img.shields.io/badge/Course%20Materials-View%20Here-success?style=for-the-badge)](https://github.com/AI-Biotechnology-Bioinformatics/AI_and_Omics_Research_Internship_2025)

*Empowering the next generation of bioinformatics researchers through hands-on learning*

</div>

### 🌐 Related Courses & Resources

Check out other amazing courses from the same organization:

- 🧬 [**R Crash Course 2025**](https://github.com/AI-Biotechnology-Bioinformatics/R_Crash_Course-2025) - Beginner-friendly R programming for biology students
- 📊 [**Microarray Data Analysis in R**](https://github.com/AI-Biotechnology-Bioinformatics/Microarray_Series_R) - Comprehensive microarray analysis workflows
- 📈 [**Logistic Regression in R**](https://github.com/AI-Biotechnology-Bioinformatics/Logistic_Regression_R) - Statistical modeling for bioinformatics

---

## 📝 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## 👥 Meet the Team

<table>
  <tr>
    <td align="center">
      <a href="https://www.linkedin.com/in/ziad-mohamed-2a956b282">
        <img src="https://img.shields.io/badge/LinkedIn-0077B5?style=for-the-badge&logo=linkedin&logoColor=white" alt="LinkedIn"/>
      </a>
      <br />
      <sub><b>Zeyad Ashraf Abdo Muhammed</b></sub>
      <br />
      <sub>🧬 Bioinformatics Pipeline</sub>
      <br />
      <a href="mailto:ziad.mohamed04@eng-st.cu.edu.eg">📧 Email</a>
    </td>
    <td align="center">
      <a href="https://www.linkedin.com/in/ahmed-abdelsalam1">
        <img src="https://img.shields.io/badge/LinkedIn-0077B5?style=for-the-badge&logo=linkedin&logoColor=white" alt="LinkedIn"/>
      </a>
      <br />
      <sub><b>Ahmed Muhammed Abdelsalam</b></sub>
      <br />
      <sub>📊 Data Analysis & Visualization</sub>
      <br />
      <a href="mailto:ahmed.mohamed0410@eng-st.cu.edu.eg">📧 Email</a>
    </td>
    <td align="center">
      <a href="https://www.linkedin.com/in/rahmafathy105">
        <img src="https://img.shields.io/badge/LinkedIn-0077B5?style=for-the-badge&logo=linkedin&logoColor=white" alt="LinkedIn"/>
      </a>
      <br />
      <sub><b>Rahma Fathy</b></sub>
      <br />
      <sub>🤖 Machine Learning Models</sub>
      <br />
      <a href="mailto:rahma.edress04@eng-st.cu.edu.eg">📧 Email</a>
    </td>
  </tr>
</table>

---

## 🌟 About This Project

This internship project represents a comprehensive exploration of AI applications in omics research, from basic data cleaning to advanced machine learning-based biomarker discovery. We've combined traditional bioinformatics approaches with modern ML techniques to create a robust pipeline for Huntington's Disease research.

### 🎯 What We've Accomplished

- ✅ Built end-to-end reproducible bioinformatics workflows
- ✅ Implemented differential gene expression analysis pipelines
- ✅ Developed ML models for disease prediction (SVM + Random Forest ensemble)
- ✅ Created modular, reusable code for omics research
- ✅ Established best practices for data preprocessing and quality control

---

## 💬 Get In Touch

We're always excited to discuss bioinformatics, machine learning, and omics research! Whether you have questions about the pipeline, suggestions for improvements, or potential collaboration ideas:

- 💡 **Open an Issue** for technical questions or bug reports
- 🤝 **Reach out directly** via email or LinkedIn for collaboration opportunities
- ⭐ **Star this repo** if you find it useful!
- 🍴 **Fork and contribute** - we welcome pull requests!

### 🏛️ Institutional Affiliation

**Cairo University**  
Faculty of Engineering  
Department of Systems and Biomedical Engineering

---

<div align="center">

**Made with ❤️ by the Systems & Biomedical Engineering Team**

[![Cairo University](https://img.shields.io/badge/Cairo%20University-Visit%20Website-darkred?style=flat-square)](https://eng.cu.edu.eg)

</div>
