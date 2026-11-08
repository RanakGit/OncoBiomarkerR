# 🧬 OncoBiomarkerR: Precision Oncology Transcriptomic Pipeline

[![Live Clinical Report](https://img.shields.io/badge/View-Live_Clinical_Report-2980b9?style=for-the-badge)](https://[your-username].github.io/OncoBiomarkerR/)
*(Click the badge above to view the interactive Quarto clinical report and heatmap)*

## Executive Summary
**OncoBiomarkerR** is an end-to-end R/Bioconductor pipeline designed to identify transcriptomic drivers of metastatic progression in Lung Adenocarcinoma (TCGA-LUAD). It programmatically ingests raw RNA-seq counts and clinical metadata, models differential expression using negative binomial generalized linear models, and maps statistically significant biomarkers to their physical genomic coordinates.

## 🛠️ Core Architecture & Tech Stack
This repository is engineered using industry-standard Bioconductor S4 object classes to ensure strict alignment between clinical metadata, count matrices, and genomic coordinates.

*   **Data Structures:** `RangedSummarizedExperiment`, `GenomicRanges`
*   **Statistical Modeling:** `DESeq2` (with `apeglm` LFC shrinkage)
*   **Data Ingestion:** `TCGAbiolinks` (GDC API integration)
*   **Visualization & Reporting:** `ComplexHeatmap`, `Quarto`
*   **Environment Management:** `renv`

## 🔬 Key Biological Findings
The pipeline successfully stratifies Stage I (localized) vs. Stage IV (metastatic) TCGA-LUAD patients. Genomic mapping isolated **123 highly significant biomarkers** (FDR < 0.05, |log2FC| > 1.0). 

Top up-regulated drivers of Stage IV progression include classic X-linked Cancer-Testis Antigens:
*   **COL2A1** (log2FC: 6.50)
*   **MAGEA12** (log2FC: 5.80)
*   **MAGEA3** (log2FC: 5.70)

## 🚀 How to Run Locally

**1. Clone the repository and restore the environment:**
```bash
git clone [https://github.com/](https://github.com/)[your-username]/OncoBiomarkerR.git
cd OncoBiomarkerR
Rscript -e "renv::restore()"
