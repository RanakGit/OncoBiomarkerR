# R/01_data_ingestion.R

library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)

#' Download and Prepare TCGA LUAD RNA-Seq Data
#'
#' This function queries the GDC for Lung Adenocarcinoma (LUAD) RNA-seq data,
#' downloads a targeted subset (to save laptop memory), and prepares a
#' RangedSummarizedExperiment object.
#'
#' @return A RangedSummarizedExperiment object
download_luad_cohort <- function() {

  message("🔍 Querying GDC for TCGA-LUAD RNA-seq data...")

  # 1. Build the Query
  query <- GDCquery(
    project = "TCGA-LUAD",
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts",
    sample.type = c("Primary Tumor")
  )

  # 2. Extract clinical data and select 20 Stage I and 20 Stage IV patients
  clinical_data <- GDCquery_clinic(project = "TCGA-LUAD", type = "clinical")

  stage_I_barcodes <- clinical_data %>%
    filter(ajcc_pathologic_stage %in% c(
      "Stage I", "Stage IA", "Stage IB", "Stage IC"
    )) %>%
    pull(submitter_id) %>%
    unique() %>%
    head(20)

  stage_IV_barcodes <- clinical_data %>%
    filter(ajcc_pathologic_stage == "Stage IV") %>%
    pull(submitter_id) %>%
    unique() %>%
    head(20)

  target_patients <- c(stage_I_barcodes, stage_IV_barcodes)
  rna_files <- getResults(query)
  target_barcodes <- rna_files %>%
    filter(.data[["cases.submitter_id"]] %in% target_patients) %>%
    pull(cases) %>%
    unique()

  if (length(target_barcodes) == 0) {
    stop("No primary-tumor RNA-seq files matched the selected clinical stages.")
  }

  # 3. Refine Query to just our targeted patients
  query_subset <- GDCquery(
    project = "TCGA-LUAD",
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts",
    barcode = target_barcodes
  )

  # 4. Download the raw files from NCI servers to your laptop
  message("📥 Downloading RNA-seq files (this may take a few minutes)...")
  GDCdownload(query_subset)

  # 5. Prepare the SummarizedExperiment Object
  message("🧬 Assembling the RangedSummarizedExperiment object...")
  se_luad <- GDCprepare(query_subset)

  sample_case_ids <- substr(colnames(se_luad), 1, 12)
  colData(se_luad)$ajcc_pathologic_stage <- clinical_data$ajcc_pathologic_stage[
    match(sample_case_ids, clinical_data$submitter_id)
  ]

  return(se_luad)
}
