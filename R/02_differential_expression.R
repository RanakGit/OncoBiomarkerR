# R/02_differential_expression.R

library(DESeq2)
library(SummarizedExperiment)

#' Run DESeq2 Pipeline on TCGA LUAD
#'
#' @param se_path Path to the saved SummarizedExperiment .rds file
#' @return A list containing the DESeqDataSet (dds) and the shrunken results
run_deseq2_pipeline <- function(se_path = "data/tcga_luad_subset.rds") {

  message("📂 Loading SummarizedExperiment object...")
  se <- readRDS(se_path)

  stage <- as.character(colData(se)$ajcc_pathologic_stage)
  clean_stage <- ifelse(
    grepl("^Stage I(A|B|C)?$", stage),
    "Stage_I",
    ifelse(stage == "Stage IV", "Stage_IV", NA_character_)
  )

  keep_samples <- !is.na(clean_stage)
  se <- se[, keep_samples]
  clean_stage <- clean_stage[keep_samples]

  if (!all(c("Stage_I", "Stage_IV") %in% clean_stage)) {
    stop("The input data must contain samples from Stage I and Stage IV.")
  }

  colData(se)$tumor_stage <- factor(
    clean_stage,
    levels = c("Stage_I", "Stage_IV")
  )

  if (!"unstranded" %in% assayNames(se)) {
    stop("The input object does not contain the raw 'unstranded' count assay.")
  }
  counts <- assay(se, "unstranded")
  if (any(!is.finite(counts)) || any(counts < 0) || any(counts != floor(counts))) {
    stop("The 'unstranded' assay must contain finite, non-negative integer counts.")
  }

  # 2. Filter out noisy, low-expressed genes to speed up math and increase power
  # Rule: A gene must have at least 10 counts in at least 10 patients
  keep <- rowSums(counts >= 10) >= 10
  counts <- counts[keep, , drop = FALSE]
  message(sprintf("✂️ Filtered low-count genes. %d genes remain for analysis.", nrow(counts)))

  # 3. Create the DESeq2 Object with raw counts
  dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = as.data.frame(colData(se)),
    design = ~ tumor_stage
  )

  # 4. Run the Core Math
  message("🚀 Running DESeq2 Generalized Linear Models (this takes a moment)...")
  dds <- DESeq(dds)

  # 5. Extract Results and Apply Log-Fold Change Shrinkage
  message("📉 Applying 'apeglm' shrinkage to effect sizes...")
  coef_name <- "tumor_stage_Stage_IV_vs_Stage_I"
  res_shrunken <- lfcShrink(dds, coef = coef_name, type = "apeglm")

  return(list(dds = dds, results = res_shrunken))
}
