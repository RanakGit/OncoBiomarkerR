# R/04_visualization.R

library(ComplexHeatmap)
library(DESeq2)
library(circlize)

#' Generate a Clinical Biomarker Heatmap
#'
#' @param dds_path Path to the DESeq2 object
#' @param sig_path Path to the significant GRanges object
#' @return A ComplexHeatmap object
plot_clinical_heatmap <- function(dds_path = "data/luad_dds.rds",
                                  sig_path = "data/sig_biomarkers.rds") {

  message("🎨 Building the clinical heatmap...")
  dds <- readRDS(dds_path)
  sig_genes <- readRDS(sig_path)

  # 1. Extract the normalized RNA-seq counts
  norm_counts <- counts(dds, normalized = TRUE)

  # 2. Subset the matrix to the significant genes
  gene_indices <- match(names(sig_genes), rownames(norm_counts))
  if (anyNA(gene_indices)) {
    stop("Some significant biomarker genes are missing from the DESeq2 object.")
  }
  mat <- norm_counts[gene_indices, , drop = FALSE]

  # 3. Apply Z-score scaling across the rows (genes)
  mat_scaled <- t(scale(t(mat)))
  if (any(!is.finite(mat_scaled))) {
    stop("Z-score scaling produced non-finite values; check for constant-expression genes.")
  }

  # 4. Extract clinical metadata for the top annotation bar
  clinical_info <- as.data.frame(colData(dds))
  if (!"tumor_stage" %in% names(clinical_info)) {
    stop("The DESeq2 object is missing the 'tumor_stage' clinical metadata.")
  }

  # 5. Build the visual top annotation
  top_anno <- HeatmapAnnotation(
    Stage = clinical_info$tumor_stage,
    col = list(Stage = c("Stage_I" = "#2980b9", "Stage_IV" = "#c0392b"))
  )

  # 6. Construct the Heatmap
  ht <- Heatmap(
    mat_scaled,
    name = "Expression\nZ-score",
    top_annotation = top_anno,
    show_row_names = FALSE,
    show_column_names = FALSE,
    cluster_columns = TRUE,
    cluster_rows = TRUE,
    row_title = sprintf("%d Significant Biomarkers", nrow(mat)),
    column_title = "TCGA-LUAD Patients (Stage I vs IV)"
  )

  return(ht)
}
