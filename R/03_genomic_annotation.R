# R/03_genomic_annotation.R

library(SummarizedExperiment)
library(GenomicRanges)

#' Annotate and Extract Significant Biomarkers
#'
#' Filters DESeq2 results for significance and maps them to their
#' physical genomic coordinates using the original SE object.
#'
#' @param res_path Path to the DESeq2 results .rds
#' @param se_path Path to the original SummarizedExperiment .rds
#' @param fdr_cutoff False Discovery Rate threshold (default 0.05)
#' @param lfc_cutoff Log2 Fold Change threshold (default 1.0)
#' @return A GRanges object of significant biomarker genes
extract_significant_ranges <- function(res_path = "data/luad_res.rds",
                                       se_path = "data/tcga_luad_subset.rds",
                                       fdr_cutoff = 0.05,
                                       lfc_cutoff = 1) {

  message("📂 Loading DESeq2 results and original SE object...")
  res <- readRDS(res_path)
  se <- readRDS(se_path)

  # 1. Convert DESeqResults to a standard data frame
  res_df <- as.data.frame(res)

  # 2. Filter for significant genes
  sig_genes <- res_df[
    !is.na(res_df$padj) &
      res_df$padj < fdr_cutoff &
      abs(res_df$log2FoldChange) > lfc_cutoff,
    ,
    drop = FALSE
  ]

  message(sprintf("🎯 Found %d highly significant biomarker genes.", nrow(sig_genes)))

  # 3. Extract the physical genomic coordinates from the original object
  gene_coords <- rowRanges(se)
  coord_indices <- match(rownames(sig_genes), names(gene_coords))
  if (anyNA(coord_indices)) {
    stop("Some significant genes do not have matching genomic coordinates.")
  }

  # 4. Subset the coordinates to only include significant genes
  sig_coords <- gene_coords[coord_indices]

  # 5. Attach the statistical results to the genomic ranges
  mcols(sig_coords)$log2FoldChange <- sig_genes$log2FoldChange
  mcols(sig_coords)$padj <- sig_genes$padj

  return(sig_coords)
}
