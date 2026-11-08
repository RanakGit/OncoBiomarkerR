# analysis/01_synthetic_data_challenge.R

# Load the core Bioconductor packages
library(SummarizedExperiment)
library(GenomicRanges)

# ---------------------------------------------------------
# 1. The Counts Matrix (assays)
# ---------------------------------------------------------
# Goal: A matrix with 100 rows (genes) and 6 columns (patients).
# Hint: Use matrix() combined with rpois(600, lambda=50) to generate random RNA-seq-like counts.

set.seed(1)
gene_ids <- paste0("G", seq_len(100))
patient_ids <- paste0("P", seq_len(6))

counts_matrix <- matrix(
  rpois(600, lambda = 50),
  nrow = 100,
  ncol = 6,
  dimnames = list(gene_ids, patient_ids)
)

# ---------------------------------------------------------
# 2. The Clinical Metadata (colData)
# ---------------------------------------------------------
# Goal: A DataFrame containing patient IDs and their drug response.
# Hint: Use DataFrame() (not base R's data.frame) with 6 rows.
# Create a 'patient_id' column ("P1", "P2", ..., "P6") and a
# 'treatment_response' column (3 "Sensitive", 3 "Resistant").

clinical_data <- DataFrame(
  patient_id = patient_ids,
  treatment_response = c(
    rep("Sensitive", 3),
    rep("Resistant", 3)
  ),
  row.names = patient_ids
)

# ---------------------------------------------------------
# 3. The Gene Coordinates (rowRanges)
# ---------------------------------------------------------
# Goal: Genomic coordinates for your 100 genes.
# Hint: Use GRanges(). Set seqnames to "chr1" for all 100,
# and use IRanges(start = 1:100, width = 1000) for dummy coordinates.

gene_coords <- GRanges(
  seqnames = rep("chr1", 100),
  ranges = IRanges(start = 1:100, width = 1000),
  names = gene_ids
)

# ---------------------------------------------------------
# 4. The Assembly
# ---------------------------------------------------------
# Goal: Combine the three components into a single RangedSummarizedExperiment.
# Note: Ensure you pass your counts_matrix inside a list, e.g., list(counts = counts_matrix).

se <- SummarizedExperiment(
  assays = list(counts = counts_matrix),
  rowRanges = gene_coords,
  colData = clinical_data
)

# ---------------------------------------------------------
# 5. The Verification Test
# ---------------------------------------------------------
# If you built the object correctly, this subsetting command will run
# without errors and print the dimensions: 100 rows, 3 columns.

resistant_se <- se[, se$treatment_response == "Resistant"]
print(dim(resistant_se))