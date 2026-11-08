# analysis/02_fetch_tcga.R

# Source the function we just wrote
source("R/01_data_ingestion.R")

# Run the pipeline
se_luad <- download_luad_cohort()

# Inspect the object you just downloaded
print(se_luad)
table(colData(se_luad)$ajcc_pathologic_stage)

# Save it as an R Data Structure (.rds) so we never have to download it again
dir.create("data", showWarnings = FALSE)
saveRDS(se_luad, file = "data/tcga_luad_subset.rds")
message("✅ Data successfully saved to data/tcga_luad_subset.rds")
