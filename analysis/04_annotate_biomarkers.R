# analysis/04_annotate_biomarkers.R

# Source the function
source("R/03_genomic_annotation.R")

# Run the extraction
sig_biomarkers <- extract_significant_ranges()

# Sort to find the most up-regulated genes driving Stage IV
top_upregulated <- sort(
  sig_biomarkers[sig_biomarkers$log2FoldChange > 0],
  by = ~ log2FoldChange,
  decreasing = TRUE
)

message("\n🔥 Top 5 Up-regulated Biomarkers (Driving Stage IV):")
print(head(top_upregulated, 5))

# Save this crucial object for the final heatmap
dir.create("data", showWarnings = FALSE)
saveRDS(sig_biomarkers, "data/sig_biomarkers.rds")
message("✅ Significant biomarkers saved to data/sig_biomarkers.rds")
