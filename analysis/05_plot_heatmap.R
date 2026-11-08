# analysis/05_plot_heatmap.R

# Source the function
source("R/04_visualization.R")

# Generate the heatmap object
ht <- plot_clinical_heatmap()

if (interactive()) {
  # Draw it to the active VS Code R graphics device
  draw(ht)
} else {
  # Save a portable copy when running without an interactive graphics device
  dir.create("data", showWarnings = FALSE)
  grDevices::pdf("data/clinical_biomarker_heatmap.pdf", width = 11, height = 9)
  ComplexHeatmap::draw(ht)
  grDevices::dev.off()
  message("✅ Heatmap saved to data/clinical_biomarker_heatmap.pdf")
}
