# analysis/03_run_deseq2.R

# Source the function
source("R/02_differential_expression.R")

# Run the pipeline
de_output <- run_deseq2_pipeline()

# Save the outputs for the next steps
dir.create("data", showWarnings = FALSE)
saveRDS(de_output$dds, "data/luad_dds.rds")
saveRDS(de_output$results, "data/luad_res.rds")

# Print a statistical summary to the terminal
message("\n✅ Pipeline Complete! Here is the summary of significant biomarkers:")
summary(de_output$results)
