# OncoBiomarkerR

[![View the clinical report](https://img.shields.io/badge/Report-GitHub%20Pages-2ea44f)](https://ranakgit.github.io/OncoBiomarkerR/)

**[Open the standalone TCGA-LUAD clinical report](https://ranakgit.github.io/OncoBiomarkerR/)**

An exploratory Bioconductor workflow that downloads a small TCGA-LUAD RNA-seq
cohort, compares Stage I with Stage IV using DESeq2, maps significant results
to genomic coordinates, and renders a clinical heatmap.

## Analysis overview

1. Select 20 Stage I and 20 Stage IV primary-tumor cases from GDC clinical
   metadata and retrieve their STAR-count RNA-seq files.
2. Filter low-count genes and run DESeq2 with the design `~ tumor_stage`;
   use `apeglm` to shrink log-fold changes.
3. Extract genes with adjusted p-value below 0.05 and absolute shrunken
   log2-fold change above 1, then attach the statistics to their `GRanges`.
4. Cluster the samples using normalized expression Z-scores for those genes
   and show tumor stage in the heatmap annotation.
5. Render a self-contained Quarto HTML report for browser viewing.

**No batch correction or batch covariate was used.** The DESeq2 model includes
tumor stage only. Results are exploratory associations from a small subset;
they are not clinically validated and should not be interpreted as evidence of
causation or treatment response.

## Results from the included report

- 40 samples: 20 Stage I and 20 Stage IV.
- 23,223 genes retained for differential-expression analysis.
- 123 genes met the stricter annotation filter (`padj < 0.05` and
  `abs(log2FoldChange) > 1`).
- Top Stage IV-upregulated genes include **COL2A1**, **MAGEA12**, **MAGEA3**,
  **SEZ6**, and **CSAG1**.
- The heatmap uses unsupervised column clustering. The stage groups do not
  form perfectly separate clusters in this cohort.

The DESeq2 summary also reports 277 upregulated and 222 downregulated genes at
adjusted p-value below 0.1; this is a different, less stringent threshold
than the 123-gene annotation filter.

## Repository layout

```text
R/                         Reusable pipeline functions
analysis/                  Milestone scripts and Quarto report source
index.html                 Standalone report served by GitHub Pages
```

Large downloaded files, serialized analysis objects, and the GDC cache are
intentionally excluded from Git. Running the pipeline downloads the inputs and
creates these local artifacts under `data/` and `GDCdata/`.

## Reproduce the analysis

Requirements: R 4.5 or compatible, internet access to the NCI GDC, and Quarto
with its CLI available on `PATH`.

Install the required R packages:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c(
  "TCGAbiolinks", "SummarizedExperiment", "dplyr",
  "DESeq2", "apeglm", "ComplexHeatmap", "circlize"
))
```

From the repository root, run the milestones in order:

```powershell
Rscript analysis/02_fetch_tcga.R
Rscript analysis/03_run_deseq2.R
Rscript analysis/04_annotate_biomarkers.R
Rscript analysis/05_plot_heatmap.R
quarto render analysis/clinical_report.qmd
Copy-Item analysis/clinical_report.html index.html
```

The ingestion step requires GDC availability and can take time. The heatmap
script saves `data/clinical_biomarker_heatmap.pdf` in non-interactive runs.
Rendering the Quarto source produces `analysis/clinical_report.html`; the
root-level `index.html` is the self-contained page deployed by GitHub Pages.

## Main artifacts

- [Live clinical report](https://ranakgit.github.io/OncoBiomarkerR/)
- [Quarto report source](analysis/clinical_report.qmd)
- [Ingestion function](R/01_data_ingestion.R)
- [Differential-expression function](R/02_differential_expression.R)
- [Genomic annotation function](R/03_genomic_annotation.R)
- [Heatmap function](R/04_visualization.R)
