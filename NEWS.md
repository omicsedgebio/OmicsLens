# OmicsLens 0.1.0

## New features

* `omicslens_load()` — unified loader for RNA-Seq, WGS (MAF/binary), and methylation beta matrices.
* `omicslens_preprocess()` — one-function preprocessing: DESeq2 VST for RNA, MAF-frequency filtering for variants, M-value conversion and variance filtering for methylation.
* `omicslens_integrate()` — MOFA2 latent factor analysis across all loaded layers.
* `omicslens_analyze()` — downstream analysis: DESeq2 differential expression, fgsea pathway enrichment, DMRcate DMR calling, Kaplan-Meier survival.
* `omicslens_app()` — interactive Shiny dashboard for exploring all results.
* `omicslens_report()` — parameterised R Markdown HTML report.
* `omicslens_plot_factors()`, `omicslens_plot_de()`, `omicslens_plot_gsea()`, `omicslens_plot_survival()` — individual ggplot2 plot helpers.
* `print.OmicsLens` and `summary.OmicsLens` S3 methods.
