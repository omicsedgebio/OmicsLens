## ── Integration tests using real public data ─────────────────────────────────
##
## Dataset 1 (airway)  : GSE52778  – Himes et al. 2014, PLoS ONE
##   SRR1039508 / SRR1039509 / SRR1039512 / SRR1039513
##   SRR1039516 / SRR1039517 / SRR1039520 / SRR1039521
##   Available via BiocManager::install("airway") — no internet required.
##
## Dataset 2 (TCGA-BRCA): TCGA-A2-A0T2 / A0CM / A0D0 / A0EQ / A0EV
##   RNA-Seq (STAR), 450k methylation, somatic MAF — downloaded via TCGAbiolinks.
##   Tests tagged skip_if_not_installed("TCGAbiolinks") and skip_if_offline().
##
## These tests verify the complete omicslens_load → preprocess pipeline on
## real data. The MOFA2 integration step is skipped (Python dependency) but
## the full analysis chain is exercised as far as preprocessing.
## ─────────────────────────────────────────────────────────────────────────────

# ── Helpers ───────────────────────────────────────────────────────────────────
.has_internet <- function() {
  tryCatch({
    con <- url("https://www.ncbi.nlm.nih.gov", open = "r")
    close(con)
    TRUE
  }, error = function(e) FALSE)
}

# ── Dataset 1: airway (GSE52778) ──────────────────────────────────────────────

test_that("airway: omicslens_download_example returns a valid OmicsLens object", {
  skip_if_not_installed("airway")
  skip_if_not_installed("SummarizedExperiment")

  obj <- omicslens_download_example(source = "airway", verbose = FALSE)

  expect_s3_class(obj, "OmicsLens")
  expect_true("rna" %in% obj$layers_present)
  expect_true("methylation" %in% obj$layers_present)
  expect_true("variants" %in% obj$layers_present)
})

test_that("airway: RNA sample IDs are SRR accessions", {
  skip_if_not_installed("airway")
  skip_if_not_installed("SummarizedExperiment")

  obj <- omicslens_download_example(source = "airway", verbose = FALSE)
  expect_true(all(grepl("^SRR", colnames(obj$rna$counts))))
})

test_that("airway: 8 samples, >10,000 genes", {
  skip_if_not_installed("airway")
  skip_if_not_installed("SummarizedExperiment")

  obj <- omicslens_download_example(source = "airway", verbose = FALSE)
  expect_equal(ncol(obj$rna$counts), 8L)
  expect_gt(nrow(obj$rna$counts), 10000L)
})

test_that("airway: methylation layer has 8000 CpGs", {
  skip_if_not_installed("airway")
  skip_if_not_installed("SummarizedExperiment")

  obj <- omicslens_download_example(source = "airway", verbose = FALSE)
  expect_equal(nrow(obj$methylation$beta), 8000L)
})

test_that("airway: metadata contains expected columns", {
  skip_if_not_installed("airway")
  skip_if_not_installed("SummarizedExperiment")

  obj <- omicslens_download_example(source = "airway", verbose = FALSE)
  expect_true(all(c("sample_id","condition","time_os","event") %in%
                    colnames(obj$metadata)))
})

test_that("airway: preprocess completes on real data", {
  skip_if_not_installed("airway")
  skip_if_not_installed("DESeq2")
  skip_if_not_installed("SummarizedExperiment")

  obj  <- omicslens_download_example(source = "airway", verbose = FALSE)
  obj2 <- omicslens_preprocess(obj, meth_var_top_n = 2000L, verbose = FALSE)

  # VST normalisation must produce a matrix with the same samples
  expect_equal(ncol(obj2$rna$normalized), 8L)
  expect_true(all(is.finite(obj2$rna$normalized)))

  # M-values must be finite and cover exactly 2000 CpGs
  expect_equal(nrow(obj2$methylation$m_values), 2000L)
  expect_true(all(is.finite(obj2$methylation$m_values)))
})

test_that("airway: preprocessed RNA has fewer genes than raw (filtering)", {
  skip_if_not_installed("airway")
  skip_if_not_installed("DESeq2")
  skip_if_not_installed("SummarizedExperiment")

  obj  <- omicslens_download_example(source = "airway", verbose = FALSE)
  obj2 <- omicslens_preprocess(obj, rna_min_counts = 10L, verbose = FALSE)
  expect_lt(nrow(obj2$rna$counts), nrow(obj$rna$counts))
})

test_that("airway: SRA sample info file is readable", {
  sra_path <- system.file("extdata", "sra_sample_info.tsv", package = "OmicsLens")
  expect_true(file.exists(sra_path))
  sra <- read.table(sra_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  expect_equal(nrow(sra), 8L)
  expect_true("SRR_accession" %in% colnames(sra))
  expect_true(all(grepl("^SRR", sra$SRR_accession)))
})

test_that("airway: variant MAF filter leaves only frequently-mutated genes", {
  skip_if_not_installed("airway")
  skip_if_not_installed("SummarizedExperiment")

  obj  <- omicslens_download_example(source = "airway", verbose = FALSE)
  # With MAF threshold above typical noise, only common events survive
  obj2 <- omicslens_preprocess(obj, variant_min_maf = 0.5, verbose = FALSE)
  expect_true(ncol(obj2$variants$binary) < ncol(obj$variants$binary))
})

# ── Dataset 2: TCGA-BRCA (network tests) ──────────────────────────────────────

test_that("TCGA metadata CSV is bundled and parses correctly", {
  meta_path <- system.file("extdata", "tcga_brca_metadata.csv",
                             package = "OmicsLens")
  expect_true(file.exists(meta_path))
  meta <- read.csv(meta_path, stringsAsFactors = FALSE)
  expect_true(all(c("sample_id","barcode","tissue","PAM50","time_os","event")
                    %in% colnames(meta)))
  expect_gte(nrow(meta), 10L)
  tumour_rows <- meta[meta$tissue == "Tumour", ]
  expect_gte(nrow(tumour_rows), 5L)
})

test_that("TCGA: omicslens_download_example(tcga_brca) requires internet", {
  skip_if_not_installed("TCGAbiolinks")
  skip_if(!.has_internet(), "No internet connection")
  skip_on_cran()

  obj <- omicslens_download_example(
    source    = "tcga_brca",
    cache_dir = file.path(tempdir(), "tcga_test"),
    verbose   = FALSE
  )
  expect_s3_class(obj, "OmicsLens")
  expect_true("rna" %in% obj$layers_present)
})

test_that("TCGA: RNA layer has Ensembl/symbol rownames after download", {
  skip_if_not_installed("TCGAbiolinks")
  skip_if(!.has_internet(), "No internet connection")
  skip_on_cran()

  obj <- omicslens_download_example(
    source    = "tcga_brca",
    cache_dir = file.path(tempdir(), "tcga_test"),
    verbose   = FALSE
  )
  rn <- rownames(obj$rna$counts)
  # Gene symbols should not look like TCGA barcodes
  expect_false(all(grepl("^TCGA", rn)))
})

test_that("TCGA: methylation beta values are in [0,1]", {
  skip_if_not_installed("TCGAbiolinks")
  skip_if(!.has_internet(), "No internet connection")
  skip_on_cran()

  obj <- omicslens_download_example(
    source    = "tcga_brca",
    cache_dir = file.path(tempdir(), "tcga_test"),
    verbose   = FALSE
  )
  skip_if(is.null(obj$methylation))
  vals <- obj$methylation$beta
  expect_true(all(vals >= 0 & vals <= 1, na.rm = TRUE))
})

test_that("TCGA: MAF contains Hugo_Symbol and Tumor_Sample_Barcode columns", {
  skip_if_not_installed("TCGAbiolinks")
  skip_if(!.has_internet(), "No internet connection")
  skip_on_cran()

  obj <- omicslens_download_example(
    source    = "tcga_brca",
    cache_dir = file.path(tempdir(), "tcga_test"),
    verbose   = FALSE
  )
  skip_if(is.null(obj$variants))
  expect_true(all(c("Hugo_Symbol","Tumor_Sample_Barcode") %in%
                    colnames(obj$variants$raw)))
})

test_that("TCGA: preprocess runs end-to-end on real TCGA data", {
  skip_if_not_installed("TCGAbiolinks")
  skip_if_not_installed("DESeq2")
  skip_if(!.has_internet(), "No internet connection")
  skip_on_cran()

  obj  <- omicslens_download_example(
    source    = "tcga_brca",
    cache_dir = file.path(tempdir(), "tcga_test"),
    verbose   = FALSE
  )
  obj2 <- omicslens_preprocess(obj, meth_var_top_n = 5000L, verbose = FALSE)

  expect_false(is.null(obj2$rna$normalized))
  if (!is.null(obj2$methylation))
    expect_false(is.null(obj2$methylation$m_values))
})
