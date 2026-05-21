#' Download and prepare public multi-omics example data
#'
#' Fetches real public data from one of two sources and returns an
#' \code{OmicsLens} object ready for \code{\link{omicslens_preprocess}}:
#'
#' \describe{
#'   \item{\code{"airway"}}{Uses the \href{https://bioconductor.org/packages/airway}{airway}
#'     Bioconductor package (GSE52778, Himes et al. 2014). Provides RNA-Seq
#'     count data for 8 airway smooth muscle cell-line samples (4 dexamethasone-
#'     treated, 4 untreated). A synthetic methylation matrix and binary mutation
#'     matrix with matching sample IDs are generated so all three layers of
#'     OmicsLens can be exercised without additional downloads.}
#'   \item{\code{"tcga_brca"}}{Downloads 10 TCGA-BRCA samples (5 tumour + 5
#'     matched normal) from the NCI Genomic Data Commons using
#'     \href{https://bioconductor.org/packages/TCGAbiolinks}{TCGAbiolinks}.
#'     Includes RNA-Seq counts (STAR), 450k methylation beta values, and somatic
#'     mutations (masked MAF). Internet access and TCGAbiolinks are required.}
#' }
#'
#' @param source Character; one of \code{"airway"} (default, no internet
#'   required) or \code{"tcga_brca"} (internet + TCGAbiolinks required).
#' @param cache_dir Directory to cache downloaded files. Defaults to a
#'   temporary directory when \code{NULL}. Only used when
#'   \code{source = "tcga_brca"}.
#' @param verbose Logical; print progress messages (default \code{TRUE}).
#'
#' @return An \code{OmicsLens} object with real data loaded and ready for
#'   \code{\link{omicslens_preprocess}}.
#'
#' @examples
#' \dontrun{
#' # Quick demo — no internet needed
#' obj <- omicslens_download_example(source = "airway")
#' obj <- omicslens_preprocess(obj)
#'
#' # Full TCGA-BRCA multi-omics (needs internet + TCGAbiolinks)
#' obj <- omicslens_download_example(source = "tcga_brca",
#'                                    cache_dir = "~/omicslens_cache")
#' }
#'
#' @seealso \code{\link{omicslens_load}}, \code{\link{omicslens_preprocess}}
#' @export
omicslens_download_example <- function(source    = c("airway", "tcga_brca"),
                                        cache_dir = NULL,
                                        verbose   = TRUE) {
  source <- match.arg(source)

  if (source == "airway") {
    .download_airway(verbose = verbose)
  } else {
    .download_tcga_brca(cache_dir = cache_dir, verbose = verbose)
  }
}

# ── airway source ────────────────────────────────────────────────────────────
.download_airway <- function(verbose) {
  if (!requireNamespace("airway", quietly = TRUE))
    stop(
      "The 'airway' package is required:\n",
      "  BiocManager::install('airway')",
      call. = FALSE
    )

  .msg(verbose, "Loading airway RNA-Seq data (GSE52778, Himes et al. 2014)...")
  utils::data("airway", package = "airway", envir = environment())
  se <- get("airway", envir = environment())

  counts <- SummarizedExperiment::assay(se, "counts")
  sample_ids <- colnames(counts)

  meta <- as.data.frame(SummarizedExperiment::colData(se))
  meta$sample_id <- rownames(meta)
  meta$condition <- as.character(meta$dex)   # "trt" | "untrt"
  meta$time_os   <- stats::runif(nrow(meta), 12, 60)  # synthetic survival
  meta$event     <- stats::rbinom(nrow(meta), 1, 0.4)

  # Synthetic 450k methylation matrix that mirrors the 8 airway samples
  set.seed(123)
  n_cpgs <- 5000L
  meth_base <- stats::rbeta(n_cpgs, 2, 5)
  meth_beta  <- matrix(
    pmin(pmax(
      outer(meth_base,
            stats::rnorm(length(sample_ids), 1, 0.05), "*"),
      0), 1),
    nrow     = n_cpgs,
    dimnames = list(paste0("cg", seq_len(n_cpgs)), sample_ids)
  )

  # Synthetic binary mutation matrix
  n_mut <- 200L
  var_bin <- matrix(
    stats::rbinom(length(sample_ids) * n_mut, 1, 0.07),
    nrow     = length(sample_ids),
    dimnames = list(sample_ids, paste0("MutGene", seq_len(n_mut)))
  )

  .msg(verbose, sprintf(
    "  Airway: %d genes x %d samples | SRR accessions in inst/extdata/sra_sample_info.tsv",
    nrow(counts), ncol(counts)
  ))

  omicslens_load(
    rna_counts     = counts,
    methylation    = meth_beta,
    variants       = var_bin,
    metadata       = meta,
    variant_format = "binary_matrix",
    verbose        = verbose
  )
}

# ── TCGA-BRCA source ─────────────────────────────────────────────────────────
.download_tcga_brca <- function(cache_dir, verbose) {
  for (pkg in c("TCGAbiolinks", "SummarizedExperiment")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      stop(sprintf("'%s' is required: BiocManager::install('%s')", pkg, pkg),
           call. = FALSE)
  }

  if (is.null(cache_dir))
    cache_dir <- file.path(tempdir(), "OmicsLens_TCGA")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  .msg(verbose, paste("TCGA download cache:", cache_dir))

  barcodes <- c(
    # Primary tumours (01A)
    "TCGA-A2-A0T2", "TCGA-A2-A0CM", "TCGA-A2-A0D0",
    "TCGA-A2-A0EQ", "TCGA-A2-A0EV",
    # Matched normals (11A)
    "TCGA-A2-A0T2", "TCGA-A2-A0CM", "TCGA-A2-A0D0",
    "TCGA-A2-A0EQ", "TCGA-A2-A0EV"
  )

  # ── RNA-Seq (STAR counts) ──────────────────────────────────────────────────
  .msg(verbose, "Querying TCGA-BRCA RNA-Seq counts...")
  q_rna <- TCGAbiolinks::GDCquery(
    project           = "TCGA-BRCA",
    data.category     = "Transcriptome Profiling",
    data.type         = "Gene Expression Quantification",
    workflow.type     = "STAR - Counts",
    barcode           = barcodes
  )
  TCGAbiolinks::GDCdownload(q_rna, directory = cache_dir, method = "api")
  se_rna <- TCGAbiolinks::GDCprepare(q_rna, directory = cache_dir)
  counts <- SummarizedExperiment::assay(se_rna, "unstranded")
  rownames(counts) <- SummarizedExperiment::rowData(se_rna)$gene_name
  .msg(verbose, sprintf("  RNA: %d genes x %d samples", nrow(counts), ncol(counts)))

  # ── Methylation (450k) ────────────────────────────────────────────────────
  .msg(verbose, "Querying TCGA-BRCA 450k methylation...")
  q_meth <- TCGAbiolinks::GDCquery(
    project       = "TCGA-BRCA",
    data.category = "DNA Methylation",
    data.type     = "Methylation Beta Value",
    platform      = "Illumina Human Methylation 450",
    barcode       = barcodes
  )
  TCGAbiolinks::GDCdownload(q_meth, directory = cache_dir, method = "api")
  se_meth <- TCGAbiolinks::GDCprepare(q_meth, directory = cache_dir)
  beta <- SummarizedExperiment::assay(se_meth)
  .msg(verbose, sprintf("  Methylation: %d CpGs x %d samples", nrow(beta), ncol(beta)))

  # ── Somatic mutations (masked MAF) ────────────────────────────────────────
  .msg(verbose, "Querying TCGA-BRCA somatic mutations...")
  q_maf <- TCGAbiolinks::GDCquery(
    project       = "TCGA-BRCA",
    data.category = "Simple Nucleotide Variation",
    data.type     = "Masked Somatic Mutation",
    barcode       = barcodes
  )
  TCGAbiolinks::GDCdownload(q_maf, directory = cache_dir, method = "api")
  maf_df <- TCGAbiolinks::GDCprepare(q_maf, directory = cache_dir)
  # Trim barcode to 12-char patient ID for joining
  maf_df$Tumor_Sample_Barcode <- substr(maf_df$Tumor_Sample_Barcode, 1, 12)
  .msg(verbose, sprintf("  Mutations: %d somatic variants", nrow(maf_df)))

  # ── Clinical metadata ─────────────────────────────────────────────────────
  meta_path <- system.file("extdata", "tcga_brca_metadata.csv",
                             package = "OmicsLens")
  meta <- utils::read.csv(meta_path, stringsAsFactors = FALSE)
  rownames(meta) <- meta$sample_id

  .msg(verbose, "Building OmicsLens object from TCGA-BRCA data...")
  omicslens_load(
    rna_counts     = counts,
    methylation    = beta,
    variants       = maf_df,
    metadata       = meta,
    variant_format = "maf",
    verbose        = verbose
  )
}
