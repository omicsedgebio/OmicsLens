#' Preprocess each data layer in an OmicsLens object
#'
#' Applies layer-specific normalisation and filtering:
#' \itemize{
#'   \item \strong{RNA-Seq}: removes low-count genes, then applies DESeq2
#'     variance-stabilising transformation (VST).
#'   \item \strong{Variants}: filters genes below the minimum mutation
#'     allele frequency (MAF) across samples.
#'   \item \strong{Methylation}: removes CpGs with excessive missing values,
#'     converts beta values to M-values (logit), and retains the top most
#'     variable CpG sites.
#' }
#'
#' @param obj An \code{OmicsLens} object produced by \code{\link{omicslens_load}}.
#' @param rna_min_counts Minimum count threshold per gene per sample. Genes
#'   with fewer than \code{rna_min_counts} counts in at least
#'   \code{rna_min_samples} samples are removed (default \code{10}).
#' @param rna_min_samples Minimum number of samples that must exceed
#'   \code{rna_min_counts}. Defaults to 10\% of samples (minimum 2).
#' @param variant_min_maf Minimum mutation frequency (proportion of samples)
#'   for a gene to be retained in the variant layer. Default \code{0.05} (5\%).
#' @param meth_max_na_frac Maximum fraction of missing values allowed per CpG
#'   site before the site is discarded (default \code{0.2}).
#' @param meth_var_top_n Number of most-variable CpG sites to retain for
#'   downstream analysis (default \code{10000}).
#' @param verbose Logical; print progress messages (default \code{TRUE}).
#'
#' @return The input \code{OmicsLens} object with preprocessed data filled in:
#' \describe{
#'   \item{\code{rna$normalized}}{VST-normalised gene expression matrix.}
#'   \item{\code{variants$binary}}{MAF-filtered binary mutation matrix.}
#'   \item{\code{methylation$m_values}}{M-value matrix (top variable CpGs).}
#' }
#'
#' @examples
#' \dontrun{
#' obj <- omicslens_load(rna_counts = rna_mat, metadata = meta_df)
#' obj <- omicslens_preprocess(obj)
#' }
#'
#' @seealso \code{\link{omicslens_load}}, \code{\link{omicslens_integrate}}
#' @export
omicslens_preprocess <- function(obj,
                                  rna_min_counts   = 10L,
                                  rna_min_samples  = NULL,
                                  variant_min_maf  = 0.05,
                                  meth_max_na_frac = 0.2,
                                  meth_var_top_n   = 10000L,
                                  verbose          = TRUE) {

  if (!inherits(obj, "OmicsLens"))
    stop("obj must be an OmicsLens object from omicslens_load().", call. = FALSE)

  # ── RNA ──────────────────────────────────────────────────────────────────────
  if ("rna" %in% obj$layers_present) {
    .msg(verbose, "Preprocessing RNA-Seq...")
    if (!requireNamespace("DESeq2", quietly = TRUE))
      stop("DESeq2 is required: BiocManager::install('DESeq2')", call. = FALSE)
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE))
      stop("SummarizedExperiment required: BiocManager::install('SummarizedExperiment')",
           call. = FALSE)

    counts    <- obj$rna$counts
    n_samples <- ncol(counts)
    min_samp  <- rna_min_samples %||% max(2L, floor(n_samples * 0.10))

    keep_rna <- rowSums(counts >= rna_min_counts) >= min_samp
    counts   <- counts[keep_rna, , drop = FALSE]
    .msg(verbose, sprintf("  RNA: kept %d / %d genes (counts >= %d in >= %d samples)",
                           sum(keep_rna), length(keep_rna),
                           rna_min_counts, min_samp))

    col_data <- if (!is.null(obj$metadata)) {
      obj$metadata[colnames(counts), , drop = FALSE]
    } else {
      data.frame(row.names = colnames(counts),
                 sample_id = colnames(counts))
    }
    col_data$sample_id <- NULL  # avoid duplicate col warning

    dds <- DESeq2::DESeqDataSetFromMatrix(
      countData = round(counts),
      colData   = col_data,
      design    = ~1            # intercept-only; purely for normalisation
    )
    dds <- DESeq2::estimateSizeFactors(dds)
    vst <- tryCatch(
      DESeq2::vst(dds, blind = TRUE),
      error = function(e) DESeq2::varianceStabilizingTransformation(dds, blind = TRUE)
    )
    vst_mat <- SummarizedExperiment::assay(vst)

    obj$rna$counts     <- counts
    obj$rna$normalized <- vst_mat
    .msg(verbose, "  RNA: DESeq2 VST normalisation complete")
  }

  # ── Variants ─────────────────────────────────────────────────────────────────
  if ("variants" %in% obj$layers_present) {
    .msg(verbose, "Preprocessing variant data...")
    bin    <- obj$variants$binary
    maf    <- colMeans(bin, na.rm = TRUE)
    keep_v <- maf >= variant_min_maf
    obj$variants$binary <- bin[, keep_v, drop = FALSE]
    .msg(verbose, sprintf("  Variants: kept %d / %d genes (MAF >= %.2f)",
                           sum(keep_v), length(keep_v), variant_min_maf))
  }

  # ── Methylation ───────────────────────────────────────────────────────────────
  if ("methylation" %in% obj$layers_present) {
    .msg(verbose, "Preprocessing methylation data...")
    beta <- obj$methylation$beta

    # Remove high-NA CpGs
    na_frac <- rowMeans(is.na(beta))
    beta    <- beta[na_frac < meth_max_na_frac, , drop = FALSE]
    .msg(verbose, sprintf("  Methylation: %d CpGs after NA filtering", nrow(beta)))

    # Beta → M-values (clip to avoid ±Inf)
    beta_c <- pmin(pmax(beta, 1e-3), 1 - 1e-3)
    m_vals <- log2(beta_c / (1 - beta_c))

    # Select top variable CpGs
    cpg_var <- apply(m_vals, 1, stats::var, na.rm = TRUE)
    top_n   <- min(as.integer(meth_var_top_n), nrow(m_vals))
    top_idx <- order(cpg_var, decreasing = TRUE)[seq_len(top_n)]
    m_vals  <- m_vals[top_idx, , drop = FALSE]

    obj$methylation$beta     <- beta[top_idx, , drop = FALSE]
    obj$methylation$m_values <- m_vals
    .msg(verbose, sprintf("  Methylation: top %d variable CpGs retained; M-values computed",
                           top_n))
  }

  .msg(verbose, "Preprocessing complete.")
  obj
}
