#' Load and validate multi-omics data
#'
#' Reads RNA-Seq count matrices, variant data (MAF or binary mutation matrix),
#' and DNA methylation beta matrices from files or in-memory objects. Validates
#' inputs and aligns sample IDs across all provided layers.
#'
#' @param rna_counts Path to a CSV/TSV/RDS file or a numeric matrix/data.frame
#'   of raw RNA-Seq counts (genes as rows, samples as columns). Row names must
#'   be gene identifiers. Integer count values are expected.
#' @param variants Path to a MAF file, CSV with binary values, or a numeric
#'   matrix (samples as rows, genes as columns, 0/1 values). Pass \code{NULL}
#'   to exclude the variant layer.
#' @param methylation Path to a CSV/TSV/RDS file or a numeric matrix/data.frame
#'   of DNA methylation beta values (CpG sites as rows, samples as columns).
#'   Values must be in \eqn{[0, 1]}. Pass \code{NULL} to exclude.
#' @param metadata A \code{data.frame} or path to a CSV file containing
#'   sample-level metadata. Must include a \code{sample_id} column. Optional
#'   columns \code{time_os} and \code{event} enable survival analysis downstream.
#' @param variant_format One of \code{"auto"} (default), \code{"maf"}, or
#'   \code{"binary_matrix"}. With \code{"auto"}, the format is inferred from
#'   the file extension and content.
#' @param verbose Logical; print loading progress messages (default \code{TRUE}).
#'
#' @return An S3 object of class \code{"OmicsLens"} with slots:
#' \describe{
#'   \item{\code{rna}}{List with element \code{counts} (matrix, genes × samples)
#'     and \code{normalized} (filled after \code{\link{omicslens_preprocess}}).}
#'   \item{\code{variants}}{List with \code{raw} (data.frame or matrix) and
#'     \code{binary} (binary matrix, samples × genes).}
#'   \item{\code{methylation}}{List with \code{beta} (matrix, CpGs × samples)
#'     and \code{m_values} (filled after preprocessing).}
#'   \item{\code{metadata}}{Sample-level metadata data.frame, or \code{NULL}.}
#'   \item{\code{layers_present}}{Character vector of loaded layer names.}
#'   \item{\code{mofa_model}}{Trained MOFA2 model, \code{NULL} until
#'     \code{\link{omicslens_integrate}} is called.}
#'   \item{\code{results}}{Named list of analysis results, empty until
#'     \code{\link{omicslens_analyze}} is called.}
#' }
#'
#' @examples
#' \dontrun{
#' # Load from files
#' obj <- omicslens_load(
#'   rna_counts  = "counts.csv",
#'   variants    = "mutations.maf",
#'   methylation = "beta_matrix.csv",
#'   metadata    = "sample_info.csv"
#' )
#' print(obj)
#'
#' # Load from in-memory matrices
#' set.seed(1)
#' rna  <- matrix(rpois(5000, 50), nrow = 100, ncol = 50,
#'                dimnames = list(paste0("Gene", 1:100),
#'                                paste0("S", 1:50)))
#' meta <- data.frame(sample_id = paste0("S", 1:50),
#'                    condition  = rep(c("Tumour","Normal"), 25))
#' obj  <- omicslens_load(rna_counts = rna, metadata = meta)
#' print(obj)
#' }
#'
#' @seealso \code{\link{omicslens_preprocess}}, \code{\link{omicslens_integrate}}
#' @export
omicslens_load <- function(rna_counts     = NULL,
                            variants       = NULL,
                            methylation    = NULL,
                            metadata       = NULL,
                            variant_format = "auto",
                            verbose        = TRUE) {

  # ── 1. RNA ─────────────────────────────────────────────────────────────────
  rna_layer <- NULL
  if (!is.null(rna_counts)) {
    .msg(verbose, "Loading RNA-Seq counts...")
    mat <- .read_matrix(rna_counts)
    mode(mat) <- "numeric"
    if (any(mat < 0, na.rm = TRUE))
      stop("rna_counts contains negative values.", call. = FALSE)
    rna_layer <- list(counts = mat, normalized = NULL)
    .msg(verbose, sprintf("  RNA: %d genes x %d samples", nrow(mat), ncol(mat)))
  }

  # ── 2. Variants ─────────────────────────────────────────────────────────────
  var_layer <- NULL
  if (!is.null(variants)) {
    .msg(verbose, "Loading variant data...")
    fmt <- if (variant_format == "auto") .detect_variant_format(variants) else variant_format
    if (fmt == "maf") {
      var_raw <- .read_maf(variants)
      var_bin <- .maf_to_binary(var_raw)
    } else {
      var_raw <- .read_matrix(variants)
      var_bin <- var_raw
      mode(var_bin) <- "numeric"
    }
    var_layer <- list(raw = var_raw, binary = var_bin)
    .msg(verbose, sprintf("  Variants: %d samples x %d genes",
                           nrow(var_bin), ncol(var_bin)))
  }

  # ── 3. Methylation ─────────────────────────────────────────────────────────
  meth_layer <- NULL
  if (!is.null(methylation)) {
    .msg(verbose, "Loading methylation beta matrix...")
    mat <- .read_matrix(methylation)
    mode(mat) <- "numeric"
    if (any(mat < 0 | mat > 1, na.rm = TRUE))
      stop(
        "Methylation values outside [0,1] detected. ",
        "Provide beta values (not M-values); OmicsLens converts internally.",
        call. = FALSE
      )
    meth_layer <- list(beta = mat, m_values = NULL)
    .msg(verbose, sprintf("  Methylation: %d CpGs x %d samples",
                           nrow(mat), ncol(mat)))
  }

  # ── 4. Metadata ─────────────────────────────────────────────────────────────
  meta <- NULL
  if (!is.null(metadata)) {
    .msg(verbose, "Loading sample metadata...")
    meta <- .read_table(metadata)
    if (!"sample_id" %in% colnames(meta))
      stop("metadata must contain a 'sample_id' column.", call. = FALSE)
    rownames(meta) <- meta$sample_id
  }

  # ── 5. Align sample IDs ─────────────────────────────────────────────────────
  layers_present <- character(0)
  sample_ids     <- list()

  if (!is.null(rna_layer)) {
    layers_present      <- c(layers_present, "rna")
    sample_ids[["rna"]] <- colnames(rna_layer$counts)
  }
  if (!is.null(var_layer)) {
    layers_present           <- c(layers_present, "variants")
    sample_ids[["variants"]] <- rownames(var_layer$binary)
  }
  if (!is.null(meth_layer)) {
    layers_present               <- c(layers_present, "methylation")
    sample_ids[["methylation"]]  <- colnames(meth_layer$beta)
  }

  if (length(layers_present) == 0)
    stop("Provide at least one data layer: rna_counts, variants, or methylation.",
         call. = FALSE)

  if (length(layers_present) > 1) {
    common  <- Reduce(intersect, sample_ids)
    dropped <- vapply(sample_ids, function(x) sum(!x %in% common), integer(1))
    if (any(dropped > 0))
      warning(
        "Sample ID mismatch detected. Retaining ", length(common),
        " common samples.\n",
        "  Dropped per layer: ",
        paste(names(dropped[dropped > 0]),
              dropped[dropped > 0], sep = "=", collapse = ", "),
        "\n  Ensure sample IDs are identical across input files.",
        call. = FALSE
      )
    if (!is.null(rna_layer))
      rna_layer$counts     <- rna_layer$counts[, common, drop = FALSE]
    if (!is.null(var_layer))
      var_layer$binary     <- var_layer$binary[common, , drop = FALSE]
    if (!is.null(meth_layer))
      meth_layer$beta      <- meth_layer$beta[, common, drop = FALSE]
    if (!is.null(meta))
      meta <- meta[meta$sample_id %in% common, , drop = FALSE]
  }

  .msg(verbose, sprintf(
    "Loaded %d layer(s): %s | %d samples",
    length(layers_present),
    paste(layers_present, collapse = ", "),
    length(unique(unlist(sample_ids)))
  ))

  structure(
    list(
      rna            = rna_layer,
      variants       = var_layer,
      methylation    = meth_layer,
      metadata       = meta,
      layers_present = layers_present,
      mofa_model     = NULL,
      results        = list()
    ),
    class = "OmicsLens"
  )
}
