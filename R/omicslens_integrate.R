#' Integrate multi-omics layers using MOFA2
#'
#' Runs Multi-Omics Factor Analysis v2 (MOFA2) on the preprocessed data layers
#' stored in an \code{OmicsLens} object. Latent factors that explain shared and
#' layer-specific variance are extracted and stored for downstream analysis.
#'
#' @param obj An \code{OmicsLens} object that has been processed with
#'   \code{\link{omicslens_preprocess}}.
#' @param n_factors Number of latent factors to learn (default \code{10}).
#'   Start with 10–15 and inspect the variance-explained plot.
#' @param scale_views Logical; scale each view to unit variance before fitting
#'   (default \code{TRUE}, recommended when combining different assay types).
#' @param seed Integer random seed for reproducibility (default \code{42}).
#' @param maxiter Maximum number of training iterations (default \code{1000}).
#' @param use_basilisk Logical; use the basilisk Python environment bundled
#'   with MOFA2 (default \code{TRUE}). Set to \code{FALSE} only if you have
#'   manually configured a Python environment with mofapy2 installed.
#' @param python_path Character; path to the Python binary that has
#'   \code{mofapy2} installed. Only used when \code{use_basilisk = FALSE}.
#'   If \code{NULL} (default), \code{reticulate} will use its default Python.
#' @param verbose Logical; print progress messages (default \code{TRUE}).
#'
#' @return The input \code{OmicsLens} object with two new entries:
#' \describe{
#'   \item{\code{mofa_model}}{The trained \code{MOFA} S4 object.}
#'   \item{\code{results$mofa}}{A list containing:
#'     \describe{
#'       \item{\code{factors}}{Numeric matrix of sample factor scores
#'         (samples × factors).}
#'       \item{\code{weights}}{Named list of feature weight matrices per view.}
#'       \item{\code{variance_explained}}{List of variance-explained values
#'         from \code{MOFA2::get_variance_explained()}.}
#'     }
#'   }
#' }
#'
#' @details
#' MOFA2 requires Python and the \code{mofapy2} package, which are managed
#' automatically via \code{basilisk} when \code{use_basilisk = TRUE}. The
#' first run may take several minutes as the Python environment is set up.
#'
#' Variant data is modelled with a Bernoulli likelihood (binary); RNA and
#' methylation layers use a Gaussian likelihood.
#'
#' @references
#' Argelaguet R et al. (2020). MOFA+: a statistical framework for
#' comprehensive integration of multi-modal single-cell data.
#' \emph{Genome Biology}, 21:111.
#'
#' @examples
#' \dontrun{
#' obj <- omicslens_load(rna_counts = rna_mat, methylation = meth_mat,
#'                        metadata = meta_df)
#' obj <- omicslens_preprocess(obj)
#' obj <- omicslens_integrate(obj, n_factors = 10)
#' head(obj$results$mofa$factors)
#' }
#'
#' @seealso \code{\link{omicslens_preprocess}}, \code{\link{omicslens_analyze}}
#' @export
omicslens_integrate <- function(obj,
                                 n_factors    = 10L,
                                 scale_views  = TRUE,
                                 seed         = 42L,
                                 maxiter      = 1000L,
                                 use_basilisk = TRUE,
                                 python_path  = NULL,
                                 verbose      = TRUE) {

  if (!inherits(obj, "OmicsLens"))
    stop("obj must be an OmicsLens object.", call. = FALSE)
  if (!requireNamespace("MOFA2", quietly = TRUE))
    stop("MOFA2 is required: BiocManager::install('MOFA2')", call. = FALSE)

  # ── Build MOFA data list ────────────────────────────────────────────────────
  mofa_data <- list()

  if ("rna" %in% obj$layers_present) {
    if (is.null(obj$rna$normalized))
      stop("RNA data is not normalised. Run omicslens_preprocess() first.",
           call. = FALSE)
    mofa_data[["RNA"]] <- obj$rna$normalized
  }

  if ("variants" %in% obj$layers_present) {
    if (is.null(obj$variants$binary))
      stop("Variant binary matrix is missing. Run omicslens_preprocess() first.",
           call. = FALSE)
    # MOFA2 expects features × samples
    mofa_data[["Variants"]] <- t(obj$variants$binary)
  }

  if ("methylation" %in% obj$layers_present) {
    if (is.null(obj$methylation$m_values))
      stop("Methylation M-values are missing. Run omicslens_preprocess() first.",
           call. = FALSE)
    mofa_data[["Methylation"]] <- obj$methylation$m_values
  }

  if (length(mofa_data) == 0)
    stop("No preprocessed layers found. Run omicslens_preprocess() first.",
         call. = FALSE)

  .msg(verbose, sprintf("Creating MOFA2 object with %d view(s): %s",
                         length(mofa_data),
                         paste(names(mofa_data), collapse = ", ")))

  # ── Create and configure MOFA object ────────────────────────────────────────
  mofa_obj <- MOFA2::create_mofa(mofa_data)

  data_opts <- MOFA2::get_default_data_options(mofa_obj)
  data_opts$scale_views <- scale_views

  model_opts <- MOFA2::get_default_model_options(mofa_obj)
  model_opts$num_factors <- as.integer(n_factors)

  # Bernoulli likelihood for binary variant layer
  if ("Variants" %in% names(mofa_data))
    model_opts$likelihoods["Variants"] <- "bernoulli"

  train_opts <- MOFA2::get_default_training_options(mofa_obj)
  train_opts$seed    <- as.integer(seed)
  train_opts$maxiter <- as.integer(maxiter)
  train_opts$verbose <- verbose

  mofa_obj <- MOFA2::prepare_mofa(
    mofa_obj,
    data_options     = data_opts,
    model_options    = model_opts,
    training_options = train_opts
  )

  if (!use_basilisk && !is.null(python_path)) {
    if (!requireNamespace("reticulate", quietly = TRUE))
      stop("reticulate is required when python_path is specified.",
           call. = FALSE)
    reticulate::use_python(python_path, required = TRUE)
    .msg(verbose, sprintf("Using Python: %s", python_path))
  }

  .msg(verbose, "Training MOFA2 model (this may take several minutes)...")
  mofa_obj <- MOFA2::run_mofa(mofa_obj, use_basilisk = use_basilisk)

  # ── Extract results ──────────────────────────────────────────────────────────
  factor_list <- MOFA2::get_factors(mofa_obj)
  # MOFA2 returns a list with one element per group; flatten for single-group data
  factors_mat <- do.call(rbind, factor_list)
  colnames(factors_mat) <- paste0("Factor", seq_len(ncol(factors_mat)))

  obj$mofa_model        <- mofa_obj
  obj$results[["mofa"]] <- list(
    factors            = factors_mat,
    weights            = MOFA2::get_weights(mofa_obj),
    variance_explained = MOFA2::get_variance_explained(mofa_obj)
  )

  .msg(verbose, sprintf(
    "MOFA2 complete. %d factors learned across %d views.",
    ncol(factors_mat), length(mofa_data)
  ))

  obj
}
