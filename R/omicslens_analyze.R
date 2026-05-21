#' Run downstream multi-omics analyses
#'
#' Uses MOFA2 factor scores to define sample groups and then runs:
#' differential expression (DESeq2), gene set enrichment (fgsea),
#' differentially methylated region calling (DMRcate), and Kaplan-Meier
#' survival analysis.
#'
#' @param obj An \code{OmicsLens} object with a trained MOFA2 model
#'   (from \code{\link{omicslens_integrate}}).
#' @param run_de Logical; run DESeq2 differential expression (default \code{TRUE}).
#' @param de_factor Integer index of the MOFA factor used to define High/Low
#'   sample groups for differential expression (default \code{1}).
#' @param run_gsea Logical; run fgsea pathway enrichment on DE results
#'   (default \code{TRUE}).
#' @param organism Character; \code{"human"} (default) or \code{"mouse"}.
#'   Controls which MSigDB gene sets are loaded (requires \pkg{msigdbr};
#'   falls back to built-in Hallmark stubs if unavailable).
#' @param run_dmr Logical; call differentially methylated regions with DMRcate
#'   (default \code{TRUE}). Silently skipped if the methylation layer is absent.
#' @param run_survival Logical; run Kaplan-Meier + Cox regression
#'   (default \code{TRUE}).
#' @param survival_time_col Name of the metadata column containing follow-up
#'   time (default \code{"time_os"}).
#' @param survival_event_col Name of the metadata column containing the event
#'   indicator (1 = event, 0 = censored; default \code{"event"}).
#' @param verbose Logical; print progress messages (default \code{TRUE}).
#'
#' @return The input \code{OmicsLens} object with analysis results appended to
#' \code{obj$results}:
#' \describe{
#'   \item{\code{results$groups}}{Named character vector mapping each sample
#'     to \code{"High"} or \code{"Low"} based on the chosen MOFA factor.}
#'   \item{\code{results$de}}{data.frame of DESeq2 results (if \code{run_de}).}
#'   \item{\code{results$gsea}}{data.table of fgsea results (if \code{run_gsea}).}
#'   \item{\code{results$dmr}}{data.frame of DMR ranges (if \code{run_dmr}).}
#'   \item{\code{results$survival}}{List with \code{fit} (survfit object),
#'     \code{cox} (coxph object), and \code{data} (data.frame used for fitting).}
#' }
#'
#' @examples
#' \dontrun{
#' obj <- omicslens_load(rna_counts = rna_mat, metadata = meta_df)
#' obj <- omicslens_preprocess(obj)
#' obj <- omicslens_integrate(obj)
#' obj <- omicslens_analyze(obj, survival_time_col = "os_months",
#'                           survival_event_col = "dead")
#' head(obj$results$de)
#' }
#'
#' @seealso \code{\link{omicslens_integrate}}, \code{\link{omicslens_app}},
#'   \code{\link{omicslens_report}}
#' @export
omicslens_analyze <- function(obj,
                               run_de              = TRUE,
                               de_factor           = 1L,
                               run_gsea            = TRUE,
                               organism            = "human",
                               run_dmr             = TRUE,
                               run_survival        = TRUE,
                               survival_time_col   = "time_os",
                               survival_event_col  = "event",
                               verbose             = TRUE) {

  if (!inherits(obj, "OmicsLens"))
    stop("obj must be an OmicsLens object.", call. = FALSE)
  if (is.null(obj$mofa_model))
    stop("No MOFA2 model found. Run omicslens_integrate() first.", call. = FALSE)

  factors <- obj$results$mofa$factors
  if (de_factor > ncol(factors))
    stop(sprintf("de_factor %d exceeds number of factors (%d).",
                 de_factor, ncol(factors)), call. = FALSE)

  # ── Define sample groups from Factor de_factor (median split) ───────────────
  f_vals <- factors[, de_factor]
  groups <- ifelse(f_vals >= stats::median(f_vals), "High", "Low")
  names(groups)           <- rownames(factors)
  obj$results[["groups"]] <- groups

  .msg(verbose, sprintf(
    "Factor %d split: %d High, %d Low",
    de_factor, sum(groups == "High"), sum(groups == "Low")
  ))

  # ── Differential Expression ─────────────────────────────────────────────────
  if (run_de && "rna" %in% obj$layers_present) {
    .msg(verbose, "Running DESeq2 differential expression...")
    obj$results[["de"]] <- tryCatch(
      .run_deseq2(obj, groups),
      error = function(e) {
        warning("DESeq2 failed: ", conditionMessage(e), call. = FALSE)
        NULL
      }
    )
    if (!is.null(obj$results$de)) {
      n_sig <- sum(obj$results$de$padj < 0.05, na.rm = TRUE)
      .msg(verbose, sprintf("  DE: %d significant genes (padj < 0.05)", n_sig))
    }
  }

  # ── Gene Set Enrichment ─────────────────────────────────────────────────────
  if (run_gsea && !is.null(obj$results$de)) {
    .msg(verbose, "Running fgsea pathway enrichment...")
    obj$results[["gsea"]] <- tryCatch(
      .run_gsea(obj$results$de, organism = organism),
      error = function(e) {
        warning("fgsea failed: ", conditionMessage(e), call. = FALSE)
        NULL
      }
    )
    if (!is.null(obj$results$gsea)) {
      n_sig <- sum(obj$results$gsea$padj < 0.05, na.rm = TRUE)
      .msg(verbose, sprintf("  GSEA: %d significant pathways (padj < 0.05)", n_sig))
    }
  }

  # ── DMR Calling ─────────────────────────────────────────────────────────────
  if (run_dmr && "methylation" %in% obj$layers_present) {
    .msg(verbose, "Calling differentially methylated regions (DMRcate)...")
    obj$results[["dmr"]] <- .run_dmr(obj, groups)
    if (!is.null(obj$results$dmr))
      .msg(verbose, sprintf("  DMR: %d regions identified", nrow(obj$results$dmr)))
  }

  # ── Survival Analysis ────────────────────────────────────────────────────────
  if (run_survival && !is.null(obj$metadata)) {
    if (all(c(survival_time_col, survival_event_col) %in% colnames(obj$metadata))) {
      .msg(verbose, "Running Kaplan-Meier and Cox regression...")
      obj$results[["survival"]] <- tryCatch(
        .run_survival(obj$metadata, groups,
                       survival_time_col, survival_event_col),
        error = function(e) {
          warning("Survival analysis failed: ", conditionMessage(e), call. = FALSE)
          NULL
        }
      )
      if (!is.null(obj$results$survival))
        .msg(verbose, "  Survival: analysis complete")
    } else {
      .msg(verbose, sprintf(
        "Skipping survival: columns '%s' or '%s' not found in metadata.",
        survival_time_col, survival_event_col
      ))
    }
  }

  .msg(verbose, "Analysis complete.")
  obj
}
