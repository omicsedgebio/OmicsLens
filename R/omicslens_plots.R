#' Plot MOFA2 factor scores
#'
#' Scatter plot of two MOFA2 latent factors, with optional group colouring.
#'
#' @param obj An \code{OmicsLens} object after \code{\link{omicslens_integrate}}.
#' @param factor_x Integer index of the factor for the x-axis (default \code{1}).
#' @param factor_y Integer index of the factor for the y-axis (default \code{2}).
#' @param color_by Name of a metadata column or \code{"group"} to colour points
#'   by the MOFA-derived High/Low group (default \code{"group"}).
#' @return A \code{ggplot2} object.
#' @export
omicslens_plot_factors <- function(obj, factor_x = 1L, factor_y = 2L,
                                    color_by = "group") {
  if (is.null(obj$results$mofa))
    stop("Run omicslens_integrate() first.", call. = FALSE)
  mat <- as.data.frame(obj$results$mofa$factors)
  mat$sample <- rownames(mat)
  fx <- paste0("Factor", factor_x)
  fy <- paste0("Factor", factor_y)
  if (color_by == "group" && !is.null(obj$results$groups))
    mat$Color <- obj$results$groups[mat$sample]
  else if (color_by %in% colnames(obj$metadata %||% data.frame()))
    mat$Color <- obj$metadata[mat$sample, color_by]
  else
    mat$Color <- "All"

  ggplot2::ggplot(mat, ggplot2::aes(x = .data[[fx]], y = .data[[fy]],
                                     color = Color, label = sample)) +
    ggplot2::geom_point(size = 3, alpha = 0.8) +
    ggplot2::scale_color_manual(values = c("High"="#e74c3c","Low"="#3498db",
                                            "All"="#7f8c8d")) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::labs(x = fx, y = fy, color = color_by,
                   title = paste("MOFA2:", fx, "vs", fy))
}

#' Plot a volcano plot of differential expression results
#'
#' @param obj An \code{OmicsLens} object after \code{\link{omicslens_analyze}}.
#' @param padj_cutoff Adjusted p-value significance threshold (default \code{0.05}).
#' @param lfc_cutoff Absolute log2 fold-change threshold (default \code{1}).
#' @param top_n_labels Number of top significant genes to label (default \code{15}).
#' @return A \code{ggplot2} object.
#' @export
omicslens_plot_de <- function(obj, padj_cutoff = 0.05, lfc_cutoff = 1,
                               top_n_labels = 15L) {
  if (is.null(obj$results$de))
    stop("Run omicslens_analyze() first.", call. = FALSE)
  de <- obj$results$de
  de$log10p <- -log10(de$padj + 1e-300)
  de$Sig    <- "NS"
  de$Sig[de$padj < padj_cutoff & de$log2FoldChange >  lfc_cutoff] <- "Up"
  de$Sig[de$padj < padj_cutoff & de$log2FoldChange < -lfc_cutoff] <- "Down"

  top_labels <- utils::head(
    de[de$Sig != "NS", ][order(de[de$Sig != "NS", "padj"]), "gene"],
    top_n_labels
  )
  de$label <- ifelse(de$gene %in% top_labels, de$gene, "")

  ggplot2::ggplot(de, ggplot2::aes(x = log2FoldChange, y = log10p,
                                    color = Sig, label = label)) +
    ggplot2::geom_point(alpha = 0.6, size = 1.5) +
    ggplot2::scale_color_manual(values = c(Up="#e74c3c", Down="#3498db", NS="#bdc3c7")) +
    ggplot2::geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff),
                         linetype = "dashed", color = "grey40") +
    ggplot2::geom_hline(yintercept = -log10(padj_cutoff),
                         linetype = "dashed", color = "grey40") +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::labs(x = "log2 Fold Change (High vs Low)",
                   y = "-log10(padj)",
                   title = "Differential Expression Volcano Plot",
                   color = "")
}

#' Plot GSEA pathway enrichment results
#'
#' @param obj An \code{OmicsLens} object after \code{\link{omicslens_analyze}}.
#' @param top_n Number of top pathways to display (default \code{20}).
#' @return A \code{ggplot2} object.
#' @export
omicslens_plot_gsea <- function(obj, top_n = 20L) {
  if (is.null(obj$results$gsea))
    stop("Run omicslens_analyze() first.", call. = FALSE)
  df <- as.data.frame(obj$results$gsea)
  df <- df[order(df$padj), ]
  df <- utils::head(df, top_n)
  df$pathway <- gsub("HALLMARK_", "", df$pathway)
  df$pathway <- factor(df$pathway, levels = rev(df$pathway))
  ggplot2::ggplot(df, ggplot2::aes(x = NES, y = pathway,
                                    size = size, color = padj)) +
    ggplot2::geom_point() +
    ggplot2::scale_color_gradient(low = "#e74c3c", high = "#95a5a6",
                                   name = "padj") +
    ggplot2::scale_size_continuous(name = "Gene set size", range = c(3, 10)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::labs(x = "Normalised Enrichment Score", y = NULL,
                   title = paste("Top", top_n, "Enriched Pathways (fgsea)"))
}

#' Plot Kaplan-Meier survival curves
#'
#' Uses \pkg{survminer} when available; falls back to base-R \code{plot.survfit}
#' if \pkg{survminer} / \pkg{survMisc} are not installed.
#'
#' @param obj An \code{OmicsLens} object after \code{\link{omicslens_analyze}}.
#' @return A \code{survminer} ggsurvplot (or invisible \code{NULL} for base plot).
#' @export
omicslens_plot_survival <- function(obj) {
  if (is.null(obj$results$survival))
    stop("No survival results. Run omicslens_analyze() with time and event columns.",
         call. = FALSE)
  surv_res <- obj$results$survival

  if (requireNamespace("survminer", quietly = TRUE) &&
      requireNamespace("survMisc",  quietly = TRUE)) {
    p <- survminer::ggsurvplot(
      surv_res$fit,
      data        = surv_res$data,
      risk.table  = TRUE,
      pval        = TRUE,
      conf.int    = TRUE,
      palette     = c("#e74c3c", "#3498db"),
      ggtheme     = ggplot2::theme_minimal(),
      title       = "Kaplan-Meier: MOFA Factor High vs Low"
    )
    return(invisible(p))
  }

  # Base-R fallback
  plot(surv_res$fit,
       col     = c("#e74c3c", "#3498db"),
       lwd     = 2,
       xlab    = "Time",
       ylab    = "Survival probability",
       main    = "Kaplan-Meier: MOFA Factor High vs Low")
  legend("topright",
         legend = levels(surv_res$data$mofa_group),
         col    = c("#e74c3c", "#3498db"),
         lwd    = 2, bty = "n")
  invisible(NULL)
}
