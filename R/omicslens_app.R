#' Launch the OmicsLens interactive Shiny dashboard
#'
#' Opens a multi-tab Shiny dashboard that lets you explore all results stored
#' in an \code{OmicsLens} object: MOFA2 factor plots, differential expression
#' volcano plots, GSEA dot plots, DMR tables, Kaplan-Meier curves, and
#' download buttons for every result table.
#'
#' @param obj An \code{OmicsLens} object. Results from
#'   \code{\link{omicslens_analyze}} are needed for most tabs, but the app
#'   will open with whatever is available.
#' @param port Integer port for the Shiny server (default \code{NULL}, which
#'   lets Shiny pick an available port).
#' @param launch_browser Logical; open a browser window automatically
#'   (default \code{TRUE}).
#' @param ... Additional arguments passed to \code{shiny::runApp()}.
#'
#' @return Invisible \code{NULL}. The function blocks until the app is closed.
#'
#' @examples
#' \dontrun{
#' # After running the full pipeline:
#' omicslens_app(obj)
#' }
#'
#' @seealso \code{\link{omicslens_analyze}}, \code{\link{omicslens_report}}
#' @export
omicslens_app <- function(obj, port = NULL, launch_browser = TRUE, ...) {
  if (!inherits(obj, "OmicsLens"))
    stop("obj must be an OmicsLens object.", call. = FALSE)
  if (!requireNamespace("shiny", quietly = TRUE))
    stop("shiny is required: install.packages('shiny')", call. = FALSE)
  if (!requireNamespace("shinydashboard", quietly = TRUE))
    stop("shinydashboard is required: install.packages('shinydashboard')",
         call. = FALSE)

  # Pass the OmicsLens object to the Shiny session via an environment.
  .omicslens_env <- new.env(parent = emptyenv())
  .omicslens_env$obj <- obj

  ui     <- .build_ui()
  server <- .build_server(.omicslens_env)

  app <- shiny::shinyApp(ui = ui, server = server)
  shiny::runApp(app, port = port, launch.browser = launch_browser, ...)
}

# ── UI ─────────────────────────────────────────────────────────────────────────
.build_ui <- function() {
  shinydashboard::dashboardPage(
    skin = "blue",

    shinydashboard::dashboardHeader(
      title = shiny::tags$span(
        shiny::tags$img(
          src = "https://img.icons8.com/fluency/48/lens.png",
          height = "30px", style = "margin-right:8px;"
        ),
        "OmicsLens"
      )
    ),

    shinydashboard::dashboardSidebar(
      shinydashboard::sidebarMenu(
        id = "tabs",
        shinydashboard::menuItem("Overview",        tabName = "overview",   icon = shiny::icon("table")),
        shinydashboard::menuItem("MOFA2 Factors",   tabName = "mofa",       icon = shiny::icon("project-diagram")),
        shinydashboard::menuItem("Diff. Expression",tabName = "de",         icon = shiny::icon("dna")),
        shinydashboard::menuItem("Pathway Enrichment", tabName = "gsea",    icon = shiny::icon("route")),
        shinydashboard::menuItem("Methylation/DMR", tabName = "dmr",        icon = shiny::icon("circle-nodes")),
        shinydashboard::menuItem("Survival",         tabName = "survival",  icon = shiny::icon("heartbeat")),
        shinydashboard::menuItem("Export",           tabName = "export",    icon = shiny::icon("download"))
      )
    ),

    shinydashboard::dashboardBody(
      shiny::tags$head(
        shiny::tags$style(shiny::HTML("
          .content-wrapper, .right-side { background-color: #f4f6f9; }
          .box { border-top: 3px solid #3c8dbc; }
          .nav-tabs > li.active > a { color: #3c8dbc; font-weight: bold; }
          h4 { color: #2c3e50; }
        "))
      ),
      shinydashboard::tabItems(
        .tab_overview(),
        .tab_mofa(),
        .tab_de(),
        .tab_gsea(),
        .tab_dmr(),
        .tab_survival(),
        .tab_export()
      )
    )
  )
}

# ── Tab definitions ────────────────────────────────────────────────────────────
.tab_overview <- function() {
  shinydashboard::tabItem("overview",
    shiny::fluidRow(
      shinydashboard::valueBoxOutput("box_layers", width = 4),
      shinydashboard::valueBoxOutput("box_samples", width = 4),
      shinydashboard::valueBoxOutput("box_genes",   width = 4)
    ),
    shiny::fluidRow(
      shinydashboard::box(title = "Data Summary", width = 6, solidHeader = TRUE,
        DT::DTOutput("tbl_summary")),
      shinydashboard::box(title = "Sample PCA (RNA, if available)", width = 6, solidHeader = TRUE,
        shiny::plotOutput("plt_pca", height = "350px"))
    )
  )
}

.tab_mofa <- function() {
  shinydashboard::tabItem("mofa",
    shiny::fluidRow(
      shinydashboard::box(title = "Variance Explained per Factor", width = 6, solidHeader = TRUE,
        shiny::plotOutput("plt_var_exp", height = "350px")),
      shinydashboard::box(title = "Factor Correlation", width = 6, solidHeader = TRUE,
        shiny::plotOutput("plt_factor_cor", height = "350px"))
    ),
    shiny::fluidRow(
      shinydashboard::box(title = "Factor Scatter Plot", width = 6, solidHeader = TRUE,
        shiny::uiOutput("factor_x_ui"),
        shiny::uiOutput("factor_y_ui"),
        shiny::plotOutput("plt_factors_scatter", height = "320px")),
      shinydashboard::box(title = "Top Feature Weights", width = 6, solidHeader = TRUE,
        shiny::uiOutput("weight_view_ui"),
        shiny::uiOutput("weight_factor_ui"),
        shiny::plotOutput("plt_weights", height = "320px"))
    )
  )
}

.tab_de <- function() {
  shinydashboard::tabItem("de",
    shiny::fluidRow(
      shinydashboard::box(title = "Volcano Plot", width = 6, solidHeader = TRUE,
        shiny::sliderInput("de_pval",  "padj threshold",  0, 0.5, 0.05, 0.01),
        shiny::sliderInput("de_lfc",   "|log2FC| cutoff", 0, 5,   1.0,  0.1),
        shiny::plotOutput("plt_volcano", height = "380px")),
      shinydashboard::box(title = "Top DE Genes Heatmap", width = 6, solidHeader = TRUE,
        shiny::numericInput("de_top_n", "Top N genes", 30, 5, 100),
        shiny::plotOutput("plt_de_heatmap", height = "380px"))
    ),
    shiny::fluidRow(
      shinydashboard::box(title = "DE Results Table", width = 12, solidHeader = TRUE,
        DT::DTOutput("tbl_de"))
    )
  )
}

.tab_gsea <- function() {
  shinydashboard::tabItem("gsea",
    shiny::fluidRow(
      shinydashboard::box(title = "GSEA Dot Plot", width = 7, solidHeader = TRUE,
        shiny::numericInput("gsea_top_n", "Top N pathways", 20, 5, 50),
        shiny::plotOutput("plt_gsea_dot", height = "520px")),
      shinydashboard::box(title = "GSEA Results", width = 5, solidHeader = TRUE,
        DT::DTOutput("tbl_gsea"))
    )
  )
}

.tab_dmr <- function() {
  shinydashboard::tabItem("dmr",
    shiny::fluidRow(
      shinydashboard::box(title = "Differentially Methylated CpGs / Regions", width = 12, solidHeader = TRUE,
        DT::DTOutput("tbl_dmr"))
    )
  )
}

.tab_survival <- function() {
  shinydashboard::tabItem("survival",
    shiny::fluidRow(
      shinydashboard::box(title = "Kaplan-Meier Curves", width = 7, solidHeader = TRUE,
        shiny::plotOutput("plt_km", height = "480px")),
      shinydashboard::box(title = "Cox Regression Summary", width = 5, solidHeader = TRUE,
        shiny::verbatimTextOutput("txt_cox"))
    )
  )
}

.tab_export <- function() {
  shinydashboard::tabItem("export",
    shiny::fluidRow(
      shinydashboard::box(title = "Download Results", width = 6, solidHeader = TRUE,
        shiny::p("Download individual analysis results as CSV files."),
        shiny::downloadButton("dl_de",      "Differential Expression"),
        shiny::br(), shiny::br(),
        shiny::downloadButton("dl_gsea",    "GSEA Pathways"),
        shiny::br(), shiny::br(),
        shiny::downloadButton("dl_dmr",     "DMR Results"),
        shiny::br(), shiny::br(),
        shiny::downloadButton("dl_factors", "MOFA2 Factor Scores")
      ),
      shinydashboard::box(title = "Generate HTML Report", width = 6, solidHeader = TRUE,
        shiny::p("Render a full reproducible HTML report of all analyses."),
        shiny::textInput("report_title",  "Report title", "OmicsLens Analysis Report"),
        shiny::textInput("report_author", "Author name",  ""),
        shiny::actionButton("btn_report",  "Generate Report",
                             icon = shiny::icon("file-alt"), class = "btn-primary"),
        shiny::br(), shiny::br(),
        shiny::uiOutput("report_download_ui")
      )
    )
  )
}

# ── Server ─────────────────────────────────────────────────────────────────────
.build_server <- function(env) {
  function(input, output, session) {
    obj <- env$obj

    # Compute DMR on-the-fly if missing and methylation data is present
    if (is.null(obj$results$dmr) &&
        "methylation" %in% obj$layers_present &&
        !is.null(obj$methylation$m_values) &&
        !is.null(obj$results$groups)) {
      obj$results$dmr <- tryCatch(
        .run_dmr(obj, obj$results$groups),
        error = function(e) NULL
      )
    }

    # ── Overview ───────────────────────────────────────────────────────────────
    output$box_layers <- shinydashboard::renderValueBox(
      shinydashboard::valueBox(
        length(obj$layers_present), "Data Layers", icon = shiny::icon("layer-group"),
        color = "blue"
      )
    )
    output$box_samples <- shinydashboard::renderValueBox({
      n <- if (!is.null(obj$rna)) ncol(obj$rna$counts) else
             if (!is.null(obj$methylation)) ncol(obj$methylation$beta) else 0
      shinydashboard::valueBox(n, "Samples", icon = shiny::icon("users"), color = "green")
    })
    output$box_genes <- shinydashboard::renderValueBox({
      n <- if (!is.null(obj$rna)) nrow(obj$rna$counts) else 0
      shinydashboard::valueBox(n, "Genes (RNA)", icon = shiny::icon("dna"), color = "orange")
    })

    output$tbl_summary <- DT::renderDT({
      rows <- list()
      if ("rna" %in% obj$layers_present)
        rows[[length(rows)+1]] <- data.frame(
          Layer="RNA-Seq", Features=nrow(obj$rna$counts), Samples=ncol(obj$rna$counts))
      if ("variants" %in% obj$layers_present)
        rows[[length(rows)+1]] <- data.frame(
          Layer="Variants", Features=ncol(obj$variants$binary), Samples=nrow(obj$variants$binary))
      if ("methylation" %in% obj$layers_present)
        rows[[length(rows)+1]] <- data.frame(
          Layer="Methylation", Features=nrow(obj$methylation$beta), Samples=ncol(obj$methylation$beta))
      DT::datatable(do.call(rbind, rows), options = list(dom = "t"), rownames = FALSE)
    })

    output$plt_pca <- shiny::renderPlot({
      if (is.null(obj$rna$normalized)) {
        graphics::plot.new()
        graphics::text(0.5, 0.5, "No RNA data available", cex = 1.4, col = "grey40")
        return(invisible(NULL))
      }
      mat <- obj$rna$normalized
      pca <- stats::prcomp(t(mat), scale. = FALSE, center = TRUE)
      df  <- as.data.frame(pca$x[, 1:2])
      df$Sample <- rownames(df)
      pct <- round(summary(pca)$importance[2, 1:2] * 100, 1)
      df$Group <- if (!is.null(obj$results$groups)) {
        obj$results$groups[df$Sample]
      } else {
        rep("All", nrow(df))
      }
      ggplot2::ggplot(df, ggplot2::aes(x = PC1, y = PC2, color = Group)) +
        ggplot2::geom_point(size = 3, alpha = 0.85) +
        ggplot2::scale_color_brewer(palette = "Set1", name = "Group") +
        ggplot2::labs(x = paste0("PC1 (", pct[1], "%)"),
                       y = paste0("PC2 (", pct[2], "%)"),
                       title = "RNA-Seq PCA") +
        ggplot2::theme_minimal(base_size = 13)
    })

    # ── MOFA2 ──────────────────────────────────────────────────────────────────
    output$plt_var_exp <- shiny::renderPlot({
      if (is.null(obj$results$mofa)) {
        graphics::plot.new()
        graphics::text(0.5, 0.5, "Run omicslens_integrate() first", cex = 1.3, col = "grey40")
        return(invisible(NULL))
      }
      ve <- obj$results$mofa$variance_explained$r2_per_factor[[1]]
      df <- as.data.frame(ve)
      df$Factor <- factor(paste0("F", seq_len(nrow(df))), levels = paste0("F", seq_len(nrow(df))))
      df_long <- tidyr::pivot_longer(df, -Factor, names_to = "View", values_to = "VarExp")
      ggplot2::ggplot(df_long, ggplot2::aes(x = Factor, y = VarExp, fill = View)) +
        ggplot2::geom_bar(stat = "identity", position = "dodge") +
        ggplot2::scale_fill_brewer(palette = "Set2") +
        ggplot2::labs(y = "Variance Explained (%)",
                       title = "Variance explained per factor per view") +
        ggplot2::theme_minimal(base_size = 13)
    })

    output$plt_factor_cor <- shiny::renderPlot({
      if (is.null(obj$results$mofa)) {
        graphics::plot.new()
        graphics::text(0.5, 0.5, "Run omicslens_integrate() first", cex = 1.3, col = "grey40")
        return(invisible(NULL))
      }
      mat <- obj$results$mofa$factors
      pheatmap::pheatmap(stats::cor(mat), display_numbers = TRUE,
                          main = "Factor Correlation",
                          color = grDevices::colorRampPalette(c("#2166ac","white","#d73027"))(50))
      invisible(NULL)
    })

    output$factor_x_ui <- shiny::renderUI({
      if (is.null(obj$results$mofa)) return(NULL)
      fn <- colnames(obj$results$mofa$factors)
      if (is.null(fn)) fn <- paste0("Factor", seq_len(ncol(obj$results$mofa$factors)))
      shiny::selectInput("factor_x", "X-axis", fn, fn[1])
    })
    output$factor_y_ui <- shiny::renderUI({
      if (is.null(obj$results$mofa)) return(NULL)
      fn <- colnames(obj$results$mofa$factors)
      if (is.null(fn)) fn <- paste0("Factor", seq_len(ncol(obj$results$mofa$factors)))
      shiny::selectInput("factor_y", "Y-axis", fn, fn[min(2, length(fn))])
    })

    output$plt_factors_scatter <- shiny::renderPlot({
      shiny::req(obj$results$mofa, input$factor_x, input$factor_y)
      mat <- as.data.frame(obj$results$mofa$factors)
      mat$Sample <- rownames(mat)
      mat$Group <- if (!is.null(obj$results$groups)) {
        obj$results$groups[mat$Sample]
      } else {
        rep("All", nrow(mat))
      }
      fx <- input$factor_x
      fy <- input$factor_y
      if (!fx %in% colnames(mat) || !fy %in% colnames(mat)) return(invisible(NULL))
      ggplot2::ggplot(mat, ggplot2::aes(x = .data[[fx]], y = .data[[fy]],
                                         color = Group)) +
        ggplot2::geom_point(size = 3, alpha = 0.85) +
        ggplot2::scale_color_manual(values = c("High" = "#e74c3c", "Low" = "#3498db",
                                                "All" = "#7f8c8d"),
                                     name = "Group") +
        ggplot2::labs(title = paste("MOFA2:", fx, "vs", fy),
                       x = fx, y = fy) +
        ggplot2::theme_minimal(base_size = 13)
    })

    output$weight_view_ui <- shiny::renderUI({
      if (is.null(obj$results$mofa)) return(NULL)
      shiny::selectInput("weight_view", "View", names(obj$results$mofa$weights))
    })
    output$weight_factor_ui <- shiny::renderUI({
      if (is.null(obj$results$mofa)) return(NULL)
      factor_names <- colnames(obj$results$mofa$factors)
      if (is.null(factor_names)) factor_names <- paste0("Factor", seq_len(ncol(obj$results$mofa$factors)))
      shiny::selectInput("weight_factor", "Factor", factor_names, factor_names[1])
    })

    output$plt_weights <- shiny::renderPlot({
      shiny::req(obj$results$mofa, input$weight_view, input$weight_factor)
      wraw <- obj$results$mofa$weights[[input$weight_view]]
      wmat <- if (is.matrix(wraw)) wraw else wraw[[1]]
      if (!input$weight_factor %in% colnames(wmat)) return(invisible(NULL))
      vals <- wmat[, input$weight_factor]
      top  <- utils::head(sort(abs(vals), decreasing = TRUE), 20)
      df   <- data.frame(Feature = names(top), Weight = vals[names(top)])
      df   <- df[order(df$Weight), ]
      df$Feature   <- factor(df$Feature, levels = df$Feature)
      df$Direction <- ifelse(df$Weight > 0, "Positive", "Negative")
      ggplot2::ggplot(df, ggplot2::aes(x = Weight, y = Feature, fill = Direction)) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::scale_fill_manual(values = c("Positive" = "#e74c3c", "Negative" = "#3498db")) +
        ggplot2::labs(x = "Weight", y = NULL,
                       title = paste("Top 20 weights:", input$weight_view, input$weight_factor)) +
        ggplot2::theme_minimal(base_size = 11) +
        ggplot2::theme(legend.position = "none")
    })

    # ── DE ─────────────────────────────────────────────────────────────────────
    output$plt_volcano <- shiny::renderPlot({
      if (is.null(obj$results$de)) {
        graphics::plot.new()
        graphics::text(0.5, 0.5, "Run omicslens_analyze() first", cex = 1.3, col = "grey40")
        return(invisible(NULL))
      }
      pval_cut <- if (is.null(input$de_pval)) 0.05 else input$de_pval
      lfc_cut  <- if (is.null(input$de_lfc))  1.0  else input$de_lfc
      de <- obj$results$de
      de$log10p <- -log10(de$padj + 1e-300)
      de$Sig    <- "NS"
      de$Sig[!is.na(de$padj) & de$padj < pval_cut & de$log2FoldChange >  lfc_cut] <- "Up"
      de$Sig[!is.na(de$padj) & de$padj < pval_cut & de$log2FoldChange < -lfc_cut] <- "Down"
      cols <- c("Up" = "#e74c3c", "Down" = "#3498db", "NS" = "#bdc3c7")
      ggplot2::ggplot(de, ggplot2::aes(x = log2FoldChange, y = log10p, color = Sig)) +
        ggplot2::geom_point(size = 1.2, alpha = 0.6) +
        ggplot2::scale_color_manual(values = cols, name = "") +
        ggplot2::geom_vline(xintercept = c(-lfc_cut, lfc_cut),
                             linetype = "dashed", color = "grey40") +
        ggplot2::geom_hline(yintercept = -log10(pval_cut),
                             linetype = "dashed", color = "grey40") +
        ggplot2::labs(x = "log2 Fold Change (High vs Low)", y = "-log10(padj)",
                       title = "Volcano Plot") +
        ggplot2::theme_minimal(base_size = 13)
    })

    output$plt_de_heatmap <- shiny::renderPlot({
      if (is.null(obj$results$de) || is.null(obj$rna$normalized)) {
        graphics::plot.new()
        graphics::text(0.5, 0.5, "No data available", cex = 1.3, col = "grey40")
        return(invisible(NULL))
      }
      top_n    <- if (is.null(input$de_top_n)) 30L else as.integer(input$de_top_n)
      top_genes <- utils::head(obj$results$de$gene[!is.na(obj$results$de$padj)], top_n)
      mat  <- obj$rna$normalized[top_genes, , drop = FALSE]
      mat  <- t(scale(t(mat)))
      ann  <- if (!is.null(obj$results$groups)) {
        data.frame(Group = obj$results$groups, row.names = names(obj$results$groups))
      } else NULL
      pheatmap::pheatmap(mat, annotation_col = ann,
                          show_colnames = FALSE,
                          color = grDevices::colorRampPalette(c("#2166ac","white","#d73027"))(100),
                          main = paste("Top", top_n, "DE genes"))
      invisible(NULL)
    })

    output$tbl_de <- DT::renderDT({
      if (is.null(obj$results$de)) return(NULL)
      df <- obj$results$de[, c("gene","baseMean","log2FoldChange","pvalue","padj"), drop = FALSE]
      df <- df[!is.na(df$padj), ]
      DT::datatable(df, rownames = FALSE,
                     options = list(pageLength = 15, scrollX = TRUE)) |>
        DT::formatRound(c("baseMean","log2FoldChange","pvalue","padj"), 4) |>
        DT::formatStyle("padj",
                          backgroundColor = DT::styleInterval(0.05, c("#abebc6","white")))
    })

    # ── GSEA ──────────────────────────────────────────────────────────────────
    output$plt_gsea_dot <- shiny::renderPlot({
      if (is.null(obj$results$gsea)) {
        graphics::plot.new()
        graphics::text(0.5, 0.5, "Run omicslens_analyze() first", cex = 1.3, col = "grey40")
        return(invisible(NULL))
      }
      top_n <- if (is.null(input$gsea_top_n)) 20L else as.integer(input$gsea_top_n)
      df   <- as.data.frame(obj$results$gsea)
      df   <- df[order(df$padj), ]
      df   <- utils::head(df, top_n)
      df$pathway <- gsub("HALLMARK_", "", df$pathway)
      df$pathway <- factor(df$pathway, levels = rev(df$pathway))
      ggplot2::ggplot(df, ggplot2::aes(x = NES, y = pathway,
                                        size = size, color = padj)) +
        ggplot2::geom_point() +
        ggplot2::scale_color_gradient(low = "#e74c3c", high = "#95a5a6",
                                       name = "padj") +
        ggplot2::scale_size_continuous(name = "Gene Set Size", range = c(3, 10)) +
        ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::labs(x = "Normalised Enrichment Score", y = NULL,
                       title = paste("Top", top_n, "Enriched Pathways"))
    })

    output$tbl_gsea <- DT::renderDT({
      if (is.null(obj$results$gsea)) return(NULL)
      df <- as.data.frame(obj$results$gsea)[, c("pathway","NES","padj","size"), drop=FALSE]
      DT::datatable(df, rownames = FALSE,
                     options = list(pageLength = 10, scrollX = TRUE)) |>
        DT::formatRound(c("NES","padj"), 4)
    })

    # ── DMR ────────────────────────────────────────────────────────────────────
    output$tbl_dmr <- DT::renderDT({
      if (is.null(obj$results$dmr))
        return(DT::datatable(
          data.frame(Message = "No methylation layer found. Load methylation data and re-run omicslens_analyze()."),
          options = list(dom = "t"), rownames = FALSE
        ))
      df  <- obj$results$dmr
      num_cols <- names(df)[vapply(df, is.numeric, logical(1))]
      tbl <- DT::datatable(df, rownames = FALSE, filter = "top",
                            options = list(pageLength = 15, scrollX = TRUE))
      if (length(num_cols) > 0)
        tbl <- DT::formatRound(tbl, num_cols, digits = 4)
      tbl
    })

    # ── Survival ───────────────────────────────────────────────────────────────
    output$plt_km <- shiny::renderPlot({
      if (is.null(obj$results$survival)) {
        graphics::plot.new()
        graphics::text(0.5, 0.5, "No survival data.\nAdd 'time_os' and 'event' columns\nto metadata.",
                        cex = 1.3, col = "grey40")
        return(invisible(NULL))
      }
      surv_res <- obj$results$survival
      if (requireNamespace("survminer", quietly = TRUE) &&
          requireNamespace("survMisc",  quietly = TRUE)) {
        p <- survminer::ggsurvplot(
          surv_res$fit, data = surv_res$data,
          risk.table    = TRUE,
          pval          = TRUE,
          conf.int      = TRUE,
          palette       = c("#e74c3c", "#3498db"),
          ggtheme       = ggplot2::theme_minimal(),
          title         = "Kaplan-Meier: MOFA Factor High vs Low",
          risk.table.height = 0.25
        )
        print(p)
      } else {
        plot(surv_res$fit,
             col  = c("#e74c3c", "#3498db"),
             lwd  = 2,
             xlab = "Time",
             ylab = "Survival probability",
             main = "Kaplan-Meier: MOFA Factor High vs Low")
        legend("topright",
               legend = levels(surv_res$data$mofa_group),
               col    = c("#e74c3c", "#3498db"),
               lwd = 2, bty = "n")
      }
    })

    output$txt_cox <- shiny::renderPrint({
      if (is.null(obj$results$survival)) {
        cat("No survival results available.")
        return(invisible(NULL))
      }
      cox <- obj$results$survival$cox
      if (is.null(cox)) {
        cat("No Cox model available.")
        return(invisible(NULL))
      }
      print(summary(cox))
    })

    # ── Export ─────────────────────────────────────────────────────────────────
    output$dl_de <- shiny::downloadHandler(
      filename = "OmicsLens_DE_results.csv",
      content  = function(f) {
        if (!is.null(obj$results$de)) utils::write.csv(obj$results$de, f, row.names = FALSE)
      }
    )
    output$dl_gsea <- shiny::downloadHandler(
      filename = "OmicsLens_GSEA_results.csv",
      content  = function(f) {
        if (!is.null(obj$results$gsea))
          utils::write.csv(as.data.frame(obj$results$gsea), f, row.names = FALSE)
      }
    )
    output$dl_dmr <- shiny::downloadHandler(
      filename = "OmicsLens_DMR_results.csv",
      content  = function(f) {
        if (!is.null(obj$results$dmr)) utils::write.csv(obj$results$dmr, f, row.names = FALSE)
      }
    )
    output$dl_factors <- shiny::downloadHandler(
      filename = "OmicsLens_MOFA_factors.csv",
      content  = function(f) {
        if (!is.null(obj$results$mofa))
          utils::write.csv(obj$results$mofa$factors, f)
      }
    )

    report_path <- shiny::reactiveVal(NULL)
    shiny::observeEvent(input$btn_report, {
      out_file <- tempfile(fileext = ".html")
      tryCatch({
        omicslens_report(obj, output_file = out_file,
                          title  = input$report_title,
                          author = input$report_author,
                          open   = FALSE)
        report_path(out_file)
        shiny::showNotification("Report generated successfully.", type = "message")
      }, error = function(e) {
        shiny::showNotification(paste("Report error:", e$message), type = "error")
      })
    })

    output$report_download_ui <- shiny::renderUI({
      if (is.null(report_path())) return(NULL)
      shiny::downloadButton("dl_report", "Download HTML Report")
    })
    output$dl_report <- shiny::downloadHandler(
      filename = "OmicsLens_report.html",
      content  = function(f) file.copy(report_path(), f)
    )
  }
}
