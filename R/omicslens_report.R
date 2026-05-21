#' Generate a reproducible HTML analysis report
#'
#' Renders a parameterised R Markdown document that summarises all analyses
#' stored in an \code{OmicsLens} object: data overview, MOFA2 factors,
#' differential expression, GSEA, DMRs, and survival curves.
#'
#' @param obj An \code{OmicsLens} object (typically after running
#'   \code{\link{omicslens_analyze}}).
#' @param output_file Path to the output HTML file. Defaults to
#'   \code{"OmicsLens_report.html"} in the current working directory.
#' @param title Character; report title (default \code{"OmicsLens Analysis Report"}).
#' @param author Character; author name to display in the report header.
#' @param open Logical; open the rendered report in the default browser
#'   (default \code{TRUE}).
#' @param ... Additional arguments passed to \code{rmarkdown::render()}.
#'
#' @return Invisibly returns the path to the generated HTML file.
#'
#' @examples
#' \dontrun{
#' omicslens_report(obj, output_file = "my_report.html",
#'                  title = "TCGA-BRCA Multi-Omics", author = "Jane Smith")
#' }
#'
#' @seealso \code{\link{omicslens_analyze}}, \code{\link{omicslens_app}}
#' @export
omicslens_report <- function(obj,
                              output_file = "OmicsLens_report.html",
                              title       = "OmicsLens Analysis Report",
                              author      = "",
                              open        = TRUE,
                              ...) {
  if (!inherits(obj, "OmicsLens"))
    stop("obj must be an OmicsLens object.", call. = FALSE)
  if (!requireNamespace("rmarkdown", quietly = TRUE))
    stop("rmarkdown is required: install.packages('rmarkdown')", call. = FALSE)

  template <- system.file("rmarkdown", "report_template.Rmd",
                            package = "OmicsLens")
  if (!nzchar(template) || !file.exists(template))
    stop("Report template not found. Reinstall OmicsLens.", call. = FALSE)

  # Write the OmicsLens object to a temp RDS so the Rmd can load it
  tmp_rds <- tempfile(fileext = ".rds")
  saveRDS(obj, tmp_rds)
  on.exit(unlink(tmp_rds), add = TRUE)

  output_file <- normalizePath(output_file, mustWork = FALSE)

  rmarkdown::render(
    input       = template,
    output_file = output_file,
    params      = list(
      obj_path = tmp_rds,
      title    = title,
      author   = author
    ),
    envir = new.env(parent = globalenv()),
    quiet = TRUE,
    ...
  )

  message("[OmicsLens] Report saved to: ", output_file)
  if (open) utils::browseURL(output_file)
  invisible(output_file)
}
