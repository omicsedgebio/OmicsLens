#' OmicsLens: Multi-Omics Integration Pipeline with Interactive Visualization
#'
#' @description
#' OmicsLens provides an end-to-end pipeline for integrating RNA-Seq,
#' whole-genome sequencing (WGS), and DNA methylation data. It wraps MOFA2
#' latent factor analysis behind a clean five-function API, adds downstream
#' analyses (DESeq2, fgsea, DMRcate, survival), and presents all results in
#' an interactive Shiny dashboard and a reproducible HTML report.
#'
#' @section Workflow:
#' \enumerate{
#'   \item \code{\link{omicslens_load}}        — load and align multi-omics data
#'   \item \code{\link{omicslens_preprocess}}  — normalize and filter each layer
#'   \item \code{\link{omicslens_integrate}}   — MOFA2 factor analysis
#'   \item \code{\link{omicslens_analyze}}     — DE, GSEA, DMR, survival
#'   \item \code{\link{omicslens_app}}         — launch the Shiny dashboard
#'   \item \code{\link{omicslens_report}}      — render an HTML report
#' }
#'
#' @section Installation:
#' \preformatted{
#' # From GitHub (recommended for latest version):
#' # install.packages("devtools")
#' devtools::install_github("priyanshpathak/OmicsLens")
#'
#' # Bioconductor dependencies must be installed first:
#' # if (!requireNamespace("BiocManager", quietly = TRUE))
#' #   install.packages("BiocManager")
#' # BiocManager::install(c("MOFA2","DESeq2","fgsea","DMRcate"))
#' }
#'
#' @author Priyansh Pathak \email{priyanshpt84@@gmail.com}
#' @references
#' Argelaguet R et al. (2020). MOFA+: a statistical framework for
#' comprehensive integration of multi-modal single-cell data.
#' \emph{Genome Biology}, 21:111.
#'
#' @keywords internal
"_PACKAGE"
