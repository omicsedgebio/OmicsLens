# ── Internal helpers ──────────────────────────────────────────────────────────

`%||%` <- function(a, b) if (!is.null(a)) a else b

.msg <- function(verbose, ...) {
  if (verbose) message("[OmicsLens] ", ...)
}

# Read a file path, matrix, or data.frame into a plain numeric matrix.
# Row names are preserved; column names are preserved.
.read_matrix <- function(x) {
  if (is.matrix(x))     return(x)
  if (is.data.frame(x)) return(as.matrix(x))
  if (!is.character(x) || length(x) != 1)
    stop("Expected a file path, matrix, or data.frame.", call. = FALSE)
  if (!file.exists(x))
    stop("File not found: ", x, call. = FALSE)
  ext <- tolower(tools::file_ext(x))
  if (ext == "rds") {
    obj <- readRDS(x)
    return(as.matrix(obj))
  }
  sep <- if (ext %in% c("tsv", "txt")) "\t" else ","
  mat <- utils::read.csv(x, sep = sep, header = TRUE, row.names = 1,
                          check.names = FALSE)
  as.matrix(mat)
}

# Read a CSV/TSV into a data.frame (first column is NOT used as row names).
.read_table <- function(x) {
  if (is.data.frame(x)) return(x)
  if (!is.character(x) || !file.exists(x))
    stop("metadata must be a data.frame or a valid file path.", call. = FALSE)
  ext <- tolower(tools::file_ext(x))
  sep <- if (ext %in% c("tsv", "txt")) "\t" else ","
  utils::read.csv(x, sep = sep, header = TRUE, stringsAsFactors = FALSE,
                   check.names = FALSE)
}

# Detect whether a variant input looks like a MAF file.
.detect_variant_format <- function(x) {
  if (!is.character(x)) return("binary_matrix")
  ext <- tolower(tools::file_ext(x))
  if (ext == "maf") return("maf")
  if (ext %in% c("csv", "tsv", "txt", "rds")) {
    # Peek at first line to check for MAF header
    first_line <- tryCatch(readLines(x, n = 5), error = function(e) "")
    has_maf_col <- any(grepl("Hugo_Symbol", first_line))
    return(if (has_maf_col) "maf" else "binary_matrix")
  }
  "binary_matrix"
}

# Parse a MAF file into a data.frame. Skips comment lines (#).
.read_maf <- function(path) {
  if (is.data.frame(path)) return(path)
  if (inherits(path, "MAF")) return(as.data.frame(path@data))
  lines <- readLines(path, n = 20)
  skip_n <- sum(startsWith(lines, "#"))
  utils::read.table(path, sep = "\t", header = TRUE, skip = skip_n,
                     comment.char = "#", stringsAsFactors = FALSE,
                     fill = TRUE, quote = "")
}

# Convert a MAF data.frame to a binary samples-x-genes matrix.
.maf_to_binary <- function(maf_df) {
  needed <- c("Hugo_Symbol", "Tumor_Sample_Barcode")
  if (!all(needed %in% colnames(maf_df)))
    stop("MAF must contain columns: Hugo_Symbol, Tumor_Sample_Barcode.",
         call. = FALSE)
  genes   <- unique(maf_df$Hugo_Symbol)
  samples <- unique(maf_df$Tumor_Sample_Barcode)
  bin <- matrix(0L, nrow = length(samples), ncol = length(genes),
                 dimnames = list(samples, genes))
  idx <- cbind(
    match(maf_df$Tumor_Sample_Barcode, samples),
    match(maf_df$Hugo_Symbol, genes)
  )
  idx <- idx[!is.na(idx[, 1]) & !is.na(idx[, 2]), , drop = FALSE]
  bin[idx] <- 1L
  bin
}

# ── DESeq2 internal wrapper ────────────────────────────────────────────────────
.run_deseq2 <- function(obj, groups) {
  if (!requireNamespace("DESeq2", quietly = TRUE))
    stop("DESeq2 is required. Install via BiocManager::install('DESeq2').",
         call. = FALSE)
  counts <- obj$rna$counts
  meta   <- obj$metadata
  if (is.null(meta)) {
    meta <- data.frame(
      row.names  = colnames(counts),
      sample_id  = colnames(counts),
      mofa_group = groups[colnames(counts)]
    )
  } else {
    meta$mofa_group <- groups[rownames(meta)]
  }
  meta$mofa_group <- factor(meta$mofa_group)
  common <- intersect(colnames(counts), rownames(meta))
  counts <- counts[, common, drop = FALSE]
  meta   <- meta[common, , drop = FALSE]
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = round(counts),
    colData   = meta,
    design    = ~ mofa_group
  )
  dds <- DESeq2::DESeq(dds, quiet = TRUE)
  res <- DESeq2::results(dds, contrast = c("mofa_group", "High", "Low"))
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  res_df[order(res_df$padj, na.last = TRUE), ]
}

# ── fgsea internal wrapper ─────────────────────────────────────────────────────
.run_gsea <- function(de_res, organism = "human") {
  if (!requireNamespace("fgsea", quietly = TRUE))
    stop("fgsea is required. Install via BiocManager::install('fgsea').",
         call. = FALSE)
  ranked <- de_res$stat
  names(ranked) <- de_res$gene
  ranked <- ranked[!is.na(ranked)]
  ranked <- sort(ranked, decreasing = TRUE)

  # Auto-detect Ensembl IDs (ENSG...) vs gene symbols
  use_ensembl <- any(grepl("^ENSG", names(ranked)))

  gene_sets <- tryCatch({
    if (requireNamespace("msigdbr", quietly = TRUE)) {
      org_name <- if (organism == "human") "Homo sapiens" else "Mus musculus"
      m_df <- msigdbr::msigdbr(species = org_name, collection = "H")
      if (use_ensembl) {
        m_df <- m_df[!is.na(m_df$ensembl_gene) & nchar(m_df$ensembl_gene) > 0, ]
        split(m_df$ensembl_gene, m_df$gs_name)
      } else {
        split(m_df$gene_symbol, m_df$gs_name)
      }
    } else {
      .minimal_hallmark_sets()
    }
  }, error = function(e) .minimal_hallmark_sets())

  fgsea::fgsea(
    pathways  = gene_sets,
    stats     = ranked,
    minSize   = 15,
    maxSize   = 500,
    eps       = 0
  )
}

# Minimal fallback hallmark gene sets when msigdbr is unavailable.
.minimal_hallmark_sets <- function() {
  list(
    HALLMARK_E2F_TARGETS      = c("E2F1","E2F2","E2F3","PCNA","MCM2","MCM3"),
    HALLMARK_MYC_TARGETS_V1   = c("MYC","NPM1","NOP56","MRPL45","PA2G4"),
    HALLMARK_INFLAMMATORY     = c("TNF","IL6","IL1B","CXCL8","NFKB1","STAT3"),
    HALLMARK_APOPTOSIS        = c("TP53","BCL2","BAX","CASP3","CASP9","PARP1"),
    HALLMARK_HYPOXIA          = c("HIF1A","VEGFA","LDHA","ENO1","BNIP3","PDK1")
  )
}

# ── DMRcate internal wrapper ───────────────────────────────────────────────────
.run_dmr <- function(obj, groups) {
  m_vals <- obj$methylation$m_values
  common <- intersect(colnames(m_vals), names(groups))
  m_vals <- m_vals[, common, drop = FALSE]
  g      <- factor(groups[common])
  lvls   <- levels(g)

  # Try DMRcate if available (requires probe names matching a known array)
  if (requireNamespace("DMRcate", quietly = TRUE) &&
      requireNamespace("limma",   quietly = TRUE)) {
    design_mat <- stats::model.matrix(~ g)
    for (atype in c("450K", "EPIC", "27K")) {
      annot <- tryCatch(
        DMRcate::cpg.annotate(
          datatype      = "array",
          object        = m_vals,
          what          = "M",
          arraytype     = atype,
          analysis.type = "differential",
          design        = design_mat,
          coef          = 2
        ),
        error   = function(e) NULL,
        warning = function(w) NULL
      )
      if (!is.null(annot)) {
        dmr_res <- tryCatch(DMRcate::dmrcate(annot, lambda = 1000, C = 2),
                            error = function(e) NULL)
        if (!is.null(dmr_res))
          return(as.data.frame(DMRcate::extractRanges(dmr_res)))
      }
    }
  }

  # Base-R fallback: Welch t-test per CpG, BH correction
  idx1 <- which(g == lvls[1])
  idx2 <- which(g == lvls[2])
  pvals <- apply(m_vals, 1, function(x) {
    tryCatch(stats::t.test(x[idx1], x[idx2])$p.value, error = function(e) NA_real_)
  })
  mean1  <- rowMeans(m_vals[, idx1, drop = FALSE], na.rm = TRUE)
  mean2  <- rowMeans(m_vals[, idx2, drop = FALSE], na.rm = TRUE)
  out <- data.frame(
    CpG      = rownames(m_vals),
    MeanM_grp1 = round(mean1, 4),
    MeanM_grp2 = round(mean2, 4),
    MeanDiff = round(mean1 - mean2, 4),
    P.Value  = pvals,
    adj.P.Val = stats::p.adjust(pvals, method = "BH"),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  colnames(out)[2] <- paste0("MeanM_", lvls[1])
  colnames(out)[3] <- paste0("MeanM_", lvls[2])
  out <- out[order(out$P.Value), ]
  sig <- out[!is.na(out$adj.P.Val) & out$adj.P.Val < 0.05, ]
  if (nrow(sig) > 0) sig else utils::head(out, 50)
}

# ── Survival internal wrapper ─────────────────────────────────────────────────
.run_survival <- function(metadata, groups, time_col, event_col) {
  if (!requireNamespace("survival", quietly = TRUE))
    stop("survival package is required.", call. = FALSE)
  common <- intersect(rownames(metadata), names(groups))
  df <- metadata[common, , drop = FALSE]
  df$mofa_group <- factor(groups[common])
  df$time  <- as.numeric(df[[time_col]])
  df$event <- as.numeric(df[[event_col]])
  df <- df[!is.na(df$time) & !is.na(df$event), ]
  surv_obj <- survival::Surv(time = df$time, event = df$event)
  fit <- survival::survfit(surv_obj ~ mofa_group, data = df)
  cox <- survival::coxph(surv_obj ~ mofa_group, data = df)
  list(fit = fit, cox = cox, data = df)
}

# ── S3 methods ────────────────────────────────────────────────────────────────

#' Print an OmicsLens object
#' @param x An OmicsLens object.
#' @param ... Ignored.
#' @export
print.OmicsLens <- function(x, ...) {
  cat("── OmicsLens object ────────────────────────────────────────\n")
  layers <- x$layers_present
  if (length(layers) == 0) {
    cat("  No layers loaded.\n")
    return(invisible(x))
  }
  cat(sprintf("  Layers (%d): %s\n", length(layers), paste(layers, collapse = ", ")))
  if (!is.null(x$rna))         cat(sprintf("  RNA        : %d genes × %d samples\n",
                                             nrow(x$rna$counts), ncol(x$rna$counts)))
  if (!is.null(x$variants))    cat(sprintf("  Variants   : %d samples × %d genes\n",
                                             nrow(x$variants$binary), ncol(x$variants$binary)))
  if (!is.null(x$methylation)) cat(sprintf("  Methylation: %d CpGs × %d samples\n",
                                             nrow(x$methylation$beta), ncol(x$methylation$beta)))
  if (!is.null(x$metadata))    cat(sprintf("  Metadata   : %d samples × %d columns\n",
                                             nrow(x$metadata), ncol(x$metadata)))
  if (!is.null(x$mofa_model)) {
    nf <- ncol(x$results$mofa$factors)
    cat(sprintf("  MOFA2      : trained (%d factors)\n", nf))
  } else {
    cat("  MOFA2      : not yet run\n")
  }
  res_names <- names(x$results)[names(x$results) != "mofa"]
  if (length(res_names) > 0)
    cat(sprintf("  Results    : %s\n", paste(res_names, collapse = ", ")))
  cat("────────────────────────────────────────────────────────────\n")
  invisible(x)
}

#' Summarise an OmicsLens object
#' @param object An OmicsLens object.
#' @param ... Ignored.
#' @export
summary.OmicsLens <- function(object, ...) {
  print(object)
  if (!is.null(object$results$de)) {
    sig <- sum(object$results$de$padj < 0.05, na.rm = TRUE)
    cat(sprintf("  DE genes (padj<0.05): %d\n", sig))
  }
  if (!is.null(object$results$gsea)) {
    sig <- sum(object$results$gsea$padj < 0.05, na.rm = TRUE)
    cat(sprintf("  GSEA pathways (padj<0.05): %d\n", sig))
  }
  if (!is.null(object$results$dmr))
    cat(sprintf("  DMRs identified: %d\n", nrow(object$results$dmr)))
  invisible(object)
}
