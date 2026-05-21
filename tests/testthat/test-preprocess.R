test_that("omicslens_preprocess filters low-count RNA genes", {
  # Set ~half the genes to counts of 0 to trigger filtering
  rna_sparse              <- .rna_mat
  rna_sparse[1:50, ]      <- 0L
  obj <- omicslens_load(rna_counts = rna_sparse, verbose = FALSE)
  skip_if_not_installed("DESeq2")
  obj2 <- omicslens_preprocess(obj, verbose = FALSE)
  expect_lt(nrow(obj2$rna$counts), nrow(obj$rna$counts))
})

test_that("omicslens_preprocess creates normalised matrix", {
  obj <- omicslens_load(rna_counts = .rna_mat, verbose = FALSE)
  skip_if_not_installed("DESeq2")
  obj2 <- omicslens_preprocess(obj, verbose = FALSE)
  expect_false(is.null(obj2$rna$normalized))
  expect_equal(ncol(obj2$rna$normalized), .n_samples)
})

test_that("omicslens_preprocess computes methylation M-values", {
  obj <- omicslens_load(methylation = .meth_mat, verbose = FALSE)
  obj2 <- omicslens_preprocess(obj, meth_var_top_n = 100L, verbose = FALSE)
  expect_false(is.null(obj2$methylation$m_values))
  expect_lte(nrow(obj2$methylation$m_values), 100L)
})

test_that("omicslens_preprocess M-values are finite", {
  obj  <- omicslens_load(methylation = .meth_mat, verbose = FALSE)
  obj2 <- omicslens_preprocess(obj, meth_var_top_n = 100L, verbose = FALSE)
  expect_true(all(is.finite(obj2$methylation$m_values)))
})

test_that("omicslens_preprocess filters variants by MAF", {
  # Force all genes to MAF < 0.05 except one
  var_sparse      <- matrix(0L, nrow = .n_samples, ncol = .n_genes_var,
                              dimnames = list(.sample_ids,
                                              paste0("Gene", seq_len(.n_genes_var))))
  var_sparse[, 1] <- 1L  # Gene1: MAF = 1.0
  obj  <- omicslens_load(variants = var_sparse,
                          variant_format = "binary_matrix",
                          verbose = FALSE)
  obj2 <- omicslens_preprocess(obj, variant_min_maf = 0.05, verbose = FALSE)
  expect_equal(ncol(obj2$variants$binary), 1L)
})

test_that("omicslens_preprocess errors if called on non-OmicsLens", {
  expect_error(omicslens_preprocess(list()), "OmicsLens")
})
