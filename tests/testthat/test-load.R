test_that("omicslens_load accepts in-memory matrices", {
  obj <- omicslens_load(rna_counts = .rna_mat, verbose = FALSE)
  expect_s3_class(obj, "OmicsLens")
  expect_true("rna" %in% obj$layers_present)
  expect_equal(dim(obj$rna$counts), c(.n_genes, .n_samples))
})

test_that("omicslens_load accepts methylation matrix", {
  obj <- omicslens_load(methylation = .meth_mat, verbose = FALSE)
  expect_s3_class(obj, "OmicsLens")
  expect_true("methylation" %in% obj$layers_present)
  expect_equal(dim(obj$methylation$beta), c(.n_cpgs, .n_samples))
})

test_that("omicslens_load accepts variant binary matrix", {
  obj <- omicslens_load(variants = .var_mat,
                         variant_format = "binary_matrix",
                         verbose = FALSE)
  expect_s3_class(obj, "OmicsLens")
  expect_true("variants" %in% obj$layers_present)
})

test_that("omicslens_load attaches metadata", {
  obj <- omicslens_load(rna_counts = .rna_mat,
                         metadata   = .meta_df,
                         verbose    = FALSE)
  expect_false(is.null(obj$metadata))
  expect_true("sample_id" %in% colnames(obj$metadata))
})

test_that("omicslens_load errors on missing sample_id column", {
  bad_meta <- data.frame(sid = .sample_ids)
  expect_error(
    omicslens_load(rna_counts = .rna_mat, metadata = bad_meta, verbose = FALSE),
    "sample_id"
  )
})

test_that("omicslens_load warns on mismatched sample IDs", {
  rna_sub <- .rna_mat[, 1:20]  # only 20 of 30 samples
  expect_warning(
    omicslens_load(rna_counts = rna_sub,
                    methylation = .meth_mat,
                    verbose = FALSE),
    "mismatch"
  )
})

test_that("omicslens_load errors when all layers are NULL", {
  expect_error(omicslens_load(verbose = FALSE), "at least one")
})

test_that("omicslens_load errors on negative RNA counts", {
  bad_rna <- .rna_mat
  bad_rna[1, 1] <- -1L
  expect_error(
    omicslens_load(rna_counts = bad_rna, verbose = FALSE),
    "negative"
  )
})

test_that("omicslens_load errors on methylation values outside [0,1]", {
  bad_meth <- .meth_mat
  bad_meth[1, 1] <- 1.5
  expect_error(
    omicslens_load(methylation = bad_meth, verbose = FALSE),
    "[0, 1]", fixed = TRUE
  )
})

test_that("print.OmicsLens produces output", {
  obj <- omicslens_load(rna_counts = .rna_mat, verbose = FALSE)
  expect_output(print(obj), "OmicsLens")
})

test_that("summary.OmicsLens produces output", {
  obj <- omicslens_load(rna_counts = .rna_mat, verbose = FALSE)
  expect_output(summary(obj), "OmicsLens")
})
