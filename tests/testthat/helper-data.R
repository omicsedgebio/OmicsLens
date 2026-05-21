# Shared synthetic data for all OmicsLens tests
set.seed(42)

.n_genes   <- 200L
.n_samples <- 30L
.n_cpgs    <- 500L
.n_genes_var <- 80L

.sample_ids <- paste0("S", seq_len(.n_samples))
.gene_ids   <- paste0("Gene", seq_len(.n_genes))
.cpg_ids    <- paste0("cg", seq_len(.n_cpgs))

# RNA count matrix (genes x samples, Poisson counts)
.rna_mat <- matrix(
  rpois(.n_genes * .n_samples, lambda = 100),
  nrow     = .n_genes,
  ncol     = .n_samples,
  dimnames = list(.gene_ids, .sample_ids)
)

# Methylation beta matrix (CpGs x samples, Beta distribution)
.meth_mat <- matrix(
  rbeta(.n_cpgs * .n_samples, 2, 5),
  nrow     = .n_cpgs,
  ncol     = .n_samples,
  dimnames = list(.cpg_ids, .sample_ids)
)

# Binary mutation matrix (samples x genes)
.var_mat <- matrix(
  rbinom(.n_samples * .n_genes_var, 1, 0.1),
  nrow     = .n_samples,
  ncol     = .n_genes_var,
  dimnames = list(.sample_ids, paste0("Gene", seq_len(.n_genes_var)))
)

# Sample metadata
.meta_df <- data.frame(
  sample_id = .sample_ids,
  condition = rep(c("Tumour","Normal"), .n_samples / 2),
  time_os   = runif(.n_samples, 6, 60),
  event     = rbinom(.n_samples, 1, 0.5),
  row.names = .sample_ids,
  stringsAsFactors = FALSE
)
