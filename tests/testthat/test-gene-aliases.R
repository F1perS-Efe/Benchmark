test_that("match_genes alias rescue: PDL1/PD1 -> CD274/PDCD1", {
  expr <- matrix(
    1:6,
    nrow = 3L,
    dimnames = list(c("CD274", "PDCD1", "GAPDH"), c("S1", "S2"))
  )

  cov_on <- match_genes(
    expr,
    c("PDL1", "PD1"),
    use_aliases = TRUE
  )
  expect_equal(cov_on$n_found, 2L)
  expect_equal(cov_on$n_total, 2L)
  expect_equal(cov_on$coverage_pct, 100)
  expect_equal(cov_on$n_alias_rescued, 2L)
  expect_true(all(c("PDL1->CD274", "PD1->PDCD1") %in% cov_on$alias_rescued))
  expect_equal(length(cov_on$missing), 0L)

  cov_off <- match_genes(
    expr,
    c("PDL1", "PD1"),
    use_aliases = FALSE
  )
  expect_equal(cov_off$n_found, 0L)
  expect_equal(cov_off$n_alias_rescued, 0L)
  expect_equal(sort(cov_off$missing), c("PD1", "PDL1"))
})

test_that("run_benchmark reports n_alias_rescued", {
  skip_if_not_installed("pROC")
  set.seed(1)
  genes <- c("CD274", "PDCD1", "GAPDH", "ACTB")
  n <- 12L
  expr <- matrix(
    rnorm(length(genes) * n),
    nrow = length(genes),
    dimnames = list(genes, paste0("S", seq_len(n)))
  )
  ds <- list(
    toy = list(
      expr = expr,
      sample_meta = data.frame(
        sample_id = colnames(expr),
        response_label = rep(c("CR", "PD"), each = 6L),
        stringsAsFactors = FALSE
      ),
      info = list(normalization = "TPM", scale = "log2")
    )
  )
  sigs <- list(
    AliasSig = list(
      genes = c("PDL1", "PD1", "GAPDH"),
      direction = "<",
      derived_from = character(0)
    )
  )
  res <- run_benchmark(
    signatures = sigs,
    datasets = ds,
    methods = "mean",
    compute_pr = FALSE,
    compute_null = FALSE,
    min_genes = 2L,
    min_n = 5L,
    min_per_class = 2L
  )
  expect_true("n_alias_rescued" %in% names(res))
  expect_equal(as.integer(res$n_alias_rescued[1]), 2L)
})
