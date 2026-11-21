.make_cov_ds <- function(gene_pool, n = 12L) {
  set.seed(42)
  expr <- matrix(
    rnorm(length(gene_pool) * n, mean = 5),
    nrow = length(gene_pool),
    dimnames = list(gene_pool, paste0("S", seq_len(n)))
  )
  list(
    cohortX = list(
      expr = expr,
      sample_meta = data.frame(
        sample_id = colnames(expr),
        response_label = rep(c("CR", "PD"), each = n / 2L),
        stringsAsFactors = FALSE
      ),
      info = list(normalization = "TPM", scale = "log2")
    )
  )
}

test_that("coverage_status invalid / low_power / ok and summarize exclusion", {
  skip_if_not_installed("pROC")

  # 10 genes in matrix; signatures choose how many overlap
  pool <- paste0("G", 1:10)
  ds <- .make_cov_ds(pool)

  sigs <- list(
    InvalidSig = list(
      genes = c(pool[1], paste0("MISS", 1:9)),
      direction = "<",
      derived_from = character(0)
    ),
    LowPowerSig = list(
      # 3/10 = 0.30 -> low_power (0.2 <= x < 0.5)
      genes = c(pool[1:3], paste0("MISS", 1:7)),
      direction = "<",
      derived_from = character(0)
    ),
    OkSig = list(
      genes = pool[1:6],
      direction = "<",
      derived_from = character(0)
    )
  )

  expect_warning(
    res <- run_benchmark(
      signatures = sigs,
      datasets = ds,
      methods = "mean",
      compute_pr = FALSE,
      compute_null = FALSE,
      min_genes = 1L,
      min_n = 5L,
      min_per_class = 2L,
      min_coverage = 0.50,
      invalid_coverage = 0.20,
      coverage_threshold = 0.70
    ),
    "Invalid gene coverage"
  )

  inv <- res[res$signature == "InvalidSig", ]
  low <- res[res$signature == "LowPowerSig", ]
  ok <- res[res$signature == "OkSig", ]

  expect_equal(inv$coverage_status[1], "invalid")
  expect_lt(inv$coverage_pct[1] / 100, 0.20)
  expect_equal(low$coverage_status[1], "low_power")
  expect_gte(low$coverage_pct[1] / 100, 0.20)
  expect_lt(low$coverage_pct[1] / 100, 0.50)
  expect_equal(ok$coverage_status[1], "ok")
  expect_gte(ok$coverage_pct[1] / 100, 0.50)

  # Build a small summary frame with AUC so summarize_auc can pool
  res2 <- res
  res2$skipped <- FALSE
  res2$auc <- c(0.7, 0.75, 0.8)
  res2$n <- 12L
  res2$auc_var <- 0.01
  res2$var_zero <- FALSE
  res2$independent <- TRUE
  res2$small_cohort <- FALSE

  sm <- summarize_auc(res2, exclude_low_coverage = TRUE)
  expect_true("n_excluded_coverage" %in% names(sm))
  # One row per signature in this toy; each signature is its own group
  sm_ok <- sm[sm$signature == "OkSig", ]
  expect_equal(sm_ok$n_datasets[1], 1L)
  expect_equal(sm_ok$n_excluded_coverage[1], 0L)

  sm_low <- sm[sm$signature == "LowPowerSig", ]
  expect_equal(sm_low$n_datasets[1], 0L)
  expect_equal(sm_low$n_excluded_coverage[1], 1L)

  sm_inv <- sm[sm$signature == "InvalidSig", ]
  expect_equal(sm_inv$n_datasets[1], 0L)
  expect_equal(sm_inv$n_excluded_coverage[1], 1L)
})
