test_that("jaccard_signatures Cytokine-5 / IFNG6 ortusmesi", {
  data(signatures_demo, package = "SigBenchR")
  jac <- jaccard_signatures(signatures_demo)
  expect_true(nrow(jac) >= 3)
  row <- jac[jac$signature1 == "Cytokine-5" & jac$signature2 == "IFNG6", ]
  if (nrow(row) == 0L) {
    row <- jac[jac$signature1 == "IFNG6" & jac$signature2 == "Cytokine-5", ]
  }
  expect_equal(nrow(row), 1)
  # IFNG ve CXCL9 ortak
  expect_equal(row$n_intersect, 2L)
  expect_gt(row$jaccard, 0)
  expect_lt(row$jaccard, 1)
})

test_that("summarize_auc n/IVW, var_zero ve independent kurallari", {
  # Oyuncak: iki bagimsiz + bir derived + bir var_zero
  res <- data.frame(
    dataset = c("A", "B", "C", "D"),
    signature = "S1",
    method = "mean",
    n = c(20, 40, 30, 25),
    auc = c(0.7, 0.8, 0.9, 1.0),
    auc_var = c(0.01, 0.04, 0.01, 0),
    var_zero = c(FALSE, FALSE, FALSE, TRUE),
    independent = c(TRUE, TRUE, FALSE, TRUE),
    small_cohort = c(FALSE, FALSE, FALSE, FALSE),
    skipped = FALSE,
    stringsAsFactors = FALSE
  )

  expect_warning(
    sm <- summarize_auc(res, method = "mean", include_nonindependent = FALSE),
    "zero or NA variance|varyansi 0"
  )
  expect_equal(nrow(sm), 1)
  # C non-indep dislandi; D var_zero n-agirlikliya dahil, IVW'den cikarildi
  expect_equal(sm$n_datasets, 3L)
  expect_equal(sm$n_excluded_nonindependent, 1L)
  expect_equal(sm$n_excluded_var_zero, 1L)
  expect_equal(sm$n_var_zero_excluded, 1L)
  expect_equal(sm$n_cohorts_ivw, 2L)
  expect_equal(sm$auc_n_weighted, (20 * 0.7 + 40 * 0.8 + 25 * 1.0) / 85)
  w <- c(1 / 0.01, 1 / 0.04)
  expect_equal(sm$auc_ivw, sum(w * c(0.7, 0.8)) / sum(w))
  expect_true(is.finite(sm$Q))
  expect_true(is.finite(sm$I2))

  sm2 <- suppressWarnings(
    summarize_auc(res, include_nonindependent = TRUE, exclude_var_zero = FALSE)
  )
  expect_equal(sm2$n_datasets, 4L)
})

test_that("auc_correlation ve adjust_null_p BH", {
  skip_if_not_installed("pROC")
  data(datasets_demo, package = "SigBenchR")
  data(signatures_demo, package = "SigBenchR")

  res <- run_benchmark(
    signatures = signatures_demo[c("Cytokine-5", "IFNG6", "Random3")],
    datasets = datasets_demo,
    methods = "mean",
    compute_pr = FALSE,
    compute_null = TRUE,
    n_perm = 15,
    seed = 3
  )

  expect_true("null_p_adj_bh" %in% names(res))
  finite_p <- is.finite(res$null_p)
  if (any(finite_p)) {
    expect_equal(
      res$null_p_adj_bh[finite_p],
      p.adjust(res$null_p[finite_p], method = "BH")
    )
  }

  corr <- auc_correlation(res, method = "mean")
  expect_true(nrow(corr) >= 3)
  expect_true(all(c("signature1", "signature2", "correlation", "n_datasets") %in% names(corr)))
})

test_that("plot_auc_forest ggplot doner", {
  skip_if_not_installed("pROC")
  data(datasets_demo, package = "SigBenchR")
  data(signatures_demo, package = "SigBenchR")
  res <- run_benchmark(
    signatures = signatures_demo["Cytokine-5"],
    datasets = datasets_demo,
    methods = "mean",
    compute_pr = FALSE,
    compute_null = FALSE
  )
  sm <- suppressWarnings(summarize_auc(res, method = "mean"))
  p <- plot_auc_forest(res, signature = "Cytokine-5", method = "mean", summary = sm)
  expect_s3_class(p, "ggplot")
})
