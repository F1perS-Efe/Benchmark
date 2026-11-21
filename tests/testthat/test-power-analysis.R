test_that("estimate_auc_power Hanley-McNeil large vs small n", {
  big <- estimate_auc_power(200L, 200L, target_auc = 0.70)
  expect_true(is.list(big))
  expect_true(is.finite(big$power))
  expect_true(is.finite(big$se_target))
  expect_gt(big$power, 0.99)

  small <- estimate_auc_power(10L, 10L, target_auc = 0.70)
  expect_true(is.finite(small$power))
  expect_lt(small$power, 0.50)
})

test_that("estimate_auc_power n_pos=1 returns NA and warns", {
  expect_warning(
    out <- estimate_auc_power(1L, 20L, target_auc = 0.70),
    "n_pos"
  )
  expect_true(is.na(out$power))
  expect_identical(out$power, NA_real_)
})

test_that("run_benchmark adds power_auc70 column", {
  skip_if_not_installed("pROC")
  data(datasets_demo, package = "SigBenchR")
  data(signatures_demo, package = "SigBenchR")
  res <- run_benchmark(
    signatures = signatures_demo["IFNG6"],
    datasets = datasets_demo["demo_A"],
    methods = "mean",
    compute_pr = FALSE,
    compute_null = FALSE,
    min_n = 5L,
    min_per_class = 2L
  )
  expect_true("power_auc70" %in% names(res))
  expect_true(any(is.finite(res$power_auc70)))
})
