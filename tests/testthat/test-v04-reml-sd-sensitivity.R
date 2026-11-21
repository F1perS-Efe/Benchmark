test_that("homogeneous cohorts: REML tau2 ~ 0 and exact logit back-transform", {
  res <- data.frame(
    dataset = c("A", "B", "C"),
    signature = "S1",
    method = "mean",
    n = c(40L, 40L, 40L),
    auc = c(0.75, 0.75, 0.75),
    auc_var = c(0.01, 0.01, 0.01),
    var_zero = FALSE,
    independent = TRUE,
    small_cohort = FALSE,
    skipped = FALSE,
    low_coverage = FALSE,
    stringsAsFactors = FALSE
  )
  sm <- suppressWarnings(summarize_auc(res, method = "mean"))
  expect_equal(nrow(sm), 1L)
  expect_lt(sm$tau2_reml[1], 1e-6)
  expect_equal(sm$auc_re[1], 0.75, tolerance = 1e-6)
  expect_equal(sm$i2_logit[1], 0)
  expect_true(sm$auc_re_ci_low[1] < 0.75)
  expect_true(0.75 < sm$auc_re_ci_high[1])
})

test_that("heterogeneous cohorts: tau2 > 0 and HKSJ not narrower than normal RE", {
  auc <- c(0.55, 0.70, 0.82, 0.90)
  auc_var <- rep(0.004, 4L)
  res <- data.frame(
    dataset = paste0("D", seq_along(auc)),
    signature = "S1",
    method = "mean",
    n = rep(50L, 4L),
    auc = auc,
    auc_var = auc_var,
    var_zero = FALSE,
    independent = TRUE,
    small_cohort = FALSE,
    skipped = FALSE,
    low_coverage = FALSE,
    stringsAsFactors = FALSE
  )
  sm <- suppressWarnings(summarize_auc(res, method = "mean"))
  expect_gt(sm$tau2_reml[1], 0)
  expect_gt(sm$i2_logit[1], 0)
  expect_gt(sm$auc_re_ci_low[1], 0)
  expect_lt(sm$auc_re_ci_high[1], 1)

  y <- stats::qlogis(auc)
  v <- auc_var / (auc * (1 - auc))^2
  tau2 <- sm$tau2_reml[1]
  w_re <- 1 / (v + tau2)
  theta_re <- sum(w_re * y) / sum(w_re)
  k <- length(y)
  q_hk <- sum(w_re * (y - theta_re)^2) / (k - 1)
  se_hk <- sqrt(max(1, q_hk) / sum(w_re))
  se_re <- sqrt(1 / sum(w_re))
  t_crit <- stats::qt(0.975, df = k - 1)
  z_crit <- stats::qnorm(0.975)
  width_hk <- 2 * t_crit * se_hk
  width_re <- 2 * z_crit * se_re
  expect_gte(width_hk, width_re)
})

test_that("edge cases: k = 1 and AUC = 1.0 excluded from REML safely", {
  res1 <- data.frame(
    dataset = "A",
    signature = "S1",
    method = "mean",
    n = 30L,
    auc = 0.80,
    auc_var = 0.02,
    var_zero = FALSE,
    independent = TRUE,
    small_cohort = FALSE,
    skipped = FALSE,
    low_coverage = FALSE,
    stringsAsFactors = FALSE
  )
  sm1 <- suppressWarnings(summarize_auc(res1, method = "mean"))
  expect_equal(sm1$auc_re[1], 0.80, tolerance = 1e-10)
  expect_true(is.na(sm1$auc_re_ci_low[1]))
  expect_true(is.na(sm1$auc_re_ci_high[1]))
  expect_true(is.na(sm1$tau2_reml[1]))

  res_vz <- data.frame(
    dataset = c("A", "B"),
    signature = "S1",
    method = "mean",
    n = c(20L, 30L),
    auc = c(1.0, 0.7),
    auc_var = c(0, 0.04),
    var_zero = c(TRUE, FALSE),
    independent = TRUE,
    small_cohort = FALSE,
    skipped = FALSE,
    low_coverage = FALSE,
    stringsAsFactors = FALSE
  )
  expect_warning(
    sm_vz <- summarize_auc(res_vz, method = "mean"),
    "varyansi 0"
  )
  expect_equal(sm_vz$auc_n_weighted[1], (20 * 1.0 + 30 * 0.7) / 50)
  expect_equal(sm_vz$auc_re[1], 0.7, tolerance = 1e-10)
  expect_true(is.na(sm_vz$auc_re_ci_low[1]))
  expect_true(is.na(sm_vz$tau2_reml[1]))
  expect_true(is.finite(sm_vz$auc_re[1]))
})

test_that("check_sd_sensitivity on riaz_2017 and hugo_2016", {
  skip_if_not_installed("pROC")
  data(datasets_ici, package = "SigBenchR")
  data(signatures_demo, package = "SigBenchR")

  tab <- suppressWarnings(check_sd_sensitivity(
    signatures = signatures_demo["IFNG6"],
    datasets = datasets_ici[c("riaz_2017", "hugo_2016")],
    methods = "mean",
    seed = 42L
  ))

  hugo <- tab[tab$dataset == "hugo_2016", , drop = FALSE]
  expect_equal(nrow(hugo), 1L)
  expect_equal(hugo$delta_orr[1], 0, tolerance = 1e-12)
  expect_equal(hugo$delta_dcr[1], 0, tolerance = 1e-12)
  expect_false(isTRUE(hugo$sd_sensitive[1]))

  riaz <- tab[tab$dataset == "riaz_2017", , drop = FALSE]
  expect_equal(nrow(riaz), 1L)
  expect_identical(as.integer(riaz$n_exclude[1]), 33L)
  expect_identical(as.integer(riaz$n_as_0[1]), 49L)
  expect_identical(as.integer(riaz$n_as_1[1]), 49L)
  expect_true(is.finite(riaz$auc_exclude[1]))
  expect_true(is.finite(riaz$auc_as_0[1]))
  expect_true(is.finite(riaz$auc_as_1[1]))
})

test_that("plot_auc_forest accepts REML summary; shiny app parses", {
  skip_if_not_installed("pROC")
  skip_if_not_installed("ggplot2")

  res <- data.frame(
    dataset = c("A", "B", "C"),
    signature = "S1",
    method = "mean",
    n = c(40L, 40L, 40L),
    auc = c(0.65, 0.72, 0.78),
    auc_ci_low = c(0.5, 0.55, 0.6),
    auc_ci_high = c(0.8, 0.85, 0.9),
    auc_var = c(0.01, 0.01, 0.01),
    var_zero = FALSE,
    independent = TRUE,
    small_cohort = FALSE,
    skipped = FALSE,
    low_coverage = FALSE,
    stringsAsFactors = FALSE
  )
  sm <- suppressWarnings(summarize_auc(res, method = "mean"))
  expect_true(is.finite(sm$auc_re[1]))
  expect_true(is.finite(sm$auc_re_ci_low[1]))
  p <- plot_auc_forest(res, signature = "S1", method = "mean", summary = sm)
  expect_s3_class(p, "ggplot")

  app_path <- system.file("shiny", "app.R", package = "SigBenchR")
  if (!nzchar(app_path)) {
    app_path <- file.path("inst", "shiny", "app.R")
  }
  expect_true(file.exists(app_path))
  expect_error(parse(file = app_path), NA)
})
