test_that("check_sensitivity fragile FALSE / TRUE / NA", {
  results <- data.frame(
    dataset = c("D1", "D1", "D1", "D1", "D2"),
    signature = c("S1", "S1", "S2", "S2", "S3"),
    method = c("mean", "zscore", "mean", "zscore", "mean"),
    auc = c(0.80, 0.82, 0.60, 0.85, 0.70),
    skipped = FALSE,
    stringsAsFactors = FALSE
  )

  sens <- check_sensitivity(results, threshold = 0.10)

  s1 <- sens[sens$signature == "S1", ]
  expect_equal(s1$n_methods, 2L)
  expect_equal(s1$auc_range, 0.02, tolerance = 1e-12)
  expect_false(s1$fragile)

  s2 <- sens[sens$signature == "S2", ]
  expect_equal(s2$auc_range, 0.25, tolerance = 1e-12)
  expect_true(s2$fragile)

  s3 <- sens[sens$signature == "S3", ]
  expect_equal(s3$n_methods, 1L)
  expect_true(is.na(s3$fragile))
})

test_that("low_coverage isi haritasinda auc_plot NA (gri) olur", {
  results <- data.frame(
    dataset = c("D1", "D2"),
    signature = c("S1", "S1"),
    method = c("mean", "mean"),
    auc = c(0.9, 0.8),
    coverage_pct = c(100, 40),
    coverage_label = c("5/5 (100%)", "2/5 (40%)"),
    low_coverage = c(FALSE, TRUE),
    skipped = FALSE,
    stringsAsFactors = FALSE
  )

  # plot_auc_heatmap ic mantigini dogrudan dogrula
  df <- results
  coverage_threshold <- 0.70
  df$low_coverage <- df$low_coverage | df$coverage_pct < (100 * coverage_threshold)
  df$auc_plot <- df$auc
  df$auc_plot[df$low_coverage] <- NA_real_

  expect_true(is.finite(df$auc_plot[1]))
  expect_true(is.na(df$auc_plot[2]))
  expect_true(df$low_coverage[2])

  p <- plot_auc_heatmap(results, method = "mean", coverage_threshold = 0.70)
  expect_s3_class(p, "ggplot")
  expect_match(p$labels$caption, "dusuk gen kapsama|Gri hucreler")
})

test_that("run_benchmark coverage_label ve median/rank uretir", {
  skip_if_not_installed("pROC")
  data(datasets_demo, package = "SigBenchR")
  data(signatures_demo, package = "SigBenchR")

  res <- run_benchmark(
    signatures = signatures_demo["Cytokine-5"],
    datasets = datasets_demo["demo_A"],
    methods = c("median", "rank"),
    coverage_threshold = 0.70
  )
  expect_true(all(c("median", "rank") %in% res$method))
  expect_true("coverage_label" %in% names(res))
  expect_match(res$coverage_label[1], "^[0-9]+/[0-9]+")
  expect_true(all(is.finite(res$auc[!res$skipped])))
})
