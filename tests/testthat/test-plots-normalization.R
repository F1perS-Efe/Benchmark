test_that("plot_auc_heatmap and plot_auc_forest add [TPM]/[FPKM] labels", {
  skip_if_not_installed("ggplot2")

  res <- data.frame(
    dataset = c("hugo_2016", "riaz_2017", "gide_2019"),
    signature = "S1",
    method = "mean",
    n = c(20L, 30L, 25L),
    auc = c(0.70, 0.75, 0.68),
    auc_ci_low = c(0.50, 0.55, 0.48),
    auc_ci_high = c(0.90, 0.95, 0.88),
    coverage_pct = c(100, 100, 100),
    coverage_label = c("5/5 (100%)", "5/5 (100%)", "5/5 (100%)"),
    low_coverage = FALSE,
    independent = TRUE,
    var_zero = FALSE,
    small_cohort = FALSE,
    skipped = FALSE,
    normalization = c("TPM", "FPKM", "TPM"),
    stringsAsFactors = FALSE
  )

  hm <- plot_auc_heatmap(res, method = "mean", show_normalization = TRUE)
  expect_s3_class(hm, "ggplot")
  expect_true("dataset_label" %in% names(hm$data))
  labs_hm <- as.character(hm$data$dataset_label)
  expect_true(any(grepl("\\[TPM\\]", labs_hm)))
  expect_true(any(grepl("\\[FPKM\\]", labs_hm)))
  expect_match(
    hm$labels$caption,
    "heterogeneous expression scales",
    fixed = TRUE
  )

  fr <- plot_auc_forest(res, signature = "S1", method = "mean", show_normalization = TRUE)
  expect_s3_class(fr, "ggplot")
  labs_fr <- as.character(fr$data$lab)
  expect_true(any(grepl("\\[TPM\\]", labs_fr)))
  expect_true(any(grepl("\\[FPKM\\]", labs_fr)))
  expect_match(
    fr$labels$caption,
    "heterogeneous expression scales",
    fixed = TRUE
  )
})
