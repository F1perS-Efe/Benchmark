test_that("demo veri ile run_benchmark AUC üretir", {
  skip_if_not_installed("pROC")
  data(datasets_demo, package = "SigBenchR")
  data(signatures_demo, package = "SigBenchR")

  res <- run_benchmark(
    signatures = signatures_demo[c("Cytokine-5", "Random3")],
    datasets = datasets_demo["demo_A"],
    methods = "mean"
  )

  expect_true(nrow(res) >= 2)
  cyto <- res[res$signature == "Cytokine-5" & !res$skipped, ]
  rand <- res[res$signature == "Random3" & !res$skipped, ]
  expect_true(nrow(cyto) == 1)
  expect_true(is.finite(cyto$auc))
  # Sinyalli imza rastgele kontrolden daha yüksek AUC bekler
  if (nrow(rand) == 1 && is.finite(rand$auc)) {
    expect_gt(cyto$auc, rand$auc)
  }
})

test_that("plot_auc_heatmap ggplot döndürür", {
  skip_if_not_installed("pROC")
  data(datasets_demo, package = "SigBenchR")
  data(signatures_demo, package = "SigBenchR")
  res <- run_benchmark(
    signatures = signatures_demo["Cytokine-5"],
    datasets = datasets_demo,
    methods = "mean"
  )
  p <- plot_auc_heatmap(res, method = "mean")
  expect_s3_class(p, "ggplot")
})
