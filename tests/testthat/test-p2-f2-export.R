# P2-F2 testleri: Yayina hazir disa aktarma paketi (.zip & klasor)

test_that("P2-F2: export_publication_bundle creates non-empty zip with figures, CSVs, manifest, and legends", {
  skip_if_not_installed("pROC")
  skip_if_not_installed("ggplot2")
  data(datasets_ici, package = "SigBenchR")
  data(signatures_demo, package = "SigBenchR")

  sigs <- signatures_demo[c("Cytokine-5", "IFNG6")]
  ds <- datasets_ici[c("van_allen_2015", "hugo_2016", "gide_2019")]
  res <- suppressWarnings(run_benchmark(
    signatures = sigs,
    datasets = ds,
    methods = "mean",
    compute_pr = FALSE,
    compute_null = FALSE,
    seed = 42L
  ))
  pl <- list(
    methods = "mean",
    direction = "<",
    sd_handling = "exclude",
    coverage_threshold = 0.70,
    seed = 42L,
    heatmap_method = "mean",
    forest_signature = "IFNG6"
  )

  zip_out <- file.path(tempdir(), "sigbenchr_pub_bundle_test.zip")
  if (file.exists(zip_out)) file.remove(zip_out)

  path <- export_publication_bundle(
    results = res,
    output_path = zip_out,
    params_list = pl,
    method = "mean",
    signature = "IFNG6",
    coverage_threshold = 0.70,
    zip = TRUE
  )

  expect_true(file.exists(path))
  listing <- utils::unzip(path, list = TRUE)
  expected_files <- c(
    "fig_auc_heatmap.png",
    "fig_auc_heatmap.pdf",
    "fig_auc_forest.png",
    "fig_auc_forest.pdf",
    "table_meta_summary.csv",
    "table_cohort_metrics.csv",
    "table_tier1_verdicts.csv",
    "manifest.json",
    "figure_legends.txt"
  )
  for (fn in expected_files) {
    expect_true(fn %in% listing$Name, info = fn)
    expect_gt(listing$Length[listing$Name == fn][1], 0, label = fn)
  }

  # zip = FALSE klasor modu da calismali
  dir_out <- file.path(tempdir(), "sigbenchr_pub_dir_test")
  if (dir.exists(dir_out)) unlink(dir_out, recursive = TRUE)
  export_publication_bundle(
    results = res,
    output_path = dir_out,
    params_list = pl,
    method = "mean",
    signature = "IFNG6",
    zip = FALSE
  )
  expect_true(file.exists(file.path(dir_out, "figure_legends.txt")))
  legends <- paste(readLines(file.path(dir_out, "figure_legends.txt"), warn = FALSE), collapse = "\n")
  expect_true(grepl("Reproducibility hash", legends, fixed = TRUE))
  expect_true(grepl("REML", legends, fixed = TRUE))
})

test_that("P2-F2: Shiny download_pub_bundle handler exists and exports from snapshot", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("pROC")

  app_path <- system.file("shiny", "app.R", package = "SigBenchR")
  if (!nzchar(app_path)) {
    app_path <- file.path("inst", "shiny", "app.R")
  }
  env <- new.env(parent = globalenv())
  source(app_path, local = env)

  shiny::testServer(env$server, {
    session$setInputs(
      sig_choice = c("Cytokine-5", "IFNG6"),
      custom_genes = "",
      datasets = c("hugo_2016", "gide_2019"),
      methods = "mean",
      direction = "<",
      sd_handling = "exclude",
      coverage_threshold = 0.70,
      sensitivity_threshold = 0.10,
      compute_pr = FALSE,
      compute_null = FALSE,
      n_boot = 0L,
      n_perm = 10L,
      match_expression = FALSE,
      seed = 42L,
      heatmap_method = "mean",
      forest_sig = "IFNG6"
    )
    session$setInputs(run = 1L)

    snap <- results()
    tmp_zip <- tempfile(fileext = ".zip")
    SigBenchR::export_publication_bundle(
      results = snap$res,
      output_path = tmp_zip,
      report_model = report_model(),
      manifest = reproducibility_manifest(),
      params_list = snap$params,
      method = snap$params$heatmap_method,
      signature = snap$params$forest_signature,
      coverage_threshold = snap$params$coverage_threshold,
      zip = TRUE
    )
    expect_true(file.exists(tmp_zip))
    expect_gt(file.info(tmp_zip)$size, 0)
  })
})
