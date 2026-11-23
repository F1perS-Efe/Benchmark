# P1-b (Adim 2/2) testleri: HTML rapor 3 katmanli cikti ve Shiny Tier-1 kartlari

test_that("P1-b-2: export_report renders 3-tier HTML with alerts and rules MD5", {
  skip_if_not_installed("rmarkdown")
  skip_if_not_installed("knitr")
  skip_if_not_installed("pROC")
  skip_if_not(rmarkdown::pandoc_available(), "pandoc yok")

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
  sm <- suppressWarnings(summarize_auc(res, method = "mean"))
  sens <- check_sensitivity(res, threshold = 0.10)
  bio_ctx <- correlate_reference_axes(ds, sigs, method = "mean")

  out <- file.path(tempdir(), "sigbenchr_p1b2_report.html")
  if (file.exists(out)) file.remove(out)

  path <- export_report(
    results = res,
    summary_tbl = sm,
    sensitivity_tbl = sens,
    params_list = list(
      methods = "mean",
      direction = "<",
      coverage_threshold = 0.70,
      sensitivity_threshold = 0.10,
      n_boot = 0L,
      n_perm = 0L,
      seed = 42L,
      compute_null = FALSE,
      heatmap_method = "mean",
      forest_signature = "IFNG6",
      bio_context = bio_ctx
    ),
    output_file = out
  )

  expect_true(file.exists(path))
  html <- paste(readLines(path, warn = FALSE), collapse = "\n")

  expect_true(grepl("Tier 1: Executive Summary", html, fixed = TRUE))
  expect_true(grepl("Tier 2: Comparative Evidence", html, fixed = TRUE))
  expect_true(grepl("Tier 3: Technical Appendix", html, fixed = TRUE))
  expect_true(grepl("R_REF_OVERLAP", html, fixed = TRUE))
  expect_true(grepl("R_SUBGROUP_POOLED", html, fixed = TRUE))

  rules <- load_interpretation_rules()
  expect_true(grepl(rules$md5, html, fixed = TRUE))
  expect_true(grepl("@media print", html, fixed = TRUE) || grepl("beforeprint", html, fixed = TRUE))
})

test_that("P1-b-2: Shiny Tier-1 cards and alerts render from snapshot without live-input drift", {
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
      datasets = c("van_allen_2015", "hugo_2016", "gide_2019"),
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

    m_before <- report_model()
    kpi_html_before <- as.character(output$kpi_strip$html)

    expect_identical(m_before$tier1$best_signature$signature, "IFNG6")
    expect_true(grepl("IFNG6", kpi_html_before, fixed = TRUE))
    expect_true(grepl("Pooled REML AUC", kpi_html_before, fixed = TRUE))
    expect_true(grepl("R_REF_OVERLAP", kpi_html_before, fixed = TRUE))

    # Canli girdileri run tetiklemeden degistir; snapshot kilidi bozulmamali
    session$setInputs(sd_handling = "as_0", direction = ">")
    m_after <- report_model()
    expect_identical(m_after$tier1$best_signature$auc_re, m_before$tier1$best_signature$auc_re)
  })
})
