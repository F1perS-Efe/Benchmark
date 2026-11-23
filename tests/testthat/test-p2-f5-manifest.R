# P2-F5 testleri: Reproducibility Manifesti

test_that("P2-F5: build_reproducibility_manifest produces deterministic hash across runs and detects changes", {
  skip_if_not_installed("pROC")
  data(datasets_ici, package = "SigBenchR")
  data(signatures_demo, package = "SigBenchR")

  sigs <- signatures_demo[c("Cytokine-5", "IFNG6")]
  ds <- datasets_ici["hugo_2016"]
  params1 <- list(methods = "mean", direction = "<", sd_handling = "exclude", coverage_threshold = 0.70, seed = 42L)

  res1 <- suppressWarnings(run_benchmark(sigs, ds, methods = "mean", compute_pr = FALSE, compute_null = FALSE, seed = 42L))
  m1 <- build_reproducibility_manifest(res1, params_list = params1, datasets = ds, signatures = sigs)

  Sys.sleep(1.1) # timestamp_utc farkli olsun
  res2 <- suppressWarnings(run_benchmark(sigs, ds, methods = "mean", compute_pr = FALSE, compute_null = FALSE, seed = 42L))
  m2 <- build_reproducibility_manifest(res2, params_list = params1, datasets = ds, signatures = sigs)

  expect_false(identical(m1$timestamp_utc, m2$timestamp_utc))
  expect_identical(m1$reproducibility_hash, m2$reproducibility_hash)
  expect_true(nzchar(m1$sigbench_version))
  expect_true(nzchar(m1$rules_meta$md5))
  expect_true("hugo_2016" %in% names(m1$cohorts_fingerprint))
  expect_true("IFNG6" %in% names(m1$signatures_fingerprint))

  # Tohum degisince hash degismeli
  params_seed <- params1; params_seed$seed <- 99L
  m_seed <- build_reproducibility_manifest(res1, params_list = params_seed, datasets = ds, signatures = sigs)
  expect_false(identical(m1$reproducibility_hash, m_seed$reproducibility_hash))

  # Esik degisince hash degismeli
  params_cov <- params1; params_cov$coverage_threshold <- 0.80
  m_cov <- build_reproducibility_manifest(res1, params_list = params_cov, datasets = ds, signatures = sigs)
  expect_false(identical(m1$reproducibility_hash, m_cov$reproducibility_hash))

  # Gen listesi degisince hash degismeli
  sigs_mod <- sigs; sigs_mod$IFNG6$genes <- c(sigs_mod$IFNG6$genes, "CXCL11")
  m_gene <- build_reproducibility_manifest(res1, params_list = params1, datasets = ds, signatures = sigs_mod)
  expect_false(identical(m1$reproducibility_hash, m_gene$reproducibility_hash))
})

test_that("P2-F5: export_manifest round-trips JSON and YAML and embeds hash in HTML report", {
  skip_if_not_installed("pROC")
  data(datasets_ici, package = "SigBenchR")
  data(signatures_demo, package = "SigBenchR")

  sigs <- signatures_demo["IFNG6"]
  ds <- datasets_ici["hugo_2016"]
  res <- suppressWarnings(run_benchmark(sigs, ds, methods = "mean", compute_pr = FALSE, compute_null = FALSE, seed = 42L))
  pl <- list(
    methods = "mean",
    direction = "<",
    sd_handling = "exclude",
    coverage_threshold = 0.70,
    sensitivity_threshold = 0.10,
    n_boot = 0L,
    n_perm = 0L,
    seed = 42L,
    compute_null = FALSE,
    heatmap_method = "mean",
    forest_signature = "IFNG6"
  )
  man <- build_reproducibility_manifest(res, params_list = pl, datasets = ds, signatures = sigs)

  json_path <- file.path(tempdir(), "sigbenchr_manifest_test.json")
  yaml_path <- file.path(tempdir(), "sigbenchr_manifest_test.yml")
  export_manifest(man, json_path, format = "json")
  export_manifest(man, yaml_path, format = "yaml")

  expect_true(file.exists(json_path))
  expect_true(file.exists(yaml_path))

  j_read <- jsonlite::fromJSON(json_path, simplifyVector = FALSE)
  y_read <- yaml::read_yaml(yaml_path)
  expect_identical(j_read$reproducibility_hash, man$reproducibility_hash)
  expect_identical(y_read$reproducibility_hash, man$reproducibility_hash)
  expect_identical(j_read$rules_meta$md5, man$rules_meta$md5)

  skip_if_not(rmarkdown::pandoc_available(), "pandoc yok")
  html_path <- file.path(tempdir(), "sigbenchr_p2f5_report.html")
  if (file.exists(html_path)) file.remove(html_path)
  export_report(
    results = res,
    summary_tbl = suppressWarnings(summarize_auc(res, method = "mean")),
    sensitivity_tbl = check_sensitivity(res),
    params_list = c(pl, list(datasets = ds, signatures = sigs)),
    output_file = html_path
  )
  html <- paste(readLines(html_path, warn = FALSE), collapse = "\n")
  expect_true(grepl(man$reproducibility_hash, html, fixed = TRUE))
})

test_that("P2-F5: Shiny builds reproducibility_manifest from snapshot", {
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
      datasets = "hugo_2016",
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
    m1 <- reproducibility_manifest()
    expect_true(nzchar(m1$reproducibility_hash))

    # Canli girdiyi run olmadan degistir; manifest snapshot hash'i sabit kalmali
    session$setInputs(seed = 999L, coverage_threshold = 0.90)
    m2 <- reproducibility_manifest()
    expect_identical(m1$reproducibility_hash, m2$reproducibility_hash)
  })
})
