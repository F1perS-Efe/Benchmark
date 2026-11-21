test_that("generate_r_code gecerli R sozdizimi ve parametreler uretir", {
  code <- generate_r_code(
    signature_name = "Cytokine-5",
    methods = c("mean", "zscore"),
    direction = "<",
    coverage_threshold = 0.70,
    n_boot = 100L,
    n_perm = 25L,
    seed = 99L,
    dataset_names = c("demo_A", "demo_B"),
    compute_pr = TRUE,
    compute_null = FALSE
  )

  expect_type(code, "character")
  expect_length(code, 1L)
  expect_silent(parse(text = code))
  expect_true(grepl("seed = 99L", code, fixed = TRUE))
  expect_true(grepl("n_boot = 100L", code, fixed = TRUE))
  expect_true(grepl("n_perm = 25L", code, fixed = TRUE))
  expect_true(grepl("Cytokine-5", code, fixed = TRUE))
  expect_true(grepl('"mean"', code, fixed = TRUE))
  expect_true(grepl('"zscore"', code, fixed = TRUE))
  expect_true(grepl("coverage_threshold = 0.7", code, fixed = TRUE))

  custom <- generate_r_code(
    signature_name = NULL,
    genes = c("IFNG", "CXCL9"),
    methods = "mean",
    seed = 7L
  )
  expect_silent(parse(text = custom))
  expect_true(grepl("IFNG", custom, fixed = TRUE))
  expect_true(grepl("seed = 7L", custom, fixed = TRUE))
})

test_that("export_report gecici dizinde HTML uretir", {
  skip_if_not_installed("rmarkdown")
  skip_if_not_installed("knitr")
  skip_if_not_installed("pROC")
  skip_if_not(rmarkdown::pandoc_available(), "pandoc yok")

  data(datasets_demo, package = "SigBenchR")
  data(signatures_demo, package = "SigBenchR")

  res <- run_benchmark(
    signatures = signatures_demo["Cytokine-5"],
    datasets = datasets_demo,
    methods = c("mean", "zscore"),
    compute_pr = FALSE,
    compute_null = FALSE,
    seed = 42L
  )
  sm <- suppressWarnings(summarize_auc(res, method = "mean"))
  sens <- check_sensitivity(res, threshold = 0.10)

  out <- file.path(tempdir(), "sigbenchr_phase3a_report.html")
  if (file.exists(out)) file.remove(out)

  path <- export_report(
    results = res,
    summary_tbl = sm,
    sensitivity_tbl = sens,
    params_list = list(
      methods = c("mean", "zscore"),
      direction = "<",
      coverage_threshold = 0.70,
      sensitivity_threshold = 0.10,
      n_boot = 0L,
      n_perm = 0L,
      seed = 42L,
      heatmap_method = "mean",
      forest_signature = "Cytokine-5"
    ),
    output_file = out
  )

  expect_true(file.exists(path))
  html <- paste(readLines(path, warn = FALSE), collapse = "\n")
  expect_true(grepl("klinik karar", html, ignore.case = TRUE))
  expect_true(grepl("sessionInfo", html, fixed = TRUE) || grepl("R version", html))
})

test_that("export_report rmarkdown yoksa anlasilir hata verir", {
  skip_if_not_installed("pROC")
  # requireNamespace mock edilemezse: en azindan imza / arguman kontrolleri
  expect_error(
    export_report(
      results = data.frame(),
      summary_tbl = data.frame(),
      sensitivity_tbl = data.frame(),
      params_list = list(),
      output_file = ""
    ),
    "output_file"
  )
})
