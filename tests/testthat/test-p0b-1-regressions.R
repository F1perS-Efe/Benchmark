# P0-b (Adim 1/2) regresyon testleri: P0-4, P0-2, P0-7

test_that("P0-4: single coverage threshold is consistent across benchmark, heatmap, and summary", {
  skip_if_not_installed("pROC")
  data(datasets_ici, package = "SigBenchR")

  # 6 genli imzanin 3 geni hugo_2016'da yok -> coverage_pct = 50%
  sigs <- list(
    HalfCov = list(
      genes = c("IFNG", "CXCL9", "CXCL10", "FAKE_G1", "FAKE_G2", "FAKE_G3"),
      direction = "positive",
      derived_from = character(0)
    )
  )
  ds <- datasets_ici["hugo_2016"]

  expect_true(exists(".classify_coverage", where = asNamespace("SigBenchR"), inherits = FALSE))

  res <- suppressWarnings(run_benchmark(
    signatures = sigs,
    datasets = ds,
    methods = "mean",
    compute_pr = FALSE,
    compute_null = FALSE
  ))

  # %50 kapsama varsayilan 0.70 esigi altinda -> low_coverage = TRUE VE coverage_status = "low_power" olmali
  expect_true(isTRUE(res$low_coverage[1]))
  expect_identical(res$coverage_status[1], "low_power")

  # Isı haritasinda gri olan hucre summarize_auc(exclude_low_coverage = TRUE) tarafindan da dislanmali
  sm <- suppressWarnings(summarize_auc(res, method = "mean", exclude_low_coverage = TRUE))
  expect_identical(sm$n_cohorts[1], 0L)
})

test_that("P0-2: Custom signature validation prevents name collision, self-DeLong, and supports independent summary toggle", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("pROC")

  app_path <- system.file("shiny", "app.R", package = "SigBenchR")
  if (!nzchar(app_path)) {
    app_path <- file.path("inst", "shiny", "app.R")
  }
  env <- new.env(parent = globalenv())
  source(app_path, local = env)

  # 1) Yerlesik imza adiyla cakisan custom imza sessizce uzerine yazmamali
  expect_true(exists(".validate_custom_signature", envir = env, inherits = FALSE))
  expect_error(
    env$.validate_custom_signature(
      genes_text = "IFNG, CXCL9, STAT1",
      custom_name = "IFNG6",
      existing_names = c("Cytokine-5", "IFNG6")
    ),
    "already exists|collision|builtin",
    ignore.case = TRUE
  )

  # 2) Bos / gecersiz gen listesi veya up/down cakismasi yakalanmali
  expect_error(
    env$.validate_custom_signature(
      genes_text = "   ",
      custom_name = "Custom",
      existing_names = "IFNG6"
    )
  )
  expect_error(
    env$.validate_custom_signature(
      genes_up = c("IFNG", "CXCL9"),
      genes_down = c("IFNG", "VEGFA"),
      custom_name = "Custom",
      existing_names = "IFNG6"
    ),
    "overlap",
    ignore.case = TRUE
  )

  # 3) Shiny server icinde Custom ozet toggle ve ayni gen setiyle kendiyle DeLong engeli
  shiny::testServer(env$server, {
    session$setInputs(
      sig_choice = "IFNG6",
      custom_genes = "IDO1, CXCL9, CXCL10, HLA-DRA, STAT1, IFNG",
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
      forest_sig = "IFNG6",
      stratify_by = FALSE,
      include_nonindependent = TRUE
    )
    session$setInputs(run = 1L)

    sm <- summary_res()
    custom_row <- sm[sm$signature == "Custom", , drop = FALSE]
    expect_identical(nrow(custom_row), 1L)
    expect_identical(custom_row$n_cohorts[1], 1L)

    # Kendiyle (birebir ayni gen listesi) DeLong karsilastirmasi engellenmeli
    session$setInputs(
      cmp_dataset = "hugo_2016",
      cmp_sig1 = "IFNG6",
      cmp_sig2 = "Custom"
    )
    session$setInputs(run_cmp = 1L)
    expect_error(cmp_res(), "identical|same gene", ignore.case = TRUE)
  })
})

test_that("P0-7: correlate_reference_axes flags circular overlap with reference axes", {
  data(datasets_ici, package = "SigBenchR")
  b <- get_builtin_signatures(include_controls = FALSE)

  bc <- correlate_reference_axes(
    datasets = datasets_ici["hugo_2016"],
    signatures = b["IFNG6"],
    method = "mean"
  )

  expect_true("jaccard_overlap" %in% names(bc))
  expect_true("overlap_flag" %in% names(bc))
  expect_true("overlap_note" %in% names(bc))

  ifn_row <- bc[bc$reference_axis == "IFN_Gamma", , drop = FALSE]
  other_rows <- bc[bc$reference_axis != "IFN_Gamma", , drop = FALSE]

  expect_equal(ifn_row$jaccard_overlap[1], 1.0, tolerance = 1e-12)
  expect_true(isTRUE(ifn_row$overlap_flag[1]))
  expect_true(!is.na(ifn_row$overlap_note[1]) && nzchar(ifn_row$overlap_note[1]))

  expect_true(all(!other_rows$overlap_flag))
})
