# P1-b (Adim 1/2) testleri: YAML kurallari ve build_report_model() 3 katmanli yapi

test_that("load_interpretation_rules loads bundled YAML with md5 provenance", {
  rules <- load_interpretation_rules()
  expect_type(rules, "list")
  expect_true(nzchar(rules$version))
  expect_true(nzchar(rules$md5))
  expect_identical(rules$coverage$ok_min_pct, 70L)
  expect_identical(rules$heterogeneity$high_i2_logit, 50L)
})

test_that("build_report_model constructs 3-tier structure and fires deterministic alerts", {
  results <- data.frame(
    dataset = c("c1", "c2", "c3"),
    signature = c("SigA", "SigA", "SigA"),
    method = c("zscore", "zscore", "zscore"),
    n = c(30L, 40L, 35L),
    n_pos = c(15L, 20L, 18L),
    n_neg = c(15L, 20L, 17L),
    n_found = c(5L, 6L, 2L),
    n_total = c(6L, 6L, 6L),
    coverage_pct = c(83.3, 100.0, 33.3),
    coverage_status = c("ok", "ok", "low_power"),
    low_coverage = c(FALSE, FALSE, TRUE),
    n_alias_rescued = c(0L, 0L, 0L),
    missing_genes = c("G6", "", "G3,G4,G5,G6"),
    auc = c(0.75, 0.78, 0.55),
    auc_var = c(0.01, 0.01, 0.02),
    skipped = c(FALSE, FALSE, FALSE),
    independent = c(TRUE, TRUE, TRUE),
    treatment = c("anti-PD1", "anti-CTLA4", "anti-PD1"),
    endpoint_type = c("RECIST", "clinical_benefit", "RECIST"),
    stringsAsFactors = FALSE
  )

  summary_tbl <- data.frame(
    signature = "SigA",
    method = "zscore",
    n_cohorts = 2L,
    auc_n_weighted = 0.767,
    auc_ivw = 0.765,
    auc_re = 0.765,
    auc_re_ci_low = 0.480,
    auc_re_ci_high = 0.920,
    I2 = 62.0,
    i2_logit = 58.4,
    tau2_reml = 0.140,
    note = "1 low-coverage cohorts excluded; Pooled across distinct subgroups (treatment=[anti-PD1, anti-CTLA4]); stratify_by recommended",
    stringsAsFactors = FALSE
  )

  ref_axes <- data.frame(
    cohort = "c1",
    signature = "SigA",
    reference_axis = "IFN_Gamma",
    cor = 0.98,
    p_value = 1e-12,
    n_samples = 30L,
    normalization = "TPM",
    jaccard_overlap = 0.83,
    overlap_flag = TRUE,
    overlap_note = "High gene overlap with IFN_Gamma",
    stringsAsFactors = FALSE
  )

  model <- build_report_model(
    results = results,
    summary_tbl = summary_tbl,
    sensitivity_tbl = data.frame(dataset = "c1", signature = "SigA", auc_range = 0.02),
    params_list = list(heatmap_method = "zscore", seed = 42L),
    ref_axes = ref_axes
  )

  # 3 katmanli yapi kontrolu
  expect_true(all(c("tier1", "tier2", "tier3") %in% names(model)))
  expect_identical(model$tier1$active_method, "zscore")
  expect_identical(model$tier1$best_signature$signature, "SigA")
  expect_identical(nrow(model$tier1$signature_cards), 1L)
  expect_true(grepl("Medium Priority", model$tier1$signature_cards$priority_tier[1]))

  # Tetiklenen deterministik uyarilar
  alert_ids <- model$tier1$alerts$rule_id
  expect_true("R_AUC_CI_CHANCE" %in% alert_ids)
  expect_true("R_HET_HIGH" %in% alert_ids)
  expect_true("R_COV_LOW" %in% alert_ids)
  expect_true("R_REF_OVERLAP" %in% alert_ids)
  expect_true("R_SUBGROUP_POOLED" %in% alert_ids)

  # Tier 2 ve Tier 3 icerigi
  expect_true("I2 logit (%)" %in% names(model$tier2$summary_formatted))
  expect_true(is.data.frame(model$tier2$treatment_strata))
  expect_true(is.data.frame(model$tier2$endpoint_strata))
  expect_identical(nrow(model$tier3$coverage_diagnostics), 3L)
  expect_true(nzchar(model$tier3$rules_meta$md5))
})

test_that("build_report_model runs end-to-end on bundled ICI cohorts", {
  skip_if_not_installed("pROC")
  data(datasets_ici, package = "SigBenchR")
  data(signatures_demo, package = "SigBenchR")

  res <- suppressWarnings(run_benchmark(
    signatures_demo[c("Cytokine-5", "IFNG6")],
    datasets_ici[c("hugo_2016", "gide_2019")],
    methods = "mean",
    compute_pr = FALSE,
    compute_null = FALSE
  ))
  model <- build_report_model(res, params_list = list(heatmap_method = "mean"))
  expect_identical(model$tier1$best_signature$signature, "IFNG6")
  expect_true("Cytokine-5" %in% model$tier1$signature_cards$signature)
})
