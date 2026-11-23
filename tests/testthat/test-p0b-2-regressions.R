# P0-b (Adim 2/2) regresyon testleri: P0-5, P0-8, P0-10

test_that("P0-5: treatment labels normalize across variants and pooled note warns", {
  df <- data.frame(
    signature = rep("SigA", 4L),
    method = rep("mean", 4L),
    dataset = c("c1", "c2", "c3", "c4"),
    n = c(30L, 35L, 28L, 32L),
    auc = c(0.70, 0.72, 0.55, 0.58),
    auc_var = c(0.01, 0.01, 0.015, 0.015),
    skipped = rep(FALSE, 4L),
    independent = rep(TRUE, 4L),
    low_coverage = rep(FALSE, 4L),
    coverage_status = rep("ok", 4L),
    treatment = c("anti-PD-1", "Anti PD1", "anti-CTLA-4", "anti-CTLA4"),
    endpoint_type = c("RECIST", "RECIST", "clinical_benefit", "clinical_benefit"),
    stringsAsFactors = FALSE
  )

  sm_tx <- suppressWarnings(summarize_auc(df, method = "mean", stratify_by = "treatment"))
  expect_identical(sort(sm_tx$treatment), c("anti-CTLA4", "anti-PD1"))
  expect_identical(nrow(sm_tx), 2L)

  sm_pool <- suppressWarnings(summarize_auc(df, method = "mean", stratify_by = NULL))
  expect_true(grepl("stratify_by|pooled|heterogen", sm_pool$note[1], ignore.case = TRUE))
})

test_that("P0-8: builtin signature panel passes scientific integrity checks", {
  sigs <- get_builtin_signatures(include_controls = TRUE)

  # 1) B2M negatif kontrollerde olmamali (MHC-I / ICI direnci)
  expect_false("B2M" %in% sigs$NegCtrl_Housekeeping$genes)
  expect_false("B2M" %in% sigs$NegCtrl_Reference_geNorm$genes)

  # 2) Birebir duplike imza (Ayers_IFNG6 == IFNG6) olmamali
  expect_false("Ayers_IFNG6" %in% names(sigs))

  # 3) Chemokine_12 Messina 2012 paneli CXCL13 icermeli
  expect_true("CXCL13" %in% sigs$Chemokine_12$genes)
  expect_false("CX3CL1" %in% sigs$Chemokine_12$genes)

  # 4) Tum imzalar >= 3 gen, gecerli direction/category/citation ve cakismasiz olmali
  for (nm in names(sigs)) {
    s <- sigs[[nm]]
    expect_true(s$direction %in% c("positive", "negative", "signed"), info = nm)
    expect_true(s$category %in% c("literature_ici", "negative_control"), info = nm)
    expect_true(is.character(s$citation) && nzchar(s$citation), info = nm)
    all_g <- if (!is.null(s$genes)) s$genes else c(s$genes_up, s$genes_down)
    req_min <- if (!is.null(s$min_genes)) as.integer(s$min_genes) else 3L
    expect_gte(length(unique(toupper(all_g))), req_min)
    expect_identical(length(all_g), length(unique(toupper(all_g))))
  }
})

test_that("P0-10: survival join matches by method, blocks cartesian duplicates, and cleans negative time", {
  skip_if_not_installed("survival")

  auc_res <- data.frame(
    dataset = c("coh1", "coh1"),
    signature = c("SigA", "SigA"),
    method = c("mean", "zscore"),
    auc = c(0.65, 0.72),
    stringsAsFactors = FALSE
  )
  surv_res <- data.frame(
    cohort = "coh1",
    signature = "SigA",
    method = "zscore",
    hr = 0.55,
    cox_p = 0.02,
    logrank_p = 0.03,
    c_index = 0.68,
    stringsAsFactors = FALSE
  )

  # 1) mean satirina zscore sagkalimi kopyalanmamali
  combo <- combine_response_and_survival(auc_res, surv_res)
  expect_identical(nrow(combo), 2L)
  expect_true(is.na(combo$hr[combo$method == "mean"]))
  expect_equal(combo$hr[combo$method == "zscore"], 0.55)

  # 2) Duplicate survival satirlari kartezyen carpima yol acmamali
  surv_dup <- rbind(surv_res, surv_res)
  expect_warning(
    combo_dup <- combine_response_and_survival(auc_res, surv_dup),
    "duplicate",
    ignore.case = TRUE
  )
  expect_identical(nrow(combo_dup), nrow(auc_res))

  # 3) Negatif survival suresi validate_dataset ve ozel time_col uzerinden engellenmeli
  set.seed(1)
  n <- 20L
  mat <- matrix(rnorm(3 * n), nrow = 3, dimnames = list(c("G1", "G2", "G3"), paste0("S", 1:n)))
  ds_bad <- list(
    expr = mat,
    sample_meta = data.frame(
      sample_id = colnames(mat),
      response_label = rep(c("CR", "PD"), each = 10L),
      os_time = c(-5, seq_len(n - 1L) * 10),
      os_event = rep(c(1L, 0L), 10L),
      stringsAsFactors = FALSE
    ),
    info = list(normalization = "TPM", scale = "log2")
  )
  sc <- stats::setNames(rnorm(n), colnames(mat))
  expect_error(evaluate_survival(ds_bad, sc), "negative")

  ds_custom <- ds_bad
  ds_custom$sample_meta$os_time <- NULL
  ds_custom$sample_meta$os_event <- NULL
  ds_custom$sample_meta$custom_time <- c(-5, seq_len(n - 1L) * 10)
  ds_custom$sample_meta$custom_event <- rep(c(1L, 0L), 10L)
  expect_warning(
    ev <- evaluate_survival(ds_custom, sc, time_col = "custom_time", event_col = "custom_event"),
    "negative",
    ignore.case = TRUE
  )
  expect_identical(ev$n_samples, 19L)
})
