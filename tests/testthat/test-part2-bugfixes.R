# Part 2 Regression Tests: BUG-12 through BUG-17 + BUG-04 UI/Rmd Labels

test_that("BUG-12: .pick_best_signature excludes negative controls, ranks by tier > ci_low > k > auc_re, and uses Exploratory Signal label", {
  rules <- load_interpretation_rules()
  expect_identical(rules$evidence_tiers$moderate$label, "Exploratory Signal")

  # Negative control has higher point AUC (0.78 on k = 1), whereas real signature
  # IFNG6 has k = 4, AUC = 0.72, CI = [0.56, 0.85]
  summary_tbl <- data.frame(
    signature = c("NegCtrl_Housekeeping_5", "IFNG6", "SingleCohortHigh"),
    method = c("zscore", "zscore", "zscore"),
    n_cohorts = c(3L, 4L, 1L),
    auc_re = c(0.78, 0.72, 0.76),
    auc_n_weighted = c(0.78, 0.72, 0.76),
    auc_re_ci_low = c(0.58, 0.56, 0.42),
    auc_re_ci_high = c(0.90, 0.85, 0.92),
    i2_logit = c(10, 12, NA_real_),
    stringsAsFactors = FALSE
  )

  cards <- .build_signature_cards(summary_tbl, rules)
  best <- .pick_best_signature(cards)
  expect_false(is.null(best))
  expect_identical(best$signature, "IFNG6")

  # Check CI <= 0.50 exploratory label
  single_card <- cards[cards$signature == "SingleCohortHigh", , drop = FALSE]
  expect_false(grepl("Medium Priority", single_card$priority_tier[1]))
  expect_true(grepl("Exploratory Signal", single_card$priority_tier[1]))
})

test_that("BUG-13: Weighted continuous signatures disallow method = 'median' and T_Cell_Inflamed_GEP scoring_compatibility omits 'median'", {
  sigs <- get_builtin_signatures(include_controls = FALSE)
  gep <- sigs[["T_Cell_Inflamed_GEP"]]
  expect_false("median" %in% gep$scoring_compatibility)
  expect_identical(gep$scoring_compatibility, c("mean", "zscore", "rank", "ss_rank"))

  set.seed(213)
  mat <- matrix(
    rnorm(length(gep$genes) * 12, mean = 5, sd = 1.2),
    nrow = length(gep$genes),
    ncol = 12,
    dimnames = list(gep$genes, paste0("S", 1:12))
  )
  expect_error(
    score_signature(mat, genes = gep$genes, weights = gep$weights, method = "median"),
    "Weighted continuous signatures do not support method = 'median'"
  )

  # In run_benchmark(), T_Cell_Inflamed_GEP + median marks only that cell skipped = TRUE
  ds <- list(
    cohort1 = list(
      expr = mat,
      sample_meta = data.frame(
        sample_id = colnames(mat),
        response_label = rep(c("CR", "PD"), each = 6),
        stringsAsFactors = FALSE
      ),
      info = list(normalization = "TPM", scale = "log2")
    )
  )
  rb <- run_benchmark(
    signatures = list(T_Cell_Inflamed_GEP = gep),
    datasets = ds,
    methods = c("mean", "median"),
    compute_pr = FALSE,
    compute_null = FALSE
  )
  med_row <- rb[rb$method == "median", , drop = FALSE]
  mean_row <- rb[rb$method == "mean", , drop = FALSE]
  expect_true(isTRUE(med_row$skipped[1]))
  expect_true(grepl("incompatible|median", med_row$skip_reason[1], ignore.case = TRUE))
  expect_false(isTRUE(mean_row$skipped[1]))
  expect_true(is.finite(mean_row$auc[1]))
})

test_that("BUG-14 & BUG-04 UI/Rmd: generate_r_code does not override demo signature directions and UI/Rmd use Repeat SD / Range (min-max)", {
  code1 <- generate_r_code(
    signature_name = "IFNG6",
    methods = c("mean", "zscore"),
    direction = ">",
    dataset_names = c("hugo_2016", "gide_2019")
  )
  expect_false(grepl("\\$direction\\s*<-", code1))

  code2 <- generate_r_code(
    signature_name = c("Cytokine-5", "IFNG6"),
    methods = c("mean", "zscore"),
    direction = ">",
    dataset_names = c("hugo_2016", "gide_2019")
  )
  expect_false(grepl("\\$direction\\s*<-", code2))

  # Custom signature still receives direction
  code_custom <- generate_r_code(
    signature_name = "IFNG6",
    genes = c("CD8A", "CD8B", "GZMA"),
    methods = "mean",
    direction = ">",
    dataset_names = "hugo_2016"
  )
  expect_true(grepl("sigs\\$Custom <- list\\(genes = .*direction = \">\"\\)", code_custom))

  app_path <- testthat::test_path("..", "..", "inst", "shiny", "app.R")
  if (!file.exists(app_path)) {
    app_path <- system.file("shiny", "app.R", package = "SigBenchR")
  }
  app_txt <- paste(readLines(app_path, warn = FALSE), collapse = "\n")
  expect_true(grepl("Custom signature direction", app_txt, fixed = TRUE))
  expect_true(grepl("Repeat SD / Range (min-max)", app_txt, fixed = TRUE))

  rmd_path <- testthat::test_path("..", "..", "inst", "rmd", "benchmark_report.Rmd")
  if (!file.exists(rmd_path)) {
    rmd_path <- system.file("rmd", "benchmark_report.Rmd", package = "SigBenchR")
  }
  rmd_txt <- paste(readLines(rmd_path, warn = FALSE), collapse = "\n")
  expect_true(grepl("Repeat SD / Range (min-max)", rmd_txt, fixed = TRUE))
})

test_that("BUG-15 & BUG-16: .annotate_immune_axis classifies CD68 as Pan_Macrophage and CD276/VTCN1/HHLA2 as B7_Immune_Exclusion", {
  expect_identical(.annotate_immune_axis("CD68"), "Pan_Macrophage")
  expect_identical(.annotate_immune_axis("CD276"), "B7_Immune_Exclusion")
  expect_identical(.annotate_immune_axis("VTCN1"), "B7_Immune_Exclusion")
  expect_identical(.annotate_immune_axis("HHLA2"), "B7_Immune_Exclusion")
  expect_identical(.annotate_immune_axis("PDCD1"), "Immune_Checkpoint")
  expect_identical(.annotate_immune_axis("ARG1"), "Myeloid_Suppression")

  vig_path <- testthat::test_path("..", "..", "vignettes", "SigBenchR-methodology.Rmd")
  if (file.exists(vig_path)) {
    vig_txt <- paste(readLines(vig_path, warn = FALSE), collapse = "\n")
    expect_true(grepl('infiltration_axis = "CD8_T_Cell"', vig_txt, fixed = TRUE))
  }
})

test_that("BUG-17: auc_correlation warns and flags low_k_unreliable when n_datasets < 5", {
  df_small <- data.frame(
    dataset = rep(c("D1", "D2", "D3", "D4"), times = 2),
    signature = rep(c("S1", "S2"), each = 4),
    method = "zscore",
    auc = c(0.61, 0.65, 0.70, 0.74, 0.60, 0.66, 0.69, 0.75),
    skipped = FALSE,
    stringsAsFactors = FALSE
  )
  expect_warning(
    res_small <- auc_correlation(df_small, method = "zscore"),
    "low_k_unreliable"
  )
  expect_true("low_k_unreliable" %in% names(res_small))
  expect_true(all(res_small$low_k_unreliable))

  df_ok <- data.frame(
    dataset = rep(paste0("D", 1:5), times = 2),
    signature = rep(c("S1", "S2"), each = 5),
    method = "zscore",
    auc = c(0.61, 0.65, 0.70, 0.74, 0.78, 0.60, 0.66, 0.69, 0.75, 0.77),
    skipped = FALSE,
    stringsAsFactors = FALSE
  )
  expect_silent(res_ok <- auc_correlation(df_ok, method = "zscore"))
  expect_false(any(res_ok$low_k_unreliable))
})
