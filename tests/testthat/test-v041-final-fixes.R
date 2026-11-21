# v0.4.1 final audit regression tests

test_that("plot_roc_curves AUC matches run_benchmark under SD exclusion (zscore/rank)", {
  skip_if_not_installed("pROC")
  skip_if_not_installed("ggplot2")
  data(datasets_ici, package = "SigBenchR")
  data(signatures_demo, package = "SigBenchR")

  ds <- datasets_ici["riaz_2017"]
  sigs <- signatures_demo[c("Cytokine-5", "IFNG6")]

  for (meth in c("zscore", "rank")) {
    res <- run_benchmark(
      signatures = sigs,
      datasets = ds,
      methods = meth,
      direction = "<",
      sd_handling = "exclude",
      compute_pr = FALSE,
      compute_null = FALSE,
      min_n = 5L,
      min_per_class = 2L
    )
    p <- plot_roc_curves(
      datasets = ds,
      signatures = sigs,
      dataset_name = "riaz_2017",
      method = meth,
      direction = "<",
      sd_handling = "exclude"
    )
    plot_auc <- unique(p$data[, c("signature", "auc")])
    for (i in seq_len(nrow(plot_auc))) {
      sig <- plot_auc$signature[i]
      expect_equal(
        plot_auc$auc[i],
        res$auc[res$signature == sig & res$method == meth][1],
        tolerance = 1e-12,
        info = paste(meth, sig)
      )
    }
  }
})

test_that("generate_r_code reproduces sd_handling, match_expression, n_perm=1000", {
  expect_identical(formals(generate_r_code)$n_perm, 1000L)
  expect_identical(formals(generate_r_code)$sd_handling, "exclude")
  expect_identical(formals(generate_r_code)$match_expression, FALSE)

  code <- generate_r_code(
    signature_name = "IFNG6",
    methods = c("mean", "zscore"),
    direction = "<",
    n_perm = 1000L,
    dataset_names = c("riaz_2017", "gide_2019"),
    compute_pr = FALSE,
    compute_null = TRUE,
    sd_handling = "as_0",
    match_expression = TRUE
  )
  expect_error(parse(text = code), NA)
  expect_true(grepl("sd_handling\\s*=\\s*[\"']as_0[\"']", code))
  expect_true(grepl("match_expression\\s*=\\s*TRUE", code))
  expect_true(grepl("n_perm\\s*=\\s*1000L", code))
})

test_that("down-only signed null keeps negative score sign and finite null_p", {
  skip_if_not_installed("pROC")
  set.seed(41)
  n <- 24L
  n_pool <- 40L
  genes_down <- paste0("D", 1:4)
  genes_missing_up <- paste0("MISSING_UP", 1:3)
  pool <- paste0("P", seq_len(n_pool))
  all_genes <- c(genes_down, pool)
  mat <- matrix(
    rnorm(length(all_genes) * n),
    nrow = length(all_genes),
    dimnames = list(all_genes, paste0("S", seq_len(n)))
  )
  # Down arm informative: higher expression in non-responders -> negative score favors responders
  mat[genes_down, 1:12] <- mat[genes_down, 1:12] - 1.5
  mat[genes_down, 13:24] <- mat[genes_down, 13:24] + 1.5
  y <- c(rep(1L, 12L), rep(0L, 12L))
  names(y) <- colnames(mat)

  ds <- list(
    synth = list(
      expr = mat,
      sample_meta = data.frame(
        sample_id = colnames(mat),
        response_label = c(rep("CR", 12L), rep("PD", 12L)),
        stringsAsFactors = FALSE
      ),
      info = list(normalization = "TPM", scale = "log2")
    )
  )
  sigs <- list(
    DownOnly = list(
      genes_up = genes_missing_up,
      genes_down = genes_down,
      direction = "<",
      derived_from = character(0)
    )
  )

  expect_warning(
    res <- run_benchmark(
      signatures = sigs,
      datasets = ds,
      methods = "mean",
      compute_pr = FALSE,
      compute_null = TRUE,
      n_perm = 40L,
      seed = 9L,
      match_expression = TRUE,
      min_genes = 3L,
      min_n = 5L,
      min_per_class = 2L
    ),
    "Signed signature"
  )
  expect_equal(as.integer(res$n_up_found[1]), 0L)
  expect_gte(as.integer(res$n_down_found[1]), 3L)
  expect_true(is.finite(res$null_p[1]))
  expect_gte(res$null_p[1], 0)
  expect_lte(res$null_p[1], 1)
  expect_true(is.finite(res$auc[1]))

  scored <- suppressWarnings(
    score_signature(
      mat,
      genes_up = genes_missing_up,
      genes_down = genes_down,
      method = "mean",
      min_genes = 3L
    )
  )
  scores_dn <- colMeans(mat[genes_down, , drop = FALSE], na.rm = TRUE)
  expect_equal(unname(scored$scores), unname(-scores_dn), tolerance = 1e-12)

  null_out <- null_permutation_test(
    mat = mat,
    response = y,
    n_genes = 0L,
    n_genes_down = 3L,
    exclude_genes = genes_down,
    method = "mean",
    direction = "<",
    n_perm = 40L,
    seed = 9L,
    match_expression = TRUE,
    signature_genes = NULL,
    signature_genes_down = genes_down[1:3],
    observed_scores = scored$scores,
    min_genes = 3L
  )
  expect_true(is.finite(null_out$obs_auc))
  expect_true(is.finite(null_out$null_p))
  expect_gte(null_out$null_p, 0)
  expect_lte(null_out$null_p, 1)
  expect_gt(null_out$obs_auc, 0.5)
})

test_that("forest uses unstratified REML summary even when stratify_by is open", {
  skip_if_not_installed("ggplot2")
  res <- data.frame(
    dataset = c("A", "B", "C", "D"),
    signature = "S1",
    method = "mean",
    n = c(40L, 40L, 40L, 40L),
    auc = c(0.62, 0.70, 0.74, 0.80),
    auc_ci_low = c(0.45, 0.52, 0.55, 0.60),
    auc_ci_high = c(0.78, 0.85, 0.90, 0.95),
    auc_var = c(0.01, 0.01, 0.012, 0.011),
    var_zero = FALSE,
    independent = TRUE,
    small_cohort = FALSE,
    skipped = FALSE,
    low_coverage = FALSE,
    treatment = c("anti-PD1", "anti-PD1", "anti-CTLA4", "anti-CTLA4"),
    endpoint_type = c("RECIST", "RECIST", "clinical_benefit", "clinical_benefit"),
    stringsAsFactors = FALSE
  )
  sm_strat <- suppressWarnings(
    summarize_auc(res, method = "mean", stratify_by = c("treatment", "endpoint_type"))
  )
  expect_gt(nrow(sm_strat), 1L)

  sm_flat <- suppressWarnings(
    summarize_auc(res, method = "mean", stratify_by = NULL)
  )
  expect_equal(nrow(sm_flat), 1L)
  expect_true(is.finite(sm_flat$auc_re[1]))
  expect_true(is.finite(sm_flat$auc_re_ci_low[1]))

  p <- plot_auc_forest(res, signature = "S1", method = "mean", summary = sm_flat)
  expect_s3_class(p, "ggplot")
  expect_true(any(as.character(p$data$lab) == "REML + HKSJ (95 pct)"))

  # Stratified summary alone would drop REML row (nrow != 1 after filter)
  p_bad <- plot_auc_forest(res, signature = "S1", method = "mean", summary = sm_strat)
  expect_false(any(as.character(p_bad$data$lab) == "REML + HKSJ (95 pct)"))

  rmd <- system.file("rmd", "benchmark_report.Rmd", package = "SigBenchR")
  if (!nzchar(rmd) || !file.exists(rmd)) {
    rmd <- file.path("inst", "rmd", "benchmark_report.Rmd")
  }
  rmd_txt <- paste(readLines(rmd, warn = FALSE), collapse = "\n")
  expect_true(grepl("summarize_auc\\(results, method = hm, stratify_by = NULL\\)", rmd_txt))
})

test_that("shiny app parses and guards input$run before results()", {
  app_path <- system.file("shiny", "app.R", package = "SigBenchR")
  if (!nzchar(app_path) || !file.exists(app_path)) {
    app_path <- file.path("inst", "shiny", "app.R")
  }
  expect_true(file.exists(app_path))
  expect_error(parse(file = app_path), NA)

  src <- paste(readLines(app_path, warn = FALSE), collapse = "\n")
  expect_true(grepl("v0\\.4\\.1", src))
  # kpi_strip: run==0 check appears before results() in that block
  kpi_block <- sub(".*output\\$kpi_strip\\s*<-\\s*shiny::renderUI\\(\\{", "", src)
  kpi_block <- sub("\\}\\)\\s*\\n\\s*sens\\s*<-.*", "", kpi_block)
  run_pos <- regexpr("identical\\(as\\.integer\\(input\\$run\\),\\s*0L\\)", kpi_block)
  res_pos <- regexpr("results\\(\\)", kpi_block)
  expect_true(run_pos[1] > 0)
  expect_true(res_pos[1] > 0)
  expect_lt(as.integer(run_pos[1]), as.integer(res_pos[1]))
})
