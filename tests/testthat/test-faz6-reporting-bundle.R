# FAZ 6 — Reporting, Export & Publication Bundle Integration Tests

test_that("6.1 & 6.2: build_report_model, export_report, and export_publication_bundle preserve Phase 4-5 statistical diagnostics and full numeric precision", {
  skip_if_not_installed("pROC")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("rmarkdown")
  skip_if_not(rmarkdown::pandoc_available(), "pandoc not available")

  set.seed(601)
  inf_genes <- c("CD8A", "CD8B", "CD3D", "CD3E")
  sig_genes <- c("IFNG", "CXCL9", "CXCL10", "STAT1", "IDO1", "HLA-DRA")
  all_genes <- unique(c(inf_genes, sig_genes, paste0("BG", 1:35)))

  # Cohort 1: small positive class (n_pos = 4 => EPV = 2 < 5, small_positive_class = TRUE) + moderate positive effect + survival
  n1 <- 20
  y1 <- c(rep(1L, 4), rep(0L, 16))
  mat1 <- matrix(rnorm(length(all_genes) * n1), nrow = length(all_genes), ncol = n1,
                 dimnames = list(all_genes, paste0("C1_", seq_len(n1))))
  sig_lat1 <- 0.6 * y1 + rnorm(n1, sd = 1.0)
  for (g in sig_genes) mat1[g, ] <- sig_lat1 + rnorm(n1, sd = 0.3)
  ds1 <- list(
    expr = mat1,
    sample_meta = data.frame(
      sample_id = colnames(mat1),
      response_label = ifelse(y1 == 1L, "CR", "PD"),
      os_time = rexp(n1, rate = 0.08) + 1.0,
      os_event = sample(c(0L, 1L), n1, replace = TRUE, prob = c(0.3, 0.7)),
      stringsAsFactors = FALSE
    ),
    info = list(normalization = "TPM", scale = "log2", treatment = "anti-PD1", endpoint_type = "RECIST")
  )

  # Cohort 2: non-separating positive effect (beta_adjusted > 0) + survival
  n2 <- 50
  y2 <- c(rep(1L, 25), rep(0L, 25))
  sig_lat2 <- 0.8 * y2 + rnorm(n2, sd = 1.0)
  mat2 <- matrix(rnorm(length(all_genes) * n2), nrow = length(all_genes), ncol = n2,
                 dimnames = list(all_genes, paste0("C2_", seq_len(n2))))
  for (g in sig_genes) mat2[g, ] <- sig_lat2 + rnorm(n2, sd = 0.3)
  ds2 <- list(
    expr = mat2,
    sample_meta = data.frame(
      sample_id = colnames(mat2),
      response_label = ifelse(y2 == 1L, "CR", "PD"),
      os_time = rexp(n2, rate = 0.08) + 1.0,
      os_event = sample(c(0L, 1L), n2, replace = TRUE, prob = c(0.3, 0.7)),
      stringsAsFactors = FALSE
    ),
    info = list(normalization = "TPM", scale = "log2", treatment = "anti-PD1", endpoint_type = "RECIST")
  )

  # Cohort 3: non-separating opposite direction effect (beta_adjusted < 0 => direction_conflict = TRUE) + survival
  n3 <- 50
  y3 <- c(rep(1L, 25), rep(0L, 25))
  sig_lat3 <- -0.8 * y3 + rnorm(n3, sd = 1.0)
  mat3 <- matrix(rnorm(length(all_genes) * n3), nrow = length(all_genes), ncol = n3,
                 dimnames = list(all_genes, paste0("C3_", seq_len(n3))))
  for (g in sig_genes) mat3[g, ] <- sig_lat3 + rnorm(n3, sd = 0.3)
  ds3 <- list(
    expr = mat3,
    sample_meta = data.frame(
      sample_id = colnames(mat3),
      response_label = ifelse(y3 == 1L, "CR", "PD"),
      os_time = rexp(n3, rate = 0.08) + 1.0,
      os_event = sample(c(0L, 1L), n3, replace = TRUE, prob = c(0.3, 0.7)),
      stringsAsFactors = FALSE
    ),
    info = list(normalization = "TPM", scale = "log2", treatment = "anti-PD1", endpoint_type = "RECIST")
  )

  # Cohort 4: high collinearity with infiltration (VIF >= 5)
  n4 <- 36
  y4 <- c(rep(1L, 18), rep(0L, 18))
  lat4 <- 0.6 * y4 + rnorm(n4, sd = 1.0)
  mat4 <- matrix(rnorm(length(all_genes) * n4), nrow = length(all_genes), ncol = n4,
                 dimnames = list(all_genes, paste0("C4_", seq_len(n4))))
  for (g in c(inf_genes, sig_genes)) mat4[g, ] <- lat4 + rnorm(n4, sd = 0.25)
  ds4 <- list(
    expr = mat4,
    sample_meta = data.frame(
      sample_id = colnames(mat4),
      response_label = ifelse(y4 == 1L, "CR", "PD"),
      os_time = rexp(n4, rate = 0.08) + 1.0,
      os_event = sample(c(0L, 1L), n4, replace = TRUE, prob = c(0.3, 0.7)),
      stringsAsFactors = FALSE
    ),
    info = list(normalization = "TPM", scale = "log2", treatment = "anti-PD1", endpoint_type = "RECIST")
  )

  datasets <- list(
    cohort_pos_low_epv = ds1,
    cohort_pos_effect = ds2,
    cohort_neg_effect = ds3,
    cohort_collinear = ds4
  )
  signatures <- list(
    IFNG6 = list(genes = sig_genes, direction = "positive")
  )

  res <- suppressWarnings(run_benchmark(
    signatures = signatures,
    datasets = datasets,
    methods = c("mean", "zscore"),
    compute_pr = TRUE,
    n_boot = 25L,
    compute_null = TRUE,
    n_perm = 20L,
    seed = 42L
  ))
  sm <- suppressWarnings(summarize_auc(res, method = "zscore", stratify_by = NULL))
  sens <- check_sensitivity(res, threshold = 0.10)
  surv_res <- suppressWarnings(run_survival_benchmark(datasets, signatures, method = "zscore"))
  bi_res <- suppressWarnings(test_beyond_infiltration(datasets, signatures$IFNG6, infiltration_axis = "CD8_T_Cell", method = "zscore"))

  params_list <- list(
    methods = c("mean", "zscore"),
    direction = "<",
    sd_handling = "exclude",
    coverage_threshold = 0.70,
    sensitivity_threshold = 0.10,
    n_boot = 25L,
    n_perm = 20L,
    seed = 42L,
    compute_pr = TRUE,
    compute_null = TRUE,
    match_expression = FALSE,
    heatmap_method = "zscore",
    forest_signature = "IFNG6",
    datasets = datasets,
    signatures = signatures
  )

  rm <- build_report_model(
    results = res,
    summary_tbl = sm,
    sensitivity_tbl = sens,
    params_list = params_list,
    surv_res = surv_res,
    beyond_inf = list(IFNG6 = bi_res)
  )

  # Verify Phase 4-6 alerts are present in report_model$tier1$alerts
  alert_ids <- rm$tier1$alerts$rule_id
  expect_true(all(c(
    "PR-SMALL-CLASS", "PR-CI-UNSTABLE", "NULL-UNIVERSE",
    "BI-DIRECTION-CONFLICT", "BI-LOW-EPV", "BI-HIGH-VIF"
  ) %in% alert_ids))

  man <- build_reproducibility_manifest(
    results = res,
    summary_tbl = sm,
    report_model = rm,
    params_list = params_list,
    datasets = datasets,
    signatures = signatures
  )
  expect_true(isTRUE(man$diagnostics_summary$has_small_positive_class))
  expect_true(isTRUE(man$diagnostics_summary$has_unreliable_pr_ci))
  expect_true(isTRUE(man$diagnostics_summary$has_survival))
  expect_true(isTRUE(man$diagnostics_summary$has_beyond_infiltration))
  expect_true("subsampled_panel_cohort" %in% unlist(man$diagnostics_summary$null_background_types))

  # 1) HTML Report generation and content verification
  html_path <- tempfile("faz6_report_", fileext = ".html")
  on.exit(unlink(html_path), add = TRUE)
  export_report(
    results = res,
    summary_tbl = sm,
    sensitivity_tbl = sens,
    params_list = params_list,
    output_file = html_path,
    surv_res = surv_res,
    report_model = rm,
    manifest = man,
    beyond_inf = rm$tier2$beyond_infiltration
  )
  expect_true(file.exists(html_path))
  html_txt <- paste(readLines(html_path, warn = FALSE), collapse = "\n")

  # Verify separate raw vs adjusted p-value headers & diagnostic explanations in HTML report
  expect_true(grepl("Raw Null p", html_txt, fixed = TRUE))
  expect_true(grepl("Within-Method BH FDR", html_txt, fixed = TRUE))
  expect_true(grepl("Global BH FDR", html_txt, fixed = TRUE))
  expect_true(grepl("Null Background Universe", html_txt, fixed = TRUE))
  expect_true(grepl("Subsampled panel cohort", html_txt, fixed = TRUE))
  expect_true(grepl("PR CI Unreliable", html_txt, fixed = TRUE))
  expect_true(grepl("Beyond Infiltration (Adjusted Logistic Regression &amp; Paired CV Delta-AUC)", html_txt, fixed = TRUE) ||
              grepl("Beyond Infiltration", html_txt, fixed = TRUE))
  expect_true(grepl("Directional Stouffer Z", html_txt, fixed = TRUE))
  expect_true(grepl("Fisher Omnibus p (Unsigned Ref)", html_txt, fixed = TRUE))
  expect_true(grepl("Raw LRT p", html_txt, fixed = TRUE))
  expect_true(grepl("BH-Adjusted LRT FDR", html_txt, fixed = TRUE))
  expect_true(grepl("Raw Cox p", html_txt, fixed = TRUE))
  expect_true(grepl("Cox BH-Adjusted FDR", html_txt, fixed = TRUE))
  expect_true(grepl("Raw Log-Rank p", html_txt, fixed = TRUE))
  expect_true(grepl("Log-Rank BH-Adjusted FDR", html_txt, fixed = TRUE))
  expect_true(grepl("PH Diagnostic p (Unadjusted)", html_txt, fixed = TRUE))
  expect_true(grepl("NanoString nCounter", html_txt, fixed = TRUE))
  expect_true(grepl(man$reproducibility_hash, html_txt, fixed = TRUE))

  # 2) Publication Bundle ZIP generation and file/column/precision verification
  zip_path <- tempfile("faz6_bundle_", fileext = ".zip")
  ex_dir <- tempfile("faz6_unzip_")
  dir.create(ex_dir)
  on.exit({
    unlink(zip_path)
    unlink(ex_dir, recursive = TRUE)
  }, add = TRUE)

  export_publication_bundle(
    results = res,
    output_path = zip_path,
    summary_tbl = sm,
    report_model = rm,
    manifest = man,
    params_list = params_list,
    method = "zscore",
    signature = "IFNG6",
    zip = TRUE
  )
  expect_true(file.exists(zip_path))
  utils::unzip(zip_path, exdir = ex_dir)

  expected_files <- c(
    "fig_auc_heatmap.png", "fig_auc_heatmap.pdf",
    "fig_auc_forest.png", "fig_auc_forest.pdf",
    "fig_roc_curves.png", "fig_roc_curves.pdf",
    "table_meta_summary.csv", "table_cohort_metrics.csv",
    "table_tier1_verdicts.csv", "table_tier1_alerts.csv",
    "table_beyond_infiltration_per_cohort.csv",
    "table_beyond_infiltration_summary.csv",
    "table_survival_benchmark.csv",
    "analysis_objects.rds", "manifest.json", "figure_legends.txt"
  )
  expect_true(all(expected_files %in% list.files(ex_dir)))

  # Check CSV columns and full numeric precision (not rounded to 3 digits)
  coh_csv <- utils::read.csv(file.path(ex_dir, "table_cohort_metrics.csv"), stringsAsFactors = FALSE)
  expect_true(all(c(
    "auc", "auc_ci_low", "auc_ci_high",
    "pr_auc", "pr_ci_low", "pr_ci_high", "tie_fraction", "small_positive_class", "n_boot_valid", "ci_unreliable",
    "null_p", "null_p_adj_bh", "null_p_adj_bh_global", "background_type", "universe_size"
  ) %in% names(coh_csv)))
  expect_equal(coh_csv$auc, res$auc, tolerance = 1e-12)
  expect_equal(coh_csv$null_p_adj_bh_global, res$null_p_adj_bh_global, tolerance = 1e-12)

  bi_pc_csv <- utils::read.csv(file.path(ex_dir, "table_beyond_infiltration_per_cohort.csv"), stringsAsFactors = FALSE)
  expect_true(all(c(
    "epv", "low_epv", "vif", "high_collinearity",
    "delta_auc_combined", "sd_delta_auc", "delta_auc_ci_low", "delta_auc_ci_high",
    "lrt_p", "lrt_p_adj_bh", "analysis_status", "adds_independent_value"
  ) %in% names(bi_pc_csv)))

  bi_sm_csv <- utils::read.csv(file.path(ex_dir, "table_beyond_infiltration_summary.csv"), stringsAsFactors = FALSE)
  expect_true(all(c(
    "combined_lrt_p", "fisher_omnibus_p", "stouffer_z", "direction_conflict"
  ) %in% names(bi_sm_csv)))
  expect_true(isTRUE(bi_sm_csv$direction_conflict[1]))

  surv_csv <- utils::read.csv(file.path(ex_dir, "table_survival_benchmark.csv"), stringsAsFactors = FALSE)
  expect_true(all(c("cox_p", "p_adj_bh", "logrank_p", "logrank_p_adj_bh", "ph_p", "ph_violated") %in% names(surv_csv)))

  rds_obj <- readRDS(file.path(ex_dir, "analysis_objects.rds"))
  expect_identical(rds_obj$manifest$reproducibility_hash, man$reproducibility_hash)
  expect_equal(rds_obj$results$auc, res$auc, tolerance = 1e-15)
})
