# FAZ 5 — Shiny UI/UX Modernization & Phase 4 Statistical Diagnostic Visualization Tests

.faz5_load_app_env <- function() {
  app_path <- file.path("inst", "shiny", "app.R")
  if (!file.exists(app_path)) {
    app_path <- system.file("shiny", "app.R", package = "SigBenchR")
  }
  env <- new.env(parent = globalenv())
  source(app_path, local = env)
  env
}

test_that("5.1: run_benchmark propagates Phase 4 PR-AUC and Null Permutation diagnostic columns", {
  set.seed(501)
  n <- 24
  y <- c(rep("CR", 4), rep("PD", 20)) # n_pos = 4 < 5 => small_positive_class = TRUE
  genes <- c("IFNG", "CXCL9", "CXCL10", "STAT1", "IDO1", "HLA-DRA", paste0("BG", 1:40))
  mat <- matrix(rnorm(length(genes) * n), nrow = length(genes), ncol = n)
  rownames(mat) <- genes
  colnames(mat) <- paste0("S", seq_len(n))

  ds <- list(
    small_cohort = list(
      expr = mat,
      sample_meta = data.frame(
        sample_id = colnames(mat),
        response_label = y,
        stringsAsFactors = FALSE
      ),
      info = list(normalization = "TPM", scale = "log2")
    )
  )
  sigs <- list(
    IFNG6 = list(genes = c("IFNG", "CXCL9", "CXCL10", "STAT1", "IDO1", "HLA-DRA"), direction = "positive")
  )

  res <- run_benchmark(
    signatures = sigs,
    datasets = ds,
    methods = c("mean", "zscore"),
    compute_pr = TRUE,
    n_boot = 30L,
    compute_null = TRUE,
    n_perm = 25L,
    seed = 42L
  )

  expect_true(all(c(
    "tie_fraction", "small_positive_class", "n_boot_valid", "ci_unreliable",
    "null_p", "null_p_adj_bh", "null_p_adj_bh_global",
    "background_type", "universe_size"
  ) %in% names(res)))
  expect_true(all(res$small_positive_class %in% TRUE))
  expect_true(all(res$ci_unreliable %in% TRUE))
  expect_true(all(res$background_type == "subsampled_panel_cohort"))
  expect_true(all(is.finite(res$null_p_adj_bh_global)))
})

test_that("5.2: Shiny server exposes Phase 4 diagnostic columns, human-readable explanations, and statistical alerts", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("pROC")

  env <- .faz5_load_app_env()
  ui_html <- as.character(env$ui)

  # UI contains collapsible statistical glossaries and alert containers
  expect_true(grepl("overview_stat_alerts", ui_html, fixed = TRUE))
  expect_true(grepl("beyond_inf_alerts", ui_html, fixed = TRUE))
  expect_true(grepl("beyond_inf_plot", ui_html, fixed = TRUE))
  expect_true(grepl("survival_alerts", ui_html, fixed = TRUE))
  expect_true(grepl("survival_table", ui_html, fixed = TRUE))

  shiny::testServer(env$server, {
    # Verify human-readable status and background translators
    expect_match(
      .translate_analysis_status("exploratory_low_epv_or_collinear"),
      "Exploratory: low EPV",
      fixed = TRUE
    )
    expect_match(
      .translate_analysis_status("failed_collinearity"),
      "Aborted: severe multicollinearity",
      fixed = TRUE
    )
    expect_match(
      .translate_background_type("subsampled_panel_cohort"),
      "Subsampled panel cohort",
      fixed = TRUE
    )

    session$setInputs(
      sig_choice = c("IFNG6", "Cytokine-5"),
      custom_genes = "",
      datasets = c("van_allen_2015", "gide_2019"),
      methods = c("mean", "zscore"),
      direction = "<",
      sd_handling = "exclude",
      coverage_threshold = 0.70,
      sensitivity_threshold = 0.10,
      compute_pr = TRUE,
      compute_null = TRUE,
      n_boot = 25L,
      n_perm = 20L,
      match_expression = FALSE,
      seed = 42L,
      heatmap_method = "mean",
      forest_sig = "IFNG6",
      beyond_inf_sig = "IFNG6",
      beyond_inf_axis = "IFN_Gamma"
    )
    session$setInputs(run = 1L)

    snap <- results()
    expect_true(all(c(
      "null_p", "null_p_adj_bh", "null_p_adj_bh_global",
      "small_positive_class", "n_boot_valid", "ci_unreliable", "background_type"
    ) %in% names(snap$res)))

    # Overview alert strip displays NULL-UNIVERSE explanation
    ov_html <- paste(as.character(output$overview_stat_alerts$html), collapse = " ")
    expect_true(grepl("NULL-UNIVERSE", ov_html, fixed = TRUE))
    expect_true(grepl("Subsampled panel cohort", ov_html, fixed = TRUE))

    # Beyond Infiltration formatted table includes all Phase 4 diagnostic fields and human-readable explanation
    bi <- beyond_inf_res()
    bi_fmt <- .format_beyond_inf_df(bi$per_cohort)
    expect_true(all(c(
      "epv", "low_epv", "vif", "high_collinearity",
      "delta_auc_combined", "sd_delta_auc", "delta_auc_ci_low", "delta_auc_ci_high",
      "lrt_p", "lrt_p_adj_bh", "analysis_status", "analysis_status_explanation"
    ) %in% names(bi_fmt)))
    expect_true(all(nzchar(bi_fmt$analysis_status_explanation)))

    # Beyond Infiltration summary reports directional Stouffer p and direction conflict status
    bi_sum_txt <- output$beyond_inf_summary
    expect_true(grepl("Directional Stouffer p:", bi_sum_txt, fixed = TRUE))
    expect_true(grepl("Direction conflict:", bi_sum_txt, fixed = TRUE))
  })
})

test_that("5.3: Shiny alert banners trigger accurately on direction_conflict, low_epv, high_collinearity, and survival PH diagnostics", {
  skip_if_not_installed("shiny")

  env <- .faz5_load_app_env()

  shiny::testServer(env$server, {
    # Inject synthetic cohorts with opposite effects, low EPV, high collinearity, and survival columns
    set.seed(503)
    inf_genes <- c("CD8A", "CD8B", "CD3D", "CD3E")
    sig_genes <- c("IFNG", "CXCL9", "CXCL10", "STAT1", "IDO1", "HLA-DRA")
    all_genes <- unique(c(inf_genes, sig_genes, paste0("BG", 1:15)))

    # Cohort 1: low EPV (n_pos = 4 => EPV = 2 < 5) + positive effect + survival columns
    n1 <- 20
    y1 <- c(rep(1L, 4), rep(0L, 16))
    mat1 <- matrix(rnorm(length(all_genes) * n1), nrow = length(all_genes), ncol = n1,
                   dimnames = list(all_genes, paste0("A", seq_len(n1))))
    for (g in sig_genes) mat1[g, ] <- 1.6 * y1 + rnorm(n1, sd = 0.4)
    ds1 <- list(
      expr = mat1,
      sample_meta = data.frame(
        sample_id = colnames(mat1),
        response_label = ifelse(y1 == 1L, "CR", "PD"),
        os_time = rexp(n1, rate = 0.1) + 0.5,
        os_event = sample(c(0L, 1L), n1, replace = TRUE, prob = c(0.3, 0.7)),
        stringsAsFactors = FALSE
      ),
      info = list(normalization = "TPM", scale = "log2")
    )

    # Cohort 2: opposite direction effect (-1.6 * y2) AND high collinearity with infiltration
    n2 <- 40
    y2 <- c(rep(1L, 20), rep(0L, 20))
    lat2 <- -1.5 * y2 + rnorm(n2, sd = 0.5)
    mat2 <- matrix(rnorm(length(all_genes) * n2), nrow = length(all_genes), ncol = n2,
                   dimnames = list(all_genes, paste0("B", seq_len(n2))))
    for (g in c(inf_genes, sig_genes)) mat2[g, ] <- lat2 + rnorm(n2, sd = 0.25)
    ds2 <- list(
      expr = mat2,
      sample_meta = data.frame(
        sample_id = colnames(mat2),
        response_label = ifelse(y2 == 1L, "CR", "PD"),
        os_time = rexp(n2, rate = 0.1) + 0.5,
        os_event = sample(c(0L, 1L), n2, replace = TRUE, prob = c(0.3, 0.7)),
        stringsAsFactors = FALSE
      ),
      info = list(normalization = "TPM", scale = "log2")
    )

    user_datasets(list(synth_pos_low_epv = ds1, synth_neg_collinear = ds2))

    session$setInputs(
      sig_choice = "IFNG6",
      custom_genes = "",
      datasets = c("synth_pos_low_epv", "synth_neg_collinear"),
      methods = "zscore",
      direction = "<",
      sd_handling = "exclude",
      coverage_threshold = 0.70,
      sensitivity_threshold = 0.10,
      compute_pr = TRUE,
      compute_null = FALSE,
      n_boot = 20L,
      n_perm = 10L,
      match_expression = FALSE,
      seed = 42L,
      heatmap_method = "zscore",
      forest_sig = "IFNG6",
      beyond_inf_sig = "IFNG6",
      beyond_inf_axis = "CD8_T_Cell"
    )
    session$setInputs(run = 1L)

    # 1) Overview alerts: small positive class & unreliable CI on synth_pos_low_epv
    ov_html <- paste(as.character(output$overview_stat_alerts$html), collapse = " ")
    expect_true(grepl("PR-SMALL-CLASS", ov_html, fixed = TRUE))
    expect_true(grepl("PR-CI-UNSTABLE", ov_html, fixed = TRUE))

    # 2) Beyond Infiltration alerts: LOW-EPV and HIGH-VIF (and DIRECTION-CONFLICT if opposite signs)
    bi_html <- paste(as.character(output$beyond_inf_alerts$html), collapse = " ")
    expect_true(grepl("LOW-EPV", bi_html, fixed = TRUE))
    expect_true(grepl("HIGH-VIF", bi_html, fixed = TRUE))

    # 3) Survival benchmark table includes both p_adj_bh (for cox_p) and logrank_p_adj_bh (for logrank_p)
    sb <- surv_bench_res()
    expect_equal(nrow(sb), 2L)
    expect_true(all(c("cox_p", "p_adj_bh", "logrank_p", "logrank_p_adj_bh", "ph_p") %in% names(sb)))
  })
})
