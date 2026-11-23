# FAZ 4 — Statistical Engine, Null Permutation & Beyond Infiltration Reliability Tests

test_that("4.1 & 4.2: test_beyond_infiltration computes fold-inside preprocessing and fold-level paired CV AUC with uncertainty", {
  set.seed(401)
  n <- 60
  y <- c(rep(1L, 30), rep(0L, 30))
  resp_lab <- ifelse(y == 1L, "CR", "PD")
  inf_genes <- c("CD8A", "CD8B", "CD3D", "CD3E")
  sig_genes <- c("SIG1", "SIG2", "SIG3")
  bg_genes <- paste0("BG", seq_len(20))
  all_genes <- c(inf_genes, sig_genes, bg_genes)

  mat <- matrix(rnorm(length(all_genes) * n), nrow = length(all_genes), ncol = n)
  rownames(mat) <- all_genes
  colnames(mat) <- paste0("S", seq_len(n))

  # Inject moderate infiltration signal and independent signature signal
  inf_lat <- rnorm(n, mean = 0.6 * y, sd = 1)
  sig_lat <- rnorm(n, mean = 1.2 * y, sd = 1)
  for (g in inf_genes) mat[g, ] <- inf_lat + rnorm(n, sd = 0.3)
  for (g in sig_genes) mat[g, ] <- sig_lat + rnorm(n, sd = 0.3)

  ds <- list(
    expr = mat,
    sample_meta = data.frame(
      sample_id = colnames(mat),
      response_label = resp_lab,
      stringsAsFactors = FALSE
    ),
    info = list(normalization = "TPM", scale = "log2")
  )

  # Test both zscore/rank (cohort-dependent) and ss_rank (single-sample)
  for (m in c("zscore", "rank", "ss_rank")) {
    res <- test_beyond_infiltration(
      datasets = list(cohort_a = ds),
      signature = list(genes = sig_genes, direction = "<"),
      method = m,
      infiltration_genes = inf_genes,
      n_folds = 5L,
      n_repeats = 3L,
      seed = 42L
    )
    pc <- res$per_cohort
    expect_equal(nrow(pc), 1L)
    expect_true(all(c(
      "sd_delta_auc", "delta_auc_ci_low", "delta_auc_ci_high",
      "lrt_p_adj_bh", "epv", "low_epv", "vif", "high_collinearity", "analysis_status"
    ) %in% names(pc)))
    expect_true(is.finite(pc$delta_auc_combined))
    expect_true(is.finite(pc$sd_delta_auc))
    expect_true(is.finite(pc$delta_auc_ci_low))
    expect_true(is.finite(pc$delta_auc_ci_high))
    expect_lte(pc$delta_auc_ci_low, pc$delta_auc_ci_high)
    expect_gt(pc$delta_auc_combined, 0)
  }
})

test_that("4.3 & 4.4: BH-adjusted LRT p-value gates adds_independent_value and directional Stouffer Z cancels opposite-sign cohorts", {
  set.seed(404)
  n <- 70
  y <- c(rep(1L, 35), rep(0L, 35))
  resp_lab <- ifelse(y == 1L, "CR", "PD")
  inf_genes <- c("CD8A", "CD8B", "CD3D", "CD3E")
  sig_genes <- c("SIG1", "SIG2", "SIG3")
  all_genes <- c(inf_genes, sig_genes, paste0("BG", 1:10))

  make_cohort <- function(beta_sign) {
    mat <- matrix(rnorm(length(all_genes) * n), nrow = length(all_genes), ncol = n)
    rownames(mat) <- all_genes
    colnames(mat) <- paste0("S", seq_len(n))
    sig_lat <- beta_sign * (1.5 * y) + rnorm(n, sd = 0.6)
    for (g in sig_genes) mat[g, ] <- sig_lat + rnorm(n, sd = 0.2)
    list(
      expr = mat,
      sample_meta = data.frame(
        sample_id = colnames(mat),
        response_label = resp_lab,
        stringsAsFactors = FALSE
      ),
      info = list(normalization = "TPM", scale = "log2")
    )
  }

  # Cohort 1 has strong positive independent effect (+); Cohort 2 has strong negative effect (-)
  ds_pos <- make_cohort(+1)
  ds_neg <- make_cohort(-1)

  res <- test_beyond_infiltration(
    datasets = list(pos_cohort = ds_pos, neg_cohort = ds_neg),
    signature = list(genes = sig_genes, direction = "<"),
    method = "zscore",
    infiltration_genes = inf_genes,
    n_folds = 5L,
    n_repeats = 2L,
    seed = 99L
  )

  expect_true(isTRUE(res$summary$direction_conflict))
  expect_true(res$per_cohort$beta_adjusted[1] > 0)
  expect_true(res$per_cohort$beta_adjusted[2] < 0)
  # Unsigned Fisher omnibus p is tiny because both LRTs are significant in opposite directions,
  # whereas directional Stouffer combined_lrt_p is non-significant because opposing effects cancel
  expect_lt(res$summary$fisher_omnibus_p, 0.01)
  expect_gt(res$summary$combined_lrt_p, 0.05)
})

test_that("4.5: Small-sample EPV and high multicollinearity (VIF) are flagged accurately", {
  set.seed(405)
  inf_genes <- c("CD8A", "CD8B", "CD3D", "CD3E")
  sig_genes <- c("SIG1", "SIG2", "SIG3")
  all_genes <- c(inf_genes, sig_genes)

  # Case A: Extreme collinearity (r > 0.99 => VIF > 10 => failed_collinearity)
  n <- 50
  y <- c(rep(1L, 25), rep(0L, 25))
  resp_lab <- ifelse(y == 1L, "CR", "PD")
  lat <- rnorm(n, mean = 0.8 * y)
  mat_col <- matrix(0, nrow = length(all_genes), ncol = n)
  rownames(mat_col) <- all_genes
  colnames(mat_col) <- paste0("C", seq_len(n))
  for (g in all_genes) mat_col[g, ] <- lat + rnorm(n, sd = 0.02)

  ds_col <- list(
    expr = mat_col,
    sample_meta = data.frame(sample_id = colnames(mat_col), response_label = resp_lab, stringsAsFactors = FALSE),
    info = list(normalization = "TPM", scale = "log2")
  )
  res_col <- test_beyond_infiltration(
    datasets = list(collinear_cohort = ds_col),
    signature = sig_genes,
    infiltration_genes = inf_genes
  )
  expect_true(isTRUE(res_col$per_cohort$high_collinearity))
  expect_gt(res_col$per_cohort$vif, 10)
  expect_equal(res_col$per_cohort$analysis_status, "failed_collinearity")
  expect_false(res_col$per_cohort$adds_independent_value)

  # Case B: Low EPV (n_pos = 6 => EPV = 3 < 5 => low_epv == TRUE)
  n_small <- 20
  y_small <- c(rep(1L, 6), rep(0L, 14))
  resp_sm <- ifelse(y_small == 1L, "CR", "PD")
  mat_sm <- matrix(rnorm(length(all_genes) * n_small), nrow = length(all_genes), ncol = n_small)
  rownames(mat_sm) <- all_genes
  colnames(mat_sm) <- paste0("M", seq_len(n_small))
  ds_sm <- list(
    expr = mat_sm,
    sample_meta = data.frame(sample_id = colnames(mat_sm), response_label = resp_sm, stringsAsFactors = FALSE),
    info = list(normalization = "TPM", scale = "log2")
  )
  res_sm <- test_beyond_infiltration(
    datasets = list(small_cohort = ds_sm),
    signature = sig_genes,
    infiltration_genes = inf_genes
  )
  expect_equal(res_sm$per_cohort$epv, 3)
  expect_true(isTRUE(res_sm$per_cohort$low_epv))
  expect_true(grepl("low_epv", res_sm$per_cohort$analysis_status))
})

test_that("4.6: null_permutation_test supports 2D mean+SD matching, continuous weights, and background_type", {
  set.seed(406)
  n_genes_tot <- 300
  n_samp <- 30
  y <- c(rep(1L, 15), rep(0L, 15))
  mat <- matrix(rnorm(n_genes_tot * n_samp), nrow = n_genes_tot, ncol = n_samp)
  rownames(mat) <- paste0("G", seq_len(n_genes_tot))
  colnames(mat) <- paste0("S", seq_len(n_samp))

  # Assign distinct mean + dispersion tiers across the matrix
  for (i in seq_len(n_genes_tot)) {
    m_shift <- (i %% 5) * 3
    s_scale <- 0.5 + (i %% 3) * 1.5
    mat[i, ] <- rnorm(n_samp, mean = m_shift, sd = s_scale)
  }

  w_sig <- c(G1 = 0.40, G2 = 0.25, G3 = 0.15, G4 = -0.30)
  np <- null_permutation_test(
    mat = mat,
    response = y,
    n_genes = 3L,
    n_genes_down = 1L,
    exclude_genes = names(w_sig),
    method = "zscore",
    direction = "<",
    n_perm = 40L,
    seed = 42L,
    match_expression = TRUE,
    signature_genes = c("G1", "G2", "G3"),
    signature_genes_down = "G4",
    weights = w_sig,
    min_genes = 2L
  )

  expect_true(is.finite(np$obs_auc))
  expect_equal(length(np$null_aucs), 40L)
  expect_true(is.finite(np$null_p))
  expect_equal(np$background_type, "subsampled_panel_cohort")

  # Reproducibility with same seed
  np2 <- null_permutation_test(
    mat = mat,
    response = y,
    n_genes = 3L,
    n_genes_down = 1L,
    exclude_genes = names(w_sig),
    method = "zscore",
    direction = "<",
    n_perm = 40L,
    seed = 42L,
    match_expression = TRUE,
    signature_genes = c("G1", "G2", "G3"),
    signature_genes_down = "G4",
    weights = w_sig,
    min_genes = 2L
  )
  expect_identical(np$null_aucs, np2$null_aucs)
  expect_identical(np$null_p, np2$null_p)
})

test_that("4.7: adjust_null_p provides within_method and global BH adjustment; run_survival_benchmark separates cox_p and logrank_p BH families", {
  df <- data.frame(
    dataset = rep("ds1", 4),
    signature = c("S1", "S2", "S1", "S2"),
    method = c("mean", "mean", "zscore", "zscore"),
    null_p = c(0.01, 0.04, 0.02, 0.20),
    stringsAsFactors = FALSE
  )

  adj_w <- adjust_null_p(df, scope = "within_method")
  expect_true(all(c("null_p_adj_bh", "null_p_adj_bh_global") %in% names(adj_w)))
  expect_equal(adj_w$null_p_adj_bh[1:2], stats::p.adjust(c(0.01, 0.04), method = "BH"))
  expect_equal(adj_w$null_p_adj_bh_global, stats::p.adjust(df$null_p, method = "BH"))

  adj_g <- adjust_null_p(df, scope = "global")
  expect_equal(adj_g$null_p_adj_bh, adj_g$null_p_adj_bh_global)

  # Survival benchmark: cox_p -> p_adj_bh, logrank_p -> logrank_p_adj_bh, ph_p remains raw diagnostic
  set.seed(407)
  n <- 40
  mat <- matrix(rnorm(10 * n), nrow = 10, ncol = n)
  rownames(mat) <- paste0("G", 1:10)
  colnames(mat) <- paste0("P", 1:n)
  ds_surv <- list(
    expr = mat,
    sample_meta = data.frame(
      sample_id = colnames(mat),
      response_label = sample(c("CR", "PD"), n, replace = TRUE),
      os_time = rexp(n, rate = 0.1) + 0.5,
      os_event = sample(c(0L, 1L), n, replace = TRUE, prob = c(0.3, 0.7)),
      stringsAsFactors = FALSE
    ),
    info = list(normalization = "TPM", scale = "log2")
  )
  sigs <- list(
    SigA = list(genes = c("G1", "G2", "G3"), direction = "positive"),
    SigB = list(genes = c("G4", "G5", "G6"), direction = "positive")
  )
  s_res <- run_survival_benchmark(
    datasets = list(c1 = ds_surv),
    signatures = sigs,
    method = "zscore"
  )
  expect_equal(nrow(s_res), 2L)
  expect_true("p_adj_bh" %in% names(s_res))
  expect_true("logrank_p_adj_bh" %in% names(s_res))
  expect_equal(s_res$p_adj_bh, stats::p.adjust(s_res$cox_p, method = "BH"))
  expect_equal(s_res$logrank_p_adj_bh, stats::p.adjust(s_res$logrank_p, method = "BH"))
})

test_that("4.8: pr_auc handles heavy ties and small positive class with explicit diagnostics", {
  # Heavy ties (only 3 distinct values across 30 samples) and small positive class (n_pos = 3 < 5)
  scores_tied <- c(rep(1.0, 2), rep(0.5, 10), rep(0.0, 18))
  resp_small <- c(rep(1L, 3), rep(0L, 27))

  res <- pr_auc(scores_tied, resp_small, direction = "<", n_boot = 50L, seed = 42L)
  expect_true(is.finite(res$pr_auc))
  expect_gt(res$tie_fraction, 0.8)
  expect_true(isTRUE(res$small_positive_class))
  expect_true(isTRUE(res$ci_unreliable))
  expect_true(is.integer(res$n_boot_valid))
})
