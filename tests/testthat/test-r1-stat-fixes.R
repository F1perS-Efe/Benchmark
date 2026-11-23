# R1 kritik istatistik ve veri sizintisi duzeltmeleri testleri

test_that("R1-1: beyond_infiltration uses repeated CV OOF AUC, detects separation, and computes BH/combined LRT p", {
  set.seed(202)
  n <- 60L
  inf_genes <- c("CD3D", "CD3E", "CD8A", "CD8B")
  noise_genes <- c("N1", "N2", "N3")
  all_genes <- c(inf_genes, noise_genes)

  inf_lat <- rnorm(n)
  y <- stats::rbinom(n, 1L, stats::plogis(1.4 * inf_lat))
  resp <- ifelse(y == 1L, "CR", "PD")

  mat <- matrix(rnorm(length(all_genes) * n, sd = 0.5), nrow = length(all_genes),
                dimnames = list(all_genes, paste0("S", seq_len(n))))
  for (g in inf_genes) mat[g, ] <- inf_lat + rnorm(n, sd = 0.25)

  ds_noise <- list(
    c1 = list(
      expr = mat,
      sample_meta = data.frame(sample_id = colnames(mat), response_label = resp, stringsAsFactors = FALSE),
      info = list(normalization = "TPM", scale = "log2")
    ),
    c2 = list(
      expr = mat + matrix(rnorm(length(all_genes) * n, sd = 0.1), nrow = length(all_genes)),
      sample_meta = data.frame(sample_id = colnames(mat), response_label = resp, stringsAsFactors = FALSE),
      info = list(normalization = "TPM", scale = "log2")
    )
  )

  # Gurultu imzasi OOF delta AUC'yi sisirmemeli ve BH/combined_lrt_p donmeli
  res_n <- test_beyond_infiltration(ds_noise, signature = list(genes = noise_genes), method = "zscore")
  expect_true(all(c("lrt_p_adj_bh", "converged", "separation_flag") %in% names(res_n$per_cohort)))
  expect_true(all(c("combined_lrt_p", "min_lrt_p_adj_bh") %in% names(res_n$summary)))
  expect_true(all(!res_n$per_cohort$adds_independent_value))
  expect_true(all(res_n$per_cohort$delta_auc_combined < 0.05))

  # Tam ayrisma (separation) durumunda separation_flag == TRUE ve adds_independent_value == FALSE
  mat_sep <- mat
  mat_sep[noise_genes, ] <- ifelse(matrix(rep(y == 1L, each = 3L), nrow = 3L), 10, -10)
  ds_sep <- list(
    sep = list(
      expr = mat_sep,
      sample_meta = data.frame(sample_id = colnames(mat_sep), response_label = resp, stringsAsFactors = FALSE),
      info = list(normalization = "TPM", scale = "log2")
    )
  )
  res_s <- test_beyond_infiltration(ds_sep, signature = list(genes = noise_genes), method = "zscore")
  expect_true(isTRUE(res_s$per_cohort$separation_flag[1]))
  expect_false(isTRUE(res_s$per_cohort$adds_independent_value[1]))
})

test_that("R1-2: benchmark and survival scoring produce identical zscore/rank on full cohort without SD label leakage", {
  data(datasets_ici, package = "SigBenchR")
  data(signatures_demo, package = "SigBenchR")

  ds <- datasets_ici$riaz_2017 # 16 SD icerir
  sig <- signatures_demo$IFNG6
  prep <- SigBenchR:::.prepare_cohort(ds, apply_log2 = TRUE, sd_handling = "exclude", filter_response = TRUE)

  for (meth in c("zscore", "rank")) {
    sc_full <- score_signature(prep$expr, genes = sig$genes, method = meth)
    sc_bench <- score_for_dataset(ds, genes = sig$genes, method = meth, sd_handling = "exclude")
    kept_ids <- names(sc_bench$scores)
    expect_equal(sc_bench$scores[kept_ids], sc_full$scores[kept_ids], tolerance = 1e-12)
  }
})

test_that("R1-3: survival orientation consistency, cox.zph PH test, tie handling, and BH FDR", {
  skip_if_not_installed("survival")
  set.seed(33)
  n <- 50L
  sc <- rnorm(n)
  # lower_better: dusuk skor = uzun sagkalim
  time <- exp(3.0 - 0.8 * sc + rnorm(n, sd = 0.4))
  event <- rep(1L, n)
  mat <- rbind(G1 = sc, G2 = sc + rnorm(n, sd = 0.1), G3 = sc + rnorm(n, sd = 0.1))
  colnames(mat) <- paste0("P", seq_len(n))
  ds <- list(
    expr = mat,
    sample_meta = data.frame(
      sample_id = colnames(mat),
      response_label = rep(c("CR", "PD"), each = 25L),
      os_time = time,
      os_event = event,
      stringsAsFactors = FALSE
    ),
    info = list(normalization = "TPM", scale = "log2")
  )

  ev <- evaluate_survival(ds, stats::setNames(sc, colnames(mat)), orient_direction = "lower_better")
  expect_true(all(c("ph_p", "ph_violated", "surv_note") %in% names(ev)))
  # lower_better ile ters cevrilmis z_score -> yuksek yonelimli skor iyi prognoz (HR < 1 ve High medyan > Low medyan)
  expect_lt(ev$hr, 1)
  expect_gt(ev$median_surv_high, ev$median_surv_low)

  # Yogun tie durumu -> logrank_p NA ve surv_note dolu
  sc_tie <- stats::setNames(c(rep(0, 48), 5, 6), colnames(mat))
  ev_tie <- evaluate_survival(ds, sc_tie)
  expect_true(is.na(ev_tie$logrank_p))
  expect_true(grepl("ties|stratum", ev_tie$surv_note, ignore.case = TRUE))

  # run_survival_benchmark BH FDR
  sigs <- list(
    S_pos = list(genes = c("G1", "G2", "G3"), direction = "positive"),
    S_neg = list(genes = c("G1", "G2", "G3"), direction = "negative")
  )
  sb <- run_survival_benchmark(list(c1 = ds), sigs, method = "mean")
  expect_true("p_adj_bh" %in% names(sb))
  expect_true(all(is.finite(sb$p_adj_bh)))
})

test_that("R1-4: power uses se_null critical bound, REML tracks n_cohorts_re, and strong tier requires k >= 5", {
  pw <- estimate_auc_power(20L, 25L, target_auc = 0.70, null_auc = 0.50, alpha = 0.05)
  expect_true(is.finite(pw$se_null))
  z_crit <- stats::qnorm(0.975)
  expected_pw <- (1 - stats::pnorm((0.50 + z_crit * pw$se_null - 0.70) / pw$se_target)) +
    stats::pnorm((0.50 - z_crit * pw$se_null - 0.70) / pw$se_target)
  expect_equal(pw$power, expected_pw, tolerance = 1e-12)

  # Boundary AUC = 1.0 uyarisi ve n_cohorts_re
  df <- data.frame(
    signature = rep("SigB", 3L),
    method = rep("mean", 3L),
    dataset = c("c1", "c2", "c3"),
    n = c(30L, 30L, 30L),
    auc = c(1.0, 0.72, 0.74),
    auc_var = c(0.01, 0.01, 0.01),
    skipped = rep(FALSE, 3L),
    independent = rep(TRUE, 3L),
    low_coverage = rep(FALSE, 3L),
    coverage_status = rep("ok", 3L),
    stringsAsFactors = FALSE
  )
  expect_warning(sm <- summarize_auc(df, method = "mean"), "AUC == 0 or 1")
  expect_identical(sm$n_cohorts_re[1], 2L)

  # k = 4 -> Strong verilmemeli; R_SMALL_K_META tetiklenmeli
  sm4 <- data.frame(
    signature = "SigStrong",
    method = "zscore",
    n_cohorts = 4L,
    auc_n_weighted = 0.78,
    auc_ivw = 0.78,
    auc_re = 0.78,
    auc_re_ci_low = 0.62,
    auc_re_ci_high = 0.89,
    i2_logit = 10.0,
    tau2_reml = 0.01,
    note = NA_character_,
    stringsAsFactors = FALSE
  )
  m4 <- build_report_model(results = data.frame(dataset = "c1", signature = "SigStrong", method = "zscore", auc = 0.78, low_coverage = FALSE), summary_tbl = sm4)
  expect_true(grepl("Medium Priority", m4$tier1$signature_cards$priority_tier[1]))
  expect_true("R_SMALL_K_META" %in% m4$tier1$alerts$rule_id)

  sm5 <- sm4; sm5$n_cohorts <- 5L
  m5 <- build_report_model(results = data.frame(dataset = "c1", signature = "SigStrong", method = "zscore", auc = 0.78, low_coverage = FALSE), summary_tbl = sm5)
  expect_true(grepl("High Priority", m5$tier1$signature_cards$priority_tier[1]))
})
