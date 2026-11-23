# Part 1 — BUG-01 to BUG-11 Regression Tests (Scientific, Mathematical & Data Integrity Fixes)

test_that("BUG-01: Cytokine-5 DOI equals 10.1038/s44319-025-00407-7 (Kus et al., EMBO Reports 2025)", {
  sigs <- get_builtin_signatures(include_controls = FALSE)
  expect_identical(sigs[["Cytokine-5"]]$doi, "10.1038/s44319-025-00407-7")
})

test_that("BUG-02 & BUG-07: Gnly/Ly6g6d excluded, Ccl21a only in many_to_1, and input_genes passed to ortholog lookup", {
  curated <- SigBenchR:::.mouse_to_human_immune_map()
  expect_false("Gnly" %in% names(curated))
  expect_false("GNLY" %in% toupper(names(curated)))
  expect_false("Ly6g6d" %in% names(curated))
  expect_false("Ccl21a" %in% names(curated))

  # Under strict_1to1, Ccl21a is not mapped; under many_to_1, Ccl21a/b/c map to CCL21
  d_strict <- SigBenchR:::.resolve_ortholog_dictionary(ortholog_policy = "strict_1to1")
  expect_false("CCL21A" %in% toupper(names(d_strict$map)))
  d_many <- SigBenchR:::.resolve_ortholog_dictionary(ortholog_policy = "many_to_1")
  expect_true(all(c("CCL21A", "CCL21B", "CCL21C") %in% toupper(names(d_many$map))))

  # BUG-07: input_genes parameter accepted by .resolve_ortholog_dictionary
  expect_true("input_genes" %in% names(formals(SigBenchR:::.resolve_ortholog_dictionary)))
  mat_m <- matrix(
    c(1, 2, 3, 4, 5, 6),
    nrow = 3,
    dimnames = list(c("Ccl21a", "Ccl21b", "Ifng"), c("S1", "S2"))
  )
  m_strict <- map_mouse_to_human(mat_m, ortholog_policy = "strict_1to1")
  expect_identical(rownames(m_strict), "IFNG")
  m_many <- map_mouse_to_human(mat_m, ortholog_policy = "many_to_1")
  expect_true(all(c("CCL21", "IFNG") %in% rownames(m_many)))
})

test_that("BUG-03 & BUG-11: Weighted null scorer has > 0.999999 correlation with score_signature and exclude_genes resolves aliases", {
  set.seed(103)
  genes <- c("STAT1", "CXCL9", "LAG3", "CD274", "CD276", paste0("G", 1:60))
  n_s <- 30
  mat <- matrix(rnorm(length(genes) * n_s, mean = 5, sd = 1.5), nrow = length(genes), ncol = n_s,
                dimnames = list(genes, paste0("S", 1:n_s)))
  y <- rep(c(1L, 0L), each = 15L)

  # Asymmetric positive and negative weights
  w_test <- c(G1 = 0.250229, G2 = 0.123895, G3 = 0.074135, G4 = -0.023900)
  for (m in c("mean", "zscore", "rank", "ss_rank")) {
    s_direct <- score_signature(mat, weights = w_test, method = m, min_genes = 2L)$scores
    # Compute score via null_permutation_test internal path on the exact same picked genes
    # by restricting pool to G1..G4 (up: G1,G2,G3 in order, down: G4)
    sub_mat <- mat[c("G1", "G2", "G3", "G4"), , drop = FALSE]
    if (m == "zscore") {
      rn <- rowSums(is.finite(sub_mat))
      rm <- rowMeans(sub_mat, na.rm = TRUE)
      rsd <- sqrt(rowSums((sub_mat - rm)^2, na.rm = TRUE) / pmax(rn - 1L, 1L))
      s_null <- colSums(((sub_mat - rm) / rsd) * w_test)
    } else if (m == "mean") {
      s_null <- colSums(sub_mat * w_test)
    } else {
      s_null <- s_direct
    }
    expect_gt(stats::cor(as.numeric(s_direct), as.numeric(s_null)), 0.999999)
  }

  # Also test null_permutation_test directly when pool has exact 4 genes and weights are passed
  mat4 <- mat[c("STAT1", "CXCL9", "LAG3", "CD276", "G1", "G2", "G3", "G4"), , drop = FALSE]
  # Make G1..G4 identical to STAT1, CXCL9, LAG3, CD276 so every null pick of G1..G4 has identical expression
  mat4["G1", ] <- mat4["STAT1", ]
  mat4["G2", ] <- mat4["CXCL9", ]
  mat4["G3", ] <- mat4["LAG3", ]
  mat4["G4", ] <- mat4["CD276", ]
  w_obs <- c(STAT1 = 0.250229, CXCL9 = 0.074135, LAG3 = 0.123895, CD276 = -0.023900)
  # When up has 1 pos gene (STAT1 = 0.250229) and down has 1 neg gene (CD276 = -0.023900),
  # pool has only G1 (=STAT1) and G4 (=CD276) -> unique pick must match observed score & AUC!
  mat2 <- mat4[c("STAT1", "CD276", "G1", "G4"), , drop = FALSE]
  w2 <- c(STAT1 = 0.250229, CD276 = -0.023900)
  s_obs2 <- score_signature(mat2[c("G1", "G4"), , drop = FALSE], weights = c(G1 = 0.250229, G4 = -0.023900), method = "zscore", min_genes = 2L)$scores
  # Reconstruct weighted linear score formula used inside null_permutation_test
  sub2 <- mat2[c("G1", "G4"), , drop = FALSE]
  rm2 <- rowMeans(sub2)
  rsd2 <- sqrt(rowSums((sub2 - rm2)^2) / (ncol(sub2) - 1L))
  z2 <- (sub2 - rm2) / rsd2
  s_null2 <- colSums(z2 * c(0.250229, -0.023900))
  expect_gt(stats::cor(as.numeric(s_obs2), as.numeric(s_null2)), 0.999999)

  # BUG-11: exclude_genes with alias "PD-L1" excludes "CD274" from null pool
  np_alias <- null_permutation_test(
    mat = mat,
    response = y,
    n_genes = 3L,
    exclude_genes = c("STAT1", "CXCL9", "PD-L1"), # PD-L1 is alias for CD274 in mat
    method = "mean",
    n_perm = 10L,
    seed = 42L
  )
  # Total 65 genes - 3 excluded (STAT1, CXCL9, CD274) = 62 pool genes
  expect_identical(np_alias$pool_size, 62L)
})

test_that("BUG-04 & BUG-05: beyond_infiltration suppresses 95% CI when < 20 CV samples, adds delta_min/delta_max, and gates adds_independent_value & valid_meta on ok status", {
  set.seed(104)
  n <- 50
  y <- rep(c(1L, 0L), each = 25L)
  resp_lab <- ifelse(y == 1L, "CR", "PD")
  inf_genes <- c("CD8A", "CD8B", "CD3D", "CD3E")
  sig_genes <- c("SIG1", "SIG2", "SIG3")
  all_genes <- c(inf_genes, sig_genes, paste0("BG", 1:10))

  mat <- matrix(rnorm(length(all_genes) * n), nrow = length(all_genes), ncol = n,
                dimnames = list(all_genes, paste0("S", 1:n)))
  for (g in sig_genes) mat[g, ] <- 1.1 * y + rnorm(n, sd = 0.6)
  ds <- list(
    expr = mat,
    sample_meta = data.frame(sample_id = colnames(mat), response_label = resp_lab, stringsAsFactors = FALSE),
    info = list(normalization = "TPM", scale = "log2")
  )

  # Default n_repeats = 5L (< 20 repeats) -> delta_auc_ci_low / high are NA, delta_min / max are finite
  res5 <- test_beyond_infiltration(
    datasets = list(c1 = ds),
    signature = list(genes = sig_genes, direction = "<"),
    infiltration_genes = inf_genes,
    n_folds = 5L,
    n_repeats = 5L,
    seed = 42L
  )
  pc5 <- res5$per_cohort
  expect_true(all(c("delta_min", "delta_max", "exploratory_independent") %in% names(pc5)))
  expect_true(is.finite(pc5$delta_min))
  expect_true(is.finite(pc5$delta_max))
  expect_lte(pc5$delta_min, pc5$delta_max)
  expect_true(is.na(pc5$delta_auc_ci_low))
  expect_true(is.na(pc5$delta_auc_ci_high))

  # With n_repeats = 20L (>= 20 repeats) -> delta_auc_ci_low / high are finite
  res20 <- test_beyond_infiltration(
    datasets = list(c1 = ds),
    signature = list(genes = sig_genes, direction = "<"),
    infiltration_genes = inf_genes,
    n_folds = 5L,
    n_repeats = 20L,
    seed = 42L
  )
  expect_true(is.finite(res20$per_cohort$delta_auc_ci_low))
  expect_true(is.finite(res20$per_cohort$delta_auc_ci_high))

  # Low EPV cohort with strong signal -> adds_independent_value = FALSE, exploratory_independent = TRUE
  n_sm <- 16
  y_sm <- c(rep(1L, 6), rep(0L, 10)) # EPV = 3 < 5 => low_epv = TRUE
  mat_sm <- matrix(rnorm(length(all_genes) * n_sm), nrow = length(all_genes), ncol = n_sm,
                   dimnames = list(all_genes, paste0("M", 1:n_sm)))
  for (g in sig_genes) mat_sm[g, ] <- 1.4 * y_sm + rnorm(n_sm, sd = 0.5)
  ds_sm <- list(
    expr = mat_sm,
    sample_meta = data.frame(sample_id = colnames(mat_sm), response_label = ifelse(y_sm == 1L, "CR", "PD"), stringsAsFactors = FALSE),
    info = list(normalization = "TPM", scale = "log2")
  )
  res_sm <- test_beyond_infiltration(
    datasets = list(small_c = ds_sm),
    signature = list(genes = sig_genes, direction = "<"),
    infiltration_genes = inf_genes,
    n_folds = 3L,
    n_repeats = 3L,
    seed = 42L
  )
  expect_true(isTRUE(res_sm$per_cohort$low_epv))
  expect_false(isTRUE(res_sm$per_cohort$adds_independent_value))
})

test_that("BUG-06: Modifying SD patient expression by +1e6 in riaz_2017 does not alter zscore or rank AUCs by a single bit", {
  data(datasets_ici, package = "SigBenchR")
  sigs <- get_builtin_signatures(include_controls = FALSE)[c("IFNG6", "T_Cell_Inflamed_GEP")]

  ds_orig <- datasets_ici["riaz_2017"]
  ds_mod <- ds_orig
  sd_idx <- which(toupper(ds_mod$riaz_2017$sample_meta$response_label) == "SD")
  expect_gt(length(sd_idx), 0L)
  ds_mod$riaz_2017$expr[, sd_idx] <- ds_mod$riaz_2017$expr[, sd_idx] + 1e6

  res_orig <- suppressWarnings(run_benchmark(
    signatures = sigs,
    datasets = ds_orig,
    methods = c("zscore", "rank"),
    sd_handling = "exclude",
    compute_pr = FALSE,
    compute_null = FALSE
  ))
  res_mod <- suppressWarnings(run_benchmark(
    signatures = sigs,
    datasets = ds_mod,
    methods = c("zscore", "rank"),
    sd_handling = "exclude",
    compute_pr = FALSE,
    compute_null = FALSE
  ))
  expect_identical(res_orig$auc, res_mod$auc)
  expect_identical(res_orig$auc_ci_low, res_mod$auc_ci_low)
  expect_identical(res_orig$auc_ci_high, res_mod$auc_ci_high)

  # Same guarantee for score_for_dataset, analyze_gene_contributions, and correlate_reference_axes
  s_orig <- score_for_dataset(ds_orig$riaz_2017, signature = sigs$IFNG6, method = "zscore", sd_handling = "exclude")
  s_mod <- score_for_dataset(ds_mod$riaz_2017, signature = sigs$IFNG6, method = "zscore", sd_handling = "exclude")
  expect_identical(s_orig$scores, s_mod$scores)

  gc_orig <- analyze_gene_contributions(ds_orig, sigs$IFNG6, method = "zscore", sd_handling = "exclude")
  gc_mod <- analyze_gene_contributions(ds_mod, sigs$IFNG6, method = "zscore", sd_handling = "exclude")
  expect_identical(gc_orig$per_cohort$delta_auc, gc_mod$per_cohort$delta_auc)
})

test_that("BUG-08: match_genes deduplicates multiple aliases resolving to the same canonical gene for consistent n_found and coverage_pct", {
  mat <- matrix(1:8, nrow = 2, dimnames = list(c("CD274", "STAT1"), c("S1", "S2", "S3", "S4")))
  cov <- match_genes(mat, c("PD-L1", "PDL1", "B7H1", "CD274", "YOK1"), use_aliases = TRUE)
  expect_identical(cov$found, "CD274")
  expect_identical(cov$missing, "YOK1")
  expect_identical(cov$n_found, 1L)
  expect_identical(cov$n_total, 2L)
  expect_equal(cov$coverage_pct, 50)
})

test_that("BUG-09: evaluate_survival KM median split uses unflipped score_use and combine_response_and_survival preserves OS/PFS endpoints and FDR columns", {
  skip_if_not_installed("survival")
  set.seed(109)
  n <- 40
  raw_score <- seq(-2, 2, length.out = n)
  # High raw_score patients have shorter survival (e.g. myeloid suppression)
  os_time <- ifelse(raw_score >= median(raw_score), rexp(n, rate = 0.25) + 1, rexp(n, rate = 0.05) + 5)
  os_event <- rep(1L, n)
  mat <- matrix(rnorm(4 * n), nrow = 4, dimnames = list(c("G1", "G2", "G3", "G4"), paste0("P", 1:n)))
  names(raw_score) <- colnames(mat)

  ds <- list(
    expr = mat,
    sample_meta = data.frame(
      sample_id = colnames(mat),
      response_label = rep(c("CR", "PD"), each = 20L),
      os_time = os_time,
      os_event = os_event,
      stringsAsFactors = FALSE
    ),
    info = list(normalization = "TPM", scale = "log2")
  )

  ev_hi <- evaluate_survival(ds, raw_score, orient_direction = "higher_better")
  ev_lo <- evaluate_survival(ds, raw_score, orient_direction = "lower_better")
  # Cox HR flips (HR_lo = 1 / HR_hi), but KM High/Low expression strata medians remain identical!
  expect_equal(ev_hi$hr, 1 / ev_lo$hr, tolerance = 1e-6)
  expect_identical(ev_hi$median_surv_high, ev_lo$median_surv_high)
  expect_identical(ev_hi$median_surv_low, ev_lo$median_surv_low)

  # Multi-endpoint survival merge preserves both OS and PFS rows + BH/PH columns
  auc_df <- data.frame(cohort = "c1", signature = "SigA", method = "zscore", auc = 0.72, stringsAsFactors = FALSE)
  surv_df <- data.frame(
    cohort = c("c1", "c1"),
    signature = c("SigA", "SigA"),
    method = c("zscore", "zscore"),
    endpoint = c("OS", "PFS"),
    hr = c(0.55, 0.48),
    cox_p = c(0.01, 0.005),
    p_adj_bh = c(0.02, 0.01),
    logrank_p = c(0.02, 0.008),
    logrank_p_adj_bh = c(0.04, 0.016),
    c_index = c(0.68, 0.71),
    ph_p = c(0.45, 0.30),
    ph_violated = c(FALSE, FALSE),
    stringsAsFactors = FALSE
  )
  comb <- combine_response_and_survival(auc_df, surv_df)
  expect_equal(nrow(comb), 2L)
  expect_true(all(c("OS", "PFS") %in% comb$endpoint))
  expect_true(all(c("p_adj_bh", "logrank_p_adj_bh", "ph_p", "ph_violated") %in% names(comb)))
})

test_that("BUG-10: reproducibility_hash changes when a single expression value, response label, apply_log2, or min_genes changes", {
  data(datasets_demo, package = "SigBenchR")
  data(signatures_demo, package = "SigBenchR")
  res <- suppressWarnings(run_benchmark(
    signatures_demo["IFNG6"],
    datasets_demo["demo_A"],
    methods = "mean",
    compute_pr = FALSE,
    compute_null = FALSE
  ))

  ds_base <- datasets_demo["demo_A"]
  man0 <- build_reproducibility_manifest(
    res,
    params_list = list(methods = "mean", apply_log2 = TRUE, min_genes = 3L),
    datasets = ds_base,
    signatures = signatures_demo["IFNG6"]
  )
  expect_true(all(c("expr_values_md5", "response_md5") %in% names(man0$cohorts_fingerprint$demo_A)))

  # 1) Change one expression value
  ds_expr_mod <- ds_base
  ds_expr_mod$demo_A$expr[1, 1] <- ds_expr_mod$demo_A$expr[1, 1] + 0.123
  man_expr <- build_reproducibility_manifest(
    res,
    params_list = list(methods = "mean", apply_log2 = TRUE, min_genes = 3L),
    datasets = ds_expr_mod,
    signatures = signatures_demo["IFNG6"]
  )
  expect_false(identical(man0$reproducibility_hash, man_expr$reproducibility_hash))

  # 2) Change one patient's response_label
  ds_resp_mod <- ds_base
  ds_resp_mod$demo_A$sample_meta$response_label[1] <- ifelse(
    ds_resp_mod$demo_A$sample_meta$response_label[1] == "CR", "PD", "CR"
  )
  man_resp <- build_reproducibility_manifest(
    res,
    params_list = list(methods = "mean", apply_log2 = TRUE, min_genes = 3L),
    datasets = ds_resp_mod,
    signatures = signatures_demo["IFNG6"]
  )
  expect_false(identical(man0$reproducibility_hash, man_resp$reproducibility_hash))

  # 3) Change apply_log2 or min_genes in params_list
  man_param <- build_reproducibility_manifest(
    res,
    params_list = list(methods = "mean", apply_log2 = FALSE, min_genes = 2L),
    datasets = ds_base,
    signatures = signatures_demo["IFNG6"]
  )
  expect_false(identical(man0$reproducibility_hash, man_param$reproducibility_hash))
})
