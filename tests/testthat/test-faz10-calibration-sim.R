# FAZ 10 — Monte Carlo Type-I Error & Null Calibration Simulation Tests

test_that("10.1: Null permutation test Type-I error rate is calibrated under pure noise", {
  n_sim <- 100L
  n_genes_tot <- 300L
  n_samp <- 40L
  sig_genes <- paste0("G", seq_len(5L))

  p_vals <- withr::with_seed(1001L, {
    vapply(seq_len(n_sim), function(i) {
      y <- rep(c(1L, 0L), each = n_samp / 2L)
      mat <- matrix(
        stats::rnorm(n_genes_tot * n_samp),
        nrow = n_genes_tot,
        ncol = n_samp,
        dimnames = list(paste0("G", seq_len(n_genes_tot)), paste0("S", seq_len(n_samp)))
      )
      res_np <- null_permutation_test(
        mat = mat,
        response = y,
        n_genes = 5L,
        exclude_genes = sig_genes,
        method = "zscore",
        direction = "<",
        n_perm = 99L,
        seed = 1000L + i,
        match_expression = TRUE,
        signature_genes = sig_genes,
        min_genes = 3L
      )
      res_np$null_p
    }, numeric(1))
  })

  expect_equal(length(p_vals), n_sim)
  expect_true(all(is.finite(p_vals)))
  fp_rate <- mean(p_vals < 0.05)
  # Exact uniform nominal alpha = 0.05 (finite B=99 => E[rate] = 0.04); allow Monte Carlo SE margin <= 0.10
  expect_lte(fp_rate, 0.10)
  # Mean null p-value under null should stay centered near 0.50 (never collapsed toward 0 or 1)
  expect_gt(mean(p_vals), 0.35)
  expect_lt(mean(p_vals), 0.65)
})

test_that("10.2: test_beyond_infiltration pure-noise signatures do not inflate fold-inside CV delta_auc or false positives", {
  n_sim <- 100L
  n_samp <- 60L
  inf_genes <- c("CD8A", "CD8B", "CD3D", "CD3E")
  noise_sig <- c("NOISE1", "NOISE2", "NOISE3", "NOISE4")
  bg_genes <- paste0("BG", seq_len(12L))
  all_genes <- c(inf_genes, noise_sig, bg_genes)

  sim_out <- withr::with_seed(1002L, {
    deltas <- numeric(n_sim)
    added <- logical(n_sim)
    lrt_ps <- numeric(n_sim)
    for (i in seq_len(n_sim)) {
      inf_lat <- stats::rnorm(n_samp)
      prob <- stats::plogis(1.1 * inf_lat)
      y <- stats::rbinom(n_samp, size = 1L, prob = prob)
      # Guarantee both classes present
      if (sum(y == 1L) < 10L || sum(y == 0L) < 10L) {
        y <- rep(c(1L, 0L), each = n_samp / 2L)
      }
      mat <- matrix(
        stats::rnorm(length(all_genes) * n_samp),
        nrow = length(all_genes),
        ncol = n_samp,
        dimnames = list(all_genes, paste0("S", seq_len(n_samp)))
      )
      for (g in inf_genes) {
        mat[g, ] <- inf_lat + stats::rnorm(n_samp, sd = 0.25)
      }
      ds <- list(
        expr = mat,
        sample_meta = data.frame(
          sample_id = colnames(mat),
          response_label = ifelse(y == 1L, "CR", "PD"),
          stringsAsFactors = FALSE
        ),
        info = list(normalization = "TPM", scale = "log2")
      )
      res_bi <- test_beyond_infiltration(
        datasets = list(sim = ds),
        signature = list(genes = noise_sig, direction = "<"),
        infiltration_genes = inf_genes,
        method = "zscore",
        n_folds = 4L,
        n_repeats = 2L,
        seed = 2000L + i
      )
      pc <- res_bi$per_cohort
      deltas[i] <- pc$delta_auc_combined[1]
      added[i] <- isTRUE(pc$adds_independent_value[1])
      lrt_ps[i] <- pc$lrt_p[1]
    }
    list(deltas = deltas, added = added, lrt_ps = lrt_ps)
  })

  ok_d <- sim_out$deltas[is.finite(sim_out$deltas)]
  expect_equal(length(ok_d), n_sim)
  # Average out-of-fold CV delta-AUC for noise added on top of true predictor is <= 0.02 (near <= 0)
  expect_lte(mean(ok_d), 0.02)
  # False positive rate of claiming independent value beyond infiltration is controlled (<= 0.05)
  expect_lte(mean(sim_out$added), 0.05)
  expect_lte(mean(sim_out$lrt_ps < 0.05, na.rm = TRUE), 0.10)
})
