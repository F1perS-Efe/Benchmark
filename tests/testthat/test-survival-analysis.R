# Survival / clinical benefit (Cox PH + Kaplan-Meier) regressions

.make_surv_dataset <- function(n = 40L,
                               seed = 42L,
                               score_effect = TRUE,
                               n_events = NULL,
                               constant_score = FALSE) {
  set.seed(seed)
  genes <- paste0("G", 1:6)
  samples <- paste0("S", seq_len(n))
  # Score driver gene: high expression -> high signature score
  mat <- matrix(
    rnorm(length(genes) * n),
    nrow = length(genes),
    dimnames = list(genes, samples)
  )
  latent <- seq_len(n) / n
  if (isTRUE(constant_score)) {
    mat["G1", ] <- 1
    mat["G2", ] <- 1
    mat["G3", ] <- 1
  } else if (isTRUE(score_effect)) {
    mat["G1", ] <- latent * 4
    mat["G2", ] <- latent * 4 + rnorm(n, sd = 0.05)
    mat["G3", ] <- latent * 4 + rnorm(n, sd = 0.05)
  }

  # High score -> longer survival (lower hazard)
  if (isTRUE(constant_score)) {
    score_proxy <- rep(1, n)
  } else {
    score_proxy <- latent
  }
  # Weibull-like times via exponential hazard
  haz <- exp(-3 * score_proxy)
  os_time <- stats::rexp(n, rate = pmax(haz, 1e-6))
  os_event <- rep(1L, n)
  if (!is.null(n_events)) {
    # Keep only n_events as events; censor the rest (set event=0, keep times)
    os_event <- rep(0L, n)
    os_event[seq_len(as.integer(n_events))] <- 1L
  }

  list(
    expr = mat,
    sample_meta = data.frame(
      sample_id = samples,
      response_label = rep(c("CR", "PD"), length.out = n),
      os_time = os_time,
      os_event = os_event,
      stringsAsFactors = FALSE
    ),
    info = list(
      normalization = "TPM",
      scale = "log2",
      dataset_name = "synth_surv"
    )
  )
}

test_that("evaluate_survival recovers protective association (hr < 1)", {
  ds <- .make_surv_dataset(n = 40L, seed = 42L, score_effect = TRUE)
  scored <- score_signature(
    ds$expr,
    genes = c("G1", "G2", "G3"),
    method = "zscore",
    min_genes = 2L
  )
  expect_false(isTRUE(scored$skipped))
  out <- evaluate_survival(
    dataset = ds,
    score = scored$scores,
    time_col = "os_time",
    event_col = "os_event",
    orient_direction = "higher_better"
  )
  expect_equal(nrow(out), 1L)
  expect_false(isTRUE(out$small_sample[1]))
  expect_lt(out$hr[1], 1)
  expect_lt(out$cox_p[1], 0.05)
  expect_lt(out$logrank_p[1], 0.05)
  expect_gt(out$c_index[1], 0.65)
})

test_that("evaluate_survival small_sample for few events and zero-variance score", {
  ds_few <- .make_surv_dataset(n = 40L, seed = 1L, n_events = 2L)
  scored <- score_signature(
    ds_few$expr,
    genes = c("G1", "G2", "G3"),
    method = "mean",
    min_genes = 2L
  )
  out_few <- evaluate_survival(ds_few, scored$scores)
  expect_true(isTRUE(out_few$small_sample[1]))
  expect_true(is.na(out_few$hr[1]))

  ds_const <- .make_surv_dataset(n = 40L, seed = 2L, constant_score = TRUE)
  scored_c <- score_signature(
    ds_const$expr,
    genes = c("G1", "G2", "G3"),
    method = "mean",
    min_genes = 2L
  )
  # Force constant score vector
  sc <- scored_c$scores
  sc[] <- 1
  out_const <- evaluate_survival(ds_const, sc)
  expect_true(isTRUE(out_const$small_sample[1]))
  expect_true(is.na(out_const$hr[1]))
})

test_that("validate_dataset rejects bad survival columns; datasets_ici still ok", {
  ds <- .make_surv_dataset(n = 12L, seed = 3L)
  expect_silent(validate_dataset(ds, name = "ok_surv"))

  ds_neg <- ds
  ds_neg$sample_meta$os_time[1] <- -1
  expect_error(validate_dataset(ds_neg, name = "bad_time"), "negatif|os_time")

  ds_ev <- ds
  ds_ev$sample_meta$os_event[1] <- 2
  expect_error(validate_dataset(ds_ev, name = "bad_event"), "0, 1|os_event")

  data(datasets_ici, package = "SigBenchR")
  for (nm in names(datasets_ici)) {
    expect_silent(validate_dataset(datasets_ici[[nm]], name = nm))
  }
})

test_that("run_survival_benchmark skips non-survival cohorts; combine merges AUC+HR", {
  set.seed(9L)
  genes <- paste0("G", 1:5)
  n <- 20L
  mat <- matrix(
    rnorm(length(genes) * n),
    nrow = length(genes),
    dimnames = list(genes, paste0("S", seq_len(n)))
  )
  ds_nosurv <- list(
    expr = mat,
    sample_meta = data.frame(
      sample_id = colnames(mat),
      response_label = rep(c("CR", "PD"), each = 10L),
      stringsAsFactors = FALSE
    ),
    info = list(normalization = "TPM", scale = "log2")
  )
  ds_surv <- .make_surv_dataset(n = 40L, seed = 11L)

  sigs <- list(
    S = list(genes = c("G1", "G2", "G3"), derived_from = character(0))
  )

  expect_message(
    empty <- run_survival_benchmark(
      datasets = list(no_surv = ds_nosurv),
      signatures = sigs,
      method = "zscore",
      min_genes = 2L
    ),
    "Skipping cohort"
  )
  expect_equal(nrow(empty), 0L)
  expect_true(all(c("cohort", "signature", "hr", "cox_p") %in% names(empty)))

  surv_res <- run_survival_benchmark(
    datasets = list(synth_surv = ds_surv, no_surv = ds_nosurv),
    signatures = sigs,
    method = "zscore",
    min_genes = 2L
  )
  expect_true(nrow(surv_res) >= 1L)
  expect_true("synth_surv" %in% surv_res$cohort)

  auc_res <- run_benchmark(
    signatures = sigs,
    datasets = list(synth_surv = ds_surv),
    methods = "zscore",
    compute_pr = FALSE,
    compute_null = FALSE,
    min_n = 5L,
    min_per_class = 2L,
    min_genes = 2L
  )
  combo <- combine_response_and_survival(auc_res, surv_res)
  expect_true("auc" %in% names(combo))
  expect_true("hr" %in% names(combo))
  expect_true(any(is.finite(combo$auc) | is.na(combo$auc)))
  expect_true("cox_p" %in% names(combo))
  expect_true("c_index" %in% names(combo))
})
