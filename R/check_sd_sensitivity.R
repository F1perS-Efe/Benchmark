#' Systematic SD-handling sensitivity for spectrum bias
#'
#' Re-runs the same signature x cohort benchmark under three response policies:
#' `"exclude"` (CR/PR vs PD), `"as_0"` (ORR: CR/PR vs SD/PD), and `"as_1"`
#' (DCR: CR/PR/SD vs PD). Comparing AUC shifts makes spectrum bias from SD
#' coding visible. The `sd_delta_threshold` default is arbitrary and is not a
#' clinical decision rule; do not tune it after inspecting results.
#'
#' @param signatures Named signature list for `run_benchmark()`.
#' @param datasets Named dataset list for `run_benchmark()`.
#' @param methods Scoring methods. Default includes all five methods.
#' @param direction Default signature direction.
#' @param min_genes Minimum genes for scoring.
#' @param coverage_threshold Coverage warning threshold (0-1).
#' @param apply_log2 Apply log2 when scale is linear.
#' @param min_n,min_per_class Small-cohort thresholds for `roc_auc()`.
#' @param sd_delta_threshold Absolute AUC delta threshold for `sd_sensitive`
#'   (default `0.05`; arbitrary).
#' @param seed Shared RNG seed for all three policy runs.
#'
#' @return One-row-per-cell `data.frame` with `n_*`, `auc_*`, `delta_orr`,
#'   `delta_dcr`, `max_abs_delta`, and `sd_sensitive`.
#' @export
check_sd_sensitivity <- function(signatures,
                                 datasets,
                                 methods = c("mean", "zscore", "median", "rank", "ss_rank"),
                                 direction = "<",
                                 min_genes = 3L,
                                 coverage_threshold = 0.70,
                                 apply_log2 = TRUE,
                                 min_n = 20L,
                                 min_per_class = 5L,
                                 sd_delta_threshold = 0.05,
                                 seed = 42L) {
  if (!is.numeric(sd_delta_threshold) || length(sd_delta_threshold) != 1L ||
      is.na(sd_delta_threshold) || sd_delta_threshold < 0) {
    stop("`sd_delta_threshold` must be a single non-negative number.", call. = FALSE)
  }

  run_one <- function(policy) {
    run_benchmark(
      signatures = signatures,
      datasets = datasets,
      methods = methods,
      direction = direction,
      min_genes = min_genes,
      coverage_threshold = coverage_threshold,
      apply_log2 = apply_log2,
      sd_handling = policy,
      min_n = min_n,
      min_per_class = min_per_class,
      compute_pr = FALSE,
      compute_null = FALSE,
      seed = seed
    )
  }

  res_ex <- run_one("exclude")
  res_0 <- run_one("as_0")
  res_1 <- run_one("as_1")

  pick <- function(df, tag) {
    out <- data.frame(
      dataset = as.character(df$dataset),
      signature = as.character(df$signature),
      method = as.character(df$method),
      n = as.integer(df$n),
      auc = as.numeric(df$auc),
      stringsAsFactors = FALSE
    )
    names(out)[names(out) == "n"] <- paste0("n_", tag)
    names(out)[names(out) == "auc"] <- paste0("auc_", tag)
    out
  }

  m <- pick(res_ex, "exclude")
  m <- merge(m, pick(res_0, "as_0"), by = c("dataset", "signature", "method"), all = TRUE)
  m <- merge(m, pick(res_1, "as_1"), by = c("dataset", "signature", "method"), all = TRUE)

  delta_orr <- m$auc_as_0 - m$auc_exclude
  delta_dcr <- m$auc_as_1 - m$auc_exclude
  abs_pair <- cbind(abs(delta_orr), abs(delta_dcr))
  max_abs_delta <- apply(abs_pair, 1L, function(x) {
    if (all(is.na(x))) {
      NA_real_
    } else {
      max(x, na.rm = TRUE)
    }
  })

  m$delta_orr <- delta_orr
  m$delta_dcr <- delta_dcr
  m$max_abs_delta <- as.numeric(max_abs_delta)
  m$sd_sensitive <- !is.na(m$max_abs_delta) & m$max_abs_delta > sd_delta_threshold

  m[order(m$dataset, m$signature, m$method), , drop = FALSE]
}
