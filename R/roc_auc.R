#' Compute ROC-AUC with fixed direction (no silent sign flipping)
#'
#' Aligns scores via [orient_score()] so higher values represent responders,
#' then invokes `pROC::roc()` with `direction = "<"`. When AUC < 0.5, values
#' are reported directly without flipping.
#'
#' Computes DeLong 95% confidence intervals (`ci.auc`) and AUC variance (`var`).
#'
#' @param scores Numeric sample score vector (optionally named).
#' @param response Integer/numeric 0/1 response vector (aligned to `scores`).
#' @param direction Signature direction: `"<"` or `">"`. Passed to [orient_score()].
#' @param min_n Minimum total sample size for `small_cohort` flag. Default `20L`.
#' @param min_per_class Minimum per-class sample size for `small_cohort` flag. Default `5L`.
#'
#' @return A list with `auc`, `auc_ci_low`, `auc_ci_high`, `auc_var`,
#'   `var_zero`, `small_cohort`, `n`, `n_pos`, `n_neg`, `roc`,
#'   `scores_oriented`, and `response`.
#' @export
roc_auc <- function(scores,
                    response,
                    direction = c("<", ">"),
                    min_n = 20L,
                    min_per_class = 5L) {
  direction <- match.arg(direction)

  scores <- as.numeric(scores)
  response <- as.integer(response)

  if (length(scores) != length(response)) {
    stop("`scores` and `response` must have the same length.", call. = FALSE)
  }

  ok <- is.finite(scores) & !is.na(response)
  scores <- scores[ok]
  response <- response[ok]

  n_pos <- sum(response == 1L)
  n_neg <- sum(response == 0L)
  n <- length(response)

  if (n_pos == 0L || n_neg == 0L) {
    stop(
      "Both response classes are required for ROC (n_pos=", n_pos, ", n_neg=", n_neg, ").",
      call. = FALSE
    )
  }

  # Single-direction alignment; pROC always receives direction = "<"
  scores_oriented <- orient_score(scores, direction = direction)

  roc_obj <- pROC::roc(
    response = response,
    predictor = scores_oriented,
    levels = c(0, 1),
    direction = "<",
    quiet = TRUE
  )

  auc_val <- as.numeric(pROC::auc(roc_obj))
  # Suppress pROC warning when AUC == 1; flagged explicitly via var_zero
  ci <- suppressWarnings(as.numeric(pROC::ci.auc(roc_obj, method = "delong")))
  # ci.auc returns c(low, median/auc, high)
  auc_ci_low <- ci[1]
  auc_ci_high <- ci[3]

  auc_var <- suppressWarnings(as.numeric(pROC::var(roc_obj, method = "delong")))
  var_zero <- is.na(auc_var) || isTRUE(auc_var == 0)

  small_cohort <- (n < min_n) || (min(n_pos, n_neg) < min_per_class)

  list(
    auc = auc_val,
    auc_ci_low = auc_ci_low,
    auc_ci_high = auc_ci_high,
    auc_var = auc_var,
    var_zero = var_zero,
    small_cohort = small_cohort,
    n = n,
    n_pos = n_pos,
    n_neg = n_neg,
    roc = roc_obj,
    scores_oriented = scores_oriented,
    response = response
  )
}
