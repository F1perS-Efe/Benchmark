#' Analytic AUC power via Hanley-McNeil (1982) variance
#'
#' Estimates the power to distinguish `target_auc` from `null_auc` using the
#' Hanley-McNeil analytic standard error and a normal approximation
#' (`stats::qnorm` / `stats::pnorm`). Deterministic; no external packages.
#'
#' Formulas (theta = target_auc):
#' - Q1 = theta / (2 - theta)
#' - Q2 = 2 * theta^2 / (1 + theta)
#' - SE(theta) = sqrt( (theta(1-theta) + (n_pos-1)(Q1 - theta^2) +
#'   (n_neg-1)(Q2 - theta^2)) / (n_pos * n_neg) )
#'
#' @param n_pos Number of positive (responder) samples.
#' @param n_neg Number of negative (non-responder) samples.
#' @param target_auc Alternative AUC to detect. Default `0.70`.
#' @param null_auc Null AUC. Default `0.50`.
#' @param alpha Type I error rate. Default `0.05`.
#' @param alternative `"two.sided"` (default) or `"greater"` / `"less"`.
#'
#' @return Named list with `power` (numeric in `[0, 1]`, or `NA_real_` if
#'   `n_pos < 2` or `n_neg < 2`) and `se_target`.
#' @export
#' @references
#' Hanley JA, McNeil BJ (1982). The meaning and use of the area under a
#' receiver operating characteristic (ROC) curve. Radiology.
estimate_auc_power <- function(n_pos,
                               n_neg,
                               target_auc = 0.70,
                               null_auc = 0.50,
                               alpha = 0.05,
                               alternative = c("two.sided", "greater", "less")) {
  alternative <- match.arg(alternative)
  n_pos <- as.integer(n_pos)[1]
  n_neg <- as.integer(n_neg)[1]
  target_auc <- as.numeric(target_auc)[1]
  null_auc <- as.numeric(null_auc)[1]
  alpha <- as.numeric(alpha)[1]

  if (!is.finite(n_pos) || !is.finite(n_neg) || n_pos < 2L || n_neg < 2L) {
    warning(
      "estimate_auc_power requires n_pos >= 2 and n_neg >= 2; returning NA.",
      call. = FALSE
    )
    return(list(power = NA_real_, se_target = NA_real_))
  }
  if (!is.finite(target_auc) || target_auc <= 0 || target_auc >= 1) {
    stop("`target_auc` must be in (0, 1).", call. = FALSE)
  }
  if (!is.finite(null_auc) || null_auc < 0 || null_auc > 1) {
    stop("`null_auc` must be in [0, 1].", call. = FALSE)
  }
  if (!is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be in (0, 1).", call. = FALSE)
  }

  hm_se <- function(theta, n1, n2) {
    q1 <- theta / (2 - theta)
    q2 <- (2 * theta * theta) / (1 + theta)
    se2 <- (
      theta * (1 - theta) +
        (n1 - 1L) * (q1 - theta * theta) +
        (n2 - 1L) * (q2 - theta * theta)
    ) / (as.numeric(n1) * as.numeric(n2))
    sqrt(max(0, se2))
  }

  se_target <- hm_se(target_auc, n_pos, n_neg)
  se_null <- hm_se(null_auc, n_pos, n_neg)

  if (!is.finite(se_target) || se_target <= 0 || !is.finite(se_null) || se_null <= 0) {
    return(list(power = NA_real_, se_target = se_target, se_null = se_null))
  }

  if (identical(alternative, "two.sided")) {
    z_crit <- stats::qnorm(1 - alpha / 2)
    c_up <- null_auc + z_crit * se_null
    c_lo <- null_auc - z_crit * se_null
    power <- (1 - stats::pnorm((c_up - target_auc) / se_target)) +
      stats::pnorm((c_lo - target_auc) / se_target)
  } else if (identical(alternative, "greater")) {
    z_crit <- stats::qnorm(1 - alpha)
    c_up <- null_auc + z_crit * se_null
    power <- 1 - stats::pnorm((c_up - target_auc) / se_target)
  } else {
    z_crit <- stats::qnorm(1 - alpha)
    c_lo <- null_auc - z_crit * se_null
    power <- stats::pnorm((c_lo - target_auc) / se_target)
  }

  power <- max(0, min(1, as.numeric(power)))
  list(power = power, se_target = as.numeric(se_target), se_null = as.numeric(se_null))
}
