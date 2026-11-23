#' Weighted AUC summary across cohorts (n, IVW, and logit-AUC REML)
#'
#' Computes cross-cohort summary statistics per `signature` x `method` (and
#' optional `stratify_by` subgroups):
#' - **n-weighted mean**: `sum(n * auc) / sum(n)`
#' - **inverse-variance weighted (IVW) mean**: `sum(w * auc) / sum(w)`, `w = 1 / auc_var`
#' - Heterogeneity: Cochran `Q` and `I^2` (IVW weights)
#' - **Logit-AUC REML + truncated HKSJ 95% CI**: on cohorts with open-interval
#'   AUC in `(0, 1)` and positive variance, estimates REML `tau2` on the
#'   `stats::qlogis` scale, Hartung-Knapp-Sidik-Jonkman standard error truncated
#'   at `max(1, q_hk)` with `stats::qt` critical values, and back-transforms
#'   estimates via `stats::plogis`; also reports logit-scale `i2_logit`.
#'
#' Default inclusion rules:
#' - `independent_only = TRUE`: includes only `independent == TRUE` (`NA` / `FALSE` excluded)
#' - `exclude_low_coverage = TRUE`: excludes rows with `coverage_status != "ok"`
#'   (if present) or `low_coverage == TRUE`; reports `n_excluded_coverage`
#' - `var_zero` / `auc_var <= 0`: excluded from IVW and REML only (retained in n-weighted mean)
#' - `skipped` or non-finite AUC rows are ignored
#'
#' With `stratify_by = NULL` (default), output preserves the unstratified schema.
#' Mixing distinct non-NA `endpoint_type` or `treatment` values within a summary
#' group emits a `warning()` and recommends `stratify_by`.
#'
#' Logit-AUC REML + HKSJ summaries are exploratory when cohort count is small
#' (`k <= 4`). Use `stratify_by` when cohorts differ in clinical endpoint or
#' treatment regimen.
#'
#' @param results Output `data.frame` from [run_benchmark()].
#' @param method Optional scoring method filter; `NULL` summarizes each method separately.
#' @param include_nonindependent If `TRUE`, retain non-independent / `NA` rows in summary.
#' @param independent_only If `TRUE` (default), retain only `independent == TRUE`.
#'   Overridden when `include_nonindependent = TRUE`.
#' @param exclude_var_zero If `TRUE` (default), exclude zero/NA variance rows from IVW and REML.
#' @param exclude_low_coverage If `TRUE` (default), exclude low-coverage rows from summary.
#' @param stratify_by Optional character vector of stratification columns (e.g.
#'   `c("treatment", "endpoint_type")`).
#'
#' @return A summary `data.frame` (including `auc_re`, `auc_re_ci_low`,
#'   `auc_re_ci_high`, `tau2_reml`, `i2_logit`).
#' @export
summarize_auc <- function(results,
                          method = NULL,
                          include_nonindependent = FALSE,
                          independent_only = TRUE,
                          exclude_var_zero = TRUE,
                          exclude_low_coverage = TRUE,
                          stratify_by = NULL) {
  if (!is.data.frame(results) || nrow(results) == 0L) {
    stop("`results` cannot be empty.", call. = FALSE)
  }
  req <- c("signature", "method", "dataset", "auc", "n", "skipped")
  miss <- setdiff(req, names(results))
  if (length(miss) > 0L) {
    stop("`results` is missing required columns: ", paste(miss, collapse = ", "), call. = FALSE)
  }

  if (!is.null(stratify_by)) {
    stratify_by <- as.character(stratify_by)
    stratify_by <- stratify_by[nzchar(stratify_by)]
    if (length(stratify_by) == 0L) {
      stratify_by <- NULL
    } else {
      miss_s <- setdiff(stratify_by, names(results))
      if (length(miss_s) > 0L) {
        stop(
          "`stratify_by` columns missing from results: ",
          paste(miss_s, collapse = ", "),
          call. = FALSE
        )
      }
    }
  }

  # include_nonindependent = TRUE overrides independent_only
  if (isTRUE(include_nonindependent)) {
    independent_only <- FALSE
  }

  df <- results
  if ("treatment" %in% names(df)) {
    df$treatment <- .normalize_treatment_label(df$treatment)
  }
  if (!is.null(method)) {
    df <- df[df$method == method, , drop = FALSE]
  }
  df <- df[!(df$skipped %in% TRUE) & is.finite(df$auc) & is.finite(df$n) & df$n > 0, , drop = FALSE]
  if (nrow(df) == 0L) {
    return(.summarize_auc_empty(stratify_by))
  }

  if (!"independent" %in% names(df)) {
    df$independent <- TRUE
  }
  if (!"var_zero" %in% names(df)) {
    df$var_zero <- if ("auc_var" %in% names(df)) {
      is.na(df$auc_var) | df$auc_var == 0
    } else {
      FALSE
    }
  }
  if (!"auc_var" %in% names(df)) {
    df$auc_var <- NA_real_
  }
  if (!"small_cohort" %in% names(df)) {
    df$small_cohort <- FALSE
  }
  if (!"low_coverage" %in% names(df)) {
    df$low_coverage <- FALSE
  }

  if (is.null(stratify_by)) {
    split_key <- paste(df$signature, df$method, sep = "\r")
  } else {
    strat_parts <- do.call(
      paste,
      c(lapply(stratify_by, function(col) {
        v <- as.character(df[[col]])
        v[is.na(v)] <- "<NA>"
        v
      }), list(sep = "\r"))
    )
    split_key <- paste(df$signature, df$method, strat_parts, sep = "\r")
  }
  groups <- split(df, split_key)

  rows <- lapply(groups, function(g) {
    row <- .summarize_auc_one_group(
      g,
      independent_only = independent_only,
      exclude_var_zero = exclude_var_zero,
      exclude_low_coverage = exclude_low_coverage,
      warn_clinical_mix = is.null(stratify_by)
    )
    if (!is.null(stratify_by)) {
      for (col in stratify_by) {
        row[[col]] <- g[[col]][1]
      }
    }
    row
  })

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' @noRd
.summarize_auc_one_group <- function(g,
                                     independent_only,
                                     exclude_var_zero,
                                     exclude_low_coverage,
                                     warn_clinical_mix = FALSE) {
  n_total_rows <- nrow(g)
  indep_flag <- as.logical(g$independent)
  n_small <- sum(as.logical(g$small_cohort), na.rm = TRUE)

  use_n <- g
  n_nonindep <- 0L
  if (isTRUE(independent_only)) {
    # Fail-closed: NA independent excluded from independent summary
    keep_indep <- !is.na(indep_flag) & indep_flag
    n_nonindep <- sum(!keep_indep)
    use_n <- use_n[keep_indep, , drop = FALSE]
  }

  n_low_cov <- 0L
  if (isTRUE(exclude_low_coverage) && nrow(use_n) > 0L) {
    keep_cov <- rep(TRUE, nrow(use_n))
    if ("coverage_status" %in% names(use_n)) {
      st <- as.character(use_n$coverage_status)
      keep_cov <- keep_cov & (is.na(st) | st == "ok")
    }
    if ("low_coverage" %in% names(use_n)) {
      low_flag <- as.logical(use_n$low_coverage)
      keep_cov <- keep_cov & (is.na(low_flag) | !low_flag)
    }
    n_low_cov <- sum(!keep_cov, na.rm = TRUE)
    use_n <- use_n[keep_cov, , drop = FALSE]
  }

  notes <- character(0)
  if (isTRUE(independent_only) && n_nonindep > 0L) {
    notes <- c(notes, sprintf("%d non-independent/NA cohorts excluded", n_nonindep))
  }
  if (isTRUE(exclude_low_coverage) && n_low_cov > 0L) {
    notes <- c(notes, sprintf("%d low-coverage cohorts excluded", n_low_cov))
  }
  if (n_small > 0L) {
    notes <- c(notes, sprintf("%d small-cohort warning(s)", n_small))
  }

  if (isTRUE(warn_clinical_mix) && nrow(use_n) > 0L) {
    mix_note <- .warn_clinical_heterogeneity(use_n, g$signature[1], g$method[1])
    if (!is.null(mix_note)) {
      notes <- c(notes, mix_note)
    }
  }

  k <- nrow(use_n)
  if (k == 0L) {
    return(data.frame(
      signature = g$signature[1],
      method = g$method[1],
      n_datasets = 0L,
      n_cohorts = 0L,
      n_cohorts_ivw = 0L,
      n_cohorts_re = 0L,
      n_total_rows = n_total_rows,
      auc_n_weighted = NA_real_,
      auc_ivw = NA_real_,
      Q = NA_real_,
      I2 = NA_real_,
      auc_re = NA_real_,
      auc_re_ci_low = NA_real_,
      auc_re_ci_high = NA_real_,
      tau2_reml = NA_real_,
      i2_logit = NA_real_,
      n_excluded_nonindependent = as.integer(n_nonindep),
      n_excluded_var_zero = 0L,
      n_var_zero_excluded = 0L,
      n_low_coverage_excluded = as.integer(n_low_cov),
      n_excluded_coverage = as.integer(n_low_cov),
      n_small_cohort = as.integer(n_small),
      note = paste(notes, collapse = "; "),
      stringsAsFactors = FALSE
    ))
  }

  n_vec <- as.numeric(use_n$n)
  auc_vec <- as.numeric(use_n$auc)
  auc_n <- sum(n_vec * auc_vec, na.rm = TRUE) / sum(n_vec, na.rm = TRUE)

  # P0-6: IVW ayri alt kume (var_zero / auc_var<=0 yalniz IVW'den)
  use_ivw <- use_n[
    !is.na(use_n$auc_var) & use_n$auc_var > 0 & is.finite(use_n$auc_var),
    ,
    drop = FALSE
  ]
  n_var_zero <- as.integer(k - nrow(use_ivw))
  if (isTRUE(exclude_var_zero) && n_var_zero > 0L) {
    notes <- c(notes, sprintf("%d zero/NA-variance rows excluded from IVW", n_var_zero))
    warning(
      "Signature '", g$signature[1], "' / method '", g$method[1],
      "': AUCs with zero or NA variance were excluded from IVW summary (still included in n-weighted mean).",
      call. = FALSE
    )
  } else if (!isTRUE(exclude_var_zero)) {
    # Even if exclude_var_zero=FALSE, IVW still needs positive variance
    use_ivw <- use_n[
      !is.na(use_n$auc_var) & use_n$auc_var > 0 & is.finite(use_n$auc_var),
      ,
      drop = FALSE
    ]
    n_var_zero <- as.integer(k - nrow(use_ivw))
  }

  if (nrow(use_ivw) == 0L) {
    auc_ivw <- NA_real_
    Q <- NA_real_
    I2 <- NA_real_
    k_ivw <- 0L
  } else {
    var_vec <- as.numeric(use_ivw$auc_var)
    auc_ivw_vec <- as.numeric(use_ivw$auc)
    w <- 1 / var_vec
    auc_ivw <- sum(w * auc_ivw_vec) / sum(w)
    k_ivw <- nrow(use_ivw)
    if (k_ivw >= 2L) {
      Q <- sum(w * (auc_ivw_vec - auc_ivw)^2)
      df_q <- k_ivw - 1L
      I2 <- if (is.finite(Q) && Q > 0) {
        max(0, (Q - df_q) / Q) * 100
      } else if (is.finite(Q) && Q == 0) {
        0
      } else {
        NA_real_
      }
    } else {
      Q <- NA_real_
      I2 <- NA_real_
    }
  }

  # Logit-AUC REML + truncated HKSJ (open-interval AUC, positive variance)
  use_re <- use_ivw[
    use_ivw$auc > 0 & use_ivw$auc < 1,
    ,
    drop = FALSE
  ]
  n_cohorts_re <- as.integer(nrow(use_re))
  n_boundary_re <- as.integer(nrow(use_ivw) - n_cohorts_re)
  if (n_boundary_re > 0L) {
    notes <- c(notes, sprintf("%d boundary AUC (0 or 1) cohort(s) excluded from logit-REML", n_boundary_re))
    warning(
      "Signature '", g$signature[1], "' / method '", g$method[1],
      "': ", n_boundary_re, " cohort(s) with AUC == 0 or 1 excluded from logit-REML.",
      call. = FALSE
    )
  }
  re_meta <- .logit_reml_meta(
    auc = as.numeric(use_re$auc),
    auc_var = as.numeric(use_re$auc_var)
  )

  data.frame(
    signature = g$signature[1],
    method = g$method[1],
    n_datasets = as.integer(k),
    n_cohorts = as.integer(k),
    n_cohorts_ivw = as.integer(k_ivw),
    n_cohorts_re = n_cohorts_re,
    n_total_rows = n_total_rows,
    auc_n_weighted = auc_n,
    auc_ivw = auc_ivw,
    Q = Q,
    I2 = I2,
    auc_re = re_meta$auc_re,
    auc_re_ci_low = re_meta$auc_re_ci_low,
    auc_re_ci_high = re_meta$auc_re_ci_high,
    tau2_reml = re_meta$tau2_reml,
    i2_logit = re_meta$i2_logit,
    n_excluded_nonindependent = as.integer(n_nonindep),
    n_excluded_var_zero = as.integer(n_var_zero),
    n_var_zero_excluded = as.integer(n_var_zero),
    n_low_coverage_excluded = as.integer(n_low_cov),
    n_excluded_coverage = as.integer(n_low_cov),
    n_small_cohort = as.integer(n_small),
    note = if (length(notes) == 0L) NA_character_ else paste(notes, collapse = "; "),
    stringsAsFactors = FALSE
  )
}

#' Logit-AUC REML + truncated HKSJ meta-analysis (base R only)
#'
#' @param auc Numeric AUC vector (open interval expected by caller).
#' @param auc_var Positive AUC variances (same length).
#' @return Named list with auc_re, auc_re_ci_low, auc_re_ci_high, tau2_reml, i2_logit.
#' @noRd
.logit_reml_meta <- function(auc, auc_var) {
  y <- stats::qlogis(auc)
  v <- auc_var / (auc * (1 - auc))^2
  ok <- is.finite(y) & is.finite(v) & v > 0
  y <- y[ok]
  v <- v[ok]
  k <- length(y)

  if (k == 0L) {
    return(list(
      auc_re = NA_real_,
      auc_re_ci_low = NA_real_,
      auc_re_ci_high = NA_real_,
      tau2_reml = NA_real_,
      i2_logit = NA_real_
    ))
  }

  if (k == 1L) {
    return(list(
      auc_re = stats::plogis(y[1]),
      auc_re_ci_low = NA_real_,
      auc_re_ci_high = NA_real_,
      tau2_reml = NA_real_,
      i2_logit = NA_real_
    ))
  }

  reml_ll <- function(tau2) {
    w <- 1 / (v + tau2)
    sw <- sum(w)
    mu <- sum(w * y) / sw
    -0.5 * (sum(log(v + tau2)) + log(sw) + sum(w * (y - mu)^2))
  }

  upper <- max(20, 10 * stats::var(y))
  opt <- stats::optimize(reml_ll, interval = c(0, upper), maximum = TRUE)
  ll0 <- reml_ll(0)
  if (is.finite(ll0) && (!is.finite(opt$objective) || ll0 >= opt$objective)) {
    tau2_reml <- 0
  } else {
    tau2_reml <- max(0, as.numeric(opt$maximum))
  }

  w_re <- 1 / (v + tau2_reml)
  theta_re <- sum(w_re * y) / sum(w_re)

  q_hk <- sum(w_re * (y - theta_re)^2) / (k - 1)
  se_hk <- sqrt(max(1, q_hk) / sum(w_re))
  t_crit <- stats::qt(0.975, df = k - 1)

  w_fe <- 1 / v
  theta_fe <- sum(w_fe * y) / sum(w_fe)
  Q_logit <- sum(w_fe * (y - theta_fe)^2)
  i2_logit <- if (is.finite(Q_logit) && Q_logit > 0) {
    max(0, (Q_logit - (k - 1)) / Q_logit) * 100
  } else {
    0
  }

  list(
    auc_re = stats::plogis(theta_re),
    auc_re_ci_low = stats::plogis(theta_re - t_crit * se_hk),
    auc_re_ci_high = stats::plogis(theta_re + t_crit * se_hk),
    tau2_reml = tau2_reml,
    i2_logit = i2_logit
  )
}

#' @noRd
.warn_clinical_heterogeneity <- function(use_n, signature, method) {
  bits <- character(0)
  for (col in c("endpoint_type", "treatment")) {
    if (!col %in% names(use_n)) {
      next
    }
    vals <- unique(as.character(use_n[[col]]))
    vals <- vals[!is.na(vals) & nzchar(vals) & vals != "NA"]
    if (length(vals) > 1L) {
      bits <- c(
        bits,
        paste0(col, "=[", paste(vals, collapse = ", "), "]")
      )
    }
  }
  if (length(bits) > 0L) {
    warning(
      "Imza '", signature, "' / method '", method,
      "': multiple distinct non-NA clinical subgroup values are pooled (",
      paste(bits, collapse = "; "),
      "). Prefer stratify_by = c(\"treatment\", \"endpoint_type\").",
      call. = FALSE
    )
    return(paste0("Pooled across distinct subgroups (", paste(bits, collapse = "; "), "); stratify_by recommended"))
  }
  NULL
}

#' Normalize treatment regimen labels across cohorts
#'
#' @noRd
.normalize_treatment_label <- function(x) {
  x <- as.character(x)
  out <- trimws(x)
  low <- tolower(gsub("[[:space:]_]+", "-", out))
  low <- gsub("-+", "-", low)
  out[!is.na(low) & low %in% c("anti-pd-1", "anti-pd1", "pd-1", "pd1")] <- "anti-PD1"
  out[!is.na(low) & low %in% c("anti-pd-l1", "anti-pdl1", "pd-l1", "pdl1")] <- "anti-PDL1"
  out[!is.na(low) & low %in% c("anti-ctla-4", "anti-ctla4", "ctla-4", "ctla4")] <- "anti-CTLA4"
  out[!is.na(low) & low %in% c("combo", "combination", "anti-pd1+anti-ctla4", "anti-pd-1+anti-ctla-4", "ipi+nivo")] <- "combo"
  out
}

.summarize_auc_empty <- function(stratify_by = NULL) {
  out <- data.frame(
    signature = character(),
    method = character(),
    n_datasets = integer(),
    n_cohorts = integer(),
    n_cohorts_ivw = integer(),
    n_cohorts_re = integer(),
    n_total_rows = integer(),
    auc_n_weighted = numeric(),
    auc_ivw = numeric(),
    Q = numeric(),
    I2 = numeric(),
    auc_re = numeric(),
    auc_re_ci_low = numeric(),
    auc_re_ci_high = numeric(),
    tau2_reml = numeric(),
    i2_logit = numeric(),
    n_excluded_nonindependent = integer(),
    n_excluded_var_zero = integer(),
    n_var_zero_excluded = integer(),
    n_low_coverage_excluded = integer(),
    n_excluded_coverage = integer(),
    n_small_cohort = integer(),
    note = character(),
    stringsAsFactors = FALSE
  )
  if (!is.null(stratify_by)) {
    for (col in stratify_by) {
      out[[col]] <- character()
    }
  }
  out
}

#' Format meta-analysis summary table for display (Shiny & HTML report)
#'
#' Pairs logit-scale REML AUC (95% CI) with logit-scale heterogeneity
#' metrics (`i2_logit` and `tau2_reml`) so heterogeneity and point estimates
#' share the same scale.
#'
#' @param sm [summarize_auc()] output `data.frame`.
#' @param stratum_col Optional single stratification column name (or character
#'   vector of stratification column names).
#' @return Formatted `data.frame` for display.
#' @export
format_meta_summary <- function(sm, stratum_col = NULL) {
  if (!is.data.frame(sm) || nrow(sm) == 0L) {
    return(sm)
  }
  em <- "\u2014"
  en <- "\u2013"
  fmt3 <- function(x) {
    x <- as.numeric(x)
    ifelse(is.finite(x), sprintf("%.3f", x), em)
  }
  fmt_i2 <- function(x) {
    x <- as.numeric(x)
    ifelse(is.finite(x), sprintf("%.1f", x), em)
  }
  fmt_auc_ci <- function(est, lo, hi) {
    est <- as.numeric(est)
    lo <- as.numeric(lo)
    hi <- as.numeric(hi)
    out <- rep(em, length(est))
    full <- is.finite(est) & is.finite(lo) & is.finite(hi)
    out[full] <- sprintf("%.3f (%.3f%s%.3f)", est[full], lo[full], en, hi[full])
    point_only <- is.finite(est) & !full
    out[point_only] <- sprintf("%.3f", est[point_only])
    out
  }
  translate_note <- function(x) {
    x <- as.character(x)
    miss <- is.na(x) | !nzchar(x)
    x[miss] <- em
    x <- gsub(
      "([0-9]+) bagimsiz-olmayan/NA satir ozetten cikarildi",
      "\\1 non-independent/NA cohorts excluded",
      x
    )
    x <- gsub(
      "([0-9]+) dusuk kapsama satir ozetten cikarildi",
      "\\1 low-coverage cohorts excluded",
      x
    )
    x <- gsub(
      "([0-9]+) kucuk kohort uyarisi",
      "\\1 small-cohort warning(s)",
      x
    )
    x <- gsub(
      "([0-9]+) var_zero/NA varyans satir IVW'den cikarildi",
      "\\1 zero/NA-variance rows excluded from IVW",
      x
    )
    x
  }

  k_col <- if ("n_cohorts" %in% names(sm)) sm$n_cohorts else sm$n_datasets
  i2_val <- if ("i2_logit" %in% names(sm)) sm$i2_logit else rep(NA_real_, nrow(sm))
  tau2_val <- if ("tau2_reml" %in% names(sm)) sm$tau2_reml else rep(NA_real_, nrow(sm))

  if (!is.null(stratum_col)) {
    stratum_cols <- intersect(as.character(stratum_col), names(sm))
    stratum <- if (length(stratum_cols) >= 2L) {
      do.call(paste, c(lapply(stratum_cols, function(c) as.character(sm[[c]])), list(sep = " / ")))
    } else if (length(stratum_cols) == 1L) {
      as.character(sm[[stratum_cols[1]]])
    } else {
      rep(em, nrow(sm))
    }
    return(data.frame(
      Signature = as.character(sm$signature),
      Method = as.character(sm$method),
      Stratum = stratum,
      k = as.integer(k_col),
      `N-Weighted AUC` = fmt3(sm$auc_n_weighted),
      `REML AUC (95% CI)` = fmt_auc_ci(sm$auc_re, sm$auc_re_ci_low, sm$auc_re_ci_high),
      `I2 logit (%)` = fmt_i2(i2_val),
      `tau2 (REML)` = fmt3(tau2_val),
      check.names = FALSE,
      stringsAsFactors = FALSE
    ))
  }

  data.frame(
    Signature = as.character(sm$signature),
    Method = as.character(sm$method),
    k = as.integer(k_col),
    `N-Weighted AUC` = fmt3(sm$auc_n_weighted),
    `IVW AUC` = fmt3(sm$auc_ivw),
    `REML AUC (95% CI)` = fmt_auc_ci(sm$auc_re, sm$auc_re_ci_low, sm$auc_re_ci_high),
    `I2 logit (%)` = fmt_i2(i2_val),
    `tau2 (REML)` = fmt3(tau2_val),
    Note = if ("note" %in% names(sm)) translate_note(sm$note) else rep(em, nrow(sm)),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

.format_meta_summary <- format_meta_summary
