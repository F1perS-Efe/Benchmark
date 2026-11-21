#' Veri setleri uzerinde agirlikli AUC ozeti (n, IVW ve logit-AUC REML)
#'
#' Her `signature` x `method` (ve istege bagli `stratify_by` alt gruplari) icin:
#' - **n-agirlikli** ortalama: `sum(n * auc) / sum(n)`
#' - **ters-varyans (IVW)** ortalama: `sum(w * auc) / sum(w)`, `w = 1 / auc_var`
#' - Heterojenlik: Cochran `Q` ve `I^2` (IVW agirliklariyla)
#' - **Logit-AUC REML + kesilmis HKSJ 95 pct GA**: AUC acik aralikta (0, 1) ve
#'   pozitif varyansli kohortlarda `stats::qlogis` olceginde REML `tau2`,
#'   Hartung-Knapp-Sidik-Jonkman `max(1, q_hk)` kesmeli standart hata ve
#'   `stats::qt` kritik degeri; nokta/GA `stats::plogis` ile geri donusur.
#'   Ayrica logit olcekli `i2_logit` raporlanir.
#'
#' Varsayilan kurallar:
#' - `independent_only = TRUE`: yalniz `independent == TRUE` (NA / FALSE dislanir)
#' - `exclude_low_coverage = TRUE`: `coverage_status != "ok"` (varsa) veya
#'   `low_coverage == TRUE` satirlar dislanir; `n_excluded_coverage` raporlanir
#' - `var_zero` / `auc_var <= 0` yalnizca IVW / REML'den cikarilir; n-agirlikli ortalamaya girer
#' - `skipped` veya sonlu olmayan AUC satirlar yok sayilir
#'
#' `stratify_by = NULL` (varsayilan) iken cikti geriye donuk uyumludur. Ayni ozet
#' grubunda birden fazla farkli non-NA `endpoint_type` veya `treatment` degeri
#' birlestiriliyorsa ASCII `warning()` uretilir ve `stratify_by` onerilir.
#' Varsayilan HTML rapor, tabakasiz ozetten hemen sonra
#' `stratify_by = "treatment"` ve `stratify_by = "endpoint_type"` tablolarini
#' ayri ayri render eder; forest grafigi tabakasiz ozeti korur.
#'
#' Logit-AUC REML + HKSJ ozeti kucuk kohort sayisinda (`k <= 4`) **kesifseldir**.
#' Farkli klinik sonlanimli kohortlarda `stratify_by` kullanin; havuzlama
#' klinik olarak anlamsiz ortalamalar uretebilir.
#'
#' @param results `run_benchmark()` cikti tablosu.
#' @param method Belirli skorlama yontemi; `NULL` ise tum yontemler ayri ozetlenir.
#' @param include_nonindependent `TRUE` ise bagimsiz-olmayan / NA satirlar ozete girer.
#' @param independent_only `TRUE` (varsayilan) ise yalniz `independent == TRUE`.
#'   `include_nonindependent = TRUE` verilirse bu bayrak gecersiz kilinir.
#' @param exclude_var_zero `TRUE` (varsayilan) ise varyansi 0/NA olanlar IVW'den cikarilir.
#' @param exclude_low_coverage `TRUE` (varsayilan) ise dusuk kapsama satirlari ozetten cikarilir.
#' @param stratify_by Opsiyonel karakter vektoru (orn. `c("treatment", "endpoint_type")`);
#'   verildiginde ozet her alt grup icinde ayri hesaplanir.
#'
#' @return `data.frame` ozet satirlari (`auc_re`, `auc_re_ci_low`, `auc_re_ci_high`,
#'   `tau2_reml`, `i2_logit` dahil).
#' @export
summarize_auc <- function(results,
                          method = NULL,
                          include_nonindependent = FALSE,
                          independent_only = TRUE,
                          exclude_var_zero = TRUE,
                          exclude_low_coverage = TRUE,
                          stratify_by = NULL) {
  if (!is.data.frame(results) || nrow(results) == 0L) {
    stop("`results` bos olamaz.", call. = FALSE)
  }
  req <- c("signature", "method", "dataset", "auc", "n", "skipped")
  miss <- setdiff(req, names(results))
  if (length(miss) > 0L) {
    stop("`results` eksik sutunlar: ", paste(miss, collapse = ", "), call. = FALSE)
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
    if ("coverage_status" %in% names(use_n)) {
      st <- as.character(use_n$coverage_status)
      keep_cov <- is.na(st) | st == "ok"
      n_low_cov <- sum(!keep_cov, na.rm = TRUE)
      use_n <- use_n[keep_cov, , drop = FALSE]
    } else if ("low_coverage" %in% names(use_n)) {
      low_flag <- as.logical(use_n$low_coverage)
      keep_cov <- is.na(low_flag) | !low_flag
      n_low_cov <- sum(!keep_cov, na.rm = TRUE)
      use_n <- use_n[keep_cov, , drop = FALSE]
    }
  }

  notes <- character(0)
  if (isTRUE(independent_only) && n_nonindep > 0L) {
    notes <- c(notes, sprintf("%d bagimsiz-olmayan/NA satir ozetten cikarildi", n_nonindep))
  }
  if (isTRUE(exclude_low_coverage) && n_low_cov > 0L) {
    notes <- c(notes, sprintf("%d dusuk kapsama satir ozetten cikarildi", n_low_cov))
  }
  if (n_small > 0L) {
    notes <- c(notes, sprintf("%d kucuk kohort uyarisi", n_small))
  }

  if (isTRUE(warn_clinical_mix) && nrow(use_n) > 0L) {
    .warn_clinical_heterogeneity(use_n, g$signature[1], g$method[1])
  }

  k <- nrow(use_n)
  if (k == 0L) {
    return(data.frame(
      signature = g$signature[1],
      method = g$method[1],
      n_datasets = 0L,
      n_cohorts = 0L,
      n_cohorts_ivw = 0L,
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
    notes <- c(notes, sprintf("%d var_zero/NA varyans satir IVW'den cikarildi", n_var_zero))
    warning(
      "Imza '", g$signature[1], "' / method '", g$method[1],
      "': varyansi 0 veya NA olan AUC'ler IVW ozetinden cikarildi (n-agirlikli ortalamaya dahil).",
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
  }
}

.summarize_auc_empty <- function(stratify_by = NULL) {
  out <- data.frame(
    signature = character(),
    method = character(),
    n_datasets = integer(),
    n_cohorts = integer(),
    n_cohorts_ivw = integer(),
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
