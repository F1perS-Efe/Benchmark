#' Empty survival result schema (0-row data.frame)
#'
#' @param with_signature Include a `signature` column.
#' @return 0-row `data.frame` with survival columns.
#' @noRd
.survival_empty_schema <- function(with_signature = FALSE) {
  out <- data.frame(
    cohort = character(),
    endpoint = character(),
    n_samples = integer(),
    n_events = integer(),
    hr = numeric(),
    hr_ci_low = numeric(),
    hr_ci_high = numeric(),
    cox_p = numeric(),
    c_index = numeric(),
    logrank_p = numeric(),
    median_surv_high = numeric(),
    median_surv_low = numeric(),
    ph_p = numeric(),
    ph_violated = logical(),
    surv_note = character(),
    small_sample = logical(),
    stringsAsFactors = FALSE
  )
  if (isTRUE(with_signature)) {
    out$signature <- character()
    out$method <- character()
    out$p_adj_bh <- numeric()
    out$logrank_p_adj_bh <- numeric()
    out <- out[, c(
      "cohort", "signature", "method", "endpoint", "n_samples", "n_events",
      "hr", "hr_ci_low", "hr_ci_high", "cox_p", "p_adj_bh", "c_index",
      "logrank_p", "logrank_p_adj_bh", "median_surv_high", "median_surv_low",
      "ph_p", "ph_violated", "surv_note", "small_sample"
    ), drop = FALSE]
  }
  out
}

#' Align a score vector to sample_meta rows
#'
#' @param score Named or unnamed numeric vector.
#' @param sample_meta Data frame with `sample_id`.
#' @return Numeric vector of length `nrow(sample_meta)` (may contain NA).
#' @noRd
.align_score_to_meta <- function(score, sample_meta) {
  if (!is.numeric(score)) {
    stop("`score` must be a numeric vector.", call. = FALSE)
  }
  sid <- as.character(sample_meta$sample_id)
  n <- length(sid)
  if (!is.null(names(score)) && any(nzchar(names(score)))) {
    out <- as.numeric(score[sid])
    names(out) <- sid
    return(out)
  }
  if (length(score) != n) {
    stop(
      "`score` length (", length(score),
      ") must equal nrow(sample_meta) (", n,
      ") when unnamed, or be named by sample_id.",
      call. = FALSE
    )
  }
  out <- as.numeric(score)
  names(out) <- sid
  out
}

#' Derive endpoint label from a time column name
#'
#' @param time_col Time column name (e.g. `"os_time"`).
#' @return Uppercase endpoint tag (e.g. `"OS"`).
#' @noRd
.survival_endpoint_label <- function(time_col) {
  tc <- as.character(time_col)[1]
  ep <- sub("_time$", "", tc, ignore.case = TRUE)
  toupper(ep)
}

#' Extract median survival for High/Low strata from survfit
#'
#' @param sf A `survfit` object stratified by High/Low.
#' @param grp Character group vector used in the fit.
#' @return Named list with `high` and `low` medians (`NA_real_` if not reached).
#' @noRd
.survfit_group_medians <- function(sf, grp) {
  med_high <- NA_real_
  med_low <- NA_real_
  tbl <- tryCatch(summary(sf)$table, error = function(e) NULL)
  if (is.null(tbl)) {
    return(list(high = med_high, low = med_low))
  }
  # Single stratum -> named numeric; multi -> matrix
  if (is.null(dim(tbl))) {
    med <- as.numeric(tbl[["median"]])
    g_unique <- unique(as.character(grp))
    if (length(g_unique) == 1L && identical(g_unique, "High")) {
      med_high <- med
    } else if (length(g_unique) == 1L && identical(g_unique, "Low")) {
      med_low <- med
    }
    return(list(high = med_high, low = med_low))
  }
  rn <- rownames(tbl)
  if (is.null(rn)) {
    return(list(high = med_high, low = med_low))
  }
  # survfit rownames look like "grp=High" or "High"
  for (i in seq_along(rn)) {
    lab <- rn[i]
    med_i <- as.numeric(tbl[i, "median"])
    if (grepl("High", lab, fixed = TRUE)) {
      med_high <- med_i
    } else if (grepl("Low", lab, fixed = TRUE)) {
      med_low <- med_i
    }
  }
  list(high = med_high, low = med_low)
}

#' Evaluate signature scores against a survival endpoint (Cox PH + KM)
#'
#' Requires the Recommended package `survival` (`Suggests`). Aligns `score`
#' to `dataset$sample_meta` by `sample_id` (or by position when lengths match).
#' Incomplete rows (`NA` in time, event, or score) are dropped. Small-sample
#' or zero-variance scores return `NA` metrics with `small_sample = TRUE`
#' without throwing.
#'
#' @param dataset Validated SigBenchR dataset list.
#' @param score Numeric score vector (named by `sample_id` preferred).
#' @param time_col Survival time column in `sample_meta`.
#' @param event_col Event indicator column (`0` = censored, `1` = event).
#' @param cutoff Grouping rule for Kaplan-Meier / log-rank. Currently only
#'   `"median"` is supported.
#' @param orient_direction `"higher_better"` keeps the scaled score;
#'   `"lower_better"` negates it before Cox PH.
#'
#' @return One-row `data.frame` with Cox and KM summary columns.
#' @export
evaluate_survival <- function(dataset,
                              score,
                              time_col = "os_time",
                              event_col = "os_event",
                              cutoff = "median",
                              orient_direction = c("higher_better", "lower_better")) {
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop(
      "Survival analysis requires the 'survival' package. ",
      "Install it with: install.packages(\"survival\")",
      call. = FALSE
    )
  }
  ds <- validate_dataset(dataset)
  orient_direction <- match.arg(orient_direction)
  cutoff <- match.arg(cutoff, choices = "median")

  time_col <- as.character(time_col)[1]
  event_col <- as.character(event_col)[1]
  if (!time_col %in% names(ds$sample_meta)) {
    stop("`sample_meta` missing time column '", time_col, "'.", call. = FALSE)
  }
  if (!event_col %in% names(ds$sample_meta)) {
    stop("`sample_meta` missing event column '", event_col, "'.", call. = FALSE)
  }

  cohort <- if (!is.null(ds$info$dataset_name) &&
                length(ds$info$dataset_name) > 0L &&
                !is.na(ds$info$dataset_name[1]) &&
                nzchar(as.character(ds$info$dataset_name)[1])) {
    as.character(ds$info$dataset_name)[1]
  } else {
    "unknown"
  }
  endpoint <- .survival_endpoint_label(time_col)

  score_vec <- .align_score_to_meta(score, ds$sample_meta)
  time <- as.numeric(ds$sample_meta[[time_col]])
  event <- as.numeric(ds$sample_meta[[event_col]])

  neg_time <- is.finite(time) & time < 0
  if (any(neg_time)) {
    warning(
      sum(neg_time), " sample(s) with negative survival time excluded.",
      call. = FALSE
    )
  }
  ok <- is.finite(time) & time >= 0 & !is.na(event) & is.finite(score_vec)
  time <- time[ok]
  event <- event[ok]
  score_use <- score_vec[ok]

  n_samples <- length(score_use)
  n_events <- as.integer(sum(event == 1, na.rm = TRUE))

  na_row <- function(small) {
    data.frame(
      cohort = cohort,
      endpoint = endpoint,
      n_samples = as.integer(n_samples),
      n_events = as.integer(n_events),
      hr = NA_real_,
      hr_ci_low = NA_real_,
      hr_ci_high = NA_real_,
      cox_p = NA_real_,
      c_index = NA_real_,
      logrank_p = NA_real_,
      median_surv_high = NA_real_,
      median_surv_low = NA_real_,
      ph_p = NA_real_,
      ph_violated = NA,
      surv_note = "Insufficient sample size or events",
      small_sample = isTRUE(small),
      stringsAsFactors = FALSE
    )
  }

  if (n_samples < 10L || n_events < 3L ||
      !is.finite(stats::sd(score_use)) || stats::sd(score_use) == 0) {
    return(na_row(TRUE))
  }

  z_score <- as.numeric(scale(score_use))
  if (identical(orient_direction, "lower_better")) {
    z_score <- -z_score
  }

  hr <- NA_real_
  hr_ci_low <- NA_real_
  hr_ci_high <- NA_real_
  cox_p <- NA_real_
  c_index <- NA_real_
  ph_p <- NA_real_
  ph_violated <- NA

  cox_out <- tryCatch(
    {
      fit <- survival::coxph(survival::Surv(time, event) ~ z_score)
      s <- summary(fit)
      zph <- tryCatch(suppressWarnings(survival::cox.zph(fit)), error = function(e) NULL)
      ph_val <- if (!is.null(zph) && !is.null(zph$table) && nrow(zph$table) >= 1L) {
        as.numeric(zph$table[1L, "p"])
      } else {
        NA_real_
      }
      list(
        hr = as.numeric(exp(stats::coef(fit)))[1],
        ci = as.numeric(exp(stats::confint(fit)))[1:2],
        cox_p = as.numeric(s$coefficients[, "Pr(>|z|)"])[1],
        c_index = as.numeric(s$concordance[1]),
        ph_p = ph_val
      )
    },
    error = function(e) NULL
  )
  if (!is.null(cox_out)) {
    hr <- cox_out$hr
    hr_ci_low <- cox_out$ci[1]
    hr_ci_high <- cox_out$ci[2]
    cox_p <- cox_out$cox_p
    c_index <- cox_out$c_index
    ph_p <- cox_out$ph_p
    ph_violated <- is.finite(ph_p) && ph_p < 0.05
  }

  logrank_p <- NA_real_
  median_surv_high <- NA_real_
  median_surv_low <- NA_real_
  surv_note <- NA_character_

  med_cut <- stats::median(score_use)
  grp <- ifelse(score_use >= med_cut, "High", "Low")
  n_high <- sum(grp == "High", na.rm = TRUE)
  n_low <- sum(grp == "Low", na.rm = TRUE)
  if (min(n_high, n_low) < 3L) {
    grp_strict <- ifelse(score_use > med_cut, "High", "Low")
    if (min(sum(grp_strict == "High"), sum(grp_strict == "Low")) >= 3L) {
      grp <- grp_strict
      n_high <- sum(grp == "High", na.rm = TRUE)
      n_low <- sum(grp == "Low", na.rm = TRUE)
    }
  }

  if (n_high >= 3L && n_low >= 3L) {
    lr <- tryCatch(
      {
        sd_fit <- survival::survdiff(survival::Surv(time, event) ~ grp)
        stats::pchisq(sd_fit$chisq, df = 1, lower.tail = FALSE)
      },
      error = function(e) NA_real_
    )
    logrank_p <- as.numeric(lr)

    meds <- tryCatch(
      {
        sf <- survival::survfit(survival::Surv(time, event) ~ grp)
        .survfit_group_medians(sf, grp)
      },
      error = function(e) list(high = NA_real_, low = NA_real_)
    )
    median_surv_high <- meds$high
    median_surv_low <- meds$low
  } else {
    surv_note <- "Too few samples in one median stratum (< 3) due to ties"
  }

  data.frame(
    cohort = cohort,
    endpoint = endpoint,
    n_samples = as.integer(n_samples),
    n_events = as.integer(n_events),
    hr = hr,
    hr_ci_low = hr_ci_low,
    hr_ci_high = hr_ci_high,
    cox_p = cox_p,
    c_index = c_index,
    logrank_p = logrank_p,
    median_surv_high = median_surv_high,
    median_surv_low = median_surv_low,
    ph_p = ph_p,
    ph_violated = ph_violated,
    surv_note = surv_note,
    small_sample = FALSE,
    stringsAsFactors = FALSE
  )
}

#' Run survival benchmark across cohorts and signatures
#'
#' Scores each signature with `score_signature()` on cohorts that contain
#' both `time_col` and `event_col`. Cohorts lacking survival columns are
#' skipped with an informative `message()`. When no suitable cohort exists,
#' returns a 0-row schema-compatible `data.frame`.
#'
#' @param datasets Named list of SigBenchR datasets.
#' @param signatures Named signature list (same contract as `run_benchmark()`).
#' @param method Scoring method passed to `score_signature()`. Default `"zscore"`.
#' @param time_col,event_col Survival columns in `sample_meta`.
#' @param cutoff,orient_direction Passed to `evaluate_survival()`.
#' @param apply_log2,min_genes Expression prep / scoring thresholds.
#'
#' @return `data.frame` with one row per cohort x signature, including separate
#'   Benjamini-Hochberg FDR columns for the continuous Cox proportional hazards
#'   family (`cox_p` -> `p_adj_bh`) and the median-dichotomized Kaplan-Meier
#'   log-rank family (`logrank_p` -> `logrank_p_adj_bh`), while keeping the
#'   Schoenfeld proportional hazards assumption check (`ph_p`, `ph_violated`)
#'   intentionally unadjusted as a model diagnostic.
#' @export
run_survival_benchmark <- function(datasets,
                                   signatures,
                                   method = "zscore",
                                   time_col = "os_time",
                                   event_col = "os_event",
                                   cutoff = "median",
                                   orient_direction = c("higher_better", "lower_better"),
                                   apply_log2 = TRUE,
                                   min_genes = 3L) {
  if (is.null(names(datasets)) || any(!nzchar(names(datasets)))) {
    stop("`datasets` must be a named list.", call. = FALSE)
  }
  if (is.null(names(signatures)) || any(!nzchar(names(signatures)))) {
    stop("`signatures` must be a named list.", call. = FALSE)
  }
  orient_direction <- match.arg(orient_direction)
  time_col <- as.character(time_col)[1]
  event_col <- as.character(event_col)[1]

  rows <- list()
  row_i <- 0L

  for (ds_name in names(datasets)) {
    prep <- tryCatch(
      .prepare_cohort(
        datasets[[ds_name]],
        name = ds_name,
        apply_log2 = apply_log2,
        filter_response = FALSE
      ),
      error = function(e) e
    )
    if (inherits(prep, "error")) {
      if (length(datasets) == 1L) {
        stop(conditionMessage(prep), call. = FALSE)
      }
      warning(
        "Skipping cohort '", ds_name, "' in run_survival_benchmark(): ",
        conditionMessage(prep),
        call. = FALSE
      )
      next
    }
    ds <- prep$ds
    has_time <- time_col %in% names(ds$sample_meta)
    has_event <- event_col %in% names(ds$sample_meta)
    if (!has_time || !has_event) {
      message(
        "Skipping cohort '", ds_name,
        "': missing survival columns '", time_col, "' / '", event_col, "'."
      )
      next
    }
    expr <- prep$expr

    for (sig_name in names(signatures)) {
      sig <- .as_signature(signatures[[sig_name]], name = sig_name)
      eff_min_genes <- .resolve_min_genes(sig, default = min_genes)
      scored <- score_signature(
        expr,
        genes = sig$genes,
        genes_up = sig$genes_up,
        genes_down = sig$genes_down,
        weights = sig$weights,
        method = method,
        min_genes = eff_min_genes
      )
      if (isTRUE(scored$skipped) || is.null(scored$scores)) {
        next
      }

      # Attach cohort label for evaluate_survival via info$dataset_name
      ds_eval <- ds
      ds_eval$info$dataset_name <- ds_name

      sig_orient <- if (!is.null(sig$direction) &&
                        as.character(sig$direction)[1] %in% c("negative", ">")) {
        "lower_better"
      } else {
        orient_direction
      }
      row <- evaluate_survival(
        dataset = ds_eval,
        score = scored$scores,
        time_col = time_col,
        event_col = event_col,
        cutoff = cutoff,
        orient_direction = sig_orient
      )
      row$signature <- sig_name
      row$method <- method
      row$p_adj_bh <- NA_real_
      row$logrank_p_adj_bh <- NA_real_
      row <- row[, c(
        "cohort", "signature", "method", "endpoint", "n_samples", "n_events",
        "hr", "hr_ci_low", "hr_ci_high", "cox_p", "p_adj_bh", "c_index",
        "logrank_p", "logrank_p_adj_bh", "median_surv_high", "median_surv_low",
        "ph_p", "ph_violated", "surv_note", "small_sample"
      ), drop = FALSE]
      row_i <- row_i + 1L
      rows[[row_i]] <- row
    }
  }

  if (length(rows) == 0L) {
    return(.survival_empty_schema(with_signature = TRUE))
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  ok_p <- is.finite(out$cox_p)
  if (any(ok_p)) {
    out$p_adj_bh[ok_p] <- stats::p.adjust(out$cox_p[ok_p], method = "BH")
  }
  ok_lr <- is.finite(out$logrank_p)
  if (any(ok_lr)) {
    out$logrank_p_adj_bh[ok_lr] <- stats::p.adjust(out$logrank_p[ok_lr], method = "BH")
  }
  out
}

#' Combine response AUC results with survival clinical-benefit metrics
#'
#' Merges on `cohort` and `signature` (`all.x = TRUE`). If `auc_res` uses
#' `dataset` instead of `cohort`, it is copied to `cohort` for the join.
#' Survival columns `hr`, `cox_p`, `logrank_p`, and `c_index` are attached
#' as clinical-benefit fields alongside response discrimination.
#'
#' @param auc_res `run_benchmark()` (or compatible) result table.
#' @param surv_res `run_survival_benchmark()` / `evaluate_survival()` table
#'   (must include `signature` for a useful join).
#'
#' @return Merged `data.frame`.
#' @export
combine_response_and_survival <- function(auc_res, surv_res) {
  if (!is.data.frame(auc_res)) {
    stop("`auc_res` must be a data.frame.", call. = FALSE)
  }
  if (!is.data.frame(surv_res)) {
    stop("`surv_res` must be a data.frame.", call. = FALSE)
  }

  auc <- auc_res
  if (!"cohort" %in% names(auc) && "dataset" %in% names(auc)) {
    auc$cohort <- as.character(auc$dataset)
  }
  if (!"cohort" %in% names(auc)) {
    stop("`auc_res` must contain `cohort` or `dataset`.", call. = FALSE)
  }
  if (!"signature" %in% names(auc)) {
    stop("`auc_res` must contain `signature`.", call. = FALSE)
  }

  if (nrow(surv_res) == 0L || !"signature" %in% names(surv_res)) {
    surv_keep <- data.frame(
      cohort = character(),
      signature = character(),
      endpoint = character(),
      hr = numeric(),
      cox_p = numeric(),
      p_adj_bh = numeric(),
      logrank_p = numeric(),
      logrank_p_adj_bh = numeric(),
      c_index = numeric(),
      ph_p = numeric(),
      ph_violated = logical(),
      stringsAsFactors = FALSE
    )
    return(merge(auc, surv_keep, by = c("cohort", "signature"), all.x = TRUE))
  }

  need <- c("cohort", "signature", "hr", "cox_p", "logrank_p", "c_index")
  miss <- setdiff(need, names(surv_res))
  if (length(miss) > 0L) {
    stop(
      "`surv_res` missing columns: ", paste(miss, collapse = ", "),
      call. = FALSE
    )
  }

  use_method <- "method" %in% names(auc) && "method" %in% names(surv_res)
  if (!use_method && "method" %in% names(auc) && length(unique(auc$method)) > 1L) {
    warning(
      "`surv_res` lacks a `method` column while `auc_res` has multiple methods; ",
      "survival metrics are joined only by (cohort, signature).",
      call. = FALSE
    )
  }

  opt_surv <- intersect(
    c("endpoint", "p_adj_bh", "logrank_p_adj_bh", "ph_p", "ph_violated"),
    names(surv_res)
  )
  surv_cols <- unique(c(need, if (use_method) "method" else character(0), opt_surv))
  surv_keep <- surv_res[, surv_cols, drop = FALSE]

  auc$.join_cohort <- tolower(trimws(as.character(auc$cohort)))
  auc$.join_sig <- tolower(trimws(as.character(auc$signature)))
  surv_keep$.join_cohort <- tolower(trimws(as.character(surv_keep$cohort)))
  surv_keep$.join_sig <- tolower(trimws(as.character(surv_keep$signature)))

  by_keys <- c(".join_cohort", ".join_sig")
  if (use_method) {
    auc$.join_method <- tolower(trimws(as.character(auc$method)))
    surv_keep$.join_method <- tolower(trimws(as.character(surv_keep$method)))
    by_keys <- c(by_keys, ".join_method")
  }

  dedup_keys <- by_keys
  if ("endpoint" %in% names(surv_keep)) {
    surv_keep$.join_endpoint <- tolower(trimws(as.character(surv_keep$endpoint)))
    dedup_keys <- c(dedup_keys, ".join_endpoint")
  }

  dup_mask <- duplicated(surv_keep[, dedup_keys, drop = FALSE])
  if (any(dup_mask)) {
    warning(
      "Duplicate survival rows detected (", sum(dup_mask),
      "); keeping first row per key to prevent Cartesian expansion.",
      call. = FALSE
    )
    surv_keep <- surv_keep[!dup_mask, , drop = FALSE]
  }
  surv_keep$.join_endpoint <- NULL

  # Report unmatched survival keys if any
  auc_keys <- do.call(paste, c(auc[, by_keys, drop = FALSE], list(sep = "\r")))
  surv_keys <- do.call(paste, c(surv_keep[, by_keys, drop = FALSE], list(sep = "\r")))
  unmatched <- setdiff(surv_keys, auc_keys)
  if (length(unmatched) > 0L) {
    warning(
      length(unmatched), " survival row(s) had no matching row in `auc_res`.",
      call. = FALSE
    )
  }

  payload_cols <- unique(c("hr", "cox_p", "logrank_p", "c_index", opt_surv))
  surv_payload <- surv_keep[, c(by_keys, payload_cols), drop = FALSE]
  out <- merge(auc, surv_payload, by = by_keys, all.x = TRUE, sort = FALSE)
  for (jk in by_keys) {
    out[[jk]] <- NULL
  }
  out
}
