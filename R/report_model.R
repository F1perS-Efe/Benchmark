#' Load deterministic interpretation rules from YAML
#'
#' Reads `inst/extdata/interpretation_rules.yml` (or a user-supplied YAML path)
#' and attaches its MD5 checksum for report provenance.
#'
#' @param rules_path Optional path to a custom YAML rules file. `NULL` uses the
#'   bundled `interpretation_rules.yml`.
#'
#' @return A list of parsed rules with `version`, `md5`, and `path` metadata.
#' @export
load_interpretation_rules <- function(rules_path = NULL) {
  if (is.null(rules_path) || !nzchar(as.character(rules_path)[1])) {
    local_cand <- file.path("inst", "extdata", "interpretation_rules.yml")
    up_cand <- file.path("..", "..", "inst", "extdata", "interpretation_rules.yml")
    if (file.exists(local_cand)) {
      rules_path <- local_cand
    } else if (file.exists(up_cand)) {
      rules_path <- up_cand
    } else {
      rules_path <- system.file("extdata", "interpretation_rules.yml", package = "SigBenchR")
    }
  }
  if (!nzchar(rules_path) || !file.exists(rules_path)) {
    stop("Interpretation rules YAML file not found in installed package (extdata/interpretation_rules.yml).", call. = FALSE)
  }
  rules <- yaml::read_yaml(rules_path)
  rules$path <- normalizePath(rules_path, mustWork = FALSE)
  rules$md5 <- unname(tools::md5sum(rules_path))
  rules
}

#' Build a 3-tier interpretive report model from benchmark outputs
#'
#' Evaluates deterministic YAML rules (`load_interpretation_rules()`) against
#' benchmark results, meta-analysis summaries, optional reference-axis
#' correlations, and optional survival results, returning a structured 3-tier
#' object (`tier1`, `tier2`, `tier3`) decoupled from R Markdown rendering.
#'
#' @param results `run_benchmark()` output `data.frame`.
#' @param summary_tbl Optional `summarize_auc()` output `data.frame`. When `NULL`,
#'   computed automatically for the active method.
#' @param sensitivity_tbl Optional `check_sensitivity()` output `data.frame`.
#' @param params_list Optional list of run parameters (e.g. `heatmap_method`).
#' @param ref_axes Optional `correlate_reference_axes()` output `data.frame`.
#' @param surv_res Optional `run_survival_benchmark()` output `data.frame`.
#' @param beyond_inf Optional `test_beyond_infiltration()` output (either a
#'   list with `per_cohort` and `summary` or a named list of such results per
#'   signature).
#' @param rules_path Optional path to a custom `interpretation_rules.yml`.
#'
#' @return Named list with `tier1` (executive summary, best signature, decision
#'   cards, deterministic alerts), `tier2` (comparative evidence, clinical
#'   strata, reference-axis overlap, beyond infiltration, survival), and `tier3`
#'   (cohort metrics, coverage diagnostics, sensitivity, run parameters, rules
#'   metadata).
#' @export
build_report_model <- function(results,
                               summary_tbl = NULL,
                               sensitivity_tbl = NULL,
                               params_list = list(),
                               ref_axes = NULL,
                               surv_res = NULL,
                               beyond_inf = NULL,
                               rules_path = NULL) {
  if (!is.data.frame(results)) {
    stop("`results` must be a data.frame.", call. = FALSE)
  }
  if (is.null(params_list)) {
    params_list <- list()
  }
  rules <- load_interpretation_rules(rules_path)

  active_method <- if (!is.null(params_list$heatmap_method) &&
                       nzchar(as.character(params_list$heatmap_method)[1])) {
    as.character(params_list$heatmap_method)[1]
  } else if (!is.null(rules$primary_method) &&
             "method" %in% names(results) &&
             rules$primary_method %in% results$method) {
    as.character(rules$primary_method)[1]
  } else if ("method" %in% names(results) && nrow(results) > 0L) {
    as.character(results$method[1])
  } else {
    "mean"
  }

  if (is.null(summary_tbl) && nrow(results) > 0L) {
    summary_tbl <- suppressWarnings(
      summarize_auc(results, method = active_method, stratify_by = NULL)
    )
  }
  if (is.null(summary_tbl)) {
    summary_tbl <- data.frame()
  }

  cards <- .build_signature_cards(summary_tbl, rules)
  best_sig <- .pick_best_signature(cards)
  bi_norm <- .normalize_beyond_inf_for_report(beyond_inf)
  alerts <- .evaluate_report_alerts(
    results = results,
    summary_tbl = summary_tbl,
    ref_axes = ref_axes,
    rules = rules,
    surv_res = surv_res,
    beyond_inf = bi_norm
  )

  # Tier 2: Comparative evidence & clinical strata
  summary_formatted <- if (nrow(summary_tbl) > 0L) {
    .format_meta_summary(summary_tbl, stratum_col = NULL)
  } else {
    summary_tbl
  }

  tx_strata <- NULL
  if (nrow(results) > 0L && "treatment" %in% names(results)) {
    tx_strata <- tryCatch(
      suppressWarnings(summarize_auc(results, method = active_method, stratify_by = "treatment")),
      error = function(e) NULL
    )
  }
  ep_strata <- NULL
  if (nrow(results) > 0L && "endpoint_type" %in% names(results)) {
    ep_strata <- tryCatch(
      suppressWarnings(summarize_auc(results, method = active_method, stratify_by = "endpoint_type")),
      error = function(e) NULL
    )
  }

  # Tier 3: Technical appendix
  cov_cols <- intersect(
    c(
      "dataset", "signature", "method", "n_found", "n_total",
      "coverage_pct", "coverage_status", "low_coverage",
      "n_alias_rescued", "missing_genes"
    ),
    names(results)
  )
  cov_diag <- if (length(cov_cols) > 0L && nrow(results) > 0L) {
    results[, cov_cols, drop = FALSE]
  } else {
    data.frame()
  }

  list(
    tier1 = list(
      active_method = active_method,
      best_signature = best_sig,
      signature_cards = cards,
      alerts = alerts
    ),
    tier2 = list(
      summary_tbl = summary_tbl,
      summary_formatted = summary_formatted,
      treatment_strata = tx_strata,
      endpoint_strata = ep_strata,
      ref_axes = ref_axes,
      beyond_infiltration = bi_norm,
      survival = surv_res
    ),
    tier3 = list(
      cohort_metrics = results,
      coverage_diagnostics = cov_diag,
      sensitivity = sensitivity_tbl,
      run_params = params_list,
      params = params_list,
      rules_meta = list(
        version = rules$version,
        md5 = rules$md5,
        path = rules$path
      )
    )
  )
}

#' @noRd
.render_rule_template <- function(tpl, vals) {
  out <- as.character(tpl)[1]
  for (nm in names(vals)) {
    out <- gsub(paste0("\\{", nm, "\\}"), as.character(vals[[nm]])[1], out)
  }
  out
}

#' @noRd
.build_signature_cards <- function(summary_tbl, rules) {
  if (!is.data.frame(summary_tbl) || nrow(summary_tbl) == 0L) {
    return(data.frame(
      signature = character(),
      method = character(),
      k = integer(),
      auc_re = numeric(),
      auc_re_ci_low = numeric(),
      auc_re_ci_high = numeric(),
      i2_logit = numeric(),
      priority_tier = character(),
      rule_id = character(),
      interpretation = character(),
      stringsAsFactors = FALSE
    ))
  }

  tiers <- rules$evidence_tiers
  k_vec <- if ("n_cohorts" %in% names(summary_tbl)) summary_tbl$n_cohorts else summary_tbl$n_datasets
  i2_vec <- if ("i2_logit" %in% names(summary_tbl)) summary_tbl$i2_logit else rep(NA_real_, nrow(summary_tbl))

  rows <- lapply(seq_len(nrow(summary_tbl)), function(i) {
    sig <- as.character(summary_tbl$signature[i])
    meth <- as.character(summary_tbl$method[i])
    k <- as.integer(k_vec[i])
    auc_est <- if (is.finite(summary_tbl$auc_re[i])) summary_tbl$auc_re[i] else summary_tbl$auc_n_weighted[i]
    ci_lo <- summary_tbl$auc_re_ci_low[i]
    ci_hi <- summary_tbl$auc_re_ci_high[i]
    i2_val <- i2_vec[i]

    fmt_num <- function(x) if (is.finite(x)) sprintf("%.3f", x) else "NA"

    if (is.na(k) || k == 0L || !is.finite(auc_est)) {
      rule <- tiers$unevaluable
    } else if (k >= tiers$strong$min_cohorts &&
               auc_est >= tiers$strong$min_auc &&
               is.finite(ci_lo) && ci_lo > tiers$strong$min_ci_low &&
               (!is.finite(i2_val) || i2_val < tiers$strong$max_i2_logit)) {
      rule <- tiers$strong
    } else if (k >= tiers$moderate$min_cohorts &&
               auc_est >= tiers$moderate$min_auc &&
               (!is.finite(ci_lo) || ci_lo > tiers$moderate$min_ci_low)) {
      rule <- tiers$moderate
    } else if (k >= tiers$moderate$min_cohorts &&
               auc_est >= tiers$moderate$min_auc) {
      # High point AUC (>= 0.60) whose 95% CI crosses below min_ci_low / 0.50:
      # flag as exploratory rather than none or Medium Priority
      rule <- tiers$moderate
    } else {
      rule <- tiers$none
    }

    interp <- .render_rule_template(rule$template, list(
      signature = sig,
      k = k,
      auc_re = fmt_num(auc_est),
      ci_low = fmt_num(ci_lo),
      ci_high = fmt_num(ci_hi)
    ))

    priority_lbl <- as.character(rule$label)
    if (identical(rule$id, tiers$moderate$id) &&
        is.finite(ci_lo) && ci_lo <= rules$auc$chance_level) {
      priority_lbl <- "Exploratory Signal (CI crosses 0.50)"
      interp <- paste0(
        interp,
        " Note: 95% CI lower bound (", fmt_num(ci_lo), ") is <= 0.50."
      )
    }

    data.frame(
      signature = sig,
      method = meth,
      k = k,
      auc_re = as.numeric(auc_est),
      auc_re_ci_low = as.numeric(ci_lo),
      auc_re_ci_high = as.numeric(ci_hi),
      i2_logit = as.numeric(i2_val),
      priority_tier = priority_lbl,
      rule_id = as.character(rule$id),
      interpretation = interp,
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)
  tier_rank <- ifelse(
    out$rule_id == "E1_STRONG", 3L,
    ifelse(out$rule_id == "E2_MODERATE", 2L,
      ifelse(out$rule_id == "E3_NONE", 1L, 0L)
    )
  )
  ci_above_chance <- is.finite(out$auc_re_ci_low) & (out$auc_re_ci_low > 0.50)
  ord <- order(
    -tier_rank,
    !ci_above_chance,
    -ifelse(is.finite(out$k), out$k, -1L),
    !is.finite(out$auc_re),
    -out$auc_re,
    out$signature
  )
  out <- out[ord, , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' @noRd
.pick_best_signature <- function(cards) {
  if (!is.data.frame(cards) || nrow(cards) == 0L) {
    return(NULL)
  }
  valid <- cards[cards$k > 0L & is.finite(cards$auc_re), , drop = FALSE]
  if (nrow(valid) == 0L) {
    return(NULL)
  }
  builtin_all <- tryCatch(
    get_builtin_signatures(include_controls = TRUE),
    error = function(e) list()
  )
  neg_ctrl_names <- names(builtin_all)[vapply(builtin_all, function(s) {
    is.list(s) && identical(s$category, "negative_control")
  }, logical(1))]
  is_neg_ctrl <- (valid$signature %in% neg_ctrl_names) |
    grepl("^NegCtrl_|^Random", valid$signature)
  if ("category" %in% names(valid)) {
    is_neg_ctrl <- is_neg_ctrl | (valid$category == "negative_control")
  }
  valid <- valid[!is_neg_ctrl, , drop = FALSE]
  if (nrow(valid) == 0L) {
    return(NULL)
  }
  tier_rank <- ifelse(
    valid$rule_id == "E1_STRONG", 3L,
    ifelse(valid$rule_id == "E2_MODERATE", 2L,
      ifelse(valid$rule_id == "E3_NONE", 1L, 0L)
    )
  )
  ci_above_chance <- is.finite(valid$auc_re_ci_low) & (valid$auc_re_ci_low > 0.50)
  ord <- order(
    -tier_rank,
    !ci_above_chance,
    -ifelse(is.finite(valid$k), valid$k, -1L),
    -valid$auc_re,
    valid$signature
  )
  valid <- valid[ord, , drop = FALSE]
  rownames(valid) <- NULL
  valid[1L, , drop = FALSE]
}

#' @noRd
.normalize_beyond_inf_for_report <- function(beyond_inf) {
  if (is.null(beyond_inf)) {
    return(NULL)
  }
  if (is.list(beyond_inf) && all(c("per_cohort", "summary") %in% names(beyond_inf))) {
    pc <- beyond_inf$per_cohort
    sm <- beyond_inf$summary
    if (is.data.frame(pc) && nrow(pc) > 0L && !"signature" %in% names(pc)) {
      pc$signature <- if (!is.null(beyond_inf$signature)) as.character(beyond_inf$signature)[1] else "Selected"
    }
    if (is.data.frame(sm) && nrow(sm) > 0L && !"signature" %in% names(sm)) {
      sm$signature <- if (!is.null(beyond_inf$signature)) as.character(beyond_inf$signature)[1] else "Selected"
    }
    return(list(per_cohort = pc, summary = sm))
  }
  if (is.list(beyond_inf) && length(beyond_inf) > 0L) {
    pc_list <- list()
    sm_list <- list()
    for (nm in names(beyond_inf)) {
      item <- beyond_inf[[nm]]
      if (is.list(item) && all(c("per_cohort", "summary") %in% names(item))) {
        if (is.data.frame(item$per_cohort) && nrow(item$per_cohort) > 0L) {
          df_p <- item$per_cohort
          df_p$signature <- nm
          pc_list[[length(pc_list) + 1L]] <- df_p
        }
        if (is.data.frame(item$summary) && nrow(item$summary) > 0L) {
          df_s <- item$summary
          df_s$signature <- nm
          sm_list[[length(sm_list) + 1L]] <- df_s
        }
      }
    }
    pc_all <- if (length(pc_list) > 0L) do.call(rbind, pc_list) else data.frame()
    sm_all <- if (length(sm_list) > 0L) do.call(rbind, sm_list) else data.frame()
    rownames(pc_all) <- NULL
    rownames(sm_all) <- NULL
    return(list(per_cohort = pc_all, summary = sm_all))
  }
  NULL
}

#' @noRd
.evaluate_report_alerts <- function(results,
                                    summary_tbl,
                                    ref_axes,
                                    rules,
                                    surv_res = NULL,
                                    beyond_inf = NULL) {
  alerts <- list()
  add_alert <- function(rule, vals, metric_val) {
    alerts[[length(alerts) + 1L]] <<- data.frame(
      rule_id = as.character(rule$id),
      severity = as.character(rule$severity),
      title = as.character(rule$title),
      message = .render_rule_template(rule$template, vals),
      metric_value = as.character(metric_val),
      stringsAsFactors = FALSE
    )
  }
  add_direct_alert <- function(rule_id, severity, title, message, metric_val) {
    alerts[[length(alerts) + 1L]] <<- data.frame(
      rule_id = as.character(rule_id),
      severity = as.character(severity),
      title = as.character(title),
      message = as.character(message),
      metric_value = as.character(metric_val),
      stringsAsFactors = FALSE
    )
  }

  # 1) CI crosses chance level (<= 0.50)
  if (is.data.frame(summary_tbl) && nrow(summary_tbl) > 0L && "auc_re_ci_low" %in% names(summary_tbl)) {
    k_vec <- if ("n_cohorts" %in% names(summary_tbl)) summary_tbl$n_cohorts else summary_tbl$n_datasets
    bad_ci <- which(k_vec > 0L & is.finite(summary_tbl$auc_re_ci_low) &
                      summary_tbl$auc_re_ci_low <= rules$auc$chance_level)
    if (length(bad_ci) > 0L) {
      min_lo <- min(summary_tbl$auc_re_ci_low[bad_ci], na.rm = TRUE)
      add_alert(rules$auc$ci_chance_rule, list(ci_low = sprintf("%.3f", min_lo)), sprintf("%.3f", min_lo))
    }
  }

  # 2) Logit-scale heterogeneity (i2_logit)
  if (is.data.frame(summary_tbl) && nrow(summary_tbl) > 0L && "i2_logit" %in% names(summary_tbl)) {
    i2_vals <- summary_tbl$i2_logit[is.finite(summary_tbl$i2_logit)]
    if (length(i2_vals) > 0L) {
      max_i2 <- max(i2_vals)
      if (max_i2 >= rules$heterogeneity$high_i2_logit) {
        add_alert(rules$heterogeneity$high_rule, list(i2_logit = sprintf("%.1f", max_i2)), sprintf("%.1f%%", max_i2))
      } else if (max_i2 >= rules$heterogeneity$moderate_i2_logit) {
        add_alert(rules$heterogeneity$moderate_rule, list(i2_logit = sprintf("%.1f", max_i2)), sprintf("%.1f%%", max_i2))
      }
    }
  }

  # 3) Low gene coverage (< ok_min_pct)
  if (is.data.frame(results) && nrow(results) > 0L) {
    low_mask <- if ("low_coverage" %in% names(results)) {
      results$low_coverage %in% TRUE
    } else if ("coverage_pct" %in% names(results)) {
      is.finite(results$coverage_pct) & results$coverage_pct < rules$coverage$ok_min_pct
    } else {
      rep(FALSE, nrow(results))
    }
    n_low <- sum(low_mask, na.rm = TRUE)
    if (n_low > 0L) {
      add_alert(rules$coverage$low_rule, list(n_low_cov = n_low), as.character(n_low))
    }
  }

  # 4) Circular reference-axis overlap
  if (is.data.frame(ref_axes) && nrow(ref_axes) > 0L) {
    ov_mask <- if ("overlap_flag" %in% names(ref_axes)) {
      ref_axes$overlap_flag %in% TRUE
    } else if ("jaccard_overlap" %in% names(ref_axes)) {
      is.finite(ref_axes$jaccard_overlap) & ref_axes$jaccard_overlap >= rules$reference_overlap$jaccard_threshold
    } else {
      rep(FALSE, nrow(ref_axes))
    }
    n_ov <- sum(ov_mask, na.rm = TRUE)
    if (n_ov > 0L) {
      add_alert(rules$reference_overlap$overlap_rule, list(n_overlap = n_ov), as.character(n_ov))
    }
  }

  # 5) Subgroup pooling warning
  if (is.data.frame(summary_tbl) && nrow(summary_tbl) > 0L && "note" %in% names(summary_tbl)) {
    pooled_hit <- any(grepl("Pooled across distinct|stratify_by", as.character(summary_tbl$note), ignore.case = TRUE), na.rm = TRUE)
    if (isTRUE(pooled_hit)) {
      add_alert(rules$subgroup_pooling$pooled_rule, list(), "pooled")
    }
  }

  # 6) Small-k exploratory meta-analysis guard (k < min_cohorts_confirmatory)
  if (!is.null(rules$small_k_meta) && is.data.frame(summary_tbl) && nrow(summary_tbl) > 0L) {
    k_vec <- if ("n_cohorts" %in% names(summary_tbl)) summary_tbl$n_cohorts else summary_tbl$n_datasets
    pos_k <- k_vec[is.finite(k_vec) & k_vec > 0L]
    if (length(pos_k) > 0L && max(pos_k) < rules$small_k_meta$min_cohorts_confirmatory) {
      add_alert(
        rules$small_k_meta$small_k_rule,
        list(k = as.integer(max(pos_k))),
        as.character(max(pos_k))
      )
    }
  }

  # 7) PR-AUC small positive/negative class & unreliable CI diagnostics
  if (is.data.frame(results) && nrow(results) > 0L) {
    if ("small_positive_class" %in% names(results) && any(results$small_positive_class %in% TRUE, na.rm = TRUE)) {
      sc_coh <- unique(as.character(results$dataset[results$small_positive_class %in% TRUE]))
      add_direct_alert(
        "PR-SMALL-CLASS",
        "warning",
        "Small positive/negative class in PR-AUC",
        sprintf(
          "Cohort(s) %s have fewer than 5 responder or non-responder samples (n_pos < 5 or n_neg < 5); PR-AUC point estimates and stratified bootstrap intervals require cautious interpretation.",
          paste(sc_coh, collapse = ", ")
        ),
        paste(sc_coh, collapse = ",")
      )
    }
    if ("ci_unreliable" %in% names(results) && any(results$ci_unreliable %in% TRUE, na.rm = TRUE)) {
      n_unr <- sum(results$ci_unreliable %in% TRUE, na.rm = TRUE)
      add_direct_alert(
        "PR-CI-UNSTABLE",
        "warning",
        "Unstable PR-AUC bootstrap confidence interval",
        sprintf(
          "%d evaluation(s) flagged ci_unreliable = TRUE due to small minority class size, heavy score ties (>50%%), or <50%% valid bootstrap resamples (n_boot_valid).",
          n_unr
        ),
        as.character(n_unr)
      )
    }
    if ("background_type" %in% names(results)) {
      bgs <- unique(stats::na.omit(as.character(results$background_type)))
      if (length(bgs) > 0L) {
        sev <- if ("partial_ortholog_subset" %in% bgs) "warning" else "info"
        add_direct_alert(
          "NULL-UNIVERSE",
          sev,
          "Null permutation background universe",
          sprintf(
            "Null permutation random gene sets were drawn from background universe type(s): %s. Raw null_p, within-method BH FDR (null_p_adj_bh), and across-method global BH FDR (null_p_adj_bh_global) are reported separately.",
            paste(bgs, collapse = ", ")
          ),
          paste(bgs, collapse = ",")
        )
      }
    }
  }

  # 8) Beyond Infiltration diagnostics (direction_conflict, low_epv, high_collinearity)
  if (is.list(beyond_inf)) {
    sm_bi <- beyond_inf$summary
    pc_bi <- beyond_inf$per_cohort
    if (is.data.frame(sm_bi) && nrow(sm_bi) > 0L && "direction_conflict" %in% names(sm_bi) &&
        any(sm_bi$direction_conflict %in% TRUE, na.rm = TRUE)) {
      add_direct_alert(
        "BI-DIRECTION-CONFLICT",
        "warning",
        "Cross-cohort effect direction conflict in Beyond Infiltration",
        "Adjusted signature effects (beta_adjusted) have opposite signs across cohorts (direction_conflict = TRUE). Directional Stouffer Z combined_lrt_p accounts for effect cancellation; do not use unsigned Fisher omnibus p to claim consistent independent value.",
        "TRUE"
      )
    }
    if (is.data.frame(pc_bi) && nrow(pc_bi) > 0L) {
      if ("low_epv" %in% names(pc_bi) && any(pc_bi$low_epv %in% TRUE, na.rm = TRUE)) {
        min_epv <- min(pc_bi$epv[pc_bi$low_epv %in% TRUE], na.rm = TRUE)
        add_direct_alert(
          "BI-LOW-EPV",
          "warning",
          "Low Events Per Variable (EPV < 5) in adjusted logistic model",
          sprintf(
            "One or more cohorts have EPV < 5 (minimum EPV = %.1f) in the 2-predictor logistic regression (response ~ sig_score_z + inf_score_z). Treat adjusted coefficients and LRT p-values as exploratory.",
            min_epv
          ),
          sprintf("%.1f", min_epv)
        )
      }
      if ("high_collinearity" %in% names(pc_bi) && any(pc_bi$high_collinearity %in% TRUE, na.rm = TRUE)) {
        max_vif <- max(pc_bi$vif[pc_bi$high_collinearity %in% TRUE], na.rm = TRUE)
        add_direct_alert(
          "BI-HIGH-VIF",
          "warning",
          "Multicollinearity with reference infiltration axis (VIF >= 5)",
          sprintf(
            "One or more evaluations exhibit high Variance Inflation Factor (max VIF = %.2f) with the reference infiltration axis. Evaluations with VIF >= 10 are marked failed_collinearity.",
            max_vif
          ),
          sprintf("%.2f", max_vif)
        )
      }
    }
  }

  # 9) Survival proportional hazards assumption diagnostic check
  if (is.data.frame(surv_res) && nrow(surv_res) > 0L && "ph_violated" %in% names(surv_res) &&
      any(surv_res$ph_violated %in% TRUE, na.rm = TRUE)) {
    n_ph <- sum(surv_res$ph_violated %in% TRUE, na.rm = TRUE)
    add_direct_alert(
      "SURV-PH-VIOLATION",
      "warning",
      "Proportional hazards assumption diagnostic violation (ph_p < 0.05)",
      sprintf(
        "%d survival evaluation(s) have unadjusted Schoenfeld residual diagnostic ph_p < 0.05 (ph_violated = TRUE). Interpret constant Cox hazard ratios with caution.",
        n_ph
      ),
      as.character(n_ph)
    )
  }

  if (length(alerts) == 0L) {
    return(data.frame(
      rule_id = character(),
      severity = character(),
      title = character(),
      message = character(),
      metric_value = character(),
      stringsAsFactors = FALSE
    ))
  }
  do.call(rbind, alerts)
}
