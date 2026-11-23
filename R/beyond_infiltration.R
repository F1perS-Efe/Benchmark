#' Test whether a signature adds predictive value beyond T-cell infiltration
#'
#' Evaluates out-of-fold discrimination via repeated stratified K-fold
#' cross-validation (`n_folds` x `n_repeats`) comparing two logistic regression
#' models per cohort:
#' - Base model: `response ~ inf_score_z`
#' - Full model: `response ~ sig_score_z + inf_score_z`
#'
#' To prevent preprocessing leakage in cross-validation (`4.1`):
#' - For cohort-dependent gene transformations (`method = "zscore"` or
#'   `"rank"`), per-gene scaling parameters (training mean/SD or training
#'   empirical rank distribution) are estimated **strictly within each training
#'   fold** and applied to the held-out test fold before aggregating gene scores.
#' - Sample-independent methods (`"mean"`, `"median"`, `"ss_rank"`) are computed
#'   per sample without cross-sample information.
#' - Predictor standardization (`sig_score_z`, `inf_score_z`) inside CV is
#'   always learned on the training fold and applied to the test fold.
#'
#' CV discrimination (`4.2`) is summarized via paired fold-level test AUCs across
#' repeats (with `sd_delta_auc`, `delta_auc_ci_low`, `delta_auc_ci_high`, and
#' pooled OOF diagnostic `delta_auc_oof`).
#'
#' Multiple testing and directional synthesis (`4.3`, `4.4`, `4.5`):
#' - `adds_independent_value` requires `lrt_p_adj_bh < 0.05` (Benjamini-Hochberg
#'   adjusted across the evaluated cohort family), positive adjusted effect
#'   (`or_adjusted > 1`), positive CV `delta_auc_combined > 0`, model
#'   convergence, no separation, and absence of severe multicollinearity
#'   (`vif < 10`).
#' - Cross-cohort synthesis uses sample-size weighted **directional Stouffer Z**
#'   (`combined_lrt_p`, `stouffer_z`, `direction_conflict`) so opposite-direction
#'   cohort effects cancel rather than reinforcing into a false-positive Fisher
#'   p-value.
#' - Reports events-per-variable (`epv = min(n_pos, n_neg) / 2`), variance
#'   inflation factor (`vif = 1 / (1 - r^2)`), `high_collinearity` (`vif >= 5`),
#'   `low_epv` (`epv < 5`), and `analysis_status`.
#'
#' @param datasets Named list of SigBenchR datasets (`validate_dataset()`).
#' @param signature A signature definition (character vector, named sign/weight
#'   vector, or list with `genes` / `genes_up` / `genes_down` / `weights`).
#' @param infiltration_axis Reference axis name: `"CD8_T_Cell"` (default, uses
#'   `CD3D/CD3E/CD8A/CD8B`), `"T_Cell_Infiltration"`, `"IFN_Gamma"`, or
#'   `"Cytotoxicity"`. Ignored when `infiltration_genes` is provided.
#' @param infiltration_genes Optional custom character vector of infiltration
#'   reference genes.
#' @param method Scoring method passed to [score_signature()]. Default `"zscore"`.
#' @param sd_handling Response policy (`"exclude"`, `"as_0"`, `"as_1"`).
#' @param direction Default ROC direction (`"<"` or `">"`). Overridden by
#'   `signature$direction` when present.
#' @param apply_log2 Passed to [prepare_expression()]. Default `TRUE`.
#' @param min_genes Minimum matched signature genes per cohort. Default `2L`.
#' @param n_folds Number of stratified CV folds for out-of-fold AUC. Default `5L`.
#' @param n_repeats Number of repeated stratified CV splits. Default `5L`.
#' @param seed RNG seed for CV fold splits (`withr::with_seed()`). Default `42L`.
#'
#' @return A list with `per_cohort` (`data.frame`) and `summary` (`data.frame`
#'   with `n_cohorts`, `n_independent`, `combined_lrt_p`, `fisher_omnibus_p`,
#'   `stouffer_z`, `direction_conflict`, `min_lrt_p_adj_bh`, `median_lrt_p`,
#'   `mean_delta_auc`).
#' @export
test_beyond_infiltration <- function(datasets,
                                     signature,
                                     infiltration_axis = "CD8_T_Cell",
                                     infiltration_genes = NULL,
                                     method = "zscore",
                                     sd_handling = "exclude",
                                     direction = "<",
                                     apply_log2 = TRUE,
                                     min_genes = 2L,
                                     n_folds = 5L,
                                     n_repeats = 5L,
                                     seed = 42L) {
  if (is.null(names(datasets)) || any(!nzchar(names(datasets)))) {
    stop("`datasets` must be a named list.", call. = FALSE)
  }
  method <- match.arg(method, c("mean", "zscore", "median", "rank", "ss_rank"))
  sig <- .as_signature(signature)
  sig_dir <- .normalize_roc_direction(sig$direction, default = direction)

  inf_genes <- .resolve_infiltration_genes(infiltration_axis, infiltration_genes)
  u_genes <- .extract_sig_genes_upper(sig)
  i_genes <- unique(toupper(inf_genes))
  jacc <- if (length(union(u_genes, i_genes)) == 0L) {
    0
  } else {
    length(intersect(u_genes, i_genes)) / length(union(u_genes, i_genes))
  }

  rows <- list()
  row_i <- 0L

  for (ds_name in names(datasets)) {
    prep <- tryCatch(
      .prepare_cohort(
        datasets[[ds_name]],
        name = ds_name,
        apply_log2 = apply_log2,
        sd_handling = sd_handling,
        filter_response = TRUE
      ),
      error = function(e) e
    )
    if (inherits(prep, "error")) {
      if (length(datasets) == 1L) {
        stop(conditionMessage(prep), call. = FALSE)
      }
      warning(
        "Skipping cohort '", ds_name, "' in test_beyond_infiltration(): ",
        conditionMessage(prep),
        call. = FALSE
      )
      next
    }
    expr_full <- prep$expr_kept
    y_kept <- prep$y_kept

    eff_min_genes <- .resolve_min_genes(sig, default = min_genes)
    s_sig <- score_signature(
      expr_full,
      genes = sig$genes,
      genes_up = sig$genes_up,
      genes_down = sig$genes_down,
      weights = sig$weights,
      method = method,
      min_genes = eff_min_genes
    )
    s_inf <- score_signature(
      expr_full,
      genes = inf_genes,
      method = method,
      min_genes = min(2L, length(inf_genes))
    )
    if (isTRUE(s_sig$skipped) || isTRUE(s_inf$skipped) ||
        is.null(s_sig$scores) || is.null(s_inf$scores)) {
      next
    }

    common <- intersect(names(s_sig$scores), names(s_inf$scores))
    common <- intersect(common, names(y_kept))
    if (length(common) < 8L) {
      next
    }

    x_sig <- orient_score(s_sig$scores[common], direction = sig_dir)
    x_inf <- as.numeric(s_inf$scores[common])
    resp <- as.integer(y_kept[common])
    ok <- is.finite(x_sig) & is.finite(x_inf) & !is.na(resp)
    if (sum(ok) < 8L || length(unique(resp[ok])) < 2L) {
      next
    }

    samp_ok <- common[ok]
    x_sig <- x_sig[ok]
    x_inf <- x_inf[ok]
    resp <- resp[ok]
    if (!is.finite(stats::sd(x_sig)) || stats::sd(x_sig) == 0 ||
        !is.finite(stats::sd(x_inf)) || stats::sd(x_inf) == 0) {
      next
    }

    n_pos <- sum(resp == 1L)
    n_neg <- sum(resp == 0L)
    epv <- min(n_pos, n_neg) / 2
    low_epv <- epv < 5.0

    df_glm <- data.frame(
      response = resp,
      sig_score_z = as.numeric(scale(x_sig)),
      inf_score_z = as.numeric(scale(x_inf))
    )
    r_inf <- stats::cor(df_glm$sig_score_z, df_glm$inf_score_z, method = "pearson")
    vif_val <- if (is.finite(r_inf) && abs(r_inf) < 1 - 1e-12) {
      1 / (1 - r_inf^2)
    } else {
      Inf
    }
    high_collinearity <- !is.finite(vif_val) || vif_val >= 5.0
    severe_collinearity <- !is.finite(vif_val) || vif_val >= 10.0

    # Full-cohort GLM with explicit warning capture & separation checks
    warn_msgs <- character(0)
    fit_res <- tryCatch(
      withCallingHandlers(
        {
          m_base <- stats::glm(response ~ inf_score_z, data = df_glm, family = stats::binomial())
          m_full <- stats::glm(response ~ sig_score_z + inf_score_z, data = df_glm, family = stats::binomial())
          cf <- summary(m_full)$coefficients
          beta_adj <- if ("sig_score_z" %in% rownames(cf)) as.numeric(cf["sig_score_z", "Estimate"]) else NA_real_
          or_adj <- if (is.finite(beta_adj)) exp(beta_adj) else NA_real_
          p_adj <- if ("sig_score_z" %in% rownames(cf)) as.numeric(cf["sig_score_z", "Pr(>|z|)"]) else NA_real_
          lrt <- stats::anova(m_base, m_full, test = "Chisq")
          lrt_p <- as.numeric(lrt[2L, "Pr(>Chi)"])
          lrt_deviance <- as.numeric(lrt[2L, "Deviance"])

          coefs_full <- stats::coef(m_full)
          fitted_full <- stats::fitted(m_full)
          conv_ok <- isTRUE(m_base$converged) && isTRUE(m_full$converged)
          sep_flag <- any(!is.finite(coefs_full) | abs(coefs_full) > 15) ||
            any(!is.finite(fitted_full) | fitted_full < 1e-6 | fitted_full > 1 - 1e-6)
          list(
            beta_adj = beta_adj,
            or_adj = or_adj,
            p_adj = p_adj,
            lrt_p = lrt_p,
            lrt_deviance = lrt_deviance,
            converged = conv_ok,
            sep_flag = sep_flag
          )
        },
        warning = function(w) {
          warn_msgs <<- c(warn_msgs, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      ),
      error = function(e) NULL
    )
    if (is.null(fit_res)) {
      next
    }
    sep_warn <- any(grepl("fitted probabilities numerically 0 or 1|did not converge", warn_msgs, ignore.case = TRUE))
    separation_flag <- isTRUE(fit_res$sep_flag) || isTRUE(sep_warn)
    converged <- isTRUE(fit_res$converged) && !any(grepl("did not converge", warn_msgs, ignore.case = TRUE))

    # Leakage-free repeated stratified K-fold CV AUC (4.1 & 4.2)
    expr_cv <- expr_full[, samp_ok, drop = FALSE]
    if (!is.null(attr(expr_full, "transcriptome_bg_expr"))) {
      bg_full <- attr(expr_full, "transcriptome_bg_expr")
      if (is.matrix(bg_full) && all(samp_ok %in% colnames(bg_full))) {
        attr(expr_cv, "transcriptome_bg_expr") <- bg_full[, samp_ok, drop = FALSE]
      }
    }
    cv_res <- .cv_glm_auc(
      df_glm = data.frame(response = resp, x_sig_raw = x_sig, x_inf_raw = x_inf),
      n_folds = as.integer(n_folds),
      n_repeats = as.integer(n_repeats),
      seed = seed,
      expr_sub = expr_cv,
      sig = sig,
      sig_dir = sig_dir,
      inf_genes = inf_genes,
      method = method,
      eff_min_genes = eff_min_genes
    )

    status_str <- if (!isTRUE(converged) || isTRUE(separation_flag)) {
      "failed_separation_or_convergence"
    } else if (isTRUE(severe_collinearity)) {
      "failed_collinearity"
    } else if (isTRUE(low_epv) || isTRUE(high_collinearity)) {
      "exploratory_low_epv_or_collinear"
    } else {
      "ok"
    }

    row_i <- row_i + 1L
    rows[[row_i]] <- data.frame(
      cohort = ds_name,
      n = as.integer(nrow(df_glm)),
      n_pos = as.integer(n_pos),
      n_neg = as.integer(n_neg),
      epv = as.numeric(epv),
      low_epv = isTRUE(low_epv),
      auc_infiltration_only = as.numeric(cv_res$auc_inf),
      auc_combined = as.numeric(cv_res$auc_comb),
      delta_auc_combined = as.numeric(cv_res$delta_auc),
      sd_delta_auc = as.numeric(cv_res$sd_delta_auc),
      delta_min = as.numeric(cv_res$delta_min),
      delta_max = as.numeric(cv_res$delta_max),
      delta_auc_ci_low = as.numeric(cv_res$delta_ci_low),
      delta_auc_ci_high = as.numeric(cv_res$delta_ci_high),
      delta_auc_oof = as.numeric(cv_res$delta_auc_oof),
      beta_adjusted = as.numeric(fit_res$beta_adj),
      or_adjusted = as.numeric(fit_res$or_adj),
      p_adjusted = as.numeric(fit_res$p_adj),
      lrt_p = as.numeric(fit_res$lrt_p),
      lrt_p_adj_bh = NA_real_,
      r_with_infiltration = as.numeric(r_inf),
      vif = as.numeric(vif_val),
      high_collinearity = isTRUE(high_collinearity),
      jaccard_overlap = as.numeric(jacc),
      converged = isTRUE(converged),
      separation_flag = isTRUE(separation_flag),
      analysis_status = status_str,
      adds_independent_value = FALSE,
      exploratory_independent = FALSE,
      stringsAsFactors = FALSE
    )
  }

  if (length(rows) == 0L) {
    empty_pc <- data.frame(
      cohort = character(), n = integer(), n_pos = integer(), n_neg = integer(),
      epv = numeric(), low_epv = logical(),
      auc_infiltration_only = numeric(), auc_combined = numeric(),
      delta_auc_combined = numeric(), sd_delta_auc = numeric(),
      delta_min = numeric(), delta_max = numeric(),
      delta_auc_ci_low = numeric(), delta_auc_ci_high = numeric(),
      delta_auc_oof = numeric(), beta_adjusted = numeric(),
      or_adjusted = numeric(), p_adjusted = numeric(),
      lrt_p = numeric(), lrt_p_adj_bh = numeric(),
      r_with_infiltration = numeric(), vif = numeric(),
      high_collinearity = logical(), jaccard_overlap = numeric(),
      converged = logical(), separation_flag = logical(),
      analysis_status = character(),
      adds_independent_value = logical(),
      exploratory_independent = logical(),
      stringsAsFactors = FALSE
    )
    empty_sm <- data.frame(
      infiltration_axis = as.character(infiltration_axis)[1],
      method = as.character(method)[1],
      n_cohorts = 0L,
      n_independent = 0L,
      combined_lrt_p = NA_real_,
      fisher_omnibus_p = NA_real_,
      stouffer_z = NA_real_,
      direction_conflict = FALSE,
      min_lrt_p_adj_bh = NA_real_,
      median_lrt_p = NA_real_,
      mean_delta_auc = NA_real_,
      stringsAsFactors = FALSE
    )
    return(list(per_cohort = empty_pc, summary = empty_sm))
  }

  pc <- do.call(rbind, rows)
  rownames(pc) <- NULL

  # 4.3: Multiple testing correction across evaluated cohorts BEFORE adds_independent_value decision
  ok_p <- is.finite(pc$lrt_p)
  if (any(ok_p)) {
    pc$lrt_p_adj_bh[ok_p] <- stats::p.adjust(pc$lrt_p[ok_p], method = "BH")
  }

  base_indep <- isTRUE_vec(
    pc$converged &
      !pc$separation_flag &
      !(pc$analysis_status %in% c("failed_collinearity", "failed_separation_or_convergence")) &
      is.finite(pc$lrt_p_adj_bh) & (pc$lrt_p_adj_bh < 0.05) &
      is.finite(pc$beta_adjusted) & (pc$beta_adjusted > 0) &
      is.finite(pc$delta_auc_combined) & (pc$delta_auc_combined > 0)
  )
  pc$adds_independent_value <- isTRUE_vec(
    base_indep & !pc$low_epv & !pc$high_collinearity & (pc$analysis_status == "ok")
  )
  pc$exploratory_independent <- isTRUE_vec(
    base_indep & (pc$low_epv | pc$high_collinearity)
  )

  # 4.4 & BUG-05: Directional Stouffer Z synthesis vs sign-blind Fisher omnibus p
  valid_meta <- is.finite(pc$lrt_p) & is.finite(pc$beta_adjusted) &
    pc$converged & !pc$separation_flag &
    !(pc$analysis_status %in% c("failed_collinearity", "failed_separation_or_convergence"))
  lrt_ok <- pc$lrt_p[valid_meta]
  bh_ok <- pc$lrt_p_adj_bh[valid_meta & is.finite(pc$lrt_p_adj_bh)]
  d_ok <- pc$delta_auc_combined[is.finite(pc$delta_auc_combined)]

  fisher_p <- if (length(lrt_ok) > 0L) {
    stat <- -2 * sum(log(pmax(lrt_ok, 1e-300)))
    stats::pchisq(stat, df = 2L * length(lrt_ok), lower.tail = FALSE)
  } else {
    NA_real_
  }

  stouffer_z <- NA_real_
  dir_p <- NA_real_
  dir_conflict <- FALSE
  if (any(valid_meta)) {
    p_m <- pmin(pmax(pc$lrt_p[valid_meta], 1e-300), 1 - 1e-16)
    b_m <- pc$beta_adjusted[valid_meta]
    n_m <- as.numeric(pc$n[valid_meta])
    signs <- ifelse(b_m >= 0, 1, -1)
    dir_conflict <- any(signs > 0) && any(signs < 0)
    # Magnitude from chi-square(1) LRT quantile, sign from adjusted beta
    z_mag <- sqrt(pmax(0, stats::qchisq(p_m, df = 1L, lower.tail = FALSE)))
    z_signed <- signs * z_mag
    w_m <- sqrt(pmax(n_m, 1))
    stouffer_z <- sum(w_m * z_signed) / sqrt(sum(w_m^2))
    # One-sided directional p-value testing positive independent value (beta > 0)
    dir_p <- stats::pnorm(stouffer_z, lower.tail = FALSE)
  } else if (length(lrt_ok) > 0L) {
    dir_p <- fisher_p
  }

  sm <- data.frame(
    infiltration_axis = as.character(infiltration_axis)[1],
    method = as.character(method)[1],
    n_cohorts = as.integer(nrow(pc)),
    n_independent = as.integer(sum(pc$adds_independent_value %in% TRUE, na.rm = TRUE)),
    combined_lrt_p = as.numeric(dir_p),
    fisher_omnibus_p = as.numeric(fisher_p),
    stouffer_z = as.numeric(stouffer_z),
    direction_conflict = isTRUE(dir_conflict),
    min_lrt_p_adj_bh = if (length(bh_ok) > 0L) min(bh_ok) else NA_real_,
    median_lrt_p = if (length(lrt_ok) > 0L) stats::median(lrt_ok) else NA_real_,
    mean_delta_auc = if (length(d_ok) > 0L) mean(d_ok) else NA_real_,
    stringsAsFactors = FALSE
  )

  list(per_cohort = pc, summary = sm)
}

#' @noRd
isTRUE_vec <- function(x) {
  !is.na(x) & as.logical(x)
}

#' Score signature and infiltration inside a CV fold without test-set leakage
#'
#' When `method %in% c("zscore", "rank")`, gene-level standardization or rank
#' scaling is learned strictly on the training columns (`tr_idx`) and applied
#' to the test columns (`te_idx`). Predictor standardization (`sig_score_z`,
#' `inf_score_z`) is always learned on `tr_idx` and applied to `te_idx`.
#'
#' @noRd
.score_fold_pair <- function(tr_idx,
                             te_idx,
                             df_glm,
                             expr_sub = NULL,
                             sig = NULL,
                             sig_dir = "<",
                             inf_genes = NULL,
                             method = "zscore",
                             eff_min_genes = 2L) {
  n_total <- nrow(df_glm)
  use_fold_gene_prep <- !is.null(expr_sub) && !is.null(sig) &&
    !is.null(inf_genes) && (method %in% c("zscore", "rank"))

  if (use_fold_gene_prep) {
    parsed_sig <- .parse_signature_genes(
      sig$genes, sig$genes_up, sig$genes_down, weights = sig$weights
    )
    cov_up <- if (length(parsed_sig$genes_up) > 0L) match_genes(expr_sub, parsed_sig$genes_up) else .empty_gene_coverage()
    cov_dn <- if (length(parsed_sig$genes_down) > 0L) match_genes(expr_sub, parsed_sig$genes_down) else .empty_gene_coverage()
    cov_inf <- match_genes(expr_sub, inf_genes)

    all_needed <- unique(c(cov_up$found, cov_dn$found, cov_inf$found))
    sub_mat <- expr_sub[all_needed, , drop = FALSE]
    tr_mat <- sub_mat[, tr_idx, drop = FALSE]
    te_mat <- sub_mat[, te_idx, drop = FALSE]

    if (identical(method, "zscore")) {
      tr_n <- rowSums(is.finite(tr_mat))
      tr_mean <- rowMeans(tr_mat, na.rm = TRUE)
      tr_sd <- sqrt(rowSums((tr_mat - tr_mean)^2, na.rm = TRUE) / pmax(tr_n - 1L, 1L))
      bad_sd <- !is.finite(tr_sd) | tr_sd == 0
      tr_scaled <- (tr_mat - tr_mean) / tr_sd
      te_scaled <- (te_mat - tr_mean) / tr_sd
      tr_scaled[bad_sd, ] <- NA_real_
      te_scaled[bad_sd, ] <- NA_real_
    } else {
      # method == "rank": learn empirical rank scaling on training fold only
      tr_scaled <- matrix(NA_real_, nrow = nrow(tr_mat), ncol = ncol(tr_mat), dimnames = dimnames(tr_mat))
      te_scaled <- matrix(NA_real_, nrow = nrow(te_mat), ncol = ncol(te_mat), dimnames = dimnames(te_mat))
      for (i in seq_len(nrow(tr_mat))) {
        x_tr <- tr_mat[i, ]
        m_tr <- sum(is.finite(x_tr))
        if (m_tr <= 1L) next
        r_tr <- rank(x_tr, ties.method = "average", na.last = "keep")
        tr_scaled[i, ] <- (r_tr - 1) / (m_tr - 1)
        x_tr_ok <- sort(x_tr[is.finite(x_tr)])
        x_te <- te_mat[i, ]
        ok_te <- is.finite(x_te)
        if (any(ok_te)) {
          # Position of each test observation relative to training distribution in [0, 1]
          less <- vapply(x_te[ok_te], function(v) sum(x_tr_ok < v), numeric(1))
          eq <- vapply(x_te[ok_te], function(v) sum(x_tr_ok == v), numeric(1))
          r_te_equiv <- 1 + less + 0.5 * pmax(0, eq - 1)
          te_scaled[i, ok_te] <- pmin(1, pmax(0, (r_te_equiv - 1) / (m_tr - 1)))
        }
      }
    }

    full_scaled <- matrix(NA_real_, nrow = nrow(sub_mat), ncol = n_total, dimnames = dimnames(sub_mat))
    full_scaled[, tr_idx] <- tr_scaled
    full_scaled[, te_idx] <- te_scaled

    # Aggregate signature and infiltration on fold-scaled matrix using "mean" (linear aggregation)
    if (isTRUE(parsed_sig$weighted)) {
      found_sig <- c(cov_up$found, cov_dn$found)
      w_found <- .resolve_found_weights(
        found_sig, parsed_sig$weights, c(cov_up$alias_rescued, cov_dn$alias_rescued)
      )
      raw_sig <- .score_weighted_on_genes(full_scaled, found_sig, w_found, method = "mean")
    } else {
      s_up <- if (length(cov_up$found) > 0L) colMeans(full_scaled[cov_up$found, , drop = FALSE], na.rm = TRUE) else NULL
      s_dn <- if (length(cov_dn$found) > 0L) colMeans(full_scaled[cov_dn$found, , drop = FALSE], na.rm = TRUE) else NULL
      raw_sig <- if (!is.null(s_up) && !is.null(s_dn)) s_up - s_dn else if (!is.null(s_up)) s_up else -s_dn
    }
    raw_sig <- orient_score(as.numeric(raw_sig), direction = sig_dir)
    raw_inf <- as.numeric(colMeans(full_scaled[cov_inf$found, , drop = FALSE], na.rm = TRUE))
  } else {
    raw_sig <- if ("x_sig_raw" %in% names(df_glm)) df_glm$x_sig_raw else df_glm$sig_score_z
    raw_inf <- if ("x_inf_raw" %in% names(df_glm)) df_glm$x_inf_raw else df_glm$inf_score_z
  }

  # Standardize predictors using ONLY training fold mean and SD
  mu_sig <- mean(raw_sig[tr_idx], na.rm = TRUE)
  sd_sig <- stats::sd(raw_sig[tr_idx], na.rm = TRUE)
  mu_inf <- mean(raw_inf[tr_idx], na.rm = TRUE)
  sd_inf <- stats::sd(raw_inf[tr_idx], na.rm = TRUE)
  if (!is.finite(sd_sig) || sd_sig == 0 || !is.finite(sd_inf) || sd_inf == 0) {
    return(NULL)
  }

  sig_z <- (raw_sig - mu_sig) / sd_sig
  inf_z <- (raw_inf - mu_inf) / sd_inf
  list(
    tr_df = data.frame(
      response = df_glm$response[tr_idx],
      sig_score_z = sig_z[tr_idx],
      inf_score_z = inf_z[tr_idx]
    ),
    te_df = data.frame(
      response = df_glm$response[te_idx],
      sig_score_z = sig_z[te_idx],
      inf_score_z = inf_z[te_idx]
    )
  )
}

#' Repeated stratified K-fold CV AUC (fold-level paired summary + OOF diagnostic)
#'
#' @noRd
.cv_glm_auc <- function(df_glm,
                        n_folds = 5L,
                        n_repeats = 5L,
                        seed = 42L,
                        expr_sub = NULL,
                        sig = NULL,
                        sig_dir = "<",
                        inf_genes = NULL,
                        method = "zscore",
                        eff_min_genes = 2L) {
  pos_idx <- which(df_glm$response == 1L)
  neg_idx <- which(df_glm$response == 0L)
  k_folds <- min(as.integer(n_folds), length(pos_idx), length(neg_idx))
  if (k_folds < 2L || n_repeats < 1L) {
    return(list(
      auc_inf = NA_real_, auc_comb = NA_real_, delta_auc = NA_real_,
      sd_delta_auc = NA_real_, delta_ci_low = NA_real_, delta_ci_high = NA_real_,
      delta_auc_oof = NA_real_, n_valid_folds = 0L
    ))
  }

  n <- nrow(df_glm)
  run_cv <- function() {
    rep_auc_b <- rep(NA_real_, n_repeats)
    rep_auc_f <- rep(NA_real_, n_repeats)
    rep_delta <- rep(NA_real_, n_repeats)
    oof_auc_b <- rep(NA_real_, n_repeats)
    oof_auc_f <- rep(NA_real_, n_repeats)
    fold_deltas_all <- numeric(0)

    for (r in seq_len(n_repeats)) {
      p_shuf <- sample(pos_idx, length(pos_idx))
      n_shuf <- sample(neg_idx, length(neg_idx))
      folds <- integer(n)
      folds[p_shuf] <- rep(seq_len(k_folds), length.out = length(p_shuf))
      folds[n_shuf] <- rep(seq_len(k_folds), length.out = length(n_shuf))

      oof_b <- rep(NA_real_, n)
      oof_f <- rep(NA_real_, n)
      f_auc_b <- rep(NA_real_, k_folds)
      f_auc_f <- rep(NA_real_, k_folds)
      f_weights <- rep(0, k_folds)

      for (f in seq_len(k_folds)) {
        tr_idx <- which(folds != f)
        te_idx <- which(folds == f)
        fold_data <- .score_fold_pair(
          tr_idx = tr_idx,
          te_idx = te_idx,
          df_glm = df_glm,
          expr_sub = expr_sub,
          sig = sig,
          sig_dir = sig_dir,
          inf_genes = inf_genes,
          method = method,
          eff_min_genes = eff_min_genes
        )
        if (is.null(fold_data)) {
          next
        }
        tr_df <- fold_data$tr_df
        te_df <- fold_data$te_df
        if (length(unique(tr_df$response)) < 2L) {
          next
        }
        mb <- tryCatch(
          suppressWarnings(stats::glm(response ~ inf_score_z, data = tr_df, family = stats::binomial())),
          error = function(e) NULL
        )
        mf <- tryCatch(
          suppressWarnings(stats::glm(response ~ sig_score_z + inf_score_z, data = tr_df, family = stats::binomial())),
          error = function(e) NULL
        )
        if (!is.null(mb) && !is.null(mf)) {
          pb_te <- as.numeric(stats::predict(mb, newdata = te_df, type = "response"))
          pf_te <- as.numeric(stats::predict(mf, newdata = te_df, type = "response"))
          oof_b[te_idx] <- pb_te
          oof_f[te_idx] <- pf_te

          y_te <- te_df$response
          n1_f <- sum(y_te == 1L)
          n0_f <- sum(y_te == 0L)
          if (n1_f >= 1L && n0_f >= 1L && all(is.finite(pb_te)) && all(is.finite(pf_te))) {
            f_auc_b[f] <- auc_mw(pb_te, y_te)
            f_auc_f[f] <- auc_mw(pf_te, y_te)
            f_weights[f] <- as.numeric(n1_f * n0_f)
          }
        }
      }

      ok_f <- is.finite(f_auc_b) & is.finite(f_auc_f) & (f_weights > 0)
      if (any(ok_f)) {
        w_f <- f_weights[ok_f] / sum(f_weights[ok_f])
        rep_auc_b[r] <- sum(w_f * f_auc_b[ok_f])
        rep_auc_f[r] <- sum(w_f * f_auc_f[ok_f])
        rep_delta[r] <- rep_auc_f[r] - rep_auc_b[r]
        fold_deltas_all <- c(fold_deltas_all, f_auc_f[ok_f] - f_auc_b[ok_f])
      }
      if (all(is.finite(oof_b))) {
        oof_auc_b[r] <- auc_mw(oof_b, df_glm$response)
      }
      if (all(is.finite(oof_f))) {
        oof_auc_f[r] <- auc_mw(oof_f, df_glm$response)
      }
    }

    auc_inf_oof <- if (any(is.finite(oof_auc_b))) mean(oof_auc_b, na.rm = TRUE) else NA_real_
    auc_comb_oof <- if (any(is.finite(oof_auc_f))) mean(oof_auc_f, na.rm = TRUE) else NA_real_
    delta_oof <- if (is.finite(auc_inf_oof) && is.finite(auc_comb_oof)) auc_comb_oof - auc_inf_oof else NA_real_

    auc_inf <- if (any(is.finite(rep_auc_b))) mean(rep_auc_b, na.rm = TRUE) else auc_inf_oof
    auc_comb <- if (any(is.finite(rep_auc_f))) mean(rep_auc_f, na.rm = TRUE) else auc_comb_oof
    delta_auc <- if (is.finite(auc_inf) && is.finite(auc_comb)) auc_comb - auc_inf else delta_oof

    delta_samples <- if (sum(is.finite(rep_delta)) >= 2L) {
      rep_delta[is.finite(rep_delta)]
    } else {
      fold_deltas_all[is.finite(fold_deltas_all)]
    }
    sd_delta <- if (length(delta_samples) >= 2L) stats::sd(delta_samples) else NA_real_
    d_min <- if (length(delta_samples) >= 1L) min(delta_samples) else NA_real_
    d_max <- if (length(delta_samples) >= 1L) max(delta_samples) else NA_real_
    ci_low <- if (length(delta_samples) >= 20L) {
      as.numeric(stats::quantile(delta_samples, probs = 0.025, names = FALSE))
    } else {
      NA_real_
    }
    ci_high <- if (length(delta_samples) >= 20L) {
      as.numeric(stats::quantile(delta_samples, probs = 0.975, names = FALSE))
    } else {
      NA_real_
    }

    list(
      auc_inf = auc_inf,
      auc_comb = auc_comb,
      delta_auc = delta_auc,
      sd_delta_auc = sd_delta,
      delta_min = d_min,
      delta_max = d_max,
      delta_ci_low = ci_low,
      delta_ci_high = ci_high,
      delta_auc_oof = delta_oof,
      n_valid_folds = as.integer(length(fold_deltas_all))
    )
  }

  if (is.null(seed) || is.na(seed)) {
    run_cv()
  } else {
    withr::with_seed(as.integer(seed), run_cv())
  }
}

#' @noRd
.resolve_infiltration_genes <- function(infiltration_axis, infiltration_genes = NULL) {
  if (!is.null(infiltration_genes) && length(infiltration_genes) > 0L) {
    g <- unique(as.character(infiltration_genes))
    return(g[!is.na(g) & nzchar(g)])
  }
  axes <- .default_reference_axes()
  ax <- as.character(infiltration_axis)[1]
  if (identical(ax, "CD8_T_Cell") || identical(ax, "T_Cell_Infiltration")) {
    return(axes$T_Cell_Infiltration$genes)
  }
  if (ax %in% names(axes)) {
    return(axes[[ax]]$genes)
  }
  warning(
    "Unknown `infiltration_axis` '", ax,
    "'; falling back to T_Cell_Infiltration (CD3D, CD3E, CD8A, CD8B).",
    call. = FALSE
  )
  axes$T_Cell_Infiltration$genes
}
