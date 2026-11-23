#' Imza x veri seti benchmark'ini calistir
#'
#' Her veri seti icinde ayri AUC hesaplar; ham ifadeleri birlestirmez.
#' Farkli normalizasyonlar yan yana gosterilirken uyari uretir.
#'
#' @param signatures Isimli liste; her eleman `list(genes = <chr>, direction = "<"|">",
#'   derived_from = <chr>)` veya isaretli `list(genes_up=..., genes_down=..., ...)`.
#'   `genes` isimli sayisal isaret vektoru de olabilir. `direction` ve
#'   `derived_from` istege bagli.
#'   `derived_from = NULL` (belirtilmemis) -> `independent = NA` (fail-closed).
#'   `derived_from = character(0)` -> bagimsiz (`TRUE`).
#' @param datasets Isimli liste; her eleman `validate_dataset()` yapisinda.
#' @param methods Skorlama yontemleri: `"mean"`, `"zscore"`, `"median"`, `"rank"`,
#'   `"ss_rank"`. Varsayilan dortlu liste geriye donuk uyumluluk icin korunur.
#' @param direction Varsayilan imza yonu (imza duzeyinde yoksa).
#' @param min_genes `score_signature()` esigi.
#' @param coverage_threshold Kapsama uyari esigi (0-1). Varsayilan 0.70.
#' @param min_coverage `coverage_status = "ok"` alt esigi (0-1). Varsayilan 0.50.
#' @param invalid_coverage `coverage_status = "invalid"` ust esigi (0-1). Varsayilan 0.20.
#' @param apply_log2 `scale = "linear"` iken log2 uygula.
#' @param sd_handling `encode_response()` politikasi.
#' @param min_n,min_per_class Kucuk kohort esikleri (`roc_auc`).
#' @param compute_pr PR-AUC hesapla. Varsayilan TRUE.
#' @param n_boot PR-AUC bootstrap tekrari. Varsayilan 1000; hizli onizleme icin dusurulebilir.
#' @param compute_null Null permutasyon testi. Varsayilan FALSE.
#' @param n_perm Null permutasyon sayisi. Varsayilan 1000.
#' @param seed Rastgelelik tohumu (PR bootstrap ve null); hucre basina turetilir.
#' @param match_expression Null icin ifade-ve-dispersiyon (2D ortalama x SD)
#'   eslesmeli ornekleme.
#'
#' @return `data.frame` sonuc tablosu:
#'   - ROC-AUC: `auc`, `auc_ci_low`, `auc_ci_high`, `auc_var`, `var_zero`, `small_cohort`
#'   - PR-AUC ve tanilari: `pr_auc`, `pr_baseline`, `pr_ci_low`, `pr_ci_high`,
#'     `tie_fraction`, `small_positive_class`, `n_boot_valid`, `ci_unreliable`
#'   - Null permutasyon ve FDR aileleri: `null_p` (ham), `null_p_adj_bh`
#'     (yontem-ici BH), `null_p_adj_bh_global` (tum yontem x imza x kohort
#'     ailesinde global BH), `background_type`, `universe_size`
#'   - Kapsama ve metaveri: `coverage_pct`, `coverage_status`, `low_coverage`,
#'     `n_alias_rescued`, `signed`, `treatment`, `endpoint_type`,
#'     `cancer_type`, `normalization`, `power_auc70`, `independent`,
#'     `skipped`, `skip_reason`.
#' @examples
#' data(datasets_demo, package = "SigBenchR")
#' data(signatures_demo, package = "SigBenchR")
#' res <- run_benchmark(
#'   signatures = signatures_demo["IFNG6"],
#'   datasets = datasets_demo["demo_A"],
#'   methods = c("mean", "zscore"),
#'   compute_pr = FALSE,
#'   compute_null = FALSE
#' )
#' head(res[, c("dataset", "signature", "method", "auc", "auc_ci_low", "auc_ci_high")])
#' @export
run_benchmark <- function(signatures,
                          datasets,
                          methods = c("mean", "zscore", "median", "rank"),
                          direction = "<",
                          min_genes = 3L,
                          coverage_threshold = 0.70,
                          min_coverage = coverage_threshold,
                          invalid_coverage = 0.20,
                          apply_log2 = TRUE,
                          sd_handling = "exclude",
                          min_n = 20L,
                          min_per_class = 5L,
                          compute_pr = TRUE,
                          n_boot = 1000L,
                          compute_null = FALSE,
                          n_perm = 1000L,
                          seed = 42L,
                          match_expression = FALSE) {
  if (is.null(names(signatures)) || any(!nzchar(names(signatures)))) {
    stop("`signatures` must be a named list.", call. = FALSE)
  }
  if (is.null(names(datasets)) || any(!nzchar(names(datasets)))) {
    stop("`datasets` must be a named list.", call. = FALSE)
  }
  method_choices <- c("mean", "zscore", "median", "rank", "ss_rank")
  methods <- match.arg(methods, choices = method_choices, several.ok = TRUE)

  norms <- vapply(datasets, function(d) {
    if (is.list(d) && is.list(d$info) && !is.null(d$info$normalization) &&
        length(d$info$normalization) > 0L && !is.na(d$info$normalization[1]) &&
        nzchar(as.character(d$info$normalization)[1])) {
      as.character(d$info$normalization)[1]
    } else {
      "unknown"
    }
  }, character(1))
  valid_norms <- unique(norms[norms != "unknown"])
  if (length(valid_norms) > 1L) {
    warning(
      "Datasets have heterogeneous normalizations: ",
      paste(valid_norms, collapse = ", "),
      ". AUCs are computed separately; raw expression is never pooled.",
      call. = FALSE
    )
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
        "Cohort '", ds_name, "' failed validation/preparation and was marked skipped: ",
        conditionMessage(prep),
        call. = FALSE
      )
      for (sig_name in names(signatures)) {
        sig <- tryCatch(.as_signature(signatures[[sig_name]], name = sig_name), error = function(e) list())
        indep_val <- if (is.null(sig$derived_from)) {
          NA
        } else {
          !(tolower(ds_name) %in% tolower(as.character(sig$derived_from)))
        }
        for (method in methods) {
          row_i <- row_i + 1L
          rows[[row_i]] <- list(
            dataset = ds_name,
            signature = sig_name,
            method = method,
            n = NA_integer_,
            n_pos = NA_integer_,
            n_neg = NA_integer_,
            coverage_pct = 0,
            n_found = 0L,
            n_total = NA_integer_,
            coverage_label = "0/0 (0%)",
            found_genes = "",
            missing_genes = "",
            low_coverage = TRUE,
            coverage_status = "invalid",
            n_alias_rescued = 0L,
            auc = NA_real_,
            auc_ci_low = NA_real_,
            auc_ci_high = NA_real_,
            auc_var = NA_real_,
            var_zero = NA,
            small_cohort = NA,
            pr_auc = NA_real_,
            pr_baseline = NA_real_,
            pr_ci_low = NA_real_,
            pr_ci_high = NA_real_,
            tie_fraction = NA_real_,
            small_positive_class = NA,
            n_boot_valid = NA_integer_,
            ci_unreliable = NA,
            null_p = NA_real_,
            background_type = NA_character_,
            universe_size = NA_integer_,
            power_auc70 = NA_real_,
            normalization = "unknown",
            scale = NA_character_,
            cancer = NA_character_,
            treatment = NA_character_,
            endpoint_type = NA_character_,
            cancer_type = NA_character_,
            signed = FALSE,
            n_up_found = 0L,
            n_down_found = 0L,
            n_all_na_dropped = 0L,
            n_excluded_sd = NA_integer_,
            skipped = TRUE,
            skip_reason = paste0("Cohort preparation failed: ", conditionMessage(prep)),
            independent = indep_val
          )
        }
      }
      next
    }
    ds <- prep$ds
    enc <- prep$enc
    expr_full <- prep$expr
    y_kept <- prep$y_kept
    y_full <- enc$response
    names(y_full) <- colnames(expr_full)

    cancer <- if (!is.null(ds$sample_meta$cancer)) {
      paste(unique(as.character(ds$sample_meta$cancer)), collapse = ";")
    } else {
      NA_character_
    }
    treatment <- if (!is.null(ds$info$treatment) &&
                     length(ds$info$treatment) > 0L &&
                     !is.na(ds$info$treatment[1])) {
      as.character(ds$info$treatment)[1]
    } else if (!is.null(ds$sample_meta$treatment)) {
      paste(unique(as.character(ds$sample_meta$treatment)), collapse = ";")
    } else {
      NA_character_
    }
    endpoint_type <- if (!is.null(ds$info$endpoint_type) &&
                         length(ds$info$endpoint_type) > 0L &&
                         !is.na(ds$info$endpoint_type[1])) {
      as.character(ds$info$endpoint_type)[1]
    } else {
      NA_character_
    }
    cancer_type <- if (!is.null(ds$info$cancer_type) &&
                       length(ds$info$cancer_type) > 0L &&
                       !is.na(ds$info$cancer_type[1])) {
      as.character(ds$info$cancer_type)[1]
    } else {
      NA_character_
    }

    for (sig_name in names(signatures)) {
      sig <- .as_signature(signatures[[sig_name]], name = sig_name)
      sig_dir <- .normalize_roc_direction(sig$direction, default = direction)

      for (method in methods) {
        cell_seed <- .cell_seed(seed, ds_name, sig_name, method)
        eff_min_genes <- .resolve_min_genes(sig, default = min_genes)

        scored <- score_signature(
          expr_full,
          genes = sig$genes,
          genes_up = sig$genes_up,
          genes_down = sig$genes_down,
          weights = sig$weights,
          method = method,
          min_genes = eff_min_genes
        )

        cov_cls <- .classify_coverage(
          scored$coverage$coverage_pct,
          ok_threshold = min_coverage,
          invalid_threshold = invalid_coverage
        )
        cov_frac <- scored$coverage$coverage_pct / 100
        low_cov <- cov_cls$low_coverage
        cov_status <- cov_cls$coverage_status
        n_alias <- as.integer(scored$coverage$n_alias_rescued)
        if (length(n_alias) != 1L || is.na(n_alias)) {
          n_alias <- 0L
        }

        if (identical(cov_status, "invalid")) {
          warning(
            "Invalid gene coverage for cohort '", ds_name, "' / signature '",
            sig_name, "': coverage=",
            sprintf("%.1f%%", 100 * cov_frac),
            " (< ", sprintf("%.0f%%", 100 * invalid_coverage),
            "). For a full gene matrix use read_user_dataset() or ",
            "load_cohort_bundle().",
            call. = FALSE
          )
        }

        base <- list(
          dataset = ds_name,
          signature = sig_name,
          method = method,
          n = NA_integer_,
          n_pos = NA_integer_,
          n_neg = NA_integer_,
          coverage_pct = scored$coverage$coverage_pct,
          n_found = scored$coverage$n_found,
          n_total = scored$coverage$n_total,
          coverage_label = sprintf(
            "%d/%d (%.0f%%)",
            scored$coverage$n_found,
            scored$coverage$n_total,
            scored$coverage$coverage_pct
          ),
          found_genes = paste(scored$coverage$found, collapse = ","),
          missing_genes = if (length(scored$coverage$missing) == 0L) {
            ""
          } else {
            paste(scored$coverage$missing, collapse = ",")
          },
          low_coverage = low_cov,
          coverage_status = cov_status,
          n_alias_rescued = n_alias,
          auc = NA_real_,
          auc_ci_low = NA_real_,
          auc_ci_high = NA_real_,
          auc_var = NA_real_,
          var_zero = NA,
          small_cohort = NA,
          pr_auc = NA_real_,
          pr_baseline = NA_real_,
          pr_ci_low = NA_real_,
          pr_ci_high = NA_real_,
          tie_fraction = NA_real_,
          small_positive_class = NA,
          n_boot_valid = NA_integer_,
          ci_unreliable = NA,
          null_p = NA_real_,
          background_type = NA_character_,
          universe_size = NA_integer_,
          power_auc70 = NA_real_,
          normalization = if (!is.null(ds$info$normalization) &&
                              length(ds$info$normalization) > 0L &&
                              !is.na(ds$info$normalization[1]) &&
                              nzchar(as.character(ds$info$normalization)[1])) {
            as.character(ds$info$normalization)[1]
          } else {
            "unknown"
          },
          scale = ds$info$scale,
          cancer = cancer,
          treatment = treatment,
          endpoint_type = endpoint_type,
          cancer_type = cancer_type,
          signed = isTRUE(scored$signed),
          n_up_found = scored$n_up_found,
          n_down_found = scored$n_down_found,
          n_all_na_dropped = scored$n_all_na_dropped,
          n_excluded_sd = enc$n_excluded_sd,
          skipped = scored$skipped,
          skip_reason = scored$skip_reason,
          independent = NA
        )

        # P0-4: NULL derived_from -> NA (fail-closed); karakter vektor -> case-insensitive
        if (is.null(sig$derived_from)) {
          base$independent <- NA
        } else {
          derived <- as.character(sig$derived_from)
          base$independent <- !(tolower(ds_name) %in% tolower(derived))
        }

        if (isTRUE(scored$skipped)) {
          row_i <- row_i + 1L
          rows[[row_i]] <- base
          next
        }

        sample_ids <- intersect(names(scored$scores), names(y_kept))
        resp <- y_kept[sample_ids]
        sc_kept <- scored$scores[sample_ids]
        keep <- is.finite(sc_kept) & !is.na(resp)

        scores_use <- sc_kept[keep]
        resp_use <- resp[keep]

        if (length(unique(resp_use)) < 2L) {
          base$skipped <- TRUE
          base$skip_reason <- "Single response class remaining (ROC cannot be computed)."
          row_i <- row_i + 1L
          rows[[row_i]] <- base
          next
        }

        roc_res <- tryCatch(
          roc_auc(
            scores_use, resp_use,
            direction = sig_dir,
            min_n = min_n,
            min_per_class = min_per_class
          ),
          error = function(e) e
        )

        if (inherits(roc_res, "error")) {
          base$skipped <- TRUE
          base$skip_reason <- conditionMessage(roc_res)
          row_i <- row_i + 1L
          rows[[row_i]] <- base
          next
        }

        base$n <- roc_res$n
        base$n_pos <- roc_res$n_pos
        base$n_neg <- roc_res$n_neg
        base$auc <- roc_res$auc
        base$auc_ci_low <- roc_res$auc_ci_low
        base$auc_ci_high <- roc_res$auc_ci_high
        base$auc_var <- roc_res$auc_var
        base$var_zero <- roc_res$var_zero
        base$small_cohort <- roc_res$small_cohort
        base$skipped <- FALSE
        base$skip_reason <- NA_character_
        base$power_auc70 <- suppressWarnings(
          estimate_auc_power(
            n_pos = roc_res$n_pos,
            n_neg = roc_res$n_neg,
            target_auc = 0.70,
            null_auc = 0.50,
            alpha = 0.05
          )$power
        )

        if (isTRUE(compute_pr)) {
          pr_res <- tryCatch(
            pr_auc(
              scores_use, resp_use,
              direction = sig_dir,
              n_boot = n_boot,
              seed = cell_seed
            ),
            error = function(e) e
          )
          if (inherits(pr_res, "error")) {
            warning(
              "pr_auc failed (", ds_name, " / ", sig_name, " / ", method, "): ",
              conditionMessage(pr_res),
              call. = FALSE
            )
          } else {
            base$pr_auc <- pr_res$pr_auc
            base$pr_baseline <- pr_res$pr_baseline
            base$pr_ci_low <- pr_res$pr_ci_low
            base$pr_ci_high <- pr_res$pr_ci_high
            base$tie_fraction <- pr_res$tie_fraction
            base$small_positive_class <- pr_res$small_positive_class
            base$n_boot_valid <- pr_res$n_boot_valid
            base$ci_unreliable <- pr_res$ci_unreliable
          }
        }

        if (isTRUE(compute_null)) {
          n_up_f <- as.integer(scored$n_up_found)
          n_dn_f <- as.integer(scored$n_down_found)
          found_all <- scored$coverage$found
          sig_up <- if (n_up_f > 0L) {
            found_all[seq_len(n_up_f)]
          } else {
            character(0)
          }
          sig_dn <- if (n_dn_f > 0L) {
            found_all[seq.int(n_up_f + 1L, length.out = n_dn_f)]
          } else {
            character(0)
          }
          null_res <- tryCatch(
            null_permutation_test(
              mat = expr_full,
              response = y_full,
              n_genes = n_up_f,
              n_genes_down = n_dn_f,
              exclude_genes = found_all,
              method = method,
              direction = sig_dir,
              n_perm = n_perm,
              seed = cell_seed,
              match_expression = match_expression,
              signature_genes = if (length(sig_up) > 0L) sig_up else NULL,
              signature_genes_down = if (length(sig_dn) > 0L) sig_dn else NULL,
              observed_scores = scores_use,
              min_genes = min_genes,
              weights = sig$weights
            ),
            error = function(e) e
          )
          if (inherits(null_res, "error")) {
            warning(
              "null_permutation_test failed (", ds_name, " / ", sig_name,
              " / ", method, "): ", conditionMessage(null_res),
              call. = FALSE
            )
          } else {
            base$null_p <- null_res$null_p
            base$background_type <- null_res$background_type
            base$universe_size <- null_res$universe_size
          }
        }

        row_i <- row_i + 1L
        rows[[row_i]] <- base
      }
    }
  }

  if (length(rows) == 0L) {
    return(data.frame())
  }

  out <- do.call(rbind.data.frame, lapply(rows, as.data.frame, stringsAsFactors = FALSE))
  rownames(out) <- NULL

  # Faz 2C & Faz 4.7: null_p icin BH (method-ici null_p_adj_bh + global null_p_adj_bh_global)
  out$null_p_adj_bh <- NA_real_
  out$null_p_adj_bh_global <- NA_real_
  if (isTRUE(compute_null) && "null_p" %in% names(out)) {
    out <- adjust_null_p(out, p_col = "null_p")
  }

  out
}

#' Map signature direction metadata to ROC direction
#'
#' Accepts legacy ROC codes (`"<"`, `">"`) and panel metadata
#' (`"positive"`, `"negative"`, `"signed"`).
#'
#' @param sig_direction Direction metadata from a signature definition (`"<"`,
#'   `">"`, `"positive"`, `"negative"`, `"signed"`, or `NULL`).
#' @param default Fallback ROC direction (`"<"` or `">"`).
#'
#' @return Normalized ROC direction string (`"<"` or `">"`).
#' @keywords internal
#' @export
normalize_roc_direction <- function(sig_direction, default = "<") {
  if (is.null(sig_direction) || length(sig_direction) == 0L ||
      is.na(sig_direction[1]) || !nzchar(as.character(sig_direction[1]))) {
    return(as.character(default)[1])
  }
  d <- as.character(sig_direction)[1]
  if (d %in% c("<", ">")) {
    return(d)
  }
  if (identical(d, "negative")) {
    return(">")
  }
  if (d %in% c("positive", "signed")) {
    return("<")
  }
  as.character(default)[1]
}

.normalize_roc_direction <- normalize_roc_direction

#' Deterministic order-independent per-cell RNG seed
#'
#' @noRd
.cell_seed <- function(seed, ds_name, sig_name, method) {
  if (is.null(seed) || length(seed) == 0L || is.na(seed[1])) {
    return(NULL)
  }
  key <- paste(
    as.character(ds_name)[1],
    as.character(sig_name)[1],
    as.character(method)[1],
    sep = "\r"
  )
  bytes <- as.integer(charToRaw(enc2utf8(key)))
  h <- 0
  mod <- 1073741823 # 2^30 - 1 (31-bit signed int safe)
  for (b in bytes) {
    h <- (h * 131 + b) %% mod
  }
  base_seed <- abs(as.numeric(seed[1])) %% mod
  as.integer((base_seed + h) %% mod)
}

#' Centralized coverage classification
#'
#' Classifies coverage percentages (0-100) into `"invalid"`, `"low_power"`,
#' or `"ok"` using a single `ok_threshold` (0-1) so heatmap, KPI, and
#' `summarize_auc()` never disagree.
#'
#' @noRd
.classify_coverage <- function(coverage_pct,
                               ok_threshold = 0.70,
                               invalid_threshold = 0.20) {
  cov_frac <- as.numeric(coverage_pct) / 100
  n <- length(cov_frac)
  cov_status <- rep("ok", n)
  bad <- !is.finite(cov_frac) | cov_frac < invalid_threshold
  low <- !bad & (cov_frac < ok_threshold)
  cov_status[bad] <- "invalid"
  cov_status[low] <- "low_power"
  list(
    coverage_status = cov_status,
    low_coverage = cov_status != "ok"
  )
}
