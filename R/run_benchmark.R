#' Imza x veri seti benchmark'ini calistir
#'
#' Her veri seti icinde ayri AUC hesaplar; ham ifadeleri birlestirmez.
#' Farkli normalizasyonlar yan yana gosterilirken uyari uretir.
#'
#' @param signatures Isimli liste; her eleman `list(genes = <chr>, direction = "<"|">",
#'   derived_from = <chr>)`. `direction` ve `derived_from` istege bagli.
#'   `derived_from = NULL` (belirtilmemis) -> `independent = NA` (fail-closed).
#'   `derived_from = character(0)` -> bagimsiz (`TRUE`).
#' @param datasets Isimli liste; her eleman `validate_dataset()` yapisinda.
#' @param methods Skorlama yontemleri: `"mean"`, `"zscore"`, `"median"`, `"rank"`.
#' @param direction Varsayilan imza yonu (imza duzeyinde yoksa).
#' @param min_genes `score_signature()` esigi.
#' @param coverage_threshold Kapsama uyari esigi (0-1). Varsayilan 0.70.
#' @param apply_log2 `scale = "linear"` iken log2 uygula.
#' @param sd_handling `encode_response()` politikasi.
#' @param min_n,min_per_class Kucuk kohort esikleri (`roc_auc`).
#' @param compute_pr PR-AUC hesapla. Varsayilan TRUE.
#' @param n_boot PR-AUC bootstrap tekrari. Varsayilan 1000; hizli onizleme icin dusurulebilir.
#' @param compute_null Null permutasyon testi. Varsayilan FALSE.
#' @param n_perm Null permutasyon sayisi. Varsayilan 1000.
#' @param seed Rastgelelik tohumu (PR bootstrap ve null); hucre basina turetilir.
#' @param match_expression Null icin ifade-eslesmeli ornekleme.
#'
#' @return `data.frame` sonuc tablosu.
#' @export
run_benchmark <- function(signatures,
                          datasets,
                          methods = c("mean", "zscore", "median", "rank"),
                          direction = "<",
                          min_genes = 3L,
                          coverage_threshold = 0.70,
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
    stop("`signatures` isimli bir liste olmalidir.", call. = FALSE)
  }
  if (is.null(names(datasets)) || any(!nzchar(names(datasets)))) {
    stop("`datasets` isimli bir liste olmalidir.", call. = FALSE)
  }
  method_choices <- c("mean", "zscore", "median", "rank")
  methods <- match.arg(methods, choices = method_choices, several.ok = TRUE)

  norms <- vapply(datasets, function(d) as.character(d$info$normalization), character(1))
  if (length(unique(norms)) > 1L) {
    warning(
      "Veri setleri farkli normalizasyonlara sahip: ",
      paste(unique(norms), collapse = ", "),
      ". AUC'ler ayri hesaplanir; ham ifadeler birlestirilmez.",
      call. = FALSE
    )
  }

  rows <- list()
  row_i <- 0L
  cell_idx <- 0L

  for (ds_name in names(datasets)) {
    ds <- validate_dataset(datasets[[ds_name]], name = ds_name)
    expr <- collapse_duplicate_genes(as.matrix(ds$expr))
    expr <- prepare_expression(
      expr,
      scale = ds$info$scale,
      apply_log2 = apply_log2
    )

    enc <- encode_response(
      ds$sample_meta$response_label,
      sd_handling = sd_handling
    )

    # P0-3: SD (veya NA yanit) orneklerini skorlama / null oncesi cikar
    expr_kept <- expr[, enc$kept, drop = FALSE]
    y_kept <- enc$response[enc$kept]
    names(y_kept) <- colnames(expr_kept)

    cancer <- if (!is.null(ds$sample_meta$cancer)) {
      paste(unique(as.character(ds$sample_meta$cancer)), collapse = ";")
    } else {
      NA_character_
    }
    treatment <- if (!is.null(ds$sample_meta$treatment)) {
      paste(unique(as.character(ds$sample_meta$treatment)), collapse = ";")
    } else {
      NA_character_
    }

    for (sig_name in names(signatures)) {
      sig <- signatures[[sig_name]]
      if (is.character(sig)) {
        sig <- list(genes = sig)
      }
      if (is.null(sig$genes)) {
        stop("Imza '", sig_name, "' icin `genes` zorunludur.", call. = FALSE)
      }
      sig_dir <- if (!is.null(sig$direction)) sig$direction else direction

      for (method in methods) {
        cell_idx <- cell_idx + 1L
        cell_seed <- if (is.null(seed)) {
          NULL
        } else {
          as.integer(seed) + as.integer(cell_idx)
        }

        scored <- score_signature(
          expr_kept,
          genes = sig$genes,
          method = method,
          min_genes = min_genes
        )

        cov_frac <- scored$coverage$coverage_pct / 100
        low_cov <- cov_frac < coverage_threshold

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
          null_p = NA_real_,
          normalization = ds$info$normalization,
          scale = ds$info$scale,
          cancer = cancer,
          treatment = treatment,
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

        sample_ids <- names(scored$scores)
        resp <- y_kept[sample_ids]
        keep <- is.finite(scored$scores) & !is.na(resp)

        scores_use <- scored$scores[keep]
        resp_use <- resp[keep]

        if (length(unique(resp_use)) < 2L) {
          base$skipped <- TRUE
          base$skip_reason <- "Yanitta tek sinif kaldi (ROC hesaplanamaz)."
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
              "pr_auc basarisiz (", ds_name, " / ", sig_name, " / ", method, "): ",
              conditionMessage(pr_res),
              call. = FALSE
            )
          } else {
            base$pr_auc <- pr_res$pr_auc
            base$pr_baseline <- pr_res$pr_baseline
            base$pr_ci_low <- pr_res$pr_ci_low
            base$pr_ci_high <- pr_res$pr_ci_high
          }
        }

        if (isTRUE(compute_null)) {
          null_res <- tryCatch(
            null_permutation_test(
              mat = expr_kept,
              response = y_kept,
              n_genes = scored$coverage$n_found,
              exclude_genes = scored$coverage$found,
              method = method,
              direction = sig_dir,
              n_perm = n_perm,
              seed = cell_seed,
              match_expression = match_expression,
              signature_genes = scored$coverage$found,
              observed_scores = scores_use,
              min_genes = min_genes
            ),
            error = function(e) e
          )
          if (inherits(null_res, "error")) {
            warning(
              "null_permutation_test basarisiz (", ds_name, " / ", sig_name,
              " / ", method, "): ", conditionMessage(null_res),
              call. = FALSE
            )
          } else {
            base$null_p <- null_res$null_p
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

  # Faz 2C: null_p icin BH (aile = method icinde imza x veri seti)
  out$null_p_adj_bh <- NA_real_
  if (isTRUE(compute_null) && "null_p" %in% names(out)) {
    out <- adjust_null_p(out, p_col = "null_p")
  }

  out
}
