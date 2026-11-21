#' Imza x veri seti benchmark'ini calistir
#'
#' Her veri seti icinde ayri AUC hesaplar; ham ifadeleri birlestirmez.
#' Farkli normalizasyonlar yan yana gosterilirken uyari uretir.
#'
#' @param signatures Isimli liste; her eleman `list(genes = <chr>, direction = "<"|">",
#'   derived_from = <chr>)`. `direction` ve `derived_from` istege bagli.
#' @param datasets Isimli liste; her eleman `validate_dataset()` yapisinda.
#' @param methods Skorlama yontemleri: `"mean"`, `"zscore"`, `"median"`, `"rank"`.
#' @param direction Varsayilan imza yonu (imza duzeyinde yoksa).
#' @param min_genes `score_signature()` esigi.
#' @param coverage_threshold Kapsama uyari esigi (0-1). Varsayilan 0.70.
#' @param apply_log2 `scale = "linear"` iken log2 uygula.
#' @param sd_handling `encode_response()` politikasi.
#'
#' @return `data.frame` sonuc tablosu (`coverage_label`, `missing_genes` dahil).
#' @export
run_benchmark <- function(signatures,
                          datasets,
                          methods = c("mean", "zscore", "median", "rank"),
                          direction = "<",
                          min_genes = 3L,
                          coverage_threshold = 0.70,
                          apply_log2 = TRUE,
                          sd_handling = "exclude") {
  if (is.null(names(signatures)) || any(!nzchar(names(signatures)))) {
    stop("`signatures` isimli bir liste olmalidir.", call. = FALSE)
  }
  if (is.null(names(datasets)) || any(!nzchar(names(datasets)))) {
    stop("`datasets` isimli bir liste olmalidir.", call. = FALSE)
  }
  method_choices <- c("mean", "zscore", "median", "rank")
  methods <- match.arg(methods, choices = method_choices, several.ok = TRUE)

  # Farkli normalizasyon uyarisi (birlestirme yok; yalnizca bilgilendirme).
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
        scored <- score_signature(
          expr,
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
          missing_genes = paste(scored$coverage$missing, collapse = ","),
          low_coverage = low_cov,
          auc = NA_real_,
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

        if (isTRUE(scored$skipped)) {
          row_i <- row_i + 1L
          rows[[row_i]] <- base
          next
        }

        # Skorlar ile yaniti ornek id uzerinden hizala.
        sample_ids <- names(scored$scores)
        meta_idx <- match(sample_ids, ds$sample_meta$sample_id)
        resp <- enc$response[meta_idx]
        keep <- enc$kept[meta_idx] & is.finite(scored$scores)

        scores_use <- scored$scores[keep]
        resp_use <- resp[keep]

        # derived_from (Faz 1: alan tasinir; ozetten cikarma Faz 2)
        derived <- sig$derived_from
        if (is.null(derived)) {
          derived <- character(0)
        }
        base$independent <- !ds_name %in% derived

        if (length(unique(resp_use)) < 2L) {
          base$skipped <- TRUE
          base$skip_reason <- "Yanitta tek sinif kaldi (ROC hesaplanamaz)."
          row_i <- row_i + 1L
          rows[[row_i]] <- base
          next
        }

        roc_res <- tryCatch(
          roc_auc(scores_use, resp_use, direction = sig_dir),
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
        base$skipped <- FALSE
        base$skip_reason <- NA_character_

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
  out
}
