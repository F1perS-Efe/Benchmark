#' AUC orman (forest) grafigi
#'
#' Her veri seti icin AUC noktasi ve DeLong 95 pct GA cubugu.
#' `summarize_auc()` ozeti verilirse IVW (yoksa n-agirlikli) referans cizgisi
#' eklenir. Ozet `auc_re` ile sonlu HKSJ GA sutunlari iceriyorsa logit-AUC
#' REML nokta tahmini ve 95 pct HKSJ yatay hata cubugu da cizilir; tek kohortta
#' (`auc_re_ci_low` NA) yalniz mevcut ozet cizgisi korunur.
#'
#' @param results `run_benchmark()` cikti.
#' @param signature Tek imza adi; `NULL` ise ilk imza.
#' @param method Skorlama yontemi; `NULL` ise ilk yontem.
#' @param summary Opsiyonel `summarize_auc()` satiri / tablosu (ayni imza+method).
#' @param show_normalization If `TRUE` (default) and a `normalization` column
#'   is present, cohort labels become `"cohort [TPM]"` style and a
#'   heterogeneous-scale caption note is appended.
#'
#' @return `ggplot` nesnesi.
#' @export
plot_auc_forest <- function(results,
                            signature = NULL,
                            method = NULL,
                            summary = NULL,
                            show_normalization = TRUE) {
  if (!is.data.frame(results) || nrow(results) == 0L) {
    stop("`results` bos olamaz.", call. = FALSE)
  }
  df <- results[!(results$skipped %in% TRUE) & is.finite(results$auc), , drop = FALSE]
  if (is.null(method)) method <- df$method[1]
  if (is.null(signature)) signature <- df$signature[1]
  df <- df[df$method == method & df$signature == signature, , drop = FALSE]
  if (nrow(df) == 0L) {
    stop("Secilen imza/yontem icin satir yok.", call. = FALSE)
  }

  if (!"auc_ci_low" %in% names(df)) df$auc_ci_low <- NA_real_
  if (!"auc_ci_high" %in% names(df)) df$auc_ci_high <- NA_real_

  indep_ok <- if ("independent" %in% names(df)) {
    is.na(df$independent) | as.logical(df$independent)
  } else {
    rep(TRUE, nrow(df))
  }
  use_norm <- isTRUE(show_normalization) && "normalization" %in% names(df)
  cohort_base <- as.character(df$dataset)
  if (use_norm) {
    norm_lab <- as.character(df$normalization)
    norm_lab[is.na(norm_lab) | !nzchar(norm_lab)] <- "unknown"
    cohort_base <- paste0(cohort_base, " [", norm_lab, "]")
  }
  df$lab <- ifelse(
    indep_ok,
    cohort_base,
    paste0(cohort_base, " (bagimsiz degil)")
  )

  caption <- "Kesikli cizgi: sans duzeyi (0.5). var_zero satirlarda GA dar/eksik olabilir."
  if (use_norm) {
    caption <- paste0(
      caption,
      " Note: Cohorts use heterogeneous expression scales; ",
      "direct cross-cohort expression pooling is not performed."
    )
  }
  has_hksj <- FALSE
  sm <- NULL
  if (!is.null(summary) && is.data.frame(summary) && nrow(summary) > 0L) {
    sm <- summary
    sm <- sm[sm$signature == signature & sm$method == method, , drop = FALSE]
    if (nrow(sm) == 1L &&
        all(c("auc_re", "auc_re_ci_low", "auc_re_ci_high") %in% names(sm)) &&
        is.finite(sm$auc_re[1]) &&
        is.finite(sm$auc_re_ci_low[1]) &&
        is.finite(sm$auc_re_ci_high[1])) {
      has_hksj <- TRUE
      re_row <- data.frame(
        auc = sm$auc_re[1],
        auc_ci_low = sm$auc_re_ci_low[1],
        auc_ci_high = sm$auc_re_ci_high[1],
        lab = "REML + HKSJ (95 pct)",
        stringsAsFactors = FALSE
      )
      df <- rbind(
        df[, c("auc", "auc_ci_low", "auc_ci_high", "lab"), drop = FALSE],
        re_row
      )
      caption <- paste0(
        caption,
        " Alt satir: logit-AUC REML nokta + kesilmis HKSJ 95 pct GA."
      )
    }
  }

  df$lab <- factor(df$lab, levels = rev(unique(df$lab)))

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = .data$auc, y = .data$lab)
  ) +
    ggplot2::geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey50") +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = .data$auc_ci_low, xmax = .data$auc_ci_high),
      orientation = "y",
      width = 0.2,
      na.rm = TRUE
    ) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::labs(
      title = paste0("AUC forest: ", signature, " (", method, ")"),
      x = "ROC-AUC (DeLong 95 pct GA)",
      y = NULL,
      caption = caption
    ) +
    ggplot2::theme_bw()

  if (!is.null(sm) && nrow(sm) == 1L) {
    ref_val <- if ("auc_re" %in% names(sm) && is.finite(sm$auc_re[1])) {
      sm$auc_re[1]
    } else {
      sm$auc_ivw[1]
    }
    ref_lab <- if ("auc_re" %in% names(sm) && is.finite(sm$auc_re[1])) {
      "auc_re (logit-AUC REML)"
    } else {
      "auc_ivw (ters-varyans)"
    }
    if (!is.finite(ref_val)) {
      ref_val <- sm$auc_n_weighted[1]
      ref_lab <- "auc_n_weighted (n-agirlikli)"
    }
    if (is.finite(ref_val)) {
      cap_extra <- paste0(
        "Mavi cizgi: ", ref_lab, " = ",
        sprintf("%.3f", ref_val),
        ". Kesikli: 0.5."
      )
      if (isTRUE(has_hksj)) {
        cap_extra <- paste0(
          cap_extra,
          " Alt satir: logit-AUC REML + HKSJ 95 pct GA."
        )
      }
      if (use_norm) {
        cap_extra <- paste0(
          cap_extra,
          " Note: Cohorts use heterogeneous expression scales; ",
          "direct cross-cohort expression pooling is not performed."
        )
      }
      p <- p + ggplot2::geom_vline(
        xintercept = ref_val,
        color = "steelblue",
        linetype = "solid"
      ) +
        ggplot2::labs(caption = cap_extra)
    }
  }

  p
}
