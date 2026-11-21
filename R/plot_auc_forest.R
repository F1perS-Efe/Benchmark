#' AUC orman (forest) grafigi
#'
#' Her veri seti icin AUC noktasi ve DeLong %95 GA cubugu.
#' `summarize_auc()` ozeti verilirse IVW (yoksa n-agirlikli) referans cizgisi eklenir.
#'
#' @param results `run_benchmark()` cikti.
#' @param signature Tek imza adi; `NULL` ise ilk imza.
#' @param method Skorlama yontemi; `NULL` ise ilk yontem.
#' @param summary Opsiyonel `summarize_auc()` satiri / tablosu (ayni imza+method).
#'
#' @return `ggplot` nesnesi.
#' @export
plot_auc_forest <- function(results,
                            signature = NULL,
                            method = NULL,
                            summary = NULL) {
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
  df$lab <- ifelse(
    indep_ok,
    as.character(df$dataset),
    paste0(df$dataset, " (bagimsiz degil)")
  )
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
      x = "ROC-AUC (DeLong %95 GA)",
      y = NULL,
      caption = "Kesikli cizgi: sans duzeyi (0.5). var_zero satirlarda GA dar/eksik olabilir."
    ) +
    ggplot2::theme_bw()

  if (!is.null(summary) && is.data.frame(summary) && nrow(summary) > 0L) {
    sm <- summary
    sm <- sm[sm$signature == signature & sm$method == method, , drop = FALSE]
    if (nrow(sm) == 1L) {
      ref_val <- sm$auc_ivw[1]
      ref_lab <- "auc_ivw (ters-varyans)"
      if (!is.finite(ref_val)) {
        ref_val <- sm$auc_n_weighted[1]
        ref_lab <- "auc_n_weighted (n-agirlikli)"
      }
      if (is.finite(ref_val)) {
        p <- p + ggplot2::geom_vline(
          xintercept = ref_val,
          color = "steelblue",
          linetype = "solid"
        ) +
          ggplot2::labs(
            caption = paste0(
              "Mavi cizgi: ", ref_lab, " = ",
              sprintf("%.3f", ref_val),
              ". Kesikli: 0.5."
            )
          )
      }
    }
  }

  p
}
