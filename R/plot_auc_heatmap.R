#' Imza x veri seti AUC isi haritasi
#'
#' Kapsama esigi altindaki (`low_coverage`) hucreler gri cizilir; bu durum
#' lejant/caption'da acikca belirtilir (sessizce gecilmez).
#'
#' @param results `run_benchmark()` ciktisi.
#' @param method Belirli bir skorlama yontemi suz. `NULL` ise ilk yontem.
#' @param coverage_threshold 0-1. Varsayilan 0.70.
#' @param show_normalization If `TRUE` (default) and a `normalization` column
#'   is present, cohort axis labels become `"cohort [TPM]"` style and a
#'   heterogeneous-scale caption note is appended.
#'
#' @return `ggplot` nesnesi.
#' @export
plot_auc_heatmap <- function(results,
                             method = NULL,
                             coverage_threshold = 0.70,
                             show_normalization = TRUE) {
  if (!is.data.frame(results) || nrow(results) == 0L) {
    stop("`results` bos olamaz.", call. = FALSE)
  }

  df <- results
  if (is.null(method)) {
    method <- df$method[1]
  }
  df <- df[df$method == method & !(df$skipped %in% TRUE), , drop = FALSE]
  if (nrow(df) == 0L) {
    stop("Secilen yontem icin cizilecek satir yok.", call. = FALSE)
  }

  # Dusuk kapsama: esigi dogrudan yeniden hesapla (onceki TRUE kilitlenmesin)
  if ("coverage_pct" %in% names(df)) {
    df$low_coverage <- df$coverage_pct < (100 * coverage_threshold)
  } else if (!"low_coverage" %in% names(df)) {
    df$low_coverage <- FALSE
  }

  df$auc_plot <- df$auc
  df$auc_plot[df$low_coverage] <- NA_real_

  cov_txt <- if ("coverage_label" %in% names(df)) {
    as.character(df$coverage_label)
  } else {
    sprintf("%.0f%%", df$coverage_pct)
  }
  df$cell_label <- ifelse(
    df$low_coverage,
    paste0("dusuk\n", cov_txt),
    sprintf("%.2f", df$auc)
  )

  use_norm <- isTRUE(show_normalization) && "normalization" %in% names(df)
  if (use_norm) {
    norm_lab <- as.character(df$normalization)
    norm_lab[is.na(norm_lab) | !nzchar(norm_lab)] <- "unknown"
    df$dataset_label <- paste0(as.character(df$dataset), " [", norm_lab, "]")
  } else {
    df$dataset_label <- as.character(df$dataset)
  }

  caption <- paste0(
    "Gri hucreler = dusuk gen kapsama (kapsama < %",
    round(100 * coverage_threshold),
    "). Gri, eksik kanit uyarisidir; sessizce yok sayilmamalidir."
  )
  if (use_norm) {
    caption <- paste0(
      caption,
      " Note: Cohorts use heterogeneous expression scales; ",
      "direct cross-cohort expression pooling is not performed."
    )
  }

  ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = .data$dataset_label,
      y = .data$signature,
      fill = .data$auc_plot
    )
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.4) +
    ggplot2::geom_text(
      ggplot2::aes(label = .data$cell_label),
      size = 3
    ) +
    ggplot2::scale_fill_gradient(
      low = "#f7fbff",
      high = "#08306b",
      limits = c(0, 1),
      na.value = "grey70",
      name = "AUC"
    ) +
    ggplot2::labs(
      title = paste0("AUC isi haritasi (", method, ")"),
      x = "Veri seti",
      y = "Imza",
      caption = caption
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 30, hjust = 1)
    )
}
