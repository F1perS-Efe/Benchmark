#' Imza x veri seti AUC isi haritasi
#'
#' Kapsama esigi altindaki hucreler gri cizilir.
#'
#' @param results `run_benchmark()` ciktisi.
#' @param method Belirli bir skorlama yontemi suz. `NULL` ise ilk yontem.
#' @param coverage_threshold 0-1. Varsayilan 0.70.
#'
#' @return `ggplot` nesnesi.
#' @export
plot_auc_heatmap <- function(results,
                             method = NULL,
                             coverage_threshold = 0.70) {
  if (!is.data.frame(results) || nrow(results) == 0L) {
    stop("`results` bos olamaz.", call. = FALSE)
  }

  df <- results
  if (is.null(method)) {
    method <- df$method[1]
  }
  df <- df[df$method == method & !isTRUE(df$skipped), , drop = FALSE]
  if (nrow(df) == 0L) {
    stop("Secilen yontem icin cizilecek satir yok.", call. = FALSE)
  }

  # Dusuk kapsama: AUC'yi NA yap -> gri kutu
  df$auc_plot <- df$auc
  low <- df$low_coverage | df$coverage_pct < (100 * coverage_threshold)
  df$auc_plot[low] <- NA_real_

  ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = .data$dataset,
      y = .data$signature,
      fill = .data$auc_plot
    )
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.4) +
    ggplot2::geom_text(
      ggplot2::aes(
        label = ifelse(
          is.na(.data$auc_plot),
          sprintf("NA\n(%.0f%%)", .data$coverage_pct),
          sprintf("%.2f", .data$auc)
        )
      ),
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
      caption = paste0(
        "Gri hucreler: kapsama < %", round(100 * coverage_threshold),
        " veya atlanan sonuc"
      )
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 30, hjust = 1)
    )
}
