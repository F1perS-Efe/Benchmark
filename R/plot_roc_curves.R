#' Secilen veri seti icin imza ROC egirlerini ciz
#'
#' @param datasets Isimli veri seti listesi (`validate_dataset` yapisi).
#' @param signatures Isimli imza listesi.
#' @param dataset_name Cizilecek veri seti adi.
#' @param method Skorlama yontemi. Varsayilan `"mean"`.
#' @param direction Varsayilan imza yonu.
#' @param apply_log2 `prepare_expression` icin.
#' @param sd_handling `encode_response` politikasi.
#' @param min_genes Minimum gen.
#'
#' @return `ggplot` nesnesi.
#' @export
plot_roc_curves <- function(datasets,
                            signatures,
                            dataset_name,
                            method = "mean",
                            direction = "<",
                            apply_log2 = TRUE,
                            sd_handling = "exclude",
                            min_genes = 3L) {
  if (!dataset_name %in% names(datasets)) {
    stop("Veri seti bulunamadi: ", dataset_name, call. = FALSE)
  }
  ds <- validate_dataset(datasets[[dataset_name]], name = dataset_name)
  expr <- collapse_duplicate_genes(as.matrix(ds$expr))
  expr <- prepare_expression(expr, scale = ds$info$scale, apply_log2 = apply_log2)
  enc <- encode_response(ds$sample_meta$response_label, sd_handling = sd_handling)

  # SD / NA yanit orneklerini skorlamadan once cikar (zscore/rank sizintisi yok)
  expr_kept <- expr[, enc$kept, drop = FALSE]
  y_kept <- enc$response[enc$kept]
  names(y_kept) <- colnames(expr_kept)

  curve_rows <- list()
  for (sig_name in names(signatures)) {
    sig <- signatures[[sig_name]]
    if (is.character(sig) || is.numeric(sig)) sig <- list(genes = sig)
    sig_dir <- if (!is.null(sig$direction)) sig$direction else direction

    scored <- score_signature(
      expr_kept,
      genes = sig$genes,
      genes_up = sig$genes_up,
      genes_down = sig$genes_down,
      method = method,
      min_genes = min_genes
    )
    if (isTRUE(scored$skipped)) next

    sample_ids <- names(scored$scores)
    resp <- y_kept[sample_ids]
    keep <- is.finite(scored$scores) & !is.na(resp)
    if (length(unique(resp[keep])) < 2L) next

    roc_res <- roc_auc(scored$scores[keep], resp[keep], direction = sig_dir)
    roc_obj <- roc_res$roc

    # pROC: specificities / sensitivities
    curve_rows[[length(curve_rows) + 1L]] <- data.frame(
      signature = sig_name,
      fpr = 1 - as.numeric(roc_obj$specificities),
      tpr = as.numeric(roc_obj$sensitivities),
      auc = roc_res$auc,
      stringsAsFactors = FALSE
    )
  }

  if (length(curve_rows) == 0L) {
    stop("Cizilecek ROC egrisi yok.", call. = FALSE)
  }
  df <- do.call(rbind, curve_rows)
  lab <- unique(df[, c("signature", "auc")])
  lab$label <- sprintf("%s (AUC=%.3f)", lab$signature, lab$auc)
  df$label <- lab$label[match(df$signature, lab$signature)]

  ggplot2::ggplot(
    df,
    ggplot2::aes(x = .data$fpr, y = .data$tpr, color = .data$label)
  ) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey60") +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title = paste0("ROC egrileri: ", dataset_name, " (", method, ")"),
      x = "1 - Ozgunluk (FPR)",
      y = "Duyarlilik (TPR)",
      color = "Imza"
    ) +
    ggplot2::theme_minimal(base_size = 12)
}
