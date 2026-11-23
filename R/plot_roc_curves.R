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
    stop("Dataset not found: ", dataset_name, call. = FALSE)
  }
  prep <- .prepare_cohort(
    datasets[[dataset_name]],
    name = dataset_name,
    apply_log2 = apply_log2,
    sd_handling = sd_handling,
    filter_response = TRUE
  )
  expr_full <- prep$expr_kept
  y_kept <- prep$y_kept

  curve_rows <- list()
  for (sig_name in names(signatures)) {
    sig <- .as_signature(signatures[[sig_name]], name = sig_name)
    sig_dir <- .normalize_roc_direction(sig$direction, default = direction)
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
    if (isTRUE(scored$skipped)) next

    sample_ids <- intersect(names(scored$scores), names(y_kept))
    resp <- y_kept[sample_ids]
    sc_kept <- scored$scores[sample_ids]
    keep <- is.finite(sc_kept) & !is.na(resp)
    if (length(unique(resp[keep])) < 2L) next

    roc_res <- roc_auc(sc_kept[keep], resp[keep], direction = sig_dir)
    roc_obj <- roc_res$roc

    # pROC: specificities / sensitivities
    curve_rows[[length(curve_rows) + 1L]] <- data.frame(
      signature = sig_name,
      fpr = 1 - as.numeric(roc_obj$specificities),
      tpr = as.numeric(roc_obj$sensitivities),
      auc = roc_res$auc,
      auc_ci_low = roc_res$auc_ci_low,
      auc_ci_high = roc_res$auc_ci_high,
      n = roc_res$n,
      n_pos = roc_res$n_pos,
      n_neg = roc_res$n_neg,
      stringsAsFactors = FALSE
    )
  }

  if (length(curve_rows) == 0L) {
    stop("No ROC curves available to plot.", call. = FALSE)
  }
  df <- do.call(rbind, curve_rows)
  lab <- unique(df[, c("signature", "auc", "auc_ci_low", "auc_ci_high")])
  lab$label <- ifelse(
    is.finite(lab$auc_ci_low) & is.finite(lab$auc_ci_high),
    sprintf("%s (AUC=%.3f [%.3f-%.3f])", lab$signature, lab$auc, lab$auc_ci_low, lab$auc_ci_high),
    sprintf("%s (AUC=%.3f)", lab$signature, lab$auc)
  )
  df$label <- lab$label[match(df$signature, lab$signature)]
  n_tot <- df$n[1]
  n_pos_tot <- df$n_pos[1]
  n_neg_tot <- df$n_neg[1]
  sub_txt <- if (is.finite(n_tot) && is.finite(n_pos_tot) && is.finite(n_neg_tot)) {
    sprintf("n = %d (responders = %d, non-responders = %d) | Dashed line: chance level (AUC = 0.50)", n_tot, n_pos_tot, n_neg_tot)
  } else {
    "Dashed line: chance level (AUC = 0.50)"
  }

  ggplot2::ggplot(
    df,
    ggplot2::aes(x = .data$fpr, y = .data$tpr, color = .data$label)
  ) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey60") +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title = paste0("ROC curves: ", dataset_name, " (", method, ")"),
      subtitle = sub_txt,
      x = "1 - Specificity (False Positive Rate)",
      y = "Sensitivity (True Positive Rate)",
      color = "Signature (95% DeLong CI)"
    ) +
    ggplot2::theme_minimal(base_size = 12)
}
