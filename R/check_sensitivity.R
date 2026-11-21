#' Skorlama dayanikliligi (sensitivity) kontrolu
#'
#' Ayni veri seti x imza cifti icin farkli skorlama yontemlerinin AUC
#' araligini hesaplar. Aralik esigi asarsa imza "kirilgan" (fragile)
#' isaretlenir.
#'
#' Esik varsayilani `0.10` keyfidir; bilimsel bir standart degildir.
#' Sonuclara bakildiktan sonra degistirilmemelidir.
#'
#' @param results `run_benchmark()` cikti tablosu.
#' @param threshold AUC araligi esigi. Varsayilan `0.10` (keyfi).
#'
#' @return `data.frame`: dataset, signature, min_auc, max_auc, auc_range,
#'   n_methods, fragile, threshold.
#' @export
check_sensitivity <- function(results, threshold = 0.10) {
  if (!is.data.frame(results) || nrow(results) == 0L) {
    stop("`results` bos olamaz.", call. = FALSE)
  }
  req <- c("dataset", "signature", "method", "auc", "skipped")
  miss <- setdiff(req, names(results))
  if (length(miss) > 0L) {
    stop("`results` eksik sutunlar: ", paste(miss, collapse = ", "), call. = FALSE)
  }
  if (!is.numeric(threshold) || length(threshold) != 1L || is.na(threshold) ||
      threshold < 0) {
    stop("`threshold` negatif olmayan tek bir sayi olmalidir.", call. = FALSE)
  }

  df <- results[!isTRUE(results$skipped) & is.finite(results$auc), , drop = FALSE]
  if (nrow(df) == 0L) {
    return(data.frame(
      dataset = character(),
      signature = character(),
      min_auc = numeric(),
      max_auc = numeric(),
      auc_range = numeric(),
      n_methods = integer(),
      fragile = logical(),
      threshold = numeric(),
      stringsAsFactors = FALSE
    ))
  }

  split_key <- paste(df$dataset, df$signature, sep = "\r")
  groups <- split(df, split_key)

  rows <- lapply(groups, function(g) {
    aucs <- g$auc
    n_methods <- length(unique(g$method))
    min_auc <- min(aucs)
    max_auc <- max(aucs)
    auc_range <- max_auc - min_auc
    fragile <- if (n_methods < 2L) {
      NA
    } else {
      auc_range > threshold
    }
    data.frame(
      dataset = g$dataset[1],
      signature = g$signature[1],
      min_auc = min_auc,
      max_auc = max_auc,
      auc_range = auc_range,
      n_methods = as.integer(n_methods),
      fragile = fragile,
      threshold = threshold,
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out[order(out$dataset, out$signature), , drop = FALSE]
}
