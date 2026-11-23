#' Veri setleri boyunca imza AUC korelasyonu
#'
#' Secilen skorlama yonteminde her imza icin veri setleri uzerindeki AUC
#' vektorleri arasinda Pearson korelasyonu hesaplar. Iki noktada Pearson
#' daima +/-1 oldugundan en az 3 ortak sonlu veri seti gerekir
#' (`n_ok < 3` ise `correlation = NA_real_`).
#'
#' @param results `run_benchmark()` cikti tablosu.
#' @param method Skorlama yontemi. `NULL` ise ilk yontem.
#'
#' @return Uzun `data.frame`: `signature1`, `signature2`, `correlation`,
#'   `n_datasets`, `method`, `low_k_unreliable`.
#' @export
auc_correlation <- function(results, method = NULL) {
  if (!is.data.frame(results) || nrow(results) == 0L) {
    stop("`results` cannot be empty.", call. = FALSE)
  }
  req <- c("dataset", "signature", "method", "auc", "skipped")
  miss <- setdiff(req, names(results))
  if (length(miss) > 0L) {
    stop("`results` is missing required columns: ", paste(miss, collapse = ", "), call. = FALSE)
  }

  df <- results[!(results$skipped %in% TRUE) & is.finite(results$auc), , drop = FALSE]
  if (is.null(method)) {
    method <- df$method[1]
  }
  df <- df[df$method == method, , drop = FALSE]
  if (nrow(df) == 0L) {
    stop("No AUC rows available for the selected method.", call. = FALSE)
  }

  sigs <- unique(as.character(df$signature))
  if (length(sigs) < 2L) {
    stop("At least two signatures are required for correlation.", call. = FALSE)
  }

  # Genis matris: satir = veri seti, sutun = imza
  datasets <- unique(as.character(df$dataset))
  mat <- matrix(NA_real_, nrow = length(datasets), ncol = length(sigs),
                dimnames = list(datasets, sigs))
  for (i in seq_len(nrow(df))) {
    mat[as.character(df$dataset[i]), as.character(df$signature[i])] <- df$auc[i]
  }

  rows <- list()
  k <- 0L
  any_low_k <- FALSE
  for (i in seq_len(length(sigs) - 1L)) {
    for (j in seq.int(i + 1L, length(sigs))) {
      x <- mat[, sigs[i]]
      y <- mat[, sigs[j]]
      ok <- is.finite(x) & is.finite(y)
      n_ok <- sum(ok)
      if (n_ok < 5L) {
        any_low_k <- TRUE
      }
      # 2 noktada Pearson daima +/-1; anlamli korelasyon icin >= 3
      cor_val <- if (n_ok >= 3L) {
        suppressWarnings(stats::cor(x[ok], y[ok], method = "pearson"))
      } else {
        NA_real_
      }
      k <- k + 1L
      rows[[k]] <- data.frame(
        signature1 = sigs[i],
        signature2 = sigs[j],
        correlation = cor_val,
        n_datasets = as.integer(n_ok),
        method = method,
        low_k_unreliable = as.logical(n_ok < 5L),
        stringsAsFactors = FALSE
      )
    }
  }
  if (any_low_k) {
    warning(
      "Cross-cohort AUC correlation computed with fewer than 5 cohorts (n_datasets < 5) has high sampling variance and is flagged with low_k_unreliable = TRUE.",
      call. = FALSE
    )
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
