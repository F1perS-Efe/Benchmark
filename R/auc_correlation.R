#' Veri setleri boyunca imza AUC korelasyonu
#'
#' Secilen skorlama yonteminde her imza icin veri setleri uzerindeki AUC
#' vektorleri arasinda Pearson korelasyonu hesaplar. En az 2 ortak veri seti
#' gerekir; aksi halde ilgili cift icin `NA`.
#'
#' @param results `run_benchmark()` cikti tablosu.
#' @param method Skorlama yontemi. `NULL` ise ilk yontem.
#'
#' @return Uzun `data.frame`: `signature1`, `signature2`, `correlation`,
#'   `n_datasets`, `method`.
#' @export
auc_correlation <- function(results, method = NULL) {
  if (!is.data.frame(results) || nrow(results) == 0L) {
    stop("`results` bos olamaz.", call. = FALSE)
  }
  req <- c("dataset", "signature", "method", "auc", "skipped")
  miss <- setdiff(req, names(results))
  if (length(miss) > 0L) {
    stop("`results` eksik sutunlar: ", paste(miss, collapse = ", "), call. = FALSE)
  }

  df <- results[!isTRUE(results$skipped) & is.finite(results$auc), , drop = FALSE]
  if (is.null(method)) {
    method <- df$method[1]
  }
  df <- df[df$method == method, , drop = FALSE]
  if (nrow(df) == 0L) {
    stop("Secilen yontem icin AUC satiri yok.", call. = FALSE)
  }

  sigs <- unique(as.character(df$signature))
  if (length(sigs) < 2L) {
    stop("Korelasyon icin en az iki imza gerekir.", call. = FALSE)
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
  for (i in seq_len(length(sigs) - 1L)) {
    for (j in seq.int(i + 1L, length(sigs))) {
      x <- mat[, sigs[i]]
      y <- mat[, sigs[j]]
      ok <- is.finite(x) & is.finite(y)
      n_ok <- sum(ok)
      cor_val <- if (n_ok >= 2L) {
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
        stringsAsFactors = FALSE
      )
    }
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
