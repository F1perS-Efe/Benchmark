#' Null p-degerlerine BH duzeltmesi (imza x veri seti ailesi)
#'
#' Coklu karsilastirma ailesi = ayni `method` icindeki tum
#' `signature` x `dataset` satirlarinin `null_p` degerleri.
#' Ham `null_p` korunur; `null_p_adj_bh` eklenir.
#'
#' @param results `run_benchmark()` cikti (veya `null_p` iceren tablo).
#' @param p_col Ham p sutunu. Varsayilan `"null_p"`.
#'
#' @return Ayni tablo + `null_p_adj_bh` sutunu.
#' @export
adjust_null_p <- function(results, p_col = "null_p") {
  if (!is.data.frame(results) || nrow(results) == 0L) {
    stop("`results` bos olamaz.", call. = FALSE)
  }
  if (!p_col %in% names(results)) {
    stop("Sutun bulunamadi: ", p_col, call. = FALSE)
  }
  if (!"method" %in% names(results)) {
    stop("`results` icinde `method` sutunu gerekli.", call. = FALSE)
  }

  out <- results
  out$null_p_adj_bh <- NA_real_
  methods <- unique(as.character(out$method))
  for (m in methods) {
    idx <- which(out$method == m & is.finite(out[[p_col]]))
    if (length(idx) > 0L) {
      out$null_p_adj_bh[idx] <- stats::p.adjust(out[[p_col]][idx], method = "BH")
    }
  }
  out
}
