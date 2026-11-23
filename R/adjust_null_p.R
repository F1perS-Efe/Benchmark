#' Null p-degerlerine BH duzeltmesi (imza x veri seti ve yontemler arasi aile)
#'
#' Coklu karsilastirma ailesi:
#' - `"within_method"` (varsayilan): Ayni `method` icindeki tum `signature` x `dataset`
#'   satirlarinin `null_p` degerleri (`null_p_adj_bh`).
#' - `"global"`: Tablodaki tum `method` x `signature` x `dataset` kombinasyonlarinin
#'   ortak ailesi uzerinden BH duzeltmesi.
#'
#' Her iki modda da yontemler-arasi global duzeltilmis p-degeri `null_p_adj_bh_global`
#' sutununda her zaman raporlanir; `scope = "global"` secildiginde `null_p_adj_bh`
#' de global duzeltilmis degere esitlenir.
#'
#' @param results `run_benchmark()` cikti (veya `null_p` iceren tablo).
#' @param p_col Ham p sutunu. Varsayilan `"null_p"`.
#' @param scope Coklu test ailesi kapsami: `"within_method"` veya `"global"`.
#'
#' @return Ayni tablo + `null_p_adj_bh` ve `null_p_adj_bh_global` sutunlari.
#' @export
adjust_null_p <- function(results,
                          p_col = "null_p",
                          scope = c("within_method", "global")) {
  scope <- match.arg(scope)
  if (!is.data.frame(results) || nrow(results) == 0L) {
    stop("`results` cannot be empty.", call. = FALSE)
  }
  if (!p_col %in% names(results)) {
    stop("Column not found: ", p_col, call. = FALSE)
  }
  if (!"method" %in% names(results)) {
    stop("`results` must contain a `method` column.", call. = FALSE)
  }

  out <- results
  out$null_p_adj_bh <- NA_real_
  out$null_p_adj_bh_global <- NA_real_

  idx_all <- which(is.finite(out[[p_col]]))
  if (length(idx_all) > 0L) {
    out$null_p_adj_bh_global[idx_all] <- stats::p.adjust(
      out[[p_col]][idx_all],
      method = "BH"
    )
  }

  if (scope == "global") {
    out$null_p_adj_bh <- out$null_p_adj_bh_global
  } else {
    methods <- unique(as.character(out$method))
    for (m in methods) {
      idx <- which(out$method == m & is.finite(out[[p_col]]))
      if (length(idx) > 0L) {
        out$null_p_adj_bh[idx] <- stats::p.adjust(out[[p_col]][idx], method = "BH")
      }
    }
  }
  out
}
