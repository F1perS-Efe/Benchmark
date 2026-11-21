#' Imza skorunu ornek basina hesapla
#'
#' Skorlama yanit etiketini (`response`) asla kullanmaz; tum ornekler
#' uzerinde etiketten bagimsiz hesaplanir.
#'
#' - `mean`: log olcekli ifadenin ornek basina ortalamasi
#' - `zscore`: her gen veri seti genelinde z-score'a cevrilir, sonra ortalanir
#'
#' @param expr Gen x ornek matris (hazirlanmis / log olcekli).
#' @param genes Imza genleri (karakter).
#' @param method `"mean"` veya `"zscore"`.
#' @param min_genes Kapsama sonrasi minimum gen sayisi. Varsayilan 3.
#'
#' @return Liste: `scores`, `coverage`, `n_all_na_dropped`, `skipped`, `skip_reason`.
#' @export
score_signature <- function(expr,
                            genes,
                            method = c("mean", "zscore"),
                            min_genes = 3L) {
  method <- match.arg(method)
  coverage <- match_genes(expr, genes)

  if (coverage$n_found < min_genes) {
    return(list(
      scores = NULL,
      coverage = coverage,
      n_all_na_dropped = 0L,
      skipped = TRUE,
      skip_reason = sprintf(
        "Bulunan gen sayisi (%d) min_genes (%d) altinda.",
        coverage$n_found, min_genes
      )
    ))
  }

  sub <- expr[coverage$found, , drop = FALSE]

  if (method == "zscore") {
    # Her gen: tum ornekler uzerinde z-score (yanit kullanilmaz).
    sub <- t(scale(t(sub)))
    # Sabit genler (sd = 0) NA uretir; skorlamada na.rm ile idare edilir.
  }

  # Ornek basina ortalama; gen NA'lari ornek icinde yok sayilir.
  scores <- apply(sub, 2L, function(x) {
    if (all(is.na(x))) {
      return(NA_real_)
    }
    mean(x, na.rm = TRUE)
  })

  all_na <- is.na(scores)
  n_all_na_dropped <- sum(all_na)
  if (n_all_na_dropped > 0L) {
    scores <- scores[!all_na]
  }

  list(
    scores = scores,
    coverage = coverage,
    n_all_na_dropped = n_all_na_dropped,
    skipped = FALSE,
    skip_reason = NA_character_
  )
}
