#' Imza skorunu ornek basina hesapla
#'
#' Skorlama yanit etiketini (`response`) asla kullanmaz; tum ornekler
#' uzerinde etiketten bagimsiz hesaplanir.
#'
#' - `mean`: log olcekli ifadenin ornek basina ortalamasi
#' - `zscore`: her gen veri seti (kohort) genelinde z-score'a cevrilir, sonra
#'   ortalanir. Basit ortalama gen z-skorudur; `sqrt(k)` olcekli GSVA degildir.
#' - `median`: ham (log olcekli) ifadenin ornek basina medyani
#' - `rank`: her gen ornekler arasi (kohort-ici) average-rank'e cevrilir, 0-1
#'   araligina olceklenir, sonra ornek basina ortalanir. Bu `singscore` degildir;
#'   kohortlar arasi ham skor birlestirilmez.
#'
#' @param expr Gen x ornek matris (hazirlanmis / log olcekli).
#' @param genes Imza genleri (karakter).
#' @param method `"mean"`, `"zscore"`, `"median"` veya `"rank"`.
#' @param min_genes Kapsama sonrasi minimum gen sayisi. Varsayilan 3.
#'
#' @return Liste: `scores`, `coverage`, `n_all_na_dropped`, `skipped`, `skip_reason`.
#' @export
score_signature <- function(expr,
                            genes,
                            method = c("mean", "zscore", "median", "rank"),
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

  if (method == "rank") {
    # Her gen (satir): ornekler arasi average-rank, sonra [0, 1] olcek.
    sub <- t(apply(sub, 1L, function(x) {
      m <- sum(!is.na(x))
      if (m <= 1L) {
        # Tek veya sifir gecerli deger: bolme guvensiz; NA birak.
        return(rep(NA_real_, length(x)))
      }
      r <- rank(x, ties.method = "average", na.last = "keep")
      (r - 1) / (m - 1)
    }))
    if (!is.null(colnames(expr))) {
      colnames(sub) <- colnames(expr)[seq_len(ncol(sub))]
    }
    if (!is.null(rownames(expr)) || length(coverage$found) > 0L) {
      rownames(sub) <- coverage$found
    }
  }

  # Ornek basina ozet; gen NA'lari ornek icinde yok sayilir.
  scores <- apply(sub, 2L, function(x) {
    if (all(is.na(x))) {
      return(NA_real_)
    }
    if (method == "median") {
      stats::median(x, na.rm = TRUE)
    } else {
      mean(x, na.rm = TRUE)
    }
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
