#' ROC-AUC hesapla (sabit yon, cift ters cevirme yok)
#'
#' Akis: `orient_score()` ile skorlari "yuksek = yanitli" yap, ardindan
#' `pROC::roc()` cagrisini her zaman `direction = "<"` ile yap.
#' AUC < 0.5 ise ters cevirmeden durustce raporlanir.
#'
#' @param scores Ornek skorlari (isimli vektor olabilir).
#' @param response 0/1 yanit vektoru (ayni uzunluk / hizali isimler).
#' @param direction Imza yonu: `"<"` veya `">"`. Yalnizca `orient_score`'a gider.
#'
#' @return Liste: `auc`, `n`, `n_pos`, `n_neg`, `roc` (pROC nesnesi).
#' @export
roc_auc <- function(scores, response, direction = c("<", ">")) {
  direction <- match.arg(direction)

  scores <- as.numeric(scores)
  response <- as.integer(response)

  if (length(scores) != length(response)) {
    stop("`scores` ve `response` ayni uzunlukta olmalidir.", call. = FALSE)
  }

  ok <- is.finite(scores) & !is.na(response)
  scores <- scores[ok]
  response <- response[ok]

  n_pos <- sum(response == 1L)
  n_neg <- sum(response == 0L)
  n <- length(response)

  if (n_pos == 0L || n_neg == 0L) {
    stop(
      "ROC icin her iki sinif da gerekli (n_pos=", n_pos, ", n_neg=", n_neg, ").",
      call. = FALSE
    )
  }

  # Tek yon hizalama; pROC her zaman direction = "<"
  scores_oriented <- orient_score(scores, direction = direction)

  roc_obj <- pROC::roc(
    response = response,
    predictor = scores_oriented,
    levels = c(0, 1),
    direction = "<",
    quiet = TRUE
  )

  list(
    auc = as.numeric(pROC::auc(roc_obj)),
    n = n,
    n_pos = n_pos,
    n_neg = n_neg,
    roc = roc_obj
  )
}
