#' ROC-AUC hesapla (sabit yon, cift ters cevirme yok)
#'
#' Akis: `orient_score()` ile skorlari "yuksek = yanitli" yap, ardindan
#' `pROC::roc()` cagrisini her zaman `direction = "<"` ile yap.
#' AUC < 0.5 ise ters cevirmeden durustce raporlanir.
#'
#' DeLong yontemiyle %95 GA (`ci.auc`) ve AUC varyansi (`var`) eklenir.
#'
#' @param scores Ornek skorlari (isimli vektor olabilir).
#' @param response 0/1 yanit vektoru (ayni uzunluk / hizali isimler).
#' @param direction Imza yonu: `"<"` veya `">"`. Yalnizca `orient_score`'a gider.
#' @param min_n Kucuk kohort icin toplam n esigi. Varsayilan 20.
#' @param min_per_class Kucuk kohort icin sinif basi minimum. Varsayilan 5.
#'
#' @return Liste: `auc`, `auc_ci_low`, `auc_ci_high`, `auc_var`, `var_zero`,
#'   `small_cohort`, `n`, `n_pos`, `n_neg`, `roc`, `scores_oriented`.
#' @export
roc_auc <- function(scores,
                    response,
                    direction = c("<", ">"),
                    min_n = 20L,
                    min_per_class = 5L) {
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

  auc_val <- as.numeric(pROC::auc(roc_obj))
  # AUC==1 iken pROC uyari verir; var_zero ile zaten isaretliyoruz
  ci <- suppressWarnings(as.numeric(pROC::ci.auc(roc_obj, method = "delong")))
  # ci.auc donusu: [low, mid/auc, high]
  auc_ci_low <- ci[1]
  auc_ci_high <- ci[3]

  auc_var <- suppressWarnings(as.numeric(pROC::var(roc_obj, method = "delong")))
  var_zero <- is.na(auc_var) || isTRUE(auc_var == 0)

  small_cohort <- (n < min_n) || (min(n_pos, n_neg) < min_per_class)

  list(
    auc = auc_val,
    auc_ci_low = auc_ci_low,
    auc_ci_high = auc_ci_high,
    auc_var = auc_var,
    var_zero = var_zero,
    small_cohort = small_cohort,
    n = n,
    n_pos = n_pos,
    n_neg = n_neg,
    roc = roc_obj,
    scores_oriented = scores_oriented,
    response = response
  )
}
