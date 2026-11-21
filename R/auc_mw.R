#' Mann-Whitney U (rank-sum) ile ROC-AUC
#'
#' Vektorlestirilmis formul (esit degerlerde average ranks):
#' `AUC = (sum(rank_pos) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)`.
#'
#' Null permutasyon hizlandirmasi icin; `pROC` ile `<= 1e-12` eslesmelidir.
#'
#' @param scores Sayisal skorlar (yuksek = yanitli yonunde olmali).
#' @param response 0/1 yanit vektoru.
#'
#' @return Skaler AUC.
#' @export
auc_mw <- function(scores, response) {
  scores <- as.numeric(scores)
  response <- as.integer(response)
  ok <- is.finite(scores) & !is.na(response)
  scores <- scores[ok]
  response <- response[ok]

  n_pos <- sum(response == 1L)
  n_neg <- sum(response == 0L)
  if (n_pos == 0L || n_neg == 0L) {
    return(NA_real_)
  }

  ranks <- rank(scores, ties.method = "average")
  r_pos <- sum(ranks[response == 1L])
  (r_pos - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}
