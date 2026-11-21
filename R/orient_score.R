#' Skor yonunu "yuksek = yanitli" olacak sekilde hizala
#'
#' Asagi-regule imzalar icin `direction = ">"` verildiginde skorlar
#' isaret degistirilir (`-scores`). Boylece sonraki ROC cagrilari
#' her zaman `direction = "<"` kullanabilir (cift ters cevirme olmaz).
#'
#' @param scores Ornek basina sayisal skor vektoru.
#' @param direction `"<"` (yuksek skor = yanitli, varsayilan) veya
#'   `">"` (dusuk skor = yanitli; skorlar ters cevrilir).
#'
#' @return Yonu hizalanmis skor vektoru.
#' @export
orient_score <- function(scores, direction = c("<", ">")) {
  direction <- match.arg(direction)
  scores <- as.numeric(scores)
  if (direction == ">") {
    return(-scores)
  }
  scores
}
