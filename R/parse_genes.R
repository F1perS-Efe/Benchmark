#' Kullanici gen listesini ayristir
#'
#' Satir, virgul, noktali virgul veya boslukla ayrilmis girdiyi
#' karakter vektorune cevirir.
#'
#' @param text Tek bir karakter dizisi.
#'
#' @return Temizlenmis gen sembolleri.
#' @export
parse_genes <- function(text) {
  if (is.null(text) || !nzchar(trimws(text))) {
    return(character(0))
  }
  parts <- unlist(strsplit(text, "[,;\\s]+", perl = TRUE), use.names = FALSE)
  parts <- trimws(parts)
  parts <- parts[nzchar(parts)]
  unique(parts)
}
