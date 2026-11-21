#' SigBenchR Shiny arayuzunu baslat
#'
#' Uygulama `inst/shiny/` altindadir. Shiny paket bagimliligi Suggests'tedir.
#'
#' @param ... `shiny::runApp()`'e iletilen ek argumanlar.
#'
#' @return `shiny::runApp()` donus degeri (gorunmez).
#' @export
run_app <- function(...) {
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop(
      "Shiny arayuzu icin 'shiny' paketini kurun: install.packages(\"shiny\")",
      call. = FALSE
    )
  }
  app_dir <- system.file("shiny", package = "SigBenchR")
  if (!nzchar(app_dir) || !dir.exists(app_dir)) {
    stop("Shiny uygulama dizini bulunamadi (inst/shiny).", call. = FALSE)
  }
  shiny::runApp(app_dir, ...)
}
