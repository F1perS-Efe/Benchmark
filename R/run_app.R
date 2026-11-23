#' SigBenchR Shiny arayuzunu baslat
#'
#' Uygulama `inst/shiny/` altindadir. Shiny paket bagimliligi Suggests'tedir.
#' Yerel dosya yukleme varsayilan olarak kapali (fail-closed) kalir; acmak icin
#' `allow_upload = TRUE` verin. Uzak sunucu / shinyapps.io giris noktasi
#' (`app.R`) `allow_upload = FALSE` ile cagirir.
#'
#' @param allow_upload Yerel CSV/TSV kohort yuklemeyi ac. Varsayilan `FALSE`.
#' @param ... `shiny::runApp()`'e iletilen ek argumanlar.
#'
#' @return `shiny::runApp()` donus degeri (gorunmez).
#' @export
run_app <- function(allow_upload = FALSE, ...) {
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop(
      "Install 'shiny' to launch the UI: install.packages(\"shiny\")",
      call. = FALSE
    )
  }
  app_dir <- system.file("shiny", package = "SigBenchR")
  if (!nzchar(app_dir) || !dir.exists(app_dir)) {
    stop("Shiny app directory not found (inst/shiny).", call. = FALSE)
  }
  old_allow <- shiny::getShinyOption("sigbenchr_allow_upload")
  old_opts <- options(shiny.maxRequestSize = 100 * 1024^2)
  on.exit(
    {
      options(old_opts)
      shiny::shinyOptions(sigbenchr_allow_upload = old_allow)
    },
    add = TRUE
  )
  shiny::shinyOptions(sigbenchr_allow_upload = isTRUE(allow_upload))
  shiny::runApp(app_dir, ...)
}
