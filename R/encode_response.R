#' Klinik yanit etiketlerini ikili 0/1'e cevir
#'
#' Varsayilan: CR/PR = 1, PD = 0. Stabil hastalik (SD) varsayilan olarak
#' dislanir (`sd_handling = "exclude"`).
#'
#' @param labels Karakter veya faktor yanit etiketleri.
#' @param positive Yanitli sayilacak etiketler. Varsayilan `c("CR", "PR")`.
#' @param negative Yanitsiz sayilacak etiketler. Varsayilan `"PD"`.
#' @param sd_label Stabil hastalik etiketi. Varsayilan `"SD"`.
#' @param sd_handling `"exclude"`, `"as_0"` veya `"as_1"`.
#'
#' @return Liste: `response` (0/1, SD dislandiysa NA), `kept` (mantiksal),
#'   `n_excluded_sd`.
#' @export
encode_response <- function(labels,
                            positive = c("CR", "PR"),
                            negative = c("PD"),
                            sd_label = "SD",
                            sd_handling = c("exclude", "as_0", "as_1")) {
  sd_handling <- match.arg(sd_handling)
  labels_up <- toupper(as.character(labels))
  positive <- toupper(positive)
  negative <- toupper(negative)
  sd_label <- toupper(sd_label)

  response <- rep(NA_integer_, length(labels_up))
  response[labels_up %in% positive] <- 1L
  response[labels_up %in% negative] <- 0L

  is_sd <- labels_up %in% sd_label
  n_excluded_sd <- 0L

  if (any(is_sd)) {
    if (sd_handling == "exclude") {
      response[is_sd] <- NA_integer_
      n_excluded_sd <- sum(is_sd)
    } else if (sd_handling == "as_0") {
      response[is_sd] <- 0L
    } else {
      response[is_sd] <- 1L
    }
  }

  known <- unique(c(positive, negative, sd_label))
  unknown_mask <- !is.na(labels_up) & nzchar(labels_up) & !(labels_up %in% known)
  if (any(unknown_mask)) {
    unk <- unique(labels_up[unknown_mask])
    warning(
      "Unknown response label(s) encoded as NA: ",
      paste(unk, collapse = ", "),
      call. = FALSE
    )
  }

  kept <- !is.na(response)
  list(
    response = response,
    kept = kept,
    n_excluded_sd = n_excluded_sd
  )
}
