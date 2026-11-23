#' Ifade matrisini skorlama oncesi hazirla
#'
#' `scale = "linear"` ise varsayilan olarak `log2(x + 1)` uygular.
#' Olcegi sessizce tahmin etmez. Negatif deger varken log donusumu hata verir.
#'
#' @param expr Gen x ornek sayisal matris.
#' @param scale `"log2"` veya `"linear"`.
#' @param apply_log2 `scale = "linear"` iken log2 uygulanacak mi? Varsayilan TRUE.
#' @param large_max_warn Buyuk maksimum deger uyarisi esigi. Varsayilan 100.
#'
#' @return Ayni boyutta hazirlanmis matris.
#' @export
prepare_expression <- function(expr,
                               scale = c("log2", "linear"),
                               apply_log2 = TRUE,
                               large_max_warn = 100) {
  scale <- match.arg(scale)
  if (!is.matrix(expr)) {
    expr <- as.matrix(expr)
  }
  if (!is.numeric(expr)) {
    stop("`expr` must be numeric.", call. = FALSE)
  }
  if (any(is.infinite(expr))) {
    stop("`expr` cannot contain infinite (Inf / -Inf) values.", call. = FALSE)
  }

  # Cok buyuk degerler log olcegi unutulmus olabilir; yalnizca uyari.
  max_val <- suppressWarnings(max(expr, na.rm = TRUE))
  if (is.finite(max_val) && max_val > large_max_warn && scale == "log2") {
    warning(
      "Expression matrix has very large values (max > ", large_max_warn,
      "). Verify that the scale is truly log2.",
      call. = FALSE
    )
  }

  if (scale == "linear" && isTRUE(apply_log2)) {
    min_val <- suppressWarnings(min(expr, na.rm = TRUE))
    if (is.finite(min_val) && min_val < 0) {
      stop(
        "Cannot apply log2(x + 1) to data containing negative values.",
        call. = FALSE
      )
    }
    expr <- log2(expr + 1)
  }

  expr
}
