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
    stop("`expr` sayisal olmalidir.", call. = FALSE)
  }

  # Cok buyuk degerler log olcegi unutulmus olabilir; yalnizca uyari.
  finite_vals <- expr[is.finite(expr)]
  if (length(finite_vals) > 0L && max(finite_vals) > large_max_warn &&
      scale == "log2") {
    warning(
      "Ifade matrisinde cok buyuk degerler var (max > ", large_max_warn,
      "). Olcek gercekten log2 mi kontrol edin.",
      call. = FALSE
    )
  }

  if (scale == "linear" && isTRUE(apply_log2)) {
    if (any(finite_vals < 0, na.rm = TRUE)) {
      stop(
        "Negatif deger iceren veriye log2(x + 1) uygulanamaz.",
        call. = FALSE
      )
    }
    expr <- log2(expr + 1)
  }

  expr
}
