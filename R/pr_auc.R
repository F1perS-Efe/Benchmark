#' PR-AUC ve tabakali bootstrap %95 GA
#'
#' Skor once `orient_score()` ile "yuksek = yanitli" yapilir, sonra
#' `PRROC::pr.curve(scores.class0 = ..., weights.class0 = response)`
#' ile Davis-Goadrich PR-AUC hesaplanir.
#'
#' Taban cizgisi (prevalans): `n_pos / n`.
#' GA: pozitif ve negatif indeksler ayri ayri iadeli orneklenir
#' (`n_pos` ve `n_neg` sabit); yuzdelik (0.025, 0.975).
#'
#' @param scores Ornek skorlari.
#' @param response 0/1 yanit.
#' @param direction Imza yonu (`orient_score` icin).
#' @param n_boot Bootstrap tekrari. Varsayilan 1000; `0` ise CI `NA`.
#' @param seed RNG tohumu (`withr::with_seed`).
#'
#' @return Liste: `pr_auc`, `pr_baseline`, `pr_ci_low`, `pr_ci_high`,
#'   `n`, `n_pos`, `n_neg`.
#' @export
pr_auc <- function(scores,
                   response,
                   direction = c("<", ">"),
                   n_boot = 1000L,
                   seed = 42L) {
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
    stop("PR-AUC icin her iki sinif da gerekli.", call. = FALSE)
  }

  oriented <- orient_score(scores, direction = direction)
  pr_obj <- PRROC::pr.curve(
    scores.class0 = oriented,
    weights.class0 = response
  )
  pr_val <- as.numeric(pr_obj$auc.davis.goadrich)
  pr_baseline <- n_pos / n

  pr_ci_low <- NA_real_
  pr_ci_high <- NA_real_

  if (n_boot > 0L) {
    pos_idx <- which(response == 1L)
    neg_idx <- which(response == 0L)

    boot_aucs <- withr::with_seed(seed, {
      vapply(seq_len(n_boot), function(i) {
        bi <- c(
          sample(pos_idx, size = n_pos, replace = TRUE),
          sample(neg_idx, size = n_neg, replace = TRUE)
        )
        pr_b <- PRROC::pr.curve(
          scores.class0 = oriented[bi],
          weights.class0 = response[bi]
        )
        as.numeric(pr_b$auc.davis.goadrich)
      }, numeric(1))
    })

    qs <- stats::quantile(boot_aucs, probs = c(0.025, 0.975), na.rm = TRUE, names = FALSE)
    pr_ci_low <- qs[1]
    pr_ci_high <- qs[2]
  }

  list(
    pr_auc = pr_val,
    pr_baseline = pr_baseline,
    pr_ci_low = pr_ci_low,
    pr_ci_high = pr_ci_high,
    n = n,
    n_pos = n_pos,
    n_neg = n_neg
  )
}
