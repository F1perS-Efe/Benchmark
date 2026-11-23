#' PR-AUC ve tabakali bootstrap %95 GA
#'
#' Skor once `orient_score()` ile "yuksek = yanitli" yapilir, sonra
#' `PRROC::pr.curve(scores.class0 = ..., weights.class0 = response)`
#' ile Davis-Goadrich PR-AUC hesaplanir. Yogun tie (bagli deger) nedeniyle
#' Davis-Goadrich interpolasyonu sonlu degilse `auc.integral` degerine
#' geri donulur (`tie_fraction` raporlanir).
#'
#' Taban cizgisi (prevalans): `n_pos / n`.
#' GA: pozitif ve negatif indeksler ayri ayri iadeli orneklenir
#' (`n_pos` ve `n_neg` sabit); yuzdelik (0.025, 0.975).
#' Kucuk pozitif/negatif sinif (`n_pos < 5` veya `n_neg < 5`) ve
#' yetersiz gecerli bootstrap tekrari (`n_boot_valid < 0.5 * n_boot`)
#' durumlarinda `small_positive_class` ve `ci_unreliable` bayraklari isaretlenir.
#'
#' @param scores Ornek skorlari.
#' @param response 0/1 yanit.
#' @param direction Imza yonu (`orient_score` icin).
#' @param n_boot Bootstrap tekrari. Varsayilan 1000; `0` ise CI `NA`.
#' @param seed RNG tohumu (`withr::with_seed`).
#'
#' @return Liste: `pr_auc`, `pr_baseline`, `pr_ci_low`, `pr_ci_high`,
#'   `n`, `n_pos`, `n_neg`, `tie_fraction`, `small_positive_class`,
#'   `n_boot_valid`, `ci_unreliable`.
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
    stop("`scores` and `response` must have the same length.", call. = FALSE)
  }

  ok <- is.finite(scores) & !is.na(response)
  scores <- scores[ok]
  response <- response[ok]

  n_pos <- sum(response == 1L)
  n_neg <- sum(response == 0L)
  n <- length(response)
  if (n_pos == 0L || n_neg == 0L) {
    stop("Both response classes are required for PR-AUC.", call. = FALSE)
  }

  oriented <- orient_score(scores, direction = direction)
  pr_baseline <- n_pos / n
  n_unique <- length(unique(oriented))
  tie_fraction <- if (n > 1L) 1 - (n_unique / n) else 0
  small_positive_class <- (n_pos < 5L) || (n_neg < 5L)

  if (n_unique < 2L || !is.finite(stats::sd(oriented)) || stats::sd(oriented) == 0) {
    warning("Constant or zero-variance scores in pr_auc(); returning NA.", call. = FALSE)
    return(list(
      pr_auc = NA_real_,
      pr_baseline = pr_baseline,
      pr_ci_low = NA_real_,
      pr_ci_high = NA_real_,
      n = n,
      n_pos = n_pos,
      n_neg = n_neg,
      tie_fraction = tie_fraction,
      small_positive_class = small_positive_class,
      n_boot_valid = 0L,
      ci_unreliable = TRUE
    ))
  }

  .extract_pr <- function(s_vec, w_vec) {
    pr_obj <- tryCatch(
      PRROC::pr.curve(scores.class0 = s_vec, weights.class0 = w_vec),
      error = function(e) NULL
    )
    if (is.null(pr_obj)) {
      return(NA_real_)
    }
    val <- as.numeric(pr_obj$auc.davis.goadrich)
    if (!is.finite(val) && !is.null(pr_obj$auc.integral)) {
      val <- as.numeric(pr_obj$auc.integral)
    }
    if (is.finite(val)) val else NA_real_
  }

  pr_val <- .extract_pr(oriented, response)
  if (!is.finite(pr_val)) {
    warning("PRROC::pr.curve returned non-finite PR-AUC; returning NA.", call. = FALSE)
    return(list(
      pr_auc = NA_real_,
      pr_baseline = pr_baseline,
      pr_ci_low = NA_real_,
      pr_ci_high = NA_real_,
      n = n,
      n_pos = n_pos,
      n_neg = n_neg,
      tie_fraction = tie_fraction,
      small_positive_class = small_positive_class,
      n_boot_valid = 0L,
      ci_unreliable = TRUE
    ))
  }

  pr_ci_low <- NA_real_
  pr_ci_high <- NA_real_
  n_boot_valid <- 0L
  ci_unreliable <- small_positive_class || (tie_fraction > 0.5)

  if (n_boot > 0L) {
    pos_idx <- which(response == 1L)
    neg_idx <- which(response == 0L)

    boot_aucs <- withr::with_seed(seed, {
      vapply(seq_len(n_boot), function(i) {
        bi <- c(
          sample(pos_idx, size = n_pos, replace = TRUE),
          sample(neg_idx, size = n_neg, replace = TRUE)
        )
        if (length(unique(oriented[bi])) < 2L) {
          return(NA_real_)
        }
        .extract_pr(oriented[bi], response[bi])
      }, numeric(1))
    })

    valid_boots <- boot_aucs[is.finite(boot_aucs)]
    n_boot_valid <- length(valid_boots)
    if (n_boot_valid > 0L) {
      qs <- stats::quantile(valid_boots, probs = c(0.025, 0.975), na.rm = TRUE, names = FALSE)
      pr_ci_low <- qs[1]
      pr_ci_high <- qs[2]
    }
    if (n_boot_valid < max(10L, as.integer(0.5 * n_boot))) {
      ci_unreliable <- TRUE
    }
  }

  list(
    pr_auc = pr_val,
    pr_baseline = pr_baseline,
    pr_ci_low = pr_ci_low,
    pr_ci_high = pr_ci_high,
    n = n,
    n_pos = n_pos,
    n_neg = n_neg,
    tie_fraction = tie_fraction,
    small_positive_class = small_positive_class,
    n_boot_valid = as.integer(n_boot_valid),
    ci_unreliable = isTRUE(ci_unreliable)
  )
}
