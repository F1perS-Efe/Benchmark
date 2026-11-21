#' Veri setleri uzerinde agirlikli AUC ozeti (n ve ters-varyans)
#'
#' Her `signature` x `method` icin:
#' - **n-agirlikli** ortalama: `sum(n * auc) / sum(n)`
#' - **ters-varyans (IVW)** ortalama: `sum(w * auc) / sum(w)`, `w = 1 / auc_var`
#' - Heterojenlik: Cochran `Q` ve `I^2` (IVW agirliklariyla)
#'
#' Varsayilan kurallar:
#' - `independent == FALSE` (`derived_from`) satirlar ozetten cikarilir
#' - `var_zero == TRUE` (veya `auc_var` NA/0) satirlar ozetten cikarilir + uyari
#' - `skipped` veya sonlu olmayan AUC satirlar yok sayilir
#'
#' DerSimonian-Laird rastgele etkili model bu fonksiyonda yoktur (opsiyonel genisleme).
#'
#' @param results `run_benchmark()` cikti tablosu.
#' @param method Belirli skorlama yontemi; `NULL` ise tum yontemler ayri ozetlenir.
#' @param include_nonindependent `TRUE` ise bagimsiz-olmayan satirlar ozete girer.
#' @param exclude_var_zero `TRUE` (varsayilan) ise varyansi 0/NA olanlar cikarilir.
#'
#' @return `data.frame` ozet satirlari.
#' @export
summarize_auc <- function(results,
                          method = NULL,
                          include_nonindependent = FALSE,
                          exclude_var_zero = TRUE) {
  if (!is.data.frame(results) || nrow(results) == 0L) {
    stop("`results` bos olamaz.", call. = FALSE)
  }
  req <- c("signature", "method", "dataset", "auc", "n", "skipped")
  miss <- setdiff(req, names(results))
  if (length(miss) > 0L) {
    stop("`results` eksik sutunlar: ", paste(miss, collapse = ", "), call. = FALSE)
  }

  df <- results
  if (!is.null(method)) {
    df <- df[df$method == method, , drop = FALSE]
  }
  df <- df[!(df$skipped %in% TRUE) & is.finite(df$auc), , drop = FALSE]
  if (nrow(df) == 0L) {
    return(.summarize_auc_empty())
  }

  if (!"independent" %in% names(df)) {
    df$independent <- TRUE
  }
  if (!"var_zero" %in% names(df)) {
    df$var_zero <- if ("auc_var" %in% names(df)) {
      is.na(df$auc_var) | df$auc_var == 0
    } else {
      FALSE
    }
  }
  if (!"auc_var" %in% names(df)) {
    df$auc_var <- NA_real_
  }
  if (!"small_cohort" %in% names(df)) {
    df$small_cohort <- FALSE
  }

  split_key <- paste(df$signature, df$method, sep = "\r")
  groups <- split(df, split_key)

  rows <- lapply(groups, function(g) {
    n_total_rows <- nrow(g)
    indep_flag <- as.logical(g$independent)
    var_zero_flag <- as.logical(g$var_zero) |
      is.na(g$auc_var) |
      (!is.na(g$auc_var) & g$auc_var == 0)
    n_nonindep <- sum(!indep_flag, na.rm = TRUE)
    n_small <- sum(as.logical(g$small_cohort), na.rm = TRUE)

    use <- g
    if (!isTRUE(include_nonindependent)) {
      keep_indep <- is.na(indep_flag) | indep_flag
      use <- use[keep_indep, , drop = FALSE]
      use_var_zero <- var_zero_flag[keep_indep]
    } else {
      use_var_zero <- var_zero_flag
    }
    # var_zero sayimi: independent filtresi sonrasi kalan satirlar
    n_var_zero <- sum(use_var_zero, na.rm = TRUE)
    if (isTRUE(exclude_var_zero) && nrow(use) > 0L) {
      use <- use[!use_var_zero, , drop = FALSE]
    }

    notes <- character(0)
    if (!isTRUE(include_nonindependent) && n_nonindep > 0L) {
      notes <- c(notes, sprintf("%d bagimsiz-olmayan satir ozetten cikarildi", n_nonindep))
    }
    if (isTRUE(exclude_var_zero) && n_var_zero > 0L) {
      notes <- c(notes, sprintf("%d var_zero/NA varyans satir ozetten cikarildi", n_var_zero))
      warning(
        "Imza '", g$signature[1], "' / method '", g$method[1],
        "': varyansi 0 veya NA olan AUC'ler agirlikli ozetten cikarildi.",
        call. = FALSE
      )
    }
    if (n_small > 0L) {
      notes <- c(notes, sprintf("%d kucuk kohort uyarisi", n_small))
    }

    k <- nrow(use)
    if (k == 0L) {
      return(data.frame(
        signature = g$signature[1],
        method = g$method[1],
        n_datasets = 0L,
        n_total_rows = n_total_rows,
        auc_n_weighted = NA_real_,
        auc_ivw = NA_real_,
        Q = NA_real_,
        I2 = NA_real_,
        n_excluded_nonindependent = as.integer(n_nonindep),
        n_excluded_var_zero = as.integer(n_var_zero),
        n_small_cohort = as.integer(n_small),
        note = paste(notes, collapse = "; "),
        stringsAsFactors = FALSE
      ))
    }

    n_vec <- as.numeric(use$n)
    auc_vec <- as.numeric(use$auc)
    var_vec <- as.numeric(use$auc_var)

    auc_n <- sum(n_vec * auc_vec, na.rm = TRUE) / sum(n_vec, na.rm = TRUE)

    # IVW yalnizca pozitif sonlu varyanslarla (1/0 = Inf sorununu onle)
    ivw_ok <- is.finite(var_vec) & var_vec > 0 & is.finite(auc_vec)
    if (sum(ivw_ok) >= 1L) {
      w <- 1 / var_vec[ivw_ok]
      auc_ivw_vec <- auc_vec[ivw_ok]
      auc_ivw <- sum(w * auc_ivw_vec) / sum(w)
      k_ivw <- sum(ivw_ok)
      if (k_ivw >= 2L) {
        Q <- sum(w * (auc_ivw_vec - auc_ivw)^2)
        df_q <- k_ivw - 1L
        I2 <- if (is.finite(Q) && Q > 0) {
          max(0, (Q - df_q) / Q) * 100
        } else if (is.finite(Q) && Q == 0) {
          0
        } else {
          NA_real_
        }
      } else {
        Q <- NA_real_
        I2 <- NA_real_
      }
    } else {
      auc_ivw <- NA_real_
      Q <- NA_real_
      I2 <- NA_real_
    }

    data.frame(
      signature = g$signature[1],
      method = g$method[1],
      n_datasets = as.integer(k),
      n_total_rows = n_total_rows,
      auc_n_weighted = auc_n,
      auc_ivw = auc_ivw,
      Q = Q,
      I2 = I2,
      n_excluded_nonindependent = as.integer(n_nonindep),
      n_excluded_var_zero = as.integer(n_var_zero),
      n_small_cohort = as.integer(n_small),
      note = if (length(notes) == 0L) NA_character_ else paste(notes, collapse = "; "),
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

.summarize_auc_empty <- function() {
  data.frame(
    signature = character(),
    method = character(),
    n_datasets = integer(),
    n_total_rows = integer(),
    auc_n_weighted = numeric(),
    auc_ivw = numeric(),
    Q = numeric(),
    I2 = numeric(),
    n_excluded_nonindependent = integer(),
    n_excluded_var_zero = integer(),
    n_small_cohort = integer(),
    note = character(),
    stringsAsFactors = FALSE
  )
}
