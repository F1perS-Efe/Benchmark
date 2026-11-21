#' Rastgele gen null permütasyon testi (vektorlestirilmis AUC)
#'
#' Bulunan imza gen sayisi kadar rastgele gen seti skorlanir; gozlenen AUC
#' ile null dagilimdan tek yonlu ampirik p hesaplanir:
#' `(sum(null_aucs >= obs_auc) + 1) / (B + 1)`.
#'
#' Null AUC'ler `auc_mw()` ile hesaplanir (dongude `pROC` yok).
#'
#' @param mat Gen x ornek ifade matrisi (hazirlanmis).
#' @param response 0/1 yanit (ornek sirasi `colnames(mat)` ile hizali veya isimli).
#' @param n_genes Cekilecek rastgele gen sayisi (kapsama sonrasi bulunan).
#' @param exclude_genes Havuzdan cikarilacak genler (imza genleri).
#' @param method Skorlama yontemi (`score_signature`).
#' @param direction Imza yonu.
#' @param n_perm Permütasyon sayisi B. Varsayilan 200.
#' @param seed RNG tohumu.
#' @param match_expression TRUE ise imza genlerinin ortalama ifadesine yakin
#'   quantile-bin eslestirilmis ornekleme yapar.
#' @param signature_genes `match_expression = TRUE` icin imza genleri
#'   (ifade ortalamasi eslestirme referansi). `NULL` ise `exclude_genes`.
#' @param observed_scores Opsiyonel: gozlenen imza skorlari (isimli).
#'   Verilmezse `signature_genes` ile skorlanir.
#' @param min_genes `score_signature` esigi.
#'
#' @return Liste: `obs_auc`, `null_aucs`, `null_p`, `n_perm`, `n_genes`,
#'   `pool_size`, `seed`.
#' @export
null_permutation_test <- function(mat,
                                  response,
                                  n_genes,
                                  exclude_genes = character(0),
                                  method = c("mean", "zscore", "median", "rank"),
                                  direction = c("<", ">"),
                                  n_perm = 200L,
                                  seed = 42L,
                                  match_expression = FALSE,
                                  signature_genes = NULL,
                                  observed_scores = NULL,
                                  min_genes = 3L) {
  method <- match.arg(method)
  direction <- match.arg(direction)
  mat <- as.matrix(mat)
  if (is.null(rownames(mat)) || is.null(colnames(mat))) {
    stop("`mat` satir/sutun adlari zorunludur.", call. = FALSE)
  }

  # Yaniti orneklerle hizala
  if (!is.null(names(response))) {
    response <- as.integer(response[colnames(mat)])
  } else {
    response <- as.integer(response)
    if (length(response) != ncol(mat)) {
      stop("`response` uzunlugu `ncol(mat)` ile eslesmeli.", call. = FALSE)
    }
  }

  exclude_up <- unique(toupper(as.character(exclude_genes)))
  rows_up <- toupper(rownames(mat))

  # Sifir varyans / tamamen NA genleri cikar
  gene_ok <- apply(mat, 1L, function(x) {
    x <- x[is.finite(x)]
    length(x) > 1L && stats::sd(x) > 0
  })
  pool_idx <- which(gene_ok & !(rows_up %in% exclude_up))
  pool_genes <- rownames(mat)[pool_idx]
  pool_size <- length(pool_genes)

  if (pool_size < n_genes) {
    warning(
      "Null havuz boyutu (", pool_size, ") n_genes (", n_genes,
      ") altinda; null testi atlandi.",
      call. = FALSE
    )
    return(list(
      obs_auc = NA_real_,
      null_aucs = numeric(0),
      null_p = NA_real_,
      n_perm = as.integer(n_perm),
      n_genes = as.integer(n_genes),
      pool_size = pool_size,
      seed = seed
    ))
  }

  # Gozlenen AUC
  if (is.null(observed_scores)) {
    sig_use <- if (is.null(signature_genes)) exclude_genes else signature_genes
    scored <- score_signature(mat, genes = sig_use, method = method, min_genes = min_genes)
    if (isTRUE(scored$skipped) || is.null(scored$scores)) {
      stop("Gozlenen imza skorlanamadi: ", scored$skip_reason, call. = FALSE)
    }
    observed_scores <- scored$scores
  }

  # Skor / yanit ortak ornekleri
  if (!is.null(names(observed_scores))) {
    common <- intersect(names(observed_scores), colnames(mat))
    obs_scores <- as.numeric(observed_scores[common])
    resp_use <- response[match(common, colnames(mat))]
  } else {
    obs_scores <- as.numeric(observed_scores)
    resp_use <- response
  }
  ok <- is.finite(obs_scores) & !is.na(resp_use)
  obs_scores <- obs_scores[ok]
  resp_use <- resp_use[ok]
  obs_oriented <- orient_score(obs_scores, direction = direction)
  obs_auc <- auc_mw(obs_oriented, resp_use)

  # Eslestirme: imza gen ortalama ifadeleri -> quantile bin
  sig_ref <- if (is.null(signature_genes)) exclude_genes else signature_genes
  sig_ref_up <- toupper(sig_ref)
  sig_means <- NULL
  if (isTRUE(match_expression)) {
    hit <- match(sig_ref_up, rows_up)
    hit <- hit[!is.na(hit)]
    if (length(hit) == 0L) {
      warning("match_expression: imza geni bulunamadi; rastgele orneklemeye dusuldu.", call. = FALSE)
      match_expression <- FALSE
    } else {
      sig_means <- rowMeans(mat[hit, , drop = FALSE], na.rm = TRUE)
    }
  }

  pool_means <- rowMeans(mat[pool_idx, , drop = FALSE], na.rm = TRUE)
  # Tum genlerin quantile binleri (havuz + referans olcegi icin)
  all_means <- rowMeans(mat, na.rm = TRUE)
  brks <- unique(stats::quantile(all_means[is.finite(all_means)],
                                 probs = seq(0, 1, length.out = 6),
                                 na.rm = TRUE, names = FALSE))
  if (length(brks) < 2L) {
    match_expression <- FALSE
  } else {
    pool_bins <- if (isTRUE(match_expression)) {
      as.integer(cut(pool_means, breaks = brks, include.lowest = TRUE, labels = FALSE))
    } else {
      NULL
    }
    sig_bins <- if (isTRUE(match_expression)) {
      as.integer(cut(sig_means, breaks = brks, include.lowest = TRUE, labels = FALSE))
    } else {
      NULL
    }
  }

  null_aucs <- withr::with_seed(seed, {
    vapply(seq_len(n_perm), function(b) {
      if (isTRUE(match_expression) && !is.null(sig_bins)) {
        picked <- character(0)
        for (bn in sig_bins) {
          cand <- pool_genes[pool_bins == bn & !(pool_genes %in% picked)]
          if (length(cand) == 0L) {
            cand <- setdiff(pool_genes, picked)
          }
          picked <- c(picked, sample(cand, size = 1L))
        }
        # Fazla/eksik uzunluk guvence
        if (length(picked) < n_genes) {
          extra <- setdiff(pool_genes, picked)
          picked <- c(picked, sample(extra, size = n_genes - length(picked)))
        }
        picked <- picked[seq_len(n_genes)]
      } else {
        picked <- sample(pool_genes, size = n_genes, replace = FALSE)
      }

      sc <- score_signature(mat, genes = picked, method = method, min_genes = min_genes)
      if (isTRUE(sc$skipped) || is.null(sc$scores)) {
        return(NA_real_)
      }
      # Ayni ornek sirasi / ortak
      if (!is.null(names(sc$scores)) && !is.null(names(obs_scores))) {
        # obs_scores isimleri common'dan geldi; mat colnames kullan
        nm <- names(sc$scores)
        idx <- match(nm, colnames(mat))
        r <- response[idx]
        keep <- is.finite(sc$scores) & !is.na(r)
        s_or <- orient_score(sc$scores[keep], direction = direction)
        auc_mw(s_or, r[keep])
      } else {
        s_or <- orient_score(sc$scores, direction = direction)
        # response mat sutun sirasinda; skor isimleri colnames alt kumesiyse hizala
        if (!is.null(names(sc$scores))) {
          r <- response[match(names(sc$scores), colnames(mat))]
        } else {
          r <- response
        }
        keep <- is.finite(s_or) & !is.na(r)
        auc_mw(s_or[keep], r[keep])
      }
    }, numeric(1))
  })

  null_aucs <- null_aucs[is.finite(null_aucs)]
  B <- length(null_aucs)
  if (B == 0L || !is.finite(obs_auc)) {
    null_p <- NA_real_
  } else {
    k <- sum(null_aucs >= obs_auc)
    null_p <- (k + 1) / (B + 1)
  }

  list(
    obs_auc = obs_auc,
    null_aucs = null_aucs,
    null_p = null_p,
    n_perm = as.integer(n_perm),
    n_genes = as.integer(n_genes),
    pool_size = pool_size,
    seed = seed
  )
}
