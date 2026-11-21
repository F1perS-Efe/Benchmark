#' Rastgele gen null permutasyon testi (vektorlestirilmis AUC)
#'
#' Bulunan imza gen sayisi kadar rastgele gen seti skorlanir; gozlenen AUC
#' ile null dagilimdan tek yonlu ampirik p hesaplanir:
#' `(sum(null_aucs >= obs_auc) + 1) / (B + 1)`.
#'
#' Null AUC'ler `auc_mw()` ile hesaplanir (dongude `pROC` yok).
#' `zscore` / `rank` icin havuz satirlari donguden once bir kez donusturulur;
#' `mean` / `median` icin dongu ici `match_genes()` tekrari yoktur.
#'
#' Not: `rank` kohort-ici ornekler arasi goreli siralamadir (singscore degil);
#' `zscore` basit ortalama gen z-skorudur (`sqrt(k)` olcekli GSVA degil).
#'
#' @param mat Gen x ornek ifade matrisi (hazirlanmis).
#' @param response 0/1 yanit (ornek sirasi `colnames(mat)` ile hizali veya isimli).
#' @param n_genes Cekilecek rastgele gen sayisi (kapsama sonrasi bulunan).
#' @param exclude_genes Havuzdan cikarilacak genler (imza genleri).
#' @param method Skorlama yontemi (`score_signature`).
#' @param direction Imza yonu.
#' @param n_perm Permutasyon sayisi B. Varsayilan 1000.
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
                                  n_perm = 1000L,
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

  # P0-3: NA yanitlari (SD vb.) olceklemeden once cikar
  keep_y <- !is.na(response)
  if (!all(keep_y)) {
    mat <- mat[, keep_y, drop = FALSE]
    response <- response[keep_y]
  }

  exclude_up <- unique(toupper(as.character(exclude_genes)))
  rows_up <- toupper(rownames(mat))

  # P1-1: vektorlestirilmis satir SD (apply yerine)
  row_n <- rowSums(is.finite(mat))
  row_mean <- rowMeans(mat, na.rm = TRUE)
  row_sd <- sqrt(
    rowSums((mat - row_mean)^2, na.rm = TRUE) / pmax(row_n - 1L, 1L)
  )
  gene_ok <- row_n > 1L & is.finite(row_sd) & row_sd > 0
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
  sig_ref_up <- unique(toupper(as.character(sig_ref)))
  sig_means <- NULL
  if (isTRUE(match_expression)) {
    hit <- match(sig_ref_up, rows_up)
    hit <- hit[!is.na(hit)]
    if (length(hit) == 0L) {
      warning("match_expression: imza geni bulunamadi; rastgele orneklemeye dusuldu.", call. = FALSE)
      match_expression <- FALSE
    } else {
      sig_means <- rowMeans(mat[hit, , drop = FALSE], na.rm = TRUE)
      sig_means <- sig_means[is.finite(sig_means)]
      if (length(sig_means) == 0L) {
        warning("match_expression: imza ortalama ifadesi NA; rastgele orneklemeye dusuldu.", call. = FALSE)
        match_expression <- FALSE
      }
    }
  }

  pool_means <- rowMeans(mat[pool_idx, , drop = FALSE], na.rm = TRUE)
  # Tum genlerin quantile binleri (havuz + referans olcegi icin)
  all_means <- rowMeans(mat, na.rm = TRUE)
  brks <- unique(stats::quantile(all_means[is.finite(all_means)],
                                 probs = seq(0, 1, length.out = 6),
                                 na.rm = TRUE, names = FALSE))
  pool_bins <- NULL
  sig_bins <- NULL
  if (length(brks) < 2L) {
    match_expression <- FALSE
  } else if (isTRUE(match_expression)) {
    pool_bins <- as.integer(cut(pool_means, breaks = brks, include.lowest = TRUE, labels = FALSE))
    sig_bins <- as.integer(cut(sig_means, breaks = brks, include.lowest = TRUE, labels = FALSE))
    sig_bins <- sig_bins[!is.na(sig_bins)]
    if (length(sig_bins) == 0L) {
      warning("match_expression: sig_bins NA; rastgele orneklemeye dusuldu.", call. = FALSE)
      match_expression <- FALSE
      pool_bins <- NULL
      sig_bins <- NULL
    }
  }

  # Havuz matrisini donguden once hazirla (zscore/rank satir donusumu bir kez)
  pool_mat <- mat[pool_idx, , drop = FALSE]
  rownames(pool_mat) <- pool_genes
  if (method == "zscore") {
    # P1-1: vektorlestirilmis z-score
    pm_n <- rowSums(is.finite(pool_mat))
    pm_mean <- rowMeans(pool_mat, na.rm = TRUE)
    pm_sd <- sqrt(
      rowSums((pool_mat - pm_mean)^2, na.rm = TRUE) / pmax(pm_n - 1L, 1L)
    )
    score_mat <- (pool_mat - pm_mean) / pm_sd
    score_mat[!is.finite(score_mat)] <- NA_real_
  } else if (method == "rank") {
    score_mat <- t(apply(pool_mat, 1L, function(x) {
      m <- sum(!is.na(x))
      if (m <= 1L) {
        return(rep(NA_real_, length(x)))
      }
      r <- rank(x, ties.method = "average", na.last = "keep")
      (r - 1) / (m - 1)
    }))
    colnames(score_mat) <- colnames(pool_mat)
    rownames(score_mat) <- rownames(pool_mat)
  } else {
    score_mat <- pool_mat
  }

  null_aucs <- withr::with_seed(seed, {
    vapply(seq_len(n_perm), function(b) {
      if (isTRUE(match_expression) && !is.null(sig_bins)) {
        picked <- character(0)
        for (bn in sig_bins) {
          cand_idx <- which(!is.na(pool_bins) & pool_bins == bn & !(pool_genes %in% picked))
          cand <- pool_genes[cand_idx]
          if (length(cand) == 0L) {
            cand <- setdiff(pool_genes, picked)
          }
          if (length(cand) == 0L) {
            return(NA_real_)
          }
          picked <- c(picked, sample(cand, size = 1L))
        }
        # Fazla/eksik uzunluk guvence
        if (length(picked) < n_genes) {
          extra <- setdiff(pool_genes, picked)
          need <- n_genes - length(picked)
          if (length(extra) < need) {
            return(NA_real_)
          }
          picked <- c(picked, sample(extra, size = need))
        }
        picked <- picked[seq_len(n_genes)]
      } else {
        picked <- sample(pool_genes, size = n_genes, replace = FALSE)
      }

      row_idx <- match(picked, rownames(score_mat))
      if (anyNA(row_idx)) {
        return(NA_real_)
      }
      sub <- score_mat[row_idx, , drop = FALSE]
      if (method == "median") {
        scores <- apply(sub, 2L, function(x) {
          if (all(is.na(x))) {
            NA_real_
          } else {
            stats::median(x, na.rm = TRUE)
          }
        })
      } else {
        scores <- colMeans(sub, na.rm = TRUE)
        scores[colSums(!is.na(sub)) == 0L] <- NA_real_
      }
      names(scores) <- colnames(score_mat)

      keep <- is.finite(scores) & !is.na(response)
      if (sum(keep) < 2L) {
        return(NA_real_)
      }
      s_or <- orient_score(scores[keep], direction = direction)
      auc_mw(s_or, response[keep])
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
