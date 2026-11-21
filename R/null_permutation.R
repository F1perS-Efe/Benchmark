#' Rastgele gen null permutasyon testi (vektorlestirilmis AUC)
#'
#' Bulunan imza gen sayisi kadar rastgele gen seti skorlanir; gozlenen AUC
#' ile null dagilimdan tek yonlu ampirik p hesaplanir:
#' `(sum(null_aucs >= obs_auc) + 1) / (B + 1)`.
#'
#' Null AUC'ler `auc_mw()` ile hesaplanir (dongude `pROC` yok).
#' `zscore` / `rank` icin havuz satirlari donguden once bir kez donusturulur;
#' `ss_rank` icin sutun bazli dikey sira matrisi (`rank_mat`) ve `N_j`
#' donguden once bir kez hesaplanir; `mean` / `median` icin dongu ici
#' `match_genes()` tekrari yoktur.
#'
#' Not: `rank` kohort-ici ornekler arasi goreli siralamadir (singscore degil);
#' `ss_rank` kohorttan bagimsiz dikey tek-ornek siralamadir (singscore
#' normalizasyonu ile uyumlu saf base R); `zscore` basit ortalama gen
#' z-skorudur (`sqrt(k)` olcekli GSVA degil).
#'
#' Isaretli imzalar icin `n_genes_down > 0L` iken her permutasyonda birbirinden
#' ayrik `n_genes` (up) + `n_genes_down` (down) rastgele gen cekilir;
#' skor `up - down` (`ss_rank` icin `s_up + s_down`) olur.
#' Yalnizca down kolu (`n_genes = 0`, `n_genes_down >= 1`) icin skor `-down`
#' (`ss_rank` icin azalan sira `s_down`) uretilir.
#'
#' @param mat Gen x ornek ifade matrisi (hazirlanmis).
#' @param response 0/1 yanit (ornek sirasi `colnames(mat)` ile hizali veya isimli).
#' @param n_genes Cekilecek rastgele up / tek yonlu gen sayisi (kapsama sonrasi).
#'   Down-only isaretli null icin `0L` olabilir.
#' @param exclude_genes Havuzdan cikarilacak genler (imza genleri).
#' @param method Skorlama yontemi (`score_signature`).
#' @param direction Imza yonu.
#' @param n_perm Permutasyon sayisi B. Varsayilan 1000.
#' @param seed RNG tohumu.
#' @param match_expression TRUE ise imza genlerinin ortalama ifadesine yakin
#'   quantile-bin eslestirilmis ornekleme yapar.
#' @param signature_genes `match_expression = TRUE` icin up / tek yonlu imza
#'   genleri (ifade ortalamasi eslestirme referansi). `NULL` ise `exclude_genes`
#'   (up kolu yoksa eslestirme atlanir).
#' @param signature_genes_down `match_expression = TRUE` ve `n_genes_down > 0`
#'   iken down kolu ifade-bin referansi. `NULL` ise down kolu guvenli rastgele
#'   fallback kullanir.
#' @param observed_scores Opsiyonel: gozlenen imza skorlari (isimli).
#'   Verilmezse `signature_genes` ile skorlanir.
#' @param min_genes `score_signature` esigi.
#' @param n_genes_down Isaretli null icin down kolu gen sayisi. Varsayilan 0.
#'
#' @return Liste: `obs_auc`, `null_aucs`, `null_p`, `n_perm`, `n_genes`,
#'   `pool_size`, `seed`.
#' @export
null_permutation_test <- function(mat,
                                  response,
                                  n_genes,
                                  exclude_genes = character(0),
                                  method = c("mean", "zscore", "median", "rank", "ss_rank"),
                                  direction = c("<", ">"),
                                  n_perm = 1000L,
                                  seed = 42L,
                                  match_expression = FALSE,
                                  signature_genes = NULL,
                                  observed_scores = NULL,
                                  min_genes = 3L,
                                  n_genes_down = 0L,
                                  signature_genes_down = NULL) {
  method <- match.arg(method)
  direction <- match.arg(direction)
  mat <- as.matrix(mat)
  n_genes <- as.integer(n_genes)
  n_genes_down <- as.integer(n_genes_down)
  if (length(n_genes) != 1L || is.na(n_genes) || n_genes < 0L) {
    stop("`n_genes` must be a non-negative integer.", call. = FALSE)
  }
  if (length(n_genes_down) != 1L || is.na(n_genes_down) || n_genes_down < 0L) {
    stop("`n_genes_down` must be a non-negative integer.", call. = FALSE)
  }
  if (n_genes < 1L && n_genes_down < 1L) {
    stop("At least one of `n_genes` and `n_genes_down` must be positive.", call. = FALSE)
  }
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
  need_genes <- n_genes + n_genes_down

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

  if (pool_size < need_genes) {
    warning(
      "Null havuz boyutu (", pool_size, ") n_genes+n_genes_down (",
      need_genes, ") altinda; null testi atlandi.",
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

  # Eslestirme: imza gen ortalama ifadeleri -> quantile bin (up ve down ayri)
  pool_means <- rowMeans(mat[pool_idx, , drop = FALSE], na.rm = TRUE)
  all_means <- rowMeans(mat, na.rm = TRUE)
  brks <- unique(stats::quantile(all_means[is.finite(all_means)],
                                 probs = seq(0, 1, length.out = 6),
                                 na.rm = TRUE, names = FALSE))
  pool_bins <- NULL
  sig_bins <- NULL
  sig_bins_down <- NULL
  match_up <- isTRUE(match_expression) && n_genes > 0L
  match_down <- isTRUE(match_expression) &&
    n_genes_down > 0L &&
    !is.null(signature_genes_down) &&
    length(signature_genes_down) > 0L

  if (length(brks) < 2L) {
    match_up <- FALSE
    match_down <- FALSE
  } else {
    pool_bins <- as.integer(cut(pool_means, breaks = brks, include.lowest = TRUE, labels = FALSE))
  }

  if (isTRUE(match_up)) {
    sig_ref <- if (is.null(signature_genes)) exclude_genes else signature_genes
    sig_ref_up <- unique(toupper(as.character(sig_ref)))
    hit <- match(sig_ref_up, rows_up)
    hit <- hit[!is.na(hit)]
    if (length(hit) == 0L) {
      warning("match_expression: up imza geni bulunamadi; up rastgele orneklemeye dusuldu.", call. = FALSE)
      match_up <- FALSE
    } else {
      sig_means <- rowMeans(mat[hit, , drop = FALSE], na.rm = TRUE)
      sig_means <- sig_means[is.finite(sig_means)]
      if (length(sig_means) == 0L) {
        warning("match_expression: up ortalama ifadesi NA; up rastgele orneklemeye dusuldu.", call. = FALSE)
        match_up <- FALSE
      } else {
        sig_bins <- as.integer(cut(sig_means, breaks = brks, include.lowest = TRUE, labels = FALSE))
        sig_bins <- sig_bins[!is.na(sig_bins)]
        if (length(sig_bins) == 0L) {
          warning("match_expression: sig_bins NA; up rastgele orneklemeye dusuldu.", call. = FALSE)
          match_up <- FALSE
          sig_bins <- NULL
        }
      }
    }
  }

  if (isTRUE(match_down)) {
    sig_ref_dn <- unique(toupper(as.character(signature_genes_down)))
    hit_dn <- match(sig_ref_dn, rows_up)
    hit_dn <- hit_dn[!is.na(hit_dn)]
    if (length(hit_dn) == 0L) {
      warning("match_expression: down imza geni bulunamadi; down rastgele orneklemeye dusuldu.", call. = FALSE)
      match_down <- FALSE
    } else {
      sig_means_dn <- rowMeans(mat[hit_dn, , drop = FALSE], na.rm = TRUE)
      sig_means_dn <- sig_means_dn[is.finite(sig_means_dn)]
      if (length(sig_means_dn) == 0L) {
        warning("match_expression: down ortalama ifadesi NA; down rastgele orneklemeye dusuldu.", call. = FALSE)
        match_down <- FALSE
      } else {
        sig_bins_down <- as.integer(
          cut(sig_means_dn, breaks = brks, include.lowest = TRUE, labels = FALSE)
        )
        sig_bins_down <- sig_bins_down[!is.na(sig_bins_down)]
        if (length(sig_bins_down) == 0L) {
          warning("match_expression: sig_bins_down NA; down rastgele orneklemeye dusuldu.", call. = FALSE)
          match_down <- FALSE
          sig_bins_down <- NULL
        }
      }
    }
  }

  # Havuz / sira matrisi donguden once
  pool_mat <- mat[pool_idx, , drop = FALSE]
  rownames(pool_mat) <- pool_genes
  rank_mat <- NULL
  N_j <- NULL
  if (method == "ss_rank") {
    rank_mat <- apply(mat, 2L, function(x) {
      rank(x, ties.method = "average", na.last = "keep")
    })
    if (is.null(dim(rank_mat))) {
      rank_mat <- matrix(rank_mat, nrow = nrow(mat), ncol = ncol(mat))
    }
    rownames(rank_mat) <- rownames(mat)
    colnames(rank_mat) <- colnames(mat)
    N_j <- colSums(!is.na(mat))
    score_mat <- NULL
  } else if (method == "zscore") {
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

  .sample_matched_arm <- function(n_need, bins, already_picked) {
    if (n_need < 1L) {
      return(character(0))
    }
    avail <- setdiff(pool_genes, already_picked)
    if (length(avail) < n_need) {
      return(NULL)
    }
    picked <- character(0)
    for (bn in bins) {
      if (length(picked) >= n_need) {
        break
      }
      cand_idx <- which(
        !is.na(pool_bins) &
          pool_bins == bn &
          pool_genes %in% avail &
          !(pool_genes %in% picked)
      )
      cand <- pool_genes[cand_idx]
      if (length(cand) == 0L) {
        cand <- setdiff(avail, picked)
      }
      if (length(cand) == 0L) {
        return(NULL)
      }
      picked <- c(picked, sample(cand, size = 1L))
    }
    if (length(picked) < n_need) {
      extra <- setdiff(avail, picked)
      need <- n_need - length(picked)
      if (length(extra) < need) {
        return(NULL)
      }
      picked <- c(picked, sample(extra, size = need))
    }
    picked[seq_len(n_need)]
  }

  .pick_null_genes <- function() {
    if (n_genes > 0L) {
      if (isTRUE(match_up) && !is.null(sig_bins)) {
        picked_up <- .sample_matched_arm(n_genes, sig_bins, character(0))
        if (is.null(picked_up)) {
          return(NULL)
        }
      } else {
        picked_up <- sample(pool_genes, size = n_genes, replace = FALSE)
      }
    } else {
      picked_up <- character(0)
    }

    if (n_genes_down > 0L) {
      if (isTRUE(match_down) && !is.null(sig_bins_down)) {
        picked_down <- .sample_matched_arm(n_genes_down, sig_bins_down, picked_up)
        if (is.null(picked_down)) {
          return(NULL)
        }
      } else {
        rem <- setdiff(pool_genes, picked_up)
        if (length(rem) < n_genes_down) {
          return(NULL)
        }
        picked_down <- sample(rem, size = n_genes_down, replace = FALSE)
      }
    } else {
      picked_down <- character(0)
    }
    list(up = picked_up, down = picked_down)
  }

  .arm_scores <- function(gene_names) {
    if (length(gene_names) == 0L) {
      return(NULL)
    }
    row_idx <- match(gene_names, rownames(score_mat))
    if (anyNA(row_idx)) {
      return(NULL)
    }
    sub <- score_mat[row_idx, , drop = FALSE]
    if (method == "median") {
      out <- apply(sub, 2L, function(x) {
        if (all(is.na(x))) NA_real_ else stats::median(x, na.rm = TRUE)
      })
    } else {
      out <- colMeans(sub, na.rm = TRUE)
      out[colSums(!is.na(sub)) == 0L] <- NA_real_
    }
    out
  }

  .scores_from_pick <- function(picked) {
    if (method == "ss_rank") {
      has_up <- length(picked$up) > 0L
      has_down <- length(picked$down) > 0L
      s_up <- if (has_up) {
        .ss_rank_arm(rank_mat, picked$up, N_j, reverse = FALSE)
      } else {
        rep(0, ncol(mat))
      }
      s_down <- if (has_down) {
        .ss_rank_arm(rank_mat, picked$down, N_j, reverse = TRUE)
      } else {
        rep(0, ncol(mat))
      }
      scores <- if (has_up && has_down) {
        s_up + s_down
      } else if (has_up) {
        s_up
      } else if (has_down) {
        s_down
      } else {
        return(NULL)
      }
      names(scores) <- colnames(mat)
      return(scores)
    }

    scores_up <- .arm_scores(picked$up)
    scores_dn <- .arm_scores(picked$down)
    if (!is.null(scores_up) && !is.null(scores_dn)) {
      scores <- scores_up - scores_dn
    } else if (!is.null(scores_up)) {
      scores <- scores_up
    } else if (!is.null(scores_dn)) {
      scores <- -scores_dn
    } else {
      return(NULL)
    }
    names(scores) <- colnames(score_mat)
    scores
  }

  null_aucs <- withr::with_seed(seed, {
    vapply(seq_len(n_perm), function(b) {
      picked <- .pick_null_genes()
      if (is.null(picked)) {
        return(NA_real_)
      }
      scores <- .scores_from_pick(picked)
      if (is.null(scores)) {
        return(NA_real_)
      }
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
