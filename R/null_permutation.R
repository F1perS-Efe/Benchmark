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
#' Surekli agirlikli (`weights`) imzalar icin her rastgele cekilis ayni pozitif ve
#' negatif katsayi buyukluk dagilimiyla (`w_up`, `w_dn`) agirliklandirilir.
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
#' @param match_expression TRUE ise imza genlerinin hem ortalama ifadesine (5 quantile)
#'   hem de gen dispersiyonuna/SD'sine (3 tertile) yakin 2D-bin eslestirilmis
#'   ornekleme yapar.
#' @param signature_genes `match_expression = TRUE` icin up / tek yonlu imza
#'   genleri (ifade ortalamasi/SD eslestirme referansi). `NULL` ise `exclude_genes`
#'   (up kolu yoksa eslestirme atlanir).
#' @param signature_genes_down `match_expression = TRUE` ve `n_genes_down > 0`
#'   iken down kolu ifade-bin referansi. `NULL` ise down kolu guvenli rastgele
#'   fallback kullanir.
#' @param observed_scores Opsiyonel: gozlenen imza skorlari (isimli).
#'   Verilmezse `signature_genes` ile skorlanir.
#' @param min_genes `score_signature` esigi.
#' @param n_genes_down Isaretli null icin down kolu gen sayisi. Varsayilan 0.
#' @param weights Opsiyonel isimli sayisal vektor: surekli gen agirliklari
#'   (orn. Ayers/Cristescu `T_Cell_Inflamed_GEP` regresyon katsayilari).
#'
#' @return Liste: `obs_auc`, `null_aucs`, `null_p`, `n_perm`, `n_genes`,
#'   `n_genes_effective`, `pool_size`, `universe_size`, `partial_background`,
#'   `background_type`, `seed`.
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
                                  signature_genes_down = NULL,
                                  weights = NULL) {
  method <- match.arg(method)
  direction <- match.arg(direction)
  if (!is.matrix(mat)) {
    mat <- as.matrix(mat)
  }
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
    stop("`mat` rownames and colnames are required.", call. = FALSE)
  }

  ostat <- attr(mat, "ortholog_stats")
  bg_expr <- attr(mat, "transcriptome_bg_expr")
  partial_background <- FALSE
  if (is.matrix(bg_expr) && ncol(bg_expr) == ncol(mat)) {
    mat <- bg_expr
  } else if (!is.null(ostat) && isTRUE(ostat$is_partial_mapping)) {
    partial_background <- TRUE
    warning(
      sprintf(
        paste0(
          "Null permutation on partial ortholog-mapped matrix (mapped=%d of %d ",
          "input genes, retention=%.1f%%) without preserved background context: ",
          "null sampling pool is restricted to mapped ortholog subset and does ",
          "not represent full transcriptome background."
        ),
        ostat$n_mapped, ostat$n_input, 100 * ostat$retention_rate
      ),
      call. = FALSE
    )
  }
  universe_size <- nrow(mat)
  background_type <- if (isTRUE(partial_background)) {
    "partial_ortholog_subset"
  } else if (universe_size < 5000L) {
    "subsampled_panel_cohort"
  } else {
    "full_transcriptome"
  }

  # Yaniti orneklerle hizala (zscore/rank tum kohortta hesaplanir, NA yanit AUC adiminda dislanir)
  if (!is.null(names(response))) {
    response <- as.integer(response[colnames(mat)])
  } else {
    response <- as.integer(response)
    if (length(response) != ncol(mat)) {
      stop("`response` length must match `ncol(mat)`.", call. = FALSE)
    }
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

  # R2: zscore/rank sirasinda sifir-SD imza genleri dusuyorsa ayni efektif sayiyi cek
  n_genes_eff <- n_genes
  n_genes_down_eff <- n_genes_down
  if (method %in% c("zscore", "rank")) {
    if (n_genes > 0L && !is.null(signature_genes) && length(signature_genes) > 0L) {
      m_up <- tryCatch(match_genes(mat, signature_genes, use_aliases = TRUE)$found, error = function(e) character(0))
      hit_up_ok <- match(m_up, rownames(mat))
      hit_up_ok <- hit_up_ok[!is.na(hit_up_ok) & gene_ok[hit_up_ok]]
      if (length(hit_up_ok) > 0L && length(hit_up_ok) < n_genes_eff) {
        n_genes_eff <- length(hit_up_ok)
      }
    }
    if (n_genes_down > 0L && !is.null(signature_genes_down) && length(signature_genes_down) > 0L) {
      m_dn <- tryCatch(match_genes(mat, signature_genes_down, use_aliases = TRUE)$found, error = function(e) character(0))
      hit_dn_ok <- match(m_dn, rownames(mat))
      hit_dn_ok <- hit_dn_ok[!is.na(hit_dn_ok) & gene_ok[hit_dn_ok]]
      if (length(hit_dn_ok) > 0L && length(hit_dn_ok) < n_genes_down_eff) {
        n_genes_down_eff <- length(hit_dn_ok)
      }
    }
  }

  need_genes <- n_genes_eff + n_genes_down_eff
  pool_idx <- which(gene_ok & !(rows_up %in% exclude_up))
  pool_genes <- rownames(mat)[pool_idx]
  pool_size <- length(pool_genes)

  if (pool_size < need_genes) {
    warning(
      "Null pool size (", pool_size, ") is below n_genes+n_genes_down (",
      need_genes, "); null permutation test skipped.",
      call. = FALSE
    )
    return(list(
      obs_auc = NA_real_,
      null_aucs = numeric(0),
      null_p = NA_real_,
      n_perm = as.integer(n_perm),
      n_genes = as.integer(n_genes),
      n_genes_effective = as.integer(n_genes_eff),
      pool_size = pool_size,
      universe_size = as.integer(universe_size),
      partial_background = partial_background,
      background_type = background_type,
      seed = seed
    ))
  }

  # Gozlenen AUC
  if (is.null(observed_scores)) {
    sig_use <- if (is.null(signature_genes)) exclude_genes else signature_genes
    scored <- score_signature(
      mat,
      genes = if (is.null(signature_genes_down) && is.null(weights)) sig_use else NULL,
      genes_up = if (!is.null(signature_genes_down) && is.null(weights)) sig_use else NULL,
      genes_down = if (!is.null(signature_genes_down) && is.null(weights)) signature_genes_down else NULL,
      weights = weights,
      method = method,
      min_genes = min_genes
    )
    if (isTRUE(scored$skipped) || is.null(scored$scores)) {
      stop("Observed signature could not be scored: ", scored$skip_reason, call. = FALSE)
    }
    observed_scores <- scored$scores
  }

  # Agirlikli imzalar icin up ve down agirlik buyukluklerini cozumle
  w_up_vec <- NULL
  w_dn_vec <- NULL
  if (!is.null(weights) && is.numeric(weights) && length(weights) > 0L && !is.null(names(weights))) {
    w_clean <- stats::setNames(as.numeric(weights), toupper(trimws(names(weights))))
    w_clean <- w_clean[nzchar(names(w_clean)) & is.finite(w_clean) & w_clean != 0]
    if (length(w_clean) > 0L) {
      pos_q <- names(w_clean)[w_clean > 0]
      neg_q <- names(w_clean)[w_clean < 0]
      if (length(pos_q) > 0L && n_genes_eff > 0L) {
        cov_p <- tryCatch(match_genes(mat, pos_q, use_aliases = TRUE), error = function(e) NULL)
        if (!is.null(cov_p) && length(cov_p$found) > 0L) {
          w_p <- abs(as.numeric(.resolve_found_weights(cov_p$found, w_clean, cov_p$alias_rescued)))
          w_p <- w_p[is.finite(w_p) & w_p > 0]
          if (length(w_p) > 0L) {
            w_up_vec <- rep_len(w_p, n_genes_eff)
          }
        }
      }
      if (length(neg_q) > 0L && n_genes_down_eff > 0L) {
        cov_n <- tryCatch(match_genes(mat, neg_q, use_aliases = TRUE), error = function(e) NULL)
        if (!is.null(cov_n) && length(cov_n$found) > 0L) {
          w_n <- abs(as.numeric(.resolve_found_weights(cov_n$found, w_clean, cov_n$alias_rescued)))
          w_n <- w_n[is.finite(w_n) & w_n > 0]
          if (length(w_n) > 0L) {
            w_dn_vec <- rep_len(w_n, n_genes_down_eff)
          }
        }
      }
    }
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

  # Faz 4.6: 2D (5 Ortalama Ifade Quantile x 3 Gen SD/Dispersiyon Tertile) eslestirme
  pool_means <- row_mean[pool_idx]
  pool_sds <- row_sd[pool_idx]
  all_means <- row_mean
  all_sds <- row_sd
  brks_mean <- unique(stats::quantile(
    all_means[is.finite(all_means)],
    probs = seq(0, 1, length.out = 6),
    na.rm = TRUE, names = FALSE
  ))
  brks_sd <- unique(stats::quantile(
    all_sds[is.finite(all_sds) & all_sds > 0],
    probs = seq(0, 1, length.out = 4),
    na.rm = TRUE, names = FALSE
  ))

  pool_bins_mean <- NULL
  pool_bins_2d <- NULL
  sig_bins_mean <- NULL
  sig_bins_2d <- NULL
  sig_bins_down_mean <- NULL
  sig_bins_down_2d <- NULL
  has_sd_bins <- length(brks_sd) >= 2L

  match_up <- isTRUE(match_expression) && n_genes > 0L
  match_down <- isTRUE(match_expression) &&
    n_genes_down > 0L &&
    !is.null(signature_genes_down) &&
    length(signature_genes_down) > 0L

  if (length(brks_mean) < 2L) {
    match_up <- FALSE
    match_down <- FALSE
  } else {
    pool_bins_mean <- as.integer(cut(pool_means, breaks = brks_mean, include.lowest = TRUE, labels = FALSE))
    if (has_sd_bins) {
      pool_bins_sd <- as.integer(cut(pool_sds, breaks = brks_sd, include.lowest = TRUE, labels = FALSE))
      pool_bins_sd[is.na(pool_bins_sd)] <- 1L
      pool_bins_2d <- (pool_bins_mean - 1L) * 3L + pool_bins_sd
    } else {
      pool_bins_2d <- pool_bins_mean
    }
  }

  if (isTRUE(match_up)) {
    sig_ref <- if (is.null(signature_genes)) exclude_genes else signature_genes
    found_ref_up <- if (length(sig_ref) > 0L) {
      tryCatch(match_genes(mat, sig_ref, use_aliases = TRUE)$found, error = function(e) character(0))
    } else {
      character(0)
    }
    hit <- match(found_ref_up, rownames(mat))
    hit <- hit[!is.na(hit)]
    if (length(hit) == 0L) {
      warning("match_expression: no up signature genes found; falling back to random sampling.", call. = FALSE)
      match_up <- FALSE
    } else {
      sig_means <- row_mean[hit]
      sig_sds <- row_sd[hit]
      keep_s <- is.finite(sig_means)
      sig_means <- sig_means[keep_s]
      sig_sds <- sig_sds[keep_s]
      if (length(sig_means) == 0L) {
        warning("match_expression: up mean expression is NA; falling back to random sampling.", call. = FALSE)
        match_up <- FALSE
      } else {
        m_bins <- as.integer(cut(sig_means, breaks = brks_mean, include.lowest = TRUE, labels = FALSE))
        ok_b <- !is.na(m_bins)
        m_bins <- m_bins[ok_b]
        if (length(m_bins) == 0L) {
          warning("match_expression: sig_bins NA; falling back to random sampling.", call. = FALSE)
          match_up <- FALSE
        } else {
          sig_bins_mean <- m_bins
          if (has_sd_bins) {
            s_bins <- as.integer(cut(sig_sds[ok_b], breaks = brks_sd, include.lowest = TRUE, labels = FALSE))
            s_bins[is.na(s_bins)] <- 1L
            sig_bins_2d <- (m_bins - 1L) * 3L + s_bins
          } else {
            sig_bins_2d <- m_bins
          }
        }
      }
    }
  }

  if (isTRUE(match_down)) {
    found_ref_dn <- if (length(signature_genes_down) > 0L) {
      tryCatch(match_genes(mat, signature_genes_down, use_aliases = TRUE)$found, error = function(e) character(0))
    } else {
      character(0)
    }
    hit_dn <- match(found_ref_dn, rownames(mat))
    hit_dn <- hit_dn[!is.na(hit_dn)]
    if (length(hit_dn) == 0L) {
      warning("match_expression: no down signature genes found; falling back to random sampling.", call. = FALSE)
      match_down <- FALSE
    } else {
      sig_means_dn <- row_mean[hit_dn]
      sig_sds_dn <- row_sd[hit_dn]
      keep_sd <- is.finite(sig_means_dn)
      sig_means_dn <- sig_means_dn[keep_sd]
      sig_sds_dn <- sig_sds_dn[keep_sd]
      if (length(sig_means_dn) == 0L) {
        warning("match_expression: down mean expression is NA; falling back to random sampling.", call. = FALSE)
        match_down <- FALSE
      } else {
        m_bins_dn <- as.integer(cut(sig_means_dn, breaks = brks_mean, include.lowest = TRUE, labels = FALSE))
        ok_bd <- !is.na(m_bins_dn)
        m_bins_dn <- m_bins_dn[ok_bd]
        if (length(m_bins_dn) == 0L) {
          warning("match_expression: sig_bins_down NA; falling back to random sampling.", call. = FALSE)
          match_down <- FALSE
        } else {
          sig_bins_down_mean <- m_bins_dn
          if (has_sd_bins) {
            s_bins_dn <- as.integer(cut(sig_sds_dn[ok_bd], breaks = brks_sd, include.lowest = TRUE, labels = FALSE))
            s_bins_dn[is.na(s_bins_dn)] <- 1L
            sig_bins_down_2d <- (m_bins_dn - 1L) * 3L + s_bins_dn
          } else {
            sig_bins_down_2d <- m_bins_dn
          }
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

  .sample_matched_arm <- function(n_need, bins_2d, bins_mean, already_picked) {
    if (n_need < 1L) {
      return(character(0))
    }
    avail <- setdiff(pool_genes, already_picked)
    if (length(avail) < n_need) {
      return(NULL)
    }
    picked <- character(0)
    for (idx_b in seq_along(bins_2d)) {
      if (length(picked) >= n_need) {
        break
      }
      bn2 <- bins_2d[idx_b]
      bn1 <- if (!is.null(bins_mean) && length(bins_mean) >= idx_b) bins_mean[idx_b] else NA_integer_
      # Oncelikle 2D (ortalama + SD) hucresinden aday sec
      cand_idx <- which(
        !is.na(pool_bins_2d) &
          pool_bins_2d == bn2 &
          pool_genes %in% avail &
          !(pool_genes %in% picked)
      )
      # 2D hucrede aday yoksa ayni ortalama ifade quantile (1D) hucresine geri don
      if (length(cand_idx) == 0L && !is.na(bn1)) {
        cand_idx <- which(
          !is.na(pool_bins_mean) &
            pool_bins_mean == bn1 &
            pool_genes %in% avail &
            !(pool_genes %in% picked)
        )
      }
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
    if (n_genes_eff > 0L) {
      if (isTRUE(match_up) && !is.null(sig_bins_2d)) {
        picked_up <- .sample_matched_arm(n_genes_eff, sig_bins_2d, sig_bins_mean, character(0))
        if (is.null(picked_up)) {
          return(NULL)
        }
      } else {
        picked_up <- sample(pool_genes, size = n_genes_eff, replace = FALSE)
      }
    } else {
      picked_up <- character(0)
    }

    if (n_genes_down_eff > 0L) {
      if (isTRUE(match_down) && !is.null(sig_bins_down_2d)) {
        picked_down <- .sample_matched_arm(n_genes_down_eff, sig_bins_down_2d, sig_bins_down_mean, picked_up)
        if (is.null(picked_down)) {
          return(NULL)
        }
      } else {
        rem <- setdiff(pool_genes, picked_up)
        if (length(rem) < n_genes_down_eff) {
          return(NULL)
        }
        picked_down <- sample(rem, size = n_genes_down_eff, replace = FALSE)
      }
    } else {
      picked_down <- character(0)
    }
    list(up = picked_up, down = picked_down)
  }

  .arm_scores <- function(gene_names, arm_weights = NULL) {
    if (length(gene_names) == 0L) {
      return(NULL)
    }
    row_idx <- match(gene_names, rownames(score_mat))
    if (anyNA(row_idx)) {
      return(NULL)
    }
    sub <- score_mat[row_idx, , drop = FALSE]
    if (!is.null(arm_weights) && length(arm_weights) == nrow(sub) && method != "median") {
      w_mat <- matrix(arm_weights, nrow = nrow(sub), ncol = ncol(sub))
      w_mat[!is.finite(sub)] <- NA_real_
      out <- colSums(sub * w_mat, na.rm = TRUE) / colSums(w_mat, na.rm = TRUE)
      out[colSums(is.finite(sub)) == 0L] <- NA_real_
      return(out)
    }
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
        .ss_rank_arm(rank_mat, picked$up, N_j, reverse = FALSE, arm_weights = w_up_vec)
      } else {
        rep(0, ncol(mat))
      }
      s_down <- if (has_down) {
        .ss_rank_arm(rank_mat, picked$down, N_j, reverse = TRUE, arm_weights = w_dn_vec)
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

    scores_up <- .arm_scores(picked$up, arm_weights = w_up_vec)
    scores_dn <- .arm_scores(picked$down, arm_weights = w_dn_vec)
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
    n_genes_effective = as.integer(n_genes_eff),
    pool_size = pool_size,
    universe_size = as.integer(universe_size),
    partial_background = partial_background,
    background_type = background_type,
    seed = seed
  )
}
