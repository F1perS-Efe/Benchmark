#' Imza skorunu ornek basina hesapla
#'
#' Skorlama yanit etiketini (`response`) asla kullanmaz; tum ornekler
#' uzerinde etiketten bagimsiz hesaplanir.
#'
#' - `mean`: log olcekli ifadenin ornek basina ortalamasi
#' - `zscore`: her gen veri seti (kohort) genelinde z-score'a cevrilir, sonra
#'   ortalanir. Basit ortalama gen z-skorudur; `sqrt(k)` olcekli GSVA degildir.
#' - `median`: ham (log olcekli) ifadenin ornek basina medyani
#' - `rank`: her gen ornekler arasi (kohort-ici) average-rank'e cevrilir, 0-1
#'   araligina olceklenir, sonra ornek basina ortalanir. Bu `singscore` degildir;
#'   kohortlar arasi ham skor birlestirilmez.
#' - `ss_rank`: her ornekte (sutunda) tum non-NA genler dikey siralanir
#'   (`rank(..., ties.method = "average", na.last = "keep")`; 1 = en dusuk ifade).
#'   Kohorta ornek eklenip cikarilsa bile tek bir hastanin skoru degismez.
#'   Tek yonlu skor `[-0.5, +0.5]`, cift yonlu (`genes_up` + `genes_down`) skor
#'   `[-1, +1]` araligindadir; formuller `singscore` tek/cift yonlu
#'   normalizasyonu ile uyumludur (saf base R).
#'
#' Isaretli (`+`/`-`) imzalar:
#' - `genes_up` / `genes_down` parametreleri, veya
#' - isimli sayisal `genes` (`c(IFNG = 1, VEGFA = -1)`), veya
#' - geriye donuk tek yonlu karakter `genes`.
#'
#' `mean` / `median` / `zscore` / `rank` icin `genes_down` doluysa
#' `score = score_up - score_down`. `ss_rank` icin `score = s_up + s_down`.
#'
#' @param expr Gen x ornek matris (hazirlanmis / log olcekli).
#' @param genes Imza genleri (karakter) veya isimli sayisal isaret vektoru.
#'   `genes_up` / `genes_down` verildiyse `NULL` olmalidir.
#' @param method `"mean"`, `"zscore"`, `"median"`, `"rank"` veya `"ss_rank"`.
#' @param min_genes Kapsama sonrasi minimum gen sayisi. Varsayilan 3.
#' @param genes_up Yukari yonlu (pozitif) imza genleri.
#' @param genes_down Asagi yonlu (negatif) imza genleri.
#'
#' @return Liste: `scores`, `coverage`, `n_all_na_dropped`, `skipped`,
#'   `skip_reason`, `signed`, `n_up_found`, `n_down_found`.
#' @export
score_signature <- function(expr,
                            genes = NULL,
                            method = c("mean", "zscore", "median", "rank", "ss_rank"),
                            min_genes = 3L,
                            genes_up = NULL,
                            genes_down = NULL) {
  method <- match.arg(method)
  parsed <- .parse_signature_genes(genes, genes_up, genes_down)

  cov_up <- if (length(parsed$genes_up) > 0L) {
    match_genes(expr, parsed$genes_up)
  } else {
    .empty_gene_coverage()
  }
  cov_down <- if (length(parsed$genes_down) > 0L) {
    match_genes(expr, parsed$genes_down)
  } else {
    .empty_gene_coverage()
  }

  n_up_found <- cov_up$n_found
  n_down_found <- cov_down$n_found
  alias_rescued <- c(cov_up$alias_rescued, cov_down$alias_rescued)
  n_total <- cov_up$n_total + cov_down$n_total
  miss_all <- c(cov_up$missing, cov_down$missing)
  n_matched <- n_total - length(miss_all)
  coverage <- list(
    found = c(cov_up$found, cov_down$found),
    missing = miss_all,
    n_found = n_up_found + n_down_found,
    n_total = n_total,
    coverage_pct = if (n_total == 0L) 0 else 100 * n_matched / n_total,
    alias_rescued = alias_rescued,
    n_alias_rescued = as.integer(length(alias_rescued))
  )

  if (isTRUE(parsed$signed) && (n_up_found == 0L || n_down_found == 0L)) {
    warning(
      "Signed signature: one arm found 0 genes (n_up_found=", n_up_found,
      ", n_down_found=", n_down_found, ").",
      call. = FALSE
    )
  }

  if (coverage$n_found < min_genes) {
    return(list(
      scores = NULL,
      coverage = coverage,
      n_all_na_dropped = 0L,
      skipped = TRUE,
      skip_reason = sprintf(
        "Bulunan gen sayisi (%d) min_genes (%d) altinda.",
        coverage$n_found, min_genes
      ),
      signed = parsed$signed,
      n_up_found = as.integer(n_up_found),
      n_down_found = as.integer(n_down_found)
    ))
  }

  if (method == "ss_rank") {
    scores <- .ss_rank_scores(expr, cov_up$found, cov_down$found)
  } else {
    scores_up <- .score_method_on_genes(expr, cov_up$found, method)
    scores_down <- .score_method_on_genes(expr, cov_down$found, method)
    if (!is.null(scores_up) && !is.null(scores_down)) {
      scores <- scores_up - scores_down
    } else if (!is.null(scores_up)) {
      scores <- scores_up
    } else {
      scores <- -scores_down
    }
  }

  all_na <- is.na(scores)
  n_all_na_dropped <- sum(all_na)
  if (n_all_na_dropped > 0L) {
    scores <- scores[!all_na]
  }

  list(
    scores = scores,
    coverage = coverage,
    n_all_na_dropped = n_all_na_dropped,
    skipped = FALSE,
    skip_reason = NA_character_,
    signed = parsed$signed,
    n_up_found = as.integer(n_up_found),
    n_down_found = as.integer(n_down_found)
  )
}

#' Parse unidirectional or signed signature gene arguments
#'
#' @noRd
.parse_signature_genes <- function(genes = NULL,
                                   genes_up = NULL,
                                   genes_down = NULL) {
  has_updown <- !is.null(genes_up) || !is.null(genes_down)
  if (has_updown && !is.null(genes)) {
    stop(
      "Provide either `genes` or `genes_up`/`genes_down`, not both.",
      call. = FALSE
    )
  }

  if (has_updown) {
    gu <- if (is.null(genes_up)) character(0) else as.character(genes_up)
    gd <- if (is.null(genes_down)) character(0) else as.character(genes_down)
    gu <- gu[!is.na(gu) & nzchar(gu)]
    gd <- gd[!is.na(gd) & nzchar(gd)]
  } else if (is.null(genes)) {
    stop("`genes` or `genes_up`/`genes_down` is required.", call. = FALSE)
  } else if (is.numeric(genes)) {
    nm <- names(genes)
    if (is.null(nm) || any(is.na(nm)) || any(!nzchar(nm))) {
      stop(
        "Numeric `genes` must be a named sign vector (e.g. c(IFNG = 1, VEGFA = -1)).",
        call. = FALSE
      )
    }
    if (anyNA(genes) || any(genes == 0)) {
      stop(
        "Numeric `genes` signs must be non-zero and non-NA.",
        call. = FALSE
      )
    }
    gu <- nm[genes > 0]
    gd <- nm[genes < 0]
  } else {
    gu <- as.character(genes)
    gu <- gu[!is.na(gu) & nzchar(gu)]
    gd <- character(0)
  }

  gu <- unique(gu)
  gd <- unique(gd)
  if (length(gu) == 0L && length(gd) == 0L) {
    stop("Signature gene list is empty.", call. = FALSE)
  }

  overlap <- intersect(toupper(gu), toupper(gd))
  if (length(overlap) > 0L) {
    stop(
      "genes_up and genes_down overlap: ",
      paste(overlap, collapse = ", "),
      call. = FALSE
    )
  }

  list(
    genes_up = gu,
    genes_down = gd,
    signed = length(gd) > 0L
  )
}

#' @noRd
.empty_gene_coverage <- function() {
  list(
    found = character(0),
    missing = character(0),
    n_found = 0L,
    n_total = 0L,
    coverage_pct = 0,
    alias_rescued = character(0),
    n_alias_rescued = 0L
  )
}

#' Cohort-relative mean/median/zscore/rank on a gene subset
#'
#' @noRd
.score_method_on_genes <- function(expr, found, method) {
  if (length(found) == 0L) {
    return(NULL)
  }
  sub <- expr[found, , drop = FALSE]

  if (method == "zscore") {
    sub <- t(scale(t(sub)))
  }

  if (method == "rank") {
    sub <- t(apply(sub, 1L, function(x) {
      m <- sum(!is.na(x))
      if (m <= 1L) {
        return(rep(NA_real_, length(x)))
      }
      r <- rank(x, ties.method = "average", na.last = "keep")
      (r - 1) / (m - 1)
    }))
    if (!is.null(colnames(expr))) {
      colnames(sub) <- colnames(expr)[seq_len(ncol(sub))]
    }
    rownames(sub) <- found
  }

  apply(sub, 2L, function(x) {
    if (all(is.na(x))) {
      return(NA_real_)
    }
    if (method == "median") {
      stats::median(x, na.rm = TRUE)
    } else {
      mean(x, na.rm = TRUE)
    }
  })
}

#' Single-sample vertical rank scores (singscore-compatible normalization)
#'
#' @noRd
.ss_rank_scores <- function(expr, found_up, found_down) {
  nc <- ncol(expr)
  rank_mat <- apply(expr, 2L, function(x) {
    rank(x, ties.method = "average", na.last = "keep")
  })
  if (is.null(dim(rank_mat))) {
    rank_mat <- matrix(rank_mat, nrow = nrow(expr), ncol = nc)
  }
  rownames(rank_mat) <- rownames(expr)
  colnames(rank_mat) <- colnames(expr)
  N <- colSums(!is.na(expr))

  has_up <- length(found_up) > 0L
  has_down <- length(found_down) > 0L

  s_up <- if (has_up) {
    .ss_rank_arm(rank_mat, found_up, N, reverse = FALSE)
  } else {
    rep(0, nc)
  }
  s_down <- if (has_down) {
    .ss_rank_arm(rank_mat, found_down, N, reverse = TRUE)
  } else {
    rep(0, nc)
  }

  scores <- if (has_up && has_down) {
    s_up + s_down
  } else if (has_up) {
    s_up
  } else {
    s_down
  }
  names(scores) <- colnames(expr)
  scores
}

#' @noRd
.ss_rank_arm <- function(rank_mat, genes, N, reverse = FALSE) {
  sub <- rank_mat[genes, , drop = FALSE]
  if (isTRUE(reverse)) {
    sub <- sweep(sub, 2L, N + 1, function(r, np1) np1 - r)
  }
  m <- colMeans(sub, na.rm = TRUE)
  k <- colSums(!is.na(sub))
  out <- (m - (k + 1) / 2) / (N - k) - 0.5
  bad <- !is.finite(out) | N <= k | k == 0L
  out[bad] <- NA_real_
  out
}
