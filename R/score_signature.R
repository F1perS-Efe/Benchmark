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
#' @param weights Optional named numeric vector of continuous gene weights
#'   (e.g. published regression coefficients such as Ayers/Cristescu 18-gene
#'   T cell-inflamed GEP where `CD276 = -0.023900` and `STAT1 = 0.250229`).
#'   When provided, scores are computed as the weighted linear combination of
#'   matched gene values (`method` transformation applied per gene first).
#'
#' @return Liste: `scores`, `coverage`, `n_all_na_dropped`, `skipped`,
#'   `skip_reason`, `signed`, `weighted`, `n_up_found`, `n_down_found`.
#' @examples
#' mat <- matrix(
#'   c(2.1, 3.5, 1.2, 4.0, 2.8, 3.9, 1.5, 4.2, 0.5, 0.2, 1.8, 0.3),
#'   nrow = 3, byrow = TRUE,
#'   dimnames = list(c("STAT1", "CXCL9", "CD276"), c("S1", "S2", "S3", "S4"))
#' )
#' # Unidirectional or continuously weighted signature scoring
#' s_unidir <- score_signature(mat, genes = c("STAT1", "CXCL9"), method = "zscore", min_genes = 2L)
#' s_weighted <- score_signature(
#'   mat,
#'   weights = c(STAT1 = 0.250229, CXCL9 = 0.074135, CD276 = -0.023900),
#'   method = "zscore",
#'   min_genes = 2L
#' )
#' @export
score_signature <- function(expr,
                            genes = NULL,
                            method = c("mean", "zscore", "median", "rank", "ss_rank"),
                            min_genes = 3L,
                            genes_up = NULL,
                            genes_down = NULL,
                            weights = NULL) {
  method <- match.arg(method)
  parsed <- .parse_signature_genes(genes, genes_up, genes_down, weights = weights)
  if (isTRUE(parsed$weighted) && identical(method, "median")) {
    stop(
      "Weighted continuous signatures do not support method = 'median'; use 'mean', 'zscore', 'rank', or 'ss_rank'.",
      call. = FALSE
    )
  }

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

  if (isTRUE(parsed$signed) && !isTRUE(parsed$weighted) &&
      (n_up_found == 0L || n_down_found == 0L)) {
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
        "Found gene count (%d) is below min_genes (%d).",
        coverage$n_found, min_genes
      ),
      signed = parsed$signed,
      weighted = isTRUE(parsed$weighted),
      n_up_found = as.integer(n_up_found),
      n_down_found = as.integer(n_down_found)
    ))
  }

  if (isTRUE(parsed$weighted)) {
    found_all <- c(cov_up$found, cov_down$found)
    w_found <- .resolve_found_weights(found_all, parsed$weights, alias_rescued)
    if (method == "ss_rank") {
      ss_res <- .resolve_ss_rank_with_context(
        expr = expr,
        found_up = cov_up$found,
        found_down = cov_down$found,
        weights = w_found
      )
      scores <- ss_res$scores
      coverage$ss_rank_universe_size <- ss_res$universe_size
      coverage$partial_ortholog_universe <- ss_res$partial_universe
    } else {
      scores <- .score_weighted_on_genes(expr, found_all, w_found, method)
    }
  } else if (method == "ss_rank") {
    ss_res <- .resolve_ss_rank_with_context(
      expr = expr,
      found_up = cov_up$found,
      found_down = cov_down$found,
      weights = NULL
    )
    scores <- ss_res$scores
    coverage$ss_rank_universe_size <- ss_res$universe_size
    coverage$partial_ortholog_universe <- ss_res$partial_universe
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
    weighted = isTRUE(parsed$weighted),
    n_up_found = as.integer(n_up_found),
    n_down_found = as.integer(n_down_found)
  )
}

#' Parse unidirectional, signed, or continuously weighted signature arguments
#'
#' @noRd
.parse_signature_genes <- function(genes = NULL,
                                   genes_up = NULL,
                                   genes_down = NULL,
                                   weights = NULL) {
  has_updown <- !is.null(genes_up) || !is.null(genes_down)
  if (is.null(weights) && has_updown && !is.null(genes)) {
    stop(
      "Provide either `genes` or `genes_up`/`genes_down`, not both.",
      call. = FALSE
    )
  }

  # Continuous weights via explicit `weights` or non-unit numeric `genes`
  if (is.null(weights) && is.numeric(genes)) {
    nm_g <- names(genes)
    if (is.null(nm_g) || any(is.na(nm_g)) || any(!nzchar(nm_g))) {
      stop(
        "Numeric `genes` must be a named sign/weight vector (e.g. c(IFNG = 1, VEGFA = -1)).",
        call. = FALSE
      )
    }
    if (any(!is.finite(genes)) || any(genes == 0)) {
      stop(
        "Numeric `genes` signs/weights must be finite and non-zero.",
        call. = FALSE
      )
    }
    if (!all(abs(as.numeric(genes)) == 1)) {
      weights <- genes
    }
  }

  if (!is.null(weights)) {
    if (!is.numeric(weights) || is.null(names(weights))) {
      stop("`weights` must be a named numeric vector of gene coefficients.", call. = FALSE)
    }
    nm <- trimws(as.character(names(weights)))
    if (any(is.na(nm)) || any(!nzchar(nm))) {
      stop("`weights` must have non-empty gene symbol names.", call. = FALSE)
    }
    w_num <- as.numeric(weights)
    if (any(!is.finite(w_num)) || any(w_num == 0)) {
      stop("`weights` values must be finite and non-zero.", call. = FALSE)
    }
    if (anyDuplicated(toupper(nm))) {
      stop("`weights` contains duplicate gene names.", call. = FALSE)
    }
    gu <- nm[w_num > 0]
    gd <- nm[w_num < 0]
    w_map <- stats::setNames(w_num, toupper(nm))
    return(list(
      genes_up = gu,
      genes_down = gd,
      signed = length(gd) > 0L,
      weighted = TRUE,
      weights = w_map
    ))
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
    signed = length(gd) > 0L,
    weighted = FALSE,
    weights = NULL
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
    rn <- rowSums(is.finite(sub))
    rm <- rowMeans(sub, na.rm = TRUE)
    rsd <- sqrt(rowSums((sub - rm)^2, na.rm = TRUE) / pmax(rn - 1L, 1L))
    sub <- (sub - rm) / rsd
    sub[!is.finite(sub)] <- NA_real_
  } else if (method == "rank") {
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

  if (method == "median") {
    apply(sub, 2L, function(x) {
      if (all(is.na(x))) NA_real_ else stats::median(x, na.rm = TRUE)
    })
  } else {
    out <- colMeans(sub, na.rm = TRUE)
    out[colSums(!is.na(sub)) == 0L] <- NA_real_
    out
  }
}

#' Map matched expression matrix row names to their signed gene weights
#'
#' Handles direct case-insensitive matches and alias-rescued symbols
#' (`"QUERY->OFFICIAL"`).
#'
#' @noRd
.resolve_found_weights <- function(found_rows, weights_upper, alias_rescued = character(0)) {
  if (length(found_rows) == 0L) {
    return(numeric(0))
  }
  alias_rev <- character(0)
  if (length(alias_rescued) > 0L) {
    for (ar in alias_rescued) {
      parts <- strsplit(as.character(ar), "->", fixed = TRUE)[[1]]
      if (length(parts) == 2L) {
        alias_rev[[toupper(parts[2])]] <- toupper(parts[1])
      }
    }
  }
  w_out <- numeric(length(found_rows))
  names(w_out) <- found_rows
  for (i in seq_along(found_rows)) {
    key <- toupper(found_rows[i])
    if (key %in% names(weights_upper)) {
      w_out[i] <- weights_upper[[key]]
    } else if (key %in% names(alias_rev) && alias_rev[[key]] %in% names(weights_upper)) {
      w_out[i] <- weights_upper[[alias_rev[[key]]]]
    } else {
      w_out[i] <- 0
    }
  }
  w_out
}

#' Compute weighted linear signature scores across samples
#'
#' Applies the per-gene `method` transformation (`"mean"` uses the
#' prepared log-scale matrix directly; `"zscore"` standardizes each gene across
#' the cohort; `"rank"` scales per-gene ranks to `[0, 1]`), then computes the
#' weighted linear predictor `sum(w_i * x_{i,j})` (re-scaled proportionally if
#' an individual sample has `NA` values in a subset of matched genes).
#'
#' @noRd
.score_weighted_on_genes <- function(expr, found, weights_found, method) {
  if (identical(method, "median")) {
    stop(
      "Weighted continuous signatures do not support method = 'median'; use 'mean', 'zscore', 'rank', or 'ss_rank'.",
      call. = FALSE
    )
  }
  if (length(found) == 0L) {
    return(NULL)
  }
  sub <- expr[found, , drop = FALSE]
  w <- as.numeric(weights_found[found])

  if (method == "zscore") {
    rn <- rowSums(is.finite(sub))
    rm <- rowMeans(sub, na.rm = TRUE)
    rsd <- sqrt(rowSums((sub - rm)^2, na.rm = TRUE) / pmax(rn - 1L, 1L))
    sub <- (sub - rm) / rsd
    sub[!is.finite(sub)] <- NA_real_
  } else if (method == "rank") {
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

  total_abs_w <- sum(abs(w))
  out <- apply(sub, 2L, function(x) {
    ok <- is.finite(x) & is.finite(w) & (w != 0)
    if (!any(ok)) {
      return(NA_real_)
    }
    w_ok <- w[ok]
    x_ok <- x[ok]
    denom <- sum(abs(w_ok))
    scale_fac <- if (denom > 0 && total_abs_w > 0) total_abs_w / denom else 1
    sum(w_ok * x_ok) * scale_fac
  })
  names(out) <- colnames(expr)
  out
}

#' Evaluate ss_rank using full transcriptome context when present and warn on partial ortholog universes
#'
#' @noRd
.resolve_ss_rank_with_context <- function(expr, found_up, found_down, weights = NULL) {
  ostat <- attr(expr, "ortholog_stats")
  bg_expr <- attr(expr, "transcriptome_bg_expr")
  rank_expr <- if (is.matrix(bg_expr) && ncol(bg_expr) == ncol(expr) &&
                   all(c(found_up, found_down) %in% rownames(bg_expr))) {
    bg_expr
  } else {
    expr
  }
  universe_size <- nrow(rank_expr)
  partial_universe <- FALSE
  if (!is.null(ostat) && isTRUE(ostat$is_partial_mapping) && universe_size <= ostat$n_mapped) {
    partial_universe <- TRUE
    warning(
      sprintf(
        paste0(
          "ss_rank on partial ortholog-mapped matrix (mapped=%d of %d input genes, ",
          "retention=%.1f%%) without preserved background context: vertical ranking ",
          "universe is restricted to %d mapped genes and does not represent full ",
          "transcriptome background."
        ),
        ostat$n_mapped, ostat$n_input, 100 * ostat$retention_rate, universe_size
      ),
      call. = FALSE
    )
  }
  scores <- .ss_rank_scores(rank_expr, found_up, found_down, weights = weights)
  list(
    scores = scores,
    universe_size = as.integer(universe_size),
    partial_universe = partial_universe
  )
}

#' Single-sample vertical rank scores (singscore-compatible normalization)
#'
#' @noRd
.ss_rank_scores <- function(expr, found_up, found_down, weights = NULL) {
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

  w_up <- if (!is.null(weights) && has_up) abs(as.numeric(weights[found_up])) else NULL
  w_down <- if (!is.null(weights) && has_down) abs(as.numeric(weights[found_down])) else NULL

  s_up <- if (has_up) {
    .ss_rank_arm(rank_mat, found_up, N, reverse = FALSE, arm_weights = w_up)
  } else {
    rep(0, nc)
  }
  s_down <- if (has_down) {
    .ss_rank_arm(rank_mat, found_down, N, reverse = TRUE, arm_weights = w_down)
  } else {
    rep(0, nc)
  }

  scores <- if (has_up && has_down) {
    if (!is.null(weights)) {
      sum_up <- sum(w_up, na.rm = TRUE)
      sum_dn <- sum(w_down, na.rm = TRUE)
      tot <- sum_up + sum_dn
      if (tot > 0) {
        (sum_up * s_up + sum_dn * s_down) / tot
      } else {
        s_up + s_down
      }
    } else {
      s_up + s_down
    }
  } else if (has_up) {
    s_up
  } else {
    s_down
  }
  names(scores) <- colnames(expr)
  scores
}

#' @noRd
.ss_rank_arm <- function(rank_mat, genes, N, reverse = FALSE, arm_weights = NULL) {
  sub <- rank_mat[genes, , drop = FALSE]
  if (isTRUE(reverse)) {
    sub <- sweep(sub, 2L, N + 1, function(r, np1) np1 - r)
  }
  if (is.null(arm_weights) || length(arm_weights) != nrow(sub)) {
    m <- colMeans(sub, na.rm = TRUE)
  } else {
    m <- apply(sub, 2L, function(col_r) {
      ok <- !is.na(col_r) & is.finite(arm_weights) & (arm_weights > 0)
      if (!any(ok)) return(NA_real_)
      sum(col_r[ok] * arm_weights[ok]) / sum(arm_weights[ok])
    })
  }
  k <- colSums(!is.na(sub))
  out <- (m - (k + 1) / 2) / (N - k) - 0.5
  bad <- !is.finite(out) | N <= k | k == 0L
  out[bad] <- NA_real_
  out
}

#' Normalize a signature definition into a standard list contract
#'
#' @noRd
.as_signature <- function(sig, name = NULL) {
  label <- if (is.null(name) || !nzchar(as.character(name)[1])) {
    "Signature"
  } else {
    paste0("Signature '", name, "'")
  }
  if (is.character(sig) || is.numeric(sig)) {
    if (length(sig) == 0L) {
      stop(label, " cannot be empty.", call. = FALSE)
    }
    if (is.numeric(sig) && !is.null(names(sig)) && !all(abs(as.numeric(sig)) == 1)) {
      return(list(
        genes = names(sig),
        genes_up = names(sig)[sig > 0],
        genes_down = names(sig)[sig < 0],
        weights = sig,
        direction = "signed"
      ))
    }
    return(list(genes = sig))
  }
  if (!is.list(sig)) {
    stop(label, " must be a character/numeric vector or a list.", call. = FALSE)
  }
  if (!is.null(sig$weights) && is.numeric(sig$weights) && !is.null(names(sig$weights))) {
    if (is.null(sig$genes) && is.null(sig$genes_up) && is.null(sig$genes_down)) {
      sig$genes <- names(sig$weights)
    }
  }
  if (is.null(sig$genes) && is.null(sig$genes_up) && is.null(sig$genes_down) && is.null(sig$weights)) {
    stop(
      label, " requires `genes`, `genes_up`/`genes_down`, or `weights`.",
      call. = FALSE
    )
  }
  sig
}

#' Resolve effective minimum gene threshold from signature metadata
#'
#' If `sig` is a list carrying a valid positive integer/numeric `min_genes`
#' field (such as 2-gene `Cytolytic_Activity`), returns `as.integer(sig$min_genes[1])`;
#' otherwise returns `as.integer(default[1])`.
#'
#' @param sig Signature object (list or vector) or `NULL`.
#' @param default Fallback minimum gene count (default `3L`).
#' @return Integer minimum gene count (`>= 1L`).
#' @noRd
.resolve_min_genes <- function(sig = NULL, default = 3L) {
  def <- suppressWarnings(as.integer(default)[1])
  if (is.na(def) || def < 1L) {
    def <- 3L
  }
  if (is.list(sig) && !is.null(sig$min_genes) && length(sig$min_genes) >= 1L) {
    mg <- suppressWarnings(as.integer(sig$min_genes[1]))
    if (!is.na(mg) && mg >= 1L) {
      return(mg)
    }
  }
  def
}
