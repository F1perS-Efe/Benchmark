## Lossless full-gene pool evaluation helpers (not run by R CMD check).
## Usage (from package root): source("data-raw/evaluate_full_gene_pool.R")
##
## Measures xz-compressed size of a full-gene cohort matrix and, when
## rounding is tested, compares mean/median/zscore/rank AUC on original vs
## rounded matrices. Rejects rounding if max_auc_diff > 1e-3 or rank ties
## increase.

#' Evaluate lossless storage and optional rounding impact on AUC
#'
#' @param expr Numeric gene x sample matrix (full gene pool).
#' @param response Named 0/1 response vector aligned to columns (optional;
#'   if NULL, a balanced synthetic response is created).
#' @param signature_genes Character vector of genes for scoring.
#' @param round_digits Integer vector of decimal places to test (e.g. 2:4).
#'   Empty / NULL skips rounding tests.
#' @param seed Integer seed for synthetic response and reproducibility.
#' @param auc_tol Maximum allowed absolute AUC difference (default 1e-3).
#'
#' @return A list with `xz_bytes`, `rounding` (data.frame or NULL), and
#'   `accepted` (logical).
evaluate_full_gene_pool <- function(expr,
                                    response = NULL,
                                    signature_genes = NULL,
                                    round_digits = c(2L, 3L, 4L),
                                    seed = 42L,
                                    auc_tol = 1e-3) {
  stopifnot(is.matrix(expr), is.numeric(expr))
  if (is.null(colnames(expr))) {
    colnames(expr) <- paste0("S", seq_len(ncol(expr)))
  }
  if (is.null(rownames(expr))) {
    stop("`expr` rownames (genes) required.", call. = FALSE)
  }

  # Lossless xz size
  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(expr, file = tmp, compress = "xz")
  xz_bytes <- file.info(tmp)$size

  if (is.null(signature_genes)) {
    signature_genes <- rownames(expr)[seq_len(min(10L, nrow(expr)))]
  }
  signature_genes <- as.character(signature_genes)

  if (is.null(response)) {
    set.seed(as.integer(seed))
    n <- ncol(expr)
    n_pos <- floor(n / 2)
    response <- c(rep(1L, n_pos), rep(0L, n - n_pos))
    names(response) <- colnames(expr)
  }
  response <- as.integer(response[colnames(expr)])
  names(response) <- colnames(expr)

  methods <- c("mean", "median", "zscore", "rank")

  score_auc <- function(mat) {
    out <- numeric(length(methods))
    names(out) <- methods
    for (m in methods) {
      scored <- score_signature(mat, genes = signature_genes, method = m, min_genes = 1L)
      if (isTRUE(scored$skipped)) {
        out[[m]] <- NA_real_
        next
      }
      keep <- is.finite(scored$scores) & !is.na(response[names(scored$scores)])
      if (sum(keep) < 4L || length(unique(response[names(scored$scores)][keep])) < 2L) {
        out[[m]] <- NA_real_
        next
      }
      roc <- roc_auc(
        scored$scores[keep],
        response[names(scored$scores)][keep],
        direction = "<",
        min_n = 1L,
        min_per_class = 1L
      )
      out[[m]] <- roc$auc
    }
    out
  }

  count_rank_ties <- function(mat) {
    genes_use <- intersect(toupper(signature_genes), toupper(rownames(mat)))
    if (length(genes_use) == 0L) {
      return(0L)
    }
    idx <- match(genes_use, toupper(rownames(mat)))
    sub <- mat[idx, , drop = FALSE]
    sum(apply(sub, 1L, function(x) {
      x <- x[is.finite(x)]
      if (length(x) < 2L) return(0L)
      length(x) - length(unique(x))
    }))
  }

  base_auc <- score_auc(expr)
  base_ties <- count_rank_ties(expr)

  rounding <- NULL
  accepted <- TRUE
  if (!is.null(round_digits) && length(round_digits) > 0L) {
    rows <- list()
    for (d in as.integer(round_digits)) {
      rounded <- round(expr, digits = d)
      auc_r <- score_auc(rounded)
      ties_r <- count_rank_ties(rounded)
      diffs <- abs(auc_r - base_auc)
      max_diff <- max(diffs, na.rm = TRUE)
      if (!is.finite(max_diff)) max_diff <- Inf
      tie_increase <- ties_r > base_ties
      reject <- isTRUE(max_diff > auc_tol) || isTRUE(tie_increase)
      if (reject) accepted <- FALSE
      rows[[length(rows) + 1L]] <- data.frame(
        digits = d,
        max_auc_diff = max_diff,
        base_ties = as.integer(base_ties),
        rounded_ties = as.integer(ties_r),
        tie_increase = tie_increase,
        reject = reject,
        stringsAsFactors = FALSE
      )
    }
    rounding <- do.call(rbind, rows)
    rownames(rounding) <- NULL
  }

  list(
    xz_bytes = as.numeric(xz_bytes),
    base_auc = base_auc,
    rounding = rounding,
    accepted = accepted,
    auc_tol = auc_tol
  )
}
