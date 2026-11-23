#' Correlate user signatures with reference immune axes
#'
#' Scores each user signature and each reference axis on the same SD-filtered
#' samples (`encode_response(...)$kept`), then reports Spearman (or other)
#' correlations. When `ref_signatures = NULL`, three built-in axes are used:
#' `IFN_Gamma` (Ayers IFNG6), `Cytotoxicity` (GZMA/GZMB/PRF1), and
#' `T_Cell_Infiltration` (CD3D/CD3E/CD8A/CD8B).
#'
#' @param datasets Named dataset list (`validate_dataset` structure).
#' @param signatures Named user signature list.
#' @param ref_signatures Optional named reference signature list. `NULL` uses
#'   the built-in axes above.
#' @param method Scoring method for `score_signature()`. Default `"zscore"`.
#' @param cor_method Correlation method for `stats::cor` / `stats::cor.test`.
#'   Default `"spearman"`.
#' @param sd_handling Passed to `encode_response()`. Default `"exclude"`.
#' @param apply_log2 Passed to `prepare_expression()`. Default `TRUE`.
#' @param min_genes Minimum genes for scoring. Default `3L`.
#' @param overlap_threshold Jaccard or overlap fraction threshold at or above
#'   which `overlap_flag` is set to `TRUE` (default `0.50`).
#'
#' @return `data.frame` with columns `cohort`, `signature`, `reference_axis`,
#'   `cor`, `p_value`, `n_samples`, `normalization`, `jaccard_overlap`,
#'   `overlap_flag`, `overlap_note`.
#' @export
correlate_reference_axes <- function(datasets,
                                     signatures,
                                     ref_signatures = NULL,
                                     method = "zscore",
                                     cor_method = "spearman",
                                     sd_handling = "exclude",
                                     apply_log2 = TRUE,
                                     min_genes = 3L,
                                     overlap_threshold = 0.50) {
  if (is.null(names(datasets)) || any(!nzchar(names(datasets)))) {
    stop("`datasets` must be a named list.", call. = FALSE)
  }
  if (is.null(names(signatures)) || any(!nzchar(names(signatures)))) {
    stop("`signatures` must be a named list.", call. = FALSE)
  }

  if (is.null(ref_signatures)) {
    ref_signatures <- .default_reference_axes()
  }
  if (is.null(names(ref_signatures)) || any(!nzchar(names(ref_signatures)))) {
    stop("`ref_signatures` must be a named list.", call. = FALSE)
  }

  rows <- list()
  row_i <- 0L

  for (ds_name in names(datasets)) {
    prep <- tryCatch(
      .prepare_cohort(
        datasets[[ds_name]],
        name = ds_name,
        apply_log2 = apply_log2,
        sd_handling = sd_handling,
        filter_response = TRUE
      ),
      error = function(e) e
    )
    if (inherits(prep, "error")) {
      if (length(datasets) == 1L) {
        stop(conditionMessage(prep), call. = FALSE)
      }
      warning(
        "Skipping cohort '", ds_name, "' in correlate_reference_axes(): ",
        conditionMessage(prep),
        call. = FALSE
      )
      next
    }
    ds <- prep$ds
    expr_full <- prep$expr
    y_kept <- prep$y_kept
    norm <- if (!is.null(ds$info$normalization) &&
                length(ds$info$normalization) > 0L &&
                !is.na(ds$info$normalization[1])) {
      as.character(ds$info$normalization)[1]
    } else {
      "unknown"
    }

    user_scores <- .score_named_signatures(
      expr_full, signatures, method = method, min_genes = min_genes
    )
    ref_scores <- .score_named_signatures(
      expr_full, ref_signatures, method = method, min_genes = min_genes
    )

    for (sig_name in names(user_scores)) {
      s_user <- user_scores[[sig_name]]
      u_genes <- .extract_sig_genes_upper(signatures[[sig_name]])
      for (axis_name in names(ref_scores)) {
        s_ref <- ref_scores[[axis_name]]
        r_genes <- .extract_sig_genes_upper(ref_signatures[[axis_name]])
        common <- intersect(intersect(names(s_user), names(s_ref)), names(y_kept))
        if (length(common) < 3L) {
          next
        }
        x <- as.numeric(s_user[common])
        y <- as.numeric(s_ref[common])
        ok <- is.finite(x) & is.finite(y)
        n_ok <- sum(ok)
        if (n_ok < 3L) {
          next
        }
        x <- x[ok]
        y <- y[ok]
        cor_val <- stats::cor(x, y, method = cor_method, use = "pairwise.complete.obs")
        p_val <- tryCatch(
          {
            ct <- stats::cor.test(x, y, method = cor_method, exact = FALSE)
            as.numeric(ct$p.value)
          },
          error = function(e) NA_real_
        )
        inter_n <- length(intersect(u_genes, r_genes))
        union_n <- length(union(u_genes, r_genes))
        jacc <- if (union_n == 0L) 0 else inter_n / union_n
        min_len <- min(length(u_genes), length(r_genes))
        min_frac <- if (min_len == 0L) 0 else inter_n / min_len
        ov_flag <- is.finite(jacc) && (jacc >= overlap_threshold || min_frac >= overlap_threshold)
        ov_note <- if (isTRUE(ov_flag)) {
          sprintf("High gene overlap with %s (Jaccard=%.2f); correlation is circular", axis_name, jacc)
        } else {
          NA_character_
        }
        row_i <- row_i + 1L
        rows[[row_i]] <- data.frame(
          cohort = ds_name,
          signature = sig_name,
          reference_axis = axis_name,
          cor = as.numeric(cor_val),
          p_value = p_val,
          n_samples = as.integer(n_ok),
          normalization = norm,
          jaccard_overlap = as.numeric(jacc),
          overlap_flag = isTRUE(ov_flag),
          overlap_note = ov_note,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (length(rows) == 0L) {
    return(data.frame(
      cohort = character(),
      signature = character(),
      reference_axis = character(),
      cor = numeric(),
      p_value = numeric(),
      p_adj_bh = numeric(),
      n_samples = integer(),
      normalization = character(),
      jaccard_overlap = numeric(),
      overlap_flag = logical(),
      overlap_note = character(),
      stringsAsFactors = FALSE
    ))
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out$p_adj_bh <- NA_real_
  ok_p <- is.finite(out$p_value)
  if (any(ok_p)) {
    out$p_adj_bh[ok_p] <- stats::p.adjust(out$p_value[ok_p], method = "BH")
  }
  col_ord <- c(
    "cohort", "signature", "reference_axis", "cor", "p_value", "p_adj_bh",
    "n_samples", "normalization", "jaccard_overlap", "overlap_flag", "overlap_note"
  )
  out[, intersect(col_ord, names(out)), drop = FALSE]
}

#' Built-in reference immune axes
#'
#' @return Named signature list.
#' @noRd
.default_reference_axes <- function() {
  list(
    IFN_Gamma = list(
      genes = c("IDO1", "CXCL9", "CXCL10", "HLA-DRA", "STAT1", "IFNG")
    ),
    Cytotoxicity = list(
      genes = c("GZMA", "GZMB", "PRF1")
    ),
    T_Cell_Infiltration = list(
      genes = c("CD3D", "CD3E", "CD8A", "CD8B")
    )
  )
}

#' Score a named signature list on an expression matrix
#'
#' @return Named list of named numeric score vectors (skipped signatures omitted).
#' @noRd
.score_named_signatures <- function(expr, signatures, method, min_genes) {
  out <- list()
  for (nm in names(signatures)) {
    sig <- .as_signature(signatures[[nm]], name = nm)
    eff_min_genes <- .resolve_min_genes(sig, default = min_genes)
    scored <- score_signature(
      expr,
      genes = sig$genes,
      genes_up = sig$genes_up,
      genes_down = sig$genes_down,
      weights = sig$weights,
      method = method,
      min_genes = eff_min_genes
    )
    if (isTRUE(scored$skipped) || is.null(scored$scores)) {
      next
    }
    out[[nm]] <- scored$scores
  }
  out
}

#' @noRd
.extract_sig_genes_upper <- function(sig) {
  if (is.character(sig)) {
    g <- sig
  } else if (is.numeric(sig)) {
    g <- names(sig)
  } else if (is.list(sig)) {
    g <- if (!is.null(sig$genes)) {
      if (is.numeric(sig$genes)) names(sig$genes) else sig$genes
    } else {
      c(sig$genes_up, sig$genes_down)
    }
  } else {
    g <- character(0)
  }
  g <- as.character(g)
  unique(toupper(g[!is.na(g) & nzchar(g)]))
}
