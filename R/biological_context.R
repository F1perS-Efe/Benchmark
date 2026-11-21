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
#'
#' @return `data.frame` with columns `cohort`, `signature`, `reference_axis`,
#'   `cor`, `p_value`, `n_samples`, `normalization`.
#' @export
correlate_reference_axes <- function(datasets,
                                     signatures,
                                     ref_signatures = NULL,
                                     method = "zscore",
                                     cor_method = "spearman",
                                     sd_handling = "exclude",
                                     apply_log2 = TRUE,
                                     min_genes = 3L) {
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
    ds <- validate_dataset(datasets[[ds_name]], name = ds_name)
    expr <- collapse_duplicate_genes(as.matrix(ds$expr))
    expr <- prepare_expression(
      expr,
      scale = ds$info$scale,
      apply_log2 = apply_log2
    )
    enc <- encode_response(
      ds$sample_meta$response_label,
      sd_handling = sd_handling
    )
    expr_kept <- expr[, enc$kept, drop = FALSE]
    norm <- if (!is.null(ds$info$normalization) &&
                length(ds$info$normalization) > 0L &&
                !is.na(ds$info$normalization[1])) {
      as.character(ds$info$normalization)[1]
    } else {
      "unknown"
    }

    user_scores <- .score_named_signatures(
      expr_kept, signatures, method = method, min_genes = min_genes
    )
    ref_scores <- .score_named_signatures(
      expr_kept, ref_signatures, method = method, min_genes = min_genes
    )

    for (sig_name in names(user_scores)) {
      s_user <- user_scores[[sig_name]]
      for (axis_name in names(ref_scores)) {
        s_ref <- ref_scores[[axis_name]]
        common <- intersect(names(s_user), names(s_ref))
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
        row_i <- row_i + 1L
        rows[[row_i]] <- data.frame(
          cohort = ds_name,
          signature = sig_name,
          reference_axis = axis_name,
          cor = as.numeric(cor_val),
          p_value = p_val,
          n_samples = as.integer(n_ok),
          normalization = norm,
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
      n_samples = integer(),
      normalization = character(),
      stringsAsFactors = FALSE
    ))
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
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
    sig <- signatures[[nm]]
    if (is.character(sig) || is.numeric(sig)) {
      sig <- list(genes = sig)
    }
    scored <- score_signature(
      expr,
      genes = sig$genes,
      genes_up = sig$genes_up,
      genes_down = sig$genes_down,
      method = method,
      min_genes = min_genes
    )
    if (isTRUE(scored$skipped) || is.null(scored$scores)) {
      next
    }
    out[[nm]] <- scored$scores
  }
  out
}
