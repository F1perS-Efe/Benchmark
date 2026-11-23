#' Veri seti listesinin zorunlu alanlarini dogrula
#'
#' Beklenen yapi:
#' `list(expr = <matris>, sample_meta = <data.frame>, info = list(...))`
#'
#' @param dataset SigBenchR veri seti listesi.
#' @param name Hata mesajlari icin istege bagli ad.
#'
#' @return Gorunmez olarak `dataset` (yan etkisiz dogrulama).
#' @export
validate_dataset <- function(dataset, name = NULL) {
  label <- if (is.null(name)) "dataset" else name
  if (!is.list(dataset)) {
    stop("`", label, "` must be a list.", call. = FALSE)
  }
  required <- c("expr", "sample_meta", "info")
  missing <- setdiff(required, names(dataset))
  if (length(missing) > 0L) {
    stop(
      "`", label, "` is missing required fields: ", paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  if (!is.matrix(dataset$expr) && !is.data.frame(dataset$expr)) {
    stop("`", label, "$expr` must be a matrix or data.frame.", call. = FALSE)
  }
  expr_mat <- if (is.matrix(dataset$expr)) dataset$expr else as.matrix(dataset$expr)
  if (!is.numeric(expr_mat)) {
    stop("`", label, "$expr` must be a numeric matrix.", call. = FALSE)
  }
  if (nrow(expr_mat) < 1L || ncol(expr_mat) < 1L) {
    stop("`", label, "$expr` must have at least one row and one column.", call. = FALSE)
  }
  if (any(is.infinite(expr_mat))) {
    stop("`", label, "$expr` cannot contain infinite (Inf / -Inf) values.", call. = FALSE)
  }
  if (!any(is.finite(expr_mat))) {
    stop("`", label, "$expr` cannot be completely NA/NaN (at least one finite expression value is required).", call. = FALSE)
  }
  if (is.null(rownames(dataset$expr)) || any(is.na(rownames(dataset$expr)) | !nzchar(trimws(rownames(dataset$expr))))) {
    stop("`", label, "$expr` requires non-empty rownames (gene symbols).", call. = FALSE)
  }
  if (is.null(colnames(dataset$expr)) || any(is.na(colnames(dataset$expr)) | !nzchar(trimws(colnames(dataset$expr))))) {
    stop("`", label, "$expr` requires non-empty colnames (sample IDs).", call. = FALSE)
  }
  if (anyDuplicated(colnames(dataset$expr))) {
    stop("`", label, "$expr` colnames (sample IDs) must be unique.", call. = FALSE)
  }
  if (anyDuplicated(toupper(rownames(dataset$expr)))) {
    warning(
      "`", label, "$expr` contains duplicate gene rownames; ",
      "merge them with collapse_duplicate_genes().",
      call. = FALSE
    )
  }
  if (!is.data.frame(dataset$sample_meta)) {
    stop("`", label, "$sample_meta` must be a data.frame.", call. = FALSE)
  }
  if (nrow(dataset$sample_meta) < 1L) {
    stop("`", label, "$sample_meta` cannot be empty (0 rows).", call. = FALSE)
  }
  if (!"sample_id" %in% names(dataset$sample_meta)) {
    stop("`", label, "$sample_meta` must contain a `sample_id` column.", call. = FALSE)
  }
  if (anyDuplicated(as.character(dataset$sample_meta$sample_id))) {
    stop("`", label, "$sample_meta$sample_id` values must be unique.", call. = FALSE)
  }
  if (!"response_label" %in% names(dataset$sample_meta)) {
    stop(
      "`", label, "$sample_meta` must contain a `response_label` column.",
      call. = FALSE
    )
  }
  info_req <- c("normalization", "scale")
  info_miss <- setdiff(info_req, names(dataset$info))
  if (length(info_miss) > 0L) {
    stop(
      "`", label, "$info` is missing required fields: ", paste(info_miss, collapse = ", "),
      call. = FALSE
    )
  }
  if (!dataset$info$scale %in% c("log2", "linear")) {
    stop("`", label, "$info$scale` must be 'log2' or 'linear'.", call. = FALSE)
  }

  # Ornek sirasi: expr sutunlari ile sample_meta$sample_id hizali olmali.
  if (!identical(as.character(dataset$sample_meta$sample_id), colnames(dataset$expr))) {
    stop(
      "`", label, "`: sample_meta$sample_id must match expr colnames in the same order.",
      call. = FALSE
    )
  }

  # Opsiyonel sagkalim sutunlari (yoksa sessiz; varsa dogrula)
  .validate_optional_survival_columns(dataset$sample_meta, label = label)

  invisible(dataset)
}

#' Opsiyonel OS/PFS sagkalim sutunlarini dogrula
#'
#' `os_time`/`os_event` ve/veya `pfs_time`/`pfs_event` varsa kurallari uygular.
#' Sutun yoksa hata/uyari uretmez (geriye donuk uyum).
#'
#' @param sample_meta Metadata data.frame.
#' @param label Hata mesaji etiketi.
#' @return Gorunmez `TRUE`.
#' @noRd
.validate_optional_survival_columns <- function(sample_meta, label = "dataset") {
  time_cols <- intersect(c("os_time", "pfs_time"), names(sample_meta))
  event_cols <- intersect(c("os_event", "pfs_event"), names(sample_meta))

  for (tc in time_cols) {
    x <- sample_meta[[tc]]
    if (!is.numeric(x)) {
      stop(
        "`", label, "$sample_meta$", tc, "` must be numeric.",
        call. = FALSE
      )
    }
    bad <- !is.na(x) & x < 0
    if (any(bad)) {
      stop(
        "`", label, "$sample_meta$", tc,
        "` cannot be negative (>= 0, excluding NA).",
        call. = FALSE
      )
    }
  }

  for (ec in event_cols) {
    x <- sample_meta[[ec]]
    # Coerce factors/characters carefully: only 0/1/NA allowed after numeric
    x_num <- suppressWarnings(as.numeric(as.character(x)))
    if (length(x_num) != length(x) || any(!is.na(x) & is.na(x_num))) {
      stop(
        "`", label, "$sample_meta$", ec,
        "` must contain only 0, 1, or NA.",
        call. = FALSE
      )
    }
    ok <- is.na(x_num) | x_num %in% c(0, 1)
    if (!all(ok)) {
      stop(
        "`", label, "$sample_meta$", ec,
        "` must contain only 0, 1, or NA.",
        call. = FALSE
      )
    }
  }

  invisible(TRUE)
}

#' Yinelenen gen satirlarini birlestir
#'
#' Varsayilan: ortalama ifadesi en yuksek olan satiri tut.
#'
#' @param expr Gen x ornek matris.
#' @param rule Simdilik yalniz `"max_mean"`.
#'
#' @return Tekillestirilmis matris.
#' @export
collapse_duplicate_genes <- function(expr, rule = c("max_mean")) {
  rule <- match.arg(rule)
  if (is.null(rownames(expr))) {
    stop("`expr` rownames are required.", call. = FALSE)
  }
  expr <- as.matrix(expr)
  keys <- toupper(rownames(expr))
  if (!anyDuplicated(keys)) {
    return(expr)
  }

  split_idx <- split(seq_len(nrow(expr)), keys)
  keep <- vapply(split_idx, function(idx) {
    if (length(idx) == 1L) {
      return(idx)
    }
    means <- rowMeans(expr[idx, , drop = FALSE], na.rm = TRUE)
    if (length(means) == 0L || all(!is.finite(means))) {
      return(idx[1L])
    }
    wi <- which.max(means)
    if (length(wi) == 0L) {
      return(idx[1L])
    }
    idx[wi]
  }, integer(1))

  expr[sort(unname(keep)), , drop = FALSE]
}

#' Validate, deduplicate, scale-prepare, and optionally response-filter a cohort
#'
#' @noRd
.prepare_cohort <- function(dataset,
                            name = NULL,
                            apply_log2 = TRUE,
                            sd_handling = "exclude",
                            filter_response = TRUE) {
  ds <- validate_dataset(dataset, name = name)
  expr <- if (is.matrix(ds$expr)) ds$expr else as.matrix(ds$expr)
  is_mouse <- !is.null(ds$info$species) &&
    identical(tolower(as.character(ds$info$species)[1]), "mouse")

  if (is_mouse) {
    orth_tbl <- ds$info$ortholog_table
    orth_pol <- if (!is.null(ds$info$ortholog_policy)) {
      as.character(ds$info$ortholog_policy)[1]
    } else {
      "strict_1to1"
    }
    keep_ctx <- if (!is.null(ds$info$keep_unmapped_context)) {
      isTRUE(ds$info$keep_unmapped_context)
    } else {
      TRUE
    }
    expr <- map_mouse_to_human(
      expr,
      ortholog_table = orth_tbl,
      ortholog_policy = orth_pol,
      keep_unmapped_context = keep_ctx
    )
  } else {
    expr <- collapse_duplicate_genes(expr)
  }
  ostat <- attr(expr, "ortholog_stats")
  bg_expr <- attr(expr, "transcriptome_bg_expr")
  expr <- prepare_expression(
    expr,
    scale = ds$info$scale,
    apply_log2 = apply_log2
  )
  if (is.matrix(bg_expr)) {
    bg_expr <- prepare_expression(
      bg_expr,
      scale = ds$info$scale,
      apply_log2 = apply_log2
    )
    attr(expr, "transcriptome_bg_expr") <- bg_expr
  }
  if (!is.null(ostat)) {
    attr(expr, "ortholog_stats") <- ostat
  }

  if (!isTRUE(filter_response)) {
    return(list(
      ds = ds,
      expr = expr,
      expr_kept = expr,
      y_kept = NULL,
      enc = NULL,
      ortholog_stats = ostat
    ))
  }

  enc <- encode_response(
    ds$sample_meta$response_label,
    sd_handling = sd_handling
  )
  expr_kept <- expr[, enc$kept, drop = FALSE]
  if (!is.null(ostat)) {
    attr(expr_kept, "ortholog_stats") <- ostat
  }
  if (is.matrix(bg_expr) && ncol(bg_expr) == ncol(expr)) {
    attr(expr_kept, "transcriptome_bg_expr") <- bg_expr[, enc$kept, drop = FALSE]
  }
  y_kept <- enc$response[enc$kept]
  names(y_kept) <- colnames(expr_kept)

  list(
    ds = ds,
    expr = expr,
    expr_kept = expr_kept,
    y_kept = y_kept,
    enc = enc,
    ortholog_stats = ostat
  )
}
