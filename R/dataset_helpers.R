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
    stop("`", label, "` bir list olmalidir.", call. = FALSE)
  }
  required <- c("expr", "sample_meta", "info")
  missing <- setdiff(required, names(dataset))
  if (length(missing) > 0L) {
    stop(
      "`", label, "` eksik alanlar: ", paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  if (!is.matrix(dataset$expr) && !is.data.frame(dataset$expr)) {
    stop("`", label, "$expr` matris veya data.frame olmalidir.", call. = FALSE)
  }
  expr_mat <- as.matrix(dataset$expr)
  if (!is.numeric(expr_mat)) {
    stop("`", label, "$expr` sayisal bir matris olmalidir.", call. = FALSE)
  }
  if (nrow(expr_mat) < 1L || ncol(expr_mat) < 1L) {
    stop("`", label, "$expr` en az bir satir ve bir sutun icermelidir.", call. = FALSE)
  }
  if (is.null(rownames(dataset$expr))) {
    stop("`", label, "$expr` satir adlari (genler) zorunludur.", call. = FALSE)
  }
  if (is.null(colnames(dataset$expr))) {
    stop("`", label, "$expr` sutun adlari (ornekler) zorunludur.", call. = FALSE)
  }
  if (anyDuplicated(colnames(dataset$expr))) {
    stop("`", label, "$expr` sutun adlari (ornek id) tekil olmalidir.", call. = FALSE)
  }
  if (anyDuplicated(rownames(dataset$expr))) {
    warning(
      "`", label, "$expr` satir adlarinda yinelenen genler var; ",
      "collapse_duplicate_genes() ile birlestirin.",
      call. = FALSE
    )
  }
  if (!is.data.frame(dataset$sample_meta)) {
    stop("`", label, "$sample_meta` data.frame olmalidir.", call. = FALSE)
  }
  if (!"sample_id" %in% names(dataset$sample_meta)) {
    stop("`", label, "$sample_meta` icinde `sample_id` sutunu zorunludur.", call. = FALSE)
  }
  if (anyDuplicated(as.character(dataset$sample_meta$sample_id))) {
    stop("`", label, "$sample_meta$sample_id` degerleri tekil olmalidir.", call. = FALSE)
  }
  if (!"response_label" %in% names(dataset$sample_meta)) {
    stop(
      "`", label, "$sample_meta` icinde `response_label` sutunu zorunludur.",
      call. = FALSE
    )
  }
  info_req <- c("normalization", "scale")
  info_miss <- setdiff(info_req, names(dataset$info))
  if (length(info_miss) > 0L) {
    stop(
      "`", label, "$info` eksik alanlar: ", paste(info_miss, collapse = ", "),
      call. = FALSE
    )
  }
  if (!dataset$info$scale %in% c("log2", "linear")) {
    stop("`", label, "$info$scale` 'log2' veya 'linear' olmalidir.", call. = FALSE)
  }

  # Ornek sirasi: expr sutunlari ile sample_meta$sample_id hizali olmali.
  if (!identical(as.character(dataset$sample_meta$sample_id), colnames(dataset$expr))) {
    stop(
      "`", label, "`: sample_meta$sample_id, expr sutun adlariyla ayni sirada olmalidir.",
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
        "`", label, "$sample_meta$", tc, "` sayisal olmalidir.",
        call. = FALSE
      )
    }
    bad <- !is.na(x) & x < 0
    if (any(bad)) {
      stop(
        "`", label, "$sample_meta$", tc,
        "` negatif olamaz (>= 0; NA haric).",
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
        "` yalnizca 0, 1 veya NA icermelidir.",
        call. = FALSE
      )
    }
    ok <- is.na(x_num) | x_num %in% c(0, 1)
    if (!all(ok)) {
      stop(
        "`", label, "$sample_meta$", ec,
        "` yalnizca 0, 1 veya NA icermelidir.",
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
    stop("`expr` satir adlari zorunludur.", call. = FALSE)
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
