#' Read a local CSV/TSV expression matrix into the SigBenchR contract
#'
#' Converts user-provided expression and clinical metadata files into
#' `list(expr, sample_meta, info)` and validates the result.
#'
#' @param expr_path Path to expression CSV/TSV (genes x samples by default).
#' @param meta_path Path to clinical metadata CSV/TSV.
#' @param sep Field separator. `NULL` (default) auto-detects from extension
#'   (`.tsv`/`.txt` -> tab, `.csv` -> comma).
#' @param genes_in_rows If `TRUE` (default), first column is gene IDs and
#'   remaining columns are samples. If `FALSE`, first column is sample IDs
#'   and the matrix is transposed with `t()`.
#' @param gene_id_type One of `"HGNC"`, `"Ensembl"`, `"other"`.
#' @param normalization One of `"TPM"`, `"FPKM"`, `"normalized_counts"`,
#'   `"log_normalized"`, `"other"`.
#' @param scale One of `"linear"`, `"log2"`. Declared scale is never
#'   re-transformed silently; `prepare_expression()` respects this flag.
#' @param sample_id_col Metadata column holding sample IDs.
#' @param response_col Metadata column holding response labels.
#' @param dataset_name Optional label stored in `info$dataset_name`.
#' @param collapse_duplicates If `TRUE`, collapse duplicate gene rows via
#'   `collapse_duplicate_genes()`.
#'
#' @return A validated SigBenchR dataset list.
#' @export
read_user_dataset <- function(expr_path,
                              meta_path,
                              sep = NULL,
                              genes_in_rows = TRUE,
                              gene_id_type = c("HGNC", "Ensembl", "other"),
                              normalization = c(
                                "TPM", "FPKM", "normalized_counts",
                                "log_normalized", "other"
                              ),
                              scale = c("linear", "log2"),
                              sample_id_col = "sample_id",
                              response_col = "response_label",
                              dataset_name = NULL,
                              collapse_duplicates = TRUE) {
  gene_id_type <- match.arg(gene_id_type)
  normalization <- match.arg(normalization)
  scale <- match.arg(scale)

  if (!file.exists(expr_path)) {
    stop("`expr_path` not found: ", expr_path, call. = FALSE)
  }
  if (!file.exists(meta_path)) {
    stop("`meta_path` not found: ", meta_path, call. = FALSE)
  }

  sep_expr <- .detect_sep(expr_path, sep)
  sep_meta <- .detect_sep(meta_path, sep)

  # Read without row.names so duplicate gene IDs are allowed (collapsed later).
  expr_raw <- utils::read.table(
    expr_path,
    header = TRUE,
    sep = sep_expr,
    check.names = FALSE,
    quote = "\"",
    comment.char = "",
    stringsAsFactors = FALSE,
    na.strings = c("NA", "")
  )
  meta_df <- utils::read.table(
    meta_path,
    header = TRUE,
    sep = sep_meta,
    check.names = FALSE,
    quote = "\"",
    comment.char = "",
    stringsAsFactors = FALSE,
    na.strings = c("NA", "")
  )

  if (ncol(expr_raw) < 2L) {
    stop("Expression file must contain at least one ID column and one value column.", call. = FALSE)
  }
  id_col <- as.character(expr_raw[[1L]])
  if (any(is.na(id_col) | !nzchar(id_col))) {
    stop("Expression matrix requires non-empty gene (or sample) ID rows.", call. = FALSE)
  }
  expr_df <- expr_raw[, -1L, drop = FALSE]

  if (!isTRUE(genes_in_rows)) {
    # samples in rows -> transpose to genes x samples
    mat_tmp <- as.matrix(expr_df)
    rownames(mat_tmp) <- id_col
    mat_tmp <- t(mat_tmp)
    id_col <- rownames(mat_tmp)
    expr_df <- as.data.frame(mat_tmp, stringsAsFactors = FALSE, check.names = FALSE)
  }

  # All expression values must be numeric
  for (j in seq_len(ncol(expr_df))) {
    if (!is.numeric(expr_df[[j]])) {
      coerced <- suppressWarnings(as.numeric(as.character(expr_df[[j]])))
      if (any(is.na(coerced) & !is.na(expr_df[[j]]) & nzchar(as.character(expr_df[[j]])))) {
        stop(
          "Expression matrix must be numeric; non-numeric value(s) found ",
          "(column: ", names(expr_df)[j], ").",
          call. = FALSE
        )
      }
      expr_df[[j]] <- coerced
    }
  }

  expr <- as.matrix(expr_df)
  storage.mode(expr) <- "double"
  rownames(expr) <- id_col
  if (!is.numeric(expr)) {
    stop("Expression matrix must be numeric.", call. = FALSE)
  }
  if (is.null(colnames(expr)) || any(!nzchar(colnames(expr)))) {
    stop("Expression matrix requires non-empty sample (or gene) ID columns.", call. = FALSE)
  }

  if (identical(gene_id_type, "Ensembl")) {
    warning(
      "gene_id_type = 'Ensembl': signature matching expects HGNC gene symbols; ",
      "Ensembl IDs may not match signature genes.",
      call. = FALSE
    )
  }

  n_dup_before <- sum(duplicated(toupper(rownames(expr))))
  if (n_dup_before > 0L && isTRUE(collapse_duplicates)) {
    n_unique_dup <- length(unique(toupper(rownames(expr))[duplicated(toupper(rownames(expr)))]))
    expr <- collapse_duplicate_genes(expr)
    warning(
      n_unique_dup, " duplicate gene(s) collapsed via collapse_duplicate_genes().",
      call. = FALSE
    )
  }

  if (!sample_id_col %in% names(meta_df)) {
    stop(
      "Column `", sample_id_col, "` not found in `meta_path`.",
      call. = FALSE
    )
  }
  if (!response_col %in% names(meta_df)) {
    stop(
      "Column `", response_col, "` not found in `meta_path`.",
      call. = FALSE
    )
  }

  meta_df$sample_id <- as.character(meta_df[[sample_id_col]])
  meta_df$response_label <- as.character(meta_df[[response_col]])

  if (any(is.na(meta_df$sample_id) | !nzchar(trimws(meta_df$sample_id)))) {
    stop("Metadata column `", sample_id_col, "` cannot contain empty or NA sample IDs.", call. = FALSE)
  }
  if (anyDuplicated(meta_df$sample_id)) {
    stop("Metadata column `", sample_id_col, "` contains duplicate sample_id values.", call. = FALSE)
  }
  if (anyDuplicated(as.character(colnames(expr)))) {
    stop("Expression matrix contains duplicate sample IDs in column headers.", call. = FALSE)
  }

  expr_ids <- as.character(colnames(expr))
  meta_ids <- as.character(meta_df$sample_id)
  common <- intersect(expr_ids, meta_ids)
  if (length(common) == 0L) {
    stop(
      "No common sample_id values between expression matrix and metadata.",
      call. = FALSE
    )
  }
  dropped_expr <- setdiff(expr_ids, common)
  dropped_meta <- setdiff(meta_ids, common)
  if (length(dropped_expr) > 0L || length(dropped_meta) > 0L) {
    warning(
      "Excluded unmatched samples: expr=", length(dropped_expr),
      ", meta=", length(dropped_meta), ".",
      call. = FALSE
    )
  }

  expr <- expr[, common, drop = FALSE]
  meta_idx <- match(common, meta_df$sample_id)
  sample_meta <- meta_df[meta_idx, , drop = FALSE]
  rownames(sample_meta) <- NULL

  info <- list(
    normalization = normalization,
    scale = scale,
    gene_id_type = gene_id_type,
    source = "user_local",
    n_genes = nrow(expr),
    n_samples = ncol(expr)
  )
  if (!is.null(dataset_name) && nzchar(as.character(dataset_name)[1])) {
    info$dataset_name <- as.character(dataset_name)[1]
  }

  ds <- list(
    expr = expr,
    sample_meta = sample_meta,
    info = info
  )
  ds <- validate_dataset(ds, name = if (!is.null(info$dataset_name)) info$dataset_name else "user")

  # Trigger unknown-label warning via encode_response
  encode_response(ds$sample_meta$response_label)

  ds
}

.detect_sep <- function(path, sep) {
  if (!is.null(sep)) {
    return(sep)
  }
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("tsv", "txt")) {
    return("\t")
  }
  if (identical(ext, "csv")) {
    return(",")
  }
  # Fallback: comma
  ","
}
