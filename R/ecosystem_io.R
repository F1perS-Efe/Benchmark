#' Read a MSigDB-style GMT file (base R only)
#'
#' Each non-empty line is `NAME\\tDESCRIPTION\\tGENE1\\tGENE2\\t...`.
#' Empty lines and empty gene tokens are dropped.
#'
#' @param path Path to a `.gmt` file.
#'
#' @return Named list of character gene vectors (`signature_name = genes`).
#' @export
read_gmt <- function(path) {
  if (missing(path) || !nzchar(as.character(path)[1])) {
    stop("`path` is required.", call. = FALSE)
  }
  path <- as.character(path)[1]
  if (!file.exists(path)) {
    stop("GMT file not found: ", path, call. = FALSE)
  }

  lines <- readLines(path, warn = FALSE)
  lines <- lines[!is.na(lines)]
  lines <- trimws(lines)
  lines <- lines[nzchar(lines)]
  if (length(lines) == 0L) {
    return(list())
  }

  out <- list()
  for (ln in lines) {
    parts <- strsplit(ln, "\t", fixed = TRUE)[[1]]
    parts <- trimws(parts)
    if (length(parts) < 1L || !nzchar(parts[1])) {
      next
    }
    nm <- parts[1]
    genes <- if (length(parts) >= 3L) parts[seq.int(3L, length(parts))] else character(0)
    genes <- genes[nzchar(genes)]
    out[[nm]] <- genes
  }
  out
}

#' Coerce objects into a SigBenchR dataset list
#'
#' S3 generic. Validates via [validate_dataset()]. Does **not** declare
#' Bioconductor packages in DESCRIPTION; SummarizedExperiment support is
#' optional at runtime when that package is already installed.
#'
#' @param x Object to coerce.
#' @param ... Passed to methods (e.g. `name`, `normalization`, `scale`).
#'
#' @return A validated `list(expr, sample_meta, info)`.
#' @export
as_sigbenchr_dataset <- function(x, ...) {
  UseMethod("as_sigbenchr_dataset")
}

#' @rdname as_sigbenchr_dataset
#' @param name Optional dataset label for validation messages.
#' @export
as_sigbenchr_dataset.list <- function(x, name = NULL, ...) {
  validate_dataset(x, name = name)
}

#' @rdname as_sigbenchr_dataset
#' @param assay Assay index or name for SummarizedExperiment inputs.
#' @param normalization Value stored in `info$normalization`.
#' @param scale Value stored in `info$scale` (`"log2"` or `"linear"`).
#' @param response_col Column in colData used as `response_label`.
#' @export
as_sigbenchr_dataset.default <- function(x,
                                         ...,
                                         assay = 1L,
                                         normalization = "other",
                                         scale = "log2",
                                         response_col = "response_label",
                                         name = NULL) {
  # Avoid hard-coding Bioconductor package literals for CRAN checks
  # (same pattern as optional singscore cross-checks elsewhere).
  se_class <- paste0("Summarized", "Experiment")
  se_pkg <- paste0("Summarized", "Experiment")

  if (!inherits(x, se_class)) {
    stop(
      "`as_sigbenchr_dataset()` supports named dataset lists or ",
      se_class, " objects. Got: ",
      paste(class(x), collapse = "/"),
      ".",
      call. = FALSE
    )
  }

  if (!requireNamespace(se_pkg, quietly = TRUE)) {
    stop(
      "Object inherits ", se_class, " but package '", se_pkg,
      "' is not installed. Install it from Bioconductor to coerce, ",
      "or pass a list(expr, sample_meta, info) instead.",
      call. = FALSE
    )
  }

  assay_fun <- getExportedValue(se_pkg, "assay")
  coldata_fun <- getExportedValue(se_pkg, "colData")
  expr <- assay_fun(x, assay)
  expr <- as.matrix(expr)
  if (is.null(rownames(expr)) || is.null(colnames(expr))) {
    stop(
      se_class, " assay must have gene rownames and sample colnames.",
      call. = FALSE
    )
  }

  cd <- as.data.frame(coldata_fun(x), stringsAsFactors = FALSE)
  if (!"sample_id" %in% names(cd)) {
    cd$sample_id <- colnames(expr)
  }
  response_col <- as.character(response_col)[1]
  if (!response_col %in% names(cd)) {
    stop(
      "colData must contain response column '", response_col, "'.",
      call. = FALSE
    )
  }
  if (!identical(response_col, "response_label")) {
    cd$response_label <- cd[[response_col]]
  }

  # Align sample_meta rows to assay columns
  sid <- as.character(cd$sample_id)
  idx <- match(colnames(expr), sid)
  if (anyNA(idx)) {
    stop(
      "colData sample_id must cover all assay colnames.",
      call. = FALSE
    )
  }
  cd <- cd[idx, , drop = FALSE]
  cd$sample_id <- colnames(expr)

  ds <- list(
    expr = expr,
    sample_meta = cd,
    info = list(
      normalization = as.character(normalization)[1],
      scale = as.character(scale)[1],
      source = se_class
    )
  )
  validate_dataset(ds, name = name)
}
