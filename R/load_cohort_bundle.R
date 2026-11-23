#' Load a local full-gene cohort bundle (.rds) with optional MD5 check
#'
#' Reads a named list of SigBenchR datasets from a local `.rds` file.
#' Does not perform any network calls; the user supplies a local path
#' (optionally after downloading from Zenodo / GitHub Releases).
#'
#' @param path Path to a local `.rds` file containing a named list of
#'   `list(expr, sample_meta, info)` cohorts.
#' @param expected_md5 Optional MD5 hex digest. When supplied, verified via
#'   `tools::md5sum(path)`; mismatch raises an ASCII `stop()`.
#'
#' @return A named list of validated cohort datasets.
#' @export
#' @importFrom tools md5sum
load_cohort_bundle <- function(path, expected_md5 = NULL) {
  if (missing(path) || !nzchar(as.character(path)[1])) {
    stop("`path` is required.", call. = FALSE)
  }
  path <- as.character(path)[1]
  if (!file.exists(path)) {
    stop("Cohort bundle not found: ", path, call. = FALSE)
  }

  if (!is.null(expected_md5)) {
    expected_md5 <- tolower(as.character(expected_md5)[1])
    got <- tolower(unname(tools::md5sum(path)))
    if (!identical(got, expected_md5)) {
      stop(
        "MD5 mismatch: expected=", expected_md5, ", got=", got, ".",
        call. = FALSE
      )
    }
  }

  bundle <- readRDS(path)
  if (!is.list(bundle) || is.null(names(bundle)) || any(!nzchar(names(bundle)))) {
    stop(
      "Cohort bundle must be a named list of dataset objects.",
      call. = FALSE
    )
  }

  out <- lapply(names(bundle), function(nm) {
    validate_dataset(bundle[[nm]], name = nm)
  })
  names(out) <- names(bundle)
  out
}
