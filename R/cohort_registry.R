#' List cohorts in the SigBenchR cohort registry
#'
#' Reads the bundled `inst/extdata/cohort_registry.yml` (or a user-supplied
#' YAML registry file) and returns a summary `data.frame`.
#'
#' @param registry_path Optional path to a custom YAML registry file. `NULL`
#'   reads the bundled builtin registry.
#'
#' @return A `data.frame` of registered cohorts.
#' @export
list_cohort_registry <- function(registry_path = NULL) {
  reg <- .read_registry_yaml(registry_path)
  cohorts <- reg$cohorts
  if (is.null(cohorts) || length(cohorts) == 0L) {
    return(data.frame(
      cohort_id = character(), source_type = character(), species = character(),
      cancer_type = character(), treatment = character(), endpoint_type = character(),
      normalization = character(), scale = character(), n_samples = integer(),
      n_genes = integer(), file_path = character(), meta_path = character(),
      md5 = character(), stringsAsFactors = FALSE
    ))
  }
  rows <- lapply(names(cohorts), function(cid) {
    x <- cohorts[[cid]]
    pick_chr <- function(v) if (is.null(v) || length(v) == 0L || is.na(v)) "" else as.character(v)[1]
    pick_int <- function(v) if (is.null(v) || length(v) == 0L || is.na(v)) NA_integer_ else as.integer(v)[1]
    data.frame(
      cohort_id = cid,
      source_type = pick_chr(x$source_type),
      species = pick_chr(x$species),
      cancer_type = pick_chr(x$cancer_type),
      treatment = pick_chr(x$treatment),
      endpoint_type = pick_chr(x$endpoint_type),
      normalization = pick_chr(x$normalization),
      scale = pick_chr(x$scale),
      n_samples = pick_int(x$n_samples),
      n_genes = pick_int(x$n_genes),
      file_path = pick_chr(x$file_path),
      meta_path = pick_chr(x$meta_path),
      md5 = pick_chr(x$md5),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Register a local full-transcriptome cohort (.rds or .csv/.tsv) in a YAML registry
#'
#' Validates the cohort dataset, computes its file MD5 checksum, normalizes the
#' treatment label via `.normalize_treatment_label()`, and stores the entry in
#' `registry_path`.
#'
#' @param cohort_id Unique identifier for the cohort.
#' @param file_path Path to local `.rds` dataset or expression `.csv`/`.tsv`.
#' @param registry_path Path to the YAML registry file to create/update.
#' @param meta_path Path to clinical metadata CSV/TSV (required when `file_path`
#'   is a CSV/TSV expression matrix).
#' @param species `"human"` (default) or `"mouse"`.
#' @param cancer_type Cancer type label (e.g. `"melanoma"`, `"NSCLC"`).
#' @param treatment Treatment regimen (normalized automatically).
#' @param endpoint_type Clinical endpoint label (e.g. `"RECIST"`).
#' @param normalization Expression normalization (for CSV/TSV inputs).
#' @param scale `"log2"` or `"linear"` (for CSV/TSV inputs).
#' @param include_builtin If `TRUE` (default) and `registry_path` does not yet
#'   exist, seed it with the bundled builtin cohorts.
#'
#' @return The newly registered cohort metadata row (invisibly).
#' @export
register_cohort <- function(cohort_id,
                            file_path,
                            registry_path,
                            meta_path = NULL,
                            species = c("human", "mouse"),
                            cancer_type = NA_character_,
                            treatment = NA_character_,
                            endpoint_type = NA_character_,
                            normalization = "TPM",
                            scale = "log2",
                            include_builtin = TRUE) {
  species <- match.arg(species)
  cohort_id <- trimws(as.character(cohort_id)[1])
  if (is.na(cohort_id) || !nzchar(cohort_id)) {
    stop("`cohort_id` is required.", call. = FALSE)
  }
  if (missing(file_path) || !file.exists(file_path)) {
    stop("Cohort file not found: ", file_path, call. = FALSE)
  }
  if (missing(registry_path) || !nzchar(as.character(registry_path)[1])) {
    stop("`registry_path` is required.", call. = FALSE)
  }
  file_path <- normalizePath(file_path, mustWork = TRUE)
  reg_dir_raw <- dirname(as.character(registry_path)[1])
  if (!dir.exists(reg_dir_raw)) {
    dir.create(reg_dir_raw, recursive = TRUE, showWarnings = FALSE)
  }
  reg_dir <- normalizePath(reg_dir_raw, mustWork = TRUE)
  registry_path <- file.path(reg_dir, basename(as.character(registry_path)[1]))

  ext <- tolower(tools::file_ext(file_path))
  if (identical(ext, "rds")) {
    ds <- readRDS(file_path)
    ds <- validate_dataset(ds, name = cohort_id)
    src_type <- "local_rds"
    norm_val <- as.character(ds$info$normalization)[1]
    scale_val <- as.character(ds$info$scale)[1]
    meta_norm <- ""
  } else {
    if (is.null(meta_path) || !file.exists(meta_path)) {
      stop("`meta_path` is required for CSV/TSV cohort registration.", call. = FALSE)
    }
    meta_norm <- normalizePath(meta_path, mustWork = TRUE)
    ds <- read_user_dataset(
      expr_path = file_path,
      meta_path = meta_norm,
      normalization = normalization,
      scale = scale,
      dataset_name = cohort_id
    )
    src_type <- "local_csv"
    norm_val <- as.character(ds$info$normalization)[1]
    scale_val <- as.character(ds$info$scale)[1]
  }

  tx_norm <- if (!is.na(treatment) && nzchar(trimws(as.character(treatment)[1]))) {
    .normalize_treatment_label(as.character(treatment)[1])
  } else if (!is.null(ds$info$treatment)) {
    .normalize_treatment_label(as.character(ds$info$treatment)[1])
  } else {
    ""
  }

  reg <- if (file.exists(registry_path)) {
    yaml::read_yaml(registry_path)
  } else if (isTRUE(include_builtin)) {
    .read_registry_yaml(NULL)
  } else {
    list(cohorts = list())
  }
  if (is.null(reg$cohorts)) {
    reg$cohorts <- list()
  }

  md5_val <- unname(tools::md5sum(file_path))
  file_path_store <- .make_relative_path(file_path, reg_dir)
  meta_path_store <- if (nzchar(meta_norm)) .make_relative_path(meta_norm, reg_dir) else ""
  reg$cohorts[[cohort_id]] <- list(
    source_type = src_type,
    species = species,
    cancer_type = if (is.na(cancer_type)) "" else as.character(cancer_type)[1],
    treatment = tx_norm,
    endpoint_type = if (is.na(endpoint_type)) "" else as.character(endpoint_type)[1],
    normalization = norm_val,
    scale = scale_val,
    n_samples = as.integer(ncol(ds$expr)),
    n_genes = as.integer(nrow(ds$expr)),
    file_path = file_path_store,
    meta_path = meta_path_store,
    md5 = md5_val
  )

  yaml::write_yaml(reg, registry_path)

  tbl <- list_cohort_registry(registry_path)
  invisible(tbl[tbl$cohort_id == cohort_id, , drop = FALSE])
}

#' Load a cohort from the registry with optional MD5 verification
#'
#' @param cohort_id Cohort identifier in the registry.
#' @param registry_path Optional path to custom YAML registry. `NULL` uses the
#'   bundled registry.
#' @param verify_md5 If `TRUE` (default), verify file MD5 checksum before
#'   loading local cohorts.
#'
#' @return A validated SigBenchR dataset list (`list(expr, sample_meta, info)`).
#' @export
load_registered_cohort <- function(cohort_id,
                                   registry_path = NULL,
                                   verify_md5 = TRUE) {
  cohort_id <- as.character(cohort_id)[1]
  tbl <- list_cohort_registry(registry_path)
  idx <- match(cohort_id, tbl$cohort_id)
  if (is.na(idx)) {
    stop("Cohort '", cohort_id, "' not found in registry.", call. = FALSE)
  }
  row <- tbl[idx, , drop = FALSE]

  if (identical(row$source_type, "builtin")) {
    env <- new.env(parent = emptyenv())
    utils::data("datasets_ici", package = "SigBenchR", envir = env)
    if (!cohort_id %in% names(env$datasets_ici)) {
      stop("Builtin cohort '", cohort_id, "' missing from datasets_ici.", call. = FALSE)
    }
    ds <- env$datasets_ici[[cohort_id]]
    ds$info$species <- if (nzchar(row$species)) row$species else "human"
    return(validate_dataset(ds, name = cohort_id))
  }

  fp <- .resolve_registry_path(row$file_path, registry_path = registry_path)
  if (!nzchar(fp) || !file.exists(fp)) {
    stop("Registered cohort file not found: ", fp, call. = FALSE)
  }
  if (isTRUE(verify_md5) && nzchar(row$md5)) {
    got_md5 <- tolower(unname(tools::md5sum(fp)))
    if (!identical(got_md5, tolower(row$md5))) {
      stop(
        "MD5 mismatch for cohort '", cohort_id, "': expected=", row$md5,
        ", got=", got_md5, ".",
        call. = FALSE
      )
    }
  }

  if (identical(row$source_type, "local_rds") || identical(tolower(tools::file_ext(fp)), "rds")) {
    ds <- readRDS(fp)
  } else {
    mp <- .resolve_registry_path(row$meta_path, registry_path = registry_path)
    ds <- read_user_dataset(
      expr_path = fp,
      meta_path = mp,
      normalization = if (nzchar(row$normalization)) row$normalization else "TPM",
      scale = if (nzchar(row$scale)) row$scale else "log2",
      dataset_name = cohort_id
    )
  }
  ds$info$species <- if (nzchar(row$species)) row$species else "human"
  if (nzchar(row$treatment) && is.null(ds$info$treatment)) {
    ds$info$treatment <- row$treatment
  }
  if (nzchar(row$endpoint_type) && is.null(ds$info$endpoint_type)) {
    ds$info$endpoint_type <- row$endpoint_type
  }
  if (nzchar(row$cancer_type) && is.null(ds$info$cancer_type)) {
    ds$info$cancer_type <- row$cancer_type
  }
  validate_dataset(ds, name = cohort_id)
}

#' @noRd
.read_registry_yaml <- function(registry_path = NULL) {
  if (is.null(registry_path) || !nzchar(as.character(registry_path)[1])) {
    registry_path <- system.file("extdata", "cohort_registry.yml", package = "SigBenchR")
  }
  if (!nzchar(registry_path) || !file.exists(registry_path)) {
    stop("Cohort registry YAML file not found in installed package (extdata/cohort_registry.yml).", call. = FALSE)
  }
  yaml::read_yaml(registry_path)
}

#' @noRd
.make_relative_path <- function(target, base_dir) {
  target <- gsub("\\\\", "/", as.character(target)[1])
  base_dir <- gsub("\\\\", "/", as.character(base_dir)[1])
  if (nchar(base_dir) > 1L) {
    base_dir <- sub("/+$", "", base_dir)
  }

  # Windows drive letter check (e.g. C:/ vs D:/): different drives cannot be relative
  t_drive <- if (grepl("^[A-Za-z]:(/|$)", target)) toupper(substr(target, 1L, 2L)) else ""
  b_drive <- if (grepl("^[A-Za-z]:(/|$)", base_dir)) toupper(substr(base_dir, 1L, 2L)) else ""
  if (nzchar(t_drive) != nzchar(b_drive)) {
    return(target)
  }
  if (nzchar(t_drive) && !identical(t_drive, b_drive)) {
    return(target)
  }

  t_parts <- strsplit(target, "/", fixed = TRUE)[[1]]
  b_parts <- strsplit(base_dir, "/", fixed = TRUE)[[1]]
  if (nzchar(t_drive) && length(t_parts) >= 1L && length(b_parts) >= 1L) {
    t_parts[1L] <- toupper(t_parts[1L])
    b_parts[1L] <- toupper(b_parts[1L])
  }

  n_min <- min(length(t_parts), length(b_parts))
  common <- 0L
  for (i in seq_len(n_min)) {
    if (identical(t_parts[i], b_parts[i])) {
      common <- i
    } else {
      break
    }
  }
  if (common < 2L || common >= length(t_parts)) {
    return(target)
  }
  up <- rep("..", length(b_parts) - common)
  down <- t_parts[seq.int(common + 1L, length(t_parts))]
  paste(c(up, down), collapse = "/")
}

#' @noRd
.resolve_registry_path <- function(path_val, registry_path = NULL) {
  path_val <- as.character(path_val)[1]
  if (is.na(path_val) || !nzchar(path_val)) {
    return("")
  }
  is_abs <- grepl("^(/|[A-Za-z]:[/\\\\])", path_val)
  if (is_abs || is.null(registry_path) || !nzchar(as.character(registry_path)[1])) {
    return(path_val)
  }
  reg_dir <- dirname(normalizePath(as.character(registry_path)[1], mustWork = FALSE))
  cand <- file.path(reg_dir, path_val)
  if (file.exists(cand)) cand else path_val
}
