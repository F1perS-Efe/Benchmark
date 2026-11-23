#' Build a deterministic reproducibility manifest for a benchmark run
#'
#' Captures package versions, sanitized run parameters, interpretation rules
#' provenance, cohort/signature fingerprints, and a deterministic
#' `reproducibility_hash` (MD5) computed from all scientific inputs and
#' parameters while excluding `timestamp_utc`.
#'
#' @param results `run_benchmark()` output `data.frame`.
#' @param summary_tbl Optional `summarize_auc()` output `data.frame`.
#' @param report_model Optional `build_report_model()` output.
#' @param params_list Optional list of run parameters.
#' @param datasets Optional named dataset list for richer cohort fingerprints.
#' @param signatures Optional named signature list for richer gene fingerprints.
#'
#' @return A named list representing the reproducibility manifest, including
#'   `sigbench_version`, `timestamp_utc`, `r_version`, `package_versions`,
#'   `run_params`, `rules_meta`, `cohorts_fingerprint`,
#'   `signatures_fingerprint`, `diagnostics_summary` (summary of evaluable/skipped
#'   rows, null permutation background universe types/sizes, PR-AUC stability
#'   flags, survival, and beyond-infiltration presence), and deterministic
#'   `reproducibility_hash` (MD5).
#' @examples
#' data(datasets_demo, package = "SigBenchR")
#' data(signatures_demo, package = "SigBenchR")
#' res <- run_benchmark(
#'   signatures_demo["IFNG6"],
#'   datasets_demo["demo_A"],
#'   methods = "mean",
#'   compute_pr = FALSE,
#'   compute_null = FALSE
#' )
#' man <- build_reproducibility_manifest(res)
#' man$reproducibility_hash
#' @export
build_reproducibility_manifest <- function(results,
                                           summary_tbl = NULL,
                                           report_model = NULL,
                                           params_list = NULL,
                                           datasets = NULL,
                                           signatures = NULL) {
  if (!is.data.frame(results)) {
    stop("`results` must be a data.frame.", call. = FALSE)
  }

  if (is.null(report_model) && nrow(results) > 0L) {
    report_model <- tryCatch(
      build_report_model(
        results = results,
        summary_tbl = summary_tbl,
        params_list = if (is.null(params_list)) list() else params_list
      ),
      error = function(e) NULL
    )
  }

  raw_params <- if (!is.null(params_list)) {
    params_list
  } else if (!is.null(report_model$tier3$run_params)) {
    report_model$tier3$run_params
  } else {
    list()
  }

  param_keys <- c(
    "methods", "direction", "sd_handling", "coverage_threshold",
    "sensitivity_threshold", "n_boot", "n_perm", "seed",
    "compute_pr", "compute_null", "match_expression",
    "heatmap_method", "forest_signature", "stratify_by"
  )
  run_params <- list()
  for (k in param_keys) {
    if (!is.null(raw_params[[k]])) {
      run_params[[k]] <- raw_params[[k]]
    }
  }
  if (is.null(run_params$methods) && "method" %in% names(results) && nrow(results) > 0L) {
    run_params$methods <- sort(unique(as.character(results$method)))
  }

  rules_meta <- if (!is.null(report_model$tier3$rules_meta)) {
    list(
      version = as.character(report_model$tier3$rules_meta$version),
      md5 = as.character(report_model$tier3$rules_meta$md5)
    )
  } else {
    r <- load_interpretation_rules()
    list(version = as.character(r$version), md5 = as.character(r$md5))
  }

  cohorts_fp <- .fingerprint_cohorts(results, datasets = datasets)
  sigs_fp <- .fingerprint_signatures(results, signatures = signatures)

  core_payload <- list(
    sigbench_version = .pkg_ver("SigBenchR"),
    run_params = run_params,
    rules_meta = rules_meta,
    cohorts_fingerprint = cohorts_fp,
    signatures_fingerprint = sigs_fp
  )
  repro_hash <- .md5_string(yaml::as.yaml(core_payload))

  pkgs <- c("pROC", "survival", "yaml", "jsonlite", "PRROC", "ggplot2")
  pkg_vers <- stats::setNames(
    vapply(pkgs, .pkg_ver, character(1)),
    pkgs
  )

  bg_types <- if ("background_type" %in% names(results)) {
    unique(stats::na.omit(as.character(results$background_type)))
  } else {
    character(0)
  }
  univ_sizes <- if ("universe_size" %in% names(results)) {
    unique(stats::na.omit(as.integer(results$universe_size)))
  } else {
    integer(0)
  }
  diag_summary <- list(
    n_rows_total = as.integer(nrow(results)),
    n_rows_evaluable = as.integer(sum(!(results$skipped %in% TRUE), na.rm = TRUE)),
    n_rows_skipped = as.integer(sum(results$skipped %in% TRUE, na.rm = TRUE)),
    null_background_types = as.list(bg_types),
    null_universe_sizes = as.list(univ_sizes),
    has_small_positive_class = isTRUE("small_positive_class" %in% names(results) && any(results$small_positive_class %in% TRUE, na.rm = TRUE)),
    has_unreliable_pr_ci = isTRUE("ci_unreliable" %in% names(results) && any(results$ci_unreliable %in% TRUE, na.rm = TRUE)),
    has_survival = isTRUE(!is.null(report_model$tier2$survival) && is.data.frame(report_model$tier2$survival) && nrow(report_model$tier2$survival) > 0L),
    has_beyond_infiltration = isTRUE(!is.null(report_model$tier2$beyond_infiltration) && is.list(report_model$tier2$beyond_infiltration) &&
                                       is.data.frame(report_model$tier2$beyond_infiltration$per_cohort) &&
                                       nrow(report_model$tier2$beyond_infiltration$per_cohort) > 0L)
  )

  list(
    sigbench_version = .pkg_ver("SigBenchR"),
    timestamp_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    r_version = R.version.string,
    package_versions = as.list(pkg_vers),
    run_params = run_params,
    rules_meta = rules_meta,
    cohorts_fingerprint = cohorts_fp,
    signatures_fingerprint = sigs_fp,
    diagnostics_summary = diag_summary,
    reproducibility_hash = repro_hash
  )
}

#' Export a reproducibility manifest to JSON or YAML
#'
#' @param manifest Output of [build_reproducibility_manifest()].
#' @param path Output file path.
#' @param format `"json"` (default) or `"yaml"`.
#'
#' @return Normalized output file path (invisibly).
#' @export
export_manifest <- function(manifest, path, format = c("json", "yaml")) {
  format <- match.arg(format)
  if (!is.list(manifest) || is.null(manifest$reproducibility_hash)) {
    stop("`manifest` must be a valid reproducibility manifest list.", call. = FALSE)
  }
  if (missing(path) || !nzchar(as.character(path)[1])) {
    stop("`path` is required.", call. = FALSE)
  }
  path <- normalizePath(as.character(path)[1], mustWork = FALSE)
  out_dir <- dirname(path)
  if (!dir.exists(out_dir)) {
    ok_dir <- dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    if (!isTRUE(ok_dir) && !dir.exists(out_dir)) {
      stop("Cannot create output directory for manifest: ", out_dir, call. = FALSE)
    }
  }
  if (!identical(as.integer(file.access(out_dir, 2)), 0L)) {
    stop("Output directory is not writable: ", out_dir, call. = FALSE)
  }

  if (identical(format, "json")) {
    jsonlite::write_json(manifest, path, auto_unbox = TRUE, pretty = TRUE, null = "null")
  } else {
    yaml::write_yaml(manifest, path)
  }
  invisible(path)
}

#' @noRd
.pkg_ver <- function(pkg) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    as.character(utils::packageVersion(pkg))
  } else {
    "not installed"
  }
}

#' @noRd
.md5_string <- function(txt) {
  raw_bytes <- charToRaw(enc2utf8(as.character(txt)[1]))
  if ("bytes" %in% names(formals(tools::md5sum))) {
    return(unname(tools::md5sum(bytes = raw_bytes)))
  }
  tmp <- tempfile(fileext = ".txt")
  on.exit(unlink(tmp), add = TRUE)
  writeBin(raw_bytes, tmp)
  unname(tools::md5sum(tmp))
}

#' @noRd
.fingerprint_cohorts <- function(results, datasets = NULL) {
  out <- list()
  if (is.list(datasets) && length(datasets) > 0L) {
    for (nm in sort(names(datasets))) {
      ds <- datasets[[nm]]
      sids <- colnames(ds$expr)
      gids <- rownames(ds$expr)
      out[[nm]] <- list(
        n_samples = as.integer(ncol(ds$expr)),
        n_genes = as.integer(nrow(ds$expr)),
        sample_ids_md5 = .md5_string(paste(sids, collapse = ",")),
        genes_md5 = .md5_string(paste(gids, collapse = ",")),
        normalization = as.character(ds$info$normalization)[1],
        scale = as.character(ds$info$scale)[1]
      )
    }
    return(out)
  }
  if (nrow(results) == 0L || !"dataset" %in% names(results)) {
    return(out)
  }
  ds_names <- sort(unique(as.character(results$dataset)))
  for (nm in ds_names) {
    sub <- results[results$dataset == nm, , drop = FALSE]
    n_val <- if ("n" %in% names(sub) && any(is.finite(sub$n))) {
      as.integer(sub$n[which(is.finite(sub$n))[1]])
    } else {
      NA_integer_
    }
    norm_val <- if ("normalization" %in% names(sub)) as.character(sub$normalization[1]) else NA_character_
    scale_val <- if ("scale" %in% names(sub)) as.character(sub$scale[1]) else NA_character_
    out[[nm]] <- list(
      n_samples = n_val,
      normalization = norm_val,
      scale = scale_val
    )
  }
  out
}

#' @noRd
.fingerprint_signatures <- function(results, signatures = NULL) {
  out <- list()
  if (is.list(signatures) && length(signatures) > 0L) {
    for (nm in sort(names(signatures))) {
      sig <- .as_signature(signatures[[nm]], name = nm)
      genes <- if (!is.null(sig$weights) && is.numeric(sig$weights) && !is.null(names(sig$weights))) {
        ord <- order(toupper(names(sig$weights)))
        w_ord <- sig$weights[ord]
        paste(toupper(names(w_ord)), sprintf("%.8g", as.numeric(w_ord)), sep = ":", collapse = ",")
      } else if (!is.null(sig$genes)) {
        if (is.numeric(sig$genes)) paste(names(sig$genes), sig$genes, sep = ":", collapse = ",") else paste(sig$genes, collapse = ",")
      } else {
        paste(c(paste0("+", sig$genes_up), paste0("-", sig$genes_down)), collapse = ",")
      }
      out[[nm]] <- list(
        direction = if (!is.null(sig$direction)) as.character(sig$direction)[1] else NA_character_,
        weighted = !is.null(sig$weights),
        genes_md5 = .md5_string(genes)
      )
    }
    return(out)
  }
  if (nrow(results) == 0L || !"signature" %in% names(results)) {
    return(out)
  }
  sig_names <- sort(unique(as.character(results$signature)))
  for (nm in sig_names) {
    sub <- results[results$signature == nm, , drop = FALSE]
    fg <- if ("found_genes" %in% names(sub)) paste(sort(unique(sub$found_genes)), collapse = "|") else ""
    mg <- if ("missing_genes" %in% names(sub)) paste(sort(unique(sub$missing_genes)), collapse = "|") else ""
    out[[nm]] <- list(
      n_total = if ("n_total" %in% names(sub)) as.integer(sub$n_total[1]) else NA_integer_,
      signed = if ("signed" %in% names(sub)) isTRUE(sub$signed[1]) else NA,
      genes_md5 = .md5_string(paste(fg, mg, sep = "||"))
    )
  }
  out
}
