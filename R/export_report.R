#' Benchmark sonuclarini HTML rapor olarak disa aktar
#'
#' `inst/rmd/benchmark_report.Rmd` sablonunu `rmarkdown::render()` ile
#' isler. Rapor "research use only / not a clinical decision support tool"
#' uyarisini icerir. `rmarkdown` Suggests bagimliligidir; yoksa anlasilir
#' hata verir.
#'
#' @param results `run_benchmark()` cikti tablosu.
#' @param summary_tbl `summarize_auc()` ozet tablosu.
#' @param sensitivity_tbl `check_sensitivity()` tablosu.
#' @param params_list Analiz parametreleri listesi. Beklenen alanlar:
#'   `methods`, `direction`, `coverage_threshold`, `sensitivity_threshold`,
#'   `n_boot`, `n_perm`, `seed`, `compute_null`; istege bagli
#'   `heatmap_method`, `forest_signature`.
#' @param output_file Cikti HTML dosya yolu.
#' @param surv_res Opsiyonel `run_survival_benchmark()` tablosu. `NULL`
#'   (varsayilan) iken eger `params_list$datasets` ve `params_list$signatures`
#'   verilmisse ve kohortlarda sagkalim sutunlari varsa otomatik hesaplanir.
#' @param beyond_inf Opsiyonel `test_beyond_infiltration()` ciktisi (veya imza
#'   bazli liste).
#' @param report_model Optional pre-built [build_report_model()] object.
#' @param manifest Optional pre-built [build_reproducibility_manifest()] object.
#'
#' @return `output_file` yolu (gorunmez).
#' @export
export_report <- function(results,
                          summary_tbl,
                          sensitivity_tbl,
                          params_list,
                          output_file,
                          surv_res = NULL,
                          report_model = NULL,
                          manifest = NULL,
                          beyond_inf = NULL) {
  if (!requireNamespace("rmarkdown", quietly = TRUE)) {
    stop(
      "Install 'rmarkdown' to render HTML reports: install.packages(\"rmarkdown\")",
      call. = FALSE
    )
  }
  if (!is.data.frame(results)) {
    stop("`results` must be a data.frame.", call. = FALSE)
  }
  if (!is.data.frame(summary_tbl)) {
    stop("`summary_tbl` must be a data.frame.", call. = FALSE)
  }
  if (!is.data.frame(sensitivity_tbl)) {
    stop("`sensitivity_tbl` must be a data.frame.", call. = FALSE)
  }
  if (!is.list(params_list)) {
    stop("`params_list` must be a list.", call. = FALSE)
  }
  if (missing(output_file) || !nzchar(as.character(output_file)[1])) {
    stop("`output_file` is required.", call. = FALSE)
  }
  output_file <- normalizePath(as.character(output_file)[1], mustWork = FALSE)

  rmd <- system.file("rmd", "benchmark_report.Rmd", package = "SigBenchR")
  if (!nzchar(rmd) || !file.exists(rmd)) {
    stop("Report template not found (inst/rmd/benchmark_report.Rmd).", call. = FALSE)
  }

  defaults <- list(
    methods = character(0),
    direction = "<",
    coverage_threshold = 0.70,
    sensitivity_threshold = 0.10,
    n_boot = NA_integer_,
    n_perm = NA_integer_,
    seed = NA_integer_,
    compute_null = FALSE,
    heatmap_method = NULL,
    forest_signature = NULL
  )
  if (is.null(params_list$compute_null)) {
    # Infer from results when callers omit the flag
    if (is.data.frame(results) && "null_p" %in% names(results) &&
        any(is.finite(results$null_p), na.rm = TRUE)) {
      params_list$compute_null <- TRUE
    }
  }
  for (nm in names(defaults)) {
    if (is.null(params_list[[nm]])) {
      params_list[[nm]] <- defaults[[nm]]
    }
  }
  if (is.null(params_list$heatmap_method) || !nzchar(as.character(params_list$heatmap_method)[1])) {
    if (length(params_list$methods) > 0L) {
      params_list$heatmap_method <- params_list$methods[1]
    } else if ("method" %in% names(results) && nrow(results) > 0L) {
      params_list$heatmap_method <- results$method[1]
    } else {
      params_list$heatmap_method <- "mean"
    }
  }
  if (is.null(params_list$forest_signature) || !nzchar(as.character(params_list$forest_signature)[1])) {
    if ("signature" %in% names(results) && nrow(results) > 0L) {
      params_list$forest_signature <- results$signature[1]
    } else {
      params_list$forest_signature <- NA_character_
    }
  }

  out_dir <- dirname(output_file)
  if (!dir.exists(out_dir)) {
    ok_dir <- dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    if (!isTRUE(ok_dir) && !dir.exists(out_dir)) {
      stop("Cannot create output directory for report: ", out_dir, call. = FALSE)
    }
  }
  if (!identical(as.integer(file.access(out_dir, 2)), 0L)) {
    stop("Output directory is not writable: ", out_dir, call. = FALSE)
  }

  if (!is.null(surv_res) && !is.data.frame(surv_res)) {
    stop("`surv_res` must be NULL or a data.frame.", call. = FALSE)
  }

  ref_axes <- if (!is.null(params_list$bio_context) && is.data.frame(params_list$bio_context)) {
    params_list$bio_context
  } else if (!is.null(params_list$datasets) && !is.null(params_list$signatures)) {
    tryCatch(
      correlate_reference_axes(
        datasets = params_list$datasets,
        signatures = params_list$signatures,
        method = params_list$heatmap_method,
        sd_handling = if (!is.null(params_list$sd_handling)) params_list$sd_handling else "exclude"
      ),
      error = function(e) NULL
    )
  } else {
    NULL
  }

  if (is.null(surv_res) && !is.null(report_model$tier2$survival) &&
      is.data.frame(report_model$tier2$survival)) {
    surv_res <- report_model$tier2$survival
  } else if (is.null(surv_res) && !is.null(params_list$surv_res) &&
             is.data.frame(params_list$surv_res)) {
    surv_res <- params_list$surv_res
  } else if (is.null(surv_res) && !is.null(params_list$datasets) && !is.null(params_list$signatures)) {
    surv_res <- tryCatch(
      run_survival_benchmark(
        datasets = params_list$datasets,
        signatures = params_list$signatures,
        method = params_list$heatmap_method
      ),
      error = function(e) NULL
    )
  }

  if (is.null(beyond_inf) && !is.null(report_model$tier2$beyond_infiltration)) {
    beyond_inf <- report_model$tier2$beyond_infiltration
  } else if (is.null(beyond_inf) && !is.null(params_list$beyond_inf)) {
    beyond_inf <- params_list$beyond_inf
  } else if (is.null(beyond_inf) && !is.null(params_list$datasets) && !is.null(params_list$signatures)) {
    beyond_inf <- tryCatch(
      {
        bi_list <- list()
        for (sn in names(params_list$signatures)) {
          bi_list[[sn]] <- test_beyond_infiltration(
            datasets = params_list$datasets,
            signature = params_list$signatures[[sn]],
            infiltration_axis = "CD8_T_Cell",
            method = params_list$heatmap_method,
            sd_handling = if (!is.null(params_list$sd_handling)) params_list$sd_handling else "exclude",
            direction = if (!is.null(params_list$direction)) params_list$direction else "<"
          )
        }
        .normalize_beyond_inf_for_report(bi_list)
      },
      error = function(e) NULL
    )
  }

  if (is.null(report_model)) {
    report_model <- build_report_model(
      results = results,
      summary_tbl = if (nrow(summary_tbl) > 0L) summary_tbl else NULL,
      sensitivity_tbl = sensitivity_tbl,
      params_list = params_list,
      ref_axes = ref_axes,
      surv_res = surv_res,
      beyond_inf = beyond_inf
    )
  }

  if (is.null(manifest)) {
    manifest <- build_reproducibility_manifest(
      results = results,
      summary_tbl = if (nrow(summary_tbl) > 0L) summary_tbl else NULL,
      report_model = report_model,
      params_list = params_list,
      datasets = params_list$datasets,
      signatures = params_list$signatures
    )
  }

  rmarkdown::render(
    input = rmd,
    output_file = basename(output_file),
    output_dir = out_dir,
    params = list(
      results = results,
      summary_tbl = summary_tbl,
      sensitivity_tbl = sensitivity_tbl,
      params_list = params_list,
      surv_res = surv_res,
      beyond_inf = report_model$tier2$beyond_infiltration,
      report_model = report_model,
      manifest = manifest
    ),
    envir = new.env(parent = asNamespace("SigBenchR")),
    quiet = TRUE
  )

  invisible(output_file)
}
