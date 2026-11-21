#' Benchmark sonuclarini HTML rapor olarak disa aktar
#'
#' `inst/rmd/benchmark_report.Rmd` sablonunu `rmarkdown::render()` ile
#' isler. `rmarkdown` Suggests bagimliligidir; yoksa anlasilir hata verir.
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
#'   (varsayilan) iken sagkalim bolumu render edilmez (geriye donuk uyum).
#'
#' @return `output_file` yolu (gorunmez).
#' @export
export_report <- function(results,
                          summary_tbl,
                          sensitivity_tbl,
                          params_list,
                          output_file,
                          surv_res = NULL) {
  if (!requireNamespace("rmarkdown", quietly = TRUE)) {
    stop(
      "HTML rapor icin 'rmarkdown' paketini kurun: install.packages(\"rmarkdown\")",
      call. = FALSE
    )
  }
  if (!is.data.frame(results)) {
    stop("`results` bir data.frame olmalidir.", call. = FALSE)
  }
  if (!is.data.frame(summary_tbl)) {
    stop("`summary_tbl` bir data.frame olmalidir.", call. = FALSE)
  }
  if (!is.data.frame(sensitivity_tbl)) {
    stop("`sensitivity_tbl` bir data.frame olmalidir.", call. = FALSE)
  }
  if (!is.list(params_list)) {
    stop("`params_list` bir liste olmalidir.", call. = FALSE)
  }
  if (missing(output_file) || !nzchar(as.character(output_file)[1])) {
    stop("`output_file` zorunludur.", call. = FALSE)
  }
  output_file <- normalizePath(as.character(output_file)[1], mustWork = FALSE)

  rmd <- system.file("rmd", "benchmark_report.Rmd", package = "SigBenchR")
  if (!nzchar(rmd) || !file.exists(rmd)) {
    stop("Rapor sablonu bulunamadi (inst/rmd/benchmark_report.Rmd).", call. = FALSE)
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
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }

  if (!is.null(surv_res) && !is.data.frame(surv_res)) {
    stop("`surv_res` NULL veya data.frame olmalidir.", call. = FALSE)
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
      surv_res = surv_res
    ),
    envir = new.env(parent = globalenv()),
    quiet = TRUE
  )

  invisible(output_file)
}
