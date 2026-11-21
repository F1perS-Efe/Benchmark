#' Benchmark sonuclarini HTML rapor olarak disa aktar
#'
#' `inst/rmd/benchmark_report.Rmd` sablonunu `rmarkdown::render()` ile
#' isler. `rmarkdown` Suggests bagimliligidır; yoksa anlasilir hata verir.
#'
#' @param results `run_benchmark()` cikti tablosu.
#' @param summary_tbl `summarize_auc()` ozet tablosu.
#' @param sensitivity_tbl `check_sensitivity()` tablosu.
#' @param params_list Analiz parametreleri listesi. Beklenen alanlar:
#'   `methods`, `direction`, `coverage_threshold`, `sensitivity_threshold`,
#'   `n_boot`, `n_perm`, `seed`; istege bagli `heatmap_method`,
#'   `forest_signature`.
#' @param output_file Cikti HTML dosya yolu.
#'
#' @return `output_file` yolu (gorunmez).
#' @export
export_report <- function(results,
                          summary_tbl,
                          sensitivity_tbl,
                          params_list,
                          output_file) {
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
    heatmap_method = NULL,
    forest_signature = NULL
  )
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

  rmarkdown::render(
    input = rmd,
    output_file = basename(output_file),
    output_dir = out_dir,
    params = list(
      results = results,
      summary_tbl = summary_tbl,
      sensitivity_tbl = sensitivity_tbl,
      params_list = params_list
    ),
    envir = new.env(parent = globalenv()),
    quiet = TRUE
  )

  invisible(output_file)
}
