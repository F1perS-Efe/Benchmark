#!/usr/bin/env Rscript
# One-command HTML demo report for SigBenchR.
# Usage: Rscript inst/scripts/render_demo_report.R [output_path]
# Default output: sigbenchr_demo_report.html

args <- commandArgs(trailingOnly = TRUE)
output_file <- if (length(args) >= 1L && nzchar(as.character(args[[1]])[1])) {
  as.character(args[[1]])[1]
} else {
  "sigbenchr_demo_report.html"
}

.find_pkg_root <- function() {
  ca <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", ca, value = TRUE)
  if (length(file_arg) >= 1L) {
    script <- sub("^--file=", "", file_arg[[1]])
    script <- normalizePath(script, mustWork = FALSE)
    root <- normalizePath(file.path(dirname(script), "..", ".."), mustWork = FALSE)
    if (file.exists(file.path(root, "DESCRIPTION"))) {
      return(root)
    }
  }
  if (file.exists("DESCRIPTION")) {
    return(normalizePath(".", mustWork = FALSE))
  }
  NULL
}

.root <- .find_pkg_root()
if (!is.null(.root) && requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(.root, quiet = TRUE)
} else if (!is.null(.root) && requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(.root, quiet = TRUE)
} else if (requireNamespace("SigBenchR", quietly = TRUE)) {
  library(SigBenchR)
} else {
  stop(
    "SigBenchR is not installed. From the package root run ",
    "devtools::install('.') or Rscript inst/scripts/render_demo_report.R ",
    "with pkgload/devtools available.",
    call. = FALSE
  )
}

if (!requireNamespace("rmarkdown", quietly = TRUE)) {
  stop(
    "Package 'rmarkdown' is required. Install with install.packages('rmarkdown').",
    call. = FALSE
  )
}

data("datasets_ici", package = "SigBenchR")
sigs <- get_builtin_signatures(include_controls = TRUE)

# Demo keeps methods compact; null permutations follow the release recipe.
methods_use <- c("mean", "zscore")
res <- run_benchmark(
  signatures = sigs,
  datasets = datasets_ici,
  methods = methods_use,
  compute_pr = FALSE,
  compute_null = TRUE,
  n_perm = 200L,
  seed = 42L
)

summary_tbl <- suppressWarnings(summarize_auc(res, method = "mean"))
sensitivity_tbl <- check_sensitivity(res, threshold = 0.10)

params_list <- list(
  methods = methods_use,
  direction = "<",
  coverage_threshold = 0.70,
  sensitivity_threshold = 0.10,
  n_boot = 0L,
  n_perm = 200L,
  seed = 42L,
  compute_null = TRUE,
  sd_handling = "exclude",
  heatmap_method = "mean",
  forest_signature = if ("IFNG6" %in% res$signature) {
    "IFNG6"
  } else if (nrow(res) > 0L) {
    as.character(res$signature[1])
  } else {
    NA_character_
  }
)

output_file <- normalizePath(output_file, mustWork = FALSE)
path <- export_report(
  results = res,
  summary_tbl = summary_tbl,
  sensitivity_tbl = sensitivity_tbl,
  params_list = params_list,
  output_file = output_file
)

message("Demo HTML report written to: ", path)
