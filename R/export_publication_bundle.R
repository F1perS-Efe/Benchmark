#' Export a publication-ready bundle (.zip or directory)
#'
#' Writes publication-resolution figures (300 DPI PNG + vector PDF), tidy CSV
#' tables (`table_meta_summary.csv`, `table_cohort_metrics.csv`,
#' `table_tier1_verdicts.csv`, `table_tier1_alerts.csv`,
#' `table_beyond_infiltration_per_cohort.csv`,
#' `table_beyond_infiltration_summary.csv`, `table_survival_benchmark.csv`),
#' machine-readable `analysis_objects.rds`, `manifest.json`, and
#' `figure_legends.txt` containing methods summary, limitations, and
#' `reproducibility_hash`.
#'
#' @param results `run_benchmark()` output `data.frame`.
#' @param output_path Destination `.zip` file path (when `zip = TRUE`) or
#'   directory path (when `zip = FALSE`).
#' @param summary_tbl Optional `summarize_auc()` `data.frame`.
#' @param report_model Optional `build_report_model()` output.
#' @param manifest Optional `build_reproducibility_manifest()` output.
#' @param params_list Optional list of run parameters.
#' @param method Scoring method for figures (`NULL` uses `params_list$heatmap_method`
#'   or the first method in `results`).
#' @param signature Signature name for the forest plot (`NULL` uses
#'   `params_list$forest_signature`, best independent signature, or first row).
#' @param coverage_threshold Coverage threshold (0-1) for heatmap shading.
#' @param zip If `TRUE` (default), archive files into a `.zip` file at
#'   `output_path`; if `FALSE`, write directly into `output_path` directory.
#' @param surv_res Optional `run_survival_benchmark()` `data.frame`.
#' @param beyond_inf Optional `test_beyond_infiltration()` result (single list
#'   or named list across signatures).
#'
#' @return Normalized `output_path` (invisibly).
#' @export
export_publication_bundle <- function(results,
                                      output_path,
                                      summary_tbl = NULL,
                                      report_model = NULL,
                                      manifest = NULL,
                                      params_list = list(),
                                      method = NULL,
                                      signature = NULL,
                                      coverage_threshold = 0.70,
                                      zip = TRUE,
                                      surv_res = NULL,
                                      beyond_inf = NULL) {
  if (!is.data.frame(results) || nrow(results) == 0L) {
    stop("`results` cannot be empty.", call. = FALSE)
  }
  if (missing(output_path) || !nzchar(as.character(output_path)[1])) {
    stop("`output_path` is required.", call. = FALSE)
  }
  output_path <- normalizePath(as.character(output_path)[1], mustWork = FALSE)

  if (is.null(params_list)) {
    params_list <- list()
  }
  if (is.null(method) || !nzchar(as.character(method)[1])) {
    method <- if (!is.null(params_list$heatmap_method) && nzchar(as.character(params_list$heatmap_method)[1])) {
      as.character(params_list$heatmap_method)[1]
    } else {
      as.character(results$method[1])
    }
  }
  if (is.null(params_list$heatmap_method)) {
    params_list$heatmap_method <- method
  }
  if (is.null(params_list$coverage_threshold)) {
    params_list$coverage_threshold <- coverage_threshold
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
        method = method
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
            method = method,
            sd_handling = if (!is.null(params_list$sd_handling)) params_list$sd_handling else "exclude",
            direction = if (!is.null(params_list$direction)) params_list$direction else "<"
          )
        }
        .normalize_beyond_inf_for_report(bi_list)
      },
      error = function(e) NULL
    )
  } else if (!is.null(beyond_inf)) {
    beyond_inf <- .normalize_beyond_inf_for_report(beyond_inf)
  }

  if (is.null(report_model)) {
    report_model <- build_report_model(
      results = results,
      summary_tbl = summary_tbl,
      params_list = params_list,
      surv_res = surv_res,
      beyond_inf = beyond_inf
    )
  } else {
    if (is.null(report_model$tier2$survival) && is.data.frame(surv_res)) {
      report_model$tier2$survival <- surv_res
    }
    if (is.null(report_model$tier2$beyond_infiltration) && !is.null(beyond_inf)) {
      report_model$tier2$beyond_infiltration <- beyond_inf
    }
  }
  if (is.null(manifest)) {
    manifest <- build_reproducibility_manifest(
      results = results,
      summary_tbl = summary_tbl,
      report_model = report_model,
      params_list = params_list,
      datasets = params_list$datasets,
      signatures = params_list$signatures
    )
  }

  if (is.null(signature) || !nzchar(as.character(signature)[1])) {
    signature <- if (!is.null(params_list$forest_signature) && nzchar(as.character(params_list$forest_signature)[1])) {
      as.character(params_list$forest_signature)[1]
    } else if (!is.null(report_model$tier1$best_signature) && nrow(report_model$tier1$best_signature) > 0L) {
      as.character(report_model$tier1$best_signature$signature[1])
    } else {
      as.character(results$signature[1])
    }
  }

  dest_dir <- if (isTRUE(zip)) dirname(output_path) else output_path
  if (!dir.exists(dest_dir)) {
    ok_dir <- dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
    if (!isTRUE(ok_dir) && !dir.exists(dest_dir)) {
      stop("Cannot create output directory for publication bundle: ", dest_dir, call. = FALSE)
    }
  }
  if (!identical(as.integer(file.access(dest_dir, 2)), 0L)) {
    stop("Output directory is not writable: ", dest_dir, call. = FALSE)
  }

  work_dir <- if (isTRUE(zip)) {
    d <- tempfile("sigbenchr_pub_")
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
    on.exit(unlink(d, recursive = TRUE, force = TRUE), add = TRUE)
    d
  } else {
    output_path
  }

  # 1) High-resolution figures (PNG 300 DPI + vector PDF)
  has_hm <- any(results$method == method & !(results$skipped %in% TRUE))
  if (isTRUE(has_hm)) {
    p_hm <- plot_auc_heatmap(results, method = method, coverage_threshold = coverage_threshold)
    ggplot2::ggsave(
      filename = file.path(work_dir, "fig_auc_heatmap.png"),
      plot = p_hm, width = 8, height = 5, dpi = 300
    )
    ggplot2::ggsave(
      filename = file.path(work_dir, "fig_auc_heatmap.pdf"),
      plot = p_hm, width = 8, height = 5
    )
  }

  has_fr <- any(results$signature == signature & results$method == method & !(results$skipped %in% TRUE))
  if (isTRUE(has_fr)) {
    p_fr <- plot_auc_forest(
      results,
      signature = signature,
      method = method,
      summary = report_model$tier2$summary_tbl
    )
    ggplot2::ggsave(
      filename = file.path(work_dir, "fig_auc_forest.png"),
      plot = p_fr, width = 8, height = 4.5, dpi = 300
    )
    ggplot2::ggsave(
      filename = file.path(work_dir, "fig_auc_forest.pdf"),
      plot = p_fr, width = 8, height = 4.5
    )
  }

  # Optional ROC curves for the first evaluable cohort when datasets & signatures are available
  if (is.list(params_list$datasets) && length(params_list$datasets) > 0L &&
      is.list(params_list$signatures) && length(params_list$signatures) > 0L) {
    roc_ds <- names(params_list$datasets)[1]
    p_roc <- tryCatch(
      plot_roc_curves(
        datasets = params_list$datasets,
        signatures = params_list$signatures,
        dataset_name = roc_ds,
        method = method,
        direction = if (!is.null(params_list$direction)) params_list$direction else "<",
        sd_handling = if (!is.null(params_list$sd_handling)) params_list$sd_handling else "exclude"
      ),
      error = function(e) NULL
    )
    if (!is.null(p_roc)) {
      ggplot2::ggsave(
        filename = file.path(work_dir, "fig_roc_curves.png"),
        plot = p_roc, width = 7.5, height = 5, dpi = 300
      )
      ggplot2::ggsave(
        filename = file.path(work_dir, "fig_roc_curves.pdf"),
        plot = p_roc, width = 7.5, height = 5
      )
    }
  }

  # 2) Publication CSV tables (full numeric precision preserved; NA written as NA)
  meta_tbl <- report_model$tier2$summary_formatted
  if (!is.data.frame(meta_tbl)) {
    meta_tbl <- data.frame()
  }
  utils::write.csv(meta_tbl, file.path(work_dir, "table_meta_summary.csv"), row.names = FALSE)

  coh_cols <- intersect(
    c(
      "dataset", "signature", "method", "n", "n_pos", "n_neg",
      "auc", "auc_ci_low", "auc_ci_high", "auc_var", "var_zero", "small_cohort",
      "pr_auc", "pr_baseline", "pr_ci_low", "pr_ci_high",
      "tie_fraction", "small_positive_class", "n_boot_valid", "ci_unreliable",
      "coverage_pct", "coverage_status", "low_coverage", "n_found", "n_total",
      "null_p", "null_p_adj_bh", "null_p_adj_bh_global",
      "background_type", "universe_size", "power_auc70",
      "independent", "normalization", "scale", "treatment", "endpoint_type",
      "skipped", "skip_reason"
    ),
    names(results)
  )
  utils::write.csv(
    results[, coh_cols, drop = FALSE],
    file.path(work_dir, "table_cohort_metrics.csv"),
    row.names = FALSE
  )

  verdict_tbl <- report_model$tier1$signature_cards
  if (!is.data.frame(verdict_tbl)) {
    verdict_tbl <- data.frame()
  }
  utils::write.csv(verdict_tbl, file.path(work_dir, "table_tier1_verdicts.csv"), row.names = FALSE)

  alerts_tbl <- report_model$tier1$alerts
  if (!is.data.frame(alerts_tbl)) {
    alerts_tbl <- data.frame(
      rule_id = character(), severity = character(), title = character(),
      message = character(), metric_value = character(), stringsAsFactors = FALSE
    )
  }
  utils::write.csv(alerts_tbl, file.path(work_dir, "table_tier1_alerts.csv"), row.names = FALSE)

  bi_obj <- report_model$tier2$beyond_infiltration
  bi_pc <- if (is.list(bi_obj) && is.data.frame(bi_obj$per_cohort)) bi_obj$per_cohort else data.frame()
  bi_sm <- if (is.list(bi_obj) && is.data.frame(bi_obj$summary)) bi_obj$summary else data.frame()
  utils::write.csv(bi_pc, file.path(work_dir, "table_beyond_infiltration_per_cohort.csv"), row.names = FALSE)
  utils::write.csv(bi_sm, file.path(work_dir, "table_beyond_infiltration_summary.csv"), row.names = FALSE)

  surv_out <- if (is.data.frame(report_model$tier2$survival)) {
    report_model$tier2$survival
  } else {
    .survival_empty_schema(with_signature = TRUE)
  }
  utils::write.csv(surv_out, file.path(work_dir, "table_survival_benchmark.csv"), row.names = FALSE)

  # Machine-readable full precision R object bundle
  saveRDS(
    list(
      results = results,
      summary_tbl = report_model$tier2$summary_tbl,
      report_model = report_model,
      manifest = manifest,
      beyond_infiltration = bi_obj,
      survival = surv_out
    ),
    file = file.path(work_dir, "analysis_objects.rds")
  )

  # 3) Manifest JSON & figure legends / methods note
  export_manifest(manifest, file.path(work_dir, "manifest.json"), format = "json")

  seed_txt <- if (!is.null(manifest$run_params$seed)) as.character(manifest$run_params$seed) else "NA"
  sd_txt <- if (!is.null(manifest$run_params$sd_handling)) as.character(manifest$run_params$sd_handling) else "exclude"
  bg_types <- if ("background_type" %in% names(results)) {
    unique(stats::na.omit(as.character(results$background_type)))
  } else {
    character(0)
  }
  bg_txt <- if (length(bg_types) > 0L) paste(bg_types, collapse = ", ") else "Not run (compute_null = FALSE)"

  legend_lines <- c(
    "SigBenchR Publication Bundle - Figure Legends & Methods Summary",
    "================================================================",
    paste0("SigBenchR version: ", manifest$sigbench_version),
    paste0("R version: ", manifest$r_version),
    paste0("Reproducibility hash (MD5): ", manifest$reproducibility_hash),
    paste0("Rules version / MD5: ", manifest$rules_meta$version, " / ", manifest$rules_meta$md5),
    paste0("Timestamp (UTC): ", manifest$timestamp_utc),
    paste0("Null permutation background universe: ", bg_txt),
    "",
    "Methods Summary:",
    sprintf(
      paste0(
        "Gene signature discrimination was evaluated per cohort without cross-platform ",
        "expression pooling (scoring method: %s; SD handling: %s; gene coverage threshold: %.0f%%; ",
        "RNG seed: %s). Cohort-level ROC-AUC and 95%% CIs were computed via DeLong's method. ",
        "PR-AUC and stratified bootstrap 95%% CIs flag small positive/negative classes (n_pos < 5 or n_neg < 5) ",
        "and tie/bootstrap instability (ci_unreliable). Random-gene null permutation tests report raw null_p, ",
        "within-method BH-adjusted FDR (null_p_adj_bh), and across-method global BH-adjusted FDR (null_p_adj_bh_global). ",
        "Beyond-infiltration adjusted logistic models report fold-inside paired CV delta-AUC (with sd_delta_auc and 95%% CI), ",
        "BH-adjusted LRT p-values (lrt_p_adj_bh), EPV, VIF, analysis_status, and directional Stouffer Z cross-cohort synthesis. ",
        "Survival benchmarks adjust continuous Cox PH (p_adj_bh) and median-dichotomized Log-rank (logrank_p_adj_bh) as separate ",
        "hypothesis families while leaving Schoenfeld PH assumption checks (ph_p) unadjusted as a model diagnostic. ",
        "Cross-cohort AUC synthesis used logit-transformed AUC random-effects meta-analysis (REML) ",
        "with truncated Hartung-Knapp-Sidik-Jonkman (HKSJ) 95%% intervals and logit-scale I2."
      ),
      method, sd_txt, 100 * as.numeric(coverage_threshold), seed_txt
    ),
    "",
    "Scientific Limitations & Scope:",
    "- Research exploration tool only; not a calibrated clinical decision support device.",
    "- Weighted signatures (e.g., T_Cell_Inflamed_GEP continuous regression coefficients) are evaluated on cohort RNA-seq scales and do not claim analytical equivalence to proprietary clinical assay platforms (e.g., NanoString nCounter housekeeping-normalized GEP).",
    "- Evaluations with low EPV (< 5), high VIF (>= 5), unreliable bootstrap CIs (ci_unreliable = TRUE), or small cohort counts (k <= 4) are exploratory.",
    "",
    "Bundle Contents & Figure Legends:",
    sprintf(
      "- fig_auc_heatmap.{png,pdf}: Signature x cohort ROC-AUC heatmap (%s). Cells below %.0f%% gene coverage are shaded grey.",
      method, 100 * as.numeric(coverage_threshold)
    ),
    sprintf(
      "- fig_auc_forest.{png,pdf}: Cohort-level ROC-AUC forest plot for %s (%s) with DeLong 95%% CIs and logit-AUC REML + HKSJ summary.",
      signature, method
    ),
    "- fig_roc_curves.{png,pdf}: Cohort ROC curves with DeLong 95% CIs and sample counts (when cohort expression data are provided).",
    "- table_cohort_metrics.csv: Full-precision cohort x signature x method ROC-AUC, PR-AUC (with tie_fraction, small_positive_class, n_boot_valid, ci_unreliable), and null permutation metrics (null_p, null_p_adj_bh, null_p_adj_bh_global, background_type, universe_size).",
    "- table_meta_summary.csv: Logit-AUC REML + HKSJ meta-analysis summary.",
    "- table_tier1_verdicts.csv & table_tier1_alerts.csv: Rule-based signature evidence tiers and statistical diagnostic alerts.",
    "- table_beyond_infiltration_per_cohort.csv & table_beyond_infiltration_summary.csv: Adjusted logistic regression, paired CV delta-AUC CIs, EPV, VIF, BH-adjusted LRT, and directional Stouffer Z synthesis.",
    "- table_survival_benchmark.csv: Cox HR (cox_p, p_adj_bh), Log-rank (logrank_p, logrank_p_adj_bh), and unadjusted Schoenfeld PH diagnostics (ph_p, ph_violated).",
    "- analysis_objects.rds & manifest.json: Machine-readable R objects and deterministic reproducibility manifest."
  )
  writeLines(legend_lines, con = file.path(work_dir, "figure_legends.txt"))

  if (isTRUE(zip)) {
    out_dir <- dirname(output_path)
    if (!dir.exists(out_dir)) {
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    }
    if (file.exists(output_path)) {
      unlink(output_path)
    }
    files_rel <- list.files(work_dir, recursive = TRUE)
    if (requireNamespace("zip", quietly = TRUE)) {
      withr::with_dir(work_dir, {
        zip::zipr(zipfile = output_path, files = files_rel)
      })
    } else {
      zip_cmd <- Sys.getenv("R_ZIPCMD", unset = "zip")
      if (!nzchar(Sys.which(zip_cmd))) {
        stop(
          "Creating a .zip publication bundle requires either the 'zip' R package ",
          "(recommended across Windows/macOS/Linux: install.packages(\"zip\")) ",
          "or a system 'zip' executable on PATH, or pass `zip = FALSE` to write ",
          "a directory bundle directly.",
          call. = FALSE
        )
      }
      status <- tryCatch(
        withr::with_dir(work_dir, {
          utils::zip(zipfile = output_path, files = files_rel, flags = "-r9Xq")
        }),
        error = function(e) 1L
      )
      if (!identical(as.integer(status), 0L) || !file.exists(output_path)) {
        stop(
          "Failed to create zip archive at '", output_path,
          "'. Install the 'zip' R package (install.packages(\"zip\")) or use `zip = FALSE`.",
          call. = FALSE
        )
      }
    }
  }

  invisible(output_path)
}
