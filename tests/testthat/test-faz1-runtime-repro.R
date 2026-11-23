# Faz 1 testleri: Publication Bundle manifest hash tutarliligi (1.1),
# Shiny snapshot degismezligi (1.2), merkezi min_genes (1.3) ve
# isaretli (genes_up/genes_down) ozel imza R kodu uretimi (1.4)

.faz1_load_app_env <- function() {
  app_path <- system.file("shiny", "app.R", package = "SigBenchR")
  if (!nzchar(app_path)) {
    app_path <- file.path("inst", "shiny", "app.R")
  }
  env <- new.env(parent = globalenv())
  source(app_path, local = env)
  env
}

test_that("Faz 1.1 & 1.2: download_manifest, download_pub_bundle, and download_report share identical reproducibility_hash under live input drift", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("pROC")
  skip_if_not_installed("rmarkdown")
  skip_if_not(rmarkdown::pandoc_available(), "pandoc yok")

  env <- .faz1_load_app_env()

  shiny::testServer(env$server, {
    session$setInputs(
      sig_choice = c("Cytokine-5", "IFNG6"),
      custom_genes = "",
      custom_genes_up = "",
      custom_genes_down = "",
      datasets = c("hugo_2016", "gide_2019"),
      methods = c("mean", "zscore"),
      direction = "<",
      sd_handling = "exclude",
      coverage_threshold = 0.70,
      sensitivity_threshold = 0.10,
      compute_pr = FALSE,
      compute_null = FALSE,
      n_boot = 0L,
      n_perm = 10L,
      match_expression = FALSE,
      seed = 42L,
      heatmap_method = "mean",
      forest_sig = "IFNG6"
    )
    session$setInputs(run = 1L)
    expect_false(stale_inputs())

    rm0 <- report_model()
    expect_identical(rm0$tier1$active_method, "mean")
    expect_identical(rm0$tier3$run_params$heatmap_method, "mean")
    expect_identical(rm0$tier3$params$heatmap_method, "mean")
    expect_identical(rm0$tier3$run_params$forest_signature, "IFNG6")

    man0 <- reproducibility_manifest()
    hash0 <- man0$reproducibility_hash
    expect_true(nzchar(hash0))

    # Changing view selectors or analysis parameters without re-running must flag stale_inputs()
    session$setInputs(heatmap_method = "zscore")
    expect_true(stale_inputs())
    session$setInputs(
      forest_sig = "Cytokine-5",
      sd_handling = "as_0",
      direction = ">",
      seed = 999L
    )
    expect_true(stale_inputs())

    # report_model() and reproducibility_manifest() must stay locked to the snapshot
    expect_identical(report_model()$tier1$active_method, "mean")
    expect_identical(reproducibility_manifest()$reproducibility_hash, hash0)

    # 1) Standalone manifest download
    man_file <- output$download_manifest
    expect_true(file.exists(man_file))
    man_standalone <- jsonlite::read_json(man_file, simplifyVector = TRUE)

    # 2) Publication bundle download
    zip_file <- output$download_pub_bundle
    expect_true(file.exists(zip_file))
    ex_dir <- tempfile("faz1_bundle_")
    dir.create(ex_dir)
    on.exit(unlink(ex_dir, recursive = TRUE), add = TRUE)
    utils::unzip(zip_file, exdir = ex_dir)
    man_bundle <- jsonlite::read_json(file.path(ex_dir, "manifest.json"), simplifyVector = TRUE)

    # 3) HTML report download
    rep_file <- output$download_report
    expect_true(file.exists(rep_file))
    rep_html <- paste(readLines(rep_file, warn = FALSE), collapse = "\n")

    # Acceptance criterion 1.1: identical reproducibility_hash and fingerprints across all exports
    expect_identical(man_standalone$reproducibility_hash, hash0)
    expect_identical(man_bundle$reproducibility_hash, hash0)
    expect_identical(man_bundle$cohorts_fingerprint, man_standalone$cohorts_fingerprint)
    expect_identical(man_bundle$signatures_fingerprint, man_standalone$signatures_fingerprint)
    expect_identical(man_bundle$run_params$heatmap_method, "mean")
    expect_identical(man_bundle$run_params$forest_signature, "IFNG6")
    expect_true(grepl(hash0, rep_html, fixed = TRUE))
  })
})

test_that("Faz 1.3: 2-gene Cytolytic_Activity (min_genes = 2L) is consistently scored across all package and Shiny entry points", {
  skip_if_not_installed("pROC")
  skip_if_not_installed("ggplot2")

  data(datasets_ici, package = "SigBenchR")
  b <- get_builtin_signatures(include_controls = FALSE)
  cyt <- b["Cytolytic_Activity"]
  expect_identical(cyt$Cytolytic_Activity$min_genes, 2L)
  expect_identical(SigBenchR:::.resolve_min_genes(cyt$Cytolytic_Activity, default = 3L), 2L)
  expect_identical(SigBenchR:::.resolve_min_genes(b$IFNG6, default = 3L), 3L)

  # Build a cohort containing GZMA and PRF1 (plus survival & reference axis genes)
  ds0 <- datasets_ici$van_allen_2015
  set.seed(101)
  extra_genes <- setdiff(c("GZMA", "PRF1", "GZMB"), rownames(ds0$expr))
  if (length(extra_genes) > 0L) {
    extra_mat <- matrix(
      stats::rnorm(length(extra_genes) * ncol(ds0$expr), mean = 4, sd = 1),
      nrow = length(extra_genes),
      ncol = ncol(ds0$expr),
      dimnames = list(extra_genes, colnames(ds0$expr))
    )
    ds0$expr <- rbind(ds0$expr, extra_mat)
    ds0$info$n_genes <- nrow(ds0$expr)
  }
  ds0$sample_meta$os_time <- stats::runif(nrow(ds0$sample_meta), min = 30, max = 900)
  ds0$sample_meta$os_event <- rep_len(c(1L, 0L, 1L), nrow(ds0$sample_meta))
  ds_pool <- list(van_allen_2015 = ds0)

  # 1) run_benchmark
  rb <- run_benchmark(cyt, ds_pool, methods = "mean", compute_pr = FALSE)
  expect_false(isTRUE(rb$skipped[1]))
  expect_identical(rb$n_found[1], 2L)

  # 2) plot_roc_curves
  p_roc <- plot_roc_curves(ds_pool, cyt, "van_allen_2015", method = "mean")
  expect_s3_class(p_roc, "ggplot")
  expect_true("Cytolytic_Activity" %in% p_roc$data$signature)

  # 3) score_for_dataset
  s_ds <- score_for_dataset(
    ds0,
    signature = cyt$Cytolytic_Activity,
    method = "mean"
  )
  expect_identical(s_ds$coverage$n_found, 2L)
  expect_true(length(s_ds$scores) > 0L)

  # 4) correlate_reference_axes
  ref_cor <- correlate_reference_axes(ds_pool, cyt, method = "mean")
  expect_true(nrow(ref_cor) > 0L)
  expect_true(all(ref_cor$signature == "Cytolytic_Activity"))

  # 5) run_survival_benchmark
  if (requireNamespace("survival", quietly = TRUE)) {
    sb <- run_survival_benchmark(ds_pool, cyt, method = "zscore")
    expect_identical(nrow(sb), 1L)
    expect_identical(sb$signature[1], "Cytolytic_Activity")
    expect_true(is.finite(sb$hr[1]))
  }

  # 6) Shiny Paired DeLong comparison with 2-gene Cytolytic_Activity
  if (requireNamespace("shiny", quietly = TRUE)) {
    env <- .faz1_load_app_env()
    env$datasets_ici$van_allen_2015 <- ds0
    shiny::testServer(env$server, {
      session$setInputs(
        sig_choice = c("Cytolytic_Activity", "IFNG6"),
        custom_genes = "",
        custom_genes_up = "",
        custom_genes_down = "",
        datasets = "van_allen_2015",
        methods = "mean",
        direction = "<",
        sd_handling = "exclude",
        coverage_threshold = 0.70,
        sensitivity_threshold = 0.10,
        compute_pr = FALSE,
        compute_null = FALSE,
        n_boot = 0L,
        n_perm = 10L,
        match_expression = FALSE,
        seed = 42L,
        heatmap_method = "mean",
        forest_sig = "Cytolytic_Activity",
        cmp_dataset = "van_allen_2015",
        cmp_sig1 = "Cytolytic_Activity",
        cmp_sig2 = "IFNG6"
      )
      session$setInputs(run = 1L)
      session$setInputs(run_cmp = 1L)

      cmp <- cmp_res()
      expect_s3_class(cmp, "data.frame")
      expect_true(is.finite(cmp$auc1[1]))
      expect_true(is.finite(cmp$auc2[1]))
      expect_true(is.finite(cmp$p_value[1]))
    })
  }
})

test_that("Faz 1.4: generate_r_code and Shiny r_code_text preserve signed custom signatures (genes_up and genes_down)", {
  skip_if_not_installed("pROC")

  data(datasets_ici, package = "SigBenchR")

  # 1) Signed custom signature only (no builtin signature)
  code_signed_only <- generate_r_code(
    signature_name = NULL,
    genes_up = c("IDO1", "CXCL9", "STAT1"),
    genes_down = c("ACTB", "GAPDH", "RPLP0"),
    methods = "mean",
    direction = "<",
    dataset_names = "hugo_2016",
    compute_pr = FALSE,
    compute_null = FALSE,
    seed = 42L
  )
  expect_silent(parse(text = code_signed_only))
  expect_true(grepl("genes_up = c(\"IDO1\", \"CXCL9\", \"STAT1\")", code_signed_only, fixed = TRUE))
  expect_true(grepl("genes_down = c(\"ACTB\", \"GAPDH\", \"RPLP0\")", code_signed_only, fixed = TRUE))

  # Evaluate generated R code and compare with direct run_benchmark()
  eval_lines <- strsplit(code_signed_only, "\n", fixed = TRUE)[[1]]
  eval_lines <- eval_lines[!grepl("^library\\(SigBenchR\\)$", eval_lines)]
  eval_lines <- eval_lines[!grepl("^(print|plot_auc_heatmap)\\(", eval_lines)]
  eval_env <- new.env(parent = asNamespace("SigBenchR"))
  utils::data("datasets_ici", package = "SigBenchR", envir = eval_env)
  suppressWarnings(eval(parse(text = paste(eval_lines, collapse = "\n")), envir = eval_env))

  direct_res <- suppressWarnings(run_benchmark(
    signatures = list(Custom = list(
      genes_up = c("IDO1", "CXCL9", "STAT1"),
      genes_down = c("ACTB", "GAPDH", "RPLP0"),
      direction = "<"
    )),
    datasets = datasets_ici["hugo_2016"],
    methods = "mean",
    direction = "<",
    compute_pr = FALSE,
    compute_null = FALSE,
    seed = 42L
  ))
  expect_equal(eval_env$res$auc, direct_res$auc, tolerance = 1e-12)
  expect_true(isTRUE(eval_env$res$signed[1]))

  # 2) Mutually exclusive flat genes and signed genes_up/genes_down
  expect_error(
    generate_r_code(
      signature_name = NULL,
      genes = c("IFNG", "STAT1"),
      genes_up = c("CXCL9"),
      genes_down = c("ACTB")
    ),
    "not both"
  )

  # 3) Shiny server passes snap$params$custom_genes_up and custom_genes_down to generate_r_code
  if (requireNamespace("shiny", quietly = TRUE)) {
    env <- .faz1_load_app_env()
    shiny::testServer(env$server, {
      session$setInputs(
        sig_choice = "IFNG6",
        custom_genes = "",
        custom_genes_up = "IDO1, CXCL9, STAT1",
        custom_genes_down = "ACTB, GAPDH, RPLP0",
        datasets = "hugo_2016",
        methods = "mean",
        direction = "<",
        sd_handling = "exclude",
        coverage_threshold = 0.70,
        sensitivity_threshold = 0.10,
        compute_pr = FALSE,
        compute_null = FALSE,
        n_boot = 0L,
        n_perm = 10L,
        match_expression = FALSE,
        seed = 42L,
        heatmap_method = "mean",
        forest_sig = "IFNG6"
      )
      session$setInputs(run = 1L)
      p <- results()$params
      code_from_snap <- generate_r_code(
        signature_name = p$sig_choice,
        genes = if (nzchar(trimws(p$custom_genes))) parse_genes(p$custom_genes) else NULL,
        genes_up = if (nzchar(trimws(p$custom_genes_up))) parse_genes(p$custom_genes_up) else NULL,
        genes_down = if (nzchar(trimws(p$custom_genes_down))) parse_genes(p$custom_genes_down) else NULL,
        methods = p$methods,
        direction = p$direction,
        coverage_threshold = p$coverage_threshold,
        n_boot = p$n_boot,
        n_perm = p$n_perm,
        seed = p$seed,
        dataset_names = p$datasets,
        compute_pr = p$compute_pr,
        compute_null = p$compute_null,
        sd_handling = p$sd_handling,
        match_expression = p$match_expression
      )
      expect_silent(parse(text = code_from_snap))
      expect_true(grepl("genes_up = c(\"IDO1\", \"CXCL9\", \"STAT1\")", code_from_snap, fixed = TRUE))
      expect_true(grepl("genes_down = c(\"ACTB\", \"GAPDH\", \"RPLP0\")", code_from_snap, fixed = TRUE))
    })
  }
})
