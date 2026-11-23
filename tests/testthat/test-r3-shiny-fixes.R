# R3 testleri: Shiny reaktif guvenligi, yukleme siniri, hata korumalari ve snapshot tutarliligi

.r3_load_app_env <- function() {
  app_path <- system.file("shiny", "app.R", package = "SigBenchR")
  if (!nzchar(app_path)) {
    app_path <- file.path("inst", "shiny", "app.R")
  }
  env <- new.env(parent = globalenv())
  source(app_path, local = env)
  env
}

test_that("R3-1: run_app sets 100 MB shiny.maxRequestSize and restores options on exit", {
  skip_if_not_installed("shiny")

  old_size <- getOption("shiny.maxRequestSize")
  old_allow <- shiny::getShinyOption("sigbenchr_allow_upload")
  on.exit(
    {
      options(shiny.maxRequestSize = old_size)
      shiny::shinyOptions(sigbenchr_allow_upload = old_allow)
    },
    add = TRUE
  )

  options(shiny.maxRequestSize = 12345)
  shiny::shinyOptions(sigbenchr_allow_upload = FALSE)

  seen_inside <- list()
  testthat::local_mocked_bindings(
    runApp = function(appDir, ...) {
      seen_inside$max_req <<- getOption("shiny.maxRequestSize")
      seen_inside$allow <<- shiny::getShinyOption("sigbenchr_allow_upload")
      invisible(TRUE)
    },
    .package = "shiny"
  )

  run_app(allow_upload = TRUE)

  expect_identical(seen_inside$max_req, 100 * 1024^2)
  expect_true(isTRUE(seen_inside$allow))
  expect_identical(getOption("shiny.maxRequestSize"), 12345)
  expect_identical(shiny::getShinyOption("sigbenchr_allow_upload"), FALSE)
})

test_that("R3-1: UI exposes get_builtin_signatures, genes_up/down, package version, and removes dead helpers", {
  skip_if_not_installed("shiny")

  env <- .r3_load_app_env()
  ui_html <- as.character(env$ui)

  builtin_names <- names(get_builtin_signatures(include_controls = TRUE))
  expect_true(all(builtin_names %in% names(env$.builtin_sigs_ui)))
  expect_true(grepl("T_Cell_Inflamed_GEP", ui_html, fixed = TRUE))
  expect_true(grepl("Myeloid_Suppression", ui_html, fixed = TRUE))
  expect_true(grepl("custom_genes_up", ui_html, fixed = TRUE))
  expect_true(grepl("custom_genes_down", ui_html, fixed = TRUE))
  expect_true(grepl(
    paste0("v", as.character(utils::packageVersion("SigBenchR"))),
    ui_html,
    fixed = TRUE
  ))

  # Dead formatters must be removed from server environment
  shiny::testServer(env$server, {
    srv_env <- environment()
    for (dead in c(".em", ".en", ".fmt3", ".fmt_i2", ".fmt_auc_ci", ".translate_note")) {
      expect_false(exists(dead, envir = srv_env, inherits = FALSE))
    }
  })
})

test_that("R3-1 & R3-2: single-signature guards, numeric validate(), sensitivity stale_inputs, and live snapshot guards work", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("pROC")

  env <- .r3_load_app_env()

  shiny::testServer(env$server, {
    # Numeric guards before run_benchmark()
    session$setInputs(
      sig_choice = c("IFNG6", "Cytokine-5"),
      custom_genes = "",
      datasets = "hugo_2016",
      methods = "mean",
      direction = "<",
      sd_handling = "exclude",
      coverage_threshold = 0.70,
      sensitivity_threshold = 0.10,
      compute_pr = TRUE,
      compute_null = FALSE,
      n_boot = 0L,
      n_perm = 10L,
      match_expression = FALSE,
      seed = 42L,
      heatmap_method = "mean",
      forest_sig = "IFNG6"
    )
    session$setInputs(run = 1L)
    expect_error(results(), class = "shiny.silent.error")

    session$setInputs(compute_pr = FALSE, n_boot = NA_integer_)
    session$setInputs(run = 2L)
    expect_error(results(), class = "shiny.silent.error")

    session$setInputs(n_boot = 0L, n_perm = 0L)
    session$setInputs(run = 3L)
    expect_error(results(), class = "shiny.silent.error")

    session$setInputs(n_perm = 10L, seed = 0L)
    session$setInputs(run = 4L)
    expect_error(results(), class = "shiny.silent.error")

    # Valid run with single signature (plus sensitivity_threshold stale check)
    session$setInputs(
      sig_choice = "IFNG6",
      datasets = c("hugo_2016", "gide_2019"),
      methods = "mean",
      seed = 42L,
      roc_dataset = "hugo_2016",
      forest_sig = "IFNG6"
    )
    session$setInputs(run = 5L)
    expect_false(stale_inputs())

    # Single signature must cleanly trigger validate(need()) instead of crashing
    expect_error(output$jaccard_fallback, class = "shiny.silent.error")
    expect_error(output$corr_fallback, class = "shiny.silent.error")
    if (requireNamespace("DT", quietly = TRUE)) {
      expect_error(output$jaccard_table, class = "shiny.silent.error")
      expect_error(output$corr_table, class = "shiny.silent.error")
    }

    # Changing sensitivity_threshold must flag stale_inputs()
    session$setInputs(sensitivity_threshold = 0.25)
    expect_true(stale_inputs())
    session$setInputs(sensitivity_threshold = 0.10)
    expect_false(stale_inputs())

    # Live input mismatches against snapshot must trigger validate(need())
    session$setInputs(heatmap_method = "rank")
    expect_error(summary_res(), class = "shiny.silent.error")
    expect_error(summary_res_forest(), class = "shiny.silent.error")
    expect_error(output$heatmap, class = "shiny.silent.error")

    session$setInputs(heatmap_method = "mean", roc_dataset = "van_allen_2015")
    expect_error(output$roc_plot, class = "shiny.silent.error")

    session$setInputs(roc_dataset = "hugo_2016", forest_sig = "Random3")
    expect_error(output$forest_plot, class = "shiny.silent.error")
  })
})

test_that("R3-2: build_sigs preserves builtin signature directions and supports signed custom genes_up/genes_down", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("pROC")

  env <- .r3_load_app_env()

  shiny::testServer(env$server, {
    # Built-in Myeloid_Suppression ("negative") and IFNG6 ("positive") must NOT be
    # overwritten by global direction = ">"
    sigs_preserved <- build_sigs(
      direction = ">",
      sig_choice = c("IFNG6", "Myeloid_Suppression"),
      custom_genes = "",
      custom_genes_up = "CXCL9, CXCL10, STAT1",
      custom_genes_down = "VEGFA, TGFB1"
    )
    expect_identical(sigs_preserved[["IFNG6"]]$direction, "positive")
    expect_identical(sigs_preserved[["Myeloid_Suppression"]]$direction, "negative")
    expect_identical(sigs_preserved[["Custom"]]$genes_up, c("CXCL9", "CXCL10", "STAT1"))
    expect_identical(sigs_preserved[["Custom"]]$genes_down, c("VEGFA", "TGFB1"))
    expect_identical(sigs_preserved[["Custom"]]$direction, ">")

    session$setInputs(
      sig_choice = c("IFNG6", "Myeloid_Suppression"),
      custom_genes = "",
      custom_genes_up = "CXCL9, CXCL10, STAT1",
      custom_genes_down = "VEGFA, TGFB1",
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

    snap <- results()
    expect_identical(snap$sigs[["IFNG6"]]$direction, "positive")
    expect_identical(snap$sigs[["Myeloid_Suppression"]]$direction, "negative")
    expect_identical(snap$sigs[["Custom"]]$genes_up, c("CXCL9", "CXCL10", "STAT1"))
    expect_identical(snap$sigs[["Custom"]]$genes_down, c("VEGFA", "TGFB1"))
    custom_row <- snap$res[snap$res$signature == "Custom", , drop = FALSE]
    expect_identical(nrow(custom_row), 1L)
    expect_false(isTRUE(custom_row$skipped[1]))
  })
})

test_that("R3-1: download_report and download_pub_bundle both use snapshot forest_signature and heatmap_method", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("pROC")
  skip_if_not_installed("rmarkdown")
  skip_if_not(rmarkdown::pandoc_available(), "pandoc yok")

  env <- .r3_load_app_env()

  shiny::testServer(env$server, {
    session$setInputs(
      sig_choice = c("Cytokine-5", "IFNG6"),
      custom_genes = "",
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

    # Drift live view selectors WITHOUT pressing Run analysis
    session$setInputs(
      heatmap_method = "zscore",
      forest_sig = "Cytokine-5"
    )

    rep_path <- output$download_report
    expect_true(file.exists(rep_path))
    html <- paste(readLines(rep_path, warn = FALSE), collapse = "\n")
    expect_true(grepl("IFNG6", html, fixed = TRUE))

    zip_path <- output$download_pub_bundle
    expect_true(file.exists(zip_path))
    expect_gt(file.info(zip_path)$size, 0)

    ex_dir <- tempfile("r3_bundle_")
    dir.create(ex_dir)
    on.exit(unlink(ex_dir, recursive = TRUE), add = TRUE)
    utils::unzip(zip_path, exdir = ex_dir)
    legends_file <- file.path(ex_dir, "figure_legends.txt")
    expect_true(file.exists(legends_file))
    legends <- paste(readLines(legends_file, warn = FALSE), collapse = "\n")
    expect_true(grepl("forest plot for IFNG6 (mean)", legends, fixed = TRUE))
    expect_true(grepl("scoring method: mean;", legends, fixed = TRUE))
    expect_false(grepl("Cytokine-5 (zscore)", legends, fixed = TRUE))

    man <- jsonlite::read_json(file.path(ex_dir, "manifest.json"), simplifyVector = TRUE)
    expect_identical(man$run_params$forest_signature, "IFNG6")
    expect_identical(man$run_params$heatmap_method, "mean")
  })
})
