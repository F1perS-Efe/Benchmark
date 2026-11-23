# R4 testleri: Bioconductor/R CMD check standartlari, izolasyon, kod uretimi ve tasinabilirlik

.r4_pkg_file <- function(...) {
  p <- system.file(..., package = "SigBenchR")
  if (!nzchar(p) || !file.exists(p)) {
    p <- file.path(...)
  }
  p
}

test_that("R4-1: inst/ files contain zero internal ::: calls", {
  inst_files <- c(
    .r4_pkg_file("shiny", "app.R"),
    .r4_pkg_file("rmd", "benchmark_report.Rmd")
  )
  for (f in inst_files) {
    expect_true(file.exists(f))
    txt <- paste(readLines(f, warn = FALSE), collapse = "\n")
    expect_false(grepl(":::", txt, fixed = TRUE), info = paste("Found ::: in", f))
  }
})

test_that("R4-1: DESCRIPTION, NAMESPACE, and ecosystem_io follow Bioconductor/CRAN standards", {
  desc_path <- .r4_pkg_file("DESCRIPTION")
  expect_true(file.exists(desc_path))
  dcf <- read.dcf(desc_path)
  expect_identical(unname(dcf[1, "LazyData"]), "false")
  expect_true("biocViews" %in% colnames(dcf))
  expect_true(nzchar(dcf[1, "biocViews"]))
  expect_true(grepl("\\brlang\\b", dcf[1, "Imports"]))
  expect_true(grepl("\\bSummarizedExperiment\\b", dcf[1, "Suggests"]))

  ns_path <- .r4_pkg_file("NAMESPACE")
  expect_true(file.exists(ns_path))
  ns_lines <- readLines(ns_path, warn = FALSE)
  expect_true(any(grepl("^importFrom\\(rlang,\\.data\\)$", ns_lines)))
  expect_false(any(grepl("globalVariables", ns_lines, fixed = TRUE)))
  expect_true("export(format_meta_summary)" %in% ns_lines)
  expect_true("export(normalize_roc_direction)" %in% ns_lines)

  # No paste0("Summarized", "Experiment") obfuscation in exported S3 method
  fn_txt <- paste(deparse(SigBenchR:::as_sigbenchr_dataset.default), collapse = "\n")
  expect_false(grepl('paste0\\("Summarized"', fn_txt))
  expect_true(grepl("SummarizedExperiment", fn_txt, fixed = TRUE))
})

test_that("R4-1: format_meta_summary and normalize_roc_direction are exported and match internal aliases", {
  expect_identical(normalize_roc_direction("negative"), ">")
  expect_identical(normalize_roc_direction("positive"), "<")
  expect_identical(normalize_roc_direction("signed"), "<")
  expect_identical(normalize_roc_direction(NULL, default = ">"), ">")

  sm <- data.frame(
    signature = "IFNG6",
    method = "mean",
    n_cohorts = 2L,
    auc_n_weighted = 0.72,
    auc_ivw = 0.73,
    auc_re = 0.74,
    auc_re_ci_low = 0.60,
    auc_re_ci_high = 0.85,
    i2_logit = 12.5,
    tau2_reml = 0.015,
    note = NA_character_,
    stringsAsFactors = FALSE
  )
  out <- format_meta_summary(sm)
  expect_s3_class(out, "data.frame")
  expect_true("REML AUC (95% CI)" %in% names(out))
})

test_that("R4-2: generate_r_code omits liu_2019 and produces valid runnable code for get_builtin_signatures()", {
  skip_if_not_installed("pROC")
  skip_if_not_installed("ggplot2")

  # liu_2019 is not in datasets_ici and must not emit datasets_ici["liu_2019"]
  code_liu <- generate_r_code(
    signature_name = "IFNG6",
    methods = "mean",
    dataset_names = c("hugo_2016", "liu_2019")
  )
  expect_silent(parse(text = code_liu))
  expect_false(grepl("datasets_ici", code_liu, fixed = TRUE))

  # get_builtin_signatures() panel support (including negative-direction signature)
  code_bi <- generate_r_code(
    signature_name = c("T_Cell_Inflamed_GEP", "Myeloid_Suppression"),
    methods = "mean",
    direction = "<",
    dataset_names = c("hugo_2016", "gide_2019"),
    compute_pr = FALSE,
    compute_null = FALSE
  )
  expect_silent(parse(text = code_bi))
  expect_true(grepl("get_builtin_signatures", code_bi, fixed = TRUE))
  expect_true(grepl("T_Cell_Inflamed_GEP", code_bi, fixed = TRUE))

  # Evaluate generated script in isolated environment (replacing library call cleanly for check/load_all)
  eval_lines <- strsplit(code_bi, "\n", fixed = TRUE)[[1]]
  eval_lines <- eval_lines[!grepl("^library\\(SigBenchR\\)$", eval_lines)]
  # Replace final plot/print calls so output captures `res` cleanly
  eval_lines <- eval_lines[!grepl("^(print|plot_auc_heatmap)\\(", eval_lines)]
  env <- new.env(parent = asNamespace("SigBenchR"))
  utils::data("datasets_ici", package = "SigBenchR", envir = env)
  suppressWarnings(eval(parse(text = paste(eval_lines, collapse = "\n")), envir = env))
  expect_true(is.data.frame(env$res))
  expect_true(all(c("T_Cell_Inflamed_GEP", "Myeloid_Suppression") %in% env$res$signature))
})

test_that("R4-2: export_report renders in isolated env decoupled from .GlobalEnv", {
  skip_if_not_installed("rmarkdown")
  skip_if_not_installed("knitr")
  skip_if_not_installed("pROC")
  skip_if_not(rmarkdown::pandoc_available(), "pandoc yok")

  data(datasets_ici, package = "SigBenchR")
  data(signatures_demo, package = "SigBenchR")
  res <- suppressWarnings(run_benchmark(
    signatures = signatures_demo[c("Cytokine-5", "IFNG6")],
    datasets = datasets_ici[c("hugo_2016", "gide_2019")],
    methods = "mean",
    compute_pr = FALSE,
    compute_null = FALSE,
    seed = 42L
  ))
  sm <- suppressWarnings(summarize_auc(res, method = "mean"))
  sens <- check_sensitivity(res, threshold = 0.10)

  # Confirm export_report() does not parent render env from globalenv()
  fn_txt <- paste(deparse(export_report), collapse = "\n")
  expect_false(grepl("globalenv\\(", fn_txt))

  # Poison .GlobalEnv with identically named variables and verify clean render
  had_results <- exists("results", envir = .GlobalEnv, inherits = FALSE)
  old_results <- if (had_results) get("results", envir = .GlobalEnv, inherits = FALSE) else NULL
  had_hm <- exists("hm", envir = .GlobalEnv, inherits = FALSE)
  old_hm <- if (had_hm) get("hm", envir = .GlobalEnv, inherits = FALSE) else NULL
  on.exit({
    if (had_results) {
      assign("results", old_results, envir = .GlobalEnv)
    } else if (exists("results", envir = .GlobalEnv, inherits = FALSE)) {
      rm("results", envir = .GlobalEnv)
    }
    if (had_hm) {
      assign("hm", old_hm, envir = .GlobalEnv)
    } else if (exists("hm", envir = .GlobalEnv, inherits = FALSE)) {
      rm("hm", envir = .GlobalEnv)
    }
  }, add = TRUE)

  assign("results", "POISON_GLOBAL_RESULTS", envir = .GlobalEnv)
  assign("hm", "POISON_GLOBAL_HM", envir = .GlobalEnv)

  out_html <- tempfile(fileext = ".html")
  on.exit(unlink(out_html), add = TRUE)
  p <- export_report(
    results = res,
    summary_tbl = sm,
    sensitivity_tbl = sens,
    params_list = list(
      methods = "mean",
      direction = "<",
      coverage_threshold = 0.70,
      sensitivity_threshold = 0.10,
      n_boot = 0L,
      n_perm = 0L,
      seed = 42L,
      heatmap_method = "mean",
      forest_signature = "IFNG6"
    ),
    output_file = out_html
  )
  expect_true(file.exists(p))
  html <- paste(readLines(p, warn = FALSE), collapse = "\n")
  expect_true(grepl("IFNG6", html, fixed = TRUE))
})

test_that("R4-2: manifest omits unused metafor and uses in-memory md5; cohort registry stores portable relative paths", {
  data(datasets_demo, package = "SigBenchR")
  data(signatures_demo, package = "SigBenchR")
  res <- suppressWarnings(run_benchmark(
    signatures_demo["IFNG6"],
    datasets_demo["demo_A"],
    methods = "mean",
    compute_pr = FALSE,
    compute_null = FALSE
  ))
  man <- build_reproducibility_manifest(res)
  expect_false("metafor" %in% names(man$package_versions))

  tf_md5 <- tempfile(fileext = ".txt")
  on.exit(unlink(tf_md5), add = TRUE)
  writeBin(charToRaw(enc2utf8("sigbenchr-test")), tf_md5)
  expect_identical(
    SigBenchR:::.md5_string("sigbenchr-test"),
    unname(tools::md5sum(tf_md5))
  )

  # Relative path portability test in cohort_registry
  root_a <- tempfile("reg_port_a_")
  dir.create(file.path(root_a, "data"), recursive = TRUE)
  on.exit(unlink(root_a, recursive = TRUE), add = TRUE)

  rds_a <- file.path(root_a, "data", "toy.rds")
  saveRDS(datasets_demo$demo_A, rds_a)
  reg_a <- file.path(root_a, "registry.yml")

  row <- register_cohort(
    cohort_id = "toy_rel",
    file_path = rds_a,
    registry_path = reg_a,
    include_builtin = FALSE
  )
  expect_identical(row$file_path, "data/toy.rds")

  # Move/copy directory to root_b and load from a different working directory
  root_b <- tempfile("reg_port_b_")
  dir.create(file.path(root_b, "data"), recursive = TRUE)
  on.exit(unlink(root_b, recursive = TRUE), add = TRUE)
  file.copy(rds_a, file.path(root_b, "data", "toy.rds"))
  file.copy(reg_a, file.path(root_b, "registry.yml"))

  ds_b <- withr::with_dir(tempdir(), {
    load_registered_cohort("toy_rel", registry_path = file.path(root_b, "registry.yml"))
  })
  expect_identical(nrow(ds_b$expr), nrow(datasets_demo$demo_A$expr))
})

test_that("R4-3: read_user_dataset and load_cohort_bundle messages are English without v0.2 reference", {
  expect_error(
    load_cohort_bundle({
      tf <- tempfile(fileext = ".rds")
      on.exit(unlink(tf), add = TRUE)
      saveRDS(list(1, 2), tf)
      tf
    }),
    "named list"
  )
  td <- tempfile("r4_read_")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  ep <- file.path(td, "e.csv")
  mp <- file.path(td, "m.csv")
  writeLines(c("gene,S1,S2,S3,S4", "IFNG,1,2,3,4", "CXCL9,2,3,4,5"), ep)
  writeLines(c("sample_id,response_label", "S1,CR", "S2,PR", "S3,PD", "S4,PD"), mp)
  w_msg <- NULL
  withCallingHandlers(
    read_user_dataset(ep, mp, gene_id_type = "Ensembl", scale = "log2"),
    warning = function(w) {
      w_msg <<- paste(w_msg, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_false(grepl("v0\\.2", w_msg))
  expect_true(grepl("signature matching expects HGNC", w_msg, fixed = TRUE))
})
