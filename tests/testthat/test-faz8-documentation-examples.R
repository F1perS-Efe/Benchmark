# FAZ 8 — Documentation, Vignettes, Scientific Transparency & Rd Integrity Tests

.faz8_pkg_root <- function() {
  for (cand in c(".", "..", "../..")) {
    if (file.exists(file.path(cand, "DESCRIPTION")) && file.exists(file.path(cand, "README.md"))) {
      return(normalizePath(cand, mustWork = TRUE))
    }
  }
  "."
}

test_that("8.1: All exported functions have synchronized man/*.Rd files with no missing argument documentation", {
  ns_exports <- getNamespaceExports("SigBenchR")
  expect_gt(length(ns_exports), 40L)

  root <- .faz8_pkg_root()
  man_dir <- file.path(root, "man")
  if (!dir.exists(man_dir)) {
    man_dir <- system.file("man", package = "SigBenchR")
  }
  skip_if_not(dir.exists(man_dir), "man/ directory not located in source tree")

  rd_files <- list.files(man_dir, pattern = "\\.Rd$", full.names = TRUE)
  expect_gt(length(rd_files), 40L)

  # Key functions updated in Phases 1-7 must document their new parameters and return fields
  key_checks <- list(
    "score_signature.Rd" = c("weights", "genes_up", "genes_down"),
    "null_permutation_test.Rd" = c("weights", "background_type", "universe_size"),
    "adjust_null_p.Rd" = c("scope", "null_p_adj_bh_global"),
    "test_beyond_infiltration.Rd" = c("sd_delta_auc", "delta_auc_ci_low", "stouffer_z", "direction_conflict", "vif", "epv"),
    "pr_auc.Rd" = c("tie_fraction", "small_positive_class", "n_boot_valid", "ci_unreliable"),
    "run_survival_benchmark.Rd" = c("logrank_p_adj_bh", "p_adj_bh", "ph_p"),
    "map_mouse_to_human.Rd" = c("ortholog_table", "ortholog_policy", "keep_unmapped_context"),
    "build_report_model.Rd" = c("beyond_inf", "surv_res"),
    "export_report.Rd" = c("beyond_inf", "surv_res"),
    "export_publication_bundle.Rd" = c("beyond_inf", "surv_res", "analysis_objects.rds"),
    "build_reproducibility_manifest.Rd" = c("diagnostics_summary", "reproducibility_hash")
  )

  for (rd_name in names(key_checks)) {
    rd_path <- file.path(man_dir, rd_name)
    expect_true(file.exists(rd_path), info = paste("Missing Rd:", rd_name))
    txt <- paste(readLines(rd_path, warn = FALSE), collapse = "\n")
    for (kw in key_checks[[rd_name]]) {
      expect_true(grepl(kw, txt, fixed = TRUE), info = paste(rd_name, "missing keyword:", kw))
    }
  }
})

test_that("8.2: README.md and NEWS.md accurately reflect Phase 1-8 architecture and scientific disclaimers", {
  root <- .faz8_pkg_root()
  readme_path <- file.path(root, "README.md")
  news_path <- file.path(root, "NEWS.md")
  skip_if_not(file.exists(readme_path), "README.md not in package root")
  skip_if_not(file.exists(news_path), "NEWS.md not in package root")

  readme_txt <- paste(readLines(readme_path, warn = FALSE), collapse = "\n")
  news_txt <- paste(readLines(news_path, warn = FALSE), collapse = "\n")

  # Stale false claims must be absent from README.md
  expect_false(grepl("İmza seçimi varsayılan olarak `signatures_demo` (Cytokine-5, IFNG6, Random3) ile sınırlıdır", readme_txt, fixed = TRUE))

  # New Phase 1-8 scientific and statistical features must be documented in README.md and NEWS.md
  for (term in c(
    "T_Cell_Inflamed_GEP", "CD276", "B7_Immune_Exclusion",
    "null_p_adj_bh_global", "background_type",
    "stouffer_z", "direction_conflict", "lrt_p_adj_bh",
    "logrank_p_adj_bh", "ci_unreliable", "small_positive_class",
    "strict_1to1", "transcriptome_bg_expr", "NanoString nCounter"
  )) {
    expect_true(grepl(term, readme_txt, fixed = TRUE), info = paste("README.md missing:", term))
    expect_true(grepl(term, news_txt, fixed = TRUE), info = paste("NEWS.md missing:", term))
  }
})

test_that("8.3: Methodology vignette (vignettes/SigBenchR-methodology.Rmd) parses and executes cleanly", {
  root <- .faz8_pkg_root()
  vig_path <- file.path(root, "vignettes", "SigBenchR-methodology.Rmd")
  skip_if_not(file.exists(vig_path), "Vignette file not in package root")
  skip_if_not_installed("knitr")

  # Extract R code chunks and execute them in an isolated environment
  tmp_r <- tempfile("faz8_vig_", fileext = ".R")
  on.exit(unlink(tmp_r), add = TRUE)
  knitr::purl(vig_path, output = tmp_r, quiet = TRUE)
  expect_true(file.exists(tmp_r))
  expect_error(parse(file = tmp_r), NA)

  vig_env <- new.env(parent = globalenv())
  expect_error(source(tmp_r, local = vig_env), NA)
  expect_true(is.data.frame(vig_env$res))
  expect_true(is.list(vig_env$bi))
  expect_true("null_p_adj_bh_global" %in% names(vig_env$res))
  expect_true("stouffer_z" %in% names(vig_env$bi$summary))
})
