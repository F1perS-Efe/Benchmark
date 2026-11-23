# FAZ 9 — R / Bioconductor Package Standards, Hygiene & Cross-Platform Portability Tests

.faz9_pkg_root <- function() {
  for (cand in c(".", "..", "../..")) {
    if (file.exists(file.path(cand, "DESCRIPTION")) && file.exists(file.path(cand, "NAMESPACE"))) {
      return(normalizePath(cand, mustWork = TRUE))
    }
  }
  "."
}

test_that("9.1: DESCRIPTION biocViews, LazyData, and roxygen2 version match package standards", {
  root <- .faz9_pkg_root()
  desc_path <- file.path(root, "DESCRIPTION")
  skip_if_not(file.exists(desc_path), "DESCRIPTION not found in package root")

  dcf <- as.list(read.dcf(desc_path)[1, ])
  expect_identical(dcf$LazyData, "false")
  expect_true("biocViews" %in% names(dcf))

  views <- trimws(strsplit(dcf$biocViews, ",", fixed = TRUE)[[1]])
  # Must not contain non-standard 'ImmunoOncology' term; must contain official Bioconductor views
  expect_false("ImmunoOncology" %in% views)
  expect_true(all(c("Software", "GeneExpression", "Transcriptomics", "Classification", "Survival") %in% views))

  roxy_desc <- system.file("DESCRIPTION", package = "roxygen2")
  if (nzchar(roxy_desc) && file.exists(roxy_desc)) {
    roxy_ver <- as.character(read.dcf(roxy_desc)[1, "Version"])
    expect_identical(dcf[["Config/roxygen2/version"]], roxy_ver)
  }
})

test_that("9.2: .Rbuildignore excludes temporary patches, scratch files, Rproj, and OS artifacts", {
  root <- .faz9_pkg_root()
  rbi_path <- file.path(root, ".Rbuildignore")
  skip_if_not(file.exists(rbi_path), ".Rbuildignore not found in package root")

  rbi_lines <- trimws(readLines(rbi_path, warn = FALSE))
  rbi_lines <- rbi_lines[nzchar(rbi_lines)]

  expect_true(any(grepl("Rproj", rbi_lines, fixed = TRUE)))
  expect_true(any(grepl("DS_Store", rbi_lines, fixed = TRUE)))
  expect_true(any(grepl("patch", rbi_lines, fixed = TRUE)))
  expect_true(any(grepl("scratch", rbi_lines, fixed = TRUE)))
  expect_true(any(grepl("renv", rbi_lines, fixed = TRUE)))
})

test_that("9.3: NAMESPACE registers S3 method as_sigbenchr_dataset.SummarizedExperiment and avoids internal helper leaks", {
  root <- .faz9_pkg_root()
  ns_path <- file.path(root, "NAMESPACE")
  skip_if_not(file.exists(ns_path), "NAMESPACE not found in package root")

  ns_lines <- trimws(readLines(ns_path, warn = FALSE))
  expect_true("S3method(as_sigbenchr_dataset,SummarizedExperiment)" %in% ns_lines)
  expect_true("S3method(as_sigbenchr_dataset,list)" %in% ns_lines)
  expect_true("S3method(as_sigbenchr_dataset,default)" %in% ns_lines)
  expect_true("importFrom(rlang,.data)" %in% ns_lines)
  expect_true("importFrom(tools,md5sum)" %in% ns_lines)

  # No dot-prefixed internal functions may be exported
  ns_exports <- getNamespaceExports("SigBenchR")
  expect_false(any(startsWith(ns_exports, ".")))

  # Test S3 dispatch on a plain list and clean error on invalid class
  data(datasets_demo, package = "SigBenchR")
  ds_ok <- as_sigbenchr_dataset(datasets_demo$demo_A)
  expect_true(is.matrix(ds_ok$expr))
  expect_error(as_sigbenchr_dataset(list(a = 1)), "missing required fields")
  expect_error(as_sigbenchr_dataset(42L), "SummarizedExperiment")
})

test_that("9.4: Cross-platform POSIX and Windows drive-letter path portability in .make_relative_path and system.file resolution", {
  # 1) POSIX relative path
  expect_identical(
    SigBenchR:::.make_relative_path("/Users/test/project/data/cohort.rds", "/Users/test/project"),
    "data/cohort.rds"
  )
  expect_identical(
    SigBenchR:::.make_relative_path("/Users/test/other/cohort.rds", "/Users/test/project/reg"),
    "../../other/cohort.rds"
  )

  # 2) Windows same-drive (case-insensitive drive letter + backslash normalization)
  expect_identical(
    SigBenchR:::.make_relative_path("C:\\Users\\test\\project\\data\\cohort.rds", "c:/Users/test/project"),
    "data/cohort.rds"
  )

  # 3) Windows different drives (C:/ vs D:/) must NOT produce invalid relative paths across drives
  expect_identical(
    SigBenchR:::.make_relative_path("D:\\cohorts\\data\\cohort.rds", "C:\\Users\\test\\project"),
    "D:/cohorts/data/cohort.rds"
  )

  # 4) system.file resolution works from any working directory without relying on local inst/extdata
  withr::with_dir(tempdir(), {
    rules <- load_interpretation_rules()
    expect_true(is.list(rules))
    expect_true(nzchar(rules$md5))

    reg_tbl <- list_cohort_registry()
    expect_true(is.data.frame(reg_tbl))
    expect_true("van_allen_2015" %in% reg_tbl$cohort_id)
  })
})
