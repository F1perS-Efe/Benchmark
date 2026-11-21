test_that("datasets_ici yapisi, boyut ve on-isleme kurallari", {
  data(datasets_ici, package = "SigBenchR")
  expect_true(all(c("van_allen_2015", "hugo_2016", "riaz_2017", "gide_2019") %in% names(datasets_ici)))
  expect_false("liu_2019" %in% names(datasets_ici))

  for (nm in names(datasets_ici)) {
    ds <- datasets_ici[[nm]]
    expect_silent(validate_dataset(ds, name = nm))
    expect_true(nrow(ds$expr) >= 2500L)
    expect_true(nrow(ds$expr) <= 3000L)
    expect_equal(colnames(ds$expr), ds$sample_meta$sample_id)
    expect_true(all(ds$sample_meta$timing == "pre"))
    expect_false(anyDuplicated(ds$sample_meta$patient_id) > 0)
    expect_true(ds$info$n == nrow(ds$sample_meta))
    expect_true(ds$info$n_pos + ds$info$n_neg + ds$info$n_excluded_sd <= ds$info$n)
    # Imza genleri paketlenmis olmali (Random3 sentetik; bulunamayabilir)
    expect_true(all(c("IFNG", "CXCL9") %in% rownames(ds$expr)))
  }
})

test_that("Cytokine-5 derived_from ICI; IFNG6 bagimsiz", {
  data(signatures_demo, package = "SigBenchR")
  data(datasets_ici, package = "SigBenchR")
  expect_true(all(names(datasets_ici) %in% signatures_demo[["Cytokine-5"]]$derived_from))
  expect_identical(signatures_demo[["IFNG6"]]$derived_from, character(0))

  skip_if_not_installed("pROC")
  res <- run_benchmark(
    signatures = signatures_demo[c("Cytokine-5", "IFNG6")],
    datasets = datasets_ici["hugo_2016"],
    methods = "mean",
    compute_pr = FALSE,
    compute_null = FALSE
  )
  cyto <- res[res$signature == "Cytokine-5", ]
  ifng <- res[res$signature == "IFNG6", ]
  expect_true(all(cyto$independent == FALSE))
  expect_true(all(ifng$independent == TRUE))
})

test_that("datasets_ici uzerinde kucuk benchmark calisir", {
  skip_if_not_installed("pROC")
  data(datasets_ici, package = "SigBenchR")
  data(signatures_demo, package = "SigBenchR")
  res <- suppressWarnings(run_benchmark(
    signatures = signatures_demo["IFNG6"],
    datasets = datasets_ici[c("hugo_2016", "van_allen_2015")],
    methods = "mean",
    compute_pr = FALSE,
    compute_null = FALSE
  ))
  expect_true(nrow(res) >= 2)
  expect_true(any(is.finite(res$auc)))
})

test_that("data/ paketi 5 MB altinda", {
  data_dir <- system.file("data", package = "SigBenchR")
  if (!nzchar(data_dir) || !dir.exists(data_dir)) {
    # gelistirme: paket kokundeki data/
    data_dir <- file.path(testthat::test_path("..", ".."), "data")
  }
  skip_if_not(dir.exists(data_dir))
  bytes <- sum(file.info(list.files(data_dir, full.names = TRUE))$size, na.rm = TRUE)
  expect_lt(bytes, 5 * 1024^2)
})
