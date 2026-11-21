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

test_that("generate_r_code datasets_ici ile uyumlu ve eval edilebilir", {
  skip_if_not_installed("pROC")
  skip_if_not_installed("ggplot2")
  data(datasets_ici, package = "SigBenchR")
  data(signatures_demo, package = "SigBenchR")

  ici_names <- names(datasets_ici)
  expect_true(all(c("van_allen_2015", "hugo_2016", "riaz_2017", "gide_2019") %in% ici_names))

  code <- generate_r_code(
    signature_name = "IFNG6",
    methods = "mean",
    direction = "<",
    dataset_names = ici_names,
    compute_pr = FALSE,
    compute_null = FALSE,
    n_boot = 10L,
    n_perm = 10L,
    seed = 42L
  )
  expect_silent(parse(text = code))
  expect_true(grepl("data(datasets_ici)", code, fixed = TRUE))
  expect_false(grepl("data(datasets_demo)", code, fixed = TRUE))
  expect_true(grepl("ds <- datasets_ici[", code, fixed = TRUE))
  for (nm in ici_names) {
    expect_true(grepl(nm, code, fixed = TRUE))
  }

  # Shiny "R kodunu kopyala" ile ayni: secilen kohort adlari datasets_ici[...] uzerinden
  code_run <- sub("^library\\(SigBenchR\\)\\r?\\n?", "", code)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  eval_ok <- TRUE
  tryCatch(
    suppressWarnings(capture.output(eval(parse(text = code_run)))),
    error = function(e) {
      eval_ok <<- FALSE
      attr(eval_ok, "message") <<- conditionMessage(e)
    }
  )
  expect_true(eval_ok, info = paste0("eval hata: ", attr(eval_ok, "message")))
})

test_that("gide_2019 anti-PD1 monoterapi; van_allen klinik yarar n_pos=19", {
  data(datasets_ici, package = "SigBenchR")
  va <- datasets_ici$van_allen_2015
  expect_identical(unique(as.character(va$sample_meta$treatment)), "anti-CTLA4")
  expect_equal(as.integer(va$info$n_pos), 19L)
  expect_equal(as.integer(va$info$n_neg), 23L)
  expect_equal(as.integer(va$info$n), 42L)

  gide <- datasets_ici$gide_2019
  expect_identical(unique(as.character(gide$sample_meta$treatment)), "anti-PD1")
  expect_true(all(grepl("^PD1_", gide$sample_meta$sample_id)))
  expect_false(any(grepl("^ipiPD1_", gide$sample_meta$sample_id)))
  expect_equal(as.integer(gide$info$n), 41L)
  expect_equal(as.integer(gide$info$n_pos), 19L)
  expect_equal(as.integer(gide$info$n_neg), 16L)
  expect_equal(as.integer(gide$info$n_excluded_sd), 6L)
})
