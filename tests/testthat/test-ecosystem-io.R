test_that("read_gmt parses temporary GMT and drops empty tokens", {
  td <- tempfile("gmt_")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  gmt_path <- file.path(td, "toy.gmt")
  writeLines(
    c(
      "SIG_A\thttp://example\tIFNG\tCXCL9\t\tSTAT1",
      "",
      "SIG_B\tdesc\tCD274\tPDCD1"
    ),
    gmt_path
  )

  gs <- read_gmt(gmt_path)
  expect_equal(sort(names(gs)), c("SIG_A", "SIG_B"))
  expect_equal(gs$SIG_A, c("IFNG", "CXCL9", "STAT1"))
  expect_equal(gs$SIG_B, c("CD274", "PDCD1"))
})

test_that("as_sigbenchr_dataset validates lists and rejects bad objects", {
  ds <- list(
    expr = matrix(
      1:4,
      nrow = 2L,
      dimnames = list(c("A", "B"), c("S1", "S2"))
    ),
    sample_meta = data.frame(
      sample_id = c("S1", "S2"),
      response_label = c("CR", "PD"),
      stringsAsFactors = FALSE
    ),
    info = list(normalization = "TPM", scale = "log2")
  )
  expect_silent(as_sigbenchr_dataset(ds))
  expect_error(as_sigbenchr_dataset(1:5), "SummarizedExperiment|named dataset lists")
  expect_error(
    as_sigbenchr_dataset(list(expr = ds$expr)),
    "missing required fields|eksik alanlar"
  )
})
