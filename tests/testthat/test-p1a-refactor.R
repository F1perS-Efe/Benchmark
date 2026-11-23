# P1-a refactor birim testleri: .prepare_cohort(), .as_signature(), tek dil (EN) mesajlar

test_that("P1-a: .prepare_cohort validates, deduplicates, scales, and optionally filters response", {
  prep_fn <- get(".prepare_cohort", envir = asNamespace("SigBenchR"))

  mat <- matrix(
    c(1, 3, 7, 15, 2, 4, 8, 16, 10, 20, 30, 40),
    nrow = 3L, ncol = 4L,
    dimnames = list(c("G1", "g1", "G2"), c("S1", "S2", "S3", "S4"))
  )
  ds <- list(
    expr = mat,
    sample_meta = data.frame(
      sample_id = colnames(mat),
      response_label = c("CR", "PD", "SD", "PR"),
      stringsAsFactors = FALSE
    ),
    info = list(normalization = "TPM", scale = "linear")
  )

  res_filt <- prep_fn(ds, name = "synth", apply_log2 = TRUE, sd_handling = "exclude", filter_response = TRUE)
  expect_identical(nrow(res_filt$expr_kept), 2L)
  expect_identical(colnames(res_filt$expr_kept), c("S1", "S2", "S4"))
  expect_identical(unname(res_filt$y_kept), c(1L, 0L, 1L))
  expect_equal(res_filt$expr_kept["g1", "S1"], log2(3 + 1), tolerance = 1e-12)

  res_unfilt <- prep_fn(ds, name = "synth", apply_log2 = TRUE, filter_response = FALSE)
  expect_identical(ncol(res_unfilt$expr_kept), 4L)
  expect_null(res_unfilt$y_kept)
})

test_that("P1-a: .as_signature normalizes character, numeric, and list inputs", {
  sig_fn <- get(".as_signature", envir = asNamespace("SigBenchR"))

  s_chr <- sig_fn(c("IFNG", "CXCL9", "STAT1"), name = "S1")
  expect_identical(s_chr, list(genes = c("IFNG", "CXCL9", "STAT1")))

  s_num <- sig_fn(c(IFNG = 1, VEGFA = -1), name = "S2")
  expect_identical(s_num, list(genes = c(IFNG = 1, VEGFA = -1)))

  s_lst <- list(genes_up = c("IFNG", "CXCL9"), genes_down = "VEGFA", direction = "signed")
  expect_identical(sig_fn(s_lst, name = "S3"), s_lst)

  expect_error(sig_fn(character(0), name = "Empty"), "empty")
  expect_error(sig_fn(list(direction = "positive"), name = "BadList"), "genes_up")
  expect_error(sig_fn(TRUE, name = "InvalidType"), "vector or a list")
})
