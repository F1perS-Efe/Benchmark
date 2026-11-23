test_that("prepare_expression linear log2 uygular; negatifte hata", {
  m <- matrix(c(1, 3, 7, 15), nrow = 2, dimnames = list(c("A", "B"), c("S1", "S2")))
  out <- prepare_expression(m, scale = "linear", apply_log2 = TRUE)
  expect_equal(out, log2(m + 1))

  m_neg <- m
  m_neg[1, 1] <- -1
  expect_error(
    prepare_expression(m_neg, scale = "linear", apply_log2 = TRUE),
    "negative|Negatif"
  )
})

test_that("match_genes büyük harf tam eşleşme yapar; alias yok", {
  expr <- matrix(1:4, nrow = 2, dimnames = list(c("Ifng", "Cxcl9"), c("S1", "S2")))
  cov <- match_genes(expr, c("IFNG", "CXCL9", "IL21"))
  expect_equal(cov$n_found, 2)
  expect_equal(cov$missing, "IL21")
  expect_equal(cov$coverage_pct, 100 * 2 / 3)
})

test_that("encode_response SD'yi varsayılan dışlar", {
  labels <- c("CR", "PR", "PD", "SD", "sd")
  enc <- encode_response(labels, sd_handling = "exclude")
  expect_equal(enc$response, c(1L, 1L, 0L, NA, NA))
  expect_equal(enc$n_excluded_sd, 2L)
  expect_equal(enc$kept, c(TRUE, TRUE, TRUE, FALSE, FALSE))
})

test_that("validate_dataset zorunlu alanları kontrol eder", {
  ds <- list(
    expr = matrix(1:4, nrow = 2, dimnames = list(c("A", "B"), c("S1", "S2"))),
    sample_meta = data.frame(
      sample_id = c("S1", "S2"),
      response_label = c("CR", "PD"),
      stringsAsFactors = FALSE
    ),
    info = list(normalization = "TPM", scale = "log2")
  )
  expect_silent(validate_dataset(ds))
  bad <- ds
  bad$info$scale <- NULL
  expect_error(validate_dataset(bad), "scale")
})
