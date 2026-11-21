test_that("mean skoru bilinen oyuncak matriste doğru", {
  expr <- matrix(
    c(1, 3, 5,
      2, 4, 6),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("G1", "G2"), c("S1", "S2", "S3"))
  )
  res <- score_signature(expr, genes = c("G1", "G2"), method = "mean", min_genes = 2)
  expect_false(res$skipped)
  expect_equal(unname(res$scores), c(1.5, 3.5, 5.5))
})

test_that("zscore etiketten bağımsız ve örnek ortalaması üretir", {
  set.seed(1)
  expr <- matrix(
    rnorm(30),
    nrow = 3,
    dimnames = list(c("A", "B", "C"), paste0("S", 1:10))
  )
  res <- score_signature(expr, genes = c("A", "B", "C"), method = "zscore", min_genes = 3)
  expect_false(res$skipped)
  expect_length(res$scores, 10)
  # Gen bazlı z-score sonrası satır ortalamaları ~0 civarı olmalı (tüm örnekler)
  expect_true(abs(mean(res$scores)) < 0.5)
})

test_that("min_genes altında atlanır", {
  expr <- matrix(1:6, nrow = 2, dimnames = list(c("A", "B"), paste0("S", 1:3)))
  res <- score_signature(expr, genes = c("A", "X", "Y", "Z"), method = "mean", min_genes = 3)
  expect_true(res$skipped)
  expect_match(res$skip_reason, "min_genes")
})
