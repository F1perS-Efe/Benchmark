test_that("withr::with_seed ile skor üretimi tekrarlanabilir", {
  skip_if_not_installed("withr")

  make_mat <- function() {
    matrix(
      rnorm(20),
      nrow = 4,
      dimnames = list(paste0("G", 1:4), paste0("S", 1:5))
    )
  }

  m1 <- withr::with_seed(123, make_mat())
  m2 <- withr::with_seed(123, make_mat())
  expect_equal(m1, m2)

  s1 <- score_signature(m1, genes = paste0("G", 1:4), method = "mean", min_genes = 3)
  s2 <- score_signature(m2, genes = paste0("G", 1:4), method = "mean", min_genes = 3)
  expect_equal(s1$scores, s2$scores)
})
