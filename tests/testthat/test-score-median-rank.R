test_that("median skoru bilinen oyuncak matriste dogru", {
  expr <- matrix(
    c(1, 3, 5,
      2, 4, 100),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("G1", "G2"), c("S1", "S2", "S3"))
  )
  # Ornek medyanlari: S1=1.5, S2=3.5, S3=52.5
  res <- score_signature(expr, genes = c("G1", "G2"), method = "median", min_genes = 2)
  expect_false(res$skipped)
  expect_equal(unname(res$scores), c(1.5, 3.5, 52.5))
})

test_that("rank skoru ties ve sabit geni guvenli yonetir", {
  # G1: 1,2,3 -> ranks 1,2,3 -> (0, 0.5, 1)
  # G2: 5,5,5 -> ranks 2,2,2 -> (0.5, 0.5, 0.5)  [m=3, (2-1)/(3-1)=0.5]
  expr <- matrix(
    c(1, 2, 3,
      5, 5, 5),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("G1", "G2"), c("S1", "S2", "S3"))
  )
  res <- score_signature(expr, genes = c("G1", "G2"), method = "rank", min_genes = 2)
  expect_false(res$skipped)
  expect_equal(unname(res$scores), c(0.25, 0.5, 0.75), tolerance = 1e-12)
})

test_that("rank m<=1 kenar durumunda NA skor uretir/dislar", {
  expr <- matrix(
    c(1, NA_real_,
      NA_real_, NA_real_),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("G1", "G2"), c("S1", "S2"))
  )
  # Her gen icin m<=1 -> olcek NA; ornek ortalamalari NA -> dislanir
  res <- score_signature(expr, genes = c("G1", "G2"), method = "rank", min_genes = 2)
  expect_false(res$skipped)
  expect_true(res$n_all_na_dropped >= 1L)
})

test_that("median response etiketini kullanmaz (yalniz ifade)", {
  # Sutunlar: pos = (10,10), neg = (1,1)
  expr <- matrix(
    c(10, 10, 1, 1),
    nrow = 2,
    dimnames = list(c("A", "B"), c("pos", "neg"))
  )
  res <- score_signature(expr, genes = c("A", "B"), method = "median", min_genes = 2)
  expect_equal(unname(res$scores["pos"]), 10)
  expect_equal(unname(res$scores["neg"]), 1)
})
