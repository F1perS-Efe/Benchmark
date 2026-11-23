test_that("skipped=TRUE satirlar vektor-guvenli dislanir", {
  skip_if_not_installed("ggplot2")

  res <- data.frame(
    dataset = c("A", "B", "C", "D"),
    signature = c("S1", "S1", "S1", "S2"),
    method = "mean",
    n = c(20L, 30L, 25L, 40L),
    auc = c(0.70, 0.80, NA_real_, 0.65),
    auc_ci_low = c(0.5, 0.6, NA_real_, 0.5),
    auc_ci_high = c(0.9, 1.0, NA_real_, 0.8),
    auc_var = c(0.01, 0.02, NA_real_, 0.01),
    coverage_pct = c(100, 100, 40, 100),
    low_coverage = c(FALSE, FALSE, TRUE, FALSE),
    independent = TRUE,
    var_zero = FALSE,
    small_cohort = FALSE,
    skipped = c(FALSE, FALSE, TRUE, FALSE),
    stringsAsFactors = FALSE
  )

  # Onceki !isTRUE(skipped) vektor bug'i burada patlar / atlananlari dislemezdi
  expect_silent(hm <- plot_auc_heatmap(res, method = "mean"))
  expect_s3_class(hm, "ggplot")

  sm <- suppressWarnings(summarize_auc(res, method = "mean"))
  expect_equal(sm$n_datasets[sm$signature == "S1"], 2L)
  expect_false(any(is.na(sm$auc_n_weighted)))

  sens <- check_sensitivity(res, threshold = 0.10)
  expect_false(any(sens$dataset == "C" & sens$signature == "S1"))
  expect_equal(nrow(sens), 3L)
  expect_setequal(paste(sens$dataset, sens$signature), c("A S1", "B S1", "D S2"))

  p <- plot_auc_forest(res, signature = "S1", method = "mean", summary = sm)
  expect_s3_class(p, "ggplot")

  # Iki imza; C atlandigi icin ortak veri seti sayisi dusuk kalabilir
  corr <- auc_correlation(res, method = "mean")
  expect_true(nrow(corr) >= 1L)
  expect_true(all(corr$n_datasets <= 2L))
})

test_that("auc_correlation 2 veri setinde NA, 3+ sonlu", {
  make_res <- function(n_ds) {
    ds <- paste0("D", seq_len(n_ds))
    data.frame(
      dataset = rep(ds, each = 2L),
      signature = rep(c("S1", "S2"), times = n_ds),
      method = "mean",
      auc = rep(c(0.6, 0.7), times = n_ds) + seq_len(2L * n_ds) * 0.01,
      skipped = FALSE,
      stringsAsFactors = FALSE
    )
  }

  corr2 <- auc_correlation(make_res(2L), method = "mean")
  expect_equal(nrow(corr2), 1L)
  expect_true(is.na(corr2$correlation[1]))
  expect_equal(corr2$n_datasets[1], 2L)

  corr3 <- auc_correlation(make_res(3L), method = "mean")
  expect_equal(nrow(corr3), 1L)
  expect_true(is.finite(corr3$correlation[1]))
  expect_equal(corr3$n_datasets[1], 3L)
  expect_gte(corr3$correlation[1], -1)
  expect_lte(corr3$correlation[1], 1)
})

test_that("generate_r_code hazir imzada kendi direction degerini korur, yalnizca Custom imzaya direction yazar", {
  code <- generate_r_code(
    signature_name = "Cytokine-5",
    direction = ">",
    methods = "mean",
    compute_pr = FALSE,
    compute_null = FALSE,
    dataset_names = "demo_A"
  )
  expect_silent(parse(text = code))
  expect_false(grepl('sigs[[1]]$direction <- ">"', code, fixed = TRUE))
  expect_true(grepl('direction = ">"', code, fixed = TRUE))
})

test_that("null_permutation_test hizli yol seed tekrarlanabilir (mean/zscore)", {
  skip_if_not_installed("withr")
  set.seed(11)
  mat <- matrix(rnorm(40 * 24), nrow = 40, ncol = 24,
                dimnames = list(paste0("G", 1:40), paste0("S", 1:24)))
  response <- c(rep(1L, 10), rep(0L, 14))
  names(response) <- colnames(mat)
  mat[1:4, 1:10] <- mat[1:4, 1:10] + 1.2

  for (meth in c("mean", "zscore", "median", "rank")) {
    a <- null_permutation_test(
      mat, response, n_genes = 4L,
      exclude_genes = paste0("G", 1:4),
      method = meth, direction = "<",
      n_perm = 25L, seed = 77L,
      signature_genes = paste0("G", 1:4),
      min_genes = 3L
    )
    b <- null_permutation_test(
      mat, response, n_genes = 4L,
      exclude_genes = paste0("G", 1:4),
      method = meth, direction = "<",
      n_perm = 25L, seed = 77L,
      signature_genes = paste0("G", 1:4),
      min_genes = 3L
    )
    expect_equal(a$null_aucs, b$null_aucs, info = meth)
    expect_equal(a$null_p, b$null_p, info = meth)
    expect_true(is.finite(a$obs_auc), info = meth)
  }
})

test_that("summarize_auc var_zero sayimi independent sonrasi", {
  res <- data.frame(
    dataset = c("A", "B", "C"),
    signature = "S1",
    method = "mean",
    n = c(20, 30, 25),
    auc = c(0.7, 0.8, 0.9),
    auc_var = c(0.01, 0.02, 0),
    var_zero = c(FALSE, FALSE, TRUE),
    independent = c(TRUE, TRUE, FALSE),
    small_cohort = FALSE,
    skipped = FALSE,
    stringsAsFactors = FALSE
  )
  # C hem non-independent hem var_zero: independent cikarildiktan sonra var_zero=0
  sm <- summarize_auc(res, method = "mean")
  expect_equal(sm$n_excluded_nonindependent, 1L)
  expect_equal(sm$n_excluded_var_zero, 0L)
  expect_equal(sm$n_datasets, 2L)
})
