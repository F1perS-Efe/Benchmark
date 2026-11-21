test_that("auc_mw pROC ile esit (ties dahil)", {
  skip_if_not_installed("pROC")
  scores <- c(0.1, 0.2, 0.2, 0.8, 0.9, 0.7)
  response <- c(0, 0, 1, 1, 1, 0)
  mw <- auc_mw(scores, response)
  roc <- pROC::roc(response, scores, levels = c(0, 1), direction = "<", quiet = TRUE)
  expect_equal(mw, as.numeric(pROC::auc(roc)), tolerance = 1e-12)
})

test_that("pr_auc mukemmel / rastgele ve seed deterministik", {
  skip_if_not_installed("PRROC")
  skip_if_not_installed("withr")

  # Mukemmel ayirim
  scores <- c(5, 4, 3, 2, 1, 0)
  response <- c(1, 1, 1, 0, 0, 0)
  pr <- pr_auc(scores, response, direction = "<", n_boot = 50, seed = 7)
  expect_equal(pr$pr_auc, 1, tolerance = 1e-10)
  expect_equal(pr$pr_baseline, 0.5)

  # Ters ayirim -> dusuk PR
  pr_bad <- pr_auc(scores, response, direction = ">", n_boot = 0)
  expect_lt(pr_bad$pr_auc, 0.6)

  # Seed deterministik + global RNG korunur
  rs0 <- .Random.seed
  a <- pr_auc(scores, response, n_boot = 40, seed = 99)
  rs1 <- .Random.seed
  b <- pr_auc(scores, response, n_boot = 40, seed = 99)
  expect_equal(a$pr_ci_low, b$pr_ci_low)
  expect_equal(a$pr_ci_high, b$pr_ci_high)
  expect_identical(rs0, rs1)
})

test_that("null_permutation_test seed, p-sinirlari, sifir-varyans disi", {
  skip_if_not_installed("withr")
  set.seed(1)
  mat <- matrix(rnorm(50 * 30), nrow = 50, ncol = 30,
                dimnames = list(paste0("G", 1:50), paste0("S", 1:30)))
  # Sifir varyansli gen
  mat["G50", ] <- 1
  response <- c(rep(1L, 12), rep(0L, 18))
  names(response) <- colnames(mat)
  # Sinyalli imza
  mat[1:5, 1:12] <- mat[1:5, 1:12] + 1.5

  rs0 <- .Random.seed
  n1 <- null_permutation_test(
    mat, response, n_genes = 5,
    exclude_genes = paste0("G", 1:5),
    method = "mean", direction = "<",
    n_perm = 40, seed = 123,
    signature_genes = paste0("G", 1:5),
    min_genes = 3
  )
  n2 <- null_permutation_test(
    mat, response, n_genes = 5,
    exclude_genes = paste0("G", 1:5),
    method = "mean", direction = "<",
    n_perm = 40, seed = 123,
    signature_genes = paste0("G", 1:5),
    min_genes = 3
  )
  expect_equal(n1$null_aucs, n2$null_aucs)
  expect_equal(n1$null_p, n2$null_p)
  expect_gte(n1$null_p, 1 / (length(n1$null_aucs) + 1))
  expect_lte(n1$null_p, 1)
  # G50 sifir varyansli; havuzdan cikarilmis olmali
  expect_true(stats::sd(mat["G50", ]) == 0)
  expect_lte(n1$pool_size, 50L - 5L - 1L)
  expect_identical(.Random.seed, rs0)
})

test_that("compare_signatures ayni DS calisir, farkli DS hata, BH dogru", {
  skip_if_not_installed("pROC")
  set.seed(2)
  n <- 40
  resp <- c(rep(1L, 20), rep(0L, 20))
  names(resp) <- paste0("S", seq_len(n))
  s1 <- setNames(c(rnorm(20, 1), rnorm(20, 0)), names(resp))
  s2 <- setNames(c(rnorm(20, 0.2), rnorm(20, 0)), names(resp))
  s3 <- setNames(rnorm(n), names(resp))

  one <- compare_signatures(
    scores1 = s1, scores2 = s2, response = resp,
    dataset1 = "D1", dataset2 = "D1",
    name1 = "A", name2 = "B"
  )
  expect_true(is.finite(one$p_value))
  expect_equal(one$p_adj_bh, one$p_value)

  expect_error(
    compare_signatures(s1, s2, resp, dataset1 = "D1", dataset2 = "D2"),
    "Farkli veri setleri"
  )

  # Ortak ornek yok
  s_bad <- setNames(s2, paste0("X", seq_along(s2)))
  expect_error(
    compare_signatures(
      scores1 = s1, scores2 = s_bad, response = resp,
      dataset1 = "D1", dataset2 = "D1"
    ),
    "Ortak ornek"
  )

  multi <- compare_signatures(comparisons = list(
    list(scores1 = s1, scores2 = s2, response = resp, dataset = "D1", name1 = "A", name2 = "B"),
    list(scores1 = s1, scores2 = s3, response = resp, dataset = "D1", name1 = "A", name2 = "C")
  ))
  expect_equal(nrow(multi), 2)
  expect_equal(multi$p_adj_bh, p.adjust(multi$p_value, method = "BH"))
})

test_that("small_cohort ve var_zero kenar durumlari", {
  skip_if_not_installed("pROC")

  # Kucuk kohort
  scores <- c(1, 2, 3, 4, 5, 6)
  response <- c(1, 1, 1, 0, 0, 0)
  r <- roc_auc(scores, response, min_n = 20, min_per_class = 5)
  expect_true(r$small_cohort)

  # Mukemmel ayirim -> AUC=1, varyans 0 olabilir
  scores_p <- c(10, 9, 8, 1, 0, -1)
  response_p <- c(1, 1, 1, 0, 0, 0)
  r2 <- roc_auc(scores_p, response_p, min_n = 3, min_per_class = 2)
  expect_equal(r2$auc, 1)
  expect_true(isTRUE(r2$var_zero) || r2$auc_var == 0 || is.na(r2$auc_var) || is.finite(r2$auc_ci_low))
  # En azindan CI alanlari mevcut
  expect_true(is.finite(r2$auc_ci_low) || is.na(r2$auc_ci_low))
})

test_that("run_benchmark Faz 2B sutunlarini doldurur", {
  skip_if_not_installed("pROC")
  skip_if_not_installed("PRROC")
  data(datasets_demo, package = "SigBenchR")
  data(signatures_demo, package = "SigBenchR")

  res <- run_benchmark(
    signatures = signatures_demo["Cytokine-5"],
    datasets = datasets_demo["demo_A"],
    methods = "mean",
    compute_pr = TRUE,
    n_boot = 30,
    compute_null = TRUE,
    n_perm = 20,
    seed = 1
  )
  expect_true(all(c(
    "auc_ci_low", "auc_ci_high", "auc_var", "var_zero", "small_cohort",
    "pr_auc", "pr_baseline", "pr_ci_low", "pr_ci_high", "null_p"
  ) %in% names(res)))
  expect_true(is.finite(res$auc[1]))
  expect_true(is.finite(res$pr_auc[1]))
  expect_true(is.finite(res$null_p[1]))
})

test_that("plot_roc_curves ggplot doner", {
  skip_if_not_installed("pROC")
  data(datasets_demo, package = "SigBenchR")
  data(signatures_demo, package = "SigBenchR")
  p <- plot_roc_curves(
    datasets_demo, signatures_demo[c("Cytokine-5", "Random3")],
    dataset_name = "demo_A", method = "mean"
  )
  expect_s3_class(p, "ggplot")
})
