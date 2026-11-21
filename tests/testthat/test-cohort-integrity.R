# Bilimsel dogrulama ve kohort butunluk testleri (v0.1.1)

test_that("dort kohort validate_dataset gecer; expr sutunlari = sample_id tekil ve sirali", {
  data(datasets_ici, package = "SigBenchR")
  expect_identical(
    sort(names(datasets_ici)),
    sort(c("van_allen_2015", "hugo_2016", "riaz_2017", "gide_2019"))
  )

  for (nm in names(datasets_ici)) {
    ds <- datasets_ici[[nm]]
    expect_silent(validate_dataset(ds, name = nm))
    sid <- as.character(ds$sample_meta$sample_id)
    expect_identical(colnames(ds$expr), sid)
    expect_false(anyDuplicated(sid) > 0L)
    expect_false(anyDuplicated(colnames(ds$expr)) > 0L)
  }
})

test_that("datasets_ici info endpoint_type and treatment are populated", {
  data(datasets_ici, package = "SigBenchR")
  expect_identical(datasets_ici$van_allen_2015$info$endpoint_type, "clinical_benefit")
  expect_identical(datasets_ici$hugo_2016$info$endpoint_type, "RECIST")
  expect_identical(datasets_ici$riaz_2017$info$endpoint_type, "RECIST")
  expect_identical(datasets_ici$gide_2019$info$endpoint_type, "RECIST")
  expect_identical(datasets_ici$van_allen_2015$info$treatment, "anti-CTLA4")
  expect_identical(datasets_ici$hugo_2016$info$treatment, "anti-PD1")
  expect_identical(datasets_ici$riaz_2017$info$treatment, "anti-PD1")
  expect_identical(datasets_ici$gide_2019$info$treatment, "anti-PD1")
})

test_that("patient_id varsa hasta basina tek ornek (pseudo-replication = 0)", {
  data(datasets_ici, package = "SigBenchR")
  for (nm in names(datasets_ici)) {
    meta <- datasets_ici[[nm]]$sample_meta
    if ("patient_id" %in% names(meta)) {
      expect_false(
        anyDuplicated(as.character(meta$patient_id)) > 0L,
        info = nm
      )
    }
  }
})

test_that("beklenen n / n_pos / n_neg / n_SD ve exclude sonrasi efektif N", {
  data(datasets_ici, package = "SigBenchR")

  expected <- data.frame(
    cohort = c("van_allen_2015", "hugo_2016", "riaz_2017", "gide_2019"),
    n = c(42L, 26L, 49L, 41L),
    n_pos = c(19L, 14L, 10L, 19L),
    n_neg = c(23L, 12L, 23L, 16L),
    n_SD = c(0L, 0L, 16L, 6L),
    n_exclude = c(42L, 26L, 33L, 35L),
    stringsAsFactors = FALSE
  )

  for (i in seq_len(nrow(expected))) {
    nm <- expected$cohort[i]
    ds <- datasets_ici[[nm]]
    lab <- toupper(as.character(ds$sample_meta$response_label))
    expect_equal(as.integer(ds$info$n), expected$n[i], info = nm)
    expect_equal(as.integer(ds$info$n_pos), expected$n_pos[i], info = nm)
    expect_equal(as.integer(ds$info$n_neg), expected$n_neg[i], info = nm)
    expect_equal(as.integer(ds$info$n_excluded_sd), expected$n_SD[i], info = nm)
    expect_equal(sum(lab %in% c("CR", "PR")), expected$n_pos[i], info = nm)
    expect_equal(sum(lab == "PD"), expected$n_neg[i], info = nm)
    expect_equal(sum(lab == "SD"), expected$n_SD[i], info = nm)

    enc <- encode_response(ds$sample_meta$response_label, sd_handling = "exclude")
    expect_equal(sum(enc$kept), expected$n_exclude[i], info = nm)
    expect_equal(enc$n_excluded_sd, expected$n_SD[i], info = nm)
  }
})

test_that("sentetik matriste mean/median/zscore/rank kesin degerler; seed=42 tekrarlanabilir", {
  skip_if_not_installed("pROC")

  # mean / zscore oyuncak matrisi
  expr_mz <- matrix(
    c(1, 3, 5,
      2, 4, 6),
    nrow = 2L,
    byrow = TRUE,
    dimnames = list(c("G1", "G2"), c("S1", "S2", "S3"))
  )
  # mean: ornek ortalamalari
  res_mean <- score_signature(expr_mz, genes = c("G1", "G2"), method = "mean", min_genes = 2L)
  expect_equal(unname(res_mean$scores), c(1.5, 3.5, 5.5), tolerance = 1e-12)

  # zscore: her gen sd=2 ile z=(-1,0,1); ornek ortalamasi ayni
  res_z <- score_signature(expr_mz, genes = c("G1", "G2"), method = "zscore", min_genes = 2L)
  expect_equal(unname(res_z$scores), c(-1, 0, 1), tolerance = 1e-12)

  # median
  expr_med <- matrix(
    c(1, 3, 5,
      2, 4, 100),
    nrow = 2L,
    byrow = TRUE,
    dimnames = list(c("G1", "G2"), c("S1", "S2", "S3"))
  )
  res_med <- score_signature(expr_med, genes = c("G1", "G2"), method = "median", min_genes = 2L)
  expect_equal(unname(res_med$scores), c(1.5, 3.5, 52.5), tolerance = 1e-12)

  # rank: G1 -> (0, 0.5, 1); G2 sabit -> (0.5, 0.5, 0.5); ortalama (0.25, 0.5, 0.75)
  expr_rk <- matrix(
    c(1, 2, 3,
      5, 5, 5),
    nrow = 2L,
    byrow = TRUE,
    dimnames = list(c("G1", "G2"), c("S1", "S2", "S3"))
  )
  res_rk <- score_signature(expr_rk, genes = c("G1", "G2"), method = "rank", min_genes = 2L)
  expect_equal(unname(res_rk$scores), c(0.25, 0.5, 0.75), tolerance = 1e-12)

  # run_benchmark seed = 42L iki cagri birebir ayni
  set.seed(1L)
  n <- 12L
  g <- paste0("G", 1:5)
  mat <- matrix(
    rnorm(length(g) * n),
    nrow = length(g),
    dimnames = list(g, paste0("S", seq_len(n)))
  )
  ds <- list(
    synth = list(
      expr = mat,
      sample_meta = data.frame(
        sample_id = colnames(mat),
        response_label = rep(c("CR", "PD"), each = 6L),
        stringsAsFactors = FALSE
      ),
      info = list(normalization = "TPM", scale = "log2")
    )
  )
  sigs <- list(
    S = list(genes = g, direction = "<", derived_from = character(0))
  )
  args <- list(
    signatures = sigs,
    datasets = ds,
    methods = c("mean", "median", "zscore", "rank"),
    seed = 42L,
    compute_pr = FALSE,
    compute_null = FALSE,
    min_n = 5L,
    min_per_class = 2L
  )
  r1 <- do.call(run_benchmark, args)
  r2 <- do.call(run_benchmark, args)
  expect_equal(r1, r2, tolerance = 1e-12)
})

test_that("sd_handling politikalari riaz_2017 ve gide_2019 n / n_pos / n_neg degistirir", {
  skip_if_not_installed("pROC")
  data(datasets_ici, package = "SigBenchR")

  expect_policy <- function(nm, policy, n_exp, n_pos_exp, n_neg_exp) {
    ds <- datasets_ici[[nm]]
    enc <- encode_response(ds$sample_meta$response_label, sd_handling = policy)
    resp <- enc$response[enc$kept]
    expect_equal(length(resp), n_exp, info = paste(nm, policy, "n"))
    expect_equal(sum(resp == 1L), n_pos_exp, info = paste(nm, policy, "n_pos"))
    expect_equal(sum(resp == 0L), n_neg_exp, info = paste(nm, policy, "n_neg"))

    # run_benchmark n / n_pos / n_neg ile tutarli
    genes <- intersect(c("IFNG", "CXCL9", "CXCL10", "STAT1", "IDO1"), rownames(ds$expr))
    expect_gte(length(genes), 3L)
    res <- run_benchmark(
      signatures = list(T = list(genes = genes, direction = "<", derived_from = character(0))),
      datasets = datasets_ici[nm],
      methods = "mean",
      sd_handling = policy,
      compute_pr = FALSE,
      compute_null = FALSE,
      min_n = 5L,
      min_per_class = 2L
    )
    expect_equal(as.integer(res$n[1]), n_exp, info = paste(nm, policy, "rb_n"))
    expect_equal(as.integer(res$n_pos[1]), n_pos_exp, info = paste(nm, policy, "rb_pos"))
    expect_equal(as.integer(res$n_neg[1]), n_neg_exp, info = paste(nm, policy, "rb_neg"))
    if (identical(policy, "exclude")) {
      expect_equal(
        as.integer(res$n_excluded_sd[1]),
        as.integer(ds$info$n_excluded_sd),
        info = paste(nm, policy, "n_excluded_sd")
      )
    }
  }

  # riaz_2017: 10 PR + 23 PD + 16 SD
  expect_policy("riaz_2017", "exclude", 33L, 10L, 23L)
  expect_policy("riaz_2017", "as_0", 49L, 10L, 39L)
  expect_policy("riaz_2017", "as_1", 49L, 26L, 23L)

  # gide_2019: 19 PR + 16 PD + 6 SD
  expect_policy("gide_2019", "exclude", 35L, 19L, 16L)
  expect_policy("gide_2019", "as_0", 41L, 19L, 22L)
  expect_policy("gide_2019", "as_1", 41L, 25L, 16L)
})
