# P0 audit regression tests

test_that("P0-2 collapse_duplicate_genes all-NA duplicates do not crash", {
  expr <- matrix(
    c(
      NA_real_, NA_real_,
      NA_real_, NA_real_,
      1, 2,
      3, 4
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(c("GENE", "GENE", "A", "B"), c("S1", "S2"))
  )
  out <- collapse_duplicate_genes(expr)
  expect_equal(nrow(out), 3L)
  expect_true("GENE" %in% toupper(rownames(out)))
})

test_that("P0-2 validate_dataset rejects non-numeric and duplicate sample_id", {
  base <- list(
    expr = matrix(1:4, nrow = 2, dimnames = list(c("A", "B"), c("S1", "S2"))),
    sample_meta = data.frame(
      sample_id = c("S1", "S2"),
      response_label = c("CR", "PD"),
      stringsAsFactors = FALSE
    ),
    info = list(normalization = "TPM", scale = "log2")
  )
  expect_silent(validate_dataset(base))

  bad_chr <- base
  bad_chr$expr <- matrix(
    c("a", "b", "c", "d"),
    nrow = 2,
    dimnames = list(c("A", "B"), c("S1", "S2"))
  )
  expect_error(validate_dataset(bad_chr), "sayisal")

  bad_id <- base
  bad_id$sample_meta$sample_id <- c("S1", "S1")
  colnames(bad_id$expr) <- c("S1", "S1")
  expect_error(validate_dataset(bad_id), "tekil")
})

test_that("P0-2 encode_response warns on unknown labels", {
  expect_warning(
    enc <- encode_response(c("CR", "PD", "UNK")),
    "Unknown response label"
  )
  expect_true(is.na(enc$response[3]))
})

test_that("P0-3 SD exclusion does not leak into zscore/rank AUC", {
  skip_if_not_installed("pROC")
  set.seed(11)
  n_keep <- 20L
  n_sd <- 8L
  genes <- paste0("G", 1:6)
  samples <- c(paste0("K", seq_len(n_keep)), paste0("SD", seq_len(n_sd)))
  expr <- matrix(
    rnorm(length(genes) * length(samples)),
    nrow = length(genes),
    dimnames = list(genes, samples)
  )
  # Make signature somewhat informative on kept samples
  expr[1:3, 1:10] <- expr[1:3, 1:10] + 1.5
  labels <- c(rep(c("CR", "PD"), each = 10L), rep("SD", n_sd))
  ds_full <- list(
    expr = expr,
    sample_meta = data.frame(
      sample_id = samples,
      response_label = labels,
      stringsAsFactors = FALSE
    ),
    info = list(normalization = "TPM", scale = "log2")
  )
  ds_pre <- list(
    expr = expr[, seq_len(n_keep), drop = FALSE],
    sample_meta = data.frame(
      sample_id = samples[seq_len(n_keep)],
      response_label = labels[seq_len(n_keep)],
      stringsAsFactors = FALSE
    ),
    info = list(normalization = "TPM", scale = "log2")
  )
  sigs <- list(S = list(genes = genes, direction = "<", derived_from = character(0)))

  res_full <- run_benchmark(
    signatures = sigs,
    datasets = list(D = ds_full),
    methods = c("zscore", "rank"),
    sd_handling = "exclude",
    compute_pr = FALSE,
    compute_null = FALSE,
    min_n = 5L,
    min_per_class = 2L
  )
  res_pre <- run_benchmark(
    signatures = sigs,
    datasets = list(D = ds_pre),
    methods = c("zscore", "rank"),
    sd_handling = "exclude",
    compute_pr = FALSE,
    compute_null = FALSE,
    min_n = 5L,
    min_per_class = 2L
  )
  for (m in c("zscore", "rank")) {
    expect_equal(
      res_full$auc[res_full$method == m],
      res_pre$auc[res_pre$method == m],
      tolerance = 1e-12
    )
  }
})

test_that("P0-4 derived_from NULL is NA and case-insensitive match", {
  skip_if_not_installed("pROC")
  data(datasets_demo, package = "SigBenchR")
  ds <- datasets_demo["demo_A"]
  genes <- rownames(ds$demo_A$expr)[1:4]

  sig_null <- list(Custom = list(genes = genes, direction = "<")) # derived_from NULL
  res_null <- run_benchmark(
    signatures = sig_null,
    datasets = ds,
    methods = "mean",
    compute_pr = FALSE,
    compute_null = FALSE,
    min_n = 5L,
    min_per_class = 2L
  )
  expect_true(is.na(res_null$independent[1]))
  sm <- suppressWarnings(summarize_auc(res_null, independent_only = TRUE))
  expect_equal(sm$n_datasets[1], 0L)

  # Case-insensitive derived_from
  set.seed(2)
  n <- 12L
  g <- paste0("G", 1:5)
  expr <- matrix(rnorm(5 * n), nrow = 5, dimnames = list(g, paste0("S", seq_len(n))))
  ds_h <- list(
    hugo_2016 = list(
      expr = expr,
      sample_meta = data.frame(
        sample_id = colnames(expr),
        response_label = rep(c("CR", "PD"), each = 6L),
        stringsAsFactors = FALSE
      ),
      info = list(normalization = "FPKM", scale = "log2")
    )
  )
  sig_case <- list(
    S = list(genes = g, direction = "<", derived_from = "Hugo_2016")
  )
  res_case <- run_benchmark(
    signatures = sig_case,
    datasets = ds_h,
    methods = "mean",
    compute_pr = FALSE,
    compute_null = FALSE,
    min_n = 5L,
    min_per_class = 2L
  )
  expect_false(isTRUE(res_case$independent[1]))
})

test_that("P0-5 low_coverage summarize and heatmap threshold recompute", {
  res <- data.frame(
    dataset = c("A", "B"),
    signature = "S1",
    method = "mean",
    n = c(20, 30),
    auc = c(0.7, 0.8),
    auc_var = c(0.01, 0.02),
    var_zero = c(FALSE, FALSE),
    independent = c(TRUE, TRUE),
    small_cohort = c(FALSE, FALSE),
    skipped = FALSE,
    coverage_pct = c(60, 90),
    low_coverage = c(TRUE, FALSE),
    coverage_label = c("3/5 (60%)", "9/10 (90%)"),
    stringsAsFactors = FALSE
  )
  sm_ex <- summarize_auc(res, exclude_low_coverage = TRUE)
  expect_equal(sm_ex$n_datasets, 1L)
  expect_equal(sm_ex$n_low_coverage_excluded, 1L)
  sm_in <- summarize_auc(res, exclude_low_coverage = FALSE)
  expect_equal(sm_in$n_datasets, 2L)

  p <- plot_auc_heatmap(res, coverage_threshold = 0.50)
  expect_false(isTRUE(p$data$low_coverage[p$data$coverage_pct == 60][1]))
  expect_true(is.finite(p$data$auc_plot[p$data$coverage_pct == 60][1]))
})

test_that("P0-6 n_weighted includes var_zero; IVW excludes", {
  res <- data.frame(
    dataset = c("A", "B"),
    signature = "S1",
    method = "mean",
    n = c(20, 30),
    auc = c(1.0, 0.7),
    auc_var = c(0, 0.04),
    var_zero = c(TRUE, FALSE),
    independent = c(TRUE, TRUE),
    small_cohort = c(FALSE, FALSE),
    skipped = FALSE,
    stringsAsFactors = FALSE
  )
  expect_warning(
    sm <- summarize_auc(res, method = "mean"),
    "varyansi 0"
  )
  expect_equal(sm$auc_n_weighted, (20 * 1.0 + 30 * 0.7) / 50)
  expect_equal(sm$auc_ivw, 0.7)
  expect_equal(sm$n_cohorts_ivw, 1L)
  expect_equal(sm$n_var_zero_excluded, 1L)
})

test_that("P0-7 n_perm default >= 1000 and match_expression stable", {
  expect_gte(as.integer(formals(run_benchmark)$n_perm), 1000L)
  expect_gte(as.integer(formals(null_permutation_test)$n_perm), 1000L)

  set.seed(3)
  # Duplicate signature gene names + low-diversity edge; keep a large variable pool
  genes <- c("G1", "G1", "G2", paste0("P", 1:20))
  mat <- matrix(
    rnorm(length(genes) * 20),
    nrow = length(genes),
    dimnames = list(genes, paste0("S", 1:20))
  )
  mat["P1", ] <- 1 # zero-variance row filtered from pool
  y <- c(rep(1L, 10), rep(0L, 10))
  names(y) <- colnames(mat)
  out <- null_permutation_test(
    mat = mat,
    response = y,
    n_genes = 3L,
    exclude_genes = c("G1", "G2"),
    method = "mean",
    direction = "<",
    n_perm = 30L,
    seed = 7L,
    match_expression = TRUE,
    signature_genes = c("G1", "G1", "G2"),
    min_genes = 2L
  )
  expect_true(is.finite(out$null_p))
  expect_gte(out$null_p, 0)
  expect_lte(out$null_p, 1)
})

test_that("P0-8 compare_signatures mixed named/unnamed errors", {
  s1 <- setNames(rnorm(10), paste0("S", 1:10))
  s2 <- rnorm(10)
  y <- setNames(rep(c(1L, 0L), each = 5), paste0("S", 1:10))
  expect_error(
    compare_signatures(
      scores1 = s1, scores2 = s2, response = y,
      dataset1 = "D1", dataset2 = "D1"
    ),
    "Mixed named/unnamed"
  )
})

test_that("P0-9 generate_r_code escapes quotes and backslashes", {
  code <- generate_r_code(
    signature_name = NULL,
    genes = c('GENE"X', "A\\B"),
    methods = "mean",
    direction = "<",
    dataset_names = c("demo_A", "demo_B"),
    compute_pr = FALSE,
    compute_null = FALSE
  )
  expect_error(parse(text = code), NA)
  expect_true(grepl("summarize_auc\\(res\\)", code))
})

test_that("P0-10 shiny app parses and ROC direction '>' is respected", {
  skip_if_not_installed("pROC")
  app_path <- system.file("shiny", "app.R", package = "SigBenchR")
  if (!nzchar(app_path) || !file.exists(app_path)) {
    app_path <- file.path("inst", "shiny", "app.R")
  }
  expect_error(parse(file = app_path), NA)

  set.seed(5)
  n <- 30L
  genes <- paste0("G", 1:4)
  expr <- matrix(rnorm(4 * n), nrow = 4, dimnames = list(genes, paste0("S", seq_len(n))))
  # Higher score -> more likely PD when using raw mean; direction ">" should flip
  expr[1:2, 1:15] <- expr[1:2, 1:15] - 2
  expr[1:2, 16:30] <- expr[1:2, 16:30] + 2
  ds <- list(
    D1 = list(
      expr = expr,
      sample_meta = data.frame(
        sample_id = colnames(expr),
        response_label = c(rep("CR", 15), rep("PD", 15)),
        stringsAsFactors = FALSE
      ),
      info = list(normalization = "TPM", scale = "log2")
    )
  )
  sigs_gt <- list(S = list(genes = genes, direction = ">"))
  p <- plot_roc_curves(
    datasets = ds,
    signatures = sigs_gt,
    dataset_name = "D1",
    method = "mean",
    direction = "<",
    sd_handling = "exclude",
    min_genes = 2L
  )
  # Direction from signature list must be ">"; oriented AUC should be high
  scored <- score_for_dataset(ds$D1, genes = genes, method = "mean", direction = ">")
  roc_gt <- roc_auc(scored$scores, scored$response, direction = ">")
  roc_lt <- roc_auc(scored$scores, scored$response, direction = "<")
  expect_equal(roc_gt$auc + roc_lt$auc, 1, tolerance = 1e-8)
  expect_gt(roc_gt$auc, 0.5)
  expect_s3_class(p, "ggplot")
})
