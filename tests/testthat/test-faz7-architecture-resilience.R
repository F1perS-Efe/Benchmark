# FAZ 7 — Software Architecture, Input Validation, Error Handling, Performance & Resilience Tests

test_that("7.1: validate_dataset and prepare_expression reject invalid matrices, infinite values, all-NA, empty IDs, and duplicate samples", {
  base_mat <- matrix(c(1.2, 2.3, 3.4, 4.5), nrow = 2, dimnames = list(c("IFNG", "STAT1"), c("S1", "S2")))
  base_meta <- data.frame(sample_id = c("S1", "S2"), response_label = c("CR", "PD"), stringsAsFactors = FALSE)
  base_info <- list(normalization = "TPM", scale = "log2")

  # 1) All-NA matrix rejected
  ds_na <- list(
    expr = matrix(NA_real_, nrow = 2, ncol = 2, dimnames = list(c("G1", "G2"), c("S1", "S2"))),
    sample_meta = base_meta,
    info = base_info
  )
  expect_error(validate_dataset(ds_na), "completely NA")

  # 2) Infinite values rejected in both validate_dataset and prepare_expression
  ds_inf <- list(expr = base_mat, sample_meta = base_meta, info = base_info)
  ds_inf$expr[1, 1] <- Inf
  expect_error(validate_dataset(ds_inf), "infinite")
  expect_error(prepare_expression(ds_inf$expr, scale = "log2"), "infinite")

  # 3) Negative values in linear scale rejected by prepare_expression
  mat_neg <- base_mat
  mat_neg[1, 2] <- -0.5
  expect_error(prepare_expression(mat_neg, scale = "linear", apply_log2 = TRUE), "negative")

  # 4) Empty gene rownames or sample colnames rejected
  ds_empty_row <- list(expr = base_mat, sample_meta = base_meta, info = base_info)
  rownames(ds_empty_row$expr)[1] <- "  "
  expect_error(validate_dataset(ds_empty_row), "non-empty rownames")

  ds_dup_col <- list(
    expr = matrix(1:4, nrow = 2, dimnames = list(c("G1", "G2"), c("S1", "S1"))),
    sample_meta = data.frame(sample_id = c("S1", "S1"), response_label = c("CR", "PD"), stringsAsFactors = FALSE),
    info = base_info
  )
  expect_error(validate_dataset(ds_dup_col), "unique")

  # 5) Case-insensitive duplicate gene warning
  ds_case_dup <- list(
    expr = matrix(1:4, nrow = 2, dimnames = list(c("Ifng", "IFNG"), c("S1", "S2"))),
    sample_meta = base_meta,
    info = base_info
  )
  expect_warning(validate_dataset(ds_case_dup), "duplicate gene rownames")
})

test_that("7.2: Invalid signatures and malformed weight vectors fail with clear user-facing errors", {
  mat <- matrix(rnorm(20), nrow = 4, dimnames = list(c("IFNG", "STAT1", "CXCL9", "CD276"), paste0("S", 1:5)))

  # Empty gene list
  expect_error(score_signature(mat, genes = character(0)), "empty")
  # Both genes and genes_up provided
  expect_error(score_signature(mat, genes = "IFNG", genes_up = "STAT1"), "not both")
  # Overlapping genes_up and genes_down
  expect_error(score_signature(mat, genes_up = c("IFNG", "STAT1"), genes_down = c("ifng")), "overlap")
  # Unnamed weights
  expect_error(score_signature(mat, weights = c(0.5, -0.2)), "named numeric vector")
  # Zero or non-finite weights
  expect_error(score_signature(mat, weights = c(IFNG = 0, STAT1 = 0.5)), "finite and non-zero")
  expect_error(score_signature(mat, weights = c(IFNG = Inf, STAT1 = 0.5)), "finite and non-zero")
  # Duplicate gene names in weights
  expect_error(score_signature(mat, weights = c(IFNG = 0.4, ifng = 0.2)), "duplicate gene names")
})

test_that("7.3: Multi-cohort failure isolation keeps valid cohorts when one cohort fails in run_benchmark, beyond_infiltration, survival, and reference axes", {
  set.seed(703)
  genes <- c("IFNG", "CXCL9", "CXCL10", "STAT1", "IDO1", "HLA-DRA", "CD8A", "CD8B", "CD3D", "CD3E", paste0("BG", 1:20))
  n <- 24
  y <- c(rep("CR", 12), rep("PD", 12))
  mat_ok <- matrix(rnorm(length(genes) * n), nrow = length(genes), ncol = n,
                   dimnames = list(genes, paste0("S", seq_len(n))))
  ds_ok <- list(
    expr = mat_ok,
    sample_meta = data.frame(
      sample_id = colnames(mat_ok),
      response_label = y,
      os_time = rexp(n, rate = 0.1) + 0.5,
      os_event = sample(c(0L, 1L), n, replace = TRUE),
      stringsAsFactors = FALSE
    ),
    info = list(normalization = "TPM", scale = "log2")
  )

  # Corrupted cohort (contains Inf in expression matrix)
  ds_bad <- ds_ok
  ds_bad$expr[1, 1] <- Inf

  # Single broken cohort still stops immediately
  sigs <- list(IFNG6 = list(genes = c("IFNG", "CXCL9", "CXCL10", "STAT1", "IDO1", "HLA-DRA"), direction = "positive"))
  expect_error(run_benchmark(sigs, list(bad_only = ds_bad), methods = "mean", compute_pr = FALSE), "infinite")

  # Multi-cohort run isolates bad_cohort as skipped while preserving good_cohort
  expect_warning(
    res_multi <- run_benchmark(
      signatures = sigs,
      datasets = list(bad_cohort = ds_bad, good_cohort = ds_ok),
      methods = "mean",
      compute_pr = FALSE,
      compute_null = FALSE
    ),
    "failed validation/preparation and was marked skipped"
  )
  expect_equal(nrow(res_multi), 2L)
  row_bad <- res_multi[res_multi$dataset == "bad_cohort", , drop = FALSE]
  row_good <- res_multi[res_multi$dataset == "good_cohort", , drop = FALSE]
  expect_true(isTRUE(row_bad$skipped[1]))
  expect_match(row_bad$skip_reason[1], "Cohort preparation failed")
  expect_false(isTRUE(row_good$skipped[1]))
  expect_true(is.finite(row_good$auc[1]))

  # Multi-cohort beyond_infiltration isolates bad_cohort and returns good_cohort
  expect_warning(
    bi_multi <- test_beyond_infiltration(
      datasets = list(bad_cohort = ds_bad, good_cohort = ds_ok),
      signature = sigs$IFNG6,
      infiltration_axis = "CD8_T_Cell"
    ),
    "Skipping cohort 'bad_cohort'"
  )
  expect_equal(bi_multi$per_cohort$cohort, "good_cohort")

  # Multi-cohort survival benchmark isolates bad_cohort and returns good_cohort
  expect_warning(
    surv_multi <- run_survival_benchmark(
      datasets = list(bad_cohort = ds_bad, good_cohort = ds_ok),
      signatures = sigs
    ),
    "Skipping cohort 'bad_cohort'"
  )
  expect_equal(surv_multi$cohort, "good_cohort")

  # Multi-cohort correlate_reference_axes isolates bad_cohort and returns good_cohort
  expect_warning(
    cor_multi <- correlate_reference_axes(
      datasets = list(bad_cohort = ds_bad, good_cohort = ds_ok),
      signatures = sigs
    ),
    "Skipping cohort 'bad_cohort'"
  )
  expect_true(all(cor_multi$cohort == "good_cohort"))
})

test_that("7.4: null_permutation_test properly resolves continuous weights and handles single-class or zero-variance edge cases", {
  set.seed(704)
  n_g <- 80
  n_s <- 24
  y <- c(rep(1L, 12), rep(0L, 12))
  mat <- matrix(rnorm(n_g * n_s), nrow = n_g, ncol = n_s,
                dimnames = list(c("STAT1", "CXCL9", "IDO1", "CD276", paste0("BG", 1:76)), paste0("S", 1:n_s)))
  w <- c(STAT1 = 0.25, CXCL9 = 0.15, IDO1 = 0.10, CD276 = -0.05)

  np <- null_permutation_test(
    mat = mat,
    response = y,
    n_genes = 3L,
    n_genes_down = 1L,
    exclude_genes = names(w),
    method = "zscore",
    direction = "<",
    n_perm = 25L,
    seed = 42L,
    signature_genes = c("STAT1", "CXCL9", "IDO1"),
    signature_genes_down = "CD276",
    weights = w,
    min_genes = 2L
  )
  expect_true(is.finite(np$obs_auc))
  expect_equal(length(np$null_aucs), 25L)
  expect_true(is.finite(np$null_p))

  # Zero-variance genes in matrix are excluded from null pool without NaN contamination
  mat_zv <- mat
  mat_zv[5:10, ] <- 5.0
  np_zv <- null_permutation_test(
    mat = mat_zv,
    response = y,
    n_genes = 3L,
    n_genes_down = 1L,
    exclude_genes = names(w),
    method = "zscore",
    direction = "<",
    n_perm = 20L,
    seed = 42L,
    signature_genes = c("STAT1", "CXCL9", "IDO1"),
    signature_genes_down = "CD276",
    weights = w,
    min_genes = 2L
  )
  expect_equal(np_zv$pool_size, 80L - 4L - 6L)
  expect_true(is.finite(np_zv$null_p))
})

test_that("7.5: Unwritable export directory and missing optional package guards fail cleanly", {
  data(datasets_demo, package = "SigBenchR")
  data(signatures_demo, package = "SigBenchR")
  res <- suppressWarnings(run_benchmark(
    signatures_demo["IFNG6"],
    datasets_demo["demo_A"],
    methods = "mean",
    compute_pr = FALSE,
    compute_null = FALSE
  ))
  man <- build_reproducibility_manifest(res)

  ro_dir <- tempfile("ro_dir_")
  dir.create(ro_dir, mode = "0500")
  on.exit({
    Sys.chmod(ro_dir, mode = "0700")
    unlink(ro_dir, recursive = TRUE)
  }, add = TRUE)

  if (!identical(as.integer(file.access(ro_dir, 2)), 0L)) {
    expect_error(
      export_manifest(man, path = file.path(ro_dir, "man.json")),
      "not writable"
    )
    expect_error(
      export_publication_bundle(res, output_path = file.path(ro_dir, "bundle.zip")),
      "not writable"
    )
  }

  # Optional SummarizedExperiment guard on non-SE object
  expect_error(
    as_sigbenchr_dataset(12345),
    "SummarizedExperiment"
  )
})
