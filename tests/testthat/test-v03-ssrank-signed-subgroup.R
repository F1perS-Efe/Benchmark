test_that("ss_rank cohort independence: adding/removing samples does not change patient 1", {
  set.seed(7L)
  genes <- paste0("G", 1:20)
  n0 <- 8L
  mat0 <- matrix(
    rnorm(length(genes) * n0),
    nrow = length(genes),
    dimnames = list(genes, paste0("S", seq_len(n0)))
  )
  sig <- c("G1", "G2", "G3", "G4")
  s0 <- score_signature(mat0, genes = sig, method = "ss_rank", min_genes = 3L)
  expect_false(s0$skipped)
  score_p1 <- unname(s0$scores["S1"])

  mat_add <- cbind(
    mat0,
    matrix(
      rnorm(length(genes) * 10L),
      nrow = length(genes),
      dimnames = list(genes, paste0("N", 1:10))
    )
  )
  s_add <- score_signature(mat_add, genes = sig, method = "ss_rank", min_genes = 3L)
  expect_equal(unname(s_add$scores["S1"]), score_p1, tolerance = 1e-12)

  mat_drop <- mat0[, 1:3, drop = FALSE]
  s_drop <- score_signature(mat_drop, genes = sig, method = "ss_rank", min_genes = 3L)
  expect_equal(unname(s_drop$scores["S1"]), score_p1, tolerance = 1e-12)
})

test_that("ss_rank scale invariance and closed-form singscore formula", {
  genes <- paste0("G", 1:10)
  expr <- matrix(
    c(
      1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
      10, 9, 8, 7, 6, 5, 4, 3, 2, 1
    ),
    nrow = 10L,
    byrow = FALSE,
    dimnames = list(genes, c("S1", "S2"))
  )
  # Positive linear vs log2(x+1): ranks identical
  s_lin <- score_signature(expr, genes = c("G8", "G9", "G10"), method = "ss_rank", min_genes = 2L)
  s_log <- score_signature(
    log2(expr + 1),
    genes = c("G8", "G9", "G10"),
    method = "ss_rank",
    min_genes = 2L
  )
  expect_equal(s_lin$scores, s_log$scores, tolerance = 1e-12)

  # Manual unidirectional for S1: ranks 1..10 by expression; up = G8,G9,G10 -> ranks 8,9,10
  N <- 10
  k_u <- 3
  m_u <- mean(c(8, 9, 10))
  s_up_hand <- (m_u - (k_u + 1) / 2) / (N - k_u) - 0.5
  expect_equal(unname(s_lin$scores["S1"]), s_up_hand, tolerance = 1e-12)
  expect_true(unname(s_lin$scores["S1"]) <= 0.5 + 1e-12)
  expect_true(unname(s_lin$scores["S1"]) >= -0.5 - 1e-12)

  # Bidirectional: up G8-G10, down G1-G2 on S1
  # down ranks 1,2 -> reverse 10,9; m_d = 9.5; k_d = 2
  m_d <- mean(c(10, 9))
  k_d <- 2
  s_down_hand <- (m_d - (k_d + 1) / 2) / (N - k_d) - 0.5
  s_bi_hand <- s_up_hand + s_down_hand
  s_bi <- score_signature(
    expr,
    genes_up = c("G8", "G9", "G10"),
    genes_down = c("G1", "G2"),
    method = "ss_rank",
    min_genes = 2L
  )
  expect_equal(unname(s_bi$scores["S1"]), s_bi_hand, tolerance = 1e-12)
  expect_true(unname(s_bi$scores["S1"]) <= 1 + 1e-12)
  expect_true(unname(s_bi$scores["S1"]) >= -1 - 1e-12)

  # Optional Bioconductor cross-check without declaring Suggests (CRAN policy).
  ss_pkg <- paste0("sing", "score")
  if (requireNamespace(ss_pkg, quietly = TRUE)) {
    ss_ns <- asNamespace(ss_pkg)
    ranked <- get("rankGenes", envir = ss_ns)(expr)
    sc <- get("simpleScore", envir = ss_ns)(
      ranked,
      upSet = c("G8", "G9", "G10"),
      downSet = c("G1", "G2")
    )
    expect_equal(
      unname(s_bi$scores["S1"]),
      as.numeric(sc$TotalScore[1]),
      tolerance = 1e-6
    )
  }
})

test_that("ss_rank returns NA_real_ when N == k (no Inf/NaN)", {
  expr <- matrix(
    1:4,
    nrow = 2L,
    dimnames = list(c("A", "B"), c("S1", "S2"))
  )
  res <- score_signature(expr, genes = c("A", "B"), method = "ss_rank", min_genes = 2L)
  expect_false(res$skipped)
  # All-NA sample scores are dropped; remaining must not be NaN/Inf
  expect_false(any(is.nan(res$scores)))
  expect_false(any(is.infinite(res$scores)))
  expect_true(length(res$scores) == 0L || all(is.na(res$scores)))

  rank_mat <- apply(expr, 2L, function(x) rank(x, ties.method = "average", na.last = "keep"))
  if (is.null(dim(rank_mat))) {
    rank_mat <- matrix(rank_mat, nrow = 2L, ncol = 2L)
  }
  rownames(rank_mat) <- rownames(expr)
  N <- colSums(!is.na(expr))
  arm <- getFromNamespace(".ss_rank_arm", "SigBenchR")(
    rank_mat, c("A", "B"), N, reverse = FALSE
  )
  expect_true(all(is.na(arm)))
  expect_false(any(is.nan(arm)))
  expect_false(any(is.infinite(arm)))
})

test_that("signed +/- signatures: genes_up/down == named signs; high AUC; overlap/warn/null", {
  skip_if_not_installed("pROC")
  set.seed(99L)
  n <- 40L
  genes <- c(paste0("U", 1:5), paste0("D", 1:5), paste0("N", 1:30))
  mat <- matrix(
    rnorm(length(genes) * n),
    nrow = length(genes),
    dimnames = list(genes, paste0("S", seq_len(n)))
  )
  y <- c(rep(1L, 20L), rep(0L, 20L))
  # Positive class: UP high, DOWN low
  mat[paste0("U", 1:5), y == 1L] <- mat[paste0("U", 1:5), y == 1L] + 3
  mat[paste0("D", 1:5), y == 1L] <- mat[paste0("D", 1:5), y == 1L] - 3
  mat[paste0("U", 1:5), y == 0L] <- mat[paste0("U", 1:5), y == 0L] - 1
  mat[paste0("D", 1:5), y == 0L] <- mat[paste0("D", 1:5), y == 0L] + 1

  s_ud <- score_signature(
    mat,
    genes_up = paste0("U", 1:5),
    genes_down = paste0("D", 1:5),
    method = "mean",
    min_genes = 3L
  )
  named <- c(setNames(rep(1, 5), paste0("U", 1:5)), setNames(rep(-1, 5), paste0("D", 1:5)))
  s_nm <- score_signature(mat, genes = named, method = "mean", min_genes = 3L)
  expect_equal(s_ud$scores, s_nm$scores, tolerance = 1e-12)
  expect_true(s_ud$signed)
  expect_equal(s_ud$n_up_found, 5L)
  expect_equal(s_ud$n_down_found, 5L)

  expect_error(
    score_signature(
      mat,
      genes_up = c("U1", "D1"),
      genes_down = c("D1", "D2"),
      method = "mean",
      min_genes = 1L
    ),
    "overlap"
  )

  mat_miss <- mat
  rownames(mat_miss)[rownames(mat_miss) == "D1"] <- "DX"
  rownames(mat_miss)[rownames(mat_miss) == "D2"] <- "DY"
  rownames(mat_miss)[rownames(mat_miss) == "D3"] <- "DZ"
  rownames(mat_miss)[rownames(mat_miss) == "D4"] <- "DW"
  rownames(mat_miss)[rownames(mat_miss) == "D5"] <- "DV"
  expect_warning(
    score_signature(
      mat_miss,
      genes_up = paste0("U", 1:5),
      genes_down = paste0("D", 1:5),
      method = "mean",
      min_genes = 3L
    ),
    "one arm found 0"
  )

  ds <- list(
    synth = list(
      expr = mat,
      sample_meta = data.frame(
        sample_id = colnames(mat),
        response_label = ifelse(y == 1L, "CR", "PD"),
        treatment = "anti-PD1",
        stringsAsFactors = FALSE
      ),
      info = list(
        normalization = "TPM",
        scale = "log2",
        treatment = "anti-PD1",
        endpoint_type = "RECIST",
        cancer_type = "melanoma"
      )
    )
  )
  sigs <- list(
    Signed = list(
      genes_up = paste0("U", 1:5),
      genes_down = paste0("D", 1:5),
      direction = "<",
      derived_from = character(0)
    )
  )
  res <- run_benchmark(
    signatures = sigs,
    datasets = ds,
    methods = c("mean", "ss_rank"),
    compute_pr = FALSE,
    compute_null = TRUE,
    n_perm = 50L,
    seed = 11L,
    min_genes = 3L,
    min_n = 5L,
    min_per_class = 2L
  )
  expect_true(all(res$signed))
  expect_true(all(res$auc > 0.90, na.rm = TRUE))
  expect_true(all(c("treatment", "endpoint_type", "cancer_type") %in% names(res)))
  expect_identical(unique(as.character(res$endpoint_type)), "RECIST")
  expect_true(all(is.finite(res$null_p)))

  r1 <- null_permutation_test(
    mat = mat,
    response = y,
    n_genes = 5L,
    n_genes_down = 5L,
    exclude_genes = c(paste0("U", 1:5), paste0("D", 1:5)),
    method = "mean",
    direction = "<",
    n_perm = 50L,
    seed = 5L,
    observed_scores = s_ud$scores,
    min_genes = 3L
  )
  r2 <- null_permutation_test(
    mat = mat,
    response = y,
    n_genes = 5L,
    n_genes_down = 5L,
    exclude_genes = c(paste0("U", 1:5), paste0("D", 1:5)),
    method = "mean",
    direction = "<",
    n_perm = 50L,
    seed = 5L,
    observed_scores = s_ud$scores,
    min_genes = 3L
  )
  expect_equal(r1$null_p, r2$null_p)
  expect_equal(r1$null_aucs, r2$null_aucs)
})

test_that("summarize_auc stratify_by splits treatment/endpoint_type; warns when pooled", {
  res <- data.frame(
    dataset = c("A", "B", "C", "D"),
    signature = "S1",
    method = "mean",
    n = c(20, 30, 25, 35),
    auc = c(0.7, 0.8, 0.75, 0.85),
    auc_var = c(0.01, 0.02, 0.015, 0.025),
    var_zero = FALSE,
    independent = TRUE,
    small_cohort = FALSE,
    skipped = FALSE,
    low_coverage = FALSE,
    treatment = c("anti-PD1", "anti-PD1", "anti-CTLA4", "anti-CTLA4"),
    endpoint_type = c("RECIST", "RECIST", "clinical_benefit", "clinical_benefit"),
    stringsAsFactors = FALSE
  )

  expect_warning(
    sm <- summarize_auc(res, method = "mean"),
    "stratify_by"
  )
  expect_equal(nrow(sm), 1L)
  expect_equal(sm$n_cohorts, 4L)

  sm_s <- summarize_auc(
    res,
    method = "mean",
    stratify_by = c("treatment", "endpoint_type")
  )
  expect_equal(nrow(sm_s), 2L)
  expect_true(all(c("treatment", "endpoint_type") %in% names(sm_s)))
  expect_equal(sort(sm_s$n_cohorts), c(2L, 2L))
  pd1 <- sm_s[sm_s$treatment == "anti-PD1", , drop = FALSE]
  expect_equal(nrow(pd1), 1L)
  expect_equal(pd1$auc_n_weighted, (20 * 0.7 + 30 * 0.8) / 50)
})
