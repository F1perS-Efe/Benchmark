# Faz 2 testleri: T_Cell_Inflamed_GEP agirlikli dogrusal model (2.1),
# Checkpoint_Axis vs B7_Immune_Exclusion ayrimi (2.2),
# Myeloid_Suppression CD68 ayrimi (2.3) ve ortak imza metadata standardi (2.4)

test_that("Faz 2.1: T_Cell_Inflamed_GEP carries exact Ayers/Cristescu 18-gene coefficients and matches reference linear predictor", {
  sigs <- get_builtin_signatures(include_controls = FALSE)
  gep <- sigs$T_Cell_Inflamed_GEP

  expected_w <- c(
    CCL5 = 0.008346,
    CD27 = 0.072293,
    CD274 = 0.042853,
    CD276 = -0.023900,
    CD8A = 0.031021,
    CMKLR1 = 0.151253,
    CXCL9 = 0.074135,
    CXCR6 = 0.004313,
    `HLA-DQA1` = 0.020091,
    `HLA-DRB1` = 0.058806,
    `HLA-E` = 0.071750,
    IDO1 = 0.060679,
    LAG3 = 0.123895,
    NKG7 = 0.075524,
    PDCD1LG2 = 0.003734,
    PSMB10 = 0.032999,
    STAT1 = 0.250229,
    TIGIT = 0.084767
  )

  expect_identical(gep$genes, names(expected_w))
  expect_equal(gep$weights, expected_w, tolerance = 1e-12)
  expect_identical(gep$direction, "signed")
  expect_identical(unname(gep$gene_directions["CD276"]), "-")
  expect_identical(unname(gep$gene_directions["STAT1"]), "+")
  expect_true(all(gep$gene_directions[names(expected_w) != "CD276"] == "+"))

  # Reference mathematical verification: exact inner product sum(w_i * X_{i,j})
  set.seed(2026)
  n_samp <- 24L
  mat <- matrix(
    stats::rnorm(18L * n_samp, mean = 6, sd = 1.5),
    nrow = 18L,
    ncol = n_samp,
    dimnames = list(names(expected_w), paste0("S", seq_len(n_samp)))
  )

  ref_scores <- as.numeric(crossprod(expected_w, mat))
  names(ref_scores) <- colnames(mat)

  sc_weighted <- score_signature(
    mat,
    genes = gep$genes,
    weights = gep$weights,
    method = "mean"
  )
  expect_true(sc_weighted$weighted)
  expect_true(sc_weighted$signed)
  expect_identical(sc_weighted$n_up_found, 17L)
  expect_identical(sc_weighted$n_down_found, 1L)
  expect_equal(sc_weighted$scores, ref_scores, tolerance = 1e-12)

  # Weighted GEP must differ from naive unweighted mean AND from naive equal-arm genes_up/genes_down
  sc_unweighted <- score_signature(mat, genes = gep$genes, method = "mean")
  sc_equal_arm <- score_signature(
    mat,
    genes_up = names(expected_w)[expected_w > 0],
    genes_down = "CD276",
    method = "mean"
  )
  expect_gt(max(abs(sc_weighted$scores - sc_unweighted$scores)), 0.1)
  expect_gt(max(abs(sc_weighted$scores - sc_equal_arm$scores)), 0.1)

  # Alias rescue preserves exact weight (query alias PDL1 -> matrix row CD274 gets 0.042853)
  w_alias <- gep$weights
  names(w_alias)[names(w_alias) == "CD274"] <- "PDL1"
  sc_alias <- score_signature(
    mat,
    genes = names(w_alias),
    weights = w_alias,
    method = "mean"
  )
  expect_identical(sc_alias$coverage$n_alias_rescued, 1L)
  expect_equal(sc_alias$scores, ref_scores, tolerance = 1e-12)
})

test_that("Faz 2.1 & 3: weighted scoring propagates across benchmark, score_for_dataset, LOO, manifest, and generate_r_code", {
  skip_if_not_installed("pROC")

  sigs <- get_builtin_signatures(include_controls = FALSE)
  gep <- sigs["T_Cell_Inflamed_GEP"]

  set.seed(42)
  n_samp <- 24L
  g_names <- gep$T_Cell_Inflamed_GEP$genes
  bg_names <- paste0("BG", seq_len(20L))
  all_rows <- c(g_names, bg_names)
  mat <- matrix(
    stats::rnorm(length(all_rows) * n_samp, mean = 5, sd = 1),
    nrow = length(all_rows),
    ncol = n_samp,
    dimnames = list(all_rows, paste0("S", seq_len(n_samp)))
  )
  # Inject responder signal aligned with positive weights and negative CD276
  mat[setdiff(g_names, "CD276"), 1:12] <- mat[setdiff(g_names, "CD276"), 1:12] + 1.5
  mat["CD276", 1:12] <- mat["CD276", 1:12] - 1.5

  ds <- list(synth_gep = list(
    expr = mat,
    sample_meta = data.frame(
      sample_id = colnames(mat),
      response_label = rep(c("CR", "PD"), each = 12L),
      stringsAsFactors = FALSE
    ),
    info = list(normalization = "TPM", scale = "log2")
  ))

  rb <- run_benchmark(gep, ds, methods = c("mean", "zscore", "ss_rank"), compute_pr = FALSE, compute_null = FALSE)
  expect_false(any(rb$skipped))
  expect_true(all(rb$auc > 0.85))

  s_ds <- score_for_dataset(ds$synth_gep, signature = gep$T_Cell_Inflamed_GEP, method = "mean")
  expect_equal(
    unname(s_ds$scores),
    as.numeric(crossprod(gep$T_Cell_Inflamed_GEP$weights, mat[g_names, , drop = FALSE])),
    tolerance = 1e-12
  )

  # LOO gene contributions use weights
  loo <- analyze_gene_contributions(ds, gep$T_Cell_Inflamed_GEP, method = "mean")
  expect_identical(nrow(loo$summary), 18L)

  # Reproducibility manifest fingerprint changes when a gene weight changes
  man1 <- build_reproducibility_manifest(rb, datasets = ds, signatures = gep)
  gep_mod <- gep
  gep_mod$T_Cell_Inflamed_GEP$weights["CD276"] <- 0.023900
  man2 <- build_reproducibility_manifest(rb, datasets = ds, signatures = gep_mod)
  expect_true(man1$signatures_fingerprint$T_Cell_Inflamed_GEP$weighted)
  expect_false(identical(
    man1$signatures_fingerprint$T_Cell_Inflamed_GEP$genes_md5,
    man2$signatures_fingerprint$T_Cell_Inflamed_GEP$genes_md5
  ))
  expect_false(identical(man1$reproducibility_hash, man2$reproducibility_hash))

  # generate_r_code supports custom continuous weights
  w_custom <- c(STAT1 = 0.25, CXCL9 = 0.15, CD276 = -0.05)
  code_w <- generate_r_code(
    signature_name = NULL,
    weights = w_custom,
    methods = "mean",
    dataset_names = "demo_A",
    compute_pr = FALSE,
    compute_null = FALSE
  )
  expect_silent(parse(text = code_w))
  expect_true(grepl("weights = c(\"STAT1\" = 0.25, \"CXCL9\" = 0.15, \"CD276\" = -0.05)", code_w, fixed = TRUE))
})

test_that("Faz 2.2, 2.3 & 2.4: Checkpoint_Axis, B7_Immune_Exclusion, Myeloid_Suppression, and metadata standard are scientifically valid", {
  sigs <- get_builtin_signatures(include_controls = TRUE)

  # 2.2: CD276 and VTCN1 separated from Checkpoint_Axis into B7_Immune_Exclusion
  expect_false("CD276" %in% sigs$Checkpoint_Axis$genes)
  expect_false("VTCN1" %in% sigs$Checkpoint_Axis$genes)
  expect_true(all(c("PDCD1", "CD274", "CTLA4", "LAG3", "TIGIT") %in% sigs$Checkpoint_Axis$genes))
  expect_identical(sigs$Checkpoint_Axis$direction, "positive")

  expect_true("B7_Immune_Exclusion" %in% names(sigs))
  expect_true(all(c("CD276", "VTCN1") %in% sigs$B7_Immune_Exclusion$genes))
  expect_identical(sigs$B7_Immune_Exclusion$direction, "negative")
  expect_identical(sigs$B7_Immune_Exclusion$min_genes, 2L)

  # 2.3: CD68 (pan-macrophage) and IDO1 excluded from Myeloid_Suppression
  expect_false("CD68" %in% sigs$Myeloid_Suppression$genes)
  expect_false("IDO1" %in% sigs$Myeloid_Suppression$genes)
  expect_true(all(c("IL10", "TGFB1", "ARG1", "CD163", "MRC1", "TREM2") %in% sigs$Myeloid_Suppression$genes))
  expect_identical(sigs$Myeloid_Suppression$direction, "negative")

  # 2.4: All built-in signatures carry complete standardized metadata fields
  req_fields <- c(
    "signature_name", "biological_interpretation", "genes", "gene_directions",
    "weights", "min_genes", "species", "direction", "category",
    "derived_from", "citation", "original_reference", "doi",
    "scoring_compatibility", "known_limitations"
  )
  for (nm in names(sigs)) {
    s <- sigs[[nm]]
    expect_true(all(req_fields %in% names(s)), info = paste("Missing metadata field in", nm))
    expect_identical(s$signature_name, nm)
    expect_true(nzchar(s$biological_interpretation), info = nm)
    expect_true(nzchar(s$original_reference), info = nm)
    expect_true(nzchar(s$doi), info = nm)
    expect_true(nzchar(s$known_limitations), info = nm)
    expect_identical(s$species, "Homo sapiens")
    expect_identical(names(s$gene_directions), s$genes)
    expect_true(all(s$gene_directions %in% c("+", "-")))
  }
})
