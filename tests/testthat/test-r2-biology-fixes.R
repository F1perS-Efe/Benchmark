# R2 kanser immunolojisi dogrulugu, null kalibrasyonu ve ortolog guvenligi testleri

test_that("R2-1: builtin panel reflects exact Ayers/Cristescu GEP, Ayers18, Rooney 2-gene CYT, Myeloid without IDO1, and omits IMPRES_Like", {
  sigs <- get_builtin_signatures(include_controls = FALSE)

  expect_identical(
    sort(sigs$T_Cell_Inflamed_GEP$genes),
    sort(c(
      "CCL5", "CD27", "CD274", "CD276", "CD8A", "CMKLR1",
      "CXCL9", "CXCR6", "HLA-DQA1", "HLA-DRB1", "HLA-E", "IDO1",
      "LAG3", "NKG7", "PDCD1LG2", "PSMB10", "STAT1", "TIGIT"
    ))
  )
  expect_true(all(c("CIITA", "GZMK", "CXCR6", "LAG3") %in% sigs$Ayers_Expanded18$genes))
  expect_identical(length(sigs$Ayers_Expanded18$genes), 18L)

  expect_identical(sigs$Cytolytic_Activity$genes, c("GZMA", "PRF1"))
  expect_identical(sigs$Cytolytic_Activity$min_genes, 2L)

  expect_false("IDO1" %in% sigs$Myeloid_Suppression$genes)
  expect_false("IMPRES_Like" %in% names(sigs))

  # 2-genli Cytolytic_Activity run_benchmark(min_genes = 3L) ile skorlanabilmeli
  set.seed(11)
  n <- 20L
  mat <- matrix(rnorm(4 * n), nrow = 4, dimnames = list(c("GZMA", "PRF1", "G3", "G4"), paste0("S", 1:n)))
  mat[c("GZMA", "PRF1"), 1:10] <- mat[c("GZMA", "PRF1"), 1:10] + 2
  ds <- list(c1 = list(
    expr = mat,
    sample_meta = data.frame(sample_id = colnames(mat), response_label = rep(c("CR", "PD"), each = 10L), stringsAsFactors = FALSE),
    info = list(normalization = "TPM", scale = "log2")
  ))
  res <- run_benchmark(sigs["Cytolytic_Activity"], ds, methods = "mean", min_genes = 3L, compute_pr = FALSE, compute_null = FALSE)
  expect_false(isTRUE(res$skipped[1]))
  expect_gt(res$auc[1], 0.8)
})

test_that("R2-2: null calibration uses effective non-zero-SD count, match_expression uses aliases, pr_auc guards constants, bio adds BH, aliases remove CD3/CD8", {
  # 1) Alias temizligi: CD3 ve CD8 belirsiz kisaltmalar olmamali
  amap <- SigBenchR:::.ici_gene_alias_map()
  expect_false("CD3" %in% names(amap))
  expect_false("CD8" %in% names(amap))

  # 2) pr_auc sabit skorlarda uyari ile NA doner
  expect_warning(
    pr_c <- pr_auc(rep(1, 10), rep(c(1L, 0L), each = 5L), n_boot = 10L),
    "Constant|zero-variance"
  )
  expect_true(is.na(pr_c$pr_auc))

  # 3) null_permutation_test sifir-SD imza genleri dusunce efektif sayiyi kullanir ve match_expression alias kurtarir
  set.seed(5)
  n <- 20L
  pool <- paste0("P", 1:20)
  all_g <- c("CD274", "STAT1", "ZERO1", "ZERO2", pool)
  mat <- matrix(rnorm(length(all_g) * n), nrow = length(all_g), dimnames = list(all_g, paste0("S", 1:n)))
  mat[c("ZERO1", "ZERO2"), ] <- 1.0 # sifir SD
  y <- rep(c(1L, 0L), each = 10L)

  np <- null_permutation_test(
    mat = mat,
    response = y,
    n_genes = 4L,
    exclude_genes = c("CD274", "STAT1", "ZERO1", "ZERO2"),
    method = "zscore",
    n_perm = 25L,
    match_expression = TRUE,
    signature_genes = c("PDL1", "STAT1", "ZERO1", "ZERO2"), # PDL1 -> CD274 alias!
    observed_scores = stats::setNames(rnorm(n), colnames(mat))
  )
  expect_identical(np$n_genes_effective, 2L)
  expect_true(is.finite(np$null_p))

  # 4) correlate_reference_axes p_adj_bh sutunu dondurur
  data(datasets_ici, package = "SigBenchR")
  data(signatures_demo, package = "SigBenchR")
  bc <- correlate_reference_axes(datasets_ici["hugo_2016"], signatures_demo["IFNG6"], method = "mean")
  expect_true("p_adj_bh" %in% names(bc))
  expect_true(all(is.finite(bc$p_adj_bh)))
})

test_that("R2-3: map_mouse_to_human drops unvalidated MHC/unknown symbols by default", {
  mat <- matrix(
    1:12, nrow = 4L, ncol = 3L,
    dimnames = list(c("Pdcd1", "H2-K1", "H2-D1", "UnknownMouseGene"), c("S1", "S2", "S3"))
  )
  mapped <- map_mouse_to_human(mat)
  expect_identical(rownames(mapped), "PDCD1")
  st <- attr(mapped, "ortholog_stats")
  expect_identical(st$n_unmapped, 3L)

  # All unmapped should error clearly
  expect_error(
    map_mouse_to_human(mat[c("H2-K1", "UnknownMouseGene"), , drop = FALSE]),
    "No genes matched"
  )
})
