# Faz 3 testleri: Fare-insan ortolog haritasi politikalari (3.1),
# ss_rank siralama evreni butunlugu (3.2) ve null permutasyon arka plani (3.3)

test_that("Faz 3.1: ortholog mapping enforces 1:1 vs many:1 policy, rejects 1:many ambiguity, and reports loss/provenance", {
  curated <- SigBenchR:::.mouse_to_human_immune_map()
  expect_gte(length(curated), 150L)

  # All housekeeping, ribosomal, and geNorm negative control genes have 1:1 mouse orthologs
  sigs <- get_builtin_signatures(include_controls = TRUE)
  expect_true(all(sigs$NegCtrl_Housekeeping$genes %in% unname(curated)))
  expect_true(all(sigs$NegCtrl_Ribosomal$genes %in% unname(curated)))
  expect_true(all(sigs$NegCtrl_Reference_geNorm$genes %in% unname(curated)))
  expect_true(all(sigs$Cytolytic_Activity$genes %in% unname(curated)))
  expect_true(all(sigs$Checkpoint_Axis$genes %in% unname(curated)))
  expect_true(all(sigs$Myeloid_Suppression$genes %in% unname(curated)))

  # Custom ortholog table with 1:many (MouseAmb -> HUM1, HUM2) and many:1 (ParaA, ParaB -> HUM_SHARED)
  custom_orth <- data.frame(
    mouse_symbol = c("Ifng", "Stat1", "MouseAmb", "MouseAmb", "ParaA", "ParaB"),
    human_symbol = c("IFNG", "STAT1", "HUM1", "HUM2", "HUM_SHARED", "HUM_SHARED"),
    stringsAsFactors = FALSE
  )
  attr(custom_orth, "source_version") <- "MGI_Test_v1"
  attr(custom_orth, "annotation_date") <- "2026-03-01"

  mat <- matrix(
    c(
      10, 20,
      30, 40,
      50, 60,
      100, 200,
      300, 400
    ),
    nrow = 5L,
    byrow = TRUE,
    dimnames = list(c("Ifng", "Stat1", "MouseAmb", "ParaA", "ParaB"), c("S1", "S2"))
  )

  # 1) strict_1to1 drops both 1:many (MouseAmb) and many:1 paralogs (ParaA, ParaB) with warning
  expect_warning(
    m_strict <- map_mouse_to_human(
      mat,
      ortholog_table = custom_orth,
      ortholog_policy = "strict_1to1"
    ),
    "ambiguous mouse gene"
  )
  expect_identical(sort(rownames(m_strict)), c("IFNG", "STAT1"))
  st_strict <- attr(m_strict, "ortholog_stats")
  expect_identical(st_strict$n_input, 5L)
  expect_identical(st_strict$n_mapped, 2L)
  expect_identical(st_strict$n_ambiguous_dropped, 3L)
  expect_identical(st_strict$source_version, "MGI_Test_v1")
  expect_identical(st_strict$annotation_date, "2026-03-01")
  expect_equal(st_strict$retention_rate, 2 / 5)
  expect_equal(st_strict$loss_rate, 3 / 5)

  # 2) many_to_1 still drops 1:many (MouseAmb) with warning, but collapses ParaA & ParaB onto HUM_SHARED
  expect_warning(
    m_m21 <- map_mouse_to_human(
      mat,
      collapse_fun = "mean",
      ortholog_table = custom_orth,
      ortholog_policy = "many_to_1"
    ),
    "ambiguous mouse gene"
  )
  expect_identical(sort(rownames(m_m21)), c("HUM_SHARED", "IFNG", "STAT1"))
  expect_equal(as.numeric(m_m21["HUM_SHARED", ]), c(200, 300))
})

test_that("Faz 3.2 & 3.3: 18,000-gene mouse matrix preserves full background for ss_rank and null_permutation_test, while 49-gene partial matrix warns and flags", {
  skip_if_not_installed("pROC")

  set.seed(303)
  n_samp <- 16L
  y <- rep(c(1L, 0L), each = 8L)

  curated <- SigBenchR:::.mouse_to_human_immune_map()
  sig_mouse <- c("Ifng", "Cxcl9", "Cxcl10", "Stat1", "Ido1")
  sig_human <- c("IFNG", "CXCL9", "CXCL10", "STAT1", "IDO1")
  mouse_49 <- unique(c(sig_mouse, names(curated)))[seq_len(49L)]
  human_49 <- unname(curated[mouse_49])
  orth_49 <- data.frame(
    mouse_symbol = mouse_49,
    human_symbol = human_49,
    stringsAsFactors = FALSE
  )

  n_bg <- 18000L - 49L
  bg_mouse <- paste0("MouseBg_", seq_len(n_bg))
  all_mouse <- c(mouse_49, bg_mouse)

  # Construct 18,000 x 16 mouse transcriptome where all 49 immune genes are highly
  # expressed (~10) relative to the 17,951 background genes (~4)
  mat_18k <- matrix(
    stats::rnorm(18000L * n_samp, mean = 4, sd = 1),
    nrow = 18000L,
    ncol = n_samp,
    dimnames = list(all_mouse, paste0("S", seq_len(n_samp)))
  )
  mat_18k[mouse_49, ] <- stats::rnorm(49L * n_samp, mean = 10, sd = 0.2)
  # Add responder separation on Ifng, Cxcl9, Cxcl10, Stat1, Ido1
  sig_mouse <- c("Ifng", "Cxcl9", "Cxcl10", "Stat1", "Ido1")
  sig_human <- c("IFNG", "CXCL9", "CXCL10", "STAT1", "IDO1")
  mat_18k[sig_mouse, 1:8] <- mat_18k[sig_mouse, 1:8] + 1.5

  # Case A: Default keep_unmapped_context = TRUE preserves the 18,000-gene universe
  mapped_full_ctx <- map_mouse_to_human(
    mat_18k,
    ortholog_table = orth_49,
    keep_unmapped_context = TRUE
  )
  st_full <- attr(mapped_full_ctx, "ortholog_stats")
  expect_identical(nrow(mapped_full_ctx), 49L)
  expect_identical(st_full$n_input, 18000L)
  expect_identical(st_full$n_background_total, 18000L)
  expect_true(st_full$is_partial_mapping)
  expect_true(st_full$context_preserved)

  sc_ctx <- score_signature(mapped_full_ctx, genes = sig_human, method = "ss_rank", min_genes = 3L)
  expect_identical(sc_ctx$coverage$ss_rank_universe_size, 18000L)
  expect_false(sc_ctx$coverage$partial_ortholog_universe)
  # Because sig_human are at the top of the 18,000-gene transcriptome, ss_rank is near +0.50
  expect_true(all(sc_ctx$scores > 0.45))

  np_ctx <- null_permutation_test(
    mat = mapped_full_ctx,
    response = y,
    n_genes = 5L,
    exclude_genes = sig_human,
    method = "ss_rank",
    n_perm = 20L,
    seed = 42L,
    observed_scores = sc_ctx$scores
  )
  expect_identical(np_ctx$universe_size, 18000L)
  expect_identical(np_ctx$pool_size, 17995L)
  expect_false(np_ctx$partial_background)

  # Case B: 49-gene partial mapping WITHOUT preserved background context (Acceptance Criterion 3.3)
  mapped_partial_49 <- map_mouse_to_human(
    mat_18k,
    ortholog_table = orth_49,
    keep_unmapped_context = FALSE
  )
  expect_warning(
    sc_partial <- score_signature(mapped_partial_49, genes = sig_human, method = "ss_rank", min_genes = 3L),
    "ss_rank on partial ortholog-mapped matrix"
  )
  expect_identical(sc_partial$coverage$ss_rank_universe_size, 49L)
  expect_true(sc_partial$coverage$partial_ortholog_universe)
  # On the truncated 49-gene immune-only matrix, non-responder ss_rank collapses near 0 instead of > 0.45
  expect_lt(mean(sc_partial$scores[9:16]), 0.15)

  expect_warning(
    np_partial <- null_permutation_test(
      mat = mapped_partial_49,
      response = y,
      n_genes = 5L,
      exclude_genes = sig_human,
      method = "ss_rank",
      n_perm = 20L,
      seed = 42L,
      observed_scores = sc_partial$scores
    ),
    "Null permutation on partial ortholog-mapped matrix"
  )
  expect_identical(np_partial$universe_size, 49L)
  expect_identical(np_partial$pool_size, 44L)
  expect_true(np_partial$partial_background)

  # Case C: 2,800-gene bundled human cohort (datasets_ici$hugo_2016)
  data(datasets_ici, package = "SigBenchR")
  hugo <- datasets_ici$hugo_2016
  sc_hugo <- score_signature(hugo$expr, genes = sig_human, method = "ss_rank", min_genes = 3L)
  expect_identical(sc_hugo$coverage$ss_rank_universe_size, nrow(hugo$expr))
  expect_false(sc_hugo$coverage$partial_ortholog_universe)

  # Case D: End-to-end run_benchmark on mouse dataset preserves background context through .prepare_cohort()
  ds_mouse <- list(mouse_18k = list(
    expr = mat_18k,
    sample_meta = data.frame(
      sample_id = colnames(mat_18k),
      response_label = rep(c("CR", "PD"), each = 8L),
      stringsAsFactors = FALSE
    ),
    info = list(
      normalization = "TPM",
      scale = "log2",
      species = "mouse",
      ortholog_table = orth_49
    )
  ))
  prep_m <- SigBenchR:::.prepare_cohort(ds_mouse$mouse_18k, name = "mouse_18k")
  expect_identical(nrow(attr(prep_m$expr, "transcriptome_bg_expr")), 18000L)

  rb_m <- run_benchmark(
    signatures = list(IFNG5 = list(genes = sig_human, direction = "positive")),
    datasets = ds_mouse,
    methods = "ss_rank",
    compute_pr = FALSE,
    compute_null = TRUE,
    n_perm = 20L,
    seed = 42L
  )
  expect_false(rb_m$skipped[1])
  expect_true(is.finite(rb_m$null_p[1]))
})
