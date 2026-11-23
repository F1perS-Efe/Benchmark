# P2-F1 testleri: Gen duzeyinde leave-one-out delta-AUC + immun anotasyon

test_that("P2-F1: analyze_gene_contributions identifies driver vs noise genes and annotates immune axes", {
  skip_if_not_installed("ggplot2")
  set.seed(123)
  n <- 40L
  genes <- c("IFNG", "CXCL9", "NOISE1", "NOISE2", "VEGFA", paste0("BG", 1:10))
  mat <- matrix(rnorm(length(genes) * n), nrow = length(genes),
                dimnames = list(genes, paste0("S", seq_len(n))))
  # IFNG guclu sinyal tasiyicisi; VEGFA negatif kolda guclu sinyal; NOISE1/NOISE2 gurultu
  mat["IFNG", 1:20] <- mat["IFNG", 1:20] + 2.5
  mat["IFNG", 21:40] <- mat["IFNG", 21:40] - 2.5
  mat["VEGFA", 1:20] <- mat["VEGFA", 1:20] - 2.5
  mat["VEGFA", 21:40] <- mat["VEGFA", 21:40] + 2.5

  ds <- list(
    synth = list(
      expr = mat,
      sample_meta = data.frame(
        sample_id = colnames(mat),
        response_label = c(rep("CR", 20L), rep("PD", 20L)),
        stringsAsFactors = FALSE
      ),
      info = list(normalization = "TPM", scale = "log2")
    )
  )

  # 1) Tek yonlu imza: IFNG (sinyal) + NOISE1 + NOISE2
  sig_uni <- list(genes = c("IFNG", "NOISE1", "NOISE2"), direction = "positive")
  res_uni <- analyze_gene_contributions(ds, sig_uni, method = "zscore", sd_handling = "na")
  sm <- res_uni$summary

  expect_identical(nrow(sm), 3L)
  expect_true(sm$mean_delta_auc[sm$gene == "IFNG"] > 0)
  expect_true(all(sm$mean_delta_auc[sm$gene %in% c("NOISE1", "NOISE2")] <= 0))
  expect_identical(sm$immune_axis[sm$gene == "IFNG"], "IFN_gamma_JAK_STAT")
  expect_identical(sm$immune_axis[sm$gene == "NOISE1"], "Other")

  # 2) Isaretli (signed) imza: genes_up + genes_down
  sig_signed <- list(
    genes_up = c("IFNG", "NOISE1"),
    genes_down = c("VEGFA", "NOISE2"),
    direction = "signed"
  )
  res_signed <- analyze_gene_contributions(ds, sig_signed, method = "zscore")
  sm_s <- res_signed$summary
  expect_true(all(c("up", "down") %in% sm_s$arm))
  expect_identical(sm_s$immune_axis[sm_s$gene == "VEGFA"], "TGFb_Stromal")
  expect_true(sm_s$mean_delta_auc[sm_s$gene == "VEGFA"] > 0)

  # 3) ggplot ciktisi
  p <- plot_gene_contributions(res_signed)
  expect_s3_class(p, "ggplot")
})

test_that("P2-F1: Shiny gene_contrib_res computes LOO delta-AUC from snapshot", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("pROC")

  app_path <- system.file("shiny", "app.R", package = "SigBenchR")
  if (!nzchar(app_path)) {
    app_path <- file.path("inst", "shiny", "app.R")
  }
  env <- new.env(parent = globalenv())
  source(app_path, local = env)

  shiny::testServer(env$server, {
    session$setInputs(
      sig_choice = c("Cytokine-5", "IFNG6"),
      custom_genes = "",
      datasets = c("hugo_2016", "gide_2019"),
      methods = "mean",
      direction = "<",
      sd_handling = "exclude",
      coverage_threshold = 0.70,
      sensitivity_threshold = 0.10,
      compute_pr = FALSE,
      compute_null = FALSE,
      n_boot = 0L,
      n_perm = 10L,
      match_expression = FALSE,
      seed = 42L,
      heatmap_method = "mean",
      forest_sig = "IFNG6",
      gene_contrib_sig = "IFNG6"
    )
    session$setInputs(run = 1L)

    gc <- gene_contrib_res()
    expect_true(is.data.frame(gc$summary))
    expect_true(nrow(gc$summary) > 0L)
    expect_true("immune_axis" %in% names(gc$summary))
    expect_true("mean_delta_auc" %in% names(gc$summary))
  })
})
