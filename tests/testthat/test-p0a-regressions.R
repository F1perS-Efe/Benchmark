# P0-a regresyon testleri: P0-1, P0-3, P0-6, P0-9

test_that("P0-1: plot_roc_curves accepts builtin positive/negative/signed directions", {
  skip_if_not_installed("pROC")
  skip_if_not_installed("ggplot2")
  data(datasets_ici, package = "SigBenchR")

  b <- get_builtin_signatures(include_controls = FALSE)
  p <- plot_roc_curves(datasets_ici, b["IFNG6"], "hugo_2016")
  expect_s3_class(p, "ggplot")

  # Sentetik veri ile negative ve signed yonlerin run_benchmark ile eslestigini dogrula
  set.seed(123)
  genes <- c("G1", "G2", "G3", "G4", "U1", "U2", "U3", "D1", "D2", "D3")
  mat <- matrix(rnorm(length(genes) * 20L), nrow = length(genes),
                dimnames = list(genes, paste0("S", 1:20)))
  mat[c("G1", "G2", "G3", "G4"), 1:10] <- mat[c("G1", "G2", "G3", "G4"), 1:10] - 1.5
  mat[c("U1", "U2", "U3"), 1:10] <- mat[c("U1", "U2", "U3"), 1:10] + 1.5
  mat[c("D1", "D2", "D3"), 1:10] <- mat[c("D1", "D2", "D3"), 1:10] - 1.5
  ds <- list(synth = list(
    expr = mat,
    sample_meta = data.frame(
      sample_id = colnames(mat),
      response_label = c(rep("CR", 10L), rep("PD", 10L)),
      stringsAsFactors = FALSE
    ),
    info = list(normalization = "TPM", scale = "log2")
  ))
  sigs <- list(
    NegSig = list(genes = c("G1", "G2", "G3", "G4"), direction = "negative"),
    SignedSig = list(genes_up = c("U1", "U2", "U3"), genes_down = c("D1", "D2", "D3"), direction = "signed")
  )

  res <- run_benchmark(sigs, ds, methods = "mean", compute_pr = FALSE, min_n = 5L, min_per_class = 2L)
  p_s <- plot_roc_curves(ds, sigs, "synth", method = "mean")
  expect_s3_class(p_s, "ggplot")
  plot_auc <- unique(p_s$data[, c("signature", "auc")])
  for (i in seq_len(nrow(plot_auc))) {
    s <- plot_auc$signature[i]
    expect_equal(plot_auc$auc[i], res$auc[res$signature == s][1], tolerance = 1e-12)
  }
})

test_that("P0-3: run_benchmark cell seed is order-independent across signature panels", {
  skip_if_not_installed("pROC")
  data(datasets_ici, package = "SigBenchR")
  data(signatures_demo, package = "SigBenchR")

  ds <- datasets_ici["hugo_2016"]
  r1 <- suppressWarnings(run_benchmark(
    signatures_demo["IFNG6"], ds,
    methods = "mean", compute_pr = TRUE, n_boot = 50L,
    compute_null = TRUE, n_perm = 100L, seed = 42L
  ))
  r2 <- suppressWarnings(run_benchmark(
    signatures_demo[c("Cytokine-5", "IFNG6")], ds,
    methods = "mean", compute_pr = TRUE, n_boot = 50L,
    compute_null = TRUE, n_perm = 100L, seed = 42L
  ))

  row1 <- r1[r1$signature == "IFNG6", , drop = FALSE]
  row2 <- r2[r2$signature == "IFNG6", , drop = FALSE]
  expect_identical(row1$null_p, row2$null_p)
  expect_identical(row1$pr_ci_low, row2$pr_ci_low)
  expect_identical(row1$pr_ci_high, row2$pr_ci_high)
})

test_that("P0-6: meta summary formatter displays logit-scale I2 and tau2_reml alongside REML AUC", {
  sm <- data.frame(
    signature = "IFNG6",
    method = "mean",
    n_cohorts = 4L,
    auc_n_weighted = 0.665,
    auc_ivw = 0.686,
    auc_re = 0.659,
    auc_re_ci_low = 0.432,
    auc_re_ci_high = 0.831,
    I2 = 42.9,
    i2_logit = 33.5,
    tau2_reml = 0.093,
    note = NA_character_,
    stringsAsFactors = FALSE
  )

  expect_true(exists(".format_meta_summary", where = asNamespace("SigBenchR"), inherits = FALSE))
  fmt_fn <- get(".format_meta_summary", envir = asNamespace("SigBenchR"))
  tbl <- fmt_fn(sm)

  expect_true("I2 logit (%)" %in% names(tbl))
  expect_true("tau2 (REML)" %in% names(tbl))
  expect_identical(tbl[["I2 logit (%)"]], "33.5")
  expect_identical(tbl[["tau2 (REML)"]], "0.093")
  expect_false("I2 (%)" %in% names(tbl))
})

test_that("P0-9: Shiny DeLong comparison uses last run snapshot rather than live inputs", {
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
      datasets = "riaz_2017",
      methods = c("mean", "zscore"),
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
      forest_sig = "Cytokine-5"
    )
    session$setInputs(run = 1L)

    # Canli girdileri snapshot sonrasinda degistir (run tetiklemeden)
    session$setInputs(
      sd_handling = "as_0",
      direction = ">",
      cmp_dataset = "riaz_2017",
      cmp_sig1 = "Cytokine-5",
      cmp_sig2 = "IFNG6"
    )
    session$setInputs(run_cmp = 1L)

    cmp <- cmp_res()
    # riaz_2017 SD dislandiginda n = 33 (as_0 ile n = 49 olurdu)
    expect_identical(cmp$n, 33L)
    # direction = "<" snapshot'i korunmali (AUC > 0.5)
    expect_gt(cmp$auc1, 0.5)
  })
})
