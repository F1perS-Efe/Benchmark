# P2-F4 testleri: Infiltrasyonun otesinde mi? (Beyond-infiltration testi)

test_that("P2-F4: test_beyond_infiltration distinguishes redundant infiltration proxy vs orthogonal signal", {
  set.seed(42)
  n <- 60L
  inf_genes <- c("CD3D", "CD3E", "CD8A", "CD8B")
  orth_genes <- c("ORTH1", "ORTH2", "ORTH3")
  proxy_genes <- c("PROXY1", "PROXY2", "PROXY3")
  bg_genes <- paste0("BG", 1:10)
  all_genes <- c(inf_genes, orth_genes, proxy_genes, bg_genes)

  # Baz infiltrasyon sinyali ve ona dik (ortogonal) gercek sinyal
  inf_latent <- rnorm(n)
  orth_latent <- rnorm(n)
  # Ortogonal sinyali infiltrasyonla iliskisiz hale getir
  orth_latent <- orth_latent - stats::coef(stats::lm(orth_latent ~ inf_latent))[2] * inf_latent
  orth_latent <- as.numeric(scale(orth_latent))

  lin_pred <- 0.8 * inf_latent + 1.8 * orth_latent
  prob <- stats::plogis(lin_pred)
  y <- stats::rbinom(n, 1L, prob)
  resp_lab <- ifelse(y == 1L, "CR", "PD")

  mat <- matrix(rnorm(length(all_genes) * n, sd = 0.4), nrow = length(all_genes),
                dimnames = list(all_genes, paste0("S", seq_len(n))))
  for (g in inf_genes) mat[g, ] <- mat[g, ] + inf_latent
  for (g in proxy_genes) {
    pure_noise <- stats::resid(stats::lm(rnorm(n) ~ y + inf_latent))
    mat[g, ] <- 0.95 * inf_latent + 0.2 * pure_noise
  }
  for (g in orth_genes) mat[g, ] <- mat[g, ] + orth_latent

  ds <- list(
    synth = list(
      expr = mat,
      sample_meta = data.frame(
        sample_id = colnames(mat),
        response_label = resp_lab,
        stringsAsFactors = FALSE
      ),
      info = list(normalization = "TPM", scale = "log2")
    )
  )

  # 1) Sadece infiltrasyonun gurultulu kopyasi -> bagimsiz ek deger yok
  res_proxy <- test_beyond_infiltration(
    datasets = ds,
    signature = list(genes = proxy_genes, direction = "positive"),
    infiltration_axis = "CD8_T_Cell",
    method = "zscore"
  )
  expect_identical(nrow(res_proxy$per_cohort), 1L)
  expect_false(isTRUE(res_proxy$per_cohort$adds_independent_value[1]))
  expect_gt(res_proxy$per_cohort$r_with_infiltration[1], 0.7)

  # 2) Infiltrasyona dik gercek ek sinyal -> LRT p < 0.05, delta > 0, adds_independent_value == TRUE
  res_orth <- test_beyond_infiltration(
    datasets = ds,
    signature = list(genes = orth_genes, direction = "positive"),
    infiltration_axis = "CD8_T_Cell",
    method = "zscore"
  )
  expect_identical(nrow(res_orth$per_cohort), 1L)
  expect_lt(res_orth$per_cohort$lrt_p[1], 0.05)
  expect_gt(res_orth$per_cohort$delta_auc_combined[1], 0)
  expect_true(isTRUE(res_orth$per_cohort$adds_independent_value[1]))
  expect_identical(res_orth$summary$n_independent[1], 1L)
})

test_that("P2-F4: Shiny beyond_inf_res evaluates from snapshot", {
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
      datasets = c("van_allen_2015", "gide_2019"),
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
      beyond_inf_sig = "IFNG6",
      beyond_inf_axis = "CD8_T_Cell"
    )
    session$setInputs(run = 1L)

    bi <- beyond_inf_res()
    expect_true(is.data.frame(bi$per_cohort))
    expect_true(is.data.frame(bi$summary))
    expect_true("lrt_p" %in% names(bi$per_cohort))
    expect_true("adds_independent_value" %in% names(bi$per_cohort))
  })
})
