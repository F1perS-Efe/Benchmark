# P2-F3 testleri: Tam transkriptom kohort kayit defteri + fare -> insan ortolog esleme

test_that("P2-F3: map_mouse_to_human maps curated 1:1 orthologs, drops unmapped MHC/unknowns, and collapses duplicates", {
  mat <- matrix(
    c(
      2, 4, 6, 8,
      10, 12, 14, 16,
      1, 2, 3, 4,
      3, 4, 5, 6,
      5, 5, 5, 5
    ),
    nrow = 5L, ncol = 4L, byrow = TRUE,
    dimnames = list(c("H2-K1", "Pdcd1", "Ifng", "IFNG", ""), c("S1", "S2", "S3", "S4"))
  )

  res_mean <- map_mouse_to_human(mat, collapse_fun = "mean")
  expect_true(all(c("PDCD1", "IFNG") %in% rownames(res_mean)))
  expect_false("HLA-A" %in% rownames(res_mean))
  expect_identical(nrow(res_mean), 2L)
  expect_equal(unname(res_mean["IFNG", "S1"]), 2.0)

  st <- attr(res_mean, "ortholog_stats")
  expect_identical(st$n_input, 4L) # bos satir cikarildi
  expect_identical(st$n_mapped, 2L)
  expect_identical(st$n_unmapped, 1L) # H2-K1 1:1 degil -> dustu
  expect_identical(st$n_collapsed, 1L)

  res_max <- map_mouse_to_human(mat, collapse_fun = "max")
  expect_equal(unname(res_max["IFNG", "S1"]), 3.0)

  res_first <- map_mouse_to_human(mat, collapse_fun = "first")
  expect_equal(unname(res_first["IFNG", "S1"]), 1.0)
})

test_that("P2-F3: run_benchmark automatically maps mouse cohorts via species='mouse'", {
  skip_if_not_installed("pROC")
  set.seed(7)
  n <- 24L
  mouse_genes <- c("Ifng", "Cxcl9", "Cxcl10", "Stat1", "Ido1", "Pdcd1", "Gzmb")
  mat <- matrix(rnorm(length(mouse_genes) * n), nrow = length(mouse_genes),
                dimnames = list(mouse_genes, paste0("M", seq_len(n))))
  mat[c("Ifng", "Cxcl9", "Stat1"), 1:12] <- mat[c("Ifng", "Cxcl9", "Stat1"), 1:12] + 2.0

  ds_mouse <- list(
    syngeneic_mc38 = list(
      expr = mat,
      sample_meta = data.frame(
        sample_id = colnames(mat),
        response_label = c(rep("CR", 12L), rep("PD", 12L)),
        stringsAsFactors = FALSE
      ),
      info = list(normalization = "TPM", scale = "log2", species = "mouse", treatment = "anti-PD1")
    )
  )
  sigs <- list(
    IFNG6 = list(
      genes = c("IDO1", "CXCL9", "CXCL10", "HLA-DRA", "STAT1", "IFNG"),
      direction = "positive",
      derived_from = character(0)
    )
  )

  res <- run_benchmark(sigs, ds_mouse, methods = "mean", compute_pr = FALSE, compute_null = FALSE)
  expect_equal(res$coverage_pct[1], 100 * 5 / 6, tolerance = 1e-6)
  expect_false(isTRUE(res$skipped[1]))
  expect_gt(res$auc[1], 0.8)
})

test_that("P2-F3: cohort registry lists builtins, registers local rds/csv with MD5 verification", {
  reg_builtin <- list_cohort_registry()
  expect_true(all(c("van_allen_2015", "hugo_2016", "riaz_2017", "gide_2019") %in% reg_builtin$cohort_id))
  expect_true(all(reg_builtin$source_type == "builtin"))

  ds_b <- load_registered_cohort("hugo_2016")
  expect_identical(ds_b$info$species, "human")

  # Yerel RDS + CSV kayit ve MD5 kontrolu
  set.seed(11)
  mat <- matrix(rnorm(6 * 12), nrow = 6, dimnames = list(c("Ifng", "Stat1", "Ido1", "Cxcl9", "Cxcl10", "H2-K1"), paste0("S", 1:12)))
  meta <- data.frame(sample_id = colnames(mat), response_label = rep(c("CR", "PD"), each = 6), stringsAsFactors = FALSE)
  ds_loc <- list(expr = mat, sample_meta = meta, info = list(normalization = "TPM", scale = "log2"))

  rds_path <- file.path(tempdir(), "local_mouse_cohort.rds")
  saveRDS(ds_loc, rds_path)
  reg_path <- file.path(tempdir(), "test_cohort_registry.yml")
  if (file.exists(reg_path)) unlink(reg_path)

  row_rds <- register_cohort(
    cohort_id = "mc38_rds",
    file_path = rds_path,
    registry_path = reg_path,
    species = "mouse",
    cancer_type = "colorectal_syngeneic",
    treatment = "anti-PD-1",
    endpoint_type = "tumor_regression"
  )
  expect_identical(row_rds$treatment, "anti-PD1")
  expect_identical(row_rds$species, "mouse")
  expect_identical(row_rds$n_samples, 12L)

  loaded_rds <- load_registered_cohort("mc38_rds", registry_path = reg_path)
  expect_identical(loaded_rds$info$species, "mouse")

  # CSV kaydi
  csv_expr <- file.path(tempdir(), "local_human_expr.csv")
  csv_meta <- file.path(tempdir(), "local_human_meta.csv")
  write.csv(cbind(gene = rownames(mat), as.data.frame(mat)), csv_expr, row.names = FALSE)
  write.csv(meta, csv_meta, row.names = FALSE)
  register_cohort(
    cohort_id = "human_csv",
    file_path = csv_expr,
    meta_path = csv_meta,
    registry_path = reg_path,
    species = "human",
    treatment = "anti-CTLA-4"
  )
  loaded_csv <- load_registered_cohort("human_csv", registry_path = reg_path)
  expect_identical(loaded_csv$info$species, "human")

  # MD5 bozulursa hata vermeli
  saveRDS(list(expr = mat + 1, sample_meta = meta, info = ds_loc$info), rds_path)
  expect_error(load_registered_cohort("mc38_rds", registry_path = reg_path), "MD5 mismatch")
})

test_that("P2-F3: Shiny registry_load loads builtin cohort from registry", {
  skip_if_not_installed("shiny")
  app_path <- system.file("shiny", "app.R", package = "SigBenchR")
  if (!nzchar(app_path)) {
    app_path <- file.path("inst", "shiny", "app.R")
  }
  env <- new.env(parent = globalenv())
  source(app_path, local = env)

  shiny::testServer(env$server, {
    session$setInputs(datasets = "hugo_2016", registry_cohort = "gide_2019")
    session$setInputs(registry_load = 1L)
    expect_true(grepl("gide_2019", output$registry_status))
  })
})
