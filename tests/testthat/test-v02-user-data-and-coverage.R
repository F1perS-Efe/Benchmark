# SigBenchR v0.2.0: local user dataset contract, bundle loader, coverage detail

.make_tiny_files <- function(dir,
                             genes_in_rows = TRUE,
                             scale_vals = NULL,
                             dup_gene = FALSE,
                             bad_numeric = FALSE,
                             response_extra = NULL) {
  genes <- c("IFNG", "CXCL9", "STAT1", "IDO1")
  if (isTRUE(dup_gene)) {
    genes <- c("IFNG", "IFNG", "CXCL9", "STAT1")
  }
  samples <- c("S1", "S2", "S3", "S4")
  if (is.null(scale_vals)) {
    mat <- matrix(
      c(
        1.1, 2.2, 3.3, 4.4,
        2.0, 2.5, 3.0, 3.5,
        0.5, 1.5, 2.5, 3.5,
        4.0, 3.0, 2.0, 1.0
      ),
      nrow = 4L, byrow = TRUE,
      dimnames = list(genes, samples)
    )
  } else {
    mat <- scale_vals
  }
  if (isTRUE(bad_numeric)) {
    mat[1, 1] <- NA
    # write as character "x" for non-numeric
    expr_path <- file.path(dir, "expr_bad.csv")
    df <- as.data.frame(mat, check.names = FALSE)
    df[1, 1] <- "not_a_number"
    utils::write.csv(cbind(gene = rownames(mat), df), expr_path, row.names = FALSE)
  } else if (isTRUE(genes_in_rows)) {
    expr_path <- file.path(dir, "expr.csv")
    utils::write.csv(cbind(gene = rownames(mat), as.data.frame(mat, check.names = FALSE)),
                     expr_path, row.names = FALSE)
  } else {
    # samples in rows
    expr_path <- file.path(dir, "expr.tsv")
    tmat <- t(mat)
    utils::write.table(
      cbind(sample = rownames(tmat), as.data.frame(tmat, check.names = FALSE)),
      expr_path,
      sep = "\t", row.names = FALSE, quote = FALSE
    )
  }

  resp <- c("CR", "PR", "PD", "PD")
  if (!is.null(response_extra)) {
    resp[4] <- response_extra
  }
  meta <- data.frame(
    sample_id = samples,
    response_label = resp,
    stringsAsFactors = FALSE
  )
  meta_path <- file.path(dir, "meta.csv")
  utils::write.csv(meta, meta_path, row.names = FALSE)
  list(expr_path = expr_path, meta_path = meta_path, genes = genes, samples = samples, mat = mat)
}

test_that("read_user_dataset basic CSV and transposed TSV load", {
  td <- tempfile("v02_basic_")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)

  f1 <- .make_tiny_files(td, genes_in_rows = TRUE)
  ds1 <- read_user_dataset(
    expr_path = f1$expr_path,
    meta_path = f1$meta_path,
    genes_in_rows = TRUE,
    gene_id_type = "HGNC",
    normalization = "TPM",
    scale = "log2"
  )
  expect_silent(validate_dataset(ds1))
  expect_equal(ds1$info$normalization, "TPM")
  expect_equal(ds1$info$scale, "log2")
  expect_equal(ds1$info$source, "user_local")
  expect_equal(ds1$info$n_genes, 4L)
  expect_equal(ds1$info$n_samples, 4L)

  td2 <- tempfile("v02_t_")
  dir.create(td2)
  on.exit(unlink(td2, recursive = TRUE), add = TRUE)
  f2 <- .make_tiny_files(td2, genes_in_rows = FALSE)
  ds2 <- read_user_dataset(
    expr_path = f2$expr_path,
    meta_path = f2$meta_path,
    genes_in_rows = FALSE,
    gene_id_type = "HGNC",
    normalization = "FPKM",
    scale = "linear"
  )
  expect_silent(validate_dataset(ds2))
  expect_equal(ds2$info$normalization, "FPKM")
  expect_equal(ds2$info$scale, "linear")
  expect_equal(nrow(ds2$expr), 4L)
  expect_equal(ncol(ds2$expr), 4L)
})

test_that("double log2 protection: scale=log2 unchanged; linear gets log2", {
  td <- tempfile("v02_log2_")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  f <- .make_tiny_files(td)

  ds_log <- read_user_dataset(
    f$expr_path, f$meta_path,
    scale = "log2", normalization = "TPM"
  )
  prepared_log <- prepare_expression(ds_log$expr, scale = ds_log$info$scale, apply_log2 = TRUE)
  expect_equal(prepared_log, ds_log$expr)

  ds_lin <- read_user_dataset(
    f$expr_path, f$meta_path,
    scale = "linear", normalization = "TPM"
  )
  prepared_lin <- prepare_expression(ds_lin$expr, scale = ds_lin$info$scale, apply_log2 = TRUE)
  expect_equal(prepared_lin, log2(ds_lin$expr + 1))
  expect_false(isTRUE(all.equal(prepared_lin, ds_lin$expr)))
})

test_that("edge cases: dup genes, Ensembl warn, unknown label, non-numeric, sample mismatch", {
  td <- tempfile("v02_edge_")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)

  # Duplicate genes
  f_dup <- .make_tiny_files(td, dup_gene = TRUE)
  expect_warning(
    ds_dup <- read_user_dataset(
      f_dup$expr_path, f_dup$meta_path,
      collapse_duplicates = TRUE, scale = "log2"
    ),
    "yinelenen gen"
  )
  expect_equal(nrow(ds_dup$expr), 3L)

  # Ensembl warning
  dir.create(file.path(td, "ens"))
  f_ens <- .make_tiny_files(file.path(td, "ens"))
  expect_warning(
    read_user_dataset(
      f_ens$expr_path, f_ens$meta_path,
      gene_id_type = "Ensembl", scale = "log2"
    ),
    "Ensembl"
  )

  # Unknown response label
  dir.create(file.path(td, "unk"))
  f_unk <- .make_tiny_files(file.path(td, "unk"), response_extra = "UNK")
  expect_warning(
    read_user_dataset(f_unk$expr_path, f_unk$meta_path, scale = "log2"),
    "Unknown response label"
  )

  # Non-numeric matrix
  dir.create(file.path(td, "bad"))
  f_bad <- .make_tiny_files(file.path(td, "bad"), bad_numeric = TRUE)
  expect_error(
    read_user_dataset(f_bad$expr_path, f_bad$meta_path, scale = "log2"),
    "sayisal"
  )

  # No overlapping sample_id
  td_miss <- file.path(td, "miss")
  dir.create(td_miss)
  f_ok <- .make_tiny_files(td_miss)
  meta_bad <- data.frame(
    sample_id = c("X1", "X2", "X3", "X4"),
    response_label = c("CR", "PR", "PD", "PD"),
    stringsAsFactors = FALSE
  )
  meta_bad_path <- file.path(td_miss, "meta_bad.csv")
  utils::write.csv(meta_bad, meta_bad_path, row.names = FALSE)
  expect_error(
    read_user_dataset(f_ok$expr_path, meta_bad_path, scale = "log2"),
    "ortak sample_id"
  )
})

test_that("load_cohort_bundle MD5 integrity (no network)", {
  td <- tempfile("v02_bundle_")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)

  set.seed(1)
  mat <- matrix(rnorm(12), nrow = 3, dimnames = list(c("A", "B", "C"), c("S1", "S2", "S3", "S4")))
  ds <- list(
    expr = mat,
    sample_meta = data.frame(
      sample_id = colnames(mat),
      response_label = c("CR", "PR", "PD", "PD"),
      stringsAsFactors = FALSE
    ),
    info = list(normalization = "TPM", scale = "log2")
  )
  bundle <- list(toy = ds)
  rds_path <- file.path(td, "bundle.rds")
  saveRDS(bundle, rds_path)
  good_md5 <- unname(tools::md5sum(rds_path))

  loaded <- load_cohort_bundle(rds_path, expected_md5 = good_md5)
  expect_true("toy" %in% names(loaded))
  expect_silent(validate_dataset(loaded$toy))

  expect_error(
    load_cohort_bundle(rds_path, expected_md5 = "00000000000000000000000000000000"),
    "MD5"
  )
})

test_that("missing_genes detail and 10 published signatures coverage", {
  skip_if_not_installed("pROC")

  # Synthetic matrix with intentional missing gene for one signature
  genes_full <- c(
    "IFNG", "CXCL9", "CXCL10", "CXCL11", "CXCL13", "CCL5", "IL21",
    "IDO1", "STAT1", "HLA-DRA", "HLA-E", "CD8A", "GZMA", "GZMB", "PRF1",
    "PDCD1", "CD274", "CTLA4", "LAG3", "TIGIT", "HAVCR2", "CD27", "CD28",
    "B2M", "TAP1", "IRF1", "NKG7", "CMKLR1", "PSMB8", "PSMB9",
    "CD80", "CD86", "HLA-A", "HLA-B"
  )
  set.seed(42)
  n <- 20L
  mat <- matrix(
    rnorm(length(genes_full) * n, mean = 5, sd = 1),
    nrow = length(genes_full),
    dimnames = list(genes_full, paste0("S", seq_len(n)))
  )
  # Drop MISSING_GENE from matrix so coverage reports it
  ds <- list(
    full = list(
      expr = mat,
      sample_meta = data.frame(
        sample_id = colnames(mat),
        response_label = rep(c("CR", "PD"), each = 10L),
        stringsAsFactors = FALSE
      ),
      info = list(normalization = "TPM", scale = "log2")
    )
  )

  sig_missing <- list(
    Partial = list(
      genes = c("IFNG", "CXCL9", "MISSING_GENE_XYZ"),
      direction = "<",
      derived_from = character(0)
    )
  )
  res_m <- run_benchmark(
    signatures = sig_missing,
    datasets = ds,
    methods = "mean",
    compute_pr = FALSE,
    compute_null = FALSE,
    min_genes = 2L,
    min_n = 5L,
    min_per_class = 2L
  )
  expect_true("found_genes" %in% names(res_m))
  expect_true("missing_genes" %in% names(res_m))
  expect_true(grepl("MISSING_GENE_XYZ", res_m$missing_genes[1], fixed = TRUE))
  expect_true(grepl("IFNG", res_m$found_genes[1], fixed = TRUE))

  # 10 published immuno-oncology signatures on full synthetic matrix (>= 85%)
  sigs10 <- list(
    IFNG6 = list(genes = c("IDO1", "CXCL9", "CXCL10", "HLA-DRA", "STAT1", "IFNG")),
    `Cytokine-5` = list(genes = c("IFNG", "CXCL9", "CXCL13", "IL21", "CCL5")),
    `Ayers-GEP-core` = list(genes = c("CXCL9", "CXCL10", "HLA-DRA", "IDO1", "STAT1", "IFNG", "CCL5", "CD8A")),
    `IMPRES-subset` = list(genes = c("PDCD1", "CD274", "CTLA4", "CD28", "CD80", "CD86")),
    Cytolytic = list(genes = c("GZMA", "PRF1", "GZMB")),
    Chemokine = list(genes = c("CXCL9", "CXCL10", "CCL5")),
    Exhaustion = list(genes = c("PDCD1", "LAG3", "HAVCR2", "TIGIT")),
    AntigenPresentation = list(genes = c("B2M", "HLA-A", "HLA-B", "TAP1")),
    InflamedGEP = list(genes = c("CXCL9", "CXCL10", "CXCL11", "STAT1", "IDO1", "HLA-E", "CMKLR1", "NKG7")),
    Interferon = list(genes = c("IFNG", "STAT1", "IRF1", "IDO1"))
  )
  for (nm in names(sigs10)) {
    sigs10[[nm]]$direction <- "<"
    sigs10[[nm]]$derived_from <- character(0)
  }

  res10 <- run_benchmark(
    signatures = sigs10,
    datasets = ds,
    methods = "mean",
    compute_pr = FALSE,
    compute_null = FALSE,
    min_genes = 2L,
    min_n = 5L,
    min_per_class = 2L,
    coverage_threshold = 0.85
  )
  expect_equal(length(unique(res10$signature)), 10L)
  expect_true(all(res10$coverage_pct >= 85))
  expect_true(all(res10$missing_genes == "" | is.na(res10$missing_genes)))

  # datasets_ici (~2800 genes): reports numeric coverage + missing_genes without crash
  data(datasets_ici, package = "SigBenchR")
  res_ici <- suppressWarnings(run_benchmark(
    signatures = sigs10,
    datasets = datasets_ici["hugo_2016"],
    methods = "mean",
    compute_pr = FALSE,
    compute_null = FALSE,
    min_genes = 1L,
    min_n = 5L,
    min_per_class = 2L
  ))
  expect_equal(nrow(res_ici), 10L)
  expect_true(all(is.finite(res_ici$coverage_pct)))
  expect_true(all(c("found_genes", "missing_genes", "n_found", "n_total") %in% names(res_ici)))
  # At least some signatures should find genes on ICI matrix
  expect_true(any(res_ici$n_found > 0L))
})

test_that("fail-closed upload security: allow_upload default FALSE; shiny parses", {
  expect_identical(formals(run_app)$allow_upload, FALSE)
  app_path <- system.file("shiny", "app.R", package = "SigBenchR")
  if (!nzchar(app_path) || !file.exists(app_path)) {
    app_path <- file.path("inst", "shiny", "app.R")
  }
  expect_true(file.exists(app_path))
  expect_error(parse(file = app_path), NA)
  # Remote entry point pins allow_upload = FALSE
  root_app <- "app.R"
  if (file.exists(root_app)) {
    txt <- paste(readLines(root_app, warn = FALSE), collapse = "\n")
    expect_true(grepl("allow_upload\\s*=\\s*FALSE", txt))
  }
})
