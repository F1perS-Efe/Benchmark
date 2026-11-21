test_that("correlate_reference_axes self-test IFNG6 vs IFN_Gamma cor > 0.99", {
  set.seed(42L)
  genes <- c(
    "IDO1", "CXCL9", "CXCL10", "HLA-DRA", "STAT1", "IFNG",
    "GZMA", "GZMB", "PRF1",
    "CD3D", "CD3E", "CD8A", "CD8B",
    "OTHER1", "OTHER2"
  )
  n <- 40L
  mat <- matrix(
    rnorm(length(genes) * n),
    nrow = length(genes),
    dimnames = list(genes, paste0("S", seq_len(n)))
  )
  ds <- list(
    cohort_A = list(
      expr = mat,
      sample_meta = data.frame(
        sample_id = colnames(mat),
        response_label = rep(c("CR", "PD"), each = 20L),
        stringsAsFactors = FALSE
      ),
      info = list(normalization = "TPM", scale = "log2")
    )
  )
  sigs <- list(
    IFNG6 = list(
      genes = c("IDO1", "CXCL9", "CXCL10", "HLA-DRA", "STAT1", "IFNG")
    )
  )
  out <- correlate_reference_axes(
    datasets = ds,
    signatures = sigs,
    method = "zscore",
    cor_method = "spearman",
    sd_handling = "exclude"
  )
  expect_true(nrow(out) >= 1L)
  self <- out[out$signature == "IFNG6" & out$reference_axis == "IFN_Gamma", , drop = FALSE]
  expect_equal(nrow(self), 1L)
  expect_gt(self$cor[1], 0.99)
  expect_equal(self$n_samples[1], 40L)
})

test_that("correlate_reference_axes sd_handling=exclude drops SD from n_samples", {
  set.seed(7L)
  genes <- c(
    "IDO1", "CXCL9", "CXCL10", "HLA-DRA", "STAT1", "IFNG",
    "GZMA", "GZMB", "PRF1",
    "CD3D", "CD3E", "CD8A", "CD8B"
  )
  n_keep <- 30L
  n_sd <- 8L
  samples <- c(paste0("K", seq_len(n_keep)), paste0("SD", seq_len(n_sd)))
  mat <- matrix(
    rnorm(length(genes) * length(samples)),
    nrow = length(genes),
    dimnames = list(genes, samples)
  )
  labels <- c(rep(c("CR", "PD"), each = 15L), rep("SD", n_sd))
  ds <- list(
    with_sd = list(
      expr = mat,
      sample_meta = data.frame(
        sample_id = samples,
        response_label = labels,
        stringsAsFactors = FALSE
      ),
      info = list(normalization = "FPKM", scale = "log2")
    )
  )
  sigs <- list(
    IFNG6 = list(
      genes = c("IDO1", "CXCL9", "CXCL10", "HLA-DRA", "STAT1", "IFNG")
    )
  )
  out <- correlate_reference_axes(
    datasets = ds,
    signatures = sigs,
    method = "mean",
    sd_handling = "exclude"
  )
  expect_true(nrow(out) >= 1L)
  expect_true(all(out$n_samples <= n_keep))
  expect_false(any(out$n_samples > n_keep))
  enc <- encode_response(labels, sd_handling = "exclude")
  expect_equal(sum(enc$kept), n_keep)
  expect_true(all(out$n_samples == sum(enc$kept)))
})
