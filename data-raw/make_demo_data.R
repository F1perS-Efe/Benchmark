## Sentetik demo veri üretici (Faz 1)
## Çalıştırma (paket kökünden): source("data-raw/make_demo_data.R")

set.seed(42)

cytokine5 <- c("IFNG", "CXCL9", "CXCL13", "IL21", "CCL5")
ifng6 <- c("IDO1", "CXCL9", "CXCL10", "HLA-DRA", "STAT1", "IFNG")
random3 <- c("GENA", "GENB", "GENC")

extra_genes <- c(
  "CD8A", "GZMB", "PRF1", "PDCD1", "CTLA4",
  "GENA", "GENB", "GENC", "GEND", "GENE",
  "GENF", "GENG", "GENH", "GENI", "GENJ"
)
all_genes <- unique(c(cytokine5, ifng6, extra_genes))

make_dataset <- function(n_pos,
                         n_neg,
                         n_sd = 2L,
                         signal = 1.2,
                         scale = c("log2", "linear"),
                         normalization = "TPM",
                         cancer = "SKCM",
                         treatment = "anti-PD1",
                         prefix = "S") {
  scale <- match.arg(scale)
  n <- n_pos + n_neg + n_sd
  sample_ids <- paste0(prefix, sprintf("%03d", seq_len(n)))

  response_label <- c(
    rep("CR", max(0L, n_pos - 1L)),
    if (n_pos > 0L) "PR" else character(0),
    rep("PD", n_neg),
    rep("SD", n_sd)
  )
  response_label <- response_label[seq_len(n)]

  y <- ifelse(response_label %in% c("CR", "PR"), 1,
              ifelse(response_label == "PD", 0, NA_real_))

  expr <- matrix(
    rnorm(length(all_genes) * n, mean = 4, sd = 1),
    nrow = length(all_genes),
    ncol = n,
    dimnames = list(all_genes, sample_ids)
  )

  for (g in intersect(cytokine5, all_genes)) {
    expr[g, ] <- expr[g, ] + ifelse(is.na(y), 0, y * signal) + rnorm(n, 0, 0.3)
  }
  for (g in setdiff(ifng6, cytokine5)) {
    if (g %in% rownames(expr)) {
      expr[g, ] <- expr[g, ] + ifelse(is.na(y), 0, y * (signal * 0.7)) +
        rnorm(n, 0, 0.3)
    }
  }

  if (scale == "linear") {
    expr <- pmax(2^expr - 1, 0)
  }

  sample_meta <- data.frame(
    sample_id = sample_ids,
    response_label = response_label,
    cancer = cancer,
    treatment = treatment,
    timing = "pre",
    stringsAsFactors = FALSE
  )

  list(
    expr = expr,
    sample_meta = sample_meta,
    info = list(
      normalization = normalization,
      scale = scale,
      description = "Sentetik Faz 1 demo verisi"
    )
  )
}

demo_ds_a <- make_dataset(
  n_pos = 20L, n_neg = 25L, n_sd = 3L,
  scale = "log2", normalization = "TPM",
  treatment = "anti-PD1", prefix = "A"
)

demo_ds_b <- make_dataset(
  n_pos = 15L, n_neg = 18L, n_sd = 2L,
  scale = "linear", normalization = "TPM",
  treatment = "anti-CTLA4", prefix = "B",
  signal = 1.0
)

# Cytokine-5 derived_from: Faz 3B gercek kohort adlari (prepare_ici_datasets.R
# ile de guncellenir). IFNG6 bagimsiz dogrulama icin bos birakilir.
ici_cohort_names <- c("van_allen_2015", "hugo_2016", "riaz_2017", "gide_2019")

signatures_demo <- list(
  `Cytokine-5` = list(
    genes = cytokine5,
    direction = "<",
    derived_from = ici_cohort_names,
    citation = "Kus et al., EMBO Reports 2025 (TCGEx)"
  ),
  IFNG6 = list(
    genes = ifng6,
    direction = "<",
    derived_from = character(0),
    citation = "Ayers et al., JCI 2017"
  ),
  Random3 = list(
    genes = random3,
    direction = "<",
    derived_from = character(0),
    citation = "Sentetik negatif kontrol"
  )
)

datasets_demo <- list(
  demo_A = demo_ds_a,
  demo_B = demo_ds_b
)

dir.create("data", showWarnings = FALSE)
save(datasets_demo, file = file.path("data", "datasets_demo.rda"), compress = "xz")
save(signatures_demo, file = file.path("data", "signatures_demo.rda"), compress = "xz")
message("Demo veri yazıldı: data/datasets_demo.rda, data/signatures_demo.rda")
