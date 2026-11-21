## Faz 3B Paket A: gercek melanoma ICI kohortlari
## Calistirma (paket kokunden): source("data-raw/prepare_ici_datasets.R")
##
## Kaynaklar:
## - Van Allen 2015: GitHub TPM (vanallenlab/VanAllen_CTLA4_Science_RNASeq_TPM)
##   Klinik: xmuyulab mel_van_cli_data.csv (clinical benefit -> n_pos=19);
##   sample_meta$treatment = "anti-CTLA4"
## - Hugo 2016: GEO GSE78220 FPKM
## - Riaz 2017: GEO GSE91061 FPKM (PRCR -> PR; UNK dislanir)
## - Gide 2019: GitHub normalized counts (miabioinformatics/Gide_Quek_CancerCell2019)
##   Klinik: xmuyulab mel_gide19_cli_data.csv pre (n=41); ENA PD1_*_PRE monoterapi;
##   sample_meta$treatment = "anti-PD1" (ipiPD1 kombinasyon HARIC)
## - Liu 2019 (cBioPortal): S3 403 / erisim engeli nedeniyle ATLANDI
##
## Tarafsiz null havuzu: tum demo imza genleri + sifir-varyanssiz genlerden
## seed = 42 ile rastgele ~2800 gen; compress = "xz".

cache_dir <- file.path("data-raw", "cache")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
dir.create("data", showWarnings = FALSE)

stopifnot(requireNamespace("readxl", quietly = TRUE))
stopifnot(requireNamespace("withr", quietly = TRUE))

TARGET_N_GENES <- 2800L
GENE_SEED <- 42L

cytokine5 <- c("IFNG", "CXCL9", "CXCL13", "IL21", "CCL5")
ifng6 <- c("IDO1", "CXCL9", "CXCL10", "HLA-DRA", "STAT1", "IFNG")
random3 <- c("GENA", "GENB", "GENC")
signature_genes <- unique(c(cytokine5, ifng6, random3))

# ---------------------------------------------------------------------------
# Yardimcilar
# ---------------------------------------------------------------------------

download_if_missing <- function(url, dest, ...) {
  if (file.exists(dest) && file.info(dest)$size > 0) {
    message("Cache: ", basename(dest))
    return(invisible(dest))
  }
  message("Indiriliyor: ", url)
  utils::download.file(url, destfile = dest, mode = "wb", quiet = FALSE, ...)
  invisible(dest)
}

parse_geo_soft_samples <- function(soft_gz) {
  lines <- readLines(gzfile(soft_gz), warn = FALSE)
  sample_starts <- which(startsWith(lines, "^SAMPLE"))
  rows <- lapply(sample_starts, function(i) {
    block_end <- length(lines)
    nxt <- sample_starts[sample_starts > i]
    if (length(nxt) > 0L) block_end <- nxt[1] - 1L
    block <- lines[i:block_end]
    get_char <- function(key) {
      hit <- grep(paste0("!Sample_characteristics_ch1 = ", key), block, value = TRUE, fixed = TRUE)
      if (length(hit) == 0L) {
        hit <- grep(paste0("!Sample_characteristics_ch1 = ", key), block, value = TRUE)
      }
      if (length(hit) == 0L) return(NA_character_)
      sub("^![^=]+=\\s*", "", hit[1])
    }
    title <- sub("^!Sample_title =\\s*", "", grep("^!Sample_title =", block, value = TRUE)[1])
    gsm <- sub("^!Sample_geo_accession =\\s*", "", grep("^!Sample_geo_accession =", block, value = TRUE)[1])
    chars <- grep("^!Sample_characteristics_ch1 =", block, value = TRUE)
    chars <- sub("^!Sample_characteristics_ch1 =\\s*", "", chars)
    kv <- strsplit(chars, ":\\s*", perl = TRUE)
    keys <- tolower(vapply(kv, `[`, character(1), 1L))
    vals <- vapply(kv, function(x) if (length(x) >= 2L) paste(x[-1], collapse = ":") else NA_character_, character(1))
    names(vals) <- keys
    data.frame(
      title = title,
      gsm = gsm,
      patient_id = unname(vals["patient id"]),
      response_raw = unname(vals["anti-pd-1 response"]),
      response_alt = unname(vals["response"]),
      visit = unname(vals["visit (pre or on treatment)"]),
      biopsy_time = unname(vals["biopsy time"]),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

load_gene_maps <- function(gene_info_gz) {
  message("NCBI gene_info okunuyor...")
  gi <- utils::read.delim(
    gzfile(gene_info_gz),
    header = TRUE,
    sep = "\t",
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  # tax_id == 9606 already in Homo_sapiens file
  entrez <- as.character(gi[[2]])
  symbol <- as.character(gi[[3]])
  keep <- nzchar(symbol) & !is.na(symbol) & symbol != "-"
  entrez_map <- setNames(symbol[keep], entrez[keep])

  dbx <- as.character(gi[[6]])
  ens_map <- character(0)
  has_ens <- grepl("Ensembl:ENSG", dbx, fixed = TRUE)
  if (any(has_ens)) {
    ens_ids <- regmatches(dbx[has_ens], regexpr("ENSG[0-9]+", dbx[has_ens], perl = TRUE))
    ens_map <- setNames(symbol[has_ens], ens_ids)
    ens_map <- ens_map[!duplicated(names(ens_map))]
  }
  list(entrez = entrez_map, ensembl = ens_map)
}

collapse_to_symbols <- function(mat, symbols) {
  stopifnot(length(symbols) == nrow(mat))
  ok <- !is.na(symbols) & nzchar(symbols)
  mat <- mat[ok, , drop = FALSE]
  symbols <- toupper(symbols[ok])
  if (!anyDuplicated(symbols)) {
    rownames(mat) <- symbols
    return(mat)
  }
  split_idx <- split(seq_len(nrow(mat)), symbols)
  keep <- vapply(split_idx, function(idx) {
    if (length(idx) == 1L) return(idx)
    means <- rowMeans(mat[idx, , drop = FALSE], na.rm = TRUE)
    idx[which.max(means)]
  }, integer(1))
  ord <- order(unname(keep))
  out <- mat[unname(keep)[ord], , drop = FALSE]
  rownames(out) <- names(keep)[ord]
  out
}

subsample_neutral_pool <- function(expr,
                                   must_keep = signature_genes,
                                   n_target = TARGET_N_GENES,
                                   seed = GENE_SEED) {
  expr <- as.matrix(expr)
  storage.mode(expr) <- "double"
  rownames(expr) <- toupper(rownames(expr))
  must_keep <- toupper(unique(must_keep))
  keep_sig <- intersect(must_keep, rownames(expr))

  var_ok <- apply(expr, 1L, function(x) {
    x <- x[is.finite(x)]
    length(x) > 1L && stats::sd(x) > 0
  })
  pool <- setdiff(rownames(expr)[var_ok], keep_sig)
  n_extra <- max(0L, as.integer(n_target) - length(keep_sig))
  n_extra <- min(n_extra, length(pool))
  picked <- withr::with_seed(seed, {
    if (n_extra == 0L) character(0) else sample(pool, size = n_extra, replace = FALSE)
  })
  keep <- unique(c(keep_sig, picked))
  list(
    expr = expr[keep, , drop = FALSE],
    n_genes = length(keep),
    n_signature_genes = length(keep_sig),
    n_null_pool = length(picked),
    missing_signature_genes = setdiff(must_keep, keep_sig)
  )
}

finalize_dataset <- function(expr,
                             sample_meta,
                             info,
                             must_keep = signature_genes) {
  # Tek hasta: sample_meta$patient_id varsa
  if ("patient_id" %in% names(sample_meta)) {
    sample_meta$patient_id <- as.character(sample_meta$patient_id)
    dup <- duplicated(sample_meta$patient_id)
    if (any(dup)) {
      message("  Pseudo-replication: ", sum(dup), " tekrarlayan hasta ornegi dusuldu.")
      sample_meta <- sample_meta[!dup, , drop = FALSE]
      expr <- expr[, sample_meta$sample_id, drop = FALSE]
    }
  }

  # Yanit sayimi (SD dislama encode_response ile ayni kural)
  lab <- toupper(as.character(sample_meta$response_label))
  n_excluded_sd <- sum(lab == "SD", na.rm = TRUE)
  n_pos <- sum(lab %in% c("CR", "PR"), na.rm = TRUE)
  n_neg <- sum(lab == "PD", na.rm = TRUE)
  n <- nrow(sample_meta)

  sub <- subsample_neutral_pool(expr, must_keep = must_keep)
  expr2 <- sub$expr
  # sutun sirasi meta ile hizali
  expr2 <- expr2[, sample_meta$sample_id, drop = FALSE]

  info$n <- n
  info$n_pos <- n_pos
  info$n_neg <- n_neg
  info$n_excluded_sd <- as.integer(n_excluded_sd)
  info$n_genes_packaged <- sub$n_genes
  info$n_signature_genes <- sub$n_signature_genes
  info$n_null_pool_genes <- sub$n_null_pool
  info$gene_subsample_seed <- GENE_SEED
  info$missing_signature_genes <- paste(sub$missing_signature_genes, collapse = ",")

  list(
    dataset = list(expr = expr2, sample_meta = sample_meta, info = info),
    summary = data.frame(
      cohort = info$cohort,
      n = n,
      n_pos = n_pos,
      n_neg = n_neg,
      n_excluded_sd = as.integer(n_excluded_sd),
      n_genes = sub$n_genes,
      normalization = info$normalization,
      scale = info$scale,
      source = info$source,
      license = info$license,
      stringsAsFactors = FALSE
    )
  )
}

map_binary_response <- function(x) {
  # xmuyulab/ims_gene_signature klinik CSV kodlamasi:
  # 1 = yanitli / klinik yarar, 0 = ilerleme, -1 = SD
  ifelse(x == 1, "PR", ifelse(x == 0, "PD", ifelse(x == -1, "SD", NA_character_)))
}

# ---------------------------------------------------------------------------
# Indirmeler
# ---------------------------------------------------------------------------

download_if_missing(
  "https://raw.githubusercontent.com/vanallenlab/VanAllen_CTLA4_Science_RNASeq_TPM/master/TPM_RSEM_VAScience2015.txt",
  file.path(cache_dir, "TPM_RSEM_VAScience2015.txt")
)
download_if_missing(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE78nnn/GSE78220/suppl/GSE78220_PatientFPKM.xlsx",
  file.path(cache_dir, "GSE78220_PatientFPKM.xlsx")
)
download_if_missing(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE78nnn/GSE78220/soft/GSE78220_family.soft.gz",
  file.path(cache_dir, "GSE78220_family.soft.gz")
)
download_if_missing(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE91nnn/GSE91061/suppl/GSE91061_BMS038109Sample.hg19KnownGene.fpkm.csv.gz",
  file.path(cache_dir, "GSE91061_fpkm.csv.gz")
)
download_if_missing(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE91nnn/GSE91061/soft/GSE91061_family.soft.gz",
  file.path(cache_dir, "GSE91061_family.soft.gz")
)
download_if_missing(
  "https://raw.githubusercontent.com/miabioinformatics/Gide_Quek_CancerCell2019/master/cancercell_normalized_counts_genenames.txt",
  file.path(cache_dir, "gide_normalized_counts.txt")
)
download_if_missing(
  "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz",
  file.path(cache_dir, "Homo_sapiens.gene_info.gz")
)
download_if_missing(
  "https://raw.githubusercontent.com/xmuyulab/ims_gene_signature/main/data/mel_van_cli_data.csv",
  file.path(cache_dir, "mel_van_cli_data.csv")
)
download_if_missing(
  "https://raw.githubusercontent.com/xmuyulab/ims_gene_signature/main/data/mel_gide19_cli_data.csv",
  file.path(cache_dir, "mel_gide19_cli_data.csv")
)
download_if_missing(
  "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJEB23709&result=read_run&fields=run_accession,sample_alias,sample_title&format=tsv",
  file.path(cache_dir, "gide_ena.tsv")
)

# Liu 2019: cBioPortal S3 403 — atlandi (Paket A, 4 kohort).
message("Liu 2019 (cBioPortal mel_dfci_2019): S3 403 Forbidden — atlandi.")

gene_maps <- load_gene_maps(file.path(cache_dir, "Homo_sapiens.gene_info.gz"))

datasets_ici <- list()
cohort_summaries <- list()

# ---------------------------------------------------------------------------
# 1) Van Allen 2015
# ---------------------------------------------------------------------------
message("=== Van Allen 2015 ===")
va <- utils::read.delim(
  file.path(cache_dir, "TPM_RSEM_VAScience2015.txt"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
ens <- sub("\\..*$", "", va[[1]])
va_mat <- as.matrix(va[, -1, drop = FALSE])
storage.mode(va_mat) <- "double"
rownames(va_mat) <- ens
colnames(va_mat) <- colnames(va)[-1]
syms <- unname(gene_maps$ensembl[ens])
va_mat <- collapse_to_symbols(va_mat, syms)

va_cli <- utils::read.csv(file.path(cache_dir, "mel_van_cli_data.csv"), row.names = 1, check.names = FALSE)
# Klinik etiket seffafligi (van_allen_2015):
# xmuyulab mel_van_cli_data.csv ikili `response` alani kullanilir (1/0; bu kohortta -1/SD yok).
# IMS paketlemesi Van Allen (Science 2015) klinik yarar (clinical benefit) siniflamasina dayanir:
#   yanit (CR/PR) veya uzun sagkalimli SD -> 1 (bizde PR / n_pos),
#   PD veya kisa sagkalimli SD -> 0 (bizde PD / n_neg).
# Sonuc: n = 42, n_pos = 19, n_neg = 23, n_excluded_sd = 0.
# sample_meta$treatment = "anti-CTLA4" (ipilimumab monoterapi).
# Sutun: MEL-IPI_Pat02 -> Pat02
va_patients <- sub("^MEL-IPI_", "", colnames(va_mat))
keep_va <- va_patients %in% rownames(va_cli)
va_mat <- va_mat[, keep_va, drop = FALSE]
va_patients <- va_patients[keep_va]
va_meta <- data.frame(
  sample_id = colnames(va_mat),
  patient_id = va_patients,
  response_label = map_binary_response(va_cli[va_patients, "response"]),
  cancer = "SKCM",
  treatment = "anti-CTLA4",
  timing = "pre",
  stringsAsFactors = FALSE
)
va_meta <- va_meta[!is.na(va_meta$response_label), , drop = FALSE]
va_mat <- va_mat[, va_meta$sample_id, drop = FALSE]

va_out <- finalize_dataset(
  va_mat,
  va_meta,
  info = list(
    cohort = "van_allen_2015",
    normalization = "TPM",
    scale = "linear",
    treatment = "anti-CTLA4",
    endpoint_type = "clinical_benefit",
    cancer_type = "melanoma",
    source = "GitHub vanallenlab/VanAllen_CTLA4_Science_RNASeq_TPM (TPM_RSEM_VAScience2015.txt); klinik: xmuyulab/ims_gene_signature mel_van_cli_data.csv",
    license = "Veri orijinal yayina / dbGaP phs000452 tabi; paket egitim/arastirma amaclidir",
    citation = "Van Allen et al., Science 2015",
    description = "Pre-treatment melanoma anti-CTLA4 (ipilimumab) RNA-seq TPM; n_pos=19 klinik yarar (IMS binary response)"
  )
)
datasets_ici$van_allen_2015 <- va_out$dataset
cohort_summaries$van_allen_2015 <- va_out$summary

# ---------------------------------------------------------------------------
# 2) Hugo 2016
# ---------------------------------------------------------------------------
message("=== Hugo 2016 ===")
hugo_df <- as.data.frame(readxl::read_xlsx(file.path(cache_dir, "GSE78220_PatientFPKM.xlsx")))
hugo_genes <- toupper(as.character(hugo_df[[1]]))
hugo_mat <- as.matrix(hugo_df[, -1, drop = FALSE])
storage.mode(hugo_mat) <- "double"
rownames(hugo_mat) <- hugo_genes
hugo_mat <- collapse_to_symbols(hugo_mat, rownames(hugo_mat))

hugo_soft <- parse_geo_soft_samples(file.path(cache_dir, "GSE78220_family.soft.gz"))
# Ornek secimi: yalniz baseline / pre-treatment; Pt16.OnTx disla; Pt27 icin tek ornek
hugo_cols <- colnames(hugo_mat)
hugo_keep_cols <- hugo_cols[!grepl("OnTx|\\.on", hugo_cols, ignore.case = TRUE)]
# Pt27A/B -> tek hasta: Pt27A tercih
if (all(c("Pt27A.baseline", "Pt27B.baseline") %in% hugo_keep_cols)) {
  hugo_keep_cols <- setdiff(hugo_keep_cols, "Pt27B.baseline")
}
hugo_mat <- hugo_mat[, hugo_keep_cols, drop = FALSE]

hugo_patient <- sub("\\..*$", "", colnames(hugo_mat))
hugo_patient <- sub("A$", "", hugo_patient) # Pt27A -> Pt27
resp_map <- setNames(
  toupper(trimws(hugo_soft$response_raw)),
  toupper(trimws(hugo_soft$patient_id))
)
resp_map <- resp_map[!is.na(names(resp_map)) & !is.na(resp_map)]
hugo_resp_raw <- unname(resp_map[toupper(hugo_patient)])
hugo_resp <- ifelse(
  grepl("COMPLETE", hugo_resp_raw), "CR",
  ifelse(grepl("PARTIAL", hugo_resp_raw), "PR",
         ifelse(grepl("PROGRESSIVE", hugo_resp_raw), "PD",
                ifelse(grepl("STABLE", hugo_resp_raw), "SD", NA_character_)))
)

hugo_meta <- data.frame(
  sample_id = colnames(hugo_mat),
  patient_id = hugo_patient,
  response_label = hugo_resp,
  cancer = "SKCM",
  treatment = "anti-PD1",
  timing = "pre",
  stringsAsFactors = FALSE
)
hugo_meta <- hugo_meta[!is.na(hugo_meta$response_label), , drop = FALSE]
hugo_mat <- hugo_mat[, hugo_meta$sample_id, drop = FALSE]

hugo_out <- finalize_dataset(
  hugo_mat,
  hugo_meta,
  info = list(
    cohort = "hugo_2016",
    normalization = "FPKM",
    scale = "linear",
    treatment = "anti-PD1",
    endpoint_type = "RECIST",
    cancer_type = "melanoma",
    source = "GEO GSE78220 (GSE78220_PatientFPKM.xlsx + family.soft)",
    license = "GEO acik erisim; orijinal yayina atif zorunlu",
    citation = "Hugo et al., Cell 2016",
    description = "Pre-treatment melanoma anti-PD1 RNA-seq FPKM; Pt16 on-tx ve Pt27 pseudo-rep temizlendi"
  )
)
datasets_ici$hugo_2016 <- hugo_out$dataset
cohort_summaries$hugo_2016 <- hugo_out$summary

# ---------------------------------------------------------------------------
# 3) Riaz 2017
# ---------------------------------------------------------------------------
message("=== Riaz 2017 ===")
riaz_df <- utils::read.csv(
  gzfile(file.path(cache_dir, "GSE91061_fpkm.csv.gz")),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
# Ilk sutun Entrez (bazen X olarak gelir)
entrez_ids <- as.character(riaz_df[[1]])
riaz_mat <- as.matrix(riaz_df[, -1, drop = FALSE])
storage.mode(riaz_mat) <- "double"
rownames(riaz_mat) <- entrez_ids
syms_r <- unname(gene_maps$entrez[entrez_ids])
riaz_mat <- collapse_to_symbols(riaz_mat, syms_r)

riaz_soft <- parse_geo_soft_samples(file.path(cache_dir, "GSE91061_family.soft.gz"))
# Soft title == FPKM colnames
riaz_soft$response_use <- ifelse(
  !is.na(riaz_soft$response_alt), riaz_soft$response_alt, riaz_soft$response_raw
)
riaz_soft$visit_use <- ifelse(
  !is.na(riaz_soft$visit), riaz_soft$visit, riaz_soft$biopsy_time
)
pre_idx <- grepl("^pre$", riaz_soft$visit_use, ignore.case = TRUE) |
  grepl("pre-treatment|\\bpre\\b", riaz_soft$visit_use, ignore.case = TRUE)
# Title ile esle
title_to_resp <- setNames(toupper(trimws(riaz_soft$response_use)), riaz_soft$title)
title_to_visit <- setNames(toupper(trimws(riaz_soft$visit_use)), riaz_soft$title)

riaz_cols <- colnames(riaz_mat)
is_pre <- grepl("_Pre_", riaz_cols, ignore.case = FALSE) |
  grepl("^PRE$", unname(title_to_visit[riaz_cols]), ignore.case = TRUE)
riaz_cols_pre <- riaz_cols[is_pre]
riaz_mat <- riaz_mat[, riaz_cols_pre, drop = FALSE]

riaz_resp_raw <- unname(title_to_resp[colnames(riaz_mat)])
# PRCR -> PR; UNK disla
riaz_resp <- ifelse(
  riaz_resp_raw %in% c("PRCR", "CR", "PR"),
  ifelse(riaz_resp_raw == "PRCR", "PR", riaz_resp_raw),
  ifelse(riaz_resp_raw %in% c("PD", "SD"), riaz_resp_raw, NA_character_)
)
# UNK ve digerleri NA
riaz_resp[riaz_resp_raw %in% c("UNK", "NA", "")] <- NA_character_

riaz_patient <- sub("_.*$", "", colnames(riaz_mat))
riaz_meta <- data.frame(
  sample_id = colnames(riaz_mat),
  patient_id = riaz_patient,
  response_label = riaz_resp,
  response_raw = riaz_resp_raw,
  cancer = "SKCM",
  treatment = "anti-PD1",
  timing = "pre",
  stringsAsFactors = FALSE
)
# UNK / NA disla
riaz_meta <- riaz_meta[!is.na(riaz_meta$response_label), , drop = FALSE]
riaz_mat <- riaz_mat[, riaz_meta$sample_id, drop = FALSE]
riaz_meta$response_raw <- NULL

riaz_out <- finalize_dataset(
  riaz_mat,
  riaz_meta,
  info = list(
    cohort = "riaz_2017",
    normalization = "FPKM",
    scale = "linear",
    treatment = "anti-PD1",
    endpoint_type = "RECIST",
    cancer_type = "melanoma",
    source = "GEO GSE91061 (BMS038 FPKM + family.soft); PRCR->PR, UNK excluded",
    license = "GEO acik erisim; orijinal yayina atif zorunlu",
    citation = "Riaz et al., Cell 2017",
    description = "Pre-treatment melanoma anti-PD1 (nivolumab) FPKM; bir hasta = bir Pre ornek"
  )
)
datasets_ici$riaz_2017 <- riaz_out$dataset
cohort_summaries$riaz_2017 <- riaz_out$summary

# ---------------------------------------------------------------------------
# 4) Gide 2019
# ---------------------------------------------------------------------------
message("=== Gide 2019 ===")
gide_df <- utils::read.delim(
  file.path(cache_dir, "gide_normalized_counts.txt"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
gide_genes <- toupper(as.character(gide_df[["Gene"]]))
# ID, Gene, Description sonrasi ornekler
gide_sample_cols <- setdiff(colnames(gide_df), c("ID", "Gene", "Description"))
gide_pre_cols <- gide_sample_cols[grepl("_PRE$", gide_sample_cols)]
gide_mat <- as.matrix(gide_df[, gide_pre_cols, drop = FALSE])
storage.mode(gide_mat) <- "double"
rownames(gide_mat) <- gide_genes
gide_mat <- collapse_to_symbols(gide_mat, rownames(gide_mat))

ena <- utils::read.delim(file.path(cache_dir, "gide_ena.tsv"), stringsAsFactors = FALSE)
alias_to_err <- setNames(as.character(ena$run_accession), as.character(ena$sample_alias))
gide_cli <- utils::read.csv(file.path(cache_dir, "mel_gide19_cli_data.csv"), row.names = 1, check.names = FALSE)
# Klinik etiket / tedavi kolu seffafligi (gide_2019):
# xmuyulab mel_gide19_cli_data.csv icinde `treatment` sutunu zamanlama (pre/on) anlamindadir.
# IMS pre alt kumesi n = 41; ENA ornek adlari bu 41 satirin tamaminda `PD1_*_PRE` (anti-PD1
# monoterapi). `ipiPD1_*_PRE` (anti-PD1 + anti-CTLA4 kombinasyon) dahil edilmez.
# Bu nedenle sample_meta$treatment = "anti-PD1" (kombinasyon degil).
# Yanit: IMS binary response (1->PR, 0->PD, -1->SD) -> n_pos=19, n_neg=16, n_SD=6.
gide_cli_pre <- gide_cli[gide_cli$treatment == "pre", , drop = FALSE]

gide_err <- unname(alias_to_err[colnames(gide_mat)])
keep_g <- !is.na(gide_err) & gide_err %in% rownames(gide_cli_pre)
gide_mat <- gide_mat[, keep_g, drop = FALSE]
gide_err <- gide_err[keep_g]
gide_resp <- map_binary_response(gide_cli_pre[gide_err, "response"])

gide_patient <- sub("_PRE$", "", colnames(gide_mat))
# Guvenlik: yalniz PD1_* monoterapi ornekleri (ipiPD1 kombinasyon disi)
stopifnot(all(grepl("^PD1_", colnames(gide_mat))))
gide_meta <- data.frame(
  sample_id = colnames(gide_mat),
  patient_id = gide_patient,
  response_label = gide_resp,
  cancer = "SKCM",
  treatment = "anti-PD1",
  timing = "pre",
  stringsAsFactors = FALSE
)
gide_meta <- gide_meta[!is.na(gide_meta$response_label), , drop = FALSE]
gide_mat <- gide_mat[, gide_meta$sample_id, drop = FALSE]

gide_out <- finalize_dataset(
  gide_mat,
  gide_meta,
  info = list(
    cohort = "gide_2019",
    normalization = "normalized_counts",
    scale = "linear",
    treatment = "anti-PD1",
    endpoint_type = "RECIST",
    cancer_type = "melanoma",
    source = "GitHub miabioinformatics/Gide_Quek_CancerCell2019 (cancercell_normalized_counts_genenames.txt); klinik: xmuyulab mel_gide19_cli_data.csv + ENA PRJEB23709",
    license = "ENA PRJEB23709 / orijinal yayina tabi; paket egitim/arastirma amaclidir",
    citation = "Gide et al., Cancer Cell 2019",
    description = "Pre-treatment melanoma anti-PD1 monoterapi (PD1_*_PRE; ipiPD1 kombinasyon haric); n=41"
  )
)
datasets_ici$gide_2019 <- gide_out$dataset
cohort_summaries$gide_2019 <- gide_out$summary

# ---------------------------------------------------------------------------
# Ozet + kaydet
# ---------------------------------------------------------------------------
summary_tbl <- do.call(rbind, cohort_summaries)
rownames(summary_tbl) <- NULL
print(summary_tbl)

utils::write.csv(
  summary_tbl,
  file = file.path("data-raw", "ici_cohort_summary.csv"),
  row.names = FALSE
)

# Dogrula
for (nm in names(datasets_ici)) {
  ds <- datasets_ici[[nm]]
  stopifnot(identical(colnames(ds$expr), ds$sample_meta$sample_id))
  stopifnot(ds$info$scale %in% c("log2", "linear"))
  message(sprintf(
    "%s: n=%d n_pos=%d n_neg=%d n_sd=%d genes=%d",
    nm, ds$info$n, ds$info$n_pos, ds$info$n_neg, ds$info$n_excluded_sd,
    nrow(ds$expr)
  ))
}

save(datasets_ici, file = file.path("data", "datasets_ici.rda"), compress = "xz")

# signatures_demo: Cytokine-5 derived_from gercek kohort adlari
load(file.path("data", "signatures_demo.rda"))
signatures_demo[["Cytokine-5"]]$derived_from <- names(datasets_ici)
signatures_demo[["IFNG6"]]$derived_from <- character(0)
signatures_demo[["Random3"]]$derived_from <- character(0)
save(signatures_demo, file = file.path("data", "signatures_demo.rda"), compress = "xz")

data_bytes <- sum(file.info(list.files("data", full.names = TRUE))$size, na.rm = TRUE)
message(sprintf(
  "data/ toplam boyutu: %.2f MB (hedef < 5 MB)",
  data_bytes / (1024^2)
))
if (data_bytes >= 5 * 1024^2) {
  warning("data/ 5 MB ustunde; TARGET_N_GENES dusurun.", call. = FALSE)
}

message("Yazildi: data/datasets_ici.rda, data/signatures_demo.rda, data-raw/ici_cohort_summary.csv")
