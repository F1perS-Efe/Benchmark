#' Map mouse gene symbols to human HGNC orthologs with transcriptome context preservation
#'
#' Maps mouse gene symbols to human HGNC symbols using either a curated 1:1
#' MGI-HGNC immuno-oncology and reference ortholog table (default), a
#' user-supplied full-genome `ortholog_table` (`data.frame` with columns
#' `mouse_symbol` and `human_symbol`), or an installed ortholog annotation
#' package (`"babelgene"`, `"homologene"`, or `"org.Mm.eg.db"`).
#'
#' Contested polymorphic MHC class I/II paralogs (`H2-K1`/`H2-D1` vs
#' `HLA-A`/`HLA-B`, `H2-Ab1`/`H2-Aa` vs `HLA-DQB1`/`HLA-DQA1`), rodent-absent
#' genes (`CCL18`, `HHLA2`, `GNLY` which has no true mouse ortholog in wild-type
#' mice, and murine `Ly6g`/`Ly6c1` vs human `LY6G6D`), multi-copy mouse
#' chemokine paralogs (`Ccl21a`/`Ccl21b`/`Ccl21c` under `strict_1to1`), and
#' blind `toupper()` conversion are disabled by default
#' (`allow_uppercase_fallback = FALSE`). Ambiguous 1:many mouse-to-human
#' mappings are never resolved silently; they are excluded and reported in
#' `attr(res, "ortholog_stats")` with an explicit warning.
#'
#' When `keep_unmapped_context = TRUE` (default), the full pre-mapping
#' transcriptome background (mapped human rows + unmapped mouse transcriptome
#' rows prefixed with `"__UNMAPPED_BG__:"`) is stored in
#' `attr(res, "transcriptome_bg_expr")` so single-sample rank scoring
#' (`method = "ss_rank"`) and `null_permutation_test()` retain the full
#' transcriptome ranking and sampling universe rather than collapsing to the
#' mapped gene subset.
#'
#' @param expr Numeric gene x sample expression matrix (or data.frame) with
#'   mouse gene symbols as rownames.
#' @param collapse_fun How to collapse duplicate mouse rows (or many:1 mouse
#'   genes when `ortholog_policy = "many_to_1"`) mapping to the same human
#'   symbol: `"mean"` (default), `"max"`, or `"first"`.
#' @param allow_uppercase_fallback If `FALSE` (default), drop symbols absent
#'   from the ortholog table instead of blindly uppercasing.
#' @param ortholog_table Optional custom ortholog mapping `data.frame` with
#'   columns `mouse_symbol` and `human_symbol`, or a character scalar naming an
#'   installed package (`"babelgene"` or `"homologene"`). `NULL` uses the
#'   built-in curated 1:1 MGI-HGNC reference table.
#' @param ortholog_policy Policy for many-to-1 mouse-to-human mappings in the
#'   ortholog dictionary: `"strict_1to1"` (default; drops many:1 paralogs from
#'   the dictionary while still collapsing duplicate rows of the same mouse
#'   symbol in `expr`) or `"many_to_1"` (collapses distinct mouse paralogs onto
#'   the shared human symbol via `collapse_fun`). Ambiguous 1:many mappings
#'   (one mouse gene -> multiple human genes) are always excluded.
#' @param keep_unmapped_context If `TRUE` (default), attach the full
#'   transcriptome matrix (`attr(res, "transcriptome_bg_expr")`) so `ss_rank`
#'   and `null_permutation_test()` preserve the full transcriptome background.
#'
#' @return Numeric matrix with human gene rownames and attributes
#'   `"ortholog_stats"` (list with `n_input`, `n_mapped`, `n_unmapped`,
#'   `n_collapsed`, `n_curated`, `n_ambiguous_dropped`, `ambiguous_mouse_genes`,
#'   `retention_rate`, `loss_rate`, `is_partial_mapping`, `mapping_scope`,
#'   `ortholog_policy`, `source_version`, `annotation_date`, `collapse_fun`,
#'   `context_preserved`) and `"transcriptome_bg_expr"` (when
#'   `keep_unmapped_context = TRUE`).
#' @export
map_mouse_to_human <- function(expr,
                               collapse_fun = c("mean", "max", "first"),
                               allow_uppercase_fallback = FALSE,
                               ortholog_table = NULL,
                               ortholog_policy = c("strict_1to1", "many_to_1"),
                               keep_unmapped_context = TRUE) {
  collapse_fun <- match.arg(collapse_fun)
  ortholog_policy <- match.arg(ortholog_policy)
  if (is.null(rownames(expr))) {
    stop("`expr` rownames (gene symbols) are required.", call. = FALSE)
  }
  mat <- as.matrix(expr)
  if (!is.numeric(mat)) {
    stop("`expr` must be numeric.", call. = FALSE)
  }

  rn <- as.character(rownames(mat))
  valid <- !is.na(rn) & nzchar(trimws(rn))
  mat <- mat[valid, , drop = FALSE]
  rn <- trimws(rn[valid])
  rownames(mat) <- rn
  n_input <- nrow(mat)
  if (n_input == 0L) {
    stop("No valid gene symbols remain in `expr`.", call. = FALSE)
  }

  orth_info <- .resolve_ortholog_dictionary(
    ortholog_table = ortholog_table,
    ortholog_policy = ortholog_policy,
    input_genes = rn
  )
  curated <- orth_info$map
  curated_keys <- toupper(names(curated))
  rn_up <- toupper(rn)

  # Check if any input genes were dropped due to 1:many or many:1 ambiguity
  ambig_in_data <- unique(rn[rn_up %in% orth_info$ambiguous_mouse_upper])
  if (length(ambig_in_data) > 0L) {
    warning(
      sprintf(
        "Excluded %d ambiguous mouse gene(s) under ortholog_policy = '%s': %s",
        length(ambig_in_data),
        ortholog_policy,
        paste(utils::head(ambig_in_data, 10L), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  hit <- match(rn_up, curated_keys)
  n_curated <- sum(!is.na(hit))

  if (isTRUE(allow_uppercase_fallback)) {
    is_ambig <- rn_up %in% orth_info$ambiguous_mouse_upper
    keep_rows <- !is_ambig
    mapped <- ifelse(!is.na(hit[keep_rows]), unname(curated[hit[keep_rows]]), rn_up[keep_rows])
    n_unmapped <- as.integer(sum(!keep_rows))
    unmapped_mat <- mat[!keep_rows, , drop = FALSE]
    mat_kept <- mat[keep_rows, , drop = FALSE]
  } else {
    keep_rows <- !is.na(hit)
    n_unmapped <- as.integer(sum(!keep_rows))
    if (!any(keep_rows)) {
      stop("No genes matched the mouse-to-human ortholog table.", call. = FALSE)
    }
    unmapped_mat <- mat[!keep_rows, , drop = FALSE]
    mat_kept <- mat[keep_rows, , drop = FALSE]
    mapped <- unname(curated[hit[keep_rows]])
  }

  out <- .collapse_rows_by_key(mat_kept, mapped, collapse_fun = collapse_fun)
  n_collapsed <- as.integer(nrow(mat_kept) - nrow(out))

  # Build full transcriptome background (mapped human rows + collapsed unmapped mouse rows)
  bg_expr <- NULL
  if (isTRUE(keep_unmapped_context)) {
    if (nrow(unmapped_mat) > 0L) {
      unmapped_keys <- paste0("__UNMAPPED_BG__:", toupper(rownames(unmapped_mat)))
      unmapped_collapsed <- .collapse_rows_by_key(
        unmapped_mat,
        unmapped_keys,
        collapse_fun = collapse_fun
      )
      bg_expr <- rbind(out, unmapped_collapsed)
    } else {
      bg_expr <- out
    }
  }

  retention_rate <- as.numeric(nrow(out)) / as.numeric(n_input)
  is_partial <- (n_unmapped > 0L) &&
    (retention_rate < 0.50 || (n_input >= 500L && nrow(out) < 500L))

  attr(out, "ortholog_stats") <- list(
    n_input = as.integer(n_input),
    n_mapped = as.integer(nrow(out)),
    n_unmapped = as.integer(n_unmapped),
    n_collapsed = n_collapsed,
    n_curated = as.integer(n_curated),
    n_ambiguous_dropped = as.integer(length(ambig_in_data)),
    ambiguous_mouse_genes = ambig_in_data,
    retention_rate = retention_rate,
    loss_rate = 1 - retention_rate,
    is_partial_mapping = isTRUE(is_partial),
    mapping_scope = orth_info$mapping_scope,
    ortholog_policy = ortholog_policy,
    source_version = orth_info$source_version,
    annotation_date = orth_info$annotation_date,
    collapse_fun = collapse_fun,
    context_preserved = isTRUE(keep_unmapped_context),
    n_background_total = if (!is.null(bg_expr)) as.integer(nrow(bg_expr)) else as.integer(nrow(out))
  )
  if (!is.null(bg_expr)) {
    attr(out, "transcriptome_bg_expr") <- bg_expr
  }
  out
}

#' Collapse matrix rows sharing the same target gene key
#'
#' @noRd
.collapse_rows_by_key <- function(mat, keys, collapse_fun = c("mean", "max", "first")) {
  collapse_fun <- match.arg(collapse_fun)
  if (nrow(mat) == 0L) {
    return(mat)
  }
  if (!anyDuplicated(keys)) {
    out <- mat
    rownames(out) <- keys
    return(out)
  }
  split_idx <- split(seq_len(nrow(mat)), keys)
  out_list <- lapply(split_idx, function(idx) {
    if (length(idx) == 1L) {
      return(as.numeric(mat[idx, ]))
    }
    sub <- mat[idx, , drop = FALSE]
    if (identical(collapse_fun, "first")) {
      return(as.numeric(sub[1L, ]))
    }
    if (identical(collapse_fun, "mean")) {
      v <- colMeans(sub, na.rm = TRUE)
      v[is.nan(v)] <- NA_real_
      return(as.numeric(v))
    }
    apply(sub, 2L, function(x) {
      if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
    })
  })
  out <- do.call(rbind, out_list)
  rownames(out) <- names(split_idx)
  colnames(out) <- colnames(mat)
  out
}

#' Resolve and validate mouse-to-human ortholog dictionary under 1:1 / many:1 policy
#'
#' @noRd
.resolve_ortholog_dictionary <- function(ortholog_table = NULL,
                                         ortholog_policy = c("strict_1to1", "many_to_1"),
                                         input_genes = NULL) {
  ortholog_policy <- match.arg(ortholog_policy)
  query_genes <- if (!is.null(input_genes) && length(input_genes) > 0L) {
    unique(trimws(as.character(input_genes)))
  } else {
    unique(names(.mouse_to_human_immune_map()))
  }
  query_genes <- query_genes[!is.na(query_genes) & nzchar(query_genes)]

  if (is.null(ortholog_table)) {
    raw_map <- .mouse_to_human_immune_map()
    if (identical(ortholog_policy, "many_to_1")) {
      raw_map <- c(
        raw_map,
        "Ccl21a" = "CCL21",
        "Ccl21b" = "CCL21",
        "Ccl21c" = "CCL21"
      )
    }
    df <- data.frame(
      mouse_symbol = names(raw_map),
      human_symbol = unname(raw_map),
      stringsAsFactors = FALSE
    )
    scope <- "curated_mgi_hgnc_reference"
    src_ver <- "MGI_HGNC_curated_1to1_v2026.03"
    ann_date <- "2026-03-01"
  } else if (is.character(ortholog_table) && length(ortholog_table) == 1L) {
    pkg <- trimws(ortholog_table[1])
    if (identical(pkg, "homologene")) {
      if (!requireNamespace("homologene", quietly = TRUE)) {
        stop("Package 'homologene' is not installed; supply a data.frame `ortholog_table` or install 'homologene'.", call. = FALSE)
      }
      hm <- homologene::mouse2human(query_genes)
      df <- data.frame(mouse_symbol = hm$mouseGene, human_symbol = hm$humanGene, stringsAsFactors = FALSE)
      scope <- "homologene_package"
      src_ver <- paste0("homologene_", as.character(utils::packageVersion("homologene")))
      ann_date <- as.character(Sys.Date())
    } else if (identical(pkg, "babelgene")) {
      if (!requireNamespace("babelgene", quietly = TRUE)) {
        stop("Package 'babelgene' is not installed; supply a data.frame `ortholog_table` or install 'babelgene'.", call. = FALSE)
      }
      bg <- babelgene::orthologs(genes = query_genes, species = "mouse", human = FALSE)
      df <- data.frame(mouse_symbol = bg$symbol, human_symbol = bg$human_symbol, stringsAsFactors = FALSE)
      scope <- "babelgene_package"
      src_ver <- paste0("babelgene_", as.character(utils::packageVersion("babelgene")))
      ann_date <- as.character(Sys.Date())
    } else {
      stop("Unsupported `ortholog_table` string: '", pkg, "'. Pass a data.frame with `mouse_symbol` and `human_symbol`.", call. = FALSE)
    }
  } else if (is.data.frame(ortholog_table)) {
    req_cols <- c("mouse_symbol", "human_symbol")
    if (!all(req_cols %in% names(ortholog_table))) {
      stop("`ortholog_table` data.frame must contain columns `mouse_symbol` and `human_symbol`.", call. = FALSE)
    }
    df <- data.frame(
      mouse_symbol = trimws(as.character(ortholog_table$mouse_symbol)),
      human_symbol = toupper(trimws(as.character(ortholog_table$human_symbol))),
      stringsAsFactors = FALSE
    )
    scope <- if (nrow(df) >= 1000L) "full_transcriptome_table" else "custom_ortholog_table"
    src_ver <- if (!is.null(attr(ortholog_table, "source_version"))) {
      as.character(attr(ortholog_table, "source_version"))[1]
    } else {
      "user_supplied_ortholog_table"
    }
    ann_date <- if (!is.null(attr(ortholog_table, "annotation_date"))) {
      as.character(attr(ortholog_table, "annotation_date"))[1]
    } else {
      "unspecified"
    }
  } else {
    stop("`ortholog_table` must be NULL, a data.frame(`mouse_symbol`, `human_symbol`), or a package name.", call. = FALSE)
  }

  ok <- !is.na(df$mouse_symbol) & nzchar(df$mouse_symbol) &
    !is.na(df$human_symbol) & nzchar(df$human_symbol)
  df <- df[ok, , drop = FALSE]
  df$mouse_upper <- toupper(df$mouse_symbol)
  df$human_upper <- toupper(df$human_symbol)
  df <- df[!duplicated(df[, c("mouse_upper", "human_upper"), drop = FALSE]), , drop = FALSE]

  # 1:many check (one mouse gene -> multiple distinct human symbols): always ambiguous
  m_counts <- table(df$mouse_upper)
  ambig_1to_many <- names(m_counts)[m_counts > 1L]

  # many:1 check (multiple distinct mouse genes -> same human symbol)
  df_no_1m <- df[!(df$mouse_upper %in% ambig_1to_many), , drop = FALSE]
  h_counts <- table(df_no_1m$human_upper)
  many_to_1_human <- names(h_counts)[h_counts > 1L]
  ambig_many_to_1 <- if (identical(ortholog_policy, "strict_1to1")) {
    df_no_1m$mouse_upper[df_no_1m$human_upper %in% many_to_1_human]
  } else {
    character(0)
  }

  all_ambig_mouse <- unique(c(ambig_1to_many, ambig_many_to_1))
  df_clean <- df[!(df$mouse_upper %in% all_ambig_mouse), , drop = FALSE]

  map_vec <- stats::setNames(df_clean$human_upper, df_clean$mouse_symbol)
  list(
    map = map_vec,
    ambiguous_mouse_upper = all_ambig_mouse,
    mapping_scope = scope,
    source_version = src_ver,
    annotation_date = ann_date
  )
}

#' Curated 1:1 MGI-HGNC immuno-oncology and reference ortholog table
#'
#' Covers all 1:1 mouse orthologs across built-in literature ICI signatures,
#' reference immune axes, and negative control panels (housekeeping, ribosomal,
#' geNorm), while excluding polymorphic MHC class I/II paralogs (`H2-K1`/`H2-D1`
#' vs `HLA-A/B`, `H2-Ab1`/`H2-Aa` vs `HLA-DQB1/DQA1`) and human genes lacking
#' 1:1 mouse orthologs (`CCL18`, `HHLA2`).
#'
#' @noRd
.mouse_to_human_immune_map <- function() {
  c(
    # Immune checkpoints & costimulatory receptors/ligands
    "Pdcd1" = "PDCD1",
    "Cd274" = "CD274",
    "Pdcd1lg2" = "PDCD1LG2",
    "Ctla4" = "CTLA4",
    "Lag3" = "LAG3",
    "Havcr2" = "HAVCR2",
    "Tigit" = "TIGIT",
    "Vsir" = "VSIR",
    "Cd276" = "CD276",
    "Vtcn1" = "VTCN1",
    "Btla" = "BTLA",
    "Icos" = "ICOS",
    "Icosl" = "ICOSLG",
    "Cd27" = "CD27",
    "Cd70" = "CD70",
    "Cd28" = "CD28",
    "Cd80" = "CD80",
    "Cd86" = "CD86",
    "Tnfrsf4" = "TNFRSF4",
    "Tnfsf4" = "TNFSF4",
    "Tnfrsf9" = "TNFRSF9",
    "Tnfsf9" = "TNFSF9",
    "Tnfrsf18" = "TNFRSF18",
    "Tnfsf18" = "TNFSF18",
    "Cd40" = "CD40",
    "Cd40lg" = "CD40LG",
    # T-cell, NK-cell & cytolytic effector markers
    "Ifng" = "IFNG",
    "Gzma" = "GZMA",
    "Gzmb" = "GZMB",
    "Gzmk" = "GZMK",
    "Prf1" = "PRF1",
    "Nkg7" = "NKG7",
    "Fasl" = "FASLG",
    "Fas" = "FAS",
    "Tnf" = "TNF",
    "Cd8a" = "CD8A",
    "Cd8b1" = "CD8B",
    "Cd4" = "CD4",
    "Cd3d" = "CD3D",
    "Cd3e" = "CD3E",
    "Cd3g" = "CD3G",
    "Cd247" = "CD247",
    "Cd2" = "CD2",
    "Il2rg" = "IL2RG",
    "Il2ra" = "IL2RA",
    "Il2rb" = "IL2RB",
    "Tagap" = "TAGAP",
    "Zap70" = "ZAP70",
    "Lck" = "LCK",
    "Tbx21" = "TBX21",
    "Eomes" = "EOMES",
    "Foxp3" = "FOXP3",
    "Tcf7" = "TCF7",
    "Tox" = "TOX",
    "Pdcd4" = "PDCD4",
    "Klrk1" = "KLRK1",
    "Klrd1" = "KLRD1",
    "Ncr1" = "NCR1",
    # Chemokines, cytokines & receptors
    "Cxcl9" = "CXCL9",
    "Cxcl10" = "CXCL10",
    "Cxcl11" = "CXCL11",
    "Cxcl12" = "CXCL12",
    "Cxcl13" = "CXCL13",
    "Cxcl16" = "CXCL16",
    "Ccl2" = "CCL2",
    "Ccl3" = "CCL3",
    "Ccl4" = "CCL4",
    "Ccl5" = "CCL5",
    "Ccl7" = "CCL7",
    "Ccl8" = "CCL8",
    "Ccl17" = "CCL17",
    "Ccl19" = "CCL19",
    "Ccl20" = "CCL20",
    "Ccl22" = "CCL22",
    "Ccr2" = "CCR2",
    "Ccr4" = "CCR4",
    "Ccr5" = "CCR5",
    "Ccr7" = "CCR7",
    "Cxcr3" = "CXCR3",
    "Cxcr4" = "CXCR4",
    "Cxcr5" = "CXCR5",
    "Cxcr6" = "CXCR6",
    "Cmklr1" = "CMKLR1",
    "Xcl1" = "XCL1",
    "Xcr1" = "XCR1",
    "Il2" = "IL2",
    "Il4" = "IL4",
    "Il6" = "IL6",
    "Il7" = "IL7",
    "Il10" = "IL10",
    "Il12a" = "IL12A",
    "Il12b" = "IL12B",
    "Il15" = "IL15",
    "Il17a" = "IL17A",
    "Il18" = "IL18",
    "Il21" = "IL21",
    # IFN signaling, antigen processing & non-polymorphic presentation factors
    "Stat1" = "STAT1",
    "Stat2" = "STAT2",
    "Stat3" = "STAT3",
    "Stat4" = "STAT4",
    "Stat5a" = "STAT5A",
    "Stat5b" = "STAT5B",
    "Stat6" = "STAT6",
    "Jak1" = "JAK1",
    "Jak2" = "JAK2",
    "Jak3" = "JAK3",
    "Tyk2" = "TYK2",
    "Irf1" = "IRF1",
    "Irf3" = "IRF3",
    "Irf4" = "IRF4",
    "Irf7" = "IRF7",
    "Irf8" = "IRF8",
    "Irf9" = "IRF9",
    "Ido1" = "IDO1",
    "Ido2" = "IDO2",
    "Ciita" = "CIITA",
    "B2m" = "B2M",
    "Tap1" = "TAP1",
    "Tap2" = "TAP2",
    "Tapbp" = "TAPBP",
    "Psmb8" = "PSMB8",
    "Psmb9" = "PSMB9",
    "Psmb10" = "PSMB10",
    "Cd74" = "CD74",
    "Nlrc5" = "NLRC5",
    # Myeloid, TAM, MDSC, DC & stromal suppression markers
    "Tgfb1" = "TGFB1",
    "Tgfb2" = "TGFB2",
    "Tgfb3" = "TGFB3",
    "Tgfbr1" = "TGFBR1",
    "Tgfbr2" = "TGFBR2",
    "Arg1" = "ARG1",
    "Arg2" = "ARG2",
    "Nos2" = "NOS2",
    "Cd68" = "CD68",
    "Cd163" = "CD163",
    "Mrc1" = "MRC1",
    "Csf1" = "CSF1",
    "Csf1r" = "CSF1R",
    "Csf2" = "CSF2",
    "Csf3r" = "CSF3R",
    "Trem2" = "TREM2",
    "Tyrobp" = "TYROBP",
    "S100a8" = "S100A8",
    "S100a9" = "S100A9",
    "Itgam" = "ITGAM",
    "Itgax" = "ITGAX",
    "Itgae" = "ITGAE",
    "Cd14" = "CD14",
    "Batf3" = "BATF3",
    "Clec9a" = "CLEC9A",
    "Vegfa" = "VEGFA",
    "Ptgs2" = "PTGS2",
    "Entpd1" = "ENTPD1",
    "Nt5e" = "NT5E",
    "Adora2a" = "ADORA2A",
    "Fap" = "FAP",
    "Acta2" = "ACTA2",
    "Col1a1" = "COL1A1",
    "Col3a1" = "COL3A1",
    # B-cell & TLS markers
    "Cd19" = "CD19",
    "Ms4a1" = "MS4A1",
    "Cd79a" = "CD79A",
    "Cd79b" = "CD79B",
    "Jchain" = "JCHAIN",
    # Negative controls: housekeeping, geNorm & ribosomal 1:1 orthologs
    "Actb" = "ACTB",
    "Gapdh" = "GAPDH",
    "Rplp0" = "RPLP0",
    "Ppia" = "PPIA",
    "Hprt" = "HPRT1",
    "Tbp" = "TBP",
    "Ywhaz" = "YWHAZ",
    "Sdha" = "SDHA",
    "Ubc" = "UBC",
    "Polr2a" = "POLR2A",
    "Rpl7" = "RPL7",
    "Rpl21" = "RPL21",
    "Rpl23a" = "RPL23A",
    "Rpl24" = "RPL24",
    "Rpl35" = "RPL35",
    "Rpl38" = "RPL38",
    "Rps2" = "RPS2",
    "Rps9" = "RPS9",
    "Rplp1" = "RPLP1"
  )
}
