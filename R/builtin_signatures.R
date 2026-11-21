#' Built-in literature ICI signature panel (+ optional negative controls)
#'
#' Returns an expanded immuno-oncology signature panel for benchmarking.
#' Demo signatures Cytokine-5 and IFNG6 keep the same gene lists as
#' `signatures_demo` (data object unchanged). Synthetic fake-gene controls
#' (`Random3`, `NegCtrl_RandomA`, `NegCtrl_RandomB`) are omitted from the
#' default panel so coverage counters reflect real gene absence only.
#'
#' Metadata per signature:
#' - `direction`: `"positive"`, `"negative"`, or `"signed"` (ROC mapping via
#'   [run_benchmark()]; legacy `"<"` / `">"` also accepted)
#' - `category`: `"literature_ici"` or `"negative_control"`
#' - `derived_from`: character vector of source cohorts, or `NULL` (fail-closed
#'   `independent = NA`) / `character(0)` (independent)
#'
#' @param include_controls If `TRUE` (default), include real negative
#'   controls (housekeeping / ribosomal / geNorm). If `FALSE`, return only
#'   literature ICI signatures.
#'
#' @return Named list of signature definitions.
#' @export
get_builtin_signatures <- function(include_controls = TRUE) {
  lit <- list(
    IFNG6 = list(
      genes = c("IDO1", "CXCL9", "CXCL10", "HLA-DRA", "STAT1", "IFNG"),
      direction = "positive",
      category = "literature_ici",
      derived_from = character(0),
      citation = "Ayers et al., JCI 2017"
    ),
    Ayers_IFNG6 = list(
      genes = c("IDO1", "CXCL9", "CXCL10", "HLA-DRA", "STAT1", "IFNG"),
      direction = "positive",
      category = "literature_ici",
      derived_from = character(0),
      citation = "Ayers et al., JCI 2017 (alias of IFNG6)"
    ),
    Ayers_Expanded18 = list(
      genes = c(
        "IDO1", "CXCL9", "CXCL10", "CXCL11", "HLA-DRA", "HLA-E",
        "STAT1", "IFNG", "CCR5", "CXCL13", "CCL5", "GZMB",
        "CD3D", "CD3E", "CD2", "IL2RG", "NKG7", "TAGAP"
      ),
      direction = "positive",
      category = "literature_ici",
      derived_from = character(0),
      citation = "Ayers et al., JCI 2017 expanded GEP core"
    ),
    T_Cell_Inflamed_GEP = list(
      genes = c(
        "CXCL9", "CXCL10", "CXCL11", "STAT1", "IDO1", "HLA-E",
        "CMKLR1", "NKG7", "CD8A", "CCL5", "CXCR6", "TIGIT"
      ),
      direction = "positive",
      category = "literature_ici",
      derived_from = character(0),
      citation = "T cell-inflamed GEP (Ayers / Cristescu literature panel)"
    ),
    Cytolytic_Activity = list(
      genes = c("GZMA", "PRF1"),
      direction = "positive",
      category = "literature_ici",
      derived_from = character(0),
      citation = "Rooney et al. cytolytic activity (GZMA, PRF1)"
    ),
    Chemokine_12 = list(
      genes = c(
        "CXCL9", "CXCL10", "CXCL11", "CCL2", "CCL3", "CCL4",
        "CCL5", "CCL8", "CCL18", "CCL19", "CCL21", "CX3CL1"
      ),
      direction = "positive",
      category = "literature_ici",
      derived_from = character(0),
      citation = "Chemokine / T-cell trafficking panel"
    ),
    Checkpoint_Axis = list(
      genes = c(
        "PDCD1", "CD274", "CTLA4", "LAG3", "HAVCR2", "TIGIT",
        "IDO1", "VSIR", "CD276", "VTCN1"
      ),
      direction = "positive",
      category = "literature_ici",
      derived_from = character(0),
      citation = "Immune checkpoint axis gene panel"
    ),
    IMPRES_Like = list(
      genes = c(
        "PDCD1", "CD274", "CTLA4", "CD28", "CD80", "CD86",
        "CD27", "CD70", "TNFRSF9", "TNFSF9", "TNFRSF4", "TNFSF4",
        "TNFRSF18", "BTLA", "CD40", "CD40LG"
      ),
      direction = "positive",
      category = "literature_ici",
      derived_from = character(0),
      citation = "IMPRES-like pairwise checkpoint gene set (Auslander et al.)"
    ),
    Myeloid_Suppression = list(
      genes = c("IL10", "TGFB1", "ARG1", "CD68", "CD163", "MRC1", "IDO1"),
      direction = "negative",
      category = "literature_ici",
      derived_from = character(0),
      citation = "Myeloid suppression / M2-like markers (higher = worse response)"
    ),
    `Cytokine-5` = list(
      genes = c("IFNG", "CXCL9", "CXCL13", "IL21", "CCL5"),
      direction = "positive",
      category = "literature_ici",
      derived_from = c(
        "van_allen_2015", "hugo_2016", "riaz_2017", "gide_2019"
      ),
      citation = "Kus et al., EMBO Reports 2025 (TCGEx)"
    )
  )

  ctrl <- list(
    NegCtrl_Housekeeping = list(
      genes = c("ACTB", "GAPDH", "RPLP0", "B2M", "PPIA"),
      direction = "positive",
      category = "negative_control",
      derived_from = character(0),
      citation = "Housekeeping negative control"
    ),
    NegCtrl_Ribosomal = list(
      genes = c(
        "RPLP0", "RPL7", "RPL21", "RPL23A", "RPL24",
        "RPL35", "RPL38", "RPS2", "RPS9", "RPLP1"
      ),
      direction = "positive",
      category = "negative_control",
      derived_from = character(0),
      citation = "Ribosomal protein negative control"
    ),
    NegCtrl_Reference_geNorm = list(
      genes = c(
        "ACTB", "GAPDH", "B2M", "HPRT1", "TBP",
        "YWHAZ", "SDHA", "UBC"
      ),
      direction = "positive",
      category = "negative_control",
      derived_from = character(0),
      citation = "geNorm-style reference gene negative control"
    )
  )

  if (isTRUE(include_controls)) {
    c(lit, ctrl)
  } else {
    lit
  }
}
