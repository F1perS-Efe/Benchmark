#' Immuno-oncology HGNC alias map (ASCII, uppercase)
#'
#' Named character vector: `ALIAS = "OFFICIAL_HGNC"`. Used by
#' [match_genes()] when `use_aliases = TRUE` and no custom `alias_map`
#' is supplied. Not exported.
#'
#' @return Named character vector (all uppercase ASCII).
#' @keywords internal
#' @noRd
.ici_gene_alias_map <- function() {
  c(
    PDL1 = "CD274",
    PD_L1 = "CD274",
    "PD-L1" = "CD274",
    B7H1 = "CD274",
    "B7-H1" = "CD274",
    PD1 = "PDCD1",
    "PD-1" = "PDCD1",
    "CTLA-4" = "CTLA4",
    CTLA4 = "CTLA4",
    TIM3 = "HAVCR2",
    "TIM-3" = "HAVCR2",
    "LAG-3" = "LAG3",
    LAG3 = "LAG3",
    TIGIT = "TIGIT",
    HLADRA = "HLA-DRA",
    "HLA-DRA" = "HLA-DRA",
    "IFN-G" = "IFNG",
    IFNG = "IFNG",
    IDO = "IDO1",
    IDO1 = "IDO1",
    OX40 = "TNFRSF4",
    "OX-40" = "TNFRSF4",
    "4-1BB" = "TNFRSF9",
    "41BB" = "TNFRSF9",
    GITR = "TNFRSF18",
    VISTA = "VSIR",
    "PD-L2" = "PDCD1LG2",
    PDL2 = "PDCD1LG2",
    CD8 = "CD8A",
    "CD8-A" = "CD8A",
    GZMB = "GZMB",
    PRF1 = "PRF1",
    CXCL9 = "CXCL9",
    CXCL10 = "CXCL10",
    STAT1 = "STAT1",
    CCR5 = "CCR5",
    "CD3" = "CD3D",
    CD3E = "CD3E",
    NKG7 = "NKG7",
    "TGF-B" = "TGFB1",
    TGFB = "TGFB1",
    "IL-10" = "IL10",
    "B7-H3" = "CD276",
    B7H3 = "CD276",
    "B7-H4" = "VTCN1",
    B7H4 = "VTCN1",
    ICOS = "ICOS",
    "ICOS-L" = "ICOSLG",
    ICOSL = "ICOSLG",
    HVEM = "TNFRSF14",
    BTLA = "BTLA",
    "CD40-L" = "CD40LG",
    CD40L = "CD40LG"
  )
}
