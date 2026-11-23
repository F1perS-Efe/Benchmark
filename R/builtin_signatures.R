#' Built-in literature ICI signature panel (+ optional negative controls)
#'
#' Returns a curated immuno-oncology gene signature panel with standardized
#' biological and methodological metadata for benchmarking. Demo signatures
#' `Cytokine-5` and `IFNG6` keep the same gene lists as `signatures_demo`.
#' Synthetic fake-gene controls (`Random3`, `NegCtrl_RandomA`,
#' `NegCtrl_RandomB`) are omitted from the default panel so coverage counters
#' reflect real gene absence only. `IMPRES_Like` is omitted because Auslander
#' et al. (2018) IMPRES uses 15 pairwise binary gene comparisons rather than
#' linear gene averaging and was trained/validated on overlapping melanoma ICI
#' cohorts (`van_allen_2015`, `hugo_2016`, `riaz_2017`).
#'
#' Standardized metadata fields per signature:
#' - `signature_name`: canonical signature identifier
#' - `biological_interpretation`: concise description of the biological axis
#' - `genes`: character vector of constituent HGNC gene symbols
#' - `genes_up`, `genes_down`: positive and negative gene subsets (when signed/weighted)
#' - `gene_directions`: named character vector (`"+"` or `"-"`) per gene
#' - `weights`: named numeric vector of published continuous coefficients (e.g.
#'   Ayers/Cristescu 18-gene `T_Cell_Inflamed_GEP` where `CD276 = -0.023900` and
#'   `STAT1 = 0.250229`), or `NULL` for equal-weight signatures
#' - `min_genes`: integer minimum matched gene threshold (`2L` or `3L`)
#' - `species`: `"Homo sapiens"`
#' - `direction`: `"positive"`, `"negative"`, or `"signed"`
#' - `category`: `"literature_ici"` or `"negative_control"`
#' - `derived_from`: character vector of source cohorts (`character(0)` when independent)
#' - `citation`, `original_reference`: bibliographic reference
#' - `doi`: persistent DOI or publication identifier
#' - `scoring_compatibility`: compatible `score_signature()` methods
#' - `known_limitations`: documented biological and methodological caveats
#'
#' @param include_controls If `TRUE` (default), include real negative
#'   controls (housekeeping / ribosomal / geNorm). If `FALSE`, return only
#'   literature ICI signatures.
#'
#' @return Named list of standardized signature definitions.
#' @export
get_builtin_signatures <- function(include_controls = TRUE) {
  gep_weights <- c(
    CCL5 = 0.008346,
    CD27 = 0.072293,
    CD274 = 0.042853,
    CD276 = -0.023900,
    CD8A = 0.031021,
    CMKLR1 = 0.151253,
    CXCL9 = 0.074135,
    CXCR6 = 0.004313,
    `HLA-DQA1` = 0.020091,
    `HLA-DRB1` = 0.058806,
    `HLA-E` = 0.071750,
    IDO1 = 0.060679,
    LAG3 = 0.123895,
    NKG7 = 0.075524,
    PDCD1LG2 = 0.003734,
    PSMB10 = 0.032999,
    STAT1 = 0.250229,
    TIGIT = 0.084767
  )

  lit <- list(
    IFNG6 = .make_builtin_sig(
      signature_name = "IFNG6",
      genes = c("IDO1", "CXCL9", "CXCL10", "HLA-DRA", "STAT1", "IFNG"),
      direction = "positive",
      category = "literature_ici",
      derived_from = character(0),
      min_genes = 3L,
      biological_interpretation = "Core interferon-gamma signaling and antigen presentation response axis",
      original_reference = "Ayers M, Lunceford J, Nebozhyn M, et al. IFN-gamma-related mRNA profile predicts clinical response to PD-1 blockade. J Clin Invest. 2017;127(8):2930-2940.",
      doi = "10.1172/JCI91190",
      known_limitations = "Measures baseline IFN-gamma inflammation; does not capture tumor-intrinsic immune evasion (e.g. JAK1/2 or B2M loss)."
    ),
    Ayers_Expanded18 = .make_builtin_sig(
      signature_name = "Ayers_Expanded18",
      genes = c(
        "CD3D", "IDO1", "CIITA", "CD3E", "CCL5", "GZMK",
        "CD2", "HLA-DRA", "CXCL13", "IL2RG", "NKG7", "HLA-E",
        "CXCR6", "LAG3", "TAGAP", "CXCL10", "STAT1", "GZMB"
      ),
      direction = "positive",
      category = "literature_ici",
      derived_from = character(0),
      min_genes = 3L,
      biological_interpretation = "Expanded 18-gene T/NK-cell infiltration, cytolytic, and MHC class I/II inflammatory profile",
      original_reference = "Ayers M, Lunceford J, Nebozhyn M, et al. IFN-gamma-related mRNA profile predicts clinical response to PD-1 blockade. J Clin Invest. 2017;127(8):2930-2940 (Expanded 18-gene immune signature).",
      doi = "10.1172/JCI91190",
      known_limitations = "Unweighted discovery panel distinct from the locked weighted 18-gene T cell-inflamed GEP."
    ),
    T_Cell_Inflamed_GEP = .make_builtin_sig(
      signature_name = "T_Cell_Inflamed_GEP",
      genes = names(gep_weights),
      weights = gep_weights,
      direction = "signed",
      category = "literature_ici",
      derived_from = character(0),
      min_genes = 3L,
      biological_interpretation = "Weighted 18-gene T cell-inflamed gene expression profile (17 positive IFN-g/adaptive resistance genes and negatively weighted CD276/B7-H3)",
      original_reference = "Ayers M, Lunceford J, Nebozhyn M, et al. J Clin Invest. 2017;127(8):2930-2940; Cristescu R, Mogg R, Ayers M, et al. Science. 2018;362(6411):eaar3593.",
      doi = "10.1172/JCI91190; 10.1126/science.aar3593",
      scoring_compatibility = c("mean", "zscore", "rank", "ss_rank", "median"),
      known_limitations = "Original clinical assay normalizes log2 expression against 11 housekeeping genes on NanoString before applying the 18 linear regression weights (CD276 weight = -0.023900)."
    ),
    Cytolytic_Activity = .make_builtin_sig(
      signature_name = "Cytolytic_Activity",
      genes = c("GZMA", "PRF1"),
      direction = "positive",
      category = "literature_ici",
      derived_from = character(0),
      min_genes = 2L,
      biological_interpretation = "CD8+ CTL and NK effector granule cytolytic activity (CYT)",
      original_reference = "Rooney MS, Shukla SA, Wu CJ, Getz G, Hacohen N. Molecular and genetic properties of tumors associated with local immune cytolytic activity. Cell. 2015;160(1-2):48-61.",
      doi = "10.1016/j.cell.2014.12.033",
      known_limitations = "Original CYT is defined as the geometric mean of GZMA and PRF1 TPM (equivalent to arithmetic mean on log2-transformed expression). Requires both genes (min_genes = 2)."
    ),
    Chemokine_12 = .make_builtin_sig(
      signature_name = "Chemokine_12",
      genes = c(
        "CXCL9", "CXCL10", "CXCL11", "CXCL13", "CCL2", "CCL3", "CCL4",
        "CCL5", "CCL8", "CCL18", "CCL19", "CCL21"
      ),
      direction = "positive",
      category = "literature_ici",
      derived_from = character(0),
      min_genes = 3L,
      biological_interpretation = "12-chemokine tertiary lymphoid structure (TLS) and lymphocyte recruitment axis",
      original_reference = "Messina JL, Fenstermacher DA, Eschrich S, et al. 12-Chemokine gene signature identifies lymph node-like structures in melanoma. Sci Rep. 2012;2:765.",
      doi = "10.1038/srep00765",
      known_limitations = "Originally derived via PCA first principal component across 12 chemokines; includes both T-cell (CXCL9/10/11/13) and myeloid-recruiting (CCL2/18) chemokines."
    ),
    Checkpoint_Axis = .make_builtin_sig(
      signature_name = "Checkpoint_Axis",
      genes = c(
        "PDCD1", "CD274", "PDCD1LG2", "CTLA4", "LAG3",
        "HAVCR2", "TIGIT", "IDO1", "VSIR"
      ),
      direction = "positive",
      category = "literature_ici",
      derived_from = character(0),
      min_genes = 3L,
      biological_interpretation = "Adaptive immune resistance and T-cell activation/exhaustion coinhibitory receptor-ligand axis (excluding non-classical cold-tumor B7 ligands CD276/VTCN1)",
      original_reference = "Pardoll DM. The blockade of immune checkpoints in cancer immunotherapy. Nat Rev Cancer. 2012;12(4):252-264; Ribas A, Wolchok JD. Science. 2018;359(6382):1350-1355.",
      doi = "10.1038/nrc3239; 10.1126/science.aar4060",
      known_limitations = "Reflects IFN-gamma-coupled adaptive checkpoint upregulation in T-cell-inflamed tumors; in non-inflamed tumors, terminal T-cell exhaustion without effector renewal can invert prognostic direction."
    ),
    B7_Immune_Exclusion = .make_builtin_sig(
      signature_name = "B7_Immune_Exclusion",
      genes = c("CD276", "VTCN1", "HHLA2"),
      direction = "negative",
      category = "literature_ici",
      derived_from = character(0),
      min_genes = 2L,
      biological_interpretation = "Non-classical B7 family coinhibitory ligands (B7-H3/CD276, B7-H4/VTCN1, B7-H7/HHLA2) enriched in immune-excluded ('cold') tumor microenvironments",
      original_reference = "Altan M, Pelekanou V, Schalper KA, et al. B7-H3 Expression in NSCLC and Its Association with B7-H4, PD-L1 and Tumor-Infiltrating Lymphocytes. Clin Cancer Res. 2017;23(17):5202-5209; Picarda E, et al. Clin Cancer Res. 2016;22(14):3425-3431.",
      doi = "10.1158/1078-0432.CCR-16-3107; 10.1158/1078-0432.CCR-15-2428",
      known_limitations = "Association with T-cell exclusion and ICI resistance is histology- and context-dependent across solid tumors; separated from PD-1/PD-L1 adaptive checkpoints."
    ),
    Myeloid_Suppression = .make_builtin_sig(
      signature_name = "Myeloid_Suppression",
      genes = c(
        "IL10", "TGFB1", "ARG1", "CD163", "MRC1",
        "CSF1R", "TREM2", "S100A8", "S100A9"
      ),
      direction = "negative",
      category = "literature_ici",
      derived_from = character(0),
      min_genes = 3L,
      biological_interpretation = "Immunosuppressive M2-like tumor-associated macrophage (TAM) and myeloid-derived suppressor cell (MDSC) axis (pan-macrophage marker CD68 and IFN-g-induced IDO1 excluded)",
      original_reference = "Cassetta L, Pollard JW. Targeting macrophages: therapeutic approaches in cancer. Nat Rev Drug Discov. 2018;17(12):887-904; Veglia F, Sanseviero E, Gabrilovich DI. Nat Rev Immunol. 2021;21(8):485-498.",
      doi = "10.1038/nrd.2018.169; 10.1038/s41577-020-00490-y",
      known_limitations = "Curated bulk RNA-seq proxy for myeloid-mediated immunosuppression; excludes pan-macrophage CD68 to avoid conflating M1 macrophage infiltration with M2/MDSC suppression."
    ),
    `Cytokine-5` = .make_builtin_sig(
      signature_name = "Cytokine-5",
      genes = c("IFNG", "CXCL9", "CXCL13", "IL21", "CCL5"),
      direction = "positive",
      category = "literature_ici",
      derived_from = c(
        "van_allen_2015", "hugo_2016", "riaz_2017", "gide_2019"
      ),
      min_genes = 3L,
      biological_interpretation = "5-gene Th1/Tfh cytokine and chemokine ICI response signature",
      original_reference = "Kus et al., EMBO Reports 2025 (TCGEx melanoma ICI cytokine panel).",
      doi = "10.1038/s44319-025-00380-x",
      known_limitations = "Derived across bundled melanoma ICI cohorts (van_allen_2015, hugo_2016, riaz_2017, gide_2019); flagged as non-independent (independent = FALSE) on those cohorts."
    )
  )

  ctrl <- list(
    NegCtrl_Housekeeping = .make_builtin_sig(
      signature_name = "NegCtrl_Housekeeping",
      genes = c("ACTB", "GAPDH", "RPLP0", "PPIA", "HPRT1"),
      direction = "positive",
      category = "negative_control",
      derived_from = character(0),
      min_genes = 3L,
      biological_interpretation = "Constitutive metabolic and cytoskeletal housekeeping negative control (B2M excluded due to MHC-I / IFN-gamma regulation)",
      original_reference = "Eisenberg E, Levanon EY. Human housekeeping genes, revisited. Trends Genet. 2013;29(10):569-574.",
      doi = "10.1016/j.tig.2013.05.010",
      known_limitations = "Negative control expected to yield near-chance discrimination (AUC ~ 0.50); metabolic shifts in highly proliferative tumors can introduce minor variance."
    ),
    NegCtrl_Ribosomal = .make_builtin_sig(
      signature_name = "NegCtrl_Ribosomal",
      genes = c(
        "RPLP0", "RPL7", "RPL21", "RPL23A", "RPL24",
        "RPL35", "RPL38", "RPS2", "RPS9", "RPLP1"
      ),
      direction = "positive",
      category = "negative_control",
      derived_from = character(0),
      min_genes = 3L,
      biological_interpretation = "Structural cytosolic ribosomal subunit negative control panel",
      original_reference = "Thorrez L, Van Deun K, Tranchevent LC, et al. Using ribosomal protein genes as reference: a tale of caution. PLoS One. 2008;3(3):e1854.",
      doi = "10.1371/journal.pone.0001854",
      known_limitations = "Negative control for specificity calibration; not intended as a biological predictor."
    ),
    NegCtrl_Reference_geNorm = .make_builtin_sig(
      signature_name = "NegCtrl_Reference_geNorm",
      genes = c(
        "ACTB", "GAPDH", "HPRT1", "TBP",
        "YWHAZ", "SDHA", "UBC", "POLR2A"
      ),
      direction = "positive",
      category = "negative_control",
      derived_from = character(0),
      min_genes = 3L,
      biological_interpretation = "geNorm reference housekeeping gene panel (B2M excluded)",
      original_reference = "Vandesompele J, De Preter K, Pattyn F, et al. Accurate normalization of real-time quantitative RT-PCR data by geometric averaging of multiple internal control genes. Genome Biol. 2002;3(7):RESEARCH0034.",
      doi = "10.1186/gb-2002-3-7-research0034",
      known_limitations = "Reference normalization control panel; shares ACTB, GAPDH, and HPRT1 with NegCtrl_Housekeeping."
    )
  )

  if (isTRUE(include_controls)) {
    c(lit, ctrl)
  } else {
    lit
  }
}

#' Construct a standardized built-in signature metadata entry
#'
#' @noRd
.make_builtin_sig <- function(signature_name,
                              genes,
                              direction = c("positive", "negative", "signed"),
                              category = c("literature_ici", "negative_control"),
                              derived_from = character(0),
                              min_genes = 3L,
                              weights = NULL,
                              species = "Homo sapiens",
                              biological_interpretation = "",
                              original_reference = "",
                              doi = "",
                              scoring_compatibility = c("mean", "zscore", "median", "rank", "ss_rank"),
                              known_limitations = "") {
  direction <- match.arg(direction)
  category <- match.arg(category)
  genes <- as.character(genes)

  if (!is.null(weights)) {
    weights <- as.numeric(weights)
    names(weights) <- genes
    genes_up <- genes[weights > 0]
    genes_down <- genes[weights < 0]
    gene_dirs <- stats::setNames(ifelse(weights > 0, "+", "-"), genes)
  } else {
    genes_up <- if (identical(direction, "negative")) character(0) else genes
    genes_down <- if (identical(direction, "negative")) genes else character(0)
    dir_sym <- if (identical(direction, "negative")) "-" else "+"
    gene_dirs <- stats::setNames(rep(dir_sym, length(genes)), genes)
  }

  list(
    signature_name = as.character(signature_name)[1],
    biological_interpretation = as.character(biological_interpretation)[1],
    genes = genes,
    genes_up = NULL,
    genes_down = NULL,
    gene_directions = gene_dirs,
    weights = weights,
    min_genes = as.integer(min_genes)[1],
    species = as.character(species)[1],
    direction = direction,
    category = category,
    derived_from = derived_from,
    citation = as.character(original_reference)[1],
    original_reference = as.character(original_reference)[1],
    doi = as.character(doi)[1],
    scoring_compatibility = as.character(scoring_compatibility),
    known_limitations = as.character(known_limitations)[1]
  )
}
