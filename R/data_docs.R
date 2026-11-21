#' Sentetik demo veri setleri (Faz 1)
#'
#' Iki kucuk sentetik ICI-benzeri kohort. Her biri sade bir listedir:
#' `list(expr, sample_meta, info)`.
#'
#' @format Isimli liste (`demo_A`, `demo_B`).
#' @examples
#' data(datasets_demo)
#' names(datasets_demo)
#' datasets_demo$demo_A$info
"datasets_demo"

#' Demo gen imzalari (Faz 1)
#'
#' Cytokine-5 (Kus et al. 2025), IFNG6 (Ayers et al. 2017) ve Random3.
#'
#' @format Isimli liste; her eleman `list(genes, direction, derived_from, citation)`.
#' @examples
#' data(signatures_demo)
#' signatures_demo$`Cytokine-5`$genes
"signatures_demo"
