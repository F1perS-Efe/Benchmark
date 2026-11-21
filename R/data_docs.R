#' Sentetik demo veri setleri (Faz 1)
#'
#' Iki kucuk sentetik ICI-benzeri kohort. Her biri sade bir listedir:
#' `list(expr, sample_meta, info)`. Birim testleri ve hizli ornekler icin;
#' Shiny varsayilani gercek kohortlar (`datasets_ici`) kullanir.
#'
#' @format Isimli liste (`demo_A`, `demo_B`).
#' @examples
#' data(datasets_demo)
#' names(datasets_demo)
#' datasets_demo$demo_A$info
"datasets_demo"

#' Gercek melanoma ICI kohortlari (Faz 3B Paket A)
#'
#' Pre-treatment, hasta basina tek ornek. Genler: demo imza genleri +
#' `seed = 42` ile orneklenmis tarafsiz null havuzu (~2800 gen/kohort).
#' Liu 2019 (cBioPortal) S3 403 nedeniyle dahil edilmedi.
#'
#' @format Isimli liste: `van_allen_2015`, `hugo_2016`, `riaz_2017`, `gide_2019`.
#' @examples
#' data(datasets_ici)
#' names(datasets_ici)
#' datasets_ici$hugo_2016$info[c("n", "n_pos", "n_neg", "normalization")]
"datasets_ici"

#' Demo gen imzalari
#'
#' Cytokine-5 (Kus et al. 2025), IFNG6 (Ayers et al. 2017) ve Random3.
#' `Cytokine-5$derived_from` gercek ICI kohort adlarini icerir
#' (`independent = FALSE`); `IFNG6$derived_from` bostur (`independent = TRUE`).
#'
#' @format Isimli liste; her eleman `list(genes, direction, derived_from, citation)`.
#' @examples
#' data(signatures_demo)
#' signatures_demo$`Cytokine-5`$genes
#' signatures_demo$`Cytokine-5`$derived_from
"signatures_demo"
