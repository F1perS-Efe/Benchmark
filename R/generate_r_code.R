#' Shiny / CLI parametrelerinden tekrarlanabilir R kodu uret
#'
#' Secili imza, yontemler ve esikler icin dogrudan calistirilabilir
#' kisa bir R betigi dondurur (genelde 10-20 satir).
#'
#' @param signature_name Hazir demo imza adi (`"Cytokine-5"`, `"IFNG6"`,
#'   `"Random3"`). Ozel gen listesi icin `NULL`.
#' @param genes Ozel imza genleri (karakter vektoru). `signature_name`
#'   `NULL` iken zorunlu.
#' @param methods Skorlama yontemleri.
#' @param direction Imza yonu (`"<"` veya `">"`).
#' @param coverage_threshold Kapsama esigi (0-1).
#' @param n_boot PR bootstrap tekrari.
#' @param n_perm Null permutasyon sayisi.
#' @param seed Rastgelelik tohumu.
#' @param dataset_names Demo veri seti adlari.
#' @param compute_pr PR-AUC hesaplansin mi.
#' @param compute_null Null permutasyon calissin mi.
#'
#' @return Tek bir karakter dizisi (R kodu; satirlar `\\n` ile ayrilir).
#' @export
generate_r_code <- function(signature_name = "Cytokine-5",
                            genes = NULL,
                            methods = c("mean", "zscore"),
                            direction = "<",
                            coverage_threshold = 0.70,
                            n_boot = 200L,
                            n_perm = 50L,
                            seed = 42L,
                            dataset_names = c("demo_A", "demo_B"),
                            compute_pr = TRUE,
                            compute_null = FALSE) {
  methods <- as.character(methods)
  if (length(methods) < 1L) {
    stop("`methods` bos olamaz.", call. = FALSE)
  }
  direction <- match.arg(as.character(direction), c("<", ">"))
  dataset_names <- as.character(dataset_names)
  if (length(dataset_names) < 1L) {
    stop("`dataset_names` bos olamaz.", call. = FALSE)
  }

  demo_names <- if (is.null(signature_name)) {
    character(0)
  } else {
    as.character(signature_name)
  }
  demo_names <- demo_names[nzchar(demo_names)]
  custom_genes <- if (is.null(genes)) {
    character(0)
  } else {
    g <- as.character(genes)
    g[nzchar(g)]
  }
  use_demo <- length(demo_names) > 0L
  use_custom <- length(custom_genes) > 0L
  if (!use_demo && !use_custom) {
    stop("Ozel imza icin `genes` gerekli (veya `signature_name` verin).", call. = FALSE)
  }

  fmt_chr_vec <- function(x) {
    paste0("c(", paste(sprintf('"%s"', x), collapse = ", "), ")")
  }
  fmt_logical <- function(x) if (isTRUE(x)) "TRUE" else "FALSE"

  ici_names <- c("van_allen_2015", "hugo_2016", "riaz_2017", "gide_2019", "liu_2019")
  use_ici_data <- any(dataset_names %in% ici_names)
  data_obj <- if (use_ici_data) "datasets_ici" else "datasets_demo"

  lines <- c(
    "library(SigBenchR)",
    paste0("data(", data_obj, ")"),
    "data(signatures_demo)",
    ""
  )

  if (use_demo && !use_custom) {
    if (length(demo_names) == 1L) {
      lines <- c(
        lines,
        paste0('sigs <- signatures_demo["', demo_names, '"]'),
        paste0('sigs[[1]]$direction <- "', direction, '"')
      )
    } else {
      lines <- c(
        lines,
        paste0("sigs <- signatures_demo[", fmt_chr_vec(demo_names), "]"),
        paste0('for (nm in names(sigs)) sigs[[nm]]$direction <- "', direction, '"')
      )
    }
  } else if (use_custom && !use_demo) {
    lines <- c(
      lines,
      paste0(
        "sigs <- list(Custom = list(genes = ",
        fmt_chr_vec(custom_genes),
        ", direction = \"", direction, "\"))"
      )
    )
  } else {
    lines <- c(
      lines,
      paste0("sigs <- signatures_demo[", fmt_chr_vec(demo_names), "]"),
      paste0('for (nm in names(sigs)) sigs[[nm]]$direction <- "', direction, '"'),
      paste0(
        "sigs$Custom <- list(genes = ",
        fmt_chr_vec(custom_genes),
        ", direction = \"", direction, "\")"
      )
    )
  }

  lines <- c(
    lines,
    paste0("ds <- ", data_obj, "[", fmt_chr_vec(dataset_names), "]"),
    "",
    "res <- run_benchmark(",
    "  signatures = sigs,",
    "  datasets = ds,",
    paste0("  methods = ", fmt_chr_vec(methods), ","),
    paste0("  direction = \"", direction, "\","),
    paste0("  coverage_threshold = ", as.numeric(coverage_threshold), ","),
    paste0("  compute_pr = ", fmt_logical(compute_pr), ","),
    paste0("  n_boot = ", as.integer(n_boot), "L,"),
    paste0("  compute_null = ", fmt_logical(compute_null), ","),
    paste0("  n_perm = ", as.integer(n_perm), "L,"),
    paste0("  seed = ", as.integer(seed), "L"),
    ")",
    "print(res)",
    paste0('plot_auc_heatmap(res, method = "', methods[1], '")')
  )

  paste(lines, collapse = "\n")
}
