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

  use_custom <- is.null(signature_name) || !nzchar(as.character(signature_name)[1])
  if (use_custom) {
    genes <- as.character(genes)
    genes <- genes[nzchar(genes)]
    if (length(genes) < 1L) {
      stop("Ozel imza icin `genes` gerekli (veya `signature_name` verin).", call. = FALSE)
    }
  } else {
    signature_name <- as.character(signature_name)[1]
  }

  fmt_chr_vec <- function(x) {
    paste0("c(", paste(sprintf('"%s"', x), collapse = ", "), ")")
  }
  fmt_logical <- function(x) if (isTRUE(x)) "TRUE" else "FALSE"

  lines <- c(
    "library(SigBenchR)",
    "data(datasets_demo)",
    "data(signatures_demo)",
    ""
  )

  if (use_custom) {
    lines <- c(
      lines,
      paste0(
        "sigs <- list(Custom = list(genes = ",
        fmt_chr_vec(genes),
        ", direction = \"", direction, "\"))"
      )
    )
  } else {
    lines <- c(
      lines,
      paste0('sigs <- signatures_demo["', signature_name, '"]')
    )
  }

  lines <- c(
    lines,
    paste0("ds <- datasets_demo[", fmt_chr_vec(dataset_names), "]"),
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
