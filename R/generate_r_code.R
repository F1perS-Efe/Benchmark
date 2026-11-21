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
#' @param dataset_names Veri seti adlari. ICI kohort adlari
#'   (`van_allen_2015`, `hugo_2016`, `riaz_2017`, `gide_2019`) verilirse
#'   uretilen kod `data(datasets_ici)` ve `datasets_ici[...]` kullanir;
#'   aksi halde `datasets_demo`.
#' @param compute_pr PR-AUC hesaplansin mi.
#' @param compute_null Null permutasyon calissin mi.
#' @param sd_handling SD politikasi (`"exclude"`, `"as_0"`, `"as_1"`).
#' @param match_expression Null icin ifade duzeyine eslestirilmis ornekleme.
#'
#' @return Tek bir karakter dizisi (R kodu; satirlar `\\n` ile ayrilir).
#' @export
generate_r_code <- function(signature_name = "Cytokine-5",
                            genes = NULL,
                            methods = c("mean", "zscore"),
                            direction = "<",
                            coverage_threshold = 0.70,
                            n_boot = 200L,
                            n_perm = 1000L,
                            seed = 42L,
                            dataset_names = c("demo_A", "demo_B"),
                            compute_pr = TRUE,
                            compute_null = FALSE,
                            sd_handling = "exclude",
                            match_expression = FALSE) {
  methods <- as.character(methods)
  if (length(methods) < 1L) {
    stop("`methods` bos olamaz.", call. = FALSE)
  }
  direction <- match.arg(as.character(direction), c("<", ">"))
  sd_handling <- match.arg(
    as.character(sd_handling),
    c("exclude", "as_0", "as_1")
  )
  match_expression <- isTRUE(match_expression)
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

  fmt_chr <- function(s) {
    paste(deparse(as.character(s)), collapse = "")
  }
  fmt_chr_vec <- function(x) {
    paste0(
      "c(",
      paste(vapply(x, fmt_chr, character(1)), collapse = ", "),
      ")"
    )
  }
  fmt_logical <- function(x) if (isTRUE(x)) "TRUE" else "FALSE"

  ici_names <- c("van_allen_2015", "hugo_2016", "riaz_2017", "gide_2019", "liu_2019")
  demo_ds_names <- c("demo_A", "demo_B")
  known_names <- c(ici_names, demo_ds_names)
  unknown_ds <- setdiff(dataset_names, known_names)
  all_ici <- length(dataset_names) > 0L && all(dataset_names %in% ici_names)
  all_demo <- length(dataset_names) > 0L && all(dataset_names %in% demo_ds_names)

  if (length(unknown_ds) > 0L) {
    data_obj <- NULL
  } else if (all_ici) {
    data_obj <- "datasets_ici"
  } else if (all_demo) {
    data_obj <- "datasets_demo"
  } else {
    data_obj <- NULL
  }

  lines <- c("library(SigBenchR)")
  if (!is.null(data_obj)) {
    lines <- c(lines, paste0("data(", data_obj, ")"))
  }
  lines <- c(lines, "data(signatures_demo)", "")

  if (use_demo && !use_custom) {
    if (length(demo_names) == 1L) {
      lines <- c(
        lines,
        paste0("sigs <- signatures_demo[", fmt_chr(demo_names), "]"),
        paste0("sigs[[1]]$direction <- ", fmt_chr(direction))
      )
    } else {
      lines <- c(
        lines,
        paste0("sigs <- signatures_demo[", fmt_chr_vec(demo_names), "]"),
        paste0("for (nm in names(sigs)) sigs[[nm]]$direction <- ", fmt_chr(direction))
      )
    }
  } else if (use_custom && !use_demo) {
    lines <- c(
      lines,
      paste0(
        "sigs <- list(Custom = list(genes = ",
        fmt_chr_vec(custom_genes),
        ", direction = ", fmt_chr(direction), "))"
      )
    )
  } else {
    lines <- c(
      lines,
      paste0("sigs <- signatures_demo[", fmt_chr_vec(demo_names), "]"),
      paste0("for (nm in names(sigs)) sigs[[nm]]$direction <- ", fmt_chr(direction)),
      paste0(
        "sigs$Custom <- list(genes = ",
        fmt_chr_vec(custom_genes),
        ", direction = ", fmt_chr(direction), ")"
      )
    )
  }

  if (is.null(data_obj)) {
    lines <- c(
      lines,
      paste0(
        "# WARNING: unrecognized dataset_names = ",
        fmt_chr_vec(dataset_names),
        "; provide a named list `ds` matching validate_dataset()."
      ),
      "ds <- list()  # TODO: user-supplied datasets"
    )
  } else {
    lines <- c(
      lines,
      paste0("ds <- ", data_obj, "[", fmt_chr_vec(dataset_names), "]")
    )
  }

  lines <- c(
    lines,
    "",
    "res <- run_benchmark(",
    "  signatures = sigs,",
    "  datasets = ds,",
    paste0("  methods = ", fmt_chr_vec(methods), ","),
    paste0("  direction = ", fmt_chr(direction), ","),
    paste0("  coverage_threshold = ", as.numeric(coverage_threshold), ","),
    paste0("  sd_handling = ", fmt_chr(sd_handling), ","),
    paste0("  compute_pr = ", fmt_logical(compute_pr), ","),
    paste0("  n_boot = ", as.integer(n_boot), "L,"),
    paste0("  compute_null = ", fmt_logical(compute_null), ","),
    paste0("  n_perm = ", as.integer(n_perm), "L,"),
    paste0("  seed = ", as.integer(seed), "L,"),
    paste0("  match_expression = ", fmt_logical(match_expression)),
    ")",
    "print(res)",
    "print(summarize_auc(res))",
    paste0("plot_auc_heatmap(res, method = ", fmt_chr(methods[1]), ")")
  )

  paste(lines, collapse = "\n")
}
