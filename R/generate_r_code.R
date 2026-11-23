#' Generate a reproducible R script from Shiny or CLI parameters
#'
#' Returns a short, directly runnable R script (typically 10-20 lines) for the
#' selected signatures, scoring methods, and thresholds.
#'
#' @param signature_name Built-in signature name(s) from [get_builtin_signatures()]
#'   or `signatures_demo` (e.g. `"Cytokine-5"`, `"IFNG6"`,
#'   `"T_Cell_Inflamed_GEP"`), or `NULL` for custom genes only.
#' @param genes Custom signature gene vector. Required when `signature_name`
#'   is `NULL`.
#' @param methods Scoring method names.
#' @param direction Default signature direction (`"<"` or `">"`).
#' @param coverage_threshold Coverage threshold in `[0, 1]`.
#' @param n_boot PR bootstrap replicates.
#' @param n_perm Null permutation replicates.
#' @param seed RNG seed.
#' @param dataset_names Dataset names. When all names belong to `datasets_ici`
#'   (`van_allen_2015`, `hugo_2016`, `riaz_2017`, `gide_2019`), generated code
#'   uses `data(datasets_ici)` and `datasets_ici[...]`; when all belong to
#'   `datasets_demo`, it uses `data(datasets_demo)`.
#' @param compute_pr Whether to compute PR-AUC.
#' @param compute_null Whether to run null permutation testing.
#' @param sd_handling Stable disease policy (`"exclude"`, `"as_0"`, `"as_1"`).
#' @param match_expression Whether to use expression-matched null sampling.
#' @param genes_up Optional up-regulated custom signature gene vector (for
#'   signed custom signatures; mutually exclusive with `genes`).
#' @param genes_down Optional down-regulated custom signature gene vector (for
#'   signed custom signatures; mutually exclusive with `genes`).
#' @param weights Optional named numeric vector of continuous gene weights for
#'   a weighted linear custom signature.
#'
#' @return A single character string of R code separated by `\n`.
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
                            match_expression = FALSE,
                            genes_up = NULL,
                            genes_down = NULL,
                            weights = NULL) {
  methods <- as.character(methods)
  if (length(methods) < 1L) {
    stop("`methods` cannot be empty.", call. = FALSE)
  }
  direction <- match.arg(as.character(direction), c("<", ">"))
  sd_handling <- match.arg(
    as.character(sd_handling),
    c("exclude", "as_0", "as_1")
  )
  match_expression <- isTRUE(match_expression)
  dataset_names <- as.character(dataset_names)
  if (length(dataset_names) < 1L) {
    stop("`dataset_names` cannot be empty.", call. = FALSE)
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
  custom_up <- if (is.null(genes_up)) {
    character(0)
  } else {
    g <- as.character(genes_up)
    g[nzchar(g)]
  }
  custom_down <- if (is.null(genes_down)) {
    character(0)
  } else {
    g <- as.character(genes_down)
    g[nzchar(g)]
  }
  use_custom_weights <- !is.null(weights)
  if (use_custom_weights) {
    if (!is.numeric(weights) || is.null(names(weights)) ||
        any(!nzchar(trimws(names(weights)))) || any(!is.finite(weights)) || any(weights == 0)) {
      stop("`weights` must be a named numeric vector of finite non-zero gene coefficients.", call. = FALSE)
    }
  }
  use_demo <- length(demo_names) > 0L
  use_custom_flat <- length(custom_genes) > 0L && !use_custom_weights
  use_custom_signed <- (length(custom_up) > 0L || length(custom_down) > 0L) && !use_custom_weights
  if (length(custom_genes) > 0L && (length(custom_up) > 0L || length(custom_down) > 0L)) {
    stop("Provide either `genes` or `genes_up`/`genes_down`, not both.", call. = FALSE)
  }
  use_custom <- use_custom_flat || use_custom_signed || use_custom_weights
  if (!use_demo && !use_custom) {
    stop(
      "`genes`, `genes_up`/`genes_down`, or `weights` is required for a custom signature (or provide `signature_name`).",
      call. = FALSE
    )
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
  fmt_named_num_vec <- function(w) {
    nms <- names(w)
    vals <- vapply(as.numeric(w), function(v) sprintf("%.8g", v), character(1))
    pairs <- paste0(vapply(nms, fmt_chr, character(1)), " = ", vals)
    paste0("c(", paste(pairs, collapse = ", "), ")")
  }
  fmt_logical <- function(x) if (isTRUE(x)) "TRUE" else "FALSE"
  fmt_custom_list <- function() {
    if (use_custom_weights) {
      w_nms <- trimws(as.character(names(weights)))
      paste0(
        "list(genes = ", fmt_chr_vec(w_nms),
        ", weights = ", fmt_named_num_vec(weights),
        ", direction = ", fmt_chr(direction), ")"
      )
    } else if (use_custom_signed) {
      up_code <- if (length(custom_up) > 0L) fmt_chr_vec(custom_up) else "NULL"
      dn_code <- if (length(custom_down) > 0L) fmt_chr_vec(custom_down) else "NULL"
      paste0(
        "list(genes_up = ", up_code,
        ", genes_down = ", dn_code,
        ", direction = ", fmt_chr(direction), ")"
      )
    } else {
      paste0(
        "list(genes = ",
        fmt_chr_vec(custom_genes),
        ", direction = ", fmt_chr(direction), ")"
      )
    }
  }

  ici_names <- c("van_allen_2015", "hugo_2016", "riaz_2017", "gide_2019")
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

  demo_sig_names <- c("Cytokine-5", "IFNG6", "Random3")
  builtin_pool <- tryCatch(
    get_builtin_signatures(include_controls = TRUE),
    error = function(e) list()
  )
  all_known_sigs <- union(demo_sig_names, names(builtin_pool))
  unknown_sigs <- if (use_demo) setdiff(demo_names, all_known_sigs) else character(0)
  use_demo_only <- use_demo && all(demo_names %in% demo_sig_names)
  use_builtin_only <- use_demo && !use_demo_only && all(demo_names %in% names(builtin_pool))
  needs_demo_sigs <- use_demo && !use_builtin_only && any(demo_names %in% demo_sig_names)

  lines <- c("library(SigBenchR)")
  if (!is.null(data_obj)) {
    lines <- c(lines, paste0("data(", data_obj, ")"))
  }
  if (needs_demo_sigs || (!use_demo && use_custom)) {
    lines <- c(lines, "data(signatures_demo)")
  }
  lines <- c(lines, "")

  if (use_demo) {
    if (length(unknown_sigs) > 0L) {
      lines <- c(
        lines,
        paste0(
          "# WARNING: unrecognized signature_name = ",
          fmt_chr_vec(unknown_sigs),
          "; supply known builtin or custom signatures."
        ),
        "sigs <- list()"
      )
    } else if (use_demo_only) {
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
    } else if (use_builtin_only) {
      lines <- c(
        lines,
        paste0("builtin_sigs <- get_builtin_signatures(include_controls = TRUE)"),
        paste0("sigs <- builtin_sigs[", fmt_chr_vec(demo_names), "]")
      )
    } else {
      lines <- c(
        lines,
        "builtin_sigs <- get_builtin_signatures(include_controls = TRUE)",
        "extra_demo <- signatures_demo[setdiff(names(signatures_demo), names(builtin_sigs))]",
        "all_sigs <- c(builtin_sigs, extra_demo)",
        paste0("sigs <- all_sigs[", fmt_chr_vec(demo_names), "]")
      )
    }
    if (use_custom) {
      lines <- c(
        lines,
        paste0("sigs$Custom <- ", fmt_custom_list())
      )
    }
  } else {
    lines <- c(
      lines,
      paste0("sigs <- list(Custom = ", fmt_custom_list(), ")")
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
