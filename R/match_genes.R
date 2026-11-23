#' Imza genlerini ifade matrisinde eslestir
#'
#' Once buyuk harfe cevrilmis tam eslesme; eslesmeyenler icin istege bagli
#' alias kurtarma (dahili `.ici_gene_alias_map()` veya ozel `alias_map`).
#'
#' @param expr Gen x ornek matris (satir adlari gen sembolleri).
#' @param genes Karakter vektoru; imza genleri.
#' @param use_aliases `TRUE` (varsayilan) ise eslesmeyen genler alias ile aranir.
#' @param alias_map Opsiyonel isimli karakter vektoru (`ALIAS = "OFFICIAL"`).
#'   `NULL` ise dahili immuno-onkoloji haritasi kullanilir.
#'
#' @return Liste: `found`, `missing`, `n_found`, `n_total`, `coverage_pct`,
#'   `alias_rescued` (`"OLD->NEW"`), `n_alias_rescued`.
#' @export
match_genes <- function(expr, genes, use_aliases = TRUE, alias_map = NULL) {
  if (is.null(rownames(expr))) {
    stop("`expr` rownames (gene symbols) are required.", call. = FALSE)
  }
  genes <- as.character(genes)
  genes <- genes[!is.na(genes) & nzchar(genes)]
  if (length(genes) == 0L) {
    stop("`genes` cannot be empty.", call. = FALSE)
  }

  genes_up <- toupper(genes)
  rows_up <- toupper(rownames(expr))

  # Yinelenen imza genlerini tekillestir (buyuk harf duzeyinde).
  genes_up <- unique(genes_up)

  hit_idx <- match(genes_up, rows_up)
  missing_up <- genes_up[is.na(hit_idx)]
  found_names <- rownames(expr)[hit_idx[!is.na(hit_idx)]]
  alias_rescued <- character(0)

  if (isTRUE(use_aliases) && length(missing_up) > 0L) {
    amap <- if (is.null(alias_map)) {
      .ici_gene_alias_map()
    } else {
      alias_map
    }
    if (is.null(names(amap)) || any(!nzchar(names(amap)))) {
      stop("`alias_map` must be a named character vector.", call. = FALSE)
    }
    amap_vals <- toupper(as.character(amap))
    names(amap_vals) <- toupper(names(amap))
    amap <- amap_vals

    still_missing <- character(0)
    for (g in missing_up) {
      if (!g %in% names(amap)) {
        still_missing <- c(still_missing, g)
        next
      }
      official <- unname(amap[[g]])
      if (length(official) != 1L || is.na(official) || !nzchar(official)) {
        still_missing <- c(still_missing, g)
        next
      }
      idx <- match(official, rows_up)
      if (is.na(idx)) {
        still_missing <- c(still_missing, g)
      } else {
        found_names <- c(found_names, rownames(expr)[idx])
        alias_rescued <- c(alias_rescued, paste0(g, "->", official))
      }
    }
    missing_up <- still_missing
  }

  # Skorlama icin tekil matris satirlari.
  if (length(found_names) > 0L) {
    found_names <- found_names[!duplicated(toupper(found_names))]
  }

  n_total <- length(genes_up)
  n_matched <- n_total - length(missing_up)
  n_found <- length(found_names)
  coverage_pct <- if (n_total == 0L) {
    0
  } else {
    100 * n_matched / n_total
  }

  list(
    found = found_names,
    missing = missing_up,
    n_found = n_found,
    n_total = n_total,
    coverage_pct = coverage_pct,
    alias_rescued = alias_rescued,
    n_alias_rescued = as.integer(length(alias_rescued))
  )
}
