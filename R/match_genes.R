#' Imza genlerini ifade matrisinde eslestir
#'
#' Once buyuk harfe cevrilmis tam eslesme; eslesmeyenler icin istege bagli
#' alias kurtarma (dahili `.ici_gene_alias_map()` veya ozel `alias_map`).
#' Yinelenen sorgular ve ayni kanonik gene cozumlenen birden fazla alias
#' tekillestirilir; boylece `n_found`, `n_total` ve `coverage_pct`
#' (`100 * n_found / n_total`) her zaman matematiksel olarak tutarlidir.
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

  genes_up <- unique(toupper(genes))
  rows_up <- toupper(rownames(expr))

  amap <- NULL
  if (isTRUE(use_aliases)) {
    amap_raw <- if (is.null(alias_map)) {
      .ici_gene_alias_map()
    } else {
      alias_map
    }
    if (is.null(names(amap_raw)) || any(!nzchar(names(amap_raw)))) {
      stop("`alias_map` must be a named character vector.", call. = FALSE)
    }
    amap <- toupper(as.character(amap_raw))
    names(amap) <- toupper(names(amap_raw))
  }

  found_names <- character(0)
  missing_up <- character(0)
  alias_rescued <- character(0)
  seen_targets <- character(0)

  for (g in genes_up) {
    direct_idx <- match(g, rows_up)
    if (!is.na(direct_idx)) {
      target_up <- rows_up[direct_idx]
      if (!target_up %in% seen_targets) {
        seen_targets <- c(seen_targets, target_up)
        found_names <- c(found_names, rownames(expr)[direct_idx])
      }
      next
    }

    if (isTRUE(use_aliases) && !is.null(amap) && g %in% names(amap)) {
      official <- unname(amap[[g]])
      if (length(official) == 1L && !is.na(official) && nzchar(official)) {
        idx <- match(official, rows_up)
        if (!is.na(idx)) {
          target_up <- rows_up[idx]
          alias_rescued <- c(alias_rescued, paste0(g, "->", official))
          if (!target_up %in% seen_targets) {
            seen_targets <- c(seen_targets, target_up)
            found_names <- c(found_names, rownames(expr)[idx])
          }
          next
        }
        # Alias points to a canonical gene that is missing from expr: deduplicate missing canonical target
        if (!official %in% seen_targets && !official %in% missing_up) {
          seen_targets <- c(seen_targets, official)
          missing_up <- c(missing_up, g)
        }
        next
      }
    }

    if (!g %in% seen_targets && !g %in% missing_up) {
      seen_targets <- c(seen_targets, g)
      missing_up <- c(missing_up, g)
    }
  }

  n_found <- length(found_names)
  n_total <- n_found + length(missing_up)
  coverage_pct <- if (n_total == 0L) {
    0
  } else {
    100 * n_found / n_total
  }

  list(
    found = found_names,
    missing = missing_up,
    n_found = as.integer(n_found),
    n_total = as.integer(n_total),
    coverage_pct = coverage_pct,
    alias_rescued = alias_rescued,
    n_alias_rescued = as.integer(length(alias_rescued))
  )
}
