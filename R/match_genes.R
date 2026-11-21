#' Imza genlerini ifade matrisinde eslestir
#'
#' Faz 1: buyuk harfe cevrilmis tam eslesme. Takma ad (alias) eslestirmesi yok.
#'
#' @param expr Gen x ornek matris (satir adlari gen sembolleri).
#' @param genes Karakter vektoru; imza genleri.
#'
#' @return Liste: `found`, `missing`, `n_found`, `n_total`, `coverage_pct`.
#' @export
match_genes <- function(expr, genes) {
  if (is.null(rownames(expr))) {
    stop("`expr` satir adlari (gen sembolleri) zorunludur.", call. = FALSE)
  }
  genes <- as.character(genes)
  genes <- genes[!is.na(genes) & nzchar(genes)]
  if (length(genes) == 0L) {
    stop("`genes` bos olamaz.", call. = FALSE)
  }

  genes_up <- toupper(genes)
  rows_up <- toupper(rownames(expr))

  # Yinelenen imza genlerini tekillestir (buyuk harf duzeyinde).
  genes_up <- unique(genes_up)

  hit_idx <- match(genes_up, rows_up)
  found_up <- genes_up[!is.na(hit_idx)]
  missing_up <- genes_up[is.na(hit_idx)]

  # Orijinal satir adlarini koru (eslesenler).
  found_names <- rownames(expr)[hit_idx[!is.na(hit_idx)]]

  n_found <- length(found_names)
  n_total <- length(genes_up)
  coverage_pct <- if (n_total == 0L) 0 else 100 * n_found / n_total

  list(
    found = found_names,
    missing = missing_up,
    n_found = n_found,
    n_total = n_total,
    coverage_pct = coverage_pct
  )
}
