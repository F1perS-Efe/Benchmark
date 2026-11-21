#' Tek veri setinde imza skorlarini uret (Shiny / DeLong yardimcisi)
#'
#' @param dataset Tek veri seti listesi.
#' @param genes Imza genleri (karakter veya isimli isaret vektoru).
#' @param method Skorlama yontemi (`"ss_rank"` dahil).
#' @param direction Yon.
#' @param apply_log2,sd_handling,min_genes Diger ayarlar.
#' @param genes_up,genes_down Isaretli imza kollari (`genes` ile birlikte verilmez).
#'
#' @return Liste: `scores` (isimli), `response` (isimli, 0/1, SD dislanmis NA),
#'   `coverage`.
#' @export
score_for_dataset <- function(dataset,
                              genes = NULL,
                              method = "mean",
                              direction = "<",
                              apply_log2 = TRUE,
                              sd_handling = "exclude",
                              min_genes = 3L,
                              genes_up = NULL,
                              genes_down = NULL) {
  ds <- validate_dataset(dataset)
  expr <- collapse_duplicate_genes(as.matrix(ds$expr))
  expr <- prepare_expression(expr, scale = ds$info$scale, apply_log2 = apply_log2)
  enc <- encode_response(ds$sample_meta$response_label, sd_handling = sd_handling)

  # P0-3: skorlamadan once SD / NA yanit orneklerini cikar
  expr_kept <- expr[, enc$kept, drop = FALSE]
  y_kept <- enc$response[enc$kept]
  names(y_kept) <- colnames(expr_kept)

  scored <- score_signature(
    expr_kept,
    genes = genes,
    genes_up = genes_up,
    genes_down = genes_down,
    method = method,
    min_genes = min_genes
  )
  if (isTRUE(scored$skipped)) {
    stop(scored$skip_reason, call. = FALSE)
  }

  sample_ids <- names(scored$scores)
  resp <- y_kept[sample_ids]
  keep <- is.finite(scored$scores) & !is.na(resp)

  list(
    scores = scored$scores[keep],
    response = resp[keep],
    coverage = scored$coverage,
    direction = direction,
    expr = expr_kept
  )
}
