#' Tek veri setinde imza skorlarini uret (Shiny / DeLong yardimcisi)
#'
#' @param dataset Tek veri seti listesi.
#' @param genes Imza genleri.
#' @param method Skorlama yontemi.
#' @param direction Yon.
#' @param apply_log2,sd_handling,min_genes Diger ayarlar.
#'
#' @return Liste: `scores` (isimli), `response` (isimli, 0/1, SD dislanmis NA),
#'   `coverage`.
#' @export
score_for_dataset <- function(dataset,
                              genes,
                              method = "mean",
                              direction = "<",
                              apply_log2 = TRUE,
                              sd_handling = "exclude",
                              min_genes = 3L) {
  ds <- validate_dataset(dataset)
  expr <- collapse_duplicate_genes(as.matrix(ds$expr))
  expr <- prepare_expression(expr, scale = ds$info$scale, apply_log2 = apply_log2)
  enc <- encode_response(ds$sample_meta$response_label, sd_handling = sd_handling)

  scored <- score_signature(expr, genes = genes, method = method, min_genes = min_genes)
  if (isTRUE(scored$skipped)) {
    stop(scored$skip_reason, call. = FALSE)
  }

  sample_ids <- names(scored$scores)
  meta_idx <- match(sample_ids, ds$sample_meta$sample_id)
  resp <- enc$response[meta_idx]
  names(resp) <- sample_ids
  keep <- enc$kept[meta_idx] & is.finite(scored$scores)

  list(
    scores = scored$scores[keep],
    response = resp[keep],
    coverage = scored$coverage,
    direction = direction,
    expr = expr
  )
}
