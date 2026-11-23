#' Tek veri setinde imza skorlarini uret (Shiny / DeLong yardimcisi)
#'
#' @param dataset Tek veri seti listesi.
#' @param genes Imza genleri (karakter veya isimli isaret vektoru).
#' @param method Skorlama yontemi (`"ss_rank"` dahil).
#' @param direction Yon.
#' @param apply_log2,sd_handling,min_genes Diger ayarlar.
#' @param genes_up,genes_down Isaretli imza kollari (`genes` ile birlikte verilmez).
#' @param signature Optional full signature definition (list or vector). When
#'   provided, `genes`, `genes_up`, `genes_down`, `weights`, and signature-level
#'   `min_genes` are resolved from `signature` if not explicitly supplied.
#' @param weights Optional named numeric vector of continuous gene weights.
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
                              genes_down = NULL,
                              signature = NULL,
                              weights = NULL) {
  if (!is.null(signature)) {
    sig_obj <- .as_signature(signature)
    if (is.null(genes) && is.null(genes_up) && is.null(genes_down)) {
      genes <- sig_obj$genes
      genes_up <- sig_obj$genes_up
      genes_down <- sig_obj$genes_down
    }
    if (is.null(weights)) {
      weights <- sig_obj$weights
    }
    min_genes <- .resolve_min_genes(sig_obj, default = min_genes)
  } else if (is.list(genes)) {
    sig_obj <- .as_signature(genes)
    genes <- sig_obj$genes
    genes_up <- sig_obj$genes_up
    genes_down <- sig_obj$genes_down
    if (is.null(weights)) {
      weights <- sig_obj$weights
    }
    min_genes <- .resolve_min_genes(sig_obj, default = min_genes)
  }

  prep <- .prepare_cohort(
    dataset,
    apply_log2 = apply_log2,
    sd_handling = sd_handling,
    filter_response = TRUE
  )
  expr_full <- prep$expr_kept
  y_kept <- prep$y_kept

  scored <- score_signature(
    expr_full,
    genes = genes,
    genes_up = genes_up,
    genes_down = genes_down,
    weights = weights,
    method = method,
    min_genes = min_genes
  )
  if (isTRUE(scored$skipped)) {
    stop(scored$skip_reason, call. = FALSE)
  }

  sample_ids <- intersect(names(scored$scores), names(y_kept))
  resp <- y_kept[sample_ids]
  sc_kept <- scored$scores[sample_ids]
  keep <- is.finite(sc_kept) & !is.na(resp)

  list(
    scores = sc_kept[keep],
    response = resp[keep],
    coverage = scored$coverage,
    direction = direction,
    expr = expr_full
  )
}
