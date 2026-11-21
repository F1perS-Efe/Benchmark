#' Imzalar arasi Jaccard gen ortusmesi
#'
#' `J = |A intersect B| / |A union B|`. Genler buyuk harfe cevrilip tekillestirilir.
#'
#' @param signatures Isimli liste; her eleman `list(genes = <chr>, ...)` veya
#'   dogrudan karakter vektoru.
#'
#' @return Uzun `data.frame`: `signature1`, `signature2`, `jaccard`,
#'   `n_intersect`, `n_union`, `n1`, `n2`.
#' @export
jaccard_signatures <- function(signatures) {
  if (is.null(names(signatures)) || any(!nzchar(names(signatures)))) {
    stop("`signatures` isimli bir liste olmalidir.", call. = FALSE)
  }
  if (length(signatures) < 2L) {
    stop("Jaccard icin en az iki imza gerekir.", call. = FALSE)
  }

  gene_sets <- lapply(signatures, function(sig) {
    if (is.character(sig) || is.numeric(sig)) {
      genes <- if (is.numeric(sig)) names(sig) else sig
    } else {
      genes <- sig$genes
      if (is.numeric(genes)) {
        genes <- names(genes)
      }
      if (is.null(genes) || length(genes) == 0L) {
        genes <- c(sig$genes_up, sig$genes_down)
      }
    }
    genes <- as.character(genes)
    genes <- genes[!is.na(genes) & nzchar(genes)]
    if (length(genes) == 0L) {
      stop("Her imzada `genes` veya `genes_up`/`genes_down` zorunludur.", call. = FALSE)
    }
    unique(toupper(genes))
  })

  nms <- names(gene_sets)
  rows <- list()
  k <- 0L
  for (i in seq_len(length(nms) - 1L)) {
    for (j in seq.int(i + 1L, length(nms))) {
      a <- gene_sets[[i]]
      b <- gene_sets[[j]]
      inter <- length(intersect(a, b))
      uni <- length(union(a, b))
      jac <- if (uni == 0L) NA_real_ else inter / uni
      k <- k + 1L
      rows[[k]] <- data.frame(
        signature1 = nms[i],
        signature2 = nms[j],
        jaccard = jac,
        n_intersect = as.integer(inter),
        n_union = as.integer(uni),
        n1 = length(a),
        n2 = length(b),
        stringsAsFactors = FALSE
      )
    }
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
