#' Ayni veri setinde iki imzayi eslestirilmis DeLong ile karsilastir
#'
#' Yalnizca ayni veri seti ve ortak gecerli ornekler uzerinde
#' `pROC::roc.test(..., method = "delong", paired = TRUE)` kullanir.
#' Farkli veri setleri arasinda karsilastirmaya izin vermez.
#'
#' Coklu cift karsilastirmalarinda karsilastirma ailesi =
#' ayni cagridaki tum imza ciftleri; `p_adj_bh = p.adjust(p, "BH")`.
#'
#' @param scores1,scores2 Isimli skor vektorleri (ornek id isimleri).
#' @param response 0/1 yanit (isimli veya `scores1` ile hizali).
#' @param dataset1,dataset2 Veri seti kimlikleri; farkliysa hata.
#' @param direction1,direction2 Imza yonleri.
#' @param name1,name2 Imza adlari (rapor icin).
#' @param comparisons Opsiyonel: coklu cift icin liste. Her eleman
#'   `scores1`, `scores2`, `response`, `dataset` (veya dataset1/2),
#'   `direction1`, `direction2`, `name1`, `name2` alanlarini tasiyabilir.
#'   Verilirse toplu BH uygulanir.
#'
#' @return `data.frame`: auc1, auc2, p_value, p_adj_bh, ...
#' @export
compare_signatures <- function(scores1 = NULL,
                               scores2 = NULL,
                               response = NULL,
                               dataset1 = NULL,
                               dataset2 = NULL,
                               direction1 = "<",
                               direction2 = "<",
                               name1 = "sig1",
                               name2 = "sig2",
                               comparisons = NULL) {
  if (!is.null(comparisons)) {
    rows <- lapply(comparisons, function(cmp) {
      ds1 <- if (!is.null(cmp$dataset1)) cmp$dataset1 else cmp$dataset
      ds2 <- if (!is.null(cmp$dataset2)) cmp$dataset2 else cmp$dataset
      compare_signatures_one(
        scores1 = cmp$scores1,
        scores2 = cmp$scores2,
        response = cmp$response,
        dataset1 = ds1,
        dataset2 = ds2,
        direction1 = if (!is.null(cmp$direction1)) cmp$direction1 else "<",
        direction2 = if (!is.null(cmp$direction2)) cmp$direction2 else "<",
        name1 = if (!is.null(cmp$name1)) cmp$name1 else "sig1",
        name2 = if (!is.null(cmp$name2)) cmp$name2 else "sig2"
      )
    })
    out <- do.call(rbind, rows)
    out$p_adj_bh <- stats::p.adjust(out$p_value, method = "BH")
    rownames(out) <- NULL
    return(out)
  }

  one <- compare_signatures_one(
    scores1 = scores1,
    scores2 = scores2,
    response = response,
    dataset1 = dataset1,
    dataset2 = dataset2,
    direction1 = direction1,
    direction2 = direction2,
    name1 = name1,
    name2 = name2
  )
  one$p_adj_bh <- one$p_value
  one
}

#' @keywords internal
compare_signatures_one <- function(scores1,
                                   scores2,
                                   response,
                                   dataset1,
                                   dataset2,
                                   direction1 = "<",
                                   direction2 = "<",
                                   name1 = "sig1",
                                   name2 = "sig2") {
  if (is.null(dataset1) || is.null(dataset2)) {
    stop("`dataset1` ve `dataset2` zorunludur.", call. = FALSE)
  }
  if (!identical(as.character(dataset1), as.character(dataset2))) {
    stop(
      "Farkli veri setleri arasinda DeLong karsilastirmasi yapilamaz (",
      dataset1, " vs ", dataset2, ").",
      call. = FALSE
    )
  }

  # Isimli hizalama
  if (!is.null(names(scores1)) && !is.null(names(scores2))) {
    common <- intersect(names(scores1), names(scores2))
    if (!is.null(names(response))) {
      common <- intersect(common, names(response))
    }
    if (length(common) == 0L) {
      stop("Ortak ornek id bulunamadi.", call. = FALSE)
    }
    scores1 <- as.numeric(scores1[common])
    scores2 <- as.numeric(scores2[common])
    if (!is.null(names(response))) {
      response <- as.integer(response[common])
    } else {
      stop("Isimli skorlar icin `response` da isimli olmalidir.", call. = FALSE)
    }
  } else {
    scores1 <- as.numeric(scores1)
    scores2 <- as.numeric(scores2)
    response <- as.integer(response)
    if (length(scores1) != length(scores2) || length(scores1) != length(response)) {
      stop(
        "scores1, scores2 ve response ayni uzunlukta olmali.",
        call. = FALSE
      )
    }
  }

  ok <- is.finite(scores1) & is.finite(scores2) & !is.na(response)
  scores1 <- scores1[ok]
  scores2 <- scores2[ok]
  response <- response[ok]

  if (length(unique(response)) < 2L) {
    stop("Karsilastirma icin her iki yanit sinifi da gerekli.", call. = FALSE)
  }

  s1 <- orient_score(scores1, direction = direction1)
  s2 <- orient_score(scores2, direction = direction2)

  roc1 <- pROC::roc(
    response = response, predictor = s1,
    levels = c(0, 1), direction = "<", quiet = TRUE
  )
  roc2 <- pROC::roc(
    response = response, predictor = s2,
    levels = c(0, 1), direction = "<", quiet = TRUE
  )

  test <- pROC::roc.test(roc1, roc2, method = "delong", paired = TRUE)

  data.frame(
    dataset = as.character(dataset1),
    signature1 = name1,
    signature2 = name2,
    n = length(response),
    n_pos = sum(response == 1L),
    n_neg = sum(response == 0L),
    auc1 = as.numeric(pROC::auc(roc1)),
    auc2 = as.numeric(pROC::auc(roc2)),
    p_value = as.numeric(test$p.value),
    stringsAsFactors = FALSE
  )
}
