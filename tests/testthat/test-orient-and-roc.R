test_that("orient_score yönü yüksek=yanıtlı yapar", {
  s <- c(1, 2, 3)
  expect_equal(orient_score(s, "<"), s)
  expect_equal(orient_score(s, ">"), -s)
})

test_that("roc_auc direction='>' sonrası pROC her zaman '<' kullanır (çift ters yok)", {
  skip_if_not_installed("pROC")

  # Yüksek skor = yanıtsız olacak şekilde kur; direction='>' ile düzeltilmeli
  scores <- c(5, 4, 3, 2, 1, 0)
  response <- c(0, 0, 0, 1, 1, 1)

  # Ham skorlarla direction='<' → AUC < 0.5
  auc_wrong_way <- roc_auc(scores, response, direction = "<")$auc
  expect_lt(auc_wrong_way, 0.5)

  # direction='>' → orient_score ters çevirir → AUC > 0.5
  auc_fixed <- roc_auc(scores, response, direction = ">")$auc
  expect_gt(auc_fixed, 0.5)

  # Çift ters olmamalı: AUC_fixed + AUC_wrong ≈ 1
  expect_equal(auc_fixed + auc_wrong_way, 1, tolerance = 1e-12)
})

test_that("elle Mann-Whitney AUC pROC ile eşleşir", {
  skip_if_not_installed("pROC")

  scores <- c(0.1, 0.2, 0.8, 0.9, 0.3, 0.7)
  response <- c(0, 0, 1, 1, 0, 1)

  # Mann-Whitney / rank-sum AUC (eşitlikte average ranks)
  ranks <- rank(scores, ties.method = "average")
  n_pos <- sum(response == 1)
  n_neg <- sum(response == 0)
  r_pos <- sum(ranks[response == 1])
  auc_mw <- (r_pos - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)

  auc_p <- roc_auc(scores, response, direction = "<")$auc
  expect_equal(auc_p, auc_mw, tolerance = 1e-12)
})

test_that("tek sınıflı yanıt hata verir", {
  skip_if_not_installed("pROC")
  expect_error(
    roc_auc(c(1, 2, 3), c(1, 1, 1), direction = "<"),
    "Both response classes|her iki sinif"
  )
})
