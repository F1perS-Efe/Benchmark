test_that("get_builtin_signatures has literature panel and real negative controls", {
  sigs <- get_builtin_signatures(include_controls = TRUE)
  cats <- vapply(sigs, function(s) as.character(s$category)[1], character(1))
  lit <- sigs[cats == "literature_ici"]
  ctrl <- sigs[cats == "negative_control"]

  expect_gte(length(lit), 8L)
  expect_gte(length(ctrl), 3L)

  # Literature demo signatures preserved by name and gene content
  expect_true(all(c("IFNG6", "Cytokine-5") %in% names(sigs)))
  expect_identical(
    sigs$IFNG6$genes,
    c("IDO1", "CXCL9", "CXCL10", "HLA-DRA", "STAT1", "IFNG")
  )
  expect_identical(
    sigs[["Cytokine-5"]]$genes,
    c("IFNG", "CXCL9", "CXCL13", "IL21", "CCL5")
  )

  # Fake synthetic gene sets excluded from default panel
  expect_false("Random3" %in% names(sigs))
  expect_false("NegCtrl_RandomA" %in% names(sigs))
  expect_false("NegCtrl_RandomB" %in% names(sigs))

  # Real negative controls present
  expect_true(all(
    c(
      "NegCtrl_Housekeeping",
      "NegCtrl_Ribosomal",
      "NegCtrl_Reference_geNorm"
    ) %in% names(sigs)
  ))

  # Neg controls should barely overlap literature (Jaccard < 0.1)
  jac <- jaccard_signatures(sigs)
  for (cn in names(ctrl)) {
    for (ln in names(lit)) {
      row <- jac[
        (jac$signature1 == cn & jac$signature2 == ln) |
          (jac$signature1 == ln & jac$signature2 == cn),
      ]
      expect_equal(nrow(row), 1L)
      expect_lt(row$jaccard[1], 0.1)
    }
  }

  lit_only <- get_builtin_signatures(include_controls = FALSE)
  cats2 <- vapply(lit_only, function(s) as.character(s$category)[1], character(1))
  expect_true(all(cats2 == "literature_ici"))
  expect_false("NegCtrl_Housekeeping" %in% names(lit_only))
})
