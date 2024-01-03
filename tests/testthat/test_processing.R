# DATA TO USE
# visium mini expression
g <- GiottoData::loadGiottoMini("visium")




# normalize ####

test_that("expression is normalized", {
  # remove normalized and scaled expr matrices
  g <- setExpression(g, NULL, spat_unit = "cell", feat_type = "rna", name = "normalized", verbose = FALSE)
  g <- setExpression(g, NULL, spat_unit = "cell", feat_type = "rna", name = "scaled", verbose = FALSE)

  expect_false(any(c("normalized", "scaled") %in% list_expression(g)$name))

  g <- normalizeGiotto(g, verbose = FALSE)

  expect_true(any(c("normalized", "scaled") %in% list_expression(g)$name))
})


# hvf ####

test_that("highly variable gene detection", {
  # reset feature metadata
  g <- setFeatureMetadata(g, NULL, spat_unit = "cell", feat_type = "rna", verbose = FALSE)
  expect_false("hvf" %in% names(fDataDT(g)))

  g <- calculateHVF(g)
  expect_true("hvf" %in% names(fDataDT(g)))

  # character "yes" and "no" expected
  checkmate::expect_character(fDataDT(g)$hvf)
})


# statistics ####

test_that("statistics are added", {
  # reset cell metadata
  g <- setCellMetadata(g, NULL, spat_unit = "cell", feat_type = "rna", verbose = FALSE)
  expect_false(any(c("nr_feats", "perc_feats", "total_expr") %in% names(pDataDT(g))))

  g <- addStatistics(g)
  expect_true(any(c("nr_feats", "perc_feats", "total_expr") %in% names(pDataDT(g))))
})

