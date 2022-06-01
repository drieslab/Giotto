# Tests for functions dependent on Giotto object

expr_path = system.file("extdata", "seqfish_field_expr.txt.gz", package = 'Giotto')
loc_path = system.file("extdata", "seqfish_field_locs.txt", package = 'Giotto')


# Tests for Giotto object creation
# ---------------------------------------------

object <- createGiottoObject(raw_exprs = expr_path,
                             spatial_locs = loc_path)

test_that("Object initialization creates Giotto object", {
  expect_is(object, "giotto")
})
