# Tests for merFISH workflow

# install Giotto environment
python_path = NULL 
if(is.null(python_path)) {
  installGiottoEnvironment()
}

expr_path = system.file("extdata", "merFISH_3D_data_expression.txt.gz", package = 'Giotto')
loc_path = system.file("extdata", "merFISH_3D_data_cell_locations.txt", package = 'Giotto')
meta_path = system.file("extdata", "merFISH_3D_metadata.txt", package = 'Giotto')

# Tests for Giotto object creation
# ---------------------------------------------

object <- createGiottoObject(raw_exprs = expr_path,
                             spatial_locs = loc_path)

test_that("Object initialization creates Giotto object", {
  expect_s4_class(object, "giotto")
})
