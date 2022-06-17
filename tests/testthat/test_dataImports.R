python_path = NULL 
if(is.null(python_path)) {
  installGiottoEnvironment()
}

### TESTS FOR DATA IMPORT FUNCTIONS
# ------------------------------------

# getSpatialDataset
getSpatialDataset(dataset = c("Mouse_brain_scRNAseq"), directory = paste0(getwd(), "/testdata/"))

test_that("Spatial dataset was downloaded", {
  expect_true(file.exists("./testdata/brain_sc_expression_matrix.txt.gz"))
  expect_true(file.exists("./testdata/brain_sc_metadata.csv"))
})

# readExprMatrix
expr_mat <- readExprMatrix(paste0(getwd(), "/testdata/brain_sc_expression_matrix.txt.gz"))

test_that("Expression matrix is read correctly", {
  expect_s4_class(expr_mat, "dgCMatrix")
  expect_equal(expr_mat@Dim, c(27998, 8039))
  
  # check a few genes
  expect_equal(expr_mat@Dimnames[[1]][20], "Sgcz")
  expect_equal(expr_mat@Dimnames[[1]][50], 'Zfp804a')
})

# get10Xmatrix_h5
mat <- get10Xmatrix_h5("./testdata/visium_brain_data/filtered_feature_bc_matrix.h5")

test_that("get10Xmatrix_h5 returns list", {
  expect_type(mat, "list")
})

# stitchFieldCoordinates
# TODO

# stitchTileCoordinates
# TODO

# -----------------------------
# remove files after testing
if (file.exists("./testdata/brain_sc_expression_matrix.txt.gz")) {
  unlink("./testdata/brain_sc_expression_matrix.txt.gz")
}

if (file.exists("./testdata/brain_sc_metadata.csv")) {
  unlink("./testdata/brain_sc_metadata.csv")
}

