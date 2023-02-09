# python_path = NULL
# if(is.null(python_path)) {
#   installGiottoEnvironment()
# }

if(!require(remotes)){
  install.packages('remotes', repos = 'http://cran.us.r-project.org')
}

if(!require(GiottoData)){
  library(remotes)
  install_github('drieslab/GiottoData')
}

### TESTS FUNCTIONS FOR CREATING/CHANGING GIOTTO INSTRUCTIONS
# --------------------------------------------------------------

instrs = createGiottoInstructions(
  #python_path = NULL,
  show_plot = TRUE,
  return_plot = NULL,
  save_plot = FALSE,
  save_dir = NULL,
  plot_format = "png",
  dpi = 300,
  units = NULL,
  height = NULL,
  width = NULL,
  is_docker = FALSE,
  plot_count = 0,
  fiji_path = NULL
)

test_that("Instructions are created", {
  # createGiottoInstructions
  
  expect_type(instrs, "list")
})


GiottoData::getSpatialDataset(dataset = "merfish_preoptic", directory = paste0(getwd(), "/testdata/"))

expr_path = "./testdata/merFISH_3D_data_expression.txt.gz"
loc_path = "./testdata/merFISH_3D_data_cell_locations.txt"
meta_path = "./testdata/merFISH_3D_metadata.txt"

# CREATE GIOTTO OBJECT FOR TESTING
object <- createGiottoObject(expression = expr_path,
                             spatial_locs = loc_path,
                             instructions = instrs,
                             verbose = FALSE)

# readGiottoInstructions
test_that("readGiottoInstructions reads a few giotto object params correctly", {
  
  expect_type(readGiottoInstructions(object, param = "show_plot"), "logical")
  expect_type(readGiottoInstructions(object, param = "plot_format"), "character")
  expect_type(readGiottoInstructions(object, param = "dpi"), "double")
})

# showGiottoInstructions
test_that("showGiottoInstructions returns expected list", {
  expect_type(showGiottoInstructions(object), "list")
})


# changeGiottoInstructions
object = changeGiottoInstructions(
  object,
  params = c("show_plot", "save_plot"),
  new_values = c(FALSE, TRUE),
  return_gobject = TRUE
)

test_that("changeGiottoInstructions changes instruction params in object", {
  expect_false(readGiottoInstructions(object, param = "show_plot"))
  expect_true(readGiottoInstructions(object, param = "save_plot"))
})


# replaceGiottoInstructions
object = replaceGiottoInstructions(object, instrs)

test_that("replaceGiottoInstructions returns object instructions to original", {
  expect_true(readGiottoInstructions(object, param = "show_plot"))
  expect_false(readGiottoInstructions(object, param = "save_plot"))
})


# ---------------------------------------------
# remove downloaded datasets after tests run
if (file.exists("./testdata/merFISH_3D_data_expression.txt.gz")) {
  unlink("./testdata/merFISH_3D_data_expression.txt.gz")
}

if (file.exists("./testdata/merFISH_3D_data_cell_locations.txt")) {
  unlink("./testdata/merFISH_3D_data_cell_locations.txt")
}

if (file.exists("./testdata/merFISH_3D_metadata.txt")) {
  unlink("./testdata/merFISH_3D_metadata.txt")
}
