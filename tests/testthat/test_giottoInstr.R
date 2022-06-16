python_path = NULL 
if(is.null(python_path)) {
  installGiottoEnvironment()
}

getSpatialDataset(dataset = "merfish_preoptic", directory = paste0(getwd(), "/testdata/"))

expr_path = "./testdata/merFISH_3D_data_expression.txt.gz"
loc_path = "./testdata/merFISH_3D_data_cell_locations.txt"
meta_path = "./testdata/merFISH_3D_metadata.txt"


### TESTS FUNCTIONS FOR CREATING/CHANGING GIOTTO INSTRUCTIONS
# --------------------------------------------------------------

# CREATE GIOTTO OBJECT FOR TESTING
object <- createGiottoObject(expression = expr_path,
                             spatial_locs = loc_path,
                             verbose = FALSE)





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
