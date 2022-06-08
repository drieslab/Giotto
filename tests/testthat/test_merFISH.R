python_path = NULL 
if(is.null(python_path)) {
  installGiottoEnvironment()
}

expr_path = system.file("extdata", "merFISH_3D_data_expression.txt.gz", package = 'Giotto')
loc_path = system.file("extdata", "merFISH_3D_data_cell_locations.txt", package = 'Giotto')
meta_path = system.file("extdata", "merFISH_3D_metadata.txt", package = 'Giotto')

### TESTS FOR MERFISH MOUSE HYPOTHALMIC PREOPTIC REGION DATASET
# --------------------------------------------------------------

# CREATE GIOTTO OBJECT
object <- createGiottoObject(expression = expr_path,
                             spatial_locs = loc_path)

test_that("Object initialization creates expected Giotto object", {
  
  # S4 object of class "giotto"
  expect_s4_class(object, "giotto")
  
  # gobject contains S4 object of class "dgCMatrix" containing raw expression
  expect_s4_class(object@expression[["cell"]][["rna"]][["raw"]], "dgCMatrix")
  expect_true(all(object@expression[["cell"]][["rna"]][["raw"]]@Dim == c(161, 73655)))
  
  # gobject contains S3 object "data.table" of dimensions 73655 x 4 containing spatial locations
  expect_s3_class(object@spatial_locs[["cell"]][["raw"]], "data.table")
  expect_length(object@spatial_locs[["cell"]][["raw"]][["sdimx"]], 73655)
  expect_length(object@spatial_locs[["cell"]][["raw"]][["sdimy"]], 73655)
  expect_length(object@spatial_locs[["cell"]][["raw"]][["sdimz"]], 73655)
  expect_length(object@spatial_locs[["cell"]][["raw"]][["cell_ID"]], 73655)
  
})

# READ IN METADATA
metadata = data.table::fread(meta_path)
object = addCellMetadata(object, new_metadata = metadata$layer_ID, vector_name = 'layer_ID')
object = addCellMetadata(object, new_metadata = metadata$orig_cell_types, vector_name = 'orig_cell_types')

test_that("Cell metadata are read and added to Giotto object", {
  
  # metadata col names
  expect_named(metadata, c("orig_cell_types", "layer_ID"))
  
  # metadata length matches number of preexisting cell_IDs
  expect_equal(nrow(metadata), length(object@cell_metadata[["cell"]][["rna"]][["cell_ID"]])) 
  
  # metadata length/types after added to object
  expect_length(object@cell_metadata[["cell"]][["rna"]][["layer_ID"]], 73655)
  expect_type(object@cell_metadata[["cell"]][["rna"]][["layer_ID"]], "integer")
  expect_length(object@cell_metadata[["cell"]][["rna"]][["orig_cell_types"]], 73655)
  expect_type(object@cell_metadata[["cell"]][["rna"]][["orig_cell_types"]], "character")
  
})

# FILTER GIOTTO OBJECT
###   kind of arbitrary threshold
###   here we just test dimensions of filtered data, is there a more robust unit test?
filtered_object <- filterGiotto(gobject = object,
                       expression_values = "raw",
                       expression_threshold = 1,
                       feat_det_in_min_cells = 50,
                       min_det_feats_per_cell = 50,
                       verbose = TRUE)

test_that("Data in filtered object is expected size", {
  
  # filtered object expression values have expected dimensions
  expect_true(all(filtered_object@expression[["cell"]][["rna"]][["raw"]]@Dim == c(153, 17814)))
  
  # filtered object spatial locations have expected length
  expect_length(filtered_object@spatial_locs[["cell"]][["raw"]][["sdimx"]], 17814)
  expect_length(filtered_object@spatial_locs[["cell"]][["raw"]][["sdimy"]], 17814)
  expect_length(filtered_object@spatial_locs[["cell"]][["raw"]][["sdimz"]], 17814)
  expect_length(filtered_object@spatial_locs[["cell"]][["raw"]][["cell_ID"]], 17814)
  
  # filtered object metadata has expected length 
  expect_length(filtered_object@cell_metadata[["cell"]][["rna"]][["layer_ID"]], 17814)
  expect_length(filtered_object@cell_metadata[["cell"]][["rna"]][["orig_cell_types"]], 17814)
  
})

# NORMALIZE GIOTTO OBJECT
object <- normalizeGiotto(gobject = object, scalefactor = 10000, verbose = T)

test_that("Normalized data added to giotto object", {
  
  # gobject now also contains S4 object of class "dgCMatrix" containing normalized expression
  expect_s4_class(object@expression[["cell"]][["rna"]][["normalized"]], "dgCMatrix")
  expect_true(all(object@expression[["cell"]][["rna"]][["normalized"]]@Dim == c(161, 73655)))
  
  # gobject now also contains S4 object of class "dgeMatrix" containing scaled expression
  expect_s4_class(object@expression[["cell"]][["rna"]][["scaled"]], "dgeMatrix")
  expect_true(all(object@expression[["cell"]][["rna"]][["scaled"]]@Dim == c(161, 73655)))
  
})

# ADD FEATURE AND CELL STATISTICS TO GIOTTO OBJECT
object <- addStatistics(gobject = object)

test_that("Feature and cell statistics are added to giotto object", {
  
  # gobject cell metadata contains nr_feats, perc_feats, total_expr
  expect_length(object@cell_metadata[["cell"]][["rna"]][["nr_feats"]], 73655)
  expect_type(object@cell_metadata[["cell"]][["rna"]][["nr_feats"]], "integer")
  expect_length(object@cell_metadata[["cell"]][["rna"]][["perc_feats"]], 73655)
  expect_type(object@cell_metadata[["cell"]][["rna"]][["perc_feats"]], "double")
  expect_length(object@cell_metadata[["cell"]][["rna"]][["total_expr"]], 73655)
  expect_type(object@cell_metadata[["cell"]][["rna"]][["total_expr"]], "double")
  
  # gobject feat metadata contains nr_cells, perc_cells, total_expr, mean_expr, mean_expr_det
  expect_length(object@feat_metadata[["cell"]][["rna"]][["nr_cells"]], 161)
  expect_type(object@feat_metadata[["cell"]][["rna"]][["nr_cells"]], "integer")
  expect_length(object@feat_metadata[["cell"]][["rna"]][["perc_cells"]], 161)
  expect_type(object@feat_metadata[["cell"]][["rna"]][["perc_cells"]], "double")
  expect_length(object@feat_metadata[["cell"]][["rna"]][["total_expr"]], 161)
  expect_type(object@feat_metadata[["cell"]][["rna"]][["total_expr"]], "double")
  expect_length(object@feat_metadata[["cell"]][["rna"]][["mean_expr"]], 161)
  expect_type(object@feat_metadata[["cell"]][["rna"]][["mean_expr"]], "double")
  expect_length(object@feat_metadata[["cell"]][["rna"]][["mean_expr_det"]], 161)
  expect_type(object@feat_metadata[["cell"]][["rna"]][["mean_expr_det"]], "double")
  
})

# ADJUST EXPRESSION VALUES 
object <- adjustGiottoMatrix(gobject = object, expression_values = c('normalized'),
                             batch_columns = NULL, covariate_columns = c('layer_ID'),
                             return_gobject = TRUE,
                             update_slot = c('custom'))