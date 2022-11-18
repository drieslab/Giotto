python_path = NULL
if(is.null(python_path)) {
  installGiottoEnvironment()
}

getSpatialDataset(dataset = "merfish_preoptic", directory = paste0(getwd(), "/testdata/"))

expr_path = "./testdata/merFISH_3D_data_expression.txt.gz"
loc_path = "./testdata/merFISH_3D_data_cell_locations.txt"
meta_path = "./testdata/merFISH_3D_metadata.txt"

### TESTS FOR MERFISH MOUSE HYPOTHALMIC PREOPTIC REGION DATASET
# --------------------------------------------------------------

# CREATE GIOTTO OBJECT
object <- createGiottoObject(expression = expr_path,
                             spatial_locs = loc_path,
                             verbose = FALSE)

test_that("Object initialization creates expected Giotto object", {

  # S4 object of class "giotto"
  expect_s4_class(object, "giotto")

  # gobject contains S4 object of class "dgCMatrix" containing raw expression
  expect_s4_class(object@expression[["cell"]][["rna"]][["raw"]], "dgCMatrix")
  expect_true(all(object@expression[["cell"]][["rna"]][["raw"]]@Dim == c(161, 73655)))

  # gobject contains S4 object "spatLocsObj" of dimensions 73655 x 4 containing spatial locations
  st = get_spatial_locations(test, spat_unit = 'cell', spat_loc_name = 'raw', return_spatlocs_Obj = TRUE)
  expect_identical(st@coordinates, object@spatial_locs[["cell"]][["raw"]]@coordinates)
  expect_s4_class(object@spatial_locs[["cell"]][["raw"]], "spatLocsObj")
  expect_length(object@spatial_locs[["cell"]][["raw"]]@coordinates[["sdimx"]], 73655)
  expect_length(object@spatial_locs[["cell"]][["raw"]]@coordinates[["sdimy"]], 73655)
  expect_length(object@spatial_locs[["cell"]][["raw"]]@coordinates[["sdimz"]], 73655)
  expect_length(object@spatial_locs[["cell"]][["raw"]]@coordinates[["cell_ID"]], 73655)

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
filtered_object <- filterGiotto(gobject = object,
                       expression_values = "raw",
                       expression_threshold = 1,
                       feat_det_in_min_cells = 50,
                       min_det_feats_per_cell = 50,
                       verbose = FALSE)

test_that("Data in filtered object is expected size", {

  # filtered object expression values have expected dimensions
  expect_true(all(filtered_object@expression[["cell"]][["rna"]][["raw"]]@Dim == c(153, 17814)))

  # filtered object spatial locations have expected length
  expect_length(filtered_object@spatial_locs[["cell"]][["raw"]]@coordinates[["sdimx"]], 17814)
  expect_length(filtered_object@spatial_locs[["cell"]][["raw"]]@coordinates[["sdimy"]], 17814)
  expect_length(filtered_object@spatial_locs[["cell"]][["raw"]]@coordinates[["sdimz"]], 17814)
  expect_length(filtered_object@spatial_locs[["cell"]][["raw"]]@coordinates[["cell_ID"]], 17814)

  # filtered object metadata has expected length
  expect_length(filtered_object@cell_metadata[["cell"]][["rna"]][["layer_ID"]], 17814)
  expect_length(filtered_object@cell_metadata[["cell"]][["rna"]][["orig_cell_types"]], 17814)

})

# NORMALIZE GIOTTO OBJECT
object <- normalizeGiotto(gobject = object, scalefactor = 10000, verbose = F)

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

# ADJUST EXPRESSION VALUES FOR BATCH/COVARIATES
object <- adjustGiottoMatrix(gobject = object, expression_values = c('normalized'),
                             batch_columns = NULL, covariate_columns = c('layer_ID'),
                             return_gobject = TRUE,
                             update_slot = c('custom'))

test_that("Adjusted values are created in 'custom' slot", {

  # expression now also contains custom object of class double
  expect_type(object@expression[["cell"]][["rna"]][["custom"]], "double")
  expect_equal(nrow(object@expression[["cell"]][["rna"]][["custom"]]), 161)
  expect_equal(ncol(object@expression[["cell"]][["rna"]][["custom"]]), 73655)

})

# RUN DIMENSION REDUCTION
object <- runPCA(gobject = object,
                 genes_to_use = NULL,
                 scale_unit = FALSE,
                 center = TRUE,
                 verbose = FALSE)

test_that("PCA S4 object is created as expected", {

  # s4 object of class "dimObj"
  expect_s4_class(object@dimension_reduction[["cells"]][["cell"]][["rna"]][["pca"]][["pca"]], "dimObj")

  # coordinates double of dims 73655 x 100
  expect_equal(nrow(object@dimension_reduction[["cells"]][["cell"]][["rna"]][["pca"]][["pca"]]@coordinates), 73655)
  expect_equal(ncol(object@dimension_reduction[["cells"]][["cell"]][["rna"]][["pca"]][["pca"]]@coordinates), 100)

  # test a few arbitrary coordinates
  expect_equal(object@dimension_reduction[["cells"]][["cell"]][["rna"]][["pca"]][["pca"]]@coordinates[5],
               -5.612915,
               tolerance = 1*10^-3)
  expect_equal(object@dimension_reduction[["cells"]][["cell"]][["rna"]][["pca"]][["pca"]]@coordinates[10],
               -24.66273,
               tolerance = 1*10^-3)

})

# UMAP
object <- runUMAP(object, dimensions_to_use = 1:8, n_components = 3, n_threads = 4, verbose = FALSE)

test_that("UMAP S4 object is created as expected", {

  # s4 object of class "dimObj"
  expect_s4_class(object@dimension_reduction[["cells"]][["cell"]][["rna"]][["umap"]][["umap"]], "dimObj")

  # coordinates double of dims 73655 x 3
  expect_equal(nrow(object@dimension_reduction[["cells"]][["cell"]][["rna"]][["umap"]][["umap"]]@coordinates), 73655)
  expect_equal(ncol(object@dimension_reduction[["cells"]][["cell"]][["rna"]][["umap"]][["umap"]]@coordinates), 3)

  # test a few arbitrary coordinates
  show_failure(expect_equal(!!object@dimension_reduction[["cells"]][["cell"]][["rna"]][["umap"]][["umap"]]@coordinates[20],
                            -3.2,
                            tolerance = 1*10^-3))
  show_failure(expect_equal(!!object@dimension_reduction[["cells"]][["cell"]][["rna"]][["umap"]][["umap"]]@coordinates[40],
                            10,
                            tolerance = 1*10^-3))

})

# CREATE NETWORK
object <- createNearestNetwork(gobject = object, dimensions_to_use = 1:8, k = 15, verbose = FALSE)

test_that("sNN S3 object is created as expected", {

  # igraph s3 object
  expect_s3_class(object@nn_network[["cell"]][["rna"]][["sNN"]][["sNN.pca"]], "igraph")

})

# LEIDEN CLUSTERING
object <- doLeidenCluster(gobject = object, resolution = 0.2, n_iterations = 200,
                          name = 'leiden_0.2')

test_that("New clusters are added to cell metadata", {

  expect_length(object@cell_metadata[["cell"]][["rna"]][["leiden_0.2"]], 73655)

  # test a few cluster assignments
  expect_equal(object@cell_metadata[["cell"]][["rna"]][["leiden_0.2"]][10], 5)
  expect_equal(object@cell_metadata[["cell"]][["rna"]][["leiden_0.2"]][80], 4)

})

# CELL TYPE MARKER GENE DETECTION
markers = findMarkers_one_vs_all(gobject = object,
                                 method = 'gini',
                                 expression_values = 'normalized',
                                 cluster_column = 'leiden_0.2',
                                 min_feats = 1, rank_score = 2,
                                 verbose = FALSE)

test_that("Cell type markers are detected", {

  # markers col names
  expect_named(markers, c("feats", "cluster", "expression", "expression_gini",
                          "detection", "detection_gini", "expression_rank",
                          "detection_rank", "comb_score", "comb_rank"))

  # number of markers
  expect_equal(nrow(markers), 584)

})

# CLUSTER ANNOTATION
selected_genes = c('Myh11', 'Klf4', 'Fn1', 'Cd24a', 'Cyr61', 'Nnat',
                   'Trh', 'Selplg', 'Pou3f2', 'Aqp4', 'Traf4',
                   'Pdgfra', 'Opalin', 'Mbp', 'Ttyh2', 'Fezf1',
                   'Cbln1', 'Slc17a6', 'Scg2', 'Isl1', 'Gad1')
cluster_order = c(6, 11, 9, 12, 4, 8, 7, 5, 13, 3, 1, 2, 10)
clusters_cell_types_hypo = c('Inhibitory', 'Inhibitory', 'Excitatory',
                             'Astrocyte','OD Mature', 'Endothelial',
                             'OD Mature', 'OD Immature', 'Ependymal', 'Ambiguous',
                             'Endothelial', 'Microglia', 'OD Mature')

names(clusters_cell_types_hypo) = as.character(sort(cluster_order))
object = annotateGiotto(gobject = object, annotation_vector = clusters_cell_types_hypo,
                        cluster_column = 'leiden_0.2', name = 'cell_types')

test_that("Cell type annotations are added to cell metadata", {

  expect_type(object@cell_metadata[["cell"]][["rna"]][["cell_types"]], "character")
  expect_length(object@cell_metadata[["cell"]][["rna"]][["cell_types"]], 73655)

  # check a few annotations
  expect_equal(object@cell_metadata[["cell"]][["rna"]][["cell_types"]][5], "Inhibitory")
  expect_equal(object@cell_metadata[["cell"]][["rna"]][["cell_types"]][250], "OD Immature")

})


# --------------------------------------------
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
