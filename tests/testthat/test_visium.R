python_path = NULL 
if(is.null(python_path)) {
  installGiottoEnvironment()
}

filtered_matrix = system.file("extdata/visium_brain_data", "filtered_feature_bc_matrix.h5", package = 'Giotto')
tissue_positions = system.file("extdata/visium_brain_data/spatial", "tissue_positions_list.csv", package = 'Giotto')
image_path = system.file("extdata/visium_brain_data/spatial", "tissue_lowres_image.png", package = 'Giotto')
scalefactors_path = system.file("extdata/visium_brain_data/spatial", "scalefactors_json.json", package = "Giotto")

### TESTS FOR 10X VISIUM MOUSE BRAIN DATASET
# --------------------------------------------------------------

# CREATE VISIUM OBJECT
object = createGiottoVisiumObject(h5_visium_path = filtered_matrix,
                                  h5_tissue_positions_path = tissue_positions,
                                  h5_image_png_path = image_path,
                                  h5_json_scalefactors_path = scalefactors_path)

test_that("Object initialization creates expected Giotto object", {
  
  # S4 object of class "giotto"
  expect_s4_class(object, "giotto")
  
  # gobject contains S4 object of class "dgCMatrix" containing raw expression
  expect_s4_class(object@expression[["cell"]][["rna"]][["raw"]], "dgCMatrix")
  expect_true(all(object@expression[["cell"]][["rna"]][["raw"]]@Dim == c(32285, 2702)))
  
  # gobject contains S3 object "data.table" containing spatial locations
  expect_s3_class(object@spatial_locs[["cell"]][["raw"]], "data.table")
  expect_length(object@spatial_locs[["cell"]][["raw"]][["sdimx"]], 2702)
  expect_length(object@spatial_locs[["cell"]][["raw"]][["sdimy"]], 2702)
  expect_length(object@spatial_locs[["cell"]][["raw"]][["cell_ID"]], 2702)
  
  # metadata length/types
  expect_length(object@cell_metadata[["cell"]][["rna"]][["cell_ID"]], 2702)
  expect_type(object@cell_metadata[["cell"]][["rna"]][["cell_ID"]], "character")
  expect_length(object@cell_metadata[["cell"]][["rna"]][["in_tissue"]], 2702)
  expect_type(object@cell_metadata[["cell"]][["rna"]][["in_tissue"]], "integer")
  expect_length(object@feat_metadata[["cell"]][["rna"]][["feat_ID"]], 32285)
  expect_type(object@feat_metadata[["cell"]][["rna"]][["feat_ID"]], "character")
  
})

# PROCESS GIOTTO VISIUM OBJECT
metadata = pDataDT(object)
in_tissue_barcodes = metadata[in_tissue == 1]$cell_ID
object = subsetGiotto(object, cell_ids = in_tissue_barcodes)

test_that("Object filtered for  in-tissue", {
  expect_true(all(object@cell_metadata[["cell"]][["rna"]][["in_tissue"]]))
})

object <- filterGiotto(gobject = object,
                       expression_threshold = 1,
                       feat_det_in_min_cells = 50,
                       min_det_feats_per_cell = 1000,
                       expression_values = c('raw'),
                       verbose = F)

object <- normalizeGiotto(gobject = object, scalefactor = 6000, verbose = F)

object <- addStatistics(gobject = object)

# DIMENSION REDUCTION
object <- calculateHVF(gobject = object, save_plot = FALSE)

gene_metadata = fDataDT(object)
featgenes = gene_metadata[hvf == 'yes' & perc_cells > 3 & mean_expr_det > 0.4]$feat_ID

object <- runPCA(gobject = object,
                 feats_to_use = featgenes)

object <- runUMAP(object, dimensions_to_use = 1:10)

object <- runtSNE(object, dimensions_to_use = 1:10)

# create network
object <- createNearestNetwork(gobject = object, dimensions_to_use = 1:10, k = 15)

# cluster
object <- doLeidenCluster(gobject = object, resolution = 0.4, n_iterations = 1000)

# SUBSET
DG_subset = subsetGiottoLocs(object,
                             x_max = 6500, x_min = 3000,
                             y_max = -2500, y_min = -5500,
                             return_gobject = TRUE)


test_that("Subset object created as expected", {
  
  # DG_subset contains S4 object of class "dgTMatrix" containing raw expression
  expect_s4_class(DG_subset@expression[["cell"]][["rna"]][["raw"]], "dgCMatrix")
  expect_true(all(DG_subset@expression[["cell"]][["rna"]][["raw"]]@Dim == c(14814, 624)))
  
  # DG_subset contains S4 object of class "dgTMatrix" containing raw expression
  expect_s4_class(DG_subset@expression[["cell"]][["rna"]][["normalized"]], "dgCMatrix")
  expect_true(all(DG_subset@expression[["cell"]][["rna"]][["normalized"]]@Dim == c(14814, 624)))
  
  # DG_subset contains S4 object of class "dgTMatrix" containing raw expression
  expect_s4_class(DG_subset@expression[["cell"]][["rna"]][["scaled"]], "dgeMatrix")
  expect_true(all(DG_subset@expression[["cell"]][["rna"]][["scaled"]]@Dim == c(14814, 624)))
  
  # DG_subset contains S3 object "data.table" containing spatial locations
  expect_s3_class(DG_subset@spatial_locs[["cell"]][["raw"]], "data.table")
  expect_length(DG_subset@spatial_locs[["cell"]][["raw"]][["sdimx"]], 624)
  expect_length(DG_subset@spatial_locs[["cell"]][["raw"]][["sdimy"]], 624)
  expect_length(DG_subset@spatial_locs[["cell"]][["raw"]][["cell_ID"]], 624)
  
  # dimObj containing dimension reduction points
  expect_s3_class(DG_subset@dimension_reduction[["cells"]][["cell"]][["rna"]][["pca"]][["pca"]], "dimObj")
  expect_s3_class(DG_subset@dimension_reduction[["cells"]][["cell"]][["rna"]][["umap"]][["umap"]], "dimObj")
  expect_s3_class(DG_subset@dimension_reduction[["cells"]][["cell"]][["rna"]][["tsne"]][["tsne"]], "dimObj")
  
})




