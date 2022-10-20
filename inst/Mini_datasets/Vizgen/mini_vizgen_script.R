
## MINI VIZGEN script and dataset preparation ##


devtools::load_all() #library(Giotto)
library(data.table)

# 0. preparation ####
# ----------------- #

## create instructions
instrs = createGiottoInstructions(save_dir = tempdir(),
                                  save_plot = FALSE,
                                  show_plot = TRUE,
                                  return_plot = FALSE)

## provide path to vizgen folder
data_path = system.file('/Mini_datasets/Vizgen/', package = 'Giotto')

## 0.1 path to images ####
# ---------------------- #

# vizgen creates data from 7 z-stack slices (z0, z1, ..., z6)
## - each z-slice has a DAPI image, polyT image
## - they also have a combined composite image, created from their segmentation kit (ab stainings)
DAPI_z0_image_path = paste0(data_path, '/', 'images/mini_dataset_dapi_z0.jpg')
DAPI_z1_image_path = paste0(data_path, '/', 'images/mini_dataset_dapi_z1.jpg')

polyT_z0_image_path = paste0(data_path, '/', 'images/mini_dataset_polyT_z0.jpg')
polyT_z1_image_path = paste0(data_path, '/', 'images/mini_dataset_polyT_z1.jpg')



## 0.2 path to transcripts ####
# --------------------------- #

## each transcript has x, y and z coordinate
tx_path = paste0(data_path, '/', 'vizgen_transcripts.gz')
tx_dt = data.table::fread(tx_path)



## 0.3 path to cell boundaries folder ####
# -------------------------------------- #

## vizgen already provides segmentation information in .hdf5 files
## the hdf5 files are organized in different tiles
## Here I have already converted the hdf5 files to a simple data.table format

boundary_path = paste0(data_path, '/cell_boundaries/')

z0_polygon_DT = fread(paste0(boundary_path, '/', 'z0_polygons.gz'))
z0_polygons = createGiottoPolygonsFromDfr(name = 'z0',
                                          segmdfr = z0_polygon_DT)

z1_polygon_DT = fread(paste0(boundary_path, '/', 'z1_polygons.gz'))
z1_polygons = createGiottoPolygonsFromDfr(name = 'z1',
                                          segmdfr = z1_polygon_DT)

# 1. create subcellular dataset with transcript and polygon information ####
# ------------------------------------------------------------------------ #
vizsubc = createGiottoObjectSubcellular(gpoints = list('rna' = tx_dt[,.(global_x, -global_y, gene, global_z)]),
                                        gpolygons = list(z0_polygons, z1_polygons),
                                        instructions = instrs)
showGiottoFeatInfo(vizsubc)
showGiottoSpatialInfo(vizsubc)


# calculate centroid for each polygon ( = cell)
# this can/will be used when aggregating for example counts to cells
vizsubc = addSpatialCentroidLocations(gobject = vizsubc,
                                      poly_info = paste0('z',0:1),
                                      return_gobject = T)
showGiottoSpatLocs(vizsubc)


# 2. add image ####
# --------------- #

# x and y information from original script
ultra_mini_extent = terra::ext(c(6400.029, 6900.037, -5150.007, -4699.967 ))

image_paths = c(DAPI_z0_image_path, DAPI_z1_image_path,
                polyT_z0_image_path, polyT_z1_image_path)
image_names = c('dapi_z0', 'dapi_z1',
                'polyT_z0', 'polyT_z1')

imagelist = createGiottoLargeImageList(raster_objects = image_paths,
                                       names = image_names,
                                       negative_y = TRUE,
                                       extent = ultra_mini_extent)

vizsubc = addGiottoImage(gobject = vizsubc,
                         largeImages = imagelist)

showGiottoImageNames(vizsubc)


# subset Giotto object based on locations
vizsubc = subsetGiottoLocsMulti(vizsubc,
                                spat_unit = c('z0', 'z1'),
                                poly_info = list(z0 = 'z0', z1 = 'z1'),
                                x_min = 6400.029,
                                x_max = 6900.037,
                                y_max = -4699.967,
                                y_min = -5150.007,
                                verbose = TRUE)

showGiottoSpatLocs(vizsubc)
showGiottoSpatialInfo(vizsubc)

# visualize
spatPlot2D(gobject = vizsubc,
           spat_unit = 'z0',
           show_image = T,
           largeImage_name = 'dapi_z0',
           point_shape = 'no_border',
           point_size = 2.5,
           point_alpha = 0.4,
           save_param = list(base_width = 7, base_height = 7))

spatPlot2D(gobject = vizsubc,
           spat_unit = 'z1',
           show_image = T,
           largeImage_name = 'polyT_z1',
           point_shape = 'no_border',
           point_size = 2.5,
           point_alpha = 0.4,
           save_param = list(base_width = 7, base_height = 7))


spatInSituPlotPoints(vizsubc,
                     feats = list('rna' = c("Htr1b", "Ackr1", "Epha7")),
                     feats_color_code = c("Htr1b" = 'green', 'Ackr1' = 'blue', 'Epha7' = 'red'),
                     show_image = TRUE,
                     largeImage_name = 'dapi_z0',
                     point_size = 0.5,
                     plot_method = 'ggplot',
                     show_polygon = TRUE,
                     use_overlap = F,
                     polygon_feat_type = 'z1',
                     polygon_color = 'yellow',
                     polygon_bg_color = 'yellow',
                     polygon_line_size = 0.2,
                     coord_fix_ratio = TRUE,
                     background_color = NA,
                     save_param = list(base_width = 7, base_height = 7))


# 3. aggregate information to matrix: polygons and transcripts ####
# --------------------------------------------------------------- #

# we will use the z1 polygon information
# we can set a global option or specify this for each command
# options('giotto.spat_unit' = 'z1') # now you don't need to think about setting spat_unit each time

vizsubc = calculateOverlapRaster(vizsubc,
                                 spatial_info = 'z0',
                                 feat_info = 'rna',
                                 feat_subset_column = 'global_z',
                                 feat_subset_ids = 0)

vizsubc = overlapToMatrix(vizsubc,
                          poly_info = 'z0',
                          feat_info = 'rna',
                          name = 'raw')

vizsubc = calculateOverlapRaster(vizsubc,
                                 spatial_info = 'z1',
                                 feat_info = 'rna',
                                 feat_subset_column = 'global_z',
                                 feat_subset_ids = 1)

vizsubc = overlapToMatrix(vizsubc,
                          poly_info = 'z1',
                          feat_info = 'rna',
                          name = 'raw')


showGiottoSpatialInfo(vizsubc)
showGiottoFeatInfo(vizsubc)

# 4. filter object on single layer #####
# ----------------------------------- ##
vizsubc <- filterGiotto(gobject = vizsubc,
                        spat_unit = 'z1',
                        expression_threshold = 1,
                        feat_det_in_min_cells = 3,
                        min_det_feats_per_cell = 5, poly_info = 'z1')



# 5. normalize on single layer #####
# ----------------------------- ##

# rna data, default.
# other feature modalities can be processed and filtered in an anologous manner
vizsubc <- normalizeGiotto(gobject = vizsubc, spat_unit = 'z1', scalefactor = 5000, verbose = T)
vizsubc <- addStatistics(gobject = vizsubc, spat_unit = 'z1')
vizsubc <- normalizeGiotto(gobject = vizsubc, spat_unit = 'z1', norm_methods = 'pearson_resid', update_slot = 'pearson')


spatPlot2D(gobject = vizsubc, spat_unit = 'z1',
           cell_color = 'total_expr', color_as_factor = F,
           largeImage_name = 'dapi_z1', show_image = TRUE,
           point_size = 3.5, point_alpha = 0.5, coord_fix_ratio = T)

spatInSituPlotPoints(vizsubc,
                     show_polygon = TRUE,
                     polygon_feat_type = 'z1',
                     polygon_color = 'white',
                     polygon_line_size = 0.1,
                     polygon_fill = 'total_expr',
                     polygon_fill_as_factor = F,
                     coord_fix_ratio = T)



# 6. highly variable genes ####
# ----------------------------- #

# typical way of calculating HVG
vizsubc <- calculateHVF(gobject = vizsubc, spat_unit = 'z1', HVFname = 'hvg_orig')

# new method based on variance of pearson residuals for each gene
vizsubc <- calculateHVF(gobject = vizsubc, spat_unit = 'z1',
                        method = 'var_p_resid', expression_values = 'pearson',
                        show_plot = T)

gene_meta = fDataDT(vizsubc, spat_unit = 'z1')


# 7. dimension reduction ####
# --------------------------- #

# ** 7.1 PCA ####

# we will run pca on the pre-scaled matrix from the pearson residual normalization
# if features are not specified it will automatically search for the hvf column in the feature metadata

vizsubc <- runPCA(gobject = vizsubc,
                  spat_unit = 'z1',
                  expression_values = 'pearson',
                  scale_unit = F, center = F)

screePlot(vizsubc,
          ncp = 20,
          spat_unit = 'z1')

showGiottoDimRed(vizsubc)

plotPCA(vizsubc,
        spat_unit = 'z1',
        dim_reduction_name = 'pca',
        dim1_to_use = 1,
        dim2_to_use = 2)



# ** 7.2 UMAP and TSN ####
vizsubc <- runUMAP(vizsubc, dimensions_to_use = 1:8, n_threads = 4, spat_unit = 'z1')
plotUMAP(gobject = vizsubc, spat_unit = 'z1')


# 8. graph-based clustering ####
# ---------------------------- #
vizsubc <- createNearestNetwork(gobject = vizsubc, dimensions_to_use = 1:8, k = 10, spat_unit = 'z1')
vizsubc <- doLeidenCluster(gobject = vizsubc, resolution = 0.05, n_iterations = 1000, spat_unit = 'z1')

# visualize UMAP cluster results
plotUMAP(gobject = vizsubc, spat_unit = 'z1',
         cell_color = 'leiden_clus',
         show_NN_network = T, point_size = 2.5)

spatInSituPlotPoints(vizsubc,
                     show_image = T,
                     largeImage_name = 'dapi_z0',
                     feats = list('rna' = c("Htr1b", "Ackr1", "Epha7")),
                     feats_color_code = c("Htr1b" = 'green', 'Ackr1' = 'blue', 'Epha7' = 'red'),
                     point_size = 0.35,
                     show_polygon = TRUE,
                     polygon_feat_type = 'z1',
                     polygon_color = 'white',
                     polygon_line_size = 0.1,
                     polygon_fill = 'leiden_clus',
                     polygon_fill_as_factor = T,
                     coord_fix_ratio = TRUE)


## 9. save Giotto object ####
# ------------------------- #

format(object.size(vizsubc), units = 'Mb')

#saveRDS(vizsubc, file = paste0(data_path, '/', 'gobject_mini_vizgen.RDS'))


# write terra files (spatvectors)
# these need to be read in again when loading the vizgen giotto object
terra::writeVector(vizsubc@feat_info$rna@spatVector, filename = paste0(data_path, '/', 'processed_data/rna_spatVector.shp'))

terra::writeVector(vizsubc@spatial_info$z0@spatVector, filename = paste0(data_path, '/', 'processed_data/z0_spatVector.shp'))
terra::writeVector(vizsubc@spatial_info$z0@spatVectorCentroids, filename = paste0(data_path, '/', 'processed_data/z0_spatVectorCentroids.shp'))

terra::writeVector(vizsubc@spatial_info$z1@spatVector, filename = paste0(data_path, '/', 'processed_data/z1_spatVector.shp'))
terra::writeVector(vizsubc@spatial_info$z1@spatVectorCentroids, filename = paste0(data_path, '/', 'processed_data/z1_spatVectorCentroids.shp'))


# change name to fix error


# set pointers to NULL to fix potential error
vizsubc@feat_info$rna@spatVector = 1

vizsubc@spatial_info$z0@spatVector = 1
vizsubc@spatial_info$z0@spatVectorCentroids = 1

vizsubc@spatial_info$z1@spatVector = 1
vizsubc@spatial_info$z1@spatVectorCentroids = 1

vizsubc@largeImages = NULL

# save object and copy object to Giotto/data folder
gobject_mini_vizgen = vizsubc
save(gobject_mini_vizgen, file = paste0(data_path, '/', 'gobject_mini_vizgen.rda'))





