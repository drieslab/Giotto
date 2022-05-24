### Jocelyn Test ##
library(Giotto)
# create data
# ST SCC ####
ST_SCC_data = list(
  dataset = 'ST_SCC',
  spatial_locs = c("https://github.com/RubD/spatial-datasets/raw/master/data/2020_ST_SCC/cell_locations/P2_1_spatial_locs.csv",
                   "https://github.com/RubD/spatial-datasets/raw/master/data/2020_ST_SCC/cell_locations/P2_2_spatial_locs.csv",
                   "https://github.com/RubD/spatial-datasets/raw/master/data/2020_ST_SCC/cell_locations/P2_3_spatial_locs.csv"),
  expr_matrix = c("https://github.com/RubD/spatial-datasets/raw/master/data/2020_ST_SCC/count_matrix/P2_1_expression.csv",
                  "https://github.com/RubD/spatial-datasets/raw/master/data/2020_ST_SCC/count_matrix/P2_2_expression.csv",
                  "https://github.com/RubD/spatial-datasets/raw/master/data/2020_ST_SCC/count_matrix/P2_3_expression.csv"),
  metadata = c("https://github.com/RubD/spatial-datasets/raw/master/data/2020_ST_SCC/scRNAseq/p2_metadata.rds",
               "https://github.com/RubD/spatial-datasets/raw/master/data/2020_ST_SCC/scRNAseq/scrna_expr.rds",
               "https://github.com/RubD/spatial-datasets/raw/master/data/2020_ST_SCC/scRNAseq/normalized_sc_matrix.RDS",
               "https://github.com/RubD/spatial-datasets/raw/master/data/2020_ST_SCC/scRNAseq/cell_type_vector.RDS",
               "https://github.com/RubD/spatial-datasets/raw/master/data/2020_ST_SCC/scRNAseq/sign_list.RDS",
               "https://github.com/RubD/spatial-datasets/raw/master/data/2020_ST_SCC/lowres_images/P2_1_0.0625.jpg",
               "https://github.com/RubD/spatial-datasets/raw/master/data/2020_ST_SCC/lowres_images/P2_2_0.0625.jpg",
               "https://github.com/RubD/spatial-datasets/raw/master/data/2020_ST_SCC/lowres_images/P2_3_0.0625.jpg",
               "https://github.com/RubD/spatial-datasets/raw/master/data/2020_ST_SCC/imgReg.zip",
               "https://github.com/RubD/spatial-datasets/raw/master/data/2020_ST_SCC/PairsLigRec.txt")
)
########### Downloading Data ################
#spatial_locations
url <- "https://github.com/RubD/spatial-datasets/raw/master/data/2020_ST_SCC/cell_locations/P2_2_spatial_locs.csv"
destfile <- "/Users/adrianasistig/Documents/GitHub/Giotto/inst/extdata/ST_SCC/P2_2_spatial_locs.csv"
download.file(url,destfile)
#expression file
url <- "https://github.com/RubD/spatial-datasets/raw/master/data/2020_ST_SCC/count_matrix/P2_2_expression.csv"
destfile <- "/Users/adrianasistig/Documents/GitHub/Giotto/inst/extdata/ST_SCC/P2_2_expression.csv"
download.file(url,destfile)

#image
url <- "https://github.com/RubD/spatial-datasets/raw/master/data/2020_ST_SCC/lowres_images/P2_2_0.0625.jpg"
destfile <- "/Users/adrianasistig/Documents/GitHub/Giotto/inst/extdata/ST_SCC/P2_2_0.0625.jpg"
download.file(url,destfile)

########## Jocelyn Test Script ###################

# set data and results directories
data_directory <- "inst/extdata/ST_SCC"
save_directory <-  "save_dir"

# create instructions
my_instructions <- createGiottoInstructions(save_plot = TRUE,
                                            save_dir = save_directory)
# create giotto object
my_giotto_object <- createGiottoObject(expression = paste0(data_directory, "/P2_2_expression.csv"),
                                       spatial_locs = paste0(data_directory, "/P2_2_spatial_locs.csv"),instructions = my_instructions)
# create image
my_giotto_image <- createGiottoImage(gobject = my_giotto_object,
                                     do_manual_adj = FALSE,
                                     scale_factor = 0.0625,
                                     mg_object = magick::image_read(paste0(data_directory, "/P2_2_0.0625.jpg")))
# add image
my_giotto_object <- addGiottoImage(gobject = my_giotto_object,
                                   images = list(my_giotto_image),
                                   spat_loc_name = "raw")
​
## filter
my_giotto_object <- filterGiotto(gobject = my_giotto_object,
                                 spat_unit = "cell",
                                 feat_type = "rna",
                                 expression_threshold = 2,
                                 feat_det_in_min_cells = 2,
                                 min_det_feats_per_cell = 100)
​
## normalize expression matrix
my_giotto_object <- normalizeGiotto(gobject = my_giotto_object,
                                    spat_unit = "cell",
                                    feat_type = "rna",
                                    norm_methods = "standard",
                                    scalefactor = 6000,
                                    verbose = TRUE)
​
## add gene statistics
my_giotto_object <- addStatistics(gobject = my_giotto_object)
​
## clustering and cell type identification
my_giotto_object <- calculateHVF(gobject = my_giotto_object,
                                 spat_unit = "cell",
                                 feat_type = "rna",
                                 expression_values = "normalized",
                                 method = "cov_groups",
                                 nr_expression_groups = 20,
                                 zscore_threshold = 1.5, save_plot = T,
                                 save_param = list(save_name = 'fig8_calculateHVG',
                                                   base_width = 6,
                                                   base_height = 6,
                                                   save_format = 'pdf'))
​
## Principal component analysis (PCA)
my_giotto_object <- runPCA(gobject = my_giotto_object,
                           spat_unit = "cell",
                           feat_type = "rna",
                           expression_values = "normalized",
                           reduction = "cells",
                           feats_to_use = "hvf")
​
## Uniform manifold approximation projection (UMAP)
my_giotto_object <- runUMAP(gobject = my_giotto_object,
                            spat_unit = "cell",
                            feat_type = "rna",
                            dimensions_to_use = 1:10)

## Create nearest network
my_giotto_object <- createNearestNetwork(gobject = my_giotto_object,
                                         spat_unit = "cell",
                                         feat_type = "rna",
                                         dimensions_to_use = 1:10)
## Leiden clustering
my_giotto_object <- doLeidenCluster(gobject = my_giotto_object,
                                    spat_unit = "cell",
                                    feat_type = "rna",
                                    name = "leiden_clus")
​
## cell type marker gene detection
gini_markers_subclusters <- findMarkers_one_vs_all(gobject = my_giotto_object,
                                                   spat_unit = "cell",
                                                   feat_type = "rna",
                                                   method = 'gini',
                                                   expression_values = 'normalized',
                                                   cluster_column = 'leiden_clus',
                                                   min_featss = 20,
                                                   min_expr_gini_score = 0.5,
                                                   min_det_gini_score = 0.5)
topgenes_gini <- gini_markers_subclusters[, head(.SD, 2), by = 'cluster']$feats
​
scran_results <- findMarkers_one_vs_all(my_giotto_object,
                                        spat_unit = "cell",
                                        feat_type = "rna",
                                        method = "scran",
                                        expression_values = 'normalized',
                                        cluster_column = 'leiden_clus',
                                        #subset_clusters = c("polygon 1", "polygon 2"),
                                        min_feats = 10
)
scran_results
topgenes_scran <- scran_results[, head(.SD, 2), by = 'cluster']$feats
topgenes_scran
# "IGFBP4"  "COL1A2"  "GJB2"    "PERP"    "S100A2"  "FABP5"   "HLA-B"   "HLA-A"
# "KRT5"    "S100A2"  "VIM"     "COL18A1" "DCD"     "MUCL1"   "KRT1"    "KRT10"

## umap plots using marker genes
dimFeatPlot2D(my_giotto_object,
              expression_values = 'scaled',
              feats = gini_markers_subclusters[, head(.SD, 1), by = 'cluster']$feats,
              cow_n_col = 3, point_size = 1,
              save_plot = FALSE)

spat_unit = set_default_spat_unit(mini_giotto_single_cell)
giotto_object_test <- createGiottoObject(mini_giotto_single_cell)

## Adriana ##
## Filter Giotto Object ##

filtered_gobject = filterGiotto(my_giotto_object,
                                spat_unit = "cell",
                                feat_type = "rna",
                                feat_det_in_min_cells = 10,
                                min_det_feats_per_cell = 10)

## Subset Giotto Function ##
random_cells = sample(slot(my_giotto_object, 'cell_ID'), 10, replace = TRUE)
random_genes = sample(slot(my_giotto_object, 'feat_ID'), 10, replace = TRUE)


subset_obj = subsetGiotto(my_giotto_object,
                          cell_ids = random_cells,
                          feat_ids = random_genes)
subset_obj

## Subset based on spatial locations ##
subset_spat_locations <- subsetGiottoLocs(my_giotto_object)
subset_spat_locations

## filter distributions ##
filtered_distr <- filterDistributions(my_giotto_object, detection = 'feats')

## filter combinations ##
filtered_gobject = filterGiotto(my_giotto_object,
  spat_unit = "cell",
  feat_type = "rna",
  feat_det_in_min_cells = 10,
  min_det_feats_per_cell = 10)
#filtered_gobject

## rna standard normalization ##
rna_standard_norm <- rna_standard_normalization(my_gioto_object)

## normalize Giotto ##
normGiotto <- normalizeGiotto(my_giotto_object)

## adjust Giotto Matrix ##

adjMat <- adjustGiottoMatrix(my_giotto_object,
                             spat_unit = 'cell',
                             feat_type = 'rna',
                             update_slot = 'normalized')
adjMat
## process Giotto ##

proGiotto <- processGiotto(my_giotto_object)

## annotate Giotto ##
#create a vector with cell type nanes as names of the vector
cluster_cell_types = c('cell_type_1','cell_type_2','cell_type_3','cell_type_4',
                       'cell_type_5','cell_type_6','cell_type_7','cell_type_8')
names(cluster_cell_types) = 1:8

# leiden clustering results
cell_metadata = pDataDT(my_giotto_object)

#show clustering
cell_metadata[['leiden_clus']]

# convert the cluster results into annotations and add the cell metadata
annGiotto <- annotateGiotto(my_giotto_object,
                            annotation_vector = cluster_cell_types,
                            cluster_column = 'leiden_clus',
                            name = 'cell_type2')
# show annotation results
spatDimPlot(my_giotto_object,
            cell_color = 'cell_type2',
            spat_point_size = 3,
            dim_point_size = 2)

## remove cell annotation ##
# cell metadata
cell_metadata = pDataDT(my_giotto_object)

# remove the cell_type column
my_giotto_object <- removeCellAnnotation(annGiotto, columns = 'cell_type2')

## remove feature annotation ##
## show feature metadata ##
######## ERROR remove feature annotation ########
annotateGiotto(my_giotto_object, columns = 'feat_ID')

## add cell metadata ##
# leiden clustering results
cell_metadata = pDataDT(my_giotto_object)

added_metadata <- addCellMetadata(my_giotto_object, new_metadata = cell_metadata)

## add gene metadata ##
gene_metadata = fDataDT(my_giotto_object)
added_metadata <- addCellMetadata(my_giotto_object, new_metadata = gene_metadata)

## add gene statistics ##
added_gene_stats <- addFeatStatistics(my_giotto_object)

## add cell statistics ##
added_cell_stats <- addCellStatistics(my_giotto_object)

## add statistics ##
add_stats <- addStatistics(my_giotto_object)

## add Genes Perc ##
# select features
random_feats = sample(slot(my_giotto_object, 'feat_ID'), 10, replace = TRUE)

#calculate percentage of those selected features per cell/spot
perc_calc <- addFeatsPerc(my_giotto_object, feats = random_feats,
                          vector_name = 'random_feat_perc')
fDataDT(perc_calc)

## show processing steps ##
steps <- showProcessingSteps(my_giotto_object)

## calculate meta table ##

## metadata
## add feat metadata ##
feat_metadata = fDataDT(my_giotto_object)
added_metadata <- addFeatMetadata(my_giotto_object, new_metadata = feat_metadata)
# show average feat expression per annotated cell type
avg_feat_expr <- calculateMetaTable(added_metadata, metadata_cols = 'feat_ID')
avg_feat_expr

## calculateMetaTableCells ##
######### ERROR calculateMetaTableCells ############
cell_metadata <- pDataDT(my_giotto_object)
colnames(cell_metadata)
# [1] "cell_ID"     "leiden_clus" "nb_cells"    "cell_type2"
added_cell_metadata <- addCellMetadata(my_giotto_object, new_metadata = cell_metadata)

# show average cell expression per annotated cell type
avg_cell_expr <- calculateMetaTableCells(added_cell_metadata, metadata_cols = 'cell_ID')
##### ERROR MESSAGE ####
# Error in calculateMetaTableCells(added_cell_metadata, metadata_cols = "cell_ID") :
#
#   You need to select one or more valid value column names from pDataDT()

.


## combine metadata ##
cell_metadata <- pDataDT(my_giotto_object)
addCellMetadata(my_giotto_object, new_metadata = cell_metadata)
combined_metadata <- combineMetadata(my_giotto_object)
combined_metadata

## combine meta feats ##
####### ERROR ###########
# get all features
all_feats <- slot(my_giotto_object, 'feat_ID')

# create 3 metagenes from the first 6 genes
clust_vector <- c(1,1,1,2,2,2,3,3,3) # 3 groups
names(clust_vector) = all_feats$rna[2:10]

meta_feats <- createMetafeats(my_giotto_object,
                              feat_clusters = clust_vector,
                              name = 'cluster_metagene')

#plot metagene expression
spatFeatPlot2D(meta_feats,
               feat_type = 'rna',
               feats = 'feat_ID',
               point_size = 3,
               cow_n_col = 2)

spatFeatPlot2D_single(meta_feats,
               spat_enr_names = 'cluster_metagene',
               cell_annotation_values = c('1','2'),
               point_size = 3,
               cow_n_col = 2)

## find network neighbours ##
# get all cells
all_cells <- slot(my_giotto_object, 'cell_ID')
# create network
SpatNet <- createSpatialDelaunayNetwork(my_giotto_object)
# make sure that network is created
showGiottoSpatNetworks(SpatNet)
# find all the spatial neighbours for the first 10 cells within the Delaunay network
find_spat_neighbor <- findNetworkNeighbors(SpatNet,
                                           spatial_network_name = 'Delaunay_network',
                                           source_cell_ids = all_cells$cell[1:10])

# > find_spat_neighbor
# cell_ID nb_cells
# 1:    8x20     both
# 2:    8x22     both
# 3:    8x32     both
# 4:    8x34     both
# 5:    8x36     both
# ---
#   629:   49x31   others
# 630:   49x33   others
# 631:   50x28   others
# 632:   50x30   others
# 633:   50x32   others

## combine spatial cell feature information ##
####### ERROR combinespatialcellfeatureinfo ######
# get spatial locations
spat_loc <- get_spatial_locations(my_giotto_object)
combine_spatial_info <- combineSpatialCellFeatureInfo(spat_loc, feat_type = "cell_ID")

## fDataDT ##
feat <- fDataDT(my_giotto_object)
feat[1:10]
# feat_ID hvf
# 1: FO538757.1  no
# 2:      NOC2L  no
# 3:     KLHL17  no
# 4:    PLEKHN1  no
# 5:       HES4 yes
# 6:      ISG15  no
# 7:       AGRN  no
# 8:     RNF223  no
# 9:   C1orf159  no
# 10:   TNFRSF18  no

## pDataDT ##
cell <- pDataDT(my_giotto_object)
cell[1:10]
# cell_ID leiden_clus nb_cells
# 1:    8x20           4     both
# 2:    8x22           3     both
# 3:    8x32           1     both
# 4:    8x34           1     both
# 5:    8x36           3     both
# 6:    9x19           5     both
# 7:    9x21           4     both
# 8:    9x23           3     both
# 9:    9x25           6     both
# 10:    9x27           1     both

## createmetafeats ##

meta_features <- createMetafeats(my_giotto_object, feat_clusters = )

## addFeatMetadata ##
## add feat metadata ##
feat_metadata = fDataDT(my_giotto_object)
added_metadata <- addFeatMetadata(my_giotto_object, new_metadata = feat_metadata)
added_metadata
