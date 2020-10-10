

library(Giotto)

# createGiottoInstructions(python_path = '/your/path')

## Giotto 0.3.5 ##
## mini-test Visium Brain Giotto 0.3.5 ##

## !! change this directory path !!:
#temp_dir = '/path/to/your/temporary/directory/results'
temp_dir = getwd()

## 1. giotto object ####
expr_path = system.file("extdata", "visium_DG_expr.txt", package = 'Giotto')
loc_path = system.file("extdata", "visium_DG_locs.txt", package = 'Giotto')

# default
brain_small <- createGiottoObject(raw_exprs = expr_path, spatial_locs = loc_path)

# with instructions (e.g. specific python path)
# myinstructions = createGiottoInstructions(python_path = '/your/path')
# brain_small <- createGiottoObject(raw_exprs = expr_path, spatial_locs = loc_path, instructions = myinstructions)
showGiottoInstructions(brain_small)


## 1. read image
png_path = system.file("extdata", "deg_image.png", package = 'Giotto')
mg_img = magick::image_read(png_path)

## 2. test and modify image alignment
mypl = spatPlot(brain_small, return_plot = T, point_alpha = 0.8)
orig_png = createGiottoImage(gobject = brain_small, mg_object = mg_img, name = 'image',
                             xmax_adj = 450, xmin_adj = 550,
                             ymax_adj = 200, ymin_adj = 200)
mypl_image = addGiottoImageToSpatPlot(mypl, orig_png)
mypl_image

## 3. add images to Giotto object ##
image_list = list(orig_png)
brain_small = addGiottoImage(gobject = brain_small,
                             images = image_list)
showGiottoImageNames(brain_small)





## 2. processing steps ####
filterDistributions(brain_small, detection = 'genes')
filterDistributions(brain_small, detection = 'cells')
filterCombinations(brain_small,
                   expression_thresholds = c(1),
                   gene_det_in_min_cells = c(20, 20, 50, 50),
                   min_det_genes_per_cell = c(100, 200, 100, 200))

brain_small <- filterGiotto(gobject = brain_small,
                            expression_threshold = 1,
                            gene_det_in_min_cells = 50,
                            min_det_genes_per_cell = 100,
                            expression_values = c('raw'),
                            verbose = T)
brain_small <- normalizeGiotto(gobject = brain_small, scalefactor = 6000, verbose = T)
brain_small <- addStatistics(gobject = brain_small)


## 3. dimension reduction ####
brain_small <- calculateHVG(gobject = brain_small)
brain_small <- runPCA(gobject = brain_small)
screePlot(brain_small, ncp = 30)
plotPCA(gobject = brain_small)
brain_small <- runUMAP(brain_small, dimensions_to_use = 1:10)
plotUMAP(gobject = brain_small)
brain_small <- runtSNE(brain_small, dimensions_to_use = 1:10)
plotTSNE(gobject = brain_small)


## 4. clustering ####
brain_small <- createNearestNetwork(gobject = brain_small, dimensions_to_use = 1:10, k = 15)
brain_small <- doLeidenCluster(gobject = brain_small, resolution = 0.4, n_iterations = 1000)
plotUMAP(gobject = brain_small, cell_color = 'leiden_clus', show_NN_network = T, point_size = 2.5)


## 5. co-visualize ####
spatDimPlot(gobject = brain_small, cell_color = 'leiden_clus',
            dim_point_size = 2, spat_point_size = 2.5)



## 6. differential expression ####
scran_markers = findMarkers_one_vs_all(gobject = brain_small,
                                       method = 'scran',
                                       expression_values = 'normalized',
                                       cluster_column = 'leiden_clus')
# violinplot
topgenes_scran = scran_markers[, head(.SD, 2), by = 'cluster']$genes
violinPlot(brain_small, genes = topgenes_scran, cluster_column = 'leiden_clus',
           strip_text = 10, strip_position = 'right')

# metadata heatmap
topgenes_scran = scran_markers[, head(.SD, 6), by = 'cluster']$genes
plotMetaDataHeatmap(brain_small, selected_genes = topgenes_scran,
                    metadata_cols = c('leiden_clus'))



## 7. annotation ####

### 7.1 spot annotation ####
clusters_cell_types = c('Gfap_cells', 'Tbr1_cells', 'Tcf7l2_cells', 'Wfs1_cells', 'Nptxr_cells')
names(clusters_cell_types) = 1:5
brain_small = annotateGiotto(gobject = brain_small, annotation_vector = clusters_cell_types,
                             cluster_column = 'leiden_clus', name = 'cell_types')
spatDimPlot(gobject = brain_small, cell_color = 'cell_types', spat_point_size = 3, dim_point_size = 3)




### 7.2 spatial enrichment: cell type distribution ####

# known markers for different mouse brain cell types:
# Zeisel, A. et al. Molecular Architecture of the Mouse Nervous System. Cell 174, 999-1014.e22 (2018).

## cell type signatures ##
## combination of all marker genes identified in Zeisel et al
sign_matrix_path = system.file("extdata", "sig_matrix.txt", package = 'Giotto')
brain_sc_markers = data.table::fread(sign_matrix_path) # file don't exist in data folder
sig_matrix = as.matrix(brain_sc_markers[,-1]); rownames(sig_matrix) = brain_sc_markers$Event

## enrichment tests
brain_small = createSpatialEnrich(brain_small, sign_matrix = sig_matrix, enrich_method = 'PAGE') #default = 'PAGE'

## heatmap of enrichment versus annotation (e.g. clustering result)
cell_types = colnames(sig_matrix)
plotMetaDataCellsHeatmap(gobject = brain_small,
                         metadata_cols = 'leiden_clus',
                         value_cols = cell_types,
                         spat_enr_names = 'PAGE',
                         x_text_size = 8, y_text_size = 8)


enrichment_results = brain_small@spatial_enrichment$PAGE
enrich_cell_types = colnames(enrichment_results)
enrich_cell_types = enrich_cell_types[enrich_cell_types != 'cell_ID']

## spatplot
spatCellPlot(gobject = brain_small, spat_enr_names = 'PAGE',
             cell_annotation_values = enrich_cell_types,
             cow_n_col = 3,coord_fix_ratio = NULL, point_size = 1)



## 8. spatial grid ####
brain_small <- createSpatialGrid(gobject = brain_small,
                                 sdimx_stepsize = 300,
                                 sdimy_stepsize = 300,
                                 minimum_padding = 50)
showGrids(brain_small)
spatPlot(gobject = brain_small, show_grid = T, point_size = 1.5)

annotated_grid = annotateSpatialGrid(brain_small)
annotated_grid_metadata = annotateSpatialGrid(brain_small,
                                              cluster_columns = c('leiden_clus', 'cell_types', 'nr_genes'))


## 9. spatial network ####
plotStatDelaunayNetwork(gobject = brain_small, maximum_distance = 300)
brain_small = createSpatialNetwork(gobject = brain_small, minimum_k = 2, maximum_distance_delaunay = 400)
brain_small = createSpatialNetwork(gobject = brain_small, minimum_k = 2, method = 'kNN', k = 10)
showNetworks(brain_small)


spatPlot(gobject = brain_small, show_network = T,
         network_color = 'blue', spatial_network_name = 'Delaunay_network',
         point_size = 2.5, cell_color = 'leiden_clus')

spatPlot(gobject = brain_small, show_network = T,
         network_color = 'blue', spatial_network_name = 'kNN_network',
         point_size = 2.5, cell_color = 'leiden_clus')



## 10. spatial genes ####
km_spatialgenes = binSpect(brain_small)
spatGenePlot(brain_small, expression_values = 'scaled', genes = km_spatialgenes[1:6]$genes,
             point_shape = 'border', point_border_stroke = 0.1,
             show_network = F, network_color = 'lightgrey', point_size = 2.5,
             cow_n_col = 2)

rank_spatialgenes = binSpect(brain_small, bin_method = 'rank')
spatGenePlot(brain_small, expression_values = 'scaled', genes = rank_spatialgenes[1:6]$genes,
             point_shape = 'border', point_border_stroke = 0.1,
             show_network = F, network_color = 'lightgrey', point_size = 2.5,
             cow_n_col = 2)

silh_spatialgenes = silhouetteRank(gobject = brain_small) # TODO: suppress print output
spatGenePlot(brain_small, expression_values = 'scaled', genes = silh_spatialgenes[1:6]$genes,
             point_shape = 'border', point_border_stroke = 0.1,
             show_network = F, network_color = 'lightgrey', point_size = 2.5,
             cow_n_col = 2)



## 11. spatial co-expression patterns ####
ext_spatial_genes = km_spatialgenes[1:100]$genes
spat_cor_netw_DT = detectSpatialCorGenes(brain_small,
                                         method = 'network', spatial_network_name = 'Delaunay_network',
                                         subset_genes = ext_spatial_genes)

spat_cor_netw_DT = clusterSpatialCorGenes(spat_cor_netw_DT, name = 'spat_netw_clus', k = 8)
heatmSpatialCorGenes(brain_small, spatCorObject = spat_cor_netw_DT, use_clus_name = 'spat_netw_clus')

netw_ranks = rankSpatialCorGroups(brain_small, spatCorObject = spat_cor_netw_DT, use_clus_name = 'spat_netw_clus')
top_netw_spat_cluster = showSpatialCorGenes(spat_cor_netw_DT, use_clus_name = 'spat_netw_clus',
                                            selected_clusters = 6, show_top_genes = 1)

cluster_genes_DT = showSpatialCorGenes(spat_cor_netw_DT, use_clus_name = 'spat_netw_clus', show_top_genes = 1)
cluster_genes = cluster_genes_DT$clus; names(cluster_genes) = cluster_genes_DT$gene_ID
brain_small = createMetagenes(brain_small, gene_clusters = cluster_genes, name = 'cluster_metagene')
spatCellPlot(brain_small,
             spat_enr_names = 'cluster_metagene',
             cell_annotation_values = netw_ranks$clusters,
             point_size = 1.5, cow_n_col = 3)




## 12. spatial HMRF domains ####
hmrf_folder = paste0(temp_dir,'/','11_HMRF/')
if(!file.exists(hmrf_folder)) dir.create(hmrf_folder, recursive = T)

# perform hmrf
my_spatial_genes = km_spatialgenes[1:100]$genes
HMRF_spatial_genes = doHMRF(gobject = brain_small,
                            expression_values = 'scaled',
                            spatial_genes = my_spatial_genes,
                            spatial_network_name = 'Delaunay_network',
                            k = 8,
                            betas = c(28,2,2),
                            output_folder = paste0(hmrf_folder, '/', 'Spatial_genes_brain/SG_top100_k8_scaled'))

# check and select hmrf
for(i in seq(28, 30, by = 2)) {
  viewHMRFresults2D(gobject = brain_small,
                    HMRFoutput = HMRF_spatial_genes,
                    k = 8, betas_to_view = i,
                    point_size = 2)
}

brain_small = addHMRF(gobject = brain_small,
                      HMRFoutput = HMRF_spatial_genes,
                      k = 8, betas_to_add = c(28),
                      hmrf_name = 'HMRF')

giotto_colors = getDistinctColors(8)
names(giotto_colors) = 1:8
spatPlot(gobject = brain_small, cell_color = 'HMRF_k8_b.28',
         point_size = 3, coord_fix_ratio = 1, cell_color_code = giotto_colors)


