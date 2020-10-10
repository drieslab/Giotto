

library(Giotto)

# createGiottoInstructions(python_path = '/your/path')

## Giotto 0.3.5 ##
## mini-test seqFish Giotto 0.3.5 ##

#temp_dir = '/your/path/'
temp_dir = getwd()
temp_dir = '~/Temp/'

## 1. giotto object ####
expr_path = system.file("extdata", "seqfish_field_expr.txt", package = 'Giotto')
loc_path = system.file("extdata", "seqfish_field_locs.txt", package = 'Giotto')

# default
VC_small <- createGiottoObject(raw_exprs = expr_path, spatial_locs = loc_path)
# with instructions (e.g. specific python path)
# myinstructions = createGiottoInstructions(python_path = '/your/path')
# VC_small <- createGiottoObject(raw_exprs = expr_path, spatial_locs = loc_path, instructions = myinstructions)
showGiottoInstructions(VC_small)

## 2. processing steps ####
VC_small <- filterGiotto(gobject = VC_small, expression_threshold = 0.5, gene_det_in_min_cells = 20, min_det_genes_per_cell = 0)
VC_small <- normalizeGiotto(gobject = VC_small, scalefactor = 6000, verbose = T)
VC_small <- addStatistics(gobject = VC_small)
VC_small <- adjustGiottoMatrix(gobject = VC_small, expression_values = c('normalized'), covariate_columns = c('nr_genes', 'total_expr'))

## 3. dimension reduction ####
VC_small <- calculateHVG(gobject = VC_small)
VC_small <- runPCA(gobject = VC_small, center = T)
?runPCA
screePlot(VC_small, ncp = 20)
jackstrawPlot(VC_small, ncp = 20)
plotPCA(VC_small)

VC_small <- runUMAP(VC_small, dimensions_to_use = 1:5, n_threads = 2)
plotUMAP(gobject = VC_small)
VC_small <- runtSNE(VC_small, dimensions_to_use = 1:5)
plotTSNE(gobject = VC_small)

## 4. clustering ####
VC_small <- createNearestNetwork(gobject = VC_small, dimensions_to_use = 1:5, k = 5)
VC_small <- doLeidenCluster(gobject = VC_small, resolution = 0.4, n_iterations = 1000)
plotUMAP(gobject = VC_small, cell_color = 'leiden_clus', show_NN_network = T, point_size = 2.5)
spatDimPlot(gobject = VC_small, cell_color = 'leiden_clus', spat_point_shape = 'voronoi')
showClusterHeatmap(gobject = VC_small, cluster_column = 'leiden_clus')
showClusterDendrogram(VC_small, h = 0.5, rotate = T, cluster_column = 'leiden_clus')

## 5. differential expression ####
gini_markers = findMarkers_one_vs_all(gobject = VC_small,
                                                  method = 'gini',
                                                  expression_values = 'normalized',
                                                  cluster_column = 'leiden_clus',
                                                  min_genes = 20,
                                                  min_expr_gini_score = 0.5,
                                                  min_det_gini_score = 0.5)
topgenes_gini = gini_markers[, head(.SD, 2), by = 'cluster']
violinPlot(VC_small, genes = topgenes_gini$genes, cluster_column = 'leiden_clus')

topgenes_gini2 = gini_markers[, head(.SD, 6), by = 'cluster']
plotMetaDataHeatmap(VC_small, selected_genes = topgenes_gini2$genes,
                    metadata_cols = c('leiden_clus'))


## 6. cell type annotation ####
clusters_cell_types = c('cell A', 'cell B', 'cell C', 'cell D',
                        'cell E', 'cell F', 'cell G')
names(clusters_cell_types) = 1:7
VC_small = annotateGiotto(gobject = VC_small, annotation_vector = clusters_cell_types,
                         cluster_column = 'leiden_clus', name = 'cell_types')
spatDimPlot(gobject = VC_small, cell_color = 'cell_types', spat_point_size = 3, dim_point_size = 3)


## 7. spatial grid ####
VC_small <- createSpatialGrid(gobject = VC_small,
                              sdimx_stepsize = 300,
                              sdimy_stepsize = 300,
                              minimum_padding = 50)
showGrids(VC_small)

spatPlot(gobject = VC_small, show_grid = T, point_size = 1.5)


## 8. spatial network ####
plotStatDelaunayNetwork(gobject = VC_small, maximum_distance = 400)
VC_small = createSpatialNetwork(gobject = VC_small, minimum_k = 2, maximum_distance_delaunay = 400)
VC_small = createSpatialNetwork(gobject = VC_small, minimum_k = 2, method = 'kNN', k = 10)
showNetworks(VC_small)

spatPlot(gobject = VC_small, show_network = T,
         network_color = 'blue', spatial_network_name = 'Delaunay_network',
         point_size = 2.5, cell_color = 'leiden_clus')

spatPlot(gobject = VC_small, show_network = T,
         network_color = 'blue', spatial_network_name = 'kNN_network',
         point_size = 2.5, cell_color = 'leiden_clus')


## 9. spatial genes ####

km_spatialgenes = binSpect(VC_small)

spatGenePlot(VC_small, expression_values = 'scaled', genes = km_spatialgenes[1:6]$genes,
             point_shape = 'border', point_border_stroke = 0.1,
             show_network = F, network_color = 'lightgrey', point_size = 2.5,
             cow_n_col = 2)


## combined approach ##
km_spatialgenes[p.value <= 0.05]

km_spatialgenes1 = km_spatialgenes
km_spatialgenes1[, k := 1]
km_spatialgenes2 = km_spatialgenes
km_spatialgenes2[, k := 2]
km_spatialgenes3 = km_spatialgenes
km_spatialgenes3[, k := 3]

km_spatialgenes_tot = data.table::rbindlist(l = list(km_spatialgenes1, km_spatialgenes2, km_spatialgenes3))

comb_km = km_spatialgenes_tot[, sum(log(p.value)), by = genes]
comb_km[, V1 := V1*-2]
comb_km[, p.val := pchisq(q = V1, df = 2*3, log.p = F, lower.tail = F)]

pchisq(q = 1, df = 2*3, log.p = T)
pchisq(q = 6, df = 2*3, log.p = T)


km_spatialgenes[p.value < 0.01 ]
comb_km[p.val < 0.01]
comb_km

## simulations ##

# create smaller object

set.seed(1234)
sample_genes = sample(VC_small@gene_ID, 100)

VC_small_subset = subsetGiotto(VC_small, gene_ids = sample_genes)
VC_small_subset <- filterGiotto(gobject = VC_small_subset, expression_threshold = 0.5, gene_det_in_min_cells = 20, min_det_genes_per_cell = 0)
VC_small_subset <- normalizeGiotto(gobject = VC_small_subset, scalefactor = 6000, verbose = T)




# pattern 1: bottom right stripe
pattern = VC_small_subset@spatial_locs[sdimx > 1500 & sdimy < -500]
pattern_ids = pattern$cell_ID

selected_genes = sample_genes[1:4]
my_dir = '/Users/rubendries/Dropbox (Personal)/Projects/GC_lab/Ruben_Dries/190225_spatial_package/Results/Paper_revisions/NatMethod_revisions/Revision_1/Spatial_sim_tests/right_lower_patch/'

right_patch_pattern = runPatternSimulation(gobject = VC_small_subset,
                                           pattern_name = 'right_patch',
                                           pattern_cell_ids = pattern_ids,
                                           gene_names = selected_genes,
                                           spatial_probs = c(0.5, 0.8, 0.9, 0.95, 0.99, 1),
                                           reps = 6,
                                           save_plot = T,
                                           max_col = 2,
                                           save_dir = my_dir)


# pattern 2: central patch
pattern = VC_small_subset@spatial_locs[sdimx > 750 & sdimx < 1250 & sdimy > -1250 & sdimy < -750]
pattern_ids = pattern$cell_ID

selected_genes = sample_genes[1:4]
my_dir = '/Users/rubendries/Dropbox (Personal)/Projects/GC_lab/Ruben_Dries/190225_spatial_package/Results/Paper_revisions/NatMethod_revisions/Revision_1/Spatial_sim_tests/center_patch/'

center_patch_pattern = runPatternSimulation(gobject = VC_small_subset,
                                           pattern_name = 'center_patch',
                                           pattern_cell_ids = pattern_ids,
                                           gene_names = selected_genes,
                                           spatial_probs = c(0.5, 0.8, 0.9, 0.95, 0.99, 1),
                                           reps = 6,
                                           save_plot = T,
                                           max_col = 2,
                                           save_dir = my_dir)

# pattern 3: stripe
pattern = VC_small_subset@spatial_locs[sdimx > 800 & sdimx < 1200]
pattern_ids = pattern$cell_ID

selected_genes = sample_genes[1:4]
my_dir = '/Users/rubendries/Dropbox (Personal)/Projects/GC_lab/Ruben_Dries/190225_spatial_package/Results/Paper_revisions/NatMethod_revisions/Revision_1/Spatial_sim_tests/stripe/'

center_patch_pattern = runPatternSimulation(gobject = VC_small_subset,
                                            pattern_name = 'stripe',
                                            pattern_cell_ids = pattern_ids,
                                            gene_names = selected_genes,
                                            spatial_probs = c(0.5, 0.8, 0.9, 0.95, 0.99, 1),
                                            reps = 6,
                                            save_plot = T,
                                            max_col = 2,
                                            save_dir = my_dir)
















rank_spatialgenes = binSpect(VC_small, bin_method = 'rank')
spatGenePlot(VC_small, expression_values = 'scaled', genes = rank_spatialgenes[1:6]$genes,
             point_shape = 'border', point_border_stroke = 0.1,
             show_network = F, network_color = 'lightgrey', point_size = 2.5,
             cow_n_col = 2)

silh_spatialgenes = silhouetteRank(gobject = VC_small) # TODO: suppress print output
spatGenePlot(VC_small, expression_values = 'scaled', genes = silh_spatialgenes[1:6]$genes,
             point_shape = 'border', point_border_stroke = 0.1,
             show_network = F, network_color = 'lightgrey', point_size = 2.5,
             cow_n_col = 2)



## 10. spatial co-expression patterns ####
ext_spatial_genes = km_spatialgenes[1:500]$genes
spat_cor_netw_DT = detectSpatialCorGenes(VC_small,
                                         method = 'network', spatial_network_name = 'Delaunay_network',
                                         subset_genes = ext_spatial_genes)

spat_cor_netw_DT = clusterSpatialCorGenes(spat_cor_netw_DT, name = 'spat_netw_clus', k = 8)
heatmSpatialCorGenes(VC_small, spatCorObject = spat_cor_netw_DT, use_clus_name = 'spat_netw_clus')

netw_ranks = rankSpatialCorGroups(VC_small, spatCorObject = spat_cor_netw_DT, use_clus_name = 'spat_netw_clus')
top_netw_spat_cluster = showSpatialCorGenes(spat_cor_netw_DT, use_clus_name = 'spat_netw_clus',
                                            selected_clusters = 6, show_top_genes = 1)

cluster_genes_DT = showSpatialCorGenes(spat_cor_netw_DT, use_clus_name = 'spat_netw_clus', show_top_genes = 1)
cluster_genes = cluster_genes_DT$clus; names(cluster_genes) = cluster_genes_DT$gene_ID
VC_small = createMetagenes(VC_small, gene_clusters = cluster_genes, name = 'cluster_metagene')
spatCellPlot(VC_small,
             spat_enr_names = 'cluster_metagene',
             cell_annotation_values = netw_ranks$clusters,
             point_size = 1.5, cow_n_col = 3)



## 11. spatial HMRF domains ####
hmrf_folder = paste0(temp_dir,'/','11_HMRF/')
if(!file.exists(hmrf_folder)) dir.create(hmrf_folder, recursive = T)

# perform hmrf
my_spatial_genes = km_spatialgenes[1:100]$genes
HMRF_spatial_genes = doHMRF(gobject = VC_small,
                            expression_values = 'scaled',
                            spatial_genes = my_spatial_genes,
                            spatial_network_name = 'Delaunay_network',
                            k = 9,
                            betas = c(28,2,2),
                            output_folder = paste0(hmrf_folder, '/', 'Spatial_genes/SG_top100_k9_scaled'))

# check and select hmrf
for(i in seq(28, 30, by = 2)) {
  viewHMRFresults2D(gobject = VC_small,
                    HMRFoutput = HMRF_spatial_genes,
                    k = 9, betas_to_view = i,
                    point_size = 2)
}

VC_small = addHMRF(gobject = VC_small,
                  HMRFoutput = HMRF_spatial_genes,
                  k = 9, betas_to_add = c(28),
                  hmrf_name = 'HMRF')

giotto_colors = Giotto:::getDistinctColors(9)
names(giotto_colors) = 1:9
spatPlot(gobject = VC_small, cell_color = 'HMRF_k9_b.28',
         point_size = 3, coord_fix_ratio = 1, cell_color_code = giotto_colors)


## 12. cell neighborhood: cell-type/cell-type interactions ####
set.seed(seed = 2841)
cell_proximities = cellProximityEnrichment(gobject = VC_small,
                                           cluster_column = 'cell_types',
                                           spatial_network_name = 'Delaunay_network',
                                           adjust_method = 'fdr',
                                           number_of_simulations = 1000)

# barplot
cellProximityBarplot(gobject = VC_small, CPscore = cell_proximities,
                     min_orig_ints = 1, min_sim_ints = 1, p_val = 0.25)

## heatmap
cellProximityHeatmap(gobject = VC_small, CPscore = cell_proximities, order_cell_types = T, scale = T,
                     color_breaks = c(-1.5, 0, 1.5), color_names = c('blue', 'white', 'red'))

# network
cellProximityNetwork(gobject = VC_small, CPscore = cell_proximities, remove_self_edges = T, only_show_enrichment_edges = T)

# network with self-edges
cellProximityNetwork(gobject = VC_small, CPscore = cell_proximities,
                     remove_self_edges = F, self_loop_strength = 0.3,
                     only_show_enrichment_edges = F,
                     rescale_edge_weights = T,
                     node_size = 8,
                     edge_weight_range_depletion = c(1, 2),
                     edge_weight_range_enrichment = c(2,5))


## visualization of specific cell types
# Option 1
spec_interaction = "cell D--cell F"
cellProximitySpatPlot2D(gobject = VC_small,
                        interaction_name = spec_interaction,
                        show_network = T,
                        cluster_column = 'cell_types',
                        cell_color = 'cell_types',
                        cell_color_code = c('cell D' = 'lightblue', 'cell F' = 'red'),
                        point_size_select = 4, point_size_other = 2)

# Option 2: create additional metadata
VC_small = addCellIntMetadata(VC_small,
                             spatial_network = 'Delaunay_network',
                             cluster_column = 'cell_types',
                             cell_interaction = spec_interaction,
                             name = 'D_F_interactions')
spatPlot(VC_small, cell_color = 'D_F_interactions', legend_symbol_size = 3,
         select_cell_groups =  c('other_cell D', 'other_cell F', 'select_cell D', 'select_cell F'))




### 13. cell neighborhood: interaction changed genes ####

## select top 25th highest expressing genes
gene_metadata = fDataDT(VC_small)
plot(gene_metadata$nr_cells, gene_metadata$mean_expr)
plot(gene_metadata$nr_cells, gene_metadata$mean_expr_det)

quantile(gene_metadata$mean_expr_det)
high_expressed_genes = gene_metadata[mean_expr_det > 4]$gene_ID

## identify genes that are associated with proximity to other cell types
CPGscoresHighGenes =  findCPG(gobject = VC_small,
                              selected_genes = high_expressed_genes,
                              spatial_network_name = 'Delaunay_network',
                              cluster_column = 'cell_types',
                              diff_test = 'permutation',
                              adjust_method = 'fdr',
                              nr_permutations = 500,
                              do_parallel = T, cores = 2)

## visualize all genes
plotCellProximityGenes(VC_small, cpgObject = CPGscoresHighGenes, method = 'dotplot')

## filter genes
CPGscoresFilt = filterCPG(CPGscoresHighGenes, min_cells = 2, min_int_cells = 2, min_fdr = 0.1,
                          min_spat_diff = 0.1, min_log2_fc = 0.1, min_zscore = 1)

## visualize subset of interaction changed genes (ICGs)
ICG_genes = c('Cpne2', 'Scg3', 'Cmtm3', 'Cplx1', 'Lingo1')
ICG_genes_types = c('cell E', 'cell D', 'cell D', 'cell G', 'cell E')
names(ICG_genes) = ICG_genes_types

plotICG(gobject = VC_small,
        cpgObject = CPGscoresHighGenes,
        source_type = 'cell A',
        source_markers = c('Csf1r', 'Laptm5'),
        ICG_genes = ICG_genes)




##### 14. cell neighborhood:  ligand-receptor cell-cell communication ####

LR_data = data.table::fread(system.file("extdata", "mouse_ligand_receptors.txt", package = 'Giotto'))

LR_data[, ligand_det := ifelse(mouseLigand %in% VC_small@gene_ID, T, F)]
LR_data[, receptor_det := ifelse(mouseReceptor %in% VC_small@gene_ID, T, F)]
LR_data_det = LR_data[ligand_det == T & receptor_det == T]
select_ligands = LR_data_det$mouseLigand
select_receptors = LR_data_det$mouseReceptor


## get statistical significance of gene pair expression changes based on expression ##
expr_only_scores = exprCellCellcom(gobject = VC_small,
                                   cluster_column = 'cell_types',
                                   random_iter = 500,
                                   gene_set_1 = select_ligands,
                                   gene_set_2 = select_receptors)

## get statistical significance of gene pair expression changes upon cell-cell interaction
spatial_all_scores = spatCellCellcom(VC_small,
                                     spatial_network_name = 'Delaunay_network',
                                     cluster_column = 'cell_types',
                                     random_iter = 500,
                                     gene_set_1 = select_ligands,
                                     gene_set_2 = select_receptors,
                                     adjust_method = 'fdr',
                                     do_parallel = T,
                                     cores = 4,
                                     verbose = 'none')


## * plot communication scores ####

## select top LR ##
selected_spat = spatial_all_scores[p.adj <= 0.5 & abs(log2fc) > 0.1 & lig_nr >= 2 & rec_nr >= 2]
data.table::setorder(selected_spat, -PI)

top_LR_ints = unique(selected_spat[order(-abs(PI))]$LR_comb)[1:33]
top_LR_cell_ints = unique(selected_spat[order(-abs(PI))]$LR_cell_comb)[1:33]

plotCCcomHeatmap(gobject = VC_small,
                 comScores = spatial_all_scores,
                 selected_LR = top_LR_ints,
                 selected_cell_LR = top_LR_cell_ints,
                 show = 'LR_expr')


plotCCcomDotplot(gobject = VC_small,
                 comScores = spatial_all_scores,
                 selected_LR = top_LR_ints,
                 selected_cell_LR = top_LR_cell_ints,
                 cluster_on = 'PI')



## * spatial vs rank ####
comb_comm = combCCcom(spatialCC = spatial_all_scores,
                      exprCC = expr_only_scores)

# top differential activity levels for ligand receptor pairs
plotRankSpatvsExpr(gobject = VC_small,
                   comb_comm,
                   expr_rnk_column = 'exprPI_rnk',
                   spat_rnk_column = 'spatPI_rnk',
                   midpoint = 10)

## * recovery ####
## predict maximum differential activity
plotRecovery(gobject = VC_small,
             comb_comm,
             expr_rnk_column = 'exprPI_rnk',
             spat_rnk_column = 'spatPI_rnk',
             ground_truth = 'spatial')


### 15. export Giotto Analyzer to Viewer ####
viewer_folder = paste0(temp_dir, '/', 'Mouse_cortex_viewer')

# select annotations, reductions and expression values to view in Giotto Viewer
exportGiottoViewer(gobject = VC_small, output_directory = viewer_folder,
                   factor_annotations = c('cell_types',
                                          'leiden_clus',
                                          'HMRF_k9_b.28'),
                   numeric_annotations = 'total_expr',
                   dim_reductions = c('umap'),
                   dim_reduction_names = c('umap'),
                   expression_values = 'scaled',
                   expression_rounding = 3,
                   overwrite_dir = T)



