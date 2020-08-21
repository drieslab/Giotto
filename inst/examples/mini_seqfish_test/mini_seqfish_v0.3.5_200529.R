

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
VC_small <- runPCA(gobject = VC_small)
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
km_spatialgenes = binSpect(VC_small, spatial_network_name = 'kNN_network')
spatGenePlot(VC_small, expression_values = 'scaled', genes = km_spatialgenes[1:6]$genes,
             point_shape = 'border', point_border_stroke = 0.1,
             show_network = F, network_color = 'lightgrey', point_size = 2.5,
             cow_n_col = 2)

km_spatialgenes[p.value <= 0.05]






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


simulateOneGenePatternGiottoObject = function(gobject,
                                       pattern_name = 'pattern',
                                       pattern_cell_ids = NULL,
                                       gene_name = NULL,
                                       spatial_prob = 0.95,
                                       gradient_direction = NULL,
                                       show_pattern = TRUE,
                                       pattern_colors = c('in' = 'green', 'out' = 'red'),
                                       ...) {

  if(is.null(pattern_cell_ids)) {
    stop('pattern_cell_ids can not be NULL \n')
  }

  ## create and add annotation for pattern
  cell_meta = pDataDT(gobject)
  cell_meta[, (pattern_name) := ifelse(cell_ID %in% pattern_cell_ids, 'in', 'out')]

  newgobject = addCellMetadata(gobject,
                               new_metadata = cell_meta[,c('cell_ID', pattern_name), with = F],
                               by_column = T,
                               column_cell_ID = 'cell_ID')

  # show pattern
  if(show_pattern == TRUE) {
    spatPlot2D(gobject = newgobject, save_plot = F, cell_color_code = pattern_colors,
               point_size = 2, cell_color = pattern_name)
  }


  ## merge cell metadata and cell coordinate data
  cell_meta = pDataDT(newgobject)
  cell_coord = newgobject@spatial_locs
  cell_meta = data.table::merge.data.table(cell_meta, cell_coord, by = 'cell_ID')

  ## get number of cells within pattern
  cell_number = nrow(cell_meta[get(pattern_name) == 'in'])


  ## normalized expression
  expr_data = newgobject@norm_expr
  result_list = list()

  ## raw expression
  raw_expr_data = newgobject@raw_exprs
  raw_result_list = list()


  ## create the spatial expression pattern for the specified gene
  # 1. rank all gene values from the cells from high to low
  # 2. move the highest expressing values to the spatial pattern using a probability
  #     - 0.5 is the control = random
  #     - 1 is perfection: all the highest values go to the pattern
  #     - 0.5 to 1 is decreasing noise levels

  if(is.null(gene_name)) stop('a gene name needs to be provided')



  # rank genes
  gene_vector = expr_data[rownames(expr_data) == gene_name, ]
  sort_expr_gene = sort(gene_vector, decreasing = T)

  # number of cells in and out the pattern
  total_cell_number = length(sort_expr_gene)
  remaining_cell_number = total_cell_number - cell_number

  # calculate outside probability
  outside_prob = 1 - spatial_prob
  prob_vector = c(rep(spatial_prob, cell_number), rep(outside_prob, remaining_cell_number))

  # first get the 'in' pattern sample values randomly
  sample_values = sample(sort_expr_gene, replace = F, size = cell_number, prob = prob_vector)

  # then take the remaining 'out' pattern values randomly
  remain_values = sort_expr_gene[!names(sort_expr_gene) %in% names(sample_values)]
  remain_values = sample(remain_values, size = length(remain_values))



  ## A. within pattern ##
  # ------------------- #
  in_cell_meta = cell_meta[get(pattern_name) == 'in']

  # if gradient is wanted
  # does not work with 0.5!! is not random!!
  if(!is.null(gradient_direction)) {
    # sort in_ids according to x, y or  xy coordinates to create gradient
    in_cell_meta[, sdimx_y := abs(sdimx)+ abs(sdimy)]
    # order according to gradient direction
    in_cell_meta = in_cell_meta[order(get(gradient_direction))]
  }
  in_ids = in_cell_meta$cell_ID

  # preparation for raw matrix
  sample_values_id_vector = names(sample_values)
  names(sample_values_id_vector) = in_ids


  ## B. outside pattern ##
  # -------------------- #
  out_ids = cell_meta[get(pattern_name) == 'out']$cell_ID

  # preparation for raw matrix
  remain_values_id_vector = names(remain_values)
  names(remain_values_id_vector) = out_ids




  ## raw matrix
  # swap the cell ids #
  raw_gene_vector = raw_expr_data[rownames(raw_expr_data) == gene_name,]

  raw_new_sample_vector = raw_gene_vector[sample_values_id_vector]
  names(raw_new_sample_vector) = names(sample_values_id_vector)

  raw_new_remain_vector = raw_gene_vector[remain_values_id_vector]
  names(raw_new_remain_vector) = names(remain_values_id_vector)

  new_sim_raw_values = c(raw_new_sample_vector, raw_new_remain_vector)
  new_sim_raw_values = new_sim_raw_values[names(raw_gene_vector)]

  # change the original matrices
  raw_expr_data[rownames(raw_expr_data) == gene_name,] = new_sim_raw_values
  newgobject@raw_exprs = raw_expr_data

  # recalculate normalized values
  newgobject <- normalizeGiotto(gobject = newgobject, ...)
  newgobject <- addStatistics(gobject = newgobject)

  return(newgobject)

}



# original
spatGenePlot2D(VC_small_subset, expression_values = 'norm', genes = 'Abca4',
             point_shape = 'border', point_border_stroke = 0.1,
             show_network = F, network_color = 'lightgrey', point_size = 2.5,
             cow_n_col = 1, show_plot = F,
             save_plot = T, save_param = list(save_dir = '~', save_folder = 'pattern_name', save_name = 'plot_prob0.9_rep1',
                                              base_width = 9, base_height = 7, units = 'cm'))

showSaveParameters()

right_gene_patch = simulateOneGenePatternGiottoObject(VC_small_subset,
                                                      pattern_name = 'right_patch',
                                                      pattern_cell_ids = pattern_ids,
                                                      gene_name = 'Abca4',
                                                      spatial_prob = 0.99,
                                                      scalefactor = 6000, verbose = T)



spatGenePlot(right_gene_patch, expression_values = 'norm', genes = 'Abca4',
             point_shape = 'border', point_border_stroke = 0.1,
             show_network = F, network_color = 'lightgrey', point_size = 2.5,
             cow_n_col = 1, show_plot = F)




run_spatial_sim_tests_one_rep = function(gobject,
                                         pattern_name = 'pattern',
                                         pattern_cell_ids = NULL,
                                         gene_name = NULL,
                                         spatial_prob = 0.95,
                                         show_pattern = FALSE,

                                         # binspect kmeans
                                         spatial_network_name = 'kNN_network',
                                         binSpect_km_param = list(nstart = 3,
                                                                  iter_max = 10,
                                                                  expression_values = 'normalized',
                                                                  get_av_expr = FALSE,
                                                                  get_high_expr = FALSE),
                                         # binspect rank
                                         binSpect_rnk_param = list(percentage_rank = 30,
                                                                  expression_values = 'normalized',
                                                                  get_av_expr = FALSE,
                                                                  get_high_expr = FALSE),

                                         # spatialDE
                                         spatialDE_param = list(expression_values = 'raw',
                                                                sig_alpha = 0.5,
                                                                unsig_alpha = 0.5),

                                         # spark
                                         spark_param = list(values_type = 'raw',
                                                            percentage = 0.1,
                                                            min_count = 10,
                                                            num_core = 5),


                                         save_plot = F,
                                         save_dir = '~',
                                         save_name = 'plot',

                                         ...) {


  simulate_patch = simulateOneGenePatternGiottoObject(gobject,
                                                      pattern_name = pattern_name,
                                                      pattern_cell_ids = pattern_cell_ids,
                                                      gene_name = gene_name,
                                                      spatial_prob = spatial_prob,
                                                      gradient_direction = NULL,
                                                      show_pattern = show_pattern,
                                                      ...)


  if(save_plot == TRUE) {

    spatGenePlot2D(simulate_patch, expression_values = 'norm', genes = gene_name,
                   point_shape = 'border', point_border_stroke = 0.1, point_size = 2.5,
                   cow_n_col = 1, show_plot = F,
                   save_plot = T, save_param = list(save_dir = save_dir, save_folder = pattern_name,
                                                    save_name = save_name,
                                                    base_width = 9, base_height = 7, units = 'cm'))

  }



  ## spatial gene detection methods ##

  ## binspect kmeans
  ## -------------- ##
  start = proc.time()
  km_spatialgenes_sim = do.call('binSpect', c(gobject =  simulate_patch,
                                              bin_method = 'kmeans',
                                              spatial_network_name = spatial_network_name,
                                              binSpect_km_param))

  km_spatialgenes_sim[, adj.p.value := p.adjust(p.value, method = 'fdr')]
  binspec_km_result = km_spatialgenes_sim[genes == gene_name]
  binspec_km_time = proc.time() - start

  binspec_km_result[, prob := spatial_prob]
  binspec_km_result[, time := binspec_km_time[['elapsed']] ]

  binspec_km_result = binspec_km_result[,.(genes, p.value, adj.p.value, prob, time)]
  colnames(binspec_km_result) = c('genes', 'p.value', 'adj.p.value', 'prob', 'time')
  binspec_km_result[, method := 'binspec_km']


  ## binspect rank
  ## ---------- ##
  start = proc.time()
  rnk_spatialgenes_sim = do.call('binSpect', c(gobject =  simulate_patch,
                                              bin_method = 'rank',
                                              spatial_network_name = spatial_network_name,
                                              binSpect_rnk_param))

  rnk_spatialgenes_sim[, adj.p.value := p.adjust(p.value, method = 'fdr')]
  binspec_rnk_result = rnk_spatialgenes_sim[genes == gene_name]
  binspec_rnk_time = proc.time() - start

  binspec_rnk_result[, prob := spatial_prob]
  binspec_rnk_result[, time := binspec_rnk_time[['elapsed']] ]

  binspec_rnk_result = binspec_rnk_result[,.(genes, p.value, adj.p.value, prob, time)]
  colnames(binspec_rnk_result) = c('genes', 'p.value', 'adj.p.value', 'prob', 'time')
  binspec_rnk_result[, method := 'binspec_rnk']



  ## spatialDE
  ## -------- ##
  start = proc.time()
  new_raw_sim_matrix = simulate_patch@raw_exprs
  sd_cells = apply(new_raw_sim_matrix, 2, sd)
  sd_non_zero_cells = names(sd_cells[sd_cells != 0])
  simulate_patch_fix = subsetGiotto(simulate_patch, cell_ids = sd_non_zero_cells)

  spatialDE_spatialgenes_sim = do.call('spatialDE', c(gobject =  simulate_patch_fix,
                                                      show_plot = FALSE,
                                                      return_plot = FALSE,
                                                      save_plot = FALSE,
                                                      spatialDE_param))

  spatialDE_spatialgenes_sim_res = spatialDE_spatialgenes_sim$results$results
  if(is.null(spatialDE_spatialgenes_sim_res)) spatialDE_spatialgenes_sim_res = spatialDE_spatialgenes_sim$results
  spatialDE_spatialgenes_sim_res = data.table::as.data.table(spatialDE_spatialgenes_sim_res)
  data.table::setorder(spatialDE_spatialgenes_sim_res, qval, pval)
  spatialDE_result = spatialDE_spatialgenes_sim_res[g == gene_name]

  spatialDE_time = proc.time() - start

  spatialDE_result[, prob := spatial_prob]
  spatialDE_result[, time := spatialDE_time[['elapsed']] ]

  spatialDE_result = spatialDE_result[,.(g, pval, qval, prob, time)]
  colnames(spatialDE_result) = c('genes', 'p.value', 'adj.p.value', 'prob', 'time')
  spatialDE_result[, method := 'spatialDE']



  ## spark
  start = proc.time()
  spark_spatialgenes_sim = do.call('spark', c(gobject =  simulate_patch,
                                              return_object = 'data.table',
                                              spark_param))

  spark_result = spark_spatialgenes_sim[genes == gene_name]
  spark_time = proc.time() - start

  spark_result[, prob := spatial_prob]
  spark_result[, time := spark_time[['elapsed']] ]

  spark_result = spark_result[,.(genes, combined_pvalue, adjusted_pvalue, prob, time)]
  colnames(spark_result) = c('genes', 'p.value', 'adj.p.value', 'prob', 'time')
  spark_result[, method := 'spark']

  results = do.call('rbind', list(binspec_km_result, binspec_rnk_result, spatialDE_result, spark_result))


  return(results)

}



test = run_spatial_sim_tests_one_rep(VC_small_subset,
                                     pattern_name = 'right_patch',
                                     pattern_cell_ids = pattern_ids,
                                     gene_name = 'Abca4',
                                     spatial_prob = 0.99,
                                     scalefactor = 6000,
                                     verbose = F,

                                     # binspect kmeans
                                     spatial_network_name = 'kNN_network',
                                     binSpect_km_param = list(nstart = 3,
                                                              iter_max = 10,
                                                              expression_values = 'normalized',
                                                              get_av_expr = FALSE,
                                                              get_high_expr = FALSE),
                                     # binspect rank
                                     binSpect_rnk_param = list(percentage_rank = 30,
                                                               expression_values = 'normalized',
                                                               get_av_expr = FALSE,
                                                               get_high_expr = FALSE),

                                     # spatialDE
                                     spatialDE_param = list(expression_values = 'raw',
                                                            sig_alpha = 0.5,
                                                            unsig_alpha = 0.5),

                                     # spark
                                     spark_param = list(values_type = 'raw',
                                                        percentage = 0.1,
                                                        min_count = 10,
                                                        num_core = 5))



run_spatial_sim_tests_multi = function(gobject,
                                       pattern_name = 'pattern',
                                       pattern_cell_ids = NULL,
                                       gene_name = NULL,
                                       spatial_probs = c(0.5, 1),
                                       reps = 2,

                                       # binspect kmeans
                                       spatial_network_name = 'kNN_network',
                                       binSpect_km_param = list(nstart = 3,
                                                                iter_max = 10,
                                                                expression_values = 'normalized',
                                                                get_av_expr = FALSE,
                                                                get_high_expr = FALSE),
                                       # binspect rank
                                       binSpect_rnk_param = list(percentage_rank = 30,
                                                                 expression_values = 'normalized',
                                                                 get_av_expr = FALSE,
                                                                 get_high_expr = FALSE),

                                       # spatialDE
                                       spatialDE_param = list(expression_values = 'raw',
                                                              sig_alpha = 0.5,
                                                              unsig_alpha = 0.5),

                                       # spark
                                       spark_param = list(values_type = 'raw',
                                                          percentage = 0.1,
                                                          min_count = 10,
                                                          num_core = 5),

                                       save_plot = F,
                                       save_dir = '~',
                                       ...) {

  prob_list = list()
  for(prob_ind in 1:length(spatial_probs)) {

    prob_i = spatial_probs[prob_ind]

    cat('\n \n start with ', prob_i, '\n \n')

    rep_list = list()
    for(rep_i in 1:reps) {


      cat('\n \n repetitiion = ', rep_i, '\n \n')


      plot_name = paste0('plot_',gene_name,'_prob', prob_i, '_rep', rep_i)


      rep_res = run_spatial_sim_tests_one_rep(gobject,
                                           pattern_name = pattern_name,
                                           pattern_cell_ids = pattern_cell_ids,
                                           gene_name = gene_name,
                                           spatial_prob = prob_i,

                                           spatial_network_name = spatial_network_name,
                                           binSpect_km_param = binSpect_km_param,
                                           binSpect_rnk_param = binSpect_rnk_param,
                                           spatialDE_param = spatialDE_param,
                                           spark_param = spark_param,

                                           save_plot = save_plot,
                                           save_dir = save_dir,
                                           save_name = plot_name,
                                           ...)

      rep_res[, rep := rep_i]
      rep_list[[rep_i]] = rep_res

    }

    rep_list_res = do.call('rbind', rep_list)
    prob_list[[prob_ind]] = rep_list_res

  }

  final_gene_results = do.call('rbind', prob_list)

  return(final_gene_results)

}


spatial_probs = c(0.5, 0.65, 0.8, 0.9, 0.95, 0.99, 1)
reps = 10

starttime = proc.time()
testmulti = run_spatial_sim_tests_multi(VC_small_subset,
                                     pattern_name = 'right_patch',
                                     pattern_cell_ids = pattern_ids,
                                     gene_name = 'Abca4',
                                     spatial_probs = c(0.5, 0.8, 1),
                                     reps = 3,
                                     spatial_network_name = 'kNN_network',
                                     binSpect_km_param = list(nstart = 3,
                                                              iter_max = 10,
                                                              expression_values = 'normalized',
                                                              get_av_expr = FALSE,
                                                              get_high_expr = FALSE),
                                     # binspect rank
                                     binSpect_rnk_param = list(percentage_rank = 30,
                                                               expression_values = 'normalized',
                                                               get_av_expr = FALSE,
                                                               get_high_expr = FALSE),

                                     # spatialDE
                                     spatialDE_param = list(expression_values = 'raw',
                                                            sig_alpha = 0.5,
                                                            unsig_alpha = 0.5),

                                     # spark
                                     spark_param = list(values_type = 'raw',
                                                        percentage = 0.1,
                                                        min_count = 10,
                                                        num_core = 5),
                                     save_plot = T,
                                     save_dir = '/Users/rubendries/Dropbox (Personal)/Projects/GC_lab/Ruben_Dries/190225_spatial_package/Results/Paper_revisions/NatMethod_revisions/Revision_1/Spatial_sim_tests/',
                                     scalefactor = 6000,
                                     verbose = F)
testmulti[, prob := as.factor(prob)]
testmulti[, method := factor(method, levels = c('binspec_km', 'binspec_rnk', 'spatialDE', 'spark'))]
total_time = proc.time() - starttime







runPatternSimulation = function(gobject,
                                pattern_name = 'pattern',
                                pattern_cell_ids = NULL,
                                gene_names = NULL,
                                save_plot = T,
                                save_dir = '~',
                                max_col = 4,
                                height = 7,
                                width = 7,
                                ...) {


  # plot pattern for first gene (the same for all)
  example_patch = simulateOneGenePatternGiottoObject(gobject,
                                                     pattern_name = pattern_name,
                                                     pattern_cell_ids = pattern_cell_ids,
                                                     gene_name = gene_names[[1]],
                                                     spatial_prob = 1,
                                                     scalefactor = 6000,
                                                     verbose = T)

  spatPlot2D(example_patch, cell_color = pattern_name,
             save_plot = T, save_param = list(save_dir = save_dir, save_folder = 'original', save_name = paste0(pattern_name,'_pattern'),
                                              base_width = 9, base_height = 7, units = 'cm'))


  all_results = list()
  for(gene_ind in 1:length(gene_names)) {

    gene = gene_names[gene_ind]

    # plot original expression
    spatGenePlot2D(gobject, expression_values = 'norm', genes = gene,
                   point_shape = 'border', point_border_stroke = 0.1,
                   show_network = F, network_color = 'lightgrey', point_size = 2.5,
                   cow_n_col = 1, show_plot = F,
                   save_plot = T, save_param = list(save_dir = save_dir, save_folder = 'original', save_name = paste0(gene,'_original'),
                                                    base_width = 9, base_height = 7, units = 'cm'))


    generesults = run_spatial_sim_tests_multi(gobject,
                                            pattern_name = pattern_name,
                                            pattern_cell_ids = pattern_cell_ids,
                                            gene_name = gene,
                                            save_plot = save_plot,
                                            save_dir = save_dir,
                                            ...)
    generesults[, prob := as.factor(prob)]
    generesults[, method := factor(method, levels = c('binspec_km', 'binspec_rnk', 'spatialDE', 'spark'))]

    all_results[[gene_ind]] = generesults

  }

  results = do.call('rbind', all_results)



  ## plot results ##

  # 4 columns max
  nr_rows = max(c(round(length(gene_names)/max_col), 1))

  # p-values
  pl = ggplot()
  pl = pl + geom_boxplot(data = results, aes(x = method, y = adj.p.value, color = prob))
  pl = pl + geom_point(data = results, aes(x = method, y = adj.p.value, color = prob), size = 2, position = position_jitterdodge())
  pl = pl + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))
  pl = pl + facet_wrap(~genes, nrow = nr_rows)
  pl = pl + geom_hline(yintercept = 0.05, color = 'red', linetype = 2)

  pdf(file = paste0(save_dir,'/',pattern_name,'_boxplot_pvalues.pdf'), width = width, height = height)
  print(pl)
  dev.off()


  # -log10 p-values
  pl = ggplot()
  pl = pl + geom_boxplot(data = results, aes(x = method, y = -log10(adj.p.value), color = prob))
  pl = pl + geom_point(data = results, aes(x = method, y = -log10(adj.p.value), color = prob), size = 2, position = position_jitterdodge())
  pl = pl + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))
  pl = pl + facet_wrap(~genes, nrow = nr_rows)

  pdf(file = paste0(save_dir,'/',pattern_name,'_boxplot_log10pvalues.pdf'), width = width, height = height)
  print(pl)
  dev.off()


  # time
  pl = ggplot()
  pl = pl + geom_boxplot(data = results, aes(x = method, y = time, color = prob))
  pl = pl + geom_point(data = results, aes(x = method, y = time, color = prob), size = 2, position = position_jitterdodge())
  pl = pl + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))

  pdf(file = paste0(save_dir,'/',pattern_name,'_boxplot_time.pdf'), width = width, height = height)
  print(pl)
  dev.off()

  # write results
  data.table::fwrite(x = results, file = paste0(save_dir,'/',pattern_name,'_results.txt'), sep = '\t', quote = F)

  return(results)

}


selected_genes = c('Abca4', 'Kcna6')
my_dir = '/Users/rubendries/Dropbox (Personal)/Projects/GC_lab/Ruben_Dries/190225_spatial_package/Results/Paper_revisions/NatMethod_revisions/Revision_1/Spatial_sim_tests/'

testpattern = runPatternSimulation(gobject = VC_small_subset,
                                   pattern_name = 'right_patch',
                                   pattern_cell_ids = pattern_ids,
                                   gene_names = selected_genes,
                                   save_plot = T,
                                   save_dir = my_dir)










## OLD ANALYSIS ##
right_patch = simulatePatternGiottoObject(VC_small, pattern_name = 'right_patch',
                                          pattern_cell_ids = pattern_ids)




# pattern 2: small patch center
pattern = VC_small@spatial_locs[sdimx < 1250 & sdimx > 750 & sdimy < -750 & sdimy > -1250]
pattern_ids = pattern$cell_ID

center_patch = simulatePatternGiottoObject(VC_small, pattern_name = 'center_patch',
                                          pattern_cell_ids = pattern_ids)






# pattern 3: center stripe
pattern = VC_small@spatial_locs[sdimx < 1200 & sdimx > 800]
pattern_ids = pattern$cell_ID


# example
center_stripe = simulatePatternGiottoObject(VC_small,
                                            spatial_prob = 0.99,
                                            pattern_name = 'center_stripe',
                                            pattern_cell_ids = pattern_ids)

km_spatialgenes_sim = binSpect(center_stripe, spatial_network_name = 'kNN_network')

spatGenePlot(center_stripe, expression_values = 'norm', genes = head(km_spatialgenes_sim$genes, 4),
             point_shape = 'border', point_border_stroke = 0.1,
             show_network = F, network_color = 'lightgrey', point_size = 2.5,
             cow_n_col = 2)

spatGenePlot(center_stripe, expression_values = 'norm', genes = tail(km_spatialgenes_sim$genes, 4),
             point_shape = 'border', point_border_stroke = 0.1,
             show_network = F, network_color = 'lightgrey', point_size = 2.5,
             cow_n_col = 2)


## simulations
test_patterns_binspect = function(gobject,
                                  pattern_ids,
                                  pattern_name,
                                  probs = c(0.5, 1),
                                  repetitions = 3,
                                  show_pattern = F,
                                  pattern_colors = c('in' = 'green', 'out' = 'red'),
                                  spatial_network_name = 'kNN_network',
                                  bin_method = 'kmeans',
                                  nstart = nstart,
                                  iter_max = iter_max,
                                  percentage_rank = percentage_rank,
                                  p_values = c(0.01, 0.05, 0.1)) {



  outer_list = list()
  for(prob_index in 1:length(probs)) {

    prob_i = probs[prob_index]

    inner_list = list()
    for(j in 1:repetitions) {

      sim_object = simulatePatternGiottoObject(gobject,
                                               spatial_prob = prob_i,
                                               pattern_name = pattern_name,
                                               pattern_cell_ids = pattern_ids,
                                               show_pattern = show_pattern,
                                               pattern_colors = pattern_colors)

      sim_results = binSpect(sim_object,
                             spatial_network_name = spatial_network_name,
                             bin_method = bin_method,
                             expression_values = 'norm',
                             nstart = nstart,
                             iter_max = iter_max,
                             percentage_rank = percentage_rank)

      sim_results[, method := bin_method]
      sim_results[, prob := prob_i]
      sim_results[, rep := j]

      inner_list[[j]] = sim_results
    }
    inner_list_res = do.call('rbind', inner_list)

    outer_list[[prob_index]] = inner_list_res

  }
  outer_list_res = do.call('rbind', outer_list)



  ## summarize per p-value
  p_results = list()
  for(p_ind in 1:length(p_values)) {

    p_i = p_values[p_ind]
    outer_list_res[, bin := ifelse(p.value <= p_i, 'yes', 'no')]
    summary_1 = outer_list_res[, table(factor(bin, levels = c('yes', 'no'))), by = .(method, prob, rep)]

    reps = nrow(summary_1)/2
    summary_1[, ident := rep(c('found', 'not_found'), reps)]
    summary_1[, frac_V1 := V1/sum(V1), by = c('method', 'prob', 'rep')]
    summary_1[, frac_V1 := as.numeric(frac_V1)]

    summary_1[, p_val := p_i]
    p_results[[p_ind]] = summary_1

  }

  summary_all = do.call('rbind', p_results)
  summary_all[, p_val := as.factor(p_val)]

  return(list(raw = outer_list_res, summary = summary_all))

}



binkmeans_sim = test_patterns_binspect(gobject = VC_small,pattern_ids,
                                       pattern_name = 'center_stripe',
                                       probs = c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99, 1),
                                       repetitions = 2,
                                       show_pattern = F,
                                       pattern_colors = c('in' = 'green', 'out' = 'red'),
                                       spatial_network_name = 'kNN_network',
                                       bin_method = 'kmeans',
                                       nstart = 3,
                                       iter_max = 10,
                                       percentage_rank = 30,
                                       p_values = c(0.01, 0.05, 0.1))


summary_all = binkmeans_sim$summary

library(ggplot2)

pl = ggplot()
pl = pl + geom_point(data = summary_all, aes(x = as.factor(prob), color = ident, y = frac_V1))
pl = pl + geom_line(data = summary_all, aes(x = as.factor(prob), color = ident, y = frac_V1, group = ident))
pl = pl + facet_grid(~ p_val)
pl = pl + theme_bw()
pl

pl = ggplot()
pl = pl + geom_point(data = summary_all, aes(x = as.factor(prob), color = p_val, y = frac_V1))
pl = pl + geom_line(data = summary_all, aes(x = as.factor(prob), color = p_val, y = frac_V1, group = p_val))
pl = pl + facet_grid(~ ident)
pl = pl + theme_bw()
pl




# spatialDE
libsizes = colSums(center_stripe@raw_exprs)
libsize_needed = median(libsizes[pattern_ids])

new_raw_sim_matrix = center_stripe@raw_exprs
sd_cells = apply(new_raw_sim_matrix, 2, sd)
sd_non_zero_cells = names(sd_cells[sd_cells != 0])
center_stripe_fix = subsetGiotto(center_stripe, cell_ids = sd_non_zero_cells)

spatialDE_spatialgenes_sim = spatialDE(gobject = center_stripe_fix)

spatialDE_spatialgenes_sim_res = spatialDE_spatialgenes_sim$results$results
if(is.null(spatialDE_spatialgenes_sim_res)) spatialDE_spatialgenes_sim_res = spatialDE_spatialgenes_sim$results
spatialDE_spatialgenes_sim_res = data.table::as.data.table(spatialDE_spatialgenes_sim_res)
data.table::setorder(spatialDE_spatialgenes_sim_res, qval, pval)

spatGenePlot(center_stripe, expression_values = 'norm', genes = head(spatialDE_spatialgenes_sim_res$g, 4),
             point_shape = 'border', point_border_stroke = 0.1,
             show_network = F, network_color = 'lightgrey', point_size = 2.5,
             cow_n_col = 2)

spatGenePlot(center_stripe, expression_values = 'norm', genes = tail(spatialDE_spatialgenes_sim_res$g, 4),
             point_shape = 'border', point_border_stroke = 0.1,
             show_network = F, network_color = 'lightgrey', point_size = 2.5,
             cow_n_col = 2)



## simulations
test_patterns_spatialDE = function(gobject,
                                  pattern_ids,
                                  pattern_name,
                                  probs = c(0.5, 1),
                                  repetitions = 3,
                                  show_pattern = F,
                                  pattern_colors = c('in' = 'green', 'out' = 'red'),
                                  p_values = c(0.01, 0.05, 0.1),
                                  ...) {



  outer_list = list()
  for(prob_index in 1:length(probs)) {

    prob_i = probs[prob_index]

    inner_list = list()
    for(j in 1:repetitions) {

      sim_object = simulatePatternGiottoObject(gobject,
                                               spatial_prob = prob_i,
                                               pattern_name = pattern_name,
                                               pattern_cell_ids = pattern_ids,
                                               show_pattern = show_pattern,
                                               pattern_colors = pattern_colors)


      new_raw_sim_matrix = sim_object@raw_exprs
      sd_cells = apply(new_raw_sim_matrix, 2, sd)
      sd_non_zero_cells = names(sd_cells[sd_cells != 0])
      sim_object_fix = subsetGiotto(sim_object, cell_ids = sd_non_zero_cells)

      sim_results = spatialDE(gobject = sim_object_fix, ...)

      sim_results_res = sim_results$results$results
      if(is.null(sim_results_res)) sim_results_res = sim_results$results
      sim_results_res = data.table::as.data.table(sim_results_res)
      data.table::setorder(sim_results_res, qval, pval)

      sim_results_res[, method := 'spatialDE']
      sim_results_res[, prob := prob_i]
      sim_results_res[, rep := j]

      inner_list[[j]] = sim_results_res
    }
    inner_list_res = do.call('rbind', inner_list)

    outer_list[[prob_index]] = inner_list_res

  }
  outer_list_res = do.call('rbind', outer_list)



  ## summarize per p-value
  p_results = list()
  for(p_ind in 1:length(p_values)) {

    p_i = p_values[p_ind]
    outer_list_res[, bin := ifelse(pval <= p_i, 'yes', 'no')]
    summary_1 = outer_list_res[, table(factor(bin, levels = c('yes', 'no'))), by = .(method, prob, rep)]

    reps = nrow(summary_1)/2
    summary_1[, ident := rep(c('found', 'not_found'), reps)]
    summary_1[, frac_V1 := V1/sum(V1), by = c('method', 'prob', 'rep')]
    summary_1[, frac_V1 := as.numeric(frac_V1)]

    summary_1[, p_val := p_i]
    p_results[[p_ind]] = summary_1

  }

  summary_all = do.call('rbind', p_results)
  summary_all[, p_val := as.factor(p_val)]

  return(list(raw = outer_list_res, summary = summary_all))

}



spatialDE_sim = test_patterns_spatialDE(gobject = VC_small,
                                        pattern_ids = pattern_ids,
                                        pattern_name = 'center_stripe',
                                        probs = c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99, 1),
                                        repetitions = 2,
                                        show_pattern = F,
                                        pattern_colors = c('in' = 'green', 'out' = 'red'),
                                        p_values = c(0.01, 0.05, 0.1))



summary_all = spatialDE_sim$summary

library(ggplot2)

pl = ggplot()
pl = pl + geom_point(data = summary_all, aes(x = as.factor(prob), color = ident, y = frac_V1))
pl = pl + geom_line(data = summary_all, aes(x = as.factor(prob), color = ident, y = frac_V1, group = ident))
pl = pl + facet_grid(~ p_val)
pl = pl + theme_bw()
pl

pl = ggplot()
pl = pl + geom_point(data = summary_all, aes(x = as.factor(prob), color = p_val, y = frac_V1))
pl = pl + geom_line(data = summary_all, aes(x = as.factor(prob), color = p_val, y = frac_V1, group = p_val))
pl = pl + facet_grid(~ ident)
pl = pl + theme_bw()
pl










# spark
testspark = spark(gobject = center_stripe,
                  return_object = 'spark')

testspark[combined_pvalue < 0.05]

spatGenePlot(center_stripe, expression_values = 'norm', genes = testspark[1:6]$genes,
             point_shape = 'border', point_border_stroke = 0.1,
             show_network = F, network_color = 'lightgrey', point_size = 2.5,
             cow_n_col = 2)

spatGenePlot(center_stripe, expression_values = 'norm', genes = tail(testspark$genes, 6),
             point_shape = 'border', point_border_stroke = 0.1,
             show_network = F, network_color = 'lightgrey', point_size = 2.5,
             cow_n_col = 2)












# silhouette
silh_spatialgenes = silhouetteRank(gobject = VC_small_sim) # TODO: suppress print output

spatGenePlot(VC_small_sim, expression_values = 'norm', genes = silh_spatialgenes[1:6]$genes,
             point_shape = 'border', point_border_stroke = 0.1,
             show_network = F, network_color = 'lightgrey', point_size = 2.5,
             cow_n_col = 2)

spatGenePlot(VC_small_sim, expression_values = 'norm', genes = tail(silh_spatialgenes$genes, 6),
             point_shape = 'border', point_border_stroke = 0.1,
             show_network = F, network_color = 'lightgrey', point_size = 2.5,
             cow_n_col = 2)















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



