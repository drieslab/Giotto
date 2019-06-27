




#' @title annotateSpatialNetwork
#' @name annotateSpatialNetwork
#' @description Annotate spatial network with cell metadata information.
#' @param gobject giotto object
#' @param spatial_network_name name of spatial network to use
#' @param cluster_column name of column to use for clusters
#' @return annotated network in data.table format
#' @export
#' @examples
#'     annotateSpatialNetwork(gobject)
annotateSpatialNetwork = function(gobject,
                                  spatial_network_name = 'spatial_network',
                                  cluster_column) {

  # get network
  if(!spatial_network_name %in% names(gobject@spatial_network)) {
    stop('\n spatial network with name: ', spatial_network_name, ' does not exist \n')
  }
  spatial_network = gobject@spatial_network[[spatial_network_name]]

  # cell metadata
  cell_metadata = pDataDT(gobject)
  if(!cluster_column %in% colnames(cell_metadata)) {
    stop('\n the cluster column does not exist in pDataDT(gobject) \n')
  }
  cluster_type_vector = cell_metadata[[cluster_column]]
  names(cluster_type_vector) = cell_metadata[['cell_ID']]


  spatial_network_annot = copy(spatial_network)
  spatial_network_annot[, to_cell_type := cluster_type_vector[to]]
  spatial_network_annot[, from_cell_type := cluster_type_vector[from]]
  spatial_network_annot[, type_int := ifelse(to_cell_type == from_cell_type, 'homo', 'hetero')]

  # specific direction
  spatial_network_annot[, from_to := paste0(from_cell_type,'-',to_cell_type)]

  # unified direction, due to 'sort'
  spatial_network_annot[, unified_int := paste(sort(c(from_cell_type, to_cell_type)), collapse = '-'), by = 1:nrow(spatial_network_annot)]

  return(spatial_network_annot)

}


#' @title make_simulated_network
#' @name make_simulated_network
#' @description Simulate random network.
#' @examples
#'     make_simulated_network(gobject)
make_simulated_network = function(gobject,
                                  spatial_network_name = 'spatial_network',
                                  cluster_column,
                                  number_of_simulations = 100) {


  spatial_network_annot = annotateSpatialNetwork(gobject = gobject,
                                                 spatial_network_name = spatial_network_name,
                                                 cluster_column = cluster_column)

  # remove double edges between same cells #
  spatial_network_annot[, unified_cells := paste(sort(c(to,from)), collapse = '-'), by = 1:nrow(spatial_network_annot)]
  spatial_network_annot = spatial_network_annot[!duplicated(unified_cells)]

  # create a simulated network
  length_ints = nrow(spatial_network_annot)
  s1_list = list()
  s2_list = list()
  for(sim in 1:number_of_simulations) {
    s1 = sample(spatial_network_annot$from_cell_type, size = length_ints)
    s1_list[[sim]] = s1

    s2 = sample(spatial_network_annot$to_cell_type, size = length_ints)
    s2_list[[sim]] = s2
  }

  s1_vector = do.call('c', s1_list)
  s2_vector = do.call('c', s2_list)
  round_vector = rep(x = 1:number_of_simulations, each = length_ints)
  round_vector = paste0('sim',round_vector)
  sample_dt = data.table(s1 = s1_vector, s2 = s2_vector, round = round_vector)

  uniq_sim_comb = unique(sample_dt[,.(s1,s2)])
  uniq_sim_comb[, unified_int := paste(sort(c(s1,s2)), collapse = '-'), by  = 1:nrow(uniq_sim_comb)]
  sample_dt[uniq_sim_comb, unified_int := unified_int, on = c(s1 = 's1', s2 = 's2')]
  sample_dt[, type_int := ifelse(s1 == s2, 'homo', 'hetero')]

  return(sample_dt)

}


#' @title cellProximityEnrichment
#' @name cellProximityEnrichment
#' @description Compute cell-cell interaction enrichment (observed vs expected)
#' @param gobject giotto object
#' @param spatial_network_name name of spatial network to use
#' @param cluster_column name of column to use for clusters
#' @param number_of_simulations number of simulations to create expected observations
#' @return Cell Proximity scores (CPscores) in data.table format
#' @export
#' @examples
#'     cellProximityEnrichment(gobject)
cellProximityEnrichment <- function(gobject,
                                    spatial_network_name = 'spatial_network',
                                    cluster_column,
                                    number_of_simulations = 100) {


  spatial_network_annot = annotateSpatialNetwork(gobject = gobject,
                                                 spatial_network_name = spatial_network_name,
                                                 cluster_column = cluster_column)

  # remove double edges between same cells #
  spatial_network_annot[, unified_cells := paste(sort(c(to,from)), collapse = '-'), by = 1:nrow(spatial_network_annot)]
  spatial_network_annot = spatial_network_annot[!duplicated(unified_cells)]

  sample_dt = make_simulated_network(gobject = gobject,
                                     spatial_network_name = spatial_network_name,
                                     cluster_column = cluster_column, number_of_simulations = number_of_simulations)

  # combine original and simulated network
  table_sim_results <- sample_dt[, table(round), by = c('unified_int', 'type_int')]
  table_sim_results[, orig := 'simulations']
  spatial_network_annot[, round := 'original']
  table_orig_results <- spatial_network_annot[, table(round), by  = c('unified_int', 'type_int')]
  table_orig_results[, orig := 'original']
  table_results <- rbind(table_orig_results, table_sim_results)


  # add missing combinations from original or simulations
  all_simulation_ints = as.character(unique(table_results[orig == 'simulations']$unified_int))
  all_original_ints = as.character(unique(table_results[orig == 'original']$unified_int))
  missing_in_original = all_simulation_ints[!all_simulation_ints %in% all_original_ints]
  missing_in_simulations = all_original_ints[!all_original_ints %in% all_simulation_ints]
  create_missing_for_original = table_results[unified_int %in% missing_in_original]
  create_missing_for_original = unique(create_missing_for_original[, c('orig', 'V1') := list('original', 0)])
  create_missing_for_simulations = table_results[unified_int %in% missing_in_simulations]
  create_missing_for_simulations = unique(create_missing_for_simulations[, c('orig', 'V1') := list('simulations', 0)])

  table_results <- do.call('rbind', list(table_results, create_missing_for_original, create_missing_for_simulations))


  ## p-values
  combo_list = rep(NA, length = length(unique(table_results$unified_int)))
  p_high = rep(NA, length = length(unique(table_results$unified_int)))
  p_low = rep(NA, length = length(unique(table_results$unified_int)))

  for(int_combo in 1:length(unique(table_results$unified_int))) {

    this_combo = as.character(unique(table_results$unified_int)[int_combo])

    sub = table_results[unified_int == this_combo]

    orig_value = sub[orig == 'original']$V1
    sim_values = sub[orig == 'simulations']$V1

    length_simulations = length(sim_values)
    if(length_simulations != number_of_simulations) {
      additional_length_needed = number_of_simulations-length_simulations
      length_simulations = c(length_simulations, rep(0, additional_length_needed))
    }

    p_orig_higher = 1 - (sum((orig_value+1) > (sim_values+1))/number_of_simulations)
    p_orig_lower = 1 - (sum((orig_value+1) < (sim_values+1))/number_of_simulations)

    combo_list[[int_combo]] = this_combo
    p_high[[int_combo]] = p_orig_higher
    p_low[[int_combo]] = p_orig_lower

  }
  res_pvalue_DT = data.table(unified_int = as.vector(combo_list), p_higher_orig = p_high, p_lower_orig = p_low)


  # depletion or enrichment in barplot format
  table_mean_results <- table_results[, .(mean(V1)), by = c('orig', 'unified_int', 'type_int')]
  table_mean_results_dc <- dcast.data.table(data = table_mean_results, formula = type_int+unified_int~orig, value.var = 'V1')
  table_mean_results_dc[, original := ifelse(is.na(original), 0, original)]
  table_mean_results_dc[, enrichm := log2((original+1)/(simulations+1))]


  table_mean_results_dc <- merge(table_mean_results_dc, res_pvalue_DT, by = 'unified_int')
  setorder(table_mean_results_dc, enrichm)
  table_mean_results_dc[, unified_int := factor(unified_int, unified_int)]
  table_mean_results_dc[, PI_value := ifelse(p_higher_orig <= p_lower_orig, -log10(p_higher_orig+(1/number_of_simulations))*enrichm, -log10(p_lower_orig+(1/number_of_simulations))*enrichm)]
  setorder(table_mean_results_dc, PI_value)

  # order
  table_mean_results_dc <- table_mean_results_dc[order(-PI_value)]
  table_mean_results_dc[, int_ranking := 1:.N]

  return(list(raw_sim_table = table_results, enrichm_res = table_mean_results_dc))

}




#' @title get_specific_interaction_gene_enrichment
#' @name get_specific_interaction_gene_enrichment
#' @description Computes gene enrichment between specified interaction
#' @examples
#'     get_specific_interaction_gene_enrichment()
get_specific_interaction_gene_enrichment <- function(sub_spatial_network,
                                                     source_col = 'source_clus', source_IDs = 'from',
                                                     neighb_col = 'neighb_clus', neighb_IDs = 'to',
                                                     expression_matrix,
                                                     interaction_name = 'to_specify',
                                                     cell_annotation,
                                                     annotation_ID = 'uniq_ID',
                                                     cell_type_col,
                                                     do_diff_test = T,
                                                     exclude_selected_cells_from_test = T) {

  int_cell_types = sort(unique(c(sub_spatial_network[[source_col]], sub_spatial_network[[neighb_col]])))

  cell_type_1 = int_cell_types[[1]]
  set1 = as.character(sub_spatial_network[get(neighb_col) == cell_type_1][[neighb_IDs]])
  set2 = as.character(sub_spatial_network[get(source_col) == cell_type_1][[source_IDs]])
  cell_type_1_ids = unique(c(set1, set2))

  # if the interaction is homotypic you can take the same cell ids
  if(length(int_cell_types) > 1) {
    cell_type_2 = int_cell_types[[2]]
    set3 = as.character(sub_spatial_network[get(neighb_col) == cell_type_2][[neighb_IDs]])
    set4 = as.character(sub_spatial_network[get(source_col) == cell_type_2][[source_IDs]])
    cell_type_2_ids = unique(c(set3, set4))
  } else {
    cell_type_2 = cell_type_1
    cell_type_2_ids = cell_type_1_ids
  }

  nr_cell_type_1 = length(cell_type_1_ids)
  nr_cell_type_2 = length(cell_type_2_ids)

  # create average scores
  if(nr_cell_type_1 == 1) {
    cell_type_1_ids_matrix = expression_matrix[, colnames(expression_matrix) %in% cell_type_1_ids]
    average_cell_type_1_ids = cell_type_1_ids_matrix
  } else {
    cell_type_1_ids_matrix = expression_matrix[, colnames(expression_matrix) %in% cell_type_1_ids]
    average_cell_type_1_ids = rowMeans(cell_type_1_ids_matrix)
  }
  if(nr_cell_type_2 == 1) {
    cell_type_2_ids_matrix = expression_matrix[, colnames(expression_matrix) %in% cell_type_2_ids]
    average_cell_type_2_ids = cell_type_2_ids_matrix
  } else {
    cell_type_2_ids_matrix = expression_matrix[, colnames(expression_matrix) %in% cell_type_2_ids]
    average_cell_type_2_ids = rowMeans(cell_type_2_ids_matrix)
  }
  combined_1_2_ids = (average_cell_type_1_ids+average_cell_type_2_ids)/2


  # global expression
  all_cell_type_1_ids = as.character(cell_annotation[get(cell_type_col) == cell_type_1][[annotation_ID]])
  all_cell_type_2_ids = as.character(cell_annotation[get(cell_type_col) == cell_type_2][[annotation_ID]])

  if(exclude_selected_cells_from_test == TRUE) {
    all_cell_type_1_ids <- all_cell_type_1_ids[!all_cell_type_1_ids %in% cell_type_1_ids]
    all_cell_type_2_ids <- all_cell_type_2_ids[!all_cell_type_2_ids %in% cell_type_2_ids]
  }

  nr_all_cell_type_1 = length(all_cell_type_1_ids)
  nr_all_cell_type_2 = length(all_cell_type_2_ids)

  all_cell_type_1_ids_matrix = expression_matrix[, colnames(expression_matrix) %in% all_cell_type_1_ids]
  all_cell_type_2_ids_matrix = expression_matrix[, colnames(expression_matrix) %in% all_cell_type_2_ids]




  # wilcox test and averages
  if(nr_cell_type_1 == 1) {
    if(do_diff_test == T) cell_type_1_test =  wilcox.test(x = as.vector(cell_type_1_ids_matrix), y = as.vector(all_cell_type_1_ids_matrix))$p.value
    average_all_cell_type_1_ids = rowMeans(all_cell_type_1_ids_matrix)
  } else {
    if(do_diff_test == T) cell_type_1_test = sapply(1:nrow(cell_type_1_ids_matrix), function(i) wilcox.test(as.vector(cell_type_1_ids_matrix[i,]), as.vector(all_cell_type_1_ids_matrix[i,]))$p.value)
    average_all_cell_type_1_ids = rowMeans(all_cell_type_1_ids_matrix)
  }

  if(nr_cell_type_2 == 1) {
    if(do_diff_test == T) cell_type_2_test =  wilcox.test(x = as.vector(cell_type_2_ids_matrix), y = as.vector(all_cell_type_2_ids_matrix))$p.value
    average_all_cell_type_2_ids = rowMeans(all_cell_type_2_ids_matrix)
  } else {
    if(do_diff_test == T) cell_type_2_test = sapply(1:nrow(cell_type_2_ids_matrix), function(i) wilcox.test(as.vector(cell_type_2_ids_matrix[i,]), as.vector(all_cell_type_2_ids_matrix[i,]))$p.value)
    average_all_cell_type_2_ids = rowMeans(all_cell_type_2_ids_matrix)
  }

  combined_all_1_2_ids = (average_all_cell_type_1_ids+average_all_cell_type_2_ids)/2


  gene_names = names(combined_1_2_ids)

  if(do_diff_test == T) {
    res_DT = as.data.table(do.call('cbind', list(genes = gene_names,
                                                 cell_expr_1 = average_cell_type_1_ids, cell_expr_2 = average_cell_type_2_ids,
                                                 comb_expr = combined_1_2_ids,
                                                 all_cell_expr_1 = average_all_cell_type_1_ids, all_cell_expr_2 = average_all_cell_type_2_ids,
                                                 all_comb_expr = combined_all_1_2_ids,
                                                 pval_1 = cell_type_1_test, pval_2 = cell_type_2_test)))
  } else {
    res_DT = as.data.table(do.call('cbind', list(genes = gene_names,
                                                 cell_expr_1 = average_cell_type_1_ids, cell_expr_2 = average_cell_type_2_ids,
                                                 comb_expr = combined_1_2_ids,
                                                 all_cell_expr_1 = average_all_cell_type_1_ids, all_cell_expr_2 = average_all_cell_type_2_ids,
                                                 all_comb_expr = combined_all_1_2_ids)))
  }



  res_DT[, c('cell_type_1', 'cell_type_2', 'interaction', 'nr_1', 'nr_2', 'all_nr_1', 'all_nr_2') := list(cell_type_1, cell_type_2, interaction_name,
                                                                                                          nr_cell_type_1, nr_cell_type_2,
                                                                                                          nr_all_cell_type_1, nr_all_cell_type_2)]

  return(res_DT)

}


#' @title get_interaction_gene_enrichment
#' @name get_interaction_gene_enrichment
#' @description Computes gene enrichment between all interactions
#' @examples
#'     get_interaction_gene_enrichment()
get_interaction_gene_enrichment <- function(spatial_network,
                                            unified_int_col = 'unified_int',
                                            source_col = 'source_clus', source_IDs = 'from',
                                            neighb_col = 'neighb_clus', neighb_IDs = 'to',
                                            expression_matrix,
                                            cell_annotation,
                                            annotation_ID = 'uniq_ID',
                                            cell_type_col,
                                            do_diff_test = T,
                                            exclude_selected_cells_from_test = T,
                                            verbose = T) {

  interactions = unique(as.character(spatial_network[[unified_int_col]]))

  result_list <- list()
  for(selected_int in interactions) {

    if(verbose == TRUE){
      cat('start ', selected_int, '\n')
    }

    sub_int_netw = spatial_network[get(unified_int_col) == selected_int]

    tempres = get_specific_interaction_gene_enrichment(sub_spatial_network = sub_int_netw,
                                                       source_col = source_col, source_IDs = source_IDs,
                                                       neighb_col = neighb_col, neighb_IDs = neighb_IDs,
                                                       expression_matrix = expression_matrix, interaction_name = selected_int,
                                                       cell_annotation = cell_annotation, annotation_ID = annotation_ID, cell_type_col = cell_type_col,
                                                       do_diff_test = do_diff_test, exclude_selected_cells_from_test = exclude_selected_cells_from_test)
    result_list[[selected_int]] <- tempres
  }

  final_res = do.call('rbind', result_list)
  return(final_res)

}


#' @title get_cell_to_cell_sorted_name_conversion
#' @name get_cell_to_cell_sorted_name_conversion
#' @description creates unified cell-cell interaction names
#' @examples
#'     get_cell_to_cell_sorted_name_conversion()
get_cell_to_cell_sorted_name_conversion <- function(all_cell_types) {

  # create all possible combinations
  all_cell_types = unique(all_cell_types)
  first_combn = as.data.table(t(combn(x = all_cell_types, m = 2)))
  sec_combn = data.table(V1 = first_combn$V2, V2 = first_combn$V1)
  self_comb = data.table(V1 = all_cell_types, V2 = all_cell_types)

  # create data.table
  mult_comb = do.call('rbind', list(first_combn, sec_combn, self_comb))
  mult_comb[, c('V1', 'V2') := list(as.character(V1), as.character(V2))]
  mult_comb[, given_name := paste(c(V1,V2), collapse = '-'), by = 1:nrow(mult_comb)]
  mult_comb[, uniq_name := paste(sort(c(V1,V2)), collapse = '-'), by = 1:nrow(mult_comb)]

  # named vector
  name_conversion = mult_comb$uniq_name
  names(name_conversion) = mult_comb$given_name

  return(name_conversion)
}



#' @title getCellProximityGeneScores
#' @name getCellProximityGeneScores
#' @description Compute cell-cell interaction enrichment (observed vs expected)
#' @param gobject giotto object
#' @param spatial_network_name name of spatial network to use
#' @param cluster_column name of column to use for clusters
#' @param expression_values expression values to use
#' @param fold_change_addendum constant to add when calculating log2 fold-change
#' @param in_two_directions shows enrichment in both directions: cell1-cell2, cell2-cell1
#' @param exclude_selected_cells_from_test exclude certain cells from test
#' @param verbose verbose
#' @return Cell Proximity Gene scores (CPGscores) in data.table format
#' @details Give more details ...
#' @export
#' @examples
#'     getCellProximityGeneScores(gobject)
getCellProximityGeneScores = function(gobject,
                                      spatial_network_name = 'spatial_network',
                                      cluster_column = 'louvain_clus.1',
                                      expression_values = c('normalized', 'scaled', 'custom'),
                                      fold_change_addendum = 0.1,
                                      in_two_directions = TRUE,
                                      exclude_selected_cells_from_test = F,
                                      verbose = T) {


  # 1. create annotated spatial network
  annot_spatial_network =  annotateSpatialNetwork(gobject = gobject,
                                                  spatial_network_name = spatial_network_name,
                                                  cluster_column = cluster_column)

  # 2. get expression values
  values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))
  expr_values = Giotto:::select_expression_values(gobject = gobject, values = values)

  # 3. get cell metadata
  cell_metadata = pDataDT(gobject)
  if(!cluster_column %in% colnames(cell_metadata)) {
    stop('\n selected cluster column does not exist \n')
  }

  # 4. calculate cell-cell interaction gene scores
  interaction_gene_scores = Giotto:::get_interaction_gene_enrichment(spatial_network = annot_spatial_network,
                                                                     unified_int_col = 'unified_int',
                                                                     source_col = 'from_cell_type', source_IDs = 'from',
                                                                     neighb_col = 'to_cell_type', neighb_IDs = 'to',
                                                                     expression_matrix = expr_values,
                                                                     cell_annotation = cell_metadata,
                                                                     annotation_ID = 'cell_ID', cell_type_col = cluster_column,
                                                                     do_diff_test = TRUE,
                                                                     exclude_selected_cells_from_test = exclude_selected_cells_from_test,
                                                                     verbose = verbose)

  # difference with spatial
  interaction_gene_scores[, diff_spat := as.numeric(comb_expr)-as.numeric(all_comb_expr)]
  interaction_gene_scores[, diff_spat_1 := as.numeric(cell_expr_1)-as.numeric(all_cell_expr_1)]
  interaction_gene_scores[, diff_spat_2 := as.numeric(cell_expr_2)-as.numeric(all_cell_expr_2)]
  setorder(interaction_gene_scores, diff_spat)

  # fold-change with spatial
  interaction_gene_scores[, log2fc_spat_1 := log2((as.numeric(cell_expr_1)+fold_change_addendum)/(as.numeric(all_cell_expr_1)+fold_change_addendum))]
  interaction_gene_scores[, log2fc_spat_2 := log2((as.numeric(cell_expr_2)+fold_change_addendum)/(as.numeric(all_cell_expr_2)+fold_change_addendum))]
  interaction_gene_scores[, log2fc_spat := log2((as.numeric(comb_expr)+fold_change_addendum)/(as.numeric(all_comb_expr)+fold_change_addendum))]



  # expand scores to all possible cell-cell combinations
  # e.g. both astrocyte-NSC and NSC-astrocyte
  if(in_two_directions == TRUE) {
    CPGscore_first_direction  = interaction_gene_scores
    CPGscore_second_direction = interaction_gene_scores
    colnames(CPGscore_second_direction) = c('genes', 'cell_expr_2', 'cell_expr_1', 'comb_expr',
                                            'all_cell_expr_2', 'all_cell_expr_1', 'all_comb_expr',
                                            'pval_2', 'pval_1',
                                            'cell_type_2', 'cell_type_1', 'interaction',
                                            'nr_2', 'nr_1', 'all_nr_2', 'all_nr_1',
                                            'diff_spat', 'diff_spat_2', 'diff_spat_1',
                                            'log2fc_spat_1', 'log2fc_spat_2', 'log2fc_spat')
    CPGscore_second_direction = CPGscore_second_direction[, colnames(CPGscore_first_direction), with = F]
    CPGscore_second_direction[, interaction := paste0(cell_type_1,'-', cell_type_2)]
    CPGscore = rbind(CPGscore_first_direction, CPGscore_second_direction)
    CPGscore = unique(CPGscore)


    # homo or hetero-typic interaction
    CPGscore[, type_int := ifelse(cell_type_1 == cell_type_2, 'homo', 'hetero')]

    # create unique cell-cell interaction names
    name_conversion = get_cell_to_cell_sorted_name_conversion(all_cell_types = unique(c(CPGscore$cell_type_1,
                                                                                         CPGscore$cell_type_2)))
    CPGscore[, unified_int := name_conversion[[interaction]], by = 1:nrow(CPGscore)]
    CPGscore[, unif_int_rank := ifelse(interaction == unified_int, 1, 2)]



  } else {

    CPGscore = interaction_gene_scores
    CPGscore[, type_int := ifelse(cell_type_1 == cell_type_2, 'homo', 'hetero')]

  }

  return(CPGscore)

}




#' @title getGeneToGeneSelection
#' @name getGeneToGeneSelection
#' @description Compute gene-gene enrichment scores.
#' @param CPGscore CPGscore, output from getCellProximityGeneScores()
#' @param selected_genes select subset of genes
#' @param specific_genes_1 specific source genes (see details)
#' @param specific_genes_2 specific target genes (see details)
#' @param min_cells min number of cells threshold
#' @param min_pval p-value threshold
#' @param min_spat_diff spatial difference threshold
#' @param min_log2_fc log2 fold-change threshold
#' @param direction up or downregulation or both
#' @param fold_change_addendum constant to add when calculating log2 fold-change
#' @param verbose verbose
#' @return Gene to gene scores in data.table format
#' @details Give more details ...
#' @export
#' @examples
#'     getGeneToGeneSelection(CPGscore)
getGeneToGeneSelection <- function(CPGscore,
                                   selected_genes = NULL,
                                   specific_genes_1 = NULL,
                                   specific_genes_2 = NULL,
                                   min_cells = 5,
                                   min_pval = 0.05,
                                   min_spat_diff = 0.2,
                                   min_log2_fc = 0.5,
                                   direction = c('both', 'up', 'down'),
                                   fold_change_addendum = 0.1,
                                   verbose = TRUE) {


  direction = match.arg(direction, choices = c('both', 'up', 'down'))

  # remove redundant data
  unif_gene_scores = CPGscore[unif_int_rank == 1]

  # make sure p-value are numeric
  unif_gene_scores[, pval_1 := as.numeric(pval_1)]
  unif_gene_scores[, pval_2 := as.numeric(pval_2)]

  ## first filtering CPGscores
  if((!is.null(specific_genes_1) & !is.null(specific_genes_2))) {

    if(length(specific_genes_1) != length(specific_genes_2)) {
      stop('\n specific_genes_1 must have the same length as specific_genes_2')
    }

    unif_gene_scores = unif_gene_scores[genes %in% c(specific_genes_1, specific_genes_2)]

  } else if(!is.null(selected_genes)) {

    unif_gene_scores = unif_gene_scores[genes %in% c(selected_genes)]
  }


  # second filtering
  all_ints_selection = unif_gene_scores[(pval_1 <= min_pval & nr_1 >= min_cells & abs(diff_spat_1) >= min_spat_diff & abs(log2fc_spat_1) >= min_log2_fc) |
                                          (pval_2 <= min_pval & nr_2 >= min_cells & abs(diff_spat_2) >= min_spat_diff  & abs(log2fc_spat_2) >= min_log2_fc)]

  unif_gene_scores = all_ints_selection
  all_ints = unique(all_ints_selection[['unified_int']])


  savelist = list()

  for(cell_int in 1:length(all_ints)) {

    sel_int = all_ints[cell_int]

    if(verbose == TRUE) {
      cat('\n', sel_int, '\n')
    }

    # cell types
    sel_type_1 = strsplit(x = sel_int, split = '-')[[1]][1]
    sel_type_2 = strsplit(x = sel_int, split = '-')[[1]][2]

    # subsets of the data
    subset_type_1 = unif_gene_scores[cell_type_1 == sel_type_1][(pval_1 <= min_pval & nr_1 >= min_cells & abs(diff_spat_1) >= min_spat_diff & abs(log2fc_spat_1) >= min_log2_fc)][cell_type_2 == sel_type_2]
    subset_type_2 = unif_gene_scores[cell_type_2 == sel_type_2][(pval_2 <= min_pval & nr_2 >= min_cells & abs(diff_spat_2) >= min_spat_diff & abs(log2fc_spat_2) >= min_log2_fc)][cell_type_1 == sel_type_1]

    subset_type_1 = subset_type_1[,.(genes, cell_type_1, nr_1, cell_expr_1, all_nr_1, all_cell_expr_1, pval_1, diff_spat_1, log2fc_spat_1, unified_int)]
    setnames(subset_type_1, 'genes', 'genes_1')
    subset_type_2 = subset_type_2[,.(genes, cell_type_2, nr_2, cell_expr_2, all_nr_2, all_cell_expr_2, pval_2, diff_spat_2, log2fc_spat_2, unified_int)]
    setnames(subset_type_2, 'genes', 'genes_2')

    # merge data again
    mergetest = merge(subset_type_1, subset_type_2, by = 'unified_int', allow.cartesian = T)


    if(nrow(mergetest) > 0) {

      # create additional columns
      mergetest[, gene_gene := paste0(genes_1,'-',genes_2), by = 1:nrow(mergetest)]
      mergetest[, unif_gene_gene := paste(sort(c(genes_1, genes_2)), collapse = '-'), by = 1:nrow(mergetest)]

      # only keep specific combinations
      if((!is.null(specific_genes_1) & !is.null(specific_genes_2))) {

        LR_combo = paste0(specific_genes_1,'-', specific_genes_2)
        RL_combo = paste0(specific_genes_2,'-', specific_genes_1)
        mergetest = mergetest[gene_gene %in% c(LR_combo, RL_combo)]
      }

      savelist[[cell_int]] = mergetest
    }



  }

  finalres = do.call('rbind', savelist)

  # remove redundant homo-typic interaction data
  if(verbose == TRUE) cat('\n calculate additional information \n')

  finalres[, type_int := ifelse(cell_type_1 == cell_type_2, 'homo', 'hetero')]
  finalres[, unif_gene_gene_int := 1:.N, by = .(unified_int, unif_gene_gene)]
  finalres = finalres[type_int == 'hetero' | (type_int == 'homo' & unif_gene_gene_int == 1)]



  if(direction == 'both') {
    finalres = finalres
  } else if(direction == 'up') {
    finalres = finalres[(pval_1 <= min_pval & nr_1 >= min_cells & diff_spat_1 >= min_spat_diff & log2fc_spat_1 >= min_log2_fc) &
                          (pval_2 <= min_pval & nr_2 >= min_cells & diff_spat_2 >= min_spat_diff  & log2fc_spat_2 >= min_log2_fc)]
  } else if(direction == 'down') {
    finalres = finalres[(pval_1 <= min_pval & nr_1 >= min_cells & diff_spat_1 <= -min_spat_diff & log2fc_spat_1 <= -min_log2_fc) &
                          (pval_2 <= min_pval & nr_2 >= min_cells & diff_spat_2 <= -min_spat_diff  & log2fc_spat_2 <= -min_log2_fc)]
  }


  # provide additional data
  finalres[, all_cell_expr := as.numeric(all_cell_expr_1)+as.numeric(all_cell_expr_2), 1:nrow(finalres)]
  finalres[, spatial_cell_expr := as.numeric(cell_expr_1)+as.numeric(cell_expr_2), 1:nrow(finalres)]
  finalres[, diff_spatial_expr := spatial_cell_expr-all_cell_expr]
  finalres[, log2fc_spatial_expr := log2((spatial_cell_expr+fold_change_addendum)/(all_cell_expr+fold_change_addendum))]
  finalres[, all_cell_rank := rank(-all_cell_expr), by = unif_gene_gene]
  finalres[, spatial_cell_rank := rank(-spatial_cell_expr), by = unif_gene_gene]

  change_values = unlist(apply(finalres, MARGIN = 1, FUN = function(x) {
    direction_test(x, min_pval = min_pval)
  }))
  finalres[, change := change_values]

  return(finalres)


}





