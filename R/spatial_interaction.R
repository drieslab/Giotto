




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



#' @title getCellProximityGeneScores
#' @name getCellProximityGeneScores
#' @description Compute cell-cell interaction enrichment (observed vs expected)
#' @param gobject giotto object
#' @param spatial_network_name name of spatial network to use
#' @param cluster_column name of column to use for clusters
#' @param expression_values expression values to use
#' @param do_diff_test do differential test
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
                                             do_diff_test = T,
                                             exclude_selected_cells_from_test = F,
                                             verbose = T) {


  # 1. create annotated spatial network
  annot_spatial_network =  annotateSpatialNetwork(gobject = gobject,
                                                  spatial_network_name = spatial_network_name,
                                                  cluster_column = cluster_column)

  # 2. get expression values
  values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))
  expr_values = select_expression_values(gobject = gobject, values = values)

  # 3. get cell metadata
  cell_metadata = pDataDT(gobject)
  if(!cluster_column %in% colnames(cell_metadata)) {
    stop('\n selected cluster column does not exist \n')
  }

  # 4. calculate cell-cell interaction gene scores
  interaction_gene_scores = get_interaction_gene_enrichment(spatial_network = annot_spatial_network,
                                                            unified_int_col = 'unified_int',
                                                            source_col = 'from_cell_type', source_IDs = 'from',
                                                            neighb_col = 'to_cell_type', neighb_IDs = 'to',
                                                            expression_matrix = expr_values,
                                                            cell_annotation = cell_metadata,
                                                            annotation_ID = 'cell_ID', cell_type_col = cluster_column,
                                                            do_diff_test = do_diff_test,
                                                            exclude_selected_cells_from_test = exclude_selected_cells_from_test,
                                                            verbose = verbose)

  interaction_gene_scores[, diff_spat := as.numeric(comb_expr)-as.numeric(all_comb_expr)]
  interaction_gene_scores[, diff_spat_1 := as.numeric(cell_expr_1)-as.numeric(all_cell_expr_1)]
  interaction_gene_scores[, diff_spat_2 := as.numeric(cell_expr_2)-as.numeric(all_cell_expr_2)]
  setorder(interaction_gene_scores, diff_spat)


  return(interaction_gene_scores)

}



#' @title getGeneToGeneScores
#' @name getGeneToGeneScores
#' @description Compute gene-gene enrichment scores.
#' @param CPGscore CPGscore, output from getAverageCellProximityGeneScores()
#' @param selected_genes subset of genes to test
#' @param selected_cell_interactions subset of cell-cell interactions to test
#' @param specific_genes_1 specific source genes (see details)
#' @param specific_genes_2 specific target genes (see details)
#' @param verbose verbose
#' @return Gene to gene scores in data.table format
#' @details Give more details ...
#' @export
#' @examples
#'     getGeneToGeneScores(CPGscore)
getGeneToGeneScores <- function(CPGscore,
                                selected_genes = NULL,
                                selected_cell_interactions,
                                specific_genes_1 = NULL,
                                specific_genes_2 = NULL,
                                verbose = FALSE) {

  all_genes = unique(CPGscore[['genes']])

  if(!is.null(specific_genes_1) & !is.null(specific_genes_2) & length(specific_genes_1) == length(specific_genes_2)) {
    selected_genes = unique(c(specific_genes_1, specific_genes_2))
    cat('\n use specific gene-gene interactions \n')
  }


  first_cell_DT = CPGscore[genes %in% selected_genes & interaction %in% selected_cell_interactions][,.(interaction, cell_type_1, genes, cell_expr_1, nr_1, all_cell_expr_1, all_nr_1 )]
  setorder(first_cell_DT, genes, interaction)
  setnames(first_cell_DT, 'genes', 'genes_cell_1')
  second_cell_DT = CPGscore[genes %in% selected_genes & interaction %in% selected_cell_interactions][,.(interaction, cell_type_2, genes, cell_expr_2, nr_2, all_cell_expr_2, all_nr_2 )]
  setorder(second_cell_DT, genes, interaction)
  setnames(second_cell_DT, 'genes', 'genes_cell_2')



  # only use specific interactions
  if(!is.null(specific_genes_1) & !is.null(specific_genes_2) & length(specific_genes_1) == length(specific_genes_2)) {


    cat('\n start specific gene-gene interactions \n')

    savelist = list()
    for(gene_int in 1:length(specific_genes_1)) {

      # make sure that both genes are present
      if(specific_genes_1[gene_int] %in% all_genes & specific_genes_2[gene_int] %in% all_genes) {

        if(verbose == TRUE) {
          cat('\n gene 1 & 2: ', specific_genes_1[gene_int], ' & ',  specific_genes_2[gene_int],'\n')
        }

        both_genes = c(specific_genes_1[gene_int], specific_genes_2[gene_int])
        subset_1 = first_cell_DT[genes_cell_1 %in% both_genes]
        subset_2 = second_cell_DT[genes_cell_2 %in% both_genes]

        merge_subset = merge(subset_1, subset_2, by = 'interaction', allow.cartesian = T)
        merge_subset = merge_subset[genes_cell_1 != genes_cell_2]

        savelist[[gene_int]] = merge_subset
      }

    }

    results = do.call('rbind', savelist)

    # calculate differences
    comb_expr = (as.numeric(results$cell_expr_1) + as.numeric(results$cell_expr_2))/2
    results[, comb_expr := comb_expr]
    all_comb_expr = (as.numeric(results$all_cell_expr_1) + as.numeric(results$all_cell_expr_2))/2
    results[, all_comb_expr := all_comb_expr]
    results[, diff_spat_1 := as.numeric(cell_expr_1) - as.numeric(all_cell_expr_1)]
    results[, diff_spat_2 := as.numeric(cell_expr_2) - as.numeric(all_cell_expr_2)]
    results[, diff_spat := comb_expr - all_comb_expr]

    results[, gene_gene := paste0(genes_cell_1,'-',genes_cell_2)]

    return(results)

  } else {
    results = merge(first_cell_DT, second_cell_DT, by = 'interaction', allow.cartesian = T)

    # calculate differences
    comb_expr = (as.numeric(results$cell_expr_1) + as.numeric(results$cell_expr_2))/2
    results[, comb_expr := comb_expr]
    all_comb_expr = (as.numeric(results$all_cell_expr_1) + as.numeric(results$all_cell_expr_2))/2
    results[, all_comb_expr := all_comb_expr]
    results[, diff_spat_1 := as.numeric(cell_expr_1) - as.numeric(all_cell_expr_1)]
    results[, diff_spat_2 := as.numeric(cell_expr_2) - as.numeric(all_cell_expr_2)]
    results[, diff_spat := comb_expr - all_comb_expr]

    results[, gene_gene := paste0(genes_cell_1,'-',genes_cell_2)]

    return(results)
  }

}


