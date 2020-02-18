


# cell type proximity enrichment ####

#' @title annotateSpatialNetwork
#' @name annotateSpatialNetwork
#' @description Annotate spatial network with cell metadata information.
#' @param gobject giotto object
#' @param spatial_network_name name of spatial network to use
#' @param cluster_column name of column to use for clusters
#' @param create_full_network convert from reduced to full network representation
#' @return annotated network in data.table format
#' @export
#' @examples
#'     annotateSpatialNetwork(gobject)
annotateSpatialNetwork = function(gobject,
                                  spatial_network_name = 'spatial_network',
                                  cluster_column,
                                  create_full_network = F) {

  # get network
  if(!spatial_network_name %in% names(gobject@spatial_network)) {
    stop('\n spatial network with name: ', spatial_network_name, ' does not exist \n')
  }
  spatial_network = gobject@spatial_network[[spatial_network_name]]

  if(create_full_network == TRUE) {
    spatial_network = Giotto:::convert_to_full_spatial_network(spatial_network)
    setnames(spatial_network,
             old = c('source', 'target', 'source_begin', 'source_end', 'target_begin', 'target_end'),
             new = c('from', 'to', 'sdimx_begin', 'sdimy_begin', 'sdimx_end', 'sdimy_end'))
  }



  # cell metadata
  cell_metadata = pDataDT(gobject)
  if(!cluster_column %in% colnames(cell_metadata)) {
    stop('\n the cluster column does not exist in pDataDT(gobject) \n')
  }
  cluster_type_vector = cell_metadata[[cluster_column]]
  names(cluster_type_vector) = cell_metadata[['cell_ID']]


  spatial_network_annot = data.table::copy(spatial_network)
  spatial_network_annot[, to_cell_type := cluster_type_vector[to]]
  spatial_network_annot[, from_cell_type := cluster_type_vector[from]]
  spatial_network_annot[, type_int := ifelse(to_cell_type == from_cell_type, 'homo', 'hetero')]

  # specific direction
  spatial_network_annot[, from_to := paste0(from_cell_type,'-',to_cell_type)]

  # unified direction, due to 'sort'
  spatial_network_annot = Giotto:::sort_combine_two_DT_columns(spatial_network_annot,
                                                               column1 = 'from_cell_type', column2 = 'to_cell_type',
                                                               myname = 'unified_int')

  #spatial_network_annot[, unified_int := paste(sort(c(from_cell_type, to_cell_type)), collapse = '--'), by = 1:nrow(spatial_network_annot)]

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
  spatial_network_annot = Giotto:::sort_combine_two_DT_columns(spatial_network_annot,
                                                               column1 = 'from', column2 = 'to',
                                                               myname = 'unified_cells')
  spatial_network_annot = spatial_network_annot[!duplicated(unified_cells)]

  # create a simulated network
  length_ints = nrow(spatial_network_annot)
  s1_list = list()
  s2_list = list()

  all_cell_type = c(spatial_network_annot$from_cell_type, spatial_network_annot$to_cell_type)
  middle_point = length(all_cell_type)/2

  for(sim in 1:number_of_simulations) {

    reshuffled_all_cell_type = sample(x = all_cell_type, size = length(all_cell_type), replace = F)

    new_from_cell_type = reshuffled_all_cell_type[1:middle_point]
    s1_list[[sim]] = new_from_cell_type

    new_to_cell_type = reshuffled_all_cell_type[(middle_point+1):length(all_cell_type)]
    s2_list[[sim]] = new_to_cell_type

  }

  s1_vector = do.call('c', s1_list)
  s2_vector = do.call('c', s2_list)
  round_vector = rep(x = 1:number_of_simulations, each = length_ints)
  round_vector = paste0('sim',round_vector)
  sample_dt = data.table::data.table(s1 = s1_vector, s2 = s2_vector, round = round_vector)

  uniq_sim_comb = unique(sample_dt[,.(s1,s2)])
  uniq_sim_comb[, unified_int := paste(sort(c(s1,s2)), collapse = '--'), by  = 1:nrow(uniq_sim_comb)]
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
#' @return List of cell Proximity scores (CPscores) in data.table format. The first
#' data.table (raw_sim_table) shows the raw observations of both the original and
#' simulated networks. The second data.table (enrichm_res) shows the enrichment results.
#' @details Spatial proximity enrichment or depletion between pairs of cell types
#' is calculated by calculating the observed over the expected frequency
#' of cell-cell proximity interactions. The expected frequency is the average frequency
#' calculated from a number of spatial network simulations. Each individual simulation is
#' obtained by reshuffling the cell type labels of each node (cell)
#' in the spatial network.
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
  spatial_network_annot[, unified_cells := paste(sort(c(to,from)), collapse = '--'), by = 1:nrow(spatial_network_annot)]
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
  res_pvalue_DT = data.table::data.table(unified_int = as.vector(combo_list), p_higher_orig = p_high, p_lower_orig = p_low)


  # depletion or enrichment in barplot format
  table_mean_results <- table_results[, .(mean(V1)), by = c('orig', 'unified_int', 'type_int')]
  table_mean_results_dc <- data.table::dcast.data.table(data = table_mean_results, formula = type_int+unified_int~orig, value.var = 'V1')
  table_mean_results_dc[, original := ifelse(is.na(original), 0, original)]
  table_mean_results_dc[, enrichm := log2((original+1)/(simulations+1))]


  table_mean_results_dc <- merge(table_mean_results_dc, res_pvalue_DT, by = 'unified_int')
  data.table::setorder(table_mean_results_dc, enrichm)
  table_mean_results_dc[, unified_int := factor(unified_int, unified_int)]
  table_mean_results_dc[, PI_value := ifelse(p_higher_orig <= p_lower_orig,
                                             -log10(p_higher_orig+(1/number_of_simulations))*enrichm,
                                             -log10(p_lower_orig+(1/number_of_simulations))*enrichm)]
  data.table::setorder(table_mean_results_dc, PI_value)

  # order
  table_mean_results_dc <- table_mean_results_dc[order(-PI_value)]
  table_mean_results_dc[, int_ranking := 1:.N]

  return(list(raw_sim_table = table_results, enrichm_res = table_mean_results_dc))

}



#' @title addCellIntMetadata
#' @name addCellIntMetadata
#' @description Creates an additional metadata column with information about interacting and non-interacting cell types of the
#' selected cell-cell interaction.
#' @param gobject giotto object
#' @param spatial_network name of spatial network to use
#' @param cluster_column column of cell types
#' @param cell_interaction cell-cell interaction to use
#' @param name name for the new metadata column
#' @param return_gobject return an updated giotto object
#' @return Giotto object
#' @details This function will create an additional metadata column which selects interacting cell types for a specific cell-cell
#' interaction. For example, if you want to color interacting astrocytes and oligodendrocytes it will create a new metadata column with
#' the values "select_astrocytes", "select_oligodendrocytes", "other_astrocytes", "other_oligodendroyctes" and "other". Where "other" is all
#' other cell types found within the selected cell type column.
#' @export
#' @examples
#'     addCellIntMetadata(gobject)
addCellIntMetadata = function(gobject,
                              spatial_network = 'spatial_network',
                              cluster_column,
                              cell_interaction,
                              name = 'select_int',
                              return_gobject = TRUE) {

  if(is.null(spatial_network)) {
    stop('spatial_network must be provided, this must be an existing spatial network \n')
  }

  if(is.null(cluster_column)) {
    stop('cluster_column must be provided, this must be an existing cell metadata column, see pData(your_giotto_object) \n')
  }

  if(is.null(cell_interaction)) {
    stop('cell_interaction must be provided, this must be cell--cell interaction between cell types in cluster_column \n')
  }

  # create spatial network
  spatial_network_annot = annotateSpatialNetwork(gobject = gobject,
                                                 spatial_network_name = spatial_network,
                                                 cluster_column = cluster_column)

  # selected vs other cells
  selected_cells = unique(c(spatial_network_annot[unified_int == cell_interaction]$to,
                            spatial_network_annot[unified_int == cell_interaction]$from))

  cell_metadata = data.table::copy(pDataDT(gobject))
  cell_type_1 = strsplit(cell_interaction, split = '--')[[1]][1]
  cell_type_2 = strsplit(cell_interaction, split = '--')[[1]][2]

  cell_metadata[, c(name) := ifelse(!get(cluster_column) %in% c(cell_type_1, cell_type_2), 'other',
                                    ifelse(get(cluster_column) == cell_type_1 & cell_ID %in% selected_cells, paste0("select_", cell_type_1),
                                           ifelse(get(cluster_column) == cell_type_2 & cell_ID %in% selected_cells, paste0("select_", cell_type_2),
                                                  ifelse(get(cluster_column) == cell_type_1, paste0("other_", cell_type_1), paste0("other_", cell_type_2)))))]

  if(return_gobject == TRUE) {
    gobject@cell_metadata = cell_metadata
    return(gobject)
  } else {
    return(cell_metadata)
  }

}




# * ####
# NEW CPGscores ####

#' @title do_ttest
#' @name do_ttest
#' @description Performs t.test on subsets of a matrix
#' @examples
#'     do_ttest()
do_ttest = function(expr_values, select_ind, other_ind, adjust_method) {

  if(length(select_ind) == 1){mean_sel = expr_values[,select_ind]} else{mean_sel = rowMeans(expr_values[,select_ind])}
  if(length(other_ind) == 1){mean_all = expr_values[,other_ind]} else{mean_all = rowMeans(expr_values[,other_ind])}

  if(length(select_ind) == 1 | length(other_ind) == 1) {
    results = NaN
  } else {
    results = apply(expr_values, MARGIN = 1, function(x) {
      p.value = stats::t.test(x[select_ind], x[other_ind])$p.value
    })
  }

  # other info
  log2fc = log2((mean_sel+1)/(mean_all+1))
  diff = mean_sel - mean_all

  resultsDT = data.table('genes' = rownames(expr_values), 'sel' = mean_sel, 'other' = mean_all, 'log2fc' = log2fc, 'diff' = diff, 'p.value' = unlist(results))
  resultsDT[, p.value := ifelse(is.nan(p.value), 1, p.value)]
  resultsDT[, p.adj := p.adjust(p.value, method = adjust_method)]
  setorder(resultsDT, p.adj)

  return(resultsDT)
}

#' @title do_limmatest
#' @name do_limmatest
#' @description Performs limma t.test on subsets of a matrix
#' @examples
#'     do_limmatest()
do_limmatest = function(expr_values, select_ind, other_ind) {

  expr_values_subset = cbind(expr_values[,select_ind], expr_values[,other_ind])
  mygroups = c(rep('sel', length(select_ind)), rep('other', length(other_ind)))
  mygroups = factor(mygroups, levels = unique(mygroups))

  design = stats::model.matrix(~0+mygroups)
  colnames(design) = levels(mygroups)
  fit = limma::lmFit(expr_values_subset, design = design)


  cont.matrix = limma::makeContrasts(
    sel_vs_other =  sel-other,
    levels = design
  )

  fitcontrast = limma::contrasts.fit(fit,  cont.matrix)
  fitc_ebayes = limma::eBayes(fitcontrast)

  # limma to DT
  limma_result = limma::topTable(fitc_ebayes, coef = 1,number = 100000, confint = T)
  limmaDT = data.table::as.data.table(limma_result); limmaDT[, genes := rownames(limma_result)]

  # other info
  if(length(select_ind) == 1){mean_sel = expr_values[,select_ind]} else{mean_sel = rowMeans(expr_values[,select_ind])}
  if(length(other_ind) == 1){mean_all = expr_values[,other_ind]} else{mean_all = rowMeans(expr_values[,other_ind])}

  log2fc = log2((mean_sel+1)/(mean_all+1))
  diff = mean_sel - mean_all

  tempDT = data.table::data.table('genes' = rownames(expr_values),
                                  'sel'= mean_sel,
                                  'other' = mean_all,
                                  'log2fc' = log2fc,
                                  'diff' =  diff)
  limmaDT = data.table:::merge.data.table(limmaDT, tempDT, by = 'genes')
  limmaDT = limmaDT[,.(genes, sel, other, log2fc, diff, P.Value, adj.P.Val)]
  colnames(limmaDT) = c('genes', 'sel', 'other', 'log2fc', 'diff', 'p.value', 'p.adj')

  setorder(limmaDT, p.adj)

  return(limmaDT)

}

#' @title do_ttest
#' @name do_ttest
#' @description Performs wilcoxon on subsets of a matrix
#' @examples
#'     do_ttest()
do_wilctest = function(expr_values, select_ind, other_ind, adjust_method) {

  if(length(select_ind) == 1){mean_sel = expr_values[,select_ind]} else{mean_sel = rowMeans(expr_values[,select_ind])}
  if(length(other_ind) == 1){mean_all = expr_values[,other_ind]} else{mean_all = rowMeans(expr_values[,other_ind])}

  if(length(select_ind) == 1 | length(other_ind) == 1) {
    results = NaN
  } else {
    results = apply(expr_values, MARGIN = 1, function(x) {
      p.value = stats::wilcox.test(x[select_ind], x[other_ind])$p.value
    })
  }

  # other info
  log2fc = log2((mean_sel+1)/(mean_all+1))
  diff = mean_sel - mean_all

  resultsDT = data.table('genes' = rownames(expr_values), 'sel' = mean_sel, 'other' = mean_all, 'log2fc' = log2fc, 'diff' = diff, 'p.value' = unlist(results))
  resultsDT[, p.value := ifelse(is.nan(p.value), 1, p.value)]
  resultsDT[, p.adj := p.adjust(p.value, method = adjust_method)]
  setorder(resultsDT, p.adj)

  return(resultsDT)

}

#' @title do_permuttest_original
#' @name do_permuttest_original
#' @description calculate original values
#' @examples
#'     do_permuttest_original()
do_permuttest_original = function(expr_values, select_ind, other_ind, name = 'orig') {

  if(length(select_ind) == 1){mean_sel = expr_values[,select_ind]} else{mean_sel = rowMeans(expr_values[,select_ind])}
  if(length(other_ind) == 1){mean_all = expr_values[,other_ind]} else{mean_all = rowMeans(expr_values[,other_ind])}

  log2fc = log2((mean_sel+1)/(mean_all+1))
  diff = mean_sel - mean_all

  resultsDT = data.table('sel' = mean_sel, 'other' = mean_all, 'log2fc' = log2fc, 'diff' =  diff)
  resultsDT[, genes := rownames(expr_values)]
  resultsDT[, name := name]

  return(resultsDT)

}

#' @title do_permuttest_random
#' @name do_permuttest_random
#' @description calculate random values
#' @examples
#'     do_permuttest_random()
do_permuttest_random = function(expr_values, select_ind, other_ind, name = 'perm_1') {

  l_select_ind = length(select_ind)
  l_other_ind = length(other_ind)

  all_ind = c(select_ind, other_ind)

  random_select = sample(all_ind, size = l_select_ind, replace = F)
  random_other = all_ind[!all_ind %in% random_select]

  # alternative
  if(length(select_ind) == 1){mean_sel = expr_values[,random_select]} else{mean_sel = rowMeans(expr_values[,random_select])}
  if(length(other_ind) == 1){mean_all = expr_values[,random_other]} else{mean_all = rowMeans(expr_values[,random_other])}

  log2fc = log2((mean_sel+1)/(mean_all+1))
  diff = mean_sel - mean_all

  resultsDT = data.table('sel' = mean_sel, 'other' = mean_all, 'log2fc' = log2fc, 'diff' = diff)
  resultsDT[, genes := rownames(expr_values)]
  resultsDT[, name := name]

  return(resultsDT)

}

#' @title do_multi_permuttest_random
#' @name do_multi_permuttest_random
#' @description calculate multiple random values
#' @examples
#'     do_multi_permuttest_random()
do_multi_permuttest_random = function(expr_values, select_ind, other_ind, n = 100, cores = 2) {

  result = parallel::mclapply(mc.cores = cores, X = 1:n, FUN = function(x) {

    perm_rand = do_permuttest_random(expr_values = expr_values, select_ind = select_ind, other_ind = other_ind, name = paste0('perm_', x))

  })

  final_result = do.call('rbind', result)

}

#' @title do_permuttest_random
#' @name do_permuttest_random
#' @description Performs permutation test on subsets of a matrix
#' @examples
#'     do_permuttest_random()
do_permuttest = function(expr_values, select_ind, other_ind, n_perm = 100, adjust_method = 'fdr', cores = 2) {


  ## original data
  original = do_permuttest_original(expr_values = expr_values, select_ind = select_ind, other_ind = other_ind, name = 'orig')

  ## random permutations
  random_perms = do_multi_permuttest_random(expr_values = expr_values, n = n_perm,
                                            select_ind = select_ind, other_ind = other_ind, cores = cores)

  ##
  random_perms[, log2fc_diff := rep(original$log2fc, n_perm) - log2fc]
  random_perms[, c('perm_sel', 'perm_other', 'perm_log2fc', 'perm_diff') := list(mean(sel), mean(other), mean(log2fc), mean(diff)), by = genes]

  ## get p-values
  random_perms[, p_higher := sum(log2fc_diff > 0), by = genes]
  random_perms[, p_higher := 1-(p_higher/n_perm)]
  random_perms[, p_lower := sum(log2fc_diff < 0), by = genes]
  random_perms[, p_lower := 1-(p_lower/n_perm)]

  ## combine results permutation and original
  random_perms_res = unique(random_perms[,.(genes, perm_sel, perm_other, perm_log2fc, perm_diff, p_higher, p_lower)])
  results_m = data.table:::merge.data.table(random_perms_res, original[,.(genes, sel, other, log2fc, diff)], by = 'genes')

  # select lowest p-value and perform p.adj
  results_m[, p.value := ifelse(p_higher <= p_lower, p_higher, p_lower)]
  results_m[, p.adj := p.adjust(p.value, method = adjust_method)]

  results_m = results_m[,.(genes, sel, other, log2fc, diff, p.value, p.adj, perm_sel, perm_other, perm_log2fc, perm_diff)]
  setorder(results_m, p.adj, -log2fc)

  return(results_m)

}


#' @title do_cell_proximity_test
#' @name do_cell_proximity_test
#' @description Performs a selected differential test on subsets of a matrix
#' @examples
#'     do_cell_proximity_test()
do_cell_proximity_test = function(expr_values,
                                  select_ind, other_ind,
                                  diff_test = c('permutation', 'limma', 't.test', 'wilcox'),
                                  n_perm = 100,
                                  adjust_method = c("bonferroni","BH", "holm", "hochberg", "hommel",
                                                    "BY", "fdr", "none"),
                                  cores = 2) {

  # get parameters
  diff_test = match.arg(diff_test, choices = c('permutation', 'limma', 't.test', 'wilcox'))
  adjust_method = match.arg(adjust_method, choices = c("bonferroni","BH", "holm", "hochberg", "hommel",
                                                       "BY", "fdr", "none"))


  if(diff_test == 'permutation') {
    result = do_permuttest(expr_values = expr_values,
                           select_ind = select_ind, other_ind = other_ind,
                           n_perm = n_perm, adjust_method = adjust_method, cores = cores)

  } else if(diff_test == 'limma') {
    result = do_limmatest(expr_values = expr_values, select_ind = select_ind, other_ind = other_ind)

  } else if(diff_test == 't.test') {
    result = do_ttest(expr_values = expr_values, select_ind = select_ind, other_ind = other_ind, adjust_method = adjust_method)

  } else if(diff_test == 'wilcox') {
    result = do_wilctest(expr_values = expr_values, select_ind = select_ind, other_ind = other_ind, adjust_method = adjust_method)

  }

  return(result)

}



#' @title findCellProximityGenes_per_interaction
#' @name findCellProximityGenes_per_interaction
#' @description Identifies genes that are differentially expressed due to proximity to other cell types.
#' @examples
#'     findCellProximityGenes_per_interaction()
findCellProximityGenes_per_interaction = function(expr_values,
                                                  cell_metadata,
                                                  annot_spatnetwork,
                                                  sel_int,
                                                  minimum_unique_cells = 1,
                                                  minimum_unique_int_cells = 1,
                                                  exclude_selected_cells_from_test = T,
                                                  diff_test = c('permutation', 'limma', 't.test', 'wilcox'),
                                                  adjust_method = 'bonferroni',
                                                  nr_permutations = 100,
                                                  cores = 1) {

  # select test to perform
  diff_test = match.arg(arg = diff_test, choices = c('permutation', 'limma', 't.test', 'wilcox'))

  # select subnetwork
  sub_spatnetwork = annot_spatnetwork[unified_int == sel_int]

  # unique cell types
  unique_cell_types = unique(c(sub_spatnetwork$to_cell_type, sub_spatnetwork$from_cell_type))

  if(length(unique_cell_types) == 2) {

    first_cell_type = unique_cell_types[1]
    second_cell_type = unique_cell_types[2]

    # first cell type ids
    to1 = sub_spatnetwork[to_cell_type == first_cell_type][['to']]
    from1 = sub_spatnetwork[from_cell_type == first_cell_type][['from']]
    cell1_ids = unique(c(to1, from1))

    # second cell type ids
    to2 = sub_spatnetwork[to_cell_type == second_cell_type][['to']]
    from2 = sub_spatnetwork[from_cell_type == second_cell_type][['from']]
    cell2_ids = unique(c(to2, from2))

    ## all cell ids
    all_cell1 = cell_metadata[cell_types == first_cell_type][['cell_ID']]
    all_cell2 = cell_metadata[cell_types == second_cell_type][['cell_ID']]

    ## exclude selected
    if(exclude_selected_cells_from_test == TRUE) {
      all_cell1 = all_cell1[!all_cell1 %in% cell1_ids]
      all_cell2 = all_cell2[!all_cell2 %in% cell2_ids]
    }

    ## FOR CELL TYPE 1
    sel_ind1 = which(colnames(expr_values) %in% cell1_ids)
    all_ind1 = which(colnames(expr_values) %in% all_cell1)

    ## FOR CELL TYPE 2
    sel_ind2 = which(colnames(expr_values) %in% cell2_ids)
    all_ind2 = which(colnames(expr_values) %in% all_cell2)


    ## do not continue if too few cells ##
    if(length(sel_ind1) < minimum_unique_cells | length(all_ind1) < minimum_unique_cells |
       length(sel_ind2) < minimum_unique_int_cells) {
      result_cell_1 = NULL
    } else {
      result_cell_1 = do_cell_proximity_test(expr_values = expr_values,
                                             select_ind = sel_ind1,
                                             other_ind = all_ind1,
                                             diff_test = diff_test,
                                             n_perm = nr_permutations,
                                             adjust_method = adjust_method,
                                             cores = cores)
      result_cell_1[, cell_type := first_cell_type]
      result_cell_1[, int_cell_type := second_cell_type]
      result_cell_1[, nr_select := length(sel_ind1)]
      result_cell_1[, int_nr_select := length(sel_ind2)]
      result_cell_1[, nr_other := length(all_ind1)]
      result_cell_1[, int_nr_other := length(all_ind2)]

    }


    ## do not continue if too few cells ##
    if(length(sel_ind2) < minimum_unique_cells | length(all_ind2) < minimum_unique_cells |
       length(sel_ind1) < minimum_unique_int_cells) {
      result_cell_2 = NULL
    } else {
      result_cell_2 = do_cell_proximity_test(expr_values = expr_values,
                                             select_ind = sel_ind2, other_ind = all_ind2,
                                             diff_test = diff_test,
                                             n_perm = nr_permutations,
                                             adjust_method = adjust_method,
                                             cores = cores)
      result_cell_2[, cell_type := second_cell_type]
      result_cell_2[, int_cell_type := first_cell_type]
      result_cell_2[, nr_select := length(sel_ind2)]
      result_cell_2[, int_nr_select := length(sel_ind1)]
      result_cell_2[, nr_other := length(all_ind2)]
      result_cell_2[, int_nr_other := length(all_ind1)]

    }


    ## COMBINE

    if(is.null(result_cell_1) & is.null(result_cell_2)) {
      return(NULL)
    } else {
      result_cells = rbind(result_cell_1, result_cell_2)
    }

  } else if(length(unique_cell_types) == 1) {

    first_cell_type = unique_cell_types[1]

    # first cell type ids
    to1 = sub_spatnetwork[to_cell_type == first_cell_type][['to']]
    from1 = sub_spatnetwork[from_cell_type == first_cell_type][['from']]
    cell1_ids = unique(c(to1, from1))

    ## all cell ids
    all_cell1 = cell_metadata[cell_types == first_cell_type][['cell_ID']]

    ## exclude selected
    if(exclude_selected_cells_from_test == TRUE) {
      all_cell1 = all_cell1[!all_cell1 %in% cell1_ids]
    }

    ## FOR CELL TYPE 1
    sel_ind1 = which(colnames(expr_values) %in% cell1_ids)
    all_ind1 = which(colnames(expr_values) %in% all_cell1)


    ## do not continue if too few cells ##
    if(length(sel_ind1) < minimum_unique_cells | length(all_ind1) < minimum_unique_cells) {
      return(NULL)
    }

    result_cells = do_cell_proximity_test(expr_values = expr_values,
                                          select_ind = sel_ind1, other_ind = all_ind1,
                                          diff_test = diff_test,
                                          n_perm = nr_permutations,
                                          adjust_method = adjust_method,
                                          cores = cores)

    result_cells[, cell_type := first_cell_type]
    result_cells[, int_cell_type := first_cell_type]
    result_cells[, nr_select := length(sel_ind1)]
    result_cells[, int_nr_select := length(sel_ind1)]
    result_cells[, nr_other := length(all_ind1)]
    result_cells[, int_nr_other := length(all_ind1)]

  }

  result_cells[, unif_int := sel_int]

  return(result_cells)

}





#' @title findCellProximityGenes
#' @name findCellProximityGenes
#' @description Identifies genes that are differentially expressed due to proximity to other cell types.
#' @param gobject giotto object
#' @param expression_values expression values to use
#' @param cluster_column name of column to use for cell types
#' @param spatial_network_name name of spatial network to use
#' @param minimum_unique_cells minimum number of target cells required
#' @param minimum_unique_int_cells minimum number of interacting cells required
#' @param diff_test which differential expression test
#' @param adjust_method which method to adjust p-values
#' @param nr_permutations number of permutations if diff_test = permutation
#' @param exclude_selected_cells_from_test exclude interacting cells other cells
#' @param do_parallel run calculations in parallel with mclapply
#' @param cores number of cores to use if do_parallel = TRUE
#' @return cpgObject that contains the differential gene scores
#' @details Function to calculate if genes are differentially expressed in cell types
#'  when they interact (approximated by physical proximity) with other cell types.
#'  The results data.table in the cpgObject contains - at least - the following columns:
#' \itemize{
#'  \item{genes:}{ All or selected list of tested genes}
#'  \item{sel:}{ average gene expression in the interacting cells from the target cell type }
#'  \item{other:}{ average gene expression in the NOT-interacting cells from the target cell type }
#'  \item{log2fc:}{ log2 fold-change between sel and other}
#'  \item{diff:}{ spatial expression difference between sel and other}
#'  \item{p.value:}{ associated p-value}
#'  \item{p.adj:}{ adjusted p-value}
#'  \item{cell_type:}{ target cell type}
#'  \item{int_cell_type:}{ interacting cell type}
#'  \item{nr_select:}{ number of cells for selected target cell type}
#'  \item{int_nr_select:}{ number of cells for interacting cell type}
#'  \item{nr_other:}{ number of other cells of selected target cell type}
#'  \item{int_nr_other:}{ number of other cells for interacting cell type}
#'  \item{unif_int:}{ cell-cell interaction}
#' }
#' @export
#' @examples
#'     findCellProximityGenes(gobject)
findCellProximityGenes = function(gobject,
                                  expression_values = 'normalized',
                                  cluster_column,
                                  spatial_network_name = 'spatial_network',
                                  minimum_unique_cells = 1,
                                  minimum_unique_int_cells = 1,
                                  diff_test = c('permutation', 'limma', 't.test', 'wilcox'),
                                  adjust_method = c("bonferroni","BH", "holm", "hochberg", "hommel",
                                                    "BY", "fdr", "none"),
                                  nr_permutations = 100,
                                  exclude_selected_cells_from_test = T,
                                  do_parallel = TRUE,
                                  cores = NA) {



  # expression values to be used
  values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))
  expr_values = Giotto:::select_expression_values(gobject = gobject, values = values)

  # difference test
  diff_test = match.arg(diff_test, choices = c('permutation', 'limma', 't.test', 'wilcox'))

  # p.adj test
  adjust_method = match.arg(adjust_method, choices = c("bonferroni","BH", "holm", "hochberg", "hommel",
                                                       "BY", "fdr", "none"))

  ## metadata
  cell_metadata = pDataDT(gobject)

  ## annotated spatial network
  annot_spatnetwork = annotateSpatialNetwork(gobject,
                                             spatial_network_name = spatial_network_name,
                                             cluster_column = cluster_column)


  all_interactions = unique(annot_spatnetwork$unified_int)

  if(do_parallel == TRUE) {

    # set number of cores
    if(is.na(cores) | !is.numeric(cores)) {
      cores = parallel::detectCores() - 1
    }

    fin_result = parallel::mclapply(mc.cores = cores, X = all_interactions, FUN = function(x) {

      tempres = findCellProximityGenes_per_interaction(expr_values = expr_values,
                                                       cell_metadata = cell_metadata,
                                                       annot_spatnetwork = annot_spatnetwork,
                                                       minimum_unique_cells = minimum_unique_cells,
                                                       minimum_unique_int_cells = minimum_unique_int_cells,
                                                       sel_int = x,
                                                       exclude_selected_cells_from_test = exclude_selected_cells_from_test,
                                                       diff_test = diff_test,
                                                       adjust_method = adjust_method,
                                                       nr_permutations = nr_permutations,
                                                       cores = 2)


    })

  } else {

    fin_result = list()

    for(i in 1:length(all_interactions)) {

      x = all_interactions[i]
      print(x)

      tempres = findCellProximityGenes_per_interaction(expr_values = expr_values,
                                                       cell_metadata = cell_metadata,
                                                       annot_spatnetwork = annot_spatnetwork,
                                                       minimum_unique_cells = minimum_unique_cells,
                                                       minimum_unique_int_cells = minimum_unique_int_cells,
                                                       sel_int = x,
                                                       exclude_selected_cells_from_test = exclude_selected_cells_from_test,
                                                       diff_test = diff_test,
                                                       adjust_method = adjust_method,
                                                       nr_permutations = nr_permutations,
                                                       cores = 2)

      fin_result[[i]] = tempres

    }


  }

  final_result = do.call('rbind', fin_result)

  final_result[, spec_int := paste0(cell_type,'--',int_cell_type)]
  final_result[, type_int := ifelse(cell_type == int_cell_type, 'homo', 'hetero')]


  #return(final_result)

  permutation_test = ifelse(diff_test == 'permutation', nr_permutations, 'no permutations')

  cpgObject = list(CPGscores = final_result,
                   Giotto_info = list('values' = values, 'cluster' = cluster_column,
                                      'spatial network' = spatial_network_name),
                   test_info = list('test' = diff_test,
                                    'p.adj' = adjust_method,
                                    'min cells' = minimum_unique_cells,
                                    'min interacting cells' = minimum_unique_int_cells,
                                    'exclude selected cells' = exclude_selected_cells_from_test,
                                    'perm' = permutation_test))
  class(cpgObject) = append(class(cpgObject), 'cpgObject')
  return(cpgObject)

}

#' @title findCPG
#' @name findCPG
#' @description Identifies genes that are differentially expressed due to proximity to other cell types.
#' @param gobject giotto object
#' @param expression_values expression values to use
#' @param cluster_column name of column to use for cell types
#' @param spatial_network_name name of spatial network to use
#' @param minimum_unique_cells minimum number of target cells required
#' @param minimum_unique_int_cells minimum number of interacting cells required
#' @param diff_test which differential expression test
#' @param adjust_method which method to adjust p-values
#' @param nr_permutations number of permutations if diff_test = permutation
#' @param exclude_selected_cells_from_test exclude interacting cells other cells
#' @param do_parallel run calculations in parallel with mclapply
#' @param cores number of cores to use if do_parallel = TRUE
#' @return cpgObject that contains the differential gene scores
#' @details Function to calculate if genes are differentially expressed in cell types
#'  when they interact (approximated by physical proximity) with other cell types.
#'  The results data.table in the cpgObject contains - at least - the following columns:
#' \itemize{
#'  \item{genes:}{ All or selected list of tested genes}
#'  \item{sel:}{ average gene expression in the interacting cells from the target cell type }
#'  \item{other:}{ average gene expression in the NOT-interacting cells from the target cell type }
#'  \item{log2fc:}{ log2 fold-change between sel and other}
#'  \item{diff:}{ spatial expression difference between sel and other}
#'  \item{p.value:}{ associated p-value}
#'  \item{p.adj:}{ adjusted p-value}
#'  \item{cell_type:}{ target cell type}
#'  \item{int_cell_type:}{ interacting cell type}
#'  \item{nr_select:}{ number of cells for selected target cell type}
#'  \item{int_nr_select:}{ number of cells for interacting cell type}
#'  \item{nr_other:}{ number of other cells of selected target cell type}
#'  \item{int_nr_other:}{ number of other cells for interacting cell type}
#'  \item{unif_int:}{ cell-cell interaction}
#' }
#' @export
#' @examples
#'     findCPG(gobject)
findCPG = function(gobject,
                   expression_values = 'normalized',
                   cluster_column,
                   spatial_network_name = 'spatial_network',
                   minimum_unique_cells = 1,
                   minimum_unique_int_cells = 1,
                   diff_test = c('permutation', 'limma', 't.test', 'wilcox'),
                   adjust_method = c("bonferroni","BH", "holm", "hochberg", "hommel",
                                     "BY", "fdr", "none"),
                   nr_permutations = 100,
                   exclude_selected_cells_from_test = T,
                   do_parallel = TRUE,
                   cores = NA) {


  findCellProximityGenes(gobject = gobject,
                         expression_values = expression_values,
                         cluster_column = cluster_column,
                         spatial_network_name = spatial_network_name,
                         minimum_unique_cells = minimum_unique_cells,
                         minimum_unique_int_cells = minimum_unique_int_cells,
                         diff_test = diff_test,
                         adjust_method = adjust_method,
                         nr_permutations = nr_permutations,
                         exclude_selected_cells_from_test = exclude_selected_cells_from_test,
                         do_parallel = do_parallel,
                         cores = cores)


}





#' @title filterCellProximityGenes
#' @name filterCellProximityGenes
#' @description Filter cell proximity gene scores.
#' @param cpgObject cell proximity gene score object
#' @param min_cells minimum number of target cell type
#' @param min_int_cells minimum number of interacting cell type
#' @param min_fdr minimum adjusted p-value
#' @param min_spat_diff minimum absolute spatial expression difference
#' @param min_log2_fc minimum absolute log2 fold-change
#' @param direction differential expression directions to keep
#' @return cpgObject that contains the filtered differential gene scores
#' @export
#' @examples
#'     filterCellProximityGenes(gobject)
filterCellProximityGenes = function(cpgObject,
                                    min_cells = 5,
                                    min_int_cells = 3,
                                    min_fdr = 0.05,
                                    min_spat_diff = 0.2,
                                    min_log2_fc = 0.5,
                                    direction = c('both', 'up', 'down')) {


  if(!'cpgObject' %in% class(cpgObject)) {
    stop('\n cpgObject needs to be the output from findCellProximityGenes() or findCPG() \n')
  }

  CPGscore = copy(cpgObject[['CPGscores']])

  # other parameters
  direction = match.arg(direction, choices = c('both', 'up', 'down'))

  selection_scores = CPGscore[(nr_select >= min_cells & int_nr_select >= min_int_cells &
                                 p.adj <= min_fdr &  abs(diff) >= min_spat_diff & abs(log2fc) >= min_log2_fc)]

  if(direction == 'both') {
    selection_scores = selection_scores
  } else if(direction == 'up') {
    selection_scores = selection_scores[log2fc >= min_log2_fc]
  } else if(direction == 'down') {
    selection_scores = selection_scores[log2fc <= -min_log2_fc]
  }


  newobj = copy(cpgObject)
  newobj[['CPGscores']] = selection_scores

  return(newobj)

}



#' @title filterCPG
#' @name filterCPG
#' @description Filter cell proximity gene scores.
#' @param cpgObject cell proximity gene score object
#' @param min_cells minimum number of target cell type
#' @param min_int_cells minimum number of interacting cell type
#' @param min_fdr minimum adjusted p-value
#' @param min_spat_diff minimum absolute spatial expression difference
#' @param min_log2_fc minimum absolute log2 fold-change
#' @param direction differential expression directions to keep
#' @return cpgObject that contains the filtered differential gene scores
#' @export
#' @examples
#'     filterCPG(gobject)
filterCPG = function(cpgObject,
                     min_cells = 5,
                     min_int_cells = 3,
                     min_fdr = 0.05,
                     min_spat_diff = 0.2,
                     min_log2_fc = 0.5,
                     direction = c('both', 'up', 'down')) {

  filterCellProximityGenes(cpgObject = cpgObject,
                           min_cells = min_cells,
                           min_int_cells = min_int_cells,
                           min_fdr = min_fdr,
                           min_spat_diff = min_spat_diff,
                           min_log2_fc = min_log2_fc,
                           direction = direction)

}



# * ####
# NEW GTGscores ####

#' @title combineCellProximityGenes_per_interaction
#' @name combineCellProximityGenes_per_interaction
#' @description Combine CPG scores per interaction
#' @examples
#'     combineCellProximityGenes_per_interaction()
combineCellProximityGenes_per_interaction =  function(cpgObject,
                                                      sel_int,
                                                      selected_genes = NULL,
                                                      specific_genes_1 = NULL,
                                                      specific_genes_2 = NULL,
                                                      min_cells = 5,
                                                      min_int_cells = 3,
                                                      min_fdr = 0.05,
                                                      min_spat_diff = 0,
                                                      min_log2_fc = 0.5) {

  if(!'cpgObject' %in% class(cpgObject)) {
    stop('\n cpgObject needs to be the output from findCellProximityGenes() or findCPG() \n')
  }

  CPGscore = copy(cpgObject[['CPGscores']])

  test_used = cpgObject[['test_info']][['test']]

  ## subset on selected interaction
  subset = CPGscore[unif_int == sel_int]

  ##  type of interactions
  type_interaction = unique(subset[['type_int']])

  ## first filtering CPGscores on genes
  if((!is.null(specific_genes_1) & !is.null(specific_genes_2))) {
    if(length(specific_genes_1) != length(specific_genes_2)) {
      stop('\n specific_genes_1 must have the same length as specific_genes_2')
    }
    subset = subset[genes %in% c(specific_genes_1, specific_genes_2)]
  } else if(!is.null(selected_genes)) {
    subset = subset[genes %in% c(selected_genes)]
  }

  ## find number of unique cell types
  unique_cell_types = unique(c(subset$cell_type, subset$int_cell_type))

  if(length(unique_cell_types) == 2) {

    ## CELL TYPE 1
    subset_cell_1 = subset[cell_type == unique_cell_types[1]]

    if(nrow(subset_cell_1) == 0) {

      if(test_used == 'permutation') {

        subset_cell_1 = data.table::data.table('genes_1' = subset[['genes']],
                                               'sel_1' = NA,
                                               'other_1' = NA,
                                               'log2fc_1' = NA,
                                               'diff_1' = NA,
                                               'p.value_1' = NA,
                                               'p.adj_1' = NA,
                                               'perm_sel_1' = NA,
                                               'perm_other_1' = NA,
                                               'perm_log2fc_1' = NA,
                                               'perm_diff_1' = NA,
                                               'cell_type_1' = unique_cell_types[1],
                                               'int_cell_type_1' = unique_cell_types[2],
                                               'nr_select_1' = NA,
                                               'nr_other_1' = NA,
                                               'unif_int' = subset[['unif_int']])

      } else {

        subset_cell_1 = data.table::data.table('genes_1' = subset[['genes']],
                                               'sel_1' = NA,
                                               'other_1' = NA,
                                               'log2fc_1' = NA,
                                               'diff_1' = NA,
                                               'p.value_1' = NA,
                                               'p.adj_1' = NA,
                                               'cell_type_1' = unique_cell_types[1],
                                               'int_cell_type_1' = unique_cell_types[2],
                                               'nr_select_1' = NA,
                                               'nr_other_1' = NA,
                                               'unif_int' = subset[['unif_int']])
      }


    } else {


      # filter on statistics
      subset_cell_1 = subset_cell_1[p.adj <= min_fdr]
      subset_cell_1 = subset_cell_1[nr_select >= min_cells]
      subset_cell_1 = subset_cell_1[int_nr_select >= min_int_cells]
      subset_cell_1 = subset_cell_1[abs(log2fc) >= min_log2_fc]
      subset_cell_1 = subset_cell_1[abs(diff) >= min_spat_diff]

      if(test_used == 'permutation') {
        # make it specific
        subset_cell_1 = subset_cell_1[,.(genes, sel, other, log2fc, diff, p.value, p.adj,
                                         perm_sel, perm_other, perm_log2fc, perm_diff,
                                         cell_type, int_cell_type, nr_select, nr_other, unif_int)]
        data.table::setnames(subset_cell_1, old = c('genes', 'sel', 'other', 'log2fc', 'diff', 'p.value', 'p.adj',
                                                    'perm_sel', 'perm_other', 'perm_log2fc', 'perm_diff',
                                                    'cell_type', 'int_cell_type', 'nr_select', 'nr_other'),
                             new = c('genes_1', 'sel_1', 'other_1', 'log2fc_1', 'diff_1', 'p.value_1', 'p.adj_1',
                                     'perm_sel_1', 'perm_other_1',  'perm_log2fc_1', 'perm_diff_1',
                                     'cell_type_1', 'int_cell_type_1', 'nr_select_1', 'nr_other_1'))

      } else {
        # make it specific
        subset_cell_1 = subset_cell_1[,.(genes, sel, other, log2fc, diff, p.value, p.adj, cell_type, int_cell_type, nr_select, nr_other, unif_int)]
        data.table::setnames(subset_cell_1, old = c('genes', 'sel', 'other', 'log2fc', 'diff', 'p.value', 'p.adj', 'cell_type', 'int_cell_type', 'nr_select', 'nr_other'),
                             new = c('genes_1', 'sel_1', 'other_1', 'log2fc_1', 'diff_1', 'p.value_1', 'p.adj_1', 'cell_type_1', 'int_cell_type_1', 'nr_select_1', 'nr_other_1'))

      }


    }




    ## CELL TYPE 2
    subset_cell_2 = subset[cell_type == unique_cell_types[2]]

    if(nrow(subset_cell_2) == 0) {

      if(test_used == 'permutation') {

        subset_cell_2 = data.table::data.table('genes_2' = subset[['genes']],
                                               'sel_2' = NA,
                                               'other_2' = NA,
                                               'log2fc_2' = NA,
                                               'diff_2' = NA,
                                               'p.value_2' = NA,
                                               'p.adj_2' = NA,
                                               'perm_sel_2' = NA,
                                               'perm_other_2' = NA,
                                               'perm_log2fc_2' = NA,
                                               'perm_diff_2' = NA,
                                               'cell_type_2' = unique_cell_types[2],
                                               'int_cell_type_2' = unique_cell_types[1],
                                               'nr_select_2' = NA,
                                               'nr_other_2' = NA,
                                               'unif_int' = subset[['unif_int']])
      } else {

        subset_cell_2 = data.table::data.table('genes_2' = subset[['genes']],
                                               'sel_2' = NA,
                                               'other_2' = NA,
                                               'log2fc_2' = NA,
                                               'diff_2' = NA,
                                               'p.value_2' = NA,
                                               'p.adj_2' = NA,
                                               'cell_type_2' = unique_cell_types[2],
                                               'int_cell_type_2' = unique_cell_types[1],
                                               'nr_select_2' = NA,
                                               'nr_other_2' = NA,
                                               'unif_int' = subset[['unif_int']])
      }



    } else {

      # filter on statistics
      subset_cell_2 = subset_cell_2[p.adj <= min_fdr]
      subset_cell_2 = subset_cell_2[nr_select >= min_cells]
      subset_cell_2 = subset_cell_2[int_nr_select >= min_int_cells]
      subset_cell_2 = subset_cell_2[abs(log2fc) >= min_log2_fc]
      subset_cell_2 = subset_cell_2[abs(diff) >= min_spat_diff]

      if(test_used == 'permutation') {

        subset_cell_2 = subset_cell_2[,.(genes, sel, other, log2fc, diff, p.value, p.adj,
                                         perm_sel, perm_other, perm_log2fc, perm_diff,
                                         cell_type, int_cell_type, nr_select, nr_other, unif_int)]
        data.table::setnames(subset_cell_2, old = c('genes', 'sel', 'other', 'log2fc', 'diff', 'p.value', 'p.adj',
                                                    'perm_sel', 'perm_other', 'perm_log2fc', 'perm_diff',
                                                    'cell_type', 'int_cell_type', 'nr_select', 'nr_other'),
                             new = c('genes_2', 'sel_2', 'other_2', 'log2fc_2', 'diff_2', 'p.value_2', 'p.adj_2',
                                     'perm_sel_2', 'perm_other_2', 'perm_log2fc_2', 'perm_diff_2',
                                     'cell_type_2', 'int_cell_type_2', 'nr_select_2', 'nr_other_2'))


      } else {
        subset_cell_2 = subset_cell_2[,.(genes, sel, other, log2fc, diff, p.value, p.adj, cell_type, int_cell_type, nr_select, nr_other, unif_int)]
        data.table::setnames(subset_cell_2, old = c('genes', 'sel', 'other', 'log2fc', 'diff', 'p.value', 'p.adj', 'cell_type', 'int_cell_type', 'nr_select', 'nr_other'),
                             new = c('genes_2', 'sel_2', 'other_2', 'log2fc_2', 'diff_2', 'p.value_2', 'p.adj_2', 'cell_type_2', 'int_cell_type_2', 'nr_select_2', 'nr_other_2'))

      }



    }

    merge_subsets = data.table:::merge.data.table(subset_cell_1, subset_cell_2, by = c('unif_int'), allow.cartesian = TRUE)

  } else if(length(unique_cell_types) == 1) {

    ## CELL TYPE 1
    subset_cell_1 = subset[cell_type == unique_cell_types[1]]

    # filter on statistics
    subset_cell_1 = subset_cell_1[p.adj <= min_fdr]
    subset_cell_1 = subset_cell_1[nr_select >= min_cells]
    subset_cell_1 = subset_cell_1[int_nr_select >= min_int_cells]
    subset_cell_1 = subset_cell_1[abs(log2fc) >= min_log2_fc]
    subset_cell_1 = subset_cell_1[abs(diff) >= min_spat_diff]

    # make it specific

    if(test_used == 'permutation') {
      subset_cell_1A = subset_cell_1[,.(genes, sel, other, log2fc, diff, p.value, p.adj,
                                        perm_sel, perm_other, perm_log2fc, perm_diff,
                                        cell_type, int_cell_type, nr_select, nr_other, unif_int)]
      data.table::setnames(subset_cell_1A, old = c('genes', 'sel', 'other', 'log2fc', 'diff', 'p.value', 'p.adj',
                                                   'perm_sel', 'perm_other', 'perm_log2fc', 'perm_diff',
                                                   'cell_type', 'int_cell_type', 'nr_select', 'nr_other'),
                           new = c('genes_1', 'sel_1', 'other_1', 'log2fc_1', 'diff_1', 'p.value_1', 'p.adj_1',
                                   'perm_sel_1', 'perm_other_1', 'perm_log2fc_1', 'perm_diff_1',
                                   'cell_type_1', 'int_cell_type_1', 'nr_select_1', 'nr_other_1'))

    } else {
      subset_cell_1A = subset_cell_1[,.(genes, sel, other, log2fc, diff, p.value, p.adj, cell_type, int_cell_type, nr_select, nr_other, unif_int)]
      data.table::setnames(subset_cell_1A, old = c('genes', 'sel', 'other', 'log2fc', 'diff', 'p.value', 'p.adj', 'cell_type', 'int_cell_type', 'nr_select', 'nr_other'),
                           new = c('genes_1', 'sel_1', 'other_1', 'log2fc_1', 'diff_1', 'p.value_1', 'p.adj_1', 'cell_type_1', 'int_cell_type_1', 'nr_select_1', 'nr_other_1'))

    }


    ## CELL TYPE 2

    if(test_used == 'permutation') {
      subset_cell_1B = subset_cell_1[,.(genes, sel, other, log2fc, diff, p.value, p.adj,
                                        perm_sel, perm_other, perm_log2fc, perm_diff,
                                        cell_type, int_cell_type, nr_select, nr_other, unif_int)]
      data.table::setnames(subset_cell_1B, old = c('genes', 'sel', 'other', 'log2fc', 'diff', 'p.value', 'p.adj',
                                                   'perm_sel', 'perm_other', 'perm_log2fc', 'perm_diff',
                                                   'cell_type', 'int_cell_type', 'nr_select', 'nr_other'),
                           new = c('genes_2', 'sel_2', 'other_2', 'log2fc_2', 'diff_2', 'p.value_2', 'p.adj_2',
                                   'perm_sel_2', 'perm_other_2', 'perm_log2fc_2', 'perm_diff_2',
                                   'cell_type_2', 'int_cell_type_2', 'nr_select_2', 'nr_other_2'))

    } else {
      subset_cell_1B = subset_cell_1[,.(genes, sel, other, log2fc, diff, p.value, p.adj, cell_type, int_cell_type, nr_select, nr_other, unif_int)]
      data.table::setnames(subset_cell_1B, old = c('genes', 'sel', 'other', 'log2fc', 'diff', 'p.value', 'p.adj', 'cell_type', 'int_cell_type', 'nr_select', 'nr_other'),
                           new = c('genes_2', 'sel_2', 'other_2', 'log2fc_2', 'diff_2', 'p.value_2', 'p.adj_2', 'cell_type_2', 'int_cell_type_2', 'nr_select_2', 'nr_other_2'))

    }

    merge_subsets = data.table:::merge.data.table(subset_cell_1A, subset_cell_1B, by = c('unif_int'), allow.cartesian = TRUE)


  }

  # restrict to gene combinations if needed
  if((!is.null(specific_genes_1) & !is.null(specific_genes_2))) {
    merge_subsets[, genes_combo := paste0(genes_1,'--',genes_2)]
    all_combos = c(paste0(specific_genes_1,'--', specific_genes_2),
                   paste0(specific_genes_2,'--', specific_genes_1))
    merge_subsets = merge_subsets[genes_combo %in% all_combos]
    merge_subsets[, genes_combo := NULL]
  }

  merge_subsets[, type_int := type_interaction]
  return(merge_subsets)

}


#' @title sort_combine_two_DT_columns
#' @name sort_combine_two_DT_columns
#' @description fast sorting and pasting of 2 character columns
#' @examples
#'     sort_combine_two_DT_columns()
sort_combine_two_DT_columns = function(DT, column1, column2, myname = 'unif_gene_gene') {


  # maybe faster with converting to factors??

  # make sure columns are character
  selected_columns = c(column1, column2)
  DT[,(selected_columns):= lapply(.SD, as.character), .SDcols = selected_columns]

  # convert characters into numeric values
  uniq_values = sort(unique(c(DT[[column1]], DT[[column2]])))
  uniq_values_num = 1:length(uniq_values)
  names(uniq_values_num) = uniq_values


  DT[,values_1_num := uniq_values_num[get(column1)]]
  DT[,values_2_num := uniq_values_num[get(column2)]]


  DT[, scolumn_1 := ifelse(values_1_num < values_2_num, get(column1), get(column2))]
  DT[, scolumn_2 := ifelse(values_1_num < values_2_num, get(column2), get(column1))]

  DT[, unif_sort_column := paste0(scolumn_1,'--',scolumn_2)]
  DT[, c('values_1_num', 'values_2_num', 'scolumn_1', 'scolumn_2') := NULL]
  data.table::setnames(DT, 'unif_sort_column', myname)

  return(DT)
}


#' @title combineCellProximityGenes
#' @name combineCellProximityGenes
#' @description Combine CPG scores in a pairwise manner.
#' @param cpgObject cell proximity gene score object
#' @param selected_ints subset of selected cell-cell interactions (optional)
#' @param selected_genes subset of selected genes (optional)
#' @param specific_genes_1 specific geneset combo (need to position match specific_genes_2)
#' @param specific_genes_2 specific geneset combo (need to position match specific_genes_1)
#' @param min_cells minimum number of target cell type
#' @param min_int_cells minimum number of interacting cell type
#' @param min_fdr minimum adjusted p-value
#' @param min_spat_diff minimum absolute spatial expression difference
#' @param min_log2_fc minimum absolute log2 fold-change
#' @param do_parallel run calculations in parallel with mclapply
#' @param cores number of cores to use if do_parallel = TRUE
#' @param verbose verbose
#' @return cpgObject that contains the filtered differential gene scores
#' @export
#' @examples
#'     combineCellProximityGenes(gobject)
combineCellProximityGenes = function(cpgObject,
                                     selected_ints = NULL,
                                     selected_genes = NULL,
                                     specific_genes_1 = NULL,
                                     specific_genes_2 = NULL,
                                     min_cells = 5,
                                     min_int_cells = 3,
                                     min_fdr = 0.05,
                                     min_spat_diff = 0,
                                     min_log2_fc = 0.5,
                                     do_parallel = TRUE,
                                     cores = NA,
                                     verbose = T) {

  ## check validity
  if(!'cpgObject' %in% class(cpgObject)) {
    stop('\n cpgObject needs to be the output from findCellProximityGenes() or findCPG() \n')
  }
  CPGscore = copy(cpgObject[['CPGscores']])

  if(!is.null(selected_ints)) {
    CPGscore = CPGscore[unif_int %in% selected_ints]
  }

  all_ints = unique(CPGscore[['unif_int']])

  # parallel
  if(do_parallel == TRUE) {

    # set number of cores
    if(is.na(cores) | !is.numeric(cores)) {
      cores = parallel::detectCores() - 1
    }

    GTGresults = parallel::mclapply(mc.cores = cores, X = all_ints, FUN = function(x) {

      tempres =  combineCellProximityGenes_per_interaction(cpgObject = cpgObject,
                                                           sel_int = x,
                                                           selected_genes = selected_genes,
                                                           specific_genes_1 = specific_genes_1,
                                                           specific_genes_2 = specific_genes_2,
                                                           min_cells = min_cells,
                                                           min_int_cells = min_int_cells,
                                                           min_fdr = min_fdr,
                                                           min_spat_diff = min_spat_diff,
                                                           min_log2_fc = min_log2_fc)

    })

  } else {
    # for loop
    GTGresults = list()

    for(i in 1:length(all_ints)) {

      x = all_ints[[i]]

      if(verbose == TRUE) print(x)

      tempres =  combineCellProximityGenes_per_interaction(cpgObject = cpgObject,
                                                           sel_int = x,
                                                           selected_genes = selected_genes,
                                                           specific_genes_1 = specific_genes_1,
                                                           specific_genes_2 = specific_genes_2,
                                                           min_cells = min_cells,
                                                           min_int_cells = min_int_cells,
                                                           min_fdr = min_fdr,
                                                           min_spat_diff = min_spat_diff,
                                                           min_log2_fc = min_log2_fc)
      GTGresults[[i]] = tempres
    }

  }

  final_results = do.call('rbind', GTGresults)

  final_results[, gene1_gene2 := paste0(genes_1,'--',genes_2)]

  final_results = sort_combine_two_DT_columns(final_results,
                                              column1 = 'genes_1', column2 = 'genes_2',
                                              myname = 'unif_gene_gene')
  #return(final_results)

  combCpgObject = list(combCPGscores = final_results,
                       Giotto_info = list('values' = cpgObject[['Giotto_info']][['values']],
                                          'cluster' = cpgObject[['Giotto_info']][['cluster']],
                                          'spatial network' = cpgObject[['Giotto_info']][['spatial network']]),
                       test_info = list('test' = cpgObject[['test_info']][['test']],
                                        'p.adj' = cpgObject[['test_info']][['p.adj']],
                                        'min cells' = cpgObject[['test_info']][['min cells']],
                                        'min interacting cells' = cpgObject[['test_info']][['min interacting cells']],
                                        'exclude selected cells' = cpgObject[['test_info']][['exclude selected cells']],
                                        'perm' = cpgObject[['test_info']][['perm']]))
  class(combCpgObject) = append(class(combCpgObject), 'combCpgObject')
  return(combCpgObject)


}


#' @title combineCPG
#' @name combineCPG
#' @description Combine CPG scores in a pairwise manner.
#' @param cpgObject cell proximity gene score object
#' @param selected_ints subset of selected cell-cell interactions (optional)
#' @param selected_genes subset of selected genes (optional)
#' @param specific_genes_1 specific geneset combo (need to position match specific_genes_2)
#' @param specific_genes_2 specific geneset combo (need to position match specific_genes_1)
#' @param min_cells minimum number of target cell type
#' @param min_int_cells minimum number of interacting cell type
#' @param min_fdr minimum adjusted p-value
#' @param min_spat_diff minimum absolute spatial expression difference
#' @param min_log2_fc minimum absolute log2 fold-change
#' @param do_parallel run calculations in parallel with mclapply
#' @param cores number of cores to use if do_parallel = TRUE
#' @param verbose verbose
#' @return cpgObject that contains the filtered differential gene scores
#' @export
#' @examples
#'     combineCPG(gobject)
combineCPG = function(cpgObject,
                      selected_ints = NULL,
                      selected_genes = NULL,
                      specific_genes_1 = NULL,
                      specific_genes_2 = NULL,
                      min_cells = 5,
                      min_int_cells = 3,
                      min_fdr = 0.05,
                      min_spat_diff = 0,
                      min_log2_fc = 0.5,
                      do_parallel = TRUE,
                      cores = NA,
                      verbose = T) {


  combineCellProximityGenes(cpgObject = cpgObject,
                            selected_ints = selected_ints,
                            selected_genes = selected_genes,
                            specific_genes_1 = specific_genes_1,
                            specific_genes_2 = specific_genes_2,
                            min_cells = min_cells,
                            min_int_cells = min_int_cells,
                            min_fdr = min_fdr,
                            min_spat_diff = min_spat_diff,
                            min_log2_fc = min_log2_fc,
                            do_parallel = do_parallel,
                            cores = cores,
                            verbose = verbose)


}











# * ####
# OLD CPGscores ####

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
                                                     diff_test = c('t.test', 'wilcox'),
                                                     minimum_unique_cells = NA,
                                                     exclude_selected_cells_from_test = T) {


  fast_check1 = sub_spatial_network[, c(neighb_IDs, neighb_col), with = F]
  colnames(fast_check1) = c('cell_ID', 'cell_type')
  fast_check2 = sub_spatial_network[, c(source_IDs, source_col), with = F]
  colnames(fast_check2) = c('cell_ID', 'cell_type')
  fast_check = unique(rbind(fast_check1, fast_check2))
  test = fast_check[, .N, by = cell_type]

  if(!is.na(minimum_unique_cells) & is.numeric(minimum_unique_cells) & any(test$N < minimum_unique_cells) == TRUE) {
    # not enough cells, do not continue
    return(NULL)
  }

  ## get unique cell IDs for all cells that belong to the cell type
  int_cell_types = sort(unique(fast_check[['cell_type']]))
  cell_type_1 = int_cell_types[[1]]
  cell_type_1_ids = fast_check[cell_type == cell_type_1][['cell_ID']]
  if(length(int_cell_types) > 1) {
    cell_type_2 = int_cell_types[[2]]
    cell_type_2_ids = fast_check[cell_type == cell_type_2][['cell_ID']]
  } else {
    cell_type_2 = cell_type_1
    cell_type_2_ids = cell_type_1_ids
  }


  nr_cell_type_1 = length(cell_type_1_ids)
  nr_cell_type_2 = length(cell_type_2_ids)

  ## create average scores to calculate fold-changes
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


  ## global expression with or without selected cell IDs
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



  ## stats test and averages
  diff_test = match.arg(diff_test, c('t.test', 'wilcox'))

  if(nr_cell_type_1 == 1) {
    if(do_diff_test == T) {
      if(diff_test == 't.test') {
        cell_type_1_test =  stats::t.test(x = as.vector(cell_type_1_ids_matrix), y = as.vector(all_cell_type_1_ids_matrix),exact = F)$p.value
      } else if(diff_test == 'wilcox') {
        cell_type_1_test =  stats::wilcox.test(x = as.vector(cell_type_1_ids_matrix), y = as.vector(all_cell_type_1_ids_matrix),exact = F)$p.value
      }
    }
    average_all_cell_type_1_ids = rowMeans(all_cell_type_1_ids_matrix)
  } else {
    if(do_diff_test == T) {

      if(diff_test == 't.test') {
        cell_type_1_test = sapply(1:nrow(cell_type_1_ids_matrix), function(i) stats::t.test(as.vector(cell_type_1_ids_matrix[i,]), as.vector(all_cell_type_1_ids_matrix[i,]),exact = F)$p.value)
      } else if(diff_test == 'wilcox') {
        cell_type_1_test = sapply(1:nrow(cell_type_1_ids_matrix), function(i) stats::wilcox.test(as.vector(cell_type_1_ids_matrix[i,]), as.vector(all_cell_type_1_ids_matrix[i,]),exact = F)$p.value)
      }

    }
    average_all_cell_type_1_ids = rowMeans(all_cell_type_1_ids_matrix)
  }

  if(nr_cell_type_2 == 1) {
    if(do_diff_test == T){
      if(diff_test == 't.test') {
        cell_type_2_test =  stats::t.test(x = as.vector(cell_type_2_ids_matrix), y = as.vector(all_cell_type_2_ids_matrix),exact = F)$p.value
      } else if(diff_test == 'wilcox') {
        cell_type_2_test =  stats::wilcox.test(x = as.vector(cell_type_2_ids_matrix), y = as.vector(all_cell_type_2_ids_matrix),exact = F)$p.value
      }
    }
    average_all_cell_type_2_ids = rowMeans(all_cell_type_2_ids_matrix)
  } else {
    if(do_diff_test == T) {
      if(diff_test == 't.test') {
        cell_type_2_test = sapply(1:nrow(cell_type_2_ids_matrix), function(i) stats::t.test(as.vector(cell_type_2_ids_matrix[i,]), as.vector(all_cell_type_2_ids_matrix[i,]),exact = F)$p.value)
      } else if(diff_test == 'wilcox') {
        cell_type_2_test = sapply(1:nrow(cell_type_2_ids_matrix), function(i) stats::wilcox.test(as.vector(cell_type_2_ids_matrix[i,]), as.vector(all_cell_type_2_ids_matrix[i,]),exact = F)$p.value)
      }
    }
    average_all_cell_type_2_ids = rowMeans(all_cell_type_2_ids_matrix)
  }

  combined_all_1_2_ids = (average_all_cell_type_1_ids+average_all_cell_type_2_ids)/2


  gene_names = names(combined_1_2_ids)

  if(do_diff_test == T) {
    res_DT = data.table::as.data.table(do.call('cbind', list(genes = gene_names,
                                                             cell_expr_1 = average_cell_type_1_ids, cell_expr_2 = average_cell_type_2_ids,
                                                             comb_expr = combined_1_2_ids,
                                                             all_cell_expr_1 = average_all_cell_type_1_ids, all_cell_expr_2 = average_all_cell_type_2_ids,
                                                             all_comb_expr = combined_all_1_2_ids,
                                                             pval_1 = cell_type_1_test, pval_2 = cell_type_2_test)))
  } else {
    res_DT = data.table::as.data.table(do.call('cbind', list(genes = gene_names,
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
                                            diff_test = c('t.test', 'wilcox'),
                                            minimum_unique_cells = NA,
                                            exclude_selected_cells_from_test = T,
                                            do_parallel = TRUE,
                                            cores = NA,
                                            verbose = T) {

  interactions = unique(as.character(spatial_network[[unified_int_col]]))

  # set number of cores
  if(is.na(cores) | !is.numeric(cores)) {
    cores = parallel::detectCores() - 1
  }


  ## parallel option
  if(do_parallel == TRUE) {
    result_list = parallel::mclapply(X = interactions, mc.cores = cores, FUN = function(selected_int) {

      sub_int_netw = spatial_network[get(unified_int_col) == selected_int]

      tempres = get_specific_interaction_gene_enrichment(sub_spatial_network = sub_int_netw,
                                                         source_col = source_col, source_IDs = source_IDs,
                                                         neighb_col = neighb_col, neighb_IDs = neighb_IDs,
                                                         expression_matrix = expression_matrix,
                                                         interaction_name = selected_int,
                                                         cell_annotation = cell_annotation,
                                                         annotation_ID = annotation_ID,
                                                         cell_type_col = cell_type_col,
                                                         do_diff_test = do_diff_test,
                                                         diff_test = diff_test,
                                                         minimum_unique_cells = minimum_unique_cells,
                                                         exclude_selected_cells_from_test = exclude_selected_cells_from_test)
    })
  } else {

    ## for-loop
    result_list <- list()
    for(selected_int in interactions) {

      if(verbose == TRUE){
        cat('start ', selected_int, '\n')
      }

      sub_int_netw = spatial_network[get(unified_int_col) == selected_int]

      tempres = get_specific_interaction_gene_enrichment(sub_spatial_network = sub_int_netw,
                                                         source_col = source_col, source_IDs = source_IDs,
                                                         neighb_col = neighb_col, neighb_IDs = neighb_IDs,
                                                         expression_matrix = expression_matrix,
                                                         interaction_name = selected_int,
                                                         cell_annotation = cell_annotation,
                                                         annotation_ID = annotation_ID,
                                                         cell_type_col = cell_type_col,
                                                         do_diff_test = do_diff_test,
                                                         diff_test = diff_test,
                                                         minimum_unique_cells = minimum_unique_cells,
                                                         exclude_selected_cells_from_test = exclude_selected_cells_from_test)
      result_list[[selected_int]] <- tempres
    }


  }

  # combine results
  final_res = do.call('rbind', result_list)

  if(verbose == TRUE){
    cat('\n interaction calculation is finished  \n')
  }

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
  first_combn = data.table::as.data.table(t(combn(x = all_cell_types, m = 2)))
  sec_combn = data.table::data.table(V1 = first_combn$V2, V2 = first_combn$V1)
  self_comb = data.table::data.table(V1 = all_cell_types, V2 = all_cell_types)

  # create data.table
  mult_comb = do.call('rbind', list(first_combn, sec_combn, self_comb))
  mult_comb[, c('V1', 'V2') := list(as.character(V1), as.character(V2))]
  mult_comb[, given_name := paste(c(V1,V2), collapse = '--'), by = 1:nrow(mult_comb)]
  mult_comb[, uniq_name := paste(sort(c(V1,V2)), collapse = '--'), by = 1:nrow(mult_comb)]

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
#' @param selected_genes selection of genes to perform calculations for
#' @param expression_values expression values to use
#' @param do_diff_test perform differential test
#' @param diff_test which differential expression test
#' @param false_discovery_test test to adjust p-values for multiple hypothesis testing
#' @param false_discovery_target adjust p-values per cell-cell pair or per gene
#' @param minimum_unique_cells minimum number of cells needed to proceed
#' @param fold_change_addendum constant to add when calculating log2 fold-change
#' @param in_two_directions shows enrichment in both directions: cell1-cell2, cell2-cell1
#' @param exclude_selected_cells_from_test exclude certain cells from test
#' @param do_parallel run enrichment calculations in parallel with mclapply
#' @param cores number of cores to use if do_parallel = TRUE
#' @param verbose verbose
#' @return Cell Proximity Gene scores (CPGscores) in data.table format
#' @details Function to calculate if genes are differentially expressed in cell types
#'  when they interact (according to physical proximity) with other cell types.
#'  The results data.table contains the following columns:
#' \itemize{
#'  \item{genes:}{ All or selected list of tested genes}
#'  \item{cell_expr_1:}{ average gene expression in cell type 1 from unified_int cell-cell interaction}
#'  \item{cell_expr_2:}{ average gene expression in cell type 2 from unified_int cell-cell interaction}
#'  \item{comb_expr:}{ combined average gene expression in cell type 1 and 2 from unified_int cell-cell interaction}
#'  \item{all_cell_expr_1:}{ average gene expression for all cells from cell type 1}
#'  \item{all_cell_expr_2:}{ average gene expression for all cells from cell type 2}
#'  \item{all_comb_expr:}{ combined average gene expression for all cells from cell type 1 and 2}
#'  \item{pval_1:}{ p-value from test between interacting cells and all cells from cell type 1}
#'  \item{pval_2:}{ p-value from test between interacting cells and all cells from cell type 2}
#'  \item{cell_type_1:}{ first cell type of cell-cell interaction}
#'  \item{cell_type_2:}{ second cell type of cell-cell interaction}
#'  \item{interaction:}{ the cell-cell interaction, based on physical proximity}
#'  \item{nr_1:}{ number of cell type 1 in the unified cell-cell interaction}
#'  \item{nr_2:}{ number of cell type 2 in the unified cell-cell interaction}
#'  \item{all_nr_1:}{ number of all cell type 1 in the whole dataset}
#'  \item{all_nr_2:}{ number of all cell type 2 in the whole dataset}
#'  \item{diff_spat:}{ difference between comb_expr and all_comb_expr}
#'  \item{diff_spat_1:}{ difference between cell_expr_1 and all_cell_expr_1}
#'  \item{diff_spat_2:}{ difference between cell_expr_1 and all_cell_expr_1}
#'  \item{log2fc_spat_1:}{ fold-change of diff_spat_1}
#'  \item{log2fc_spat_2:}{ fold-change of diff_spat_2}
#'  \item{log2fc_spat:}{ fold-change of diff_spat}
#'  \item{type_int:}{ type of interaction}
#'  \item{unified_int:}{ interaction with alphabetically sorted cell type 1 and cell type 2}
#'  \item{unif_int_rank:}{ 1 or 2}
#'  \item{fdr_1:}{ fdr from test between interacting cells and all cells from cell type 1}
#'  \item{fdr_2:}{ fdr from test between interacting cells and all cells from cell type 2}
#' }
#' @export
#' @examples
#'     getCellProximityGeneScores(gobject)
getCellProximityGeneScores = function(gobject,
                                      spatial_network_name = 'spatial_network',
                                      cluster_column = 'louvain_clus.1',
                                      selected_genes = NULL,
                                      expression_values = c('normalized', 'scaled', 'custom'),
                                      do_diff_test = TRUE,
                                      diff_test = c('t.test', 'wilcox'),
                                      false_discovery_test = c("holm", "hochberg", "hommel", "bonferroni",
                                                               "BH", "BY", "fdr", "none"),
                                      false_discovery_target = c('cell_interactions', 'genes'),
                                      minimum_unique_cells = NA,
                                      fold_change_addendum = 0.1,
                                      in_two_directions = TRUE,
                                      exclude_selected_cells_from_test = F,
                                      do_parallel = TRUE,
                                      cores = NA,
                                      verbose = T) {


  # 1. create annotated spatial network
  annot_spatial_network =  annotateSpatialNetwork(gobject = gobject,
                                                  spatial_network_name = spatial_network_name,
                                                  cluster_column = cluster_column)

  # 2. get expression values
  values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))
  expr_values = Giotto:::select_expression_values(gobject = gobject, values = values)

  # subset expression values
  if(!is.null(selected_genes)) {
    non_detected_genes = selected_genes[!selected_genes %in% rownames(expr_values)]
    if(length(non_detected_genes) > 0) {
      cat('\t the following genes were not detected ', non_detected_genes, '\t')
    }
    expr_values = expr_values[rownames(expr_values) %in% selected_genes, ]
  }

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
                                                                     annotation_ID = 'cell_ID',
                                                                     cell_type_col = cluster_column,
                                                                     do_diff_test = do_diff_test,
                                                                     minimum_unique_cells = minimum_unique_cells,
                                                                     diff_test = diff_test,
                                                                     exclude_selected_cells_from_test = exclude_selected_cells_from_test,
                                                                     do_parallel = do_parallel,
                                                                     cores = cores,
                                                                     verbose = verbose)

  # difference with spatial
  interaction_gene_scores[, diff_spat := as.numeric(comb_expr)-as.numeric(all_comb_expr)]
  interaction_gene_scores[, diff_spat_1 := as.numeric(cell_expr_1)-as.numeric(all_cell_expr_1)]
  interaction_gene_scores[, diff_spat_2 := as.numeric(cell_expr_2)-as.numeric(all_cell_expr_2)]
  data.table::setorder(interaction_gene_scores, diff_spat)

  # fold-change with spatial
  interaction_gene_scores[, log2fc_spat_1 := log2((as.numeric(cell_expr_1)+fold_change_addendum)/(as.numeric(all_cell_expr_1)+fold_change_addendum))]
  interaction_gene_scores[, log2fc_spat_2 := log2((as.numeric(cell_expr_2)+fold_change_addendum)/(as.numeric(all_cell_expr_2)+fold_change_addendum))]
  interaction_gene_scores[, log2fc_spat := log2((as.numeric(comb_expr)+fold_change_addendum)/(as.numeric(all_comb_expr)+fold_change_addendum))]



  # expand scores to all possible cell-cell combinations
  # e.g. both astrocyte-NSC and NSC-astrocyte
  if(in_two_directions == TRUE & do_diff_test == TRUE) {
    CPGscore_first_direction  = interaction_gene_scores
    CPGscore_second_direction = interaction_gene_scores
    colnames(CPGscore_second_direction) = c('genes', 'cell_expr_2', 'cell_expr_1', 'comb_expr',
                                            'all_cell_expr_2', 'all_cell_expr_1', 'all_comb_expr',
                                            'pval_2', 'pval_1',
                                            'cell_type_2', 'cell_type_1', 'interaction',
                                            'nr_2', 'nr_1', 'all_nr_2', 'all_nr_1',
                                            'diff_spat', 'diff_spat_2', 'diff_spat_1',
                                            'log2fc_spat_2', 'log2fc_spat_1', 'log2fc_spat')
    CPGscore_second_direction = CPGscore_second_direction[, colnames(CPGscore_first_direction), with = F]
    CPGscore_second_direction[, interaction := paste0(cell_type_1,'--', cell_type_2)]
    CPGscore = rbind(CPGscore_first_direction, CPGscore_second_direction)
    CPGscore = unique(CPGscore)


    # homo or hetero-typic interaction
    CPGscore[, type_int := ifelse(cell_type_1 == cell_type_2, 'homo', 'hetero')]

    # create unique cell-cell interaction names
    name_conversion = Giotto:::get_cell_to_cell_sorted_name_conversion(all_cell_types = unique(c(CPGscore$cell_type_1,
                                                                                                 CPGscore$cell_type_2)))
    CPGscore[, unified_int := name_conversion[[interaction]], by = 1:nrow(CPGscore)]
    CPGscore[, unif_int_rank := ifelse(interaction == unified_int, 1, 2)]



  } else {

    CPGscore = interaction_gene_scores
    CPGscore[, type_int := ifelse(cell_type_1 == cell_type_2, 'homo', 'hetero')]

  }

  ## add multiple hypothesis testing and set colomn classes manually
  if(do_diff_test == TRUE) {

    false_discovery_test = match.arg(false_discovery_test,  c("holm", "hochberg", "hommel", "bonferroni",
                                                              "BH", "BY", "fdr", "none"))

    false_discovery_target = match.arg(false_discovery_target, c('cell_interactions', 'genes'))
    false_discovery_other = ifelse(false_discovery_target == 'cell_interactions', 'genes', 'unified_int')

    # don't calculate fdr in two ways
    # set N to include missing values
    rank1_scores = unique(CPGscore[unif_int_rank == 1])
    rank1_scores[,fdr_1 := p.adjust(p = pval_1, method = false_discovery_test, n = .N), by = false_discovery_other]
    rank1_scores[,fdr_2 := p.adjust(p = pval_2, method = false_discovery_test, n = .N), by = false_discovery_other]


    part1 = rank1_scores[,.(genes, cell_type_1, cell_type_2, fdr_1, fdr_2)]
    part2 = rank1_scores[,.(genes, cell_type_2, cell_type_1, fdr_2, fdr_1)]
    setnames(part2, c('cell_type_2', 'cell_type_1', 'fdr_2', 'fdr_1'), c('cell_type_1', 'cell_type_2', 'fdr_1', 'fdr_2'))
    part_12 = unique(rbind(part1, part2))

    CPGscore_1 = data.table:::merge.data.table(CPGscore[unif_int_rank == 1], part_12[,.(genes, cell_type_1, cell_type_2, fdr_1, fdr_2)],
                                  by.x = c('genes', 'cell_type_1', 'cell_type_2'),
                                  by.y = c('genes', 'cell_type_1', 'cell_type_2'))
    CPGscore_2 = data.table:::merge.data.table(CPGscore[unif_int_rank == 2], part_12[,.(genes, cell_type_1, cell_type_2, fdr_1, fdr_2)],
                                  by.x = c('genes', 'cell_type_1', 'cell_type_2'),
                                  by.y = c('genes', 'cell_type_1', 'cell_type_2'))
    CPGscore = rbind(CPGscore_1, CPGscore_2)

    #return(list(CPGscore, rank1_scores))

    # bonferroni
    #total_comparisons_per_gene = length(unique(CPGscore[['unified_int']]))
    #CPGscore[, fdr_1_a := as.numeric(pval_1) * total_comparisons_per_gene]
    #CPGscore[, fdr_2_a := as.numeric(pval_2) * total_comparisons_per_gene]

    #CPGscore[, fdr_1_a := ifelse(fdr_1_a > 1, 1, fdr_1_a)]
    #CPGscore[, fdr_2_a := ifelse(fdr_2_a > 1, 1, fdr_2_a)]


    # set classes manually
    changeCols = c('cell_expr_1', 'cell_expr_2', 'comb_expr',
                   'all_cell_expr_1' , 'all_cell_expr_2',  'all_comb_expr',
                   'pval_1', 'pval_2', 'diff_spat', 'diff_spat_1', 'diff_spat_2',
                   'log2fc_spat_1', 'log2fc_spat_2', 'log2fc_spat',
                   'fdr_1', 'fdr_2')
    CPGscore[,(changeCols):= lapply(.SD, as.numeric), .SDcols = changeCols]

  } else {
    # set classes manually
    changeCols = c('cell_expr_1', 'cell_expr_2', 'comb_expr',
                   'all_cell_expr_1' , 'all_cell_expr_2',  'all_comb_expr',
                   'diff_spat', 'diff_spat_1', 'diff_spat_2',
                   'log2fc_spat_1', 'log2fc_spat_2', 'log2fc_spat')
    CPGscore[,(changeCols):= lapply(.SD, as.numeric), .SDcols = changeCols]
  }
  changeCols = c('nr_1', 'nr_2', 'all_nr_1', 'all_nr_2')
  CPGscore[,(changeCols):= lapply(.SD, as.integer), .SDcols = changeCols]

  return(CPGscore)

}



# * ####
# OLD GTGscores ####

#' @title combine_ints_f
#' @name combine_ints_f
#' @description function to combine gene enrichment interactions
#' @param cell_int selected cell interaction
#' @param all_ints all interactions
#' @param unif_gene_scores unif_gene_scores results
#' @param specific_genes_1 specific source genes (see details)
#' @param specific_genes_2 specific target genes (see details)
#' @param min_cells min number of cells threshold
#' @param min_pval p-value threshold
#' @param min_spat_diff spatial difference threshold
#' @param min_log2_fc log2 fold-change threshold
#' @return Gene to gene scores in data.table format
combine_ints_f = function(cell_int,
                          all_ints,
                          unif_gene_scores,
                          specific_genes_1 = NULL,
                          specific_genes_2 = NULL,
                          min_cells = 5,
                          min_fdr = 0.05,
                          min_spat_diff = 0.2,
                          min_log2_fc = 0.5) {


  sel_int = all_ints[cell_int]

  # cell types
  sel_type_1 = strsplit(x = sel_int, split = '--')[[1]][1]
  sel_type_2 = strsplit(x = sel_int, split = '--')[[1]][2]

  # subsets of the data
  subset_type_1 = unif_gene_scores[cell_type_1 == sel_type_1][(fdr_1 <= min_fdr & nr_1 >= min_cells & abs(diff_spat_1) >= min_spat_diff & abs(log2fc_spat_1) >= min_log2_fc)][cell_type_2 == sel_type_2]
  subset_type_2 = unif_gene_scores[cell_type_2 == sel_type_2][(fdr_2 <= min_fdr & nr_2 >= min_cells & abs(diff_spat_2) >= min_spat_diff & abs(log2fc_spat_2) >= min_log2_fc)][cell_type_1 == sel_type_1]

  subset_type_1 = subset_type_1[,.(genes, cell_type_1, nr_1, cell_expr_1, all_nr_1, all_cell_expr_1, fdr_1, diff_spat_1, log2fc_spat_1, unified_int)]
  setnames(subset_type_1, 'genes', 'genes_1')
  subset_type_2 = subset_type_2[,.(genes, cell_type_2, nr_2, cell_expr_2, all_nr_2, all_cell_expr_2, fdr_2, diff_spat_2, log2fc_spat_2, unified_int)]
  setnames(subset_type_2, 'genes', 'genes_2')

  # merge data again
  mergetest = merge(subset_type_1, subset_type_2, by = 'unified_int', allow.cartesian = T)

  if(nrow(mergetest) > 0) {


    # create additional columns
    mergetest[, gene_gene := paste0(genes_1,'--',genes_2), by = 1:nrow(mergetest)]

    # only keep specific combinations
    if((!is.null(specific_genes_1) & !is.null(specific_genes_2))) {

      LR_combo = paste0(specific_genes_1,'--', specific_genes_2)
      RL_combo = paste0(specific_genes_2,'--', specific_genes_1)
      mergetest = mergetest[gene_gene %in% c(LR_combo, RL_combo)]

    }

    if(nrow(mergetest) > 0) {
      mergetest[, unif_gene_gene := paste(sort(c(genes_1, genes_2)), collapse = '--'), by = 1:nrow(mergetest)]
    }
  }

  return(mergetest)
}



#' @title getGeneToGeneScores
#' @name getGeneToGeneScores
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
#' @param do_parallel run enrichment calculations in parallel with mclapply
#' @param cores number of cores to use if do_parallel = TRUE
#' @param verbose verbose
#' @return Gene to gene scores in data.table format
#' @details This converts the single gene cell proximityscores into pairwise combinations
#' of genes, which allows you to determine if 2 genes are differentially expressed in interacting
#' cell types.
#' @export
#' @examples
#'     getGeneToGeneScores(CPGscore)
getGeneToGeneScores <- function(CPGscore,
                                selected_genes = NULL,
                                specific_genes_1 = NULL,
                                specific_genes_2 = NULL,
                                min_cells = 5,
                                min_fdr = 0.05,
                                min_spat_diff = 0.2,
                                min_log2_fc = 0.5,
                                direction = c('both', 'up', 'down'),
                                fold_change_addendum = 0.1,
                                do_parallel = TRUE,
                                cores = NA,
                                verbose = TRUE) {


  direction = match.arg(direction, choices = c('both', 'up', 'down'))

  # remove redundant data
  unif_gene_scores = CPGscore[unif_int_rank == 1]

  # make sure fdr's are numeric
  unif_gene_scores[, fdr_1 := as.numeric(fdr_1)]
  unif_gene_scores[, fdr_2 := as.numeric(fdr_2)]

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
  all_ints_selection = unif_gene_scores[(fdr_1 <= min_fdr & nr_1 >= min_cells & abs(diff_spat_1) >= min_spat_diff & abs(log2fc_spat_1) >= min_log2_fc) |
                                          (fdr_2 <= min_fdr & nr_2 >= min_cells & abs(diff_spat_2) >= min_spat_diff  & abs(log2fc_spat_2) >= min_log2_fc)]

  unif_gene_scores = all_ints_selection
  all_ints = unique(all_ints_selection[['unified_int']])

  # parallel option
  if(do_parallel == TRUE) {

    # set number of cores
    if(is.na(cores) | !is.numeric(cores)) {
      cores = parallel::detectCores() - 1
    }

    savelist = parallel::mclapply(X = 1:length(all_ints), mc.cores = cores, FUN = function(x) {

      combine_ints_f(cell_int = x,
                     all_ints = all_ints,
                     unif_gene_scores = unif_gene_scores,
                     min_cells = min_cells,
                     min_fdr = min_fdr,
                     min_spat_diff = min_spat_diff,
                     min_log2_fc = min_log2_fc)

    })




  } else {

    # for loop
    savelist = list()

    for(cell_int in 1:length(all_ints)) {

      sel_int = all_ints[cell_int]

      if(verbose == TRUE) {
        cat('\n', sel_int, '\n')
      }

      # cell types
      sel_type_1 = strsplit(x = sel_int, split = '--')[[1]][1]
      sel_type_2 = strsplit(x = sel_int, split = '--')[[1]][2]

      # subsets of the data
      subset_type_1 = unif_gene_scores[cell_type_1 == sel_type_1][(fdr_1 <= min_fdr & nr_1 >= min_cells & abs(diff_spat_1) >= min_spat_diff & abs(log2fc_spat_1) >= min_log2_fc)][cell_type_2 == sel_type_2]
      subset_type_2 = unif_gene_scores[cell_type_2 == sel_type_2][(fdr_2 <= min_fdr & nr_2 >= min_cells & abs(diff_spat_2) >= min_spat_diff & abs(log2fc_spat_2) >= min_log2_fc)][cell_type_1 == sel_type_1]

      subset_type_1 = subset_type_1[,.(genes, cell_type_1, nr_1, cell_expr_1, all_nr_1, all_cell_expr_1, fdr_1, diff_spat_1, log2fc_spat_1, unified_int)]
      setnames(subset_type_1, 'genes', 'genes_1')
      subset_type_2 = subset_type_2[,.(genes, cell_type_2, nr_2, cell_expr_2, all_nr_2, all_cell_expr_2, fdr_2, diff_spat_2, log2fc_spat_2, unified_int)]
      setnames(subset_type_2, 'genes', 'genes_2')

      # merge data again
      mergetest = merge(subset_type_1, subset_type_2, by = 'unified_int', allow.cartesian = T)

      if(nrow(mergetest) > 0) {


        # create additional columns
        mergetest[, gene_gene := paste0(genes_1,'--',genes_2), by = 1:nrow(mergetest)]

        # only keep specific combinations
        if((!is.null(specific_genes_1) & !is.null(specific_genes_2))) {

          LR_combo = paste0(specific_genes_1,'--', specific_genes_2)
          RL_combo = paste0(specific_genes_2,'--', specific_genes_1)
          mergetest = mergetest[gene_gene %in% c(LR_combo, RL_combo)]

        }

        if(nrow(mergetest) > 0) {
          mergetest[, unif_gene_gene := paste(sort(c(genes_1, genes_2)), collapse = '--'), by = 1:nrow(mergetest)]
          savelist[[cell_int]] = mergetest
        }
      }

    }

  }





  finalres = do.call('rbind', savelist)

  # remove redundant homo-typic interaction data
  if(verbose == TRUE) cat('\n calculate additional information \n')

  finalres[, type_int := ifelse(cell_type_1 == cell_type_2, 'homo', 'hetero')]
  finalres[, unif_gene_gene_int := 1:.N, by = .(unified_int, unif_gene_gene)]
  finalres = finalres[type_int == 'hetero' | (type_int == 'homo' & unif_gene_gene_int == 1)]


  changeCols = c('cell_expr_1', 'cell_expr_2',
                 'all_cell_expr_1' , 'all_cell_expr_2',
                 'diff_spat_1', 'diff_spat_2',
                 'log2fc_spat_1', 'log2fc_spat_2',
                 'fdr_1', 'fdr_2')
  finalres[,(changeCols):= lapply(.SD, as.numeric), .SDcols = changeCols]

  if(direction == 'both') {
    finalres = finalres[(fdr_1 <= min_fdr & nr_1 >= min_cells & abs(diff_spat_1) >= min_spat_diff & abs(log2fc_spat_1) >= min_log2_fc) &
                          (fdr_2 <= min_fdr & nr_2 >= min_cells & abs(diff_spat_2) >= min_spat_diff  & abs(log2fc_spat_2) >= min_log2_fc)]
  } else if(direction == 'up') {
    finalres = finalres[(fdr_1 <= min_fdr & nr_1 >= min_cells & diff_spat_1 >= min_spat_diff & log2fc_spat_1 >= min_log2_fc) &
                          (fdr_2 <= min_fdr & nr_2 >= min_cells & diff_spat_2 >= min_spat_diff  & log2fc_spat_2 >= min_log2_fc)]
  } else if(direction == 'down') {
    finalres = finalres[(fdr_1 <= min_fdr & nr_1 >= min_cells & diff_spat_1 <= -min_spat_diff & log2fc_spat_1 <= -min_log2_fc) &
                          (fdr_2 <= min_fdr & nr_2 >= min_cells & diff_spat_2 <= -min_spat_diff  & log2fc_spat_2 <= -min_log2_fc)]
  }


  # provide additional data
  finalres[, all_cell_expr := as.numeric(all_cell_expr_1)+as.numeric(all_cell_expr_2), 1:nrow(finalres)]
  finalres[, spatial_cell_expr := as.numeric(cell_expr_1)+as.numeric(cell_expr_2), 1:nrow(finalres)]
  finalres[, diff_spatial_expr := spatial_cell_expr-all_cell_expr]
  finalres[, log2fc_spatial_expr := log2((spatial_cell_expr+fold_change_addendum)/(all_cell_expr+fold_change_addendum))]
  finalres[, all_cell_rank := rank(-all_cell_expr), by = unif_gene_gene]
  finalres[, spatial_cell_rank := rank(-spatial_cell_expr), by = unif_gene_gene]

  change_values = unlist(apply(finalres, MARGIN = 1, FUN = function(x) {
    Giotto:::direction_test(x, min_fdr = min_fdr)
  }))
  finalres[, change := change_values]

  data.table::setorder(finalres, -log2fc_spatial_expr)

  return(finalres)


}





# * ####
# cell communication ####


#' @title average_gene_gene_expression_in_groups
#' @name average_gene_gene_expression_in_groups
#' @description calculate average expression per cluster
#' @param gobject giotto object to use
#' @param cluster_column cluster column with cell type information
#' @param gene_set_1 first specific gene set from gene pairs
#' @param gene_set_2 second specific gene set from gene pairs
#' @return data.table with average expression scores for each cluster
#' @details Details will follow soon.
#' @examples
#'     average_gene_gene_expression_in_groups(gobject)
average_gene_gene_expression_in_groups = function(gobject,
                                                  cluster_column = 'cell_types',
                                                  gene_set_1,
                                                  gene_set_2) {

  average_DT = create_average_DT(gobject = gobject, meta_data_name = cluster_column)

  # change column names back to original
  new_colnames = gsub(pattern = 'cluster_', replacement = '', colnames(average_DT))
  colnames(average_DT) = new_colnames

  # keep order of colnames
  colnames_order = new_colnames

  # gene_set_1 and gene_set_2 need to have same length and all genes need to be present in data
  if(length(gene_set_1) != length(gene_set_2)) {
    stop('\n length of set1 needs to be the same as that of set2 \n')
  }

  if(!all(c(gene_set_1, gene_set_2) %in% rownames(average_DT) == T)) {
    stop('\n all selected genes from set 1 and 2 need to be present \n')
  }

  LR_pairs = paste0(gene_set_1,'-',gene_set_2)

  # get ligand and receptor information
  ligand_match = average_DT[match(gene_set_1, rownames(average_DT)), ,drop = F]
  receptor_match = average_DT[match(gene_set_2, rownames(average_DT)), ,drop = F]

  all_ligand_cols = colnames(ligand_match)
  lig_test = data.table::as.data.table(reshape2::melt(ligand_match, measure.vars = all_ligand_cols))
  lig_test[, ligand := rep(rownames(ligand_match), ncol(ligand_match))]
  lig_test[, ligand := strsplit(ligand,'\\.')[[1]][1] , by = 1:nrow(lig_test)]
  lig_test[, LR_comb := rep(LR_pairs, ncol(ligand_match))]
  setnames(lig_test, 'value', 'lig_expr')
  setnames(lig_test, 'variable', 'lig_cell_type')

  all_receptor_cols = colnames(receptor_match)
  rec_test = data.table::as.data.table(reshape2::melt(receptor_match, measure.vars = all_receptor_cols))
  rec_test[, receptor := rep(rownames(receptor_match), ncol(receptor_match))]
  rec_test[, receptor := strsplit(receptor,'\\.')[[1]][1] , by = 1:nrow(rec_test)]
  rec_test[, LR_comb := rep(LR_pairs, ncol(receptor_match))]
  setnames(rec_test, 'value', 'rec_expr')
  setnames(rec_test, 'variable', 'rec_cell_type')

  lig_rec_test = merge(lig_test, rec_test, by = 'LR_comb', allow.cartesian = T)
  lig_rec_test[, LR_expr := lig_expr+rec_expr]


  lig_rec_test[, lig_cell_type := factor(lig_cell_type, levels = colnames_order)]
  lig_rec_test[, rec_cell_type := factor(rec_cell_type, levels = colnames_order)]
  setorder(lig_rec_test, LR_comb, lig_cell_type, rec_cell_type)

  return(lig_rec_test)

}


#' @title exprCellCellcom
#' @name exprCellCellcom
#' @description Cell-Cell communication scores based on expression only
#' @param gobject giotto object to use
#' @param cluster_column cluster column with cell type information
#' @param random_iter number of iterations
#' @param gene_set_1 first specific gene set from gene pairs
#' @param gene_set_2 second specific gene set from gene pairs
#' @param log2FC_addendum addendum to add when calculating log2FC
#' @param adjust_method which method to adjust p-values
#' @param adjust_target adjust multiple hypotheses at the cell or gene level
#' @param verbose verbose
#' @return Cell-Cell communication scores for gene pairs based on expression only
#' @details Statistical framework to identify if pairs of genes (such as ligand-receptor combinations)
#' are expressed at higher levels than expected based on a reshuffled null distribution of gene expression values,
#' without considering the spatial position of cells.
#' More details will follow soon.
#' @export
#' @examples
#'     exprCellCellcom(gobject)
exprCellCellcom = function(gobject,
                           cluster_column = 'cell_types',
                           random_iter = 100,
                           gene_set_1,
                           gene_set_2,
                           log2FC_addendum = 0.1,
                           adjust_method = c("fdr", "bonferroni","BH", "holm", "hochberg", "hommel",
                                             "BY", "none"),
                           adjust_target = c('genes', 'cells'),
                           verbose = T) {


  # get parameters
  adjust_method = match.arg(adjust_method, choices = c("fdr", "bonferroni","BH", "holm", "hochberg", "hommel",
                                                       "BY", "none"))
  adjust_target = match.arg(adjust_target, choices = c('genes', 'cells'))

    # get information about number of cells
  cell_metadata = pDataDT(gobject)
  nr_cell_types = cell_metadata[,.N, by = c(cluster_column)]
  nr_cells = nr_cell_types$N
  names(nr_cells) = nr_cell_types$cluster_column


  comScore = average_gene_gene_expression_in_groups(gobject = gobject,
                                                    cluster_column = cluster_column,
                                                    gene_set_1 = gene_set_1, gene_set_2 = gene_set_2)

  comScore[, lig_nr := nr_cells[lig_cell_type]]
  comScore[, rec_nr := nr_cells[rec_cell_type]]

  # prepare for randomized scores
  total_av = rep(0, nrow(comScore))
  total_sum = rep(0, nrow(comScore))
  total_bool = rep(0, nrow(comScore))


  ## parallel option ##
  # not yet available


  for(sim in 1:random_iter) {

    if(verbose == TRUE) cat('simulation ', sim, '\n')


    # create temporary giotto
    tempGiotto = subsetGiotto(gobject = gobject)

    # randomize annoation
    cell_types = cell_metadata[[cluster_column]]
    random_cell_types = sample(x = cell_types, size = length(cell_types))
    tempGiotto = addCellMetadata(gobject = tempGiotto, new_metadata = random_cell_types, by_column = F)

    # get random communication scores
    randomScore = average_gene_gene_expression_in_groups(gobject = tempGiotto,
                                                         cluster_column = 'random_cell_types',
                                                         gene_set_1 = gene_set_1, gene_set_2 = gene_set_2)

    # average random score
    total_av = total_av + randomScore[['LR_expr']]

    # difference between observed and random
    difference = comScore[['LR_expr']] - randomScore[['LR_expr']]

    # calculate total difference
    total_sum = total_sum+difference

    # calculate p-values
    difference[difference > 0] = 1
    difference[difference < 0] = -1
    total_bool = total_bool + difference

  }

  comScore[, rand_expr := total_av/random_iter]
  comScore[, av_diff := total_sum/random_iter]
  comScore[, log2fc := log2((LR_expr+log2FC_addendum)/(rand_expr+log2FC_addendum))]
  comScore[, pvalue := total_bool/random_iter]
  comScore[, pvalue := ifelse(pvalue > 0, 1-pvalue, 1+pvalue)]
  comScore[, LR_cell_comb := paste0(lig_cell_type,'--',rec_cell_type)]

  if(adjust_target == 'genes') {
    comScore[, p.adj := p.adjust(pvalue, method = adjust_method), by = .(LR_cell_comb)]
  } else if(adjust_target == 'cells'){
    comScore[, p.adj := p.adjust(pvalue, method = adjust_method), by = .(LR_comb)]
  }
  comScore[, PI := log2fc*(1-p.adj)]
  data.table::setorder(comScore, LR_comb, -LR_expr)

  return(comScore)

}



#' @title create_cell_type_random_cell_IDs
#' @name create_cell_type_random_cell_IDs
#' @description creates randomized cell ids within a selection of cell types
#' @param gobject giotto object to use
#' @param cluster_column cluster column with cell type information
#' @param needed_cell_types vector of cell type names for which a random id will be found
#' @return list of randomly sampled cell ids with same cell type composition
#' @details Details will follow.
#' @examples
#'     create_cell_type_random_cell_IDs(gobject)
create_cell_type_random_cell_IDs = function(gobject,
                                            cluster_column = 'cell_types',
                                            needed_cell_types) {

  # subset metadata to choose from
  full_metadata = pDataDT(gobject)
  possible_metadata = full_metadata[get(cluster_column) %in% unique(needed_cell_types)]

  sample_ids = list()

  uniq_types = unique(needed_cell_types)

  for(i in 1:length(uniq_types)) {

    uniq_type = uniq_types[i]
    length_random = length(needed_cell_types[needed_cell_types == uniq_type])
    sub_sample_ids = possible_metadata[get(cluster_column) == uniq_type][sample(x = 1:.N, size = length_random)][['cell_ID']]
    sample_ids[[i]] = sub_sample_ids

  }
  return(unlist(sample_ids))
}




#' @title specificCellCellcommunicationScores
#' @name specificCellCellcommunicationScores
#' @description Specific Cell-Cell communication scores based on spatial expression of interacting cells
#' @param gobject giotto object to use
#' @param spatial_network_name spatial network to use for identifying interacting cells
#' @param cluster_column cluster column with cell type information
#' @param random_iter number of iterations
#' @param cell_type_1 first cell type
#' @param cell_type_2 second cell type
#' @param gene_set_1 first specific gene set from gene pairs
#' @param gene_set_2 second specific gene set from gene pairs
#' @param log2FC_addendum addendum to add when calculating log2FC
#' @param min_observations minimum number of interactions needed to be considered
#' @param adjust_method which method to adjust p-values
#' @param adjust_target adjust multiple hypotheses at the cell or gene level
#' @param verbose verbose
#' @return Cell-Cell communication scores for gene pairs based on spatial interaction
#' @details Statistical framework to identify if pairs of genes (such as ligand-receptor combinations)
#' are expressed at higher levels than expected based on a reshuffled null distribution
#' of gene expression values in cells that are spatially in proximity to eachother..
#' More details will follow soon.
#' @export
#' @examples
#'     specificCellCellcommunicationScores(gobject)
specificCellCellcommunicationScores = function(gobject,
                                               spatial_network_name = 'spatial_network',
                                               cluster_column = 'cell_types',
                                               random_iter = 100,
                                               cell_type_1 = 'astrocyte',
                                               cell_type_2 = 'endothelial',
                                               gene_set_1,
                                               gene_set_2,
                                               log2FC_addendum = 0.1,
                                               min_observations = 2,
                                               adjust_method = c("fdr", "bonferroni","BH", "holm", "hochberg", "hommel",
                                                                 "BY", "none"),
                                               adjust_target = c('genes', 'cells'),
                                               verbose = T) {

  # get parameters
  adjust_method = match.arg(adjust_method, choices = c("fdr", "bonferroni","BH", "holm", "hochberg", "hommel",
                                                       "BY", "none"))
  adjust_target = match.arg(adjust_target, choices = c('genes', 'cells'))

  # metadata
  cell_metadata = pDataDT(gobject = gobject)

  # get annotated spatial network
  annot_network = annotateSpatialNetwork(gobject, spatial_network_name = spatial_network_name, cluster_column = cluster_column)

  cell_direction_1 = paste0(cell_type_1,'-',cell_type_2)
  cell_direction_2 = paste0(cell_type_2,'-',cell_type_1)

  subset_annot_network = annot_network[from_to %in% c(cell_direction_1, cell_direction_1)]

  # make sure that there are sufficient observations
  if(nrow(subset_annot_network) <= min_observations) {

    return(NULL)

  } else {


    # subset giotto object to only interacting cells
    subset_ids = unique(c(subset_annot_network$to, subset_annot_network$from))
    subsetGiotto = subsetGiotto(gobject = gobject, cell_ids = subset_ids)

    # get information about number of cells
    temp_meta = pDataDT(subsetGiotto)
    nr_cell_types = temp_meta[cell_ID %in% subset_ids][,.N, by = c(cluster_column)]
    nr_cells = nr_cell_types$N
    names(nr_cells) = nr_cell_types$cell_types

    # get average communication scores
    comScore = average_gene_gene_expression_in_groups(gobject = subsetGiotto,
                                                      cluster_column = cluster_column,
                                                      gene_set_1 = gene_set_1, gene_set_2 = gene_set_2)
    comScore = comScore[(lig_cell_type == cell_type_1 & rec_cell_type == cell_type_2) |
                          (lig_cell_type == cell_type_2 & rec_cell_type == cell_type_1)]

    comScore[, lig_nr := nr_cells[lig_cell_type]]
    comScore[, rec_nr := nr_cells[rec_cell_type]]

    # prepare for randomized scores
    total_av = rep(0, nrow(comScore))
    total_sum = rep(0, nrow(comScore))
    total_bool = rep(0, nrow(comScore))

    # identify which cell types you need
    subset_metadata = cell_metadata[cell_ID %in% subset_ids]
    needed_cell_types = subset_metadata[[cluster_column]]

    # simulations
    for(sim in 1:random_iter) {

      if(verbose == TRUE) cat('simulation ', sim, '\n')

      # get random ids and subset
      random_ids = create_cell_type_random_cell_IDs(gobject = gobject, cluster_column = cluster_column,
                                                    needed_cell_types = needed_cell_types)
      tempGiotto = subsetGiotto(gobject = gobject, cell_ids = random_ids)

      # get random communication scores
      randomScore = average_gene_gene_expression_in_groups(gobject = tempGiotto,
                                                           cluster_column = cluster_column,
                                                           gene_set_1 = gene_set_1, gene_set_2 = gene_set_2)
      randomScore = randomScore[(lig_cell_type == cell_type_1 & rec_cell_type == cell_type_2) |
                                  (lig_cell_type == cell_type_2 & rec_cell_type == cell_type_1)]


      # average random score
      total_av = total_av + randomScore[['LR_expr']]

      # difference between observed and random
      difference = comScore[['LR_expr']] - randomScore[['LR_expr']]

      # calculate total difference
      total_sum = total_sum+difference

      # calculate p-values
      difference[difference > 0] = 1
      difference[difference < 0] = -1
      total_bool = total_bool + difference

    }

    comScore[, rand_expr := total_av/random_iter]
    comScore[, av_diff := total_sum/random_iter]
    comScore[, log2fc := log2((LR_expr+log2FC_addendum)/(rand_expr+log2FC_addendum))]
    comScore[, pvalue := total_bool/random_iter]
    comScore[, pvalue := ifelse(pvalue > 0, 1-pvalue, 1+pvalue)]
    comScore[, LR_cell_comb := paste0(lig_cell_type,'--',rec_cell_type)]

    if(adjust_target == 'genes') {
      comScore[, p.adj := p.adjust(pvalue, method = adjust_method), by = .(LR_cell_comb)]
    } else if(adjust_target == 'cells'){
      comScore[, p.adj := p.adjust(pvalue, method = adjust_method), by = .(LR_comb)]
    }

    comScore[, PI := log2fc*(1-p.adj)]

    return(comScore)

  }
}


#' @title spatCellCellcom
#' @name spatCellCellcom
#' @description Spatial Cell-Cell communication scores based on spatial expression of interacting cells
#' @param gobject giotto object to use
#' @param spatial_network_name spatial network to use for identifying interacting cells
#' @param cluster_column cluster column with cell type information
#' @param random_iter number of iterations
#' @param gene_set_1 first specific gene set from gene pairs
#' @param gene_set_2 second specific gene set from gene pairs
#' @param log2FC_addendum addendum to add when calculating log2FC
#' @param min_observations minimum number of interactions needed to be considered
#' @param adjust_method which method to adjust p-values
#' @param adjust_target adjust multiple hypotheses at the cell or gene level
#' @param do_parallel run calculations in parallel with mclapply
#' @param cores number of cores to use if do_parallel = TRUE
#' @param verbose verbose
#' @return Cell-Cell communication scores for gene pairs based on spatial interaction
#' @details Statistical framework to identify if pairs of genes (such as ligand-receptor combinations)
#' are expressed at higher levels than expected based on a reshuffled null distribution
#' of gene expression values in cells that are spatially in proximity to eachother..
#' More details will follow soon.
#' @export
#' @examples
#'     spatCellCellcom(gobject)
spatCellCellcom = function(gobject,
                           spatial_network_name = 'spatial_network',
                           cluster_column = 'cell_types',
                           random_iter = 100,
                           gene_set_1,
                           gene_set_2,
                           log2FC_addendum = 0.1,
                           min_observations = 2,
                           adjust_method = c("fdr", "bonferroni","BH", "holm", "hochberg", "hommel",
                                             "BY", "none"),
                           adjust_target = c('genes', 'cells'),
                           do_parallel = TRUE,
                           cores = NA,
                           verbose = c('a little', 'a lot', 'none')) {

  verbose = match.arg(verbose, choices = c('a little', 'a lot', 'none'))

  cell_metadata = pDataDT(gobject)

  ## get all combinations between cell types
  all_uniq_values = unique(cell_metadata[[cluster_column]])
  same_DT = data.table(V1 = all_uniq_values, V2 = all_uniq_values)
  combn_DT = as.data.table(t(combn(all_uniq_values, m = 2)))
  combn_DT = rbind(same_DT, combn_DT)

  ## parallel option ##
  if(do_parallel == TRUE) {

    # set number of cores
    if(is.na(cores) | !is.numeric(cores)) {
      cores = parallel::detectCores() - 1
    }

    savelist = parallel::mclapply(mc.cores = cores, X = 1:nrow(combn_DT), FUN = function(row) {

      cell_type_1 = combn_DT[row][['V1']]
      cell_type_2 = combn_DT[row][['V2']]

      specific_scores = specificCellCellcommunicationScores(gobject = gobject,
                                                            cluster_column = cluster_column,
                                                            random_iter = random_iter,
                                                            cell_type_1 = cell_type_1,
                                                            cell_type_2 = cell_type_2,
                                                            gene_set_1 = gene_set_1,
                                                            gene_set_2 = gene_set_2,
                                                            spatial_network_name = spatial_network_name,
                                                            log2FC_addendum = log2FC_addendum,
                                                            min_observations = min_observations,
                                                            adjust_method = adjust_method,
                                                            adjust_target = adjust_target)

    })

  } else {

    ## for loop over all combinations ##
    savelist = list()
    countdown = nrow(combn_DT)

    for(row in 1:nrow(combn_DT)) {

      cell_type_1 = combn_DT[row][['V1']]
      cell_type_2 = combn_DT[row][['V2']]

      if(verbose == 'a little' | verbose == 'a lot') cat('\n\n PROCESS nr ', countdown,': ', cell_type_1, ' and ', cell_type_2, '\n\n')

      if(verbose %in% c('a little', 'none')) {
        specific_verbose = F
      } else {
        specific_verbose = T
      }

      specific_scores = specificCellCellcommunicationScores(gobject = gobject,
                                                            cluster_column = cluster_column,
                                                            random_iter = random_iter,
                                                            cell_type_1 = cell_type_1,
                                                            cell_type_2 = cell_type_2,
                                                            gene_set_1 = gene_set_1,
                                                            gene_set_2 = gene_set_2,
                                                            spatial_network_name = spatial_network_name,
                                                            log2FC_addendum = log2FC_addendum,
                                                            min_observations = min_observations,
                                                            adjust_method = adjust_method,
                                                            adjust_target = adjust_target,
                                                            verbose = specific_verbose)
      savelist[[row]] = specific_scores
      countdown = countdown - 1
    }

  }

  finalDT = do.call('rbind', savelist)
  data.table::setorder(finalDT, LR_comb, -LR_expr)

  return(finalDT)
}


#' @title combCCcom
#' @name combCCcom
#' @description Combine spatial and expression based cell-cell communication data.tables
#' @param spatialCC spatial cell-cell communication scores
#' @param exprCC expression cell-cell communication scores
#' @param min_lig_nr minimum number of ligand cells
#' @param min_rec_nr minimum number of receptor cells
#' @param min_padj_value minimum adjusted p-value
#' @param min_log2fc minimum log2 fold-change
#' @param min_av_diff minimum average expression difference
#' @return combined data.table with spatial and expression communication data
#' @export
#' @examples
#'     combCCcom(gobject)
combCCcom = function(spatialCC,
                     exprCC,
                     min_lig_nr = 3,
                     min_rec_nr = 3,
                     min_padj_value = 1,
                     min_log2fc = 0,
                     min_av_diff = 0) {


  spatialCC = spatialCC[lig_nr >= min_lig_nr & rec_nr >= min_rec_nr &
                          p.adj <= min_padj_value & abs(log2fc) >= min_log2fc & abs(av_diff) >= min_av_diff]

  data.table::setnames(x = spatialCC,
                       old = c('lig_expr', 'rec_expr', 'LR_expr', 'lig_nr', 'rec_nr',
                               'rand_expr', 'av_diff', 'log2fc', 'pvalue', 'p.adj', 'PI'),
                       new = c('lig_expr_spat', 'rec_expr_spat', 'LR_expr_spat', 'lig_nr_spat', 'rec_nr_spat',
                               'rand_expr_spat', 'av_diff_spat', 'log2fc_spat', 'pvalue_spat', 'p.adj_spat', 'PI_spat'))

  merge_DT = data.table:::merge.data.table(spatialCC, exprCC, by = c('LR_comb', 'LR_cell_comb',
                                                                    'lig_cell_type', 'rec_cell_type',
                                                                    'ligand', 'receptor'))

  # rank for expression levels
  merge_DT[, LR_expr_rnk := rank(-LR_expr), by = LR_comb]
  merge_DT[, LR_spat_rnk := rank(-LR_expr_spat), by = LR_comb]

  # rank for differential activity levels
  merge_DT[, exprPI_rnk := rank(-PI), by = LR_comb]
  merge_DT[, spatPI_rnk := rank(-PI_spat), by = LR_comb]

  return(merge_DT)

}

