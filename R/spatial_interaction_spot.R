## Giotto cell type enrichment for spots functions ####

# * ####
# cell type proximity for spots ####

#' @title cell_proximity_spots_internal
#' @name cell_proximity_spots_internal
#' @description Compute cell-cell interactions observed value inner each spot
#' @param cell_IDs cell_IDs get from gobject@cell_ID
#' @param dwls_values data.table of cell type enrichment in each spot and multiply
#' by cell number in each spot
#' @return List of cell proximity observed value in data.table format. Columns:
#' unified_int, internal value, type_int.
cell_proximity_spots_internal <- function(cell_IDs,
                                          dwls_values){
  
  proximity_dt = data.table::data.table()
  # calculate proximity for each spot
  for (cell_i in 1:length(cell_IDs)){
    cell_ID = cell_IDs[cell_i]
    # dwls value for one spot and remove 0 cell type
    dwls_spot = dwls_values[cell_ID,]
    dwls_spot = dwls_spot[dwls_spot > 0]
    
    # calculate proximity of same cell type (A==B)
    same_ct = data.table::data.table()
    if (length(dwls_spot) >= 1){
      same_ct = (dwls_spot-1) * dwls_spot / 2
      # transfer format
      unified_int_same = names(same_ct)
      unified_int_same = paste0(unified_int_same,'--',unified_int_same)
      same_ct = data.table::data.table('unified_int' = unified_int_same,'internal' = same_ct)
    }
    
    # calculate proximity of different cell type (A==B) 
    diff_ct = data.table::data.table()
    if (length(dwls_spot) >= 2){
      diff_ct = dwls_spot %o% dwls_spot
      #modifiy duplicate value
      diag(diff_ct) = NA
      diff_ct[lower.tri(diff_ct)] = NA
      # transfer format to data.table
      diff_ct = data.table::as.data.table(reshape2::melt(diff_ct))
      diff_ct = diff_ct[value != 'NA' ]
      diff_ct[, c('Var1', 'Var2') := lapply(.SD, as.character),.SDcols = c('Var1', 'Var2')]
      diff_ct[, unified_int := ifelse(Var1 < Var2, paste0(Var1,'--',Var2), paste0(Var2,'--',Var1))]
      diff_ct = diff_ct[, c('unified_int', 'value')]
      data.table::setnames(diff_ct, old = c('value'), new = c('internal'))
    }
    
    # merge spot proximity to proximity data.table
    proximity_dt = rbind(proximity_dt, same_ct, diff_ct)
  }
  
  proximity_dt = proximity_dt[internal > 0]
  proximity_dt[, internal := sum(internal), by=c('unified_int')]
  proximity_dt = unique(proximity_dt)
  
  return(proximity_dt)
}


#' @title cell_proximity_spots_external
#' @name cell_proximity_spots_external
#' @description Compute cell-cell interactions observed value for interacted spots
#' @param pairs data.table of paired spots. Format: cell_ID1, cell_ID2, N
#' @param dwls_values data.table of cell type enrichment in each spot and multiply
#' by cell number in each spot
#' @return List of cell proximity observed value in data.table format. Columns:
#' unified_int, internal value, type_int.
cell_proximity_spots_external <- function(pairs,
                                          dwls_values){
  
  cell_IDs = unique(c(pairs$from, pairs$to))
  pairs = pairs[, .N, by = c('from','to')]
  # add internal pairs to make full matirx 
  pairs_spots = data.table::data.table(from = cell_IDs, to = cell_IDs, N = 0)
  pairs_balance = data.table::data.table(from = pairs$to, to = pairs$from, N = pairs$N)
  pairs_for_mat = rbind(pairs_spots, pairs, pairs_balance)
  pairs_for_mat = pairs_for_mat[, .N, by = c('from','to')]
  
  # make square matrix of interaction between spots
  pairs_mat = reshape2::acast(pairs_for_mat, from ~ to, value.var = 'N' ,fill = 0)
  pairs_mat = pairs_mat[cell_IDs,cell_IDs]
  
  #calculate cell-tyep/cell-type interactions
  dwls_sub = dwls_values[cell_IDs,]
  proximity_dt = data.table::data.table()
  cts = colnames(dwls_sub)
  cts = sort(cts)
  for (i in 1:length(cts)){
    ct1 = cts[i]
    dwls_ct1 = dwls_sub[, ct1]
    
    for (j in i:length(cts)){
      ct2 = cts[j]
      dwls_ct2 = dwls_sub[, ct2]
      if (i == j ){f = 0.5}else{f=1}
      proximity_2cts = dwls_ct1 %o% dwls_ct2 * pairs_mat * f
      proximity_2cts = sum(proximity_2cts)
      proximity_2cts = data.table::data.table(unified_int = paste0(ct1,'--',ct2),
                                              external = proximity_2cts)
      proximity_dt = rbind(proximity_dt, proximity_2cts)
    }
  }
  return(proximity_dt)
}


#' @title cell_proximity_spots
#' @name cell_proximity_spots
#' @description Compute cell-cell interactions observed value for internal and external spots
#' @param cell_IDs cell_IDs to calculate internal cell-type/cell-type interactions
#' @param pairs data.table of paired spots. Format: cell_ID1, cell_ID2, N
#' @param dwls_values data.table of cell type enrichment in each spot and multiply
#' by cell number in each spot
#' @return List of cell proximity observed value in data.table format. Columns:
#' unified_int, type_int, V1, external, internal.
cell_proximity_spots <- function(cell_IDs,
                                 pairs_external,
                                 dwls_values){
  
  # compute cell-type/cell-type interactions in each spot (internal)
  if (length(cell_IDs) > 0){
    proximity_in = cell_proximity_spots_internal(cell_IDs = cell_IDs,
                                                 dwls_values = dwls_values)
  }
  
  # compute cell-type/cell-type interactions between spots (external)
  # get paired spots barcodes
  proximity_ex = cell_proximity_spots_external(pairs = pairs_external,
                                               dwls_values = dwls_values)
  
  if (length(cell_IDs) > 0) {
    proximity_dt = merge(proximity_ex, proximity_in, by= 'unified_int', all=TRUE)
  }else{
    proximity_dt = proximity_ex[, 'internal' := 0]
  }
  proximity_dt[is.na(proximity_dt)] = 0
  proximity_dt[, V1 := internal + external]
  
  proximity_dt[, s1 := strsplit(as.character(unified_int), split = '--')[[1]][1], by = 1:nrow(proximity_dt)]
  proximity_dt[, s2 := strsplit(as.character(unified_int), split = '--')[[1]][2], by = 1:nrow(proximity_dt)]
  proximity_dt[, type_int := ifelse(s1 == s2, 'homo', 'hetero')]
  proximity_dt = proximity_dt[, c('unified_int', 'type_int', 'V1', 'external', 'internal')]
  return(proximity_dt)
}


#' @title cellProximityEnrichmentSpots
#' @name cellProximityEnrichmentSpots
#' @description Compute cell-cell interaction enrichment for spots (observed vs expected)
#' @param gobject giotto object
#' @param spatial_network_name name of spatial network to use
#' @param cluster_column name of column to use for clusters
#' @param cells_in_spot cell number in each spot
#' @param number_of_simulations number of simulations to create expected observations
#' @param adjust_method method to adjust p.values
#' @param set_seed use of seed
#' @param seed_number seed number to use
#' @return List of cell Proximity scores (CPscores) in data.table format. The first
#' data.table (raw_sim_table) shows the raw observations of both the original and
#' simulated networks. The second data.table (enrichm_res) shows the enrichment results.
#' @details Spatial proximity enrichment or depletion between pairs of cell types
#' is calculated by calculating the observed over the expected frequency
#' of cell-cell proximity interactions. The expected frequency is the average frequency
#' calculated from a number of spatial network simulations. Each individual simulation is
#' obtained by reshuffling the cell type labels of each node (spot)
#' in the spatial network.
#' @export
cellProximityEnrichmentSpots <- function(gobject,
                                         spatial_network_name = 'spatial_network',
                                         cluster_column = 'cell_ID',
                                         cells_in_spot = 1,
                                         number_of_simulations = 100,
                                         adjust_method = c("none", "fdr", "bonferroni","BH",
                                                           "holm", "hochberg", "hommel",
                                                           "BY"),
                                         set_seed = TRUE,
                                         seed_number = 1234) {
  
  # p.adj test
  sel_adjust_method = match.arg(adjust_method, choices = c("none", "fdr", "bonferroni","BH",
                                                           "holm", "hochberg", "hommel",
                                                           "BY"))
  
  spatial_network_annot = annotateSpatialNetwork(gobject = gobject,
                                                 spatial_network_name = spatial_network_name,
                                                 cluster_column = cluster_column)
  
  # data.table variables
  unified_cells = type_int = N = NULL
  
  spatial_network_annot = Giotto:::sort_combine_two_DT_columns(spatial_network_annot, 'to', 'from', 'unified_cells')
  spatial_network_annot = spatial_network_annot[!duplicated(unified_cells)]
  
  # exact spatial_enrichment matrix
  dwls_values = gobject@spatial_enrichment$DWLS
  data.table::setDF(dwls_values)
  rownames_dwls = dwls_values[,'cell_ID']
  dwls_values = as.matrix(dwls_values[,-1])
  rownames(dwls_values) = rownames_dwls
  dwls_values_adjust = dwls_values * cells_in_spot
  
  # compute cell-type/cell-type interactions
  orig_pairs_external = spatial_network_annot[, .N, by = c('from', 'to')]
  table_orig_results = cell_proximity_spots(cell_IDs = gobject@cell_ID,
                                            pairs_external = orig_pairs_external, 
                                            dwls_values = dwls_values_adjust)
  table_orig_results[, orig := 'original']
  table_orig_results[, round := 'original']  
  
  # make simulated network
  sample_dt = make_simulated_network(gobject = gobject,
                                     spatial_network_name = spatial_network_name,
                                     cluster_column = cluster_column,
                                     number_of_simulations = number_of_simulations,
                                     set_seed = set_seed,
                                     seed_number = seed_number)
  
  # method for get simulation cell-type/cell-type interaction for each round
  data.table::setnames(sample_dt, old = c('s1', 's2'), new = c('from', 'to'))
  table_sim_results = NULL
  for(sim in 1:number_of_simulations) {
    r = paste0('sim',sim)
    sim_pairs = sample_dt[round == r, c("from","to")]    
    
    sim_cell_IDs = unique(sim_pairs[from == to, from])
    sim_pairs_ex = sim_pairs[from != to, ]
    sim_pairs_ex[, N :=1]
    
    sim_dt_round = cell_proximity_spots(cell_IDs = sim_cell_IDs,
                                        pairs_external = sim_pairs_ex, 
                                        dwls_values = dwls_values_adjust)
    
    sim_dt_round[, orig := 'simulations']
    sim_dt_round[, round := r]
    table_sim_results = rbind(table_sim_results, sim_dt_round)
  }
  
  
  table_results = rbind(table_orig_results, table_sim_results)
  
  # add missing combinations from original or simulations
  # probably not needed anymore
  all_sim_ints = as.character(unique(table_results[orig == 'simulations']$unified_int))
  all_orig_ints = as.character(unique(table_results[orig == 'original']$unified_int))
  missing_in_orig = all_sim_ints[!all_sim_ints %in% all_orig_ints]
  missing_in_sim = all_orig_ints[!all_orig_ints %in% all_sim_ints]
  create_missing_for_orig = table_results[unified_int %in% missing_in_orig]
  create_missing_for_orig = unique(create_missing_for_orig[, c('orig', 'V1') := list('original', 0)])
  create_missing_for_sim = table_results[unified_int %in% missing_in_sim]
  create_missing_for_sim = unique(create_missing_for_sim[, c('orig', 'V1') := list('simulations', 0)])
  
  table_results <- do.call('rbind', list(table_results, create_missing_for_orig, create_missing_for_sim))
  
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
      sim_values = c(sim_values, rep(0, additional_length_needed))
      #length_simulations = c(length_simulations, rep(0, additional_length_needed))
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
  
  # adjust p-values for mht
  
  # data.table variables
  p.adj_higher = p.adj_lower = p_lower_orig = p_higher_orig = PI_value = int_ranking = NULL
  
  table_mean_results_dc[, p.adj_higher := stats::p.adjust(p_higher_orig, method = sel_adjust_method)]
  table_mean_results_dc[, p.adj_lower := stats::p.adjust(p_lower_orig, method = sel_adjust_method)]
  
  
  table_mean_results_dc[, PI_value := ifelse(p.adj_higher <= p.adj_lower,
                                             -log10(p.adj_higher+(1/number_of_simulations))*enrichm,
                                             -log10(p.adj_lower+(1/number_of_simulations))*enrichm)]
  data.table::setorder(table_mean_results_dc, PI_value)
  
  # order
  table_mean_results_dc <- table_mean_results_dc[order(-PI_value)]
  table_mean_results_dc[, int_ranking := 1:.N]
  
  return(list(raw_sim_table = table_results, enrichm_res = table_mean_results_dc))
}

# * ####
# cell proximity with gene expression ####
#' @title geneExpDWLS
#' @name geneExpDWLS
#' @description Compute predicted gene expression value by spacialDWSL results and 
#' average gene expression for cell type
#' @param gobject gottio object
#' @param ave_celltype_exp data.table of gene expression in each cell type
#' @return matrix
geneExpDWLS = function(gobject,
                       ave_celltype_exp){
  
  # exact spatial_enrichment matrix
  dwls_values = gobject@spatial_enrichment$DWLS
  
  # 1. check if cell_type_vector and matrix are compatible
  if(ncol(ave_celltype_exp) != ncol(dwls_values) - 1) {
    stop('ncol(ave_celltype_exp) needs to be the same as ncol(dwls_values) - 1')
  }
  
  cell_types = colnames(ave_celltype_exp)
  data.table::setcolorder(dwls_values,c('cell_ID',cell_types))
  
  # 2. for each spot
  # calculate dwls predicted expression for genes 
  expMatrixDWLS = matrix(data = NA,
                         nrow = nrow(ave_celltype_exp),
                         ncol = nrow(dwls_values))
  
  average_exp = as.matrix(ave_celltype_exp)
  for(spot_i in 1:nrow(dwls_values)){
    spot = dwls_values[spot_i,1]
    spot_dwls = dwls_values[spot_i, -1]
    data.table::setDF(spot_dwls)
    spot_dwls = as.vector(t(spot_dwls)[,1])
    
    spot_exp = average_exp %*% spot_dwls
    
    expMatrixDWLS[, spot_i] = spot_exp
  }
  rownames(expMatrixDWLS) = rownames(ave_celltype_exp)
  colnames(expMatrixDWLS) = dwls_values[,'cell_ID']
  
  return(expMatrixDWLS)
}

#' @title cellProximityEnrichmentEachSpot
#' @name cellProximityEnrichmentEachSpot
#' @description Compute cell-cell interaction enrichment for each spot with its 
#' interacted spots (observed)
#' @param gobject giotto object
#' @param feat_type feature type
#' @param spatial_network_name name of spatial network to use
#' @param cluster_column name of column to use for clusters
#' @return matrix that rownames are cell-cell interaction pairs and colnames are cell_IDs
#' @export
cellProximityEnrichmentEachSpot <- function(gobject,
                                            spatial_network_name = 'spatial_network',
                                            cluster_column = 'cell_ID') {
  
  
  spatial_network_annot = annotateSpatialNetwork(gobject = gobject,
                                                 spatial_network_name = spatial_network_name,
                                                 cluster_column = cluster_column)
  
  # data.table variables
  unified_cells = type_int = N = NULL
  
  spatial_network_annot = Giotto:::sort_combine_two_DT_columns(spatial_network_annot, 'to', 'from', 'unified_cells')
  spatial_network_annot = spatial_network_annot[!duplicated(unified_cells)]
  
  # exact spatial_enrichment matrix
  dwls_values = gobject@spatial_enrichment$DWLS
  data.table::setDF(dwls_values)
  rownames_dwls = dwls_values[,'cell_ID']
  dwls_values = as.matrix(dwls_values[,-1])
  rownames(dwls_values) = rownames_dwls
  
  # calculate cell-cell interaction with gene expression in each spot
  # proximity_spots is the 
  # get cell-cell types pairs
  cts = colnames(dwls_values)
  ct_pairs = data.table::as.data.table(t(combn(cts,m=2)))
  self_pairs = data.table::data.table(V1 = cts, V2 = cts)
  ct_pairs = rbind(ct_pairs, self_pairs)
  ct_pairs[, unified_int := ifelse(V1 < V2, paste0(V1,'--',V2), paste0(V2,'--',V1)), by = 1:nrow(ct_pairs)]
  unified_int = ct_pairs$unified_int
  
  # get paired spots barcodes
  orig_pairs = spatial_network_annot[, .N, by = c('from', 'to')]  
  cell_IDs = unique(c(orig_pairs$from, orig_pairs$to))
  
  # make matrix that rows are cell-cell types and columns are cell_IDs 
  proximityMat = matrix(data = 0,
                        nrow = length(unified_int),
                        ncol = length(cell_IDs))
  
  rownames(proximityMat) = unified_int
  colnames(proximityMat) = cell_IDs
  
  # for each spot, calculate cell type proximity to it 
  for (cell_i in 1:length(cell_IDs)){
    cell_ID = cell_IDs[cell_i]
    spot_pairs = orig_pairs[from == cell_ID | to == cell_ID]
    
    spot_proximity = cell_proximity_spots(cell_IDs = cell_ID,
                                          pairs_external = spot_pairs, 
                                          dwls_values = dwls_values)
    
    spot_proximity = unique(spot_proximity[,c('unified_int','V1')])
    
    # add to proximityMat(matrix)
    proximityMat[spot_proximity$unified_int, cell_i] = spot_proximity$V1
  }
  return(proximityMat)
}


#' @title findCellProximityGenesForSpots
#' @name findCellProximityGenesForSpots
#' @description Compute the correlation between cell-cell interaction enrichment
#' and gene expression
#' @param gobject giotto object
#' @param feat_type feature type
#' @param expression_values expression values to use
#' @param ave_celltype_exp average gene expresion in each cell type
#' @param spatial_network_name name of spatial network to use
#' @param pcc_threshold pcc cut off 
#' @return data.table
#' @export

findCellProximityGenesForSpots <- function(gobject,
                                           expression_values = c('normalized', 'scaled', 'custom'),
                                           ave_celltype_exp,
                                           spatial_network_name = 'spatial_network',
                                           pcc_threshold = 0.1) {
  
  # expression data
  values = match.arg(expression_values, choices = unique(c('normalized', 'scaled', 'custom', expression_values)))
  expr_observed = select_expression_values(gobject = gobject, values = values)
  
  # Compute predicted gene expression value
  expr_predicted = geneExpDWLS(gobject = gobject,
                               ave_celltype_exp = ave_celltype_exp)
  
  # Get the difference expssion matrix between observed and predicted expression
  intersect_gene  = intersect(rownames(expr_predicted), rownames(expr_observed))
  expr_diff = expr_observed[intersect_gene,] - expr_predicted[intersect_gene,]
  expr_diff = as.matrix(expr_diff)
  
  # compute cell proximity for each spot
  proximityMat = cellProximityEnrichmentEachSpot(gobject = gobject, 
                                                 spatial_network_name = spatial_network_name)
  
  # select overlapped spots
  intersect_cell_IDs = intersect(colnames(expr_diff), colnames(proximityMat))
  expr_diff = expr_diff[, intersect_cell_IDs]
  proximityMat = proximityMat[, intersect_cell_IDs]
  
  # compute correlation between genes and cell-types to find ICGs
  pcc_mat = cor(t(expr_diff), t(proximityMat))
  pcc_dt = data.table::as.data.table(reshape2::melt(pcc_mat))
  data.table::setnames(pcc_dt, old = c('Var1', 'Var2', 'value'), new = c('feat', 'unified_int', 'pcc'))
  
  # select by pcc_threshold
  filter_table = pcc_dt[pcc >= pcc_threshold | pcc <= -pcc_threshold]
  
  # rank for pcc value
  filter_table[, 'rank_pos' := rank(pcc)]
  filter_table[, 'rank_neg' := rank(-pcc)]
  
  # add interacted cell types
  filter_table[, cell_type_1 := strsplit(x = as.character(unified_int), split = '--')[[1]][1], by = 1:nrow(filter_table)]
  filter_table[, cell_type_2 := strsplit(x = as.character(unified_int), split = '--')[[1]][2], by = 1:nrow(filter_table)]
  
  # add type int
  filter_table[, type_int := ifelse(cell_type_1 == cell_type_2, 'homo', 'hetero')]
  
  return(filter_table)
}

#' @title plotCCgeneDotplot
#' @name plotCCgeneDotplot
#' @description Plots dotplot for correlation scores between cell-cell interactions and gene expression
#' @param gobject giotto object
#' @param pccScores correlation scores from \code{\link{findCellProximityFeatsForSpots}}
#' @param select_feats selected genes
#' @param selected_cell_LR selected cell-cell combinations for genes
#' @param show_feat_names show gene names
#' @param show_cell_LR_names show cell-cell names
#' @param cluster_on values to use for clustering of cell-cell and gene pairs
#' @param cor_method correlation method used for clustering
#' @param aggl_method agglomeration method used by hclust
#' @param show_plot show plots
#' @param return_plot return plotting object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters from \code{\link{all_plots_save_function}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return ggplot
#' @export
plotCCgeneDotplot = function(gobject,
                             pccScores,
                             select_feats = NULL,
                             selected_cell_LR = NULL,
                             cluster_on = 'pcc',
                             cor_method = c("pearson", "kendall", "spearman"),
                             aggl_method = c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid"),
                             show_feat_names = TRUE,
                             show_cell_LR_names = TRUE,
                             show_plot = NA,
                             return_plot = NA,
                             save_plot = NA,
                             save_param =  list(),
                             default_save_name = 'plotCCgeneDotplot'){
  # get parameters
  cor_method = match.arg(cor_method, choices = c("pearson", "kendall", "spearman"))
  aggl_method = match.arg(aggl_method, choices = c("ward.D", "ward.D2", "single", "complete",
                                                   "average", "mcquitty", "median", "centroid"))
  
  
  # data.table variables
  feat = unified_int = sd = NULL
  
  # plot method
  if(!is.null(select_feats) & !is.null(selected_cell_LR)) {
    selDT = pccScores[feat %in% select_feats & unified_int %in% selected_cell_LR]
  } else if(!is.null(select_feats)) {
    selDT = pccScores[feat %in% select_feats]
  } else if(!is.null(selected_cell_LR)) {
    selDT = pccScores[unified_int %in% selected_cell_LR]
  } else {
    selDT = pccScores
  }
  selDT = selDT[, pcc_abs := abs(pcc), by = 1:nrow(selDT)]
  
  
  
  # creat matrix
  cluster_on = match.arg(cluster_on, choices = c('pcc'))
  selDT_d = data.table::dcast.data.table(selDT, unified_int~feat, value.var = cluster_on, fill = 0)
  selDT_m = Giotto:::dt_to_matrix(selDT_d)
  
  # remove zero variance
  sd_rows = apply(selDT_m, 1, sd)
  sd_rows_zero = names(sd_rows[sd_rows == 0])
  if(length(sd_rows_zero) > 0) selDT_m = selDT_m[!rownames(selDT_m) %in% sd_rows_zero, ]
  
  sd_cols = apply(selDT_m, 2, sd)
  sd_cols_zero = names(sd_cols[sd_cols == 0])
  if(length(sd_cols_zero) > 0) selDT_m = selDT_m[, !colnames(selDT_m) %in% sd_cols_zero]
  
  
  ## cells
  corclus_cells_dist = stats::as.dist(1-Giotto:::cor_flex(x = Giotto:::t_flex(selDT_m), method = cor_method))
  hclusters_cells = stats::hclust(d = corclus_cells_dist, method = aggl_method)
  clus_names = rownames(selDT_m)
  names(clus_names) = 1:length(clus_names)
  clus_sort_names = clus_names[hclusters_cells$order]
  selDT[, unified_int := factor(unified_int, clus_sort_names)]
  
  ## genes
  corclus_genes_dist = stats::as.dist(1-Giotto:::cor_flex(x = selDT_m, method = cor_method))
  hclusters_genes = stats::hclust(d = corclus_genes_dist, method = aggl_method)
  clus_names = colnames(selDT_m)
  names(clus_names) = 1:length(clus_names)
  clus_sort_names = clus_names[hclusters_genes$order]
  selDT[, feat := factor(feat, clus_sort_names)]
  
  
  
  pl = ggplot2::ggplot()
  pl = pl + ggplot2::geom_point(data = selDT, Giotto:::aes_string2(x = 'unified_int',
                                                                   y = 'feat', size = 'pcc_abs', color = 'pcc'))
  pl = pl + ggplot2::theme_classic()
  if(show_feat_names == TRUE) pl = pl + ggplot2::theme(axis.text.y = element_text(),
                                                       axis.ticks.y = element_line())
  if(show_cell_LR_names == TRUE) pl = pl + ggplot2::theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
                                                          axis.ticks.x = element_line())
  pl = pl + ggplot2::scale_fill_gradientn(colours = c('darkblue', 'blue', 'white', 'red', 'darkred'))
  pl = pl + ggplot2::scale_size_continuous(range = c(0.5, 5)) + scale_color_gradientn(colours = c('darkblue', 'blue', 'white', 'red', 'darkred'))
  pl = pl + ggplot2::labs(x = 'cell-cell', y = 'gene')
  
  
  # print, return and save parameters
  show_plot = ifelse(is.na(show_plot), readGiottoInstructions(gobject, param = 'show_plot'), show_plot)
  save_plot = ifelse(is.na(save_plot), readGiottoInstructions(gobject, param = 'save_plot'), save_plot)
  return_plot = ifelse(is.na(return_plot), readGiottoInstructions(gobject, param = 'return_plot'), return_plot)
  
  ## print plot
  if(show_plot == TRUE) {
    print(pl)
  }
  
  ## save plot
  if(save_plot == TRUE) {
    do.call('all_plots_save_function', c(list(gobject = gobject, plot_object = pl, default_save_name = default_save_name), save_param))
  }
  
  ## return plot
  if(return_plot == TRUE) {
    return(pl)
  }
  
}










