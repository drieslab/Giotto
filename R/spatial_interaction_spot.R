## Giotto cell type enrichment for spots functions ####

# * ####
# cell type proximity for spots ####

#' @title cell_proximity_spots_internal
#' @name cell_proximity_spots_internal
#' @description Compute cell-cell interactions observed value inner each spot
#' @param cell_IDs cell_IDs 
#' @param dwls_values data.table of cell type enrichment in each spot and multiply
#' by cell number in each spot
#' @return List of cell proximity observed value in data.table format. Columns:
#' unified_int, internal value, type_int.
#' @keywords internal
cell_proximity_spots_internal = function(cell_IDs,
                                         dwls_values){
  
  # data.table variables
  value = unified_int = Var1 = Var2 = internal = NULL

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
  proximity_dt[, internal := sum(internal), by = c('unified_int')]
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
#' @keywords internal
cell_proximity_spots_external = function(pairs,
                                         dwls_values){

  cell_IDs = unique(c(pairs$from, pairs$to))
  pairs = pairs[, .N, by = c('from','to')]
  # add internal pairs to make full matrix
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
#' @keywords internal
cell_proximity_spots = function(cell_IDs,
                                pairs_external,
                                dwls_values){
  
  # data.table variables
  V1 =  internal = external = s1 = s2 = unified_int = type_int = NULL
  
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
#'
#' @param gobject giotto object
#' @param spat_unit spatial unit (e.g. 'cell')
#' @param feat_type feature type (e.g. 'rna')
#' @param spatial_network_name name of spatial network to use
#' @param cluster_column name of column to use for clusters
#' @param cells_in_spot cell number in each spot
#' @param number_of_simulations number of simulations to create expected observations
#' @param adjust_method method to adjust p.values (e.g. "none", "fdr", "bonferroni","BH","holm", "hochberg", "hommel","BY")
#' @param set_seed use of seed. Default = TRUE
#' @param seed_number seed number to use. Default = 1234
#'
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
                                         spat_unit = NULL,
                                         feat_type = NULL,
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
                                                 spat_unit = spat_unit,
                                                 feat_type = feat_type,
                                                 spatial_network_name = spatial_network_name,
                                                 cluster_column = cluster_column)

  # data.table variables
  orig = from = to = unified_int = unified_cells = type_int = N = V1 = original = enrichm = simulations = NULL

  spatial_network_annot = sort_combine_two_DT_columns(spatial_network_annot, 'to', 'from', 'unified_cells')
  spatial_network_annot = spatial_network_annot[!duplicated(unified_cells)]

  # exact spatial_enrichment matrix
  dwls_values = get_spatial_enrichment(gobject = gobject,
                                       spat_unit = spat_unit,
                                       feat_type = feat_type,
                                       enrichm_name = 'DWLS',
                                       output = 'data.table')
  data.table::setDF(dwls_values)
  rownames_dwls = dwls_values[,'cell_ID']
  dwls_values = as.matrix(dwls_values[,-1])
  rownames(dwls_values) = rownames_dwls
  dwls_values_adjust = dwls_values * cells_in_spot

  # compute cell-type/cell-type interactions
  print("1/5 Computing cell-type/cell-type interactions")

  orig_pairs_external = spatial_network_annot[, .N, by = c('from', 'to')]
  table_orig_results = cell_proximity_spots(cell_IDs = pDataDT(gobject)$cell_ID,
                                            pairs_external = orig_pairs_external,
                                            dwls_values = dwls_values_adjust)
  table_orig_results[, orig := 'original']
  table_orig_results[, round := 'original']

  # make simulated network
  print("2/5 Make simulated network")

  sample_dt = make_simulated_network(gobject = gobject,
                                     spat_unit = spat_unit,
                                     feat_type = feat_type,
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
  print("3/5 Calculating p-values")

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
  print("4/5 Depletion or enrichment in barplot format")

  table_mean_results <- table_results[, .(mean(V1)), by = c('orig', 'unified_int', 'type_int')]
  table_mean_results_dc <- data.table::dcast.data.table(data = table_mean_results, formula = type_int+unified_int~orig, value.var = 'V1')
  table_mean_results_dc[, original := ifelse(is.na(original), 0, original)]
  table_mean_results_dc[, enrichm := log2((original+1)/(simulations+1))]


  table_mean_results_dc <- merge(table_mean_results_dc, res_pvalue_DT, by = 'unified_int')
  data.table::setorder(table_mean_results_dc, enrichm)
  table_mean_results_dc[, unified_int := factor(unified_int, unified_int)]

  # adjust p-values for mht

  print("5/5 Calculating adjust p-values for mht")

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
# cell proximity with feature expression ####

#' @title featExpDWLS
#' @name featExpDWLS
#' @description Compute predicted feature expression value by spatialDWSL results and
#' average feature expression for cell type
#'
#' @param gobject giotto object
#' @param spat_unit spatial unit (e.g. 'cell')
#' @param feat_type feature type (e.g. 'rna')
#' @param ave_celltype_exp data.table of feature expression in each cell type
#'
#' @return matrix
#' @export
featExpDWLS = function(gobject,
                       spat_unit = NULL,
                       feat_type = NULL,
                       ave_celltype_exp){

  # exact spatial_enrichment matrix
  dwls_values = get_spatial_enrichment(gobject,
                                       spat_unit = spat_unit,
                                       feat_type = feat_type,
                                       enrichm_name = 'DWLS',
                                       output = 'data.table')

  # 1. check if cell_type_vector and matrix are compatible
  if(ncol(ave_celltype_exp) != ncol(dwls_values) - 1) {
    stop('ncol(ave_celltype_exp) needs to be the same as ncol(dwls_values) - 1')
  }

  cell_types = colnames(ave_celltype_exp)
  data.table::setcolorder(dwls_values,c('cell_ID',cell_types))

  # 2. for each spot
  # calculate dwls predicted expression for features
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
  colnames(expMatrixDWLS) = dwls_values$cell_ID

  return(expMatrixDWLS)
}


#' @title cal_expr_residual
#' @name cal_expr_residual
#' @description Calculate feature expression residual (observed_exp - DWLS_predicted)
#'
#' @param gobject giotto object
#' @param spat_unit spatial unit (e.g. 'cell')
#' @param feat_type feature type (e.g. 'rna')
#' @param expression_values expression values to use (e.g. 'normalized', 'scaled', 'custom')
#' @param ave_celltype_exp average expression matrix in cell types
#'
#' @keywords internal
cal_expr_residual <- function(gobject,
                              spat_unit = NULL,
                              feat_type = NULL,
                              expression_values = c('normalized', 'scaled', 'custom'),
                              ave_celltype_exp) {

  # expression data
  values = match.arg(expression_values, choices = c('normalized', 'scaled', 'custom'))
  # expr_observed = get_expression_values(gobject = gobject,
  #                                       values = expression_values,
  #                                       output = 'matrix')
  expr_observed = slot(gobject@expression[[spat_unit]][[feat_type]][[values]], 'exprMat')

  # Compute predicted feature expression value
  expr_predicted = featExpDWLS(gobject = gobject,
                               spat_unit = spat_unit,
                               feat_type = feat_type,
                               ave_celltype_exp = ave_celltype_exp)

  # Get the difference expression matrix between observed and predicted expression
  intersect_feature  = intersect(rownames(expr_predicted), rownames(expr_observed))
  expr_residual = expr_observed[intersect_feature,] - expr_predicted[intersect_feature,]
  expr_residual = as.matrix(expr_residual)

  return(expr_residual)
}


#' @title cellProximityEnrichmentEachSpot
#' @name cellProximityEnrichmentEachSpot
#' @description Compute cell-cell interaction enrichment for each spot with its
#' interacted spots (observed)
#'
#' @param gobject giotto object
#' @param spat_unit spatial unit (e.g. 'cell')
#' @param feat_type feature type (e.g. 'rna')
#' @param spatial_network_name name of spatial network to use
#' @param cluster_column name of column to use for clusters
#'
#' @return matrix that rownames are cell-cell interaction pairs and colnames are cell_IDs
#' @export
cellProximityEnrichmentEachSpot <- function(gobject,
                                            spat_unit = NULL,
                                            feat_type = NULL,
                                            spatial_network_name = 'spatial_network',
                                            cluster_column = 'cell_ID') {
  
  spatial_network_annot = annotateSpatialNetwork(gobject = gobject,
                                                 spat_unit = spat_unit,
                                                 feat_type = feat_type,
                                                 spatial_network_name = spatial_network_name,
                                                 cluster_column = cluster_column)

  # data.table variables
  V1 = V2 = from = to = int_cell_IDS = Var1 = Var2 = unified_cells = type_int = N = NULL

  spatial_network_annot = sort_combine_two_DT_columns(spatial_network_annot, 'to', 'from', 'unified_cells')
  spatial_network_annot = spatial_network_annot[!duplicated(unified_cells)]

  # exact spatial_enrichment matrix
  dwls_values = get_spatial_enrichment(gobject = gobject,
                                       spat_unit = spat_unit,
                                       feat_type = feat_type,
                                       enrichm_name = 'DWLS',
                                       output = 'data.table')
  data.table::setDF(dwls_values)
  rownames_dwls = dwls_values[,'cell_ID']
  dwls_values = as.matrix(dwls_values[,-1])
  rownames(dwls_values) = rownames_dwls

  # calculate cell-cell interaction with feature expression in each spot
  # proximity_spots is the
  # get cell-cell types pairs
  cts = colnames(dwls_values)
  ct_pairs = data.table::data.table(V1 = rep(cts,each = length(cts)), V2 = rep(cts,length(cts)))
  ct_pairs[, unified_int := paste0(V1,'--',V2), by = 1:nrow(ct_pairs)]
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
    spot_pairs[, int_cell_IDS := ifelse(from==cell_ID, to, from)]
    int_num = spot_pairs$N

    dwls_target_cell = dwls_values[cell_ID,]
    dwls_int_cells = dwls_values[spot_pairs$int_cell_IDS,]

    # filter 0 and kept the data type
    # rowSum(dwls) = c(1,1,1,1.....)
    idx1 = which(dwls_target_cell > 0) #length(idx) must > 0
    dwls_target_cell = dwls_target_cell[idx1]

    if (length(int_num) > 1){
      idx2 = which(colSums(dwls_int_cells) > 0)
      dwls_int_cells = dwls_int_cells[, idx2]
      
      # all the interacted cells dwls have same cell type with proportion=1
      if (length(idx2) == 1){
        dwls_int_cells = matrix(dwls_int_cells, ncol = 1,
                                dimnames = list(spot_pairs$int_cell_IDS,names(idx2)))
      }
      
    } else{
      # target cell only contain 1 inteacted cell
      idx2 = which(dwls_int_cells > 0)
      dwls_int_cells = dwls_int_cells[idx2]
      dwls_int_cells = matrix(dwls_int_cells,nrow=1,byrow = TRUE,
                              dimnames = list(spot_pairs$int_cell_IDS,names(dwls_int_cells)))
    }


    spot_proximity = dwls_target_cell %o% (dwls_int_cells * int_num)
    spot_proximity = apply(spot_proximity, 3, rowSums)
    if (length(dwls_target_cell) == 1){
      # change to the right data class
      spot_proximity = matrix(spot_proximity,nrow=1,byrow = TRUE,
                              dimnames = list(names(dwls_target_cell),names(spot_proximity)))

    }
    spot_proximity = reshape2::melt(spot_proximity)
    spot_proximity = data.table::data.table(spot_proximity)
    spot_proximity[, c('Var1', 'Var2') := lapply(.SD, as.character),.SDcols = c('Var1', 'Var2')]
    spot_proximity[, unified_int := paste0(Var1,'--',Var2)]

    # add to proximityMat(matrix)
    proximityMat[spot_proximity$unified_int, cell_i] = spot_proximity$value
  }
  return(proximityMat)
}

#' @title cal_diff_per_interaction
#' @name cal_diff_per_interaction
#' @description calculate correlation between expression residual and
#' cell proximity score of selected cell for spots
#' @keywords internal
cal_diff_per_interaction <- function(sel_int,
                                     other_ints,
                                     select_ind,
                                     other_ind,
                                     proximityMat,
                                     expr_residual){
  
  pcc_diff <- sel <- other <- NULL

  # get data

  prox_sel = proximityMat[sel_int, select_ind]
  prox_sel = as.matrix(prox_sel)
  expr_sel = expr_residual[, select_ind]

  prox_other = proximityMat[other_ints,other_ind]
  prox_other = prox_other[rowSums(prox_other) != 0, ]
  expr_other = expr_residual[, other_ind]

  # calculate pcc between expresidiual and proximity
  pcc_sel = stats::cor(t(expr_sel), prox_sel)

  pcc_other = stats::cor(t(expr_other), t(prox_other))
  pcc_other = rowMeans(pcc_other)

  features = rownames(pcc_sel)
  pcc_dt = data.table::data.table(features = features,
                                  pcc_sel = as.vector(pcc_sel),
                                  pcc_other = pcc_other[features])

  pcc_dt[, pcc_diff := pcc_sel - pcc_other]

  # calculate mean exression residual
  expr_sel_mean = rowMeans(expr_sel)
  expr_other_mean = rowMeans(expr_other)
  expr_residual_dt = data.table::data.table(features = features,
                                            sel = expr_sel_mean[features],
                                            other = expr_other_mean[features])
  expr_residual_dt[, diff := sel - other]

  results_dt = data.table::merge.data.table(expr_residual_dt, pcc_dt, by = 'features')

  return(results_dt)
}

#' @title do_permuttest_original_spot
#' @name do_permuttest_original_spot
#' @description calculate original values for spots
#' @keywords internal
do_permuttest_original_spot <- function(sel_int,
                                        other_ints,
                                        select_ind,
                                        other_ind,
                                        name = 'orig',
                                        proximityMat,
                                        expr_residual) {

  resultsDT = cal_diff_per_interaction(sel_int = sel_int,
                                       other_ints = other_ints,
                                       select_ind = select_ind,
                                       other_ind = other_ind,
                                       proximityMat = proximityMat,
                                       expr_residual = expr_residual)
  resultsDT[, name := name]
  return(resultsDT)
}

#' @title do_permuttest_random_spot
#' @name do_permuttest_random_spot
#' @description calculate random values for spots
#' @keywords internal
do_permuttest_random_spot <- function(sel_int,
                                      other_ints,
                                      select_ind,
                                      other_ind,
                                      name = 'perm_1',
                                      proximityMat,
                                      expr_residual,
                                      set_seed = TRUE,
                                      seed_number = 1234) {

  # data.table variables
  features = NULL

  l_sel_int = length(sel_int)
  l_other_ints = length(other_ints)
  l_select_ind = length(select_ind)
  l_other_ind = length(other_ind)

  all_IDs = colnames(proximityMat)
  all_ints = rownames(proximityMat)
  all_ints = all_ints[!rownames(proximityMat) %in% sel_int]

  if(set_seed == TRUE) {
    set.seed(seed = seed_number)
  }
  random_sel_int = sample(all_ints, size = l_sel_int, replace = F)
  random_other_ints = sample(all_ints, size = l_other_ints, replace = F)

  # keep the random selete not all the zeros
  prox = proximityMat[random_sel_int,]
  prox = prox[prox>0]
  random_select = c(sample(all_IDs, size = l_select_ind - 1, replace = F), names(prox[1]))
  random_other = c(sample(all_IDs, size = l_other_ind, replace = F), names(prox[length(prox)]))

  resultsDT = cal_diff_per_interaction(sel_int = random_sel_int,
                                       other_ints = random_other_ints,
                                       select_ind = random_select,
                                       other_ind = random_other,
                                       proximityMat = proximityMat,
                                       expr_residual = expr_residual)
  resultsDT[, name := name]

  return(resultsDT)
}


#' @title do_multi_permuttest_random_spot
#' @name do_multi_permuttest_random_spot
#' @description calculate multiple random values for spots
#' @keywords internal
do_multi_permuttest_random_spot = function(sel_int,
                                           other_ints,
                                           select_ind,
                                           other_ind,
                                           proximityMat,
                                           expr_residual,
                                           n = 100,
                                           cores = 2,
                                           set_seed = TRUE,
                                           seed_number = 1234) {

  if(set_seed == TRUE) {
    seed_number_list = seed_number:(seed_number + (n-1))
  }

  result = giotto_lapply(X = 1:n, cores = cores, fun = function(x) {

    seed_number = seed_number_list[x]

    perm_rand = do_permuttest_random_spot(sel_int = sel_int,
                                          other_ints = other_ints,
                                          select_ind = select_ind,
                                          other_ind = other_ind,
                                          name = paste0('perm_', x),
                                          proximityMat =  proximityMat,
                                          expr_residual = expr_residual,
                                          set_seed = set_seed,
                                          seed_number = seed_number)

  })

  final_result = do.call('rbind', result)

}

#' @title do_permuttest_spot
#' @name do_permuttest_spot
#' @description Performs permutation test on subsets of a matrix for spots
#' @keywords internal
do_permuttest_spot = function(sel_int,
                              other_ints,
                              select_ind,
                              other_ind,
                              proximityMat,
                              expr_residual,
                              n_perm = 100,
                              adjust_method = 'fdr',
                              cores = 2,
                              set_seed = TRUE,
                              seed_number = 1234) {
  
  

  # data.table variables
  log2fc_diff = log2fc = sel = other = features = p_higher = p_lower = perm_sel = NULL
  perm_other = perm_log2fc = perm_diff = p.value = p.adj = pcc_sel = pcc_diff = NULL
  perm_pcc_sel = perm_pcc_diff = pcc_other = NULL

  ## original data
  original = do_permuttest_original_spot(sel_int = sel_int,
                                         other_ints = other_ints ,
                                         select_ind = select_ind,
                                         other_ind = other_ind,
                                         name = 'orig',
                                         proximityMat = proximityMat,
                                         expr_residual = expr_residual)

  ## random permutations
  random_perms = do_multi_permuttest_random_spot(sel_int = sel_int,
                                                 other_ints = other_ints,
                                                 select_ind = select_ind,
                                                 other_ind = other_ind,
                                                 proximityMat = proximityMat,
                                                 expr_residual = expr_residual,
                                                 n = n_perm,
                                                 cores = cores,
                                                 set_seed = set_seed,
                                                 seed_number = seed_number)

  ##
  #random_perms[, log2fc_diff := rep(original$log2fc, n_perm) - log2fc]
  random_perms[, c('perm_sel', 'perm_other', 'perm_pcc_sel', 'perm_pcc_diff') := list(mean(sel), mean(other), mean(pcc_sel), mean(pcc_diff)), by = features]

  ## get p-values
  random_perms[, p_higher := sum(pcc_diff > 0), by = features]
  random_perms[, p_higher := 1-(p_higher/n_perm)]
  random_perms[, p_lower := sum(pcc_diff < 0), by = features]
  random_perms[, p_lower := 1-(p_lower/n_perm)]

  ## combine results permutation and original
  random_perms_res = unique(random_perms[,.(features, perm_sel, perm_other, perm_pcc_sel, perm_pcc_diff, p_higher, p_lower)])
  results_m = data.table::merge.data.table(random_perms_res, original[,.(features, sel, other, diff, pcc_sel, pcc_other, pcc_diff)], by = 'features')

  # select lowest p-value and perform p.adj
  results_m[, p.value := ifelse(p_higher <= p_lower, p_higher, p_lower)]
  results_m[, p.adj := stats::p.adjust(p.value, method = adjust_method)]

  results_m = results_m[,.(features, sel, other, pcc_sel, pcc_other, pcc_diff, p.value, p.adj, perm_sel, perm_other, perm_pcc_sel, perm_pcc_diff)]
  setorder(results_m, p.adj, -pcc_diff)

  return(results_m)

}


#' @title do_cell_proximity_test_spot
#' @name do_cell_proximity_test_spot
#' @description Performs a selected differential test on subsets of a matrix for spots
#' @keywords internal
do_cell_proximity_test_spot = function(sel_int,
                                       other_ints,
                                       select_ind,
                                       other_ind,
                                       proximityMat,
                                       expr_residual,
                                       diff_test,
                                       n_perm = 100,
                                       adjust_method = 'fdr',
                                       cores = 2,
                                       set_seed = TRUE,
                                       seed_number = 1234) {

  # get parameters
  diff_test = match.arg(diff_test, choices = c('permutation', 'limma', 't.test', 'wilcox'))
  adjust_method = match.arg(adjust_method, choices = c("bonferroni","BH", "holm", "hochberg", "hommel",
                                                       "BY", "fdr", "none"))


  if(diff_test == 'permutation') {
    result = do_permuttest_spot(sel_int = sel_int,
                                other_ints = other_ints,
                                select_ind = select_ind,
                                other_ind = other_ind,
                                proximityMat = proximityMat,
                                expr_residual = expr_residual,
                                n_perm = n_perm,
                                adjust_method = adjust_method,
                                cores = cores,
                                set_seed = set_seed,
                                seed_number = seed_number)

  }
  return(result)

}

#' @title findICF_per_interaction_spot
#' @name findICF_per_interaction_spot
#' @description Identifies features that are differentially expressed due to proximity to other cell types for spots.
#' @keywords internal
findICF_per_interaction_spot <- function(sel_int,
                                         all_ints,
                                         proximityMat,
                                         expr_residual,
                                         dwls_values,
                                         dwls_cutoff = 0.001,
                                         CCI_cell_score = 0.01,
                                         minimum_unique_cells = 1,
                                         minimum_unique_int_cells = 1,
                                         diff_test = 'permutation',
                                         n_perm = 100,
                                         adjust_method = 'fdr',
                                         cores = 2,
                                         set_seed = TRUE,
                                         seed_number = 1234){
  
  # data.table variables
  cell_type = int_cell_type = nr_select = int_nr_select = unif_int = unified_int = NULL

  sel_ct = strsplit(sel_int, '--')[[1]][1]
  int_ct = strsplit(sel_int, '--')[[1]][2]

  # filter out cells that without these two cellsltype
  prox_sel = proximityMat[sel_int,]
  prox_sel = prox_sel[which(prox_sel != 0)]
  prox_sel = prox_sel[which(prox_sel > CCI_cell_score)]
  spec_IDs = names(prox_sel)

  # find other cells contribution to cell type
  dwls_all_cell = dwls_values[, sel_ct]
  dwls_all_cell = dwls_all_cell[dwls_all_cell > dwls_cutoff]
  all_IDs = intersect(names(dwls_all_cell), colnames(proximityMat))
  other_IDs = setdiff(all_IDs, spec_IDs)

  other_ints = all_ints[cell_type == sel_ct]$unified_int

  other_ints = other_ints[-which(other_ints == sel_int)]

  ## do not continue if too few cells ##
  if(length(spec_IDs) < minimum_unique_cells | length(other_IDs) < minimum_unique_cells) {
    result = NULL
  } else {
    result = do_cell_proximity_test_spot(sel_int = sel_int,
                                         other_ints = other_ints,
                                         select_ind = spec_IDs,
                                         other_ind = other_IDs,
                                         proximityMat = proximityMat,
                                         expr_residual = expr_residual,
                                         diff_test = diff_test,
                                         n_perm = n_perm,
                                         adjust_method = adjust_method,
                                         cores = cores,
                                         set_seed = set_seed,
                                         seed_number = seed_number)

    result[, cell_type := sel_ct]
    result[, int_cell_type := int_ct]
    result[, nr_select := length(spec_IDs)]
    result[, int_nr_select := length(other_IDs)]
    result[, unif_int := sel_int]
  }

  return(result)

}

#' @title giotto_lapply
#' @keywords internal
giotto_lapply = function(X, cores = NA, fun, ...) {

  # get type of os
  os = .Platform$OS.type

  # set number of cores automatically, but with limit of 10
  cores = determine_cores(cores)

  if(os == 'unix') {
    save_list = parallel::mclapply(X = X, mc.cores = cores,
                                   FUN = fun, ...)
  } else if(os == 'windows') {
    save_list = parallel::mclapply(X = X, mc.cores = 1,
                                   FUN = fun, ...)

    # !! unexplainable errors are returned for some nodes !! #
    # currently disabled #
    #cl <- parallel::makeCluster(cores)
    #save_list = parallel::parLapply(cl = cl, X = X,
    #                                fun = fun, ...)
  }

  return(save_list)
}


#' @title findICFSpot
#' @name findICFSpot
#' @description Identifies cell-to-cell Interaction Changed Features (ICF) for spots,
#' i.e. features expression residual that are different due to proximity to other cell types.
#'
#' @param gobject A giotto object
#' @param spat_unit spatial unit (e.g. 'cell')
#' @param feat_type feature type (e.g. 'rna')
#' @param expression_values expression values to use
#' @param ave_celltype_exp average feature expression in each cell type
#' @param selected_features subset of selected features (optional)
#' @param spatial_network_name name of spatial network to use
#' @param minimum_unique_cells minimum number of target cells required
#' @param minimum_unique_int_cells minimum number of interacting cells required
#' @param CCI_cell_score cell proximity score to filter no interacted cell
#' @param dwls_cutoff cell type proportion cutoff to label the cell
#' @param diff_test which differential expression test
#' @param adjust_method which method to adjust p-values
#' @param nr_permutations number of permutations if diff_test = permutation
#' @param do_parallel run calculations in parallel with mclapply
#' @param cores number of cores to use if do_parallel = TRUE
#' @param set_seed set a seed for reproducibility
#' @param seed_number seed number
#'
#' @return icfObject that contains the differential feat scores
#' @details Function to calculate if features expression residual are differentially expressed in cell types
#'  when they interact (approximated by physical proximity) with other cell types.
#'  Feature expression residual calculated as:
#'  (observed expression in spot - cell_type_proportion * average_expressed_in_cell_type)
#'  The results data.table in the icfObject contains - at least - the following columns:
#' \itemize{
#'  \item{features:}{ All or selected list of tested features}
#'  \item{sel:}{ average feature expression residual in the interacting cells from the target cell type }
#'  \item{other:}{ average feature expression residual in the NOT-interacting cells from the target cell type }
#'  \item{pcc_sel:}{ correlation between cell proximity score and expression residual in the interacting cells from the target cell type}
#'  \item{pcc_other:}{ correlation between cell proximity score and expression residual in the NOT-interacting cells from the target cell type }
#'  \item{pcc_diff:}{ correlation difference between sel and other}
#'  \item{p.value:}{ associated p-value}
#'  \item{p.adj:}{ adjusted p-value}
#'  \item{cell_type:}{ target cell type}
#'  \item{int_cell_type:}{ interacting cell type}
#'  \item{nr_select:}{ number of cells for selected target cell type}
#'  \item{int_nr_select:}{ number of cells for interacting cell type}
#'  \item{unif_int:}{ cell-cell interaction}
#' }
#' @export
findICFSpot <- function(gobject,
                        spat_unit = NULL,
                        feat_type = NULL,
                        expression_values = c('normalized', 'scaled', 'custom'),
                        ave_celltype_exp,
                        selected_features = NULL,
                        spatial_network_name = 'Delaunay_network',
                        minimum_unique_cells = 5,
                        minimum_unique_int_cells = 5,
                        CCI_cell_score = 0.1,
                        dwls_cutoff = 0.001,
                        diff_test = 'permutation',
                        nr_permutations = 100,
                        adjust_method = 'fdr',
                        do_parallel = TRUE,
                        cores = 2,
                        set_seed = TRUE,
                        seed_number = 1234) {
  
  # data.table variables
  unified_int = NULL

  # expression data
  values = match.arg(expression_values, choices = c('normalized', 'scaled', 'custom'))
  features_overlap = intersect(slot(gobject, "feat_ID")[[feat_type]], rownames(ave_celltype_exp))
  ave_celltype_exp_sel = ave_celltype_exp[features_overlap,]
  expr_residual = cal_expr_residual(gobject = gobject,
                                    spat_unit = spat_unit,
                                    feat_type = feat_type,
                                    ave_celltype_exp = ave_celltype_exp_sel)

  ## test selected features ##
  if(!is.null(selected_features)) {
    expr_residual = expr_residual[rownames(expr_residual) %in% selected_features, ]
  }

  # compute cell proximity for each spot
  proximityMat = cellProximityEnrichmentEachSpot(gobject = gobject,
                                                 spatial_network_name = spatial_network_name)
  # select overlapped spots
  #intersect_cell_IDs = intersect(colnames(expr_residual), colnames(proximityMat))
  #expr_residual = expr_residual[, intersect_cell_IDs]
  #proximityMat = proximityMat[, intersect_cell_IDs]

  # compute correlation between features and cell-types to find ICFs
  all_ints = data.table::data.table(unified_int = rownames(proximityMat))
  all_ints[, cell_type := strsplit(as.character(unified_int), '--')[[1]][1], by = 1:nrow(all_ints)]
  all_ints[, int_cell_type := strsplit(as.character(unified_int), '--')[[1]][2], by = 1:nrow(all_ints)]

  # exact spatial_enrichment matrix
  dwls_values = get_spatial_enrichment(gobject = gobject,
                                       spat_unit = spat_unit,
                                       feat_type = feat_type,
                                       output = 'data.table')
  data.table::setDF(dwls_values)
  rownames_dwls = dwls_values[,'cell_ID']
  dwls_values = as.matrix(dwls_values[,-1])
  rownames(dwls_values) = rownames_dwls


  if(do_parallel == TRUE) {


    fin_result = giotto_lapply(X = all_ints$unified_int, cores = cores, fun = function(x) {

      tempres = findICF_per_interaction_spot(sel_int = x,
                                             all_ints = all_ints,
                                             proximityMat = proximityMat,
                                             expr_residual = expr_residual,
                                             dwls_values = dwls_values,
                                             dwls_cutoff = dwls_cutoff,
                                             CCI_cell_score = CCI_cell_score,
                                             minimum_unique_cells = minimum_unique_cells,
                                             minimum_unique_int_cells = minimum_unique_int_cells,
                                             n_perm = nr_permutations,
                                             adjust_method = adjust_method,
                                             cores = cores,
                                             set_seed = set_seed,
                                             seed_number = seed_number)
    })


  } else {

    fin_result = list()

    for(i in 1:length(all_ints$unified_int)) {

      x = all_ints$unified_int[i]
      print(x)

      tempres = findICF_per_interaction_spot(sel_int = x,
                                             all_ints = all_ints,
                                             proximityMat = proximityMat,
                                             expr_residual = expr_residual,
                                             dwls_values = dwls_values,
                                             dwls_cutoff = dwls_cutoff,
                                             CCI_cell_score = CCI_cell_score,
                                             minimum_unique_cells = minimum_unique_cells,
                                             minimum_unique_int_cells = minimum_unique_int_cells,
                                             n_perm = nr_permutations,
                                             adjust_method = adjust_method,
                                             cores = 2,
                                             set_seed = set_seed,
                                             seed_number = seed_number)

      fin_result[[i]] = tempres

    }

  }
  
  final_result = do.call('rbind', fin_result)

  # data.table variables
  spec_int = cell_type = int_cell_type = type_int = NULL

  final_result[, spec_int := paste0(cell_type,'--',int_cell_type)]
  final_result[, type_int := ifelse(cell_type == int_cell_type, 'homo', 'hetero')]


  #return(final_result)

  permutation_test = ifelse(diff_test == 'permutation', nr_permutations, 'no permutations')

  icfObject = list(ICFscores = final_result,
                   Giotto_info = list('values' = values,
                                      'cluster' = 'cell_ID',
                                      'spatial network' = spatial_network_name),
                   test_info = list('test' = diff_test,
                                    'p.adj' = adjust_method,
                                    'min cells' = minimum_unique_cells,
                                    'min interacting cells' = minimum_unique_int_cells,
                                    'perm' = permutation_test))
  class(icfObject) = append(class(icfObject), 'icfObject')
  return(icfObject)

}


#' @title filterICFSpot
#' @name filterICFSpot
#' @description Filter Interaction Changed Feature scores for spots.
#'
#' @param icfObject ICF (interaction changed feature) score object
#' @param min_cells minimum number of source cell type
#' @param min_cells_expr_resi minimum expression residual level for source cell type
#' @param min_int_cells minimum number of interacting neighbor cell type
#' @param min_int_cells_expr_resi minimum expression residual level for interacting neighbor cell type
#' @param min_fdr minimum adjusted p-value
#' @param min_pcc_diff minimum absolute pcc difference difference
#' @param min_zscore minimum z-score change
#' @param zscores_column calculate z-scores over cell types or features
#' @param direction differential expression directions to keep
#'
#' @return icfObject that contains the filtered differential feature scores
#' @export
filterICFSpot = function(icfObject,
                         min_cells = 4,
                         min_cells_expr_resi = 0.05,
                         min_int_cells = 4,
                         min_int_cells_expr_resi = 0.05,
                         min_fdr = 0.5,
                         min_pcc_diff = 0.05,
                         min_zscore = 0.05,
                         zscores_column = c('cell_type', 'features'),
                         direction = c('both', 'up', 'down')) {

  # data.table variables
  nr_select = int_nr_select = zscores = pcc_diff = sel = other = p.adj = NULL
  log2fc = min_log2_fc = NULL

  if(!'icfObject' %in% class(icfObject)) {
    stop('\n icfObject needs to be the output from findInteractionChangedFeats() or findICF() \n')
  }

  zscores_column = match.arg(zscores_column, choices = c('cell_type', 'features'))

  ICFscore = copy(icfObject[['ICFscores']])

  # other parameters
  direction = match.arg(direction, choices = c('both', 'up', 'down'))


  ## sequential filter steps ##
  # 1. minimum number of source and target cells
  selection_scores = ICFscore[nr_select >= min_cells & int_nr_select >= min_int_cells]

  # 2. create z-scores for log2fc per cell type
  selection_scores[, zscores := scale(pcc_diff), by = c(zscores_column)]

  # 3. filter based on z-scores and minimum levels
  comb_DT = rbind(selection_scores[zscores >= min_zscore & abs(pcc_diff) >= min_pcc_diff & sel >= min_cells_expr_resi],
                  selection_scores[zscores <= -min_zscore & abs(pcc_diff) >= min_pcc_diff & other >= min_int_cells_expr_resi])

  # 4. filter based on adjusted p-value (fdr)
  comb_DT = comb_DT[p.adj < min_fdr]


  if(direction == 'both') {
    selection_scores = selection_scores
  } else if(direction == 'up') {
    selection_scores = selection_scores[log2fc >= min_log2_fc]
  } else if(direction == 'down') {
    selection_scores = selection_scores[log2fc <= -min_log2_fc]
  }


  newobj = copy(icfObject)
  newobj[['ICFscores']] = comb_DT

  return(newobj)

}

#' @title plotICFSpot
#' @name plotICFSpot
#' @description Create barplot to visualize interaction changed features
#'
#' @param gobject giotto object
#' @param icfObject ICF (interaction changed feature) score object
#' @param source_type cell type of the source cell
#' @param source_markers markers for the source cell type
#' @param ICF_features named character vector of ICF features
#' @param cell_color_code cell color code for the interacting cell types
#' @param show_plot show plots
#' @param return_plot return plotting object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters from \code{\link{all_plots_save_function}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#'
#' @return plot
#' @export
plotICFSpot <- function(gobject,
                        icfObject,
                        source_type,
                        source_markers,
                        ICF_features,
                        cell_color_code = NULL,
                        show_plot = NA,
                        return_plot = NA,
                        save_plot = NA,
                        save_param =  list(),
                        default_save_name = 'plotICFSpot') {


  # data.table variables
  cell_type = int_cell_type = pcc_diff = NULL


  if(!'icfObject' %in% class(icfObject)) {
    stop('\n icfObject needs to be the output from findInteractionChangedFeats() or findICF() \n')
  }

  ICFscores = icfObject[['ICFscores']]

  # combine features
  names(source_markers) = rep('marker', length(source_markers))
  neighbor_types = names(ICF_features)
  all_features = c(source_markers, ICF_features)

  # warning if there are features selected that are not detected
  detected_features = unique(ICFscores[['features']])
  not_detected_features = all_features[!all_features %in% detected_features]
  if(length(not_detected_features) > 0) {
    cat('These selected features are not in the icfObject: \n',
        not_detected_features, '\n')
  }

  # data.table set column names
  features = group = NULL

  tempDT = ICFscores[features %in% all_features][cell_type == source_type][int_cell_type %in% neighbor_types]
  tempDT[, features := factor(features, levels = all_features)]
  tempDT[, group := names(all_features[all_features == features]), by = 1:nrow(tempDT)]


  if(is.null(cell_color_code)) {
    mycolors = getDistinctColors(n = length(unique(tempDT$int_cell_type)))
    names(mycolors) = unique(tempDT$int_cell_type)
  } else {
    mycolors = cell_color_code
  }


  pl = ggplot2::ggplot()
  pl = pl + ggplot2::theme_classic() + ggplot2::theme(axis.text.x = ggplot2::element_text(size = 14, angle = 45, vjust = 1, hjust = 1),
                                                      axis.text.y = ggplot2::element_text(size = 14),
                                                      axis.title = ggplot2::element_text(size = 14))
  pl = pl + ggplot2::geom_bar(data = tempDT, ggplot2::aes(x = features, y = pcc_diff, fill = int_cell_type), stat = 'identity', position = ggplot2::position_dodge())
  pl = pl + ggplot2::scale_fill_manual(values = mycolors)
  pl = pl + ggplot2::labs(x = '', title = paste0('fold-change z-scores in ' ,source_type))



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

#' @title plotCellProximityFeatSpot
#' @name plotCellProximityFeatSpot
#' @description Create visualization for cell proximity feature scores
#'
#' @param gobject giotto object
#' @param icfObject ICF (interaction changed feature) score object
#' @param method plotting method to use
#' @param min_cells minimum number of source cell type
#' @param min_cells_expr_resi Default = 0.05
#' @param min_int_cells minimum number of interacting neighbor cell type
#' @param min_int_cells_expr_resi Default = 0.05
#' @param min_pcc_diff Default = 0.05
#' @param min_fdr minimum adjusted p-value
#' @param min_zscore minimum z-score change
#' @param zscores_column calculate z-scores over cell types or features
#' @param direction differential expression directions to keep
#' @param cell_color_code vector of colors with cell types as names
#' @param show_plot show plots
#' @param return_plot return plotting object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters from \code{\link{all_plots_save_function}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#'
#' @return plot
#' @export
plotCellProximityFeatSpot = function(gobject,
                                      icfObject,
                                      method = c('volcano', 'cell_barplot', 'cell-cell', 'cell_sankey', 'heatmap', 'dotplot'),
                                      min_cells = 4,
                                      min_cells_expr_resi = 0.05,
                                      min_int_cells = 4,
                                      min_int_cells_expr_resi = 0.05,
                                      min_fdr = 0.5,
                                      min_pcc_diff = 0.05,
                                      min_zscore = 0.05,
                                      zscores_column = c('cell_type', 'features'),
                                      direction = c('both', 'up', 'down'),
                                      cell_color_code = NULL,
                                      show_plot = NA,
                                      return_plot = NA,
                                      save_plot = NA,
                                      save_param =  list(),
                                      default_save_name = 'plotCellProximityFeats') {


  if(!'icfObject' %in% class(icfObject)) {
    stop('\n icfObject needs to be the output from findInteractionChangedFeats() or findICF() \n')
  }

  # print, return and save parameters
  show_plot = ifelse(is.na(show_plot), readGiottoInstructions(gobject, param = 'show_plot'), show_plot)
  save_plot = ifelse(is.na(save_plot), readGiottoInstructions(gobject, param = 'save_plot'), save_plot)
  return_plot = ifelse(is.na(return_plot), readGiottoInstructions(gobject, param = 'return_plot'), return_plot)


  ## first filter
  filter_icf = filterICFSpot(icfObject,
                             min_cells = min_cells,
                             min_cells_expr_resi = min_cells_expr_resi,
                             min_int_cells = min_int_cells,
                             min_int_cells_expr_resi = min_int_cells_expr_resi,
                             min_fdr = min_fdr,
                             min_pcc_diff = min_pcc_diff,
                             min_zscore = min_zscore,
                             zscores_column = c('cell_type', 'features'),
                             direction = c('both', 'up', 'down'))

  complete_part = filter_icf[['ICFscores']]

  ## other parameters
  method = match.arg(method, choices = c('volcano', 'cell_barplot', 'cell-cell', 'cell_sankey', 'heatmap', 'dotplot'))


  # variables
  pcc_diff = p.adj =  unif_int = N = cell_type = int_cell_type = NULL

  ## create data.table for visualization
  if(method == 'volcano') {

    ## volcanoplot
    pl = ggplot2::ggplot()
    pl = pl + ggplot2::geom_point(data = complete_part, ggplot2::aes(x = pcc_diff, y = ifelse(is.infinite(-log10(p.adj)), 1000, -log10(p.adj))))
    pl = pl + ggplot2::theme_classic()
    pl = pl + ggplot2::geom_vline(xintercept = 0, linetype = 2)
    pl = pl + ggplot2::labs(x = 'pcc diff', y = '-log10(p.adjusted)')


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


  } else if(method == 'cell-cell') {

    nr_int_selection_scores = complete_part[, .N, by = unif_int]
    order_interactions = nr_int_selection_scores[order(N)]$unif_int

    complete_part[, unif_int := factor(unif_int, order_interactions)]

    pl <- ggplot2::ggplot()
    pl <- pl + ggplot2::geom_bar(data = complete_part, ggplot2::aes(x = unif_int, fill = unif_int))
    pl <- pl + ggplot2::theme_classic() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 1))
    pl <- pl + ggplot2::coord_flip()

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


  } else if(method == 'cell_barplot') {


    # by source cell type plot
    nr_source_selection_scores = complete_part[, .N, by = cell_type]
    order_source = nr_source_selection_scores[order(N)]$cell_type

    complete_part[, cell_type := factor(cell_type, order_source)]

    pl <- ggplot2::ggplot()
    pl <- pl + ggplot2::geom_bar(data = complete_part, ggplot2::aes(x = cell_type, fill = int_cell_type))
    if(!is.null(cell_color_code)) {
      pl <- pl + ggplot2::scale_fill_manual(values = cell_color_code)
    }
    pl <- pl + ggplot2::theme_classic() + ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
    pl <- pl + ggplot2::labs(x = '', y = '# of features influenced by cell neighborhood')


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

  } else if(method == 'cell_sankey') {


    # package check for ggalluvial
    package_check(pkg_name = 'ggalluvial', repository = 'CRAN')


    testalluv = complete_part[, .N, by = c('int_cell_type', 'cell_type')]

    # library(ggalluvial) # this is needed for it to work, why??
    # maybe use requireNamespace() instead?

    pl <- ggplot2::ggplot(testalluv,
                          ggplot2::aes(y = N, axis1 = cell_type, axis2 = int_cell_type)) +
      ggalluvial::geom_alluvium(aes(fill = cell_type), width = 1/12) +
      ggalluvial::geom_stratum(width = 1/12, fill = "black", color = "grey") +
      ggplot2::scale_x_discrete(limits = c("cell type", "neighbours"), expand = c(.05, .05)) +
      ggplot2::geom_label(stat = "stratum", label.strata = TRUE, size = 3) +
      ggplot2::theme_classic() + ggplot2::labs(x = '', y = '# of features influenced by cell neighborhood')

    if(!is.null(cell_color_code)) {
      pl <- pl + ggplot2::scale_fill_manual(values = cell_color_code)
    }



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

  } else if(method == 'dotplot') {

    changed_features = complete_part[, .N, by = c('cell_type', 'int_cell_type')]

    changed_features[, cell_type := factor(cell_type, unique(cell_type))]
    changed_features[, int_cell_type := factor(int_cell_type, unique(int_cell_type))]

    pl = ggplot2::ggplot()
    pl = pl + ggplot2::theme_classic()
    pl = pl + ggplot2::geom_point(data = changed_features, ggplot2::aes(x = cell_type, y = int_cell_type, size = N))
    pl = pl + ggplot2::scale_size_continuous(guide=guide_legend(title = '# of ICFs'))
    pl = pl + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 1, hjust = 1))
    pl = pl + ggplot2::labs(x = 'source cell type', y = 'neighbor cell type')

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

  } else if(method == 'heatmap') {

    changed_features = complete_part[, .N, by = c('cell_type', 'int_cell_type')]

    changed_features[, cell_type := factor(cell_type, unique(cell_type))]
    changed_features[, int_cell_type := factor(int_cell_type, unique(int_cell_type))]

    changed_features_d = data.table::dcast.data.table(changed_features, cell_type~int_cell_type, value.var = 'N', fill = 0)
    changed_features_m = dt_to_matrix(changed_features_d)

    col_fun = circlize::colorRamp2(breaks = stats::quantile(log2(changed_features_m+1)),
                                   colors =  c("white", 'white', "blue", "yellow", "red"))

    heatm = ComplexHeatmap::Heatmap(as.matrix(log2(changed_features_m+1)), col = col_fun,
                                    row_title = 'cell_type', column_title = 'int_cell_type', heatmap_legend_param = list(title = 'log2(# DEGs)'))

    ## print plot
    if(show_plot == TRUE) {
      print(heatm)
    }

    ## save plot
    if(save_plot == TRUE) {
      do.call('all_plots_save_function', c(list(gobject = gobject, plot_object = heatm, default_save_name = default_save_name), save_param))
    }

    ## return plot
    if(return_plot == TRUE) {
      return(heatm)
    }

  }
}


# * ####
# cell communication spots ####


#' @title specific_CCCScores_spots
#' @name specific_CCCScores_spots
#' @description Specific Cell-Cell communication scores based on spatial expression of interacting cells at spots resulation
#'
#' @param gobject giotto object to use
#' @param spat_unit spatial unit (e.g. 'cell')
#' @param feat_type feature type (e.g. 'rna')
#' @param expr_residual spatial network to use for identifying interacting cells
#' @param dwls_values dwls matrix
#' @param proximityMat cell cell communication score matrix
#' @param random_iter number of iterations
#' @param cell_type_1 first cell type
#' @param cell_type_2 second cell type
#' @param feature_set_1 first specific feature set from feature pairs
#' @param feature_set_2 second specific feature set from feature pairs
#' @param min_observations minimum number of interactions needed to be considered
#' @param detailed provide more detailed information (random variance and z-score)
#' @param adjust_method which method to adjust p-values
#' @param adjust_target adjust multiple hypotheses at the cell or feature level
#' @param set_seed set a seed for reproducibility
#' @param seed_number seed number
#' @param verbose verbose
#'
#' @return Cell-Cell communication scores for feature pairs based on spatial interaction
#' @details Statistical framework to identify if pairs of features (such as ligand-receptor combinations)
#' are expressed at higher levels than expected based on a reshuffled null distribution
#' of feature expression values in cells that are spatially in proximity to eachother.
#' \itemize{
#'  \item{LR_comb:}{Pair of ligand and receptor}
#'  \item{lig_cell_type:}{ cell type to assess expression level of ligand }
#'  \item{lig_expr:}{ average expressionresidual(observed - DWLS_predicted) of ligand in lig_cell_type }
#'  \item{ligand:}{ ligand name }
#'  \item{rec_cell_type:}{ cell type to assess expression level of receptor }
#'  \item{rec_expr:}{ average expression residual(observed - DWLS_predicted) of receptor in rec_cell_type}
#'  \item{receptor:}{ receptor name }
#'  \item{LR_expr:}{ combined average ligand and receptor expression }
#'  \item{lig_nr:}{ total number of cells from lig_cell_type that spatially interact with cells from rec_cell_type }
#'  \item{rec_nr:}{ total number of cells from rec_cell_type that spatially interact with cells from lig_cell_type }
#'  \item{rand_expr:}{ average combined ligand and receptor expression residual from random spatial permutations }
#'  \item{av_diff:}{ average difference between LR_expr and rand_expr over all random spatial permutations }
#'  \item{sd_diff:}{ (optional) standard deviation of the difference between LR_expr and rand_expr over all random spatial permutations }
#'  \item{z_score:}{ (optinal) z-score }
#'  \item{log2fc:}{ LR_expr - rand_expr }
#'  \item{pvalue:}{ p-value }
#'  \item{LR_cell_comb:}{ cell type pair combination }
#'  \item{p.adj:}{ adjusted p-value }
#'  \item{PI:}{ significanc score: log2fc * -log10(p.adj) }
#' }
#' @keywords internal
specific_CCCScores_spots = function(gobject,
                                    spat_unit = NULL,
                                    feat_type = NULL,
                                    expr_residual,
                                    dwls_values,
                                    proximityMat,
                                    random_iter = 1000,
                                    cell_type_1 = 'astrocytes',
                                    cell_type_2 = 'endothelial',
                                    feature_set_1,
                                    feature_set_2,
                                    min_observations = 2,
                                    detailed = FALSE,
                                    adjust_method = c('fdr', 'bonferroni',' BH', 'holm', 'hochberg', 'hommel',
                                                      'BY','none'),
                                    adjust_target = c('features', 'cells'),
                                    set_seed = FALSE,
                                    seed_number = 1234,
                                    verbose = T){

  # data.table variables
  from_to = cell_ID = lig_cell_type = rec_cell_type = lig_nr = rec_nr = rand_expr = NULL
  av_diff = log2fc = LR_expr = pvalue = LR_cell_comb = p.adj = LR_comb = PI = NULL
  sd_diff = z_score = ligand = receptor = NULL

  # get parameters
  adjust_method = match.arg(adjust_method, choices = c("fdr", "bonferroni","BH", "holm", "hochberg", "hommel",
                                                       "BY", "none"))
  adjust_target = match.arg(adjust_target, choices = c('features', 'cells'))

  # select cell_ids with cell-types
  cell_direction_1 = paste0(cell_type_1,'--',cell_type_2)
  cell_direction_2 = paste0(cell_type_2,'--',cell_type_1)

  print(paste0('Processing specific CCC Scores: ', cell_direction_1))

  proxi_1 = proximityMat[cell_direction_1,]
  proxi_2 = proximityMat[cell_direction_2,]

  ct1_cell_ids = names(proxi_1[proxi_1 > 0])
  ct2_cell_ids = names(proxi_2[proxi_2 > 0])

  # dwls value for cell types
  dwls_ct1 = dwls_values[,cell_type_1]
  dwls_ct2 = dwls_values[,cell_type_2]

  # make sure that there are sufficient observations
  if(length(ct1_cell_ids) <= min_observations | length(ct2_cell_ids) <= min_observations) {

    return(NULL)

  } else {

    # get feature expression residual for ligand and receptor
    expr_res_L = expr_residual[feature_set_1, ct1_cell_ids]
    expr_res_R = expr_residual[feature_set_2, ct2_cell_ids]

    # compute Ligand value
    lig_expr = t(t(expr_res_L) * dwls_ct1[ct1_cell_ids])
    rec_expr = t(t(expr_res_R) * dwls_ct2[ct2_cell_ids])

    lig_expr = round(rowMeans(lig_expr), 7)
    rec_expr = round(rowMeans(rec_expr), 7)

    comScore = data.table::data.table(LR_comb = paste0(feature_set_1, '-', feature_set_2),
                                      lig_cell_type = rep(cell_type_1, length(feature_set_1)),
                                      lig_expr = lig_expr,
                                      ligand = feature_set_1,
                                      rec_cell_type = rep(cell_type_2, length(feature_set_2)),
                                      rec_expr = rec_expr,
                                      receptor = feature_set_2,
                                      lig_nr = rep(length(ct1_cell_ids), length(feature_set_1)),
                                      rec_nr = rep(length(ct2_cell_ids), length(feature_set_1))
    )

    comScore[, LR_expr := lig_expr + rec_expr]
    comScore = comScore[, .(LR_comb, lig_cell_type, lig_expr, ligand,
                            rec_cell_type, rec_expr, receptor, LR_expr, lig_nr, rec_nr)]

    # prepare for randomized scores
    total_av = rep(0, nrow(comScore))

    if(detailed == FALSE) {
      total_sum = rep(0, nrow(comScore))
    } else {
      total_sum = matrix(nrow = nrow(comScore), ncol = random_iter)
    }

    total_bool = rep(0, nrow(comScore))

    # all_cell_ids = pDataDT(gobject = gobject,
    #                        spat_unit = spat_unit,
    #                        feat_type = feat_type)$cell_ID
    
    all_cell_ids = colnames(expr_residual)
    
    ## simulations ##
    for(sim in 1:random_iter) {

      if(verbose == TRUE) cat('simulation ', sim, '\n')

      # get random ids and subset
      if(set_seed == TRUE) {
        seed_number = seed_number+sim
        set.seed(seed = seed_number)
      }
      #random_ids_1 = all_cell_ids[sample(length(all_cell_ids), size = length(ct1_cell_ids))]
      #random_ids_2 = all_cell_ids[sample(length(all_cell_ids), size = length(ct2_cell_ids))]
      
      random_ids_1 = sample(all_cell_ids, size = length(ct1_cell_ids), replace = FALSE)
      random_ids_2 = sample(all_cell_ids, size = length(ct2_cell_ids), replace = FALSE)
      
      # get feature expression residual for ligand and receptor
      random_expr_res_L = expr_residual[feature_set_1, random_ids_1]
      random_expr_res_R = expr_residual[feature_set_2, random_ids_2]

      # compute Ligand value
      random_lig_expr = t(t(random_expr_res_L) * dwls_ct1[random_ids_1])
      random_rec_expr = t(t(random_expr_res_R) * dwls_ct2[random_ids_2])

      random_lig_expr = round(rowMeans(random_lig_expr), 7)
      random_rec_expr = round(rowMeans(random_rec_expr), 7)

      randomScore = data.table::data.table(lig_expr = random_lig_expr,
                                           rec_expr = random_rec_expr)
      randomScore = randomScore[, LR_expr := lig_expr + rec_expr]

      # average random score
      total_av = total_av + randomScore[['LR_expr']]

      # difference between observed and random
      difference = comScore[['LR_expr']] - randomScore[['LR_expr']]

      # calculate total difference
      if(detailed == FALSE) {
        total_sum = total_sum+difference
      } else {
        total_sum[,sim] = difference
      }

      # calculate p-values
      difference[difference > 0] = 1
      difference[difference < 0] = -1
      total_bool = total_bool + difference

    }

    comScore[, rand_expr := total_av/random_iter]

    if(detailed == TRUE) {
      av_difference_scores = rowMeans_flex(total_sum)
      sd_difference_scores = apply(total_sum, MARGIN = 1, FUN = stats::sd)

      comScore[, av_diff := av_difference_scores]
      comScore[, sd_diff := sd_difference_scores]
      comScore[, z_score := (LR_expr - rand_expr)/sd_diff]

    } else {
      comScore[, av_diff := total_sum/random_iter]
    }

    comScore[, log2fc := LR_expr - rand_expr]
    comScore[, pvalue := total_bool/random_iter]
    comScore[, pvalue := ifelse(pvalue > 0, 1-pvalue, 1+pvalue)]
    comScore[, LR_cell_comb := paste0(lig_cell_type,'--',rec_cell_type)]

    if(adjust_target == 'features') {
      comScore[, p.adj := stats::p.adjust(pvalue, method = adjust_method), by = .(LR_cell_comb)]
    } else if(adjust_target == 'cells'){
      comScore[, p.adj := stats::p.adjust(pvalue, method = adjust_method), by = .(LR_comb)]
    }

    # get minimum adjusted p.value that is not zero
    all_p.adj = comScore[['p.adj']]
    lowest_p.adj = min(all_p.adj[all_p.adj != 0])
    comScore[, PI := ifelse(p.adj == 0, log2fc*(-log10(lowest_p.adj)), log2fc*(-log10(p.adj)))]

    return(comScore)

  }
}



#' @title spatCellCellcomSpots
#' @name spatCellCellcomSpots
#' @description Spatial Cell-Cell communication scores based on spatial expression of interacting cells at spots resolution
#'
#' @param gobject giotto object to use
#' @param spat_unit spatial unit (e.g. 'cell')
#' @param feat_type feature type (e.g. 'rna')
#' @param ave_celltype_exp Matrix with average expression per cell type
#' @param expression_values (e.g. 'normalized', 'scaled', 'custom')
#' @param spatial_network_name spatial network to use for identifying interacting cells
#' @param cluster_column cluster column with cell type information
#' @param random_iter number of iterations
#' @param feature_set_1 first specific feature set from feature pairs
#' @param feature_set_2 second specific feature set from feature pairs
#' @param min_observations minimum number of interactions needed to be considered
#' @param detailed provide more detailed information (random variance and z-score)
#' @param adjust_method which method to adjust p-values
#' @param adjust_target adjust multiple hypotheses at the cell or feature level
#' @param do_parallel run calculations in parallel with mclapply
#' @param cores number of cores to use if do_parallel = TRUE
#' @param set_seed set a seed for reproducibility
#' @param seed_number seed number
#' @param verbose verbose (e.g. 'a little', 'a lot', 'none')
#'
#' @return Cell-Cell communication scores for feature pairs based on spatial interaction
#' @details Statistical framework to identify if pairs of features (such as ligand-receptor combinations)
#' are expressed at higher levels than expected based on a reshuffled null distribution
#' of feature expression values in cells that are spatially in proximity to eachother..
#' \itemize{
#'  \item{LR_comb:}{Pair of ligand and receptor}
#'  \item{lig_cell_type:}{ cell type to assess expression level of ligand }
#'  \item{lig_expr:}{ average expression residual(observed - DWLS_predicted) of ligand in lig_cell_type }
#'  \item{ligand:}{ ligand name }
#'  \item{rec_cell_type:}{ cell type to assess expression level of receptor }
#'  \item{rec_expr:}{ average expression residual(observed - DWLS_predicted) of receptor in rec_cell_type}
#'  \item{receptor:}{ receptor name }
#'  \item{LR_expr:}{ combined average ligand and receptor expression residual}
#'  \item{lig_nr:}{ total number of cells from lig_cell_type that spatially interact with cells from rec_cell_type }
#'  \item{rec_nr:}{ total number of cells from rec_cell_type that spatially interact with cells from lig_cell_type }
#'  \item{rand_expr:}{ average combined ligand and receptor expression residual from random spatial permutations }
#'  \item{av_diff:}{ average difference between LR_expr and rand_expr over all random spatial permutations }
#'  \item{sd_diff:}{ (optional) standard deviation of the difference between LR_expr and rand_expr over all random spatial permutations }
#'  \item{z_score:}{ (optinal) z-score }
#'  \item{log2fc:}{ LR_expr - rand_expr }
#'  \item{pvalue:}{ p-value }
#'  \item{LR_cell_comb:}{ cell type pair combination }
#'  \item{p.adj:}{ adjusted p-value }
#'  \item{PI:}{ significanc score: log2fc * -log10(p.adj) }
#' }
#' @export
spatCellCellcomSpots = function(gobject,
                                spat_unit = NULL,
                                feat_type = NULL,
                                ave_celltype_exp,
                                spatial_network_name = 'Delaunay_network',
                                cluster_column = 'cell_ID',
                                random_iter = 1000,
                                feature_set_1,
                                feature_set_2,
                                min_observations = 2,
                                expression_values = c('normalized', 'scaled', 'custom'),
                                detailed = FALSE,
                                adjust_method = c('fdr', 'bonferroni', 'BH', 'holm', 'hochberg', 'hommel',
                                                  'BY', 'none'),
                                adjust_target = c('features', 'cells'),
                                do_parallel = TRUE,
                                cores = NA,
                                set_seed = TRUE,
                                seed_number = 1234,
                                verbose = c('a little', 'a lot', 'none')){
  
  # data.table vars
  V1 = V2 = LR_cell_comb = NULL

  # code start
  verbose = match.arg(verbose, choices = c('a little', 'a lot', 'none'))

  ## check if spatial network exists ##
  spat_networks = names(gobject@spatial_network[[spat_unit]])
  
  if(!spatial_network_name %in% spat_networks) {
    stop(spatial_network_name, ' is not an existing spatial network \n',
         'use showGiottoSpatNetworks() to see the available networks \n',
         'or create a new spatial network with createSpatialNetwork() \n')
  }


  # expression data
  values = match.arg(expression_values, choices = c('normalized', 'scaled', 'custom'))
  expr_residual = cal_expr_residual(gobject = gobject,
                                    spat_unit = spat_unit,
                                    feat_type = feat_type,
                                    ave_celltype_exp = ave_celltype_exp)

  # compute cell proximity for each spot
  proximityMat = cellProximityEnrichmentEachSpot(gobject = gobject,
                                                 spat_unit = spat_unit,
                                                 feat_type = feat_type,
                                                 spatial_network_name = spatial_network_name)

  # select overlapped spots
  intersect_cell_IDs = intersect(colnames(expr_residual), colnames(proximityMat))
  expr_residual = expr_residual[, intersect_cell_IDs]
  proximityMat = proximityMat[, intersect_cell_IDs]

  # exact spatial_enrichment matrix
  dwls_values = get_spatial_enrichment(gobject = gobject,
                                       spat_unit = spat_unit,
                                       feat_type = feat_type,
                                       output = 'data.table')
  data.table::setDF(dwls_values)
  rownames_dwls = dwls_values[,'cell_ID']
  dwls_values = as.matrix(dwls_values[,-1])
  rownames(dwls_values) = rownames_dwls

  # check feature list
  LR_comb = data.table::data.table(ligand = feature_set_1, receptor = feature_set_2)
  # check LR pair not captured in giotto object
  LR_out = LR_comb[!LR_comb$ligand %in% rownames(expr_residual) | !LR_comb$receptor %in% rownames(expr_residual)]

  if (dim(LR_out)[1] > 0){
    print('Ligand or receptor were removed after computing expresion residual.')
    print(LR_out)
    LR_comb = LR_comb[LR_comb$ligand %in% rownames(expr_residual) &  LR_comb$receptor %in% rownames(expr_residual) ]
    feature_set_1 = LR_comb$ligand
    feature_set_2 = LR_comb$receptor
  }

  ## get all combinations between cell types
  combn_DT = data.table::data.table(LR_cell_comb = rownames(proximityMat))
  combn_DT[, V1 := strsplit(LR_cell_comb, '--')[[1]][1], by = 1:nrow(combn_DT)]
  combn_DT[, V2 := strsplit(LR_cell_comb, '--')[[1]][2], by = 1:nrow(combn_DT)]

  ## parallel option ##
  if(do_parallel == TRUE) {


    savelist = giotto_lapply(X = 1:nrow(combn_DT), cores = cores, fun = function(row) {

      cell_type_1 = combn_DT[row][['V1']]
      cell_type_2 = combn_DT[row][['V2']]


      specific_scores = specific_CCCScores_spots(gobject = gobject,
                                                 spat_unit = spat_unit,
                                                 feat_type = feat_type,
                                                 expr_residual = expr_residual,
                                                 dwls_values = dwls_values,
                                                 proximityMat = proximityMat,
                                                 random_iter = random_iter,
                                                 cell_type_1 = cell_type_1,
                                                 cell_type_2 = cell_type_2,
                                                 feature_set_1 = feature_set_1,
                                                 feature_set_2 = feature_set_2,
                                                 min_observations = min_observations,
                                                 detailed = detailed,
                                                 adjust_method = adjust_method,
                                                 adjust_target = adjust_target,
                                                 set_seed = set_seed,
                                                 seed_number = seed_number)
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

      specific_scores = specific_CCCScores_spots(gobject = gobject,
                                                 spat_unit = spat_unit,
                                                 feat_type = feat_type,
                                                 expr_residual = expr_residual,
                                                 dwls_values = dwls_values,
                                                 proximityMat = proximityMat,
                                                 random_iter = random_iter,
                                                 cell_type_1 = cell_type_1,
                                                 cell_type_2 = cell_type_2,
                                                 feature_set_1 = feature_set_1,
                                                 feature_set_2 = feature_set_2,
                                                 min_observations = min_observations,
                                                 detailed = detailed,
                                                 adjust_method = adjust_method,
                                                 adjust_target = adjust_target,
                                                 set_seed = set_seed,
                                                 seed_number = seed_number,
                                                 verbose = specific_verbose)
      savelist[[row]] = specific_scores
      countdown = countdown - 1
    }

  }

  finalDT = do.call('rbind', savelist)

  # data.table variables
  LR_comb = LR_expr = NULL

  data.table::setorder(finalDT, LR_comb, -LR_expr)

  return(finalDT)

}








