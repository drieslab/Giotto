

## spatial gene detection ####

#' @title spat_fish_func
#' @name spat_fish_func
#' @description performs fisher exact test
#' @keywords internal
spat_fish_func = function(gene,
                          bin_matrix,
                          spat_mat,
                          calc_hub = F,
                          hub_min_int = 3) {

  gene_vector = bin_matrix[rownames(bin_matrix) == gene,]

  gene_vectorA = gene_vector[names(gene_vector) %in% rownames(spat_mat)]
  gene_vectorA = gene_vectorA[match(rownames(spat_mat), names(gene_vectorA))]

  gene_vectorB = gene_vector[names(gene_vector) %in% colnames(spat_mat)]
  gene_vectorB = gene_vectorB[match(colnames(spat_mat), names(gene_vectorB))]

  test1 = spat_mat*gene_vectorA
  test2 = t_giotto(t_giotto(spat_mat)*gene_vectorB)

  sourcevalues = test1[spat_mat == 1]
  targetvalues = test2[spat_mat == 1]

  # option 1
  test = paste0(sourcevalues,'-',targetvalues)


  if(length(unique(test)) < 4) {

    possibs = c("1-1","0-1","1-0","0-0")
    missings_possibs = possibs[!possibs %in% unique(test)]
    test = c(test, missings_possibs)

    table_test = table(test)
    table_test[names(table_test) %in% missings_possibs] = 0
    table_matrix = matrix(table_test, byrow = T, nrow = 2)

  } else {
    table_matrix = matrix(table(test), byrow = T, nrow = 2)
  }

  if(calc_hub == TRUE) {
    high_cells = names(gene_vector[gene_vector == 1])
    subset_spat_mat = spat_mat[rownames(spat_mat) %in% high_cells, colnames(spat_mat) %in% high_cells]

    if(length(subset_spat_mat) == 1) {
      hub_nr = 0
    } else {
      subset_spat_mat = spat_mat[rownames(spat_mat) %in% high_cells, colnames(spat_mat) %in% high_cells]
      rowhubs = rowSums_giotto(subset_spat_mat)
      colhubs = colSums_giotto(subset_spat_mat)
      hub_nr = length(unique(c(names(colhubs[colhubs > hub_min_int]), names(rowhubs[colhubs > hub_min_int]))))
    }

    fish_res = stats::fisher.test(table_matrix)[c('p.value','estimate')]
    return(c(genes = list(gene), fish_res, hubs = list(hub_nr)))

  } else {

    fish_res = stats::fisher.test(table_matrix)[c('p.value','estimate')]
    return(c(genes = list(gene), fish_res))
  }

}

#' @title spat_fish_func_DT
#' @name spat_fish_func_DT
#' @description performs fisher exact test with data.table implementation
#' @keywords internal
spat_fish_func_DT = function(bin_matrix_DTm,
                             spat_netw_min,
                             calc_hub = F,
                             hub_min_int = 3,
                             cores = NA) {

  # set number of cores automatically, but with limit of 10
  cores = determine_cores(cores)
  data.table::setDTthreads(threads = cores)

  # data.table variables
  from_value = to_value = gene_ID = N = to = from = cell_ID = V1 = NULL

  # get binarized expression values for the neighbors
  spatial_network_min_ext = data.table::merge.data.table(spat_netw_min, bin_matrix_DTm, by.x = 'from', by.y = 'variable', allow.cartesian = T)
  data.table::setnames(spatial_network_min_ext, 'value', 'from_value')

  spatial_network_min_ext = data.table::merge.data.table(spatial_network_min_ext, by.x = c('to', 'gene_ID'), bin_matrix_DTm, by.y = c('variable', 'gene_ID'))
  data.table::setnames(spatial_network_min_ext, 'value', 'to_value')


  # summarize the different combinations
  spatial_network_min_ext[, combn := paste0(from_value,'-',to_value)]
  freq_summary = spatial_network_min_ext[, .N, by = .(gene_ID, combn)]
  data.table::setorder(freq_summary, gene_ID, combn)

  genes = unique(freq_summary$gene_ID)
  all_combn = c('0-0', '0-1', '1-0', '1-1')

  # create a zeroes DT to add missing observations
  freq_summary_zeroes = data.table::data.table(gene_ID = rep(genes, each = 4),
                                   combn = rep(all_combn, length(genes)),
                                   N = 0)
  freq_summary2 = rbind(freq_summary, freq_summary_zeroes)
  freq_summary2[, N := sum(N), by = .(gene_ID, combn)]
  freq_summary2 = unique(freq_summary2)

  # sort the combinations and run fisher test
  data.table::setorder(freq_summary2, gene_ID, combn, -N)
  fish_results = freq_summary2[, stats::fisher.test(matrix(N, nrow = 2))[c(1,3)], by = gene_ID]


  ## hubs ##
  if(calc_hub == TRUE) {

    double_pos = spatial_network_min_ext[combn == '1-1']

    double_pos_to = double_pos[, .N, by = .(gene_ID, to)]
    data.table::setnames(double_pos_to, 'to', 'cell_ID')
    double_pos_from = double_pos[, .N, by = .(gene_ID, from)]
    data.table::setnames(double_pos_from, 'from', 'cell_ID')

    double_pos_both = rbind(double_pos_to, double_pos_from)
    double_pos_both = double_pos_both[, sum(N), by = .(gene_ID, cell_ID)]
    data.table::setorder(double_pos_both, gene_ID, -V1)

    # get hubs and add 0's
    hub_DT = double_pos_both[V1 > hub_min_int, .N, by = gene_ID]
    hub_DT_zeroes = data.table::data.table(gene_ID = unique(spatial_network_min_ext$gene_ID), N = 0)
    hub_DT2 = rbind(hub_DT, hub_DT_zeroes)

    hub_DT2 = hub_DT2[, sum(N), by = gene_ID]
    data.table::setnames(hub_DT2, 'V1', 'hub_nr')

    fish_results = data.table::merge.data.table(fish_results, hub_DT2, by = 'gene_ID')

  }

  return(fish_results)

}



#' @title spat_OR_func
#' @name spat_OR_func
#' @description calculate odds-ratio
#' @keywords internal
spat_OR_func = function(gene,
                        bin_matrix,
                        spat_mat,
                        calc_hub = F,
                        hub_min_int = 3) {

  gene_vector = bin_matrix[rownames(bin_matrix) == gene,]

  gene_vectorA = gene_vector[names(gene_vector) %in% rownames(spat_mat)]
  gene_vectorA = gene_vectorA[match(rownames(spat_mat), names(gene_vectorA))]

  gene_vectorB = gene_vector[names(gene_vector) %in% colnames(spat_mat)]
  gene_vectorB = gene_vectorB[match(colnames(spat_mat), names(gene_vectorB))]

  test1 = spat_mat*gene_vectorA
  test2 = t_giotto(t_giotto(spat_mat)*gene_vectorB)

  sourcevalues = test1[spat_mat == 1]
  targetvalues = test2[spat_mat == 1]

  # option 1
  test = paste0(sourcevalues,'-',targetvalues)


  if(length(unique(test)) < 4) {

    possibs = c("1-1","0-1","1-0","0-0")
    missings_possibs = possibs[!possibs %in% unique(test)]
    test = c(test, missings_possibs)

    table_test = table(test)
    table_test[names(table_test) %in% missings_possibs] = 0
    table_matrix = matrix(table_test, byrow = T, nrow = 2)

  } else {
    table_matrix = matrix(table(test), byrow = T, nrow = 2)
  }


  if(calc_hub == TRUE) {
    high_cells = names(gene_vector[gene_vector == 1])
    subset_spat_mat = spat_mat[rownames(spat_mat) %in% high_cells, colnames(spat_mat) %in% high_cells]

    if(length(subset_spat_mat) == 1) {
      hub_nr = 0
    } else {
      rowhubs = rowSums_giotto(subset_spat_mat)
      colhubs = colSums_giotto(subset_spat_mat)
      hub_nr = length(unique(c(names(colhubs[colhubs > hub_min_int]), names(rowhubs[colhubs > hub_min_int]))))
    }

    fish_matrix = table_matrix
    fish_matrix = fish_matrix/1000
    OR = ((fish_matrix[1]*fish_matrix[4]) / (fish_matrix[2]*fish_matrix[3]))

    return(c(genes = list(gene), OR, hubs = list(hub_nr)))

  }

  fish_matrix = table_matrix
  fish_matrix = fish_matrix/1000
  OR = ((fish_matrix[1]*fish_matrix[4]) / (fish_matrix[2]*fish_matrix[3]))
  return(c(genes = list(gene), OR))

}


#' @title spat_OR_func_DT
#' @name spat_OR_func_DT
#' @description calculate odds-ratio with data.table implementation
#' @keywords internal
spat_OR_func_DT = function(bin_matrix_DTm,
                           spat_netw_min,
                           calc_hub = F,
                           hub_min_int = 3,
                           cores = NA) {

  # set number of cores automatically, but with limit of 10
  cores = determine_cores(cores)
  data.table::setDTthreads(threads = cores)

  # data.table variables
  from_value = to_value = gene_ID = N = to = from = cell_ID = V1 = NULL

  # get binarized expression values for the neighbors
  spatial_network_min_ext = data.table::merge.data.table(spat_netw_min, bin_matrix_DTm, by.x = 'from', by.y = 'variable', allow.cartesian = T)
  data.table::setnames(spatial_network_min_ext, 'value', 'from_value')

  spatial_network_min_ext = data.table::merge.data.table(spatial_network_min_ext, by.x = c('to', 'gene_ID'), bin_matrix_DTm, by.y = c('variable', 'gene_ID'))
  data.table::setnames(spatial_network_min_ext, 'value', 'to_value')


  # summarize the different combinations
  spatial_network_min_ext[, combn := paste0(from_value,'-',to_value)]
  freq_summary = spatial_network_min_ext[, .N, by = .(gene_ID, combn)]
  data.table::setorder(freq_summary, gene_ID, combn)

  genes = unique(freq_summary$gene_ID)
  all_combn = c('0-0', '0-1', '1-0', '1-1')

  # create a zeroes DT to add missing observations
  freq_summary_zeroes = data.table::data.table(gene_ID = rep(genes, each = 4),
                                   combn = rep(all_combn, length(genes)),
                                   N = 0)
  freq_summary2 = rbind(freq_summary, freq_summary_zeroes)
  freq_summary2[, N := sum(N), by = .(gene_ID, combn)]
  freq_summary2 = unique(freq_summary2)

  # sort the combinations and run fisher test
  setorder(freq_summary2, gene_ID, combn, -N)
  or_results = freq_summary2[, OR_test_fnc(matrix(N, nrow = 2)), by = gene_ID]


  ## hubs ##
  if(calc_hub == TRUE) {

    double_pos = spatial_network_min_ext[combn == '1-1']

    double_pos_to = double_pos[, .N, by = .(gene_ID, to)]
    data.table::setnames(double_pos_to, 'to', 'cell_ID')
    double_pos_from = double_pos[, .N, by = .(gene_ID, from)]
    data.table::setnames(double_pos_from, 'from', 'cell_ID')

    double_pos_both = rbind(double_pos_to, double_pos_from)
    double_pos_both = double_pos_both[, sum(N), by = .(gene_ID, cell_ID)]
    data.table::setorder(double_pos_both, gene_ID, -V1)

    # get hubs and add 0's
    hub_DT = double_pos_both[V1 > hub_min_int, .N, by = gene_ID]
    hub_DT_zeroes = data.table::data.table(gene_ID = unique(spatial_network_min_ext$gene_ID), N = 0)
    hub_DT2 = rbind(hub_DT, hub_DT_zeroes)

    hub_DT2 = hub_DT2[, sum(N), by = gene_ID]
    data.table::setnames(hub_DT2, 'V1', 'hub_nr')

    or_results = data.table::merge.data.table(or_results, hub_DT2, by = 'gene_ID')

  }

  return(or_results)

}


#' @title OR_test_fnc
#' @name OR_test_fnc
#' @description calculate odds-ratio from a 2x2 matrix
#' @keywords internal
OR_test_fnc = function(matrix) {
  OR = ((matrix[1]*matrix[4]) / (matrix[2]*matrix[3]))
  list('estimate' = OR)
}


#' @title calc_spatial_enrichment_minimum
#' @name calc_spatial_enrichment_minimum
#' @description calculate spatial enrichment using a simple and efficient for loop
#' @keywords internal
calc_spatial_enrichment_minimum = function(spatial_network,
                                           bin_matrix,
                                           adjust_method = 'fdr',
                                           do_fisher_test = TRUE) {

  # data.table variables
  from = to = genes = variable = value = p.value = adj.p.value = score = estimate = NULL

  spatial_network_min = spatial_network[,.(from, to)]

  all_colindex = 1:ncol(bin_matrix)
  names(all_colindex) = colnames(bin_matrix)

  # code for possible combinations
  convert_code = c(1, 2, 3, 4)
  names(convert_code) = c('0-0', '0-1', '1-0', '1-1')

  # preallocate final matrix for results
  matrix_res = matrix(data = NA, nrow = nrow(bin_matrix), ncol = nrow(spatial_network_min))

  ## 1. summarize results for each edge in the network
  for(row_i in 1:nrow(spatial_network_min)) {

    from_id = spatial_network_min[row_i][['from']]
    to_id = spatial_network_min[row_i][['to']]

    sumres = data.table::as.data.table(bin_matrix[, all_colindex[c(from_id, to_id)]])
    sumres[, combn := paste0(get(from_id),'-',get(to_id))]

    ## maybe a slightly faster alternative ##
    #sumres[, sum := get(from_id)+get(to_id)]
    #sumres[, combn := ifelse(sum == 0, 1,
    #                          ifelse(sum == 2, 4,
    #                                 ifelse(get(from_id) == 1, 3, 2)))]
    #code_res = sumres[['combn']]

    code_res = convert_code[sumres$combn]
    matrix_res[, row_i] = code_res
  }

  rownames(matrix_res) = rownames(bin_matrix)


  # preallocate matrix for table results
  table_res = matrix(data = NA, nrow(matrix_res), ncol = 4)

  ## 2. calculate the frequencies of possible combinations ##
  # '0-0' = 1, '0-1' = 2, '1-0' = 3 and '1-1' = 4
  for(row_i in 1:nrow(matrix_res)) {

    x = matrix_res[row_i,]
    x = factor(x, levels = c(1,2,3,4))
    tabres = as.vector(table(x))

    table_res[row_i,] = tabres
  }

  rownames(table_res) = rownames(matrix_res)
  colnames(table_res) = 1:4

  rable_resDT = data.table::as.data.table(table_res)
  rable_resDT[, genes := rownames(table_res)]

  rable_resDTm = data.table::melt.data.table(rable_resDT, id.vars = 'genes')
  data.table::setorder(rable_resDTm, genes, variable)

  ## run fisher test ##
  if(do_fisher_test == TRUE) {
    results = rable_resDTm[, stats::fisher.test(matrix(value, nrow = 2))[c(1,3)], by = genes]

    # replace zero p-values with lowest p-value
    min_pvalue = min(results$p.value[results$p.value > 0])
    results[, p.value := ifelse(p.value == 0, min_pvalue, p.value)]
    results[, adj.p.value := stats::p.adjust(p.value, method = adjust_method)]

    # sort genes based on p-value and estimate
    results[, score := -log(p.value) * estimate]
    data.table::setorder(results, -score)


  } else {

    results = rable_resDTm[,  OR_test_fnc(matrix(value, nrow = 2)), by = genes]
    data.table::setorder(results, -estimate)

  }

  return(results)

}

#' @title calc_spatial_enrichment_matrix
#' @name calc_spatial_enrichment_matrix
#' @description calculate spatial enrichment using a matrix approach
#' @keywords internal
calc_spatial_enrichment_matrix = function(spatial_network,
                                          bin_matrix,
                                          adjust_method = 'fdr',
                                          do_fisher_test = TRUE,
                                          do_parallel = TRUE,
                                          cores = NA,
                                          calc_hub = FALSE,
                                          hub_min_int = 3,
                                          verbose = TRUE) {


  # data.table variables
  verbose = genes = p.value = estimate = adj.p.value = score = NULL

  # convert spatial network data.table to spatial matrix
  dc_spat_network = data.table::dcast.data.table(spatial_network, formula = to~from, value.var = 'distance', fill = 0)
  spat_mat = dt_to_matrix(dc_spat_network)
  spat_mat[spat_mat > 0] = 1


  ## parallel
  if(do_parallel == TRUE) {

    if(do_fisher_test == TRUE) {

      save_list = suppressMessages(giotto_lapply(X = rownames(bin_matrix), cores = cores, fun = spat_fish_func,
                                                          bin_matrix = bin_matrix, spat_mat = spat_mat,
                                                          calc_hub = calc_hub, hub_min_int = hub_min_int))

    } else {
      save_list =  suppressMessages(giotto_lapply(X = rownames(bin_matrix), cores = cores, fun = spat_OR_func,
                                                           bin_matrix = bin_matrix, spat_mat = spat_mat,
                                                           calc_hub = calc_hub, hub_min_int = hub_min_int))

    }

  } else {

    ## serial
    save_list = list()

    if(do_fisher_test == TRUE) {
      for(gene in rownames(bin_matrix)) {
        if(verbose == TRUE) print(gene)

        save_list[[gene]] = suppressMessages(spat_fish_func(gene = gene, bin_matrix = bin_matrix, spat_mat = spat_mat,
                                                                     calc_hub = calc_hub, hub_min_int = hub_min_int))

      }
    } else {
      for(gene in rownames(bin_matrix)) {
        if(verbose == TRUE) print(gene)

        save_list[[gene]] = suppressMessages(spat_OR_func(gene = gene, bin_matrix = bin_matrix, spat_mat = spat_mat,
                                                                   calc_hub = calc_hub, hub_min_int = hub_min_int))

      }
    }

  }

  result = data.table::as.data.table(do.call('rbind', save_list))
  result[, genes := unlist(genes)]


  if(do_fisher_test == TRUE) {
    result[, c('p.value', 'estimate') := list(as.numeric(p.value), as.numeric(estimate))]

    # convert p.value = 0 to lowest p-value
    min_pvalue = min(result$p.value[result$p.value > 0])
    result[, p.value := ifelse(p.value == 0, min_pvalue, p.value)]
    result[, adj.p.value := stats::p.adjust(p.value, method = adjust_method)]

    result[, score := -log(p.value) * estimate]
    data.table::setorder(result, -score)

  } else {

    data.table::setnames(result, 'V1', 'estimate')
    data.table::setorder(result, -estimate)
  }

  return(result)

}


#' @title calc_spatial_enrichment_DT
#' @name calc_spatial_enrichment_DT
#' @description calculate spatial enrichment using the data.table implementation
#' @keywords internal
calc_spatial_enrichment_DT = function(bin_matrix,
                                      spatial_network,
                                      calc_hub = F,
                                      hub_min_int = 3,
                                      group_size = 'automatic',
                                      do_fisher_test = TRUE,
                                      adjust_method = 'fdr',
                                      cores = NA) {


  # set number of cores automatically, but with limit of 10
  cores = determine_cores(cores)
  data.table::setDTthreads(threads = cores)

  # data.table variables
  from = to = gene_ID = p.value = adj.p.value = score = estimate = NULL

  # create minimum spatial network
  spat_netw_min = spatial_network[,.(from, to)]

  # divide matrix in groups
  if(!is.na(group_size) & is.numeric(group_size)) {
    group_size = group_size
    if(group_size > nrow(bin_matrix)) {
      stop('group_size is too big, it can not be greater than the number of genes')
    }
  } else if(group_size == 'automatic') {
    test_number = ceiling(nrow(bin_matrix)/10)
    test_number = max(2, test_number)
    group_size = min(200, test_number)
  }

  groups = ceiling(nrow(bin_matrix)/group_size)
  cut_groups = cut(1:nrow(bin_matrix), breaks = groups, labels = 1:groups)
  indexes = 1:nrow(bin_matrix)
  names(indexes) = cut_groups


  total_list = list()
  for(group in unique(cut_groups)) {

    sel_indices = indexes[names(indexes) == group]

    bin_matrix_DT = data.table::as.data.table(bin_matrix[sel_indices,])
    bin_matrix_DT[, gene_ID := rownames(bin_matrix[sel_indices,])]
    bin_matrix_DTm = data.table::melt.data.table(bin_matrix_DT, id.vars = 'gene_ID')

    if(do_fisher_test == TRUE) {
      test = spat_fish_func_DT(bin_matrix_DTm = bin_matrix_DTm,
                               spat_netw_min = spat_netw_min,
                               calc_hub = calc_hub,
                               hub_min_int = hub_min_int,
                               cores = cores)
    } else {
      test = spat_OR_func_DT(bin_matrix_DTm = bin_matrix_DTm,
                             spat_netw_min = spat_netw_min,
                             calc_hub = calc_hub,
                             hub_min_int = hub_min_int,
                             cores = cores)
    }


    total_list[[group]] = test

  }

  result = do.call('rbind', total_list)

  if(do_fisher_test == TRUE) {
    min_pvalue = min(result$p.value[result$p.value > 0])
    result[, p.value := ifelse(p.value == 0, min_pvalue, p.value)]
    result[, adj.p.value := stats::p.adjust(p.value, method = adjust_method)]

    result[, score := -log(p.value) * estimate]
    data.table::setorder(result, -score)
    data.table::setnames(result, 'gene_ID', 'genes')
  } else {
    data.table::setorder(result, -estimate)
    data.table::setnames(result, 'gene_ID', 'genes')
  }

  return(result)

}







#' @title binSpectSingle
#' @name binSpectSingle
#' @description binSpect for a single spatial network
#' @param gobject giotto object
#' @param bin_method method to binarize gene expression
#' @param expression_values expression values to use
#' @param subset_genes only select a subset of genes to test
#' @param spatial_network_name name of spatial network to use (default = 'spatial_network')
#' @param reduce_network default uses the full network
#' @param kmeans_algo kmeans algorithm to use (kmeans, kmeans_arma, kmeans_arma_subset)
#' @param nstart kmeans: nstart parameter
#' @param iter_max kmeans: iter.max parameter
#' @param extreme_nr number of top and bottom cells (see details)
#' @param sample_nr total number of cells to sample (see details)
#' @param percentage_rank percentage of top cells for binarization
#' @param do_fisher_test perform fisher test
#' @param adjust_method p-value adjusted method to use (see \code{\link[stats]{p.adjust}})
#' @param calc_hub calculate the number of hub cells
#' @param hub_min_int minimum number of cell-cell interactions for a hub cell
#' @param get_av_expr calculate the average expression per gene of the high expressing cells
#' @param get_high_expr calculate the number of high expressing cells  per gene
#' @param implementation enrichment implementation (data.table, simple, matrix)
#' @param group_size number of genes to process together with data.table implementation (default = automatic)
#' @param do_parallel run calculations in parallel with mclapply
#' @param cores number of cores to use if do_parallel = TRUE
#' @param verbose be verbose
#' @param set.seed set a seed before kmeans binarization
#' @param bin_matrix a binarized matrix, when provided it will skip the binarization process
#' @return data.table with results (see details)
#' @details We provide two ways to identify spatial genes based on gene expression binarization.
#' Both methods are identicial except for how binarization is performed.
#' \itemize{
#'   \item{1. binarize: }{Each gene is binarized (0 or 1) in each cell with \bold{kmeans} (k = 2) or based on \bold{rank} percentile}
#'   \item{2. network: }{Alll cells are connected through a spatial network based on the physical coordinates}
#'   \item{3. contingency table: }{A contingency table is calculated based on all edges of neighboring cells and the binarized expression (0-0, 0-1, 1-0 or 1-1)}
#'   \item{4. For each gene an odds-ratio (OR) and fisher.test (optional) is calculated}
#' }
#' Three different kmeans algorithmes have been implemented:
#' \itemize{
#'   \item{1. kmeans: }{default, see \code{\link[stats]{kmeans}} }
#'   \item{2. kmeans_arma: }{from ClusterR, see \code{\link[ClusterR]{KMeans_arma}} }
#'   \item{3. kmeans_arma_subst: }{from ClusterR, see \code{\link[ClusterR]{KMeans_arma}},
#'    but random subsetting the vector for each gene to increase speed. Change extreme_nr and sample_nr for control.  }
#' }
#' Other statistics are provided (optional):
#' \itemize{
#'   \item{Number of cells with high expression (binary = 1)}
#'   \item{Average expression of each gene within high expressing cells }
#'   \item{Number of hub cells, these are high expressing cells that have a user defined number of high expressing neighbors}
#' }
#' By selecting a subset of likely spatial genes (e.g. soft thresholding highly variable genes) can accelerate the speed.
#' The simple implementation is usually faster, but lacks the possibility to run in parallel and to calculate hub cells.
#' The data.table implementation might be more appropriate for large datasets by setting the group_size (number of genes) parameter to divide the workload.
#' @export
binSpectSingle = function(gobject,
                          bin_method = c('kmeans', 'rank'),
                          expression_values = c('normalized', 'scaled', 'custom'),
                          subset_genes = NULL,
                          spatial_network_name = 'Delaunay_network',
                          reduce_network = FALSE,
                          kmeans_algo = c('kmeans', 'kmeans_arma', 'kmeans_arma_subset'),
                          nstart = 3,
                          iter_max = 10,
                          extreme_nr = 50,
                          sample_nr = 50,
                          percentage_rank = 30,
                          do_fisher_test = TRUE,
                          adjust_method = 'fdr',
                          calc_hub = FALSE,
                          hub_min_int = 3,
                          get_av_expr = TRUE,
                          get_high_expr = TRUE,
                          implementation = c('data.table', 'simple', 'matrix'),
                          group_size = 'automatic',
                          do_parallel = TRUE,
                          cores = NA,
                          verbose = T,
                          set.seed = NULL,
                          bin_matrix = NULL) {

  if(verbose == TRUE) cat('\n This is the single parameter version of binSpect')


  # set number of cores automatically, but with limit of 10
  cores = determine_cores(cores)
  data.table::setDTthreads(threads = cores)

  # data.table: set global variable
  genes = p.value = estimate = score = NULL

  # set binarization method
  bin_method = match.arg(bin_method, choices = c('kmeans', 'rank'))

  # kmeans algorithm
  kmeans_algo = match.arg(kmeans_algo, choices = c('kmeans', 'kmeans_arma', 'kmeans_arma_subset'))

  # implementation
  implementation = match.arg(implementation, choices = c('data.table', 'simple', 'matrix'))


  # spatial network
  spatial_network = select_spatialNetwork(gobject,name = spatial_network_name,return_network_Obj = FALSE)
  if(is.null(spatial_network)) {
    stop('spatial_network_name: ', spatial_network_name, ' does not exist, create a spatial network first')
  }

  # convert to full network
  if(reduce_network == FALSE) {
    spatial_network = convert_to_full_spatial_network(spatial_network)
    setnames(spatial_network, c('source', 'target'), c('from', 'to'))
  }




  ## start binarization ##
  ## ------------------ ##

  if(!is.null(bin_matrix)) {
    bin_matrix = bin_matrix
  } else {
    if(bin_method == 'kmeans') {

      bin_matrix = kmeans_binarize_wrapper(gobject = gobject,
                                           expression_values = expression_values,
                                           subset_genes = subset_genes,
                                           kmeans_algo = kmeans_algo,
                                           nstart = nstart,
                                           iter_max = iter_max,
                                           extreme_nr = extreme_nr,
                                           sample_nr = sample_nr,
                                           set.seed = set.seed)

    } else if(bin_method == 'rank') {

      # expression
      values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))
      expr_values = select_expression_values(gobject = gobject, values = values)

      if(!is.null(subset_genes)) {
        expr_values = expr_values[rownames(expr_values) %in% subset_genes, ]
      }

      max_rank = (ncol(expr_values)/100)*percentage_rank
      bin_matrix = t_giotto(apply(X = expr_values, MARGIN = 1, FUN = rank_binarize, max_rank = max_rank))
    }
  }

  if(verbose == TRUE) cat('\n 1. matrix binarization complete \n')

  ## start with enrichment ##
  ## --------------------- ##

  if(implementation == 'simple') {
    if(do_parallel == TRUE) {
      warning('Parallel not yet implemented for simple. Enrichment will default to serial.')
    }

    if(calc_hub == TRUE) {
      warning('Hub calculation is not possible with the simple implementation, change to matrix if requird.')
    }


    result = calc_spatial_enrichment_minimum(spatial_network = spatial_network,
                                             bin_matrix = bin_matrix,
                                             adjust_method = adjust_method,
                                             do_fisher_test = do_fisher_test)


  } else if(implementation == 'matrix') {

    result = calc_spatial_enrichment_matrix(spatial_network = spatial_network,
                                            bin_matrix = bin_matrix,
                                            adjust_method = adjust_method,
                                            do_fisher_test = do_fisher_test,
                                            do_parallel = do_parallel,
                                            cores = cores,
                                            calc_hub = calc_hub,
                                            hub_min_int = hub_min_int)

  } else if(implementation == 'data.table') {

    result = calc_spatial_enrichment_DT(bin_matrix = bin_matrix,
                                        spatial_network = spatial_network,
                                        calc_hub = calc_hub,
                                        hub_min_int = hub_min_int,
                                        group_size = group_size,
                                        do_fisher_test = do_fisher_test,
                                        adjust_method = adjust_method,
                                        cores = cores)
  }

  if(verbose == TRUE) cat('\n 2. spatial enrichment test completed \n')





  ## start with average high expression ##
  ## ---------------------------------- ##

  if(get_av_expr == TRUE) {

    # expression
    values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))
    expr_values = select_expression_values(gobject = gobject, values = values)

    if(!is.null(subset_genes)) {
      expr_values = expr_values[rownames(expr_values) %in% subset_genes, ]
    }

    sel_expr_values = expr_values * bin_matrix
    av_expr = apply(sel_expr_values, MARGIN = 1, FUN = function(x) {
      mean(x[x > 0])
    })
    av_expr_DT = data.table::data.table(genes = names(av_expr), av_expr = av_expr)
    result = merge(result, av_expr_DT, by = 'genes')

    if(verbose == TRUE) cat('\n 3. (optional) average expression of high expressing cells calculated \n')
  }



  ## start with number of high expressing cells ##
  ## ------------------------------------------ ##

  if(get_high_expr == TRUE) {
    high_expr = rowSums(bin_matrix)
    high_expr_DT = data.table::data.table(genes = names(high_expr), high_expr = high_expr)
    result = merge(result, high_expr_DT, by = 'genes')

    if(verbose == TRUE) cat('\n 4. (optional) number of high expressing cells calculated \n')
  }


  # sort
  if(do_fisher_test == TRUE) {
    data.table::setorder(result, -score)
  } else {
    data.table::setorder(result, -estimate)
  }

  return(result)

}





#' @title binSpectMulti
#' @name binSpectMulti
#' @description binSpect for multiple spatial kNN networks
#' @param gobject giotto object
#' @param bin_method method to binarize gene expression
#' @param expression_values expression values to use
#' @param subset_genes only select a subset of genes to test
#' @param spatial_network_k different k's for a spatial kNN to evaluate
#' @param reduce_network default uses the full network
#' @param kmeans_algo kmeans algorithm to use (kmeans, kmeans_arma, kmeans_arma_subset)
#' @param nstart kmeans: nstart parameter
#' @param iter_max kmeans: iter.max parameter
#' @param extreme_nr number of top and bottom cells (see details)
#' @param sample_nr total number of cells to sample (see details)
#' @param percentage_rank percentage of top cells for binarization
#' @param do_fisher_test perform fisher test
#' @param adjust_method p-value adjusted method to use (see \code{\link[stats]{p.adjust}})
#' @param calc_hub calculate the number of hub cells
#' @param hub_min_int minimum number of cell-cell interactions for a hub cell
#' @param get_av_expr calculate the average expression per gene of the high expressing cells
#' @param get_high_expr calculate the number of high expressing cells  per gene
#' @param implementation enrichment implementation (data.table, simple, matrix)
#' @param group_size number of genes to process together with data.table implementation (default = automatic)
#' @param do_parallel run calculations in parallel with mclapply
#' @param cores number of cores to use if do_parallel = TRUE
#' @param verbose be verbose
#' @param knn_params list of parameters to create spatial kNN network
#' @param set.seed set a seed before kmeans binarization
#' @param summarize summarize the p-values or adjusted p-values
#' @return data.table with results (see details)
#' @details We provide two ways to identify spatial genes based on gene expression binarization.
#' Both methods are identicial except for how binarization is performed.
#' \itemize{
#'   \item{1. binarize: }{Each gene is binarized (0 or 1) in each cell with \bold{kmeans} (k = 2) or based on \bold{rank} percentile}
#'   \item{2. network: }{Alll cells are connected through a spatial network based on the physical coordinates}
#'   \item{3. contingency table: }{A contingency table is calculated based on all edges of neighboring cells and the binarized expression (0-0, 0-1, 1-0 or 1-1)}
#'   \item{4. For each gene an odds-ratio (OR) and fisher.test (optional) is calculated}
#' }
#' Three different kmeans algorithmes have been implemented:
#' \itemize{
#'   \item{1. kmeans: }{default, see \code{\link[stats]{kmeans}} }
#'   \item{2. kmeans_arma: }{from ClusterR, see \code{\link[ClusterR]{KMeans_arma}} }
#'   \item{3. kmeans_arma_subst: }{from ClusterR, see \code{\link[ClusterR]{KMeans_arma}},
#'    but random subsetting the vector for each gene to increase speed. Change extreme_nr and sample_nr for control.  }
#' }
#' Other statistics are provided (optional):
#' \itemize{
#'   \item{Number of cells with high expression (binary = 1)}
#'   \item{Average expression of each gene within high expressing cells }
#'   \item{Number of hub cells, these are high expressing cells that have a user defined number of high expressing neighbors}
#' }
#' By selecting a subset of likely spatial genes (e.g. soft thresholding highly variable genes) can accelerate the speed.
#' The simple implementation is usually faster, but lacks the possibility to run in parallel and to calculate hub cells.
#' The data.table implementation might be more appropriate for large datasets by setting the group_size (number of genes) parameter to divide the workload.
#' @export
binSpectMulti = function(gobject,
                         bin_method = c('kmeans', 'rank'),
                         expression_values = c('normalized', 'scaled', 'custom'),
                         subset_genes = NULL,
                         spatial_network_k = c(5, 10, 20),
                         reduce_network = FALSE,
                         kmeans_algo = c('kmeans', 'kmeans_arma', 'kmeans_arma_subset'),
                         nstart = 3,
                         iter_max = 10,
                         extreme_nr = 50,
                         sample_nr = 50,
                         percentage_rank = c(10, 30),
                         do_fisher_test = TRUE,
                         adjust_method = 'fdr',
                         calc_hub = FALSE,
                         hub_min_int = 3,
                         get_av_expr = TRUE,
                         get_high_expr = TRUE,
                         implementation = c('data.table', 'simple', 'matrix'),
                         group_size = 'automatic',
                         do_parallel = TRUE,
                         cores = NA,
                         verbose = T,
                         knn_params = NULL,
                         set.seed = NULL,
                         summarize = c('adj.p.value', 'p.value')
                         ) {


  if(verbose == TRUE) cat('\n This is the multi parameter version of binSpect')

  # set number of cores automatically, but with limit of 10
  cores = determine_cores(cores)
  data.table::setDTthreads(threads = cores)

  # check bin_method
  bin_method = match.arg(bin_method, choices = c('kmeans', 'rank'))

  # summarization level
  summarize = match.arg(summarize, choices = c('adj.p.value', 'p.value'))

  ## bin method rank
  if(bin_method == 'rank') {

    total_trials = length(spatial_network_k)*length(percentage_rank)
    result_list = vector(mode = 'list', length = total_trials)
    i = 1

    for(k in spatial_network_k) {

      if(is.null(knn_params)) {
        knn_params = list(minimum_k = 1)
      }
      temp_gobject = do.call('createSpatialKNNnetwork', c(gobject = gobject,
                                                          name = 'temp_knn_network',
                                                          k = k, knn_params))

      for(rank_i in percentage_rank) {

        if(verbose == TRUE) cat('\n Run for k = ', k, ' and rank % = ', rank_i,'\n')

        result = binSpectSingle(gobject = temp_gobject,
                                bin_method = bin_method,
                                expression_values = expression_values,
                                subset_genes = subset_genes,
                                spatial_network_name = 'temp_knn_network',
                                reduce_network = reduce_network,
                                kmeans_algo = kmeans_algo,
                                percentage_rank = rank_i,
                                do_fisher_test = do_fisher_test,
                                adjust_method = adjust_method,
                                calc_hub = calc_hub,
                                hub_min_int = hub_min_int,
                                get_av_expr = get_av_expr,
                                get_high_expr = get_high_expr,
                                implementation = implementation,
                                group_size = group_size,
                                do_parallel = do_parallel,
                                cores = cores,
                                verbose = verbose,
                                set.seed = set.seed)

        result_list[[i]] = result
        i = i+1
      }
    }
    combined_result = data.table::rbindlist(result_list)

  } else if(bin_method == 'kmeans') {

    ## bin method kmeans
    total_trials = length(spatial_network_k)
    result_list = vector(mode = 'list', length = total_trials)
    i = 1


    # pre-calculate bin_matrix once
    bin_matrix = kmeans_binarize_wrapper(gobject = gobject,
                                         expression_values = expression_values,
                                         subset_genes = subset_genes,
                                         kmeans_algo = kmeans_algo,
                                         nstart = nstart,
                                         iter_max = iter_max,
                                         extreme_nr = extreme_nr,
                                         sample_nr = sample_nr,
                                         set.seed = set.seed)

    for(k in spatial_network_k) {

      if(is.null(knn_params)) {
        knn_params = list(minimum_k = 1)
      }
      temp_gobject = do.call('createSpatialKNNnetwork', c(gobject = gobject,
                                                          name = 'temp_knn_network',
                                                          k = k, knn_params))

      if(verbose == TRUE) cat('\n Run for k = ', k,'\n')

      result = binSpectSingle(gobject = temp_gobject,
                              bin_method = bin_method,
                              expression_values = expression_values,
                              subset_genes = subset_genes,
                              spatial_network_name = 'temp_knn_network',
                              reduce_network = reduce_network,
                              kmeans_algo = kmeans_algo,
                              nstart = nstart,
                              iter_max = iter_max,
                              extreme_nr = extreme_nr,
                              sample_nr = sample_nr,
                              do_fisher_test = do_fisher_test,
                              adjust_method = adjust_method,
                              calc_hub = calc_hub,
                              hub_min_int = hub_min_int,
                              get_av_expr = get_av_expr,
                              get_high_expr = get_high_expr,
                              implementation = implementation,
                              group_size = group_size,
                              do_parallel = do_parallel,
                              cores = cores,
                              verbose = verbose,
                              set.seed = set.seed,
                              bin_matrix = bin_matrix)

      result_list[[i]] = result
      i = i+1

    }

    combined_result = data.table::rbindlist(result_list)

  }


  # data.table variables
  genes = V1 = p.val = NULL

  ## merge results into 1 p-value per gene ##
  simple_result = combined_result[, sum(log(get(summarize))), by = genes]
  simple_result[, V1 := V1*-2]
  simple_result[, p.val := stats::pchisq(q = V1, df = total_trials, log.p = F, lower.tail = F)]

  return(list(combined = combined_result, simple = simple_result[,.(genes, p.val)]))

}



#' @title binSpect
#' @name binSpect
#' @description Previously: binGetSpatialGenes. BinSpect (Binary Spatial Extraction of genes) is a fast computational method
#' that identifies genes with a spatially coherent expression pattern.
#' @param gobject giotto object
#' @param bin_method method to binarize gene expression
#' @param expression_values expression values to use
#' @param subset_genes only select a subset of genes to test
#' @param spatial_network_name name of spatial network to use (default = 'spatial_network')
#' @param spatial_network_k different k's for a spatial kNN to evaluate
#' @param reduce_network default uses the full network
#' @param kmeans_algo kmeans algorithm to use (kmeans, kmeans_arma, kmeans_arma_subset)
#' @param nstart kmeans: nstart parameter
#' @param iter_max kmeans: iter.max parameter
#' @param extreme_nr number of top and bottom cells (see details)
#' @param sample_nr total number of cells to sample (see details)
#' @param percentage_rank percentage of top cells for binarization
#' @param do_fisher_test perform fisher test
#' @param adjust_method p-value adjusted method to use (see \code{\link[stats]{p.adjust}})
#' @param calc_hub calculate the number of hub cells
#' @param hub_min_int minimum number of cell-cell interactions for a hub cell
#' @param get_av_expr calculate the average expression per gene of the high expressing cells
#' @param get_high_expr calculate the number of high expressing cells  per gene
#' @param implementation enrichment implementation (data.table, simple, matrix)
#' @param group_size number of genes to process together with data.table implementation (default = automatic)
#' @param do_parallel run calculations in parallel with mclapply
#' @param cores number of cores to use if do_parallel = TRUE
#' @param verbose be verbose
#' @param knn_params list of parameters to create spatial kNN network
#' @param set.seed set a seed before kmeans binarization
#' @param bin_matrix a binarized matrix, when provided it will skip the binarization process
#' @param summarize summarize the p-values or adjusted p-values
#' @return data.table with results (see details)
#' @details We provide two ways to identify spatial genes based on gene expression binarization.
#' Both methods are identicial except for how binarization is performed.
#' \itemize{
#'   \item{1. binarize: }{Each gene is binarized (0 or 1) in each cell with \bold{kmeans} (k = 2) or based on \bold{rank} percentile}
#'   \item{2. network: }{Alll cells are connected through a spatial network based on the physical coordinates}
#'   \item{3. contingency table: }{A contingency table is calculated based on all edges of neighboring cells and the binarized expression (0-0, 0-1, 1-0 or 1-1)}
#'   \item{4. For each gene an odds-ratio (OR) and fisher.test (optional) is calculated}
#' }
#' Three different kmeans algorithmes have been implemented:
#' \itemize{
#'   \item{1. kmeans: }{default, see \code{\link[stats]{kmeans}} }
#'   \item{2. kmeans_arma: }{from ClusterR, see \code{\link[ClusterR]{KMeans_arma}} }
#'   \item{3. kmeans_arma_subst: }{from ClusterR, see \code{\link[ClusterR]{KMeans_arma}},
#'    but random subsetting the vector for each gene to increase speed. Change extreme_nr and sample_nr for control.  }
#' }
#' Other statistics are provided (optional):
#' \itemize{
#'   \item{Number of cells with high expression (binary = 1)}
#'   \item{Average expression of each gene within high expressing cells }
#'   \item{Number of hub cells, these are high expressing cells that have a user defined number of high expressing neighbors}
#' }
#' By selecting a subset of likely spatial genes (e.g. soft thresholding highly variable genes) can accelerate the speed.
#' The simple implementation is usually faster, but lacks the possibility to run in parallel and to calculate hub cells.
#' The data.table implementation might be more appropriate for large datasets by setting the group_size (number of genes) parameter to divide the workload.
#' @export
binSpect = function(gobject,
                    bin_method = c('kmeans', 'rank'),
                    expression_values = c('normalized', 'scaled', 'custom'),
                    subset_genes = NULL,
                    spatial_network_name = 'Delaunay_network',
                    spatial_network_k = NULL,
                    reduce_network = FALSE,
                    kmeans_algo = c('kmeans', 'kmeans_arma', 'kmeans_arma_subset'),
                    nstart = 3,
                    iter_max = 10,
                    extreme_nr = 50,
                    sample_nr = 50,
                    percentage_rank = 30,
                    do_fisher_test = TRUE,
                    adjust_method = 'fdr',
                    calc_hub = FALSE,
                    hub_min_int = 3,
                    get_av_expr = TRUE,
                    get_high_expr = TRUE,
                    implementation = c('data.table', 'simple', 'matrix'),
                    group_size = 'automatic',
                    do_parallel = TRUE,
                    cores = NA,
                    verbose = T,
                    knn_params = NULL,
                    set.seed = NULL,
                    bin_matrix = NULL,
                    summarize = c('p.value', 'adj.p.value')) {


  if(!is.null(spatial_network_k)) {

    output = binSpectMulti(gobject = gobject,
                           bin_method = bin_method,
                           expression_values = expression_values,
                           subset_genes = subset_genes,
                           spatial_network_k = spatial_network_k,
                           reduce_network = reduce_network,
                           kmeans_algo = kmeans_algo,
                           nstart = nstart,
                           iter_max = iter_max,
                           extreme_nr = extreme_nr,
                           sample_nr = sample_nr,
                           percentage_rank = percentage_rank,
                           do_fisher_test = do_fisher_test,
                           adjust_method = adjust_method,
                           calc_hub = calc_hub,
                           hub_min_int = hub_min_int,
                           get_av_expr = get_av_expr,
                           get_high_expr = get_high_expr,
                           implementation = implementation,
                           group_size = group_size,
                           do_parallel = do_parallel,
                           cores = cores,
                           verbose = verbose,
                           knn_params = knn_params,
                           set.seed = set.seed,
                           summarize = summarize)

  } else {

    output = binSpectSingle(gobject = gobject,
                            bin_method = bin_method,
                            expression_values = expression_values,
                            subset_genes = subset_genes,
                            spatial_network_name = spatial_network_name,
                            reduce_network = reduce_network,
                            kmeans_algo = kmeans_algo,
                            nstart = nstart,
                            iter_max = iter_max,
                            extreme_nr = extreme_nr,
                            sample_nr = sample_nr,
                            percentage_rank = percentage_rank,
                            do_fisher_test = do_fisher_test,
                            adjust_method = adjust_method,
                            calc_hub = calc_hub,
                            hub_min_int = hub_min_int,
                            get_av_expr = get_av_expr,
                            get_high_expr = get_high_expr,
                            implementation = implementation,
                            group_size = group_size,
                            do_parallel = do_parallel,
                            cores = cores,
                            verbose = verbose,
                            set.seed = set.seed,
                            bin_matrix = bin_matrix)

  }

  return(output)

}




#' @title silhouetteRank
#' @name silhouetteRank
#' @description Previously: calculate_spatial_genes_python. This method computes a silhouette score per gene based on the
#' spatial distribution of two partitions of cells (expressed L1, and non-expressed L0).
#' Here, rather than L2 Euclidean norm, it uses a rank-transformed, exponentially weighted
#' function to represent the local physical distance between two cells.
#' New multi aggregator implementation can be found at \code{\link{silhouetteRankTest}}
#' @param gobject giotto object
#' @param expression_values expression values to use
#' @param metric distance metric to use
#' @param subset_genes only run on this subset of genes
#' @param rbp_p fractional binarization threshold
#' @param examine_top top fraction to evaluate with silhouette
#' @param python_path specify specific path to python if required
#' @return data.table with spatial scores
#' @export
silhouetteRank <- function(gobject,
                           expression_values = c('normalized', 'scaled', 'custom'),
                           metric = "euclidean",
                           subset_genes = NULL,
                           rbp_p = 0.95,
                           examine_top = 0.3,
                           python_path = NULL) {


  # expression values
  values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))
  expr_values = select_expression_values(gobject = gobject, values = values)

  # subset genes
  if(!is.null(subset_genes)) {

    subset_genes = subset_genes[subset_genes %in% gobject@gene_ID]
    expr_values = expr_values[rownames(expr_values) %in% subset_genes, ]

  }


  # data.table variables
  sdimx = sdimy = NULL

  # spatial locations
  spatlocs = as.matrix(gobject@spatial_locs[,.(sdimx, sdimy)])

  # python path
  if(is.null(python_path)) {
    python_path = readGiottoInstructions(gobject, param = "python_path")
  }

  ## prepare python path and louvain script
  reticulate::use_python(required = T, python = python_path)
  python_silh_function = system.file("python", "python_spatial_genes.py", package = 'Giotto')
  reticulate::source_python(file = python_silh_function)


  output_python = python_spatial_genes(spatial_locations = spatlocs,
                                       expression_matrix = as.data.frame(expr_values),
                                       metric = metric,
                                       rbp_p = rbp_p,
                                       examine_top = examine_top)

  # unlist output
  genes = unlist(lapply(output_python, FUN = function(x) {
    y = x[1][[1]]
  }))
  scores = unlist(lapply(output_python, FUN = function(x) {
    y = x[2][[1]]
  }))

  spatial_python_DT = data.table::data.table(genes = genes, scores = scores)

  return(spatial_python_DT)


}




#' @title silhouetteRankTest
#' @name silhouetteRankTest
#' @description Multi parameter aggregator version of \code{\link{silhouetteRank}}
#' @param gobject giotto object
#' @param expression_values expression values to use
#' @param subset_genes only run on this subset of genes
#' @param overwrite_input_bin overwrite input bin
#' @param rbp_ps fractional binarization thresholds
#' @param examine_tops top fractions to evaluate with silhouette
#' @param matrix_type type of matrix
#' @param num_core number of cores to use
#' @param parallel_path path to GNU parallel function
#' @param output output directory
#' @param query_sizes size of query
#' @param verbose be verbose
#' @return data.table with spatial scores
#' @export
silhouetteRankTest = function(gobject,
                               expression_values = c('normalized', 'scaled', 'custom'),
                               subset_genes = NULL,
                               overwrite_input_bin = TRUE,
                               rbp_ps = c(0.95, 0.99),
                               examine_tops = c(0.005, 0.010, 0.050, 0.100, 0.300),
                               matrix_type = "dissim",
                               num_core = 4,
                               parallel_path = "/usr/bin",
                               output = NULL,
                               query_sizes = 10L,
                               verbose = FALSE) {


  # data.table variables
  cell_ID = sdimx = sdimy = sdimz = NULL

  ## test if R packages are installed
  # check envstats
  package_check(pkg_name = 'EnvStats', repository = c('CRAN'))

  # check eva
  if ('eva' %in% rownames(installed.packages()) == FALSE) {
    stop("\n package ", 'eva', " is not yet installed \n",
         "To install: \n",
         "install.packages('eva_0.2.5.tar.gz', repos=NULL, type='source')",
         "see https://cran.r-project.org/src/contrib/Archive/eva/")
  }

  ## test if python package is installed
  module_test = reticulate::py_module_available('silhouetteRank')
  if(module_test == FALSE) {
    warning("silhouetteRank python module is not installed:
            install in the right environment or python path with:

            'pip install silhouetteRank'

            or from within R in the Giotto environment with:

            conda_path = reticulate::miniconda_path()
            conda_full_path = paste0(conda_path,'/','bin/conda')
            full_envname = paste0(conda_path,'/envs/giotto_env')
            reticulate::py_install(packages = 'silhouetteRank',
                       envname = full_envname,
                       method = 'conda',
                       conda = conda_full_path,
                       pip = TRUE,
                       python_version = '3.6')")
  }



  # expression values
  values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))
  expr_values = select_expression_values(gobject = gobject, values = values)

  # subset genes
  if(!is.null(subset_genes)) {

    subset_genes = subset_genes[subset_genes %in% gobject@gene_ID]
    expr_values = expr_values[rownames(expr_values) %in% subset_genes, ]

  }

  # spatial locations
  spatlocs = gobject@spatial_locs

  ## save dir and log
  if(is.null(output)) {

    save_dir = readGiottoInstructions(gobject, param = "save_dir")
    silh_output_dir = paste0(save_dir, '/', 'silhouetteRank_output/')
    if(!file.exists(silh_output_dir)) dir.create(silh_output_dir, recursive = TRUE)

  } else if(file.exists(output)) {

    silh_output_dir = paste0(output, '/', 'silhouetteRank_output/')
    if(!file.exists(silh_output_dir)) dir.create(silh_output_dir, recursive = TRUE)

  } else {

    silh_output_dir = paste0(output, '/', 'silhouetteRank_output/')
    if(!file.exists(silh_output_dir)) dir.create(silh_output_dir, recursive = TRUE)

  }

  # log directory
  log_dir  = paste0(silh_output_dir, '/', 'logs/')
  if(!file.exists(log_dir)) dir.create(log_dir, recursive = TRUE)


  ## write spatial locations to .txt file
  if(ncol(spatlocs) == 3) {
    format_spatlocs = spatlocs[,.(cell_ID, sdimx, sdimy)]
    colnames(format_spatlocs) = c('ID', 'x', 'y')
  } else {
    format_spatlocs = spatlocs[,.(cell_ID, sdimx, sdimy, sdimz)]
    colnames(format_spatlocs) = c('ID', 'x', 'y', 'z')
  }

  write.table(x = format_spatlocs, row.names = F,
              file = paste0(silh_output_dir,'/', 'format_spatlocs.txt'),
              quote = F, sep = '\t')

  spatlocs_path = paste0(silh_output_dir,'/', 'format_spatlocs.txt')

  ## write expression to .txt file
  write.table(x = as.matrix(expr_values),
              file = paste0(silh_output_dir,'/', 'expression.txt'),
              quote = F, sep = '\t', col.names=NA)

  expr_values_path = paste0(silh_output_dir,'/', 'expression.txt')



  ## prepare python path and louvain script
  python_path = readGiottoInstructions(gobject, param = 'python_path')
  reticulate::use_python(required = T, python = python_path)
  python_silh_function = system.file("python", "silhouette_rank_wrapper.py", package = 'Giotto')
  reticulate::source_python(file = python_silh_function)


  output_silh = silhouette_rank(expr = expr_values_path,
                                centroid = spatlocs_path,
                                overwrite_input_bin = overwrite_input_bin,
                                rbp_ps = rbp_ps,
                                examine_tops = examine_tops,
                                matrix_type = matrix_type,
                                verbose = verbose,
                                num_core = num_core,
                                parallel_path = parallel_path,
                                output = silh_output_dir,
                                query_sizes = as.integer(query_sizes))

  return(output_silh)

}






#' @title spatialDE
#' @name spatialDE
#' @description Compute spatial variable genes with spatialDE method
#' @param gobject Giotto object
#' @param expression_values gene expression values to use
#' @param size size of plot
#' @param color low/medium/high color scheme for plot
#' @param sig_alpha alpha value for significance
#' @param unsig_alpha alpha value for unsignificance
#' @param python_path specify specific path to python if required
#' @param show_plot show plot
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters, see \code{\link{showSaveParameters}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return a list of data.frames with results and plot (optional)
#' @details This function is a wrapper for the SpatialDE method implemented in the ...
#' @export
spatialDE <- function(gobject = NULL,
                      expression_values = c('raw', 'normalized', 'scaled', 'custom'),
                      size = c(4,2,1),
                      color = c("blue", "green", "red"),
                      sig_alpha = 0.5,
                      unsig_alpha = 0.5,
                      python_path = NULL,
                      show_plot = NA,
                      return_plot = NA,
                      save_plot = NA,
                      save_param = list(),
                      default_save_name = 'SpatialDE'){



  # test if SPARK is installed ##

  module_test = reticulate::py_module_available('SpatialDE')
  if(module_test == FALSE) {
    warning("SpatialDE python module is not installed:
            install in the right environment or python path with:

            'pip install spatialde'

            or from within R in the Giotto environment with:

            conda_path = reticulate::miniconda_path()
            conda_full_path = paste0(conda_path,'/','bin/conda')
            full_envname = paste0(conda_path,'/envs/giotto_env')
            reticulate::py_install(packages = c('NaiveDE', 'patsy', 'SpatialDE'),
                                   envname = full_envname,
                                   method = 'conda',
                                   conda = conda_full_path,
                                   pip = TRUE,
                                   python_version = '3.6')")
  }


  # print message with information #
  message("using 'SpatialDE' for spatial gene/pattern detection. If used in published research, please cite:
  Svensson, Valentine, Sarah A. Teichmann, and Oliver Stegle. 'SpatialDE: Identification of Spatially Variable Genes.'
          Nature Methods 15, no. 5 (May 2018): 343-46. https://doi.org/10.1038/nmeth.4636.")



  # data.table variables
  cell_ID = NULL

  # expression
  values = match.arg(expression_values, c('raw', 'normalized', 'scaled', 'custom'))
  expr_values = select_expression_values(gobject = gobject, values = values)

  ## python path
  if(is.null(python_path)) {
    python_path = readGiottoInstructions(gobject, param = "python_path")
  }

  ## source python file
  reticulate::use_python(required = T, python = python_path)
  reader_path = system.file("python", "SpatialDE_wrapper.py", package = 'Giotto')
  reticulate::source_python(file = reader_path)

  ## get spatial locations
  spatial_locs <- as.data.frame(gobject@spatial_locs)
  rownames(spatial_locs) <- spatial_locs$cell_ID
  spatial_locs <- subset(spatial_locs, select = -cell_ID)

  ## run spatialDE
  Spatial_DE_results = Spatial_DE(as.data.frame(t(as.matrix(expr_values))), spatial_locs)

  results <- as.data.frame(reticulate::py_to_r(Spatial_DE_results[[1]]))

  if(length(Spatial_DE_results) == 2){
    ms_results = as.data.frame(reticulate::py_to_r(Spatial_DE_results[[2]]))
    spatial_genes_results = list(results, ms_results)
    names(spatial_genes_results) = c("results", "ms_results")
  } else{
    spatial_genes_results =  results
    ms_results = NULL
  }


  # print, return and save parameters
  show_plot = ifelse(is.na(show_plot), readGiottoInstructions(gobject, param = 'show_plot'), show_plot)
  save_plot = ifelse(is.na(save_plot), readGiottoInstructions(gobject, param = 'save_plot'), save_plot)
  return_plot = ifelse(is.na(return_plot), readGiottoInstructions(gobject, param = 'return_plot'), return_plot)

  ## create plot
  if(show_plot == TRUE | save_plot == TRUE | return_plot == TRUE) {
    FSV_plot = FSV_show(results = results,
                        ms_results = ms_results,
                        size =size,
                        color = color,
                        sig_alpha = sig_alpha,
                        unsig_alpha = unsig_alpha)
  }

  ## print plot
  if(show_plot == TRUE) {
    print(FSV_plot)
  }

  ## save plot
  if(save_plot == TRUE) {
    do.call('all_plots_save_function', c(list(gobject = gobject, plot_object = FSV_plot, default_save_name = default_save_name), save_param))
  }

  ## return results and plot (optional)
  if(return_plot == TRUE) {
    return(list(results = spatial_genes_results, plot = FSV_plot))
  } else {
    return(list(results =  spatial_genes_results))
  }

}


#' @title spatialAEH
#' @name spatialAEH
#' @description Compute spatial variable genes with spatialDE method
#' @param gobject Giotto object
#' @param SpatialDE_results results of \code{\link{spatialDE}} function
#' @param name_pattern name for the computed spatial patterns
#' @param expression_values gene expression values to use
#' @param pattern_num number of spatial patterns to look for
#' @param l lengthscale
#' @param python_path specify specific path to python if required
#' @param return_gobject show plot
#' @return An updated giotto object
#' @details This function is a wrapper for the SpatialAEH method implemented in the ...
#' @export
spatialAEH <- function(gobject = NULL,
                       SpatialDE_results = NULL,
                       name_pattern = 'AEH_patterns',
                       expression_values = c('raw', 'normalized', 'scaled', 'custom'),
                       pattern_num = 6,
                       l = 1.05,
                       python_path = NULL,
                       return_gobject = TRUE) {

  # data.table variables
  cell_ID = NULL

  # expression
  values = match.arg(expression_values, c('raw', 'normalized', 'scaled', 'custom'))
  expr_values = select_expression_values(gobject = gobject, values = values)

  ## python path
  if(is.null(python_path)) {
    python_path = readGiottoInstructions(gobject, param = "python_path")
  }

  ## source python file
  reticulate::use_python(required = T, python = python_path)
  reader_path = system.file("python", "SpatialDE_wrapper.py", package = 'Giotto')
  reticulate::source_python(file = reader_path)


  ## spatial locations
  spatial_locs <- as.data.frame(gobject@spatial_locs)
  rownames(spatial_locs) <- spatial_locs$cell_ID
  spatial_locs <- subset(spatial_locs, select = -cell_ID)

  # extract results you need
  results = SpatialDE_results[['results']][['results']]

  ## automatic expression histology
  AEH_results = Spatial_DE_AEH(filterd_exprs = as.data.frame(t_giotto(as.matrix(expr_values))),
                               coordinates = spatial_locs,
                               results = as.data.frame(results),
                               pattern_num = pattern_num,
                               l = l)
  histology_results <- as.data.frame(reticulate::py_to_r(AEH_results[[1]]))
  cell_pattern_score <- as.data.frame((reticulate::py_to_r(AEH_results[[2]])))

  spatial_pattern_results <- list(histology_results, cell_pattern_score)
  names(spatial_pattern_results) <- c("histology_results","cell_pattern_score")


  if(return_gobject == TRUE) {

    dt_res = as.data.table(spatial_pattern_results[['cell_pattern_score']])
    dt_res[['cell_ID']] = rownames(spatial_pattern_results[['cell_pattern_score']])
    gobject@spatial_enrichment[[name_pattern]] = dt_res
    return(gobject)

  } else {

    return(list(results = spatial_pattern_results))

  }
}


#' @title FSV_show
#' @name FSV_show
#' @description Visualize spatial varible genes caculated by spatial_DE
#' @param results results caculated by spatial_DE
#' @param ms_results ms_results caculated by spatial_DE
#' @param size indicate different levels of qval
#' @param color indicate different SV features
#' @param sig_alpha transparency of significant genes
#' @param unsig_alpha transparency of unsignificant genes
#' @return ggplot object
#' @details Description of parameters.
#' @keywords internal
FSV_show <- function(results,
                     ms_results = NULL,
                     size = c(4,2,1),
                     color = c("blue", "green", "red"),
                     sig_alpha = 0.5,
                     unsig_alpha = 0.5){

  results$FSV95conf = 2 * sqrt(results$s2_FSV)
  results$intervals <- cut(results$FSV95conf,c(0, 1e-1, 1e0, Inf),label = F)
  results$log_pval <- log10(results$pval)

  if(is.null(ms_results)){
    results$model_bic = results$model
  }
  else{
    results= merge(results,ms_results[,c("g","model")],by.x = "g",by.y = "g",all.x = T,
                   suffixes=(c(" ",'_bic')))
  }

  results$model_bic <- factor(results$model_bic)
  results$intervals <- factor(results$intervals)


  pl <- ggplot2::ggplot()
  pl <- pl + ggplot2::theme_bw()
  pl <- pl + ggplot2::geom_point(data = results[results$qval < 0.05,],
                                 ggplot2::aes_string(x = "FSV", y = "log_pval",fill = "model_bic",size = "intervals"),
                                 show.legend = T, shape = 21,alpha = sig_alpha,
                                 #size = size[results_cp_s$inftervals],
                                 stroke = 0.1, color = "black") +
    ggplot2::geom_point(data = results[results$qval > 0.05,],
                        ggplot2::aes_string(x = "FSV", y = "log_pval",size = "intervals"),
                        show.legend = T, shape = 21,alpha = unsig_alpha,
                        fill = "black", #size = size[results_cp_ns$inftervals],
                        stroke = 0.1, color = "black") +
    ggplot2::scale_size_manual(values = size,guide=FALSE)+
    ggplot2::scale_color_manual(values = color)+
    ggplot2::scale_fill_discrete(name="Spatial Patterns",
                                 breaks=c("linear", "PER", "SE"),
                                 labels=c("linear", "periodical", "general"))+
    ggplot2::geom_hline(yintercept = max(results[results$qval < 0.05,]$log_pval),linetype = "dashed")+
    ggplot2::geom_text(ggplot2::aes(0.9,max(results[results$qval < 0.05,]$log_pval),
                                    label = "FDR = 0.05", vjust = -1))+
    ggplot2::scale_y_reverse()

  print(pl)
}




#' @title trendSceek
#' @name trendSceek
#' @description Compute spatial variable genes with trendsceek method
#' @param gobject Giotto object
#' @param expression_values gene expression values to use
#' @param subset_genes subset of genes to run trendsceek on
#' @param nrand An integer specifying the number of random resamplings of the mark distribution as to create the null-distribution.
#' @param ncores An integer specifying the number of cores to be used by BiocParallel
#' @param \dots Additional parameters to the \code{\link[trendsceek]{trendsceek_test}} function
#' @return data.frame with trendsceek spatial genes results
#' @details This function is a wrapper for the trendsceek_test method implemented in the trendsceek package
#' @export
trendSceek <- function(gobject,
                       expression_values = c("normalized", "raw"),
                       subset_genes = NULL,
                       nrand = 100,
                       ncores = 8,
                       ...) {

  # verify if optional package is installed
  package_check(pkg_name = 'trendsceek',
                repository = c('github'),
                github_repo = 'edsgard/trendsceek')

  # print message with information #
  message("using 'trendsceek' for spatial gene/pattern detection. If used in published research, please cite:
  Edsgard, Daniel, Per Johnsson, and Rickard Sandberg. 'Identification of Spatial Expression Trends in Single-Cell Gene Expression Data.'
          Nature Methods 15, no. 5 (May 2018): 339-42. https://doi.org/10.1038/nmeth.4634.")

  ## expression data
  values = match.arg(expression_values, c("normalized", "raw"))
  expr_values = select_expression_values(gobject = gobject, values = values)

  ## normalization function
  if (values == "normalized") {
    log.fcn = NA
  }
  else if (values == "raw") {
    log.fcn = log10
  }

  ## subset genes
  if (!is.null(subset_genes)) {
    subset_genes = subset_genes[subset_genes %in% gobject@gene_ID]
    expr_values = expr_values[rownames(expr_values) %in% subset_genes, ]
  }


  ## initial locations

  # data.table variables
  cell_ID = NULL

  spatial_locations = copy(gobject@spatial_locs)
  spatial_locations[, cell_ID := NULL]
  pp = trendsceek::pos2pp(spatial_locations)

  ## initial gene counts
  pp = trendsceek::set_marks(pp, as.matrix(expr_values), log.fcn = log.fcn)

  # eliminates running errors caused by too many zeros
  pp[["marks"]] = pp[["marks"]] + 1e-7

  ## run trendsceek
  trendsceektest = trendsceek::trendsceek_test(pp, nrand = nrand, ncores = ncores, ...)

  ## get final results
  trendsceektest = trendsceektest$supstats_wide

  return(trendsceektest)
}




#' @title spark
#' @name spark
#' @description Compute spatially expressed genes with SPARK method
#' @param gobject giotto object
#' @param percentage The percentage of cells that are expressed for analysis
#' @param min_count minimum number of counts for a gene to be included
#' @param expression_values type of values to use (raw by default)
#' @param num_core number of cores to use
#' @param covariates The covariates in experiments, i.e. confounding factors/batch effect. Column name of giotto cell metadata.
#' @param return_object type of result to return (data.table or spark object)
#' @param \dots Additional parameters to the \code{\link[SPARK]{spark.vc}} function
#' @return data.table with SPARK spatial genes results or the SPARK object
#' @details This function is a wrapper for the method implemented in the SPARK package:
#' \itemize{
#'  \item{1. CreateSPARKObject }{create a SPARK object from a Giotto object}
#'  \item{2. spark.vc }{ Fits the count-based spatial model to estimate the parameters,
#'  see \code{\link[SPARK]{spark.vc}} for additional parameters}
#'  \item{3. spark.test }{ Testing multiple kernel matrices}
#' }
#' @export
spark = function(gobject,
                 percentage = 0.1,
                 min_count = 10,
                 expression_values = 'raw',
                 num_core = 5,
                 covariates = NULL,
                 return_object = c('data.table', 'spark'),
                 ...) {


  # determine parameter
  return_object = match.arg(return_object, c('data.table', 'spark'))

  # data.table variables
  genes =  adjusted_pvalue = combined_pvalue = NULL

  ## test if SPARK is installed ##
  package_check(pkg_name = 'SPARK',
                repository = c('github'),
                github_repo = 'xzhoulab/SPARK')


  # print message with information #
  message("using 'SPARK' for spatial gene/pattern detection. If used in published research, please cite:
  Sun, Shiquan, Jiaqiang Zhu, and Xiang Zhou. 'Statistical Analysis of Spatial Expression Pattern for Spatially Resolved Transcriptomic Studies.'
          BioRxiv, October 21, 2019, 810903. https://doi.org/10.1101/810903.")


  ## extract expression values from gobject
  expr = select_expression_values(gobject = gobject, values = expression_values)

  ## extract coordinates from gobject
  locs = as.data.frame(gobject@spatial_locs)
  rownames(locs) = colnames(expr)

  ## create SPARK object for analysis and filter out lowly expressed genes
  sobject = SPARK::CreateSPARKObject(counts = expr,
                                     location = locs[,1:2],
                                     percentage = percentage,
                                     min_total_counts = min_count)

  ## total counts for each cell
  sobject@lib_size = apply(sobject@counts, 2, sum)

  ## extract covariates ##
  if(!is.null(covariates)) {

    # first filter giotto object based on spark object
    filter_cell_ids = colnames(sobject@counts)
    filter_gene_ids = rownames(sobject@counts)
    tempgobject = subsetGiotto(gobject, cell_ids = filter_cell_ids, gene_ids = filter_gene_ids)

    metadata = pDataDT(tempgobject)

    if(!covariates %in% colnames(metadata)) {
      warning(covariates, ' was not found in the cell metadata of the giotto object, will be set to NULL \n')
      covariates = NULL
    } else {
      covariates = metadata[[covariates]]
    }
  }

  ## Fit statistical model under null hypothesis
  sobject = SPARK::spark.vc(sobject,
                            covariates = covariates,
                            lib_size = sobject@lib_size,
                            num_core = num_core,
                            verbose = F,
                            ...)

  ## test spatially expressed pattern genes
  ## calculating pval
  sobject = SPARK::spark.test(sobject,
                              check_positive = T,
                              verbose = F)

  ## return results ##
  if(return_object == 'spark'){
    return(sobject)
  }else if(return_object == 'data.table'){
    DT_results = data.table::as.data.table(sobject@res_mtest)
    gene_names = rownames(sobject@counts)
    DT_results[, genes := gene_names]
    data.table::setorder(DT_results, adjusted_pvalue, combined_pvalue)
    return(DT_results)
  }
}





# * ####
## PCA spatial patterns ####

#' @title detectSpatialPatterns
#' @name detectSpatialPatterns
#' @description Identify spatial patterns through PCA on average expression in a spatial grid.
#' @param gobject giotto object
#' @param expression_values expression values to use
#' @param spatial_grid_name name of spatial grid to use (default = 'spatial_grid')
#' @param min_cells_per_grid minimum number of cells in a grid to be considered
#' @param scale_unit scale features
#' @param ncp number of principal components to calculate
#' @param show_plot show plots
#' @param PC_zscore minimum z-score of variance explained by a PC
#' @return spatial pattern object 'spatPatObj'
#' @details
#' Steps to identify spatial patterns:
#' \itemize{
#'   \item{1. average gene expression for cells within a grid, see createSpatialGrid}
#'   \item{2. perform PCA on the average grid expression profiles}
#'   \item{3. convert variance of principlal components (PCs) to z-scores and select PCs based on a z-score threshold}
#' }
#' @export
detectSpatialPatterns <- function(gobject,
                                  expression_values = c('normalized', 'scaled', 'custom'),
                                  spatial_grid_name = 'spatial_grid',
                                  min_cells_per_grid = 4,
                                  scale_unit = F,
                                  ncp = 100,
                                  show_plot = T,
                                  PC_zscore = 1.5) {


  # expression values to be used
  values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))
  expr_values = select_expression_values(gobject = gobject, values = values)


  # spatial grid and spatial locations
  if(is.null(gobject@spatial_grid)) {
    stop("\n you need to create a spatial grid, see createSpatialGrid(), for this function to work \n")
  }
  if(!spatial_grid_name %in% names(gobject@spatial_grid)) {
    stop("\n you need to provide an existing spatial grid name for this function to work \n")
  }

  #spatial_grid = gobject@spatial_grid[[spatial_grid_name]]
  spatial_grid = select_spatialGrid(gobject, spatial_grid_name)

  # annotate spatial locations with spatial grid information
  spatial_locs = copy(gobject@spatial_locs)

  if(all(c('sdimx', 'sdimy', 'sdimz') %in% colnames(spatial_locs))) {
    spatial_locs = annotate_spatlocs_with_spatgrid_3D(spatloc = spatial_locs, spatgrid = spatial_grid)
  } else if(all(c('sdimx', 'sdimy') %in% colnames(spatial_locs))) {
    spatial_locs = annotate_spatlocs_with_spatgrid_2D(spatloc = spatial_locs, spatgrid = spatial_grid)
  }


  # data.table variables
  gr_loc = zscore = variance.percent = loc_ID = gene_ID = NULL

  # filter grid, minimum number of cells per grid
  cells_per_grid = sort(table(spatial_locs$gr_loc))
  cells_per_grid = cells_per_grid[cells_per_grid >= min_cells_per_grid]
  loc_names = names(cells_per_grid)

  # average expression per grid
  loc_av_expr_list <- list()
  for(loc_name in loc_names) {

    loc_cell_IDs = spatial_locs[gr_loc == loc_name]$cell_ID
    subset_expr = expr_values[, colnames(expr_values) %in% loc_cell_IDs]
    if(is.vector(subset_expr) == TRUE) {
      loc_av_expr = subset_expr
    } else {
      loc_av_expr = rowMeans(subset_expr)
    }
    loc_av_expr_list[[loc_name]] <- loc_av_expr
  }
  loc_av_expr_matrix = do.call('cbind', loc_av_expr_list)

  # START TEST
  loc_av_expr_matrix = as.matrix(loc_av_expr_matrix)
  # STOP

  # perform pca on grid matrix
  mypca <- FactoMineR::PCA(X = t(loc_av_expr_matrix), scale.unit = scale_unit, ncp = ncp, graph = F)

  # screeplot
  screeplot = factoextra::fviz_eig(mypca, addlabels = T, ylim = c(0, 50))
  if(show_plot == TRUE) {
    print(screeplot)
  }

  # select variable PCs
  eig.val <- factoextra::get_eigenvalue(mypca)
  eig.val_DT <- data.table::as.data.table(eig.val)
  eig.val_DT$names = rownames(eig.val)
  eig.val_DT[, zscore := scale(variance.percent)]
  eig.val_DT[, rank := rank(variance.percent)]
  dims_to_keep = eig.val_DT[zscore > PC_zscore]$names


  # if no dimensions are kept, return message
  if(is.null(dims_to_keep) | length(dims_to_keep) < 1) {
    return(cat('\n no PC dimensions retained, lower the PC zscore \n'))
  }

  # coordinates for cells
  pca_matrix <- mypca$ind$coord
  if(length(dims_to_keep) == 1) {
    pca_matrix_DT = data.table::data.table('dimkeep' = pca_matrix[,1],
                                           loc_ID = colnames(loc_av_expr_matrix))
    data.table::setnames(pca_matrix_DT, old = 'dimkeep', dims_to_keep)
  } else {
    pca_matrix_DT <- data.table::as.data.table(pca_matrix[,1:length(dims_to_keep)])
    pca_matrix_DT[, loc_ID := colnames(loc_av_expr_matrix)]
  }


  # correlation of genes with PCs
  feat_matrix <- mypca$var$cor
  if(length(dims_to_keep) == 1) {
    feat_matrix_DT = data.table::data.table('featkeep' = feat_matrix[,1],
                                            gene_ID = rownames(loc_av_expr_matrix))
    data.table::setnames(feat_matrix_DT, old = 'featkeep', dims_to_keep)
  } else {
    feat_matrix_DT <- data.table::as.data.table(feat_matrix[,1:length(dims_to_keep)])
    feat_matrix_DT[, gene_ID := rownames(loc_av_expr_matrix)]
  }


  spatPatObject = list(pca_matrix_DT = pca_matrix_DT,
                       feat_matrix_DT = feat_matrix_DT,
                       spatial_grid = spatial_grid)

  class(spatPatObject) <- append(class(spatPatObject), 'spatPatObj')

  return(spatPatObject)
}



#' @title showPattern2D
#' @name showPattern2D
#' @description show patterns for 2D spatial data
#' @param gobject giotto object
#' @param spatPatObj Output from detectSpatialPatterns
#' @param dimension dimension to plot
#' @param trim Trim ends of the PC values.
#' @param background_color background color for plot
#' @param grid_border_color color for grid
#' @param show_legend show legend of ggplot
#' @param point_size size of points
#' @param show_plot show plot
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters, see \code{\link{showSaveParameters}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return ggplot
#' @export
showPattern2D <- function(gobject,
                          spatPatObj,
                          dimension = 1,
                          trim = c(0.02, 0.98),
                          background_color = 'white',
                          grid_border_color = 'grey',
                          show_legend = T,
                          point_size = 1,
                          show_plot = NA,
                          return_plot = NA,
                          save_plot = NA,
                          save_param =  list(),
                          default_save_name = 'showPattern2D') {

  if(!'spatPatObj' %in% class(spatPatObj)) {
    stop('\n spatPatObj needs to be the output from detectSpatialPatterns \n')
  }

  # select PC and subset data
  selected_PC = paste0('Dim.', dimension)
  PC_DT = spatPatObj$pca_matrix_DT
  if(!selected_PC %in% colnames(PC_DT)) {
    stop('\n This dimension was not found in the spatial pattern object \n')
  }
  PC_DT = PC_DT[,c(selected_PC, 'loc_ID'), with = F]

  # annotate grid with PC values
  annotated_grid = merge(spatPatObj$spatial_grid, by.x = 'gr_name', PC_DT, by.y = 'loc_ID')

  # trim PC values
  if(!is.null(trim)) {
    boundaries = stats::quantile(annotated_grid[[selected_PC]], probs = trim)
    annotated_grid[[selected_PC]][annotated_grid[[selected_PC]] < boundaries[1]] = boundaries[1]
    annotated_grid[[selected_PC]][annotated_grid[[selected_PC]] > boundaries[2]] = boundaries[2]

  }

  # 2D-plot
  #


  dpl <- ggplot2::ggplot()
  dpl <- dpl + ggplot2::theme_bw()
  dpl <- dpl + ggplot2::geom_tile(data = annotated_grid,
                                  aes_string(x = 'x_start', y = 'y_start', fill = selected_PC),
                                  color = grid_border_color, show.legend = show_legend)
  dpl <- dpl + ggplot2::scale_fill_gradient2('low' = 'darkblue', mid = 'white', high = 'darkred', midpoint = 0,
                                             guide = guide_legend(title = ''))
  dpl <- dpl + ggplot2::theme(axis.text.x = element_text(size = 8, angle = 45, vjust = 1, hjust = 1),
                              panel.background = element_rect(fill = background_color),
                              panel.grid = element_blank(),
                              plot.title = element_text(hjust = 0.5))
  dpl <- dpl + ggplot2::labs(x = 'x coordinates', y = 'y coordinates')


  # print, return and save parameters
  show_plot = ifelse(is.na(show_plot), readGiottoInstructions(gobject, param = 'show_plot'), show_plot)
  save_plot = ifelse(is.na(save_plot), readGiottoInstructions(gobject, param = 'save_plot'), save_plot)
  return_plot = ifelse(is.na(return_plot), readGiottoInstructions(gobject, param = 'return_plot'), return_plot)

  ## print plot
  if(show_plot == TRUE) {
    print(dpl)
  }

  ## save plot
  if(save_plot == TRUE) {
    do.call('all_plots_save_function', c(list(gobject = gobject, plot_object = dpl, default_save_name = default_save_name), save_param))
  }

  ## return plot
  if(return_plot == TRUE) {
    return(dpl)
  }

}

#' @title showPattern
#' @name showPattern
#' @description show patterns for 2D spatial data
#' @param gobject giotto object
#' @param spatPatObj Output from detectSpatialPatterns
#' @inheritDotParams showPattern2D -gobject -spatPatObj
#' @return ggplot
#' @seealso \code{\link{showPattern2D}}
#' @export
showPattern = function(gobject, spatPatObj, ...) {

  showPattern2D(gobject = gobject, spatPatObj = spatPatObj, ...)

}

#' @title showPattern3D
#' @name showPattern3D
#' @description show patterns for 3D spatial data
#' @param gobject giotto object
#' @param spatPatObj Output from detectSpatialPatterns
#' @param dimension dimension to plot
#' @param trim Trim ends of the PC values.
#' @param background_color background color for plot
#' @param grid_border_color color for grid
#' @param show_legend show legend of plot
#' @param point_size adjust the point size
#' @param axis_scale scale the axis
#' @param custom_ratio cutomize the scale of the axis
#' @param x_ticks the tick number of x_axis
#' @param y_ticks the tick number of y_axis
#' @param z_ticks the tick number of z_axis
#' @param show_plot show plot
#' @param return_plot return plot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters, see \code{\link{showSaveParameters}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return plotly
#' @export
showPattern3D <- function(gobject,
                          spatPatObj,
                          dimension = 1,
                          trim = c(0.02, 0.98),
                          background_color = 'white',
                          grid_border_color = 'grey',
                          show_legend = T,
                          point_size = 1,
                          axis_scale = c("cube","real","custom"),
                          custom_ratio = NULL,
                          x_ticks = NULL,
                          y_ticks = NULL,
                          z_ticks = NULL,
                          show_plot = NA,
                          return_plot = NA,
                          save_plot = NA,
                          save_param =  list(),
                          default_save_name = 'showPattern3D') {

  # data.table variables
  center_x = x_start = x_end = center_y = y_start = y_end = center_z = z_start = z_end = NULL

  if(!'spatPatObj' %in% class(spatPatObj)) {
    stop('\n spatPatObj needs to be the output from detectSpatialPatterns \n')
  }

  # select PC and subset data
  selected_PC = paste0('Dim.', dimension)
  PC_DT = spatPatObj$pca_matrix_DT
  if(!selected_PC %in% colnames(PC_DT)) {
    stop('\n This dimension was not found in the spatial pattern object \n')
  }
  PC_DT = PC_DT[,c(selected_PC, 'loc_ID'), with = F]

  # annotate grid with PC values
  annotated_grid = merge(spatPatObj$spatial_grid, by.x = 'gr_name', PC_DT, by.y = 'loc_ID')

  # trim PC values
  if(!is.null(trim)) {
    boundaries = stats::quantile(annotated_grid[[selected_PC]], probs = trim)
    annotated_grid[[selected_PC]][annotated_grid[[selected_PC]] < boundaries[1]] = boundaries[1]
    annotated_grid[[selected_PC]][annotated_grid[[selected_PC]] > boundaries[2]] = boundaries[2]

  }


  annotated_grid <- data.table(annotated_grid)
  annotated_grid[,center_x:=(x_start+x_end)/2]
  annotated_grid[,center_y:=(y_start+y_end)/2]
  annotated_grid[,center_z:=(z_start+z_end)/2]


  axis_scale = match.arg(axis_scale, c("cube","real","custom"))

  ratio = plotly_axis_scale_3D(annotated_grid,sdimx = "center_x",sdimy = "center_y",sdimz = "center_z",
                               mode = axis_scale,custom_ratio = custom_ratio)

  dpl <- plotly::plot_ly(type = 'scatter3d',
                         x = annotated_grid$center_x, y = annotated_grid$center_y, z = annotated_grid$center_z,
                         color = annotated_grid[[selected_PC]],marker = list(size = point_size),
                         mode = 'markers', colors = c( 'darkblue','white','darkred'))
  dpl <- dpl %>% plotly::layout(scene = list(
    xaxis = list(title = "X",nticks = x_ticks),
    yaxis = list(title = "Y",nticks = y_ticks),
    zaxis = list(title = "Z",nticks = z_ticks),
    aspectmode='manual',
    aspectratio = list(x=ratio[[1]],
                       y=ratio[[2]],
                       z=ratio[[3]])))
  dpl <- dpl %>% plotly::colorbar(title = paste(paste("dim.",dimension,sep = ""),"genes", sep = " "))

  # print, return and save parameters
  show_plot = ifelse(is.na(show_plot), readGiottoInstructions(gobject, param = 'show_plot'), show_plot)
  save_plot = ifelse(is.na(save_plot), readGiottoInstructions(gobject, param = 'save_plot'), save_plot)
  return_plot = ifelse(is.na(return_plot), readGiottoInstructions(gobject, param = 'return_plot'), return_plot)

  ## print plot
  if(show_plot == TRUE) {
    print(dpl)
  }

  ## save plot
  if(save_plot == TRUE) {
    do.call('all_plots_save_function', c(list(gobject = gobject, plot_object = dpl, default_save_name = default_save_name), save_param))
  }

  ## return plot
  if(return_plot == TRUE) {
    return(dpl)
  }

}




#' @title showPatternGenes
#' @name showPatternGenes
#' @description show genes correlated with spatial patterns
#' @param gobject giotto object
#' @param spatPatObj Output from detectSpatialPatterns
#' @param dimension dimension to plot genes for.
#' @param top_pos_genes Top positively correlated genes.
#' @param top_neg_genes Top negatively correlated genes.
#' @param point_size size of points
#' @param return_DT if TRUE, it will return the data.table used to generate the plots
#' @param show_plot show plot
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters, see \code{\link{showSaveParameters}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return ggplot
#' @export
showPatternGenes <- function(gobject,
                             spatPatObj,
                             dimension = 1,
                             top_pos_genes = 5,
                             top_neg_genes = 5,
                             point_size = 1,
                             return_DT = FALSE,
                             show_plot = NA,
                             return_plot = NA,
                             save_plot = NA,
                             save_param =  list(),
                             default_save_name = 'showPatternGenes') {

  # data.table variables
  gene_ID = NULL

  if(!'spatPatObj' %in% class(spatPatObj)) {
    stop('\n spatPatObj needs to be the output from detectSpatialPatterns \n')
  }


  # select PC to use
  selected_PC = paste0('Dim.', dimension)

  gene_cor_DT = spatPatObj$feat_matrix_DT
  if(!selected_PC %in% colnames(gene_cor_DT)) {
    stop('\n This dimension was not found in the spatial pattern object \n')
  }
  gene_cor_DT = gene_cor_DT[,c(selected_PC, 'gene_ID'), with = F]

  # order and subset
  gene_cor_DT = gene_cor_DT[!is.na(get(selected_PC))][order(get(selected_PC))]

  subset = gene_cor_DT[c(1:top_neg_genes, (nrow(gene_cor_DT)-top_pos_genes):nrow(gene_cor_DT))]
  subset[, gene_ID := factor(gene_ID, gene_ID)]

  ## return DT and make not plot ##
  if(return_DT == TRUE) {
    return(subset)
  }

  pl <- ggplot()
  pl <- pl + ggplot2::theme_classic()
  pl <- pl + ggplot2::geom_point(data = subset, aes_string(x = selected_PC, y = 'gene_ID'), size = point_size)
  pl <- pl + ggplot2::geom_vline(xintercept = 0, linetype = 2)
  pl <- pl + ggplot2::labs(x = 'correlation', y = '', title = selected_PC)
  pl <- pl + ggplot2::theme(plot.title = element_text(hjust = 0.5))


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


#' @title selectPatternGenes
#' @name selectPatternGenes
#' @description Select genes correlated with spatial patterns
#' @param spatPatObj Output from detectSpatialPatterns
#' @param dimensions dimensions to identify correlated genes for.
#' @param top_pos_genes Top positively correlated genes.
#' @param top_neg_genes Top negatively correlated genes.
#' @param min_pos_cor Minimum positive correlation score to include a gene.
#' @param min_neg_cor Minimum negative correlation score to include a gene.
#' @param return_top_selection only return selection based on correlation criteria (boolean)
#' @return Data.table with genes associated with selected dimension (PC).
#' @details Description.
#' @export
selectPatternGenes <- function(spatPatObj,
                               dimensions = 1:5,
                               top_pos_genes = 10,
                               top_neg_genes = 10,
                               min_pos_cor = 0.5,
                               min_neg_cor = -0.5,
                               return_top_selection = FALSE) {


  if(!'spatPatObj' %in% class(spatPatObj)) {
    stop('\n spatPatObj needs to be the output from detectSpatialPatterns \n')
  }

  # data.table variables
  top_pos_rank = value = top_neg_rank = topvalue = gene_ID = variable = NULL


  # select PC to use
  selected_PCs = paste0('Dim.', dimensions)
  gene_cor_DT = spatPatObj$feat_matrix_DT
  if(any(selected_PCs %in% colnames(gene_cor_DT) == F)) {
    stop('\n not all dimensions were found back \n')
  }
  gene_cor_DT = gene_cor_DT[,c(selected_PCs, 'gene_ID'), with = FALSE]

  # melt and select
  gene_cor_DT_m = data.table::melt.data.table(gene_cor_DT, id.vars = 'gene_ID')
  gene_cor_DT_m[, top_pos_rank := rank(value), by = 'variable']
  gene_cor_DT_m[, top_neg_rank := rank(-value), by = 'variable']
  selection = gene_cor_DT_m[top_pos_rank %in% 1:top_pos_genes | top_neg_rank %in% 1:top_neg_genes]

  # filter on min correlation
  selection = selection[value > min_pos_cor | value < min_neg_cor]

  # return all the top correlated genes + information
  if(return_top_selection == TRUE) {
    return(selection)
  }

  # remove duplicated genes by only retaining the most correlated dimension
  selection[, topvalue := max(abs(value)), by = 'gene_ID']
  uniq_selection = selection[value == topvalue]

  # add other genes back
  output_selection = uniq_selection[,.(gene_ID, variable)]
  other_genes = gene_cor_DT[!gene_ID %in% output_selection$gene_ID][['gene_ID']]
  other_genes_DT = data.table::data.table(gene_ID = other_genes, variable = 'noDim')

  comb_output_genes = rbind(output_selection, other_genes_DT)
  setnames(comb_output_genes, 'variable', 'patDim')

  return(comb_output_genes)

}







# ** ####
## Spatial co-expression ####
## ----------- ##

#' @title do_spatial_knn_smoothing
#' @name do_spatial_knn_smoothing
#' @description smooth gene expression over a kNN spatial network
#' @param gobject giotto object
#' @param expression_values gene expression values to use
#' @param subset_genes subset of genes to use
#' @param spatial_network_name name of spatial network to use
#' @param b smoothing factor beteen 0 and 1 (default: automatic)
#' @return matrix with smoothened gene expression values based on kNN spatial network
#' @details This function will smoothen the gene expression values per cell according to
#' its neighbors in the selected spatial network. \cr
#' b is a smoothening factor that defaults to 1 - 1/k, where k is the median number of
#' k-neighbors in the selected spatial network. Setting b = 0 means no smoothing and b = 1
#' means no contribution from its own expression.
#' @keywords internal
do_spatial_knn_smoothing = function(gobject,
                                    expression_values = c('normalized', 'scaled', 'custom'),
                                    subset_genes = NULL,
                                    spatial_network_name = 'Delaunay_network',
                                    b = NULL) {

  # checks
  if(!is.null(b)) {
    if(b > 1 | b < 0) {
      stop('b needs to be between 0 (no spatial contribution) and 1 (only spatial contribution)')
    }
  }

  # get spatial network
  spatial_network = select_spatialNetwork(gobject,name = spatial_network_name,return_network_Obj = FALSE)

  # get expression matrix
  values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))
  expr_values = select_expression_values(gobject = gobject, values = values)

  if(!is.null(subset_genes)) {
    expr_values = expr_values[rownames(expr_values) %in% subset_genes,]
  }

  # data.table variables
  gene_ID = value = NULL

  # merge spatial network with expression data
  expr_values_dt = data.table::as.data.table(expr_values); expr_values_dt[, gene_ID := rownames(expr_values)]
  expr_values_dt_m = data.table::melt.data.table(expr_values_dt, id.vars = 'gene_ID', variable.name = 'cell_ID')

  ## test ##
  spatial_network = convert_to_full_spatial_network(spatial_network)
  ## stop test ##

  #print(spatial_network)

  spatial_network_ext = data.table::merge.data.table(spatial_network, expr_values_dt_m, by.x = 'target', by.y = 'cell_ID', allow.cartesian = T)

  #print(spatial_network_ext)

  # calculate mean over all k-neighbours
  # exclude 0's?
  # trimmed mean?
  spatial_network_ext_smooth = spatial_network_ext[, mean(value), by = c('source', 'gene_ID')]

  # convert back to matrix
  spatial_smooth_dc = dcast.data.table(data = spatial_network_ext_smooth, formula = gene_ID~source, value.var = 'V1')
  spatial_smooth_matrix = dt_to_matrix(spatial_smooth_dc)


  # if network was not fully connected, some cells might be missing and are not smoothed
  # add the original values for those cells back
  all_cells = colnames(expr_values)
  smoothed_cells = colnames(spatial_smooth_matrix)
  missing_cells = all_cells[!all_cells %in% smoothed_cells]
  if(length(missing_cells) > 0) {
    missing_matrix = expr_values[, missing_cells]
    spatial_smooth_matrix = cbind(spatial_smooth_matrix[rownames(expr_values),], missing_matrix)
  }

  spatial_smooth_matrix = spatial_smooth_matrix[rownames(expr_values), colnames(expr_values)]

  # combine original and smoothed values according to smoothening b
  # create best guess for b if not given
  if(is.null(b)) {
    k = stats::median(table(spatial_network$source))
    smooth_b = 1 - 1/k
  } else {
    smooth_b = b
  }

  expr_b = 1 - smooth_b
  spatsmooth_expr_values = ((smooth_b*spatial_smooth_matrix) + (expr_b*expr_values))

  return(spatsmooth_expr_values)

}


#' @title do_spatial_grid_averaging
#' @name do_spatial_grid_averaging
#' @description smooth gene expression over a defined spatial grid
#' @param gobject giotto object
#' @param expression_values gene expression values to use
#' @param subset_genes subset of genes to use
#' @param spatial_grid_name name of spatial grid to use
#' @param min_cells_per_grid minimum number of cells to consider a grid
#' @return matrix with smoothened gene expression values based on spatial grid
#' @keywords internal
do_spatial_grid_averaging = function(gobject,
                                     expression_values = c('normalized', 'scaled', 'custom'),
                                     subset_genes = NULL,
                                     spatial_grid_name = 'spatial_grid',
                                     min_cells_per_grid = 4) {


  # get expression matrix
  values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))
  expr_values = select_expression_values(gobject = gobject, values = values)

  if(!is.null(subset_genes)) {
    expr_values = expr_values[rownames(expr_values) %in% subset_genes,]
  }

  # spatial grid and spatial locations
  if(is.null(gobject@spatial_grid)) {
    stop("\n you need to create a spatial grid, see createSpatialGrid(), for this function to work \n")
  }
  if(!spatial_grid_name %in% names(gobject@spatial_grid)) {
    stop("\n you need to provide an existing spatial grid name for this function to work \n")
  }

  #spatial_grid = gobject@spatial_grid[[spatial_grid_name]]
  spatial_grid = select_spatialGrid(gobject, spatial_grid_name)


  # annotate spatial locations with spatial grid information
  spatial_locs = copy(gobject@spatial_locs)

  if(all(c('sdimx', 'sdimy', 'sdimz') %in% colnames(spatial_locs))) {
    spatial_locs = annotate_spatlocs_with_spatgrid_3D(spatloc = spatial_locs, spatgrid = spatial_grid)
  } else if(all(c('sdimx', 'sdimy') %in% colnames(spatial_locs))) {
    spatial_locs = annotate_spatlocs_with_spatgrid_2D(spatloc = spatial_locs, spatgrid = spatial_grid)
  }


  # data.table variables
  gr_loc = NULL

  # filter grid, minimum number of cells per grid
  cells_per_grid = sort(table(spatial_locs$gr_loc))
  cells_per_grid = cells_per_grid[cells_per_grid >= min_cells_per_grid]
  loc_names = names(cells_per_grid)

  # average expression per grid
  loc_av_expr_list <- list()
  for(loc_name in loc_names) {

    loc_cell_IDs = spatial_locs[gr_loc == loc_name]$cell_ID
    subset_expr = expr_values[, colnames(expr_values) %in% loc_cell_IDs]
    if(is.vector(subset_expr) == TRUE) {
      loc_av_expr = subset_expr
    } else {
      loc_av_expr = rowMeans(subset_expr)
    }
    loc_av_expr_list[[loc_name]] <- loc_av_expr
  }
  loc_av_expr_matrix = do.call('cbind', loc_av_expr_list)
  loc_av_expr_matrix = as.matrix(loc_av_expr_matrix)

  return(loc_av_expr_matrix)
}


#' @title detectSpatialCorGenes
#' @name detectSpatialCorGenes
#' @description Detect genes that are spatially correlated
#' @param gobject giotto object
#' @param method method to use for spatial averaging
#' @param expression_values gene expression values to use
#' @param subset_genes subset of genes to use
#' @param spatial_network_name name of spatial network to use
#' @param network_smoothing  smoothing factor beteen 0 and 1 (default: automatic)
#' @param spatial_grid_name name of spatial grid to use
#' @param min_cells_per_grid minimum number of cells to consider a grid
#' @param cor_method correlation method
#' @return returns a spatial correlation object: "spatCorObject"
#' @details
#' For method = network, it expects a fully connected spatial network. You can make sure to create a
#' fully connected network by setting minimal_k > 0 in the \code{\link{createSpatialNetwork}} function.
#' \itemize{
#'  \item{1. grid-averaging: }{average gene expression values within a predefined spatial grid}
#'  \item{2. network-averaging: }{smoothens the gene expression matrix by averaging the expression within one cell
#'  by using the neighbours within the predefined spatial network. b is a smoothening factor
#'  that defaults to 1 - 1/k, where k is the median number of  k-neighbors in the
#'  selected spatial network. Setting b = 0 means no smoothing and b = 1 means no contribution
#'  from its own expression.}
#' }
#' The spatCorObject can be further explored with showSpatialCorGenes()
#' @seealso \code{\link{showSpatialCorGenes}}
#' @export
detectSpatialCorGenes <- function(gobject,
                                  method = c('grid', 'network'),
                                  expression_values = c('normalized', 'scaled', 'custom'),
                                  subset_genes = NULL,
                                  spatial_network_name = 'Delaunay_network',
                                  network_smoothing = NULL,
                                  spatial_grid_name = 'spatial_grid',
                                  min_cells_per_grid = 4,
                                  cor_method = c('pearson', 'kendall', 'spearman')) {


  ## correlation method to be used
  cor_method = match.arg(cor_method, choices = c('pearson', 'kendall', 'spearman'))

  ## method to be used
  method = match.arg(method, choices = c('grid', 'network'))

  # get expression matrix
  values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))
  expr_values = select_expression_values(gobject = gobject, values = values)

  if(!is.null(subset_genes)) {
    expr_values = expr_values[rownames(expr_values) %in% subset_genes,]
  }


  ## spatial averaging or smoothing
  if(method == 'grid') {

    loc_av_expr_matrix = do_spatial_grid_averaging(gobject = gobject,
                                                            expression_values = expression_values,
                                                            subset_genes = subset_genes,
                                                            spatial_grid_name = spatial_grid_name,
                                                            min_cells_per_grid = min_cells_per_grid)

    # data.table variables
    gene_ID = variable = NULL

    cor_spat_matrix = cor_giotto(t_giotto(as.matrix(loc_av_expr_matrix)), method = cor_method)
    cor_spat_matrixDT = data.table::as.data.table(cor_spat_matrix)
    cor_spat_matrixDT[, gene_ID := rownames(cor_spat_matrix)]
    cor_spat_DT = data.table::melt.data.table(data = cor_spat_matrixDT,
                                              id.vars = 'gene_ID', value.name = 'spat_cor')
  }

  if(method == 'network') {

    knn_av_expr_matrix = do_spatial_knn_smoothing(gobject = gobject,
                                                  expression_values = expression_values,
                                                  subset_genes = subset_genes,
                                                  spatial_network_name = spatial_network_name,
                                                  b = network_smoothing)

    #print(knn_av_expr_matrix[1:4, 1:4])

    cor_spat_matrix = cor_giotto(t_giotto(as.matrix(knn_av_expr_matrix)), method = cor_method)
    cor_spat_matrixDT = data.table::as.data.table(cor_spat_matrix)
    cor_spat_matrixDT[, gene_ID := rownames(cor_spat_matrix)]
    cor_spat_DT = data.table::melt.data.table(data = cor_spat_matrixDT,
                                              id.vars = 'gene_ID', value.name = 'spat_cor')


  }



  # data.table variables
  cordiff = spat_cor = expr_cor = spatrank= exprrank = rankdiff = NULL

  ## 2. perform expression correlation at single-cell level without spatial information
  cor_matrix = cor_giotto(t_giotto(expr_values), method = cor_method)
  cor_matrixDT = data.table::as.data.table(cor_matrix)
  cor_matrixDT[, gene_ID := rownames(cor_matrix)]
  cor_DT = data.table::melt.data.table(data = cor_matrixDT,
                                       id.vars = 'gene_ID', value.name = 'expr_cor')

  ## 3. merge spatial and expression correlation
  data.table::setorder(cor_spat_DT, gene_ID, variable)
  data.table::setorder(cor_DT, gene_ID, variable)
  doubleDT = cbind(cor_spat_DT, expr_cor = cor_DT[['expr_cor']])

  # difference in correlation scores
  doubleDT[, cordiff := spat_cor - expr_cor]

  # difference in rank scores
  doubleDT[, spatrank := frank(-spat_cor, ties.method = 'first'), by = gene_ID]
  doubleDT[, exprrank := frank(-expr_cor, ties.method = 'first'), by = gene_ID]
  doubleDT[, rankdiff := spatrank - exprrank]

  # sort data
  data.table::setorder(doubleDT, gene_ID, -spat_cor)

  spatCorObject = list(cor_DT = doubleDT,
                       gene_order = rownames(cor_spat_matrix),
                       cor_hclust = list(),
                       cor_clusters = list())

  class(spatCorObject) = append(class(spatCorObject), 'spatCorObject')

  return(spatCorObject)

}





#' @title showSpatialCorGenes
#' @name showSpatialCorGenes
#' @description Shows and filters spatially correlated genes
#' @param spatCorObject spatial correlation object
#' @param use_clus_name cluster information to show
#' @param selected_clusters subset of clusters to show
#' @param genes subset of genes to show
#' @param min_spat_cor filter on minimum spatial correlation
#' @param min_expr_cor filter on minimum single-cell expression correlation
#' @param min_cor_diff filter on minimum correlation difference (spatial vs expression)
#' @param min_rank_diff filter on minimum correlation rank difference (spatial vs expression)
#' @param show_top_genes show top genes per gene
#' @return data.table with filtered information
#' @export
showSpatialCorGenes = function(spatCorObject,
                               use_clus_name = NULL,
                               selected_clusters = NULL,
                               genes = NULL,
                               min_spat_cor = 0.5,
                               min_expr_cor = NULL,
                               min_cor_diff = NULL,
                               min_rank_diff = NULL,
                               show_top_genes = NULL) {

  # data.table variables
  clus = gene_ID = spat_cor = cor_diff = rankdiff = NULL

  if(!'spatCorObject' %in% class(spatCorObject)) {
    stop('\n spatCorObject needs to be the output from detectSpatialCorGenes() \n')
  }

  filter_DT = copy(spatCorObject[['cor_DT']])

  if(!is.null(use_clus_name)) {

    clusters_part = spatCorObject[['cor_clusters']][[use_clus_name]]

    # combine spatial correlation info and clusters
    clusters = clusters_part
    names_clusters = names(clusters_part)
    clusters_DT = data.table::data.table('gene_ID' = names_clusters, 'clus' = clusters)
    filter_DT = data.table::merge.data.table(filter_DT, clusters_DT, by = 'gene_ID')
  }

  ## 0. subset clusters
  if(!is.null(selected_clusters)) {
    filter_DT = filter_DT[clus %in% selected_clusters]
  }


  ## 1. subset genes
  if(!is.null(genes)) {
    filter_DT = filter_DT[gene_ID %in% genes]
  }

  ## 2. select spatial correlation
  if(!is.null(min_spat_cor)) {
    filter_DT = filter_DT[spat_cor >= min_spat_cor]
  }

  ## 3. minimum expression correlation
  if(!is.null(min_expr_cor)) {
    filter_DT = filter_DT[spat_cor >= min_expr_cor]
  }

  ## 4. minimum correlation difference
  if(!is.null(min_cor_diff)) {
    filter_DT = filter_DT[cor_diff >= min_cor_diff]
  }

  ## 5. minimum correlation difference
  if(!is.null(min_rank_diff)) {
    filter_DT = filter_DT[rankdiff >= min_rank_diff]
  }

  ## 6. show only top genes
  if(!is.null(show_top_genes)) {
    filter_DT = filter_DT[, head(.SD, show_top_genes), by = gene_ID]
  }

  return(filter_DT)

}



#' @title clusterSpatialCorGenes
#' @name clusterSpatialCorGenes
#' @description Cluster based on spatially correlated genes
#' @param spatCorObject spatial correlation object
#' @param name name for spatial clustering results
#' @param hclust_method method for hierarchical clustering
#' @param k number of clusters to extract
#' @param return_obj return spatial correlation object (spatCorObject)
#' @return spatCorObject or cluster results
#' @export
clusterSpatialCorGenes = function(spatCorObject,
                                  name = 'spat_clus',
                                  hclust_method = 'ward.D',
                                  k = 10,
                                  return_obj = TRUE) {


  # check input
  if(!'spatCorObject' %in% class(spatCorObject)) {
    stop('\n spatCorObject needs to be the output from detectSpatialCorGenes() \n')
  }

  # create correlation matrix
  cor_DT = spatCorObject[['cor_DT']]
  cor_DT_dc = data.table::dcast.data.table(cor_DT, formula = gene_ID~variable, value.var = 'spat_cor')
  cor_matrix = dt_to_matrix(cor_DT_dc)

  # re-ordering matrix
  my_gene_order = spatCorObject[['gene_order']]
  cor_matrix = cor_matrix[my_gene_order, my_gene_order]

  # cluster
  cor_dist = stats::as.dist(1-cor_matrix)
  cor_h = stats::hclust(d = cor_dist, method = hclust_method)
  cor_clus = stats::cutree(cor_h, k = k)

  if(return_obj == TRUE) {
    spatCorObject[['cor_hclust']][[name]] = cor_h
    spatCorObject[['cor_clusters']][[name]] = cor_clus
    spatCorObject[['cor_coexpr_groups']][[name]] = NA

    return(spatCorObject)

  } else {
    return(list('hclust' = cor_h, 'clusters' = cor_clus))
  }

}



#' @title heatmSpatialCorGenes
#' @name heatmSpatialCorGenes
#' @description Create heatmap of spatially correlated genes
#' @param gobject giotto object
#' @param spatCorObject spatial correlation object
#' @param use_clus_name name of clusters to visualize (from clusterSpatialCorGenes())
#' @param show_cluster_annot show cluster annotation on top of heatmap
#' @param show_row_dend show row dendrogram
#' @param show_column_dend show column dendrogram
#' @param show_row_names show row names
#' @param show_column_names show column names
#' @param show_plot show plot
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters, see \code{\link{showSaveParameters}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @param \dots additional parameters to the \code{\link[ComplexHeatmap]{Heatmap}} function from ComplexHeatmap
#' @return Heatmap generated by ComplexHeatmap
#' @export
heatmSpatialCorGenes = function(gobject,
                                spatCorObject,
                                use_clus_name = NULL,
                                show_cluster_annot = TRUE,
                                show_row_dend = T,
                                show_column_dend = F,
                                show_row_names = F,
                                show_column_names = F,
                                show_plot = NA,
                                return_plot = NA,
                                save_plot = NA,
                                save_param =  list(),
                                default_save_name = 'heatmSpatialCorGenes',
                                ...) {

  ## check input
  if(!'spatCorObject' %in% class(spatCorObject)) {
    stop('\n spatCorObject needs to be the output from detectSpatialCorGenes() \n')
  }

  ## create correlation matrix
  cor_DT = spatCorObject[['cor_DT']]
  cor_DT_dc = data.table::dcast.data.table(cor_DT, formula = gene_ID~variable, value.var = 'spat_cor')
  cor_matrix = dt_to_matrix(cor_DT_dc)

  # re-ordering matrix
  my_gene_order = spatCorObject[['gene_order']]
  cor_matrix = cor_matrix[my_gene_order, my_gene_order]


  ## fix row and column names
  cor_matrix = cor_matrix[rownames(cor_matrix), rownames(cor_matrix)]

  ## default top annotation
  ha = NULL

  if(!is.null(use_clus_name)) {
    hclust_part = spatCorObject[['cor_hclust']][[use_clus_name]]

    if(is.null(hclust_part)) {
      cat(use_clus_name, ' does not exist, make one with spatCorCluster \n')
      hclust_part = TRUE

    } else {
      clusters_part = spatCorObject[['cor_clusters']][[use_clus_name]]

      if(show_cluster_annot) {
        uniq_clusters = unique(clusters_part)

        # color vector
        mycolors = getDistinctColors(length(uniq_clusters))
        names(mycolors) = uniq_clusters
        ha = ComplexHeatmap::HeatmapAnnotation(bar = as.vector(clusters_part),
                                               col = list(bar = mycolors),
                                               annotation_legend_param = list(title = NULL))
      }

    }
  } else {
    hclust_part = TRUE
  }


  ## create heatmap
  heatm = ComplexHeatmap::Heatmap(matrix = as.matrix(cor_matrix),
                                  cluster_rows = hclust_part,
                                  cluster_columns = hclust_part,
                                  show_row_dend = show_row_dend,
                                  show_column_dend = show_column_dend,
                                  show_row_names = show_row_names,
                                  show_column_names = show_column_names,
                                  top_annotation = ha, ...)

  # print, return and save parameters
  show_plot = ifelse(is.na(show_plot), readGiottoInstructions(gobject, param = 'show_plot'), show_plot)
  save_plot = ifelse(is.na(save_plot), readGiottoInstructions(gobject, param = 'save_plot'), save_plot)
  return_plot = ifelse(is.na(return_plot), readGiottoInstructions(gobject, param = 'return_plot'), return_plot)

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





#' @title rankSpatialCorGroups
#' @name rankSpatialCorGroups
#' @description Rank spatial correlated clusters according to correlation structure
#' @param gobject giotto object
#' @param spatCorObject spatial correlation object
#' @param use_clus_name name of clusters to visualize (from clusterSpatialCorGenes())
#' @param show_plot show plot
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters, see \code{\link{showSaveParameters}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return data.table with positive (within group) and negative (outside group) scores
#' @export
rankSpatialCorGroups = function(gobject,
                                spatCorObject,
                                use_clus_name = NULL,
                                show_plot = NA,
                                return_plot = FALSE,
                                save_plot = NA,
                                save_param =  list(),
                                default_save_name = 'rankSpatialCorGroups') {


  ## check input
  if(!'spatCorObject' %in% class(spatCorObject)) {
    stop('\n spatCorObject needs to be the output from detectSpatialCorGenes() \n')
  }

  ## check if cluster exist
  if(is.null(use_clus_name)) {
    stop('use_clus_name does not exist \n')
  }
  clusters_part = spatCorObject[['cor_clusters']][[use_clus_name]]

  if(is.null(clusters_part)) {
    stop('use_clus_name does not exist \n')
  }

  ## create correlation matrix
  cor_DT = spatCorObject[['cor_DT']]
  cor_DT_dc = data.table::dcast.data.table(cor_DT, formula = gene_ID~variable, value.var = 'spat_cor')
  cor_matrix = dt_to_matrix(cor_DT_dc)

  # re-ordering matrix
  my_gene_order = spatCorObject[['gene_order']]
  cor_matrix = cor_matrix[my_gene_order, my_gene_order]



  res_cor_list = list()
  res_neg_cor_list = list()
  nr_genes_list = list()

  for(id in 1:length(unique(clusters_part))) {

    clus_id = unique(clusters_part)[id]
    selected_genes = names(clusters_part[clusters_part == clus_id])
    nr_genes_list[[id]] = length(selected_genes)

    sub_cor_matrix = cor_matrix[rownames(cor_matrix) %in% selected_genes, colnames(cor_matrix) %in% selected_genes]
    mean_score = mean_giotto(sub_cor_matrix)
    res_cor_list[[id]] = mean_score

    sub_neg_cor_matrix = cor_matrix[rownames(cor_matrix) %in% selected_genes, !colnames(cor_matrix) %in% selected_genes]
    mean_neg_score = mean_giotto(sub_neg_cor_matrix)
    res_neg_cor_list[[id]] = mean_neg_score
  }


  # data.table variables
  cor_neg_adj = cor_neg_score = adj_cor_score = cor_score = clusters = nr_genes = NULL

  res_cor_DT = data.table::data.table('clusters' = unique(clusters_part),
                                      cor_score = unlist(res_cor_list),
                                      cor_neg_score = unlist(res_neg_cor_list),
                                      nr_genes = unlist(nr_genes_list))

  res_cor_DT[, cor_neg_adj := 1-(cor_neg_score-min(cor_neg_score))]
  res_cor_DT[, adj_cor_score := cor_neg_adj * cor_score]
  data.table::setorder(res_cor_DT, -adj_cor_score)
  res_cor_DT[, clusters := factor(x = clusters, levels = rev(clusters))]

  # print, return and save parameters
  show_plot = ifelse(is.na(show_plot), readGiottoInstructions(gobject, param = 'show_plot'), show_plot)
  save_plot = ifelse(is.na(save_plot), readGiottoInstructions(gobject, param = 'save_plot'), save_plot)
  return_plot = ifelse(is.na(return_plot), readGiottoInstructions(gobject, param = 'return_plot'), return_plot)


  pl = ggplot2::ggplot()
  pl = pl + ggplot2::geom_point(data = res_cor_DT, ggplot2::aes(x = clusters, y = adj_cor_score, size = nr_genes))
  pl = pl + ggplot2::theme_classic()
  pl = pl + ggplot2::labs(x = 'clusters', y = 'pos r x (1 - (neg_r - min(neg_r)))')


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
  } else {
    return(res_cor_DT)
  }

}





# ** ####
## Simulate single-gene spatial patterns ####
## --------------------------------------- ##



#' @title simulateOneGenePatternGiottoObject
#' @name simulateOneGenePatternGiottoObject
#' @description Create a simulated spatial pattern for one selected gnee
#' @param gobject giotto object
#' @param pattern_name name of spatial pattern
#' @param pattern_cell_ids cell ids that make up the spatial pattern
#' @param gene_name selected gene
#' @param spatial_prob probability for a high expressing gene value to be part of the spatial pattern
#' @param gradient_direction direction of gradient
#' @param show_pattern show the discrete spatial pattern
#' @param pattern_colors 2 color vector for the spatial pattern
#' @param \dots additional parameters for (re-)normalizing
#' @return Reprocessed Giotto object for which one gene has a forced spatial pattern
#' @export
simulateOneGenePatternGiottoObject = function(gobject,
                                              pattern_name = 'pattern',
                                              pattern_cell_ids = NULL,
                                              gene_name = NULL,
                                              spatial_prob = 0.95,
                                              gradient_direction = NULL,
                                              show_pattern = TRUE,
                                              pattern_colors = c('in' = 'green', 'out' = 'red'),
                                              ...) {

  # data.table variables
  cell_ID = sdimx_y = sdimx = sdimy = NULL

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






#' @title run_spatial_sim_tests_one_rep
#' @name run_spatial_sim_tests_one_rep
#' @description runs all spatial tests for 1 probability and 1 rep
#' @keywords internal
run_spatial_sim_tests_one_rep = function(gobject,
                                         pattern_name = 'pattern',
                                         pattern_cell_ids = NULL,
                                         gene_name = NULL,
                                         spatial_prob = 0.95,
                                         show_pattern = FALSE,
                                         spatial_network_name = 'kNN_network',
                                         spat_methods = c('binSpect_single', 'binSpect_multi', 'spatialDE', 'spark', 'silhouetteRank'),
                                         spat_methods_params = list(NA, NA, NA, NA, NA),
                                         spat_methods_names = c('binSpect_single', 'binSpect_multi', 'spatialDE', 'spark', 'silhouetteRank'),
                                         save_plot = FALSE,
                                         save_raw = FALSE,
                                         save_norm = FALSE,
                                         save_dir = '~',
                                         save_name = 'plot',
                                         run_simulations = TRUE,
                                         ...) {


  # data.table variables
  genes = prob = time = adj.p.value = method = p.val = sd = qval = pval = g = adjusted_pvalue = NULL

  ## test if spat_methods, params and names have the same length
  if(length(spat_methods) != length(spat_methods_params)) {
    stop('number of spatial detection methods to test need to be equal to number of spatial methods parameters \n')
  }
  if(length(spat_methods) != length(spat_methods_names)) {
    stop('number of spatial detection methods to test need to be equal to number of spatial methods names \n')
  }


  ## simulate pattern ##
  simulate_patch = simulateOneGenePatternGiottoObject(gobject,
                                                      pattern_name = pattern_name,
                                                      pattern_cell_ids = pattern_cell_ids,
                                                      gene_name = gene_name,
                                                      spatial_prob = spatial_prob,
                                                      gradient_direction = NULL,
                                                      show_pattern = show_pattern,
                                                      ...)

  # save plot
  if(save_plot == TRUE) {

    spatGenePlot2D(simulate_patch, expression_values = 'norm', genes = gene_name,
                   point_shape = 'border', point_border_stroke = 0.1, point_size = 2.5,
                   cow_n_col = 1, show_plot = F,
                   save_plot = T, save_param = list(save_dir = save_dir, save_folder = pattern_name,
                                                    save_name = save_name,
                                                    base_width = 9, base_height = 7, units = 'cm'))

  }

  # save raw data
  if(save_raw == TRUE) {

    folder_path = paste0(save_dir, '/', pattern_name)
    if(!file.exists(folder_path)) dir.create(folder_path, recursive = TRUE)

    write.table(x = as.matrix(simulate_patch@raw_exprs),
                file = paste0(save_dir, '/', pattern_name,'/',  save_name, '_raw_data.txt'),
                sep = '\t')
  }

  # save normalized data
  if(save_norm == TRUE) {

    folder_path = paste0(save_dir, '/', pattern_name)
    if(!file.exists(folder_path)) dir.create(folder_path, recursive = TRUE)

    write.table(x = as.matrix(simulate_patch@norm_expr),
                file = paste0(save_dir, '/', pattern_name,'/', save_name, '_norm_data.txt'),
                sep = '\t')
  }



  ## do simulations ##
  if(run_simulations == TRUE) {

    result_list = list()
    for(test in 1:length(spat_methods)) {

      # method
      selected_method = spat_methods[test]
      if(!selected_method %in% c('binSpect_single', 'binSpect_multi', 'spatialDE', 'spark', 'silhouetteRank')) {
        stop(selected_method, ' is not a know spatial method \n')
      }

      # params
      selected_params = spat_methods_params[[test]]

      if(length(selected_params) == 1) {

        if(is.na(selected_params)) {

          if(selected_method == 'binSpect_single') {
            selected_params = list(bin_method = 'kmeans',
                                   nstart = 3,
                                   iter_max = 10,
                                   expression_values = 'normalized',
                                   get_av_expr = FALSE,
                                   get_high_expr = FALSE)

          } else if(selected_method == 'binSpect_multi') {
            selected_params = list(bin_method = 'kmeans',
                                   spatial_network_k = c(5, 10, 20),
                                   nstart = 3,
                                   iter_max = 10,
                                   expression_values = 'normalized',
                                   get_av_expr = FALSE,
                                   get_high_expr = FALSE,
                                   summarize = 'adj.p.value')

          } else if(selected_method == 'spatialDE') {
            selected_params = list(expression_values = 'raw',
                                   sig_alpha = 0.5,
                                   unsig_alpha = 0.5,
                                   show_plot = FALSE,
                                   return_plot = FALSE,
                                   save_plot = FALSE)


          } else if(selected_method == 'spark') {
            selected_params = list(expression_values = 'raw',
                                   return_object = 'data.table',
                                   percentage = 0.1,
                                   min_count = 10,
                                   num_core = 5)

          }  else if(selected_method == 'silhouetteRank') {
            selected_params = list(expression_values = 'normalized',
                                   overwrite_input_bin = FALSE,
                                   rbp_ps = c(0.95, 0.99),
                                   examine_tops = c(0.005, 0.010),
                                   matrix_type = "dissim",
                                   num_core = 4,
                                   parallel_path = "/usr/bin",
                                   output = NULL,
                                   query_sizes = 10L)

          }

        }

      }

      # name
      selected_name = spat_methods_names[test]


      ## RUN Spatial Analysis ##
      if(selected_method == 'binSpect_single') {

        start = proc.time()
        spatial_gene_results = do.call('binSpectSingle', c(gobject =  simulate_patch,
                                                           selected_params))

        spatial_gene_results = spatial_gene_results[genes == gene_name]
        total_time = proc.time() - start

        spatial_gene_results[, prob := spatial_prob]
        spatial_gene_results[, time := total_time[['elapsed']] ]

        spatial_gene_results = spatial_gene_results[,.(genes, adj.p.value, prob, time)]
        colnames(spatial_gene_results) = c('genes', 'adj.p.value', 'prob', 'time')

        spatial_gene_results[, method := selected_name]


      } else if(selected_method == 'binSpect_multi') {

        start = proc.time()
        spatial_gene_results = do.call('binSpectMulti', c(gobject =  simulate_patch,
                                                          selected_params))

        spatial_gene_results = spatial_gene_results$simple
        spatial_gene_results = spatial_gene_results[genes == gene_name]
        total_time = proc.time() - start

        spatial_gene_results[, prob := spatial_prob]
        spatial_gene_results[, time := total_time[['elapsed']] ]

        spatial_gene_results = spatial_gene_results[,.(genes, p.val, prob, time)]
        colnames(spatial_gene_results) = c('genes', 'adj.p.value', 'prob', 'time')

        spatial_gene_results[, method := selected_name]


      } else if(selected_method == 'spatialDE') {

        start = proc.time()
        new_raw_sim_matrix = simulate_patch@raw_exprs
        sd_cells = apply(new_raw_sim_matrix, 2, sd)
        sd_non_zero_cells = names(sd_cells[sd_cells != 0])
        simulate_patch_fix = subsetGiotto(simulate_patch, cell_ids = sd_non_zero_cells)

        spatial_gene_results = do.call('spatialDE', c(gobject =  simulate_patch_fix,
                                                      selected_params))

        spatialDE_spatialgenes_sim_res = spatial_gene_results$results$results
        if(is.null(spatialDE_spatialgenes_sim_res)) spatialDE_spatialgenes_sim_res = spatial_gene_results$results

        spatialDE_spatialgenes_sim_res = data.table::as.data.table(spatialDE_spatialgenes_sim_res)
        data.table::setorder(spatialDE_spatialgenes_sim_res, qval, pval)
        spatialDE_result = spatialDE_spatialgenes_sim_res[g == gene_name]

        spatialDE_time = proc.time() - start

        spatialDE_result[, prob := spatial_prob]
        spatialDE_result[, time := spatialDE_time[['elapsed']] ]

        spatial_gene_results = spatialDE_result[,.(g, qval, prob, time)]
        colnames(spatial_gene_results) = c('genes', 'adj.p.value', 'prob', 'time')
        spatial_gene_results[, method := 'spatialDE']


      } else if(selected_method == 'spark') {

        ## spark
        start = proc.time()
        spark_spatialgenes_sim = do.call('spark', c(gobject =  simulate_patch,
                                                    selected_params))

        spark_result = spark_spatialgenes_sim[genes == gene_name]
        spark_time = proc.time() - start

        spark_result[, prob := spatial_prob]
        spark_result[, time := spark_time[['elapsed']] ]

        spatial_gene_results = spark_result[,.(genes, adjusted_pvalue, prob, time)]
        colnames(spatial_gene_results) = c('genes', 'adj.p.value', 'prob', 'time')
        spatial_gene_results[, method := 'spark']


      } else if(selected_method == 'silhouetteRank') {

        ## silhouetterank
        start = proc.time()

        spatial_gene_results = do.call('silhouetteRankTest', c(gobject = simulate_patch,
                                                           selected_params))

        data.table::setnames(spatial_gene_results, old = 'gene', new = 'genes')
        spatial_gene_results = spatial_gene_results[genes == gene_name]
        silh_time = proc.time() - start

        spatial_gene_results[, prob := spatial_prob]
        spatial_gene_results[, time := silh_time[['elapsed']] ]

        # silhrank uses qval by default
        spatial_gene_results = spatial_gene_results[,.(genes, qval, prob, time)]
        colnames(spatial_gene_results) = c('genes', 'adj.p.value', 'prob', 'time')
        spatial_gene_results[, method := 'silhouette']

      }

      result_list[[test]] = spatial_gene_results

    }

    results = data.table::rbindlist(l = result_list)
    return(results)

  } else {
    return(NULL)
  }

}





#' @title run_spatial_sim_tests_multi
#' @name run_spatial_sim_tests_multi
#' @description runs all spatial tests for multiple probabilities and repetitions
#' @keywords internal
run_spatial_sim_tests_multi = function(gobject,
                                       pattern_name = 'pattern',
                                       pattern_cell_ids = NULL,
                                       gene_name = NULL,
                                       spatial_probs = c(0.5, 1),
                                       reps = 2,

                                       spatial_network_name = 'kNN_network',
                                       spat_methods = c('binSpect_single', 'binSpect_multi', 'spatialDE', 'spark', 'silhouetteRank'),
                                       spat_methods_params = list(NA, NA, NA, NA, NA),
                                       spat_methods_names = c('binSpect_single', 'binSpect_multi', 'spatialDE', 'spark', 'silhouetteRank'),

                                       save_plot = FALSE,
                                       save_raw = FALSE,
                                       save_norm = FALSE,
                                       save_dir = '~',
                                       verbose = TRUE,
                                       run_simulations = TRUE,
                                       ... ) {


  prob_list = list()
  for(prob_ind in 1:length(spatial_probs)) {

    prob_i = spatial_probs[prob_ind]

    if(verbose) cat('\n \n start with ', prob_i, '\n \n')

    rep_list = list()
    for(rep_i in 1:reps) {


      if(verbose) cat('\n \n repetitiion = ', rep_i, '\n \n')


      plot_name = paste0('plot_',gene_name,'_prob', prob_i, '_rep', rep_i)


      rep_res = run_spatial_sim_tests_one_rep(gobject,
                                              pattern_name = pattern_name,
                                              pattern_cell_ids = pattern_cell_ids,
                                              gene_name = gene_name,
                                              spatial_prob = prob_i,

                                              spatial_network_name = spatial_network_name,

                                              spat_methods = spat_methods,
                                              spat_methods_params = spat_methods_params,
                                              spat_methods_names = spat_methods_names,

                                              save_plot = save_plot,
                                              save_raw = save_raw,
                                              save_norm = save_norm,

                                              save_dir = save_dir,
                                              save_name = plot_name,
                                              run_simulations = run_simulations,
                                              ...)

      if(run_simulations == TRUE) {
        rep_res[, rep := rep_i]
        rep_list[[rep_i]] = rep_res
      }


    }

    if(run_simulations == TRUE) {
      rep_list_res = do.call('rbind', rep_list)
      prob_list[[prob_ind]] = rep_list_res
    }


  }

  if(run_simulations == TRUE) {
    final_gene_results = do.call('rbind', prob_list)
    return(final_gene_results)
  }


}




#' @title runPatternSimulation
#' @name runPatternSimulation
#' @description Creates a known spatial pattern for selected genes one-by-one and runs the different spatial gene detection tests
#' @param gobject giotto object
#' @param pattern_name name of spatial pattern
#' @param pattern_colors 2 color vector for the spatial pattern
#' @param pattern_cell_ids cell ids that make up the spatial pattern
#' @param gene_names selected genes
#' @param spatial_probs probabilities to test for a high expressing gene value to be part of the spatial pattern
#' @param reps number of random simulation repetitions
#' @param spatial_network_name which spatial network to use for binSpectSingle
#' @param spat_methods vector of spatial methods to test
#' @param spat_methods_params list of parameters list for each element in the vector of spatial methods to test
#' @param spat_methods_names name for each element in the vector of spatial elements to test
#' @param scalefactor library size scaling factor when re-normalizing dataset
#' @param save_plot save intermediate random simulation plots or not
#' @param save_raw save the raw expression matrix of the simulation
#' @param save_norm save the normalized expression matrix of the simulation
#' @param save_dir directory to save results to
#' @param max_col maximum number of columns for final plots
#' @param height height of final plots
#' @param width width of final plots
#' @param run_simulations run simulations (default = TRUE)
#' @param \dots additional parameters for renormalization
#' @return data.table with results
#' @export
runPatternSimulation = function(gobject,
                                pattern_name = 'pattern',
                                pattern_colors = c('in' = 'green', 'out' = 'red'),
                                pattern_cell_ids = NULL,
                                gene_names = NULL,
                                spatial_probs = c(0.5, 1),
                                reps = 2,
                                spatial_network_name = 'kNN_network',
                                spat_methods = c('binSpect_single', 'binSpect_multi', 'spatialDE', 'spark', 'silhouetteRank'),
                                spat_methods_params = list(NA, NA, NA, NA, NA),
                                spat_methods_names = c('binSpect_single', 'binSpect_multi', 'spatialDE', 'spark', 'silhouetteRank'),
                                scalefactor = 6000,
                                save_plot = T,
                                save_raw = T,
                                save_norm = T,
                                save_dir = '~',
                                max_col = 4,
                                height = 7,
                                width = 7,
                                run_simulations = TRUE,
                                ...) {


  # data.table variables
  prob = method = adj.p.value = time = NULL


  # plot pattern for first gene (the same for all)
  example_patch = simulateOneGenePatternGiottoObject(gobject,
                                                     pattern_name = pattern_name,
                                                     pattern_cell_ids = pattern_cell_ids,
                                                     gene_name = gene_names[[1]],
                                                     spatial_prob = 1,
                                                     scalefactor = scalefactor,
                                                     verbose = T)

  spatPlot2D(example_patch, cell_color = pattern_name, cell_color_code = pattern_colors,
             save_plot = save_plot, save_param = list(save_dir = save_dir, save_folder = 'original', save_name = paste0(pattern_name,'_pattern'),
                                                      base_width = 9, base_height = 7, units = 'cm'))


  all_results = list()
  for(gene_ind in 1:length(gene_names)) {

    gene = gene_names[gene_ind]

    # plot original expression
    spatGenePlot2D(gobject, expression_values = 'norm', genes = gene,
                   point_shape = 'border', point_border_stroke = 0.1,
                   show_network = F, network_color = 'lightgrey', point_size = 2.5,
                   cow_n_col = 1, show_plot = F,
                   save_plot = save_plot, save_param = list(save_dir = save_dir, save_folder = 'original', save_name = paste0(gene,'_original'),
                                                            base_width = 9, base_height = 7, units = 'cm'))


    generesults = run_spatial_sim_tests_multi(gobject,
                                              pattern_name = pattern_name,
                                              pattern_cell_ids = pattern_cell_ids,
                                              gene_name = gene,
                                              spatial_network_name = spatial_network_name,
                                              spat_methods = spat_methods,
                                              spat_methods_params = spat_methods_params,
                                              spat_methods_names = spat_methods_names,
                                              save_plot = save_plot,
                                              save_raw = save_raw,
                                              save_norm = save_norm,
                                              save_dir = save_dir,
                                              spatial_probs = spatial_probs,
                                              reps = reps,
                                              run_simulations = run_simulations,
                                              ...)

    if(run_simulations == TRUE) {
      generesults[, prob := as.factor(prob)]
      uniq_methods = sort(unique(generesults$method))
      generesults[, method := factor(method, levels = uniq_methods)]

      if(save_plot == TRUE) {

        subdir = paste0(save_dir,'/',pattern_name,'/')
        if(!file.exists(subdir)) dir.create(path = subdir, recursive = TRUE)
        # write results
        data.table::fwrite(x = generesults, file = paste0(subdir,'/',gene,'_results.txt'), sep = '\t', quote = F)

      }

      all_results[[gene_ind]] = generesults

    }

  }


  ## create combined results and visuals
  if(run_simulations == TRUE) {

    results = do.call('rbind', all_results)

    ## plot results ##

    if(save_plot == TRUE) {
      # 4 columns max
      nr_rows = max(c(round(length(gene_names)/max_col), 1))

      # p-values
      pl = ggplot2::ggplot()
      pl = pl + ggplot2::geom_boxplot(data = results, ggplot2::aes(x = method, y = adj.p.value, color = prob))
      pl = pl + ggplot2::geom_point(data = results, ggplot2::aes(x = method, y = adj.p.value, color = prob), size = 2, position = ggplot2::position_jitterdodge())
      pl = pl + ggplot2::theme_bw() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 1, hjust = 1))
      pl = pl + ggplot2::facet_wrap(~genes, nrow = nr_rows)
      pl = pl + ggplot2::geom_hline(yintercept = 0.05, color = 'red', linetype = 2)

      grDevices::pdf(file = paste0(save_dir,'/',pattern_name,'_boxplot_pvalues.pdf'), width = width, height = height)
      print(pl)
      grDevices::dev.off()



      # -log10 p-values
      pl = ggplot2::ggplot()
      pl = pl + ggplot2::geom_boxplot(data = results, ggplot2::aes(x = method, y = -log10(adj.p.value), color = prob))
      pl = pl + ggplot2::geom_point(data = results, ggplot2::aes(x = method, y = -log10(adj.p.value), color = prob), size = 2, position = ggplot2::position_jitterdodge())
      pl = pl + ggplot2::theme_bw() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 1, hjust = 1))
      pl = pl + ggplot2::facet_wrap(~genes, nrow = nr_rows)

      grDevices::pdf(file = paste0(save_dir,'/',pattern_name,'_boxplot_log10pvalues.pdf'), width = width, height = height)
      print(pl)
      grDevices::dev.off()


      # time
      pl = ggplot2::ggplot()
      pl = pl + ggplot2::geom_boxplot(data = results, ggplot2::aes(x = method, y = time, color = prob))
      pl = pl + ggplot2::geom_point(data = results, ggplot2::aes(x = method, y = time, color = prob), size = 2, position = ggplot2::position_jitterdodge())
      pl = pl + ggplot2::theme_bw() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 1, hjust = 1))

      grDevices::pdf(file = paste0(save_dir,'/',pattern_name,'_boxplot_time.pdf'), width = width, height = height)
      print(pl)
      grDevices::dev.off()
    }


    # write results
    data.table::fwrite(x = results, file = paste0(save_dir,'/',pattern_name,'_results.txt'), sep = '\t', quote = F)
    return(results)

  } else {
    return(NULL)
  }

}


