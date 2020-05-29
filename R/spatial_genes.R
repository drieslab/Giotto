

## spatial gene detection ####

#' @title spat_fish_func
#' @name spat_fish_func
#' @description performs fisher exact test
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
  test2 = Giotto:::t_giotto(Giotto:::t_giotto(spat_mat)*gene_vectorB)

  sourcevalues = test1[spat_mat == 1]
  targetvalues = test2[spat_mat == 1]

  # option 1
  test = paste0(sourcevalues,'-',targetvalues)


  if(length(unique(test)) < 4) {

    possibs = c("1-1","0-1","1-0","0-0")
    missings_possibs = possibs[!possibs %in% unique(test)]
    test = c(test, missings_possibs)
  }

  if(calc_hub == TRUE) {
    high_cells = names(gene_vector[gene_vector == 1])
    subset_spat_mat = spat_mat[rownames(spat_mat) %in% high_cells, colnames(spat_mat) %in% high_cells]

    if(length(subset_spat_mat) == 1) {
      hub_nr = 0
    } else {
      subset_spat_mat = spat_mat[rownames(spat_mat) %in% high_cells, colnames(spat_mat) %in% high_cells]
      rowhubs = Giotto:::rowSums_giotto(subset_spat_mat)
      colhubs = Giotto:::colSums_giotto(subset_spat_mat)
      hub_nr = length(unique(c(names(colhubs[colhubs > hub_min_int]), names(rowhubs[colhubs > hub_min_int]))))
    }
    fish_res = fisher.test(matrix(table(test), byrow = T, nrow = 2))[c('p.value','estimate')]
    return(c(genes = list(gene), fish_res, hubs = list(hub_nr)))

  } else {

    fish_res = fisher.test(matrix(table(test), byrow = T, nrow = 2))[c('p.value','estimate')]
    return(c(genes = list(gene), fish_res))
  }

}


#' @title spat_OR_func
#' @name spat_OR_func
#' @description calculate odds-ratio
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

    fish_matrix = matrix(table(test), byrow = T, nrow = 2)
    fish_matrix = fish_matrix/1000
    OR = ((fish_matrix[1]*fish_matrix[4]) / (fish_matrix[2]*fish_matrix[3]))

    return(c(genes = list(gene), OR, hubs = list(hub_nr)))

  }

  fish_matrix = matrix(table(test), byrow = T, nrow = 2)
  fish_matrix = fish_matrix/1000
  OR = ((fish_matrix[1]*fish_matrix[4]) / (fish_matrix[2]*fish_matrix[3]))
  return(c(genes = list(gene), OR))

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
#' @param nstart kmeans: nstart parameter
#' @param iter_max kmeans: iter.max parameter
#' @param percentage_rank percentage of top cells for binarization
#' @param do_fisher_test perform fisher test
#' @param calc_hub calculate the number of hub cells
#' @param hub_min_int minimum number of cell-cell interactions for a hub cell
#' @param get_av_expr calculate the average expression per gene of the high expressing cells
#' @param get_high_expr calculate the number of high expressing cells  per gene
#' @param do_parallel run calculations in parallel with mclapply
#' @param cores number of cores to use if do_parallel = TRUE
#' @param verbose be verbose
#' @return data.table with results (see details)
#' @details We provide two ways to identify spatial genes based on gene expression binarization.
#' Both methods are identicial except for how binarization is performed.
#' \itemize{
#'   \item{1. binarize: }{Each gene is binarized (0 or 1) in each cell with \bold{kmeans} (k = 2) or based on \bold{rank} percentile}
#'   \item{2. network: }{Alll cells are connected through a spatial network based on the physical coordinates}
#'   \item{3. contingency table: }{A contingency table is calculated based on all edges of neighboring cells and the binarized expression (0-0, 0-1, 1-0 or 1-1)}
#'   \item{4. For each gene an odds-ratio (OR) and fisher.test (optional) is calculated}
#' }
#' Other statistics are provided (optional):
#' \itemize{
#'   \item{Number of cells with high expression (binary = 1)}
#'   \item{Average expression of each gene within high expressing cells }
#'   \item{Number of hub cells, these are high expressing cells that have a user defined number of high expressing neighbors}
#' }
#' By selecting a subset of likely spatial genes (e.g. soft thresholding highly variable genes) or using multiple cores can accelerate the speed.
#' @export
#' @examples
#'     binSpect(gobject)
binSpect = function(gobject,
                    bin_method = c('kmeans', 'rank'),
                    expression_values = c('normalized', 'scaled', 'custom'),
                    subset_genes = NULL,
                    spatial_network_name = 'Delaunay_network',
                    nstart = 3,
                    iter_max = 10,
                    percentage_rank = 30,
                    do_fisher_test = TRUE,
                    calc_hub = FALSE,
                    hub_min_int = 3,
                    get_av_expr = TRUE,
                    get_high_expr = TRUE,
                    do_parallel = TRUE,
                    cores = NA,
                    verbose = T) {


  # set binarization method
  bin_method = match.arg(bin_method, choices = c('kmeans', 'rank'))

  # spatial network
  spatial_network = select_spatialNetwork(gobject,name = spatial_network_name,return_network_Obj = FALSE)
  if(is.null(spatial_network)) {
    stop('spatial_network_name: ', spatial_network_name, ' does not exist, create a spatial network first')
  }

  # expression
  values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))
  expr_values = Giotto:::select_expression_values(gobject = gobject, values = values)

  if(!is.null(subset_genes)) {
    expr_values = expr_values[rownames(expr_values) %in% subset_genes, ]
  }

  # binarize matrix
  if(bin_method == 'kmeans') {
    bin_matrix = t_giotto(apply(X = expr_values, MARGIN = 1, FUN = Giotto:::kmeans_binarize, nstart = nstart, iter.max = iter_max))
  } else if(bin_method == 'rank') {
    max_rank = (ncol(expr_values)/100)*percentage_rank
    bin_matrix = t_giotto(apply(X = expr_values, MARGIN = 1, FUN = Giotto:::rank_binarize, max_rank = max_rank))
  }

  if(verbose == TRUE) cat('\n 1. matrix binarization complete \n')


  # spatial matrix
  dc_spat_network = dcast.data.table(spatial_network, formula = to~from, value.var = 'distance', fill = 0)
  spat_mat = Giotto:::dt_to_matrix(dc_spat_network)
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

  if(verbose == TRUE) cat('\n 2. spatial enrichment test completed \n')

  result = as.data.table(do.call('rbind', save_list))
  result[, genes := unlist(genes)]
  #result[, genes := rownames(bin_matrix)]

  ## extra info: average expression of high expression group
  if(get_av_expr == TRUE) {
    sel_expr_values = expr_values * bin_matrix
    av_expr = apply(sel_expr_values, MARGIN = 1, FUN = function(x) {
      mean(x[x > 0])
    })
    av_expr_DT = data.table::data.table(genes = names(av_expr), av_expr = av_expr)
    result = merge(result, av_expr_DT, by = 'genes')

    if(verbose == TRUE) cat('\n 3. average expression of high expressing cells calculated \n')
  }

  ## extra info: number of high expressing cells
  if(get_high_expr == TRUE) {
    high_expr = rowSums(bin_matrix)
    high_expr_DT = data.table::data.table(genes = names(high_expr), high_expr = high_expr)
    result = merge(result, high_expr_DT, by = 'genes')

    if(verbose == TRUE) cat('\n 4. number of high expressing cells calculated \n')
  }


  ## order data.table
  if(do_fisher_test == TRUE) {
    result[, c('p.value', 'estimate') := list(as.numeric(p.value), as.numeric(estimate))]

    # convert p.value = 0 to lowest p-value
    min_pvalue = min(result$p.value[result$p.value > 0])
    result[, p.value := ifelse(p.value == 0, min_pvalue, p.value)]

    result[, score := -log(p.value) * estimate]
    data.table::setorder(result, -score)



  } else {
    data.table::setnames(result, 'V1', 'estimate')
    data.table::setorder(result, -estimate)
  }

  return(result)

}






#' @title silhouetteRank
#' @name silhouetteRank
#' @description Previously: calculate_spatial_genes_python. This method computes a silhouette score per gene based on the
#' spatial distribution of two partitions of cells (expressed L1, and non-expressed L0).
#' Here, rather than L2 Euclidean norm, it uses a rank-transformed, exponentially weighted
#' function to represent the local physical distance between two cells.
#' @param gobject giotto object
#' @param expression_values expression values to use
#' @param metric distance metric to use
#' @param subset_genes only run on this subset of genes
#' @param rbp_p fractional binarization threshold
#' @param examine_top top fraction to evaluate with silhouette
#' @param python_path specify specific path to python if required
#' @return data.table with spatial scores
#' @export
#' @examples
#'     silhouetteRank(gobject)
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


  # spatial locations
  spatlocs = as.matrix(gobject@spatial_locs[,.(sdimx, sdimy)])

  # python path
  if(is.null(python_path)) {
    python_path = readGiottoInstructions(gobject, param = "python_path")
  }

  ## prepare python path and louvain script
  reticulate::use_python(required = T, python = python_path)
  python_leiden_function = system.file("python", "python_spatial_genes.py", package = 'Giotto')
  reticulate::source_python(file = python_leiden_function)


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
#' @param save_param list of saving parameters from all_plots_save_function()
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return a list of data.frames with results and plot (optional)
#' @details This function is a wrapper for the SpatialDE method implemented in the ...
#' @export
#' @examples
#'     spatialDE(gobject)
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


  # expression
  values = match.arg(expression_values, c('raw', 'normalized', 'scaled', 'custom'))
  expr_values = Giotto:::select_expression_values(gobject = gobject, values = values)

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
  Spatial_DE_results = Spatial_DE(as.data.frame(t(expr_values)), spatial_locs)

  results <- as.data.frame(reticulate::py_to_r(Spatial_DE_results[[1]]))

  if(length(Spatial_DE_results) == 2){
    ms_results = as.data.frame(reticulate::py_to_r(Spatial_DE_results[[2]]))
    spatial_genes_results = list(results, ms_results)
    names(spatial_genes_results) = c("results", "ms_results")
  } else{
    spatial_genes_results =  results
    ms_results = NULL
  }


  ## create plot
  FSV_plot = FSV_show(results = results,
                      ms_results = ms_results,
                      size =size,
                      color = color,
                      sig_alpha = sig_alpha,
                      unsig_alpha = unsig_alpha)


  # print, return and save parameters
  show_plot = ifelse(is.na(show_plot), readGiottoInstructions(gobject, param = 'show_plot'), show_plot)
  save_plot = ifelse(is.na(save_plot), readGiottoInstructions(gobject, param = 'save_plot'), save_plot)
  return_plot = ifelse(is.na(return_plot), readGiottoInstructions(gobject, param = 'return_plot'), return_plot)


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
#' @param SpatialDE_results results of \code{\link{SpatialDE}} function
#' @param name_pattern name for the computed spatial patterns
#' @param expression_values gene expression values to use
#' @param pattern_num number of spatial patterns to look for
#' @param l lengthscale
#' @param python_path specify specific path to python if required
#' @param return_gobject show plot
#' @return An updated giotto object
#' @details This function is a wrapper for the SpatialAEH method implemented in the ...
#' @export
#' @examples
#'     spatialAEH(gobject)
spatialAEH <- function(gobject = NULL,
                       SpatialDE_results = NULL,
                       name_pattern = 'AEH_patterns',
                       expression_values = c('raw', 'normalized', 'scaled', 'custom'),
                       pattern_num = 6,
                       l = 1.05,
                       python_path = NULL,
                       return_gobject = TRUE) {

  # expression
  values = match.arg(expression_values, c('raw', 'normalized', 'scaled', 'custom'))
  expr_values = Giotto:::select_expression_values(gobject = gobject, values = values)

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
  AEH_results = Spatial_DE_AEH(filterd_exprs = as.data.frame(t(expr_values)),
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



#' @title trendSceek
#' @name trendSceek
#' @description Compute spatial variable genes with trendsceek method
#' @param gobject Giotto object
#' @param expression_values gene expression values to use
#' @param subset_genes subset of genes to run trendsceek on
#' @param nrand An integer specifying the number of random resamplings of the mark distribution as to create the null-distribution.
#' @param ncores An integer specifying the number of cores to be used by BiocParallel
#' @param ... Additional parameters to the \code{\link[dbscan]{trendsceek_test}} function
#' @return data.frame with trendsceek spatial genes results
#' @details This function is a wrapper for the trendsceek_test method implemented in the trendsceek package
#' @export
#' @examples
#'     trendSceek(gobject)
trendSceek <- function(gobject,
                       expression_values = c("normalized", "raw"),
                       subset_genes = NULL,
                       nrand = 100,
                       ncores = 8,
                       ...) {

  if("trendsceek" %in% rownames(installed.packages()) == FALSE) {
    stop("\n package 'trendsceek' is not yet installed \n",
         "To install: \n",
         "See https://github.com/edsgard/trendsceek"
    )
  }

  ## expression data
  values = match.arg(expression_values, c("normalized", "raw"))
  expr_values = Giotto:::select_expression_values(gobject = gobject, values = values)

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
  spatial_locations = copy(gobject@spatial_locs)
  spatial_locations[, cell_ID := NULL]
  pp = trendsceek::pos2pp(spatial_locations)

  ## initial gene counts
  pp = trendsceek::set_marks(pp, expr_values, log.fcn = log.fcn)

  # eliminates running errors caused by too many zeros
  pp[["marks"]] = pp[["marks"]] + 1e-7

  ## run trendsceek
  trendsceektest = trendsceek::trendsceek_test(pp, nrand = nrand, ncores = ncores, ...)

  ## get final results
  trendsceektest = trendsceektest$supstats_wide

  return(trendsceektest)
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
#' @examples
#'     detectSpatialPatterns(gobject)
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
  spatial_grid = gobject@spatial_grid[[spatial_grid_name]]


  # annotate spatial locations with spatial grid information
  spatial_locs = copy(gobject@spatial_locs)

  if(all(c('sdimx', 'sdimy', 'sdimz') %in% colnames(spatial_locs))) {
    spatial_locs = annotate_spatlocs_with_spatgrid_3D(spatloc = spatial_locs, spatgrid = spatial_grid)
  } else if(all(c('sdimx', 'sdimy') %in% colnames(spatial_locs))) {
    spatial_locs = annotate_spatlocs_with_spatgrid_2D(spatloc = spatial_locs, spatgrid = spatial_grid)
  }



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
#' @param show_plot show plot
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters from \code{\link{all_plots_save_function}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return ggplot
#' @export
#' @examples
#'     showPattern2D(gobject)
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
#' @param dimension dimension to plot
#' @param trim Trim ends of the PC values.
#' @param background_color background color for plot
#' @param grid_border_color color for grid
#' @param show_legend show legend of ggplot
#' @param show_plot show plot
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters from \code{\link{all_plots_save_function}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return ggplot
#' @seealso \code{\link{showPattern2D}}
#' @export
#' @examples
#'     showPattern(gobject)
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
#' @param save_param list of saving parameters from \code{\link{all_plots_save_function}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return plotly
#' @export
#' @examples
#'     showPattern3D(gobject)
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
#' @param save_param list of saving parameters from all_plots_save_function()
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return ggplot
#' @export
#' @examples
#'     showPatternGenes(gobject)
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
#' @return Data.table with genes associated with selected dimension (PC).
#' @details Description.
#' @export
#' @examples
#'     selectPatternGenes(gobject)
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


#' @title Spatial_DE
#' @name Spatial_DE
#' @description calculate spatial varible genes with spatialDE method
#' @param gobject Giotto object
#' @param show_plot show FSV plot
#' @param python_path specify specific path to python if required
#' @return a list or a dataframe of SVs
#' @details Description.
#' @export
#' @examples
#'     Spatial_DE(gobject)
Spatial_DE <- function(gobject = NULL,
                       show_plot = T,
                       size = c(4,2,1),
                       color = c("blue", "green", "red"),
                       sig_alpha = 0.5,
                       unsig_alpha = 0.5,
                       python_path = NULL){

  ## python path
  if(is.null(python_path)) {
    python_path = readGiottoInstructions(gobject, param = "python_path")
    #python_path = system('which python', intern = T)
  }

  ## source python file
  reticulate::use_python(required = T, python = python_path)
  reader_path = system.file("python", "SpatialDE_wrapper.py", package = 'Giotto')
  reticulate::source_python(file = reader_path)

  #reader_path = system.file("python", "SpatialDE_wrapper.py", package = 'Giotto')
  #source_python(reader_path)

  spatial_locs <- as.data.frame(gobject@spatial_locs)
  rownames(spatial_locs) <- spatial_locs$cell_ID
  spatial_locs <- subset(spatial_locs, select = -cell_ID)

  Spatial_DE_results = Spatial_DE(as.data.frame(t(gobject@raw_exprs)), spatial_locs)

  results <- as.data.frame(py_to_r(Spatial_DE_results[[1]]))

  if(length(Spatial_DE_results) == 2){
    ms_results <- as.data.frame(py_to_r(Spatial_DE_results[[2]]))
    spatial_genes_results <- list(results,ms_results)
    names(spatial_genes_results) <- c("results","ms_results")
  }

  else{
    spatial_genes_results =  results
    ms_results = NULL
  }


  if(show_plot == T){
    FSV_show(results = results,
             ms_results = ms_results,
             size =size,
             color = color,
             sig_alpha = sig_alpha,
             unsig_alpha = unsig_alpha)
  }

  return(spatial_genes_results)
}



#' @title Spatial_AEH
#' @name Spatial_AEH
#' @description calculate automatic expression histology with spatialDE method
#' @param gobject Giotto object
#' @param results output from spatial_DE
#' @param pattern_num the number of gene expression patterns
#' @param show_AEH show AEH plot
#' @param python_path specify specific path to python if required
#' @return a list or a dataframe of SVs
#' @details Description.
#' @export
#' @examples
#'     Spatial_AEH(gobject)
Spatial_AEH <- function(gobject = NULL,
                        results = NULL,
                        pattern_num = 5,
                        l = 1.05,
                        show_AEH = T,
                        sdimx = NULL,
                        sdimy = NULL,
                        point_size = 3,
                        point_alpha = 1,
                        low_color = "blue",
                        mid_color = "white",
                        high_color = "red",
                        midpoint = 0,
                        python_path = NULL){


  ## python path
  if(is.null(python_path)) {
    python_path = readGiottoInstructions(gobject, param = "python_path")
  }

  ## source python file
  reticulate::use_python(required = T, python = python_path)
  reader_path = system.file("python", "SpatialDE_wrapper.py", package = 'Giotto')
  reticulate::source_python(file = reader_path)



  if(is.null(sdimx)|is.null(sdimy)){
    sdimx = "sdimx"
    sdimy = "sdimy"
  }

  spatial_locs <- as.data.frame(gobject@spatial_locs)
  rownames(spatial_locs) <- spatial_locs$cell_ID
  spatial_locs <- subset(spatial_locs, select = -cell_ID)

  AEH_results = Spatial_DE_AEH(filterd_exprs = as.data.frame(t(gobject@raw_exprs)),
                               coordinates = spatial_locs,
                               results = as.data.frame(results),
                               pattern_num = pattern_num,
                               l = l)
  histology_results <- as.data.frame(py_to_r(AEH_results[[1]]))
  cell_pattern_score <- as.data.frame((py_to_r(AEH_results[[2]])))

  spatial_pattern_results <- list(histology_results,cell_pattern_score)
  names(spatial_pattern_results) <- c("histology_results","cell_pattern_score")

  if(show_AEH){
    GenePattern_show(gobject = gobject,
                     AEH_results = spatial_pattern_results,
                     sdimx = sdimx,
                     sdimy = sdimy,
                     point_size = point_size,
                     point_alpha = point_alpha,
                     low_color = low_color,
                     mid_color = mid_color,
                     high_color = high_color,
                     midpoint = 0)
  }

  return(spatial_pattern_results)
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
#' @examples
#'     do_spatial_knn_smoothing(gobject)
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
  expr_values = Giotto:::select_expression_values(gobject = gobject, values = values)

  if(!is.null(subset_genes)) {
    expr_values = expr_values[rownames(expr_values) %in% subset_genes,]
  }

  # merge spatial network with expression data
  expr_values_dt = data.table::as.data.table(expr_values); expr_values_dt[, gene_ID := rownames(expr_values)]
  expr_values_dt_m = data.table::melt.data.table(expr_values_dt, id.vars = 'gene_ID', variable.name = 'cell_ID')

  ## test ##
  spatial_network = Giotto:::convert_to_full_spatial_network(spatial_network)
  ## stop test ##

  #print(spatial_network)

  spatial_network_ext = data.table:::merge.data.table(spatial_network, expr_values_dt_m, by.x = 'target', by.y = 'cell_ID', allow.cartesian = T)

  #print(spatial_network_ext)

  # calculate mean over all k-neighbours
  # exclude 0's?
  # trimmed mean?
  spatial_network_ext_smooth = spatial_network_ext[, mean(value), by = c('source', 'gene_ID')]

  # convert back to matrix
  spatial_smooth_dc = dcast.data.table(data = spatial_network_ext_smooth, formula = gene_ID~source, value.var = 'V1')
  spatial_smooth_matrix = Giotto:::dt_to_matrix(spatial_smooth_dc)
  spatial_smooth_matrix = spatial_smooth_matrix[rownames(expr_values), colnames(expr_values)]

  # combine original and smoothed values according to smoothening b
  # create best guess for b if not given
  if(is.null(b)) {
    k = median(table(spatial_network$source))
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
#' @examples
#'     do_spatial_grid_averaging(gobject)
do_spatial_grid_averaging = function(gobject,
                                     expression_values = c('normalized', 'scaled', 'custom'),
                                     subset_genes = NULL,
                                     spatial_grid_name = 'spatial_grid',
                                     min_cells_per_grid = 4) {


  # get expression matrix
  values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))
  expr_values = Giotto:::select_expression_values(gobject = gobject, values = values)

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
  spatial_grid = gobject@spatial_grid[[spatial_grid_name]]


  # annotate spatial locations with spatial grid information
  spatial_locs = copy(gobject@spatial_locs)

  if(all(c('sdimx', 'sdimy', 'sdimz') %in% colnames(spatial_locs))) {
    spatial_locs = Giotto:::annotate_spatlocs_with_spatgrid_3D(spatloc = spatial_locs, spatgrid = spatial_grid)
  } else if(all(c('sdimx', 'sdimy') %in% colnames(spatial_locs))) {
    spatial_locs = Giotto:::annotate_spatlocs_with_spatgrid_2D(spatloc = spatial_locs, spatgrid = spatial_grid)
  }


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
#' @param b smoothing factor beteen 0 and 1 (default: automatic)
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
#' @examples
#'     detectSpatialCorGenes(gobject)
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
  expr_values = Giotto:::select_expression_values(gobject = gobject, values = values)

  if(!is.null(subset_genes)) {
    expr_values = expr_values[rownames(expr_values) %in% subset_genes,]
  }


  ## spatial averaging or smoothing
  if(method == 'grid') {

    loc_av_expr_matrix = Giotto:::do_spatial_grid_averaging(gobject = gobject,
                                                            expression_values = expression_values,
                                                            subset_genes = subset_genes,
                                                            spatial_grid_name = spatial_grid_name,
                                                            min_cells_per_grid = min_cells_per_grid)

    cor_spat_matrix = cor_giotto(t_giotto(as.matrix(loc_av_expr_matrix)), method = cor_method)
    cor_spat_matrixDT = data.table::as.data.table(cor_spat_matrix)
    cor_spat_matrixDT[, gene_ID := rownames(cor_spat_matrix)]
    cor_spat_DT = data.table::melt.data.table(data = cor_spat_matrixDT,
                                              id.vars = 'gene_ID', value.name = 'spat_cor')
  }

  if(method == 'network') {

    knn_av_expr_matrix = Giotto:::do_spatial_knn_smoothing(gobject = gobject,
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
#' @examples
#'     showSpatialCorGenes(gobject)
showSpatialCorGenes = function(spatCorObject,
                               use_clus_name = NULL,
                               selected_clusters = NULL,
                               genes = NULL,
                               min_spat_cor = 0.5,
                               min_expr_cor = NULL,
                               min_cor_diff = NULL,
                               min_rank_diff = NULL,
                               show_top_genes = NULL) {

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
    filter_DT = data.table:::merge.data.table(filter_DT, clusters_DT, by = 'gene_ID')
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
#' @examples
#'     clusterSpatialCorGenes(gobject)
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
  cor_matrix = Giotto:::dt_to_matrix(cor_DT_dc)

  # re-ordering matrix
  my_gene_order = spatCorObject[['gene_order']]
  cor_matrix = cor_matrix[my_gene_order, my_gene_order]

  # cluster
  cor_dist = as.dist(1-cor_matrix)
  cor_h = hclust(d = cor_dist, method = hclust_method)
  cor_clus = cutree(cor_h, k = k)

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
#' @param save_param list of saving parameters from \code{\link{all_plots_save_function}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @param ... additional parameters to the \code{\link[ComplexHeatmap]{Heatmap}} function from ComplexHeatmap
#' @return Heatmap generated by ComplexHeatmap
#' @export
#' @examples
#'     heatmSpatialCorGenes(gobject)
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
  cor_matrix = Giotto:::dt_to_matrix(cor_DT_dc)

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
        mycolors = Giotto:::getDistinctColors(length(uniq_clusters))
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
#' @param save_param list of saving parameters from \code{\link{all_plots_save_function}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return data.table with positive (within group) and negative (outside group) scores
#' @export
#' @examples
#'     rankSpatialCorGroups(gobject)
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
  cor_matrix = Giotto:::dt_to_matrix(cor_DT_dc)

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

  res_cor_DT = data.table('clusters' = unique(clusters_part),
                          cor_score = unlist(res_cor_list),
                          cor_neg_score = unlist(res_neg_cor_list),
                          nr_genes = unlist(nr_genes_list))
  res_cor_DT[, cor_neg_adj := 1-(cor_neg_score-min(cor_neg_score))]
  res_cor_DT[, adj_cor_score := cor_neg_adj * cor_score]
  setorder(res_cor_DT, -adj_cor_score)
  res_cor_DT[, clusters := factor(x = clusters, levels = rev(clusters))]

  # print, return and save parameters
  show_plot = ifelse(is.na(show_plot), readGiottoInstructions(gobject, param = 'show_plot'), show_plot)
  save_plot = ifelse(is.na(save_plot), readGiottoInstructions(gobject, param = 'save_plot'), save_plot)
  return_plot = ifelse(is.na(return_plot), readGiottoInstructions(gobject, param = 'return_plot'), return_plot)


  pl = ggplot()
  pl = pl + geom_point(data = res_cor_DT, aes(x = clusters, y = adj_cor_score, size = nr_genes))
  pl = pl + theme_classic()
  pl = pl + labs(x = 'clusters', y = 'pos r x (1 - (neg_r - min(neg_r)))')


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





