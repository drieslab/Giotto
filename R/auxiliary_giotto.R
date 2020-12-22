
## Giotto stat functions ####

#' @title mean_giotto
#' @description mean function that works with multiple matrix representations
#' @param x vector
#' @param \dots additional parameters
#' @return numeric
#' @export
mean_giotto = function(x, ...) {

  if(methods::is(x, 'dgCMatrix')) {
    return(Matrix::mean(x, ...)) # replace with sparseMatrixStats
  } else if(methods::is(x, 'Matrix')) {
    return(Matrix::mean(x, ...))
  } else {
    return(base::mean(x, ...))
  }
}


#' @title rowSums_giotto
#' @description rowSums function that works with multiple matrix representations
#' @param mymatrix matrix object
#' @return numeric vector
#' @export
rowSums_giotto = function(mymatrix) {

  if(methods::is(mymatrix, 'dgCMatrix')) {
    return(Matrix::rowSums(mymatrix)) # replace with sparseMatrixStats
  } else if(methods::is(mymatrix, 'Matrix')) {
    return(Matrix::rowSums(mymatrix))
  } else {
    temp_matrix = as.matrix(mymatrix)
    temp_res = matrixStats::rowSums2(temp_matrix)
    names(temp_res) = rownames(temp_matrix)
    return(temp_res)
  }
}


#' @title rowMeans_giotto
#' @description rowMeans function that works with multiple matrix representations
#' @param mymatrix matrix object
#' @return numeric vector
#' @export
rowMeans_giotto = function(mymatrix) {

  if(methods::is(mymatrix, 'dgCMatrix')) {
    return(Matrix::rowMeans(mymatrix)) # replace with sparseMatrixStats
  } else if(methods::is(mymatrix, 'Matrix')) {
    return(Matrix::rowMeans(mymatrix))
  } else {
    temp_matrix = as.matrix(mymatrix)
    temp_res = matrixStats::rowMeans2(temp_matrix)
    names(temp_res) = rownames(temp_matrix)
    return(temp_res)

  }
}

#' @title colSums_giotto
#' @description colSums function that works with multiple matrix representations
#' @param mymatrix matrix object
#' @return numeric vector
#' @export
colSums_giotto = function(mymatrix) {

  if(methods::is(mymatrix, 'dgCMatrix')) {
    return(Matrix::colSums(mymatrix)) # replace with sparseMatrixStats
  } else if(methods::is(mymatrix, 'Matrix')) {
    return(Matrix::colSums(mymatrix))
  } else {
    temp_matrix = as.matrix(mymatrix)
    temp_res = matrixStats::colSums2(temp_matrix)
    names(temp_res) = colnames(temp_matrix)
    return(temp_res)
  }
}

#' @title colMeans_giotto
#' @description colMeans function that works with multiple matrix representations
#' @param mymatrix matrix object
#' @return numeric vector
#' @export
colMeans_giotto = function(mymatrix) {

  if(methods::is(mymatrix, 'dgCMatrix')) {
    return(Matrix::colMeans(mymatrix)) # replace with sparseMatrixStats
  } else if(methods::is(mymatrix, 'Matrix')) {
    return(Matrix::colMeans(mymatrix))
  } else {
    temp_matrix = as.matrix(mymatrix)
    temp_res = matrixStats::colMeans2(temp_matrix)
    names(temp_res) = colnames(temp_matrix)
    return(temp_res)
  }
}

#' @title t_giotto
#' @description t function that works with multiple matrix representations
#' @param mymatrix matrix object
#' @return transposed matrix
#' @export
t_giotto = function(mymatrix) {

  if(methods::is(mymatrix, 'dgCMatrix')) {
    return(Matrix::t(mymatrix)) # replace with sparseMatrixStats
  } else if(methods::is(mymatrix, 'Matrix')) {
    return(Matrix::t(mymatrix))
  } else {
    mymatrix = as.matrix(mymatrix)
    mymatrix = base::t(mymatrix)
    return(mymatrix)
  }
}



#' @title cor_sparse adapted from wydr package
#' @keywords internal
cor_sparse <- function(x) {
  n = nrow(x)
  covmat = (as.matrix(Matrix::crossprod(x)) - n * Matrix::tcrossprod(Matrix::colMeans(x))) / (n - 1)
  cormat = covmat / base::tcrossprod(base::sqrt(base::diag(covmat)))
  cormat
}

#' @title cor_giotto
#' @keywords internal
cor_giotto = function(x, ...) {
  x = as.matrix(x)
  return(stats::cor(x, ...))
}


## * ####
## Giotto auxiliary functions ####

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


#' @title mean_expr_det_test
#' @keywords internal
mean_expr_det_test = function(mymatrix, detection_threshold = 1) {
  mean_expr_detected = unlist(apply(X = mymatrix, MARGIN = 1, FUN = function(x) {
    detected_x = x[x > detection_threshold]
    mean(detected_x)
  }))
}

#' @title libNorm_giotto
#' @keywords internal
libNorm_giotto <- function(mymatrix, scalefactor){
  libsizes = colSums_giotto(mymatrix)

  if(methods::is(mymatrix, 'dgCMatrix')) {
    norm_expr = Matrix::t(Matrix::t(mymatrix)/ libsizes)*scalefactor # replace with sparseMatrixStats
  } else if(methods::is(mymatrix, 'Matrix')) {
    norm_expr = Matrix::t(Matrix::t(mymatrix)/ libsizes)*scalefactor
  } else {
    norm_expr = t(t(as.matrix(mymatrix))/ libsizes)*scalefactor
  }
}

#' @title logNorm_giotto
#' @keywords internal
logNorm_giotto = function(mymatrix, base, offset) {

  if(methods::is(mymatrix, 'dgCMatrix')) {
    mymatrix@x = log(mymatrix@x + offset)/log(base) # replace with sparseMatrixStats
  } else if(methods::is(mymatrix, 'Matrix')) {
    mymatrix@x = log(mymatrix@x + offset)/log(base)
  } else {
    mymatrix = log(as.matrix(mymatrix) + offset)/log(base)
  }

  return(mymatrix)
}

#' @title pDataDT
#' @description show cell metadata
#' @param gobject giotto object
#' @return data.table with cell metadata
#' @export
#' @examples
#'
#' data(mini_giotto_single_cell) # loads existing Giotto object
#' pDataDT(mini_giotto_single_cell)
#'
pDataDT <- function(gobject) {

  if(!inherits(gobject, c('ExpressionSet', 'SCESet', 'seurat', 'giotto'))) {
    stop('only works with ExpressionSet (-like) objects')
  }

  if(inherits(gobject, c('ExpressionSet', 'SCESet'))) {
    return(data.table::as.data.table(Biobase::pData(gobject)))
  }
  else if(inherits(gobject, 'giotto')) {
    return(gobject@cell_metadata)
  }
  else if(inherits(gobject, 'seurat')) {
    return(data.table::as.data.table(gobject@meta.data))
  }

}

#' @title fDataDT
#' @description show gene metadata
#' @param gobject giotto object
#' @return data.table with gene metadata
#' @export
#' @examples
#'
#' data(mini_giotto_single_cell) # loads existing Giotto object
#' fDataDT(mini_giotto_single_cell)
#'
fDataDT <- function(gobject) {

  if(!inherits(gobject, c('ExpressionSet', 'SCESet', 'giotto'))) {
    stop('only works with ExpressionSet (-like) objects')
  }
  else if(inherits(gobject,'giotto')) {
    return(gobject@gene_metadata)
  }
  return(data.table::as.data.table(Biobase::fData(gobject)))

}


#' @title select_expression_values
#' @description helper function to select expression values
#' @param gobject giotto object
#' @param values expression values to extract
#' @return expression matrix
#' @keywords internal
select_expression_values <- function(gobject, values) {

  if(values == 'scaled' & is.null(gobject@norm_scaled_expr)) {
    stop('run first scaling step')
  } else if(values == 'scaled') {
    expr_values = gobject@norm_scaled_expr
  } else if(values == 'normalized' & is.null(gobject@norm_expr)) {
    stop('run first normalization step')
  } else if(values == 'normalized') {
    expr_values = gobject@norm_expr
  } else if(values == 'custom' & is.null(gobject@custom_expr)) {
    stop('first add custom expression matrix')
  } else if(values == 'custom') {
    expr_values = gobject@custom_expr
  } else if(values == 'raw') {
    expr_values = gobject@raw_exprs
  }

  return(expr_values)

}


#' @title create_average_DT
#' @description calculates average gene expression for a cell metadata factor (e.g. cluster)
#' @param gobject giotto object
#' @param meta_data_name name of metadata column to use
#' @param expression_values which expression values to use
#' @return data.table with average gene epression values for each factor
#' @keywords internal
create_average_DT <- function(gobject, meta_data_name,
                              expression_values = c('normalized', 'scaled', 'custom')) {


  # expression values to be used
  values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))
  expr_data = select_expression_values(gobject = gobject, values = values)

  # metadata
  cell_metadata = pDataDT(gobject)
  myrownames = rownames(expr_data)

  savelist <- list()
  for(group in unique(cell_metadata[[meta_data_name]])) {

    name = paste0('cluster_', group)

    temp = expr_data[, cell_metadata[[meta_data_name]] == group]
    temp_DT = rowMeans_giotto(temp)

    savelist[[name]] <- temp_DT
  }

  finalDF = do.call('cbind', savelist)
  rownames(finalDF) = myrownames

  return(as.data.frame(finalDF))
}

#' @title create_average_detection_DT
#' @description calculates average gene detection for a cell metadata factor (e.g. cluster)
#' @param gobject giotto object
#' @param meta_data_name name of metadata column to use
#' @param expression_values which expression values to use
#' @param detection_threshold detection threshold to consider a gene detected
#' @return data.table with average gene epression values for each factor
#' @keywords internal
create_average_detection_DT <- function(gobject, meta_data_name,
                                        expression_values = c('normalized', 'scaled', 'custom'),
                                        detection_threshold = 0) {

  # expression values to be used
  values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))
  expr_data = select_expression_values(gobject = gobject, values = values)

  # metadata
  cell_metadata <- pDataDT(gobject)
  myrownames <- rownames(expr_data)

  savelist <- list()
  for(group in unique(cell_metadata[[meta_data_name]])) {

    name = paste0('cluster_', group)

    temp = expr_data[, cell_metadata[[meta_data_name]] == group]
    temp = as.matrix(temp)

    if(is.matrix(temp)) {
      temp_DT = rowSums_giotto(temp > detection_threshold)/ncol(temp)
    } else {
      temp_DT = as.numeric(temp > detection_threshold)
    }

    savelist[[name]] <- temp_DT
  }

  finalDF = do.call('cbind', savelist)
  rownames(finalDF) = myrownames

  return(as.data.frame(finalDF))
}





#' @title subsetGiotto
#' @description subsets Giotto object including previous analyses.
#' @param gobject giotto object
#' @param cell_ids cell IDs to keep
#' @param gene_ids gene IDs to keep
#' @param verbose be verbose
#' @return giotto object
#' @export
#' @examples
#' \donttest{
#'
#'data(mini_giotto_single_cell)
#'
#'random_cells = sample(slot(mini_giotto_single_cell, 'cell_ID'), 10)
#'random_genes = sample(slot(mini_giotto_single_cell, 'gene_ID'), 10)
#'
#'subset_obj = subsetGiotto(mini_giotto_single_cell,
#'                          cell_ids = random_cells,
#'                          gene_ids = random_genes)
#'
#' }
#'
subsetGiotto <- function(gobject,
                         cell_ids = NULL,
                         gene_ids = NULL,
                         verbose = FALSE) {


  g_cell_IDs = gobject@cell_ID
  g_gene_IDs = gobject@gene_ID

  ## filter index
  if(!is.null(cell_ids)) {
    filter_bool_cells = g_cell_IDs %in% cell_ids
  } else filter_bool_cells = g_cell_IDs %in% g_cell_IDs
  if(!is.null(gene_ids)) {
    filter_bool_genes = g_gene_IDs %in% gene_ids
  } else filter_bool_genes = g_gene_IDs %in% g_gene_IDs

  cells_to_keep = g_cell_IDs[filter_bool_cells]
  genes_to_keep = g_gene_IDs[filter_bool_genes]

  ## FILTER ##
  # filter raw data
  gobject@raw_exprs = gobject@raw_exprs[filter_bool_genes, filter_bool_cells]

  # filter spatial locations
  gobject@spatial_locs = gobject@spatial_locs[filter_bool_cells]

  # filter cell_ID and gene_ID
  gobject@cell_ID = colnames(gobject@raw_exprs)
  gobject@gene_ID = rownames(gobject@raw_exprs)



  ## FILTER optional slots ##

  ## expression data ##
  # normalized expression
  if(!is.null(gobject@norm_expr)) {
    gobject@norm_expr = gobject@norm_expr[filter_bool_genes, filter_bool_cells]
  }
  # (normalized) rescaled expression
  if(!is.null(gobject@norm_scaled_expr)) {
    gobject@norm_scaled_expr = gobject@norm_scaled_expr[filter_bool_genes, filter_bool_cells]
  }
  # custom expression
  if(!is.null(gobject@custom_expr)) {
    gobject@custom_expr = gobject@custom_expr[filter_bool_genes, filter_bool_cells]
  }

  ## cell & gene metadata ##
  # cell metadata
  if(!is.null(gobject@cell_metadata)) {
    gobject@cell_metadata = gobject@cell_metadata[filter_bool_cells,]
  }
  # gene metadata
  if(!is.null(gobject@gene_metadata)) {
    gobject@gene_metadata = gobject@gene_metadata[filter_bool_genes,]
  }

  # data.table variables
  to = from = V = NULL

  ## spatial network & grid ##
  # cell spatial network
  if(!is.null(gobject@spatial_network)) {
    for(network in names(gobject@spatial_network)) {
      gobject@spatial_network[[network]]$networkDT =   gobject@spatial_network[[network]]$networkDT[to %in% cells_to_keep & from %in% cells_to_keep]
    }
  }

  # spatial grid
  # need to be recomputed


  ## dimension reduction ##
  # cell dim reduction
  if(!is.null(gobject@dimension_reduction$cells)) {
    if(verbose == TRUE) print(' subset dimensions reductions ')

    # for pca
    for(pca_name in names(gobject@dimension_reduction[['cells']][['pca']]) ) {
      old_coord = gobject@dimension_reduction[['cells']][['pca']][[pca_name]][['coordinates']]
      new_coord = old_coord[rownames(old_coord) %in% cells_to_keep,]
      gobject@dimension_reduction[['cells']][['pca']][[pca_name]][['coordinates']] = new_coord
    }

    # for umap
    for(umap_name in names(gobject@dimension_reduction[['cells']][['umap']]) ) {
      old_coord = gobject@dimension_reduction[['cells']][['umap']][[umap_name]][['coordinates']]
      new_coord = old_coord[rownames(old_coord) %in% cells_to_keep,]
      gobject@dimension_reduction[['cells']][['umap']][[umap_name]][['coordinates']] = new_coord
    }

    # for tsne
    for(tsne_name in names(gobject@dimension_reduction[['cells']][['tsne']]) ) {
      old_coord = gobject@dimension_reduction[['cells']][['tsne']][[tsne_name]][['coordinates']]
      new_coord = old_coord[rownames(old_coord) %in% cells_to_keep,]
      gobject@dimension_reduction[['cells']][['tsne']][[tsne_name]][['coordinates']] = new_coord
    }

  }


  ## nn network ##
  if(!is.null(gobject@nn_network$cells)) {
    if(verbose == TRUE) print(' subset networks ')

    for(knn_name in names(gobject@nn_network[['kNN']])) {

      # layout
      old_layout = gobject@nn_network[['cells']][['kNN']][[knn_name]][['layout']]
      if(!is.null(old_layout)) {
        new_layout = old_layout[filter_bool_cells,]
        gobject@nn_network[['cells']][['kNN']][[knn_name]][['layout']] = new_layout
      }
      # igraph object
      old_graph = gobject@nn_network[['cells']][['kNN']][[knn_name]][['igraph']]
      vertices_to_keep = V(old_graph)[filter_bool_cells]
      new_subgraph = igraph::subgraph(graph = old_graph, v = vertices_to_keep)
      gobject@nn_network[['cells']][['kNN']][[knn_name]][['igraph']] = new_subgraph
    }

    for(snn_name in names(gobject@nn_network[['sNN']])) {

      # layout
      old_layout = gobject@nn_network[['cells']][['sNN']][[snn_name]][['layout']]
      if(!is.null(old_layout)) {
        new_layout = old_layout[filter_bool_cells,]
        gobject@nn_network[['cells']][['sNN']][[snn_name]][['layout']] = new_layout
      }
      # igraph object
      old_graph = gobject@nn_network[['cells']][['sNN']][[snn_name]][['igraph']]
      vertices_to_keep = V(old_graph)[filter_bool_cells]
      new_subgraph = igraph::subgraph(graph = old_graph, v = vertices_to_keep)
      gobject@nn_network[['cells']][['sNN']][[snn_name]][['igraph']] = new_subgraph
    }

  }


  ## spatial enrichment ##
  if(!is.null(gobject@spatial_enrichment)) {
    if(verbose == TRUE) print(' subset spatial enrichment results ')
    for(spat_enrich_name in names(gobject@spatial_enrichment)) {
      gobject@spatial_enrichment[[spat_enrich_name]] = gobject@spatial_enrichment[[spat_enrich_name]][filter_bool_cells]
    }
  }


  ## update parameters used ##
  parameters_list = gobject@parameters
  number_of_rounds = length(parameters_list)
  update_name = paste0(number_of_rounds,'_subset')
  # parameters to include
  cells_removed = length(filter_bool_cells[filter_bool_cells==FALSE])
  genes_removed = length(filter_bool_genes[filter_bool_genes==FALSE])
  parameters_list[[update_name]] = c('cells removed' = cells_removed,
                                     'genes removed' = genes_removed)
  gobject@parameters = parameters_list


  return(gobject)

}





#' @title subsetGiottoLocs
#' @description subsets Giotto object based on spatial locations
#' @param gobject giotto object
#' @param x_max maximum x-coordinate
#' @param x_min minimum x-coordinate
#' @param y_max maximum y-coordinate
#' @param y_min minimum y-coordinate
#' @param z_max maximum z-coordinate
#' @param z_min minimum z-coordinate
#' @param return_gobject return Giotto object
#' @param verbose be verbose
#' @return giotto object
#' @details if return_gobject = FALSE, then a filtered combined metadata data.table will be returned
#' @export
#' @examples
#' \donttest{
#'
#' data(mini_giotto_single_cell)
#'
#' # spatial plot
#' spatPlot(mini_giotto_single_cell)
#'
#' # subset giotto object based on spatial locations
#' subset_obj = subsetGiottoLocs(mini_giotto_single_cell,
#' x_max = 1500, x_min = 1000,
#' y_max = -500, y_min = -1000)
#'
#' # spatial plot of subset giotto object
#' spatPlot(subset_obj)
#'
#' }
subsetGiottoLocs = function(gobject,
                            x_max = NULL,
                            x_min = NULL,
                            y_max = NULL,
                            y_min = NULL,
                            z_max = NULL,
                            z_min = NULL,
                            return_gobject = T,
                            verbose = FALSE) {

  comb_metadata = combineMetadata(gobject = gobject)
  comb_colnames =  colnames(comb_metadata)

  # x spatial dimension
  if('sdimx' %in% comb_colnames) {
    if(is.null(x_max)) x_max = max(comb_colnames[['sdimx']])
    if(is.null(x_min)) x_min = min(comb_colnames[['sdimx']])

    comb_metadata = comb_metadata[get('sdimx') < x_max & get('sdimx') > x_min]
  }

  # y spatial dimension
  if('sdimy' %in% comb_colnames) {
    if(is.null(y_max)) y_max = max(comb_colnames[['sdimy']])
    if(is.null(y_min)) y_min = min(comb_colnames[['sdimy']])

    comb_metadata = comb_metadata[get('sdimy') < y_max & get('sdimy') > y_min]
  }

  # z spatial dimension
  if('sdimz' %in% comb_colnames) {
    if(is.null(z_max)) z_max = max(comb_colnames[['sdimz']])
    if(is.null(z_min)) z_min = min(comb_colnames[['sdimz']])

    comb_metadata = comb_metadata[get('sdimz') < z_max & get('sdimz') > z_min]
  }

  if(return_gobject == TRUE) {

    filtered_cell_IDs = comb_metadata[['cell_ID']]

    subset_object = subsetGiotto(gobject = gobject, cell_ids = filtered_cell_IDs, verbose = verbose)

    return(subset_object)

  } else {
    return(comb_metadata)
  }

}







#' @title filterDistributions
#' @description show gene or cell distribution after filtering on expression threshold
#' @param gobject giotto object
#' @param expression_values expression values to use
#' @param expression_threshold threshold to consider a gene expressed
#' @param detection consider genes or cells
#' @param plot_type type of plot
#' @param nr_bins number of bins for histogram plot
#' @param fill_color fill color for plots
#' @param scale_axis ggplot transformation for axis (e.g. log2)
#' @param axis_offset offset to be used together with the scaling transformation
#' @param show_plot show plot
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters from \code{\link{all_plots_save_function}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return ggplot object
#' @export
#' @examples
#'
#' data(mini_giotto_single_cell)
#'
#' # distribution plot of genes
#' filterDistributions(mini_giotto_single_cell, detection = 'genes')
#'
#' # distribution plot of cells
#' filterDistributions(mini_giotto_single_cell, detection = 'cells')
#'
filterDistributions <- function(gobject,
                                expression_values = c('raw', 'normalized', 'scaled', 'custom'),
                                expression_threshold = 1,
                                detection = c('genes', 'cells'),
                                plot_type = c('histogram', 'violin'),
                                nr_bins = 30,
                                fill_color = 'lightblue',
                                scale_axis = 'identity',
                                axis_offset = 0,
                                show_plot = NA,
                                return_plot = NA,
                                save_plot = NA,
                                save_param =  list(),
                                default_save_name = 'filterDistributions') {

  # expression values to be used
  values = match.arg(expression_values, c('raw', 'normalized', 'scaled', 'custom'))
  expr_values = select_expression_values(gobject = gobject, values = values)

  # plot distribution for genes or cells
  detection = match.arg(detection, c('genes', 'cells'))

  # plot type
  plot_type = match.arg(plot_type, c('histogram', 'violin'))

  # variables
  V1 = NULL

  # for genes
  if(detection == 'genes') {

    gene_detection_levels = data.table::as.data.table(rowSums_giotto(expr_values >= expression_threshold))

    if(plot_type == 'violin') {

      pl <- ggplot2::ggplot()
      pl <- pl + ggplot2::theme_classic()
      pl <- pl + ggplot2::geom_violin(data = gene_detection_levels, ggplot2::aes(x = 'genes', y = V1+axis_offset),
                                      fill = fill_color)
      pl <- pl + ggplot2::scale_y_continuous(trans = scale_axis)
      pl <- pl + ggplot2::labs(y = 'gene detected in # of cells', x = '')

    } else if(plot_type == 'histogram') {

      pl <- ggplot2::ggplot()
      pl <- pl + ggplot2::theme_classic()
      pl <- pl + ggplot2::geom_histogram(data = gene_detection_levels, ggplot2::aes(x = V1+axis_offset),
                                         color = 'white', bins = nr_bins, fill = fill_color)
      pl <- pl + ggplot2::scale_x_continuous(trans = scale_axis)
      pl <- pl + ggplot2::labs(x = 'gene detected in # of cells')

    }

    # for cells
  } else if(detection == 'cells') {

    cell_detection_levels = data.table::as.data.table(colSums_giotto(expr_values >= expression_threshold))

    if(plot_type == 'violin') {

      pl <- ggplot2::ggplot()
      pl <- pl + ggplot2::theme_classic()
      pl <- pl + ggplot2::geom_violin(data = cell_detection_levels, ggplot2::aes(x = 'cells', y = V1+axis_offset),
                                      fill = fill_color)
      pl <- pl + ggplot2::scale_y_continuous(trans = scale_axis)
      pl <- pl + ggplot2::labs(y = 'genes detected per cell', x = '')

    } else if(plot_type == 'histogram') {

      pl <- ggplot2::ggplot()
      pl <- pl + ggplot2::theme_classic()
      pl <- pl + ggplot2::geom_histogram(data = cell_detection_levels, ggplot2::aes(x = V1+axis_offset),
                                         color = 'white', bins = nr_bins, fill = fill_color)
      pl <- pl + ggplot2::scale_x_continuous(trans = scale_axis)
      pl <- pl + ggplot2::labs(x = 'genes detected per cell')

    }
  }

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



#' @title filterCombinations
#' @description Shows how many genes and cells are lost with combinations of thresholds.
#' @param gobject giotto object
#' @param expression_values expression values to use
#' @param expression_thresholds all thresholds to consider a gene expressed
#' @param gene_det_in_min_cells minimum number of cells that should express a gene to consider that gene further
#' @param min_det_genes_per_cell minimum number of expressed genes per cell to consider that cell further
#' @param scale_x_axis ggplot transformation for x-axis (e.g. log2)
#' @param x_axis_offset x-axis offset to be used together with the scaling transformation
#' @param scale_y_axis ggplot transformation for y-axis (e.g. log2)
#' @param y_axis_offset y-axis offset to be used together with the scaling transformation
#' @param show_plot show plot
#' @param return_plot return only ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters from \code{\link{all_plots_save_function}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return list of data.table and ggplot object
#' @details Creates a scatterplot that visualizes the number of genes and cells that are
#' lost with a specific combination of a gene and cell threshold given an arbitrary cutoff
#' to call a gene expressed. This function can be used to make an informed decision at the
#' filtering step with filterGiotto.
#' @export
#' @examples
#'
#' data(mini_giotto_single_cell)
#'
#' # assess the effect of multiple filter criteria
#' filterCombinations(mini_giotto_single_cell,
#' gene_det_in_min_cells = c(2, 4, 8),
#' min_det_genes_per_cell = c(5, 10, 20))
#'
filterCombinations <- function(gobject,
                               expression_values = c('raw', 'normalized', 'scaled', 'custom'),
                               expression_thresholds = c(1, 2),
                               gene_det_in_min_cells = c(5, 50),
                               min_det_genes_per_cell = c(200, 400),
                               scale_x_axis = 'identity',
                               x_axis_offset = 0,
                               scale_y_axis = 'identity',
                               y_axis_offset = 0,
                               show_plot = TRUE,
                               return_plot = FALSE,
                               save_plot = NA,
                               save_param =  list(),
                               default_save_name = 'filterCombinations') {



  # expression values to be used
  values = match.arg(expression_values, c('raw', 'normalized', 'scaled', 'custom'))
  expr_values = select_expression_values(gobject = gobject, values = values)

  # gene and cell minimums need to have the same length
  if(length(gene_det_in_min_cells) != length(min_det_genes_per_cell)) {
    stop('\n gene_det_in_min_cells and min_det_genes_per_cell need to be the same size \n')
  }

  # compute the number of removed genes and cells
  result_list = list()
  for(thresh_i in 1:length(expression_thresholds)) {

    threshold = expression_thresholds[thresh_i]

    det_genes_res = list()
    det_cells_res = list()
    for(combn_i in 1:length(gene_det_in_min_cells)) {

      min_cells_for_gene = gene_det_in_min_cells[combn_i]
      min_genes_per_cell = min_det_genes_per_cell[combn_i]


      # first remove genes
      filter_index_genes = rowSums_giotto(expr_values >= threshold) >= min_cells_for_gene
      removed_genes = length(filter_index_genes[filter_index_genes == FALSE])
      det_cells_res[[combn_i]] = removed_genes

      # then remove cells
      filter_index_cells = colSums_giotto(expr_values[filter_index_genes, ] >= threshold) >= min_genes_per_cell
      removed_cells = length(filter_index_cells[filter_index_cells == FALSE])
      det_genes_res[[combn_i]] = removed_cells
    }

    temp_dt = data.table::data.table('threshold' = threshold,
                                     removed_genes = unlist(det_cells_res),
                                     removed_cells = unlist(det_genes_res))

    result_list[[thresh_i]] = temp_dt

  }

  result_DT = do.call('rbind', result_list)

  # data.table variables
  # gene_detected_in_min_cells = min_detected_genes_per_cell = combination = NULL

  # data.table variables
  gene_detected_in_min_cells = min_detected_genes_per_cell = combination = NULL

  result_DT[['gene_detected_in_min_cells']] = gene_det_in_min_cells
  result_DT[['min_detected_genes_per_cell']] = min_det_genes_per_cell
  result_DT[['combination']] = paste0(result_DT$gene_detected_in_min_cells,'-',result_DT$min_detected_genes_per_cell)

  result_DT = result_DT[,.(threshold,
                           gene_detected_in_min_cells,
                           min_detected_genes_per_cell,
                           combination,
                           removed_genes,
                           removed_cells)]

  maximum_x_value = max(result_DT[['removed_cells']], na.rm = T)
  maximum_y_value = max(result_DT[['removed_genes']], na.rm = T)

  pl <- ggplot2::ggplot()
  pl <- pl + ggplot2::theme_classic()
  pl <- pl + ggplot2::geom_line(data = result_DT, aes(x = removed_cells+x_axis_offset,
                                                      y = removed_genes+y_axis_offset,
                                                      group = as.factor(threshold)), linetype = 2)
  pl <- pl + ggplot2::geom_point(data = result_DT, aes(x = removed_cells+x_axis_offset,
                                                       y = removed_genes+y_axis_offset,
                                                       color = as.factor(threshold)))
  pl <- pl + scale_color_discrete(guide = guide_legend(title = 'threshold(s)'))
  pl <- pl + ggrepel::geom_text_repel(data = result_DT, aes(x = removed_cells+x_axis_offset,
                                                            y = removed_genes+y_axis_offset,
                                                            label = combination))
  pl <- pl + ggplot2::scale_x_continuous(trans = scale_x_axis, limits = c(0, maximum_x_value))
  pl <- pl + ggplot2::scale_y_continuous(trans = scale_y_axis, limits = c(0, maximum_y_value))
  pl <- pl + ggplot2::labs(x = 'number of removed cells', y = 'number of removed genes')


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
  } else {
    return(list(results = result_DT, ggplot = pl))
  }

}


#' @title filterGiotto
#' @description filter Giotto object based on expression threshold
#' @param gobject giotto object
#' @param expression_values expression values to use
#' @param expression_threshold threshold to consider a gene expressed
#' @param gene_det_in_min_cells minimum # of cells that need to express a gene
#' @param min_det_genes_per_cell minimum # of genes that need to be detected in a cell
#' @param verbose verbose
#' @return giotto object
#' @details The function \code{\link{filterCombinations}} can be used to explore the effect of different parameter values.
#' @export
#' @examples
#'
#' data(mini_giotto_single_cell)
#'
#' filtered_gobject = filterGiotto(mini_giotto_single_cell,
#'                                 gene_det_in_min_cells = 10,
#'                                 min_det_genes_per_cell = 10)
#'
#'
filterGiotto <- function(gobject,
                         expression_values = c('raw', 'normalized', 'scaled', 'custom'),
                         expression_threshold = 1,
                         gene_det_in_min_cells = 100,
                         min_det_genes_per_cell = 100,
                         verbose = F) {

  # expression values to be used
  values = match.arg(expression_values, c('raw', 'normalized', 'scaled', 'custom'))
  expr_values = select_expression_values(gobject = gobject, values = values)

  # approach:
  # 1. first remove genes that are not frequently detected
  # 2. then remove cells that do not have sufficient detected genes

  ## filter genes
  filter_index_genes = rowSums_giotto(expr_values >= expression_threshold) >= gene_det_in_min_cells
  selected_gene_ids = gobject@gene_ID[filter_index_genes]

  ## filter cells
  filter_index_cells = colSums_giotto(expr_values[filter_index_genes, ] >= expression_threshold) >= min_det_genes_per_cell
  selected_cell_ids = gobject@cell_ID[filter_index_cells]

  newGiottoObject = subsetGiotto(gobject = gobject,
                                 cell_ids = selected_cell_ids,
                                 gene_ids = selected_gene_ids)

  ## print output ##
  removed_genes = length(filter_index_genes[filter_index_genes == FALSE])
  total_genes   = length(filter_index_genes)

  removed_cells = length(filter_index_cells[filter_index_cells == FALSE])
  total_cells   = length(filter_index_cells)

  if(verbose == TRUE) {
    cat('Number of cells removed: ', removed_cells, ' out of ', total_cells, '\n')
    cat('Number of genes removed: ', removed_genes, ' out of ', total_genes, '\n')
  }


  ## update parameters used ##
  parameters_list  = newGiottoObject@parameters
  number_of_rounds = length(parameters_list)
  update_name      = paste0(number_of_rounds,'_filter')
  # parameters to include
  parameters_list[[update_name]] = c('used expression values' = values,
                                     'gene expression threshold' = expression_threshold,
                                     'minimum # of genes detected per cell' = min_det_genes_per_cell,
                                     'minimum times a gene is detected over all cells' = gene_det_in_min_cells)
  newGiottoObject@parameters = parameters_list

  return(newGiottoObject)


}




#' @title normalizeGiotto
#' @description fast normalize and/or scale expresion values of Giotto object
#' @param gobject giotto object
#' @param norm_methods normalization method to use
#' @param library_size_norm normalize cells by library size
#' @param scalefactor scale factor to use after library size normalization
#' @param log_norm transform values to log-scale
#' @param log_offset offset value to add to expression matrix, default = 1
#' @param logbase log base to use to log normalize expression values
#' @param scale_genes z-score genes over all cells
#' @param scale_cells z-score cells over all genes
#' @param scale_order order to scale genes and cells
#' @param verbose be verbose
#' @return giotto object
#' @details Currently there are two 'methods' to normalize your raw counts data.
#'
#' A. The standard method follows the standard protocol which can be adjusted using
#' the provided parameters and follows the following order: \cr
#' \itemize{
#'   \item{1. Data normalization for total library size and scaling by a custom scale-factor.}
#'   \item{2. Log transformation of data.}
#'   \item{3. Z-scoring of data by genes and/or cells.}
#' }
#' B. The normalization method as provided by the osmFISH paper is also implemented: \cr
#' \itemize{
#'   \item{1. First normalize genes, for each gene divide the counts by the total gene count and
#' multiply by the total number of genes.}
#'   \item{2. Next normalize cells, for each cell divide the normalized gene counts by the total
#' counts per cell and multiply by the total number of cells.}
#' }
#' This data will be saved in the Giotto slot for custom expression.
#' @export
#' @examples
#'
#'
#' data(mini_giotto_single_cell)
#'
#' norm_gobject = normalizeGiotto(mini_giotto_single_cell)
#'
#'
normalizeGiotto <- function(gobject,
                             norm_methods = c('standard', 'osmFISH'),
                             library_size_norm = TRUE,
                             scalefactor = 6e3,
                             log_norm = TRUE,
                             log_offset = 1,
                             logbase = 2,
                             scale_genes = T,
                             scale_cells = T,
                             scale_order = c('first_genes', 'first_cells'),
                             verbose = F) {

  raw_expr = gobject@raw_exprs
  gene_names = rownames(raw_expr)
  col_names = colnames(raw_expr)

  norm_methods = match.arg(arg = norm_methods, choices = c('standard', 'osmFISH'))

  # normalization according to standard methods
  if(norm_methods == 'standard') {

    ## 1. library size normalize
    if(library_size_norm == TRUE) {
      norm_expr = libNorm_giotto(mymatrix = raw_expr, scalefactor = scalefactor)
    } else {
      norm_expr = raw_expr
    }

    ## 2. lognormalize
    if(log_norm == TRUE) {
      norm_expr = logNorm_giotto(mymatrix = norm_expr,  base = logbase, offset = log_offset)
    } else {
      norm_expr = norm_expr
    }

    ## 3. scale
    if(scale_genes == TRUE & scale_cells == TRUE) {

      scale_order = match.arg(arg = scale_order, choices = c('first_genes', 'first_cells'))

      if(scale_order == 'first_genes') {
        if(verbose == TRUE) cat('\n first scale genes and then cells \n')
        if(!methods::is(norm_expr, class2 = 'matrix')) norm_expr = as.matrix(norm_expr)
        norm_scaled_expr = armaScaleRow(Z = norm_expr)
        norm_scaled_expr = armaScaleCol(Z = norm_scaled_expr)
      } else if(scale_order == 'first_cells') {
        if(verbose == TRUE) cat('\n first scale cells and then genes \n')
        if(!methods::is(norm_expr, class2 = 'matrix')) norm_expr = as.matrix(norm_expr)
        norm_scaled_expr = armaScaleCol(Z = norm_expr)
        norm_scaled_expr = armaScaleRow(Z = norm_scaled_expr)
      } else {
        stop('\n scale order must be given \n')
      }

    } else if(scale_genes == TRUE) {
      if(!methods::is(norm_expr, class2 = 'matrix')) norm_expr = as.matrix(norm_expr)
      norm_scaled_expr = armaScaleRow(Z = norm_expr)
    } else if(scale_cells == TRUE) {
      if(!methods::is(norm_expr, class2 = 'matrix')) norm_expr = as.matrix(norm_expr)
      norm_scaled_expr = armaScaleCol(Z = norm_expr)
    } else {
      norm_scaled_expr = NULL
    }


    ## 4. add cell and gene names back
    if(!is.null(norm_expr)) {
      rownames(norm_expr) = gene_names
      colnames(norm_expr) = col_names
    }
    if(!is.null(norm_scaled_expr)) {
      rownames(norm_scaled_expr) = gene_names
      colnames(norm_scaled_expr) = col_names
    }

    # return Giotto object
    gobject@norm_expr = norm_expr
    gobject@norm_scaled_expr = norm_scaled_expr

  }

  # normalization according to osmFISH method
  else if(norm_methods == 'osmFISH') {

    # 1. normalize per gene with scale-factor equal to number of genes
    norm_genes = (raw_expr/rowSums_giotto(raw_expr)) * nrow(raw_expr)
    # 2. normalize per cells with scale-factor equal to number of cells
    norm_genes_cells = t((t(norm_genes)/colSums_giotto(norm_genes)) * ncol(raw_expr))

    # return results to Giotto object
    cat('\n osmFISH-like normalized data will be returned to the custom Giotto slot \n')
    gobject@custom_expr = norm_genes_cells

  }




  ## update parameters used ##
  parameters_list  = gobject@parameters
  number_of_rounds = length(parameters_list)
  update_name      = paste0(number_of_rounds,'_normalize')

  # parameters to include
  if(norm_methods == 'standard') {
    parameters_list[[update_name]] = c('normalization method' = norm_methods,
                                       'normalized to library size' = ifelse(library_size_norm == T, 'yes', 'no'),
                                       'scalefactor' = scalefactor,
                                       'log-normalized' =  ifelse(log_norm == T, 'yes', 'no'),
                                       'logbase' = ifelse(is.null(logbase), NA, logbase),
                                       'log offset' = log_offset,
                                       'genes scaled' = ifelse(scale_genes == T, 'yes', 'no'),
                                       'cell scaled' = ifelse(scale_cells == T, 'yes', 'no'),
                                       'if both, order of scaling' = scale_order)
  }

  if(norm_methods == 'osmFISH') {
    parameters_list[[update_name]] = c('normalization method' = norm_methods)
  }

  gobject@parameters = parameters_list

  return(gobject)
}



#' @title adjustGiottoMatrix
#' @description Adjust expression values to account for known batch effects or technological covariates.
#' @param gobject giotto object
#' @param expression_values expression values to use
#' @param batch_columns metadata columns that represent different batch (max = 2)
#' @param covariate_columns metadata columns that represent covariates to regress out
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @param update_slot expression slot that will be updated (default = custom)
#' @return giotto object
#' @details This function implements the \code{\link[limma]{removeBatchEffect}} function to
#' remove known batch effects and to adjust expression values according to provided covariates.
#' @export
#' @examples
#'
#' data(mini_giotto_single_cell)
#'
#' adjust_gobject = adjustGiottoMatrix(mini_giotto_single_cell)
#'
adjustGiottoMatrix <- function(gobject,
                               expression_values = c('normalized', 'scaled', 'custom'),
                               batch_columns = NULL,
                               covariate_columns = NULL,
                               return_gobject = TRUE,
                               update_slot = c('custom')) {

  # metadata
  cell_metadata = pDataDT(gobject)

  if(!is.null(batch_columns)) {
    if(!all(batch_columns %in% colnames(cell_metadata))) {
      stop('\n batch column name(s) were not found in the cell metadata \n')
    }
  }

  if(!is.null(covariate_columns)) {
    if(!all(covariate_columns %in% colnames(cell_metadata))) {
      stop('\n covariate column name(s) were not found in the cell metadata \n')
    }
  }

  update_slot = match.arg(update_slot, c('normalized', 'scaled', 'custom'))

  values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))
  expr_data = select_expression_values(gobject = gobject, values = values)

  # batch columns
  if(!is.null(batch_columns)) {
    batch_column_1 = cell_metadata[[ batch_columns[1] ]]
    if(length(batch_columns) > 1) {
      batch_column_2 = cell_metadata[[ batch_columns[2] ]]
    } else {
      batch_column_2 = NULL
    }
  } else {
    batch_column_1 = NULL
    batch_column_2 = NULL
  }

  # covariate columns
  if(!is.null(covariate_columns)) {
    covariates = as.matrix(cell_metadata[, covariate_columns, with = F])
  } else {
    covariates = NULL
  }


  adjusted_matrix = limma::removeBatchEffect(x = expr_data,
                                             batch = batch_column_1,
                                             batch2 =  batch_column_2,
                                             covariates = covariates)

  if(return_gobject == TRUE) {
    if(update_slot == 'normalized') {
      gobject@norm_expr = adjusted_matrix
    } else if(update_slot == 'scaled') {
      gobject@norm_scaled_expr = adjusted_matrix
    } else if(update_slot == 'custom') {
      gobject@custom_expr = adjusted_matrix
    }

    return(gobject)

  } else {

    return(adjusted_matrix)

  }

}



#' @title processGiotto
#' @description Wrapper for the different Giotto object processing functions
#' @param gobject giotto object
#' @param filter_params additional parameters to filterGiotto
#' @param norm_params additional parameters to normalizeGiotto
#' @param stat_params additional parameters to addStatistics
#' @param adjust_params additional parameters to adjustGiottoMatrix
#' @param verbose be verbose (default is TRUE)
#' @return giotto object
#' @details See \code{\link{filterGiotto}}, \code{\link{normalizeGiotto}},
#' \code{\link{addStatistics}} and \code{\link{adjustGiottoMatrix}} for more
#' information about the different parameters in each step. If you do not provide
#' them it will use the default values.
#' @export
#' @examples
#'
#' data(mini_giotto_single_cell)
#'
#' processed_object = processGiotto(mini_giotto_single_cell,
#'                                  filter_params = list(gene_det_in_min_cells = 10,
#'                                  min_det_genes_per_cell = 10))
#'
processGiotto = function(gobject,
                         filter_params = list(),
                         norm_params = list(),
                         stat_params = list(),
                         adjust_params = list(),
                         verbose = TRUE){

  # filter Giotto
  if(verbose == TRUE) cat('1. start filter step \n')
  if(class(filter_params) != 'list') stop('filter_params need to be a list of parameters for filterGiotto \n')
  gobject = do.call('filterGiotto', c(gobject = gobject, filter_params))

  # normalize Giotto
  if(verbose == TRUE) cat('2. start normalization step \n')
  if(class(norm_params) != 'list') stop('norm_params need to be a list of parameters for normalizeGiotto \n')
  gobject = do.call('normalizeGiotto', c(gobject = gobject, norm_params))

  # add Statistics
  if(verbose == TRUE) cat('3. start cell and gene statistics step \n')
  if(class(stat_params) != 'list') stop('stat_params need to be a list of parameters for addStatistics \n')
  stat_params[['return_gobject']] = TRUE # force this to be true
  gobject = do.call('addStatistics', c(gobject = gobject, stat_params))

  # adjust Giotto
  if(verbose == TRUE) cat('3. start adjusted matrix step \n')
  if(class(adjust_params) != 'list') stop('adjust_params need to be a list of parameters for adjustGiottoMatrix \n')
  adjust_params[['return_gobject']] = TRUE # force this to be true
  gobject = do.call('adjustGiottoMatrix', c(gobject = gobject, adjust_params))


  return(gobject)

}




## * ####
## Gene & Cell metadata functions ####


#' @title annotateGiotto
#' @description Converts cluster results into a user provided annotation.
#' @param gobject giotto object
#' @param annotation_vector named annotation vector (names = cluster ids)
#' @param cluster_column cluster column to convert to annotation names
#' @param name new name for annotation column
#' @return giotto object
#' @details You need to specifify which (cluster) column you want to annotate
#' and you need to provide an annotation vector like this:
#' \itemize{
#'   \item{1. identify the cell type of each cluster}
#'   \item{2. create a vector of these cell types, e.g. cell_types =  c('T-cell', 'B-cell', 'Stromal')}
#'   \item{3. provide original cluster names to previous vector, e.g. names(cell_types) = c(2, 1, 3)}
#' }
#' @export
#' @examples
#'
#' data(mini_giotto_single_cell)
#'
#' # show leiden clustering results
#' cell_metadata = pDataDT(mini_giotto_single_cell)
#' cell_metadata[['leiden_clus']]
#'
#' # create vector with cell type names as names of the vector
#' clusters_cell_types = c('cell_type_1', 'cell_type_2', 'cell_type_3')
#' names(clusters_cell_types) = 1:3
#'
#' # convert cluster results into annotations and add to cell metadata
#' mini_giotto_single_cell = annotateGiotto(gobject = mini_giotto_single_cell,
#'                                          annotation_vector = clusters_cell_types,
#'                                          cluster_column = 'leiden_clus', name = 'cell_types2')
#'
#' # visualize annotation results
#' spatDimPlot(gobject = mini_giotto_single_cell,
#'             cell_color = 'cell_types2',
#'             spat_point_size = 3, dim_point_size = 3)
#'
#'
annotateGiotto <- function(gobject,
                           annotation_vector = NULL,
                           cluster_column = NULL,
                           name = 'cell_types') {


  # data.table: set global variable
  temp_cluster_name = NULL

  if(is.null(annotation_vector) | is.null(cluster_column)) {
    stop('\n You need to provide both a named annotation vector and the corresponding cluster column  \n')
  }

  cell_metadata = pDataDT(gobject)

  # 1. verify if cluster column exist
  if(!cluster_column %in% colnames(cell_metadata)) {
    stop('\n Cluster column is not found in cell metadata \n')
  }

  # 2. verify if each cluster has an annotation
  uniq_names = names(annotation_vector)
  uniq_clusters = unique(cell_metadata[[cluster_column]])
  missing_annotations = uniq_clusters[!uniq_clusters %in% uniq_names]
  no_matching_annotations = uniq_names[!uniq_names %in% uniq_clusters]

  if(length(missing_annotations) > 0) {
    cat('Not all clusters have an accompanying annotation in the annotation_vector: \n',
        'These names are missing: ', as.character(missing_annotations), '\n',
        'These annotations have no match: ', as.character(no_matching_annotations), '\n')
    stop('Annotation interrupted \n')
  }


  # 3. remove previous annotation name if it's the same
  # but only if new name is not the same as cluster to be used
  if(name %in% colnames(cell_metadata)) {
    cat('\n annotation name ', name,' was already used \n',
        'and will be overwritten \n')

    cell_metadata[, temp_cluster_name := annotation_vector[[as.character(get(cluster_column))]], by = 1:nrow(cell_metadata)]
    cell_metadata[, (name) := NULL]

  } else {

    cell_metadata[, temp_cluster_name := annotation_vector[[as.character(get(cluster_column))]], by = 1:nrow(cell_metadata)]
  }

  data.table::setnames(cell_metadata, old = 'temp_cluster_name', new = name)
  gobject@cell_metadata = cell_metadata

  return(gobject)


}


#' @title removeCellAnnotation
#' @description removes cell annotation of giotto object
#' @param gobject giotto object
#' @param columns names of columns to remove
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @return giotto object
#' @details if return_gobject = FALSE, it will return the cell metadata
#' @export
#' @examples
#'
#' data(mini_giotto_single_cell) # load full mini giotto object
#'
#' # show cell metadata
#' pDataDT(mini_giotto_single_cell)
#'
#' # remove cell_types column
#' mini_giotto_single_cell = removeCellAnnotation(mini_giotto_single_cell,
#'                                                columns = 'cell_types')
#'
removeCellAnnotation <- function(gobject,
                                 columns = NULL,
                                 return_gobject = TRUE) {

  if(is.null(columns)) {
    stop('\t You need to provide a vector of metadata column names to remove \t')
  }

  gobject@cell_metadata[, (columns) := NULL]

  if(return_gobject == TRUE) {
    return(gobject)
  } else {
    gobject@cell_metadata
  }

}


#' @title removeGeneAnnotation
#' @description removes gene annotation of giotto object
#' @param gobject giotto object
#' @param columns names of columns to remove
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @return giotto object
#' @details if return_gobject = FALSE, it will return the gene metadata
#' @export
#' @examples
#'
#' data(mini_giotto_single_cell) # load full mini giotto object
#'
#' # show gene metadata
#' fDataDT(mini_giotto_single_cell)
#'
#' # remove nr_cells column
#' mini_giotto_single_cell = removeGeneAnnotation(mini_giotto_single_cell,
#'                                                columns = 'nr_cells')
#'
removeGeneAnnotation <- function(gobject,
                                 columns = NULL,
                                 return_gobject = TRUE) {

  if(is.null(columns)) {
    stop('\t You need to provide a vector of metadata column names to remove \t')
  }

  gobject@gene_metadata[, (columns) := NULL]

  if(return_gobject == TRUE) {
    return(gobject)
  } else {
    gobject@gene_metadata
  }

}


#' @title addCellMetadata
#' @description adds cell metadata to the giotto object
#' @param gobject giotto object
#' @param new_metadata new cell metadata to use (data.table, data.frame, ...)
#' @param vector_name (optional) custom name if you provide a single vector
#' @param by_column merge metadata based on cell_ID column in pDataDT (default = FALSE)
#' @param column_cell_ID column name of new metadata to use if by_column = TRUE
#' @return giotto object
#' @details You can add additional cell metadata in two manners:
#' \itemize{
#'   \item{1. Provide a data.table or data.frame with cell annotations in the same order as the cell_ID column in pDataDT(gobject) }
#'   \item{2. Provide a data.table or data.frame with cell annotations and specificy which column contains the cell IDs, these cell IDs need to match with the cell_ID column in pDataDT(gobject)}
#' }
#' @export
addCellMetadata <- function(gobject,
                            new_metadata,
                            vector_name = NULL,
                            by_column = FALSE,
                            column_cell_ID = NULL) {

  # data.table variables
  cell_ID = NULL

  cell_metadata = gobject@cell_metadata
  ordered_cell_IDs = gobject@cell_ID

  if(is.vector(new_metadata) | is.factor(new_metadata)) {
    original_name = deparse(substitute(new_metadata))
    new_metadata = data.table::as.data.table(new_metadata)

    if(!is.null(vector_name) & is.character(vector_name)) {
      colnames(new_metadata) = vector_name
    } else {
      colnames(new_metadata) = original_name
    }

  } else {
    new_metadata = data.table::as.data.table(new_metadata)
  }

  if(is.null(column_cell_ID)) {
    column_cell_ID = 'cell_ID'
  }

  # overwrite columns with same name
  new_col_names = colnames(new_metadata)
  new_col_names = new_col_names[new_col_names != column_cell_ID]
  old_col_names = colnames(cell_metadata)
  old_col_names = old_col_names[old_col_names != 'cell_ID']
  same_col_names = new_col_names[new_col_names %in% old_col_names]


  if(length(same_col_names) >= 1) {
    cat('\n these column names were already used: ', same_col_names, '\n',
        'and will be overwritten \n')
    cell_metadata[, (same_col_names) := NULL]
  }



  if(by_column == FALSE) {
    cell_metadata = cbind(cell_metadata, new_metadata)
  } else {
    if(is.null(column_cell_ID)) stop('You need to provide cell_ID column')
    cell_metadata <- data.table::merge.data.table(cell_metadata, by.x = 'cell_ID',
                                                   new_metadata, by.y = column_cell_ID,
                                                   all.x = T)
  }

  # reorder
  cell_metadata = cell_metadata[match(ordered_cell_IDs, cell_ID)]

  gobject@cell_metadata <- cell_metadata
  return(gobject)
}


#' @title addGeneMetadata
#' @description adds gene metadata to the giotto object
#' @param gobject giotto object
#' @param new_metadata new metadata to use
#' @param by_column merge metadata based on gene_ID column in fDataDT
#' @param column_gene_ID column name of new metadata to use if by_column = TRUE
#' @return giotto object
#' @details You can add additional gene metadata in two manners:
#' 1. Provide a data.table or data.frame with gene annotations in the same order as the gene_ID column in fDataDT(gobject)
#' 2. Provide a data.table or data.frame with gene annotations and specificy which column contains the gene IDs,
#' these gene IDs need to match with the gene_ID column in fDataDT(gobject)
#' @export
addGeneMetadata <- function(gobject,
                            new_metadata,
                            by_column = F,
                            column_gene_ID = NULL) {

  # data.table variables
  gene_ID = NULL

  gene_metadata = gobject@gene_metadata
  ordered_gene_IDs = gobject@gene_ID

  if(by_column == FALSE) {
    gene_metadata = cbind(gene_metadata, new_metadata)
  } else {
    if(is.null(column_gene_ID)) stop('You need to provide gene_ID column')
    gene_metadata <- data.table::merge.data.table(gene_metadata, by.x = 'gene_ID',
                                                   new_metadata, by.y = column_gene_ID,
                                                   all.x = T)
  }

  # reorder
  gene_metadata = gene_metadata[match(ordered_gene_IDs, gene_ID)]

  gobject@gene_metadata <- gene_metadata
  return(gobject)
}



#' @title addGeneStatistics
#' @description adds gene statistics to the giotto object
#' @param gobject giotto object
#' @param expression_values expression values to use
#' @param detection_threshold detection threshold to consider a gene detected
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @return giotto object if return_gobject = TRUE
#' @details
#' This function will add the following statistics to gene metadata:
#' \itemize{
#'   \item{nr_cells: }{Denotes in how many cells the gene is detected}
#'   \item{per_cells: }{Denotes in what percentage of cells the gene is detected}
#'   \item{total_expr: }{Shows the total sum of gene expression in all cells}
#'   \item{mean_expr: }{Average gene expression in all cells}
#'   \item{mean_expr_det: }{Average gene expression in cells with detectable levels of the gene}
#' }
#' @export
#' @examples
#'
#' data(mini_giotto_single_cell)
#'
#' updated_giotto_object = addGeneStatistics(mini_giotto_single_cell)
#'
addGeneStatistics <- function(gobject,
                              expression_values = c('normalized', 'scaled', 'custom'),
                              detection_threshold = 0,
                              return_gobject = TRUE) {

  # expression values to be used
  expression_values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))
  expr_data = select_expression_values(gobject = gobject, values = expression_values)

  # calculate stats
  gene_stats = data.table::data.table(genes = rownames(expr_data),
                          nr_cells = rowSums_giotto(expr_data > detection_threshold),
                          perc_cells = (rowSums_giotto(expr_data > detection_threshold)/ncol(expr_data))*100,
                          total_expr = rowSums_giotto(expr_data),
                          mean_expr = rowMeans_giotto(expr_data))

  # data.table variables
  mean_expr_det = NULL

  mean_expr_detected = mean_expr_det_test(expr_data, detection_threshold = detection_threshold)
  gene_stats[, mean_expr_det := mean_expr_detected]


  if(return_gobject == TRUE) {

    # remove previous statistics
    gene_metadata = fDataDT(gobject)
    metadata_names = colnames(gene_metadata)
    if('nr_cells' %in% metadata_names) {
      cat('\n gene statistics has already been applied once, will be overwritten \n')
      gene_metadata[, c('nr_cells', 'perc_cells', 'total_expr', 'mean_expr', 'mean_expr_det') := NULL]
      gobject@gene_metadata = gene_metadata
    }

    gobject = addGeneMetadata(gobject = gobject, new_metadata = gene_stats,
                                      by_column = T, column_gene_ID = 'genes')

    ## update parameters used ##
    parameters_list = gobject@parameters
    number_of_rounds = length(parameters_list)
    update_name = paste0(number_of_rounds,'_gene_stats')
    # parameters to include
    parameters_list[[update_name]] = c('expression values used' = expression_values,
                                       'detection_threshold' = detection_threshold)
    gobject@parameters = parameters_list

    return(gobject)


  } else {
    return(gene_stats)
  }

}


#' @title addCellStatistics
#' @description adds cells statistics to the giotto object
#' @param gobject giotto object
#' @param expression_values expression values to use
#' @param detection_threshold detection threshold to consider a gene detected
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @return giotto object if return_gobject = TRUE
#' @details
#' This function will add the following statistics to cell metadata:
#' \itemize{
#'   \item{nr_genes: }{Denotes in how many genes are detected per cell}
#'   \item{perc_genes: }{Denotes what percentage of genes is detected per cell}
#'   \item{total_expr: }{Shows the total sum of gene expression per cell}
#' }
#' @export
#' @examples
#'
#' data(mini_giotto_single_cell)
#'
#' updated_giotto_object = addCellStatistics(mini_giotto_single_cell)
#'
addCellStatistics <- function(gobject,
                              expression_values = c('normalized', 'scaled', 'custom'),
                              detection_threshold = 0,
                              return_gobject = TRUE) {

  # expression values to be used
  expression_values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))
  expr_data = select_expression_values(gobject = gobject, values = expression_values)

  # calculate stats
  cell_stats = data.table::data.table(cells = colnames(expr_data),
                          nr_genes = colSums_giotto(expr_data > detection_threshold),
                          perc_genes = (colSums_giotto(expr_data > detection_threshold)/nrow(expr_data))*100,
                          total_expr = colSums_giotto(expr_data))



  if(return_gobject == TRUE) {

    # remove previous statistics
    cell_metadata = pDataDT(gobject)
    metadata_names = colnames(cell_metadata)
    if('nr_genes' %in% metadata_names) {
      cat('\n cells statistics has already been applied once, will be overwritten \n')
      cell_metadata[, c('nr_genes', 'perc_genes', 'total_expr') := NULL]
      gobject@cell_metadata = cell_metadata
    }

    gobject = addCellMetadata(gobject = gobject, new_metadata = cell_stats,
                                      by_column = T, column_cell_ID = 'cells')

    ## update parameters used ##
    parameters_list = gobject@parameters
    number_of_rounds = length(parameters_list)
    update_name = paste0(number_of_rounds,'_cell_stats')
    # parameters to include
    parameters_list[[update_name]] = c('expression values used' = expression_values,
                                       'detection_threshold' = detection_threshold)
    gobject@parameters = parameters_list

    return(gobject)


  } else {
    return(cell_stats)
  }

}



#' @title addStatistics
#' @description adds genes and cells statistics to the giotto object
#' @param gobject giotto object
#' @param expression_values expression values to use
#' @param detection_threshold detection threshold to consider a gene detected
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @return giotto object if return_gobject = TRUE, else a list with results
#' @details See \code{\link{addGeneStatistics}} and \code{\link{addCellStatistics}}
#' @export
#' @examples
#'
#' data(mini_giotto_single_cell)
#'
#' updated_giotto_object = addStatistics(mini_giotto_single_cell)
#'
addStatistics <- function(gobject,
                          expression_values = c('normalized', 'scaled', 'custom'),
                          detection_threshold = 0,
                          return_gobject = TRUE) {

  # get gene statistics
  gene_stats = addGeneStatistics(gobject = gobject,
                                 expression_values = expression_values,
                                 detection_threshold = detection_threshold,
                                 return_gobject = return_gobject)

  if(return_gobject == TRUE) {
    gobject = gene_stats
  }

  # get cell statistics
  cell_stats = addCellStatistics(gobject = gobject,
                                 expression_values = expression_values,
                                 detection_threshold = detection_threshold,
                                 return_gobject = return_gobject)

  if(return_gobject == TRUE) {
    gobject = cell_stats
    return(gobject)
  } else {
    return(gene_stats = gene_stats, cell_stats = cell_stats)
  }

}


#' @title addGenesPerc
#' @description calculates the total percentage of (normalized) counts for a subset of selected genes
#' @param gobject giotto object
#' @param expression_values expression values to use
#' @param genes vector of selected genes
#' @param vector_name column name as seen in pDataDT()
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @return giotto object if return_gobject = TRUE, else a vector with % results
#' @export
#' @examples
#'
#' data(mini_giotto_single_cell)
#'
#' # select genes (e.g. Rpl or mitochondrial)
#' random_genes = sample(slot(mini_giotto_single_cell, 'gene_ID'), 5)
#'
#' # calculate percentage of those selected genes per cells/spot
#' updated_giotto_object = addGenesPerc(mini_giotto_single_cell,
#'                                      genes = random_genes,
#'                                      vector_name = 'random_gene_perc')
#'
#' # visualize result in data.table format
#' pDataDT(updated_giotto_object)
#'
#'
addGenesPerc = function(gobject,
                        expression_values = c('normalized', 'scaled', 'custom'),
                        genes = NULL,
                        vector_name = 'gene_perc',
                        return_gobject = TRUE) {

  # tests
  if(is.null(genes)) {
    stop('You need to provide a vector of gene names \n')
  }

  if(!methods::is(gobject, 'giotto')) {
    stop('You need to provide a giotto object \n')
  }


  # expression values to be used
  expression_values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))
  expr_data = select_expression_values(gobject = gobject, values = expression_values)

  totalsum = colSums_giotto(expr_data)
  gene_sum = colSums_giotto(expr_data[rownames(expr_data) %in% genes,])
  perc_genes = round((gene_sum/totalsum)*100, 2)

  if(return_gobject == TRUE) {
    temp_gobj = addCellMetadata(gobject = gobject,
                                new_metadata = perc_genes,
                                vector_name = vector_name,
                                by_column = F)
    return(temp_gobj)
  } else {
    return(perc_genes)
  }

}





## * ####
## Giotto auxiliary functions ####


#' @title showProcessingSteps
#' @description shows the sequential processing steps that were performed
#' on a Giotto object in a summarized format
#' @param gobject giotto object
#' @return list of processing steps and names
#' @export
#' @examples
#'
#' data(mini_giotto_single_cell)
#'
#' showProcessingSteps(mini_giotto_single_cell)
#'
showProcessingSteps <- function(gobject) {

  parameters = gobject@parameters

  cat('Processing steps: \n \n')

  for(step in names(parameters)) {
    cat('\n', step, '\n')

    sub_step = parameters[[step]]

    if(any(grepl('name', names(sub_step)) == TRUE)) {

      selected_names = grep('name', names(sub_step), value = T)
      cat('\t name info: ', sub_step[selected_names], '\n')

    }

  }


}



#' @title create_cluster_matrix
#' @description creates aggregated matrix for a given clustering column
#' @keywords internal
create_cluster_matrix <- function(gobject,
                                  expression_values = c('normalized', 'scaled', 'custom'),
                                  cluster_column,
                                  gene_subset = NULL) {

  # data.table variables
  genes = NULL

  values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))

  # average expression per cluster
  aggr_sc_clusters <- create_average_DT(gobject = gobject,
                                        meta_data_name = cluster_column,
                                        expression_values = values)
  aggr_sc_clusters_DT <- data.table::as.data.table(aggr_sc_clusters)
  aggr_sc_clusters_DT[, genes := rownames(aggr_sc_clusters)]
  aggr_sc_clusters_DT_melt <- data.table::melt.data.table(aggr_sc_clusters_DT,
                                                          variable.name = 'cluster',
                                                          id.vars = 'genes',
                                                          value.name = 'expression')

  # create matrix
  testmat = data.table::dcast.data.table(aggr_sc_clusters_DT_melt,
                                         formula = genes~cluster,
                                         value.var = 'expression')
  testmatrix = dt_to_matrix(testmat)
  #testmatrix = as.matrix(testmat[,-1])
  #rownames(testmatrix) = testmat[['genes']]

  # create subset if required
  if(!is.null(gene_subset)) {
    gene_subset_detected = gene_subset[gene_subset %in% rownames(testmatrix)]
    testmatrix = testmatrix[rownames(testmatrix) %in% gene_subset_detected, ]
  }

  return(testmatrix)

}





#' @title calculateMetaTable
#' @description calculates the average gene expression for one or more (combined) annotation columns.
#' @param gobject giotto object
#' @param expression_values expression values to use
#' @param metadata_cols annotation columns found in pDataDT(gobject)
#' @param selected_genes subset of genes to use
#' @return data.table with average expression values for each gene per (combined) annotation
#' @export
#' @examples
#'
#' data(mini_giotto_single_cell)
#'
#' # show cell metadata
#' pDataDT(mini_giotto_single_cell)
#'
#' # show average gene expression per annotated cell type
#' calculateMetaTable(mini_giotto_single_cell,
#'                    metadata_cols = 'cell_types')
#'
calculateMetaTable = function(gobject,
                              expression_values =  c("normalized", "scaled", "custom"),
                              metadata_cols = NULL,
                              selected_genes = NULL) {

  if(is.null(metadata_cols)) stop('\n You need to select one or more valid column names from pDataDT() \n')

  # data.table variables
  uniq_ID = NULL

  ## get metadata and create unique groups
  metadata = data.table::copy(pDataDT(gobject))
  if(length(metadata_cols) > 1) {
    metadata[, uniq_ID := paste(.SD, collapse = '-'), by = 1:nrow(metadata), .SDcols = metadata_cols]
  } else {
    metadata[, uniq_ID := get(metadata_cols)]
  }

  ## possible groups
  possible_groups = unique(metadata[,metadata_cols, with = F])
  if(length(metadata_cols) > 1) {
    possible_groups[, uniq_ID := paste(.SD, collapse = '-'), by = 1:nrow(possible_groups), .SDcols = metadata_cols]
  } else {
    possible_groups[, uniq_ID := get(metadata_cols)]
  }

  ## get expression data
  values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))
  expr_values = select_expression_values(gobject = gobject, values = values)
  if(!is.null(selected_genes)) {
    expr_values = expr_values[rownames(expr_values) %in% selected_genes, ]
  }

  ## summarize unique groups (average)
  result_list = list()

  for(row in 1:nrow(possible_groups)) {

    uniq_identifiier = possible_groups[row][['uniq_ID']]
    selected_cell_IDs = metadata[uniq_ID == uniq_identifiier][['cell_ID']]
    sub_expr_values = expr_values[, colnames(expr_values) %in% selected_cell_IDs]
    if(is.vector(sub_expr_values) == FALSE) {
      subvec = rowMeans_giotto(sub_expr_values)
    } else {
      subvec = sub_expr_values
    }
    result_list[[row]] = subvec
  }
  finaldt = data.table::as.data.table(do.call('rbind', result_list))
  possible_groups_res = cbind(possible_groups, finaldt)
  possible_groups_res_melt = data.table::melt.data.table(possible_groups_res, id.vars = c(metadata_cols, 'uniq_ID'))

  return(possible_groups_res_melt)

}


#' @title calculateMetaTableCells
#' @description calculates the average metadata values for one or more (combined) annotation columns.
#' @param gobject giotto object
#' @param value_cols metadata or enrichment value columns to use
#' @param metadata_cols annotation columns found in pDataDT(gobject)
#' @param spat_enr_names which spatial enrichment results to include
#' @return data.table with average metadata values per (combined) annotation
#' @export
calculateMetaTableCells = function(gobject,
                                   value_cols = NULL,
                                   metadata_cols = NULL,
                                   spat_enr_names = NULL) {


  if(is.null(metadata_cols)) stop('\n You need to select one or more valid column names from pDataDT() \n')
  if(is.null(value_cols)) stop('\n You need to select one or more valid value column names from pDataDT() \n')

  cell_metadata = combineMetadata(gobject = gobject,
                                  spat_enr_names = spat_enr_names)

  ## only keep columns that exist
  cell_metadata_cols = colnames(cell_metadata)

  if(!all(value_cols %in% cell_metadata_cols)) {
    missing_value_cols = value_cols[!value_cols %in% cell_metadata_cols]
    cat('These value columns were not found: ', missing_value_cols)
  }
  value_cols = value_cols[value_cols %in% cell_metadata_cols]

  if(!all(metadata_cols %in% cell_metadata_cols)) {
    missing_metadata_cols = metadata_cols[!metadata_cols %in% cell_metadata_cols]
    cat('These metadata columns were not found: ', missing_metadata_cols)
  }
  metadata_cols = metadata_cols[metadata_cols %in% cell_metadata_cols]

  if(!length(metadata_cols) > 0 | !length(value_cols) > 0) {
    stop('\n missing sufficient metadata or value columns \n')
  }

  workdt = cell_metadata[, lapply(.SD, mean), by = metadata_cols, .SDcols = value_cols]
  workdtmelt = data.table::melt.data.table(workdt, measure.vars = value_cols)

  return(workdtmelt)

}




#' @title combineMetadata
#' @description This function combines the cell metadata with spatial locations and
#' enrichment results from \code{\link{runSpatialEnrich}}
#' @param gobject Giotto object
#' @param spat_enr_names names of spatial enrichment results to include
#' @return Extended cell metadata in data.table format.
#' @export
combineMetadata = function(gobject,
                           spat_enr_names = NULL) {

  # cell metadata
  metadata = pDataDT(gobject)

  # spatial locations
  spatial_locs = copy(gobject@spatial_locs)

  # data.table variables
  cell_ID = NULL

  metadata = cbind(metadata, spatial_locs[, cell_ID := NULL])

  # cell/spot enrichment data
  available_enr = names(gobject@spatial_enrichment)

  # output warning if not found
  not_available = spat_enr_names[!spat_enr_names %in% available_enr]
  if(length(not_available) > 0) {
    cat('These spatial enrichment results have not been found: \n',
        not_available)
  }

  spat_enr_names = spat_enr_names[spat_enr_names %in% available_enr]

  if(!is.null(spat_enr_names) & length(spat_enr_names) > 0) {

    result_list = list()
    for(spatenr in 1:length(spat_enr_names)) {

      spatenr_name = spat_enr_names[spatenr]
      temp_spat = copy(gobject@spatial_enrichment[[spatenr_name]])
      temp_spat[, 'cell_ID' := NULL]

      result_list[[spatenr]] = temp_spat
    }
    final_meta = do.call('cbind', c(list(metadata), result_list))

    duplicates = sum(duplicated(colnames(final_meta)))
    if(duplicates > 0) cat('Some column names are not unique.
                           If you add results from multiple enrichments,
                           consider giving the signatures unique names')

  } else {

    final_meta = metadata

  }

  return(final_meta)

}





#' @title createMetagenes
#' @description This function creates an average metagene for gene clusters.
#' @param gobject Giotto object
#' @param expression_values expression values to use
#' @param gene_clusters numerical vector with genes as names
#' @param name name of the metagene results
#' @param return_gobject return giotto object
#' @return giotto object
#' @details An example for the 'gene_clusters' could be like this:
#' cluster_vector = c(1, 1, 2, 2); names(cluster_vector) = c('geneA', 'geneB', 'geneC', 'geneD')
#' @export
#' @examples
#'
#' data(mini_giotto_single_cell)
#'
#' # get all genes
#' all_genes = slot(mini_giotto_single_cell, 'gene_ID')
#'
#' # create 2 metagenes from the first 6 genes
#' cluster_vector = c(1, 1, 1, 2, 2, 2) # 2 groups
#' names(cluster_vector) = all_genes[1:6]
#'
#' mini_giotto_single_cell = createMetagenes(mini_giotto_single_cell,
#'                                           gene_clusters = cluster_vector,
#'                                           name = 'cluster_metagene')
#'
#' # show metagene expression
#' spatCellPlot(mini_giotto_single_cell,
#'             spat_enr_names = 'cluster_metagene',
#'             cell_annotation_values = c('1', '2'),
#'             point_size = 3.5, cow_n_col = 2)
#'
createMetagenes = function(gobject,
                           expression_values = c('normalized', 'scaled', 'custom'),
                           gene_clusters,
                           name = 'metagene',
                           return_gobject = TRUE) {

  # expression values to be used
  values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))
  expr_values = select_expression_values(gobject = gobject, values = values)


  ## calculate metagene ##
  res_list = list()

  for(id in sort(unique(gene_clusters))) {

    clus_id = id

    selected_genes = names(gene_clusters[gene_clusters == clus_id])
    sub_mat = expr_values[rownames(expr_values) %in% selected_genes,]

    # calculate mean
    if(length(selected_genes) == 1) {
      mean_score = sub_mat
    } else{
      mean_score = colMeans_giotto(sub_mat)
    }

    res_list[[id]] = mean_score
  }

  res_final = data.table::as.data.table(t(do.call('rbind', res_list)))
  colnames(res_final) = as.character(sort(unique(gene_clusters)))

  # data.table variables
  cell_ID = NULL

  res_final[, cell_ID := colnames(expr_values)]


  if(return_gobject == TRUE) {

    ## enrichment scores
    spenr_names = names(gobject@spatial_enrichment)

    if(name %in% spenr_names) {
      cat('\n ', name, ' has already been used, will be overwritten \n')
    }

    ## update parameters used ##
    parameters_list = gobject@parameters
    number_of_rounds = length(parameters_list)
    update_name = paste0(number_of_rounds,'_create_metagene')

    # parameters to include
    parameters_list[[update_name]] = c('expression values' = values,
                                       'metagene name' = name)
    gobject@parameters = parameters_list

    gobject@spatial_enrichment[[name]] = res_final

    return(gobject)



  } else {
    return(res_final)
  }

}


#' @title findNetworkNeighbors
#' @description Find the spatial neighbors for a selected group of cells within the selected spatial network.
#' @param gobject Giotto object
#' @param spatial_network_name name of spatial network
#' @param source_cell_ids cell ids for which you want to know the spatial neighbors
#' @param name name of the results
#' @return data.table
#' @export
#' @examples
#'
#' data(mini_giotto_single_cell)
#'
#' # get all cells
#' all_cells = slot(mini_giotto_single_cell, 'cell_ID')
#'
#' # find all the spatial neighbours for the first 5 cells
#' # within the Delaunay network
#' findNetworkNeighbors(mini_giotto_single_cell,
#'                      spatial_network_name = 'Delaunay_network',
#'                      source_cell_ids = all_cells[1:5])
#'
findNetworkNeighbors = function(gobject,
                                spatial_network_name,
                                source_cell_ids = NULL,
                                name = 'nb_cells') {

  # get spatial network
  if(!is.null(spatial_network_name)) {
    spatial_network = select_spatialNetwork(gobject, name = spatial_network_name, return_network_Obj = FALSE)
  } else {
    stop('You need to select a spatial network')
  }

  # source cell ids that are found back
  all_cell_ids = gobject@cell_ID
  source_cells = all_cell_ids[all_cell_ids %in% source_cell_ids]

  if(length(source_cells) == 0) {
    stop('No source cell ids were selected or found')
  }


  full_network_DT = convert_to_full_spatial_network(spatial_network)
  potential_target_cells = full_network_DT[source %in% source_cells][['target']]
  source_and_target_cells = potential_target_cells[potential_target_cells %in% source_cells]
  target_cells = potential_target_cells[!potential_target_cells %in% source_and_target_cells]

  cell_meta = pDataDT(gobject)

  # data.table variables
  nb_cells = cell_ID = NULL

  cell_meta[, nb_cells := ifelse(cell_ID %in% source_and_target_cells, 'both',
                                 ifelse(cell_ID %in% source_cells, 'source',
                                        ifelse(cell_ID %in% target_cells, 'neighbor', 'others')))]
  nb_annot = cell_meta[, c('cell_ID', 'nb_cells'), with = FALSE]
  data.table::setnames(nb_annot, 'nb_cells', name)

  return(nb_annot)

}






