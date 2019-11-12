
#' @title pDataDT
#' @description show cell metadata
#' @param gobject giotto object
#' @return data.table
#' @export
#' @examples
#'     pDataDT(gobject)
pDataDT <- function(gobject) {

  if(!class(gobject) %in% c('ExpressionSet', 'SCESet', 'seurat', 'giotto')) {
    stop('only works with ExpressionSet (-like) objects')
  }

  if(class(gobject) %in% c('ExpressionSet', 'SCESet')) {
    return(as.data.table(pData(gobject)))
  }
  else if(class(gobject) == 'giotto') {
    return(gobject@cell_metadata)
  }
  else if(class(gobject) == 'seurat') {
    return(as.data.table(gobject@meta.data))
  }

}

#' @title fDataDT
#' @description show gene metadata
#' @param gobject giotto object
#' @return data.table
#' @export
#' @examples
#'     pDataDT(gobject)
fDataDT <- function(gobject) {

  if(!class(gobject) %in% c('ExpressionSet', 'SCESet', 'giotto')) {
    stop('only works with ExpressionSet (-like) objects')
  }
  else if(class(gobject) == 'giotto') {
    return(gobject@gene_metadata)
  }
  return(data.table::as.data.table(fData(gobject)))

}



#' @title select_expression_values
#' @description helper function to select expression values
#' @param gobject giotto object
#' @param values expression values to extract
#' @return expression matrix
#' @export
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
create_average_DT <- function(gobject, meta_data_name,
                              expression_values = c('normalized', 'scaled', 'custom')) {


  # expression values to be used
  values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))
  expr_data = select_expression_values(gobject = gobject, values = values)

  # metadata
  cell_metadata <- pDataDT(gobject)
  myrownames <- rownames(expr_data)

  savelist <- list()
  for(group in unique(cell_metadata[[meta_data_name]])) {

    name = paste0('cluster_', group)

    temp <- expr_data[, cell_metadata[[meta_data_name]] == group]
    temp_DT <- rowMeans(as.matrix(temp))

    savelist[[name]] <- temp_DT
  }

  finalDF <- do.call('cbind', savelist)
  rownames(finalDF) <- myrownames

  return(as.data.frame(finalDF))
}

#' @title create_average_detection_DT
#' @description calculates average gene detection for a cell metadata factor (e.g. cluster)
#' @param gobject giotto object
#' @param meta_data_name name of metadata column to use
#' @param expression_values which expression values to use
#' @param detection_threshold detection threshold to consider a gene detected
#' @return data.table with average gene epression values for each factor
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

    temp <- expr_data[, cell_metadata[[meta_data_name]] == group]

    if(is.matrix(temp)) {
      temp_DT = rowSums(as.matrix(temp) > detection_threshold)/ncol(temp)
    } else {
      temp_DT = as.numeric(temp > detection_threshold)
    }

    savelist[[name]] <- temp_DT
  }

  finalDF <- do.call('cbind', savelist)
  rownames(finalDF) <- myrownames

  return(as.data.frame(finalDF))
}





#' @title subsetGiotto
#' @description subsets Giotto object including previous calculations
#' @param gobject giotto object
#' @param cell_ids cell IDs to keep
#' @param gene_ids gene IDs to keep
#' @return giotto object
#' @export
#' @examples
#'     subsetGiotto(gobject)
subsetGiotto <- function(gobject, cell_ids = NULL, gene_ids = NULL) {


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

  ## spatial network & grid ##
  # cell spatial network
  if(!is.null(gobject@spatial_network)) {
    for(network in names(gobject@spatial_network)) {
      gobject@spatial_network[[network]] =   gobject@spatial_network[[network]][to %in% cells_to_keep & from %in% cells_to_keep]
    }
  }

  # spatial grid
  # need to be recomputed


  ## dimension reduction ##
  # gene dim reduction
  # TODO

  # cell dim reduction
  if(!is.null(gobject@dimension_reduction$cells)) {

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

  }


  ## nn network ##
  if(!is.null(gobject@nn_network$cells)) {
    print('\n subset networks \n')

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



#' @title filterDistributions
#' @description show gene or cell filter distributions
#' @param gobject giotto object
#' @param expression_values expression values to use
#' @param expression_threshold threshold to consider a gene expressed
#' @param detection look at genes or cells
#' @param plot_type type of plot
#' @param nr_bins number of bins for histogram plot
#' @param fill_color fill color for plots
#' @param scale_axis ggplot transformation for axis (e.g. log2)
#' @param axis_offset offset to be used together with the scaling transformation
#' @param show_plot show plot
#' @return ggplot object
#' @export
#' @examples
#'     filterDistributions(gobject)
filterDistributions <- function(gobject,
                                expression_values = c('raw', 'normalized', 'scaled', 'custom'),
                                expression_threshold = 1,
                                detection = c('genes', 'cells'),
                                plot_type = c('histogram', 'violin'),
                                nr_bins = 30,
                                fill_color = 'lightblue',
                                scale_axis = 'identity',
                                axis_offset = 0,
                                show_plot = TRUE) {

  # expression values to be used
  values = match.arg(expression_values, c('raw', 'normalized', 'scaled', 'custom'))
  expr_values = select_expression_values(gobject = gobject, values = values)

  # plot distribution for genes or cells
  detection = match.arg(detection, c('genes', 'cells'))

  # plot type
  plot_type = match.arg(plot_type, c('histogram', 'violin'))

  # for genes
  if(detection == 'genes') {

    gene_detection_levels = data.table::as.data.table(rowSums(expr_values >= expression_threshold))

    if(plot_type == 'violin') {

      pl <- ggplot2::ggplot()
      pl <- pl + ggplot2::theme_classic()
      pl <- pl + ggplot2::geom_violin(data = gene_detection_levels, aes(x = 'genes', y = V1+axis_offset),
                                      fill = fill_color)
      pl <- pl + ggplot2::scale_y_continuous(trans = scale_axis)
      pl <- pl + ggplot2::labs(y = 'gene detected in # of cells', x = '')

    } else if(plot_type == 'histogram') {

      pl <- ggplot2::ggplot()
      pl <- pl + ggplot2::theme_classic()
      pl <- pl + ggplot2::geom_histogram(data = gene_detection_levels, aes(x = V1+axis_offset),
                                         color = 'white', bins = nr_bins, fill = fill_color)
      pl <- pl + ggplot2::scale_x_continuous(trans = scale_axis)
      pl <- pl + ggplot2::labs(x = 'gene detected in # of cells')

    }

    # for cells
  } else if(detection == 'cells') {

    cell_detection_levels = data.table::as.data.table(colSums(expr_values >= expression_threshold))

    if(plot_type == 'violin') {

      pl <- ggplot2::ggplot()
      pl <- pl + ggplot2::theme_classic()
      pl <- pl + ggplot2::geom_violin(data = cell_detection_levels, aes(x = 'cells', y = V1+axis_offset),
                                      fill = fill_color)
      pl <- pl + ggplot2::scale_y_continuous(trans = scale_axis)
      pl <- pl + ggplot2::labs(y = 'genes detected per cell', x = '')

    } else if(plot_type == 'histogram') {

      pl <- ggplot2::ggplot()
      pl <- pl + ggplot2::theme_classic()
      pl <- pl + ggplot2::geom_histogram(data = cell_detection_levels, aes(x = V1+axis_offset),
                                         color = 'white', bins = nr_bins, fill = fill_color)
      pl <- pl + ggplot2::scale_x_continuous(trans = scale_axis)
      pl <- pl + ggplot2::labs(x = 'genes detected per cell')

    }
  }

  if(show_plot == TRUE) {
    print(pl)
  }

  return(pl)

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
#' @return list of data.table and ggplot object
#' @details Creates a scatterplot that visualizes the number of genes and cells that are
#' lost with a specific combination of a gene and cell threshold given an arbitrary cutoff
#' to call a gene expressed. This function can be used to make an informed decision at the
#' filtering step with filterGiotto.
#' @export
#' @examples
#'     filterCombinations(gobject)
filterCombinations <- function(gobject,
                               expression_values = c('raw', 'normalized', 'scaled', 'custom'),
                               expression_thresholds = c(1, 2),
                               gene_det_in_min_cells = c(5, 50),
                               min_det_genes_per_cell = c(200, 400),
                               scale_x_axis = 'identity',
                               x_axis_offset = 0,
                               scale_y_axis = 'identity',
                               y_axis_offset = 0,
                               show_plot = TRUE) {


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
      filter_index_genes = rowSums(expr_values >= threshold) >= min_cells_for_gene
      removed_genes = length(filter_index_genes[filter_index_genes == FALSE])
      det_cells_res[[combn_i]] = removed_genes

      # then remove cells
      filter_index_cells = colSums(expr_values[filter_index_genes, ] >= threshold) >= min_genes_per_cell
      removed_cells = length(filter_index_cells[filter_index_cells == FALSE])
      det_genes_res[[combn_i]] = removed_cells
    }

    temp_dt = data.table::data.table('threshold' = threshold,
                                     removed_genes = unlist(det_cells_res),
                                     removed_cells = unlist(det_genes_res))

    result_list[[thresh_i]] = temp_dt

  }

  result_DT = do.call('rbind', result_list)
  result_DT[['gene_detected_in_min_cells']] = gene_det_in_min_cells
  result_DT[['min_detected_genes_per_cell']] = min_det_genes_per_cell
  result_DT[['combination']] = paste0(result_DT$gene_detected_in_min_cells,'-',result_DT$min_detected_genes_per_cell)

  result_DT = result_DT[,.(threshold,
                           gene_detected_in_min_cells, min_detected_genes_per_cell,
                           combination, removed_genes, removed_cells)]

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
  if(show_plot == TRUE) {
    print(pl)
  }

  return(list(results = result_DT, ggplot = pl))

}


#' @title filterGiotto
#' @description filter Giotto object
#' @param gobject giotto object
#' @param expression_values expression values to use
#' @param expression_threshold threshold to consider a gene expressed
#' @param gene_det_in_min_cells minimum # of cells that need to express a gene
#' @param min_det_genes_per_cell minimum # of genes that need to be detected in a cell
#' @param verbose verbose
#' @return giotto object
#' @export
#' @examples
#'     filterGiotto(gobject)
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
  filter_index_genes = rowSums(expr_values >= expression_threshold) >= gene_det_in_min_cells
  selected_gene_ids = gobject@gene_ID[filter_index_genes]

  ## filter cells
  filter_index_cells = colSums(expr_values[filter_index_genes, ] >= expression_threshold) >= min_det_genes_per_cell
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
  parameters_list[[update_name]] = c('used expression values' = expression_values,
                                     'gene expression threshold' = expression_threshold,
                                     'minimum # of genes detected per cell' = min_det_genes_per_cell,
                                     'minimum times a gene is detected over all cells' = gene_det_in_min_cells)
  newGiottoObject@parameters = parameters_list

  return(newGiottoObject)


}


#' @title normalizeGiotto
#' @description normalize and/or scale expresion values of Giotto object
#' @param gobject giotto object
#' @param norm_methods normalization method to use
#' @param library_size_norm normalize cells by library size
#' @param scalefactor scale factor to use after library size normalization
#' @param log_norm transform values to log-scale
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
#' 1. Data normalization for total library size and scaling by a custom scale-factor. \cr
#' 2. Log transformation of data. \cr
#' 3. Z-scoring of data by genes and/or cells.
#'
#' B. The normalization method as provided by the osmFISH paper is also implemented: \cr
#' 1. First normalize genes, for each gene divide the counts by the total gene count and
#' multiply by the total number of genes. \cr
#' 2. Next normalize cells, for each cell divide the normalized gene counts by the total
#' counts per cell and multiply by the total number of cells. \cr
#' This data will be saved in the Giotto slot for custom expression.
#' @export
#' @examples
#'     normalizeGiotto(gobject)
normalizeGiotto <- function(gobject,
                            norm_methods = c('standard', 'osmFISH'),
                            library_size_norm = TRUE,
                            scalefactor = 6e3,
                            log_norm = TRUE,
                            logbase = 2,
                            scale_genes = T,
                            scale_cells = T,
                            scale_order = c('first_genes', 'first_cells'),
                            verbose = F) {

  raw_expr = gobject@raw_exprs

  norm_methods = match.arg(arg = norm_methods, choices = c('standard', 'osmFISH'))

  # normalization according to standard methods
  if(norm_methods == 'standard') {

    ## 1. library size normalize
    if(library_size_norm == TRUE) {
      norm_expr = t((t(raw_expr)/colSums(raw_expr))*scalefactor)
    } else {
      norm_expr = raw_expr
    }

    ## 2. lognormalize
    if(log_norm == TRUE) {
      norm_expr = log(x = norm_expr+1, base = logbase)
    } else {
      norm_expr = norm_expr
    }

    ## 3. scale
    if(scale_genes == TRUE & scale_cells == TRUE) {

      scale_order = match.arg(arg = scale_order, choices = c('first_genes', 'first_cells'))

      if(scale_order == 'first_genes') {
        if(verbose == TRUE) cat('\n first scale genes and then cells \n')
        norm_scaled_expr = t(scale(x = t(norm_expr)))
        norm_scaled_expr = scale(x = norm_scaled_expr)
      } else if(scale_order == 'first_cells') {
        if(verbose == TRUE) cat('\n first scale cells and then genes \n')
        norm_scaled_expr = scale(x = norm_expr)
        norm_scaled_expr = t(scale(x = t(norm_scaled_expr)))
      } else {
        stop('\n scale order must be given \n')
      }

    } else if(scale_genes == TRUE) {
      norm_scaled_expr = t(scale(x = t(norm_expr)))
    } else if(scale_cells == TRUE) {
      norm_scaled_expr = scale(x = norm_expr)
    } else {
      norm_scaled_expr = NULL
    }


    ## 4. reverse log-scale
    # only when data have been logged
    # and when data have been scaled
    # not implemented


    # return Giotto object
    gobject@norm_expr = norm_expr
    gobject@norm_scaled_expr = norm_scaled_expr

  }

  # normalization according to osmFISH method
  else if(norm_methods == 'osmFISH') {
    # 1. normalize per gene with scale-factor equal to number of genes
    norm_genes = (raw_expr/rowSums(raw_expr)) * nrow(raw_expr)
    # 2. normalize per cells with scale-factor equal to number of cells
    norm_genes_cells = t((t(norm_genes)/colSums(norm_genes)) * ncol(raw_expr))

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
#' @description normalize and/or scale expresion values of Giotto object
#' @param gobject giotto object
#' @param expression_values expression values to use
#' @param batch_columns metadata columns that represent different batch
#' @param covariate_columns metadata columns that represent covariates to regress out
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @param update_slot expression slot that will be updated (default = custom)
#' @return giotto object
#' @details Description of adjusting steps ...
#' @export
#' @examples
#'     adjustGiottoMatrix(gobject)
adjustGiottoMatrix <- function(gobject,
                               expression_values = c('normalized', 'scaled', 'custom'),
                               batch_columns = NULL,
                               covariate_columns = NULL,
                               return_gobject = TRUE,
                               update_slot = c('custom')) {

  # metadata
  cell_metadata = Giotto::pDataDT(gobject)

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



#' @title annotateGiotto
#' @description adds cell annotation to giotto object based on clustering
#' @param gobject giotto object
#' @param annotation_vector named annotation vector (names = cluster ids)
#' @param cluster_column cluster column to convert to annotation names
#' @param name new name for annotation column
#' @return giotto object
#' @details Description of how to add cell metadata ...
#' @export
#' @examples
#'     annotateGiotto(gobject)
annotateGiotto <- function(gobject, annotation_vector = NULL, cluster_column = NULL, name = 'cell_types') {

  if(is.null(annotation_vector) | is.null(cluster_column)) {
    stop('\n You need to provide both a named annotation vector and the corresponding cluster column  \n')
  }


  cell_metadata = pDataDT(gobject)

  # 1. verify if cluster column exist
  if(!cluster_column %in% colnames(cell_metadata)) {
    stop('\n Cluster column is not found in cell metadata \n')
  }

  # 2. remove previous annotation name if it's the same
  # but only if new name is not the same as cluster to be used
  if(name %in% colnames(cell_metadata)) {
    cat('\n annotation name ', name,' was already used \n',
        'and will be overwritten \n')

    cell_metadata[, temp_cluster_name := annotation_vector[[get(cluster_column)]], by = 1:nrow(cell_metadata)]
    cell_metadata[, (name) := NULL]

  } else {

    cell_metadata[, temp_cluster_name := annotation_vector[[get(cluster_column)]], by = 1:nrow(cell_metadata)]
  }

  setnames(cell_metadata, old = 'temp_cluster_name', new = name)
  gobject@cell_metadata = cell_metadata

  return(gobject)


}



#' @title addCellMetadata
#' @description adds cell metadata to the giotto object
#' @param gobject giotto object
#' @param new_metadata new metadata to use
#' @param by_column merge metadata based on cell_ID column in pDataDT
#' @param column_cell_ID column name of new metadata to use if by_column = TRUE
#' @return giotto object
#' @details Description of how to add cell metadata ...
#' @export
#' @examples
#'     addCellMetadata(gobject)
addCellMetadata <- function(gobject,
                            new_metadata,
                            by_column = F,
                            column_cell_ID = NULL) {

  cell_metadata = gobject@cell_metadata
  ordered_cell_IDs = gobject@cell_ID

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
    cell_metadata <- data.table:::merge.data.table(cell_metadata, by.x = 'cell_ID',
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
#' @param column_cell_ID column name of new metadata to use if by_column = TRUE
#' @return giotto object
#' @details Description of how to add gene metadata ...
#' @export
#' @examples
#'     addGeneMetadata(gobject)
addGeneMetadata <- function(gobject,
                            new_metadata,
                            by_column = F,
                            column_gene_ID = NULL) {

  gene_metadata = gobject@gene_metadata
  ordered_gene_IDs = gobject@gene_ID

  if(by_column == FALSE) {
    gene_metadata = cbind(gene_metadata, new_metadata)
  } else {
    if(is.null(column_gene_ID)) stop('You need to provide gene_ID column')
    gene_metadata <- data.table:::merge.data.table(gene_metadata, by.x = 'gene_ID',
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
#' @details Details about gene statistics that are returned.
#' @export
#' @examples
#'     addGeneStatistics(gobject)
addGeneStatistics <- function(gobject,
                              expression_values = c('normalized', 'scaled', 'custom'),
                              detection_threshold = 0,
                              return_gobject = TRUE) {

  # expression values to be used
  expression_values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))
  expr_data = select_expression_values(gobject = gobject, values = expression_values)

  # calculate stats
  gene_stats = data.table::data.table(genes = rownames(expr_data),
                          nr_cells = rowSums(expr_data > detection_threshold),
                          perc_cells = (rowSums(expr_data > detection_threshold)/ncol(expr_data))*100,
                          total_expr = rowSums(expr_data),
                          mean_expr = rowMeans(expr_data))

    mean_expr_detected = unlist(apply(X = expr_data, MARGIN = 1, FUN = function(x) {
    detected_x = x[x > detection_threshold]
    mean(detected_x)
  }))
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
#' @details Details about cell statistics that are returned.
#' @export
#' @examples
#'     addCellStatistics(gobject)
addCellStatistics <- function(gobject,
                              expression_values = c('normalized', 'scaled', 'custom'),
                              detection_threshold = 0,
                              return_gobject = TRUE) {

  # expression values to be used
  expression_values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))
  expr_data = select_expression_values(gobject = gobject, values = expression_values)

  # calculate stats
  cell_stats = data.table::data.table(cells = colnames(expr_data),
                          nr_genes = colSums(expr_data > detection_threshold),
                          perc_genes = (colSums(expr_data > detection_threshold)/nrow(expr_data))*100,
                          total_expr = colSums(expr_data))



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
#' @details Details about gene and cell statistics that are returned.
#' @export
#' @examples
#'     addStatistics(gobject)
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



#' @title showProcessingSteps
#' @description shows the sequential processing steps that were performed
#' @param gobject giotto object
#' @return list of processing steps and names
#' @export
#' @examples
#'     showProcessingSteps(gobject)
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
#' @description creates aggregated matrix for a given clustering
#' @examples
#'     create_cluster_matrix(gobject)
create_cluster_matrix <- function(gobject,
                                  expression_values = c('normalized', 'scaled', 'custom'),
                                  cluster_column,
                                  gene_subset = NULL) {

  values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))

  # average expression per cluster
  aggr_sc_clusters <- create_average_DT(gobject = gobject, meta_data_name = cluster_column,
                                        expression_values = values)
  aggr_sc_clusters_DT <- data.table::as.data.table(aggr_sc_clusters)
  aggr_sc_clusters_DT[, genes := rownames(aggr_sc_clusters)]
  aggr_sc_clusters_DT_melt <- data.table::melt.data.table(aggr_sc_clusters_DT,
                                                          variable.name = 'cluster',
                                                          id.vars = 'genes',
                                                          value.name = 'expression')

  # create matrix
  testmat = data.table::dcast.data.table(aggr_sc_clusters_DT_melt,
                                         formula = genes~cluster, value.var = 'expression')
  testmatrix = as.matrix(testmat[,-1])
  rownames(testmatrix) = testmat[['genes']]

  # create subset if required
  if(!is.null(gene_subset)) {
    gene_subset_detected = gene_subset[gene_subset %in% rownames(testmatrix)]
    testmatrix = testmatrix[rownames(testmatrix) %in% gene_subset_detected, ]
  }

  return(testmatrix)

}


#' @title dt_to_matrix
#' @description converts data.table to matrix
#' @examples
#'     dt_to_matrix(x)
dt_to_matrix <- function(x) {
  rownames = as.character(x[[1]])
  mat = as.matrix(x[,-1])
  rownames(mat) = rownames
  return(mat)
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
#'     calculateMetaTable(gobject)
calculateMetaTable = function(gobject,
                              expression_values =  c("normalized", "scaled", "custom"),
                              metadata_cols = NULL,
                              selected_genes = NULL) {

  if(is.null(metadata_cols)) stop('\n You need to select one or more valid column names from pDataDT() \n')

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
      subvec = rowMeans(sub_expr_values)
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

