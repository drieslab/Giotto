
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
  return(as.data.table(fData(gobject)))

}


#' @title select_expression_values
#' @description helper function to select expression values
#' @param gobject giotto object
#' @param values expression values to extract
#' @return expression matrix
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


#' @title filterGiotto
#' @description filter Giotto object
#' @param gobject giotto object
#' @param expression_values expression values to use
#' @param expression_threshold threshold to consider a gene expressed
#' @param minimum_detected_genes minimum # of genes that need to be detected in a cell
#' @param minimum_expression_in_cell minimum # of cells that need to express a gene
#' @param verbose verbose
#' @return giotto object
#' @export
#' @examples
#'     filterGiotto(gobject)
filterGiotto <- function(gobject,
                         expression_values = c('raw', 'normalized', 'scaled', 'custom'),
                         expression_threshold = 1,
                         minimum_detected_genes = 100,
                         minimum_expression_in_cell = 100,
                         verbose = F) {

  # expression values to be used
  values = base::match.arg(expression_values, c('raw', 'normalized', 'scaled', 'custom'))
  expr_values = Giotto:::select_expression_values(gobject = gobject, values = values)

  ## filter cells
  filter_index_cells = base::colSums(expr_values >= expression_threshold) >= minimum_detected_genes
  selected_cell_ids = gobject@cell_ID[filter_index_cells]
  ## filter genes
  filter_index_genes = base::rowSums(expr_values >= expression_threshold) >= minimum_expression_in_cell
  selected_gene_ids = gobject@gene_ID[filter_index_genes]

  newGiottoObject = Giotto::subsetGiotto(gobject = gobject,
                                 cell_ids = selected_cell_ids,
                                 gene_ids = selected_gene_ids)

  ## print output ##
  removed_cells = base::length(filter_index_cells[filter_index_cells == FALSE])
  total_cells   = base::length(filter_index_cells)

  removed_genes = base::length(filter_index_genes[filter_index_genes == FALSE])
  total_genes   = base::length(filter_index_genes)

  if(verbose == TRUE) {
    cat('Number of cells removed: ', removed_cells, ' out of ', total_cells, '\n')
    cat('Number of genes removed: ', removed_genes, ' out of ', total_genes, '\n')
  }


  ## update parameters used ##
  parameters_list  = gobject@parameters
  number_of_rounds = base::length(parameters_list)
  update_name      = base::paste0(number_of_rounds,'_filter')
  # parameters to include
  parameters_list[[update_name]] = c('used expression values' = expression_values,
                                     'gene expression threshold' = expression_threshold,
                                     'minimum genes detected per cell' = minimum_detected_genes,
                                     'minimum times a gene is detected over all cells' = minimum_expression_in_cell)
  gobject@parameters = parameters_list

  return(newGiottoObject)


}


#' @title normalizeGiotto
#' @description normalize and/or scale expresion values of Giotto object
#' @param gobject giotto object
#' @param library_size_norm normalize cells by library size
#' @param scalefactor scale factor to use after library size normalization
#' @param logbase log base to use to log normalize expression values
#' @param scale_genes z-score genes over all cells
#' @param scale_cells z-score cells over all genes
#' @param scale_order order to scale genes and cells
#' @param verbose be verbose
#' @return giotto object
#' @details Description of normalization steps ...
#' @export
#' @examples
#'     normalizeGiotto(gobject)
normalizeGiotto <- function(gobject,
                            library_size_norm = T,
                            scalefactor = 6e3,
                            logbase = 2,
                            scale_genes = T,
                            scale_cells = T,
                            scale_order = c('first_genes', 'first_cells'),
                            verbose = F) {

  raw_expr = gobject@raw_exprs

    ## 1. library size normalize
  if(library_size_norm == TRUE) {
    norm_expr = base::t((base::t(raw_expr)/base::colSums(raw_expr))*scalefactor)

    if(!is.null(logbase)) {
      ## 2. lognormalize
      norm_expr = base::log(x = norm_expr+1, base = logbase)
    }
  }



  ## 3. scale
  if(scale_genes == TRUE & scale_cells == TRUE) {

    scale_order = match.arg(arg = scale_order, choices = c('first_genes', 'first_cells'))

    if(scale_order == 'first_genes') {
      if(verbose == TRUE) cat('\n first scale genes and then cells \n')
      norm_scaled_expr = base::t(scale(x = base::t(norm_expr)))
      norm_scaled_expr = base::scale(x = norm_scaled_expr)
    } else if(scale_order == 'first_cells') {
      if(verbose == TRUE) cat('\n first scale cells and then genes \n')
      norm_scaled_expr = base::scale(x = norm_expr)
      norm_scaled_expr = base::t(scale(x = base::t(norm_scaled_expr)))
    } else {
      stop('\n scale order must be given \n')
    }

  } else if(scale_genes == TRUE) {
    norm_scaled_expr = base::t(scale(x = base::t(norm_expr)))
  } else if(scale_cells == TRUE) {
    norm_scaled_expr = base::scale(x = norm_expr)
  } else {
    norm_scaled_expr = NULL
  }

  # return Giotto object
  gobject@norm_expr = norm_expr
  gobject@norm_scaled_expr = norm_scaled_expr

  ## update parameters used ##
  parameters_list  = gobject@parameters
  number_of_rounds = base::length(parameters_list)
  update_name      = base::paste0(number_of_rounds,'_normalize')
  # parameters to include
  parameters_list[[update_name]] = c('normalized to library size' = ifelse(library_size_norm == T, 'yes', 'no'),
                                     'scalefactor' = scalefactor,
                                     'logbase' = ifelse(is.null(logbase), NA, logbase),
                                     'genes scaled' = ifelse(scale_genes == T, 'yes', 'no'),
                                     'cell scaled' = ifelse(scale_cells == T, 'yes', 'no'))
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
    if(!all(batch_columns %in% base::colnames(cell_metadata))) {
      stop('\n batch column name(s) were not found in the cell metadata \n')
    }
  }

  if(!is.null(covariate_columns)) {
    if(!all(covariate_columns %in% base::colnames(cell_metadata))) {
      stop('\n covariate column name(s) were not found in the cell metadata \n')
    }
  }

  update_slot = base::match.arg(update_slot, c('normalized', 'scaled', 'custom'))

  values = base::match.arg(expression_values, c('normalized', 'scaled', 'custom'))
  expr_data = Giotto:::select_expression_values(gobject = gobject, values = values)

  # batch columns
  if(!base::is.null(batch_columns)) {
    batch_column_1 = cell_metadata[[ batch_columns[1] ]]
    if(base::length(batch_columns) > 1) {
      batch_column_2 = cell_metadata[[ batch_columns[2] ]]
    } else {
      batch_column_2 = NULL
    }
  } else {
    batch_column_1 = NULL
    batch_column_2 = NULL
  }

  # covariate columns
  if(!base::is.null(covariate_columns)) {
    covariates = base::as.matrix(cell_metadata[, covariate_columns, with = F])
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
  if(name %in% colnames(cell_metadata)) {
    cat('\n annotation name ', name,' was already used \n',
        'and will be overwritten \n')
    cell_metadata[, (name) := NULL]
  }

  cell_metadata[, temp_cluster_name := annotation_vector[[get(cluster_column)]], by = 1:nrow(cell_metadata)]
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
    cell_metadata = base::cbind(cell_metadata, new_metadata)
  } else {
    if(base::is.null(column_cell_ID)) stop('You need to provide cell_ID column')
    cell_metadata <- data.table:::merge.data.table(cell_metadata, by.x = 'cell_ID',
                                                   new_metadata, by.y = column_cell_ID,
                                                   all.x = T)
  }

  # reorder
  cell_metadata = cell_metadata[base::match(ordered_cell_IDs, cell_ID)]

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
    gene_metadata = base::cbind(gene_metadata, new_metadata)
  } else {
    if(base::is.null(column_gene_ID)) stop('You need to provide gene_ID column')
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
  values = base::match.arg(expression_values, c('normalized', 'scaled', 'custom'))
  expr_data = Giotto:::select_expression_values(gobject = gobject, values = values)

  # calculate stats
  gene_stats = data.table(genes = base::rownames(expr_data),
                          nr_cells = base::rowSums(expr_data > detection_threshold),
                          perc_cells = (base::rowSums(expr_data > detection_threshold)/base::ncol(expr_data))*100,
                          total_expr = base::rowSums(expr_data),
                          mean_expr = base::rowMeans(expr_data))

    mean_expr_detected = base::unlist(base::apply(X = expr_data, MARGIN = 1, FUN = function(x) {
    detected_x = x[x > detection_threshold]
    mean(detected_x)
  }))
  gene_stats[, mean_expr_det := mean_expr_detected]


  if(return_gobject == TRUE) {

    # remove previous statistics
    gene_metadata = Giotto::fDataDT(gobject)
    metadata_names = base::colnames(gene_metadata)
    if('nr_cells' %in% metadata_names) {
      base::cat('\n gene statistics has already been applied once, will be overwritten \n')
      gene_metadata[, c('nr_cells', 'perc_cells', 'total_expr', 'mean_expr', 'mean_expr_det') := NULL]
      gobject@gene_metadata = gene_metadata
    }

    gobject = Giotto::addGeneMetadata(gobject = gobject, new_metadata = gene_stats,
                                      by_column = T, column_gene_ID = 'genes')

    ## update parameters used ##
    parameters_list = gobject@parameters
    number_of_rounds = base::length(parameters_list)
    update_name = base::paste0(number_of_rounds,'_gene_stats')
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
  values = base::match.arg(expression_values, c('normalized', 'scaled', 'custom'))
  expr_data = Giotto:::select_expression_values(gobject = gobject, values = values)

  # calculate stats
  cell_stats = data.table(cells = base::colnames(expr_data),
                          nr_genes = base::colSums(expr_data > detection_threshold),
                          perc_genes = (base::colSums(expr_data > detection_threshold)/base::nrow(expr_data))*100,
                          total_expr = base::colSums(expr_data))



  if(return_gobject == TRUE) {

    # remove previous statistics
    cell_metadata = Giotto::pDataDT(gobject)
    metadata_names = base::colnames(cell_metadata)
    if('nr_genes' %in% metadata_names) {
      cat('\n cells statistics has already been applied once, will be overwritten \n')
      cell_metadata[, c('nr_genes', 'perc_genes', 'total_expr') := NULL]
      gobject@cell_metadata = cell_metadata
    }

    gobject = Giotto::addCellMetadata(gobject = gobject, new_metadata = cell_stats,
                                      by_column = T, column_cell_ID = 'cells')

    ## update parameters used ##
    parameters_list = gobject@parameters
    number_of_rounds = base::length(parameters_list)
    update_name = base::paste0(number_of_rounds,'_cell_stats')
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

  base::cat('Processing steps: \n \n')

  for(step in base::names(parameters)) {
    base::cat('\n', step, '\n')

    sub_step = parameters[[step]]

    if(base::any(base::grepl('name', base::names(sub_step)) == TRUE)) {

      selected_names = base::grep('name', base::names(sub_step), value = T)
      base::cat('\t name info: ', sub_step[selected_names], '\n')

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

  values = base::match.arg(expression_values, c('normalized', 'scaled', 'custom'))

  # average expression per cluster
  aggr_sc_clusters <- Giotto:::create_average_DT(gobject = gobject, meta_data_name = cluster_column,
                                        expression_values = values)
  aggr_sc_clusters_DT <- data.table::as.data.table(aggr_sc_clusters)
  aggr_sc_clusters_DT[, genes := base::rownames(aggr_sc_clusters)]
  aggr_sc_clusters_DT_melt <- data.table::melt.data.table(aggr_sc_clusters_DT,
                                                          variable.name = 'cluster',
                                                          id.vars = 'genes',
                                                          value.name = 'expression')

  # create matrix
  testmat = data.table::dcast.data.table(aggr_sc_clusters_DT_melt,
                                         formula = genes~cluster, value.var = 'expression')
  testmatrix = base::as.matrix(testmat[,-1])
  base::rownames(testmatrix) = testmat[['genes']]

  # create subset if required
  if(!base::is.null(gene_subset)) {
    gene_subset_detected = gene_subset[gene_subset %in% base::rownames(testmatrix)]
    testmatrix = testmatrix[base::rownames(testmatrix) %in% gene_subset_detected, ]
  }

  return(testmatrix)

}



