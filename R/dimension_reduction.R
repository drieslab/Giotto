



#' @title create_dimObject
#' @name create_dimObject
#' @description Creates an object that stores a dimension reduction output
#' @param name arbitrary name for object
#' @param reduction_method method used to reduce dimensions
#' @param coordinates accepts the coordinates after dimension reduction
#' @param misc any additional information will be added to this slot
#' @return number of distinct colors
create_dimObject = function(name = 'test',
                            reduction_method = NULL,
                            coordinates = NULL,
                            misc = NULL) {


  number_of_dimensions = ncol(coordinates)
  colnames(coordinates) <- paste0('Dim.',1:number_of_dimensions)


  dimObj = list(name = name,
                reduction_method = reduction_method,
                coordinates = coordinates,
                misc = misc)

  class(dimObj) <- append(class(dimObj), 'dimObj')
  return(dimObj)

}


#' @title runPCA
#' @name runPCA
#' @description run PCA
#' @param gobject giotto object
#' @param expression_values expression values to use
#' @param reduction cells or genes
#' @param name arbitrary name for PCA run
#' @param genes_to_use subset of genes to use for PCA
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @param scale.unit scale features before PCA
#' @param ncp number of principal components to calculate
#' @param ... additional parameters for PCA
#' @return giotto object with updated PCA dimension recuction
#' @details Description of PCA steps...
#' @export
#' @examples
#'     runPCA(gobject)
runPCA <- function(gobject,
                   expression_values = c('normalized', 'scaled', 'custom'),
                   reduction = c('cells', 'genes'),
                   name = 'pca',
                   genes_to_use = NULL,
                   return_gobject = TRUE,
                   scale.unit = F,
                   ncp = 200,
                   ...) {


  # expression values to be used
  values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))
  expr_values = select_expression_values(gobject = gobject, values = values)


  # subset expression matrix
  if(!is.null(genes_to_use)) {
    expr_values = expr_values[rownames(expr_values) %in% genes_to_use, ]
  }

  # do PCA dimension reduction
  reduction = match.arg(reduction, c('cells', 'genes'))
  if(reduction == 'cells') {
    pca_object <- FactoMineR::PCA(X = t(expr_values),
                                  scale.unit = scale.unit, ncp = ncp, graph = F, ...)
  } else {
    pca_object <- FactoMineR::PCA(X = expr_values,
                                  scale.unit = scale.unit, ncp = ncp, graph = F, ...)
  }

  if(return_gobject == TRUE) {

    pca_names = names(gobject@dimension_reduction[[reduction]][['pca']])

    if(name %in% pca_names) {
      cat('\n ', name, ' has already been used, will be overwritten \n')

    }

    dimObject = create_dimObject(name = name, reduction_method = 'pca',
                                 coordinates = pca_object$ind$coord,
                                 misc = pca_object)

    gobject@dimension_reduction[[reduction]][['pca']][[name]] <- dimObject



    ## update parameters used ##
    parameters_list = gobject@parameters
    number_of_rounds = length(parameters_list)
    update_name = paste0(number_of_rounds,'_pca')
    # parameters to include
    parameters_list[[update_name]] = c('reduction type:' = reduction,
                                       'expression values' = expression_values,
                                       'number of genes used:' = length(genes_to_use),
                                       'ncp' = ncp,
                                       'scale.unit' = scale.unit,
                                       'name for pca' = name)
    gobject@parameters = parameters_list

    return(gobject)


  } else {
    return(pca_object)
  }
}



#' @title runUMAP
#' @name runUMAP
#' @description run UMAP
#' @param gobject giotto object
#' @param expression_values expression values to use
#' @param reduction cells or genes
#' @param dim_reduction_to_use use another dimension reduction set as input
#' @param dim_reduction_name name of dimension reduction set to use
#' @param dimensions_to_use number of dimensions to use as input
#' @param name arbitrary name for UMAP run
#' @param genes_to_use if dim_reduction_to_use = NULL, which genes to use
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @param n_neighbors UMAP param: number of neighbors
#' @param n_components UMAP param: number of components
#' @param n_epochs UMAP param: number of epochs
#' @param min_dist UMAP param: minimum distance
#' @param n_threads UMAP param: threads to use
#' @param spread UMAP param: spread
#' @param set_seed use of seed
#' @param seed_number seed number to use
#' @param ... additional UMAP parameters
#' @return giotto object with updated UMAP dimension recuction
#' @details Description of UMAP steps...
#' @export
#' @examples
#'     runUMAP(gobject)
runUMAP <- function(gobject,
                    expression_values = c('normalized', 'scaled', 'custom'),
                    reduction = c('cells', 'genes'),
                    dim_reduction_to_use = 'pca',
                    dim_reduction_name = 'pca',
                    dimensions_to_use = 1:10,
                    name = 'umap',
                    genes_to_use = NULL,
                    return_gobject = TRUE,
                    n_neighbors = 40,
                    n_components = 2,
                    n_epochs = 400,
                    min_dist = 0.01,
                    n_threads = 4,
                    spread = 5,
                    set_seed = T,
                    seed_number = 1234,
                    ...) {

  reduction = match.arg(reduction, choices = c('cells', 'genes'))

  ## umap on cells ##
  if(reduction == 'cells') {

    ## using dimension reduction ##
    if(!is.null(dim_reduction_to_use)) {

      ## TODO: check if reduction exists

      matrix_to_use = gobject@dimension_reduction[['cells']][[dim_reduction_to_use]][[dim_reduction_name]][['coordinates']][, dimensions_to_use]


    } else {
      ## using original matrix ##
      values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))
      expr_values = select_expression_values(gobject = gobject, values = values)

      # subset expression matrix
      if(!is.null(genes_to_use)) {
        expr_values = expr_values[rownames(expr_values) %in% genes_to_use, ]
      }

      matrix_to_use = t(expr_values)

    }


    ## run umap ##
    if(set_seed == TRUE) {
      set.seed(seed = seed_number)
    }

    uwot_clus <- uwot::umap(X = matrix_to_use, n_neighbors = n_neighbors, n_components = n_components,
                            n_epochs = n_epochs, min_dist = min_dist, n_threads = n_threads, spread = spread, ...)
    uwot_clus_pos_DT <- as.data.table(uwot_clus)
    uwot_clus_pos_DT[, cell_ID := rownames(matrix_to_use)]


    if(set_seed == TRUE) {
      set.seed(Sys.time())
    }


    if(return_gobject == TRUE) {

      umap_names = names(gobject@dimension_reduction[[reduction]][['umap']])

      if(name %in% umap_names) {
        cat('\n ', name, ' has already been used, will be overwritten \n')

      }


      coordinates = uwot_clus
      rownames(coordinates) = rownames(matrix_to_use)

      dimObject = create_dimObject(name = name,
                                   reduction_method = 'umap',
                                   coordinates = coordinates,
                                   misc = NULL)

      gobject@dimension_reduction[[reduction]][['umap']][[name]] <- dimObject



      ## update parameters used ##
      parameters_list = gobject@parameters
      number_of_rounds = length(parameters_list)
      update_name = paste0(number_of_rounds,'_umap')
      # parameters to include
      parameters_list[[update_name]] = c('reduction type' = reduction,
                                         'dimension reduction used' = dim_reduction_to_use,
                                         'name for dimension reduction' = dim_reduction_name,
                                         'dimensions used' = length(dimensions_to_use),
                                         'expression values' = expression_values,
                                         'number of genes used' = length(genes_to_use),
                                         'n_neighbors' = n_neighbors,
                                         'n_components' = n_components,
                                         'n_epochs' = n_epochs,
                                         'min_dist' = min_dist,
                                         'spread' = spread,
                                         'name for umap' = name)
      gobject@parameters = parameters_list


      return(gobject)


    } else {
      return(uwot_clus_pos_DT)
    }




  } else if(reduction == 'genes') {

    cat('To be completed')

  }

}



#' @title runtSNE
#' @name runtSNE
#' @description run tSNE
#' @param gobject giotto object
#' @param expression_values expression values to use
#' @param reduction cells or genes
#' @param dim_reduction_to_use use another dimension reduction set as input
#' @param dim_reduction_name name of dimension reduction set to use
#' @param dimensions_to_use number of dimensions to use as input
#' @param name arbitrary name for tSNE run
#' @param genes_to_use if dim_reduction_to_use = NULL, which genes to use
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @param dims tSNE param: number of dimensions to return
#' @param perplexity tSNE param: perplexity
#' @param theta tSNE param: theta
#' @param do_PCA_first tSNE param: do PCA before tSNE (default = FALSE)
#' @param set_seed use of seed
#' @param seed_number seed number to use
#' @param ... additional tSNE parameters
#' @return giotto object with updated tSNE dimension recuction
#' @details Description of tSNE steps and params ...
#' @export
#' @examples
#'     runtSNE(gobject)
runtSNE <- function(gobject,
                    expression_values = c('normalized', 'scaled', 'custom'),
                    reduction = c('cells', 'genes'),
                    dim_reduction_to_use = 'pca',
                    dim_reduction_name = 'pca',
                    dimensions_to_use = 1:10,
                    name = 'tsne',
                    genes_to_use = NULL,
                    return_gobject = TRUE,
                    dims = 2,
                    perplexity = 30,
                    theta = 0.5,
                    do_PCA_first = F,
                    set_seed = T,
                    seed_number = 1234,
                    ...) {

  reduction = match.arg(reduction, choices = c('cells', 'genes'))

  ## umap on cells ##
  if(reduction == 'cells') {

    ## using dimension reduction ##
    if(!is.null(dim_reduction_to_use)) {

      ## TODO: check if reduction exists

      matrix_to_use = gobject@dimension_reduction[['cells']][[dim_reduction_to_use]][[dim_reduction_name]][['coordinates']][, dimensions_to_use]


    } else {
      ## using original matrix ##
      values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))
      expr_values = select_expression_values(gobject = gobject, values = values)

      # subset expression matrix
      if(!is.null(genes_to_use)) {
        expr_values = expr_values[rownames(expr_values) %in% genes_to_use, ]
      }

      matrix_to_use = t(expr_values)

    }


    if(set_seed == TRUE) {
      set.seed(seed = seed_number)
    }

    ## run umap ##
    tsne_clus = Rtsne::Rtsne(X = matrix_to_use,
                             dims = dims,
                             perplexity = perplexity,
                             theta = theta,
                             pca = do_PCA_first, ...)

    tsne_clus_pos_DT <- as.data.table(tsne_clus$Y)
    tsne_clus_pos_DT[, cell_ID := rownames(matrix_to_use)]

    if(set_seed == TRUE) {
      set.seed(Sys.time())
    }


    if(return_gobject == TRUE) {

      tsne_names = names(gobject@dimension_reduction[[reduction]][['tsne']])

      if(name %in% tsne_names) {
        cat('\n ', name, ' has already been used, will be overwritten \n')

      }


      coordinates = tsne_clus$Y
      rownames(coordinates) = rownames(matrix_to_use)

      dimObject = create_dimObject(name = name,
                                   reduction_method = 'tsne',
                                   coordinates = coordinates,
                                   misc = tsne_clus)

      gobject@dimension_reduction[[reduction]][['tsne']][[name]] <- dimObject



      ## update parameters used ##
      parameters_list = gobject@parameters
      number_of_rounds = length(parameters_list)
      update_name = paste0(number_of_rounds,'_tsne')
      # parameters to include
      parameters_list[[update_name]] = c('reduction type' = reduction,
                                         'dimension reduction used' = dim_reduction_to_use,
                                         'name for dimension reduction' = dim_reduction_name,
                                         'dimensions used' = length(dimensions_to_use),
                                         'expression values' = expression_values,
                                         'number of genes used' = length(genes_to_use),
                                         'output dimensions' = dims,
                                         'perplexity' = perplexity,
                                         'theta' = theta,
                                         'perform PCA first' = do_PCA_first,
                                         'name for tsne' = name)
      gobject@parameters = parameters_list

      return(gobject)


    } else {
      return(tsne_clus_pos_DT)
    }




  } else if(reduction == 'genes') {

    cat('To be completed')

  }

}

