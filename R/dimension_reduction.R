



#' @title create_dimObject
#' @name create_dimObject
#' @description Creates an object that stores a dimension reduction output
#' @param name arbitrary name for object
#' @param reduction_method method used to reduce dimensions
#' @param coordinates accepts the coordinates after dimension reduction
#' @param misc any additional information will be added to this slot
#' @param my_rownames rownames
#' @keywords internal
#' @return number of distinct colors
create_dimObject = function(name = 'test',
                            reduction_method = NULL,
                            coordinates = NULL,
                            misc = NULL,
                            my_rownames = NULL) {


  number_of_dimensions = ncol(coordinates)
  colnames(coordinates) = paste0('Dim.',1:number_of_dimensions)

  if(!is.null(my_rownames)) {
    rownames(coordinates) = my_rownames
  }

  dimObj = list(name = name,
                reduction_method = reduction_method,
                coordinates = coordinates,
                misc = misc)

  class(dimObj) <- append(class(dimObj), 'dimObj')
  return(dimObj)

}



#' @title select_dimReduction
#' @name select_dimReduction
#' @description Creates an object that stores a dimension reduction output
#' @keywords internal
#' @return dim reduction coordinates (default) or dim reduction object
select_dimReduction = function(gobject,
                               reduction = c('cells', 'genes'),
                               reduction_method = c('pca', 'umap', 'tsne'),
                               name = 'pca',
                               return_dimObj = FALSE) {


  ## check parameters
  reduction = match.arg(arg = reduction, choices = c('cells', 'genes'))
  reduction_method = match.arg(arg = reduction_method, choices = c('pca', 'umap', 'tsne'))

  ## check reduction
  reduction_res = gobject@dimension_reduction[[reduction]]
  if(is.null(reduction_res)) {
    stop('No dimension reduction for ', reduction, ' has been applied \n')
  }

  ## check method
  reduction_res = reduction_res[[reduction_method]]
  if(is.null(reduction_res)) {
    stop(reduction_method, ' has not been performed on this dataset \n')
  }

  ## check name for method
  reduction_res = reduction_res[[name]]
  if(is.null(reduction_res)) {
    stop(name, ': this name is not available for method: ', reduction_method, '\n')
  }

  ## return object or coordinates
  if(return_dimObj == TRUE) {
    return(reduction_res)
  } else {
    return(reduction_res$coordinates)
  }

}




#' @title standardise_giotto
#' @name standardise_giotto
#' @description standardises a matrix
#' @param x matrix
#' @param center center data
#' @param scale scale data
#' @keywords internal
#' @return standardized matrix
standardise_giotto = function (x, center = TRUE, scale = TRUE)
{
  if (center & scale) {
    y <- t_giotto(x) - Rfast::colmeans(x)
    y <- y/sqrt(Rfast::rowsums(y^2)) * sqrt((dim(x)[1] -
                                               1))
    y <- t_giotto(y)
  }
  else if (center & !scale) {
    m <- Rfast::colmeans(x)
    y <- Rfast::eachrow(x, m, oper = "-")
  }
  else if (!center & scale) {
    s <- Rfast::colVars(x, std = TRUE)
    y <- Rfast::eachrow(x, s, oper = "/")
  } else {
    y = x
  }
}

#' @title pca_giotto
#' @name pca_giotto
#' @description performs PCA based on Rfast
#' @param mymatrix matrix or object that can be converted to matrix
#' @param center center data
#' @param scale scale features
#' @param k number of principal components to calculate
#' @keywords internal
#' @return list of eigenvalues, eigenvectors and pca coordinates
pca_giotto = function(mymatrix, center = T, scale = T, k = 50) {

  if(!is.null(k) & k > ncol(mymatrix)) {
    warning('k > ncol(matrix), will be set to ncol(matrix)')
    k = ncol(mymatrix)
  }

  if(!is.matrix(mymatrix)) mymatrix = as.matrix(mymatrix)
  my_t_matrix = t_giotto(mymatrix)
  pca_f = Rfast::hd.eigen(x = my_t_matrix, center = center, scale = scale, k = k, vectors = TRUE)

  # calculate pca coordinates
  rotated_mat = standardise_giotto(x = my_t_matrix, center = center, scale = scale)
  coords = rotated_mat %*% pca_f$vectors
  colnames(coords) = paste0('Dim.', 1:ncol(coords))

  return(list(eigenvalues = pca_f$values, eigenvectors = pca_f$vectors, coords = coords))

}


#' @title runPCA_prcomp_irlba
#' @name runPCA_prcomp_irlba
#' @description performs PCA based on the irlba package
#' @param x matrix or object that can be converted to matrix
#' @param ncp number of principal components to calculate
#' @param center center data
#' @param scale scale features
#' @param rev reverse PCA
#' @param set_seed use of seed
#' @param seed_number seed number to use
#' @keywords internal
#' @return list of eigenvalues, loadings and pca coordinates
runPCA_prcomp_irlba = function(x,
                               ncp = 100,
                               center = TRUE,
                               scale = TRUE,
                               rev = FALSE,
                               set_seed = TRUE,
                               seed_number = 1234,
                               ...) {

  min_ncp = min(dim(x))

  if(ncp >= min_ncp) {
    warning("ncp >= minimum dimension of x, will be set to minimum dimension of x - 1")
    ncp = min_ncp-1
  }

  if(rev == TRUE) {

    x = t_giotto(x)

    if(set_seed == TRUE) {
      set.seed(seed = seed_number)
    }
    pca_res = irlba::prcomp_irlba(x = x, n = ncp, center = center, scale. = scale, ...)
    # eigenvalues
    eigenvalues = pca_res$sdev^2
    # PC loading
    loadings = pca_res$x
    rownames(loadings) = rownames(x)
    colnames(loadings) = paste0('Dim.', 1:ncol(loadings))
    # coordinates
    coords = pca_res$rotation
    rownames(coords) = colnames(x)
    colnames(coords) = paste0('Dim.', 1:ncol(coords))
    result = list(eigenvalues = eigenvalues, loadings = loadings, coords = coords)

  } else {

    if(set_seed == TRUE) {
      set.seed(seed = seed_number)
    }
    pca_res = irlba::prcomp_irlba(x = x, n = ncp, center = center, scale. = scale, ...)
    # eigenvalues
    eigenvalues = pca_res$sdev^2
    # PC loading
    loadings = pca_res$rotation
    rownames(loadings) = colnames(x)
    colnames(loadings) = paste0('Dim.', 1:ncol(loadings))
    # coordinates
    coords = pca_res$x
    rownames(coords) = rownames(x)
    colnames(coords) = paste0('Dim.', 1:ncol(coords))
    result = list(eigenvalues = eigenvalues, loadings = loadings, coords = coords)

  }

  return(result)

}



#' @title runPCA_factominer
#' @name runPCA_factominer
#' @description performs PCA based on the factominer package
#' @param x matrix or object that can be converted to matrix
#' @param ncp number of principal components to calculate
#' @param scale scale features
#' @param rev reverse PCA
#' @param set_seed use of seed
#' @param seed_number seed number to use
#' @keywords internal
#' @return list of eigenvalues, loadings and pca coordinates
runPCA_factominer = function(x,
                             ncp = 100,
                             scale = TRUE,
                             rev = FALSE,
                             set_seed = TRUE,
                             seed_number = 1234,
                             ...) {

  if(!methods::is(x, 'matrix')) {
    x = as.matrix(x)
  }

  if(rev == TRUE) {

    x = t_giotto(x)

    if(ncp > nrow(x)) {
      warning("ncp > nrow(x), will be set to nrow(x)")
      ncp = nrow(x)
    }

    if(set_seed == TRUE) {
      set.seed(seed = seed_number)
    }
    pca_res = FactoMineR::PCA(X = x, ncp = ncp, scale.unit = scale, graph = F, ...)

    # eigenvalues
    eigenvalues = pca_res$eig[,1]

    # PC loading
    loadings = pca_res$ind$coord
    rownames(loadings) = rownames(x)
    colnames(loadings) = paste0('Dim.', 1:ncol(loadings))

    # coordinates
    coords = sweep(pca_res$var$coord, 2, sqrt(eigenvalues[1:ncp]), FUN = "/")
    rownames(coords) = colnames(x)
    colnames(coords) = paste0('Dim.', 1:ncol(coords))

    result = list(eigenvalues = eigenvalues, loadings = loadings, coords = coords)

  } else {

    if(ncp > ncol(x)) {
      warning("ncp > ncol(x), will be set to ncol(x)")
      ncp = ncol(x)
    }

    if(set_seed == TRUE) {
      set.seed(seed = seed_number)
    }
    pca_res = FactoMineR::PCA(X = x, ncp = ncp, scale.unit = scale, graph = F, ...)

    # eigenvalues
    eigenvalues = pca_res$eig[,1]

    # PC loading
    loadings = sweep(pca_res$var$coord, 2, sqrt(eigenvalues[1:ncp]), FUN = "/")
    rownames(loadings) = colnames(x)
    colnames(loadings) = paste0('Dim.', 1:ncol(loadings))

    # coordinates
    coords = pca_res$ind$coord
    rownames(coords) = rownames(x)
    colnames(coords) = paste0('Dim.', 1:ncol(coords))

    result = list(eigenvalues = eigenvalues, loadings = loadings, coords = coords)

  }

  return(result)

}



#' @title create_genes_to_use_matrix
#' @name create_genes_to_use_matrix
#' @description subsets matrix based on vector of genes or hvg column
#' @param gobject giotto object
#' @param sel_matrix selected expression matrix
#' @param genes_to_use genes to use, character or vector of genes
#' @param verbose verbosity
#' @keywords internal
#' @return subsetted matrix based on selected genes
create_genes_to_use_matrix = function(gobject,
                                      sel_matrix,
                                      genes_to_use,
                                      verbose = TRUE) {

  # cell metadata
  gene_metadata = fDataDT(gobject)

  # for hvg genes
  if(is.character(genes_to_use) & length(genes_to_use) == 1) {
    if(genes_to_use %in% colnames(gene_metadata)) {
      if(verbose == TRUE) cat(genes_to_use, ' was found in the gene metadata information and will be used to select highly variable genes \n')
      genes_to_use = gene_metadata[get(genes_to_use) == 'yes'][['gene_ID']]
      sel_matrix = sel_matrix[rownames(sel_matrix) %in% genes_to_use, ]
    } else {
      if(verbose == TRUE) cat(genes_to_use, ' was not found in the gene metadata information, all genes will be used \n')
    }
  } else {
    if(verbose == TRUE) cat('a custom vector of genes will be used to subset the matrix \n')
    sel_matrix = sel_matrix[rownames(sel_matrix) %in% genes_to_use, ]
  }

  return(sel_matrix)

}


#' @title runPCA
#' @name runPCA
#' @description runs a Principal Component Analysis
#' @param gobject giotto object
#' @param expression_values expression values to use
#' @param reduction cells or genes
#' @param name arbitrary name for PCA run
#' @param genes_to_use subset of genes to use for PCA
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @param center center data first (default = TRUE)
#' @param scale_unit scale features before PCA (default = TRUE)
#' @param ncp number of principal components to calculate
#' @param method which implementation to use
#' @param rev do a reverse PCA
#' @param set_seed use of seed
#' @param seed_number seed number to use
#' @param verbose verbosity of the function
#' @param ... additional parameters for PCA (see details)
#' @return giotto object with updated PCA dimension recuction
#' @details See \code{\link[irlba]{prcomp_irlba}} and \code{\link[FactoMineR]{PCA}} for more information about other parameters.
#' \itemize{
#'   \item genes_to_use = NULL: will use all genes from the selected matrix
#'   \item genes_to_use = <hvg name>: can be used to select a column name of
#'   highly variable genes, created by (see \code{\link{calculateHVG}})
#'   \item genes_to_use = c('geneA', 'geneB', ...): will use all manually provided genes
#' }
#' @export
#' @examples
#'
#' data(mini_giotto_single_cell)
#'
#' # run PCA
#' mini_giotto_single_cell <- runPCA(gobject = mini_giotto_single_cell,
#'                                   center = TRUE, scale_unit = TRUE)
#'
#' # plot PCA results
#' plotPCA(mini_giotto_single_cell)
#'
runPCA <- function(gobject,
                   expression_values = c('normalized', 'scaled', 'custom'),
                   reduction = c('cells', 'genes'),
                   name = 'pca',
                   genes_to_use = 'hvg',
                   return_gobject = TRUE,
                   center = TRUE,
                   scale_unit = TRUE,
                   ncp = 100,
                   method = c('irlba','factominer'),
                   rev = FALSE,
                   set_seed = TRUE,
                   seed_number = 1234,
                   verbose = TRUE,
                   ...) {


  # expression values to be used
  values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))
  expr_values = select_expression_values(gobject = gobject, values = values)

  ## subset matrix
  if(!is.null(genes_to_use)) {
    expr_values = create_genes_to_use_matrix(gobject = gobject,
                                                      sel_matrix = expr_values,
                                                      genes_to_use = genes_to_use,
                                                      verbose = verbose)
  }


  # do PCA dimension reduction
  reduction = match.arg(reduction, c('cells', 'genes'))

  # PCA implementation
  method = match.arg(method, c('irlba','factominer'))

  if(reduction == 'cells') {
    # PCA on cells
    if(method == 'irlba') {
      pca_object = runPCA_prcomp_irlba(x = t_giotto(expr_values),
                                       center = center, scale = scale_unit, ncp = ncp,
                                       rev = rev, set_seed = set_seed, seed_number = seed_number, ...)
    } else if(method == 'factominer') {
      pca_object = runPCA_factominer(x = t_giotto(expr_values),
                                     scale = scale_unit, ncp = ncp, rev = rev,
                                     set_seed = set_seed, seed_number = seed_number, ...)
    } else {
      stop('only PCA methods from the irlba and factominer package have been implemented \n')
    }

  } else {
    # PCA on genes
    if(method == 'irlba') {
      pca_object = runPCA_prcomp_irlba(x = expr_values,
                                       center = center, scale = scale_unit, ncp = ncp,
                                       rev = rev, set_seed = set_seed, seed_number = seed_number, ...)
    } else if(method == 'factominer') {
      pca_object = runPCA_factominer(x = expr_values,
                                     scale = scale_unit, ncp = ncp, rev = rev,
                                     set_seed = set_seed, seed_number = seed_number, ...)
    } else {
      stop('only PCA methods from the irlba and factominer package have been implemented \n')
    }

  }


  if(return_gobject == TRUE) {

    pca_names = names(gobject@dimension_reduction[[reduction]][['pca']])

    if(name %in% pca_names) {
      cat('\n ', name, ' has already been used, will be overwritten \n')

    }

    dimObject = create_dimObject(name = name,
                                          reduction_method = 'pca',
                                          coordinates = pca_object$coords,
                                          misc = list(eigenvalues = pca_object$eigenvalues,
                                                      loadings = pca_object$loadings),
                                          my_rownames = colnames(expr_values))

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
                                       'package' = method,
                                       'center' = center,
                                       'scale_unit' = scale_unit,
                                       'name for pca' = name)
    gobject@parameters = parameters_list

    return(gobject)


  } else {
    return(pca_object)
  }
}



#' @title create_screeplot
#' @name create_screeplot
#' @description create screeplot with ggplot
#' @param pca_obj pca dimension reduction object
#' @param ncp number of principal components to calculate
#' @param ylim y-axis limits on scree plot
#' @return ggplot
#' @keywords internal
create_screeplot = function(pca_obj, ncp = 20, ylim = c(0, 20)) {


  # data.table: set global variable
  PC = NULL

  eigs = pca_obj$misc$eigenvalues

  # variance explained
  var_expl = eigs/sum(eigs)*100
  var_expl_cum = cumsum(eigs)/sum(eigs)*100

  # create data.table
  screeDT = data.table::data.table('PC' = paste0('PC.', 1:length(var_expl)),
                                   'var_expl' = var_expl,
                                   'var_expl_cum' = var_expl_cum)
  screeDT[, PC := factor(PC, levels = PC)]

  max_ncp = length(eigs)
  ncp = ifelse(ncp > max_ncp, max_ncp, ncp)

  pl = ggplot2::ggplot()
  pl = pl + ggplot2::theme_bw()
  pl = pl + ggplot2::geom_bar(data = screeDT[1:ncp], ggplot2::aes(x = PC, y = var_expl), stat = 'identity')
  pl = pl + ggplot2::coord_cartesian(ylim = ylim)
  pl = pl + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1))
  pl = pl + ggplot2::labs(x = '', y = '% of variance explained per PC')

  cpl = ggplot2::ggplot()
  cpl = cpl + ggplot2::theme_bw()
  cpl = cpl + ggplot2::geom_bar(data = screeDT[1:ncp], ggplot2::aes(x = PC, y = var_expl_cum), stat = 'identity')
  cpl = cpl + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1))
  cpl = cpl + ggplot2::labs(x = '', y = 'cumulative % of variance explained')

  savelist = list(pl, cpl)

  ## combine plots with cowplot
  combo_plot <- cowplot::plot_grid(plotlist = savelist,
                                   ncol = 1,
                                   rel_heights = c(1),
                                   rel_widths = c(1),
                                   align = 'v')

  return(combo_plot)

}



#' @title screePlot
#' @name screePlot
#' @description identify significant prinicipal components (PCs) using an screeplot (a.k.a. elbowplot)
#' @param gobject giotto object
#' @param name name of PCA object if available
#' @param expression_values expression values to use
#' @param reduction cells or genes
#' @param method which implementation to use
#' @param rev do a reverse PCA
#' @param genes_to_use subset of genes to use for PCA
#' @param center center data before PCA
#' @param scale_unit scale features before PCA
#' @param ncp number of principal components to calculate
#' @param ylim y-axis limits on scree plot
#' @param verbose verobsity
#' @param show_plot show plot
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters from all_plots_save_function()
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @param ... additional arguments to pca function, see \code{\link{runPCA}}
#' @return ggplot object for scree method
#' @details
#'  Screeplot works by plotting the explained variance of each
#'  individual PC in a barplot allowing you to identify which PC provides a significant
#'  contribution (a.k.a 'elbow method'). \cr
#'  Screeplot will use an available pca object, based on the parameter 'name', or it will
#'  create it if it's not available (see \code{\link{runPCA}})
#' @export
#' @examples
#'
#' data(mini_giotto_single_cell)
#'
#' screePlot(mini_giotto_single_cell, ncp = 10)
#'
screePlot = function(gobject,
                     name = 'pca',
                     expression_values = c('normalized', 'scaled', 'custom'),
                     reduction = c('cells', 'genes'),
                     method = c('irlba','factominer'),
                     rev = FALSE,
                     genes_to_use = NULL,
                     center = F,
                     scale_unit = F,
                     ncp = 100,
                     ylim = c(0, 20),
                     verbose = T,
                     show_plot = NA,
                     return_plot = NA,
                     save_plot = NA,
                     save_param = list(),
                     default_save_name = 'screePlot',
                     ...) {


  # select direction of reduction
  reduction = match.arg(reduction, c('cells', 'genes'))
  pca_obj = gobject@dimension_reduction[[reduction]]$pca[[name]]

  # print, return and save parameters
  show_plot = ifelse(is.na(show_plot), readGiottoInstructions(gobject, param = 'show_plot'), show_plot)
  save_plot = ifelse(is.na(save_plot), readGiottoInstructions(gobject, param = 'save_plot'), save_plot)
  return_plot = ifelse(is.na(return_plot), readGiottoInstructions(gobject, param = 'return_plot'), return_plot)


  # if pca already exists plot
  if(!is.null(pca_obj)) {
    if(verbose == TRUE) cat('PCA with name: ', name, ' already exists and will be used for the screeplot \n')

    screeplot = create_screeplot(pca_obj = pca_obj, ncp = ncp, ylim = ylim)

  } else {
    # if pca doesn't exists, then create pca and then plot
    if(verbose == TRUE) cat('PCA with name: ', name, ' does NOT exists, PCA will be done first \n')

    # expression values to be used
    values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))
    expr_values = select_expression_values(gobject = gobject, values = values)

    # PCA implementation
    method = match.arg(method, c('irlba','factominer'))

    ## subset matrix
    if(!is.null(genes_to_use)) {
      expr_values = create_genes_to_use_matrix(gobject = gobject,
                                               sel_matrix = expr_values,
                                               genes_to_use = genes_to_use,
                                               verbose = verbose)
    }

    # reduction of cells
    if(reduction == 'cells') {

      # PCA on cells
      if(method == 'irlba') {
        pca_object = runPCA_prcomp_irlba(x = t_giotto(expr_values), center = center, scale = scale_unit, ncp = ncp, rev = rev, ...)
      } else if(method == 'factominer') {
        pca_object = runPCA_factominer(x = t_giotto(expr_values), scale = scale_unit, ncp = ncp, rev = rev, ...)
      } else {
        stop('only PCA methods from the irlba and factominer package have been implemented \n')
      }

      dimObject = create_dimObject(name = name,
                                            reduction_method = 'pca',
                                            coordinates = pca_object$coords,
                                            misc = list(eigenvalues = pca_object$eigenvalues,
                                                        loadings = pca_object$loadings),
                                            my_rownames = colnames(expr_values))

      screeplot = create_screeplot(pca_obj = dimObject, ncp = ncp, ylim = ylim)
    }

  }

  ## print plot
  if(show_plot == TRUE) {
    print(screeplot)
  }

  ## save plot
  if(save_plot == TRUE) {
    do.call('all_plots_save_function', c(list(gobject = gobject, plot_object = screeplot, default_save_name = default_save_name), save_param))
  }

  ## return plot
  if(return_plot == TRUE) {
    return(screeplot)
  }
}


#' @title create_jackstrawplot
#' @name create_jackstrawplot
#' @description create jackstrawplot with ggplot
#' @param jackstraw_data result from jackstraw function
#' @param ncp number of principal components to calculate
#' @param ylim y-axis limits on jackstraw plot
#' @param threshold p.value threshold to call a PC significant
#' @keywords internal
#' @return ggplot
create_jackstrawplot = function(jackstraw_data,
                                ncp = 20,
                                ylim = c(0, 1),
                                threshold = 0.01) {

  # data.table variables
  PC = p.val = NULL

  testDT = data.table(PC = paste0('PC.', 1:length(jackstraw_data)),
                      p.val = jackstraw_data)
  testDT[, PC := factor(PC, levels = PC)]
  testDT[, sign := ifelse(p.val <= threshold, 'sign', 'n.s.')]

  pl = ggplot2::ggplot()
  pl = pl + ggplot2::theme_bw()
  pl = pl + ggplot2::geom_point(data = testDT[1:ncp], ggplot2::aes(x = PC, y = p.val, fill = sign), shape = 21)
  pl = pl + ggplot2::scale_fill_manual(values  = c('n.s.' = 'lightgrey', 'sign' = 'darkorange'))
  pl = pl + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1))
  pl = pl + ggplot2::coord_cartesian(ylim = ylim)
  pl = pl + ggplot2::theme(panel.grid.major.x = ggplot2::element_blank())
  pl = pl + ggplot2::labs(x = '', y = 'p-value per PC')

  return(pl)

}


#' @title jackstrawPlot
#' @name jackstrawPlot
#' @description identify significant prinicipal components (PCs)
#' @param gobject giotto object
#' @param expression_values expression values to use
#' @param reduction cells or genes
#' @param genes_to_use subset of genes to use for PCA
#' @param center center data before PCA
#' @param scale_unit scale features before PCA
#' @param ncp number of principal components to calculate
#' @param ylim y-axis limits on jackstraw plot
#' @param iter number of interations for jackstraw
#' @param threshold p-value threshold to call a PC significant
#' @param verbose show progress of jackstraw method
#' @param show_plot show plot
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters from all_plots_save_function()
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return ggplot object for jackstraw method
#' @details
#'  The Jackstraw method uses the \code{\link[jackstraw]{permutationPA}} function. By
#'  systematically permuting genes it identifies robust, and thus significant, PCs.
#'  \cr
#' @export
#' @examples
#'
#' \donttest{
#'
#' data(mini_giotto_single_cell)
#'
#' # jackstraw package is required to run
#' jackstrawPlot(mini_giotto_single_cell, ncp = 10)
#'
#' }
#'
jackstrawPlot = function(gobject,
                         expression_values = c('normalized', 'scaled', 'custom'),
                         reduction = c('cells', 'genes'),
                         genes_to_use = NULL,
                         center = FALSE,
                         scale_unit = FALSE,
                         ncp = 20,
                         ylim = c(0, 1),
                         iter = 10,
                         threshold = 0.01,
                         verbose = TRUE,
                         show_plot = NA,
                         return_plot = NA,
                         save_plot = NA,
                         save_param = list(),
                         default_save_name = 'jackstrawPlot') {


  package_check(pkg_name = "jackstraw", repository = "CRAN")

  # print message with information #
  if(verbose) message("using 'jackstraw' to identify significant PCs If used in published research, please cite:
  Neo Christopher Chung and John D. Storey (2014).
  'Statistical significance of variables driving systematic variation in high-dimensional data. Bioinformatics")

  # select direction of reduction
  reduction = match.arg(reduction, c('cells', 'genes'))

  # print, return and save parameters
  show_plot = ifelse(is.na(show_plot), readGiottoInstructions(gobject, param = 'show_plot'), show_plot)
  save_plot = ifelse(is.na(save_plot), readGiottoInstructions(gobject, param = 'save_plot'), save_plot)
  return_plot = ifelse(is.na(return_plot), readGiottoInstructions(gobject, param = 'return_plot'), return_plot)

  # expression values to be used
  values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))
  expr_values = select_expression_values(gobject = gobject, values = values)

  ## subset matrix
  if(!is.null(genes_to_use)) {
    expr_values = create_genes_to_use_matrix(gobject = gobject,
                                             sel_matrix = expr_values,
                                             genes_to_use = genes_to_use,
                                             verbose = verbose)
  }

  # reduction of cells
  if(reduction == 'cells') {

    if(scale_unit == TRUE | center == TRUE) {
      expr_values = t_giotto(scale(t_giotto(expr_values), center = center, scale = scale_unit))
    }

    jtest = jackstraw::permutationPA(dat = expr_values, B = iter, threshold = threshold, verbose = verbose)

    ## results ##
    nr_sign_components = jtest$r
    cat('number of estimated significant components: ', nr_sign_components, '\n')
    final_results = jtest$p
    jackplot = create_jackstrawplot(jackstraw_data = final_results, ncp = ncp, ylim = ylim, threshold = threshold)

  }

  ## print plot
  if(show_plot == TRUE) {
    print(jackplot)
  }

  ## save plot
  if(save_plot == TRUE) {
    do.call('all_plots_save_function', c(list(gobject = gobject, plot_object = jackplot, default_save_name = default_save_name), save_param))
  }

  ## return plot
  if(return_plot == TRUE) {
    return(jackplot)
  }

}



#' @title signPCA
#' @name signPCA
#' @description identify significant prinicipal components (PCs)
#' @param gobject giotto object
#' @param name name of PCA object if available
#' @param method method to use to identify significant PCs
#' @param expression_values expression values to use
#' @param reduction cells or genes
#' @param pca_method which implementation to use
#' @param rev do a reverse PCA
#' @param genes_to_use subset of genes to use for PCA
#' @param center center data before PCA
#' @param scale_unit scale features before PCA
#' @param ncp number of principal components to calculate
#' @param scree_ylim y-axis limits on scree plot
#' @param jack_iter number of interations for jackstraw
#' @param jack_threshold p-value threshold to call a PC significant
#' @param jack_ylim y-axis limits on jackstraw plot
#' @param verbose verbosity
#' @param show_plot show plot
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters from all_plots_save_function()
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return ggplot object for scree method and maxtrix of p-values for jackstraw
#' @details Two different methods can be used to assess the number of relevant or significant
#'  prinicipal components (PC's). \cr
#'  1. Screeplot works by plotting the explained variance of each
#'  individual PC in a barplot allowing you to identify which PC provides a significant
#'  contribution  (a.k.a. 'elbow method'). \cr
#'  2. The Jackstraw method uses the \code{\link[jackstraw]{permutationPA}} function. By
#'  systematically permuting genes it identifies robust, and thus significant, PCs.
#'  \cr
#' @export
signPCA <- function(gobject,
                    name = 'pca',
                    method = c('screeplot', 'jackstraw'),
                    expression_values = c("normalized", "scaled", "custom"),
                    reduction = c("cells", "genes"),
                    pca_method = c('irlba', 'factominer'),
                    rev = FALSE,
                    genes_to_use = NULL,
                    center = T,
                    scale_unit = T,
                    ncp = 50,
                    scree_ylim = c(0,10),
                    jack_iter = 10,
                    jack_threshold = 0.01,
                    jack_ylim = c(0, 1),
                    verbose = TRUE,
                    show_plot = NA,
                    return_plot = NA,
                    save_plot = NA,
                    save_param = list(),
                    default_save_name = 'signPCA') {

  # select method
  method = match.arg(method, choices = c('screeplot', 'jackstraw'))

  # select direction of reduction
  reduction = match.arg(reduction, c('cells', 'genes'))

  # expression values to be used
  values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))
  expr_values = select_expression_values(gobject = gobject, values = values)

  # print, return and save parameters
  show_plot = ifelse(is.na(show_plot), readGiottoInstructions(gobject, param = 'show_plot'), show_plot)
  save_plot = ifelse(is.na(save_plot), readGiottoInstructions(gobject, param = 'save_plot'), save_plot)
  return_plot = ifelse(is.na(return_plot), readGiottoInstructions(gobject, param = 'return_plot'), return_plot)


  ## subset matrix
  if(!is.null(genes_to_use)) {
    expr_values = create_genes_to_use_matrix(gobject = gobject,
                                             sel_matrix = expr_values,
                                             genes_to_use = genes_to_use,
                                             verbose = verbose)
  }

  # reduction of cells
  if(reduction == 'cells') {

    if(method == 'screeplot') {

      screeplot = screePlot(gobject = gobject,
                            name = name,
                            expression_values = values,
                            reduction = reduction,
                            genes_to_use = genes_to_use,
                            center = center,
                            scale_unit = scale_unit,
                            ncp = ncp,
                            rev = rev,
                            method = pca_method,
                            ylim = scree_ylim,
                            verbose = verbose,
                            show_plot = FALSE,
                            return_plot = TRUE,
                            save_plot = FALSE,
                            save_param = list(),
                            default_save_name = 'screePlot')

      ## print plot
      if(show_plot == TRUE) {
        print(screeplot)
      }

      ## save plot
      if(save_plot == TRUE) {
        do.call('all_plots_save_function', c(list(gobject = gobject, plot_object = screeplot, default_save_name = default_save_name), save_param))
      }

      ## return plot
      if(return_plot == TRUE) {
        return(screeplot)
      }


    } else if(method == 'jackstraw') {


      jackplot = jackstrawPlot(gobject = gobject,
                               expression_values = values,
                               reduction = reduction,
                               genes_to_use = genes_to_use,
                               center = center,
                               scale_unit = scale_unit,
                               ncp = ncp,
                               ylim = jack_ylim,
                               iter = jack_iter,
                               threshold = jack_threshold,
                               verbose = verbose,
                               show_plot = FALSE,
                               return_plot = TRUE,
                               save_plot = FALSE,
                               save_param = list(),
                               default_save_name = 'jackstrawPlot')

      ## print plot
      if(show_plot == TRUE) {
        print(jackplot)
      }

      ## save plot
      if(save_plot == TRUE) {
        do.call('all_plots_save_function', c(list(gobject = gobject, plot_object = jackplot, default_save_name = default_save_name), save_param))
      }

      ## return plot
      if(return_plot == TRUE) {
        return(jackplot)
      } else {
        return(jackplot) # poentially return all results instead
      }

    }

  } else {
    cat('gene reduction not yet implemented')
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
#' @param n_threads UMAP param: threads/cores to use
#' @param spread UMAP param: spread
#' @param set_seed use of seed
#' @param seed_number seed number to use
#' @param verbose verbosity of function
#' @param ... additional UMAP parameters
#' @return giotto object with updated UMAP dimension recuction
#' @details See \code{\link[uwot]{umap}} for more information about these and other parameters.
#' \itemize{
#'   \item Input for UMAP dimension reduction can be another dimension reduction (default = 'pca')
#'   \item To use gene expression as input set dim_reduction_to_use = NULL
#'   \item If dim_reduction_to_use = NULL, genes_to_use can be used to select a column name of
#'   highly variable genes (see \code{\link{calculateHVG}}) or simply provide a vector of genes
#'   \item multiple UMAP results can be stored by changing the \emph{name} of the analysis
#' }
#' @export
#' @examples
#'
#' data(mini_giotto_single_cell)
#'
#' mini_giotto_single_cell <- runUMAP(mini_giotto_single_cell,
#'                                    dimensions_to_use = 1:3,
#'                                    n_threads = 1,
#'                                    n_neighbors = 3)
#'
#' plotUMAP(gobject = mini_giotto_single_cell)
#'
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
                    n_threads = NA,
                    spread = 5,
                    set_seed = TRUE,
                    seed_number = 1234,
                    verbose = T,
                    ...) {

  reduction = match.arg(reduction, choices = c('cells', 'genes'))

  # set cores to use
  n_threads = determine_cores(cores = n_threads)

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

      ## subset matrix
      if(!is.null(genes_to_use)) {
        expr_values = create_genes_to_use_matrix(gobject = gobject,
                                                 sel_matrix = expr_values,
                                                 genes_to_use = genes_to_use,
                                                 verbose = verbose)
      }

      matrix_to_use = t_giotto(expr_values)
    }


    ## run umap ##
    if(set_seed == TRUE) {
      set.seed(seed = seed_number)
    }

    uwot_clus <- uwot::umap(X = as.matrix(matrix_to_use), n_neighbors = n_neighbors, n_components = n_components,
                            n_epochs = n_epochs, min_dist = min_dist, n_threads = n_threads, spread = spread, ...)
    uwot_clus_pos_DT <- data.table::as.data.table(uwot_clus)

    # data.table variables
    cell_ID = NULL
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

    cat('\n Gene reduction is not yet implemented \n')

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
#' @param verbose verbosity of the function
#' @param ... additional tSNE parameters
#' @return giotto object with updated tSNE dimension recuction
#' @details See \code{\link[Rtsne]{Rtsne}} for more information about these and other parameters. \cr
#' \itemize{
#'   \item Input for tSNE dimension reduction can be another dimension reduction (default = 'pca')
#'   \item To use gene expression as input set dim_reduction_to_use = NULL
#'   \item If dim_reduction_to_use = NULL, genes_to_use can be used to select a column name of
#'   highly variable genes (see \code{\link{calculateHVG}}) or simply provide a vector of genes
#'   \item multiple tSNE results can be stored by changing the \emph{name} of the analysis
#' }
#' @export
#' @examples
#'
#' data(mini_giotto_single_cell)
#'
#' mini_giotto_single_cell <- runtSNE(mini_giotto_single_cell,
#'                                    dimensions_to_use = 1:3,
#'                                    n_threads = 1,
#'                                    n_neighbors = 3,
#'                                    perplexity = 1)
#'
#' plotTSNE(gobject = mini_giotto_single_cell)
#'
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
                    verbose = TRUE,
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

      ## subset matrix
      if(!is.null(genes_to_use)) {
        expr_values = create_genes_to_use_matrix(gobject = gobject,
                                                 sel_matrix = expr_values,
                                                 genes_to_use = genes_to_use,
                                                 verbose = verbose)
      }

      matrix_to_use = t_giotto(expr_values)

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

    tsne_clus_pos_DT <- data.table::as.data.table(tsne_clus$Y)

    # data.table variables
    cell_ID = NULL
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

    cat('\n Not yet implemented \n')

  }

}

