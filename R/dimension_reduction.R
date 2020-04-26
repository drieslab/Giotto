



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


#' @title standardise_giotto
#' @name standardise_giotto
#' @description standardises a matrix
#' @param x matrix 
#' @param center center data
#' @param scale scale data
#' @return standardized matrix
standardise_giotto = function (x, center = TRUE, scale = TRUE) 
{
  if (center & scale) {
    y <- t(x) - Rfast::colmeans(x)
    y <- y/sqrt(Rfast::rowsums(y^2)) * sqrt((dim(x)[1] - 
                                               1))
    y <- t(y)
  }
  else if (center & !scale) {
    m <- Rfast::colmeans(x)
    y <- eachrow(x, m, oper = "-")
  }
  else if (!center & scale) {
    s <- Rfast::colVars(x, std = TRUE)
    y <- eachrow(x, s, oper = "/")
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
#' @return list of eigenvalues, eigenvectors and pca coordinates
pca_giotto = function(mymatrix, center = T, scale = T, k = 50) {
  
  if(!is.null(k) & k > ncol(mymatrix)) {
    warning('k > ncol(matrix), will be set to ncol(matrix)')
    k = ncol(mymatrix)
  }
  
  if(!is.matrix(mymatrix)) mymatrix = as.matrix(mymatrix)
  my_t_matrix = t(mymatrix)
  pca_f = Rfast::hd.eigen(x = my_t_matrix, center = center, scale = scale, k = k, vectors = TRUE)
  
  # calculate pca coordinates
  rotated_mat = standardise_giotto(x = my_t_matrix, center = center, scale = scale)
  coords = rotated_mat %*% pca_f$vectors
  colnames(coords) = paste0('Dim.', 1:ncol(coords))
  
  return(list(eigenvalues = pca_f$values, eigenvectors = pca_f$vectors, coords = coords))
  
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
#' @param center center data first (default = FALSE)
#' @param scale_unit scale features before PCA
#' @param ncp number of principal components to calculate
#' @param ... additional parameters for PCA (see details)
#' @return giotto object with updated PCA dimension recuction
#' @details See \code{\link[FactoMineR]{PCA}} for more information about other parameters.
#' @export
#' @examples
#'     runPCA(gobject)
runPCA <- function(gobject,
                   expression_values = c('normalized', 'scaled', 'custom'),
                   reduction = c('cells', 'genes'),
                   name = 'pca',
                   genes_to_use = NULL,
                   return_gobject = TRUE,
                   center = F,
                   scale_unit = F,
                   ncp = 100,
                   ...) {


  # expression values to be used
  values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))
  expr_values = Giotto:::select_expression_values(gobject = gobject, values = values)


  # subset expression matrix
  if(!is.null(genes_to_use)) {
    expr_values = expr_values[rownames(expr_values) %in% genes_to_use, ]
  }

  # do PCA dimension reduction
  reduction = match.arg(reduction, c('cells', 'genes'))
  if(reduction == 'cells') {
    pca_object = pca_giotto(mymatrix = expr_values, center = center, scale = scale_unit, k = ncp)
  } else {
    pca_object = pca_giotto(mymatrix = t(expr_values), center = center, scale = scale_unit, k = ncp)
  }

  if(return_gobject == TRUE) {

    pca_names = names(gobject@dimension_reduction[[reduction]][['pca']])

    if(name %in% pca_names) {
      cat('\n ', name, ' has already been used, will be overwritten \n')

    }

    dimObject = create_dimObject(name = name, reduction_method = 'pca',
                                 coordinates = pca_object$coords,
                                 misc = pca_object, my_rownames = colnames(expr_values))

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
create_screeplot = function(pca_obj, ncp = 20, ylim = c(0, 20)) {
  
  eigs = pca_obj$misc$eigenvalues
  
  # variance explained
  var_expl = eigs/sum(eigs)*100
  var_expl_cum = cumsum(eigs)/sum(eigs)*100
  
  # create data.table
  screeDT = data.table('PC' = paste0('PC.', 1:length(var_expl)),
                       'var_expl' = var_expl,
                       'var_expl_cum' = var_expl_cum)
  screeDT[, PC := factor(PC, levels = PC)]
  
  pl = ggplot()
  pl = pl + theme_bw()
  pl = pl + geom_bar(data = screeDT[1:ncp], aes(x = PC, y = var_expl), stat = 'identity')
  pl = pl + coord_cartesian(ylim = ylim)
  pl = pl + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  pl = pl + labs(x = '', y = '% of variance explained per PC')
  
  cpl = ggplot()
  cpl = cpl + theme_bw()
  cpl = cpl + geom_bar(data = screeDT[1:ncp], aes(x = PC, y = var_expl_cum), stat = 'identity')
  cpl = cpl + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  cpl = cpl + labs(x = '', y = 'cumulative % of variance explained')
  
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
#' @return ggplot object for scree method 
#' @details 
#'  Screeplot works by plotting the explained variance of each
#'  individual PC in a barplot allowing you to identify which PC does not show a significant
#'  contribution anymore ( = 'elbow method'). \cr
#' @export
#' @examples
#'     screePlot(gobject)
screePlot = function(gobject,
                     name = 'pca',
                     expression_values = c('normalized', 'scaled', 'custom'),
                     reduction = c('cells', 'genes'),
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
                     default_save_name = 'screePlot') {
  
  
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
    
    # subset expression matrix
    if(!is.null(genes_to_use)) {
      expr_values = expr_values[rownames(expr_values) %in% genes_to_use, ]
    }
    
    # reduction of cells
    if(reduction == 'cells') {
      pca_object = pca_giotto(mymatrix = expr_values, center = center, scale = scale_unit, k = ncp)
      
      dimObject = create_dimObject(name = name, reduction_method = 'pca',
                                   coordinates = pca_object$coords,
                                   misc = pca_object, my_rownames = colnames(expr_values))
      
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
#' @param p-value threshold to call a PC significant
#' @return ggplot
create_jackstrawplot = function(jackstraw_data, ncp = 20, ylim = c(0, 1), threshold = 0.01) {
  
  
  testDT = data.table(PC = paste0('PC.', 1:length(jackstraw_data)),
                      p.val = jackstraw_data)
  testDT[, PC := factor(PC, levels = PC)]
  testDT[, sign := ifelse(p.val <= threshold, 'sign', 'n.s.')]
  
  pl = ggplot()
  pl = pl + theme_bw()
  pl = pl + geom_point(data = testDT[1:ncp], aes(x = PC, y = p.val, fill = sign), shape = 21)
  pl = pl + scale_fill_manual(values  = c('n.s.' = 'lightgrey', 'sign' = 'darkorange'))
  pl = pl + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  pl = pl + coord_cartesian(ylim = ylim)
  pl = pl + theme(panel.grid.major.x = element_blank())
  pl = pl + labs(x = '', y = 'p-value per PC')
  
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
#'  \cr multiple PCA results can be stored by changing the \emph{name} parameter
#' @export
#' @examples
#'     jackstrawPlot(gobject)
jackstrawPlot = function(gobject,
                         expression_values = c('normalized', 'scaled', 'custom'),
                         reduction = c('cells', 'genes'),
                         genes_to_use = NULL,
                         center = F,
                         scale_unit = F,
                         ncp = 20,
                         ylim = c(0, 1),
                         iter = 10,
                         threshold = 0.01,
                         verbose = T,
                         show_plot = NA,
                         return_plot = NA,
                         save_plot = NA,
                         save_param = list(),
                         default_save_name = 'jackstrawPlot') {
  
  
  # select direction of reduction
  reduction = match.arg(reduction, c('cells', 'genes'))
  
  # print, return and save parameters
  show_plot = ifelse(is.na(show_plot), readGiottoInstructions(gobject, param = 'show_plot'), show_plot)
  save_plot = ifelse(is.na(save_plot), readGiottoInstructions(gobject, param = 'save_plot'), save_plot)
  return_plot = ifelse(is.na(return_plot), readGiottoInstructions(gobject, param = 'return_plot'), return_plot)
  
  # expression values to be used
  values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))
  expr_values = select_expression_values(gobject = gobject, values = values)
  
  # subset expression matrix
  if(!is.null(genes_to_use)) {
    expr_values = expr_values[rownames(expr_values) %in% genes_to_use, ]
  }
  
  # reduction of cells
  if(reduction == 'cells') {
    
    if(scale_unit == TRUE | center == TRUE) {
      expr_values = t(scale(t(expr_values), center = center, scale = scale_unit))
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
#' @param genes_to_use subset of genes to use for PCA
#' @param center center data before PCA
#' @param scale_unit scale features before PCA
#' @param ncp number of principal components to calculate
#' @param scree_ylim y-axis limits on scree plot
#' @param jack_iter number of interations for jackstraw
#' @param jack_threshold p-value threshold to call a PC significant
#' @param jack_ylim y-axis limits on jackstraw plot
#' @param jack_verbose show progress of jackstraw method
#' @param show_plot show plot
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters from all_plots_save_function()
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return ggplot object for scree method and maxtrix of p-values for jackstraw
#' @details Two different methods can be used to assess the number of relevant or significant
#'  prinicipal components (PC's). \cr
#'  1. Screeplot works by plotting the explained variance of each
#'  individual PC in a barplot allowing you to identify which PC does not show a significant
#'  contribution anymore ( = 'elbow method'). \cr
#'  2. The Jackstraw method uses the \code{\link[jackstraw]{permutationPA}} function. By
#'  systematically permuting genes it identifies robust, and thus significant, PCs.
#'  \cr multiple PCA results can be stored by changing the \emph{name} parameter
#' @export
#' @examples
#'     signPCA(gobject)
signPCA <- function(gobject,
                    name = 'pca',
                    method = c('screeplot', 'jackstraw'),
                    expression_values = c("normalized", "scaled", "custom"),
                    reduction = c("cells", "genes"),
                    genes_to_use = NULL,
                    center = T,
                    scale_unit = T,
                    ncp = 50,
                    scree_ylim = c(0,10),
                    jack_iter = 10,
                    jack_threshold = 0.01,
                    jack_ylim = c(0, 1),
                    jack_verbose = T,
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


  # subset expression matrix
  if(!is.null(genes_to_use)) {
    expr_values = expr_values[rownames(expr_values) %in% genes_to_use, ]
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
#' @param n_threads UMAP param: threads to use
#' @param spread UMAP param: spread
#' @param set_seed use of seed
#' @param seed_number seed number to use
#' @param ... additional UMAP parameters
#' @return giotto object with updated UMAP dimension recuction
#' @details See \code{\link[uwot]{umap}} for more information about these and other parameters.
#' \itemize{
#'   \item Input for UMAP dimension reduction can be another dimension reduction (default = 'pca')
#'   \item To use gene expression as input set dim_reduction_to_use = NULL
#'   \item multiple UMAP results can be stored by changing the \emph{name} of the analysis
#' }
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
                    n_threads = 1,
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
    uwot_clus_pos_DT <- data.table::as.data.table(uwot_clus)
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
#' @param ... additional tSNE parameters
#' @return giotto object with updated tSNE dimension recuction
#' @details See \code{\link[Rtsne]{Rtsne}} for more information about these and other parameters. \cr
#' \itemize{
#'   \item Input for tSNE dimension reduction can be another dimension reduction (default = 'pca')
#'   \item To use gene expression as input set dim_reduction_to_use = NULL
#'   \item multiple tSNE results can be stored by changing the \emph{name} of the analysis
#' }
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

    tsne_clus_pos_DT <- data.table::as.data.table(tsne_clus$Y)
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

