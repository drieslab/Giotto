#' @title kmeans_binarize
#' @name kmeans_binarize
#' @description create binarized scores using kmeans
kmeans_binarize = function(x) {

  sel_gene_km = stats::kmeans(x, centers = 2)$cluster
  mean_1 = mean(x[sel_gene_km == 1])
  mean_2 = mean(x[sel_gene_km == 2])

  if(mean_1 > mean_2) {
    mean_1_value = 1
    mean_2_value = 0
  } else {
    mean_1_value = 0
    mean_2_value = 1
  }

  sel_gene_bin = x
  sel_gene_bin[sel_gene_km == 1] = mean_1_value
  sel_gene_bin[sel_gene_km == 2] = mean_2_value

  return(sel_gene_bin)

}

#' @title rank_binarize
#' @name rank_binarize
#' @description create binarized scores using arbitrary rank of top genes
rank_binarize = function(x, max_rank = 200) {

  sel_gene_rank = rank(-x, ties.method = 'average')

  sel_gene_bin = x
  sel_gene_bin[sel_gene_rank <= max_rank] = 1
  sel_gene_bin[sel_gene_rank > max_rank] = 0

  return(sel_gene_bin)

}

#' @title fish_function
#' @name fish_function
#' @description perform fisher exact test
fish_function = function(x_to, x_from) {

  fish_table = table(x_to == '1',
                     x_from == '1')

  fish_res = stats::fisher.test(fish_table)

  return(list(pval = fish_res$p.value, OR = fish_res$estimate))
}

#' @title calculate_binarized_spatial_genes
#' @name calculate_binarized_spatial_genes
#' @description compute spatial genes on binarized gene values
calculate_binarized_spatial_genes <- function(gobject,
                                              expression_values = c('normalized', 'scaled', 'custom'),
                                              spatial_network_name = 'spatial_network',
                                              bin_method = c('kmeans', 'rank'),
                                              percentage_rank = 10) {

  # spatial network
  spatial_network = gobject@spatial_network[[spatial_network_name]]

  # expression
  values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))
  expr_values = select_expression_values(gobject = gobject, values = values)

  # binarization method
  bin_method = match.arg(bin_method, choices = c('kmeans', 'rank'))

  # binarize matrix
  if(bin_method == 'kmeans') {
    bin_matrix = t(apply(X = expr_values, MARGIN = 1, FUN = kmeans_binarize))
  } else if(bin_method == 'rank') {

    max_rank = (ncol(expr_values)/100)*percentage_rank
    bin_matrix = t(apply(X = expr_values, MARGIN = 1, FUN = function(x) rank_binarize(x = x, max_rank = max_rank)))
  }
  bin_matrix_DT = data.table::as.data.table(melt(bin_matrix, varnames = c('genes', 'cells'), value.name = 'value'))

  # extra info: nr of cells with high expression
  nr_high_cells = bin_matrix_DT[, .N, by = .(genes, value)][value == 1]
  nr_high_cells = nr_high_cells[,.(genes,N)]

  # combine binarized matrix with spatial network
  spatial_network_min = spatial_network[,.(to, from)]
  spatial_network_min = data.table:::merge.data.table(x = spatial_network_min, by.x = 'to', y = bin_matrix_DT, by.y = 'cells', allow.cartesian = T)
  setnames(spatial_network_min, 'value', 'to_value')
  spatial_network_min[bin_matrix_DT, from_value := value, on = c(genes = 'genes', from = 'cells')]
  setorder(spatial_network_min, 'genes')

  # perform fisher test to identify spatial genes
  spatial_results = spatial_network_min[, fish_function(x_to = to_value, x_from = from_value), by = genes]
  spatial_results = merge(spatial_results, nr_high_cells, by = 'genes')
  setorder(spatial_results, -OR)

  return(spatial_results)

}

#' @title calculateSpatialGenes
#' @name calculate_spatial_genes
#' @description compute spatial genes on gene expression products between neighbors
calculate_spatial_genes <- function(gobject,
                                    expression_values = c('normalized', 'scaled', 'custom'),
                                    spatial_network_name = 'spatial_network',
                                    method = c('gini'),
                                    simulations = 0,
                                    detection_threshold = 0) {


  values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))
  expr_values = select_expression_values(gobject = gobject, values = values)

  # calculate detection
  detection = rowSums(expr_values > detection_threshold)/ncol(expr_values)
  detection_DT = data.table::data.table(genes = names(detection), detection = detection)


  # transpose expression data
  t_expression_data = t(expr_values)
  t_expression_data_DT <- data.table::as.data.table(t_expression_data)
  t_expression_data_DT[, cell_ID := rownames(t_expression_data)]

  # spatial network
  spatial_network = gobject@spatial_network[[spatial_network_name]]

  # merge spatial network and expression data
  geneset = colnames(t_expression_data)
  spatial_netw_ext_to <- merge(spatial_network, by.x = 'to', t_expression_data_DT, by.y = 'cell_ID')
  spatial_netw_ext_to_m <- data.table::melt.data.table(spatial_netw_ext_to, measure.vars = geneset, variable.name = 'to_gene', value.name = 'to_expression')
  spatial_netw_ext_from <- merge(spatial_network, by.x = 'from', t_expression_data_DT, by.y = 'cell_ID')
  spatial_netw_ext_from_m <- data.table::melt.data.table(spatial_netw_ext_from, measure.vars = geneset, variable.name = 'from_gene', value.name = 'from_expression')
  spatial_netw_ext_to_m[spatial_netw_ext_from_m, from_expression := from_expression, on = c(from='from',to_gene='from_gene')]


  # score expression
  # edge score == product
  comp_edge_score = spatial_netw_ext_to_m$to_expression * spatial_netw_ext_to_m$from_expression
  spatial_netw_ext_to_m[, edge_score := comp_edge_score]

  # gini scores
  gini_scores_DT = spatial_netw_ext_to_m[, .(spatial_gini_score = mygini_fun(rank(edge_score, ties.method = 'average'))), by = 'to_gene']

  # create copy for spatial random network
  spatial_random = copy(spatial_netw_ext_to_m)
  mygenes = unique(spatial_netw_ext_to_m$to_gene)

  # create original index
  spatial_netw_ext_to_m[, original_index := 1:.N]


  ## simulations
  if(simulations > 0) {

    for(sim in 1:simulations) {

      cat('start simulation ', sim, '\n')

      # OPTIM
      # create random index per gene and use to create random expression
      spatial_netw_ext_to_m[, random_index := sample(original_index, size = .N), by = 'to_gene']
      random_expr = spatial_netw_ext_to_m$from_expression[spatial_netw_ext_to_m$random_index]
      spatial_random[, from_expression := random_expr]
      cat('random expression is done \n')

      # average expression between random pairs and gini score
      new_comp_edge_score = spatial_random$to_expression * spatial_random$from_expression
      spatial_random[, edge_score := new_comp_edge_score]
      random_gini_score = spatial_random[, .(spatial_gini_score = mygini_fun(rank(edge_score, ties.method = 'average'))), by = 'to_gene']
      cat('random_gini score calculation \n')

      random_name = paste0('rand_gini_', sim)
      gini_scores_DT[, eval(random_name) := random_gini_score$spatial_gini_score]

    }

    random_columns = paste0('rand_gini_', 1:simulations)

    # calculate p-value, mean and sd of observed and random gini scores
    gini_scores_DT_m = data.table::melt.data.table(gini_scores_DT, measure.vars = random_columns, variable.name = 'random', value.name = 'random_scores')
    gini_scores_DT_m[, prob := ifelse(spatial_gini_score > random_scores, 0, 1)]

    summary_gini_scores = gini_scores_DT_m[, .(mean_score = mean(random_scores), sd_score = sd(random_scores), prob = sum(prob)/.N), by = .(to_gene, spatial_gini_score)]
    summary_gini_scores[, gini_ratio := spatial_gini_score/mean_score]
    summary_gini_scores <- merge(summary_gini_scores, detection_DT, by.x = 'to_gene', by.y = 'genes')

    # order
    setorder(summary_gini_scores, -gini_ratio, -detection)

  } else {

    summary_gini_scores <- merge(gini_scores_DT, detection_DT, by.x = 'to_gene', by.y = 'genes')

  }

  return(summary_gini_scores)

}


#' @title calculate_spatial_genes_python
#' @name calculate_spatial_genes_python
#' @description Calculate spatial genes using distance matrix.
#' @param gobject giotto object
#' @param expression_values expression values to use
#' @param metric distance metric to use
#' @param subset_genes only run on this subset of genes
#' @param rbp_p fractional binarization threshold
#' @param examine_top top fraction to evaluate with silhouette
#' @param python_path specify specific path to python if required
#' @return data.table with spatial scores
#' @details Description of how we compute spatial pattern genes.
#' @export
#' @examples
#'     calculate_spatial_genes_python(gobject)
calculate_spatial_genes_python <- function(gobject,
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
    python_path = system('which python', intern = T)
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







#' @title calculateSpatialGenes
#' @name calculateSpatialGenes
#' @description compute genes that are spatially clustered
#' @param gobject giotto object
#' @param expression_values expression values to use
#' @param method method to calculate spatial genes
#' @param spatial_network_name name of spatial network to use (default = 'spatial_network')
#' @param detection_threshold detection threshold to consider a gene detected
#' @param loess_span loess span for loess regression
#' @param pred_difference minimum difference between observed and predicted
#' @param split_gene_groups number of groups to split genes in
#' @param show_plot show plots
#' @param rank_percentage percentage of top cells for binarization
#' @param pvalue minimum p-value
#' @param OddsRatio minimum odds ratio
#' @param min_N minimum number of cells that need to display high expression upon binarization
#' @param max_N maximum number of cells that can display high expression upon binarization
#' @param SVname name for identified spatial genes (default = 'SV')
#' @param show_genes show top genes on plot
#' @param nr_genes # of genes to plot if show_genes = TRUE
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @return giotto object spatial genes appended to fDataDT
#' @details Description of how we compute spatial genes.
#' @export
#' @examples
#'     calculateSpatialGenes(gobject)
calculateSpatialGenes <- function(gobject,
                                  expression_values = c('normalized', 'scaled', 'custom'),
                                  method = c('kmeans',  'gini', 'rank'),
                                  spatial_network_name = 'spatial_network',
                                  simulations = 10,
                                  detection_threshold = 0,
                                  loess_span = 0.2,
                                  pred_difference = 0.01,
                                  split_gene_groups = 10,
                                  show_plot = T,
                                  rank_percentage = 10,
                                  pvalue = 0.01,
                                  OddsRatio = 2,
                                  min_N = 20,
                                  max_N = 5000,
                                  SVname = 'SV',
                                  show_genes = T,
                                  nr_genes = 20,
                                  return_gobject = T) {


  # spatial network is required to run this function
  if(is.null(gobject@spatial_network)) {
    stop('\n This function requires a spatial network, run createSpatialNetwork() first \n')
  } else if(is.null(gobject@spatial_network[[spatial_network_name]])) {
    stop('\n This function requires an existing spatial network name, run names(your_giotto_object@spatial_network) or make it with createSpatialNetwork() \n')
  }


  method = match.arg(method, choices = c('kmeans', 'gini', 'rank'))

  # split genes to run function in for loop OR in parallel
  # TODO: write parallel function

  # gini method
  # high score between pairs of cells #
  if(method == 'gini') {

    geneids = gobject@gene_ID
    mycut = cut(1:length(geneids), breaks = split_gene_groups, labels = paste0('group_', 1:split_gene_groups))
    names(geneids) <- mycut


    savelist <- list()
    for(group in unique(names(geneids))) {

      cat('\n \n START group ', group, '\n')

      selected_genes = geneids[names(geneids) == group]
      temp_giotto = subsetGiotto(gobject = gobject, gene_ids = selected_genes, cell_ids = NULL)
      temp_spatial = calculate_spatial_genes(gobject = temp_giotto, expression_values = expression_values,
                                             spatial_network_name = spatial_network_name,
                                             method = method, simulations = simulations, detection_threshold = detection_threshold)

      savelist[[group]] = temp_spatial

    }

    spatial_results = do.call('rbind', savelist)

    # TEMP
    #return(spatial_results)

    if(simulations > 0) {
      spatial_results[, gini_ratio := log2((spatial_gini_score+1)/(mean_score+1))]
    }

    # calculate prediction
    loess_model = stats::loess(formula = spatial_gini_score~detection, data = spatial_results, span = loess_span)

    # predict spatial gini ratio score based on detection
    predicted_gini = stats::predict(object = loess_model, newdata = spatial_results$detection, se = T)
    pred_gini = predicted_gini$fit
    pred_se_gini = predicted_gini$se.fit
    spatial_results[, c('pred_gini_score', 'pred_se_gini_score') := list(pred_gini, pred_se_gini)]
    spatial_results[, pred_diff := spatial_gini_score - pred_gini_score]
    spatial_results[, SVgenes := ifelse(spatial_gini_score > pred_gini_score + pred_difference, 'yes', 'no')]
    setorder(spatial_results, -pred_diff)

    if(show_plot == TRUE) {

      pl <- ggplot2::ggplot()
      pl <- pl + ggplot2::geom_point(data = spatial_results, aes(x = detection, y = spatial_gini_score, color = SVgenes))
      pl <- pl + ggplot2::geom_line(data = spatial_results, aes(x = detection, y = pred_gini_score, group = 1), color = 'blue', size = 1.5)
      pl <- pl + ggplot2::geom_line(data = spatial_results, aes(x = detection, y = pred_gini_score+pred_difference, group = 1), color = 'red', size = 1.5)
      if(show_genes == TRUE) {
        pl <- pl + ggrepel::geom_text_repel(data = spatial_results[SVgenes == 'yes'][1:nr_genes], aes(x = detection, y = spatial_gini_score, label = to_gene))
      }
      pl <- pl + ggplot2::labs(x = 'gene detection fraction', y = 'gini-score')
      print(pl)
    }




  }

  else if(method == 'kmeans') {
    spatial_results = calculate_binarized_spatial_genes(gobject = gobject, spatial_network_name = spatial_network_name,
                                                        bin_method = 'kmeans',
                                                        expression_values = expression_values)
    spatial_results[, SVgenes := ifelse(pval <= pvalue & OR > OddsRatio & N > min_N & N < max_N, 'yes', 'no')]
  }

  else if(method == 'rank') {
    spatial_results = calculate_binarized_spatial_genes(gobject = gobject, spatial_network_name = spatial_network_name,
                                                        bin_method = 'rank',
                                                        percentage_rank = rank_percentage,
                                                        expression_values = expression_values)
    spatial_results[, SVgenes := ifelse(pval <= pvalue & OR > OddsRatio & N > min_N & N < max_N, 'yes', 'no')]
  }





  if(return_gobject == TRUE) {

    gene_metadata = fDataDT(gobject)
    column_names_gene_metadata = colnames(gene_metadata)
    if(SVname %in% column_names_gene_metadata) {
      cat('\n ', SVname, ' has already been used, will be overwritten \n')
      gene_metadata[, eval(SVname) := NULL]
      gobject@gene_metadata = gene_metadata
    }

    if(method == 'gini') {
      SVgenes = spatial_results[,.(to_gene, SVgenes)]
      setnames(SVgenes, 'SVgenes', SVname)
      gobject <- addGeneMetadata(gobject = gobject, new_metadata = SVgenes, by_column = T, column_gene_ID = 'to_gene')
    } else if(method %in% c('kmeans', 'rank')) {
      SVgenes = spatial_results[,.(genes, SVgenes)]
      setnames(SVgenes, 'SVgenes', SVname)
      gobject <- addGeneMetadata(gobject = gobject, new_metadata = SVgenes, by_column = T, column_gene_ID = 'genes')
    }



    ## update parameters used ##
    parameters_list = gobject@parameters
    number_of_rounds = length(parameters_list)
    update_name = paste0(number_of_rounds,'_sv')

    # parameters to include
    parameters_list[[update_name]] = c('method used' = method,
                                       'expression values' = expression_values,
                                       'simulations' = simulations,
                                       'detection threshold' = detection_threshold,
                                       'loess span' = loess_span,
                                       'loess prediction difference' = pred_difference,
                                       'rank percentage' = rank_percentage,
                                       'pvalue' = pvalue,
                                       'Odds-Ratio' = OddsRatio,
                                       'min # of cells with high expression' = min_N,
                                       'max # of cells with high expression' = max_N,
                                       'Spatial Variable name' = SVname)
    gobject@parameters = parameters_list

    return(gobject)
  } else {
    return(spatial_results)
  }

}




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
#' @details Description of how we compute spatial pattern genes.
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
  spatial_locs = gobject@spatial_locs

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



#' @title showPattern
#' @name showPattern
#' @description create a spatial grid
#' @param spatPatObj Output from detectSpatialPatterns
#' @param dimension dimension to plot
#' @param trim Trim ends of the PC values.
#' @param background_color background color for plot
#' @param grid_border_color color for grid
#' @param show_legend show legend of ggplot
#' @param show_plot Show the plot.
#' @return ggplot
#' @details Description.
#' @export
#' @examples
#'     showPattern(gobject)
showPattern <- function(spatPatObj,
                        dimension = 1,
                        trim = c(0.02, 0.98),
                        background_color = 'white',
                        grid_border_color = 'grey',
                        show_legend = T,
                        show_plot = F) {

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

  # plot
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

  if(show_plot == TRUE) {
    print(dpl)
  }
  return(dpl)

}




#' @title showPatternGenes
#' @name showPatternGenes
#' @description create a spatial grid
#' @param spatPatObj Output from detectSpatialPatterns
#' @param dimension dimension to plot genes for.
#' @param top_pos_genes Top positively correlated genes.
#' @param top_neg_genes Top negatively correlated genes.
#' @param point_size size of points
#' @param show_plot Show the plot.
#' @return ggplot
#' @details Description.
#' @export
#' @examples
#'     showPatternGenes(gobject)
showPatternGenes <- function(spatPatObj,
                             dimension = 1,
                             top_pos_genes = 5,
                             top_neg_genes = 5,
                             point_size = 1,
                             show_plot = F) {

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
  gene_cor_DT = gene_cor_DT[order(get(selected_PC))]

  subset = gene_cor_DT[c(1:top_neg_genes, (nrow(gene_cor_DT)-top_pos_genes):nrow(gene_cor_DT))]
  subset[, gene_ID := factor(gene_ID, gene_ID)]

  pl <- ggplot()
  pl <- pl + ggplot2::theme_classic()
  pl <- pl + ggplot2::geom_point(data = subset, aes_string(x = selected_PC, y = 'gene_ID'), size = point_size)
  pl <- pl + ggplot2::geom_vline(xintercept = 0, linetype = 2)
  pl <- pl + ggplot2::labs(x = 'correlation', y = '', title = selected_PC)
  pl <- pl + ggplot2::theme(plot.title = element_text(hjust = 0.5))
  if(show_plot == TRUE) {
    print(pl)
  }
  return(pl)
}






#' @title selectPatternGenes
#' @name selectPatternGenes
#' @description create a spatial grid
#' @param spatPatObj Output from detectSpatialPatterns
#' @param dimensions dimensions to identify correlated genes for.
#' @param top_pos_genes Top positively correlated genes.
#' @param top_neg_genes Top negatively correlated genes.
#' @param min_pos_cor Minimum positive correlation score to include a gene.
#' @param min_neg_cor Minimum negative correlation score to include a gene.
#' @return ggplot
#' @details Description.
#' @export
#' @examples
#'     selectPatternGenes(gobject)
selectPatternGenes <- function(spatPatObj,
                               dimensions = 1:5,
                               top_pos_genes = 10,
                               top_neg_genes = 10,
                               min_pos_cor = 0.5,
                               min_neg_cor = -0.5) {


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











