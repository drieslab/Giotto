

#' @title cellProximityBarplot
#' @name cellProximityBarplot
#' @description Create barplot from cell-cell proximity scores
#' @param gobject giotto object
#' @param CPscore CPscore, output from cellProximityEnrichment()
#' @param min_orig_ints filter on minimum original cell-cell interactions
#' @param min_sim_ints filter on minimum simulated cell-cell interactions
#' @param p_val p-value
#' @param show_plot show plot
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters from \code{\link{all_plots_save_function}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return ggplot barplot
#' @details This function creates a barplot that shows the  spatial proximity
#'  enrichment or depletion of cell type pairs.
#' @export
cellProximityBarplot = function(gobject,
                                CPscore,
                                min_orig_ints = 5,
                                min_sim_ints = 5,
                                p_val = 0.05,
                                show_plot = NA,
                                return_plot = NA,
                                save_plot = NA,
                                save_param =  list(),
                                default_save_name = 'cellProximityBarplot') {


  table_mean_results_dc = CPscore$enrichm_res

  ## filter to remove low number of cell-cell proximity interactions ##

  # data.table variables
  original = simulations = p_higher_orig = p_lower_orig = enrichm = type_int = unified_int = NULL

  table_mean_results_dc_filter = table_mean_results_dc[original >= min_orig_ints & simulations >= min_sim_ints]
  table_mean_results_dc_filter = table_mean_results_dc_filter[p_higher_orig <= p_val | p_lower_orig <= p_val]

  pl <- ggplot2::ggplot()
  pl <- pl + ggplot2::geom_bar(data = table_mean_results_dc_filter, ggplot2::aes(x = unified_int, y = enrichm, fill = type_int), stat = 'identity', show.legend = F)
  pl <- pl + ggplot2::coord_flip()
  pl <- pl + ggplot2::theme_bw()
  pl <- pl + ggplot2::labs(y = 'enrichment/depletion')
  pl

  bpl <- ggplot2::ggplot()
  bpl <- bpl + ggplot2::geom_bar(data = table_mean_results_dc_filter, ggplot2::aes(x = unified_int, y = original, fill = type_int), stat = 'identity', show.legend = T)
  bpl <- bpl + ggplot2::coord_flip()
  bpl <- bpl + ggplot2::theme_bw() + ggplot2::theme(axis.text.y = element_blank())
  bpl <- bpl + ggplot2::labs(y = '# of interactions')
  bpl

  combo_plot <- cowplot::plot_grid(pl, bpl, ncol = 2, rel_heights = c(1), rel_widths = c(3,1.5), align = 'h')


  # print, return and save parameters
  show_plot = ifelse(is.na(show_plot), readGiottoInstructions(gobject, param = 'show_plot'), show_plot)
  save_plot = ifelse(is.na(save_plot), readGiottoInstructions(gobject, param = 'save_plot'), save_plot)
  return_plot = ifelse(is.na(return_plot), readGiottoInstructions(gobject, param = 'return_plot'), return_plot)

  ## print plot
  if(show_plot == TRUE) {
    print(combo_plot)
  }

  ## save plot
  if(save_plot == TRUE) {
    do.call('all_plots_save_function', c(list(gobject = gobject, plot_object = combo_plot, default_save_name = default_save_name), save_param))
  }

  ## return plot
  if(return_plot == TRUE) {
    return(combo_plot)
  }

}

#' @title cellProximityHeatmap
#' @name cellProximityHeatmap
#' @description Create heatmap from cell-cell proximity scores
#' @param gobject giotto object
#' @param CPscore CPscore, output from cellProximityEnrichment()
#' @param scale scale cell-cell proximity interaction scores
#' @param order_cell_types order cell types based on enrichment correlation
#' @param color_breaks numerical vector of length 3 to represent min, mean and maximum
#' @param color_names character color vector of length 3
#' @param show_plot show plot
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters from \code{\link{all_plots_save_function}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return ggplot heatmap
#' @details This function creates a heatmap that shows the  spatial proximity
#'  enrichment or depletion of cell type pairs.
#' @export
cellProximityHeatmap = function(gobject,
                                CPscore,
                                scale = T,
                                order_cell_types = T,
                                color_breaks = NULL,
                                color_names = NULL,
                                show_plot = NA,
                                return_plot = NA,
                                save_plot = NA,
                                save_param =  list(),
                                default_save_name = 'cellProximityHeatmap') {


  enrich_res = CPscore$enrichm_res

  # data.table variables
  first_type = second_type = unified_int = NULL

  enrich_res[, first_type := strsplit(x = as.character(unified_int), split = '--')[[1]][1], by = 1:nrow(enrich_res)]
  enrich_res[, second_type := strsplit(x = as.character(unified_int), split = '--')[[1]][2], by = 1:nrow(enrich_res)]

  # create matrix
  enrich_mat = data.table::dcast.data.table(data = enrich_res,formula = first_type~second_type, value.var = 'enrichm')
  matrix_d <- as.matrix(enrich_mat[,-1]); rownames(matrix_d) = as.vector(enrich_mat[[1]])
  t_matrix_d <- t_giotto(matrix_d)

  # fill in NAs based on values in upper and lower matrix triangle
  t_matrix_d[upper.tri(t_matrix_d)][is.na(t_matrix_d[upper.tri(t_matrix_d)])] = matrix_d[upper.tri(matrix_d)][is.na(t_matrix_d[upper.tri(t_matrix_d)])]
  t_matrix_d[lower.tri(t_matrix_d)][is.na(t_matrix_d[lower.tri(t_matrix_d)])] = matrix_d[lower.tri(matrix_d)][is.na(t_matrix_d[lower.tri(t_matrix_d)])]
  t_matrix_d[is.na(t_matrix_d)] = 0
  final_matrix = t_matrix_d

  # scale data
  if(scale == TRUE) {
    final_matrix <- t_giotto(scale(t_giotto(final_matrix)))
    final_matrix <- t_giotto(final_matrix)
    final_matrix[lower.tri(final_matrix)] <- t_giotto(final_matrix)[lower.tri(final_matrix)]
  }

  # # if NA values, impute as mean
  #if(any(is.na(final_matrix)) == TRUE) {
  #  myrowmeans = apply(X = final_matrix, MARGIN = 1, FUN = function(x) mean(na.omit(x)))
  #  mymatrixmeans = matrix(data = rep(myrowmeans, ncol(final_matrix)), nrow = nrow(final_matrix), ncol = ncol(final_matrix))
  #  final_matrix[is.na(final_matrix)] = mymatrixmeans[which(is.na(final_matrix))]
  #}

  # order cell types
  if(order_cell_types == TRUE) {

    cordist = stats::as.dist(1-cor_giotto(final_matrix))
    clus = stats::hclust(cordist)
    myorder = clus$order
    mylabels = clus$labels
    names(mylabels) = 1:length(mylabels)
    sample_order = mylabels[myorder]

    final_matrix = final_matrix[sample_order, sample_order]
  }

  # create custom colors or not
  if(!is.null(color_breaks) & !is.null(color_names)) {

    if(length(color_breaks) != 3 | !is.numeric(color_breaks)) {
      stop('\n color_breaks needs to be a numerical vector of length 3 \n')
    }

    if(length(color_names) != 3 | !is.character(color_names)) {
      stop('\n color_names needs to be a character vector of length 3 \n')
    }

    heatm = ComplexHeatmap::Heatmap(matrix = final_matrix, cluster_rows = F, cluster_columns = F,
                                    col = circlize::colorRamp2(breaks = color_breaks, colors = color_names))
  } else {
    heatm = ComplexHeatmap::Heatmap(matrix = final_matrix, cluster_rows = F, cluster_columns = F)
  }



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


#' @title cellProximityNetwork
#' @name cellProximityNetwork
#' @description Create network from cell-cell proximity scores
#' @param gobject giotto object
#' @param CPscore CPscore, output from cellProximityEnrichment()
#' @param remove_self_edges remove enrichment/depletion edges with itself
#' @param self_loop_strength size of self-loops
#' @param color_depletion color for depleted cell-cell interactions
#' @param color_enrichment color for enriched cell-cell interactions
#' @param rescale_edge_weights rescale edge weights (boolean)
#' @param edge_weight_range_depletion numerical vector of length 2 to rescale depleted edge weights
#' @param edge_weight_range_enrichment numerical vector of length 2 to rescale enriched edge weights
#' @param layout layout algorithm to use to draw nodes and edges
#' @param only_show_enrichment_edges show only the enriched pairwise scores
#' @param edge_width_range range of edge width
#' @param node_size size of nodes
#' @param node_text_size size of node labels
#' @param show_plot show plot
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters from \code{\link{all_plots_save_function}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return igraph plot
#' @details This function creates a network that shows the  spatial proximity
#'  enrichment or depletion of cell type pairs.
#' @export
cellProximityNetwork = function(gobject,
                                CPscore,
                                remove_self_edges = FALSE,
                                self_loop_strength = 0.1,
                                color_depletion = 'lightgreen',
                                color_enrichment = 'red',
                                rescale_edge_weights = TRUE,
                                edge_weight_range_depletion = c(0.1, 1),
                                edge_weight_range_enrichment = c(1, 5),
                                layout = c('Fruchterman', 'DrL', 'Kamada-Kawai'),
                                only_show_enrichment_edges = F,
                                edge_width_range = c(0.1, 2),
                                node_size = 4,
                                node_text_size = 6,
                                show_plot = NA,
                                return_plot = NA,
                                save_plot = NA,
                                save_param =  list(),
                                default_save_name = 'cellProximityNetwork') {

  # extract scores

  # data.table variables
  cell_1 = cell_2 = unified_int = color = size = name = NULL

  CPscores = CPscore[['enrichm_res']]
  CPscores[, cell_1 := strsplit(as.character(unified_int), split = '--')[[1]][1], by = 1:nrow(CPscores)]
  CPscores[, cell_2 := strsplit(as.character(unified_int), split = '--')[[1]][2], by = 1:nrow(CPscores)]

  # create igraph with enrichm as weight edges
  igd = igraph::graph_from_data_frame(d = CPscores[,c('cell_1', 'cell_2', 'enrichm')], directed = F)

  if(remove_self_edges == TRUE) {
    igd = igraph::simplify(graph = igd, remove.loops = TRUE, remove.multiple = FALSE)
  }

  edges_sizes = igraph::get.edge.attribute(igd, 'enrichm')
  post_edges_sizes = edges_sizes[edges_sizes > 0]
  neg_edges_sizes = edges_sizes[edges_sizes <= 0]

  # rescale if wanted
  if(rescale_edge_weights == TRUE) {
    pos_edges_sizes_resc = scales::rescale(x = post_edges_sizes, to = edge_weight_range_enrichment)
    neg_edges_sizes_resc = scales::rescale(x = neg_edges_sizes, to = edge_weight_range_depletion)
    edges_sizes_resc = c(pos_edges_sizes_resc, neg_edges_sizes_resc)
  } else {
    edges_sizes_resc = c(post_edges_sizes, neg_edges_sizes)
  }

  # colors
  edges_colors = ifelse(edges_sizes > 0, 'enriched', 'depleted')


  # create coordinates for layout
  if(class(layout) %in% c('data.frame', 'data.table')) {
    if(ncol(layout) < 2) {
      stop('custom layout needs to have at least 2 columns')
    }

    if(nrow(layout) != length(igraph::E(igd))) {
      stop('rows of custom layout need to be the same as number of edges')
    }

  } else {
    layout = match.arg(arg = layout, choices = c('Fruchterman', 'DrL', 'Kamada-Kawai'))
  }




  #iplot = igraph::plot.igraph(igd, edge.color = edges_colors, edge.width = edges_sizes_resc, layout = coords)

  igd = igraph::set.edge.attribute(graph = igd, index = igraph::E(igd), name = 'color', value = edges_colors)
  igd = igraph::set.edge.attribute(graph = igd, index = igraph::E(igd), name = 'size', value = as.numeric(edges_sizes_resc))

  ## only show attractive edges
  if(only_show_enrichment_edges == TRUE) {
    colors = igraph::get.edge.attribute(igd, name = 'color')
    subvertices_ids = which(colors == 'enriched')
    igd = igraph::subgraph.edges(graph = igd, eids = subvertices_ids)

    # get new rescale vector (in case vector id is lost)
    edges_sizes_resc = igraph::E(igd)$size
  }

  ## get coordinates layouts
  if(layout == 'Fruchterman') {
    coords = igraph::layout_with_fr(graph = igd, weights = edges_sizes_resc)
  } else if(layout == 'DrL') {
    coords = igraph::layout_with_drl(graph = igd, weights = edges_sizes_resc)
  } else if(layout == 'Kamada-Kawai') {
    coords = igraph::layout_with_kk(graph = igd, weights = edges_sizes_resc)
  } else {
    stop('\n Currently no other layouts have been implemented \n')
  }


  #longDT = as.data.table(igraph::as_long_data_frame(igd))
  #return(longDT)
  #return(list(igd, coords))

  ## create plot
  gpl = ggraph::ggraph(graph = igd, layout = coords)
  gpl = gpl + ggraph::geom_edge_link(ggplot2::aes(color = factor(color), edge_width = size, edge_alpha = size), show.legend = F)
  gpl = gpl + ggraph::geom_edge_loop(ggplot2::aes(color = factor(color), edge_width = size, edge_alpha = size, strength = self_loop_strength), show.legend = F)
  gpl = gpl + ggraph::scale_edge_color_manual(values = c('enriched' = color_enrichment, 'depleted' = color_depletion))
  gpl = gpl + ggraph::scale_edge_width(range = edge_width_range)
  gpl = gpl + ggraph::scale_edge_alpha(range = c(0.1,1))
  gpl = gpl + ggraph::geom_node_text(ggplot2::aes(label = name), repel = TRUE, size = node_text_size)
  gpl = gpl + ggraph::geom_node_point(size = node_size)
  gpl = gpl + ggplot2::theme_bw() + ggplot2::theme(panel.grid = ggplot2::element_blank(),
                                                   panel.border = ggplot2::element_blank(),
                                                   axis.title = ggplot2::element_blank(),
                                                   axis.text = ggplot2::element_blank(),
                                                   axis.ticks = ggplot2::element_blank())


  # print, return and save parameters
  show_plot = ifelse(is.na(show_plot), readGiottoInstructions(gobject, param = 'show_plot'), show_plot)
  save_plot = ifelse(is.na(save_plot), readGiottoInstructions(gobject, param = 'save_plot'), save_plot)
  return_plot = ifelse(is.na(return_plot), readGiottoInstructions(gobject, param = 'return_plot'), return_plot)

  ## print plot
  if(show_plot == TRUE) {
    print(gpl)
  }

  ## save plot
  if(save_plot == TRUE) {
    do.call('all_plots_save_function', c(list(gobject = gobject, plot_object = gpl, default_save_name = default_save_name), save_param))
  }

  ## return plot
  if(return_plot == TRUE) {
    return(gpl)
  }


}


#' @title cellProximityVisPlot_2D_ggplot
#' @name cellProximityVisPlot_2D_ggplot
#' @description Visualize 2D cell-cell interactions according to spatial coordinates in ggplot mode
#' @return ggplot
#' @details Description of parameters.
#' @keywords internal
cellProximityVisPlot_2D_ggplot <- function(gobject,
                                           interaction_name = NULL,
                                           cluster_column = NULL,
                                           sdimx = NULL,
                                           sdimy = NULL,
                                           cell_color = NULL,
                                           cell_color_code = NULL,
                                           color_as_factor = T,
                                           show_other_cells = F,
                                           show_network = F,
                                           show_other_network = F,
                                           network_color = NULL,
                                           spatial_network_name = 'Delaunay_network',
                                           show_grid = F,
                                           grid_color = NULL,
                                           spatial_grid_name = 'spatial_grid',
                                           coord_fix_ratio = 1,
                                           show_legend = T,
                                           point_size_select = 2,
                                           point_select_border_col = 'black',
                                           point_select_border_stroke = 0.05,
                                           point_size_other = 1,
                                           point_alpha_other = 0.3,
                                           point_other_border_col = 'lightgrey',
                                           point_other_border_stroke = 0.01,
                                           ...){

  # data.table variables
  unified_int = sdimx_begin = sdimy_begin = sdimx_end = sdimy_end = x_start = x_end = NULL
  y_start = y_end = cell_ID = NULL

  if(is.null(interaction_name)) {
    stop('\n you need to specific at least one interaction name, run cellProximityEnrichment \n')
  }


  cell_locations  = gobject@spatial_locs
  spatial_grid    = gobject@spatial_grid[[spatial_grid_name]]
  cell_metadata   = gobject@cell_metadata



  spatial_network = annotateSpatialNetwork(gobject = gobject,
                                           spatial_network_name = spatial_network_name,
                                           cluster_column = cluster_column)

  cell_IDs_to_keep = unique(c(spatial_network[unified_int %in% interaction_name]$to,
                              spatial_network[unified_int %in% interaction_name]$from))

  print(cell_IDs_to_keep)

  if(show_other_cells){
    CellType <- strsplit(interaction_name,"--")
    all_cell_IDs = cell_metadata[cell_metadata[[cluster_column]] == CellType[[1]][1] |
                                   cell_metadata[[cluster_column]] == CellType[[1]][2],]$cell_ID
    other_cell_IDs <- setdiff(all_cell_IDs, cell_IDs_to_keep)
  }


  # annotated cell data
  if(nrow(cell_metadata) == 0) {
    cell_locations_metadata = cell_locations
  } else {
    cell_locations_metadata <- merge(cell_locations, cell_metadata,by = "cell_ID")
  }


  # first 2 dimensions need to be defined
  if(is.null(sdimx) | is.null(sdimy)) {
    cat('first and second dimenion need to be defined, default is first 2 \n')
    sdimx = 'sdimx'
    sdimy = 'sdimy'
  }

  pl <- ggplot2::ggplot()
  pl <- pl + ggplot2::theme_classic()

  if(!is.null(spatial_network) & show_network == TRUE) {
    if(is.null(network_color)) network_color = 'red'
    if(show_other_network){
      pl <- pl + ggplot2::geom_segment(data = spatial_network[!unified_int %in% interaction_name],
                                       aes(x = sdimx_begin, y = sdimy_begin, xend = sdimx_end, yend = sdimy_end),
                                       color = 'lightgrey', size = 0.5, alpha = 0.5)
    }
    pl <- pl + ggplot2::geom_segment(data = spatial_network[unified_int %in% interaction_name],
                                     aes(x = sdimx_begin, y = sdimy_begin, xend = sdimx_end, yend = sdimy_end),
                                     color = network_color, size = 0.5, alpha = 0.5)
  }

  if(!is.null(spatial_grid) & show_grid == TRUE) {
    if(is.null(grid_color)) grid_color = 'black'
    pl <- pl + ggplot2::geom_rect(data = spatial_grid, aes(xmin = x_start, xmax = x_end, ymin = y_start, ymax = y_end),
                                  color = grid_color, fill = NA)
  }

  # cell color default
  if(is.null(cell_color)) {
    cell_color = 'lightblue'
    pl <- pl + ggplot2::geom_point(data = cell_locations[!cell_ID %in% cell_IDs_to_keep], aes_string(x = sdimx, y = sdimy),
                                   show.legend = show_legend, shape = 21, fill = 'lightgrey', size = point_size_other)
    pl <- pl + ggplot2::geom_point(data = cell_locations[cell_ID %in% cell_IDs_to_keep], aes_string(x = sdimx, y = sdimy),
                                   show.legend = show_legend, shape = 21, fill = cell_color, size = point_size_select)
    if(show_other_cells){
      pl <- pl + ggplot2::geom_point(data = cell_locations[cell_ID %in% other_cell_IDs], aes_string(x = sdimx, y = sdimy),
                                     show.legend = show_legend, shape = 21, fill = cell_color, alpha = point_alpha_other,
                                     size = point_size_select * 0.5)
    }
  }
  else if (is.character(cell_color)) {
    if(cell_color %in% colnames(cell_locations_metadata)) {

      if(color_as_factor == TRUE) {
        factor_data = factor(cell_locations_metadata[[cell_color]])
        cell_locations_metadata[[cell_color]] <- factor_data
      }

      pl <- pl + ggplot2::geom_point(data = cell_locations_metadata[!cell_ID %in% cell_IDs_to_keep], aes_string(x = sdimx, y = sdimy),
                                     fill = 'lightgrey', shape = 21, size = point_size_other,
                                     color = point_other_border_col, stroke = point_other_border_stroke)
      pl <- pl + ggplot2::geom_point(data = cell_locations_metadata[cell_ID %in% cell_IDs_to_keep], aes_string(x = sdimx, y = sdimy, fill = cell_color),
                                     show.legend = show_legend, shape = 21, size = point_size_select,
                                     color = point_select_border_col, stroke = point_select_border_stroke)
      if(show_other_cells){
        pl <- pl + ggplot2::geom_point(data = cell_locations_metadata[cell_ID %in% other_cell_IDs], aes_string(x = sdimx, y = sdimy,fill = cell_color),
                                       show.legend = show_legend, shape = 21, alpha = point_alpha_other,
                                       size = point_size_select * 0.5)
      }



      if(!is.null(cell_color_code)) {
        pl <- pl + ggplot2::scale_fill_manual(values = cell_color_code)
      } else if(color_as_factor == T) {
        number_colors = length(unique(factor_data))
        cell_color_code = getDistinctColors(n = number_colors)
        names(cell_color_code) = unique(factor_data)
        pl <- pl + ggplot2::scale_fill_manual(values = cell_color_code)
      } else if(color_as_factor == F){
        pl <- pl + ggplot2::scale_fill_gradient(low = 'blue', high = 'red')
      }

    } else {
      pl <- pl + ggplot2::geom_point(data = cell_locations_metadata[!cell_ID %in% cell_IDs_to_keep], aes_string(x = sdimx, y = sdimy),
                                     show.legend = show_legend, shape = 21, fill = 'lightgrey', size = point_size_other,
                                     color = point_other_border_col, stroke = point_other_border_stroke)
      pl <- pl + ggplot2::geom_point(data = cell_locations_metadata[cell_ID %in% cell_IDs_to_keep], aes_string(x = sdimx, y = sdimy),
                                     show.legend = show_legend, shape = 21, fill = cell_color, size = point_size_select,
                                     color = point_select_border_col, stroke = point_select_border_stroke)
    }

  }

  pl <- pl + ggplot2::theme_bw() + ggplot2::theme(plot.title = element_text(hjust = 0.5),
                                                  legend.title = element_text(size = 10),
                                                  legend.text = element_text(size = 10))

  # fix coord ratio
  if(!is.null(coord_fix_ratio)) {
    pl <- pl + ggplot2::coord_fixed(ratio = coord_fix_ratio)
  }

  pl <- pl + ggplot2::labs(x = 'x coordinates', y = 'y coordinates')

  return(pl)
}


#' @title cellProximityVisPlot_2D_plotly
#' @name cellProximityVisPlot_2D_plotly
#' @description Visualize 2D cell-cell interactions according to spatial coordinates in plotly mode
#' @keywords internal
cellProximityVisPlot_2D_plotly <- function(gobject,
                                           interaction_name = NULL,
                                           cluster_column = NULL,
                                           sdimx = NULL,
                                           sdimy = NULL,
                                           cell_color = NULL,
                                           cell_color_code = NULL,
                                           color_as_factor = T,
                                           show_other_cells = F,
                                           show_network = F,
                                           show_other_network = F,
                                           network_color = NULL,
                                           spatial_network_name = 'Delaunay_network',
                                           show_grid = F,
                                           grid_color = NULL,
                                           spatial_grid_name = 'spatial_grid',
                                           show_legend = T,
                                           point_size_select = 2,
                                           point_size_other = 1,
                                           point_alpha_other = 0.3,
                                           axis_scale = c("cube","real","custom"),
                                           custom_ratio = NULL,
                                           x_ticks = NULL,
                                           y_ticks = NULL,
                                           ...){


  # data.table variables
  cell_ID = unified_int = NULL

  if(is.null(interaction_name)) {
    stop('\n you need to specific at least one interaction name, run cellProximityEnrichment \n')
  }


  cell_locations  = gobject@spatial_locs
  spatial_grid    = gobject@spatial_grid[[spatial_grid_name]]
  cell_metadata   = gobject@cell_metadata


  spatial_network = annotateSpatialNetwork(gobject = gobject, spatial_network_name = spatial_network_name, cluster_column = cluster_column)

  cell_IDs_to_keep = unique(c(spatial_network[unified_int %in% interaction_name]$to, spatial_network[unified_int %in% interaction_name]$from))

  if(show_other_cells){
    CellType <- strsplit(interaction_name,"-")
    all_cell_IDs = cell_metadata[cell_metadata[[cluster_column]] == CellType[[1]][1] |
                                   cell_metadata[[cluster_column]] == CellType[[1]][2],]$cell_ID
    other_cell_IDs <- setdiff(all_cell_IDs, cell_IDs_to_keep)
  }

  # annotated cell data
  if(nrow(cell_metadata) == 0) {
    cell_locations_metadata = cell_locations
  } else {
    cell_locations_metadata <- merge(cell_locations, cell_metadata, by = "cell_ID")
  }




  # first 2 dimensions need to be defined
  if(is.null(sdimx) | is.null(sdimy)) {
    cat('first and second dimenion need to be defined, default is first 2 \n')
    sdimx = 'sdimx'
    sdimy = 'sdimy'
  }


  cat('create 2D plotly plot')

  axis_scale = match.arg(axis_scale, c("cube","real","custom"))

  ratio = plotly_axis_scale_2D(cell_locations_metadata,sdimx = sdimx,sdimy = sdimy,
                               mode = axis_scale,custom_ratio = custom_ratio)

  pl <- plotly::plot_ly()

  if(show_network == TRUE) {
    if(is.null(network_color)){
      network_color = "red"
    }
    if(show_other_network){
      pl <- pl %>% plotly::add_segments(name = paste("unselected",spatial_network_name,sep = " "),
                                        type = "scatter",
                                        data = spatial_network[!unified_int %in% interaction_name],
                                        x = ~sdimx_begin,
                                        y =~sdimy_begin,
                                        xend = ~sdimx_end,
                                        yend = ~sdimy_end,
                                        line = list(color = "lightgrey",
                                                    width = 0.5),
                                        opacity=0.3)
    }
    pl <- pl %>% plotly::add_segments(name = spatial_network_name,
                                      type = "scatter",
                                      data = spatial_network[unified_int %in% interaction_name],
                                      x = ~sdimx_begin,
                                      y = ~sdimy_begin,
                                      xend = ~sdimx_end,
                                      yend = ~sdimy_end,
                                      line = list(color = network_color,
                                                  width = 0.5),
                                      opacity=0.8)

  }
  if(show_grid == TRUE){
    if(is.null(grid_color)) {
      grid_color = 'black'
    }
    edges <- plotly_grid(spatial_grid)
    spl <- spl %>% plotly::add_segments(name = "spatial_grid",
                                        type = "scatter",
                                        data = edges,
                                        x = ~x,
                                        y = ~y,
                                        xend = ~x_end,
                                        yend = ~y_end,
                                        line = list(color = grid_color,
                                                    width = 1),
                                        opacity=1)

  }


  if(!is.null(cell_color)) {
    if(cell_color %in% colnames(cell_locations_metadata)){
      if(is.null(cell_color_code)) {
        number_colors=length(unique(cell_locations_metadata[[cell_color]]))
        cell_color_code = getDistinctColors(n = number_colors)
      }
      cell_locations_metadata[[cell_color]] <- as.factor(cell_locations_metadata[[cell_color]])

      pl <- pl %>% plotly::add_trace(type = 'scatter',mode = 'markers',
                                     #name = "selected cells",
                                     data=cell_locations_metadata[cell_ID %in% cell_IDs_to_keep],
                                     x = ~sdimx, y = ~sdimy,
                                     color = cell_locations_metadata[cell_ID %in% cell_IDs_to_keep][[cell_color]],
                                     colors = cell_color_code,
                                     marker = list(size = point_size_select))
      if(show_other_cells){
        pl <- pl %>% plotly::add_trace(type = 'scatter',mode = 'markers',
                                       #name = "selected cells outside network",
                                       data=cell_locations_metadata[cell_ID %in% other_cell_IDs],
                                       x = ~sdimx, y = ~sdimy,
                                       color = cell_locations_metadata[cell_ID %in% other_cell_IDs][[cell_color]],
                                       colors = cell_color_code,
                                       opacity = point_alpha_other,
                                       marker = list(size = point_size_select * 0.7))
      }
      pl <- pl %>%  plotly::add_trace(type = 'scatter',mode = "markers",
                                      name = "unselected cells",
                                      data=cell_locations_metadata[!cell_ID %in% cell_IDs_to_keep],
                                      x = ~sdimx, y = ~sdimy,
                                      marker = list(size = point_size_other,color = "lightgray",colors = "lightgray"),
                                      opacity = point_alpha_other)

    }
    else{
      cat('cell_color not exist!\n')
    }
  }
  else {
    pl <- pl %>% plotly::add_trace(type = 'scatter',mode = 'markers',
                                   name = "selected cells",
                                   data=cell_locations_metadata[cell_ID %in% cell_IDs_to_keep],
                                   x = ~sdimx, y = ~sdimy,
                                   marker = list(size = point_size_select,color = "lightblue",colors = "lightblue"))
    if(show_other_cells){
      pl <- pl %>% plotly::add_trace(type = 'scatter',mode = 'markers',
                             data=cell_locations_metadata[cell_ID %in% other_cell_IDs],
                             x = ~sdimx, y = ~sdimy,
                             name = "selected cells outside network",
                             marker = list(size = point_size_select*0.7,color = "lightblue",colors = "lightblue"),
                             opacity = point_alpha_other)
    }
    pl <-  pl %>%  plotly::add_trace(type = 'scatter',mode = "markers",
                                     name = "unselected cells",
                                     data=cell_locations_metadata[!cell_ID %in% cell_IDs_to_keep],
                                     x = ~sdimx, y = ~sdimy,
                                     marker = list(size = point_size_other,color = "lightgray",colors = "lightgray"),
                                     opacity = point_alpha_other)

  }

  pl <- pl %>%
    plotly::layout(list(xaxis = list(title = 'X',nticks = x_ticks),
                        yaxis = list(title = 'Y',nticks = y_ticks)),
                   legend = list(x = 100, y = 0.5,
                                 font = list(family = "sans-serif",size = 12)))
  return((pl))

}

#' @title cellProximityVisPlot_3D_plotly
#' @name cellProximityVisPlot_3D_plotly
#' @description Visualize 3D cell-cell interactions according to spatial coordinates in plotly mode
#' @keywords internal
cellProximityVisPlot_3D_plotly <- function(gobject,
                                           interaction_name = NULL,
                                           cluster_column = NULL,
                                           sdimx = NULL,
                                           sdimy = NULL,
                                           sdimz = NULL,
                                           cell_color = NULL,
                                           cell_color_code = NULL,
                                           color_as_factor = T,
                                           show_other_cells = F,
                                           show_network = F,
                                           show_other_network = F,
                                           network_color = NULL,
                                           spatial_network_name = 'Delaunay_network',
                                           show_grid = F,
                                           grid_color = NULL,
                                           spatial_grid_name = 'spatial_grid',
                                           show_legend = T,
                                           point_size_select = 2,
                                           point_size_other = 1,
                                           point_alpha_other = 0.5,
                                           axis_scale = c("cube","real","custom"),
                                           custom_ratio = NULL,
                                           x_ticks = NULL,
                                           y_ticks = NULL,
                                           z_ticks = NULL,
                                           ...){

  # data.table variables
  cell_ID = unified_int = NULL

  if(is.null(interaction_name)) {
    stop('\n you need to specific at least one interaction name, run cellProximityEnrichment \n')
  }


  cell_locations  = gobject@spatial_locs
  spatial_grid    = gobject@spatial_grid[[spatial_grid_name]]
  cell_metadata   = gobject@cell_metadata


  spatial_network = annotateSpatialNetwork(gobject = gobject, spatial_network_name = spatial_network_name, cluster_column = cluster_column)

  cell_IDs_to_keep = unique(c(spatial_network[unified_int %in% interaction_name]$to, spatial_network[unified_int %in% interaction_name]$from))

  if(show_other_cells){
    CellType <- strsplit(interaction_name,"-")
    all_cell_IDs = cell_metadata[cell_metadata[[cluster_column]] == CellType[[1]][1] |
                                   cell_metadata[[cluster_column]] == CellType[[1]][2],]$cell_ID
    other_cell_IDs <- setdiff(all_cell_IDs, cell_IDs_to_keep)
  }

  # annotated cell data
  if(nrow(cell_metadata) == 0) {
    cell_locations_metadata = cell_locations
  } else {
    cell_locations_metadata <- merge(cell_locations, cell_metadata, by = "cell_ID")
  }




  # first 2 dimensions need to be defined
  if(is.null(sdimx) | is.null(sdimy)) {
    cat('first and second dimenion need to be defined, default is first 2 \n')
    sdimx = 'sdimx'
    sdimy = 'sdimy'
  }





  # if 3 dimensions are defined create a 3D plot

  cat('create 3D plotly plot')

  pl <- plotly::plot_ly()

  axis_scale = match.arg(axis_scale, c("cube","real","custom"))

  ratio = plotly_axis_scale_3D(cell_locations_metadata,sdimx = sdimx,sdimy = sdimy,sdimz = sdimz,
                               mode = axis_scale,custom_ratio = custom_ratio)

  if(!is.null(cell_color)) {
    if(cell_color %in% colnames(cell_locations_metadata)){
      if(is.null(cell_color_code)) {
        number_colors=length(unique(cell_locations_metadata[[cell_color]]))
        cell_color_code = getDistinctColors(n = number_colors)
      }
      cell_locations_metadata[[cell_color]] <- as.factor(cell_locations_metadata[[cell_color]])

      pl <- pl %>% plotly::add_trace(type = 'scatter3d',mode = 'markers',
                                     #name = "selected cells",
                                     data=cell_locations_metadata[cell_ID %in% cell_IDs_to_keep],
                                     x = ~sdimx, y = ~sdimy, z = ~sdimz,
                                     color = cell_locations_metadata[cell_ID %in% cell_IDs_to_keep][[cell_color]],
                                     colors = cell_color_code,
                                     marker = list(size = point_size_select))%>%
        plotly::add_trace(type = 'scatter3d',mode = "markers",
                          name = "unselected cells",
                          data=cell_locations_metadata[!cell_ID %in% cell_IDs_to_keep],
                          x = ~sdimx, y = ~sdimy, z = ~sdimz,
                          marker = list(size = point_size_other,color = "lightgray",colors = "lightgray"),
                          opacity = point_alpha_other)
      if(show_other_cells){
        pl <- pl %>% plotly::add_trace(type = 'scatter3d',mode = 'markers',
                                       #name = "selected cells outside network",
                                       data=cell_locations_metadata[cell_ID %in% other_cell_IDs],
                                       x = ~sdimx, y = ~sdimy, z = ~sdimz,
                                       color = cell_locations_metadata[cell_ID %in% other_cell_IDs][[cell_color]],
                                       colors = cell_color_code,
                                       opacity = point_alpha_other,
                                       marker = list(size = point_size_select*0.7))
      }
    }
    else{
      cat('cell_color not exist!\n')
    }
  }
  else {
    pl <- pl %>% plotly::add_trace(type = 'scatter3d',mode = 'markers',
                                   name = "selected cells",
                                   data=cell_locations_metadata[cell_ID %in% cell_IDs_to_keep],
                                   x = ~sdimx, y = ~sdimy, z = ~sdimz,
                                   marker = list(size = point_size_select,color = "lightblue",colors = "lightblue"))%>%
      plotly::add_trace(type = 'scatter3d',mode = "markers",
                        name = "unselected cells",
                        data=cell_locations_metadata[!cell_ID %in% cell_IDs_to_keep],
                        x = ~sdimx, y = ~sdimy, z = ~sdimz,
                        marker = list(size = point_size_other,color = "lightgray",colors = "lightgray"),
                        opacity = point_alpha_other)
    if(show_other_cells){
      pl <- pl %>% plotly::add_trace(type = 'scatter3d',mode = 'markers',
                             data=cell_locations_metadata[cell_ID %in% other_cell_IDs],
                             x = ~sdimx, y = ~sdimy, z = ~sdimz,
                             name = "selected cells outside network",
                             marker = list(size = point_size_select*0.7,color = "lightblue",colors = "lightblue"),
                             opacity = point_alpha_other)
    }
  }
  if(!is.null(spatial_network) & show_network == TRUE) {
    if(is.null(network_color)) {
      network_color = 'red'
    }
    unselect_network <- spatial_network[!unified_int %in% interaction_name]
    select_network <- spatial_network[unified_int %in% interaction_name]
    pl <- pl %>% plotly::add_trace(name = "sptial network",mode = "lines", type = "scatter3d",opacity=0.5,
                                   data = plotly_network(select_network),
                                   x = ~x,y=~y,z=~z,inherit = F,line=list(color=network_color))
    if(show_other_network == T){
      pl <- pl %>% plotly::add_trace(name = "unselected sptial network",mode = "lines", type = "scatter3d",opacity=0.1,
                                     data = plotly_network(unselect_network),
                                     x = ~x,y=~y,z=~z,inherit = F,line=list(color="lightgray"))
    }


  }

  pl <- pl %>% plotly::layout(scene = list(
    xaxis = list(title = "X",nticks = x_ticks),
    yaxis = list(title = "Y",nticks = y_ticks),
    zaxis = list(title = "Z",nticks = z_ticks),
    aspectmode='manual',
    aspectratio = list(x=ratio[[1]],
                       y=ratio[[2]],
                       z=ratio[[3]])))
  return(pl)


}

#' @title cellProximityVisPlot
#' @name cellProximityVisPlot
#' @description Visualize cell-cell interactions according to spatial coordinates
#' @param gobject giotto object
#' @param interaction_name cell-cell interaction name
#' @param cluster_column cluster column with cell clusters
#' @param sdimx x-axis dimension name (default = 'sdimx')
#' @param sdimy y-axis dimension name (default = 'sdimy')
#' @param sdimz z-axis dimension name (default = 'sdimz')
#' @param cell_color color for cells (see details)
#' @param cell_color_code named vector with colors
#' @param color_as_factor convert color column to factor
#' @param show_other_cells show not selected cells
#' @param show_network show underlying spatial network
#' @param show_other_network show underlying spatial network of other cells
#' @param network_color color of spatial network
#' @param spatial_network_name name of spatial network to use
#' @param show_grid show spatial grid
#' @param grid_color color of spatial grid
#' @param spatial_grid_name name of spatial grid to use
#' @param coord_fix_ratio fix ratio between x and y-axis
#' @param show_legend show legend
#' @param point_size_select size of selected points
#' @param point_select_border_col border color of selected points
#' @param point_select_border_stroke stroke size of selected points
#' @param point_size_other size of other points
#' @param point_alpha_other alpha of other points
#' @param point_other_border_col border color of other points
#' @param point_other_border_stroke stroke size of other points
#' @param axis_scale scale of axis
#' @param custom_ratio custom ratio of scales
#' @param x_ticks x ticks
#' @param y_ticks y ticks
#' @param z_ticks z ticks
#' @param plot_method method to plot
#' @param \dots additional parameters
#' @return ggplot or plotly
#' @details Description of parameters.
#' @export
cellProximityVisPlot <- function(gobject,
                                 interaction_name = NULL,
                                 cluster_column = NULL,
                                 sdimx = NULL,
                                 sdimy = NULL,
                                 sdimz = NULL,
                                 cell_color = NULL,
                                 cell_color_code = NULL,
                                 color_as_factor = T,
                                 show_other_cells = F,
                                 show_network = F,
                                 show_other_network = F,
                                 network_color = NULL,
                                 spatial_network_name = 'Delaunay_network',
                                 show_grid = F,
                                 grid_color = NULL,
                                 spatial_grid_name = 'spatial_grid',
                                 coord_fix_ratio = 1,
                                 show_legend = T,
                                 point_size_select = 2,
                                 point_select_border_col = 'black',
                                 point_select_border_stroke = 0.05,
                                 point_size_other = 1,
                                 point_alpha_other = 0.3,
                                 point_other_border_col = 'lightgrey',
                                 point_other_border_stroke = 0.01,
                                 axis_scale = c("cube","real","custom"),
                                 custom_ratio = NULL,
                                 x_ticks = NULL,
                                 y_ticks = NULL,
                                 z_ticks = NULL,
                                 plot_method = c('ggplot', 'plotly'),
                                 ...) {


  ## decide plot method
  plot_method = match.arg(plot_method, choices = c('ggplot', 'plotly'))
  axis_scale = match.arg(axis_scale, c("cube","real","custom"))


  if(plot_method == 'ggplot') {

    if(is.null(sdimx) | is.null(sdimy)) {

      warning("plot_method = ggplot, but spatial dimensions for sdimx and sdimy for 2D plotting are not given. \n
              It will default to the 'sdimx' and 'sdimy' ")
      sdimx = 'sdimx'
      sdimy = 'sdimy'
      #stop('\n ggplot is in 2D and you need to define sdimx and sdimy \n')
    }

    if(length(c(sdimx, sdimy, sdimz)) == 3){
      warning("ggplot is not able to produce 3D plot! Please choose plotly method\n")
    }
    result = cellProximityVisPlot_2D_ggplot(gobject = gobject,
                                            interaction_name = interaction_name,
                                            cluster_column = cluster_column,
                                            sdimx = sdimx,
                                            sdimy = sdimy,
                                            cell_color = cell_color,
                                            cell_color_code = cell_color_code,
                                            color_as_factor = color_as_factor,
                                            show_other_cells = show_other_cells,
                                            show_network = show_network,
                                            show_other_network = show_other_network,
                                            network_color = network_color,
                                            spatial_network_name = spatial_network_name,
                                            show_grid = show_grid,
                                            grid_color = grid_color,
                                            spatial_grid_name = spatial_grid_name,
                                            coord_fix_ratio = coord_fix_ratio,
                                            show_legend = show_legend,
                                            point_size_select = point_size_select,
                                            point_select_border_col = point_select_border_col,
                                            point_select_border_stroke = point_select_border_stroke,
                                            point_size_other = point_size_other,
                                            point_alpha_other =point_alpha_other,
                                            point_other_border_col = point_other_border_col,
                                            point_other_border_stroke = point_other_border_stroke,
                                            ...)

  }
  else if(plot_method == 'plotly') {

    if(length(c(sdimx, sdimy, sdimz)) == 3) {

      result = cellProximityVisPlot_3D_plotly(gobject = gobject,
                                              interaction_name = interaction_name,
                                              cluster_column = cluster_column,
                                              sdimx = sdimx,
                                              sdimy = sdimy,
                                              sdimz = sdimz,
                                              cell_color = cell_color,
                                              cell_color_code = cell_color_code,
                                              color_as_factor = color_as_factor,
                                              show_other_cells = show_other_cells,
                                              show_network = show_network,
                                              show_other_network = show_other_network,
                                              network_color = network_color,
                                              spatial_network_name = spatial_network_name,
                                              show_grid = show_grid,
                                              grid_color = grid_color,
                                              spatial_grid_name = spatial_grid_name,
                                              show_legend = show_legend,
                                              point_size_select = point_size_select,
                                              point_size_other = point_size_other,
                                              point_alpha_other = point_alpha_other,
                                              axis_scale = axis_scale,
                                              custom_ratio = custom_ratio,
                                              x_ticks = x_ticks,
                                              y_ticks = y_ticks,
                                              z_ticks = z_ticks,
                                              ...)

    }
    else {

      if(is.null(sdimx) | is.null(sdimy)) {
        stop('\n plotly in 2D requires you to define sdimx and sdimy \n')
      }

      ## run: visPlot_2D_plotly
      result = cellProximityVisPlot_2D_plotly(gobject = gobject,
                                              interaction_name = interaction_name,
                                              cluster_column = cluster_column,
                                              sdimx = sdimx,
                                              sdimy = sdimy,
                                              cell_color = cell_color,
                                              cell_color_code = cell_color_code,
                                              color_as_factor = color_as_factor,
                                              show_other_cells = show_other_cells,
                                              show_network = show_network,
                                              show_other_network = show_other_network,
                                              network_color = network_color,
                                              spatial_network_name = spatial_network_name,
                                              show_grid = show_grid,
                                              grid_color = grid_color,
                                              spatial_grid_name = spatial_grid_name,
                                              show_legend = show_legend,
                                              point_size_select = point_size_select,
                                              point_size_other = point_size_other,
                                              point_alpha_other = point_alpha_other,
                                              axis_scale = axis_scale,
                                              custom_ratio = custom_ratio,
                                              x_ticks = x_ticks,
                                              y_ticks = y_ticks,
                                              ...)


    }

  }
  return(result)

}







# * ####
# NEW plots ####

#' @title plotCellProximityGenes
#' @name plotCellProximityGenes
#' @description Create visualization for cell proximity gene scores
#' @param gobject giotto object
#' @param cpgObject ICG (interaction changed gene) score object
#' @param method plotting method to use
#' @param min_cells minimum number of source cell type
#' @param min_cells_expr minimum expression level for source cell type
#' @param min_int_cells minimum number of interacting neighbor cell type
#' @param min_int_cells_expr minimum expression level for interacting neighbor cell type
#' @param min_fdr minimum adjusted p-value
#' @param min_spat_diff minimum absolute spatial expression difference
#' @param min_log2_fc minimum log2 fold-change
#' @param min_zscore minimum z-score change
#' @param zscores_column calculate z-scores over cell types or genes
#' @param direction differential expression directions to keep
#' @param cell_color_code vector of colors with cell types as names
#' @param show_plot show plots
#' @param return_plot return plotting object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters from \code{\link{all_plots_save_function}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return plot
#' @export
plotCellProximityGenes = function(gobject,
                                  cpgObject,
                                  method = c('volcano', 'cell_barplot', 'cell-cell', 'cell_sankey', 'heatmap', 'dotplot'),
                                  min_cells = 4,
                                  min_cells_expr = 1,
                                  min_int_cells = 4,
                                  min_int_cells_expr = 1,
                                  min_fdr = 0.1,
                                  min_spat_diff = 0.2,
                                  min_log2_fc = 0.2,
                                  min_zscore = 2,
                                  zscores_column = c('cell_type', 'genes'),
                                  direction = c('both', 'up', 'down'),
                                  cell_color_code = NULL,
                                  show_plot = NA,
                                  return_plot = NA,
                                  save_plot = NA,
                                  save_param =  list(),
                                  default_save_name = 'plotCellProximityGenes') {


  if(!'cpgObject' %in% class(cpgObject)) {
    stop('\n cpgObject needs to be the output from findCellProximityGenes() or findCPG() \n')
  }

  # print, return and save parameters
  show_plot = ifelse(is.na(show_plot), readGiottoInstructions(gobject, param = 'show_plot'), show_plot)
  save_plot = ifelse(is.na(save_plot), readGiottoInstructions(gobject, param = 'save_plot'), save_plot)
  return_plot = ifelse(is.na(return_plot), readGiottoInstructions(gobject, param = 'return_plot'), return_plot)


  ## first filter
  filter_cpg = filterInteractionChangedGenes(cpgObject = cpgObject,
                                             min_cells = min_cells,
                                             min_cells_expr = min_cells_expr,
                                             min_int_cells = min_int_cells,
                                             min_int_cells_expr = min_int_cells_expr,
                                             min_fdr = min_fdr,
                                             min_spat_diff = min_spat_diff,
                                             min_log2_fc = min_log2_fc,
                                             min_zscore = min_zscore,
                                             zscores_column = zscores_column,
                                             direction = direction)

  complete_part = filter_cpg[['CPGscores']]

  ## other parameters
  method = match.arg(method, choices = c('volcano', 'cell_barplot', 'cell-cell', 'cell_sankey', 'heatmap', 'dotplot'))


  # variables
  log2fc = p.adj =  unif_int = N = cell_type = int_cell_type = NULL

  ## create data.table for visualization
  if(method == 'volcano') {

    ## volcanoplot
    pl = ggplot2::ggplot()
    pl = pl + ggplot2::geom_point(data = complete_part, ggplot2::aes(x = log2fc, y = ifelse(is.infinite(-log10(p.adj)), 1000, -log10(p.adj))))
    pl = pl + ggplot2::theme_classic()
    pl = pl + ggplot2::geom_vline(xintercept = 0, linetype = 2)
    pl = pl + ggplot2::labs(x = 'log2 fold-change', y = '-log10(p.adjusted)')


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
    pl <- pl + ggplot2::labs(x = '', y = '# of genes influenced by cell neighborhood')


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
      ggplot2::theme_classic() + ggplot2::labs(x = '', y = '# of genes influenced by cell neighborhood')

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

    changed_genes = complete_part[, .N, by = c('cell_type', 'int_cell_type')]

    changed_genes[, cell_type := factor(cell_type, unique(cell_type))]
    changed_genes[, int_cell_type := factor(int_cell_type, unique(int_cell_type))]

    pl = ggplot2::ggplot()
    pl = pl + ggplot2::theme_classic()
    pl = pl + ggplot2::geom_point(data = changed_genes, ggplot2::aes(x = cell_type, y = int_cell_type, size = N))
    pl = pl + ggplot2::scale_size_continuous(guide=guide_legend(title = '# of ICGs'))
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

    changed_genes = complete_part[, .N, by = c('cell_type', 'int_cell_type')]

    changed_genes[, cell_type := factor(cell_type, unique(cell_type))]
    changed_genes[, int_cell_type := factor(int_cell_type, unique(int_cell_type))]

    changed_genes_d = data.table::dcast.data.table(changed_genes, cell_type~int_cell_type, value.var = 'N', fill = 0)
    changed_genes_m = dt_to_matrix(changed_genes_d)

    col_fun = circlize::colorRamp2(breaks = stats::quantile(log2(changed_genes_m+1)),
                                   colors =  c("white", 'white', "blue", "yellow", "red"))

    heatm = ComplexHeatmap::Heatmap(log2(changed_genes_m+1), col = col_fun,
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


#' @title plotCPG
#' @name plotCPG
#' @description Create visualization for cell proximity gene scores
#' @param gobject giotto object
#' @param cpgObject ICG (interaction changed gene) score object
#' @param method plotting method to use
#' @param min_cells minimum number of source cell type
#' @param min_cells_expr minimum expression level for source cell type
#' @param min_int_cells minimum number of interacting neighbor cell type
#' @param min_int_cells_expr minimum expression level for interacting neighbor cell type
#' @param min_fdr minimum adjusted p-value
#' @param min_spat_diff minimum absolute spatial expression difference
#' @param min_log2_fc minimum log2 fold-change
#' @param min_zscore minimum z-score change
#' @param zscores_column calculate z-scores over cell types or genes
#' @param direction differential expression directions to keep
#' @param cell_color_code vector of colors with cell types as names
#' @param show_plot show plots
#' @param return_plot return plotting object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters from \code{\link{all_plots_save_function}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return plot
#' @export
plotCPG = function(gobject,
                   cpgObject,
                   method = c('volcano', 'cell_barplot', 'cell-cell', 'cell_sankey', 'heatmap', 'dotplot'),
                   min_cells = 5,
                   min_cells_expr = 1,
                   min_int_cells = 3,
                   min_int_cells_expr = 1,
                   min_fdr = 0.05,
                   min_spat_diff = 0.2,
                   min_log2_fc = 0.2,
                   min_zscore = 2,
                   zscores_column = c('cell_type', 'genes'),
                   direction = c('both', 'up', 'down'),
                   cell_color_code = NULL,
                   show_plot = NA,
                   return_plot = NA,
                   save_plot = NA,
                   save_param =  list(),
                   default_save_name = 'plotCPG') {


  plotCellProximityGenes(gobject = gobject,
                         cpgObject = cpgObject,
                         method = method,
                         min_cells = min_cells,
                         min_cells_expr = min_cells_expr,
                         min_int_cells = min_int_cells,
                         min_int_cells_expr = min_int_cells_expr,
                         min_fdr = min_fdr,
                         min_spat_diff = min_spat_diff,
                         min_log2_fc = min_log2_fc,
                         min_zscore = min_zscore,
                         zscores_column = zscores_column,
                         direction = direction,
                         cell_color_code = cell_color_code,
                         show_plot = show_plot,
                         return_plot = return_plot,
                         save_plot = save_plot,
                         save_param =  save_param,
                         default_save_name = default_save_name)


}





#' @title plotInteractionChangedGenes
#' @name plotInteractionChangedGenes
#' @description Create barplot to visualize interaction changed genes
#' @param gobject giotto object
#' @param cpgObject ICG (interaction changed gene) score object
#' @param source_type cell type of the source cell
#' @param source_markers markers for the source cell type
#' @param ICG_genes named character vector of ICG genes
#' @param cell_color_code cell color code for the interacting cell types
#' @param show_plot show plots
#' @param return_plot return plotting object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters from \code{\link{all_plots_save_function}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return plot
#' @export
plotInteractionChangedGenes = function(gobject,
                                       cpgObject,
                                       source_type,
                                       source_markers,
                                       ICG_genes,
                                       cell_color_code = NULL,
                                       show_plot = NA,
                                       return_plot = NA,
                                       save_plot = NA,
                                       save_param =  list(),
                                       default_save_name = 'plotInteractionChangedGenes') {


  # data.table variables
  cell_type = int_cell_type = log2fc = NULL


  if(!'cpgObject' %in% class(cpgObject)) {
    stop('\n cpgObject needs to be the output from findCellProximityGenes() or findCPG() \n')
  }

  CPGscores = cpgObject[['CPGscores']]

  # combine genes
  names(source_markers) = rep('marker', length(source_markers))
  neighbor_types = names(ICG_genes)
  all_genes = c(source_markers, ICG_genes)

  # warning if there are genes selected that are not detected
  detected_genes = unique(CPGscores[['genes']])
  not_detected_genes = all_genes[!all_genes %in% detected_genes]
  if(length(not_detected_genes) > 0) {
    cat('These selected genes are not in the cpgObject: \n',
        not_detected_genes, '\n')
  }

  # data.table set column names
  genes = group = NULL

  tempDT = CPGscores[genes %in% all_genes][cell_type == source_type][int_cell_type %in% neighbor_types]
  tempDT[, genes := factor(genes, levels = all_genes)]
  tempDT[, group := names(all_genes[all_genes == genes]), by = 1:nrow(tempDT)]


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
  pl = pl + ggplot2::geom_bar(data = tempDT, ggplot2::aes(x = genes, y = log2fc, fill = int_cell_type), stat = 'identity', position = ggplot2::position_dodge())
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


#' @title plotICG
#' @name plotICG
#' @description Create barplot to visualize interaction changed genes
#' @param gobject giotto object
#' @param cpgObject ICG (interaction changed gene) score object
#' @param source_type cell type of the source cell
#' @param source_markers markers for the source cell type
#' @param ICG_genes named character vector of ICG genes
#' @param cell_color_code cell color code for the interacting cell types
#' @param show_plot show plots
#' @param return_plot return plotting object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters from \code{\link{all_plots_save_function}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return plot
#' @export
plotICG = function(gobject,
                   cpgObject,
                   source_type,
                   source_markers,
                   ICG_genes,
                   cell_color_code = NULL,
                   show_plot = NA,
                   return_plot = NA,
                   save_plot = NA,
                   save_param =  list(),
                   default_save_name = 'plotICG') {


  plotInteractionChangedGenes(gobject = gobject,
                             cpgObject = cpgObject,
                             source_type = source_type,
                             source_markers = source_markers,
                             ICG_genes = ICG_genes,
                             cell_color_code = cell_color_code,
                             show_plot = show_plot,
                             return_plot = return_plot,
                             save_plot = save_plot,
                             save_param =  save_param,
                             default_save_name = default_save_name)

}






#' @title plotCombineInteractionChangedGenes
#' @name plotCombineInteractionChangedGenes
#' @description Create visualization for combined (pairwise) ICG scores
#' @param gobject giotto object
#' @param combCpgObject ICGscores, output from combineInteractionChangedGenes()
#' @param selected_interactions interactions to show
#' @param selected_gene_to_gene pairwise gene combinations to show
#' @param detail_plot show detailed info in both interacting cell types
#' @param simple_plot show a simplified plot
#' @param simple_plot_facet facet on interactions or genes with simple plot
#' @param facet_scales ggplot facet scales paramter
#' @param facet_ncol ggplot facet ncol parameter
#' @param facet_nrow ggplot facet nrow parameter
#' @param colors vector with two colors to use
#' @param show_plot show plots
#' @param return_plot return plotting object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters from \code{\link{all_plots_save_function}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return ggplot
#' @export
plotCombineInteractionChangedGenes <- function(gobject,
                                               combCpgObject,
                                               selected_interactions = NULL,
                                               selected_gene_to_gene = NULL,
                                               detail_plot = T,
                                               simple_plot = F,
                                               simple_plot_facet = c('interaction', 'genes'),
                                               facet_scales = 'fixed',
                                               facet_ncol = length(selected_gene_to_gene),
                                               facet_nrow = length(selected_interactions),
                                               colors = c('#9932CC', '#FF8C00'),
                                               show_plot = NA,
                                               return_plot = NA,
                                               save_plot = NA,
                                               save_param =  list(),
                                               default_save_name = 'plotCombineICG') {



  ## check validity
  if(!'combCpgObject' %in% class(combCpgObject)) {
    stop('\n combCpgObject needs to be the output from combineCellProximityGenes() or combineCPG() \n')
  }
  combCPGscore = copy(combCpgObject[['combCPGscores']])

  if(is.null(selected_interactions) | is.null(selected_gene_to_gene)) {
    stop('\n You need to provide a selection of cell-cell interactions and genes-genes to plot \n')
  }


  # data.table variables
  unif_gene_gene = unif_int = other_2 = sel_2 = other_1 = sel_1 = cols = NULL


  subDT = combCPGscore[unif_gene_gene %in% selected_gene_to_gene & unif_int %in% selected_interactions]

  # order interactions and gene-to-gene according to input
  subDT[, unif_gene_gene := factor(unif_gene_gene, levels = selected_gene_to_gene)]
  subDT[, unif_int := factor(unif_int, levels = selected_interactions)]

  if(simple_plot == F) {

    pl <- ggplot2::ggplot()
    pl <- pl + ggplot2::theme_bw()

    if(detail_plot == TRUE) {
      pl <- pl + ggplot2::geom_point(data = subDT, aes(x = 0, y = other_2, colour = "other cell expression"),shape = 1)
      pl <- pl + ggplot2::geom_point(data = subDT, aes(x = 0, y = sel_2, colour = "selected cell expression"), shape = 1)
      pl <- pl + ggplot2::geom_point(data = subDT, aes(x = other_1, y = 0, colour = "other cell expression"), shape = 1)
      pl <- pl + ggplot2::geom_point(data = subDT, aes(x = sel_1, y = 0,colour ="selected cell expression"), shape = 1)
    }

    pl <- pl + ggplot2::geom_point(data = subDT, aes(x = other_1, y = other_2, colour = "other cell expression"),size = 2)
    pl <- pl + ggplot2::geom_point(data = subDT, aes(x = sel_1, y = sel_2, colour ="selected cell expression"), size = 2)
    pl <- pl + ggplot2::geom_segment(data = subDT, aes(x = other_1, xend = sel_1,
                                                       y = other_2, yend = sel_2), linetype = 2)
    #pl <- pl + ggplot2::labs(x = 'gene 1 in celltype 1', y = 'gene 2 in celltype 2')
    pl <- pl + ggplot2::labs(x = paste(subDT$genes_1, subDT$cell_type_1, sep = " in ")
                             , y = paste(subDT$genes_2, subDT$cell_type_2, sep = " in "))
    pl <- pl + ggplot2::scale_colour_manual(name="expression source",values = colors)
    pl <- pl + ggplot2::facet_wrap(~unif_gene_gene+unif_int, nrow = facet_nrow, ncol = facet_ncol,
                                   scales = facet_scales)

  }else {

    simple_plot_facet = match.arg(arg = simple_plot_facet, choices = c('interaction', 'genes'))

    if(simple_plot_facet == 'interaction') {
      pl <- ggplot2::ggplot()
      pl <- pl + ggplot2::theme_bw()
      pl <- pl + ggplot2::geom_segment(data = subDT, aes(x = sum(c(other_1, other_2)), xend = sum(c(sel_1, sel_2)),
                                                         y = unif_gene_gene, yend = unif_gene_gene), linetype = 2)
      pl <- pl + ggplot2::geom_point(data = subDT, aes(x = sum(c(other_1, other_2)), y = unif_gene_gene,colour = "other cell expression"))
      pl <- pl + ggplot2::geom_point(data = subDT, aes(x = sum(c(sel_1, sel_2)), y = unif_gene_gene,colour ="selected cell expression"))
      pl <- pl + ggplot2::scale_colour_manual(name="expression source",values=cols)
      pl <- pl + ggplot2::facet_wrap(~unif_int, scales = facet_scales)
      pl <- pl + ggplot2::labs(x = 'interactions', y = 'gene-gene')
    } else {
      pl <- ggplot2::ggplot()
      pl <- pl + ggplot2::theme_bw()
      pl <- pl + ggplot2::geom_segment(data = subDT, aes(x = sum(c(other_1, other_2)), xend = sum(c(sel_1, sel_2)),
                                                         y = unif_int, yend = unif_int), linetype = 2)
      pl <- pl + ggplot2::geom_point(data = subDT, aes(x = sum(c(other_1, other_2)), y = unif_int, colour = "other cell expression"))
      pl <- pl + ggplot2::geom_point(data = subDT, aes(x = sum(c(sel_1, sel_2)), y = unif_int, colour ="selected cell expression"))
      pl <- pl + ggplot2::scale_colour_manual(name="expression source",values=cols)
      pl <- pl + ggplot2::facet_wrap(~unif_gene_gene, scales = facet_scales)
      pl <- pl + ggplot2::labs(x = 'gene-gene', y = 'interactions')
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



#' @title plotCombineCellProximityGenes
#' @description Create visualization for combined (pairwise) ICG scores
#' @inheritDotParams plotCombineInteractionChangedGenes
#' @seealso \code{\link{plotCombineInteractionChangedGenes}}
#' @export
plotCombineCellProximityGenes <- function(...) {

  .Deprecated(new = "plotCombineInteractionChangedGenes")

  plotCombineInteractionChangedGenes(...)

}




#' @title plotCombineICG
#' @name plotCombineICG
#' @description Create visualization for combined (pairwise) ICG scores
#' @param gobject giotto object
#' @param combCpgObject ICGscores, output from combineInteractionChangedGenes()
#' @param selected_interactions interactions to show
#' @param selected_gene_to_gene pairwise gene combinations to show
#' @param detail_plot show detailed info in both interacting cell types
#' @param simple_plot show a simplified plot
#' @param simple_plot_facet facet on interactions or genes with simple plot
#' @param facet_scales ggplot facet scales paramter
#' @param facet_ncol ggplot facet ncol parameter
#' @param facet_nrow ggplot facet nrow parameter
#' @param colors vector with two colors to use
#' @param show_plot show plots
#' @param return_plot return plotting object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters from \code{\link{all_plots_save_function}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return ggplot
#' @export
plotCombineICG <- function(gobject,
                           combCpgObject,
                           selected_interactions = NULL,
                           selected_gene_to_gene = NULL,
                           detail_plot = T,
                           simple_plot = F,
                           simple_plot_facet = c('interaction', 'genes'),
                           facet_scales = 'fixed',
                           facet_ncol = length(selected_gene_to_gene),
                           facet_nrow = length(selected_interactions),
                           colors = c('#9932CC', '#FF8C00'),
                           show_plot = NA,
                           return_plot = NA,
                           save_plot = NA,
                           save_param =  list(),
                           default_save_name = 'plotCombineICG') {


  plotCombineCellProximityGenes(gobject = gobject,
                                combCpgObject = combCpgObject,
                                selected_interactions = selected_interactions,
                                selected_gene_to_gene = selected_gene_to_gene,
                                detail_plot = detail_plot,
                                simple_plot = simple_plot,
                                simple_plot_facet = simple_plot_facet,
                                facet_scales = facet_scales,
                                facet_ncol = facet_ncol,
                                facet_nrow = facet_nrow,
                                colors = colors,
                                show_plot = show_plot,
                                return_plot = return_plot,
                                save_plot = save_plot,
                                save_param =  save_param,
                                default_save_name = default_save_name)


}



#' @title plotCombineCPG
#' @description Create visualization for combined (pairwise) ICG scores
#' @inheritDotParams plotCombineICG
#' @seealso \code{\link{plotCombineICG}}
#' @export
plotCombineCPG <- function(...) {

  .Deprecated(new = "plotCombineICG")

  plotCombineICG(...)

}




# * ####
# cell communication plots ####

#' @title plotCombineCellCellCommunication
#' @name plotCombineCellCellCommunication
#' @description Create visualization for combined (pairwise) cell proximity gene scores
#' @param gobject giotto object
#' @param combCCcom combined communcation scores, output from combCCcom()
#' @param selected_LR selected ligand-receptor pair
#' @param selected_cell_LR selected cell-cell interaction pair for ligand-receptor pair
#' @param detail_plot show detailed info in both interacting cell types
#' @param simple_plot show a simplified plot
#' @param simple_plot_facet facet on interactions or genes with simple plot
#' @param facet_scales ggplot facet scales paramter
#' @param facet_ncol ggplot facet ncol parameter
#' @param facet_nrow ggplot facet nrow parameter
#' @param colors vector with two colors to use
#' @param show_plot show plots
#' @param return_plot return plotting object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters from \code{\link{all_plots_save_function}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return ggplot
#' @export
plotCombineCellCellCommunication <- function(gobject,
                                             combCCcom,
                                             selected_LR = NULL,
                                             selected_cell_LR = NULL,
                                             detail_plot = T,
                                             simple_plot = F,
                                             simple_plot_facet = c('interaction', 'genes'),
                                             facet_scales = 'fixed',
                                             facet_ncol = length(selected_LR),
                                             facet_nrow = length(selected_cell_LR),
                                             colors = c('#9932CC', '#FF8C00'),
                                             show_plot = NA,
                                             return_plot = NA,
                                             save_plot = NA,
                                             save_param =  list(),
                                             default_save_name = 'plotCombineCellCellCommunication') {



  # data.table variables
  LR_comb = LR_cell_comb = lig_expr = lig_expr_spat = rec_expr = rec_expr_spat = LR_expr = LR_expr_spat = NULL

  ## check validity
  if(is.null(selected_cell_LR) | is.null(selected_LR)) {
    stop('\n You need to provide a selection of cell-cell interactions and genes-genes to plot \n')
  }

  subDT = combCCcom[LR_comb %in% selected_LR & LR_cell_comb %in% selected_cell_LR]

  # order interactions and gene-to-gene according to input
  subDT[, LR_comb := factor(LR_comb, levels = selected_LR)]
  subDT[, LR_cell_comb := factor(LR_cell_comb, levels = selected_cell_LR)]

  if(simple_plot == F) {

    pl <- ggplot2::ggplot()
    pl <- pl + ggplot2::theme_bw()

    if(detail_plot == TRUE) {
      pl <- pl + ggplot2::geom_point(data = subDT, aes(x = 0, y = lig_expr, colour = "overall cell expression"),shape = 1)
      pl <- pl + ggplot2::geom_point(data = subDT, aes(x = 0, y = lig_expr_spat, colour = "spatial cell expression"), shape = 1)
      pl <- pl + ggplot2::geom_point(data = subDT, aes(x = rec_expr, y = 0, colour = "overall cell expression"), shape = 1)
      pl <- pl + ggplot2::geom_point(data = subDT, aes(x = rec_expr_spat, y = 0,colour ="spatial cell expression"), shape = 1)
    }

    pl <- pl + ggplot2::geom_point(data = subDT, aes(x = rec_expr, y = lig_expr, colour = "overall cell expression"),size = 2)
    pl <- pl + ggplot2::geom_point(data = subDT, aes(x = rec_expr_spat, y = lig_expr_spat, colour ="spatial cell expression"), size = 2)
    pl <- pl + ggplot2::geom_segment(data = subDT, aes(x = rec_expr, xend = rec_expr_spat,
                                                       y = lig_expr, yend = lig_expr_spat), linetype = 2)
    #pl <- pl + ggplot2::labs(x = 'gene 1 in celltype 1', y = 'gene 2 in celltype 2')
    pl <- pl + ggplot2::labs(x = paste(subDT$receptor, subDT$rec_cell_type, sep = " in ")
                             , y = paste(subDT$ligand, subDT$lig_cell_type, sep = " in "))
    pl <- pl + ggplot2::scale_colour_manual(name="expression source",values = colors)
    pl <- pl + ggplot2::facet_wrap(~LR_comb+LR_cell_comb, nrow = facet_nrow, ncol = facet_ncol,
                                   scales = facet_scales)

  }else {

    simple_plot_facet = match.arg(arg = simple_plot_facet, choices = c('interaction', 'genes'))

    if(simple_plot_facet == 'interaction') {
      pl <- ggplot2::ggplot()
      pl <- pl + ggplot2::theme_bw()
      pl <- pl + ggplot2::geom_segment(data = subDT, aes(x = LR_expr, xend = LR_expr_spat,
                                                         y = LR_comb, yend = LR_comb), linetype = 2)
      pl <- pl + ggplot2::geom_point(data = subDT, aes(x = LR_expr, y = LR_comb, colour = "overall cell expression"))
      pl <- pl + ggplot2::geom_point(data = subDT, aes(x = LR_expr_spat, y = LR_comb,colour ="spatial cell expression"))
      pl <- pl + ggplot2::scale_colour_manual(name="expression source",values=colors)
      pl <- pl + ggplot2::facet_wrap(~LR_cell_comb, scales = 'fixed')
      pl <- pl + ggplot2::labs(x = 'interactions', y = 'gene-gene')
      pl

    } else {
      pl <- ggplot2::ggplot()
      pl <- pl + ggplot2::theme_bw()
      pl <- pl + ggplot2::geom_segment(data = subDT, aes(x = LR_expr, xend = LR_expr_spat,
                                                         y = LR_cell_comb, yend = LR_cell_comb), linetype = 2)
      pl <- pl + ggplot2::geom_point(data = subDT, aes(x = LR_expr, y = LR_cell_comb, colour = "overall cell expression"))
      pl <- pl + ggplot2::geom_point(data = subDT, aes(x = LR_expr_spat, y = LR_cell_comb, colour ="spatial cell expression"))
      pl <- pl + ggplot2::scale_colour_manual(name="expression source",values=colors)
      pl <- pl + ggplot2::facet_wrap(~LR_comb, scales = facet_scales)
      pl <- pl + ggplot2::labs(x = 'gene-gene', y = 'interactions')
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



#' @title plotCombineCCcom
#' @name plotCombineCCcom
#' @description Create visualization for combined (pairwise) cell proximity gene scores
#' @param gobject giotto object
#' @param combCCcom combined communcation scores, output from combCCcom()
#' @param selected_LR selected ligand-receptor pair
#' @param selected_cell_LR selected cell-cell interaction pair for ligand-receptor pair
#' @param detail_plot show detailed info in both interacting cell types
#' @param simple_plot show a simplified plot
#' @param simple_plot_facet facet on interactions or genes with simple plot
#' @param facet_scales ggplot facet scales paramter
#' @param facet_ncol ggplot facet ncol parameter
#' @param facet_nrow ggplot facet nrow parameter
#' @param colors vector with two colors to use
#' @param show_plot show plots
#' @param return_plot return plotting object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters from \code{\link{all_plots_save_function}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return ggplot
#' @export
plotCombineCCcom = function(gobject,
                            combCCcom,
                            selected_LR = NULL,
                            selected_cell_LR = NULL,
                            detail_plot = T,
                            simple_plot = F,
                            simple_plot_facet = c('interaction', 'genes'),
                            facet_scales = 'fixed',
                            facet_ncol = length(selected_LR),
                            facet_nrow = length(selected_cell_LR),
                            colors = c('#9932CC', '#FF8C00'),
                            show_plot = NA,
                            return_plot = NA,
                            save_plot = NA,
                            save_param =  list(),
                            default_save_name = 'plotCombineCCcom') {


  plotCombineCellCellCommunication(gobject = gobject,
                                   combCCcom = combCCcom,
                                   selected_LR = selected_LR,
                                   selected_cell_LR = selected_cell_LR,
                                   detail_plot = detail_plot,
                                   simple_plot = simple_plot,
                                   simple_plot_facet = simple_plot_facet,
                                   facet_scales = facet_scales,
                                   facet_ncol = facet_ncol,
                                   facet_nrow = facet_nrow,
                                   colors = colors,
                                   show_plot = show_plot,
                                   return_plot = return_plot,
                                   save_plot = save_plot,
                                   save_param =  save_param,
                                   default_save_name = default_save_name)

}



#' @title plotCCcomHeatmap
#' @name plotCCcomHeatmap
#' @description Plots heatmap for ligand-receptor communication scores in cell-cell interactions
#' @param gobject giotto object
#' @param comScores communinication scores from \code{\link{exprCellCellcom}} or \code{\link{spatCellCellcom}}
#' @param selected_LR selected ligand-receptor combinations
#' @param selected_cell_LR selected cell-cell combinations for ligand-receptor combinations
#' @param show_LR_names show ligand-receptor names
#' @param show_cell_LR_names show cell-cell names
#' @param show values to show on heatmap
#' @param cor_method correlation method used for clustering
#' @param aggl_method agglomeration method used by hclust
#' @param show_plot show plots
#' @param return_plot return plotting object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters from \code{\link{all_plots_save_function}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return ggplot
#' @export
plotCCcomHeatmap = function(gobject,
                            comScores,
                            selected_LR = NULL,
                            selected_cell_LR = NULL,
                            show_LR_names = TRUE,
                            show_cell_LR_names = TRUE,
                            show = c('PI', 'LR_expr', 'log2fc'),
                            cor_method = c("pearson", "kendall", "spearman"),
                            aggl_method = c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid"),
                            show_plot = NA,
                            return_plot = NA,
                            save_plot = NA,
                            save_param =  list(),
                            default_save_name = 'plotCCcomHeatmap') {


  # get parameters
  cor_method = match.arg(cor_method, choices = c("pearson", "kendall", "spearman"))
  aggl_method = match.arg(aggl_method, choices = c("ward.D", "ward.D2", "single", "complete",
                                                   "average", "mcquitty", "median", "centroid"))


  # data.table variables
  LR_comb = LR_cell_comb = NULL

  # plot method
  if(!is.null(selected_LR) & !is.null(selected_cell_LR)) {
    selDT = comScores[LR_comb %in% selected_LR & LR_cell_comb %in% selected_cell_LR]
  } else if(!is.null(selected_LR)) {
    selDT = comScores[LR_comb %in% selected_LR]
  } else if(!is.null(selected_cell_LR)) {
    selDT = comScores[LR_cell_comb %in% selected_cell_LR]
  } else {
    selDT = comScores
  }

  # creat matrix
  show = match.arg(show, choices = c('PI', 'LR_expr', 'log2fc'))
  selDT_d = data.table::dcast.data.table(selDT, LR_cell_comb~LR_comb, value.var = show, fill = 0)
  selDT_m = dt_to_matrix(selDT_d)

  ## cells
  corclus_cells_dist = stats::as.dist(1-cor_giotto(x = t_giotto(selDT_m), method = cor_method))
  hclusters_cells = stats::hclust(d = corclus_cells_dist, method = aggl_method)
  clus_names = rownames(selDT_m)
  names(clus_names) = 1:length(clus_names)
  clus_sort_names = clus_names[hclusters_cells$order]
  selDT[, LR_cell_comb := factor(LR_cell_comb, clus_sort_names)]

  ## genes
  corclus_genes_dist = stats::as.dist(1-cor_giotto(x = selDT_m, method = cor_method))
  hclusters_genes = stats::hclust(d = corclus_genes_dist, method = aggl_method)
  clus_names = colnames(selDT_m)
  names(clus_names) = 1:length(clus_names)
  clus_sort_names = clus_names[hclusters_genes$order]
  selDT[, LR_comb := factor(LR_comb, clus_sort_names)]



  pl = ggplot2::ggplot()
  pl = pl + ggplot2::geom_raster(data = selDT, aes_string(x = 'LR_cell_comb',
                                                  y = 'LR_comb', fill = show))

  pl = pl + ggplot2::theme_classic() + ggplot2::theme(axis.text.x = element_blank(),
                                     axis.ticks = element_blank(),
                                     axis.text.y = element_blank())
  if(show_LR_names == TRUE) pl <- pl + ggplot2::theme(axis.text.y = element_text(),
                                             axis.ticks.y = element_line())
  if(show_cell_LR_names == TRUE) pl <- pl + ggplot2::theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1),
                                                  axis.ticks.x = element_line())
  pl = pl + ggplot2::scale_fill_gradientn(colours = c('darkblue', 'blue', 'white', 'red', 'darkred'))
  pl = pl + ggplot2::labs(x = 'cell-cell', y = 'ligand-receptor')


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



#' @title plotCCcomDotplot
#' @name plotCCcomDotplot
#' @description Plots dotplot for ligand-receptor communication scores in cell-cell interactions
#' @param gobject giotto object
#' @param comScores communinication scores from \code{\link{exprCellCellcom}} or \code{\link{spatCellCellcom}}
#' @param selected_LR selected ligand-receptor combinations
#' @param selected_cell_LR selected cell-cell combinations for ligand-receptor combinations
#' @param show_LR_names show ligand-receptor names
#' @param show_cell_LR_names show cell-cell names
#' @param cluster_on values to use for clustering of cell-cell and ligand-receptor pairs
#' @param cor_method correlation method used for clustering
#' @param aggl_method agglomeration method used by hclust
#' @param show_plot show plots
#' @param return_plot return plotting object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters from \code{\link{all_plots_save_function}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return ggplot
#' @export
plotCCcomDotplot = function(gobject,
                            comScores,
                            selected_LR = NULL,
                            selected_cell_LR = NULL,
                            show_LR_names = TRUE,
                            show_cell_LR_names = TRUE,
                            cluster_on = c('PI', 'LR_expr', 'log2fc'),
                            cor_method = c("pearson", "kendall", "spearman"),
                            aggl_method = c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid"),
                            show_plot = NA,
                            return_plot = NA,
                            save_plot = NA,
                            save_param =  list(),
                            default_save_name = 'plotCCcomDotplot') {

  # get parameters
  cor_method = match.arg(cor_method, choices = c("pearson", "kendall", "spearman"))
  aggl_method = match.arg(aggl_method, choices = c("ward.D", "ward.D2", "single", "complete",
                                                   "average", "mcquitty", "median", "centroid"))


  # data.table variables
  LR_comb = LR_cell_comb = sd = NULL

  # plot method
  if(!is.null(selected_LR) & !is.null(selected_cell_LR)) {
    selDT = comScores[LR_comb %in% selected_LR & LR_cell_comb %in% selected_cell_LR]
  } else if(!is.null(selected_LR)) {
    selDT = comScores[LR_comb %in% selected_LR]
  } else if(!is.null(selected_cell_LR)) {
    selDT = comScores[LR_cell_comb %in% selected_cell_LR]
  } else {
    selDT = comScores
  }

  # creat matrix
  cluster_on = match.arg(cluster_on, choices = c('PI', 'LR_expr', 'log2fc'))
  selDT_d = data.table::dcast.data.table(selDT, LR_cell_comb~LR_comb, value.var = cluster_on, fill = 0)
  selDT_m = dt_to_matrix(selDT_d)

  # remove zero variance
  sd_rows = apply(selDT_m, 1, sd)
  sd_rows_zero = names(sd_rows[sd_rows == 0])
  if(length(sd_rows_zero) > 0) selDT_m = selDT_m[!rownames(selDT_m) %in% sd_rows_zero, ]

  sd_cols = apply(selDT_m, 2, sd)
  sd_cols_zero = names(sd_cols[sd_cols == 0])
  if(length(sd_cols_zero) > 0) selDT_m = selDT_m[, !colnames(selDT_m) %in% sd_cols_zero]



  ## cells
  corclus_cells_dist = stats::as.dist(1-cor_giotto(x = t_giotto(selDT_m), method = cor_method))
  hclusters_cells = stats::hclust(d = corclus_cells_dist, method = aggl_method)
  clus_names = rownames(selDT_m)
  names(clus_names) = 1:length(clus_names)
  clus_sort_names = clus_names[hclusters_cells$order]
  selDT[, LR_cell_comb := factor(LR_cell_comb, clus_sort_names)]

  ## genes
  corclus_genes_dist = stats::as.dist(1-cor_giotto(x = selDT_m, method = cor_method))
  hclusters_genes = stats::hclust(d = corclus_genes_dist, method = aggl_method)
  clus_names = colnames(selDT_m)
  names(clus_names) = 1:length(clus_names)
  clus_sort_names = clus_names[hclusters_genes$order]
  selDT[, LR_comb := factor(LR_comb, clus_sort_names)]



  pl = ggplot2::ggplot()
  pl = pl + ggplot2::geom_point(data = selDT, aes_string(x = 'LR_cell_comb',
                                                 y = 'LR_comb', size = 'pvalue', color = 'log2fc'))
  pl = pl + ggplot2::theme_classic()
  if(show_LR_names == TRUE) pl = pl + ggplot2::theme(axis.text.y = element_text(),
                                             axis.ticks.y = element_line())
  if(show_cell_LR_names == TRUE) pl = pl + ggplot2::theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1),
                                                  axis.ticks.x = element_line())
  pl = pl + ggplot2::scale_fill_gradientn(colours = c('darkblue', 'blue', 'white', 'red', 'darkred'))
  pl = pl + ggplot2::scale_size_continuous(range = c(5, 0.5)) + scale_color_gradientn(colours = c('darkblue', 'blue', 'white', 'red', 'darkred'))
  pl = pl + ggplot2::labs(x = 'cell-cell', y = 'ligand-receptor')

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



#' @title plotRankSpatvsExpr
#' @name plotRankSpatvsExpr
#' @description Plots dotplot to compare ligand-receptor rankings from spatial and expression information
#' @param gobject giotto object
#' @param combCC combined communinication scores from \code{\link{combCCcom}}
#' @param expr_rnk_column column with expression rank information to use
#' @param spat_rnk_column column with spatial rank information to use
#' @param midpoint midpoint of colors
#' @param size_range size ranges of dotplot
#' @param xlims x-limits, numerical vector of 2
#' @param ylims y-limits, numerical vector of 2
#' @param selected_ranks numerical vector, will be used to print out the percentage of top spatial ranks are recovered
#' @param show_plot show plots
#' @param return_plot return plotting object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters from \code{\link{all_plots_save_function}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return ggplot
#' @export
plotRankSpatvsExpr = function(gobject,
                              combCC,
                              expr_rnk_column = 'LR_expr_rnk',
                              spat_rnk_column = 'LR_spat_rnk',
                              midpoint = 10,
                              size_range = c(0.01, 1.5),
                              xlims = NULL,
                              ylims = NULL,
                              selected_ranks = c(1, 10, 20),
                              show_plot = NA,
                              return_plot = NA,
                              save_plot = NA,
                              save_param =  list(),
                              default_save_name = 'plotRankSpatvsExpr') {


  # data.table variables
  spt_rank = variable = value = NULL

  total_rnks = max(unique(combCC[[expr_rnk_column]]))

  rnk_list = list()
  spt_list = list()
  for(rnk in 1:total_rnks) {

    mytab = table(cut(sort(combCC[get(expr_rnk_column) == rnk][[spat_rnk_column]]), breaks = seq(0,total_rnks,1), labels = c(1:total_rnks)))
    rnk_list[[rnk]] = mytab
    spt_list[[rnk]] = names(mytab)
  }

  rnk_res = as.data.table(do.call('rbind', rnk_list))
  rnk_res[, spt_rank := 1:total_rnks]

  rnk_res_m = melt.data.table(rnk_res, id.vars = 'spt_rank')
  rnk_res_m[, spt_rank := as.numeric(spt_rank)]
  rnk_res_m[, variable := as.numeric(variable)]

  rnk_res_m[, diff := variable - spt_rank]

  for(i in selected_ranks) {
    perc_recovered = 100*(sum(rnk_res_m[abs(diff) < i]$value)/sum(rnk_res_m$value))
    cat('for top ', i, ' expression ranks, you recover ', round(perc_recovered, 2), '% of the highest spatial rank \n')
  }


  # full plot
  pl = ggplot2::ggplot()
  pl = pl + ggplot2::theme_classic() + ggplot2::theme(axis.text = element_blank())
  pl = pl + ggplot2::geom_point(data = rnk_res_m, ggplot2::aes(x = variable, y = spt_rank, size = value, color = value))
  pl = pl + ggplot2::scale_color_gradient2(low = 'blue', mid = 'yellow', high = 'red', midpoint = midpoint, guide = guide_legend(title = ''))
  pl = pl + ggplot2::scale_size_continuous(range = size_range, guide = F)
  pl = pl + ggplot2::labs(x = 'expression rank', y = 'spatial rank')

  if(!is.null(xlims)) {
    pl = pl + xlim(xlims)
  }

  if(!is.null(ylims)) {
    pl = pl + ylim(ylims)
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




#' @title plotRecovery_sub
#' @name plotRecovery_sub
#' @description Plots recovery plot to compare ligand-receptor rankings from spatial and expression information
#' @param combCC combined communinication scores from \code{\link{combCCcom}}
#' @param first_col first column to use
#' @param second_col second column to use
plotRecovery_sub = function(combCC,
                            first_col = 'LR_expr_rnk',
                            second_col = 'LR_spat_rnk') {


  # data.table variables
  concord = perc = not_concord = secondrank = secondrank_perc = NULL

  mergeDT_filt = combCC[get(first_col) == 1]

  mymat = matrix(data = NA, nrow = max(combCC[[second_col]]), ncol = 2)
  for(i in 1:max(combCC[[second_col]])) {

    mergeDT_filt[, concord := ifelse(get(second_col) <= i, 'yes', 'no')]
    mytable = table(mergeDT_filt$concord)

    matching = mytable['yes']
    if(is.na(matching)) matching = 0
    mymat[i, 1] = matching

    non_matching = mytable['no']
    if(is.na(non_matching)) non_matching = 0
    mymat[i, 2] = non_matching

  }

  mymatDT = data.table::as.data.table(mymat); colnames(mymatDT) = c('concord', 'not_concord')
  mymatDT[, perc := 100*(concord / (concord+not_concord))]
  mymatDT[, secondrank := 1:nrow(mymatDT)]
  mymatDT[, secondrank_perc := (secondrank/max(secondrank))*100]

  # percentage explained
  perc_explained = mymatDT[, sum(perc)]/(100*nrow(mymat))
  cat('percentage explained = ', perc_explained)


  pl = ggplot2::ggplot()
  pl = pl + ggplot2::theme_classic()
  pl = pl + ggplot2::geom_point(data = mymatDT, aes(x = secondrank_perc, y = perc))
  pl = pl + ggplot2::scale_x_continuous(expand = c(0,0), limits = c(0,100))
  pl = pl + ggplot2::scale_y_continuous(expand = c(0,0), limits = c(0, 100))
  pl = pl + ggplot2::geom_abline(slope = 1, intercept = 0, color = 'blue')

  return(pl)

}





#' @title plotRecovery
#' @name plotRecovery
#' @description Plots recovery plot to compare ligand-receptor rankings from spatial and expression information
#' @param gobject giotto object
#' @param combCC combined communinication scores from \code{\link{combCCcom}}
#' @param expr_rnk_column column with expression rank information to use
#' @param spat_rnk_column column with spatial rank information to use
#' @param ground_truth what to consider as ground truth (default: spatial)
#' @param show_plot show plots
#' @param return_plot return plotting object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters from \code{\link{all_plots_save_function}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return ggplot
#' @export
plotRecovery = function(gobject,
                        combCC,
                        expr_rnk_column = 'exprPI_rnk',
                        spat_rnk_column = 'spatPI_rnk',
                        ground_truth = c('spatial', 'expression'),
                        show_plot = NA,
                        return_plot = NA,
                        save_plot = NA,
                        save_param =  list(),
                        default_save_name = 'plotRecovery') {

  ground_truth = match.arg(ground_truth, choices = c('spatial', 'expression'))


  if(ground_truth == 'spatial') {

    pl = plotRecovery_sub(combCC = combCC,
                          first_col = spat_rnk_column,
                          second_col = expr_rnk_column)
    pl = pl + ggplot2::labs(x = '% expression rank included', y = '% highest spatial rank recovered')

  } else if(ground_truth == 'expression') {

    pl = plotRecovery_sub(combCC = combCC,
                          first_col = expr_rnk_column,
                          second_col = spat_rnk_column)
    pl = pl + ggplot2::labs(x = '% spatial rank included', y = '% highest expression rank recovered')

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









# * ####
# cell proximity spatplots ####

#' @title cellProximitySpatPlot2D
#' @name cellProximitySpatPlot2D
#' @description Visualize 2D cell-cell interactions according to spatial coordinates in ggplot mode
#' @param gobject giotto object
#' @param interaction_name cell-cell interaction name
#' @param cluster_column cluster column with cell clusters
#' @param sdimx x-axis dimension name (default = 'sdimx')
#' @param sdimy y-axis dimension name (default = 'sdimy')
#' @param cell_color color for cells (see details)
#' @param cell_color_code named vector with colors
#' @param color_as_factor convert color column to factor
#' @param show_other_cells decide if show cells not in network
#' @param show_network show spatial network of selected cells
#' @param show_other_network show spatial network of not selected cells
#' @param network_color color of spatial network
#' @param spatial_network_name name of spatial network to use
#' @param show_grid show spatial grid
#' @param grid_color color of spatial grid
#' @param spatial_grid_name name of spatial grid to use
#' @param coord_fix_ratio fix ratio between x and y-axis
#' @param show_legend show legend
#' @param point_size_select size of selected points
#' @param point_select_border_col border color of selected points
#' @param point_select_border_stroke stroke size of selected points
#' @param point_size_other size of other points
#' @param point_alpha_other opacity of other points
#' @param point_other_border_col border color of other points
#' @param point_other_border_stroke stroke size of other points
#' @param show_plot show plots
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters from \code{\link{all_plots_save_function}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return ggplot
#' @details Description of parameters.
#' @export
cellProximitySpatPlot2D <- function(gobject,
                                    interaction_name = NULL,
                                    cluster_column = NULL,
                                    sdimx = 'sdimx',
                                    sdimy = 'sdimy',
                                    cell_color = NULL,
                                    cell_color_code = NULL,
                                    color_as_factor = T,
                                    show_other_cells = F,
                                    show_network = F,
                                    show_other_network = F,
                                    network_color = NULL,
                                    spatial_network_name = 'Delaunay_network',
                                    show_grid = F,
                                    grid_color = NULL,
                                    spatial_grid_name = 'spatial_grid',
                                    coord_fix_ratio = 1,
                                    show_legend = T,
                                    point_size_select = 2,
                                    point_select_border_col = 'black',
                                    point_select_border_stroke = 0.05,
                                    point_size_other = 1,
                                    point_alpha_other = 0.3,
                                    point_other_border_col = 'lightgrey',
                                    point_other_border_stroke = 0.01,
                                    show_plot = NA,
                                    return_plot = NA,
                                    save_plot = NA,
                                    save_param =  list(),
                                    default_save_name = 'cellProximitySpatPlot2D') {
  if(is.null(interaction_name)) {
    stop('\n you need to specific at least one interaction name, run cellProximityEnrichment \n')
  }


  cell_locations  = gobject@spatial_locs
  spatial_grid    = gobject@spatial_grid[[spatial_grid_name]]
  cell_metadata   = gobject@cell_metadata

  spatial_network = annotateSpatialNetwork(gobject = gobject,
                                           spatial_network_name = spatial_network_name,
                                           cluster_column = cluster_column)


  # data.table variables
  unified_int = sdimx_begin = sdimy_begin = sdimx_end = sdimy_end = x_start = x_end = y_start = y_end = cell_ID = NULL

  cell_IDs_to_keep = unique(c(spatial_network[unified_int %in% interaction_name]$to,
                              spatial_network[unified_int %in% interaction_name]$from))

  #print(cell_IDs_to_keep)

  if(show_other_cells){
    CellType <- strsplit(interaction_name,"--")
    all_cell_IDs = cell_metadata[cell_metadata[[cluster_column]] == CellType[[1]][1] |
                                   cell_metadata[[cluster_column]] == CellType[[1]][2],]$cell_ID
    other_cell_IDs <- setdiff(all_cell_IDs, cell_IDs_to_keep)
  }


  # annotated cell data
  if(nrow(cell_metadata) == 0) {
    cell_locations_metadata = cell_locations
  } else {
    cell_locations_metadata <- merge(cell_locations, cell_metadata,by = "cell_ID")
  }


  # first 2 dimensions need to be defined
  if(is.null(sdimx) | is.null(sdimy)) {
    cat('first and second dimenion need to be defined, default is first 2 \n')
    sdimx = 'sdimx'
    sdimy = 'sdimy'
  }

  pl <- ggplot2::ggplot()
  pl <- pl + ggplot2::theme_classic()

  if(!is.null(spatial_network) & show_network == TRUE) {
    if(is.null(network_color)) network_color = 'red'
    if(show_other_network){
      pl <- pl + ggplot2::geom_segment(data = spatial_network[!unified_int %in% interaction_name],
                                       aes(x = sdimx_begin, y = sdimy_begin, xend = sdimx_end, yend = sdimy_end),
                                       color = 'lightgrey', size = 0.5, alpha = 0.5)
    }
    pl <- pl + ggplot2::geom_segment(data = spatial_network[unified_int %in% interaction_name],
                                     aes(x = sdimx_begin, y = sdimy_begin, xend = sdimx_end, yend = sdimy_end),
                                     color = network_color, size = 0.5, alpha = 0.5)
  }

  if(!is.null(spatial_grid) & show_grid == TRUE) {
    if(is.null(grid_color)) grid_color = 'black'
    pl <- pl + ggplot2::geom_rect(data = spatial_grid, aes(xmin = x_start, xmax = x_end, ymin = y_start, ymax = y_end),
                                  color = grid_color, fill = NA)
  }

  # cell color default
  if(is.null(cell_color)) {
    cell_color = 'lightblue'
    pl <- pl + ggplot2::geom_point(data = cell_locations[!cell_ID %in% cell_IDs_to_keep], aes_string(x = sdimx, y = sdimy),
                                   show.legend = show_legend, shape = 21, fill = 'lightgrey', size = point_size_other)
    pl <- pl + ggplot2::geom_point(data = cell_locations[cell_ID %in% cell_IDs_to_keep], aes_string(x = sdimx, y = sdimy),
                                   show.legend = show_legend, shape = 21, fill = cell_color, size = point_size_select)
    if(show_other_cells){
      pl <- pl + ggplot2::geom_point(data = cell_locations[cell_ID %in% other_cell_IDs], aes_string(x = sdimx, y = sdimy),
                                     show.legend = show_legend, shape = 21, fill = cell_color, alpha = point_alpha_other,
                                     size = point_size_select * 0.5)
    }
  }
  else if (is.character(cell_color)) {
    if(cell_color %in% colnames(cell_locations_metadata)) {

      if(color_as_factor == TRUE) {
        factor_data = factor(cell_locations_metadata[[cell_color]])
        cell_locations_metadata[[cell_color]] <- factor_data
      }

      pl <- pl + ggplot2::geom_point(data = cell_locations_metadata[!cell_ID %in% cell_IDs_to_keep], aes_string(x = sdimx, y = sdimy),
                                     fill = 'lightgrey', shape = 21, size = point_size_other,
                                     color = point_other_border_col, stroke = point_other_border_stroke)
      pl <- pl + ggplot2::geom_point(data = cell_locations_metadata[cell_ID %in% cell_IDs_to_keep], aes_string(x = sdimx, y = sdimy, fill = cell_color),
                                     show.legend = show_legend, shape = 21, size = point_size_select,
                                     color = point_select_border_col, stroke = point_select_border_stroke)
      if(show_other_cells){
        pl <- pl + ggplot2::geom_point(data = cell_locations_metadata[cell_ID %in% other_cell_IDs], aes_string(x = sdimx, y = sdimy,fill = cell_color),
                                       show.legend = show_legend, shape = 21, alpha = point_alpha_other,
                                       size = point_size_select * 0.5)
      }



      if(!is.null(cell_color_code)) {
        pl <- pl + ggplot2::scale_fill_manual(values = cell_color_code)
      } else if(color_as_factor == T) {
        number_colors = length(unique(factor_data))
        cell_color_code = getDistinctColors(n = number_colors)
        names(cell_color_code) = unique(factor_data)
        pl <- pl + ggplot2::scale_fill_manual(values = cell_color_code)
      } else if(color_as_factor == F){
        pl <- pl + ggplot2::scale_fill_gradient(low = 'blue', high = 'red')
      }

    } else {
      pl <- pl + ggplot2::geom_point(data = cell_locations_metadata[!cell_ID %in% cell_IDs_to_keep], aes_string(x = sdimx, y = sdimy),
                                     show.legend = show_legend, shape = 21, fill = 'lightgrey', size = point_size_other,
                                     color = point_other_border_col, stroke = point_other_border_stroke)
      pl <- pl + ggplot2::geom_point(data = cell_locations_metadata[cell_ID %in% cell_IDs_to_keep], aes_string(x = sdimx, y = sdimy),
                                     show.legend = show_legend, shape = 21, fill = cell_color, size = point_size_select,
                                     color = point_select_border_col, stroke = point_select_border_stroke)
    }

  }

  pl <- pl + ggplot2::theme_bw() + ggplot2::theme(plot.title = element_text(hjust = 0.5),
                                                  legend.title = element_text(size = 10),
                                                  legend.text = element_text(size = 10))

  # fix coord ratio
  if(!is.null(coord_fix_ratio)) {
    pl <- pl + ggplot2::coord_fixed(ratio = coord_fix_ratio)
  }

  pl <- pl + ggplot2::labs(x = 'x coordinates', y = 'y coordinates')



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



#' @title cellProximitySpatPlot
#' @name cellProximitySpatPlot
#' @description Visualize 2D cell-cell interactions according to spatial coordinates in ggplot mode
#' @param gobject giotto object
#' @inheritDotParams cellProximitySpatPlot2D -gobject
#' @return ggplot
#' @details Description of parameters.
#' @export
#' @seealso  \code{\link{cellProximitySpatPlot2D}} and \code{\link{cellProximitySpatPlot3D}} for 3D
cellProximitySpatPlot = function(gobject, ...) {

  cellProximitySpatPlot2D(gobject = gobject, ...)

}


#' @title cellProximitySpatPlot2D
#' @name cellProximitySpatPlot3D
#' @description Visualize 3D cell-cell interactions according to spatial coordinates in plotly mode
#' @param gobject giotto object
#' @param interaction_name cell-cell interaction name
#' @param cluster_column cluster column with cell clusters
#' @param sdimx x-axis dimension name (default = 'sdimx')
#' @param sdimy y-axis dimension name (default = 'sdimy')
#' @param sdimz z-axis dimension name (default = 'sdimz')
#' @param cell_color color for cells (see details)
#' @param cell_color_code named vector with colors
#' @param color_as_factor convert color column to factor
#' @param show_other_cells decide if show cells not in network
#' @param show_network show spatial network of selected cells
#' @param show_other_network show spatial network of not selected cells
#' @param network_color color of spatial network
#' @param spatial_network_name name of spatial network to use
#' @param show_grid show spatial grid
#' @param grid_color color of spatial grid
#' @param spatial_grid_name name of spatial grid to use
#' @param show_legend show legend
#' @param point_size_select size of selected points
#' @param point_size_other size of other points
#' @param point_alpha_other opacity of other points
#' @param axis_scale scale of axis
#' @param x_ticks ticks on x-axis
#' @param y_ticks ticks on y-axis
#' @param z_ticks ticks on z-axis
#' @param custom_ratio custom ratio of axes
#' @param show_plot show plots
#' @param return_plot return plotly object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters from \code{\link{all_plots_save_function}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @param \dots additional parameters
#' @return plotly
#' @details Description of parameters.
#' @export
cellProximitySpatPlot3D = function(gobject,
                                   interaction_name = NULL,
                                   cluster_column = NULL,
                                   sdimx = "sdimx",
                                   sdimy = "sdimy",
                                   sdimz = "sdimz",
                                   cell_color = NULL,
                                   cell_color_code = NULL,
                                   color_as_factor = T,
                                   show_other_cells = T,
                                   show_network = T,
                                   show_other_network = F,
                                   network_color = NULL,
                                   spatial_network_name = 'Delaunay_network',
                                   show_grid = F,
                                   grid_color = NULL,
                                   spatial_grid_name = 'spatial_grid',
                                   show_legend = T,
                                   point_size_select = 4,
                                   point_size_other = 2,
                                   point_alpha_other = 0.5,
                                   axis_scale = c("cube","real","custom"),
                                   custom_ratio = NULL,
                                   x_ticks = NULL,
                                   y_ticks = NULL,
                                   z_ticks = NULL,
                                   show_plot = NA,
                                   return_plot = NA,
                                   save_plot = NA,
                                   save_param =  list(),
                                   default_save_name = 'cellProximitySpatPlot3D',
                                   ...){
  if(is.null(sdimz)){
    pl = cellProximityVisPlot_2D_plotly(gobject = gobject,
                                        interaction_name = interaction_name,
                                        cluster_column = cluster_column,
                                        sdimx = sdimx,
                                        sdimy = sdimy,
                                        cell_color = cell_color,
                                        cell_color_code = cell_color_code,
                                        color_as_factor = color_as_factor,
                                        show_other_cells = show_other_cells,
                                        show_network = show_network,
                                        show_other_network = show_other_network,
                                        network_color = network_color,
                                        spatial_network_name = spatial_network_name,
                                        show_grid = show_grid,
                                        grid_color = grid_color,
                                        spatial_grid_name = spatial_grid_name,
                                        show_legend = show_legend,
                                        point_size_select = point_size_select,
                                        point_size_other = point_size_other,
                                        point_alpha_other = point_alpha_other,
                                        axis_scale = axis_scale,
                                        custom_ratio = custom_ratio,
                                        x_ticks = x_ticks,
                                        y_ticks = y_ticks,
                                        ...)
  }
  else{
    pl = cellProximityVisPlot_3D_plotly(gobject = gobject,
                                        interaction_name = interaction_name,
                                        cluster_column = cluster_column,
                                        sdimx = sdimx,
                                        sdimy = sdimy,
                                        sdimz = sdimz,
                                        cell_color = cell_color,
                                        cell_color_code = cell_color_code,
                                        color_as_factor = color_as_factor,
                                        show_other_cells = show_other_cells,
                                        show_network = show_network,
                                        show_other_network = show_other_network,
                                        network_color = network_color,
                                        spatial_network_name = spatial_network_name,
                                        show_grid = show_grid,
                                        grid_color = grid_color,
                                        spatial_grid_name = spatial_grid_name,
                                        show_legend = show_legend,
                                        point_size_select = point_size_select,
                                        point_size_other = point_size_other,
                                        point_alpha_other = point_alpha_other,
                                        axis_scale = axis_scale,
                                        custom_ratio = custom_ratio,
                                        x_ticks = x_ticks,
                                        y_ticks = y_ticks,
                                        z_ticks = z_ticks,
                                        ...)
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

