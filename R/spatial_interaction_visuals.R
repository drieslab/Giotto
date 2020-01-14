

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
#' @examples
#'     cellProximityBarplot(CPscore)
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
  table_mean_results_dc_filter = table_mean_results_dc[original >= min_orig_ints & simulations >= min_sim_ints]
  table_mean_results_dc_filter = table_mean_results_dc_filter[p_higher_orig <= p_val | p_lower_orig <= p_val]

  pl <- ggplot2::ggplot()
  pl <- pl + ggplot2::geom_bar(data = table_mean_results_dc_filter, aes(x = unified_int, y = enrichm, fill = type_int), stat = 'identity', show.legend = F)
  pl <- pl + ggplot2::coord_flip()
  pl <- pl + ggplot2::theme_bw()
  pl <- pl + ggplot2::labs(y = 'enrichment/depletion')
  pl

  bpl <- ggplot2::ggplot()
  bpl <- bpl + ggplot2::geom_bar(data = table_mean_results_dc_filter, aes(x = unified_int, y = original, fill = type_int), stat = 'identity', show.legend = T)
  bpl <- bpl + ggplot2::coord_flip()
  bpl <- bpl + ggplot2::theme_bw() + theme(axis.text.y = element_blank())
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
#' @examples
#'     cellProximityHeatmap(CPscore)
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
  enrich_res[, first_type := strsplit(x = as.character(unified_int), split = '--')[[1]][1], by = 1:nrow(enrich_res)]
  enrich_res[, second_type := strsplit(x = as.character(unified_int), split = '--')[[1]][2], by = 1:nrow(enrich_res)]

  # create matrix
  enrich_mat = data.table::dcast.data.table(data = enrich_res,formula = first_type~second_type, value.var = 'enrichm')
  matrix_d <- as.matrix(enrich_mat[,-1]); rownames(matrix_d) = as.vector(enrich_mat[[1]])
  t_matrix_d <- t(matrix_d)

  # fill in NAs based on values in upper and lower matrix triangle
  t_matrix_d[upper.tri(t_matrix_d)][is.na(t_matrix_d[upper.tri(t_matrix_d)])] = matrix_d[upper.tri(matrix_d)][is.na(t_matrix_d[upper.tri(t_matrix_d)])]
  t_matrix_d[lower.tri(t_matrix_d)][is.na(t_matrix_d[lower.tri(t_matrix_d)])] = matrix_d[lower.tri(matrix_d)][is.na(t_matrix_d[lower.tri(t_matrix_d)])]
  t_matrix_d[is.na(t_matrix_d)] = 0
  final_matrix = t_matrix_d

  # scale data
  if(scale == TRUE) {
    final_matrix <- t(scale(t(final_matrix)))
    final_matrix <- t(final_matrix)
    final_matrix[lower.tri(final_matrix)] <- t(final_matrix)[lower.tri(final_matrix)]
  }

  # # if NA values, impute as mean
  #if(any(is.na(final_matrix)) == TRUE) {
  #  myrowmeans = apply(X = final_matrix, MARGIN = 1, FUN = function(x) mean(na.omit(x)))
  #  mymatrixmeans = matrix(data = rep(myrowmeans, ncol(final_matrix)), nrow = nrow(final_matrix), ncol = ncol(final_matrix))
  #  final_matrix[is.na(final_matrix)] = mymatrixmeans[which(is.na(final_matrix))]
  #}

  # order cell types
  if(order_cell_types == TRUE) {

    cordist = stats::as.dist(1-cor(final_matrix))
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
#' @examples
#'     cellProximityNetwork(CPscore)
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


  if(layout == 'Fruchterman') {
    coords = igraph::layout_with_fr(graph = igd, weights = edges_sizes_resc)
  } else if(layout == 'DrL') {
    coords = igraph::layout_with_drl(graph = igd, weights = edges_sizes_resc)
  } else if(layout == 'Kamada-Kawai') {
    coords = igraph::layout_with_kk(graph = igd, weights = edges_sizes_resc)
  } else {
    stop('\n Currently no other layouts have been implemented \n')
  }

  #iplot = igraph::plot.igraph(igd, edge.color = edges_colors, edge.width = edges_sizes_resc, layout = coords)

  igd = igraph::set.edge.attribute(graph = igd, index = igraph::E(igd), name = 'color', value = edges_colors)
  igd = igraph::set.edge.attribute(graph = igd, index = igraph::E(igd), name = 'size', value = as.numeric(edges_sizes_resc))

  ## only show attractive edges
  if(only_show_enrichment_edges == TRUE) {
    colors = igraph::get.edge.attribute(igd, name = 'color')
    subvertices_ids = which(colors == 'enriched')
    igd = igraph::subgraph.edges(graph = igd, eids = subvertices_ids)
  }

  #longDT = as.data.table(igraph::as_long_data_frame(igd))
  #return(longDT)
  #return(list(igd, coords))

  ## create plot
  gpl = ggraph::ggraph(graph = igd, layout = coords)
  gpl = gpl + ggraph::geom_edge_link(aes(color = factor(color), edge_width = size, edge_alpha = size), show.legend = F)
  gpl = gpl + ggraph::geom_edge_loop(aes(color = factor(color), edge_width = size, edge_alpha = size, strength = self_loop_strength), show.legend = F)
  gpl = gpl + ggraph::scale_edge_color_manual(values = c('enriched' = color_enrichment, 'depleted' = color_depletion))
  gpl = gpl + ggraph::scale_edge_width(range = edge_width_range)
  gpl = gpl + ggraph::scale_edge_alpha(range = c(0.1,1))
  gpl = gpl + ggraph::geom_node_text(aes(label = name), repel = TRUE, size = node_text_size)
  gpl = gpl + ggraph::geom_node_point(size = node_size)
  gpl = gpl + ggplot2::theme_bw() + ggplot2::theme(panel.grid = element_blank(),
                                                   panel.border = element_blank(),
                                                   axis.title = element_blank(),
                                                   axis.text = element_blank(),
                                                   axis.ticks = element_blank())


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
#' @param gobject giotto object
#' @param interaction_name cell-cell interaction name
#' @param cluster_column cluster column with cell clusters
#' @param sdimx x-axis dimension name (default = 'sdimx')
#' @param sdimy y-axis dimension name (default = 'sdimy')
#' @param cell_color color for cells (see details)
#' @param cell_color_code named vector with colors
#' @param color_as_factor convert color column to factor
#' @param show_other_cells decide if show cells not in network
#' @param show_network show underlying spatial network
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
#' @param point_other_border_col border color of other points
#' @param point_other_border_stroke stroke size of other points
#' @return ggplot
#' @details Description of parameters.
#' @export
#' @examples
#'     cellProximityVisPlot_2D_ggplot(gobject)
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
                                           spatial_network_name = 'spatial_network',
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
        cell_color_code = Giotto:::getDistinctColors(n = number_colors)
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
#' @param gobject giotto object
#' @param interaction_name cell-cell interaction name
#' @param cluster_column cluster column with cell clusters
#' @param sdimx x-axis dimension name (default = 'sdimx')
#' @param sdimy y-axis dimension name (default = 'sdimy')
#' @param cell_color color for cells (see details)
#' @param cell_color_code named vector with colors
#' @param color_as_factor convert color column to factor
#' @param show_other_cells decide if show cells not in network
#' @param show_network show underlying spatial network
#' @param network_color color of spatial network
#' @param spatial_network_name name of spatial network to use
#' @param show_grid show spatial grid
#' @param grid_color color of spatial grid
#' @param spatial_grid_name name of spatial grid to use
#' @param coord_fix_ratio fix ratio between x and y-axis
#' @param show_legend show legend
#' @param point_size_select size of selected points
#' @return plotly
#' @details Description of parameters.
#' @export
#' @examples
#'     cellProximityVisPlot_2D_plotly(gobject)

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
                                           spatial_network_name = 'spatial_network',
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
        cell_color_code = Giotto:::getDistinctColors(n = number_colors)
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
      pl <- pl %>% add_trace(type = 'scatter',mode = 'markers',
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
#' @param show_network show underlying spatial network
#' @param network_color color of spatial network
#' @param spatial_network_name name of spatial network to use
#' @param show_grid show spatial grid
#' @param grid_color color of spatial grid
#' @param spatial_grid_name name of spatial grid to use
#' @param coord_fix_ratio fix ratio between x and y-axis
#' @param show_legend show legend
#' @param point_size_select size of selected points
#' @return plotly
#' @details Description of parameters.
#' @export
#' @examples
#'     cellProximityVisPlot_3D_plotly(gobject)

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
                                           spatial_network_name = 'spatial_network',
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
        cell_color_code = Giotto:::getDistinctColors(n = number_colors)
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
      pl <- pl %>% add_trace(type = 'scatter3d',mode = 'markers',
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
#' @param show_network show underlying spatial network
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
#' @param point_other_border_col border color of other points
#' @param point_other_border_stroke stroke size of other points
#' @return ggplot or plotly
#' @details Description of parameters.
#' @export
#' @examples
#'     cellProximityVisPlot(gobject)
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
                                 spatial_network_name = 'spatial_network',
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





#' @title showGeneExpressionProximityScore
#' @name showGeneExpressionProximityScore
#' @description Create heatmap from cell-cell proximity scores
#' @param scores CPscore, output from getAverageCellProximityGeneScores()
#' @param selected_gene gene to show
#' @param sort_column column name to use for sorting
#' @return ggplot barplot
#' @details Give more details ...
#' @export
#' @examples
#'     showGeneExpressionProximityScore(scores)
showGeneExpressionProximityScore <- function(scores,
                                             selected_gene,
                                             sort_column = 'diff_spat') {

  # gene specific
  subset_scores = scores[genes == selected_gene]

  order_dt = subset_scores[order(get(sort_column))]
  order_ints = order_dt$interaction

  specific_int_scores = subset_scores[,.(genes, interaction, cell_expr_1, cell_expr_2, comb_expr,
                                         all_cell_expr_1, all_cell_expr_2, all_comb_expr,
                                         diff_spat, diff_spat_1, diff_spat_2)]
  specific_int_scores = data.table::melt.data.table(data = specific_int_scores, id.vars = c('genes', 'interaction'))
  specific_int_scores[, value := as.numeric(value)]
  specific_int_scores[, interaction := factor(interaction, order_ints)]
  specific_int_scores[, score := ifelse(grepl('spat', variable), 'difference',
                                        ifelse(grepl('all', variable), 'all cells', 'int cells'))]
  specific_int_scores[, variable := factor(variable, levels = c('cell_expr_1', 'cell_expr_2', 'comb_expr',
                                                                'all_cell_expr_1', 'all_cell_expr_2', 'all_comb_expr',
                                                                'diff_spat_1', 'diff_spat_2', 'diff_spat'))]

  pl <- ggplot2::ggplot()
  pl <- pl + ggplot2::theme_classic()
  pl <- pl + ggplot2::geom_bar(data = specific_int_scores, aes(x = interaction, y = as.numeric(value), fill = score), stat = 'identity')
  pl <- pl + ggplot2::facet_wrap(~ variable, ncol = 9)
  pl <- pl + ggplot2::coord_flip()
  pl <- pl + ggplot2::labs(x = '', y = 'normalized expression', title = selected_gene)
  pl

  print(pl)


}




#' @title showIntExpressionProximityScore
#' @name showIntExpressionProximityScore
#' @description Create heatmap from cell-cell proximity scores
#' @param scores scores, output from getAverageCellProximityGeneScores()
#' @param selected_interaction interaction to show
#' @param sort_column column name to use for sorting
#' @param show_enriched_n show top enriched interactions
#' @param show_depleted_n show top depleted interactions
#' @return ggplot barplot
#' @details Give more details ...
#' @export
#' @examples
#'     showIntExpressionProximityScore(scores)
showIntExpressionProximityScore <- function(scores,
                                            selected_interaction,
                                            sort_column = 'diff_spat',
                                            show_enriched_n = 5,
                                            show_depleted_n = 5) {

  # interaction specific
  subset_scores = scores[interaction == selected_interaction]

  order_dt = subset_scores[order(get(sort_column))]
  order_genes = order_dt$genes
  head_genes = head(order_genes, show_enriched_n)
  tail_genes = tail(order_genes, show_depleted_n)

  specific_gene_scores = subset_scores[,.(genes, interaction, cell_expr_1, cell_expr_2, comb_expr,
                                          all_cell_expr_1, all_cell_expr_2, all_comb_expr,
                                          diff_spat, diff_spat_1, diff_spat_2)]
  specific_gene_scores = data.table::melt.data.table(data = specific_gene_scores, id.vars = c('genes', 'interaction'))
  specific_gene_scores[, value := as.numeric(value)]
  specific_gene_scores[, genes := factor(genes, order_genes)]
  specific_gene_scores[, score := ifelse(grepl('spat', variable), 'difference',
                                         ifelse(grepl('all', variable), 'all cells', 'int cells'))]
  specific_gene_scores[, variable := factor(variable, levels = c('cell_expr_1', 'cell_expr_2', 'comb_expr',
                                                                 'all_cell_expr_1', 'all_cell_expr_2', 'all_comb_expr',
                                                                 'diff_spat_1', 'diff_spat_2', 'diff_spat'))]

  specific_gene_scores = specific_gene_scores[genes %in% c(head_genes, tail_genes)]

  pl <- ggplot2::ggplot()
  pl <- pl + ggplot2::theme_classic()
  pl <- pl + ggplot2::geom_bar(data = specific_gene_scores, aes(x = genes, y = as.numeric(value), fill = score), stat = 'identity')
  pl <- pl + ggplot2::facet_wrap(~ variable, ncol = 9)
  pl <- pl + ggplot2::coord_flip()
  pl <- pl + ggplot2::labs(x = '', y = 'normalized expression', title = '1-3')
  pl

  print(pl)

}


#' @title showTopGeneToGene
#' @name showTopGeneToGene
#' @description Show enriched/depleted gene-gene enrichments
#' @param GTGscore GTGscore, output from getGeneToGeneScores()
#' @param top_interactions number of top gene-gene enrichments to show
#' @param direction show top increased or decreased gene-gene enrichments
#' @param complement_data include non-enriched gene-gene scores from other cell-cell interactions
#' @param subset_cell_ints subset cell-cell interactions to show
#' @param subset_genes subset genes to show
#' @return ggplot barplot
#' @details Give more details ...
#' @export
#' @examples
#'     showTopGeneToGene(scores)
showTopGeneToGene = function(GTGscore,
                             top_interactions = 10,
                             direction = c('increased', 'decreased'),
                             complement_data = T,
                             subset_cell_ints = NULL,
                             subset_genes = NULL) {

  direction = match.arg(direction, choices = c('increased', 'decreased'))
  data.table::setorder(GTGscore, interaction, -diff_spat)


  # subset genes or cell cell interactions
  if(!is.null(subset_cell_ints)) {
    GTGscore = GTGscore[interaction %in% subset_cell_ints]
  }

  if(!is.null(subset_genes)) {
    GTGscore = GTGscore[genes_cell_1 %in% subset_genes | genes_cell_2 %in% subset_genes]
  }


  ## increased ##
  if(direction == 'increased') {

    topDT = GTGscore[, head(.SD, top_interactions), by = interaction]
    if(complement_data == TRUE) {
      topDT = GTGscore[gene_gene %in% topDT$gene_gene]
    }
    topDT[, direction := 'increased']

  }

  ## decreased ##
  else if(direction == 'decreased') {

    topDT = GTGscore[, tail(.SD, top_interactions), by = interaction]
    if(complement_data == TRUE) {
      topDT = GTGscore[gene_gene %in% topDT$gene_gene]
    }
    topDT[, direction := 'decreased']

  }

  # sort gene_gene
  setorder(topDT, interaction, -diff_spat)
  gene_gene_order = unique(topDT[['gene_gene']])
  topDT[, gene_gene := factor(gene_gene, level = gene_gene_order)]

  pl <- ggplot2::ggplot()
  pl <- pl + ggplot2::theme_classic()
  pl <- pl + ggplot2::geom_segment(data = topDT, aes(x = comb_expr, xend = all_comb_expr, y = gene_gene, yend = gene_gene), linetype = 2)
  pl <- pl + ggplot2::geom_point(data = topDT, aes_string(x = 'all_comb_expr', y = 'gene_gene'), color = 'blue', size = 2)
  pl <- pl + ggplot2::geom_point(data = topDT, aes_string(x = 'comb_expr', y = 'gene_gene'), color = 'red', size = 2)
  pl <- pl + ggplot2::facet_wrap( ~ interaction, ncol = length(unique(topDT$interaction)))
  pl

}



#' @title direction_test_CPG
#' @name direction_test_CPG
#' @description shows direction of change
#' @examples
#'     direction_test_CPG()
direction_test = function(x, min_fdr = 0.05) {

  fdr_1 = as.numeric(x[['fdr_1']])
  if(is.na(fdr_1) == TRUE) fdr_1 = 1
  fdr_2 = as.numeric(x[['fdr_2']])
  if(is.na(fdr_2) == TRUE) fdr_2 = 1

  log2fc_1 = as.numeric(x[['log2fc_spat_1']])
  log2fc_2 = as.numeric(x[['log2fc_spat_2']])

  if(fdr_1 > min_fdr & fdr_2 > min_fdr) {
    return('n.s.')
  } else if(fdr_1 <= min_fdr & fdr_2 <= min_fdr) {

    if(log2fc_1 > 0 & log2fc_2 > 0) {
      return('both_up')
    } else if(log2fc_1 < 0 & log2fc_2 < 0) {
      return('both_down')
    } else {
      return('mixed')
    }

  } else if((fdr_1 <= min_fdr & log2fc_1 > 0) | (fdr_2 <= min_fdr & log2fc_2 > 0)) {
    return('up')
  } else if((fdr_1 <= min_fdr & log2fc_1 < 0) | (fdr_2 <= min_fdr & log2fc_2 < 0)) {
    return('down')
  } else {
    return('missing')
  }

}



#' @title filterCPGscores
#' @name filterCPGscores
#' @description visualize Cell Proximity Gene enrichment scores
#' @param method visualization method
#' @param min_cells min number of cells threshold
#' @param min_fdr false_discovery threshold
#' @param min_spat_diff spatial difference threshold
#' @param min_log2_fc min log2 fold-change
#' @param keep_int_duplicates keep both cell_A-cell_B and cell_B-cell_A
#' @param direction expression changes to keep
#' @return Gene to gene scores in data.table format
#' @details This function filters the output from getCellProximityGeneScores based on
#' false-discovery rate, minimum absolute difference, minimum log fold-change and
#' direction of change.
#' @export
#' @examples
#'     filterCPGscores(CPGscore)
filterCPGscores = function(CPGscore,
                           min_cells = 5,
                           min_fdr = 0.05,
                           min_spat_diff = 0.2,
                           min_log2_fc = 0.5,
                           keep_int_duplicates = TRUE,
                           direction = c('both', 'up', 'down')) {



  # other parameters
  direction = match.arg(direction, choices = c('both', 'up', 'down'))

  if(keep_int_duplicates == FALSE) {
    selection_scores = data.table::copy(CPGscore[unif_int_rank == 1])
  } else {
    selection_scores = data.table::copy(CPGscore)
  }

  selection_scores = CPGscore[(nr_1 >= min_cells & fdr_1 <= min_fdr &
                                 abs(diff_spat_1) >= min_spat_diff &
                                 abs(log2fc_spat_1) >= min_log2_fc) |
                                (nr_2 >= min_cells & fdr_2 <= min_fdr &
                                   abs(diff_spat_2) >= min_spat_diff &
                                   abs(log2fc_spat_2) >= min_log2_fc)]

  if(direction == 'both') {
    selection_scores = selection_scores
  } else if(direction == 'up') {
    selection_scores = selection_scores[(log2fc_spat_1 >= min_log2_fc & fdr_1 <= min_fdr & diff_spat_1 >= min_spat_diff) |
                                          (log2fc_spat_2 >= min_log2_fc & fdr_2 <= min_fdr & diff_spat_2 >= min_spat_diff)]
  } else if(direction == 'down') {
    selection_scores = selection_scores[(log2fc_spat_1 <= -min_log2_fc & fdr_1 <= min_fdr & diff_spat_1 <= -min_spat_diff) |
                                          (log2fc_spat_2 <= -min_log2_fc & fdr_2 <= min_fdr & diff_spat_2 <= -min_spat_diff)]
  }


  change_values = unlist(apply(selection_scores, MARGIN = 1, FUN = function(x) {
    direction_test(x, min_fdr = min_fdr)
  }))
  selection_scores[, change := change_values]
  return(selection_scores)

}



#' @title showCPGscores
#' @name showCPGscores
#' @description visualize Cell Proximity Gene enrichment scores
#' @param CPGscore CPGscore, output from getCellProximityGeneScores()
#' @param method visualization method
#' @param min_cells min number of cells threshold
#' @param min_fdr fdr threshold
#' @param min_spat_diff spatial difference threshold
#' @param min_log2_fc min log2 fold-change
#' @param keep_int_duplicates keep both cell_A-cell_B and cell_B-cell_A
#' @param direction up or downregulation or both
#' @param cell_color_code color code for cell types
#' @param show_plot show plot
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters from \code{\link{all_plots_save_function}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return Gene to gene scores in data.table format
#' @details Give more details ...
#' @export
#' @examples
#'     showCPGscores(CPGscore)
showCPGscores = function(gobject,
                         CPGscore,
                         method = c('volcano', 'cell_barplot', 'cell-cell', 'cell_sankey'),
                         min_cells = 5,
                         min_fdr = 0.05,
                         min_spat_diff = 0.2,
                         min_log2_fc = 0.5,
                         keep_int_duplicates = TRUE,
                         direction = c('both', 'up', 'down'),
                         cell_color_code = NULL,
                         show_plot = NA,
                         return_plot = NA,
                         save_plot = NA,
                         save_param =  list(),
                         default_save_name = 'showCPGscores'
) {


  # print, return and save parameters
  show_plot = ifelse(is.na(show_plot), readGiottoInstructions(gobject, param = 'show_plot'), show_plot)
  save_plot = ifelse(is.na(save_plot), readGiottoInstructions(gobject, param = 'save_plot'), save_plot)
  return_plot = ifelse(is.na(return_plot), readGiottoInstructions(gobject, param = 'return_plot'), return_plot)


  ## first filter
  filter_set = filterCPGscores(CPGscore = CPGscore,
                               min_cells = min_cells,
                               min_fdr = min_fdr,
                               min_spat_diff = min_spat_diff,
                               min_log2_fc = min_log2_fc,
                               keep_int_duplicates = keep_int_duplicates,
                               direction = direction)


  ## other parameters
  method = match.arg(method, choices = c('volcano', 'cell_barplot', 'cell-cell', 'cell_sankey'))

  ## create data.table for visualization
  subset = filter_set[(fdr_1 <= min_fdr & unif_int_rank == 1) | (fdr_1 <= min_fdr & fdr_2 <= min_fdr & unif_int_rank == 2)]
  part1 = subset[unif_int_rank == 1,.(genes, interaction, cell_type_1, cell_type_2, log2fc_spat_1, fdr_1)]
  colnames(part1) = c('genes', 'interaction', 'source', 'neighbor', 'log2fc', 'fdr')
  part2 = subset[unif_int_rank == 2,.(genes, interaction, cell_type_2, cell_type_1, log2fc_spat_2, fdr_2)]
  colnames(part2) = c('genes', 'interaction', 'source', 'neighbor', 'log2fc', 'fdr')
  complete_part = rbind(part1, part2)



  if(method == 'volcano') {

    ## volcanoplot
    pl = ggplot()
    pl = pl + geom_point(data = complete_part, aes(x = log2fc, y = -log10(fdr)))
    pl = pl + theme_classic()
    pl = pl + geom_vline(xintercept = 0, linetype = 2)


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

    nr_int_selection_scores = complete_part[, .N, by = interaction]
    order_interactions = nr_int_selection_scores[order(N)]$interaction

    complete_part[, interaction := factor(interaction, order_interactions)]

    pl <- ggplot2::ggplot()
    pl <- pl + ggplot2::geom_bar(data = complete_part, aes(x = interaction, fill = interaction))
    pl <- pl + ggplot2::theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1))
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
    nr_source_selection_scores = complete_part[, .N, by = source]
    order_source = nr_source_selection_scores[order(N)]$source

    complete_part[, source := factor(source, order_source)]

    pl <- ggplot2::ggplot()
    pl <- pl + ggplot2::geom_bar(data = complete_part, aes(x = source, fill = neighbor))
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

    testalluv = complete_part[, .N, by = c('neighbor', 'source')]

    library(ggalluvial)

    pl <- ggplot2::ggplot(testalluv,
                          aes(y = N, axis1 = source, axis2 = neighbor)) +
      ggalluvial::geom_alluvium(aes(fill = source), width = 1/12) +
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

  }

}



#' @title showGTGscores
#' @name showGTGscores
#' @description visualize Cell Proximity Gene enrichment scores
#' @param CPGscore CPGscore, output from getCellProximityGeneScores()
#' @param method visualization method
#' @param min_cells min number of cells threshold
#' @param min_pval p-value threshold
#' @param min_spat_diff spatial difference threshold
#' @param min_log2_fc log2 fold-change threshold
#' @param direction up or downregulation or both
#' @param cell_color_code color code for cell types
#' @param show_plot print plot
#' @param specific_genes_1 subset of genes, matched with specific_genes_2
#' @param specific_genes_2 subset of genes, matched with specific_genes_1
#' @param first_cell_name name for first cells
#' @param second_cell_name name for second cells
#' @return ggplot
#' @details Give more details ...
#' @export
#' @examples
#'     showGTGscores(CPGscore)
showGTGscores = function(GTGscore,
                         method = c('cell_barplot', 'cell-cell', 'cell_sankey'),
                         min_cells = 5,
                         min_pval = 0.05,
                         min_spat_diff = 0.2,
                         min_log2_fc = 0.5,
                         direction = c('both', 'up', 'down'),
                         cell_color_code = NULL,
                         show_plot = T,
                         specific_genes_1 = NULL,
                         specific_genes_2 = NULL,
                         first_cell_name = 'ligand cell',
                         second_cell_name = 'receptor cell',
                         return_DT = F) {


  direction = match.arg(direction, choices = c('both', 'up', 'down'))
  method = match.arg(method, choices = c('cell_barplot', 'cell-cell', 'cell_sankey'))


  ## filter ##

  ## p-value, nr of cells and spatial difference
  if(direction == 'both') {
    GTGscore = GTGscore[(nr_1 >= min_cells & pval_1 < min_pval & abs(diff_spat_1) >= min_spat_diff & abs(log2fc_spat_1) >= min_log2_fc) &
                          (nr_2 >= min_cells & pval_2 < min_pval & abs(diff_spat_2) >= min_spat_diff & abs(log2fc_spat_2) >= min_log2_fc)]
  } else if(direction == 'up') {
    GTGscore = GTGscore[(nr_1 >= min_cells & pval_1 < min_pval & diff_spat_1 >= min_spat_diff & log2fc_spat_1 >= min_log2_fc) &
                          (nr_2 >= min_cells & pval_2 < min_pval & diff_spat_2 >= min_spat_diff & log2fc_spat_2 >= min_log2_fc)]
  } else if(direction == 'down') {
    GTGscore = GTGscore[(nr_1 >= min_cells & pval_1 < min_pval & diff_spat_1 <= -min_spat_diff & log2fc_spat_1 <= -min_log2_fc) &
                          (nr_2 >= min_cells & pval_2 < min_pval & diff_spat_2 <= -min_spat_diff & log2fc_spat_2 <= -min_log2_fc)]
  }



  ## specific genes
  if((!is.null(specific_genes_1) & !is.null(specific_genes_2))) {

    if(length(specific_genes_1) != length(specific_genes_2)) {
      stop('\n specific_genes_1 must have the same length as specific_genes_2')
    }

    LR_combo = paste0(specific_genes_1,'-', specific_genes_2)
    RL_combo = paste0(specific_genes_2,'-', specific_genes_1)
    GTGscore = GTGscore[gene_gene %in% c(LR_combo, RL_combo)]

  }


  # return data.table
  if(return_DT == TRUE) {
    return(GTGscore)
  }



  if(method == 'cell-cell') {

    nr_int_GTGscore = GTGscore[, .N, by = unified_int]
    order_interactions = nr_int_GTGscore[order(N)]$unified_int

    GTGscore[, unified_int := factor(unified_int, order_interactions)]
    GTGscore[, type_int := ifelse(cell_type_1 == cell_type_2, 'homo', 'hetero')]

    pl <- ggplot2::ggplot()
    pl <- pl + ggplot2::geom_bar(data = GTGscore, aes(x = unified_int, fill = type_int))
    pl <- pl + ggplot2::theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1))
    pl <- pl + ggplot2::coord_flip()

    if(show_plot == TRUE) {
      print(pl)
    }

  } else if(method == 'cell_barplot') {

    if((!is.null(specific_genes_1) & !is.null(specific_genes_2))) {

      LR_GTGscore = GTGscore[(genes_1 %in% specific_genes_1 & genes_2 %in% specific_genes_2)]
      RL_GTGscore = GTGscore[(genes_1 %in% specific_genes_2 & genes_2 %in% specific_genes_1)]

      hipLR = LR_GTGscore[, .N, by = .(cell_type_1, cell_type_2)]
      setnames(hipLR, c('cell_type_1', 'cell_type_2'), c('first_cell', 'second_cell'))
      hipRL = RL_GTGscore[, .N, by = .(cell_type_2, cell_type_1)]
      setnames(hipRL, c('cell_type_2', 'cell_type_1'), c('first_cell', 'second_cell'))
      hip = rbind(hipLR, hipRL)

      hip[, type_int := ifelse(first_cell == second_cell, 'homo', 'hetero')]

    } else {

      hip = GTGscore[, .N, by = .(cell_type_1, cell_type_2)]
      setnames(hip, c('cell_type_1', 'cell_type_2'), c('first_cell', 'second_cell'))

    }

    first_order = hip[, sum(N), by = first_cell]
    first_cell_order = first_order[order(V1)][['first_cell']]

    second_order = hip[, sum(N), by = second_cell]
    second_cell_order = second_order[order(V1)][['second_cell']]

    hip[, first_cell := factor(first_cell, levels = first_cell_order)]
    hip[, second_cell := factor(second_cell, levels = second_cell_order)]


    pl <- ggplot2::ggplot()
    pl <- pl + ggplot2::geom_bar(data = hip, aes(x = first_cell, fill = second_cell, y = N), stat = 'identity')
    pl <- pl + ggplot2::theme_classic()  + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

    if(!is.null(cell_color_code)) {
      pl <- pl + ggplot2::scale_fill_manual(values = cell_color_code)
    }

    if(show_plot == TRUE) {
      print(pl)
    }


  } else if(method == 'cell_sankey') {


    if((!is.null(specific_genes_1) & !is.null(specific_genes_2))) {

      LR_GTGscore = GTGscore[(genes_1 %in% specific_genes_1 & genes_2 %in% specific_genes_2)]
      RL_GTGscore = GTGscore[(genes_1 %in% specific_genes_2 & genes_2 %in% specific_genes_1)]

      hipLR = LR_GTGscore[, .N, by = .(cell_type_1, cell_type_2)]
      setnames(hipLR, c('cell_type_1', 'cell_type_2'), c('first_cell', 'second_cell'))
      hipRL = RL_GTGscore[, .N, by = .(cell_type_2, cell_type_1)]
      setnames(hipRL, c('cell_type_2', 'cell_type_1'), c('first_cell', 'second_cell'))
      hip = rbind(hipLR, hipRL)

      hip[, type_int := ifelse(first_cell == second_cell, 'homo', 'hetero')]

    } else {

      hip = GTGscore[, .N, by = .(cell_type_1, cell_type_2)]
      setnames(hip, c('cell_type_1', 'cell_type_2'), c('first_cell', 'second_cell'))

    }

    first_order = hip[, sum(N), by = first_cell]
    first_cell_order = first_order[order(V1)][['first_cell']]

    second_order = hip[, sum(N), by = second_cell]
    second_cell_order = second_order[order(V1)][['second_cell']]

    hip[, first_cell := factor(first_cell, levels = first_cell_order)]
    hip[, second_cell := factor(second_cell, levels = second_cell_order)]


    pl <- ggplot2::ggplot(hip,
                 aes(y = N, axis1 = first_cell, axis2 = second_cell)) +
      ggalluvial::geom_alluvium(aes(fill = first_cell), width = 1/12) +
      ggalluvial::geom_stratum(width = 1/12, fill = "black", color = "grey") +
      ggplot2::scale_x_discrete(limits = c(first_cell_name, second_cell_name), expand = c(.05, .05)) +
      ggplot2::geom_label(stat = "stratum", label.strata = TRUE, size = 3) +
      ggplot2::theme_classic() + ggplot2::labs(x = '', y = '# spatial L-R')
    if(!is.null(cell_color_code)) {
      pl <- pl + ggplot2::scale_fill_manual(values = cell_color_code)
    }
    if(show_plot == TRUE) {
      print(pl)
    }

  }

  return(pl)

}




#' @title plotGTGscores
#' @name plotGTGscores
#' @description Create heatmap from cell-cell proximity scores
#' @param gobject giotto object
#' @param GTGscore GTGscore, output from getGeneToGeneScores()
#' @param selected_interactions interactions to show
#' @param selected_genes genes to show
#' @param detail_plot show detailed info in both interacting cell types
#' @param simple_plot show a simplified plot
#' @param simple_plot_facet facet on interactions or genes with simple plot
#' @param facet_scales ggplot facet scales paramter
#' @param facet_ncol ggplot facet ncol parameter
#' @param facet_nrow ggplot facet nrow parameter
#' @param colors vector with 2 colors to represent respectively all and selected cells
#' @param show_plot show plots
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters from \code{\link{all_plots_save_function}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return ggplot barplot
#' @details Give more details ...
#' @export
#' @examples
#'     plotGTGscores(GTGscore)
plotGTGscores <- function(gobject,
                          GTGscore,
                          selected_interactions = NULL,
                          selected_gene_to_gene = NULL,
                          detail_plot = T,
                          simple_plot = F,
                          simple_plot_facet = c('interaction', 'genes'),
                          facet_scales = 'fixed',
                          facet_ncol = length(selected_gene_to_gene),
                          facet_nrow = length(selected_interactions),
                          colors = c('blue', 'red'),
                          show_plot = NA,
                          return_plot = NA,
                          save_plot = NA,
                          save_param =  list(),
                          default_save_name = 'plotGTGscores') {

  if(is.null(selected_interactions) | is.null(selected_gene_to_gene)) {
    stop('\n You need to provide a selection of cell-cell interactions and genes-genes to plot \n')
  }


  subDT = GTGscore[unif_gene_gene %in% selected_gene_to_gene & unified_int %in% selected_interactions]

  subDT[, all_cell_expr_1 := as.numeric(all_cell_expr_1)]
  subDT[, cell_expr_1 := as.numeric(cell_expr_1)]
  subDT[, all_cell_expr_2 := as.numeric(all_cell_expr_2)]
  subDT[, cell_expr_2 := as.numeric(cell_expr_2)]
  subDT[, all_cell_expr := as.numeric(all_cell_expr)]
  subDT[, spatial_cell_expr := as.numeric(spatial_cell_expr)]


  # order interactions and gene-to-gene according to input
  subDT[, gene_gene := factor(gene_gene, levels = selected_gene_to_gene)]
  subDT[, unified_int := factor(unified_int, levels = selected_interactions)]


    if(simple_plot == F) {

    pl <- ggplot2::ggplot()
    pl <- pl + ggplot2::theme_bw()

    if(detail_plot == TRUE) {
      pl <- pl + ggplot2::geom_point(data = subDT, aes(x = 0, y = all_cell_expr_2, colour = "all cell expression"),shape = 1)
      pl <- pl + ggplot2::geom_point(data = subDT, aes(x = 0, y = cell_expr_2,colour = "selected cell expression"), shape = 1)
      pl <- pl + ggplot2::geom_point(data = subDT, aes(x = all_cell_expr_1, y = 0, colour = "all cell expression"), shape = 1)
      pl <- pl + ggplot2::geom_point(data = subDT, aes(x = cell_expr_1, y = 0,colour ="selected cell expression"), shape = 1)
    }

    pl <- pl + ggplot2::geom_point(data = subDT, aes(x = all_cell_expr_1, y = all_cell_expr_2,colour = "all cell expression"),size = 2)
    pl <- pl + ggplot2::geom_point(data = subDT, aes(x = cell_expr_1, y = cell_expr_2,colour ="selected cell expression"), size = 2)
    pl <- pl + ggplot2::geom_segment(data = subDT, aes(x = all_cell_expr_1, xend = cell_expr_1,
                                              y = all_cell_expr_2, yend = cell_expr_2), linetype = 2)
    #pl <- pl + ggplot2::labs(x = 'gene 1 in celltype 1', y = 'gene 2 in celltype 2')
    pl <- pl + ggplot2::labs(x = paste(subDT$genes_1,subDT$cell_type_1,sep = " in ")
                             , y = paste(subDT$genes_2,subDT$cell_type_2,sep = " in "))
    pl <- pl + ggplot2::scale_colour_manual(name="expression source",values = colors)
    pl <- pl + ggplot2::facet_wrap(~unif_gene_gene+unified_int, nrow = facet_nrow, ncol = facet_ncol,
                          scales = facet_scales)

  }else {

    simple_plot_facet = match.arg(arg = simple_plot_facet, choices = c('interaction', 'genes'))

    if(simple_plot_facet == 'interaction') {
      pl <- ggplot2::ggplot()
      pl <- pl + ggplot2::theme_bw()
      pl <- pl + ggplot2::geom_segment(data = subDT, aes(x = all_cell_expr, xend = spatial_cell_expr,
                                                y = unif_gene_gene, yend = unif_gene_gene), linetype = 2)
      pl <- pl + ggplot2::geom_point(data = subDT, aes(x = all_cell_expr, y = unif_gene_gene,colour = "all cell expression"))
      pl <- pl + ggplot2::geom_point(data = subDT, aes(x = spatial_cell_expr, y = unif_gene_gene,colour ="selected cell expression"))
      pl <- pl + ggplot2::scale_colour_manual(name="expression source",values=cols)
      pl <- pl + ggplot2::facet_wrap(~unified_int, scales = facet_scales)
      pl <- pl + ggplot2::labs(x = 'interactions', y = 'gene-gene')
    } else {
      pl <- ggplot2::ggplot()
      pl <- pl + ggplot2::theme_bw()
      pl <- pl + ggplot2::geom_segment(data = subDT, aes(x = all_cell_expr, xend = spatial_cell_expr,
                                                y = unified_int, yend = unified_int), linetype = 2)
      pl <- pl + ggplot2::geom_point(data = subDT, aes(x = all_cell_expr, y = unified_int,colour = "all cell expression"))
      pl <- pl + ggplot2::geom_point(data = subDT, aes(x = spatial_cell_expr, y = unified_int,colour ="selected cell expression"))
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

  return(pl)

}




#' @title plotCPGscores
#' @name plotCPGscores
#' @description Create heatmap from cell-cell proximity scores
#' @param CPGscores CPGscores, output from getCellProximityGeneScores()
#' @param selected_interactions interactions to show
#' @param selected_genes genes to show
#' @param detail_plot show detailed info in both interacting cell types
#' @param simple_plot show a simplified plot
#' @param simple_plot_facet facet on interactions or genes with simple plot
#' @param facet_scales ggplot facet scales paramter
#' @param facet_ncol ggplot facet ncol parameter
#' @param facet_nrow ggplot facet nrow parameter
#' @param show_plot show plot
#' @return ggplot barplot
#' @details Give more details ...
#' @export
#' @examples
#'     plotCPGscores(CPGscores)
plotCPGscores <- function(CPGscores,
                          selected_interactions = NULL,
                          selected_genes = NULL,
                          detail_plot = T,
                          simple_plot = F,
                          simple_plot_facet = c('interaction', 'genes'),
                          facet_scales = 'fixed',
                          facet_ncol = length(selected_genes),
                          facet_nrow = length(selected_interactions),
                          show_plot = F) {

  if(is.null(selected_interactions) | is.null(selected_genes)) {
    stop('\n You need to provide a selection of cell-cell interactions and genes to plot \n')
  }


  subDT = CPGscores[genes %in% selected_genes & interaction %in% selected_interactions]

  subDT[, all_cell_expr_1 := as.numeric(all_cell_expr_1)]
  subDT[, cell_expr_1 := as.numeric(cell_expr_1)]
  subDT[, all_cell_expr_2 := as.numeric(all_cell_expr_2)]
  subDT[, cell_expr_2 := as.numeric(cell_expr_2)]
  subDT[, all_comb_expr := as.numeric(all_comb_expr)]
  subDT[, comb_expr := as.numeric(comb_expr)]


  # order interactions and gene-to-gene according to input
  subDT[, genes := factor(genes, levels = selected_genes)]
  subDT[, interaction := factor(interaction, levels = selected_interactions)]


  if(simple_plot == F) {

    pl <- ggplot2::ggplot()
    pl <- pl + ggplot2::theme_bw()

    if(detail_plot == TRUE) {
      pl <- pl + ggplot2::geom_point(data = subDT, aes(x = 0, y = all_cell_expr_2), color = 'blue', shape = 1)
      pl <- pl + ggplot2::geom_point(data = subDT, aes(x = 0, y = cell_expr_2), color = 'red', shape = 1)
      pl <- pl + ggplot2::geom_point(data = subDT, aes(x = all_cell_expr_1, y = 0), color = 'blue', shape = 1)
      pl <- pl + ggplot2::geom_point(data = subDT, aes(x = cell_expr_1, y = 0), color = 'red', shape = 1)
    }

    pl <- pl + ggplot2::geom_point(data = subDT, aes(x = all_cell_expr_1, y = all_cell_expr_2), color = 'blue', size = 2)
    pl <- pl + ggplot2::geom_point(data = subDT, aes(x = cell_expr_1, y = cell_expr_2), color = 'red', size = 2)
    pl <- pl + ggplot2::geom_segment(data = subDT, aes(x = all_cell_expr_1, xend = cell_expr_1,
                                              y = all_cell_expr_2, yend = cell_expr_2), linetype = 2)
    pl <- pl + ggplot2::labs(x = 'gene in celltype 1', y = 'gene in celltype 2')
    pl <- pl + ggplot2::facet_wrap(~genes+interaction, nrow = facet_nrow, ncol = facet_ncol,
                          scales = facet_scales)

  } else {

    simple_plot_facet = match.arg(arg = simple_plot_facet, choices = c('interaction', 'genes'))

    if(simple_plot_facet == 'interaction') {
      pl <- ggplot2::ggplot()
      pl <- pl + ggplot2::theme_bw()
      pl <- pl + ggplot2::geom_segment(data = subDT, aes(x = all_comb_expr, xend = comb_expr,
                                                y = genes, yend = genes), linetype = 2)
      pl <- pl + ggplot2::geom_point(data = subDT, aes(x = all_comb_expr, y = genes), color = 'blue')
      pl <- pl + ggplot2::geom_point(data = subDT, aes(x = comb_expr, y = genes), color = 'red')
      pl <- pl + ggplot2::facet_wrap(~interaction, scales = facet_scales)
      pl <- pl + ggplot2::labs(x = 'interactions', y = 'genes')
    } else {
      pl <- ggplot2::ggplot()
      pl <- pl + ggplot2::theme_bw()
      pl <- pl + ggplot2::geom_segment(data = subDT, aes(x = all_comb_expr, xend = comb_expr,
                                                y = interaction, yend = interaction), linetype = 2)
      pl <- pl + ggplot2::geom_point(data = subDT, aes(x = all_comb_expr, y = interaction), color = 'blue')
      pl <- pl + ggplot2::geom_point(data = subDT, aes(x = comb_expr, y = interaction), color = 'red')
      pl <- pl + ggplot2::facet_wrap(~genes, scales = facet_scales)
      pl <- pl + ggplot2::labs(x = 'genes', y = 'interactions')
    }
  }

  if(show_plot == TRUE) {
    print(pl)
  }

  return(pl)

}














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
#' @param show_network show underlying spatial network
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
#' @examples
#'     cellProximitySpatPlot2D(gobject)
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
                                    spatial_network_name = 'spatial_network',
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
        cell_color_code = Giotto:::getDistinctColors(n = number_colors)
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
#' @param interaction_name cell-cell interaction name
#' @param cluster_column cluster column with cell clusters
#' @param sdimx x-axis dimension name (default = 'sdimx')
#' @param sdimy y-axis dimension name (default = 'sdimy')
#' @param cell_color color for cells (see details)
#' @param cell_color_code named vector with colors
#' @param color_as_factor convert color column to factor
#' @param show_other_cells decide if show cells not in network
#' @param show_network show underlying spatial network
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
#' @seealso  \code{\link{cellProximitySpatPlot2D}} and \code{\link{cellProximitySpatPlot3D}} for 3D
#' @examples
#'     cellProximitySpatPlot(gobject)
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
#' @param show_network show underlying spatial network
#' @param network_color color of spatial network
#' @param spatial_network_name name of spatial network to use
#' @param show_grid show spatial grid
#' @param grid_color color of spatial grid
#' @param spatial_grid_name name of spatial grid to use
#' @param show_legend show legend
#' @param point_size_select size of selected points
#' @param point_size_other size of other points
#' @param show_plot show plots
#' @param return_plot return plotly object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters from \code{\link{all_plots_save_function}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return plotly
#' @details Description of parameters.
#' @export
#' @examples
#'     cellProximitySpatPlot3D(gobject)

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
                                   spatial_network_name = 'spatial_network',
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

