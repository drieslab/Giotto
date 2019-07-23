

#' @title visPlot
#' @name visPlot
#' @description Visualize cells according to spatial coordinates
#' @param gobject giotto object
#' @param sdimx x-axis dimension name (default = 'sdimx')
#' @param sdimy y-axis dimension name (default = 'sdimy')
#' @param sdimz z-axis dimension name (default = 'sdimz')
#' @param point_size size of point (cell)
#' @param point.border.col color of border around points
#' @param point.border.stroke stroke size of border around points
#' @param cell_color color for cells (see details)
#' @param cell_color_code named vector with colors
#' @param color_as_factor convert color column to factor
#' @param select_cell_groups select a subset of the groups from cell_color
#' @param show_network show underlying spatial network
#' @param network_color color of spatial network
#' @param spatial_network_name name of spatial network to use
#' @param show_grid show spatial grid
#' @param grid_color color of spatial grid
#' @param spatial_grid_name name of spatial grid to use
#' @param coord_fix_ratio fix ratio between x and y-axis
#' @param title title of plot
#' @param show.legend show legend
#' @param show.plot show plot
#' @return ggplot
#' @details Description of parameters.
#' @export
#' @examples
#'     visPlot(gobject)
visPlot <- function(gobject,
                    sdimx = NULL,
                    sdimy = NULL,
                    sdimz = NULL,
                    point_size = 1,
                    point.border.col = 'black',
                    point.border.stroke = 0.1,
                    cell_color = NULL,
                    cell_color_code = NULL,
                    color_as_factor = T,
                    select_cell_groups = NULL,
                    show_network = F,
                    network_color = NULL,
                    spatial_network_name = 'spatial_network',
                    show_grid = F,
                    grid_color = NULL,
                    spatial_grid_name = 'spatial_grid',
                    coord_fix_ratio = 0.6,
                    title = '',
                    show.legend = T,
                    show.plot = F) {

  cell_locations  = gobject@spatial_locs

  ## extract spatial network
  if(!is.null(spatial_network_name)) {
    spatial_network = gobject@spatial_network[[spatial_network_name]]
  } else {
    spatial_network = NULL
  }

  ## extract spatial grid
  if(!is.null(spatial_grid_name)) {
    spatial_grid    = gobject@spatial_grid[[spatial_grid_name]]
  } else {
    spatial_grid = NULL
  }

  cell_metadata   = gobject@cell_metadata
  cell_metadata   = cell_metadata[, !grepl('cell_ID', colnames(cell_metadata)), with = F]

  if(nrow(cell_metadata) == 0) {
    cell_locations_metadata = cell_locations
  } else {
    cell_locations_metadata <- cbind(cell_locations, cell_metadata)
  }

  # create subsets of needed
  if(!is.null(select_cell_groups)) {
    cell_locations_metadata_selected = cell_locations_metadata[get(cell_color) %in% select_cell_groups]
    cell_locations_metadata_other = cell_locations_metadata[!get(cell_color) %in% select_cell_groups]
  }

  # first 2 dimensions need to be defined
  if(is.null(sdimx) | is.null(sdimy)) {
    # cat('first and second dimenion need to be defined, default is first 2 \n')
    sdimx = 'sdimx'
    sdimy = 'sdimy'
  }



  # if 3 dimensions are defined create a 3D plot
  if(!is.null(sdimx) & !is.null(sdimy) & !is.null(sdimz)) {

    cat('create 3D plot')

    if(!is.null(cell_color) & cell_color %in% colnames(cell_locations_metadata)) {
      if(is.null(cell_color_code)) cell_color_code <- 'lightblue'
      p <- plotly::plot_ly(type = 'scatter3d',
                   x = cell_locations_metadata$sdimx, y = cell_locations_metadata$sdimy, z = cell_locations_metadata$sdimz,
                   color = cell_locations_metadata[[cell_color]],
                   mode = 'markers', colors = cell_color_code)
      print(p)
    } else {
      p <- plotly::plot_ly(type = 'scatter3d',
                   x = cell_locations_metadata$sdimx, y = cell_locations_metadata$sdimy, z = cell_locations_metadata$sdimz,
                   mode = 'markers', colors = 'lightblue')

      if(show.plot == TRUE) {
        print(pl)
      }
      return(pl)
    }



  } else {

    pl <- ggplot2::ggplot()
    pl <- pl + ggplot2::theme_bw()

    ## plot spatial network
    if(!is.null(spatial_network) & show_network == TRUE) {
      if(is.null(network_color)) network_color = 'red'
      pl <- pl + ggplot2::geom_segment(data = spatial_network, aes(x = sdimx_begin, y = sdimy_begin,
                                                          xend = sdimx_end, yend = sdimy_end),
                              color = network_color, size = 0.5, alpha = 0.5)
    }

    ## plot spatial grid
    if(!is.null(spatial_grid) & show_grid == TRUE) {
      if(is.null(grid_color)) grid_color = 'black'
      pl <- pl + ggplot2::geom_rect(data = spatial_grid, aes(xmin = x_start, xmax = x_end,
                                                    ymin = y_start, ymax = y_end),
                           color = grid_color, fill = NA)
    }

    # cell color default
    if(is.null(cell_color)) {

      cell_color = 'lightblue'
      pl <- pl + ggplot2::geom_point(data = cell_locations_metadata, aes_string(x = sdimx, y = sdimy),
                            show.legend = show.legend, shape = 21,
                            fill = cell_color, size = point_size,
                            stroke = point.border.stroke, color = point.border.col)

    }

    else if (is.character(cell_color)) {

      if(cell_color %in% colnames(cell_locations_metadata)) {

        if(color_as_factor == TRUE) {
          if(is.null(select_cell_groups)) {
            factor_data = factor(cell_locations_metadata[[cell_color]])
            cell_locations_metadata[[cell_color]] <- factor_data
          } else {
            factor_data_selected = factor(cell_locations_metadata_selected[[cell_color]])
            cell_locations_metadata_selected[[cell_color]] <- factor_data_selected
            factor_data_other = factor(cell_locations_metadata_other[[cell_color]])
            cell_locations_metadata_other[[cell_color]] <- factor_data_other
          }

        }

        if(is.null(select_cell_groups)) {
          pl <- pl + ggplot2::geom_point(data = cell_locations_metadata, aes_string(x = sdimx, y = sdimy, fill = cell_color),
                                show.legend = show.legend, shape = 21, size = point_size,
                                stroke = point.border.stroke, color = point.border.col)
        } else {
          cell_color_other = 'grey'
          pl <- pl + ggplot2::geom_point(data = cell_locations_metadata_other, aes_string(x = sdimx, y = sdimy),
                                fill = cell_color_other,
                                show.legend = show.legend, shape = 21, size = point_size/2,
                                stroke = point.border.stroke, color = point.border.col)

          pl <- pl + ggplot2::geom_point(data = cell_locations_metadata_selected, aes_string(x = sdimx, y = sdimy, fill = cell_color),
                                show.legend = show.legend, shape = 21, size = point_size,
                                stroke = point.border.stroke, color = point.border.col)
        }




        if(!is.null(cell_color_code)) {
          pl <- pl + ggplot2::scale_fill_manual(values = cell_color_code)
        } else if(color_as_factor == T) {
          if(is.null(select_cell_groups)) {
            number_colors = length(unique(factor_data))
            cell_color_code = getDistinctColors(n = number_colors)
            names(cell_color_code) = unique(factor_data)
          } else {
            number_colors = length(unique(factor_data_selected))
            cell_color_code = getDistinctColors(n = number_colors)
            names(cell_color_code) = unique(factor_data_selected)
          }
          pl <- pl + ggplot2::scale_fill_manual(values = cell_color_code)
        } else if(color_as_factor == F){
          pl <- pl + ggplot2::scale_fill_gradient(low = 'blue', high = 'red')
        }

      } else {


        if(is.null(select_cell_groups)) {
          pl <- pl + ggplot2::geom_point(data = cell_locations_metadata, aes_string(x = sdimx, y = sdimy),
                                show.legend = show.legend, shape = 21, fill = cell_color,
                                size = point_size,
                                stroke = point.border.stroke, color = point.border.col)
        } else {
          cell_color_other = 'grey'
          pl <- pl + ggplot2::geom_point(data = cell_locations_metadata_other, aes_string(x = sdimx, y = sdimy),
                                show.legend = show.legend, shape = 21, fill = cell_color_other,
                                size = point_size/2,
                                stroke = point.border.stroke, color = point.border.col)

          pl <- pl + ggplot2::geom_point(data = cell_locations_metadata_selected, aes_string(x = sdimx, y = sdimy),
                                show.legend = show.legend, shape = 21, fill = cell_color,
                                size = point_size,
                                stroke = point.border.stroke, color = point.border.col)
        }

      }

    }

    pl <- pl + ggplot2::theme(plot.title = element_text(hjust = 0.5),
                     legend.title = element_text(size = 10),
                     legend.text = element_text(size = 10))

    # fix coord ratio
    if(!is.null(coord_fix_ratio)) {
      pl <- pl + ggplot2::coord_fixed(ratio = coord_fix_ratio)
    }

    pl <- pl + ggplot2::labs(x = 'x coordinates', y = 'y coordinates', title = title)


    if(show.plot == TRUE) {
      print(pl)
    }
    return(pl)
  }

}




#' @title visGenePlot
#' @name visGenePlot
#' @description Visualize cells and gene expression according to spatial coordinates
#' @param gobject giotto object
#' @param expression_values gene expression values to use
#' @param genes genes to show
#' @param show_network show underlying spatial network
#' @param network_color color of spatial network
#' @param spatial_network_name name of spatial network to use
#' @param show_grid show spatial grid
#' @param grid_color color of spatial grid
#' @param spatial_grid_name name of spatial grid to use
#' @param midpoint expression midpoint
#' @param scale_alpha_with_expression scale expression with ggplot alpha parameter
#' @param point_size size of point (cell)
#' @param point.border.col color of border around points
#' @param point.border.stroke stroke size of border around points
#' @param cow_n_col cowplot param: how many columns
#' @param cow_rel_h cowplot param: relative height
#' @param cow_rel_w cowplot param: relative width
#' @param cow_align cowplot param: how to align
#' @param show.legend show legend
#' @param show_plots show plots
#' @return ggplot
#' @details Description of parameters.
#' @export
#' @examples
#'     visGenePlot(gobject)
visGenePlot <- function(gobject,
                        expression_values = c('normalized', 'scaled', 'custom'),
                        genes,
                        show_network = F,
                        network_color = NULL,
                        spatial_network_name = 'spatial_network',
                        show_grid = F,
                        grid_color = NULL,
                        spatial_grid_name = 'spatial_grid',
                        midpoint = 0,
                        scale_alpha_with_expression = TRUE,
                        point_size = 1,
                        point.border.col = 'black',
                        point.border.stroke = 0.1,
                        show.legend = T,
                        cow_n_col = 2,
                        cow_rel_h = 1,
                        cow_rel_w = 1,
                        cow_align = 'h',
                        show_plots = F) {

  selected_genes = genes

  values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))
  expr_values = select_expression_values(gobject = gobject, values = values)

  # only keep genes that are in the dataset
  selected_genes = selected_genes[selected_genes %in% rownames(expr_values) ]

  # get selected gene expression values in data.table format
  if(length(selected_genes) == 1) {
    subset_expr_data = expr_values[rownames(expr_values) %in% selected_genes, ]
    t_sub_expr_data_DT = data.table::data.table('selected_gene' = subset_expr_data, 'cell_ID' = colnames(expr_values))
    data.table::setnames(t_sub_expr_data_DT, 'selected_gene', selected_genes)
  } else {
    subset_expr_data = expr_values[rownames(expr_values) %in% selected_genes, ]
    t_sub_expr_data = t(subset_expr_data)
    t_sub_expr_data_DT = data.table::as.data.table(t_sub_expr_data)
    t_sub_expr_data_DT[, cell_ID := rownames(t_sub_expr_data)]
  }


  ## extract cell locations
  cell_locations  = gobject@spatial_locs

  ## extract spatial network
  if(!is.null(spatial_network_name)) {
    spatial_network = gobject@spatial_network[[spatial_network_name]]
  } else {
    spatial_network = NULL
  }

  ## extract spatial grid
  if(!is.null(spatial_grid_name)) {
    spatial_grid    = gobject@spatial_grid[[spatial_grid_name]]
  } else {
    spatial_grid = NULL
  }

  ## extract cell metadata
  cell_metadata   = gobject@cell_metadata
  cell_metadata   = cell_metadata[, !grepl('cell_ID', colnames(cell_metadata)), with = F]

  if(nrow(cell_metadata) == 0) {
    cell_locations_metadata = cell_locations
  } else {
    cell_locations_metadata <- cbind(cell_locations, cell_metadata)
  }

  cell_locations_metadata_genes <- merge(cell_locations_metadata, t_sub_expr_data_DT, by = 'cell_ID')

  ## plotting ##
  savelist <- list()
  for(gene in selected_genes) {

    pl <- ggplot2::ggplot()
    pl <- pl + ggplot2::theme_classic()

    ## plot spatial network
    if(!is.null(spatial_network) & show_network == TRUE) {
      if(is.null(network_color)) network_color = 'red'
      pl <- pl + ggplot2::geom_segment(data = spatial_network, aes(x = sdimx_begin, y = sdimy_begin,
                                                          xend = sdimx_end, yend = sdimy_end),
                              color = network_color, size = 0.5, alpha = 0.5)
    }

    ## plot spatial grid
    if(!is.null(spatial_grid) & show_grid == TRUE) {
      if(is.null(grid_color)) grid_color = 'black'
      pl <- pl + ggplot2::geom_rect(data = spatial_grid, aes(xmin = x_start, xmax = x_end,
                                                    ymin = y_start, ymax = y_end),
                           color = grid_color, fill = NA)
    }


    if(scale_alpha_with_expression == TRUE) {
      pl <- pl + ggplot2::geom_point(data = cell_locations_metadata_genes, aes_string(x = 'sdimx', y = 'sdimy',
                                                                             fill = gene, alpha = gene),
                            shape = 21,
                            color = point.border.col, size = point_size, stroke = point.border.stroke,
                            show.legend = show.legend)
    } else {
      pl <- pl + ggplot2::geom_point(data = cell_locations_metadata_genes, aes_string(x = 'sdimx', y = 'sdimy',
                                                                             fill = gene),
                            shape = 21,
                            color = point.border.col, size = point_size, stroke = point.border.stroke,
                            show.legend = show.legend)
    }
    pl <- pl + ggplot2::scale_alpha_continuous(guide = 'none')
    pl <- pl + ggplot2::scale_fill_gradient2(low = 'darkblue', mid = 'white', high = 'darkred',
                                    midpoint = midpoint, guide = guide_colorbar(title = ''))
    pl <- pl + ggplot2::labs(x = 'coord x', y = 'coord y', title = gene)
    pl <- pl + ggplot2::theme(plot.title = element_text(hjust = 0.5))

    if(show_plots == TRUE) {
      print(pl)
    }

    savelist[[gene]] <- pl
  }

  # combine plots with cowplot
  combo_plot <- cowplot::plot_grid(plotlist = savelist,
                                   ncol = cow_n_col,
                                   rel_heights = cow_rel_h, rel_widths = cow_rel_w, align = cow_align)
  combined_cowplot = cowplot::plot_grid(combo_plot)

  return(combined_cowplot)

}




#' @title visDimPlot
#' @name visDimPlot
#' @description Visualize cells according to dimension reduction coordinates
#' @param gobject giotto object
#' @param dim_reduction_to_use dimension reduction to use
#' @param dim_reduction_name dimension reduction name
#' @param dim1_to_use dimension to use on x-axis
#' @param dim2_to_use dimension to use on y-axis
#' @param show_NN_network show underlying NN network
#' @param nn_network_to_use type of NN network to use (kNN vs sNN)
#' @param network_name name of NN network to use, if show_NN_network = TRUE
#' @param cell_color color for cells (see details)
#' @param color_as_factor convert color column to factor
#' @param cell_color_code named vector with colors
#' @param show_cluster_center plot center of selected clusters
#' @param show_center_label plot label of selected clusters
#' @param center.point.size size of center points
#' @param label.size  size of labels
#' @param label.fontface font of labels
#' @param edge_alpha column to use for alpha of the edges
#' @param point_size size of point (cell)
#' @param point.border.col color of border around points
#' @param point.border.stroke stroke size of border around points
#' @param show.legend show legend
#' @return ggplot
#' @details Description of parameters.
#' @export
#' @examples
#'     visDimPlot(gobject)
visDimPlot <- function(gobject,
                       dim_reduction_to_use = 'umap',
                       dim_reduction_name = 'umap',
                       dim1_to_use = 1,
                       dim2_to_use = 2,
                       show_NN_network = F,
                       nn_network_to_use = 'sNN',
                       network_name = 'sNN.pca',
                       cell_color = NULL,
                       color_as_factor = T,
                       cell_color_code = NULL,
                       show_cluster_center = F,
                       show_center_label = T,
                       center.point.size = 4,
                       center.point.border.col = 'black',
                       center.point.border.stroke = 0.1,
                       label.size = 4,
                       label.fontface = 'bold',
                       edge_alpha = NULL,
                       point_size = 1,
                       point.border.col = 'black',
                       point.border.stroke = 0.1,
                       show.legend = T) {


  ## dimension reduction ##
  dim_dfr = gobject@dimension_reduction$cells[[dim_reduction_to_use]][[dim_reduction_name]]$coordinates[,c(dim1_to_use, dim2_to_use)]
  dim_names = colnames(dim_dfr)
  dim_DT = data.table::as.data.table(dim_dfr); dim_DT[, cell_ID := rownames(dim_dfr)]

  ## annotated cell metadata
  cell_metadata = gobject@cell_metadata
  annotated_DT = merge(cell_metadata, dim_DT, by = 'cell_ID')


  # create input for network
  if(show_NN_network == TRUE) {

    # nn_network
    selected_nn_network = gobject@nn_network[[nn_network_to_use]][[network_name]][['igraph']]
    network_DT = data.table::as.data.table(igraph::as_data_frame(selected_nn_network, what = 'edges'))

    # annotated network
    old_dim_names = dim_names

    annotated_network_DT <- merge(network_DT, dim_DT, by.x = 'from', by.y = 'cell_ID')
    from_dim_names = paste0('from_', old_dim_names)
    data.table::setnames(annotated_network_DT, old = old_dim_names, new = from_dim_names)

    annotated_network_DT <- merge(annotated_network_DT, dim_DT, by.x = 'to', by.y = 'cell_ID')
    to_dim_names = paste0('to_', old_dim_names)
    data.table::setnames(annotated_network_DT, old = old_dim_names, new = to_dim_names)

  }


  # visualize

  pl <- ggplot2::ggplot()
  pl <- pl + ggplot2::theme_classic()

  # network layer
  if(show_NN_network == TRUE) {

    if(is.null(edge_alpha)) {
      edge_alpha = 0.5
      pl <- pl + ggplot2::geom_segment(data = annotated_network_DT, aes_string(x = from_dim_names[1], y = from_dim_names[2],
                                                                      xend = to_dim_names[1], yend = to_dim_names[2]), alpha = edge_alpha, show.legend = show.legend)
    } else if(is.numeric(edge_alpha)) {
      pl <- pl + ggplot2::geom_segment(data = annotated_network_DT, aes_string(x = from_dim_names[1], y = from_dim_names[2],
                                                                      xend = to_dim_names[1], yend = to_dim_names[2]), alpha = edge_alpha, show.legend = show.legend)
    } else if(is.character(edge_alpha)) {

      if(edge_alpha %in% colnames(annotated_network_DT)) {
        pl <- pl + ggplot2::geom_segment(data = annotated_network_DT, aes_string(x = from_dim_names[1], y = from_dim_names[2],
                                                                        xend = to_dim_names[1], yend = to_dim_names[2], alpha = edge_alpha), show.legend = show.legend)
      }
    }
  }

  # point layer
  if(is.null(cell_color)) {
    cell_color = 'lightblue'
    pl <- pl + ggplot2::geom_point(data = annotated_DT, aes_string(x = dim_names[1], dim_names[2]),
                          color = cell_color, show.legend = show.legend, size = point_size)

  } else if (is.character(cell_color)) {

    if(cell_color %in% colnames(annotated_DT)) {

      class_cell_color = class(annotated_DT[[cell_color]])



      # convert numericals to factors
      if(color_as_factor == TRUE) {
        factor_data = factor(annotated_DT[[cell_color]])
        annotated_DT[[cell_color]] <- factor_data
        # for centers
        if(show_cluster_center == TRUE | show_center_label == TRUE) {
          annotated_DT_centers = annotated_DT[, .(center_1 = median(get(dim_names[1])), center_2 = median(get(dim_names[2]))), by = cell_color]
          factor_center_data = factor(annotated_DT_centers[[cell_color]])
          annotated_DT_centers[[cell_color]] <- factor_center_data
        }
      } else {

        # TEST: centers can only be shown for factors that are part of the metadata
        if((show_cluster_center == TRUE | show_center_label == TRUE) & class_cell_color %in% c('character', 'factor')) {
          annotated_DT_centers = annotated_DT[, .(center_1 = median(get(dim_names[1])), center_2 = median(get(dim_names[2]))), by = cell_color]
        }

      }

      pl <- pl + ggplot2::geom_point(data = annotated_DT, aes_string(x = dim_names[1], y = dim_names[2], fill = cell_color),
                            show.legend = show.legend, shape = 21, size = point_size,
                            color = point.border.col, stroke = point.border.stroke)

      # plot centers
      if(show_cluster_center == TRUE & (color_as_factor == TRUE | class_cell_color %in% c('character', 'factor'))) {

        pl <- pl + ggplot2::geom_point(data = annotated_DT_centers,
                              aes_string(x = 'center_1', y = 'center_2', fill = cell_color),
                              color = center.point.border.col, stroke = center.point.border.stroke,
                              size = center.point.size, shape = 21)
      }

      # plot labels
      if(show_center_label == TRUE) {
        pl <- pl + ggrepel::geom_text_repel(data = annotated_DT_centers,
                                            aes_string(x = 'center_1', y = 'center_2', label = cell_color),
                                            size = label.size, fontface = label.fontface)
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

    }

  } else {
    pl <- pl + ggplot2::geom_point(data = annotated_DT, aes_string(x = dim_names[1], y = dim_names[2]),
                          show.legend = show.legend, shape = 21, fill = cell_color,
                          size = point_size,
                          color = point.border.col, stroke = point.border.stroke)
  }


  ## add % variance explained to names of plot for PCA ##
  if(dim_reduction_to_use == 'pca') {

    eigenvaluesDT = data.table::as.data.table(gobject@dimension_reduction$cells[[dim_reduction_to_use]][[dim_reduction_name]]$misc$eig)
    var_expl_vec = eigenvaluesDT[c(dim1_to_use, dim2_to_use)][['percentage of variance']]
    dim1_x_variance = var_expl_vec[1]
    dim2_y_variance = var_expl_vec[2]

    x_name = paste0('pca','-',dim_names[1])
    y_name = paste0('pca','-',dim_names[2])

    x_title = sprintf('%s explains %.02f%% of variance', x_name, var_expl_vec[1])
    y_title = sprintf('%s explains %.02f%% of variance', y_name, var_expl_vec[2])

    pl <- pl + ggplot2::labs(x = x_title, y = y_title)

  } else {

    x_title = paste0(dim_reduction_to_use,'-',dim_names[1])
    y_title = paste0(dim_reduction_to_use,'-',dim_names[2])

    pl <- pl + ggplot2::labs(x = x_title, y = y_title)

  }


  return(pl)

}



#' @title plotUMAP
#' @name plotUMAP
#' @description Short wrapper for UMAP visualization
#' @param gobject giotto object
#' @param ... other parameters that are part of visDimPlot()
#' @return ggplot
#' @details Description of parameters.
#' @export
#' @seealso \code{\link{visDimPlot}}
#' @examples
#'     plotUMAP(gobject)
plotUMAP = function(gobject, ...) {

  visDimPlot(gobject = gobject, dim_reduction_to_use = 'umap', dim_reduction_name = 'umap', ...)

}

#' @title plotTSNE
#' @name plotTSNE
#' @description Short wrapper for tSNE visualization
#' @param gobject giotto object
#' @param ... other parameters that are part of visDimPlot()
#' @return ggplot
#' @details Description of parameters.
#' @export
#' @seealso \code{\link{visDimPlot}}
#' @examples
#'     plotTSNE(gobject)
plotTSNE = function(gobject, ...) {

  visDimPlot(gobject = gobject, dim_reduction_to_use = 'tsne', dim_reduction_name = 'tsne', ...)

}

#' @title plotPCA
#' @name plotPCA
#' @description Short wrapper for PCA visualization
#' @param gobject giotto object
#' @param ... other parameters that are part of visDimPlot()
#' @return ggplot
#' @details Description of parameters.
#' @export
#' @seealso \code{\link{visDimPlot}}
#' @examples
#'     plotPCA(gobject)
plotPCA = function(gobject, ...) {

  visDimPlot(gobject = gobject, dim_reduction_to_use = 'pca', dim_reduction_name = 'pca', ...)

}




#' @title visForceLayoutPlot
#' @name visForceLayoutPlot
#' @description Visualize cells according to forced layout algorithm coordinates
#' @param gobject giotto object
#' @param nn_network_to_use type of NN network to use (kNN vs sNN)
#' @param network_name NN network to use
#' @param layout_name name of layout to use
#' @param dim1_to_use dimension to use on x-axis
#' @param dim2_to_use dimension to use on y-axis
#' @param show_NN_network show underlying NN network
#' @param cell_color color for cells (see details)
#' @param color_as_factor convert color column to factor
#' @param cell_color_code named vector with colors
#' @param edge_alpha column to use for alpha of the edges
#' @param point_size size of point (cell)
#' @param point.border.col color of border around points
#' @param point.border.stroke stroke size of border around points
#' @param show.legend show legend
#' @return ggplot
#' @details Description of parameters.
#' @export
#' @examples
#'     visForceLayoutPlot(gobject)
visForceLayoutPlot <- function(gobject,
                               nn_network_to_use = 'sNN',
                               network_name = 'sNN.pca',
                               layout_name = 'layout',
                               dim1_to_use = 1,
                               dim2_to_use = 2,
                               show_NN_network = T,
                               cell_color = NULL,
                               color_as_factor = F,
                               cell_color_code = NULL,
                               edge_alpha = NULL,
                               point_size = 1,
                               point.border.col = 'black',
                               point.border.stroke = 0.1,
                               show.legend = T) {


  ## layout ##
  co = gobject@nn_network[[nn_network_to_use]][[network_name]][['layout']][,c(dim1_to_use, dim2_to_use)]
  dim_dfr = as.data.frame(co)
  colnames(dim_dfr) = paste0('Dim.', 1:ncol(dim_dfr))
  dim_names = colnames(dim_dfr)
  dim_DT = data.table::as.data.table(dim_dfr); dim_DT[, cell_ID := paste0('cell_', rownames(dim_dfr))]

  ## annotated cell metadata
  cell_metadata = gobject@cell_metadata
  annotated_DT = merge(cell_metadata, dim_DT, by = 'cell_ID')

  # create input for network
  if(show_NN_network == TRUE) {

    # nn_network
    selected_nn_network = gobject@nn_network[[nn_network_to_use]][[network_name]][['igraph']]
    network_DT = data.table::as.data.table(igraph::as_data_frame(selected_nn_network, what = 'edges'))
    #network_DT[, from := paste0('cell_', from)]
    #network_DT[, to := paste0('cell_', to)]

    # annotated network
    old_dim_names = dim_names

    annotated_network_DT <- merge(network_DT, dim_DT, by.x = 'from', by.y = 'cell_ID')
    from_dim_names = paste0('from_', old_dim_names)
    data.table::setnames(annotated_network_DT, old = old_dim_names, new = from_dim_names)

    annotated_network_DT <- merge(annotated_network_DT, dim_DT, by.x = 'to', by.y = 'cell_ID')
    to_dim_names = paste0('to_', old_dim_names)
    data.table::setnames(annotated_network_DT, old = old_dim_names, new = to_dim_names)

  }


  # visualize

  pl <- ggplot2::ggplot()
  pl <- pl + ggplot2::theme_classic()

  # network layer
  if(show_NN_network == TRUE) {

    if(is.null(edge_alpha)) {
      edge_alpha = 0.5
      pl <- pl + ggplot2::geom_segment(data = annotated_network_DT, aes_string(x = from_dim_names[1], y = from_dim_names[2],
                                                                      xend = to_dim_names[1], yend = to_dim_names[2]), alpha = edge_alpha, show.legend = show.legend)
    } else if(is.numeric(edge_alpha)) {
      pl <- pl + ggplot2::geom_segment(data = annotated_network_DT, aes_string(x = from_dim_names[1], y = from_dim_names[2],
                                                                      xend = to_dim_names[1], yend = to_dim_names[2]), alpha = edge_alpha, show.legend = show.legend)
    } else if(is.character(edge_alpha)) {

      if(edge_alpha %in% colnames(annotated_network_DT)) {
        pl <- pl + ggplot2::geom_segment(data = annotated_network_DT, aes_string(x = from_dim_names[1], y = from_dim_names[2],
                                                                        xend = to_dim_names[1], yend = to_dim_names[2], alpha = edge_alpha), show.legend = show.legend)
      }
    }
  }

  # point layer
  if(is.null(cell_color)) {
    cell_color = 'lightblue'
    pl <- pl + ggplot2::geom_point(data = annotated_DT, aes_string(x = dim_names[1], dim_names[2]),
                          color = cell_color, show.legend = show.legend, size = point_size)

  } else if (is.character(cell_color)) {

    if(cell_color %in% colnames(annotated_DT)) {

      if(color_as_factor == TRUE) {
        factor_data = factor(annotated_DT[[cell_color]])
        annotated_DT[[cell_color]] <- factor_data
      }

      pl <- pl + ggplot2::geom_point(data = annotated_DT, aes_string(x = dim_names[1], y = dim_names[2], fill = cell_color),
                            show.legend = show.legend, shape = 21, size = point_size,
                            color = point.border.col, stroke = point.border.stroke)


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

    }

  } else {
    pl <- pl + ggplot2::geom_point(data = annotated_DT, aes_string(x = dim_names[1], y = dim_names[2]),
                          show.legend = show.legend, shape = 21, fill = cell_color,
                          size = point_size,
                          color = point.border.col, stroke = point.border.stroke)
  }

  return(pl)

}




#' @title visSpatDimPlot
#' @name visSpatDimPlot
#' @description Visualize cells according to spatial AND dimension reduction coordinates
#' @param gobject giotto object
#' @param plot_alignment direction to align plot
#' @param dim_reduction_to_use dimension reduction to use
#' @param dim_reduction_name dimension reduction name
#' @param dim1_to_use dimension to use on x-axis
#' @param dim2_to_use dimension to use on y-axis
#' @param show_NN_network show underlying NN network
#' @param nn_network_to_use type of NN network to use (kNN vs sNN)
#' @param network_name name of NN network to use, if show_NN_network = TRUE
#' @param cell_color color for cells (see details)
#' @param color_as_factor convert color column to factor
#' @param cell_color_code named vector with colors
#' @param edge_alpha column to use for alpha of the edges
#' @param show_spatial_network show spatial network
#' @param spatial_network_name name of spatial network to use
#' @param spatial_network_color color of spatial network
#' @param show_spatial_grid show spatial grid
#' @param spatial_grid_name name of spatial grid to use
#' @param spatial_grid_color color of spatial grid
#' @param show.legend show legend
#' @param show.plot show plot
#' @return ggplot
#' @details Description of parameters.
#' @export
#' @examples
#'     visSpatDimPlot(gobject)
visSpatDimPlot <- function(gobject,
                           plot_alignment = c('vertical', 'horizontal'),
                           dim_reduction_to_use = 'umap',
                           dim_reduction_name = 'umap',
                           dim1_to_use = 1,
                           dim2_to_use = 2,
                           show_NN_network = F,
                           nn_network_to_use = 'sNN',
                           network_name = 'sNN.pca',
                           show_cluster_center = F,
                           show_center_label = T,
                           center.point.size = 4,
                           label.size = 4,
                           label.fontface = 'bold',
                           cell_color = NULL,
                           color_as_factor = T,
                           cell_color_code = NULL,
                           dim_point_size = 1,
                           dim.point.border.col = 'black',
                           dim.point.border.stroke = 0.1,
                           edge_alpha = 0.05,
                           show_spatial_network = F,
                           spatial_network_name = 'spatial_network',
                           spatial_network_color = NULL,
                           show_spatial_grid = F,
                           spatial_grid_name = 'spatial_grid',
                           spatial_grid_color = NULL,
                           spatial_point_size = 1,
                           spatial.point.border.col = 'black',
                           spatial.point.border.stroke = 0.1,
                           show.legend = T,
                           show.plot = F) {

  plot_alignment = match.arg(plot_alignment, choices = c('vertical', 'horizontal'))


  # create matching cell_color_code
  if(is.null(cell_color_code)) {
    if(is.character(cell_color)) {

      cell_metadata = pDataDT(gobject)
      if(cell_color %in% colnames(cell_metadata)) {

        if(color_as_factor == TRUE) {
          number_colors = length(unique(cell_metadata[[cell_color]]))
          cell_color_code = getDistinctColors(n = number_colors)
          names(cell_color_code) = unique(cell_metadata[[cell_color]])
          cell_color_code = cell_color_code
        }
      }
    }
  }


  # dimension reduction plot
  dmpl = visDimPlot(gobject = gobject, dim_reduction_to_use = dim_reduction_to_use, dim_reduction_name = dim_reduction_name,
                    dim1_to_use = dim1_to_use, dim2_to_use = dim2_to_use,
                    show_NN_network = show_NN_network,
                    nn_network_to_use = nn_network_to_use, network_name = network_name,
                    cell_color = cell_color, color_as_factor = color_as_factor, cell_color_code = cell_color_code,
                    edge_alpha = edge_alpha, show.legend = show.legend,
                    point_size = dim_point_size,
                    point.border.col = dim.point.border.col,
                    point.border.stroke = dim.point.border.stroke,
                    show_cluster_center = show_cluster_center,
                    show_center_label = show_center_label,
                    center.point.size = center.point.size,
                    label.size = label.size,
                    label.fontface = label.fontface)

  # spatial plot
  spl = visPlot(gobject = gobject,
                show_network = show_spatial_network, spatial_network_name = spatial_network_name,
                show_grid = show_spatial_grid, spatial_grid_name = spatial_grid_name,
                cell_color = cell_color,
                cell_color_code = cell_color_code, color_as_factor = color_as_factor, coord_fix_ratio = NULL,
                network_color = spatial_network_color, grid_color = spatial_grid_color,
                show.legend = show.legend, show.plot = show.plot,
                point_size = spatial_point_size, point.border.col = spatial.point.border.col,
                point.border.stroke = spatial.point.border.stroke)


  if(plot_alignment == 'vertical') {
    combo_plot <- cowplot::plot_grid(dmpl, spl, ncol = 1, rel_heights = c(1), rel_widths = c(1), align = 'v')
    return(cowplot::plot_grid(combo_plot))
  } else {
    combo_plot <- cowplot::plot_grid(dmpl, spl, ncol = 2, rel_heights = c(1), rel_widths = c(1), align = 'h')
    return(cowplot::plot_grid(combo_plot))
  }
}



#' @title visDimGenePlot
#' @name visDimGenePlot
#' @description Visualize cells and gene expression according to dimension reduction coordinates
#' @param gobject giotto object
#' @param expression_values gene expression values to use
#' @param genes genes to show
#' @param dim_reduction_to_use dimension reduction to use
#' @param dim_reduction_name dimension reduction name
#' @param dim1_to_use dimension to use on x-axis
#' @param dim2_to_use dimension to use on y-axis
#' @param show_NN_network show underlying NN network
#' @param nn_network_to_use type of NN network to use (kNN vs sNN)
#' @param network_name name of NN network to use, if show_NN_network = TRUE
#' @param edge_alpha column to use for alpha of the edges
#' @param scale_alpha_with_expression scale expression with ggplot alpha parameter
#' @param point_size size of point (cell)
#' @param point.border.col color of border around points
#' @param point.border.stroke stroke size of border around points
#' @param midpoint size of point (cell)
#' @param cow_n_col cowplot param: how many columns
#' @param cow_rel_h cowplot param: relative height
#' @param cow_rel_w cowplot param: relative width
#' @param cow_align cowplot param: how to align
#' @param show.legend show legend
#' @param show_plots show plots
#' @return ggplot
#' @details Description of parameters.
#' @export
#' @examples
#'     visDimGenePlot(gobject)
visDimGenePlot <- function(gobject,
                           expression_values = c('normalized', 'scaled', 'custom'),
                           genes = NULL,
                           dim_reduction_to_use = 'umap',
                           dim_reduction_name = 'umap',
                           dim1_to_use = 1,
                           dim2_to_use = 2,
                           show_NN_network = T,
                           nn_network_to_use = 'sNN',
                           network_name = 'sNN.pca',
                           edge_alpha = NULL,
                           scale_alpha_with_expression = TRUE,
                           point_size = 1,
                           point.border.col = 'black',
                           point.border.stroke = 0.1,
                           midpoint = 0,
                           cow_n_col = 2,
                           cow_rel_h = 1,
                           cow_rel_w = 1,
                           cow_align = 'h',
                           show.legend = T,
                           show_plots = T) {


  ## select genes ##
  selected_genes = genes
  values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))
  expr_values = select_expression_values(gobject = gobject, values = values)

  # only keep genes that are in the dataset
  selected_genes = selected_genes[selected_genes %in% rownames(expr_values) ]

  #
  if(length(selected_genes) == 1) {
    subset_expr_data = expr_values[rownames(expr_values) %in% selected_genes, ]
    t_sub_expr_data_DT = data.table('selected_gene' = subset_expr_data, 'cell_ID' = colnames(expr_values))
    data.table::setnames(t_sub_expr_data_DT, 'selected_gene', selected_genes)
  } else {
    subset_expr_data = expr_values[rownames(expr_values) %in% selected_genes, ]
    t_sub_expr_data = t(subset_expr_data)
    t_sub_expr_data_DT = data.table::as.data.table(t_sub_expr_data)
    t_sub_expr_data_DT[, cell_ID := rownames(t_sub_expr_data)]
  }


  ## dimension reduction ##
  dim_dfr = gobject@dimension_reduction$cells[[dim_reduction_to_use]][[dim_reduction_name]]$coordinates[,c(dim1_to_use, dim2_to_use)]
  dim_names = colnames(dim_dfr)
  dim_DT = data.table::as.data.table(dim_dfr); dim_DT[, cell_ID := rownames(dim_dfr)]

  ## annotated cell metadata
  cell_metadata = gobject@cell_metadata
  annotated_DT = merge(cell_metadata, dim_DT, by = 'cell_ID')

  ## merge gene info
  annotated_gene_DT = merge(annotated_DT, t_sub_expr_data_DT, by = 'cell_ID')

  # create input for network
  if(show_NN_network == TRUE) {

    # nn_network
    selected_nn_network = gobject@nn_network[[nn_network_to_use]][[network_name]]
    network_DT = data.table::as.data.table(igraph::as_data_frame(selected_nn_network, what = 'edges'))
    network_DT[, from := paste0('cell_', from)]
    network_DT[, to := paste0('cell_', to)]

    # annotated network
    old_dim_names = dim_names

    annotated_network_DT <- merge(network_DT, dim_DT, by.x = 'from', by.y = 'cell_ID')
    from_dim_names = paste0('from_', old_dim_names)
    data.table::setnames(annotated_network_DT, old = old_dim_names, new = from_dim_names)

    annotated_network_DT <- merge(annotated_network_DT, dim_DT, by.x = 'to', by.y = 'cell_ID')
    to_dim_names = paste0('to_', old_dim_names)
    data.table::setnames(annotated_network_DT, old = old_dim_names, new = to_dim_names)

  }



  ## visualize multipe plots ##
  savelist <- list()
  for(gene in selected_genes) {


    ## OLD need to be combined ##
    pl <- ggplot2::ggplot()
    pl <- pl + ggplot2::theme_classic()

    # network layer
    if(show_NN_network == TRUE) {

      if(is.null(edge_alpha)) {
        edge_alpha = 0.5
        pl <- pl + ggplot2::geom_segment(data = annotated_network_DT, aes_string(x = from_dim_names[1], y = from_dim_names[2],
                                                                        xend = to_dim_names[1], yend = to_dim_names[2]), alpha = edge_alpha, show.legend = show.legend)
      } else if(is.numeric(edge_alpha)) {
        pl <- pl + ggplot2::geom_segment(data = annotated_network_DT, aes_string(x = from_dim_names[1], y = from_dim_names[2],
                                                                        xend = to_dim_names[1], yend = to_dim_names[2]), alpha = edge_alpha, show.legend = show.legend)
      } else if(is.character(edge_alpha)) {

        if(edge_alpha %in% colnames(annotated_network_DT)) {
          pl <- pl + ggplot2::geom_segment(data = annotated_network_DT, aes_string(x = from_dim_names[1], y = from_dim_names[2],
                                                                          xend = to_dim_names[1], yend = to_dim_names[2], alpha = edge_alpha), show.legend = show.legend)
        }
      }
    }


    # point layer
    if(is.null(genes)) {
      cell_color = 'lightblue'
      pl <- pl + ggplot2::geom_point(data = annotated_gene_DT, aes_string(x = dim_names[1], dim_names[2]),
                            fill = cell_color, show.legend = show.legend, size =  point_size)

    } else {
      if(scale_alpha_with_expression == TRUE) {
        pl <- pl + ggplot2::geom_point(data = annotated_gene_DT, aes_string(x = dim_names[1], y = dim_names[2], fill = gene, alpha = gene),
                              show.legend = show.legend, shape = 21, size = point_size,
                              color = point.border.col, stroke = point.border.stroke)
      } else {
        pl <- pl + ggplot2::geom_point(data = annotated_gene_DT, aes_string(x = dim_names[1], y = dim_names[2], fill = gene),
                              show.legend = show.legend, shape = 21,
                              size =  point_size,
                              color = point.border.col, stroke = point.border.stroke)
      }

      pl <- pl + ggplot2::scale_fill_gradient2(low = 'grey', mid = 'lightgrey', high = 'red', midpoint = midpoint)
    }

    pl <- pl + ggplot2::labs(x = 'coord x', y = 'coord y')

    if(show_plots == TRUE) {
      print(pl)
    }

    savelist[[gene]] <- pl
  }

  # combine plots with cowplot
  combo_plot <- cowplot::plot_grid(plotlist = savelist,
                          ncol = cow_n_col,
                          rel_heights = cow_rel_h, rel_widths = cow_rel_w, align = cow_align)
  combined_cowplot = cowplot::plot_grid(combo_plot)

  return(combined_cowplot)

}



#' @title visSpatDimGenePlot
#' @name visSpatDimGenePlot
#' @description Visualize cells according to spatial AND dimension reduction coordinates
#' @param gobject giotto object
#' @param expression_values gene expression values to use
#' @param plot_alignment direction to align plot
#' @param genes genes to show
#' @param dim_reduction_to_use dimension reduction to use
#' @param dim_reduction_name dimension reduction name
#' @param dim1_to_use dimension to use on x-axis
#' @param dim2_to_use dimension to use on y-axis
#' @param dim_point_size dim reduction plot: point size
#' @param dim.point.border.col color of border around points
#' @param dim.point.border.stroke stroke size of border around points
#' @param show_NN_network show underlying NN network
#' @param nn_network_to_use type of NN network to use (kNN vs sNN)
#' @param network_name name of NN network to use, if show_NN_network = TRUE
#' @param edge_alpha_dim dim reduction plot: column to use for alpha of the edges
#' @param scale_alpha_with_expression scale expression with ggplot alpha parameter
#' @param spatial_network_name name of spatial network to use
#' @param spatial_grid_name name of spatial grid to use
#' @param spatial_point_size spatial plot: point size
#' @param spatial.point.border.col color of border around points
#' @param spatial.point.border.stroke stroke size of border around points
#' @param midpoint size of point (cell)
#' @param point_size size of point (cell)
#' @param cow_n_col cowplot param: how many columns
#' @param cow_rel_h cowplot param: relative height
#' @param cow_rel_w cowplot param: relative width
#' @param cow_align cowplot param: how to align
#' @param show.legend show legend
#' @param show.plot show plot
#' @return ggplot
#' @details Description of parameters.
#' @export
#' @examples
#'     visSpatDimGenePlot(gobject)
visSpatDimGenePlot <- function(gobject,
                               expression_values = c('normalized', 'scaled', 'custom'),
                               plot_alignment = c('horizontal', 'vertical'),
                               genes,
                               dim_reduction_to_use = 'umap',
                               dim_reduction_name = 'umap',
                               dim1_to_use = 1,
                               dim2_to_use = 2,
                               dim_point_size = 1,
                               dim.point.border.col = 'black',
                               dim.point.border.stroke = 0.1,
                               show_NN_network = F,
                               nn_network_to_use = 'sNN',
                               network_name = 'sNN.pca',
                               edge_alpha_dim = NULL,
                               scale_alpha_with_expression = TRUE,
                               spatial_network_name = 'spatial_network',
                               spatial_grid_name = 'spatial_grid',
                               spatial_point_size = 1,
                               spatial.point.border.col = 'black',
                               spatial.point.border.stroke = 0.1,
                               midpoint = 0,
                               point_size = 1,
                               cow_n_col = 2,
                               cow_rel_h = 1,
                               cow_rel_w = 1,
                               cow_align = 'h',
                               show.legend = T,
                               show_plots = F) {

  {

    plot_alignment = match.arg(plot_alignment, choices = c('vertical', 'horizontal'))

    # dimension reduction plot
    dmpl = visDimGenePlot(gobject = gobject, genes = genes, expression_values = expression_values,
                          dim_reduction_to_use = dim_reduction_to_use, dim_reduction_name = dim_reduction_name,
                          dim1_to_use = dim1_to_use, dim2_to_use = dim2_to_use,
                          show_NN_network = show_NN_network,
                          nn_network_to_use = nn_network_to_use, network_name = network_name,
                          edge_alpha = edge_alpha_dim, midpoint = midpoint,
                          show.legend = show.legend,scale_alpha_with_expression = scale_alpha_with_expression,
                          point_size = dim_point_size,
                          point.border.col = dim.point.border.col,
                          point.border.stroke = dim.point.border.stroke,
                          cow_n_col = cow_n_col, cow_rel_h = cow_rel_h, cow_rel_w = cow_rel_w, cow_align = cow_align, show_plots = show_plots)


    # spatial plot
    spl = visGenePlot(gobject = gobject, genes = genes, expression_values = expression_values,
                      spatial_network_name = spatial_network_name, midpoint = midpoint,
                      show_plots = show_plots, scale_alpha_with_expression = scale_alpha_with_expression,
                      point_size = spatial_point_size,
                      point.border.col = spatial.point.border.col,
                      point.border.stroke = spatial.point.border.stroke,
                      show.legend = show.legend,
                      cow_n_col = cow_n_col, cow_rel_h = cow_rel_h, cow_rel_w = cow_rel_w, cow_align = cow_align)

    if(plot_alignment == 'vertical') {
      combo_plot <- cowplot::plot_grid(dmpl, spl, ncol = 1, rel_heights = c(1), rel_widths = c(1), align = 'v')
      return(cowplot::plot_grid(combo_plot))
    } else {
      combo_plot <- cowplot::plot_grid(dmpl, spl, ncol = 2, rel_heights = c(1), rel_widths = c(1), align = 'h')
      return(cowplot::plot_grid(combo_plot))
    }
  }



}



