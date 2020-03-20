## original functions ####
## -------------------- ##

#' @title visPlot_3D_plotly
#' @name visPlot_3D_plotly
#' @description Visualize cells according to spatial coordinates
#' @param gobject giotto object
#' @param sdimx x-axis dimension name (default = 'sdimx')
#' @param sdimy y-axis dimension name (default = 'sdimy')
#' @param sdimz z-axis dimension name (default = 'sdimz')
#' @param point_size size of point (cell)
#' @param point_border_col color of border around points
#' @param point_border_stroke stroke size of border around points
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
#' @param show_legend show legend
#' @param show_plot show plot
#' @return ggplot
#' @details Description of parameters.
#' @export
#' @examples
#'     visPlot_3D_plotly(gobject)
visPlot_3D_plotly = function(gobject,
                             sdimx = NULL,
                             sdimy = NULL,
                             sdimz = NULL,
                             point_size = 3,
                             cell_color = NULL,
                             cell_color_code = NULL,
                             select_cell_groups = NULL,
                             select_cells = NULL,
                             show_other_cells = T,
                             other_cell_color = "lightgrey",
                             other_point_size = 0.5,
                             show_network = F,
                             network_color = NULL,
                             network_alpha = 1,
                             other_cell_alpha = 0.5,
                             spatial_network_name = 'spatial_network',
                             spatial_grid_name = 'spatial_grid',
                             title = '',
                             show_legend = T,
                             axis_scale = c("cube","real","custom"),
                             custom_ratio = NULL,
                             x_ticks = NULL,
                             y_ticks = NULL,
                             z_ticks = NULL,
                             show_plot = F) {
  
  
  ## get spatial cell locations
  cell_locations  = gobject@spatial_locs
  
  ## extract spatial network
  if(!is.null(spatial_network_name)) {
    #spatial_network = gobject@spatial_network[[spatial_network_name]]
    spatial_network = select_spatialNetwork(gobject, name = spatial_network_name, return_network_Obj = FALSE)
  } else {
    spatial_network = NULL
  }
  
  ## extract spatial grid
  if(!is.null(spatial_grid_name)) {
    spatial_grid    = gobject@spatial_grid[[spatial_grid_name]]
  } else {
    spatial_grid = NULL
  }
  
  ## get cell metadata
  cell_metadata   = gobject@cell_metadata
  
  
  if(nrow(cell_metadata) == 0) {
    cell_locations_metadata = cell_locations
  } else {
    cell_locations_metadata <- merge(cell_locations, cell_metadata,by = 'cell_ID')
  }
  
  ## create subsets if needed
  if(!is.null(select_cells) & !is.null(select_cell_groups)) {
    cat('You have selected both individual cell IDs and a group of cells \n')
    group_cell_IDs = cell_locations_metadata[get(cell_color) %in% select_cell_groups][['cell_ID']]
    select_cells = unique(c(select_cells, group_cell_IDs))
  } else if(!is.null(select_cell_groups)) {
    select_cells = cell_locations_metadata[get(cell_color) %in% select_cell_groups][['cell_ID']]
  }
  
  if(!is.null(select_cells)) {
    cell_locations_metadata_other = cell_locations_metadata[!cell_locations_metadata$cell_ID %in% select_cells]
    cell_locations_metadata_selected = cell_locations_metadata[cell_locations_metadata$cell_ID %in% select_cells]
    spatial_network <- spatial_network[spatial_network$to %in% select_cells & spatial_network$from %in% select_cells]
    
    # if specific cells are selected
    cell_locations_metadata = cell_locations_metadata_selected
  }
  
  
  
  ### set scale
  axis_scale = match.arg(axis_scale, c("cube","real","custom"))
  
  ### set ratio
  ratio = plotly_axis_scale_3D(cell_locations,
                               sdimx = sdimx, sdimy = sdimy, sdimz = sdimz,
                               mode = axis_scale, custom_ratio = custom_ratio)
  
  pl <- plotly::plot_ly()
  if(!is.null(cell_color)) {
    if(cell_color %in% colnames(cell_locations_metadata)){
      if(is.null(cell_color_code)) {
        number_colors=length(unique(cell_locations_metadata[[cell_color]]))
        cell_color_code = Giotto:::getDistinctColors(n = number_colors)
      }
      cell_locations_metadata[[cell_color]] <- as.factor(cell_locations_metadata[[cell_color]])
      pl <- pl %>% plotly::add_trace(type = 'scatter3d',mode = "markers",data = cell_locations_metadata,
                                     x = ~sdimx, y = ~sdimy, z = ~sdimz,
                                     color = cell_locations_metadata[[cell_color]],
                                     colors = cell_color_code,
                                     marker = list(size = point_size))
      
      
      if(!is.null(select_cells) & show_other_cells){
        pl <- pl %>% plotly::add_trace(type = "scatter3d",mode="markers",
                                       data=cell_locations_metadata_other,
                                       name = "unselected cells",
                                       x=~sdimx,y=~sdimy,z=~sdimz,
                                       marker = list(size = other_point_size,color = other_cell_color),
                                       opacity=other_cell_alpha)
      }
    }
    else{
      cat('cell_color does not exist! \n')
    }
  } else {
    pl <- pl %>% plotly::add_trace(type = 'scatter3d',
                                   data = cell_locations_metadata,
                                   x = ~sdimx,
                                   y = ~sdimy,
                                   z = ~sdimz,
                                   mode = 'markers',
                                   marker = list(size = point_size),
                                   colors = 'lightblue',name = "selected cells")
    if(!is.null(select_cells) & show_other_cells){
      pl <- pl %>% plotly::add_trace(type = "scatter3d",
                                     mode="markers",
                                     data=cell_locations_metadata_other,
                                     name = "unselected cells",
                                     x=~sdimx,y=~sdimy,z=~sdimz,
                                     marker = list(size = other_point_size,color=other_cell_color),
                                     opacity = other_cell_alpha)
    }
  }
  
  ## plot spatial network
  if(!is.null(spatial_network) & show_network == TRUE) {
    if(is.null(network_color)) {
      network_color = 'red'
    }
    edges <- plotly_network(spatial_network)
    
    pl <- pl %>% plotly::add_trace(name = "sptial network",
                                   mode = "lines",
                                   type = "scatter3d",
                                   data = edges,
                                   x = ~x,y=~y,z=~z,
                                   line=list(color=network_color,width = 0.5),opacity=network_alpha)
  }
  ## plot spatial grid
  # 3D grid is not clear to view
  
  
  pl <- pl %>%
    plotly::layout(scene = list(xaxis = list(title = 'X',nticks = x_ticks),
                                yaxis = list(title = 'Y',nticks = y_ticks),
                                zaxis = list(title = 'Z',nticks = z_ticks),
                                aspectmode='manual',
                                aspectratio = list(x=ratio[[1]], y=ratio[[2]], z=ratio[[3]])),
                   legend = list(x = 100, y = 0.5,
                                 font = list(family = "sans-serif",size = 12)))
  
  
  return((pl))
}


#' @title visPlot_2D_ggplot
#' @name visPlot_2D_ggplot
#' @description Visualize cells according to spatial coordinates
#' @param gobject giotto object
#' @param sdimx x-axis dimension name (default = 'sdimx')
#' @param sdimy y-axis dimension name (default = 'sdimy')
#' @param point_size size of point (cell)
#' @param point_border_col color of border around points
#' @param point_border_stroke stroke size of border around points
#' @param cell_color color for cells (see details)
#' @param cell_color_code named vector with colors
#' @param color_as_factor convert color column to factor
#' @param select_cell_groups select subset of cells/clusters based on cell_color parameter
#' @param select_cells select subset of cells based on cell IDs
#' @param show_other_cells display not selected cells
#' @param other_cell_color color of not selected cells
#' @param show_network show underlying spatial network
#' @param network_color color of spatial network
#' @param spatial_network_name name of spatial network to use
#' @param show_grid show spatial grid
#' @param grid_color color of spatial grid
#' @param spatial_grid_name name of spatial grid to use
#' @param coord_fix_ratio fix ratio between x and y-axis
#' @param title title of plot
#' @param show_legend show legend
#' @param show_plot show plot
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_dir directory to save the plot
#' @param save_folder (optional) folder in directory to save the plot
#' @param save_name name of plot
#' @param save_format format of plot (e.g. tiff, png, pdf, ...)
#' @param show_saved_plot load & display the saved plot
#' @return ggplot
#' @details Description of parameters.
#' @export
#' @examples
#'     visPlot_2D_ggplot(gobject)
visPlot_2D_ggplot = function(gobject,
                             sdimx = NULL,
                             sdimy = NULL,
                             point_size = 3,
                             point_border_col = 'black',
                             point_border_stroke = 0.1,
                             cell_color = NULL,
                             cell_color_code = NULL,
                             color_as_factor = T,
                             select_cell_groups = NULL,
                             select_cells = NULL,
                             show_other_cells = T,
                             other_cell_color = 'lightgrey',
                             show_network = F,
                             network_color = NULL,
                             network_alpha = 1,
                             other_cells_alpha = 0.1,
                             spatial_network_name = 'spatial_network',
                             show_grid = F,
                             grid_color = NULL,
                             spatial_grid_name = 'spatial_grid',
                             coord_fix_ratio = 0.6,
                             title = '',
                             show_legend = T,
                             axis_scale = c("cube","real","custom"),
                             custom_ratio = NULL,
                             x_ticks = NULL,
                             y_ticks = NULL,
                             z_ticks = NULL,
                             show_plot = F,
                             return_plot = TRUE,
                             save_plot = F,
                             save_dir = NULL,
                             save_folder = NULL,
                             save_name = NULL,
                             save_format = NULL,
                             show_saved_plot = F,
                             ...) {
  
  
  ## get spatial cell locations
  cell_locations  = gobject@spatial_locs
  
  ## extract spatial network
  if(!is.null(spatial_network_name)) {
    spatial_network = select_spatialNetwork(gobject, name = spatial_network_name, return_network_Obj = FALSE)
  } else {
    spatial_network = NULL
  }
  
  ## extract spatial grid
  if(!is.null(spatial_grid_name)) {
    spatial_grid    = gobject@spatial_grid[[spatial_grid_name]]
  } else {
    spatial_grid = NULL
  }
  
  ## get cell metadata
  cell_metadata   = gobject@cell_metadata
  cell_metadata   = cell_metadata[, !grepl('cell_ID', colnames(cell_metadata)), with = F]
  
  if(nrow(cell_metadata) == 0) {
    cell_locations_metadata = cell_locations
  } else {
    cell_locations_metadata <- cbind(cell_locations, cell_metadata)
  }
  
  ## create subsets if needed
  if(!is.null(select_cells) & !is.null(select_cell_groups)) {
    cat('You have selected both individual cell IDs and a group of cells \n')
    group_cell_IDs = cell_locations_metadata[get(cell_color) %in% select_cell_groups][['cell_ID']]
    select_cells = unique(c(select_cells, group_cell_IDs))
  } else if(!is.null(select_cell_groups)) {
    select_cells = cell_locations_metadata[get(cell_color) %in% select_cell_groups][['cell_ID']]
  }
  
  if(!is.null(select_cells)) {
    cell_locations_metadata_other = cell_locations_metadata[!cell_locations_metadata$cell_ID %in% select_cells]
    cell_locations_metadata_selected = cell_locations_metadata[cell_locations_metadata$cell_ID %in% select_cells]
    spatial_network <- spatial_network[spatial_network$to %in% select_cells & spatial_network$from %in% select_cells]
    
    # if specific cells are selected
    cell_locations_metadata = cell_locations_metadata_selected
  }
  
  
  
  ### create 2D plot with ggplot ###
  cat('create 2D plot with ggplot \n')
  
  
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
    pl <- pl + ggplot2::geom_point(data = cell_locations_metadata,
                                   aes_string(x = sdimx, y = sdimy),
                                   show.legend = show_legend, shape = 21,
                                   fill = cell_color, size = point_size,
                                   stroke = point_border_stroke, color = point_border_col)
    
  }
  
  else if (is.character(cell_color)) {
    
    if(cell_color %in% colnames(cell_locations_metadata)) {
      
      if(color_as_factor == TRUE) {
        if(is.null(select_cells)) {
          factor_data = factor(cell_locations_metadata[[cell_color]])
          cell_locations_metadata[[cell_color]] <- factor_data
        } else {
          factor_data_selected = factor(cell_locations_metadata_selected[[cell_color]])
          cell_locations_metadata_selected[[cell_color]] <- factor_data_selected
          factor_data_other = factor(cell_locations_metadata_other[[cell_color]])
          cell_locations_metadata_other[[cell_color]] <- factor_data_other
        }
        
      }
      
      if(is.null(select_cells)) {
        pl <- pl + ggplot2::geom_point(data = cell_locations_metadata,
                                       aes_string(x = sdimx, y = sdimy, fill = cell_color),
                                       show.legend = show_legend, shape = 21, size = point_size,
                                       stroke = point_border_stroke, color = point_border_col)
      } else {
        if(show_other_cells){
          cell_color_other = other_cell_color
          pl <- pl + ggplot2::geom_point(data = cell_locations_metadata_other,
                                         aes_string(x = sdimx, y = sdimy), fill = cell_color_other,
                                         show.legend = show_legend, shape = 21, size = point_size/2,
                                         stroke = point_border_stroke, color = point_border_col)
        }
        
        pl <- pl + ggplot2::geom_point(data = cell_locations_metadata_selected,
                                       aes_string(x = sdimx, y = sdimy, fill = cell_color),
                                       show.legend = show_legend, shape = 21, size = point_size,
                                       stroke = point_border_stroke, color = point_border_col)
      }
      
      
      
      
      if(!is.null(cell_color_code)) {
        pl <- pl + ggplot2::scale_fill_manual(values = cell_color_code)
      } else if(color_as_factor == T) {
        if(is.null(select_cells)) {
          number_colors = length(unique(factor_data))
          cell_color_code = Giotto:::getDistinctColors(n = number_colors)
          names(cell_color_code) = unique(factor_data)
        } else {
          number_colors = length(unique(factor_data_selected))
          cell_color_code = Giotto:::getDistinctColors(n = number_colors)
          names(cell_color_code) = unique(factor_data_selected)
        }
        pl <- pl + ggplot2::scale_fill_manual(values = cell_color_code)
      } else if(color_as_factor == F){
        pl <- pl + ggplot2::scale_fill_gradient(low = 'blue', high = 'red')
      }
      
    } else {
      
      
      if(is.null(select_cells)) {
        pl <- pl + ggplot2::geom_point(data = cell_locations_metadata, aes_string(x = sdimx, y = sdimy),
                                       show.legend = show_legend, shape = 21, fill = cell_color,
                                       size = point_size, stroke = point_border_stroke, color = point_border_col)
      } else {
        if(show_other_cells){
          cell_color_other = 'grey'
          pl <- pl + ggplot2::geom_point(data = cell_locations_metadata_other, aes_string(x = sdimx, y = sdimy),
                                         show.legend = show_legend, shape = 21, fill = cell_color_other,
                                         size = point_size/2,
                                         stroke = point_border_stroke, color = point_border_col)
          
          pl <- pl + ggplot2::geom_point(data = cell_locations_metadata_selected, aes_string(x = sdimx, y = sdimy),
                                         show.legend = show_legend, shape = 21, fill = cell_color,
                                         size = point_size,
                                         stroke = point_border_stroke, color = point_border_col)
        }
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
  
  
  ## print plot
  if(show_plot == TRUE) {
    print(pl)
  }
  
  
  ## save plot
  if(save_plot == TRUE) {
    
    ggplot_save_function(gobject = gobject,
                         plot_object = pl,
                         save_dir = save_dir,
                         save_folder = save_folder,
                         save_name = save_name,
                         save_format = save_format,
                         show_saved_plot = show_saved_plot,
                         ...)
  }
  
  ## return plot
  if(return_plot == TRUE) return(pl)
  
}



#' @title visPlot_2D_plotly
#' @name visPlot_2D_plotly
#' @description Visualize cells according to spatial coordinates
#' @param gobject giotto object
#' @param sdimx x-axis dimension name (default = 'sdimx')
#' @param sdimy y-axis dimension name (default = 'sdimy')
#' @param point_size size of point (cell)
#' @param cell_color color for cells (see details)
#' @param cell_color_code named vector with colors
#' @param color_as_factor convert color column to factor
#' @param select_cell_groups select a subset of the groups from cell_color
#' @param show_network show underlying spatial network
#' @param network_color color of spatial network
#' @param spatial_network_name name of spatial network to use
#' @param show_grid show spatial grid
#' @param grid_color color of spatial grid
#' @param grid_alpha alpha of spatial grid
#' @param spatial_grid_name name of spatial grid to use
#' @param show_legend show legend
#' @param show_plot show plot
#' @return plotly
#' @details Description of parameters.
#' @export
#' @examples
#'     visPlot_2D_plotly(gobject)
visPlot_2D_plotly = function(gobject,
                             sdimx = NULL,
                             sdimy = NULL,
                             point_size = 3,
                             cell_color = NULL,
                             cell_color_code = NULL,
                             color_as_factor = T,
                             select_cell_groups = NULL,
                             select_cells = NULL,
                             show_other_cells = T,
                             other_cell_color = "lightgrey",
                             other_point_size = 0.5,
                             show_network = F,
                             network_color = "lightgray",
                             network_alpha = 1,
                             other_cell_alpha = 0.5,
                             spatial_network_name = 'spatial_network',
                             show_grid = F,
                             grid_color = NULL,
                             grid_alpha = 1,
                             spatial_grid_name = 'spatial_grid',
                             show_legend = T,
                             axis_scale = c("cube","real","custom"),
                             custom_ratio = NULL,
                             x_ticks = NULL,
                             y_ticks = NULL,
                             show_plot = F) {
  
  
  ## get spatial cell locations
  cell_locations  = gobject@spatial_locs
  
  ##if(is.null(sdimx)|is.null(sdimy)){
  ##   warning("sdimx or sdimy is null! Default value will be used!\n")
  ##   sdimx = "sdimx"
  ##   sdimy = "sdimy"
  ##}
  
  ## extract spatial network
  if(!is.null(spatial_network_name)) {
    spatial_network = select_spatialNetwork(gobject, name = spatial_network_name, return_network_Obj = FALSE)
  } else {
    spatial_network = NULL
  }
  
  ## extract spatial grid
  if(!is.null(spatial_grid_name)) {
    spatial_grid    = gobject@spatial_grid[[spatial_grid_name]]
  } else {
    spatial_grid = NULL
  }
  
  ## get cell metadata
  cell_metadata   = gobject@cell_metadata
  cell_metadata   = cell_metadata[, !grepl('cell_ID', colnames(cell_metadata)), with = F]
  
  if(nrow(cell_metadata) == 0) {
    cell_locations_metadata = cell_locations
  } else {
    cell_locations_metadata <- cbind(cell_locations, cell_metadata)
  }
  
  ## create subsets if needed
  if(!is.null(select_cells) & !is.null(select_cell_groups)) {
    cat('You have selected both individual cell IDs and a group of cells \n')
    group_cell_IDs = cell_locations_metadata[get(cell_color) %in% select_cell_groups][['cell_ID']]
    select_cells = unique(c(select_cells, group_cell_IDs))
  }
  else if(!is.null(select_cell_groups)) {
    select_cells = cell_locations_metadata[get(cell_color) %in% select_cell_groups][['cell_ID']]
  }
  
  if(!is.null(select_cells)) {
    cell_locations_metadata_other = cell_locations_metadata[!cell_locations_metadata$cell_ID %in% select_cells]
    cell_locations_metadata_selected = cell_locations_metadata[cell_locations_metadata$cell_ID %in% select_cells]
    spatial_network <- spatial_network[spatial_network$to %in% select_cells & spatial_network$from %in% select_cells]
    
    # if specific cells are selected
    cell_locations_metadata = cell_locations_metadata_selected
  }
  
  
  
  ### set scale
  axis_scale = match.arg(axis_scale, c("cube","real","custom"))
  
  ### set ratio
  ratio = plotly_axis_scale_2D(cell_locations,
                               sdimx = sdimx,
                               sdimy = sdimy,
                               mode = axis_scale,
                               custom_ratio = custom_ratio)
  pl <- plotly::plot_ly()
  
  if(show_network == TRUE) {
    if(is.null(spatial_network)){
      stop("No usable spatial network specified! Please choose a network with spatial_network_name=xxx")
    }
    else{
      if(is.null(network_alpha)) {
        network_alpha = 0.5
      }
      else if(is.character(network_alpha)){
        warning("Edge_alpha for plotly mode is not adjustable yet. Default 0.5 will be set\n")
        network_alpha = 0.5
      }
      pl <- pl %>% plotly::add_segments(name = spatial_network_name,
                                        type = "scatter",
                                        x = spatial_network[["sdimx_begin"]],
                                        y = spatial_network[["sdimy_begin"]],
                                        xend = spatial_network[["sdimx_end"]],
                                        yend = spatial_network[["sdimy_end"]],
                                        line = list(color = network_color,
                                                    width = 0.5),
                                        opacity=network_alpha)
    }
  }
  if(show_grid == TRUE){
    if(is.null(spatial_grid)){
      stop("No usable spatial grid specified! Please choose a network with spatial_grid_name=xxx")
    }
    else{
      if(is.null(grid_color)) {
        grid_color = 'black'
      }
      edges <- plotly_grid(spatial_grid)
      pl <- pl %>% plotly::add_segments(name = "spatial_grid",
                                        type = "scatter",
                                        data = edges,
                                        x = ~x,
                                        y = ~y,
                                        xend = ~x_end,
                                        yend = ~y_end,
                                        line = list(color = grid_color,
                                                    width = 1),
                                        opacity=grid_alpha)
      
    }
  }
  
  if(!is.null(cell_color)) {
    if(cell_color %in% colnames(cell_locations_metadata)){
      if(is.null(cell_color_code)) {
        number_colors=length(unique(cell_locations_metadata[[cell_color]]))
        cell_color_code = Giotto:::getDistinctColors(n = number_colors)
      }
      cell_locations_metadata[[cell_color]] <- as.factor(cell_locations_metadata[[cell_color]])
      pl <- pl %>% plotly::add_trace(type = 'scatter',
                                     mode = 'markers',
                                     x = cell_locations_metadata[[sdimx]],
                                     y = cell_locations_metadata[[sdimy]],
                                     color = cell_locations_metadata[[cell_color]],
                                     colors = cell_color_code,
                                     marker = list(size = point_size))
      
      
      if(!is.null(select_cells) & show_other_cells){
        pl <- pl %>% plotly::add_trace(type = "scatter",
                                       mode="markers",
                                       data=cell_locations_metadata_other,
                                       name = "unselected cells",
                                       x=~sdimx,
                                       y=~sdimy,
                                       marker = list(size = other_point_size,color = other_cell_color),
                                       opacity=other_cell_alpha)
      }
    }
    else{
      cat('cell_color does not exist! \n')
    }
  } else {
    pl <- pl %>% plotly::add_trace(type = 'scatter',
                                   mode = 'markers',
                                   name = "selected cells",
                                   x = cell_locations_metadata[[sdimx]],
                                   y = cell_locations_metadata[[sdimy]],
                                   colors = 'lightblue',
                                   marker = list(size = point_size))
    if(!is.null(select_cells) & show_other_cells){
      pl <- pl %>% plotly::add_trace(type = "scatter",
                                     mode="markers",
                                     data=cell_locations_metadata_other,
                                     name = "unselected cells",
                                     x=~sdimx,y=~sdimy,
                                     marker = list(size = other_point_size,color = other_cell_color),
                                     opacity = other_cell_alpha)
    }
  }
  
  
  pl <- pl %>%
    plotly::layout(list(xaxis = list(title = 'X',nticks = x_ticks),
                        yaxis = list(title = 'Y',nticks = y_ticks)),
                   legend = list(x = 100, y = 0.5,
                                 font = list(family = "sans-serif",size = 12)))
  
  
  return((pl))
  
  
}



#' @title visPlot
#' @name visPlot
#' @description Visualize cells according to spatial coordinates
#' @param gobject giotto object
#' @param sdimx x-axis dimension name (default = 'sdimx')
#' @param sdimy y-axis dimension name (default = 'sdimy')
#' @param sdimz z-axis dimension name (default = 'sdimz')
#' @param point_size size of point (cell)
#' @param point_border_col color of border around points
#' @param point_border_stroke stroke size of border around points
#' @param cell_color color for cells (see details)
#' @param cell_color_code named vector with colors
#' @param color_as_factor convert color column to factor
#' @param select_cell_groups select subset of cells/clusters based on cell_color parameter
#' @param select_cells select subset of cells based on cell IDs
#' @param show_other_cells display not selected cells
#' @param other_cell_color color of not selected cells
#' @param show_network show underlying spatial network
#' @param network_color color of spatial network
#' @param spatial_network_name name of spatial network to use
#' @param show_grid show spatial grid
#' @param grid_color color of spatial grid
#' @param spatial_grid_name name of spatial grid to use
#' @param coord_fix_ratio fix ratio between x and y-axis
#' @param title title of plot
#' @param show_legend show legend
#' @param show_plot show plot
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_dir directory to save the plot
#' @param save_folder (optional) folder in directory to save the plot
#' @param save_name name of plot
#' @param save_format format of plot (e.g. tiff, png, pdf, ...)
#' @param show_saved_plot load & display the saved plot
#' @return ggplot
#' @details Description of parameters.
#' @export
#' @examples
#'     visPlot(gobject)
visPlot <- function(gobject,
                    sdimx = NULL,
                    sdimy = NULL,
                    sdimz = NULL,
                    point_size = 3,
                    point_border_col = 'black',
                    point_border_stroke = 0.1,
                    cell_color = NULL,
                    cell_color_code = NULL,
                    color_as_factor = T,
                    select_cell_groups = NULL,
                    select_cells = NULL,
                    show_other_cells = T,
                    other_cell_color = 'lightgrey',
                    show_network = F,
                    network_color = NULL,
                    network_alpha = 1,
                    other_cell_alpha = 0.1,
                    spatial_network_name = 'spatial_network',
                    show_grid = F,
                    grid_color = NULL,
                    grid_alpha = 1,
                    spatial_grid_name = 'spatial_grid',
                    coord_fix_ratio = 0.6,
                    title = '',
                    show_legend = T,
                    axis_scale = c("cube","real","custom"),
                    custom_ratio = NULL,
                    x_ticks = NULL,
                    y_ticks = NULL,
                    z_ticks = NULL,
                    plot_method = c('ggplot', 'plotly'),
                    show_plot = F,
                    return_plot = TRUE,
                    save_plot = F,
                    save_dir = NULL,
                    save_folder = NULL,
                    save_name = NULL,
                    save_format = NULL,
                    show_saved_plot = F,
                    ...) {
  
  
  ## decide plot method
  plot_method = match.arg(plot_method, choices = c('ggplot', 'plotly'))
  
  
  if(plot_method == 'ggplot') {
    
    if(is.null(sdimx) | is.null(sdimy)) {
      
      warning("plot_method = ggplot, but spatial dimensions for sdimx and sdimy for 2D plotting are not given. \n
              It will default to the 'sdimx' and 'sdimy' ")
      sdimx = 'sdimx'
      sdimy = 'sdimy'
      #stop('\n ggplot is in 2D and you need to define sdimx and sdimy \n')
    }
    
    result = visPlot_2D_ggplot(gobject = gobject,
                               sdimx = sdimx,
                               sdimy = sdimy,
                               point_size = point_size,
                               point_border_col = point_border_col,
                               point_border_stroke = point_border_stroke,
                               cell_color = cell_color,
                               cell_color_code = cell_color_code,
                               color_as_factor = color_as_factor,
                               select_cell_groups = select_cell_groups,
                               select_cells = select_cells,
                               show_network = show_network,
                               other_cell_color = other_cell_color,
                               network_color = network_color,
                               network_alpha = network_alpha,
                               other_cells_alpha = other_cell_alpha,
                               show_other_cells = show_other_cells,
                               spatial_network_name = spatial_network_name,
                               show_grid = show_grid,
                               grid_color = grid_color,
                               spatial_grid_name = spatial_grid_name,
                               coord_fix_ratio = coord_fix_ratio,
                               title = title,
                               show_legend = show_legend,
                               axis_scale = axis_scale,
                               custom_ratio = custom_ratio,
                               x_ticks = x_ticks,
                               y_ticks = y_ticks,
                               z_ticks = z_ticks,
                               show_plot = show_plot,
                               return_plot = return_plot,
                               save_plot = save_plot,
                               save_dir = save_dir,
                               save_folder = save_folder,
                               save_name = save_name,
                               save_format = save_format,
                               show_saved_plot = show_saved_plot,
                               ...)
    
    
  }
  else if(plot_method == 'plotly') {
    
    if(length(c(sdimx, sdimy, sdimz)) == 3) {
      
      result = visPlot_3D_plotly(gobject = gobject,
                                 sdimx = sdimx,
                                 sdimy = sdimy,
                                 sdimz = sdimz,
                                 point_size = point_size,
                                 cell_color = cell_color,
                                 cell_color_code = cell_color_code,
                                 select_cell_groups = select_cell_groups,
                                 select_cells = select_cells,
                                 show_other_cells = show_other_cells,
                                 show_network = show_network,
                                 other_cell_alpha = other_cell_alpha,
                                 network_color = network_color,
                                 network_alpha = network_alpha,
                                 spatial_network_name = spatial_network_name,
                                 spatial_grid_name = spatial_grid_name,
                                 title = title,
                                 show_legend = show_legend,
                                 axis_scale = axis_scale,
                                 custom_ratio = custom_ratio,
                                 x_ticks = x_ticks,
                                 y_ticks = y_ticks,
                                 z_ticks = z_ticks,
                                 show_plot = show_plot)
      
    }
    else {
      
      if(is.null(sdimx) | is.null(sdimy)) {
        
        warning("plot_method = plotly, but spatial dimensions for sdimx and sdimy for 2D plotting are not given. \n
                It will default to the 'sdimx' and 'sdimy'. For 3D visualization also specify sdimz.")
        sdimx = 'sdimx'
        sdimy = 'sdimy'
        #stop('\n plotly in 2D requires you to define sdimx and sdimy \n')
        
      }
      
      ## run: visPlot_2D_plotly
      result = visPlot_2D_plotly(gobject = gobject,
                                 sdimx = sdimx,
                                 sdimy = sdimy,
                                 point_size = point_size,
                                 cell_color = cell_color,
                                 cell_color_code = cell_color_code,
                                 color_as_factor = color_as_factor,
                                 select_cell_groups = select_cell_groups,
                                 select_cells = select_cells,
                                 show_other_cells = show_other_cells,
                                 show_network = show_network,
                                 network_color =network_color,
                                 network_alpha = network_alpha,
                                 other_cell_alpha = other_cell_alpha,
                                 spatial_network_name = spatial_network_name,
                                 show_grid = show_grid,
                                 grid_color = grid_color,
                                 grid_alpha = grid_alpha,
                                 spatial_grid_name = spatial_grid_name,
                                 show_legend = show_legend,
                                 axis_scale = axis_scale,
                                 custom_ratio = custom_ratio,
                                 x_ticks = x_ticks,
                                 y_ticks = y_ticks,
                                 show_plot = show_plot)
      
      
    }
    
  }
  return(result)
  
}




#' @title visGenePlot_2D_ggplot
#' @name visGenePlot_2D_ggplot
#' @description Visualize cells and gene expression according to spatial coordinates
#' @param gobject giotto object
#' @param expression_values gene expression values to use
#' @param genes genes to show
#' @param genes_high_color color represents high gene expression
#' @param genes_mid_color color represents middle gene expression
#' @param genes_low_color color represents low gene expression
#' @param show_network show underlying spatial network
#' @param network_color color of spatial network
#' @param spatial_network_name name of spatial network to use
#' @param show_grid show spatial grid
#' @param grid_color color of spatial grid
#' @param spatial_grid_name name of spatial grid to use
#' @param midpoint expression midpoint
#' @param scale_alpha_with_expression scale expression with ggplot alpha parameter
#' @param point_size size of point (cell)
#' @param point_border_col color of border around points
#' @param point_border_stroke stroke size of border around points
#' @param cow_n_col cowplot param: how many columns
#' @param cow_rel_h cowplot param: relative height
#' @param cow_rel_w cowplot param: relative width
#' @param cow_align cowplot param: how to align
#' @param show_legend show legend
#' @param show_plots show plots
#' @return ggplot
#' @details Description of parameters.
#' @export
#' @examples
#'     visGenePlot_2D_ggplot(gobject)
visGenePlot_2D_ggplot <- function(gobject,
                                  expression_values = c('normalized', 'scaled', 'custom'),
                                  genes,
                                  genes_high_color = 'darkred',
                                  genes_mid_color = "white",
                                  genes_low_color = "darkblue",
                                  show_network = F,
                                  network_color = NULL,
                                  spatial_network_name = 'spatial_network',
                                  edge_alpha = NULL,
                                  show_grid = F,
                                  grid_color = NULL,
                                  spatial_grid_name = 'spatial_grid',
                                  midpoint = 0,
                                  scale_alpha_with_expression = FALSE,
                                  point_size = 1,
                                  point_border_col = 'black',
                                  point_border_stroke = 0.1,
                                  show_legend = T,
                                  cow_n_col = 2,
                                  cow_rel_h = 1,
                                  cow_rel_w = 1,
                                  cow_align = 'h',
                                  show_plots = F){
  selected_genes = genes
  
  values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))
  expr_values = Giotto:::select_expression_values(gobject = gobject, values = values)
  
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
    spatial_network = select_spatialNetwork(gobject, name = spatial_network_name, return_network_Obj = FALSE)
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
      if(is.null(network_color)) {
        network_color = 'red'
      }
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
                                     color = point_border_col, size = point_size, stroke = point_border_stroke,
                                     show.legend = show_legend)
    } else {
      pl <- pl + ggplot2::geom_point(data = cell_locations_metadata_genes, aes_string(x = 'sdimx', y = 'sdimy',
                                                                                      fill = gene),
                                     shape = 21,
                                     color = point_border_col, size = point_size, stroke = point_border_stroke,
                                     show.legend = show_legend)
    }
    pl <- pl + ggplot2::scale_alpha_continuous(guide = 'none')
    pl <- pl + ggplot2::scale_fill_gradient2(low = genes_low_color, mid = genes_mid_color, high = genes_high_color,
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



#' @title visGenePlot_3D_plotly
#' @name visGenePlot_3D_plotly
#' @description Visualize cells and gene expression according to spatial coordinates
#' @param gobject giotto object
#' @param expression_values gene expression values to use
#' @param genes genes to show
#' @param genes_high_color color represents high gene expression
#' @param genes_mid_color color represents middle gene expression
#' @param genes_low_color color represents low gene expression
#' @param show_network show underlying spatial network
#' @param network_color color of spatial network
#' @param spatial_network_name name of spatial network to use
#' @param show_grid show spatial grid
#' @param grid_color color of spatial grid
#' @param spatial_grid_name name of spatial grid to use
#' @param point_size size of point (cell)
#' @param cow_n_col cowplot param: how many columns
#' @param cow_rel_h cowplot param: relative height
#' @param cow_rel_w cowplot param: relative width
#' @param cow_align cowplot param: how to align
#' @param show_legend show legend
#' @param axis_scale three mode to adjust axis scale
#' @param x_ticks number of ticks on x axis
#' @param y_ticks number of ticks on y axis
#' @param z_ticks number of ticks on z axis
#' @param show_plots show plots
#' @return plotly
#' @details Description of parameters.
#' @export
#' @examples
#'     visGenePlot_3D_plotly(gobject)
visGenePlot_3D_plotly <- function(gobject,
                                  expression_values = c('normalized', 'scaled', 'custom'),
                                  genes,
                                  show_network = F,
                                  network_color = NULL,
                                  spatial_network_name = 'spatial_network',
                                  edge_alpha = NULL,
                                  show_grid = F,
                                  genes_high_color = NULL,
                                  genes_mid_color = "white",
                                  genes_low_color = "blue",
                                  spatial_grid_name = 'spatial_grid',
                                  point_size = 1,
                                  show_legend = T,
                                  axis_scale = c("cube","real","custom"),
                                  custom_ratio = NULL,
                                  x_ticks = NULL,
                                  y_ticks = NULL,
                                  z_ticks = NULL,
                                  show_plots = F){
  selected_genes = genes
  
  values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))
  expr_values = Giotto:::select_expression_values(gobject = gobject, values = values)
  
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
    spatial_network = select_spatialNetwork(gobject, name = spatial_network_name, return_network_Obj = FALSE)
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
  
  
  
  axis_scale = match.arg(axis_scale, c("cube","real","custom"))
  
  ratio = plotly_axis_scale_3D(cell_locations_metadata_genes,sdimx = "sdimx",sdimy = "sdimy",sdimz = "sdimz",
                               mode = axis_scale,custom_ratio = custom_ratio)
  
  
  ## spatial network data
  if(!is.null(spatial_network) & show_network == TRUE){
    edges <- plotly_network(spatial_network)
  }
  ##Point layer
  if(length(selected_genes) > 4){
    stop("\n The max number of genes showed together is 4.Otherwise it will be too small to see\n
         \n If you have more genes to show, please divide them into groups\n")
  }
  savelist <- list()
  for(i in 1:length(selected_genes)){
    gene = selected_genes[i]
    if(!is.null(genes_high_color)){
      if(length(genes_high_color)!=length(selected_genes) & length(genes_high_color)!=1){
        stop('\n The number of genes and their corresbonding do not match\n')
      }
      else if(length(genes_high_color) == 1){
        genes_high_color = rep(genes_high_color,length(selected_genes))
      }
    }
    else{
      genes_high_color = rep("red",length(selected_genes))
    }
    pl <- plotly::plot_ly(name = gene,
                          
                          scene=paste("scene",i,sep = "")) %>%
      plotly::add_trace(data = cell_locations_metadata_genes,
                        type = 'scatter3d',mode = "markers",
                        x = ~sdimx, y = ~sdimy, z = ~sdimz,
                        marker = list(size = point_size),
                        color = cell_locations_metadata_genes[[gene]],
                        colors = c(genes_low_color,genes_mid_color,genes_high_color[i]))
    
    ## plot spatial network
    if(show_network == TRUE) {
      if(is.null(network_color)) {
        network_color = 'lightblue'
      }
      if(is.null(edge_alpha)) {
        edge_alpha = 0.5
      }
      else if (is.character(edge_alpha)){
        edge_alpha = 0.5
        cat("\nEdge_alpha for plotly mode is not adjustable yet. Default 0.5 will be set\n")
      }
      pl <- pl %>% plotly::add_trace(name = "sptial network",
                                     mode = "lines",
                                     type = "scatter3d",
                                     data = edges,
                                     x = ~x,y=~y,z=~z,
                                     line=list(color=network_color,width = 0.5),
                                     opacity = edge_alpha,
                                     showlegend = F)
    }
    ##plot spatial grid
    if(!is.null(spatial_grid) & show_grid == TRUE){
      cat("\n spatial grid is not clear in 3D plot \n")
    }
    pl <- pl %>% plotly::colorbar(title = gene)
    savelist[[gene]] <- pl
  }
  
  if(length(savelist) == 1){
    savelist[[1]] <- savelist[[1]] %>% plotly::layout(scene = list(xaxis = list(title = "X",nticks = x_ticks),
                                                                   yaxis = list(title = "Y",nticks = y_ticks),
                                                                   zaxis = list(title = "Z",nticks = z_ticks),
                                                                   aspectmode='manual',
                                                                   aspectratio = list(x=ratio[[1]],
                                                                                      y=ratio[[2]],
                                                                                      z=ratio[[3]])))
    if(show_plots){
      print(savelist[[1]])
    }
    return (savelist[[1]])
  }
  else if(length(savelist)==2){
    cowplot <- suppressWarnings(plotly::subplot(savelist)%>% plotly::layout(scene = list(xaxis = list(title = "X",nticks = x_ticks),
                                                                                         yaxis = list(title = "Y",nticks = y_ticks),
                                                                                         zaxis = list(title = "Z",nticks = z_ticks),
                                                                                         aspectmode='manual',
                                                                                         aspectratio = list(x=ratio[[1]],
                                                                                                            y=ratio[[2]],
                                                                                                            z=ratio[[3]])),
                                                                            scene2 = list(xaxis = list(title = "X",nticks = x_ticks),
                                                                                          yaxis = list(title = "Y",nticks = y_ticks),
                                                                                          zaxis = list(title = "Z",nticks = z_ticks),
                                                                                          aspectmode='manual',
                                                                                          aspectratio = list(x=ratio[[1]],
                                                                                                             y=ratio[[2]],
                                                                                                             z=ratio[[3]])),
                                                                            #annotations = annotations,
                                                                            legend = list(x = 100, y = 0)))
  }
  else if(length(savelist)==3){
    cowplot <- suppressWarnings(plotly::subplot(savelist)%>% plotly::layout(scene = list(xaxis = list(title = "X",nticks = x_ticks),
                                                                                         yaxis = list(title = "Y",nticks = y_ticks),
                                                                                         zaxis = list(title = "Z",nticks = z_ticks),
                                                                                         aspectmode='manual',
                                                                                         aspectratio = list(x=ratio[[1]],
                                                                                                            y=ratio[[2]],
                                                                                                            z=ratio[[3]])),
                                                                            scene2 = list(xaxis = list(title = "X",nticks = x_ticks),
                                                                                          yaxis = list(title = "Y",nticks = y_ticks),
                                                                                          zaxis = list(title = "Z",nticks = z_ticks),
                                                                                          aspectmode='manual',
                                                                                          aspectratio = list(x=ratio[[1]],
                                                                                                             y=ratio[[2]],
                                                                                                             z=ratio[[3]])),
                                                                            scene3 = list(xaxis = list(title = "X",nticks = x_ticks),
                                                                                          yaxis = list(title = "Y",nticks = y_ticks),
                                                                                          zaxis = list(title = "Z",nticks = z_ticks),
                                                                                          aspectmode='manual',
                                                                                          aspectratio = list(x=ratio[[1]],
                                                                                                             y=ratio[[2]],
                                                                                                             z=ratio[[3]])),
                                                                            legend = list(x = 100, y = 0)))
  }
  else if(length(savelist)==4){
    
    
    cowplot <- suppressWarnings(plotly::subplot(savelist)%>% plotly::layout(scene = list(xaxis = list(title = "X",nticks = x_ticks),
                                                                                         yaxis = list(title = "Y",nticks = y_ticks),
                                                                                         zaxis = list(title = "Z",nticks = z_ticks),
                                                                                         aspectmode='manual',
                                                                                         aspectratio = list(x=ratio[[1]],
                                                                                                            y=ratio[[2]],
                                                                                                            z=ratio[[3]])),
                                                                            scene2 = list(xaxis = list(title = "X",nticks = x_ticks),
                                                                                          yaxis = list(title = "Y",nticks = y_ticks),
                                                                                          zaxis = list(title = "Z",nticks = z_ticks),
                                                                                          aspectmode='manual',
                                                                                          aspectratio = list(x=ratio[[1]],
                                                                                                             y=ratio[[2]],
                                                                                                             z=ratio[[3]])),
                                                                            scene3 = list(xaxis = list(title = "X",nticks = x_ticks),
                                                                                          yaxis = list(title = "Y",nticks = y_ticks),
                                                                                          zaxis = list(title = "Z",nticks = z_ticks),
                                                                                          aspectmode='manual',
                                                                                          aspectratio = list(x=ratio[[1]],
                                                                                                             y=ratio[[2]],
                                                                                                             z=ratio[[3]])),
                                                                            scene4 = list(xaxis = list(title = "X",nticks = x_ticks),
                                                                                          yaxis = list(title = "Y",nticks = y_ticks),
                                                                                          zaxis = list(title = "Z",nticks = z_ticks),
                                                                                          aspectmode='manual',
                                                                                          aspectratio = list(x=ratio[[1]],
                                                                                                             y=ratio[[2]],
                                                                                                             z=ratio[[3]])),
                                                                            legend = list(x = 100, y = 0)))
  }
  
  return(cowplot)
  
  }


#' @title visGenePlot
#' @name visGenePlot
#' @description Visualize cells and gene expression according to spatial coordinates
#' @param gobject giotto object
#' @param expression_values gene expression values to use
#' @param genes genes to show
#' @param genes_high_color color represents high gene expression
#' @param genes_mid_color color represents middle gene expression
#' @param genes_low_color color represents low gene expression
#' @param show_network show underlying spatial network
#' @param network_color color of spatial network
#' @param spatial_network_name name of spatial network to use
#' @param show_grid show spatial grid
#' @param grid_color color of spatial grid
#' @param spatial_grid_name name of spatial grid to use
#' @param midpoint expression midpoint
#' @param scale_alpha_with_expression scale expression with ggplot alpha parameter
#' @param point_size size of point (cell)
#' @param point_border_col color of border around points
#' @param point_border_stroke stroke size of border around points
#' @param show_legend show legend
#' @param cow_n_col cowplot param: how many columns
#' @param cow_rel_h cowplot param: relative height
#' @param cow_rel_w cowplot param: relative width
#' @param cow_align cowplot param: how to align
#' @param axis_scale three mode to adjust axis scale
#' @param x_ticks number of ticks on x axis
#' @param y_ticks number of ticks on y axis
#' @param z_ticks number of ticks on z axis
#' @param plot_method two methods of plot
#' @param show_plots show plots
#' @return ggplot or plotly
#' @details Description of parameters.
#' @export
#' @examples
#'     visGenePlot(gobject)
visGenePlot <- function(gobject,
                        expression_values = c('normalized', 'scaled', 'custom'),
                        genes,
                        genes_high_color = NULL,
                        genes_mid_color = "white",
                        genes_low_color = "blue",
                        show_network = F,
                        network_color = NULL,
                        spatial_network_name = 'spatial_network',
                        edge_alpha = NULL,
                        show_grid = F,
                        grid_color = NULL,
                        spatial_grid_name = 'spatial_grid',
                        midpoint = 0,
                        scale_alpha_with_expression = FALSE,
                        point_size = 1,
                        point_border_col = 'black',
                        point_border_stroke = 0.1,
                        show_legend = T,
                        cow_n_col = 2,
                        cow_rel_h = 1,
                        cow_rel_w = 1,
                        cow_align = 'h',
                        axis_scale = c("cube","real","custom"),
                        custom_ratio = NULL,
                        x_ticks = NULL,
                        y_ticks = NULL,
                        z_ticks = NULL,
                        plot_method = c( 'ggplot','plotly'),
                        show_plots = F){
  
  
  plot_method = match.arg(plot_method, choices = c('ggplot','plotly'))
  
  if(plot_method == 'ggplot'){
    
    if(is.null(genes_high_color)){
      genes_high_color = "red"
    }
    result = visGenePlot_2D_ggplot(gobject = gobject,
                                   expression_values = expression_values,
                                   genes = genes,
                                   genes_high_color = genes_high_color,
                                   genes_mid_color = genes_mid_color,
                                   genes_low_color = genes_low_color,
                                   show_network = show_network,
                                   network_color = network_color,
                                   spatial_network_name = spatial_network_name,
                                   edge_alpha = edge_alpha,
                                   show_grid = show_grid,
                                   grid_color = grid_color,
                                   spatial_grid_name = spatial_grid_name,
                                   midpoint = midpoint,
                                   scale_alpha_with_expression = scale_alpha_with_expression,
                                   point_size = point_size,
                                   point_border_col = point_border_col,
                                   point_border_stroke = point_border_stroke,
                                   show_legend = show_legend,
                                   cow_n_col = cow_n_col,
                                   cow_rel_h = cow_rel_h,
                                   cow_rel_w = cow_rel_w,
                                   cow_align = cow_align,
                                   show_plots = show_plots)
  }
  else{
    
    result = visGenePlot_3D_plotly(gobject,
                                   expression_values = expression_values,
                                   genes = genes,
                                   show_network = show_network,
                                   network_color = network_color,
                                   spatial_network_name = spatial_network_name,
                                   edge_alpha = edge_alpha,
                                   show_grid = show_grid,
                                   genes_high_color = genes_high_color,
                                   genes_mid_color = genes_mid_color,
                                   genes_low_color = genes_low_color,
                                   spatial_grid_name = spatial_grid_name,
                                   point_size = point_size,
                                   show_legend = show_legend,
                                   axis_scale = axis_scale,
                                   custom_ratio = custom_ratio,
                                   x_ticks = x_ticks,
                                   y_ticks = y_ticks,
                                   z_ticks = z_ticks,
                                   show_plots = show_plots)
  }
  return(result)
}








#' @title visDimPlot_2D_ggplot
#' @name visDimPlot_2D_ggplot
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
#' @param select_cell_groups select subset of cells/clusters based on cell_color parameter
#' @param select_cells select subset of cells based on cell IDs
#' @param show_other_cells display not selected cells
#' @param other_cell_color color of not selected cells
#' @param other_point_size size of not selected cells
#' @param show_cluster_center plot center of selected clusters
#' @param show_center_label plot label of selected clusters
#' @param center_point_size size of center points
#' @param label_size  size of labels
#' @param label_fontface font of labels
#' @param edge_alpha column to use for alpha of the edges
#' @param point_size size of point (cell)
#' @param point_border_col color of border around points
#' @param point_border_stroke stroke size of border around points
#' @param show_legend show legend
#' @return ggplot
#' @details Description of parameters.
#' @export
#' @examples
#'     visDimPlot_2D_ggplot(gobject)
visDimPlot_2D_ggplot <- function(gobject,
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
                                 select_cell_groups = NULL,
                                 select_cells = NULL,
                                 show_other_cells = T,
                                 other_cell_color = 'lightgrey',
                                 other_point_size = 0.5,
                                 show_cluster_center = F,
                                 show_center_label = T,
                                 center_point_size = 4,
                                 center_point_border_col = 'black',
                                 center_point_border_stroke = 0.1,
                                 label_size = 4,
                                 label_fontface = 'bold',
                                 edge_alpha = NULL,
                                 point_size = 1,
                                 point_border_col = 'black',
                                 point_border_stroke = 0.1,
                                 show_legend = T,
                                 show_plot = F,
                                 return_plot = TRUE,
                                 save_plot = F,
                                 save_dir = NULL,
                                 save_folder = NULL,
                                 save_name = NULL,
                                 save_format = NULL,
                                 show_saved_plot = F,
                                 ...){
  
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
  
  if(dim_reduction_to_use == "pca"){
    eigenvaluesDT = data.table::as.data.table(gobject@dimension_reduction$cells[[dim_reduction_to_use]][[dim_reduction_name]]$misc$eig)
    var_expl_vec = eigenvaluesDT[c(dim1_to_use, dim2_to_use)][['percentage of variance']]
    dim1_x_variance = var_expl_vec[1]
    dim2_y_variance = var_expl_vec[2]
  }
  
  
  
  ## create subsets if needed
  if(!is.null(select_cells) & !is.null(select_cell_groups)) {
    if(is.null(cell_color)) {
      stop('\n selection of cells is based on cell_color paramter, which is a metadata column \n')
    }
    cat('You have selected both individual cell IDs and a group of cells \n')
    group_cell_IDs = annotated_DT[get(cell_color) %in% select_cell_groups][['cell_ID']]
    select_cells = unique(c(select_cells, group_cell_IDs))
  } else if(!is.null(select_cell_groups)) {
    select_cells = annotated_DT[get(cell_color) %in% select_cell_groups][['cell_ID']]
  }
  
  if(!is.null(select_cells)) {
    annotated_DT_other = annotated_DT[!annotated_DT$cell_ID %in% select_cells]
    annotated_DT_selected = annotated_DT[annotated_DT$cell_ID %in% select_cells]
    
    if(show_NN_network == TRUE) {
      annotated_network_DT <- annotated_network_DT[annotated_network_DT$to %in% select_cells & annotated_network_DT$from %in% select_cells]
    }
    
    # if specific cells are selected
    annotated_DT = annotated_DT_selected
  }
  
  ## if no subsets are required
  if(is.null(select_cells) & is.null(select_cell_groups)) {
    annotated_DT_selected = annotated_DT
    annotated_DT_other    = NULL
  }
  
  
  
  pl <- ggplot2::ggplot()
  pl <- pl + ggplot2::theme_classic()
  
  
  
  ## add network layer
  if(show_NN_network == TRUE) {
    pl = plot_network_layer_ggplot(ggobject = pl,
                                   annotated_network_DT = annotated_network_DT,
                                   edge_alpha = edge_alpha,
                                   show_legend = show_legend)
  }
  
  
  
  ## add point layer
  pl = plot_point_layer_ggplot(ggobject = pl,
                               annotated_DT_selected = annotated_DT_selected,
                               annotated_DT_other = annotated_DT_other,
                               cell_color = cell_color,
                               color_as_factor = color_as_factor,
                               cell_color_code = cell_color_code,
                               select_cell_groups = select_cell_groups,
                               select_cells = select_cells,
                               show_other_cells = show_other_cells,
                               other_cell_color = other_cell_color,
                               other_point_size = other_point_size,
                               show_cluster_center = show_cluster_center,
                               show_center_label = show_center_label,
                               center_point_size = center_point_size,
                               center_point_border_col = center_point_border_col,
                               center_point_border_stroke = center_point_border_stroke,
                               label_size = label_size,
                               label_fontface = label_fontface,
                               edge_alpha = edge_alpha,
                               point_size = point_size,
                               point_border_col = point_border_col,
                               point_border_stroke = point_border_stroke,
                               show_legend = show_legend)
  
  
  ## add % variance explained to names of plot for PCA ##
  if(dim_reduction_to_use == 'pca') {
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
  
  
  ## print plot
  if(show_plot == TRUE) {
    print(pl)
  }
  
  
  ## save plot
  if(save_plot == TRUE) {
    
    ggplot_save_function(gobject = gobject,
                         plot_object = pl,
                         save_dir = save_dir,
                         save_folder = save_folder,
                         save_name = save_name,
                         save_format = save_format,
                         show_saved_plot = show_saved_plot,
                         ...)
  }
  
  ## return plot
  if(return_plot == TRUE) return(pl)
  
}


#' @title visDimPlot_2D_plotly
#' @name visDimPlot_2D_plotly
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
#' @param center_point_size size of center points
#' @param label_size  size of labels
#' @param edge_alpha column to use for alpha of the edges
#' @param point_size size of point (cell)
#' @return plotly
#' @details Description of parameters.
#' @export
#' @examples
#'     visDimPlot_2D_plotly(gobject)
visDimPlot_2D_plotly <- function(gobject,
                                 dim_reduction_to_use = 'umap',
                                 dim_reduction_name = 'umap',
                                 dim1_to_use = 1,
                                 dim2_to_use = 2,
                                 select_cell_groups = NULL,
                                 select_cells = NULL,
                                 show_other_cells = T,
                                 other_cell_color = 'lightgrey',
                                 other_point_size = 0.5,
                                 show_NN_network = F,
                                 nn_network_to_use = 'sNN',
                                 network_name = 'sNN.pca',
                                 color_as_factor = T,
                                 cell_color = NULL,
                                 cell_color_code = NULL,
                                 show_cluster_center = F,
                                 show_center_label = T,
                                 center_point_size = 4,
                                 label_size = 4,
                                 edge_alpha = NULL,
                                 point_size = 5){
  
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
  
  if(dim_reduction_to_use == "pca"){
    eigenvaluesDT = data.table::as.data.table(gobject@dimension_reduction$cells[[dim_reduction_to_use]][[dim_reduction_name]]$misc$eig)
    var_expl_vec = eigenvaluesDT[c(dim1_to_use, dim2_to_use)][['percentage of variance']]
    dim1_x_variance = var_expl_vec[1]
    dim2_y_variance = var_expl_vec[2]
  }
  
  if(!is.null(select_cells) & !is.null(select_cell_groups)) {
    if(is.null(cell_color)) {
      stop('\n selection of cells is based on cell_color paramter, which is a metadata column \n')
    }
    cat('You have selected both individual cell IDs and a group of cells \n')
    group_cell_IDs = annotated_DT[get(cell_color) %in% select_cell_groups][['cell_ID']]
    select_cells = unique(c(select_cells, group_cell_IDs))
  } else if(!is.null(select_cell_groups)) {
    select_cells = annotated_DT[get(cell_color) %in% select_cell_groups][['cell_ID']]
  }
  
  if(!is.null(select_cells)) {
    annotated_DT_other = annotated_DT[!annotated_DT$cell_ID %in% select_cells]
    annotated_DT_selected = annotated_DT[annotated_DT$cell_ID %in% select_cells]
    
    if(show_NN_network == TRUE) {
      annotated_network_DT <- annotated_network_DT[annotated_network_DT$to %in% select_cells & annotated_network_DT$from %in% select_cells]
    }
    
    # if specific cells are selected
    annotated_DT = annotated_DT_selected
  }
  
  ## if no subsets are required
  if(is.null(select_cells) & is.null(select_cell_groups)) {
    annotated_DT_selected = annotated_DT
    annotated_DT_other    = NULL
  }
  
  pl <- plotly::plot_ly()
  if(show_NN_network == TRUE) {
    if(is.null(edge_alpha)) {
      edge_alpha = 0.5
    }
    else if(is.character(edge_alpha)){
      warning("Edge_alpha for plotly mode is not adjustable yet. Default 0.5 will be set\n")
      edge_alpha = 0.5
    }
    pl <- pl %>% plotly::add_segments(name = network_name,
                                      type = "scatter",
                                      x = annotated_network_DT[[from_dim_names[1]]],
                                      y = annotated_network_DT[[from_dim_names[2]]],
                                      xend = annotated_network_DT[[to_dim_names[1]]],
                                      yend = annotated_network_DT[[to_dim_names[2]]],
                                      line = list(color = "lightgray",
                                                  width = 0.5),
                                      opacity=edge_alpha)
  }
  
  if(is.null(cell_color)){
    cell_color = "lightblue"
    pl <- pl %>% plotly::add_trace(type = "scatter",mode = "markers",
                                   x = annotated_DT[[dim_names[1]]],
                                   y = annotated_DT[[dim_names[2]]],
                                   color = cell_color,
                                   colors = cell_color,
                                   marker = list(size = point_size))
  }
  
  else if(cell_color %in% colnames(annotated_DT)){
    if(is.null(cell_color_code)){
      number_colors=length(unique(annotated_DT[[cell_color]]))
      cell_color_code = Giotto:::getDistinctColors(n = number_colors)
    }
    if(color_as_factor){
      annotated_DT[[cell_color]] <- as.factor(annotated_DT[[cell_color]])
    }
    
    
    pl <- pl %>% plotly::add_trace(type = "scatter",mode = "markers",
                                   x = annotated_DT[[dim_names[1]]],
                                   y = annotated_DT[[dim_names[2]]],
                                   color = annotated_DT[[cell_color]],
                                   colors = cell_color_code,
                                   legendgroup = annotated_DT[[cell_color]],
                                   marker = list(size = point_size))
    
    if(!is.null(select_cells)&show_other_cells){
      pl <- pl %>% plotly::add_trace(type = "scatter",mode = "markers",
                                     x = annotated_DT_other[[dim_names[1]]],
                                     y = annotated_DT_other[[dim_names[2]]],
                                     #legendgroup = annotated_DT[[cell_color]],
                                     marker = list(size = other_point_size,color = other_cell_color),
                                     showlegend = F)
    }
    
    if(show_cluster_center == TRUE | show_center_label == TRUE) {
      annotated_DT_centers = annotated_DT[, .(center_1 = median(get(dim_names[1])),
                                              center_2 = median(get(dim_names[2]))),
                                          by = cell_color]
      annotated_DT_centers[[cell_color]] <- as.factor(annotated_DT_centers[[cell_color]])
      if(show_cluster_center == TRUE){
        pl <- pl %>% plotly::add_trace(type = "scatter",mode = "markers",
                                       x = annotated_DT_centers[["center_1"]],
                                       y = annotated_DT_centers[["center_2"]],
                                       color = annotated_DT_centers[[cell_color]],
                                       colors = cell_color_code,
                                       legendgroup = annotated_DT_centers[[cell_color]],
                                       marker = list(size = center_point_size,symbol = "x",symbols = "x"),
                                       showlegend = F)
      }
      
      if(show_center_label == TRUE){
        pl <- pl %>%  plotly::add_text(x = annotated_DT_centers[["center_1"]],
                                       y = annotated_DT_centers[["center_2"]],
                                       type = 'scatter',mode = 'text',
                                       text = annotated_DT_centers[[cell_color]],
                                       textposition = 'middle right',
                                       textfont = list(color = '#000000', size = 16),
                                       showlegend = F)
      }
      
    }
  }
  
  else{
    stop("cell_color does not exist!\n")
  }
  
  
  
  if(dim_reduction_to_use == 'pca') {
    x_name = paste0('pca','-',dim_names[1])
    y_name = paste0('pca','-',dim_names[2])
    x_title = sprintf('%s explains %.02f%% of variance', x_name, var_expl_vec[1])
    y_title = sprintf('%s explains %.02f%% of variance', y_name, var_expl_vec[2])
  }
  else{
    x_title = paste(dim_reduction_to_use,dim_names[1],sep = " ")
    y_title = paste(dim_reduction_to_use,dim_names[2],sep = " ")
  }
  pl <- pl %>% plotly::layout(xaxis = list(title = x_title),
                              yaxis = list(title = y_title),
                              legend = list(x = 100, y = 0.5,font = list(family = "sans-serif",size = 12)))
  
  return (pl)
}


#' @title visDimPlot_3D_plotly
#' @name visDimPlot_3D_plotly
#' @description Visualize cells according to dimension reduction coordinates
#' @param gobject giotto object
#' @param dim_reduction_to_use dimension reduction to use
#' @param dim_reduction_name dimension reduction name
#' @param dim1_to_use dimension to use on x-axis
#' @param dim2_to_use dimension to use on y-axis
#' @param dim3_to_use dimension to use on z-axis
#' @param show_NN_network show underlying NN network
#' @param nn_network_to_use type of NN network to use (kNN vs sNN)
#' @param network_name name of NN network to use, if show_NN_network = TRUE
#' @param cell_color color for cells (see details)
#' @param color_as_factor convert color column to factor
#' @param cell_color_code named vector with colors
#' @param show_cluster_center plot center of selected clusters
#' @param show_center_label plot label of selected clusters
#' @param center_point_size size of center points
#' @param label_size  size of labels
#' @param edge_alpha column to use for alpha of the edges
#' @param point_size size of point (cell)
#' @return plotly
#' @details Description of parameters.
#' @export
#' @examples
#'     visDimPlot_3D_plotly(gobject)

visDimPlot_3D_plotly <- function(gobject,
                                 dim_reduction_to_use = 'umap',
                                 dim_reduction_name = 'umap',
                                 dim1_to_use = 1,
                                 dim2_to_use = 2,
                                 dim3_to_use = 3,
                                 
                                 select_cell_groups = NULL,
                                 select_cells = NULL,
                                 show_other_cells = T,
                                 other_cell_color = 'lightgrey',
                                 other_point_size = 0.5,
                                 
                                 show_NN_network = F,
                                 nn_network_to_use = 'sNN',
                                 network_name = 'sNN.pca',
                                 color_as_factor = T,
                                 cell_color = NULL,
                                 cell_color_code = NULL,
                                 show_cluster_center = F,
                                 show_center_label = T,
                                 center_point_size = 4,
                                 label_size = 4,
                                 edge_alpha = NULL,
                                 point_size = 1){
  ## dimension reduction ##
  dim_dfr = gobject@dimension_reduction$cells[[dim_reduction_to_use]][[dim_reduction_name]]$coordinates[,c(dim1_to_use, dim2_to_use,dim3_to_use)]
  dim_names = colnames(dim_dfr)
  dim_DT = data.table::as.data.table(dim_dfr); dim_DT[, cell_ID := rownames(dim_dfr)]
  
  ## annotated cell metadata
  cell_metadata = gobject@cell_metadata
  annotated_DT = data.table::merge.data.table(cell_metadata, dim_DT, by = 'cell_ID')
  
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
  
  if(dim_reduction_to_use == "pca"){
    eigenvaluesDT = data.table::as.data.table(gobject@dimension_reduction$cells[[dim_reduction_to_use]][[dim_reduction_name]]$misc$eig)
    var_expl_vec = eigenvaluesDT[c(dim1_to_use, dim2_to_use,dim3_to_use)][['percentage of variance']]
    dim1_x_variance = var_expl_vec[1]
    dim2_y_variance = var_expl_vec[2]
    dim3_z_variance = var_expl_vec[3]
  }
  
  ## create subsets if needed
  if(!is.null(select_cells) & !is.null(select_cell_groups)) {
    if(is.null(cell_color)) {
      stop('\n selection of cells is based on cell_color paramter, which is a metadata column \n')
    }
    cat('You have selected both individual cell IDs and a group of cells \n')
    group_cell_IDs = annotated_DT[get(cell_color) %in% select_cell_groups][['cell_ID']]
    select_cells = unique(c(select_cells, group_cell_IDs))
  } else if(!is.null(select_cell_groups)) {
    select_cells = annotated_DT[get(cell_color) %in% select_cell_groups][['cell_ID']]
  }
  
  if(!is.null(select_cells)) {
    annotated_DT_other = annotated_DT[!annotated_DT$cell_ID %in% select_cells]
    annotated_DT_selected = annotated_DT[annotated_DT$cell_ID %in% select_cells]
    
    if(show_NN_network == TRUE) {
      annotated_network_DT <- annotated_network_DT[annotated_network_DT$to %in% select_cells & annotated_network_DT$from %in% select_cells]
    }
    
    # if specific cells are selected
    annotated_DT = annotated_DT_selected
  }
  
  ## if no subsets are required
  if(is.null(select_cells) & is.null(select_cell_groups)) {
    annotated_DT_selected = annotated_DT
    annotated_DT_other    = NULL
  }
  
  pl <- plotly::plot_ly()
  if(is.null(cell_color)){
    cell_color = "lightblue"
    pl <- pl %>% plotly::add_trace(type = 'scatter3d',mode = "markers",
                                   x = annotated_DT[[dim_names[1]]],
                                   y = annotated_DT[[dim_names[2]]],
                                   z = annotated_DT[[dim_names[3]]],
                                   color = cell_color,
                                   colors = cell_color,
                                   marker = list(size = 2),
                                   legendgroup = annotated_DT[[cell_color]])
  }
  else{
    if(cell_color %in% colnames(annotated_DT)){
      if(is.null(cell_color_code)) {
        number_colors=length(unique(annotated_DT[[cell_color]]))
        cell_color_code = Giotto:::getDistinctColors(n = number_colors)
      }
      if(color_as_factor){
        annotated_DT[[cell_color]] <- as.factor(annotated_DT[[cell_color]])
      }
      
      pl <- pl %>% plotly::add_trace(type = 'scatter3d',mode = "markers",
                                     x = annotated_DT[[dim_names[1]]],
                                     y = annotated_DT[[dim_names[2]]],
                                     z = annotated_DT[[dim_names[3]]],
                                     color = annotated_DT[[cell_color]],
                                     colors = cell_color_code,
                                     marker = list(size = point_size),
                                     legendgroup = annotated_DT[[cell_color]])
      
      if(!is.null(select_cells)&show_other_cells){
        pl <- pl %>% plotly::add_trace(type = 'scatter3d',mode = "markers",
                                       x = annotated_DT_other[[dim_names[1]]],
                                       y = annotated_DT_other[[dim_names[2]]],
                                       z = annotated_DT_other[[dim_names[3]]],
                                       #colors = other_cell_color,
                                       marker = list(size = other_point_size,color = other_cell_color),
                                       showlegend = F)
      }
      
      
      if(show_cluster_center == TRUE | show_center_label == TRUE){
        annotated_DT_centers = annotated_DT[, .(center_1 = median(get(dim_names[1])),
                                                center_2 = median(get(dim_names[2])),
                                                center_3 = median(get(dim_names[3]))),
                                            by = cell_color]
        annotated_DT_centers[[cell_color]] <- as.factor(annotated_DT_centers[[cell_color]])
        if(show_cluster_center == TRUE){
          pl <- pl %>% plotly::add_trace(mode = "markers",
                                         type = "scatter3d",
                                         data = annotated_DT_centers,
                                         x = ~center_1,
                                         y = ~center_2,
                                         z = ~center_3,
                                         color = annotated_DT_centers[[cell_color]],
                                         colors = cell_color_code,
                                         inherit = F,
                                         marker=list(size = 2,symbol = "x",symbols = "x"),
                                         legendgroup = annotated_DT_centers[[cell_color]],
                                         showlegend = F)
        }
        if(show_center_label == TRUE){
          cat(" center label is not clear to see in 3D plot\n You can shut it down with show_center_label = F\n")
          pl <- pl %>% plotly::add_trace(mode = "text",
                                         type = "scatter3d",
                                         data = annotated_DT_centers,
                                         x = ~center_1,
                                         y = ~center_2,
                                         z = ~center_3,
                                         text = annotated_DT_centers[[cell_color]],
                                         legendgroup = annotated_DT_centers[[cell_color]],
                                         inherit = F,
                                         showlegend = F)
        }
        
      }
    }
    
    else{
      stop("cell_color does not exist!\n")
    }
  }
  
  if(show_NN_network){
    edges <- plotly_network(annotated_network_DT,
                            "from_Dim.1","from_Dim.2","from_Dim.3",
                            "to_Dim.1","to_Dim.2","to_Dim.3")
    if(is.null(edge_alpha)){
      edge_alpha = 0.5
    }
    else if(is.character(edge_alpha)){
      warning("Edge_alpha for plotly mode is not adjustable yet. Default 0.5 will be set\n")
      edge_alpha = 0.5
    }
    
    pl <- pl %>% plotly::add_trace(name = network_name,
                                   mode = "lines",
                                   type = "scatter3d",
                                   data = edges,
                                   x = ~x,y=~y,z=~z,
                                   inherit = F,
                                   line=list(color="lightgray", width = 0.5),
                                   opacity=edge_alpha)
  }
  
  if(dim_reduction_to_use == 'pca') {
    x_name = paste0('pca','-',dim_names[1])
    y_name = paste0('pca','-',dim_names[2])
    z_name = paste0('pca','-',dim_names[3])
    x_title = sprintf('%s explains %.02f%% of variance', x_name, var_expl_vec[1])
    y_title = sprintf('%s explains %.02f%% of variance', y_name, var_expl_vec[2])
    z_title = sprintf('%s explains %.02f%% of variance', z_name, var_expl_vec[3])
  }
  else{
    x_title = paste(dim_reduction_to_use,dim_names[1],sep = " ")
    y_title = paste(dim_reduction_to_use,dim_names[2],sep = " ")
    z_title = paste(dim_reduction_to_use,dim_names[3],sep = " ")
  }
  pl <- pl %>%  plotly::layout(scene = list(xaxis = list(title = x_title),
                                            yaxis = list(title = y_title),
                                            zaxis = list(title = z_title)),
                               legend = list(x = 100, y = 0.5,font = list(family = "sans-serif",size = 12)))
  return(pl)
}



#' @title visDimPlot
#' @name visDimPlot
#' @description Visualize cells according to dimension reduction coordinates
#' @param gobject giotto object
#' @param dim_reduction_to_use dimension reduction to use
#' @param dim_reduction_name dimension reduction name
#' @param dim1_to_use dimension to use on x-axis
#' @param dim2_to_use dimension to use on y-axis
#' @param dim3_to_use dimension to use on z-axis
#' @param show_NN_network show underlying NN network
#' @param nn_network_to_use type of NN network to use (kNN vs sNN)
#' @param network_name name of NN network to use, if show_NN_network = TRUE
#' @param cell_color color for cells (see details)
#' @param color_as_factor convert color column to factor
#' @param cell_color_code named vector with colors
#' @param show_cluster_center plot center of selected clusters
#' @param show_center_label plot label of selected clusters
#' @param center_point_size size of center points
#' @param label_size  size of labels
#' @param label_fontface font of labels
#' @param edge_alpha column to use for alpha of the edges
#' @param point_size size of point (cell)
#' @param point_border_col color of border around points
#' @param point_border_stroke stroke size of border around points
#' @param show_legend show legend
#' @param show_plot show plot
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_dir directory to save the plot
#' @param save_folder (optional) folder in directory to save the plot
#' @param save_name name of plot
#' @param save_format format of plot (e.g. tiff, png, pdf, ...)
#' @param show_saved_plot load & display the saved plot
#' @return ggplot or plotly
#' @details Description of parameters.
#' @export
#' @examples
#'     visDimPlot(gobject)
visDimPlot <- function(gobject,
                       dim_reduction_to_use = 'umap',
                       dim_reduction_name = 'umap',
                       dim1_to_use = 1,
                       dim2_to_use = 2,
                       dim3_to_use = NULL,
                       show_NN_network = F,
                       nn_network_to_use = 'sNN',
                       network_name = 'sNN.pca',
                       cell_color = NULL,
                       color_as_factor = T,
                       cell_color_code = NULL,
                       select_cell_groups = NULL,
                       select_cells = NULL,
                       show_other_cells = T,
                       other_cell_color = 'lightgrey',
                       other_point_size = 0.5,
                       show_cluster_center = F,
                       show_center_label = T,
                       center_point_size = 4,
                       center_point_border_col = 'black',
                       center_point_border_stroke = 0.1,
                       label_size = 4,
                       label_fontface = 'bold',
                       edge_alpha = NULL,
                       point_size = 3,
                       point_border_col = 'black',
                       point_border_stroke = 0.1,
                       plot_method = c('ggplot', 'plotly'),
                       show_legend = T,
                       show_plot = F,
                       return_plot = TRUE,
                       save_plot = F,
                       save_dir = NULL,
                       save_folder = NULL,
                       save_name = NULL,
                       save_format = NULL,
                       show_saved_plot = F,
                       ...){
  
  plot_method = match.arg(plot_method, choices = c('ggplot', 'plotly'))
  
  
  if(plot_method == "ggplot"){
    if(!is.null(dim3_to_use)) {
      warning('\n ggplot is only in 2D. If you want to plot the third dim, please choose plotly method with \"plot_method = \"plotly\"\"\n')
    }
    
    result = visDimPlot_2D_ggplot(gobject = gobject,
                                  dim_reduction_to_use = dim_reduction_to_use,
                                  dim_reduction_name = dim_reduction_name,
                                  dim1_to_use = dim1_to_use,
                                  dim2_to_use = dim2_to_use,
                                  show_NN_network = show_NN_network,
                                  nn_network_to_use = nn_network_to_use,
                                  network_name = network_name,
                                  cell_color = cell_color,
                                  color_as_factor = color_as_factor,
                                  cell_color_code = cell_color_code,
                                  select_cell_groups = select_cell_groups,
                                  select_cells = select_cells,
                                  show_other_cells = show_other_cells,
                                  other_cell_color = other_cell_color,
                                  other_point_size = other_point_size,
                                  show_cluster_center = show_cluster_center,
                                  show_center_label = show_center_label,
                                  center_point_size = center_point_size,
                                  center_point_border_col = center_point_border_col,
                                  center_point_border_stroke = center_point_border_stroke,
                                  label_size = label_size,
                                  label_fontface = label_fontface,
                                  edge_alpha = edge_alpha,
                                  point_size = point_size,
                                  point_border_col =point_border_col,
                                  point_border_stroke = point_border_stroke,
                                  show_legend = show_legend,
                                  show_plot = show_plot,
                                  return_plot = return_plot,
                                  save_plot = save_plot,
                                  save_dir = save_dir,
                                  save_folder = save_folder,
                                  save_name = save_name,
                                  save_format = save_format,
                                  show_saved_plot = show_saved_plot,
                                  ...)
  }
  
  
  else if(plot_method == "plotly"){
    if(is.null(dim3_to_use)){
      result = visDimPlot_2D_plotly(gobject = gobject,
                                    dim_reduction_to_use = dim_reduction_to_use,
                                    dim_reduction_name = dim_reduction_name,
                                    dim1_to_use = dim1_to_use,
                                    dim2_to_use = dim2_to_use,
                                    show_NN_network = show_NN_network,
                                    nn_network_to_use = nn_network_to_use,
                                    network_name = network_name,
                                    cell_color = cell_color,
                                    color_as_factor = color_as_factor,
                                    cell_color_code = cell_color_code,
                                    show_cluster_center = show_cluster_center,
                                    show_center_label = show_center_label,
                                    center_point_size = center_point_size,
                                    label_size = label_size,
                                    edge_alpha = edge_alpha,
                                    point_size = point_size)
    }
    
    else{
      result = visDimPlot_3D_plotly(gobject = gobject,
                                    dim_reduction_to_use = dim_reduction_to_use,
                                    dim_reduction_name = dim_reduction_name,
                                    dim1_to_use = dim1_to_use,
                                    dim2_to_use = dim2_to_use,
                                    dim3_to_use = dim3_to_use,
                                    show_NN_network = show_NN_network,
                                    nn_network_to_use = nn_network_to_use,
                                    network_name =network_name,
                                    cell_color = cell_color,
                                    color_as_factor = color_as_factor,
                                    cell_color_code = cell_color_code,
                                    show_cluster_center = show_cluster_center,
                                    show_center_label = show_center_label,
                                    center_point_size = center_point_size,
                                    label_size = label_size,
                                    edge_alpha = edge_alpha,
                                    point_size = point_size)
    }
  }
  return(result)
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
#' @param point_border_col color of border around points
#' @param point_border_stroke stroke size of border around points
#' @param show_legend show legend
#' @param show_plot show plot
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_dir directory to save the plot
#' @param save_folder (optional) folder in directory to save the plot
#' @param save_name name of plot
#' @param save_format format of plot (e.g. tiff, png, pdf, ...)
#' @param show_saved_plot load & display the saved plot
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
                               color_as_factor = TRUE,
                               cell_color_code = NULL,
                               edge_alpha = NULL,
                               point_size = 1,
                               point_border_col = 'black',
                               point_border_stroke = 0.1,
                               show_legend = T,
                               show_plot = F,
                               return_plot = TRUE,
                               save_plot = F,
                               save_dir = NULL,
                               save_folder = NULL,
                               save_name = NULL,
                               save_format = NULL,
                               show_saved_plot = F,
                               ...) {
  
  
  ## layout ##
  co = gobject@nn_network[[nn_network_to_use]][[network_name]][['layout']][,c(dim1_to_use, dim2_to_use)]
  dim_dfr = as.data.frame(co)
  colnames(dim_dfr) = paste0('Dim.', 1:ncol(dim_dfr))
  dim_names = colnames(dim_dfr)
  dim_DT = data.table::as.data.table(dim_dfr); dim_DT[, cell_ID := gobject@cell_ID]
  
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
                                                                               xend = to_dim_names[1], yend = to_dim_names[2]), alpha = edge_alpha, show.legend = show_legend)
    } else if(is.numeric(edge_alpha)) {
      pl <- pl + ggplot2::geom_segment(data = annotated_network_DT, aes_string(x = from_dim_names[1], y = from_dim_names[2],
                                                                               xend = to_dim_names[1], yend = to_dim_names[2]), alpha = edge_alpha, show.legend = show_legend)
    } else if(is.character(edge_alpha)) {
      
      if(edge_alpha %in% colnames(annotated_network_DT)) {
        pl <- pl + ggplot2::geom_segment(data = annotated_network_DT, aes_string(x = from_dim_names[1], y = from_dim_names[2],
                                                                                 xend = to_dim_names[1], yend = to_dim_names[2], alpha = edge_alpha), show.legend = show_legend)
      }
    }
  }
  
  # point layer
  if(is.null(cell_color)) {
    cell_color = 'lightblue'
    pl <- pl + ggplot2::geom_point(data = annotated_DT, aes_string(x = dim_names[1], dim_names[2]),
                                   color = cell_color, show.legend = show_legend, size = point_size)
    
  } else if (is.character(cell_color)) {
    
    if(cell_color %in% colnames(annotated_DT)) {
      
      if(color_as_factor == TRUE) {
        factor_data = factor(annotated_DT[[cell_color]])
        annotated_DT[[cell_color]] <- factor_data
      }
      
      pl <- pl + ggplot2::geom_point(data = annotated_DT, aes_string(x = dim_names[1], y = dim_names[2], fill = cell_color),
                                     show.legend = show_legend, shape = 21, size = point_size,
                                     color = point_border_col, stroke = point_border_stroke)
      
      
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
      
    }
    
  } else {
    pl <- pl + ggplot2::geom_point(data = annotated_DT, aes_string(x = dim_names[1], y = dim_names[2]),
                                   show.legend = show_legend, shape = 21, fill = cell_color,
                                   size = point_size,
                                   color = point_border_col, stroke = point_border_stroke)
  }
  
  
  
  ## print plot
  if(show_plot == TRUE) {
    print(pl)
  }
  
  ## save plot
  if(save_plot == TRUE) {
    
    ggplot_save_function(gobject = gobject,
                         plot_object = pl,
                         save_dir = save_dir,
                         save_folder = save_folder,
                         save_name = save_name,
                         save_format = save_format,
                         show_saved_plot = show_saved_plot,
                         ...)
  }
  
  ## return plot
  if(return_plot == TRUE) return(pl)
  
}




#' @title visSpatDimPlot_2D
#' @name visSpatDimPlot_2D
#' @description Visualize cells according to spatial AND dimension reduction coordinates in ggplot2 mode
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
#' @param select_cell_groups select subset of cells/clusters based on cell_color parameter
#' @param select_cells select subset of cells based on cell IDs
#' @param show_other_cells display not selected cells
#' @param other_cell_color color of not selected cells
#' @param nn_network_alpha column to use for alpha of the edges
#' @param show_spatial_network show spatial network
#' @param spatial_network_name name of spatial network to use
#' @param spatial_network_color color of spatial network
#' @param show_spatial_grid show spatial grid
#' @param spatial_grid_name name of spatial grid to use
#' @param spatial_grid_color color of spatial grid
#' @param show_legend show legend
#' @param show_plot show plot
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_dir directory to save the plot
#' @param save_folder (optional) folder in directory to save the plot
#' @param save_name name of plot
#' @param save_format format of plot (e.g. tiff, png, pdf, ...)
#' @param show_saved_plot load & display the saved plot
#' @return ggplot
#' @details Description of parameters.
#' @export
#' @examples
#'     visSpatDimPlot_2D(gobject)
visSpatDimPlot_2D <- function(gobject,
                              plot_alignment = c('vertical', 'horizontal'),
                              dim_reduction_to_use = 'umap',
                              dim_reduction_name = 'umap',
                              dim1_to_use = 1,
                              dim2_to_use = 2,
                              sdimx=NULL,
                              sdimy=NULL,
                              show_NN_network = F,
                              nn_network_to_use = 'sNN',
                              network_name = 'sNN.pca',
                              show_cluster_center = F,
                              show_center_label = T,
                              center_point_size = 4,
                              label_size = 4,
                              label_fontface = 'bold',
                              cell_color = NULL,
                              color_as_factor = T,
                              cell_color_code = NULL,
                              select_cell_groups = NULL,
                              select_cells = NULL,
                              show_other_cells = T,
                              other_cell_color = 'lightgrey',
                              dim_plot_mode = NULL,
                              dim_point_size = 1,
                              dim_point_border_col = 'black',
                              dim_point_border_stroke = 0.1,
                              nn_network_alpha = 0.05,
                              show_spatial_network = F,
                              spatial_network_name = 'spatial_network',
                              spatial_network_color = NULL,
                              show_spatial_grid = F,
                              spatial_grid_name = 'spatial_grid',
                              spatial_grid_color = NULL,
                              spatial_point_size = 1,
                              spatial_point_border_col = 'black',
                              spatial_point_border_stroke = 0.1,
                              show_legend = T,
                              show_plot = F,
                              plot_method = 'ggplot'){
  
  plot_alignment = match.arg(plot_alignment, choices = c( 'vertical','horizontal'))
  
  
  # create matching cell_color_code
  if(is.null(cell_color_code)) {
    if(is.character(cell_color)) {
      
      cell_metadata = pDataDT(gobject)
      if(cell_color %in% colnames(cell_metadata)) {
        
        if(color_as_factor == TRUE) {
          number_colors = length(unique(cell_metadata[[cell_color]]))
          cell_color_code = Giotto:::getDistinctColors(n = number_colors)
          names(cell_color_code) = unique(cell_metadata[[cell_color]])
          cell_color_code = cell_color_code
        }
      }
    }
  }
  
  # dimension reduction plot
  dmpl = visDimPlot(gobject = gobject,
                    dim_reduction_to_use = dim_reduction_to_use, dim_reduction_name = dim_reduction_name,
                    dim1_to_use = dim1_to_use,
                    dim2_to_use = dim2_to_use,
                    show_NN_network = show_NN_network,
                    nn_network_to_use = nn_network_to_use,
                    network_name = network_name,
                    cell_color = cell_color,
                    color_as_factor = color_as_factor,
                    cell_color_code = cell_color_code,
                    select_cell_groups = select_cell_groups,
                    select_cells = select_cells,
                    show_other_cells = show_other_cells,
                    other_cell_color = other_cell_color,
                    edge_alpha = nn_network_alpha, show_legend = show_legend,
                    point_size = dim_point_size,
                    point_border_col = dim_point_border_col,
                    point_border_stroke = dim_point_border_stroke,
                    show_cluster_center = show_cluster_center,
                    show_center_label = show_center_label,
                    center_point_size = center_point_size,
                    label_size = label_size,
                    label_fontface = label_fontface,
                    plot_method = plot_method)
  
  # spatial plot
  spl = visPlot(gobject = gobject,
                sdimx = sdimx,
                sdimy = sdimy,
                show_network = show_spatial_network,
                spatial_network_name = spatial_network_name,
                show_grid = show_spatial_grid,
                spatial_grid_name = spatial_grid_name,
                cell_color = cell_color,
                cell_color_code = cell_color_code,
                color_as_factor = color_as_factor,
                select_cell_groups = select_cell_groups,
                select_cells = select_cells,
                show_other_cells = show_other_cells,
                other_cell_color = other_cell_color,
                coord_fix_ratio = NULL,
                network_color = spatial_network_color,
                grid_color = spatial_grid_color,
                show_legend = show_legend,
                show_plot = show_plot,
                point_size = spatial_point_size,
                point_border_col = spatial_point_border_col,
                point_border_stroke = spatial_point_border_stroke)
  
  if(plot_alignment == 'vertical') {
    combo_plot <- cowplot::plot_grid(dmpl, spl, ncol = 1, rel_heights = c(1), rel_widths = c(1), align = 'v')
    return(cowplot::plot_grid(combo_plot))
  } else {
    combo_plot <- cowplot::plot_grid(dmpl, spl, ncol = 2, rel_heights = c(1), rel_widths = c(1), align = 'h')
    return(cowplot::plot_grid(combo_plot))
  }
  
}


#' @title visSpatDimPlot_3D
#' @name visSpatDimPlot_3D
#' @description Visualize cells according to spatial AND dimension reduction coordinates in plotly mode
#' @param gobject giotto object
#' @param plot_alignment direction to align plot
#' @param dim_reduction_to_use dimension reduction to use
#' @param dim_reduction_name dimension reduction name
#' @param dim1_to_use dimension to use on x-axis
#' @param dim2_to_use dimension to use on y-axis
#' @param dim3_to_use dimension to use on z-axis
#' @param show_NN_network show underlying NN network
#' @param nn_network_to_use type of NN network to use (kNN vs sNN)
#' @param network_name name of NN network to use, if show_NN_network = TRUE
#' @param cell_color color for cells (see details)
#' @param color_as_factor convert color column to factor
#' @param cell_color_code named vector with colors
#' @param nn_network_alpha column to use for alpha of the edges
#' @param show_spatial_network show spatial network
#' @param spatial_network_name name of spatial network to use
#' @param spatial_network_color color of spatial network
#' @param spatial_network_alpha alpha of spatial network
#' @param show_spatial_grid show spatial grid
#' @param spatial_grid_name name of spatial grid to use
#' @param spatial_grid_color color of spatial grid
#' @param spatial_grid_alpha alpha of spatial grid
#' @param show_legend show legend
#' @param show_plot show plot
#' @param legend_text_size text size of legend
#' @return plotly
#' @details Description of parameters.
#' @export
#' @examples
#'     visSpatDimPlot_3D(gobject)
visSpatDimPlot_3D <- function(gobject,
                              plot_alignment = c('horizontal','vertical'),
                              dim_reduction_to_use = 'umap',
                              dim_reduction_name = 'umap',
                              dim1_to_use = 1,
                              dim2_to_use = 2,
                              dim3_to_use = NULL,
                              sdimx=NULL,
                              sdimy=NULL,
                              sdimz=NULL,
                              show_NN_network = F,
                              nn_network_to_use = 'sNN',
                              network_name = 'sNN.pca',
                              show_cluster_center = F,
                              show_center_label = T,
                              center_point_size = 4,
                              label_size = 16,
                              cell_color = NULL,
                              color_as_factor = T,
                              cell_color_code = NULL,
                              dim_point_size = 3,
                              nn_network_alpha = 0.5,
                              show_spatial_network = F,
                              spatial_network_name = 'spatial_network',
                              network_color = "lightgray",
                              spatial_network_alpha = 0.5,
                              show_spatial_grid = F,
                              spatial_grid_name = 'spatial_grid',
                              spatial_grid_color = NULL,
                              spatial_grid_alpha = 0.5,
                              spatial_point_size = 3,
                              axis_scale = c("cube","real","custom"),
                              custom_ratio = NULL,
                              x_ticks = NULL,
                              y_ticks = NULL,
                              z_ticks = NULL,
                              legend_text_size = 12){
  
  plot_alignment = match.arg(plot_alignment, choices = c( 'horizontal','vertical'))
  
  # ********data prepare********#
  ## dimension reduction ##
  dim_dfr = gobject@dimension_reduction$cells[[dim_reduction_to_use]][[dim_reduction_name]]$coordinates[,c(dim1_to_use, dim2_to_use,dim3_to_use)]
  dim_names = colnames(dim_dfr)
  dim_DT = data.table::as.data.table(dim_dfr); dim_DT[, cell_ID := rownames(dim_dfr)]
  
  ## annotated cell metadata
  cell_metadata = gobject@cell_metadata
  annotated_DT = merge(cell_metadata, dim_DT, by = 'cell_ID')
  annotated_DT = merge(annotated_DT, gobject@spatial_locs,by = 'cell_ID')
  
  if(dim_reduction_to_use == "pca"){
    eigenvaluesDT = data.table::as.data.table(gobject@dimension_reduction$cells[[dim_reduction_to_use]][[dim_reduction_name]]$misc$eig)
    var_expl_vec = eigenvaluesDT[c(dim1_to_use, dim2_to_use,dim3_to_use)][['percentage of variance']]
    dim1_x_variance = var_expl_vec[1]
    dim2_y_variance = var_expl_vec[2]
    if(!is.null(dim3_to_use)){
      dim3_z_variance = var_expl_vec[3]
    }
  }
  
  ## nn network
  if(show_NN_network){
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
  
  
  ## extract spatial network
  if(!is.null(spatial_network_name)) {
    spatial_network = select_spatialNetwork(gobject, name = spatial_network_name, return_network_Obj = FALSE)
  } else {
    spatial_network = NULL
  }
  
  ## extract spatial grid
  if(!is.null(spatial_grid_name)) {
    spatial_grid = gobject@spatial_grid[[spatial_grid_name]]
  } else {
    spatial_grid = NULL
  }
  
  # create matching cell_color_code
  if(is.null(cell_color_code)) {
    if(is.character(cell_color)) {
      
      cell_metadata = pDataDT(gobject)
      if(cell_color %in% colnames(cell_metadata)) {
        
        if(color_as_factor == TRUE) {
          number_colors = length(unique(cell_metadata[[cell_color]]))
          cell_color_code = Giotto:::getDistinctColors(n = number_colors)
          names(cell_color_code) = unique(cell_metadata[[cell_color]])
        }
      }
    }
  }
  
  
  #********** dim plot ***********#
  #2D plot
  if(is.null(dim3_to_use)){
    dpl <- plotly::plot_ly()
    if(show_NN_network == TRUE) {
      if(is.null(nn_network_alpha)) {
        nn_network_alpha = 0.5
      }
      else if(is.character(nn_network_alpha)){
        warning("Edge_alpha for plotly mode is not adjustable yet. Default 0.5 will be set\n")
        nn_network_alpha = 0.5
      }
      dpl <- dpl %>% plotly::add_segments(name = network_name,
                                          type = "scatter",
                                          x = annotated_network_DT[[from_dim_names[1]]],
                                          y = annotated_network_DT[[from_dim_names[2]]],
                                          xend = annotated_network_DT[[to_dim_names[1]]],
                                          yend = annotated_network_DT[[to_dim_names[2]]],
                                          line = list(color = network_color,
                                                      width = 0.5),
                                          opacity=nn_network_alpha)
    }
    
    if(is.null(cell_color)){
      #cell_color = "lightblue"
      dpl <- dpl %>% plotly::add_trace(type = "scatter",mode = "markers",
                                       x = annotated_DT[[dim_names[1]]],
                                       y = annotated_DT[[dim_names[2]]],
                                       #color = "lightblue",
                                       #colors ="lightblue",
                                       marker = list(size = dim_point_size,
                                                     color = "lightblue"),
                                       showlegend = F)
    }
    
    else if(cell_color %in% colnames(annotated_DT)){
      if(color_as_factor){
        annotated_DT[[cell_color]] <- as.factor(annotated_DT[[cell_color]])
      }
      
      
      dpl <- dpl %>% plotly::add_trace(type = "scatter",mode = "markers",
                                       x = annotated_DT[[dim_names[1]]],
                                       y = annotated_DT[[dim_names[2]]],
                                       color = annotated_DT[[cell_color]],
                                       colors = cell_color_code,
                                       legendgroup = annotated_DT[[cell_color]],
                                       marker = list(size = dim_point_size))
    }
    
    else{
      stop("cell_color does not exist!\n")
    }
    
    
    if((show_cluster_center == TRUE | show_center_label == TRUE)&!is.null(cell_color)) {
      annotated_DT_centers = annotated_DT[, .(center_1 = median(get(dim_names[1])),
                                              center_2 = median(get(dim_names[2]))),
                                          by = cell_color]
      annotated_DT_centers[[cell_color]] <- as.factor(annotated_DT_centers[[cell_color]])
      if(show_cluster_center == TRUE){
        dpl <- dpl %>% plotly::add_trace(type = "scatter",mode = "markers",
                                         x = annotated_DT_centers[["center_1"]],
                                         y = annotated_DT_centers[["center_2"]],
                                         color = annotated_DT_centers[[cell_color]],
                                         colors = cell_color_code,
                                         legendgroup = annotated_DT_centers[[cell_color]],
                                         marker = list(size = center_point_size,symbol = "x",symbols = "x"),
                                         showlegend = F)
      }
      
      if(show_center_label == TRUE){
        dpl <- dpl %>%  plotly::add_text(x = annotated_DT_centers[["center_1"]],
                                         y = annotated_DT_centers[["center_2"]],
                                         type = 'scatter',mode = 'text',
                                         text = annotated_DT_centers[[cell_color]],
                                         textposition = 'middle right',
                                         textfont = list(color = '#000000', size = label_size),showlegend = F)
      }
      
    }
    if(dim_reduction_to_use == 'pca') {
      x_name = paste0('pca','-',dim_names[1])
      y_name = paste0('pca','-',dim_names[2])
      x_title = sprintf('%s explains %.02f%% of variance', x_name, var_expl_vec[1])
      y_title = sprintf('%s explains %.02f%% of variance', y_name, var_expl_vec[2])
    }
    else{
      x_title = paste(dim_reduction_to_use,dim_names[1],sep = " ")
      y_title = paste(dim_reduction_to_use,dim_names[2],sep = " ")
    }
    dpl <- dpl %>% plotly::layout(xaxis = list(title = x_title),
                                  yaxis = list(title = y_title),
                                  legend = list(x = 100, y = 0.5,font = list(family = "sans-serif",size = legend_text_size)))
  }
  #3D plot
  else if(!is.null(dim3_to_use)){
    dpl <- plotly::plot_ly(scene = "scene1")
    if(is.null(cell_color)){
      #cell_color = "lightblue"
      dpl <- dpl %>% plotly::add_trace(type = 'scatter3d',mode = "markers",
                                       x = annotated_DT[[dim_names[1]]],
                                       y = annotated_DT[[dim_names[2]]],
                                       z = annotated_DT[[dim_names[3]]],
                                       color = "lightblue",
                                       colors = "lightblue",
                                       marker = list(size = dim_point_size),
                                       showlegend = F)
      #legendgroup = annotated_DT[[cell_color]])
    }
    else{
      if(cell_color %in% colnames(annotated_DT)){
        if(is.null(cell_color_code)) {
          number_colors=length(unique(annotated_DT[[cell_color]]))
          cell_color_code = Giotto:::getDistinctColors(n = number_colors)
        }
        if(color_as_factor){
          annotated_DT[[cell_color]] <- as.factor(annotated_DT[[cell_color]])
        }
        dpl <- dpl %>% plotly::add_trace(type = 'scatter3d',mode = "markers",
                                         x = annotated_DT[[dim_names[1]]],
                                         y = annotated_DT[[dim_names[2]]],
                                         z = annotated_DT[[dim_names[3]]],
                                         color = annotated_DT[[cell_color]],
                                         colors = cell_color_code,
                                         marker = list(size = dim_point_size),
                                         legendgroup = annotated_DT[[cell_color]])
      }
      
      else{
        stop("cell_color does not exist!\n")
      }
    }
    
    if(show_NN_network){
      edges <- plotly_network(annotated_network_DT,
                              "from_Dim.1","from_Dim.2","from_Dim.3",
                              "to_Dim.1","to_Dim.2","to_Dim.3")
      if(is.null(nn_network_alpha)){
        nn_network_alpha = 0.5
      }
      else if(is.character(nn_network_alpha)){
        warning("Edge_alpha for plotly mode is not adjustable yet. Default 0.5 will be set\n")
        nn_network_alpha = 0.5
      }
      
      dpl <- dpl %>% plotly::add_trace(name = network_name,
                                       mode = "lines",
                                       type = "scatter3d",
                                       data = edges,
                                       x = ~x,y=~y,z=~z,
                                       line=list(color=network_color),
                                       opacity=nn_network_alpha)
    }
    if((show_cluster_center == TRUE | show_center_label == TRUE)& !is.null(cell_color)){
      annotated_DT_centers = annotated_DT[, .(center_1 = median(get(dim_names[1])),
                                              center_2 = median(get(dim_names[2])),
                                              center_3 = median(get(dim_names[3]))),
                                          by = cell_color]
      annotated_DT_centers[[cell_color]] <- as.factor(annotated_DT_centers[[cell_color]])
      if(show_cluster_center == TRUE){
        dpl <- dpl %>% plotly::add_trace(mode = "markers",
                                         type = "scatter3d",
                                         data = annotated_DT_centers,
                                         x = ~center_1,
                                         y = ~center_2,
                                         z = ~center_3,
                                         color = annotated_DT_centers[[cell_color]],
                                         colors = cell_color_code,
                                         marker=list(size = 2,symbol = "x",symbols = "x"),
                                         legendgroup = annotated_DT_centers[[cell_color]],
                                         showlegend = F)
      }
      if(show_center_label == TRUE){
        cat(" center label is not clear to see in 3D plot\n You can shut it down with show_center_label = F\n")
        dpl <- dpl %>% plotly::add_trace(mode = "text",
                                         type = "scatter3d",
                                         data = annotated_DT_centers,
                                         x = ~center_1,
                                         y = ~center_2,
                                         z = ~center_3,
                                         text = annotated_DT_centers[[cell_color]],
                                         legendgroup = annotated_DT_centers[[cell_color]],
                                         showlegend = F)
      }
      
    }
    if(dim_reduction_to_use == 'pca') {
      x_name = paste0('pca','-',dim_names[1])
      y_name = paste0('pca','-',dim_names[2])
      z_name = paste0('pca','-',dim_names[3])
      x_title = sprintf('%s explains %.02f%% of variance', x_name, var_expl_vec[1])
      y_title = sprintf('%s explains %.02f%% of variance', y_name, var_expl_vec[2])
      z_title = sprintf('%s explains %.02f%% of variance', z_name, var_expl_vec[3])
    }
    else{
      x_title = paste(dim_reduction_to_use,dim_names[1],sep = " ")
      y_title = paste(dim_reduction_to_use,dim_names[2],sep = " ")
      z_title = paste(dim_reduction_to_use,dim_names[3],sep = " ")
    }
    
  }
  
  
  
  #********** spatial plot ***********#
  if(is.null(sdimx) | is.null(sdimy)) {
    # cat('first and second dimenion need to be defined, default is first 2 \n')
    sdimx = 'sdimx'
    sdimy = 'sdimy'
  }
  
  # 2D plot
  if(is.null(sdimz)){
    spl <- plotly::plot_ly()
    
    if(show_spatial_network == TRUE) {
      if(is.null(spatial_network)){
        stop("No usable spatial network specified! Please choose a network with spatial_network_name=xxx")
      }
      else{
        if(is.null(spatial_network_alpha)) {
          spatial_network_alpha = 0.5
        }
        else if(is.character(spatial_network_alpha)){
          warning("Edge_alpha for plotly mode is not adjustable yet. Default 0.5 will be set\n")
          spatial_network_alpha = 0.5
        }
        spl <- spl %>% plotly::add_segments(name = spatial_network_name,
                                            type = "scatter",
                                            x = spatial_network[["sdimx_begin"]],
                                            y = spatial_network[["sdimy_begin"]],
                                            xend = spatial_network[["sdimx_end"]],
                                            yend = spatial_network[["sdimy_end"]],
                                            line = list(color = network_color,
                                                        width = 0.5),
                                            opacity=spatial_network_alpha)
      }
    }
    if(show_spatial_grid == TRUE){
      if(is.null(spatial_grid)){
        stop("No usable spatial grid specified! Please choose a network with spatial_grid_name=xxx")
      }
      else{
        if(is.null(spatial_grid_color)) {
          spatial_grid_color = 'black'
        }
        edges <- plotly_grid(spatial_grid)
        spl <- spl %>% plotly::add_segments(name = "spatial_grid",
                                            type = "scatter",
                                            data = edges,
                                            x = ~x,
                                            y = ~y,
                                            xend = ~x_end,
                                            yend = ~y_end,
                                            line = list(color = spatial_grid_color,
                                                        width = 1),
                                            opacity=spatial_grid_alpha)
        
      }
    }
    if(is.null(cell_color)){
      #cell_color = "lightblue"
      spl <- spl %>% plotly::add_trace(type = "scatter",mode = "markers",
                                       x = annotated_DT[[sdimx]],
                                       y = annotated_DT[[sdimy]],
                                       #color = "lightblue",
                                       #colors = "lightblue",
                                       marker = list(size = spatial_point_size,
                                                     color = "lightblue"),
                                       showlegend = F)
    }
    
    else if(cell_color %in% colnames(annotated_DT)){
      if(color_as_factor){
        annotated_DT[[cell_color]] <- as.factor(annotated_DT[[cell_color]])
      }
      
      
      spl <- spl %>% plotly::add_trace(type = "scatter",mode = "markers",
                                       x = annotated_DT[[sdimx]],
                                       y = annotated_DT[[sdimy]],
                                       color = annotated_DT[[cell_color]],
                                       colors = cell_color_code,
                                       legendgroup = annotated_DT[[cell_color]],
                                       marker = list(size = spatial_point_size),
                                       showlegend = F)
    }
    else{
      stop("cell_color doesn't exist!\n")
    }
    spl <- spl %>% plotly::layout(xaxis = list(title = "X"),
                                  yaxis = list(title = "Y"),
                                  legend = list(x = 100, y = 0.5,font = list(family = "sans-serif",size = legend_text_size)))
    
  }
  
  
  # 3D plot
  else{
    axis_scale = match.arg(axis_scale, c("cube","real","custom"))
    
    ratio = plotly_axis_scale_3D(annotated_DT,sdimx = sdimx,sdimy = sdimy,sdimz = sdimz,
                                 mode = axis_scale,custom_ratio = custom_ratio)
    spl <- plotly::plot_ly(scene = "scene2")
    if(!is.null(cell_color)) {
      if(cell_color %in% colnames(annotated_DT)){
        annotated_DT[[cell_color]] <- as.factor(annotated_DT[[cell_color]])
        spl <- spl %>% plotly::add_trace(type = 'scatter3d',mode = 'markers',
                                         x = annotated_DT[[sdimx]],
                                         y = annotated_DT[[sdimy]],
                                         z = annotated_DT[[sdimz]],
                                         color = annotated_DT[[cell_color]],
                                         colors = cell_color_code,
                                         legendgroup = annotated_DT[[cell_color]],
                                         marker = list(size = spatial_point_size),
                                         showlegend = F)
      }
      else{
        stop("cell_color doesn't exist!\n")
      }
    }
    else{
      spl <- spl %>% plotly::add_trace(type = 'scatter3d',mode = 'markers',
                                       x = annotated_DT$sdimx,
                                       y = annotated_DT$sdimy,
                                       z = annotated_DT$sdimz,
                                       color = "lightblue",
                                       colors = "lightblue",
                                       #legendgroup = annotated_DT[[cell_color]],
                                       marker = list(size = spatial_point_size),
                                       showlegend = F)
    }
    if(show_spatial_network == TRUE) {
      if(is.null(spatial_network)){
        stop("No usable spatial network specified! Please choose a network with spatial_network_name=xxx")
      }
      else{
        if(is.null(spatial_network_alpha)) {
          spatial_network_alpha = 0.5
        }
        else if(is.character(spatial_network_alpha)){
          warning("Edge_alpha for plotly mode is not adjustable yet. Default 0.5 will be set\n")
          spatial_network_alpha = 0.5
        }
        edges <- plotly_network(spatial_network)
        
        spl <- spl %>% plotly::add_trace(name = "sptial network",
                                         mode = "lines",
                                         type = "scatter3d",
                                         data = edges,
                                         x = ~x,y=~y,z=~z,
                                         line=list(color=network_color),
                                         opacity=spatial_network_alpha)
      }
    }
    
    if(show_spatial_grid == TRUE){
      cat("3D grid is not clear to view\n")
    }
    
  }
  
  
  
  
  if(is.null(dim3_to_use) & is.null(sdimz)){
    if(plot_alignment == 'vertical'){
      combo_plot <- plotly::subplot(dpl,spl,nrows = 2,titleX = TRUE,titleY = TRUE)
    }
    else{
      combo_plot <- plotly::subplot(dpl,spl,titleX = TRUE,titleY = TRUE)
    }
  }
  
  else if(!is.null(dim3_to_use) & is.null(sdimz)){
    if(plot_alignment == 'vertical'){
      combo_plot <- plotly::subplot(dpl,spl,nrows = 2,titleX = TRUE,titleY = TRUE)%>%
        plotly::layout(scene = list(domain = list(x = c(0, 1), y = c(0,0.5)),
                                    xaxis = list(title = x_title),
                                    yaxis = list(title = y_title),
                                    zaxis = list(title = z_title)))
    }
    else{
      combo_plot <- plotly::subplot(dpl,spl,titleX = TRUE,titleY = TRUE) %>%
        plotly::layout(scene = list(domain = list(x = c(0, 0.5), y = c(0,1)),
                                    xaxis = list(title = x_title),
                                    yaxis = list(title = y_title),
                                    zaxis = list(title = z_title)))
    }
  }
  
  else if(is.null(dim3_to_use) & !is.null(sdimz)){
    if(plot_alignment == 'vertical'){
      combo_plot <- plotly::subplot(dpl,spl,nrows = 2,titleX = TRUE,titleY = TRUE)%>%
        plotly::layout(scene2 = list(domain = list(x = c(0, 1), y = c(0.5,1)),
                                     xaxis = list(title = "X",nticks = x_ticks),
                                     yaxis = list(title = "Y",nticks = y_ticks),
                                     zaxis = list(title = "Z",nticks = z_ticks),
                                     aspectmode='manual',
                                     aspectratio = list(x=ratio[[1]],
                                                        y=ratio[[2]],
                                                        z=ratio[[3]])))
    }
    else{
      combo_plot <- plotly::subplot(dpl,spl,titleX = TRUE,titleY = TRUE) %>%
        plotly::layout(scene2 = list(domain = list(x = c(0.5, 1), y = c(0,1)),
                                     xaxis = list(title = "X",nticks = x_ticks),
                                     yaxis = list(title = "Y",nticks = y_ticks),
                                     zaxis = list(title = "Z",nticks = z_ticks),
                                     aspectmode='manual',
                                     aspectratio = list(x=ratio[[1]],
                                                        y=ratio[[2]],
                                                        z=ratio[[3]])))
    }
  }
  
  else if(!is.null(dim3_to_use) & !is.null(sdimz)){
    if(plot_alignment == 'vertical'){
      combo_plot <- plotly::subplot(dpl,spl,nrows = 2,titleX = TRUE,titleY = TRUE)%>%
        plotly::layout(scene = list(domain = list(x = c(0, 1), y = c(0,0.5)),
                                    xaxis = list(title = x_title),
                                    yaxis = list(title = y_title),
                                    zaxis = list(title = z_title)),
                       scene2 = list(domain = list(x = c(0, 1), y = c(0.5,1)),
                                     xaxis = list(title = "X",nticks = x_ticks),
                                     yaxis = list(title = "Y",nticks = y_ticks),
                                     zaxis = list(title = "Z",nticks = z_ticks),
                                     aspectmode='manual',
                                     aspectratio = list(x=ratio[[1]],
                                                        y=ratio[[2]],
                                                        z=ratio[[3]])))
    }
    else{
      combo_plot <- plotly::subplot(dpl,spl,titleX = TRUE,titleY = TRUE) %>%
        plotly::layout(scene = list(domain = list(x = c(0, 0.5), y = c(0,1)),
                                    xaxis = list(title = x_title),
                                    yaxis = list(title = y_title),
                                    zaxis = list(title = z_title)),
                       scene2 = list(domain = list(x = c(0.5, 1), y = c(0,1)),
                                     xaxis = list(title = "X",nticks = x_ticks),
                                     yaxis = list(title = "Y",nticks = y_ticks),
                                     zaxis = list(title = "Z",nticks = z_ticks),
                                     aspectmode='manual',
                                     aspectratio = list(x=ratio[[1]],
                                                        y=ratio[[2]],
                                                        z=ratio[[3]])))
    }
  }
  
  
  
  return(combo_plot)
  
}




#' @title visSpatDimPlot
#' @name visSpatDimPlot
#' @description integration of visSpatDimPlot_2D and visSpatDimPlot_3D
#' @param gobject giotto object
#' @param plot_mode choose the mode to draw plot : ggplot or plotly
#' @param plot_alignment direction to align plot
#' @param dim_reduction_to_use dimension reduction to use
#' @param dim_reduction_name dimension reduction name
#' @param dim1_to_use dimension to use on x-axis
#' @param dim2_to_use dimension to use on y-axis
#' @param dim3_to_use dimension to use on z-axis
#' @param show_NN_network show underlying NN network
#' @param nn_network_to_use type of NN network to use (kNN vs sNN)
#' @param network_name name of NN network to use, if show_NN_network = TRUE
#' @param cell_color color for cells (see details)
#' @param color_as_factor convert color column to factor
#' @param cell_color_code named vector with colors
#' @param select_cell_groups select subset of cells/clusters based on cell_color parameter
#' @param select_cells select subset of cells based on cell IDs
#' @param show_other_cells display not selected cells
#' @param other_cell_color color of not selected cells
#' @param nn_network_alpha column to use for alpha of the edges
#' @param show_spatial_network show spatial network
#' @param spatial_network_name name of spatial network to use
#' @param spatial_network_color color of spatial network
#' @param spatial_network_alpha alpha of spatial network
#' @param show_spatial_grid show spatial grid
#' @param spatial_grid_name name of spatial grid to use
#' @param spatial_grid_color color of spatial grid
#' @param spatial_grid_alpha alpha of spatial grid
#' @param show_legend show legend
#' @param show_plot show plot
#' @param legend_text_size text size of legend
#' @return ggplot or plotly
#' @details Description of parameters.
#' @export
#' @examples
#'     visSpatDimPlot(gobject)
visSpatDimPlot <- function(gobject,
                           plot_method = c("ggplot", "plotly"),
                           plot_alignment = NULL,
                           dim_reduction_to_use = 'umap',
                           dim_reduction_name = 'umap',
                           dim1_to_use = 1,
                           dim2_to_use = 2,
                           dim3_to_use = NULL,
                           sdimx=NULL,
                           sdimy=NULL,
                           sdimz=NULL,
                           show_NN_network = F,
                           nn_network_to_use = 'sNN',
                           network_name = 'sNN.pca',
                           show_cluster_center = F,
                           show_center_label = T,
                           center_point_size = 4,
                           label_size = NULL,          #ggplot 4, plotly 16
                           label_fontface = 'bold',
                           cell_color = NULL,
                           color_as_factor = T,
                           cell_color_code = NULL,
                           select_cell_groups = NULL,
                           select_cells = NULL,
                           show_other_cells = T,
                           other_cell_color = 'lightgrey',
                           dim_point_size = 3,
                           dim_point_border_col = 'black',
                           dim_point_border_stroke = 0.1,
                           nn_network_alpha = NULL,    #ggplot 0.05, plotly 0.5
                           show_spatial_network = F,
                           spatial_network_name = 'spatial_network',
                           network_color = "lightgray",
                           spatial_network_alpha = 0.5,
                           show_spatial_grid = F,
                           spatial_grid_name = 'spatial_grid',
                           spatial_grid_color = NULL,
                           spatial_grid_alpha = 0.5,
                           spatial_point_size = 3,
                           legend_text_size = 12,
                           spatial_point_border_col = 'black',
                           spatial_point_border_stroke = 0.1,
                           show_legend = T,
                           axis_scale = c("cube","real","custom"),
                           custom_ratio = NULL,
                           x_ticks = NULL,
                           y_ticks = NULL,
                           z_ticks = NULL,
                           show_plot = F){
  
  plot_method = match.arg(plot_method, c("ggplot", "plotly"))
  
  axis_scale = match.arg(axis_scale,c("cube","real","custom"))
  if(plot_method == "ggplot"){
    if(!is.null(sdimz) | !is.null(dim3_to_use)){
      stop("3D plots can not be drawn in ggplot2 mode, please choose plotly mode with plot_mode = \"plotly\"\n")
    }
    if(is.null(plot_alignment)){
      plot_alignment = "vertical"
    }
    if(is.null(label_size)){
      label_size = 4
    }
    if(is.null(nn_network_alpha)){
      nn_network_alpha = 0.05
    }
    pl <- visSpatDimPlot_2D(gobject = gobject,
                            plot_alignment = plot_alignment,
                            dim_reduction_to_use = dim_reduction_to_use,
                            dim_reduction_name =dim_reduction_name,
                            dim1_to_use = dim1_to_use,
                            dim2_to_use = dim2_to_use,
                            sdimx=sdimx,
                            sdimy=sdimy,
                            show_NN_network = show_NN_network,
                            nn_network_to_use = nn_network_to_use,
                            network_name =network_name,
                            show_cluster_center = show_cluster_center,
                            show_center_label = show_center_label,
                            center_point_size = center_point_size,
                            label_size = label_size,
                            label_fontface = label_fontface,
                            cell_color = cell_color,
                            color_as_factor = color_as_factor,
                            cell_color_code = cell_color_code,
                            select_cell_groups = select_cell_groups,
                            select_cells = select_cells,
                            show_other_cells = show_other_cells,
                            other_cell_color = other_cell_color,
                            dim_plot_mode = dim_plot_mode,
                            dim_point_size = dim_point_size,
                            dim_point_border_col =dim_point_border_col,
                            dim_point_border_stroke = dim_point_border_stroke,
                            nn_network_alpha = nn_network_alpha,
                            show_spatial_network = show_spatial_network,
                            spatial_network_name = spatial_network_name,
                            spatial_network_color = network_color,
                            show_spatial_grid = show_spatial_grid,
                            spatial_grid_name = spatial_grid_name,
                            spatial_grid_color = spatial_grid_color,
                            spatial_point_size = spatial_point_size,
                            spatial_point_border_col = spatial_point_border_col,
                            spatial_point_border_stroke = spatial_point_border_stroke,
                            show_legend = show_legend,
                            show_plot = show_plot,
                            plot_method = plot_method)
  }
  
  else if(plot_method == "plotly"){
    if(is.null(plot_alignment)){
      plot_alignment = "horizontal"
    }
    if(is.null(label_size)){
      label_size = 16
    }
    if(is.null(nn_network_alpha)){
      nn_network_alpha = 0.5
    }
    pl <- visSpatDimPlot_3D(gobject = gobject,
                            plot_alignment = plot_alignment,
                            dim_reduction_to_use = dim_reduction_to_use,
                            dim_reduction_name = dim_reduction_name,
                            dim1_to_use = dim1_to_use,
                            dim2_to_use = dim2_to_use,
                            dim3_to_use = dim3_to_use,
                            sdimx=sdimx,
                            sdimy=sdimy,
                            sdimz=sdimz,
                            show_NN_network = show_NN_network,
                            nn_network_to_use = nn_network_to_use,
                            network_name = network_name,
                            show_cluster_center = show_cluster_center,
                            show_center_label = show_center_label,
                            center_point_size = center_point_size,
                            label_size = label_size,
                            cell_color = cell_color,
                            color_as_factor = color_as_factor,
                            cell_color_code = cell_color_code,
                            dim_point_size = dim_point_size,
                            nn_network_alpha = nn_network_alpha,
                            show_spatial_network = show_spatial_network,
                            spatial_network_name = spatial_network_name,
                            network_color = network_color,
                            spatial_network_alpha = spatial_network_alpha,
                            show_spatial_grid = show_spatial_grid,
                            spatial_grid_name = spatial_grid_name,
                            spatial_grid_color = spatial_grid_color,
                            spatial_grid_alpha = spatial_grid_alpha,
                            spatial_point_size = spatial_point_size,
                            axis_scale = axis_scale,
                            custom_ratio = custom_ratio,
                            x_ticks = x_ticks,
                            y_ticks = y_ticks,
                            z_ticks = z_ticks,
                            legend_text_size = legend_text_size)
  }
  return(pl)
}





#' @title visDimGenePlot_2D_ggplot
#' @name visDimGenePlot_2D_ggplot
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
#' @param point_border_col color of border around points
#' @param point_border_stroke stroke size of border around points
#' @param midpoint size of point (cell)
#' @param cow_n_col cowplot param: how many columns
#' @param cow_rel_h cowplot param: relative height
#' @param cow_rel_w cowplot param: relative width
#' @param cow_align cowplot param: how to align
#' @param show_legend show legend
#' @param show_plots show plots
#' @return ggplot
#' @details Description of parameters.
#' @export
#' @examples
#'     visDimGenePlot_2D_ggplot(gobject)
visDimGenePlot_2D_ggplot <- function(gobject,
                                     expression_values = c('normalized', 'scaled', 'custom'),
                                     genes = NULL,
                                     dim_reduction_to_use = 'umap',
                                     dim_reduction_name = 'umap',
                                     dim1_to_use = 1,
                                     dim2_to_use = 2,
                                     show_NN_network = F,
                                     nn_network_to_use = 'sNN',
                                     network_name = 'sNN.pca',
                                     network_color = "lightgray",
                                     edge_alpha = NULL,
                                     scale_alpha_with_expression = FALSE,
                                     point_size = 1,
                                     genes_high_color = "red",
                                     genes_mid_color = "white",
                                     genes_low_color = "blue",
                                     point_border_col = 'black',
                                     point_border_stroke = 0.1,
                                     midpoint = 0,
                                     cow_n_col = 2,
                                     cow_rel_h = 1,
                                     cow_rel_w = 1,
                                     cow_align = 'h',
                                     show_legend = T,
                                     show_plots = F){
  ## select genes ##
  selected_genes = genes
  values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))
  expr_values = Giotto:::select_expression_values(gobject = gobject, values = values)
  
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
  
  
  
  ## visualize multipe plots ##
  ## 2D plots ##
  savelist <- list()
  
  for(gene in selected_genes) {
    
    
    ## OLD need to be combined ##
    pl <- ggplot2::ggplot()
    pl <- pl + ggplot2::theme_classic()
    
    # network layer
    if(show_NN_network == TRUE) {
      
      if(is.null(edge_alpha)) {
        edge_alpha = 0.5
        pl <- pl + ggplot2::geom_segment(data = annotated_network_DT,
                                         aes_string(x = from_dim_names[1], y = from_dim_names[2],
                                                    xend = to_dim_names[1], yend = to_dim_names[2]),
                                         alpha = edge_alpha, color=network_color,size = 0.1,
                                         show.legend = F)
      } else if(is.numeric(edge_alpha)) {
        pl <- pl + ggplot2::geom_segment(data = annotated_network_DT,
                                         aes_string(x = from_dim_names[1], y = from_dim_names[2],
                                                    xend = to_dim_names[1], yend = to_dim_names[2]),
                                         alpha = edge_alpha, color=network_color,size = 0.1,
                                         show.legend = F)
      } else if(is.character(edge_alpha)) {
        
        if(edge_alpha %in% colnames(annotated_network_DT)) {
          pl <- pl + ggplot2::geom_segment(data = annotated_network_DT,
                                           aes_string(x = from_dim_names[1], y = from_dim_names[2],
                                                      xend = to_dim_names[1],
                                                      yend = to_dim_names[2], alpha = edge_alpha),
                                           color=network_color,
                                           show.legend = F)
        }
      }
    }
    
    
    # point layer
    if(is.null(genes)) {
      cell_color = 'lightblue'
      pl <- pl + ggplot2::geom_point(data = annotated_gene_DT, aes_string(x = dim_names[1], dim_names[2]),
                                     fill = cell_color, show.legend = show_legend, size =  point_size)
      
    } else {
      if(scale_alpha_with_expression == TRUE) {
        pl <- pl + ggplot2::geom_point(data = annotated_gene_DT, aes_string(x = dim_names[1], y = dim_names[2], fill = gene, alpha = gene),
                                       show.legend = show_legend, shape = 21, size = point_size,
                                       color = point_border_col, stroke = point_border_stroke)
      } else {
        pl <- pl + ggplot2::geom_point(data = annotated_gene_DT, aes_string(x = dim_names[1], y = dim_names[2], fill = gene),
                                       show.legend = show_legend, shape = 21,
                                       size =  point_size,
                                       color = point_border_col, stroke = point_border_stroke)
      }
      
      pl <- pl + ggplot2::scale_fill_gradient2(low = genes_low_color, mid = genes_mid_color, high = genes_high_color, midpoint = midpoint)
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


#' @title visDimGenePlot_3D_plotly
#' @name visDimGenePlot_3D_plotly
#' @description Visualize cells and gene expression according to dimension reduction coordinates
#' @param gobject giotto object
#' @param expression_values gene expression values to use
#' @param genes genes to show
#' @param dim_reduction_to_use dimension reduction to use
#' @param dim_reduction_name dimension reduction name
#' @param dim1_to_use dimension to use on x-axis
#' @param dim2_to_use dimension to use on y-axis
#' @param dim3_to_use dimension to use on z-axis
#' @param show_NN_network show underlying NN network
#' @param nn_network_to_use type of NN network to use (kNN vs sNN)
#' @param network_name name of NN network to use, if show_NN_network = TRUE
#' @param edge_alpha column to use for alpha of the edges
#' @param point_size size of point (cell)
#' @param show_legend show legend
#' @param show_plots show plots
#' @return ggplot
#' @details Description of parameters.
#' @export
#' @examples
#'     visDimGenePlot_3D_plotly(gobject)


visDimGenePlot_3D_plotly <- function(gobject,
                                     expression_values = c('normalized', 'scaled', 'custom'),
                                     genes = NULL,
                                     dim_reduction_to_use = 'umap',
                                     dim_reduction_name = 'umap',
                                     dim1_to_use = 1,
                                     dim2_to_use = 2,
                                     dim3_to_use = 3,
                                     show_NN_network = F,
                                     nn_network_to_use = 'sNN',
                                     network_name = 'sNN.pca',
                                     network_color = "lightgray",
                                     edge_alpha = NULL,
                                     point_size = 1,
                                     genes_high_color = NULL,
                                     genes_mid_color = "white",
                                     genes_low_color = "blue",
                                     show_legend = T,
                                     show_plots = F){
  ## select genes ##
  selected_genes = genes
  values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))
  expr_values = Giotto:::select_expression_values(gobject = gobject, values = values)
  
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
  dim_dfr = gobject@dimension_reduction$cells[[dim_reduction_to_use]][[dim_reduction_name]]$coordinates[,c(dim1_to_use, dim2_to_use,dim3_to_use)]
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
  
  
  
  ## visualize multipe plots ##
  ## 3D plots ##
  
  
  if(show_NN_network == TRUE){
    edges <- plotly_network(annotated_network_DT,
                            "from_Dim.1","from_Dim.2","from_Dim.3",
                            "to_Dim.1","to_Dim.2","to_Dim.3")
  }
  ##Point layer
  if(length(selected_genes) > 4){
    stop("\n The max number of genes showed together is 4.Otherwise it will be too small to see\n
         \n If you have more genes to show, please divide them into groups\n")
  }
  if(!is.null(genes_high_color)){
    if(length(genes_high_color)!=length(selected_genes)&length(genes_high_color) != 1){
      stop('\n The number of genes and their corresbonding do not match\n')
    }
  }
  else if (is.null(genes_high_color)){
    genes_high_color = rep("red",length(selected_genes))
  }
  else{
    genes_high_color = rep(genes_high_color,length(selected_genes))
  }
  
  titleX = title = paste(dim_reduction_to_use,dim_names[1],sep = " ")
  titleY = title = paste(dim_reduction_to_use,dim_names[2],sep = " ")
  titleZ = title = paste(dim_reduction_to_use,dim_names[3],sep = " ")
  savelist <- list()
  for(i in 1:length(selected_genes)){
    
    gene = selected_genes[i]
    
    pl <- plotly::plot_ly(name = gene,scene=paste("scene",i,sep = ""))
    pl <- pl %>%  plotly::add_trace(data = annotated_gene_DT,type = 'scatter3d',mode = "markers",
                                    x = annotated_gene_DT[[dim_names[1]]],
                                    y = annotated_gene_DT[[dim_names[2]]],
                                    z = annotated_gene_DT[[dim_names[3]]],
                                    color = annotated_gene_DT[[gene]],
                                    colors = c(genes_low_color,genes_mid_color,genes_high_color[i]),
                                    marker = list(size = point_size))
    
    ## plot spatial network
    if(show_NN_network == TRUE) {
      pl <- pl %>% plotly::add_trace(name = "sptial network",mode = "lines",
                                     type = "scatter3d",opacity = edge_alpha,
                                     showlegend = F,
                                     data = edges,
                                     x = ~x,y=~y,z=~z,
                                     line=list(color=network_color, width = 0.5))
    }
    pl <- pl %>% plotly::colorbar(title = gene)
    savelist[[gene]] <- pl
  }
  
  if(length(savelist) == 1){
    savelist[[1]] <- savelist[[1]] %>% plotly::layout(scene = list(
      xaxis = list(title = titleX),
      yaxis = list(title = titleY),
      zaxis = list(title = titleZ)))
    if(show_plots){
      print(savelist[[1]])
    }
    return (savelist[[1]])
  }
  else if(length(savelist)==2){
    cowplot <- suppressWarnings(plotly::subplot(savelist,titleX = TRUE,titleY = TRUE)%>%
                                  plotly::layout(scene = list(domain = list(x = c(0, 0.5), y = c(0,1)),
                                                              xaxis = list(title = titleX),
                                                              yaxis = list(title = titleY),
                                                              zaxis = list(title = titleZ)),
                                                 scene2 = list(domain = list(x = c(0.5, 1), y = c(0,1)),
                                                               xaxis = list(title = titleX),
                                                               yaxis = list(title = titleY),
                                                               zaxis = list(title = titleZ)),
                                                 legend = list(x = 100, y = 0)))
  }
  else if(length(savelist)==3){
    cowplot <- suppressWarnings(plotly::subplot(savelist,titleX = TRUE,titleY = TRUE)%>%
                                  plotly::layout(scene = list(domain = list(x = c(0, 0.5), y = c(0,0.5)),
                                                              xaxis = list(title = titleX),
                                                              yaxis = list(title = titleY),
                                                              zaxis = list(title = titleZ)),
                                                 scene2 = list(domain = list(x = c(0.5, 1), y = c(0,0.5)),
                                                               xaxis = list(title = titleX),
                                                               yaxis = list(title = titleY),
                                                               zaxis = list(title = titleZ)),
                                                 scene3 = list(domain = list(x = c(0, 0.5), y = c(0.5,1)),
                                                               xaxis = list(title = titleX),
                                                               yaxis = list(title = titleY),
                                                               zaxis = list(title = titleZ)),
                                                 legend = list(x = 100, y = 0)))
  }
  else if(length(savelist)==4){
    
    cowplot <- suppressWarnings(plotly::subplot(savelist)%>% plotly::layout(scene = list(domain = list(x = c(0, 0.5), y = c(0,0.5)),
                                                                                         xaxis = list(title = titleX),
                                                                                         yaxis = list(title = titleY),
                                                                                         zaxis = list(title = titleZ)),
                                                                            scene2 = list(domain = list(x = c(0.5, 1), y = c(0,0.5)),
                                                                                          xaxis = list(title = titleX),
                                                                                          yaxis = list(title = titleY),
                                                                                          zaxis = list(title = titleZ)),
                                                                            scene3 = list(domain = list(x = c(0, 0.5), y = c(0.5,1)),
                                                                                          xaxis = list(title = titleX),
                                                                                          yaxis = list(title = titleY),
                                                                                          zaxis = list(title = titleZ)),
                                                                            scene4 = list(domain = list(x = c(0.5, 1), y = c(0.5,1)),
                                                                                          xaxis = list(title = titleX),
                                                                                          yaxis = list(title = titleY),
                                                                                          zaxis = list(title = titleZ)),
                                                                            legend = list(x = 100, y = 0)))
  }
  return(cowplot)
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
#' @param dim3_to_use dimension to use on z-axis
#' @param show_NN_network show underlying NN network
#' @param nn_network_to_use type of NN network to use (kNN vs sNN)
#' @param network_name name of NN network to use, if show_NN_network = TRUE
#' @param edge_alpha column to use for alpha of the edges
#' @param scale_alpha_with_expression scale expression with ggplot alpha parameter
#' @param point_size size of point (cell)
#' @param point_border_col color of border around points
#' @param point_border_stroke stroke size of border around points
#' @param midpoint size of point (cell)
#' @param cow_n_col cowplot param: how many columns
#' @param cow_rel_h cowplot param: relative height
#' @param cow_rel_w cowplot param: relative width
#' @param cow_align cowplot param: how to align
#' @param show_legend show legend
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
                           dim3_to_use = NULL,
                           show_NN_network = F,
                           nn_network_to_use = 'sNN',
                           network_name = 'sNN.pca',
                           network_color = "lightgray",
                           edge_alpha = NULL,
                           scale_alpha_with_expression = FALSE,
                           point_size = 1,
                           genes_high_color = NULL,
                           genes_mid_color = "white",
                           genes_low_color = "blue",
                           point_border_col = 'black',
                           point_border_stroke = 0.1,
                           midpoint = 0,
                           cow_n_col = 2,
                           cow_rel_h = 1,
                           cow_rel_w = 1,
                           cow_align = 'h',
                           show_legend = T,
                           plot_method = c('ggplot','plotly'),
                           show_plots = F){
  
  plot_method = match.arg(plot_method, choices = c('ggplot','plotly'))
  
  
  if(plot_method == 'ggplot'){
    if(!is.null(dim3_to_use)){
      warning("\n ggplot is only in 2D. If you want to plot the third dim, please choose plotly method with \"plot_method = \"plotly\"\"\n")
    }
    if(is.null(genes_high_color)){
      genes_high_color = "red"
    }
    result <- visDimGenePlot_2D_ggplot(gobject = gobject,
                                       expression_values =expression_values ,
                                       genes = genes,
                                       dim_reduction_to_use = dim_reduction_to_use,
                                       dim_reduction_name = dim_reduction_name,
                                       dim1_to_use = dim1_to_use,
                                       dim2_to_use = dim2_to_use,
                                       show_NN_network = show_NN_network,
                                       nn_network_to_use = nn_network_to_use,
                                       network_name = network_name,
                                       network_color = network_color,
                                       edge_alpha = edge_alpha,
                                       scale_alpha_with_expression = scale_alpha_with_expression,
                                       point_size = point_size,
                                       genes_high_color = genes_high_color,
                                       genes_mid_color = genes_mid_color,
                                       genes_low_color = genes_low_color,
                                       point_border_col = point_border_col,
                                       point_border_stroke = point_border_stroke,
                                       midpoint = midpoint,
                                       cow_n_col = cow_n_col,
                                       cow_rel_h = cow_rel_h,
                                       cow_rel_w = cow_rel_w,
                                       cow_align = cow_align,
                                       show_legend = show_legend,
                                       show_plots = show_plots)
    
  }
  else{
    if(is.null(dim3_to_use)){
      stop("\n plotly is in 3D and you need to define dim3_to_use \n")
    }
    result <- visDimGenePlot_3D_plotly(gobject =gobject,
                                       expression_values = expression_values,
                                       genes = genes,
                                       dim_reduction_to_use = dim_reduction_to_use,
                                       dim_reduction_name = dim_reduction_name,
                                       dim1_to_use = dim1_to_use,
                                       dim2_to_use = dim2_to_use,
                                       dim3_to_use = dim3_to_use,
                                       show_NN_network = show_NN_network,
                                       nn_network_to_use = nn_network_to_use,
                                       network_name = network_name,
                                       network_color = network_color,
                                       edge_alpha = edge_alpha,
                                       point_size = point_size,
                                       genes_high_color = genes_high_color,
                                       genes_mid_color = genes_mid_color,
                                       genes_low_color = genes_low_color,
                                       show_legend = show_legend,
                                       show_plots = show_plots)
  }
  return(result)
}



#' @title visSpatDimGenePlot_2D
#' @name visSpatDimGenePlot_2D
#' @description Visualize cells according to spatial AND dimension reduction coordinates in ggplot mode
#' @param gobject giotto object
#' @param expression_values gene expression values to use
#' @param plot_alignment direction to align plot
#' @param genes genes to show
#' @param dim_reduction_to_use dimension reduction to use
#' @param dim_reduction_name dimension reduction name
#' @param dim1_to_use dimension to use on x-axis
#' @param dim2_to_use dimension to use on y-axis
#' @param dim_point_size dim reduction plot: point size
#' @param dim_point_border_col color of border around points
#' @param dim_point_border_stroke stroke size of border around points
#' @param show_NN_network show underlying NN network
#' @param nn_network_to_use type of NN network to use (kNN vs sNN)
#' @param network_name name of NN network to use, if show_NN_network = TRUE
#' @param edge_alpha_dim dim reduction plot: column to use for alpha of the edges
#' @param scale_alpha_with_expression scale expression with ggplot alpha parameter
#' @param spatial_network_name name of spatial network to use
#' @param spatial_grid_name name of spatial grid to use
#' @param spatial_point_size spatial plot: point size
#' @param spatial_point_border_col color of border around points
#' @param spatial_point_border_stroke stroke size of border around points
#' @param midpoint size of point (cell)
#' @param point_size size of point (cell)
#' @param cow_n_col cowplot param: how many columns
#' @param cow_rel_h cowplot param: relative height
#' @param cow_rel_w cowplot param: relative width
#' @param cow_align cowplot param: how to align
#' @param show_legend show legend
#' @param show_plot show plot
#' @return ggplot
#' @details Description of parameters.
#' @export
#' @examples
#'     visSpatDimGenePlot_2D(gobject)
visSpatDimGenePlot_2D <- function(gobject,
                                  expression_values = c('normalized', 'scaled', 'custom'),
                                  plot_alignment = c('horizontal', 'vertical'),
                                  genes,
                                  dim_reduction_to_use = 'umap',
                                  dim_reduction_name = 'umap',
                                  dim1_to_use = 1,
                                  dim2_to_use = 2,
                                  point_size = 1,
                                  dim_point_border_col = 'black',
                                  dim_point_border_stroke = 0.1,
                                  show_NN_network = F,
                                  show_spatial_network = F,
                                  show_spatial_grid = F,
                                  nn_network_to_use = 'sNN',
                                  network_name = 'sNN.pca',
                                  edge_alpha_dim = NULL,
                                  scale_alpha_with_expression = FALSE,
                                  spatial_network_name = 'spatial_network',
                                  spatial_grid_name = 'spatial_grid',
                                  spatial_point_size = 1,
                                  spatial_point_border_col = 'black',
                                  spatial_point_border_stroke = 0.1,
                                  midpoint = 0,
                                  genes_high_color = "red",
                                  genes_mid_color="white",
                                  genes_low_color = "blue",
                                  cow_n_col = 2,
                                  cow_rel_h = 1,
                                  cow_rel_w = 1,
                                  cow_align = 'h',
                                  axis_scale = c("cube","real","custom"),
                                  custom_ratio = NULL,
                                  x_ticks = NULL,
                                  y_ticks = NULL,
                                  show_legend = T,
                                  show_plots = F) {
  
  {
    
    plot_alignment = match.arg(plot_alignment, choices = c('vertical', 'horizontal'))
    
    # dimension reduction plot
    dmpl = visDimGenePlot(gobject = gobject,
                          expression_values = expression_values,
                          genes = genes,
                          dim_reduction_to_use = dim_reduction_to_use,
                          dim_reduction_name = dim_reduction_name,
                          dim1_to_use = dim1_to_use,
                          dim2_to_use = dim2_to_use,
                          show_NN_network = show_NN_network,
                          nn_network_to_use = nn_network_to_use,
                          network_name =network_name,
                          network_color = network_color,
                          edge_alpha = edge_alpha,
                          scale_alpha_with_expression = scale_alpha_with_expression,
                          point_size = point_size,
                          genes_high_color = genes_high_color,
                          genes_mid_color=genes_mid_color,
                          genes_low_color = genes_low_color,
                          point_border_col =dim_point_border_col,
                          point_border_stroke =dim_point_border_stroke,
                          midpoint = midpoint,
                          cow_n_col = cow_n_col,
                          cow_rel_h = cow_rel_h,
                          cow_rel_w = cow_rel_w,
                          cow_align = cow_align,
                          show_legend = show_legend,
                          plot_method = 'ggplot',
                          show_plots = show_plots)
    
    # spatial plot
    spl = visGenePlot(gobject=gobject,
                      expression_values = expression_values,
                      genes = genes,
                      genes_high_color = genes_high_color,
                      genes_mid_color = genes_mid_color,
                      genes_low_color = genes_low_color,
                      show_network = show_spatial_network,
                      network_color = network_color,
                      spatial_network_name = spatial_network_name,
                      edge_alpha = edge_alpha,
                      show_grid = show_spatial_grid,
                      grid_color = grid_color,
                      spatial_grid_name = spatial_grid_name,
                      midpoint = midpoint,
                      scale_alpha_with_expression = scale_alpha_with_expression,
                      point_size = point_size,
                      point_border_col =spatial_point_border_col,
                      point_border_stroke = spatial_point_border_stroke,
                      show_legend = show_legend,
                      cow_n_col = cow_n_col,
                      cow_rel_h = cow_rel_h,
                      cow_rel_w = cow_rel_w,
                      cow_align = cow_align,
                      axis_scale = axis_scale,
                      custom_ratio = custom_ratio,
                      x_ticks = x_ticks,
                      y_ticks = y_ticks,
                      plot_method = 'ggplot',
                      show_plots = F)
    
    if(plot_alignment == 'vertical') {
      combo_plot <- cowplot::plot_grid(dmpl, spl, ncol = 1, rel_heights = c(1), rel_widths = c(1), align = 'v')
      return(cowplot::plot_grid(combo_plot))
    } else {
      combo_plot <- cowplot::plot_grid(dmpl, spl, ncol = 2, rel_heights = c(1), rel_widths = c(1), align = 'h')
      return(cowplot::plot_grid(combo_plot))
    }
  }
  
}

#' @title visSpatDimGenePlot_3D
#' @name visSpatDimGenePlot_3D
#' @description Visualize cells according to spatial AND dimension reduction coordinates in plotly mode
#' @param gobject giotto object
#' @param plot_alignment direction to align plot
#' @param dim_reduction_to_use dimension reduction to use
#' @param dim_reduction_name dimension reduction name
#' @param dim1_to_use dimension to use on x-axis
#' @param dim2_to_use dimension to use on y-axis
#' @param dim3_to_use dimension to use on z-axis
#' @param show_NN_network show underlying NN network
#' @param nn_network_to_use type of NN network to use (kNN vs sNN)
#' @param network_name name of NN network to use, if show_NN_network = TRUE
#' @param genes_high_color color represent high gene expression (see details)
#' @param genes_low_color color represent high gene expression (see details)
#' @param nn_network_alpha column to use for alpha of the edges
#' @param show_spatial_network show spatial network
#' @param spatial_network_name name of spatial network to use
#' @param network_color color of spatial/nn network
#' @param spatial_network_alpha alpha of spatial network
#' @param show_spatial_grid show spatial grid
#' @param spatial_grid_name name of spatial grid to use
#' @param spatial_grid_color color of spatial grid
#' @param spatial_grid_alpha alpha of spatial grid
#' @param show_legend show legend
#' @param show_plot show plot
#' @param legend_text_size text size of legend
#' @return plotly
#' @details Description of parameters.
#' @export
#' @examples
#'     visSpatDimPlot_3D(gobject)
visSpatDimGenePlot_3D <- function(gobject,
                                  expression_values = c('normalized', 'scaled', 'custom'),
                                  plot_alignment = c('horizontal','vertical'),
                                  dim_reduction_to_use = 'umap',
                                  dim_reduction_name = 'umap',
                                  dim1_to_use = 1,
                                  dim2_to_use = 2,
                                  dim3_to_use = NULL,
                                  sdimx=NULL,
                                  sdimy=NULL,
                                  sdimz=NULL,
                                  genes,
                                  show_NN_network = F,
                                  nn_network_to_use = 'sNN',
                                  network_name = 'sNN.pca',
                                  label_size = 16,
                                  genes_low_color = "blue",
                                  genes_mid_color = "white",
                                  genes_high_color = "red",
                                  dim_point_size = 3,
                                  nn_network_alpha = 0.5,
                                  show_spatial_network = F,
                                  spatial_network_name = 'spatial_network',
                                  network_color = "lightgray",
                                  spatial_network_alpha = 0.5,
                                  show_spatial_grid = F,
                                  spatial_grid_name = 'spatial_grid',
                                  spatial_grid_color = NULL,
                                  spatial_grid_alpha = 0.5,
                                  spatial_point_size = 3,
                                  legend_text_size = 12,
                                  axis_scale = c("cube","real","custom"),
                                  custom_ratio = NULL,
                                  x_ticks = NULL,
                                  y_ticks = NULL,
                                  z_ticks = NULL){
  
  plot_alignment = match.arg(plot_alignment, choices = c( 'horizontal','vertical'))
  
  # ********data prepare********#
  ## select genes ##
  selected_genes = genes
  values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))
  expr_values = Giotto:::select_expression_values(gobject = gobject, values = values)
  
  # only keep genes that are in the dataset
  selected_genes = selected_genes[selected_genes %in% rownames(expr_values) ]
  subset_expr_data = expr_values[rownames(expr_values) %in% selected_genes, ]
  t_sub_expr_data_DT = data.table('selected_gene' = subset_expr_data, 'cell_ID' = colnames(expr_values))
  data.table::setnames(t_sub_expr_data_DT, 'selected_gene', selected_genes)
  
  
  ## dimension reduction ##
  dim_dfr = gobject@dimension_reduction$cells[[dim_reduction_to_use]][[dim_reduction_name]]$coordinates[,c(dim1_to_use, dim2_to_use,dim3_to_use)]
  dim_names = colnames(dim_dfr)
  dim_DT = data.table::as.data.table(dim_dfr); dim_DT[, cell_ID := rownames(dim_dfr)]
  
  ## annotated cell metadata
  cell_metadata = gobject@cell_metadata
  annotated_DT = merge(cell_metadata, dim_DT, by = 'cell_ID')
  annotated_DT = merge(annotated_DT, gobject@spatial_locs,by = 'cell_ID')
  annotated_DT = merge(annotated_DT, t_sub_expr_data_DT,by = 'cell_ID')
  
  
  ## nn network
  if(show_NN_network){
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
  
  
  ## extract spatial network
  if(!is.null(spatial_network_name)) {
    spatial_network = select_spatialNetwork(gobject, name = spatial_network_name, return_network_Obj = FALSE)
  } else {
    spatial_network = NULL
  }
  
  ## extract spatial grid
  if(!is.null(spatial_grid_name)) {
    spatial_grid = gobject@spatial_grid[[spatial_grid_name]]
  } else {
    spatial_grid = NULL
  }
  
  
  
  #********** dim plot ***********#
  #2D plot
  if(is.null(dim3_to_use)){
    dpl <- plotly::plot_ly()
    if(show_NN_network == TRUE) {
      if(is.null(nn_network_alpha)) {
        nn_network_alpha = 0.5
      }
      else if(is.character(nn_network_alpha)){
        warning("Edge_alpha for plotly mode is not adjustable yet. Default 0.5 will be set\n")
        nn_network_alpha = 0.5
      }
      dpl <- dpl %>% plotly::add_segments(name = network_name,
                                          type = "scatter",
                                          x = annotated_network_DT[[from_dim_names[1]]],
                                          y = annotated_network_DT[[from_dim_names[2]]],
                                          xend = annotated_network_DT[[to_dim_names[1]]],
                                          yend = annotated_network_DT[[to_dim_names[2]]],
                                          line = list(color = network_color,
                                                      width = 0.5),
                                          opacity=nn_network_alpha)
    }
    
    dpl <- dpl %>% plotly::add_trace(type = "scatter",mode = "markers",
                                     x = annotated_DT[[dim_names[1]]],
                                     y = annotated_DT[[dim_names[2]]],
                                     color = annotated_DT[[selected_genes]],
                                     colors =c(genes_low_color,genes_mid_color,genes_high_color),
                                     marker = list(size = dim_point_size),
                                     showlegend = F)
    
    
    x_title = paste(dim_reduction_to_use,dim_names[1],sep = " ")
    y_title = paste(dim_reduction_to_use,dim_names[2],sep = " ")
    
    dpl <- dpl %>% plotly::layout(xaxis = list(title = x_title),
                                  yaxis = list(title = y_title),
                                  legend = list(x = 100, y = 0.5,font = list(family = "sans-serif",size = legend_text_size)))
  }
  #3D plot
  else if(!is.null(dim3_to_use)){
    dpl <- plotly::plot_ly(scene = "scene1")
    
    dpl <- dpl %>% plotly::add_trace(type = 'scatter3d',mode = "markers",
                                     x = annotated_DT[[dim_names[1]]],
                                     y = annotated_DT[[dim_names[2]]],
                                     z = annotated_DT[[dim_names[3]]],
                                     color = annotated_DT[[selected_genes]],
                                     colors = c(genes_low_color,genes_mid_color,genes_high_color),
                                     marker = list(size = dim_point_size),
                                     showlegend = F)
    #legendgroup = annotated_DT[[cell_color]])
    
    
    if(show_NN_network){
      edges <- plotly_network(annotated_network_DT,
                              "from_Dim.1","from_Dim.2","from_Dim.3",
                              "to_Dim.1","to_Dim.2","to_Dim.3")
      if(is.null(nn_network_alpha)){
        nn_network_alpha = 0.5
      }
      else if(is.character(nn_network_alpha)){
        warning("Edge_alpha for plotly mode is not adjustable yet. Default 0.5 will be set\n")
        nn_network_alpha = 0.5
      }
      
      dpl <- dpl %>% plotly::add_trace(name = network_name,
                                       mode = "lines",
                                       type = "scatter3d",
                                       data = edges,
                                       x = ~x,y=~y,z=~z,
                                       line=list(color=network_color),
                                       opacity=nn_network_alpha)
    }
    
    
    x_title = paste(dim_reduction_to_use,dim_names[1],sep = " ")
    y_title = paste(dim_reduction_to_use,dim_names[2],sep = " ")
    z_title = paste(dim_reduction_to_use,dim_names[3],sep = " ")
    
  }
  dpl <- dpl %>% plotly::colorbar(title = selected_genes)
  
  
  #********** spatial plot ***********#
  if(is.null(sdimx) | is.null(sdimy)) {
    # cat('first and second dimenion need to be defined, default is first 2 \n')
    sdimx = 'sdimx'
    sdimy = 'sdimy'
  }
  
  # 2D plot
  if(is.null(sdimz)){
    spl <- plotly::plot_ly()
    
    if(show_spatial_network == TRUE) {
      if(is.null(spatial_network)){
        stop("No usable spatial network specified! Please choose a network with spatial_network_name=xxx")
      }
      else{
        if(is.null(spatial_network_alpha)) {
          spatial_network_alpha = 0.5
        }
        else if(is.character(spatial_network_alpha)){
          warning("Edge_alpha for plotly mode is not adjustable yet. Default 0.5 will be set\n")
          spatial_network_alpha = 0.5
        }
        spl <- spl %>% plotly::add_segments(name = spatial_network_name,
                                            type = "scatter",
                                            x = spatial_network[["sdimx_begin"]],
                                            y = spatial_network[["sdimy_begin"]],
                                            xend = spatial_network[["sdimx_end"]],
                                            yend = spatial_network[["sdimy_end"]],
                                            line = list(color = network_color,
                                                        width = 0.5),
                                            opacity=spatial_network_alpha)
      }
    }
    if(show_spatial_grid == TRUE){
      if(is.null(spatial_grid)){
        stop("No usable spatial grid specified! Please choose a network with spatial_grid_name=xxx")
      }
      else{
        if(is.null(spatial_grid_color)) {
          spatial_grid_color = 'black'
        }
        edges <- plotly_grid(spatial_grid)
        spl <- spl %>% plotly::add_segments(name = "spatial_grid",
                                            type = "scatter",
                                            data = edges,
                                            x = ~x,
                                            y = ~y,
                                            xend = ~x_end,
                                            yend = ~y_end,
                                            line = list(color = spatial_grid_color,
                                                        width = 1),
                                            opacity=spatial_grid_alpha)
        
      }
    }
    
    spl <- spl %>% plotly::add_trace(type = "scatter",mode = "markers",
                                     x = annotated_DT[[sdimx]],
                                     y = annotated_DT[[sdimy]],
                                     color = annotated_DT[[selected_genes]],
                                     colors = c(genes_low_color,genes_mid_color,genes_high_color),
                                     marker = list(size = spatial_point_size),
                                     showlegend = F)
    
    
    spl <- spl %>% plotly::layout(xaxis = list(title = "X"),
                                  yaxis = list(title = "Y"),
                                  legend = list(x = 100, y = 0.5,font = list(family = "sans-serif",size = legend_text_size)))
    
  }
  
  
  # 3D plot
  else{
    axis_scale = match.arg(axis_scale, c("cube","real","custom"))
    ratio = plotly_axis_scale_3D(annotated_DT,sdimx = sdimx,sdimy = sdimy,sdimz = sdimz,
                                 mode = axis_scale,custom_ratio = custom_ratio)
    
    
    spl <- plotly::plot_ly(scene = "scene2")
    
    spl <- spl %>% plotly::add_trace(type = 'scatter3d',mode = 'markers',
                                     x = annotated_DT[[sdimx]],
                                     y = annotated_DT[[sdimy]],
                                     z = annotated_DT[[sdimz]],
                                     color = annotated_DT[[selected_genes]],
                                     colors = c(genes_low_color,genes_mid_color,genes_high_color),
                                     #legendgroup = annotated_DT[[cell_color]],
                                     marker = list(size = spatial_point_size),
                                     showlegend = F)
    
    if(show_spatial_network == TRUE) {
      if(is.null(spatial_network)){
        stop("No usable spatial network specified! Please choose a network with spatial_network_name=xxx")
      }
      else{
        if(is.null(spatial_network_alpha)) {
          spatial_network_alpha = 0.5
        }
        else if(is.character(spatial_network_alpha)){
          warning("Edge_alpha for plotly mode is not adjustable yet. Default 0.5 will be set\n")
          spatial_network_alpha = 0.5
        }
        edges <- plotly_network(spatial_network)
        
        spl <- spl %>% plotly::add_trace(name = "sptial network",
                                         mode = "lines",
                                         type = "scatter3d",
                                         data = edges,
                                         x = ~x,y=~y,z=~z,
                                         line=list(color=network_color),
                                         opacity=spatial_network_alpha)
      }
    }
    
    if(show_spatial_grid == TRUE){
      cat("3D grid is not clear to view\n")
    }
    
  }
  
  
  
  spl <- plotly::hide_colorbar(spl)
  if(is.null(dim3_to_use) & is.null(sdimz)){
    if(plot_alignment == 'vertical'){
      combo_plot <- plotly::subplot(dpl,spl,nrows = 2,titleX = TRUE,titleY = TRUE)
    }
    else{
      combo_plot <- plotly::subplot(dpl,spl,titleX = TRUE,titleY = TRUE)
    }
  }
  
  else if(!is.null(dim3_to_use) & is.null(sdimz)){
    if(plot_alignment == 'vertical'){
      combo_plot <- plotly::subplot(dpl,spl,nrows = 2,titleX = TRUE,titleY = TRUE)%>%
        plotly::layout(scene = list(domain = list(x = c(0, 1), y = c(0,0.5)),
                                    xaxis = list(title = x_title),
                                    yaxis = list(title = y_title),
                                    zaxis = list(title = z_title)))
    }
    else{
      combo_plot <- plotly::subplot(dpl,spl,titleX = TRUE,titleY = TRUE) %>%
        plotly::layout(scene = list(domain = list(x = c(0, 0.5), y = c(0,1)),
                                    xaxis = list(title = x_title),
                                    yaxis = list(title = y_title),
                                    zaxis = list(title = z_title)))
    }
  }
  
  else if(is.null(dim3_to_use) & !is.null(sdimz)){
    if(plot_alignment == 'vertical'){
      combo_plot <- plotly::subplot(dpl,spl,nrows = 2,titleX = TRUE,titleY = TRUE)%>%
        plotly::layout(scene2 = list(
          xaxis = list(title = "X",nticks = x_ticks),
          yaxis = list(title = "Y",nticks = y_ticks),
          zaxis = list(title = "Z",nticks = z_ticks),
          aspectmode='manual',
          aspectratio = list(x=ratio[[1]],
                             y=ratio[[2]],
                             z=ratio[[3]])))
    }
    else{
      combo_plot <- plotly::subplot(dpl,spl,titleX = TRUE,titleY = TRUE) %>%
        plotly::layout(scene2 = list(
          xaxis = list(title = "X",nticks = x_ticks),
          yaxis = list(title = "Y",nticks = y_ticks),
          zaxis = list(title = "Z",nticks = z_ticks),
          aspectmode='manual',
          aspectratio = list(x=ratio[[1]],
                             y=ratio[[2]],
                             z=ratio[[3]])))
    }
  }
  
  else if(!is.null(dim3_to_use) & !is.null(sdimz)){
    if(plot_alignment == 'vertical'){
      combo_plot <- plotly::subplot(dpl,spl,nrows = 2,titleX = TRUE,titleY = TRUE)%>%
        plotly::layout(scene = list(domain = list(x = c(0, 1), y = c(0,0.5)),
                                    xaxis = list(title = x_title),
                                    yaxis = list(title = y_title),
                                    zaxis = list(title = z_title)),
                       scene2 = list(
                         xaxis = list(title = "X",nticks = x_ticks),
                         yaxis = list(title = "Y",nticks = y_ticks),
                         zaxis = list(title = "Z",nticks = z_ticks),
                         aspectmode='manual',
                         aspectratio = list(x=ratio[[1]],
                                            y=ratio[[2]],
                                            z=ratio[[3]])))
    }
    else{
      combo_plot <- plotly::subplot(dpl,spl,titleX = TRUE,titleY = TRUE) %>%
        plotly::layout(scene = list(domain = list(x = c(0, 0.5), y = c(0,1)),
                                    xaxis = list(title = x_title),
                                    yaxis = list(title = y_title),
                                    zaxis = list(title = z_title)),
                       scene2 = list(
                         xaxis = list(title = "X",nticks = x_ticks),
                         yaxis = list(title = "Y",nticks = y_ticks),
                         zaxis = list(title = "Z",nticks = z_ticks),
                         aspectmode='manual',
                         aspectratio = list(x=ratio[[1]],
                                            y=ratio[[2]],
                                            z=ratio[[3]])))
    }
  }
  
  
  
  return(combo_plot)
  
}


#' @title visSpatDimGenePlot
#' @name visSpatDimGenePlot
#' @description integration of visSpatDimGenePlot_2D(ggplot) and visSpatDimGenePlot_3D(plotly)
#' @param gobject giotto object
#' @param expression_values gene expression values to use
#' @param plot_alignment direction to align plot
#' @param genes genes to show
#' @param dim_reduction_to_use dimension reduction to use
#' @param dim_reduction_name dimension reduction name
#' @param dim1_to_use dimension to use on x-axis
#' @param dim2_to_use dimension to use on y-axis
#' @param dim3_to_use dimension to use on z-axis
#' @param sdimx x-axis dimension name (default = 'sdimx')
#' @param sdimy y-axis dimension name (default = 'sdimy')
#' @param sdimz z-axis dimension name (default = 'sdimz')
#' @param dim_point_size dim reduction plot: point size
#' @param dim_point_border_col color of border around points
#' @param dim_point_border_stroke stroke size of border around points
#' @param show_NN_network show underlying NN network
#' @param nn_network_to_use type of NN network to use (kNN vs sNN)
#' @param network_name name of NN network to use, if show_NN_network = TRUE
#' @param edge_alpha_dim dim reduction plot: column to use for alpha of the edges
#' @param scale_alpha_with_expression scale expression with ggplot alpha parameter
#' @param label_size size for the label
#' @param genes_high_color color to represent high expression of gene
#' @param genes_low_color color to represent low expression of gene
#' @param axis_scale three modes to adjust axis scale ratio
#' @param x_ticks number of ticks on x axis
#' @param y_ticks number of ticks on y axis
#' @param z_ticks number of ticks on z axis
#' @param custom_ratio set the axis scale ratio on custom
#' @param legend_text_size the size of the text in legend
#' @param spatial_network_name name of spatial network to use
#' @param spatial_grid_name name of spatial grid to use
#' @param spatial_point_size spatial plot: point size
#' @param spatial_point_border_col color of border around points
#' @param spatial_point_border_stroke stroke size of border around points
#' @param midpoint size of point (cell)
#' @param point_size size of point (cell)
#' @param cow_n_col cowplot param: how many columns
#' @param cow_rel_h cowplot param: relative height
#' @param cow_rel_w cowplot param: relative width
#' @param cow_align cowplot param: how to align
#' @param show_legend show legend
#' @param show_plot show plot
#' @return ggplot or plotly
#' @details Description of parameters.
#' @export
#' @examples
#'     visSpatDimGenePlot(gobject)
visSpatDimGenePlot <- function(gobject,
                               plot_method = c("ggplot", "plotly"),
                               expression_values = c('normalized', 'scaled', 'custom'),
                               plot_alignment = c('horizontal','vertical'),
                               dim_reduction_to_use = 'umap',
                               dim_reduction_name = 'umap',
                               dim1_to_use = 1,
                               dim2_to_use = 2,
                               dim3_to_use = NULL,
                               sdimx=NULL,
                               sdimy=NULL,
                               sdimz=NULL,
                               genes,
                               dim_point_border_col = 'black',
                               dim_point_border_stroke = 0.1,
                               show_NN_network = F,
                               nn_network_to_use = 'sNN',
                               network_name = 'sNN.pca',
                               edge_alpha_dim = NULL,
                               scale_alpha_with_expression = FALSE,
                               label_size = 16,
                               genes_low_color = "blue",
                               genes_mid_color = "white",
                               genes_high_color = "red",
                               dim_point_size = 3,
                               nn_network_alpha = 0.5,
                               show_spatial_network = F,
                               spatial_network_name = 'spatial_network',
                               network_color = "lightgray",
                               spatial_network_alpha = 0.5,
                               show_spatial_grid = F,
                               spatial_grid_name = 'spatial_grid',
                               spatial_grid_color = NULL,
                               spatial_grid_alpha = 0.5,
                               spatial_point_size = 3,
                               spatial_point_border_col = 'black',
                               spatial_point_border_stroke = 0.1,
                               legend_text_size = 12,
                               axis_scale = c("cube","real","custom"),
                               custom_ratio = NULL,
                               x_ticks = NULL,
                               y_ticks = NULL,
                               z_ticks = NULL,
                               midpoint = 0,
                               point_size = 1,
                               cow_n_col = 2,
                               cow_rel_h = 1,
                               cow_rel_w = 1,
                               cow_align = 'h',
                               show_legend = T,
                               show_plots = F){
  plot_method = match.arg(plot_method, c("ggplot", "plotly"))
  plot_alignment = match.arg(plot_alignment, c( 'horizontal','vertical'))
  
  if(plot_method == "ggplot"){
    if(!is.null(sdimz) | !is.null(dim3_to_use)){
      stop("3D plots couldn't be drawn in ggplot mode,please choose plotly mode with plot_mode = \"plotly\"\n")
    }
    
    
    if(is.null(label_size)){
      label_size = 4
    }
    if(is.null(nn_network_alpha)){
      nn_network_alpha = 0.05
    }
    pl <- visSpatDimGenePlot_2D (gobject = gobject,
                                 expression_values = expression_values,
                                 plot_alignment = plot_alignment,
                                 genes = genes,
                                 dim_reduction_to_use = dim_reduction_to_use,
                                 dim_reduction_name = dim_reduction_name,
                                 dim1_to_use = dim1_to_use,
                                 dim2_to_use = dim2_to_use,
                                 point_size = point_size,
                                 dim_point_border_col =dim_point_border_col,
                                 dim_point_border_stroke = dim_point_border_stroke,
                                 show_NN_network = show_NN_network,
                                 show_spatial_network = show_spatial_network,
                                 show_spatial_grid = show_spatial_grid,
                                 nn_network_to_use = nn_network_to_use,
                                 network_name =network_name,
                                 edge_alpha_dim = edge_alpha_dim,
                                 scale_alpha_with_expression = scale_alpha_with_expression,
                                 spatial_network_name = spatial_network_name,
                                 spatial_grid_name = spatial_grid_name,
                                 spatial_point_size = spatial_point_size,
                                 spatial_point_border_col = spatial_point_border_col,
                                 spatial_point_border_stroke = spatial_point_border_stroke,
                                 midpoint = midpoint,
                                 genes_high_color = genes_high_color,
                                 genes_mid_color=genes_mid_color,
                                 genes_low_color = genes_low_color,
                                 cow_n_col = cow_n_col,
                                 cow_rel_h = cow_rel_h,
                                 cow_rel_w = cow_rel_w,
                                 cow_align = cow_align,
                                 axis_scale = axis_scale,
                                 custom_ratio = custom_ratio,
                                 x_ticks = x_ticks,
                                 y_ticks = y_ticks,
                                 show_legend = show_legend,
                                 show_plots = show_plots)
  }
  
  else if(plot_method == "plotly"){
    if(is.null(label_size)){
      label_size = 16
    }
    if(is.null(nn_network_alpha)){
      nn_network_alpha = 0.5
    }
    if(length(genes) > 1){
      warning("\n Now plotly mode can just accept one gene, only the first gene will be plot\n")
      genes = genes[1]
    }
    pl <- visSpatDimGenePlot_3D(gobject = gobject,
                                expression_values = expression_values,
                                plot_alignment = plot_alignment,
                                dim_reduction_to_use = dim_reduction_to_use,
                                dim_reduction_name = dim_reduction_name,
                                dim1_to_use = dim1_to_use,
                                dim2_to_use = dim2_to_use,
                                dim3_to_use = dim3_to_use,
                                sdimx=sdimx,
                                sdimy=sdimy,
                                sdimz=sdimz,
                                genes = genes,
                                show_NN_network = show_NN_network,
                                nn_network_to_use = nn_network_to_use,
                                network_name = network_name,
                                label_size = label_size,
                                genes_low_color = genes_low_color,
                                genes_mid_color = genes_mid_color,
                                genes_high_color = genes_high_color,
                                dim_point_size = dim_point_size,
                                nn_network_alpha = nn_network_alpha,
                                show_spatial_network = show_spatial_network,
                                spatial_network_name = spatial_network_name,
                                network_color = network_color,
                                spatial_network_alpha = spatial_network_alpha,
                                show_spatial_grid = show_spatial_grid,
                                spatial_grid_name = spatial_grid_name,
                                spatial_grid_color = spatial_grid_color,
                                spatial_grid_alpha = spatial_grid_alpha,
                                spatial_point_size = spatial_point_size,
                                legend_text_size = legend_text_size,
                                axis_scale = axis_scale,
                                custom_ratio = custom_ratio,
                                x_ticks = x_ticks,
                                y_ticks = y_ticks,
                                z_ticks = z_ticks)
  }
  return(pl)
}




## 2-D ggplots ####
## ----------- ##


#' @title plot_network_layer_ggplot
#' @name plot_network_layer_ggplot
#' @description Visualize cells in network layer according to dimension reduction coordinates
#' @param gobject giotto object
#' @param annotated_network_DT annotated network data.table of selected cells
#' @param edge_alpha alpha of network edges
#' @param show_legend show legend
#' @return ggplot
#' @details Description of parameters.
#' @export
#' @examples
#'     plot_network_layer_ggplot(gobject)
plot_network_layer_ggplot = function(ggobject,
                                     annotated_network_DT,
                                     edge_alpha = NULL,
                                     show_legend = T) {


  from_dims = grep('from_Dim', colnames(annotated_network_DT), value = T)
  to_dims = grep('to_Dim', colnames(annotated_network_DT), value = T)



  pl <- ggobject

  if(is.null(edge_alpha)) {
    edge_alpha = 0.05
    pl <- pl + ggplot2::geom_segment(data = annotated_network_DT, aes_string2(x = from_dims[1], y = from_dims[2],
                                                                             xend = to_dims[1], yend = to_dims[2]),
                                     alpha = edge_alpha, show.legend = show_legend)

  } else if(is.numeric(edge_alpha)) {
    pl <- pl + ggplot2::geom_segment(data = annotated_network_DT, aes_string2(x = from_dims[1], y = from_dims[2],
                                                                             xend = to_dims[1], yend = to_dims[2]),
                                     alpha = edge_alpha, show.legend = show_legend)
  } else if(is.character(edge_alpha)) {

    if(edge_alpha %in% colnames(annotated_network_DT)) {
      pl <- pl + ggplot2::geom_segment(data = annotated_network_DT, aes_string2(x = from_dims[1], y = from_dims[2],
                                                                               xend = to_dims[1], yend = to_dims[2],
                                                                               alpha = edge_alpha),
                                       show.legend = show_legend)
    }
  }

  return(pl)

}




#' @title plot_point_layer_ggplot
#' @name plot_point_layer_ggplot
#' @description Visualize cells in point layer according to dimension reduction coordinates
#' @param gobject giotto object
#' @param annotated_DT_selected annotated data.table of selected cells
#' @param annotated_DT_other annotated data.table of not selected cells
#' @param cell_color color for cells (see details)
#' @param color_as_factor convert color column to factor
#' @param cell_color_code named vector with colors
#' @param cell_color_gradient vector with 3 colors for numeric data
#' @param gradient_midpoint midpoint for color gradient
#' @param gradient_limits vector with lower and upper limits
#' @param select_cell_groups select subset of cells/clusters based on cell_color parameter
#' @param select_cells select subset of cells based on cell IDs
#' @param point_size size of point (cell)
#' @param point_border_col color of border around points
#' @param point_border_stroke stroke size of border around points
#' @param show_cluster_center plot center of selected clusters
#' @param show_center_label plot label of selected clusters
#' @param center_point_size size of center points
#' @param label_size  size of labels
#' @param label_fontface font of labels
#' @param edge_alpha column to use for alpha of the edges
#' @param show_other_cells display not selected cells
#' @param other_cell_color color of not selected cells
#' @param other_point_size size of not selected cells
#' @param show_legend show legend
#' @return ggplot
#' @details Description of parameters.
#' @export
#' @examples
#'     plot_point_layer_ggplot(gobject)
plot_point_layer_ggplot = function(ggobject,
                                   annotated_DT_selected,
                                   annotated_DT_other,
                                   cell_color = NULL,
                                   color_as_factor = T,
                                   cell_color_code = NULL,
                                   cell_color_gradient = c('blue', 'white', 'red'),
                                   gradient_midpoint = 0,
                                   gradient_limits = NULL,
                                   select_cell_groups = NULL,
                                   select_cells = NULL,
                                   point_size = 1,
                                   point_border_col = 'black',
                                   point_border_stroke = 0.1,
                                   show_cluster_center = F,
                                   show_center_label = T,
                                   center_point_size = 4,
                                   center_point_border_col = 'black',
                                   center_point_border_stroke = 0.1,
                                   label_size = 4,
                                   label_fontface = 'bold',
                                   edge_alpha = NULL,
                                   show_other_cells = T,
                                   other_cell_color = 'lightgrey',
                                   other_point_size = 0.5,
                                   show_legend = T
) {


  pl = ggobject



  ## first plot other non-selected cells
  if((!is.null(select_cells) | !is.null(select_cell_groups)) & show_other_cells == TRUE) {

    dims = grep('Dim.', colnames(annotated_DT_other), value = T)
    pl = pl + ggplot2::geom_point(data = annotated_DT_other, aes_string(x = dims[1], dims[2]),
                                  color = other_cell_color, show.legend = F, size = other_point_size)

  }


  ## order of color
  # 1. if NULL then default to lightblue
  # 2. if character vector
  # 2.1 if length of cell_color is longer than 1 and has colors
  # 2.2 if not part of metadata then suppose its color
  # 2.3 part of metadata
  # 2.3.1 numerical column
  # 2.3.2 factor column or character to factor


  ## point layer
  dims = grep('Dim.', colnames(annotated_DT_selected), value = T)

  if(is.null(cell_color)) {

    cell_color = 'lightblue'
    pl <- pl + ggplot2::geom_point(data = annotated_DT_selected, aes_string(x = dims[1], dims[2]),
                                   color = cell_color, show.legend = show_legend, size = point_size)


  } else if(length(cell_color) > 1) {

    if(is.numeric(cell_color) | is.factor(cell_color)) {
      if(nrow(annotated_DT_selected) != length(cell_color)) stop('\n vector needs to be the same lengths as number of cells \n')
      annotated_DT_selected[['temp_color']] = cell_color

      pl <- pl + ggplot2::geom_point(data = annotated_DT_selected, aes_string2(x = dims[1], y = dims[2], fill = 'temp_color'),
                                     show.legend = show_legend, shape = 21,
                                     size = point_size,
                                     color = point_border_col, stroke = point_border_stroke)

    } else if(is.character(cell_color)) {
      if(!all(cell_color %in% grDevices::colors())) stop('cell_color is not numeric, a factor or vector of colors \n')
      pl <- pl + ggplot2::geom_point(data = annotated_DT_selected, aes_string2(x = dims[1], y = dims[2]),
                                     show.legend = show_legend, shape = 21, fill = cell_color,
                                     size = point_size,
                                     color = point_border_col, stroke = point_border_stroke)

    }

  } else if (is.character(cell_color)) {

    if(!cell_color %in% colnames(annotated_DT_selected)) {
      if(!cell_color %in% grDevices::colors()) stop(cell_color,' is not a color or a column name \n')
      pl <- pl + ggplot2::geom_point(data = annotated_DT_selected, aes_string(x = dims[1], y = dims[2]),
                                     show.legend = show_legend, shape = 21, fill = cell_color,
                                     size = point_size,
                                     color = point_border_col, stroke = point_border_stroke)

    } else {

      class_cell_color = class(annotated_DT_selected[[cell_color]])

      if((class_cell_color == 'integer' | class_cell_color == 'numeric') & color_as_factor == FALSE) {

        # set upper and lower limits
        if(!is.null(gradient_limits) & is.vector(gradient_limits) & length(gradient_limits) == 2) {
          lower_lim = gradient_limits[[1]]
          upper_lim = gradient_limits[[2]]

          numeric_data = annotated_DT_selected[[cell_color]]
          limit_numeric_data = ifelse(numeric_data > upper_lim, upper_lim,
                                      ifelse(numeric_data < lower_lim, lower_lim, numeric_data))
          annotated_DT_selected[[cell_color]] = limit_numeric_data
        }

        pl <- pl + ggplot2::geom_point(data = annotated_DT_selected,
                                       aes_string2(x = dims[1], y = dims[2], fill = cell_color),
                                       show.legend = show_legend, shape = 21, size = point_size,
                                       color = point_border_col, stroke = point_border_stroke)

      } else {

        # convert character or numeric to factor
        if(color_as_factor == TRUE) {
          factor_data = factor(annotated_DT_selected[[cell_color]])
          annotated_DT_selected[[cell_color]] <- factor_data
        }

        # if you want to show centers or labels then calculate centers
        if(show_cluster_center == TRUE | show_center_label == TRUE) {
          annotated_DT_centers = annotated_DT_selected[, .(center_1 = median(get(dims[1])), center_2 = median(get(dims[2]))), by = cell_color]
          factor_center_data = factor(annotated_DT_centers[[cell_color]])
          annotated_DT_centers[[cell_color]] <- factor_center_data
        }

        pl <- pl + ggplot2::geom_point(data = annotated_DT_selected,
                                       aes_string2(x = dims[1], y = dims[2], fill = cell_color),
                                       show.legend = show_legend, shape = 21, size = point_size,
                                       color = point_border_col, stroke = point_border_stroke)


        ## plot centers
        if(show_cluster_center == TRUE & (color_as_factor == TRUE | class_cell_color %in% c('character', 'factor'))) {

          pl <- pl + ggplot2::geom_point(data = annotated_DT_centers,
                                         aes_string2(x = 'center_1', y = 'center_2', fill = cell_color),
                                         color = center_point_border_col, stroke = center_point_border_stroke,
                                         size = center_point_size, shape = 21)
        }

        ## plot labels
        if(show_center_label == TRUE) {
          pl <- pl + ggrepel::geom_text_repel(data = annotated_DT_centers,
                                              aes_string2(x = 'center_1', y = 'center_2', label = cell_color),
                                              size = label_size, fontface = label_fontface)
        }

      }


      ## specificy colors to use
      if(!is.null(cell_color_code)) {

        pl <- pl + ggplot2::scale_fill_manual(values = cell_color_code)

      } else if(color_as_factor == T) {

        number_colors = length(unique(factor_data))
        cell_color_code = Giotto:::getDistinctColors(n = number_colors)
        names(cell_color_code) = unique(factor_data)
        pl <- pl + ggplot2::scale_fill_manual(values = cell_color_code)

      } else if(color_as_factor == F){

        if(is.null(gradient_midpoint)) {
          gradient_midpoint = median(annotated_DT_selected[[cell_color]])
        }
        pl <- pl + ggplot2::scale_fill_gradient2(low = cell_color_gradient[[1]],
                                                 mid = cell_color_gradient[[2]],
                                                 high = cell_color_gradient[[3]],
                                                 midpoint = gradient_midpoint)

      }
    }
  }
  return(pl)
}




#' @title plot_point_layer_ggplot_noFILL
#' @name plot_point_layer_ggplot_noFILL
#' @description Visualize cells in point layer according to dimension reduction coordinates without borders
#' @param gobject giotto object
#' @param annotated_DT_selected annotated data.table of selected cells
#' @param annotated_DT_other annotated data.table of not selected cells
#' @param cell_color color for cells (see details)
#' @param color_as_factor convert color column to factor
#' @param cell_color_code named vector with colors
#' @param cell_color_gradient vector with 3 colors for numeric data
#' @param gradient_midpoint midpoint for color gradient
#' @param gradient_limits vector with lower and upper limits
#' @param select_cell_groups select subset of cells/clusters based on cell_color parameter
#' @param select_cells select subset of cells based on cell IDs
#' @param point_size size of point (cell)
#' @param show_cluster_center plot center of selected clusters
#' @param show_center_label plot label of selected clusters
#' @param center_point_size size of center points
#' @param label_size  size of labels
#' @param label_fontface font of labels
#' @param edge_alpha column to use for alpha of the edges
#' @param show_other_cells display not selected cells
#' @param other_cell_color color of not selected cells
#' @param other_point_size size of not selected cells
#' @param show_legend show legend
#' @return ggplot
#' @details Description of parameters.
#' @export
#' @examples
#'     plot_point_layer_ggplot_noFILL(gobject)
plot_point_layer_ggplot_noFILL = function(ggobject,
                                   annotated_DT_selected,
                                   annotated_DT_other,
                                   cell_color = NULL,
                                   color_as_factor = T,
                                   cell_color_code = NULL,
                                   cell_color_gradient = c('blue', 'white', 'red'),
                                   gradient_midpoint = 0,
                                   gradient_limits = NULL,
                                   select_cell_groups = NULL,
                                   select_cells = NULL,
                                   point_size = 1,
                                   show_cluster_center = F,
                                   show_center_label = T,
                                   center_point_size = 4,
                                   label_size = 4,
                                   label_fontface = 'bold',
                                   edge_alpha = NULL,
                                   show_other_cells = T,
                                   other_cell_color = 'lightgrey',
                                   other_point_size = 0.5,
                                   show_legend = T
) {


  pl = ggobject



  ## first plot other non-selected cells
  if((!is.null(select_cells) | !is.null(select_cell_groups)) & show_other_cells == TRUE) {

    dims = grep('Dim.', colnames(annotated_DT_other), value = T)
    pl = pl + ggplot2::geom_point(data = annotated_DT_other, aes_string(x = dims[1], dims[2]),
                                  color = other_cell_color, show.legend = F, size = other_point_size)

  }


  ## order of color
  # 1. if NULL then default to lightblue
  # 2. if character vector
  # 2.1 if length of cell_color is longer than 1 and has colors
  # 2.2 if not part of metadata then suppose its color
  # 2.3 part of metadata
  # 2.3.1 numerical column
  # 2.3.2 factor column or character to factor


  ## point layer
  dims = grep('Dim.', colnames(annotated_DT_selected), value = T)

  if(is.null(cell_color)) {

    cell_color = 'lightblue'
    pl <- pl + ggplot2::geom_point(data = annotated_DT_selected, aes_string(x = dims[1], dims[2]),
                                   color = cell_color, show.legend = show_legend, size = point_size)


  } else if(length(cell_color) > 1) {

    if(is.numeric(cell_color) | is.factor(cell_color)) {
      if(nrow(annotated_DT_selected) != length(cell_color)) stop('\n vector needs to be the same lengths as number of cells \n')
      annotated_DT_selected[['temp_color']] = cell_color

      pl <- pl + ggplot2::geom_point(data = annotated_DT_selected, aes_string2(x = dims[1], y = dims[2], color = 'temp_color'),
                                     show.legend = show_legend, shape = 19, size = point_size)

    } else if(is.character(cell_color)) {
      if(!all(cell_color %in% grDevices::colors())) stop('cell_color is not numeric, a factor or vector of colors \n')
      pl <- pl + ggplot2::geom_point(data = annotated_DT_selected, aes_string2(x = dims[1], y = dims[2]),
                                     show.legend = show_legend, shape = 19, fill = cell_color, size = point_size)

    }

  } else if (is.character(cell_color)) {

    if(!cell_color %in% colnames(annotated_DT_selected)) {
      if(!cell_color %in% grDevices::colors()) stop(cell_color,' is not a color or a column name \n')
      pl <- pl + ggplot2::geom_point(data = annotated_DT_selected, aes_string(x = dims[1], y = dims[2]),
                                     show.legend = show_legend, shape = 19, color = cell_color, size = point_size)

    } else {

      class_cell_color = class(annotated_DT_selected[[cell_color]])

      if((class_cell_color == 'integer' | class_cell_color == 'numeric') & color_as_factor == FALSE) {

        # set upper and lower limits
        if(!is.null(gradient_limits) & is.vector(gradient_limits) & length(gradient_limits) == 2) {
          lower_lim = gradient_limits[[1]]
          upper_lim = gradient_limits[[2]]

          numeric_data = annotated_DT_selected[[cell_color]]
          limit_numeric_data = ifelse(numeric_data > upper_lim, upper_lim,
                                      ifelse(numeric_data < lower_lim, lower_lim, numeric_data))
          annotated_DT_selected[[cell_color]] = limit_numeric_data
        }

        pl <- pl + ggplot2::geom_point(data = annotated_DT_selected,
                                       aes_string2(x = dims[1], y = dims[2], color = cell_color),
                                       show.legend = show_legend, shape = 19, size = point_size)

      } else {

        # convert character or numeric to factor
        if(color_as_factor == TRUE) {
          factor_data = factor(annotated_DT_selected[[cell_color]])
          annotated_DT_selected[[cell_color]] <- factor_data
        }

        # if you want to show centers or labels then calculate centers
        if(show_cluster_center == TRUE | show_center_label == TRUE) {
          annotated_DT_centers = annotated_DT_selected[, .(center_1 = median(get(dims[1])), center_2 = median(get(dims[2]))), by = cell_color]
          factor_center_data = factor(annotated_DT_centers[[cell_color]])
          annotated_DT_centers[[cell_color]] <- factor_center_data
        }

        pl <- pl + ggplot2::geom_point(data = annotated_DT_selected,
                                       aes_string2(x = dims[1], y = dims[2], color = cell_color),
                                       show.legend = show_legend, shape = 19, size = point_size)


        ## plot centers
        if(show_cluster_center == TRUE & (color_as_factor == TRUE | class_cell_color %in% c('character', 'factor'))) {

          pl <- pl + ggplot2::geom_point(data = annotated_DT_centers,
                                         aes_string2(x = 'center_1', y = 'center_2', color = cell_color),
                                         size = center_point_size, shape = 19)
        }

        ## plot labels
        if(show_center_label == TRUE) {
          pl <- pl + ggrepel::geom_text_repel(data = annotated_DT_centers,
                                              aes_string2(x = 'center_1', y = 'center_2', label = cell_color),
                                              size = label_size, fontface = label_fontface)
        }

      }


      ## specificy colors to use
      if(!is.null(cell_color_code)) {

        pl <- pl + ggplot2::scale_color_manual(values = cell_color_code)

      } else if(color_as_factor == T) {

        number_colors = length(unique(factor_data))
        cell_color_code = Giotto:::getDistinctColors(n = number_colors)
        names(cell_color_code) = unique(factor_data)
        pl <- pl + ggplot2::scale_color_manual(values = cell_color_code)

      } else if(color_as_factor == F){

        if(is.null(gradient_midpoint)) {
          gradient_midpoint = median(annotated_DT_selected[[cell_color]])
        }
        pl <- pl + ggplot2::scale_color_gradient2(low = cell_color_gradient[[1]],
                                                 mid = cell_color_gradient[[2]],
                                                 high = cell_color_gradient[[3]],
                                                 midpoint = gradient_midpoint)

      }
    }
  }
  return(pl)
}




#' @title dimPlot2D_single
#' @name dimPlot2D_single
#' @description Visualize cells according to dimension reduction coordinates
#' @param gobject giotto object
#' @param dim_reduction_to_use dimension reduction to use
#' @param dim_reduction_name dimension reduction name
#' @param dim1_to_use dimension to use on x-axis
#' @param dim2_to_use dimension to use on y-axis
#' @param spat_enr_names names of spatial enrichment results to include
#' @param show_NN_network show underlying NN network
#' @param nn_network_to_use type of NN network to use (kNN vs sNN)
#' @param network_name name of NN network to use, if show_NN_network = TRUE
#' @param cell_color color for cells (see details)
#' @param color_as_factor convert color column to factor
#' @param cell_color_code named vector with colors
#' @param cell_color_gradient vector with 3 colors for numeric data
#' @param gradient_midpoint midpoint for color gradient
#' @param gradient_limits vector with lower and upper limits
#' @param select_cell_groups select subset of cells/clusters based on cell_color parameter
#' @param select_cells select subset of cells based on cell IDs
#' @param show_other_cells display not selected cells
#' @param other_cell_color color of not selected cells
#' @param other_point_size size of not selected cells
#' @param show_cluster_center plot center of selected clusters
#' @param show_center_label plot label of selected clusters
#' @param center_point_size size of center points
#' @param label_size  size of labels
#' @param label_fontface font of labels
#' @param edge_alpha column to use for alpha of the edges
#' @param point_shape point with border or not (border or no_border)
#' @param point_size size of point (cell)
#' @param point_border_col color of border around points
#' @param point_border_stroke stroke size of border around points
#' @param title title for plot, defaults to cell_color parameter
#' @param show_legend show legend
#' @param legend_text size of legend text
#' @param legend_symbol_size size of legend symbols
#' @param background_color color of plot background
#' @param axis_text size of axis text
#' @param axis_title size of axis title
#' @param show_plot show plot
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters from \code{\link{all_plots_save_function}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return ggplot
#' @details Description of parameters. For 3D plots see \code{\link{dimPlot3D}}
#' @export
#' @examples
#'     dimPlot2D_single(gobject)
dimPlot2D_single <- function(gobject,
                             dim_reduction_to_use = 'umap',
                             dim_reduction_name = 'umap',
                             dim1_to_use = 1,
                             dim2_to_use = 2,
                             spat_enr_names = NULL,
                             show_NN_network = F,
                             nn_network_to_use = 'sNN',
                             network_name = 'sNN.pca',
                             cell_color = NULL,
                             color_as_factor = T,
                             cell_color_code = NULL,
                             cell_color_gradient = c('blue', 'white', 'red'),
                             gradient_midpoint = NULL,
                             gradient_limits = NULL,
                             select_cell_groups = NULL,
                             select_cells = NULL,
                             show_other_cells = T,
                             other_cell_color = 'lightgrey',
                             other_point_size = 0.5,
                             show_cluster_center = F,
                             show_center_label = T,
                             center_point_size = 4,
                             center_point_border_col = 'black',
                             center_point_border_stroke = 0.1,
                             label_size = 4,
                             label_fontface = 'bold',
                             edge_alpha = NULL,
                             point_shape = c('border', 'no_border'),
                             point_size = 1,
                             point_border_col = 'black',
                             point_border_stroke = 0.1,
                             title = NULL,
                             show_legend = T,
                             legend_text = 8,
                             legend_symbol_size = 1,
                             background_color = 'white',
                             axis_text = 8,
                             axis_title = 8,
                             show_plot = NA,
                             return_plot = NA,
                             save_plot = NA,
                             save_param = list(),
                             default_save_name = 'dimPlot2D_single'
){


  ## point shape ##
  point_shape = match.arg(point_shape, c('border', 'no_border'))

  ## dimension reduction ##
  # test if dimension reduction was performed
  if(is.null(gobject@dimension_reduction$cells[[dim_reduction_to_use]][[dim_reduction_name]])) {
    stop('\n dimension reduction: ', dim_reduction_to_use, ' or dimension reduction name: ',dim_reduction_name,' is not available \n')
  }
  dim_dfr = gobject@dimension_reduction$cells[[dim_reduction_to_use]][[dim_reduction_name]]$coordinates[,c(dim1_to_use, dim2_to_use)]
  dim_names = colnames(dim_dfr)
  dim_DT = data.table::as.data.table(dim_dfr); dim_DT[, cell_ID := rownames(dim_dfr)]

  ## annotated cell metadata
  cell_metadata = combineMetadata(gobject = gobject,
                                  spat_enr_names = spat_enr_names)
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

  # add % variance information if reduction is PCA
  if(dim_reduction_to_use == "pca"){
    eigenvaluesDT = data.table::as.data.table(gobject@dimension_reduction$cells[[dim_reduction_to_use]][[dim_reduction_name]]$misc$eig)
    var_expl_vec = eigenvaluesDT[c(dim1_to_use, dim2_to_use)][['percentage of variance']]
    dim1_x_variance = var_expl_vec[1]
    dim2_y_variance = var_expl_vec[2]
  }



  ## create subsets if needed
  if(!is.null(select_cells) & !is.null(select_cell_groups)) {
    if(is.null(cell_color)) {
      stop('\n selection of cells is based on cell_color paramter, which is a metadata column \n')
    }
    cat('You have selected both individual cell IDs and a group of cells \n')
    group_cell_IDs = annotated_DT[get(cell_color) %in% select_cell_groups][['cell_ID']]
    select_cells = unique(c(select_cells, group_cell_IDs))
  } else if(!is.null(select_cell_groups)) {
    select_cells = annotated_DT[get(cell_color) %in% select_cell_groups][['cell_ID']]
  }

  if(!is.null(select_cells)) {
    annotated_DT_other = annotated_DT[!annotated_DT$cell_ID %in% select_cells]
    annotated_DT_selected = annotated_DT[annotated_DT$cell_ID %in% select_cells]

    if(show_NN_network == TRUE) {
      annotated_network_DT <- annotated_network_DT[annotated_network_DT$to %in% select_cells & annotated_network_DT$from %in% select_cells]
    }

    # if specific cells are selected
    annotated_DT = annotated_DT_selected
  }

  ## if no subsets are required
  if(is.null(select_cells) & is.null(select_cell_groups)) {
    annotated_DT_selected = annotated_DT
    annotated_DT_other    = NULL
  }



  pl <- ggplot2::ggplot()
  pl <- pl + ggplot2::theme_classic()

  ## add network layer
  if(show_NN_network == TRUE) {
    pl = plot_network_layer_ggplot(ggobject = pl,
                                   annotated_network_DT = annotated_network_DT,
                                   edge_alpha = edge_alpha,
                                   show_legend = show_legend)
  }

  #return(list(pl, annotated_DT_selected, annotated_DT_other))

  if(point_shape == 'border') {
    ## add point layer
    pl = plot_point_layer_ggplot(ggobject = pl,
                                 annotated_DT_selected = annotated_DT_selected,
                                 annotated_DT_other = annotated_DT_other,
                                 cell_color = cell_color,
                                 color_as_factor = color_as_factor,
                                 cell_color_code = cell_color_code,
                                 cell_color_gradient = cell_color_gradient,
                                 gradient_midpoint = gradient_midpoint,
                                 gradient_limits = gradient_limits,
                                 select_cell_groups = select_cell_groups,
                                 select_cells = select_cells,
                                 show_other_cells = show_other_cells,
                                 other_cell_color = other_cell_color,
                                 other_point_size = other_point_size,
                                 show_cluster_center = show_cluster_center,
                                 show_center_label = show_center_label,
                                 center_point_size = center_point_size,
                                 center_point_border_col = center_point_border_col,
                                 center_point_border_stroke = center_point_border_stroke,
                                 label_size = label_size,
                                 label_fontface = label_fontface,
                                 edge_alpha = edge_alpha,
                                 point_size = point_size,
                                 point_border_col = point_border_col,
                                 point_border_stroke = point_border_stroke,
                                 show_legend = show_legend)

  } else if(point_shape == 'no_border') {

    pl = plot_point_layer_ggplot_noFILL(ggobject = pl,
                                        annotated_DT_selected = annotated_DT_selected,
                                        annotated_DT_other = annotated_DT_other,
                                        cell_color = cell_color,
                                        color_as_factor = color_as_factor,
                                        cell_color_code = cell_color_code,
                                        cell_color_gradient = cell_color_gradient,
                                        gradient_midpoint = gradient_midpoint,
                                        gradient_limits = gradient_limits,
                                        select_cell_groups = select_cell_groups,
                                        select_cells = select_cells,
                                        show_other_cells = show_other_cells,
                                        other_cell_color = other_cell_color,
                                        other_point_size = other_point_size,
                                        show_cluster_center = show_cluster_center,
                                        show_center_label = show_center_label,
                                        center_point_size = center_point_size,
                                        label_size = label_size,
                                        label_fontface = label_fontface,
                                        edge_alpha = edge_alpha,
                                        point_size = point_size,
                                        show_legend = show_legend)

  }




  ## add % variance explained to names of plot for PCA ##
  if(dim_reduction_to_use == 'pca') {
    x_name = paste0('pca','-',dim_names[1])
    y_name = paste0('pca','-',dim_names[2])

    # provide x, y and plot titles
    x_title = sprintf('%s explains %.02f%% of variance', x_name, var_expl_vec[1])
    y_title = sprintf('%s explains %.02f%% of variance', y_name, var_expl_vec[2])

    if(is.null(title)) title = cell_color
    pl <- pl + ggplot2::labs(x = x_title, y = y_title, title = title)

  } else {

    # provide x, y and plot titles
    x_title = paste0(dim_reduction_to_use,'-',dim_names[1])
    y_title = paste0(dim_reduction_to_use,'-',dim_names[2])

    if(is.null(title)) title = cell_color
    pl <- pl + ggplot2::labs(x = x_title, y = y_title, title = title)

  }

  ## adjust titles
  pl <- pl + ggplot2::theme(plot.title = element_text(hjust = 0.5),
                            legend.title = element_blank(),
                            legend.text = element_text(size = legend_text),
                            axis.text = element_text(size = axis_text),
                            axis.title = element_text(size = axis_title),
                            panel.grid = element_blank(),
                            panel.background = element_rect(fill = background_color))

  ## change symbol size of legend
  if(color_as_factor == TRUE) {
    if(point_shape == 'border') {
      pl = pl + guides(fill = guide_legend(override.aes = list(size = legend_symbol_size)))
    } else if(point_shape == 'no_border') {
      pl = pl + guides(color = guide_legend(override.aes = list(size = legend_symbol_size)))
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




#' @title dimPlot2D
#' @name dimPlot2D
#' @description Visualize cells according to dimension reduction coordinates
#' @param gobject giotto object
#' @param groub_by create multiple plots based on cell annotation column
#' @param group_by_subset subset the group_by factor column
#' @param dim_reduction_to_use dimension reduction to use
#' @param dim_reduction_name dimension reduction name
#' @param dim1_to_use dimension to use on x-axis
#' @param dim2_to_use dimension to use on y-axis
#' @param spat_enr_names names of spatial enrichment results to include
#' @param show_NN_network show underlying NN network
#' @param nn_network_to_use type of NN network to use (kNN vs sNN)
#' @param network_name name of NN network to use, if show_NN_network = TRUE
#' @param cell_color color for cells (see details)
#' @param color_as_factor convert color column to factor
#' @param cell_color_code named vector with colors
#' @param cell_color_gradient vector with 3 colors for numeric data
#' @param gradient_midpoint midpoint for color gradient
#' @param gradient_limits vector with lower and upper limits
#' @param select_cell_groups select subset of cells/clusters based on cell_color parameter
#' @param select_cells select subset of cells based on cell IDs
#' @param show_other_cells display not selected cells
#' @param other_cell_color color of not selected cells
#' @param other_point_size size of not selected cells
#' @param show_cluster_center plot center of selected clusters
#' @param show_center_label plot label of selected clusters
#' @param center_point_size size of center points
#' @param label_size  size of labels
#' @param label_fontface font of labels
#' @param edge_alpha column to use for alpha of the edges
#' @param point_shape point with border or not (border or no_border)
#' @param point_size size of point (cell)
#' @param point_border_col color of border around points
#' @param point_border_stroke stroke size of border around points
#' @param title title for plot, defaults to cell_color parameter
#' @param show_legend show legend
#' @param legend_text size of legend text
#' @param legend_symbol_size size of legend symbols
#' @param background_color color of plot background
#' @param axis_text size of axis text
#' @param axis_title size of axis title
#' @param cow_n_col cowplot param: how many columns
#' @param cow_rel_h cowplot param: relative height
#' @param cow_rel_w cowplot param: relative width
#' @param cow_align cowplot param: how to align
#' @param show_plot show plot
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters from \code{\link{all_plots_save_function}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return ggplot
#' @details Description of parameters. For 3D plots see \code{\link{dimPlot3D}}
#' @export
#' @examples
#'     dimPlot2D(gobject)
dimPlot2D = function(gobject,
                     group_by = NULL,
                     group_by_subset = NULL,
                     dim_reduction_to_use = 'umap',
                     dim_reduction_name = 'umap',
                     dim1_to_use = 1,
                     dim2_to_use = 2,
                     spat_enr_names = NULL,
                     show_NN_network = F,
                     nn_network_to_use = 'sNN',
                     network_name = 'sNN.pca',
                     cell_color = NULL,
                     color_as_factor = T,
                     cell_color_code = NULL,
                     cell_color_gradient = c('blue', 'white', 'red'),
                     gradient_midpoint = NULL,
                     gradient_limits = NULL,
                     select_cell_groups = NULL,
                     select_cells = NULL,
                     show_other_cells = T,
                     other_cell_color = 'lightgrey',
                     other_point_size = 0.5,
                     show_cluster_center = F,
                     show_center_label = T,
                     center_point_size = 4,
                     center_point_border_col = 'black',
                     center_point_border_stroke = 0.1,
                     label_size = 4,
                     label_fontface = 'bold',
                     edge_alpha = NULL,
                     point_shape = c('border', 'no_border'),
                     point_size = 1,
                     point_border_col = 'black',
                     point_border_stroke = 0.1,
                     title = NULL,
                     show_legend = T,
                     legend_text = 8,
                     legend_symbol_size = 1,
                     background_color = 'white',
                     axis_text = 8,
                     axis_title = 8,
                     cow_n_col = 2,
                     cow_rel_h = 1,
                     cow_rel_w = 1,
                     cow_align = 'h',
                     show_plot = NA,
                     return_plot = NA,
                     save_plot = NA,
                     save_param = list(),
                     default_save_name = 'dimPlot2D') {


  ## check group_by
  if(is.null(group_by)) {

    dimPlot2D_single(gobject = gobject,
                     dim_reduction_to_use = dim_reduction_to_use,
                     dim_reduction_name = dim_reduction_name,
                     dim1_to_use = dim1_to_use,
                     dim2_to_use = dim2_to_use,
                     spat_enr_names = spat_enr_names,
                     show_NN_network = show_NN_network,
                     nn_network_to_use = nn_network_to_use,
                     network_name = network_name,
                     cell_color = cell_color,
                     color_as_factor = color_as_factor,
                     cell_color_code = cell_color_code,
                     cell_color_gradient = cell_color_gradient,
                     gradient_midpoint = gradient_midpoint,
                     gradient_limits = gradient_limits,
                     select_cell_groups = select_cell_groups,
                     select_cells = select_cells,
                     show_other_cells = show_other_cells,
                     other_cell_color = other_cell_color,
                     other_point_size = other_point_size,
                     show_cluster_center = show_cluster_center,
                     show_center_label = show_center_label,
                     center_point_size = center_point_size,
                     center_point_border_col = center_point_border_col,
                     center_point_border_stroke = center_point_border_stroke,
                     label_size = label_size,
                     label_fontface = label_fontface,
                     edge_alpha = edge_alpha,
                     point_shape = point_shape,
                     point_size = point_size,
                     point_border_col = point_border_col,
                     point_border_stroke = point_border_stroke,
                     title = title,
                     show_legend = show_legend,
                     legend_text = legend_text,
                     legend_symbol_size = legend_symbol_size,
                     background_color = background_color,
                     axis_text = axis_text,
                     axis_title = axis_title,
                     show_plot = show_plot,
                     return_plot = return_plot,
                     save_plot = save_plot,
                     save_param = save_param,
                     default_save_name = default_save_name)



  } else {

    comb_metadata = combineMetadata(gobject = gobject,
                                    spat_enr_names = spat_enr_names)
    possible_meta_groups = colnames(comb_metadata)

    ## check if group_by is found
    if(!group_by %in% possible_meta_groups) {
      stop("group_by ", group_by, " was not found in pDataDT()")
    }

    unique_groups = unique(comb_metadata[[group_by]])

    # subset unique_groups
    if(!is.null(group_by_subset)) {
      not_found = group_by_subset[!group_by_subset %in% unique_groups]

      if(length(not_found) > 0) {
        cat('the following subset was not found: ', not_found)
      }
      unique_groups = unique_groups[unique_groups %in% group_by_subset]
    }


    # create matching cell_color_code
    if(is.null(cell_color_code)) {
      if(is.character(cell_color)) {

        if(cell_color %in% colnames(comb_metadata)) {

          if(color_as_factor == TRUE) {
            number_colors = length(unique(comb_metadata[[cell_color]]))
            cell_color_code = Giotto:::getDistinctColors(n = number_colors)
            names(cell_color_code) = unique(comb_metadata[[cell_color]])
            cell_color_code = cell_color_code
          }
        }
      }
    }




    # print, return and save parameters
    show_plot = ifelse(is.na(show_plot), readGiottoInstructions(gobject, param = 'show_plot'), show_plot)
    save_plot = ifelse(is.na(save_plot), readGiottoInstructions(gobject, param = 'save_plot'), save_plot)
    return_plot = ifelse(is.na(return_plot), readGiottoInstructions(gobject, param = 'return_plot'), return_plot)

    ## plotting ##
    savelist <- list()


    for(group_id in 1:length(unique_groups)) {

      group = unique_groups[group_id]

      subset_cell_IDs = comb_metadata[get(group_by) == group][['cell_ID']]
      temp_gobject = subsetGiotto(gobject = gobject, cell_ids = subset_cell_IDs)

      pl = dimPlot2D_single(gobject = temp_gobject,
                            dim_reduction_to_use = dim_reduction_to_use,
                            dim_reduction_name = dim_reduction_name,
                            dim1_to_use = dim1_to_use,
                            dim2_to_use = dim2_to_use,
                            spat_enr_names = spat_enr_names,
                            show_NN_network = show_NN_network,
                            nn_network_to_use = nn_network_to_use,
                            network_name = network_name,
                            cell_color = cell_color,
                            cell_color_code = cell_color_code,
                            color_as_factor = color_as_factor,
                            cell_color_gradient = cell_color_gradient,
                            gradient_midpoint = gradient_midpoint,
                            gradient_limits = gradient_limits,
                            select_cell_groups = select_cell_groups,
                            select_cells = select_cells,
                            show_other_cells = show_other_cells,
                            other_cell_color = other_cell_color,
                            other_point_size = other_point_size,
                            show_cluster_center = show_cluster_center,
                            show_center_label = show_center_label,
                            center_point_size = center_point_size,
                            center_point_border_col = center_point_border_col,
                            center_point_border_stroke = center_point_border_stroke,
                            label_size = label_size,
                            label_fontface = label_fontface,
                            edge_alpha = edge_alpha,
                            point_shape = point_shape,
                            point_size = point_size,
                            point_border_col = point_border_col,
                            point_border_stroke = point_border_stroke,
                            title = group,
                            show_legend = show_legend,
                            legend_text = legend_text,
                            legend_symbol_size = legend_symbol_size,
                            background_color = background_color,
                            axis_text = axis_text,
                            axis_title = axis_title,
                            show_plot = FALSE,
                            return_plot = TRUE,
                            save_plot = FALSE,
                            save_param = list(),
                            default_save_name = 'dimPlot2D_single')


      savelist[[group_id]] <- pl


    }

    # combine plots with cowplot
    combo_plot <- cowplot::plot_grid(plotlist = savelist,
                                     ncol = cow_n_col,
                                     rel_heights = cow_rel_h,
                                     rel_widths = cow_rel_w,
                                     align = cow_align)

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

}







#' @title dimPlot
#' @name dimPlot
#' @description Visualize cells according to dimension reduction coordinates
#' @param gobject giotto object
#' @param groub_by create multiple plots based on cell annotation column
#' @param group_by_subset subset the group_by factor column
#' @param dim_reduction_to_use dimension reduction to use
#' @param dim_reduction_name dimension reduction name
#' @param dim1_to_use dimension to use on x-axis
#' @param dim2_to_use dimension to use on y-axis
#' @param spat_enr_names names of spatial enrichment results to include
#' @param show_NN_network show underlying NN network
#' @param nn_network_to_use type of NN network to use (kNN vs sNN)
#' @param network_name name of NN network to use, if show_NN_network = TRUE
#' @param cell_color color for cells (see details)
#' @param color_as_factor convert color column to factor
#' @param cell_color_code named vector with colors
#' @param cell_color_gradient vector with 3 colors for numeric data
#' @param gradient_midpoint midpoint for color gradient
#' @param gradient_limits vector with lower and upper limits
#' @param select_cell_groups select subset of cells/clusters based on cell_color parameter
#' @param select_cells select subset of cells based on cell IDs
#' @param show_other_cells display not selected cells
#' @param other_cell_color color of not selected cells
#' @param other_point_size size of not selected cells
#' @param show_cluster_center plot center of selected clusters
#' @param show_center_label plot label of selected clusters
#' @param center_point_size size of center points
#' @param label_size  size of labels
#' @param label_fontface font of labels
#' @param edge_alpha column to use for alpha of the edges
#' @param point_shape point with border or not (border or no_border)
#' @param point_size size of point (cell)
#' @param point_border_col color of border around points
#' @param point_border_stroke stroke size of border around points
#' @param title title for plot, defaults to cell_color parameter
#' @param show_legend show legend
#' @param legend_text size of legend text
#' @param legend_symbol_size size of legend symbols
#' @param background_color color of plot background
#' @param axis_text size of axis text
#' @param axis_title size of axis title
#' @param cow_n_col cowplot param: how many columns
#' @param cow_rel_h cowplot param: relative height
#' @param cow_rel_w cowplot param: relative width
#' @param cow_align cowplot param: how to align
#' @param show_plot show plot
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters from \code{\link{all_plots_save_function}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return ggplot
#' @details Description of parameters, see \code{\link{dimPlot2D}}. For 3D plots see \code{\link{dimPlot3D}}
#' @export
#' @examples
#'     dimPlot(gobject)
dimPlot = function(gobject,
                   group_by = NULL,
                   group_by_subset = NULL,
                   dim_reduction_to_use = 'umap',
                   dim_reduction_name = 'umap',
                   dim1_to_use = 1,
                   dim2_to_use = 2,
                   spat_enr_names = NULL,
                   show_NN_network = F,
                   nn_network_to_use = 'sNN',
                   network_name = 'sNN.pca',
                   cell_color = NULL,
                   color_as_factor = T,
                   cell_color_code = NULL,
                   cell_color_gradient = c('blue', 'white', 'red'),
                   gradient_midpoint = NULL,
                   gradient_limits = NULL,
                   select_cell_groups = NULL,
                   select_cells = NULL,
                   show_other_cells = T,
                   other_cell_color = 'lightgrey',
                   other_point_size = 0.5,
                   show_cluster_center = F,
                   show_center_label = T,
                   center_point_size = 4,
                   center_point_border_col = 'black',
                   center_point_border_stroke = 0.1,
                   label_size = 4,
                   label_fontface = 'bold',
                   edge_alpha = NULL,
                   point_shape = c('border', 'no_border'),
                   point_size = 1,
                   point_border_col = 'black',
                   point_border_stroke = 0.1,
                   show_legend = T,
                   legend_text = 8,
                   legend_symbol_size = 1,
                   background_color = 'white',
                   axis_text = 8,
                   axis_title = 8,
                   title = NULL,
                   cow_n_col = 2,
                   cow_rel_h = 1,
                   cow_rel_w = 1,
                   cow_align = 'h',
                   show_plot = NA,
                   return_plot = NA,
                   save_plot = NA,
                   save_param = list(),
                   default_save_name = 'dimPlot') {

  dimPlot2D(gobject = gobject,
            group_by = group_by,
            group_by_subset = group_by_subset,
            dim_reduction_to_use = dim_reduction_to_use,
            dim_reduction_name = dim_reduction_name,
            dim1_to_use = dim1_to_use,
            dim2_to_use = dim2_to_use,
            spat_enr_names = spat_enr_names,
            show_NN_network = show_NN_network,
            nn_network_to_use = nn_network_to_use,
            network_name = network_name,
            cell_color = cell_color,
            color_as_factor = color_as_factor,
            cell_color_code = cell_color_code,
            cell_color_gradient = cell_color_gradient,
            gradient_midpoint = gradient_midpoint,
            gradient_limits = gradient_limits,
            select_cell_groups = select_cell_groups,
            select_cells = select_cells,
            show_other_cells = show_other_cells,
            other_cell_color = other_cell_color,
            other_point_size = other_point_size,
            show_cluster_center = show_cluster_center,
            show_center_label = show_center_label,
            center_point_size = center_point_size,
            center_point_border_col = center_point_border_col,
            center_point_border_stroke = center_point_border_stroke,
            label_size = label_size,
            label_fontface = label_fontface,
            edge_alpha = edge_alpha,
            point_shape = point_shape,
            point_size = point_size,
            point_border_col = point_border_col,
            point_border_stroke = point_border_stroke,
            title = title,
            show_legend = show_legend,
            legend_text = legend_text,
            legend_symbol_size = legend_symbol_size,
            background_color = background_color,
            axis_text = axis_text,
            axis_title = axis_title,
            cow_n_col = cow_n_col,
            cow_rel_h = cow_rel_h,
            cow_rel_w = cow_rel_w,
            cow_align = cow_align,
            show_plot = show_plot,
            return_plot = return_plot,
            save_plot = save_plot,
            save_param = save_param,
            default_save_name = default_save_name)

}






#' @title plotUMAP_2D
#' @name plotUMAP_2D
#' @description Short wrapper for UMAP visualization
#' @param gobject giotto object
#' @param groub_by create multiple plots based on cell annotation column
#' @param group_by_subset subset the group_by factor column
#' @param dim_reduction_name dimension reduction name
#' @param dim1_to_use dimension to use on x-axis
#' @param dim2_to_use dimension to use on y-axis
#' @param spat_enr_names names of spatial enrichment results to include
#' @param show_NN_network show underlying NN network
#' @param nn_network_to_use type of NN network to use (kNN vs sNN)
#' @param network_name name of NN network to use, if show_NN_network = TRUE
#' @param cell_color color for cells (see details)
#' @param color_as_factor convert color column to factor
#' @param cell_color_code named vector with colors
#' @param cell_color_gradient vector with 3 colors for numeric data
#' @param gradient_midpoint midpoint for color gradient
#' @param gradient_limits vector with lower and upper limits
#' @param select_cell_groups select subset of cells/clusters based on cell_color parameter
#' @param select_cells select subset of cells based on cell IDs
#' @param show_other_cells display not selected cells
#' @param other_cell_color color of not selected cells
#' @param other_point_size size of not selected cells
#' @param show_cluster_center plot center of selected clusters
#' @param show_center_label plot label of selected clusters
#' @param center_point_size size of center points
#' @param label_size  size of labels
#' @param label_fontface font of labels
#' @param edge_alpha column to use for alpha of the edges
#' @param point_shape point with border or not (border or no_border)
#' @param point_size size of point (cell)
#' @param point_border_col color of border around points
#' @param point_border_stroke stroke size of border around points
#' @param title title for plot, defaults to cell_color parameter
#' @param show_legend show legend
#' @param legend_text size of legend text
#' @param legend_symbol_size size of legend symbols
#' @param background_color color of plot background
#' @param axis_text size of axis text
#' @param axis_title size of axis title
#' @param cow_n_col cowplot param: how many columns
#' @param cow_rel_h cowplot param: relative height
#' @param cow_rel_w cowplot param: relative width
#' @param cow_align cowplot param: how to align
#' @param show_plot show plot
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters from \code{\link{all_plots_save_function}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return ggplot
#' @details Description of parameters, see \code{\link{dimPlot2D}}. For 3D plots see \code{\link{plotUMAP_3D}}
#' @export
#' @examples
#'     plotUMAP_2D(gobject)
plotUMAP_2D = function(gobject, dim_reduction_name = 'umap', default_save_name = 'UMAP_2D', ...) {

  dimPlot2D(gobject = gobject,
            dim_reduction_to_use = 'umap',  dim_reduction_name = dim_reduction_name,
            default_save_name = default_save_name, ...)

}


#' @title plotUMAP
#' @name plotUMAP
#' @description Short wrapper for UMAP visualization
#' @param gobject giotto object
#' @param groub_by create multiple plots based on cell annotation column
#' @param group_by_subset subset the group_by factor column
#' @param dim_reduction_name dimension reduction name
#' @param dim1_to_use dimension to use on x-axis
#' @param dim2_to_use dimension to use on y-axis
#' @param spat_enr_names names of spatial enrichment results to include
#' @param show_NN_network show underlying NN network
#' @param nn_network_to_use type of NN network to use (kNN vs sNN)
#' @param network_name name of NN network to use, if show_NN_network = TRUE
#' @param cell_color color for cells (see details)
#' @param color_as_factor convert color column to factor
#' @param cell_color_code named vector with colors
#' @param cell_color_gradient vector with 3 colors for numeric data
#' @param gradient_midpoint midpoint for color gradient
#' @param gradient_limits vector with lower and upper limits
#' @param select_cell_groups select subset of cells/clusters based on cell_color parameter
#' @param select_cells select subset of cells based on cell IDs
#' @param show_other_cells display not selected cells
#' @param other_cell_color color of not selected cells
#' @param other_point_size size of not selected cells
#' @param show_cluster_center plot center of selected clusters
#' @param show_center_label plot label of selected clusters
#' @param center_point_size size of center points
#' @param label_size  size of labels
#' @param label_fontface font of labels
#' @param edge_alpha column to use for alpha of the edges
#' @param point_shape point with border or not (border or no_border)
#' @param point_size size of point (cell)
#' @param point_border_col color of border around points
#' @param point_border_stroke stroke size of border around points
#' @param title title for plot, defaults to cell_color parameter
#' @param show_legend show legend
#' @param legend_text size of legend text
#' @param legend_symbol_size size of legend symbols
#' @param background_color color of plot background
#' @param axis_text size of axis text
#' @param axis_title size of axis title
#' @param cow_n_col cowplot param: how many columns
#' @param cow_rel_h cowplot param: relative height
#' @param cow_rel_w cowplot param: relative width
#' @param cow_align cowplot param: how to align
#' @param show_plot show plot
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters from \code{\link{all_plots_save_function}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return ggplot
#' @details Description of parameters, see \code{\link{dimPlot2D}}. For 3D plots see \code{\link{plotUMAP_3D}}
#' @export
#' @examples
#'     plotUMAP(gobject)
plotUMAP = function(gobject, dim_reduction_name = 'umap', default_save_name = 'UMAP', ...) {

  dimPlot2D(gobject = gobject,
            dim_reduction_to_use = 'umap',  dim_reduction_name = dim_reduction_name,
            default_save_name = default_save_name, ...)

}





#' @title plotTSNE_2D
#' @name plotTSNE_2D
#' @description Short wrapper for tSNE visualization
#' @param gobject giotto object
#' @param groub_by create multiple plots based on cell annotation column
#' @param group_by_subset subset the group_by factor column
#' @param dim_reduction_name dimension reduction name
#' @param dim1_to_use dimension to use on x-axis
#' @param dim2_to_use dimension to use on y-axis
#' @param spat_enr_names names of spatial enrichment results to include
#' @param show_NN_network show underlying NN network
#' @param nn_network_to_use type of NN network to use (kNN vs sNN)
#' @param network_name name of NN network to use, if show_NN_network = TRUE
#' @param cell_color color for cells (see details)
#' @param color_as_factor convert color column to factor
#' @param cell_color_code named vector with colors
#' @param cell_color_gradient vector with 3 colors for numeric data
#' @param gradient_midpoint midpoint for color gradient
#' @param gradient_limits vector with lower and upper limits
#' @param select_cell_groups select subset of cells/clusters based on cell_color parameter
#' @param select_cells select subset of cells based on cell IDs
#' @param show_other_cells display not selected cells
#' @param other_cell_color color of not selected cells
#' @param other_point_size size of not selected cells
#' @param show_cluster_center plot center of selected clusters
#' @param show_center_label plot label of selected clusters
#' @param center_point_size size of center points
#' @param label_size  size of labels
#' @param label_fontface font of labels
#' @param edge_alpha column to use for alpha of the edges
#' @param point_shape point with border or not (border or no_border)
#' @param point_size size of point (cell)
#' @param point_border_col color of border around points
#' @param point_border_stroke stroke size of border around points
#' @param title title for plot, defaults to cell_color parameter
#' @param show_legend show legend
#' @param legend_text size of legend text
#' @param legend_symbol_size size of legend symbols
#' @param background_color color of plot background
#' @param axis_text size of axis text
#' @param axis_title size of axis title
#' @param cow_n_col cowplot param: how many columns
#' @param cow_rel_h cowplot param: relative height
#' @param cow_rel_w cowplot param: relative width
#' @param cow_align cowplot param: how to align
#' @param show_plot show plot
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters from \code{\link{all_plots_save_function}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return ggplot
#' @details Description of parameters, see \code{\link{dimPlot2D}}. For 3D plots see \code{\link{plotTSNE_3D}}
#' @export
#' @examples
#'     plotTSNE_2D(gobject)
plotTSNE_2D = function(gobject, dim_reduction_name = 'tsne', default_save_name = 'tSNE_2D', ...) {

  dimPlot2D(gobject = gobject,
            dim_reduction_to_use = 'tsne', dim_reduction_name = dim_reduction_name,
            default_save_name = default_save_name, ...)

}

#' @title plotTSNE
#' @name plotTSNE
#' @description Short wrapper for tSNE visualization
#' @param gobject giotto object
#' @param groub_by create multiple plots based on cell annotation column
#' @param group_by_subset subset the group_by factor column
#' @param dim_reduction_name dimension reduction name
#' @param dim1_to_use dimension to use on x-axis
#' @param dim2_to_use dimension to use on y-axis
#' @param spat_enr_names names of spatial enrichment results to include
#' @param show_NN_network show underlying NN network
#' @param nn_network_to_use type of NN network to use (kNN vs sNN)
#' @param network_name name of NN network to use, if show_NN_network = TRUE
#' @param cell_color color for cells (see details)
#' @param color_as_factor convert color column to factor
#' @param cell_color_code named vector with colors
#' @param cell_color_gradient vector with 3 colors for numeric data
#' @param gradient_midpoint midpoint for color gradient
#' @param gradient_limits vector with lower and upper limits
#' @param select_cell_groups select subset of cells/clusters based on cell_color parameter
#' @param select_cells select subset of cells based on cell IDs
#' @param show_other_cells display not selected cells
#' @param other_cell_color color of not selected cells
#' @param other_point_size size of not selected cells
#' @param show_cluster_center plot center of selected clusters
#' @param show_center_label plot label of selected clusters
#' @param center_point_size size of center points
#' @param label_size  size of labels
#' @param label_fontface font of labels
#' @param edge_alpha column to use for alpha of the edges
#' @param point_shape point with border or not (border or no_border)
#' @param point_size size of point (cell)
#' @param point_border_col color of border around points
#' @param point_border_stroke stroke size of border around points
#' @param title title for plot, defaults to cell_color parameter
#' @param show_legend show legend
#' @param legend_text size of legend text
#' @param legend_symbol_size size of legend symbols
#' @param background_color color of plot background
#' @param axis_text size of axis text
#' @param axis_title size of axis title
#' @param cow_n_col cowplot param: how many columns
#' @param cow_rel_h cowplot param: relative height
#' @param cow_rel_w cowplot param: relative width
#' @param cow_align cowplot param: how to align
#' @param show_plot show plot
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters from \code{\link{all_plots_save_function}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return ggplot
#' @details Description of parameters, see \code{\link{dimPlot2D}}. For 3D plots see \code{\link{plotTSNE_3D}}
#' @export
#' @examples
#'     plotTSNE(gobject)
plotTSNE = function(gobject, dim_reduction_name = 'tsne', default_save_name = 'tSNE', ...) {

  dimPlot2D(gobject = gobject,
            dim_reduction_to_use = 'tsne', dim_reduction_name = dim_reduction_name,
            default_save_name = default_save_name, ...)

}



#' @title plotPCA_2D
#' @name plotPCA_2D
#' @description Short wrapper for PCA visualization
#' @param gobject giotto object
#' @param groub_by create multiple plots based on cell annotation column
#' @param group_by_subset subset the group_by factor column
#' @param dim_reduction_name dimension reduction name
#' @param dim1_to_use dimension to use on x-axis
#' @param dim2_to_use dimension to use on y-axis
#' @param spat_enr_names names of spatial enrichment results to include
#' @param show_NN_network show underlying NN network
#' @param nn_network_to_use type of NN network to use (kNN vs sNN)
#' @param network_name name of NN network to use, if show_NN_network = TRUE
#' @param cell_color color for cells (see details)
#' @param color_as_factor convert color column to factor
#' @param cell_color_code named vector with colors
#' @param cell_color_gradient vector with 3 colors for numeric data
#' @param gradient_midpoint midpoint for color gradient
#' @param gradient_limits vector with lower and upper limits
#' @param select_cell_groups select subset of cells/clusters based on cell_color parameter
#' @param select_cells select subset of cells based on cell IDs
#' @param show_other_cells display not selected cells
#' @param other_cell_color color of not selected cells
#' @param other_point_size size of not selected cells
#' @param show_cluster_center plot center of selected clusters
#' @param show_center_label plot label of selected clusters
#' @param center_point_size size of center points
#' @param label_size  size of labels
#' @param label_fontface font of labels
#' @param edge_alpha column to use for alpha of the edges
#' @param point_shape point with border or not (border or no_border)
#' @param point_size size of point (cell)
#' @param point_border_col color of border around points
#' @param point_border_stroke stroke size of border around points
#' @param title title for plot, defaults to cell_color parameter
#' @param show_legend show legend
#' @param legend_text size of legend text
#' @param legend_symbol_size size of legend symbols
#' @param background_color color of plot background
#' @param axis_text size of axis text
#' @param axis_title size of axis title
#' @param cow_n_col cowplot param: how many columns
#' @param cow_rel_h cowplot param: relative height
#' @param cow_rel_w cowplot param: relative width
#' @param cow_align cowplot param: how to align
#' @param show_plot show plot
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters from \code{\link{all_plots_save_function}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return ggplot
#' @details Description of parameters, see \code{\link{dimPlot2D}}. For 3D plots see \code{\link{plotPCA_3D}}
#' @export
#' @examples
#'     plotPCA_2D(gobject)
plotPCA_2D = function(gobject, dim_reduction_name = 'pca', default_save_name = 'PCA_2D', ...) {

  dimPlot2D(gobject = gobject,
            dim_reduction_to_use = 'pca', dim_reduction_name = dim_reduction_name,
            default_save_name = default_save_name, ...)

}



#' @title plotPCA
#' @name plotPCA
#' @description Short wrapper for PCA visualization
#' @param gobject giotto object
#' @param groub_by create multiple plots based on cell annotation column
#' @param group_by_subset subset the group_by factor column
#' @param dim_reduction_name dimension reduction name
#' @param dim1_to_use dimension to use on x-axis
#' @param dim2_to_use dimension to use on y-axis
#' @param spat_enr_names names of spatial enrichment results to include
#' @param show_NN_network show underlying NN network
#' @param nn_network_to_use type of NN network to use (kNN vs sNN)
#' @param network_name name of NN network to use, if show_NN_network = TRUE
#' @param cell_color color for cells (see details)
#' @param color_as_factor convert color column to factor
#' @param cell_color_code named vector with colors
#' @param cell_color_gradient vector with 3 colors for numeric data
#' @param gradient_midpoint midpoint for color gradient
#' @param gradient_limits vector with lower and upper limits
#' @param select_cell_groups select subset of cells/clusters based on cell_color parameter
#' @param select_cells select subset of cells based on cell IDs
#' @param show_other_cells display not selected cells
#' @param other_cell_color color of not selected cells
#' @param other_point_size size of not selected cells
#' @param show_cluster_center plot center of selected clusters
#' @param show_center_label plot label of selected clusters
#' @param center_point_size size of center points
#' @param label_size  size of labels
#' @param label_fontface font of labels
#' @param edge_alpha column to use for alpha of the edges
#' @param point_shape point with border or not (border or no_border)
#' @param point_size size of point (cell)
#' @param point_border_col color of border around points
#' @param point_border_stroke stroke size of border around points
#' @param show_legend show legend
#' @param title title for plot, defaults to cell_color parameter
#' @param legend_text size of legend text
#' @param legend_symbol_size size of legend symbols
#' @param background_color color of plot background
#' @param axis_text size of axis text
#' @param axis_title size of axis title
#' @param cow_n_col cowplot param: how many columns
#' @param cow_rel_h cowplot param: relative height
#' @param cow_rel_w cowplot param: relative width
#' @param cow_align cowplot param: how to align
#' @param show_plot show plot
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters from \code{\link{all_plots_save_function}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return ggplot
#' @details Description of parameters, see \code{\link{dimPlot2D}}. For 3D plots see \code{\link{plotPCA_3D}}
#' @export
#' @examples
#'     plotPCA(gobject)
plotPCA = function(gobject, dim_reduction_name = 'pca', default_save_name = 'PCA', ...) {

  dimPlot2D(gobject = gobject,
            dim_reduction_to_use = 'pca', dim_reduction_name = dim_reduction_name,
            default_save_name = default_save_name, ...)
}




#' @title plot_spat_point_layer_ggplot
#' @name plot_spat_point_layer_ggplot
#' @description creat ggplot point layer for spatial coordinates
#' @param gobject giotto object
#' @param sdimx x-axis dimension name (default = 'sdimx')
#' @param sdimy y-axis dimension name (default = 'sdimy')
#' @param cell_locations_metadata_selected annotated location from selected cells
#' @param cell_locations_metadata_other annotated location from non-selected cells
#' @param cell_color color for cells (see details)
#' @param color_as_factor convert color column to factor
#' @param cell_color_code named vector with colors
#' @param cell_color_gradient vector with 3 colors for numeric data
#' @param gradient_midpoint midpoint for color gradient
#' @param gradient_limits vector with lower and upper limits
#' @param select_cell_groups select subset of cells/clusters based on cell_color parameter
#' @param select_cells select subset of cells based on cell IDs
#' @param point_size size of point (cell)
#' @param point_border_col color of border around points
#' @param point_border_stroke stroke size of border around points
#' @param show_cluster_center plot center of selected clusters
#' @param show_center_label plot label of selected clusters
#' @param center_point_size size of center points
#' @param label_size  size of labels
#' @param label_fontface font of labels
#' @param show_other_cells display not selected cells
#' @param other_cell_color color for not selected cells
#' @param other_point_size point size for not selected cells
#' @param show_legend show legend
#' @return ggplot
#' @details Description of parameters.
#' @export
#' @examples
#'     plot_spat_point_layer_ggplot(gobject)
plot_spat_point_layer_ggplot = function(ggobject,
                                        sdimx = NULL,
                                        sdimy = NULL,
                                        cell_locations_metadata_selected,
                                        cell_locations_metadata_other,
                                        cell_color = NULL,
                                        color_as_factor = T,
                                        cell_color_code = NULL,
                                        cell_color_gradient = c('blue', 'white', 'red'),
                                        gradient_midpoint = NULL,
                                        gradient_limits = NULL,
                                        select_cell_groups = NULL,
                                        select_cells = NULL,
                                        point_size = 2,
                                        point_border_col = 'lightgrey',
                                        point_border_stroke = 0.1,
                                        show_cluster_center = F,
                                        show_center_label = T,
                                        center_point_size = 4,
                                        center_point_border_col = 'black',
                                        center_point_border_stroke = 0.1,
                                        label_size = 4,
                                        label_fontface = 'bold',
                                        show_other_cells = T,
                                        other_cell_color = 'lightgrey',
                                        other_point_size = 1,
                                        show_legend = TRUE

) {

  ## specify spatial dimensions first
  if(is.null(sdimx) | is.null(sdimy)) {

    warning("plot_method = ggplot, but spatial dimensions for sdimx and/or sdimy are not specified. \n
            It will default to the 'sdimx' and 'sdimy' ")
    sdimx = 'sdimx'
    sdimy = 'sdimy'
  }

  ## ggplot object
  pl = ggobject

  ## first plot other non-selected cells
  if((!is.null(select_cells) | !is.null(select_cell_groups)) & show_other_cells == TRUE) {
    pl <- pl + ggplot2::geom_point(data = cell_locations_metadata_other, aes_string(x = sdimx, sdimy),
                                   color = other_cell_color, show.legend = F, size = other_point_size)
  }


  ## order of color
  # 1. if NULL then default to lightblue
  # 2. if character vector
  # 2.1 if length of cell_color is longer than 1 and has colors
  # 2.2 if not part of metadata then suppose its color
  # 2.3 part of metadata
  # 2.3.1 numerical column
  # 2.3.2 factor column or character to factor


  # cell color default
  if(is.null(cell_color)) {

    cell_color = 'lightblue'
    pl <- pl + ggplot2::geom_point(data = cell_locations_metadata_selected,
                                   aes_string(x = sdimx, y = sdimy),
                                   show.legend = show_legend, shape = 21,
                                   fill = cell_color, size = point_size,
                                   stroke = point_border_stroke, color = point_border_col)


  } else if(length(cell_color) > 1) {

    if(is.numeric(cell_color) | is.factor(cell_color)) {
      if(nrow(cell_locations_metadata_selected) != length(cell_color)) stop('\n vector needs to be the same lengths as number of cells \n')
      cell_locations_metadata_selected[['temp_color']] = cell_color

      pl <- pl + ggplot2::geom_point(data = cell_locations_metadata_selected,
                                     aes_string2(x = sdimx, y = sdimy, fill = 'temp_color'),
                                     show.legend = show_legend, shape = 21,
                                     size = point_size,
                                     color = point_border_col, stroke = point_border_stroke)

    } else if(is.character(cell_color)) {
      if(!all(cell_color %in% grDevices::colors())) stop('cell_color is not numeric, a factor or vector of colors \n')
      pl <- pl + ggplot2::geom_point(data = cell_locations_metadata_selected,
                                     aes_string2(x = sdimx, y = sdimy),
                                     show.legend = show_legend, shape = 21, fill = cell_color,
                                     size = point_size,
                                     color = point_border_col, stroke = point_border_stroke)

    }

  } else if(is.character(cell_color)) {
    if(!cell_color %in% colnames(cell_locations_metadata_selected)) {
      if(!cell_color %in% grDevices::colors()) stop(cell_color,' is not a color or a column name \n')
      pl <- pl + ggplot2::geom_point(data = cell_locations_metadata_selected,
                                     aes_string2(x = sdimx, y = sdimy),
                                     show.legend = show_legend, shape = 21, fill = cell_color,
                                     size = point_size,
                                     color = point_border_col, stroke = point_border_stroke)

    } else {

      class_cell_color = class(cell_locations_metadata_selected[[cell_color]])

      if((class_cell_color == 'integer' | class_cell_color == 'numeric') & color_as_factor == FALSE) {
        # set upper and lower limits
        if(!is.null(gradient_limits) & is.vector(gradient_limits) & length(gradient_limits) == 2) {
          lower_lim = gradient_limits[[1]]
          upper_lim = gradient_limits[[2]]

          numeric_data = cell_locations_metadata_selected[[cell_color]]
          limit_numeric_data = ifelse(numeric_data > upper_lim, upper_lim,
                                      ifelse(numeric_data < lower_lim, lower_lim, numeric_data))
          cell_locations_metadata_selected[[cell_color]] = limit_numeric_data
        }

        pl <- pl + ggplot2::geom_point(data = cell_locations_metadata_selected,
                                       aes_string2(x = sdimx, y = sdimy, fill = cell_color),
                                       show.legend = show_legend, shape = 21,
                                       size = point_size,
                                       color = point_border_col, stroke = point_border_stroke)



      } else {

        # convert character or numeric to factor
        if(color_as_factor == TRUE) {
          factor_data = factor(cell_locations_metadata_selected[[cell_color]])
          cell_locations_metadata_selected[[cell_color]] <- factor_data
        }

        # if you want to show centers or labels then calculate centers
        if(show_cluster_center == TRUE | show_center_label == TRUE) {
          annotated_DT_centers = cell_locations_metadata_selected[, .(center_1 = median(get('sdimx')), center_2 = median(get('sdimy'))), by = cell_color]
          factor_center_data = factor(annotated_DT_centers[[cell_color]])
          annotated_DT_centers[[cell_color]] <- factor_center_data
        }

        pl <- pl + ggplot2::geom_point(data = cell_locations_metadata_selected,
                                       aes_string2(x = sdimx, y = sdimy, fill = cell_color),
                                       show.legend = show_legend, shape = 21, size = point_size,
                                       color = point_border_col, stroke = point_border_stroke)


        ## plot centers
        if(show_cluster_center == TRUE & (color_as_factor == TRUE | class_cell_color %in% c('character', 'factor'))) {

          pl <- pl + ggplot2::geom_point(data = annotated_DT_centers,
                                         aes_string2(x = 'center_1', y = 'center_2', fill = cell_color),
                                         color = center_point_border_col, stroke = center_point_border_stroke,
                                         size = center_point_size, shape = 21)
        }

        ## plot labels
        if(show_center_label == TRUE) {
          pl <- pl + ggrepel::geom_text_repel(data = annotated_DT_centers,
                                              aes_string2(x = 'center_1', y = 'center_2', label = cell_color),
                                              size = label_size, fontface = label_fontface)
        }

      }

      ## specificy colors to use
      if(!is.null(cell_color_code)) {

        pl <- pl + ggplot2::scale_fill_manual(values = cell_color_code)

      } else if(color_as_factor == T) {

        number_colors = length(unique(factor_data))
        cell_color_code = Giotto:::getDistinctColors(n = number_colors)
        names(cell_color_code) = unique(factor_data)
        pl <- pl + ggplot2::scale_fill_manual(values = cell_color_code)

      } else if(color_as_factor == F){

        if(is.null(gradient_midpoint)) {
          gradient_midpoint = median(cell_locations_metadata_selected[[cell_color]])
        }

        pl <- pl + ggplot2::scale_fill_gradient2(low = cell_color_gradient[[1]],
                                                 mid = cell_color_gradient[[2]],
                                                 high = cell_color_gradient[[3]],
                                                 midpoint = gradient_midpoint)

      }
    }
  }
  return(pl)
}


#' @title plot_spat_point_layer_ggplot_noFILL
#' @name plot_spat_point_layer_ggplot_noFILL
#' @description creat ggplot point layer for spatial coordinates without borders
#' @param gobject giotto object
#' @param sdimx x-axis dimension name (default = 'sdimx')
#' @param sdimy y-axis dimension name (default = 'sdimy')
#' @param cell_locations_metadata_selected annotated location from selected cells
#' @param cell_locations_metadata_other annotated location from non-selected cells
#' @param cell_color color for cells (see details)
#' @param color_as_factor convert color column to factor
#' @param cell_color_code named vector with colors
#' @param cell_color_gradient vector with 3 colors for numeric data
#' @param gradient_midpoint midpoint for color gradient
#' @param gradient_limits vector with lower and upper limits
#' @param select_cell_groups select subset of cells/clusters based on cell_color parameter
#' @param select_cells select subset of cells based on cell IDs
#' @param point_size size of point (cell)
#' @param show_cluster_center plot center of selected clusters
#' @param show_center_label plot label of selected clusters
#' @param center_point_size size of center points
#' @param label_size  size of labels
#' @param label_fontface font of labels
#' @param show_other_cells display not selected cells
#' @param other_cell_color color for not selected cells
#' @param other_point_size point size for not selected cells
#' @param show_legend show legend
#' @return ggplot
#' @details Description of parameters.
#' @export
#' @examples
#'     plot_spat_point_layer_ggplot_noFILL(gobject)
plot_spat_point_layer_ggplot_noFILL = function(ggobject,
                                               sdimx = NULL,
                                               sdimy = NULL,
                                               cell_locations_metadata_selected,
                                               cell_locations_metadata_other,
                                               cell_color = NULL,
                                               color_as_factor = T,
                                               cell_color_code = NULL,
                                               cell_color_gradient = c('blue', 'white', 'red'),
                                               gradient_midpoint = NULL,
                                               gradient_limits = NULL,
                                               select_cell_groups = NULL,
                                               select_cells = NULL,
                                               point_size = 2,
                                               show_cluster_center = F,
                                               show_center_label = T,
                                               center_point_size = 4,
                                               label_size = 4,
                                               label_fontface = 'bold',
                                               show_other_cells = T,
                                               other_cell_color = 'lightgrey',
                                               other_point_size = 1,
                                               show_legend = TRUE

) {

  ## specify spatial dimensions first
  if(is.null(sdimx) | is.null(sdimy)) {

    warning("plot_method = ggplot, but spatial dimensions for sdimx and/or sdimy are not specified. \n
            It will default to the 'sdimx' and 'sdimy' ")
    sdimx = 'sdimx'
    sdimy = 'sdimy'
  }

  ## ggplot object
  pl = ggobject

  ## first plot other non-selected cells
  if((!is.null(select_cells) | !is.null(select_cell_groups)) & show_other_cells == TRUE) {
    pl <- pl + ggplot2::geom_point(data = cell_locations_metadata_other, aes_string(x = sdimx, sdimy),
                                   color = other_cell_color, show.legend = F, size = other_point_size)
  }


  ## order of color
  # 1. if NULL then default to lightblue
  # 2. if character vector
  # 2.1 if length of cell_color is longer than 1 and has colors
  # 2.2 if not part of metadata then suppose its color
  # 2.3 part of metadata
  # 2.3.1 numerical column
  # 2.3.2 factor column or character to factor


  # cell color default
  if(is.null(cell_color)) {

    cell_color = 'lightblue'
    pl <- pl + ggplot2::geom_point(data = cell_locations_metadata_selected,
                                   aes_string(x = sdimx, y = sdimy),
                                   show.legend = show_legend, shape = 19,
                                   color = cell_color, size = point_size)


  } else if(length(cell_color) > 1) {

    if(is.numeric(cell_color) | is.factor(cell_color)) {
      if(nrow(cell_locations_metadata_selected) != length(cell_color)) stop('\n vector needs to be the same lengths as number of cells \n')
      cell_locations_metadata_selected[['temp_color']] = cell_color

      pl <- pl + ggplot2::geom_point(data = cell_locations_metadata_selected, aes_string2(x = sdimx, y = sdimy, color = 'temp_color'),
                                     show.legend = show_legend, shape = 19, size = point_size)

    } else if(is.character(cell_color)) {
      if(!all(cell_color %in% grDevices::colors())) stop('cell_color is not numeric, a factor or vector of colors \n')
      pl <- pl + ggplot2::geom_point(data = cell_locations_metadata_selected, aes_string2(x = sdimx, y = sdimy),
                                     show.legend = show_legend, shape = 19, color = cell_color, size = point_size)

    }

  } else if(is.character(cell_color)) {
    if(!cell_color %in% colnames(cell_locations_metadata_selected)) {
      if(!cell_color %in% grDevices::colors()) stop(cell_color,' is not a color or a column name \n')
      pl <- pl + ggplot2::geom_point(data = cell_locations_metadata_selected,
                                     aes_string2(x = sdimx, y = sdimy),
                                     show.legend = show_legend, shape = 19, color = cell_color, size = point_size)

    } else {

      class_cell_color = class(cell_locations_metadata_selected[[cell_color]])

      if((class_cell_color == 'integer' | class_cell_color == 'numeric') & color_as_factor == FALSE) {
        # set upper and lower limits
        if(!is.null(gradient_limits) & is.vector(gradient_limits) & length(gradient_limits) == 2) {
          lower_lim = gradient_limits[[1]]
          upper_lim = gradient_limits[[2]]

          numeric_data = cell_locations_metadata_selected[[cell_color]]
          limit_numeric_data = ifelse(numeric_data > upper_lim, upper_lim,
                                      ifelse(numeric_data < lower_lim, lower_lim, numeric_data))
          cell_locations_metadata_selected[[cell_color]] = limit_numeric_data
        }

        pl <- pl + ggplot2::geom_point(data = cell_locations_metadata_selected,
                                       aes_string2(x = sdimx, y = sdimy, color = cell_color),
                                       show.legend = show_legend, shape = 19, size = point_size)



      } else {

        # convert character or numeric to factor
        if(color_as_factor == TRUE) {
          factor_data = factor(cell_locations_metadata_selected[[cell_color]])
          cell_locations_metadata_selected[[cell_color]] <- factor_data
        }

        # if you want to show centers or labels then calculate centers
        if(show_cluster_center == TRUE | show_center_label == TRUE) {
          annotated_DT_centers = cell_locations_metadata_selected[, .(center_1 = median(get('sdimx')), center_2 = median(get('sdimy'))), by = cell_color]
          factor_center_data = factor(annotated_DT_centers[[cell_color]])
          annotated_DT_centers[[cell_color]] <- factor_center_data
        }

        pl <- pl + ggplot2::geom_point(data = cell_locations_metadata_selected,
                                       aes_string2(x = sdimx, y = sdimy, color = cell_color),
                                       show.legend = show_legend, shape = 19, size = point_size)


        ## plot centers
        if(show_cluster_center == TRUE & (color_as_factor == TRUE | class_cell_color %in% c('character', 'factor'))) {

          pl <- pl + ggplot2::geom_point(data = annotated_DT_centers,
                                         aes_string2(x = 'center_1', y = 'center_2', color = cell_color),
                                         size = center_point_size, shape = 19)
        }

        ## plot labels
        if(show_center_label == TRUE) {
          pl <- pl + ggrepel::geom_text_repel(data = annotated_DT_centers,
                                              aes_string2(x = 'center_1', y = 'center_2', label = cell_color),
                                              size = label_size, fontface = label_fontface)
        }

      }

      ## specificy colors to use
      if(!is.null(cell_color_code)) {

        pl <- pl + ggplot2::scale_color_manual(values = cell_color_code)

      } else if(color_as_factor == T) {

        number_colors = length(unique(factor_data))
        cell_color_code = Giotto:::getDistinctColors(n = number_colors)
        names(cell_color_code) = unique(factor_data)
        pl <- pl + ggplot2::scale_color_manual(values = cell_color_code)

      } else if(color_as_factor == F){

        if(is.null(gradient_midpoint)) {
          gradient_midpoint = median(cell_locations_metadata_selected[[cell_color]])
        }

        pl <- pl + ggplot2::scale_color_gradient2(low = cell_color_gradient[[1]],
                                                  mid = cell_color_gradient[[2]],
                                                  high = cell_color_gradient[[3]],
                                                  midpoint = gradient_midpoint)

      }
    }
  }
  return(pl)
}






#' @title spatPlot2D_single
#' @name spatPlot2D_single
#' @description Visualize cells according to spatial coordinates
#' @param gobject giotto object
#' @param sdimx x-axis dimension name (default = 'sdimx')
#' @param sdimy y-axis dimension name (default = 'sdimy')
#' @param spat_enr_names names of spatial enrichment results to include
#' @param cell_color color for cells (see details)
#' @param color_as_factor convert color column to factor
#' @param cell_color_code named vector with colors
#' @param cell_color_gradient vector with 3 colors for numeric data
#' @param gradient_midpoint midpoint for color gradient
#' @param gradient_limits vector with lower and upper limits
#' @param select_cell_groups select subset of cells/clusters based on cell_color parameter
#' @param select_cells select subset of cells based on cell IDs
#' @param point_shape point with border or not (border or no_border)
#' @param point_size size of point (cell)
#' @param point_border_col color of border around points
#' @param point_border_stroke stroke size of border around points
#' @param show_cluster_center plot center of selected clusters
#' @param show_center_label plot label of selected clusters
#' @param center_point_size size of center points
#' @param label_size  size of labels
#' @param label_fontface font of labels
#' @param show_network show underlying spatial network
#' @param spatial_network_name name of spatial network to use
#' @param network_color color of spatial network
#' @param network_alpha alpha of spatial network
#' @param show_grid show spatial grid
#' @param spatial_grid_name name of spatial grid to use
#' @param grid_color color of spatial grid
#' @param show_other_cells display not selected cells
#' @param other_cell_color color of not selected cells
#' @param other_point_size point size of not selected cells
#' @param other_cells_alpha alpha of not selected cells
#' @param coord_fix_ratio fix ratio between x and y-axis
#' @param title title of plot
#' @param show_legend show legend
#' @param legend_text size of legend text
#' @param legend_symbol_size size of legend symbols
#' @param background_color color of plot background
#' @param axis_text size of axis text
#' @param axis_title size of axis title
#' @param show_plot show plot
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters from \code{\link{all_plots_save_function}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return ggplot
#' @details Description of parameters.
#' @export
#' @seealso \code{\link{spatPlot3D}}
#' @examples
#'     spatPlot2D_single(gobject)
spatPlot2D_single = function(gobject,
                             sdimx = 'sdimx',
                             sdimy = 'sdimy',
                             spat_enr_names = NULL,
                             cell_color = NULL,
                             color_as_factor = T,
                             cell_color_code = NULL,
                             cell_color_gradient = c('blue', 'white', 'red'),
                             gradient_midpoint = NULL,
                             gradient_limits = NULL,
                             select_cell_groups = NULL,
                             select_cells = NULL,
                             point_shape = c('border', 'no_border'),
                             point_size = 3,
                             point_border_col = 'black',
                             point_border_stroke = 0.1,
                             show_cluster_center = F,
                             show_center_label = F,
                             center_point_size = 4,
                             center_point_border_col = 'black',
                             center_point_border_stroke = 0.1,
                             label_size = 4,
                             label_fontface = 'bold',
                             show_network = F,
                             spatial_network_name = NULL,
                             network_color = NULL,
                             network_alpha = 1,
                             show_grid = F,
                             spatial_grid_name = 'spatial_grid',
                             grid_color = NULL,
                             show_other_cells = T,
                             other_cell_color = 'lightgrey',
                             other_point_size = 1,
                             other_cells_alpha = 0.1,
                             coord_fix_ratio = NULL,
                             title = NULL,
                             show_legend = T,
                             legend_text = 8,
                             legend_symbol_size = 1,
                             background_color = 'white',
                             axis_text = 8,
                             axis_title = 8,
                             show_plot = NA,
                             return_plot = NA,
                             save_plot = NA,
                             save_param =  list(),
                             default_save_name = 'spatPlot2D_single'
) {


  ## point shape ##
  point_shape = match.arg(point_shape, choices = c('border', 'no_border'))

  ## get spatial cell locations
  cell_locations  = gobject@spatial_locs

  ## extract spatial network
  if(!is.null(spatial_network_name)) {
    spatial_network = select_spatialNetwork(gobject, name = spatial_network_name, return_network_Obj = FALSE)
  } else {
    spatial_network = NULL
  }

  ## extract spatial grid
  if(!is.null(spatial_grid_name)) {
    spatial_grid    = gobject@spatial_grid[[spatial_grid_name]]
  } else {
    spatial_grid = NULL
  }

  ## get cell metadata
  cell_metadata = combineMetadata(gobject = gobject,
                                  spat_enr_names = spat_enr_names)
  #cell_metadata   = cell_metadata[, !grepl('cell_ID', colnames(cell_metadata)), with = F]

  if(nrow(cell_metadata) == 0) {
    cell_locations_metadata = cell_locations
  } else {
    #cell_locations_metadata <- cbind(cell_locations, cell_metadata)
    cell_locations_metadata <- cell_metadata
  }

  ## create subsets if needed
  if(!is.null(select_cells) & !is.null(select_cell_groups)) {
    cat('You have selected both individual cell IDs and a group of cells \n')
    group_cell_IDs = cell_locations_metadata[get(cell_color) %in% select_cell_groups][['cell_ID']]
    select_cells = unique(c(select_cells, group_cell_IDs))
  } else if(!is.null(select_cell_groups)) {
    select_cells = cell_locations_metadata[get(cell_color) %in% select_cell_groups][['cell_ID']]
  }

  if(!is.null(select_cells)) {
    cell_locations_metadata_other = cell_locations_metadata[!cell_locations_metadata$cell_ID %in% select_cells]
    cell_locations_metadata_selected = cell_locations_metadata[cell_locations_metadata$cell_ID %in% select_cells]
    spatial_network <- spatial_network[spatial_network$to %in% select_cells & spatial_network$from %in% select_cells]

    # if specific cells are selected
    # cell_locations_metadata = cell_locations_metadata_selected

  } else if(is.null(select_cells)) {

    cell_locations_metadata_selected = cell_locations_metadata
    cell_locations_metadata_other = NULL

  }



  ### create 2D plot with ggplot ###
  #cat('create 2D plot with ggplot \n')


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

  ## plot point layer
  if(point_shape == 'border') {
    pl = plot_spat_point_layer_ggplot(ggobject = pl,
                                      sdimx = sdimx,
                                      sdimy = sdimy,
                                      cell_locations_metadata_selected = cell_locations_metadata_selected,
                                      cell_locations_metadata_other = cell_locations_metadata_other,
                                      cell_color = cell_color,
                                      color_as_factor = color_as_factor,
                                      cell_color_code = cell_color_code,
                                      cell_color_gradient = cell_color_gradient,
                                      gradient_midpoint = gradient_midpoint,
                                      gradient_limits = gradient_limits,
                                      select_cell_groups = select_cell_groups,
                                      select_cells = select_cells,
                                      point_size = point_size,
                                      point_border_stroke = point_border_stroke,
                                      point_border_col = point_border_col,
                                      show_cluster_center = show_cluster_center,
                                      show_center_label = show_center_label,
                                      center_point_size = center_point_size,
                                      center_point_border_col = center_point_border_col,
                                      center_point_border_stroke = center_point_border_stroke,
                                      label_size = label_size,
                                      label_fontface = label_fontface,
                                      show_other_cells = show_other_cells,
                                      other_cell_color = other_cell_color,
                                      other_point_size = other_point_size,
                                      show_legend = show_legend)
  } else if(point_shape == 'no_border') {
    pl = plot_spat_point_layer_ggplot_noFILL(ggobject = pl,
                                             sdimx = sdimx,
                                             sdimy = sdimy,
                                             cell_locations_metadata_selected = cell_locations_metadata_selected,
                                             cell_locations_metadata_other = cell_locations_metadata_other,
                                             cell_color = cell_color,
                                             color_as_factor = color_as_factor,
                                             cell_color_code = cell_color_code,
                                             cell_color_gradient = cell_color_gradient,
                                             gradient_midpoint = gradient_midpoint,
                                             gradient_limits = gradient_limits,
                                             select_cell_groups = select_cell_groups,
                                             select_cells = select_cells,
                                             point_size = point_size,
                                             show_cluster_center = show_cluster_center,
                                             show_center_label = show_center_label,
                                             center_point_size = center_point_size,
                                             label_size = label_size,
                                             label_fontface = label_fontface,
                                             show_other_cells = show_other_cells,
                                             other_cell_color = other_cell_color,
                                             other_point_size = other_point_size,
                                             show_legend = show_legend)

  }



  ## adjust theme settings
  pl <- pl + ggplot2::theme(plot.title = element_text(hjust = 0.5),
                            legend.title = element_blank(),
                            legend.text = element_text(size = legend_text),
                            axis.title = element_text(size = axis_title),
                            axis.text = element_text(size = axis_text),
                            panel.grid = element_blank(),
                            panel.background = element_rect(fill = background_color))

  ## change symbol size of legend
  if(color_as_factor == TRUE) {
    if(point_shape == 'border') {
      pl = pl + guides(fill = guide_legend(override.aes = list(size = legend_symbol_size)))
    } else if(point_shape == 'no_border') {
      pl = pl + guides(color = guide_legend(override.aes = list(size = legend_symbol_size)))
    }
  }


  # fix coord ratio
  if(!is.null(coord_fix_ratio)) {
    pl <- pl + ggplot2::coord_fixed(ratio = coord_fix_ratio)
  }

  # provide x, y and plot titles
  if(is.null(title)) title = cell_color
  pl <- pl + ggplot2::labs(x = 'x coordinates', y = 'y coordinates', title = title)


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






#' @title spatPlot2D
#' @name spatPlot2D
#' @description Visualize cells according to spatial coordinates
#' @param gobject giotto object
#' @param groub_by create multiple plots based on cell annotation column
#' @param group_by_subset subset the group_by factor column
#' @param sdimx x-axis dimension name (default = 'sdimx')
#' @param sdimy y-axis dimension name (default = 'sdimy')
#' @param spat_enr_names names of spatial enrichment results to include
#' @param cell_color color for cells (see details)
#' @param color_as_factor convert color column to factor
#' @param cell_color_code named vector with colors
#' @param cell_color_gradient vector with 3 colors for numeric data
#' @param gradient_midpoint midpoint for color gradient
#' @param gradient_limits vector with lower and upper limits
#' @param select_cell_groups select subset of cells/clusters based on cell_color parameter
#' @param select_cells select subset of cells based on cell IDs
#' @param point_shape point with border or not (border or no_border)
#' @param point_size size of point (cell)
#' @param point_border_col color of border around points
#' @param point_border_stroke stroke size of border around points
#' @param show_cluster_center plot center of selected clusters
#' @param show_center_label plot label of selected clusters
#' @param center_point_size size of center points
#' @param label_size  size of labels
#' @param label_fontface font of labels
#' @param show_network show underlying spatial network
#' @param spatial_network_name name of spatial network to use
#' @param network_color color of spatial network
#' @param network_alpha alpha of spatial network
#' @param show_grid show spatial grid
#' @param spatial_grid_name name of spatial grid to use
#' @param grid_color color of spatial grid
#' @param show_other_cells display not selected cells
#' @param other_cell_color color of not selected cells
#' @param other_point_size point size of not selected cells
#' @param other_cells_alpha alpha of not selected cells
#' @param coord_fix_ratio fix ratio between x and y-axis
#' @param title title of plot
#' @param show_legend show legend
#' @param legend_text size of legend text
#' @param legend_symbol_size size of legend symbols
#' @param background_color color of plot background
#' @param axis_text size of axis text
#' @param axis_title size of axis title
#' @param cow_n_col cowplot param: how many columns
#' @param cow_rel_h cowplot param: relative height
#' @param cow_rel_w cowplot param: relative width
#' @param cow_align cowplot param: how to align
#' @param show_plot show plot
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters from \code{\link{all_plots_save_function}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return ggplot
#' @details Description of parameters.
#' @export
#' @seealso \code{\link{spatPlot3D}}
#' @examples
#'     spatPlot2D(gobject)
spatPlot2D = function(gobject,
                      group_by = NULL,
                      group_by_subset = NULL,
                      sdimx = 'sdimx',
                      sdimy = 'sdimy',
                      spat_enr_names = NULL,
                      cell_color = NULL,
                      color_as_factor = T,
                      cell_color_code = NULL,
                      cell_color_gradient = c('blue', 'white', 'red'),
                      gradient_midpoint = NULL,
                      gradient_limits = NULL,
                      select_cell_groups = NULL,
                      select_cells = NULL,
                      point_shape = c('border', 'no_border'),
                      point_size = 3,
                      point_border_col = 'black',
                      point_border_stroke = 0.1,
                      show_cluster_center = F,
                      show_center_label = F,
                      center_point_size = 4,
                      center_point_border_col = 'black',
                      center_point_border_stroke = 0.1,
                      label_size = 4,
                      label_fontface = 'bold',
                      show_network = F,
                      spatial_network_name = NULL,
                      network_color = NULL,
                      network_alpha = 1,
                      show_grid = F,
                      spatial_grid_name = 'spatial_grid',
                      grid_color = NULL,
                      show_other_cells = T,
                      other_cell_color = 'lightgrey',
                      other_point_size = 1,
                      other_cells_alpha = 0.1,
                      coord_fix_ratio = NULL,
                      title = NULL,
                      show_legend = T,
                      legend_text = 8,
                      legend_symbol_size = 1,
                      background_color = 'white',
                      axis_text = 8,
                      axis_title = 8,
                      cow_n_col = 2,
                      cow_rel_h = 1,
                      cow_rel_w = 1,
                      cow_align = 'h',
                      show_plot = NA,
                      return_plot = NA,
                      save_plot = NA,
                      save_param =  list(),
                      default_save_name = 'spatPlot2D') {


  ## check group_by
  if(is.null(group_by)) {

    spatPlot2D_single(gobject = gobject,
                      sdimx = sdimx,
                      sdimy = sdimy,
                      spat_enr_names = spat_enr_names,
                      cell_color = cell_color,
                      color_as_factor = color_as_factor,
                      cell_color_code = cell_color_code,
                      cell_color_gradient = cell_color_gradient,
                      gradient_midpoint = gradient_midpoint,
                      gradient_limits = gradient_limits,
                      select_cell_groups = select_cell_groups,
                      select_cells = select_cells,
                      point_shape = point_shape,
                      point_size = point_size,
                      point_border_col = point_border_col,
                      point_border_stroke = point_border_stroke,
                      show_cluster_center = show_cluster_center,
                      show_center_label = show_center_label,
                      center_point_size = center_point_size,
                      center_point_border_col = center_point_border_col,
                      center_point_border_stroke = center_point_border_stroke,
                      label_size = label_size,
                      label_fontface = label_fontface,
                      show_network = show_network,
                      spatial_network_name = spatial_network_name,
                      network_color = network_color,
                      network_alpha = network_alpha,
                      show_grid = show_grid,
                      spatial_grid_name = spatial_grid_name,
                      grid_color = grid_color,
                      show_other_cells = show_other_cells,
                      other_cell_color = other_cell_color,
                      other_point_size = other_point_size,
                      other_cells_alpha = other_cells_alpha,
                      coord_fix_ratio = coord_fix_ratio,
                      show_legend = show_legend,
                      legend_text = legend_text,
                      legend_symbol_size = legend_symbol_size,
                      background_color = background_color,
                      axis_text = axis_text,
                      axis_title = axis_title,
                      title = title,
                      show_plot = show_plot,
                      return_plot = return_plot,
                      save_plot = save_plot,
                      save_param =  save_param,
                      default_save_name = default_save_name)


  } else {

    ## metadata
    comb_metadata = combineMetadata(gobject = gobject,
                                    spat_enr_names = spat_enr_names)
    possible_meta_groups = colnames(comb_metadata)

    ## check if group_by is found
    if(!group_by %in% possible_meta_groups) {
      stop("group_by ", group_by, " was not found in pDataDT()")
    }

    unique_groups = unique(comb_metadata[[group_by]])

    # subset unique_groups
    if(!is.null(group_by_subset)) {
      not_found = group_by_subset[!group_by_subset %in% unique_groups]
      if(length(not_found) > 0) {
        cat('the following subset was not found: ', not_found)
      }
      unique_groups = unique_groups[unique_groups %in% group_by_subset]
    }

    # create matching cell_color_code
    if(is.null(cell_color_code)) {
      if(is.character(cell_color)) {

        if(cell_color %in% colnames(comb_metadata)) {

          if(color_as_factor == TRUE) {
            number_colors = length(unique(comb_metadata[[cell_color]]))
            cell_color_code = Giotto:::getDistinctColors(n = number_colors)
            names(cell_color_code) = unique(comb_metadata[[cell_color]])
            cell_color_code = cell_color_code
          }
        }
      }
    }


    # print, return and save parameters
    show_plot = ifelse(is.na(show_plot), readGiottoInstructions(gobject, param = 'show_plot'), show_plot)
    save_plot = ifelse(is.na(save_plot), readGiottoInstructions(gobject, param = 'save_plot'), save_plot)
    return_plot = ifelse(is.na(return_plot), readGiottoInstructions(gobject, param = 'return_plot'), return_plot)

    ## plotting ##
    savelist <- list()


    for(group_id in 1:length(unique_groups)) {

      group = unique_groups[group_id]

      subset_cell_IDs = comb_metadata[get(group_by) == group][['cell_ID']]
      temp_gobject = subsetGiotto(gobject = gobject, cell_ids = subset_cell_IDs)

      pl = spatPlot2D_single(gobject = temp_gobject,
                             sdimx = sdimx,
                             sdimy = sdimy,
                             spat_enr_names = spat_enr_names,
                             cell_color = cell_color,
                             cell_color_code = cell_color_code,
                             color_as_factor = color_as_factor,
                             cell_color_gradient = cell_color_gradient,
                             gradient_midpoint = gradient_midpoint,
                             gradient_limits = gradient_limits,
                             select_cell_groups = select_cell_groups,
                             select_cells = select_cells,
                             point_shape = point_shape,
                             point_size = point_size,
                             point_border_col = point_border_col,
                             point_border_stroke = point_border_stroke,
                             show_cluster_center = show_cluster_center,
                             show_center_label = show_center_label,
                             center_point_size = center_point_size,
                             center_point_border_col = center_point_border_col,
                             center_point_border_stroke = center_point_border_stroke,
                             label_size = label_size,
                             label_fontface = label_fontface,
                             show_network = show_network,
                             spatial_network_name = spatial_network_name,
                             network_color = network_color,
                             network_alpha = network_alpha,
                             show_grid = show_grid,
                             spatial_grid_name = spatial_grid_name,
                             grid_color = grid_color,
                             show_other_cells = show_other_cells,
                             other_cell_color = other_cell_color,
                             other_point_size = other_point_size,
                             other_cells_alpha = other_cells_alpha,
                             coord_fix_ratio = coord_fix_ratio,
                             title = group,
                             show_legend = show_legend,
                             legend_text = legend_text,
                             legend_symbol_size = legend_symbol_size,
                             background_color = background_color,
                             axis_text = axis_text,
                             axis_title = axis_title,
                             show_plot = FALSE,
                             return_plot = TRUE,
                             save_plot = FALSE,
                             save_param =  list(),
                             default_save_name = 'spatPlot2D')


      savelist[[group_id]] <- pl

    }

    # combine plots with cowplot
    combo_plot <- cowplot::plot_grid(plotlist = savelist,
                                     ncol = cow_n_col,
                                     rel_heights = cow_rel_h,
                                     rel_widths = cow_rel_w,
                                     align = cow_align)


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

}






#' @title spatPlot
#' @name spatPlot
#' @description Visualize cells according to spatial coordinates
#' @param gobject giotto object
#' @param groub_by create multiple plots based on cell annotation column
#' @param group_by_subset subset the group_by factor column
#' @param sdimx x-axis dimension name (default = 'sdimx')
#' @param sdimy y-axis dimension name (default = 'sdimy')
#' @param spat_enr_names names of spatial enrichment results to include
#' @param cell_color color for cells (see details)
#' @param color_as_factor convert color column to factor
#' @param cell_color_code named vector with colors
#' @param cell_color_gradient vector with 3 colors for numeric data
#' @param gradient_midpoint midpoint for color gradient
#' @param gradient_limits vector with lower and upper limits
#' @param select_cell_groups select subset of cells/clusters based on cell_color parameter
#' @param select_cells select subset of cells based on cell IDs
#' @param point_shape point with border or not (border or no_border)
#' @param point_size size of point (cell)
#' @param point_border_col color of border around points
#' @param point_border_stroke stroke size of border around points
#' @param show_cluster_center plot center of selected clusters
#' @param show_center_label plot label of selected clusters
#' @param center_point_size size of center points
#' @param label_size  size of labels
#' @param label_fontface font of labels
#' @param show_network show underlying spatial network
#' @param spatial_network_name name of spatial network to use
#' @param network_color color of spatial network
#' @param network_alpha alpha of spatial network
#' @param show_grid show spatial grid
#' @param spatial_grid_name name of spatial grid to use
#' @param grid_color color of spatial grid
#' @param show_other_cells display not selected cells
#' @param other_cell_color color of not selected cells
#' @param other_point_size point size of not selected cells
#' @param other_cells_alpha alpha of not selected cells
#' @param coord_fix_ratio fix ratio between x and y-axis
#' @param title title of plot
#' @param show_legend show legend
#' @param legend_text size of legend text
#' @param legend_symbol_size size of legend symbols
#' @param background_color color of plot background
#' @param axis_text size of axis text
#' @param axis_title size of axis title
#' @param cow_n_col cowplot param: how many columns
#' @param cow_rel_h cowplot param: relative height
#' @param cow_rel_w cowplot param: relative width
#' @param cow_align cowplot param: how to align
#' @param show_plot show plot
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters from \code{\link{all_plots_save_function}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return ggplot
#' @details Description of parameters.
#' @export
#' @seealso \code{\link{spatPlot3D}}
#' @examples
#'     spatPlot(gobject)
spatPlot = function(gobject,
                    group_by = NULL,
                    group_by_subset = NULL,
                    sdimx = 'sdimx',
                    sdimy = 'sdimy',
                    spat_enr_names = NULL,
                    cell_color = NULL,
                    color_as_factor = T,
                    cell_color_code = NULL,
                    cell_color_gradient = c('blue', 'white', 'red'),
                    gradient_midpoint = NULL,
                    gradient_limits = NULL,
                    select_cell_groups = NULL,
                    select_cells = NULL,
                    point_shape = c('border', 'no_border'),
                    point_size = 3,
                    point_border_col = 'black',
                    point_border_stroke = 0.1,
                    show_cluster_center = F,
                    show_center_label = F,
                    center_point_size = 4,
                    center_point_border_col = 'black',
                    center_point_border_stroke = 0.1,
                    label_size = 4,
                    label_fontface = 'bold',
                    show_network = F,
                    spatial_network_name = NULL,
                    network_color = NULL,
                    network_alpha = 1,
                    show_grid = F,
                    spatial_grid_name = 'spatial_grid',
                    grid_color = NULL,
                    show_other_cells = T,
                    other_cell_color = 'lightgrey',
                    other_point_size = 1,
                    other_cells_alpha = 0.1,
                    coord_fix_ratio = NULL,
                    title = NULL,
                    show_legend = T,
                    legend_text = 8,
                    legend_symbol_size = 1,
                    background_color = 'white',
                    axis_text = 8,
                    axis_title = 8,
                    cow_n_col = 2,
                    cow_rel_h = 1,
                    cow_rel_w = 1,
                    cow_align = 'h',
                    show_plot = NA,
                    return_plot = NA,
                    save_plot = NA,
                    save_param =  list(),
                    default_save_name = 'spatPlot') {

  spatPlot2D(gobject = gobject,
             group_by = group_by,
             group_by_subset = group_by_subset,
             sdimx = sdimx,
             sdimy = sdimy,
             spat_enr_names = spat_enr_names,
             cell_color = cell_color,
             color_as_factor = color_as_factor,
             cell_color_code = cell_color_code,
             cell_color_gradient = cell_color_gradient,
             gradient_midpoint = gradient_midpoint,
             gradient_limits = gradient_limits,
             select_cell_groups = select_cell_groups,
             select_cells = select_cells,
             point_shape = point_shape,
             point_size = point_size,
             point_border_col = point_border_col,
             point_border_stroke = point_border_stroke,
             show_cluster_center = show_cluster_center,
             show_center_label = show_center_label,
             center_point_size = center_point_size,
             center_point_border_col = center_point_border_col,
             center_point_border_stroke = center_point_border_stroke,
             label_size = label_size,
             label_fontface = label_fontface,
             show_network = show_network,
             spatial_network_name = spatial_network_name,
             network_color = network_color,
             network_alpha = network_alpha,
             show_grid = show_grid,
             spatial_grid_name = spatial_grid_name,
             grid_color = grid_color,
             show_other_cells = show_other_cells,
             other_cell_color = other_cell_color,
             other_point_size = other_point_size,
             other_cells_alpha = other_cells_alpha,
             coord_fix_ratio = coord_fix_ratio,
             title = title,
             show_legend = show_legend,
             legend_text = legend_text,
             legend_symbol_size = legend_symbol_size,
             background_color = background_color,
             axis_text = axis_text,
             axis_title = axis_title,
             cow_n_col = cow_n_col,
             cow_rel_h = cow_rel_h,
             cow_rel_w = cow_rel_w,
             cow_align = cow_align,
             show_plot = show_plot,
             return_plot = return_plot,
             save_plot = save_plot,
             save_param =  save_param,
             default_save_name = default_save_name)

}








#' @title spatDimPlot2D
#' @name spatDimPlot2D
#' @description Visualize cells according to spatial AND dimension reduction coordinates 2D
#' @param gobject giotto object
#' @param plot_alignment direction to align plot
#' @param dim_reduction_to_use dimension reduction to use
#' @param dim_reduction_name dimension reduction name
#' @param dim1_to_use dimension to use on x-axis
#' @param dim2_to_use dimension to use on y-axis
#' @param sdimx = spatial dimension to use on x-axis
#' @param sdimy = spatial dimension to use on y-axis
#' @param spat_enr_names names of spatial enrichment results to include
#' @param cell_color color for cells (see details)
#' @param color_as_factor convert color column to factor
#' @param cell_color_code named vector with colors
#' @param cell_color_gradient vector with 3 colors for numeric data
#' @param gradient_midpoint midpoint for color gradient
#' @param gradient_limits vector with lower and upper limits
#' @param select_cell_groups select subset of cells/clusters based on cell_color parameter
#' @param select_cells select subset of cells based on cell IDs
#' @param dim_point_shape point with border or not (border or no_border)
#' @param dim_point_size size of points in dim. reduction space
#' @param dim_point_border_col border color of points in dim. reduction space
#' @param dim_point_border_stroke border stroke of points in dim. reduction space
#' @param spat_point_shape point with border or not (border or no_border)
#' @param spat_point_size size of spatial points
#' @param spat_point_border_col border color of spatial points
#' @param spat_point_border_stroke border stroke of spatial points
#' @param dim_show_cluster_center show the center of each cluster
#' @param dim_show_center_label provide a label for each cluster
#' @param dim_center_point_size size of the center point
#' @param dim_center_point_border_col border color of center point
#' @param dim_center_point_border_stroke stroke size of center point
#' @param dim_label_size size of the center label
#' @param dim_label_fontface font of the center label
#' @param spat_show_cluster_center show the center of each cluster
#' @param spat_show_center_label provide a label for each cluster
#' @param spat_center_point_size size of the center point
#' @param spat_label_size size of the center label
#' @param spat_label_fontface font of the center label
#' @param show_NN_network show underlying NN network
#' @param nn_network_to_use type of NN network to use (kNN vs sNN)
#' @param network_name name of NN network to use, if show_NN_network = TRUE
#' @param nn_network_alpha column to use for alpha of the edges
#' @param show_spatial_network show spatial network
#' @param spat_network_name name of spatial network to use
#' @param spat_network_color color of spatial network
#' @param show_spatial_grid show spatial grid
#' @param spat_grid_name name of spatial grid to use
#' @param spat_grid_color color of spatial grid
#' @param show_other_cells display not selected cells
#' @param other_cell_color color of not selected cells
#' @param dim_other_point_size size of not selected dim cells
#' @param spat_other_point_size size of not selected spat cells
#' @param spat_other_cells_alpha alpha of not selected spat cells
#' @param dim_show_legend show legend of dimension reduction plot
#' @param spat_show_legend show legend of spatial plot
#' @param legend_text size of legend text
#' @param legend_symbol_size size of legend symbols
#' @param dim_background_color background color of points in dim. reduction space
#' @param spat_background_color background color of spatial points
#' @param axis_text size of axis text
#' @param axis_title size of axis title
#' @param show_plot show plot
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters from \code{\link{all_plots_save_function}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return ggplot
#' @details Description of parameters.
#' @export
#' @seealso \code{\link{spatDimPlot3D}}
#' @examples
#'     spatDimPlot2D(gobject)
spatDimPlot2D <- function(gobject,
                          plot_alignment = c('vertical', 'horizontal'),
                          dim_reduction_to_use = 'umap',
                          dim_reduction_name = 'umap',
                          dim1_to_use = 1,
                          dim2_to_use = 2,
                          sdimx = 'sdimx',
                          sdimy = 'sdimy',
                          spat_enr_names = NULL,
                          cell_color = NULL,
                          color_as_factor = T,
                          cell_color_code = NULL,
                          cell_color_gradient = c('blue', 'white', 'red'),
                          gradient_midpoint = NULL,
                          gradient_limits = NULL,
                          select_cell_groups = NULL,
                          select_cells = NULL,
                          dim_point_shape = c('border', 'no_border'),
                          dim_point_size = 1,
                          dim_point_border_col = 'black',
                          dim_point_border_stroke = 0.1,
                          spat_point_shape = c('border', 'no_border'),
                          spat_point_size = 1,
                          spat_point_border_col = 'black',
                          spat_point_border_stroke = 0.1,
                          dim_show_cluster_center = F,
                          dim_show_center_label = T,
                          dim_center_point_size = 4,
                          dim_center_point_border_col = 'black',
                          dim_center_point_border_stroke = 0.1,
                          dim_label_size = 4,
                          dim_label_fontface = 'bold',
                          spat_show_cluster_center = F,
                          spat_show_center_label = F,
                          spat_center_point_size = 4,
                          spat_label_size = 4,
                          spat_label_fontface = 'bold',
                          show_NN_network = F,
                          nn_network_to_use = 'sNN',
                          network_name = 'sNN.pca',
                          nn_network_alpha = 0.05,
                          show_spatial_network = F,
                          spat_network_name = 'spatial_network',
                          spat_network_color = 'blue',
                          spat_network_alpha = 0.5,
                          show_spatial_grid = F,
                          spat_grid_name = 'spatial_grid',
                          spat_grid_color = 'blue',
                          show_other_cells = T,
                          other_cell_color = 'lightgrey',
                          dim_other_point_size = 1,
                          spat_other_point_size = 1,
                          spat_other_cells_alpha = 0.5,
                          dim_show_legend = F,
                          spat_show_legend = F,
                          legend_text = 8,
                          legend_symbol_size = 1,
                          dim_background_color = 'white',
                          spat_background_color = 'white',
                          axis_text = 8,
                          axis_title = 8,
                          show_plot = NA,
                          return_plot = NA,
                          save_plot = NA,
                          save_param =  list(),
                          default_save_name = 'spatDimPlot2D'
){

  plot_alignment = match.arg(plot_alignment, choices = c( 'vertical','horizontal'))


  # create matching cell_color_code
  if(is.null(cell_color_code)) {
    if(is.character(cell_color)) {

      cell_metadata = pDataDT(gobject)
      if(cell_color %in% colnames(cell_metadata)) {

        if(color_as_factor == TRUE) {
          number_colors = length(unique(cell_metadata[[cell_color]]))
          cell_color_code = Giotto:::getDistinctColors(n = number_colors)
          names(cell_color_code) = unique(cell_metadata[[cell_color]])
          cell_color_code = cell_color_code
        }
      }
    }
  }

  # dimension reduction plot
  dmpl = dimPlot2D(gobject = gobject,
                   group_by = NULL,
                   group_by_subset = NULL,
                   dim_reduction_to_use = dim_reduction_to_use,
                   dim_reduction_name = dim_reduction_name,
                   dim1_to_use = dim1_to_use,
                   dim2_to_use = dim2_to_use,
                   spat_enr_names = spat_enr_names,
                   cell_color = cell_color,
                   color_as_factor = color_as_factor,
                   cell_color_code = cell_color_code,
                   cell_color_gradient = cell_color_gradient,
                   gradient_midpoint = gradient_midpoint,
                   gradient_limits = gradient_limits,
                   select_cell_groups = select_cell_groups,
                   select_cells = select_cells,
                   point_shape = dim_point_shape,
                   point_size = dim_point_size,
                   point_border_col = dim_point_border_col,
                   point_border_stroke = dim_point_border_stroke,
                   show_cluster_center = dim_show_cluster_center,
                   show_center_label = dim_show_center_label,
                   center_point_size = dim_center_point_size,
                   center_point_border_col = dim_center_point_border_col,
                   center_point_border_stroke = dim_center_point_border_stroke,
                   label_size = dim_label_size,
                   label_fontface = dim_label_fontface,
                   show_NN_network = show_NN_network,
                   nn_network_to_use = nn_network_to_use,
                   network_name = network_name,
                   edge_alpha = nn_network_alpha,
                   show_other_cells = show_other_cells,
                   other_cell_color = other_cell_color,
                   other_point_size = dim_other_point_size,
                   show_legend = dim_show_legend,
                   legend_text = legend_text,
                   legend_symbol_size = legend_symbol_size,
                   background_color = dim_background_color,
                   axis_text = axis_text,
                   axis_title = axis_title,
                   show_plot = FALSE,
                   return_plot = TRUE,
                   save_plot = FALSE
  )

  # spatial plot
  spl = spatPlot2D(gobject = gobject,
                   group_by = NULL,
                   group_by_subset = NULL,
                   sdimx = sdimx,
                   sdimy = sdimy,
                   spat_enr_names = spat_enr_names,
                   cell_color = cell_color,
                   cell_color_code = cell_color_code,
                   color_as_factor = color_as_factor,
                   cell_color_gradient = cell_color_gradient,
                   gradient_midpoint = gradient_midpoint,
                   gradient_limits = gradient_limits,
                   select_cell_groups = select_cell_groups,
                   select_cells = select_cells,
                   point_shape = spat_point_shape,
                   point_size = spat_point_size,
                   point_border_col = spat_point_border_col,
                   point_border_stroke = spat_point_border_stroke,
                   show_cluster_center = spat_show_cluster_center,
                   show_center_label = spat_show_center_label,
                   center_point_size = spat_center_point_size,
                   center_point_border_col = spat_center_point_border_col,
                   center_point_border_stroke = spat_center_point_border_stroke,
                   label_size = spat_label_size,
                   label_fontface = spat_label_fontface,
                   show_network = show_spatial_network,
                   spatial_network_name = spat_network_name,
                   network_color = spat_network_color,
                   network_alpha = spat_network_alpha,
                   show_grid = show_spatial_grid,
                   spatial_grid_name = spat_grid_name,
                   grid_color = spat_grid_color,
                   show_other_cells = show_other_cells,
                   other_cell_color = other_cell_color,
                   other_point_size = spat_other_point_size,
                   other_cells_alpha = spat_other_cells_alpha,
                   coord_fix_ratio = NULL,
                   title = '',
                   show_legend = spat_show_legend,
                   legend_text = legend_text,
                   legend_symbol_size = legend_symbol_size,
                   background_color = spat_background_color,
                   axis_text = axis_text,
                   axis_title = axis_title,
                   show_plot = FALSE,
                   return_plot = TRUE,
                   save_plot = FALSE)


  if(plot_alignment == 'vertical') {
    ncol = 1
    nrow = 2
    combo_plot = cowplot::plot_grid(dmpl, spl, ncol = ncol, nrow = nrow, rel_heights = c(1), rel_widths = c(1), align = 'v')
  } else {
    ncol = 2
    nrow = 1
    combo_plot = cowplot::plot_grid(dmpl, spl, ncol = ncol, nrow = nrow, rel_heights = c(1), rel_widths = c(1), align = 'h')
  }


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




#' @title spatDimPlot
#' @name spatDimPlot
#' @description Visualize cells according to spatial AND dimension reduction coordinates 2D
#' @param gobject giotto object
#' @param plot_alignment direction to align plot
#' @param dim_reduction_to_use dimension reduction to use
#' @param dim_reduction_name dimension reduction name
#' @param dim1_to_use dimension to use on x-axis
#' @param dim2_to_use dimension to use on y-axis
#' @param sdimx = spatial dimension to use on x-axis
#' @param sdimy = spatial dimension to use on y-axis
#' @param spat_enr_names names of spatial enrichment results to include
#' @param cell_color color for cells (see details)
#' @param color_as_factor convert color column to factor
#' @param cell_color_code named vector with colors
#' @param cell_color_gradient vector with 3 colors for numeric data
#' @param gradient_midpoint midpoint for color gradient
#' @param gradient_limits vector with lower and upper limits
#' @param select_cell_groups select subset of cells/clusters based on cell_color parameter
#' @param select_cells select subset of cells based on cell IDs
#' @param dim_point_shape point with border or not (border or no_border)
#' @param dim_point_size size of points in dim. reduction space
#' @param dim_point_border_col border color of points in dim. reduction space
#' @param dim_point_border_stroke border stroke of points in dim. reduction space
#' @param spat_point_shape point with border or not (border or no_border)
#' @param spat_point_size size of spatial points
#' @param spat_point_border_col border color of spatial points
#' @param spat_point_border_stroke border stroke of spatial points
#' @param dim_show_cluster_center show the center of each cluster
#' @param dim_show_center_label provide a label for each cluster
#' @param dim_center_point_size size of the center point
#' @param dim_center_point_border_col border color of center point
#' @param dim_center_point_border_stroke stroke size of center point
#' @param dim_label_size size of the center label
#' @param dim_label_fontface font of the center label
#' @param spat_show_cluster_center show the center of each cluster
#' @param spat_show_center_label provide a label for each cluster
#' @param spat_center_point_size size of the center point
#' @param spat_label_size size of the center label
#' @param spat_label_fontface font of the center label
#' @param show_NN_network show underlying NN network
#' @param nn_network_to_use type of NN network to use (kNN vs sNN)
#' @param network_name name of NN network to use, if show_NN_network = TRUE
#' @param nn_network_alpha column to use for alpha of the edges
#' @param show_spatial_network show spatial network
#' @param spat_network_name name of spatial network to use
#' @param spat_network_color color of spatial network
#' @param show_spatial_grid show spatial grid
#' @param spat_grid_name name of spatial grid to use
#' @param spat_grid_color color of spatial grid
#' @param show_other_cells display not selected cells
#' @param other_cell_color color of not selected cells
#' @param dim_other_point_size size of not selected dim cells
#' @param spat_other_point_size size of not selected spat cells
#' @param spat_other_cells_alpha alpha of not selected spat cells
#' @param dim_show_legend show legend of dimension reduction plot
#' @param spat_show_legend show legend of spatial plot
#' @param legend_text size of legend text
#' @param legend_symbol_size size of legend symbols
#' @param dim_background_color background color of points in dim. reduction space
#' @param spat_background_color background color of spatial points
#' @param axis_text size of axis text
#' @param axis_title size of axis title
#' @param show_plot show plot
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters from \code{\link{all_plots_save_function}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return ggplot
#' @details Description of parameters.
#' @export
#' @seealso \code{\link{spatDimPlot2D}} and \code{\link{spatDimPlot3D}} for 3D visualization.
#' @examples
#'     spatDimPlot(gobject)
spatDimPlot = function(gobject,
                       plot_alignment = c('vertical', 'horizontal'),
                       dim_reduction_to_use = 'umap',
                       dim_reduction_name = 'umap',
                       dim1_to_use = 1,
                       dim2_to_use = 2,
                       sdimx = 'sdimx',
                       sdimy = 'sdimy',
                       spat_enr_names = NULL,
                       cell_color = NULL,
                       color_as_factor = T,
                       cell_color_code = NULL,
                       cell_color_gradient = c('blue', 'white', 'red'),
                       gradient_midpoint = NULL,
                       gradient_limits = NULL,
                       select_cell_groups = NULL,
                       select_cells = NULL,
                       dim_point_shape = c('border', 'no_border'),
                       dim_point_size = 1,
                       dim_point_border_col = 'black',
                       dim_point_border_stroke = 0.1,
                       spat_point_shape = c('border', 'no_border'),
                       spat_point_size = 1,
                       spat_point_border_col = 'black',
                       spat_point_border_stroke = 0.1,
                       dim_show_cluster_center = F,
                       dim_show_center_label = T,
                       dim_center_point_size = 4,
                       dim_center_point_border_col = 'black',
                       dim_center_point_border_stroke = 0.1,
                       dim_label_size = 4,
                       dim_label_fontface = 'bold',
                       spat_show_cluster_center = F,
                       spat_show_center_label = F,
                       spat_center_point_size = 4,
                       spat_label_size = 4,
                       spat_label_fontface = 'bold',
                       show_NN_network = F,
                       nn_network_to_use = 'sNN',
                       network_name = 'sNN.pca',
                       nn_network_alpha = 0.05,
                       show_spatial_network = F,
                       spat_network_name = 'spatial_network',
                       spat_network_color = 'blue',
                       spat_network_alpha = 0.5,
                       show_spatial_grid = F,
                       spat_grid_name = 'spatial_grid',
                       spat_grid_color = 'blue',
                       show_other_cells = T,
                       other_cell_color = 'lightgrey',
                       dim_other_point_size = 1,
                       spat_other_point_size = 1,
                       spat_other_cells_alpha = 0.5,
                       dim_show_legend = F,
                       spat_show_legend = F,
                       legend_text = 8,
                       legend_symbol_size = 1,
                       dim_background_color = 'white',
                       spat_background_color = 'white',
                       axis_text = 8,
                       axis_title = 8,
                       show_plot = NA,
                       return_plot = NA,
                       save_plot = NA,
                       save_param =  list(),
                       default_save_name = 'spatDimPlot') {

  spatDimPlot2D(gobject = gobject,
                plot_alignment = plot_alignment,
                dim_reduction_to_use = dim_reduction_to_use,
                dim_reduction_name = dim_reduction_name,
                dim1_to_use = dim1_to_use,
                dim2_to_use = dim2_to_use,
                sdimx = sdimx,
                sdimy = sdimy,
                spat_enr_names = spat_enr_names,
                cell_color = cell_color,
                color_as_factor = color_as_factor,
                cell_color_code = cell_color_code,
                cell_color_gradient = cell_color_gradient,
                gradient_midpoint = gradient_midpoint,
                gradient_limits = gradient_limits,
                select_cell_groups = select_cell_groups,
                select_cells = select_cells,
                dim_point_shape = dim_point_shape,
                dim_point_size = dim_point_size,
                dim_point_border_col = dim_point_border_col,
                dim_point_border_stroke = dim_point_border_stroke,
                spat_point_shape = spat_point_shape,
                spat_point_size = spat_point_size,
                spat_point_border_col = spat_point_border_col,
                spat_point_border_stroke = spat_point_border_stroke,
                dim_show_cluster_center = dim_show_cluster_center,
                dim_show_center_label = dim_show_center_label,
                dim_center_point_size = dim_center_point_size,
                dim_center_point_border_col = dim_center_point_border_col,
                dim_center_point_border_stroke = dim_center_point_border_stroke,
                dim_label_size = dim_label_size,
                dim_label_fontface = dim_label_fontface,
                spat_show_cluster_center = spat_show_cluster_center,
                spat_show_center_label = spat_show_center_label,
                spat_center_point_size = spat_center_point_size,
                spat_label_size = spat_label_size,
                spat_label_fontface = spat_label_fontface,
                show_NN_network = show_NN_network,
                nn_network_to_use = nn_network_to_use,
                network_name = network_name,
                nn_network_alpha = nn_network_alpha,
                show_spatial_network = show_spatial_network,
                spat_network_name = spat_network_name,
                spat_network_color = spat_network_color,
                spat_network_alpha = spat_network_alpha,
                show_spatial_grid = show_spatial_grid,
                spat_grid_name = spat_grid_name,
                spat_grid_color = spat_grid_color,
                show_other_cells = show_other_cells,
                other_cell_color = other_cell_color,
                dim_other_point_size = dim_other_point_size,
                spat_other_point_size = spat_other_point_size,
                spat_other_cells_alpha = spat_other_cells_alpha,
                dim_show_legend = dim_show_legend,
                spat_show_legend = spat_show_legend,
                legend_text = legend_text,
                legend_symbol_size = legend_symbol_size,
                dim_background_color = dim_background_color,
                spat_background_color = spat_background_color,
                axis_text = axis_text,
                axis_title = axis_title,
                show_plot = show_plot,
                return_plot = return_plot,
                save_plot = save_plot,
                save_param =  save_param,
                default_save_name = default_save_name)

}



#' @title spatGenePlot2D
#' @name spatGenePlot2D
#' @description Visualize cells and gene expression according to spatial coordinates
#' @param gobject giotto object
#' @param expression_values gene expression values to use
#' @param genes genes to show
#' @param genes_high_color color represents high gene expression
#' @param genes_mid_color color represents middle gene expression
#' @param genes_low_color color represents low gene expression
#' @param show_network show underlying spatial network
#' @param network_color color of spatial network
#' @param spatial_network_name name of spatial network to use
#' @param show_grid show spatial grid
#' @param grid_color color of spatial grid
#' @param spatial_grid_name name of spatial grid to use
#' @param midpoint expression midpoint
#' @param scale_alpha_with_expression scale expression with ggplot alpha parameter
#' @param point_shape point with border or not (border or no_border)
#' @param point_size size of point (cell)
#' @param point_border_col color of border around points
#' @param point_border_stroke stroke size of border around points
#' @param cow_n_col cowplot param: how many columns
#' @param cow_rel_h cowplot param: relative height
#' @param cow_rel_w cowplot param: relative width
#' @param cow_align cowplot param: how to align
#' @param show_legend show legend
#' @param legend_text size of legend text
#' @param background_color color of plot background
#' @param axis_text size of axis text
#' @param axis_title size of axis title
#' @param show_plot show plots
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters from \code{\link{all_plots_save_function}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @param ... parameters for cowplot::save_plot()
#' @return ggplot
#' @details Description of parameters.
#' @export
#' @seealso \code{\link{spatGenePlot3D}}
#' @examples
#'     spatGenePlot2D(gobject)
spatGenePlot2D <- function(gobject,
                           expression_values = c('normalized', 'scaled', 'custom'),
                           genes,
                           genes_high_color = 'darkred',
                           genes_mid_color = "white",
                           genes_low_color = "darkblue",
                           show_network = F,
                           network_color = NULL,
                           spatial_network_name = 'Delaunay_network',
                           edge_alpha = NULL,
                           show_grid = F,
                           grid_color = NULL,
                           spatial_grid_name = 'spatial_grid',
                           midpoint = 0,
                           scale_alpha_with_expression = FALSE,
                           point_shape = c('border', 'no_border'),
                           point_size = 1,
                           point_border_col = 'black',
                           point_border_stroke = 0.1,
                           show_legend = T,
                           legend_text = 8,
                           background_color = 'white',
                           axis_text = 8,
                           axis_title = 8,
                           cow_n_col = 2,
                           cow_rel_h = 1,
                           cow_rel_w = 1,
                           cow_align = 'h',
                           show_plot = NA,
                           return_plot = NA,
                           save_plot = NA,
                           save_param =  list(),
                           default_save_name = 'spatGenePlot2D') {


  # print, return and save parameters
  show_plot = ifelse(is.na(show_plot), readGiottoInstructions(gobject, param = 'show_plot'), show_plot)
  save_plot = ifelse(is.na(save_plot), readGiottoInstructions(gobject, param = 'save_plot'), save_plot)
  return_plot = ifelse(is.na(return_plot), readGiottoInstructions(gobject, param = 'return_plot'), return_plot)

  # point shape
  point_shape = match.arg(point_shape, choices = c('border', 'no_border'))

  # expression values
  values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))
  expr_values = Giotto:::select_expression_values(gobject = gobject, values = values)

  # only keep genes that are in the dataset
  selected_genes = genes
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
    spatial_network = select_spatialNetwork(gobject,name = spatial_network_name,return_network_Obj = FALSE)
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
  cell_metadata = combineMetadata(gobject = gobject)
  #cell_metadata   = gobject@cell_metadata
  #cell_metadata   = cell_metadata[, !grepl('cell_ID', colnames(cell_metadata)), with = F]

  if(nrow(cell_metadata) == 0) {
    cell_locations_metadata = cell_locations
  } else {
    cell_locations_metadata = cell_metadata
    #cell_locations_metadata <- cbind(cell_locations, cell_metadata)
  }

  cell_locations_metadata_genes <- merge(cell_locations_metadata, t_sub_expr_data_DT, by = 'cell_ID')

  ## plotting ##

  savelist <- list()

  for(gene in selected_genes) {

    pl <- ggplot2::ggplot()
    pl <- pl + ggplot2::theme_classic()

    ## plot spatial network
    if(!is.null(spatial_network) & show_network == TRUE) {
      if(is.null(network_color)) {
        network_color = 'red'
      }
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


    ### plot cells ###

    ## with border ##
    if(point_shape == 'border') {

      if(scale_alpha_with_expression == TRUE) {
        pl <- pl + ggplot2::geom_point(data = cell_locations_metadata_genes, aes_string2(x = 'sdimx', y = 'sdimy',
                                                                                         fill = gene, alpha = gene),
                                       shape = 21,
                                       color = point_border_col, size = point_size, stroke = point_border_stroke,
                                       show.legend = show_legend)
      } else {
        pl <- pl + ggplot2::geom_point(data = cell_locations_metadata_genes, aes_string2(x = 'sdimx', y = 'sdimy',
                                                                                         fill = gene),
                                       shape = 21,
                                       color = point_border_col, size = point_size, stroke = point_border_stroke,
                                       show.legend = show_legend)
      }
      pl <- pl + ggplot2::scale_alpha_continuous(guide = 'none')
      pl <- pl + ggplot2::scale_fill_gradient2(low = genes_low_color, mid = genes_mid_color, high = genes_high_color,
                                               midpoint = midpoint, guide = guide_colorbar(title = ''))
      pl <- pl + ggplot2::labs(x = 'coord x', y = 'coord y', title = gene)

    }


    ## no border ##
    if(point_shape == 'no_border') {

      if(scale_alpha_with_expression == TRUE) {
        pl <- pl + ggplot2::geom_point(data = cell_locations_metadata_genes, aes_string2(x = 'sdimx', y = 'sdimy',
                                                                                         color = gene, alpha = gene),
                                       shape = 19, size = point_size,  show.legend = show_legend)
      } else {
        pl <- pl + ggplot2::geom_point(data = cell_locations_metadata_genes, aes_string2(x = 'sdimx', y = 'sdimy',
                                                                                         color = gene),
                                       shape = 19, size = point_size, show.legend = show_legend)
      }
      pl <- pl + ggplot2::scale_alpha_continuous(guide = 'none')
      pl <- pl + ggplot2::scale_color_gradient2(low = genes_low_color, mid = genes_mid_color, high = genes_high_color,
                                               midpoint = midpoint, guide = guide_colorbar(title = ''))
      pl <- pl + ggplot2::labs(x = 'coord x', y = 'coord y', title = gene)

    }




    pl <- pl + ggplot2::theme(plot.title = element_text(hjust = 0.5),
                              legend.title = element_blank(),
                              legend.text = element_text(size = legend_text),
                              axis.title = element_text(size = axis_title),
                              axis.text = element_text(size = axis_text),
                              panel.grid = element_blank(),
                              panel.background = element_rect(fill = background_color))


    savelist[[gene]] <- pl
  }

  # combine plots with cowplot
  combo_plot <- cowplot::plot_grid(plotlist = savelist,
                                   ncol = cow_n_col,
                                   rel_heights = cow_rel_h, rel_widths = cow_rel_w, align = cow_align)


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





#' @title spatGenePlot
#' @name spatGenePlot
#' @description Visualize cells and gene expression according to spatial coordinates
#' @param gobject giotto object
#' @param expression_values gene expression values to use
#' @param genes genes to show
#' @param genes_high_color color represents high gene expression
#' @param genes_mid_color color represents middle gene expression
#' @param genes_low_color color represents low gene expression
#' @param show_network show underlying spatial network
#' @param network_color color of spatial network
#' @param spatial_network_name name of spatial network to use
#' @param show_grid show spatial grid
#' @param grid_color color of spatial grid
#' @param spatial_grid_name name of spatial grid to use
#' @param midpoint expression midpoint
#' @param scale_alpha_with_expression scale expression with ggplot alpha parameter
#' @param point_shape point with border or not (border or no_border)
#' @param point_size size of point (cell)
#' @param point_border_col color of border around points
#' @param point_border_stroke stroke size of border around points
#' @param cow_n_col cowplot param: how many columns
#' @param cow_rel_h cowplot param: relative height
#' @param cow_rel_w cowplot param: relative width
#' @param cow_align cowplot param: how to align
#' @param show_legend show legend
#' @param legend_text size of legend text
#' @param background_color color of plot background
#' @param axis_text size of axis text
#' @param axis_title size of axis title
#' @param show_plot show plots
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters from \code{\link{all_plots_save_function}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @param ... parameters for cowplot::save_plot()
#' @return ggplot
#' @details Description of parameters.
#' @export
#' @seealso \code{\link{spatGenePlot3D}} and \code{\link{spatGenePlot2D}}
#' @examples
#'     spatGenePlot(gobject)
spatGenePlot = function(gobject,
                        expression_values = c('normalized', 'scaled', 'custom'),
                        genes,
                        genes_high_color = 'darkred',
                        genes_mid_color = "white",
                        genes_low_color = "darkblue",
                        show_network = F,
                        network_color = NULL,
                        spatial_network_name = 'Delaunay_network',
                        edge_alpha = NULL,
                        show_grid = F,
                        grid_color = NULL,
                        spatial_grid_name = 'spatial_grid',
                        midpoint = 0,
                        scale_alpha_with_expression = FALSE,
                        point_shape = c('border', 'no_border'),
                        point_size = 1,
                        point_border_col = 'black',
                        point_border_stroke = 0.1,
                        show_legend = T,
                        legend_text = 8,
                        background_color = 'white',
                        axis_text = 8,
                        axis_title = 8,
                        cow_n_col = 2,
                        cow_rel_h = 1,
                        cow_rel_w = 1,
                        cow_align = 'h',
                        show_plot = NA,
                        return_plot = NA,
                        save_plot = NA,
                        save_param =  list(),
                        default_save_name = 'spatGenePlot') {

  spatGenePlot2D(gobject = gobject,
                 expression_values = expression_values,
                 genes = genes,
                 genes_high_color = genes_high_color,
                 genes_mid_color = genes_mid_color,
                 genes_low_color = genes_low_color,
                 show_network = show_network,
                 network_color = network_color,
                 spatial_network_name = spatial_network_name,
                 edge_alpha = edge_alpha,
                 show_grid = show_grid,
                 grid_color = grid_color,
                 spatial_grid_name = spatial_grid_name,
                 midpoint = midpoint,
                 scale_alpha_with_expression = scale_alpha_with_expression,
                 point_shape = point_shape,
                 point_size = point_size,
                 point_border_col = point_border_col,
                 point_border_stroke = point_border_stroke,
                 show_legend = show_legend,
                 legend_text = legend_text,
                 background_color = background_color,
                 axis_text = axis_text,
                 axis_title = axis_title,
                 cow_n_col = cow_n_col,
                 cow_rel_h = cow_rel_h,
                 cow_rel_w = cow_rel_w,
                 cow_align = cow_align,
                 show_plot = show_plot,
                 return_plot = return_plot,
                 save_plot = save_plot,
                 save_param =  save_param,
                 default_save_name = default_save_name)

}



#' @title dimGenePlot2D
#' @name dimGenePlot2D
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
#' @param point_shape point with border or not (border or no_border)
#' @param point_size size of point (cell)
#' @param point_border_col color of border around points
#' @param point_border_stroke stroke size of border around points
#' @param midpoint size of point (cell)
#' @param show_legend show legend
#' @param legend_text size of legend text
#' @param background_color color of plot background
#' @param axis_text size of axis text
#' @param axis_title size of axis title
#' @param cow_n_col cowplot param: how many columns
#' @param cow_rel_h cowplot param: relative height
#' @param cow_rel_w cowplot param: relative width
#' @param cow_align cowplot param: how to align
#' @param show_plot show plots
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters from \code{\link{all_plots_save_function}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @param ... parameters for cowplot::save_plot()
#' @return ggplot
#' @details Description of parameters.
#' @export
#' @seealso \code{\link{dimGenePlot3D}}
#' @examples
#'     dimGenePlot2D(gobject)
dimGenePlot2D <- function(gobject,
                          expression_values = c('normalized', 'scaled', 'custom'),
                          genes = NULL,
                          dim_reduction_to_use = 'umap',
                          dim_reduction_name = 'umap',
                          dim1_to_use = 1,
                          dim2_to_use = 2,
                          show_NN_network = F,
                          nn_network_to_use = 'sNN',
                          network_name = 'sNN.pca',
                          network_color = "lightgray",
                          edge_alpha = NULL,
                          scale_alpha_with_expression = FALSE,
                          point_shape = c('border', 'no_border'),
                          point_size = 1,
                          genes_high_color = "red",
                          genes_mid_color = "white",
                          genes_low_color = "blue",
                          point_border_col = 'black',
                          point_border_stroke = 0.1,
                          midpoint = 0,
                          show_legend = T,
                          legend_text = 8,
                          background_color = 'white',
                          axis_text = 8,
                          axis_title = 8,
                          cow_n_col = 2,
                          cow_rel_h = 1,
                          cow_rel_w = 1,
                          cow_align = 'h',
                          show_plot = NA,
                          return_plot = NA,
                          save_plot = NA,
                          save_param =  list(),
                          default_save_name = 'dimGenePlot2D') {


  # print, return and save parameters
  show_plot = ifelse(is.na(show_plot), readGiottoInstructions(gobject, param = 'show_plot'), show_plot)
  save_plot = ifelse(is.na(save_plot), readGiottoInstructions(gobject, param = 'save_plot'), save_plot)
  return_plot = ifelse(is.na(return_plot), readGiottoInstructions(gobject, param = 'return_plot'), return_plot)

  # point shape
  point_shape = match.arg(point_shape, choices = c('border', 'no_border'))

  ## select genes ##
  selected_genes = genes
  values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))
  expr_values = Giotto:::select_expression_values(gobject = gobject, values = values)

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



  ## visualize multipe plots ##
  ## 2D plots ##
  savelist <- list()

  for(gene in selected_genes) {


    ## OLD need to be combined ##
    pl <- ggplot2::ggplot()
    pl <- pl + ggplot2::theme_classic()

    # network layer
    if(show_NN_network == TRUE) {

      if(is.null(edge_alpha)) {
        edge_alpha = 0.5
        pl <- pl + ggplot2::geom_segment(data = annotated_network_DT,
                                         aes_string(x = from_dim_names[1], y = from_dim_names[2],
                                                    xend = to_dim_names[1], yend = to_dim_names[2]),
                                         alpha = edge_alpha, color=network_color,size = 0.1,
                                         show.legend = F)
      } else if(is.numeric(edge_alpha)) {
        pl <- pl + ggplot2::geom_segment(data = annotated_network_DT,
                                         aes_string(x = from_dim_names[1], y = from_dim_names[2],
                                                    xend = to_dim_names[1], yend = to_dim_names[2]),
                                         alpha = edge_alpha, color=network_color,size = 0.1,
                                         show.legend = F)
      } else if(is.character(edge_alpha)) {

        if(edge_alpha %in% colnames(annotated_network_DT)) {
          pl <- pl + ggplot2::geom_segment(data = annotated_network_DT,
                                           aes_string(x = from_dim_names[1], y = from_dim_names[2],
                                                      xend = to_dim_names[1],
                                                      yend = to_dim_names[2], alpha = edge_alpha),
                                           color=network_color,
                                           show.legend = F)
        }
      }
    }



    # point layer
    if(is.null(genes)) {
      cell_color = 'lightblue'
      pl <- pl + ggplot2::geom_point(data = annotated_gene_DT, aes_string(x = dim_names[1], dim_names[2]),
                                     fill = cell_color, show.legend = show_legend, size =  point_size)

    } else {

      # with border
      if(point_shape == 'border') {

        if(scale_alpha_with_expression == TRUE) {
          pl <- pl + ggplot2::geom_point(data = annotated_gene_DT, aes_string2(x = dim_names[1], y = dim_names[2], fill = gene, alpha = gene),
                                         show.legend = show_legend, shape = 21, size = point_size,
                                         color = point_border_col, stroke = point_border_stroke)
        } else {
          pl <- pl + ggplot2::geom_point(data = annotated_gene_DT, aes_string2(x = dim_names[1], y = dim_names[2], fill = gene),
                                         show.legend = show_legend, shape = 21,
                                         size =  point_size,
                                         color = point_border_col, stroke = point_border_stroke)
        }

        pl <- pl + ggplot2::scale_fill_gradient2(low = genes_low_color, mid = genes_mid_color, high = genes_high_color, midpoint = midpoint)

      }


      # without border
      # with border
      if(point_shape == 'no_border') {

        if(scale_alpha_with_expression == TRUE) {
          pl <- pl + ggplot2::geom_point(data = annotated_gene_DT, aes_string2(x = dim_names[1], y = dim_names[2], color = gene, alpha = gene),
                                         show.legend = show_legend, shape = 19, size = point_size)
        } else {
          pl <- pl + ggplot2::geom_point(data = annotated_gene_DT, aes_string2(x = dim_names[1], y = dim_names[2], color = gene),
                                         show.legend = show_legend, shape = 19, size =  point_size)
        }

        pl <- pl + ggplot2::scale_color_gradient2(low = genes_low_color, mid = genes_mid_color, high = genes_high_color, midpoint = midpoint)

      }


    }


    ## add title
    pl <- pl + ggplot2::labs(x = 'coord x', y = 'coord y', title = gene)

    ## aesthetics
    pl <- pl + ggplot2::theme(plot.title = element_text(hjust = 0.5),
                              legend.title = element_blank(),
                              legend.text = element_text(size = legend_text),
                              axis.title = element_text(size = axis_title),
                              axis.text = element_text(size = axis_text),
                              panel.grid = element_blank(),
                              panel.background = element_rect(fill = background_color))

    savelist[[gene]] <- pl
  }




  # combine plots with cowplot
  combo_plot <- cowplot::plot_grid(plotlist = savelist,
                                   ncol = cow_n_col,
                                   rel_heights = cow_rel_h, rel_widths = cow_rel_w, align = cow_align)


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





#' @title dimGenePlot
#' @name dimGenePlot
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
#' @param point_border_col color of border around points
#' @param point_border_stroke stroke size of border around points
#' @param midpoint size of point (cell)
#' @param cow_n_col cowplot param: how many columns
#' @param cow_rel_h cowplot param: relative height
#' @param cow_rel_w cowplot param: relative width
#' @param cow_align cowplot param: how to align
#' @param show_legend show legend
#' @param show_plot show plots
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters from \code{\link{all_plots_save_function}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @param ... parameters for cowplot::save_plot()
#' @return ggplot
#' @details Description of parameters.
#' @export
#' @seealso \code{\link{dimGenePlot3D}}
#' @examples
#'     dimGenePlot(gobject)
dimGenePlot = function(gobject,
                       expression_values = c('normalized', 'scaled', 'custom'),
                       genes = NULL,
                       dim_reduction_to_use = 'umap',
                       dim_reduction_name = 'umap',
                       dim1_to_use = 1,
                       dim2_to_use = 2,
                       show_NN_network = F,
                       nn_network_to_use = 'sNN',
                       network_name = 'sNN.pca',
                       network_color = "lightgray",
                       edge_alpha = NULL,
                       scale_alpha_with_expression = FALSE,
                       point_shape = c('border', 'no_border'),
                       point_size = 1,
                       genes_high_color = "red",
                       genes_mid_color = "white",
                       genes_low_color = "blue",
                       point_border_col = 'black',
                       point_border_stroke = 0.1,
                       midpoint = 0,
                       show_legend = T,
                       legend_text = 8,
                       background_color = 'white',
                       axis_text = 8,
                       axis_title = 8,
                       cow_n_col = 2,
                       cow_rel_h = 1,
                       cow_rel_w = 1,
                       cow_align = 'h',
                       show_plot = NA,
                       return_plot = NA,
                       save_plot = NA,
                       save_param =  list(),
                       default_save_name = 'dimGenePlot') {

  dimGenePlot2D(gobject = gobject,
                expression_values = expression_values,
                genes = genes,
                dim_reduction_to_use = dim_reduction_to_use,
                dim_reduction_name = dim_reduction_name,
                dim1_to_use = dim1_to_use,
                dim2_to_use = dim2_to_use,
                show_NN_network = show_NN_network,
                nn_network_to_use = nn_network_to_use,
                network_name = network_name,
                network_color = network_color,
                edge_alpha = edge_alpha,
                scale_alpha_with_expression = scale_alpha_with_expression,
                point_shape = point_shape,
                point_size = point_size,
                genes_high_color = genes_high_color,
                genes_mid_color = genes_mid_color,
                genes_low_color = genes_low_color,
                point_border_col = point_border_col,
                point_border_stroke = point_border_stroke,
                midpoint = midpoint,
                show_legend = show_legend,
                legend_text = legend_text,
                background_color = background_color,
                axis_text = axis_text,
                axis_title = axis_title,
                cow_n_col = cow_n_col,
                cow_rel_h = cow_rel_h,
                cow_rel_w = cow_rel_w,
                cow_align = cow_align,
                show_plot = show_plot,
                return_plot = return_plot,
                save_plot = save_plot,
                save_param = save_param,
                default_save_name = default_save_name)

}






#' @title spatDimGenePlot2D
#' @name spatDimGenePlot2D
#' @description Visualize cells according to spatial AND dimension reduction coordinates in ggplot mode
#' @param gobject giotto object
#' @param expression_values gene expression values to use
#' @param plot_alignment direction to align plot
#' @param genes genes to show
#' @param dim_reduction_to_use dimension reduction to use
#' @param dim_reduction_name dimension reduction name
#' @param dim1_to_use dimension to use on x-axis
#' @param dim2_to_use dimension to use on y-axis
#' @param dim_point_shape dim reduction points with border or not (border or no_border)
#' @param dim_point_size dim reduction plot: point size
#' @param dim_point_border_col color of border around points
#' @param dim_point_border_stroke stroke size of border around points
#' @param show_NN_network show underlying NN network
#' @param nn_network_to_use type of NN network to use (kNN vs sNN)
#' @param network_name name of NN network to use, if show_NN_network = TRUE
#' @param edge_alpha_dim dim reduction plot: column to use for alpha of the edges
#' @param scale_alpha_with_expression scale expression with ggplot alpha parameter
#' @param spatial_network_name name of spatial network to use
#' @param spatial_grid_name name of spatial grid to use
#' @param spat_point_shape spatial points with border or not (border or no_border)
#' @param spat_point_size spatial plot: point size
#' @param spat_point_border_col color of border around points
#' @param spat_point_border_stroke stroke size of border around points
#' @param midpoint size of point (cell)
#' @param show_legend show legend
#' @param legend_text size of legend text
#' @param dim_background_color color of plot background for dimension plot
#' @param spat_background_color color of plot background for spatial plot
#' @param axis_text size of axis text
#' @param axis_title size of axis title
#' @param cow_n_col cowplot param: how many columns
#' @param cow_rel_h cowplot param: relative height
#' @param cow_rel_w cowplot param: relative width
#' @param cow_align cowplot param: how to align
#' @param show_plot show plots
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters from \code{\link{all_plots_save_function}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return ggplot
#' @details Description of parameters.
#' @export
#' @seealso \code{\link{spatDimGenePlot3D}}
#' @examples
#'     spatDimGenePlot2D(gobject)
spatDimGenePlot2D <- function(gobject,
                              expression_values = c('normalized', 'scaled', 'custom'),
                              plot_alignment = c('vertical', 'horizontal'),
                              genes,
                              dim_reduction_to_use = 'umap',
                              dim_reduction_name = 'umap',
                              dim1_to_use = 1,
                              dim2_to_use = 2,
                              dim_point_shape = c('border', 'no_border'),
                              dim_point_size = 1,
                              dim_point_border_col = 'black',
                              dim_point_border_stroke = 0.1,
                              show_NN_network = F,
                              show_spatial_network = F,
                              show_spatial_grid = F,
                              nn_network_to_use = 'sNN',
                              network_name = 'sNN.pca',
                              edge_alpha_dim = NULL,
                              scale_alpha_with_expression = FALSE,
                              spatial_network_name = 'Delaunay_network',
                              spatial_grid_name = 'spatial_grid',
                              spat_point_shape = c('border', 'no_border'),
                              spat_point_size = 1,
                              spat_point_border_col = 'black',
                              spat_point_border_stroke = 0.1,
                              midpoint = 0,
                              genes_high_color = "red",
                              genes_mid_color= "white",
                              genes_low_color = "blue",
                              cow_n_col = 2,
                              cow_rel_h = 1,
                              cow_rel_w = 1,
                              cow_align = 'h',
                              show_legend = T,
                              legend_text = 8,
                              dim_background_color = 'white',
                              spat_background_color = 'white',
                              axis_text = 8,
                              axis_title = 8,
                              show_plot = NA,
                              return_plot = NA,
                              save_plot = NA,
                              save_param =  list(),
                              default_save_name = 'spatDimGenePlot2D') {

  plot_alignment = match.arg(plot_alignment, choices = c('vertical', 'horizontal'))

  # dimension reduction plot
  dmpl = dimGenePlot2D(gobject = gobject,
                       expression_values = expression_values,
                       genes = genes,
                       dim_reduction_to_use = dim_reduction_to_use,
                       dim_reduction_name = dim_reduction_name,
                       dim1_to_use = dim1_to_use,
                       dim2_to_use = dim2_to_use,
                       show_NN_network = show_NN_network,
                       nn_network_to_use = nn_network_to_use,
                       network_name =network_name,
                       network_color = network_color,
                       edge_alpha = edge_alpha,
                       scale_alpha_with_expression = scale_alpha_with_expression,
                       point_shape = dim_point_shape,
                       point_size = dim_point_size,
                       genes_high_color = genes_high_color,
                       genes_mid_color= genes_mid_color,
                       genes_low_color = genes_low_color,
                       point_border_col = dim_point_border_col,
                       point_border_stroke = dim_point_border_stroke,
                       midpoint = midpoint,
                       show_legend = show_legend,
                       legend_text = legend_text,
                       background_color = dim_background_color,
                       axis_text = axis_text,
                       axis_title = axis_title,
                       cow_n_col = cow_n_col,
                       cow_rel_h = cow_rel_h,
                       cow_rel_w = cow_rel_w,
                       cow_align = cow_align,
                       show_plot = FALSE,
                       return_plot = TRUE,
                       save_plot = FALSE)

  # spatial plot
  spl = spatGenePlot2D(gobject=gobject,
                       expression_values = expression_values,
                       genes = genes,
                       genes_high_color = genes_high_color,
                       genes_mid_color = genes_mid_color,
                       genes_low_color = genes_low_color,
                       show_network = show_spatial_network,
                       network_color = network_color,
                       spatial_network_name = spatial_network_name,
                       edge_alpha = edge_alpha,
                       show_grid = show_spatial_grid,
                       grid_color = grid_color,
                       spatial_grid_name = spatial_grid_name,
                       midpoint = midpoint,
                       scale_alpha_with_expression = scale_alpha_with_expression,
                       point_shape = spat_point_shape,
                       point_size = spat_point_size,
                       point_border_col =spat_point_border_col,
                       point_border_stroke = spat_point_border_stroke,
                       show_legend = show_legend,
                       legend_text = legend_text,
                       background_color = spat_background_color,
                       axis_text = axis_text,
                       axis_title = axis_title,
                       cow_n_col = cow_n_col,
                       cow_rel_h = cow_rel_h,
                       cow_rel_w = cow_rel_w,
                       cow_align = cow_align,
                       show_plot = FALSE,
                       return_plot = TRUE,
                       save_plot = FALSE)


  if(plot_alignment == 'vertical') {
    ncol = 1
    nrow = 2
    combo_plot = cowplot::plot_grid(dmpl, spl, ncol = ncol, nrow = nrow, rel_heights = c(1), rel_widths = c(1), align = 'v')
  } else {
    ncol = 2
    nrow = 1
    combo_plot = cowplot::plot_grid(dmpl, spl, ncol = ncol, nrow = nrow, rel_heights = c(1), rel_widths = c(1), align = 'h')
  }

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







#' @title spatDimGenePlot
#' @name spatDimGenePlot
#' @description Visualize cells according to spatial AND dimension reduction coordinates in ggplot mode
#' @param gobject giotto object
#' @param expression_values gene expression values to use
#' @param plot_alignment direction to align plot
#' @param genes genes to show
#' @param dim_reduction_to_use dimension reduction to use
#' @param dim_reduction_name dimension reduction name
#' @param dim1_to_use dimension to use on x-axis
#' @param dim2_to_use dimension to use on y-axis
#' @param dim_point_shape dimension points with border or not (border or no_border)
#' @param dim_point_size dim reduction plot: point size
#' @param dim_point_border_col color of border around points
#' @param dim_point_border_stroke stroke size of border around points
#' @param show_NN_network show underlying NN network
#' @param nn_network_to_use type of NN network to use (kNN vs sNN)
#' @param network_name name of NN network to use, if show_NN_network = TRUE
#' @param edge_alpha_dim dim reduction plot: column to use for alpha of the edges
#' @param scale_alpha_with_expression scale expression with ggplot alpha parameter
#' @param spatial_network_name name of spatial network to use
#' @param spatial_grid_name name of spatial grid to use
#' @param spat_point_shape spatial points with border or not (border or no_border)
#' @param spat_point_size spatial plot: point size
#' @param spat_point_border_col color of border around points
#' @param spat_point_border_stroke stroke size of border around points
#' @param midpoint size of point (cell)
#' @param show_legend show legend
#' @param legend_text size of legend text
#' @param dim_background_color color of plot background for dimension plot
#' @param spat_background_color color of plot background for spatial plot
#' @param axis_text size of axis text
#' @param axis_title size of axis title
#' @param cow_n_col cowplot param: how many columns
#' @param cow_rel_h cowplot param: relative height
#' @param cow_rel_w cowplot param: relative width
#' @param cow_align cowplot param: how to align
#' @param show_plot show plots
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters from \code{\link{all_plots_save_function}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return ggplot
#' @details Description of parameters.
#' @export
#' @seealso \code{\link{spatDimGenePlot3D}}
#' @examples
#'     spatDimGenePlot(gobject)
spatDimGenePlot = function(gobject,
                           expression_values = c('normalized', 'scaled', 'custom'),
                           plot_alignment = c('vertical', 'horizontal'),
                           genes,
                           dim_reduction_to_use = 'umap',
                           dim_reduction_name = 'umap',
                           dim1_to_use = 1,
                           dim2_to_use = 2,
                           dim_point_shape = c('border', 'no_border'),
                           dim_point_size = 1,
                           dim_point_border_col = 'black',
                           dim_point_border_stroke = 0.1,
                           show_NN_network = F,
                           show_spatial_network = F,
                           show_spatial_grid = F,
                           nn_network_to_use = 'sNN',
                           network_name = 'sNN.pca',
                           edge_alpha_dim = NULL,
                           scale_alpha_with_expression = FALSE,
                           spatial_network_name = 'Delaunay_network',
                           spatial_grid_name = 'spatial_grid',
                           spat_point_shape = c('border', 'no_border'),
                           spat_point_size = 1,
                           spat_point_border_col = 'black',
                           spat_point_border_stroke = 0.1,
                           midpoint = 0,
                           genes_high_color = "red",
                           genes_mid_color="white",
                           genes_low_color = "blue",
                           show_legend = T,
                           legend_text = 8,
                           dim_background_color = 'white',
                           spat_background_color = 'white',
                           axis_text = 8,
                           axis_title = 8,
                           cow_n_col = 2,
                           cow_rel_h = 1,
                           cow_rel_w = 1,
                           cow_align = 'h',
                           show_plot = NA,
                           return_plot = NA,
                           save_plot = NA,
                           save_param =  list(),
                           default_save_name = 'spatDimGenePlot') {

  spatDimGenePlot2D(gobject = gobject,
                    expression_values = expression_values,
                    plot_alignment = plot_alignment,
                    genes = genes,
                    dim_reduction_to_use = dim_reduction_to_use,
                    dim_reduction_name = dim_reduction_name,
                    dim1_to_use = dim1_to_use,
                    dim2_to_use = dim2_to_use,
                    dim_point_shape = dim_point_shape,
                    dim_point_size = dim_point_size,
                    dim_point_border_col = dim_point_border_col,
                    dim_point_border_stroke = dim_point_border_stroke,
                    show_NN_network = show_NN_network,
                    show_spatial_network = show_spatial_network,
                    show_spatial_grid = show_spatial_grid,
                    nn_network_to_use = nn_network_to_use,
                    network_name = network_name,
                    edge_alpha_dim = edge_alpha_dim,
                    scale_alpha_with_expression = scale_alpha_with_expression,
                    spatial_network_name = spatial_network_name,
                    spatial_grid_name = spatial_grid_name,
                    spat_point_shape = spat_point_shape,
                    spat_point_size = spat_point_size,
                    spat_point_border_col = spat_point_border_col,
                    spat_point_border_stroke = spat_point_border_stroke,
                    midpoint = midpoint,
                    genes_high_color = genes_high_color,
                    genes_mid_color = genes_mid_color,
                    genes_low_color = genes_low_color,
                    show_legend = show_legend,
                    legend_text = legend_text,
                    dim_background_color = dim_background_color,
                    spat_background_color = spat_background_color,
                    axis_text = axis_text,
                    axis_title = axis_title,
                    cow_n_col = cow_n_col,
                    cow_rel_h = cow_rel_h,
                    cow_rel_w = cow_rel_w,
                    cow_align = cow_align,
                    show_plot = show_plot,
                    return_plot = return_plot,
                    save_plot = save_plot,
                    save_param =  save_param,
                    default_save_name = default_save_name)

}







#' @title spatCellPlot2D
#' @name spatCellPlot2D
#' @description Visualize cells according to spatial coordinates
#' @param gobject giotto object
#' @param sdimx x-axis dimension name (default = 'sdimx')
#' @param sdimy y-axis dimension name (default = 'sdimy')
#' @param spat_enr_names names of spatial enrichment results to include
#' @param cell_annotation_values numeric cell annotation columns
#' @param cell_color_gradient vector with 3 colors for numeric data
#' @param gradient_midpoint midpoint for color gradient
#' @param gradient_limits vector with lower and upper limits
#' @param select_cell_groups select subset of cells/clusters based on cell_color parameter
#' @param select_cells select subset of cells based on cell IDs
#' @param point_shape point with border or not (border or no_border)
#' @param point_size size of point (cell)
#' @param point_border_col color of border around points
#' @param point_border_stroke stroke size of border around points
#' @param show_cluster_center plot center of selected clusters
#' @param show_center_label plot label of selected clusters
#' @param center_point_size size of center points
#' @param label_size  size of labels
#' @param label_fontface font of labels
#' @param show_network show underlying spatial network
#' @param spatial_network_name name of spatial network to use
#' @param network_color color of spatial network
#' @param network_alpha alpha of spatial network
#' @param show_grid show spatial grid
#' @param spatial_grid_name name of spatial grid to use
#' @param grid_color color of spatial grid
#' @param show_other_cells display not selected cells
#' @param other_cell_color color of not selected cells
#' @param other_point_size point size of not selected cells
#' @param other_cells_alpha alpha of not selected cells
#' @param coord_fix_ratio fix ratio between x and y-axis
#' @param show_legend show legend
#' @param legend_text size of legend text
#' @param legend_symbol_size size of legend symbols
#' @param background_color color of plot background
#' @param axis_text size of axis text
#' @param axis_title size of axis title
#' @param show_plot show plot
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters from \code{\link{all_plots_save_function}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return ggplot
#' @details Description of parameters.
#' @export
#' @examples
#'     spatCellPlot2D(gobject)
spatCellPlot2D = function(gobject,
                          sdimx = 'sdimx',
                          sdimy = 'sdimy',
                          spat_enr_names = NULL,
                          cell_annotation_values = NULL,
                          cell_color_gradient = c('blue', 'white', 'red'),
                          gradient_midpoint = NULL,
                          gradient_limits = NULL,
                          select_cell_groups = NULL,
                          select_cells = NULL,
                          point_shape = c('border', 'no_border'),
                          point_size = 3,
                          point_border_col = 'black',
                          point_border_stroke = 0.1,
                          show_cluster_center = F,
                          show_center_label = F,
                          center_point_size = 4,
                          center_point_border_col = 'black',
                          center_point_border_stroke = 0.1,
                          label_size = 4,
                          label_fontface = 'bold',
                          show_network = F,
                          spatial_network_name = 'Delaunay_network',
                          network_color = NULL,
                          network_alpha = 1,
                          show_grid = F,
                          spatial_grid_name = 'spatial_grid',
                          grid_color = NULL,
                          show_other_cells = T,
                          other_cell_color = 'lightgrey',
                          other_point_size = 1,
                          other_cells_alpha = 0.1,
                          coord_fix_ratio = NULL,
                          show_legend = T,
                          legend_text = 8,
                          legend_symbol_size = 1,
                          background_color = 'white',
                          axis_text = 8,
                          axis_title = 8,
                          cow_n_col = 2,
                          cow_rel_h = 1,
                          cow_rel_w = 1,
                          cow_align = 'h',
                          show_plot = NA,
                          return_plot = NA,
                          save_plot = NA,
                          save_param =  list(),
                          default_save_name = 'spatCellPlot2D'
) {


  comb_metadata = combineMetadata(gobject = gobject,
                                  spat_enr_names = spat_enr_names)

  # keep only available columns
  possible_value_cols = colnames(comb_metadata)
  if(is.null(cell_annotation_values)) {
    stop('you need to choose which continuous/numerical cell annotations or enrichments you want to visualize')
  }
  cell_annotation_values = cell_annotation_values[cell_annotation_values %in% possible_value_cols]


  # print, return and save parameters
  show_plot = ifelse(is.na(show_plot), readGiottoInstructions(gobject, param = 'show_plot'), show_plot)
  save_plot = ifelse(is.na(save_plot), readGiottoInstructions(gobject, param = 'save_plot'), save_plot)
  return_plot = ifelse(is.na(return_plot), readGiottoInstructions(gobject, param = 'return_plot'), return_plot)

  ## plotting ##
  savelist <- list()

  for(annot in cell_annotation_values) {

    pl = spatPlot2D(gobject = gobject,
                    group_by = NULL,
                    group_by_subset = NULL,
                    sdimx = sdimx,
                    sdimy = sdimy,
                    spat_enr_names = spat_enr_names,
                    cell_color = annot,
                    color_as_factor = F,
                    cell_color_gradient = cell_color_gradient,
                    gradient_midpoint = gradient_midpoint,
                    gradient_limits = gradient_limits,
                    select_cell_groups = select_cell_groups,
                    select_cells = select_cells,
                    point_shape = point_shape,
                    point_size = point_size,
                    point_border_col = point_border_col,
                    point_border_stroke = point_border_stroke,
                    show_cluster_center = show_cluster_center,
                    show_center_label = show_center_label,
                    center_point_size = center_point_size,
                    center_point_border_col = center_point_border_col,
                    center_point_border_stroke = center_point_border_stroke,
                    label_size = label_size,
                    label_fontface = label_fontface,
                    show_network = show_network,
                    spatial_network_name = spatial_network_name,
                    network_color = network_color,
                    network_alpha = network_alpha,
                    show_grid = show_grid,
                    spatial_grid_name = spatial_grid_name,
                    grid_color = grid_color,
                    show_other_cells = show_other_cells,
                    other_cell_color = other_cell_color,
                    other_point_size = other_point_size,
                    other_cells_alpha = other_cells_alpha,
                    coord_fix_ratio = coord_fix_ratio,
                    title = annot,
                    show_legend = show_legend,
                    legend_text = legend_text,
                    legend_symbol_size = legend_symbol_size,
                    background_color = background_color,
                    axis_text = axis_text,
                    axis_title = axis_title,
                    show_plot = FALSE,
                    return_plot = TRUE,
                    save_plot = FALSE,
                    save_param =  list(),
                    default_save_name = 'spatPlot2D')


    savelist[[annot]] <- pl

  }


  # combine plots with cowplot
  combo_plot <- cowplot::plot_grid(plotlist = savelist,
                                   ncol = cow_n_col,
                                   rel_heights = cow_rel_h,
                                   rel_widths = cow_rel_w,
                                   align = cow_align)


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


#' @title spatCellPlot
#' @name spatCellPlot
#' @description Visualize cells according to spatial coordinates
#' @param gobject giotto object
#' @param sdimx x-axis dimension name (default = 'sdimx')
#' @param sdimy y-axis dimension name (default = 'sdimy')
#' @param spat_enr_names names of spatial enrichment results to include
#' @param cell_annotation_values numeric cell annotation columns
#' @param cell_color_gradient vector with 3 colors for numeric data
#' @param gradient_midpoint midpoint for color gradient
#' @param gradient_limits vector with lower and upper limits
#' @param select_cell_groups select subset of cells/clusters based on cell_color parameter
#' @param select_cells select subset of cells based on cell IDs
#' @param point_shape point with border or not (border or no_border)
#' @param point_size size of point (cell)
#' @param point_border_col color of border around points
#' @param point_border_stroke stroke size of border around points
#' @param show_cluster_center plot center of selected clusters
#' @param show_center_label plot label of selected clusters
#' @param center_point_size size of center points
#' @param label_size  size of labels
#' @param label_fontface font of labels
#' @param show_network show underlying spatial network
#' @param spatial_network_name name of spatial network to use
#' @param network_color color of spatial network
#' @param network_alpha alpha of spatial network
#' @param show_grid show spatial grid
#' @param spatial_grid_name name of spatial grid to use
#' @param grid_color color of spatial grid
#' @param show_other_cells display not selected cells
#' @param other_cell_color color of not selected cells
#' @param other_point_size point size of not selected cells
#' @param other_cells_alpha alpha of not selected cells
#' @param coord_fix_ratio fix ratio between x and y-axis
#' @param show_legend show legend
#' @param legend_text size of legend text
#' @param legend_symbol_size size of legend symbols
#' @param background_color color of plot background
#' @param axis_text size of axis text
#' @param axis_title size of axis title
#' @param show_plot show plot
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters from \code{\link{all_plots_save_function}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return ggplot
#' @details Description of parameters.
#' @export
#' @examples
#'     spatCellPlot(gobject)
spatCellPlot = function(gobject,
                        sdimx = 'sdimx',
                        sdimy = 'sdimy',
                        spat_enr_names = NULL,
                        cell_annotation_values = NULL,
                        cell_color_gradient = c('blue', 'white', 'red'),
                        gradient_midpoint = NULL,
                        gradient_limits = NULL,
                        select_cell_groups = NULL,
                        select_cells = NULL,
                        point_shape = c('border', 'no_border'),
                        point_size = 3,
                        point_border_col = 'black',
                        point_border_stroke = 0.1,
                        show_cluster_center = F,
                        show_center_label = F,
                        center_point_size = 4,
                        center_point_border_col = 'black',
                        center_point_border_stroke = 0.1,
                        label_size = 4,
                        label_fontface = 'bold',
                        show_network = F,
                        spatial_network_name = 'Delaunay_network',
                        network_color = NULL,
                        network_alpha = 1,
                        show_grid = F,
                        spatial_grid_name = 'spatial_grid',
                        grid_color = NULL,
                        show_other_cells = T,
                        other_cell_color = 'lightgrey',
                        other_point_size = 1,
                        other_cells_alpha = 0.1,
                        coord_fix_ratio = NULL,
                        show_legend = T,
                        legend_text = 8,
                        legend_symbol_size = 1,
                        background_color = 'white',
                        axis_text = 8,
                        axis_title = 8,
                        cow_n_col = 2,
                        cow_rel_h = 1,
                        cow_rel_w = 1,
                        cow_align = 'h',
                        show_plot = NA,
                        return_plot = NA,
                        save_plot = NA,
                        save_param =  list(),
                        default_save_name = 'spatCellPlot') {

  spatCellPlot2D(gobject = gobject,
                 sdimx = sdimx,
                 sdimy = sdimy,
                 spat_enr_names = spat_enr_names,
                 cell_annotation_values = cell_annotation_values,
                 cell_color_gradient = cell_color_gradient,
                 gradient_midpoint = gradient_midpoint,
                 gradient_limits = gradient_limits,
                 select_cell_groups = select_cell_groups,
                 select_cells = select_cells,
                 point_shape = point_shape,
                 point_size = point_size,
                 point_border_col = point_border_col,
                 point_border_stroke = point_border_stroke,
                 show_cluster_center = show_cluster_center,
                 show_center_label = show_center_label,
                 center_point_size = center_point_size,
                 center_point_border_col = center_point_border_col,
                 center_point_border_stroke = center_point_border_stroke,
                 label_size = label_size,
                 label_fontface = label_fontface,
                 show_network = show_network,
                 spatial_network_name = spatial_network_name,
                 network_color = network_color,
                 network_alpha = network_alpha,
                 show_grid = show_grid,
                 spatial_grid_name = spatial_grid_name,
                 grid_color = grid_color,
                 show_other_cells = show_other_cells,
                 other_cell_color = other_cell_color,
                 other_point_size = other_point_size,
                 other_cells_alpha = other_cells_alpha,
                 coord_fix_ratio = coord_fix_ratio,
                 show_legend = show_legend,
                 legend_text = legend_text,
                 legend_symbol_size = legend_symbol_size,
                 background_color = background_color,
                 axis_text = axis_text,
                 axis_title = axis_title,
                 cow_n_col = cow_n_col,
                 cow_rel_h = cow_rel_h,
                 cow_rel_w = cow_rel_w,
                 cow_align = cow_align,
                 show_plot = show_plot,
                 return_plot = return_plot,
                 save_plot = save_plot,
                 save_param =  save_param,
                 default_save_name = default_save_name)

}





#' @title dimCellPlot2D
#' @name dimCellPlot2D
#' @description Visualize cells according to dimension reduction coordinates
#' @param gobject giotto object
#' @param dim_reduction_to_use dimension reduction to use
#' @param dim_reduction_name dimension reduction name
#' @param dim1_to_use dimension to use on x-axis
#' @param dim2_to_use dimension to use on y-axis
#' @param spat_enr_names names of spatial enrichment results to include
#' @param cell_annotation_values numeric cell annotation columns
#' @param show_NN_network show underlying NN network
#' @param nn_network_to_use type of NN network to use (kNN vs sNN)
#' @param network_name name of NN network to use, if show_NN_network = TRUE
#' @param cell_color color for cells (see details)
#' @param color_as_factor convert color column to factor
#' @param cell_color_code named vector with colors
#' @param cell_color_gradient vector with 3 colors for numeric data
#' @param gradient_midpoint midpoint for color gradient
#' @param gradient_limits vector with lower and upper limits
#' @param select_cell_groups select subset of cells/clusters based on cell_color parameter
#' @param select_cells select subset of cells based on cell IDs
#' @param show_other_cells display not selected cells
#' @param other_cell_color color of not selected cells
#' @param other_point_size size of not selected cells
#' @param show_cluster_center plot center of selected clusters
#' @param show_center_label plot label of selected clusters
#' @param center_point_size size of center points
#' @param label_size  size of labels
#' @param label_fontface font of labels
#' @param edge_alpha column to use for alpha of the edges
#' @param point_shape point with border or not (border or no_border)
#' @param point_size size of point (cell)
#' @param point_border_col color of border around points
#' @param point_border_stroke stroke size of border around points
#' @param title title for plot, defaults to cell_color parameter
#' @param show_legend show legend
#' @param legend_text size of legend text
#' @param legend_symbol_size size of legend symbols
#' @param background_color color of plot background
#' @param axis_text size of axis text
#' @param axis_title size of axis title
#' @param show_plot show plot
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters from \code{\link{all_plots_save_function}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return ggplot
#' @details Description of parameters. For 3D plots see \code{\link{dimPlot3D}}
#' @export
#' @examples
#'     dimCellPlot2D(gobject)
dimCellPlot2D = function(gobject,
                         dim_reduction_to_use = 'umap',
                         dim_reduction_name = 'umap',
                         dim1_to_use = 1,
                         dim2_to_use = 2,
                         spat_enr_names = NULL,
                         cell_annotation_values = NULL,
                         show_NN_network = F,
                         nn_network_to_use = 'sNN',
                         network_name = 'sNN.pca',
                         cell_color_gradient = c('blue', 'white', 'red'),
                         gradient_midpoint = NULL,
                         gradient_limits = NULL,
                         select_cell_groups = NULL,
                         select_cells = NULL,
                         show_other_cells = T,
                         other_cell_color = 'lightgrey',
                         other_point_size = 0.5,
                         show_cluster_center = F,
                         show_center_label = T,
                         center_point_size = 4,
                         center_point_border_col = 'black',
                         center_point_border_stroke = 0.1,
                         label_size = 4,
                         label_fontface = 'bold',
                         edge_alpha = NULL,
                         point_shape = c('border', 'no_border'),
                         point_size = 1,
                         point_border_col = 'black',
                         point_border_stroke = 0.1,
                         show_legend = T,
                         legend_text = 8,
                         legend_symbol_size = 1,
                         background_color = 'white',
                         axis_text = 8,
                         axis_title = 8,
                         cow_n_col = 2,
                         cow_rel_h = 1,
                         cow_rel_w = 1,
                         cow_align = 'h',
                         show_plot = NA,
                         return_plot = NA,
                         save_plot = NA,
                         save_param = list(),
                         default_save_name = 'dimCellPlot2D') {


  comb_metadata = combineMetadata(gobject = gobject,
                                  spat_enr_names = spat_enr_names)

  # keep only available columns
  possible_value_cols = colnames(comb_metadata)
  if(is.null(cell_annotation_values)) {
    stop('you need to choose which continuous/numerical cell annotations or enrichments you want to visualize')
  }
  cell_annotation_values = cell_annotation_values[cell_annotation_values %in% possible_value_cols]

  # print, return and save parameters
  show_plot = ifelse(is.na(show_plot), readGiottoInstructions(gobject, param = 'show_plot'), show_plot)
  save_plot = ifelse(is.na(save_plot), readGiottoInstructions(gobject, param = 'save_plot'), save_plot)
  return_plot = ifelse(is.na(return_plot), readGiottoInstructions(gobject, param = 'return_plot'), return_plot)

  ## plotting ##
  savelist <- list()

  for(annot in cell_annotation_values) {

    pl = dimPlot2D(gobject = gobject,
                   group_by = NULL,
                   group_by_subset = NULL,
                   dim_reduction_to_use = dim_reduction_to_use,
                   dim_reduction_name = dim_reduction_name,
                   dim1_to_use = dim1_to_use,
                   dim2_to_use = dim2_to_use,
                   spat_enr_names = spat_enr_names,
                   show_NN_network = show_NN_network,
                   nn_network_to_use = nn_network_to_use,
                   network_name = network_name,
                   cell_color = annot,
                   color_as_factor = FALSE,
                   cell_color_gradient = cell_color_gradient,
                   gradient_midpoint = gradient_midpoint,
                   gradient_limits = gradient_limits,
                   select_cell_groups = select_cell_groups,
                   select_cells = select_cells,
                   show_other_cells = show_other_cells,
                   other_cell_color = other_cell_color,
                   other_point_size = other_point_size,
                   show_cluster_center = show_cluster_center,
                   show_center_label = show_center_label,
                   center_point_size = center_point_size,
                   center_point_border_col = center_point_border_col,
                   center_point_border_stroke = center_point_border_stroke,
                   label_size = label_size,
                   label_fontface = label_fontface,
                   edge_alpha = edge_alpha,
                   point_shape = point_shape,
                   point_size = point_size,
                   point_border_col = point_border_col,
                   point_border_stroke = point_border_stroke,
                   title = annot,
                   show_legend = show_legend,
                   legend_text = legend_text,
                   legend_symbol_size = legend_symbol_size,
                   background_color = background_color,
                   axis_text = axis_text,
                   axis_title = axis_title,
                   show_plot = FALSE,
                   return_plot = TRUE,
                   save_plot = FALSE,
                   save_param = list(),
                   default_save_name = 'dimPlot2D')


    savelist[[annot]] <- pl

  }


  # combine plots with cowplot
  combo_plot <- cowplot::plot_grid(plotlist = savelist,
                                   ncol = cow_n_col,
                                   rel_heights = cow_rel_h,
                                   rel_widths = cow_rel_w,
                                   align = cow_align)


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




#' @title dimCellPlot
#' @name dimCellPlot
#' @description Visualize cells according to dimension reduction coordinates
#' @param gobject giotto object
#' @param dim_reduction_to_use dimension reduction to use
#' @param dim_reduction_name dimension reduction name
#' @param dim1_to_use dimension to use on x-axis
#' @param dim2_to_use dimension to use on y-axis
#' @param spat_enr_names names of spatial enrichment results to include
#' @param cell_annotation_values numeric cell annotation columns
#' @param show_NN_network show underlying NN network
#' @param nn_network_to_use type of NN network to use (kNN vs sNN)
#' @param network_name name of NN network to use, if show_NN_network = TRUE
#' @param cell_color color for cells (see details)
#' @param color_as_factor convert color column to factor
#' @param cell_color_code named vector with colors
#' @param cell_color_gradient vector with 3 colors for numeric data
#' @param gradient_midpoint midpoint for color gradient
#' @param gradient_limits vector with lower and upper limits
#' @param select_cell_groups select subset of cells/clusters based on cell_color parameter
#' @param select_cells select subset of cells based on cell IDs
#' @param show_other_cells display not selected cells
#' @param other_cell_color color of not selected cells
#' @param other_point_size size of not selected cells
#' @param show_cluster_center plot center of selected clusters
#' @param show_center_label plot label of selected clusters
#' @param center_point_size size of center points
#' @param label_size  size of labels
#' @param label_fontface font of labels
#' @param edge_alpha column to use for alpha of the edges
#' @param point_shape point with border or not (border or no_border)
#' @param point_size size of point (cell)
#' @param point_border_col color of border around points
#' @param point_border_stroke stroke size of border around points
#' @param title title for plot, defaults to cell_color parameter
#' @param show_legend show legend
#' @param legend_text size of legend text
#' @param legend_symbol_size size of legend symbols
#' @param background_color color of plot background
#' @param axis_text size of axis text
#' @param axis_title size of axis title
#' @param show_plot show plot
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters from \code{\link{all_plots_save_function}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return ggplot
#' @details Description of parameters. For 3D plots see \code{\link{dimCellPlot2D}}
#' @export
#' @examples
#'     dimCellPlot(gobject)
dimCellPlot = function(gobject,
                       dim_reduction_to_use = 'umap',
                       dim_reduction_name = 'umap',
                       dim1_to_use = 1,
                       dim2_to_use = 2,
                       spat_enr_names = NULL,
                       cell_annotation_values = NULL,
                       show_NN_network = F,
                       nn_network_to_use = 'sNN',
                       network_name = 'sNN.pca',
                       cell_color_gradient = c('blue', 'white', 'red'),
                       gradient_midpoint = NULL,
                       gradient_limits = NULL,
                       select_cell_groups = NULL,
                       select_cells = NULL,
                       show_other_cells = T,
                       other_cell_color = 'lightgrey',
                       other_point_size = 0.5,
                       show_cluster_center = F,
                       show_center_label = T,
                       center_point_size = 4,
                       center_point_border_col = 'black',
                       center_point_border_stroke = 0.1,
                       label_size = 4,
                       label_fontface = 'bold',
                       edge_alpha = NULL,
                       point_shape = c('border', 'no_border'),
                       point_size = 1,
                       point_border_col = 'black',
                       point_border_stroke = 0.1,
                       show_legend = T,
                       legend_text = 8,
                       legend_symbol_size = 1,
                       background_color = 'white',
                       axis_text = 8,
                       axis_title = 8,
                       cow_n_col = 2,
                       cow_rel_h = 1,
                       cow_rel_w = 1,
                       cow_align = 'h',
                       show_plot = NA,
                       return_plot = NA,
                       save_plot = NA,
                       save_param = list(),
                       default_save_name = 'dimCellPlot') {

  dimCellPlot2D(gobject = gobject,
                dim_reduction_to_use = dim_reduction_to_use,
                dim_reduction_name = dim_reduction_name,
                dim1_to_use = dim1_to_use,
                dim2_to_use = dim2_to_use,
                spat_enr_names = spat_enr_names,
                cell_annotation_values = cell_annotation_values,
                show_NN_network = show_NN_network,
                nn_network_to_use = nn_network_to_use,
                network_name = network_name,
                cell_color_gradient = cell_color_gradient,
                gradient_midpoint = gradient_midpoint,
                gradient_limits = gradient_limits,
                select_cell_groups = select_cell_groups,
                select_cells = select_cells,
                show_other_cells = show_other_cells,
                other_cell_color = other_cell_color,
                other_point_size = other_point_size,
                show_cluster_center = show_cluster_center,
                show_center_label = show_center_label,
                center_point_size = center_point_size,
                center_point_border_col = center_point_border_col,
                center_point_border_stroke = center_point_border_stroke,
                label_size = label_size,
                label_fontface = label_fontface,
                edge_alpha = edge_alpha,
                point_shape = point_shape,
                point_size = point_size,
                point_border_col = point_border_col,
                point_border_stroke = point_border_stroke,
                show_legend = show_legend,
                legend_text = legend_text,
                legend_symbol_size = legend_symbol_size,
                background_color = background_color,
                axis_text = axis_text,
                axis_title = axis_title,
                cow_n_col = cow_n_col,
                cow_rel_h = cow_rel_h,
                cow_rel_w = cow_rel_w,
                cow_align = cow_align,
                show_plot = show_plot,
                return_plot = return_plot,
                save_plot = save_plot,
                save_param = save_param,
                default_save_name = default_save_name)

}




#' @title spatDimCellPlot2D
#' @name spatDimCellPlot2D
#' @description Visualize numerical features of cells according to spatial AND dimension reduction coordinates in 2D
#' @param gobject giotto object
#' @param plot_alignment direction to align plot
#' @param spat_enr_names names of spatial enrichment results to include
#' @param cell_annotation_values numeric cell annotation columns
#' @param dim_reduction_to_use dimension reduction to use
#' @param dim_reduction_name dimension reduction name
#' @param dim1_to_use dimension to use on x-axis
#' @param dim2_to_use dimension to use on y-axis
#' @param sdimx = spatial dimension to use on x-axis
#' @param sdimy = spatial dimension to use on y-axis
#' @param cell_color_gradient vector with 3 colors for numeric data
#' @param gradient_midpoint midpoint for color gradient
#' @param gradient_limits vector with lower and upper limits
#' @param select_cell_groups select subset of cells/clusters based on cell_color parameter
#' @param select_cells select subset of cells based on cell IDs
#' @param dim_point_shape dim reduction points with border or not (border or no_border)
#' @param dim_point_size size of points in dim. reduction space
#' @param dim_point_border_col border color of points in dim. reduction space
#' @param dim_point_border_stroke border stroke of points in dim. reduction space
#' @param spat_point_shape spatial points with border or not (border or no_border)
#' @param spat_point_size size of spatial points
#' @param spat_point_border_col border color of spatial points
#' @param spat_point_border_stroke border stroke of spatial points
#' @param dim_show_cluster_center show the center of each cluster
#' @param dim_show_center_label provide a label for each cluster
#' @param dim_center_point_size size of the center point
#' @param dim_center_point_border_col border color of center point
#' @param dim_center_point_border_stroke stroke size of center point
#' @param dim_label_size size of the center label
#' @param dim_label_fontface font of the center label
#' @param spat_show_cluster_center show the center of each cluster
#' @param spat_show_center_label provide a label for each cluster
#' @param spat_center_point_size size of the center point
#' @param spat_label_size size of the center label
#' @param spat_label_fontface font of the center label
#' @param show_NN_network show underlying NN network
#' @param nn_network_to_use type of NN network to use (kNN vs sNN)
#' @param nn_network_name name of NN network to use, if show_NN_network = TRUE
#' @param dim_edge_alpha column to use for alpha of the edges
#' @param spat_show_network show spatial network
#' @param spatial_network_name name of spatial network to use
#' @param spat_network_color color of spatial network
#' @param spat_show_grid show spatial grid
#' @param spatial_grid_name name of spatial grid to use
#' @param spat_grid_color color of spatial grid
#' @param show_other_cells display not selected cells
#' @param other_cell_color color of not selected cells
#' @param dim_other_point_size size of not selected dim cells
#' @param spat_other_point_size size of not selected spat cells
#' @param spat_other_cells_alpha alpha of not selected spat cells
#' @param coord_fix_ratio ratio for coordinates
#' @param cow_n_col cowplot param: how many columns
#' @param cow_rel_h cowplot param: relative height
#' @param cow_rel_w cowplot param: relative width
#' @param cow_align cowplot param: how to align
#' @param show_legend show legend
#' @param legend_text size of legend text
#' @param legend_symbol_size size of legend symbols
#' @param dim_background_color background color of points in dim. reduction space
#' @param spat_background_color background color of spatial points
#' @param axis_text size of axis text
#' @param axis_title size of axis title
#' @param show_plot show plot
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters from \code{\link{all_plots_save_function}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return ggplot
#' @details Description of parameters.
#' @export
#' @examples
#'     spatDimCellPlot2D(gobject)
spatDimCellPlot2D <- function(gobject,
                              plot_alignment = c('vertical', 'horizontal'),
                              spat_enr_names = NULL,
                              cell_annotation_values = NULL,
                              dim_reduction_to_use = 'umap',
                              dim_reduction_name = 'umap',
                              dim1_to_use = 1,
                              dim2_to_use = 2,
                              sdimx = 'sdimx',
                              sdimy = 'sdimy',
                              cell_color_gradient = c('blue', 'white', 'red'),
                              gradient_midpoint = NULL,
                              gradient_limits = NULL,
                              select_cell_groups = NULL,
                              select_cells = NULL,
                              dim_point_shape = c('border', 'no_border'),
                              dim_point_size = 1,
                              dim_point_border_col = 'black',
                              dim_point_border_stroke = 0.1,
                              spat_point_shape = c('border', 'no_border'),
                              spat_point_size = 1,
                              spat_point_border_col = 'black',
                              spat_point_border_stroke = 0.1,
                              dim_show_cluster_center = F,
                              dim_show_center_label = T,
                              dim_center_point_size = 4,
                              dim_center_point_border_col = 'black',
                              dim_center_point_border_stroke = 0.1,
                              dim_label_size = 4,
                              dim_label_fontface = 'bold',
                              spat_show_cluster_center = F,
                              spat_show_center_label = F,
                              spat_center_point_size = 4,
                              spat_center_point_border_col = 'black',
                              spat_center_point_border_stroke = 0.1,
                              spat_label_size = 4,
                              spat_label_fontface = 'bold',
                              show_NN_network = F,
                              nn_network_to_use = 'sNN',
                              nn_network_name = 'sNN.pca',
                              dim_edge_alpha = 0.5,
                              spat_show_network = F,
                              spatial_network_name = 'Delaunay_network',
                              spat_network_color = 'red',
                              spat_network_alpha = 0.5,
                              spat_show_grid = F,
                              spatial_grid_name = 'spatial_grid',
                              spat_grid_color = 'green',
                              show_other_cells = TRUE,
                              other_cell_color = 'grey',
                              dim_other_point_size = 0.5,
                              spat_other_point_size = 0.5,
                              spat_other_cells_alpha = 0.5,
                              show_legend = T,
                              legend_text = 8,
                              legend_symbol_size = 1,
                              dim_background_color = 'white',
                              spat_background_color = 'white',
                              axis_text = 8,
                              axis_title = 8,
                              coord_fix_ratio = NULL,
                              cow_n_col = 2,
                              cow_rel_h = 1,
                              cow_rel_w = 1,
                              cow_align = 'h',
                              show_plot = NA,
                              return_plot = NA,
                              save_plot = NA,
                              save_param =  list(),
                              default_save_name = 'spatDimCellPlot2D') {

  plot_alignment = match.arg(plot_alignment, choices = c('vertical', 'horizontal'))

  # dimension reduction plot
  dmpl = dimCellPlot2D(gobject = gobject,
                       dim_reduction_to_use = dim_reduction_to_use,
                       dim_reduction_name = dim_reduction_name,
                       dim1_to_use = dim1_to_use,
                       dim2_to_use = dim2_to_use,
                       spat_enr_names = spat_enr_names,
                       cell_annotation_values = cell_annotation_values,
                       cell_color_gradient = cell_color_gradient,
                       gradient_midpoint = gradient_midpoint,
                       gradient_limits = gradient_limits,
                       select_cell_groups = select_cell_groups,
                       select_cells = select_cells,
                       point_shape = dim_point_shape,
                       point_size = dim_point_size,
                       point_border_col = dim_point_border_col,
                       point_border_stroke = dim_point_border_stroke,
                       show_cluster_center = dim_show_cluster_center,
                       show_center_label = dim_show_center_label,
                       center_point_size = dim_center_point_size,
                       center_point_border_col = dim_center_point_border_col,
                       center_point_border_stroke = dim_center_point_border_stroke,
                       label_size = dim_label_size,
                       label_fontface = dim_label_fontface,
                       show_NN_network = show_NN_network,
                       nn_network_to_use = nn_network_to_use,
                       network_name = nn_network_name,
                       edge_alpha = dim_edge_alpha,
                       show_other_cells = show_other_cells,
                       other_cell_color = other_cell_color,
                       other_point_size = dim_other_point_size,
                       show_legend = show_legend,
                       legend_text = legend_text,
                       legend_symbol_size = legend_symbol_size,
                       background_color = dim_background_color,
                       axis_text = axis_text,
                       axis_title = axis_title,
                       cow_n_col = cow_n_col,
                       cow_rel_h = cow_rel_h,
                       cow_rel_w = cow_rel_w,
                       cow_align = cow_align,
                       show_plot = FALSE,
                       return_plot = TRUE,
                       save_plot = FALSE)

  # spatial plot
  spl = spatCellPlot2D(gobject = gobject,
                       sdimx = sdimx,
                       sdimy = sdimy,
                       spat_enr_names = spat_enr_names,
                       cell_annotation_values = cell_annotation_values,
                       cell_color_gradient = cell_color_gradient,
                       gradient_midpoint = gradient_midpoint,
                       gradient_limits = gradient_limits,
                       select_cell_groups = select_cell_groups,
                       select_cells = select_cells,
                       point_shape = spat_point_shape,
                       point_size = spat_point_size,
                       point_border_col = spat_point_border_col,
                       point_border_stroke = spat_point_border_stroke,
                       show_cluster_center = spat_show_cluster_center,
                       show_center_label = spat_show_center_label,
                       center_point_size = spat_center_point_size,
                       center_point_border_col = spat_center_point_border_col,
                       center_point_border_stroke = spat_center_point_border_stroke,
                       label_size = spat_label_size,
                       label_fontface = spat_label_fontface,
                       show_network = spat_show_network,
                       spatial_network_name = spatial_network_name,
                       network_color = spat_network_color,
                       network_alpha = spat_network_alpha,
                       show_grid = spat_show_grid,
                       spatial_grid_name = spatial_grid_name,
                       grid_color = spat_grid_color,
                       show_other_cells = show_other_cells,
                       other_cell_color = other_cell_color,
                       other_point_size = spat_other_point_size,
                       other_cells_alpha = spat_other_cells_alpha,
                       coord_fix_ratio = coord_fix_ratio,
                       show_legend = show_legend,
                       legend_text = legend_text,
                       legend_symbol_size = legend_symbol_size,
                       background_color = spat_background_color,
                       axis_text = axis_text,
                       axis_title = axis_title,
                       cow_n_col = cow_n_col,
                       cow_rel_h = cow_rel_h,
                       cow_rel_w = cow_rel_w,
                       cow_align = cow_align,
                       show_plot = FALSE,
                       return_plot = TRUE,
                       save_plot = FALSE)


  if(plot_alignment == 'vertical') {
    ncol = 1
    nrow = 2
    combo_plot = cowplot::plot_grid(dmpl, spl, ncol = ncol, nrow = nrow, rel_heights = c(1), rel_widths = c(1), align = 'v')
  } else {
    ncol = 2
    nrow = 1
    combo_plot = cowplot::plot_grid(dmpl, spl, ncol = ncol, nrow = nrow, rel_heights = c(1), rel_widths = c(1), align = 'h')
  }

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




#' @title spatDimCellPlot
#' @name spatDimCellPlot
#' @description Visualize numerical features of cells according to spatial AND dimension reduction coordinates in 2D
#' @param gobject giotto object
#' @param plot_alignment direction to align plot
#' @param spat_enr_names names of spatial enrichment results to include
#' @param cell_annotation_values numeric cell annotation columns
#' @param dim_reduction_to_use dimension reduction to use
#' @param dim_reduction_name dimension reduction name
#' @param dim1_to_use dimension to use on x-axis
#' @param dim2_to_use dimension to use on y-axis
#' @param sdimx = spatial dimension to use on x-axis
#' @param sdimy = spatial dimension to use on y-axis
#' @param cell_color_gradient vector with 3 colors for numeric data
#' @param gradient_midpoint midpoint for color gradient
#' @param gradient_limits vector with lower and upper limits
#' @param select_cell_groups select subset of cells/clusters based on cell_color parameter
#' @param select_cells select subset of cells based on cell IDs
#' @param dim_point_shape spatial points with border or not (border or no_border)
#' @param dim_point_size size of points in dim. reduction space
#' @param dim_point_border_col border color of points in dim. reduction space
#' @param dim_point_border_stroke border stroke of points in dim. reduction space
#' @param spat_point_shape spatial points with border or not (border or no_border)
#' @param spat_point_size size of spatial points
#' @param spat_point_border_col border color of spatial points
#' @param spat_point_border_stroke border stroke of spatial points
#' @param dim_show_cluster_center show the center of each cluster
#' @param dim_show_center_label provide a label for each cluster
#' @param dim_center_point_size size of the center point
#' @param dim_center_point_border_col border color of center point
#' @param dim_center_point_border_stroke stroke size of center point
#' @param dim_label_size size of the center label
#' @param dim_label_fontface font of the center label
#' @param spat_show_cluster_center show the center of each cluster
#' @param spat_show_center_label provide a label for each cluster
#' @param spat_center_point_size size of the center point
#' @param spat_label_size size of the center label
#' @param spat_label_fontface font of the center label
#' @param show_NN_network show underlying NN network
#' @param nn_network_to_use type of NN network to use (kNN vs sNN)
#' @param nn_network_name name of NN network to use, if show_NN_network = TRUE
#' @param dim_edge_alpha column to use for alpha of the edges
#' @param spat_show_network show spatial network
#' @param spatial_network_name name of spatial network to use
#' @param spat_network_color color of spatial network
#' @param spat_show_grid show spatial grid
#' @param spatial_grid_name name of spatial grid to use
#' @param spat_grid_color color of spatial grid
#' @param show_other_cells display not selected cells
#' @param other_cell_color color of not selected cells
#' @param dim_other_point_size size of not selected dim cells
#' @param spat_other_point_size size of not selected spat cells
#' @param spat_other_cells_alpha alpha of not selected spat cells
#' @param show_legend show legend
#' @param legend_text size of legend text
#' @param legend_symbol_size size of legend symbols
#' @param dim_background_color background color of points in dim. reduction space
#' @param spat_background_color background color of spatial points
#' @param axis_text size of axis text
#' @param axis_title size of axis title
#' @param coord_fix_ratio ratio for coordinates
#' @param cow_n_col cowplot param: how many columns
#' @param cow_rel_h cowplot param: relative height
#' @param cow_rel_w cowplot param: relative width
#' @param cow_align cowplot param: how to align
#' @param show_plot show plot
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters from \code{\link{all_plots_save_function}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return ggplot
#' @details Description of parameters.
#' @export
#' @examples
#'     spatDimCellPlot(gobject)
spatDimCellPlot = function(gobject,
                           plot_alignment = c('vertical', 'horizontal'),
                           spat_enr_names = NULL,
                           cell_annotation_values = NULL,
                           dim_reduction_to_use = 'umap',
                           dim_reduction_name = 'umap',
                           dim1_to_use = 1,
                           dim2_to_use = 2,
                           sdimx = 'sdimx',
                           sdimy = 'sdimy',
                           cell_color_gradient = c('blue', 'white', 'red'),
                           gradient_midpoint = NULL,
                           gradient_limits = NULL,
                           select_cell_groups = NULL,
                           select_cells = NULL,
                           dim_point_shape = c('border', 'no_border'),
                           dim_point_size = 1,
                           dim_point_border_col = 'black',
                           dim_point_border_stroke = 0.1,
                           spat_point_shape = c('border', 'no_border'),
                           spat_point_size = 1,
                           spat_point_border_col = 'black',
                           spat_point_border_stroke = 0.1,
                           dim_show_cluster_center = F,
                           dim_show_center_label = T,
                           dim_center_point_size = 4,
                           dim_center_point_border_col = 'black',
                           dim_center_point_border_stroke = 0.1,
                           dim_label_size = 4,
                           dim_label_fontface = 'bold',
                           spat_show_cluster_center = F,
                           spat_show_center_label = F,
                           spat_center_point_size = 4,
                           spat_center_point_border_col = 'black',
                           spat_center_point_border_stroke = 0.1,
                           spat_label_size = 4,
                           spat_label_fontface = 'bold',
                           show_NN_network = F,
                           nn_network_to_use = 'sNN',
                           nn_network_name = 'sNN.pca',
                           dim_edge_alpha = 0.5,
                           spat_show_network = F,
                           spatial_network_name = 'Delaunay_network',
                           spat_network_color = 'red',
                           spat_network_alpha = 0.5,
                           spat_show_grid = F,
                           spatial_grid_name = 'spatial_grid',
                           spat_grid_color = 'green',
                           show_other_cells = TRUE,
                           other_cell_color = 'grey',
                           dim_other_point_size = 0.5,
                           spat_other_point_size = 0.5,
                           spat_other_cells_alpha = 0.5,
                           coord_fix_ratio = NULL,
                           cow_n_col = 2,
                           cow_rel_h = 1,
                           cow_rel_w = 1,
                           cow_align = 'h',
                           show_legend = T,
                           legend_text = 8,
                           legend_symbol_size = 1,
                           dim_background_color = 'white',
                           spat_background_color = 'white',
                           axis_text = 8,
                           axis_title = 8,
                           show_plot = NA,
                           return_plot = NA,
                           save_plot = NA,
                           save_param =  list(),
                           default_save_name = 'spatDimCellPlot') {

  spatDimCellPlot2D(gobject = gobject,
                    plot_alignment = plot_alignment,
                    spat_enr_names = spat_enr_names,
                    cell_annotation_values = cell_annotation_values,
                    dim_reduction_to_use = dim_reduction_to_use,
                    dim_reduction_name = dim_reduction_name,
                    dim1_to_use = dim1_to_use,
                    dim2_to_use = dim2_to_use,
                    sdimx = sdimx,
                    sdimy = sdimy,
                    cell_color_gradient = cell_color_gradient,
                    gradient_midpoint = gradient_midpoint,
                    gradient_limits = gradient_limits,
                    select_cell_groups = select_cell_groups,
                    select_cells = select_cells,
                    dim_point_shape = dim_point_shape,
                    dim_point_size = dim_point_size,
                    dim_point_border_col = dim_point_border_col,
                    dim_point_border_stroke = dim_point_border_stroke,
                    spat_point_shape = spat_point_shape,
                    spat_point_size = spat_point_size,
                    spat_point_border_col = spat_point_border_col,
                    spat_point_border_stroke = spat_point_border_stroke,
                    dim_show_cluster_center = dim_show_cluster_center,
                    dim_show_center_label = dim_show_center_label,
                    dim_center_point_size = dim_center_point_size,
                    dim_center_point_border_col = dim_center_point_border_col,
                    dim_center_point_border_stroke = dim_center_point_border_stroke,
                    dim_label_size = dim_label_size,
                    dim_label_fontface = dim_label_fontface,
                    spat_show_cluster_center = spat_show_cluster_center,
                    spat_show_center_label = spat_show_center_label,
                    spat_center_point_size = spat_center_point_size,
                    spat_center_point_border_col = spat_center_point_border_col,
                    spat_center_point_border_stroke = spat_center_point_border_stroke,
                    spat_label_size = spat_label_size,
                    spat_label_fontface = spat_label_fontface,
                    show_NN_network = show_NN_network,
                    nn_network_to_use = nn_network_to_use,
                    nn_network_name = nn_network_name,
                    dim_edge_alpha = dim_edge_alpha,
                    spat_show_network = spat_show_network,
                    spatial_network_name = spatial_network_name,
                    spat_network_color = spat_network_color,
                    spat_network_alpha = spat_network_alpha,
                    spat_show_grid = spat_show_grid,
                    spatial_grid_name = spatial_grid_name,
                    spat_grid_color = spat_grid_color,
                    show_other_cells = show_other_cells,
                    other_cell_color = other_cell_color,
                    dim_other_point_size = dim_other_point_size,
                    spat_other_point_size = spat_other_point_size,
                    spat_other_cells_alpha = spat_other_cells_alpha,
                    coord_fix_ratio = coord_fix_ratio,
                    cow_n_col = cow_n_col,
                    cow_rel_h = cow_rel_h,
                    cow_rel_w = cow_rel_w,
                    cow_align = cow_align,
                    show_legend = show_legend,
                    legend_text = legend_text,
                    legend_symbol_size = legend_symbol_size,
                    dim_background_color = dim_background_color,
                    spat_background_color = spat_background_color,
                    axis_text = axis_text,
                    axis_title = axis_title,
                    show_plot = show_plot,
                    return_plot = return_plot,
                    save_plot = save_plot,
                    save_param =  save_param,
                    default_save_name = default_save_name)

}







# ** ####
## 3-D plotly ####
## ----------- ##
#' @title dimPlot3D
#' @name dimPlot3D
#' @description Visualize cells according to dimension reduction coordinates
#' @param gobject giotto object
#' @param dim_reduction_to_use dimension reduction to use
#' @param dim_reduction_name dimension reduction name
#' @param dim1_to_use dimension to use on x-axis
#' @param dim2_to_use dimension to use on y-axis
#' @param dim3_to_use dimension to use on z-axis
#' @param show_NN_network show underlying NN network
#' @param nn_network_to_use type of NN network to use (kNN vs sNN)
#' @param network_name name of NN network to use, if show_NN_network = TRUE
#' @param cell_color color for cells (see details)
#' @param color_as_factor convert color column to factor
#' @param cell_color_code named vector with colors
#' @param select_cell_groups select subset of cells/clusters based on cell_color parameter
#' @param select_cells select subset of cells based on cell IDs
#' @param show_other_cells display not selected cells
#' @param other_cell_color color of not selected cells
#' @param other_point_size size of not selected cells
#' @param show_cluster_center plot center of selected clusters
#' @param show_center_label plot label of selected clusters
#' @param center_point_size size of center points
#' @param label_size  size of labels
#' @param edge_alpha column to use for alpha of the edges
#' @param point_size size of point (cell)
#' @param show_legend show legend
#' @param show_plot show plot
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters from \code{\link{all_plots_save_function}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return plotly
#' @details Description of parameters.
#' @export
#' @examples
#'     dimPlot3D(gobject)
#'
dimPlot3D = function(gobject,
                     dim_reduction_to_use = 'umap',
                     dim_reduction_name = 'umap',
                     dim1_to_use = 1,
                     dim2_to_use = 2,
                     dim3_to_use = 3,

                     select_cell_groups = NULL,
                     select_cells = NULL,
                     show_other_cells = T,
                     other_cell_color = 'lightgrey',
                     other_point_size = 2,

                     show_NN_network = F,
                     nn_network_to_use = 'sNN',
                     network_name = 'sNN.pca',
                     color_as_factor = T,
                     cell_color = NULL,
                     cell_color_code = NULL,
                     show_cluster_center = F,
                     show_center_label = T,
                     center_point_size = 4,
                     label_size = 4,
                     edge_alpha = NULL,
                     point_size = 3,

                     show_plot = NA,
                     return_plot = NA,
                     save_plot = NA,
                     save_param =  list(),
                     default_save_name = "dim3D"){
  if(is.null(dim3_to_use)){
    cat('create 2D plot\n')

    pl = visDimPlot_2D_plotly(gobject = gobject,
                              dim_reduction_to_use = dim_reduction_to_use,
                              dim_reduction_name = dim_reduction_name,
                              dim1_to_use = dim1_to_use,
                              dim2_to_use = dim2_to_use,

                              select_cell_groups = select_cell_groups,
                              select_cells = select_cells,
                              show_other_cells = show_other_cells,
                              other_cell_color = other_cell_color,
                              other_point_size = other_point_size,

                              show_NN_network = show_NN_network,
                              nn_network_to_use = nn_network_to_use,
                              network_name = network_name,
                              color_as_factor = color_as_factor,
                              cell_color = cell_color,
                              cell_color_code = cell_color_code,
                              show_cluster_center = show_cluster_center,
                              show_center_label = show_center_label,
                              center_point_size = center_point_size,
                              label_size = label_size,
                              edge_alpha = edge_alpha,
                              point_size = point_size)
  }

  else{
    cat('create 3D plot\n')
    pl = visDimPlot_3D_plotly(gobject = gobject,
                              dim_reduction_to_use = dim_reduction_to_use,
                              dim_reduction_name = dim_reduction_name,
                              dim1_to_use = dim1_to_use,
                              dim2_to_use = dim2_to_use,
                              dim3_to_use = dim3_to_use,

                              select_cell_groups = select_cell_groups,
                              select_cells = select_cells,
                              show_other_cells = show_other_cells,
                              other_cell_color = other_cell_color,
                              other_point_size = other_point_size,

                              show_NN_network = show_NN_network,
                              nn_network_to_use = nn_network_to_use,
                              network_name = network_name,
                              color_as_factor = color_as_factor,
                              cell_color = cell_color,
                              cell_color_code = cell_color_code,
                              show_cluster_center = show_cluster_center,
                              show_center_label = show_center_label,
                              center_point_size = center_point_size,
                              label_size = label_size,
                              edge_alpha = edge_alpha,
                              point_size = point_size)
  }


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


#' @title plotUMAP_3D
#' @name plotUMAP_3D
#' @description Visualize cells according to dimension reduction coordinates
#' @param gobject giotto object
#' @param dim1_to_use dimension to use on x-axis
#' @param dim2_to_use dimension to use on y-axis
#' @param dim3_to_use dimension to use on z-axis
#' @param dim_reduction_name umap dimension reduction name
#' @param show_NN_network show underlying NN network
#' @param nn_network_to_use type of NN network to use (kNN vs sNN)
#' @param network_name name of NN network to use, if show_NN_network = TRUE
#' @param cell_color color for cells (see details)
#' @param color_as_factor convert color column to factor
#' @param cell_color_code named vector with colors
#' @param select_cell_groups select subset of cells/clusters based on cell_color parameter
#' @param select_cells select subset of cells based on cell IDs
#' @param show_other_cells display not selected cells
#' @param other_cell_color color of not selected cells
#' @param other_point_size size of not selected cells
#' @param show_cluster_center plot center of selected clusters
#' @param show_center_label plot label of selected clusters
#' @param center_point_size size of center points
#' @param label_size  size of labels
#' @param edge_alpha column to use for alpha of the edges
#' @param point_size size of point (cell)
#' @param show_legend show legend
#' @param show_plot show plot
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters from \code{\link{all_plots_save_function}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return plotly
#' @details Description of parameters.
#' @export
#' @examples
#'     plotUMAP_3D(gobject)
#'
plotUMAP_3D = function(gobject, dim_reduction_name = 'umap', default_save_name = 'UMAP_3D', ...) {

  dimPlot3D(gobject = gobject,
            dim_reduction_to_use = 'umap', dim_reduction_name = dim_reduction_name,
            default_save_name = default_save_name, ...)

}


#' @title plotTSNE_3D
#' @name plotTSNE_3D
#' @description Visualize cells according to dimension reduction coordinates
#' @param gobject giotto object
#' @param dim1_to_use dimension to use on x-axis
#' @param dim2_to_use dimension to use on y-axis
#' @param dim3_to_use dimension to use on z-axis
#' @param dim_reduction_name tsne dimension reduction name
#' @param show_NN_network show underlying NN network
#' @param nn_network_to_use type of NN network to use (kNN vs sNN)
#' @param network_name name of NN network to use, if show_NN_network = TRUE
#' @param cell_color color for cells (see details)
#' @param color_as_factor convert color column to factor
#' @param cell_color_code named vector with colors
#' @param select_cell_groups select subset of cells/clusters based on cell_color parameter
#' @param select_cells select subset of cells based on cell IDs
#' @param show_other_cells display not selected cells
#' @param other_cell_color color of not selected cells
#' @param other_point_size size of not selected cells
#' @param show_cluster_center plot center of selected clusters
#' @param show_center_label plot label of selected clusters
#' @param center_point_size size of center points
#' @param label_size  size of labels
#' @param edge_alpha column to use for alpha of the edges
#' @param point_size size of point (cell)
#' @param show_legend show legend
#' @param show_plot show plot
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters from \code{\link{all_plots_save_function}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return plotly
#' @details Description of parameters.
#' @export
#' @examples
#'     plotTSNE_3D(gobject)
#'
plotTSNE_3D = function(gobject, dim_reduction_name = 'tsne', default_save_name = 'TSNE_3D', ...) {

  dimPlot3D(gobject = gobject,
            dim_reduction_to_use = 'tsne', dim_reduction_name = dim_reduction_name,
            default_save_name = default_save_name, ...)

}


#' @title plotPCA_3D
#' @name plotPCA_3D
#' @description Visualize cells according to 3D PCA dimension reduction
#' @param gobject giotto object
#' @param dim1_to_use dimension to use on x-axis
#' @param dim2_to_use dimension to use on y-axis
#' @param dim3_to_use dimension to use on z-axis
#' @param dim_reduction_name pca dimension reduction name
#' @param show_NN_network show underlying NN network
#' @param nn_network_to_use type of NN network to use (kNN vs sNN)
#' @param network_name name of NN network to use, if show_NN_network = TRUE
#' @param cell_color color for cells (see details)
#' @param color_as_factor convert color column to factor
#' @param cell_color_code named vector with colors
#' @param select_cell_groups select subset of cells/clusters based on cell_color parameter
#' @param select_cells select subset of cells based on cell IDs
#' @param show_other_cells display not selected cells
#' @param other_cell_color color of not selected cells
#' @param other_point_size size of not selected cells
#' @param show_cluster_center plot center of selected clusters
#' @param show_center_label plot label of selected clusters
#' @param center_point_size size of center points
#' @param label_size  size of labels
#' @param edge_alpha column to use for alpha of the edges
#' @param point_size size of point (cell)
#' @param show_legend show legend
#' @param show_plot show plot
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters from \code{\link{all_plots_save_function}}
#' @param default_save_name default save name for saving, ideally change save_name in save_param
#' @return plotly
#' @details Description of parameters.
#' @export
#' @examples
#'     plotPCA_3D(gobject)
#'
plotPCA_3D = function(gobject, dim_reduction_name = 'pca', default_save_name = 'PCA_3D', ...) {

  dimPlot3D(gobject = gobject,
            dim_reduction_to_use = 'pca', dim_reduction_name = dim_reduction_name,
            default_save_name = default_save_name, ...)

}


#' @title spatPlot3D
#' @name spatPlot3D
#' @description Visualize cells according to spatial coordinates
#' @param gobject giotto object
#' @param sdimx x-axis dimension name (default = 'sdimx')
#' @param sdimy y-axis dimension name (default = 'sdimy')
#' @param sdimz z-axis dimension name (default = 'sdimy')
#' @param point_size size of point (cell)
#' @param cell_color color for cells (see details)
#' @param cell_color_code named vector with colors
#' @param select_cell_groups select subset of cells/clusters based on cell_color parameter
#' @param select_cells select subset of cells based on cell IDs
#' @param show_other_cells display not selected cells
#' @param other_cell_color color of not selected cells
#' @param show_network show underlying spatial network
#' @param network_color color of spatial network
#' @param spatial_network_name name of spatial network to use
#' @param show_grid show spatial grid
#' @param grid_color color of spatial grid
#' @param spatial_grid_name name of spatial grid to use
#' @param title title of plot
#' @param axis_scale the way to scale the axis
#' @param custom_ratio customize the scale of the plot
#' @param x_ticks set the number of ticks on the x-axis
#' @param y_ticks set the number of ticks on the y-axis
#' @param z_ticks set the number of ticks on the z-axis
#' @param show_legend show legend
#' @param show_plot show plot
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters from \code{\link{all_plots_save_function}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return ggplot
#' @details Description of parameters.
#' @export
#' @examples
#'     spatPlot3D(gobject)
#'
spatPlot3D = function(gobject,
                      sdimx = "sdimx",
                      sdimy = "sdimy",
                      sdimz = "sdimz",
                      point_size = 3,
                      cell_color = NULL,
                      cell_color_code = NULL,
                      select_cell_groups = NULL,
                      select_cells = NULL,
                      show_other_cells = T,
                      other_cell_color = "lightgrey",
                      other_point_size = 0.5,
                      show_network = F,
                      network_color = NULL,
                      network_alpha = 1,
                      other_cell_alpha = 0.5,
                      spatial_network_name = 'Delaunay_network',
                      show_grid = F,
                      grid_color = NULL,
                      spatial_grid_name = 'spatial_grid',
                      title = '',
                      show_legend = T,
                      axis_scale = c("cube","real","custom"),
                      custom_ratio = NULL,
                      x_ticks = NULL,
                      y_ticks = NULL,
                      z_ticks = NULL,
                      show_plot = NA,
                      return_plot = NA,
                      save_plot = NA,
                      save_param =  list(),
                      default_save_name = "spat3D") {
  if(is.null(sdimz)){
    cat('create 2D plot\n')

    pl = visPlot_2D_plotly(gobject = gobject,
                           sdimx = sdimx,
                           sdimy = sdimy,
                           point_size = point_size,
                           cell_color = cell_color,
                           cell_color_code = cell_color_code,
                           select_cell_groups = select_cell_groups,
                           select_cells = select_cells,
                           show_other_cells = show_other_cells,
                           other_cell_color = other_cell_color,
                           other_point_size =other_point_size,
                           show_network = show_network,
                           network_color = network_color,
                           network_alpha = network_alpha,
                           other_cell_alpha =other_cell_alpha,
                           spatial_network_name = spatial_network_name,
                           show_grid = show_grid,
                           grid_color = grid_color,
                           grid_alpha = grid_alpha,
                           spatial_grid_name = spatial_grid_name,
                           show_legend = show_legend,
                           axis_scale = axis_scale,
                           custom_ratio = custom_ratio,
                           x_ticks = x_ticks,
                           y_ticks = y_ticks,
                           show_plot = F)
  }
  else{
    cat('create 3D plot\n')
    pl = visPlot_3D_plotly(gobject = gobject,
                           sdimx = sdimx,
                           sdimy = sdimy,
                           sdimz = sdimz,
                           point_size = point_size,
                           cell_color = cell_color,
                           cell_color_code = cell_color_code,
                           select_cell_groups = select_cell_groups,
                           select_cells = select_cells,
                           show_other_cells = show_other_cells,
                           other_cell_color = other_cell_color,
                           other_point_size =other_point_size,
                           show_network = show_network,
                           network_color = network_color,
                           network_alpha = network_alpha,
                           other_cell_alpha =other_cell_alpha,
                           spatial_network_name = spatial_network_name,
                           spatial_grid_name = spatial_grid_name,
                           show_legend = show_legend,
                           axis_scale = axis_scale,
                           custom_ratio = custom_ratio,
                           x_ticks = x_ticks,
                           y_ticks = y_ticks,
                           z_ticks = z_ticks,
                           show_plot = F)
  }

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


#' @title spatDimPlot3D
#' @name spatDimPlot3D
#' @description Visualize cells according to spatial AND dimension reduction coordinates in plotly mode
#' @param gobject giotto object
#' @param plot_alignment direction to align plot
#' @param dim_reduction_to_use dimension reduction to use
#' @param dim_reduction_name dimension reduction name
#' @param dim1_to_use dimension to use on x-axis
#' @param dim2_to_use dimension to use on y-axis
#' @param dim3_to_use dimension to use on z-axis
#' @param sdimx = spatial dimension to use on x-axis
#' @param sdimy = spatial dimension to use on y-axis
#' @param sdimz = spatial dimension to use on z-axis
#' @param show_NN_network show underlying NN network
#' @param nn_network_to_use type of NN network to use (kNN vs sNN)
#' @param network_name name of NN network to use, if show_NN_network = TRUE
#' @param show_cluster_center show the center of each cluster
#' @param show_center_label provide a label for each cluster
#' @param center_point_size size of the center point
#' @param label_size size of the center label
#' @param cell_color color for cells (see details)
#' @param color_as_factor convert color column to factor
#' @param cell_color_code named vector with colors
#' @param select_cell_groups select subset of cells/clusters based on cell_color parameter
#' @param select_cells select subset of cells based on cell IDs
#' @param show_other_cells display not selected cells
#' @param other_cell_color color of not selected cells
#' @param other_point_size size of not selected cells
#' @param dim_point_size size of points in dim. reduction space
#' @param dim_point_border_col border color of points in dim. reduction space
#' @param dim_point_border_stroke border stroke of points in dim. reduction space
#' @param nn_network_alpha column to use for alpha of the edges
#' @param show_spatial_network show spatial network
#' @param spatial_network_name name of spatial network to use
#' @param spatial_network_color color of spatial network
#' @param show_spatial_grid show spatial grid
#' @param spatial_grid_name name of spatial grid to use
#' @param spatial_grid_color color of spatial grid
#' @param spatial_point_size size of spatial points
#' @param spatial_other_point_size size of not selected spatial points
#' @param spatial_network_alpha alpha of spatial network
#' @param spatial_other_cells_alpha alpha of not selected spatial points
#' @param dim_other_point_size size of not selected dim. reduction points
#' @param show_legend show legend
#' @param show_plot show plot
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters from \code{\link{all_plots_save_function}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return plotly
#' @details Description of parameters.
#' @export
#' @examples
#'     spatDimPlot3D(gobject)
spatDimPlot3D <- function(gobject,
                          plot_alignment = c('horizontal','vertical'),
                          dim_reduction_to_use = 'umap',
                          dim_reduction_name = 'umap',
                          dim1_to_use = 1,
                          dim2_to_use = 2,
                          dim3_to_use = 3,
                          sdimx="sdimx",
                          sdimy="sdimy",
                          sdimz="sdimz",
                          show_NN_network = F,
                          nn_network_to_use = 'sNN',
                          network_name = 'sNN.pca',
                          show_cluster_center = F,
                          show_center_label = T,
                          center_point_size = 4,
                          label_size = 16,

                          select_cell_groups = NULL,
                          select_cells = NULL,
                          show_other_cells = T,
                          other_cell_color = 'lightgrey',
                          other_point_size = 1.5,

                          cell_color = NULL,
                          color_as_factor = T,
                          cell_color_code = NULL,
                          dim_point_size = 3,
                          nn_network_alpha = 0.5,
                          show_spatial_network = F,
                          spatial_network_name = 'Delaunay_network',
                          network_color = "lightgray",
                          spatial_network_alpha = 0.5,
                          show_spatial_grid = F,
                          spatial_grid_name = 'spatial_grid',
                          spatial_grid_color = NULL,
                          spatial_grid_alpha = 0.5,
                          spatial_point_size = 3,
                          axis_scale = c("cube","real","custom"),
                          custom_ratio = NULL,
                          x_ticks = NULL,
                          y_ticks = NULL,
                          z_ticks = NULL,
                          legend_text_size = 12,

                          show_plot = NA,
                          return_plot = NA,
                          save_plot = NA,
                          save_param =  list(),
                          default_save_name = "spatDimPlot3D"){

  plot_alignment = match.arg(plot_alignment, choices = c( 'horizontal','vertical'))

  # ********data prepare********#
  ## dimension reduction ##
  dim_dfr = gobject@dimension_reduction$cells[[dim_reduction_to_use]][[dim_reduction_name]]$coordinates[,c(dim1_to_use, dim2_to_use,dim3_to_use)]
  dim_names = colnames(dim_dfr)
  dim_DT = data.table::as.data.table(dim_dfr); dim_DT[, cell_ID := rownames(dim_dfr)]

  ## annotated cell metadata
  cell_metadata = gobject@cell_metadata
  annotated_DT = merge(cell_metadata, dim_DT, by = 'cell_ID')
  annotated_DT = merge(annotated_DT, gobject@spatial_locs,by = 'cell_ID')

  if(dim_reduction_to_use == "pca"){
    eigenvaluesDT = data.table::as.data.table(gobject@dimension_reduction$cells[[dim_reduction_to_use]][[dim_reduction_name]]$misc$eig)
    var_expl_vec = eigenvaluesDT[c(dim1_to_use, dim2_to_use,dim3_to_use)][['percentage of variance']]
    dim1_x_variance = var_expl_vec[1]
    dim2_y_variance = var_expl_vec[2]
    if(!is.null(dim3_to_use)){
      dim3_z_variance = var_expl_vec[3]
    }
  }

  ## nn network
  if(show_NN_network){
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


  ## extract spatial network
  if(!is.null(spatial_network_name)) {
    spatial_network = select_spatialNetwork(gobject,name = spatial_network_name,return_network_Obj = FALSE)
  } else {
    spatial_network = NULL
  }

  ## extract spatial grid
  if(!is.null(spatial_grid_name)) {
    spatial_grid = gobject@spatial_grid[[spatial_grid_name]]
  } else {
    spatial_grid = NULL
  }

  # create matching cell_color_code
  if(is.null(cell_color_code)) {
    if(is.character(cell_color)) {

      cell_metadata = pDataDT(gobject)
      if(cell_color %in% colnames(cell_metadata)) {

        if(color_as_factor == TRUE) {
          number_colors = length(unique(cell_metadata[[cell_color]]))
          cell_color_code = Giotto:::getDistinctColors(n = number_colors)
          names(cell_color_code) = unique(cell_metadata[[cell_color]])
        }
      }
    }
  }

  if(!is.null(select_cells) & !is.null(select_cell_groups)) {
    if(is.null(cell_color)) {
      stop('\n selection of cells is based on cell_color paramter, which is a metadata column \n')
    }
    cat('You have selected both individual cell IDs and a group of cells \n')
    group_cell_IDs = annotated_DT[get(cell_color) %in% select_cell_groups][['cell_ID']]
    select_cells = unique(c(select_cells, group_cell_IDs))
  } else if(!is.null(select_cell_groups)) {
    select_cells = annotated_DT[get(cell_color) %in% select_cell_groups][['cell_ID']]
  }

  if(!is.null(select_cells)) {
    annotated_DT_other = annotated_DT[!annotated_DT$cell_ID %in% select_cells]
    annotated_DT_selected = annotated_DT[annotated_DT$cell_ID %in% select_cells]

    if(show_NN_network == TRUE) {
      annotated_network_DT <- annotated_network_DT[annotated_network_DT$to %in% select_cells & annotated_network_DT$from %in% select_cells]
    }
    if(show_spatial_network == TRUE){
      spatial_network <- spatial_network[spatial_network$to %in% select_cells & spatial_network$from %in% select_cells]
    }

    # if specific cells are selected
    annotated_DT = annotated_DT_selected
  }

  ## if no subsets are required
  if(is.null(select_cells) & is.null(select_cell_groups)) {
    annotated_DT_selected = annotated_DT
    annotated_DT_other    = NULL
  }

  #********** dim plot ***********#
  #2D plot
  if(is.null(dim3_to_use)){
    dpl <- plotly::plot_ly()
    if(show_NN_network == TRUE) {
      if(is.null(nn_network_alpha)) {
        nn_network_alpha = 0.5
      }
      else if(is.character(nn_network_alpha)){
        warning("Edge_alpha for plotly mode is not adjustable yet. Default 0.5 will be set\n")
        nn_network_alpha = 0.5
      }
      dpl <- dpl %>% plotly::add_segments(name = network_name,
                                          type = "scatter",
                                          x = annotated_network_DT[[from_dim_names[1]]],
                                          y = annotated_network_DT[[from_dim_names[2]]],
                                          xend = annotated_network_DT[[to_dim_names[1]]],
                                          yend = annotated_network_DT[[to_dim_names[2]]],
                                          line = list(color = network_color,
                                                      width = 0.5),
                                          opacity=nn_network_alpha)
    }

    if(is.null(cell_color)){
      #cell_color = "lightblue"
      dpl <- dpl %>% plotly::add_trace(type = "scatter",mode = "markers",
                                       x = annotated_DT[[dim_names[1]]],
                                       y = annotated_DT[[dim_names[2]]],
                                       #color = "lightblue",
                                       #colors ="lightblue",
                                       marker = list(size = dim_point_size,
                                                     color = "lightblue"),
                                       showlegend = F)
    }

    else if(cell_color %in% colnames(annotated_DT)){
      if(color_as_factor){
        annotated_DT[[cell_color]] <- as.factor(annotated_DT[[cell_color]])
      }


      dpl <- dpl %>% plotly::add_trace(type = "scatter",mode = "markers",
                                       x = annotated_DT[[dim_names[1]]],
                                       y = annotated_DT[[dim_names[2]]],
                                       color = annotated_DT[[cell_color]],
                                       colors = cell_color_code,
                                       legendgroup = annotated_DT[[cell_color]],
                                       marker = list(size = dim_point_size))
    }

    else{
      stop("cell_color does not exist!\n")
    }


    if((show_cluster_center == TRUE | show_center_label == TRUE)&!is.null(cell_color)) {
      annotated_DT_centers = annotated_DT[, .(center_1 = median(get(dim_names[1])),
                                              center_2 = median(get(dim_names[2]))),
                                          by = cell_color]
      annotated_DT_centers[[cell_color]] <- as.factor(annotated_DT_centers[[cell_color]])
      if(show_cluster_center == TRUE){
        dpl <- dpl %>% plotly::add_trace(type = "scatter",mode = "markers",
                                         x = annotated_DT_centers[["center_1"]],
                                         y = annotated_DT_centers[["center_2"]],
                                         color = annotated_DT_centers[[cell_color]],
                                         colors = cell_color_code,
                                         legendgroup = annotated_DT_centers[[cell_color]],
                                         marker = list(size = center_point_size,symbol = "x",symbols = "x"),
                                         showlegend = F)
      }

      if(show_center_label == TRUE){
        dpl <- dpl %>%  plotly::add_text(x = annotated_DT_centers[["center_1"]],
                                         y = annotated_DT_centers[["center_2"]],
                                         type = 'scatter',mode = 'text',
                                         text = annotated_DT_centers[[cell_color]],
                                         textposition = 'middle right',
                                         textfont = list(color = '#000000', size = label_size),showlegend = F)
      }

    }
    if(show_other_cells == TRUE){
      dpl <- dpl %>% plotly::add_trace(type = "scatter",mode = "markers",
                                       x = annotated_DT_other[[dim_names[1]]],
                                       y = annotated_DT_other[[dim_names[2]]],
                                       marker = list(size = other_point_size,color = other_cell_color),
                                       showlegend = FALSE)
    }
    if(dim_reduction_to_use == 'pca') {
      x_name = paste0('pca','-',dim_names[1])
      y_name = paste0('pca','-',dim_names[2])
      x_title = sprintf('%s explains %.02f%% of variance', x_name, var_expl_vec[1])
      y_title = sprintf('%s explains %.02f%% of variance', y_name, var_expl_vec[2])
    }
    else{
      x_title = paste(dim_reduction_to_use,dim_names[1],sep = " ")
      y_title = paste(dim_reduction_to_use,dim_names[2],sep = " ")
    }
    dpl <- dpl %>% plotly::layout(xaxis = list(title = x_title),
                                  yaxis = list(title = y_title),
                                  legend = list(x = 100, y = 0.5,font = list(family = "sans-serif",size = legend_text_size)))
  }
  #3D plot
  else if(!is.null(dim3_to_use)){
    dpl <- plotly::plot_ly(scene = "scene1")
    if(is.null(cell_color)){
      #cell_color = "lightblue"
      dpl <- dpl %>% plotly::add_trace(type = 'scatter3d',mode = "markers",
                                       x = annotated_DT[[dim_names[1]]],
                                       y = annotated_DT[[dim_names[2]]],
                                       z = annotated_DT[[dim_names[3]]],
                                       color = "lightblue",
                                       colors = "lightblue",
                                       marker = list(size = dim_point_size),
                                       showlegend = F)
      #legendgroup = annotated_DT[[cell_color]])
    }
    else{
      if(cell_color %in% colnames(annotated_DT)){
        if(is.null(cell_color_code)) {
          number_colors=length(unique(annotated_DT[[cell_color]]))
          cell_color_code = Giotto:::getDistinctColors(n = number_colors)
        }
        if(color_as_factor){
          annotated_DT[[cell_color]] <- as.factor(annotated_DT[[cell_color]])
        }
        dpl <- dpl %>% plotly::add_trace(type = 'scatter3d',mode = "markers",
                                         x = annotated_DT[[dim_names[1]]],
                                         y = annotated_DT[[dim_names[2]]],
                                         z = annotated_DT[[dim_names[3]]],
                                         color = annotated_DT[[cell_color]],
                                         colors = cell_color_code,
                                         marker = list(size = dim_point_size),
                                         legendgroup = annotated_DT[[cell_color]])
      }

      else{
        stop("cell_color does not exist!\n")
      }
    }
    if(show_other_cells == TRUE){
      dpl <- dpl %>% plotly::add_trace(type = "scatter3d",mode = "markers",
                                       x = annotated_DT_other[[dim_names[1]]],
                                       y = annotated_DT_other[[dim_names[2]]],
                                       z = annotated_DT_other[[dim_names[3]]],
                                       marker = list(size = other_point_size,color = other_cell_color),
                                       showlegend = FALSE)
    }

    if(show_NN_network){
      edges <- plotly_network(annotated_network_DT,
                              "from_Dim.1","from_Dim.2","from_Dim.3",
                              "to_Dim.1","to_Dim.2","to_Dim.3")
      if(is.null(nn_network_alpha)){
        nn_network_alpha = 0.5
      }
      else if(is.character(nn_network_alpha)){
        warning("Edge_alpha for plotly mode is not adjustable yet. Default 0.5 will be set\n")
        nn_network_alpha = 0.5
      }

      dpl <- dpl %>% plotly::add_trace(name = network_name,
                                       mode = "lines",
                                       type = "scatter3d",
                                       data = edges,
                                       x = ~x,y=~y,z=~z,
                                       line=list(color=network_color),
                                       opacity=nn_network_alpha)
    }
    if((show_cluster_center == TRUE | show_center_label == TRUE)& !is.null(cell_color)){
      annotated_DT_centers = annotated_DT[, .(center_1 = median(get(dim_names[1])),
                                              center_2 = median(get(dim_names[2])),
                                              center_3 = median(get(dim_names[3]))),
                                          by = cell_color]
      annotated_DT_centers[[cell_color]] <- as.factor(annotated_DT_centers[[cell_color]])
      if(show_cluster_center == TRUE){
        dpl <- dpl %>% plotly::add_trace(mode = "markers",
                                         type = "scatter3d",
                                         data = annotated_DT_centers,
                                         x = ~center_1,
                                         y = ~center_2,
                                         z = ~center_3,
                                         color = annotated_DT_centers[[cell_color]],
                                         colors = cell_color_code,
                                         marker=list(size = 2,symbol = "x",symbols = "x"),
                                         legendgroup = annotated_DT_centers[[cell_color]],
                                         showlegend = F)
      }
      if(show_center_label == TRUE){
        cat(" center label is not clear to see in 3D plot\n You can shut it down with show_center_label = F\n")
        dpl <- dpl %>% plotly::add_trace(mode = "text",
                                         type = "scatter3d",
                                         data = annotated_DT_centers,
                                         x = ~center_1,
                                         y = ~center_2,
                                         z = ~center_3,
                                         text = annotated_DT_centers[[cell_color]],
                                         legendgroup = annotated_DT_centers[[cell_color]],
                                         showlegend = F)
      }

    }
    if(dim_reduction_to_use == 'pca') {
      x_name = paste0('pca','-',dim_names[1])
      y_name = paste0('pca','-',dim_names[2])
      z_name = paste0('pca','-',dim_names[3])
      x_title = sprintf('%s explains %.02f%% of variance', x_name, var_expl_vec[1])
      y_title = sprintf('%s explains %.02f%% of variance', y_name, var_expl_vec[2])
      z_title = sprintf('%s explains %.02f%% of variance', z_name, var_expl_vec[3])
    }
    else{
      x_title = paste(dim_reduction_to_use,dim_names[1],sep = " ")
      y_title = paste(dim_reduction_to_use,dim_names[2],sep = " ")
      z_title = paste(dim_reduction_to_use,dim_names[3],sep = " ")
    }

  }



  #********** spatial plot ***********#
  if(is.null(sdimx) | is.null(sdimy)) {
    # cat('first and second dimenion need to be defined, default is first 2 \n')
    sdimx = 'sdimx'
    sdimy = 'sdimy'
  }

  # 2D plot
  if(is.null(sdimz)){
    spl <- plotly::plot_ly()

    if(show_spatial_network == TRUE) {
      if(is.null(spatial_network)){
        stop("No usable spatial network specified! Please choose a network with spatial_network_name=xxx")
      }
      else{
        if(is.null(spatial_network_alpha)) {
          spatial_network_alpha = 0.5
        }
        else if(is.character(spatial_network_alpha)){
          warning("Edge_alpha for plotly mode is not adjustable yet. Default 0.5 will be set\n")
          spatial_network_alpha = 0.5
        }
        spl <- spl %>% plotly::add_segments(name = spatial_network_name,
                                            type = "scatter",
                                            x = spatial_network[["sdimx_begin"]],
                                            y = spatial_network[["sdimy_begin"]],
                                            xend = spatial_network[["sdimx_end"]],
                                            yend = spatial_network[["sdimy_end"]],
                                            line = list(color = network_color,
                                                        width = 0.5),
                                            opacity=spatial_network_alpha)
      }
    }
    if(show_spatial_grid == TRUE){
      if(is.null(spatial_grid)){
        stop("No usable spatial grid specified! Please choose a network with spatial_grid_name=xxx")
      }
      else{
        if(is.null(spatial_grid_color)) {
          spatial_grid_color = 'black'
        }
        edges <- plotly_grid(spatial_grid)
        spl <- spl %>% plotly::add_segments(name = "spatial_grid",
                                            type = "scatter",
                                            data = edges,
                                            x = ~x,
                                            y = ~y,
                                            xend = ~x_end,
                                            yend = ~y_end,
                                            line = list(color = spatial_grid_color,
                                                        width = 1),
                                            opacity=spatial_grid_alpha)

      }
    }
    if(is.null(cell_color)){
      #cell_color = "lightblue"
      spl <- spl %>% plotly::add_trace(type = "scatter",mode = "markers",
                                       x = annotated_DT[[sdimx]],
                                       y = annotated_DT[[sdimy]],
                                       #color = "lightblue",
                                       #colors = "lightblue",
                                       marker = list(size = spatial_point_size,
                                                     color = "lightblue"),
                                       showlegend = F)
    }

    else if(cell_color %in% colnames(annotated_DT)){
      if(color_as_factor){
        annotated_DT[[cell_color]] <- as.factor(annotated_DT[[cell_color]])
      }


      spl <- spl %>% plotly::add_trace(type = "scatter",mode = "markers",
                                       x = annotated_DT[[sdimx]],
                                       y = annotated_DT[[sdimy]],
                                       color = annotated_DT[[cell_color]],
                                       colors = cell_color_code,
                                       legendgroup = annotated_DT[[cell_color]],
                                       marker = list(size = spatial_point_size),
                                       showlegend = F)
    }
    else{
      stop("cell_color doesn't exist!\n")
    }
    if(show_other_cells == TRUE){
      spl <- spl %>% plotly::add_trace(type = "scatter",mode = "markers",
                                       x = annotated_DT_other[[sdimx]],
                                       y = annotated_DT_other[[sdimy]],
                                       marker = list(size = other_point_size,color = other_cell_color),
                                       showlegend = FALSE)
    }
    spl <- spl %>% plotly::layout(xaxis = list(title = "X"),
                                  yaxis = list(title = "Y"),
                                  legend = list(x = 100, y = 0.5,font = list(family = "sans-serif",size = legend_text_size)))

  }


  # 3D plot
  else{
    axis_scale = match.arg(axis_scale, c("cube","real","custom"))

    ratio = plotly_axis_scale_3D(annotated_DT,sdimx = sdimx,sdimy = sdimy,sdimz = sdimz,
                                 mode = axis_scale,custom_ratio = custom_ratio)
    spl <- plotly::plot_ly(scene = "scene2")
    if(!is.null(cell_color)) {
      if(cell_color %in% colnames(annotated_DT)){
        annotated_DT[[cell_color]] <- as.factor(annotated_DT[[cell_color]])
        spl <- spl %>% plotly::add_trace(type = 'scatter3d',mode = 'markers',
                                         x = annotated_DT[[sdimx]],
                                         y = annotated_DT[[sdimy]],
                                         z = annotated_DT[[sdimz]],
                                         color = annotated_DT[[cell_color]],
                                         colors = cell_color_code,
                                         legendgroup = annotated_DT[[cell_color]],
                                         marker = list(size = spatial_point_size),
                                         showlegend = F)
      }
      else{
        stop("cell_color doesn't exist!\n")
      }
    }
    else{
      spl <- spl %>% plotly::add_trace(type = 'scatter3d',mode = 'markers',
                                       x = annotated_DT$sdimx,
                                       y = annotated_DT$sdimy,
                                       z = annotated_DT$sdimz,
                                       color = "lightblue",
                                       colors = "lightblue",
                                       #legendgroup = annotated_DT[[cell_color]],
                                       marker = list(size = spatial_point_size),
                                       showlegend = F)
    }
    if(show_other_cells == TRUE){
      spl <- spl %>% plotly::add_trace(type = "scatter3d",mode = "markers",
                                       x = annotated_DT_other[[sdimx]],
                                       y = annotated_DT_other[[sdimy]],
                                       z = annotated_DT_other[[sdimz]],
                                       marker = list(size = other_point_size,color = other_cell_color),
                                       showlegend = FALSE)
    }
    if(show_spatial_network == TRUE) {
      if(is.null(spatial_network)){
        stop("No usable spatial network specified! Please choose a network with spatial_network_name=xxx")
      }
      else{
        if(is.null(spatial_network_alpha)) {
          spatial_network_alpha = 0.5
        }
        else if(is.character(spatial_network_alpha)){
          warning("Edge_alpha for plotly mode is not adjustable yet. Default 0.5 will be set\n")
          spatial_network_alpha = 0.5
        }
        edges <- plotly_network(spatial_network)

        spl <- spl %>% plotly::add_trace(name = "sptial network",
                                         mode = "lines",
                                         type = "scatter3d",
                                         data = edges,
                                         x = ~x,y=~y,z=~z,
                                         line=list(color=network_color),
                                         opacity=spatial_network_alpha)
      }
    }

    if(show_spatial_grid == TRUE){
      cat("3D grid is not clear to view\n")
    }

  }




  if(is.null(dim3_to_use) & is.null(sdimz)){
    if(plot_alignment == 'vertical'){
      combo_plot <- plotly::subplot(dpl,spl,nrows = 2,titleX = TRUE,titleY = TRUE)
    }
    else{
      combo_plot <- plotly::subplot(dpl,spl,titleX = TRUE,titleY = TRUE)
    }
  }

  else if(!is.null(dim3_to_use) & is.null(sdimz)){
    if(plot_alignment == 'vertical'){
      combo_plot <- plotly::subplot(dpl,spl,nrows = 2,titleX = TRUE,titleY = TRUE)%>%
        plotly::layout(scene = list(domain = list(x = c(0, 1), y = c(0,0.5)),
                                    xaxis = list(title = x_title),
                                    yaxis = list(title = y_title),
                                    zaxis = list(title = z_title)))
    }
    else{
      combo_plot <- plotly::subplot(dpl,spl,titleX = TRUE,titleY = TRUE) %>%
        plotly::layout(scene = list(domain = list(x = c(0, 0.5), y = c(0,1)),
                                    xaxis = list(title = x_title),
                                    yaxis = list(title = y_title),
                                    zaxis = list(title = z_title)))
    }
  }

  else if(is.null(dim3_to_use) & !is.null(sdimz)){
    if(plot_alignment == 'vertical'){
      combo_plot <- plotly::subplot(dpl,spl,nrows = 2,titleX = TRUE,titleY = TRUE)%>%
        plotly::layout(scene2 = list(domain = list(x = c(0, 1), y = c(0.5,1)),
                                     xaxis = list(title = "X",nticks = x_ticks),
                                     yaxis = list(title = "Y",nticks = y_ticks),
                                     zaxis = list(title = "Z",nticks = z_ticks),
                                     aspectmode='manual',
                                     aspectratio = list(x=ratio[[1]],
                                                        y=ratio[[2]],
                                                        z=ratio[[3]])))
    }
    else{
      combo_plot <- plotly::subplot(dpl,spl,titleX = TRUE,titleY = TRUE) %>%
        plotly::layout(scene2 = list(domain = list(x = c(0.5, 1), y = c(0,1)),
                                     xaxis = list(title = "X",nticks = x_ticks),
                                     yaxis = list(title = "Y",nticks = y_ticks),
                                     zaxis = list(title = "Z",nticks = z_ticks),
                                     aspectmode='manual',
                                     aspectratio = list(x=ratio[[1]],
                                                        y=ratio[[2]],
                                                        z=ratio[[3]])))
    }
  }

  else if(!is.null(dim3_to_use) & !is.null(sdimz)){
    if(plot_alignment == 'vertical'){
      combo_plot <- plotly::subplot(dpl,spl,nrows = 2,titleX = TRUE,titleY = TRUE)%>%
        plotly::layout(scene = list(domain = list(x = c(0, 1), y = c(0,0.5)),
                                    xaxis = list(title = x_title),
                                    yaxis = list(title = y_title),
                                    zaxis = list(title = z_title)),
                       scene2 = list(domain = list(x = c(0, 1), y = c(0.5,1)),
                                     xaxis = list(title = "X",nticks = x_ticks),
                                     yaxis = list(title = "Y",nticks = y_ticks),
                                     zaxis = list(title = "Z",nticks = z_ticks),
                                     aspectmode='manual',
                                     aspectratio = list(x=ratio[[1]],
                                                        y=ratio[[2]],
                                                        z=ratio[[3]])))
    }
    else{
      combo_plot <- plotly::subplot(dpl,spl,titleX = TRUE,titleY = TRUE) %>%
        plotly::layout(scene = list(domain = list(x = c(0, 0.5), y = c(0,1)),
                                    xaxis = list(title = x_title),
                                    yaxis = list(title = y_title),
                                    zaxis = list(title = z_title)),
                       scene2 = list(domain = list(x = c(0.5, 1), y = c(0,1)),
                                     xaxis = list(title = "X",nticks = x_ticks),
                                     yaxis = list(title = "Y",nticks = y_ticks),
                                     zaxis = list(title = "Z",nticks = z_ticks),
                                     aspectmode='manual',
                                     aspectratio = list(x=ratio[[1]],
                                                        y=ratio[[2]],
                                                        z=ratio[[3]])))
    }
  }

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




#' @title spatGenePlot3D
#' @name spatGenePlot3D
#' @description Visualize cells and gene expression according to spatial coordinates
#' @param gobject giotto object
#' @param expression_values gene expression values to use
#' @param genes genes to show
#' @param genes_high_color color represents high gene expression
#' @param genes_mid_color color represents middle gene expression
#' @param genes_low_color color represents low gene expression
#' @param show_network show underlying spatial network
#' @param network_color color of spatial network
#' @param spatial_network_name name of spatial network to use
#' @param show_grid show spatial grid
#' @param grid_color color of spatial grid
#' @param spatial_grid_name name of spatial grid to use
#' @param midpoint expression midpoint
#' @param scale_alpha_with_expression scale expression with ggplot alpha parameter
#' @param point_size size of point (cell)
#' @param show_legend show legend
#' @param show_plot show plots
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters from \code{\link{all_plots_save_function}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @param ... parameters for cowplot::save_plot()
#' @return ggplot
#' @details Description of parameters.
#' @export
#' @examples
#'     spatGenePlot3D(gobject)
spatGenePlot3D <- function(gobject,
                           expression_values = c('normalized', 'scaled', 'custom'),
                           genes,
                           show_network = F,
                           network_color = NULL,
                           spatial_network_name = 'Delaunay_network',
                           edge_alpha = NULL,
                           show_grid = F,

                           cluster_column = NULL,
                           select_cell_groups = NULL,
                           select_cells = NULL,
                           show_other_cells = T,
                           other_cell_color = 'lightgrey',
                           other_point_size = 1,

                           genes_high_color = NULL,
                           genes_mid_color = "white",
                           genes_low_color = "blue",
                           spatial_grid_name = 'spatial_grid',
                           point_size = 2,
                           show_legend = T,
                           axis_scale = c("cube","real","custom"),
                           custom_ratio = NULL,
                           x_ticks = NULL,
                           y_ticks = NULL,
                           z_ticks = NULL,
                           show_plot = NA,
                           return_plot = NA,
                           save_plot = NA,
                           save_param =  list(),
                           default_save_name = "spatGenePlot3D"){
  selected_genes = genes

  values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))
  expr_values = Giotto:::select_expression_values(gobject = gobject, values = values)

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
    spatial_network = select_spatialNetwork(gobject,name = spatial_network_name,return_network_Obj = FALSE)
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

  if(!is.null(select_cells) & !is.null(select_cell_groups)) {
    cat('You have selected both individual cell IDs and a group of cells \n')
    group_cell_IDs = cell_locations_metadata[get(cluster_column) %in% select_cell_groups][['cell_ID']]
    select_cells = unique(c(select_cells, group_cell_IDs))
  } else if(!is.null(select_cell_groups)) {
    select_cells = cell_locations_metadata[get(cluster_column) %in% select_cell_groups][['cell_ID']]
  }

  if(!is.null(select_cells)) {
    cell_locations_metadata_other = cell_locations_metadata[!cell_locations_metadata$cell_ID %in% select_cells]
    cell_locations_metadata_selected = cell_locations_metadata[cell_locations_metadata$cell_ID %in% select_cells]
    spatial_network <- spatial_network[spatial_network$to %in% select_cells & spatial_network$from %in% select_cells]

    # if specific cells are selected
    cell_locations_metadata = cell_locations_metadata_selected
  }

  cell_locations_metadata_genes <- merge(cell_locations_metadata, t_sub_expr_data_DT, by = 'cell_ID')



  ## plotting ##
  axis_scale = match.arg(axis_scale, c("cube","real","custom"))

  ratio = plotly_axis_scale_3D(cell_locations_metadata_genes,sdimx = "sdimx",sdimy = "sdimy",sdimz = "sdimz",
                               mode = axis_scale,custom_ratio = custom_ratio)


  ## spatial network data
  if(!is.null(spatial_network) & show_network == TRUE){
    edges <- plotly_network(spatial_network)
  }
  ##Point layer
  if(length(selected_genes) > 4){
    stop("\n The max number of genes showed together is 4.Otherwise it will be too small to see\n
              \n If you have more genes to show, please divide them into groups\n")
  }
  savelist <- list()
  for(i in 1:length(selected_genes)){
    gene = selected_genes[i]
    if(!is.null(genes_high_color)){
      if(length(genes_high_color)!=length(selected_genes) & length(genes_high_color)!=1){
        stop('\n The number of genes and their corresbonding do not match\n')
      }
      else if(length(genes_high_color) == 1){
        genes_high_color = rep(genes_high_color,length(selected_genes))
      }
    }
    else{
      genes_high_color = rep("red",length(selected_genes))
    }
    pl <- plotly::plot_ly(name = gene,

                          scene=paste("scene",i,sep = "")) %>%
      plotly::add_trace(data = cell_locations_metadata_genes,
                        type = 'scatter3d',mode = "markers",
                        x = ~sdimx, y = ~sdimy, z = ~sdimz,
                        marker = list(size = point_size),
                        color = cell_locations_metadata_genes[[gene]],
                        colors = c(genes_low_color,genes_mid_color,genes_high_color[i]))

    if(show_other_cells == T){
      pl <- pl %>% plotly::add_trace(name = "unselected cells",
                                     data = cell_locations_metadata_other,
                                     type = 'scatter3d',mode = "markers",
                                     x = ~sdimx, y = ~sdimy, z = ~sdimz,
                                     marker = list(size = other_point_size,color = other_cell_color))
    }

    ## plot spatial network
    if(show_network == TRUE) {
      if(is.null(network_color)) {
        network_color = 'lightblue'
      }
      if(is.null(edge_alpha)) {
        edge_alpha = 0.5
      }
      else if (is.character(edge_alpha)){
        edge_alpha = 0.5
        cat("\nEdge_alpha for plotly mode is not adjustable yet. Default 0.5 will be set\n")
      }
      pl <- pl %>% plotly::add_trace(name = "sptial network",
                                     mode = "lines",
                                     type = "scatter3d",
                                     data = edges,
                                     x = ~x,y=~y,z=~z,
                                     line=list(color=network_color,width = 0.5),
                                     opacity = edge_alpha,
                                     showlegend = F)
    }
    ##plot spatial grid
    if(!is.null(spatial_grid) & show_grid == TRUE){
      cat("\n spatial grid is not clear in 3D plot \n")
    }
    pl <- pl %>% plotly::colorbar(title = gene)
    savelist[[gene]] <- pl
  }


  if(length(savelist) == 1){
    cowplot <- savelist[[1]] %>% plotly::layout(scene = list(xaxis = list(title = "X",nticks = x_ticks),
                                                             yaxis = list(title = "Y",nticks = y_ticks),
                                                             zaxis = list(title = "Z",nticks = z_ticks),
                                                             aspectmode='manual',
                                                             aspectratio = list(x=ratio[[1]],
                                                                                y=ratio[[2]],
                                                                                z=ratio[[3]])))

  }
  else if(length(savelist)==2){
    cowplot <- suppressWarnings(plotly::subplot(savelist)%>% plotly::layout(scene = list(xaxis = list(title = "X",nticks = x_ticks),
                                                                                         yaxis = list(title = "Y",nticks = y_ticks),
                                                                                         zaxis = list(title = "Z",nticks = z_ticks),
                                                                                         aspectmode='manual',
                                                                                         aspectratio = list(x=ratio[[1]],
                                                                                                            y=ratio[[2]],
                                                                                                            z=ratio[[3]])),
                                                                            scene2 = list(xaxis = list(title = "X",nticks = x_ticks),
                                                                                          yaxis = list(title = "Y",nticks = y_ticks),
                                                                                          zaxis = list(title = "Z",nticks = z_ticks),
                                                                                          aspectmode='manual',
                                                                                          aspectratio = list(x=ratio[[1]],
                                                                                                             y=ratio[[2]],
                                                                                                             z=ratio[[3]])),
                                                                            #annotations = annotations,
                                                                            legend = list(x = 100, y = 0)))
  }
  else if(length(savelist)==3){
    cowplot <- suppressWarnings(plotly::subplot(savelist)%>% plotly::layout(scene = list(xaxis = list(title = "X",nticks = x_ticks),
                                                                                         yaxis = list(title = "Y",nticks = y_ticks),
                                                                                         zaxis = list(title = "Z",nticks = z_ticks),
                                                                                         aspectmode='manual',
                                                                                         aspectratio = list(x=ratio[[1]],
                                                                                                            y=ratio[[2]],
                                                                                                            z=ratio[[3]])),
                                                                            scene2 = list(xaxis = list(title = "X",nticks = x_ticks),
                                                                                          yaxis = list(title = "Y",nticks = y_ticks),
                                                                                          zaxis = list(title = "Z",nticks = z_ticks),
                                                                                          aspectmode='manual',
                                                                                          aspectratio = list(x=ratio[[1]],
                                                                                                             y=ratio[[2]],
                                                                                                             z=ratio[[3]])),
                                                                            scene3 = list(xaxis = list(title = "X",nticks = x_ticks),
                                                                                          yaxis = list(title = "Y",nticks = y_ticks),
                                                                                          zaxis = list(title = "Z",nticks = z_ticks),
                                                                                          aspectmode='manual',
                                                                                          aspectratio = list(x=ratio[[1]],
                                                                                                             y=ratio[[2]],
                                                                                                             z=ratio[[3]])),
                                                                            legend = list(x = 100, y = 0)))
  }
  else if(length(savelist)==4){


    cowplot <- suppressWarnings(plotly::subplot(savelist)%>% plotly::layout(scene = list(xaxis = list(title = "X",nticks = x_ticks),
                                                                                         yaxis = list(title = "Y",nticks = y_ticks),
                                                                                         zaxis = list(title = "Z",nticks = z_ticks),
                                                                                         aspectmode='manual',
                                                                                         aspectratio = list(x=ratio[[1]],
                                                                                                            y=ratio[[2]],
                                                                                                            z=ratio[[3]])),
                                                                            scene2 = list(xaxis = list(title = "X",nticks = x_ticks),
                                                                                          yaxis = list(title = "Y",nticks = y_ticks),
                                                                                          zaxis = list(title = "Z",nticks = z_ticks),
                                                                                          aspectmode='manual',
                                                                                          aspectratio = list(x=ratio[[1]],
                                                                                                             y=ratio[[2]],
                                                                                                             z=ratio[[3]])),
                                                                            scene3 = list(xaxis = list(title = "X",nticks = x_ticks),
                                                                                          yaxis = list(title = "Y",nticks = y_ticks),
                                                                                          zaxis = list(title = "Z",nticks = z_ticks),
                                                                                          aspectmode='manual',
                                                                                          aspectratio = list(x=ratio[[1]],
                                                                                                             y=ratio[[2]],
                                                                                                             z=ratio[[3]])),
                                                                            scene4 = list(xaxis = list(title = "X",nticks = x_ticks),
                                                                                          yaxis = list(title = "Y",nticks = y_ticks),
                                                                                          zaxis = list(title = "Z",nticks = z_ticks),
                                                                                          aspectmode='manual',
                                                                                          aspectratio = list(x=ratio[[1]],
                                                                                                             y=ratio[[2]],
                                                                                                             z=ratio[[3]])),
                                                                            legend = list(x = 100, y = 0)))

  }


  show_plot = ifelse(is.na(show_plot), readGiottoInstructions(gobject, param = 'show_plot'), show_plot)
  save_plot = ifelse(is.na(save_plot), readGiottoInstructions(gobject, param = 'save_plot'), save_plot)
  return_plot = ifelse(is.na(return_plot), readGiottoInstructions(gobject, param = 'return_plot'), return_plot)
  ## print plot
  if(show_plot == TRUE) {
    print(cowplot)
  }

  ## save plot
  if(save_plot == TRUE) {
    do.call('all_plots_save_function', c(list(gobject = gobject, plot_object = cowplot, default_save_name = default_save_name), save_param))
  }

  ## return plot
  if(return_plot == TRUE) {
    return(cowplot)
  }



}



#' @title dimGenePlot3D
#' @name dimGenePlot3D
#' @description Visualize cells and gene expression according to dimension reduction coordinates
#' @param gobject giotto object
#' @param expression_values gene expression values to use
#' @param genes genes to show
#' @param dim_reduction_to_use dimension reduction to use
#' @param dim_reduction_name dimension reduction name
#' @param dim1_to_use dimension to use on x-axis
#' @param dim2_to_use dimension to use on y-axis
#' @param dim3_to_use dimension to use on z-axis
#' @param show_NN_network show underlying NN network
#' @param nn_network_to_use type of NN network to use (kNN vs sNN)
#' @param network_name name of NN network to use, if show_NN_network = TRUE
#' @param edge_alpha column to use for alpha of the edges
#' @param point_size size of point (cell)
#' @param show_legend show legend
#' @param show_plot show plots
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters from \code{\link{all_plots_save_function}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @param ... parameters for cowplot::save_plot()
#' @return ggplot
#' @details Description of parameters.
#' @export
#' @examples
#'     dimGenePlot3D(gobject)
dimGenePlot3D <- function(gobject,
                          expression_values = c('normalized', 'scaled', 'custom'),
                          genes = NULL,
                          dim_reduction_to_use = 'umap',
                          dim_reduction_name = 'umap',
                          dim1_to_use = 1,
                          dim2_to_use = 2,
                          dim3_to_use = 3,
                          show_NN_network = F,
                          nn_network_to_use = 'sNN',
                          network_name = 'sNN.pca',
                          network_color = "lightgray",

                          cluster_column = NULL,
                          select_cell_groups = NULL,
                          select_cells = NULL,
                          show_other_cells = T,
                          other_cell_color = 'lightgrey',
                          other_point_size = 1,

                          edge_alpha = NULL,
                          point_size = 2,
                          genes_high_color = NULL,
                          genes_mid_color = "white",
                          genes_low_color = "blue",
                          show_legend = T,
                          show_plot = NA,
                          return_plot = NA,
                          save_plot = NA,
                          save_param =  list(),
                          default_save_name = "dimGenePlot3D"){
  ## select genes ##
  selected_genes = genes
  values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))
  expr_values = Giotto:::select_expression_values(gobject = gobject, values = values)

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
  dim_dfr = gobject@dimension_reduction$cells[[dim_reduction_to_use]][[dim_reduction_name]]$coordinates[,c(dim1_to_use, dim2_to_use,dim3_to_use)]
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


  ## create subsets if needed
  if(!is.null(select_cells) & !is.null(select_cell_groups)) {
    if(is.null(cluster_column)) {
      stop('\n selection of cells is based on cell_color paramter, which is a metadata column \n')
    }
    cat('You have selected both individual cell IDs and a group of cells \n')
    group_cell_IDs = annotated_DT[get(cluster_column) %in% select_cell_groups][['cell_ID']]
    select_cells = unique(c(select_cells, group_cell_IDs))
  } else if(!is.null(select_cell_groups)) {
    select_cells = annotated_DT[get(cluster_column) %in% select_cell_groups][['cell_ID']]
  }

  if(!is.null(select_cells)) {
    annotated_DT_other = annotated_DT[!annotated_DT$cell_ID %in% select_cells]
    annotated_DT_selected = annotated_DT[annotated_DT$cell_ID %in% select_cells]

    if(show_NN_network == TRUE) {
      annotated_network_DT <- annotated_network_DT[annotated_network_DT$to %in% select_cells & annotated_network_DT$from %in% select_cells]
    }

    # if specific cells are selected
    annotated_DT = annotated_DT_selected
  }

  ## if no subsets are required
  if(is.null(select_cells) & is.null(select_cell_groups)) {
    annotated_DT_selected = annotated_DT
    annotated_DT_other    = NULL
  }

  ## merge gene info
  annotated_gene_DT = merge(annotated_DT, t_sub_expr_data_DT, by = 'cell_ID')


  ## visualize multipe plots ##
  ## 3D plots ##


  if(show_NN_network == TRUE){
    edges <- plotly_network(annotated_network_DT,
                            "from_Dim.1","from_Dim.2","from_Dim.3",
                            "to_Dim.1","to_Dim.2","to_Dim.3")
  }
  ##Point layer
  if(length(selected_genes) > 4){
    stop("\n The max number of genes showed together is 4.Otherwise it will be too small to see\n
              \n If you have more genes to show, please divide them into groups\n")
  }
  if(!is.null(genes_high_color)){
    if(length(genes_high_color)!=length(selected_genes)&length(genes_high_color) != 1){
      stop('\n The number of genes and their corresbonding do not match\n')
    }
  }
  else if (is.null(genes_high_color)){
    genes_high_color = rep("red",length(selected_genes))
  }
  else{
    genes_high_color = rep(genes_high_color,length(selected_genes))
  }

  titleX = title = paste(dim_reduction_to_use,dim_names[1],sep = " ")
  titleY = title = paste(dim_reduction_to_use,dim_names[2],sep = " ")
  titleZ = title = paste(dim_reduction_to_use,dim_names[3],sep = " ")
  savelist <- list()
  for(i in 1:length(selected_genes)){

    gene = selected_genes[i]

    pl <- plotly::plot_ly(name = gene,scene=paste("scene",i,sep = ""))
    pl <- pl %>%  plotly::add_trace(data = annotated_gene_DT,type = 'scatter3d',mode = "markers",
                                    x = annotated_gene_DT[[dim_names[1]]],
                                    y = annotated_gene_DT[[dim_names[2]]],
                                    z = annotated_gene_DT[[dim_names[3]]],
                                    color = annotated_gene_DT[[gene]],
                                    colors = c(genes_low_color,genes_mid_color,genes_high_color[i]),
                                    marker = list(size = point_size))
    if(show_other_cells == T){
      pl <- pl %>% plotly::add_trace(name = "unselected cells",
                                     data = annotated_DT_other,
                                     type = 'scatter3d',mode = "markers",
                                     x = annotated_DT_other[[dim_names[1]]],
                                     y = annotated_DT_other[[dim_names[2]]],
                                     z = annotated_DT_other[[dim_names[3]]],
                                     marker = list(size = other_point_size,color = other_cell_color))
    }

    ## plot spatial network
    if(show_NN_network == TRUE) {
      pl <- pl %>% plotly::add_trace(name = "sptial network",mode = "lines",
                                     type = "scatter3d",opacity = edge_alpha,
                                     showlegend = F,
                                     data = edges,
                                     x = ~x,y=~y,z=~z,
                                     line=list(color=network_color, width = 0.5))
    }
    pl <- pl %>% plotly::colorbar(title = gene)
    savelist[[gene]] <- pl
  }

  if(length(savelist) == 1){
    cowplot <- savelist[[1]] %>% plotly::layout(scene = list(
      xaxis = list(title = titleX),
      yaxis = list(title = titleY),
      zaxis = list(title = titleZ)))
  }
  else if(length(savelist)==2){
    cowplot <- suppressWarnings(plotly::subplot(savelist,titleX = TRUE,titleY = TRUE)%>%
                                  plotly::layout(scene = list(domain = list(x = c(0, 0.5), y = c(0,1)),
                                                              xaxis = list(title = titleX),
                                                              yaxis = list(title = titleY),
                                                              zaxis = list(title = titleZ)),
                                                 scene2 = list(domain = list(x = c(0.5, 1), y = c(0,1)),
                                                               xaxis = list(title = titleX),
                                                               yaxis = list(title = titleY),
                                                               zaxis = list(title = titleZ)),
                                                 legend = list(x = 100, y = 0)))
  }
  else if(length(savelist)==3){
    cowplot <- suppressWarnings(plotly::subplot(savelist,titleX = TRUE,titleY = TRUE)%>%
                                  plotly::layout(scene = list(domain = list(x = c(0, 0.5), y = c(0,0.5)),
                                                              xaxis = list(title = titleX),
                                                              yaxis = list(title = titleY),
                                                              zaxis = list(title = titleZ)),
                                                 scene2 = list(domain = list(x = c(0.5, 1), y = c(0,0.5)),
                                                               xaxis = list(title = titleX),
                                                               yaxis = list(title = titleY),
                                                               zaxis = list(title = titleZ)),
                                                 scene3 = list(domain = list(x = c(0, 0.5), y = c(0.5,1)),
                                                               xaxis = list(title = titleX),
                                                               yaxis = list(title = titleY),
                                                               zaxis = list(title = titleZ)),
                                                 legend = list(x = 100, y = 0)))
  }
  else if(length(savelist)==4){

    cowplot <- suppressWarnings(plotly::subplot(savelist)%>% plotly::layout(scene = list(domain = list(x = c(0, 0.5), y = c(0,0.5)),
                                                                                         xaxis = list(title = titleX),
                                                                                         yaxis = list(title = titleY),
                                                                                         zaxis = list(title = titleZ)),
                                                                            scene2 = list(domain = list(x = c(0.5, 1), y = c(0,0.5)),
                                                                                          xaxis = list(title = titleX),
                                                                                          yaxis = list(title = titleY),
                                                                                          zaxis = list(title = titleZ)),
                                                                            scene3 = list(domain = list(x = c(0, 0.5), y = c(0.5,1)),
                                                                                          xaxis = list(title = titleX),
                                                                                          yaxis = list(title = titleY),
                                                                                          zaxis = list(title = titleZ)),
                                                                            scene4 = list(domain = list(x = c(0.5, 1), y = c(0.5,1)),
                                                                                          xaxis = list(title = titleX),
                                                                                          yaxis = list(title = titleY),
                                                                                          zaxis = list(title = titleZ)),
                                                                            legend = list(x = 100, y = 0)))
  }
  show_plot = ifelse(is.na(show_plot), readGiottoInstructions(gobject, param = 'show_plot'), show_plot)
  save_plot = ifelse(is.na(save_plot), readGiottoInstructions(gobject, param = 'save_plot'), save_plot)
  return_plot = ifelse(is.na(return_plot), readGiottoInstructions(gobject, param = 'return_plot'), return_plot)
  ## print plot
  if(show_plot == TRUE) {
    print(cowplot)
  }

  ## save plot
  if(save_plot == TRUE) {
    do.call('all_plots_save_function', c(list(gobject = gobject, plot_object = cowplot, default_save_name = default_save_name), save_param))
  }

  ## return plot
  if(return_plot == TRUE) {
    return(cowplot)
  }
}




#' @title spatDimGenePlot3D
#' @name spatDimGenePlot3D
#' @description Visualize cells according to spatial AND dimension reduction coordinates in ggplot mode
#' @param gobject giotto object
#' @param expression_values gene expression values to use
#' @param plot_alignment direction to align plot
#' @param genes genes to show
#' @param dim_reduction_to_use dimension reduction to use
#' @param dim_reduction_name dimension reduction name
#' @param dim1_to_use dimension to use on x-axis
#' @param dim2_to_use dimension to use on y-axis
#' @param dim3_to_use dimension to use on z-axis
#' @param dim_point_size dim reduction plot: point size
#' @param show_NN_network show underlying NN network
#' @param nn_network_to_use type of NN network to use (kNN vs sNN)
#' @param network_name name of NN network to use, if show_NN_network = TRUE
#' @param edge_alpha_dim dim reduction plot: column to use for alpha of the edges
#' @param scale_alpha_with_expression scale expression with ggplot alpha parameter
#' @param spatial_network_name name of spatial network to use
#' @param spatial_grid_name name of spatial grid to use
#' @param spatial_point_size spatial plot: point size
#' @param point_size size of point (cell)
#' @param show_legend show legend
#' @param show_plot show plots
#' @param return_plot return plotly object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters from \code{\link{all_plots_save_function}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return plotly
#' @details Description of parameters.
#' @export
#' @examples
#'     spatDimGenePlot3D(gobject)

spatDimGenePlot3D <- function(gobject,
                              expression_values = c('normalized', 'scaled', 'custom'),
                              plot_alignment = c('horizontal','vertical'),
                              dim_reduction_to_use = 'umap',
                              dim_reduction_name = 'umap',
                              dim1_to_use = 1,
                              dim2_to_use = 2,
                              dim3_to_use = NULL,
                              sdimx="sdimx",
                              sdimy="sdimy",
                              sdimz="sdimz",
                              genes,

                              cluster_column = NULL,
                              select_cell_groups = NULL,
                              select_cells = NULL,
                              show_other_cells = T,
                              other_cell_color = 'lightgrey',
                              other_point_size = 1.5,

                              show_NN_network = F,
                              nn_network_to_use = 'sNN',
                              network_name = 'sNN.pca',
                              label_size = 16,
                              genes_low_color = "blue",
                              genes_mid_color = "white",
                              genes_high_color = "red",
                              dim_point_size = 3,
                              nn_network_alpha = 0.5,
                              show_spatial_network = F,
                              spatial_network_name = 'Delaunay_network',
                              network_color = "lightgray",
                              spatial_network_alpha = 0.5,
                              show_spatial_grid = F,
                              spatial_grid_name = 'spatial_grid',
                              spatial_grid_color = NULL,
                              spatial_grid_alpha = 0.5,
                              spatial_point_size = 3,
                              legend_text_size = 12,
                              axis_scale = c("cube","real","custom"),
                              custom_ratio = NULL,
                              x_ticks = NULL,
                              y_ticks = NULL,
                              z_ticks = NULL,

                              show_plot = NA,
                              return_plot = NA,
                              save_plot = NA,
                              save_param =  list(),
                              default_save_name = "spatDimGenePlot3D"){

  plot_alignment = match.arg(plot_alignment, choices = c( 'horizontal','vertical'))

  # ********data prepare********#
  ## select genes ##
  if(length(genes) > 1){
    warning("\n Now 3D mode can just accept one gene, only the first gene will be plot\n")
    genes = genes[1]
  }
  selected_genes = genes
  values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))
  expr_values = Giotto:::select_expression_values(gobject = gobject, values = values)

  # only keep genes that are in the dataset
  selected_genes = selected_genes[selected_genes %in% rownames(expr_values) ]
  subset_expr_data = expr_values[rownames(expr_values) %in% selected_genes, ]
  t_sub_expr_data_DT = data.table('selected_gene' = subset_expr_data, 'cell_ID' = colnames(expr_values))
  data.table::setnames(t_sub_expr_data_DT, 'selected_gene', selected_genes)


  ## dimension reduction ##
  dim_dfr = gobject@dimension_reduction$cells[[dim_reduction_to_use]][[dim_reduction_name]]$coordinates[,c(dim1_to_use, dim2_to_use,dim3_to_use)]
  dim_names = colnames(dim_dfr)
  dim_DT = data.table::as.data.table(dim_dfr); dim_DT[, cell_ID := rownames(dim_dfr)]

  ## annotated cell metadata
  cell_metadata = gobject@cell_metadata
  annotated_DT = merge(cell_metadata, dim_DT, by = 'cell_ID')
  annotated_DT = merge(annotated_DT, gobject@spatial_locs,by = 'cell_ID')
  annotated_DT = merge(annotated_DT, t_sub_expr_data_DT,by = 'cell_ID')


  ## nn network
  if(show_NN_network){
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


  ## extract spatial network
  if(!is.null(spatial_network_name)) {
    spatial_network = select_spatialNetwork(gobject,name = spatial_network_name,return_network_Obj = FALSE)
  } else {
    spatial_network = NULL
  }

  ## extract spatial grid
  if(!is.null(spatial_grid_name)) {
    spatial_grid = gobject@spatial_grid[[spatial_grid_name]]
  } else {
    spatial_grid = NULL
  }

  if(!is.null(select_cells) & !is.null(select_cell_groups)) {
    if(is.null(cluster_column)) {
      stop('\n selection of cells is based on cell_color paramter, which is a metadata column \n')
    }
    cat('You have selected both individual cell IDs and a group of cells \n')
    group_cell_IDs = annotated_DT[get(cluster_column) %in% select_cell_groups][['cell_ID']]
    select_cells = unique(c(select_cells, group_cell_IDs))
  } else if(!is.null(select_cell_groups)) {
    select_cells = annotated_DT[get(cluster_column) %in% select_cell_groups][['cell_ID']]
  }

  if(!is.null(select_cells)) {
    annotated_DT_other = annotated_DT[!annotated_DT$cell_ID %in% select_cells]
    annotated_DT_selected = annotated_DT[annotated_DT$cell_ID %in% select_cells]

    if(show_NN_network == TRUE) {
      annotated_network_DT <- annotated_network_DT[annotated_network_DT$to %in% select_cells & annotated_network_DT$from %in% select_cells]
    }
    if(show_spatial_network == TRUE){
      spatial_network <- spatial_network[spatial_network$to %in% select_cells & spatial_network$from %in% select_cells]
    }

    # if specific cells are selected
    annotated_DT = annotated_DT_selected
  }

  ## if no subsets are required
  if(is.null(select_cells) & is.null(select_cell_groups)) {
    annotated_DT_selected = annotated_DT
    annotated_DT_other    = NULL
  }

  #********** dim plot ***********#
  #2D plot
  if(is.null(dim3_to_use)){
    dpl <- plotly::plot_ly()
    if(show_NN_network == TRUE) {
      if(is.null(nn_network_alpha)) {
        nn_network_alpha = 0.5
      }
      else if(is.character(nn_network_alpha)){
        warning("Edge_alpha for plotly mode is not adjustable yet. Default 0.5 will be set\n")
        nn_network_alpha = 0.5
      }
      dpl <- dpl %>% plotly::add_segments(name = network_name,
                                          type = "scatter",
                                          x = annotated_network_DT[[from_dim_names[1]]],
                                          y = annotated_network_DT[[from_dim_names[2]]],
                                          xend = annotated_network_DT[[to_dim_names[1]]],
                                          yend = annotated_network_DT[[to_dim_names[2]]],
                                          line = list(color = network_color,
                                                      width = 0.5),
                                          opacity=nn_network_alpha)
    }

    dpl <- dpl %>% plotly::add_trace(type = "scatter",mode = "markers",
                                     x = annotated_DT[[dim_names[1]]],
                                     y = annotated_DT[[dim_names[2]]],
                                     color = annotated_DT[[selected_genes]],
                                     colors =c(genes_low_color,genes_mid_color,genes_high_color),
                                     marker = list(size = dim_point_size),
                                     showlegend = F)

    if(show_other_cells == TRUE){
      dpl <- dpl %>% plotly::add_trace(type = "scatter",mode = "markers",
                                       x = annotated_DT_other[[dim_names[1]]],
                                       y = annotated_DT_other[[dim_names[2]]],
                                       marker = list(size = other_point_size,color = other_cell_color),
                                       showlegend = FALSE)
    }

    x_title = paste(dim_reduction_to_use,dim_names[1],sep = " ")
    y_title = paste(dim_reduction_to_use,dim_names[2],sep = " ")

    dpl <- dpl %>% plotly::layout(xaxis = list(title = x_title),
                                  yaxis = list(title = y_title),
                                  legend = list(x = 100, y = 0.5,font = list(family = "sans-serif",size = legend_text_size)))
  }
  #3D plot
  else if(!is.null(dim3_to_use)){
    dpl <- plotly::plot_ly(scene = "scene1")

    dpl <- dpl %>% plotly::add_trace(type = 'scatter3d',mode = "markers",
                                     x = annotated_DT[[dim_names[1]]],
                                     y = annotated_DT[[dim_names[2]]],
                                     z = annotated_DT[[dim_names[3]]],
                                     color = annotated_DT[[selected_genes]],
                                     colors = c(genes_low_color,genes_mid_color,genes_high_color),
                                     marker = list(size = dim_point_size),
                                     showlegend = F)
    #legendgroup = annotated_DT[[cell_color]])
    if(show_other_cells == TRUE){
      dpl <- dpl %>% plotly::add_trace(type = "scatter3d",mode = "markers",
                                       x = annotated_DT_other[[dim_names[1]]],
                                       y = annotated_DT_other[[dim_names[2]]],
                                       z = annotated_DT_other[[dim_names[3]]],
                                       marker = list(size = other_point_size,color = other_cell_color),
                                       showlegend = FALSE)
    }

    if(show_NN_network){
      edges <- plotly_network(annotated_network_DT,
                              "from_Dim.1","from_Dim.2","from_Dim.3",
                              "to_Dim.1","to_Dim.2","to_Dim.3")
      if(is.null(nn_network_alpha)){
        nn_network_alpha = 0.5
      }
      else if(is.character(nn_network_alpha)){
        warning("Edge_alpha for plotly mode is not adjustable yet. Default 0.5 will be set\n")
        nn_network_alpha = 0.5
      }

      dpl <- dpl %>% plotly::add_trace(name = network_name,
                                       mode = "lines",
                                       type = "scatter3d",
                                       data = edges,
                                       x = ~x,y=~y,z=~z,
                                       line=list(color=network_color),
                                       opacity=nn_network_alpha)
    }


    x_title = paste(dim_reduction_to_use,dim_names[1],sep = " ")
    y_title = paste(dim_reduction_to_use,dim_names[2],sep = " ")
    z_title = paste(dim_reduction_to_use,dim_names[3],sep = " ")

  }
  dpl <- dpl %>% plotly::colorbar(title = selected_genes)


  #********** spatial plot ***********#
  if(is.null(sdimx) | is.null(sdimy)) {
    # cat('first and second dimenion need to be defined, default is first 2 \n')
    sdimx = 'sdimx'
    sdimy = 'sdimy'
  }

  # 2D plot
  if(is.null(sdimz)){
    spl <- plotly::plot_ly()

    if(show_spatial_network == TRUE) {
      if(is.null(spatial_network)){
        stop("No usable spatial network specified! Please choose a network with spatial_network_name=xxx")
      }
      else{
        if(is.null(spatial_network_alpha)) {
          spatial_network_alpha = 0.5
        }
        else if(is.character(spatial_network_alpha)){
          warning("Edge_alpha for plotly mode is not adjustable yet. Default 0.5 will be set\n")
          spatial_network_alpha = 0.5
        }
        spl <- spl %>% plotly::add_segments(name = spatial_network_name,
                                            type = "scatter",
                                            x = spatial_network[["sdimx_begin"]],
                                            y = spatial_network[["sdimy_begin"]],
                                            xend = spatial_network[["sdimx_end"]],
                                            yend = spatial_network[["sdimy_end"]],
                                            line = list(color = network_color,
                                                        width = 0.5),
                                            opacity=spatial_network_alpha)
      }
    }
    if(show_spatial_grid == TRUE){
      if(is.null(spatial_grid)){
        stop("No usable spatial grid specified! Please choose a network with spatial_grid_name=xxx")
      }
      else{
        if(is.null(spatial_grid_color)) {
          spatial_grid_color = 'black'
        }
        edges <- plotly_grid(spatial_grid)
        spl <- spl %>% plotly::add_segments(name = "spatial_grid",
                                            type = "scatter",
                                            data = edges,
                                            x = ~x,
                                            y = ~y,
                                            xend = ~x_end,
                                            yend = ~y_end,
                                            line = list(color = spatial_grid_color,
                                                        width = 1),
                                            opacity=spatial_grid_alpha)

      }
    }

    spl <- spl %>% plotly::add_trace(type = "scatter",mode = "markers",
                                     x = annotated_DT[[sdimx]],
                                     y = annotated_DT[[sdimy]],
                                     color = annotated_DT[[selected_genes]],
                                     colors = c(genes_low_color,genes_mid_color,genes_high_color),
                                     marker = list(size = spatial_point_size),
                                     showlegend = F)
    if(show_other_cells == TRUE){
      spl <- spl %>% plotly::add_trace(type = "scatter",mode = "markers",
                                       x = annotated_DT_other[[sdimx]],
                                       y = annotated_DT_other[[sdimy]],
                                       marker = list(size = other_point_size,color = other_cell_color),
                                       showlegend = FALSE)
    }

    spl <- spl %>% plotly::layout(xaxis = list(title = "X"),
                                  yaxis = list(title = "Y"),
                                  legend = list(x = 100, y = 0.5,font = list(family = "sans-serif",size = legend_text_size)))

  }


  # 3D plot
  else{
    axis_scale = match.arg(axis_scale, c("cube","real","custom"))
    ratio = plotly_axis_scale_3D(annotated_DT,sdimx = sdimx,sdimy = sdimy,sdimz = sdimz,
                                 mode = axis_scale,custom_ratio = custom_ratio)


    spl <- plotly::plot_ly(scene = "scene2")

    spl <- spl %>% plotly::add_trace(type = 'scatter3d',mode = 'markers',
                                     x = annotated_DT[[sdimx]],
                                     y = annotated_DT[[sdimy]],
                                     z = annotated_DT[[sdimz]],
                                     color = annotated_DT[[selected_genes]],
                                     colors = c(genes_low_color,genes_mid_color,genes_high_color),
                                     #legendgroup = annotated_DT[[cell_color]],
                                     marker = list(size = spatial_point_size),
                                     showlegend = F)
    if(show_other_cells == TRUE){
      spl <- spl %>% plotly::add_trace(type = "scatter3d",mode = "markers",
                                       x = annotated_DT_other[[sdimx]],
                                       y = annotated_DT_other[[sdimy]],
                                       z = annotated_DT_other[[sdimz]],
                                       marker = list(size = other_point_size,color = other_cell_color),
                                       showlegend = FALSE)
    }

    if(show_spatial_network == TRUE) {
      if(is.null(spatial_network)){
        stop("No usable spatial network specified! Please choose a network with spatial_network_name=xxx")
      }
      else{
        if(is.null(spatial_network_alpha)) {
          spatial_network_alpha = 0.5
        }
        else if(is.character(spatial_network_alpha)){
          warning("Edge_alpha for plotly mode is not adjustable yet. Default 0.5 will be set\n")
          spatial_network_alpha = 0.5
        }
        edges <- plotly_network(spatial_network)

        spl <- spl %>% plotly::add_trace(name = "sptial network",
                                         mode = "lines",
                                         type = "scatter3d",
                                         data = edges,
                                         x = ~x,y=~y,z=~z,
                                         line=list(color=network_color),
                                         opacity=spatial_network_alpha)
      }
    }

    if(show_spatial_grid == TRUE){
      cat("3D grid is not clear to view\n")
    }

  }



  spl <- plotly::hide_colorbar(spl)
  if(is.null(dim3_to_use) & is.null(sdimz)){
    if(plot_alignment == 'vertical'){
      combo_plot <- plotly::subplot(dpl,spl,nrows = 2,titleX = TRUE,titleY = TRUE)
    }
    else{
      combo_plot <- plotly::subplot(dpl,spl,titleX = TRUE,titleY = TRUE)
    }
  }

  else if(!is.null(dim3_to_use) & is.null(sdimz)){
    if(plot_alignment == 'vertical'){
      combo_plot <- plotly::subplot(dpl,spl,nrows = 2,titleX = TRUE,titleY = TRUE)%>%
        plotly::layout(scene = list(domain = list(x = c(0, 1), y = c(0,0.5)),
                                    xaxis = list(title = x_title),
                                    yaxis = list(title = y_title),
                                    zaxis = list(title = z_title)))
    }
    else{
      combo_plot <- plotly::subplot(dpl,spl,titleX = TRUE,titleY = TRUE) %>%
        plotly::layout(scene = list(domain = list(x = c(0, 0.5), y = c(0,1)),
                                    xaxis = list(title = x_title),
                                    yaxis = list(title = y_title),
                                    zaxis = list(title = z_title)))
    }
  }

  else if(is.null(dim3_to_use) & !is.null(sdimz)){
    if(plot_alignment == 'vertical'){
      combo_plot <- plotly::subplot(dpl,spl,nrows = 2,titleX = TRUE,titleY = TRUE)%>%
        plotly::layout(scene2 = list(
          xaxis = list(title = "X",nticks = x_ticks),
          yaxis = list(title = "Y",nticks = y_ticks),
          zaxis = list(title = "Z",nticks = z_ticks),
          aspectmode='manual',
          aspectratio = list(x=ratio[[1]],
                             y=ratio[[2]],
                             z=ratio[[3]])))
    }
    else{
      combo_plot <- plotly::subplot(dpl,spl,titleX = TRUE,titleY = TRUE) %>%
        plotly::layout(scene2 = list(
          xaxis = list(title = "X",nticks = x_ticks),
          yaxis = list(title = "Y",nticks = y_ticks),
          zaxis = list(title = "Z",nticks = z_ticks),
          aspectmode='manual',
          aspectratio = list(x=ratio[[1]],
                             y=ratio[[2]],
                             z=ratio[[3]])))
    }
  }

  else if(!is.null(dim3_to_use) & !is.null(sdimz)){
    if(plot_alignment == 'vertical'){
      combo_plot <- plotly::subplot(dpl,spl,nrows = 2,titleX = TRUE,titleY = TRUE)%>%
        plotly::layout(scene = list(domain = list(x = c(0, 1), y = c(0,0.5)),
                                    xaxis = list(title = x_title),
                                    yaxis = list(title = y_title),
                                    zaxis = list(title = z_title)),
                       scene2 = list(
                         xaxis = list(title = "X",nticks = x_ticks),
                         yaxis = list(title = "Y",nticks = y_ticks),
                         zaxis = list(title = "Z",nticks = z_ticks),
                         aspectmode='manual',
                         aspectratio = list(x=ratio[[1]],
                                            y=ratio[[2]],
                                            z=ratio[[3]])))
    }
    else{
      combo_plot <- plotly::subplot(dpl,spl,titleX = TRUE,titleY = TRUE) %>%
        plotly::layout(scene = list(domain = list(x = c(0, 0.5), y = c(0,1)),
                                    xaxis = list(title = x_title),
                                    yaxis = list(title = y_title),
                                    zaxis = list(title = z_title)),
                       scene2 = list(
                         xaxis = list(title = "X",nticks = x_ticks),
                         yaxis = list(title = "Y",nticks = y_ticks),
                         zaxis = list(title = "Z",nticks = z_ticks),
                         aspectmode='manual',
                         aspectratio = list(x=ratio[[1]],
                                            y=ratio[[2]],
                                            z=ratio[[3]])))
    }
  }

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
