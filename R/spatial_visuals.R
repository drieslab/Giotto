






## * ####
## 2-D ggplots ####
## ----------- ##

## ** image object compatibility ####

#' @title Optimized largeImage resampling
#' @name plot_auto_largeImage_resample
#' @description Downsample \code{largeImage} for plotting. Uses \code{\link[terra]{spatSample}}
#' to load only a portion of the original image, speeding up plotting and lowering memory
#' footprint.
#'
#' Default behavior of \code{spatSample} is to crop if only a smaller ROI is needed for plotting
#' followed by the sampling process in order to reduce wasted sampling by focusing the
#' sample space. For very large ROIs, this crop can be time intensive and require
#' writing to disk.
#'
#' This function examines the ROI dimensions as defined through the limits of the spatial
#' locations to be plotted, and decides between the following two methods in order to
#' avoid this issue:
#' \itemize{
#'   \item{\strong{Method A.} First crop original image and then sample n values where
#'   n = 500,000 to generate final image}
#'   \item{\strong{Method B.} First oversample n values and then crop, where n = 500,000
#'   scaled by a value >1. Scaling factor increases the smaller the ROI is and
#'   is defined by: original dimensions/crop dimensions where the larger ratio between
#'   x and y dims is chosen. Scale factor is capped by \code{max_resample_scale}}
#' }
#' Control points for this function are set by \code{max_crop} which decides the max
#' ROI area after which switchover to method B happens in order to avoid laborious crops
#' and \code{max_resample_scale} which determines the maximum scale factor for number
#' of values to sample. Both values can be adjusted depending on system resources.
#' Additionally, \code{flex_resample} determines if this switching behavior happens.
#' When set to \code{FALSE}, only method A is used.
#' @param gobject \code{gobject} containing \code{largeImage} object
#' @param giottoLargeImage \code{largeImage} object to resample if not provided through
#' \code{gobject} and \code{largeImage_name}
#' @param largeImage_name name of \code{largeImage} in \code{gobject}
#' @param spat_unit spatial unit
#' @param spat_loc_name name of spatial locations to plot
#' @param polygon_feat_type name of polygon/spatial_info to plot
#' @param include_image_in_border [boolean] expand the extent sampled to also show image in
#' border regions not included in spatlocs. This prevents images in plots from
#' being sharply cut off around the furthest spatial locations. (default is \code{TRUE})
#' @param flex_resample [boolean] Whether to allow automatic selection of sampling
#' workflow as defined in details sections. (default is \code{TRUE})
#' @param max_crop maximum crop size allowed for \strong{method A} before switching to
#' \strong{method B} (see description)
#' @param max_resample_scale maximum cells allowed to resample to compensate for
#' decreased resolution when cropping after sampling
#' @return a \code{giottoLargeImage} cropped and resampled properly for plotting
#' @seealso \code{\link[terra]{spatSample}}
#' @keywords internal
plot_auto_largeImage_resample = function(gobject,
                                         giottoLargeImage = NULL,
                                         largeImage_name = NULL,
                                         spat_unit = NULL,
                                         spat_loc_name = NULL,
                                         polygon_feat_type = NULL,
                                         include_image_in_border = TRUE,
                                         flex_resample = TRUE,
                                         max_crop = 1e+08,
                                         max_resample_scale = 100) {

  # If no giottoLargeImage, select specified giottoLargeImage. If none specified, select first one.
  if(is.null(giottoLargeImage)) {
    giottoLargeImage = get_giottoLargeImage(gobject = gobject,
                                            name = largeImage_name)
  }

  # Set spat_unit
  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = spat_unit)

  # Get spatial locations
  cell_locations = get_spatial_locations(gobject = gobject,
                                         spat_unit = spat_unit,
                                         spat_loc_name = spat_loc_name)

  # If no spatial locations are available, rely on first existing polygon extent
  if(is.null(cell_locations)) {
    sub_obj = get_polygon_info(gobject = gobject,
                               polygon_name = polygon_feat_type)
    # No point in checking feat_info right now. Subcellular gobj creation requires both poly and points
    # if(is.null(sub_obj)) {
    #   sub_obj = get_feature_info(gobject = gobject,
    #                              feat_type = feat_type)
    # }
    sub_ext = terra::ext(sub_obj)[1:4]
    cell_locations = data.table::data.table(sdimx = c(sub_ext[['xmin']], sub_ext[['xmax']]),
                                            sdimy = c(sub_ext[['ymin']], sub_ext[['ymax']]))
  }

  if(is.null(cell_locations)) stop('No spatial locations or polygons discovered.\n Cannot determine largeImage resample extent\n')

  # Get image extent minmax
  im_minmax = terra::ext(giottoLargeImage@raster_object)[1:4]
  # Determine crop
  if(isTRUE(include_image_in_border)) {
    # with crop padding
    x_halfPaddedRange = diff(range(cell_locations$sdimx))*0.625
    y_halfPaddedRange = diff(range(cell_locations$sdimy))*0.625
    x_mean = mean(cell_locations$sdimx)
    y_mean = mean(cell_locations$sdimy)

    xmax_crop = x_mean + x_halfPaddedRange
    xmin_crop = x_mean - x_halfPaddedRange
    ymax_crop = y_mean + y_halfPaddedRange
    ymin_crop = y_mean - y_halfPaddedRange

    if(xmin_crop < im_minmax[['xmin']]) xmin_crop = im_minmax[['xmin']]
    if(xmax_crop > im_minmax[['xmax']]) xmax_crop = im_minmax[['xmax']]
    if(ymin_crop < im_minmax[['ymin']]) ymin_crop = im_minmax[['ymin']]
    if(ymax_crop > im_minmax[['ymax']]) ymax_crop = im_minmax[['ymax']]
  } else {
    # no crop padding
    x_range = range(cell_locations$sdimx)
    y_range = range(cell_locations$sdimy)
    xmin_crop = x_range[1]
    xmax_crop = x_range[2]
    ymin_crop = y_range[1]
    ymax_crop = y_range[2]
  }

  # setup crop extent object
  crop_ext = terra::ext(xmin_crop, xmax_crop,
                        ymin_crop, ymax_crop)

  # zoom and resample giottoLargeImage
  crop_xdim = abs(xmax_crop - xmin_crop)
  crop_ydim = abs(ymax_crop - ymin_crop)
  crop_area_px = crop_xdim * giottoLargeImage@scale_factor[['x']] * giottoLargeImage@scale_factor[['y']] * crop_ydim
  im_xdim = abs(im_minmax[2] - im_minmax[1])
  im_ydim = abs(im_minmax[2] - im_minmax[1])
  crop_relative_scale = max(im_xdim/crop_xdim, im_ydim/crop_ydim)

  if(!isTRUE(flex_resample) | crop_area_px <= max_crop) {
    # METHOD A: Crop if needed then resample to final image
    if(!isTRUE(flex_resample) & crop_area_px > max_crop) {
      warning('Plotting large regions with flex_resample == FALSE will be costly in time and drive space.')
    }
    # For ROIs with area smaller than max_crop OR if flex_resample is FALSE
    giottoLargeImage@raster_object = terra::spatSample(giottoLargeImage@raster_object,
                                                       size = 500000,
                                                       method = 'regular',
                                                       as.raster = TRUE,
                                                       ext = crop_ext)
  } else {
    # METHOD B: Resample then crop to final image
    # Sample n values where n = default val scaled by a value >1
    if(crop_relative_scale <= max_resample_scale) {
      # Scale factor is fullsize image dim/crop dim. Larger of the two ratios is chosen
      resample_image = terra::spatSample(giottoLargeImage@raster_object,
                                         size = round(500000 * crop_relative_scale),
                                         method = 'regular',
                                         as.raster = TRUE)
    } else {
      # For scale factors that are too large, scaling is capped by max_resample_scale
      resample_image = terra::spatSample(giottoLargeImage@raster_object,
                                         size = 500000 * max_resample_scale,
                                         method = 'regular',
                                         as.raster = TRUE)
    }
    giottoLargeImage@raster_object = terra::crop(x = resample_image,
                                                 y = crop_ext)
  }

  return(giottoLargeImage)
}


## ** dim reduction plotting ####

#' @title plot_network_layer_ggplot
#' @name plot_network_layer_ggplot
#' @description Visualize cells in network layer according to dimension reduction coordinates
#' @param gobject giotto object
#' @param annotated_network_DT annotated network data.table of selected cells
#' @param edge_alpha alpha of network edges
#' @param show_legend show legend
#' @return ggplot
#' @details Description of parameters.
#' @keywords internal
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
#' @param point_alpha transparency of point
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
#' @keywords internal
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
                                   point_alpha = 1,
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
    pl = pl + ggplot2::geom_point(data = annotated_DT_other,
                                  ggplot2::aes_string(x = dims[1], dims[2]),
                                  color = other_cell_color,
                                  show.legend = F,
                                  size = other_point_size,
                                  alpha = point_alpha)

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
    pl <- pl + ggplot2::geom_point(data = annotated_DT_selected,
                                   ggplot2::aes_string(x = dims[1], dims[2]),
                                   color = cell_color,
                                   show.legend = show_legend,
                                   size = point_size,
                                   alpha = point_alpha)


  } else if(length(cell_color) > 1) {

    if(is.numeric(cell_color) | is.factor(cell_color)) {
      if(nrow(annotated_DT_selected) != length(cell_color)) stop('\n vector needs to be the same lengths as number of cells \n')
      annotated_DT_selected[['temp_color']] = cell_color

      pl <- pl + ggplot2::geom_point(data = annotated_DT_selected,
                                     aes_string2(x = dims[1], y = dims[2], fill = 'temp_color'),
                                     show.legend = show_legend, shape = 21,
                                     size = point_size,
                                     color = point_border_col,
                                     stroke = point_border_stroke,
                                     alpha = point_alpha)

    } else if(is.character(cell_color)) {
      if(!all(cell_color %in% grDevices::colors())) stop('cell_color is not numeric, a factor or vector of colors \n')
      pl <- pl + ggplot2::geom_point(data = annotated_DT_selected, aes_string2(x = dims[1], y = dims[2]),
                                     show.legend = show_legend, shape = 21, fill = cell_color,
                                     size = point_size,
                                     color = point_border_col,
                                     stroke = point_border_stroke,
                                     alpha = point_alpha)

    }

  } else if (is.character(cell_color)) {

    if(!cell_color %in% colnames(annotated_DT_selected)) {
      if(!cell_color %in% grDevices::colors()) stop(cell_color,' is not a color or a column name \n')
      pl <- pl + ggplot2::geom_point(data = annotated_DT_selected, ggplot2::aes_string(x = dims[1], y = dims[2]),
                                     show.legend = show_legend, shape = 21, fill = cell_color,
                                     size = point_size,
                                     color = point_border_col,
                                     stroke = point_border_stroke,
                                     alpha = point_alpha)

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
                                       color = point_border_col, stroke = point_border_stroke,
                                       alpha = point_alpha)

      } else {

        # convert character or numeric to factor
        if(color_as_factor == TRUE) {
          factor_data = factor(annotated_DT_selected[[cell_color]])
          annotated_DT_selected[[cell_color]] <- factor_data
        }

        # if you want to show centers or labels then calculate centers
        if(show_cluster_center == TRUE | show_center_label == TRUE) {
          annotated_DT_centers = annotated_DT_selected[, .(center_1 = stats::median(get(dims[1])),
                                                           center_2 = stats::median(get(dims[2]))), by = cell_color]
          factor_center_data = factor(annotated_DT_centers[[cell_color]])
          annotated_DT_centers[[cell_color]] <- factor_center_data
        }

        pl <- pl + ggplot2::geom_point(data = annotated_DT_selected,
                                       aes_string2(x = dims[1], y = dims[2], fill = cell_color),
                                       show.legend = show_legend, shape = 21, size = point_size,
                                       color = point_border_col, stroke = point_border_stroke,
                                       alpha = point_alpha)


        ## plot centers
        if(show_cluster_center == TRUE & (color_as_factor == TRUE | class_cell_color %in% c('character', 'factor'))) {

          pl <- pl + ggplot2::geom_point(data = annotated_DT_centers,
                                         aes_string2(x = 'center_1', y = 'center_2', fill = cell_color),
                                         color = center_point_border_col, stroke = center_point_border_stroke,
                                         size = center_point_size, shape = 21,
                                         alpha = point_alpha)
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
        cell_color_code = getDistinctColors(n = number_colors)
        names(cell_color_code) = unique(factor_data)
        pl <- pl + ggplot2::scale_fill_manual(values = cell_color_code)

      } else if(color_as_factor == F){

        if(is.null(gradient_midpoint)) {
          gradient_midpoint = stats::median(annotated_DT_selected[[cell_color]])
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
#' @param point_alpha transparancy of point
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
#' @keywords internal
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
                                   point_alpha = 1,
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
                                  color = other_cell_color, show.legend = F, size = other_point_size,
                                  alpha = point_alpha)

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
                                   color = cell_color, show.legend = show_legend, size = point_size,
                                   alpha = point_alpha)


  } else if(length(cell_color) > 1) {

    if(is.numeric(cell_color) | is.factor(cell_color)) {
      if(nrow(annotated_DT_selected) != length(cell_color)) stop('\n vector needs to be the same lengths as number of cells \n')
      annotated_DT_selected[['temp_color']] = cell_color

      pl <- pl + ggplot2::geom_point(data = annotated_DT_selected, aes_string2(x = dims[1], y = dims[2], color = 'temp_color'),
                                     show.legend = show_legend, shape = 19, size = point_size,
                                     alpha = point_alpha)

    } else if(is.character(cell_color)) {
      if(!all(cell_color %in% grDevices::colors())) stop('cell_color is not numeric, a factor or vector of colors \n')
      pl <- pl + ggplot2::geom_point(data = annotated_DT_selected, aes_string2(x = dims[1], y = dims[2]),
                                     show.legend = show_legend, shape = 19, fill = cell_color, size = point_size,
                                     alpha = point_alpha)

    }

  } else if (is.character(cell_color)) {

    if(!cell_color %in% colnames(annotated_DT_selected)) {
      if(!cell_color %in% grDevices::colors()) stop(cell_color,' is not a color or a column name \n')
      pl <- pl + ggplot2::geom_point(data = annotated_DT_selected, aes_string(x = dims[1], y = dims[2]),
                                     show.legend = show_legend, shape = 19, color = cell_color, size = point_size,
                                     alpha = point_alpha)

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
                                       show.legend = show_legend, shape = 19, size = point_size,
                                       alpha = point_alpha)

      } else {

        # convert character or numeric to factor
        if(color_as_factor == TRUE) {
          factor_data = factor(annotated_DT_selected[[cell_color]])
          annotated_DT_selected[[cell_color]] <- factor_data
        }

        # if you want to show centers or labels then calculate centers
        if(show_cluster_center == TRUE | show_center_label == TRUE) {
          annotated_DT_centers = annotated_DT_selected[, .(center_1 = stats::median(get(dims[1])),
                                                           center_2 = stats::median(get(dims[2]))), by = cell_color]
          factor_center_data = factor(annotated_DT_centers[[cell_color]])
          annotated_DT_centers[[cell_color]] <- factor_center_data
        }

        pl <- pl + ggplot2::geom_point(data = annotated_DT_selected,
                                       aes_string2(x = dims[1], y = dims[2], color = cell_color),
                                       show.legend = show_legend, shape = 19, size = point_size,
                                       alpha = point_alpha)


        ## plot centers
        if(show_cluster_center == TRUE & (color_as_factor == TRUE | class_cell_color %in% c('character', 'factor'))) {

          pl <- pl + ggplot2::geom_point(data = annotated_DT_centers,
                                         aes_string2(x = 'center_1', y = 'center_2', color = cell_color),
                                         size = center_point_size, shape = 19,
                                         alpha = point_alpha)
        }

        ## plot labels
        if(show_center_label == TRUE) {
          pl <- pl + ggrepel::geom_text_repel(data = annotated_DT_centers,
                                              aes_string2(x = 'center_1', y = 'center_2', label = cell_color),
                                              size = label_size, fontface = label_fontface,
                                              alpha = point_alpha)
        }

      }


      ## specificy colors to use
      if(!is.null(cell_color_code)) {

        pl <- pl + ggplot2::scale_color_manual(values = cell_color_code)

      } else if(color_as_factor == T) {

        number_colors = length(unique(factor_data))
        cell_color_code = getDistinctColors(n = number_colors)
        names(cell_color_code) = unique(factor_data)
        pl <- pl + ggplot2::scale_color_manual(values = cell_color_code)

      } else if(color_as_factor == F){

        if(is.null(gradient_midpoint)) {
          gradient_midpoint = stats::median(annotated_DT_selected[[cell_color]])
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
#' @param spat_unit spatial unit
#' @param feat_type feature type
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
#' @param point_alpha transparancy of point
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
#' @keywords internal
dimPlot2D_single <- function(gobject,
                             spat_unit = NULL,
                             feat_type = NULL,
                             dim_reduction_to_use = 'umap',
                             dim_reduction_name = NULL,
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
                             point_alpha = 1,
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


  # Set feat_type and spat_unit
  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = spat_unit)
  feat_type = set_default_feat_type(gobject = gobject,
                                    spat_unit = spat_unit,
                                    feat_type = feat_type)

  # specify dim_reduction_name according to provided feat_type
  if(!is.null(dim_reduction_to_use)) {
    if(is.null(dim_reduction_name)) {
      if(feat_type == 'rna') {
        dim_reduction_name = dim_reduction_to_use
      } else {
        dim_reduction_name = paste0(feat_type,'.',dim_reduction_to_use)
      }
    }
  }

  ## point shape ##
  point_shape = match.arg(point_shape, c('border', 'no_border'))

  ## dimension reduction ##
  # test if dimension reduction was performed

  dim_red_names = list_dim_reductions_names(gobject = gobject, data_type = 'cells',
                                            spat_unit = spat_unit, feat_type = feat_type,
                                            dim_type = dim_reduction_to_use)

  if(!dim_reduction_name %in% dim_red_names) {
    stop('\n dimension reduction: ', dim_reduction_to_use, ' or dimension reduction name: ',dim_reduction_name,' is not available \n')
  }

  #if(is.null(gobject@dimension_reduction$cells[[spat_unit]][[feat_type]][[dim_reduction_to_use]][[dim_reduction_name]])) {
  #  stop('\n dimension reduction: ', dim_reduction_to_use, ' or dimension reduction name: ',dim_reduction_name,' is not available \n')
  #}

  dim_dfr = get_dimReduction(gobject = gobject,
                             spat_unit = spat_unit,
                             reduction = 'cells',
                             reduction_method = dim_reduction_to_use,
                             name = dim_reduction_name,
                             return_dimObj = F)
  dim_dfr = dim_dfr[,c(dim1_to_use, dim2_to_use)]

  #dim_dfr = gobject@dimension_reduction$cells[[dim_reduction_to_use]][[dim_reduction_name]]$coordinates[,c(dim1_to_use, dim2_to_use)]

  dim_names = colnames(dim_dfr)

  # data.table variables
  cell_ID = NULL

  dim_DT = data.table::as.data.table(dim_dfr); dim_DT[, cell_ID := as.character(rownames(dim_dfr))]

  ## annotated cell metadata
  cell_metadata = combineMetadata(gobject = gobject,
                                  spat_unit = spat_unit,
                                  feat_type = feat_type,
                                  spat_enr_names = spat_enr_names,
                                  spat_loc_name = NULL)

  cell_metadata[, cell_ID := as.character(cell_ID)]

  annotated_DT = data.table::merge.data.table(cell_metadata, dim_DT, by = 'cell_ID')


  # create input for network
  if(show_NN_network == TRUE) {

    # nn_network
    selected_nn_network = get_NearestNetwork(gobject = gobject,
                                             spat_unit = spat_unit,
                                             nn_network_to_use = nn_network_to_use,
                                             network_name = network_name,
                                             output = 'igraph')

    #selected_nn_network = gobject@nn_network[[spat_unit]][[nn_network_to_use]][[network_name]][['igraph']]
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

    eigenvalues = gobject@dimension_reduction$cells[[spat_unit]][[dim_reduction_to_use]][[dim_reduction_name]]$misc$eigenvalues

    if(!is.null(eigenvalues)) {
      total = sum(eigenvalues)
      var_expl_vec = (eigenvalues/total) * 100
      dim1_x_variance = var_expl_vec[dim1_to_use]
      dim2_y_variance = var_expl_vec[dim2_to_use]

    }

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
                                 point_alpha = point_alpha,
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
                                        point_alpha = point_alpha,
                                        show_legend = show_legend)

  }




  ## add % variance explained to names of plot for PCA ##
  if(dim_reduction_to_use == 'pca') {

    if(!is.null(eigenvalues)) {
      x_name = paste0('pca','-',dim_names[1])
      y_name = paste0('pca','-',dim_names[2])

      # provide x, y and plot titles
      x_title = sprintf('%s explains %.02f%% of variance', x_name, var_expl_vec[dim1_to_use])
      y_title = sprintf('%s explains %.02f%% of variance', y_name, var_expl_vec[dim2_to_use])

      if(is.null(title)) title = cell_color
      pl <- pl + ggplot2::labs(x = x_title, y = y_title, title = title)
    }



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



#' @title Plot dimension reduction 2D
#' @name dimPlot2D
#' @description Visualize cells according to dimension reduction coordinates
#' @param gobject giotto object
#' @param spat_unit spatial unit
#' @param feat_type feature type
#' @param group_by create multiple plots based on cell annotation column
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
#' @param center_point_border_col border color of center points
#' @param center_point_border_stroke border stroke size of center points
#' @param label_size  size of labels
#' @param label_fontface font of labels
#' @param edge_alpha column to use for alpha of the edges
#' @param point_shape point with border or not (border or no_border)
#' @param point_size size of point (cell)
#' @param point_alpha transparancy of point
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
#' @param cow_rel_h cowplot param: relative heights of rows (e.g. c(1,2))
#' @param cow_rel_w cowplot param: relative widths of columns (e.g. c(1,2))
#' @param cow_align cowplot param: how to align
#' @param show_plot show plot
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters, see \code{\link{showSaveParameters}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return ggplot
#' @details Description of parameters. For 3D plots see \code{\link{dimPlot3D}}
#' @family reduced dimension visualizations
#' @export
dimPlot2D = function(gobject,
                     spat_unit = NULL,
                     feat_type = NULL,
                     group_by = NULL,
                     group_by_subset = NULL,
                     dim_reduction_to_use = 'umap',
                     dim_reduction_name = NULL,
                     dim1_to_use = 1,
                     dim2_to_use = 2,
                     spat_enr_names = NULL,
                     show_NN_network = FALSE,
                     nn_network_to_use = 'sNN',
                     network_name = 'sNN.pca',
                     cell_color = NULL,
                     color_as_factor = TRUE,
                     cell_color_code = NULL,
                     cell_color_gradient = c('blue', 'white', 'red'),
                     gradient_midpoint = NULL,
                     gradient_limits = NULL,
                     select_cell_groups = NULL,
                     select_cells = NULL,
                     show_other_cells = TRUE,
                     other_cell_color = 'lightgrey',
                     other_point_size = 0.5,
                     show_cluster_center = FALSE,
                     show_center_label = TRUE,
                     center_point_size = 4,
                     center_point_border_col = 'black',
                     center_point_border_stroke = 0.1,
                     label_size = 4,
                     label_fontface = 'bold',
                     edge_alpha = NULL,
                     point_shape = c('border', 'no_border'),
                     point_size = 1,
                     point_alpha = 1,
                     point_border_col = 'black',
                     point_border_stroke = 0.1,
                     title = NULL,
                     show_legend = TRUE,
                     legend_text = 10,
                     legend_symbol_size = 2,
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
                     spat_unit = spat_unit,
                     feat_type = feat_type,
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
                     point_alpha = point_alpha,
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
                                    spat_unit = spat_unit,
                                    feat_type = feat_type,
                                    spat_enr_names = spat_enr_names,
                                    spat_loc_name = NULL)
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
            cell_color_code = getDistinctColors(n = number_colors)
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
    savelist = list()


    for(group_id in 1:length(unique_groups)) {

      group = unique_groups[group_id]

      subset_cell_IDs = comb_metadata[get(group_by) == group][['cell_ID']]
      temp_gobject = subsetGiotto(gobject = gobject,
                                  spat_unit = spat_unit,
                                  feat_type = feat_type,
                                  cell_ids = subset_cell_IDs)

      pl = dimPlot2D_single(gobject = temp_gobject,
                            spat_unit = spat_unit,
                            feat_type = feat_type,
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
                            point_alpha = point_alpha,
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
#' @inheritDotParams dimPlot2D
#' @return ggplot
#' @details Description of parameters, see \code{\link{dimPlot2D}}. For 3D plots see \code{\link{dimPlot3D}}
#' @family reduced dimension visualizations
#' @export

dimPlot = function(...) {

  dimPlot2D(...)

}






#' @title plotUMAP_2D
#' @name plotUMAP_2D
#' @description Short wrapper for UMAP visualization
#' @param gobject giotto object
#' @param dim_reduction_name name of UMAP
#' @param default_save_name default save name of UMAP plot
#' @inheritDotParams dimPlot2D -gobject -dim_reduction_to_use -dim_reduction_name -default_save_name
#' @return ggplot
#' @details Description of parameters, see \code{\link{dimPlot2D}}. For 3D plots see \code{\link{plotUMAP_3D}}
#' @family reduced dimension visualizations
#' @export
plotUMAP_2D = function(gobject,
                       dim_reduction_name = NULL,
                       default_save_name = 'UMAP_2D',
                       ...) {

  dimPlot2D(gobject = gobject,
            dim_reduction_to_use = 'umap',
            dim_reduction_name = dim_reduction_name,
            default_save_name = default_save_name,
            ...)

}


#' @title plotUMAP
#' @name plotUMAP
#' @description Short wrapper for UMAP visualization
#' @param gobject giotto object
#' @param dim_reduction_name name of UMAP
#' @param default_save_name default save name of UMAP plot
#' @inheritDotParams dimPlot2D -gobject -dim_reduction_to_use -dim_reduction_name -default_save_name
#' @return ggplot
#' @details Description of parameters, see \code{\link{dimPlot2D}}. For 3D plots see \code{\link{plotUMAP_3D}}
#' @family reduced dimension visualizations
#' @export
plotUMAP = function(gobject,
                    dim_reduction_name = NULL,
                    default_save_name = 'UMAP',
                    ...) {

  dimPlot2D(gobject = gobject,
            dim_reduction_to_use = 'umap',
            dim_reduction_name = dim_reduction_name,
            default_save_name = default_save_name,
            ...)

}





#' @title plotTSNE_2D
#' @name plotTSNE_2D
#' @description Short wrapper for tSNE visualization
#' @param gobject giotto object
#' @param dim_reduction_name name of TSNE
#' @param default_save_name default save name of TSNE plot
#' @inheritDotParams dimPlot2D -gobject -dim_reduction_to_use -dim_reduction_name -default_save_name
#' @return ggplot
#' @details Description of parameters, see \code{\link{dimPlot2D}}. For 3D plots see \code{\link{plotTSNE_3D}}
#' @family reduced dimension visualizations
#' @export
plotTSNE_2D = function(gobject,
                       dim_reduction_name = NULL,
                       default_save_name = 'tSNE_2D',
                       ...) {

  dimPlot2D(gobject = gobject,
            dim_reduction_to_use = 'tsne',
            dim_reduction_name = dim_reduction_name,
            default_save_name = default_save_name,
            ...)

}

#' @title plotTSNE
#' @name plotTSNE
#' @description Short wrapper for tSNE visualization
#' @param gobject giotto object
#' @param dim_reduction_name name of TSNE
#' @param default_save_name default save name of TSNE plot
#' @inheritDotParams dimPlot2D -gobject -dim_reduction_to_use -dim_reduction_name -default_save_name
#' @return ggplot
#' @details Description of parameters, see \code{\link{dimPlot2D}}. For 3D plots see \code{\link{plotTSNE_3D}}
#' @family reduced dimension visualizations
#' @export
plotTSNE = function(gobject,
                    dim_reduction_name = NULL,
                    default_save_name = 'tSNE',
                    ...) {

  dimPlot2D(gobject = gobject,
            dim_reduction_to_use = 'tsne',
            dim_reduction_name = dim_reduction_name,
            default_save_name = default_save_name,
            ...)

}



#' @title plotPCA_2D
#' @name plotPCA_2D
#' @description Short wrapper for PCA visualization
#' @param gobject giotto object
#' @param dim_reduction_name name of PCA
#' @param default_save_name default save name of PCA plot
#' @inheritDotParams dimPlot2D -gobject -dim_reduction_to_use -dim_reduction_name -default_save_name
#' @return ggplot
#' @details Description of parameters, see \code{\link{dimPlot2D}}. For 3D plots see \code{\link{plotPCA_3D}}
#' @family reduced dimension visualizations
#' @export
plotPCA_2D = function(gobject,
                      dim_reduction_name = NULL,
                      default_save_name = 'PCA_2D',
                      ...) {

  dimPlot2D(gobject = gobject,
            dim_reduction_to_use = 'pca',
            dim_reduction_name = dim_reduction_name,
            default_save_name = default_save_name,
            ...)

}



#' @title plotPCA
#' @name plotPCA
#' @description Short wrapper for PCA visualization
#' @param gobject giotto object
#' @param dim_reduction_name name of PCA
#' @param default_save_name default save name of PCA plot
#' @inheritDotParams dimPlot2D -gobject -dim_reduction_to_use -dim_reduction_name -default_save_name
#' @return ggplot
#' @details Description of parameters, see \code{\link{dimPlot2D}}. For 3D plots see \code{\link{plotPCA_3D}}
#' @family reduced dimension visualizations
#' @export
plotPCA = function(gobject,
                   dim_reduction_name = NULL,
                   default_save_name = 'PCA',
                   ...) {

  dimPlot2D(gobject = gobject,
            dim_reduction_to_use = 'pca',
            dim_reduction_name = dim_reduction_name,
            default_save_name = default_save_name,
            ...)
}


## ** spatial plotting ####

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
#' @param point_alpha transparancy of point
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
#' @keywords internal
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
                                        point_alpha = 1,
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
                                        show_legend = TRUE) {

  ## specify spatial dimensions first
  if(is.null(sdimx) | is.null(sdimy)) {

    warning("plot_method = ggplot, but spatial dimensions for sdimx and/or sdimy are not specified. \n
            It will default to the 'sdimx' and 'sdimy' ")
    sdimx = 'sdimx'
    sdimy = 'sdimy'
  }


  ## point parameters ##
  point_size = as.numeric(point_size)
  point_alpha = as.numeric(point_alpha)
  point_border_stroke = as.numeric(point_border_stroke)
  center_point_size = as.numeric(center_point_size)
  center_point_border_stroke = as.numeric(center_point_border_stroke)
  label_size = as.numeric(label_size)
  other_point_size = as.numeric(other_point_size)

  ## ggplot object
  pl = ggobject

  ## first plot other non-selected cells
  if((!is.null(select_cells) | !is.null(select_cell_groups)) & show_other_cells == TRUE) {

    #print('OTHER CELLS WILL BE PLOTTED')

    pl <- pl + ggplot2::geom_point(data = cell_locations_metadata_other,
                                   aes_string(x = sdimx, sdimy),
                                   color = other_cell_color,
                                   show.legend = F,
                                   size = other_point_size, alpha = point_alpha)
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
                                   stroke = point_border_stroke, color = point_border_col,
                                   alpha = point_alpha)


  } else if(length(cell_color) > 1) {

    if(is.numeric(cell_color) | is.factor(cell_color)) {
      if(nrow(cell_locations_metadata_selected) != length(cell_color)) stop('\n vector needs to be the same lengths as number of cells \n')
      cell_locations_metadata_selected[['temp_color']] = cell_color

      pl <- pl + ggplot2::geom_point(data = cell_locations_metadata_selected,
                                     aes_string2(x = sdimx, y = sdimy, fill = 'temp_color'),
                                     show.legend = show_legend, shape = 21,
                                     size = point_size,
                                     color = point_border_col, stroke = point_border_stroke,
                                     alpha = point_alpha)

    } else if(is.character(cell_color)) {
      if(!all(cell_color %in% grDevices::colors())) stop('cell_color is not numeric, a factor or vector of colors \n')

      pl <- pl + ggplot2::geom_point(data = cell_locations_metadata_selected,
                                     aes_string2(x = sdimx, y = sdimy),
                                     show.legend = show_legend, shape = 21, fill = cell_color,
                                     size = point_size,
                                     color = point_border_col, stroke = point_border_stroke,
                                     alpha = point_alpha)
    }

  } else if(is.character(cell_color)) {
    if(!cell_color %in% colnames(cell_locations_metadata_selected)) {
      if(!cell_color %in% grDevices::colors()) stop(cell_color,' is not a color or a column name \n')
      pl <- pl + ggplot2::geom_point(data = cell_locations_metadata_selected,
                                     aes_string2(x = sdimx, y = sdimy),
                                     show.legend = show_legend, shape = 21, fill = cell_color,
                                     size = point_size,
                                     color = point_border_col, stroke = point_border_stroke,
                                     alpha = point_alpha)

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
                                       color = point_border_col,
                                       stroke = point_border_stroke,
                                       alpha = point_alpha)



      } else {

        # convert character or numeric to factor
        if(color_as_factor == TRUE) {
          factor_data = factor(cell_locations_metadata_selected[[cell_color]])
          cell_locations_metadata_selected[[cell_color]] <- factor_data
        }

        # if you want to show centers or labels then calculate centers
        if(show_cluster_center == TRUE | show_center_label == TRUE) {
          annotated_DT_centers = cell_locations_metadata_selected[, .(center_1 = stats::median(get('sdimx')),
                                                                      center_2 = stats::median(get('sdimy'))), by = cell_color]
          factor_center_data = factor(annotated_DT_centers[[cell_color]])
          annotated_DT_centers[[cell_color]] <- factor_center_data
        }

        pl <- pl + ggplot2::geom_point(data = cell_locations_metadata_selected,
                                       aes_string2(x = sdimx, y = sdimy, fill = cell_color),
                                       show.legend = show_legend, shape = 21, size = point_size,
                                       color = point_border_col, stroke = point_border_stroke,
                                       alpha = point_alpha)


        ## plot centers
        if(show_cluster_center == TRUE & (color_as_factor == TRUE | class_cell_color %in% c('character', 'factor'))) {

          pl <- pl + ggplot2::geom_point(data = annotated_DT_centers,
                                         aes_string2(x = 'center_1', y = 'center_2', fill = cell_color),
                                         color = center_point_border_col, stroke = center_point_border_stroke,
                                         size = center_point_size, shape = 21,
                                         alpha = point_alpha)
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
        cell_color_code = getDistinctColors(n = number_colors)
        names(cell_color_code) = unique(factor_data)
        pl <- pl + ggplot2::scale_fill_manual(values = cell_color_code)

      } else if(color_as_factor == F){

        if(is.null(gradient_midpoint)) {
          gradient_midpoint = stats::median(cell_locations_metadata_selected[[cell_color]])
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
#' @param point_alpha transparancy of point
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
#' @keywords internal
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
                                               point_alpha = 1,
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

  ## point parameters ##
  point_size = as.numeric(point_size)
  point_alpha = as.numeric(point_alpha)
  center_point_size = as.numeric(center_point_size)
  label_size = as.numeric(label_size)
  other_point_size = as.numeric(other_point_size)


  ## ggplot object
  pl = ggobject

  ## first plot other non-selected cells
  if((!is.null(select_cells) | !is.null(select_cell_groups)) & show_other_cells == TRUE) {
    pl <- pl + ggplot2::geom_point(data = cell_locations_metadata_other,
                                   aes_string(x = sdimx, sdimy),
                                   color = other_cell_color,
                                   show.legend = F, size = other_point_size, alpha = point_alpha)
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
                                   color = cell_color, size = point_size,
                                   alpha = point_alpha)


  } else if(length(cell_color) > 1) {

    if(is.numeric(cell_color) | is.factor(cell_color)) {
      if(nrow(cell_locations_metadata_selected) != length(cell_color)) stop('\n vector needs to be the same lengths as number of cells \n')
      cell_locations_metadata_selected[['temp_color']] = cell_color

      pl = pl + ggplot2::geom_point(data = cell_locations_metadata_selected, aes_string2(x = sdimx, y = sdimy, color = 'temp_color'),
                                     show.legend = show_legend, shape = 19, size = point_size, alpha = point_alpha)

    } else if(is.character(cell_color)) {
      if(!all(cell_color %in% grDevices::colors())) stop('cell_color is not numeric, a factor or vector of colors \n')

       pl = pl + ggplot2::geom_point(data = cell_locations_metadata_selected,
                                      aes_string2(x = sdimx, y = sdimy),
                                      show.legend = show_legend, shape = 19,
                                      color = cell_color, size = point_size,
                                      alpha = point_alpha)
    }

  } else if(is.character(cell_color)) {
    if(!cell_color %in% colnames(cell_locations_metadata_selected)) {
      if(!cell_color %in% grDevices::colors()) stop(cell_color,' is not a color or a column name \n')
      pl = pl + ggplot2::geom_point(data = cell_locations_metadata_selected,
                                     aes_string2(x = sdimx, y = sdimy),
                                     show.legend = show_legend, shape = 19, color = cell_color, size = point_size,
                                     alpha = point_alpha)

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
                                       show.legend = show_legend, shape = 19, size = point_size,
                                       alpha = point_alpha)



      } else {

        # convert character or numeric to factor
        if(color_as_factor == TRUE) {
          factor_data = factor(cell_locations_metadata_selected[[cell_color]])
          cell_locations_metadata_selected[[cell_color]] <- factor_data
        }

        # if you want to show centers or labels then calculate centers
        if(show_cluster_center == TRUE | show_center_label == TRUE) {
          annotated_DT_centers = cell_locations_metadata_selected[, .(center_1 = stats::median(get('sdimx')),
                                                                      center_2 = stats::median(get('sdimy'))), by = cell_color]
          factor_center_data = factor(annotated_DT_centers[[cell_color]])
          annotated_DT_centers[[cell_color]] <- factor_center_data
        }

        pl <- pl + ggplot2::geom_point(data = cell_locations_metadata_selected,
                                       aes_string2(x = sdimx, y = sdimy, color = cell_color),
                                       show.legend = show_legend, shape = 19, size = point_size,
                                       alpha = point_alpha)


        ## plot centers
        if(show_cluster_center == TRUE & (color_as_factor == TRUE | class_cell_color %in% c('character', 'factor'))) {

          pl <- pl + ggplot2::geom_point(data = annotated_DT_centers,
                                         aes_string2(x = 'center_1', y = 'center_2', color = cell_color),
                                         size = center_point_size, shape = 19, alpha = point_alpha)
        }

        ## plot labels
        if(show_center_label == TRUE) {
          pl <- pl + ggrepel::geom_text_repel(data = annotated_DT_centers,
                                              aes_string2(x = 'center_1', y = 'center_2', label = cell_color),
                                              size = label_size, fontface = label_fontface, alpha = point_alpha)
        }

      }

      ## specificy colors to use
      if(!is.null(cell_color_code)) {

        pl <- pl + ggplot2::scale_color_manual(values = cell_color_code)

      } else if(color_as_factor == T) {

        number_colors = length(unique(factor_data))
        cell_color_code = getDistinctColors(n = number_colors)
        names(cell_color_code) = unique(factor_data)
        pl <- pl + ggplot2::scale_color_manual(values = cell_color_code)

      } else if(color_as_factor == F){

        if(is.null(gradient_midpoint)) {
          gradient_midpoint = stats::median(cell_locations_metadata_selected[[cell_color]])
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



#' @title plot_spat_voronoi_layer_ggplot
#' @name plot_spat_voronoi_layer_ggplot
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
#' @param point_alpha transparancy of point
#' @param show_cluster_center plot center of selected clusters
#' @param show_center_label plot label of selected clusters
#' @param center_point_size size of center points
#' @param label_size  size of labels
#' @param label_fontface font of labels
#' @param show_other_cells display not selected cells
#' @param other_cell_color color for not selected cells
#' @param other_point_size point size for not selected cells
#' @param background_color background color
#' @param vor_border_color borde colorr of voronoi plot
#' @param vor_max_radius maximum radius for voronoi 'cells'
#' @param vor_alpha transparancy of voronoi 'cells'
#' @param show_legend show legend
#' @return ggplot
#' @details Description of parameters.
#' @keywords internal
#' @export
plot_spat_voronoi_layer_ggplot = function(ggobject,
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
                                          point_alpha = 1,
                                          show_cluster_center = F,
                                          show_center_label = T,
                                          center_point_size = 4,
                                          label_size = 4,
                                          label_fontface = 'bold',
                                          show_other_cells = T,
                                          other_cell_color = 'lightgrey',
                                          other_point_size = 1,
                                          background_color = 'white',
                                          vor_border_color = 'white',
                                          vor_max_radius = 200,
                                          vor_alpha = 1,
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



  ## order of color
  # 1. if NULL then default to lightblue
  # 2. if character vector
  # 2.1 if length of cell_color is longer than 1 and has colors
  # 2.2 if not part of metadata then suppose its color
  # 2.3 part of metadata
  # 2.3.1 numerical column
  # 2.3.2 factor column or character to factor

  # data.table variables
  temp_color = NULL

  # cell color default
  if(is.null(cell_color)) {

    ## 1. default colors when no colors are assigned ##

    cell_color = 'lightblue'
    cell_locations_metadata_selected[, 'temp_color' := 'selected']

    if(!is.null(cell_locations_metadata_other)) cell_locations_metadata_other[, 'temp_color' := 'other']


    combn_cell_locations_metadata = rbind(cell_locations_metadata_selected, cell_locations_metadata_other)

    pl = pl + ggforce::geom_voronoi_tile(data = combn_cell_locations_metadata,
                                         aes(x = sdimx, y = sdimy, group = -1L, fill = as.factor(temp_color)),
                                         colour = vor_border_color, max.radius = vor_max_radius, show.legend = show_legend,
                                         alpha = vor_alpha)

    if(show_other_cells == TRUE) {
      pl = pl + ggplot2::scale_fill_manual(values = c(selected = cell_color, other = other_cell_color))
    } else {
      pl = pl + ggplot2::scale_fill_manual(values = c(selected = cell_color, other = background_color))
    }

    # theme specific changes
    pl = pl + theme(legend.title = element_blank())



  } else if(length(cell_color) > 1) {

    ## 2. continuous vector to convert to colors ##
    if(is.numeric(cell_color) | is.factor(cell_color)) {
      if(nrow(cell_locations_metadata_selected) != length(cell_color)) stop('\n vector needs to be the same lengths as number of cells \n')

      cell_locations_metadata_selected[['temp_color']] = cell_color
      if(!is.null(cell_locations_metadata_other)) cell_locations_metadata_other[['temp_color']] = NA
      combn_cell_locations_metadata = rbind(cell_locations_metadata_selected, cell_locations_metadata_other)

      pl = pl + ggforce::geom_voronoi_tile(data = combn_cell_locations_metadata,
                                           aes(x = sdimx, y = sdimy, group = -1L, fill = temp_color),
                                           colour = vor_border_color, max.radius = vor_max_radius, show.legend = show_legend,
                                           alpha = vor_alpha)

      if(is.null(gradient_midpoint)) {
        gradient_midpoint = stats::median(cell_locations_metadata_selected[['temp_color']])
      }

      mybg_color = ifelse(show_other_cells == TRUE, other_cell_color, background_color)

      pl <- pl + ggplot2::scale_fill_gradient2(low = cell_color_gradient[[1]],
                                               mid = cell_color_gradient[[2]],
                                               high = cell_color_gradient[[3]],
                                               midpoint = gradient_midpoint,
                                               na.value = mybg_color)

      # theme specific changes
      pl = pl + theme(legend.title = element_blank())


    } else if(is.character(cell_color)) {

      ## 3. character vector to convert to colors ##

      if(!all(cell_color %in% grDevices::colors())) stop('cell_color is not numeric, a factor or vector of colors \n')

      if(nrow(cell_locations_metadata_selected) != length(cell_color)) stop('\n vector needs to be the same lengths as number of cells \n')

      other_cell_color = ifelse(show_other_cells == TRUE, other_cell_color, background_color)

      cell_locations_metadata_selected[['temp_color']] = cell_color
      if(!is.null(cell_locations_metadata_other)) cell_locations_metadata_other[['temp_color']] = other_cell_color
      combn_cell_locations_metadata = rbind(cell_locations_metadata_selected, cell_locations_metadata_other)

      pl = pl + ggforce::geom_voronoi_tile(data = combn_cell_locations_metadata,
                                           aes(x = sdimx, y = sdimy, group = -1L, fill = temp_color),
                                           colour = vor_border_color, max.radius = vor_max_radius, show.legend = show_legend,
                                           alpha = vor_alpha)

      my_color_code = unique(combn_cell_locations_metadata[['temp_color']])
      names(my_color_code) = my_color_code

      pl <- pl + ggplot2::scale_fill_manual(values = my_color_code)

      # theme specific changes
      pl = pl + theme(legend.title = element_blank())

    }




  } else if(is.character(cell_color)) {
    if(!cell_color %in% colnames(cell_locations_metadata_selected)) {
      if(!cell_color %in% grDevices::colors()) stop(cell_color,' is not a color or a column name \n')

      ## 4. use a specific color ##
      other_cell_color = ifelse(show_other_cells == TRUE, other_cell_color, background_color)

      cell_locations_metadata_selected[['temp_color']] = 'selected'
      if(!is.null(cell_locations_metadata_other)) cell_locations_metadata_other[['temp_color']] = 'other'
      combn_cell_locations_metadata = rbind(cell_locations_metadata_selected, cell_locations_metadata_other)

      pl = pl + ggforce::geom_voronoi_tile(data = combn_cell_locations_metadata,
                                           aes(x = sdimx, y = sdimy, group = -1L, fill = temp_color),
                                           colour = vor_border_color, max.radius = vor_max_radius, show.legend = show_legend,
                                           alpha = vor_alpha)

      my_color_code = unique(combn_cell_locations_metadata[['temp_color']])
      names(my_color_code) = my_color_code
      pl = pl + ggplot2::scale_fill_manual(values = c(selected = cell_color, other = other_cell_color))

      # theme specific changes
      pl = pl + theme(legend.title = element_blank())

    } else {

      class_cell_color = class(cell_locations_metadata_selected[[cell_color]])

      if((class_cell_color == 'integer' | class_cell_color == 'numeric') & color_as_factor == FALSE) {

        ## 5. use continuous column from metadata ##

        # set upper and lower limits
        if(!is.null(gradient_limits) & is.vector(gradient_limits) & length(gradient_limits) == 2) {
          lower_lim = gradient_limits[[1]]
          upper_lim = gradient_limits[[2]]

          numeric_data = cell_locations_metadata_selected[[cell_color]]

          limit_numeric_data = ifelse(numeric_data > upper_lim, upper_lim,
                                      ifelse(numeric_data < lower_lim, lower_lim, numeric_data))
          cell_locations_metadata_selected[[cell_color]] = limit_numeric_data

        }

        cell_locations_metadata_selected[['temp_color']] = cell_locations_metadata_selected[[cell_color]]
        if(!is.null(cell_locations_metadata_other)) cell_locations_metadata_other[['temp_color']] = NA
        combn_cell_locations_metadata = rbind(cell_locations_metadata_selected, cell_locations_metadata_other)

        pl = pl + ggforce::geom_voronoi_tile(data = combn_cell_locations_metadata,
                                             aes(x = sdimx, y = sdimy, group = -1L, fill = temp_color),
                                             colour = vor_border_color, max.radius = vor_max_radius, show.legend = show_legend,
                                             alpha = vor_alpha)

        mybg_color = ifelse(show_other_cells == TRUE, other_cell_color, background_color)

        if(is.null(gradient_midpoint)) {
          gradient_midpoint = stats::median(cell_locations_metadata_selected[['temp_color']])
        }

        pl = pl + ggplot2::scale_fill_gradient2(low = cell_color_gradient[[1]],
                                                mid = cell_color_gradient[[2]],
                                                high = cell_color_gradient[[3]],
                                                midpoint = gradient_midpoint,
                                                na.value = mybg_color,
                                                name = cell_color)



      } else {


        ## 6. use factor or character column from metadata ##
        # convert character or numeric to factor
        if(color_as_factor == TRUE) {
          factor_data = factor(cell_locations_metadata_selected[[cell_color]])
          cell_locations_metadata_selected[[cell_color]] <- factor_data
        }

        # if you want to show centers or labels then calculate centers
        if(show_cluster_center == TRUE | show_center_label == TRUE) {
          annotated_DT_centers = cell_locations_metadata_selected[, .(center_1 = stats::median(get('sdimx')),
                                                                      center_2 = stats::median(get('sdimy'))), by = cell_color]
          factor_center_data = factor(annotated_DT_centers[[cell_color]])
          annotated_DT_centers[[cell_color]] <- factor_center_data
        }

        cell_locations_metadata_selected[['temp_color']] = cell_locations_metadata_selected[[cell_color]]
        if(!is.null(cell_locations_metadata_other)) cell_locations_metadata_other[['temp_color']] = 'other'
        combn_cell_locations_metadata = rbind(cell_locations_metadata_selected, cell_locations_metadata_other)

        pl = pl + ggforce::geom_voronoi_tile(data = combn_cell_locations_metadata,
                                             aes(x = sdimx, y = sdimy, group = -1L, fill = temp_color),
                                             colour = vor_border_color, max.radius = vor_max_radius, show.legend = show_legend,
                                             alpha = vor_alpha)


        other_cell_color = ifelse(show_other_cells == TRUE, other_cell_color, background_color)

        ## specificy colors to use
        if(!is.null(cell_color_code)) {

          cell_color_code[['other']] = other_cell_color
          pl = pl + ggplot2::scale_fill_manual(values = cell_color_code,
                                               name = cell_color)

        } else if(color_as_factor == T) {

          number_colors = length(unique(factor_data))
          cell_color_code = getDistinctColors(n = number_colors)
          names(cell_color_code) = unique(factor_data)

          cell_color_code[['other']] = other_cell_color
          pl = pl + ggplot2::scale_fill_manual(values = cell_color_code, name = cell_color)

        }

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


    }
  }



  ## lastly overlay POINTS ##
  ## first plot other non-selected cells
  if((!is.null(select_cells) | !is.null(select_cell_groups)) & show_other_cells == TRUE) {

    pl <- pl + ggplot2::geom_point(data = cell_locations_metadata_other,
                                   aes_string(x = sdimx, sdimy),
                                   color = 'black', show.legend = F, size = other_point_size,
                                   alpha = point_alpha)
  }

  ## plot selected cells
  pl <- pl + ggplot2::geom_point(data = cell_locations_metadata_selected,
                                 aes_string(x = sdimx, y = sdimy),
                                 show.legend = F, color = 'black', size = point_size,
                                 alpha = point_alpha)


  return(pl)
}





#' @title plot_spat_image_layer_ggplot
#' @name plot_spat_image_layer_ggplot
#' @description create background image in ggplot
#' @param gg_obj ggplot2 object
#' @param gobject giotto object
#' @param gimage a giotto image or a list/vector of giotto images
#' @param feat_type feature type
#' @param spat_unit spatial unit
#' @param spat_loc_name name for spatial locations
#' @param polygon_feat_type name for feature type associated with polygon information
#' @param sdimx x-axis dimension name (default = 'sdimx')
#' @param sdimy y-axis dimension name (default = 'sdimy')
#' @return ggplot
#' @keywords internal
plot_spat_image_layer_ggplot = function(gg_obj,
                                        gobject,
                                        gimage,
                                        feat_type = NULL,
                                        spat_unit = NULL,
                                        spat_loc_name = NULL,
                                        polygon_feat_type = NULL,
                                        sdimx = NULL,
                                        sdimy = NULL) {


  if(is.null(gobject) | is.null(gimage)) {
    stop('A giotto object and a giotto image need to be provided')
  }

  if(is.null(sdimx) | is.null(sdimy)) {
    warning("plot_method = ggplot, but spatial dimensions for sdimx and/or sdimy are not specified. \n
            It will default to the 'sdimx' and 'sdimy' ")
    sdimx = 'sdimx'
    sdimy = 'sdimy'
  }

  # Set feat_type and spat_unit
  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = spat_unit)
  feat_type = set_default_feat_type(gobject = gobject,
                                    spat_unit = spat_unit,
                                    feat_type = feat_type)

  # spatial locations
  spatlocs = get_spatial_locations(gobject = gobject,
                                   spat_unit = spat_unit,
                                   spat_loc_name = spat_loc_name)

  # Get spatial extent for positioning purposes
  spat_ext = spatlocs[, c('sdimx', 'sdimy'), with = FALSE]

  # When spatial locations are missing but subcellular info is present
  # Pull plot extent info from polygon info if present

  if(is.null(spat_ext)) {
    gpoly = get_polygon_info(gobject = gobject,
                             polygon_name = polygon_feat_type)

    poly_ext = terra::ext(gpoly)[1:4]
    spat_ext = data.table::data.table(sdimx = c(poly_ext[['xmin']], poly_ext[['xmax']]),
                                      sdimy = c(poly_ext[['ymin']], poly_ext[['ymax']]))
  }

  # If still missing, send warning
  if(is.null(spat_ext)) {
    warning('No spatial locations or polygon info found.\n Plot spatial extent may be incorrect\n')
  }

  # Assign region to plot
  gg_obj = gg_obj + geom_blank(data = spat_ext, aes_string(sdimx, sdimy))

  if((is.list(gimage) | is.vector(gimage)) & length(gimage) > 1) {

    for(i in 1:length(gimage)) {

      if(inherits(gimage[[i]], 'giottoImage')) {
        # extract min and max from object
        my_xmax = gimage[[i]]@minmax[1]
        my_xmin = gimage[[i]]@minmax[2]
        my_ymax = gimage[[i]]@minmax[3]
        my_ymin = gimage[[i]]@minmax[4]

        # convert giotto image object into array
        img_array = as.numeric(gimage[[i]]@mg_object[[1]])

        # extract adjustments from object
        xmax_b = gimage[[i]]@boundaries[1]
        xmin_b = gimage[[i]]@boundaries[2]
        ymax_b = gimage[[i]]@boundaries[3]
        ymin_b = gimage[[i]]@boundaries[4]

        gg_obj = gg_obj + annotation_raster(img_array,
                                            xmin = my_xmin-xmin_b, xmax = my_xmax+xmax_b,
                                            ymin = my_ymin-ymin_b, ymax = my_ymax+ymax_b)

        # TODO geom_raster to accommodate single-channel

      } else if(inherits(gimage[[i]], 'giottoLargeImage')) {
        # get plotting minmax
        extent = terra::ext(gimage[[i]]@raster_object)[1:4]
        xmin = extent[['xmin']]
        xmax = extent[['xmax']]
        ymin = extent[['ymin']]
        ymax = extent[['ymax']]

        # convert raster object into array with 3 channels
        img_array = terra::as.array(gimage[[i]]@raster_object)
        img_array = img_array/max(img_array)
        if(dim(img_array)[3] == 1) {
          img_array_RGB = array(NA, dim = c(dim(img_array)[1:2],3))
          img_array_RGB[,,1:3] = img_array
        } else {
          img_array_RGB = img_array
        }

        gg_obj = gg_obj + annotation_raster(img_array_RGB,
                                            xmin = xmin, xmax = xmax,
                                            ymin = ymin, ymax = ymax)

        # TODO geom_raster to accommodate single-channel
      }

    }

  } else {

    if(methods::is(gimage, 'giottoImage')) {
      # extract min and max from object
      my_xmax = gimage@minmax[1]
      my_xmin = gimage@minmax[2]
      my_ymax = gimage@minmax[3]
      my_ymin = gimage@minmax[4]

      # convert giotto image object into array
      img_array = as.numeric(gimage@mg_object[[1]])

      # extract adjustments from object
      xmax_b = gimage@boundaries[1]
      xmin_b = gimage@boundaries[2]
      ymax_b = gimage@boundaries[3]
      ymin_b = gimage@boundaries[4]

      gg_obj = gg_obj + annotation_raster(img_array,
                                          xmin = my_xmin-xmin_b, xmax = my_xmax+xmax_b,
                                          ymin = my_ymin-ymin_b, ymax = my_ymax+ymax_b)

    } else if(methods::is(gimage, 'giottoLargeImage')) {
      # get plotting minmax
      extent = terra::ext(gimage@raster_object)[1:4]
      xmin = extent[['xmin']]
      xmax = extent[['xmax']]
      ymin = extent[['ymin']]
      ymax = extent[['ymax']]

      # convert raster object into array with 3 channels
      img_array = terra::as.array(gimage@raster_object)
      img_array = img_array/max(img_array)
      if(dim(img_array)[3] == 1) {
        img_array_RGB = array(NA, dim = c(dim(img_array)[1:2],3))
        img_array_RGB[,,1:3] = img_array
      } else {
        img_array_RGB = img_array
      }

      gg_obj = gg_obj + annotation_raster(img_array_RGB,
                                          xmin = xmin, xmax = xmax,
                                          ymin = ymin, ymax = ymax)
    }

  }

  if(!is.null(spatlocs))  gg_obj = gg_obj + geom_point(data = spatlocs, aes_string(sdimx, sdimy), alpha = 0.5, size = 0.4)

  return(gg_obj)

}




#' @title plot_spat_scatterpie_layer_ggplot
#' @name plot_spat_scatterpie_layer_ggplot
#' @description create scatterpie in ggplot
#' @param gobject giotto object
#' @param sdimx x-axis dimension name (default = 'sdimx')
#' @param sdimy y-axis dimension name (default = 'sdimy')
#' @param spatial_locations spatial locations
#' @param spatial_enrichment spatial enrichment results
#' @param radius radius of scatterpie
#' @param color color of lines within pie chart
#' @param cell_color_code color code for the cell types
#' @return ggplot
#' @keywords internal
plot_spat_scatterpie_layer_ggplot = function(ggobject,
                                             sdimx = 'sdimx',
                                             sdimy = 'sdimy',
                                             spatial_locations = NULL,
                                             spatial_enrichment = NULL,
                                             radius = 10,
                                             color = NA,
                                             alpha = 1,
                                             cell_color_code = NULL) {


  # get cell names
  cell_names = colnames(spatial_enrichment)[-1]

  # combine spatial locations and enrichment results
  combined_spat_enrichm = data.table::merge.data.table(x = spatial_enrichment,
                                                       y = spatial_locations,
                                                       by = 'cell_ID')

  # plot scatterpie
  pl = ggobject
  pl = pl + scatterpie::geom_scatterpie(data = combined_spat_enrichm,
                                        aes(x = sdimx, y = sdimy, r = radius),
                                        cols = cell_names,
                                        color = color,
                                        alpha = alpha)

  ## specificy colors to use
  if(!is.null(cell_color_code)) {
    pl = pl + ggplot2::scale_fill_manual(values = cell_color_code)
  } else {

    number_colors = length(unique(cell_names))
    cell_color_code = getDistinctColors(n = number_colors)
    names(cell_color_code) = unique(cell_names)
    pl = pl + ggplot2::scale_fill_manual(values = cell_color_code)

  }

  return(pl)

}








#' @title spatPlot2D_single
#' @name spatPlot2D_single
#' @description Visualize cells according to spatial coordinates
#' @param gobject giotto object
#' @param feat_type feature type
#' @param spat_unit spatial unit
#' @param show_image show a tissue background image
#' @param gimage a giotto image
#' @param image_name name of a giotto image
#' @param largeImage_name name of a giottoLargeImage
#' @param spat_loc_name name of spatial locations
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
#' @param point_shape shape of points (border, no_border or voronoi)
#' @param point_size size of point (cell)
#' @param point_alpha transparancy of point
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
#' @param vor_border_color border colorr for voronoi plot
#' @param vor_max_radius maximum radius for voronoi 'cells'
#' @param vor_alpha transparancy of voronoi 'cells'
#' @param axis_text size of axis text
#' @param axis_title size of axis title
#' @param show_plot show plot
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param verbose be verbose
#' @param save_param list of saving parameters from \code{\link{all_plots_save_function}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return ggplot
#' @details Description of parameters.
#' @keywords internal
#' @seealso \code{\link{spatPlot3D}}
spatPlot2D_single = function(gobject,
                             feat_type = NULL,
                             spat_unit = NULL,
                             show_image = F,
                             gimage = NULL,
                             image_name = NULL,
                             largeImage_name = NULL,
                             spat_loc_name = NULL,
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
                             point_shape = c('border', 'no_border', 'voronoi'),
                             point_size = 3,
                             point_alpha = 1,
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
                             title = NULL,
                             show_legend = T,
                             legend_text = 8,
                             legend_symbol_size = 1,
                             background_color = 'white',
                             vor_border_color = 'white',
                             vor_max_radius = 200,
                             vor_alpha = 1,
                             axis_text = 8,
                             axis_title = 8,
                             show_plot = NA,
                             return_plot = NA,
                             save_plot = NA,
                             verbose = FALSE,
                             save_param = list(),
                             default_save_name = 'spatPlot2D_single') {


  if(verbose == TRUE) {
    cat('\n verbose == TRUE \n')
  }

  # Check params
  if(!is.null(image_name) && !is.null(largeImage_name)) stop('Only one type of image can be used at a time')

  # Set feat_type and spat_unit
  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = spat_unit)
  feat_type = set_default_feat_type(gobject = gobject,
                                    spat_unit = spat_unit,
                                    feat_type = feat_type)

  ## giotto image ##
  if(show_image == TRUE) {

    if(!is.null(gimage)) {
      gimage = gimage
    } else if(!is.null(image_name)) {
      # If there is input to image_name arg

      if(length(image_name) == 1) {
        gimage = gobject@images[[image_name]]
        if(is.null(gimage)) warning('image_name: ', image_name, ' does not exist \n')
      } else {
        gimage = list()
        for(gim in 1:length(image_name)) {
          gimage[[gim]] = gobject@images[[gim]]
          if(is.null(gimage[[gim]])) warning('image_name: ', gim, ' does not exists \n')
        }
      }

    } else if(!is.null(largeImage_name)) {
      # If there is input to largeImage_name arg

      if(length(largeImage_name) == 1) {
        gimage = plot_auto_largeImage_resample(gobject = gobject,
                                               largeImage_name = largeImage_name,
                                               spat_unit = spat_unit,
                                               spat_loc_name = spat_loc_name,
                                               include_image_in_border = TRUE)
      } else {
        gimage = list()
        for(gim in 1:length(largeImage_name)) {
          gimage[[gim]] = plot_auto_largeImage_resample(gobject = gobject,
                                                        largeImage_name = largeImage_name[[gim]],
                                                        spat_unit = spat_unit,
                                                        spat_loc_name = spat_loc_name,
                                                        include_image_in_border = TRUE)
        }
      }

    } else {
      # Default to first image available in images if no input given to image_name or largeImage_name args
      image_name = names(gobject@images)[1]
      gimage = gobject@images[[image_name]]
      if(is.null(gimage)) warning('image_name: ', image_name, ' does not exist \n')
    }
  }


  ## point shape ##
  point_shape = match.arg(point_shape, choices = c('border', 'no_border', 'voronoi'))

  ## get spatial cell locations
  cell_locations = get_spatial_locations(gobject = gobject,
                                         spat_unit = spat_unit,
                                         spat_loc_name = spat_loc_name)


  ## extract spatial network
  if(show_network == TRUE) {
    spatial_network = get_spatialNetwork(gobject,
                                         spat_unit = spat_unit,
                                         name = spatial_network_name,
                                         return_network_Obj = FALSE)
  } else {
    spatial_network = NULL
  }

  ## extract spatial grid
  if(show_grid == TRUE) {
    spatial_grid = get_spatialGrid(gobject = gobject,
                                   spat_unit = spat_unit,
                                   feat_type = feat_type,
                                   name = spatial_grid_name)
  } else {
    spatial_grid = NULL
  }


  ## get cell metadata

  if(is.null(spat_loc_name)) {
    if(!is.null(gobject@spatial_locs)) {
      spat_loc_name = names(gobject@spatial_locs[[spat_unit]])[[1]]
      # cat('No spatial locations have been selected, the first one -',spat_loc_name, '- will be used \n')
    } else {
      spat_loc_name = NULL
      cat('No spatial locations have been found \n')
      return(NULL)
    }
  }

  cell_metadata = combineMetadata(gobject = gobject,
                                  feat_type = feat_type,
                                  spat_unit = spat_unit,
                                  spat_loc_name = spat_loc_name,
                                  spat_enr_names = spat_enr_names)

  if(nrow(cell_metadata) == 0) {
    cell_locations_metadata = cell_locations
  } else {
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
    spatial_network = spatial_network[spatial_network$to %in% select_cells & spatial_network$from %in% select_cells]

    # if specific cells are selected
    # cell_locations_metadata = cell_locations_metadata_selected

  } else if(is.null(select_cells)) {

    cell_locations_metadata_selected = cell_locations_metadata
    cell_locations_metadata_other = NULL

  }

  # update cell_color_code
  # only keep names from selected groups
  if(!is.null(select_cell_groups) & !is.null(cell_color_code)) {
    cell_color_code = cell_color_code[names(cell_color_code) %in% select_cell_groups]
  }

  # data.table and ggplot variables
  sdimx_begin = sdimy_begin = sdimx_end = sdimy_end = x_start = x_end = y_start = y_end = NULL


  ### create 2D plot with ggplot ###
  #cat('create 2D plot with ggplot \n')

  if(verbose == TRUE) {
    cat('Data table with selected information (e.g. cells): \n')
    print(cell_locations_metadata_selected[1:5,])

    cat('Data table with non-selected information (e.g. cells): \n')
    print(cell_locations_metadata_other[1:5,])
  }




  pl = ggplot2::ggplot()
  pl = pl + ggplot2::theme_bw()

  ## plot image ##
  if(show_image == TRUE & !is.null(gimage)) {
    pl = plot_spat_image_layer_ggplot(gg_obj = pl,
                                      gobject = gobject,
                                      feat_type = feat_type,
                                      spat_unit = spat_unit,
                                      spat_loc_name = spat_loc_name,
                                      gimage = gimage,
                                      sdimx = sdimx,
                                      sdimy = sdimy)
  }

  ## plot spatial network
  if(!is.null(spatial_network) & show_network == TRUE) {
    if(is.null(network_color)) network_color = 'red'
    pl <- pl + ggplot2::geom_segment(data = spatial_network, aes(x = sdimx_begin,
                                                                 y = sdimy_begin,
                                                                 xend = sdimx_end,
                                                                 yend = sdimy_end),
                                     color = network_color, size = 0.5, alpha = network_alpha)
  }


  ## plot spatial grid
  if(!is.null(spatial_grid) & show_grid == TRUE) {
    if(is.null(grid_color)) grid_color = 'black'
    pl <- pl + ggplot2::geom_rect(data = spatial_grid, aes(xmin = x_start,
                                                           xmax = x_end,
                                                           ymin = y_start,
                                                           ymax = y_end),
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
                                      point_alpha = point_alpha,
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
                                             point_alpha = point_alpha,
                                             show_cluster_center = show_cluster_center,
                                             show_center_label = show_center_label,
                                             center_point_size = center_point_size,
                                             label_size = label_size,
                                             label_fontface = label_fontface,
                                             show_other_cells = show_other_cells,
                                             other_cell_color = other_cell_color,
                                             other_point_size = other_point_size,
                                             show_legend = show_legend)

  } else if(point_shape == 'voronoi') {

    pl = plot_spat_voronoi_layer_ggplot(ggobject = pl,
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
                                        point_alpha = point_alpha,
                                        show_cluster_center = show_cluster_center,
                                        show_center_label = show_center_label,
                                        center_point_size = center_point_size,
                                        label_size = label_size,
                                        label_fontface = label_fontface,
                                        show_other_cells = show_other_cells,
                                        other_cell_color = other_cell_color,
                                        other_point_size = other_point_size,
                                        background_color = background_color,
                                        vor_border_color = vor_border_color,
                                        vor_max_radius = vor_max_radius,
                                        vor_alpha = vor_alpha,
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
    if(point_shape %in% c('border', 'voronoi')) {
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
#' @param spat_unit spatial unit
#' @param feat_type feature type
#' @param show_image show a tissue background image
#' @param gimage a giotto image
#' @param image_name name of a giotto image or multiple images with group_by
#' @param largeImage_name name of a giottoLargeImage or multiple images with group_by
#' @param group_by create multiple plots based on cell annotation column
#' @param group_by_subset subset the group_by factor column
#' @param spat_loc_name name of spatial locations
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
#' @param point_shape shape of points (border, no_border or voronoi)
#' @param point_size size of point (cell)
#' @param point_alpha transparancy of point
#' @param point_border_col color of border around points
#' @param point_border_stroke stroke size of border around points
#' @param show_cluster_center plot center of selected clusters
#' @param show_center_label plot label of selected clusters
#' @param center_point_size size of center points
#' @param center_point_border_col border color of center points
#' @param center_point_border_stroke border stroke size of center points
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
#' @param coord_fix_ratio fix ratio between x and y-axis (default = 1)
#' @param title title of plot
#' @param show_legend show legend
#' @param legend_text size of legend text
#' @param legend_symbol_size size of legend symbols
#' @param background_color color of plot background
#' @param vor_border_color border colorr for voronoi plot
#' @param vor_max_radius maximum radius for voronoi 'cells'
#' @param vor_alpha transparancy of voronoi 'cells'
#' @param axis_text size of axis text
#' @param axis_title size of axis title
#' @param cow_n_col cowplot param: how many columns
#' @param cow_rel_h cowplot param: relative height
#' @param cow_rel_w cowplot param: relative width
#' @param cow_align cowplot param: how to align
#' @param show_plot show plot
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters, see \code{\link{showSaveParameters}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return ggplot
#' @details Description of parameters.
#' coord_fix_ratio: set to NULL to use default ggplot parameters
#' @family spatial visualizations
#' @export
#' @seealso \code{\link{spatPlot3D}}
spatPlot2D = function(gobject,
                      spat_unit = NULL,
                      feat_type = NULL,
                      show_image = F,
                      gimage = NULL,
                      image_name = NULL,
                      largeImage_name = NULL,
                      group_by = NULL,
                      group_by_subset = NULL,
                      spat_loc_name = NULL,
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
                      point_shape = c('border', 'no_border', 'voronoi'),
                      point_size = 3,
                      point_alpha = 1,
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
                      coord_fix_ratio = 1,
                      title = NULL,
                      show_legend = T,
                      legend_text = 10,
                      legend_symbol_size = 2,
                      background_color = 'white',
                      vor_border_color = 'white',
                      vor_max_radius = 200,
                      vor_alpha = 1,
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
                      spat_unit = spat_unit,
                      feat_type = feat_type,
                      show_image = show_image,
                      gimage = gimage,
                      image_name = image_name,
                      largeImage_name = largeImage_name,
                      spat_loc_name = spat_loc_name,
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
                      point_alpha = point_alpha,
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
                      vor_border_color = vor_border_color,
                      vor_max_radius = vor_max_radius,
                      vor_alpha = vor_alpha,
                      axis_text = axis_text,
                      axis_title = axis_title,
                      title = title,
                      show_plot = show_plot,
                      return_plot = return_plot,
                      save_plot = save_plot,
                      save_param =  save_param,
                      default_save_name = default_save_name)


  } else {

    # Set feat_type and spat_unit
    spat_unit = set_default_spat_unit(gobject = gobject,
                                      spat_unit = spat_unit)
    feat_type = set_default_feat_type(gobject = gobject,
                                      spat_unit = spat_unit,
                                      feat_type = feat_type)

    ## metadata
    comb_metadata = combineMetadata(gobject = gobject,
                                    spat_loc_name = spat_loc_name,
                                    feat_type = feat_type,
                                    spat_unit = spat_unit,
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
            cell_color_code = getDistinctColors(n = number_colors)
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
      temp_gobject = subsetGiotto(gobject = gobject,
                                  spat_unit = spat_unit,
                                  feat_type = feat_type,
                                  cell_ids = subset_cell_IDs,
                                  verbose = FALSE)


      if(length(unique_groups) == length(image_name)) {
        spec_image_name = image_name[group_id]
      } else {
        spec_image_name = image_name
      }

      pl = spatPlot2D_single(gobject = temp_gobject,
                             spat_unit = spat_unit,
                             feat_type = feat_type,
                             show_image = show_image,
                             gimage = gimage,
                             image_name = spec_image_name,
                             spat_loc_name = spat_loc_name,
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
                             point_alpha = point_alpha,
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
                             vor_border_color = vor_border_color,
                             vor_max_radius = vor_max_radius,
                             vor_alpha = vor_alpha,
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
#' @inheritDotParams spatPlot2D
#' @return ggplot
#' @details Description of parameters.
#' @family spatial visualizations
#' @export
#' @seealso \code{\link{spatPlot3D}}
spatPlot = function(...) {

  spatPlot2D(...)

}








## ** spatial deconvolution plotting ####


#' @title spatDeconvPlot
#' @name spatDeconvPlot
#' @description Visualize cell type enrichment / deconvolution results in a scatterpie
#' @param gobject giotto object
#' @param spat_unit spatial unit
#' @param feat_type feature type
#' @param deconv_name name of deconvolution results to use
#' @param show_image show a tissue background image
#' @param gimage a giotto image
#' @param image_name name of a giotto image
#' @param largeImage_name name of a giottoLargeImage
#' @param spat_loc_name name of spatial locations
#' @param sdimx x-axis dimension name (default = 'sdimx')
#' @param sdimy y-axis dimension name (default = 'sdimy')
#' @param cell_color_code named vector with colors
#' @param line_color color of line within pie charts
#' @param radius radios of pie charts
#' @param alpha alpha of pie charts
#' @param coord_fix_ratio fix ratio between x and y-axis
#' @param title title of plot
#' @param legend_text size of legend text
#' @param background_color color of plot background
#' @param title title for plot (default = deconv_name)
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
spatDeconvPlot = function(gobject,
                          spat_unit = NULL,
                          feat_type = NULL,
                          deconv_name = 'DWLS',
                          show_image = F,
                          gimage = NULL,
                          image_name = NULL,
                          largeImage_name = NULL,
                          spat_loc_name = NULL,
                          sdimx = 'sdimx',
                          sdimy = 'sdimy',
                          cell_color_code = NULL,
                          line_color = NA,
                          radius = 10,
                          alpha = 1,
                          legend_text = 8,
                          background_color = 'white',
                          title = NULL,
                          axis_text = 8,
                          axis_title = 8,
                          coord_fix_ratio = TRUE,
                          show_plot = NA,
                          return_plot = NA,
                          save_plot = NA,
                          save_param =  list(),
                          default_save_name = 'spatDeconvPlot') {


  # check for installed packages
  package_check(pkg_name = "scatterpie", repository = "CRAN")

  # Set feat_type and spat_unit
  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = spat_unit)
  feat_type = set_default_feat_type(gobject = gobject,
                                    spat_unit = spat_unit,
                                    feat_type = feat_type)

  ## giotto image ##
  if(show_image == TRUE) {

    if(!is.null(gimage)) {
      gimage = gimage
    } else if(!is.null(image_name)) {
      # If there is input to image_name arg

      if(length(image_name) == 1) {
        gimage = gobject@images[[image_name]]
        if(is.null(gimage)) warning('image_name: ', image_name, ' does not exist \n')
      } else {
        gimage = list()
        for(gim in 1:length(image_name)) {
          gimage[[gim]] = gobject@images[[gim]]
          if(is.null(gimage[[gim]])) warning('image_name: ', gim, ' does not exists \n')
        }
      }

    } else if(!is.null(largeImage_name)) {
      # If there is input to largeImage_name arg

      if(length(largeImage_name) == 1) {
        gimage = plot_auto_largeImage_resample(gobject = gobject,
                                               largeImage_name = largeImage_name,
                                               spat_unit = spat_unit,
                                               spat_loc_name = spat_loc_name,
                                               include_image_in_border = TRUE)
      } else {
        gimage = list()
        for(gim in 1:length(largeImage_name)) {
          gimage[[gim]] = plot_auto_largeImage_resample(gobject = gobject,
                                                        largeImage_name = largeImage_name[[gim]],
                                                        spat_unit = spat_unit,
                                                        spat_loc_name = spat_loc_name,
                                                        include_image_in_border = TRUE)
        }
      }

    } else {
      # Default to first image available in images if no input given to image_name or largeImage_name args
      image_name = names(gobject@images)[1]
      gimage = gobject@images[[image_name]]
      if(is.null(gimage)) warning('image_name: ', image_name, ' does not exist \n')
    }
  }


  ## get spatial cell locations
  spatial_locations = get_spatial_locations(gobject = gobject,
                                            spat_unit = spat_unit,
                                            spat_loc_name = spat_loc_name)

  ## deconvolution results
  spatial_enrichment = get_spatial_enrichment(gobject = gobject,
                                              spat_unit = spat_unit,
                                              feat_type = feat_type,
                                              enrichm_name = deconv_name)




  ### create 2D plot with ggplot ###
  #cat('create 2D plot with ggplot \n')


  pl <- ggplot2::ggplot()
  pl <- pl + ggplot2::theme_bw()

  ## plot image ##
  if(show_image == TRUE & !is.null(gimage)) {
    pl = plot_spat_image_layer_ggplot(gg_obj = pl,
                                      gobject = gobject,
                                      feat_type = feat_type,
                                      spat_unit = spat_unit,
                                      spat_loc_name = spat_loc_name,
                                      gimage = gimage,
                                      sdimx = sdimx,
                                      sdimy = sdimy)
  }


  ## plot scatterpie ##
  pl = plot_spat_scatterpie_layer_ggplot(ggobject = pl,
                                         sdimx = sdimx,
                                         sdimy = sdimy,
                                         spatial_locations = spatial_locations,
                                         spatial_enrichment = spatial_enrichment,
                                         radius = radius,
                                         color = line_color,
                                         alpha = alpha,
                                         cell_color_code = cell_color_code)


  ## adjust theme settings
  pl = pl + ggplot2::theme(plot.title = element_text(hjust = 0.5),
                           legend.title = element_blank(),
                           legend.text = element_text(size = legend_text),
                           axis.title = element_text(size = axis_title),
                           axis.text = element_text(size = axis_text),
                           panel.grid = element_blank(),
                           panel.background = element_rect(fill = background_color))

  # fix coord ratio
  if(!is.null(coord_fix_ratio)) {
    pl <- pl + ggplot2::coord_fixed(ratio = coord_fix_ratio)
  }

  # provide x, y and plot titles
  if(is.null(title)) title = deconv_name
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
    do.call('all_plots_save_function', c(list(gobject = gobject,
                                              plot_object = pl,
                                              default_save_name = default_save_name), save_param))
  }

  ## return plot
  if(return_plot == TRUE) {
    return(pl)
  }


}




## ** spatial and dim reduction plotting ####


#' @title spatDimPlot2D
#' @name spatDimPlot2D
#' @description Visualize cells according to spatial AND dimension reduction coordinates 2D
#' @param gobject giotto object
#' @param spat_unit spatial unit
#' @param feat_type feature type
#' @param show_image show a tissue background image
#' @param gimage a giotto image
#' @param image_name name of a giotto image
#' @param largeImage_name name of a giottoLargeImage
#' @param spat_loc_name name of spatial locations
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
#' @param dim_point_alpha transparancy of point in dim. reduction space
#' @param dim_point_border_col border color of points in dim. reduction space
#' @param dim_point_border_stroke border stroke of points in dim. reduction space
#' @param spat_point_shape shape of points (border, no_border or voronoi)
#' @param spat_point_size size of spatial points
#' @param spat_point_alpha transparancy of spatial points
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
#' @param spat_center_point_border_col border color of spatial center points
#' @param spat_center_point_border_stroke border strike size of spatial center points
#' @param spat_label_size size of the center label
#' @param spat_label_fontface font of the center label
#' @param show_NN_network show underlying NN network
#' @param nn_network_to_use type of NN network to use (kNN vs sNN)
#' @param network_name name of NN network to use, if show_NN_network = TRUE
#' @param nn_network_alpha column to use for alpha of the edges
#' @param show_spatial_network show spatial network
#' @param spat_network_name name of spatial network to use
#' @param spat_network_color color of spatial network
#' @param spat_network_alpha alpha of spatial network
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
#' @param vor_border_color border colorr for voronoi plot
#' @param vor_max_radius maximum radius for voronoi 'cells'
#' @param vor_alpha transparancy of voronoi 'cells'
#' @param axis_text size of axis text
#' @param axis_title size of axis title
#' @param show_plot show plot
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters, see \code{\link{showSaveParameters}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return ggplot
#' @details Description of parameters.
#' @family spatial and dimension reduction visualizations
#' @export
#' @seealso \code{\link{spatDimPlot3D}}
spatDimPlot2D <- function(gobject,
                          spat_unit = NULL,
                          feat_type = NULL,
                          show_image = F,
                          gimage = NULL,
                          image_name = NULL,
                          largeImage_name = NULL,
                          spat_loc_name = NULL,
                          plot_alignment = c('vertical', 'horizontal'),
                          dim_reduction_to_use = 'umap',
                          dim_reduction_name = NULL,
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
                          dim_point_alpha = 1,
                          dim_point_border_col = 'black',
                          dim_point_border_stroke = 0.1,
                          spat_point_shape = c('border', 'no_border', 'voronoi'),
                          spat_point_size = 1,
                          spat_point_alpha = 1,
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
                          spat_center_point_border_col = 'blue',
                          spat_center_point_border_stroke = 0.1,
                          spat_label_size = 4,
                          spat_label_fontface = 'bold',
                          show_NN_network = F,
                          nn_network_to_use = 'sNN',
                          network_name = 'sNN.pca',
                          nn_network_alpha = 0.05,
                          show_spatial_network = F,
                          spat_network_name = 'Delaunay_network',
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
                          legend_text = 10,
                          legend_symbol_size = 2,
                          dim_background_color = 'white',
                          spat_background_color = 'white',
                          vor_border_color = 'white',
                          vor_max_radius = 200,
                          vor_alpha = 1,
                          axis_text = 8,
                          axis_title = 8,
                          show_plot = NA,
                          return_plot = NA,
                          save_plot = NA,
                          save_param =  list(),
                          default_save_name = 'spatDimPlot2D'){

  # Set feat_type and spat_unit
  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = spat_unit)
  feat_type = set_default_feat_type(gobject = gobject,
                                    spat_unit = spat_unit,
                                    feat_type = feat_type)

  plot_alignment = match.arg(plot_alignment, choices = c( 'vertical','horizontal'))


  # create matching cell_color_code
  if(is.null(cell_color_code)) {
    if(is.character(cell_color)) {

      cell_metadata = pDataDT(gobject,
                              spat_unit = spat_unit,
                              feat_type = feat_type)
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
  dmpl = dimPlot2D(gobject = gobject,
                   spat_unit = spat_unit,
                   feat_type = feat_type,
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
                   point_alpha = dim_point_alpha,
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
                   spat_unit = spat_unit,
                   feat_type = feat_type,
                   show_image = show_image,
                   gimage = gimage,
                   image_name = image_name,
                   largeImage_name = largeImage_name,
                   spat_loc_name = spat_loc_name,
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
                   point_alpha = spat_point_alpha,
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
                   vor_border_color = vor_border_color,
                   vor_max_radius = vor_max_radius,
                   vor_alpha = vor_alpha,
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
#' @inheritDotParams spatDimPlot2D
#' @return ggplot
#' @details Description of parameters.
#' @family spatial and dimension reduction visualizations
#' @export
#' @seealso \code{\link{spatDimPlot2D}} and \code{\link{spatDimPlot3D}} for 3D visualization.
spatDimPlot = function(...) {

  spatDimPlot2D(...)

}



## ** spatial feature plotting ####

#' @title spatFeatPlot2D_single
#' @name spatFeatPlot2D_single
#' @description Visualize cells and feature expression according to spatial coordinates
#' @param gobject giotto object
#' @param spat_unit spatial unit
#' @param feat_type feature type
#' @param show_image show a tissue background image
#' @param gimage a giotto image
#' @param image_name name of a giotto image
#' @param largeImage_name name of a giottoLargeImage
#' @param spat_loc_name name of spatial locations
#' @param sdimx x-axis dimension name (default = 'sdimx')
#' @param sdimy y-axis dimension name (default = 'sdimy')
#' @param expression_values gene expression values to use
#' @param feats features to show
#' @param cell_color_gradient vector with 3 colors for numeric data
#' @param gradient_midpoint midpoint for color gradient
#' @param gradient_limits vector with lower and upper limits
#' @param show_network show underlying spatial network
#' @param network_color color of spatial network
#' @param spatial_network_name name of spatial network to use
#' @param edge_alpha alpha of edge
#' @param show_grid show spatial grid
#' @param grid_color color of spatial grid
#' @param spatial_grid_name name of spatial grid to use
#' @param midpoint expression midpoint
#' @param scale_alpha_with_expression scale expression with ggplot alpha parameter
#' @param point_shape shape of points (border, no_border or voronoi)
#' @param point_size size of point (cell)
#' @param point_alpha transparancy of points
#' @param point_border_col color of border around points
#' @param point_border_stroke stroke size of border around points
#' @param cow_n_col cowplot param: how many columns
#' @param cow_rel_h cowplot param: relative height
#' @param cow_rel_w cowplot param: relative width
#' @param cow_align cowplot param: how to align
#' @param show_legend show legend
#' @param legend_text size of legend text
#' @param background_color color of plot background
#' @param vor_border_color border colorr for voronoi plot
#' @param vor_max_radius maximum radius for voronoi 'cells'
#' @param vor_alpha transparency of voronoi 'cells'
#' @param axis_text size of axis text
#' @param axis_title size of axis title
#' @param show_plot show plots
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters, see \code{\link{showSaveParameters}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return ggplot
#' @details Description of parameters.
#' @family spatial feature expression visualizations
#' @export
#' @seealso \code{\link{spatGenePlot3D}}
spatFeatPlot2D_single <- function(gobject,
                                  spat_unit = NULL,
                                  feat_type = NULL,
                                  show_image = F,
                                  gimage = NULL,
                                  image_name = NULL,
                                  largeImage_name = NULL,
                                  spat_loc_name = NULL,
                                  sdimx = 'sdimx',
                                  sdimy = 'sdimy',
                                  expression_values = c('normalized', 'scaled', 'custom'),
                                  feats,
                                  cell_color_gradient = c('blue', 'white', 'red'),
                                  gradient_midpoint = NULL,
                                  gradient_limits = NULL,
                                  show_network = F,
                                  network_color = NULL,
                                  spatial_network_name = 'Delaunay_network',
                                  edge_alpha = NULL,
                                  show_grid = F,
                                  grid_color = NULL,
                                  spatial_grid_name = 'spatial_grid',
                                  midpoint = 0,
                                  scale_alpha_with_expression = FALSE,
                                  point_shape = c('border', 'no_border', 'voronoi'),
                                  point_size = 1,
                                  point_alpha = 1,
                                  point_border_col = 'black',
                                  point_border_stroke = 0.1,
                                  show_legend = T,
                                  legend_text = 8,
                                  background_color = 'white',
                                  vor_border_color = 'white',
                                  vor_alpha = 1,
                                  vor_max_radius = 200,
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
                                  default_save_name = 'spatFeatPlot2D_single') {


  # data.table variables
  cell_ID = NULL

  # print, return and save parameters
  show_plot = ifelse(is.na(show_plot), readGiottoInstructions(gobject, param = 'show_plot'), show_plot)
  save_plot = ifelse(is.na(save_plot), readGiottoInstructions(gobject, param = 'save_plot'), save_plot)
  return_plot = ifelse(is.na(return_plot), readGiottoInstructions(gobject, param = 'return_plot'), return_plot)

  # Set feat_type and spat_unit
  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = spat_unit)
  feat_type = set_default_feat_type(gobject = gobject,
                                    spat_unit = spat_unit,
                                    feat_type = feat_type)

  ## giotto image ##
  if(show_image == TRUE) {

    if(!is.null(gimage)) {
      gimage = gimage
    } else if(!is.null(image_name)) {
      # if there is input to image_name arg
      gimage = gobject@images[[image_name]]
      if(is.null(gimage)) warning('image_name: ', image_name, ' does not exist \n')
    } else if (!is.null(largeImage_name)) {
      # if there is input to largeImage_name arg

      gimage = plot_auto_largeImage_resample(gobject = gobject,
                                             largeImage_name = largeImage_name,
                                             spat_unit = spat_unit,
                                             spat_loc_name = spat_loc_name,
                                             include_image_in_border = TRUE)
    } else {
      # Default to first image available in images if no input given to image_name or largeImage_name args
      image_name = names(gobject@images)[1]
      gimage = gobject@images[[image_name]]
      if(is.null(gimage)) warning('image_name: ', image_name, ' does not exist \n')
    }
  }

  # point shape
  point_shape = match.arg(point_shape, choices = c('border', 'no_border', 'voronoi'))

  # expression values
  values = match.arg(expression_values, unique(c('normalized', 'scaled', 'custom', expression_values)))
  expr_values = get_expression_values(gobject = gobject,
                                      spat_unit = spat_unit,
                                      feat_type = feat_type,
                                      values = values)

  # only keep feats that are in the dataset
  selected_feats = feats
  selected_feats = selected_feats[selected_feats %in% rownames(expr_values) ]


  # get selected feat expression values in data.table format
  if(length(selected_feats) == 1) {
    subset_expr_data = expr_values[rownames(expr_values) %in% selected_feats, ]
    t_sub_expr_data_DT = data.table::data.table('selected_feat' = subset_expr_data, 'cell_ID' = colnames(expr_values))
    data.table::setnames(t_sub_expr_data_DT, 'selected_feat', selected_feats)
  } else {
    subset_expr_data = expr_values[rownames(expr_values) %in% selected_feats, ]
    t_sub_expr_data = t_flex(subset_expr_data)
    t_sub_expr_data_DT = data.table::as.data.table(as.matrix(t_sub_expr_data))
    t_sub_expr_data_DT[, cell_ID := rownames(t_sub_expr_data)]
  }


  ## extract cell locations
  if(is.null(spat_loc_name)) {
    if(!is.null(gobject@spatial_locs)) {
      spat_loc_name = names(gobject@spatial_locs[[spat_unit]])[[1]]
      # cat('No spatial locations have been selected, the first one -',spat_loc_name, '- will be used \n')
    } else {
      spat_loc_name = NULL
      cat('No spatial locations have been found \n')
      return(NULL)
    }
  }

  cell_locations  = get_spatial_locations(gobject,
                                          spat_unit = spat_unit,
                                          spat_loc_name = spat_loc_name)

  ## extract spatial network
  if(show_network == TRUE) {
    spatial_network = get_spatialNetwork(gobject,
                                         spat_unit = spat_unit,
                                         name = spatial_network_name,
                                         return_network_Obj = FALSE)
  } else {
    spatial_network = NULL
  }

  ## extract spatial grid
  if(show_grid == TRUE) {
    spatial_grid = get_spatialGrid(gobject,
                                   spat_unit = spat_unit,
                                   feat_type = feat_type,
                                   name = spatial_grid_name)

  } else {
    spatial_grid = NULL
  }

  ## extract cell metadata
  cell_metadata = combineMetadata(gobject = gobject,
                                  spat_unit = spat_unit,
                                  spat_loc_name = spat_loc_name,
                                  feat_type = feat_type)

  if(nrow(cell_metadata) == 0) {
    cell_locations_metadata = cell_locations
  } else {
    cell_locations_metadata = cell_metadata
  }

  cell_locations_metadata_feats <- merge(cell_locations_metadata,
                                         t_sub_expr_data_DT,
                                         by = 'cell_ID')


  ## plotting ##
  savelist <- list()

  for(feat in selected_feats) {

    pl <- ggplot2::ggplot()
    pl <- pl + ggplot2::theme_classic()


    ## plot image ## TODO
    ## plot image ##
    if(show_image == TRUE & !is.null(gimage)) {
      pl = plot_spat_image_layer_ggplot(gg_obj = pl,
                                        gobject = gobject,
                                        feat_type = feat_type,
                                        spat_unit = spat_unit,
                                        spat_loc_name = spat_loc_name,
                                        gimage = gimage,
                                        sdimx = sdimx,
                                        sdimy = sdimy)
    }

    ## plot network or grid first if point_shape is border or no_border point
    if(point_shape %in% c('border', 'no_border')) {

      ## plot spatial network
      if(!is.null(spatial_network) & show_network == TRUE) {
        if(is.null(network_color)) {
          network_color = 'red'
        }
        xbegin = paste0(sdimx, '_begin')
        ybegin = paste0(sdimy, '_begin')
        xend = paste0(sdimx, '_end')
        yend = paste0(sdimy, '_end')
        pl <- pl + ggplot2::geom_segment(data = spatial_network, aes_string(x = xbegin,
                                                                            y = ybegin,
                                                                            xend = xend,
                                                                            yend = yend),
                                         color = network_color, size = 0.5, alpha = 0.5)
      }

      ## plot spatial grid
      if(!is.null(spatial_grid) & show_grid == TRUE) {
        if(is.null(grid_color)) grid_color = 'black'

        xmin = paste0(gsub(pattern = 'sdim', replacement = '', x = sdimx), '_start')
        ymin = paste0(gsub(pattern = 'sdim', replacement = '', x = sdimy), '_start')
        xmax = paste0(gsub(pattern = 'sdim', replacement = '', x = sdimx), '_end')
        ymax = paste0(gsub(pattern = 'sdim', replacement = '', x = sdimy), '_end')

        pl <- pl + ggplot2::geom_rect(data = spatial_grid, aes_string(xmin = xmin,
                                                                      xmax = xmax,
                                                                      ymin = ymin,
                                                                      ymax = ymax),
                                      color = grid_color, fill = NA)
      }

    }



    ### plot cells ###

    ## set gradient limits if needed ##
    if(!is.null(gradient_limits) & is.vector(gradient_limits) & length(gradient_limits) == 2) {
      lower_lim = gradient_limits[[1]]
      upper_lim = gradient_limits[[2]]
      numeric_data = cell_locations_metadata_feats[[feat]]
      limit_numeric_data = ifelse(numeric_data > upper_lim, upper_lim,
                                  ifelse(numeric_data < lower_lim, lower_lim, numeric_data))
      cell_locations_metadata_feats[[feat]] = limit_numeric_data
    }

    if(is.null(gradient_midpoint)) {
      gradient_midpoint = stats::median(cell_locations_metadata_feats[[feat]])
    }


    ## with border ##
    if(point_shape == 'border') {

      if(scale_alpha_with_expression == TRUE) {
        pl <- pl + ggplot2::geom_point(data = cell_locations_metadata_feats,
                                       aes_string2(x = sdimx,
                                                   y = sdimy,
                                                   fill = feat,
                                                   alpha = feat),
                                       shape = 21,
                                       color = point_border_col, size = point_size,
                                       stroke = point_border_stroke,
                                       show.legend = show_legend)
      } else {
        pl <- pl + ggplot2::geom_point(data = cell_locations_metadata_feats,
                                       aes_string2(x = sdimx,
                                                   y = sdimy,
                                                   fill = feat),
                                       shape = 21,
                                       color = point_border_col, size = point_size,
                                       stroke = point_border_stroke,
                                       show.legend = show_legend, alpha = point_alpha)
      }


      ## scale and labs ##
      pl <- pl + ggplot2::scale_alpha_continuous(guide = 'none')
      pl <- pl + ggplot2::scale_fill_gradient2(low = cell_color_gradient[[1]],
                                               mid = cell_color_gradient[[2]],
                                               high = cell_color_gradient[[3]],
                                               midpoint = gradient_midpoint,
                                               guide = guide_colorbar(title = ''))
      pl <- pl + ggplot2::labs(x = 'coord x', y = 'coord y', title = feat)


    }



    ## no border ##
    if(point_shape == 'no_border') {

      if(scale_alpha_with_expression == TRUE) {
        pl <- pl + ggplot2::geom_point(data = cell_locations_metadata_feats,
                                       aes_string2(x = sdimx,
                                                   y = sdimy,
                                                   color = feat,
                                                   alpha = feat),
                                       shape = 19, size = point_size,
                                       show.legend = show_legend)
      } else {
        pl <- pl + ggplot2::geom_point(data = cell_locations_metadata_feats,
                                       aes_string2(x = sdimx,
                                                   y = sdimy,
                                                   color = feat),
                                       shape = 19, size = point_size,
                                       show.legend = show_legend, alpha = point_alpha)
      }


      ## scale and labs ##
      pl <- pl + ggplot2::scale_alpha_continuous(guide = 'none')
      pl <- pl + ggplot2::scale_color_gradient2(low = cell_color_gradient[[1]],
                                                mid = cell_color_gradient[[2]],
                                                high = cell_color_gradient[[3]],
                                                midpoint = gradient_midpoint,
                                                guide = guide_colorbar(title = ''))
      pl <- pl + ggplot2::labs(x = 'coord x', y = 'coord y', title = feat)

    }


    ## voronoi ##
    if(point_shape == 'voronoi') {

      if(scale_alpha_with_expression == TRUE) {
        pl = pl + ggforce::geom_voronoi_tile(data = cell_locations_metadata_feats,
                                             aes_string(x = sdimx, y = sdimy,
                                                        group = '-1L',
                                                        fill = feat, alpha = feat),
                                             colour = vor_border_color,
                                             max.radius = vor_max_radius,
                                             show.legend = show_legend)
      } else {
        pl = pl + ggforce::geom_voronoi_tile(data = cell_locations_metadata_feats,
                                             aes_string(x = sdimx, y = sdimy,
                                                        group = '-1L',
                                                        fill = feat),
                                             colour = vor_border_color,
                                             max.radius = vor_max_radius,
                                             show.legend = show_legend,
                                             alpha = vor_alpha)
      }


      ## plot spatial network
      if(!is.null(spatial_network) & show_network == TRUE) {
        if(is.null(network_color)) {
          network_color = 'red'
        }
        xbegin = paste0(sdimx, '_begin')
        ybegin = paste0(sdimy, '_begin')
        xend = paste0(sdimx, '_end')
        yend = paste0(sdimy, '_end')
        pl <- pl + ggplot2::geom_segment(data = spatial_network, aes_string(x = xbegin, y = ybegin,
                                                                            xend = xend, yend = yend),
                                         color = network_color, size = 0.5, alpha = 0.5)
      }

      ## plot spatial grid
      if(!is.null(spatial_grid) & show_grid == TRUE) {
        if(is.null(grid_color)) grid_color = 'black'

        xmin = paste0(gsub(pattern = 'sdim', replacement = '', x = sdimx), '_start')
        ymin = paste0(gsub(pattern = 'sdim', replacement = '', x = sdimy), '_start')
        xmax = paste0(gsub(pattern = 'sdim', replacement = '', x = sdimx), '_end')
        ymax = paste0(gsub(pattern = 'sdim', replacement = '', x = sdimy), '_end')

        pl <- pl + ggplot2::geom_rect(data = spatial_grid, aes_string(xmin = xmin, xmax = xmax,
                                                                      ymin = ymin, ymax = ymax),
                                      color = grid_color, fill = NA)
      }


      ## scale and labs ##
      pl <- pl + ggplot2::scale_alpha_continuous(guide = 'none')
      pl <- pl + ggplot2::scale_fill_gradient2(low = cell_color_gradient[[1]],
                                               mid = cell_color_gradient[[2]],
                                               high = cell_color_gradient[[3]],
                                               midpoint = gradient_midpoint,
                                               guide = guide_colorbar(title = ''))
      pl <- pl + ggplot2::labs(x = 'coord x', y = 'coord y', title = feat)


    }

    ## theme ##
    pl <- pl + ggplot2::theme(plot.title = element_text(hjust = 0.5),
                              legend.title = element_blank(),
                              legend.text = element_text(size = legend_text),
                              axis.title = element_text(size = axis_title),
                              axis.text = element_text(size = axis_text),
                              panel.grid = element_blank(),
                              panel.background = element_rect(fill = background_color))


    savelist[[feat]] <- pl
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


#' @title Plot data in physical space 2D
#' @name spatFeatPlot2D
#' @description Visualize cells and feature expression according to spatial coordinates
#' @param gobject giotto object
#' @param feat_type feature type
#' @param spat_unit spatial unit
#' @param show_image show a tissue background image
#' @param gimage a giotto image
#' @param image_name name of a giotto image or multiple images if group_by
#' @param largeImage_name name of a giottoLargeImage or multiple images if group by
#' @param spat_loc_name name of spatial locations
#' @param group_by create multiple plots based on cell annotation column
#' @param group_by_subset subset the group_by factor column
#' @param sdimx x-axis dimension name (default = 'sdimx')
#' @param sdimy y-axis dimension name (default = 'sdimy')
#' @param expression_values gene expression values to use
#' @param feats features to show
#' @param cell_color_gradient vector with 3 colors for numeric data
#' @param gradient_midpoint midpoint for color gradient
#' @param gradient_limits vector with lower and upper limits
#' @param show_network show underlying spatial network
#' @param network_color color of spatial network
#' @param spatial_network_name name of spatial network to use
#' @param edge_alpha alpha of edge
#' @param show_grid show spatial grid
#' @param grid_color color of spatial grid
#' @param spatial_grid_name name of spatial grid to use
#' @param midpoint expression midpoint
#' @param scale_alpha_with_expression scale expression with ggplot alpha parameter
#' @param point_shape shape of points (border, no_border or voronoi)
#' @param point_size size of point (cell)
#' @param point_alpha transparancy of points
#' @param point_border_col color of border around points
#' @param point_border_stroke stroke size of border around points
#' @param cow_n_col cowplot param: how many columns
#' @param cow_rel_h cowplot param: relative height
#' @param cow_rel_w cowplot param: relative width
#' @param cow_align cowplot param: how to align
#' @param show_legend show legend
#' @param legend_text size of legend text
#' @param background_color color of plot background
#' @param vor_border_color border colorr for voronoi plot
#' @param vor_max_radius maximum radius for voronoi 'cells'
#' @param vor_alpha transparancy of voronoi 'cells'
#' @param axis_text size of axis text
#' @param axis_title size of axis title
#' @param show_plot show plots
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters, see \code{\link{showSaveParameters}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return ggplot
#' @details Description of parameters.
#' @family spatial feature expression visualizations
#' @export
#' @seealso \code{\link{spatGenePlot3D}}
spatFeatPlot2D <- function(gobject,
                           feat_type = NULL,
                           spat_unit = NULL,
                           show_image = F,
                           gimage = NULL,
                           image_name = NULL,
                           largeImage_name = NULL,
                           spat_loc_name = NULL,
                           group_by = NULL,
                           group_by_subset = NULL,
                           sdimx = 'sdimx',
                           sdimy = 'sdimy',
                           expression_values = c('normalized', 'scaled', 'custom'),
                           feats,
                           cell_color_gradient = c('blue', 'white', 'red'),
                           gradient_midpoint = NULL,
                           gradient_limits = NULL,
                           show_network = F,
                           network_color = NULL,
                           spatial_network_name = 'Delaunay_network',
                           edge_alpha = NULL,
                           show_grid = F,
                           grid_color = NULL,
                           spatial_grid_name = 'spatial_grid',
                           midpoint = 0,
                           scale_alpha_with_expression = FALSE,
                           point_shape = c('border', 'no_border', 'voronoi'),
                           point_size = 1,
                           point_alpha = 1,
                           point_border_col = 'black',
                           point_border_stroke = 0.1,
                           show_legend = T,
                           legend_text = 8,
                           background_color = 'white',
                           vor_border_color = 'white',
                           vor_alpha = 1,
                           vor_max_radius = 200,
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
                           default_save_name = 'spatFeatPlot2D') {

  ## check group_by
  if(is.null(group_by)) {

    spatFeatPlot2D_single(gobject = gobject,
                          feat_type = feat_type,
                          spat_unit = spat_unit,
                          show_image = show_image,
                          gimage = gimage,
                          spat_loc_name = spat_loc_name,
                          image_name = image_name,
                          largeImage_name = largeImage_name,
                          sdimx = sdimx,
                          sdimy = sdimy,
                          expression_values = expression_values,
                          feats = feats,
                          cell_color_gradient = cell_color_gradient,
                          gradient_midpoint = gradient_midpoint,
                          gradient_limits = gradient_limits,
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
                          point_alpha = point_alpha,
                          point_border_col = point_border_col,
                          point_border_stroke = point_border_stroke,
                          show_legend = show_legend,
                          legend_text = legend_text,
                          background_color = background_color,
                          vor_border_color = vor_border_color,
                          vor_alpha = vor_alpha,
                          vor_max_radius = vor_max_radius,
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


  } else {

    # Set feat_type and spat_unit
    spat_unit = set_default_spat_unit(gobject = gobject,
                                      spat_unit = spat_unit)
    feat_type = set_default_feat_type(gobject = gobject,
                                      feat_type = feat_type,
                                      spat_unit = spat_unit)

    ## metadata
    comb_metadata = combineMetadata(gobject = gobject,
                                    spat_unit = spat_unit,
                                    spat_loc_name = spat_loc_name,
                                    feat_type = feat_type)
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


    # print, return and save parameters
    show_plot = ifelse(is.na(show_plot), readGiottoInstructions(gobject, param = 'show_plot'), show_plot)
    save_plot = ifelse(is.na(save_plot), readGiottoInstructions(gobject, param = 'save_plot'), save_plot)
    return_plot = ifelse(is.na(return_plot), readGiottoInstructions(gobject, param = 'return_plot'), return_plot)

    ## plotting ##
    savelist <- list()


    for(group_id in 1:length(unique_groups)) {

      group = unique_groups[group_id]

      subset_cell_IDs = comb_metadata[get(group_by) == group][['cell_ID']]
      temp_gobject = subsetGiotto(gobject = gobject,
                                  feat_type = feat_type,
                                  spat_unit = spat_unit,
                                  cell_ids = subset_cell_IDs)


      if(length(unique_groups) == length(image_name)) {
        spec_image_name = image_name[group_id]
      } else {
        spec_image_name = image_name
      }

      pl = spatFeatPlot2D_single(gobject = temp_gobject,
                                 feat_type = feat_type,
                                 spat_unit = spat_unit,
                                 show_image = show_image,
                                 gimage = gimage,
                                 image_name = spec_image_name,
                                 largeImage_name = largeImage_name,
                                 spat_loc_name = spat_loc_name,
                                 sdimx = sdimx,
                                 sdimy = sdimy,
                                 expression_values = expression_values,
                                 feats = feats,
                                 cell_color_gradient = cell_color_gradient,
                                 gradient_midpoint = gradient_midpoint,
                                 gradient_limits = gradient_limits,
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
                                 point_alpha = point_alpha,
                                 point_border_col = point_border_col,
                                 point_border_stroke = point_border_stroke,
                                 show_legend = show_legend,
                                 legend_text = legend_text,
                                 background_color = background_color,
                                 vor_border_color = vor_border_color,
                                 vor_alpha = vor_alpha,
                                 vor_max_radius = vor_max_radius,
                                 axis_text = axis_text,
                                 axis_title = axis_title,
                                 cow_n_col = 1,
                                 cow_rel_h = cow_rel_h,
                                 cow_rel_w = cow_rel_w,
                                 cow_align = cow_align,
                                 show_plot = FALSE,
                                 return_plot = TRUE,
                                 save_plot = FALSE,
                                 save_param =  save_param,
                                 default_save_name = 'spatFeatPlot2D')


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



#' @title spatGenePlot2D
#' @name spatGenePlot2D
#' @description Visualize cells and gene expression according to spatial coordinates
#' @param gobject giotto object
#' @param genes genes to show
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @inheritDotParams spatFeatPlot2D -gobject -feats -default_save_name
#' @return ggplot
#' @details Description of parameters, see \code{\link{spatFeatPlot2D}}
#' @family spatial gene expression visualizations
#' @export
spatGenePlot2D <- function(gobject,
                           genes,
                           default_save_name = 'spatGenePlot2D',
                           ...) {

  .Deprecated(new = "spatFeatPlot2D")

  spatFeatPlot2D(gobject = gobject,
                 feat_type = 'rna',
                 feats = genes,
                 default_save_name = default_save_name,
                 ...)

}


#' @title spatGenePlot
#' @name spatGenePlot
#' @description Visualize cells and gene expression according to spatial coordinates
#' @inheritDotParams spatGenePlot2D
#' @return ggplot
#' @details Description of parameters.
#' @family spatial gene expression visualizations
#' @export
#' @seealso \code{\link{spatGenePlot3D}} and \code{\link{spatGenePlot2D}}
spatGenePlot = function(...) {

  .Deprecated(new = "spatFeatPlot2D")

  spatGenePlot2D(...)

}



## ** dim reduction feature plotting ####

#' @title dimFeatPlot2D
#' @name dimFeatPlot2D
#' @description Visualize gene expression according to dimension reduction coordinates
#' @param gobject giotto object
#' @param feat_type  feature type
#' @param spat_unit spatial unit
#' @param expression_values gene expression values to use
#' @param feats features to show
#' @param dim_reduction_to_use dimension reduction to use
#' @param dim_reduction_name dimension reduction name
#' @param dim1_to_use dimension to use on x-axis
#' @param dim2_to_use dimension to use on y-axis
#' @param show_NN_network show underlying NN network
#' @param nn_network_to_use type of NN network to use (kNN vs sNN)
#' @param network_name name of NN network to use, if show_NN_network = TRUE
#' @param network_color color of NN network
#' @param edge_alpha column to use for alpha of the edges
#' @param scale_alpha_with_expression scale expression with ggplot alpha parameter
#' @param point_shape point with border or not (border or no_border)
#' @param point_size size of point (cell)
#' @param point_alpha transparancy of points
#' @param cell_color_gradient vector with 3 colors for numeric data
#' @param gradient_midpoint midpoint for color gradient
#' @param gradient_limits vector with lower and upper limits
#' @param point_border_col color of border around points
#' @param point_border_stroke stroke size of border around points
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
#' @param save_param list of saving parameters, see \code{\link{showSaveParameters}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return ggplot
#' @details Description of parameters.
#' @family dimension reduction feature expression visualizations
#' @export
dimFeatPlot2D <- function(gobject,
                          spat_unit = NULL,
                          feat_type = NULL,
                          expression_values = c('normalized', 'scaled', 'custom'),
                          feats = NULL,
                          dim_reduction_to_use = 'umap',
                          dim_reduction_name = NULL,
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
                          point_alpha = 1,
                          cell_color_gradient = c('blue', 'white', 'red'),
                          gradient_midpoint = NULL,
                          gradient_limits = NULL,
                          point_border_col = 'black',
                          point_border_stroke = 0.1,
                          show_legend = T,
                          legend_text = 10,
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
                          default_save_name = 'dimFeatPlot2D') {


  # print, return and save parameters
  show_plot = ifelse(is.na(show_plot), readGiottoInstructions(gobject, param = 'show_plot'), show_plot)
  save_plot = ifelse(is.na(save_plot), readGiottoInstructions(gobject, param = 'save_plot'), save_plot)
  return_plot = ifelse(is.na(return_plot), readGiottoInstructions(gobject, param = 'return_plot'), return_plot)

  # point shape
  point_shape = match.arg(point_shape, choices = c('border', 'no_border'))

  # Set feat_type and spat_unit
  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = spat_unit)
  feat_type = set_default_feat_type(gobject = gobject,
                                    spat_unit = spat_unit,
                                    feat_type = feat_type)

  # specify dim_reduction_name according to provided feat_type
  if(!is.null(dim_reduction_to_use)) {
    if(is.null(dim_reduction_name)) {
      if(feat_type == 'rna') {
        dim_reduction_name = dim_reduction_to_use
      } else {
        dim_reduction_name = paste0(feat_type,'.',dim_reduction_to_use)
      }
    }
  }


  # expression values
  values = match.arg(expression_values, unique(c('normalized', 'scaled', 'custom', expression_values)))
  expr_values = get_expression_values(gobject = gobject,
                                      spat_unit = spat_unit,
                                      feat_type = feat_type,
                                      values = values)

  # only keep feats that are in the dataset
  selected_feats = feats
  selected_feats = selected_feats[selected_feats %in% rownames(expr_values) ]

  #
  if(length(selected_feats) == 1) {
    subset_expr_data = expr_values[rownames(expr_values) %in% selected_feats, ]
    t_sub_expr_data_DT = data.table('selected_feat' = subset_expr_data, 'cell_ID' = colnames(expr_values))
    data.table::setnames(t_sub_expr_data_DT, 'selected_feat', selected_feats)
  } else {
    subset_expr_data = expr_values[rownames(expr_values) %in% selected_feats, ]
    t_sub_expr_data = t_flex(subset_expr_data)
    t_sub_expr_data_DT = data.table::as.data.table(as.matrix(t_sub_expr_data))

    # data.table variables
    cell_ID = NULL

    t_sub_expr_data_DT[, cell_ID := rownames(t_sub_expr_data)]
  }


  ## dimension reduction ##
  dim_dfr = get_dimReduction(gobject = gobject,
                             feat_type = feat_type,
                             spat_unit = spat_unit,
                             reduction = 'cells',
                             reduction_method = dim_reduction_to_use,
                             name = dim_reduction_name,
                             return_dimObj = FALSE)

  #dim_dfr = gobject@dimension_reduction$cells[[dim_reduction_to_use]][[dim_reduction_name]]$coordinates[,c(dim1_to_use, dim2_to_use)]
  dim_names = colnames(dim_dfr)
  dim_DT = data.table::as.data.table(dim_dfr); dim_DT[, cell_ID := rownames(dim_dfr)]

  ## annotated cell metadata
  cell_metadata = pDataDT(gobject,
                          spat_unit = spat_unit,
                          feat_type = feat_type)

  annotated_DT = data.table::merge.data.table(cell_metadata, dim_DT, by = 'cell_ID')

  ## merge feat info
  annotated_feat_DT = data.table::merge.data.table(annotated_DT, t_sub_expr_data_DT, by = 'cell_ID')

  # create input for network
  if(show_NN_network == TRUE) {

    # nn_network
    selected_nn_network = get_NearestNetwork(gobject = gobject,
                                             spat_unit = spat_unit,
                                             feat_type = feat_type,
                                             nn_network_to_use = nn_network_to_use,
                                             network_name = network_name,
                                             output = 'igraph')
    # gobject@nn_network[[nn_network_to_use]][[network_name]][['igraph']]

    network_DT = data.table::as.data.table(igraph::as_data_frame(selected_nn_network, what = 'edges'))

    # annotated network
    old_dim_names = dim_names

    annotated_network_DT = data.table::merge.data.table(network_DT, dim_DT, by.x = 'from', by.y = 'cell_ID')
    from_dim_names = paste0('from_', old_dim_names)
    data.table::setnames(annotated_network_DT, old = old_dim_names, new = from_dim_names)

    annotated_network_DT = data.table::merge.data.table(annotated_network_DT, dim_DT, by.x = 'to', by.y = 'cell_ID')
    to_dim_names = paste0('to_', old_dim_names)
    data.table::setnames(annotated_network_DT, old = old_dim_names, new = to_dim_names)

  }

  ## visualize multiple plots ##
  ## 2D plots ##
  savelist <- list()


  for(feat in selected_feats) {


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


    ## point layer ##
    if(is.null(feats)) {
      cell_color = 'lightblue'
      cat('no feats selected')
      pl <- pl + ggplot2::geom_point(data = annotated_feat_DT,
                                     aes_string(x = dim_names[1], dim_names[2]),
                                     fill = cell_color, show.legend = show_legend,
                                     size =  point_size, alpha = point_alpha)

    } else {


      ## set gradient limits if needed ##
      if(!is.null(gradient_limits) & is.vector(gradient_limits) & length(gradient_limits) == 2) {
        lower_lim = gradient_limits[[1]]
        upper_lim = gradient_limits[[2]]
        numeric_data = annotated_feat_DT[[feat]]
        limit_numeric_data = ifelse(numeric_data > upper_lim, upper_lim,
                                    ifelse(numeric_data < lower_lim, lower_lim, numeric_data))
        annotated_feat_DT[[feat]] = limit_numeric_data
      }

      if(is.null(gradient_midpoint)) {
        gradient_midpoint = stats::median(annotated_feat_DT[[feat]])
      }



      ## with border ##
      if(point_shape == 'border') {

        if(scale_alpha_with_expression == TRUE) {
          pl <- pl + ggplot2::geom_point(data = annotated_feat_DT, aes_string2(x = dim_names[1],
                                                                               y = dim_names[2],
                                                                               fill = feat, alpha = feat),
                                         show.legend = show_legend, shape = 21, size = point_size,
                                         color = point_border_col, stroke = point_border_stroke)
        } else {
          pl <- pl + ggplot2::geom_point(data = annotated_feat_DT, aes_string2(x = dim_names[1],
                                                                               y = dim_names[2],
                                                                               fill = feat),
                                         show.legend = show_legend, shape = 21,
                                         size =  point_size,
                                         color = point_border_col, stroke = point_border_stroke,
                                         alpha = point_alpha)
        }

        ## scale and labs ##
        pl <- pl + ggplot2::scale_alpha_continuous(guide = 'none')
        pl <- pl + ggplot2::scale_fill_gradient2(low = cell_color_gradient[[1]],
                                                 mid = cell_color_gradient[[2]],
                                                 high = cell_color_gradient[[3]],
                                                 midpoint = gradient_midpoint,
                                                 guide = guide_colorbar(title = ''))
      }


      ## without border ##
      if(point_shape == 'no_border') {

        if(scale_alpha_with_expression == TRUE) {
          pl <- pl + ggplot2::geom_point(data = annotated_feat_DT, aes_string2(x = dim_names[1],
                                                                               y = dim_names[2],
                                                                               color = feat, alpha = feat),
                                         show.legend = show_legend, shape = 19, size = point_size)
        } else {
          pl <- pl + ggplot2::geom_point(data = annotated_feat_DT, aes_string2(x = dim_names[1],
                                                                               y = dim_names[2],
                                                                               color = feat),
                                         show.legend = show_legend, shape = 19, size =  point_size,
                                         alpha = point_alpha)
        }

        ## scale and labs ##
        pl <- pl + ggplot2::scale_alpha_continuous(guide = 'none')
        pl <- pl + ggplot2::scale_color_gradient2(low = cell_color_gradient[[1]],
                                                  mid = cell_color_gradient[[2]],
                                                  high = cell_color_gradient[[3]],
                                                  midpoint = gradient_midpoint,
                                                  guide = guide_colorbar(title = ''))

      }
    }

    ## add title
    pl <- pl + ggplot2::labs(x = 'coord x', y = 'coord y', title = feat)

    ## aesthetics
    pl <- pl + ggplot2::theme(plot.title = element_text(hjust = 0.5),
                              legend.title = element_blank(),
                              legend.text = element_text(size = legend_text),
                              axis.title = element_text(size = axis_title),
                              axis.text = element_text(size = axis_text),
                              panel.grid = element_blank(),
                              panel.background = element_rect(fill = background_color))

    savelist[[feat]] <- pl
  }




  # combine plots with cowplot
  combo_plot <- cowplot::plot_grid(plotlist = savelist,
                                   ncol = cow_n_col,
                                   rel_heights = cow_rel_h, rel_widths = cow_rel_w,
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



#' @title dimGenePlot2D
#' @name dimGenePlot2D
#' @description Visualize gene expression according to dimension reduction coordinates
#' @param gobject giotto object
#' @param genes genes to show
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @inheritDotParams dimFeatPlot2D -gobject -feats -default_save_name
#' @return ggplot
#' @details Description of parameters.
#' @family dimension reduction gene expression visualizations
#' @export
dimGenePlot2D <- function(gobject,
                          genes = NULL,
                          default_save_name = 'dimGenePlot2D',
                          ...) {

  .Deprecated(new = "dimFeatPlot2D")

  dimFeatPlot2D(gobject = gobject,
                feat_type = 'rna',
                feats = genes,
                default_save_name = default_save_name,
                ...)
}



#' @title dimGenePlot
#' @name dimGenePlot
#' @description Visualize gene expression according to dimension reduction coordinates
#' @inheritDotParams dimGenePlot2D
#' @return ggplot
#' @details Description of parameters.
#' @family dimension reduction gene expression visualizations
#' @export
#' @seealso \code{\link{dimGenePlot3D}}
dimGenePlot = function(...) {

  .Deprecated(new = "dimGenePlot2D")

  dimGenePlot2D(...)

}


## ** spatial and dim reduction feature plotting ####


#' @title spatDimFeatPlot2D
#' @name spatDimFeatPlot2D
#' @description Visualize cells according to spatial AND dimension reduction coordinates in ggplot mode
#' @param gobject giotto object
#' @param spat_unit spatial unit
#' @param feat_type feature type
#' @param show_image show a tissue background image
#' @param gimage a giotto image
#' @param image_name name of a giotto image
#' @param largeImage_name name of a giottoLargeImage
#' @param expression_values feat expression values to use
#' @param plot_alignment direction to align plot
#' @param feats features to show
#' @param dim_reduction_to_use dimension reduction to use
#' @param dim_reduction_name dimension reduction name
#' @param dim1_to_use dimension to use on x-axis
#' @param dim2_to_use dimension to use on y-axis
#' @param dim_point_shape dim reduction points with border or not (border or no_border)
#' @param dim_point_size dim reduction plot: point size
#' @param dim_point_alpha transparancy of dim. reduction points
#' @param dim_point_border_col color of border around points
#' @param dim_point_border_stroke stroke size of border around points
#' @param show_NN_network show underlying NN network
#' @param show_spatial_network show underlying spatial netwok
#' @param nn_network_to_use type of NN network to use (kNN vs sNN)
#' @param network_name name of NN network to use, if show_NN_network = TRUE
#' @param dim_network_color color of NN network
#' @param dim_edge_alpha dim reduction plot: column to use for alpha of the edges
#' @param scale_alpha_with_expression scale expression with ggplot alpha parameter
#' @param sdimx spatial x-axis dimension name (default = 'sdimx')
#' @param sdimy spatial y-axis dimension name (default = 'sdimy')
#' @param spatial_network_name name of spatial network to use
#' @param spatial_network_color color of spatial network
#' @param show_spatial_grid show spatial grid
#' @param grid_color color of spatial grid
#' @param spatial_grid_name name of spatial grid to use
#' @param spat_point_shape spatial points with border or not (border or no_border)
#' @param spat_point_size spatial plot: point size
#' @param spat_point_alpha transparancy of spatial points
#' @param spat_point_border_col color of border around points
#' @param spat_point_border_stroke stroke size of border around points
#' @param spat_edge_alpha edge alpha
#' @param cell_color_gradient vector with 3 colors for numeric data
#' @param gradient_midpoint midpoint for color gradient
#' @param gradient_limits vector with lower and upper limits
#' @param show_legend show legend
#' @param legend_text size of legend text
#' @param dim_background_color color of plot background for dimension plot
#' @param spat_background_color color of plot background for spatial plot
#' @param vor_border_color border colorr for voronoi plot
#' @param vor_max_radius maximum radius for voronoi 'cells'
#' @param vor_alpha transparancy of voronoi 'cells'
#' @param axis_text size of axis text
#' @param axis_title size of axis title
#' @param cow_n_col cowplot param: how many columns
#' @param cow_rel_h cowplot param: relative height
#' @param cow_rel_w cowplot param: relative width
#' @param cow_align cowplot param: how to align
#' @param show_plot show plots
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters, see \code{\link{showSaveParameters}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return ggplot
#' @details Description of parameters.
#' @family spatial and dimension reduction feature expression visualizations
#' @export
spatDimFeatPlot2D <- function(gobject,
                              spat_unit = NULL,
                              feat_type = NULL,
                              show_image = F,
                              gimage = NULL,
                              image_name = NULL,
                              largeImage_name = NULL,
                              expression_values = c('normalized', 'scaled', 'custom'),
                              plot_alignment = c('vertical', 'horizontal'),
                              feats,
                              dim_reduction_to_use = 'umap',
                              dim_reduction_name = 'umap',
                              dim1_to_use = 1,
                              dim2_to_use = 2,
                              dim_point_shape = c('border', 'no_border'),
                              dim_point_size = 1,
                              dim_point_alpha = 1,
                              dim_point_border_col = 'black',
                              dim_point_border_stroke = 0.1,
                              show_NN_network = F,
                              show_spatial_network = F,
                              dim_network_color = 'gray',
                              nn_network_to_use = 'sNN',
                              network_name = 'sNN.pca',
                              dim_edge_alpha = NULL,
                              scale_alpha_with_expression = FALSE,
                              sdimx = 'sdimx',
                              sdimy = 'sdimy',
                              spatial_network_name = 'Delaunay_network',
                              spatial_network_color = NULL,
                              show_spatial_grid = F,
                              grid_color = NULL,
                              spatial_grid_name = 'spatial_grid',
                              spat_point_shape = c('border', 'no_border', 'voronoi'),
                              spat_point_size = 1,
                              spat_point_alpha = 1,
                              spat_point_border_col = 'black',
                              spat_point_border_stroke = 0.1,
                              spat_edge_alpha = NULL,
                              cell_color_gradient = c('blue', 'white', 'red'),
                              gradient_midpoint = NULL,
                              gradient_limits = NULL,
                              cow_n_col = 2,
                              cow_rel_h = 1,
                              cow_rel_w = 1,
                              cow_align = 'h',
                              show_legend = TRUE,
                              legend_text = 10,
                              dim_background_color = 'white',
                              spat_background_color = 'white',
                              vor_border_color = 'white',
                              vor_max_radius = 200,
                              vor_alpha = 1,
                              axis_text = 8,
                              axis_title = 8,
                              show_plot = NA,
                              return_plot = NA,
                              save_plot = NA,
                              save_param =  list(),
                              default_save_name = 'spatDimFeatPlot2D') {

  plot_alignment = match.arg(plot_alignment, choices = c('vertical', 'horizontal'))

  # dimension reduction plot
  dmpl = dimFeatPlot2D(gobject = gobject,
                       feat_type = feat_type,
                       spat_unit = spat_unit,
                       expression_values = expression_values,
                       feats = feats,
                       dim_reduction_to_use = dim_reduction_to_use,
                       dim_reduction_name = dim_reduction_name,
                       dim1_to_use = dim1_to_use,
                       dim2_to_use = dim2_to_use,
                       show_NN_network = show_NN_network,
                       nn_network_to_use = nn_network_to_use,
                       network_name = network_name,
                       network_color = dim_network_color,
                       edge_alpha = dim_edge_alpha,
                       scale_alpha_with_expression = scale_alpha_with_expression,
                       point_shape = dim_point_shape,
                       point_size = dim_point_size,
                       point_alpha = dim_point_alpha,
                       cell_color_gradient = cell_color_gradient,
                       gradient_midpoint = gradient_midpoint,
                       gradient_limits = gradient_limits,
                       point_border_col = dim_point_border_col,
                       point_border_stroke = dim_point_border_stroke,
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
  spl = spatFeatPlot2D(gobject = gobject,
                       feat_type = feat_type,
                       spat_unit = spat_unit,
                       show_image = show_image,
                       gimage = gimage,
                       image_name = image_name,
                       largeImage_name = largeImage_name,
                       sdimx = sdimx,
                       sdimy = sdimy,
                       expression_values = expression_values,
                       feats = feats,
                       cell_color_gradient = cell_color_gradient,
                       gradient_midpoint = gradient_midpoint,
                       gradient_limits = gradient_limits,
                       show_network = show_spatial_network,
                       network_color = spatial_network_color,
                       spatial_network_name = spatial_network_name,
                       edge_alpha = spat_edge_alpha,
                       show_grid = show_spatial_grid,
                       grid_color = grid_color,
                       spatial_grid_name = spatial_grid_name,
                       scale_alpha_with_expression = scale_alpha_with_expression,
                       point_shape = spat_point_shape,
                       point_size = spat_point_size,
                       point_alpha = spat_point_alpha,
                       point_border_col = spat_point_border_col,
                       point_border_stroke = spat_point_border_stroke,
                       show_legend = show_legend,
                       legend_text = legend_text,
                       background_color = spat_background_color,
                       vor_border_color = vor_border_color,
                       vor_max_radius = vor_max_radius,
                       vor_alpha = vor_alpha,
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






#' @title spatDimGenePlot2D
#' @name spatDimGenePlot2D
#' @description Visualize cells according to spatial AND dimension reduction coordinates in ggplot mode
#' @param gobject giotto object
#' @param genes genes to show
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @inheritDotParams spatDimFeatPlot2D -gobject -feats -default_save_name
#' @return ggplot
#' @details Description of parameters.
#' @family spatial and dimension reduction gene expression visualizations
#' @export
#' @seealso \code{\link{spatDimGenePlot3D}}
spatDimGenePlot2D <- function(gobject,
                              genes,
                              default_save_name = 'spatDimGenePlot2D',
                              ...) {

  .Deprecated(new = "spatDimFeatPlot2D")

  spatDimFeatPlot2D(gobject = gobject,
                    feat_type = 'rna',
                    feats = genes,
                    default_save_name = default_save_name,
                    ...)

}






#' @title spatDimGenePlot
#' @name spatDimGenePlot
#' @description Visualize cells according to spatial AND dimension reduction coordinates in ggplot mode
#' @inheritDotParams spatDimGenePlot2D
#' @return ggplot
#' @details Description of parameters.
#' @family spatial and dimension reduction gene expression visualizations
#' @export
#' @seealso \code{\link{spatDimGenePlot3D}}
spatDimGenePlot = function(...) {

  .Deprecated(new = "spatDimGenePlot2D")

  spatDimGenePlot2D(...)
}






#' @title spatCellPlot2D
#' @name spatCellPlot2D
#' @description Visualize cells according to spatial coordinates
#' @param gobject giotto object
#' @param feat_type feature type
#' @param spat_unit spatial unit
#' @param show_image show a tissue background image
#' @param gimage a giotto image
#' @param image_name name of a giotto image
#' @param largeImage_name name of a giottoLargeImage
#' @param sdimx x-axis dimension name (default = 'sdimx')
#' @param sdimy y-axis dimension name (default = 'sdimy')
#' @param spat_enr_names names of spatial enrichment results to include
#' @param cell_annotation_values numeric cell annotation columns
#' @param cell_color_gradient vector with 3 colors for numeric data
#' @param gradient_midpoint midpoint for color gradient
#' @param gradient_limits vector with lower and upper limits
#' @param select_cell_groups select subset of cells/clusters based on cell_color parameter
#' @param select_cells select subset of cells based on cell IDs
#' @param point_shape shape of points (border, no_border or voronoi)
#' @param point_size size of point (cell)
#' @param point_alpha transparancy of spatial points
#' @param point_border_col color of border around points
#' @param point_border_stroke stroke size of border around points
#' @param show_cluster_center plot center of selected clusters
#' @param show_center_label plot label of selected clusters
#' @param center_point_size size of center points
#' @param center_point_border_col border color of center points
#' @param center_point_border_stroke border stroke size of center points
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
#' @param vor_border_color border colorr for voronoi plot
#' @param vor_max_radius maximum radius for voronoi 'cells'
#' @param vor_alpha transparancy of voronoi 'cells'
#' @param axis_text size of axis text
#' @param axis_title size of axis title
#' @param cow_n_col cowplot param: how many columns
#' @param cow_rel_h cowplot param: relative height
#' @param cow_rel_w cowplot param: relative width
#' @param cow_align cowplot param: how to align
#' @param show_plot show plot
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters, see \code{\link{showSaveParameters}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return ggplot
#' @details Description of parameters.
#' @family spatial cell annotation visualizations
#' @export
spatCellPlot2D = function(gobject,
                          spat_unit = NULL,
                          feat_type = NULL,
                          show_image = F,
                          gimage = NULL,
                          image_name = NULL,
                          largeImage_name = NULL,
                          sdimx = 'sdimx',
                          sdimy = 'sdimy',
                          spat_enr_names = NULL,
                          cell_annotation_values = NULL,
                          cell_color_gradient = c('blue', 'white', 'red'),
                          gradient_midpoint = NULL,
                          gradient_limits = NULL,
                          select_cell_groups = NULL,
                          select_cells = NULL,
                          point_shape = c('border', 'no_border', 'voronoi'),
                          point_size = 3,
                          point_alpha = 1,
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
                          vor_border_color = 'white',
                          vor_max_radius = 200,
                          vor_alpha = 1,
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
                          default_save_name = 'spatCellPlot2D') {


  # Set feat_type and spat_unit
  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = spat_unit)
  feat_type = set_default_feat_type(gobject = gobject,
                                    spat_unit = spat_unit,
                                    feat_type = feat_type)

  comb_metadata = combineMetadata(gobject = gobject,
                                  spat_unit = spat_unit,
                                  feat_type = feat_type,
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
                    spat_unit = spat_unit,
                    feat_type = feat_type,
                    show_image = show_image,
                    gimage = gimage,
                    image_name = image_name,
                    largeImage_name = largeImage_name,
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
                    point_alpha = point_alpha,
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
                    vor_border_color = vor_border_color,
                    vor_max_radius = vor_max_radius,
                    vor_alpha = vor_alpha,
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
#' @inheritDotParams spatCellPlot2D
#' @return ggplot
#' @details Description of parameters.
#' @family spatial cell annotation visualizations
#' @export
spatCellPlot = function(...) {

  spatCellPlot2D(...)

}





#' @title dimCellPlot2D
#' @name dimCellPlot2D
#' @description Visualize cells according to dimension reduction coordinates
#' @param gobject giotto object
#' @param spat_unit spatial unit
#' @param feat_type feature type
#' @param dim_reduction_to_use dimension reduction to use
#' @param dim_reduction_name dimension reduction name
#' @param dim1_to_use dimension to use on x-axis
#' @param dim2_to_use dimension to use on y-axis
#' @param spat_enr_names names of spatial enrichment results to include
#' @param cell_annotation_values numeric cell annotation columns
#' @param show_NN_network show underlying NN network
#' @param nn_network_to_use type of NN network to use (kNN vs sNN)
#' @param network_name name of NN network to use, if show_NN_network = TRUE
#' @param cell_color_code named vector with colors for cell annotation values
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
#' @param center_point_border_col border color of center points
#' @param center_point_border_stroke border stroke size of center points
#' @param label_size  size of labels
#' @param label_fontface font of labels
#' @param edge_alpha column to use for alpha of the edges
#' @param point_shape point with border or not (border or no_border)
#' @param point_size size of point (cell)
#' @param point_alpha transparancy of dim. reduction points
#' @param point_border_col color of border around points
#' @param point_border_stroke stroke size of border around points
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
#' @param save_param list of saving parameters, see \code{\link{showSaveParameters}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return ggplot
#' @details Description of parameters. For 3D plots see \code{\link{dimPlot3D}}
#' @family dimension reduction cell annotation visualizations
#' @export
dimCellPlot2D = function(gobject,
                         spat_unit = NULL,
                         feat_type = NULL,
                         dim_reduction_to_use = 'umap',
                         dim_reduction_name = 'umap',
                         dim1_to_use = 1,
                         dim2_to_use = 2,
                         spat_enr_names = NULL,
                         cell_annotation_values = NULL,
                         show_NN_network = F,
                         nn_network_to_use = 'sNN',
                         network_name = 'sNN.pca',
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
                         point_alpha = 1,
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


  # Set feat_type and spat_unit
  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = spat_unit)
  feat_type = set_default_feat_type(gobject = gobject,
                                    spat_unit = spat_unit,
                                    feat_type = feat_type)

  comb_metadata = combineMetadata(gobject = gobject,
                                  spat_unit = spat_unit,
                                  feat_type = feat_type,
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
                   spat_unit = spat_unit,
                   feat_type = feat_type,
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
                   point_alpha = point_alpha,
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
#' @inheritDotParams dimCellPlot2D -gobject
#' @return ggplot
#' @details Description of parameters. For 3D plots see \code{\link{dimCellPlot2D}}
#' @family dimension reduction cell annotation visualizations
#' @export
dimCellPlot = function(gobject, ...) {

  dimCellPlot2D(gobject = gobject, ...)

}




#' @title spatDimCellPlot2D
#' @name spatDimCellPlot2D
#' @description Visualize numerical features of cells according to spatial AND dimension reduction coordinates in 2D
#' @param gobject giotto object
#' @param spat_unit spatial unit
#' @param feat_type feature type
#' @param show_image show a tissue background image
#' @param gimage a giotto image
#' @param image_name name of a giotto image
#' @param largeImage_name name of a giottoLargeImage
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
#' @param dim_point_alpha transparancy of dim. reduction points
#' @param dim_point_border_col border color of points in dim. reduction space
#' @param dim_point_border_stroke border stroke of points in dim. reduction space
#' @param spat_point_shape shape of points (border, no_border or voronoi)
#' @param spat_point_size size of spatial points
#' @param spat_point_alpha transparancy of spatial points
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
#' @param spat_center_point_size size of the spatial center points
#' @param spat_center_point_border_col border color of the spatial center points
#' @param spat_center_point_border_stroke stroke size of the spatial center points
#' @param spat_label_size size of the center label
#' @param spat_label_fontface font of the center label
#' @param show_NN_network show underlying NN network
#' @param nn_network_to_use type of NN network to use (kNN vs sNN)
#' @param nn_network_name name of NN network to use, if show_NN_network = TRUE
#' @param dim_edge_alpha column to use for alpha of the edges
#' @param spat_show_network show spatial network
#' @param spatial_network_name name of spatial network to use
#' @param spat_network_color color of spatial network
#' @param spat_network_alpha alpha of spatial network
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
#' @param vor_border_color border colorr for voronoi plot
#' @param vor_max_radius maximum radius for voronoi 'cells'
#' @param vor_alpha transparancy of voronoi 'cells'
#' @param axis_text size of axis text
#' @param axis_title size of axis title
#' @param show_plot show plot
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters, see \code{\link{showSaveParameters}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return ggplot
#' @details Description of parameters.
#' @family spatial and dimension reduction cell annotation visualizations
#' @export
spatDimCellPlot2D <- function(gobject,
                              feat_type = NULL,
                              spat_unit = NULL,
                              show_image = F,
                              gimage = NULL,
                              image_name = NULL,
                              largeImage_name = NULL,
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
                              dim_point_alpha = 1,
                              dim_point_border_col = 'black',
                              dim_point_border_stroke = 0.1,
                              spat_point_shape = c('border', 'no_border', 'voronoi'),
                              spat_point_size = 1,
                              spat_point_alpha = 1,
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
                              vor_border_color = 'white',
                              vor_max_radius = 200,
                              vor_alpha = 1,
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
                       spat_unit = spat_unit,
                       feat_type = feat_type,
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
                       point_alpha = dim_point_alpha,
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
                       spat_unit = spat_unit,
                       feat_type = feat_type,
                       show_image = show_image,
                       gimage = gimage,
                       image_name = image_name,
                       largeImage_name = largeImage_name,
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
                       point_alpha = spat_point_alpha,
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
                       vor_border_color = vor_border_color,
                       vor_max_radius = vor_max_radius,
                       vor_alpha = vor_alpha,
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
#' @inheritDotParams spatDimCellPlot2D
#' @return ggplot
#' @details Description of parameters.
#' @family spatial and dimension reduction cell annotation visualizations
#' @export
spatDimCellPlot = function(...) {

  spatDimCellPlot2D(...)

}






# * ####
## 3-D plotly ####
## ----------- ##

# ** dimension plot ####


#' @title dimPlot_2D_plotly
#' @name dimPlot_2D_plotly
#' @description Visualize cells at their 2D dimension reduction coordinates with plotly
#' @return plotly object
#' @keywords internal
dimPlot_2D_plotly <- function(gobject,
                              spat_unit = NULL,
                              feat_type = NULL,
                              dim_reduction_to_use = 'umap',
                              dim_reduction_name = 'umap',
                              dim1_to_use = 1,
                              dim2_to_use = 2,
                              spat_enr_names = NULL,
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


  # Set feat_type and spat_unit
  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = spat_unit)
  feat_type = set_default_feat_type(gobject = gobject,
                                    spat_unit = spat_unit,
                                    feat_type = feat_type)

  # data.table variables
  cell_ID = NULL

  ## dimension reduction ##
  dim_dfr = get_dimReduction(gobject = gobject,
                             spat_unit = spat_unit,
                             feat_type = feat_type,
                             reduction = 'cells',
                             reduction_method = dim_reduction_to_use,
                             name = dim_reduction_name,
                             return_dimObj = FALSE)

  dim_dfr = dim_dfr[,c(dim1_to_use, dim2_to_use)]
  dim_names = colnames(dim_dfr)
  dim_DT = data.table::as.data.table(dim_dfr)
  dim_DT[, cell_ID := rownames(dim_dfr)]


  ## annotated cell metadata
  cell_metadata = combineMetadata(gobject = gobject,
                                  spat_unit = spat_unit,
                                  feat_type = feat_type,
                                  spat_enr_names = spat_enr_names)
  annotated_DT = merge(cell_metadata, dim_DT, by = 'cell_ID')


  # create input for network
  if(show_NN_network == TRUE) {

    # nn_network
    selected_nn_network = get_NearestNetwork(gobject = gobject,
                                             spat_unit = spat_unit,
                                             feat_type = feat_type,
                                             nn_network_to_use = nn_network_to_use,
                                             network_name = network_name,
                                             output = 'igraph')
    network_DT = data.table::as.data.table(igraph::as_data_frame(selected_nn_network, what = 'edges'))

    # annotated network
    old_dim_names = dim_names

    annotated_network_DT = merge(network_DT, dim_DT, by.x = 'from', by.y = 'cell_ID')
    from_dim_names = paste0('from_', old_dim_names)
    data.table::setnames(annotated_network_DT, old = old_dim_names, new = from_dim_names)

    annotated_network_DT <- merge(annotated_network_DT, dim_DT, by.x = 'to', by.y = 'cell_ID')
    to_dim_names = paste0('to_', old_dim_names)
    data.table::setnames(annotated_network_DT, old = old_dim_names, new = to_dim_names)

  }


  if(dim_reduction_to_use == "pca"){

    pca_object = get_dimReduction(gobject = gobject,
                                  spat_unit = spat_unit,
                                  feat_type = feat_type,
                                  reduction = 'cells',
                                  reduction_method = dim_reduction_to_use,
                                  name = dim_reduction_name,
                                  return_dimObj = TRUE)
    eigenvalues = slot(pca_object, 'misc')$eigenvalues

    if(!is.null(eigenvalues)) {
      total = sum(eigenvalues)
      var_expl_vec = (eigenvalues/total) * 100
      dim1_x_variance = var_expl_vec[dim1_to_use]
      dim2_y_variance = var_expl_vec[dim2_to_use]
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
      annotated_network_DT = annotated_network_DT[annotated_network_DT$to %in% select_cells & annotated_network_DT$from %in% select_cells]
    }

    # if specific cells are selected
    #annotated_DT = annotated_DT_selected

  }


  ## if no subsets are required
  if(is.null(select_cells) & is.null(select_cell_groups)) {
    annotated_DT_selected = annotated_DT
    annotated_DT_other    = NULL
  }


  ## annotated_DT_selected = all selected cells or all cells if no selection
  ## annotated_DT_other = all not selected cells or NULL if no selection


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
                                      opacity = edge_alpha)
  }

  if(is.null(cell_color)){
    cell_color = "lightblue"
    pl <- pl %>% plotly::add_trace(type = "scatter",mode = "markers",
                                   x = annotated_DT_selected[[dim_names[1]]],
                                   y = annotated_DT_selected[[dim_names[2]]],
                                   color = cell_color,
                                   colors = cell_color,
                                   marker = list(size = point_size))
  }

  else if(cell_color %in% colnames(annotated_DT_selected)){
    if(is.null(cell_color_code)){
      number_colors = length(unique(annotated_DT[[cell_color]]))
      cell_color_code = getDistinctColors(n = number_colors)
    }
    if(color_as_factor){
      annotated_DT_selected[[cell_color]] <- as.factor(annotated_DT_selected[[cell_color]])
    }


    pl <- pl %>% plotly::add_trace(type = "scatter",mode = "markers",
                                   x = annotated_DT_selected[[dim_names[1]]],
                                   y = annotated_DT_selected[[dim_names[2]]],
                                   color = annotated_DT_selected[[cell_color]],
                                   colors = cell_color_code,
                                   legendgroup = annotated_DT_selected[[cell_color]],
                                   marker = list(size = point_size))

    if(!is.null(select_cells) & show_other_cells){
      pl <- pl %>% plotly::add_trace(type = "scatter",mode = "markers",
                                     x = annotated_DT_other[[dim_names[1]]],
                                     y = annotated_DT_other[[dim_names[2]]],
                                     #legendgroup = annotated_DT[[cell_color]],
                                     marker = list(size = other_point_size,color = other_cell_color),
                                     showlegend = F)
    }

    if(show_cluster_center == TRUE | show_center_label == TRUE) {
      annotated_DT_centers = annotated_DT_selected[, .(center_1 = stats::median(get(dim_names[1])),
                                                       center_2 = stats::median(get(dim_names[2]))),
                                                   by = cell_color]
      annotated_DT_centers[[cell_color]] = as.factor(annotated_DT_centers[[cell_color]])
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

    if(!is.null(eigenvalues)) {
      x_name = paste0('pca','-',dim_names[1])
      y_name = paste0('pca','-',dim_names[2])
      x_title = sprintf('%s explains %.02f%% of variance', x_name, var_expl_vec[1])
      y_title = sprintf('%s explains %.02f%% of variance', y_name, var_expl_vec[2])
    }


  } else {
    x_title = paste(dim_reduction_to_use, dim_names[1],sep = " ")
    y_title = paste(dim_reduction_to_use, dim_names[2],sep = " ")
  }
  pl <- pl %>% plotly::layout(xaxis = list(title = x_title),
                              yaxis = list(title = y_title),
                              legend = list(x = 100, y = 0.5,font = list(family = "sans-serif",size = 12)))

  return (pl)
}


#' @title dimPlot_3D_plotly
#' @name dimPlot_3D_plotly
#' @description Visualize cells at their 3D dimension reduction coordinates with plotly
#' @return plotly object
#' @keywords internal
dimPlot_3D_plotly <- function(gobject,
                              spat_unit = NULL,
                              feat_type = NULL,
                              dim_reduction_to_use = 'umap',
                              dim_reduction_name = 'umap',
                              dim1_to_use = 1,
                              dim2_to_use = 2,
                              dim3_to_use = 3,
                              spat_enr_names = NULL,

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



  # Set feat_type and spat_unit
  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = spat_unit)
  feat_type = set_default_feat_type(gobject = gobject,
                                    spat_unit = spat_unit,
                                    feat_type = feat_type)

  # data.table variables
  cell_ID = NULL

  ## dimension reduction ##
  dim_dfr = get_dimReduction(gobject = gobject,
                             spat_unit = spat_unit,
                             feat_type = feat_type,
                             reduction = 'cells',
                             reduction_method = dim_reduction_to_use,
                             name = dim_reduction_name,
                             return_dimObj = FALSE)
  dim_dfr = dim_dfr[,c(dim1_to_use, dim2_to_use, dim3_to_use)]
  dim_names = colnames(dim_dfr)
  dim_DT = data.table::as.data.table(dim_dfr)
  dim_DT[, cell_ID := rownames(dim_dfr)]


  ## annotated cell metadata
  cell_metadata = combineMetadata(gobject = gobject,
                                  spat_unit = spat_unit,
                                  feat_type = feat_type,
                                  spat_enr_names = spat_enr_names)
  annotated_DT = merge(cell_metadata, dim_DT, by = 'cell_ID')


  # create input for network
  if(show_NN_network == TRUE) {

    # nn_network
    selected_nn_network = get_NearestNetwork(gobject = gobject,spat_unit = spat_unit,
                                             feat_type = feat_type,
                                             nn_network_to_use = nn_network_to_use,
                                             network_name = network_name,
                                             output = 'igraph')
    network_DT = data.table::as.data.table(igraph::as_data_frame(selected_nn_network, what = 'edges'))

    # annotated network
    old_dim_names = dim_names

    annotated_network_DT = merge(network_DT, dim_DT, by.x = 'from', by.y = 'cell_ID')
    from_dim_names = paste0('from_', old_dim_names)
    data.table::setnames(annotated_network_DT, old = old_dim_names, new = from_dim_names)

    annotated_network_DT = merge(annotated_network_DT, dim_DT, by.x = 'to', by.y = 'cell_ID')
    to_dim_names = paste0('to_', old_dim_names)
    data.table::setnames(annotated_network_DT, old = old_dim_names, new = to_dim_names)

  }

  if(dim_reduction_to_use == "pca"){
    pca_object = get_dimReduction(gobject = gobject,
                                  spat_unit = spat_unit,
                                  feat_type = feat_type,
                                  reduction = 'cells',
                                  reduction_method = dim_reduction_to_use,
                                  name = dim_reduction_name,
                                  return_dimObj = TRUE)

    eigenvalues = slot(pca_object, 'misc')$eigenvalues
    if(!is.null(eigenvalues)) {
      total = sum(eigenvalues)
      var_expl_vec = (eigenvalues/total) * 100
      dim1_x_variance = var_expl_vec[dim1_to_use]
      dim2_y_variance = var_expl_vec[dim2_to_use]
    }

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
      annotated_network_DT = annotated_network_DT[annotated_network_DT$to %in% select_cells & annotated_network_DT$from %in% select_cells]
    }

    # if specific cells are selected
    annotated_DT = annotated_DT_selected
  }

  ## if no subsets are required
  if(is.null(select_cells) & is.null(select_cell_groups)) {
    annotated_DT_selected = annotated_DT
    annotated_DT_other    = NULL
  }

  ## annotated_DT_selected = all selected cells or all cells if no selection
  ## annotated_DT_other = all not selected cells or NULL if no selection


  pl <- plotly::plot_ly()
  if(is.null(cell_color)){
    cell_color = "lightblue"
    pl <- pl %>% plotly::add_trace(type = 'scatter3d',mode = "markers",
                                   x = annotated_DT_selected[[dim_names[1]]],
                                   y = annotated_DT_selected[[dim_names[2]]],
                                   z = annotated_DT_selected[[dim_names[3]]],
                                   color = cell_color,
                                   colors = cell_color,
                                   marker = list(size = 2),
                                   legendgroup = annotated_DT_selected[[cell_color]])
  }
  else{
    if(cell_color %in% colnames(annotated_DT_selected)){
      if(is.null(cell_color_code)) {
        number_colors=length(unique(annotated_DT_selected[[cell_color]]))
        cell_color_code = getDistinctColors(n = number_colors)
      }
      if(color_as_factor){
        annotated_DT_selected[[cell_color]] <- as.factor(annotated_DT_selected[[cell_color]])
      }

      pl <- pl %>% plotly::add_trace(type = 'scatter3d',mode = "markers",
                                     x = annotated_DT_selected[[dim_names[1]]],
                                     y = annotated_DT_selected[[dim_names[2]]],
                                     z = annotated_DT_selected[[dim_names[3]]],
                                     color = annotated_DT_selected[[cell_color]],
                                     colors = cell_color_code,
                                     marker = list(size = point_size),
                                     legendgroup = annotated_DT_selected[[cell_color]])

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
        annotated_DT_centers = annotated_DT_selected[, .(center_1 = stats::median(get(dim_names[1])),
                                                         center_2 = stats::median(get(dim_names[2])),
                                                         center_3 = stats::median(get(dim_names[3]))),
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

    if(!is.null(eigenvalues)) {
      x_name = paste0('pca','-',dim_names[1])
      y_name = paste0('pca','-',dim_names[2])
      z_name = paste0('pca','-',dim_names[3])
      x_title = sprintf('%s explains %.02f%% of variance', x_name, var_expl_vec[1])
      y_title = sprintf('%s explains %.02f%% of variance', y_name, var_expl_vec[2])
      z_title = sprintf('%s explains %.02f%% of variance', z_name, var_expl_vec[3])
    }


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






#' @title dimPlot3D
#' @name dimPlot3D
#' @description Visualize cells according to dimension reduction coordinates
#' @param gobject giotto object
#' @param spat_unit spatial unit
#' @param feat_type feature type
#' @param dim_reduction_to_use dimension reduction to use
#' @param dim_reduction_name dimension reduction name
#' @param dim1_to_use dimension to use on x-axis
#' @param dim2_to_use dimension to use on y-axis
#' @param dim3_to_use dimension to use on z-axis
#' @param spat_enr_names names of spatial enrichment results to include
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
#' @param show_plot show plot
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters, see \code{\link{showSaveParameters}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return plotly
#' @details Description of parameters.
#' @family reduced dimension visualizations
#' @export
dimPlot3D = function(gobject,
                     spat_unit = NULL,
                     feat_type = NULL,
                     dim_reduction_to_use = 'umap',
                     dim_reduction_name = 'umap',
                     dim1_to_use = 1,
                     dim2_to_use = 2,
                     dim3_to_use = 3,
                     spat_enr_names = NULL,

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

    pl = dimPlot_2D_plotly(gobject = gobject,
                           spat_unit = spat_unit,
                           feat_type = feat_type,
                           dim_reduction_to_use = dim_reduction_to_use,
                           dim_reduction_name = dim_reduction_name,
                           dim1_to_use = dim1_to_use,
                           dim2_to_use = dim2_to_use,
                           spat_enr_names = spat_enr_names,

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
    pl = dimPlot_3D_plotly(gobject = gobject,
                           spat_unit = spat_unit,
                           feat_type = feat_type,
                           dim_reduction_to_use = dim_reduction_to_use,
                           dim_reduction_name = dim_reduction_name,
                           dim1_to_use = dim1_to_use,
                           dim2_to_use = dim2_to_use,
                           dim3_to_use = dim3_to_use,
                           spat_enr_names = spat_enr_names,

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
#' @param dim_reduction_name name of UMAP
#' @param default_save_name default save name of UMAP plot
#' @inheritDotParams dimPlot3D -gobject -dim_reduction_to_use -dim_reduction_name -default_save_name
#' @return plotly
#' @details Description of parameters.
#' @family reduced dimension visualizations
#' @export
plotUMAP_3D = function(gobject,
                       dim_reduction_name = 'umap',
                       default_save_name = 'UMAP_3D',
                       ...) {

  dimPlot3D(gobject = gobject,
            dim_reduction_to_use = 'umap',
            dim_reduction_name = dim_reduction_name,
            default_save_name = default_save_name,
            ...)

}


#' @title plotTSNE_3D
#' @name plotTSNE_3D
#' @description Visualize cells according to dimension reduction coordinates
#' @param gobject giotto object
#' @param dim_reduction_name name of TSNE
#' @param default_save_name default save name of TSNE plot
#' @inheritDotParams dimPlot3D -gobject -dim_reduction_to_use -dim_reduction_name -default_save_name
#' @return plotly
#' @details Description of parameters.
#' @family reduced dimension visualizations
#' @export
plotTSNE_3D = function(gobject,
                       dim_reduction_name = 'tsne',
                       default_save_name = 'TSNE_3D',
                       ...) {

  dimPlot3D(gobject = gobject,
            dim_reduction_to_use = 'tsne',
            dim_reduction_name = dim_reduction_name,
            default_save_name = default_save_name,
            ...)

}


#' @title plotPCA_3D
#' @name plotPCA_3D
#' @description Visualize cells according to 3D PCA dimension reduction
#' @param gobject giotto object
#' @param dim_reduction_name name of PCA
#' @param default_save_name default save name of PCA plot
#' @inheritDotParams dimPlot3D -gobject -dim_reduction_to_use -dim_reduction_name -default_save_name
#' @return plotly
#' @details Description of parameters.
#' @family reduced dimension visualizations
#' @export
plotPCA_3D = function(gobject,
                      dim_reduction_name = 'pca',
                      default_save_name = 'PCA_3D',
                      ...) {

  dimPlot3D(gobject = gobject,
            dim_reduction_to_use = 'pca',
            dim_reduction_name = dim_reduction_name,
            default_save_name = default_save_name,
            ...)

}






# ** ####
# ** spatial 3D plot ####

#' @title spatPlot_2D_plotly
#' @name spatPlot_2D_plotly
#' @description Visualize cells at their 2D spatial locations with plotly
#' @return plotly object
#' @keywords internal
spatPlot_2D_plotly = function(gobject,
                              spat_unit = NULL,
                              feat_type = NULL,
                              spat_loc_name = NULL,
                              sdimx = NULL,
                              sdimy = NULL,
                              spat_enr_names = NULL,
                              point_size = 3,
                              cell_color = NULL,
                              cell_color_code = NULL,
                              color_as_factor = T,
                              select_cell_groups = NULL,
                              select_cells = NULL,
                              show_other_cells = T,
                              other_cell_color = "lightgrey",
                              other_point_size = 0.5,
                              show_network = FALSE,
                              spatial_network_name = 'spatial_network',
                              network_color = "lightgray",
                              network_alpha = 1,
                              other_cell_alpha = 0.5,
                              show_grid = FALSE,
                              spatial_grid_name = 'spatial_grid',
                              grid_color = NULL,
                              grid_alpha = 1,
                              show_legend = T,
                              axis_scale = c("cube","real","custom"),
                              custom_ratio = NULL,
                              x_ticks = NULL,
                              y_ticks = NULL,
                              show_plot = F) {

  # Set feat_type and spat_unit
  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = spat_unit)
  feat_type = set_default_feat_type(gobject = gobject,
                                    spat_unit = spat_unit,
                                    feat_type = feat_type)

  ## get spatial cell locations
  cell_locations  = get_spatial_locations(gobject,
                                          spat_unit = spat_unit,
                                          spat_loc_name = spat_loc_name)


  ## extract spatial network
  if(show_network == TRUE) {
    spatial_network = get_spatialNetwork(gobject,
                                         spat_unit = spat_unit,
                                         name = spatial_network_name,
                                         return_network_Obj = FALSE)
  } else {
    spatial_network = NULL
  }

  ## extract spatial grid
  if(show_grid == TRUE) {
    spatial_grid = get_spatialGrid(gobject,
                                   spat_unit = spat_unit,
                                   feat_type = feat_type,
                                   spatial_grid_name)
  } else {
    spatial_grid = NULL
  }

  ## get cell metadata
  cell_metadata = combineMetadata(gobject = gobject,
                                  feat_type = feat_type,
                                  spat_unit = spat_unit,
                                  spat_enr_names = spat_enr_names)

  if(nrow(cell_metadata) == 0) {
    cell_locations_metadata = cell_locations
  } else {
    cell_locations_metadata = cell_metadata
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



  ### set scale
  axis_scale = match.arg(axis_scale, c("cube","real","custom"))

  ### set ratio
  ratio = plotly_axis_scale_2D(cell_locations,
                               sdimx = sdimx,
                               sdimy = sdimy,
                               mode = axis_scale,
                               custom_ratio = custom_ratio)



  pl <- plotly::plot_ly()

  ## create network
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
                                        opacity = network_alpha)
    }
  }

  ## create grid
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
    if(cell_color %in% colnames(cell_locations_metadata_selected)){
      if(is.null(cell_color_code)) {
        number_colors=length(unique(cell_locations_metadata_selected[[cell_color]]))
        cell_color_code = getDistinctColors(n = number_colors)
      }
      cell_locations_metadata_selected[[cell_color]] <- as.factor(cell_locations_metadata_selected[[cell_color]])
      pl <- pl %>% plotly::add_trace(type = 'scatter',
                                     mode = 'markers',
                                     x = cell_locations_metadata_selected[[sdimx]],
                                     y = cell_locations_metadata_selected[[sdimy]],
                                     color = cell_locations_metadata_selected[[cell_color]],
                                     colors = cell_color_code,
                                     marker = list(size = point_size))


      if(!is.null(select_cells) & show_other_cells){
        pl <- pl %>% plotly::add_trace(type = "scatter",
                                       mode="markers",
                                       data= cell_locations_metadata_other,
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
                                   x = cell_locations_metadata_selected[[sdimx]],
                                   y = cell_locations_metadata_selected[[sdimy]],
                                   colors = 'lightblue',
                                   marker = list(size = point_size))

    if(!is.null(select_cells) & show_other_cells){
      pl <- pl %>% plotly::add_trace(type = "scatter",
                                     mode = "markers",
                                     data = cell_locations_metadata_other,
                                     name = "unselected cells",
                                     x =~ sdimx,
                                     y =~ sdimy,
                                     marker = list(size = other_point_size,
                                                   color = other_cell_color),
                                     opacity = other_cell_alpha)
    }
  }


  pl <- pl %>%
    plotly::layout(list(xaxis = list(title = 'X', nticks = x_ticks),
                        yaxis = list(title = 'Y', nticks = y_ticks)),
                   legend = list(x = 100, y = 0.5,
                                 font = list(family = "sans-serif",
                                             size = 12)))


  return((pl))


}



#' @title spatPlot_3D_plotly
#' @name spatPlot_3D_plotly
#' @description Visualize cells at their 3D spatial locations with plotly
#' @return plotly object
#' @keywords internal
spatPlot_3D_plotly = function(gobject,
                              spat_unit = NULL,
                              feat_type = NULL,
                              spat_loc_name = NULL,
                              sdimx = NULL,
                              sdimy = NULL,
                              sdimz = NULL,
                              spat_enr_names = NULL,
                              point_size = 3,
                              cell_color = NULL,
                              cell_color_code = NULL,
                              select_cell_groups = NULL,
                              select_cells = NULL,
                              show_other_cells = T,
                              other_cell_color = "lightgrey",
                              other_point_size = 0.5,
                              show_network = FALSE,
                              spatial_network_name = 'spatial_network',
                              network_color = NULL,
                              network_alpha = 1,
                              other_cell_alpha = 0.5,
                              show_grid = FALSE,
                              spatial_grid_name = 'spatial_grid',
                              title = '',
                              show_legend = TRUE,
                              axis_scale = c("cube","real","custom"),
                              custom_ratio = NULL,
                              x_ticks = NULL,
                              y_ticks = NULL,
                              z_ticks = NULL,
                              show_plot = FALSE) {


  # Set feat_type and spat_unit
  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = spat_unit)
  feat_type = set_default_feat_type(gobject = gobject,
                                    spat_unit = spat_unit,
                                    feat_type = feat_type)

  ## get spatial cell locations
  cell_locations  = get_spatial_locations(gobject,
                                          spat_unit = spat_unit,
                                          spat_loc_name = spat_loc_name)

  ## extract spatial network
  if(show_network == TRUE) {
    spatial_network = get_spatialNetwork(gobject,
                                         spat_unit = spat_unit,
                                         name = spatial_network_name,
                                         return_network_Obj = FALSE)
  } else {
    spatial_network = NULL
  }

  ## extract spatial grid
  if(show_grid == TRUE) {
    spatial_grid = get_spatialGrid(gobject,
                                   spat_unit = spat_unit,
                                   feat_type = feat_type,
                                   spatial_grid_name)
  } else {
    spatial_grid = NULL
  }

  ## get cell metadata
  cell_metadata = combineMetadata(gobject = gobject,
                                  feat_type = feat_type,
                                  spat_unit = spat_unit,
                                  spat_enr_names = spat_enr_names)

  if(nrow(cell_metadata) == 0) {
    cell_locations_metadata = cell_locations
  } else {
    cell_locations_metadata = cell_metadata
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
    #cell_locations_metadata = cell_locations_metadata_selected

  } else if(is.null(select_cells)) {

    cell_locations_metadata_selected = cell_locations_metadata
    cell_locations_metadata_other = NULL

  }



  ### set scale
  axis_scale = match.arg(axis_scale, c("cube","real","custom"))

  ### set ratio
  ratio = plotly_axis_scale_3D(cell_locations,
                               sdimx = sdimx,
                               sdimy = sdimy,
                               sdimz = sdimz,
                               mode = axis_scale,
                               custom_ratio = custom_ratio)



  pl <- plotly::plot_ly()
  if(!is.null(cell_color)) {
    if(cell_color %in% colnames(cell_locations_metadata_selected)){
      if(is.null(cell_color_code)) {
        number_colors=length(unique(cell_locations_metadata_selected[[cell_color]]))
        cell_color_code = getDistinctColors(n = number_colors)
      }
      cell_locations_metadata_selected[[cell_color]] <- as.factor(cell_locations_metadata_selected[[cell_color]])
      pl <- pl %>% plotly::add_trace(type = 'scatter3d',mode = "markers",data = cell_locations_metadata_selected,
                                     x = ~sdimx, y = ~sdimy, z = ~sdimz,
                                     color = cell_locations_metadata_selected[[cell_color]],
                                     colors = cell_color_code,
                                     marker = list(size = point_size))


      if(!is.null(select_cells) & show_other_cells){
        pl <- pl %>% plotly::add_trace(type = "scatter3d",mode="markers",
                                       data = cell_locations_metadata_other,
                                       name = "unselected cells",
                                       x = ~sdimx,
                                       y = ~sdimy,
                                       z = ~sdimz,
                                       marker = list(size = other_point_size,
                                                     color = other_cell_color),
                                       opacity=other_cell_alpha)
      }
    }
    else{
      cat('cell_color does not exist! \n')
    }
  } else {
    pl <- pl %>% plotly::add_trace(type = 'scatter3d',
                                   data = cell_locations_metadata_selected,
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
                                     marker = list(size = other_point_size,
                                                   color=other_cell_color),
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
                                   x = ~x,
                                   y = ~y,
                                   z = ~z,
                                   line = list(color=network_color,width = 0.5),
                                   opacity = network_alpha)
  }

  ## plot spatial grid
  # 3D grid is not clear to view


  pl <- pl %>%
    plotly::layout(scene = list(xaxis = list(title = 'X',nticks = x_ticks),
                                yaxis = list(title = 'Y',nticks = y_ticks),
                                zaxis = list(title = 'Z',nticks = z_ticks),
                                aspectmode='manual',
                                aspectratio = list(x=ratio[[1]],
                                                   y=ratio[[2]],
                                                   z=ratio[[3]])),
                   legend = list(x = 100, y = 0.5,
                                 font = list(family = "sans-serif",size = 12)))


  return(pl)

}




#' @title spatPlot3D
#' @name spatPlot3D
#' @description Visualize cells according to spatial coordinates
#' @param gobject giotto object
#' @param spat_unit spatial unit
#' @param feat_type feature type
#' @param sdimx x-axis dimension name (default = 'sdimx')
#' @param sdimy y-axis dimension name (default = 'sdimy')
#' @param sdimz z-axis dimension name (default = 'sdimy')
#' @param spat_enr_names names of spatial enrichment results to include
#' @param point_size size of point (cell)
#' @param cell_color color for cells (see details)
#' @param cell_color_code named vector with colors
#' @param select_cell_groups select subset of cells/clusters based on cell_color parameter
#' @param select_cells select subset of cells based on cell IDs
#' @param show_other_cells display not selected cells
#' @param other_cell_color color of not selected cells
#' @param other_point_size size of not selected cells
#' @param other_cell_alpha alpha of not selected cells
#' @param show_network show underlying spatial network
#' @param network_color color of spatial network
#' @param network_alpha opacity of spatial network
#' @param spatial_network_name name of spatial network to use
#' @param show_grid show spatial grid
#' @param grid_color color of spatial grid
#' @param grid_alpha opacity of spatial grid
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
#' @param save_param list of saving parameters, see \code{\link{showSaveParameters}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return ggplot
#' @family spatial visualizations
#' @export
spatPlot3D = function(gobject,
                      spat_unit = NULL,
                      feat_type = NULL,
                      sdimx = "sdimx",
                      sdimy = "sdimy",
                      sdimz = "sdimz",
                      spat_enr_names = NULL,
                      point_size = 3,
                      cell_color = NULL,
                      cell_color_code = NULL,
                      select_cell_groups = NULL,
                      select_cells = NULL,
                      show_other_cells = T,
                      other_cell_color = "lightgrey",
                      other_point_size = 0.5,
                      other_cell_alpha = 0.5,
                      show_network = F,
                      spatial_network_name = 'Delaunay_network',
                      network_color = NULL,
                      network_alpha = 1,
                      show_grid = F,
                      spatial_grid_name = 'spatial_grid',
                      grid_color = NULL,
                      grid_alpha = 1,
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

    pl = spatPlot_2D_plotly(gobject = gobject,
                            spat_unit = spat_unit,
                            feat_type = feat_type,
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
    pl = spatPlot_3D_plotly(gobject = gobject,
                            spat_unit = spat_unit,
                            feat_type = feat_type,
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









# ** ####
# ** spatial & dimension 3D plot ####

#' @title spatDimPlot3D
#' @name spatDimPlot3D
#' @description Visualize cells according to spatial AND dimension reduction coordinates in plotly mode
#' @param gobject giotto object
#' @param spat_unit spatial unit
#' @param feat_type feature type
#' @param plot_alignment direction to align plot
#' @param dim_reduction_to_use dimension reduction to use
#' @param dim_reduction_name dimension reduction name
#' @param dim1_to_use dimension to use on x-axis
#' @param dim2_to_use dimension to use on y-axis
#' @param dim3_to_use dimension to use on z-axis
#'
#' @param spat_loc_name name for spatial locations
#' @param sdimx = spatial dimension to use on x-axis
#' @param sdimy = spatial dimension to use on y-axis
#' @param sdimz = spatial dimension to use on z-axis
#'
#' @param spat_enr_names names of spatial enrichment results to include
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
#'
#' @param select_cell_groups select subset of cells/clusters based on cell_color parameter
#' @param select_cells select subset of cells based on cell IDs
#' @param show_other_cells display not selected cells
#' @param other_cell_color color of not selected cells
#' @param other_point_size size of not selected cells
#'
#' @param dim_point_size size of points in dim. reduction space
#' @param nn_network_color color of nn network
#' @param nn_network_alpha column to use for alpha of the edges
#' @param show_spatial_network show spatial network
#' @param spatial_network_name name of spatial network to use
#' @param spatial_network_color color of spatial network
#'
#' @param show_spatial_grid show spatial grid
#' @param spatial_grid_name name of spatial grid to use
#' @param spatial_grid_color color of spatial grid
#' @param spatial_grid_alpha alpha of spatial grid
#' @param spatial_point_size size of spatial points
#' @param spatial_network_color color of spatial network
#' @param spatial_network_alpha alpha of spatial network
#'
#' @param axis_scale the way to scale the axis
#' @param custom_ratio customize the scale of the plot
#' @param x_ticks set the number of ticks on the x-axis
#' @param y_ticks set the number of ticks on the y-axis
#' @param z_ticks set the number of ticks on the z-axis
#' @param legend_text_size size of legend
#'
#' @param show_plot show plot
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters, see \code{\link{showSaveParameters}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return plotly
#' @details Description of parameters.
#' @family spatial and dimension reduction visualizations
#' @export
spatDimPlot3D <- function(gobject,
                          spat_unit = NULL,
                          feat_type = NULL,
                          plot_alignment = c('horizontal','vertical'),
                          dim_reduction_to_use = 'umap',
                          dim_reduction_name = 'umap',
                          dim1_to_use = 1,
                          dim2_to_use = 2,
                          dim3_to_use = 3,

                          spat_loc_name = NULL,
                          sdimx="sdimx",
                          sdimy="sdimy",
                          sdimz="sdimz",
                          spat_enr_names = NULL,

                          show_NN_network = FALSE,
                          nn_network_to_use = 'sNN',
                          network_name = 'sNN.pca',
                          nn_network_color = 'lightgray',
                          nn_network_alpha = 0.5,

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


                          show_spatial_network = F,
                          spatial_network_name = 'Delaunay_network',
                          spatial_network_color = "lightgray",
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

  # Set feat_type and spat_unit
  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = spat_unit)
  feat_type = set_default_feat_type(gobject = gobject,
                                    spat_unit = spat_unit,
                                    feat_type = feat_type)

  # data.table variables
  cell_ID = NULL

  plot_alignment = match.arg(plot_alignment, choices = c( 'horizontal','vertical'))

  # ********data prepare********#
  ## dimension reduction ##
  dim_dfr = get_dimReduction(gobject = gobject,
                             spat_unit = spat_unit,
                             feat_type = feat_type,
                             reduction = 'cells',
                             reduction_method = dim_reduction_to_use,
                             name = dim_reduction_name,
                             return_dimObj = FALSE)
  dim_dfr = dim_dfr[,c(dim1_to_use, dim2_to_use, dim3_to_use)]
  dim_names = colnames(dim_dfr)
  dim_DT = data.table::as.data.table(dim_dfr)
  dim_DT[, cell_ID := rownames(dim_dfr)]


  ## annotated cell metadata
  cell_metadata = combineMetadata(gobject = gobject,
                                  feat_type = feat_type,
                                  spat_unit = spat_unit,
                                  spat_enr_names = spat_enr_names)
  annotated_DT = merge(cell_metadata, dim_DT, by = 'cell_ID')
  spatial_locations = get_spatial_locations(gobject,
                                            spat_unit = spat_unit,
                                            spat_loc_name = spat_loc_name)
  annotated_DT = merge(annotated_DT, spatial_locations, by = 'cell_ID')


  if(dim_reduction_to_use == "pca"){
    pca_object = get_dimReduction(gobject = gobject,
                                  spat_unit = spat_unit,
                                  feat_type = feat_type,
                                  reduction = 'cells',
                                  reduction_method = dim_reduction_to_use,
                                  name = dim_reduction_name,
                                  return_dimObj = TRUE)
    eigenvalues = slot(pca_object, 'misc')$eigenvalues

    if(!is.null(eigenvalues)) {
      total = sum(eigenvalues)
      var_expl_vec = (eigenvalues/total) * 100
      dim1_x_variance = var_expl_vec[dim1_to_use]
      dim2_y_variance = var_expl_vec[dim2_to_use]
      if(!is.null(dim3_to_use)){
        dim3_z_variance = var_expl_vec[3]
      }
    }


  }



  ## nn network
  if(show_NN_network){

    # nn_network
    selected_nn_network = get_NearestNetwork(gobject = gobject,
                                             feat_type = feat_type,
                                             spat_unit = spat_unit,
                                             nn_network_to_use = nn_network_to_use,
                                             network_name = network_name,
                                             output = 'igraph')
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
  if(show_spatial_network == TRUE) {
    spatial_network = get_spatialNetwork(gobject,
                                         spat_unit = spat_unit,
                                         name = spatial_network_name,
                                         return_network_Obj = FALSE)
  } else {
    spatial_network = NULL
  }


  ## extract spatial grid
  if(show_spatial_grid == TRUE) {
    spatial_grid = get_spatialGrid(gobject,
                                   spat_unit = spat_unit,
                                   feat_type = feat_type,
                                   spatial_grid_name)
  } else {
    spatial_grid = NULL
  }


  # create matching cell_color_code
  if(is.null(cell_color_code)) {
    if(is.character(cell_color)) {

      cell_metadata = pDataDT(gobject,
                              feat_type = feat_type,
                              spat_unit = spat_unit)
      if(cell_color %in% colnames(cell_metadata)) {

        if(color_as_factor == TRUE) {
          number_colors = length(unique(cell_metadata[[cell_color]]))
          cell_color_code = getDistinctColors(n = number_colors)
          names(cell_color_code) = unique(cell_metadata[[cell_color]])
        }
      }
    }
  }


  ## subset cell selection ##
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
    # annotated_DT = annotated_DT_selected
  }


  ## if no subsets are required
  if(is.null(select_cells) & is.null(select_cell_groups)) {
    annotated_DT_selected = annotated_DT
    annotated_DT_other    = NULL
  }

  ## annotated_DT_selected = all selected cells or all cells if no selection
  ## annotated_DT_other = all not selected cells or NULL if no selection



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
                                          line = list(color = nn_network_color,
                                                      width = 0.5),
                                          opacity=nn_network_alpha)
    }

    if(is.null(cell_color)){
      #cell_color = "lightblue"
      dpl <- dpl %>% plotly::add_trace(type = "scatter",mode = "markers",
                                       x = annotated_DT_selected[[dim_names[1]]],
                                       y = annotated_DT_selected[[dim_names[2]]],
                                       #color = "lightblue",
                                       #colors ="lightblue",
                                       marker = list(size = dim_point_size,
                                                     color = "lightblue"),
                                       showlegend = F)
    }

    else if(cell_color %in% colnames(annotated_DT_selected)){
      if(color_as_factor){
        annotated_DT_selected[[cell_color]] <- as.factor(annotated_DT_selected[[cell_color]])
      }


      dpl <- dpl %>% plotly::add_trace(type = "scatter",mode = "markers",
                                       x = annotated_DT_selected[[dim_names[1]]],
                                       y = annotated_DT_selected[[dim_names[2]]],
                                       color = annotated_DT_selected[[cell_color]],
                                       colors = cell_color_code,
                                       legendgroup = annotated_DT_selected[[cell_color]],
                                       marker = list(size = dim_point_size))
    }

    else{
      stop("cell_color does not exist!\n")
    }


    if((show_cluster_center == TRUE | show_center_label == TRUE) & !is.null(cell_color)) {
      annotated_DT_centers = annotated_DT_selected[, .(center_1 = stats::median(get(dim_names[1])),
                                                       center_2 = stats::median(get(dim_names[2]))),
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

      if(!is.null(eigenvalues)) {
        x_name = paste0('pca','-',dim_names[1])
        y_name = paste0('pca','-',dim_names[2])
        x_title = sprintf('%s explains %.02f%% of variance', x_name, var_expl_vec[1])
        y_title = sprintf('%s explains %.02f%% of variance', y_name, var_expl_vec[2])
      }

    } else {
      x_title = paste(dim_reduction_to_use, dim_names[1],sep = " ")
      y_title = paste(dim_reduction_to_use, dim_names[2],sep = " ")
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
                                       x = annotated_DT_selected[[dim_names[1]]],
                                       y = annotated_DT_selected[[dim_names[2]]],
                                       z = annotated_DT_selected[[dim_names[3]]],
                                       color = "lightblue",
                                       colors = "lightblue",
                                       marker = list(size = dim_point_size),
                                       showlegend = F)
      #legendgroup = annotated_DT_selected[[cell_color]])
    }
    else{
      if(cell_color %in% colnames(annotated_DT_selected)){
        if(is.null(cell_color_code)) {
          number_colors=length(unique(annotated_DT_selected[[cell_color]]))
          cell_color_code = getDistinctColors(n = number_colors)
        }
        if(color_as_factor){
          annotated_DT_selected[[cell_color]] <- as.factor(annotated_DT_selected[[cell_color]])
        }
        dpl <- dpl %>% plotly::add_trace(type = 'scatter3d',mode = "markers",
                                         x = annotated_DT_selected[[dim_names[1]]],
                                         y = annotated_DT_selected[[dim_names[2]]],
                                         z = annotated_DT_selected[[dim_names[3]]],
                                         color = annotated_DT_selected[[cell_color]],
                                         colors = cell_color_code,
                                         marker = list(size = dim_point_size),
                                         legendgroup = annotated_DT_selected[[cell_color]])
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
                                       line=list(color = nn_network_color),
                                       opacity=nn_network_alpha)
    }
    if((show_cluster_center == TRUE | show_center_label == TRUE)& !is.null(cell_color)){
      annotated_DT_centers = annotated_DT_selected[, .(center_1 = stats::median(get(dim_names[1])),
                                                       center_2 = stats::median(get(dim_names[2])),
                                                       center_3 = stats::median(get(dim_names[3]))),
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

  ## 2D plot ##
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
                                            line = list(color = spatial_network_color,
                                                        width = 0.5),
                                            opacity = spatial_network_alpha)
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
                                       x = annotated_DT_selected[[sdimx]],
                                       y = annotated_DT_selected[[sdimy]],
                                       #color = "lightblue",
                                       #colors = "lightblue",
                                       marker = list(size = spatial_point_size,
                                                     color = "lightblue"),
                                       showlegend = F)
    }

    else if(cell_color %in% colnames(annotated_DT_selected)){
      if(color_as_factor){
        annotated_DT_selected[[cell_color]] <- as.factor(annotated_DT_selected[[cell_color]])
      }


      spl <- spl %>% plotly::add_trace(type = "scatter",mode = "markers",
                                       x = annotated_DT_selected[[sdimx]],
                                       y = annotated_DT_selected[[sdimy]],
                                       color = annotated_DT_selected[[cell_color]],
                                       colors = cell_color_code,
                                       legendgroup = annotated_DT_selected[[cell_color]],
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


  ## 3D plot ##
  else{
    axis_scale = match.arg(axis_scale, c("cube","real","custom"))

    ratio = plotly_axis_scale_3D(annotated_DT_selected, sdimx = sdimx,sdimy = sdimy,sdimz = sdimz,
                                 mode = axis_scale,custom_ratio = custom_ratio)
    spl <- plotly::plot_ly(scene = "scene2")
    if(!is.null(cell_color)) {
      if(cell_color %in% colnames(annotated_DT_selected)){
        annotated_DT_selected[[cell_color]] <- as.factor(annotated_DT_selected[[cell_color]])
        spl <- spl %>% plotly::add_trace(type = 'scatter3d',mode = 'markers',
                                         x = annotated_DT_selected[[sdimx]],
                                         y = annotated_DT_selected[[sdimy]],
                                         z = annotated_DT_selected[[sdimz]],
                                         color = annotated_DT_selected[[cell_color]],
                                         colors = cell_color_code,
                                         legendgroup = annotated_DT_selected[[cell_color]],
                                         marker = list(size = spatial_point_size),
                                         showlegend = F)
      }
      else{
        stop("cell_color doesn't exist!\n")
      }
    }
    else{
      spl <- spl %>% plotly::add_trace(type = 'scatter3d',mode = 'markers',
                                       x = annotated_DT_selected$sdimx,
                                       y = annotated_DT_selected$sdimy,
                                       z = annotated_DT_selected$sdimz,
                                       color = "lightblue",
                                       colors = "lightblue",
                                       #legendgroup = annotated_DT_selected[[cell_color]],
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
                                         line=list(color = spatial_network_color),
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




# ** ####
# ** gene 3D plot ####


#' @title spatGenePlot3D
#' @name spatGenePlot3D
#' @description Visualize cells and gene expression according to spatial coordinates
#' @param gobject giotto object
#' @param feat_type feature type
#' @param spat_unit spatial unit
#' @param spat_loc_name name of spatial locations to use
#' @param expression_values gene expression values to use
#' @param genes genes to show
#'
#' @param cluster_column cluster column to select groups
#' @param select_cell_groups select subset of cells/clusters based on cell_color parameter
#' @param select_cells select subset of cells based on cell IDs
#' @param show_other_cells display not selected cells
#' @param other_cell_color color of not selected cells
#' @param other_point_size size of not selected cells
#'
#' @param genes_high_color color represents high gene expression
#' @param genes_mid_color color represents middle gene expression
#' @param genes_low_color color represents low gene expression
#' @param show_network show underlying spatial network
#' @param network_color color of spatial network
#' @param spatial_network_name name of spatial network to use
#' @param edge_alpha alpha of edges
#' @param show_grid show spatial grid
#' @param spatial_grid_name name of spatial grid to use
#'
#' @param point_size size of point (cell)
#' @param show_legend show legend
#'
#' @param axis_scale the way to scale the axis
#' @param custom_ratio customize the scale of the plot
#' @param x_ticks set the number of ticks on the x-axis
#' @param y_ticks set the number of ticks on the y-axis
#' @param z_ticks set the number of ticks on the z-axis
#'
#' @param show_plot show plots
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters, see \code{\link{showSaveParameters}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return ggplot
#' @details Description of parameters.
#' @family spatial gene expression visualizations
#' @export
spatGenePlot3D <- function(gobject,
                           spat_unit = NULL,
                           feat_type = NULL,
                           spat_loc_name = NULL,
                           expression_values = c('normalized', 'scaled', 'custom'),
                           genes,

                           show_network = FALSE,
                           network_color = NULL,
                           spatial_network_name = 'Delaunay_network',
                           edge_alpha = NULL,

                           cluster_column = NULL,
                           select_cell_groups = NULL,
                           select_cells = NULL,
                           show_other_cells = T,
                           other_cell_color = 'lightgrey',
                           other_point_size = 1,

                           genes_high_color = NULL,
                           genes_mid_color = "white",
                           genes_low_color = "blue",

                           show_grid = FALSE,
                           spatial_grid_name = 'spatial_grid',
                           point_size = 2,
                           show_legend = TRUE,
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

  # Set feat_type and spat_unit
  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = spat_unit)
  feat_type = set_default_feat_type(gobject = gobject,
                                    spat_unit = spat_unit,
                                    feat_type = feat_type)

  # data.table variables
  cell_ID = NULL

  selected_genes = genes

  values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))
  expr_values = get_expression_values(gobject = gobject,
                                      spat_unit = spat_unit,
                                      feat_type = feat_type,
                                      values = values)

  # only keep genes that are in the dataset
  selected_genes = selected_genes[selected_genes %in% rownames(expr_values) ]

  # get selected gene expression values in data.table format
  if(length(selected_genes) == 1) {
    subset_expr_data = expr_values[rownames(expr_values) %in% selected_genes, ]
    t_sub_expr_data_DT = data.table::data.table('selected_gene' = subset_expr_data, 'cell_ID' = colnames(expr_values))
    data.table::setnames(t_sub_expr_data_DT, 'selected_gene', selected_genes)
  } else {
    subset_expr_data = expr_values[rownames(expr_values) %in% selected_genes, ]
    t_sub_expr_data = t_flex(subset_expr_data)
    t_sub_expr_data_DT = data.table::as.data.table(as.matrix(t_sub_expr_data))
    t_sub_expr_data_DT[, cell_ID := rownames(t_sub_expr_data)]
  }


  ## extract cell locations
  cell_locations  = get_spatial_locations(gobject = gobject,
                                          spat_unit = spat_unit,
                                          spat_loc_name = spat_loc_name)


  ## extract spatial network
  if(show_network == TRUE) {
    spatial_network = get_spatialNetwork(gobject,
                                         spat_unit = spat_unit,
                                         name = spatial_network_name,
                                         return_network_Obj = FALSE)
  } else {
    spatial_network = NULL
  }

  ## extract spatial grid
  if(show_grid == TRUE) {
    spatial_grid = get_spatialGrid(gobject,
                                   spat_unit = spat_unit,
                                   feat_type = feat_type,
                                   spatial_grid_name)
  } else {
    spatial_grid = NULL
  }

  ## extract cell metadata
  cell_metadata = pDataDT(gobject,
                          feat_type = feat_type,
                          spat_unit = spat_unit)
  cell_metadata = cell_metadata[, !grepl('cell_ID', colnames(cell_metadata)), with = F]


  if(nrow(cell_metadata) == 0) {
    cell_locations_metadata = cell_locations
  } else {
    cell_locations_metadata = cbind(cell_locations, cell_metadata)
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
    spatial_network = spatial_network[spatial_network$to %in% select_cells & spatial_network$from %in% select_cells]

    # if specific cells are selected
    cell_locations_metadata = cell_locations_metadata_selected
  }

  cell_locations_metadata_genes = merge(cell_locations_metadata, t_sub_expr_data_DT, by = 'cell_ID')



  ## plotting ##
  axis_scale = match.arg(axis_scale, c("cube", "real", "custom"))

  ratio = plotly_axis_scale_3D(cell_locations_metadata_genes,
                               sdimx = "sdimx",sdimy = "sdimy",sdimz = "sdimz",
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
#' @param feat_type feature type
#' @param spat_unit spatial unit
#' @param expression_values gene expression values to use
#' @param genes genes to show
#' @param dim_reduction_to_use dimension reduction to use
#' @param dim_reduction_name dimension reduction name
#' @param dim1_to_use dimension to use on x-axis
#' @param dim2_to_use dimension to use on y-axis
#' @param dim3_to_use dimension to use on z-axis
#'
#' @param show_NN_network show underlying NN network
#' @param nn_network_to_use type of NN network to use (kNN vs sNN)
#' @param network_name name of NN network to use, if show_NN_network = TRUE
#' @param network_color color of NN network
#'
#' @param cluster_column cluster column to select groups
#' @param select_cell_groups select subset of cells/clusters based on cell_color parameter
#' @param select_cells select subset of cells based on cell IDs
#' @param show_other_cells display not selected cells
#' @param other_cell_color color of not selected cells
#' @param other_point_size size of not selected cells
#'
#' @param edge_alpha column to use for alpha of the edges
#' @param point_size size of point (cell)
#'
#' @param genes_high_color color for high expression levels
#' @param genes_mid_color color for medium expression levels
#' @param genes_low_color color for low expression levels
#'
#' @param show_legend show legend
#' @param show_plot show plots
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters, see \code{\link{showSaveParameters}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return ggplot
#' @details Description of parameters.
#' @family dimension reduction gene expression visualizations
#' @export
dimGenePlot3D <- function(gobject,
                          feat_type = NULL,
                          spat_unit = NULL,
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

  # Set feat_type and spat_unit
  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = spat_unit)
  feat_type = set_default_feat_type(gobject = gobject,
                                    spat_unit = spat_unit,
                                    feat_type = feat_type)

  ## select genes ##
  selected_genes = genes
  values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))
  expr_values = get_expression_values(gobject = gobject,
                                      spat_unit = spat_unit,
                                      feat_type = feat_type,
                                      values = values)

  # only keep genes that are in the dataset
  selected_genes = selected_genes[selected_genes %in% rownames(expr_values) ]

  #
  if(length(selected_genes) == 1) {
    subset_expr_data = expr_values[rownames(expr_values) %in% selected_genes, ]
    t_sub_expr_data_DT = data.table('selected_gene' = subset_expr_data, 'cell_ID' = colnames(expr_values))
    data.table::setnames(t_sub_expr_data_DT, 'selected_gene', selected_genes)
  } else {
    subset_expr_data = expr_values[rownames(expr_values) %in% selected_genes, ]
    t_sub_expr_data = t_flex(subset_expr_data)
    t_sub_expr_data_DT = data.table::as.data.table(as.matrix(t_sub_expr_data))

    # data.table variables
    cell_ID = NULL

    t_sub_expr_data_DT[, cell_ID := rownames(t_sub_expr_data)]
  }


  ## dimension reduction ##
  dim_dfr = select_dimReduction(gobject,
                                reduction = 'cells',
                                reduction_method = dim_reduction_to_use,
                                name = dim_reduction_name,
                                return_dimObj = FALSE)
  dim_dfr = dim_dfr[,c(dim1_to_use, dim2_to_use, dim3_to_use)]
  dim_names = colnames(dim_dfr)
  dim_DT = data.table::as.data.table(dim_dfr)
  dim_DT[, cell_ID := rownames(dim_dfr)]

  ## annotated cell metadata
  cell_metadata = pDataDT(gobject,
                          feat_type = feat_type,
                          spat_unit = spat_unit)
  annotated_DT = merge(cell_metadata, dim_DT, by = 'cell_ID')



  # create input for network
  if(show_NN_network == TRUE) {

    # nn_network
    selected_nn_network = get_NearestNetwork(gobject = gobject,
                                             feat_type = feat_type,
                                             spat_unit = spat_unit,
                                             nn_network_to_use = nn_network_to_use,
                                             network_name = network_name,
                                             output = 'igraph')
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
#' @param feat_type feature type
#' @param spat_unit spatial unit
#' @param spat_loc_name name of spatial locations to use
#' @param expression_values gene expression values to use
#' @param plot_alignment direction to align plot
#' @param dim_reduction_to_use dimension reduction to use
#' @param dim_reduction_name dimension reduction name
#' @param dim1_to_use dimension to use on x-axis
#' @param dim2_to_use dimension to use on y-axis
#' @param dim3_to_use dimension to use on z-axis
#' @param sdimx spatial dimension to use on x-axis
#' @param sdimy spatial dimension to use on y-axis
#' @param sdimz spatial dimension to use on z-axis
#' @param genes genes to show
#'
#' @param cluster_column cluster column to select groups
#' @param select_cell_groups select subset of cells/clusters based on cell_color parameter
#' @param select_cells select subset of cells based on cell IDs
#' @param show_other_cells display not selected cells
#' @param other_cell_color color of not selected cells
#' @param other_point_size size of not selected cells
#'
#' @param dim_point_size dim reduction plot: point size
#' @param show_NN_network show underlying NN network
#' @param nn_network_to_use type of NN network to use (kNN vs sNN)
#' @param nn_network_color color of NN network
#' @param nn_network_alpha alpha of NN network
#' @param network_name name of NN network to use, if show_NN_network = TRUE
#'
#' @param label_size size of labels
#' @param genes_high_color color for high expression levels
#' @param genes_mid_color color for medium expression levels
#' @param genes_low_color color for low expression levels
#'
#' @param show_spatial_network show spatial network (boolean)
#' @param spatial_network_name name of spatial network to use
#' @param spatial_network_color color of spatial network
#' @param spatial_network_alpha alpha of spatial network
#'
#' @param show_spatial_grid show spatial grid (boolean)
#' @param spatial_grid_name name of spatial grid to use
#' @param spatial_grid_color color of spatial grid
#' @param spatial_grid_alpha alpha of spatial grid
#'
#' @param spatial_point_size spatial plot: point size
#' @param legend_text_size size of legend
#'
#' @param axis_scale the way to scale the axis
#' @param custom_ratio customize the scale of the plot
#' @param x_ticks set the number of ticks on the x-axis
#' @param y_ticks set the number of ticks on the y-axis
#' @param z_ticks set the number of ticks on the z-axis
#'
#' @param show_plot show plots
#' @param return_plot return plotly object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters, see \code{\link{showSaveParameters}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return plotly
#' @details Description of parameters.
#' @family spatial and dimension reduction gene expression visualizations
#' @export
spatDimGenePlot3D <- function(gobject,
                              feat_type = NULL,
                              spat_unit = NULL,
                              spat_loc_name = NULL,
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

                              show_NN_network = FALSE,
                              nn_network_to_use = 'sNN',
                              nn_network_color = 'lightgrey',
                              nn_network_alpha = 0.5,
                              network_name = 'sNN.pca',
                              label_size = 16,
                              genes_low_color = "blue",
                              genes_mid_color = "white",
                              genes_high_color = "red",
                              dim_point_size = 3,

                              show_spatial_network = FALSE,
                              spatial_network_name = 'Delaunay_network',
                              spatial_network_color = "lightgray",
                              spatial_network_alpha = 0.5,

                              show_spatial_grid = FALSE,
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

  # Set feat_type and spat_unit
  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = spat_unit)
  feat_type = set_default_feat_type(gobject = gobject,
                                    spat_unit = spat_unit,
                                    feat_type = feat_type)

  # data.table variables
  cell_ID = NULL

  plot_alignment = match.arg(plot_alignment, choices = c( 'horizontal','vertical'))

  # ********data prepare********#
  ## select genes ##
  if(length(genes) > 1){
    warning("\n Now 3D mode can just accept one gene, only the first gene will be plot\n")
    genes = genes[1]
  }
  selected_genes = genes
  values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))
  expr_values = get_expression_values(gobject = gobject,
                                      spat_unit = spat_unit,
                                      feat_type = feat_type,
                                      values = values)

  # only keep genes that are in the dataset
  selected_genes = selected_genes[selected_genes %in% rownames(expr_values) ]
  subset_expr_data = expr_values[rownames(expr_values) %in% selected_genes, ]
  t_sub_expr_data_DT = data.table('selected_gene' = subset_expr_data, 'cell_ID' = colnames(expr_values))
  data.table::setnames(t_sub_expr_data_DT, 'selected_gene', selected_genes)


  ## dimension reduction ##
  dim_dfr = select_dimReduction(gobject,
                                reduction = 'cells',
                                reduction_method = dim_reduction_to_use,
                                name = dim_reduction_name,
                                return_dimObj = FALSE)
  dim_dfr = dim_dfr[,c(dim1_to_use, dim2_to_use, dim3_to_use)]
  dim_names = colnames(dim_dfr)
  dim_DT = data.table::as.data.table(dim_dfr)
  dim_DT[, cell_ID := rownames(dim_dfr)]


  ## annotated cell metadata
  cell_metadata = pDataDT(gobject,
                          feat_type = feat_type,
                          spat_unit = spat_unit)
  cell_locations = get_spatial_locations(gobject = gobject,
                                         spat_unit = spat_unit,
                                         spat_loc_name = spat_loc_name)
  annotated_DT = merge(cell_metadata, dim_DT, by = 'cell_ID')
  annotated_DT = merge(annotated_DT, cell_locations, by = 'cell_ID')
  annotated_DT = merge(annotated_DT, t_sub_expr_data_DT,by = 'cell_ID')


  ## nn network
  if(show_NN_network){

    # nn_network
    selected_nn_network = get_NearestNetwork(gobject = gobject,
                                             feat_type = feat_type,
                                             spat_unit = spat_unit,
                                             nn_network_to_use = nn_network_to_use,
                                             network_name = network_name,
                                             output = 'igraph')
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
  if(show_spatial_network == TRUE) {
    spatial_network = get_spatialNetwork(gobject,
                                         spat_unit = spat_unit,
                                         name = spatial_network_name,
                                         return_network_Obj = FALSE)
  } else {
    spatial_network = NULL
  }

  ## extract spatial grid
  if(show_spatial_grid == TRUE) {
    spatial_grid = get_spatialGrid(gobject,
                                   spat_unit = spat_unit,
                                   feat_type = feat_type,
                                   spatial_grid_name)
  } else {
    spatial_grid = NULL
  }


  ## select subset of cells ##
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
                                          line = list(color = nn_network_color,
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
                                       line=list(color = nn_network_color),
                                       opacity= nn_network_alpha)
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
                                            line = list(color = spatial_network_color,
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
                                         line=list(color = spatial_network_color),
                                         opacity = spatial_network_alpha)
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


# ** ####
# ** interactive plotting ####

#' Select image regions by plotting interactive polygons
#'
#' @description Plot interactive polygons on an image and retrieve the polygons coordinates.
#' @param x A `ggplot` or `rast` plot object to draw polygons on
#' @param width,height An integer, defining the width/height in pixels.
#' @param ... Graphical parameters passed on to `polygon` or `geom_point`.
#'
#' @return A `data.table` containing x,y coordinates from the plotted polygons.
#'
#' @import shiny miniUI data.table magrittr ggplot2
#'
#' @examples
#' \dontrun{
#' # Using a ggplot2 plot
#' library(ggplot2)
#' df <- data.frame(x = 1:5, y = 1:5)
#' my_plot <- ggplot(df, aes(x,y)) + geom_point()
#' plotInteractivePolygons(my_plot)
#'
#' # Using a terra rast image
#' library(terra)
#' r = rast(system.file("ex/elev.tif", package="terra"))
#' plotInteractivePolygons(r)
#' plotInteractivePolygons(r, border = "red", lwd = 2)
#'
#' # Using an image contained in Giotto object
#' library(Giotto)
#' my_spatPlot <- spatPlot2D(gobject = my_giotto_object,
#'                           show_image = TRUE,
#'                           point_alpha = 0.75,
#'                           save_plot = FALSE)
#' plotInteractivePolygons(my_spatPlot, height = 500)
#' my_polygon_coordinates <- plotInteractivePolygons(my_spatPlot, height = 500)
#' }
#' @export
plotInteractivePolygons <- function(x, width = "auto", height = "auto", ...) {

  # data.table variables
  y = name = NULL

  if(is.null(x)) stop('plot object is empty')

  ui <- miniPage(
    gadgetTitleBar("Plot Interactive Polygons"),
    miniContentPanel(
      textInput("polygon_name", label = "Polygon name", value = "polygon 1"),
      sliderInput("xrange", label = "x coordinates",
                  min = min(x[["layers"]][[1]]$data$sdimx),
                  max = max(x[["layers"]][[1]]$data$sdimx),
                  value = c(min(x[["layers"]][[1]]$data$sdimx),
                            max(x[["layers"]][[1]]$data$sdimx))) ,
      sliderInput("yrange", label = "y coordinates",
                  min = min(x[["layers"]][[1]]$data$sdimy),
                  max = max(x[["layers"]][[1]]$data$sdimy),
                  value = c(min(x[["layers"]][[1]]$data$sdimy),
                            max(x[["layers"]][[1]]$data$sdimy))) ,
      plotOutput("plot", click = "plot_click")
    )
  )

  server <- function(input, output,session) {
    output$plot <- renderPlot({
      if ("ggplot" %in% class(x)) {
        x$coordinates$default <- TRUE
        x +
          geom_polygon(data = clicklist(), aes(x,y, color = name, fill = name),
                       alpha = 0, ...) +
          coord_fixed(xlim = c(input$xrange[1], input$xrange[2]),
                      ylim = c(input$yrange[1], input$yrange[2])) +
          theme(legend.position = 'none')
      } else {
        terra::plot(x)
        lapply(split(clicklist(), by = "name"), function (x) graphics::polygon(x$x, x$y, ...) )
      }
    }, res = 96, width = width, height = height)

    clicklist <- reactiveVal(data.table(x = numeric(), y = numeric(), name = character())) # empty table
    observeEvent(input$plot_click, {
      click_x <- input$plot_click$x
      click_y <- input$plot_click$y
      polygon_name <- input$polygon_name
      temp <- clicklist() # get the table of past clicks
      temp <- rbind(temp,data.table(x = click_x, y = click_y, name = polygon_name))
      clicklist(temp)
    })


    output$info <- renderTable(clicklist())

    observeEvent(input$done, {
      returnValue <- clicklist()
      stopApp(returnValue)
    })
  }
  runGadget(ui, server)
}


#' Get cells located within the polygons area
#'
#' @param gobject A Giotto object
#' @param polygon_slot Slot name where polygon coordinates are stored in Giotto object
#' @param cells_loc_slot Slot name where cell coordinates are stored in Giotto object
#'
#' @return A `SpatVector` with cell IDs, x,y coordinates, and polygon name where
#'  each cell is located in.
#'
#' @export
#'
#' @examples
#'
#' \dontrun{
#' ## Plot interactive polygons
#' my_spatPlot <- spatPlot2D(gobject = my_giotto_object,
#'                           show_image = TRUE,
#'                           point_alpha = 0.75,
#'                           save_plot = FALSE)
#' my_polygon_coords <- plotInteractivePolygons(my_spatPlot)
#'
#' ## Add polygon coordinates to Giotto object
#' my_giotto_polygons <- createGiottoPolygonsFromDfr(my_polygon_coords)
#' my_giotto_object <- addGiottoPolygons(gobject = my_giotto_object,
#'                                       gpolygons = list(my_giotto_polygons))
#'
#' ## Get cells located within polygons area
#' getCellsFromPolygon(my_giotto_object)
#' }
#'

getCellsFromPolygon <- function(gobject,
                                polygon_slot = "spatial_info",
                                cells_loc_slot = "spatial_locs") {
  ## verify Giotto object
  if (!inherits(gobject, "giotto")) {
    stop("gobject needs to be a giotto object")
  }

  ## get polygon spatvector
  my_polygon_spatplot <- methods::slot(methods::slot(gobject, polygon_slot)$cell,"spatVector")

  ## get spatial locs from cells
  my_spatial_locs <- methods::slot(gobject, cells_loc_slot)$cell$raw

  ## create spatvector from spatial locs
  my_cells_spatplot <- terra::vect(as.matrix(my_spatial_locs[,1:2]),
                                   type = "points",
                                   atts = my_spatial_locs)

  ## get the insersect of polygons and cells
  my_intersect <- terra::intersect(my_cells_spatplot, my_polygon_spatplot)

  return(my_intersect)
}


#' Add corresponding polygon IDs to cell metadata
#'
#' @param gobject A Giotto object
#' @param cellsFromPolygon A `SpatVector` with cell IDs located inside each polygon
#' @param feat_type feature name where metadata will be added
#'
#' @return A Giotto object with a modified cell_metadata slot that includes the
#' polygon name where each cell is located or NA if the cell is not located
#' within a polygon area
#'
#' @export
#'
#' @examples
#'
#' \dontrun{
#' ## Plot interactive polygons
#' my_polygon_coords <- plotInteractivePolygons(my_spatPlot)
#'
#' ## Add polygon coordinates to Giotto object
#' my_giotto_polygons <- createGiottoPolygonsFromDfr(my_polygon_coords)
#' my_giotto_object <- addGiottoPolygons(gobject = my_giotto_object,
#'                                       gpolygons = list(my_giotto_polygons))
#'
#' ## Get cells located within polygons area
#' my_polygon_cells <- getCellsFromPolygon(my_giotto_object)
#'
#' my_giotto_object <- addCellsFromPolygon(my_giotto_object, my_polygon_cells)
#' }
#'

addCellsFromPolygon <- function(gobject,
                                cellsFromPolygon,
                                feat_type = "rna") {

  ## verify Giotto object
  if (!inherits(gobject, "giotto")) {
    stop("gobject needs to be a giotto object")
  }

  ## get original metadas
  cell_metadata <- methods::slot(gobject, "cell_metadata")$cell[[feat_type]]

  ## convert cellsFromPolygon to data frame
  cellsFromPolygondata <- as.data.frame(cellsFromPolygon)

  ## merge metadata as data.table
  new_cell_metadata <- merge(cell_metadata, cellsFromPolygondata[,c("cell_ID", "poly_ID")],
                             by = "cell_ID", all.x = TRUE)

  ## replace NA's with a defult name
  new_cell_metadata[is.na(new_cell_metadata$poly_ID),"poly_ID"] <- "no_polygon"

  ## set rows order
  new_cell_metadata <- new_cell_metadata[match(cell_metadata$cell_ID, new_cell_metadata$cell_ID),]

  ## add cell metadata to Giotto object
  addCellMetadata(gobject = gobject, new_metadata = new_cell_metadata[,-1])

}


#' Compare gene expression between polygon areas
#'
#' @param gobject A Giotto object
#' @param spat_unit spatial unit (e.g. "cell")
#' @param feat_type feature type (e.g. "rna", "dna", "protein")
#' @param selected_feats vector of selected features to plot
#' @param expression_values gene expression values to use ("normalized", "scaled", "custom")
#' @param method method to use to detect differentially expressed feats ("scran", "gini", "mast")
#' @param \dots Arguments passed to \link[ComplexHeatmap]{Heatmap}
#'
#' @return A ComplexHeatmap::Heatmap object
#' @export
comparePolygonExpression <- function(gobject,
                                     spat_unit = "cell",
                                     feat_type = "rna",
                                     selected_feats = "top_genes",
                                     expression_values = "normalized",
                                     method = "scran",
                                     ...) {

  # verify gobject
  if (!inherits(gobject, "giotto")) {
    stop("gobject needs to be a giotto object")
  }

  # get expression
  my_expression <- gobject@expression[[spat_unit]][[feat_type]][[expression_values]]

  # get cell_ID and poly_ID from metadata
  my_metadata <- gobject@cell_metadata[[spat_unit]][[feat_type]]
  my_metadata <- my_metadata[,c("cell_ID", "poly_ID")]

  if (length(selected_feats) == 1 && selected_feats == "top_genes") {
    # find top features
    scran_results <- findMarkers_one_vs_all(gobject,
                                            spat_unit = "cell",
                                            feat_type = "rna",
                                            method = method,
                                            expression_values = "normalized",
                                            cluster_column = "poly_ID",
                                            min_feats = 10)

    selected_feats <- scran_results[, head(.SD, 2), by = 'cluster']$feats
  }
  # select features
  my_expression <- my_expression[selected_feats,]

  # convert to data frame
  my_rownames <- rownames(my_expression)

  # calculate zscore

  my_zscores <- my_expression

  for (gene in my_rownames) {
    mean_expression_gene <- mean(my_expression[gene,])
    sd_expression_gene <- stats::sd(my_expression[gene,])
    for (cell in colnames(my_expression)) {
      my_zscores[gene, cell] <- (my_expression[gene, cell]-mean_expression_gene)/sd_expression_gene
    }

  }

  # calculate mean zscore per polygon
  my_zscores_mean <- data.table::data.table(feat_ID = my_rownames)

  for(i in unique(my_metadata$poly_ID)) {
    my_cells <- my_metadata[my_metadata$poly_ID == i, "cell_ID" ]
    my_sub_zscores <- my_zscores[,my_cells$cell_ID]
    mean_zscores <- Matrix::rowMeans(my_sub_zscores)
    my_zscores_mean <- cbind(my_zscores_mean, mean_zscores)
  }

  my_zscores_mean <- as.matrix(my_zscores_mean[,-1])
  colnames(my_zscores_mean) <- unique(my_metadata$poly_ID)
  rownames(my_zscores_mean) <- my_rownames

  # plot heatmap
  my_heatmap <- ComplexHeatmap::Heatmap(my_zscores_mean,
                                        heatmap_legend_param = list(title = "Normalized mean z score"),
                                        cluster_rows = FALSE,
                                        cluster_columns = FALSE,
                                        column_order = sort(colnames(my_zscores_mean)),
                                        ...)
  return(my_heatmap)
}

