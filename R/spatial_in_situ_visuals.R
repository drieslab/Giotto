




#' @title Plot cell polygon layer
#' @name plot_cell_polygon_layer
#' @description Low level function to plot a polygon
#' @return ggplot
#' @details This functions plots a polygon based on spatial cell information.
#' This is most likely a polygon that corresponds to the cell shape.
#' @keywords internal
plot_cell_polygon_layer = function(ggobject = NULL,
                                   polygon_dt,
                                   polygon_grouping = 'poly_ID',
                                   sdimx = 'x',
                                   sdimy = 'y',
                                   fill = NULL,
                                   poly_fill_gradient = c('blue', 'white', 'red'),
                                   fill_gradient_midpoint =  NULL,
                                   fill_as_factor = TRUE,
                                   fill_code = NULL,
                                   bg_color = 'black',
                                   color = 'black',
                                   alpha = 0.5,
                                   size = 2) {

  # data.table variables
  final_fill = NULL

  # check fill column
  if(!is.null(fill)) {
    if(fill_as_factor == TRUE) {
      polygon_dt[, final_fill := as.factor(get(fill))]
    } else {
      polygon_dt[, final_fill := get(fill)]
    }
  }

  # create layer
  if(!is.null(ggobject) & methods::is(ggobject, 'ggplot')) {
    pl = ggobject
  } else {
    pl = ggplot2::ggplot()
  }


  # specific fill color for polygon shapes
  if(!is.null(fill)) {
    pl = pl + ggplot2::geom_polygon(data = polygon_dt,
                                    ggplot2::aes_string(x = 'x',
                                                        y = 'y',
                                                        group = polygon_grouping,
                                                        fill = 'final_fill'),
                                    alpha = alpha,
                                    color = color,
                                    size = size)

    # manual fill colors for factor values
    if(fill_as_factor == TRUE) {
      if(!is.null(fill_code)) {
        pl = pl + ggplot2::scale_fill_manual(values = fill_code)
      } else {
        fill_values_names = unique(polygon_dt[['final_fill']])
        fill_code = getDistinctColors(length(fill_values_names))
        names(fill_code) = fill_values_names
        pl = pl + ggplot2::scale_fill_manual(values = fill_code)
      }
    }

    # gradient fill colors for numerical values
    if(fill_as_factor == FALSE) {

      if(is.null(fill_gradient_midpoint)) {
        fill_gradient_midpoint = stats::median(polygon_dt[['final_fill']])
      }

      pl = pl + ggplot2::scale_fill_gradient2(low = poly_fill_gradient[[1]],
                                              mid = poly_fill_gradient[[2]],
                                              high = poly_fill_gradient[[3]],
                                              midpoint = fill_gradient_midpoint,
                                              guide = guide_colorbar(title = ''))
    }


  } else {
    pl = pl + ggplot2::geom_polygon(data = polygon_dt,
                                    ggplot2::aes_string(x = 'x',
                                                        y = 'y',
                                                        group = 'poly_ID'),
                                    fill = bg_color,
                                    alpha = alpha,
                                    color = color,
                                    size = size)
  }

  return(pl)

}



#' @title select_gimage
#' @name select_gimage
#' @description selects and creates giotto images for plotting
#' @keywords internal
select_gimage = function(gobject,
                         gimage = NULL,
                         image_name = NULL,
                         largeImage_name = NULL,
                         spat_unit = NULL,
                         spat_loc_name = NULL,
                         feat_type = NULL,
                         polygon_feat_type = NULL) {


  if(!is.null(gimage)) gimage = gimage


  else if(!is.null(image_name)) {

    if(length(image_name) == 1) {
      gimage = gobject@images[[image_name]]
      if(is.null(gimage)) warning('image_name: ', image_name, ' does not exists')
    } else {
      gimage = list()
      for(gim in 1:length(image_name)) {
        gimage[[gim]] = gobject@images[[gim]]
        if(is.null(gimage[[gim]])) warning('image_name: ', gim, ' does not exists')
      }
    }

  } else if(!is.null(largeImage_name)) {
    # If there is input to largeImage_name arg

    if(length(largeImage_name) == 1) {
      gimage = plot_auto_largeImage_resample(gobject = gobject,
                                             largeImage_name = largeImage_name,
                                             spat_unit = spat_unit,
                                             spat_loc_name = spat_loc_name,
                                             polygon_feat_type = polygon_feat_type,
                                             include_image_in_border = TRUE)
    } else {
      gimage = list()
      for(gim in 1:length(largeImage_name)) {
        gimage[[gim]] = plot_auto_largeImage_resample(gobject = gobject,
                                                      largeImage_name = largeImage_name[[gim]],
                                                      spat_unit = spat_unit,
                                                      spat_loc_name = spat_loc_name,
                                                      polygon_feat_type = polygon_feat_type,
                                                      include_image_in_border = TRUE)
      }
    }

  } else {
    # Default to first image available in images if no input given to image_name or largeImage_name args
    image_name = names(gobject@images)[1]
    gimage = gobject@images[[image_name]]

    if(is.null(gimage)) warning('image_name: ', image_name, ' does not exist \n')
  }

  return(gimage)

}




#' @title expand_feature_info
#' @name expand_feature_info
#' @description low level function to expand feature coordinates
#' @return data.table
#' @keywords internal
expand_feature_info = function(spatial_feat_info,
                               expand_counts = FALSE,
                               count_info_column = 'count',
                               jitter = c(0,0),
                               verbose = TRUE) {


  # 1. expand feature locations with multiple counts (e.g. in seq-Scope or Stereo-seq)
  if(isTRUE(expand_counts)) {

    if(!count_info_column %in% colnames(spatial_feat_info)) stop('count_info_column ', count_info_column, ' does not exist')

    if(isTRUE(verbose)) {wrap_msg('Start expanding feature information based on count column')}
    #print(spatial_feat_info)

    extra_feats = spatial_feat_info[get(count_info_column) > 1]
    extra_feats = extra_feats[,rep(get(count_info_column), get(count_info_column)), by = .(feat_ID, x, y, feat, spat_unit)]
    spatial_feat_info = rbind(extra_feats[,.(feat_ID, x, y, feat, spat_unit)], spatial_feat_info[get(count_info_column) == 1, .(feat_ID, x, y, feat, spat_unit)])

    #print(spatial_feat_info)

  }

  # 2. add jitter to x and y coordinates

  if(!identical(c(0,0), jitter)) {

    if(isTRUE(verbose)) {wrap_msg('Start adding jitter to x and y based on provided max jitter information')}

    # create jitter for x and y coordinates: from 0 to max-x or max-y
    tx_number = nrow(spatial_feat_info)
    x_jitter = sample(0:jitter[[1]], size = tx_number, replace = TRUE)
    y_jitter = sample(0:jitter[[2]], size = tx_number, replace = TRUE)

    spatial_feat_info[, c('x', 'y') := list(x+x_jitter, y+y_jitter)]
  }

  return(spatial_feat_info)

}



#' @title plot_feature_points_layer
#' @name plot_feature_points_layer
#' @description low level function to plot a points at the spatial in situ level
#' @return ggplot
#' @details This function can plot multiple features over multiple modalities. These plots can get very big very fast.
#' @keywords internal
plot_feature_points_layer = function(ggobject,
                                     spatial_feat_info,
                                     feats,
                                     feats_color_code = NULL,
                                     feat_shape_code = NULL,
                                     sdimx = 'x',
                                     sdimy = 'y',
                                     color = 'feat_ID',
                                     shape = 'feat',
                                     point_size = 1.5,
                                     show_legend = TRUE,
                                     plot_method = c('ggplot', 'scattermore', 'scattermost'),
                                     expand_counts = FALSE,
                                     count_info_column = 'count',
                                     jitter = c(0,0),
                                     verbose = TRUE) {

  # data.table variables
  feat_ID = NULL

  spatial_feat_info_subset = spatial_feat_info[feat_ID %in% unlist(feats)]

  # expand feature coordinates and/or add jitter to coordiantes
  if(isTRUE(expand_counts) | !identical(c(0,0), jitter)) {
    spatial_feat_info_subset = expand_feature_info(spatial_feat_info = spatial_feat_info_subset,
                                                   expand_counts = expand_counts,
                                                   count_info_column = count_info_column,
                                                   jitter = jitter,
                                                   verbose = verbose)
  }

  cat(' --| Plotting ', nrow(spatial_feat_info_subset), ' feature points\n')

  if(!is.null(ggobject) & inherits(ggobject, 'ggplot')) {
    pl = ggobject
  } else {
    pl = ggplot2::ggplot()
  }

  pl = pl + giotto_point(plot_method = plot_method,
                         data = spatial_feat_info_subset,
                         ggplot2::aes_string(x = sdimx,
                                             y = sdimy,
                                             color = color,
                                             shape = shape),
                         size = point_size,
                         show.legend = show_legend)



  if(!is.null(feats_color_code)) {
    pl = pl + ggplot2::scale_color_manual(values = feats_color_code)
  } else {
    feats_names = unique(spatial_feat_info_subset[[color]])
    feats_color_code = getDistinctColors(length(feats_names))
    names(feats_color_code) = feats_names
    pl = pl + ggplot2::scale_color_manual(values = feats_color_code)
  }

  return(pl)
}


#' @title spatInSituPlotPoints
#' @name spatInSituPlotPoints
#' @description Function to plot multiple features for multiple modalities at the spatial in situ level
#' @param gobject giotto object
#' @param show_image show a tissue background image
#' @param gimage a giotto image
#' @param image_name name of a giotto image
#' @param largeImage_name name of a giottoLargeImage
#' @param spat_unit spatial unit
#' @param spat_loc_name name of spatial locations
#' @param feats features to plot
#' @param feat_type feature types of the feats
#' @param feats_color_code code to color the provided features
#' @param feat_shape_code code to shape the provided feature types
#' @param sdimx spatial dimension x
#' @param sdimy spatial dimension y
#' @param point_size size of the points
#' @param expand_counts expand feature coordinate counts (see details)
#' @param count_info_column column name with count information (if expand_counts = TRUE)
#' @param jitter maximum x,y jitter provided as c(x, y)
#' @param show_polygon overlay polygon information (e.g. cell shape)
#' @param use_overlap use polygon and feature coordinates overlap results
#' @param polygon_feat_type feature type associated with polygon information
#' @param polygon_color color for polygon border
#' @param polygon_bg_color color for polygon background (overruled by polygon_fill)
#' @param polygon_fill fill color or column for polygon
#' @param polygon_fill_gradient polygon fill gradient colors given in order from low to high
#' @param polygon_fill_gradient_midpoint value to set as gradient midpoint (optional). If
#'   left as \code{NULL}, the median value detected will be chosen
#' @param polygon_fill_as_factor is fill color a factor
#' @param polygon_fill_code code to color the fill column
#' @param polygon_alpha alpha of polygon
#' @param polygon_line_size line width of the polygon's outline
#' @param axis_text axis text size
#' @param axis_title title text size
#' @param legend_text legend text size
#' @param coord_fix_ratio fix ratio of coordinates
#' @param background_color background color
#' @param show_legend show legend
#' @param plot_method method to plot points
#' @param show_plot show plots
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters, see \code{\link{showSaveParameters}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @param verbose verbosity
#' @return ggplot
#' @details TODO
#' @family In Situ visualizations
#' @export
spatInSituPlotPoints = function(gobject,
                                show_image = F,
                                gimage = NULL,
                                image_name = NULL,
                                largeImage_name = NULL,
                                spat_unit = NULL,
                                spat_loc_name = NULL,
                                feats = NULL,
                                feat_type = 'rna',
                                feats_color_code = NULL,
                                feat_shape_code = NULL,
                                sdimx = 'x',
                                sdimy = 'y',
                                point_size = 1.5,
                                expand_counts = FALSE,
                                count_info_column = 'count',
                                jitter = c(0,0),
                                show_polygon = TRUE,
                                use_overlap = TRUE,
                                polygon_feat_type = 'cell',
                                polygon_color = 'black',
                                polygon_bg_color = 'black',
                                polygon_fill = NULL,
                                polygon_fill_gradient = c('blue', 'white', 'red'),
                                polygon_fill_gradient_midpoint =  NULL,
                                polygon_fill_as_factor = NULL,
                                polygon_fill_code = NULL,
                                polygon_alpha = 0.5,
                                polygon_line_size = 2,
                                axis_text = 8,
                                axis_title = 8,
                                legend_text = 6,
                                coord_fix_ratio = 1,
                                background_color = 'black',
                                show_legend = TRUE,
                                plot_method = c('ggplot', 'scattermore', 'scattermost'),
                                show_plot = NA,
                                return_plot = NA,
                                save_plot = NA,
                                save_param =  list(),
                                default_save_name = 'spatInSituPlotPoints',
                                verbose = TRUE) {


  if(is.null(feats)) {
    warning('You need to select features (feats) and modify feature types (feat_type) if you want to show individual features (e.g. transcripts) \n')
  }


  # print, return and save parameters
  show_plot = ifelse(is.na(show_plot), readGiottoInstructions(gobject, param = 'show_plot'), show_plot)
  save_plot = ifelse(is.na(save_plot), readGiottoInstructions(gobject, param = 'save_plot'), save_plot)
  return_plot = ifelse(is.na(return_plot), readGiottoInstructions(gobject, param = 'return_plot'), return_plot)


  ## giotto image ##
  if(show_image == TRUE) {

    gimage = select_gimage(gobject = gobject,
                           gimage = gimage,
                           image_name = image_name,
                           largeImage_name = largeImage_name,
                           spat_unit = spat_unit,
                           spat_loc_name = spat_loc_name,
                           feat_type = feat_type,
                           polygon_feat_type = polygon_feat_type)

    if(isTRUE(verbose)) wrap_msg('select image done')

  }






  # start plotting
  plot = ggplot2::ggplot()

  ## 0. plot image ##
  if(show_image == TRUE & !is.null(gimage)) {
    plot = plot_spat_image_layer_ggplot(gg_obj = plot,
                                        gobject = gobject,
                                        spat_unit = spat_unit,
                                        feat_type = feat_type,
                                        spat_loc_name = spat_loc_name,
                                        polygon_feat_type = polygon_feat_type,
                                        gimage = gimage,
                                        sdimx = 'sdimx',
                                        sdimy = 'sdimy')

    if(isTRUE(verbose)) wrap_msg('plot image layer done')
  }


  ## 1. plot morphology first
  if(show_polygon == TRUE) {

    # Set feat_type and spat_unit
    polygon_feat_type = set_default_spat_unit(gobject = gobject,
                                      spat_unit = polygon_feat_type)
    feat_type = set_default_feat_type(gobject = gobject,
                                      spat_unit = polygon_feat_type,
                                      feat_type = feat_type)

    #feat_type = set_default_feat_type(gobject = gobject, feat_type = feat_type)
    #if(is.null(polygon_feat_type)) {
    #  polygon_feat_type = gobject@expression_feat[[1]]
    #}

    polygon_combo = combineCellData(gobject = gobject,
                                    spat_loc_name = spat_loc_name,
                                    feat_type = feat_type,
                                    include_poly_info = TRUE,
                                    poly_info = polygon_feat_type)
    polygon_dt = polygon_combo[[feat_type]]
    data.table::setnames(polygon_dt, old = 'cell_ID', new = 'poly_ID')

    #polygon_info = get_polygon_info(gobject = gobject,
    #                                polygon_name = polygon_feat_type)
    #polygon_dt = spatVector_to_dt(polygon_info)

    plot = plot_cell_polygon_layer(ggobject = plot,
                                   polygon_dt = polygon_dt,
                                   polygon_grouping = 'poly_ID',
                                   sdimx = sdimx,
                                   sdimy = sdimy,
                                   fill = polygon_fill,
                                   poly_fill_gradient = polygon_fill_gradient,
                                   fill_gradient_midpoint = polygon_fill_gradient_midpoint,
                                   fill_as_factor = polygon_fill_as_factor,
                                   fill_code = polygon_fill_code,
                                   bg_color = polygon_bg_color,
                                   color = polygon_color,
                                   alpha = polygon_alpha,
                                   size = polygon_line_size)

    if(isTRUE(verbose)) wrap_msg('plot polygon layer done')


  }


  ## 2. plot features second

  if(!is.null(feats)) {
    # use_overlap = TRUE will use the overlap results
    # use_overlap = FALSE will use the raw tx coordinate results
    if(use_overlap == TRUE) {

      # TODO: check if overlap exists, if not print warning message and default to non-overlap results
      spatial_feat_info = combineFeatureOverlapData(gobject = gobject,
                                                    feat_type = feat_type,
                                                    sel_feats = feats,
                                                    poly_info = polygon_feat_type)
    } else {

      spatial_feat_info = combineFeatureData(gobject = gobject,
                                             spat_unit =  polygon_feat_type,
                                             feat_type = feat_type,
                                             sel_feats = feats)
    }

    spatial_feat_info = do.call('rbind', spatial_feat_info)

    plot = plot_feature_points_layer(ggobject = plot,
                                     spatial_feat_info = spatial_feat_info,
                                     feats = feats,
                                     feats_color_code = feats_color_code,
                                     feat_shape_code = feat_shape_code,
                                     expand_counts = expand_counts,
                                     count_info_column = count_info_column,
                                     jitter = jitter,
                                     sdimx = 'x',
                                     sdimy = 'y',
                                     color = 'feat_ID',
                                     shape = 'feat',
                                     point_size = point_size,
                                     show_legend = show_legend,
                                     plot_method = plot_method)

    if(isTRUE(verbose)) wrap_msg('plot feature points layer done')

  }


  ## 3. adjust theme settings
  plot <- plot + ggplot2::theme(plot.title = element_text(hjust = 0.5),
                                legend.title = element_blank(),
                                legend.text = element_text(size = legend_text),
                                axis.title = element_text(size = axis_title),
                                axis.text = element_text(size = axis_text),
                                panel.grid = element_blank(),
                                panel.background = element_rect(fill = background_color))



  if(!is.null(coord_fix_ratio)) {
    plot = plot + ggplot2::coord_fixed(ratio = coord_fix_ratio)
  }


  ## print plot
  if(show_plot == TRUE) {
    print(plot)
  }

  ## save plot
  if(save_plot == TRUE) {
    do.call('all_plots_save_function', c(list(gobject = gobject,
                                              plot_object = plot,
                                              default_save_name = default_save_name),
                                         save_param))
  }

  ## return plot
  if(return_plot == TRUE) {
    return(plot)
  }

}



#' @title plot_feature_hexbin_layer
#' @name plot_feature_hexbin_layer
#' @description low level function to plot hexbins at the spatial in situ level
#' @return ggplot
#' @details This function can plot one feature for one modality.
#' @keywords internal
plot_feature_hexbin_layer = function(ggobject = NULL,
                                     spatial_feat_info,
                                     sel_feat,
                                     sdimx = 'x',
                                     sdimy = 'y',
                                     binwidth = NULL,
                                     min_axis_bins = 10L,
                                     alpha = 0.5) {

  # data.table variables
  feat_ID = NULL

  spatial_feat_info_subset = spatial_feat_info[feat_ID %in% sel_feat]

  # set default binwidth to 1/10 of minor axis
  if(is.null(binwidth)) {
    minorRange = spatial_feat_info_subset[, min(diff(sapply(.SD, range))), .SDcols = c('x','y')]
    binwidth = as.integer(minorRange/min_axis_bins)
  }

  if(!is.null(ggobject) & methods::is(ggobject, 'ggplot')) {
    pl = ggobject
  } else {
    pl = ggplot2::ggplot()
  }

  pl = pl + ggplot2::geom_hex(data = spatial_feat_info_subset,
                              ggplot2::aes_string(x = sdimx,
                                                  y = sdimy),
                              binwidth = binwidth,
                              alpha = alpha)
  pl = pl + labs(title = sel_feat)
  return(pl)

}



#' @title spatInSituPlotHex_single
#' @name spatInSituPlotHex_single
#' @description function to plot hexbins at the spatial in situ level
#' @return ggplot
#' @details This function can plot one feature for one modality.
#' @keywords internal
spatInSituPlotHex_single = function(gobject,
                                    feat = NULL,
                                    feat_type = 'rna',
                                    sdimx = 'x',
                                    sdimy = 'y',
                                    binwidth = NULL,
                                    min_axis_bins = NULL,
                                    alpha = 0.5,
                                    show_polygon = TRUE,
                                    polygon_feat_type = 'cell',
                                    polygon_color = 'black',
                                    polygon_fill = NULL,
                                    polygon_fill_as_factor = NULL,
                                    polygon_alpha = 0.5,
                                    polygon_size = 0.5,
                                    coord_fix_ratio = NULL,
                                    axis_text = 8,
                                    axis_title = 8,
                                    legend_text = 6,
                                    background_color = 'black') {


  if(is.null(feat)) {
    stop('You need to select a feature (feat) and modify feature types (feat_type) if needed \n')
  }

  plot = ggplot2::ggplot()

  ## polygon layer ##
  if(show_polygon == TRUE) {

    if(is.null(polygon_feat_type)) {
      polygon_feat_type = gobject@expression_feat[[1]]
    }


    #polygon_dt = combineSpatialCellMetadataInfo(gobject, feat_type = polygon_feat_type)
    #polygon_dt = polygon_dt[[polygon_feat_type]]

    polygon_info = get_polygon_info(gobject = gobject,
                                       polygon_name = polygon_feat_type)
    polygon_dt = spatVector_to_dt(polygon_info)

    plot = plot_cell_polygon_layer(ggobject = gobject,
                                   polygon_dt,
                                   polygon_grouping = 'poly_ID',
                                   sdimx = sdimx,
                                   sdimy = sdimy,
                                   fill = polygon_fill,
                                   fill_as_factor = polygon_fill_as_factor,
                                   color = polygon_color,
                                   alpha = polygon_alpha,
                                   size = polygon_size)

  }



  ## hexbin layer ##

  form_feat = list(feat_type = c(feat))
  spatial_feat_info = combineFeatureOverlapData(gobject = gobject,
                                                feat_type = feat_type,
                                                sel_feats = form_feat,
                                                poly_info = polygon_feat_type)

  #spatial_feat_info = combineSpatialCellFeatureInfo(gobject = gobject,
  #                                                  feat_type = feat_type,
  #                                                  selected_features = feat)
  spatial_feat_info = do.call('rbind', spatial_feat_info)

  plot = plot_feature_hexbin_layer(ggobject = plot,
                                   spatial_feat_info = spatial_feat_info,
                                   sel_feat = feat,
                                   sdimx = sdimx,
                                   sdimy = sdimy,
                                   binwidth = binwidth,
                                   min_axis_bins = min_axis_bins)


  ## adjust theme settings
  plot <- plot + ggplot2::theme(plot.title = element_text(hjust = 0.5),
                                legend.title = element_blank(),
                                legend.text = element_text(size = legend_text),
                                axis.title = element_text(size = axis_title),
                                axis.text = element_text(size = axis_text),
                                panel.grid = element_blank(),
                                panel.background = element_rect(fill = background_color))

  # fix coord ratio
  if(!is.null(coord_fix_ratio)) {
    plot = plot + ggplot2::coord_fixed(ratio = coord_fix_ratio)
  }

  return(plot)

}



#' @title spatInSituPlotHex
#' @name spatInSituPlotHex
#' @description Function to plot hexbins for features for multiple modalities at the spatial in situ level
#' @param gobject giotto object
#' @param feats features to plot
#' @param feat_type feature types of the feats
#' @param sdimx spatial dimension x
#' @param sdimy spatial dimension y
#' @param binwidth numeric vector for x and y width of bins (default is minor axis
#' range/10, where the 10 is from \code{min_axis_bins})
#' @param min_axis_bins number of bins to create per range defined by minor axis.
#' (default value is 10)
#' @param alpha alpha of hexbin plot
#' @param show_polygon overlay polygon information (cell shape)
#' @param polygon_feat_type feature type associated with polygon information
#' @param polygon_color color for polygon border
#' @param polygon_fill fill color or column for polygon
#' @param polygon_fill_as_factor is fill color a factor
#' @param polygon_alpha alpha of polygon
#' @param polygon_size size of polygon border
#' @param coord_fix_ratio fix ratio between x and y-axis
#' @param axis_text axis text size
#' @param axis_title title text size
#' @param legend_text legend text size
#' @param background_color background color
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
#' @details TODO
#' @family In Situ visualizations
#' @export
spatInSituPlotHex = function(gobject,
                             feats = NULL,
                             feat_type = 'rna',
                             sdimx = 'x',
                             sdimy = 'y',
                             binwidth = NULL,
                             min_axis_bins = 10,
                             alpha = 0.5,
                             show_polygon = TRUE,
                             polygon_feat_type = 'cell',
                             polygon_color = 'black',
                             polygon_fill = NULL,
                             polygon_fill_as_factor = NULL,
                             polygon_alpha = 0.5,
                             polygon_size = 0.5,
                             coord_fix_ratio = 1,
                             axis_text = 8,
                             axis_title = 8,
                             legend_text = 6,
                             background_color = 'white',
                             cow_n_col = NULL,
                             cow_rel_h = 1,
                             cow_rel_w = 1,
                             cow_align = 'h',
                             show_plot = NA,
                             return_plot = NA,
                             save_plot = NA,
                             save_param =  list(),
                             default_save_name = 'spatInSituPlotHex') {


  if(is.null(feats)) {
    stop('You need to select features (feats) and modify feature types (feat_type) if needed \n')
  }

  # print, return and save parameters
  show_plot = ifelse(is.na(show_plot), readGiottoInstructions(gobject, param = 'show_plot'), show_plot)
  save_plot = ifelse(is.na(save_plot), readGiottoInstructions(gobject, param = 'save_plot'), save_plot)
  return_plot = ifelse(is.na(return_plot), readGiottoInstructions(gobject, param = 'return_plot'), return_plot)

  ## plotting ##
  savelist <- list()

  for(sel_feat in feats) {

    pl = spatInSituPlotHex_single(gobject = gobject,
                                  feat = sel_feat,
                                  feat_type = feat_type,
                                  sdimx = sdimx,
                                  sdimy = sdimy,
                                  binwidth = binwidth,
                                  min_axis_bins = min_axis_bins,
                                  alpha = alpha,
                                  show_polygon = show_polygon,
                                  polygon_feat_type = polygon_feat_type,
                                  polygon_color = polygon_color,
                                  polygon_fill = polygon_fill,
                                  polygon_fill_as_factor = polygon_fill_as_factor,
                                  polygon_alpha = polygon_alpha,
                                  polygon_size = polygon_size,
                                  coord_fix_ratio = coord_fix_ratio,
                                  axis_text = axis_text,
                                  axis_title = axis_title,
                                  legend_text = legend_text,
                                  background_color = background_color)

    savelist[[sel_feat]] = pl

  }

  if(length(savelist) == 1) {
    combo_plot = savelist[[1]]
  } else {
    # combine plots with cowplot
    combo_plot <- cowplot::plot_grid(plotlist = savelist,
                                     ncol = set_default_cow_n_col(cow_n_col = cow_n_col,
                                                                  nr_plots = length(savelist)),
                                     rel_heights = cow_rel_h,
                                     rel_widths = cow_rel_w,
                                     align = cow_align)
  }


  ## print plot
  if(show_plot == TRUE) {
    print(combo_plot)
  }

  ## save plot
  if(save_plot == TRUE) {
    do.call('all_plots_save_function', c(list(gobject = gobject,
                                              plot_object = combo_plot,
                                              default_save_name = default_save_name),
                                         save_param))
  }

  ## return plot
  if(return_plot == TRUE) {
    return(combo_plot)
  }


}



#' @title plot_feature_raster_density_layer
#' @name plot_feature_raster_density_layer
#' @description low level function to plot density plots at the spatial in situ level
#' @return ggplot
#' @details This function can plot one feature for one modality.
#' @keywords internal
plot_feature_raster_density_layer = function(ggobject = NULL,
                                              spatial_feat_info,
                                              sel_feat,
                                              sdimx = 'x',
                                              sdimy = 'y',
                                              alpha = 0.5) {

  # data.table variable
  feat_ID = NULL

  spatial_feat_info_subset = spatial_feat_info[feat_ID %in% unlist(sel_feat)]

  if(!is.null(ggobject) & methods::is(ggobject, 'ggplot')) {
    pl = ggobject
  } else {
    pl = ggplot2::ggplot()
  }

  pl = pl + ggplot2::stat_density_2d(data = spatial_feat_info_subset,
                                     ggplot2::aes_string(x = sdimx,
                                                         y = sdimy,
                                                         fill = 'after_stat(density)'),
                                     geom = "raster",
                                     alpha = alpha,
                                     contour = FALSE)
  pl = pl + labs(title = sel_feat)
  pl = pl + ggplot2::scale_fill_continuous(type = "viridis")

  return(pl)

}



#' @title spatInSituPlotDensity_single
#' @name spatInSituPlotDensity_single
#' @description low level function to plot density plots at the spatial in situ level
#' @return ggplot
#' @details This function can plot one feature for one modality.
#' @keywords internal
spatInSituPlotDensity_single = function(gobject,
                                        feat = NULL,
                                        feat_type = 'rna',
                                        sdimx = 'x',
                                        sdimy = 'y',
                                        alpha = 0.95,
                                        show_polygon = TRUE,
                                        polygon_feat_type = 'cell',
                                        polygon_color = 'black',
                                        polygon_fill = NULL,
                                        polygon_fill_as_factor = NULL,
                                        polygon_alpha = 0.5,
                                        polygon_size = 0.5,
                                        coord_fix_ratio = NULL,
                                        axis_text = 8,
                                        axis_title = 8,
                                        legend_text = 6,
                                        background_color = 'black') {


  if(is.null(feat)) {
    stop('You need to select a feature (feat) and modify feature types (feat_type) if needed \n')
  }

  plot = ggplot2::ggplot()

  ## polygon layer ##
  if(show_polygon == TRUE) {

    if(is.null(polygon_feat_type)) {
      polygon_feat_type = gobject@expression_feat[[1]]
    }


    polygon_info = get_polygon_info(gobject = gobject,
                                    polygon_name = polygon_feat_type)
    polygon_dt = spatVector_to_dt(polygon_info)

    #polygon_dt = combineSpatialCellMetadataInfo(gobject, feat_type = polygon_feat_type)
    #polygon_dt = polygon_dt[[polygon_feat_type]]

    plot = plot_cell_polygon_layer(ggobject = gobject,
                                   polygon_dt,
                                   polygon_grouping = 'poly_ID',
                                   sdimx = sdimx,
                                   sdimy = sdimy,
                                   fill = polygon_fill,
                                   fill_as_factor = polygon_fill_as_factor,
                                   color = polygon_color,
                                   alpha = polygon_alpha,
                                   size = polygon_size)
  }



  ## density layer ##
  form_feat = list(feat_type = c(feat))
  spatial_feat_info = combineFeatureOverlapData(gobject = gobject,
                                                feat_type = feat_type,
                                                sel_feats = form_feat,
                                                poly_info = polygon_feat_type)

  #spatial_feat_info = combineSpatialCellFeatureInfo(gobject = gobject,
  #                                                  feat_type = feat_type,
  #                                                  selected_features = feat)
  spatial_feat_info = do.call('rbind', spatial_feat_info)

  plot = plot_feature_raster_density_layer(ggobject = plot,
                                           spatial_feat_info = spatial_feat_info,
                                           sel_feat = feat,
                                           sdimx = sdimx,
                                           sdimy = sdimy,
                                           alpha = alpha)


  ## adjust theme settings
  plot <- plot + ggplot2::theme(plot.title = element_text(hjust = 0.5),
                                legend.title = element_blank(),
                                legend.text = element_text(size = legend_text),
                                axis.title = element_text(size = axis_title),
                                axis.text = element_text(size = axis_text),
                                panel.grid = element_blank(),
                                panel.background = element_rect(fill = background_color))

  # fix coord ratio
  if(!is.null(coord_fix_ratio)) {
    plot = plot + ggplot2::coord_fixed(ratio = coord_fix_ratio)
  }

  return(plot)

}



#' @title spatInSituPlotDensity
#' @name spatInSituPlotDensity
#' @description Function for density plots for features for multiple modalities at the spatial in situ level
#' @param gobject giotto object
#' @param feats features to plot
#' @param feat_type feature types of the feats
#' @param sdimx spatial dimension x
#' @param sdimy spatial dimension y
#' @param alpha alpha of density plot
#' @param show_polygon overlay polygon information (cell shape)
#' @param polygon_feat_type feature type associated with polygon information
#' @param polygon_color color for polygon border
#' @param polygon_fill fill color or column for polygon
#' @param polygon_fill_as_factor is fill color a factor
#' @param polygon_alpha alpha of polygon
#' @param polygon_size size of polygon border
#' @param coord_fix_ratio fix ratio between x and y-axis
#' @param axis_text axis text size
#' @param axis_title title text size
#' @param legend_text legend text size
#' @param background_color background color
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
#' @details TODO
#' @family In Situ visualizations
#' @export
spatInSituPlotDensity = function(gobject,
                                 feats = NULL,
                                 feat_type = 'rna',
                                 sdimx = 'x',
                                 sdimy = 'y',
                                 alpha = 0.95,
                                 show_polygon = TRUE,
                                 polygon_feat_type = 'cell',
                                 polygon_color = 'black',
                                 polygon_fill = NULL,
                                 polygon_fill_as_factor = NULL,
                                 polygon_alpha = 0.5,
                                 polygon_size = 0.5,
                                 coord_fix_ratio = 1,
                                 axis_text = 8,
                                 axis_title = 8,
                                 legend_text = 6,
                                 background_color = 'black',
                                 cow_n_col = NULL,
                                 cow_rel_h = 1,
                                 cow_rel_w = 1,
                                 cow_align = 'h',
                                 show_plot = NA,
                                 return_plot = NA,
                                 save_plot = NA,
                                 save_param =  list(),
                                 default_save_name = 'spatInSituPlotDensity') {


  if(is.null(feats)) {
    stop('You need to select features (feat) and modify feature types (feat_type) if needed \n')
  }

  # print, return and save parameters
  show_plot = ifelse(is.na(show_plot), readGiottoInstructions(gobject, param = 'show_plot'), show_plot)
  save_plot = ifelse(is.na(save_plot), readGiottoInstructions(gobject, param = 'save_plot'), save_plot)
  return_plot = ifelse(is.na(return_plot), readGiottoInstructions(gobject, param = 'return_plot'), return_plot)




  ## plotting ##
  savelist <- list()

  for(sel_feat in feats) {

    pl = spatInSituPlotDensity_single(gobject = gobject,
                                      feat = sel_feat,
                                      feat_type = feat_type,
                                      sdimx = sdimx,
                                      sdimy = sdimy,
                                      alpha = alpha,
                                      show_polygon = show_polygon,
                                      polygon_feat_type = polygon_feat_type,
                                      polygon_color = polygon_color,
                                      polygon_fill = polygon_fill,
                                      polygon_fill_as_factor = polygon_fill_as_factor,
                                      polygon_alpha = polygon_alpha,
                                      polygon_size = polygon_size,
                                      coord_fix_ratio = coord_fix_ratio,
                                      axis_text = axis_text,
                                      axis_title = axis_title,
                                      legend_text = legend_text,
                                      background_color = background_color)

    savelist[[sel_feat]] = pl

  }

  if(length(savelist) == 1) {
    combo_plot = savelist[[1]]
  } else {
    # combine plots with cowplot
    combo_plot <- cowplot::plot_grid(plotlist = savelist,
                                     ncol = set_default_cow_n_col(cow_n_col = cow_n_col,
                                                                  nr_plots = length(savelist)),
                                     rel_heights = cow_rel_h,
                                     rel_widths = cow_rel_w,
                                     align = cow_align)
  }



  ## print plot
  if(show_plot == TRUE) {
    print(combo_plot)
  }

  ## save plot
  if(save_plot == TRUE) {
    do.call('all_plots_save_function', c(list(gobject = gobject,
                                              plot_object = combo_plot,
                                              default_save_name = default_save_name),
                                         save_param))
  }

  ## return plot
  if(return_plot == TRUE) {
    return(combo_plot)
  }

}

