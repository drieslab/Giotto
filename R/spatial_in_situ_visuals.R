




#' @name plot_cell_polygon_layer
#' @description low level function to plot a polygon
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
                                   fill_as_factor = TRUE,
                                   color = 'black',
                                   alpha = 0.5) {


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

  if(!is.null(fill)) {
    pl = pl + ggplot2::geom_polygon(data = polygon_dt,
                                    ggplot2::aes_string(x = sdimx,
                                                        y = sdimy,
                                                        group = polygon_grouping,
                                                        fill = 'final_fill'),
                                    alpha = alpha,
                                    color = color)
  } else {
    pl = pl + ggplot2::geom_polygon(data = polygon_dt,
                                    ggplot2::aes_string(x = sdimx,
                                                        y = sdimy,
                                                        group = polygon_grouping),
                                    fill = 'lightblue',
                                    alpha = alpha,
                                    color = color)
  }

  return(pl)

}







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
                                     plot_method = c('ggplot', 'scattermore', 'scattermost')) {


  spatial_feat_info_subset = spatial_feat_info[feat_ID %in% unlist(feats)]

  if(!is.null(ggobject) & methods::is(ggobject, 'ggplot')) {
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
                         size = point_size, show.legend = show_legend)

  if(!is.null(feats_color_code)) {
    pl = pl + ggplot2::scale_color_manual(values = feats_color_code)
  }

  return(pl)
}



#' @name spatInSituPlotPoints
#' @description Function to plot multiple features for multiple modalities at the spatial in situ level
#' @param gobject giotto object
#' @param feats features to plot
#' @param feat_type feature types of the feats
#' @param feats_color_code code to color the provided features
#' @param feat_shape_code code to shape the provided feature types
#' @param sdimx spatial dimension x
#' @param sdimy spatial dimension y
#' @param point_size size of the points
#' @param show_polygon overlay polygon information (cell shape)
#' @param polygon_feat_type feature type associated with polygon information
#' @param polygon_color color for polygon border
#' @param polygon_fill fill color or column for polygon
#' @param polygon_fill_as_factor is fill color a factor
#' @param polygon_alpha alpha of polygon
#' @param axis_text axis text size
#' @param axis_title title text size
#' @param legend_text legend text size
#' @param background_color background color
#' @param show_legend show legend
#' @param plot_method method to plot points
#' @param show_plot show plots
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters, see \code{\link{showSaveParameters}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return ggplot
#' @details TODO
#' @family In Situ visualizations
#' @export
spatInSituPlotPoints = function(gobject,
                                feats = NULL,
                                feat_type = 'rna',
                                feats_color_code = NULL,
                                feat_shape_code = NULL,
                                sdimx = 'x',
                                sdimy = 'y',
                                point_size = 1.5,
                                show_polygon = TRUE,
                                polygon_feat_type = 'cell',
                                polygon_color = 'black',
                                polygon_fill = NULL,
                                polygon_fill_as_factor = NULL,
                                polygon_alpha = 0.5,
                                axis_text = 8,
                                axis_title = 8,
                                legend_text = 6,
                                background_color = 'black',
                                show_legend = TRUE,
                                plot_method = c('ggplot', 'scattermore', 'scattermost'),
                                show_plot = NA,
                                return_plot = NA,
                                save_plot = NA,
                                save_param =  list(),
                                default_save_name = 'spatInSituPlotPoints') {


  if(is.null(feats)) {
    stop('You need to select features (feats) and modify feature types (feat_type) if needed \n')
  }


  # print, return and save parameters
  show_plot = ifelse(is.na(show_plot), readGiottoInstructions(gobject, param = 'show_plot'), show_plot)
  save_plot = ifelse(is.na(save_plot), readGiottoInstructions(gobject, param = 'save_plot'), save_plot)
  return_plot = ifelse(is.na(return_plot), readGiottoInstructions(gobject, param = 'return_plot'), return_plot)


  # start plotting
  plot = ggplot2::ggplot()

  if(show_polygon == TRUE) {

    if(is.null(polygon_feat_type)) {
      polygon_feat_type = gobject@expression_feat[[1]]
    }


    #testobj@spatial_info$cell@spatVector

    polygon_info = get_polygon_info(gobject = gobject,
                                       polygon_name = polygon_feat_type)
    polygon_dt = spatVector_to_dt(polygon_info)

    #polygon_dt = spatVector_to_dt(gobject@spatial_info[[polygon_feat_type]]@spatVector)

    #polygon_dt = combineCellData(gobject = gobject,
    #                             feat_type = polygon_feat_type)
    #polygon_dt = polygon_dt[[polygon_feat_type]]

    plot = plot_cell_polygon_layer(ggobject = gobject,
                                    polygon_dt,
                                    polygon_grouping = 'poly_ID',
                                    sdimx = sdimx,
                                    sdimy = sdimy,
                                    fill = polygon_fill,
                                    fill_as_factor = polygon_fill_as_factor,
                                    color = polygon_color,
                                    alpha = polygon_alpha)

  }

  spatial_feat_info = combineFeatureOverlapData(gobject = gobject,
                                                feat_type = feat_type,
                                                sel_feats = feats)

  #spatial_feat_info = combineSpatialCellFeatureInfo(gobject = gobject, feat_type = feat_type)
  spatial_feat_info = do.call('rbind', spatial_feat_info)

  plot = plot_feature_points_layer(ggobject = plot,
                                   spatial_feat_info = spatial_feat_info,
                                   feats = feats,
                                   feats_color_code = feats_color_code,
                                   feat_shape_code = feat_shape_code,
                                   sdimx = 'x',
                                   sdimy = 'y',
                                   color = 'feat_ID',
                                   shape = 'feat',
                                   point_size = point_size,
                                   show_legend = show_legend,
                                   plot_method = plot_method)

  ## adjust theme settings
  plot <- plot + ggplot2::theme(plot.title = element_text(hjust = 0.5),
                                legend.title = element_blank(),
                                legend.text = element_text(size = legend_text),
                                axis.title = element_text(size = axis_title),
                                axis.text = element_text(size = axis_text),
                                panel.grid = element_blank(),
                                panel.background = element_rect(fill = background_color))




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
                                      bins = 10,
                                      alpha = 0.5) {



  spatial_feat_info_subset = spatial_feat_info[feat_ID %in% sel_feat]

  if(!is.null(ggobject) & methods::is(ggobject, 'ggplot')) {
    pl = ggobject
  } else {
    pl = ggplot2::ggplot()
  }

  pl = pl + ggplot2::geom_hex(data = spatial_feat_info_subset,
                              ggplot2::aes_string(x = sdimx,
                                                  y = sdimy),
                              bins = bins,
                              alpha = alpha)
  pl = pl + labs(title = sel_feat)
  return(pl)

}



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
                                     bins = 10,
                                     alpha = 0.5,
                                     show_polygon = TRUE,
                                     polygon_feat_type = 'cell',
                                     polygon_color = 'black',
                                     polygon_fill = NULL,
                                     polygon_fill_as_factor = NULL,
                                     polygon_alpha = 0.5,
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
                                    alpha = polygon_alpha)

  }



  ## hexbin layer ##

  form_feat = list(feat_type = c(feat))
  spatial_feat_info = combineFeatureOverlapData(gobject = gobject,
                                                feat_type = feat_type,
                                                sel_feats = form_feat)

  #spatial_feat_info = combineSpatialCellFeatureInfo(gobject = gobject,
  #                                                  feat_type = feat_type,
  #                                                  selected_features = feat)
  spatial_feat_info = do.call('rbind', spatial_feat_info)

  plot = plot_feature_hexbin_layer(ggobject = plot,
                                    spatial_feat_info = spatial_feat_info,
                                    sel_feat = feat,
                                    sdimx = sdimx,
                                    sdimy = sdimy,
                                    bins = bins)


  ## adjust theme settings
  plot <- plot + ggplot2::theme(plot.title = element_text(hjust = 0.5),
                                legend.title = element_blank(),
                                legend.text = element_text(size = legend_text),
                                axis.title = element_text(size = axis_title),
                                axis.text = element_text(size = axis_text),
                                panel.grid = element_blank(),
                                panel.background = element_rect(fill = background_color))



  return(plot)

}




#' @name spatInSituPlotHex
#' @description Function to plot hexbins for features for multiple modalities at the spatial in situ level
#' @param gobject giotto object
#' @param feats features to plot
#' @param feat_type feature types of the feats
#' @param sdimx spatial dimension x
#' @param sdimy spatial dimension y
#' @param bins number of hexbins in one direction
#' @param show_polygon overlay polygon information (cell shape)
#' @param polygon_feat_type feature type associated with polygon information
#' @param polygon_color color for polygon border
#' @param polygon_fill fill color or column for polygon
#' @param polygon_fill_as_factor is fill color a factor
#' @param polygon_alpha alpha of polygon
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
                             bins = 10,
                             alpha = 0.5,
                             show_polygon = TRUE,
                             polygon_feat_type = 'cell',
                             polygon_color = 'black',
                             polygon_fill = NULL,
                             polygon_fill_as_factor = NULL,
                             polygon_alpha = 0.5,
                             axis_text = 8,
                             axis_title = 8,
                             legend_text = 6,
                             background_color = 'white',
                             cow_n_col = 2,
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
                                   bins = bins,
                                   alpha = alpha,
                                   show_polygon = show_polygon,
                                   polygon_feat_type = polygon_feat_type,
                                   polygon_color = polygon_color,
                                   polygon_fill = polygon_fill,
                                   polygon_fill_as_factor = polygon_fill_as_factor,
                                   polygon_alpha = polygon_alpha,
                                   axis_text = axis_text,
                                   axis_title = axis_title,
                                   legend_text = legend_text,
                                   background_color = background_color)

    savelist[[sel_feat]] = pl

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

  spatial_feat_info_subset = spatial_feat_info[feat_ID %in% unlist(sel_feat)]

  if(!is.null(ggobject) & methods::is(ggobject, 'ggplot')) {
    pl = ggobject
  } else {
    pl = ggplot2::ggplot()
  }

  pl = pl + ggplot2::stat_density_2d(data = spatial_feat_info_subset,
                                     ggplot2::aes_string(x = sdimx,
                                                         y = sdimy,
                                                         fill = '..density..'),
                                     geom = "raster",
                                     alpha = alpha,
                                     contour = FALSE)
  pl = pl + labs(title = sel_feat)
  pl = pl + ggplot2::scale_fill_continuous(type = "viridis")

  return(pl)

}



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
                                    alpha = polygon_alpha)
  }



  ## density layer ##
  form_feat = list(feat_type = c(feat))
  spatial_feat_info = combineFeatureOverlapData(gobject = gobject,
                                                feat_type = feat_type,
                                                sel_feats = form_feat)

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



  return(plot)

}




#' @name spatInSituPlotDensity
#' @description Function for density plots for features for multiple modalities at the spatial in situ level
#' @param gobject giotto object
#' @param feats features to plot
#' @param feat_type feature types of the feats
#' @param sdimx spatial dimension x
#' @param sdimy spatial dimension y
#' @param show_polygon overlay polygon information (cell shape)
#' @param polygon_feat_type feature type associated with polygon information
#' @param polygon_color color for polygon border
#' @param polygon_fill fill color or column for polygon
#' @param polygon_fill_as_factor is fill color a factor
#' @param polygon_alpha alpha of polygon
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
                                  axis_text = 8,
                                  axis_title = 8,
                                  legend_text = 6,
                                  background_color = 'black',
                                  cow_n_col = 2,
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
                                       axis_text = axis_text,
                                       axis_title = axis_title,
                                       legend_text = legend_text,
                                       background_color = background_color)

    savelist[[sel_feat]] = pl

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

