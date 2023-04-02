
#' @title aes_string2
#' @name aes_string2
#' @param \dots aes_string parameters
#' @keywords internal
#' @description makes sure aes_string can also be used with names that start with numeric values
aes_string2 <- function(...){
  args <- lapply(list(...), function(x) sprintf("`%s`", x))
  do.call(ggplot2::aes_string, args)
}


#' @title giotto_point
#' @name giotto_point
#' @param \dots geom_point parameters
#' @keywords internal
#' @description uses ggplot::geom_point, scattermore::geom_scattermore or scattermore::geom_scattermost
giotto_point = function(plot_method = c('ggplot', 'scattermore', 'scattermost'),
                        ...) {

  plot_method = match.arg(arg = plot_method, choices = c('ggplot', 'scattermore', 'scattermost'))

  if(plot_method == 'ggplot') {
    ggplot2::geom_point(...)
  } else if(plot_method == 'scattermore') {
    package_check(pkg_name = "scattermore",
                  repository = 'CRAN')
    scattermore::geom_scattermore(...)
  } else if(plot_method == 'scattermost') {
    package_check(pkg_name = "scattermore",
                  repository = 'CRAN')
    scattermore::geom_scattermost(...)
  }

}







#' @title ggplot_save_function
#' @name ggplot_save_function
#' @description Function to automatically save plots to directory of interest
#' @param gobject giotto object
#' @param plot_object ggplot object to plot
#' @param save_dir directory to save to
#' @param save_folder folder in save_dir to save to
#' @param save_name name of plot
#' @param save_format format (e.g. png, tiff, pdf, ...)
#' @param show_saved_plot load & display the saved plot
#' @param ncol number of columns
#' @param nrow number of rows
#' @param scale scale
#' @param base_width width
#' @param base_height height
#' @param base_aspect_ratio aspect ratio
#' @param units units
#' @param dpi Plot resolution
#' @param limitsize When TRUE (the default), ggsave will not save images larger than 50x50 inches, to prevent the common error of specifying dimensions in pixels.
#' @param plot_count count number for plot
#' @param \dots additional parameters to cowplot::save_plot
#' @seealso \code{\link[cowplot]{save_plot}}
#' @keywords internal
ggplot_save_function = function(gobject,
                                plot_object,
                                save_dir = NULL,
                                save_folder = NULL,
                                save_name = NULL,
                                default_save_name = 'giotto_plot',
                                save_format = NULL,
                                show_saved_plot = F,
                                ncol = 1,
                                nrow = 1,
                                scale = 1,
                                base_width = NULL,
                                base_height = NULL,
                                base_aspect_ratio = NULL,
                                units = NULL,
                                dpi = NULL,
                                limitsize = TRUE,
                                plot_count = NULL,
                                ...) {

  if(is.null(plot_object)) {
    stop('\t there is no object to plot \t')
  }

  ## get save information and set defaults
  if(is.null(save_dir)) save_dir = readGiottoInstructions(gobject, param = 'save_dir')
  if(is.null(save_folder)) save_folder = NULL
  if(is.null(plot_count)) plot_count = getOption('giotto.plot_count')
  if(is.null(save_name)) {
    save_name = default_save_name
    save_name = paste0(plot_count,'-', save_name)
    options('giotto.plot_count' = plot_count + 1)
  }
  if(is.null(save_format)) save_format = readGiottoInstructions(gobject, param = 'plot_format')
  if(is.null(dpi)) dpi = readGiottoInstructions(gobject, param = 'dpi')
  if(is.null(base_width)) base_width = readGiottoInstructions(gobject, param = 'width')
  if(is.null(base_height)) base_height = readGiottoInstructions(gobject, param = 'height')
  if(is.null(base_aspect_ratio)) base_aspect_ratio = 1.1
  if(is.null(units)) units = 'in'

  ## checking
  dpi = as.numeric(dpi)
  base_width = as.numeric(base_width)
  base_height = as.numeric(base_height)
  base_aspect_ratio = as.numeric(base_aspect_ratio)

  # create saving location
  if(!is.null(save_folder)) {
    file_location = paste0(save_dir,'/', save_folder)
  } else {
    file_location = save_dir
  }
  if(!file.exists(file_location)) dir.create(file_location, recursive = T)
  file_name = paste0(save_name, ".", save_format)

  cowplot::save_plot(plot = plot_object,
                     filename = file_name,
                     path = file_location,
                     device = save_format,
                     ncol = ncol,
                     nrow = nrow,
                     scale = scale,
                     base_width = base_width,
                     base_height = base_height,
                     base_aspect_ratio = base_aspect_ratio,
                     units = units,
                     dpi = dpi,
                     limitsize = limitsize,
                     ...)

  # show saved plot if requested
  if(show_saved_plot == TRUE) {

    if(save_format == 'png') {
      if("png" %in% rownames(installed.packages()) == FALSE) {
        cat("\n package 'png' is not yet installed \n")
      } else {
        img = png::readPNG(source = paste0(file_location, '/', file_name))
        grid::grid.raster(img)
        }

    } else if(save_format == 'tiff') {
      if("tiff" %in% rownames(installed.packages()) == FALSE) {
        cat("\n package 'tiff' is not yet installed \n")
      } else {
        img = tiff::readTIFF(source =  paste0(file_location, '/', file_name))
        grid::grid.raster(img)
      }
    } else {
      cat('\t only png & tiff are currently supported \t')
    }
  }
}



#' @title general_save_function
#' @name general_save_function
#' @description Function to automatically save plots to directory of interest
#' @param gobject giotto object
#' @param plot_object non-ggplot object to plot
#' @param save_dir directory to save to
#' @param save_folder folder in save_dir to save to
#' @param save_name name of plot
#' @param save_format format (e.g. png, tiff, pdf, ...)
#' @param show_saved_plot load & display the saved plot
#' @param base_width width
#' @param base_height height
#' @param base_aspect_ratio aspect ratio
#' @param units units
#' @param dpi Plot resolution
#' @param plot_count count number for plot
#' @keywords internal
general_save_function = function(gobject,
                                 plot_object,
                                 save_dir = NULL,
                                 save_folder = NULL,
                                 save_name = NULL,
                                 default_save_name = 'giotto_plot',
                                 save_format = c('png', 'tiff', 'pdf', 'svg'),
                                 show_saved_plot = F,
                                 base_width = NULL,
                                 base_height = NULL,
                                 base_aspect_ratio = NULL,
                                 units = NULL,
                                 dpi = NULL,
                                 plot_count = NULL,
                                 ...) {




  if(is.null(plot_object)) {
    stop('\t there is no object to plot \t')
  }
  save_format = match.arg(save_format, choices = c('png', 'tiff', 'pdf', 'svg'))

  if(any('plotly' %in% class(plot_object)) == TRUE){
    save_format = "html"
  }

  ## get save information and set defaults
  if(is.null(save_dir)) save_dir = readGiottoInstructions(gobject, param = 'save_dir')
  if(is.null(save_folder)) save_folder = NULL
  if(is.null(plot_count)) plot_count = getOption('giotto.plot_count')
  if(is.null(save_name)) {
    save_name = default_save_name
    save_name = paste0(plot_count,'-', save_name)
    options('giotto.plot_count' = plot_count + 1)
  }
  if(is.null(save_format)) save_format = readGiottoInstructions(gobject, param = 'plot_format')
  if(is.null(dpi)) dpi = readGiottoInstructions(gobject, param = 'dpi')
  if(is.null(base_width)) base_width = readGiottoInstructions(gobject, param = 'width')
  if(is.null(base_height)) base_height = readGiottoInstructions(gobject, param = 'height')
  if(is.null(base_aspect_ratio)) base_aspect_ratio = 1.1
  if(is.null(units)) units = 'px'

  ## checking
  dpi = as.numeric(dpi)
  base_width = as.numeric(base_width)
  base_height = as.numeric(base_height)
  base_aspect_ratio = as.numeric(base_aspect_ratio)

  # create saving location
  if(!is.null(save_folder)) {
    file_location = paste0(save_dir,'/', save_folder)
  } else {
    file_location = save_dir
  }
  if(!file.exists(file_location)) dir.create(file_location, recursive = T)
  file_name = paste0(save_name, ".", save_format)
  full_location = paste0(file_location,'/', file_name)

  if(any('plotly' %in% class(plot_object)) == TRUE){
    htmlwidgets::saveWidget(plotly::as_widget(plot_object), file = full_location)
  }

  else{
    if(save_format == 'png') {
      grDevices::png(filename = full_location, width = base_width, height = base_height, res = dpi, units = units, ...)
      print(plot_object)
      grDevices::dev.off()
    }

    if(save_format == 'tiff') {
      grDevices::tiff(filename = full_location, width = base_width, height = base_height, units = units, ...)
      print(plot_object)
      grDevices::dev.off()
    }

    if(save_format == 'pdf') {
      grDevices::pdf(file = full_location, width = base_width, height = base_height, useDingbats = F, ...)
      print(plot_object)
      grDevices::dev.off()
    }

    if(save_format == 'svg') {
      grDevices::svg(filename = full_location, width = base_width, height = base_height, ...)
      print(plot_object)
      grDevices::dev.off()
    }


    # show saved plot if requested
    if(show_saved_plot == TRUE) {

      if(save_format == 'png') {
        if("png" %in% rownames(installed.packages()) == FALSE) {
          cat("\n package 'png' is not yet installed \n")
        } else {
          img <- png::readPNG(source = paste0(file_location, '/', file_name))
          grid::grid.raster(img)
        }

      } else if(save_format == 'tiff') {
        if("tiff" %in% rownames(installed.packages()) == FALSE) {
          cat("\n package 'tiff' is not yet installed \n")
        } else {
          img <- tiff::readTIFF(source =  paste0(file_location, '/', file_name))
          grid::grid.raster(img)
        }
      } else {
        cat('\t only png & tiff are currently supported \t')
      }
    }

  }

}

#' @title all_plots_save_function
#' @name all_plots_save_function
#' @description Function to automatically save plots to directory of interest
#' @param gobject giotto object
#' @param plot_object object to plot
#' @param save_dir directory to save to
#' @param save_folder folder in save_dir to save to
#' @param save_name name of plot
#' @param default_save_name default name to save a plot
#' @param save_format format (e.g. png, tiff, pdf, ...)
#' @param show_saved_plot load & display the saved plot
#' @param ncol number of columns
#' @param nrow number of rows
#' @param scale scale
#' @param base_width width
#' @param base_height height
#' @param base_aspect_ratio aspect ratio
#' @param units units
#' @param dpi Plot resolution
#' @param limitsize When TRUE (the default), ggsave will not save images larger than 50x50 inches, to prevent the common error of specifying dimensions in pixels.
#' @param plot_count count number for plot
#' @param \dots additional parameters to ggplot_save_function or general_save_function
#' @seealso \code{\link{general_save_function}}
#' @keywords internal
all_plots_save_function = function(gobject,
                                   plot_object,
                                   save_dir = NULL,
                                   save_folder = NULL,
                                   save_name = NULL,
                                   default_save_name = 'giotto_plot',
                                   save_format = NULL,
                                   show_saved_plot = F,
                                   ncol = 1,
                                   nrow = 1,
                                   scale = 1,
                                   base_width = NULL,
                                   base_height = NULL,
                                   base_aspect_ratio = NULL,
                                   units = NULL,
                                   dpi = NULL,
                                   limitsize = TRUE,
                                   plot_count = NULL,
                                   ...) {


  if(any('ggplot' %in% class(plot_object)) == TRUE) {

    ggplot_save_function(gobject = gobject,
                         plot_object = plot_object,
                         save_dir = save_dir,
                         save_folder = save_folder,
                         save_name = save_name,
                         default_save_name = default_save_name,
                         save_format = save_format,
                         show_saved_plot = show_saved_plot,
                         ncol = ncol,
                         nrow = nrow,
                         scale = scale,
                         base_width = base_width,
                         base_height = base_height,
                         base_aspect_ratio = base_aspect_ratio,
                         units = units,
                         dpi = dpi,
                         limitsize = limitsize,
                         plot_count = plot_count,
                         ...)

  } else {

    general_save_function(gobject = gobject,
                          plot_object = plot_object,
                          save_dir = save_dir,
                          save_folder = save_folder,
                          save_name = save_name,
                          default_save_name = default_save_name,
                          save_format = save_format,
                          show_saved_plot = show_saved_plot,
                          base_width = base_width,
                          base_height = base_height,
                          base_aspect_ratio = base_aspect_ratio,
                          units = units,
                          dpi = dpi,
                          plot_count = plot_count,
                          ...)

  }

}




#' @title showSaveParameters
#' @name showSaveParameters
#' @description Description of Giotto saving options, links to \code{\link{all_plots_save_function}}
#' @return Instruction on how to use the automatic plot saving options within Giotto
#' @export
#' @examples
#' \dontrun{
#'   showSaveParameters()
#' }
showSaveParameters = function() {

  cat("This is a simple guide to help you with automatically saving plots. \n")
  cat("Importantly, defaults for all these parameters can be set at the beginning with createGiottoInstructions() \n")
  cat("See https://rubd.github.io/Giotto/articles/instructions_and_plotting.html for more information and examples \n \n")

  cat("Each plotting function in Giotto has 4 important parameters for showing and/or saving a plot: \n
      - show_plot: TRUE or FALSE, show the plot to the console
      - return_plot: TRUE or FALSE, return the plot to the console (e.g. to further modify or save the plot
      - save_plot: TRUE or FALSE, automatically save the plot
      - save_param: a list of parameters that can be set \n")

  cat('\n')

  cat("The following list of parameters can be provided to save_param: \n
      - save_dir: directory to save the plot to
      - save_folder: if not NULL, a subfolder within save_dir that will be created to save the plot to
      - save_name: name of the plot (no extension needed, see save_format)
      - save_format: picture format to use, default is .png
      - ncol: number of columns for multiplots
      - nrow: number of rows for multiplot
      - scale: scale of plots
      - base_width: width of plot
      - base_height: height of plot
      - base_aspect_ratio: ratio of plot
      - units: plotting units (e.g. in)
      - dpi: dpi for each plot if plot is in raster format\n")

  cat('\n')

  cat("Example: \n
      plotfunction(...,
                   save_plot = TRUE,
                   save_param = list(save_name = 'favorite_name', units = 'png'))")

}




#' @title showClusterHeatmap
#' @name showClusterHeatmap
#' @description Creates heatmap based on identified clusters
#' @param gobject giotto object
#' @param spat_unit spatial unit (e.g. "cell")
#' @param feat_type feature type (e.g. "rna", "dna", "protein")
#' @param expression_values expression values to use (e.g. "normalized", "scaled", "custom")
#' @param feats vector of features to use, default to 'all'
#' @param genes deprecated. Replaced by \code{feats} param
#' @param cluster_column name of column to use for clusters (e.g. "leiden_clus")
#' @param cor correlation score to calculate distance (e.g. "pearson", "spearman")
#' @param distance distance method to use for hierarchical clustering, default to "ward.D"
#' @param show_plot show plot. TRUE or FALSE
#' @param return_plot return ggplot object. TRUE or FALSE
#' @param save_plot directly save the plot. TRUE or FALSE
#' @param save_param list of saving parameters, see \code{\link{showSaveParameters}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @param ... additional parameters passed to \code{\link[ComplexHeatmap]{Heatmap}} function
#' @return ggplot
#' @details Correlation heatmap of selected clusters.
#' @export
showClusterHeatmap <- function(gobject,
                               spat_unit = NULL,
                               feat_type = NULL,
                               expression_values = c('normalized', 'scaled', 'custom'),
                               feats = 'all',
                               genes = NULL,
                               cluster_column,
                               cor = c('pearson', 'spearman'),
                               distance = 'ward.D',
                               show_plot = NA,
                               return_plot = NA,
                               save_plot = NA,
                               save_param =  list(),
                               default_save_name = 'showClusterHeatmap',
                               ...) {


  ## deprecated arguments
  if(!is.null(genes)) {
    feats = genes
    warning('genes is deprecated, use feats in the future \n')
  }


  # Set feat_type and spat_unit
  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = spat_unit)
  feat_type = set_default_feat_type(gobject = gobject,
                                    spat_unit = spat_unit,
                                    feat_type = feat_type)


  ## correlation
  cor = match.arg(cor, c('pearson', 'spearman'))
  values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))

  ## subset expression data
  if(feats[1] != 'all') {
    all_feats = gobject@feat_ID[[feat_type]]
    detected_feats = feats[feats %in% all_feats]
  } else {
    # NULL = all feats in calculateMetaTable()
    detected_feats = NULL
  }

  metatable = calculateMetaTable(gobject = gobject,
                                 spat_unit = spat_unit,
                                 feat_type = feat_type,
                                 expression_values = values,
                                 metadata_cols = cluster_column,
                                 selected_feats = detected_feats)
  dcast_metatable = data.table::dcast.data.table(metatable, formula = variable~uniq_ID, value.var = 'value')
  testmatrix = dt_to_matrix(x = dcast_metatable)

  # correlation
  cormatrix = cor_flex(x = testmatrix, method = cor)
  cordist = stats::as.dist(1 - cormatrix, diag = T, upper = T)
  corclus = stats::hclust(d = cordist, method = distance)

  hmap = ComplexHeatmap::Heatmap(matrix = cormatrix,
                                 cluster_rows = corclus,
                                 cluster_columns = corclus,
                                 ...)


  # print, return and save parameters
  show_plot = ifelse(is.na(show_plot), readGiottoInstructions(gobject, param = 'show_plot'), show_plot)
  save_plot = ifelse(is.na(save_plot), readGiottoInstructions(gobject, param = 'save_plot'), save_plot)
  return_plot = ifelse(is.na(return_plot), readGiottoInstructions(gobject, param = 'return_plot'), return_plot)

  ## print plot
  if(show_plot == TRUE) {
    print(hmap)
  }

  ## save plot
  if(save_plot == TRUE) {
    do.call('all_plots_save_function', c(list(gobject = gobject, plot_object = hmap, default_save_name = default_save_name), save_param))
  }

  ## return plot
  if(return_plot == TRUE) {
    return(hmap)
  }
}



#' @title showClusterDendrogram
#' @name showClusterDendrogram
#' @description Creates dendrogram for selected clusters.
#' @param gobject giotto object
#' @param spat_unit spatial unit (e.g. "cell")
#' @param feat_type feature type (e.g. "rna", "dna", "protein")
#' @param expression_values expression values to use (e.g. "normalized", "scaled", "custom")
#' @param cluster_column name of column to use for clusters (e.g. "leiden_clus")
#' @param cor correlation score to calculate distance (e.g. "pearson", "spearman")
#' @param distance distance method to use for hierarchical clustering, default to "ward.D"
#' @param h height of horizontal lines to plot
#' @param h_color color of horizontal lines
#' @param rotate rotate dendrogram 90 degrees
#' @param show_plot show plot. TRUE or FALSE
#' @param return_plot return ggplot object. TRUE or FALSE
#' @param save_plot directly save the plot. TRUE or FALSE
#' @param save_param list of saving parameters, see \code{\link{showSaveParameters}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @param ... additional parameters passed to \code{\link[ggdendro]{ggdendrogram}}
#' @return ggplot
#' @details Expression correlation dendrogram for selected clusters.
#' @export
showClusterDendrogram <- function(gobject,
                                  spat_unit = NULL,
                                  feat_type = NULL,
                                  expression_values = c('normalized', 'scaled', 'custom'),
                                  cluster_column,
                                  cor = c('pearson', 'spearman'),
                                  distance = 'ward.D',
                                  h = NULL,
                                  h_color = 'red',
                                  rotate = FALSE,
                                  show_plot = NA,
                                  return_plot = NA,
                                  save_plot = NA,
                                  save_param =  list(),
                                  default_save_name = 'showClusterDendrogram',
                                  ...) {


  # verify if optional package is installed
  package_check(pkg_name = "ggdendro", repository = "CRAN")

  cor = match.arg(cor, c('pearson', 'spearman'))
  values = match.arg(expression_values, unique(c('normalized', 'scaled', 'custom', expression_values)))

  # Set feat_type and spat_unit
  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = spat_unit)
  feat_type = set_default_feat_type(gobject = gobject,
                                    spat_unit = spat_unit,
                                    feat_type = feat_type)

  metatable = calculateMetaTable(gobject = gobject,
                                 spat_unit = spat_unit,
                                 feat_type = feat_type,
                                 expression_values = values,
                                 metadata_cols = cluster_column)
  dcast_metatable = data.table::dcast.data.table(metatable, formula = variable~uniq_ID, value.var = 'value')
  testmatrix = dt_to_matrix(x = dcast_metatable)

  # correlation
  cormatrix = cor_flex(x = testmatrix, method = cor)
  cordist = stats::as.dist(1 - cormatrix, diag = T, upper = T)
  corclus = stats::hclust(d = cordist, method = distance)

  cordend = stats::as.dendrogram(object = corclus)

  # plot dendrogram
  pl = ggdendro::ggdendrogram(cordend, rotate = rotate, ...)

  # add horizontal or vertical lines
  if(!is.null(h)) {
    pl = pl + ggplot2::geom_hline(yintercept = h, col = h_color)
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



#' @title decide_cluster_order
#' @name decide_cluster_order
#' @description creates order for clusters
#' @param gobject giotto object
#' @param spat_unit spatial unit (e.g. "cell")
#' @param feat_type feature type (e.g. "rna", "dna", "protein")
#' @param expression_values expression values to use (e.g. "normalized", "scaled", "custom")
#' @param feats features to use (e.g. genes)
#' @param genes deprecated, use feats
#' @param cluster_column name of column to use for clusters (e.g. "leiden_clus")
#' @param cluster_order method to determine cluster order (e.g. "size", "correlation", "custom")
#' @param cluster_custom_order custom order for clusters
#' @param cor_method method for correlation, default to 'pearson'
#' @param hclust_method method for hierarchical clustering, default to 'ward.D'
#' @return custom
#' @details Calculates order for clusters.
#' @keywords internal
decide_cluster_order = function(gobject,
                                spat_unit = NULL,
                                feat_type = NULL,
                                expression_values = c('normalized', 'scaled', 'custom'),
                                feats,
                                genes = NULL,
                                cluster_column = NULL,
                                cluster_order = c('size', 'correlation', 'custom'),
                                cluster_custom_order = NULL,
                                cor_method = 'pearson',
                                hclust_method = 'ward.D') {

  ## deprecated arguments
  if(!is.null(genes)) {
    feats = genes
  }

  # Set feat_type and spat_unit
  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = spat_unit)
  feat_type = set_default_feat_type(gobject = gobject,
                                    spat_unit = spat_unit,
                                    feat_type = feat_type)

  # epxression data
  values = match.arg(expression_values, unique(c('normalized', 'scaled', 'custom', expression_values)))
  expr_values = get_expression_values(gobject = gobject,
                                      spat_unit = spat_unit,
                                      feat_type = feat_type,
                                      values = values,
                                      output = 'matrix')

  # subset expression data
  detected_feats = feats[feats %in% rownames(expr_values)]
  subset_values = expr_values[rownames(expr_values) %in% detected_feats, ]

  # metadata
  cell_metadata = pDataDT(gobject,
                          spat_unit = spat_unit,
                          feat_type = feat_type)

  ## check parameters
  if(is.null(cluster_column)) stop('\n cluster column must be selected \n')
  if(!cluster_column %in% colnames(cell_metadata)) stop('\n cluster column is not found \n')

  ## cluster order ##
  cluster_order = match.arg(cluster_order, c('size', 'correlation', 'custom'))


  if(cluster_order == 'size') {
    ## sorts clusters from big to small (# of cells per cluster)
    clus_sort_sizes = sort(table(cell_metadata[[cluster_column]]))
    clus_sort_names = names(clus_sort_sizes)


  } else if(cluster_order == 'correlation') {
    ## sorts clusters based on their correlation
    subset_matrix = create_cluster_matrix(gobject = gobject,
                                          spat_unit = spat_unit,
                                          feat_type = feat_type,
                                          cluster_column = cluster_column,
                                          feat_subset = detected_feats,
                                          expression_values = values)

    cormatrix = cor_flex(x = subset_matrix, method = cor_method)
    cordist = stats::as.dist(1 - cormatrix, diag = T, upper = T)
    corclus = stats::hclust(d = cordist, method = hclust_method)
    clus_names = rownames(cormatrix)
    names(clus_names) = 1:length(clus_names)
    clus_sort_names = clus_names[corclus$order]
    clus_sort_names = gsub(pattern = 'cluster_', replacement = '', x = clus_sort_names)


  } else if(cluster_order == 'custom') {
    ## sorts based on a given custom order
    if(is.null(cluster_custom_order)) {
      stop('\n custom order parameter is not given \n')
    }
    clus_sort_names = cluster_custom_order
  }
  return(clus_sort_names)
}


#' @title createHeatmap_DT
#' @name createHeatmap_DT
#' @description creates order for clusters
#' @param gobject giotto object
#' @param spat_unit spatial unit (e.g. "cell")
#' @param feat_type feature type (e.g. "rna", "dna", "protein")
#' @param expression_values expression values to use (e.g. "normalized", "scaled", "custom")
#' @param feats features to use
#' @param genes deprecated, use feats
#' @param cluster_column name of column to use for clusters (e.g. "leiden_clus")
#' @param cluster_order method to determine cluster order (e.g. "size", "correlation", "custom")
#' @param cluster_custom_order custom order for clusters
#' @param cluster_cor_method method for cluster correlation, default to "pearson"
#' @param cluster_hclust_method method for hierarchical clustering of clusters, default to "ward.D"
#' @param feat_order method to determine features order (e.g. "correlation", "custom")
#' @param gene_order deprecated, use feat_order in the future
#' @param feat_custom_order custom order for features
#' @param gene_custom_order deprecated, use feat_custom_order in the future
#' @param feat_cor_method method for features correlation, default to "pearson"
#' @param gene_cor_method deprecated, use feat_cor_method in the future
#' @param feat_hclust_method method for hierarchical clustering of features, default to "complete"
#' @param gene_hclust_method deprecated, use feat_hclust_method in the future
#' @return list
#' @details Creates input data.tables for plotHeatmap function.
#' @keywords internal
createHeatmap_DT = function(gobject,
                            spat_unit = NULL,
                            feat_type = NULL,
                            expression_values = c('normalized', 'scaled', 'custom'),
                            feats,
                            genes = NULL,
                            cluster_column = NULL,
                            cluster_order = c('size', 'correlation', 'custom'),
                            cluster_custom_order = NULL,
                            cluster_cor_method = 'pearson',
                            cluster_hclust_method = 'ward.D',
                            feat_order = c('correlation', 'custom'),
                            gene_order = NULL,
                            feat_custom_order = NULL,
                            gene_custom_order = NULL,
                            feat_cor_method = 'pearson',
                            gene_cor_method = NULL,
                            feat_hclust_method = 'complete',
                            gene_hclust_method = NULL) {


  ## deprecated arguments
  if(!is.null(genes)) {
    feats = genes
    warning('genes is deprecated, use feats in the future \n')
  }
  if(!is.null(gene_order)) {
    feat_order = gene_order
    warning('gene_order is deprecated, use feat_order in the future \n')
  }
  if(!is.null(gene_custom_order)) {
    feat_custom_order = gene_custom_order
    warning('gene_custom_order is deprecated, use feat_custom_order in the future \n')
  }
  if(!is.null(gene_cor_method)) {
    feat_cor_method = gene_cor_method
    warning('gene_cor_method is deprecated, use feat_cor_method in the future \n')
  }
  if(!is.null(gene_hclust_method)) {
    feat_hclust_method = gene_hclust_method
    warning('gene_hclust_method is deprecated, use feat_hclust_method in the future \n')
  }


  # Set feat_type and spat_unit
  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = spat_unit)
  feat_type = set_default_feat_type(gobject = gobject,
                                    spat_unit = spat_unit,
                                    feat_type = feat_type)


  # epxression data
  values = match.arg(expression_values, unique(c('normalized', 'scaled', 'custom', expression_values)))
  expr_values = get_expression_values(gobject = gobject,
                                      spat_unit = spat_unit,
                                      feat_type = feat_type,
                                      values = values,
                                      output = 'matrix')

  # subset expression data
  detected_feats = feats[feats %in% rownames(expr_values)]
  subset_values = expr_values[rownames(expr_values) %in% detected_feats, ]

  # metadata
  cell_metadata = pDataDT(gobject,
                          spat_unit = spat_unit,
                          feat_type = feat_type)

  # feat order
  feat_order = match.arg(feat_order, c('correlation', 'custom'))


  ### cluster order ###
  clus_sort_names = decide_cluster_order(gobject = gobject,
                                         spat_unit = spat_unit,
                                         feat_type = feat_type,
                                         expression_values = expression_values,
                                         feats = feats,
                                         cluster_column = cluster_column,
                                         cluster_order = cluster_order,
                                         cluster_custom_order = cluster_custom_order,
                                         cor_method = cluster_cor_method,
                                         hclust_method = cluster_hclust_method)

  ## data.table ##
  subset_values_DT <- data.table::as.data.table(reshape2::melt(as.matrix(subset_values),
                                                               varnames = c('feats', 'cells'),
                                                               value.name = 'expression'))
  subset_values_DT <- merge(subset_values_DT,
                            by.x = 'cells',
                            cell_metadata[, c('cell_ID', cluster_column), with = F],
                            by.y = 'cell_ID')
  subset_values_DT[[cluster_column]] <- factor(subset_values_DT[[cluster_column]], levels = clus_sort_names)

  # data.table variables
  z_scores = scale_scores = V1 = cells = NULL

  subset_values_DT[, feats := factor(feats, unique(detected_feats))]
  subset_values_DT[, z_scores := scale(expression), by = feats]
  subset_values_DT[, scale_scores := scales::rescale(x = expression, to = c(0,1)), by = feats]


  ## order cells by mean expression ##
  cell_order_DT = subset_values_DT[, mean(expression), by = c('cells', cluster_column)]
  cell_order_DT = cell_order_DT[order(get(cluster_column), V1)]
  subset_values_DT[, cells := factor(cells, cell_order_DT$cells)]

  ## get x-coordines for vertical lines in heatmap
  x_lines = cumsum(as.vector(table(cell_order_DT[[cluster_column]])))

  ## order feats ##
  if(feat_order == 'correlation') {
    featsum_per_clus = subset_values_DT[, sum(expression), by = c('feats', cluster_column)]

    my_formula = paste0('feats~',cluster_column)
    test_mat = data.table::dcast.data.table(data = featsum_per_clus, formula = my_formula, value.var = 'V1')
    test_matrix = as.matrix(test_mat[,-1]); rownames(test_matrix) = test_mat$feats

    feat_dist = stats::as.dist(1 - cor_flex(t_flex(test_matrix), method = feat_cor_method))
    feat_clus = stats::hclust(feat_dist, method = feat_hclust_method)

    feat_labels = rownames(test_matrix)
    feat_index = 1:length(feat_labels)
    names(feat_index) = feat_labels

    final_feat_order = names(feat_index[match(feat_clus$order, feat_index)])
    subset_values_DT[, feats := factor(feats, final_feat_order)]

  } else if(feat_order == 'custom') {

    if(is.null(feat_custom_order)) {
      stop('\n with custom feat order the feat_custom_order parameter needs to be provided \n')
    }
    subset_values_DT[, feats := factor(feats, feat_custom_order)]

  }

  cell_order_DT[['cells']] = factor(cell_order_DT[['cells']], levels = as.character(cell_order_DT[['cells']]))

  return(list(DT = subset_values_DT, x_lines = x_lines, cell_DT = cell_order_DT))
}



#' @title plotHeatmap
#' @name plotHeatmap
#' @description Creates heatmap for genes and clusters.
#' @param gobject giotto object
#' @param spat_unit spatial unit (e.g. "cell")
#' @param feat_type feature type (e.g. "rna", "dna", "protein")
#' @param expression_values expression values to use (e.g. "normalized", "scaled", "custom")
#' @param feats features to use
#' @param genes deprecated, use feats
#' @param cluster_column name of column to use for clusters (e.g. "leiden_clus")
#' @param cluster_order method to determine cluster order (e.g. "size", "correlation", "custom")
#' @param cluster_custom_order custom order for clusters
#' @param cluster_color_code color code for clusters
#' @param cluster_cor_method method for cluster correlation, default to "pearson"
#' @param cluster_hclust_method method for hierarchical clustering of clusters, default to "ward.D"
#' @param feat_order method to determine features order (e.g. "correlation", "custom")
#' @param feat_custom_order custom order for features
#' @param feat_cor_method method for features correlation, default to "pearson"
#' @param feat_hclust_method method for hierarchical clustering of features, default to "complete"
#' @param gene_order deprecated, use feat_order
#' @param gene_custom_order deprecated, use feat_custom_order
#' @param gene_cor_method deprecated, use feat_cor_method
#' @param gene_hclust_method deprecated, use feat_hclust_method
#' @param show_values which values to show on heatmap (e.g. "rescaled", "z-scaled", "original")
#' @param size_vertical_lines sizes for vertical lines
#' @param gradient_colors colors for heatmap gradient
#' @param feat_label_selection subset of features to show on y-axis
#' @param gene_label_selection deprecated, use feat_label_selection
#' @param axis_text_y_size size for y-axis text
#' @param legend_nrows number of rows for the cluster legend
#' @param show_plot show plot. TRUE or FALSE
#' @param return_plot return ggplot object. TRUE or FALSE
#' @param save_plot directly save the plot. TRUE or FALSE
#' @param save_param list of saving parameters, see \code{\link{showSaveParameters}}
#' @param default_save_name default save name
#' @return ggplot
#' @details If you want to display many genes there are 2 ways to proceed:
#' \itemize{
#'   \item{1. set axis_text_y_size to a really small value and show all features}
#'   \item{2. provide a subset of features to display to feat_label_selection}
#' }
#' @export
plotHeatmap <- function(gobject,
                        spat_unit = NULL,
                        feat_type = NULL,
                        expression_values = c('normalized', 'scaled', 'custom'),
                        feats,
                        genes = NULL,
                        cluster_column = NULL,
                        cluster_order = c('size', 'correlation', 'custom'),
                        cluster_custom_order = NULL,
                        cluster_color_code = NULL,
                        cluster_cor_method = 'pearson',
                        cluster_hclust_method = 'ward.D',
                        feat_order = c('correlation', 'custom'),
                        gene_order = NULL,
                        feat_custom_order = NULL,
                        gene_custom_order = NULL,
                        feat_cor_method = 'pearson',
                        gene_cor_method = NULL,
                        feat_hclust_method = 'complete',
                        gene_hclust_method = NULL,
                        show_values = c('rescaled', 'z-scaled', 'original'),
                        size_vertical_lines = 1.1,
                        gradient_colors = c('blue', 'yellow', 'red'),
                        feat_label_selection = NULL,
                        gene_label_selection = NULL,
                        axis_text_y_size = NULL,
                        legend_nrows = 1,
                        show_plot = NA,
                        return_plot = NA,
                        save_plot = NA,
                        save_param =  list(),
                        default_save_name = 'plotHeatmap') {

  show_values = match.arg(show_values, choices = c('rescaled', 'z-scaled', 'original'))


  heatmap_data = createHeatmap_DT(gobject = gobject,
                                  spat_unit = spat_unit,
                                  feat_type = feat_type,
                                  expression_values = expression_values,
                                  feats = feats,
                                  genes = genes,
                                  cluster_column = cluster_column,
                                  cluster_order = cluster_order,
                                  cluster_custom_order = cluster_custom_order,
                                  cluster_cor_method = cluster_cor_method,
                                  cluster_hclust_method = cluster_hclust_method,
                                  feat_order = feat_order,
                                  gene_order = gene_order,
                                  feat_custom_order = feat_custom_order,
                                  gene_custom_order = gene_custom_order,
                                  feat_cor_method = feat_cor_method,
                                  gene_cor_method = gene_cor_method,
                                  feat_hclust_method = feat_hclust_method,
                                  gene_hclust_method = gene_hclust_method)

  cell_order_DT = heatmap_data[['cell_DT']]
  subset_values_DT = heatmap_data[['DT']]
  x_lines = heatmap_data[['x_lines']]

  ## assign colors to each cluster
  if(is.null(cluster_color_code)) {
    clus_values = unique(cell_order_DT[[cluster_column]])
    clus_colors = getDistinctColors(n = length(clus_values))
    names(clus_colors) = clus_values
  } else {
    clus_colors = cluster_color_code
  }

  ## bar on top ##
  clus_pl <- ggplot2::ggplot()
  clus_pl <- clus_pl + ggplot2::geom_raster(data = cell_order_DT, ggplot2::aes_string(x = 'cells', y = '1', fill = cluster_column))
  clus_pl <- clus_pl + ggplot2::geom_vline(xintercept = x_lines, color = 'white', size = size_vertical_lines)
  clus_pl <- clus_pl + ggplot2::scale_fill_manual(values = clus_colors, guide = ggplot2::guide_legend(title = '', nrow = legend_nrows))
  clus_pl <- clus_pl + ggplot2::theme(axis.text = ggplot2::element_blank(),
                                      axis.ticks = ggplot2::element_blank(),
                                      axis.line = ggplot2::element_blank(),
                                      legend.position = 'top',
                                      plot.margin = ggplot2::margin(0, 0, 0, 5.5, "pt"))
  clus_pl <- clus_pl + ggplot2::labs(x = '', y = 'clusters')

  ## rescale values ##
  if(show_values == 'original') {
    value_column = 'expression'
    midpoint = max(subset_values_DT$expression)/2
  } else if(show_values == 'z-scaled') {
    value_column = 'z_scores'
    midpoint = max(subset_values_DT$z_scores)/2
  } else {
    value_column = 'scale_scores'
    midpoint = 0.5
  }

  # create empty plot
  empty = ggplot2::ggplot()
  empty = empty + ggplot2::theme_classic()
  empty = empty  + ggplot2::theme(axis.text = ggplot2::element_blank(),
                                  axis.ticks = ggplot2::element_blank(),
                                  axis.line = ggplot2::element_blank(),
                                  plot.margin = ggplot2::margin(0, 0, 0, 0, "cm"))

  ## parse gradient colors
  if(length(gradient_colors) > 3) cat('\n only first 3 colors will be used for gradient \n')
  low_color = gradient_colors[[1]]
  mid_color = gradient_colors[[2]]
  high_color = gradient_colors[[3]]

  ### heatmap ###
  hmap <- ggplot2::ggplot()
  hmap <- hmap + ggplot2::geom_tile(data = subset_values_DT, aes_string(x = 'cells', y = 'feats', fill = value_column))
  hmap <- hmap + ggplot2::geom_vline(xintercept = x_lines, color = 'white', size = size_vertical_lines)
  hmap <- hmap + ggplot2::scale_fill_gradient2(low = low_color, mid = mid_color, high = high_color,
                                               midpoint = midpoint,
                                               guide = ggplot2::guide_colorbar(title = ''))

  if(is.null(feat_label_selection)) {

    hmap <- hmap + ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                                  axis.ticks.x = ggplot2::element_blank(),
                                  axis.text.y = ggplot2::element_text(size = axis_text_y_size),
                                  plot.margin = ggplot2::margin(0, 0, 5.5, 5.5, "pt"))

    ## align and combine
    aligned <- cowplot::align_plots(clus_pl, empty, hmap + ggplot2::theme(legend.position = "none"),  align = "v", axis = "l")
    aligned <- append(aligned, list(cowplot::get_legend(hmap)))
    combplot = cowplot::plot_grid(plotlist = aligned,
                                  ncol = 2, rel_widths = c(1, 0.2),
                                  nrow = 2, rel_heights = c(0.2, 1))



    # print, return and save parameters
    show_plot = ifelse(is.na(show_plot), readGiottoInstructions(gobject, param = 'show_plot'), show_plot)
    save_plot = ifelse(is.na(save_plot), readGiottoInstructions(gobject, param = 'save_plot'), save_plot)
    return_plot = ifelse(is.na(return_plot), readGiottoInstructions(gobject, param = 'return_plot'), return_plot)

    ## print plot
    if(show_plot == TRUE) {
      print(combplot)
    }

    ## save plot
    if(save_plot == TRUE) {
      do.call('all_plots_save_function', c(list(gobject = gobject, plot_object = combplot, default_save_name = default_save_name), save_param))
    }

    ## return plot
    if(return_plot == TRUE) {
      return(combplot)
    }


  } else {

    # set defaults
    if(is.null(axis_text_y_size)) axis_text_y_size = 0.8*11/ggplot2::.pt

    # finish heatmap
    hmap <- hmap + ggplot2::theme(axis.text = ggplot2::element_blank(),
                                  axis.ticks = ggplot2::element_blank(),
                                  plot.margin = ggplot2::margin(0, 0, 5.5, 5.5, "pt"))
    ### axis ###
    featDT = subset_values_DT[,c('feats'), with = F]
    featDT = unique(setorder(featDT, feats))

    # data.table variables
    featOrder = subset_feats = NULL

    featDT[, featOrder := 1:.N]
    featDT[, subset_feats := ifelse(feats %in% feat_label_selection, as.character(feats), '')]

    axis <- ggplot2::ggplot(data = featDT, aes(x = 0, y = featOrder, label = subset_feats))
    axis <- axis + ggrepel::geom_text_repel(min.segment.length = grid::unit(0, "pt"),
                                            color = "grey30",  ## ggplot2 theme_grey() axis text
                                            size = axis_text_y_size  ## default ggplot2 theme_grey() axis text
    )
    axis <- axis + ggplot2::scale_x_continuous(limits = c(0, 1), expand = c(0, 0),
                                      breaks = NULL, labels = NULL, name = NULL)
    axis <- axis + ggplot2::scale_y_continuous(limits = c(0, nrow(featDT)), expand = c(0, 0),
                                      breaks = NULL, labels = NULL, name = NULL)
    axis <- axis + ggplot2::theme(panel.background = ggplot2::element_blank(),
                         plot.margin = ggplot2::margin(0, 0, 0, 0, "cm"))


    #return(list(hmap, clus_pl, axis, empty))

    ## align and combine
    aligned <- cowplot::align_plots(clus_pl, empty, empty, hmap + theme(legend.position = "none"), axis, align = "h", axis = "b")
    aligned <- append(aligned, list(cowplot::get_legend(hmap)))
    combplot = cowplot::plot_grid(plotlist = aligned,
                                  ncol = 3, rel_widths = c(1, 0.2, 0.1),
                                  nrow = 2, rel_heights = c(0.2, 1))


    # print, return and save parameters
    show_plot = ifelse(is.na(show_plot), readGiottoInstructions(gobject, param = 'show_plot'), show_plot)
    save_plot = ifelse(is.na(save_plot), readGiottoInstructions(gobject, param = 'save_plot'), save_plot)
    return_plot = ifelse(is.na(return_plot), readGiottoInstructions(gobject, param = 'return_plot'), return_plot)


    ## print plot
    if(show_plot == TRUE) {
      print(combplot)
    }

    ## save plot
    if(save_plot == TRUE) {
      do.call('all_plots_save_function', c(list(gobject = gobject, plot_object = combplot, default_save_name = default_save_name), save_param))
    }

    ## return plot
    if(return_plot == TRUE) {
      return(combplot)
    }

  }

}






#' @title plotMetaDataHeatmap
#' @name plotMetaDataHeatmap
#' @description Creates heatmap for genes within aggregated clusters.
#' @param gobject giotto object
#' @param spat_unit spatial unit (e.g. "cell")
#' @param feat_type feature type (e.g. "rna", "dna", "protein")
#' @param expression_values expression values to use (e.g. "normalized", "scaled", "custom")
#' @param metadata_cols annotation columns found in pDataDT(gobject)
#' @param selected_feats subset of features to use
#' @param selected_genes deprecated. See \code{selected_feats} param
#' @param first_meta_col if more than 1 metadata column, select the x-axis factor
#' @param second_meta_col if more than 1 metadata column, select the facetting factor
#' @param show_values which values to show on heatmap (e.g. "zscores", "original", "zscores_rescaled")
#' @param custom_cluster_order custom cluster order (default = NULL)
#' @param clus_cor_method correlation method for clusters, default to "pearson"
#' @param clus_cluster_method hierarchical cluster method for the clusters, default to "complete"
#' @param custom_feat_order custom feature order (default = NULL)
#' @param custom_gene_order deprecated. See \code{custom_feat_order} param
#' @param feat_cor_method correlation method for features, default to "pearson"
#' @param gene_cor_method deprecated. See \code{feat_cor_method} param
#' @param feat_cluster_method hierarchical cluster method for the features, default to "complete"
#' @param gene_cluster_method deprecated. See \code{feat_cluster_method} param
#' @param gradient_color vector with 3 colors for numeric data
#' @param gradient_midpoint midpoint for color gradient
#' @param gradient_limits vector with lower and upper limits
#' @param x_text_size size of x-axis text
#' @param x_text_angle angle of x-axis text
#' @param y_text_size size of y-axis text
#' @param strip_text_size size of strip text
#' @param show_plot show plot. TRUE or FALSE
#' @param return_plot return ggplot object. TRUE or FALSE
#' @param save_plot directly save the plot. TRUE or FALSE
#' @param save_param list of saving parameters, see \code{\link{showSaveParameters}}
#' @param default_save_name default save name
#' @return ggplot or data.table
#' @details Creates heatmap for the average expression of selected genes in the different annotation/cluster groups.
#' Calculation of cluster or gene order is done on the provided expression values, but visualization
#' is by default on the z-scores. Other options are the original values or z-scores rescaled per gene (-1 to 1).
#' @seealso \code{\link{plotMetaDataCellsHeatmap}} for numeric cell annotation instead of gene expression.
#' @export
plotMetaDataHeatmap = function(gobject,
                               spat_unit = NULL,
                               feat_type = NULL,
                               expression_values =  c("normalized", "scaled", "custom"),
                               metadata_cols = NULL,
                               selected_feats = NULL,
                               selected_genes = NULL,
                               first_meta_col = NULL,
                               second_meta_col = NULL,
                               show_values = c('zscores', 'original', 'zscores_rescaled'),
                               custom_cluster_order = NULL,
                               clus_cor_method = 'pearson',
                               clus_cluster_method = 'complete',
                               custom_feat_order = NULL,
                               custom_gene_order = NULL,
                               feat_cor_method = 'pearson',
                               gene_cor_method = NULL,
                               feat_cluster_method = 'complete',
                               gene_cluster_method = NULL,
                               gradient_color = c('blue', 'white', 'red'),
                               gradient_midpoint = 0,
                               gradient_limits = NULL,
                               x_text_size = 10,
                               x_text_angle = 45,
                               y_text_size = 10,
                               strip_text_size = 8,
                               show_plot = NA,
                               return_plot = NA,
                               save_plot = NA,
                               save_param =  list(),
                               default_save_name = 'plotMetaDataHeatmap') {


  ## deprecated arguments
  if(!is.null(selected_genes)) {
    selected_feats = selected_genes
    warning('selected_genes is deprecated, use selected_feats in the future \n')
  }
  if(!is.null(custom_gene_order)) {
    custom_feat_order = custom_gene_order
    warning('custom_gene_order is deprecated, use custom_feat_order in the future \n')
  }
  if(!is.null(gene_cor_method)) {
    feat_cor_method = gene_cor_method
    warning('gene_cor_method is deprecated, use feat_cor_method in the future \n')
  }
  if(!is.null(gene_cluster_method)) {
    feat_cluster_method = gene_cluster_method
    warning('gene_cluster_method is deprecated, use feat_cluster_method in the future \n')
  }

  # Set feat_type and spat_unit
  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = spat_unit)
  feat_type = set_default_feat_type(gobject = gobject,
                                    spat_unit = spat_unit,
                                    feat_type = feat_type)

  metaDT = calculateMetaTable(gobject = gobject,
                              spat_unit = spat_unit,
                              feat_type = feat_type,
                              expression_values = expression_values,
                              metadata_cols = metadata_cols,
                              selected_feats = selected_feats)

  # data.table variables
  zscores = value = zscores_rescaled_per_feat = NULL

  metaDT[, zscores := scale(value), by = c('variable')]
  metaDT[, zscores_rescaled_per_feat := scales::rescale(zscores, to = c(-1,1)), by = c('variable')]

  show_values = match.arg(show_values, choices = c('zscores', 'original', 'zscores_rescaled'))
  if(show_values == 'zscores') {
    show_values = 'zscores'
  } else if(show_values == 'original') {
    show_values = 'value'
  } else {
    show_values = 'zscores_rescaled_per_feat'
  }

  ## visualization
  if(length(metadata_cols) > 2) {
    cat('\n visualization is only possible for 1 or 2 metadata annotations, data.table is returned \n')
    return(metaDT)
  }


  ## order of feats and clusters

  main_factor = ifelse(length(metadata_cols) == 1, metadata_cols, first_meta_col)
  testmain = metaDT[, mean(value), by = c('variable', main_factor)]
  dfunction <- function(d, col_name1, col_name2, value.var) {
    data.table::dcast.data.table(d, paste(col_name1, "~", col_name2), value.var = value.var)
  }
  testmain_matrix = dfunction(d = testmain, col_name1 = 'variable', col_name2 = main_factor, value.var = 'V1')
  testmain_mat = as.matrix(testmain_matrix[,-1]); rownames(testmain_mat) = testmain_matrix$variable

  # for clusters
  if(is.null(custom_cluster_order)) {
    cormatrix = cor_flex(x = testmain_mat, method = clus_cor_method)
    cordist = stats::as.dist(1 - cormatrix, diag = T, upper = T)
    corclus = stats::hclust(d = cordist, method = clus_cluster_method)
    clus_names = rownames(cormatrix)
    names(clus_names) = 1:length(clus_names)
    clus_sort_names = clus_names[corclus$order]
  } else {
    clus_sort_names = unique(as.character(custom_cluster_order))
    if(all(colnames(testmain_mat) %in% clus_sort_names) == FALSE) {
      stop('\n custom cluster order is given, but not all clusters are represented \n')
    }
  }


  # for feats
  if(is.null(custom_feat_order)) {
    feat_cormatrix = cor_flex(x = t(testmain_mat), method = feat_cor_method)
    feat_cordist = stats::as.dist(1 - feat_cormatrix, diag = T, upper = T)
    feat_corclus = stats::hclust(d = feat_cordist, method = feat_cluster_method)
    feat_names = rownames(feat_cormatrix)
    names(feat_names) = 1:length(feat_names)
    feat_sort_names = feat_names[feat_corclus$order]
  }  else {
    feat_sort_names = unique(as.character(custom_feat_order))
    if(all(rownames(testmain_mat) %in% feat_sort_names) == FALSE) {
      stop('\n custom feat order is given, but not all feats are represented \n')
    }
  }

  if(length(metadata_cols) == 1) {

    # data.table variables
    factor_column = variable = NULL

    metaDT[, factor_column := factor(get(metadata_cols), levels = clus_sort_names)]
    metaDT[, variable := factor(get('variable'), levels = feat_sort_names)]

    # set gradient information
    if(!is.null(gradient_limits) & is.vector(gradient_limits) & length(gradient_limits) == 2) {
      lower_lim = gradient_limits[[1]]
      upper_lim = gradient_limits[[2]]

      numeric_data = metaDT[[show_values]]
      limit_numeric_data = ifelse(numeric_data > upper_lim, upper_lim,
                                  ifelse(numeric_data < lower_lim, lower_lim, numeric_data))
      metaDT[[show_values]] = limit_numeric_data
    }


    pl <- ggplot2::ggplot()
    pl <- pl + ggplot2::geom_tile(data = metaDT, ggplot2::aes_string(x = 'factor_column', y = 'variable', fill = show_values), color = 'black')
    pl <- pl + ggplot2::scale_fill_gradient2(low = gradient_color[[1]],
                                    mid = gradient_color[[2]],
                                    high = gradient_color[[3]],
                                    midpoint = gradient_midpoint)
    pl <- pl + ggplot2::theme_classic()
    pl <- pl + ggplot2::theme(axis.text.x = ggplot2::element_text(size = x_text_size, angle = x_text_angle, hjust = 1, vjust = 1),
                     axis.text.y = ggplot2::element_text(size = y_text_size),
                     legend.title = ggplot2::element_blank())
    pl <- pl + ggplot2::labs(x = metadata_cols, y = 'feats')


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


  else {

    if(is.null(first_meta_col) | is.null(second_meta_col)) {
      cat('\n both first_meta_col and second_meta_col need to be defined, return data.table \n')
      return(metaDT)
    } else {

      # data.table variables
      factor_1_column = factor_2_column = variable = NULL

      metaDT[, factor_1_column := factor(get(first_meta_col), clus_sort_names)]
      metaDT[, factor_2_column := as.factor(get(second_meta_col))]
      metaDT[, variable := factor(get('variable'), levels = feat_sort_names)]

      # set gradient information
      if(!is.null(gradient_limits) & is.vector(gradient_limits) & length(gradient_limits) == 2) {
        lower_lim = gradient_limits[[1]]
        upper_lim = gradient_limits[[2]]

        numeric_data = metaDT[[show_values]]
        limit_numeric_data = ifelse(numeric_data > upper_lim, upper_lim,
                                    ifelse(numeric_data < lower_lim, lower_lim, numeric_data))
        metaDT[[show_values]] = limit_numeric_data
      }

      pl <- ggplot2::ggplot()
      pl <- pl + ggplot2::geom_tile(data = metaDT, ggplot2::aes_string(x = 'factor_1_column', y = 'variable', fill = show_values), color = 'black')
      pl <- pl + ggplot2::scale_fill_gradient2(low = gradient_color[[1]],
                                      mid = gradient_color[[2]],
                                      high = gradient_color[[3]],
                                      midpoint = gradient_midpoint)
      pl <- pl + ggplot2::facet_grid(stats::reformulate('factor_2_column'))
      pl <- pl + ggplot2::theme_classic()
      pl <- pl + ggplot2::theme(axis.text.x = ggplot2::element_text(size = x_text_size, angle = x_text_angle, hjust = 1, vjust = 1),
                       axis.text.y = ggplot2::element_text(size = y_text_size),
                       strip.text = ggplot2::element_text(size = strip_text_size),
                       legend.title= ggplot2::element_blank())
      pl <- pl + ggplot2::labs(x = first_meta_col, y = 'feats', title = second_meta_col)


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

  }

}


#' @title plotMetaDataCellsHeatmap
#' @name plotMetaDataCellsHeatmap
#' @description Creates heatmap for numeric cell metadata within aggregated clusters.
#' @param gobject giotto object
#' @param spat_unit spatial unit (e.g. "cell")
#' @param feat_type feature type (e.g. "rna", "dna", "protein")
#' @param metadata_cols annotation columns found in pDataDT(gobject)
#' @param spat_enr_names spatial enrichment results to include
#' @param value_cols value columns to use
#' @param first_meta_col if more than 1 metadata column, select the x-axis factor
#' @param second_meta_col if more than 1 metadata column, select the facetting factor
#' @param show_values which values to show on heatmap (e.g. "zscores", "original", "zscores_rescaled")
#' @param custom_cluster_order custom cluster order (default = NULL)
#' @param clus_cor_method correlation method for clusters, default to "pearson"
#' @param clus_cluster_method hierarchical cluster method for the clusters, default to "complete"
#' @param custom_values_order custom values order (default = NULL)
#' @param values_cor_method correlation method for values, default to "pearson"
#' @param values_cluster_method hierarchical cluster method for the values, default to "complete"
#' @param midpoint midpoint of show_values
#' @param x_text_size size of x-axis text
#' @param x_text_angle angle of x-axis text
#' @param y_text_size size of y-axis text
#' @param strip_text_size size of strip text
#' @param show_plot show plot. TRUE or FALSE
#' @param return_plot return ggplot object. TRUE or FALSE
#' @param save_plot directly save the plot. TRUE or FALSE
#' @param save_param list of saving parameters, see \code{\link{showSaveParameters}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return ggplot or data.table
#' @details Creates heatmap for the average values of selected value columns in the different annotation groups.
#' @seealso \code{\link{plotMetaDataHeatmap}} for gene expression instead of numeric cell annotation data.
#' @export
plotMetaDataCellsHeatmap = function(gobject,
                                    spat_unit = NULL,
                                    feat_type = NULL,
                                    metadata_cols = NULL,
                                    spat_enr_names = NULL,
                                    value_cols = NULL,
                                    first_meta_col = NULL,
                                    second_meta_col = NULL,
                                    show_values = c('zscores', 'original', 'zscores_rescaled'),
                                    custom_cluster_order = NULL,
                                    clus_cor_method = 'pearson',
                                    clus_cluster_method = 'complete',
                                    custom_values_order = NULL,
                                    values_cor_method = 'pearson',
                                    values_cluster_method = 'complete',
                                    midpoint = 0,
                                    x_text_size = 8,
                                    x_text_angle = 45,
                                    y_text_size = 8,
                                    strip_text_size = 8,
                                    show_plot = NA,
                                    return_plot = NA,
                                    save_plot = NA,
                                    save_param =  list(),
                                    default_save_name = 'plotMetaDataCellsHeatmap') {


  metaDT = calculateMetaTableCells(gobject = gobject,
                                   spat_unit = spat_unit,
                                   feat_type = feat_type,
                                   value_cols = value_cols,
                                   metadata_cols = metadata_cols,
                                   spat_enr_names = spat_enr_names)

  # data.table variables
  zscores = zscores_rescaled_per_gene = value = NULL

  metaDT[, zscores := scale(value), by = c('variable')]
  metaDT[, zscores_rescaled_per_gene := scales::rescale(zscores, to = c(-1,1)), by = c('variable')]

  show_values = match.arg(show_values, choices = c('zscores', 'original', 'zscores_rescaled'))
  if(show_values == 'zscores') {
    show_values = 'zscores'
  } else if(show_values == 'original') {
    show_values = 'value'
  } else {
    show_values = 'zscores_rescaled_per_gene'
  }

  ## visualization
  if(length(metadata_cols) > 2) {
    cat('\n visualization is only possible for 1 or 2 metadata annotations, data.table is returned \n')
    return(metaDT)
  }


  ## order of genes and clusters

  main_factor = ifelse(length(metadata_cols) == 1, metadata_cols, first_meta_col)
  testmain = metaDT[, mean(value), by = c('variable', main_factor)]
  dfunction <- function(d, col_name1, col_name2, value.var) {
    data.table::dcast.data.table(d, paste(col_name1, "~", col_name2), value.var = value.var)
  }
  testmain_matrix = dfunction(d = testmain, col_name1 = 'variable', col_name2 = main_factor, value.var = 'V1')
  testmain_mat = as.matrix(testmain_matrix[,-1]); rownames(testmain_mat) = testmain_matrix$variable

  # for clusters
  if(is.null(custom_cluster_order)) {
    cormatrix = cor_flex(x = testmain_mat, method = clus_cor_method)
    cordist = stats::as.dist(1 - cormatrix, diag = T, upper = T)
    corclus = stats::hclust(d = cordist, method = clus_cluster_method)
    clus_names = rownames(cormatrix)
    names(clus_names) = 1:length(clus_names)
    clus_sort_names = clus_names[corclus$order]
  } else {
    clus_sort_names = unique(as.character(custom_cluster_order))
    if(all(colnames(testmain_mat) %in% clus_sort_names) == FALSE) {
      stop('\n custom cluster order is given, but not all clusters are represented \n')
    }
  }


  # for genes
  if(is.null(custom_values_order)) {
    values_cormatrix = cor_flex(x = t(testmain_mat), method = values_cor_method)
    values_cordist = stats::as.dist(1 - values_cormatrix, diag = T, upper = T)
    values_corclus = stats::hclust(d = values_cordist, method = values_cluster_method)
    values_names = rownames(values_cormatrix)
    names(values_names) = 1:length(values_names)
    values_sort_names = values_names[values_corclus$order]
  }  else {
    values_sort_names = unique(as.character(custom_values_order))
    if(all(rownames(testmain_mat) %in% values_sort_names) == FALSE) {
      stop('\n custom values order is given, but not all values are represented \n')
    }
  }

  if(length(metadata_cols) == 1) {

    # data.table variables
    factor_column = variable = NULL

    metaDT[, factor_column := factor(get(metadata_cols), levels = clus_sort_names)]
    metaDT[, variable := factor(get('variable'), levels = values_sort_names)]

    pl <- ggplot2::ggplot()
    pl <- pl + ggplot2::geom_tile(data = metaDT, ggplot2::aes_string(x = 'factor_column', y = 'variable', fill = show_values), color = 'black')
    pl <- pl + ggplot2::scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = midpoint)
    pl <- pl + ggplot2::theme_classic()
    pl <- pl + ggplot2::theme(axis.text.x = ggplot2::element_text(size = x_text_size, angle = x_text_angle, hjust = 1, vjust = 1),
                     axis.text.y = ggplot2::element_text(size = y_text_size),
                     legend.title = ggplot2::element_blank())
    pl <- pl + ggplot2::labs(x = metadata_cols, y = 'genes')


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


  else {

    if(is.null(first_meta_col) | is.null(second_meta_col)) {
      cat('\n both first_meta_col and second_meta_col need to be defined, return data.table \n')
      return(metaDT)
    } else {

      # data.table variables
      factor_1_column = factor_2_column = variable = NULL

      metaDT[, factor_1_column := factor(get(first_meta_col), clus_sort_names)]
      metaDT[, factor_2_column := as.factor(get(second_meta_col))]
      metaDT[, variable := factor(get('variable'), levels = values_sort_names)]

      pl <- ggplot2::ggplot()
      pl <- pl + ggplot2::geom_tile(data = metaDT, ggplot2::aes_string(x = 'factor_1_column', y = 'variable', fill = show_values), color = 'black')
      pl <- pl + ggplot2::scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = midpoint)
      pl <- pl + ggplot2::facet_grid(stats::reformulate('factor_2_column'))
      pl <- pl + ggplot2::theme_classic()
      pl <- pl + ggplot2::theme(axis.text.x = ggplot2::element_text(size = x_text_size, angle = x_text_angle, hjust = 1, vjust = 1),
                       axis.text.y = ggplot2::element_text(size = y_text_size),
                       strip.text = ggplot2::element_text(size = strip_text_size),
                       legend.title = ggplot2::element_blank())
      pl <- pl + ggplot2::labs(x = first_meta_col, y = 'genes', title = second_meta_col)


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

  }

}





#' @title violinPlot
#' @name violinPlot
#' @description Creates violinplot for selected clusters
#' @param gobject giotto object
#' @param spat_unit spatial unit (e.g. "cell")
#' @param feat_type feature type (e.g. "rna", "dna", "protein")
#' @param expression_values expression values to use (e.g. "normalized", "scaled", "custom")
#' @param feats features to plot
#' @param genes deprecated, use feats argument
#' @param cluster_column name of column to use for clusters (e.g. "leiden_clus")
#' @param color_violin color violin according to "genes" or "clusters"
#' @param cluster_custom_order custom order of clusters
#' @param cluster_color_code color code for clusters
#' @param strip_position position of gene labels (e.g. "top", "right", "left", "bottom")
#' @param strip_text size of strip text
#' @param axis_text_x_size size of x-axis text
#' @param axis_text_y_size size of y-axis text
#' @param show_plot show plot. TRUE or FALSE
#' @param return_plot return ggplot object. TRUE or FALSE
#' @param save_plot directly save the plot. TRUE or FALSE
#' @param save_param list of saving parameters, see \code{\link{showSaveParameters}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return ggplot
#' @export
violinPlot <- function(gobject,
                       spat_unit = NULL,
                       feat_type = NULL,
                       expression_values = c('normalized', 'scaled', 'custom'),
                       feats = NULL,
                       genes = NULL,
                       cluster_column,
                       cluster_custom_order = NULL,
                       color_violin = c('feats', 'cluster'),
                       cluster_color_code = NULL,
                       strip_position = c('top', 'right', 'left', 'bottom'),
                       strip_text = 7,
                       axis_text_x_size = 10,
                       axis_text_y_size = 6,
                       show_plot = NA,
                       return_plot = NA,
                       save_plot = NA,
                       save_param =  list(),
                       default_save_name = 'violinPlot') {


  ## deprecated arguments
  if(!is.null(genes)) {
    feats = genes
    warning('genes is deprecated, use feats in the future \n')
  }

  # Set feat_type and spat_unit
  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = spat_unit)
  feat_type = set_default_feat_type(gobject = gobject,
                                    spat_unit = spat_unit,
                                    feat_type = feat_type)

  ## strip position
  strip_position = match.arg(strip_position, c('top', 'right', 'left', 'bottom'))

  ## color of violin plots
  color_violin = match.arg(color_violin, c('feats', 'cluster'))

  ## expression data ##
  values = match.arg(expression_values, unique(c('normalized', 'scaled', 'custom', expression_values)))
  expr_data = get_expression_values(gobject = gobject,
                                    feat_type = feat_type,
                                    spat_unit = spat_unit,
                                    values = values,
                                    output = 'matrix')

  # only keep feats that are in the dataset
  selected_feats = feats[feats %in% rownames(expr_data)]
  if(length(selected_feats[duplicated(selected_feats)]) != 0) {
    cat('These feats have duplicates: \n',
        selected_feats[duplicated(selected_feats)])
    selected_feats = unique(selected_feats)
  }

  # stop and provide warning if no feats have been found
  if(length(selected_feats) == 0) {
    stop("No overlapping features have been found, check inputer for parameter 'feats'")
  }

  subset_data = as.matrix(expr_data[rownames(expr_data) %in% selected_feats, ])

  if(length(feats) == 1) {
    t_subset_data = subset_data
  } else {
    t_subset_data = t_flex(subset_data)
  }

  # metadata
  metadata = pDataDT(gobject,
                     feat_type = feat_type,
                     spat_unit = spat_unit)

  if(length(feats) == 1) {
    metadata_expr <- cbind(metadata,  t_subset_data)
    setnames(metadata_expr, 'V1', feats)
  } else {
    metadata_expr <- cbind(metadata, t_subset_data)
  }


  metadata_expr_m = data.table::melt.data.table(metadata_expr, measure.vars = unique(selected_feats), variable.name = 'feats')
  metadata_expr_m[, feats := factor(feats, selected_feats)]
  metadata_expr_m[[cluster_column]] = as.factor(metadata_expr_m[[cluster_column]])

  if(!is.null(cluster_custom_order)) {
    cluster_custom_order = unique(cluster_custom_order)
    metadata_expr_m[[cluster_column]] = factor(x = metadata_expr_m[[cluster_column]], levels = cluster_custom_order)
  }


  pl <- ggplot2::ggplot()
  pl <- pl + ggplot2::theme_classic()

  if(color_violin == 'feats') {
    pl <- pl + ggplot2::geom_violin(data = metadata_expr_m, aes_string(x = cluster_column, y = 'value', fill = 'feats'), width = 1, scale = 'width', show.legend = F)
  } else {
    pl <- pl + ggplot2::geom_violin(data = metadata_expr_m,
                                    aes_string(x = cluster_column, y = 'value',
                                               fill = cluster_column),
                                    width = 1, scale = 'width',
                                    show.legend = F)

    # provide own color scheme for clusters
    if(!is.null(cluster_color_code)) {
      pl <- pl + ggplot2::scale_fill_manual(values = cluster_color_code)
    }
  }


  pl <- pl + ggplot2::facet_wrap(.~feats, ncol = 1, strip.position = strip_position)
  pl <- pl + ggplot2::theme(strip.text = element_text(size = strip_text),
                            axis.text.x = element_text(size = axis_text_x_size, angle = 45, hjust = 1, vjust = 1),
                            axis.text.y = element_text(size = axis_text_y_size))
  pl <- pl + ggplot2::labs(x = '', y = 'normalized expression')



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

#' @title plotly_network
#' @name plotly_network
#' @description provide network segment to draw in 3D plot_ly()
#' @param network network object
#' @param x default to "sdimx_begin"
#' @param y default to "sdimy_begin"
#' @param z default to "sdimz_begin"
#' @param x_end default to "sdimx_end"
#' @param y_end default to "sdimy_end"
#' @param z_end default to "sdimz_end"
#'
#' @return edges in network as data.table()
#' @keywords internal
plotly_network <- function(network,
                           x = "sdimx_begin",
                           y = "sdimy_begin",
                           z = "sdimz_begin",
                           x_end = "sdimx_end",
                           y_end="sdimy_end",
                           z_end="sdimz_end"){

  edges = data.table::data.table(edge_id = 1:(3*dim(network)[1]),
                                 x = 0,
                                 y = 0,
                                 z = 0)

  edges[edges$edge_id%%3 == 1]$x = as.double(network[[x]])
  edges[edges$edge_id%%3 == 1]$y = as.double(network[[y]])
  edges[edges$edge_id%%3 == 1]$z = as.double(network[[z]])

  edges[edges$edge_id%%3 == 2]$x = as.double(network[[x_end]])
  edges[edges$edge_id%%3 == 2]$y = as.double(network[[y_end]])
  edges[edges$edge_id%%3 == 2]$z = as.double(network[[z_end]])

  edges[edges$edge_id%%3 == 0]$x = NA
  edges[edges$edge_id%%3 == 0]$y = NA
  edges[edges$edge_id%%3 == 0]$z = NA

  return(edges)
}


#' @title plotly_grid
#' @name plotly_grid
#' @description provide grid segment to draw in plot_ly()
#'
#' @param x_start default to "x_start"
#' @param y_start default to "y_start"
#' @param x_end default to "x_end"
#' @param y_end default to "y_end"
#' @param spatial_grid spatial_grid in giotto object
#'
#' @return edges in spatial grid as data.table()
#' @keywords internal
plotly_grid <- function(spatial_grid,
                        x_start = "x_start",
                        y_start = "y_start",
                        x_end = "x_end",
                        y_end = "y_end"){

  edge_num <- length(unique(spatial_grid[[x_start]])) + length(unique(spatial_grid[[y_start]])) + 2
  x_line <- unique(as.numeric(unlist(spatial_grid[,c(x_start,x_end)])))
  y_line <- unique(as.numeric(unlist(spatial_grid[,c(y_start,y_end)])))

  x_min <- min(spatial_grid[[x_start]])
  x_max <- max(spatial_grid[[x_end]])

  y_min <- min(spatial_grid[[y_start]])
  y_max <- max(spatial_grid[[y_end]])

  edges <- data.table::data.table(edge_id = 1:edge_num,x = 0,y = 0,x_end = 0,y_end = 0)

  edges[1:length(x_line),]$x <- x_line
  edges[1:length(x_line),]$x_end <- x_line
  edges[1:length(x_line),]$y <- y_min
  edges[1:length(x_line),]$y_end <- y_max

  edges[(length(x_line)+1):edge_num,]$x <- x_min
  edges[(length(x_line)+1):edge_num,]$x_end <- x_max
  edges[(length(x_line)+1):edge_num,]$y <- y_line
  edges[(length(x_line)+1):edge_num,]$y_end <- y_line

  return(edges)
}


#' @title plotly_axis_scale_3D
#' @name plotly_axis_scale_3D
#' @description adjust the axis scale in 3D plotly plot
#' @param cell_locations spatial_loc in giotto object
#' @param sdimx x axis of cell spatial location
#' @param sdimy y axis of cell spatial location
#' @param sdimz z axis of cell spatial location
#' @param mode axis adjustment mode
#' @param custom_ratio set the ratio artificially
#' @return edges in spatial grid as data.table()
#' @keywords internal
plotly_axis_scale_3D <- function(cell_locations,
                                 sdimx = NULL,
                                 sdimy = NULL,
                                 sdimz = NULL,
                                 mode = c("cube","real","custom"),
                                 custom_ratio = NULL){
  mode = match.arg(mode, c("cube","real","custom"))
  if(mode == "real"){
    x_ratio = max(cell_locations[[sdimx]]) - min(cell_locations[[sdimx]])
    y_ratio = max(cell_locations[[sdimy]]) - min(cell_locations[[sdimy]])
    z_ratio = max(cell_locations[[sdimz]]) - min(cell_locations[[sdimz]])

    min_size = min(x_ratio,y_ratio,z_ratio)
    x_ratio = round(x_ratio/min_size)
    y_ratio = round(y_ratio/min_size)
    z_ratio = round(z_ratio/min_size)
  }
  else if(mode == "cube"){
    x_ratio = 1
    y_ratio = 1
    z_ratio = 1
  }
  else{
    if(is.null(custom_ratio) | length(custom_ratio) < 3){
      stop("\n Please specify your costom axis scale, or choose axis_scale = \"real\"/\"cube\"\n")
    }
    else{
      x_ratio = custom_ratio[1]
      y_ratio = custom_ratio[2]
      z_ratio = custom_ratio[3]
    }
  }
  ratio <- list(x_ratio,y_ratio,z_ratio)
  return (ratio)
}


#' @title plotly_axis_scale_2D
#' @name plotly_axis_scale_2D
#' @description adjust the axis scale in 2D plotly plot
#' @param cell_locations spatial_loc in giotto object
#' @param sdimx x axis of cell spatial location
#' @param sdimy y axis of cell spatial location
#' @param mode axis adjustment mode
#' @param custom_ratio set the ratio artificially
#' @return edges in spatial grid as data.table()
#' @keywords internal
plotly_axis_scale_2D <- function(cell_locations,
                                 sdimx = NULL,
                                 sdimy = NULL,
                                 mode = c("cube","real","custom"),
                                 custom_ratio = NULL){

  mode = match.arg(mode, c("cube","real","custom"))
  if(mode == "real"){
    x_ratio = max(cell_locations[[sdimx]]) - min(cell_locations[[sdimx]])
    y_ratio = max(cell_locations[[sdimy]]) - min(cell_locations[[sdimy]])

    min_size = min(x_ratio,y_ratio)
    x_ratio = round(x_ratio/min_size)
    y_ratio = round(y_ratio/min_size)
  }
  else if(mode == "cube"){
    x_ratio = 1
    y_ratio = 1
  }
  else{
    if(is.null(custom_ratio) | length(custom_ratio) < 2){
      stop("\n Please specify your costom axis scale, or choose axis_scale = \"real\"/\"cube\"\n")
    }
    else{
      x_ratio = custom_ratio[1]
      y_ratio = custom_ratio[2]
    }
  }
  ratio <- list(x_ratio,y_ratio)
  return (ratio)
}











