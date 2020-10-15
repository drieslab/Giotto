
#' @title aes_string2
#' @name aes_string2
#' @param \dots aes_string parameters
#' @keywords internal
#' @description makes sure aes_string can also be used with names that start with numeric values
aes_string2 <- function(...){
  args <- lapply(list(...), function(x) sprintf("`%s`", x))
  do.call(ggplot2::aes_string, args)
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
                                ...) {

  if(is.null(plot_object)) {
    stop('\t there is no object to plot \t')
  }

  ## get save information and set defaults
  if(is.null(save_dir)) save_dir = readGiottoInstructions(gobject, param = 'save_dir')
  if(is.null(save_folder)) save_folder = NULL
  if(is.null(save_name)) save_name = default_save_name
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
  if(is.null(save_name)) save_name = default_save_name
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
                          ...)

  }

}




#' @title showSaveParameters
#' @name showSaveParameters
#' @description Description of Giotto saving options, links to \code{\link{all_plots_save_function}}
#' @return Instruction on how to use the automatic plot saving options within Giotto
#' @export
#' @examples
#'     showSaveParameters()
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
#' @param expression_values expression values to use
#' @param genes vector of genes to use, default to 'all'
#' @param cluster_column name of column to use for clusters
#' @param cor correlation score to calculate distance
#' @param distance distance method to use for hierarchical clustering
#' @param show_plot show plot
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters, see \code{\link{showSaveParameters}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @param ... additional parameters for the Heatmap function from ComplexHeatmap
#' @return ggplot
#' @details Correlation heatmap of selected clusters.
#' @export
#' @examples
#'
#' data(mini_giotto_single_cell)
#'
#' # cell metadata
#' cell_metadata = pDataDT(mini_giotto_single_cell)
#'
#' # create heatmap
#' showClusterHeatmap(mini_giotto_single_cell,
#'                    cluster_column = 'cell_types')
#'
showClusterHeatmap <- function(gobject,
                               expression_values = c('normalized', 'scaled', 'custom'),
                               genes = 'all',
                               cluster_column,
                               cor = c('pearson', 'spearman'),
                               distance = 'ward.D',
                               show_plot = NA,
                               return_plot = NA,
                               save_plot = NA,
                               save_param =  list(),
                               default_save_name = 'showClusterHeatmap',
                               ...) {

  ## correlation
  cor = match.arg(cor, c('pearson', 'spearman'))
  values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))

  ## subset expression data
  if(genes[1] != 'all') {
    all_genes = gobject@gene_ID
    detected_genes = genes[genes %in% all_genes]
  } else {
    # NULL = all genes in calculateMetaTable()
    detected_genes = NULL
  }

  metatable = calculateMetaTable(gobject = gobject,
                                 expression_values = values,
                                 metadata_cols = cluster_column,
                                 selected_genes = detected_genes)
  dcast_metatable = data.table::dcast.data.table(metatable, formula = variable~uniq_ID, value.var = 'value')
  testmatrix = dt_to_matrix(x = dcast_metatable)

  # correlation
  cormatrix = cor_giotto(x = testmatrix, method = cor)
  cordist = stats::as.dist(1 - cormatrix, diag = T, upper = T)
  corclus = stats::hclust(d = cordist, method = distance)

  hmap = ComplexHeatmap::Heatmap(matrix = cormatrix,
                                 cluster_rows = corclus,
                                 cluster_columns = corclus, ...)


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
#' @param expression_values expression values to use
#' @param cluster_column name of column to use for clusters
#' @param cor correlation score to calculate distance
#' @param distance distance method to use for hierarchical clustering
#' @param h height of horizontal lines to plot
#' @param h_color color of horizontal lines
#' @param rotate rotate dendrogram 90 degrees
#' @param show_plot show plot
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters, see \code{\link{showSaveParameters}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @param ... additional parameters for ggdendrogram()
#' @return ggplot
#' @details Expression correlation dendrogram for selected clusters.
#' @export
#' @examples
#'
#' data(mini_giotto_single_cell)
#'
#' # cell metadata
#' cell_metadata = pDataDT(mini_giotto_single_cell)
#'
#' # create heatmap
#' showClusterDendrogram(mini_giotto_single_cell,
#'                       cluster_column = 'cell_types')
#'
showClusterDendrogram <- function(gobject,
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

  cor = match.arg(cor, c('pearson', 'spearman'))
  values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))


  metatable = calculateMetaTable(gobject = gobject, expression_values = values, metadata_cols = cluster_column)
  dcast_metatable = data.table::dcast.data.table(metatable, formula = variable~uniq_ID, value.var = 'value')
  testmatrix = dt_to_matrix(x = dcast_metatable)

  # correlation
  cormatrix = cor_giotto(x = testmatrix, method = cor)
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
#' @param expression_values expression values to use
#' @param genes genes to use
#' @param cluster_column name of column to use for clusters
#' @param cluster_order method to determine cluster order
#' @param cluster_custom_order custom order for clusters
#' @param cor_method method for correlation
#' @param hclust_method method for hierarchical clustering
#' @return custom
#' @details Calculates order for clusters.
#' @keywords internal
decide_cluster_order = function(gobject,
                                expression_values = c('normalized', 'scaled', 'custom'),
                                genes,
                                cluster_column = NULL,
                                cluster_order = c('size', 'correlation', 'custom'),
                                cluster_custom_order = NULL,
                                cor_method = 'pearson',
                                hclust_method = 'ward.D') {


  # epxression data
  values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))
  expr_values = select_expression_values(gobject = gobject, values = values)

  # subset expression data
  detected_genes = genes[genes %in% rownames(expr_values)]
  subset_values = expr_values[rownames(expr_values) %in% detected_genes, ]

  # metadata
  cell_metadata = pDataDT(gobject)

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
                                          cluster_column = cluster_column,
                                          gene_subset = detected_genes,
                                          expression_values = values)

    cormatrix = cor_giotto(x = subset_matrix, method = cor_method)
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
#' @param expression_values expression values to use
#' @param genes genes to use
#' @param cluster_column name of column to use for clusters
#' @param cluster_order method to determine cluster order
#' @param cluster_custom_order custom order for clusters
#' @param cluster_cor_method method for cluster correlation
#' @param cluster_hclust_method method for hierarchical clustering of clusters
#' @param gene_order method to determine gene order
#' @param gene_custom_order custom order for genes
#' @param gene_cor_method method for gene correlation
#' @param gene_hclust_method method for hierarchical clustering of genes
#' @return list
#' @details Creates input data.tables for plotHeatmap function.
#' @keywords internal
createHeatmap_DT <- function(gobject,
                             expression_values = c('normalized', 'scaled', 'custom'),
                             genes,
                             cluster_column = NULL,
                             cluster_order = c('size', 'correlation', 'custom'),
                             cluster_custom_order = NULL,
                             cluster_cor_method = 'pearson',
                             cluster_hclust_method = 'ward.D',
                             gene_order = c('correlation', 'custom'),
                             gene_custom_order = NULL,
                             gene_cor_method = 'pearson',
                             gene_hclust_method = 'complete') {

  # epxression data
  values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))
  expr_values = select_expression_values(gobject = gobject, values = values)

  # subset expression data
  detected_genes = genes[genes %in% rownames(expr_values)]
  subset_values = expr_values[rownames(expr_values) %in% detected_genes, ]

  # metadata
  cell_metadata = pDataDT(gobject)

  # gene order
  gene_order = match.arg(gene_order, c('correlation', 'custom'))


  ### cluster order ###
  clus_sort_names = decide_cluster_order(gobject = gobject,
                                         expression_values = expression_values,
                                         genes = genes,
                                         cluster_column = cluster_column,
                                         cluster_order = cluster_order,
                                         cluster_custom_order = cluster_custom_order,
                                         cor_method = cluster_cor_method,
                                         hclust_method = cluster_hclust_method)

  ## data.table ##
  subset_values_DT <- data.table::as.data.table(reshape2::melt(as.matrix(subset_values), varnames = c('genes', 'cells'), value.name = 'expression'))
  subset_values_DT <- merge(subset_values_DT, by.x = 'cells', cell_metadata[, c('cell_ID', cluster_column), with = F], by.y = 'cell_ID')
  subset_values_DT[[cluster_column]] <- factor(subset_values_DT[[cluster_column]], levels = clus_sort_names)

  # data.table variables
  z_scores = scale_scores = V1 = cells = NULL

  subset_values_DT[, genes := factor(genes, unique(detected_genes))]
  subset_values_DT[, z_scores := scale(expression), by = genes]
  subset_values_DT[, scale_scores := scales::rescale(x = expression, to = c(0,1)), by = genes]


  ## order cells by mean expression ##
  cell_order_DT = subset_values_DT[, mean(expression), by = c('cells', cluster_column)]
  cell_order_DT = cell_order_DT[order(get(cluster_column), V1)]
  subset_values_DT[, cells := factor(cells, cell_order_DT$cells)]

  ## get x-coordines for vertical lines in heatmap
  x_lines = cumsum(as.vector(table(cell_order_DT[[cluster_column]])))

  ## order genes ##
  if(gene_order == 'correlation') {
    genesum_per_clus = subset_values_DT[, sum(expression), by = c('genes', cluster_column)]

    my_formula = paste0('genes~',cluster_column)
    test_mat = data.table::dcast.data.table(data = genesum_per_clus, formula = my_formula, value.var = 'V1')
    test_matrix = as.matrix(test_mat[,-1]); rownames(test_matrix) = test_mat$genes

    gene_dist = stats::as.dist(1 - cor_giotto(t_giotto(test_matrix), method = gene_cor_method))
    gene_clus = stats::hclust(gene_dist, method = gene_hclust_method)

    gene_labels = rownames(test_matrix)
    gene_index = 1:length(gene_labels)
    names(gene_index) = gene_labels

    final_gene_order = names(gene_index[match(gene_clus$order, gene_index)])
    subset_values_DT[, genes := factor(genes, final_gene_order)]

  } else if(gene_order == 'custom') {

    if(is.null(gene_custom_order)) {
      stop('\n with custom gene order the gene_custom_order parameter needs to be provided \n')
    }
    subset_values_DT[, genes := factor(genes, gene_custom_order)]

  }

  cell_order_DT[['cells']] = factor(cell_order_DT[['cells']], levels = as.character(cell_order_DT[['cells']]))

  return(list(DT = subset_values_DT, x_lines = x_lines, cell_DT = cell_order_DT))
}



#' @title plotHeatmap
#' @name plotHeatmap
#' @description Creates heatmap for genes and clusters.
#' @param gobject giotto object
#' @param expression_values expression values to use
#' @param genes genes to use
#' @param cluster_column name of column to use for clusters
#' @param cluster_order method to determine cluster order
#' @param cluster_custom_order custom order for clusters
#' @param cluster_color_code color code for clusters
#' @param cluster_cor_method method for cluster correlation
#' @param cluster_hclust_method method for hierarchical clustering of clusters
#' @param gene_order method to determine gene order
#' @param gene_custom_order custom order for genes
#' @param gene_cor_method method for gene correlation
#' @param gene_hclust_method method for hierarchical clustering of genes
#' @param show_values which values to show on heatmap
#' @param size_vertical_lines sizes for vertical lines
#' @param gradient_colors colors for heatmap gradient
#' @param gene_label_selection subset of genes to show on y-axis
#' @param axis_text_y_size size for y-axis text
#' @param legend_nrows number of rows for the cluster legend
#' @param show_plot show plot
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters, see \code{\link{showSaveParameters}}
#' @param default_save_name default save name
#' @return ggplot
#' @details If you want to display many genes there are 2 ways to proceed:
#' \itemize{
#'   \item{1. set axis_text_y_size to a really small value and show all genes}
#'   \item{2. provide a subset of genes to display to gene_label_selection}
#' }
#' @export
#' @examples
#' \dontrun{
#'
#' data(mini_giotto_single_cell)
#'
#' # get all genes
#' all_genes = slot(mini_giotto_single_cell, 'gene_ID')
#'
#' # plot heatmap
#' plotHeatmap(mini_giotto_single_cell,
#'             genes = all_genes[1:10])
#'
#' # look at cell metadata
#' cell_metadata = pDataDT(mini_giotto_single_cell)
#'
#' # plot heatmap per cell type, a column name from cell_metadata
#' plotHeatmap(mini_giotto_single_cell,
#'             genes = all_genes[1:10],
#'             cluster_column = 'cell_types')
#'
#' }
plotHeatmap <- function(gobject,
                        expression_values = c('normalized', 'scaled', 'custom'),
                        genes,
                        cluster_column = NULL,
                        cluster_order = c('size', 'correlation', 'custom'),
                        cluster_custom_order = NULL,
                        cluster_color_code = NULL,
                        cluster_cor_method = 'pearson',
                        cluster_hclust_method = 'ward.D',
                        gene_order = c('correlation', 'custom'),
                        gene_custom_order = NULL,
                        gene_cor_method = 'pearson',
                        gene_hclust_method = 'complete',
                        show_values = c('rescaled', 'z-scaled', 'original'),
                        size_vertical_lines = 1.1,
                        gradient_colors = c('blue', 'yellow', 'red'),
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
                                  expression_values = expression_values,
                                  genes = genes,
                                  cluster_column = cluster_column,
                                  cluster_order = cluster_order,
                                  cluster_custom_order = cluster_custom_order,
                                  cluster_cor_method = cluster_cor_method,
                                  cluster_hclust_method = cluster_hclust_method,
                                  gene_order = gene_order,
                                  gene_custom_order = gene_custom_order,
                                  gene_cor_method = gene_cor_method,
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
  hmap <- hmap + ggplot2::geom_tile(data = subset_values_DT, aes_string(x = 'cells', y = 'genes', fill = value_column))
  hmap <- hmap + ggplot2::geom_vline(xintercept = x_lines, color = 'white', size = size_vertical_lines)
  hmap <- hmap + ggplot2::scale_fill_gradient2(low = low_color, mid = mid_color, high = high_color,
                                               midpoint = midpoint,
                                               guide = ggplot2::guide_colorbar(title = ''))

  if(is.null(gene_label_selection)) {

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
    geneDT = subset_values_DT[,c('genes'), with = F]
    geneDT = unique(setorder(geneDT, genes))

    # data.table variables
    geneOrder = subset_genes = NULL

    geneDT[, geneOrder := 1:.N]
    geneDT[, subset_genes := ifelse(genes %in% gene_label_selection, as.character(genes), '')]

    axis <- ggplot(data = geneDT, aes(x = 0, y = geneOrder, label = subset_genes))
    axis <- axis + ggrepel::geom_text_repel(min.segment.length = grid::unit(0, "pt"),
                                            color = "grey30",  ## ggplot2 theme_grey() axis text
                                            size = axis_text_y_size  ## default ggplot2 theme_grey() axis text
    )
    axis <- axis + ggplot2::scale_x_continuous(limits = c(0, 1), expand = c(0, 0),
                                      breaks = NULL, labels = NULL, name = NULL)
    axis <- axis + ggplot2::scale_y_continuous(limits = c(0, nrow(geneDT)), expand = c(0, 0),
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
#' @param expression_values expression values to use
#' @param metadata_cols annotation columns found in pDataDT(gobject)
#' @param selected_genes subset of genes to use
#' @param first_meta_col if more than 1 metadata column, select the x-axis factor
#' @param second_meta_col if more than 1 metadata column, select the facetting factor
#' @param show_values which values to show on heatmap
#' @param custom_cluster_order custom cluster order (default = NULL)
#' @param clus_cor_method correlation method for clusters
#' @param clus_cluster_method hierarchical cluster method for the clusters
#' @param custom_gene_order custom gene order (default = NULL)
#' @param gene_cor_method correlation method for genes
#' @param gene_cluster_method hierarchical cluster method for the genes
#' @param gradient_color vector with 3 colors for numeric data
#' @param gradient_midpoint midpoint for color gradient
#' @param gradient_limits vector with lower and upper limits
#' @param x_text_size size of x-axis text
#' @param x_text_angle angle of x-axis text
#' @param y_text_size size of y-axis text
#' @param strip_text_size size of strip text
#' @param show_plot show plot
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters, see \code{\link{showSaveParameters}}
#' @param default_save_name default save name
#' @return ggplot or data.table
#' @details Creates heatmap for the average expression of selected genes in the different annotation/cluster groups.
#' Calculation of cluster or gene order is done on the provided expression values, but visualization
#' is by default on the z-scores. Other options are the original values or z-scores rescaled per gene (-1 to 1).
#' @seealso \code{\link{plotMetaDataCellsHeatmap}} for numeric cell annotation instead of gene expression.
#' @export
#' @examples
#' \dontrun{
#'
#' data(mini_giotto_single_cell)
#'
#' # get all genes
#' all_genes = slot(mini_giotto_single_cell, 'gene_ID')
#'
#' # look at cell metadata
#' cell_metadata = pDataDT(mini_giotto_single_cell)
#'
#' # plot heatmap per cell type, a column name from cell_metadata
#' plotMetaDataHeatmap(mini_giotto_single_cell,
#'                     selected_genes = all_genes[1:10],
#'                     metadata_cols = 'cell_types')
#'
#' }
plotMetaDataHeatmap = function(gobject,
                               expression_values =  c("normalized", "scaled", "custom"),
                               metadata_cols = NULL,
                               selected_genes = NULL,
                               first_meta_col = NULL,
                               second_meta_col = NULL,
                               show_values = c('zscores', 'original', 'zscores_rescaled'),
                               custom_cluster_order = NULL,
                               clus_cor_method = 'pearson',
                               clus_cluster_method = 'complete',
                               custom_gene_order = NULL,
                               gene_cor_method = 'pearson',
                               gene_cluster_method = 'complete',
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


  metaDT = calculateMetaTable(gobject = gobject,
                              expression_values = expression_values,
                              metadata_cols = metadata_cols,
                              selected_genes = selected_genes)

  # data.table variables
  zscores = value = zscores_rescaled_per_gene = NULL

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
    cormatrix = cor_giotto(x = testmain_mat, method = clus_cor_method)
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
  if(is.null(custom_gene_order)) {
    gene_cormatrix = cor_giotto(x = t(testmain_mat), method = gene_cor_method)
    gene_cordist = stats::as.dist(1 - gene_cormatrix, diag = T, upper = T)
    gene_corclus = stats::hclust(d = gene_cordist, method = gene_cluster_method)
    gene_names = rownames(gene_cormatrix)
    names(gene_names) = 1:length(gene_names)
    gene_sort_names = gene_names[gene_corclus$order]
  }  else {
    gene_sort_names = unique(as.character(custom_gene_order))
    if(all(rownames(testmain_mat) %in% gene_sort_names) == FALSE) {
      stop('\n custom gene order is given, but not all genes are represented \n')
    }
  }

  if(length(metadata_cols) == 1) {

    # data.table variables
    factor_column = variable = NULL

    metaDT[, factor_column := factor(get(metadata_cols), levels = clus_sort_names)]
    metaDT[, variable := factor(get('variable'), levels = gene_sort_names)]

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
      metaDT[, variable := factor(get('variable'), levels = gene_sort_names)]

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


#' @title plotMetaDataCellsHeatmap
#' @name plotMetaDataCellsHeatmap
#' @description Creates heatmap for numeric cell metadata within aggregated clusters.
#' @param gobject giotto object
#' @param metadata_cols annotation columns found in pDataDT(gobject)
#' @param spat_enr_names spatial enrichment results to include
#' @param value_cols value columns to use
#' @param first_meta_col if more than 1 metadata column, select the x-axis factor
#' @param second_meta_col if more than 1 metadata column, select the facetting factor
#' @param show_values which values to show on heatmap
#' @param custom_cluster_order custom cluster order (default = NULL)
#' @param clus_cor_method correlation method for clusters
#' @param clus_cluster_method hierarchical cluster method for the clusters
#' @param custom_values_order custom values order (default = NULL)
#' @param values_cor_method correlation method for values
#' @param values_cluster_method hierarchical cluster method for the values
#' @param midpoint midpoint of show_values
#' @param x_text_size size of x-axis text
#' @param x_text_angle angle of x-axis text
#' @param y_text_size size of y-axis text
#' @param strip_text_size size of strip text
#' @param show_plot show plot
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters, see \code{\link{showSaveParameters}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return ggplot or data.table
#' @details Creates heatmap for the average values of selected value columns in the different annotation groups.
#' @seealso \code{\link{plotMetaDataHeatmap}} for gene expression instead of numeric cell annotation data.
#' @export
plotMetaDataCellsHeatmap = function(gobject,
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
    cormatrix = cor_giotto(x = testmain_mat, method = clus_cor_method)
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
    values_cormatrix = cor_giotto(x = t(testmain_mat), method = values_cor_method)
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
#' @param expression_values expression values to use
#' @param genes genes to plot
#' @param cluster_column name of column to use for clusters
#' @param color_violin color violinplots according to genes or clusters
#' @param cluster_custom_order custom order of clusters
#' @param color_violin color violin according to genes or clusters
#' @param cluster_color_code color code for clusters
#' @param strip_position position of gene labels
#' @param strip_text size of strip text
#' @param axis_text_x_size size of x-axis text
#' @param axis_text_y_size size of y-axis text
#' @param show_plot show plot
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters, see \code{\link{showSaveParameters}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @return ggplot
#' @export
#' @examples
#' \dontrun{
#'
#' data(mini_giotto_single_cell)
#'
#' # get all genes
#' all_genes = slot(mini_giotto_single_cell, 'gene_ID')
#'
#' # look at cell metadata
#' cell_metadata = pDataDT(mini_giotto_single_cell)
#'
#' # plot violinplot with selected genes and stratified for identified cell types
#' violinPlot(mini_giotto_single_cell,
#'            genes = all_genes[1:10],
#'            cluster_column = 'cell_types')
#'
#' }
violinPlot <- function(gobject,
                       expression_values = c('normalized', 'scaled', 'custom'),
                       genes,
                       cluster_column,
                       cluster_custom_order = NULL,
                       color_violin = c('genes', 'cluster'),
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

  ## strip position
  strip_position = match.arg(strip_position, c('top', 'right', 'left', 'bottom'))

  ## color of violin plots
  color_violin = match.arg(color_violin, c('genes', 'cluster'))

  ## expression data ##
  values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))
  expr_data = select_expression_values(gobject = gobject, values = values)

  # only keep genes that are in the dataset
  selected_genes = genes[genes %in% rownames(expr_data)]
  if(length(selected_genes[duplicated(selected_genes)]) != 0) {
    cat('These genes have duplicates: \n',
        selected_genes[duplicated(selected_genes)])

    selected_genes = unique(selected_genes)
  }
  subset_data = as.matrix(expr_data[rownames(expr_data) %in% selected_genes, ])

  if(length(genes) == 1) {
    t_subset_data = subset_data
  } else {
    t_subset_data = t_giotto(subset_data)
  }

  # metadata
  metadata = pDataDT(gobject)

  if(length(genes) == 1) {
    metadata_expr <- cbind(metadata,  t_subset_data)
    setnames(metadata_expr, 'V1', genes)
  } else {
    metadata_expr <- cbind(metadata, t_subset_data)
  }


  metadata_expr_m <- data.table::melt.data.table(metadata_expr, measure.vars = unique(selected_genes), variable.name = 'genes')
  metadata_expr_m[, genes := factor(genes, selected_genes)]
  metadata_expr_m[[cluster_column]] = as.factor(metadata_expr_m[[cluster_column]])

  if(!is.null(cluster_custom_order)) {
    cluster_custom_order = unique(cluster_custom_order)
    metadata_expr_m[[cluster_column]] = factor(x = metadata_expr_m[[cluster_column]], levels = cluster_custom_order)
  }


  pl <- ggplot2::ggplot()
  pl <- pl + ggplot2::theme_classic()

  if(color_violin == 'genes') {
    pl <- pl + ggplot2::geom_violin(data = metadata_expr_m, aes_string(x = cluster_column, y = 'value', fill = 'genes'), width = 1, scale = 'width', show.legend = F)
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


  pl <- pl + ggplot2::facet_wrap(.~genes, ncol = 1, strip.position = strip_position)
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
#' @param gobject network in giotto object
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
#' @param spatial_grid spatial_grid in giotto object
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











