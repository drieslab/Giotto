




#' @title showClusterHeatmap
#' @name showClusterHeatmap
#' @description Creates heatmap based on identified clusters
#' @param gobject giotto object
#' @param expression_values expression values to use
#' @param cluster_column name of column to use for clusters
#' @param distance correlation score to calculate distance
#' @return ggplot
#' @details Correlation heatmap of clusters.
#' @export
#' @examples
#'     showClusterHeatmap(gobject)
showClusterHeatmap <- function(gobject,
                               expression_values = c('normalized', 'scaled', 'custom'),
                               cluster_column,
                               distance = c('pearson', 'spearman')) {

  distance = match.arg(distance, c('pearson', 'spearman'))
  values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))


  testmatrix = create_cluster_matrix(gobject = gobject,
                                     cluster_column = cluster_column,
                                     expression_values = values)

  # correlation
  cormatrix = stats::cor(x = testmatrix, method = distance)
  cordist = stats::as.dist(1 - cormatrix, diag = T, upper = T)
  corclus = stats::hclust(d = cordist, method = 'ward.D')

  hmap = ComplexHeatmap::Heatmap(matrix = cormatrix,
                 cluster_rows = corclus,
                 cluster_columns = corclus)

  return(hmap)

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
#' @examples
#'     decide_cluster_order(gobject)
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
    cormatrix = stats::cor(x = subset_matrix, method = cor_method)
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
#' @param gene_cor_method method for gene correlation
#' @param gene_hclust_method method for hierarchical clustering of genes
#' @return list
#' @details Creates input data.tables for plotHeatmap function.
#' @examples
#'     createHeatmap_DT(gobject)
createHeatmap_DT <- function(gobject,
                             expression_values = c('normalized', 'scaled', 'custom'),
                             genes,
                             cluster_column = NULL,
                             cluster_order = c('size', 'correlation', 'custom'),
                             cluster_custom_order = NULL,
                             cluster_cor_method = 'pearson',
                             cluster_hclust_method = 'ward.D',
                             gene_order = c('custom', 'correlation'),
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
  gene_order = match.arg(gene_order, c('custom', 'correlation'))


  ### cluster order ###
  clus_sort_names = decide_cluster_order(gobject = gobject,
                                         expression_values = expression_values,
                                         genes = genes,
                                         cluster_column = cluster_column,
                                         cluster_order = cluster_order,
                                         cor_method = cluster_cor_method,
                                         hclust_method = cluster_hclust_method)

  ## data.table ##
  subset_values_DT <- data.table::as.data.table(melt(subset_values, varnames = c('genes', 'cells'), value.name = 'expression'))
  subset_values_DT <- merge(subset_values_DT, by.x = 'cells', cell_metadata[, c('cell_ID', cluster_column), with = F], by.y = 'cell_ID')
  subset_values_DT[[cluster_column]] <- factor(subset_values_DT[[cluster_column]], levels = clus_sort_names)

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

    gene_dist = stats::as.dist(1-cor(t(test_matrix), method = gene_cor_method))
    gene_clus = stats::hclust(gene_dist, method = gene_hclust_method)

    gene_labels = rownames(test_matrix)
    gene_index = 1:length(gene_labels)
    names(gene_index) = gene_labels

    final_gene_order = names(gene_index[match(gene_clus$order, gene_index)])
    subset_values_DT[, genes := factor(genes, final_gene_order)]
  }

  cell_order_DT[['cells']] = factor(cell_order_DT[['cells']], levels = as.character(cell_order_DT[['cells']]))

  return(list(DT = subset_values_DT, x_lines = x_lines, cell_DT = cell_order_DT))
}



#' @title plotHeatmap
#' @name plotHeatmap
#' @description creates order for clusters
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
#' @param gene_cor_method method for gene correlation
#' @param gene_hclust_method method for hierarchical clustering of genes
#' @param show_values which values to show on heatmap
#' @param size_vertical_lines sizes for vertical lines
#' @param gradient_colors colors for heatmap gradient
#' @param gene_label_selection subset of genes to show on y-axis
#' @param axis_text_y_size size for y-axis text
#' @param legend_nrows number of rows for the cluster legend
#' @return ggplot
#' @details Creates heatmap for genes and clusters.
#' @export
#' @examples
#'     plotHeatmap(gobject)
plotHeatmap <- function(gobject,
                        expression_values = c('normalized', 'scaled', 'custom'),
                        genes,
                        cluster_column = NULL,
                        cluster_order = c('size', 'correlation', 'custom'),
                        cluster_custom_order = NULL,
                        cluster_color_code = NULL,
                        cluster_cor_method = 'pearson',
                        cluster_hclust_method = 'ward.D',
                        gene_order = c('custom', 'correlation'),
                        gene_cor_method = 'pearson',
                        gene_hclust_method = 'complete',
                        show_values = c('rescaled', 'z-scaled', 'original'),
                        size_vertical_lines = 1.1,
                        gradient_colors = c('blue', 'yellow', 'red'),
                        gene_label_selection = NULL,
                        axis_text_y_size = NULL,
                        legend_nrows = 1) {

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
                                  gene_cor_method = gene_cor_method,
                                  gene_hclust_method = gene_hclust_method)

  cell_order_DT = heatmap_data[['cell_DT']]
  subset_values_DT = heatmap_data[['DT']]
  x_lines = heatmap_data[['x_lines']]

  ## assign colors to each cluster
  if(is.null(cluster_color_code)) {
    clus_values = unique(cell_order_DT[[cluster_column]])
    clus_colors = Giotto:::getDistinctColors(n = length(clus_values))
    names(clus_colors) = clus_values
  } else {
    clus_colors = cluster_color_code
  }

  ## bar on top ##
  clus_pl <- ggplot2::ggplot()
  clus_pl <- clus_pl + ggplot2::geom_raster(data = cell_order_DT, aes_string(x = 'cells', y = '1', fill = cluster_column))
  clus_pl <- clus_pl + ggplot2::geom_vline(xintercept = x_lines, color = 'white', size = size_vertical_lines)
  clus_pl <- clus_pl + ggplot2::scale_fill_manual(values = clus_colors, guide = guide_legend(title = '', nrow = legend_nrows))
  clus_pl <- clus_pl + ggplot2::theme(axis.text = element_blank(),
                                      axis.ticks = element_blank(),
                                      axis.line = element_blank(),
                                      legend.position = 'top',
                                      plot.margin = margin(0, 0, 0, 5.5, "pt"))
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
  empty = ggplot()
  empty = empty + theme_classic()
  empty = empty  + ggplot2::theme(axis.text = element_blank(),
                                  axis.ticks = element_blank(),
                                  axis.line = element_blank(),
                                  plot.margin = margin(0, 0, 0, 0, "cm"))

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
                                               guide = guide_colorbar(title = ''))

  if(is.null(gene_label_selection)) {

    hmap <- hmap + ggplot2::theme(axis.text.x = element_blank(),
                                  axis.ticks.x = element_blank(),
                                  axis.text.y = element_text(size = axis_text_y_size),
                                  plot.margin = margin(0, 0, 5.5, 5.5, "pt"))

    ## align and combine
    aligned <- cowplot::align_plots(clus_pl, empty, hmap + theme(legend.position = "none"),  align = "v", axis = "l")
    aligned <- append(aligned, list(cowplot::get_legend(hmap)))
    cowplot::plot_grid(plotlist = aligned,
                       ncol = 2, rel_widths = c(1, 0.2),
                       nrow = 2, rel_heights = c(0.2, 1))


  } else {

    # set defaults
    if(is.null(axis_text_y_size)) axis_text_y_size = 0.8*11/ggplot2::.pt

    # finish heatmap
    hmap <- hmap + ggplot2::theme(axis.text = element_blank(),
                                  axis.ticks = element_blank(),
                                  plot.margin = margin(0, 0, 5.5, 5.5, "pt"))
    ### axis ###
    geneDT = subset_values_DT[,c('genes'), with = F]
    geneDT = unique(setorder(geneDT, genes))
    geneDT[, geneOrder := 1:.N]
    geneDT[, subset_genes := ifelse(genes %in% gene_label_selection, as.character(genes), '')]

    axis <- ggplot(data = geneDT, aes(x = 0, y = geneOrder, label = subset_genes))
    axis <- axis + ggrepel::geom_text_repel(min.segment.length = grid::unit(0, "pt"),
                                            color = "grey30",  ## ggplot2 theme_grey() axis text
                                            size = axis_text_y_size  ## default ggplot2 theme_grey() axis text
    )
    axis <- axis + scale_x_continuous(limits = c(0, 1), expand = c(0, 0),
                                      breaks = NULL, labels = NULL, name = NULL)
    axis <- axis + scale_y_continuous(limits = c(0, nrow(geneDT)), expand = c(0, 0),
                                      breaks = NULL, labels = NULL, name = NULL)
    axis <- axis + theme(panel.background = element_blank(),
                         plot.margin = margin(0, 0, 0, 0, "cm"))


    #return(list(hmap, clus_pl, axis, empty))

    ## align and combine
    aligned <- cowplot::align_plots(clus_pl, empty, empty, hmap + theme(legend.position = "none"), axis, align = "h", axis = "b")
    aligned <- append(aligned, list(cowplot::get_legend(hmap)))
    cowplot::plot_grid(plotlist = aligned,
                       ncol = 3, rel_widths = c(1, 0.2, 0.1),
                       nrow = 2, rel_heights = c(0.2, 1))



  }

}



#' @title violinPlot
#' @name violinPlot
#' @description Creates heatmap based on identified clusters
#' @param gobject giotto object
#' @param expression_values expression values to use
#' @param genes genes to plot
#' @param cluster_column name of column to use for clusters
#' @param color_violin color violinplots according to genes or clusters
#' @return ggplot
#' @details Correlation heatmap of clusters vs genes.
#' @export
#' @examples
#'     violinPlot(gobject)
violinPlot <- function(gobject,
                       expression_values = c('normalized', 'scaled', 'custom'),
                       genes,
                       cluster_column,
                       cluster_custom_order = NULL,
                       color_violin = c('genes', 'cluster'),
                       cluster_color_code = NULL,
                       strip_text = 7,
                       axis_text_x_size = 10,
                       axis_text_y_size = 6) {


  ## color of violin plots
  color_violin = match.arg(color_violin, c('genes', 'cluster'))

  ## expression data ##
  values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))
  expr_data = select_expression_values(gobject = gobject, values = values)

  # only keep genes that are in the dataset
  selected_genes = genes[genes %in% rownames(expr_data) ]
  subset_data = as.matrix(expr_data[rownames(expr_data) %in% selected_genes, ])

  if(length(genes) == 1) {
    t_subset_data = subset_data
  } else {
    t_subset_data = t(subset_data)
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
  pl <- pl + ggplot2::facet_wrap(~genes, ncol = 1)
  pl <- pl + ggplot2::theme(strip.text = element_text(size = strip_text),
                   axis.text.x = element_text(size = axis_text_x_size, angle = 45, hjust = 1, vjust = 1),
                   axis.text.y = element_text(size = axis_text_y_size))
  pl <- pl + ggplot2::labs(x = '', y = 'normalized expression')
  pl

  return(pl)


}


#' @title plotly_network
#' @name plotly_network
#' @description provide network segment to draw in 3D plot_ly()
#' @param gobject network in giotto object
#' @return edges in network as data.table()
#' @export
#' @examples
#'     plotly_network(gobject)
plotly_network <- function(network,
                           x = "sdimx_begin",y = "sdimy_begin",z = "sdimz_begin",
                           x_end = "sdimx_end",y_end="sdimy_end",z_end="sdimz_end"){
  edges <- data.table(edge_id = 1:(3*dim(network)[1]),x = 0,y = 0,z = 0)
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
#' @export
#' @examples
#'     plotly_grid(gobject)
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

  edges <- data.table(edge_id = 1:edge_num,x = 0,y = 0,x_end = 0,y_end = 0)

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
