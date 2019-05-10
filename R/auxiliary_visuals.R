



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
  cormatrix = cor(x = testmatrix, method = distance)
  cordist = as.dist(1 - cormatrix, diag = T, upper = T)
  corclus = hclust(d = cordist, method = 'ward.D')

  hmap = ComplexHeatmap::Heatmap(matrix = cormatrix,
                 cluster_rows = corclus,
                 cluster_columns = corclus)

  return(hmap)

}


#' @title plotHeatmap
#' @name plotHeatmap
#' @description Creates heatmap based on identified clusters
#' @param gobject giotto object
#' @param expression_values expression values to use
#' @param genes genes to plot
#' @param cluster_column name of column to use for clusters
#' @param cluster_order method to order the clusters (see details)
#' @param cluster_custom_order vector with custom order of clusters
#' @param gene_order order of genes (see details)
#' @param show_values values to plot on heatmap (see details)
#' @return ggplot
#' @details Correlation heatmap of clusters vs genes.
#' @export
#' @examples
#'     plotHeatmap(gobject)
plotHeatmap <- function(gobject,
                        expression_values = c('normalized', 'scaled', 'custom'),
                        genes,
                        cluster_column = NULL,
                        cluster_order = c('size', 'correlation', 'custom'),
                        cluster_custom_order = NULL,
                        gene_order = c('custom', 'correlation'),
                        show_values = c('rescaled', 'z-scaled', 'original')) {

  show_values = match.arg(show_values, choices = c('rescaled', 'z-scaled', 'original'))

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


  # cluster order
  if(!is.null(cluster_column)) {

    if(!cluster_column %in% colnames(cell_metadata)) {
      stop('\n cluster column is not found \n')
    }

    ## cluster order ##
    cluster_order = match.arg(cluster_order, c('size', 'correlation', 'custom'))

    if(cluster_order == 'size') {
      clus_sort_sizes = sort(table(cell_metadata[[cluster_column]]))
      clus_sort_names = names(clus_sort_sizes)

    } else if(cluster_order == 'correlation') {

      subset_matrix = create_cluster_matrix(gobject = gobject,
                                            cluster_column = cluster_column,
                                            gene_subset = detected_genes,
                                            expression_values = values)

      cormatrix = cor(x = subset_matrix, method = 'pearson')
      cordist = as.dist(1 - cormatrix, diag = T, upper = T)
      corclus = hclust(d = cordist, method = 'ward.D')

      clus_names = rownames(cormatrix)
      names(clus_names) = 1:length(clus_names)
      clus_sort_names = clus_names[corclus$order]

      clus_sort_names = gsub(pattern = 'cluster_', replacement = '', x = clus_sort_names)

    } else if(cluster_order == 'custom') {
      if(is.null(cluster_custom_order)) {
        stop('\n custom order parameter is not given \n')
      }
      clus_sort_names = cluster_custom_order
    }

    ## data.table ##
    subset_values_DT <- as.data.table(melt(subset_values, varnames = c('genes', 'cells'), value.name = 'expression'))
    subset_values_DT <- merge(subset_values_DT, by.x = 'cells', cell_metadata[, c('cell_ID', cluster_column), with = F], by.y = 'cell_ID')
    subset_values_DT[[cluster_column]] <- factor(subset_values_DT[[cluster_column]], levels = clus_sort_names)

    subset_values_DT[, genes := factor(genes, unique(detected_genes))]
    subset_values_DT[, z_scores := base::scale(expression), by = genes]
    subset_values_DT[, scale_scores := scales::rescale(x = expression, to = c(0,1)), by = genes]

    # order cells
    cell_order_DT = subset_values_DT[, mean(expression), by = c('cells', cluster_column)]
    cell_order_DT = cell_order_DT[order(get(cluster_column), V1)]
    subset_values_DT[, cells := factor(cells, cell_order_DT$cells)]

    clus_values = unique(cell_order_DT[[cluster_column]])
    clus_colors = getDistinctColors(n = length(clus_values))
    names(clus_colors) = clus_values

    x_lines = cumsum(as.vector(table(cell_order_DT[[cluster_column]])))



    # order genes
    if(gene_order == 'correlation') {
      genesum_per_clus = subset_values_DT[, sum(expression), by = c('genes', cluster_column)]

      my_formula = paste0('genes~',cluster_column)
      test_mat = dcast.data.table(data = genesum_per_clus, formula = my_formula, value.var = 'V1')
      test_matrix = as.matrix(test_mat[,-1]); rownames(test_matrix) = test_mat$genes

      gene_dist = dist(test_matrix, method = 'euclidean')
      gene_dist = as.dist(1-cor(t(test_matrix), method = 'pearson'))
      gene_clus = hclust(gene_dist, method = 'complete')

      gene_labels = rownames(test_matrix)
      gene_index = 1:length(gene_labels)
      names(gene_index) = gene_labels

      final_gene_order = names(gene_index[match(gene_clus$order, gene_index)])
      subset_values_DT[, genes := factor(genes, final_gene_order)]
    }

    cell_order_DT[['cells']] = factor(cell_order_DT[['cells']], levels = as.character(cell_order_DT[['cells']]))

    clus_pl <- ggplot()
    clus_pl <- clus_pl + geom_raster(data = cell_order_DT, aes_string(x = 'cells', y = '1', fill = cluster_column))
    clus_pl <- clus_pl + geom_vline(xintercept = x_lines, color = 'white', size = 1.1)
    clus_pl <- clus_pl + scale_fill_manual(values = clus_colors, guide = guide_legend(title = '', nrow = 1))
    clus_pl <- clus_pl + theme(axis.text = element_blank(),
                               axis.ticks = element_blank(),
                               axis.line = element_blank(),
                               legend.position = 'top')
    clus_pl <- clus_pl + labs(x = '', y = '')
    clus_pl


    # rescale values
    pl <- ggplot()
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
    pl <- pl + geom_raster(data = subset_values_DT, aes_string(x = 'cells', y = 'genes', fill = value_column))
    pl <- pl + geom_vline(xintercept = x_lines, color = 'white', size = 1.1)
    pl <- pl + theme(axis.text.x = element_blank(),
                     axis.ticks.x = element_blank())
    pl <- pl + scale_fill_gradient2(low = 'blue', mid = 'yellow', high = 'red', midpoint = midpoint)
    pl


    combo_plot <- cowplot::plot_grid(clus_pl, pl, ncol = 1, rel_heights = c(0.5, 2), rel_widths = c(1), align = 'v', axis = 'rlbt')
    print(cowplot::plot_grid(combo_plot))


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
                       color_violin = c('genes', 'cluster'),
                       strip.text = 7,
                       axis.text.x.size = 10,
                       axis.text.y.size = 6) {


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


  metadata_expr_m <- melt.data.table(metadata_expr, measure.vars = unique(selected_genes), variable.name = 'genes')
  metadata_expr_m[, genes := factor(genes, selected_genes)]
  metadata_expr_m[[cluster_column]] = as.factor(metadata_expr_m[[cluster_column]])



  pl <- ggplot()
  pl <- pl + theme_classic()

  if(color_violin == 'genes') {
    pl <- pl + geom_violin(data = metadata_expr_m, aes_string(x = cluster_column, y = 'value', fill = 'genes'), width = 1, scale = 'width', show.legend = F)
  } else {
    pl <- pl + geom_violin(data = metadata_expr_m, aes_string(x = cluster_column, y = 'value', fill = cluster_column), width = 1, scale = 'width', show.legend = F)
  }
  pl <- pl + facet_wrap(~genes, ncol = 1)
  pl <- pl + theme(strip.text = element_text(size = strip.text),
                   axis.text.x = element_text(size = axis.text.x.size, angle = 45, hjust = 1, vjust = 1),
                   axis.text.y = element_text(size = axis.text.y.size))
  pl <- pl + labs(x = '', y = 'normalized expression')
  pl

  return(pl)


}

