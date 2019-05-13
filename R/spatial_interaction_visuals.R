



#' @title cellProximityBarplot
#' @name cellProximityBarplot
#' @description Create barplot from cell-cell proximity scores
#' @param CPscore CPscore, output from cellProximityEnrichment()
#' @param min_orig_ints filter on minimum original cell-cell interactions
#' @param min_sim_ints filter on minimum simulated cell-cell interactions
#' @param p_val p-value
#' @return ggplot barplot
#' @details Give more details ...
#' @export
#' @examples
#'     cellProximityBarplot(CPscore)
cellProximityBarplot = function(CPscore,
                                min_orig_ints = 5,
                                min_sim_ints = 5,
                                p_val = 0.05) {


  table_mean_results_dc = CPscore$enrichm_res

  ## filter to remove low number of cell-cell proximity interactions ##
  table_mean_results_dc_filter = table_mean_results_dc[original >= min_orig_ints & simulations >= min_sim_ints]
  table_mean_results_dc_filter = table_mean_results_dc_filter[p_higher_orig <= p_val | p_lower_orig <= p_val]

  pl <- ggplot()
  pl <- pl + geom_bar(data = table_mean_results_dc_filter, aes(x = unified_int, y = enrichm, fill = type_int), stat = 'identity', show.legend = F)
  pl <- pl + coord_flip()
  pl <- pl + theme_bw()
  pl <- pl + labs(y = 'enrichment/depletion')
  pl

  bpl <- ggplot()
  bpl <- bpl + geom_bar(data = table_mean_results_dc_filter, aes(x = unified_int, y = original, fill = type_int), stat = 'identity', show.legend = T)
  bpl <- bpl + coord_flip()
  bpl <- bpl + theme_bw() + theme(axis.text.y = element_blank())
  bpl <- bpl + labs(y = '# of interactions')
  bpl

  combo_plot <- cowplot::plot_grid(pl, bpl, ncol = 2, rel_heights = c(1), rel_widths = c(3,1.5), align = 'h')
  print(cowplot::plot_grid(combo_plot))

}

#' @title cellProximityBarplot
#' @name cellProximityHeatmap
#' @description Create heatmap from cell-cell proximity scores
#' @param CPscore CPscore, output from cellProximityEnrichment()
#' @param scale scale cell-cell proximity interaction scores
#' @param order_cell_types order
#' @return ggplot heatmap
#' @details Give more details ...
#' @export
#' @examples
#'     cellProximityHeatmap(CPscore)
cellProximityHeatmap = function(CPscore,
                                scale = T,
                                order_cell_types = T) {


  enrich_res = CPscore$enrichm_res
  enrich_res[, first_type := strsplit(x = as.character(unified_int), split = '-')[[1]][1], by = 1:nrow(enrich_res)]
  enrich_res[, second_type := strsplit(x = as.character(unified_int), split = '-')[[1]][2], by = 1:nrow(enrich_res)]

  # create matrix
  enrich_mat = dcast.data.table(data = enrich_res,formula = first_type~second_type, value.var = 'enrichm')
  matrix_d <- as.matrix(enrich_mat[,-1]); rownames(matrix_d) = as.vector(enrich_mat[[1]])
  t_matrix_d <- t(matrix_d)
  t_matrix_d[upper.tri(t_matrix_d)] <- matrix_d[upper.tri(matrix_d)]
  final_matrix = t_matrix_d

  # scale data
  if(scale == TRUE) {
    final_matrix <- t(scale(t(final_matrix)))
    final_matrix <- t(final_matrix)
    final_matrix[upper.tri(final_matrix)] <- final_matrix[upper.tri(final_matrix)]
  }

  # if NA values, impute as mean
  if(any(is.na(final_matrix)) == TRUE) {
    myrowmeans = apply(X = final_matrix, MARGIN = 1, FUN = function(x) mean(na.omit(x)))
    mymatrixmeans = matrix(data = rep(myrowmeans, ncol(final_matrix)), nrow = nrow(final_matrix), ncol = ncol(final_matrix))
    final_matrix[is.na(final_matrix)] = mymatrixmeans[which(is.na(final_matrix))]
  }

  # order cell types
  if(order_cell_types == TRUE) {

    cordist = as.dist(1-cor(final_matrix))
    clus = hclust(cordist)
    myorder = clus$order
    mylabels = clus$labels
    names(mylabels) = 1:length(mylabels)
    sample_order = mylabels[myorder]

    final_matrix = final_matrix[sample_order, sample_order]
  }

  ComplexHeatmap::Heatmap(matrix = final_matrix, cluster_rows = F, cluster_columns = F)

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
  specific_int_scores = melt.data.table(data = specific_int_scores, id.vars = c('genes', 'interaction'))
  specific_int_scores[, value := as.numeric(value)]
  specific_int_scores[, interaction := factor(interaction, order_ints)]
  specific_int_scores[, score := ifelse(grepl('spat', variable), 'difference',
                                        ifelse(grepl('all', variable), 'all cells', 'int cells'))]
  specific_int_scores[, variable := factor(variable, levels = c('cell_expr_1', 'cell_expr_2', 'comb_expr',
                                                                'all_cell_expr_1', 'all_cell_expr_2', 'all_comb_expr',
                                                                'diff_spat_1', 'diff_spat_2', 'diff_spat'))]

  pl <- ggplot()
  pl <- pl + theme_classic()
  pl <- pl + geom_bar(data = specific_int_scores, aes(x = interaction, y = as.numeric(value), fill = score), stat = 'identity')
  pl <- pl + facet_wrap(~ variable, ncol = 9)
  pl <- pl + coord_flip()
  pl <- pl + labs(x = '', y = 'normalized expression', title = selected_gene)
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
  specific_gene_scores = melt.data.table(data = specific_gene_scores, id.vars = c('genes', 'interaction'))
  specific_gene_scores[, value := as.numeric(value)]
  specific_gene_scores[, genes := factor(genes, order_genes)]
  specific_gene_scores[, score := ifelse(grepl('spat', variable), 'difference',
                                         ifelse(grepl('all', variable), 'all cells', 'int cells'))]
  specific_gene_scores[, variable := factor(variable, levels = c('cell_expr_1', 'cell_expr_2', 'comb_expr',
                                                                 'all_cell_expr_1', 'all_cell_expr_2', 'all_comb_expr',
                                                                 'diff_spat_1', 'diff_spat_2', 'diff_spat'))]

  specific_gene_scores = specific_gene_scores[genes %in% c(head_genes, tail_genes)]

  pl <- ggplot()
  pl <- pl + theme_classic()
  pl <- pl + geom_bar(data = specific_gene_scores, aes(x = genes, y = as.numeric(value), fill = score), stat = 'identity')
  pl <- pl + facet_wrap(~ variable, ncol = 9)
  pl <- pl + coord_flip()
  pl <- pl + labs(x = '', y = 'normalized expression', title = '1-3')
  pl

  print(pl)

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
#' @param show.legend show legend
#' @return ggplot
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
                                 show_network = F,
                                 network_color = NULL,
                                 spatial_network_name = 'spatial_network',
                                 show_grid = F,
                                 grid_color = NULL,
                                 spatial_grid_name = 'spatial_grid',
                                 coord_fix_ratio = 0.6,
                                 show.legend = T,
                                 point_size_select = 2,
                                 point_size_other = 1,
                                 ...) {


  if(is.null(interaction_name)) {
    stop('\n you need to specific at least one interaction name, run cellProximityEnrichment \n')
  }


  cell_locations  = gobject@spatial_locs
  spatial_grid    = gobject@spatial_grid[[spatial_grid_name]]
  cell_metadata   = gobject@cell_metadata
  cell_metadata   = cell_metadata[, !grepl('cell_ID', colnames(cell_metadata)), with = F]


  spatial_network = annotateSpatialNetwork(gobject = gobject, spatial_network_name = spatial_network_name, cluster_column = cluster_column)

  cell_IDs_to_keep = unique(c(spatial_network[unified_int %in% interaction_name]$to, spatial_network[unified_int %in% interaction_name]$from))


  # annotated cell data
  if(nrow(cell_metadata) == 0) {
    cell_locations_metadata = cell_locations
  } else {
    cell_locations_metadata <- cbind(cell_locations, cell_metadata)
  }




  # first 2 dimensions need to be defined
  if(is.null(sdimx) | is.null(sdimy)) {
    cat('first and second dimenion need to be defined, default is first 2 \n')
    sdimx = 'sdimx'
    sdimy = 'sdimy'
  }



  # if 3 dimensions are defined create a 3D plot
  if(!is.null(sdimx) & !is.null(sdimy) & !is.null(sdimz)) {

    cat('create 3D plot')

    if(!is.null(cell_color) & cell_color %in% colnames(cell_locations_metadata)) {
      if(is.null(cell_color_code)) cell_color_code <- 'lightblue'
      p <- plot_ly(type = 'scatter3d',
                   x = cell_locations_metadata$sdimx, y = cell_locations_metadata$sdimy, z = cell_locations_metadata$sdimz,
                   color = cell_locations_metadata[[cell_color]],
                   mode = 'markers', colors = cell_color_code)
      print(p)
    } else {
      p <- plot_ly(type = 'scatter3d',
                   x = cell_locations_metadata$sdimx, y = cell_locations_metadata$sdimy, z = cell_locations_metadata$sdimz,
                   mode = 'markers', colors = 'lightblue')
      print(p)
    }



  } else {

    pl <- ggplot()
    pl <- pl + theme_classic()

    if(!is.null(spatial_network) & show_network == TRUE) {
      if(is.null(network_color)) network_color = 'red'
      pl <- pl + geom_segment(data = spatial_network[!unified_int %in% interaction_name], aes(x = sdimx_begin, y = sdimy_begin, xend = sdimx_end, yend = sdimy_end), color = 'lightgrey', size = 0.5, alpha = 0.5)
      pl <- pl + geom_segment(data = spatial_network[unified_int %in% interaction_name], aes(x = sdimx_begin, y = sdimy_begin, xend = sdimx_end, yend = sdimy_end), color = network_color, size = 0.5, alpha = 0.5)
    }

    if(!is.null(spatial_grid) & show_grid == TRUE) {
      if(is.null(grid_color)) grid_color = 'black'
      pl <- pl + geom_rect(data = spatial_grid, aes(xmin = x_start, xmax = x_end, ymin = y_start, ymax = y_end), color = grid_color, fill = NA)
    }

    # cell color default
    if(is.null(cell_color)) {
      cell_color = 'lightblue'
      pl <- pl + geom_point(data = cell_locations[!cell_ID %in% cell_IDs_to_keep], aes_string(x = sdimx, y = sdimy), show.legend = show.legend, shape = 21, fill = 'lightgrey', size = point_size_other)
      pl <- pl + geom_point(data = cell_locations[cell_ID %In% cell_IDs_to_keep], aes_string(x = sdimx, y = sdimy), show.legend = show.legend, shape = 21, fill = cell_color, size = point_size_select)
    }
    else if (is.character(cell_color)) {
      if(cell_color %in% colnames(cell_locations_metadata)) {

        if(color_as_factor == TRUE) {
          factor_data = factor(cell_locations_metadata[[cell_color]])
          cell_locations_metadata[[cell_color]] <- factor_data
        }

        pl <- pl + geom_point(data = cell_locations_metadata[!cell_ID %in% cell_IDs_to_keep], aes_string(x = sdimx, y = sdimy), fill = 'lightgrey', shape = 21, size = point_size_other)
        pl <- pl + geom_point(data = cell_locations_metadata[cell_ID %in% cell_IDs_to_keep], aes_string(x = sdimx, y = sdimy, fill = cell_color), show.legend = show.legend, shape = 21, size = point_size_select)



        if(!is.null(cell_color_code)) {
          pl <- pl + scale_fill_manual(values = cell_color_code)
        } else if(color_as_factor == T) {
          number_colors = length(unique(factor_data))
          cell_color_code = Giotto:::getDistinctColors(n = number_colors)
          names(cell_color_code) = unique(factor_data)
          pl <- pl + scale_fill_manual(values = cell_color_code)
        } else if(color_as_factor == F){
          pl <- pl + scale_fill_gradient(low = 'blue', high = 'red')
        }

      } else {
        pl <- pl + geom_point(data = cell_locations_metadata[!cell_ID %in% cell_IDs_to_keep], aes_string(x = sdimx, y = sdimy), show.legend = show.legend, shape = 21, fill = 'lightgrey', size = point_size_other)
        pl <- pl + geom_point(data = cell_locations_metadata[cell_ID %in% cell_IDs_to_keep], aes_string(x = sdimx, y = sdimy), show.legend = show.legend, shape = 21, fill = cell_color, size = point_size_select)
      }

    }

    pl <- pl + theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                                  legend.title = element_text(size = 10),
                                  legend.text = element_text(size = 10))

    # fix coord ratio
    if(!is.null(coord_fix_ratio)) {
      pl <- pl + coord_fixed(ratio = coord_fix_ratio)
    }

    pl <- pl + labs(x = 'x coordinates', y = 'y coordinates')
    pl

    print(pl)
  }

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
  setorder(GTGscore, interaction, -diff_spat)


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

  pl <- ggplot()
  pl <- pl + theme_classic()
  pl <- pl + geom_segment(data = topDT, aes(x = comb_expr, xend = all_comb_expr, y = gene_gene, yend = gene_gene), linetype = 2)
  pl <- pl + geom_point(data = topDT, aes_string(x = 'all_comb_expr', y = 'gene_gene'), color = 'blue', size = 2)
  pl <- pl + geom_point(data = topDT, aes_string(x = 'comb_expr', y = 'gene_gene'), color = 'red', size = 2)
  pl <- pl + facet_wrap( ~ interaction, ncol = length(unique(topDT$interaction)))
  pl

}


