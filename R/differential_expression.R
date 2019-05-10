



#' @title findMarkers
#' @name findMarkers
#' @description Identify marker genes for selected clusters
#' @param gobject giotto object
#' @param expression_values gene expression values to use
#' @param cluster_column clusters to use
#' @param subset_clusters selection of clusters to compare
#' @param group_1 group 1 cluster IDs from cluster_column for pairwise comparison
#' @param group_2 group 2 cluster IDs from cluster_column for pairwise comparison
#' @param method method to use
#' @param ... additional parameters
#' @return differentially expressed genes (DEG) output
#' @details Description of different DEG methods.
#' @export
#' @examples
#'     findMarkers(gobject)
findMarkers <- function(gobject,
                        expression_values = c('normalized', 'scaled', 'custom'),
                        cluster_column,
                        subset_clusters = NULL,
                        group_1 = NULL,
                        group_2 = NULL,
                        method = c('scran'),
                        ...) {


  # expression data
  values = match.arg(expression_values, choices = c('normalized', 'scaled', 'custom'))
  expr_data = Giotto:::select_expression_values(gobject = gobject, values = values)

  # select method
  method = match.arg(method, choices = c('scran'))

  # cluster column
  cell_metadata = pDataDT(gobject)
  if(!cluster_column %in% colnames(cell_metadata)) {
    stop('\n cluster column not found \n')
  }



  # subset expression data
  if(!is.null(subset_clusters)) {

    cell_metadata = cell_metadata[get(cluster_column) %in% subset_clusters]

    subset_cell_IDs = cell_metadata[['cell_ID']]
    expr_data = expr_data[, colnames(expr_data) %in% subset_cell_IDs]


  } else if(!is.null(group_1) & !is.null(group_2)) {

    cell_metadata = cell_metadata[get(cluster_column) %in% c(group_1, group_2)]

    # create new pairwise group
    group_1_name = paste0(group_1, collapse = '_')
    group_2_name = paste0(group_2, collapse = '_')
    cell_metadata[, pairwise_select_comp := ifelse(get(cluster_column) %in% group_1, group_1_name, group_2_name)]

    cluster_column = 'pairwise_select_comp'

    # expression data
    subset_cell_IDs = cell_metadata[['cell_ID']]
    expr_data = expr_data[, colnames(expr_data) %in% subset_cell_IDs]


  }



  if(method == 'scran') {

    marker_results = scran::findMarkers(x = expr_data, clusters = cell_metadata[[cluster_column]], ...)

    savelist = lapply(names(marker_results), FUN = function(x) {

      dfr = marker_results[[x]]

      DT = as.data.table(dfr)
      DT[, gene_ID := rownames(dfr)]
      DT[, cluster_ID := x]

    })

    return(savelist)

  } else if(method != 'scran') {

    cat('\n no other DEG methods implemented yet \n')

  }



}



#' @title findGiniMarkers
#' @name findGiniMarkers
#' @description Identify marker genes for selected clusters
#' @param gobject giotto object
#' @param expression_values gene expression values to use
#' @param cluster_column clusters to use
#' @param subset_clusters selection of clusters to compare
#' @param group_1 group 1 cluster IDs from cluster_column for pairwise comparison
#' @param group_2 group 2 cluster IDs from cluster_column for pairwise comparison
#' @param min_expr_gini_score minimum gini coefficient on expression
#' @param min_det_gini_score minimum gini coefficient on detection
#' @param detection_threshold detection threshold for gene expression
#' @param rank_score rank scores to include
#' @return gini genes output
#' @details Description of parameters.
#' @export
#' @examples
#'     findGiniMarkers(gobject)
findGiniMarkers <- function(gobject,
                            expression_values = c('normalized', 'scaled', 'custom'),
                            cluster_column,
                            subset_clusters = NULL,
                            group_1 = NULL,
                            group_2 = NULL,
                            min_expr_gini_score = 0.5,
                            min_det_gini_score = 0.5,
                            detection_threshold = 0,
                            rank_score = 1) {

  values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))


  # cluster column
  cell_metadata = pDataDT(gobject)
  if(!cluster_column %in% colnames(cell_metadata)) {
    stop('\n cluster column not found \n')
  }


  # subset
  if(!is.null(subset_clusters)) {

    cell_metadata = cell_metadata[get(cluster_column) %in% subset_clusters]
    subset_cell_IDs = cell_metadata[['cell_ID']]
    gobject = subsetGiotto(gobject = gobject, cell_ids = subset_cell_IDs)

  } else if(!is.null(group_1) & !is.null(group_2)) {

    cell_metadata = cell_metadata[get(cluster_column) %in% c(group_1, group_2)]

    # create new pairwise group
    group_1_name = paste0(group_1, collapse = '_')
    group_2_name = paste0(group_2, collapse = '_')
    cell_metadata[, pairwise_select_comp := ifelse(get(cluster_column) %in% group_1, group_1_name, group_2_name)]

    cluster_column = 'pairwise_select_comp'

    # expression data
    subset_cell_IDs = cell_metadata[['cell_ID']]
    gobject = subsetGiotto(gobject = gobject, cell_ids = subset_cell_IDs)

  }




  # average expression per cluster
  aggr_sc_clusters <- create_average_DT(gobject = gobject, meta_data_name = cluster_column,
                                        expression_values = values)
  aggr_sc_clusters_DT <- as.data.table(aggr_sc_clusters)
  aggr_sc_clusters_DT[, genes := rownames(aggr_sc_clusters)]
  aggr_sc_clusters_DT_melt <- melt.data.table(aggr_sc_clusters_DT, variable.name = 'cluster', id.vars = 'genes', value.name = 'expression')

  # detection per cluster
  aggr_detection_sc_clusters <- create_average_detection_DT(gobject = gobject, meta_data_name = cluster_column,
                                                            expression_values = values, detection_threshold = detection_threshold)
  aggr_detection_sc_clusters_DT <- as.data.table(aggr_detection_sc_clusters)
  aggr_detection_sc_clusters_DT[, genes := rownames(aggr_detection_sc_clusters)]
  aggr_detection_sc_clusters_DT_melt <- melt.data.table(aggr_detection_sc_clusters_DT, variable.name = 'cluster', id.vars = 'genes', value.name = 'detection')

  # gini
  aggr_sc_clusters_DT_melt[, expression_gini := mygini_fun(expression), by = genes]
  aggr_detection_sc_clusters_DT_melt[, detection_gini := mygini_fun(detection), by = genes]

  # combine
  aggr_sc = cbind(aggr_sc_clusters_DT_melt, aggr_detection_sc_clusters_DT_melt[,.(detection, detection_gini)])
  aggr_sc[, expression_rank := rank(-expression), by = genes]
  aggr_sc[, detection_rank := rank(-detection), by = genes]

  top_genes_scores = aggr_sc[expression_rank <= rank_score & detection_rank <= rank_score]
  top_genes_scores_filtered = top_genes_scores[expression > min_expr_gini_score & detection > min_det_gini_score]
  top_genes_scores_filtered[, comb_gini := expression_gini*detection_gini]
  top_genes_scores_filtered[, comb_rank := rank(-comb_gini), by = cluster]
  setorder(top_genes_scores_filtered, cluster, comb_rank)

  return(top_genes_scores_filtered)

}
