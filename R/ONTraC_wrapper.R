#' @title getONTraCv1Input
#' @name getONTraCv1Input
#' @description generate the input data for ONTraC v1
#' @inheritParams data_access_params
#' @inheritParams read_data_params
#' @param output_path the path to save the output file
#' @param cell_type the cell type column name in the metadata
#' @returns data.table with columns: Cell_ID, Sample, x, y, Cell_Type
#' @details This function generate the input data for ONTraC v1
#' @examples
#' g <- GiottoData::loadGiottoMini("visium")
#'
#' getONTraCv1Input(
#'   gobject = g,
#'   cell_type = "custom_leiden"
#' )
#' @export
getONTraCv1Input <- function(gobject, # nolint: object_name_linter.
                             cell_type,
                             output_path = getwd(),
                             spat_unit = NULL,
                             feat_type = NULL,
                             verbose = TRUE) {
  # Set feat_type and spat_unit
  spat_unit <- set_default_spat_unit(
    gobject = gobject,
    spat_unit = spat_unit
  )
  feat_type <- set_default_feat_type(
    gobject = gobject,
    spat_unit = spat_unit,
    feat_type = feat_type
  )

  pos_df <- getSpatialLocations(
    gobject = gobject,
    spat_unit = spat_unit,
    output = "data.table"
  )
  meta_df <- pDataDT(
    gobject = gobject,
    spat_unit = spat_unit,
    feat_type = feat_type
  )
  output_df <- merge(x = pos_df, y = meta_df, by = "cell_ID")

  # check if the cell_type column exits
  if (!cell_type %in% colnames(output_df)) {
    vmsg(.v = verbose, paste(
      "Given",
      cell_type,
      "do not exist in giotto object's metadata!"
    ))
    return(NULL)
  }

  # add default sample name for one sample obj
  if (!"list_ID" %in% colnames(output_df)) {
    output_df$list_ID <- "ONTraC"
  }

  output_df <- output_df[, .SD, .SDcols = c(
    "cell_ID",
    "list_ID",
    "sdimx",
    "sdimy",
    cell_type
  )]
  colnames(output_df) <- c("Cell_ID", "Sample", "x", "y", "Cell_Type")
  file_path <- file.path(output_path, "ONTraC_dataset_input.csv")
  write.csv(output_df, file = file_path, quote = FALSE, row.names = FALSE)
  vmsg(.v = verbose, paste("ONTraC input file was saved as", file_path))

  return(output_df)
}


#' @title load_cell_bin_niche_cluster
#' @name load_cell_bin_niche_cluster
#' @description load cell-level binarized niche cluster
#' @inheritParams data_access_params
#' @inheritParams read_data_params
#' @param ontrac_results_dir the directory where the ONTraC results are saved
#' @returns gobject with cell-level binarized niche cluster
#' @details This function loads the ONTraC outputed cell-level binarized niche
#' cluster into the giotto object.
load_cell_bin_niche_cluster <- function(gobject,
                                        ontrac_results_dir = getwd()) {
  bin_niche_cluster_df <- read.csv(file = file.path(
    ontrac_results_dir,
    "GNN_dir", "cell_level_max_niche_cluster.csv.gz"
  ))
  colnames(bin_niche_cluster_df) <- c("cell_ID", "NicheCluster")
  gobject <- GiottoClass::addCellMetadata(gobject,
    new_metadata = bin_niche_cluster_df,
    by_column = TRUE,
    column_cell_ID = "cell_ID"
  )
  return(gobject)
}


#' @title load_cell_NT_score
#' @name load_cell_NT_score
#' @description load cell-level NT score
#' @inheritParams data_access_params
#' @inheritParams read_data_params
#' @param ontrac_results_dir the directory where the ONTraC results are saved
#' @returns gobject with cell-level NT score
#' @details This function loads the ONTraC outputed cell-level NT score
load_cell_NT_score <- function(gobject, # nolint: object_name_linter.
                               ontrac_results_dir = getwd()) {
  NT_score_df <- read.csv(file = file.path( # nolint: object_name_linter.
    ontrac_results_dir,
    "NTScore_dir", "NTScore.csv.gz"
  ))[c("Cell_ID", "Cell_NTScore")]
  colnames(NT_score_df) <- c("cell_ID", "NTScore") # nolint: object_name_linter.
  gobject <- addCellMetadata(gobject,
    new_metadata = NT_score_df, # nolint: object_name_linter.
    by_column = TRUE,
    column_cell_ID = "cell_ID"
  )

  return(gobject)
}


#' @title load_cell_niche_cluster_prob
#' @name load_cell_niche_cluster_prob
#' @description load cell-niche cluster probability
#' @inheritParams data_access_params
#' @inheritParams read_data_params
#' @param ontrac_results_dir the directory where the ONTraC results are saved
#' @param name name for the probability matrix
#' @returns gobject with cell-niche cluster probability matrix
#' @details This function loads the ONTraC outputed cell-niche cluster
#' probability as an exprObj into the giotto object.
load_cell_niche_cluster_prob <- function(gobject,
                                         ontrac_results_dir = getwd(),
                                         spat_unit = "cell",
                                         feat_type = "niche cluster",
                                         name = "prob") {
  niche_cluster_prob_df <- read.csv(file = file.path(
    ontrac_results_dir,
    "GNN_dir", "cell_level_niche_cluster.csv.gz"
  ))
  rownames(niche_cluster_prob_df) <- niche_cluster_prob_df$Cell_ID
  niche_cluster_prob_df$Cell_ID <- NULL
  expobj <- createExprObj(t(niche_cluster_prob_df),
    spat_unit = spat_unit,
    feat_type = feat_type
  )
  gobject <- GiottoClass::setExpression(
    gobject = gobject,
    x = expobj,
    spat_unit = spat_unit,
    feat_type = feat_type,
    name = name
  )

  return(gobject)
}


#' @title load_nc_connectivity
#' @name load_nc_connectivity
#' @description load niche cluster connectivity
#' @inheritParams data_access_params
#' @inheritParams read_data_params
#' @param ontrac_results_dir the directory where the ONTraC results are saved
#' @param name name for the connectivity matrix
#' @returns gobject with niche cluster connectivity matrix
#' @details This function loads the ONTraC outputed niche cluster connectivity
#' matrix as an exprObj into the giotto object.
load_nc_connectivity <- function(gobject,
                                 ontrac_results_dir = getwd(),
                                 spat_unit = "niche cluster",
                                 feat_type = "connectivity",
                                 name = "normalized") {
  connectivity_df <- read.csv(file = file.path(
    ontrac_results_dir,
    "GNN_dir", "consolidate_out_adj.csv.gz"
  ), header = FALSE)
  rownames(connectivity_df) <- paste0(
    "NicheCluster_",
    seq_len(dim(connectivity_df)[1]) - 1
  )
  colnames(connectivity_df) <- paste0(
    "NicheCluster_",
    seq_len(dim(connectivity_df)[2]) - 1
  )
  expobj <- createExprObj(t(connectivity_df),
    spat_unit = spat_unit,
    feat_type = feat_type
  )
  gobject <- GiottoClass::setExpression(
    gobject = gobject,
    x = expobj,
    spat_unit = spat_unit,
    feat_type = feat_type,
    name = name
  )

  return(gobject)
}


#' @title loadOntraCResults
#' @name loadOntraCResults
#' @description load ONTraC results
#' @inheritParams data_access_params
#' @inheritParams read_data_params
#' @param ontrac_results_dir the directory where the ONTraC results are saved
#' @returns gobject with ONTraC results
#' @details This function loads the ONTraC results into the giotto object.
#' @export
loadOntraCResults <- function(gobject, # nolint: object_name_linter.
                              ontrac_results_dir = getwd()) {
  gobject <- load_cell_bin_niche_cluster(gobject, ontrac_results_dir)
  gobject <- load_cell_NT_score(gobject, ontrac_results_dir)
  gobject <- load_cell_niche_cluster_prob(gobject, ontrac_results_dir)
  gobject <- GiottoClass::addCellMetadata(
    gobject = gobject,
    spat_unit = "cell",
    feat_type = "niche cluster",
    new_metadata = pDataDT(gobject),
    by_column = TRUE,
    column_cell_ID = "cell_ID"
  )

  gobject <- load_nc_connectivity(gobject, ontrac_results_dir)

  return(gobject)
}

#' @title plotNicheClusterConnectivity
#' @name plotNicheClusterConnectivity
#' @description plot niche cluster connectivity
#' @inheritParams data_access_params
#' @inheritParams plot_output_params
#' @param spat_unit name of spatial unit niche stored cluster features
#' @param feat_type name of the feature type stored niche cluster connectivities
#' @param values name of the expression matrix stored connectivity values
#' @details This function plots the niche cluster connectivity matrix
#' @export
plotNicheClusterConnectivity <- function( # nolint: object_name_linter.
    gobject,
    spat_unit = "niche cluster",
    feat_type = "connectivity",
    values = "normalized",
    show_plot = NULL,
    return_plot = NULL,
    save_plot = NULL,
    save_param = list(),
    theme_param = list(),
    default_save_name = "NicheClusterConnectivity") {
  # load `guide_edge_colourbar` function in ggraph,
  # otherwise it will raise an error when using `scale_edge_colour_gradientn`
  library(ggraph)

  # get the niche cluster connectivity matrix
  niche_cluster_connectivites <- getExpression(
    gobject = gobject,
    values = "normalized",
    spat_unit = "niche cluster",
    feat_type = "connectivity",
    output = "matrix"
  )

  # transform the matrix to data.frame for constructing igraph object
  niche_cluster_connectivites <- cbind(
    expand.grid(dimnames(niche_cluster_connectivites)),
    value = as.vector(as.matrix(
      niche_cluster_connectivites
    ))
  )
  colnames(niche_cluster_connectivites) <- c("from", "to", "connectivites")

  # construct igraph object
  igd <- igraph::graph_from_data_frame(
    d = niche_cluster_connectivites[, c("from", "to", "connectivites")],
    directed = FALSE
  )
  igd <- igraph::simplify(
    graph = igd,
    remove.loops = TRUE,
    remove.multiple = FALSE
  )
  edges_sizes <- igraph::edge_attr(igd, "connectivites")
  edges_colors <- edges_sizes
  igd <- igraph::set.edge.attribute(
    graph = igd,
    index = igraph::E(igd),
    name = "color",
    value = edges_colors
  )
  igd <- igraph::set.edge.attribute(
    graph = igd,
    index = igraph::E(igd),
    name = "size",
    value = edges_sizes
  )

  # plot
  ## layout
  coords <- igraph::layout_with_drl(
    graph = igd,
    weights = edges_sizes,
    use.seed = TRUE
  )
  gpl <- ggraph::ggraph(graph = igd, layout = coords)

  ## edges
  gpl <- gpl + ggraph::geom_edge_link(
    ggplot2::aes(
      colour = edges_sizes,
      edge_width = 5,
      edge_alpha = size # nolint: object_usage_linter.
    ),
    show.legend = FALSE
  )
  gpl <- gpl + ggraph::scale_edge_alpha(range = c(0.1, 1))
  gpl <- gpl + ggraph::scale_edge_colour_gradientn(
    colours = RColorBrewer::brewer.pal(9, "Reds"),
    name = "Value"
  )

  ## node
  gpl <- gpl + ggraph::geom_node_point(
    ggplot2::aes(colour = name), # nolint: object_usage_linter.
    size = 10
  )
  gpl <- gpl + ggplot2::scale_fill_gradientn(colours = viridis::turbo(100))
  gpl <- gpl + ggraph::geom_node_text(
    ggplot2::aes(label = name), # nolint: object_usage_linter.
    repel = TRUE
  )

  ## theme
  gpl <- gpl + ggplot2::theme_bw() + ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    panel.border = ggplot2::element_blank(),
    axis.title = ggplot2::element_blank(),
    axis.text = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank()
  )
  gpl

  # return or save
  return(GiottoVisuals::plot_output_handler(
    gobject = gobject,
    plot_object = gpl,
    save_plot = save_plot,
    return_plot = return_plot,
    show_plot = show_plot,
    default_save_name = default_save_name,
    save_param = save_param,
    else_return = NULL
  ))
}

#' @title plotCTCompositionInNicheCluster
#' @name plotCTCompositionInNicheCluster
#' @description plot cell type composition within each niche cluster
#' @param cell_type the cell type column name in the metadata
#' @inheritParams data_access_params
#' @inheritParams plot_output_params
#' @param spat_unit name of spatial unit niche stored cluster features
#' @param feat_type name of the feature type stored niche cluster connectivities
#' @param values name of the expression matrix stored connectivity values
#' @details This function plots the niche cluster connectivity matrix
#' @export
plotCTCompositionInNicheCluster <- function( # nolint: object_name_linter.
    gobject,
    cell_type,
    values = "prob",
    spat_unit = "cell",
    feat_type = "niche cluster",
    show_plot = NULL,
    return_plot = NULL,
    save_plot = NULL,
    save_param = list(),
    theme_param = list(),
    default_save_name = "CellTypeCompositionInNicheCluster") {
  # Get the cell type composition within each niche cluster
  ## extract the cell-level niche cluster probability matrix
  exp <- getExpression(
    gobject = gobject,
    values = values,
    spat_unit = spat_unit,
    feat_type = feat_type,
    output = "exprObj"
  )
  prob_df <- as.data.frame(t(as.matrix(exp@exprMat)))
  prob_df$cell_ID <- rownames(prob_df)
  ## combine the cell type and niche cluster probability matrix
  combined_df <- merge(
    pDataDT(gobject, feat_type = feat_type)[, c("cell_ID", cell_type)],
    prob_df,
    by = "cell_ID"
  )

  # Calculate the normalized cell type composition within each niche cluster
  cell_type_counts_df <- combined_df %>%
    tidyr::pivot_longer(
      cols = dplyr::starts_with("NicheCluster_"),
      names_to = "Cluster",
      values_to = "Probability"
    ) %>%
    dplyr::group_by(Cell_Type, Cluster) %>% # nolint: object_usage_linter.
    dplyr::summarise(Sum = sum(Probability, # nolint: object_usage_linter.
      na.rm = TRUE
    )) %>%
    tidyr::spread(key = "Cluster", value = "Sum", fill = 0)
  cell_type_counts_df <- as.data.frame(cell_type_counts_df)
  rownames(cell_type_counts_df) <- cell_type_counts_df$Cell_Type
  cell_type_counts_df$Cell_Type <- NULL
  normalized_df <- as.data.frame(t(
    t(cell_type_counts_df) / colSums(cell_type_counts_df)
  ))


  # Reshape the data frame into long format
  normalized_df$Cell_Type <- rownames(normalized_df)
  df_long <- normalized_df %>%
    tidyr::pivot_longer(
      cols = -Cell_Type, # nolint: object_usage_linter.
      names_to = "Cluster",
      values_to = "Composition"
    )

  # Create the heatmap using ggplot2
  pl <- ggplot(df_long, aes(
    x = Cell_Type, # nolint: object_usage_linter.
    y = Cluster, # nolint: object_usage_linter.
    fill = Composition # nolint: object_usage_linter.
  )) +
    geom_tile() +
    viridis::scale_fill_viridis(option = "inferno", limits = c(0, 1)) +
    theme_minimal() +
    labs(
      title = "Normalized cell type compositions within each niche cluster",
      x = "Cell_Type",
      y = "Cluster"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  # return or save
  return(GiottoVisuals::plot_output_handler(
    gobject = gobject,
    plot_object = pl,
    save_plot = save_plot,
    return_plot = return_plot,
    show_plot = show_plot,
    default_save_name = default_save_name,
    save_param = save_param,
    else_return = NULL
  ))
}


#' @title plotCellTypeNTScore
#' @name plotCellTypeNTScore
#' @description plot NTScore by cell type
#' @param cell_type the cell type column name in the metadata
#' @inheritParams data_access_params
#' @inheritParams plot_output_params
#' @export
plotCellTypeNTScore <- function(gobject, # nolint: object_name_linter.
                                cell_type,
                                values = "NTScore",
                                spat_unit = "cell",
                                feat_type = "rna",
                                show_plot = NULL,
                                return_plot = NULL,
                                save_plot = NULL,
                                save_param = list(),
                                theme_param = list(),
                                default_save_name = "CellTypeNTScore") {
  # Get the cell type composition within each niche cluster
  data_df <- pDataDT(
    gobject = gobject,
    spat_unit = spat_unit,
    feat_type = feat_type
  )
  avg_scores <- data_df %>%
    dplyr::group_by(Cell_Type) %>% # nolint: object_usage_linter.
    dplyr::summarise(Avg_NTScore = mean(NTScore)) # nolint: object_usage_linter.
  data_df$Cell_Type <- factor(data_df$Cell_Type,
    levels = avg_scores$Cell_Type[order(avg_scores$Avg_NTScore)]
  )

  pl <- ggplot(data_df, aes(
    x = NTScore, # nolint: object_usage_linter.
    y = Cell_Type, # nolint: object_usage_linter.
    fill = Cell_Type
  )) +
    geom_violin() +
    theme_minimal() +
    labs(
      title = "Violin Plot of NTScore by Cell Type",
      x = "NTScore",
      y = "Cell Type"
    ) +
    ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1))

  # return or save
  return(GiottoVisuals::plot_output_handler(
    gobject = gobject,
    plot_object = pl,
    save_plot = save_plot,
    return_plot = return_plot,
    show_plot = show_plot,
    default_save_name = default_save_name,
    save_param = save_param,
    else_return = NULL
  ))
}
