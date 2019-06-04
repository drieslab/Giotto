


#' @title createNearestNetwork
#' @name createNearestNetwork
#' @description create a nearest neighbour network based on previously computed \cr
#' dimension reductions
#' @param gobject giotto object
#' @param expression_values expression values to use
#' @param type kNN or sNN
#' @param dim_reduction_to_use dimension reduction method to use
#' @param dim_reduction_name name of dimension reduction set to use
#' @param dimensions_to_use number of dimensions to use as input
#' @param name arbitrary name for NN network
#' @param genes_to_use if dim_reduction_to_use = NULL, which genes to use
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @param k number of k neighbors to use
#' @param minimum_shared minimum shared neighbors
#' @param top_shared keep at ...
#' @param ... additional parameters
#' @return giotto object with updated NN network
#' @details Description of nearest neighbor network creation and filter steps.
#' @export
#' @examples
#'     createNearestNetwork(gobject)
createNearestNetwork <- function(gobject,
                                 expression_values = c('normalized', 'scaled', 'custom'),
                                 type = c('sNN', 'kNN'),
                                 dim_reduction_to_use = 'pca',
                                 dim_reduction_name = 'pca',
                                 dimensions_to_use = 1:10,
                                 genes_to_use = NULL,
                                 name = 'sNN.pca',
                                 return_gobject = TRUE,
                                 k = 30,
                                 minimum_shared = 5,
                                 top_shared = 3,
                                 ...) {

  # type of NN network
  type = match.arg(type, c('sNN', 'kNN'))

  ## using dimension reduction ##
  if(!is.null(dim_reduction_to_use)) {

    ## TODO: check if reduction exists

    # use only available dimensions if dimensions < dimensions_to_use
    dim_coord = gobject@dimension_reduction[['cells']][[dim_reduction_to_use]][[dim_reduction_name]][['coordinates']]
    dimensions_to_use = dimensions_to_use[dimensions_to_use %in% 1:ncol(dim_coord)]
    matrix_to_use = dim_coord[, dimensions_to_use]


  } else {
    ## using original matrix ##
    values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))
    expr_values = select_expression_values(gobject = gobject, values = values)

    # subset expression matrix
    if(!is.null(genes_to_use)) {
      expr_values = expr_values[rownames(expr_values) %in% genes_to_use, ]
    }

    # features as columns
    # cells as rows
    matrix_to_use = t(expr_values)

  }

  # vector for cell_ID
  cell_names = rownames(matrix_to_use)
  names(cell_names) = 1:nrow(matrix_to_use)


  ## run nearest-neighbour algorithm ##
  nn_network = dbscan::kNN(x = matrix_to_use, k = k, sort = TRUE, ...)
  nn_network_dt = data.table(from = rep(1:nrow(nn_network$id), k),
                             to = as.vector(nn_network$id),
                             weight = 1/(1 + as.vector(nn_network$dist)),
                             distance = as.vector(nn_network$dist))
  nn_network_dt[, from_cell_ID := cell_names[from]]
  nn_network_dt[, to_cell_ID := cell_names[to]]


  if(type == 'sNN') {

    snn_network = dbscan::sNN(x = nn_network, k = k, kt = NULL, ...)
    snn_network_dt = data.table(from = rep(1:nrow(snn_network$id), k),
                                to = as.vector(snn_network$id),
                                weight = 1/(1 + as.vector(snn_network$dist)),
                                distance = as.vector(snn_network$dist),
                                shared = as.vector(snn_network$shared))
    snn_network_dt = snn_network_dt[complete.cases(snn_network_dt)]
    snn_network_dt[, from_cell_ID := cell_names[from]]
    snn_network_dt[, to_cell_ID := cell_names[to]]

    # rank snn
    setorder(snn_network_dt, from, -shared)
    snn_network_dt[, rank := 1:.N, by = from]

    # filter snn
    snn_network_dt = snn_network_dt[rank <= top_shared | shared >= minimum_shared]

  }

  ## convert to igraph object
  all_index = unique(x = c(nn_network_dt$from_cell_ID, nn_network_dt$to_cell_ID))


  if(type == 'kNN') {
    nn_network_igraph = igraph::graph_from_data_frame(nn_network_dt[,.(from_cell_ID, to_cell_ID, weight, distance)], directed = TRUE, vertices = all_index)
  } else if(type == 'sNN') {
    missing_indices = all_index[!all_index %in% unique(snn_network_dt$from)]
    nn_network_igraph = igraph::graph_from_data_frame(snn_network_dt[,.(from_cell_ID, to_cell_ID, weight, distance, shared, rank)], directed = TRUE, vertices = all_index)
  }




  if(return_gobject == TRUE) {

    nn_names = names(gobject@nn_network[[type]])

    if(name %in% nn_names) {
      cat('\n ', name, ' has already been used, will be overwritten \n')

    }

    gobject@nn_network[[type]][[name]][['igraph']] <- nn_network_igraph

    ## update parameters used ##
    parameters_list = gobject@parameters
    number_of_rounds = length(parameters_list)
    update_name = paste0(number_of_rounds,'_nn_network')
    # parameters to include
    parameters_list[[update_name]] = c('nearest-neighbour type' = type,
                                       'dimension reduction used' = dim_reduction_to_use,
                                       'name for dimension reduction' = dim_reduction_name,
                                       'dimensions used' = length(dimensions_to_use),
                                       'expression values' = expression_values,
                                       'number of genes used' = length(genes_to_use),
                                       'k neighbours' = k,
                                       'minimum shared edges for edge in snn' = minimum_shared,
                                       'top edges for edge in snn' = top_shared,
                                       'name for network' = name)
    gobject@parameters = parameters_list

    return(gobject)

  } else {
    return(nn_network_igraph)
  }

}


#' @title addNetworkLayout
#' @name addNetworkLayout
#' @description Add a network layout for a select nearest neighbor network
#' @param gobject giotto object
#' @param nn_network_to_use kNN or sNN
#' @param network_name name of NN network to be used
#' @param layout_type layout algorithm to use
#' @param options_list list of options for selected layout
#' @param layout_name name for layout
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @return giotto object with updated layout for selected NN network
#' @details Description of layouts and options.
#' @export
#' @examples
#'     addNetworkLayout(gobject)
addNetworkLayout = function(gobject,
                            nn_network_to_use = NULL,
                            network_name = NULL,
                            layout_type = c('drl'),
                            options_list = NULL,
                            layout_name = 'layout',
                            return_gobject = TRUE) {

  if(is.null(nn_network_to_use) | is.null(network_name)) {
    stop('\n first create a nearest network \n')
  }
  ig_object = gobject@nn_network[[nn_network_to_use]][[network_name]][['igraph']]

  layout_type = match.arg(arg = layout_type, c('drl'))

  if(layout_type == 'drl') {
    if(is.null(options_list)) {
      layout_options = igraph::drl_defaults$default
    } else {
      layout_options = options_list
    }
    layout_coord = igraph::layout_with_drl(graph = ig_object, options = layout_options)
  }


  if(return_gobject == TRUE) {

    nn_names = names(gobject@nn_network[[nn_network_to_use]])

    if(layout_name %in% nn_names) {
      cat('\n ', layout_name, ' has already been used, will be overwritten \n')

    }

    gobject@nn_network[[nn_network_to_use]][[network_name]][['layout']] <- layout_coord

    ## update parameters used ##
    parameters_list = gobject@parameters
    number_of_rounds = length(parameters_list)
    update_name = paste0(number_of_rounds,'_nn_network_layout')
    # parameters to include
    parameters_list[[update_name]] = c('nn network to use' = nn_network_to_use,
                                       'network name' = network_name,
                                       'layout type' = layout_type,
                                       'layout name ' = layout_name)
    gobject@parameters = parameters_list

    return(gobject)


  } else {
    return(layout_coord)
  }



}


#' @title nnDT_to_kNN
#' @name nnDT_to_kNN
#' @description Convert a nearest network data.table to a kNN object
#' @param nnDT nearest neighbor network in data.table format
#' @return kNN object
nnDT_to_kNN <- function(nnDT) {

  k = unique(table(nnDT$from))

  if(length(k) > 1) {
    stop('\n k is not the same for all cells \n')
  }

  nnDT[, rank := 1:.N, by = from]

  # distance matrix
  dist_prep = dcast.data.table(nnDT, formula = from~rank, value.var = 'distance')
  dist_prep[, from := NULL]
  dist_matrix = as.matrix(dist_prep)

  # id matrix
  id_prep = dcast.data.table(nnDT, formula = from~rank, value.var = 'to')
  id_prep[, from := NULL]
  id_matrix = as.matrix(id_prep)

  return(structure(list(dist = dist_matrix, id = id_matrix, k = k, sort = TRUE),
                   class = c("kNN", "NN")))
}

