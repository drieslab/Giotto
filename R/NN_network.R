


#' @title createNearestNetwork
#' @name createNearestNetwork
#' @description create a nearest neighbour (NN) network
#' @param gobject giotto object
#' @param type sNN or kNN
#' @param dim_reduction_to_use dimension reduction method to use
#' @param dim_reduction_name name of dimension reduction set to use
#' @param dimensions_to_use number of dimensions to use as input
#' @param name arbitrary name for NN network
#' @param genes_to_use if dim_reduction_to_use = NULL, which genes to use
#' @param expression_values expression values to use
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @param k number of k neighbors to use
#' @param minimum_shared minimum shared neighbors
#' @param top_shared keep at ...
#' @param verbose be verbose
#' @param ... additional parameters for kNN and sNN functions from dbscan
#' @return giotto object with updated NN network
#' @details This function creates a k-nearest neighbour (kNN) or shared nearest neighbour (sNN) network
#' based on the provided dimension reduction space. To run it directly on the gene expression matrix
#' set \emph{dim_reduction_to_use = NULL}.
#'
#' See also \code{\link[dbscan]{kNN}} and \code{\link[dbscan]{sNN}} for more information about
#' how the networks are created.
#'
#' Output for kNN:
#' \itemize{
#'   \item{from: }{cell_ID for source cell}
#'   \item{to: }{cell_ID for target cell}
#'   \item{distance: }{distance between cells}
#'   \item{weight: }{weight = 1/(1 + distance)}
#' }
#'
#' Output for sNN:
#' \itemize{
#'   \item{from: }{cell_ID for source cell}
#'   \item{to: }{cell_ID for target cell}
#'   \item{distance: }{distance between cells}
#'   \item{weight: }{1/(1 + distance)}
#'   \item{shared: }{number of shared neighbours}
#'   \item{rank: }{ranking of pairwise cell neighbours}
#' }
#' For sNN networks two additional parameters can be set:
#' \itemize{
#'   \item{minimum_shared: }{minimum number of shared neighbours needed}
#'   \item{top_shared: }{keep this number of the top shared neighbours, irrespective of minimum_shared setting}
#' }
#'
#' @export
#' @examples
#'
#' data(mini_giotto_single_cell)
#'
#' mini_giotto_single_cell <- createNearestNetwork(gobject = mini_giotto_single_cell,
#'                                                 dimensions_to_use = 1:3, k = 3)
#'
createNearestNetwork <- function(gobject,
                                 type = c('sNN', 'kNN'),
                                 dim_reduction_to_use = 'pca',
                                 dim_reduction_name = 'pca',
                                 dimensions_to_use = 1:10,
                                 genes_to_use = NULL,
                                 expression_values = c('normalized', 'scaled', 'custom'),
                                 name = 'sNN.pca',
                                 return_gobject = TRUE,
                                 k = 30,
                                 minimum_shared = 5,
                                 top_shared = 3,
                                 verbose = T,
                                 ...) {

  # type of NN network
  type = match.arg(type, c('sNN', 'kNN'))

  ## using dimension reduction ##
  if(!is.null(dim_reduction_to_use)) {

    ## check if reduction exists
    check = gobject@dimension_reduction[['cells']][[dim_reduction_to_use]][[dim_reduction_name]]
    if(is.null(check)) stop('dimension reduction does not exist, check if you did ', dim_reduction_to_use,
                            ' and if ', dim_reduction_name, ' was the name used')

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

    # features as columns & cells as rows
    matrix_to_use = t_giotto(expr_values)

  }

  # vector for cell_ID
  cell_names = rownames(matrix_to_use)
  names(cell_names) = 1:nrow(matrix_to_use)

  ## run nearest-neighbour algorithm ##
  if(k >= nrow(matrix_to_use)) {
    k = (nrow(matrix_to_use)-1)
    if(verbose == TRUE) cat('\n k is higher than total number of cells, adjusted to (total number of cells - 1) \n')
  }

  nn_network = dbscan::kNN(x = matrix_to_use, k = k, sort = TRUE, ...)

  # data.table variables
  from = to = weight = distance = from_cell_ID = to_cell_ID = shared = NULL

  nn_network_dt = data.table::data.table(from = rep(1:nrow(nn_network$id), k),
                                         to = as.vector(nn_network$id),
                                         weight = 1/(1 + as.vector(nn_network$dist)),
                                         distance = as.vector(nn_network$dist))
  nn_network_dt[, from_cell_ID := cell_names[from]]
  nn_network_dt[, to_cell_ID := cell_names[to]]


  if(type == 'sNN') {

    snn_network = dbscan::sNN(x = nn_network, k = k, kt = NULL, ...)
    snn_network_dt = data.table::data.table(from = rep(1:nrow(snn_network$id), k),
                                            to = as.vector(snn_network$id),
                                            weight = 1/(1 + as.vector(snn_network$dist)),
                                            distance = as.vector(snn_network$dist),
                                            shared = as.vector(snn_network$shared))
    snn_network_dt = snn_network_dt[stats::complete.cases(snn_network_dt)]
    snn_network_dt[, from_cell_ID := cell_names[from]]
    snn_network_dt[, to_cell_ID := cell_names[to]]

    # rank snn
    data.table::setorder(snn_network_dt, from, -shared)
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
#' @description Add a network layout for a selected nearest neighbor network
#' @param gobject giotto object
#' @param nn_network_to_use kNN or sNN
#' @param network_name name of NN network to be used
#' @param layout_type layout algorithm to use
#' @param options_list list of options for selected layout
#' @param layout_name name for layout
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @return giotto object with updated layout for selected NN network
#' @details This function creates layout coordinates based on the provided kNN or sNN.
#' Currently only the force-directed graph layout "drl", see \code{\link[igraph]{layout_with_drl}},
#' is implemented. This provides an alternative to tSNE or UMAP based visualizations.
#' @export
addNetworkLayout = function(gobject,
                            nn_network_to_use = "sNN",
                            network_name = "sNN.pca",
                            layout_type = c('drl'),
                            options_list = NULL,
                            layout_name = 'layout',
                            return_gobject = TRUE) {

  ## checks
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
#' @keywords internal
#' @return kNN object
nnDT_to_kNN <- function(nnDT) {

  # data.table variable
  from = NULL

  k = unique(table(nnDT$from))

  if(length(k) > 1) {
    stop('\n k is not the same for all cells \n')
  }

  nnDT[, rank := 1:.N, by = from]

  # distance matrix
  dist_prep = data.table::dcast.data.table(nnDT, formula = from~rank, value.var = 'distance')
  dist_prep[, from := NULL]
  dist_matrix = as.matrix(dist_prep)

  # id matrix
  id_prep = data.table::dcast.data.table(nnDT, formula = from~rank, value.var = 'to')
  id_prep[, from := NULL]
  id_matrix = as.matrix(id_prep)

  return(structure(list(dist = dist_matrix, id = id_matrix, k = k, sort = TRUE),
                   class = c("kNN", "NN")))
}



#' @title select_NearestNetwork
#' @name select_NearestNetwork
#' @description Extracts a NN-network from a Giotto object
#' @param gobject giotto object
#' @param nn_network_to_use kNN or sNN
#' @param network_name name of NN network to be used
#' @param output return a igraph or data.table object
#' @return igraph or data.table object
#' @keywords internal
select_NearestNetwork = function(gobject,
                                 nn_network_to_use = 'sNN',
                                 network_name = 'sNN.pca',
                                 output = c('igraph', 'data.table')) {

  output = match.arg(arg = output, choices = c('igraph', 'data.table'))

  ## select network to use
  if(is.null(nn_network_to_use) | is.null(network_name)) {
    stop('\n you need to select network type: knn or snn \n
         and you need to select the network name you created\n')
  } else {
    igraph_object = gobject@nn_network[[nn_network_to_use]][[network_name]][['igraph']]
    if(is.null(igraph_object)) {
      cat('\n nn_network_to_use or network_name does not exist, \n
           create a nearest-neighbor network first \n')
    }
  }

  ## convert igraph to data.table
  if(output == 'data.table') {
    igraph_object = data.table::as.data.table(igraph::get.data.frame(x = igraph_object))
    return(igraph_object)
  }

  return(igraph_object)

}


#' @title extractNearestNetwork
#' @name extractNearestNetwork
#' @description Function to calculate gene signature enrichment scores per spatial position using a rank based approach.
#' @inheritDotParams select_NearestNetwork
#' @seealso \code{\link{select_NearestNetwork}}
#' @keywords internal
extractNearestNetwork = function(...) {

  .Deprecated(new = "select_NearestNetwork")

  select_NearestNetwork(...)

}




