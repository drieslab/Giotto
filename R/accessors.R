


## ACCESSOR functions to access and retrieve values in one of the giotto class slots ##


#' @name  select_expression_values
#' @description helper function to select expression values
#' @param gobject giotto object
#' @param feat_type feature type
#' @param values expression values to extract
#' @return expression matrix
#' @keywords internal
select_expression_values <- function(gobject,
                                     feat_type = 'rna',
                                     values) {

  potential_values = names(gobject@expression[[feat_type]])

  ## special cases for giotto standard pipeline
  if(values == 'scaled' & is.null(gobject@expression[[feat_type]][[values]])) {
    stop('run first scaling (& normalization) step')
  } else if(values == 'normalized' & is.null(gobject@expression[[feat_type]][[values]])) {
    stop('run first normalization step')
  } else if(values == 'custom' & is.null(gobject@expression[[feat_type]][[values]])) {
    stop('first add custom expression matrix')
  }


  if(values %in% potential_values) {
    expr_values = gobject@expression[[feat_type]][[values]]
    return(expr_values)
  } else {
    stop("The ", feat_type ," expression matrix with name ","'", values, "'"," can not be found \n")
  }

}



#' @name select_spatial_locations
#' @description helper function to select spatial location data.table
#' @param gobject giotto object
#' @param spat_loc_name name of spatial locations
#' @return data.table with coordinates
#' @keywords internal
select_spatial_locations <- function(gobject,
                                     spat_loc_name = 'raw') {


  # spatial locations
  # if NULL (not given) and spatial locations have been added, then use first one
  # if NULL (not given) and spatial loactions have NOT been added, then keep NULL
  if(is.null(spat_loc_name)) {
    if(!is.null(gobject@spatial_locs)) {
      spat_loc_name = names(gobject@spatial_locs)[[1]]
      cat('No spatial locations have been selected, the first one -',spat_loc_name, '- will be used \n')
    } else {
      spat_loc_name = NULL
      cat('No spatial locations have been found \n')
      return(NULL)
    }
  }

  potential_names = names(gobject@spatial_locs)

  if(spat_loc_name %in% potential_names) {
    spatloc = data.table::copy(gobject@spatial_locs[[spat_loc_name]])
    return(spatloc)
  } else {
    stop("The spatial locations with name ","'", spat_loc_name, "'"," can not be found \n")
  }
}




#' @name select_dimReduction
#' @description helper function to select a dimension reduction object
#' @keywords internal
#' @return dim reduction coordinates (default) or dim reduction object
select_dimReduction = function(gobject,
                               reduction = c('cells', 'genes'),
                               reduction_method = c('pca', 'umap', 'tsne'),
                               name = 'pca',
                               return_dimObj = FALSE) {


  ## check parameters
  reduction = match.arg(arg = reduction, choices = c('cells', 'genes'))
  reduction_method = match.arg(arg = reduction_method, choices = unique(c('pca', 'umap', 'tsne', reduction_method)))


  ## check reduction
  reduction_res = gobject@dimension_reduction[[reduction]]
  if(is.null(reduction_res)) {
    stop('No dimension reduction for ', reduction, ' has been applied \n')
  }

  ## check method
  reduction_res = reduction_res[[reduction_method]]
  if(is.null(reduction_res)) {
    stop(reduction_method, ' has not been performed on this dataset \n')
  }

  ## check name for method
  reduction_res = reduction_res[[name]]
  if(is.null(reduction_res)) {
    stop(name, ': this name is not available for method: ', reduction_method, '\n')
  }

  ## return object or coordinates
  if(return_dimObj == TRUE) {
    return(reduction_res)
  } else {
    return(reduction_res$coordinates)
  }

}



#' @name select_NearestNetwork
#' @description selects a NN-network from a Giotto object
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

#' @name extractNearestNetwork
#' @description Function to calculate gene signature enrichment scores per spatial position using a rank based approach.
#' @inheritDotParams select_NearestNetwork
#' @seealso \code{\link{select_NearestNetwork}}
#' @keywords internal
extractNearestNetwork = function(...) {

  .Deprecated(new = "select_NearestNetwork")

  select_NearestNetwork(...)

}



#' @name select_spatialNetwork
#' @description function to select a spatial network
#' @param gobject giotto object
#' @param name name of spatial network
#' @param return_network_Obj return network object (default = FALSE)
#' @keywords internal
select_spatialNetwork <- function(gobject,
                                  name = NULL,
                                  return_network_Obj = FALSE) {

  if (!is.element(name, names(gobject@spatial_network))){
    message = sprintf("spatial network %s has not been created. Returning NULL.
                      check which spatial networks exist with showNetworks() \n", name)
    warning(message)
    return(NULL)
  }else{
    networkObj = gobject@spatial_network[[name]]
    networkDT = networkObj$networkDT
  }

  if (return_network_Obj == TRUE){
    return(networkObj)
  }else{
    return(networkDT)
  }
}


#' @name select_spatialGrid
#' @description accessor function to select spatial grid
#' @param gobject giotto object
#' @param name name of spatial grid
#' @param return_network_Obj return grid object (default = FALSE)
#' @keywords internal
select_spatialGrid <- function(gobject,
                               name = NULL,
                               return_grid_Obj = FALSE) {

  if (!is.element(name, names(gobject@spatial_grid))){
    message = sprintf("spatial grid %s has not been created. Returning NULL.
                      check which spatial grids exist with showGrids() \n", name)
    warning(message)
    return(NULL)
  }else{
    gridObj = gobject@spatial_grid[[name]]
    gridDT = gridObj$gridDT
  }

  if (return_grid_Obj == TRUE){
    return(gridObj)
  }else{
    return(gridDT)
  }
}




#' @name select_polygon_info
#' @description  selection giotto polygon spatVector
#' @param gobject giotto object
#' @param polygon_name name of polygons
#' @param polygon_overlap include polygon overlap information
#' @keywords internal
select_polygon_info = function(gobject,
                               polygon_name = 'cell',
                               polygon_overlap = NULL) {

  potential_names = names(gobject@spatial_info)

  if(!polygon_name %in% potential_names) {
    stop('There is no polygon information with name ', polygon_name, '\n')
  } else {

    if(is.null(polygon_overlap)) {
      poly_info = gobject@spatial_info[[polygon_name]]@spatVector
    } else {
      potential_overlaps = names(gobject@spatial_info[[polygon_name]]@overlaps)

      if(!polygon_overlap %in% potential_overlaps) {
        stop('There is no polygon overlap information with name ', polygon_overlap, '\n')
      } else {
        poly_info = gobject@spatial_info[[polygon_name]]@overlaps[[polygon_overlap]]
      }
    }
    return(poly_info)
  }
}

#' @name select_feature_info
#' @description  selection giotto points spatVector
#' @param gobject giotto object
#' @param feat_name name of feature
#' @keywords internal
select_feature_info = function(gobject,
                               feat_name = 'rna') {

  potential_names = names(gobject@feat_info)

  if(!feat_name %in% potential_names) {
    stop('There is no feature information with name ', feat_name, '\n')
  } else {
    feat_info = gobject@feat_info[[feat_name]]@spatVector
    return(feat_info)
  }
}



