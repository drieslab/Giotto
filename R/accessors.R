


## Get and set functions to get and set values in one of the giotto class slots ##


## expression values slot ####

#' @title  Get expression values
#' @name  get_expression_values
#' @description Function to get expression values from giotto object
#' @param gobject giotto object
#' @param spat_unit spatial unit (e.g. "cell")
#' @param feat_type feature type (e.g. "rna", "dna", "protein")
#' @param values expression values to extract (e.g. "raw", "normalized", "scaled")
#' @return expression matrix
#' @family expression accessor functions
#' @family functions to get data from giotto object
#' @export
get_expression_values <- function(gobject,
                                  feat_type = NULL,
                                  spat_unit = NULL,
                                  values) {


  # Set feat_type and spat_unit
  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = spat_unit)
  feat_type = set_default_feat_type(gobject = gobject,
                                    spat_unit = spat_unit,
                                    feat_type = feat_type)

  potential_values = names(gobject@expression[[spat_unit]][[feat_type]])

  ## special cases for giotto standard pipeline
  if(values == 'scaled' & is.null(gobject@expression[[spat_unit]][[feat_type]][[values]])) {
    stop('run first scaling (& normalization) step')
  } else if(values == 'normalized' & is.null(gobject@expression[[spat_unit]][[feat_type]][[values]])) {
    stop('run first normalization step')
  } else if(values == 'custom' & is.null(gobject@expression[[spat_unit]][[feat_type]][[values]])) {
    stop('first add custom expression matrix')
  }


  if(values %in% potential_values) {
    expr_values = gobject@expression[[spat_unit]][[feat_type]][[values]]
    return(expr_values)
  } else {
    stop("The spatial unit ", spat_unit ," for expression matrix ", feat_type, " and with name ","'", values, "'"," can not be found \n")
  }

}


#' @title select_expression_values
#' @name select_expression_values
#' @inheritDotParams get_expression_values
#' @seealso \code{\link{get_expression_values}}
#' @keywords internal
select_expression_values = function(...) {

  .Deprecated(new = "get_expression_values")

  get_expression_values(...)

}


#' @title  Set expression values
#' @name  set_expression_values
#' @description Function to set expression values for giotto object
#' @param gobject giotto object
#' @param spat_unit spatial unit  (e.g. "cell")
#' @param feat_type feature type (e.g. "rna", "dna", "protein")
#' @param name name for the expression slot
#' @param values matrix of expression values
#' @return giotto object
#' @family expression accessor functions
#' @family functions to set data in giotto object
#' @export
set_expression_values <- function(gobject,
                                  spat_unit = NULL,
                                  feat_type = NULL,
                                  name = 'test',
                                  values) {


  # Set feat_type and spat_unit
  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = spat_unit)
  feat_type = set_default_feat_type(gobject = gobject,
                                    spat_unit = spat_unit,
                                    feat_type = feat_type)


  ## 1. check if specified name has already been used
  potential_names = names(gobject@expression[[spat_unit]][[feat_type]])
  if(name %in% potential_names) {
    cat(name, ' already exist and will be replaced with new values \n')
  }

  ## TODO: 2. check input for values matrix


  ## 3. update and return giotto object
  gobject@expression[[spat_unit]][[feat_type]][[name]] = values
  return(gobject)

}





## spatial locations slot ####

#' @title Get spatial locations
#' @name get_spatial_locations
#' @description Function to get a spatial location data.table
#' @param gobject giotto object
#' @param spat_unit spatial unit (e.g. "cell")
#' @param spat_loc_name name of spatial locations (defaults to first name in spatial_locs slot, e.g. "raw")
#' @return data.table with coordinates
#' @family spatial location data accessor functions
#' @family functions to get data from giotto object
#' @export
get_spatial_locations <- function(gobject,
                                  spat_unit = NULL,
                                  spat_loc_name = NULL) {


  # spatial locations
  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = spat_unit)

  # if NULL (not given) and spatial locations have been added, then use first one
  # if NULL (not given) and spatial locations have NOT been added, then keep NULL
  if(is.null(spat_loc_name)) {
    if(!is.null(gobject@spatial_locs)) {
      spat_loc_name = names(gobject@spatial_locs[[spat_unit]])[[1]]
      # cat('No spatial locations have been selected, the first one -',spat_loc_name, '- will be used \n')
    } else {
      spat_loc_name = NULL
      cat('No spatial locations have been found \n')
      return(NULL)
    }
  }

  potential_names = names(gobject@spatial_locs[[spat_unit]])

  if(spat_loc_name %in% potential_names) {
    spatloc = data.table::copy(gobject@spatial_locs[[spat_unit]][[spat_loc_name]])
    return(spatloc)
  } else {
    stop("The spatial locations with name ","'", spat_loc_name, "'"," can not be found \n")
  }
}



#' @title select_spatial_locations
#' @name select_spatial_locations
#' @inheritDotParams get_spatial_locations
#' @seealso \code{\link{get_spatial_locations}}
#' @keywords internal
select_spatial_locations = function(...) {

  .Deprecated(new = "get_spatial_locations")

  get_spatial_locations(...)

}




#' @title Set spatial locations
#' @name set_spatial_locations
#' @description Function to set a spatial location slot
#' @param gobject giotto object
#' @param spat_unit spatial unit (e.g. "cell")
#' @param spat_loc_name name of spatial locations, default "raw"
#' @param spatlocs spatial locations
#' @param verbose be verbose
#' @return giotto object
#' @family spatial location data accessor functions
#' @family functions to set data in giotto object
#' @export
set_spatial_locations <- function(gobject,
                                  spat_unit = NULL,
                                  spat_loc_name = 'raw',
                                  spatlocs,
                                  verbose = TRUE) {

  # spatial locations
  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = spat_unit)

  ## 1. check if specified name has already been used
  potential_names = names(gobject@spatial_locs[[spat_unit]][[spat_loc_name]])
  if(spat_loc_name %in% potential_names) {
    if(verbose == TRUE) {
    cat(spat_loc_name, ' already exists and will be replaced with new spatial locations \n')
    }
  }

  ## TODO: 2. check input for spatial locations


  ## 3. update and return giotto object
  gobject@spatial_locs[[spat_unit]][[spat_loc_name]] = spatlocs
  return(gobject)

}






## dimension reduction slot ####

#' @title Get dimension reduction
#' @name get_dimReduction
#' @param gobject giotto object
#' @param spat_unit spatial unit (e.g. "cell")
#' @param feat_type feature type (e.g. "rna", "dna", "protein")
#' @param reduction reduction on cells or features (e.g. "cells", "feats")
#' @param reduction_method reduction method (e.g. "pca", "umap", "tsne")
#' @param name name of reduction results
#' @param return_dimObj return full dimension object result. Default = FALSE
#' @description Function to get a dimension reduction object
#' @return dim reduction coordinates (default) or dim reduction object
#' @family dimensional reduction data accessor functions
#' @family functions to get data from giotto object
#' @export
get_dimReduction = function(gobject,
                            spat_unit = NULL,
                            feat_type = NULL,
                            reduction = c('cells', 'feats'),
                            reduction_method = c('pca', 'umap', 'tsne'),
                            name = 'pca',
                            return_dimObj = FALSE) {


  # Set feat_type and spat_unit
  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = spat_unit)
  feat_type = set_default_feat_type(gobject = gobject,
                                    spat_unit = spat_unit,
                                    feat_type = feat_type)

  ## check parameters
  reduction = match.arg(arg = reduction, choices = c('cells', 'feats'))
  reduction_method = match.arg(arg = reduction_method, choices = unique(c('pca', 'umap', 'tsne', reduction_method)))

  ## check reduction
  reduction_res = gobject@dimension_reduction[[reduction]][[spat_unit]][[feat_type]]
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


#' @title select_dimReduction
#' @name select_dimReduction
#' @inheritDotParams get_dimReduction
#' @seealso \code{\link{get_dimReduction}}
#' @keywords internal
select_dimReduction = function(...) {

  .Deprecated(new = "get_dimReduction")

  get_dimReduction(...)

}


#' @title Set dimension reduction
#' @name set_dimReduction
#' @description Function to set a dimension reduction slot
#' @param gobject giotto object
#' @param spat_unit spatial unit (e.g. "cell")
#' @param feat_type feature type (e.g. "rna", "dna", "protein")
#' @param reduction reduction on cells or features
#' @param reduction_method reduction method (e.g. "pca")
#' @param name name of reduction results
#' @param dimObject dimension object result to set
#' @return giotto object
#' @family dimensional reduction data accessor functions
#' @family functions to set data in giotto object                                    dimObject = pca_object)
#' @export
set_dimReduction <- function(gobject,
                             spat_unit = NULL,
                             feat_type = NULL,
                             reduction = c('cells', 'genes'),
                             reduction_method = c('pca', 'umap', 'tsne'),
                             name = 'pca',
                             dimObject) {


  # Set feat_type and spat_unit
  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = spat_unit)
  feat_type = set_default_feat_type(gobject = gobject,
                                    spat_unit = spat_unit,
                                    feat_type = feat_type)

  ## 1. check if specified name has already been used
  potential_names = names(gobject@dimension_reduction[[reduction]][[spat_unit]][[feat_type]][[reduction_method]])
  if(name %in% potential_names) {
    cat(name, ' already exist and will be replaced with new dimension reduction object \n')
  }

  ## TODO: 2. check input for dimension reduction object


  ## 3. update and return giotto object
  gobject@dimension_reduction[[reduction]][[spat_unit]][[feat_type]][[reduction_method]][[name]] = dimObject

  return(gobject)

}







## nearest neighbor network slot ####

#' @title Get nearest network
#' @name get_NearestNetwork
#' @description Get a NN-network from a Giotto object
#' @param gobject giotto object
#' @param spat_unit spatial unit (e.g. "cell")
#' @param feat_type feature type (e.g. "rna", "dna", "protein")
#' @param nn_network_to_use "kNN" or "sNN"
#' @param network_name name of NN network to be used
#' @param output return a igraph or data.table object. Default 'igraph'
#' @return igraph or data.table object
#' @family expression space nearest network accessor functions
#' @family functions to get data from giotto object
#' @export
get_NearestNetwork = function(gobject,
                              spat_unit = NULL,
                              feat_type = NULL,
                              nn_network_to_use = NULL,
                              network_name = NULL,
                              output = c('igraph', 'data.table')) {

  output = match.arg(arg = output, choices = c('igraph', 'data.table'))

  # Set feat_type and spat_unit
  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = spat_unit)
  feat_type = set_default_feat_type(gobject = gobject,
                                    spat_unit = spat_unit,
                                    feat_type = feat_type)


  # automatic nearest network selection
  if(is.null(nn_network_to_use)) {
    nn_network_to_use = names(gobject@nn_network[[spat_unit]])[[1]]
    if(is.null(nn_network_to_use)) {
      stop('There is currently no nearest-neighbor network build for spatial unit: ', spat_unit)
    } else {
      cat('The NN network type was not specified, default to the first: ', nn_network_to_use)
    }
  }

  if(is.null(network_name)) {
    network_name = names(gobject@nn_network[[spat_unit]][[nn_network_to_use]])[[1]]
    if(is.null(network_name)) {
      stop('There is currently no nearest-neighbor network build for spatial unit: ', spat_unit,
           'and network type: ', nn_network_to_use)
    }else {
      cat('The NN network name was not specified, default to the first: ', network_name)
    }
  }

  # get igraph object
  igraph_object = gobject@nn_network[[spat_unit]][[nn_network_to_use]][[network_name]]
  if(is.null(igraph_object)) {
    stop('nn_network_to_use: ', nn_network_to_use, ' or network_name: ', network_name, 'does not exist.
          Create a nearest-neighbor network first')
  }


  ## convert igraph to data.table
  if(output == 'data.table') {
    igraph_object = data.table::as.data.table(igraph::get.data.frame(x = igraph_object))
    return(igraph_object)
  }

  return(igraph_object)

}


#' @title Extract nearest network
#' @name extractNearestNetwork
#' @inheritDotParams get_NearestNetwork
#' @seealso \code{\link{get_NearestNetwork}}
#' @keywords internal
extractNearestNetwork = function(...) {

  .Deprecated(new = "get_NearestNetwork")

  get_NearestNetwork(...)

}


#' @title Select nearest network
#' @name select_NearestNetwork
#' @inheritDotParams get_NearestNetwork
#' @seealso \code{\link{get_NearestNetwork}}
#' @keywords internal
select_NearestNetwork = function(...) {

  .Deprecated(new = "get_NearestNetwork")

  get_NearestNetwork(...)

}


#' @title Set nearest network
#' @name set_NearestNetwork
#' @description Set a NN-network for a Giotto object
#' @param gobject giotto object
#' @param spat_unit spatial unit (e.g. "cell")
#' @param feat_type feature type (e.g. "rna", "dna", "protein")
#' @param nn_network_to_use "kNN" or "sNN"
#' @param network_name name of NN network to be used
#' @param nn_network nearest network
#' @return giotto object
#' @family expression space nearest network accessor functions
#' @family functions to set data in giotto object
#' @export
set_NearestNetwork = function(gobject,
                              spat_unit = NULL,
                              feat_type = NULL,
                              nn_network_to_use = 'sNN',
                              network_name = 'sNN.pca',
                              nn_network) {

  # Set feat_type and spat_unit
  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = spat_unit)
  feat_type = set_default_feat_type(gobject = gobject,
                                    spat_unit = spat_unit,
                                    feat_type = feat_type)

  ## 1. check if specified name has already been used
  potential_names = names(gobject@nn_network[[spat_unit]][[nn_network_to_use]])
  if(network_name %in% potential_names) {
    cat(network_name, ' already exist and will be replaced with nearest neighbor network \n')
  }

  ## TODO: 2. check input for nearest neighbor network
  # convert to igraph if data.table class


  ## 3. update and return giotto object
  gobject@nn_network[[spat_unit]][[nn_network_to_use]][[network_name]] = nn_network

  return(gobject)

}






## spatial network slot ####

#' @title Get spatial network
#' @name get_spatialNetwork
#' @description Function to get a spatial network
#' @param gobject giotto object
#' @param spat_unit spatial unit (e.g. "cell")
#' @param feat_type feature type (e.g. "rna", "dna", "protein")
#' @param name name of spatial network
#' @param return_network_Obj return network object (default = FALSE)
#' @family spatial network data accessor functions
#' @family functions to get data from giotto object
#' @export
get_spatialNetwork <- function(gobject,
                               spat_unit = NULL,
                               feat_type = NULL,
                               name = NULL,
                               return_network_Obj = FALSE) {



  # Set feat_type and spat_unit
  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = spat_unit)
  feat_type = set_default_feat_type(gobject = gobject,
                                    spat_unit = spat_unit,
                                    feat_type = feat_type)

  if (!is.element(name, names(gobject@spatial_network[[spat_unit]]))){
    message = sprintf("spatial network %s has not been created. Returning NULL.
                      check which spatial networks exist with showNetworks() \n", name)
    warning(message)
    return(NULL)
  }else{
    networkObj = gobject@spatial_network[[spat_unit]][[name]]
    networkDT = networkObj$networkDT
  }

  if (return_network_Obj == TRUE){
    return(networkObj)
  }else{
    return(networkDT)
  }
}


#' @title Select spatial network
#' @name select_spatialNetwork
#' @inheritDotParams get_spatialNetwork
#' @seealso \code{\link{get_spatialNetwork}}
#' @keywords internal
select_spatialNetwork = function(...) {

  .Deprecated(new = "get_spatialNetwork")

  get_spatialNetwork(...)

}


#' @title Set spatial network
#' @name set_spatialNetwork
#' @description Function to set a spatial network
#' @param gobject giotto object
#' @param spat_unit spatial unit  (e.g. "cell")
#' @param feat_type feature type (e.g. "rna", "dna", "protein")
#' @param name name of spatial network
#' @param spatial_network spatial network
#' @return giotto object
#' @family spatial network data accessor functions
#' @family functions to set data in giotto object
#' @export
set_spatialNetwork <- function(gobject,
                               spat_unit = NULL,
                               feat_type = NULL,
                               name = NULL,
                               spatial_network) {

  # Set feat_type and spat_unit
  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = spat_unit)
  feat_type = set_default_feat_type(gobject = gobject,
                                    spat_unit = spat_unit,
                                    feat_type = feat_type)

  ## 1. check if specified name has already been used
  potential_names = names(gobject@spatial_network[[spat_unit]])
  if(name %in% potential_names) {
    cat(name, ' already exist and will be replaced with new spatial network \n')
  }

  ## TODO: 2. check input for spatial network


  ## 3. update and return giotto object
  gobject@spatial_network[[spat_unit]][[name]] = spatial_network
  return(gobject)

}




## spatial grid slot ####

#' @title Get spatial grid
#' @name get_spatialGrid
#' @description Function to get spatial grid
#' @param gobject giotto object
#' @param spat_unit spatial unit (e.g. "cell")
#' @param feat_type feature type (e.g. "rna", "dna", "protein")
#' @param name name of spatial grid
#' @param return_grid_Obj return grid object (default = FALSE)
#' @family spatial grid data accessor functions
#' @family functions to get data from giotto object
#' @export
get_spatialGrid <- function(gobject,
                            spat_unit = NULL,
                            feat_type = NULL,
                            name = NULL,
                            return_grid_Obj = FALSE) {

  # Set feat_type and spat_unit
  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = spat_unit)
  feat_type = set_default_feat_type(gobject = gobject,
                                    spat_unit = spat_unit,
                                    feat_type = feat_type)

  if (!is.element(name, names(gobject@spatial_grid[[spat_unit]]))){
    message = sprintf("spatial grid %s has not been created. Returning NULL.
                      check which spatial grids exist with showGrids() \n", name)
    warning(message)
    return(NULL)
  }else{
    gridObj = gobject@spatial_grid[[spat_unit]][[name]]
    gridDT = gridObj$gridDT
  }

  if (return_grid_Obj == TRUE){
    return(gridObj)
  }else{
    return(gridDT)
  }
}


#' @title Select spatial grid
#' @name select_spatialGrid
#' @inheritDotParams get_spatialGrid
#' @seealso \code{\link{get_spatialGrid}}
#' @keywords internal
select_spatialGrid = function(...) {

  .Deprecated(new = "get_spatialGrid")

  get_spatialGrid(...)

}


#' @title Set spatial grid
#' @name set_spatialGrid
#' @description Function to set a spatial grid
#' @param gobject giotto object
#' @param spat_unit spatial unit (e.g. "cell")
#' @param feat_type feature type (e.g. "rna", "dna", "protein")
#' @param name name of spatial grid
#' @param spatial_grid spatial grid object
#' @return giotto object
#' @family spatial grid data accessor functions
#' @family functions to set data in giotto object
#' @export
set_spatialGrid <- function(gobject,
                            spat_unit = NULL,
                            feat_type = NULL,
                            name = NULL,
                            spatial_grid) {


  # Set feat_type and spat_unit
  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = spat_unit)
  feat_type = set_default_feat_type(gobject = gobject,
                                    spat_unit = spat_unit,
                                    feat_type = feat_type)

  ## 1. check if specified name has already been used
  potential_names = names(gobject@spatial_grid[[spat_unit]])
  if(name %in% potential_names) {
    cat(name, ' already exist and will be replaced with new spatial grid \n')
  }

  ## TODO: 2. check input for spatial grid


  ## 3. update and return giotto object
  gobject@spatial_grid[[spat_unit]][[name]] = spatial_grid

  return(gobject)

}








## polygon cell info ####

#' @title Get polygon info
#' @name get_polygon_info
#' @description Get giotto polygon spatVector
#' @param gobject giotto object
#' @param polygon_name name of polygons. Default "cell"
#' @param polygon_overlap include polygon overlap information
#' @family polygon info data accessor functions
#' @family functions to get data from giotto object
#' @export
get_polygon_info = function(gobject,
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


#' @title Select polygon info
#' @name select_polygon_info
#' @inheritDotParams get_polygon_info
#' @seealso \code{\link{get_polygon_info}}
#' @keywords internal
select_polygon_info = function(...) {

  .Deprecated(new = "get_polygon_info")

  get_polygon_info(...)

}



#' @title Set polygon info
#' @name set_polygon_info
#' @description Set giotto polygon spatVector
#' @param gobject giotto object
#' @param polygon_name name of polygons. Default "cell"
#' @param gpolygon giotto polygon
#' @return giotto object
#' @family polygon info data accessor functions
#' @family functions to set data in giotto object
#' @export
set_polygon_info = function(gobject,
                            polygon_name = 'cell',
                            gpolygon) {



  ## 1. check if specified name has already been used
  potential_names = names(gobject@spatial_info)
  if(polygon_name %in% potential_names) {
    cat(polygon_name, ' already exist and will be replaced with new giotto polygon \n')
  }

  ## TODO: 2. check input for giotto polygon


  ## 3. update and return giotto object
  gobject@spatial_info[[polygon_name]] = gpolygon
  return(gobject)

}






## feature info ####

#' @title Get feature info
#' @name get_feature_info
#' @description Get giotto points spatVector
#' @param gobject giotto object
#' @param feat_type name of feature (e.g. "rna", "dna", "protein")
#' @family feature info data accessor functions
#' @family functions to get data from giotto object
#' @export
get_feature_info = function(gobject,
                            feat_type = NULL) {

  # specify feat_type
  feat_type = set_default_feat_type(gobject = gobject,
                                    feat_type = feat_type)

  potential_names = names(gobject@feat_info)

  if(!feat_type %in% potential_names) {
    stop('There is no feature information with name ', feat_type, '\n')
  } else {
    feat_info = gobject@feat_info[[feat_type]]@spatVector
    return(feat_info)
  }
}

#' @title Select feature info
#' @name select_feature_info
#' @inheritDotParams get_feature_info
#' @seealso \code{\link{get_feature_info}}
#' @keywords internal
select_feature_info = function(...) {

  .Deprecated(new = "get_feature_info")

  get_feature_info(...)

}


#' @title Set feature info
#' @name set_feature_info
#' @description Set giotto polygon spatVector for features
#' @param gobject giotto object
#' @param feat_type name of feat (e.g. "rna", "dna", "protein")
#' @param gpolygon giotto polygon
#' @return giotto object
#' @family feature info data accessor functions
#' @family functions to set data in giotto object
#' @export
set_feature_info = function(gobject,
                            feat_type = NULL,
                            gpolygon) {

  # specify feat_type
  if(is.null(feat_type)) {
    feat_type = gobject@expression_feat[[1]]
  }

  ## 1. check if specified name has already been used
  potential_names = names(gobject@feat_info)
  if(feat_type %in% potential_names) {
    cat(feat_type, ' already exist and will be replaced with new giotto polygon \n')
  }

  ## TODO: 2. check input for giotto polygon


  ## 3. update and return giotto object
  gobject@feat_info[[feat_type]] = gpolygon
  return(gobject)

}



## spatial enrichment slot ####


#' @title Get spatial enrichment
#' @name get_spatial_enrichment
#' @description Function to get a spatial enrichment data.table
#' @param gobject giotto object
#' @param spat_unit spatial unit (e.g."cell")
#' @param feat_type feature type (e.g. "rna", "dna", "protein")
#' @param enrichm_name name of spatial enrichment results. Default "DWLS"
#' @return data.table with fractions
#' @family spatial enrichment data accessor functions
#' @family functions to get data from giotto object
#' @export
get_spatial_enrichment <- function(gobject,
                                   spat_unit = NULL,
                                   feat_type = NULL,
                                   enrichm_name = 'DWLS') {


  # Set feat_type and spat_unit
  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = spat_unit)
  feat_type = set_default_feat_type(gobject = gobject,
                                    spat_unit = spat_unit,
                                    feat_type = feat_type)

  # spatial locations
  # if NULL (not given) and spatial locations have been added, then use first one
  # if NULL (not given) and spatial loactions have NOT been added, then keep NULL
  if(is.null(enrichm_name)) {
    if(!is.null(gobject@spatial_enrichment)) {
      enrichm_name = names(gobject@spatial_enrichment[[spat_unit]][[feat_type]])[[1]]
      # cat('No spatial locations have been selected, the first one -',spat_loc_name, '- will be used \n')
    } else {
      enrichm_name = NULL
      cat('No spatial enrichment results have been found \n')
      return(NULL)
    }
  }



  potential_names = names(gobject@spatial_enrichment[[spat_unit]][[feat_type]])

  if(enrichm_name %in% potential_names) {
    enr_res = data.table::copy(gobject@spatial_enrichment[[spat_unit]][[feat_type]][[enrichm_name]])
    return(enr_res)
  } else {
    stop("The spatial enrichment result with name ","'", enrichm_name, "'"," can not be found \n")
  }
}


#' @title Set spatial enrichment
#' @name set_spatial_enrichment
#' @description Function to set a spatial enrichment slot
#' @param gobject giotto object
#' @param spat_unit spatial unit (e.g. "cell")
#' @param feat_type feature type (e.g. "rna", "dna", "protein")
#' @param enrichm_name name of spatial enrichment results. Default "DWLS"
#' @param spatenrichment spatial enrichment results
#' @return giotto object
#' @family spatial enrichment data accessor functions
#' @family functions to set data in giotto object
#' @export
set_spatial_enrichment <- function(gobject,
                                   spat_unit = NULL,
                                   feat_type = NULL,
                                   enrichm_name = 'enrichment',
                                   spatenrichment) {


  # Set feat_type and spat_unit
  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = spat_unit)
  feat_type = set_default_feat_type(gobject = gobject,
                                    spat_unit = spat_unit,
                                    feat_type = feat_type)

  ## 1. check if specified name has already been used
  potential_names = names(gobject@spatial_enrichment[[spat_unit]][[feat_type]][[enrichm_name]])
  if(enrichm_name %in% potential_names) {
    cat(enrichm_name, ' already exist and will be replaced with new spatial enrichment results \n')
  }

  ## TODO: 2. check input for spatial locations


  ## 3. update and return giotto object
  gobject@spatial_enrichment[[spat_unit]][[feat_type]][[enrichm_name]] = spatenrichment
  return(gobject)

}



## MG image slot ####



#' @title Get \emph{magick}-based giotto \code{image}
#' @name get_giottoImage_MG
#' @description Get a giottoImage from a giotto object
#' @param gobject giotto object
#' @param name name of giottoImage \code{\link{showGiottoImageNames}}
#' @return a giottoImage
#' @keywords internal
get_giottoImage_MG = function(gobject = NULL,
                              name = NULL) {

  if(is.null(gobject)) stop('The giotto object holding the giottoImage needs to be provided \n')
  g_image_names = names(gobject@images)
  if(is.null(g_image_names)) stop('No giottoImages have been found \n')

  if(is.null(name)) {
    name = g_image_names[1]
  }

  if(!name %in% g_image_names) stop(name, ' was not found among the image names, see showGiottoImageNames()')

  g_image = gobject@images[[name]]

  return(g_image)
}



#' @title Set \emph{magick}-based giotto \code{image}
#' @name set_giottoImage_MG
#' @description Set a giottoImage for a giotto object with no additional modifications
#' @param gobject giotto object
#' @param image_object a giottoImage object
#' @param name name to assign giottoImage
#' @param verbose be verbose
#' @return giotto object
#' @keywords internal
set_giottoImage_MG = function(gobject,
                              image_object,
                              name = NULL,
                              verbose = TRUE) {

  # Check params
  if(is.null(gobject)) stop('gobject must be given \n')
  if(is.null(image_object)) stop('image_object to be attached must be given \n')

  # Default to name present in image object name slot
  if(is.null(name)) name = image_object@name

  # Find existing names
  potential_names = list_images_names(gobject = gobject, img_type = 'image')

  if(verbose == TRUE) {
    if(name %in% potential_names) cat(name, ' already exists and will be replaced with new image object \n')
  }

  gobject@images[[name]] = image_object

  return(gobject)

}



## large image slot ####



#' @title Get \emph{terra}-based giotto \code{largeImage}
#' @name get_giottoLargeImage
#' @description Set a giottoLargeImage from a giottoObject
#' @param gobject giotto object
#' @param name name of giottoLargeImage \code{\link{showGiottoImageNames}}
#' @return a giottoLargeImage
#' @keywords internal
get_giottoLargeImage = function(gobject = NULL,
                                name = NULL) {

  if(is.null(gobject)) stop('The giotto object holding the giottoLargeImage needs to be provided \n')
  g_image_names = names(gobject@largeImages)
  if(is.null(g_image_names)) stop('No giottoLargeImages have been found \n')

  if(is.null(name)) {
    name = g_image_names[1]
  }

  if(!name %in% g_image_names) stop(name,' was not found among the largeImage names. See showGiottoImageNames() \n')

  g_imageL = gobject@largeImages[[name]]

  return(g_imageL)
}




#' @title Set \emph{terra}-based giotto \code{largeImage}
#' @name set_giottoLargeImage
#' @description Set a giottoLargeImage for a giotto object with no additional modifications
#' @param gobject giotto object
#' @param largeImage_object a giottoLargeImage object
#' @param name name to assign giottoLargeImage
#' @param verbose be verbose
#' @return giotto object
#' @keywords internal
set_giottoLargeImage = function(gobject,
                                largeImage_object,
                                name = NULL,
                                verbose = TRUE) {

  # Check params
  if(is.null(gobject)) stop('gobject must be given \n')
  if(is.null(largeImage_object)) stop('largeImage_object to be attached must be given \n')

  # Default to name present in image object name slot
  if(is.null(name)) name = largeImage_object@name

  # Find existing names
  potential_names = list_images_names(gobject = gobject, img_type = 'largeImage')

  if(verbose == TRUE) {
    if(name %in% potential_names) cat(name, 'already exists and will be replaced with new image object \n')
  }

  gobject@largeImages[[name]] = largeImage_object

  return(gobject)

}



## all image slots ####



#' @title Get giotto image object
#' @name get_giottoImage
#' @description Get giotto image object from gobject
#' @param gobject giotto object
#' @param image_type type of giotto image object. Either "image" or "largeImage"
#' @param name name of a giotto image object \code{\link{showGiottoImageNames}}
#' @return a giotto image object
#' @family image data accessor functions
#' @family functions to get data from giotto object
#' @export
get_giottoImage = function(gobject = NULL,
                           image_type = c('image','largeImage'),
                           name = NULL) {

  # Check image type
  image_type = match.arg(image_type, choices = c('image','largeImage'))

  # Select get function
  if(image_type == 'image') {
    g_img = get_giottoImage_MG(gobject = gobject,
                               name = name)
  }
  if(image_type == 'largeImage') {
    g_img = get_giottoLargeImage(gobject = gobject,
                                 name = name)
  }
  return(g_img)
}



#' @title Set giotto image object
#' @name set_giottoImage
#' @description Directly attach a giotto image to giotto object
#' @details \emph{\strong{Use with care!}} This function directly attaches giotto image
#'   objects to the gobject without further modifications of spatial positioning values
#'   within the image object that are generally needed in order for them to
#'   plot in the correct location relative to the other modalities of spatial data. \cr
#'   For the more general-purpose method of attaching image objects, see \code{\link{addGiottoImage}}
#' @param gobject giotto object
#' @param image giotto image object to be attached without modification to the
#'   giotto object
#' @param image_type type of giotto image object. Either "image" or "largeImage"
#' @param name name of giotto image object
#' @param verbose be verbose
#' @return giotto object
#' @family image data accessor functions
#' @family functions to set data in giotto object
#' @seealso \code{\link{addGiottoImage}}
#' @export
set_giottoImage = function(gobject = NULL,
                           image = NULL,
                           image_type = NULL,
                           name = NULL,
                           verbose = TRUE) {

  # Check image type
  image_type = match.arg(image_type, choices = c('image','largeImage'))

  # Select set function
  if(image_type == 'image') {
    gobject = set_giottoImage_MG(gobject = gobject,
                                 image_object = image,
                                 name = name,
                                 verbose = verbose)
  }
  if(image_type == 'largeImage') {
    gobject = set_giottoLargeImage(gobject = gobject,
                                   largeImage_object = image,
                                   name = name,
                                   verbose = verbose)
  }
  return(gobject)
}



## Show functions ####

#' @title showGiottoExpression
#' @name showGiottoExpression
#' @description shows the available matrices
#' @param gobject giotto object
#' @param nrows number of rows to print for each matrix
#' @param ncols number of columns to print for each matrix
#' @return prints the name and small subset of available matrices
#' @family functions to show data in giotto object
#' @keywords show
#' @export
showGiottoExpression = function(gobject, nrows = 4, ncols = 4) {

  if(is.null(gobject)) stop('A giotto object needs to be provided \n')

  available_data = list_expression(gobject = gobject)
  if(is.null(available_data)) cat('No expression data available \n')

  for(spatial_unit in unique(available_data$spat_unit)) {

    cat('Spatial unit: ', spatial_unit, ' \n\n')

    for(feature_type in unique(available_data[available_data$spat_unit == spatial_unit,]$feat_type)) {

      cat('--> Feature: ', feature_type, ' \n\n')

      for(mat_i in available_data[available_data$spat_unit == spat_unit & available_data$feat_type == feature_type,]$name) {

        cat('----> Name: ', mat_i, 'matrix: \n')

        print(gobject@expression[[spatial_unit]][[feature_type]][[mat_i]][1:nrows, 1:ncols])
        cat('\n')
      }
    }
  }
}

#' @title showGiottoSpatLocs
#' @name showGiottoSpatLocs
#' @description shows the available spatial locations
#' @param gobject giotto object
#' @param nrows number of rows to print for each spatial location data.table
#' @return prints the name and small subset of available data.table
#' @family functions to show data in giotto object
#' @keywords show
#' @export
showGiottoSpatLocs = function(gobject, nrows = 4) {

  if(is.null(gobject)) stop('A giotto object needs to be provided \n')

  available_data = list_spatial_locations(gobject = gobject)
  if(is.null(available_data)) cat('No spatial locations available \n')

  for(spatial_unit in unique(available_data$spat_unit)) {

    cat('Spatial unit: ', spatial_unit, ' \n\n')

    for(spatlocname in available_data[available_data$spat_unit == spatial_unit,]$name) {
      cat('--> Name: ', spatlocname, ' \n\n')
      print(gobject@spatial_locs[[spatial_unit]][[spatlocname]][1:nrows,])
      cat('\n')
    }
  }

}


#' @title showGiottoSpatEnrichments
#' @name showGiottoSpatEnrichments
#' @description shows the available spatial enrichment results
#' @param gobject giotto object
#' @param nrows number of rows to print for each spatial enrichment data.table
#' @return prints the name and small subset of available data.table
#' @family functions to show data in giotto object
#' @keywords show
#' @export
showGiottoSpatEnrichments = function(gobject,
                                     nrows = 4) {

  if(is.null(gobject)) stop('A giotto object needs to be provided \n')

  available_data = list_spatial_enrichments(gobject = gobject)

  if(is.null(available_data)) cat('No spatial enrichments available \n')

  for(spatial_unit in unique(available_data$spat_unit)) {

    cat('Spatial unit: ', spatial_unit, ' \n\n')

    for(feature_type in available_data[spat_unit == spatial_unit][['feat_type']]) {

      cat('--> Feature type: ', feature_type, ' \n\n')

      for(spatenrichname in available_data[spat_unit == spatial_unit][feat_type == feature_type][['name']]) {

        cat('----> Name ', spatenrichname, ': \n\n')
        print(gobject@spatial_enrichment[[spatial_unit]][[feature_type]][[spatenrichname]][1:nrows,])

      }

    }

  }


}



#' @title showGiottoDimRed
#' @name showGiottoDimRed
#' @description shows the available dimension reductions
#' @param gobject giotto object
#' @param nrows number of coordinates rows to print
#' @param ncols number of coordinates columns to print
#' @return prints the name and small subset of available dimension reduction coordinates
#' @family functions to show data in giotto object
#' @keywords show
#' @export
showGiottoDimRed = function(gobject,
                            nrows = 3,
                            ncols = 2) {

  if(is.null(gobject)) stop('A giotto object needs to be provided \n')

  available_data = list_dim_reductions(gobject)
  if(is.null(available_data)) cat('No dimensional reductions available \n')

  for(data_type_red in unique(available_data$data_type)) {
    data_type_subset = available_data$data_type == data_type_red

    if(data_type_red == 'feats') cat('Dim reduction on features:')
    if(data_type_red == 'cells') cat('Dim reduction on cells:')

    cat('\n',
        '-------------------------',
        '\n\n\n')

    for(spatial_unit in unique(available_data[data_type_subset,]$spat_unit)) {
      spat_unit_subset = available_data$spat_unit == spatial_unit

      cat('Spatial unit ', spatial_unit, ': \n\n')

      for(feature_type in unique(available_data[data_type_subset & spat_unit_subset,]$feat_type)) {
        feat_type_subset = available_data$feat_type == feature_type

        cat('--> Feature type ', feature_type, ': \n\n')

        for(dimRed_type in unique(available_data[data_type_subset & spat_unit_subset & feat_type_subset,]$dim_type)) {
          dim_type_subset = available_data$dim_type == dimRed_type

          cat('----> Dim reduction type, ', dimRed_type, ': \n\n')

          for(dimRed_name in available_data[data_type_subset & spat_unit_subset & feat_type_subset & dim_type_subset,]$name) {

            cat('------> ', dimRed_name, 'coordinates: \n')
            print(gobject@dimension_reduction[[data_type_red]][[spatial_unit]][[feature_type]][[dimRed_type]][[dimRed_name]][['coordinates']][1:nrows, 1:ncols])
            cat('\n')

          }
        }
      }
    }
  }
}





#' @title showGiottoSpatialInfo
#' @name showGiottoSpatialInfo
#' @description show the available giotto spatial polygon information
#' @param gobject giotto object
#' @family functions to show data in giotto object
#' @keywords show
#' @export
showGiottoSpatialInfo = function(gobject) {

  if(is.null(gobject)) stop('A giotto object needs to be provided \n')

  available_data = list_spatial_info(gobject = gobject)
  if(is.null(available_data)) cat('No spatial info available \n')

  for(info in available_data$spat_info) {

    cat("For Spatial info: ", info, "\n\n")
    print(gobject@spatial_info[[info]])
    cat("-----------------------------")
    cat("\n \n")
  }

}


#' @title showGiottoFeatInfo
#' @name showGiottoFeatInfo
#' @description show the available giotto spatial feature information
#' @param gobject giotto object
#' @family functions to show data in giotto object
#' @keywords show
#' @export
showGiottoFeatInfo = function(gobject) {

  if(is.null(gobject)) stop('A giotto object needs to be provided \n')

  available_data = list_feature_info(gobject = gobject)
  if(is.null(available_data)) cat('No feature info available \n')

  for(info in available_data$feat_info) {

    cat("For Feature info: ", info, "\n\n")
    print(gobject@feat_info[[info]])
    cat("-----------------------------")
    cat("\n \n")
  }

}




#' @title showGiottoSpatNetworks
#' @name showGiottoSpatNetworks
#' @description Prints the available spatial networks that are attached to the Giotto object
#' @param gobject a giotto object
#' @param nrows number of rows to print
#' @return prints names and small subset of available spatial network info
#' @family functions to show data in giotto object
#' @keywords show
#' @export
showGiottoSpatNetworks = function(gobject,
                                  nrows = 4) {

  if(is.null(gobject)) stop('A giotto object needs to be provided \n')

  available_data = list_spatial_networks(gobject = gobject)
  if(is.null(available_data)) cat('No spatial networks are available \n')

  for(spatial_unit in unique(available_data$spat_unit)) {

    cat('Spatial unit:', spatial_unit, '\n\n')

    for(network_name in available_data[available_data$spat_unit == spatial_unit,]$name) {

      cat('--> Name:', network_name, '\n\n')

      print(gobject@spatial_network[[spatial_unit]][[network_name]][['networkDT']][1:nrows,])
      cat('\n')

    }
  }
}


#' @title Show networks
#' @name showNetworks
#' @inheritDotParams showGiottoSpatNetworks
#' @seealso \code{\link{showGiottoSpatNetworks}}
#' @export
showNetworks = function(...) {

  .Deprecated(new = "showGiottoSpatNetworks")

  showGiottoSpatNetworks(...)

}


#' @title showGiottoSpatGrids
#' @name showGiottoSpatGrids
#' @description Prints the available spatial grids that are attached to the Giotto object
#' @param gobject giotto object
#' @param nrows number of rows to print
#' @return prints name of available spatial grids
#' @family functions to show data in giotto object
#' @keywords show
#' @export
showGiottoSpatGrids = function(gobject,
                               nrows = 4) {

  if(is.null(gobject)) stop('A giotto object needs to be provided \n')

  available_data = list_spatial_grids(gobject = gobject)
  if(is.null(available_data)) cat('No available spatial grids \n')

  for(spatial_unit in unique(available_data$spat_unit)) {

    cat('Spatial grid:', spatial_unit, '\n\n')

    for(grid_name in available_data[available_data$spat_unit == spatial_unit,]$name) {

      cat('--> Name:', grid_name, '\n\n')

      print(gobject@spatial_grid[[spatial_unit]][[grid_name]][['gridDT']][1:nrows,])
      cat('\n')

    }
  }
}


#' @title Show Spatial Grids
#' @name showGrids
#' @inheritDotParams showGiottoSpatGrids
#' @seealso \code{\link{showGiottoSpatGrids}}
#' @export
showGrids = function(...) {

  .Deprecated(new = "showGiottoSpatGrids")

  showGiottoSpatGrids(...)

}


#' @title showGiottoImageNames
#' @name showGiottoImageNames
#' @description Prints the available giotto images that are attached to the Giotto object
#' @param gobject a giotto object
#' @return prints names of available giotto image objects
#' @family functions to show data in giotto object
#' @keywords show
#' @export
showGiottoImageNames = function(gobject) {

  if(is.null(gobject)) stop('A giotto object needs to be provided \n')

  available_data = list_images(gobject = gobject)
  if(is.null(available_data)) cat('No available images \n')

  for(image_type in unique(available_data$img_type)) {

    cat('Image type:', image_type, '\n\n')

    for(image_name in available_data[available_data$img_type == image_type,]$name) {

      cat('--> Name:', image_name, '\n')

    }
    cat('\n')
  }
}



# List functions ####


#' @title list_giotto_data
#' @name list_giotto_data
#' @description list the available data within specified giotto object slot
#' @param gobject giotto object
#' @param slot giotto object slot of interest (e.g. "expression", "spatial_locs", etc.)
#' @param ... additional params to pass
#' @return names and locations of data within giotto object slot
#' @keywords internal
list_giotto_data = function(gobject = NULL,
                            slot = NULL,
                            ...) {

  if(slot == 'expression') return(list_expression(gobject = gobject,...))
  if(slot == 'cell_metadata') return(list_cell_metadata(gobject = gobject,...))
  if(slot == 'feat_metadata') return(list_feat_metadata(gobject = gobject,...))
  if(slot == 'spatial_locs') return(list_spatial_locations(gobject = gobject,...))
  if(slot == 'spatial_enrichment') return(list_spatial_enrichments(gobject = gobject,...))
  if(slot == 'dimension_reduction') return(list_dim_reductions(gobject = gobject,...))
  if(slot == 'spatial_info') return(list_spatial_info(gobject = gobject))
  if(slot == 'feat_info') return(list_feature_info(gobject = gobject))
  if(slot == 'spatial_network') return(list_spatial_networks(gobject = gobject,...))
  if(slot == 'spatial_grid') return(list_spatial_grids(gobject = gobject,...))
  if(slot == 'images') return(list_images_names(gobject = gobject, img_type = 'image'))
  if(slot == 'largeImages') return(list_images_names(gobject = gobject, img_type = 'largeImage'))

}


#' @title list_expression
#' @name list_expression
#' @description lists the available matrices
#' @param gobject giotto object
#' @param spat_unit spatial unit (e.g. "cell")
#' @param feat_type feature type (e.g. "rna", "dna", "protein")
#' @return names and locations of available matrices as data.table. col order matters.
list_expression = function(gobject,
                           spat_unit = NULL,
                           feat_type = NULL) {

  availableExpr = data.table()
  for(spatial_unit in names(gobject@expression)) {
    for(feature_type in names(gobject@expression[[spatial_unit]])) {
      for(mat_i in names(gobject@expression[[spatial_unit]][[feature_type]])) {
        availableExpr = rbind(availableExpr,
                              list(spat_unit = spatial_unit,
                                   feat_type = feature_type,
                                   name = mat_i))
      }
    }
  }

  # check if a specific category is desired
  if(!is.null(spat_unit)) spat_unit_subset = availableExpr$spat_unit == spat_unit else spat_unit_subset = TRUE
  if(!is.null(feat_type)) feat_type_subset = availableExpr$feat_type == feat_type else feat_type_subset = TRUE

  availableExpr = availableExpr[spat_unit_subset & feat_type_subset,]

  # return data.table (NULL if empty)
  if(nrow(availableExpr) == 0) return(NULL)
  else return(availableExpr)
}



#' @title list_expression_names
#' @name list_expression_names
#' @description lists the available matrices names for a given spatial unit and feature type
#' @param gobject giotto object
#' @param spat_unit spatial unit (e.g. "cell")
#' @param feat_type feature type (e.g. "rna", "dna", "protein")
#' @return vector with names of available matrices
list_expression_names = function(gobject,
                                 spat_unit,
                                 feat_type) {

  expression_names = names(gobject@expression[[spat_unit]][[feat_type]])

  return(expression_names)
}



#' @title list_cell_metadata
#' @name list_cell_metadata
#' @description lists the available cell metadata
#' @param gobject giotto object
#' @param spat_unit spatial unit (e.g. "cell")
#' @param feat_type feature type (e.g. "rna", "dna", "protein")
#' @return names and locations of available cell metadata as data.table
list_cell_metadata = function(gobject,
                              spat_unit = NULL,
                              feat_type = NULL) {

  availableCMet = data.table()
  for(spatial_unit in names(gobject@cell_metadata)) {
    for(feature_type in names(gobject@cell_metadata[[spatial_unit]])) {
      availableCMet = rbind(availableCMet,
                            list(spat_unit = spatial_unit,
                                 feat_type = feature_type))
    }
  }

  # check if a specific category is desired
  if(!is.null(spat_unit)) spat_unit_subset = availableCMet$spat_unit == spat_unit else spat_unit_subset = TRUE

  availableCMet = availableCMet[spat_unit_subset,]

  # return data.table (NULL if empty)
  if(nrow(availableCMet) == 0) return(NULL)
  else return(availableCMet)
}



#' @title list_feat_metadata
#' @name list_feat_metadata
#' @description lists the available feature metadata
#' @param gobject giotto object
#' @param spat_unit spatial unit (e.g. "cell")
#' @param feat_type feature type (e.g. "rna", "dna", "protein")
#' @return names and locations of available feature metadata as data.table
list_feat_metadata = function(gobject,
                              spat_unit = NULL,
                              feat_type = NULL) {

  availableFMet = data.table()
  for(spatial_unit in names(gobject@feat_metadata)) {
    for(feature_type in names(gobject@feat_metadata[[spatial_unit]])) {
      availableFMet = rbind(availableFMet,
                            list(spat_unit = spatial_unit,
                                 feat_type = feature_type))
    }
  }

  # check if a specific category is desired
  if(!is.null(spat_unit)) spat_unit_subset = availableFMet$spat_unit == spat_unit else spat_unit_subset = TRUE

  availableFMet = availableFMet[spat_unit_subset,]

  # return data.table (NULL if empty)
  if(nrow(availableFMet) == 0) return(NULL)
  else return(availableFMet)
}



#' @title list_spatial_locations
#' @name list_spatial_locations
#' @description shows the available spatial locations
#' @param gobject giotto object
#' @param spat_unit spatial unit (e.g. "cell")
#' @return names and locations of available data.table as data.table
list_spatial_locations = function(gobject,
                                  spat_unit = NULL) {

  availableSpatLocs = data.table()
  for(spatial_unit in names(gobject@spatial_locs)) {
    for(spatloc_name in names(gobject@spatial_locs[[spatial_unit]]))
      availableSpatLocs = rbind(availableSpatLocs,
                                list(spat_unit = spatial_unit,
                                     name = spatloc_name))
  }

  # check if a specific category is desired
  if(!is.null(spat_unit)) {
    availableSpatLocs = availableSpatLocs[availableSpatLocs$spat_unit == spat_unit,]
  }

  if(nrow(availableSpatLocs) == 0) return(NULL)
  else return(availableSpatLocs)
}




#' @title list_spatial_locations_names
#' @name list_spatial_locations_names
#' @description lists the available spatial location names for a given spatial unit
#' @param gobject giotto object
#' @param spat_unit spatial unit (e.g. "cell")
#' @return vector with names of available spatial locations
list_spatial_locations_names = function(gobject,
                                        spat_unit) {

  spatlocs_names = names(gobject@spatial_locs[[spat_unit]])

  return(spatlocs_names)
}



#' @title list_spatial_enrichments
#' @name list_spatial_enrichments
#' @description return the available spatial enrichment results
#' @param gobject giotto object
#' @param spat_unit spatial unit (e.g. "cell")
#' @param feat_type feature type (e.g. "rna", "dna", "protein")
#' @return names and locations of available data as data.table
list_spatial_enrichments = function(gobject,
                                    spat_unit = NULL,
                                    feat_type = NULL) {

  availableSpatEnr = data.table()

  for(spatial_unit in names(gobject@spatial_enrichment)) {

    for(feature_type in names(gobject@spatial_enrichment[[spatial_unit]])) {

      for(spatenr_name in names(gobject@spatial_enrichment[[spatial_unit]][[feature_type]])) {

        availableSpatEnr = rbind(availableSpatEnr,
                                 list(spat_unit = spatial_unit,
                                      feat_type = feature_type,
                                      name = spatenr_name))
      }

    }

  }

  # check if a specific category is desired
  if(!is.null(spat_unit)) {
    availableSpatEnr = availableSpatEnr[availableSpatEnr$spat_unit == spat_unit,]
  }

  if(nrow(availableSpatEnr) == 0) return(NULL)
  else return(availableSpatEnr)
}





#' @title list_spatial_enrichments_names
#' @name list_spatial_enrichments_names
#' @description returns the available spatial enrichment names for a given spatial unit
#' @param gobject giotto object
#' @param spat_unit spatial unit (e.g. "cell")
#' @param feat_type feature type (e.g. "rna", "dna", "protein")
#' @return vector of names for available spatial enrichments
list_spatial_enrichments_names = function(gobject,
                                          spat_unit,
                                          feat_type) {

  spatenr_names = names(gobject@spatial_enrichment[[spat_unit]][[feat_type]])

  return(spatenr_names)
}





#' @title list_dim_reductions
#' @name list_dim_reductions
#' @description return the available dimension reductions
#' @param gobject giotto object
#' @param data_type "cells" or "feats" data used in dim reduction
#' @param spat_unit spatial unit (e.g. "cell")
#' @param feat_type feature type (e.g. "rna", "dna", "protein")
#' @param dim_type dimensional reduction method (e.g. "pca", "umap")
#' @return names and locations of dimension reduction as a data.table
list_dim_reductions = function(gobject,
                               data_type = NULL,
                               spat_unit = NULL,
                               feat_type = NULL,
                               dim_type = NULL) {

  availableDimRed = data.table()
  for(dataType in names(gobject@dimension_reduction)) {
    for(spatUnit in names(gobject@dimension_reduction[[dataType]])) {
      for(featType in names(gobject@dimension_reduction[[dataType]][[spatUnit]])) {
        for(dimType in names(gobject@dimension_reduction[[dataType]][[spatUnit]][[featType]])) {
          for(subType in names(gobject@dimension_reduction[[dataType]][[spatUnit]][[featType]][[dimType]])) {
            availableDimRed = rbind(availableDimRed,
                                    list(data_type = dataType,
                                         spat_unit = spatUnit,
                                         feat_type = featType,
                                         dim_type = dimType,
                                         name = subType))
          }
        }
      }
    }
  }

  # check if a specific category is desired
  if(!is.null(data_type)) data_type_subset = availableDimRed$data_type == data_type else data_type_subset = TRUE
  if(!is.null(spat_unit)) spat_unit_subset = availableDimRed$spat_unit == spat_unit else spat_unit_subset = TRUE
  if(!is.null(feat_type)) feat_type_subset = availableDimRed$feat_type == feat_type else feat_type_subset = TRUE
  if(!is.null(dim_type)) dimred_type_subset = availableDimRed$dim_type == dim_type else dimred_type_subset = TRUE

  availableDimRed = availableDimRed[data_type_subset & spat_unit_subset & feat_type_subset & dimred_type_subset,]

  # NULL if there is no data
  if(nrow(availableDimRed) == 0) return(NULL)
  else return(availableDimRed)
}



#' @title list_dim_reductions_names
#' @name list_dim_reductions_names
#' @description return the available dimension reductions object names
#' @param gobject giotto object
#' @param data_type cells or feats dim reduction
#' @param spat_unit spatial unit (e.g. "cell")
#' @param feat_type feature type (e.g. "rna", "dna", "protein")
#' @param dim_type dimensional reduction type (method)
#' @return names pf dimension reduction object
#' @details function that can be used to find which names have been used
list_dim_reductions_names = function(gobject,
                                     data_type = 'cells',
                                     spat_unit,
                                     feat_type,
                                     dim_type) {

  dim_red_object_names = names(gobject@dimension_reduction[[data_type]][[spat_unit]][[feat_type]][[dim_type]])

  return(dim_red_object_names)
}





#' @title list_spatial_info
#' @name list_spatial_info
#' @description return the available giotto spatial polygon information
#' @param gobject giotto object
#' @return names of available spatial polygon information
list_spatial_info = function(gobject) {

  availableSpatInfo = data.table()
  for(info in names(gobject@spatial_info)) {
    availableSpatInfo = rbind(availableSpatInfo,
                              list(spat_info = info))
  }

  if(nrow(availableSpatInfo) == 0) return(NULL)
  else return(availableSpatInfo)
}




#' @title list_spatial_info_names
#' @name list_spatial_info_names
#' @description return the available names for giotto spatial polygon information
#' @param gobject giotto object
#' @return vector with names of available spatial polygon information
list_spatial_info_names = function(gobject) {

  spat_info_names = names(gobject@spatial_info)

  return(spat_info_names)
}



#' @title list_feature_info
#' @name list_feature_info
#' @description return the available giotto spatial feature information
#' @param gobject giotto object
#' @return names of available feature information
list_feature_info = function(gobject) {

  availableFeatInfo = data.table()
  for(info in names(gobject@feat_info)) {
    availableFeatInfo = rbind(availableFeatInfo,
                              list(feat_info = info))
  }

  if(nrow(availableFeatInfo) == 0) return(NULL)
  else return(availableFeatInfo)
}


#' @title list_feature_info_names
#' @name list_feature_info_names
#' @description return the available names for giotto feature information
#' @param gobject giotto object
#' @return vector with names of available feature information
list_feature_info_names = function(gobject) {

  feat_info_names = names(gobject@feat_info)

  return(feat_info_names)
}



#' @title list_spatial_networks
#' @name list_spatial_networks
#' @description return the available spatial networks that are attached to the Giotto object
#' @param gobject giotto object
#' @param spat_unit spatial unit (e.g. "cell")
#' @return data.table of names and locations of available spatial networks. col order matters
list_spatial_networks = function(gobject,
                                 spat_unit = NULL) {

  availableSpatNetworks = data.table()
  for(spatial_unit in names(gobject@spatial_network)) {
    for(spat_network_name in names(gobject@spatial_network[[spatial_unit]])) {
      availableSpatNetworks = rbind(availableSpatNetworks,
                                    list(spat_unit = spatial_unit,
                                         name = spat_network_name))
    }
  }

  # check if a specific category is desired
  if(!is.null(spat_unit)) spat_unit_subset = availableSpatNetworks$spat_unit == spat_unit else spat_unit_subset = TRUE

  availableSpatNetworks = availableSpatNetworks[spat_unit_subset,]

  if(nrow(availableSpatNetworks) == 0) return(NULL)
  else return(availableSpatNetworks)
}


#' @title list_spatial_networks_names
#' @name list_spatial_networks_names
#' @description return the available names for giotto feature information
#' @param gobject giotto object
#' @param spat_unit spatial unit (e.g. "cell")
#' @param feat_type feature type (e.g. "rna", "dna", "protein")
#' @return vector with names of available feature information
list_spatial_networks_names = function(gobject,
                                       spat_unit = NULL,
                                       feat_type = NULL) {

  spat_network_names = names(gobject@spatial_network[[spat_unit]][[feat_type]])

  return(spat_network_names)
}




#' @title list_spatial_grids
#' @name list_spatial_grids
#' @description return the available spatial grids that are attached to the Giotto object
#' @param gobject giotto object
#' @param spat_unit spatial unit (e.g. "cell")
#' @return data.table of names and locations of available spatial grids. col order matters
list_spatial_grids = function(gobject,
                              spat_unit = NULL) {

  availableSpatGrids = data.table()
  for(spatial_unit in names(gobject@spatial_grid)) {
    for(grid_names in names(gobject@spatial_grid[[spatial_unit]])) {
      availableSpatGrids = rbind(availableSpatGrids,
                                 list(spat_unit = spatial_unit,
                                      name = grid_names))
    }
  }

  # check if a specific category is desired
  if(!is.null(spat_unit)) {
    availableSpatGrids = availableSpatGrids[availableSpatGrids$spat_unit == spat_unit,]
  }

  if(nrow(availableSpatGrids) == 0) return(NULL)
  else return(availableSpatGrids)
}




#' @title list_spatial_grids_names
#' @name list_spatial_grids_names
#' @description return the available spatial grids name for a given spatial unit that are attached to the Giotto object
#' @param gobject giotto object
#' @param spat_unit spatial unit (e.g. "cell")
#' @return vector with names of available spatial grids names
list_spatial_grids_names = function(gobject,
                                    spat_unit = NULL) {

  spat_grid_names = names(gobject@spatial_grid[[spat_unit]])

  return(spat_grid_names)
}


#' @title list_images
#' @name list_images
#' @description Prints the available giotto images that are attached to the Giotto object
#' @param gobject giotto object
#' @param img_type "image" or "largeImage"
#' @return data.table of giotto image names attached to the giotto object
list_images = function(gobject,
                       img_type = NULL) {

  availableImages = data.table()

  g_image_names = names(gobject@images)
  g_limage_names = names(gobject@largeImages)

  for(image_type in c('image', 'largeImage')) {
    if(image_type == 'image') {
      for(name in g_image_names) {
        availableImages = rbind(availableImages,
                                list(img_type = image_type,
                                     name = name))
      }
    }
    if(image_type == 'largeImage') {
      for(name in g_limage_names) {
        availableImages = rbind(availableImages,
                                list(img_type = image_type,
                                     name = name))
      }
    }
  }

  # check if a specific category is desired
  if(!is.null(img_type)) img_type_subset = availableImages$img_type == img_type else img_type_subset = TRUE

  availableImages = availableImages[img_type_subset,]

  # NULL if there is no data
  if(nrow(availableImages) == 0) return(NULL)
  else return(availableImages)
}



#' @title list_images_names
#' @name list_images_names
#' @description return the available image names for a given image type that are attached to the Giotto object
#' @param gobject a giotto object
#' @param img_type "image" or "largeImage"
#' @return vector with names of available image names
list_images_names = function(gobject,
                             img_type) {

  if(!img_type %in% c('image', 'largeImage')) stop('img_type must be either "image" or "largeImage\n"')

  if(img_type == 'image') img_names = names(gobject@images)
  if(img_type == 'largeImage') img_names = names(gobject@largeImages)

  return(img_names)
}



