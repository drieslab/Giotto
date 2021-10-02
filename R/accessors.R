


## Get and set functions to get and set values in one of the giotto class slots ##


## expression values slot ####

#' @name  get_expression_values
#' @description function to get expression values from giotto object
#' @param gobject giotto object
#' @param feat_type feature type
#' @param values expression values to extract
#' @return expression matrix
#' @keywords internal
get_expression_values <- function(gobject,
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


#' @name select_expression_values
#' @inheritDotParams get_expression_values
#' @seealso \code{\link{get_expression_values}}
#' @keywords internal
select_expression_values = function(...) {

  .Deprecated(new = "get_expression_values")

  get_expression_values(...)

}



#' @name  set_expression_values
#' @description function to set expression values for giotto object
#' @param gobject giotto object
#' @param feat_type feature type
#' @param name name for the expression slot
#' @param values expression values
#' @return giotto object
#' @keywords internals
set_expression_values <- function(gobject,
                                  feat_type = 'rna',
                                  name = 'test',
                                  values) {


  ## 1. check if specified name has already been used
  potential_names = names(gobject@expression[[feat_type]])
  if(name %in% potential_names) {
    cat(name, ' already exist and will be replaced with new values \n')
  }

  ## TODO: 2. check input for values matrix


  ## 3. update and return giotto object
  gobject@expression[[feat_type]][[name]] = values
  return(gobject)

}





## spatial locations slot ####


#' @name get_spatial_locations
#' @description function to get a spatial location data.table
#' @param gobject giotto object
#' @param spat_loc_name name of spatial locations
#' @return data.table with coordinates
#' @keywords internal
get_spatial_locations <- function(gobject,
                                  spat_loc_name = 'raw') {


  # spatial locations
  # if NULL (not given) and spatial locations have been added, then use first one
  # if NULL (not given) and spatial loactions have NOT been added, then keep NULL
  if(is.null(spat_loc_name)) {
    if(!is.null(gobject@spatial_locs)) {
      spat_loc_name = names(gobject@spatial_locs)[[1]]
      # cat('No spatial locations have been selected, the first one -',spat_loc_name, '- will be used \n')
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




#' @name select_expression_values
#' @inheritDotParams get_spatial_locations
#' @seealso \code{\link{get_spatial_locations}}
#' @keywords internal
select_spatial_locations = function(...) {

  .Deprecated(new = "get_spatial_locations")

  get_spatial_locations(...)

}




#' @name set_spatial_locations
#' @description function to set a spatial location slot
#' @param gobject giotto object
#' @param spat_loc_name name of spatial locations
#' @param spatlocs spatial locations
#' @param verbose be verbose
#' @return giotto object
#' @keywords internal
set_spatial_locations <- function(gobject,
                                  spat_loc_name = 'raw',
                                  spatlocs,
                                  verbose = TRUE) {


  ## 1. check if specified name has already been used
  potential_names = names(gobject@spatial_locs[[spat_loc_name]])
  if(spat_loc_name %in% potential_names) {
    if(verbose == TRUE) {
    cat(spat_loc_name, ' already exists and will be replaced with new spatial locations \n')
    }
  }

  ## TODO: 2. check input for spatial locations


  ## 3. update and return giotto object
  gobject@spatial_locs[[spat_loc_name]] = spatlocs
  return(gobject)

}






## dimension reduction slot ####


#' @name get_dimReduction
#' @description function to get a dimension reduction object
#' @keywords internal
#' @return dim reduction coordinates (default) or dim reduction object
get_dimReduction = function(gobject,
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


#' @name select_dimReduction
#' @inheritDotParams get_dimReduction
#' @seealso \code{\link{get_dimReduction}}
#' @keywords internal
select_dimReduction = function(...) {

  .Deprecated(new = "get_dimReduction")

  get_dimReduction(...)

}



#' @name set_dimReduction
#' @description function to set a dimension reduction slot
#' @param gobject giotto object
#' @param reduction reduction with cells or genes
#' @param reduction_method reduction name
#' @param name reduction name
#' @param dimObject dimension object
#' @return giotto object
#' @keywords internal
set_dimReduction <- function(gobject,
                             reduction = c('cells', 'genes'),
                             reduction_method = c('pca', 'umap', 'tsne'),
                             name = 'pca',
                             dimObject) {


  ## 1. check if specified name has already been used
  potential_names = names(gobject@dimension_reduction[[reduction]][[reduction_method]])
  if(name %in% potential_names) {
    cat(name, ' already exist and will be replaced with new dimension reduction object \n')
  }

  ## TODO: 2. check input for dimension reduction object


  ## 3. update and return giotto object
  gobject@dimension_reduction[[reduction]][[reduction_method]][[name]] = dimObject
  return(gobject)

}







## nearest neighbor network slot ####


#' @name get_NearestNetwork
#' @description get a NN-network from a Giotto object
#' @param gobject giotto object
#' @param nn_network_to_use kNN or sNN
#' @param network_name name of NN network to be used
#' @param output return a igraph or data.table object
#' @return igraph or data.table object
#' @keywords internal
get_NearestNetwork = function(gobject,
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
#' @inheritDotParams get_NearestNetwork
#' @seealso \code{\link{get_NearestNetwork}}
#' @keywords internal
extractNearestNetwork = function(...) {

  .Deprecated(new = "get_NearestNetwork")

  get_NearestNetwork(...)

}


#' @name select_NearestNetwork
#' @inheritDotParams get_NearestNetwork
#' @seealso \code{\link{get_NearestNetwork}}
#' @keywords internal
select_NearestNetwork = function(...) {

  .Deprecated(new = "get_NearestNetwork")

  get_NearestNetwork(...)

}



#' @name set_NearestNetwork
#' @description set a NN-network for a Giotto object
#' @param gobject giotto object
#' @param nn_network_to_use kNN or sNN
#' @param network_name name of NN network to be used
#' @param nn_network nearest network
#' @return giotto object
#' @keywords internal
set_NearestNetwork = function(gobject,
                              nn_network_to_use = 'sNN',
                              network_name = 'sNN.pca',
                              nn_network) {


  ## 1. check if specified name has already been used
  potential_names = names(gobject@nn_network[[nn_network_to_use]])
  if(network_name %in% potential_names) {
    cat(network_name, ' already exist and will be replaced with nearest neighbor network \n')
  }

  ## TODO: 2. check input for nearest neighbor network
  # convert to igraph if data.table class


  ## 3. update and return giotto object
  gobject@nn_network[[nn_network_to_use]][[network_name]] = nn_network
  return(gobject)

}






## spatial network slot ####


#' @name get_spatialNetwork
#' @description function to get a spatial network
#' @param gobject giotto object
#' @param name name of spatial network
#' @param return_network_Obj return network object (default = FALSE)
#' @keywords internal
get_spatialNetwork <- function(gobject,
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



#' @name select_spatialNetwork
#' @inheritDotParams get_spatialNetwork
#' @seealso \code{\link{get_spatialNetwork}}
#' @keywords internal
select_spatialNetwork = function(...) {

  .Deprecated(new = "get_spatialNetwork")

  get_spatialNetwork(...)

}



#' @name set_spatialNetwork
#' @description function to set a spatial network
#' @param gobject giotto object
#' @param name name of spatial network
#' @param spatial_network spatial network
#' @keywords internal
set_spatialNetwork <- function(gobject,
                               name = NULL,
                               spatial_network) {



  ## 1. check if specified name has already been used
  potential_names = names(gobject@spatial_network)
  if(name %in% potential_names) {
    cat(name, ' already exist and will be replaced with new spatial network \n')
  }

  ## TODO: 2. check input for spatial network


  ## 3. update and return giotto object
  gobject@spatial_network[[name]] = spatial_network
  return(gobject)

}




## spatial grid slot ####


#' @name get_spatialGrid
#' @description function to get spatial grid
#' @param gobject giotto object
#' @param name name of spatial grid
#' @param return_network_Obj return grid object (default = FALSE)
#' @keywords internal
get_spatialGrid <- function(gobject,
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



#' @name select_spatialGrid
#' @inheritDotParams get_spatialGrid
#' @seealso \code{\link{get_spatialGrid}}
#' @keywords internal
select_spatialGrid = function(...) {

  .Deprecated(new = "get_spatialGrid")

  get_spatialGrid(...)

}



#' @name set_spatialGrid
#' @description function to set a spatial grid
#' @param gobject giotto object
#' @param name name of spatial grid
#' @param spatial_grid spatial grid
#' @keywords internal
set_spatialGrid <- function(gobject,
                            name = NULL,
                            spatial_grid) {



  ## 1. check if specified name has already been used
  potential_names = names(gobject@spatial_grid)
  if(name %in% potential_names) {
    cat(name, ' already exist and will be replaced with new spatial grid \n')
  }

  ## TODO: 2. check input for spatial grid


  ## 3. update and return giotto object
  gobject@spatial_grid[[name]] = spatial_grid
  return(gobject)

}








## polygon cell info ####


#' @name get_polygon_info
#' @description  get giotto polygon spatVector
#' @param gobject giotto object
#' @param polygon_name name of polygons
#' @param polygon_overlap include polygon overlap information
#' @keywords internal
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



#' @name select_polygon_info
#' @inheritDotParams get_polygon_info
#' @seealso \code{\link{get_polygon_info}}
#' @keywords internal
select_polygon_info = function(...) {

  .Deprecated(new = "get_polygon_info")

  get_polygon_info(...)

}




#' @name set_polygon_info
#' @description  set giotto polygon spatVector
#' @param gobject giotto object
#' @param polygon_name name of polygons
#' @param gpolygon giotto polygon
#' @keywords internal
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

#' @name get_feature_info
#' @description  get giotto points spatVector
#' @param gobject giotto object
#' @param feat_name name of feature
#' @keywords internal
get_feature_info = function(gobject,
                               feat_name = 'rna') {

  potential_names = names(gobject@feat_info)

  if(!feat_name %in% potential_names) {
    stop('There is no feature information with name ', feat_name, '\n')
  } else {
    feat_info = gobject@feat_info[[feat_name]]@spatVector
    return(feat_info)
  }
}


#' @name select_feature_info
#' @inheritDotParams get_feature_info
#' @seealso \code{\link{get_feature_info}}
#' @keywords internal
select_feature_info = function(...) {

  .Deprecated(new = "get_feature_info")

  get_feature_info(...)

}



#' @name set_feature_info
#' @description  set giotto polygon spatVector for features
#' @param gobject giotto object
#' @param feat_name name of feat
#' @param gpolygon giotto polygon
#' @keywords internal
set_feature_info = function(gobject,
                            feat_name = 'cell',
                            gpolygon) {



  ## 1. check if specified name has already been used
  potential_names = names(gobject@feat_info)
  if(feat_name %in% potential_names) {
    cat(feat_name, ' already exist and will be replaced with new giotto polygon \n')
  }

  ## TODO: 2. check input for giotto polygon


  ## 3. update and return giotto object
  gobject@feat_name[[feat_name]] = gpolygon
  return(gobject)

}





## spatial enrichment slot ####


#' @name get_spatial_enrichment
#' @description function to get a spatial enrichment data.table
#' @param gobject giotto object
#' @param enrichm_name name of spatial enrichment results
#' @return data.table with fractions
#' @keywords internal
get_spatial_enrichment <- function(gobject,
                                   enrichm_name = 'DWLS') {


  # spatial locations
  # if NULL (not given) and spatial locations have been added, then use first one
  # if NULL (not given) and spatial loactions have NOT been added, then keep NULL
  if(is.null(enrichm_name)) {
    if(!is.null(gobject@spatial_enrichment)) {
      enrichm_name = names(gobject@spatial_enrichment)[[1]]
      # cat('No spatial locations have been selected, the first one -',spat_loc_name, '- will be used \n')
    } else {
      enrichm_name = NULL
      cat('No spatial enrichment results have been found \n')
      return(NULL)
    }
  }

  potential_names = names(gobject@spatial_enrichment)

  if(enrichm_name %in% potential_names) {
    enr_res = data.table::copy(gobject@spatial_enrichment[[enrichm_name]])
    return(enr_res)
  } else {
    stop("The spatial enrichment result with name ","'", enrichm_name, "'"," can not be found \n")
  }
}


#' @name set_spatial_enrichment
#' @description function to set a spatial enrichment slot
#' @param gobject giotto object
#' @param enrichm_name name of spatial enrichment results
#' @param spatenrichment spatial enrichment results
#' @return giotto object
#' @keywords internal
set_spatial_enrichment <- function(gobject,
                                   enrichm_name = 'enrichment',
                                   spatenrichment) {


  ## 1. check if specified name has already been used
  potential_names = names(gobject@spatial_enrichment[[enrichm_name]])
  if(enrichm_name %in% potential_names) {
    cat(enrichm_name, ' already exist and will be replaced with new spatial enrichment results \n')
  }

  ## TODO: 2. check input for spatial locations


  ## 3. update and return giotto object
  gobject@spatial_enrichment[[enrichm_name]] = spatenrichment
  return(gobject)

}






## Show functions ####


#' @name showGiottoExpression
#' @description shows the available matrices
#' @param gobject giotto object
#' @param nrows number of rows to print for each matrix
#' @param ncols number of columns to print for each matrix
#' @return prints the name and small subset of available matrices
#' @export
showGiottoExpression = function(gobject, nrows = 4, ncols = 4) {

  for(feat_type in names(gobject@expression)) {

    cat('Feature ', feat_type, ': \n\n')

    for(mat_i in names(gobject@expression[[feat_type]])) {

      cat('---> ', mat_i, 'matrix: \n')

      print(gobject@expression[[feat_type]][[mat_i]][1:nrows, 1:ncols])
      cat('\n')
    }

  }

}


#' @name showGiottoSpatLocs
#' @description shows the available spatial locations
#' @param gobject giotto object
#' @param nrows number of rows to print for each spatial location data.table
#' @return prints the name and small subset of available data.table
#' @export
showGiottoSpatLocs = function(gobject, nrows = 4) {

  for(spatlocname in names(gobject@spatial_locs)) {
    cat('Name ', spatlocname, ': \n\n')
    print(gobject@spatial_locs[[spatlocname]][1:nrows,])
  }
}




#' @name showGiottoSpatEnrichments
#' @description shows the available spatial enrichment results
#' @param gobject giotto object
#' @param nrows number of rows to print for each spatial enrichment data.table
#' @return prints the name and small subset of available data.table
#' @export
showGiottoSpatEnrichments = function(gobject, nrows = 4) {

  for(spatenrichname in names(gobject@spatial_enrichment)) {
    cat('Name ', spatenrichname, ': \n\n')
    print(gobject@spatial_enrichment[[spatenrichname]][1:nrows,])
  }
}




#' @name showGiottoDimRed
#' @description shows the available dimension reductions
#' @param gobject giotto object
#' @param nrows number of coordinates rows to print
#' @param ncols number of coordinates columns to print
#' @return prints the name and small subset of available dimension reduction coordinates
#' @export
showGiottoDimRed = function(gobject,
                            nrows = 3,
                            ncols = 2) {


  # for features
  cat('Dim reduction on features:',
      '\n',
      '-------------------------',
      '\n\n\n')

  for(dim_type in names(gobject@dimension_reduction[['feats']])) {

    cat('Dim reduction ', dim_type, ': \n\n')

    for(sub_type in names(gobject@dimension_reduction[['feats']][[dim_type]])) {

      cat('---> ', sub_type, 'coordinates: \n')

      print(gobject@dimension_reduction[['feats']][[dim_type]][[sub_type]][['coordinates']][1:nrows, 1:ncols])
      cat('\n')
    }

  }

  # for cells
  cat('Dim reduction on cells:',
      '\n',
      '----------------------',
      '\n\n\n')
  for(dim_type in names(gobject@dimension_reduction[['cells']])) {

    cat('Dim reduction ', dim_type, ': \n\n')

    for(sub_type in names(gobject@dimension_reduction[['cells']][[dim_type]])) {

      cat('---> ', sub_type, 'coordinates: \n')

      print(gobject@dimension_reduction[['cells']][[dim_type]][[sub_type]][['coordinates']][1:nrows, 1:ncols])
      cat('\n')
    }

  }

}






#' @name showGiottoSpatialInfo
#' @description show the available giotto spatial polygon information
#' @param gobject giotto object
#' @keywords show
#' @export
showGiottoSpatialInfo = function(gobject) {

  for(info in names(gobject@spatial_info)) {

    cat("For Spatial info: ", info, "\n\n")
    print(gobject@spatial_info[[info]])
    cat("-----------------------------")
    cat("\n \n")
  }

}


#' @name showGiottoFeatInfo
#' @description show the available giotto spatial feature information
#' @param gobject giotto object
#' @keywords show
#' @export
showGiottoFeatInfo = function(gobject) {

  for(info in names(gobject@feat_info)) {

    cat("For Feature info: ", info, "\n\n")
    print(gobject@feat_info[[info]])
    cat("-----------------------------")
    cat("\n \n")
  }

}



#' @title showNetworks
#' @description Prints the available spatial networks that are attached to the Giotto object
#' @param gobject a giotto object
#' @param verbose verbosity of function#'
#' @return vector
#' @export
showNetworks = function(gobject,
                        verbose = TRUE) {

  if(is.null(gobject)) stop('A giotto object needs to be provided \n')
  g_network_names = names(gobject@spatial_network)

  if(verbose == TRUE) {
    cat('The following images are available: ',
        g_network_names, '\n')
  }

  return(g_network_names)
}


#' @title showGrids
#' @description Prints the available spatial grids that are attached to the Giotto object
#' @param gobject a giotto object
#' @param verbose verbosity of function#'
#' @return vector
#' @export
showGrids = function(gobject,
                     verbose = TRUE) {

  if(is.null(gobject)) stop('A giotto object needs to be provided \n')
  g_grid_names = names(gobject@spatial_grid)

  if(verbose == TRUE) {
    cat('The following grids are available: ',
        g_grid_names, '\n')
  }

  return(g_grid_names)
}






