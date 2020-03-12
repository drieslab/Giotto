### update spatial network object structure

create_spatialNetworkObject <- function(name = NULL,
                             method = NULL,
                             parameters = NULL,
                             outputObj = NULL,
                             networkDT = NULL,
                             cellShapeObj = NULL,
                             crossSectionObjects = NULL,
                             misc = NULL) {

  networkObj = list(name = name,
                    method = method,
                    parameters = parameters,
                    outputObj = outputObj,
                    networkDT = networkDT,
                    cellShapeObj = cellShapeObj,
                    crossSectionObjects = crossSectionObjects,
                    misc = misc)

  class(networkObj) <- append(class(networkObj), "spatialNetworkObj")
  return(networkObj)

}

select_spatialNetwork <- function(gobject,
                                  name = NULL,
                                  return_network_Obj = FALSE) {
  if (!is.element(name,names(gobject@spatial_network))){
    stop(sprintf("spatial network %s has not been created.",name))
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

createSpatialNetwork <- function (gobject, k = 4, dimensions = "all", maximum_distance = NULL,
          minimum_k = 0, name = "spatial_network", verbose = F, return_gobject = TRUE,method="knn")
{
  spatial_locations = gobject@spatial_locs
  spatial_locations = spatial_locations[, grepl("sdim", colnames(spatial_locations)),
                                        with = F]
  if (dimensions != "all") {
    spatial_locations = spatial_locations[, dimensions, with = FALSE]
  }
  spatial_locations <- as.matrix(spatial_locations)
  rownames(spatial_locations) <- gobject@cell_ID
  cell_ID_vec <- c(1:nrow(spatial_locations))
  names(cell_ID_vec) <- rownames(spatial_locations)
  knn_spatial <- dbscan::kNN(x = spatial_locations, k = k)
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  outputObj = knn_spatial
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  knn_sptial.norm = data.frame(from = rep(1:nrow(knn_spatial$id),
                                          k), to = as.vector(knn_spatial$id), weight = 1/(1 + as.vector(knn_spatial$dist)),
                               distance = as.vector(knn_spatial$dist))
  nw_sptial.norm = igraph::graph_from_data_frame(knn_sptial.norm,
                                                 directed = FALSE)
  spatial_network_DT = data.table::as.data.table(knn_sptial.norm)
  spatial_network_DT[, `:=`(from, names(cell_ID_vec[from]))]
  spatial_network_DT[, `:=`(to, names(cell_ID_vec[to]))]
  cel_metadata = gobject@cell_metadata
  spatial_locs_DT = gobject@spatial_locs
  spatial_locations_annot <- spatial_locs_DT[, `:=`(cell_ID,
                                                    gobject@cell_ID)]
  spatial_network_DT <- merge(spatial_network_DT, by.x = "from",
                              spatial_locations_annot, by.y = "cell_ID")
  coord_names = colnames(spatial_locations)
  coord_begin = paste0(coord_names, "_begin")
  setnames(spatial_network_DT, coord_names, coord_begin)
  spatial_network_DT <- merge(spatial_network_DT, by.x = "to",
                              spatial_locations_annot, by.y = "cell_ID")
  coord_names = colnames(spatial_locations)
  coord_end = paste0(coord_names, "_end")
  setnames(spatial_network_DT, coord_names, coord_end)
  temp_fullnetwork = convert_to_full_spatial_network(spatial_network_DT)
  if (!is.null(maximum_distance)) {
    temp_fullnetwork = temp_fullnetwork[distance <= maximum_distance |
                                          rank_int <= minimum_k]
  }
  spatial_network_DT = convert_to_reduced_spatial_network(temp_fullnetwork)

  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

  parameters = list("maximum_distance" = maximum_distance,
                    "minimum_k" = minimum_k,
                    "k" = k,
                    "dimensions" = dimensions)

  spatial_network_Obj = create_spatialNetworkObject(name = name,
                                                     method = method,
                                                     parameters = parameters,
                                                     outputObj = outputObj,
                                                     networkDT = spatial_network_DT,
                                                     misc = NULL)

  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

  if (return_gobject == TRUE) {
    spn_names = names(gobject@spatial_network[[name]])
    if (name %in% spn_names) {
      cat("\n ", name, " has already been used, will be overwritten \n")
    }
    parameters_list = gobject@parameters
    number_of_rounds = length(parameters_list)
    update_name = paste0(number_of_rounds, "_spatial_network")
    parameters_list[[update_name]] = c(`k neighbours` = k,
                                       `dimensions used` = dimensions, `maximum distance threshold` = ifelse(is.null(maximum_distance),
                                                                                                             NA, maximum_distance), `name of spatial network` = name)
    gobject@parameters = parameters_list
    #gobject@spatial_network[[name]] = spatial_network_DT

    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
    gobject@spatial_network[[name]] = spatial_network_Obj
    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

    return(gobject)
  }
  else {
    return(spatial_network_DT)
  }
}

create_delaunayNetwork_geometry <- function(spatial_locations){

  delaunay_triangle = geometry::delaunayn(spatial_locations[,c("sdimx","sdimy")], options = "Pp")
  geometry_obj = list("delaunay_triangle"=delaunay_triangle)
  delaunay_edges <- rbind(delaunay_triangle[,c(1,2)],
                          delaunay_triangle[,c(1,3)],
                          delaunay_triangle[,c(2,3)])

  ### making sure of no duplication ###
  delaunay_edges_dedup = dplyr::distinct(as.data.frame(delaunay_edges))
  igraph_obj = graph_from_edgelist(as.matrix(delaunay_edges_dedup))
  adj_obj = as_adjacency_matrix(igraph_obj)
  igraph_obj2 = graph.adjacency(adj_obj)
  delaunay_edges_dedup2 = get.data.frame(igraph_obj2)
  delaunay_edges_dedup = delaunay_edges_dedup2
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  delaunay_network_DT = data.table::data.table(from = names(cell_ID_vec)[delaunay_edges_dedup$from],
                                               to = names(cell_ID_vec)[delaunay_edges_dedup$to],
                                               sdimx_begin = spatial_locations[delaunay_edges_dedup$from,"sdimx"],
                                               sdimy_begin = spatial_locations[delaunay_edges_dedup$from,"sdimy"],
                                               sdimx_end = spatial_locations[delaunay_edges_dedup$to,"sdimx"],
                                               sdimy_end = spatial_locations[delaunay_edges_dedup$to,"sdimy"])
  out_object = list("geometry_obj"=geometry_obj,
                    "delaunay_network_DT"=delaunay_network_DT)
  return(out_object)

}

create_delaunayNetwork_geometry_3D <- function(spatial_locations){

  delaunay_tetrahedra <- geometry::delaunayn(spatial_locations)
  geometry_obj = list("delaunay_tetrahedra"=delaunay_tetrahedra)
  delaunay_edges <- rbind(delaunay_tetrahedra[,c(1,2)],
                          delaunay_tetrahedra[,c(1,3)],
                          delaunay_tetrahedra[,c(1,4)],
                          delaunay_tetrahedra[,c(2,3)],
                          delaunay_tetrahedra[,c(2,4)],
                          delaunay_tetrahedra[,c(3,4)])

  ### making sure of no duplication ###
  delaunay_edges_dedup = dplyr::distinct(as.data.frame(delaunay_edges))
  igraph_obj = graph_from_edgelist(as.matrix(delaunay_edges_dedup))
  adj_obj = as_adjacency_matrix(igraph_obj)
  igraph_obj2 = graph.adjacency(adj_obj)
  delaunay_edges_dedup2 = get.data.frame(igraph_obj2)
  delaunay_edges_dedup = delaunay_edges_dedup2
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  delaunay_network_DT = data.table::data.table(from = names(cell_ID_vec)[delaunay_edges_dedup$from],
                                               to = names(cell_ID_vec)[delaunay_edges_dedup$to],
                                               sdimx_begin = spatial_locations[delaunay_edges_dedup$from,"sdimx"],
                                               sdimy_begin = spatial_locations[delaunay_edges_dedup$from,"sdimy"],
                                               sdimz_begin = spatial_locations[delaunay_edges_dedup$from,"sdimz"],
                                               sdimx_end = spatial_locations[delaunay_edges_dedup$to,"sdimx"],
                                               sdimy_end = spatial_locations[delaunay_edges_dedup$to,"sdimy"],
                                               sdimz_end = spatial_locations[delaunay_edges_dedup$to,"sdimz"])

  out_object = list("geometry_obj"=geometry_obj,
                    "delaunay_network_DT"=delaunay_network_DT)
  return(out_object)

}

create_delaunayNetwork_RTriangle <- function(spatial_locations,Y=TRUE,j=TRUE,S=0){

  RTriangle_obj = RTriangle::triangulate(RTriangle::pslg(spatial_locations),
                                         Y = Y, j = j, S = S, ...)

  delaunay_network_DT = data.table::data.table(from = names(cell_ID_vec)[RTriangle_obj$E[,
                                                                                         1]], to = names(cell_ID_vec)[RTriangle_obj$E[, 2]], sdimx_begin = RTriangle_obj$P[RTriangle_obj$E[,
                                                                                                                                                                                           1], 1], sdimy_begin = RTriangle_obj$P[RTriangle_obj$E[,
                                                                                                                                                                                                                                                 1], 2], sdimx_end = RTriangle_obj$P[RTriangle_obj$E[,
                                                                                                                                                                                                                                                                                                     2], 1], sdimy_end = RTriangle_obj$P[RTriangle_obj$E[,
                                                                                                                                                                                                                                                                                                                                                         2], 2])


  out_object = list("RTriangle_obj"=RTriangle_obj,
                    "delaunay_network_DT"=delaunay_network_DT)
  return(out_object)
}

create_delaunayNetwork_deldir <- function(spatial_locations){

  deldir_obj = deldir::deldir(x = spatial_locations[,"sdimx"],y = spatial_locations[,"sdimy"])

  delaunay_network_DT = data.table::data.table(from = names(cell_ID_vec)[deldir_obj$delsgs$ind1],
                                               to = names(cell_ID_vec)[deldir_obj$delsgs$ind2],
                                               sdimx_begin = deldir_obj$delsgs$x1,
                                               sdimy_begin = deldir_obj$delsgs$y1,
                                               sdimx_end = deldir_obj$delsgs$x2,
                                               sdimy_end = deldir_obj$delsgs$y2)
  out_object = list("deldir_obj"=deldir_obj,
                    "delaunay_network_DT"=delaunay_network_DT)
  return(out_object)
}

createDelaunayNetwork <- function (gobject, method = "delaunayn_geometry", dimensions = c("sdimx", "sdimy"), name = "delaunay_network",
          maximum_distance = "auto", minimum_k = 0, Y = TRUE, j = TRUE,
          S = 0, verbose = T, return_gobject = TRUE, ...)
{
  spatial_locations = gobject@spatial_locs
  spatial_locations = spatial_locations[, grepl("sdim", colnames(spatial_locations)),
                                        with = F]
  if (length(dimensions) != 2) {
    stop("A Delaunay network can only be computed on 2 dimensions \n")
  }
  else {
    spatial_locations = spatial_locations[, dimensions, with = F]
  }
  spatial_locations = as.matrix(spatial_locations)
  rownames(spatial_locations) = gobject@cell_ID
  cell_ID_vec = c(1:nrow(spatial_locations))
  names(cell_ID_vec) = rownames(spatial_locations)

  if (method == "RTriangle"){

    delaunay_output = create_delaunayNetwork_RTriangle(spatial_locations,Y = Y, j = j, S = S)

    RTriangle_obj = delaunay_output$RTriangle_obj
    delaunay_network_DT = delaunay_output$delaunay_network_DT
    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
    parameters = list("maximum_distance" = maximum_distance,
                      "minimum_k" = minimum_k,
                      "Y" = Y,
                      "j" = j,
                      "S" = S)
    outputObj = RTriangle_obj
    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###


  }else if (method == "deldir"){

    delaunay_output = create_delaunayNetwork_RTriangle(spatial_locations,Y = Y, j = j, S = S)

    deldir_obj = delaunay_output$deldir_obj
    delaunay_network_DT = delaunay_output$delaunay_network_DT
    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
    parameters = list("maximum_distance" = maximum_distance,
                      "minimum_k" = minimum_k)
    outputObj = deldir_obj
    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

  }else if (method == "delaunayn_geometry"){

    delaunay_output = create_delaunayNetwork_geometry(spatial_locations)

    outputObj = delaunay_output$geometry_obj
    delaunay_network_DT = delaunay_output$delaunay_network_DT
  }


  delaunay_network_DT[, `:=`(distance, dist(x = matrix(c(sdimx_begin,
                                                         sdimy_begin, sdimx_end, sdimy_end), nrow = 2, byrow = T))),
                      by = 1:nrow(delaunay_network_DT)]
  delaunay_network_DT[, `:=`(distance, as.numeric(distance))]
  delaunay_network_DT[, `:=`(weight, 1/distance)]
  data.table::setorder(delaunay_network_DT, from, distance)
  delaunay_network_DT = delaunay_network_DT[, .(to, from, weight,
                                                distance, sdimx_begin, sdimy_begin, sdimx_end, sdimy_end)]
  temp_fullnetwork = Giotto:::convert_to_full_spatial_network(delaunay_network_DT)
  if (maximum_distance == "auto") {
    temp_fullnetwork = temp_fullnetwork[distance <= boxplot.stats(temp_fullnetwork$distance)$stats[5]]
  }
  else if (!is.null(maximum_distance)) {
    temp_fullnetwork = temp_fullnetwork[distance <= maximum_distance |
                                          rank_int <= minimum_k]
  }

  delaunay_network_DT = Giotto:::convert_to_reduced_spatial_network(temp_fullnetwork)

  ### calcualte cell diameter ###
  cellDiameter_mean = mean(delaunay_network_DT$distance)
  cellDiameter_median = median(delaunay_network_DT$distance)

  cellShapeObj = list("cellDiameter_mean" = cellDiameter_mean,
                      "cellDiameter_median" = cellDiameter_median
  )


  ###
  ###
  delaunay_network_Obj = create_spatialNetworkObject(name = name,
                                                     method = method,
                                                     parameters = parameters,
                                                     outputObj = outputObj,
                                                     networkDT = delaunay_network_DT,
                                                     cellShapeObj = cellShapeObj,
                                                     misc = NULL)
  ###
  ###

  if (return_gobject == TRUE) {
    spn_names = names(gobject@spatial_network[[name]])
    if (name %in% spn_names) {
      cat("\n ", name, " has already been used, will be overwritten \n")
    }
    parameters_list = gobject@parameters
    number_of_rounds = length(parameters_list)
    update_name = paste0(number_of_rounds, "_delaunay_spatial_network")
    parameters_list[[update_name]] = c(`dimensions used` = paste(dimensions,
                                                                 collapse = "-"), `maximum distance threshold` = ifelse(is.null(maximum_distance),
                                                                                                                        NA, maximum_distance), `RTriangle Y:` = Y, `RTriangle j:` = j,
                                       `RTriangle S:` = S, `name of spatial network` = name)
    gobject@parameters = parameters_list
    # gobject@spatial_network[[name]] = delaunay_network_DT

    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
    gobject@spatial_network[[name]] = delaunay_network_Obj
    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###


    return(gobject)
  }
  else {
    return(delaunay_network_DT)
  }
}

createDelaunayNetwork3D <- function (gobject, method = "delaunayn_geometry", dimensions = c("sdimx", "sdimy","sdimz"), name = "delaunay_network_3D",
                                   maximum_distance = "auto", return_gobject = TRUE, ...)
{
  spatial_locations = gobject@spatial_locs
  spatial_locations = spatial_locations[, grepl("sdim", colnames(spatial_locations)),
                                        with = F]
  spatial_locations = as.matrix(spatial_locations)
  rownames(spatial_locations) = gobject@cell_ID
  cell_ID_vec = c(1:nrow(spatial_locations))
  names(cell_ID_vec) = rownames(spatial_locations)

##
  if (method == "delaunayn_geometry"){

    delaunay_output = create_delaunayNetwork_geometry_3D(spatial_locations)

    outputObj = delaunay_output$geometry_obj
    delaunay_network_DT = delaunay_output$delaunay_network_DT

  }
  delaunay_network_DT[, `:=`(distance, dist(x = matrix(c(sdimx_begin,
                                                         sdimy_begin, sdimx_end, sdimy_end), nrow = 2, byrow = T))),
                      by = 1:nrow(delaunay_network_DT)]
  delaunay_network_DT[, `:=`(distance, as.numeric(distance))]
  delaunay_network_DT[, `:=`(weight, 1/distance)]
  data.table::setorder(delaunay_network_DT, from, distance)
  delaunay_network_DT = delaunay_network_DT[, .(to, from, weight,
                                                distance, sdimx_begin, sdimy_begin, sdimx_end, sdimy_end)]
  temp_fullnetwork = Giotto:::convert_to_full_spatial_network(delaunay_network_DT)
  if (maximum_distance == "auto") {
    temp_fullnetwork = temp_fullnetwork[distance <= boxplot.stats(temp_fullnetwork$distance)$stats[5]]
  }
  else if (!is.null(maximum_distance)) {
    temp_fullnetwork = temp_fullnetwork[distance <= maximum_distance |
                                          rank_int <= minimum_k]
  }
  delaunay_network_DT = Giotto:::convert_to_reduced_spatial_network(temp_fullnetwork)
  cellDiameter_mean = mean(delaunay_network_DT$distance)
  cellDiameter_median = median(delaunay_network_DT$distance)

  cellShapeObj = list("cellDiameter_mean" = cellDiameter_mean,
                      "cellDiameter_median" = cellDiameter_median
                      )

  if (return_gobject == TRUE) {
    spn_names = names(gobject@spatial_network[[name]])
    if (name %in% spn_names) {
      cat("\n ", name, " has already been used, will be overwritten \n")
    }
    parameters_list = gobject@parameters
    number_of_rounds = length(parameters_list)
    update_name = paste0(number_of_rounds, "_delaunay_spatial_network_3D")
    parameters_list[[update_name]] = c(`dimensions used` = paste(dimensions,
                                                                 collapse = "-"), `maximum distance threshold` = ifelse(is.null(maximum_distance),
                                                                                                                        NA, maximum_distance), `name of spatial network` = name)
    gobject@parameters = parameters_list
    # gobject@spatial_network[[name]] = delaunay_network_DT

    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
    ###
    ###
    delaunay_network_Obj = create_spatialNetworkObject(name = name,
                                                       method = method,
                                                       #parameters = parameters,
                                                       outputObj = outputObj,
                                                       networkDT = delaunay_network_DT,
                                                       cellShapeObj = cellShapeObj,
                                                       misc = NULL)
    ###
    ###
    gobject@spatial_network[[name]] = delaunay_network_Obj
    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###


    return(gobject)
  }
  else {
    return(delaunay_network_DT)
  }
}
