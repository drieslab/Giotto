
## Spatial structure helper functions ####



#' @title spatNetwDistributionsDistance
#' @description This function return histograms displaying the distance distribution for each spatial k-neighbor
#' @param gobject Giotto object
#' @param spatial_network_name name of spatial network
#' @param hist_bins number of binds to use for the histogram
#' @param test_distance_limit effect of different distance threshold on k-neighbors
#' @param ncol number of columns to visualize the histograms in
#' @param show_plot show plot
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters from \code{\link{all_plots_save_function}}
#' @param default_save_name default save name for saving, alternatively change save_name in save_param
#' @return ggplot plot
#' @export
spatNetwDistributionsDistance <- function(gobject,
                                          spatial_network_name = 'spatial_network',
                                          hist_bins = 30,
                                          test_distance_limit =  NULL,
                                          ncol = 1,
                                          show_plot = NA,
                                          return_plot = NA,
                                          save_plot = NA,
                                          save_param =  list(),
                                          default_save_name = 'spatNetwDistributionsDistance') {


  # data.table variables
  distance = rank_int = status = label = keep = NULL

  ## spatial network
  #spatial_network = gobject@spatial_network[[spatial_network_name]]
  spatial_network = select_spatialNetwork(gobject,name = spatial_network_name,return_network_Obj = FALSE)

  ## convert to full network with rank_int column
  spatial_network = convert_to_full_spatial_network(spatial_network)

  if(is.null(spatial_network)) {
    stop('spatial network ', spatial_network_name, ' was not found')
  }

  if(!is.null(test_distance_limit)) {
    removed_neighbors = spatial_network[distance > test_distance_limit, .N, by = rank_int]
    removed_neighbors[, status := 'remove']
    keep_neighbors = spatial_network[distance <= test_distance_limit, .N, by = rank_int]
    keep_neighbors[, status := 'keep']

    dist_removal_dt = rbind(removed_neighbors, keep_neighbors)
    setorder(dist_removal_dt, rank_int)

    dist_removal_dt_dcast = dcast.data.table(data = dist_removal_dt, rank_int~status, value.var = 'N', fill = 0)
    dist_removal_dt_dcast[, label := paste0('keep:',keep, '\n remove:',remove)]
  }

  # text location coordinates
  middle_distance = max(spatial_network$distance)/(3/2)
  freq_dt = spatial_network[, table(cut(distance, breaks = 30)), by = rank_int]
  middle_height = max(freq_dt$V1)/(3/2)

  pl = ggplot()
  pl = pl + labs(title = 'distance distribution per k-neighbor')
  pl = pl + theme_classic()
  pl = pl + geom_histogram(data = spatial_network, aes(x = distance), color = 'white', fill = 'black', bins = hist_bins)
  pl = pl + facet_wrap(~rank_int, ncol = ncol)
  if(!is.null(test_distance_limit)) {
    pl = pl + geom_vline(xintercept = test_distance_limit, color = 'red')
    pl = pl + geom_text(data = dist_removal_dt_dcast, aes(x = middle_distance, y = middle_height, label = label))
  }

  # print, return and save parameters
  show_plot = ifelse(is.na(show_plot), readGiottoInstructions(gobject, param = 'show_plot'), show_plot)
  save_plot = ifelse(is.na(save_plot), readGiottoInstructions(gobject, param = 'save_plot'), save_plot)
  return_plot = ifelse(is.na(return_plot), readGiottoInstructions(gobject, param = 'return_plot'), return_plot)

  ## print plot
  if(show_plot == TRUE) {
    print(pl)
  }

  ## save plot
  if(save_plot == TRUE) {
    do.call('all_plots_save_function', c(list(gobject = gobject, plot_object = pl, default_save_name = default_save_name), save_param))
  }

  ## return plot
  if(return_plot == TRUE) {
    return(pl)
  }


}




#' @title spatNetwDistributionsKneighbors
#' @description This function returns a histogram displaying the number of k-neighbors distribution for each cell
#' @param gobject Giotto object
#' @param spatial_network_name name of spatial network
#' @param hist_bins number of binds to use for the histogram
#' @param show_plot show plot
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters from \code{\link{all_plots_save_function}}
#' @param default_save_name default save name for saving, alternatively change save_name in save_param
#' @return ggplot plot
#' @export
spatNetwDistributionsKneighbors = function(gobject,
                                           spatial_network_name = 'spatial_network',
                                           hist_bins = 30,
                                           show_plot = NA,
                                           return_plot = NA,
                                           save_plot = NA,
                                           save_param =  list(),
                                           default_save_name = 'spatNetwDistributionsKneighbors') {

  # data.table variables
  N = NULL

  ## spatial network
  #spatial_network = gobject@spatial_network[[spatial_network_name]]
  spatial_network = select_spatialNetwork(gobject,name = spatial_network_name,return_network_Obj = FALSE)

  ## convert to full network with rank_int column
  spatial_network = convert_to_full_spatial_network(spatial_network)

  if(is.null(spatial_network)) {
    stop('spatial network ', spatial_network_name, ' was not found')
  }

  spatial_network_dt = as.data.table(spatial_network[, table(source)])

  pl = ggplot()
  pl = pl + labs(title = 'k-neighbor distribution for all cells', x = 'k-neighbors/cell')
  pl = pl + theme_classic()
  pl = pl + geom_histogram(data = spatial_network_dt, aes(x = N), color = 'white', fill = 'black', bins = hist_bins)


  # print, return and save parameters
  show_plot = ifelse(is.na(show_plot), readGiottoInstructions(gobject, param = 'show_plot'), show_plot)
  save_plot = ifelse(is.na(save_plot), readGiottoInstructions(gobject, param = 'save_plot'), save_plot)
  return_plot = ifelse(is.na(return_plot), readGiottoInstructions(gobject, param = 'return_plot'), return_plot)

  ## print plot
  if(show_plot == TRUE) {
    print(pl)
  }

  ## save plot
  if(save_plot == TRUE) {
    do.call('all_plots_save_function', c(list(gobject = gobject, plot_object = pl, default_save_name = default_save_name), save_param))
  }

  ## return plot
  if(return_plot == TRUE) {
    return(pl)
  }


}



#' @title spatNetwDistributionsDistance
#' @description This function return histograms displaying the distance distribution for each spatial k-neighbor
#' @param gobject Giotto object
#' @param spatial_network_name name of spatial network
#' @param distribution show the distribution of cell-to-cell distance or number of k neighbors
#' @param hist_bins number of binds to use for the histogram
#' @param test_distance_limit effect of different distance threshold on k-neighbors
#' @param ncol number of columns to visualize the histograms in
#' @param show_plot show plot
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters from \code{\link{all_plots_save_function}}
#' @param default_save_name default save name for saving, alternatively change save_name in save_param
#' @details The \strong{distance} option shows the spatial distance distribution for each nearest neighbor rank (1st, 2nd, 3th, ... neigbor).
#' With this option the user can also test the effect of a distance limit on the spatial network. This distance limit can be used to remove neigbor
#' cells that are considered to far away. \cr
#' The \strong{k_neighbors} option shows the number of k neighbors distribution over all cells.
#' @return ggplot plot
#' @export
spatNetwDistributions <- function(gobject,
                                  spatial_network_name = 'spatial_network',
                                  distribution = c('distance', 'k_neighbors'),
                                  hist_bins = 30,
                                  test_distance_limit =  NULL,
                                  ncol = 1,
                                  show_plot = NA,
                                  return_plot = NA,
                                  save_plot = NA,
                                  save_param =  list(),
                                  default_save_name = 'spatNetwDistributions') {

  ## histogram to show
  distribution = match.arg(distribution, choices = distribution)

  ## spatial network
  #spatial_network = gobject@spatial_network[[spatial_network_name]]
  spatial_network = select_spatialNetwork(gobject,name = spatial_network_name,return_network_Obj = FALSE)
  if(is.null(spatial_network)) {
    stop('spatial network ', spatial_network_name, ' was not found')
  }


  if(distribution == 'distance') {

    spatNetwDistributionsDistance(gobject = gobject,
                                  spatial_network_name = spatial_network_name,
                                  hist_bins = hist_bins,
                                  test_distance_limit =  test_distance_limit,
                                  ncol = ncol,
                                  show_plot = show_plot,
                                  return_plot = return_plot,
                                  save_plot = save_plot,
                                  save_param =  save_param,
                                  default_save_name = default_save_name)

  } else if(distribution == 'k_neighbors') {

    spatNetwDistributionsKneighbors(gobject = gobject,
                                    spatial_network_name = spatial_network_name,
                                    hist_bins = hist_bins,
                                    show_plot = show_plot,
                                    return_plot = return_plot,
                                    save_plot = save_plot,
                                    save_param =  save_param,
                                    default_save_name = default_save_name)

  }

}






#' @title convert_to_full_spatial_network
#' @name convert_to_full_spatial_network
#' @param reduced_spatial_network_DT reduced spatial network in data.table format
#' @keywords internal
#' @description convert to a full spatial network
convert_to_full_spatial_network =  function(reduced_spatial_network_DT) {

  # data.table variables
  distance = rank_int = NULL

  # find location coordinates
  coordinates = grep('sdim', colnames(reduced_spatial_network_DT), value = T)

  begin_coordinates = grep('begin', coordinates, value = T)
  new_begin_coordinates = gsub(x = begin_coordinates, pattern = '_begin', replacement = '')
  new_begin_coordinates = gsub(x = new_begin_coordinates, pattern = 'sdim', replacement = 'source_')

  end_coordinates = grep('end', coordinates, value = T)
  new_end_coordinates = gsub(x = end_coordinates, pattern = '_end', replacement = '')
  new_end_coordinates = gsub(x = new_end_coordinates, pattern = 'sdim', replacement = 'target_')

  # create normal source --> target
  part1 = data.table::copy(reduced_spatial_network_DT)
  part1 = part1[, c('from', 'to', begin_coordinates, end_coordinates, 'distance', 'weight'), with = F]
  colnames(part1) = c('source', 'target', new_begin_coordinates, new_end_coordinates, 'distance', 'weight')

  # revert order target (now source) --> source (now target)
  part2 = data.table::copy(reduced_spatial_network_DT[, c('to', 'from', end_coordinates, begin_coordinates, 'distance', 'weight'), with = F])
  colnames(part2) = c('source', 'target', new_begin_coordinates, new_end_coordinates, 'distance', 'weight')

  # combine and remove duplicates
  full_spatial_network_DT = rbind(part1, part2)
  full_spatial_network_DT = unique(full_spatial_network_DT)

  # create ranking of interactions
  data.table::setorder(full_spatial_network_DT, source, distance)
  full_spatial_network_DT[, rank_int := 1:.N, by = 'source']

  # create unified column
  full_spatial_network_DT = sort_combine_two_DT_columns(full_spatial_network_DT, 'source', 'target', 'rnk_src_trgt')

  return(full_spatial_network_DT)

}

#' @title convert_to_reduced_spatial_network
#' @name convert_to_reduced_spatial_network
#' @param full_spatial_network_DT full spatial network in data.table format
#' @keywords internal
#' @description convert to a reduced spatial network
convert_to_reduced_spatial_network =  function(full_spatial_network_DT) {


  # data.table variables
  rnk_src_trgt = NULL

  # remove duplicates
  reduced_spatial_network_DT = full_spatial_network_DT[!duplicated(rnk_src_trgt)]
  reduced_spatial_network_DT[, c('rank_int', 'rnk_src_trgt') := NULL] # don't make sense in a reduced network

  # convert to names for a reduced network
  source_coordinates = grep('source_', colnames(reduced_spatial_network_DT), value = T)
  new_source_coordinates = gsub(x = source_coordinates, pattern = 'source_', replacement = 'sdim')
  new_source_coordinates = paste0(new_source_coordinates,'_begin')

  target_coordinates = grep('target_', colnames(reduced_spatial_network_DT), value = T)
  new_target_coordinates = gsub(x = target_coordinates, pattern = 'target_', replacement = 'sdim')
  new_target_coordinates = paste0(new_target_coordinates,'_end')

  reduced_spatial_network_DT = reduced_spatial_network_DT[, c('source', 'target', source_coordinates, target_coordinates, 'distance', 'weight'), with = F]
  colnames(reduced_spatial_network_DT) = c('from', 'to', new_source_coordinates, new_target_coordinates, 'distance', 'weight')
  return(reduced_spatial_network_DT)

}


#' @title create_spatialNetworkObject
#' @name create_spatialNetworkObject
#' @param name name
#' @param method method
#' @param parameters parameters
#' @param outputObj outputObj
#' @param networkDT networkDT
#' @param cellShapeObj cellShapeObj
#' @param networkDT_before_filter networkDT_before_filter
#' @param crossSectionObjects crossSectionObjects
#' @keywords internal
#' @description creates a spatial network object to store the created spatial network and additional information
create_spatialNetworkObject <- function(name = NULL,
                                        method = NULL,
                                        parameters = NULL,
                                        outputObj = NULL,
                                        networkDT = NULL,
                                        cellShapeObj = NULL,
                                        networkDT_before_filter = NULL,
                                        crossSectionObjects = NULL,
                                        misc = NULL) {

  networkObj = list(name = name,
                    method = method,
                    parameters = parameters,
                    outputObj = outputObj,
                    networkDT = networkDT,
                    networkDT_before_filter = networkDT_before_filter,
                    cellShapeObj = cellShapeObj,
                    crossSectionObjects = crossSectionObjects,
                    misc = misc)

  class(networkObj) <- append(class(networkObj), "spatialNetworkObj")
  return(networkObj)

}

#' @title select_spatialNetwork
#' @name select_spatialNetwork
#' @description function to select a spatial network
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

#' @title calculate_distance_and_weight
#' @name calculate_distance_and_weight
#' @param networkDT spatial network as data.table
#' @param sdimx spatial dimension x
#' @param sdimy spatial dimension y
#' @param sdimz spatial dimension z
#' @param d2_or_d3 number of dimensions
#' @description calculate_distance_and_weight
#' @keywords internal
calculate_distance_and_weight <- function(networkDT = NULL,
                                          sdimx = "sdimx",
                                          sdimy = "sdimy",
                                          sdimz = "sdimz",
                                          d2_or_d3=c(2,3)){

  # data.table variables
  distance = weight = from = NULL

  if(is.null(networkDT)) {
    stop('parameter networkDT can not be NULL \n')
  }

  # number of spatial dimensions TODO: chech with Huipeng!
  # d2_or_d3 = match.arg(d2_or_d3, choices = c(2,3))

  if (d2_or_d3==3){
    ## make it dynamic for all possible coordinates combinations ##
    xbegin_name = paste0(sdimx,'_begin')
    ybegin_name = paste0(sdimy,'_begin')
    zbegin_name =  paste0(sdimz,'_begin')
    xend_name = paste0(sdimx,'_end')
    yend_name = paste0(sdimy,'_end')
    zend_name = paste0(sdimz,'_end')
    mycols = c(xbegin_name, ybegin_name, zbegin_name,
               xend_name, yend_name, zend_name)
  }else if (d2_or_d3==2){
    xbegin_name = paste0(sdimx,'_begin')
    ybegin_name = paste0(sdimy,'_begin')
    xend_name = paste0(sdimx,'_end')
    yend_name = paste0(sdimy,'_end')
    mycols = c(xbegin_name, ybegin_name,
               xend_name, yend_name)
  }

  ## calculate distance and weight + filter ##
  networkDT[, `:=`(distance, stats::dist(x = matrix(.SD, nrow = 2, byrow = T))),
            by = 1:nrow(networkDT), .SDcols = mycols]

  networkDT[, `:=`(distance, as.numeric(distance))]
  networkDT[, `:=`(weight, 1/distance)]
  data.table::setorder(networkDT, from, distance)

  networkDT = networkDT[, c('to', 'from', 'weight',
                            'distance', mycols), with = F]

  return(networkDT)
}

#' @title filter_network
#' @name filter_network
#' @description function to filter a spatial network
#' @param networkDT spatial network in data.table format
#' @param maximum_distance maximum distance between cell centroids
#' @param minimum_k minimum number of neighbors
#' @keywords internal
filter_network <- function(networkDT = NULL,
                           maximum_distance = NULL,
                           minimum_k = NULL){

  # data.table variables
  distance = rank_int = NULL

  temp_fullnetwork = convert_to_full_spatial_network(networkDT)

  ## filter based on distance or minimum number of neighbors
  if (maximum_distance == "auto") {
    temp_fullnetwork = temp_fullnetwork[distance <= grDevices::boxplot.stats(temp_fullnetwork$distance)$stats[5] | rank_int <= minimum_k]
  }
  else if (!is.null(maximum_distance)) {
    temp_fullnetwork = temp_fullnetwork[distance <= maximum_distance | rank_int <= minimum_k]
  }
  networkDT = convert_to_reduced_spatial_network(temp_fullnetwork)

  return(networkDT)
}

## Delaunay network ####

#' @title create_delaunayNetwork_geometry
#' @description Create a spatial Delaunay network.
#' @keywords internal
create_delaunayNetwork_geometry <- function(spatial_locations,
                                               sdimx = 'sdimx',
                                               sdimy = 'sdimy',
                                               options = "Pp",
                                               ...) {


  # verify if optional package is installed
  package_check(pkg_name = "geometry", repository = "CRAN")

  # data.table variables
  from = to = NULL

  ## vector with original cell names ##
  cell_ID_vec = spatial_locations$cell_ID
  names(cell_ID_vec) = c(1:nrow(spatial_locations))

  ## create delaunay network
  delaunay_triangle = geometry::delaunayn(p = spatial_locations[, c(sdimx, sdimy), with = F],
                                          options = options, ...)

  ## save delaunay network object
  geometry_obj = list("delaunay_triangle" = delaunay_triangle)

  ## prepare delaunay network data.table results
  delaunay_edges <- as.data.table(rbind(delaunay_triangle[ ,c(1,2)],
                                        delaunay_triangle[ ,c(1,3)],
                                        delaunay_triangle[ ,c(2,3)]))

  delaunay_edges_dedup = unique(delaunay_edges)
  igraph_obj = igraph::graph_from_edgelist(as.matrix(delaunay_edges_dedup))
  adj_obj = igraph::as_adjacency_matrix(igraph_obj)
  igraph_obj2 = igraph::graph.adjacency(adj_obj)
  delaunay_edges_dedup2 = igraph::get.data.frame(igraph_obj2)
  delaunay_edges_dedup = data.table::as.data.table(delaunay_edges_dedup2)


  xbegin_name = paste0(sdimx,'_begin')
  ybegin_name = paste0(sdimy,'_begin')
  xend_name = paste0(sdimx,'_end')
  yend_name = paste0(sdimy,'_end')

  delaunay_network_DT = data.table::data.table(from = cell_ID_vec[delaunay_edges_dedup$from],
                                               to = cell_ID_vec[delaunay_edges_dedup$to],
                                               xbegin_name = spatial_locations[delaunay_edges_dedup$from, sdimx],
                                               ybegin_name = spatial_locations[delaunay_edges_dedup$from, sdimy],
                                               xend_name = spatial_locations[delaunay_edges_dedup$to, sdimx],
                                               yend_name = spatial_locations[delaunay_edges_dedup$to, sdimy])
  data.table::setnames(delaunay_network_DT,
                       old = c('xbegin_name', 'ybegin_name', 'xend_name', 'yend_name'),
                       new = c(xbegin_name, ybegin_name, xend_name, yend_name))
  data.table::setorder(delaunay_network_DT, from, to)

  out_object = list("geometry_obj" = geometry_obj,
                    "delaunay_network_DT" = delaunay_network_DT)

  return(out_object)
}

#' @title create_delaunayNetwork_geometry_3D
#' @description Create a spatial 3D Delaunay network with geometry
#' @keywords internal
create_delaunayNetwork_geometry_3D <- function(spatial_locations,
                                                  sdimx = 'sdimx',
                                                  sdimy = 'sdimy',
                                                  sdimz = 'sdimz',
                                                  options = options,
                                                  ...){


  # verify if optional package is installed
  package_check(pkg_name = "geometry", repository = "CRAN")


  # data.table variables
  from = to = NULL

  ## vector with original cell names ##
  cell_ID_vec = spatial_locations$cell_ID
  names(cell_ID_vec) = c(1:nrow(spatial_locations))


  delaunay_tetrahedra <- geometry::delaunayn(p = spatial_locations[, c(sdimx, sdimy, sdimz), with = F],
                                             options = options, ...)

  geometry_obj = list("delaunay_tetrahedra"=delaunay_tetrahedra)
  delaunay_edges <- as.data.table(rbind(delaunay_tetrahedra[,c(1,2)],
                                        delaunay_tetrahedra[,c(1,3)],
                                        delaunay_tetrahedra[,c(1,4)],
                                        delaunay_tetrahedra[,c(2,3)],
                                        delaunay_tetrahedra[,c(2,4)],
                                        delaunay_tetrahedra[,c(3,4)]))


  ### making sure of no duplication ###
  delaunay_edges_dedup = unique(delaunay_edges)
  igraph_obj = igraph::graph_from_edgelist(as.matrix(delaunay_edges_dedup))
  adj_obj = igraph::as_adjacency_matrix(igraph_obj)
  igraph_obj2 = igraph::graph.adjacency(adj_obj)
  delaunay_edges_dedup2 = igraph::get.data.frame(igraph_obj2)
  delaunay_edges_dedup = data.table::as.data.table(delaunay_edges_dedup2)
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

  xbegin_name = paste0(sdimx,'_begin')
  ybegin_name = paste0(sdimy,'_begin')
  zbegin_name = paste0(sdimz,'_begin')
  xend_name = paste0(sdimx,'_end')
  yend_name = paste0(sdimy,'_end')
  zend_name = paste0(sdimz,'_end')

  delaunay_network_DT = data.table::data.table(from = cell_ID_vec[delaunay_edges_dedup$from],
                                               to = cell_ID_vec[delaunay_edges_dedup$to],
                                               xbegin_name = spatial_locations[delaunay_edges_dedup$from, sdimx],
                                               ybegin_name = spatial_locations[delaunay_edges_dedup$from, sdimy],
                                               zbegin_name = spatial_locations[delaunay_edges_dedup$from, sdimz],
                                               xend_name = spatial_locations[delaunay_edges_dedup$to, sdimx],
                                               yend_name = spatial_locations[delaunay_edges_dedup$to, sdimy],
                                               zend_name = spatial_locations[delaunay_edges_dedup$to, sdimz])

  data.table::setnames(delaunay_network_DT,
                       old = c('xbegin_name', 'ybegin_name', 'zbegin_name', 'xend_name', 'yend_name', 'zend_name'),
                       new = c(xbegin_name, ybegin_name, zbegin_name, xend_name, yend_name, zend_name))
  data.table::setorder(delaunay_network_DT, from, to)

  out_object = list("geometry_obj"=geometry_obj,
                    "delaunay_network_DT"=delaunay_network_DT)
  return(out_object)

}

#' @title create_delaunayNetwork_RTriangle
#' @description Create a spatial Delaunay network with RTriangle
#' @keywords internal
create_delaunayNetwork_RTriangle <- function(spatial_locations,
                                                sdimx = 'sdimx',
                                                sdimy = 'sdimy',
                                                Y=TRUE,
                                                j=TRUE,
                                                S=0,
                                                ...){


  # verify if optional package is installed
  package_check(pkg_name = "RTriangle", repository = "CRAN")

  # data.table variables
  from = to = NULL

  ## vector with original cell names ##
  cell_ID_vec = spatial_locations$cell_ID
  names(cell_ID_vec) = c(1:nrow(spatial_locations))

  spatial_matrix = as.matrix(spatial_locations[, c(sdimx, sdimy), with = F])
  RTriangle_obj = RTriangle::triangulate(RTriangle::pslg(spatial_matrix),
                                         Y = Y,
                                         j = j,
                                         S = S,
                                         ...)


  ## prepare delaunay network data.table results
  xbegin_name = paste0(sdimx,'_begin')
  ybegin_name = paste0(sdimy,'_begin')
  xend_name = paste0(sdimx,'_end')
  yend_name = paste0(sdimy,'_end')

  delaunay_network_DT = data.table::data.table(from = cell_ID_vec[RTriangle_obj$E[,1]],
                                               to = cell_ID_vec[RTriangle_obj$E[, 2]],
                                               xbegin_name = RTriangle_obj$P[RTriangle_obj$E[, 1],1],
                                               ybegin_name = RTriangle_obj$P[RTriangle_obj$E[, 1], 2],
                                               xend_name = RTriangle_obj$P[RTriangle_obj$E[, 2], 1],
                                               yend_name = RTriangle_obj$P[RTriangle_obj$E[, 2], 2])

  data.table::setnames(delaunay_network_DT,
                       old = c('xbegin_name', 'ybegin_name', 'xend_name', 'yend_name'),
                       new = c(xbegin_name, ybegin_name, xend_name, yend_name))
  data.table::setorder(delaunay_network_DT, from, to)

  out_object = list("RTriangle_obj"=RTriangle_obj,
                    "delaunay_network_DT"=delaunay_network_DT)
  return(out_object)
}


#' @title create_delaunayNetwork_deldir
#' @description Create a spatial Delaunay network with deldir
#' @keywords internal
create_delaunayNetwork_deldir <- function(spatial_locations,
                                             sdimx = 'sdimx',
                                             sdimy = 'sdimy',
                                             ...){


  # data.table variables
  from = to = NULL

  ## vector with original cell names ##
  cell_ID_vec = spatial_locations$cell_ID
  names(cell_ID_vec) = c(1:nrow(spatial_locations))


  deldir_obj = deldir::deldir(x = spatial_locations[[sdimx]],
                              y = spatial_locations[[sdimy]],
                              ...)


  ## prepare delaunay network data.table results
  xbegin_name = paste0(sdimx,'_begin')
  ybegin_name = paste0(sdimy,'_begin')
  xend_name = paste0(sdimx,'_end')
  yend_name = paste0(sdimy,'_end')

  delaunay_network_DT = data.table::data.table(from = cell_ID_vec[deldir_obj$delsgs$ind1],
                                               to = cell_ID_vec[deldir_obj$delsgs$ind2],
                                               xbegin_name = deldir_obj$delsgs$x1,
                                               ybegin_name = deldir_obj$delsgs$y1,
                                               xend_name = deldir_obj$delsgs$x2,
                                               yend_name = deldir_obj$delsgs$y2)

  data.table::setnames(delaunay_network_DT,
                       old = c('xbegin_name', 'ybegin_name', 'xend_name', 'yend_name'),
                       new = c(xbegin_name, ybegin_name, xend_name, yend_name))
  data.table::setorder(delaunay_network_DT, from, to)

  out_object = list("deldir_obj"=deldir_obj,
                    "delaunay_network_DT"=delaunay_network_DT)
  return(out_object)
}








#' @title create_delaunayNetwork2D
#' @description Create a spatial 2D Delaunay network.
#' @keywords internal
create_delaunayNetwork2D <- function (gobject,
                                         method = c("delaunayn_geometry", "RTriangle", "deldir"),
                                         sdimx = 'sdimx',
                                         sdimy = 'sdimy',
                                         name = "delaunay_network",
                                         maximum_distance = "auto", # all
                                         minimum_k = 0, # all
                                         options = "Pp", # geometry
                                         Y = TRUE, # RTriange
                                         j = TRUE, # RTriange
                                         S = 0, # RTriange
                                         verbose = T,
                                         return_gobject = TRUE,
                                         ...)
{


  # get parameter values
  method = match.arg(method, c("delaunayn_geometry", "RTriangle", "deldir"))

  spatial_locations = gobject@spatial_locs
  spatial_locations = spatial_locations[, c('cell_ID', sdimx, sdimy), with = F]

  # 1. default is all dimensions as presented by spatial locations
  # 2. otherwise try to grab spatial coordinates
  # 3. stop if final result is not two columns


  if (method == "RTriangle"){

    delaunay_output = create_delaunayNetwork_RTriangle(spatial_locations = spatial_locations,
                                                          sdimx = sdimx,
                                                          sdimy = sdimy,
                                                          Y = Y,
                                                          j = j,
                                                          S = S,
                                                          ...)

    outputObj = delaunay_output$geometry_obj
    delaunay_network_DT = delaunay_output$delaunay_network_DT

    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
    parameters = list("maximum_distance" = maximum_distance,
                      "minimum_k" = minimum_k,
                      "Y" = Y,
                      "j" = j,
                      "S" = S)
    outputObj = outputObj
    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###


  }else if (method == "deldir"){

    delaunay_output = create_delaunayNetwork_deldir(spatial_locations = spatial_locations,
                                                       sdimx = sdimx,
                                                       sdimy = sdimy,
                                                       ...)

    outputObj = delaunay_output$geometry_obj
    delaunay_network_DT = delaunay_output$delaunay_network_DT

    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
    parameters = list("maximum_distance" = maximum_distance,
                      "minimum_k" = minimum_k)
    outputObj = outputObj
    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

  } else if (method == "delaunayn_geometry"){

    delaunay_output = create_delaunayNetwork_geometry(spatial_locations = spatial_locations,
                                                         sdimx = sdimx,
                                                         sdimy = sdimy,
                                                         options = options,
                                                         ...)

    outputObj = delaunay_output$geometry_obj
    delaunay_network_DT = delaunay_output$delaunay_network_DT

    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
    parameters = list("options" = options)
    outputObj = outputObj
    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###


  }


  ## calculate distance and weight + filter ##
  delaunay_network_DT = calculate_distance_and_weight(delaunay_network_DT,
                                                      sdimx = sdimx,
                                                      sdimy = sdimy,
                                                      d2_or_d3=2)
  networkDT_before_filter = delaunay_network_DT
  delaunay_network_DT = filter_network(delaunay_network_DT,
                                       maximum_distance = maximum_distance,
                                       minimum_k = minimum_k)

  ## calculate cell shape parameters ##
  meanCellDistance = get_distance(delaunay_network_DT,method="mean")
  medianCellDistance = get_distance(delaunay_network_DT,method="median")

  cellShapeObj = list("meanCellDistance" = meanCellDistance,
                      "medianCellDistance" = medianCellDistance)

  ###
  ###
  delaunay_network_Obj = create_spatialNetworkObject(name = name,
                                                              method = method,
                                                              parameters = parameters,
                                                              outputObj = outputObj,
                                                              networkDT = delaunay_network_DT,
                                                              networkDT_before_filter = networkDT_before_filter,
                                                              cellShapeObj = cellShapeObj,
                                                              misc = NULL)
  ###
  ###

  if (return_gobject == TRUE) {

    spn_names = names(gobject@spatial_network)
    if (name %in% spn_names) {
      cat("\n ", name, " has already been used, will be overwritten \n")
    }
    parameters_list = gobject@parameters
    number_of_rounds = length(parameters_list)
    update_name = paste0(number_of_rounds, "_delaunay_spatial_network")

    if(method == "delaunayn_geometry") {
      parameters_list[[update_name]] = c(`dimensions used` = paste0('dimensions: ', sdimx, ' and ', sdimy),
                                         `method` = method,
                                         `maximum distance threshold` = ifelse(is.null(maximum_distance),  NA, maximum_distance),
                                         `name of spatial network` = name)
    } else if(method == "RTriangle") {
      parameters_list[[update_name]] = c(`dimensions used` = paste0('dimensions: ', sdimx, ' and ', sdimy),
                                         `method` = method,
                                         `maximum distance threshold` = ifelse(is.null(maximum_distance),  NA, maximum_distance),
                                         `RTriangle Y:` = Y,
                                         `RTriangle j:` = j,
                                         `RTriangle S:` = S,
                                         `name of spatial network` = name)
    } else if(method == 'deldir') {
      parameters_list[[update_name]] = c(`dimensions used` = paste0('dimensions: ', sdimx, ' and ', sdimy),
                                         `method` = method,
                                         `maximum distance threshold` = ifelse(is.null(maximum_distance),  NA, maximum_distance),
                                         `name of spatial network` = name)
    }

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




#' @title create_delaunayNetwork3D
#' @description Create a spatial 3D Delaunay network.
#' @keywords internal
create_delaunayNetwork3D <- function (gobject,
                                      method = "delaunayn_geometry",
                                      sdimx = 'sdimx',
                                      sdimy = 'sdimy',
                                      sdimz = 'sdimz',
                                      name = "delaunay_network_3D",
                                      maximum_distance = "auto",
                                      minimum_k = 0, # all
                                      options = "Pp", # geometry
                                      return_gobject = TRUE,
                                      ...)
{

  # get parameter values
  method = match.arg(method, c("delaunayn_geometry", "RTriangle", "deldir"))

  spatial_locations = gobject@spatial_locs
  spatial_locations = spatial_locations[, c('cell_ID', sdimx, sdimy, sdimz), with = F]


  ## delaunay geometry method ##
  if (method == "delaunayn_geometry"){

    delaunay_output = create_delaunayNetwork_geometry_3D(spatial_locations = spatial_locations,
                                                            sdimx = sdimx,
                                                            sdimy = sdimy,
                                                            sdimz = sdimz,
                                                            options = options,
                                                            ...)

    outputObj = delaunay_output$geometry_obj
    delaunay_network_DT = delaunay_output$delaunay_network_DT

    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
    parameters = list("options" = options)
    outputObj = outputObj
    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

  }

  ## calculate distance and weight + filter ##
  networkDT_before_filter = calculate_distance_and_weight(delaunay_network_DT,
                                                          sdimx = sdimx,
                                                          sdimy = sdimy,
                                                          sdimz = sdimz,
                                                          d2_or_d3=3)
  delaunay_network_DT = filter_network(networkDT_before_filter,
                                       maximum_distance=maximum_distance,
                                       minimum_k=minimum_k)

  ## calculate cell shape parameters ##
  meanCellDistance = get_distance(delaunay_network_DT,method="mean")
  medianCellDistance = get_distance(delaunay_network_DT,method="median")

  cellShapeObj = list("meanCellDistance" = meanCellDistance,
                      "medianCellDistance" = medianCellDistance
  )

  if (return_gobject == TRUE) {
    spn_names = names(gobject@spatial_network)
    if (name %in% spn_names) {
      cat("\n ", name, " has already been used, will be overwritten \n")
    }
    parameters_list = gobject@parameters
    number_of_rounds = length(parameters_list)
    update_name = paste0(number_of_rounds, "_delaunay_spatial_network_3D")

    parameters_list[[update_name]] = c(`dimensions used` = paste0('dimensions: ', sdimx, ', ', sdimy, ' and ', sdimz),
                                       `method` = method,
                                       `maximum distance threshold` = ifelse(is.null(maximum_distance),  NA, maximum_distance),
                                       `minimum k` = minimum_k,
                                       `name of spatial network` = name)

    gobject@parameters = parameters_list

    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
    ###
    ###
    delaunay_network_Obj = create_spatialNetworkObject(name = name,
                                                       method = method,
                                                       parameters = parameters,
                                                       outputObj = outputObj,
                                                       networkDT = delaunay_network_DT,
                                                       networkDT_before_filter = networkDT_before_filter,
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








#' @title createSpatialDelaunayNetwork
#' @description Create a spatial Delaunay network based on cell centroid physical distances.
#' @param gobject giotto object
#' @param method package to use to create a Delaunay network
#' @param dimensions which spatial dimensions to use. Use "sdimx" (spatial dimension x), "sdimy", "sdimz" respectively to refer to X (or the 1st), Y (or the 2nd) and Z(or the 3rd) dimension, see details. (default = all)
#' @param name name for spatial network (default = 'delaunay_network')
#' @param maximum_distance distance cuttof for Delaunay neighbors to consider. If "auto", "upper wisker" value of the distance vector between neighbors is used; see the boxplot{graphics} documentation for more details.(default = "auto")
#' @param minimum_k minimum number of neigbhours if maximum_distance != NULL
#' @param options (geometry) String containing extra control options for the underlying Qhull command; see the Qhull documentation (../doc/qhull/html/qdelaun.html) for the available options. (default = 'Pp', do not report precision problems)
#' @param Y (RTriangle) If TRUE prohibits the insertion of Steiner points on the mesh boundary.
#' @param j (RTriangle) If TRUE jettisons vertices that are not part of the final triangulation from the output.
#' @param S (RTriangle) Specifies the maximum number of added Steiner points.
#' @param verbose verbose
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @param \dots Other additional parameters
#' @return giotto object with updated spatial network slot
#' @details Creates a spatial Delaunay network as explained in \code{\link[geometry]{delaunayn}} (default), \code{\link[deldir]{deldir}}, or \code{\link[RTriangle]{triangulate}}.
#' @export
createSpatialDelaunayNetwork <- function(gobject,
                                         method = c("deldir", "delaunayn_geometry", "RTriangle"),
                                         dimensions = "all",
                                         name = "Delaunay_network",
                                         maximum_distance = "auto", # all
                                         minimum_k = 0, # all
                                         options = "Pp", # geometry
                                         Y = TRUE, # RTriange
                                         j = TRUE, # RTriange
                                         S = 0, # RTriange
                                         verbose = T,
                                         return_gobject = TRUE,
                                         ...) {


  # get parameter values
  method = match.arg(method, c("deldir", "delaunayn_geometry", "RTriangle"))

  # determine the network dimesions
  spatial_locations = gobject@spatial_locs
  spatial_locations = spatial_locations[, grepl("sdim", colnames(spatial_locations)),  with = F]

  if (dimensions != "all") {
    spatial_locations = spatial_locations[, dimensions, with = FALSE]
  }
  spatial_locations <- as.matrix(spatial_locations)
  d2_or_d3 = dim(spatial_locations)[2]


  # create 2D or 3D delaunay network
  if (d2_or_d3 == 2){

    first_dimension = colnames(spatial_locations)[[1]]
    second_dimension = colnames(spatial_locations)[[2]]

    out = create_delaunayNetwork2D(gobject=gobject,
                                      method = method,
                                      sdimx = first_dimension,
                                      sdimy = second_dimension,
                                      name = name,
                                      maximum_distance = maximum_distance,
                                      minimum_k = minimum_k,
                                      options = options,
                                      Y = Y,
                                      j = j,
                                      S = S,
                                      verbose = verbose,
                                      return_gobject = return_gobject,
                                      ...)
  }else if(d2_or_d3 == 3){

    if (method!="delaunayn_geometry"){
      stop(method, ' method only applies to 2D data, use delaunayn_geometry, see details \n')
    }else{

      first_dimension = colnames(spatial_locations)[[1]]
      second_dimension = colnames(spatial_locations)[[2]]
      third_dimension = colnames(spatial_locations)[[3]]

      out = create_delaunayNetwork3D(gobject=gobject,
                                        method = method,
                                        sdimx = first_dimension,
                                        sdimy = second_dimension,
                                        sdimz = third_dimension,
                                        name = name,
                                        maximum_distance = maximum_distance,
                                        minimum_k = minimum_k,
                                        options = options,
                                        return_gobject = return_gobject,
                                        ...)
    }
  }

  return(out)

}




#' @title plotStatDelaunayNetwork
#' @description Plots network statistics for a Delaunay network..
#' @param gobject giotto object
#' @param method package to use to create a Delaunay network
#' @param dimensions which spatial dimensions to use (maximum 2 dimensions)
#' @param maximum_distance distance cuttof for Delaunay neighbors to consider
#' @param minimum_k minimum neigbhours if maximum_distance != NULL
#' @param options (geometry) String containing extra control options for the underlying Qhull command; see the Qhull documentation (../doc/qhull/html/qdelaun.html) for the available options. (default = 'Pp', do not report precision problems)
#' @param Y (RTriangle) If TRUE prohibits the insertion of Steiner points on the mesh boundary.
#' @param j (RTriangle) If TRUE jettisons vertices that are not part of the final triangulation from the output.
#' @param S (RTriangle) Specifies the maximum number of added Steiner points.
#' @param show_plot show plots
#' @param return_plot return ggplot object
#' @param save_plot directly save the plot [boolean]
#' @param save_param list of saving parameters, see \code{\link{showSaveParameters}}
#' @param default_save_name default save name for saving, don't change, change save_name in save_param
#' @param \dots Other parameters
#' @return giotto object with updated spatial network slot
#' @export
plotStatDelaunayNetwork = function(gobject,
                                   method = c("deldir", "delaunayn_geometry", "RTriangle"),
                                   dimensions = "all",
                                   maximum_distance = "auto", # all
                                   minimum_k = 0, # all
                                   options = "Pp", # geometry
                                   Y = TRUE, # RTriange
                                   j = TRUE, # RTriange
                                   S = 0, # RTriange
                                   show_plot = NA,
                                   return_plot = NA,
                                   save_plot = NA,
                                   save_param =  list(),
                                   default_save_name = 'plotStatDelaunayNetwork',
                                   ...) {


  # data.table variables
  distance = rank_int = N = NULL

  delaunay_network_DT = createSpatialDelaunayNetwork(gobject = gobject,
                                                        method = method,
                                                        dimensions = dimensions,
                                                        name = 'temp_network',
                                                        maximum_distance = maximum_distance, # all
                                                        minimum_k = minimum_k, # all
                                                        options = options, # geometry
                                                        Y = Y, # RTriange
                                                        j = j, # RTriange
                                                        S = S, # RTriange
                                                        return_gobject = F,
                                                        ...)

  delaunay_network_DT_c = convert_to_full_spatial_network(reduced_spatial_network_DT = delaunay_network_DT)


  ## create visuals
  pl1 = ggplot(delaunay_network_DT, aes(x=factor(""), y=distance))
  pl1 = pl1 + theme_classic() + theme(plot.title = element_text(hjust=0.5))
  pl1 = pl1 + geom_boxplot(outlier.colour = "red", outlier.shape = 1)
  pl1 = pl1 + labs(title = 'Delaunay network', y = 'cell-cell distances', x = '')

  pl2 = ggplot(delaunay_network_DT_c, aes(x=factor(rank_int), y=distance))
  pl2 = pl2 + theme_classic() + theme(plot.title = element_text(hjust=0.5))
  pl2 = pl2 + geom_boxplot(outlier.colour = "red", outlier.shape = 1)
  pl2 = pl2 + labs(title = 'Delaunay network by neigbor ranking', y = 'cell-cell distances', x = '')

  neighbors = delaunay_network_DT_c[, .N, by = source]
  pl3 = ggplot()
  pl3 = pl3 + theme_classic() + theme(plot.title = element_text(hjust=0.5))
  pl3 = pl3 + geom_histogram(data = neighbors, aes(x = as.factor(N)), stat = 'count')
  pl3 = pl3 + labs(title = 'Delaunay network neigbors per cell', y = 'count', x = '')
  pl3

  savelist = list(pl1, pl2, pl3)



  ## combine plots with cowplot
  combo_plot <- cowplot::plot_grid(pl1, pl2, NULL, pl3,
                                   ncol = 2,
                                   rel_heights = c(1, 1), rel_widths = c(1, 2), align = 'v')


  ## print, return and save parameters
  show_plot = ifelse(is.na(show_plot), readGiottoInstructions(gobject, param = 'show_plot'), show_plot)
  save_plot = ifelse(is.na(save_plot), readGiottoInstructions(gobject, param = 'save_plot'), save_plot)
  return_plot = ifelse(is.na(return_plot), readGiottoInstructions(gobject, param = 'return_plot'), return_plot)

  ## print plot
  if(show_plot == TRUE) {
    print(combo_plot)
  }

  ## save plot
  if(save_plot == TRUE) {
    do.call('all_plots_save_function', c(list(gobject = gobject, plot_object = combo_plot, default_save_name = default_save_name), save_param))
  }

  ## return plot
  if(return_plot == TRUE) {
    return(combo_plot)
  }


}




## kNN network ####

#' @title create_KNNnetwork_dbscan
#' @description Create a spatial knn network with dbscan
#' @keywords internal
create_KNNnetwork_dbscan = function(spatial_locations,
                                    sdimx = 'sdimx',
                                    sdimy = 'sdimy',
                                    sdimz = 'sdimz',
                                    k = 4,
                                    ...) {

  # data.table variables
  from = to = NULL

  ## vector with original cell names ##
  cell_ID_vec = spatial_locations$cell_ID
  names(cell_ID_vec) = c(1:nrow(spatial_locations))


  ## set dimension coordinates to NULL if they don't exist
  if(!is.null(sdimx)) {
    if(sdimx %in% colnames(spatial_locations)) {
      sdimx = sdimx
    } else {
      sdimx = NULL
    }
  }

  if(!is.null(sdimy)) {
    if(sdimy %in% colnames(spatial_locations)) {
      sdimy = sdimy
    } else {
      sdimy = NULL
    }
  }

  if(!is.null(sdimz)) {
    if(sdimz %in% colnames(spatial_locations)) {
      sdimz = sdimz
    } else {
      sdimz = NULL
    }
  }


  ## create knn network
  spatial_locations_matrix = as.matrix(spatial_locations[, c(sdimx, sdimy, sdimz), with = F])

  knn_spatial <- dbscan::kNN(x = spatial_locations_matrix,
                             k = k,
                             ...)

  knn_sptial.norm = data.frame(from = rep(1:nrow(knn_spatial$id), k),
                               to = as.vector(knn_spatial$id),
                               weight = 1/(1 + as.vector(knn_spatial$dist)),
                               distance = as.vector(knn_spatial$dist))
  nw_sptial.norm = igraph::graph_from_data_frame(knn_sptial.norm, directed = FALSE)
  network_DT = data.table::as.data.table(knn_sptial.norm)


  #spatial_network_DT[, `:=`(from, cell_ID_vec[from])]
  #spatial_network_DT[, `:=`(to, cell_ID_vec[to])]


  xbegin_name = paste0(sdimx,'_begin')
  ybegin_name = paste0(sdimy,'_begin')
  zbegin_name = paste0(sdimz,'_begin')
  xend_name = paste0(sdimx,'_end')
  yend_name = paste0(sdimy,'_end')
  zend_name = paste0(sdimz,'_end')

  if(!is.null(sdimz)) {
    spatial_network_DT = data.table::data.table(from = cell_ID_vec[network_DT$from],
                                                to = cell_ID_vec[network_DT$to],
                                                xbegin_name = spatial_locations[network_DT$from, sdimx],
                                                ybegin_name = spatial_locations[network_DT$from, sdimy],
                                                zbegin_name = spatial_locations[network_DT$from, sdimz],
                                                xend_name = spatial_locations[network_DT$to, sdimx],
                                                yend_name = spatial_locations[network_DT$to, sdimy],
                                                zend_name = spatial_locations[network_DT$to, sdimz],
                                                distance = network_DT$distance,
                                                weight = network_DT$weight)

    data.table::setnames(spatial_network_DT,
                         old = c('xbegin_name', 'ybegin_name', 'zbegin_name', 'xend_name', 'yend_name', 'zend_name'),
                         new = c(xbegin_name, ybegin_name, zbegin_name, xend_name, yend_name, zend_name))
    data.table::setorder(spatial_network_DT, from, to)


  } else {
    spatial_network_DT = data.table::data.table(from = cell_ID_vec[network_DT$from],
                                                to = cell_ID_vec[network_DT$to],
                                                xbegin_name = spatial_locations[network_DT$from, sdimx],
                                                ybegin_name = spatial_locations[network_DT$from, sdimy],
                                                xend_name = spatial_locations[network_DT$to, sdimx],
                                                yend_name = spatial_locations[network_DT$to, sdimy],
                                                distance = network_DT$distance,
                                                weight = network_DT$weight)
    data.table::setnames(spatial_network_DT,
                         old = c('xbegin_name', 'ybegin_name', 'xend_name', 'yend_name'),
                         new = c(xbegin_name, ybegin_name, xend_name, yend_name))
    data.table::setorder(spatial_network_DT, from, to)

  }

  out_object = list("knn_obj" = knn_spatial,
                    "spatial_network_DT"= spatial_network_DT)
  return(out_object)

}



#' @title createSpatialKNNnetwork
#' @description Create a spatial knn network.
#' @param gobject giotto object
#' @param method method to create kNN network
#' @param dimensions which spatial dimensions to use (default = all)
#' @param name name for spatial network (default = 'spatial_network')
#' @param k number of nearest neighbors based on physical distance
#' @param maximum_distance distance cuttof for nearest neighbors to consider for kNN network
#' @param minimum_k minimum nearest neigbhours if maximum_distance != NULL
#' @param verbose verbose
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @param \dots additional arguments to the selected method function
#' @return giotto object with updated spatial network slot
#'
#' \strong{dimensions: } default = 'all' which takes all possible dimensions.
#' Alternatively you can provide a character vector that specififies the spatial dimensions to use, e.g. c("sdimx', "sdimy")
#' or a numerical vector, e.g. 2:3
#'
#' \strong{maximum_distance: } to create a network based on maximum distance only, you also need to set k to a very high value, e.g. k = 100
#'
#'
#' @export
createSpatialKNNnetwork <- function (gobject,
                                     method = "dbscan",
                                     dimensions = "all",
                                     name = "knn_network",
                                     k = 4,
                                     maximum_distance = NULL,
                                     minimum_k = 0,
                                     verbose = F,
                                     return_gobject = TRUE,
                                     ...)
{


  # data.table variables
  distance = rank_int = NULL

  # get parameter values
  method = match.arg(method, c("dbscan"))

  spatial_locations = gobject@spatial_locs


  if (dimensions != "all") {
    temp_spatial_locations = spatial_locations[, dimensions, with = FALSE]
  } else {
    temp_spatial_locations = spatial_locations[, grepl('sdim', colnames(spatial_locations)), with = FALSE]
  }
  temp_spatial_locations <- as.matrix(temp_spatial_locations)

  first_dimension = colnames(temp_spatial_locations)[[1]]
  second_dimension = colnames(temp_spatial_locations)[[2]]
  if(ncol(temp_spatial_locations) > 2) {
    third_dimension = colnames(temp_spatial_locations)[[3]]
  } else {
    third_dimension = NULL
  }

  if (method == "dbscan"){

    spatial_locations = spatial_locations[, c('cell_ID', first_dimension, second_dimension, third_dimension), with = F]


    knn_output = create_KNNnetwork_dbscan(spatial_locations = spatial_locations,
                                          k = k,
                                          sdimx = first_dimension,
                                          sdimy = second_dimension,
                                          sdimz = third_dimension,
                                          ...)

    outputObj = knn_output$knn_obj
    spatial_network_DT = knn_output$spatial_network_DT

    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
    parameters = list("neighbors" = k,
                      "maximum_distance" = maximum_distance,
                      "minimum_k" = minimum_k)
    outputObj = outputObj
    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###


  } else {

    stop('no other methods to create kNN spatial networks have been implemented')

  }


  temp_fullnetwork = convert_to_full_spatial_network(spatial_network_DT)
  if (!is.null(maximum_distance)) {
    temp_fullnetwork = temp_fullnetwork[distance <= maximum_distance | rank_int <= minimum_k]
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
    spn_names = names(gobject@spatial_network)
    if (name %in% spn_names) {
      cat("\n ", name, " has already been used, will be overwritten \n")
    }
    parameters_list = gobject@parameters
    number_of_rounds = length(parameters_list)
    update_name = paste0(number_of_rounds, "_spatial_network")

    parameters_list[[update_name]] = c(`k neighbours` = k,
                                       `dimensions used` = dimensions,
                                       `maximum distance threshold` = ifelse(is.null(maximum_distance), NA, maximum_distance),
                                       `name of spatial network` = name)
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









## spatial network ####

#' @title createSpatialNetwork
#' @description Create a spatial network based on cell centroid physical distances.
#' @param gobject giotto object
#' @param dimensions which spatial dimensions to use (default = all)
#' @param method which method to use to create a spatial network. (default = Delaunay)
#' @param delaunay_method Delaunay method to use
#' @param maximum_distance_delaunay distance cuttof for nearest neighbors to consider for Delaunay network
#' @param options (geometry) String containing extra control options for the underlying Qhull command; see the Qhull documentation (../doc/qhull/html/qdelaun.html) for the available options. (default = 'Pp', do not report precision problems)
#' @param Y (RTriangle) If TRUE prohibits the insertion of Steiner points on the mesh boundary.
#' @param j (RTriangle) If TRUE jettisons vertices that are not part of the final triangulation from the output.
#' @param S (RTriangle) Specifies the maximum number of added Steiner points.
#' @param name name for spatial network (default = 'spatial_network')
#' @param knn_method method to create kNN network
#' @param k number of nearest neighbors based on physical distance
#' @param minimum_k minimum nearest neigbhours if maximum_distance != NULL
#' @param maximum_distance_knn distance cuttof for nearest neighbors to consider for kNN network
#' @param verbose verbose
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @param \dots Additional parameters for the selected function
#' @return giotto object with updated spatial network slot
#' @details Creates a spatial network connecting single-cells based on their physical distance to each other.
#' For Delaunay method, neighbors will be decided by delaunay triangulation and a maximum distance criteria. For kNN method, number of neighbors can be determined by k, or maximum distance from each cell with or without
#' setting a minimum k for each cell.
#'
#' \strong{dimensions: } default = 'all' which takes all possible dimensions.
#' Alternatively you can provide a character vector that specififies the spatial dimensions to use, e.g. c("sdimx', "sdimy")
#' or a numerical vector, e.g. 2:3
#'
#' @export
createSpatialNetwork <- function(gobject,
                                 name = NULL,
                                 dimensions = "all",
                                 method = c('Delaunay', 'kNN'),
                                 delaunay_method = c("deldir", "delaunayn_geometry", "RTriangle"),
                                 maximum_distance_delaunay = "auto",
                                 options = "Pp",
                                 Y = TRUE,
                                 j = TRUE,
                                 S = 0,
                                 minimum_k = 0,
                                 knn_method = "dbscan",
                                 k = 4,
                                 maximum_distance_knn = NULL,
                                 verbose = F,
                                 return_gobject = TRUE,
                                 ...){

  # get paramters
  method = match.arg(method, c('Delaunay', 'kNN'))


  if(method=="kNN"){
    if(is.null(name)){
      name = paste0(method,"_","network")
    }

    knn_method = match.arg(knn_method,c("dbscan"))

    out = createSpatialKNNnetwork(gobject = gobject,
                                  method = knn_method,
                                  dimensions = dimensions,
                                  k = k,
                                  maximum_distance = maximum_distance_knn,
                                  minimum_k = minimum_k,
                                  name = name,
                                  verbose = verbose,
                                  return_gobject = return_gobject,
                                  ...)

  } else if (method=="Delaunay"){

    delaunay_method = match.arg(delaunay_method, c("deldir", "delaunayn_geometry", "RTriangle"))
    if(is.null(name)){
      name = paste0(method,"_","network")
    }
    out = createSpatialDelaunayNetwork(gobject=gobject,
                                       method = delaunay_method,
                                       dimensions = dimensions,
                                       name = name,
                                       maximum_distance = maximum_distance_delaunay,
                                       options = options,
                                       minimum_k = minimum_k,
                                       Y = Y,
                                       j = j,
                                       S = S,
                                       verbose = verbose,
                                       return_gobject = return_gobject,
                                       ...)
  }

  return(out)
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



#' @title annotateSpatialNetwork
#' @name annotateSpatialNetwork
#' @description Annotate spatial network with cell metadata information.
#' @param gobject giotto object
#' @param spatial_network_name name of spatial network to use
#' @param cluster_column name of column to use for clusters
#' @param create_full_network convert from reduced to full network representation
#' @return annotated network in data.table format
#' @export
annotateSpatialNetwork = function(gobject,
                                  spatial_network_name = 'Delaunay_network',
                                  cluster_column,
                                  create_full_network = F) {

  # get network
  if(!spatial_network_name %in% names(gobject@spatial_network)) {
    stop('\n spatial network with name: ', spatial_network_name, ' does not exist \n')
  }
  spatial_network = select_spatialNetwork(gobject,name = spatial_network_name,return_network_Obj = FALSE)



  if(create_full_network == TRUE) {

    spatial_network = convert_to_full_spatial_network(spatial_network)

    # convert to names for a reduced network
    source_coordinates = grep('source_', colnames(spatial_network), value = T)
    new_source_coordinates = gsub(x = source_coordinates, pattern = 'source_', replacement = 'sdim')
    new_source_coordinates = paste0(new_source_coordinates,'_begin')

    target_coordinates = grep('target_', colnames(spatial_network), value = T)
    new_target_coordinates = gsub(x = target_coordinates, pattern = 'target_', replacement = 'sdim')
    new_target_coordinates = paste0(new_target_coordinates,'_end')

    data.table::setnames(spatial_network,
                         old = c('source', 'target', source_coordinates, target_coordinates),
                         new = c('from', 'to', new_source_coordinates, new_target_coordinates))
  }



  # cell metadata
  cell_metadata = pDataDT(gobject)
  if(!cluster_column %in% colnames(cell_metadata)) {
    stop('\n the cluster column does not exist in pDataDT(gobject) \n')
  }
  cluster_type_vector = cell_metadata[[cluster_column]]
  names(cluster_type_vector) = cell_metadata[['cell_ID']]

  # data.table variables
  to_cell_type = to = from_cell_type = from = type_int = from_to = NULL

  spatial_network_annot = data.table::copy(spatial_network)
  spatial_network_annot[, to_cell_type := cluster_type_vector[to]]
  spatial_network_annot[, from_cell_type := cluster_type_vector[from]]
  spatial_network_annot[, type_int := ifelse(to_cell_type == from_cell_type, 'homo', 'hetero')]

  # specific direction
  spatial_network_annot[, from_to := paste0(from_cell_type,'-',to_cell_type)]

  # unified direction, due to 'sort'
  spatial_network_annot = sort_combine_two_DT_columns(spatial_network_annot,
                                                      column1 = 'from_cell_type',
                                                      column2 = 'to_cell_type',
                                                      myname = 'unified_int')

  return(spatial_network_annot)

}





## Spatial grid ####

#' @title find_grid_3D
#' @name find_grid_3D
#' @description find grid location in 3D
#' @keywords internal
find_grid_3D <- function(grid_DT, x_loc, y_loc, z_loc) {

  # data.table variables
  x_start = x_end = y_start = y_end = z_start = z_end = NULL

  name = grid_DT[x_loc > x_start & x_loc < x_end & y_loc > y_start & y_loc < y_end & z_loc > z_start & z_loc < z_end]$gr_name
  return(name)
}

#' @title find_grid_2D
#' @name find_grid_2D
#' @description find grid location in 2D
#' @keywords internal
find_grid_2D <- function(grid_DT, x_loc, y_loc) {

  # data.table variables
  x_start = x_end = y_start = y_end = NULL

  name = grid_DT[x_loc > x_start & x_loc < x_end & y_loc > y_start & y_loc < y_end]$gr_name
  return(name)
}

#' @title find_grid_x
#' @name find_grid_x
#' @description find grid location on x-axis
#' @keywords internal
find_grid_x <- function(grid_DT, x_loc) {

  # data.table variables
  x_start = x_end = gr_x_name = NULL

  grid_DT_x = unique(grid_DT[,.(x_start, x_end, gr_x_name)])
  name_x = grid_DT_x[x_loc > x_start & x_loc < x_end]$gr_x_name
  return(name_x)
}

#' @title find_grid_y
#' @name find_grid_y
#' @description find grid location on y-axis
#' @keywords internal
find_grid_y <- function(grid_DT, y_loc) {

  # data.table variables
  y_start = y_end = gr_y_name = NULL

  grid_DT_y = unique(grid_DT[,.(y_start, y_end, gr_y_name)])
  name_y = grid_DT_y[y_loc > y_start & y_loc < y_end]$gr_y_name
  return(name_y)
}

#' @title find_grid_z
#' @name find_grid_z
#' @description find grid location on z-axis
#' @keywords internal
find_grid_z <- function(grid_DT, z_loc) {

  # data.table variables
  z_start = z_end = gr_z_name = NULL

  grid_DT_z = unique(grid_DT[,.(z_start, z_end, gr_z_name)])
  name_z = grid_DT_z[z_loc > z_start & z_loc < z_end]$gr_z_name
  return(name_z)
}



#' @title create_spatialGridObject
#' @description create a spatial grid object
#' @keywords internal
create_spatialGridObject <- function(name = NULL,
                                     method = NULL,
                                     parameters = NULL,
                                     gridDT = NULL,
                                     outputObj = NULL,
                                     misc = NULL) {

  gridObj = list(name = name,
                 method = method,
                 parameters = parameters,
                 gridDT = gridDT,
                 misc = misc)

  class(gridObj) <- append(class(gridObj), "spatialGridObj")
  return(gridObj)
  }


#' @title create_spatialGrid_default_2D
#' @description create a 2D spatial grid
#' @keywords internal
create_spatialGrid_default_2D <- function(gobject,
                                          sdimx_stepsize = NULL,
                                          sdimy_stepsize = NULL,
                                          minimum_padding = 1) {


  # data.table variables
  gr_name = gr_x_name = gr_y_name = gr_x_loc = gr_y_loc = gr_loc = NULL

  spatlocs = data.table::copy(gobject@spatial_locs)

  if(is.null(spatlocs)) stop('\n spatial locations are needed to create a spatial grid \n')

  ## calculate sequences for desired stepsize
  # x-axis
  x_range = range(spatlocs$sdimx)
  x_start = x_range[[1]] - minimum_padding
  x_end = x_range[[2]] + minimum_padding
  dimx_steps = ceiling( (x_end-x_start) / sdimx_stepsize)
  dimx_start = mean(c(x_start, x_end))-((dimx_steps/2)*sdimx_stepsize)
  dimx_end = mean(c(x_start, x_end))+((dimx_steps/2)*sdimx_stepsize)
  my_x_seq = seq(from = dimx_start, to = dimx_end, by = sdimx_stepsize)

  # y-axis
  y_range = range(spatlocs$sdimy)
  y_start = y_range[[1]] - minimum_padding
  y_end = y_range[[2]] + minimum_padding
  dimy_steps = ceiling( (y_end-y_start) / sdimy_stepsize)
  dimy_start = mean(c(y_start, y_end))-((dimy_steps/2)*sdimy_stepsize)
  dimy_end = mean(c(y_start, y_end))+((dimy_steps/2)*sdimy_stepsize)
  my_y_seq = seq(from = dimy_start, to = dimy_end, by = sdimy_stepsize)


  ## create grid with starts and ends
  grid_starts = data.table::as.data.table(expand.grid(my_x_seq[-length(my_x_seq)],
                                                      my_y_seq[-length(my_y_seq)]))
  colnames(grid_starts) = c('x_start', 'y_start')
  grid_ends = data.table::as.data.table(expand.grid(my_x_seq[-1],
                                                    my_y_seq[-1]))
  colnames(grid_ends) = c('x_end', 'y_end')
  spatgrid = cbind(grid_starts, grid_ends)


  ## first label the grid itself ##
  spatgrid[, gr_name := paste0('gr_', 1:.N)]

  # x-axis
  x_labels = sort(unique(spatgrid$x_start))
  x_gr_names = paste0('gr_x_', 1:length(x_labels))
  names(x_gr_names) = x_labels
  x_gr_names_vector = x_gr_names[as.character(spatgrid$x_start)]
  spatgrid[, gr_x_name := x_gr_names_vector]

  # y-axis
  y_labels = sort(unique(spatgrid$y_start))
  y_gr_names = paste0('gr_y_', 1:length(y_labels))
  names(y_gr_names) = y_labels
  y_gr_names_vector = y_gr_names[as.character(spatgrid$y_start)]
  spatgrid[, gr_y_name := y_gr_names_vector]

  ## for all dimensions ##
  # converter
  gr_dim_names = spatgrid$gr_name
  names(gr_dim_names) = paste0(spatgrid$gr_x_name,'-', spatgrid$gr_y_name)


  return(spatgrid)

}


#' @title create_spatialGrid_default_3D
#' @description create a 3D spatial grid
#' @keywords internal
create_spatialGrid_default_3D <- function(gobject,
                                          sdimx_stepsize = NULL,
                                          sdimy_stepsize = NULL,
                                          sdimz_stepsize = NULL,
                                          minimum_padding = 1) {


  # data.table variables
  gr_name = gr_x_name = gr_y_name = gr_z_name = gr_x_loc = gr_y_loc = gr_z_loc = gr_loc = NULL

  spatlocs = data.table::copy(gobject@spatial_locs)

  if(is.null(spatlocs)) stop('\n spatial locations are needed to create a spatial grid \n')

  ## calculate sequences for desired stepsize
  # x-axis
  x_range = range(spatlocs$sdimx)
  x_start = x_range[[1]] - minimum_padding
  x_end = x_range[[2]] + minimum_padding
  dimx_steps = ceiling( (x_end-x_start) / sdimx_stepsize)
  dimx_start = mean(c(x_start, x_end))-((dimx_steps/2)*sdimx_stepsize)
  dimx_end = mean(c(x_start, x_end))+((dimx_steps/2)*sdimx_stepsize)
  my_x_seq = seq(from = dimx_start, to = dimx_end, by = sdimx_stepsize)

  # y-axis
  y_range = range(spatlocs$sdimy)
  y_start = y_range[[1]] - minimum_padding
  y_end = y_range[[2]] + minimum_padding
  dimy_steps = ceiling( (y_end-y_start) / sdimy_stepsize)
  dimy_start = mean(c(y_start, y_end))-((dimy_steps/2)*sdimy_stepsize)
  dimy_end = mean(c(y_start, y_end))+((dimy_steps/2)*sdimy_stepsize)
  my_y_seq = seq(from = dimy_start, to = dimy_end, by = sdimy_stepsize)

  # z-axis
  z_range = range(spatlocs$sdimz)
  z_start = z_range[[1]] - minimum_padding
  z_end = z_range[[2]] + minimum_padding
  dimz_steps = ceiling( (z_end-z_start) / sdimz_stepsize)
  dimz_start = mean(c(z_start, z_end))-((dimz_steps/2)*sdimz_stepsize)
  dimz_end = mean(c(z_start, z_end))+((dimz_steps/2)*sdimz_stepsize)
  my_z_seq = seq(from = dimz_start, to = dimz_end, by = sdimz_stepsize)

  ## create grid with starts and ends
  grid_starts = data.table::as.data.table(expand.grid(my_x_seq[-length(my_x_seq)],
                                                      my_y_seq[-length(my_y_seq)],
                                                      my_z_seq[-length(my_z_seq)]))
  colnames(grid_starts) = c('x_start', 'y_start', 'z_start')
  grid_ends = data.table::as.data.table(expand.grid(my_x_seq[-1],
                                                    my_y_seq[-1],
                                                    my_z_seq[-1]))
  colnames(grid_ends) = c('x_end', 'y_end', 'z_end')
  spatgrid = cbind(grid_starts, grid_ends)


  ## first label the grid itself ##
  spatgrid[, gr_name := paste0('gr_', 1:.N)]

  # x-axis
  x_labels = sort(unique(spatgrid$x_start))
  x_gr_names = paste0('gr_x_', 1:length(x_labels))
  names(x_gr_names) = x_labels
  x_gr_names_vector = x_gr_names[as.character(spatgrid$x_start)]
  spatgrid[, gr_x_name := x_gr_names_vector]

  # y-axis
  y_labels = sort(unique(spatgrid$y_start))
  y_gr_names = paste0('gr_y_', 1:length(y_labels))
  names(y_gr_names) = y_labels
  y_gr_names_vector = y_gr_names[as.character(spatgrid$y_start)]
  spatgrid[, gr_y_name := y_gr_names_vector]

  # z-axis
  z_labels = sort(unique(spatgrid$z_start))
  z_gr_names = paste0('gr_z_', 1:length(z_labels))
  names(z_gr_names) = z_labels
  z_gr_names_vector = z_gr_names[as.character(spatgrid$z_start)]
  spatgrid[, gr_z_name := z_gr_names_vector]

  ## for all dimensions ##
  # converter
  gr_dim_names = spatgrid$gr_name
  names(gr_dim_names) = paste0(spatgrid$gr_x_name,'-', spatgrid$gr_y_name, '-', spatgrid$gr_z_name)

  return(spatgrid)

}



#' @title createSpatialDefaultGrid
#' @description Create a spatial grid using the default method
#' @param gobject giotto object
#' @param sdimx_stepsize stepsize along the x-axis
#' @param sdimy_stepsize stepsize along the y-axis
#' @param sdimz_stepsize stepsize along the z-axis
#' @param minimum_padding minimum padding on the edges
#' @param name name for spatial grid (default = 'spatial_grid')
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @return giotto object with updated spatial grid slot
#' @details Creates a spatial grid with defined x, y (and z) dimensions.
#' The dimension units are based on the provided spatial location units.
#' @export
createSpatialDefaultGrid <- function(gobject,
                                     sdimx_stepsize = NULL,
                                     sdimy_stepsize = NULL,
                                     sdimz_stepsize = NULL,
                                     minimum_padding = 1,
                                     name = NULL,
                                     return_gobject = TRUE) {

  # check parameters
  if(is.null(name)) {
    name = 'spatial_grid'
  }

  if(length(c(sdimx_stepsize, sdimy_stepsize, sdimz_stepsize)) == 3) {

    resultgrid = create_spatialGrid_default_3D(gobject = gobject,
                                               sdimx_stepsize = sdimx_stepsize,
                                               sdimy_stepsize = sdimy_stepsize,
                                               sdimz_stepsize = sdimz_stepsize,
                                               minimum_padding = minimum_padding)

  } else if(!is.null(sdimx_stepsize) & !is.null(sdimy_stepsize)) {

    resultgrid = create_spatialGrid_default_2D(gobject = gobject,
                                               sdimx_stepsize = sdimx_stepsize,
                                               sdimy_stepsize = sdimy_stepsize,
                                               minimum_padding = minimum_padding)

  } else {
    cat('\n the stepsize for the x-axis (sdimx) and y-axis (sdimy) is the minimally required \n')
    cat('\n Additionally for a 3D spatial grid the z-axis (sdimz) is also required \n')
  }


  if(return_gobject == TRUE) {

    # 1. check if name has already been used
    spg_names = names(gobject@spatial_grid)

    if(name %in% spg_names) {
      cat('\n ', name, ' has already been used, will be overwritten \n')
    }

    # 2. create spatial grid object
    parameters = list("sdimx_stepsize" = sdimx_stepsize,
                      "sdimy_stepsize" = sdimy_stepsize,
                      "sdimz_stepsize" = sdimz_stepsize,
                      "minimum_padding" = minimum_padding)

    spatgridobj = create_spatialGridObject(name = name,
                                           method = 'default',
                                           parameters = parameters,
                                           gridDT = resultgrid,
                                           outputObj = NULL, # NULL with default
                                           misc = NULL)

    # 3. assign spatial grid object
    gobject@spatial_grid[[name]] <- spatgridobj


    # 4. update log
    ## update parameters used ##
    parameters_list = gobject@parameters
    number_of_rounds = length(parameters_list)
    update_name = paste0(number_of_rounds,'_grid')

    # parameters to include
    parameters_list[[update_name]] = c('name' = name,
                                       'method' = 'default',
                                       'x stepsize' = sdimx_stepsize,
                                       'y stepsize' = sdimy_stepsize,
                                       'z stepsize' = sdimz_stepsize,
                                       'minimum padding' = minimum_padding)

    gobject@parameters = parameters_list

    return(gobject)

  } else {

    return(resultgrid)

  }

}



#' @title select_spatialGrid
#' @description accessor function to select spatial grid
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




#' @title createSpatialGrid
#' @description Create a spatial grid using the default method
#' @param gobject giotto object
#' @param name name for spatial grid
#' @param method method to create a spatial grid
#' @param sdimx_stepsize stepsize along the x-axis
#' @param sdimy_stepsize stepsize along the y-axis
#' @param sdimz_stepsize stepsize along the z-axis
#' @param minimum_padding minimum padding on the edges
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @return giotto object with updated spatial grid slot
#' @details Creates a spatial grid with defined x, y (and z) dimensions.
#' The dimension units are based on the provided spatial location units.
#' \itemize{
#'   \item{default method: }{\code{\link{createSpatialDefaultGrid}}}
#' }
#' @export
createSpatialGrid <- function(gobject,
                              name = NULL,
                              method = c('default'),
                              sdimx_stepsize = NULL,
                              sdimy_stepsize = NULL,
                              sdimz_stepsize = NULL,
                              minimum_padding = 1,
                              return_gobject = TRUE) {

  # get paramters
  method = match.arg(method, c('default'))

  if(method == 'default') {

    out = createSpatialDefaultGrid(gobject = gobject,
                                   sdimx_stepsize = sdimx_stepsize,
                                   sdimy_stepsize = sdimy_stepsize,
                                   sdimz_stepsize = sdimz_stepsize,
                                   minimum_padding = minimum_padding,
                                   name = name,
                                   return_gobject = return_gobject)

  }

  return(out)

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



#' @title annotate_spatlocs_with_spatgrid_2D
#' @description annotate spatial locations with 2D spatial grid information
#' @param spatloc spatial_locs slot from giotto object
#' @param spatgrid selected spatial_grid slot from giotto object
#' @return annotated spatial location data.table
#' @keywords internal
annotate_spatlocs_with_spatgrid_2D = function(spatloc,
                                              spatgrid) {

  ## second label the spatial locations ##
  spatlocs = data.table::copy(spatloc)

  # data.table variables
  gr_x_loc = gr_y_loc = gr_loc = NULL

  x_vector = spatlocs$sdimx
  x_breaks = sort(unique(spatgrid$x_end))
  x_breaks_labels = paste0('gr_x_', 1:length(x_breaks))
  minimum_x = min(spatgrid$x_start)
  my_x_gr = cut(x = x_vector, breaks = c(minimum_x, x_breaks), include.lowest = T, right = T, labels = x_breaks_labels)
  spatlocs[, gr_x_loc := as.character(my_x_gr)]

  y_vector = spatlocs$sdimy
  y_breaks = sort(unique(spatgrid$y_end))
  y_breaks_labels = paste0('gr_y_', 1:length(y_breaks))
  minimum_y = min(spatgrid$y_start)
  my_y_gr = cut(x = y_vector, breaks = c(minimum_y, y_breaks), include.lowest = T, right = T, labels = y_breaks_labels)
  spatlocs[, gr_y_loc := as.character(my_y_gr)]


  ## for all dimensions ##
  # converter
  gr_dim_names = spatgrid$gr_name
  names(gr_dim_names) = paste0(spatgrid$gr_x_name,'-', spatgrid$gr_y_name)

  indiv_dim_names = paste0(spatlocs$gr_x_loc,'-', spatlocs$gr_y_loc)
  my_gr = gr_dim_names[indiv_dim_names]
  spatlocs[, gr_loc := as.character(my_gr)]

  return(spatlocs)

}


#' @title annotate_spatlocs_with_spatgrid_3D
#' @description annotate spatial locations with 3D spatial grid information
#' @param spatloc spatial_locs slot from giotto object
#' @param spatgrid selected spatial_grid slot from giotto object
#' @return annotated spatial location data.table
#' @keywords internal
annotate_spatlocs_with_spatgrid_3D = function(spatloc,
                                              spatgrid) {

  ## second label the spatial locations ##
  spatlocs = data.table::copy(spatloc)

  # data.table variables
  gr_x_loc = gr_y_loc = gr_z_loc = gr_loc = NULL

  x_vector = spatlocs$sdimx
  x_breaks = sort(unique(spatgrid$x_end))
  x_breaks_labels = paste0('gr_x_', 1:length(x_breaks))
  minimum_x = min(spatgrid$x_start)
  my_x_gr = cut(x = x_vector, breaks = c(minimum_x, x_breaks), include.lowest = T, right = T, labels = x_breaks_labels)
  spatlocs[, gr_x_loc := as.character(my_x_gr)]

  y_vector = spatlocs$sdimy
  y_breaks = sort(unique(spatgrid$y_end))
  y_breaks_labels = paste0('gr_y_', 1:length(y_breaks))
  minimum_y = min(spatgrid$y_start)
  my_y_gr = cut(x = y_vector, breaks = c(minimum_y, y_breaks), include.lowest = T, right = T, labels = y_breaks_labels)
  spatlocs[, gr_y_loc := as.character(my_y_gr)]

  z_vector = spatlocs$sdimz
  z_breaks = sort(unique(spatgrid$z_end))
  z_breaks_labels = paste0('gr_z_', 1:length(z_breaks))
  minimum_z = min(spatgrid$z_start)
  my_z_gr = cut(x = z_vector, breaks = c(minimum_z, z_breaks), include.lowest = T, right = T, labels = z_breaks_labels)
  spatlocs[, gr_z_loc := as.character(my_z_gr)]


  ## for all dimensions ##
  # converter
  gr_dim_names = spatgrid$gr_name
  names(gr_dim_names) = paste0(spatgrid$gr_x_name,'-', spatgrid$gr_y_name, '-', spatgrid$gr_z_name)

  indiv_dim_names = paste0(spatlocs$gr_x_loc,'-', spatlocs$gr_y_loc, '-', spatlocs$gr_z_loc)
  my_gr = gr_dim_names[indiv_dim_names]
  spatlocs[, gr_loc := as.character(my_gr)]

  return(spatlocs)

}




#' @title annotateSpatialGrid
#' @description annotate spatial grid with cell ID and cell metadata (optional)
#' @param gobject Giotto object
#' @param spatial_grid_name name of spatial grid, see \code{\link{showGrids}}
#' @param cluster_columns names of cell metadata, see \code{\link{pDataDT}}
#' @return annotated spatial grid data.table
#' @export
annotateSpatialGrid = function(gobject,
                               spatial_grid_name = 'spatial_grid',
                               cluster_columns = NULL) {


  # get grid
  spatial_grid = select_spatialGrid(gobject = gobject,
                                    name = spatial_grid_name)
  spatial_locs = data.table::copy(gobject@spatial_locs)

  # 1. annotate spatial grid with spatial locations
  if(all(c('sdimx', 'sdimy', 'sdimz') %in% colnames(spatial_locs))) {
    annotgrid_locs = annotate_spatlocs_with_spatgrid_3D(spatloc = spatial_locs, spatgrid = spatial_grid)
  } else if(all(c('sdimx', 'sdimy') %in% colnames(spatial_locs))) {
    annotgrid_locs = annotate_spatlocs_with_spatgrid_2D(spatloc = spatial_locs, spatgrid = spatial_grid)
  }

  # 2.select metadata
  cell_metadata = pDataDT(gobject)

  if(!is.null(cluster_columns)) {

    annotation_vector = cluster_columns
    possible_annotations = colnames(cell_metadata)

    missing_annotation = annotation_vector[!annotation_vector %in% possible_annotations]
    if(length(missing_annotation) > 0) {
      cat('These annotations were not found back in the cell metadata (pDataDT): \n',
          missing_annotation, '\n')
    }

    annotation_vector_found = annotation_vector[annotation_vector %in% possible_annotations]
    cell_meta_selected = cell_metadata[, c('cell_ID', annotation_vector_found), with = F]

    annotated_grid = data.table::merge.data.table(x = annotgrid_locs, y = cell_meta_selected, by = 'cell_ID')

    return(annotated_grid)

  } else {

    return(annotgrid_locs)

  }
}







