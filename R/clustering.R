
#' @title doLeidenCluster
#' @name doLeidenCluster
#' @description cluster cells using a NN-network and the Leiden community detection algorithm
#' @param gobject giotto object
#' @param name name for cluster
#' @param nn_network_to_use type of NN network to use (kNN vs sNN)
#' @param network_name name of NN network to use
#' @param python_path specify specific path to python if required
#' @param resolution resolution
#' @param weight_col weight column to use for edges
#' @param partition_type The type of partition to use for optimisation.
#' @param init_membership initial membership of cells for the partition
#' @param n_iterations number of interations to run the Leiden algorithm.
#' If the number of iterations is negative, the Leiden algorithm is run until
#' an iteration in which there was no improvement.
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @param set_seed set seed
#' @param seed_number number for seed
#' @return giotto object with new clusters appended to cell metadata
#' @details
#' This function is a wrapper for the Leiden algorithm implemented in python,
#' which can detect communities in graphs of millions of nodes (cells),
#' as long as they can fit in memory. See the \url{https://github.com/vtraag/leidenalg}{leidenalg}
#' github page or the \url{https://leidenalg.readthedocs.io/en/stable/index.html}{readthedocs}
#' page for more information.
#'
#' Partition types available and information:
#' \itemize{
#'   \item{RBConfigurationVertexPartition: }{Implements Reichardt and Bornholdt’s Potts model
#'    with a configuration null model. This quality function is well-defined only for positive edge weights.
#'     This quality function uses a linear resolution parameter.}
#'   \item{ModularityVertexPartition: }{Implements modularity.
#'    This quality function is well-defined only for positive edge weights. It does \emph{not} use the resolution parameter}
#' }
#'
#' Set \emph{weight_col = NULL} to give equal weight (=1) to each edge.
#'
#' @export
#'
doLeidenCluster = function(gobject,
                           name = 'leiden_clus',
                           nn_network_to_use = 'sNN',
                           network_name = 'sNN.pca',
                           python_path = NULL,
                           resolution = 1,
                           weight_col = 'weight',
                           partition_type = c('RBConfigurationVertexPartition',
                                              'ModularityVertexPartition'),
                           init_membership = NULL,
                           n_iterations = 1000,
                           return_gobject = TRUE,
                           set_seed = T,
                           seed_number = 1234) {

  ## get cell IDs ##
  cell_ID_vec = gobject@cell_ID

  ## select network to use
  igraph_object = select_NearestNetwork(gobject,
                                        nn_network_to_use = nn_network_to_use,
                                        network_name = network_name)

  ## select partition type
  partition_type = match.arg(partition_type,
                             choices = c('RBConfigurationVertexPartition', 'ModularityVertexPartition'))

  ## check or make paths
  # python path
  if(is.null(python_path)) {
    python_path = readGiottoInstructions(gobject, param = "python_path")
  }

  ## prepare python path and louvain script
  reticulate::use_python(required = T, python = python_path)
  python_leiden_function = system.file("python", "python_leiden.py", package = 'Giotto')
  reticulate::source_python(file = python_leiden_function)

  ## start seed
  if(set_seed == TRUE) {
    seed_number = as.integer(seed_number)
  } else {
    seed_number = as.integer(sample(x = 1:10000, size = 1))
  }

  ## extract NN network
  network_edge_dt = data.table::as.data.table(igraph::as_data_frame(x = igraph_object, what = 'edges'))

  # data.table variables
  weight = NULL

  # add weight for edges or set to 1 for all
  if(!is.null(weight_col)) {
    if(!weight_col %in% colnames(network_edge_dt)) {
      stop('\n weight column is not an igraph attribute \n')
    } else {
      # weight is defined by attribute of igraph object
      network_edge_dt = network_edge_dt[,c('from', 'to', weight_col), with = F]
      setnames(network_edge_dt, weight_col, 'weight')
    }
  } else {
    # weight is the same
    network_edge_dt = network_edge_dt[,c('from', 'to'), with = F]
    network_edge_dt[, weight := 1]
  }


  ## do python leiden clustering
  reticulate::py_set_seed(seed = seed_number, disable_hash_randomization = TRUE)
  pyth_leid_result = python_leiden(df = network_edge_dt,
                                   partition_type = partition_type,
                                   initial_membership = init_membership,
                                   weights = 'weight',
                                   n_iterations = n_iterations,
                                   seed = seed_number,
                                   resolution_parameter = resolution)

  ident_clusters_DT = data.table::data.table(cell_ID = pyth_leid_result[[1]], 'name' = pyth_leid_result[[2]])
  data.table::setnames(ident_clusters_DT, 'name', name)


  ## add clusters to metadata ##
  if(return_gobject == TRUE) {

    cluster_names = names(gobject@cell_metadata)
    if(name %in% cluster_names) {
      cat('\n ', name, ' has already been used, will be overwritten \n')
      cell_metadata = gobject@cell_metadata
      cell_metadata[, eval(name) := NULL]
      gobject@cell_metadata = cell_metadata
    }

    gobject = addCellMetadata(gobject = gobject, new_metadata = ident_clusters_DT[, c('cell_ID', name), with = F],
                              by_column = T, column_cell_ID = 'cell_ID')

    ## update parameters used ##
    parameters_list = gobject@parameters
    number_of_rounds = length(parameters_list)
    update_name = paste0(number_of_rounds,'_cluster')

    parameters_list[[update_name]] = c('cluster algorithm' = 'Leiden',
                                       'nn network' = nn_network_to_use,
                                       'network name' = network_name,
                                       'name for clusters' = name,
                                       'pyth leiden resolution' = resolution,
                                       'pyth leiden weight' = weight_col,
                                       'pyth leiden partition' = partition_type,
                                       'pyth leiden iterations' = n_iterations)

    gobject@parameters = parameters_list
    return(gobject)


  } else {

    # else return clustering result
    return(ident_clusters_DT)
  }


}


#' @title doLouvainCluster_community
#' @name doLouvainCluster_community
#' @description cluster cells using a NN-network and the Louvain algorithm from the community module in Python
#' @param gobject giotto object
#' @param name name for cluster
#' @param nn_network_to_use type of NN network to use (kNN vs sNN)
#' @param network_name name of NN network to use
#' @param python_path specify specific path to python if required
#' @param resolution resolution
#' @param weight_col weight column to use for edges
#' @param louv_random Will randomize the node evaluation order and the community evaluation
#' order to get different partitions at each call
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @param set_seed set seed
#' @param seed_number number for seed
#' @return giotto object with new clusters appended to cell metadata
#' @details This function is a wrapper for the Louvain algorithm implemented in Python,
#' which can detect communities in graphs of nodes (cells).
#' See the \url{https://python-louvain.readthedocs.io/en/latest/index.html}{readthedocs}
#' page for more information.
#'
#' Set \emph{weight_col = NULL} to give equal weight (=1) to each edge.
#'
#' @keywords internal
#'
doLouvainCluster_community <- function(gobject,
                                       name = 'louvain_clus',
                                       nn_network_to_use = 'sNN',
                                       network_name = 'sNN.pca',
                                       python_path = NULL,
                                       resolution = 1,
                                       weight_col = NULL,
                                       louv_random = F,
                                       return_gobject = TRUE,
                                       set_seed = F,
                                       seed_number = 1234) {


  ## get cell IDs ##
  cell_ID_vec = gobject@cell_ID

  ## select network to use
  igraph_object = select_NearestNetwork(gobject,
                                        nn_network_to_use = nn_network_to_use,
                                        network_name = network_name)

  ## check or make paths
  # python path
  if(is.null(python_path)) {
    python_path = readGiottoInstructions(gobject, param = "python_path")
  }

  # prepare python path and louvain script
  reticulate::use_python(required = T, python = python_path)
  python_louvain_function = system.file("python", "python_louvain.py", package = 'Giotto')
  reticulate::source_python(file = python_louvain_function)

  # start seed
  if(set_seed == TRUE) {
    seed_number = as.integer(seed_number)
  } else {
    seed_number = as.integer(sample(x = 1:10000, size = 1))
  }

  network_edge_dt = data.table::as.data.table(igraph::as_data_frame(x = igraph_object, what = 'edges'))

  # data.table variables
  weight = NULL

  if(!is.null(weight_col)) {

    if(!weight_col %in% colnames(network_edge_dt)) {
      stop('\n weight column is not an igraph attribute \n')
    } else {
      # weight is defined by attribute of igraph object
      network_edge_dt = network_edge_dt[,c('from', 'to', weight_col), with = F]
      setnames(network_edge_dt, weight_col, 'weight')
    }
  } else {
    # weight is the same
    network_edge_dt = network_edge_dt[,c('from', 'to'), with = F]
    network_edge_dt[, weight := 1]
  }

  # do python louvain clustering
  if(louv_random == FALSE) {
    reticulate::py_set_seed(seed = seed_number, disable_hash_randomization = TRUE)
    pyth_louv_result = python_louvain(df = network_edge_dt, resolution = resolution, randomize = F)
  } else {
    reticulate::py_set_seed(seed = seed_number, disable_hash_randomization = TRUE)
    pyth_louv_result = python_louvain(df = network_edge_dt, resolution = resolution, random_state = seed_number)
  }
  ident_clusters_DT = data.table::data.table(cell_ID = rownames(pyth_louv_result), 'name' = pyth_louv_result[[1]])
  data.table::setnames(ident_clusters_DT, 'name', name)


  ## return
  if(return_gobject == TRUE) {

    cluster_names = names(gobject@cell_metadata)
    if(name %in% cluster_names) {
      cat('\n ', name, ' has already been used, will be overwritten \n')
      cell_metadata = gobject@cell_metadata
      cell_metadata[, eval(name) := NULL]
      gobject@cell_metadata = cell_metadata
    }

    gobject = addCellMetadata(gobject = gobject, new_metadata = ident_clusters_DT[, c('cell_ID', name), with = F],
                              by_column = T, column_cell_ID = 'cell_ID')


    ## update parameters used ##
    parameters_list = gobject@parameters
    number_of_rounds = length(parameters_list)
    update_name = paste0(number_of_rounds,'_cluster')


    parameters_list[[update_name]] = c('cluster algorithm' = 'louvain_community',
                                       'nn network' = nn_network_to_use,
                                       'network name' = network_name,
                                       'name for clusters' = name,
                                       'louvain resolution' = resolution,
                                       'louvain weight' = weight_col)

    gobject@parameters = parameters_list
    return(gobject)


  } else {

    # else return clustering result
    return(ident_clusters_DT)

  }

}


#' @title doLouvainCluster_multinet
#' @name doLouvainCluster_multinet
#' @description cluster cells using a NN-network and the Louvain algorithm from the multinet package in R.
#' @param gobject giotto object
#' @param name name for cluster
#' @param nn_network_to_use type of NN network to use (kNN vs sNN)
#' @param network_name name of NN network to use
#' @param gamma Resolution parameter for modularity in the generalized louvain method.
#' @param omega Inter-layer weight parameter in the generalized louvain method.
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @param set_seed set seed
#' @param seed_number number for seed
#' @return giotto object with new clusters appended to cell metadata
#' @details See \code{\link[multinet]{glouvain_ml}} from the multinet package in R for
#' more information.
#'
#' @keywords internal
doLouvainCluster_multinet <- function(gobject,
                                      name = 'louvain_clus',
                                      nn_network_to_use = 'sNN',
                                      network_name = 'sNN.pca',
                                      gamma = 1,
                                      omega = 1,
                                      return_gobject = TRUE,
                                      set_seed = F,
                                      seed_number = 1234) {


  if("multinet" %in% rownames(installed.packages()) == FALSE) {
    stop("\n package 'multinet' is not yet installed \n",
         "To install: \n",
         "install.packages('multinet')"
    )
  }


  ## get cell IDs ##
  cell_ID_vec = gobject@cell_ID

  ## select network to use
  igraph_object = select_NearestNetwork(gobject,
                                        nn_network_to_use = nn_network_to_use,
                                        network_name = network_name)

  # create mlnetworkobject
  mln_object <- multinet::ml_empty()
  multinet::add_actors_ml(mlnetwork = mln_object, actors = igraph::V(igraph_object))
  multinet::add_igraph_layer_ml(mlnetwork = mln_object, g = igraph_object, name = name)

  # start seed
  if(set_seed == TRUE) {
    set.seed(seed = seed_number)
  }

  # data.table variables
  cell_ID = actor = weight_col = NULL

  louvain_clusters = multinet::glouvain_ml(mlnetwork = mln_object, gamma = gamma, omega = omega)
  ident_clusters_DT = data.table::as.data.table(louvain_clusters)
  ident_clusters_DT[, cell_ID := actor]
  data.table::setnames(ident_clusters_DT, 'cid', name)

  # exit seed
  if(set_seed == TRUE) {
    set.seed(Sys.time())
  }


  ## return
  if(return_gobject == TRUE) {

    cluster_names = names(gobject@cell_metadata)
    if(name %in% cluster_names) {
      cat('\n ', name, ' has already been used, will be overwritten \n')
      cell_metadata = gobject@cell_metadata
      cell_metadata[, eval(name) := NULL]
      gobject@cell_metadata = cell_metadata
    }

    gobject = addCellMetadata(gobject = gobject, new_metadata = ident_clusters_DT[, c('cell_ID', name), with = F],
                              by_column = T, column_cell_ID = 'cell_ID')


    ## update parameters used ##
    parameters_list = gobject@parameters
    number_of_rounds = length(parameters_list)
    update_name = paste0(number_of_rounds,'_cluster')


    parameters_list[[update_name]] = c('cluster algorithm' = 'louvain_multinet',
                                       'nn network' = nn_network_to_use,
                                       'network name' = network_name,
                                       'name for clusters' = name,
                                       'louvain gamma' = gamma,
                                       'louvain omega' = omega,
                                       'louvain weight' = weight_col)

    gobject@parameters = parameters_list
    return(gobject)


  } else {

    # else return clustering result
    return(ident_clusters_DT)

  }

}


#' @title doLouvainCluster
#' @name doLouvainCluster
#' @description cluster cells using a NN-network and the Louvain algorithm.
#' @param gobject giotto object
#' @param version implemented version of Louvain clustering to use
#' @param name name for cluster
#' @param nn_network_to_use type of NN network to use (kNN vs sNN)
#' @param network_name name of NN network to use
#' @param python_path [community] specify specific path to python if required
#' @param resolution [community] resolution
#' @param louv_random [community] Will randomize the node evaluation order and the community evaluation
#' order to get different partitions at each call
#' @param weight_col weight column name
#' @param gamma [multinet] Resolution parameter for modularity in the generalized louvain method.
#' @param omega [multinet] Inter-layer weight parameter in the generalized louvain method
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @param set_seed set seed
#' @param seed_number number for seed
#' @param \dots additional parameters
#' @return giotto object with new clusters appended to cell metadata
#' @details Louvain clustering using the community or multinet implementation of the louvain clustering algorithm.
#' @seealso \code{\link{doLouvainCluster_community}} and \code{\link{doLouvainCluster_multinet}}
#' @export
doLouvainCluster = function(gobject,
                            version = c('community', 'multinet'),
                            name = 'louvain_clus',
                            nn_network_to_use = 'sNN',
                            network_name = 'sNN.pca',
                            python_path = NULL,
                            resolution = 1,
                            weight_col = NULL,
                            gamma = 1,
                            omega = 1,
                            louv_random = F,
                            return_gobject = TRUE,
                            set_seed = F,
                            seed_number = 1234,
                            ...) {

  ## louvain clustering version to use
  version = match.arg(version, c('community', 'multinet'))


  # python community implementation
  if(version == 'community') {

    result = doLouvainCluster_community(gobject = gobject,
                                        name = name,
                                        nn_network_to_use = nn_network_to_use,
                                        network_name = network_name,
                                        python_path = python_path,
                                        resolution = resolution,
                                        weight_col = weight_col,
                                        louv_random = louv_random,
                                        return_gobject = return_gobject,
                                        set_seed = set_seed,
                                        seed_number = seed_number,
                                        ...)

    return(result)

    ## r multinet implementation
  } else if(version == 'multinet') {

    result = doLouvainCluster_multinet(gobject = gobject,
                                       name = name,
                                       nn_network_to_use = nn_network_to_use,
                                       network_name = network_name,
                                       gamma = gamma,
                                       omega = omega,
                                       return_gobject = return_gobject,
                                       set_seed = set_seed,
                                       seed_number = seed_number,
                                       ...)

    return(result)
  }


}



#' @title doRandomWalkCluster
#' @name doRandomWalkCluster
#' @description Cluster cells using a random walk approach.
#' @param gobject giotto object
#' @param name name for cluster
#' @param nn_network_to_use type of NN network to use (kNN vs sNN)
#' @param network_name name of NN network to use
#' @param walk_steps number of walking steps
#' @param walk_clusters number of final clusters
#' @param walk_weights cluster column defining the walk weights
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @param set_seed set seed
#' @param seed_number number for seed
#' @return giotto object with new clusters appended to cell metadata
#' @details See \code{\link[igraph]{cluster_walktrap}} function from the igraph
#' package in R for more information.
#' @export
doRandomWalkCluster <- function(gobject,
                                name = 'random_walk_clus',
                                nn_network_to_use = 'sNN',
                                network_name = 'sNN.pca',
                                walk_steps = 4,
                                walk_clusters = 10,
                                walk_weights = NA,
                                return_gobject = TRUE,
                                set_seed = F,
                                seed_number = 1234) {

  ## get cell IDs ##
  cell_ID_vec = gobject@cell_ID

  ## select network to use
  igraph_object = select_NearestNetwork(gobject,
                                        nn_network_to_use = nn_network_to_use,
                                        network_name = network_name)


  # start seed
  if(set_seed == TRUE) {
    set.seed(seed = seed_number)
  }

  randomwalk_clusters <- igraph::cluster_walktrap(graph = igraph_object, steps = walk_steps, weights = walk_weights)
  randomwalk_clusters <- as.factor(igraph::cut_at(communities = randomwalk_clusters, no =  walk_clusters))

  ident_clusters_DT <- data.table::data.table('cell_ID' = igraph::V(igraph_object)$name, 'name' = randomwalk_clusters)
  data.table::setnames(ident_clusters_DT, 'name', name)

  # exit seed
  if(set_seed == TRUE) {
    set.seed(Sys.time())
  }


  ## return
  if(return_gobject == TRUE) {

    cluster_names = names(gobject@cell_metadata)
    if(name %in% cluster_names) {
      cat('\n ', name, ' has already been used, will be overwritten \n')
      cell_metadata = gobject@cell_metadata
      cell_metadata[, eval(name) := NULL]
      gobject@cell_metadata = cell_metadata
    }

    gobject = addCellMetadata(gobject = gobject, new_metadata = ident_clusters_DT[, c('cell_ID', name), with = F],
                              by_column = T, column_cell_ID = 'cell_ID')


    ## update parameters used ##
    parameters_list = gobject@parameters
    number_of_rounds = length(parameters_list)
    update_name = paste0(number_of_rounds,'_cluster')


    parameters_list[[update_name]] = c('cluster algorithm' = 'random walk',
                                       'nn network' = nn_network_to_use,
                                       'network name' = network_name,
                                       'name for clusters' = name,
                                       'random walk steps' = walk_steps,
                                       'random walk clusters' = walk_clusters,
                                       'random walk weights' = walk_weights)

    gobject@parameters = parameters_list
    return(gobject)


  } else {

    # else return clustering result
    return(ident_clusters_DT)

  }
}


#' @title doSNNCluster
#' @name doSNNCluster
#' @description Cluster cells using a SNN cluster approach.
#' @param gobject giotto object
#' @param name name for cluster
#' @param nn_network_to_use type of NN network to use (only works on kNN)
#' @param network_name name of kNN network to use
#' @param k Neighborhood size for nearest neighbor sparsification to create the shared NN graph.
#' @param eps Two objects are only reachable from each other if they share at least eps nearest neighbors.
#' @param minPts minimum number of points that share at least eps nearest neighbors for a point to be considered a core points.
#' @param borderPoints should borderPoints be assigned to clusters like in DBSCAN?
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @param set_seed set seed
#' @param seed_number number for seed
#' @return giotto object with new clusters appended to cell metadata
#' @details See \code{\link[dbscan]{sNNclust}} from dbscan package
#' @export
doSNNCluster <- function(gobject,
                         name = 'sNN_clus',
                         nn_network_to_use = 'kNN',
                         network_name = 'kNN.pca',
                         k = 20,
                         eps = 4,
                         minPts = 16,
                         borderPoints = TRUE,
                         return_gobject = TRUE,
                         set_seed = F,
                         seed_number = 1234) {


  ## get cell IDs ##
  cell_ID_vec = gobject@cell_ID

  ## select network to use
  igraph_object = select_NearestNetwork(gobject,
                                        nn_network_to_use = nn_network_to_use,
                                        network_name = network_name)


  if(nn_network_to_use == 'sNN') {
    stop('\n sNNclust can only be used with kNN-network \n')
  }


  # start seed
  if(set_seed == TRUE) {
    set.seed(seed = seed_number)
  }

  # data.table variables
  from = from_T = to_T = to = weight = distance = NULL

  ## SNN clust
  igraph_DT = data.table::as.data.table(igraph::as_data_frame(igraph_object, what = 'edges'))
  igraph_DT = igraph_DT[order(from)]

  cell_id_numeric = unique(x = c(igraph_DT$from, igraph_DT$to))
  names(cell_id_numeric) <- 1:length(cell_id_numeric)
  igraph_DT[, from_T := as.numeric(names(cell_id_numeric[cell_id_numeric == from])), by = 1:nrow(igraph_DT)]
  igraph_DT[, to_T := as.numeric(names(cell_id_numeric[cell_id_numeric == to])), by = 1:nrow(igraph_DT)]
  temp_igraph_DT = igraph_DT[,.(from_T, to_T, weight, distance)]
  data.table::setnames(temp_igraph_DT, old = c('from_T', 'to_T'), new = c('from', 'to'))

  kNN_object = nnDT_to_kNN(nnDT = temp_igraph_DT)
  sNN_clusters = dbscan::sNNclust(x = kNN_object, k = k, eps = eps,
                                  minPts = minPts, borderPoints = borderPoints)

  ident_clusters_DT <- data.table::data.table('cell_ID' = cell_id_numeric[1:nrow(kNN_object$dist)], 'name' = sNN_clusters$cluster)
  data.table::setnames(ident_clusters_DT, 'name', name)

  # exit seed
  if(set_seed == TRUE) {
    set.seed(Sys.time())
  }

  ## add clusters to metadata ##
  if(return_gobject == TRUE) {

    cluster_names = names(gobject@cell_metadata)
    if(name %in% cluster_names) {
      cat('\n ', name, ' has already been used, will be overwritten \n')
      cell_metadata = gobject@cell_metadata
      cell_metadata[, eval(name) := NULL]
      gobject@cell_metadata = cell_metadata
    }

    gobject = addCellMetadata(gobject = gobject, new_metadata = ident_clusters_DT[, c('cell_ID', name), with = F],
                              by_column = T, column_cell_ID = 'cell_ID')


    ## update parameters used ##
    parameters_list = gobject@parameters
    number_of_rounds = length(parameters_list)
    update_name = paste0(number_of_rounds,'_cluster')

    parameters_list[[update_name]] = c('cluster algorithm' = 'SNN cluster',
                                       'nn network' = nn_network_to_use,
                                       'network name' = network_name,
                                       'name for clusters' = name,
                                       'k for sNNclust' = k,
                                       'eps for sNNclust' = eps,
                                       'minPts for sNNcluster' = minPts,
                                       'assign borderPoints' = borderPoints)

    gobject@parameters = parameters_list
    return(gobject)


  } else {

    # else return clustering result
    return(ident_clusters_DT)

  }

}





#' @title doKmeans
#' @name doKmeans
#' @description cluster cells using kmeans algorithm
#' @param gobject giotto object
#' @param expression_values expression values to use
#' @param genes_to_use subset of genes to use
#' @param dim_reduction_to_use dimension reduction to use
#' @param dim_reduction_name dimensions reduction name
#' @param dimensions_to_use dimensions to use
#' @param distance_method distance method
#' @param centers number of final clusters
#' @param iter_max kmeans maximum iterations
#' @param nstart kmeans nstart
#' @param algorithm kmeans algorithm
#' @param name name for kmeans clustering
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @param set_seed set seed
#' @param seed_number number for seed
#' @return giotto object with new clusters appended to cell metadata
#' @details Description on how to use Kmeans clustering method.
#' @seealso  \code{\link[stats]{kmeans}}
#' @export
#' @examples
#'
#' data(mini_giotto_single_cell)
#'
#' mini_giotto_single_cell = doKmeans(mini_giotto_single_cell, centers = 4, name = 'kmeans_clus')
#' plotUMAP_2D(mini_giotto_single_cell, cell_color = 'kmeans_clus', point_size = 3)
#'
doKmeans <- function(gobject,
                     expression_values = c('normalized', 'scaled', 'custom'),
                     genes_to_use = NULL,
                     dim_reduction_to_use = c('cells', 'pca', 'umap', 'tsne'),
                     dim_reduction_name = 'pca',
                     dimensions_to_use = 1:10,
                     distance_method = c("original", "pearson", "spearman",
                                         "euclidean", "maximum", "manhattan",
                                         "canberra", "binary", "minkowski"),
                     centers = 10,
                     iter_max = 100,
                     nstart = 1000,
                     algorithm = "Hartigan-Wong",
                     name = 'kmeans',
                     return_gobject = TRUE,
                     set_seed = T,
                     seed_number = 1234) {



  dim_reduction_to_use = match.arg(dim_reduction_to_use, choices = c('cells', 'pca', 'umap', 'tsne'))
  distance_method = match.arg(distance_method, choices = c("original", "pearson", "spearman",
                                                           "euclidean", "maximum", "manhattan",
                                                           "canberra", "binary", "minkowski"))
  values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))

  ## using dimension reduction ##
  if(dim_reduction_to_use != 'cells' & !is.null(dim_reduction_to_use)) {

    ## TODO: check if reduction exists

    # use only available dimensions if dimensions < dimensions_to_use
    dim_coord = gobject@dimension_reduction[['cells']][[dim_reduction_to_use]][[dim_reduction_name]][['coordinates']]
    dimensions_to_use = dimensions_to_use[dimensions_to_use %in% 1:ncol(dim_coord)]
    matrix_to_use = dim_coord[, dimensions_to_use]


  } else {
    ## using original matrix ##
    expr_values = select_expression_values(gobject = gobject, values = values)

    # subset expression matrix
    if(!is.null(genes_to_use)) {
      expr_values = expr_values[rownames(expr_values) %in% genes_to_use, ]
    }

    # features as columns
    # cells as rows
    matrix_to_use = t_giotto(expr_values)

  }

  ## distance
  if(distance_method == 'original') {
    celldist = matrix_to_use
  } else if(distance_method %in% c('spearman', 'pearson')) {
    celldist = stats::as.dist(1-cor_giotto(x = t_giotto(matrix_to_use), method = distance_method))
  } else if(distance_method %in% c("euclidean", "maximum", "manhattan",
                                   "canberra", "binary", "minkowski")) {
    celldist = stats::dist(x = matrix_to_use, method = distance_method)
  }

  ## kmeans clustering
  # set seed
  if(set_seed == TRUE) {
    seed_number = as.integer(seed_number)
  } else {
    seed_number = as.integer(sample(x = 1:10000, size = 1))
  }
  set.seed(seed = seed_number)

  # start clustering
  kclusters = stats::kmeans(x = celldist, centers = centers,
                            iter.max = iter_max, nstart = nstart,
                            algorithm =  algorithm)


  ident_clusters_DT = data.table::data.table(cell_ID = names(kclusters[['cluster']]),
                                             'name' = kclusters[['cluster']])
  data.table::setnames(ident_clusters_DT, 'name', name)


  ## add clusters to metadata ##
  if(return_gobject == TRUE) {

    cluster_names = names(gobject@cell_metadata)
    if(name %in% cluster_names) {
      cat('\n ', name, ' has already been used, will be overwritten \n')
      cell_metadata = gobject@cell_metadata
      cell_metadata[, eval(name) := NULL]
      gobject@cell_metadata = cell_metadata
    }

    gobject = addCellMetadata(gobject = gobject, new_metadata = ident_clusters_DT[, c('cell_ID', name), with = F],
                              by_column = T, column_cell_ID = 'cell_ID')


    ## update parameters used ##
    parameters_list = gobject@parameters
    number_of_rounds = length(parameters_list)
    update_name = paste0(number_of_rounds,'_kmeans_cluster')
    # parameters to include


    parameters_list[[update_name]] = c('expression values' = values,
                                       'dim reduction used' = dim_reduction_to_use,
                                       'dim reduction name' = dim_reduction_name,
                                       'name for clusters' = name,
                                       'distance method' = distance_method,
                                       'centers' = centers,
                                       'iter_max' = iter_max,
                                       'nstart' =  nstart
    )

    gobject@parameters = parameters_list

    return(gobject)

  } else {

    return(ident_clusters_DT)

  }

}





#' @title doHclust
#' @name doHclust
#' @description cluster cells using hierarchical clustering algorithm
#' @param gobject giotto object
#' @param expression_values expression values to use
#' @param genes_to_use subset of genes to use
#' @param dim_reduction_to_use dimension reduction to use
#' @param dim_reduction_name dimensions reduction name
#' @param dimensions_to_use dimensions to use
#' @param distance_method distance method
#' @param agglomeration_method agglomeration method for hclust
#' @param k number of final clusters
#' @param h cut hierarchical tree at height = h
#' @param name name for hierarchical clustering
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @param set_seed set seed
#' @param seed_number number for seed
#' @return giotto object with new clusters appended to cell metadata
#' @details Description on how to use Kmeans clustering method.
#' @seealso  \code{\link[stats]{hclust}}
#' @export
#' @examples
#'
#' data(mini_giotto_single_cell)
#'
#' mini_giotto_single_cell = doHclust(mini_giotto_single_cell, k = 4, name = 'hier_clus')
#' plotUMAP_2D(mini_giotto_single_cell, cell_color = 'hier_clus', point_size = 3)
#'
doHclust <- function(gobject,
                     expression_values = c('normalized', 'scaled', 'custom'),
                     genes_to_use = NULL,
                     dim_reduction_to_use = c('cells', 'pca', 'umap', 'tsne'),
                     dim_reduction_name = 'pca',
                     dimensions_to_use = 1:10,
                     distance_method = c("pearson", "spearman", "original",
                                         "euclidean", "maximum", "manhattan",
                                         "canberra", "binary", "minkowski"),
                     agglomeration_method = c("ward.D2","ward.D", "single",
                                              "complete", "average", "mcquitty",
                                              "median", "centroid" ),
                     k = 10,
                     h = NULL,
                     name = 'hclust',
                     return_gobject = TRUE,
                     set_seed = T,
                     seed_number = 1234) {



  dim_reduction_to_use = match.arg(dim_reduction_to_use, choices = c('cells', 'pca', 'umap', 'tsne'))
  distance_method = match.arg(distance_method, choices = c("pearson", "spearman",  "original",
                                                           "euclidean", "maximum", "manhattan",
                                                           "canberra", "binary", "minkowski"))
  agglomeration_method = match.arg(agglomeration_method, choices = c("ward.D2","ward.D", "single",
                                                                     "complete", "average", "mcquitty",
                                                                     "median", "centroid" ))
  values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))


  ## using dimension reduction ##
  if(dim_reduction_to_use != 'cells' & !is.null(dim_reduction_to_use)) {

    ## TODO: check if reduction exists

    # use only available dimensions if dimensions < dimensions_to_use
    dim_coord = gobject@dimension_reduction[['cells']][[dim_reduction_to_use]][[dim_reduction_name]][['coordinates']]
    dimensions_to_use = dimensions_to_use[dimensions_to_use %in% 1:ncol(dim_coord)]
    matrix_to_use = dim_coord[, dimensions_to_use]


  } else {
    ## using original matrix ##
    expr_values = select_expression_values(gobject = gobject, values = values)

    # subset expression matrix
    if(!is.null(genes_to_use)) {
      expr_values = expr_values[rownames(expr_values) %in% genes_to_use, ]
    }

    # features as columns
    # cells as rows
    matrix_to_use = t_giotto(expr_values)

  }

  ## distance
  if(distance_method == 'original') {
    celldist = matrix_to_use
  } else if(distance_method %in% c('spearman', 'pearson')) {
    celldist = stats::as.dist(1-cor_giotto(x = t_giotto(matrix_to_use), method = distance_method))
  } else if(distance_method %in% c("euclidean", "maximum", "manhattan",
                                   "canberra", "binary", "minkowski")) {
    celldist = stats::dist(x = matrix_to_use, method = distance_method)
  }

  ## hierarchical clustering
  # set seed
  if(set_seed == TRUE) {
    seed_number = as.integer(seed_number)
  } else {
    seed_number = as.integer(sample(x = 1:10000, size = 1))
  }
  set.seed(seed = seed_number)

  # start clustering
  hclusters = stats::hclust(d = celldist, method = agglomeration_method)
  hclusters_cut = stats::cutree(tree = hclusters, k = k, h = h)

  ident_clusters_DT = data.table::data.table(cell_ID = names(hclusters_cut),
                                             'name' = hclusters_cut)
  data.table::setnames(ident_clusters_DT, 'name', name)


  ## add clusters to metadata ##
  if(return_gobject == TRUE) {

    cluster_names = names(gobject@cell_metadata)
    if(name %in% cluster_names) {
      cat('\n ', name, ' has already been used, will be overwritten \n')
      cell_metadata = gobject@cell_metadata
      cell_metadata[, eval(name) := NULL]
      gobject@cell_metadata = cell_metadata
    }

    gobject = addCellMetadata(gobject = gobject, new_metadata = ident_clusters_DT[, c('cell_ID', name), with = F],
                              by_column = T, column_cell_ID = 'cell_ID')


    ## update parameters used ##
    parameters_list = gobject@parameters
    number_of_rounds = length(parameters_list)
    update_name = paste0(number_of_rounds,'_hierarchical_cluster')
    # parameters to include


    parameters_list[[update_name]] = c('expression values' = values,
                                       'dim reduction used' = dim_reduction_to_use,
                                       'dim reduction name' = dim_reduction_name,
                                       'name for clusters' = name,
                                       'distance method' = distance_method,
                                       'k' = k,
                                       'h' = h
    )

    gobject@parameters = parameters_list

    return(gobject)

  } else {

    return(list('hclust' = hclusters, 'DT' = ident_clusters_DT))

  }

}





#' @title clusterCells
#' @name clusterCells
#' @description cluster cells using a variety of different methods
#' @param gobject giotto object
#' @param name name for new clustering result
#' @param nn_network_to_use type of NN network to use (kNN vs sNN)
#' @param network_name name of NN network to use
#' @param cluster_method community cluster method to use
#' @param pyth_leid_resolution resolution for leiden
#' @param pyth_leid_weight_col column to use for weights
#' @param pyth_leid_part_type partition type to use
#' @param pyth_leid_init_memb initial membership
#' @param pyth_leid_iterations number of iterations
#' @param pyth_louv_resolution resolution for louvain
#' @param pyth_louv_weight_col python louvain param: weight column
#' @param python_louv_random python louvain param: random
#' @param python_path specify specific path to python if required
#' @param louvain_gamma louvain param: gamma or resolution
#' @param louvain_omega louvain param: omega
#' @param walk_steps randomwalk: number of steps
#' @param walk_clusters randomwalk: number of clusters
#' @param walk_weights randomwalk: weight column
#' @param sNNclust_k SNNclust: k neighbors to use
#' @param sNNclust_eps SNNclust: epsilon
#' @param sNNclust_minPts SNNclust: min points
#' @param borderPoints SNNclust: border points
#' @param expression_values expression values to use
#' @param genes_to_use = NULL,
#' @param dim_reduction_to_use dimension reduction to use
#' @param dim_reduction_name name of reduction 'pca',
#' @param dimensions_to_use dimensions to use
#' @param distance_method distance method
#' @param km_centers kmeans centers
#' @param km_iter_max kmeans iterations
#' @param km_nstart kmeans random starting points
#' @param km_algorithm kmeans algorithm
#' @param hc_agglomeration_method hierarchical clustering method
#' @param hc_k hierachical number of clusters
#' @param hc_h hierarchical tree cutoff
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @param set_seed set seed
#' @param seed_number number for seed
#' @return giotto object with new clusters appended to cell metadata
#' @details Wrapper for the different clustering methods.
#' @seealso \code{\link{doLeidenCluster}}, \code{\link{doLouvainCluster_community}}, \code{\link{doLouvainCluster_multinet}},
#' \code{\link{doLouvainCluster}}, \code{\link{doRandomWalkCluster}}, \code{\link{doSNNCluster}},
#' \code{\link{doKmeans}}, \code{\link{doHclust}}
#' @export
clusterCells <- function(gobject,
                         cluster_method = c('leiden',
                                            'louvain_community', 'louvain_multinet',
                                            'randomwalk', 'sNNclust',
                                            'kmeans', 'hierarchical'),
                         name = 'cluster_name',

                         nn_network_to_use = 'sNN',
                         network_name = 'sNN.pca',

                         pyth_leid_resolution = 1,
                         pyth_leid_weight_col = 'weight',
                         pyth_leid_part_type = c('RBConfigurationVertexPartition', 'ModularityVertexPartition'),
                         pyth_leid_init_memb = NULL,
                         pyth_leid_iterations = 1000,

                         pyth_louv_resolution = 1,
                         pyth_louv_weight_col = NULL,
                         python_louv_random = F,

                         python_path = NULL,

                         louvain_gamma = 1,
                         louvain_omega = 1,

                         walk_steps = 4,
                         walk_clusters = 10,
                         walk_weights = NA,

                         sNNclust_k = 20,
                         sNNclust_eps = 4,
                         sNNclust_minPts = 16,
                         borderPoints = TRUE,

                         expression_values = c('normalized', 'scaled', 'custom'),
                         genes_to_use = NULL,
                         dim_reduction_to_use = c('cells', 'pca', 'umap', 'tsne'),
                         dim_reduction_name = 'pca',
                         dimensions_to_use = 1:10,
                         distance_method = c("original", "pearson", "spearman",
                                             "euclidean", "maximum", "manhattan",
                                             "canberra", "binary", "minkowski"),
                         km_centers = 10,
                         km_iter_max = 100,
                         km_nstart = 1000,
                         km_algorithm = "Hartigan-Wong",

                         hc_agglomeration_method = c("ward.D2","ward.D", "single",
                                                     "complete", "average", "mcquitty",
                                                     "median", "centroid" ),
                         hc_k = 10,
                         hc_h = NULL,


                         return_gobject = TRUE,
                         set_seed = T,
                         seed_number = 1234) {

  ## select cluster method
  cluster_method = match.arg(arg = cluster_method, choices = c('leiden',
                                                               'louvain_community', 'louvain_multinet',
                                                               'randomwalk', 'sNNclust',
                                                               'kmeans', 'hierarchical'))


  if(cluster_method == 'leiden') {

    result = doLeidenCluster(gobject = gobject,
                             name = name,
                             nn_network_to_use = nn_network_to_use,
                             network_name = network_name,
                             python_path = python_path,
                             resolution = pyth_leid_resolution,
                             weight_col = pyth_leid_weight_col,
                             partition_type = pyth_leid_part_type,
                             init_membership = pyth_leid_init_memb,
                             n_iterations = pyth_leid_iterations,
                             return_gobject = return_gobject,
                             set_seed = set_seed,
                             seed_number = seed_number)

  } else if(cluster_method == 'louvain_community') {

    result = doLouvainCluster_community(gobject = gobject,
                                        name = name,
                                        nn_network_to_use = nn_network_to_use,
                                        network_name = network_name,
                                        python_path = python_path,
                                        resolution = pyth_louv_resolution,
                                        weight_col = pyth_louv_weight_col,
                                        louv_random = python_louv_random,
                                        return_gobject = return_gobject,
                                        set_seed = set_seed,
                                        seed_number = seed_number)

  } else if(cluster_method == 'louvain_multinet') {

    result = doLouvainCluster_multinet(gobject = gobject,
                                       name = name,
                                       nn_network_to_use = nn_network_to_use,
                                       network_name = network_name,
                                       gamma = louvain_gamma,
                                       omega = louvain_omega,
                                       return_gobject = return_gobject,
                                       set_seed = set_seed,
                                       seed_number = seed_number)

  } else if(cluster_method == 'randomwalk') {

    result = doRandomWalkCluster(gobject = gobject,
                                 name = name,
                                 nn_network_to_use = nn_network_to_use,
                                 network_name = network_name,
                                 walk_steps = walk_steps,
                                 walk_clusters = walk_clusters,
                                 walk_weights = walk_weights,
                                 return_gobject = return_gobject,
                                 set_seed = set_seed,
                                 seed_number = seed_number)

  } else if(cluster_method == 'sNNclust') {

    result = doSNNCluster(gobject = gobject,
                          name = name,
                          nn_network_to_use = nn_network_to_use,
                          network_name = network_name,
                          k = sNNclust_k,
                          eps = sNNclust_eps,
                          minPts = sNNclust_minPts,
                          borderPoints = borderPoints,
                          return_gobject = return_gobject,
                          set_seed = set_seed,
                          seed_number = seed_number)

  } else if(cluster_method == 'kmeans') {

    result = doKmeans(gobject = gobject,
                      name = name,
                      expression_values = expression_values,
                      genes_to_use = genes_to_use,
                      dim_reduction_to_use = dim_reduction_to_use,
                      dim_reduction_name = dim_reduction_name,
                      dimensions_to_use = dimensions_to_use,
                      distance_method = distance_method,
                      centers = km_centers,
                      iter_max = km_iter_max,
                      nstart = km_nstart,
                      algorithm = km_algorithm,
                      return_gobject = return_gobject,
                      set_seed = set_seed,
                      seed_number = seed_number)

  } else if(cluster_method == 'hierarchical') {

    result = doHclust(gobject = gobject,
                      name = name,
                      expression_values = expression_values,
                      genes_to_use = genes_to_use,
                      dim_reduction_to_use = dim_reduction_to_use,
                      dim_reduction_name = dim_reduction_name,
                      dimensions_to_use = dimensions_to_use,
                      distance_method = distance_method,
                      agglomeration_method = hc_agglomeration_method,
                      k = hc_k,
                      h = hc_h,
                      return_gobject = return_gobject,
                      set_seed = set_seed,
                      seed_number = seed_number)

  }


  return(result)

}






#' @title doLeidenSubCluster
#' @name doLeidenSubCluster
#' @description Further subcluster cells using a NN-network and the Leiden algorithm
#' @param gobject giotto object
#' @param name name for new clustering result
#' @param cluster_column cluster column to subcluster
#' @param selected_clusters only do subclustering on these clusters
#' @param hvg_param parameters for calculateHVG
#' @param hvg_min_perc_cells threshold for detection in min percentage of cells
#' @param hvg_mean_expr_det threshold for mean expression level in cells with detection
#' @param use_all_genes_as_hvg forces all genes to be HVG and to be used as input for PCA
#' @param min_nr_of_hvg minimum number of HVG, or all genes will be used as input for PCA
#' @param pca_param parameters for runPCA
#' @param nn_param parameters for parameters for createNearestNetwork
#' @param k_neighbors number of k for createNearestNetwork
#' @param resolution resolution of Leiden clustering
#' @param n_iterations number of interations to run the Leiden algorithm.
#' @param python_path specify specific path to python if required
#' @param nn_network_to_use type of NN network to use (kNN vs sNN)
#' @param network_name name of NN network to use
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @param verbose verbose
#' @return giotto object with new subclusters appended to cell metadata
#' @details This function performs subclustering using the Leiden algorithm on selected clusters.
#' The systematic steps are:
#' \itemize{
#'   \item{1. subset Giotto object}
#'   \item{2. identify highly variable genes}
#'   \item{3. run PCA}
#'   \item{4. create nearest neighbouring network}
#'   \item{5. do Leiden clustering}
#' }
#' @seealso \code{\link{doLeidenCluster}}
#' @export
doLeidenSubCluster = function(gobject,
                              name = 'sub_pleiden_clus',
                              cluster_column = NULL,
                              selected_clusters = NULL,
                              hvg_param = list(reverse_log_scale = T, difference_in_cov = 1, expression_values = 'normalized'),
                              hvg_min_perc_cells = 5,
                              hvg_mean_expr_det = 1,
                              use_all_genes_as_hvg = FALSE,
                              min_nr_of_hvg = 5,
                              pca_param = list(expression_values = 'normalized', scale_unit = T),
                              nn_param = list(dimensions_to_use = 1:20),
                              k_neighbors = 10,
                              resolution = 0.5,
                              n_iterations = 500,
                              python_path = NULL,
                              nn_network_to_use = 'sNN',
                              network_name = 'sNN.pca',
                              return_gobject = TRUE,
                              verbose = T) {


  iter_list = list()

  cell_metadata = pDataDT(gobject)

  if(is.null(cluster_column)) {
    stop('\n You need to provide a cluster column to subcluster on \n')
  }
  unique_clusters = sort(unique(cell_metadata[[cluster_column]]))


  # data.table variables
  hvg = perc_cells = mean_expr_det = parent_cluster = comb = tempclus = NULL


  for(cluster in unique_clusters) {

    if(verbose == TRUE) cat('\n start with cluster: ', cluster, '\n')

    ## get subset
    subset_cell_IDs = cell_metadata[get(cluster_column) == cluster][['cell_ID']]
    temp_giotto = subsetGiotto(gobject = gobject, cell_ids = subset_cell_IDs)

    ## if cluster is not selected
    if(!is.null(selected_clusters) & !cluster %in% selected_clusters) {

      temp_cluster = data.table('cell_ID' = subset_cell_IDs, 'tempclus' = 1, 'parent_cluster' = cluster)
      iter_list[[cluster]] = temp_cluster

    } else {
      # continue for selected clusters or all clusters if there is no selection

      ## calculate stats
      temp_giotto <- addStatistics(gobject = temp_giotto)

      ## calculate variable genes
      temp_giotto = do.call('calculateHVG', c(gobject = temp_giotto, hvg_param))

      ## get hvg
      gene_metadata = fDataDT(temp_giotto)
      featgenes     = gene_metadata[hvg == 'yes' & perc_cells >= hvg_min_perc_cells & mean_expr_det >= hvg_mean_expr_det]$gene_ID

      ## catch too low number of hvg
      if(use_all_genes_as_hvg == TRUE) {
        featgenes == gene_metadata$gene_ID
      } else {
        if(verbose == TRUE) cat('\n', length(featgenes), 'highly variable genes have been selected \n')
        if(length(featgenes) <= min_nr_of_hvg) {
          cat('\n too few genes, will continue with all genes instead \n')
          featgenes = gene_metadata$gene_ID
        }
      }

      ## run PCA
      temp_giotto = do.call('runPCA', c(gobject =  temp_giotto, genes_to_use = list(featgenes), pca_param))

      ## nearest neighbor and clustering
      temp_giotto = do.call('createNearestNetwork', c(gobject = temp_giotto, k = k_neighbors, nn_param))

      ## Leiden Cluster
      ## TO DO: expand to all clustering options
      temp_cluster = doLeidenCluster(gobject = temp_giotto,
                                     resolution = resolution,
                                     n_iterations = n_iterations,
                                     python_path = python_path,
                                     name = 'tempclus',
                                     return_gobject = F)

      temp_cluster[, parent_cluster := cluster]

      iter_list[[cluster]] = temp_cluster



    }

  }

  together = do.call('rbind', iter_list)
  together[, comb := paste0(parent_cluster,'.',tempclus)]

  # rename with subcluster of original name
  #new_cluster_column = paste0(cluster_column,'_sub')
  setnames(together, 'comb', name)


  if(return_gobject == TRUE) {

    cluster_names = names(gobject@cell_metadata)
    if(name %in% cluster_names) {
      cat('\n ', name, ' has already been used, will be overwritten \n')
      cell_metadata = gobject@cell_metadata
      cell_metadata[, eval(name) := NULL]
      gobject@cell_metadata = cell_metadata
    }

    gobject <- addCellMetadata(gobject, new_metadata = together[, c('cell_ID', name), with = F],
                               by_column = T, column_cell_ID = 'cell_ID')

    ## update parameters used ##
    parameters_list = gobject@parameters
    number_of_rounds = length(parameters_list)
    update_name = paste0(number_of_rounds,'_sub_cluster')

    # parameters to include
    parameters_list[[update_name]] = c('subclus name' = name,
                                       'resolution ' = resolution,
                                       'k neighbors ' = k_neighbors)

    gobject@parameters = parameters_list

    return(gobject)

  } else {
    return(together)
  }

}


#' @title doLouvainSubCluster_community
#' @name doLouvainSubCluster_community
#' @description subcluster cells using a NN-network and the Louvain community detection algorithm
#' @param gobject giotto object
#' @param name name for new clustering result
#' @param cluster_column cluster column to subcluster
#' @param selected_clusters only do subclustering on these clusters
#' @param hvg_param parameters for calculateHVG
#' @param hvg_min_perc_cells threshold for detection in min percentage of cells
#' @param hvg_mean_expr_det threshold for mean expression level in cells with detection
#' @param use_all_genes_as_hvg forces all genes to be HVG and to be used as input for PCA
#' @param min_nr_of_hvg minimum number of HVG, or all genes will be used as input for PCA
#' @param pca_param parameters for runPCA
#' @param nn_param parameters for parameters for createNearestNetwork
#' @param k_neighbors number of k for createNearestNetwork
#' @param resolution resolution
#' @param python_path specify specific path to python if required
#' @param nn_network_to_use type of NN network to use (kNN vs sNN)
#' @param network_name name of NN network to use
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @param verbose verbose
#' @return giotto object with new subclusters appended to cell metadata
#' @details This function performs subclustering using the Louvain community algorithm on selected clusters.
#' The systematic steps are:
#' \itemize{
#'   \item{1. subset Giotto object}
#'   \item{2. identify highly variable genes}
#'   \item{3. run PCA}
#'   \item{4. create nearest neighbouring network}
#'   \item{5. do Louvain community clustering}
#' }
#' @seealso \code{\link{doLouvainCluster_community}}
#' @keywords internal
doLouvainSubCluster_community = function(gobject,
                                         name = 'sub_louvain_comm_clus',
                                         cluster_column = NULL,
                                         selected_clusters = NULL,
                                         hvg_param = list(reverse_log_scale = T, difference_in_cov = 1, expression_values = 'normalized'),
                                         hvg_min_perc_cells = 5,
                                         hvg_mean_expr_det = 1,
                                         use_all_genes_as_hvg = FALSE,
                                         min_nr_of_hvg = 5,
                                         pca_param = list(expression_values = 'normalized', scale_unit = T),
                                         nn_param = list(dimensions_to_use = 1:20),
                                         k_neighbors = 10,
                                         resolution = 0.5,
                                         python_path = NULL,
                                         nn_network_to_use = 'sNN',
                                         network_name = 'sNN.pca',
                                         return_gobject = TRUE,
                                         verbose = T) {


  iter_list = list()

  cell_metadata = pDataDT(gobject)

  if(is.null(cluster_column)) {
    stop('\n You need to provide a cluster column to subcluster on \n')
  }
  unique_clusters = sort(unique(cell_metadata[[cluster_column]]))

  ## if clusters start at 0, then add +1 for the index ##
  index_offset = ifelse(0 %in% unique_clusters, 1, 0)

  for(cluster in unique_clusters) {

    if(verbose == TRUE) cat('\n start with cluster: ', cluster, '\n')

    ## get subset
    subset_cell_IDs = cell_metadata[get(cluster_column) == cluster][['cell_ID']]
    temp_giotto = subsetGiotto(gobject = gobject, cell_ids = subset_cell_IDs)

    ## if cluster is not selected
    if(!is.null(selected_clusters) & !cluster %in% selected_clusters) {

      temp_cluster = data.table('cell_ID' = subset_cell_IDs, 'tempclus' = 1, 'parent_cluster' = cluster)
      iter_list[[cluster+index_offset]] = temp_cluster

    } else {
      # continue for selected clusters or all clusters if there is no selection

      ## calculate stats
      temp_giotto <- addStatistics(gobject = temp_giotto)

      ## calculate variable genes
      temp_giotto = do.call('calculateHVG', c(gobject = temp_giotto, hvg_param))

      ## get hvg
      gene_metadata = fDataDT(temp_giotto)

      # data.table variables
      hvg = perc_cells = mean_expr_det = NULL

      featgenes     = gene_metadata[hvg == 'yes' & perc_cells >= hvg_min_perc_cells & mean_expr_det >= hvg_mean_expr_det]$gene_ID

      ## catch too low number of hvg
      if(use_all_genes_as_hvg == TRUE) {
        featgenes == gene_metadata$gene_ID
      } else {
        if(verbose == TRUE) cat('\n', length(featgenes), 'highly variable genes have been selected \n')
        if(length(featgenes) <= min_nr_of_hvg) {
          cat('\n too few genes, will continue with all genes instead \n')
          featgenes = gene_metadata$gene_ID
        }
      }

      ## run PCA
      temp_giotto = do.call('runPCA', c(gobject =  temp_giotto, genes_to_use = list(featgenes), pca_param))

      ## nearest neighbor and clustering
      temp_giotto = do.call('createNearestNetwork', c(gobject = temp_giotto, k = k_neighbors, nn_param))

      ## Leiden Cluster
      ## TO DO: expand to all clustering options
      temp_cluster = doLouvainCluster_community(gobject = temp_giotto,
                                                resolution = resolution,
                                                python_path = python_path,
                                                name = 'tempclus',
                                                return_gobject = F)

      # data.table variables
      parent_cluster = NULL

      temp_cluster[, parent_cluster := cluster]

      iter_list[[cluster+index_offset]] = temp_cluster



    }

  }

  together = do.call('rbind', iter_list)

  # data.table variables
  comb = tempclus = NULL

  together[, comb := paste0(parent_cluster,'.',tempclus)]

  # rename with subcluster of original name
  #new_cluster_column = paste0(cluster_column,'_sub')
  setnames(together, 'comb', name)


  if(return_gobject == TRUE) {

    cluster_names = names(gobject@cell_metadata)
    if(name %in% cluster_names) {
      cat('\n ', name, ' has already been used, will be overwritten \n')
      cell_metadata = gobject@cell_metadata
      cell_metadata[, eval(name) := NULL]
      gobject@cell_metadata = cell_metadata
    }

    gobject <- addCellMetadata(gobject, new_metadata = together[, c('cell_ID', name), with = F],
                               by_column = T, column_cell_ID = 'cell_ID')

    ## update parameters used ##
    parameters_list = gobject@parameters
    number_of_rounds = length(parameters_list)
    update_name = paste0(number_of_rounds,'_sub_cluster')

    # parameters to include
    parameters_list[[update_name]] = c('subclus name' = name,
                                       'resolution ' = resolution,
                                       'k neighbors ' = k_neighbors)

    gobject@parameters = parameters_list

    return(gobject)

  } else {
    return(together)
  }

}




#' @title doLouvainSubCluster_multinet
#' @name doLouvainSubCluster_multinet
#' @description subcluster cells using a NN-network and the Louvain multinet detection algorithm
#' @param gobject giotto object
#' @param name name for new clustering result
#' @param cluster_column cluster column to subcluster
#' @param selected_clusters only do subclustering on these clusters
#' @param hvg_param parameters for calculateHVG
#' @param hvg_min_perc_cells threshold for detection in min percentage of cells
#' @param hvg_mean_expr_det threshold for mean expression level in cells with detection
#' @param use_all_genes_as_hvg forces all genes to be HVG and to be used as input for PCA
#' @param min_nr_of_hvg minimum number of HVG, or all genes will be used as input for PCA
#' @param pca_param parameters for runPCA
#' @param nn_param parameters for parameters for createNearestNetwork
#' @param k_neighbors number of k for createNearestNetwork
#' @param gamma gamma
#' @param omega omega
#' @param nn_network_to_use type of NN network to use (kNN vs sNN)
#' @param network_name name of NN network to use
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @param verbose verbose
#' @return giotto object with new subclusters appended to cell metadata
#' @details This function performs subclustering using the Louvain multinet algorithm on selected clusters.
#' The systematic steps are:
#' \itemize{
#'   \item{1. subset Giotto object}
#'   \item{2. identify highly variable genes}
#'   \item{3. run PCA}
#'   \item{4. create nearest neighbouring network}
#'   \item{5. do Louvain multinet clustering}
#' }
#' @seealso \code{\link{doLouvainCluster_multinet}}
#' @keywords internal
doLouvainSubCluster_multinet =  function(gobject,
                                         name = 'sub_louvain_mult_clus',
                                         cluster_column = NULL,
                                         selected_clusters = NULL,
                                         hvg_param = list(reverse_log_scale = T, difference_in_cov = 1, expression_values = 'normalized'),
                                         hvg_min_perc_cells = 5,
                                         hvg_mean_expr_det = 1,
                                         use_all_genes_as_hvg = FALSE,
                                         min_nr_of_hvg = 5,
                                         pca_param = list(expression_values = 'normalized', scale_unit = T),
                                         nn_param = list(dimensions_to_use = 1:20),
                                         k_neighbors = 10,
                                         gamma = 1,
                                         omega = 1,
                                         nn_network_to_use = 'sNN',
                                         network_name = 'sNN.pca',
                                         return_gobject = TRUE,
                                         verbose = T) {


  if("multinet" %in% rownames(installed.packages()) == FALSE) {
    stop("\n package 'multinet' is not yet installed \n",
         "To install: \n",
         "install.packages('multinet')"
    )
  }

  iter_list = list()

  cell_metadata = pDataDT(gobject)

  if(is.null(cluster_column)) {
    stop('\n You need to provide a cluster column to subcluster on \n')
  }
  unique_clusters = sort(unique(cell_metadata[[cluster_column]]))

  ## if clusters start at 0, then add +1 for the index ##
  index_offset = ifelse(0 %in% unique_clusters, 1, 0)


  # data.table variables
  hvg = perc_cells = mean_expr_det = parent_cluster = cell_ID = comb = tempclus = NULL


  for(cluster in unique_clusters) {

    if(verbose == TRUE) cat('\n start with cluster: ', cluster, '\n')

    ## get subset
    subset_cell_IDs = cell_metadata[get(cluster_column) == cluster][['cell_ID']]
    temp_giotto = subsetGiotto(gobject = gobject, cell_ids = subset_cell_IDs)

    ## if cluster is not selected
    if(!is.null(selected_clusters) & !cluster %in% selected_clusters) {

      temp_cluster = data.table('cell_ID' = subset_cell_IDs, 'tempclus' = 1, 'parent_cluster' = cluster)
      iter_list[[cluster+index_offset]] = temp_cluster

    } else {
      # continue for selected clusters or all clusters if there is no selection

      ## calculate stats
      temp_giotto <- addStatistics(gobject = temp_giotto)

      ## calculate variable genes
      temp_giotto = do.call('calculateHVG', c(gobject = temp_giotto, hvg_param))

      ## get hvg
      gene_metadata = fDataDT(temp_giotto)
      featgenes     = gene_metadata[hvg == 'yes' & perc_cells >= hvg_min_perc_cells & mean_expr_det >= hvg_mean_expr_det]$gene_ID

      ## catch too low number of hvg
      if(use_all_genes_as_hvg == TRUE) {
        featgenes == gene_metadata$gene_ID
      } else {
        if(verbose == TRUE) cat('\n', length(featgenes), 'highly variable genes have been selected \n')
        if(length(featgenes) <= min_nr_of_hvg) {
          cat('\n too few genes, will continue with all genes instead \n')
          featgenes = gene_metadata$gene_ID
        }
      }

      ## run PCA
      temp_giotto = do.call('runPCA', c(gobject =  temp_giotto, genes_to_use = list(featgenes), pca_param))

      ## nearest neighbor and clustering
      temp_giotto = do.call('createNearestNetwork', c(gobject = temp_giotto, k = k_neighbors, nn_param))

      ## Leiden Cluster
      ## TO DO: expand to all clustering options
      temp_cluster = doLouvainCluster_multinet(gobject = temp_giotto,
                                               gamma = gamma,
                                               omega = omega,
                                               name = 'tempclus',
                                               return_gobject = F)

      temp_cluster[, parent_cluster := cluster]
      temp_cluster = temp_cluster[,.(cell_ID, tempclus, parent_cluster)]

      iter_list[[cluster+index_offset]] = temp_cluster



    }

  }

  together = do.call('rbind', iter_list)
  together[, comb := paste0(parent_cluster,'.',tempclus)]

  # rename with subcluster of original name
  #new_cluster_column = paste0(cluster_column,'_sub')
  setnames(together, 'comb', name)


  if(return_gobject == TRUE) {

    cluster_names = names(gobject@cell_metadata)
    if(name %in% cluster_names) {
      cat('\n ', name, ' has already been used, will be overwritten \n')
      cell_metadata = gobject@cell_metadata
      cell_metadata[, eval(name) := NULL]
      gobject@cell_metadata = cell_metadata
    }

    gobject <- addCellMetadata(gobject, new_metadata = together[, c('cell_ID', name), with = F],
                               by_column = T, column_cell_ID = 'cell_ID')

    ## update parameters used ##
    parameters_list = gobject@parameters
    number_of_rounds = length(parameters_list)
    update_name = paste0(number_of_rounds,'_sub_cluster')

    # parameters to include
    parameters_list[[update_name]] = c('subclus name' = name,
                                       'k neighbors ' = k_neighbors,
                                       'gamma' = gamma,
                                       'omega' = omega)

    gobject@parameters = parameters_list

    return(gobject)

  } else {
    return(together)
  }

}



#' @title doLouvainSubCluster
#' @name doLouvainSubCluster
#' @description subcluster cells using a NN-network and the Louvain algorithm
#' @param gobject giotto object
#' @param name name for new clustering result
#' @param version version of Louvain algorithm to use
#' @param cluster_column cluster column to subcluster
#' @param selected_clusters only do subclustering on these clusters
#' @param hvg_param parameters for calculateHVG
#' @param hvg_min_perc_cells threshold for detection in min percentage of cells
#' @param hvg_mean_expr_det threshold for mean expression level in cells with detection
#' @param use_all_genes_as_hvg forces all genes to be HVG and to be used as input for PCA
#' @param min_nr_of_hvg minimum number of HVG, or all genes will be used as input for PCA
#' @param pca_param parameters for runPCA
#' @param nn_param parameters for parameters for createNearestNetwork
#' @param k_neighbors number of k for createNearestNetwork
#' @param resolution resolution for community algorithm
#' @param gamma gamma
#' @param omega omega
#' @param python_path specify specific path to python if required
#' @param nn_network_to_use type of NN network to use (kNN vs sNN)
#' @param network_name name of NN network to use
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @param verbose verbose
#' @return giotto object with new subclusters appended to cell metadata
#' @details This function performs subclustering using the Louvain algorithm on selected clusters.
#' The systematic steps are:
#' \itemize{
#'   \item{1. subset Giotto object}
#'   \item{2. identify highly variable genes}
#'   \item{3. run PCA}
#'   \item{4. create nearest neighbouring network}
#'   \item{5. do Louvain clustering}
#' }
#' @seealso \code{\link{doLouvainCluster_multinet}} and \code{\link{doLouvainCluster_community}}
#' @export
doLouvainSubCluster =  function(gobject,
                                name = 'sub_louvain_clus',
                                version = c('community', 'multinet'),
                                cluster_column = NULL,
                                selected_clusters = NULL,
                                hvg_param = list(reverse_log_scale = T, difference_in_cov = 1, expression_values = 'normalized'),
                                hvg_min_perc_cells = 5,
                                hvg_mean_expr_det = 1,
                                use_all_genes_as_hvg = FALSE,
                                min_nr_of_hvg = 5,
                                pca_param = list(expression_values = 'normalized', scale_unit = T),
                                nn_param = list(dimensions_to_use = 1:20),
                                k_neighbors = 10,
                                resolution = 0.5,
                                gamma = 1,
                                omega = 1,
                                python_path = NULL,
                                nn_network_to_use = 'sNN',
                                network_name = 'sNN.pca',
                                return_gobject = TRUE,
                                verbose = T) {

  ## louvain clustering version to use
  version = match.arg(version, c('community', 'multinet'))


  # python community implementation
  if(version == 'community') {

    result = doLouvainSubCluster_community(gobject = gobject,
                                           cluster_column = cluster_column,
                                           selected_clusters = selected_clusters,
                                           hvg_param = hvg_param,
                                           hvg_mean_expr_det = hvg_mean_expr_det,
                                           pca_param = pca_param,
                                           nn_param = nn_param,
                                           k_neighbors = k_neighbors,
                                           resolution = resolution,
                                           python_path = python_path,
                                           nn_network_to_use = nn_network_to_use,
                                           network_name = network_name,
                                           name = name,
                                           return_gobject = return_gobject,
                                           verbose = verbose)

  } else if(version == 'multinet') {

    result = doLouvainSubCluster_multinet(gobject = gobject,
                                          cluster_column = cluster_column,
                                          selected_clusters = selected_clusters,
                                          hvg_param = hvg_param,
                                          hvg_mean_expr_det = hvg_mean_expr_det,
                                          pca_param = pca_param,
                                          nn_param = nn_param,
                                          k_neighbors = k_neighbors,
                                          gamma = gamma,
                                          omega = omega,
                                          nn_network_to_use = nn_network_to_use,
                                          network_name = network_name,
                                          name = name,
                                          return_gobject = return_gobject,
                                          verbose = verbose)

  }

  return(result)

}





#' @title subClusterCells
#' @name subClusterCells
#' @description subcluster cells
#' @param gobject giotto object
#' @param name name for new clustering result
#' @param cluster_method clustering method to use
#' @param cluster_column cluster column to subcluster
#' @param selected_clusters only do subclustering on these clusters
#' @param hvg_param parameters for calculateHVG
#' @param hvg_min_perc_cells threshold for detection in min percentage of cells
#' @param hvg_mean_expr_det threshold for mean expression level in cells with detection
#' @param use_all_genes_as_hvg forces all genes to be HVG and to be used as input for PCA
#' @param min_nr_of_hvg minimum number of HVG, or all genes will be used as input for PCA
#' @param pca_param parameters for runPCA
#' @param nn_param parameters for parameters for createNearestNetwork
#' @param k_neighbors number of k for createNearestNetwork
#' @param resolution resolution
#' @param n_iterations number of interations to run the Leiden algorithm.
#' @param gamma gamma
#' @param omega omega
#' @param python_path specify specific path to python if required
#' @param nn_network_to_use type of NN network to use (kNN vs sNN)
#' @param network_name name of NN network to use
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @param verbose verbose
#' @return giotto object with new subclusters appended to cell metadata
#' @details This function performs subclustering on selected clusters.
#' The systematic steps are:
#' \itemize{
#'   \item{1. subset Giotto object}
#'   \item{2. identify highly variable genes}
#'   \item{3. run PCA}
#'   \item{4. create nearest neighbouring network}
#'   \item{5. do clustering}
#' }
#' @seealso \code{\link{doLouvainCluster_multinet}}, \code{\link{doLouvainCluster_community}}
#' and  @seealso \code{\link{doLeidenCluster}}
#' @export
subClusterCells <- function(gobject,
                            name = 'sub_clus',
                            cluster_method = c('leiden',
                                               'louvain_community',
                                               'louvain_multinet'),
                            cluster_column = NULL,
                            selected_clusters = NULL,
                            hvg_param = list(reverse_log_scale = T, difference_in_cov = 1, expression_values = 'normalized'),
                            hvg_min_perc_cells = 5,
                            hvg_mean_expr_det = 1,
                            use_all_genes_as_hvg = FALSE,
                            min_nr_of_hvg = 5,
                            pca_param = list(expression_values = 'normalized', scale_unit = T),
                            nn_param = list(dimensions_to_use = 1:20),
                            k_neighbors = 10,
                            resolution = 1,
                            n_iterations = 1000,
                            gamma = 1,
                            omega = 1,
                            python_path = NULL,
                            nn_network_to_use = 'sNN',
                            network_name = 'sNN.pca',
                            return_gobject = TRUE,
                            verbose = T) {

  ## select cluster method
  cluster_method = match.arg(arg = cluster_method, choices = c('leiden',
                                                               'louvain_community',
                                                               'louvain_multinet'))


  if(cluster_method == 'leiden') {

    result = doLeidenSubCluster(gobject = gobject,
                                cluster_column = cluster_column,
                                selected_clusters = selected_clusters,
                                hvg_param = hvg_param,
                                hvg_min_perc_cells = hvg_min_perc_cells,
                                hvg_mean_expr_det = hvg_mean_expr_det,
                                use_all_genes_as_hvg = use_all_genes_as_hvg,
                                min_nr_of_hvg = min_nr_of_hvg,
                                pca_param = pca_param,
                                nn_param = nn_param,
                                k_neighbors = k_neighbors,
                                resolution = resolution,
                                n_iterations = n_iterations,
                                python_path = python_path,
                                nn_network_to_use = nn_network_to_use,
                                network_name = network_name,
                                name = name,
                                return_gobject = return_gobject,
                                verbose = verbose)

  } else if(cluster_method == 'louvain_community') {

    result = doLouvainSubCluster_community(gobject = gobject,
                                        cluster_column = cluster_column,
                                        selected_clusters = selected_clusters,
                                        hvg_param = hvg_param,
                                        hvg_min_perc_cells = hvg_min_perc_cells,
                                        hvg_mean_expr_det = hvg_mean_expr_det,
                                        use_all_genes_as_hvg = use_all_genes_as_hvg,
                                        min_nr_of_hvg = min_nr_of_hvg,
                                        pca_param = pca_param,
                                        nn_param = nn_param,
                                        k_neighbors = k_neighbors,
                                        resolution = resolution,
                                        python_path = python_path,
                                        nn_network_to_use = nn_network_to_use,
                                        network_name = network_name,
                                        name = name,
                                        return_gobject = return_gobject,
                                        verbose = verbose)

  } else if(cluster_method == 'louvain_multinet') {

    result = doLouvainSubCluster_multinet(gobject = gobject,
                                       cluster_column = cluster_column,
                                       selected_clusters = selected_clusters,
                                       hvg_param = hvg_param,
                                       hvg_min_perc_cells = hvg_min_perc_cells,
                                       hvg_mean_expr_det = hvg_mean_expr_det,
                                       use_all_genes_as_hvg = use_all_genes_as_hvg,
                                       min_nr_of_hvg = min_nr_of_hvg,
                                       pca_param = pca_param,
                                       nn_param = nn_param,
                                       k_neighbors = k_neighbors,
                                       gamma = gamma,
                                       omega = omega,
                                       nn_network_to_use = nn_network_to_use,
                                       network_name = network_name,
                                       name = name,
                                       return_gobject = return_gobject,
                                       verbose = verbose)

  }

  return(result)

}








#' @title getClusterSimilarity
#' @name getClusterSimilarity
#' @description Creates data.table with pairwise correlation scores between each cluster.
#' @param gobject giotto object
#' @param expression_values expression values to use
#' @param cluster_column name of column to use for clusters
#' @param cor correlation score to calculate distance
#' @return data.table
#' @details Creates data.table with pairwise correlation scores between each cluster and
#' the group size (# of cells) for each cluster. This information can be used together
#' with mergeClusters to combine very similar or small clusters into bigger clusters.
#' @export
#' @examples
#'
#' data("mini_giotto_single_cell")
#'
#' cluster_similarities = getClusterSimilarity(mini_giotto_single_cell,
#'                                             cluster_column = 'leiden_clus')
#'
getClusterSimilarity <- function(gobject,
                                 expression_values = c('normalized', 'scaled', 'custom'),
                                 cluster_column,
                                 cor = c('pearson', 'spearman')) {


  # data.table variables
  group1 = group2 = unified_group = value = NULL

  cor = match.arg(cor, c('pearson', 'spearman'))
  values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))

  metadata = pDataDT(gobject)

  # get clustersize
  clustersize = metadata[, .N, by = cluster_column]
  colnames(clustersize) = c('clusters', 'size')

  # data.table variables
  clusters = NULL

  clustersize[, clusters := as.character(clusters)]

  # scores per cluster
  metatable = calculateMetaTable(gobject = gobject, expression_values = values, metadata_cols = cluster_column)
  dcast_metatable = data.table::dcast.data.table(metatable, formula = variable~uniq_ID, value.var = 'value')
  testmatrix = dt_to_matrix(x = dcast_metatable)

  # correlation matrix
  cormatrix = cor_giotto(x = testmatrix, method = cor)
  cor_table = data.table::as.data.table(reshape2::melt(cormatrix))
  data.table::setnames(cor_table, old = c('Var1', 'Var2'), c('group1', 'group2'))
  cor_table[, c('group1', 'group2') := list(as.character(group1), as.character(group2))]
  cor_table[, unified_group := paste(sort(c(group1, group2)), collapse = '--'), by = 1:nrow(cor_table)]
  cor_table = cor_table[!duplicated(cor_table[, .(value, unified_group)])]

  cor_table = merge(cor_table, by.x = 'group1', clustersize, by.y = 'clusters')
  setnames(cor_table, 'size', 'group1_size')
  cor_table = merge(cor_table, by.x = 'group2', clustersize, by.y = 'clusters')
  setnames(cor_table, 'size', 'group2_size')

  return(cor_table)


}




#' @title mergeClusters
#' @name mergeClusters
#' @description Merge selected clusters based on pairwise correlation scores and size of cluster.
#' @param gobject giotto object
#' @param expression_values expression values to use
#' @param cluster_column name of column to use for clusters
#' @param cor correlation score to calculate distance
#' @param new_cluster_name new name for merged clusters
#' @param min_cor_score min correlation score to merge pairwise clusters
#' @param max_group_size max cluster size that can be merged
#' @param force_min_group_size size of clusters that will be merged with their most similar neighbor(s)
#' @param max_sim_clusters maximum number of clusters to potentially merge to reach force_min_group_size
#' @param return_gobject return giotto object
#' @param verbose be verbose
#' @return Giotto object
#' @details Merge selected clusters based on pairwise correlation scores and size of cluster.
#' To avoid large clusters to merge the max_group_size can be lowered. Small clusters can
#' be forcibly merged with their most similar pairwise cluster by adjusting the
#' force_min_group_size parameter. Clusters smaller than this value will be merged
#' independent on the provided min_cor_score value. The force_min_group_size might not always
#' be reached if clusters have already been merged before \cr
#' A giotto object is returned by default, if FALSE then the merging vector will be returned.
#' @export
#' @examples
#'
#' data("mini_giotto_single_cell")
#'
#' pDataDT(mini_giotto_single_cell)
#' mini_giotto_single_cell = mergeClusters(mini_giotto_single_cell,
#'                                         cluster_column = 'leiden_clus',
#'                                         min_cor_score = 0.7,
#'                                         force_min_group_size = 4)
#' pDataDT(mini_giotto_single_cell)
#' plotUMAP_2D(mini_giotto_single_cell, cell_color = 'merged_cluster', point_size = 3)
#'
mergeClusters <- function(gobject,
                          expression_values = c('normalized', 'scaled', 'custom'),
                          cluster_column,
                          cor = c('pearson', 'spearman'),
                          new_cluster_name = 'merged_cluster',
                          min_cor_score = 0.8,
                          max_group_size = 20,
                          force_min_group_size = 10,
                          max_sim_clusters = 10,
                          return_gobject = TRUE,
                          verbose = TRUE) {

  # expression values to be used
  values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))

  # correlation score to be used
  cor = match.arg(cor, c('pearson', 'spearman'))

  # calculate similarity data.table
  similarityDT = getClusterSimilarity(gobject = gobject,
                                      expression_values = values,
                                      cluster_column = cluster_column,
                                      cor = cor)

  ## get clusters that can be merged
  # 1. clusters with high correlation

  # data.table variables
  group1 = group2 = group1_size = value = cumsum_val = group2_size = min_reached = cumsum_reached = NULL

  filter_set_first = similarityDT[group1 != group2][group1_size < max_group_size][value >= min_cor_score]

  # 2. small clusters
  minimum_set = similarityDT[group1 != group2][group1_size < force_min_group_size][order(-value)][, head(.SD, max_sim_clusters), by = group1]

  # 2.1 take all clusters necessary to reach force_min_group_size
  minimum_set[, cumsum_val := cumsum(group2_size) + group1_size, by = group1]
  minimum_set[, min_reached := ifelse(cumsum_val > max_group_size, 1, 0)]
  minimum_set[, cumsum_reached := cumsum(min_reached), by = group1]
  minimum_set = minimum_set[cumsum_reached <= 1]
  minimum_set[, c('cumsum_val', 'min_reached', 'cumsum_reached') := NULL]

  filter_set = unique(do.call('rbind', list(filter_set_first, minimum_set)))

  ## get list of correlated groups
  finallist = list()
  start_i = 1
  for(row in 1:nrow(filter_set)) {

    first_clus = filter_set[row][['group1']]
    second_clus = filter_set[row][['group2']]

    res = lapply(finallist, function(x) {any(x %in% c(first_clus, second_clus))})

    if(all(res == F)) {

      #print('not in list yet')
      finallist[[start_i]] = c(first_clus, second_clus)
      start_i = start_i + 1

    } else if(length(res[res == T]) == 2) {

      #print('do nothing')
      NULL

    } else {

      #print('already in list')
      who = which(res == TRUE)[[1]]
      finallist[[who]] = unique(c(finallist[[who]], first_clus, second_clus))

    }

  }


  ## update metadata
  metadata = data.table::copy(pDataDT(gobject))

  finalvec = NULL
  for(ll in 1:length(finallist)) {
    tempvec = finallist[[ll]]; names(tempvec) = rep(paste0('m_', ll), length(tempvec))
    finalvec = c(finalvec, tempvec)
  }

  metadata[, eval(new_cluster_name) := ifelse(as.character(get(cluster_column)) %in% finalvec,
                                              names(finalvec[finalvec == as.character(get(cluster_column))]),
                                              as.character(get(cluster_column))), by = 1:nrow(metadata)]




  if(return_gobject == TRUE) {

    cluster_names = names(gobject@cell_metadata)
    if(new_cluster_name %in% cluster_names) {
      cat('\n ', new_cluster_name, ' has already been used, will be overwritten \n')
      cell_metadata = gobject@cell_metadata
      cell_metadata[, eval(new_cluster_name) := NULL]
      gobject@cell_metadata = cell_metadata
    }

    gobject = addCellMetadata(gobject = gobject,
                              new_metadata = metadata[, c('cell_ID', new_cluster_name), with = F],
                              by_column = T, column_cell_ID = 'cell_ID')


    ## update parameters used ##
    parameters_list = gobject@parameters
    number_of_rounds = length(parameters_list)
    update_name = paste0(number_of_rounds,'_merge_cluster')

    parameters_list[[update_name]] = c('expression values' = values,
                                       'cluster column to merge' = cluster_column,
                                       'cor score' = cor,
                                       'minimum cor score' = min_cor_score,
                                       'max group size to merge' = max_group_size,
                                       'group size that will be forcibly merged' = force_min_group_size)

    gobject@parameters = parameters_list
    return(gobject)



  } else {

    return(list(mergevector = finalvec, metadata = metadata))

  }
}





#' @title split_dendrogram_in_two
#' @name split_dendrogram_in_two
#' @description Merge selected clusters based on pairwise correlation scores and size of cluster.
#' @param dend dendrogram object
#' @return list of two dendrograms and height of node
#' @keywords internal
split_dendrogram_in_two = function(dend) {

  top_height = attributes(dend)$height
  divided_leaves_labels = dendextend::cut_lower_fun(dend, h = top_height)

  # this works for both numericala nd character leave names
  all_leaves = dendextend::get_leaves_attr(dend = dend, attribute = 'label')
  selected_labels_ind_1 = all_leaves %in% divided_leaves_labels[[1]]
  selected_labels_ind_2 = all_leaves %in% divided_leaves_labels[[2]]
  numerical_leaves = unlist(dend)
  names(numerical_leaves) = all_leaves

  dend_1 = dendextend::find_dendrogram(dend = dend, selected_labels = names(numerical_leaves[selected_labels_ind_1]))
  dend_2 = dendextend::find_dendrogram(dend = dend, selected_labels = names(numerical_leaves[selected_labels_ind_2]))

  #dend_1 = dendextend::find_dendrogram(dend = dend, selected_labels = divided_leaves_labels[[1]])
  #dend_2 = dendextend::find_dendrogram(dend = dend, selected_labels = divided_leaves_labels[[2]])

  return(list(theight = top_height, dend1 =  dend_1, dend2 = dend_2))
}



#' @title node_clusters
#' @name node_clusters
#' @description Merge selected clusters based on pairwise correlation scores and size of cluster.
#' @param hclus_obj hclus object
#' @param verbose be verbose
#' @return list of splitted dendrogram nodes from high to low node height
#' @keywords internal
node_clusters = function(hclus_obj, verbose = TRUE) {

  heights = sort(hclus_obj[['height']], decreasing = T)
  mydend = stats::as.dendrogram(hclus_obj)


  result_list = list()
  j = 1

  dend_list = list()
  i = 1
  dend_list[[i]] = mydend

  ## create split at each height ##
  for(n_height in heights) {

    if(verbose == TRUE) cat('height ', n_height, '\n')

    # only use dendrogram objects
    ind = lapply(dend_list, FUN = function(x) class(x) == 'dendrogram')
    dend_list = dend_list[unlist(ind)]

    # check which heights are available
    available_h = as.numeric(unlist(lapply(dend_list, FUN = function(x) attributes(x)$height)))

    # get dendrogram associated with height and split in two
    # select_dend_ind = which(available_h == n_height)
    select_dend_ind = which.min(abs(available_h - n_height))
    select_dend = dend_list[[select_dend_ind]]
    tempres = split_dendrogram_in_two(dend = select_dend)

    # find leave labels
    toph = tempres[[1]]
    first_group = dendextend::get_leaves_attr(tempres[[2]], attribute = 'label')
    second_group = dendextend::get_leaves_attr(tempres[[3]], attribute = 'label')

    result_list[[j]] = list('height' = toph, 'first' = first_group, 'sec' = second_group)
    j = j+1



    ## add dendrograms to list
    ind = lapply(tempres, FUN = function(x) class(x) == 'dendrogram')
    tempres_dend = tempres[unlist(ind)]

    dend_list[[i+1]] = tempres_dend[[1]]
    dend_list[[i+2]] = tempres_dend[[2]]

    i = i+2


  }

  return(list(dend_list, result_list))

}



#' @title getDendrogramSplits
#' @name getDendrogramSplits
#' @description Split dendrogram at each node and keep the leave (label) information..
#' @param gobject giotto object
#' @param expression_values expression values to use
#' @param cluster_column name of column to use for clusters
#' @param cor correlation score to calculate distance
#' @param distance distance method to use for hierarchical clustering
#' @param h height of horizontal lines to plot
#' @param h_color color of horizontal lines
#' @param show_dend show dendrogram
#' @param verbose be verbose
#' @return data.table object
#' @details Creates a data.table with three columns and each row represents a node in the
#' dendrogram. For each node the height of the node is given together with the two
#' subdendrograms. This information can be used to determine in a hierarchical manner
#' differentially expressed marker genes at each node.
#' @export
#' @examples
#' \dontrun{
#' data("mini_giotto_single_cell")
#'
#' splits = getDendrogramSplits(mini_giotto_single_cell, cluster_column = 'leiden_clus')
#' }

getDendrogramSplits = function(gobject,
                               expression_values = c('normalized', 'scaled', 'custom'),
                               cluster_column,
                               cor = c('pearson', 'spearman'),
                               distance = 'ward.D',
                               h = NULL,
                               h_color = 'red',
                               show_dend = TRUE,
                               verbose = TRUE) {


  # package check for dendextend
  package_check(pkg_name = "dendextend", repository = "CRAN")

  # data.table variables
  nodeID = NULL

  cor = match.arg(cor, c('pearson', 'spearman'))
  values = match.arg(expression_values, c('normalized', 'scaled', 'custom'))

  # create average expression matrix per cluster
  metatable = calculateMetaTable(gobject = gobject, expression_values = values, metadata_cols = cluster_column)
  dcast_metatable = data.table::dcast.data.table(metatable, formula = variable~uniq_ID, value.var = 'value')
  testmatrix = dt_to_matrix(x = dcast_metatable)

  # correlation
  cormatrix = cor_giotto(x = testmatrix, method = cor)
  cordist = stats::as.dist(1 - cormatrix, diag = T, upper = T)
  corclus = stats::hclust(d = cordist, method = distance)

  cordend = stats::as.dendrogram(object = corclus)

  ## print dendrogram ##
  if(show_dend == TRUE) {
    # plot dendrogram
    graphics::plot(cordend)

    # add horizontal lines
    if(!is.null(h)) {
      graphics::abline(h = h, col = h_color)
    }
  }


  splitList = node_clusters(hclus_obj = corclus, verbose = verbose)

  splitDT = data.table::as.data.table(t_giotto(data.table::as.data.table(splitList[[2]])))
  colnames(splitDT) = c('node_h', 'tree_1', 'tree_2')
  splitDT[, nodeID := paste0('node_', 1:.N)]

  return(splitDT)

}



