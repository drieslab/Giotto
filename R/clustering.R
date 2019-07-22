


#' @title clusterCells
#' @name clusterCells
#' @description cluster cells using a NN-network and community detection algorithms
#' @param gobject giotto object
#' @param nn_network_to_use type of NN network to use (kNN vs sNN)
#' @param network_name name of NN network to use
#' @param cluster_method community cluster method to use
#' @param pyth_louv_resolution python louvain param: resolution
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
#' @param name name for cluster
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @param set_seed set seed
#' @param seed_number number for seed
#' @param ... additional parameters
#' @return giotto object appended with new cluster
#' @details Description of different clustering methods.
#' @export
#' @examples
#'     clusterCells(gobject)
clusterCells <- function(gobject,
                         nn_network_to_use = 'sNN',
                         network_name = 'sNN.pca',
                         cluster_method = c('python_leiden', 'python_louvain', 'louvain', 'randomwalk', 'sNNclust'),

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

                         name = 'plouvain_clus',
                         return_gobject = TRUE,
                         set_seed = T,
                         seed_number = 1234,
                         ...) {


  cluster_method = match.arg(arg = cluster_method, choices = c('python_leiden', 'python_louvain', 'louvain', 'randomwalk', 'sNNclust'))
  cell_ID_vec = gobject@cell_ID

  ## select network to use
  if(is.null(nn_network_to_use) | is.null(network_name)) {
    cat('you need to select network type: knn or snn \n
        and you need to select the network name you created\n')
  } else if(nn_network_to_use == 'sNN' & cluster_method == 'sNNclust') {
    stop('\n sNNclust can already be used with a kNN, not a pre-existing sNN \n')
  } else {
    igraph_object = gobject@nn_network[[nn_network_to_use]][[network_name]][['igraph']]
  }


  # perform clustering
  if(cluster_method == 'louvain') {
    ## louvain clustering using the R library 'multinet' ##

    # create mlnetworkobject
    mln_object <- multinet::ml_empty()
    multinet::add_actors_ml(mlnetwork = mln_object, actors = V(igraph_object))
    multinet::add_igraph_layer_ml(mlnetwork = mln_object, g = igraph_object, name = name)

    # start seed
    if(set_seed == TRUE) {
      set.seed(seed = seed_number)
    }

    louvain_clusters = multinet::glouvain_ml(mlnetwork = mln_object, gamma = louvain_gamma, omega = louvain_omega, ...)
    ident_clusters_DT = data.table::as.data.table(louvain_clusters)
    ident_clusters_DT[, cell_ID := actor]
    setnames(ident_clusters_DT, 'cid', name)

    # exit seed
    if(set_seed == TRUE) {
      set.seed(Sys.time())
    }




  } else if(cluster_method == 'python_leiden') {
    ## louvain clustering using the python module 'leidenalg' ##

    pyth_leid_part_type = match.arg(pyth_leid_part_type,
                                    choices = c('RBConfigurationVertexPartition', 'ModularityVertexPartition'))

    ## check or make paths
    # python path
    if(is.null(python_path)) {
      python_path = system('which python', intern = T)
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

    # add weight for edges or set to 1 for all
    if(!is.null(pyth_louv_weight_col)) {
      if(!pyth_louv_weight_col %in% colnames(network_edge_dt)) {
        stop('\n weight column is not an igraph attribute \n')
      } else {
        # weight is defined by attribute of igraph object
        network_edge_dt = network_edge_dt[,c('from', 'to', pyth_louv_weight_col), with = F]
        setnames(network_edge_dt, pyth_louv_weight_col, 'weight')
      }
    } else {
      # weight is the same
      network_edge_dt = network_edge_dt[,c('from', 'to'), with = F]
      network_edge_dt[, weight := 1]
    }


    ## do python leiden clustering
    reticulate::py_set_seed(seed = seed_number, disable_hash_randomization = TRUE)
    pyth_leid_result = python_leiden(df = network_edge_dt,
                                     partition_type = pyth_leid_part_type,
                                     initial_membership = pyth_leid_init_memb,
                                     weights = pyth_leid_weight_col,
                                     n_iterations = pyth_leid_iterations,
                                     seed = seed_number,
                                     resolution_parameter = pyth_leid_resolution)

    ident_clusters_DT = data.table::data.table(cell_ID = pyth_leid_result[[1]], 'name' = pyth_leid_result[[2]])
    data.table::setnames(ident_clusters_DT, 'name', name)




  } else if(cluster_method == 'python_louvain') {
    ## louvain clustering using the python module 'community' ##

    ## check or make paths
    # python path
    if(is.null(python_path)) {
      python_path = system('which python', intern = T)
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

    if(!is.null(pyth_louv_weight_col)) {

      if(!pyth_louv_weight_col %in% colnames(network_edge_dt)) {
        stop('\n weight column is not an igraph attribute \n')
      } else {
        # weight is defined by attribute of igraph object
        network_edge_dt = network_edge_dt[,c('from', 'to', pyth_louv_weight_col), with = F]
        setnames(network_edge_dt, pyth_louv_weight_col, 'weight')
      }
    } else {
      # weight is the same
      network_edge_dt = network_edge_dt[,c('from', 'to'), with = F]
      network_edge_dt[, weight := 1]
    }

    # do python louvain clustering
    #print(pyth_louv_resolution)
    #print(seed_number)

    if(python_louv_random == FALSE) {
      reticulate::py_set_seed(seed = seed_number, disable_hash_randomization = TRUE)
      pyth_louv_result = python_louvain(df = network_edge_dt, resolution = pyth_louv_resolution, randomize = F)
    } else {
      reticulate::py_set_seed(seed = seed_number, disable_hash_randomization = TRUE)
      pyth_louv_result = python_louvain(df = network_edge_dt, resolution = pyth_louv_resolution, random_state = seed_number)
    }
    ident_clusters_DT = data.table::data.table(cell_ID = rownames(pyth_louv_result), 'name' = pyth_louv_result[[1]])
    data.table::setnames(ident_clusters_DT, 'name', name)



  } else if(cluster_method == 'random_walk') {
    ## random walk from the R package 'igraph'

    # start seed
    if(set_seed == TRUE) {
      set.seed(seed = seed_number)
    }

    randomwalk_clusters <- igraph::cluster_walktrap(graph = igraph_object, steps = walk_steps, weights = walk_weights)
    randomwalk_clusters <- as.factor(igraph::cut_at(communities = randomwalk_clusters, no =  walk_clusters))

    ident_clusters_DT <- data.table::data.table('cell_ID' = V(igraph_object)$name, 'name' = randomwalk_clusters)
    # TODO: remove if tested
    #ident_clusters_DT <- data.table('cell_ID' = cell_ID_vec, 'name' = randomwalk_clusters)
    setnames(ident_clusters_DT, 'name', name)

    # exit seed
    if(set_seed == TRUE) {
      set.seed(Sys.time())
    }




  } else if(cluster_method == 'sNNclust') {
  ## SNNclust from the R package 'dbscan'

    if(nn_network_to_use == 'sNN') {
      stop('\n sNNclust can only be used with kNN-network \n')
    }


    # start seed
    if(set_seed == TRUE) {
      set.seed(seed = seed_number)
    }

    # TODO: NOT TESTED YET #
    igraph_DT = data.table::as.data.table(igraph::as_data_frame(igraph_object, what = 'edges'))
    igraph_DT = igraph_DT[order(from)]

    cell_id_numeric = unique(x = c(igraph_DT$from, igraph_DT$to))
    names(cell_id_numeric) <- 1:length(cell_id_numeric)
    igraph_DT[, from_T := as.numeric(names(cell_id_numeric[cell_id_numeric == from])), by = 1:nrow(igraph_DT)]
    igraph_DT[, to_T := as.numeric(names(cell_id_numeric[cell_id_numeric == to])), by = 1:nrow(igraph_DT)]
    temp_igraph_DT = igraph_DT[,.(from_T, to_T, weight, distance)]
    data.table::setnames(temp_igraph_DT, old = c('from_T', 'to_T'), new = c('from', 'to'))

    kNN_object = nnDT_to_kNN(nnDT = temp_igraph_DT)
    sNN_clusters = dbscan::sNNclust(x = kNN_object, k = sNNclust_k, eps = sNNclust_eps,
                                    minPts = sNNclust_minPts, borderPoints = borderPoints)

    ident_clusters_DT <- data.table::data.table('cell_ID' = cell_id_numeric[1:nrow(kNN_object$dist)], 'name' = sNN_clusters$cluster)
    data.table::setnames(ident_clusters_DT, 'name', name)

    # exit seed
    if(set_seed == TRUE) {
      set.seed(Sys.time())
    }

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
    # parameters to include
    if(cluster_method == 'louvain') {
      parameters_list[[update_name]] = c('nn network' = nn_network_to_use,
                                         'network name' = network_name,
                                         'name for clusters' = name,
                                         'cluster algorithm' = cluster_method,
                                         'pyth leiden resolution' = NA,
                                         'pyth leiden weight' = NA,
                                         'pyth leiden partition' = NA,
                                         'pyth leiden iterations' = NA,
                                         'pyth louvain resolution' = NA,
                                         'pyth louvain weight' = NA,
                                         'louvain gamma' = louvain_gamma,
                                         'louvain omega' = louvain_omega,
                                         'random walk steps' = NA,
                                         'random walk clusters' = NA,
                                         'random walk weights' = NA,
                                         'k for sNNclust' = NA,
                                         'eps for sNNclust' = NA,
                                         'minPts for sNNcluster' = NA,
                                         'assign borderPoints' = NA
      )

    } else if(cluster_method == 'python_leiden') {

      parameters_list[[update_name]] = c('nn network' = nn_network_to_use,
                                         'network name' = network_name,
                                         'name for clusters' = name,
                                         'cluster algorithm' = cluster_method,
                                         'pyth leiden resolution' = pyth_leid_resolution,
                                         'pyth leiden weight' = pyth_leid_weight_col,
                                         'pyth leiden partition' = pyth_leid_part_type,
                                         'pyth leiden iterations' = pyth_leid_iterations,
                                         'pyth louvain resolution' = NA,
                                         'pyth louvain weight' = NA,
                                         'louvain gamma' = NA,
                                         'louvain omega' = NA,
                                         'random walk steps' = NA,
                                         'random walk clusters' = NA,
                                         'random walk weights' = NA,
                                         'k for sNNclust' = NA,
                                         'eps for sNNclust' = NA,
                                         'minPts for sNNcluster' = NA,
                                         'assign borderPoints' = NA
      )


    } else if(cluster_method == 'python_louvain') {

      parameters_list[[update_name]] = c('nn network' = nn_network_to_use,
                                         'network name' = network_name,
                                         'name for clusters' = name,
                                         'cluster algorithm' = cluster_method,
                                         'pyth leiden resolution' = NA,
                                         'pyth leiden weight' = NA,
                                         'pyth leiden partition' = NA,
                                         'pyth leiden iterations' = NA,
                                         'pyth louvain resolution' = pyth_louv_resolution,
                                         'pyth louvain weight' = pyth_louv_weight_col,
                                         'louvain gamma' = NA,
                                         'louvain omega' = NA,
                                         'random walk steps' = NA,
                                         'random walk clusters' = NA,
                                         'random walk weights' = NA,
                                         'k for sNNclust' = NA,
                                         'eps for sNNclust' = NA,
                                         'minPts for sNNcluster' = NA,
                                         'assign borderPoints' = NA
      )

    }
    else if(cluster_method == 'random_walk') {
      parameters_list[[update_name]] = c('nn network' = nn_network_to_use,
                                         'network name' = network_name,
                                         'name for clusters' = name,
                                         'cluster algorithm' = cluster_method,
                                         'pyth leiden resolution' = NA,
                                         'pyth leiden weight' = NA,
                                         'pyth leiden partition' = NA,
                                         'pyth leiden iterations' = NA,
                                         'pyth louvain resolution' = NA,
                                         'pyth louvain weight' = NA,
                                         'louvain gamma' = NA,
                                         'louvain omega' = NA,
                                         'random walk steps' = walk_steps,
                                         'random walk clusters' = walk_clusters,
                                         'random walk weights' = walk_weights,
                                         'k for sNNclust' = NA,
                                         'eps for sNNclust' = NA,
                                         'minPts for sNNcluster' = NA,
                                         'assign borderPoints' = NA
      )

    } else if(cluster_method == 'sNNclust') {
      parameters_list[[update_name]] = c('nn network' = nn_network_to_use,
                                         'network name' = network_name,
                                         'name for clusters' = name,
                                         'cluster algorithm' = cluster_method,
                                         'pyth leiden resolution' = NA,
                                         'pyth leiden weight' = NA,
                                         'pyth leiden partition' = NA,
                                         'pyth leiden iterations' = NA,
                                         'pyth louvain resolution' = NA,
                                         'pyth louvain weight' = NA,
                                         'louvain gamma' = NA,
                                         'louvain omega' = NA,
                                         'random walk steps' = NA,
                                         'random walk clusters' = NA,
                                         'random walk weights' = NA,
                                         'k for sNNclust' = sNNclust_k,
                                         'eps for sNNclust' = sNNclust_eps,
                                         'minPts for sNNcluster' = sNNclust_minPts,
                                         'assign borderPoints' = borderPoints
      )

    }

    gobject@parameters = parameters_list

    return(gobject)
  } else {
    return(ident_clusters_DT)
  }
}





#' @title doLeidenCluster
#' @name doLeidenCluster
#' @description cluster cells using a NN-network and the Leiden community detection algorithm
#' @param gobject giotto object
#' @param nn_network_to_use type of NN network to use (kNN vs sNN)
#' @param network_name name of NN network to use
#' @param python_path specify specific path to python if required
#' @param resolution resolution
#' @param weight_col weight column
#' @param partition_type partition type to use
#' @param init_membership initial membership of cells
#' @param n_iterations number of interations
#' @param name name for cluster
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @param set_seed set seed
#' @param seed_number number for seed
#' @param ... additional parameters
#' @return giotto object appended with new cluster
#' @details Description of Leiden clustering method.
#' @export
#' @examples
#'     doLeidenCluster(gobject)
doLeidenCluster <- function(gobject,
                            nn_network_to_use = 'sNN',
                            network_name = 'sNN.pca',
                            python_path = NULL,
                            resolution = 1,
                            weight_col = 'weight',
                            partition_type = c('RBConfigurationVertexPartition', 'ModularityVertexPartition'),
                            init_membership = NULL,
                            n_iterations = 1000,
                            name = 'pleiden_clus',
                            return_gobject = TRUE,
                            set_seed = T,
                            seed_number = 1234,
                            ...) {


  result = clusterCells(cluster_method = 'python_leiden',
                        gobject = gobject,
                        nn_network_to_use = nn_network_to_use,
                        network_name = network_name,
                        python_path = python_path,
                        pyth_leid_resolution = resolution,
                        pyth_leid_weight_col = weight_col,
                        pyth_leid_part_type = partition_type,
                        pyth_leid_init_memb = init_membership,
                        pyth_leid_iterations = n_iterations,
                        name = name,
                        return_gobject = return_gobject,
                        set_seed = set_seed,
                        seed_number = seed_number,
                        ...)

  return(result)

}



#' @title doLouvainPCluster
#' @name doLouvainPCluster
#' @description cluster cells using a NN-network and the python Louvain community detection algorithm
#' @param gobject giotto object
#' @param nn_network_to_use type of NN network to use (kNN vs sNN)
#' @param network_name name of NN network to use
#' @param python_path specify specific path to python if required
#' @param resolution resolution
#' @param weight_col weight column
#' @param louv_random random boolean
#' @param name name for cluster
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @param set_seed set seed
#' @param seed_number number for seed
#' @param ... additional parameters
#' @return giotto object appended with new cluster
#' @details Description of python Louvain clustering method.
#' @export
#' @examples
#'     doLouvainPCluster(gobject)
doLouvainPCluster <- function(gobject,
                              nn_network_to_use = 'sNN',
                              network_name = 'sNN.pca',
                              python_path = NULL,
                              resolution = 1,
                              weight_col = NULL,
                              louv_random = F,
                              name = 'plouvain_clus',
                              return_gobject = TRUE,
                              set_seed = F,
                              seed_number = 1234,
                              ...) {


  result = clusterCells(cluster_method = 'python_louvain',
                        gobject = gobject,
                        nn_network_to_use = nn_network_to_use,
                        network_name = network_name,
                        pyth_louv_resolution = resolution,
                        pyth_louv_weight_col = weight_col,
                        python_louv_random = louv_random,
                        name = name,
                        return_gobject = return_gobject,
                        set_seed = set_seed,
                        seed_number = seed_number,
                        ...)

  return(result)

}


#' @title doLouvainRCluster
#' @name doLouvainRCluster
#' @description cluster cells using a NN-network and the R Louvain community detection algorithm
#' @param gobject giotto object
#' @param nn_network_to_use type of NN network to use (kNN vs sNN)
#' @param network_name name of NN network to use
#' @param gamma gamma parameter
#' @param omega omega parameter
#' @param name name for cluster
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @param set_seed set seed
#' @param seed_number number for seed
#' @param ... additional parameters
#' @return giotto object appended with new cluster
#' @details Description of R Louvain clustering method.
#' @export
#' @examples
#'     doLouvainRCluster(gobject)
doLouvainRCluster <- function(gobject,
                              nn_network_to_use = 'sNN',
                              network_name = 'sNN.pca',
                              gamma = 1,
                              omega = 1,
                              name = 'rlouvain_clus',
                              return_gobject = TRUE,
                              set_seed = F,
                              seed_number = 1234,
                              ...) {

  result = clusterCells(cluster_method = 'louvain',
                        gobject = gobject,
                        nn_network_to_use = nn_network_to_use,
                        network_name = network_name,
                        louvain_gamma = gamma,
                        louvain_omega = omega,
                        name = name,
                        return_gobject = return_gobject,
                        set_seed = set_seed,
                        seed_number = seed_number,
                        ...)

  return(result)
}


#' @title doRandomWalkCluster
#' @name doRandomWalkCluster
#' @description cluster cells using a NN-network and the R Louvain community detection algorithm
#' @param gobject giotto object
#' @param nn_network_to_use type of NN network to use (kNN vs sNN)
#' @param network_name name of NN network to use
#' @param walk_steps number of random steps
#' @param walk_clusters number of clusters
#' @param walk_weights weight column
#' @param name name for cluster
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @param set_seed set seed
#' @param seed_number number for seed
#' @param ... additional parameters
#' @return giotto object appended with new cluster
#' @details Description of igraph Random Walk clustering method.
#' @export
#' @examples
#'     doRandomWalkCluster(gobject)
doRandomWalkCluster <- function(gobject,
                                nn_network_to_use = 'sNN',
                                network_name = 'sNN.pca',
                                walk_steps = 4,
                                walk_clusters = 10,
                                walk_weights = NA,
                                name = 'random_walk_clus',
                                return_gobject = TRUE,
                                set_seed = F,
                                seed_number = 1234,
                                ...) {

  result = clusterCells(cluster_method = 'randomwalk',
                        gobject = gobject,
                        nn_network_to_use = nn_network_to_use,
                        network_name = network_name,
                        walk_steps = walk_steps,
                        walk_clusters = walk_clusters,
                        walk_weights = walk_weights,
                        name = name,
                        return_gobject = return_gobject,
                        set_seed = set_seed,
                        seed_number = seed_number,
                        ...)

  return(result)

}




#' @title doSNNCluster
#' @name doSNNCluster
#' @description cluster cells using a NN-network and the SNNclust community detection algorithm
#' @param gobject giotto object
#' @param nn_network_to_use type of NN network to use (kNN vs sNN)
#' @param network_name name of NN network to use
#' @param k number of clusters
#' @param eps eps parameter dbscan
#' @param minPts minimum Points
#' @param borderPoints include points at border
#' @param name name for cluster
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @param set_seed set seed
#' @param seed_number number for seed
#' @param ... additional parameters
#' @return giotto object appended with new cluster
#' @details Description of SNNclust clustering method.
#' @export
#' @examples
#'     doSNNCluster(gobject)
doSNNCluster <- function(gobject,
                         nn_network_to_use = 'kNN',
                         network_name = 'kNN.pca',
                         k = 20,
                         eps = 4,
                         minPts = 16,
                         borderPoints = TRUE,
                         name = 'sNN_clus',
                         return_gobject = TRUE,
                         set_seed = F,
                         seed_number = 1234,
                         ...) {


  result = clusterCells(cluster_method = 'sNNclust',
                        gobject = gobject,
                        nn_network_to_use = nn_network_to_use,
                        network_name = network_name,
                        sNNclust_k = k,
                        sNNclust_eps = eps,
                        sNNclust_minPts = minPts,
                        borderPoints = borderPoints,
                        name = name,
                        return_gobject = return_gobject,
                        set_seed = set_seed,
                        seed_number = seed_number,
                        ...)

  return(result)

}







#' @title iterCluster
#' @name iterCluster
#' @description cluster cells iteratively
#' @param gobject giotto object
#' @param nr_rounds number of iterative rounds
#' @param hvg_param parameters for calculateHVG
#' @param hvg_min_perc_cells threshold for detection in min percentage of cells
#' @param hvg_mean_expr_det threshold for mean expression level in cells with detection
#' @param pca_param parameters for runPCA
#' @param nn_param parameters for parameters for runPCA
#' @param resolution resolution for Leiden clustering
#' @param n_iterations number of iterations for Leiden clustering
#' @param python_path python path to use for Leiden clustering
#' @param nn_network_to_use NN network to use
#' @param network_name NN network name
#' @param name name of clustering
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @param ... additional parameters
#' @return giotto object appended with new cluster
#' @details Description of iterative clustering.
#' @export
#' @examples
#'     iterCluster(gobject)
iterCluster <- function(gobject,
                        nr_rounds = 5,
                        hvg_param = list(reverse_log_scale = T, difference_in_variance = 1, expression_values = 'normalized'),
                        hvg_min_perc_cells = 5,
                        hvg_mean_expr_det = 1,
                        pca_param = list(expression_values = 'custom', scale.unit = T),
                        nn_param = list(dimensions_to_use = 1:20),
                        k_neighbors = 20,
                        resolution = 1,
                        n_iterations = 1000,
                        python_path = "/Users/rubendries/Bin/anaconda3/envs/py36/bin/python",
                        nn_network_to_use = 'sNN',
                        network_name = 'sNN.pca',
                        name = 'iter_clus',
                        return_gobject = TRUE,
                        ...) {


  final_cluster_list = list()
  final_groups_list = list()


  # create temporary giotto object
  temp_giotto = gobject

  # iterations
  for(round in 1:nr_rounds) {

    ## calculate stats
    temp_giotto <- addStatistics(gobject = temp_giotto)

    ## calculate variable genes
    temp_giotto = do.call('calculateHVG', c(gobject = temp_giotto, hvg_param))

    ## get hvg
    gene_metadata = fDataDT(temp_giotto)
    featgenes     = gene_metadata[hvg == 'yes' & perc_cells >= hvg_min_perc_cells & mean_expr_det >= hvg_mean_expr_det]$gene_ID


    ## run PCA
    #temp_giotto   = runPCA(gobject = temp_giotto, genes_to_use = featgenes, expression_values = 'custom', scale.unit = T)
    temp_giotto = do.call('runPCA', c(gobject =  temp_giotto, genes_to_use = list(featgenes), pca_param))

    ## nearest neighbor and clustering
    #temp_giotto = createNearestNetwork(gobject = temp_giotto, dimensions_to_use = 1:20)
    temp_giotto = do.call('createNearestNetwork', c(gobject = temp_giotto, k = k_neighbors, nn_param))

    ## Leiden Cluster
    ## TO DO: expand to all clustering options
    temp_giotto = doLeidenCluster(gobject = temp_giotto,
                                  resolution = resolution,
                                  n_iterations = n_iterations,
                                  python_path = python_path,
                                  name = 'tempclus',
                                  ...)

    ## get network and add cluster
    my_nn_network = temp_giotto@nn_network[[nn_network_to_use]][[network_name]][['igraph']]
    cell_metadata = pDataDT(temp_giotto)
    my_nn_network = igraph::set_vertex_attr(my_nn_network, name = 'tempcluster', value = cell_metadata[['tempclus']])

    ## convert network into data.table
    edgeDT = data.table::as.data.table(igraph::as_data_frame(x = my_nn_network, what = 'edges'))
    vertexDT = data.table::as.data.table(igraph::as_data_frame(x = my_nn_network, what = 'vertices'))

    ## identify if edge goes to inside or outside of cluster
    edgeDT <- merge(x = edgeDT, by.x = 'from', y = vertexDT, by.y = 'name')
    data.table::setnames(edgeDT, 'tempcluster', 'from_clus')
    edgeDT <- merge(x = edgeDT, by.x = 'to', y = vertexDT, by.y = 'name')
    data.table::setnames(edgeDT, 'tempcluster', 'to_clus')
    edgeDT[, type_edge := ifelse(from_clus == to_clus, 'same', 'other')]


    ## calculate ratio of inside edges / outside edges for each cluster
    vcountlist = c()
    ratiolist = c()
    indexlist = c()


    ## number of clusters ##
    nr_clusters = length(unique( cell_metadata[['tempclus']] ))

    for(i in sort(unique( cell_metadata[['tempclus']] ))) {

      #cat('\n \n for :', i, '\n')

      # subset edgeDT for each cluster
      sub_edgeDT = edgeDT[from_clus == i | to_clus == i]
      # number of unique vertices
      total_v = length(unique(c(sub_edgeDT$to, sub_edgeDT$from)))
      # calculate ratio
      mytable = table(sub_edgeDT$type_edge)
      within_edges = mytable['same']
      within_edges = ifelse(is.na(within_edges), 0, within_edges)
      outside_edges = mytable['other']
      outside_edges = ifelse(is.na(outside_edges), 0, outside_edges)

      ratio = (within_edges+1)/(outside_edges+1)

      indexlist[[i]] = i
      vcountlist[[i]] = total_v
      ratiolist[[i]] = ratio

    }


    mytempDT = data.table::data.table(index = indexlist, ratio = ratiolist, vc = vcountlist)
    print(mytempDT)

    # identify cluster with the maximum ratio (most tx coherent cluster)
    cell_cluster_index = mytempDT[ratio == max(ratio)][['index']]
    cat('\n index to keep: ', cell_cluster_index, '\n')

    # identify ids from selected cluster and all others outside cluster
    store_ids = cell_metadata[tempclus == cell_cluster_index][['cell_ID']]
    remain_ids = cell_metadata[tempclus != cell_cluster_index][['cell_ID']]

    final_cluster_list[[round]] = store_ids
    final_groups_list[[round]] = rep(round, length(store_ids))

    if(round == nr_rounds | nr_clusters == 1) {

      final_cluster_list[[round]] = cell_metadata[['cell_ID']]
      final_groups_list[[round]] = rep(round, length(cell_metadata[['cell_ID']]))

      if(nr_clusters == 1) {
        break
        print(i)
        print('end \n')
      }

    } else {
      temp_giotto = subsetGiotto(temp_giotto, cell_ids = remain_ids)
    }

  }


  # combine information
  cell_ids    = unlist(final_cluster_list)
  cell_groups = unlist(final_groups_list)
  annot_DT = data.table(cell_ids = cell_ids, cell_groups = cell_groups)
  setnames(annot_DT, 'cell_groups', name)



  ## add clusters to metadata ##
  if(return_gobject == TRUE) {

    cluster_names = names(gobject@cell_metadata)
    if(name %in% cluster_names) {
      cat('\n ', name, ' has already been used, will be overwritten \n')
      cell_metadata = gobject@cell_metadata
      cell_metadata[, eval(name) := NULL]
      gobject@cell_metadata = cell_metadata
    }

    gobject = addCellMetadata(gobject = gobject, new_metadata = annot_DT[, c('cell_ids', name), with = F],
                              by_column = T, column_cell_ID = 'cell_ids')


    ## update parameters used ##
    parameters_list = gobject@parameters
    number_of_rounds = length(parameters_list)
    update_name = paste0(number_of_rounds,'_iter_cluster')

    # parameters to include
    parameters_list[[update_name]] = c('iterclus name' = name,
                                       'number of iterative rounds ' = nr_rounds,
                                       'number of iterations leiden ' = n_iterations,
                                       'resolution ' = resolution)

    gobject@parameters = parameters_list

    return(gobject)
  } else {
    return(annot_DT)
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
#' @param iter.max kmeans maximum iterations
#' @param nstart kmeans nstart
#' @param algorithm kmeans algorithm
#' @param name name for kmeans clustering
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @param set_seed set seed
#' @param seed_number number for seed
#' @param ... additional parameters
#' @return giotto object appended with new cluster
#' @details Description on how to use Kmeans clustering method.
#' @export
#' @examples
#'     doKmeans(gobject)
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
                     iter.max = 100,
                     nstart = 1000,
                     algorithm = "Hartigan-Wong",
                     name = 'kmeans',
                     return_gobject = TRUE,
                     set_seed = T,
                     seed_number = 1234) {



  dim_reduction_to_use = match.arg(dim_reduction_to_use, choices = c('cells', 'pca', 'umap', 'tsne'))
  distance_method = match.arg(distance_method, choices = c("pearson", "spearman",  "original",
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
    matrix_to_use = t(expr_values)

  }

  ## distance
  if(distance_method == 'original') {
    celldist = matrix_to_use
  } else if(distance_method %in% c('spearman', 'pearson')) {
    celldist = stats::as.dist(1-cor(x = t(matrix_to_use), method = distance_method))
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
                            iter.max = iter.max, nstart = nstart,
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
                                       'iter.max' = iter.max,
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
#' @param ... additional parameters
#' @return giotto object appended with new cluster
#' @details Description on how to use Kmeans clustering method.
#' @export
#' @examples
#'     doHclust(gobject)
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
    matrix_to_use = t(expr_values)

  }

  ## distance
  if(distance_method == 'original') {
    celldist = matrix_to_use
  } else if(distance_method %in% c('spearman', 'pearson')) {
    celldist = stats::as.dist(1-cor(x = t(matrix_to_use), method = distance_method))
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





#' @title doLeidenSubCluster
#' @name doLeidenSubCluster
#' @description cluster cells iteratively
#' @param gobject giotto object
#' @param cluster_column cluster column to subcluster
#' @param selected_clusters only do subclustering on these clusters
#' @param hvg_param parameters for calculateHVG
#' @param hvg_min_perc_cells threshold for detection in min percentage of cells
#' @param hvg_mean_expr_det threshold for mean expression level in cells with detection
#' @param pca_param parameters for runPCA
#' @param nn_param parameters for parameters for runPCA
#' @param k_neighbors number of k for sNN
#' @param resolution resolution for Leiden clustering
#' @param n_iterations number of iterations for Leiden clustering
#' @param python_path python path to use for Leiden clustering
#' @param nn_network_to_use NN network to use
#' @param network_name NN network name
#' @param name name of clustering
#' @param return_gobject boolean: return giotto object (default = TRUE)
#' @param verbose be verbose
#' @param ... additional parameters
#' @return giotto object appended with new cluster
#' @details Description of iterative clustering.
#' @export
#' @examples
#'     doLeidenSubCluster(gobject)
doLeidenSubCluster = function(gobject,
                              cluster_column = NULL,
                              selected_clusters = NULL,
                              hvg_param = list(reverse_log_scale = T, difference_in_variance = 1, expression_values = 'normalized'),
                              hvg_min_perc_cells = 5,
                              hvg_mean_expr_det = 1,
                              pca_param = list(expression_values = 'custom', scale.unit = T),
                              nn_param = list(dimensions_to_use = 1:20),
                              k_neighbors = 10,
                              resolution = 0.5,
                              n_iterations = 500,
                              python_path = "/Users/rubendries/Bin/anaconda3/envs/py36/bin/python",
                              nn_network_to_use = 'sNN',
                              network_name = 'sNN.pca',
                              name = 'sub_pleiden_clus',
                              return_gobject = TRUE,
                              verbose = T,
                              ...) {




  iter_list = list()

  cell_metadata = pDataDT(gobject)

  if(is.null(cluster_column)) {
    stop('\n You need to provide a cluster column to subcluster on \n')
  }
  unique_clusters = sort(unique(cell_metadata[[cluster_column]]))


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
                                     return_gobject = F,
                                     ...)

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




