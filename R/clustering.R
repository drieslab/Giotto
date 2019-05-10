


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
                         pyth_leid_iterations = 200,

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
                         set_seed = F,
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
    add_igraph_layer_ml(mlnetwork = mln_object, g = igraph_object, name = name)

    # start seed
    if(set_seed == TRUE) {
      set.seed(seed = seed_number)
    }

    louvain_clusters = multinet::glouvain_ml(mlnetwork = mln_object, gamma = louvain_gamma, omega = louvain_omega, ...)
    ident_clusters_DT = as.data.table(louvain_clusters)
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
    network_edge_dt = as.data.table(igraph::as_data_frame(x = igraph_object, what = 'edges'))

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


    network_edge_dt = as.data.table(igraph::as_data_frame(x = igraph_object, what = 'edges'))

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

    ident_clusters_DT <- data.table('cell_ID' = V(igraph_object)$name, 'name' = randomwalk_clusters)
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
    igraph_DT = as.data.table(igraph::as_data_frame(igraph_object, what = 'edges'))
    igraph_DT = igraph_DT[order(from)]

    cell_id_numeric = unique(x = c(igraph_DT$from, igraph_DT$to))
    names(cell_id_numeric) <- 1:length(cell_id_numeric)
    igraph_DT[, from_T := as.numeric(names(cell_id_numeric[cell_id_numeric == from])), by = 1:nrow(igraph_DT)]
    igraph_DT[, to_T := as.numeric(names(cell_id_numeric[cell_id_numeric == to])), by = 1:nrow(igraph_DT)]
    temp_igraph_DT = igraph_DT[,.(from_T, to_T, weight, distance)]
    setnames(temp_igraph_DT, old = c('from_T', 'to_T'), new = c('from', 'to'))

    kNN_object = nnDT_to_kNN(nnDT = temp_igraph_DT)
    sNN_clusters = dbscan::sNNclust(x = kNN_object, k = sNNclust_k, eps = sNNclust_eps,
                                    minPts = sNNclust_minPts, borderPoints = borderPoints)

    ident_clusters_DT <- data.table('cell_ID' = cell_id_numeric[1:nrow(kNN_object$dist)], 'name' = sNN_clusters$cluster)
    setnames(ident_clusters_DT, 'name', name)

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
                            n_iterations = 200,
                            name = 'pleiden_clus',
                            return_gobject = TRUE,
                            set_seed = F,
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










