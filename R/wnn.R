#' Multi omics integration with WNN
#'
#' @param gobject A Giotto object with individual PCA modalities pre-calculated
#' @param k k number, default = 20
#' @param spat_unit spatial unit
#' @param modality_1 modality 1 name. Default = "rna"
#' @param modality_2 modality 2 name. Default = "protein"
#'
#' @return A Giotto object with integrated UMAP (integrated.umap) within the dimension_reduction slot and Leiden clusters (integrated_leiden_clus) in the cellular metadata.
#' @export
runWNN <- function(gobject,
                   spat_unit = "cell",
                   modality_1 = "modality1",
                   modality_2 = "modality2",
                   k = 20) {

  # validate Giotto object
  if (!inherits(gobject, "giotto")) {
    stop("gobject needs to be a giotto object")
  }

  # validate modalities
  if(is.null())

    # extract PCA

    ## modality 1
    kNN_1 <- get_NearestNetwork(gobject,
                                spat_unit = spat_unit,
                                feat_type = modality_1,
                                nn_network_to_use = "kNN")
  #kNN_1 <- gobject@nn_network$cell$modality1$kNN$modality1_kNN.pca
  # pca_1 <- get_dimReduction(gobject,
  #                             spat_unit = "cell",
  #                             feat_type = modality_1,
  #                             reduction = "cells",
  #                             reduction_method = "pca",
  #                             name = "pca")
  pca_1 <- get_dimReduction(gobject,
                            spat_unit = "cell",
                            feat_type = modality_1,
                            reduction = "cells",
                            reduction_method = "pca")

  ## modality 2
  kNN_2 <- get_NearestNetwork(gobject,
                              spat_unit = spat_unit,
                              feat_type = modality_2,
                              nn_network_to_use = "kNN")
  #kNN_2 <- gobject@nn_network$cell$modality2$kNN$modality2_kNN.pca
  # pca_2 <- get_dimReduction(gobject,
  #                                 spat_unit = "cell",
  #                                 feat_type = modality_2,
  #                                 reduction = "cells",
  #                                 reduction_method = "pca",
  #                                 name = "modality2.pca")
  pca_2 <- get_dimReduction(gobject,
                            spat_unit = "cell",
                            feat_type = modality_2,
                            reduction = "cells",
                            reduction_method = "pca")

  ## get cell names
  cell_names <- unique(igraph::get.edgelist(kNN_1)[,1])

  ######################## distances calculation ############################

  ### distances modality1 modality1
  cell_distances_1_1 <- list()

  for (cell_a in cell_names) {

    my_kNN <- kNN_1[[cell_a]][[cell_a]]

    cell_distances_1_1[[cell_a]] <- rep(0, k)
    names(cell_distances_1_1[[cell_a]]) <- names(my_kNN)

    for (cell_i in names(my_kNN)) {
      dimensions_cell_a_i <- pca_1[c(cell_a, cell_i),]
      cell_distances_1_1[[cell_a]][cell_i] <- sqrt(sum((dimensions_cell_a_i[1,] - dimensions_cell_a_i[2,])^2))
    }
  }

  ### distances modality2 modality2
  cell_distances_2_2 <- list()

  for (cell_a in cell_names) {
    my_kNN <- kNN_2[[cell_a]][[cell_a]]

    cell_distances_2_2[[cell_a]] <- rep(0, k)
    names(cell_distances_2_2[[cell_a]]) <- names(my_kNN)

    for (cell_i in names(my_kNN)) {
      dimensions_cell_a_i <- pca_2[c(cell_a, cell_i),]
      cell_distances_2_2[[cell_a]][cell_i] <- sqrt(sum((dimensions_cell_a_i[1,] - dimensions_cell_a_i[2,])^2))
    }
  }

  ########################### all cell-cell distances ############################

  ## modality1 modality1

  print("Calculating low dimensional cell-cell distances for modality 1")

  calculate_all_cell_distances <- function(cell_b) {
    dimensions_cell_a_b <- pca_1[c(cell_a, cell_b),]
    result <- sqrt(sum((dimensions_cell_a_b[1,] - dimensions_cell_a_b[2,])^2))
    return(result)
  }

  calculate_all_cell_distances_1_1 <- function(cell_a) {
    result_b <- sapply(cell_names,
                       calculate_all_cell_distances)
    return(result_b)
  }

  all_cell_distances_1_1 <- lapply(cell_names,
                                   calculate_all_cell_distances_1_1)
  names(all_cell_distances_1_1) <- cell_names

  ## modality2 modality2

  print("Calculating low dimensional cell-cell distances for modality 2")

  calculate_all_cell_distances <- function(cell_b) {
    dimensions_cell_a_b <- pca_2[c(cell_a, cell_b),]
    result <- sqrt(sum((dimensions_cell_a_b[1,] - dimensions_cell_a_b[2,])^2))
    return(result)
  }

  calculate_all_cell_distances_2_2 <- function(cell_a) {
    result_b <- sapply(cell_names,
                       calculate_all_cell_distances)
    return(result_b)
  }

  all_cell_distances_2_2 <- lapply(cell_names,
                                   calculate_all_cell_distances_2_2)
  names(all_cell_distances_2_2) <- cell_names

  ######################## within-modality prediction ############################

  print("Calculating within-modality prediction")

  ### predicted modality1 modality1
  predicted_1_1 <- list()

  for (cell_a in cell_names) {
    dimensions_cell_a <- pca_1[kNN_1[[cell_a]][[cell_a]],]

    predicted_1_1[[cell_a]] <- colSums(dimensions_cell_a)/k
  }

  ### predicted modality2 modality2
  predicted_2_2 <- list()

  for (cell_a in cell_names) {
    dimensions_cell_a <- pca_2[kNN_2[[cell_a]][[cell_a]],]

    predicted_2_2[[cell_a]] <- colSums(dimensions_cell_a)/k
  }

  ######################## cross-modality prediction ############################

  print("Calculating cross-modality prediction")

  ## predicted modality1 modality2
  predicted_1_2 <- list()

  for (cell_a in cell_names) {
    dimensions_cell_a <- pca_1[kNN_2[[cell_a]][[cell_a]],]

    predicted_1_2[[cell_a]] <- colSums(dimensions_cell_a)/k
  }

  ## predicted modality2 modality1
  predicted_2_1 <- list()

  for (cell_a in cell_names) {
    dimensions_cell_a <- pca_2[kNN_1[[cell_a]][[cell_a]],]

    predicted_2_1[[cell_a]] <- colSums(dimensions_cell_a)/k
  }

  ###################### calculate jaccard similarities ##########################

  print("Calculating Jaccard similarities")

  ## modality1 modality1
  sNN_1 <- createNearestNetwork(gobject,
                                spat_unit = "cell",
                                feat_type = modality_1,
                                type = "sNN",
                                dim_reduction_to_use = "pca",
                                dimensions_to_use = 1:100,
                                return_gobject = FALSE,
                                minimum_shared = 1,
                                k = 20)

  sNN_1 <- igraph::as_data_frame(sNN_1)

  ## modality2 modality2

  sNN_2 <- createNearestNetwork(gobject,
                                spat_unit = "cell",
                                feat_type = modality_2,
                                type = "sNN",
                                dim_reduction_to_use = "pca",
                                dimensions_to_use = 1:100,
                                return_gobject = FALSE,
                                minimum_shared = 1,
                                k = 20)

  sNN_2 <- igraph::as_data_frame(sNN_2)

  print("Calculating kernel bandwidths")

  # cell-specific kernel bandwidth.

  ## modality1
  modality1_sigma_i <- numeric()

  for(cell_a in cell_names) {
    ### 20 small jaccard values
    jaccard_values <- sNN_1[sNN_1$from == cell_a,]

    if (nrow(jaccard_values == 20)) {
      further_cell_cell_distances <- all_cell_distances_1_1[[cell_a]][jaccard_values$to]
    } else {
      further_cell_cell_distances <- tail(sort(all_cell_distances_1_1[[cell_a]]), 20)
    }

    modality1_sigma_i[cell_a] <- mean(further_cell_cell_distances) #  cell-specific kernel bandwidth.
  }

  ## modality2

  modality2_sigma_i <- numeric()

  for(cell_a in cell_names) {
    ### 20 small jaccard values
    jaccard_values <- sNN_2[sNN_2$from == cell_a,]

    if (nrow(jaccard_values == 20)) {
      further_cell_cell_distances <- all_cell_distances_2_2[[cell_a]][jaccard_values$to]
    } else {
      further_cell_cell_distances <- tail(sort(all_cell_distances_2_2[[cell_a]]), 20)
    }

    modality2_sigma_i[cell_a] <- mean(further_cell_cell_distances) #  cell-specific kernel bandwidth.
  }


  ###################### cell-specific modality weights ##########################

  print("Calculating modality weights")

  ## modality1 modality1
  theta_1_1 <- list()

  for (cell_a in cell_names) {
    modality1_i <- pca_1[cell_a,] # profile of current cell
    d_modality1_i_modality2_predicted <- sqrt(sum((modality1_i - predicted_1_1[[cell_a]])^2))

    first_knn <- names(sort(cell_distances_1_1[[cell_a]]))[1]
    modality1_knn1 <- pca_1[first_knn,] # profile of the nearest neighbor
    d_modality1_i_modality1_knn1 <- sqrt(sum((modality1_i - modality1_knn1)^2))

    difference_distances <- d_modality1_i_modality2_predicted - d_modality1_i_modality1_knn1
    max_value <- max(c(difference_distances, 0))

    theta_1_1[[cell_a]] <- exp( (-max_value)/(modality1_sigma_i[cell_a] - d_modality1_i_modality1_knn1) )
  }

  ## modality2 modality2
  theta_modality2_modality2 <- list()

  for (cell_a in cell_names) {
    modality2_i <- pca_2[cell_a,] # profile of current cell
    d_modality2_i_modality2_predicted <- sqrt(sum((modality2_i - predicted_2_2[[cell_a]])^2))

    first_knn <- names(sort(cell_distances_2_2[[cell_a]]))[1]
    modality2_knn1 <- pca_2[first_knn,] # profile of the nearest neighbor
    d_modality2_i_modality2_knn1 <- sqrt(sum((modality2_i - modality2_knn1)^2))

    difference_distances <- d_modality2_i_modality2_predicted - d_modality2_i_modality2_knn1
    max_value <- max(c(difference_distances, 0))

    theta_modality2_modality2[[cell_a]] <- exp( (-max_value)/(modality2_sigma_i[cell_a] - d_modality2_i_modality2_knn1) )
  }


  ## modality1 modality2
  theta_modality1_modality2 <- list()

  for (cell_a in cell_names) {
    modality1_i <- pca_1[cell_a,] # profile of current cell
    d_modality1_i_modality2_predicted <- sqrt(sum((modality1_i - predicted_1_2[[cell_a]])^2))

    first_knn <- names(sort(cell_distances_1_1[[cell_a]]))[1]
    modality1_knn1 <- pca_1[first_knn,] # profile of the nearest neighbor
    d_modality1_i_modality1_knn1 <- sqrt(sum((modality1_i - modality1_knn1)^2))

    difference_distances <- d_modality1_i_modality2_predicted - d_modality1_i_modality1_knn1
    max_value <- max(c(difference_distances, 0))

    theta_modality1_modality2[[cell_a]] <- exp( (-max_value)/(modality1_sigma_i[cell_a] - d_modality1_i_modality1_knn1) )
  }


  ## modality2 modality1
  theta_modality2_modality1 <- list()

  for (cell_a in cell_names) {
    modality2_i <- pca_2[cell_a,] # profile of current cell
    d_modality2_i_modality1_predicted <- sqrt(sum((modality2_i - predicted_2_1[[cell_a]])^2))

    first_knn <- names(sort(cell_distances_2_2[[cell_a]]))[1]
    modality2_knn1 <- pca_2[first_knn,] # profile of the nearest neighbor
    d_modality2_i_modality2_knn1 <- sqrt(sum((modality2_i - modality2_knn1)^2))

    difference_distances <- d_modality2_i_modality1_predicted - d_modality2_i_modality2_knn1
    max_value <- max(c(difference_distances, 0))

    theta_modality2_modality1[[cell_a]] <- exp( (-max_value)/(modality2_sigma_i[cell_a] - d_modality2_i_modality2_knn1) )
  }


  ##################### ratio of affinities ######################################

  print("Calculating WNN")

  epsilon = 10^-4

  ## modality1
  ratio_modality1 <- list()

  for (cell_a in cell_names) {
    ratio_modality1[[cell_a]] <- theta_1_1[[cell_a]]/(theta_modality1_modality2[[cell_a]] + epsilon)
  }


  ## modality2
  ratio_modality2 <- list()

  for (cell_a in cell_names) {
    ratio_modality2[[cell_a]] <- theta_modality2_modality2[[cell_a]]/(theta_modality2_modality1[[cell_a]] + epsilon)
  }


  ########################### normalization ######################################

  w_modality1 <- list()

  for (cell_a in cell_names) {
    w_modality1[[cell_a]] <- exp(ratio_modality1[[cell_a]])/(exp(ratio_modality1[[cell_a]]) + exp(ratio_modality2[[cell_a]]))
  }

  w_modality2 <- list()

  for (cell_a in cell_names) {
    w_modality2[[cell_a]] <- exp(ratio_modality2[[cell_a]])/(exp(ratio_modality1[[cell_a]]) + exp(ratio_modality2[[cell_a]]))
  }

  gobject@dimension_reduction$cells$cell[['WNN']][[paste0("weight_",modality_1)]] <- w_modality1
  gobject@dimension_reduction$cells$cell[['WNN']][[paste0("weight_",modality_2)]] <- w_modality2

  ######################### Calculating a WNN graph ##############################

  theta_weighted <- matrix(rep(0,length(cell_names)*length(cell_names)),
                           ncol = length(cell_names),
                           nrow = length(cell_names))

  colnames(theta_weighted) <- cell_names
  rownames(theta_weighted) <- cell_names

  kernelpower <- 1

  for (cell_a in cell_names) {

    for (cell_b in cell_names) {

      ## theta_modality1
      theta_modality1_cella_cellb <- exp(-1*(all_cell_distances_1_1[[cell_a]][cell_b] / modality1_sigma_i[cell_a] ) ** kernelpower)

      ## theta_modality2
      theta_modality2_cella_cellb <- exp(-1*(all_cell_distances_2_2[[cell_a]][cell_b] / modality2_sigma_i[cell_a] ) ** kernelpower)

      ## theta_weighted
      theta_weighted[cell_a,cell_b] <-  w_modality1[[cell_a]]*theta_modality1_cella_cellb + w_modality2[[cell_a]]*theta_modality2_cella_cellb
    }
  }

  gobject@dimension_reduction$cells$cell[['WNN']][['theta_weighted']] <- theta_weighted

  return(gobject)
}


#' Run integrated UMAP
#'
#' @param gobject A giotto object
#' @param spat_unit spatial unit
#' @param modality1 modality 1 name. Default = "rna"
#' @param modality2 modality 2 name. Default = "protein"
#' @param k k number
#' @param spread UMAP param: spread
#' @param min_dist UMAP param: min_dist
#' @param ... additional UMAP parameters
#'
#' @return A Giotto object with integrated UMAP
#' @export
runIntegratedUMAP <- function(gobject,
                              spat_unit = "cell",
                              modality1 = "rna",
                              modality2 = "protein",
                              k = 20,
                              spread = 5,
                              min_dist = 0.01,
                              ...) {

  ################# Calculate integrated Nearest Neighbors #######################

  theta_weighted <- get_dimReduction(gobject = gobject,
                                     spat_unit = spat_unit,
                                     feat_type = "WNN",
                                     name = "theta_weighted")

  #theta_weighted <- gobject@dimension_reduction$cells$cell$WNN$theta_weighted
  theta_weighted[is.na(theta_weighted)] <- 0

  cell_names <- colnames(theta_weighted)

  nn_network = dbscan::kNN(x = theta_weighted, k = k, sort = TRUE)
  from = to = weight = distance = from_cell_ID = to_cell_ID = shared = NULL
  nn_network_dt = data.table::data.table(from = rep(1:nrow(nn_network$id),
                                                    k),
                                         to = as.vector(nn_network$id),
                                         weight = 1/(1 + as.vector(nn_network$dist)),
                                         distance = as.vector(nn_network$dist))
  nn_network_dt[, `:=`(from_cell_ID, cell_names[from])]
  nn_network_dt[, `:=`(to_cell_ID, cell_names[to])]
  all_index = unique(x = c(nn_network_dt$from_cell_ID, nn_network_dt$to_cell_ID))

  ################################ Create igraph #################################

  nn_network_igraph = igraph::graph_from_data_frame(nn_network_dt[,.(from_cell_ID, to_cell_ID, weight, distance)],
                                                    directed = TRUE,
                                                    vertices = all_index)

  #nn_network_igraph
  gobject@nn_network[[spat_unit]][[modality1]]$kNN$integrated_kNN <- nn_network_igraph


  ######################### Calculate integrated UMAP ############################

  print("Calculating integrated UMAP")

  #### using nn_network pre-calculation
  set.seed(4567)
  integrated_umap <- uwot::umap(X = theta_weighted,
                                n_neighbors = k,
                                nn_method = list(idx = nn_network$id,
                                                 dist = nn_network$dist),
                                spread = spread,
                                min_dist = min_dist,
                                ...)

  colnames(integrated_umap) <- c("Dim.1", "Dim.2")

  ## add umap
  gobject@dimension_reduction$cells[[spat_unit]][[modality1]][["umap"]][["integrated.umap"]] <- list(name = "integrated.umap",
                                                                                                     feat_type = modality1,
                                                                                                     spat_unit = spat_unit,
                                                                                                     reduction_method = "umap",
                                                                                                     coordinates = integrated_umap,
                                                                                                     misc = NULL)

  gobject@dimension_reduction$cells[[spat_unit]][[modality2]][["umap"]][["integrated.umap"]] <- list(name = "integrated.umap",
                                                                                                     feat_type = modality2,
                                                                                                     spat_unit = spat_unit,
                                                                                                     reduction_method = "umap",
                                                                                                     coordinates = integrated_umap,
                                                                                                     misc = NULL)

  return(gobject)
}
