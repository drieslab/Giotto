#' Multi omics integration with WNN
#'
#' @param gobject A Giotto object with individual PCA modalities pre-calculated
#' @param k k number, default = 20
#'
#' @return A Giotto object with integrated UMAP (integrated.umap) within the dimension_reduction slot and Leiden clusters (integrated_leiden_clus) in the cellular metadata.
#' @export
#'
runWNN <- function(gobject,
                   k = 20) {

  # extract PCA

  ## rna
  kNN_rna <- gobject@nn_network$cell$rna$kNN$rna_kNN.pca
  pca_rna <- get_dimReduction(gobject,
                              spat_unit = "cell",
                              feat_type = "rna",
                              reduction = "cells",
                              reduction_method = "pca",
                              name = "pca")

  ## protein
  kNN_protein <- gobject@nn_network$cell$protein$kNN$protein_kNN.pca
  pca_protein <- get_dimReduction(gobject,
                                  spat_unit = "cell",
                                  feat_type = "protein",
                                  reduction = "cells",
                                  reduction_method = "pca",
                                  name = "protein.pca")

  ## get cell names
  cell_names <- unique(igraph::get.edgelist(kNN_rna)[,1])

  ######################## distances calculation ############################

  ### distances rna rna
  cell_distances_rna_rna <- list()

  for (cell_a in cell_names) {

    my_kNN <- kNN_rna[[cell_a]][[cell_a]]

    cell_distances_rna_rna[[cell_a]] <- rep(0, k)
    names(cell_distances_rna_rna[[cell_a]]) <- names(my_kNN)

    for (cell_i in names(my_kNN)) {
      dimensions_cell_a_i <- pca_rna[c(cell_a, cell_i),]
      cell_distances_rna_rna[[cell_a]][cell_i] <- sqrt(sum((dimensions_cell_a_i[1,] - dimensions_cell_a_i[2,])^2))
    }
  }

  ### distances protein protein
  cell_distances_protein_protein <- list()

  for (cell_a in cell_names) {
    my_kNN <- kNN_protein[[cell_a]][[cell_a]]

    cell_distances_protein_protein[[cell_a]] <- rep(0, k)
    names(cell_distances_protein_protein[[cell_a]]) <- names(my_kNN)

    for (cell_i in names(my_kNN)) {
      dimensions_cell_a_i <- pca_protein[c(cell_a, cell_i),]
      cell_distances_protein_protein[[cell_a]][cell_i] <- sqrt(sum((dimensions_cell_a_i[1,] - dimensions_cell_a_i[2,])^2))
    }
  }

  ########################### all cell-cell distances ############################

  ## rna rna

  print("Calculating low dimensional cell-cell distances for rna")

  calculate_all_cell_distances <- function(cell_b) {
    dimensions_cell_a_b <- pca_rna[c(cell_a, cell_b),]
    result <- sqrt(sum((dimensions_cell_a_b[1,] - dimensions_cell_a_b[2,])^2))
    return(result)
  }

  calculate_all_cell_distances_rna_rna <- function(cell_a) {
    result_b <- sapply(cell_names,
                       calculate_all_cell_distances)
    return(result_b)
  }

  all_cell_distances_rna_rna <- lapply(cell_names,
                                       calculate_all_cell_distances_rna_rna)
  names(all_cell_distances_rna_rna) <- cell_names

  ## protein protein

  print("Calculating low dimensional cell-cell distances for protein")

  calculate_all_cell_distances <- function(cell_b) {
    dimensions_cell_a_b <- pca_protein[c(cell_a, cell_b),]
    result <- sqrt(sum((dimensions_cell_a_b[1,] - dimensions_cell_a_b[2,])^2))
    return(result)
  }

  calculate_all_cell_distances_protein_protein <- function(cell_a) {
    result_b <- sapply(cell_names,
                       calculate_all_cell_distances)
    return(result_b)
  }

  all_cell_distances_protein_protein <- lapply(cell_names,
                                               calculate_all_cell_distances_protein_protein)
  names(all_cell_distances_protein_protein) <- cell_names

  ######################## within-modality prediction ############################

  print("Calculating within-modality prediction")

  ### predicted rna rna
  predicted_rna_rna <- list()

  for (cell_a in cell_names) {
    dimensions_cell_a <- pca_rna[kNN_rna[[cell_a]][[cell_a]],]

    predicted_rna_rna[[cell_a]] <- colSums(dimensions_cell_a)/k
  }

  ### predicted protein protein
  predicted_protein_protein <- list()

  for (cell_a in cell_names) {
    dimensions_cell_a <- pca_protein[kNN_protein[[cell_a]][[cell_a]],]

    predicted_protein_protein[[cell_a]] <- colSums(dimensions_cell_a)/k
  }

  ######################## cross-modality prediction ############################

  print("Calculating cross-modality prediction")

  ## predicted rna protein
  predicted_rna_protein <- list()

  for (cell_a in cell_names) {
    dimensions_cell_a <- pca_rna[kNN_protein[[cell_a]][[cell_a]],]

    predicted_rna_protein[[cell_a]] <- colSums(dimensions_cell_a)/k
  }

  ## predicted protein rna
  predicted_protein_rna <- list()

  for (cell_a in cell_names) {
    dimensions_cell_a <- pca_protein[kNN_rna[[cell_a]][[cell_a]],]

    predicted_protein_rna[[cell_a]] <- colSums(dimensions_cell_a)/k
  }

  ###################### calculate jaccard similarities ##########################

  print("Calculating Jaccard similarities")

  ## rna rna
  sNN_rna <- createNearestNetwork(my_giotto_object,
                                  spat_unit = "cell",
                                  feat_type = "rna",
                                  type = "sNN",
                                  dim_reduction_to_use = "pca",
                                  dimensions_to_use = 1:100,
                                  return_gobject = FALSE,
                                  minimum_shared = 1,
                                  k = 20)

  sNN_rna <- igraph::as_data_frame(sNN_rna)

  ## protein protein

  sNN_protein <- createNearestNetwork(my_giotto_object,
                                      spat_unit = "cell",
                                      feat_type = "protein",
                                      type = "sNN",
                                      dim_reduction_to_use = "pca",
                                      dimensions_to_use = 1:100,
                                      return_gobject = FALSE,
                                      minimum_shared = 1,
                                      k = 20)

  sNN_protein <- igraph::as_data_frame(sNN_protein)

  print("Calculating kernel bandwidths")

  # cell-specific kernel bandwidth.

  ## rna
  rna_sigma_i <- numeric()

  for(cell_a in cell_names) {
    ### 20 small jaccard values
    jaccard_values <- sNN_rna[sNN_rna$from == cell_a,]

    if (nrow(jaccard_values == 20)) {
      further_cell_cell_distances <- all_cell_distances_rna_rna[[cell_a]][jaccard_values$to]
    } else {
      further_cell_cell_distances <- tail(sort(all_cell_distances_rna_rna[[cell_a]]), 20)
    }

    rna_sigma_i[cell_a] <- mean(further_cell_cell_distances) #  cell-specific kernel bandwidth.
  }

  ## protein

  protein_sigma_i <- numeric()

  for(cell_a in cell_names) {
    ### 20 small jaccard values
    jaccard_values <- sNN_protein[sNN_protein$from == cell_a,]

    if (nrow(jaccard_values == 20)) {
      further_cell_cell_distances <- all_cell_distances_protein_protein[[cell_a]][jaccard_values$to]
    } else {
      further_cell_cell_distances <- tail(sort(all_cell_distances_protein_protein[[cell_a]]), 20)
    }

    protein_sigma_i[cell_a] <- mean(further_cell_cell_distances) #  cell-specific kernel bandwidth.
  }


  ###################### cell-specific modality weights ##########################

  print("Calculating modality weights")

  ## rna rna
  theta_rna_rna <- list()

  for (cell_a in cell_names) {
    rnai <- pca_rna[cell_a,] # profile of current cell
    d_rnai_rnapredicted <- sqrt(sum((rnai - predicted_rna_rna[[cell_a]])^2))

    first_knn <- names(sort(cell_distances_rna_rna[[cell_a]]))[1]
    rknn1 <- pca_rna[first_knn,] # profile of the nearest neighbor
    d_rnai_rknn1 <- sqrt(sum((rnai - rknn1)^2))

    difference_distances <- d_rnai_rnapredicted - d_rnai_rknn1
    max_value <- max(c(difference_distances, 0))

    theta_rna_rna[[cell_a]] <- exp( (-max_value)/(rna_sigma_i[cell_a] - d_rnai_rknn1) )
  }

  ## protein protein
  theta_protein_protein <- list()

  for (cell_a in cell_names) {
    proteini <- pca_protein[cell_a,] # profile of current cell
    d_proteini_proteinpredicted <- sqrt(sum((proteini - predicted_protein_protein[[cell_a]])^2))

    first_knn <- names(sort(cell_distances_protein_protein[[cell_a]]))[1]
    pknn1 <- pca_protein[first_knn,] # profile of the nearest neighbor
    d_proteini_pknn1 <- sqrt(sum((proteini - pknn1)^2))

    difference_distances <- d_proteini_proteinpredicted - d_proteini_pknn1
    max_value <- max(c(difference_distances, 0))

    theta_protein_protein[[cell_a]] <- exp( (-max_value)/(protein_sigma_i[cell_a] - d_proteini_pknn1) )
  }


  ## rna protein
  theta_rna_protein <- list()

  for (cell_a in cell_names) {
    rnai <- pca_rna[cell_a,] # profile of current cell
    d_rnai_proteinpredicted <- sqrt(sum((rnai - predicted_rna_protein[[cell_a]])^2))

    first_knn <- names(sort(cell_distances_rna_rna[[cell_a]]))[1]
    rknn1 <- pca_rna[first_knn,] # profile of the nearest neighbor
    d_rnai_rknn1 <- sqrt(sum((rnai - rknn1)^2))

    difference_distances <- d_rnai_proteinpredicted - d_rnai_rknn1
    max_value <- max(c(difference_distances, 0))

    theta_rna_protein[[cell_a]] <- exp( (-max_value)/(rna_sigma_i[cell_a] - d_rnai_rknn1) )
  }


  ## protein rna
  theta_protein_rna <- list()

  for (cell_a in cell_names) {
    proteini <- pca_protein[cell_a,] # profile of current cell
    d_proteini_rnapredicted <- sqrt(sum((proteini - predicted_protein_rna[[cell_a]])^2))

    first_knn <- names(sort(cell_distances_protein_protein[[cell_a]]))[1]
    pknn1 <- pca_protein[first_knn,] # profile of the nearest neighbor
    d_proteini_pknn1 <- sqrt(sum((proteini - pknn1)^2))

    difference_distances <- d_proteini_rnapredicted - d_proteini_pknn1
    max_value <- max(c(difference_distances, 0))

    theta_protein_rna[[cell_a]] <- exp( (-max_value)/(protein_sigma_i[cell_a] - d_proteini_pknn1) )
  }


  ##################### ratio of affinities ######################################

  print("Calculating WNN")

  epsilon = 10^-4

  ## rna
  ratio_rna <- list()

  for (cell_a in cell_names) {
    ratio_rna[[cell_a]] <- theta_rna_rna[[cell_a]]/(theta_rna_protein[[cell_a]] + epsilon)
  }


  ## protein
  ratio_protein <- list()

  for (cell_a in cell_names) {
    ratio_protein[[cell_a]] <- theta_protein_protein[[cell_a]]/(theta_protein_rna[[cell_a]] + epsilon)
  }


  ########################### normalization ######################################

  w_rna <- list()

  for (cell_a in cell_names) {
    w_rna[[cell_a]] <- exp(ratio_rna[[cell_a]])/(exp(ratio_rna[[cell_a]]) + exp(ratio_protein[[cell_a]]))
  }

  w_protein <- list()

  for (cell_a in cell_names) {
    w_protein[[cell_a]] <- exp(ratio_protein[[cell_a]])/(exp(ratio_rna[[cell_a]]) + exp(ratio_protein[[cell_a]]))
  }


  ######################### Calculating a WNN graph ##############################

  theta_weighted <- matrix(rep(0,length(cell_names)*length(cell_names)),
                           ncol = length(cell_names),
                           nrow = length(cell_names))

  colnames(theta_weighted) <- cell_names
  rownames(theta_weighted) <- cell_names

  kernelpower <- 1

  for (cell_a in cell_names) {

    for (cell_b in cell_names) {

      ## theta_rna
      theta_rna_cella_cellb <- exp(-1*(all_cell_distances_rna_rna[[cell_a]][cell_b] / rna_sigma_i[cell_a] ) ** kernelpower)

      ## theta_protein
      theta_protein_cella_cellb <- exp(-1*(all_cell_distances_protein_protein[[cell_a]][cell_b] / protein_sigma_i[cell_a] ) ** kernelpower)

      ## theta_weighted
      theta_weighted[cell_a,cell_b] <-  w_rna[[cell_a]]*theta_rna_cella_cellb + w_protein[[cell_a]]*theta_protein_cella_cellb
    }
  }

  gobject@dimension_reduction$cells$cell[['WNN']][['theta_weighted']] <- theta_weighted

  return(gobject)
}


#' Run integrated UMAP
#'
#' @param gobject A giotto object
#' @param k k number
#' @param spread UMAP param: spread
#' @param min_dist UMAP param: min_dist
#' @param ... additional UMAP parameters
#'
#' @return A Giotto object with integrated UMAP
#' @export
#'
runIntegratedUMAP <- function(gobject,
                              k = 20,
                              spread = 5,
                              min_dist = 0.01,
                              ...) {
  ################# Calculate integrated Nearest Neighbors #######################

  theta_weighted <- my_giotto_object@dimension_reduction$cells$cell$WNN$theta_weighted
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

  gobject@nn_network$cell$rna$kNN$integrated_kNN <- nn_network_igraph


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
  gobject@dimension_reduction$cells$cell$rna[["umap"]][["integrated.umap"]] <- list(name = "integrated.umap",
                                                                                    feat_type = "rna",
                                                                                    spat_unit = "cell",
                                                                                    reduction_method = "umap",
                                                                                    coordinates = integrated_umap,
                                                                                    misc = NULL)

  gobject@dimension_reduction$cells$cell$protein[["umap"]][["integrated.umap"]] <- list(name = "integrated.umap",
                                                                                        feat_type = "protein",
                                                                                        spat_unit = "cell",
                                                                                        reduction_method = "umap",
                                                                                        coordinates = integrated_umap,
                                                                                        misc = NULL)

  return(gobject)
}
