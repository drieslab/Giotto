#' Multi omics integration with WNN
#'
#' @param gobject A Giotto object with individual PCA modalities pre-calculated
#' @param k k number, default = 20
#'
#' @return A Giotto object with integrated UMAP (integrated.umap) within the dimension_reduction slot and Leiden clusters (integrated_leiden_clus) in the cellular metadata.
#' @export
#'
#' @examples
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
    result_b <- lapply(cell_names,
                       calculate_all_cell_distances)
    result_b <- sapply(result_b, paste, collapse = " ")
    names(result_b) <- cell_names
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
    result_b <- lapply(cell_names,
                       calculate_all_cell_distances)
    result_b <- sapply(result_b, paste, collapse = " ")
    names(result_b) <- cell_names
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

  # ###################### calculate jaccard similarities ##########################
  #
  # ## rna rna
  #
  # function_b <- function(cell_b) {
  #   combined_neighbors <- c(kNN_rna[[cell_a]][[cell_a]], kNN_rna[[cell_b]][[cell_b]])
  #   duplicated_neighbors <- combined_neighbors[duplicated(combined_neighbors)]
  #   result <- length(duplicated_neighbors)/(k)
  #   return(result)
  # }
  #
  # calculate_jaccard_rna_rna <- function(cell_a) {
  #   result <- lapply(cell_names,
  #                    function_b)
  #   result <- sapply(result, paste, collapse = " ")
  #   result <- as.numeric(result)
  #   names(result) <- cell_names
  #   return(result)
  # }
  #
  # jaccard_rna_rna <- lapply(cell_names,
  #                           calculate_jaccard_rna_rna)
  #
  # names(jaccard_rna_rna) <- cell_names
  #
  #
  # ## protein protein
  #
  # function_b <- function(cell_b) {
  #   combined_neighbors <- c(kNN_protein[[cell_a]][[cell_a]], kNN_protein[[cell_b]][[cell_b]])
  #   duplicated_neighbors <- combined_neighbors[duplicated(combined_neighbors)]
  #   result <- length(duplicated_neighbors)/(k)
  #   return(result)
  # }
  #
  # calculate_jaccard_protein_protein <- function(cell_a) {
  #   result <- lapply(cell_names,
  #                    function_b)
  #   result <- sapply(result, paste, collapse = " ")
  #   result <- as.numeric(result)
  #   names(result) <- cell_names
  #   return(result)
  # }
  #
  # jaccard_protein_protein <- lapply(cell_names,
  #                                   calculate_jaccard_protein_protein)
  #
  # names(jaccard_protein_protein) <- cell_names

  ###################### cell-specific modality weights ##########################

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

    ## 20 small jaccard values
    # non_zero_jaccard_values <- jaccard_rna_rna[[cell_a]][jaccard_rna_rna[[cell_a]] > 0]
    # jaccard_values <- head(sort(non_zero_jaccard_values), 20)
    #
    # if (length(unique(jaccard_values)) < 20) {
      further_cell_cell_distances <- tail(sort(all_cell_distances_rna_rna[[cell_a]]), 20)
    # } else {
    #   further_cell_cell_distances <- all_cell_distances_rna_rna[[cell_a]][names(jaccard_values)]
    # }

    sigma_i <- mean(further_cell_cell_distances) #  cell-specific kernel bandwidth.

    theta_rna_rna[[cell_a]] <- exp( (-max_value)/(sigma_i - d_rnai_rknn1) )
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

    ## 20 small jaccard values
    # non_zero_jaccard_values <- jaccard_protein_protein[[cell_a]][jaccard_protein_protein[[cell_a]] > 0]
    # jaccard_values <- head(sort(non_zero_jaccard_values), 20)
    #
    # if (length(unique(jaccard_values)) < 20) {
      further_cell_cell_distances <- tail(sort(all_cell_distances_protein_protein[[cell_a]]), 20)
    # } else {
    #   further_cell_cell_distances <- all_cell_distances_protein_protein[[cell_a]][names(jaccard_values)]
    # }

    sigma_i <- mean(further_cell_cell_distances) #  cell-specific kernel bandwidth.

    theta_protein_protein[[cell_a]] <- exp( (-max_value)/(sigma_i - d_proteini_pknn1) )
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

    ## 20 small jaccard values
    # non_zero_jaccard_values <- jaccard_rna_rna[[cell_a]][jaccard_rna_rna[[cell_a]] > 0]
    # jaccard_values <- head(sort(non_zero_jaccard_values), 20)
    #
    # if (length(unique(jaccard_values)) < 20) {
      further_cell_cell_distances <- tail(sort(all_cell_distances_rna_rna[[cell_a]]), 20)
    # } else {
    #   further_cell_cell_distances <- all_cell_distances_rna_rna[[cell_a]][names(jaccard_values)]
    # }

    sigma_i <- mean(further_cell_cell_distances) #  cell-specific kernel bandwidth.

    theta_rna_protein[[cell_a]] <- exp( (-max_value)/(sigma_i - d_rnai_rknn1) )
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

    ## kernel width
    ### evaluate 20 small jaccard values
    # non_zero_jaccard_values <- jaccard_protein_protein[[cell_a]][jaccard_protein_protein[[cell_a]] > 0]
    # jaccard_values <- head(sort(non_zero_jaccard_values), 20)
    #
    # if (length(unique(jaccard_values)) < 20) {
      further_cell_cell_distances <- tail(sort(all_cell_distances_protein_protein[[cell_a]]), 20)
    # } else {
    #   further_cell_cell_distances <- all_cell_distances_protein_protein[[cell_a]][names(jaccard_values)]
    # }

    sigma_i <- mean(further_cell_cell_distances) #  cell-specific kernel bandwidth.

    theta_protein_rna[[cell_a]] <- exp( (-max_value)/(sigma_i - d_proteini_pknn1) )
  }


  ##################### ratio of affinities ######################################
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

    ## kernel width
    ### evaluate 20 small jaccard values
    # non_zero_jaccard_values <- jaccard_rna_rna[[cell_a]][jaccard_rna_rna[[cell_a]] > 0]
    # jaccard_values <- head(sort(non_zero_jaccard_values), 20)
    #
    # if (length(unique(jaccard_values)) < 20) {
      further_cell_cell_distances <- tail(sort(all_cell_distances_rna_rna[[cell_a]]), 20)
    # } else {
    #   further_cell_cell_distances <- all_cell_distances_rna_rna[[cell_a]][names(jaccard_values)]
    # }

    sigma_i_rna <- mean(further_cell_cell_distances) #  cell-specific kernel bandwidth.

    ## kernel width
    ### evaluate 20 small jaccard values
    # non_zero_jaccard_values <- jaccard_protein_protein[[cell_a]][jaccard_protein_protein[[cell_a]] > 0]
    # jaccard_values <- head(sort(non_zero_jaccard_values), 20)
    #
    # if (length(unique(jaccard_values)) < 20) {
      further_cell_cell_distances <- tail(sort(all_cell_distances_protein_protein[[cell_a]]), 20)
    # } else {
    #   further_cell_cell_distances <- all_cell_distances_protein_protein[[cell_a]][names(jaccard_values)]
    # }

    sigma_i_protein <- mean(further_cell_cell_distances) #  cell-specific kernel bandwidth.

    for (cell_b in cell_names) {

      ## theta_rna
      theta_rna_cella_cellb <- exp(-1*(all_cell_distances_rna_rna[[cell_a]][cell_b] / sigma_i_rna ) ** kernelpower)

      ## theta_protein
      theta_protein_cella_cellb <- exp(-1*(all_cell_distances_protein_protein[[cell_a]][cell_b] / sigma_i_protein ) ** kernelpower)

      ## theta_weighted
      theta_weighted[cell_a,cell_b] <-  w_rna[[cell_a]]*theta_rna_cella_cellb + w_protein[[cell_a]]*theta_protein_cella_cellb
    }
  }


  ################# Calculate integrated Nearest Neighbors #######################

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

  gobject@nn_network$cell$kNN$integrated_kNN <- nn_network_igraph

  ######################## calculate leiden clustering ###########################

  gobject <- doLeidenCluster(gobject = gobject,
                             spat_unit = "cell",
                             feat_type = "rna",
                             nn_network_to_use = "kNN",
                             network_name = "integrated_kNN",
                             name = "integrated_leiden_clus")

  ######################### Calculate integrated UMAP ############################

  #### using nn_network pre-calculation
  integrated_umap <- uwot::umap(X = theta_weighted, n_neighbors = 20,
                                nn_method = list(idx = nn_network$id,
                                                 dist = nn_network$dist),
                                spread = 10,
                                min_dist = 0.05,
                                a = 1.5,
                                b = 1.1
  )

  #### using matrix directly
  # integrated_umap <- uwot::umap(X = theta_weighted,
  #                               n_neighbors = 20,
  #                               spread = 10,
  #                               min_dist = 0.05,
  #                               a = 1.5,
  #                               b = 1.1
  # )

  colnames(integrated_umap) <- c("Dim.1", "Dim.2")
  head(integrated_umap)

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
