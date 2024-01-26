#' Multi omics integration with WNN
#'
#' @param gobject A Giotto object with individual PCA modalities pre-calculated
#' @param k k number, default = 20
#' @param spat_unit spatial unit
#' @param modality_1 modality 1 name. Default = "rna"
#' @param modality_2 modality 2 name. Default = "protein"
#' @param pca_name_modality_1 Default = 'rna.pca'
#' @param pca_name_modality_2 Default = 'protein.pca'
#' @param integrated_feat_type integrated feature type (e.g. 'rna_protein')
#' @param matrix_result_name Default = 'theta_weighted_matrix'
#' @param w_name_modality_1 name for modality 1 weights
#' @param w_name_modality_2 name for modality 2 weights
#' @param verbose be verbose
#'
#' @return A Giotto object with integrated UMAP (integrated.umap) within the dimension_reduction slot and Leiden clusters (integrated_leiden_clus) in the cellular metadata.
#' @export
runWNN <- function(gobject,
                   spat_unit = 'cell',
                   modality_1 = 'rna',
                   modality_2 = 'protein',
                   pca_name_modality_1 = 'rna.pca',
                   pca_name_modality_2 = 'protein.pca',
                   k = 20,
                   integrated_feat_type = NULL,
                   matrix_result_name = NULL,
                   w_name_modality_1 = NULL,
                   w_name_modality_2 = NULL,
                   verbose = FALSE) {

  # validate Giotto object
  if (!inherits(gobject, "giotto")) {
    stop("gobject needs to be a giotto object")
  }

  # validate modalities
  if(!modality_1 %in% names(gobject@dimension_reduction$cells[[spat_unit]]) || !modality_2 %in% names(gobject@dimension_reduction$cells[[spat_unit]]) ) {
    stop(paste(modality_1, "and", modality_2, " pca must exist"))
  }

  # extract PCA

  ## modality 1
  kNN_1 <- get_NearestNetwork(gobject,
                              spat_unit = spat_unit,
                              feat_type = modality_1,
                              nn_network_to_use = "kNN")
  kNN_1 <- slot(kNN_1, "igraph")

  pca_1 <- get_dimReduction(gobject,
                            spat_unit = spat_unit,
                            feat_type = modality_1,
                            reduction = "cells",
                            reduction_method = "pca",
                            name = pca_name_modality_1)
  pca_1 <- slot(pca_1, "coordinates")

  ## modality 2
  kNN_2 <- get_NearestNetwork(gobject,
                              spat_unit = spat_unit,
                              feat_type = modality_2,
                              nn_network_to_use = "kNN")
  kNN_2 <- slot(kNN_2, "igraph")

  pca_2 <- get_dimReduction(gobject,
                            spat_unit = "cell",
                            feat_type = modality_2,
                            reduction = "cells",
                            reduction_method = "pca",
                            name = pca_name_modality_2)
  pca_2 <- slot(pca_2, "coordinates")

  ## get cell names
  cell_names <- unique(igraph::get.edgelist(kNN_1)[,1])

  ######################## distances calculation ############################

  print(paste("Calculating",modality_1,"-",modality_1,"distance"))

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

  print(paste("Calculating",modality_2,"-",modality_2,"distance"))

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

  if(verbose) print(paste("Calculating low dimensional cell-cell distances for", modality_1))

  all_cell_distances_1_1 = dist(pca_1)
  all_cell_distances_1_1 = as.matrix(all_cell_distances_1_1)

  ## modality2 modality2

  if(verbose) print(paste("Calculating low dimensional cell-cell distances for", modality_2))


  all_cell_distances_2_2 = dist(pca_2)
  all_cell_distances_2_2 = as.matrix(all_cell_distances_2_2)


  ######################## within-modality prediction ############################

  if(verbose) print("Calculating within-modality prediction")

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

  if(verbose) print("Calculating cross-modality prediction")

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

  if(verbose) print("Calculating Jaccard similarities")

  ## modality1 modality1
  sNN_1 <- createNearestNetwork(gobject,
                                spat_unit = "cell",
                                feat_type = modality_1,
                                type = "sNN",
                                dim_reduction_to_use = "pca",
                                dim_reduction_name = pca_name_modality_1,
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
                                dim_reduction_name = pca_name_modality_2,
                                dimensions_to_use = 1:100,
                                return_gobject = FALSE,
                                minimum_shared = 1,
                                k = 20)

  sNN_2 <- igraph::as_data_frame(sNN_2)

  if(verbose) print("Calculating kernel bandwidths")

  # cell-specific kernel bandwidth.

  ## modality1
  modality1_sigma_i <- numeric()

  for(cell_a in cell_names) {
    ### 20 small jaccard values
    jaccard_values <- sNN_1[sNN_1$from == cell_a,]

    if (nrow(jaccard_values == 20)) {

      further_cell_cell_distances <- all_cell_distances_1_1[cell_a,jaccard_values$to]

    } else {
      further_cell_cell_distances <- tail(sort(all_cell_distances_1_1[cell_a,]), 20)
    }

    modality1_sigma_i[cell_a] <- mean(further_cell_cell_distances) #  cell-specific kernel bandwidth.
  }

  ## modality2

  modality2_sigma_i <- numeric()

  for(cell_a in cell_names) {
    ### 20 small jaccard values
    jaccard_values <- sNN_2[sNN_2$from == cell_a,]

    if (nrow(jaccard_values == 20)) {
      further_cell_cell_distances <- all_cell_distances_2_2[cell_a,jaccard_values$to]
    } else {
      further_cell_cell_distances <- tail(sort(all_cell_distances_2_2[cell_a,]), 20)
    }

    modality2_sigma_i[cell_a] <- mean(further_cell_cell_distances) #  cell-specific kernel bandwidth.
  }


  ###################### cell-specific modality weights ##########################

  if(verbose) print("Calculating modality weights")

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

  if(verbose) print("Calculating WNN")

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

  if(verbose) print("Calculating WNN normalization")

  w_modality1 <- rep(0, length(cell_names))
  names(w_modality1) = cell_names

  for (cell_a in cell_names) {
    w_modality1[cell_a] <- exp(ratio_modality1[[cell_a]])/(exp(ratio_modality1[[cell_a]]) + exp(ratio_modality2[[cell_a]]))
  }

  w_modality2 <- rep(0, length(cell_names))
  names(w_modality2) = cell_names

  for (cell_a in cell_names) {
    w_modality2[cell_a] <- exp(ratio_modality2[[cell_a]])/(exp(ratio_modality1[[cell_a]]) + exp(ratio_modality2[[cell_a]]))
  }


  ######################### Calculating a WNN graph ##############################

  if(verbose) print("Calculating WNN graph")

  theta_weighted <- matrix(rep(0,length(cell_names)*length(cell_names)),
                           ncol = length(cell_names),
                           nrow = length(cell_names))

  colnames(theta_weighted) <- cell_names
  rownames(theta_weighted) <- cell_names

  kernelpower <- 1

  ## theta_modality1

  theta_modality1_cella_cellb <- exp(-1*(all_cell_distances_1_1/ modality1_sigma_i) ** kernelpower)

  ## theta_modality2
  theta_modality2_cella_cellb <- exp(-1*(all_cell_distances_2_2/ modality2_sigma_i) ** kernelpower)

  ## theta_weighted
  theta_weighted <-  w_modality1*theta_modality1_cella_cellb + w_modality2*theta_modality2_cella_cellb


  # save theta_weighted

  if(verbose) print("Saving WNN results")


  ## set integrated feat_type and result name if not provided
  if(is.null(integrated_feat_type)) {integrated_feat_type = paste0(modality_1,'_',modality_2)}

  if(is.null(matrix_result_name)) {matrix_result_name = 'theta_weighted_matrix'}

  gobject <- set_multiomics(gobject = gobject,
                            result = theta_weighted,
                            spat_unit = spat_unit,
                            feat_type = integrated_feat_type,
                            integration_method = 'WNN',
                            result_name = matrix_result_name,
                            verbose = TRUE)


  # save modalities weight

  ## modality 1
  if(is.null(w_name_modality_1)) {w_name_modality_1 = paste0('w_',modality_1)}

  gobject <- set_multiomics(gobject = gobject,
                            result = w_modality1,
                            spat_unit = spat_unit,
                            feat_type = integrated_feat_type,
                            integration_method = 'WNN',
                            result_name = w_name_modality_1,
                            verbose = TRUE)

  ## modality 2
  if(is.null(w_name_modality_2)) {w_name_modality_2 = paste0('w_',modality_2)}

  gobject <- set_multiomics(gobject = gobject,
                            result = w_modality2,
                            spat_unit = spat_unit,
                            feat_type = integrated_feat_type,
                            integration_method = 'WNN',
                            result_name = w_name_modality_2,
                            verbose = TRUE)

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
#' @param integrated_feat_type integrated feature type (e.g. 'rna_protein')
#' @param integration_method multiomics integration method used. Default = 'WNN'
#' @param matrix_result_name Default = 'theta_weighted_matrix'
#' @param force force calculation of integrated kNN. Default = FALSE
#' @param ... additional UMAP parameters
#'
#' @return A Giotto object with integrated UMAP
#' @export
runIntegratedUMAP <- function(gobject,
                              spat_unit = "cell",
                              modality1 = "rna",
                              modality2 = "protein",
                              integrated_feat_type = NULL,
                              integration_method = 'WNN',
                              matrix_result_name = 'theta_weighted_matrix',
                              k = 20,
                              spread = 5,
                              min_dist = 0.01,
                              force = FALSE,
                              ...) {

  if(is.null(integrated_feat_type)) {
    integrated_feat_type = paste0(modality1,'_',modality2)}

  theta_weighted <- get_multiomics(gobject,
                                   spat_unit = spat_unit,
                                   feat_type = integrated_feat_type,
                                   integration_method = integration_method,
                                   result_name = matrix_result_name)

  #theta_weighted <- gobject@dimension_reduction$cells$cell$WNN$theta_weighted
  theta_weighted[is.na(theta_weighted)] <- 0

  if(is.null(gobject@nn_network[[spat_unit]][[modality1]]$kNN$integrated_kNN) || force == TRUE) {

    ################# Calculate integrated Nearest Neighbors #######################

    print('Calculating integrated Nearest Neighbors')

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

    print('Creating igraph')

    nn_network_igraph = igraph::graph_from_data_frame(nn_network_dt[,.(from_cell_ID, to_cell_ID, weight, distance)],
                                                      directed = TRUE,
                                                      vertices = all_index)

    ## store igraph
    nnNetObj = create_nn_net_obj(name = 'integrated_kNN',
                                 nn_type = 'kNN',
                                 igraph = nn_network_igraph,
                                 spat_unit = spat_unit,
                                 feat_type = modality1)

    gobject = set_NearestNetwork(gobject = gobject,
                                 nn_network = nnNetObj,
                                 spat_unit = spat_unit,
                                 feat_type = modality1,
                                 nn_network_to_use = 'kNN',
                                 network_name = 'integrated_kNN')

    ## store nn_network id
    gobject <- set_multiomics(gobject = gobject,
                              result = nn_network$id,
                              spat_unit = spat_unit,
                              feat_type = integrated_feat_type,
                              integration_method = 'WNN',
                              result_name = 'integrated_kNN_id',
                              verbose = TRUE)

    ## store nn_network dist
    gobject <- set_multiomics(gobject = gobject,
                              result = nn_network$dist,
                              spat_unit = spat_unit,
                              feat_type = integrated_feat_type,
                              integration_method = 'WNN',
                              result_name = 'integrated_kNN_dist',
                              verbose = TRUE)
  }

  ######################### Calculate integrated UMAP ############################

  print('Calculating integrated UMAP')

  nn_network_id = get_multiomics(gobject,
                                 spat_unit = spat_unit,
                                 feat_type = integrated_feat_type,
                                 integration_method = integration_method,
                                 result_name = 'integrated_kNN_id')

  nn_network_dist = get_multiomics(gobject,
                                   spat_unit = spat_unit,
                                   feat_type = integrated_feat_type,
                                   integration_method = integration_method,
                                   result_name = 'integrated_kNN_dist')


  #### using nn_network pre-calculation
  set.seed(4567)
  integrated_umap <- uwot::umap(X = theta_weighted,
                                n_neighbors = k,
                                nn_method = list(idx = nn_network_id,
                                                 dist = nn_network_dist),
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

