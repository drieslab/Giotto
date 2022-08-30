
## sp classes ####




## anndata object ####

#' @title Convert anndata to Giotto
#' @name anndataToGiotto
#' @description Converts a spatial anndata (e.g. scanpy) .h5ad file into a Giotto object
#' @param anndata_path path to the .h5ad file
#' @param metadata_cols metadata columns to include
#' @param instructions giotto instructions
#' @param \dots additional parameters to \code{\link{createGiottoObject}}
#' @return Giotto object
#' @details Function in beta. Converts a .h5ad file into a Giotto object.
#' @export
anndataToGiotto = function(anndata_path,
                           metadata_cols = c("total_counts", "pct_counts_mt"),
                           instructions = NULL,
                           ...) {



  # test if scanpy is found
  module_test = reticulate::py_module_available('scanpy')
  if(module_test == FALSE) {
    warning("scanpy python module is not installed:
            install in the right environment or python path with:

            'pip install scanpy'

            or try from within R in the Giotto environment with:

            conda_path = reticulate::miniconda_path()
            conda_full_path = paste0(conda_path,'/','bin/conda')
            full_envname = paste0(conda_path,'/envs/giotto_env')
            reticulate::py_install(packages = c('scanpy'),
                                   envname = full_envname,
                                   method = 'conda',
                                   conda = conda_full_path,
                                   pip = TRUE,
                                   python_version = '3.6')")
  }

  # load python modules
  sc <- reticulate::import("scanpy")
  pd <- reticulate::import("pandas")

  if(!file.exists(anndata_path)) stop('path to anndata does not exist \n')
  adata <- sc$read(anndata_path)

  ## get count data
  exprs <- t(adata$X)
  colnames(exprs) <- adata$obs_names$to_list()
  rownames(exprs) <- adata$var_names$to_list()

  ## get spatial data
  spatial <- as.data.frame(adata$obsm["spatial"])
  spatial_names = c('X', 'Y', 'Z')
  colnames(spatial) <- spatial_names[1:ncol(spatial)]
  row.names(spatial) <- colnames(exprs)

  ## get metadata
  obs <- adata$obs
  metadata <- obs[, metadata_cols]

  # create giotto object
  giotto_object <- createGiottoObject(raw_exprs = exprs,
                                      spatial_locs = spatial,
                                      instructions = instructions,
                                      cell_metadata = metadata,
                                      ...)

  return(giotto_object)

}



## Seurat object ####


#' @title Convert Giotto to Seurat
#' @name giottoToSeurat
#' @description Converts Giotto object into a Seurat object
#' @param obj_use Giotto object
#' @param ... additional params to pass to \code{\link{get_spatial_locations}}
#' @return Seurat object
#' @export
giottoToSeurat <- function(obj_use = NULL,
                           ...){

  requireNamespace('Seurat')

  # verify if optional package is installed
  package_check(pkg_name = "Seurat", repository = "CRAN")

  # check whether any raw data exist -- required for Seurat
  raw_exist <- sapply(obj_use@expression_feat,function(x)
    'raw' %in% names(obj_use@expression[[x]]))
  if (length(raw_exist) > 0){
    assays_all <- names(raw_exist[1])
    assays_all <- union(assays_all,obj_use@expression_feat)
  } else {
    stop("Raw count data not found. Required for Seurat object.")
  }

  # create Seurat object when at least one raw data is available
  for (i in 1:length(assays_all)){
    assay_use <- assays_all[i]
    expr_use <- obj_use@expression[[assay_use]]
    slot_use <- names(expr_use)
    if (i == 1){
      data_raw <- expr_use[['raw']]
      sobj <- Seurat::CreateSeuratObject(counts = data_raw,
                                         assay = assay_use)
      if ('normalized' %in% slot_use){
        sobj <- Seurat::SetAssayData(sobj,slot = 'data',
                                     new.data = expr_use[['normalized']],
                                     assay = assay_use)
      }
      if ('scaled' %in% slot_use){
        sobj <- Seurat::SetAssayData(sobj,slot = 'scale.data',
                                     new.data = expr_use[['scaled']],
                                     assay = assay_use)
      }
    } else {
      if ('raw' %in% slot_use){
        data_raw <- expr_use[['raw']]
        flag_raw <- 1
      } else {
        flag_raw <- 0
      }
      if ('normalized' %in% slot_use){
        data_norm <- expr_use[['normalized']]
        flag_norm <- 1
      } else {
        flag_norm <- 0
      }
      if (flag_raw == 1){
        assay_obj <- Seurat::CreateAssayObject(counts = data_raw)
      } else if (flag_raw == 0 & flag_norm == 1){
        assay_obj <- Seurat::CreateAssayObject(data = data_norm)
      } else {
        stop(paste0('Raw and normalized data not found for assay ',assay_use))
      }
      sobj[[assay_use]] <- assay_obj
      if ('scaled' %in% slot_use){
        data_scale <- expr_use[['scaled']]
        sobj <- Seurat::SetAssayData(sobj,slot = "scale.data",
                                     new.data = data_scale,
                                     assay = assay_use)
      }
    }

    # add cell metadata
    meta_cells <- as.data.frame(pDataDT(obj_use,feat_type = assay_use))
    rownames(meta_cells) <- meta_cells$cell_ID
    meta_cells <- meta_cells[,-which(colnames(meta_cells) == 'cell_ID')]
    colnames(meta_cells) <- paste0(assay_use,"_",colnames(meta_cells))
    sobj <- Seurat::AddMetaData(sobj,metadata = meta_cells[Seurat::Cells(sobj),])

    # add feature metadata
    meta_genes <- as.data.frame(fDataDT(obj_use,feat_type = assay_use))
    rownames(meta_genes) <- meta_genes$feat_ID
    sobj[[assay_use]]@meta.features <- cbind(sobj[[assay_use]]@meta.features,meta_genes)
  }

  # spatial coordinates
  loc_use <- as.data.frame(select_spatial_locations(gobject = obj_use,...))
  rownames(loc_use) <- loc_use$cell_ID
  sobj <- Seurat::AddMetaData(sobj,metadata = loc_use)
  # add spatial coordinates as new dim reduct object
  loc_2 <- loc_use[,c('sdimx','sdimy')]
  colnames(loc_2) <- c('spatial_1','spatial_2')
  sobj[['spatial']] <- Seurat::CreateDimReducObject(embeddings = as.matrix(loc_2),
                                            assay = names(sobj@assays)[1],
                                            key = 'spatial_')

  # dim reduction
  # note: Seurat requires assay name specification for each dim reduc
  if (!is.null(obj_use@dimension_reduction)){
    dr_use <- names(obj_use@dimension_reduction[[1]])
    for (i in dr_use){
      dr_name <- i
      dr_obj <- select_dimReduction(gobject = obj_use,return_dimObj = T,reduction_method = dr_name,name = dr_name)
      emb_use <- dr_obj$coordinates[Seurat::Cells(sobj),]
      if (sum(c('loadings','eigenvalues') %in% names(dr_obj$misc)) == 2){
        loadings_use <- dr_obj$misc$loadings
        stdev_use <- dr_obj$misc$eigenvalues
        sobj[[dr_name]] <- Seurat::CreateDimReducObject(embeddings = as.matrix(emb_use),
                                                        loadings = loadings_use,
                                                        key = dr_name,
                                                        stdev = stdev_use,
                                                        assay = names(sobj@assays)[1])
      } else {
        sobj[[dr_name]] <- Seurat::CreateDimReducObject(embeddings = as.matrix(emb_use),
                                                        key = dr_name,
                                                        assay = names(sobj@assays)[1])
      }
    }
  }

  # network objects
  # expression network
  nn_all <- names(obj_use@nn_network)
  if (!is.null(nn_all)){
    for (i in nn_all){
      nn_use <- select_NearestNetwork(gobject = obj_use,
                                      nn_network_to_use = i,
                                      output = 'data.table',
                                      network_name = names(obj_use@nn_network[[i]]))
      idx1 <- match(nn_use$from,Seurat::Cells(sobj))
      idx2 <- match(nn_use$to,Seurat::Cells(sobj))
      edge_weight <- nn_use$weight
      nn_mtx <- Matrix::sparseMatrix(i = idx1,j = idx2,x = edge_weight,dims = c(ncol(sobj),ncol(sobj)))
      rownames(nn_mtx) <- colnames(nn_mtx) <- Seurat::Cells(sobj)
      nn_name <- paste0('expr_',i)
      sobj[[nn_name]] <- Seurat::as.Graph(nn_mtx)
      sobj[[nn_name]]@assay.used <- names(sobj@assays)[1]
    }
  }

  # spatial network
  sn_all <- names(obj_use@spatial_network)
  if (!is.null(sn_all)){
    for (i in sn_all){
      snt_use <- select_spatialNetwork(gobject = obj_use,name = i)
      idx1 <- match(snt_use$from,Seurat::Cells(sobj))
      idx2 <- match(snt_use$to,Seurat::Cells(sobj))
      edge_weight <- snt_use$weight
      nn_mtx <- Matrix::sparseMatrix(i = idx1,j = idx2,x = edge_weight,dims = c(ncol(sobj),ncol(sobj)))
      rownames(nn_mtx) <- colnames(nn_mtx) <- Seurat::Cells(sobj)
      nn_name <- paste0('spatial_',i)
      sobj[[nn_name]] <- Seurat::as.Graph(nn_mtx)
      sobj[[nn_name]]@assay.used <- names(sobj@assays)[1]
    }
  }

  return (sobj)
}


#' @title seuratToGiotto_OLD
#' @name seuratToGiotto_OLD
#' @description Converts Seurat object into a Giotto object. Deprecated, see \code{\link{giottoToSeurat}}
#' @param obj_use Seurat object
#' @param ... additional params to pass
#' @return Giotto object
#' @export
seuratToGiotto_OLD <- function(obj_use = NULL,
                               ...){
  requireNamespace('Seurat')
  requireNamespace('Giotto')

  # get general info in basic seurat structures
  obj_assays <- names(obj_use@assays)
  if ('Spatial' %in% obj_assays){
    obj_assays <- c('Spatial',obj_assays[-which(obj_assays == 'Spatial')])
  }

  obj_dimReduc <- names(obj_use@reductions)
  obj_dimReduc_assay <- sapply(obj_dimReduc,function(x)
    obj_use[[x]]@assay.used)

  obj_graph_expr <- names(obj_use@graphs)
  obj_graph_expr_assay <- sapply(obj_graph_expr,function(x)
    obj_use[[x]]@assay.used)

  obj_meta_cells <- obj_use@meta.data
  obj_meta_genes <- lapply(obj_assays,function(x)
    obj_use[[x]]@meta.features)
  names(obj_meta_genes) <- obj_assays

  obj_img <- obj_use@images
  obj_img_names <- names(obj_img)
  loc_use <- lapply(obj_img_names,function(x){
    temp <- obj_img[[x]]@coordinates
    temp <- as.data.frame(temp[,c('col','row')])
    # temp$region <- x
    return (temp)
  })
  loc_use <- Reduce(rbind,loc_use)

  # add assay data: raw, normalized & scaled
  for (i in 1:length(obj_assays)){
    data_raw <- Seurat::GetAssayData(obj_use,slot = 'counts',assay = obj_assays[i])
    data_norm <- Seurat::GetAssayData(obj_use,slot = 'data',assay = obj_assays[i])
    data_scale <- Seurat::GetAssayData(obj_use,slot = 'scale.data',assay = obj_assays[i])

    if (i == 1 & obj_assays[i] == 'Spatial'){
      feat_use <- 'rna'
      test <- createGiottoObject(expression = obj_use[[obj_assays[i]]]@counts,
                                 spatial_locs = loc_use,
                                 expression_feat = 'rna')
      test <- addCellMetadata(test,feat_type = feat_use,new_metadata = obj_meta_cells)
    } else {
      feat_use <- obj_assays[i]
      test@expression[[feat_use]][['raw']] <- data_raw
      test@feat_ID[[feat_use]] = rownames(data_raw)
      test@feat_metadata[[feat_use]] = data.table::data.table(feat_ID = test@feat_ID[[feat_use]])
    }
    if (nrow(data_norm) > 0){
      test@expression[[feat_use]][['normalized']] <- data_norm
    }
    if (nrow(data_scale) > 0){
      test@expression[[feat_use]][['scaled']] <- data_scale
    }

    # gene metadata
    if (length(obj_meta_genes[[i]]) > 0){
      test <- addFeatMetadata(test,feat_type = feat_use,
                              new_metadata = obj_meta_genes[[i]])
    }
  }

  # add dim reduction
  for (i in obj_dimReduc){
    if (!i %in% c('pca','umap','tsne')){
      next
    } else {
      dimReduc_name <- i
      dimReduc_method <- i
      dimReduc_coords <- obj_use[[i]]@cell.embeddings
      dimReduc_misc <- list(obj_use[[i]]@stdev,
                            obj_use[[i]]@feature.loadings,
                            obj_use[[i]]@feature.loadings.projected)
      names(dimReduc_misc) <- c('eigenvalues','loadings','loadings_projected')
      dimObject <- create_dimObject(name = dimReduc_name,reduction_method = dimReduc_method,
                                             coordinates = dimReduc_coords,misc = dimReduc_misc)
      test@dimension_reduction[['cells']][[dimReduc_method]][[dimReduc_name]] <- dimObject
    }
  }

  # add expr nearest neighbors
  for (i in obj_graph_expr){
    mtx_use <- obj_use[[i]]
    ig_use <- igraph::graph_from_adjacency_matrix(mtx_use,weighted = T)
    g_type <- unlist(strsplit(i,split = "_"))[2]
    g_val <- i
    test@nn_network[[g_type]][[g_val]][['igraph']] <- ig_use
  }
  return (test)

}



#' @title seuratToGiotto
#' @name seuratToGiotto
#' @description Converts Seurat object into a Giotto object
#' @param sobject Seurat object
#' @return Giotto object
#' @export
seuratToGiotto = function(sobject, spatial_assay = 'Spatial',
                          dim_reduction = c('pca','umap'), subcellular_assay = 'Vizgen'){
  requireNamespace('Seurat')
  requireNamespace('Giotto')

  if(is.null(Seurat::GetAssayData(object = sobject, slot = "counts", assay = spatial_assay))) {
    cat('No raw expression values are provided in spatial_assay\n')
    return(sobject)

  } else {

    exp = Seurat::GetAssayData(object = sobject, slot = "counts", assay = spatial_assay)
    if(!is.null(sobject@assays$SCT)){
        normexp = Seurat::GetAssayData(object = sobject, slot = "counts", assay = 'SCT')
      }

      if(!is.null(sobject@assays$Spatial@data)){
        normexp = Seurat::GetAssayData(object = sobject, slot = "data", assay = 'Spatial')
      }

    # Cell Metadata
    cell_metadata = sobject@meta.data

    # Dimension Reduction
    if(sum(sapply(dim_reduction,function(x) length(sobject@reductions[[x]]))) == 0) {
      dim_reduc = NULL
    } else {
      dimReduc_list = lapply(dim_reduction,function(x){
        dim_coord = as.matrix(Seurat::Embeddings(object = sobject,reduction = x))
        dim_load = as.matrix(Seurat::Loadings(object = sobject, reduction = x))
        dim_eig = Seurat::Stdev(sobject, reduction = x)
        if (length(dim_eig) > 0){
          dim_eig = sapply(dim_eig, function(x) x ^ 2)
        }
        colnames(dim_coord) = paste0('Dim.',1:ncol(dim_coord))
        if (length(dim_load) > 0){
          colnames(dim_load) = paste0('Dim.',1:ncol(dim_load))
        }
        dimReduc = list(type = "cells",
                        spat_unit = "cell",
                        name = x,
                        reduction_method = x,
                        coordinates = dim_coord,
                        feat_type = 'rna',
                        misc =list(eigenvalues = dim_eig, loadings = dim_load))
        return (dimReduc)
      })
      # names(dimReduc_list) <- dim_reduction
    }

    # Spatial Locations
    if(length(sobject@assays[[spatial_assay]]) == 0) {
      spat_loc = NULL
    } else {

      spat_coord = Seurat::GetTissueCoordinates(sobject)
      spat_coord = cbind(rownames(spat_coord), data.frame(spat_coord, row.names=NULL))
      colnames(spat_coord) = c("cell_ID", "sdimy", "sdimx")
      spat_loc = spat_coord
    }


    # Subcellular
    name = names(sobject@images)
    if(length(sobject@assays[[subcellular_assay]]) == 1) {

      spat_coord = Seurat::GetTissueCoordinates(sobject)
      colnames(spat_coord) = c("sdimx", "sdimy", "cell_ID")
      exp = exp[  , c(intersect(spat_coord$cell_ID, colnames(exp)))]
      spat_loc = spat_coord
    }
    if (!length(sobject@images) == 0) {
      if ("molecules" %in% methods::slotNames(sobject@images[[name]]) == TRUE) {
        if(!length(sobject@images[[name]][["molecules"]]) == 0) {

          assay = names(sobject@assays)
          featnames = rownames(sobject@assays[[assay]]@meta.features)
          mol_spatlocs = data.table::data.table()

          for (x in featnames) {
            df = (Seurat::FetchData(sobject[[name]][["molecules"]], vars = x))
            mol_spatlocs = rbind(mol_spatlocs, df)
          }
          gpoints = createGiottoPoints(mol_spatlocs, feat_type = "rna")

        }
      }
    }
  }

  gobject = createGiottoObject(exp,
                               spatial_locs = spat_loc,
                               dimension_reduction = dimReduc_list)
   if (exists('normexp') == TRUE) {
    gobject@expression$cell$rna$normalized = normexp
    }
  gobject = addCellMetadata(gobject = gobject, new_metadata = cell_metadata)


  if (exists('gpoints') == TRUE) {
    gobject = addGiottoPoints(gobject = gobject,
                              gpoints = list(gpoints))
  }

  return (gobject)
}




## SpatialExperiment object ####

#' Utility function to convert a Giotto object to a SpatialExperiment object.
#'
#' @param giottoObj Input Giotto object to convert to a SpatialExperiment object.
#'
#' @return A SpatialExperiment object that contains data from the input Giotto object.
#' @export
giottoToSpatialExperiment <- function(giottoObj){

  requireNamespace(c(
    "SummarizedExperiment",
    "SingleCellExperiment",
    "SpatialExperiment",
    "S4Vectors"))

  #expression matrices
  giottoExpr <- Giotto:::list_expression(giottoObj)
  # check if has > 0 rows
  if(!is.null(giottoExpr)){
    message("Copying expression matrix: ", giottoExpr[1]$name)
    exprMat <- get_expression_values(
      gobject = giottoObj,
      spat_unit = giottoExpr[1]$spat_unit,
      feat_type = giottoExpr[1]$feat_type,
      values = giottoExpr[1]$name)
    names(rownames(exprMat)) <- NULL #sp doesnt allow
    names(colnames(exprMat)) <- NULL
    exprMat <- list(exprMat)
    names(exprMat)[1] <- giottoExpr[1]$name
    spe <- SpatialExperiment::SpatialExperiment(assays = exprMat)
    giottoExpr <- giottoExpr[-1, ]
  } else{
    stop("The input Giotto object must contain atleast one expression matrix.")
  }

  if(nrow(giottoExpr) > 0){
    for(i in seq(nrow(giottoExpr))){
      message("Copying expression matrix: ", giottoExpr[i]$name)
      SummarizedExperiment::assay(
        spe,
        paste0(giottoExpr[i]$name, "_",
               giottoExpr[i]$feat_type, "_",
               giottoExpr[i]$spat_unit),
        withDimnames = FALSE) <- get_expression_values(
        gobject = giottoObj,
        spat_unit = giottoExpr[i]$spat_unit,
        feat_type = giottoExpr[i]$feat_type,
        values = giottoExpr[i]$name)
    }
  }

  # metadata to coldata
  # check if exists
  if(nrow(pDataDT(giottoObj)) > 0){
    message("Copying phenotype data")
    pData <- pDataDT(giottoObj)
    SummarizedExperiment::colData(spe) <- S4Vectors::DataFrame(pData, row.names = pData$cell_ID)
  } else{
    message("No phenotype data found in input Giotto object")
  }

  # check if exists
  if(nrow(fDataDT(giottoObj)) > 0){
    message("Copying feature metadata")
    SummarizedExperiment::rowData(spe) <- fDataDT(giottoObj)
  } else{
    message("No feature metadata found in input Giotto object")
  }

  # check if exists
  if(!is.null(get_spatial_locations(giottoObj))){
    message("Copying spatial locations")
    SpatialExperiment::spatialCoords(spe) <- data.matrix(get_spatial_locations(giottoObj)[, 1:2])
  } else{
    message("No spatial locations found in the input Giotto object")
  }


  giottoReductions <- Giotto:::list_dim_reductions(giottoObj)
  if(!is.null(giottoReductions)){
    message("Copying reduced dimensions")
    for(i in seq(nrow(giottoReductions))){
      SingleCellExperiment::reducedDim(spe, giottoReductions[i]$name) <- get_dimReduction(
        gobject = giottoObj,
        reduction = "cells",
        spat_unit = giottoReductions[i]$spat_unit,
        feat_type = giottoReductions[i]$feat_type,
        reduction_method = giottoReductions[i]$dim_type,
        name = giottoReductions[i]$name)
    }
  } else{
    message("No reduced dimensions found in the input Giotto object")
  }

  #also where to store spatialgrid? metadata?

 ## NN graph
  giottoNearestNetworks <- Giotto:::list_nearest_networks(giottoObj)
  if(!is.null(giottoNearestNetworks)){
    message("Copying nearest networks")
    for(i in seq(nrow(giottoNearestNetworks))){
      nn_network <- get_NearestNetwork(
        gobject = giottoObj,
        spat_unit = giottoNearestNetworks[i]$spat_unit,
        output = "data.table",
        nn_network_to_use = giottoNearestNetworks[i]$type,
        network_name = giottoNearestNetworks[i]$name)

      # spe stores in colpairs, with col indices instead of colnames
      cell1 <- match(nn_network$from, pDataDT(giottoObj)$cell_ID)
      cell2 <- match(nn_network$to, pDataDT(giottoObj)$cell_ID)

      SingleCellExperiment::colPair(spe, giottoNearestNetworks[i]$name) <- S4Vectors::SelfHits(
        cell1, cell2, nnode=ncol(spe), nn_network[, -1:-2]) #removing from and to
    }
  } else{
    message("No nearest networks found in the input Giotto object")
  }

  # spatial network from giotto
  giottoSpatialNetworks <- Giotto:::list_spatial_networks(giottoObj)
  if(!is.null(giottoSpatialNetworks)){
    message("Copying spatial networks")
    for(i in seq(nrow(giottoSpatialNetworks))){
      sp_network <- get_spatialNetwork(gobject = giottoObj, spat_unit = giottoSpatialNetworks[i]$spat_unit, name = giottoSpatialNetworks[i]$name)

      # spe stores in colpairs, with col indices instead of colnames
      cell1 <- match(sp_network$from, pDataDT(giottoObj)$cell_ID)
      cell2 <- match(sp_network$to, pDataDT(giottoObj)$cell_ID)

      SingleCellExperiment::colPair(spe, giottoSpatialNetworks[i]$name) <- S4Vectors::SelfHits(
        cell1, cell2, nnode=ncol(spe), sp_network[, -1:-2]) #removing from and to
    }
  } else{
    message("No spatial networks found in the input Giotto object")
  }

  # images from giotto?
  giottoImages <- Giotto:::list_images(giottoObj)
  if(!is.null(giottoImages)){
    message("Copying spatial images")
    for(i in seq(nrow(giottoImages))){
      img <- get_giottoImage(
        gobject = giottoObj,
        image_type = giottoImages[i]$img_type,
        name = giottoImages[i]$name)

      if(!is.null(img@file_path)){
        spe <- SpatialExperiment::addImg(spe,
                                         sample_id = "sample01", # ask how to find sample? different samples get appended to cell_ids
                                         image_id = img@name,
                                         imageSource = img@file_path,
                                         scaleFactor = NA_real_,
                                         load = TRUE)
      }
      else{
        message("\t - Skipping image with NULL file path")
      }
      S4Vectors::metadata(spe)[[img@name]] <- img
    }
  } else{
    message("No spatial images found in the input Giotto object")
  }

  # return SPE
  return(spe)
}



