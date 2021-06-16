
## sp classes ####




## anndata object ####

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





#' @name giottoToSeurat
#' @description Converts Giotto object into a Seurat object
#' @param obj_use Giotto object
#' @return Seurat object
#' @export
giottoToSeurat <- function(obj_use = NULL,...){
  require(Seurat)
  require(Giotto)
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
      sobj <- CreateSeuratObject(counts = data_raw,assay = assay_use)
      if ('normalized' %in% slot_use){
        sobj <- SetAssayData(sobj,slot = 'data',new.data = expr_use[['normalized']],assay = assay_use)
      }
      if ('scaled' %in% slot_use){
        sobj <- SetAssayData(sobj,slot = 'scale.data',new.data = expr_use[['scaled']],assay = assay_use)
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
        assay_obj <- CreateAssayObject(counts = data_raw)
      } else if (flag_raw == 0 & flag_norm == 1){
        assay_obj <- CreateAssayObject(data = data_norm)
      } else {
        stop(paste0('Raw and normalized data not found for assay ',assay_use))
      }
      sobj[[assay_use]] <- assay_obj
      if ('scaled' %in% slot_use){
        data_scale <- expr_use[['scaled']]
        sobj <- SetAssayData(sobj,slot = "scale.data",new.data = data_scale,assay = assay_use)
      }
    }

    # add cell metadata
    meta_cells <- as.data.frame(pDataDT(obj_use,feat_type = assay_use))
    rownames(meta_cells) <- meta_cells$cell_ID
    meta_cells <- meta_cells[,-which(colnames(meta_cells) == 'cell_ID')]
    colnames(meta_cells) <- paste0(assay_use,"_",colnames(meta_cells))
    sobj <- AddMetaData(sobj,metadata = meta_cells[Cells(sobj),])

    # add feature metadata
    meta_genes <- as.data.frame(fDataDT(obj_use,feat_type = assay_use))
    rownames(meta_genes) <- meta_genes$feat_ID
    sobj[[assay_use]]@meta.features <- cbind(sobj[[assay_use]]@meta.features,meta_genes)
  }

  # spatial coordinates
  loc_use <- as.data.frame(select_spatial_locations(gobject = obj_use,...))
  rownames(loc_use) <- loc_use$cell_ID
  sobj <- AddMetaData(sobj,metadata = loc_use)
  # add spatial coordinates as new dim reduct object
  loc_2 <- loc_use[,c('sdimx','sdimy')]
  colnames(loc_2) <- c('spatial_1','spatial_2')
  sobj[['spatial']] <- CreateDimReducObject(embeddings = as.matrix(loc_2),
                                            assay = names(sobj@assays)[1],key = 'spatial_')

  # dim reduction
  # note: Seurat requires assay name specification for each dim reduc
  if (!is.null(obj_use@dimension_reduction)){
    dr_use <- names(obj_use@dimension_reduction[[1]])
    for (i in dr_use){
      dr_name <- i
      dr_obj <- select_dimReduction(gobject = obj_use,return_dimObj = T,reduction_method = dr_name,name = dr_name)
      emb_use <- dr_obj$coordinates[Cells(sobj),]
      if (sum(c('loadings','eigenvalues') %in% names(dr_obj$misc)) == 2){
        loadings_use <- dr_obj$misc$loadings
        stdev_use <- dr_obj$misc$eigenvalues
        sobj[[dr_name]] <- CreateDimReducObject(embeddings = as.matrix(emb_use),loadings = loadings_use,
                                                key = dr_name,stdev = stdev_use,assay = names(sobj@assays)[1])
      } else {
        sobj[[dr_name]] <- CreateDimReducObject(embeddings = as.matrix(emb_use),
                                                key = dr_name,assay = names(sobj@assays)[1])
      }
    }
  }

  # network objects
  # expression network
  nn_all <- names(obj_use@nn_network)
  if (!is.null(nn_all)){
    for (i in nn_all){
      nn_use <- select_NearestNetwork(gobject = obj_use,nn_network_to_use = i,
                                      output = 'data.table',network_name = names(obj_use@nn_network[[i]]))
      idx1 <- match(nn_use$from,Cells(sobj))
      idx2 <- match(nn_use$to,Cells(sobj))
      edge_weight <- nn_use$weight
      nn_mtx <- Matrix::sparseMatrix(i = idx1,j = idx2,x = edge_weight,dims = c(ncol(sobj),ncol(sobj)))
      rownames(nn_mtx) <- colnames(nn_mtx) <- Cells(sobj)
      nn_name <- paste0('expr_',i)
      sobj[[nn_name]] <- as.Graph(nn_mtx)
      sobj[[nn_name]]@assay.used <- names(sobj@assays)[1]
    }
  }

  # spatial network
  sn_all <- names(obj_use@spatial_network)
  if (!is.null(sn_all)){
    for (i in sn_all){
      snt_use <- select_spatialNetwork(gobject = obj_use,name = i)
      idx1 <- match(snt_use$from,Cells(sobj))
      idx2 <- match(snt_use$to,Cells(sobj))
      edge_weight <- snt_use$weight
      nn_mtx <- Matrix::sparseMatrix(i = idx1,j = idx2,x = edge_weight,dims = c(ncol(sobj),ncol(sobj)))
      rownames(nn_mtx) <- colnames(nn_mtx) <- Cells(sobj)
      nn_name <- paste0('spatial_',i)
      sobj[[nn_name]] <- as.Graph(nn_mtx)
      sobj[[nn_name]]@assay.used <- names(sobj@assays)[1]
    }
  }

  return (sobj)
}
