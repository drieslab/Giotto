


#' @title anndataToGiotto
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


  if(!file.exists(anndata_path)) stop('path to anndata does not exist \n')
  adata <- sc$read(anndata_path)

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
