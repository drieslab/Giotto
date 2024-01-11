#' Perform library normalization on a dbSparseMatrix object
#' @name libNormDB
#' @details
#' Each column is divided by its sum and multiplied by the scalefactor.
#' @param dbMatrix tbl_duckdb_connection (required)
#' @param scalefactor numeric value to scale sparse matrix (required)
#'
#' @return dbMatrix
#' @export
#'
#' @examples
#' # using a dbSparseMatrix object
#' dbsm <- dbMatrix::sim_dbSparseMatrix()
#' scaled_dbsm <- lib_norm_db(dbMatrix = dbsm, scalefactor = 1000)
libNormDB <- function(dbMatrix, scalefactor) {
  # check inputs
  if (!is.numeric(scalefactor)) {
    stop("scalefactor must be numeric")
  }
  
  if (!inherits(dbMatrix, "dbMatrix")) {
    stop("Invalid `dbMatrix`. `dbMatrix` must be a `dbMatrix` object")
  }

  res <- dbMatrix[] |>
    dplyr::group_by(j) |>
    dplyr::mutate(x = (as.numeric(x) / sum(x, na.rm = TRUE)) * scalefactor) |>
    dplyr::ungroup() |>
    dplyr::collapse()

  dbMatrix[] <- res

  return(dbMatrix)
}
