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

#' Perform log normalization on a dbSparseMatrix object
#' @name logNormDB
#' @details
#' An offset is added to each non-zero value in dbSparseMatrix, the natural log is taken
#' and then divided by log(base).
#' @param dbMatrix tbl_duckdb_connection (required)
#' @param base base of log divisor to scale sparse matrix (required)
#' @param offset numeric value to add to each value in the sparse matrix (required)
#'
#' @return dbMatrix
#' @export
#'
#' @examples
#' # using a dbSparseMatrix object
#' dbsm <- dbMatrix::sim_dbSparseMatrix()
#' log_dbsm <- logNorm(dbMatrix = dbsm, base = 2, offset = 1000)
logNormDB <- function(dbMatrix, base, offset){
  # check inputs
  if(!is.numeric(base) | base == 0 | base == 1){
    stop("base must be numeric, non-zero, nor equal to 1")
  }

  if (!is.numeric(offset)) {
    stop("offset must be numeric")
  }

  if (!inherits(dbMatrix, "dbMatrix")) {
    stop("Invalid `dbMatrix`. `dbMatrix` must be a `dbMatrix` object")
  }

  res = dbMatrix[] |>
    dplyr::mutate(x = log(x + offset)/log(base)) |>
    dplyr::collapse()

  dbMatrix[] <- res

  return(dbMatrix)
}