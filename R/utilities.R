

## General utility functionss ##

#' @name mean_flex
#' @keywords internal
mean_flex = function(x, ...) {

  if(methods::is(x, 'dgCMatrix')) {
    return(Matrix::mean(x, ...)) # replace with sparseMatrixStats
  } else if(methods::is(x, 'Matrix')) {
    return(Matrix::mean(x, ...))
  } else {
    return(base::mean(x, ...))
  }
}

#' @name  rowSums_flex
#' @keywords internal
rowSums_flex = function(mymatrix) {

  if(methods::is(mymatrix, 'DelayedMatrix')) {
    return(DelayedMatrixStats::rowSums2(mymatrix))
  } else if(methods::is(mymatrix, 'dgCMatrix')) {
    return(Matrix::rowSums(mymatrix)) # replace with sparseMatrixStats
  } else if(methods::is(mymatrix, 'Matrix')) {
    return(Matrix::rowSums(mymatrix))
  } else {
    temp_matrix = as.matrix(mymatrix)
    temp_res = matrixStats::rowSums2(temp_matrix)
    names(temp_res) = rownames(temp_matrix)
    return(temp_res)
  }
}

#' @name rowMeans_flex
#' @keywords internal
rowMeans_flex = function(mymatrix) {

  # replace by MatrixGenerics?

  if(methods::is(mymatrix, 'DelayedMatrix')) {
    return(DelayedMatrixStats::rowMeans2(mymatrix))
  } else if(methods::is(mymatrix, 'dgCMatrix')) {
    return(Matrix::rowMeans(mymatrix)) # replace with sparseMatrixStats
  } else if(methods::is(mymatrix, 'Matrix')) {
    return(Matrix::rowMeans(mymatrix))
  } else {
    temp_matrix = as.matrix(mymatrix)
    temp_res = matrixStats::rowMeans2(temp_matrix)
    names(temp_res) = rownames(temp_matrix)
    return(temp_res)

  }
}

#' @name  colSums_flex
#' @keywords internal
colSums_flex = function(mymatrix) {

  if(methods::is(mymatrix, 'DelayedMatrix')) {
    return(DelayedMatrixStats::colSums2(mymatrix))
  } else if(methods::is(mymatrix, 'dgCMatrix')) {
    return(Matrix::colSums(mymatrix)) # replace with sparseMatrixStats
  } else if(methods::is(mymatrix, 'Matrix')) {
    return(Matrix::colSums(mymatrix))
  } else {
    temp_matrix = as.matrix(mymatrix)
    temp_res = matrixStats::colSums2(temp_matrix)
    names(temp_res) = colnames(temp_matrix)
    return(temp_res)
  }
}

#' @name  colMeans_flex
#' @keywords internal
colMeans_flex = function(mymatrix) {

  if(methods::is(mymatrix, 'DelayedMatrix')) {
    return(DelayedMatrixStats::colMeans2(mymatrix))
  } else if(methods::is(mymatrix, 'dgCMatrix')) {
    return(Matrix::colMeans(mymatrix)) # replace with sparseMatrixStats
  } else if(methods::is(mymatrix, 'Matrix')) {
    return(Matrix::colMeans(mymatrix))
  } else {
    temp_matrix = as.matrix(mymatrix)
    temp_res = matrixStats::colMeans2(temp_matrix)
    names(temp_res) = colnames(temp_matrix)
    return(temp_res)
  }
}

#' @name  t_flex
#' @keywords internal
t_flex = function(mymatrix) {

  if(methods::is(mymatrix, 'DelayedMatrix')) {
    return(t(mymatrix))
  } else if(methods::is(mymatrix, 'dgCMatrix')) {
    return(Matrix::t(mymatrix)) # replace with sparseMatrixStats
  } else if(methods::is(mymatrix, 'Matrix')) {
    return(Matrix::t(mymatrix))
  } else {
    mymatrix = as.matrix(mymatrix)
    mymatrix = base::t(mymatrix)
    return(mymatrix)
  }
}

#' @name  cor_flex
#' @keywords internal
cor_flex = function(x, ...) {
  x = as.matrix(x)
  return(stats::cor(x, ...))
}


#' @name  flex_lapply
#' @keywords internal
flex_lapply = function(X, FUN, cores = NA, fun = NULL, ...) {

  # a simple wrapper for future.apply::future_lapply
  # probably does not need any additional changes

  # potential addition:
  # check if future::plan() was already set by user
  # if not, set plan(multisession, workers = cores) by default


  # backwards compatible with previous version
  if(!is.null(fun)) {
    FUN = fun
  }


  # get type of os
  os = .Platform$OS.type

  # set number of cores automatically, but with limit of 10
  cores = determine_cores(cores)


  # future_lapply call
  save_list = future.apply::future_lapply(X = X, FUN = FUN, ...)

  #if(os == 'unix') {
  #  save_list = parallel::mclapply(X = X, mc.cores = cores,
  #                                 FUN = fun, ...)
  #} else if(os == 'windows') {
  #  save_list = parallel::mclapply(X = X, mc.cores = 1,
  #                                 FUN = fun, ...)
  #}

  return(save_list)
}

#' @name  my_arowMeans
#' @keywords internal
my_arowMeans = function(x) {
  if(is.null(nrow(x))) {
    x # if only one column is selected
    #mean(x)
  } else {
    rowMeans_flex(x)
  }
}

#' @name  my_growMeans
#' @keywords internal
my_growMeans = function(x, offset = 0.1) {
  if(is.null(nrow(x))) {
    x # if only one column is selected
    #exp(mean(log(x+offset)))-offset
  } else {
    exp(rowMeans_flex(log(x+offset)))-offset
  }
}

#' @name  my_rowMeans
#' @keywords internal
my_rowMeans = function(x, method = c('arithmic', 'geometric'), offset = 0.1) {
  method = match.arg(method, c('arithmic', 'geometric'))
  if(method == 'arithmic') return(my_arowMeans(x = x))
  if(method == 'geometric') return(my_growMeans(x = x, offset = offset))
}


#' @name standardise_flex
#' @description standardizes a matrix
#' @param x matrix
#' @param center center data
#' @param scale scale data
#' @keywords internal
#' @return standardized matrix
standardise_flex = function (x, center = TRUE, scale = TRUE) {

  if (center & scale) {
    y = Giotto:::t_flex(x) - matrixStats:::colMeans2(x)
    y = y/sqrt(matrixStats::rowSums2(y^2)) * sqrt((dim(x)[1] - 1))
    y = Giotto:::t_flex(y)
  }
  else if (center & !scale) {
    y = Giotto:::t_flex(x) - matrixStats:::colMeans2(x)
    y = Giotto:::t_flex(y)
  }
  else if (!center & scale) {
    csd = matrixStats::colSds(x)
    y = Giotto:::t_flex(Giotto:::t_flex(x) / csd )
  } else {
    y = x
  }
}

