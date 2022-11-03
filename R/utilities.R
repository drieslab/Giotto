

## General utility functions ##



#' @title mean_flex
#' @name mean_flex
#' @param x data to use
#' @param ... other arguments to pass
#' @keywords internal
mean_flex = function(x, ...) {

  if(inherits(x, 'dgCMatrix')) {
    return(Matrix::mean(x, ...)) # replace with sparseMatrixStats
  } else if(inherits(x, 'Matrix')) {
    return(Matrix::mean(x, ...))
  } else {
    return(base::mean(x, ...))
  }
}



#' @title rowSums_flex
#' @name rowSums_flex
#' @param mymatrix matrix to use
#' @keywords internal
rowSums_flex = function(mymatrix) {

  if(inherits(mymatrix, 'DelayedMatrix')) {
    return(DelayedMatrixStats::rowSums2(mymatrix))
  } else if(inherits(mymatrix, 'dgCMatrix')) {
    return(Matrix::rowSums(mymatrix)) # replace with sparseMatrixStats
  } else if(inherits(mymatrix, 'Matrix')) {
    return(Matrix::rowSums(mymatrix))
  } else {
    temp_matrix = as.matrix(mymatrix)
    temp_res = matrixStats::rowSums2(temp_matrix)
    names(temp_res) = rownames(temp_matrix)
    return(temp_res)
  }
}



#' @title rowMeans_flex
#' @name rowMeans_flex
#' @param mymatrix matrix to use
#' @keywords internal
rowMeans_flex = function(mymatrix) {

  # replace by MatrixGenerics?

  if(inherits(mymatrix, 'DelayedMatrix')) {
    return(DelayedMatrixStats::rowMeans2(mymatrix))
  } else if(inherits(mymatrix, 'dgCMatrix')) {
    return(Matrix::rowMeans(mymatrix)) # replace with sparseMatrixStats
  } else if(inherits(mymatrix, 'Matrix')) {
    return(Matrix::rowMeans(mymatrix))
  } else {
    temp_matrix = as.matrix(mymatrix)
    temp_res = matrixStats::rowMeans2(temp_matrix)
    names(temp_res) = rownames(temp_matrix)
    return(temp_res)

  }
}



#' @title colSums_flex
#' @name colSums_flex
#' @param mymatrix matrix to use
#' @keywords internal
colSums_flex = function(mymatrix) {

  if(inherits(mymatrix, 'DelayedMatrix')) {
    return(DelayedMatrixStats::colSums2(mymatrix))
  } else if(inherits(mymatrix, 'dgCMatrix')) {
    return(Matrix::colSums(mymatrix)) # replace with sparseMatrixStats
  } else if(inherits(mymatrix, 'Matrix')) {
    return(Matrix::colSums(mymatrix))
  } else {
    temp_matrix = as.matrix(mymatrix)
    temp_res = matrixStats::colSums2(temp_matrix)
    names(temp_res) = colnames(temp_matrix)
    return(temp_res)
  }
}



#' @title colMeans_flex
#' @name colMeans_flex
#' @param mymatrix matrix to use
#' @keywords internal
colMeans_flex = function(mymatrix) {

  if(inherits(mymatrix, 'DelayedMatrix')) {
    return(DelayedMatrixStats::colMeans2(mymatrix))
  } else if(inherits(mymatrix, 'dgCMatrix')) {
    return(Matrix::colMeans(mymatrix)) # replace with sparseMatrixStats
  } else if(inherits(mymatrix, 'Matrix')) {
    return(Matrix::colMeans(mymatrix))
  } else {
    temp_matrix = as.matrix(mymatrix)
    temp_res = matrixStats::colMeans2(temp_matrix)
    names(temp_res) = colnames(temp_matrix)
    return(temp_res)
  }
}



#' @title t_flex
#' @name t_flex
#' @param mymatrix matrix to use
#' @keywords internal
t_flex = function(mymatrix) {

  if(inherits(mymatrix, 'DelayedMatrix')) {
    return(t(mymatrix))
  } else if(inherits(mymatrix, 'dgCMatrix')) {
    return(Matrix::t(mymatrix)) # replace with sparseMatrixStats
  } else if(inherits(mymatrix, 'Matrix')) {
    return(Matrix::t(mymatrix))
  } else {
    mymatrix = as.matrix(mymatrix)
    mymatrix = base::t(mymatrix)
    return(mymatrix)
  }
}



#' @title cor_flex
#' @name cor_flex
#' @param x data to use
#' @param ... other arguments passed to stats::cor()
#' @keywords internal
cor_flex = function(x, ...) {
  x = as.matrix(x)
  return(stats::cor(x, ...))
}



#' @title lapply_flex
#' @name lapply_flex
#' @param X list to use
#' @param FUN function to be performed
#' @param cores cores to use
#' @param fun deprecated. Backwards compatibility for FUN
#' @param ... other arguments to pass
#' @keywords internal
lapply_flex = function(X,
                       FUN,
                       cores = NA,
                       fun = NULL,
                       ...) {

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


#' @title flex_lapply
#' @name  flex_lapply
#' @inheritDotParams lapply_flex
#' @seealso \code{\link{lapply_flex}}
#' @keywords internal
flex_lapply = function(X,
                       FUN,
                       cores = NA,
                       fun = NULL,
                       ...) {


  .Deprecated(new = "lapply_flex")

  lapply_flex(...)

}



#' @title my_arowMeans
#' @name  my_arowMeans
#' @param x data to use
#' @keywords internal
my_arowMeans = function(x) {
  if(is.null(nrow(x))) {
    x # if only one column is selected
    #mean(x)
  } else {
    rowMeans_flex(x)
  }
}



#' @title my_growMeans
#' @name  my_growMeans
#' @param x data to use
#' @param offset offset
#' @keywords internal
my_growMeans = function(x, offset = 0.1) {
  if(is.null(nrow(x))) {
    x # if only one column is selected
    #exp(mean(log(x+offset)))-offset
  } else {
    exp(rowMeans_flex(log(x+offset)))-offset
  }
}

#' @title my_rowMeans
#' @name  my_rowMeans
#' @param x data to use
#' @param method method is either "arithmic" or "geometric"
#' @param offset offset
#' @keywords internal
my_rowMeans = function(x, method = c('arithmic', 'geometric'), offset = 0.1) {
  method = match.arg(method, c('arithmic', 'geometric'))
  if(method == 'arithmic') return(my_arowMeans(x = x))
  if(method == 'geometric') return(my_growMeans(x = x, offset = offset))
}



#' @title standardise_flex
#' @name standardise_flex
#' @description standardizes a matrix
#' @param x matrix
#' @param center center data
#' @param scale scale data
#' @keywords internal
#' @return standardized matrix
standardise_flex = function (x, center = TRUE, scale = TRUE) {

  if (center & scale) {
    y = t_flex(x) - colMeans_flex(x)
    y = y/sqrt(rowSums_flex(y^2)) * sqrt((dim(x)[1] - 1))
    y = t_flex(y)
  }
  else if (center & !scale) {
    y = t_flex(x) - colMeans_flex(x)
    y = t_flex(y)
  }
  else if (!center & scale) {
    csd = matrixStats::colSds(x)
    y = t_flex(t_flex(x) / csd )
  } else {
    y = x
  }
}



#' @title Test if list element exists
#' @name list_element_exists
#' @description Test if nth element of list exists
#' @param x list
#' @param index element index
#' @keywords internal
#' @return boolean
list_element_exists = function(x, index) {
  tryCatch({
    if(length(x[[index]]) > -1)
      return(TRUE)
  }, error = function(e) {
    return(FALSE)
  })
}



# methods/operator imports ####


# Additional methods functions that might also be imported later...

#  getSlots slotNames .hasSlot
# Usage tag content:

#' .hasSlot(object, name)
#'
#' slotNames(x)
#' getSlots(x)
#' slot(object, name, check = TRUE) <- value


#' See \code{methods::\link[methods]{slot}} for details.
#'
#' @rdname slot
#' @name slot
#' @keywords internal
#' @export
#' @importFrom methods slot slot<-
#' @usage slot(object, name)
NULL


#' See \code{methods::\link[methods]{validObject}} for details.
#'
#' @rdname validObject
#' @name validObject
#' @keywords internal
#' @export
#' @importFrom methods validObject
#' @usage validObject(object, test = FALSE, complete = FALSE)
NULL



# Show function helpers ####

#' @title Hierarchical tree printing
#' @name print_leaf
#' @param level_index Which col of availability matrix to start print from
#' @param availableDT availability matrix given as data.table
#' @param inherit_last (boolean) determine behavior from previous level for last level (intended for values print)
#' Only TRUE behavior defined. #TODO
#' @param indent ident characters to print for this level (top level is '')
#' @keywords internal
#' @details Much inspiration taken from https://rdrr.io/cran/fs/src/R/tree.R
print_leaf = function(level_index,
                      availableDT,
                      inherit_last = TRUE,
                      indent) {


  ch = box_chars()

  leafs = unique(unlist(availableDT[, level_index, with = FALSE])) # Read unique items for this level into 'leafs'
  for (i in seq_along(leafs)) {
    if(isTRUE(inherit_last) & level_index == ncol(availableDT)) { # values layer has no special indent and ends
      writeLines(paste0(indent, capture.output(cat(leafs[[i]]))))
      cat(indent ,"\n", sep = "")
    } else {

      if (i == length(leafs)) {
        cat(indent, ch$b, leafs[[i]], "\n", sep = "")
        print_leaf(level_index = level_index + 1, # increment level_index
                   availableDT = availableDT[as.vector(availableDT[, level_index, with = FALSE] == leafs[[i]])], # pass subset
                   inherit_last = inherit_last,
                   indent = paste0(indent, "   "))
      } else {
        cat(indent, ch$t, leafs[[i]], "\n", sep = "")
        print_leaf(level_index = level_index + 1,
                   availableDT = availableDT[as.vector(availableDT[, level_index, with = FALSE] == leafs[[i]])], # pass subset
                   inherit_last = inherit_last,
                   indent = paste0(indent, paste0(ch$v, "  ")))

      }
    }
  }
}


#' @title Box characters
#' @name box_chars
#' @keywords internal
#' @details Much inspiration taken from https://rdrr.io/cran/fs/src/R/tree.R
# These are derived from https://github.com/r-lib/cli/blob/e9acc82b0d20fa5c64dd529400b622c0338374ed/R/tree.R#L111
box_chars <- function() {
  list(
    "h" = "\u2500",                   # horizontal
    "v" = "\u2502",                   # vertical
    "l" = "\u2514",
    "j" = "\u251C",
    "b" = "\u2514\u2500\u2500",       # branch
    "t" = "\u251C\u2500\u2500",       # T
    "i" = "\u2502  ",                 # layer
    "s" = "   "                       # spaces
  )
}


# Print Formatting ####

#' @title Wrap message
#' @name wrap_msg
#' @keywords internal
#' @param ... additional strings and/or elements to pass to cat
#' @param collapse how to join elements of string (default is no space)
wrap_msg = function(..., collapse = '') {
  cat(...) %>%
    capture.output() %>%
    strwrap() %>%
    paste(., collapse = '\n') %>%
    message()
}





