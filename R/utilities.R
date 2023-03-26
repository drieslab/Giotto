

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
#' @include generics.R
#' @keywords internal
t_flex = function(mymatrix) {

  if(inherits(mymatrix, 'DelayedMatrix')) {
    return(t(mymatrix))
  } else if(inherits(mymatrix, 'dgCMatrix')) {
    return(Matrix::t(mymatrix)) # replace with sparseMatrixStats
  } else if(inherits(mymatrix, 'Matrix')) {
    return(Matrix::t(mymatrix))
  } else if(inherits(mymatrix, 'spatLocsObj')){
    return(t(mymatrix))
  } else if(inherits(mymatrix, 'spatialNetworkObj')) {
    return(t(mymatrix))
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
#' @param future.seed whether to set a seed when using future_lapply
#' @param fun deprecated. Backwards compatibility for FUN
#' @param ... other arguments to pass
#' @keywords internal
lapply_flex = function(X,
                       FUN,
                       cores = NA,
                       future.seed = TRUE,
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
  save_list = future.apply::future_lapply(X = X, FUN = FUN, future.seed = future.seed, ...)

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

#' @title standardise_db
#' @name standardise_db
#' @description standardizes a database table
#' @param x expression database tbl
#' @param center center data
#' @param scale scale data
#' @param by feat or cell
#' @keywords internal
#' @return standardized database table
standardise_db = function (x, center = TRUE, scale = TRUE, by = c('feat', 'cell')) {
  
  i = match.arg(by, c('feat', 'cell')) 
  
  if (center & scale) {
    colMeans <- x |>
       group_by(i) |> #scale feats
       summarise(mean_i = mean(xlognorm, na.rm = TRUE), .groups = "drop")

    centered_expr <- x |>
       inner_join(colMeans, by = i) |>
       group_by(i) |> #scale feats
       mutate(xlognorm_centered = (xlognorm - mean_i)) |>
       select('feat', 'cell', xlognorm_centered) |>
       ungroup()

    scaled_expr <- centered_expr |>
       group_by(i) |> # scale feats
       summarise(sqrtsumsq = sqrt(sum(xlognorm_centered^2, na.rm = TRUE)), #sqrt(rowSums_flex(y^2))
                 sqrtNrow = sqrt(n() - 1), # sqrt((dim(x)[1] - 1))
                 .groups = "drop") |>
       ungroup()
    
    res <- centered_expr |>
       inner_join(scaled_expr, by = i) |>
       group_by(i) |> #scale feats
       mutate(xlognorm_centered_scaled = (xlognorm_centered / sqrtsumsq * sqrtNrow)) |>
       select('feat', 'cell', xlognorm_centered_scaled)

  } else if (center & !scale) {
    colMeans <- x |>
       group_by(i) |> #scale feats
       summarise(mean_i = mean(xlognorm, na.rm = TRUE), .groups = "drop")

    res <- x |>
       inner_join(colMeans, by = i) |>
       group_by(i) |> #scale feats
       mutate(xlognorm_centered = (xlognorm - mean_i)) |>
       select('feat', 'cell', xlognorm_centered) |>
       ungroup()
    
  } else if (!center & scale) {
    scaled_expr <- x |>
       group_by(i) |> # scale feats
       summarise(sqrtsumsq = sqrt(sum(xlognorm_centered^2, na.rm = TRUE)), #sqrt(rowSums_flex(y^2))
                 sqrtNrow = sqrt(n() - 1), # sqrt((dim(x)[1] - 1))
                 .groups = "drop") |>
       ungroup()

    res <- x |>
          inner_join(scaled_expr, by = i) |>
          group_by(i) |> #scale feats
          mutate(xlognorm_centered_scaled = (xlognorm_centered / sqrtsumsq * sqrtNrow)) |>
          select('feat', 'cell', xlognorm_centered_scaled)
  } else {
    res = x
  }

  return(res)
}


#' @title Test if list element exists
#' @name list_element_exists
#' @description Test if nth element of list exists
#' @param x list
#' @param index element index
#' @keywords internal
#' @return boolean
#' @noRd
list_element_exists = function(x, index) {
  tryCatch({
    if(length(x[[index]]) > -1)
      return(TRUE)
  }, error = function(e) {
    return(FALSE)
  })
}



# Based on https://stackoverflow.com/questions/37878620/reorder-rows-in-data-table-in-a-specific-order
#' @title Set specific data.table row order
#' @param x data.table
#' @param neworder numerical vector to reorder rows
#' @keywords internal
#' @noRd
set_row_order_dt = function(x, neworder) {
  if('.r' %in% colnames(x)) {
    temp_r = x[, .SD, .SDcols = '.r']
    data.table::setorderv(temp_r[, eval(call(":=", as.name(".r_alt"), call("order", neworder)))], ".r_alt")[, ".r_alt" := NULL]
    data.table::setorderv(x[, eval(call(":=", as.name(".r"), call("order", neworder)))], ".r")[, ".r" := NULL]
    x[, eval(call(':=', as.name('.r'), temp_r$.r))]
  } else {
    data.table::setorderv(x[, eval(call(":=", as.name(".r"), call("order", neworder)))], ".r")[, ".r" := NULL]
  }

}



# Determine if character is missing, NULL or an empty string.
# Returns TRUE for all cases
is_empty_char = function(x) {
  if(is.null(x)) return(TRUE)
  if(is.na(x)) return(TRUE)
  if(x == '') return(TRUE)

  FALSE
}







# list nesting depth ####

#' @title Find depth of subnesting
#' @name depth
#' @param this object to evaluate
#' @param method max (default) or min nesting to detect
#' @param sig signature or class to check for. Default is 'data.frame'
#' @description Recursively determines how many max or min layers of subnesting
#' there is, with the end object (defined by param sig or a list of length 0)
#' being layer 0
#' @details https://stackoverflow.com/questions/13432863/determine-level-of-nesting-in-r
#' @keywords internal
depth <- function(this,
                  method = c('max', 'min'),
                  sig = 'data.frame') {

  method = match.arg(arg = method, choices = c('max', 'min'))

  # Stop conditions:

  # Stop if matches signature to search for
  if(inherits(this, sig)) {
    return(0L)
  }
  # Stop if an empty list is discovered
  if(inherits(this, 'list') && length(this) == 0L) {
    return(0L)
  }
  # Stop if object is not a list AND recurse if it is.
  # Report minimum or maximum depth depending on method
  if(method == 'max') {
    ifelse(inherits(this, 'list'), 1L + max(sapply(this, function(x) depth(x, method = method, sig = sig))), 0L)
  } else if(method == 'min') {
    ifelse(inherits(this, 'list'), 1L + min(sapply(this, function(x) depth(x, method = method, sig = sig))), 0L)
  }

}








# Import from methods ####

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


#' See \code{methods::\link[methods]{slotNames}} for details.
#'
#' @rdname slotNames
#' @name slotNames
#' @keywords internal
#' @importFrom methods slotNames
#' @usage slotNames(object)
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
#' @description Helper function to print unicode characters using escape codes.
#' @keywords internal
#' @details Much inspiration taken from \pkg{fs} \href{https://rdrr.io/cran/fs/src/R/tree.R}{tree.R}
# These are derived from https://github.com/r-lib/cli/blob/e9acc82b0d20fa5c64dd529400b622c0338374ed/R/tree.R#L111
box_chars = function() {
  if(is_utf8_output()) {
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
  } else {
    list(
      "h" = "-",                        # horizontal
      "v" = "|",                        # vertical
      "l" = "\\",
      "j" = "+",
      "b" = "\\--",                     # branch
      "t" = "+--",                      # T
      "i" = "|  ",                      # layer
      "s" = "   "                       # spaces
    )
  }
}



#' @describeIn box_chars Determine if print is latex output
#' @keywords internal
is_latex_output = function() {
  if(!('knitr' %in% loadedNamespaces())) return(FALSE)
  get('is_latex_output', asNamespace('knitr'))()
}



#' @describeIn box_chars Determine if system is using UTF-8 encoding
#' @keywords internal
is_utf8_output = function() {
  opt = getOption('cli.unicode', default = NULL)
  if(!is.null(opt)) return(isTRUE(opt))
  opt = getOption('giotto.unicode', default = NULL)
  if(!is.null(opt)) return(isTRUE(opt))

  is_utf8 = (l10n_info()$`UTF-8` & !is_latex_output())
  options('giotto.unicode' = is_utf8)
  return(is_utf8)

}





#' @title Print abbreviated matrix
#' @name abb_mat
#' @description print abbreviated matrix exprObj. Works for Matrix pkg denseMatrix,
#' matrix, data.frame and classes that inherit them.
#' @keywords internal
abb_mat = function(exprObj, nrows, ncols, header = TRUE) {
  mat = as.matrix(exprObj[])
  four_names = head(colnames(mat), 4)
  mat_cols = ncol(mat)
  mat_rows = nrow(mat)

  # suppress colnames
  mat = mat[1:if(nrows <= mat_rows) nrows else mat_rows,
            1:if(ncols <= mat_cols) ncols else mat_cols]
  colnames(mat) = NULL

  # prints
  if(isTRUE(header)) {
    cat('An object of class', class(exprObj), '\n')
    cat(paste0('for spatial unit: "', exprObj@spat_unit, '" and feature type: "', exprObj@feat_type, '"\n'))
  }
  cat('  Provenance:', exprObj@provenance)
  if(isTRUE(header)) cat('\n\ncontains:\n')
  else cat('\n')
  cat(paste0(mat_rows, ' x ', mat_cols, ' dense matrix of class "', class(exprObj[]), '"\n\n'))
  print(mat)
  cat('\n First four colnames:')
  cat('\n', wrap_txt(four_names, strWidth = 40), '\n')
}



#' @title Print abbreviated spatlocs
#' @name abb_spatlocs
#' @description print abbreviated spatLocsObj
#' @keywords internal
abb_spatlocs = function(spatLocsObj, nrows) {
  # data.table vars
  sdimx = sdimy = NULL

  DT_rows = nrow(spatLocsObj[])
  spatlocs = spatLocsObj[][1:if(nrows <= DT_rows) nrows else DT_rows,]

  # prints
  cat('An object of class', class(spatLocsObj), '\n')
  cat('provenance:', slot(spatLocsObj, 'provenance'))
  cat('\n    ------------------------\n')
  print(spatlocs)
  cat('\nranges:\n')
  try(expr = print(sapply(slot(spatLocsObj, 'coordinates')[,.(sdimx,sdimy)], range)),
      silent = TRUE)

  cat('\n\n')
}



# Print Formatting ####

#' @title Wrap message
#' @name wrap_msg
#' @param ... additional strings and/or elements to pass to wrap_txt
#' @param sep how to join elements of string (default is one space)
#' @keywords internal
wrap_msg = function(..., sep = ' ') {
  message(wrap_txt(..., sep = sep))
}

#' @title Wrap text
#' @name wrap_txt
#' @param ... additional params to pass
#' @param sep how to join elements of string (default is one space)
#' @param strWidth externally set wrapping width. (default value of 100 is not effected)
#' @param errWidth default = FALSE. Set strWidth to be compatible with error printout
#' @keywords internal
wrap_txt = function(..., sep = ' ', strWidth = 100, errWidth = FALSE) {
  custom_width = ifelse(is.null(match.call()$strWidth), yes = FALSE, no = TRUE)
  if(!isTRUE(custom_width)) {
    if(isTRUE(errWidth)) strWidth = getOption('width') - 6
  }

  cat(..., sep = sep) %>%
    capture.output() %>%
    strwrap(., prefix =  ' ', initial = '', # indent later lines, no indent first line
            width = min(80, getOption("width"), strWidth)) %>%
    paste(., collapse = '\n')
}




# Color text (8 colors) ####

#' @title Colorize print text
#' @name color_tag
#' @details supported colors checking is modified from \pkg{cli}
#' \href{https://github.com/r-lib/cli/blob/HEAD/R/num-ansi-colors.R}{aab-num-ansi-colors.R}
#' @keywords internal
color_tag = function() {
  list(
    r = '\u001b[31m', # red
    g = '\u001b[32m', # green
    y = '\u001b[33m', # yellow
    b = '\u001b[34m', # blue
    p = '\u001b[35m', # purple
    t = '\u001b[36m', # teal
    x = '\u001b[39m'  # none (return)
  )
}



#' @describeIn color_tag Determine if system should print color
#' @keywords internal
use_color_text = function() {
  opt = getOption('giotto.color_show', default = NULL)
  ansi8_color = ansi_colors() >= 8L
  if(!is.null(opt)) {
    if(!isTRUE(opt)) return(opt)
    if(isTRUE(opt) & isTRUE(ansi8_color)) return(opt)
    if(isTRUE(opt) & !isTRUE(ansi8_color)) {
      wrap_msg('Color text not supported on this system.
               Set options("giotto.color_show" = FALSE)')
    }
  } else ansi8_color
}



#' @describeIn color_tag Determine if system can print at least 8 colors
#' @keywords internal
ansi_colors = function() {

  # options
  opt = getOption('cli.num_colors', default = NULL)
  if(!is.null(opt)) return(as.integer(opt))
  opt = getOption('giotto.num_colors', default = NULL)
  if(!is.null(opt)) return(as.integer(opt))

  if((env = Sys.getenv('R_CLI_NUM_COLORS', '')) != '') {
    return(as.integer(env))
  }

  # crayon compatibility (allow color disabling through crayon)
  cray_opt_has = getOption('crayon.enabled', NULL)
  cray_opt_num = getOption('crayon.colors', NULL)
  if(!is.null(cray_opt_has) & !isTRUE(cray_opt_has)) return(1L) # disable
  if(isTRUE(cray_opt_has) & !is.null(cray_opt_num)) {
    return(as.integer(cray_opt_num))
  }
  if(isTRUE(cray_opt_has) & is.null(cray_opt_num)) {
    return(8L)
  }

  # 'NO_COLOR env setting disabling
  if(!is.na(Sys.getenv('NO_COLOR', NA_character_))) return(1L)

  # if knitr then no color in .Rmd chunks
  if(isTRUE(getOption('knitr.in.progress'))) return(1L)

  if(.Platform$GUI == 'AQUA') return(1L)

  # No specific usage cases needed
  # Colors only used in show functions
  if(Sys.getenv('RSTUDIO') == 1) return(8L) # at least

  # Must be placed after 'RSTUDIO' check
  if(.Platform$GUI == 'Rgui') return(1L)

  # Windows Emacs
  if(.Platform$OS.type == 'windows' &
     '--ess' %in% commandArgs() &
     is_emacs_with_color()) {
    return(8L)
  }

  # end catch
  return(8L)

}



#' @describeIn color_tag Determine if emacs can print color
#' @keywords internal
is_emacs_with_color = function() {
  (Sys.getenv('EMACS') != '' || Sys.getenv('INSIDE_EMACS') !=
     '') & !is.na(emacs_version()[1]) & emacs_version()[1] >=
    23
}



#' @describeIn color_tag Determine emacs version
#' @keywords internal
emacs_version = function() {
  ver <- Sys.getenv('INSIDE_EMACS')
  ver <- gsub('[^0-9\\.]+', '', ver)
  if (ver == '')
    return(NA_integer_)
  ver <- strsplit(ver, '.', fixed = TRUE)[[1]]
  as.numeric(ver)
}


# time formatting ####

#' @title Format time for printing
#' @name time_format
#' @keywords internal
#' @details Code from \code{\link[data.table]{timetaken}}
time_format = function(secs) {
  if(secs > 60) {
    secs = as.integer(secs)
    sprintf("%02d:%02d:%02d", secs%/%3600L, (secs%/%60L)%%60L,
            secs%%60L)
  }
  else {
    sprintf(if(secs >= 10)
      "%.1fs"
      else "%.3fs", secs)
  }
}

# radians and degrees ####
#' @title Radian/degree conversions
#' @name degrees
#' @description Convert radians to degrees and vice versa
#' @param deg degrees
#' @param rad radians
#' @keywords internal
NULL

#' @describeIn degrees Degrees to radians
radians = function(deg) {
  deg * pi / 180
}

#' @describeIn degrees Radians to degrees
#' @keywords internal
degrees = function(rad) {
  rad * 180 / pi
}















# guard functions ####

#' @inheritParams data_access_params
#' @param n Frames back in which to evaluate the gobject param
#' @keywords internal
#' @noRd
guard_against_notgiotto = function(gobject, n = 1L, ...) {
  fn_name = deparse(sys.calls()[[sys.nframe() - n]])
  orig_name = deparse(eval(call('substitute', as.name(substitute(gobject)), parent.frame())))
  if(!hasArg(gobject)) stop(wrap_txt(fn_name, ':\ngiotto object must be given',
                                     errWidth = TRUE),
                            call. = FALSE)
  if(!inherits(gobject, 'giotto')) {
    stop(wrap_txt(fn_name, ':\n', orig_name, 'is not a giotto object',
                  errWidth = TRUE),
         call. = FALSE)
  }
}
















