## CLASS ####
# ------- ###


setClass(
  "VisiumHDReader",
  slots = list(
    visiumHD_dir = "character",
    expression_source = "character",
    gene_column_index = "numeric",
    barcodes = "character",
    array_subset_row = "numeric",
    array_subset_col = "numeric",
    pxl_subset_row = "numeric",
    pxl_subset_col = "numeric",
    calls = "list"
  ),
  prototype = list(
    expression_source = 'raw',
    gene_column_index = 2,
    barcodes = NULL,
    array_subset_row = NULL,
    array_subset_col = NULL,
    pxl_subset_row = NULL,
    pxl_subset_col = NULL,
    calls = list()
  )
)



# * show ####
setMethod("show", signature("VisiumHDReader"), function(object) {
  cat(sprintf("Giotto <%s>\n", "VisiumHDReader"))
  print_slots <- c("dir", "expression_source", "gene_column_index",
                   "barcodes", "array_subset_row", "array_subset_col",
                   "pxl_subset_row", "pxl_subset_col",
                   "funs")
  pre <- sprintf(
    "%s :", format(print_slots)
  )
  names(pre) <- print_slots
  
  # dir
  d <- object@visiumHD_dir
  if (length(d) > 0L) {
    nch <- nchar(d)
    d <- abbrev_path(d)
    cat(pre["dir"], d, "\n")
  } else {
    cat(pre["dir"], "\n")
  }
  
  # expression_source
  expression_source <- object@expression_source
  cat(pre["expression_source"], expression_source, "\n")
  
  # gene_column_index
  gene_column_index <- object@gene_column_index
  cat(pre["gene_column_index"], gene_column_index, "\n")
  
  # barcodes
  barcodes <- ifelse(!is.null(object@barcodes), "found", "none")
  cat(pre["barcodes"], barcodes, "\n")
  
  # array_subset_row
  array_subset_row <- ifelse(!is.null(object@array_subset_row), "found", "none")
  cat(pre["array_subset_row"], array_subset_row, "\n")
  
  # array_subset_col
  array_subset_col <- ifelse(!is.null(object@array_subset_col), "found", "none")
  cat(pre["array_subset_col"], array_subset_col, "\n")
  
  # pxl_subset_row
  pxl_subset_row <- ifelse(!is.null(object@pxl_subset_row), "found", "none")
  cat(pre["pxl_subset_row"], pxl_subset_row, "\n")
  
  # pxl_subset_col
  pxl_subset_col <- ifelse(!is.null(object@pxl_subset_col), "found", "none")
  cat(pre["pxl_subset_col"], pxl_subset_col, "\n")
  
  # funs
  .reader_fun_prints(x = object, pre = pre["funs"])
})

# * print ####
setMethod("print", signature("VisiumHDReader"), function(x, ...) show(x))



#' @title Import a Visium HD assay
#' @name importVisiumHD
#' @description
#' Giotto import functionalities for Visium HD datasets. This function generates
#' a `VisiumHDReader` instance that has convenient reader functions for converting
#' individual pieces of Visium HD data into Giotto-compatible representations when
#' the param `visiumHD_dir` is provided.
#' A function that creates the full `giotto` object is also available.
#' These functions should have all param values provided as defaults, but
#' can be flexibly modified to do things such as look in alternative
#' directories or paths.
#' @param visiumHD_dir Visium HD output directory (e.g. square_016um)
#' @param expression_source character. Raw or filter expression data. Defaults to 'raw'
#' @param gene_column_index numeric. Expression column to use for gene names 
#' 1 = Ensembl and 2 = gene symbols
#' @param barcodes character vector. (optional) Use if you only want to load
#' a subset of the pixel barcodes
#' @param array_subset_row numeric vector. (optional) Vector with min and max values
#' to subset based on array rows
#' @param array_subset_col numeric vector. (optional) Vector with min and max values
#' to subset based on array columns
#' @param pxl_subset_row numeric vector. (optional) Vector with min and max values
#' to subset based on row pixels
#' @param pxl_subset_col numeric vector. (optional) Vector with min and max values
#' to subset based on column pixels
#' @details
#' Loading functions are generated after the `visiumHD_dir` is added.
#' @returns VisiumHDReader object
#' @examples
#' # Create a `VisiumHDReader` object
#' reader <- importVisiumHD()
#'
#' \dontrun{
#' # Set the visiumHD_dir
#' reader$visiumHD_dir <- "path to visium HD dir"
#' readerHD$visiumHD_dir <- visiumHD_dir
#' 
#' # Load tissue positions or create cell metadata
#' tissue_pos = readerHD$load_tissue_position()
#' metadata <- readerHD$load_metadata()
#' 
#' Load matrix or create expression object
#' matrix <- readerHD$load_matrix()
#' expression_obj = readerHD$load_expression()
#' 
#' Load transcript data (cell metadata, expression object, and transcripts per pixel)
#' my_transcripts = readerHD$load_transcripts(array_subset_row = c(500, 1000),
#'                                            array_subset_col = c(500, 1000))
#'                                            
#' # Create a `giotto` object and add the loaded data
#' TODO
#' }
#' @export
importVisiumHD <- function(
    visiumHD_dir = NULL, 
    expression_source = 'raw', 
    gene_column_index = 2,
    barcodes = NULL,
    array_subset_row = NULL,
    array_subset_col = NULL,
    pxl_subset_row = NULL,
    pxl_subset_col = NULL) {
  
  # get params
  a <- list(Class = "VisiumHDReader")
  
  if (!is.null(visiumHD_dir)) {
    a$visiumHD_dir <- visiumHD_dir
  }
  
  a$expression_source <- expression_source
  a$gene_column_index <- gene_column_index
  
  if (!is.null(barcodes)) {
    a$barcodes <- barcodes
  }
  
  if (!is.null(array_subset_row)) {
    a$array_subset_row <- array_subset_row
  }
  
  if (!is.null(array_subset_col)) {
    a$array_subset_col <- array_subset_col
  }
  
  if (!is.null(pxl_subset_row)) {
    a$pxl_subset_row <- pxl_subset_row
  }
  
  if (!is.null(pxl_subset_col)) {
    a$pxl_subset_col <- pxl_subset_col
  }
  
  do.call(new, args = a)
}


# * init ####
setMethod("initialize", signature("VisiumHDReader"), function(
    .Object, visiumHD_dir, 
    expression_source, 
    gene_column_index, 
    barcodes, 
    array_subset_row,
    array_subset_col,
    pxl_subset_row,
    pxl_subset_col
) {
  
  # provided params (if any)
  if (!missing(visiumHD_dir)) {
    checkmate::assert_directory_exists(visiumHD_dir)
    .Object@visiumHD_dir <- visiumHD_dir
  }
  
  if (!missing(expression_source)) {
    .Object@expression_source <- expression_source
  }
  
  if (!missing(gene_column_index)) {
    .Object@gene_column_index <- gene_column_index
  }
  
  if (!missing(barcodes)) {
    .Object@barcodes <- barcodes
  }
  
  if (!missing(array_subset_row)) {
    .Object@array_subset_row <- array_subset_row
  }
  
  if (!missing(array_subset_col)) {
    .Object@array_subset_col <- array_subset_col
  }
  
  if (!missing(pxl_subset_row)) {
    .Object@pxl_subset_row <- pxl_subset_row
  }
  
  if (!missing(pxl_subset_col)) {
    .Object@pxl_subset_col <- pxl_subset_col
  }
  
  # NULL case
  if (length(.Object@visiumHD_dir) == 0) {
    return(.Object) # return early if no path given
  }
  
  
  # detect paths and subdirs
  p <- .Object@visiumHD_dir
  
  
  .visiumHD_detect <- function(pattern, path = p, recursive = FALSE) {
    .detect_in_dir(pattern = pattern, path = path, recursive = recursive, platform = "visiumHD")
  }
  
  
  filter_expr_dir <- .visiumHD_detect(pattern = "filtered_feature_bc_matrix", path = p)
  raw_expr_dir    <- .visiumHD_detect(pattern = "raw_feature_bc_matrix", path = p)
  
  s <- .Object@expression_source
  if(s == 'raw') {
    expr_dir = raw_expr_dir
  } else if(s == 'filter') {
    expr_dir = filter_expr_dir
  } else {
    stop('expression source for visiumHD can only be raw or filter')
  }
  
  spatial_dir <- .visiumHD_detect(pattern = "spatial", path = p)
  
  
  c_index <- .Object@gene_column_index
  if(!c_index %in% c(1, 2)) {
    stop('gene column index can only be 1 (Ensembl) or 2 (gene symbols)')
  }
  
  
  ## matrix load call
  matrix_fun <- function(
    path = expr_dir,
    gene_column_index = c_index,
    remove_zero_rows = TRUE,
    split_by_type = TRUE,
    verbose = NULL
  ) {
    .visiumHD_matrix(
      path = path,
      gene_column_index = gene_column_index,
      remove_zero_rows = remove_zero_rows,
      split_by_type = split_by_type,
      verbose = verbose
    )
  }
  .Object@calls$load_matrix <- matrix_fun
  
  
  
  ## expression load call
  expression_fun <- function(
    path = expr_dir,
    gene_column_index = c_index,
    remove_zero_rows = TRUE,
    split_by_type = TRUE,
    verbose = NULL
  ) {
    
    .visiumHD_expression(
      path = path,
      gene_column_index = gene_column_index,
      remove_zero_rows = remove_zero_rows,
      split_by_type = split_by_type,
      verbose = verbose
    )
  }
  .Object@calls$load_expression <- expression_fun
  
  
  
  ## tissue position load call
  tissue_position_fun <- function(
    path = spatial_dir,
    verbose = NULL
  ) {
    .visiumHD_tissue_positions(
      path = path,
      verbose = verbose
    )
  }
  .Object@calls$load_tissue_position <- tissue_position_fun
  
  
  
  ## metadata load call
  meta_fun <- function(
    path = spatial_dir,
    verbose = NULL) {
    
    .visiumHD_meta(
      path = path,
      verbose = verbose
    )
  }
  .Object@calls$load_metadata <- meta_fun
  
  
  
  ## transcript load call
  transcript_fun <- function(expr_path = expr_dir,
                             tissue_positions_path = spatial_dir,
                             barcodes = .Object@barcodes,
                             array_subset_row = .Object@array_subset_row,
                             array_subset_col = .Object@array_subset_col,
                             pxl_subset_row = .Object@pxl_subset_row,
                             pxl_subset_col = .Object@pxl_subset_col) {
    
    .visiumHD_transcript(expr_path = expr_path,
                         tissue_positions_path = tissue_positions_path,
                         barcodes = barcodes,
                         array_subset_row = array_subset_row,
                         array_subset_col = array_subset_col,
                         pxl_subset_row = pxl_subset_row,
                         pxl_subset_col = pxl_subset_col,
                         verbose = TRUE)
    
  }
  .Object@calls$load_transcripts <- transcript_fun
  
  return(.Object)
})


# * access ####

#' @export
setMethod("$", signature("VisiumHDReader"), function(x, name) {
  basic_info <- c("visiumHD_dir", "expression_source", "gene_column_index", "barcodes",
                  "array_subset_row", "array_subset_col",
                  "pxl_subset_row", "pxl_subset_col")
  if (name %in% basic_info) return(methods::slot(x, name))
  
  return(x@calls[[name]])
})

#' @export
setMethod("$<-", signature("VisiumHDReader"), function(x, name, value) {
  basic_info <- c("visiumHD_dir", "expression_source", "gene_column_index", "barcodes",
                  "array_subset_row", "array_subset_col",
                  "pxl_subset_row", "pxl_subset_col")
  if (name %in% basic_info) {
    methods::slot(x, name) <- value
    return(initialize(x))
  }
  
  stop(sprintf("Only items in '%s' can be set",
               paste0(basic_info, collapse = "', '")))
})

#' @export
`.DollarNames.VisiumHDReader` <- function(x, pattern) {
  dn <- c("visiumHD_dir", "expression_source", "gene_column_index", "barcodes",
          "array_subset_row", "array_subset_col",
          "pxl_subset_row", "pxl_subset_col")
  if (length(methods::slot(x, "calls")) > 0) {
    dn <- c(dn, paste0(names(methods::slot(x, "calls")), "()"))
  }
  return(dn)
}





.visiumHD_matrix = function(path,
                            gene_column_index = 2,
                            remove_zero_rows = TRUE,
                            split_by_type = TRUE,
                            verbose = TRUE) {
  
  # check if path is provided
  if (missing(path)) {
    stop(wrap_txt(
      "No path to matrix file provided or auto-detected"
    ), call. = FALSE)
  }
  
  # check existence and access rights of files
  checkmate::assert_directory_exists(path)
  
  vmsg(.v = verbose, "loading expression matrix ...")
  vmsg(.v = verbose, .is_debug = TRUE, path)
  
  # load expression results with the 10X default matrix function
  matrix_results <- get10Xmatrix(path_to_data = path, 
                                 gene_column_index = gene_column_index,
                                 remove_zero_rows = remove_zero_rows, 
                                 split_by_type = split_by_type)
  
  return(matrix_results)
  
}





.visiumHD_expression = function(path,
                                gene_column_index = 2,
                                remove_zero_rows = TRUE,
                                split_by_type = TRUE,
                                verbose = TRUE) {
  
  # check if path is provided
  if (missing(path)) {
    stop(wrap_txt(
      "No path to matrix file provided or auto-detected"
    ), call. = FALSE)
  }
  
  # check existence and access rights of files
  checkmate::assert_directory_exists(path)
  
  vmsg(.v = verbose, "loading expression matrix ...")
  vmsg(.v = verbose, .is_debug = TRUE, path)
  
  # load expression results with the 10X default matrix function
  matrix_results <- get10Xmatrix(path_to_data = path, 
                                 gene_column_index = gene_column_index,
                                 remove_zero_rows = remove_zero_rows, 
                                 split_by_type = split_by_type)
  
  
  exprObj = createExprObj(expression_data = matrix_results,
                          spat_unit = "pixel",
                          feat_type = 'rna',
                          name = "raw",
                          provenance = "pixel")
  
  return(list('rna' = exprObj))
  
}




.visiumHD_tissue_positions = function(path,
                                      verbose = TRUE) {
  
  # check if path is provided
  if (missing(path)) {
    stop(wrap_txt(
      "No path to tissue positions file provided or auto-detected"
    ), call. = FALSE)
  }
  
  # check existence and access rights of files
  checkmate::assert_directory_exists(path)
  
  vmsg(.v = verbose, "loading tissue positions file ...")
  vmsg(.v = verbose, .is_debug = TRUE, path)
  
  # check existence and access rights of files
  tissue_positions_path = file.path(path, 'tissue_positions.parquet')
  checkmate::assert_file_exists(tissue_positions_path)
  
  # read with parquet and data.table
  tissue_positions = data.table::as.data.table(x = arrow::read_parquet(tissue_positions_path))
  
  return(tissue_positions)
  
}




.visiumHD_meta = function(
    path,
    verbose = TRUE) {
  
  # check if path is provided
  if (missing(path)) {
    stop(wrap_txt(
      "No path to tissue positions file provided or auto-detected"
    ), call. = FALSE)
  }
  
  # check existence and access rights of files
  checkmate::assert_directory_exists(path)
  
  vmsg(.v = verbose, "loading tissue positions file ...")
  vmsg(.v = verbose, .is_debug = TRUE, path)
  
  # check existence and access rights of files
  tissue_positions_path = file.path(path, 'tissue_positions.parquet')
  checkmate::assert_file_exists(tissue_positions_path)
  
  # read with parquet and data.table
  tissue_positions = data.table::as.data.table(x = arrow::read_parquet(tissue_positions_path))
  
  vmsg(.v = verbose, "creating metadata ...")
  
  data.table::setnames(tissue_positions, 'barcode', 'cell_ID')
  
  cx <- createCellMetaObj(
    metadata = tissue_positions,
    spat_unit = "pixel",
    feat_type = "rna",
    provenance = "pixel",
    verbose = verbose
  )
  return(cx)
  
}



.visiumHD_transcript = function(expr_path,
                                gene_column_index = 2,
                                remove_zero_rows = TRUE,
                                split_by_type = TRUE,
                                tissue_positions_path,
                                barcodes = NULL,
                                array_subset_row = NULL,
                                array_subset_col = NULL,
                                pxl_subset_row = NULL,
                                pxl_subset_col = NULL,
                                verbose = TRUE) {
  
  
  # function to create expression matrix
  matrix = .visiumHD_matrix(
    path = expr_path,
    gene_column_index = gene_column_index,
    remove_zero_rows = remove_zero_rows,
    split_by_type = split_by_type,
    verbose = verbose
  )
  
  
  # function to create tissue position data.table
  tissue_positions = .visiumHD_tissue_positions(
    path = tissue_positions_path,
    verbose = verbose
  )
  
  
  
  vmsg(.v = verbose, "creating visiumHD tissue position x expression data file ...")
  
  # subset data
  if(!is.null(barcodes)) {
    vmsg(.v = verbose, "subsetting visiumHD on barcodes")
    tissue_positions = tissue_positions[barcode %in% barcodes]
  }  
  
  if(!is.null(array_subset_row)) {
    if(is.vector(array_subset_row) & length(array_subset_row) == 2) {
      vmsg(.v = verbose, "subsetting visiumHD on array rows")
      tissue_positions = tissue_positions[array_row > array_subset_row[1] & array_row < array_subset_row[2]]
    } else {
      stop('array_subset_row was provided but is not a vector with length 2') 
    }
  }
  
  if(!is.null(array_subset_col)) {
    if(is.vector(array_subset_col) & length(array_subset_col) == 2) {
      vmsg(.v = verbose, "subsetting visiumHD on array columns")
      tissue_positions = tissue_positions[array_col > array_subset_col[1] & array_col < array_subset_col[2]]
    } else {
      stop('array_subset_col was provided but is not a vector with length 2') 
    }
  }
  
  if(!is.null(pxl_subset_row)) {
    if(is.vector(pxl_subset_row) & length(pxl_subset_row) == 2) {
      vmsg(.v = verbose, "subsetting visiumHD on row pixels")
      tissue_positions = tissue_positions[pxl_row_in_fullres > pxl_subset_row[1] & pxl_row_in_fullres < pxl_subset_row[2]]
    } else {
      cat('pxl_subset_row is ', pxl_subset_row)
      stop('pxl_subset_row was provided but is not a vector with length 2') 
    }
  }
  
  if(!is.null(pxl_subset_col)) {
    if(is.vector(pxl_subset_col) & length(pxl_subset_col) == 2) {
      vmsg(.v = verbose, "subsetting visiumHD on column pixels")
      tissue_positions = tissue_positions[pxl_col_in_fullres > pxl_subset_col[1] & pxl_col_in_fullres < pxl_subset_col[2]]
    } else {
      cat(pxl_subset_col)
      stop('pxl_subset_col was provided but is not a vector with length 2') 
    }
  }
  
  # also subset matrix if needed
  if(any(!is.null(c(barcodes, 
                    array_subset_row, array_subset_col, 
                    pxl_subset_row, pxl_subset_col)))) {
    vmsg(.v = verbose, "subsetting visiumHD on expression matrix")
    matrix = matrix[, colnames(matrix) %in% tissue_positions$barcode]
  }
  
  
  
  
  
  
  # convert expression matrix to minimal data.table object
  matrix_tile_dt = data.table::as.data.table(Matrix::summary(matrix))
  genes = matrix@Dimnames[[1]]
  samples = matrix@Dimnames[[2]]
  matrix_tile_dt[, gene := genes[i]]
  matrix_tile_dt[, pixel := samples[j]]
  
  
  # merge data.table matrix and spatial coordinates to create input for Giotto Polygons
  gpoints = data.table::merge.data.table(matrix_tile_dt, tissue_positions, by.x = 'pixel', by.y = 'barcode')
  gpoints = gpoints[,.(pixel, pxl_row_in_fullres, pxl_col_in_fullres, gene, x)]
  colnames(gpoints) = c('pixel', 'x', 'y', 'gene', 'counts')
  
  gpoints = createGiottoPoints(x = gpoints[,.(x, y, gene, pixel, counts)])
  
  # ensure output is always a list
  if (!is.list(gpoints)) {
    gpoints <- list(gpoints)
    names(gpoints) <- objName(gpoints[[1L]])
  }
  
  return(list('matrix' = matrix, 'tissue_positions' = tissue_positions, 'gpoints' = gpoints))
  
}

