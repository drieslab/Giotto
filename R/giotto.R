
#### Giotto class ####

#' @title S4 giotto Class
#' @description Framework of giotto object to store and work with spatial expression data
#' @keywords giotto, object
#' @slot expression expression information
#' @slot expression_feat available features (e.g. rna, protein, ...)
#' @slot spatial_locs spatial location coordinates for cells/spots/grids
#' @slot spatial_info information about spatial units (Giotto spatVector)
#' @slot cell_metadata metadata for cells
#' @slot feat_metadata metadata for available features
#' @slot feat_info information about features (Giotto spatVector)
#' @slot cell_ID unique cell IDs
#' @slot feat_ID unique feature IDs for all features or modalities
#' @slot spatial_network spatial network in data.table/data.frame format
#' @slot spatial_grid spatial grid in data.table/data.frame format
#' @slot spatial_enrichment slot to save spatial enrichment-like results
#' @slot dimension_reduction slot to save dimension reduction coordinates
#' @slot nn_network nearest neighbor network in igraph format
#' @slot images slot to store giotto images
#' @slot parameters slot to save parameters that have been used
#' @slot instructions slot for global function instructions
#' @slot offset_file offset file used to stitch together image fields
#' @slot OS_platform Operating System to run Giotto analysis on
#' @slot join_info information about joined Giotto objects
#' @details
#' [\strong{expression}] There are several ways to provide expression information:
#'
#' [\strong{expression_feat}] The different featurs or modalities such as rna, protein, metabolites, ...
#' that are provided in the expression slot.
#'
#'
#' @export
giotto <- setClass(
  "giotto",
  slots = c(
    expression = "ANY",
    expression_feat = "ANY",
    spatial_locs = "ANY",
    spatial_info = "ANY",
    cell_metadata = "ANY",
    feat_metadata = "ANY",
    feat_info = "ANY",
    cell_ID = "ANY",
    feat_ID = "ANY",
    spatial_network = "ANY",
    spatial_grid = "ANY",
    spatial_enrichment = "ANY",
    dimension_reduction = 'ANY',
    nn_network = "ANY",
    images = "ANY",
    parameters = "ANY",
    instructions = "ANY",
    offset_file = "ANY",
    OS_platform = "ANY",
    join_info = "ANY"

  ),

  prototype = list(
    expression = NULL,
    expression_feat = NULL,
    spatial_locs = NULL,
    spatial_info = NULL,
    cell_metadata = NULL,
    feat_metadata = NULL,
    feat_info = NULL,
    cell_ID = NULL,
    feat_ID = NULL,
    spatial_network = NULL,
    spatial_grid = NULL,
    spatial_enrichment = NULL,
    dimension_reduction = NULL,
    nn_network = NULL,
    images = NULL,
    parameters = NULL,
    instructions = NULL,
    offset_file = NULL,
    OS_platform = NULL,
    join_info = NULL
  )
)



#' show method for giotto class
#' @param object giotto object
#' @aliases show,giotto-method
#' @docType methods
#' @rdname show-methods

setMethod(
  f = "show",
  signature = "giotto",
  definition = function(object) {

    cat("An object of class",  class(object), "\n")

    for(feat_type in unique(object@expression_feat)) {
      cat("features = ", feat_type, "\n")

      cat(
        nrow(x = object@expression[[feat_type]][['raw']]),
        "features across",
        ncol(x = object@expression[[feat_type]][['raw']]),
        "samples.\n \n"
      )
    }

    cat('Steps and parameters used: \n \n')
    print(object@parameters)
    invisible(x = NULL)
  }
)


#### Giotto python path ####

#' @title set_giotto_python_path
#' @name set_giotto_python_path
#' @description sets the python path and/or installs miniconda and the python modules
set_giotto_python_path = function(python_path = NULL,
                                  packages_to_install = c('pandas', 'networkx', 'python-igraph',
                                                          'leidenalg', 'python-louvain', 'python.app',
                                                          'scikit-learn')) {

  ## if a python path is provided, use that path
  if(!is.null(python_path)) {
    cat('\n python path provided \n')
    python_path = as.character(python_path)
    reticulate::use_python(required = T, python = python_path)

  } else {

    ## identify operating system and adjust the necessary packages
    os_specific_system = get_os()

    # exclude python.app for windows and linux
    if(os_specific_system != 'osx') {
      packages_to_install = packages_to_install[packages_to_install != 'python.app']
    }



    ## check if giotto environment is already installed
    conda_path = reticulate::miniconda_path()
    if(os_specific_system == 'osx') {
      full_path = paste0(conda_path, "/envs/giotto_env/bin/pythonw")
    } else if(os_specific_system == 'windows') {
      full_path = paste0(conda_path, "\\envs\\giotto_env\\python.exe")
    } else if(os_specific_system == 'linux') {
      full_path = paste0(conda_path, "/envs/giotto_env/bin/python")
    }

    if(file.exists(full_path)) {
      cat('\n giotto environment found \n')
      python_path = full_path
      reticulate::use_python(required = T, python = python_path)


    } else {

      ## if giotto environment is not found, ask to install miniconda and giotto environment
      install_giotto = utils::askYesNo(msg = 'Install a miniconda Python environment for Giotto?')

      if(install_giotto == TRUE) {

        ## check and install miniconda locally if necessary
        conda_path = reticulate::miniconda_path()

        if(!file.exists(conda_path)) {
          cat('\n |---- install local miniconda ----| \n')
          reticulate::install_miniconda()
        }

        cat('\n |---- install giotto environment ----| \n')
        conda_path = reticulate::miniconda_path()

        ## for unix-like systems ##
        if(.Platform[['OS.type']] == 'unix') {

          conda_full_path = paste0(conda_path,'/','bin/conda')
          reticulate::conda_create(envname = 'giotto_env',
                                   conda = conda_full_path)


          full_envname = paste0(conda_path,'/envs/giotto_env')

          if(os_specific_system == 'osx') {
            python_full_path = paste0(conda_path, "/envs/giotto_env/bin/pythonw")
          } else if(os_specific_system == 'linux') {
            python_full_path = paste0(conda_path, "/envs/giotto_env/bin/python")
          }


          reticulate::py_install(packages = packages_to_install,
                                 envname = 'giotto_env',
                                 method = 'conda',
                                 conda = conda_full_path,
                                 python_version = '3.6')

          reticulate::py_install(packages = 'smfishhmrf',
                                 envname = full_envname,
                                 method = 'conda',
                                 conda = conda_full_path,
                                 pip = TRUE,
                                 python_version = '3.6')


          ## for windows systems ##
        } else if(.Platform[['OS.type']] == 'windows') {

          conda_full_path = paste0(conda_path,'/','condabin/conda.bat')
          reticulate::conda_create(envname = 'giotto_env',
                                   conda = conda_full_path)


          full_envname = paste0(conda_path,'/envs/giotto_env')
          python_full_path = paste0(conda_path, "/envs/giotto_env/python.exe")

          reticulate::py_install(packages = packages_to_install,
                                 envname = 'giotto_env',
                                 method = 'conda',
                                 conda = conda_full_path,
                                 python_version = '3.6',
                                 channel = c('conda-forge', 'vtraag'))

          reticulate::py_install(packages = 'smfishhmrf',
                                 envname = full_envname,
                                 method = 'conda',
                                 conda = conda_full_path,
                                 pip = TRUE,
                                 python_version = '3.6')
        }


        python_path = python_full_path
        reticulate::use_python(required = T, python = python_path)

      } else {

        if(.Platform[['OS.type']] == 'unix') {
          python_path = try(system('which python', intern = T))
        } else if(.Platform[['OS.type']] == 'windows') {
          python_path = try(system('where python', intern = T))
        }

        if(class(python_path) == "try-error") {
          cat('\n no python path found, set it manually when needed \n')
          python_path = '/need/to/set/path/to/python'
        } else {
          cat('\n try default python \n')
          python_path = python_path
          reticulate::use_python(required = T, python = python_path)
        }
      }
    }
  }
  return(python_path)
}





#### Giotto instructions ####


#' @title createGiottoInstructions
#' @description Function to set global instructions for giotto functions
#' @param python_path path to python binary to use
#' @param show_plot print plot to console, default = TRUE
#' @param return_plot return plot as object, default = TRUE
#' @param save_plot automatically save plot, dafault = FALSE
#' @param save_dir path to directory where to save plots
#' @param plot_format format of plots (defaults to png)
#' @param dpi resolution for raster images
#' @param units units of format (defaults to in)
#' @param height height of plots
#' @param width width of  plots
#' @param is_docker using docker implementation of Giotto (defaults to FALSE)
#' @param plot_count [global option] start count for creating automatic unique plots
#' @return named vector with giotto instructions
#' @seealso More online information can be found here \url{https://rubd.github.io/Giotto_site/articles/instructions_and_plotting.html}
#' @export
createGiottoInstructions <- function(python_path =  NULL,
                                     show_plot = NULL,
                                     return_plot = NULL,
                                     save_plot = NULL,
                                     save_dir = NULL,
                                     plot_format = NULL,
                                     dpi = NULL,
                                     units = NULL,
                                     height = NULL,
                                     width = NULL,
                                     is_docker = FALSE,
                                     plot_count = 0,
                                     fiji_path = NULL) {

  # pyton path to use
  if(is_docker){
    python_path = set_giotto_python_path(python_path = "/usr/bin/python3") # fixed in docker version
  }  else{
    python_path = set_giotto_python_path(python_path = python_path)
  }

  # print plot to console
  if(is.null(show_plot)) {
    show_plot = TRUE
  }

  # print plot to console
  if(is.null(return_plot)) {
    return_plot = TRUE
  }

  # print plot to console
  if(is.null(save_plot)) {
    save_plot = FALSE
  }

  # directory to save results to
  if(is.null(save_dir)) {
    save_dir = getwd()
  }
  save_dir = as.character(save_dir)

  # plot format
  if(is.null(plot_format)) {
    plot_format = "png"
  }
  plot_format = as.character(plot_format)

  # dpi of raster images
  if(is.null(dpi)) {
    dpi = 300
  }
  dpi = as.numeric(dpi)

  # units for height and width
  if(is.null(units)) {
    units = 'in'
  }
  units = as.character(units)

  # height of plot
  if(is.null(height)) {
    height = 9
  }
  height = as.numeric(height)

  # width of plot
  if(is.null(width)) {
    width = 9
  }
  width = as.numeric(width)


  ## global options ##
  # ---------------- #

  # plot count
  options('giotto.plot_count' = plot_count)

  # fiji path
  options('giotto.fiji' = fiji_path)


  # return instructions list
  instructions_list = list(python_path, show_plot, return_plot,
                           save_plot, save_dir, plot_format, dpi,
                           units, height, width, is_docker)
  names(instructions_list) = c('python_path', 'show_plot', 'return_plot',
                               'save_plot', 'save_dir', 'plot_format', 'dpi',
                               'units', 'height', 'width', 'is_docker')
  return(instructions_list)

}


#' @title readGiottoInstrunctions
#' @description Retrieves the instruction associated with the provided parameter
#' @param giotto_instructions giotto object or result from createGiottoInstructions()
#' @param param parameter to retrieve
#' @return specific parameter
#' @export
readGiottoInstructions <- function(giotto_instructions,
                                   param = NULL) {

  # get instructions if provided the giotto object
  if(class(giotto_instructions) == 'giotto') {
    giotto_instructions = giotto_instructions@instructions
  }

  # stop if parameter is not found
  if(is.null(param)) {
    stop('\t readGiottoInstructions needs a parameter to work \t')
  } else if(!param %in% names(giotto_instructions)) {
    stop('\t parameter ', param, ' is not part of Giotto parameters \t')
  } else {
    specific_instruction = giotto_instructions[[param]]
  }
  return(specific_instruction)
}


#' @title showGiottoInstructions
#' @description Function to display all instructions from giotto object
#' @param gobject giotto object
#' @return named vector with giotto instructions
#' @export
showGiottoInstructions = function(gobject) {

  instrs = gobject@instructions
  return(instrs)
}


#' @title changeGiottoInstructions
#' @description Function to change one or more instructions from giotto object
#' @param gobject giotto object
#' @param params parameter(s) to change
#' @param new_values new value(s) for parameter(s)
#' @param return_gobject (boolean) return giotto object
#' @return giotto object with one or more changed instructions
#' @export
changeGiottoInstructions = function(gobject,
                                    params = NULL,
                                    new_values = NULL,
                                    return_gobject = TRUE) {

  instrs = gobject@instructions

  if(is.null(params) | is.null(new_values)) {
    stop('\t params and new_values can not be NULL \t')
  }

  if(length(params) != length(new_values)) {
    stop('\t length of params need to be the same as new values \t')
  }

  if(!all(params %in% names(instrs))) {
    stop('\t all params need to be part of Giotto instructions \t')
  }

  ## swap with new values
  instrs[params] = new_values

  ## make sure that classes remain consistent
  new_instrs = lapply(1:length(instrs), function(x) {

    if(names(instrs[x]) %in% c('dpi', 'height', 'width')) {
      instrs[[x]] = as.numeric(instrs[[x]])
    } else {
      instrs[[x]] = instrs[[x]]
    }

  })

  names(new_instrs) = names(instrs)



  if(return_gobject == TRUE) {
    gobject@instructions = new_instrs
    return(gobject)
  } else {
    return(new_instrs)
  }

}




#' @title replaceGiottoInstructions
#' @description Function to replace all instructions from giotto object
#' @param gobject giotto object
#' @param instructions new instructions (e.g. result from createGiottoInstructions)
#' @return giotto object with replaces instructions
#' @export
replaceGiottoInstructions = function(gobject,
                                     instructions = NULL) {

  old_instrs_list = gobject@instructions

  # validate new instructions
  if(!all(names(instructions) %in% names(old_instrs_list)) | is.null(instructions)) {
    stop('\n You need to provide a named list for all instructions, like the outcome of createGiottoInstructions \n')
  } else {
    gobject@instructions = instructions
    return(gobject)
  }

}



#### Giotto matrices ####

#' @title readExprMatrix
#' @description Function to read an expression matrix into a sparse matrix.
#' @param path path to the expression matrix
#' @param cores number of cores to use
#' @param transpose transpose matrix
#' @return sparse matrix
#' @details The expression matrix needs to have both unique column names and row names
#' @export
readExprMatrix = function(path,
                          cores = NA,
                          transpose = FALSE) {

  # check if path is a character vector and exists
  if(!is.character(path)) stop('path needs to be character vector')
  if(!file.exists(path)) stop('the path: ', path, ' does not exist')

  # set number of cores automatically, but with limit of 10
  cores = determine_cores(cores)
  data.table::setDTthreads(threads = cores)

  # read and convert
  DT = suppressWarnings(data.table::fread(input = path, nThread = cores))
  spM = Matrix::Matrix(as.matrix(DT[,-1]), dimnames = list(DT[[1]], colnames(DT[,-1])), sparse = T)

  if(transpose == TRUE) {
    spM = t_flex(spM)
  }

  return(spM)
}



#' @title evaluate_expr_matrix
#' @description Evaluate expression matrices.
#' @param inputmatrix inputmatrix to evaluate
#' @param sparse create sparse matrix (default = TRUE)
#' @param cores how many cores to use
#' @return sparse matrix
#' @details The inputmatrix can be a matrix, sparse matrix, data.frame, data.table or path to any of these.
#' @keywords internal
evaluate_expr_matrix = function(inputmatrix,
                                sparse = TRUE,
                                cores = NA) {


  if(methods::is(inputmatrix, 'character')) {
    mymatrix = readExprMatrix(inputmatrix, cores =  cores)
  } else if(methods::is(inputmatrix, 'Matrix')) {
    mymatrix = inputmatrix
  } else if(methods::is(inputmatrix, 'DelayedMatrix')) {
    mymatrix = inputmatrix
  } else if(methods::is(inputmatrix, 'data.table')) {
    if(sparse == TRUE) {
      # force sparse class
      mymatrix = Matrix::Matrix(as.matrix(inputmatrix[,-1]),
                                dimnames = list(inputmatrix[[1]],
                                                colnames(inputmatrix[,-1])), sparse = TRUE)
    } else {
      # let Matrix decide
      mymatrix = Matrix::Matrix(as.matrix(inputmatrix[,-1]),
                                dimnames = list(inputmatrix[[1]],
                                                colnames(inputmatrix[,-1])))
    }

  } else if(inherits(inputmatrix, what = c('data.frame', 'matrix'))) {
    mymatrix = methods::as(as.matrix(inputmatrix), "sparseMatrix")
  } else {
    stop("raw_exprs needs to be a path or an object of class 'Matrix', 'data.table', 'data.frame' or 'matrix'")
  }

  return(mymatrix)
}



#' @name depth
#' @keywords internal
depth <- function(this) {
  if(is.list(this) && length(this) == 0) {
    return(0)
  }
  ifelse(is.list(this), 1L + max(sapply(this, depth)), 0L)
}


#' @name extr_single_list
#' @keywords internal
extr_single_list = function(gobject,
                            expr_s_list,
                            expression_feat = 'rna',
                            cores = 1) {



  # to make it compatible with previous version
  #if(length(expr_s_list) == 1 & !list(expr_s_list)) {
  #  expr_s_list = list('raw' = expr_s_list)
  #}

  expression_name = names(expr_s_list)

  ## 1. raw expression data is required
  if(!'raw' %in% expression_name) {
    stop("\n the raw expression matrix with name 'raw' must be provided \n")
  } else {

    exprname = 'raw'
    expr_i = match('raw', expression_name)
    expr     = expr_s_list[[expr_i]]
    expr_res = evaluate_expr_matrix(inputmatrix = expr, cores = cores, sparse = TRUE)

    # check rownames and colnames
    if(any(duplicated(rownames(expr_res)))) {
      stop("row names for ", exprname, " contains duplicates, please remove or rename")
    }

    if(any(duplicated(colnames(expr_res)))) {
      stop("column names for ", exprname, " contains duplicates, please remove or rename")
    }


    # update giotto object
    gobject@expression[[expression_feat]] = list()
    gobject@expression[[expression_feat]][['raw']] = expr_res

    # prepare other slots related to raw matrix

    # check existin cell_IDs; need to be consistent across features
    if(!is.null(gobject@cell_ID)) {
      test = identical(gobject@cell_ID, colnames(expr_res))
      if(test == FALSE) stop('column names (cell ids) between feature expression matrices need to be identical and in the same order')
    } else {
      gobject@cell_ID = colnames(expr_res)
    }

    # can be unique for each feature
    gobject@feat_ID[[expression_feat]] = rownames(expr_res)

    raw_exprs_dim = dim(expr_res)
  }

  ## 2. load all other provided matrices
  if(length(expression_name) > 1) {
    expression_name_other = expression_name[expression_name != 'raw']

    for(exprname in expression_name_other) {

      expr_i   = match(exprname, expression_name)
      expr     = expr_s_list[[expr_i]]
      expr_res = evaluate_expr_matrix(inputmatrix = expr, cores = cores, sparse = TRUE)

      # check rownames and colnames
      if(any(duplicated(rownames(expr_res)))) {
        stop("row names for ", exprname, " contains duplicates, please remove or rename")
      }

      if(any(duplicated(colnames(expr_res)))) {
        stop("column names for ", exprname, " contains duplicates, please remove or rename")
      }

      # compare processed matrix dimensions with raw matrix
      if(all(dim(expr_res) == raw_exprs_dim) &
         all(colnames(expr_res) == gobject@cell_ID) &
         all(rownames(expr_res) == gobject@feat_ID[[expression_feat]])) {

        gobject@expression[[expression_feat]][[exprname]] = expr_res

      } else {
        stop('\n dimensions, row or column names are not the same between ', exprname, ' and raw expression \n')
      }

    }

  }

  return(gobject)

}


#' @name extract_expression_list
#' @keywords internal
extract_expression_list = function(gobject,
                                   expr_list,
                                   expression_feat = 'rna',
                                   cores = 1,
                                   verbose = TRUE) {


  ## to make it compatible with previous version

  # single matrix
  if(inherits(expr_list, c('matrix', 'Matrix', 'DelayedMatrix'))) {
    expr_list = list('raw' = expr_list)
  }

  # single path to matrix
  if(length(expr_list) == 1 & !is.list(expr_list)) {
    expr_list = list('raw' = expr_list)
  }

  # 1. get depth of list
  if(verbose == TRUE) print(str(expr_list))
  list_depth = depth(expr_list)

  # nothing provided
  if(list_depth == 0) {
    stop('Depth of expression list is 0, no expression information is provided \n')
  }

  # only one set provided
  if(list_depth == 1) {
    if(verbose == TRUE) message('Depth of provided expression list is 1, working with one type of molecular feature \n')

    if(length(expression_feat) > list_depth) {
      stop('More expression feature types provided than expected \n')
    }

    gobject = extr_single_list(gobject = gobject,
                               expr_s_list = expr_list,
                               expression_feat = expression_feat,
                               cores = cores)



  } else {
    if(verbose == TRUE) message('Depth of provided expression list is more than 1, working with multiple types of molecular features \n')


    if(length(expression_feat) > length(expr_list)) {
      stop('More expression feature types provided than expected \n')
    }

    if(length(expression_feat) < length(expr_list)) {
      stop('Too few expression feature types provided than expected \n')
    }

    for(feat_type in names(expr_list)) {

      gobject = extr_single_list(gobject = gobject,
                                 expr_s_list = expr_list[[feat_type]],
                                 expression_feat = feat_type,
                                 cores = cores)

    }

  }

  return(gobject)

}




#### Giotto locations ####

#' @name evaluate_spatial_locations
#' @description Evaluate spatial location input
#' @param spatial_locs spatial locations to evaluate
#' @param cores how many cores to use
#' @param dummy_n number of rows to create dummy spaial locations
#' @param expr_matrix expression matrix to compare the cell IDs with
#' @return data.table
#' @keywords internal
evaluate_spatial_locations = function(spatial_locs,
                                      cores = 1,
                                      dummy_n = 2,
                                      expr_matrix = NULL) {

  if(is.null(spatial_locs)) {
    warning('\n spatial locations are not given, dummy 2D data will be created \n')

    # create 2D rectangular dummy positions
    ceil_value = ceiling(sqrt(dummy_n))
    dummy_matrix = t(utils::combn(x = ceil_value, m = 2))
    final_dummy = rbind(matrix(data = rep(1:ceil_value, 2), ncol = 2),
                        dummy_matrix,
                        dummy_matrix[, c(2,1)])
    final_dummy = final_dummy[1:dummy_n,]

    spatial_locs = data.table::data.table(sdimx = final_dummy[,1],
                                          sdimy = final_dummy[,2])

  } else {

    if(!any(class(spatial_locs) %in% c('data.table', 'data.frame', 'matrix', 'character'))) {
      stop('spatial_locs needs to be a data.table or data.frame-like object or a path to one of these')
    }
    if(methods::is(spatial_locs, 'character')) {
      if(!file.exists(spatial_locs)) stop('path to spatial locations does not exist')
      spatial_locs = data.table::fread(input = spatial_locs, nThread = cores)
    } else {
      spatial_locs = data.table::as.data.table(spatial_locs)
    }


    # check if all columns are numeric
    column_classes = lapply(spatial_locs, FUN = class)
    #non_numeric_classes = column_classes[column_classes != 'numeric']
    non_numeric_classes = column_classes[!column_classes %in% c('numeric','integer')]

    if(length(non_numeric_classes) > 0) {

      non_numeric_indices = which(!column_classes %in% c('numeric','integer'))

      warning('There are non numeric or integer columns for the spatial location input at column position(s): ', non_numeric_indices,
              '\n The first non-numeric column will be considered as a cell ID to test for consistency with the expression matrix',
              '\n Other non numeric columns will be removed')

      if(!is.null(expr_matrix)) {
        potential_cell_IDs = spatial_locs[[ non_numeric_indices[1] ]]
        expr_matrix_IDs = colnames(expr_matrix)

        if(!identical(potential_cell_IDs, expr_matrix_IDs)) {
          warning('The cell IDs from the expression matrix and spatial locations do not seem to be identical')
        }

      }


      spatial_locs = spatial_locs[,-non_numeric_indices, with = F]

    }

    # check number of columns: too few
    if(ncol(spatial_locs) < 2) {
      stop('There need to be at least 2 numeric columns for spatial locations \n')
    }

    # check number of columns: too many
    if(ncol(spatial_locs) > 3) {
      warning('There are more than 3 columns for spatial locations, only the first 3 will be used \n')
      spatial_locs = spatial_locs[, 1:3]
    }

  }


  return(spatial_locs)
}



#' @name extract_spatial_locations_list
#' @description Extract spatial locations
#' @param gobject giotto object
#' @param spat_loc_list list of spatial locations
#' @param cores how many cores to use
#' @param dummy_n number of rows to create dummy spaial locations
#' @param expr_matrix expression matrix to compare the cell IDs with
#' @param verbose be verbose
#' @return updated giotto object
#' @keywords internal
extract_spatial_locations_list = function(gobject,
                                          spat_loc_list,
                                          cores = 1,
                                          dummy_n = 1,
                                          expr_matrix = NULL,
                                          verbose = TRUE) {


  ## to make it compatible with previous version

  # single data.table-like
  if(inherits(spat_loc_list, c('data.frame', 'data.table', 'matrix', 'character'))) {
    spat_loc_list = list('raw' = spat_loc_list)
  }


  spat_locs_names = names(spat_loc_list)

  ## 1. raw expression data is required
  if(!'raw' %in% spat_locs_names) {
    stop("\n the raw / real / physical coordinates must be provided with name 'raw' \n")
  }


  gobject@spatial_locs = list()

  for(spatlocname in spat_locs_names) {

    cat('spatlocname: ', spatlocname, '\n')

    spatial_locs = spat_loc_list[[spatlocname]]

    spatial_locs = evaluate_spatial_locations(spatial_locs = spatial_locs,
                                              cores = cores,
                                              dummy_n = dummy_n,
                                              expr_matrix = expr_matrix)

    # check if dimensions agree
    # dummy_n = number of columns in expression matrix, also dummy number if no spatial location are being provided
    if(nrow(spatial_locs) != dummy_n) {
      stop('\n Number of rows of spatial location must equal number of columns of raw expression matrix \n')
    }

    # force dimension names
    spatial_dimensions = c('x', 'y', 'z')
    colnames(spatial_locs) = paste0('sdim', spatial_dimensions[1:ncol(spatial_locs)])

    # add cell_ID column
    spatial_locs[, cell_ID := gobject@cell_ID]

    gobject@spatial_locs[[spatlocname]] = spatial_locs

  }

  return(gobject)

}





#### Giotto spatial info ####


#' @name evaluate_spatial_info
#' @description Evaluate spatial information input
#' @param spatial_info spatial information to evaluate
#' @param cores how many cores to use
#' @param spatial_locs spatial location data.table to compare the cell IDs with
#' @return data.table
#' @keywords internal
evaluate_spatial_info = function(spatial_info,
                                 cores = 1,
                                 spatial_locs) {


  ## 1. load or read spatial information data ##
  if(inherits(spatial_info, 'character')) {

    if(!file.exists(spatial_info)) stop('path to spatial information does not exist')
    spatial_info = data.table::fread(input = spatial_info, nThread = cores)

  } else if(inherits(spatial_info, 'data.frame')) {

    spatial_info = data.table::as.data.table(spatial_info)

  } else {

    stop('If spatial information is provided then it needs to be a file path or a data.frame-like object')

  }

  ## 2. check and name columns ##
  nr_cols = ncol(spatial_info)

  if(nr_cols < 4) stop('Spatial information needs to have at least 4 columns: \n',
                       'x, y, (z) information columns \n',
                       'cell ID and polygon point column \n')

  column_classes = lapply(spatial_info, FUN = class)

  # 3D or 2D data
  if(all(column_classes[1:3] == 'numeric')) {
    colnames(spatial_info)[1:5] = c('sdimx', 'sdimy', 'sdimz', 'cell_ID', 'point')
  } else if(all(column_classes[1:2] == 'numeric')){
    colnames(spatial_info)[1:4] = c('sdimx', 'sdimy', 'cell_ID', 'point')
  } else {
    stop('First 3 or 2 columns need to be numeric for 3D and 2D data respectively')
  }



  ## 3. check cell ID ##
  spatial_locs_cell_IDs = spatial_locs[['cell_ID']]

  spatial_info_cell_IDs = spatial_info[['cell_ID']]

  if(all(spatial_info_cell_IDs %in% spatial_locs_cell_IDs)) {
    return(spatial_info)
  } else {
    stop('cell IDs in spatial information are missing in the spatial locations slot')
  }

}




#### Giotto spatial feature info ####


#' @name evaluate_feat_info
#' @description Evaluate spatial feature information input
#' @param spatial_feat_info spatial feature information to evaluate
#' @param cores how many cores to use
#' @param feat_ID feature IDs to check with
#' @return data.table
#' @keywords internal
evaluate_feat_info = function(spatial_feat_info,
                              feat_type,
                              cores = 1,
                              feat_ID) {


  ## 1. load or read spatial information data ##
  if(inherits(spatial_feat_info, 'character')) {

    if(!file.exists(spatial_feat_info)) stop('path to spatial information does not exist')
    spatial_feat_info = data.table::fread(input = spatial_feat_info, nThread = cores)

  } else if(inherits(spatial_feat_info, 'data.frame')) {

    spatial_feat_info = data.table::as.data.table(spatial_feat_info)

  } else {

    stop('If spatial feature information is provided then it needs to be a file path or a data.frame-like object')

  }


  ## 2. check and name columns ##
  nr_cols = ncol(spatial_feat_info)

  if(nr_cols < 3) stop('Spatial feature information needs to have at least 3 columns: \n',
                       'x, y, (z) information columns \n',
                       'and feature ID column \n')

  column_classes = lapply(spatial_feat_info, FUN = class)


  # 3D or 2D data
  if(all(column_classes[1:3] == 'numeric')) {
    colnames(spatial_feat_info)[1:4] = c('sdimx', 'sdimy', 'sdimz', 'feat_ID')
  } else if(all(column_classes[1:2] == 'numeric')){
    colnames(spatial_feat_info)[1:3] = c('sdimx', 'sdimy', 'feat_ID')
  } else {
    stop('First 3 or 2 columns need to be numeric for 3D and 2D data respectively')
  }



  ## 3. check cell ID ##

  spatial_feature_info_feat_IDs = spatial_feat_info[['feat_ID']]

  if(all(spatial_feature_info_feat_IDs %in% feat_ID)) {
    return(spatial_feat_info)
  } else {
    stop('feat IDs in spatial feature information are missing in the feature ID slot')
  }

}


#### creating Giotto objects ####

#' @title create Giotto object
#' @description Function to create a giotto object
#' @param expression expression information
#' @param raw_exprs deprecated, use expression
#' @param expression_feat available features (e.g. rna, protein, ...)
#' @param spatial_locs data.table or data.frame with coordinates for cell centroids
#' @param spatial_info information about spatial units
#' @param spatial_info list of giotto polygon objects with spatial information,
#' see \code{\link{createGiottoPolygonsFromMask}} and \code{\link{createGiottoPolygonsFromDfr}}
#' @param cell_metadata cell annotation metadata
#' @param feat_metadata feature annotation metadata for each unique feature
#' @param feat_info list of giotto point objects with feature info,
#' see \code{\link{createGiottoPoints}}
#' @param spatial_network list of spatial network(s)
#' @param spatial_network_name list of spatial network name(s)
#' @param spatial_grid list of spatial grid(s)
#' @param spatial_grid_name list of spatial grid name(s)
#' @param spatial_enrichment list of spatial enrichment score(s) for each spatial region
#' @param spatial_enrichment_name list of spatial enrichment name(s)
#' @param dimension_reduction list of dimension reduction(s)
#' @param nn_network list of nearest neighbor network(s)
#' @param images list of images
#' @param offset_file file used to stitch fields together (optional)
#' @param instructions list of instructions or output result from \code{\link{createGiottoInstructions}}
#' @param cores how many cores or threads to use to read data if paths are provided
#' @param verbose be verbose when building Giotto object
#' @return giotto object
#' @details
#'
#' See \url{https://rubd.github.io/Giotto_site/articles/howto_giotto_class.html} for more details
#'
#' [\strong{Requirements}] To create a giotto object you need to provide at least a matrix with genes as
#' row names and cells as column names. This matrix can be provided as a base matrix, sparse Matrix, data.frame,
#' data.table or as a path to any of those.
#' To include spatial information about cells (or regions) you need to provide a matrix, data.table or data.frame (or path to them)
#' with coordinates for all spatial dimensions. This can be 2D (x and y) or 3D (x, y, x).
#' The row order for the cell coordinates should be the same as the column order for the provided expression data.
#'
#' [\strong{Instructions}] Additionally an instruction file, generated manually or with \code{\link{createGiottoInstructions}}
#' can be provided to instructions, if not a default instruction file will be created
#' for the Giotto object.
#'
#' [\strong{Multiple fields}] In case a dataset consists of multiple fields, like seqFISH+ for example,
#' an offset file can be provided to stitch the different fields together. \code{\link{stitchFieldCoordinates}}
#' can be used to generate such an offset file.
#'
#' [\strong{Processed data}] Processed count data, such as normalized data, can be provided using
#' one of the different expression slots (norm_expr, norm_scaled_expr, custom_expr).
#'
#' [\strong{Metadata}] Cell and gene metadata can be provided using the cell and gene metadata slots.
#' This data can also be added afterwards using the \code{\link{addGeneMetadata}} or \code{\link{addCellMetadata}} functions.
#'
#' [\strong{Other information}] Additional information can be provided through the appropriate slots:
#' \itemize{
#'   \item{spatial networks}
#'   \item{spatial girds}
#'   \item{spatial enrichments}
#'   \item{dimensions reduction}
#'   \item{nearest neighbours networks}
#'   \item{images}
#' }
#'
#' @keywords giotto
#' @importFrom methods new
#' @export
createGiottoObject <- function(expression,
                               raw_exprs = NULL,
                               expression_feat = 'rna',
                               spatial_locs = NULL,
                               spatial_info = NULL,
                               cell_metadata = NULL,
                               feat_metadata = NULL,
                               feat_info = NULL,
                               spatial_network = NULL,
                               spatial_network_name = NULL,
                               spatial_grid = NULL,
                               spatial_grid_name = NULL,
                               spatial_enrichment = NULL,
                               spatial_enrichment_name = NULL,
                               dimension_reduction = NULL,
                               nn_network = NULL,
                               images = NULL,
                               offset_file = NULL,
                               instructions = NULL,
                               cores = NA,
                               verbose = TRUE) {

  # create minimum giotto
  gobject = giotto(expression = list(),
                   expression_feat = expression_feat,
                   spatial_locs = spatial_locs,
                   spatial_info = NULL,
                   cell_metadata = cell_metadata,
                   feat_metadata = feat_metadata,
                   feat_info = feat_info,
                   cell_ID = NULL,
                   feat_ID = NULL,
                   spatial_network = NULL,
                   spatial_grid = NULL,
                   spatial_enrichment = NULL,
                   dimension_reduction = NULL,
                   nn_network = NULL,
                   images = NULL,
                   parameters = NULL,
                   offset_file = offset_file,
                   instructions = instructions,
                   OS_platform = .Platform[['OS.type']],
                   join_info = NULL)


  ## data.table: set global variable
  cell_ID = feat_ID = NULL

  ## check if all optional packages are installed
  # TODO: update at the end
  # TODO: extract from suggest field of DESCRIPTION
  extra_packages = c("scran", "MAST", "png", "tiff", "biomaRt", "trendsceek", "multinet", "RTriangle", "FactoMiner")

  pack_index = extra_packages %in% rownames(utils::installed.packages())
  extra_installed_packages = extra_packages[pack_index]
  extra_not_installed_packages = extra_packages[!pack_index]

  if(any(pack_index == FALSE) == TRUE) {
    cat("Consider to install these (optional) packages to run all possible Giotto commands for spatial analyses: ",
        extra_not_installed_packages)
    cat("\n Giotto does not automatically install all these packages as they are not absolutely required and this reduces the number of dependencies \n")
  }


  ## if cores is not set, then set number of cores automatically, but with limit of 10
  cores = determine_cores(cores)
  data.table::setDTthreads(threads = cores)


  ## expression data ##
  ## --------------- ##


  ## deprecated arguments
  if(!is.null(raw_exprs)) {
    expression = raw_exprs
    warning('raw_exprs argument is deprecated, use expression argument in the future \n')
  }

  if(!is.null(expression)) {

    gobject = extract_expression_list(gobject = gobject,
                                      expr_list = expression,
                                      expression_feat = expression_feat,
                                      cores = cores,
                                      verbose = verbose)
  }




  ## parameters ##
  ## ---------- ##
  gobject@parameters = list()


  ## set instructions ##
  ## ---------------- ##
  if(is.null(instructions)) {
    # create all default instructions
    gobject@instructions = createGiottoInstructions()
  }


  ## test if python modules are available
  python_modules = c('pandas', 'igraph', 'leidenalg', 'community', 'networkx', 'sklearn')
  my_python_path = gobject@instructions$python_path
  for(module in python_modules) {
    if(reticulate::py_module_available(module) == FALSE) {
      warning('module: ', module, ' was not found with python path: ', my_python_path, '\n')
    }
  }


  ## spatial locations ##
  ## ----------------- ##
  raw_cell_dim = ncol(gobject@expression[[1]][[1]]) # number of columns

  # list of spatial location data.table, each with a unique name
  # the default name = 'raw' and correspond to the real physical coordinates
  # additional spatial locations can be provided

  gobject = extract_spatial_locations_list(gobject,
                                           spat_loc_list = spatial_locs,
                                           cores = cores,
                                           dummy_n = raw_cell_dim,
                                           expr_matrix = gobject@expression[[1]][['raw']],
                                           verbose = TRUE)



  ## spatial info ##
  ## ------------ ##
  ## place to store segmentation info in polygon format style


  if(is.null(spatial_info)) {

    gobject@spatial_info = NULL

  } else {

    gobject = addGiottoPolygons(gobject = gobject,
                                gpolygons = spatial_info)

  }






  ## cell metadata ##
  ## ------------- ##
  if(is.null(cell_metadata)) {

    for(feat_type in expression_feat) {
      gobject@cell_metadata[[feat_type]] = data.table::data.table(cell_ID = gobject@cell_ID)
    }


  } else {

    if(length(cell_metadata) != length(expression_feat)) {
      stop('Number of different molecular features need to correspond with the cell_metadata list length \n')
    }

    for(feat_type in expression_feat) {
      gobject@cell_metadata[[feat_type]] = data.table::as.data.table(gobject@cell_metadata[[feat_type]])
      gobject@cell_metadata[[feat_type]][, cell_ID := gobject@cell_ID]

      # put cell_ID first
      all_colnames = colnames(gobject@cell_metadata[[feat_type]])
      other_colnames = grep('cell_ID', all_colnames, invert = T, value = T)
      gobject@cell_metadata[[feat_type]] = gobject@cell_metadata[[feat_type]][, c('cell_ID', other_colnames), with = FALSE]
    }
  }



  ## feat metadata ##
  ## ------------- ##
  if(is.null(feat_metadata)) {

    for(feat_type in expression_feat) {
      gobject@feat_metadata[[feat_type]] = data.table::data.table(feat_ID = gobject@feat_ID[[feat_type]])
    }

  } else {

    if(length(feat_metadata) != length(expression_feat)) {
      stop('Number of different molecular features need to correspond with the feat_metadata list length \n')
    }

    for(feat_type in expression_feat) {
      gobject@feat_metadata[[feat_type]] = data.table::as.data.table(gobject@feat_metadata[[feat_type]])
      gobject@feat_metadata[[feat_type]][, feat_ID := gobject@feat_ID[[feat_type]]]
    }

  }

  ## feature info ##
  ## ------------ ##
  ## place to store individual feature info
  if(is.null(feat_info)) {

    gobject@feat_info = NULL

  } else {

    gobject = addGiottoPoints(gobject = gobject,
                              gpoints = feat_info)

    #feat_info_names = names(feat_info)
    #for(feat_type in feat_info_names) {
    #  if(!feat_type %in% expression_feat) {
    #    warning('The feat info for ', feat_type, ' was not found back in ', expression_feat, ' and will not be used')
    #  } else {
    #    feat_info_object = feat_info[[feat_type]]
    #    feat_ids = gobject@feat_ID[[feat_type]]
    #    gobject@feat_info[[feat_type]] = feat_info_object
    #  }
    #}

  }




  ### OPTIONAL:
  ## spatial network
  if(!is.null(spatial_network)) {
    if(is.null(spatial_network_name) | length(spatial_network) != length(spatial_network_name)) {
      stop('\n each spatial network must be given a unique name \n')
    } else {

      for(network_i in 1:length(spatial_network)) {

        networkname = spatial_network_name[[network_i]]
        network     = spatial_network[[network_i]]

        if(any(c('data.frame', 'data.table') %in% class(network))) {
          if(all(c('to', 'from', 'weight', 'sdimx_begin', 'sdimy_begin', 'sdimx_end', 'sdimy_end') %in% colnames(network))) {
            spatial_network_Obj = create_spatialNetworkObject(name = networkname, networkDT = network)
            gobject@spatial_network[[networkname]] = spatial_network_Obj
          } else {
            stop('\n network ', networkname, ' does not have all necessary column names, see details \n')
          }
        } else {
          stop('\n network ', networkname, ' is not a data.frame or data.table \n')
        }
      }
    }
  }


  ## spatial grid
  if(!is.null(spatial_grid)) {
    if(is.null(spatial_grid_name) | length(spatial_grid) != length(spatial_grid_name)) {
      stop('\n each spatial grid must be given a unique name \n')
    } else {

      for(grid_i in 1:length(spatial_grid)) {

        gridname = spatial_grid_name[[grid_i]]
        grid     = spatial_grid[[grid_i]]

        if(any(c('data.frame', 'data.table') %in% class(grid))) {
          if(all(c('x_start', 'y_start', 'x_end', 'y_end', 'gr_name') %in% colnames(grid))) {
            gobject@spatial_grid[[gridname]] = grid
          } else {
            stop('\n grid ', gridname, ' does not have all necessary column names, see details \n')
          }
        } else {
          stop('\n grid ', gridname, ' is not a data.frame or data.table \n')
        }
      }
    }
  }

  ## spatial enrichment
  if(!is.null(spatial_enrichment)) {
    if(is.null(spatial_enrichment_name) | length(spatial_enrichment) != length(spatial_enrichment_name)) {
      stop('\n each spatial enrichment data.table or data.frame must be given a unique name \n')
    } else {

      for(spat_enrich_i in 1:length(spatial_enrichment)) {

        spatenrichname = spatial_enrichment_name[[spat_enrich_i]]
        spatenrich     = spatial_enrichment[[spat_enrich_i]]

        if(nrow(spatenrich) != nrow(gobject@cell_metadata)) {
          stop('\n spatial enrichment ', spatenrichname, ' does not have the same number of rows as spots/cells, see details \n')
        } else {

          gobject@spatial_enrichment[[spatenrichname]] = spatenrich

        }
      }
    }
  }


  ## dimension reduction
  if(!is.null(dimension_reduction)) {

    for(dim_i in 1:length(dimension_reduction)) {

      dim_red = dimension_reduction[[dim_i]]

      if(all(c('type', 'name', 'reduction_method', 'coordinates', 'misc') %in% names(dim_red))) {

        coord_data = dim_red[['coordinates']]

        if(all(rownames(coord_data) %in% gobject@cell_ID)) {

          type_value = dim_red[['type']] # cells or genes
          reduction_meth_value = dim_red[['reduction_method']] # e.g. umap, tsne, ...
          name_value = dim_red[['name']]  # uniq name
          misc_value = dim_red[['misc']]  # additional data

          gobject@dimension_reduction[[type_value]][[reduction_meth_value]][[name_value]] = dim_red[c('name', 'reduction_method', 'coordinates', 'misc')]
        } else {
          stop('\n rownames for coordinates are not found in gobject IDs \n')
        }

      } else {
        stop('\n each dimension reduction list must contain all required slots, see details. \n')
      }

    }

  }

  # NN network
  if(!is.null(nn_network)) {

    for(nn_i in 1:length(nn_network)) {

      nn_netw = nn_network[[nn_i]]

      if(all(c('type', 'name', 'igraph') %in% names(nn_netw))) {

        igraph_data = nn_netw[['igraph']]

        if(all(names(igraph::V(igraph_data)) %in% gobject@cell_ID)) {

          type_value = nn_netw[['type']] # sNN or kNN
          name_value = nn_netw[['name']]  # uniq name

          gobject@nn_network[[type_value]][[name_value]][['igraph']] = igraph_data
        } else {
          stop('\n igraph vertex names are not found in gobject IDs \n')
        }

      } else {
        stop('\n each nn network list must contain all required slots, see details. \n')
      }

    }

  }

  ## images ##
  # expect a list of giotto object images
  if(!is.null(images)) {

    if(is.null(names(images))) {
      names(images) = paste0('image.', 1:length(images))
    }

    for(image_i in 1:length(images)) {

      im = images[[image_i]]
      im_name = names(images)[[image_i]]

      if(methods::is(im, 'giottoImage')) {
        gobject@images[[im_name]] = im
      } else {
        warning('image: ', im, ' is not a giotto image object')
      }

    }

  }



  # other information
  # TODO

  return(gobject)

}


#' @name get_img_minmax
#' @keywords internal
get_img_minmax = function(mg_img) {
  #Get magick object dimensions. xmin and ymax assumed to be 0.
  info = magick::image_info(mg_img)
  img_xmax = info$width     #width
  img_xmin = 0              #x origin
  img_ymax = 0              #y origin
  img_ymin = -(info$height) #height

  return(list('img_xmax' = img_xmax,
              'img_xmin' = img_xmin,
              'img_ymax' = img_ymax,
              'img_ymin' = img_ymin))
}


#' @name get_adj_rescale_img
#' @keywords internal
get_adj_rescale_img = function(img_minmax,
                               spatial_locs,
                               scale_factor = 1) {

  #Spatial minmax
  my_xmin = min(spatial_locs$sdimx)
  my_xmax = max(spatial_locs$sdimx)
  my_ymin = min(spatial_locs$sdimy)
  my_ymax = max(spatial_locs$sdimy)

  #Find scaled image adjustments based on scaled spatlocs
  xmin_adj_scaled = (my_xmin*scale_factor) - (img_minmax$img_xmin)
  xmin_adj_orig = xmin_adj_scaled/scale_factor

  xmax_adj_scaled = (img_minmax$img_xmax) - (my_xmax*scale_factor)
  xmax_adj_orig = xmax_adj_scaled/scale_factor

  ymin_adj_scaled = (my_ymin*scale_factor) - (img_minmax$img_ymin)
  ymin_adj_orig = ymin_adj_scaled/scale_factor

  ymax_adj_scaled = (img_minmax$img_ymax) - (my_ymax*scale_factor)
  ymax_adj_orig = ymax_adj_scaled/scale_factor

  #return scaled adjustments
  return(c('xmin_adj_orig' = xmin_adj_orig,
           'xmax_adj_orig' = xmax_adj_orig,
           'ymin_adj_orig' = ymin_adj_orig,
           'ymax_adj_orig' = ymax_adj_orig))

}





#' @title createGiottoVisiumObject
#' @description creates Giotto object directly from a 10X visium folder
#' @param visium_dir path to the 10X visium directory [required]
#' @param expr_data raw or filtered data (see details)
#' @param gene_column_index which column index to select (see details)
#' @param h5_visium_path path to visium 10X .h5 file
#' @param h5_gene_ids gene names as symbols (default) or ensemble gene ids
#' @param h5_tissue_positions_path path to tissue locations (.csv file)
#' @param h5_image_png_path path to tissue .png file (optional)
#' @param h5_json_scalefactors_path path to .json scalefactors (optional)
#' @param png_name select name of png to use (see details)
#' @param do_manual_adj flag to use manual adj values instead of automatic image alignment
#' @param xmax_adj adjustment of the maximum x-value to align the image
#' @param xmin_adj adjustment of the minimum x-value to align the image
#' @param ymax_adj adjustment of the maximum y-value to align the image
#' @param ymin_adj adjustment of the minimum y-value to align the image
#' @param instructions list of instructions or output result from \code{\link{createGiottoInstructions}}
#' @param cores how many cores or threads to use to read data if paths are provided
#' @param verbose be verbose
#' @return giotto object
#' @details
#' If starting from a Visium 10X directory:
#' \itemize{
#'   \item{expr_data: raw will take expression data from raw_feature_bc_matrix and filter from filtered_feature_bc_matrix}
#'   \item{gene_column_index: which gene identifiers (names) to use if there are multiple columns (e.g. ensemble and gene symbol)}
#'   \item{png_name: by default the first png will be selected, provide the png name to override this (e.g. myimage.png)}
#'   \item{the file scalefactors_json.json will be detected automaticated and used to attempt to align the data}
#' }
#'
#' If starting from a Visium 10X .h5 file
#' \itemize{
#'   \item{h5_visium_path: full path to .h5 file: /your/path/to/visium_file.h5}
#'   \item{h5_tissue_positions_path: full path to spatial locations file: /you/path/to/tissue_positions_list.csv}
#'   \item{h5_image_png_path: full path to png: /your/path/to/images/tissue_lowres_image.png}
#'   \item{h5_json_scalefactors_path: full path to .json file: /your/path/to/scalefactors_json.json}
#' }
#'
#' @export
createGiottoVisiumObject = function(visium_dir = NULL,
                                    expr_data = c('raw', 'filter'),
                                    gene_column_index = 1,
                                    h5_visium_path = NULL,
                                    h5_gene_ids = c('symbols', 'ensembl'),
                                    h5_tissue_positions_path = NULL,
                                    h5_image_png_path = NULL,
                                    h5_json_scalefactors_path = NULL,
                                    png_name = NULL,
                                    do_manual_adj = FALSE,
                                    xmax_adj = 0,
                                    xmin_adj = 0,
                                    ymax_adj = 0,
                                    ymin_adj = 0,
                                    instructions = NULL,
                                    cores = NA,
                                    verbose = TRUE) {

  if(!is.null(h5_visium_path)) {

    if(verbose) cat("A path to an .h5 10X file was provided and will be used \n")

    if(!file.exists(h5_visium_path)) stop("The provided path ", h5_visium_path, " does not exist \n")

    # spatial locations
    if(is.null(h5_tissue_positions_path)) stop("A path to the tissue positions (.csv) needs to be provided to h5_tissue_positions_path \n")
    if(!file.exists(h5_tissue_positions_path)) stop("The provided path ", h5_tissue_positions_path, " does not exist \n")

    # get matrix counts
    h5_results = get10Xmatrix_h5(path_to_data = h5_visium_path, gene_ids = h5_gene_ids)
    raw_matrix = h5_results[['Gene Expression']]

    # spatial locations
    spatial_results = data.table::fread(h5_tissue_positions_path)
    spatial_results = spatial_results[match(colnames(raw_matrix), V1)]
    colnames(spatial_results) = c('barcode', 'in_tissue', 'array_row', 'array_col', 'col_pxl', 'row_pxl')
    spatial_locs = spatial_results[,.(row_pxl,-col_pxl)]
    colnames(spatial_locs) = c('sdimx', 'sdimy')

    # image (optional)
    if(!is.null(h5_image_png_path)) {


      if(!file.exists(h5_image_png_path)) stop("The provided path ", h5_image_png_path,
                                               " does not exist. Set to NULL to exclude or provide the correct path. \n")

      mg_img = magick::image_read(h5_image_png_path)

      ## check if automatic alignment can be done
      png_name = basename(h5_image_png_path)

      if(png_name == 'tissue_lowres_image.png') {
        if(file.exists(h5_json_scalefactors_path)) {
          if(verbose == TRUE && do_manual_adj == FALSE) cat('png and scalefactors paths are found and automatic alignment for the lowres image will be attempted \n')

          json_info = jsonlite::read_json(h5_json_scalefactors_path)
          scale_factor = json_info[['tissue_lowres_scalef']]

          visium_png = createGiottoImage(gobject = NULL,
                                         spatial_locs = spatial_locs,
                                         mg_object = mg_img,
                                         name = 'image',
                                         scale_factor = scale_factor,
                                         do_manual_adj = do_manual_adj,
                                         xmax_adj = xmax_adj,
                                         xmin_adj = xmin_adj,
                                         ymax_adj = ymax_adj,
                                         ymin_adj = ymin_adj)

        }
      } else if(png_name == 'tissue_hires_image.png') {
        if(file.exists(h5_json_scalefactors_path)) {
          if(verbose == TRUE && do_manual_adj == FALSE) cat('png and scalefactors paths are found and automatic alignment for the hires image will be attempted \n')

          json_info = jsonlite::read_json(h5_json_scalefactors_path)
          scale_factor = json_info[['tissue_hires_scalef']]

          visium_png = createGiottoImage(gobject = NULL,
                                         spatial_locs = spatial_locs,
                                         mg_object = mg_img,
                                         name = 'image',
                                         scale_factor = scale_factor,
                                         do_manual_adj = do_manual_adj,
                                         xmax_adj = xmax_adj,
                                         xmin_adj = xmin_adj,
                                         ymax_adj = ymax_adj,
                                         ymin_adj = ymin_adj)

        }
      } else {
        visium_png = createGiottoImage(gobject = NULL, spatial_locs =  spatial_locs,
                                       mg_object = mg_img, name = 'image',
                                       xmax_adj = xmax_adj, xmin_adj = xmin_adj,
                                       ymax_adj = ymax_adj, ymin_adj = ymin_adj)
      }

      visium_png_list = list(visium_png)
      names(visium_png_list) = c('image')
    } else {
      visium_png_list = NULL
    }

    # create Giotto object
    giotto_object = createGiottoObject(expression = raw_matrix,
                                       expression_feat = 'rna',
                                       spatial_locs = spatial_locs,
                                       instructions = instructions,
                                       cell_metadata = list('rna' = spatial_results[,.(in_tissue, array_row, array_col)]),
                                       images = visium_png_list)
    return(giotto_object)


  } else {

    if(verbose) cat("A structured visium directory will be used \n")

    # data.table: set global variable
    V1 = row_pxl = col_pxl = in_tissue = array_row = array_col = NULL

    ## check arguments
    if(is.null(visium_dir)) stop('visium_dir needs to be a path to a visium directory \n')
    if(!file.exists(visium_dir)) stop(visium_dir, ' does not exist \n')
    expr_data = match.arg(expr_data, choices = c('raw', 'filter'))

    # set number of cores automatically, but with limit of 10
    cores = determine_cores(cores)
    data.table::setDTthreads(threads = cores)

    ## matrix
    if(expr_data == 'raw') {
      data_path = paste0(visium_dir, '/', 'raw_feature_bc_matrix/')
      raw_matrix = get10Xmatrix(path_to_data = data_path, gene_column_index = gene_column_index)
    } else if(expr_data == 'filter') {
      data_path = paste0(visium_dir, '/', 'filtered_feature_bc_matrix/')
      raw_matrix = get10Xmatrix(path_to_data = data_path, gene_column_index = gene_column_index)
    }

    ## spatial locations and image
    spatial_path = paste0(visium_dir, '/', 'spatial/')
    spatial_results = data.table::fread(paste0(spatial_path, '/','tissue_positions_list.csv'))
    spatial_results = spatial_results[match(colnames(raw_matrix), V1)]
    colnames(spatial_results) = c('barcode', 'in_tissue', 'array_row', 'array_col', 'col_pxl', 'row_pxl')
    spatial_locs = spatial_results[,.(row_pxl,-col_pxl)]
    colnames(spatial_locs) = c('sdimx', 'sdimy')

    ## spatial image
    if(is.null(png_name)) {
      png_list = list.files(spatial_path, pattern = "*.png")
      png_name = png_list[1]
    }
    png_path = paste0(spatial_path,'/',png_name)
    if(!file.exists(png_path)) stop(png_path, ' does not exist! \n')

    mg_img = magick::image_read(png_path)


    if(png_name == 'tissue_lowres_image.png') {

      scalefactors_path = paste0(spatial_path,'/','scalefactors_json.json')

      if(file.exists(scalefactors_path)) {
        if(verbose == TRUE && do_manual_adj == FALSE) cat('png and scalefactors paths are found and automatic alignment for the lowres image will be attempted \n')

        json_info = jsonlite::read_json(scalefactors_path)
        scale_factor = json_info[['tissue_lowres_scalef']]

        visium_png = createGiottoImage(gobject = NULL,
                                       spatial_locs = spatial_locs,
                                       mg_object = mg_img,
                                       name = 'image',
                                       scale_factor = scale_factor,
                                       do_manual_adj = do_manual_adj,
                                       xmax_adj = xmax_adj,
                                       xmin_adj = xmin_adj,
                                       ymax_adj = ymax_adj,
                                       ymin_adj = ymin_adj)

      }
    } else if(png_name == 'tissue_hires_image.png') {

      scalefactors_path = paste0(spatial_path,'/','scalefactors_json.json')

       if(file.exists(scalefactors_path)) {
        if(verbose == TRUE && do_manual_adj == FALSE) cat('png and scalefactors paths are found and automatic alignment for the hires image will be attempted \n')

        json_info = jsonlite::read_json(scalefactors_path)
        scale_factor = json_info[['tissue_hires_scalef']]

        visium_png = createGiottoImage(gobject = NULL,
                                       spatial_locs = spatial_locs,
                                       mg_object = mg_img,
                                       name = 'image',
                                       scale_factor = scale_factor,
                                       do_manual_adj = do_manual_adj,
                                       xmax_adj = xmax_adj,
                                       xmin_adj = xmin_adj,
                                       ymax_adj = ymax_adj,
                                       ymin_adj = ymin_adj)

      }
    } else {
      visium_png = createGiottoImage(gobject = NULL, spatial_locs =  spatial_locs,
                                     mg_object = mg_img, name = 'image',
                                     xmax_adj = xmax_adj, xmin_adj = xmin_adj,
                                     ymax_adj = ymax_adj, ymin_adj = ymin_adj)
    }

    visium_png_list = list(visium_png)
    names(visium_png_list) = c('image')

    giotto_object = createGiottoObject(expression = raw_matrix,
                                       expression_feat = 'rna',
                                       spatial_locs = spatial_locs,
                                       instructions = instructions,
                                       cell_metadata = list('rna' = spatial_results[,.(in_tissue, array_row, array_col)]),
                                       images = visium_png_list)
    return(giotto_object)

  }

}




#' @name createGiottoObjectSubcellular
#' @description Function to create a giotto object starting from subcellular polygon (e.g. cell) and points (e.g. transcripts) information
#' @param gpoints giotto points
#' @param gpolygons giotto polygons
#' @param cell_metadata cell annotation metadata
#' @param feat_metadata feature annotation metadata for each unique feature
#' @param spatial_network list of spatial network(s)
#' @param spatial_network_name list of spatial network name(s)
#' @param spatial_grid list of spatial grid(s)
#' @param spatial_grid_name list of spatial grid name(s)
#' @param spatial_enrichment list of spatial enrichment score(s) for each spatial region
#' @param spatial_enrichment_name list of spatial enrichment name(s)
#' @param dimension_reduction list of dimension reduction(s)
#' @param nn_network list of nearest neighbor network(s)
#' @param images list of images
#' @param instructions list of instructions or output result from \code{\link{createGiottoInstructions}}
#' @param cores how many cores or threads to use to read data if paths are provided
#' @param verbose be verbose when building Giotto object
#' @return giotto object
#' @keywords giotto
#' @export
createGiottoObjectSubcellular = function(gpoints = NULL,
                                         gpolygons = NULL,
                                         cell_metadata = NULL,
                                         feat_metadata = NULL,
                                         spatial_network = NULL,
                                         spatial_network_name = NULL,
                                         spatial_grid = NULL,
                                         spatial_grid_name = NULL,
                                         spatial_enrichment = NULL,
                                         spatial_enrichment_name = NULL,
                                         dimension_reduction = NULL,
                                         nn_network = NULL,
                                         images = NULL,
                                         instructions = NULL,
                                         cores = NA,
                                         verbose = TRUE) {


  # create minimum giotto
  gobject = giotto(expression = NULL,
                   expression_feat = NULL,
                   spatial_locs = NULL,
                   spatial_info = NULL,
                   cell_metadata = NULL,
                   feat_metadata = NULL,
                   feat_info = NULL,
                   cell_ID = NULL,
                   feat_ID = NULL,
                   spatial_network = NULL,
                   spatial_grid = NULL,
                   spatial_enrichment = NULL,
                   dimension_reduction = NULL,
                   nn_network = NULL,
                   images = NULL,
                   parameters = NULL,
                   offset_file = NULL,
                   instructions = instructions,
                   OS_platform = .Platform[['OS.type']])


  ## if cores is not set, then set number of cores automatically, but with limit
  cores = determine_cores(cores)
  data.table::setDTthreads(threads = cores)



  # gpolygons and gpoints need to be provided
  if(is.null(gpolygons)) {
    stop('gpolygons = NULL, cell polygon information needs to be given (e.g. cell boundary, nucleus, ...)')
  }

  if(is.null(gpoints)) {
    stop('gpoints = NULL, feature information needs to be given (e.g. transcript or protein location)')
  }


  ## extract polygon information ##
  ## --------------------------- ##
  polygon_res = extract_polygon_list(polygonlist = gpolygons)
  gobject@spatial_info = polygon_res

  ## cell ID ##
  ## ------- ##
  gobject@cell_ID = gobject@spatial_info[['cell']]@spatVector$poly_ID

  ## extract points information ##
  ## -------------------------- ##
  points_res = extract_points_list(pointslist = gpoints)
  gobject@feat_info = points_res

  ## expression features ##
  ## ------------------- ##
  gobject@expression_feat = names(points_res)
  expression_feat = gobject@expression_feat

  ## feat ID ##
  ## ------- ##
  for(feat in gobject@expression_feat) {
    unique_feats = unique(gobject@feat_info[[feat]]@spatVector$feat_ID)
    gobject@feat_ID[[feat]] = unique_feats
  }


  ## parameters ##
  ## ---------- ##
  gobject@parameters = list()

  ## set instructions ##
  ## ---------------- ##
  if(is.null(instructions)) {
    # create all default instructions
    gobject@instructions = createGiottoInstructions()
  }



  ## cell metadata ##
  ## ------------- ##
  if(is.null(cell_metadata)) {

    for(feat_type in expression_feat) {
      gobject@cell_metadata[[feat_type]] = data.table::data.table(cell_ID = gobject@cell_ID)
    }


  } else {

    if(length(cell_metadata) != length(expression_feat)) {
      stop('Number of different molecular features need to correspond with the cell_metadata list length \n')
    }

    for(feat_type in expression_feat) {
      gobject@cell_metadata[[feat_type]] = data.table::as.data.table(gobject@cell_metadata[[feat_type]])
      gobject@cell_metadata[[feat_type]][, cell_ID := gobject@cell_ID]

      # put cell_ID first
      all_colnames = colnames(gobject@cell_metadata[[feat_type]])
      other_colnames = grep('cell_ID', all_colnames, invert = T, value = T)
      gobject@cell_metadata[[feat_type]] = gobject@cell_metadata[[feat_type]][, c('cell_ID', other_colnames), with = FALSE]
    }
  }


  ## feat metadata ##
  ## ------------- ##
  if(is.null(feat_metadata)) {

    for(feat_type in expression_feat) {
      gobject@feat_metadata[[feat_type]] = data.table::data.table(feat_ID = gobject@feat_ID[[feat_type]])
    }

  } else {

    if(length(feat_metadata) != length(expression_feat)) {
      stop('Number of different molecular features need to correspond with the feat_metadata list length \n')
    }

    for(feat_type in expression_feat) {
      gobject@feat_metadata[[feat_type]] = data.table::as.data.table(gobject@feat_metadata[[feat_type]])
      gobject@feat_metadata[[feat_type]][, feat_ID := gobject@feat_ID[[feat_type]]]
    }

  }

  ### OPTIONAL:
  ## spatial network
  if(!is.null(spatial_network)) {
    if(is.null(spatial_network_name) | length(spatial_network) != length(spatial_network_name)) {
      stop('\n each spatial network must be given a unique name \n')
    } else {

      for(network_i in 1:length(spatial_network)) {

        networkname = spatial_network_name[[network_i]]
        network     = spatial_network[[network_i]]

        if(any(c('data.frame', 'data.table') %in% class(network))) {
          if(all(c('to', 'from', 'weight', 'sdimx_begin', 'sdimy_begin', 'sdimx_end', 'sdimy_end') %in% colnames(network))) {
            spatial_network_Obj = create_spatialNetworkObject(name = networkname, networkDT = network)
            gobject@spatial_network[[networkname]] = spatial_network_Obj
          } else {
            stop('\n network ', networkname, ' does not have all necessary column names, see details \n')
          }
        } else {
          stop('\n network ', networkname, ' is not a data.frame or data.table \n')
        }
      }
    }
  }


  ## spatial grid
  if(!is.null(spatial_grid)) {
    if(is.null(spatial_grid_name) | length(spatial_grid) != length(spatial_grid_name)) {
      stop('\n each spatial grid must be given a unique name \n')
    } else {

      for(grid_i in 1:length(spatial_grid)) {

        gridname = spatial_grid_name[[grid_i]]
        grid     = spatial_grid[[grid_i]]

        if(any(c('data.frame', 'data.table') %in% class(grid))) {
          if(all(c('x_start', 'y_start', 'x_end', 'y_end', 'gr_name') %in% colnames(grid))) {
            gobject@spatial_grid[[gridname]] = grid
          } else {
            stop('\n grid ', gridname, ' does not have all necessary column names, see details \n')
          }
        } else {
          stop('\n grid ', gridname, ' is not a data.frame or data.table \n')
        }
      }
    }
  }

  ## spatial enrichment
  if(!is.null(spatial_enrichment)) {
    if(is.null(spatial_enrichment_name) | length(spatial_enrichment) != length(spatial_enrichment_name)) {
      stop('\n each spatial enrichment data.table or data.frame must be given a unique name \n')
    } else {

      for(spat_enrich_i in 1:length(spatial_enrichment)) {

        spatenrichname = spatial_enrichment_name[[spat_enrich_i]]
        spatenrich     = spatial_enrichment[[spat_enrich_i]]

        if(nrow(spatenrich) != nrow(gobject@cell_metadata)) {
          stop('\n spatial enrichment ', spatenrichname, ' does not have the same number of rows as spots/cells, see details \n')
        } else {

          gobject@spatial_enrichment[[spatenrichname]] = spatenrich

        }
      }
    }
  }


  ## dimension reduction
  if(!is.null(dimension_reduction)) {

    for(dim_i in 1:length(dimension_reduction)) {

      dim_red = dimension_reduction[[dim_i]]

      if(all(c('type', 'name', 'reduction_method', 'coordinates', 'misc') %in% names(dim_red))) {

        coord_data = dim_red[['coordinates']]

        if(all(rownames(coord_data) %in% gobject@cell_ID)) {

          type_value = dim_red[['type']] # cells or genes
          reduction_meth_value = dim_red[['reduction_method']] # e.g. umap, tsne, ...
          name_value = dim_red[['name']]  # uniq name
          misc_value = dim_red[['misc']]  # additional data

          gobject@dimension_reduction[[type_value]][[reduction_meth_value]][[name_value]] = dim_red[c('name', 'reduction_method', 'coordinates', 'misc')]
        } else {
          stop('\n rownames for coordinates are not found in gobject IDs \n')
        }

      } else {
        stop('\n each dimension reduction list must contain all required slots, see details. \n')
      }

    }

  }

  # NN network
  if(!is.null(nn_network)) {

    for(nn_i in 1:length(nn_network)) {

      nn_netw = nn_network[[nn_i]]

      if(all(c('type', 'name', 'igraph') %in% names(nn_netw))) {

        igraph_data = nn_netw[['igraph']]

        if(all(names(igraph::V(igraph_data)) %in% gobject@cell_ID)) {

          type_value = nn_netw[['type']] # sNN or kNN
          name_value = nn_netw[['name']]  # uniq name

          gobject@nn_network[[type_value]][[name_value]][['igraph']] = igraph_data
        } else {
          stop('\n igraph vertex names are not found in gobject IDs \n')
        }

      } else {
        stop('\n each nn network list must contain all required slots, see details. \n')
      }

    }

  }

  ## images ##
  # expect a list of giotto object images
  if(!is.null(images)) {

    if(is.null(names(images))) {
      names(images) = paste0('image.', 1:length(images))
    }

    for(image_i in 1:length(images)) {

      im = images[[image_i]]
      im_name = names(images)[[image_i]]

      if(methods::is(im, 'giottoImage')) {
        gobject@images[[im_name]] = im
      } else {
        warning('image: ', im, ' is not a giotto image object')
      }

    }

  }

  return(gobject)


}






#### logging of giotto functions ####

#' @name get_args
#' @keywords internal
get_args <- function(toplevel = 2, verbose = FALSE) {

  nframes = sys.nframe()

  if(verbose == TRUE) {
    cat('\n number of frames: ')
    print(nframes)
    cat('\n')
  }


  cl = sys.call(-toplevel)

  if(verbose == TRUE) {
    cat('\n system call: ')
    print(cl)
    cat('\n')
  }


  # function name
  fname = as.character(cl[[1]])

  if(length(fname) > 1) {
    fname = fname[[3]]
  }

  if(verbose == TRUE) {
    cat('\n function name: ')
    print(fname)
    cat('\n')
  }


  # function
  #f = get(x = fname, mode = "function", pos = 'package:Giotto')
  f = get(x = fname, mode = "function", pos = sys.frame(-2))

  # get used arguments
  cl = match.call(definition=f, call=cl)
  user_args = as.list(cl)[-1]

  # all fun arguments
  fun_args <- formals(fun = fname)
  fun_args[names(user_args)] = user_args

  unl_args = unlist(fun_args)
  final_args = as.character(unl_args)
  names(final_args) = names(unl_args)

  # select first from vector
  bool_det = grepl("c\\(", final_args)
  if(any(bool_det) == TRUE) {

    for(bool_name in names(final_args[bool_det])) {

      bool_vec = final_args[bool_name]
      new_vec = strsplit(bool_vec, split = "\"")[[1]][2]

      final_args[bool_name] = new_vec

    }
  }

  return(final_args)

}



#' @name update_giotto_params
#' @keywords internal
update_giotto_params = function(gobject,
                                description = '_test',
                                return_gobject = TRUE,
                                toplevel = 2) {

  parameters_list = gobject@parameters
  number_of_rounds = length(parameters_list)
  update_name = paste0(number_of_rounds, description)

  parameters_list[[update_name]] = get_args(toplevel = toplevel)

  if(return_gobject == TRUE) {
    gobject@parameters = parameters_list
    return(gobject)
  } else {
    return(list(plist = parameters_list, newname = update_name))
  }
}



#### joining giotto object ####

#' @name join_expression_matrices
#' @keywords internal
join_expression_matrices = function(matrix_list) {

  # find all features
  final_feats = list()
  for(matr_i in 1:length(matrix_list)) {
    rowfeats = rownames(matrix_list[[matr_i]])
    final_feats[[matr_i]] = rowfeats
  }

  final_feats = unique(unlist(final_feats))
  final_feats = sort(final_feats)



  # extend matrices with missing ids
  final_mats = list()
  for(matr_i in 1:length(matrix_list)) {
    matr = matrix_list[[matr_i]]

    missing_feats = final_feats[!final_feats %in% rownames(matr)]

    missing_mat = Matrix::Matrix(data = 0,
                                 nrow = length(missing_feats), ncol = ncol(matr),
                                 dimnames = list(missing_feats, colnames(matr)))

    mat_ext = rbind(matr, missing_mat)
    mat_ext = mat_ext[match(final_feats, rownames(mat_ext)), ]

    final_mats[[matr_i]] = mat_ext

  }

  combined_matrix = do.call('cbind', final_mats)
  return(list(matrix = combined_matrix, sort_all_feats = final_feats))
}

#' @name join_spatlocs
#' @keywords internal
join_spatlocs = function(dt_list) {

  final_list = do.call('rbind', dt_list)
  return(final_list)
}

#' @name join_cell_meta
#' @keywords internal
join_cell_meta = function(dt_list) {

  final_list = do.call('rbind', dt_list)
  return(final_list)

}


#' @name joinGiottoObjects
#' @description Function to join multiple giotto objects together
#' @param gobject_list list of giotto objects
#' @param gobject_names unique giotto names for each giotto object
#' @param join_method method to join giotto objects
#' @param z_step distance along z-axis if method is z-stack
#' @param x_shift shift along x-axis if method is x-shift
#' @param x_padding padding bewteen datasets/images if method is x-shift
#' @param verbose be verbose
#' @return giotto object
#' @keywords giotto
#' @export
joinGiottoObjects = function(gobject_list,
                             gobject_names = NULL,
                             join_method = c('x_shift', 'z_stack'),
                             z_step = 1000,
                             x_shift = NULL,
                             x_padding = 0,
                             verbose = TRUE) {


  ## check params
  if(!is.vector(gobject_names) | !is.character(gobject_names)) {
    stop('gobject_names need to be a vector with unique names for the giotto objects')
  }

  if(length(gobject_list) != length(gobject_names)) {
    stop('each giotto object in the list needs to have a unique (short) name')
  }

  join_method = match.arg(arg = join_method, choices = c('x_shift', 'z_stack'))

  # keep instructions from first giotto object
  first_instructions = gobject_list[[1]]@instructions

  # keep features from first giotto object
  first_features = gobject_list[[1]]@expression_feat


  updated_object_list = list()




  ## 0. re-scale spatial locations ##
  ## ----------------------------- ##






  ## 1. update giotto objects ##
  ## ------------------------ ##
  if(verbose == TRUE) cat('start updating objects \n')

  all_cell_ID_list = list()
  all_image_list = list()
  xshift_list = list()

  for(gobj_i in 1:length(gobject_list)) {

    gobj = gobject_list[[gobj_i]]
    gname = gobject_names[[gobj_i]]

    ## 0. update cell ID
    gobj@cell_ID = paste0(gname,'-',gobj@cell_ID)

    all_cell_ID_list[[gobj_i]] = gobj@cell_ID


    ## 1. update expression
    # provide unique cell ID name
    for(feat in gobj@expression_feat) {

      for(matr in names(gobj@expression[[feat]])) {
        colnames(gobj@expression[[feat]][[matr]]) = gobj@cell_ID
      }

    }


    ## 2. update images
    # change individual names

    images_found = !is.null(gobj@images)

    if(images_found) {

      names(gobj@images) = paste0(gname, '-', names(gobj@images))
      for(imname in names(gobj@images)) {

        gobj@images[[imname]]@name = paste0(gname,'-', gobj@images[[imname]]@name)


        if(join_method == 'x_shift') {

          if(is.null(x_shift)) {

            # estimate x_shift step directly from giotto image
            gimage = gobj@images[[imname]]

            my_xmax = gimage@minmax[1]
            my_xmin = gimage@minmax[2]
            xmax_b = gimage@boundaries[1]
            xmin_b = gimage@boundaries[2]
            xmin = my_xmin-xmin_b
            xmax = my_xmax+xmax_b

            add_to_x = ((gobj_i - 1) * (xmax-xmin)) + ((gobj_i - 1) * x_padding)

          } else {
            add_to_x = ((gobj_i - 1) * x_shift) + ((gobj_i - 1) * x_padding)
          }


          gobj@images[[imname]]@minmax[c("xmax_sloc", "xmin_sloc")] =  gobj@images[[imname]]@minmax[c("xmax_sloc", "xmin_sloc")] + add_to_x

          xshift_list[[gobj_i]] = add_to_x
        }

        all_image_list[[imname]] = gobj@images[[imname]]

      }


    }



    ## 3. update spatial location
    # add padding to x-axis
    # update cell ID
    for(locs in names(gobj@spatial_locs)) {
      myspatlocs = gobj@spatial_locs[[locs]]

      if(join_method == 'z_stack') {
        myspatlocs[, sdimz := (gobj_i - 1) * z_step]
        myspatlocs[, cell_ID := gobj@cell_ID]
        myspatlocs = myspatlocs[,.(sdimx, sdimy, sdimz, cell_ID)]
      } else if(join_method == 'x_shift') {

        if(is.null(x_shift)) {
          add_to_x = xshift_list[[gobj_i]]

        } else {
          add_to_x = (gobj_i - 1) * x_shift + ((gobj_i - 1) * x_padding)
        }

        myspatlocs[, sdimx := sdimx+add_to_x]
        myspatlocs[, cell_ID := gobj@cell_ID]
      }

      gobj@spatial_locs[[locs]] = myspatlocs
    }

    # cell metadata
    # rbind metadata
    # create capture area specific names
    for(feat in names(gobj@cell_metadata)) {
      gobj@cell_metadata[[feat]][['cell_ID']] = gobj@cell_ID
      gobj@cell_metadata[[feat]][['list_ID']] = gname
    }



    updated_object_list[[gobj_i]] = gobj

  }

  #return(updated_object_list)




  ## 2. prepare for new giotto object ##
  ## -------------------------------- ##
  comb_gobject = Giotto:::giotto(expression = list(),
                                 expression_feat = first_features,
                                 spatial_locs = NULL,
                                 spatial_info = NULL,
                                 cell_metadata = NULL,
                                 feat_metadata = NULL,
                                 feat_info = NULL,
                                 cell_ID = NULL,
                                 feat_ID = NULL,
                                 spatial_network = NULL,
                                 spatial_grid = NULL,
                                 spatial_enrichment = NULL,
                                 dimension_reduction = NULL,
                                 nn_network = NULL,
                                 images = NULL,
                                 parameters = NULL,
                                 offset_file = NULL,
                                 instructions = first_instructions,
                                 OS_platform = .Platform[['OS.type']],
                                 join_info = NULL)




  ## 3. merge updated data  ##
  ## ------------------------ ##

  first_obj = updated_object_list[[1]]

  ## cell IDs
  combined_cell_ID = unlist(all_cell_ID_list)
  comb_gobject@cell_ID = combined_cell_ID


  ## expression
  if(verbose == TRUE) cat('start expression combination \n')
  expr_names = names(first_obj@expression)

  for(name in expr_names) {

    for(mode in names(first_obj@expression[[name]])) {

      savelist = list()
      for(gobj_i in 1:length(updated_object_list)) {

        mat = updated_object_list[[gobj_i]]@expression[[name]][[mode]]
        savelist[[gobj_i]] = mat
      }

      combmat = join_expression_matrices(matrix_list = savelist)
      comb_gobject@expression[[name]][[mode]] = combmat$matrix

      comb_gobject@feat_ID[[name]] = combmat$sort_all_feats
      comb_gobject@feat_metadata[[name]] = data.table::data.table(feat_ID = combmat$sort_all_feats)
    }

  }


  ## spatial locations
  if(verbose == TRUE) cat('start spatial location combination \n')
  spatloc_names = names(first_obj@spatial_locs)

  for(name in spatloc_names) {

    savelist = list()
    for(gobj_i in 1:length(updated_object_list)) {
      spatlocs = updated_object_list[[gobj_i]]@spatial_locs[[name]]
      savelist[[gobj_i]] = spatlocs
    }

    combspatlocs = join_spatlocs(dt_list = savelist)
    comb_gobject@spatial_locs[[name]] = combspatlocs
  }


  ## cell metadata
  if(verbose == TRUE) cat('start cell metadata combination \n')
  metanames = names(first_obj@cell_metadata)

  for(name in metanames) {

    savelist = list()
    for(gobj_i in 1:length(updated_object_list)) {
      cellmeta = updated_object_list[[gobj_i]]@cell_metadata[[name]]
      savelist[[gobj_i]] = cellmeta
    }

    combcellmeta = join_cell_meta(dt_list = savelist)
    comb_gobject@cell_metadata[[name]] = combcellmeta

  }


  ## images
  if(verbose == TRUE) cat('start image')

  # keep individual images
  # each individual image has updated x and y locations
  # so all images can be viewed together by plotting them one-by-one
  # but images can also be easify viewed separately by grouping them
  comb_gobject@images = all_image_list


  ## TODO:
  # update giotto object with join-information
  # - list ID names
  # - xshift values

  # add option to perform yshift

  comb_gobject@join_info = list(list_IDs = gobject_names,
                                join_method = join_method,
                                z_step = z_step,
                                x_shift = x_shift,
                                x_padding = x_padding)

  return(comb_gobject)


}



