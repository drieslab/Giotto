

#' @title S4 giotto Class
#' @description Framework of giotto object to store and work with spatial expression data
#' @keywords giotto, object
#' @slot raw_exprs raw expression counts
#' @slot norm_expr normalized expression counts
#' @slot norm_scaled_expr normalized and scaled expression counts
#' @slot custom_expr custom normalized counts
#' @slot spatial_locs spatial location coordinates for cells
#' @slot cell_metadata metadata for cells
#' @slot gene_metadata metadata for genes
#' @slot cell_ID unique cell IDs
#' @slot gene_ID unique gene IDs
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
#' @useDynLib Giotto
#' @importFrom Rcpp sourceCpp
#' @export
giotto <- setClass(
  "giotto",
  slots = c(
    raw_exprs = "ANY",
    norm_expr = "ANY",
    norm_scaled_expr = "ANY",
    custom_expr = "ANY",
    spatial_locs = "ANY",
    cell_metadata = "ANY",
    gene_metadata = "ANY",
    cell_ID = "ANY",
    gene_ID = "ANY",
    spatial_network = "ANY",
    spatial_grid = "ANY",
    spatial_enrichment = "ANY",
    dimension_reduction = 'ANY',
    nn_network = "ANY",
    images = "ANY",
    parameters = "ANY",
    instructions = "ANY",
    offset_file = "ANY",
    OS_platform = "ANY"
  ),

  prototype = list(
    raw_exprs = NULL,
    norm_expr = NULL,
    norm_scaled_expr = NULL,
    custom_expr = NULL,
    spatial_locs = NULL,
    cell_metadata = NULL,
    gene_metadata = NULL,
    cell_ID = NULL,
    gene_ID = NULL,
    spatial_network = NULL,
    spatial_grid = NULL,
    spatial_enrichment = NULL,
    dimension_reduction = NULL,
    nn_network = NULL,
    images = NULL,
    parameters = NULL,
    instructions = NULL,
    offset_file = NULL,
    OS_platform = NULL
  ),

  validity = function(object) {

    if(any(lapply(list(object@raw_exprs), is.null) == TRUE)) {
      return('expression and spatial locations slots need to be filled in')
    }

    # check validity of instructions vector if provided by the user
    if(!is.null(object@instructions)) {

      instr_names = c("python_path", "save_dir", "plot_format", "dpi", "units", "height", "width")
      missing_names = instr_names[!instr_names %in% names(object@instructions)]
      if(length(missing_names) != 0) {
        return('\t Instruction parameters are missing for: ', missing_names, '\t')
      }
    }

    #if(any(lapply(list(object@raw_exprs, object@spatial_locs), is.null) == TRUE)) {
    #  return('expression and spatial locations slots need to be filled in')
    #}

    #if(ncol(object@raw_exprs) != nrow(object@spatial_locs)) {
    #  return('number of cells do not correspond between expression matrix and spatial locations')
    #}
    return(TRUE)
  }
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
    cat(
      "An object of class",
      class(object),
      "\n",
      nrow(x = object@raw_exprs),
      "genes across",
      ncol(x = object@raw_exprs),
      "samples.\n \n"
    )
    cat('Steps and parameters used: \n \n')
    print(object@parameters)
    invisible(x = NULL)
  }
)





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
                                     is_docker = FALSE) {

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

  instructions_list = list(python_path, show_plot, return_plot,
                           save_plot, save_dir, plot_format, dpi, units, height, width, is_docker)
  names(instructions_list) = c('python_path', 'show_plot', 'return_plot',
                               'save_plot', 'save_dir', 'plot_format', 'dpi', 'units', 'height', 'width', 'is_docker')
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
    spM = t_giotto(spM)
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



#' @title evaluate_spatial_locations
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
    warning('\n spatial locations are not given, dummy 3D data will be created \n')
    spatial_locs = data.table::data.table(sdimx = 1:dummy_n,
                                          sdimy = 1:dummy_n,
                                          sdimz = 1:dummy_n)

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




#' @title create Giotto object
#' @description Function to create a giotto object
#' @param raw_exprs matrix with raw expression counts [required]
#' @param spatial_locs data.table or data.frame with coordinates for cell centroids
#' @param norm_expr normalized expression values
#' @param norm_scaled_expr scaled expression values
#' @param custom_expr custom expression values
#' @param cell_metadata cell annotation metadata
#' @param gene_metadata gene annotation metadata
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
#' @return giotto object
#' @details
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
#'   \item{dimensions reductions}
#'   \item{nearest neighbours networks}
#'   \item{images}
#' }
#'
#' @keywords giotto
#' @importFrom methods new
#' @export
createGiottoObject <- function(raw_exprs,
                               spatial_locs = NULL,
                               norm_expr = NULL,
                               norm_scaled_expr = NULL,
                               custom_expr = NULL,
                               cell_metadata = NULL,
                               gene_metadata = NULL,
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
                               cores = NA) {

  # create minimum giotto
  gobject = giotto(raw_exprs = raw_exprs,
                   spatial_locs = spatial_locs,
                   norm_expr = NULL,
                   norm_scaled_expr = NULL,
                   custom_expr = NULL,
                   cell_metadata = cell_metadata,
                   gene_metadata = gene_metadata,
                   cell_ID = NULL,
                   gene_ID = NULL,
                   spatial_network = NULL,
                   spatial_grid = NULL,
                   spatial_enrichment = NULL,
                   dimension_reduction = NULL,
                   nn_network = NULL,
                   images = NULL,
                   parameters = NULL,
                   offset_file = offset_file,
                   instructions = instructions,
                   OS_platform = .Platform[['OS.type']])


  # data.table: set global variable
  cell_ID = gene_ID = NULL

  # check if all optional packages are installed
  extra_packages = c("scran", "MAST", "smfishHmrf", "trendsceek", "SPARK", "multinet", "RTriangle", "FactoMiner")
  pack_index = extra_packages %in% rownames(utils::installed.packages())
  extra_installed_packages = extra_packages[pack_index]
  extra_not_installed_packages = extra_packages[!pack_index]

  if(any(pack_index == FALSE) == TRUE) {
    cat("Consider to install these (optional) packages to run all possible Giotto commands for spatial analyses: ",
        extra_not_installed_packages)
    cat("\n Giotto does not automatically install all these packages as they are not absolutely required and this reduces the number of dependencies")
  }


  # if cores is not set, then set number of cores automatically, but with limit of 10
  cores = determine_cores(cores)
  data.table::setDTthreads(threads = cores)

  ## read expression matrix
  raw_exprs = evaluate_expr_matrix(raw_exprs, cores = cores, sparse = TRUE)
  gobject@raw_exprs = raw_exprs

  # check rownames and colnames
  if(any(duplicated(rownames(raw_exprs)))) {
    stop("row names contain duplicates, please remove or rename")
  }

  if(any(duplicated(colnames(raw_exprs)))) {
    stop("column names contain duplicates, please remove or rename")
  }

  # prepare other slots related to matrix
  gobject@cell_ID = colnames(raw_exprs)
  gobject@gene_ID = rownames(raw_exprs)
  gobject@parameters = list()


  ## set instructions
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


  ## spatial locations
  dummy_n = ncol(raw_exprs)
  spatial_locs = evaluate_spatial_locations(spatial_locs = spatial_locs,
                                            cores = cores,
                                            dummy_n = dummy_n,
                                            expr_matrix = raw_exprs)


  # check if dimensions agree
  if(nrow(spatial_locs) != ncol(raw_exprs)) {
    stop('\n Number of rows of spatial location must equal number of columns of expression matrix \n')
  }

  ## force dimension names
  spatial_dimensions = c('x', 'y', 'z')
  colnames(spatial_locs) = paste0('sdim', spatial_dimensions[1:ncol(spatial_locs)])

  # add cell_ID column
  spatial_locs[, cell_ID := colnames(raw_exprs)]
  gobject@spatial_locs = spatial_locs
  #gobject@spatial_locs[, cell_ID := colnames(raw_exprs)]




  ## OPTIONAL:
  # add other normalized expression data
  if(!is.null(norm_expr)) {

    norm_expr = evaluate_expr_matrix(norm_expr, cores = cores, sparse = F)

    if(all(dim(norm_expr) == dim(raw_exprs)) &
       all(colnames(norm_expr) == colnames(raw_exprs)) &
       all(rownames(norm_expr) == rownames(raw_exprs))) {

      gobject@norm_expr = as.matrix(norm_expr)
    } else {
      stop('\n dimensions, row or column names are not the same between normalized and raw expression \n')
    }
  }

  # add other normalized and scaled expression data
  if(!is.null(norm_scaled_expr)) {

    norm_scaled_expr = evaluate_expr_matrix(norm_scaled_expr, cores = cores, sparse = F)

    if(all(dim(norm_scaled_expr) == dim(raw_exprs)) &
       all(colnames(norm_scaled_expr) == colnames(raw_exprs)) &
       all(rownames(norm_scaled_expr) == rownames(raw_exprs))) {

      gobject@norm_scaled_expr = as.matrix(norm_scaled_expr)
    } else {
      stop('\n dimensions, row or column names are not the same between normalized + scaled and raw expression \n')
    }
  }

  # add other custom normalized expression data
  if(!is.null(custom_expr)) {

    custom_expr = evaluate_expr_matrix(custom_expr, cores = cores, sparse = F)

    if(all(dim(custom_expr) == dim(raw_exprs)) &
       all(colnames(custom_expr) == colnames(raw_exprs)) &
       all(rownames(custom_expr) == rownames(raw_exprs))) {

      gobject@custom_expr = as.matrix(custom_expr)
    } else {
      stop('\n dimensions, row or column names are not the same between custom normalized and raw expression \n')
    }
  }



  ## cell metadata
  if(is.null(cell_metadata)) {
    gobject@cell_metadata = data.table::data.table(cell_ID = colnames(raw_exprs))
  } else {
    gobject@cell_metadata = data.table::as.data.table(gobject@cell_metadata)
    gobject@cell_metadata[, cell_ID := colnames(raw_exprs)]
    # put cell_ID first
    all_colnames = colnames(gobject@cell_metadata)
    other_colnames = grep('cell_ID', all_colnames, invert = T, value = T)
    gobject@cell_metadata = gobject@cell_metadata[, c('cell_ID', other_colnames), with = FALSE]
  }

  ## gene metadata
  if(is.null(gene_metadata)) {
    gobject@gene_metadata = data.table::data.table(gene_ID = rownames(raw_exprs))
  } else {
    gobject@gene_metadata = data.table::as.data.table(gobject@gene_metadata)
    gobject@gene_metadata[, gene_ID := rownames(raw_exprs)]
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
            spatial_network_Obj = create_spatialNetworkObject(name = networkname,networkDT = network)
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

      if(methods::is(im, 'imageGiottoObj')) {
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






#' @title createGiottoVisiumObject
#' @description creates Giotto object directly from a 10X visium folder
#' @param visium_dir path to the 10X visium directory [required]
#' @param expr_data raw or filtered data (see details)
#' @param gene_column_index which column index to select (see details)
#' @param h5_visium_path path to visium 10X .h5 file
#' @param h5_gene_ids gene names as symbols (default) or ensemble gene ids
#' @param h5_tissue_positions_path path to tissue locations (.csv file)
#' @param h5_image_png_path path to tissue .png file (optional)
#' @param png_name select name of png to use (see details)
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
#' }
#'
#' If starting from a Visium 10X .h5 file
#' \itemize{
#'   \item{h5_visium_path: full path to .h5 file: /your/path/to/visium_file.h5}
#'   \item{h5_tissue_positions_path: full path to spatial locations file: /you/path/to/tissue_positions_list.csv}
#'   \item{h5_image_png_path: full path to png: /your/path/to/images/tissue_lowres_image.png}
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
                                    png_name = NULL,
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

      visium_png = createGiottoImage(gobject = NULL, spatial_locs =  spatial_locs,
                                     mg_object = mg_img, name = 'image',
                                     xmax_adj = xmax_adj, xmin_adj = xmin_adj,
                                     ymax_adj = ymax_adj, ymin_adj = ymin_adj)
      visium_png_list = list(visium_png)
      names(visium_png_list) = c('image')
    } else {
      visium_png_list = NULL
    }

    # create Giotto object
    giotto_object = createGiottoObject(raw_exprs = raw_matrix,
                                       spatial_locs = spatial_locs,
                                       instructions = instructions,
                                       cell_metadata = spatial_results[,.(in_tissue, array_row, array_col)],
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


    visium_png = createGiottoImage(gobject = NULL, spatial_locs =  spatial_locs,
                                   mg_object = mg_img, name = 'image',
                                   xmax_adj = xmax_adj, xmin_adj = xmin_adj,
                                   ymax_adj = ymax_adj, ymin_adj = ymin_adj)

    visium_png_list = list(visium_png)
    names(visium_png_list) = c('image')

    giotto_object = createGiottoObject(raw_exprs = raw_matrix,
                                       spatial_locs = spatial_locs,
                                       instructions = instructions,
                                       cell_metadata = spatial_results[,.(in_tissue, array_row, array_col)],
                                       images = visium_png_list)
    return(giotto_object)

  }

}




