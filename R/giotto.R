
#### Giotto class ####


##### Giotto Class Check Function ####
#TODO
# @title Check Giotto Object
# @name checkGiottoObj
# @description check function for S4 giotto object
# @param object giotto object to check
# @keywords internal
# checkGiottoObj = function(object) {
#
#   # !!! (validity check still under construction) !!!
#   # Add check class definition when finished
#
#   errors = character()
#   ch = box_chars()
#
#   gDataSlots = list(expression = 'expression',
#                     expression_feat = 'expression_feat',
#                     spatial_locs = 'spatial_locs',
#                     spatial_info = 'spatial_info',
#                     feat_info = 'feat_info',
#                     cell_metadata = 'cell_metadata',
#                     feat_metadata = 'feat_metadata',
#                     cell_ID = 'cell_ID',
#                     feat_ID = 'feat_ID',
#                     spatial_network = 'spatial_network',
#                     spatial_grid = 'spatial_grid',
#                     spatial_enrichment = 'spatial_enrichment',
#                     dimension_reduction = 'dimension_reduction',
#                     nn_network = 'nn_network')
#
#
#   if(!is.null(slot(object, 'expression'))) {
#     avail_expr = list_expression(object)
#     for(expr in seq(nrow(avail_expr)))
#     if(!inherits(slot(object, 'expression')[[avail_expr[[expr]]]][[avail_expr[[expr]]]][[avail_expr[[expr]]]], 'data.table')) {
#       print('uh oh')
#     }
#   }
#
#
#
#
#
#   # Check for nesting issues (due to updates) - to be deprecated
#   available_nn = list_nearest_networks(object)
#   available_grids = list_spatial_grids(object)
#
#
#
#
#   ## Match spatial units and feature types ##
#
#   # Find existing spat_units in source data
#   uniqueSrcSU = list_expression(object)$spat_unit
#   uniqueSrcSU = unique(uniqueSrcSU, list_spatial_info_names(object))
#
#   # Find existing feat_types in source data
#   uniqueSrcFT = list_expression(object)$feat_type
#   uniqueSrcFT = unique(uniqueSrcFT, list_feature_info_names(object))
#
#   # Find existing spat_units and feat_types within other data slots
#   uniqueDatSU = lapply(gDataSlots, function(x) list_giotto_data(gobject = object, slot = x)$spat_unit)
#   uniqueDatFT = lapply(gDataSlots, function(x) list_giotto_data(gobject = object, slot = x)$feat_type)
#
#   # # If any spat_units exist, ensure that associated objects have identical cell_IDs
#   # if(length(uniqueSU) > 0) {
#   #   for(SU in seq_along(uniqueSU)) {
#   #
#   #     all(lapply(gDataSlots, function(x) {
#   #       list_giotto_data(gobject = object,
#   #                        slot = x)
#   #     }))
#   #   }
#   # }
#
#
#   # cell_ID slot
#   slot_cell_ID = slot(object, 'cell_ID')
#
#   ## check nesting depth
#   if(depth(slot_cell_ID) != 1) {
#     msg = paste0('Nesting depth of cell_ID is ', depth(slot_cell_ID), '.
#         Should be 1. Ex:
#         \t.
#         \t', ch$b,'spat_unit
#         \t', ch$s, ch$b,'cell IDs\n')
#     errors = c(errors, msg)
#   }
#
#   for(spat_unit in names(slot_cell_ID)) {
#     if(!inherits(slot_cell_ID[[spat_unit]], 'character')) {
#       msg = paste0('cell_ID slot info for spat_unit "', spat_unit,'" is of class "',
#                    class(slot_cell_ID[[spat_unit]]), '".\n    Should be "character"\n')
#       errors = c(errors, msg)
#     }
#   }
#
#
#   if(length(errors) == 0) TRUE else errors
# }




##### Giotto Class Definition ####

#' @title S4 giotto Class
#' @description Framework of giotto object to store and work with spatial expression data
#' @concept giotto object
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
#' @slot largeImages slot to store giottoLargeImage objects
#' @slot parameters slot to save parameters that have been used
#' @slot instructions slot for global function instructions
#' @slot offset_file offset file used to stitch together image fields
#' @slot OS_platform Operating System to run Giotto analysis on
#' @slot join_info information about joined Giotto objects
#' @details
#' [\strong{expression}] There are several ways to provide expression information:
#'
#' [\strong{expression_feat}] The different features or modalities such as rna, protein, metabolites, ...
#' that are provided in the expression slot.
#'
#'
#' @export
giotto <- setClass(
  "giotto",
  slots = c(
    expression = "nullOrList",
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
    largeImages = "ANY",
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
    largeImages = NULL,
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

    for(spat_unit in names(object@expression_feat)) {
      cat('spatial units = ', spat_unit, '\n')
      for(feat_type in unique(object@expression_feat)) {
        cat("features = ", feat_type, "\n")

        cat(
          nrow(x = object@expression[[spat_unit]][[feat_type]][['raw']]),
          "features across",
          ncol(x = object@expression[[spat_unit]][[feat_type]][['raw']]),
          "samples.\n \n"
        )
      }
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
                                  packages_to_install = c('pandas==1.1.5',
                                                          'networkx==2.6.3',
                                                          'python-igraph==0.9.6',
                                                          'leidenalg==0.8.7',
                                                          'python-louvain==0.15',
                                                          'python.app==2',
                                                          'scikit-learn==0.24.2')) {

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

        if(inherits(python_path, 'try-error')) {
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


#' @title Create instructions for giotto functions
#' @name createGiottoInstructions
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
#' @param fiji_path path to fiji executable
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

#' @title Read giotto instructions associated with giotto object
#' @name readGiottoInstructions
#' @description Retrieves the instruction associated with the provided parameter
#' @param giotto_instructions giotto object or result from createGiottoInstructions()
#' @param param parameter to retrieve
#' @return specific parameter
#' @export
readGiottoInstructions <- function(giotto_instructions,
                                   param = NULL) {

  # get instructions if provided the giotto object
  if(inherits(giotto_instructions, 'giotto')) {
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


#' @title Show giotto instructions associated with giotto object
#' @name showGiottoInstructions
#' @description Function to display all instructions from giotto object
#' @param gobject giotto object
#' @return named vector with giotto instructions
#' @export
showGiottoInstructions = function(gobject) {

  instrs = gobject@instructions
  return(instrs)
}


#' @title Change giotto instruction(s) associated with giotto object
#' @name changeGiottoInstructions
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



#' @title Replace all giotto instructions in giotto object
#' @name replaceGiottoInstructions
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




#### aggregateExprObj Class ####

##### aggregateExprObj Class Check Function ####

#' @title Check aggregateExprObj
#' @name checkAggExprObj
#' @description Check function for S4 aggregateExprObj
#' @param object S4 aggregateExprObj to check
#' @keywords internal
checkAggExprObj = function(object) {
  errors = character()

  # Check for expr info
  if(is.null(slot(object, 'expressionMat'))) {
    msg = paste0('No expression matrix found in expressionMat slot.\nThis object contains no expression information.\n')
    errors = c(errors, msg)
  }

  # Check expr matrix class
  if(!inherits(slot(object, 'expressionMat'), c('dgeMatrix', 'dgCMatrix'))) {
    msg = paste0('Expression matrix should be provided as either dgCmatrix (sparse) or dgeMatrix (dense) from the Matrix package.\n')
    errors = c(errors, msg)
  }
  if(inherits(slot(object, 'expressionMat'), 'dgCMatrix') & !isTRUE(slot(object, 'sparse'))) {
    msg = paste0('expressionMat slot contains "dgCMatrix".\nsparse slot should be "TRUE"\n')
    errors = c(errors, msg)
  }
  if(inherits(slot(object, 'expressionMat'), 'dgeMatrix') & isTRUE(slot(object, 'sparse'))) {
    msg = paste0('expressionMat slot contains "dgeMatrix".\nsparse slot should be "FALSE"\n')
    errors = c(errors, msg)
  }

  if(length(errors) == 0) TRUE else errors
}



##### aggregateExprObj Class Definition ####

#' @title S4 aggregateExprObj
#' @description Framework to store aggregated expression information
#' @slot name name of aggregateExprObj
#' @slot sparse (boolean) if matrix is sparse
#' @slot expressionMat matrix of expression information
#' @slot provenance origin data of aggregated expression information (if applicable)
#' @slot spat_unit spatial unit of aggregated expression (e.g. 'cell')
#' @slot feat_type feature type of aggregated expression (e.g. 'rna', 'protein')
#' @slot misc misc
#' @export
setClass('aggregateExprObj',
         slots = c(name = 'nullOrChar',
                   sparse = 'logical',
                   expressionMat = 'ANY',
                   provenance = 'ANY',
                   spat_unit = 'character',
                   feat_type = 'character',
                   misc = 'ANY'),
         prototype = list(name = NULL,
                          sparse = NULL,
                          expressionMat = NULL,
                          provenance = NULL,
                          spat_unit = NULL,
                          feat_type = NULL,
                          misc = NULL),
         validity = checkAggExprObj)

##### aggregateExprObj Show Function ####

#' show method for aggregateExprObj class
#' @param object aggregated expression object
#' @aliases show,aggregateExprObj-method
#' @docType methods
#' @rdname show-methods
setMethod(
  f = "show", signature('aggregateExprObj'), function(object) {

    cat("An object of class",  class(object), "\n\n")
    if(!is.null(slot(object, 'sparse'))) {
      if(isTRUE(slot(object, 'sparse'))) {
        cat('Contains sparse matrix with aggregated expression information\n')
      } else if(!isTRUE(slot(object, 'sparse'))) {
        cat('Contains dense matrix with aggregated expression information\n')
      }
    }

    cat(paste0('for spatial unit: "', slot(object, 'spat_unit'), '" and feature type: "', slot(object, 'feat_type'),'" \n'))
    if(!is.null(slot(object, 'provenance'))) cat('  Provenance: ', unlist(slot(object, 'provenance')),'\n')

    print(slot(object, 'expressionMat'))

    cat('\n')

  })






#### Giotto matrices ####

#' @title Read expression matrix
#' @name readExprMatrix
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



#' @title Evaluate expression matrix
#' @name evaluate_expr_matrix
#' @description Evaluate expression matrices that are provided as input and converts
#' them to preferred format for Giotto object. A filepath can also be provided through
#' \code{inputmatrix} param. If this is done, the function will attempt to read the
#' matrix file in using \code{\link{readExprMatrix}}.
#' @param inputmatrix inputmatrix to evaluate
#' @param sparse create sparse matrix (default = TRUE)
#' @param cores how many cores to use
#' @return sparse matrix
#' @details The inputmatrix can be a matrix, sparse matrix, data.frame, data.table or path to any of these.
#' @keywords internal
evaluate_expr_matrix = function(inputmatrix,
                                sparse = TRUE,
                                cores = NA) {


  if(inherits(inputmatrix, 'character')) {
    mymatrix = readExprMatrix(inputmatrix, cores = cores)
  } else if(inherits(inputmatrix, 'Matrix')) {
    mymatrix = inputmatrix
  } else if(inherits(inputmatrix, 'DelayedMatrix')) {
    mymatrix = inputmatrix
  } else if(inherits(inputmatrix, 'data.table')) {
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
    stop("expression input needs to be a path to matrix-like data or an object of class 'Matrix', 'data.table', 'data.frame' or 'matrix'")
  }


  # check rownames and colnames
  if(any(duplicated(rownames(mymatrix)))) {
    stop("row names contains duplicates, please remove or rename")
  }

  if(any(duplicated(colnames(mymatrix)))) {
    stop("column names contains duplicates, please remove or rename")
  }


  return(mymatrix)
}



#' @title Find depth of subnesting
#' @name depth
#' @param this object to evaluate
#' @param method max (default) or min nesting to detect
#' @param sig signature or class to check for. Default is 'data.frame'
#' @description Determines how many max or min layers of subnesting there is, with the end
#' object (defined by param sig or a list of length 0) being layer 0
#' @details https://stackoverflow.com/questions/13432863/determine-level-of-nesting-in-r
#' @keywords internal
depth <- function(this,
                  method = c('max', 'min'),
                  sig = 'data.frame') {

  method = match.arg(arg = method, choices = c('max', 'min'))

  if(inherits(this, sig)) {
    return(0)
  }
  if(is.list(this) && length(this) == 0) {
    return(0)
  }
  if(method == 'max') {
    ifelse(is.list(this), 1L + max(sapply(this, function(x) depth(x, method = method, sig = sig))), 0L)
  } else if(method == 'min') {
    ifelse(is.list(this), 1L + min(sapply(this, function(x) depth(x, method = method, sig = sig))), 0L)
  }

}




#' @title Read expression data
#' @name read_expression_data
#' @param expr_list (nested) list with expression data
#' @param sparse read matrix data in a sparse manner
#' @param cores number of cores to use
#' @param default_feat_type default feature type if nothing is provided
#' @param verbose be verbose
#' @details
#'
#' mylistA = list('a' = matrix(1:5), 'b' = matrix(1:5))
#' depth(mylistA)
#'
#' mylistB = list(A = list('a' = matrix(1:5), 'b' = matrix(1:5)),
#'                B = list('c' = matrix(1:5),'d' = matrix(1:5)))
#' depth(mylistB)
#'
#' mylistC = list('RNA' = list('RAW' = list('cell' = matrix(1:5), 'nucleus' = matrix(6:10)),
#'                             'NORM' = list('cell' = matrix(11:15),'nucleus' = matrix(20:25))),
#'                'PROT' = list('RAW' = list('cell' = matrix(16:20))))
#' depth(mylistC)
#'
#' mymatD = matrix(data = 1:4)
#'
#' @keywords internal
read_expression_data = function(expr_list = NULL,
                                sparse = TRUE,
                                cores = NA,
                                default_feat_type = NULL,
                                verbose = TRUE) {

  # Check
  if(is.null(expr_list)) return(NULL)

  # import box characters
  ch = box_chars()

  # Set default feature type if missing
  if(is.null(default_feat_type)) default_feat_type = 'rna'

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

  # no expression information
  if(list_depth == 0) {
    stop('Depth of expression list is 0, no expression information is provided \n')
  }

  # too much information
  if(list_depth > 3) {
    stop('Depth of expression list is more than 3, only 3 levels are possible:
       0)', ch$s, '.
       1)', ch$s, ch$b, 'spatial unit (e.g. cell)
       2)', ch$s, ch$s, ch$b, 'feature (e.g. RNA)
       3)', ch$s, ch$s, ch$s, ch$b, 'data type (e.g. raw)\n')
  }


  return_list = list()


  # 2. for list with 1 depth
  if(list_depth == 1) {

    cat('list depth of 1 \n')


    for(data in names(expr_list)) {

      # print(expr_list[[data]])
      res_mat = evaluate_expr_matrix(inputmatrix = expr_list[[data]],
                                     sparse = sparse,
                                     cores = cores)
      # add default feat == 'rna'
      # add default region == 'cell'
      return_list[['cell']][[default_feat_type]][[data]] = res_mat

    }


  } else if(list_depth == 2) {

    cat('list depth of 2 \n')
    # add default region == 'cell'

    for(feat in names(expr_list)) {
      for(data in names(expr_list[[feat]])) {

        res_mat = evaluate_expr_matrix(inputmatrix = expr_list[[feat]][[data]],
                                       sparse = sparse,
                                       cores = cores)
        # add default region == 'cell'
        return_list[['cell']][[feat]][[data]] = res_mat

      }
    }

  } else if(list_depth == 3) {

    cat('list depth of 3 \n')

    for(region in names(expr_list)) {
      for(feat in names(expr_list[[region]])) {
        for(data in names(expr_list[[region]][[feat]])) {

          res_mat = evaluate_expr_matrix(inputmatrix = expr_list[[region]][[feat]][[data]],
                                         sparse = sparse,
                                         cores = cores)
          # add default region == 'cell'
          return_list[[region]][[feat]][[data]] = res_mat

        }
      }
    }



  } else {
    stop('unexpected list_depth error')
  }

  return(return_list)

}



#' @title Set cell and feature IDs
#' @name set_cell_and_feat_IDs
#' @description sets cell and feature IDs based on provided expression data
#' @keywords internal
set_cell_and_feat_IDs = function(gobject) {



  # 1. set cell_ID for each region
  # each regions can have multiple features, but the cell_IDs (spatial units) should be the same
  for(spat_unit in names(gobject@expression)) {
    gobject@cell_ID[[spat_unit]] =  colnames(gobject@expression[[spat_unit]][[1]][[1]])
  }

  # 2. ensure cell_ID and colnames for each matrix are the same
  for(spat_unit in names(gobject@expression)) {

    for(feat in names(gobject@expression[[spat_unit]])) {
      for(data in names(gobject@expression[[spat_unit]][[feat]])) {
        colnames_matrix = colnames(gobject@expression[[spat_unit]][[feat]][[data]])
        if(!identical(colnames_matrix, gobject@cell_ID[[spat_unit]])) {
          stop('Colnames are not the same for feat: ', feat,', spatial unit: ', spat_unit ,', and data: ', data)
        }
      }
    }
  }


  # 3. set feat_ID for each feature
  for(spat_unit in names(gobject@expression)) {
    feat_types_covered = list()
    for(feat_type in names(gobject@expression[[spat_unit]])) {
      if(!feat_type %in% feat_types_covered) {
        gobject@feat_ID[[feat_type]] =  rownames(gobject@expression[[spat_unit]][[feat_type]][[1]])
        feat_types_covered[[feat_type]] = feat_type
      }
    }
  }

  return(gobject)

}




#### Giotto metadata ####

#' @title Set cell metadata
#' @name set_cell_metadata
#' @description Set cell metadata
#' @param gobject giotto object
#' @param cell_metadata cell_metadata supplied as a nested list with spat_unit$feat_type$name
#' @keywords internal
set_cell_metadata = function(gobject,
                             cell_metadata) {

  # if metadata is not provided, then:
  # create metadata for each spatial unit and feature type combination

  # define for data.table :=
  cell_ID = NULL

  if(is.null(cell_metadata)) {

    for(spat_unit in names(gobject@expression)) {

      for(feat_type in names(gobject@expression[[spat_unit]])) {

        if(is.null(gobject@spatial_info)) {
          gobject@cell_metadata[[spat_unit]][[feat_type]] = data.table::data.table(cell_ID = gobject@cell_ID[[spat_unit]])
        } else {
          for(poly in names(gobject@spatial_info)) {
            gobject@cell_metadata[[poly]][[feat_type]] = data.table::data.table(cell_ID = gobject@cell_ID[[poly]])
          }
        }
      }
    }

  } else {

    # extract all metadata information
    # need to be nested list (feature type and spatial unit)
    for(spat_unit in names(cell_metadata)) {

      for(feat_type in names(cell_metadata[[spat_unit]])) {


        gobject@cell_metadata[[spat_unit]][[feat_type]] = data.table::as.data.table(cell_metadata[[spat_unit]][[feat_type]])
        gobject@cell_metadata[[spat_unit]][[feat_type]][, cell_ID := gobject@cell_ID[[spat_unit]]]

        # put cell_ID first
        all_colnames = colnames(gobject@cell_metadata[[spat_unit]][[feat_type]])
        other_colnames = grep('cell_ID', all_colnames, invert = T, value = T)
        gobject@cell_metadata[[spat_unit]][[feat_type]] = gobject@cell_metadata[[spat_unit]][[feat_type]][, c('cell_ID', other_colnames), with = FALSE]

      }
    }
  }

  return(gobject)

}


#' @title Set feature metadata
#' @name set_feature_metadata
#' @description Set feature metadata
#' @keywords internal
set_feature_metadata = function(gobject,
                                feat_metadata) {

  # define for data.table :=
  feat_ID = NULL

  if(is.null(feat_metadata)) {

    for(spat_unit in names(gobject@expression)) {
      for(feat_type in names(gobject@expression[[spat_unit]])) {
        gobject@feat_metadata[[spat_unit]][[feat_type]] = data.table::data.table(feat_ID = gobject@feat_ID[[feat_type]])
      }
    }

  } else {

    for(spat_unit in names(gobject@expression)) {
      for(feat_type in names(gobject@expression[[spat_unit]])) {
        gobject@feat_metadata[[spat_unit]][[feat_type]] = data.table::as.data.table(feat_metadata[[spat_unit]][[feat_type]])
        gobject@feat_metadata[[spat_unit]][[feat_type]][, feat_ID := gobject@feat_ID[[feat_type]]]

        # put feat_ID first
        all_colnames = colnames(gobject@feat_metadata[[spat_unit]][[feat_type]])
        other_colnames = grep('feat_ID', all_colnames, invert = T, value = T)
        gobject@feat_metadata[[spat_unit]][[feat_type]] = gobject@feat_metadata[[spat_unit]][[feat_type]][, c('feat_ID', other_colnames), with = FALSE]

      }
    }
  }

  return(gobject)

}





#### Giotto locations ####

#' @title evaluate_spatial_locations_OLD
#' @name evaluate_spatial_locations_OLD
#' @description Evaluate spatial location input
#' @param spatial_locs spatial locations to evaluate
#' @param cores how many cores to use
#' @param dummy_n number of rows to create dummy spaial locations
#' @param expr_matrix expression matrix to compare the cell IDs with
#' @return data.table
#' @keywords internal
evaluate_spatial_locations_OLD = function(spatial_locs,
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
              '\n Other non numeric columns will be removed\n\n')

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


#' @title Evaluate spatial locations
#' @name evaluate_spatial_locations
#' @description Evaluate spatial location input
#' @param spatial_locs spatial locations to evaluate
#' @param cores how many cores to use
#' @return data.table
#' @keywords internal
evaluate_spatial_locations = function(spatial_locs,
                                      cores = 1) {

  # data.table variables
  cell_ID = NULL

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


  potential_cell_IDs = NULL

  if(length(non_numeric_classes) > 0) {

    non_numeric_indices = which(!column_classes %in% c('numeric','integer'))

    message('There are non numeric or integer columns for the spatial location input at column position(s): ', non_numeric_indices,
            '\n The first non-numeric column will be considered as a cell ID to test for consistency with the expression matrix',
            '\n Other non numeric columns will be removed')


    potential_cell_IDs = spatial_locs[[names(non_numeric_classes)[[1]]]]

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

  # for spatial dimension names
  spatial_dimensions = c('x', 'y', 'z')
  colnames(spatial_locs) = paste0('sdim', spatial_dimensions[1:ncol(spatial_locs)])

  if(!is.null(potential_cell_IDs)) {
    spatial_locs[, cell_ID := potential_cell_IDs]
  }

  return(spatial_locs)

}


#' @title Read spatial location data
#' @name read_spatial_location_data
#' @description read spatial locations
#' @param gobject giotto object
#' @param spat_loc_list list of spatial locations
#' @param cores how many cores to use
#' @param verbose be verbose
#' @return updated giotto object
#' @keywords internal
read_spatial_location_data = function(gobject,
                                      spat_loc_list,
                                      cores = 1,
                                      verbose = TRUE) {

  if(is.null(spat_loc_list)) return(NULL)


  ## to make it compatible with previous version

  # single matrix
  if(inherits(spat_loc_list, c('data.frame', 'data.table', 'matrix', 'character'))) {
    spat_loc_list = list('raw' = spat_loc_list)
  }

  # single path to matrix
  if(length(spat_loc_list) == 1 & !is.list(spat_loc_list)) {
    spat_loc_list = list('raw' = spat_loc_list)
  }



  # 1. get depth of list
  if(verbose == TRUE) print(str(spat_loc_list))
  list_depth = depth(spat_loc_list)

  # no expression information
  if(list_depth == 0) {
    stop('Depth of spatial location list is 0, no expression information is provided \n')
  }

  # too much information
  if(list_depth > 2) {
    stop('Depth of spatial location list is more than 2, only 2 levels are possible:
         1) spatial unit (e.g. cell) --> 2) coordinate (e.g. raw) \n')
  }


  return_list = list()


  # 2. for list with 1 depth
  if(list_depth == 1) {

    cat('list depth of 1 \n')

    for(coord in names(spat_loc_list)) {

      res_spatlocs = evaluate_spatial_locations(spatial_locs = spat_loc_list[[coord]],
                                                cores = cores)

      # add default region == 'cell'
      return_list[['cell']][[coord]] = res_spatlocs

    }


  } else if(list_depth == 2) {

    cat('list depth of 2 \n')
    # add default region == 'cell'

    for(spat_unit in names(spat_loc_list)) {
      for(coord in names(spat_loc_list[[spat_unit]])) {

        res_spatlocs = evaluate_spatial_locations(spatial_locs = spat_loc_list[[spat_unit]][[coord]],
                                                     cores = cores)

        # add default region == 'cell'
        return_list[[spat_unit]][[coord]] = res_spatlocs

      }
    }

  } else {
    stop('unexpected list_depth error')
  }

  return(return_list)

}


#' @title Check spatial location data
#' @name check_spatial_location_data
#' @description check cell ID (spatial unit) names between spatial location and expression data
#' @keywords internal
check_spatial_location_data = function(gobject) {


  for(spat_unit in names(gobject@spatial_locs)) {

    expected_cell_ID_names = gobject@cell_ID[[spat_unit]]

    for(coord in names(gobject@spatial_locs[[spat_unit]])) {

      # 1. get colnames
      spatial_colnames = colnames(gobject@spatial_locs[[spat_unit]][[coord]])

      # if cell_ID column is provided then compare with expected cell_IDs
      if('cell_ID' %in% spatial_colnames) {

        spatial_cell_id_names = gobject@spatial_locs[[spat_unit]][[coord]][['cell_ID']]

        if(!identical(spatial_cell_id_names, expected_cell_ID_names)) {
          message('spatloc cell_IDs: ')
          cat('  ', head(spatial_cell_id_names,3), '...', tail(spatial_cell_id_names,3), '\n')
          message('expression cell_IDs: ')
          cat('  ', head(expected_cell_ID_names,3), '...', tail(expected_cell_ID_names,3), '\n')

          stop('cell_IDs between spatial and expression information are not the same for: \n
                 spatial unit: ', spat_unit, ' and coordinates: ', coord, ' \n')
        }

      } else {

        # if cell_ID column is not provided then add expected cell_IDs

        if(nrow(gobject@spatial_locs[[spat_unit]][[coord]]) != length(expected_cell_ID_names)) {
          stop('Number of rows of spatial locations do not match with cell IDs for: \n
                 spatial unit: ', spat_unit, ' and coordinates: ', coord, ' \n')
        }

        gobject@spatial_locs[[spat_unit]][[coord]][['cell_ID']] = expected_cell_ID_names

      }

    }
  }

  return(gobject)

}


#  TODO Following 3 functions need wrapper function for working with spatlocs and spat_info
#  Additionally add way to refer to subsets of spatial locations by list_ID


#' @title Scale spatial locations
#' @name scale_spatial_locations
#' @description Simple scaling of spatial locations by given \code{scale_factor}.
#' Values will be scaled from the coordinate origin or coordinates provided through
#' \code{scenter} param.
#' @param spatlocs spatial locations information to scale
#' @param scale_factor scaling factor to apply to coordinates.
#' @param scenter center from which to scale spatial coordinates. Given as vector
#' of xy(z) coordinates.
#' @details \code{scale_factor} either given as a single value where it will be applied to
#' x, y, and z (if available) dimensions or as a vector of named values for 'x',
#' 'y', (and 'z').
#' @keywords internal
scale_spatial_locations = function(spatlocs,
                                   scale_factor = c(x=1,y=1,z=1),
                                   scenter = c(x=0,y=0,z=0)) {

  hasZ = 'sdimz' %in% names(spatlocs)

  if(length(scale_factor) == 1) scale_factor = c(x = scale_factor, y = scale_factor, z = scale_factor)
  if(!all(names(scenter) %in% c('x','y','z'))) stop('scenter value names not recognized')
  if(!all(names(scale_factor) %in% c('x','y','z'))) stop('scale_factor value names not recognized')

  # Adjust for scaling center
  spatlocs$sdimx = spatlocs$sdimx - scenter[['x']]
  spatlocs$sdimy = spatlocs$sdimy - scenter[['y']]

  # Perform scale
  spatlocs$sdimx = spatlocs$sdimx*scale_factor[['x']]
  spatlocs$sdimy = spatlocs$sdimy*scale_factor[['y']]

  if(isTRUE(hasZ)) {
    # Adjust for scaling z center
    spatlocs$sdimz = spatlocs$sdimz - scenter[['z']]

    # Perform z scale
    spatlocs$sdimz = spatlocs$sdimx*scale_factor[['z']]

    # Revert z scaling center adjustments
    spatlocs$sdimz = spatlocs$sdimz + scenter[['z']]
  }

  # Revert scaling center adjustments
  spatlocs$sdimx = spatlocs$sdimx + scenter[['x']]
  spatlocs$sdimy = spatlocs$sdimy + scenter[['y']]

  return(spatlocs)
}


# TODO deprecate in favor of name change since z axis can be involved$
#' @title xy_translate_spatial_locations
#' @name xy_translate_spatial_locations
#' @description Translate given X Y coordinates by given x and y translation values
#' @param spatlocs spatial locations to use
#' @param xtranslate value to translate coordinates in the positive x direction
#' @param ytranslate value to translate coordinates in the positive y direction
#' @param ztranslate value to translate coordinates in the positive z direction
#' @keywords internal
xy_translate_spatial_locations = function(spatlocs,
                                          xtranslate = 0,
                                          ytranslate = 0,
                                          ztranslate = 0) {

  spatlocs$sdimx = spatlocs$sdimx + xtranslate
  spatlocs$sdimy = spatlocs$sdimy + ytranslate
  if('sdimz' %in% names(spatlocs)) spatlocs$sdimz = spatlocs$sdimz + ztranslate

  return(spatlocs)
}



#' @title Rotate spatial locations
#' @name rotate_spatial_locations
#' @description Rotate given spatlocs by given radians
#' @param spatlocs spatial locations to use
#' @param rotateradians Named vector of radians for rotation along each of the 3 coordinate
#' axes. If only a single value is provided, it will be treated as xy rotation.
#' @param rcenter center of rotation given as vector xy(z) coordinates (defaults to coordinate center)
#' @details Radians are provided through \code{rotateradians} param as a named vector
#' with values for \code{xy} (yaw), \code{zy} (pitch), \code{xz} (roll)
#' @keywords internal
rotate_spatial_locations = function(spatlocs,
                                    rotateradians = c(xy=0,zy=0,xz=0),
                                    rcenter = c(x=0,y=0,z=0)) {

  if(length(rotateradians) == 1) rotateradians = c(xy=rotateradians,zy=0,xz=0)
  if(!all(names(rotateradians) %in% c('xy','zy','xz'))) stop('rotateradians value names not recognized')
  if(!all(names(rcenter) %in% c('x','y','z'))) stop('rcenter value names not recognized')
  hasZ = 'sdimz' %in% names(spatlocs)

  # xy center of rotation adjustment
  spatlocs$sdimx = spatlocs$sdimx - rcenter[['x']]
  spatlocs$sdimy = spatlocs$sdimy - rcenter[['y']]

  xvals = spatlocs$sdimx
  yvals = spatlocs$sdimy

  # Perform rotation XY
  if(rotateradians[['xy']] != 0) {
    spatlocs$sdimx = xvals*cos(rotateradians[['xy']]) + yvals*sin(rotateradians[['xy']])
    spatlocs$sdimy = -xvals*sin(rotateradians[['xy']]) + yvals*cos(rotateradians[['xy']])
  }

  # if z values are available
  if(isTRUE(hasZ)) {
    # z center of rotation adjustment
    spatlocs$sdimz = spatlocs$sdimz - rcenter[['z']]

    zvals = spatlocs$sdimz

    # Perform rotations
    if(rotateradians[['zy']] != 0) {
      spatlocs$sdimz = zvals*cos(rotateradians[['zy']]) + yvals*sin(rotateradians[['zy']])
      spatlocs$sdimy = -zvals*sin(rotateradians[['zy']]) + yvals*cos(rotateradians[['zy']])
    }

    if(rotateradians[['xz']] != 0) {
      spatlocs$sdimx = xvals*cos(rotateradians[['xz']]) + zvals*sin(rotateradians[['xz']])
      spatlocs$sdimz = -xvals*sin(rotateradians[['xz']]) + zvals*cos(rotateradians[['xz']])
    }

    # Revert z center of rotation adjustment
    spatlocs$sdimz = spatlocs$sdimz + rcenter[['z']]
  }

  # Revert xy center of rotation adjustment
  spatlocs$sdimx = spatlocs$sdimx + rcenter[['x']]
  spatlocs$sdimy = spatlocs$sdimy + rcenter[['y']]

  return(spatlocs)
}




#### Giotto spatial network ####

#' @title Read spatial networks
#' @name read_spatial_networks
#' @description read spatial networks from list
#' @keywords internal
read_spatial_networks = function(gobject,
                                 spatial_network) {

  if(is.null(spatial_network)) {
    message('No spatial networks are provided')
    return(gobject)

  } else {

    for(spat_unit in names(spatial_network)) {
      for(name in names(spatial_network[[spat_unit]])) {

        # first check if corresponding expression matrix exists
        if(!is.null(slot(gobject, 'expression')[[spat_unit]])) {


          # TODO: use fread if it's an existing path

          network = spatial_network[[spat_unit]][[name]]

          if(any(c('data.frame', 'data.table') %in% class(network))) {
            if(all(c('to', 'from', 'weight', 'sdimx_begin', 'sdimy_begin', 'sdimx_end', 'sdimy_end') %in% colnames(network))) {
              spatial_network_Obj = new('spatialNetworkObj',
                                        name = name,
                                        networkDT = network)
              gobject = set_spatialNetwork(gobject = gobject,
                                           spat_unit = spat_unit,
                                           name = name,
                                           spatial_network = spatial_network_Obj)
            } else {
              warning('\n spatial unit: ', spat_unit, ' with network name: ', name, ' does not have all necessary column names, see details
                      and will not be added to the Giotto object \n')
            }
          } else {
            warning('\n spatial unit: ', spat_unit, ' with network name: ', name, ' is not a data.frame or data.table
                    and will not be added to the Giotto object \n')
          }

        }
      }
    }
  }

  return(gobject)

}


#### Giotto spatial enrichment ####


#' @title Read spatial enrichment
#' @name read_spatial_enrichment
#' @description read spatial enrichment results from list
#' @keywords internal
read_spatial_enrichment = function(gobject,
                                   spatial_enrichment) {

  if(is.null(spatial_enrichment)) {
    message('No spatial enrichment results are provided')
    return(gobject)

  } else {

    for(spat_unit in names(spatial_enrichment)) {
      for(name in names(spatial_enrichment[[spat_unit]])) {

        # first check if corresponding expression matrix exists
        if(!is.null(gobject@expression[[spat_unit]])) {


          # TODO: use fread if it's an existing path

          spat_enrich = spatial_enrichment[[spat_unit]][[name]]

          if(nrow(spat_enrich) != ncol(gobject@expression[[spat_unit]][[1]][[1]])) {
            stop('\n spatial enrichment for spatial unit: ', spat_unit, ' and name: ', name, ' does not match the corresponding expression data \n')
          } else {
            gobject@spatial_enrichment[[spat_unit]][[name]] = spat_enrich
          }

        }

      }

    }

  }

  return(gobject)

}



#### Giotto dimension reduction ####


#' @title Read dimensional reduction
#' @name read_dimension_reduction
#' @description read dimension reduction results from list
#' @keywords internal
read_dimension_reduction <- function(gobject,
                                     dimension_reduction) {


  if(is.null(dimension_reduction)) {
    message('No dimension reduction results are provided')
    return(gobject)

  } else {

    for(dim_i in 1:length(dimension_reduction)) {

      dim_red = dimension_reduction[[dim_i]]

      if(all(c('type', 'spat_unit', 'name', 'reduction_method', 'coordinates', 'misc', 'feat_type') %in% names(dim_red))) {

        coord_data = dim_red[['coordinates']]
        spat_unit = dim_red[['spat_unit']]

        if(all(rownames(coord_data) %in% gobject@cell_ID[[spat_unit]])) {

          type_value = dim_red[['type']] # cells or genes
          reduction_meth_value = dim_red[['reduction_method']] # e.g. umap, tsne, ...
          name_value = dim_red[['name']]  # uniq name
          misc_value = dim_red[['misc']]  # additional data
          feat_type = dim_red[['feat_type']]

          gobject@dimension_reduction[[type_value]][[spat_unit]][[feat_type]][[reduction_meth_value]][[name_value]] = dim_red[c('name', 'reduction_method', 'coordinates', 'misc')]
        } else {
          stop('\n rownames for coordinates are not found in gobject IDs \n')
        }

      } else {
        stop('\n each dimension reduction list must contain all required slots, see details. \n')
      }

    }

  }

  return(gobject)

}


#### Giotto nearest network ####

#' @title Read nearest neighbor networks
#' @name read_nearest_networks
#' @description read nearest network results from list
#' @keywords internal
read_nearest_networks = function(gobject,
                                 nn_network) {

  if(is.null(nn_network)) {
    message('No nearest network results are provided')
    return(gobject)

  } else {

    for(nn_i in 1:length(nn_network)) {

      nn_netw = nn_network[[nn_i]]

      if(all(c('spat_unit', 'type', 'name', 'igraph') %in% names(nn_netw))) {

        igraph_data = nn_netw[['igraph']]
        spat_unit = nn_netw[['spat_unit']]

        if(all(names(igraph::V(igraph_data)) %in% gobject@cell_ID[[spat_unit]])) {

          type_value = nn_netw[['type']] # sNN or kNN
          name_value = nn_netw[['name']]  # uniq name

          gobject@nn_network[[spat_unit]][[type_value]][[name_value]][['igraph']] = igraph_data
        } else {
          stop('\n igraph vertex names are not found in gobject IDs \n')
        }

      } else {
        stop('\n each nn network list must contain all required slots, see details. \n')
      }

    }

  }

  return(gobject)

}








#### Giotto spatial info ####


#' @title Evaluate spatial info
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


#' @title Evaluate feature info
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

#' @title Create a giotto object
#' @name createGiottoObject
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
#' @param spatial_grid list of spatial grid(s)
#' @param spatial_grid_name list of spatial grid name(s)
#' @param spatial_enrichment list of spatial enrichment score(s) for each spatial region
#' @param dimension_reduction list of dimension reduction(s)
#' @param nn_network list of nearest neighbor network(s)
#' @param images list of images
#' @param largeImages list of large images
#' @param offset_file file used to stitch fields together (optional)
#' @param instructions list of instructions or output result from \code{\link{createGiottoInstructions}}
#' @param cores how many cores or threads to use to read data if paths are provided
#' @param verbose be verbose when building Giotto object
#' @return giotto object
#' @details
#'
#' See \url{http://giottosuite.com/articles/getting_started_gobject.html} for more details
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
#' This data can also be added afterwards using the \code{\link{addFeatMetadata}} or \code{\link{addCellMetadata}} functions.
#'
#' [\strong{Other information}] Additional information can be provided through the appropriate slots:
#' \itemize{
#'   \item{spatial networks}
#'   \item{spatial grids}
#'   \item{spatial enrichments}
#'   \item{dimensions reduction}
#'   \item{nearest neighbours networks}
#'   \item{images}
#' }
#'
#' @concept giotto
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
                               spatial_grid = NULL,
                               spatial_grid_name = NULL,
                               spatial_enrichment = NULL,
                               dimension_reduction = NULL,
                               nn_network = NULL,
                               images = NULL,
                               largeImages = NULL,
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
                   largeImages = NULL,
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
  extra_packages = c("scran", "MAST", "png", "tiff", "biomaRt", "trendsceek", "multinet", "RTriangle", "FactoMineR")

  pack_index = extra_packages %in% rownames(utils::installed.packages())
  extra_installed_packages = extra_packages[pack_index]
  extra_not_installed_packages = extra_packages[!pack_index]

  if(any(pack_index == FALSE) == TRUE) {
    wrap_msg("Consider to install these (optional) packages to run all possible Giotto commands for spatial analyses: ", extra_not_installed_packages)
    wrap_msg(" Giotto does not automatically install all these packages as they are not absolutely required and this reduces the number of dependencies")
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

    expression_data = read_expression_data(expr_list = expression,
                                           sparse = TRUE,
                                           cores = cores,
                                           default_feat_type = expression_feat,
                                           verbose = verbose)
    gobject@expression = expression_data

    gobject = set_cell_and_feat_IDs(gobject)

  }

  if(verbose) message('finished expression data')


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
  raw_cell_dim_list = list()

  for(spat_unit in names(gobject@expression)) {
    for(feat_type in names(gobject@expression[[spat_unit]])) {
      raw_cell_dim_list[[spat_unit]][[feat_type]] = ncol(gobject@expression[[spat_unit]][[feat_type]][[1]])
    }
  }

  #raw_cell_dim = ncol(gobject@expression[[1]][[1]][[1]]) # number of columns

  # list of spatial location data.table, each with a unique name
  # the default name = 'raw' and correspond to the real physical coordinates
  # additional spatial locations can be provided


  if(!is.null(spatial_locs)) {

    spatial_location_data = read_spatial_location_data(spat_loc_list = spatial_locs,
                                                       cores = cores,
                                                       verbose = verbose)
    gobject@spatial_locs = spatial_location_data

    # TODO: ensure spatial locations and expression matrices have the same cell IDs
    # TODO: give cell IDs if not provided

    gobject =  check_spatial_location_data(gobject)

  } else {

    if(verbose == TRUE) warning('\n No spatial locations have been provided, dummy locations will be created \n')

    # for each spatial unit create a dummy raw spatial location matrix

    for(spat_unit in names(raw_cell_dim_list)) {

      # create square dummy coordinates
      nr_cells = raw_cell_dim_list[[spat_unit]][[1]]
      x = ceiling(sqrt(nr_cells))
      first_col  = rep(1:x, each = x)[1:nr_cells]
      second_col = rep(1:x, times = x)[1:nr_cells]

      spatial_locs = data.table::data.table(cell_ID = gobject@cell_ID[[spat_unit]],
                                            sdimx = first_col,
                                            sdimy = second_col)

      gobject@spatial_locs[[spat_unit]][['raw']] = spatial_locs

    }
  }

  if(verbose) message('finished spatial location data')


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
  gobject = set_cell_metadata(gobject = gobject,
                              cell_metadata = cell_metadata)

  if(verbose) message('finished cell metadata')

  ## feat metadata ##
  ## ------------- ##
  gobject = set_feature_metadata(gobject = gobject,
                                 feat_metadata = feat_metadata)



  ## feature info ##
  ## ------------ ##
  ## place to store individual feature info
  if(is.null(feat_info)) {
    gobject@feat_info = NULL
  } else {
    gobject = addGiottoPoints(gobject = gobject,
                              gpoints = feat_info)
  }




  ### OPTIONAL:
  ## spatial network
  gobject = read_spatial_networks(gobject = gobject,
                                  spatial_network = spatial_network)



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
            if(!inherits(grid, 'data.table')) data.table::as.data.table(grid)
            grid = new('spatialGridObj',
                       name = gridname,
                       gridDT = grid)
            # TODO Assign grid as the first spat_unit and feat_type. Assigment process will need to be improved later
            avail_spat_feats = list_expression(gobject)
            gobject = set_spatialGrid(gobject = gobject,
                                      spat_unit = avail_spat_feats$spat_unit[[1]],
                                      feat_type = avail_spat_feats$feat_type[[1]],
                                      name = gridname,
                                      spatial_grid = grid)
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
  gobject = read_spatial_enrichment(gobject = gobject,
                                    spatial_enrichment = spatial_enrichment)



  ## dimension reduction
  gobject = read_dimension_reduction(gobject = gobject,
                                     dimension_reduction = dimension_reduction)


  # NN network
  gobject = read_nearest_networks(gobject = gobject,
                                  nn_network = nn_network)


  ## images ##
  # expect a list of giotto object images
  # prefer to make giottoImage creation separate from this function
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




#' @title Create a giotto object from 10x visium data
#' @name createGiottoVisiumObject
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

  # data.table: set global variable
  barcode = row_pxl = col_pxl = in_tissue = array_row = array_col = NULL

  if(!is.null(h5_visium_path)) {

    if(verbose) wrap_msg("A path to an .h5 10X file was provided and will be used \n")

    if(!file.exists(h5_visium_path)) stop("The provided path ", h5_visium_path, " does not exist \n")

    # spatial locations
    if(is.null(h5_tissue_positions_path)) stop("A path to the tissue positions (.csv) needs to be provided to h5_tissue_positions_path \n")
    if(!file.exists(h5_tissue_positions_path)) stop("The provided path ", h5_tissue_positions_path, " does not exist \n")

    # get matrix counts
    h5_results = get10Xmatrix_h5(path_to_data = h5_visium_path, gene_ids = h5_gene_ids)
    raw_matrix = h5_results[['Gene Expression']]

    # spatial locations
    spatial_results = data.table::fread(h5_tissue_positions_path)
    colnames(spatial_results) = c('barcode', 'in_tissue', 'array_row', 'array_col', 'col_pxl', 'row_pxl')
    spatial_results = spatial_results[match(colnames(raw_matrix), barcode)]
    spatial_locs = spatial_results[,.(row_pxl,-col_pxl)]
    colnames(spatial_locs) = c('sdimx', 'sdimy')

    # image (optional)
    if(!is.null(h5_image_png_path)) {


      if(!file.exists(h5_image_png_path)) stop("The provided path ", h5_image_png_path,
                                               " does not exist. Set to NULL to exclude or provide the correct path. \n")


      ## check if automatic alignment can be done
      png_name = basename(h5_image_png_path)

      if(png_name == 'tissue_lowres_image.png') {
        if(file.exists(h5_json_scalefactors_path)) {
          if(verbose == TRUE && do_manual_adj == FALSE) wrap_msg('png and scalefactors paths are found and automatic alignment for the lowres image will be attempted\n\n')

          json_info = jsonlite::read_json(h5_json_scalefactors_path)
          scale_factor = json_info[['tissue_lowres_scalef']]

          visium_png = createGiottoImage(gobject = NULL,
                                         spatial_locs = spatial_locs,
                                         mg_object = h5_image_png_path,
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
          if(verbose == TRUE && do_manual_adj == FALSE) wrap_msg('png and scalefactors paths are found and automatic alignment for the hires image will be attempted\n\n')

          json_info = jsonlite::read_json(h5_json_scalefactors_path)
          scale_factor = json_info[['tissue_hires_scalef']]

          visium_png = createGiottoImage(gobject = NULL,
                                         spatial_locs = spatial_locs,
                                         mg_object = h5_image_png_path,
                                         name = 'image',
                                         scale_factor = scale_factor,
                                         do_manual_adj = do_manual_adj,
                                         xmax_adj = xmax_adj,
                                         xmin_adj = xmin_adj,
                                         ymax_adj = ymax_adj,
                                         ymin_adj = ymin_adj)

        }
      } else {
        visium_png = createGiottoImage(gobject = NULL,
                                       spatial_locs =  spatial_locs,
                                       mg_object = h5_image_png_path,
                                       name = 'image',
                                       xmax_adj = xmax_adj,
                                       xmin_adj = xmin_adj,
                                       ymax_adj = ymax_adj,
                                       ymin_adj = ymin_adj)
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
                                       cell_metadata = list('cell' = list('rna' = spatial_results[,.(in_tissue, array_row, array_col)])),
                                       images = visium_png_list)
    return(giotto_object)


  } else {

    if(verbose) message("A structured visium directory will be used \n")

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
    # spatial_results = data.table::fread(paste0(spatial_path, '/','tissue_positions_list.csv'))
    spatial_results = data.table::fread(Sys.glob(paths = file.path(spatial_path, 'tissue_positions*')))
    colnames(spatial_results) = c('barcode', 'in_tissue', 'array_row', 'array_col', 'col_pxl', 'row_pxl')
    spatial_results = spatial_results[match(colnames(raw_matrix), barcode)]
    spatial_locs = spatial_results[,.(row_pxl,-col_pxl)]
    colnames(spatial_locs) = c('sdimx', 'sdimy')

    ## spatial image
    if(is.null(png_name)) {
      png_list = list.files(spatial_path, pattern = "*.png")
      png_name = png_list[1]
    }
    png_path = paste0(spatial_path,'/',png_name)
    if(!file.exists(png_path)) stop(png_path, ' does not exist! \n')





    if(png_name == 'tissue_lowres_image.png') {

      scalefactors_path = paste0(spatial_path,'/','scalefactors_json.json')

      if(file.exists(scalefactors_path)) {
        if(verbose == TRUE && do_manual_adj == FALSE) wrap_msg('png and scalefactors paths are found and automatic alignment for the lowres image will be attempted\n\n')

        json_info = jsonlite::read_json(scalefactors_path)
        scale_factor = json_info[['tissue_lowres_scalef']]

        visium_png = createGiottoImage(gobject = NULL,
                                       spatial_locs = spatial_locs,
                                       mg_object = png_path,
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
        if(verbose == TRUE && do_manual_adj == FALSE) wrap_msg('png and scalefactors paths are found and automatic alignment for the hires image will be attempted\n\n')

        json_info = jsonlite::read_json(scalefactors_path)
        scale_factor = json_info[['tissue_hires_scalef']]

        visium_png = createGiottoImage(gobject = NULL,
                                       spatial_locs = spatial_locs,
                                       mg_object = png_path,
                                       name = 'image',
                                       scale_factor = scale_factor,
                                       do_manual_adj = do_manual_adj,
                                       xmax_adj = xmax_adj,
                                       xmin_adj = xmin_adj,
                                       ymax_adj = ymax_adj,
                                       ymin_adj = ymin_adj)

      }
    } else {
      visium_png = createGiottoImage(gobject = NULL,
                                     spatial_locs =  spatial_locs,
                                     mg_object = png_path,
                                     name = 'image',
                                     xmax_adj = xmax_adj,
                                     xmin_adj = xmin_adj,
                                     ymax_adj = ymax_adj,
                                     ymin_adj = ymin_adj)
    }

    visium_png_list = list(visium_png)
    names(visium_png_list) = c('image')

    giotto_object = createGiottoObject(expression = raw_matrix,
                                       expression_feat = 'rna',
                                       spatial_locs = spatial_locs,
                                       instructions = instructions,
                                       cell_metadata = list('cell' = list('rna' = spatial_results[,.(in_tissue, array_row, array_col)])),
                                       images = visium_png_list)
    return(giotto_object)

  }

}




#' @title Create a giotto object from subcellular data
#' @name createGiottoObjectSubcellular
#' @description Function to create a giotto object starting from subcellular polygon (e.g. cell) and points (e.g. transcripts) information
#' @param gpolygons giotto polygons
#' @param polygon_mask_list_params list parameters for \code{\link{createGiottoPolygonsFromMask}}
#' @param polygon_dfr_list_params list parameters for \code{\link{createGiottoPolygonsFromDfr}}
#' @param gpoints giotto points
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
#' @param largeImages list of large images
#' @param instructions list of instructions or output result from \code{\link{createGiottoInstructions}}
#' @param cores how many cores or threads to use to read data if paths are provided
#' @param verbose be verbose when building Giotto object
#' @return giotto object
#' @details There are two different ways to create a Giotto Object with subcellular information:
#' - Starting from polygons (spatial units e.g. cell) represented by a mask or dataframe file and giotto points (analyte coordinates e.g. transcripts)
#' - Starting from polygons (spatial units e.g. cell) represented by a mask or dataframe file and raw intensity images (e.g. protein stains)
#' @concept giotto
#' @export
createGiottoObjectSubcellular = function(gpolygons = NULL,
                                         polygon_mask_list_params = NULL,
                                         polygon_dfr_list_params = NULL,
                                         gpoints = NULL,
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
                                         largeImages = NULL,
                                         largeImages_list_params = NULL,
                                         instructions = NULL,
                                         cores = NA,
                                         verbose = TRUE) {

  # define for data.table :=
  cell_ID = NULL
  feat_ID = NULL

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
                   largeImages = NULL,
                   parameters = NULL,
                   offset_file = NULL,
                   instructions = instructions,
                   OS_platform = .Platform[['OS.type']])


  ## if cores is not set, then set number of cores automatically, but with limit
  cores = determine_cores(cores)
  data.table::setDTthreads(threads = cores)



  # gpolygons and gpoints need to be provided
  if(is.null(gpolygons)) {
    stop('gpolygons = NULL, spatial polygon information needs to be given (e.g. cell boundary, nucleus, ...)')
  }

  if(is.null(gpoints) & is.null(largeImages)) {
    stop('both gpoints = NULL and largeImages = NULL: \n
         Some kind of feature information needs to be provided (e.g. transcript location or protein intensities)')
  }


  ## extract polygon information ##
  ## --------------------------- ##

  if(is.null(polygon_mask_list_params)) {
    polygon_mask_list_params = list(mask_method = 'guess',
                                    remove_background_polygon = TRUE,
                                    background_algo = c('range'),
                                    fill_holes = TRUE,
                                    poly_IDs = NULL,
                                    flip_vertical = TRUE,
                                    shift_vertical_step = TRUE,
                                    flip_horizontal = TRUE,
                                    shift_horizontal_step = TRUE,
                                    calc_centroids = FALSE,
                                    fix_multipart = TRUE)
  }

  if(is.null(polygon_dfr_list_params)) {
    polygon_dfr_list_params = list(calc_centroids = FALSE)
  }

  if(verbose) cat("1. Start extracting polygon information \n")

  polygon_res = extract_polygon_list(polygonlist = gpolygons,
                                     polygon_mask_list_params = polygon_mask_list_params,
                                     polygon_dfr_list_params = polygon_dfr_list_params)
  gobject@spatial_info = polygon_res

  if(verbose) cat("2. Finished extracting polygon information \n")


  if(verbose) cat("3. Add centroid / spatial locations if available \n")
  for(polygon_info in list_spatial_info_names(gobject)) {

    centroidsDT = gobject@spatial_info[[polygon_info]]@spatVectorCentroids
    if(!is.null(centroidsDT)) {

      if(verbose) cat(" - Add centroid / spatial locations for ", polygon_info, " \n")

      centroidsDT = spatVector_to_dt(centroidsDT)
      centroidsDT_loc = centroidsDT[, .(poly_ID, x, y)]
      colnames(centroidsDT_loc) = c('cell_ID', 'sdimx', 'sdimy')

      gobject = set_spatial_locations(gobject = gobject,
                                      spat_unit = polygon_info,
                                      spat_loc_name = 'raw',
                                      spatlocs = centroidsDT_loc,
                                      verbose = FALSE)
    }

  }
  if(verbose) cat("3. Finish adding centroid / spatial locations \n")

  ## cell ID ##
  ## ------- ##
  for(poly in names(gobject@spatial_info)) {
    unique_poly_names = unique(gobject@spatial_info[[poly]]@spatVector$poly_ID)
    gobject@cell_ID[[poly]] = unique_poly_names
  }

  #gobject@cell_ID = gobject@spatial_info[['cell']]@spatVector$poly_ID
  #gobject@cell_ID = gobject@spatial_info[[1]]@spatVector$poly_ID

  ## extract points information ##
  ## -------------------------- ##


  if(!is.null(gpoints)) {
    if(verbose) cat("3. Start extracting spatial feature information \n")

    points_res = extract_points_list(pointslist = gpoints)
    gobject@feat_info = points_res

    if(verbose) cat("4. Finished extracting spatial feature information \n")

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



  if(!is.null(gpoints)) {

    ## cell metadata ##
    ## ------------- ##
    if(is.null(cell_metadata)) {

      for(feat_type in expression_feat) {
        for(poly in names(gobject@spatial_info)) {
          gobject@cell_metadata[[poly]][[feat_type]] = data.table::data.table(cell_ID = gobject@cell_ID[[poly]])
        }
      }


    } else {

      if(length(cell_metadata) != length(expression_feat)) {
        stop('Number of different molecular features need to correspond with the cell_metadata list length \n')
      }

      for(feat_type in expression_feat) {
        for(poly in names(gobject@spatial_info)) {

          gobject@cell_metadata[[poly]][[feat_type]] = data.table::as.data.table(gobject@cell_metadata[[poly]][[feat_type]])
          gobject@cell_metadata[[poly]][[feat_type]][, cell_ID := gobject@cell_ID[[poly]]]

          # put cell_ID first
          all_colnames = colnames(gobject@cell_metadata[[poly]][[feat_type]])
          other_colnames = grep('cell_ID', all_colnames, invert = T, value = T)
          gobject@cell_metadata[[poly]][[feat_type]] = gobject@cell_metadata[[poly]][[feat_type]][, c('cell_ID', other_colnames), with = FALSE]
        }

      }
    }


    ## feat metadata ##
    ## ------------- ##
    if(is.null(feat_metadata)) {

      for(feat_type in expression_feat) {
        for(poly in names(gobject@spatial_info)) {
          gobject@feat_metadata[[poly]][[feat_type]] = data.table::data.table(feat_ID = gobject@feat_ID[[feat_type]])
        }
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
            spatial_network_Obj = new('spatialNetworkObj',
                                      name = networkname,
                                      networkDT = network)
            gobject = set_spatialNetwork(gobject = gobject,
                                         spat_unit = names(slot(gobject, 'spatial_info'))[[1]],
                                         name = networkname,
                                         spatial_network = spatial_network_Obj)
          } else {
            stop('\n network ', networkname, ' does not have all necessary column names, see details\n')
          }
        } else {
          stop('\n network ', networkname, ' is not a data.frame or data.table\n')
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

        if(inherits(grid, c('data.table', 'data.frame'))) {
          if(!inherits(grid, 'data.table')) grid = as.data.table(grid)
          if(all(c('x_start', 'y_start', 'x_end', 'y_end', 'gr_name') %in% colnames(grid))) {
            # Assume first spat_info and first expression_feat as spat_unit and feat_type respectively
            gobject = set_spatialGrid(gobject = gobject,
                                      spat_unit = names(slot(gobject, 'spatial_info'))[[1]],
                                      feat_type = expression_feat[[1]],
                                      name = gridname,
                                      spatial_grid = grid)
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

  ## largeImages ##
  # expect a list of giotto largeImages or file paths to such images
  if(!is.null(largeImages)) {

    if(verbose) cat("3. Start loading large images \n")



    if(is.null(names(largeImages))) {
      names(largeImages) = paste0('largeImage.', 1:length(largeImages))
    }


    for(largeImage_i in 1:length(largeImages)) {

      im = largeImages[[largeImage_i]]
      im_name = names(largeImages)[[largeImage_i]]

      if(inherits(im, 'giottoLargeImage')) { # giotto largeImage
        gobject@largeImages[[im_name]] = im
      } else if(inherits(im, 'character') & file.exists(im)) { # file path


        if(is.null(largeImages_list_params)) {
          largeImages_list_params = list(negative_y = TRUE,
                                         extent = NULL,
                                         use_rast_ext = FALSE,
                                         image_transformations = NULL,
                                         xmax_bound = NULL,
                                         xmin_bound = NULL,
                                         ymax_bound = NULL,
                                         ymin_bound = NULL,
                                         scale_factor = 1,
                                         verbose = TRUE)
        }

        glargeImage = do.call('createGiottoLargeImage', c(raster_object = im,
                                                          name = im_name,
                                                          largeImages_list_params))

        glargeImage = createGiottoLargeImage(raster_object = im,
                                             negative_y = FALSE,
                                             name = im_name)

        gobject = addGiottoImage(gobject = gobject,
                                 largeImages = list(glargeImage))

      } else {
        warning('large image: ', im, ' is not an existing file path or a giotto largeImage object')
      }

    }

    if(verbose) cat("4. Finished loading large images \n")

  }

  return(gobject)

}




#' @title Create Nanostring CosMx Giotto Object
#' @name createGiottoCosMxObject
#' @description Given the path to a CosMx experiment directory, creates a Giotto
#' object.
#' @param cosmx_dir full path to the exported cosmx directory
#' @param data_to_use which type(s) of expression data to build the gobject with
#' Default is \code{'all'} information available. \code{'subcellular'} loads the transcript
#' coordinates only. \code{'aggregate'} loads the provided aggregated expression matrix.
#' @param FOVs field of views to load (only affects subcellular data and images)
#' @inheritParams createGiottoObjectSubcellular
#' @return a giotto object
#' @export
#' @details
#' [\strong{Expected Directory}] This function generates a giotto object when given a
#' link to a cosmx output directory. It expects the following items within the directory
#' where the \strong{bolded} portions are what this function matches against:
#' \itemize{
#'   \item{\strong{CellComposite} (folder of images)}
#'   \item{\strong{CellLabels} (folder of images)}
#'   \item{\strong{CellOverlay} (folder of images)}
#'   \item{\strong{CompartmentLabels} (folder of images)}
#'   \item{experimentname_\strong{exprMat_file}.csv (file)}
#'   \item{experimentname_\strong{fov_positions_file}.csv (file)}
#'   \item{experimentname_\strong{metadata_file}.csv (file)}
#'   \item{experimentname_\strong{tx_file}.csv (file)}
#' }
#'
#' [\strong{Workflows}] Workflow to use is accessed through the data_to_use param
#' \itemize{
#'   \item{'all' - loads and requires subcellular information from tx_file and fov_positions_file
#'   and also the existing aggregated information (expression, spatial locations, and metadata)
#'   from exprMat_file and metadata_file.}
#'   \item{'subcellular' - loads and requires subcellular information from tx_file and
#'   fov_positions_file only.}
#'   \item{'aggregate' - loads and requires the existing aggregate information (expression,
#'   spatial locations, and metadata) from exprMat_file and metadata_file.}
#' }
#'
#' [\strong{Images}] Images in the default CellComposite, CellLabels, CompartmentLabels, and CellOverlay
#' folders will be loaded as giotto largeImage objects in all workflows as long as they are available.
#' Additionally, CellComposite images will be converted to giotto image objects, making plotting with
#' these image objects more responsive when accessing them from a server.
#' \code{\link{showGiottoImageNames}} can be used to see the available images.
#'
#'
createGiottoCosMxObject = function(cosmx_dir = NULL,
                                   data_to_use = c('all','subcellular','aggregate'),
                                   FOVs = NULL,
                                   instructions = NULL,
                                   cores = NA,
                                   verbose = TRUE) {

  # Define for data.table
  fov = target = x_local_px = y_local_px = z = cell_ID = CenterX_global_px = CenterY_global_px =
    CenterX_local_px = CenterY_local_px = NULL

  # 0. test if folder structure exists and is as expected
  if(is.null(cosmx_dir) | !dir.exists(cosmx_dir)) stop('The full path to a cosmx directory must be given.\n')
  if(isTRUE(verbose)) message('A structured CosMx directory will be used\n')

  # find directories (length = 1 if present, length = 0 if missing)
  dir_items = list(`CellLabels folder` = '*CellLabels',
                   `CompartmentLabels folder` = '*CompartmentLabels',
                   `CellComposite folder` = '*CellComposite',
                   `CellOverlay folder` = '*CellOverlay',
                   `transcript locations file` = '*tx_file*',
                   `fov positions file` = '*fov_positions_file*',
                   `expression matrix file` = '*exprMat_file*',
                   `metadata file` = '*metadata_file*')
  dir_items = lapply(dir_items, function(x) Sys.glob(paths = file.path(cosmx_dir, x)))
  dir_items_lengths = lengths(dir_items)

  if(isTRUE(verbose)) {
    message('Checking directory contents...')
    for(item in names(dir_items)) {
      if(dir_items_lengths[[item]] > 0) {
        message('--| ' ,item, ' found')
      } else {
        warning(item, ' is missing\n')
      }
    }
  }

  # select first directory in list if multiple are detected
  if(any(dir_items_lengths > 1)) {
    warning('Multiple matches for expected subdirectory item(s).\n First matching item selected')

    multiples = which(dir_items_lengths > 1)
    for(mult_i in multiples) {
      message(names(dir_items)[[mult_i]], 'multiple matches found:')
      print(dir_items[[mult_i]])
      dir_items[[mult_i]] = dir_items[[mult_i]][[1]]
    }
  }
  if(isTRUE(verbose)) message('Directory check done')


  # 1. read in data
  # set number of cores automatically, but with limit of 10
  cores = determine_cores(cores)
  data.table::setDTthreads(threads = cores)

  # determine data to use
  data_to_use = match.arg(arg = data_to_use, choices = c('all','subcellular','aggregate'))

  # load in subcellular information, subcellular FOV objects, then join
  if(data_to_use == 'all' | data_to_use == 'subcellular') {

    if(isTRUE(verbose)) message('Loading subcellular information...')

    # subcellular checks
    if(!file.exists(dir_items$`transcript locations file`)) stop('No transcript locations file (.csv) detected')
    if(!file.exists(dir_items$`fov positions file`)) stop('No fov positions file (.csv) detected')

    # FOVs to load
    fov_offset_file = fread(input = dir_items$`fov positions file`, nThread = cores)
    if(is.null(FOVs)) FOVs = fov_offset_file$fov
    FOV_ID = as.list(sprintf('%03d', FOVs))

    #TODO Load only relevant portions of file?

    tx_coord_all = fread(input = dir_items$`transcript locations file`, nThread = cores)
    if(isTRUE(verbose)) message('Subcellular load done')


    fov_gobjects_list = lapply(FOV_ID, function(x) {

      if(isTRUE(verbose)) message('Starting FOV ', x)

      # Build image paths
      if(isTRUE(verbose)) message('Loading image information...')

      composite_dir = Sys.glob(paths = file.path(dir_items$`CellComposite folder`, paste0('*',x, '*')))
      cellLabel_dir = Sys.glob(paths = file.path(dir_items$`CellLabels folder`, paste0('*',x, '*')))
      compartmentLabel_dir = Sys.glob(paths = file.path(dir_items$`CompartmentLabels folder`, paste0('*',x, '*')))
      cellOverlay_dir = Sys.glob(paths = file.path(dir_items$`CellOverlay folder`, paste0('*',x, '*')))
      # Missing warnings
      if(length(composite_dir) == 0) {warning('No composite images found') ; composite_dir = NULL}
      if(length(cellLabel_dir) == 0) {stop('No cell mask images found')} # cell masks are necessary
      if(length(compartmentLabel_dir) == 0) {warning('No compartment label images found') ; compartmentLabel_dir = NULL}
      if(length(cellOverlay_dir) == 0) {warning('No cell polygon overlay images found') ; cellOverlay_dir = NULL}

      if(isTRUE(verbose)) message('Image load done')

      # get FOV specific tx locations
      tx_coord = tx_coord_all[fov == as.numeric(x)]
      tx_coord = tx_coord[,.(target, x_local_px, y_local_px, z)]
      colnames(tx_coord) = c('feat_ID','x','y','z')

      # build giotto object
      if(isTRUE(verbose)) message('Building subcellular giotto object...')
      fov_subset = createGiottoObjectSubcellular(gpoints = list('rna' = tx_coord),
                                                 gpolygons = list('cell' = cellLabel_dir),
                                                 polygon_mask_list_params = list(mask_method = 'guess',
                                                                                 flip_vertical = TRUE,
                                                                                 flip_horizontal = FALSE,
                                                                                 shift_horizontal_step = FALSE),
                                                 instructions = instructions,
                                                 cores = cores)



      # find centroids as spatial locations
      if(isTRUE(verbose)) message('Finding polygon centroids as cell spatial locations...')
      fov_subset = addSpatialCentroidLocations(fov_subset,
                                               poly_info = 'cell',
                                               spat_loc_name = 'raw')



      # create and add giotto image objects
      if(isTRUE(verbose)) message('Attaching image files...')
      print(composite_dir)
      print(cellOverlay_dir)
      print(compartmentLabel_dir)

      gImage_list = list()

      # load image if files are found
      if(!is.null(composite_dir)) gImage_list$composite = createGiottoLargeImage(raster_object = composite_dir, negative_y = F, name = 'composite')
      if(!is.null(cellOverlay_dir)) gImage_list$overlay = createGiottoLargeImage(raster_object = cellOverlay_dir, negative_y = F, name = 'overlay')
      if(!is.null(compartmentLabel_dir)) gImage_list$compartment = createGiottoLargeImage(raster_object = compartmentLabel_dir, negative_y = F, name = 'compartment') #TODO



      if(length(gImage_list) > 0) {
        fov_subset = addGiottoImage(gobject = fov_subset,
                                    largeImages = gImage_list)

        # convert to MG for faster loading (particularly relevant for pulling from server)
        fov_subset = convertGiottoLargeImageToMG(giottoLargeImage = gImage_list$composite, gobject = fov_subset, return_gobject = TRUE)
        # fov_subset = convertGiottoLargeImageToMG(giottoLargeImage = gImage_list$overlay, gobject = fov_subset, return_gobject = TRUE)
        # fov_subset = convertGiottoLargeImageToMG(giottoLargeImage = gImage_list$compartment, gobject = fov_subset, return_gobject = TRUE)
      } else {
        message('No images found for fov')
      }


    }) #lapply end

    # join giotto objects according to FOV positions file
    if(isTRUE(verbose)) message('Joining FOV gobjects...')
    new_gobj_names = paste0('fov', FOV_ID)

    id_match = match(as.numeric(FOV_ID), fov_offset_file$fov)
    x_shifts = fov_offset_file[id_match]$x_global_px
    y_shifts = fov_offset_file[id_match]$y_global_px

    # Join giotto objects
    cosmx_gobject = joinGiottoObjects(gobject_list = fov_gobjects_list,
                                      gobject_names = new_gobj_names,
                                      join_method = 'shift',
                                      x_shift = x_shifts,
                                      y_shift = y_shifts)

  }

  # load in pre-generated aggregated expression matrix
  if(data_to_use == 'aggregate' | data_to_use == 'all') {

    # load aggregate information
    message('Loading provided aggregated information...')

    # aggregate checks
    if(!file.exists(dir_items$`expression matrix file`)) stop('No expression matrix file (.csv) detected')
    if(!file.exists(dir_items$`metadata file`)) stop('No metadata file (.csv) detected. Needed for cell spatial locations.')

    # read in aggregate data
    expr_mat = fread(input = dir_items$`expression matrix file`, nThread = cores)
    metadata = fread(input = dir_items$`metadata file`, nThread = cores)

    # setorder expression and spatlocs
    data.table::setorder(metadata, fov, cell_ID)
    data.table::setorder(expr_mat, fov, cell_ID)


    # generate unique cell IDs
    expr_mat$cell_ID = paste0('fov', sprintf('%03d', expr_mat$fov), '-', 'cell_', expr_mat$cell_ID)
    expr_mat = expr_mat[,-'fov']

    metadata$fov_cell_ID = metadata$cell_ID
    metadata$cell_ID = paste0('fov', sprintf('%03d', metadata$fov), '-', 'cell_', metadata$cell_ID)
    # reorder
    metadata = metadata[,.SD, c('cell_ID','fov','fov_cell_ID')]


    # extract spatial locations
    spatlocs = metadata[,.(CenterX_global_px, CenterY_global_px, cell_ID)]
    spatlocs_fov = metadata[,.(CenterX_local_px, CenterY_local_px, cell_ID)]
    #TODO FOV shifts can be regenerated from this if necessary

    # rename spatloc column names
    spatloc_colnames = c('sdimx', 'sdimy', 'cell_ID')
    colnames(spatlocs) = spatloc_colnames
    colnames(spatlocs_fov) = spatloc_colnames

    # cleanup metadata and spatlocs
    metadata = metadata[,c('CenterX_global_px', 'CenterY_global_px', 'CenterX_local_px', 'CenterY_local_px'):=NULL]
    # find unique cell_IDs present in both expression and metadata
    giotto_cell_ID = unique(intersect(expr_mat$cell_ID, metadata$cell_ID))

    # subset to only unique cell_IDs
    expr_mat = expr_mat[cell_ID %in% giotto_cell_ID,]
    metadata = metadata[cell_ID %in% giotto_cell_ID,]

    # convert expression to sparse matrix
    spM = Matrix::Matrix(as.matrix(expr_mat[,-1]), dimnames = list(expr_mat[[1]], colnames(expr_mat[,-1])), sparse = TRUE)
    spM = t_flex(spM)

    ## Ready for downstream aggregate gobject creation or appending into existing subcellular Giotto object ##



    # create standard gobject from aggregate matrix
    if(data_to_use == 'aggregate') {

      # Create aggregate gobject
      if(isTRUE(verbose)) message('Building giotto object...')
      cosmx_gobject = createGiottoObject(expression = list('raw' = spM),
                                         cell_metadata = list('cell' = list('rna' = metadata)),
                                         spatial_locs = spatlocs,
                                         instructions = instructions,
                                         cores = cores)


      # # load in images if FOV positions are available
      # if(file.exists(dir_items$`fov positions file`)) {
      #   fov_offset_file = fread(input = dir_items$`fov positions file`, nThread = cores)
      #
      #   # load images
      #   if(isTRUE(verbose)) message('Attaching image files...')
      #   composite_dir = Sys.glob(paths = file.path(dir_items$`CellComposite folder`, paste0('/*')))
      #   cellLabel_dir = Sys.glob(paths = file.path(dir_items$`CellLabels folder`, paste0('/*')))
      #   compartmentLabel_dir = Sys.glob(paths = file.path(dir_items$`CompartmentLabels folder`, paste0('/*')))
      #   overlay_dir = Sys.glob(paths = file.path(dir_items$`CellOverlay folder`, paste0('/*')))
      #
      #   if(length(cellLabel_imgList) > 0) cellLabel_imgList = lapply(cellLabel_dir, function(x) {createGiottoLargeImage(x,name = 'cellLabel',negative_y = TRUE)})
      #   if(length(composite_imgList) > 0) composite_imgList = lapply(composite_dir, function(x) {createGiottoLargeImage(x,name = 'composite',negative_y = TRUE)})
      #   if(length(compartmentLabel_dir) > 0) compartmentLabel_imgList = lapply(compartmentLabel_dir, function(x) {createGiottoLargeImage(x,name = 'composite',negative_y = TRUE)})
      #   if(length(overlay_dir) > 0) overlay_imgList = lapply(overlay_dir, function(x) {createGiottoLargeImage(x,name = 'composite',negative_y = TRUE)})
      #
      #
      # } else {warning('No FOV positions file (.csv) detected\n No images or polygons will be loaded.')}

    }

  }

  # add in pre-generated aggregated expression matrix information for 'all' workflow
  if(data_to_use == 'all') {

    # Add aggregate expression information
    if(isTRUE(verbose)) message('Appending provided aggregate expression data as...\n spat_unit: "cell_agg"\n feat_type: "rna"\n name: "raw"')
    # add expression data to expression slot
    cosmx_gobject = set_expression_values(cosmx_gobject,
                                          spat_unit = 'cell_agg',
                                          feat_type = 'rna',
                                          name = 'raw',
                                          values = spM)

    # Add spatial locations
    if(isTRUE(verbose)) message('Appending metadata provided spatial locations data as...\n --> spat_unit: "cell_agg" name: "raw"\n --> spat_unit: "cell" name: "raw_fov"')
    if(isTRUE(verbose)) message('Polygon centroid derived spatial locations assigned as...\n --> spat_unit: "cell" name: "raw" (default)')
    cosmx_gobject = set_spatial_locations(cosmx_gobject,
                                          spat_unit = 'cell_agg',
                                          spat_loc_name = 'raw',
                                          spatlocs = spatlocs)
    cosmx_gobject = set_spatial_locations(cosmx_gobject,
                                          spat_unit = 'cell_agg',
                                          spat_loc_name = 'raw_fov',
                                          spatlocs = spatlocs_fov)

    # initialize cell and feat IDs and metadata slots for 'cell_agg' spat_unit
    cosmx_gobject@cell_ID[['cell_agg']] = colnames(cosmx_gobject@expression[['cell_agg']][[1]][[1]])
    cosmx_gobject@feat_ID[['rna']] = unique(c(cosmx_gobject@feat_ID, rownames(cosmx_gobject@expression[['cell_agg']][['rna']][[1]])))
    cosmx_gobject@cell_metadata[['cell_agg']][['rna']] = data.table::data.table(cell_ID = cosmx_gobject@cell_ID[['cell_agg']])
    cosmx_gobject@feat_metadata[['cell_agg']][['rna']] = data.table::data.table(feat_ID = cosmx_gobject@feat_ID[['rna']])

    # Add metadata to both the given and the poly spat_units
    if(isTRUE(verbose)) message('Appending provided cell metadata...')
    cosmx_gobject = addCellMetadata(cosmx_gobject,
                                    spat_unit = 'cell',
                                    feat_type = 'rna',
                                    new_metadata = metadata,
                                    by_column = TRUE,
                                    column_cell_ID = 'cell_ID')
    cosmx_gobject = addCellMetadata(cosmx_gobject,
                                    spat_unit = 'cell_agg',
                                    feat_type = 'rna',
                                    new_metadata = metadata,
                                    by_column = TRUE,
                                    column_cell_ID = 'cell_ID')
  }

  message('done')
  return(cosmx_gobject)

}




#' @title Create 10x Xenium Giotto Object
#' @name createGiottoXeniumObject
#' @description Given the path to a Xenium experiment output folder, creates a Giotto
#' object
#' @param xenium_dir full path to the exported xenium directory
#' @param data_to_use which type(s) of expression data to build the gobject with
#' (e.g. default: \strong{'subcellular'}, 'aggregate', or 'all')
#' @param load_format files formats from which to load the data. Either `csv` or
#' `parquet` currently supported.
#' @param h5_expression (boolean) whether to load cell_feature_matrix from .h5 file.
#' Default is \code{TRUE}
#' @param h5_gene_ids use gene symbols (default) or ensembl ids for the .h5 gene
#' expression matrix
#' @param bounds_to_load vector of boundary information to load (e.g. \code{'cell'}
#' or \code{'nucleus'} by themselves or \code{c('cell', 'nucleus')} to load both
#' at the same time.)
#' @inheritParams get10Xmatrix
#' @inheritParams createGiottoObjectSubcellular
#' @export
createGiottoXeniumObject = function(xenium_dir,
                                    data_to_use = c('subcellular','aggregate'),
                                    load_format = 'csv',
                                    h5_expression = TRUE,
                                    h5_gene_ids = c('symbols', 'ensembl'),
                                    gene_column_index = 1,
                                    bounds_to_load = c('cell'),
                                    # include_analysis = FALSE,
                                    instructions = NULL,
                                    cores = NA,
                                    verbose = TRUE) {

  # Determine data to load
  data_to_use = match.arg(arg = data_to_use, choices = c('subcellular','aggregate'))

  # Determine load formats
  load_format = 'csv' # TODO Remove this and add as param once other options are available
  load_format = match.arg(arg = load_format, choices = c('csv', 'parquet', 'zarr'))

  # Check needed packages
  if(load_format == 'parquet') {
    package_check(pkg_name = 'arrow', repository = 'CRAN')
    package_check(pkg_name = 'dplyr', repository = 'CRAN')
  }
  if(isTRUE(h5_expression)) {
    package_check(pkg_name = 'hdf5r', repository = 'CRAN')
  }

  ch = box_chars()


  # 0. test if folder structure exists and is as expected


  if(is.null(xenium_dir) | !dir.exists(xenium_dir)) stop('The full path to a xenium directory must be given.\n')
  if(isTRUE(verbose)) message('A structured Xenium directory will be used\n')

  # find items (length = 1 if present, length = 0 if missing)
  dir_items = list(`analysis info` = '*analysis*',
                   `boundary info` = '*bound*',
                   `cell feature matrix` = '*cell_feature_matrix*',
                   `cell metadata` = '*cells*',
                   `image info` = '*tif',
                   `panel metadata` = '*panel*',
                   `raw transcript info` = '*transcripts*',
                   `experiment info (.xenium)` = '*.xenium')

  dir_items = lapply(dir_items, function(x) Sys.glob(paths = file.path(xenium_dir, x)))
  dir_items_lengths = lengths(dir_items)

  if(isTRUE(verbose)) {
    message('Checking directory contents...')
    for(item in names(dir_items)) {
      if(dir_items_lengths[[item]] > 0) {
        message(ch$s, '> ' ,item, ' found')
        for(item_i in seq_along(dir_items[[item]])) { # print found item names
          subItem = gsub(pattern = '.*/', replacement = '', x = dir_items[[item]][[item_i]])
          message(ch$s, ch$s, ch$l,ch$h,ch$h, subItem)
        }
      } else {
        if(data_to_use == 'subcellular') {
          # necessary items
          if(item %in% c('boundary info', 'raw transcript info')) stop(item, ' is missing\n')
          # optional items
          if(item %in% c('image info', 'experiment info (.xenium)', 'panel metadata')) warning(item, ' is missing (optional)\n')
          # items to ignore: analysis info, cell feature matrix, cell metadata
        } else if(data_to_use == 'aggregate') {
          # necessary items
          if(item %in% c('cell feature matrix', 'cell metadata')) stop(item, ' is missing\n')
          # optional items
          if(item %in% c('image info', 'experiment info (.xenium)', 'panel metadata', 'analysis info')) warning(item, ' is missing (optional)\n')
          # items to ignore: boundary info, raw transcript info
        }
      }
    }
  }


  # 1. Select data to load


  # **** transcript info ****
  tx_path = dir_items$`raw transcript info`[grepl(pattern = load_format, dir_items$`raw transcript info`)]
  # **** cell metadata ****
  cell_meta_path = dir_items$`cell metadata`[grepl(pattern = load_format, dir_items$`cell metadata`)]

  # **** boundary info ****
  # Select bound load format
  if(load_format != 'zarr') { # No zarr available for boundary info
    dir_items$`boundary info` = dir_items$`boundary info`[grepl(pattern = load_format, dir_items$`boundary info`)]
  } else dir_items$`boundary info` = dir_items$`boundary info`[grepl(pattern = 'csv', dir_items$`boundary info`)]

  # Organize bound paths by type of bound (bounds_to_load param)
  bound_names = bounds_to_load
  bounds_to_load = as.list(bounds_to_load)
  bound_paths = lapply(bounds_to_load, function(x) dir_items$`boundary info`[grepl(pattern = x, dir_items$`boundary info`)])
  names(bound_paths) = bound_names

  # **** aggregated expression info ****
  if(isTRUE(h5_expression)) { # h5 expression matrix loading is default
    agg_expr_path = dir_items$`cell feature matrix`[grepl(pattern = 'h5', dir_items$`cell feature matrix`)]
    h5_gene_ids = match.arg(arg = h5_gene_ids, choices = c('symbols', 'ensembl'))
  } else if(load_format == 'zarr') {
    agg_expr_path = dir_items$`cell feature matrix`[grepl(pattern = 'zarr', dir_items$`cell feature matrix`)]
  } else { # No parquet for aggregated expression - default to normal 10x loading
    agg_expr_path = dir_items$`cell feature matrix`[sapply(dir_items$`cell feature matrix`, function(x) file_test(op = '-d', x))]
    if(length(agg_expr_path) == 0) stop('Expression matrix cannot be loaded.\nHas cell_feature_matrix(.tar.gz) been unpacked into a directory?')
  }


  if(isTRUE(verbose)) message('Directory check done')


  # 2. read in data


  # set number of cores automatically, but with limit of 10
  cores = determine_cores(cores)
  data.table::setDTthreads(threads = cores)

  if(isTRUE(verbose)) message('Loading feature metadata...')
  feat_meta = data.table::fread(dir_items$`panel metadata`[[1]], nThread = cores)
  colnames(feat_meta)[[1]] = 'feat_ID'

  if(load_format == 'csv') {
    # ---------------------------------------------------------------------------- #
    # **** subcellular info ****
    if(data_to_use == 'subcellular') {
      # append missing QC probe info to feat_meta
      if(isTRUE(h5_expression)) {
        h5 = hdf5r::H5File$new(agg_expr_path)
        tryCatch({
          root = names(h5)
          feature_id = h5[[paste0(root, "/features/id")]][]
          feature_info = h5[[paste0(root,"/features/feature_type")]][]
          feature_names = h5[[paste0(root, "/features/name")]][]
          features_dt = data.table::data.table(
            'id' = feature_id,
            'name' = feature_names,
            'feature_type' = feature_info
          )
        }, finally = {
          h5$close_all()
        })
      } else {
        features_dt = data.table::fread(paste0(agg_expr_path, '/features.tsv.gz'), header = F)
      }
      colnames(features_dt) = c('id', 'feat_ID', 'feat_class')
      feat_meta = merge(features_dt[,c(2,3)], feat_meta, all.x = TRUE, by = 'feat_ID')

      if(isTRUE(verbose)) message('Loading transcript level info...')
      tx_dt = data.table::fread(tx_path[[1]], nThread = cores)
      if(isTRUE(verbose)) message('Loading boundary info...')
      bound_dt_list = lapply(bound_paths, function(x) data.table::fread(x[[1]], nThread = cores))
    }
    # **** aggregate info ****
    if(isTRUE(verbose)) message('Loading cell metadata...')
    cell_meta = data.table::fread(cell_meta_path[[1]], nThread = cores)

    if(data_to_use == 'aggregate') {
      if(isTRUE(verbose)) message('Loading aggregated expression...')
      if(length(agg_expr_path) == 0) stop('Aggregated expression not found.\nPlease confirm h5_expression and load_format params are correct\n')
      if(isTRUE(h5_expression)) agg_expr = get10Xmatrix_h5(path_to_data = agg_expr_path,
                                                           gene_ids = h5_gene_ids,
                                                           remove_zero_rows = TRUE,
                                                           split_by_type = TRUE)
      else agg_expr = get10Xmatrix(path_to_data = agg_expr_path,
                                   gene_column_index = gene_column_index,
                                   remove_zero_rows = TRUE,
                                   split_by_type = TRUE)
    }
    # ---------------------------------------------------------------------------- #
  } else if(load_format == 'parquet') {
    # ---------------------------------------------------------------------------- #
    # **** subcellular info ****
    if(data_to_use == 'subcellular') {
      # append missing QC probe info to feat_meta
      if(isTRUE(h5_expression)) {
        h5 = hdf5r::H5File$new(agg_expr_path)
        tryCatch({
          root = names(h5)
          feature_id = h5[[paste0(root, "/features/id")]][]
          feature_info = h5[[paste0(root,"/features/feature_type")]][]
          feature_names = h5[[paste0(root, "/features/name")]][]
          features_dt = data.table::data.table(
            'id' = feature_id,
            'name' = feature_names,
            'feature_type' = feature_info
          )
        }, finally = {
          h5$close_all()
        })
      } else {
        features_dt = arrow::read_tsv_arrow(paste0(agg_expr_path, '/features.tsv.gz'),
                                            col_names = FALSE) %>%
          data.table::setDT()
      }
      colnames(features_dt) = c('id', 'feat_ID', 'feat_class')
      feat_meta = merge(features_dt[,c(2,3)], feat_meta, all.x = TRUE, by = 'feat_ID')

      if(isTRUE(verbose)) message('Loading transcript level info...')
      tx_dt = arrow::read_parquet(file = tx_path[[1]], as_data_frame = FALSE) %>%
        dplyr::mutate(transcript_id = cast(transcript_id, arrow::string())) %>%
        dplyr::mutate(cell_id = cast(cell_id, arrow::string())) %>%
        dplyr::mutate(feature_name = cast(feature_name, arrow::string())) %>%
        as.data.frame() %>%
        data.table::setDT()
      if(isTRUE(verbose)) message('Loading boundary info...')
      bound_dt_list = lapply(bound_paths, function(x) {
        arrow::read_parquet(file = x[[1]], as_data_frame = FALSE) %>%
          dplyr::mutate(cell_id = cast(cell_id, arrow::string())) %>%
          as.data.frame() %>%
          data.table::setDT()})
    }
    # **** aggregate info ****
    if(data_to_use == 'aggregate') {
      if(isTRUE(verbose)) message('Loading cell metadata...')
      cell_meta = arrow::read_parquet(file = cell_meta_path[[1]], as_data_frame = FALSE) %>%
        dplyr::mutate(cell_id = cast(cell_id, arrow::string())) %>%
        as.data.frame() %>%
        data.table::setDT()

      if(length(agg_expr_path) == 0) stop('Aggregated expression not found.\nPlease confirm h5_expression and load_format params are correct\n')
      # NOTE: no parquet for agg_expr.
      if(isTRUE(h5_expression)) agg_expr = get10Xmatrix_h5(path_to_data = agg_expr_path,
                                                           gene_ids = h5_gene_ids,
                                                           remove_zero_rows = TRUE,
                                                           split_by_type = TRUE)
      else agg_expr = get10Xmatrix(path_to_data = agg_expr_path,
                                   gene_column_index = gene_column_index,
                                   remove_zero_rows = TRUE,
                                   split_by_type = TRUE)
    }
    # ---------------------------------------------------------------------------- #
  } else if(load_format == 'zarr') { # TODO
    # ---------------------------------------------------------------------------- #
    # **** subcellular info ****
    if(data_to_use == 'subcellular') {
      # NOTE: no zarr for boundaries
    }
    # **** aggregate info ****
    if(data_to_use == 'aggregate') {
      if(length(agg_expr_path) == 0) stop('Aggregated expression not found.\nPlease confirm h5_expression and load_format params are correct\n')
    }
    # ---------------------------------------------------------------------------- #
  }

  # TODO load images


  # 3. Create giotto objects

  # define for data.table
  cell_id = cell_ID = x_centroid = y_centroid = .SD = NULL

  if(data_to_use == 'subcellular') {
    if(isTRUE(verbose)) message('Building subcellular giotto object...')
    # Giotto points object
    if(isTRUE(verbose)) message('> points data...')
    gpoints = createGiottoPoints(x = tx_dt[, !'cell_id'], feat_type = 'rna')
    # Giotto polygons object
    if(isTRUE(verbose)) message('> polygons data...')
    gpolys = lapply(names(bound_dt_list),
                    function(bound_type) {
                      if(isTRUE(verbose)) message('  [', bound_type, '] bounds...')
                      bound_dt_list[[bound_type]][, cell_id := as.character(cell_id)]
                      createGiottoPolygonsFromDfr(segmdfr = bound_dt_list[[bound_type]],
                                                  name = bound_type)
                    })

    xenium_gobject = createGiottoObjectSubcellular(gpoints = list(rna = gpoints),
                                                   gpolygons = gpolys,
                                                   instructions = instructions,
                                                   cores = cores,
                                                   verbose = verbose)

    # generate centroids
    if(isTRUE(verbose)) message('Calculating polygon centroids...')
    xenium_gobject = addSpatialCentroidLocations(xenium_gobject,
                                                 poly_info = c(names(bound_dt_list)))

    # add in feature metadata
    xenium_gobject = addFeatMetadata(gobject = xenium_gobject,
                                     new_metadata = feat_meta,
                                     by_column = TRUE,
                                     column_feat_ID = 'feat_ID')
  }

  if(data_to_use == 'aggregate') {
    # clean up names for aggregate matrices
    names(agg_expr) = gsub(pattern = ' ', replacement = '_' ,names(agg_expr))
    geneExpMat = which(names(agg_expr) == 'Gene_Expression')
    names(agg_expr)[[geneExpMat]] = 'raw'

    # set cell_id as character
    cell_meta = cell_meta[, data.table::setnames(.SD, 'cell_id', 'cell_ID')]
    cell_meta = cell_meta[, cell_ID := as.character(cell_ID)]

    # set up spatial locations
    agg_spatlocs = cell_meta[, .(x_centroid, y_centroid, cell_ID)]

    # set up metadata
    agg_meta = cell_meta[, !c('x_centroid','y_centroid')]

    if(isTRUE(verbose)) message('Building aggregate giotto object...')
    xenium_gobject = createGiottoObject(expression = agg_expr,
                                        spatial_locs = agg_spatlocs,
                                        instructions = instructions,
                                        cores = cores,
                                        verbose = verbose)

    # append aggregate metadata
    xenium_gobject = addCellMetadata(gobject = xenium_gobject,
                                     new_metadata = agg_meta,
                                     by_column = TRUE,
                                     column_cell_ID = 'cell_ID')
    xenium_gobject = addFeatMetadata(gobject = xenium_gobject,
                                     new_metadata = feat_meta,
                                     by_column = TRUE,
                                     column_feat_ID = 'feat_ID')

  }

  return(xenium_gobject)

}









#### logging of giotto functions ####

#' @title Log args used
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



#' @title Update giotto parameters
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



#' @title Giotto object history
#' @name objHistory
#' @description Print and return giotto object history
#' @param object giotto object
#' @export
objHistory = function(object) {
  cat('Steps and parameters used: \n \n')
  print(object@parameters)
  cat('\n\n')
  invisible(x = object@parameters)
}



#### joining giotto object ####

#' @title join_expression_matrices
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

#' @title join_spatlocs
#' @name join_spatlocs
#' @keywords internal
join_spatlocs = function(dt_list) {

  final_list = do.call('rbind', dt_list)
  return(final_list)
}

#' @title join_cell_meta
#' @name join_cell_meta
#' @keywords internal
join_cell_meta = function(dt_list) {

  final_list = do.call('rbind', dt_list)
  return(final_list)

}


#' @title Join giotto objects
#' @name joinGiottoObjects
#' @description Function to join multiple giotto objects together
#' @param gobject_list list of giotto objects
#' @param gobject_names unique giotto names for each giotto object
#' @param join_method method to join giotto objects, see details
#' @param z_vals distance(s) along z-axis if method is z-stack (default is step of 1000)
#' @param x_shift list of values to shift along x-axis if method is shift
#' @param y_shift list of values to shift along y-axis if method is shift
#' @param x_padding padding between datasets/images if method is shift
#' @param y_padding padding between datasets/images if method is shift
#' @param verbose be verbose
#' @return giotto object
#' @details This function joins both the expression and spatial information of
#' multiple giotto objects into a single one. Giotto supports multiple ways of
#' joining spatial information as selected through param \code{join_method}:
#'
#' \itemize{
#'   \item{\strong{\code{"shift"}}} {
#'     (default) Spatial locations of different datasets are shifted
#'     by numeric vectors of values supplied through \code{x_shift}, \code{y_shift},
#'     \code{x_padding}, and \code{y_padding}. This is particularly useful for data
#'     that is provided as tiles or ROIs or when analyzing multiple spatial datasets
#'     together and keeping their spatial data separate.
#'
#'     \strong{If shift values are given then a value is needed for each giotto object
#'     to be joined in \code{gobject_list}. Order matters.}
#'
#'     If a regular step value is desired instead of a specific list of values, use
#'     \code{x_padding} and \code{y_padding}. Both shift and padding values can be used
#'     at the same time.
#'
#'     Leaving \code{x_shift} and \code{y_shift} values as \code{NULL} will have Giotto
#'     estimate an appropriate \code{x_shift} value based on the x dimension of
#'     available image objects. If no image objects are available, a default behavior of
#'     \code{x_padding = 1000} will be applied.
#'   }
#'   \item{\strong{\code{"z_stack"}}} {
#'     Datasets are spatially combined with no change to x and y
#'     spatial locations, but a z value is incorporated for each dataset based on input
#'     supplied through param \code{z_vals}. To specify a z value for each dataset to
#'     join, a numeric vector must be given with a value for each element in \code{gobject_list}.
#'     Order matters.
#'
#'     Alternatively, a single numeric value can be supplied to \code{z_vals} in which
#'     case this input will be treated as a z step value.
#'   }
#'   \item{\strong{\code{"no_change"}}} {
#'     No changes are applied to the spatial locations of the datasets when joining.
#'   }
#' }
#'
#' @concept giotto
#' @export
joinGiottoObjects = function(gobject_list,
                             gobject_names = NULL,
                             join_method = c('shift', 'z_stack', 'no_change'),
                             z_vals = 1000,
                             x_shift = NULL,
                             y_shift = NULL,
                             x_padding = NULL,
                             y_padding = NULL,
                             verbose = TRUE) {

  # define for data.table := and .()
  sdimz = NULL
  cell_ID = NULL
  sdimx = NULL
  sdimy = NULL

  n_gobjects = length(gobject_list)

  ## check general input params
  if(n_gobjects == 0) stop('A list of Giotto objects to be joined must be provided.')
  if(n_gobjects == 1) stop('Only one gobject provided in gobject_list.')
  if(!is.vector(gobject_names) | !is.character(gobject_names)) stop('gobject_names need to be a vector with unique names for the giotto objects')
  if(n_gobjects != length(gobject_names)) stop('each giotto object in the list needs to have a unique (short) name')
  if(is.null(join_method)) message('No join_method given. Defaulting to "shift"')


  ## determine join method
  join_method = match.arg(arg = join_method, choices = c('shift', 'z_stack', 'no_change'))
  if(isTRUE(verbose)) message('Join method:', join_method)



  # **** For shift workflow ****
  if(join_method == 'shift') {
    # Make sure enough x_shift and y_shift values are given to cover all gobjects
    if(!is.null(x_shift)) if(length(x_shift) != n_gobjects) stop('A numeric vector with an x_shift value for each gobject in gobject_list must be given.\n')
    if(!is.null(y_shift)) if(length(y_shift) != n_gobjects) stop('A numeric vector with a y_shift value for each gobject in gobject_list must be given.\n')

    # Set defaults if no shift params are given
    if(is.null(x_shift) & is.null(y_shift) & is.null(x_padding) & is.null(y_padding)) {
      message('No xy shift or specific padding values given. Using defaults')
      x_padding = 1000
    }
    # Assign default padding values if NULL
    if(is.null(x_padding)) x_padding = 0
    if(is.null(y_padding)) y_padding = 0
  }



  # **** For no_change workflow ****
  if(join_method == 'no_change') {
    join_method = 'shift'
    x_shift = rep(0,n_gobjects)
    y_shift = rep(0,n_gobjects)
    x_padding = 0
    y_padding = 0
  }



  # **** For z_stack workflow ****
  if(join_method == 'z_stack') {
    if(!(is.atomic(z_vals) && is.numeric(z_vals))) {
      stop('z_vals requires either a single numeric or an atomic vector of numerics with one value for each z slice (Giotto object). \n')
    }
    if((length(z_vals) != n_gobjects) && (length(z_vals) != 1)) {
      stop('If more than one z_value is given, there must be one for each Giotto object to be joined. \n')
    }

    # expand z_vals if given as a step value
    if(length(z_vals) == 1) {
      if(isTRUE(verbose)) message('Only one value given through z_vals param\n Treating this value as a z step')
      z_vals = ((1:n_gobjects) - 1) * z_vals # Find z vals stepwise
    }
  }



  # keep instructions from first giotto object
  first_instructions = gobject_list[[1]]@instructions

  # keep features from all giotto objects
  vect = vector()
  for(obj_i in 1:n_gobjects) {
    obj_i_features = gobject_list[[obj_i]]@expression_feat
    vect = c(vect, obj_i_features)
  }
  first_features = unique(vect)



  updated_object_list = list()

  ## 0. re-scale spatial locations ##
  ## ----------------------------- ##






  ## 1. update giotto objects ##
  ## ------------------------ ##
  if(verbose == TRUE) message('start updating objects')

  all_feat_ID_list = list()
  all_cell_ID_list = list()
  all_image_list = list()
  all_largeImage_list = list()

  xshift_list = list()
  yshift_list = list()

  all_spatinfo_list = list()


  for(gobj_i in 1:n_gobjects) {

    gobj = gobject_list[[gobj_i]]
    gname = gobject_names[[gobj_i]]


    ## 0. update cell ID and feat ID
    for(spat_unit in names(gobj@cell_ID)) {
      gobj@cell_ID[[spat_unit]] = paste0(gname,'-',gobj@cell_ID[[spat_unit]])
      all_cell_ID_list[[spat_unit]][[gobj_i]] = gobj@cell_ID[[spat_unit]]
    }



    for(feat_type in names(gobj@feat_ID)) {
      all_feat_ID_list[[feat_type]][[gobj_i]] = gobj@feat_ID[[feat_type]]
    }



    ## 1. update expression and all feature IDs
    # provide unique cell ID name
    for(spat_unit in names(gobj@expression)) {

      for(feat_type in names(gobj@expression[[spat_unit]])) {

        for(matr in names(gobj@expression[[spat_unit]][[feat_type]])) {
          colnames(gobj@expression[[spat_unit]][[feat_type]][[matr]]) = gobj@cell_ID[[spat_unit]]
        }

        #all_feat_ID_list[[feat_type]][[gobj_i]] = gobj@feat_ID[[feat_type]]
      }

    }




    ## 2. update images
    # change individual names

    images_found = !is.null(gobj@images)

    if(images_found) {

      names(gobj@images) = paste0(gname, '-', names(gobj@images))
      for(imname in names(gobj@images)) {

        gobj@images[[imname]]@name = paste0(gname,'-', gobj@images[[imname]]@name)


        if(join_method == 'shift') {


          ## shift in x-direction
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
            x_shift_i = x_shift[[gobj_i]]
            add_to_x = x_shift_i + (x_padding * (gobj_i - 1))
          }

          if(verbose) cat('Image: for ',imname, ' add_to_x = ', add_to_x, '\n')

          gobj@images[[imname]]@minmax[c("xmax_sloc", "xmin_sloc")] =  gobj@images[[imname]]@minmax[c("xmax_sloc", "xmin_sloc")] + add_to_x
          xshift_list[[gobj_i]] = add_to_x


          ## shift in y-direction
          if(!is.null(y_shift)) {
            y_shift_i = y_shift[[gobj_i]]
            add_to_y = y_shift_i + (y_padding * (gobj_i - 1))

            if(verbose) cat('Image: for ',imname, ' add_to_y = ', add_to_y, '\n')

            gobj@images[[imname]]@minmax[c("ymax_sloc", "ymin_sloc")] =  gobj@images[[imname]]@minmax[c("ymax_sloc", "ymin_sloc")] + add_to_y
            yshift_list[[gobj_i]] = add_to_y
          }



        }

        all_image_list[[imname]] = gobj@images[[imname]]

      }


    }


    ## 2.2 update largeImages
    # change individual names

    images_found = !is.null(gobj@largeImages)

    if(images_found) {

      names(gobj@largeImages) = paste0(gname, '-', names(gobj@largeImages))
      for(imname in names(gobj@largeImages)) {

        gobj@largeImages[[imname]]@name = paste0(gname,'-', gobj@largeImages[[imname]]@name)


        if(join_method == 'shift') {


          ## shift in x-direction (always happens if not already defined during giottoImage section)
          if(!list_element_exists(xshift_list, gobj_i)) {
            if(is.null(x_shift)) {

              # estimate x_shift step directly from giotto image
              extent = terra::ext(gobj@largeImages[[imname]]@raster_object)

              xmax = extent$xmax[[1]]
              xmin = extent$xmin[[1]]

              add_to_x = ((gobj_i - 1) * (xmax-xmin)) + ((gobj_i - 1) * x_padding)

            } else {
              x_shift_i = x_shift[[gobj_i]]
              add_to_x = x_shift_i + (x_padding * (gobj_i - 1))
            }

            # record xshift (if not already done)
            xshift_list[[gobj_i]] = add_to_x
          }


          if(verbose) cat('largeImage: for ',imname, ' add_to_x = ', add_to_x, '\n')



          gobj@largeImages[[imname]]@raster_object =
            terra::shift(gobj@largeImages[[imname]]@raster_object, dx = xshift_list[[gobj_i]])


          ## shift in y-direction (only happens when y_shift is provided)
          if(!is.null(y_shift)) {
            if(!list_element_exists(yshift_list, gobj_i)) {
              y_shift_i = y_shift[[gobj_i]]
              add_to_y = y_shift_i + (y_padding * (gobj_i - 1))

              yshift_list[[gobj_i]] = add_to_y
            }


            if(verbose) cat('largeImage: for ',imname, ' add_to_y = ', add_to_y, '\n')


            gobj@largeImages[[imname]]@raster_object =
              terra::shift(gobj@largeImages[[imname]]@raster_object, dy = yshift_list[[gobj_i]])

          }



        }

        all_largeImage_list[[imname]] = gobj@largeImages[[imname]]

      }


    }




    ## 3. update spatial location
    # add padding to x-axis
    # update cell ID

    # If no images were present
    if(length(xshift_list) == 0) xshift_list = ((seq_along(gobject_list) - 1) * x_padding)

    for(spat_unit in names(gobj@spatial_locs)) {

      for(locs in names(gobj@spatial_locs[[spat_unit]])) {
        myspatlocs = gobj@spatial_locs[[spat_unit]][[locs]]

        if(join_method == 'z_stack') {
          myspatlocs[, sdimz := z_vals[gobj_i]]
          myspatlocs[, cell_ID := gobj@cell_ID[[spat_unit]] ]
          myspatlocs = myspatlocs[,.(sdimx, sdimy, sdimz, cell_ID)]

        } else if(join_method == 'shift') {


          # shift for x-axis
          if(is.null(x_shift)) {
            add_to_x = xshift_list[[gobj_i]]
          } else {
            x_shift_i = x_shift[[gobj_i]]
            add_to_x = x_shift_i + (x_padding * (gobj_i - 1))
          }

          if(verbose) cat('Spatial locations: for ',locs, ' add_to_x = ', add_to_x, '\n')


          myspatlocs[, sdimx := sdimx + add_to_x]
          myspatlocs[, cell_ID := gobj@cell_ID[[spat_unit]] ]
        }

        # shift for y-axis
        if(!is.null(y_shift)) {
          y_shift_i = y_shift[[gobj_i]]
          add_to_y = y_shift_i + (y_padding * (gobj_i - 1))

          if(verbose) cat('Spatial locations: for ',locs, ' add_to_y = ', add_to_y, '\n')

          myspatlocs[, sdimy := sdimy + add_to_y]
        }


        gobj@spatial_locs[[spat_unit]][[locs]] = myspatlocs
      }
    }





    ## 4. cell metadata
    # rbind metadata
    # create capture area specific names
    for(spat_unit in names(gobj@cell_metadata)) {
      for(feat_type in names(gobj@cell_metadata[[spat_unit]])) {
        gobj@cell_metadata[[spat_unit]][[feat_type]][['cell_ID']] = gobj@cell_ID[[spat_unit]]
        gobj@cell_metadata[[spat_unit]][[feat_type]][['list_ID']] = gname
      }
    }



    ## 5. prepare spatial information
    spatinfo_vector = vector()
    for(spat_info in names(gobj@spatial_info)) {

      # spatVector
      poly_ids = gobj@spatial_info[[spat_info]]@spatVector$poly_ID
      gobj@spatial_info[[spat_info]]@spatVector$poly_ID = paste0(gname,'-',poly_ids)

      # spatVectorCentroids
      if(!is.null(gobj@spatial_info[[spat_info]]@spatVectorCentroids)) {
        poly_ids = gobj@spatial_info[[spat_info]]@spatVectorCentroids$poly_ID
        gobj@spatial_info[[spat_info]]@spatVectorCentroids$poly_ID = paste0(gname,'-',poly_ids)
      }


      # overlaps??
      # TODO


      if(join_method == 'shift') {

        ## for x-axis
        if(is.null(x_shift)) {
          add_to_x = xshift_list[[gobj_i]]
        } else {
          x_shift_i = x_shift[[gobj_i]]
          add_to_x = x_shift_i + (x_padding * (gobj_i - 1))
        }

        ## for y-axis
        if(!is.null(y_shift)) {
          y_shift_i = y_shift[[gobj_i]]
          add_to_y = y_shift_i + (y_padding * (gobj_i - 1))
        } else {
          add_to_y = 0
        }

        if(verbose) cat('Spatial info: for ',spat_info, ' add_to_x = ', add_to_x, '\n')
        if(verbose) cat('Spatial info: for ',spat_info, ' add_to_y = ', add_to_y, '\n')

        gobj@spatial_info[[spat_info]]@spatVector = terra::shift(x = gobj@spatial_info[[spat_info]]@spatVector, dx = add_to_x, dy = add_to_y)
        if(!is.null(gobj@spatial_info[[spat_info]]@spatVectorCentroids)) {
          gobj@spatial_info[[spat_info]]@spatVectorCentroids = terra::shift(x = gobj@spatial_info[[spat_info]]@spatVectorCentroids, dx = add_to_x, dy = add_to_y)
        }


      }

      spatinfo_vector = c(spatinfo_vector, spat_info)
      all_spatinfo_list[[gobj_i]] = spatinfo_vector
    }


    ## 5. prepare spatial information
    for(feat_info in names(gobj@feat_info)) {

      # spatVector
      feat_ids_uniq = gobj@feat_info[[feat_info]]@spatVector$feat_ID_uniq
      gobj@feat_info[[feat_info]]@spatVector$feat_ID_uniq = paste0(gname,'-',feat_ids_uniq)

      # networks??
      # TODO


      if(join_method == 'shift') {

        ## for x-axis
        if(is.null(x_shift)) {
          add_to_x = xshift_list[[gobj_i]]
        } else {
          x_shift_i = x_shift[[gobj_i]]
          add_to_x = x_shift_i + (x_padding * (gobj_i - 1))
        }

        ## for y-axis
        if(!is.null(y_shift)) {
          y_shift_i = y_shift[[gobj_i]]
          add_to_y = y_shift_i + (y_padding * (gobj_i - 1))
        } else {
          add_to_y = 0
        }

        if(verbose) cat('Feature info: for ',feat_info, ' add_to_x = ', add_to_x, '\n')
        if(verbose) cat('Feature info: for ',feat_info, ' add_to_y = ', add_to_y, '\n')

        gobj@feat_info[[feat_info]]@spatVector = terra::shift(x = gobj@feat_info[[feat_info]]@spatVector, dx = add_to_x, dy = add_to_y)
      }




    }


    updated_object_list[[gobj_i]] = gobj

  }

  #return(updated_object_list)




  ## 2. prepare for new giotto object ##
  ## -------------------------------- ##
  comb_gobject = new('giotto',
                     expression = list(),
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

  # comb_gobject = Giotto:::giotto(expression = list(),
  #                                expression_feat = first_features,
  #                                spatial_locs = NULL,
  #                                spatial_info = NULL,
  #                                cell_metadata = NULL,
  #                                feat_metadata = NULL,
  #                                feat_info = NULL,
  #                                cell_ID = NULL,
  #                                feat_ID = NULL,
  #                                spatial_network = NULL,
  #                                spatial_grid = NULL,
  #                                spatial_enrichment = NULL,
  #                                dimension_reduction = NULL,
  #                                nn_network = NULL,
  #                                images = NULL,
  #                                parameters = NULL,
  #                                offset_file = NULL,
  #                                instructions = first_instructions,
  #                                OS_platform = .Platform[['OS.type']],
  #                                join_info = NULL)




  ## 3. merge updated data  ##
  ## ------------------------ ##

  first_obj = updated_object_list[[1]]

  ## cell IDs
  for(spat_unit in names(all_cell_ID_list)) {
    combined_cell_ID = unlist(all_cell_ID_list[[spat_unit]])
    comb_gobject@cell_ID[[spat_unit]] = combined_cell_ID
  }

  ## feat IDs
  for(feat_type in names(all_feat_ID_list)) {
    combined_feat_ID = unique(unlist(all_cell_ID_list[[feat_type]]))
    comb_gobject@feat_ID[[feat_type]] = combined_feat_ID
  }



  ## expression and feat IDs
  ## if no expression matrices are provided, then just combine all feature IDs
  if(verbose == TRUE) message('start expression combination')

  expr_names = names(first_obj@expression)

  if(is.null(expr_names)) {
    ## feat IDS

    for(feat in first_features) {
      combined_feat_ID = unique(unlist(all_feat_ID_list[[feat]]))
      comb_gobject@feat_ID[[feat]] = combined_feat_ID

      # TODO: figure out best way
      comb_gobject@feat_metadata[[spat_unit]][[feat]] = data.table::data.table(feat_ID = combined_feat_ID)

    }



  } else {

    for(spat_unit in names(first_obj@expression)) {

      for(feat_type in names(first_obj@expression[[spat_unit]])) {

        for(mode in names(first_obj@expression[[spat_unit]][[feat_type]])) {

          savelist = list()
          for(gobj_i in 1:length(updated_object_list)) {

            mat = updated_object_list[[gobj_i]]@expression[[spat_unit]][[feat_type]][[mode]]
            savelist[[gobj_i]] = mat
          }

          combmat = join_expression_matrices(matrix_list = savelist)
          comb_gobject@expression[[spat_unit]][[feat_type]][[mode]] = combmat$matrix

          comb_gobject@feat_ID[[feat_type]] = combmat$sort_all_feats
          comb_gobject@feat_metadata[[spat_unit]][[feat_type]] = data.table::data.table(feat_ID = combmat$sort_all_feats)

        }


      }

    }

  }





  ## spatial locations
  if(verbose == TRUE) cat('start spatial location combination \n')


  for(spat_unit in names(first_obj@spatial_locs)) {

    for(name in names(first_obj@spatial_locs[[spat_unit]])) {

      savelist = list()
      for(gobj_i in 1:length(updated_object_list)) {
        spatlocs = updated_object_list[[gobj_i]]@spatial_locs[[spat_unit]][[name]]
        savelist[[gobj_i]] = spatlocs
      }

      combspatlocs = join_spatlocs(dt_list = savelist)
      comb_gobject@spatial_locs[[spat_unit]][[name]] = combspatlocs
    }

  }





  ## cell metadata
  if(verbose == TRUE) cat('start cell metadata combination \n')

  for(spat_unit in names(first_obj@cell_metadata)) {

     for(feat_type in names(first_obj@cell_metadata[[spat_unit]])) {

      savelist = list()
      for(gobj_i in 1:length(updated_object_list)) {
        cellmeta = updated_object_list[[gobj_i]]@cell_metadata[[spat_unit]][[feat_type]]
        savelist[[gobj_i]] = cellmeta
      }

      combcellmeta = join_cell_meta(dt_list = savelist)
      comb_gobject@cell_metadata[[spat_unit]][[feat_type]] = combcellmeta

    }
  }






  ## spatial info
  if(verbose == TRUE) cat('start spatial polygon combination \n')
  available_spat_info = unique(unlist(all_spatinfo_list))

  cat('available_spat_info: \n')
  print(available_spat_info)

  for(spat_info in available_spat_info) {

    savelist_vector = list()
    savelist_centroids = list()
    for(gobj_i in 1:length(updated_object_list)) {

      print(gobj_i)
      print(updated_object_list[[gobj_i]]@spatial_info[[spat_info]])


      spat_information_vector = updated_object_list[[gobj_i]]@spatial_info[[spat_info]]@spatVector
      savelist_vector[[gobj_i]] = spat_information_vector

      spat_information_centroids = updated_object_list[[gobj_i]]@spatial_info[[spat_info]]@spatVectorCentroids
      savelist_centroids[[gobj_i]] = spat_information_centroids
    }

    comb_spatvectors = do.call('rbind', savelist_vector)
    comb_spatcentroids = do.call('rbind', savelist_centroids)

    comb_polygon = create_giotto_polygon_object(name = spat_info,
                                                spatVector = comb_spatvectors,
                                                spatVectorCentroids = comb_spatcentroids,
                                                overlaps = NULL)


    comb_gobject@spatial_info[[spat_info]] = comb_polygon

  }



  ## feature info
  if(verbose == TRUE) cat('start spatial points combination \n')


  for(feat in comb_gobject@expression_feat) {

    savelist_vector = list()

    for(gobj_i in 1:length(updated_object_list)) {

      print(updated_object_list[[gobj_i]]@feat_info[[feat]])
      print(gobj_i)

      if(is.null(updated_object_list[[gobj_i]]@feat_info)) {
        spat_point_vector = NULL
      } else {
        spat_point_vector = updated_object_list[[gobj_i]]@feat_info[[feat]]@spatVector
      }


      savelist_vector[[gobj_i]] = spat_point_vector

      # add network

    }

    comb_spatvectors = do.call('rbind', savelist_vector)

    if(is.null(comb_spatvectors)) {
      comb_points = NULL
    } else {
      comb_points = create_giotto_points_object(feat_type = feat,
                                                spatVector = comb_spatvectors,
                                                networks = NULL)
    }

    comb_gobject@feat_info[[feat]] = comb_points

  }






  ## images
  if(verbose == TRUE) cat('start image\n')

  # keep individual images
  # each individual image has updated x and y locations
  # so all images can be viewed together by plotting them one-by-one
  # but images can also be easily viewed separately by grouping them
  comb_gobject@images = all_image_list
  comb_gobject@largeImages = all_largeImage_list


  ## TODO:
  # update giotto object with join-information
  # - list ID names
  # - xshift values

  # add option to perform yshift

  comb_gobject@join_info = list(list_IDs = gobject_names,
                                join_method = join_method,
                                z_vals = z_vals,
                                x_shift = x_shift,
                                y_shift = y_shift,
                                x_padding = x_padding,
                                y_padding = y_padding)

  return(comb_gobject)


}



