



#### Giotto python path ####

#' @title set_giotto_python_path
#' @name set_giotto_python_path
#' @description sets the python path and/or installs miniconda and the python modules
set_giotto_python_path = function(python_path = NULL,
                                  packages_to_install = c('pandas==1.5.1',
                                                          'networkx==2.8.8',
                                                          'python-igraph==0.10.2',
                                                          'leidenalg==0.9.0',
                                                          'python-louvain==0.16',
                                                          'python.app==1.4',
                                                          'scikit-learn==1.1.3')) {

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

  instructions_list = create_giotto_instructions(
    python_path = python_path,
    show_plot = show_plot,
    return_plot = return_plot,
    save_plot = save_plot,
    save_dir = save_dir,
    plot_format = plot_format,
    dpi = dpi,
    units = units,
    height = height,
    width = width,
    is_docker = is_docker
  )

  return(instructions_list)

}


#' @keywords internal
create_giotto_instructions = function(python_path = NULL,
                                      show_plot = NULL,
                                      return_plot = NULL,
                                      save_plot = NULL,
                                      save_dir = NULL,
                                      plot_format = NULL,
                                      dpi = NULL,
                                      units = NULL,
                                      height = NULL,
                                      width = NULL,
                                      is_docker = NULL) {
  instructions_list = list(python_path = python_path,
                           show_plot = show_plot,
                           return_plot = return_plot,
                           save_plot = save_plot,
                           save_dir = save_dir,
                           plot_format = plot_format,
                           dpi = dpi,
                           units = units,
                           height = height,
                           width = width,
                           is_docker = is_docker)
  class(instructions_list) = c('giottoInstructions', 'list')
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
    stop('\t parameter ', param, ' is not in Giotto instructions \t')
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
#' @description Function to change one or more instructions from giotto object.
#' If more than one item is supplied to \code{params} and \code{new_values}, use
#' a vector of values. Does not call \code{initialize} on the giotto object
#' @param gobject giotto object
#' @param params parameter(s) to change
#' @param new_values new value(s) for parameter(s)
#' @param return_gobject (boolean, default = TRUE) return giotto object
#' @param init_gobject (boolean, default = TRUE) initialize gobject if returning
#' @return giotto object with one or more changed instructions
#' @export
changeGiottoInstructions = function(gobject,
                                    params = NULL,
                                    new_values = NULL,
                                    return_gobject = TRUE,
                                    init_gobject = TRUE) {

  instrs = gobject@instructions

  if(is.null(params) | is.null(new_values)) {
    stop('\t params and new_values can not be NULL \t')
  }

  if(length(params) != length(new_values)) {
    stop('\t length of params need to be the same as new values \t')
  }

  # if(!all(params %in% names(instrs))) {
  #   stop('\t all params need to be part of Giotto instructions \t')
  # }

  ## swap with new values
  instrs[params] = new_values

  ## make sure that classes remain consistent
  new_instrs = lapply(1:length(instrs), function(x) {

    if(names(instrs[x]) %in% c('dpi', 'height', 'width')) {
      instrs[[x]] = as.numeric(instrs[[x]])
    } else if(names(instrs[x]) %in% c('show_plot', 'return_plot', 'save_plot', 'is_docker')) {
      instrs[[x]] = as.logical(instrs[[x]])
    } else if(names(instrs[x]) %in% c('active_spat_unit', 'active_feat_type', 'plot_format', 'units')) {
      instrs[[x]] = as.character(instrs[[x]])
    } else {
      instrs[[x]] = instrs[[x]]
    }

  })

  names(new_instrs) = names(instrs)



  if(isTRUE(return_gobject)) {
    gobject@instructions = new_instrs
    if(isTRUE(init_gobject)) gobject = initialize(gobject)
    return(gobject)
  } else {
    return(new_instrs)
  }

}



#' @title Replace all giotto instructions in giotto object
#' @name replaceGiottoInstructions
#' @description Function to replace all instructions from giotto object. Does
#' not call \code{initialize} on the giotto object
#' @param gobject giotto object
#' @param instructions new instructions (e.g. result from createGiottoInstructions)
#' @param init_gobject (boolean, default = TRUE) initialize gobject when returning
#' @return giotto object with replaces instructions
#' @export
replaceGiottoInstructions = function(gobject,
                                     instructions = NULL,
                                     init_gobject = TRUE) {

  instrs_needed = names(create_giotto_instructions())

  # validate new instructions
  if(!all(instrs_needed %in% names(instructions)) | is.null(instructions)) {
    stop(wrap_txt('You need to provide a named list for all instructions,',
                  'like the outcome of createGiottoInstructions',
                  errWidth = TRUE))
  } else {
    gobject@instructions = instructions
    if(isTRUE(init_gobject)) gobject = initialize(gobject)
    return(gobject)
  }

}






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
                          cores = determine_cores(),
                          transpose = FALSE) {

  # check if path is a character vector and exists
  if(!is.character(path)) stop('path needs to be character vector')
  if(!file.exists(path)) stop('the path: ', path, ' does not exist')

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
                                cores = determine_cores()) {


  if(inherits(inputmatrix, 'character')) {
    inputmatrix = path.expand(inputmatrix)
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

  } else if(inherits(inputmatrix, 'exprObj')) {

    inputmatrix[] = evaluate_expr_matrix(inputmatrix[], sparse = sparse, cores = cores)
    mymatrix = inputmatrix

  }
  else {
    stop(wrap_txt("expression input needs to be a path to matrix-like data or an",
                  "object of class 'Matrix', 'data.table', 'data.frame' or 'matrix'",
                  errWidth = TRUE))
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
  if(is.list(this) && length(this) == 0L) {
    return(0L)
  }
  # Stop if object is not a list AND recurse if it is.
  # Report minimum or maximum depth depending on method
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
#' @param provenance provenance information
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
                                cores = determine_cores(),
                                default_feat_type = NULL,
                                verbose = TRUE,
                                provenance = NULL) {

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
  if(list_depth == 0L) {
    stop('Depth of expression list is 0, no expression information is provided \n')
  }

  # too much information
  if(list_depth > 3L) {
    stop('Depth of expression list is more than 3, only 3 levels are possible:
       0)', ch$s, '.
       1)', ch$s, ch$b, 'spatial unit (e.g. cell)
       2)', ch$s, ch$s, ch$b, 'feature (e.g. RNA)
       3)', ch$s, ch$s, ch$s, ch$b, 'data type (e.g. raw)\n')
  }


  return_list = list()


  # 2. for list with 1 depth
  if(list_depth == 1L) {

    if(isTRUE(verbose)) message('list depth of 1')


    # assign exprObj tagged name info if not named
    expr_list = assign_objnames_2_list(expr_list, force_replace = FALSE)
    expr_list = assign_listnames_2_obj(expr_list)


    for(data in names(expr_list)) {

      # print(expr_list[[data]])
      res_mat = evaluate_expr_matrix(inputmatrix = expr_list[[data]],
                                     sparse = sparse,
                                     cores = cores)
      # add default feat == 'rna'
      # add default region == 'cell'
      if(!inherits(res_mat, 'exprObj')) {
        exprObj = create_expr_obj(name = data,
                                  exprMat = res_mat,
                                  spat_unit = 'cell',
                                  provenance = if(is.null(provenance)) 'cell' else provenance,
                                  feat_type = default_feat_type,
                                  misc = NULL)
      } else {
        exprObj = res_mat
      }

      return_list = append(return_list, exprObj)

    }


  } else if(list_depth == 2L) {

    cat('list depth of 2 \n')
    # add default region == 'cell'

    for(feat in names(expr_list)) {


      # assign exprObj tagged name info if not named
      expr_list = assign_objnames_2_list(expr_list, force_replace = FALSE)
      expr_list = assign_listnames_2_obj(expr_list)


      for(data in names(expr_list[[feat]])) {

        res_mat = evaluate_expr_matrix(inputmatrix = expr_list[[feat]][[data]],
                                       sparse = sparse,
                                       cores = cores)
        # add default region == 'cell'
        if(!inherits(res_mat, 'exprObj')) {
          exprObj = create_expr_obj(name = data,
                                    exprMat = res_mat,
                                    spat_unit = 'cell',
                                    provenance = if(is.null(provenance)) 'cell' else provenance,
                                    feat_type = feat,
                                    misc = NULL)
        } else {
          exprObj = res_mat
          featType(exprObj) = feat
        }


        return_list = append(return_list, exprObj)

      }
    }

  } else if(list_depth == 3L) {

    cat('list depth of 3 \n')

    for(region in names(expr_list)) {
      for(feat in names(expr_list[[region]])) {


        # assign exprObj tagged name info if not named
        expr_list = assign_objnames_2_list(expr_list, force_replace = FALSE)
        expr_list = assign_listnames_2_obj(expr_list)


        for(data in names(expr_list[[region]][[feat]])) {

          res_mat = evaluate_expr_matrix(inputmatrix = expr_list[[region]][[feat]][[data]],
                                         sparse = sparse,
                                         cores = cores)
          # add default region == 'cell'
          if(!inherits(res_mat, 'exprObj')) {
            exprObj = create_expr_obj(name = data,
                                      exprMat = res_mat,
                                      spat_unit = region,
                                      provenance = if(is.null(provenance)) region else provenance,
                                      feat_type = feat,
                                      misc = NULL)
          } else {
            exprObj = res_mat
            featType(exprObj) = feat
            spatUnit(exprObj) = region
          }


          return_list = append(return_list, exprObj)

        }
      }
    }



  } else {
    stop('unexpected list_depth error')
  }

  return(return_list)

}



#' @title Initialize cell and feature IDs
#' @name init_cell_and_feat_IDs
#' @description sets cell and feature IDs based on provided expression data.
#' Enforces that across a single spatial unit, all expression matrices MUST
#' have the same set of cell_IDs
#' @keywords internal
#' @noRd
init_cell_and_feat_IDs = function(gobject) {

  spat_unit = feat_type = name = NULL

  # wipe values
  slot(gobject, 'cell_ID') = NULL
  slot(gobject, 'feat_ID') = NULL

  # find available expr and info
  avail_expr = list_expression(gobject)
  avail_si = list_spatial_info(gobject)
  avail_fi = list_feature_info(gobject)

  used_spat_units = unique(c(avail_expr$spat_unit, avail_si$spat_info))
  used_feat_types = unique(c(avail_expr$feat_type, avail_fi$feat_info))

  # 1. set cell_ID for each region
  # each regions can have multiple features, but the cell_IDs (spatial units) should be the same
  # Select spatial unit to initialize then pass to set_cell_id
  # set_cell_id decides which data to initialize from
  for(spatial_unit in used_spat_units) {
    gobject = set_cell_id(gobject = gobject,
                          spat_unit = spatial_unit,
                          cell_IDs = 'initialize',
                          verbose = FALSE,
                          set_defaults = TRUE)
  }

  # 2. ensure cell_ID and colnames for each matrix are the same
  # for(expr_i in seq(avail_expr[, .N])) {
  #   spatial_unit = avail_expr[expr_i, spat_unit]
  #   feature_type = avail_expr[expr_i, feat_type]
  #   data = avail_expr[expr_i, name]
  #
  #   colnames_matrix = colnames(get_expression_values(gobject,
  #                                                    spat_unit = spatial_unit,
  #                                                    feat_type = feature_type,
  #                                                    values = data,
  #                                                    output = 'matrix'))
  #   gobj_cell_ID = get_cell_id(gobject,
  #                              spat_unit = spatial_unit)
  #   if(!identical(colnames_matrix, gobj_cell_ID)) {
  #     stop('Colnames are not the same for feat: ', feature_type,', spatial unit: ', spatial_unit ,', and data: ', data)
  #   }
  # }


  # 2. set feat_ID for each feature
  for(feature_type in used_feat_types) {
    gobject = set_feat_id(gobject = gobject,
                          feat_type = feature_type,
                          feat_IDs = 'initialize',
                          verbose = FALSE,
                          set_defaults = TRUE)
  }


  return(gobject)

}


#### Giotto metadata ####

#' @title Read cell metadata
#' @name read_cell_metadata
#' @description read cell metadata from list
#' @param gobject giotto object
#' @param metadata nested list of cell metadata information
#' @param provenance provenance information (optional)
#' @param verbose be verbose
#' @keywords internal
read_cell_metadata = function(gobject,
                              metadata,
                              provenance = NULL,
                              verbose = TRUE) {

  # data.table vars
  cell_ID = NULL

  if(is.null(metadata)) {
    return(NULL)
  }

  if(inherits(metadata, 'cellMetaObj')) {
    return(list(metadata))
  }


  # list reading



  cellMetaObj_list = list()

  # extract all metadata information
  # need to be nested list (feature type and spatial unit)
  for(spat_unit in names(metadata)) {

    if(!spat_unit %in% list_cell_id_names(gobject)) {
      if(isTRUE(verbose)) warning('spat_unit "', spat_unit, '" not found in gobject cell_ID slot.\nPlease check metadata list nesting.\n')
    }

    for(feat_type in names(metadata[[spat_unit]])) {

      if(!feat_type %in% list_feat_id_names(gobject)) {
        if(isTRUE(verbose)) warning('feat_type "', feat_type, '" not found in gobject feat_ID slot.\nPlease check metadata list nesting.\n')
      }

      # TODO load with fread if path given as character

      metaDT = data.table::as.data.table(metadata[[spat_unit]][[feat_type]])

      # if cell ID col is missing, try to automatically set
      if(is.null(metaDT[['cell_ID']])) {
        id_error = try(metaDT[, cell_ID := get_cell_id(gobject, spat_unit = spat_unit)], silent = TRUE)
        if(inherits(id_error, 'try-error')) stop('cannot automatically set metadata cell_ID based on gobject cell_ID slot.')
      } else if(spat_unit %in% list_cell_id_names(gobject)) {

        # if cell ID col is present in both, try to match
        if(!identical(metaDT[, cell_ID], get_cell_id(gobject, spat_unit = spat_unit))) {
          stop('metadata cell_ID does not match that in gobject cell_ID slot for spat_unit "', spat_unit, '".\n')
        }

      }

      # put cell_ID first
      setcolorder(metaDT, 'cell_ID')

      metaObj = new('cellMetaObj',
                    metaDT = metaDT,
                    col_desc = NA_character_, # unknown
                    spat_unit = spat_unit,
                    provenance = if(is.null(provenance)) spat_unit else provenance,
                    feat_type = feat_type)

      cellMetaObj_list = append(cellMetaObj_list, metaObj)

    }
  }

  return(cellMetaObj_list)
}



#' @title Read feature metadata
#' @name read_feature_metadata
#' @description read feature metadata from list
#' @param gobject giotto object
#' @param metadata nested list of feature metadata information
#' @param provenance provenance information (optional)
#' @param verbose be verbose
#' @keywords internal
read_feature_metadata = function(gobject,
                                 metadata,
                                 provenance = NULL,
                                 verbose = TRUE) {

  # data.table vars
  feat_ID = NULL

  featMetaObj_list = list()

  # extract all metadata information
  # need to be nested list (feature type and spatial unit)
  for(spat_unit in names(metadata)) {

    if(!spat_unit %in% list_cell_id_names(gobject)) {
      if(isTRUE(verbose)) warning('spat_unit "', spat_unit, '" not found in gobject cell_ID slot.\nPlease check metadata list nesting.\n')
    }

    for(feat_type in names(metadata[[spat_unit]])) {

      if(!feat_type %in% list_feat_id_names(gobject)) {
        if(isTRUE(verbose)) warning('feat_type "', feat_type, '" not found in gobject feat_ID slot.\nPlease check metadata list nesting.\n')
      }

      # TODO load with fread if path given as character

      metaDT = data.table::as.data.table(metadata[[spat_unit]][[feat_type]])

      # if feat ID col is missing, try to automatically set
      if(is.null(metaDT[['feat_ID']])) {
        id_error = try(metaDT[, feat_ID := get_feat_id(gobject, feat_type = feat_type)], silent = TRUE)
        if(inherits(id_error, 'try-error')) stop('cannot automatically set metadata feat_ID based on gobject feat_ID slot.')
      } else if(spat_unit %in% list_feat_id_names(gobject)) {

        # if feat ID col is present in both, try to match
        if(!identical(metaDT[, feat_ID], get_feat_id(gobject, feat_type = feat_type))) {
          stop('metadata feat_ID does not match that in gobject feat_ID slot for feat_type "', feat_type, '".\n')
        }

      }

      # put feat_ID first
      all_colnames = colnames(metaDT)
      other_colnames = grep('feat_ID', all_colnames, invert = TRUE, value = TRUE)
      metaDT = metaDT[, c('feat_ID', other_colnames), with = FALSE]

      metaObj = new('featMetaObj',
                    metaDT = metaDT,
                    col_desc = NA_character_, # unknown
                    spat_unit = spat_unit,
                    provenance = if(is.null(provenance)) spat_unit else provenance,
                    feat_type = feat_type)

      featMetaObj_list = append(featMetaObj_list, metaObj)

    }
  }

  return(featMetaObj_list)
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
                                      cores = determine_cores()) {

  # data.table variables
  cell_ID = NULL

  if(!any(class(spatial_locs) %in% c('data.table', 'data.frame', 'matrix', 'character'))) {
    stop('spatial_locs needs to be a data.table or data.frame-like object or a path to one of these')
  }
  if(inherits(spatial_locs, 'character')) {
    if(!file.exists(spatial_locs)) stop('path to spatial locations does not exist')
    spatial_locs = data.table::fread(input = spatial_locs, nThread = cores)
  } else {
    spatial_locs = tryCatch(
      data.table::setDT(spatial_locs),
      error = function(e) data.table::as.data.table(spatial_locs)
    )
  }

  # check if all columns are numeric
  column_classes = lapply(spatial_locs, FUN = class)
  non_numeric_classes = column_classes[!column_classes %in% c('numeric','integer')]


  potential_cell_IDs = NULL

  # find non-numeric cols (possible cell_ID col)
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

  # Assign first non-numeric as cell_ID
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
#' @param provenance provenance information (optional)
#' @param verbose be verbose
#' @return updated giotto object
#' @keywords internal
#' @noRd
read_spatial_location_data = function(gobject,
                                      spat_loc_list,
                                      cores = determine_cores(),
                                      provenance = NULL,
                                      verbose = TRUE) {

  # data.table vars
  cell_ID = NULL

  if(is.null(spat_loc_list)) return(NULL)

  try_val = try(as.list(spat_loc_list), silent = TRUE)
  if(inherits(try_val, 'try-error')) {
    return(list(spat_loc_list))
  } else {
    spat_loc_list = try_val
  }

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
    stop(wrap_txt('Depth of spatial location list is more than 2, only 2 levels are possible:
                  1) spatial unit (e.g. cell) --> 2) coordinate (e.g. raw) \n',
                  errWidth = TRUE))
  }

  # 2. Based on depth of nesting expect related info then eval, check, and assemble return list
  ### 2.1 evaluate spatlocs - (read) and find col classes and accordingly assign DT and colnames
  ### 2.2 check spatlocs - compare guessed cell_ID col vs gobject cell_ID slot (from expr)
  ### 2.3 create spatloc objects
  return_list = list()

  obj_list = list()
  spat_unit_list = c()
  name_list = c()

  # for list with 1 depth, expect name info
  if(list_depth == 1L) {

    if(isTRUE(verbose)) (message('list depth of 1'))

    obj_names = names(spat_loc_list)
    for(obj_i in seq_along(spat_loc_list)) {
      spat_unit_list = rep('cell', length(spat_loc_list)) # add default region = 'cell'

      spatlocs = spat_loc_list[[obj_i]]
      name = if(is.null(obj_names[[obj_i]])) paste0('coord_', obj_i) else obj_names[[obj_i]]

      obj_list = append(obj_list, spatlocs)
      name_list = c(name_list, name)

      # # add default region == 'cell'
      # return_list = append(return_list,
      #                      create_spat_locs_obj(
      #                        name = coord,
      #                        coordinates = res_spatlocs,
      #                        spat_unit = 'cell',
      #                        provenance = if(is.null(provenance)) 'cell' else provenance
      #                      ))
    }

    # for list with 2 depth, expect name info and spat_unit info
  } else if(list_depth == 2L) {

    if(isTRUE(verbose)) message('list depth of 2')

    spat_unit_names = names(spat_loc_list)
    for(unit_i in seq_along(spat_loc_list)) {

      obj_names = names(spat_loc_list[[unit_i]])
      for(obj_i in seq_along(spat_loc_list[[unit_i]])) {

        spatlocs = spat_loc_list[[obj_i]]
        spat_unit = if(is.null(spat_unit_names[[unit_i]])) paste0('unit_', unit_i) else spat_unit_names[[unit_i]]
        name = if(is.null(obj_names[[obj_i]])) paste0('coord_', obj_i) else obj_names[[obj_i]]

        obj_list = append(obj_list, spatlocs)
        spat_unit_list = c(spat_unit_list, spat_unit)
        name_list = c(name_list, name)

        # return_list = append(return_list,
        #                      create_spat_locs_obj(
        #                        name = coord,
        #                        coordinates = res_spatlocs,
        #                        spat_unit = spat_unit,
        #                        provenance = if(is.null(provenance)) spat_unit else provenance
        #                      ))
      }
    }

  } else {
    stop('unexpected list_depth error')
  }


  # create spatLocObj return lst
  if(length(obj_list > 0)) {

    return_list = lapply(seq_along(obj_list), function(obj_i) {

      if(inherits(obj_list[[obj_i]], 'spatLocsObj')) {
        return(obj_list[[obj_i]])
      } else {

        return(
          create_spat_locs_obj(
            name = name_list[[obj_i]],
            coordinates = obj_list[[obj_i]],
            spat_unit = spat_unit_list[[obj_i]],
            provenance = if(is.null(provenance)) spat_unit_list[[obj_i]] else provenance, # assumed
            misc = NULL
          )
        )
      }
    })
    return(return_list)
  } else {
    warning('No objects found in spatial locations input list')
  }

}







#' @title Check spatial location data
#' @name check_spatial_location_data
#' @description check cell ID (spatial unit) names between spatial location and expression data.
#' It will look for identical IDs after sorting.
#' @keywords internal
check_spatial_location_data = function(gobject) {

  # define for data.table
  cell_ID = spat_unit = name = NULL

  # find available spatial locations
  avail_sl = list_spatial_locations(gobject)
  avail_ex = list_expression(gobject)
  avail_si = list_spatial_info(gobject)

  missing_unit = !(avail_sl$spat_unit) %in% c(avail_ex$spat_unit, avail_si$spat_info)
  if(any(missing_unit)) {
    stop(wrap_txt('No expression or polygon information discovered for spat_unit:',
                  avail_sl$spat_unit[missing_unit],
                  'Please add expression or polygon information for this spatial',
                  'unit first'))
  }

  for(spat_unit_i in avail_sl[['spat_unit']]) {

    expected_cell_ID_names = get_cell_id(gobject = gobject,
                                         spat_unit = spat_unit_i)

    for(coord_i in avail_sl[spat_unit == spat_unit_i, name]) {

      # 1. get colnames
      spatlocsDT = get_spatial_locations(gobject,
                                         spat_unit = spat_unit_i,
                                         spat_loc_name = coord_i,
                                         output = 'data.table',
                                         copy_obj = FALSE)
      missing_cell_IDs = spatlocsDT[, all(is.na(cell_ID))]

      # if cell_ID column is provided then compare with expected cell_IDs
      if(!isTRUE(missing_cell_IDs)) {

        spatial_cell_id_names = spatlocsDT[['cell_ID']]

        spatial_cell_id_names = sort(spatial_cell_id_names)
        expected_cell_ID_names = sort(expected_cell_ID_names)

        if(!identical(spatial_cell_id_names, expected_cell_ID_names)) {
          message('spatloc cell_IDs: ')
          cat('  ', head(spatial_cell_id_names,3), '...', tail(spatial_cell_id_names,3), '\n')
          message('expression cell_IDs: ')
          cat('  ', head(expected_cell_ID_names,3), '...', tail(expected_cell_ID_names,3), '\n')

          stop('cell_IDs between spatial and expression information are not the same for: \n
                 spatial unit: ', spat_unit_i, ' and coordinates: ', coord_i, ' \n')
        }

      } else {

        # if cell_ID column is not provided then add expected cell_IDs

        ## error if spatlocs and cell_ID do not match in length
        if(spatlocsDT[,.N] != length(expected_cell_ID_names)) {
          stop('Number of rows of spatial locations do not match with cell IDs for: \n
                 spatial unit: ', spat_unit_i, ' and coordinates: ', coord_i, ' \n')
        }

        ## ! modify coords within gobject by reference
        spatlocsDT = spatlocsDT[, cell_ID := expected_cell_ID_names]

      }

    }
  }

  return(invisible())

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


# internal for deprecation
xy_translate_spatial_locations = function(...) {
  .Deprecated(new = 'shift_spatial_locations')

  shift_spatial_locations(...)
}


#' @title Shift spatial locations
#' @name shift_spatial_locations
#' @description Shift given coordinates by given translation values
#' @param spatlocs spatial locations to use
#' @param dx value to shift coordinates in the positive x direction
#' @param dy value to shift coordinates in the positive y direction
#' @param dz value to shift coordinates in the positive z direction
#' @param xtranslate deprecated. use dx
#' @param ytranslate deprecated. use dy
#' @param ztranslate deprecated. use dz
#' @param copy_obj copy/duplicate object (default = TRUE)
#' @keywords internal
shift_spatial_locations = function(spatlocs,
                                   dx = 0,
                                   dy = 0,
                                   dz = 0,
                                   xtranslate = NULL,
                                   ytranslate = NULL,
                                   ztranslate = NULL,
                                   copy_obj = TRUE) {
  sdimx = sdimy = sdimz = NULL

  if(!is.null(xtranslate)) {
    warning(wrap_txt('xtranslate is deprecated. use dx'))
    dx = xtranslate
  }
  if(!is.null(ytranslate)) {
    warning(wrap_txt('ytranslate is deprecated. use dy'))
    dy = ytranslate
  }
  if(!is.null(ztranslate)) {
    warning(wrap_txt('ztranslate is deprecated. use dz'))
    dz = ztranslate
  }

  spatlocs[, sdimx := sdimx + dx]
  spatlocs[, sdimy := sdimy + dy]
  if('sdimz' %in% names(spatlocs)) spatlocs[, sdimz := sdimz + dz]

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





#' @title Evaluate spatial network
#' @name evaluate_spatial_network
#' @description function to evaluate a spatial network
#' @keywords internal
evaluate_spatial_network = function(spatial_network) {

  if(inherits(spatial_network, 'spatialNetworkObj')) {
    spatial_network[] = evaluate_spatial_network(spatial_network[])
    return(spatial_network)
  }

  if(!inherits(spatial_network, 'data.frame')) {
    stop(wrap_txt('The spatial network must be a data.frame(-like) object',
                  errWidth = TRUE))
  }
  if(!inherits(spatial_network, 'data.table')) {
    spatial_network = data.table::setDT(spatial_network)
  }

  netw_names = colnames(spatial_network)
  required_cols = c('from', 'to',
                    'sdimx_begin', 'sdimy_begin',
                    'sdimx_end', 'sdimy_end',
                    'distance', 'weight')

  missing_cols = required_cols[!required_cols %in% netw_names]

  if(length(missing_cols) > 0) {
    stop('missing columns: ', list(missing_cols))
  }

  return(spatial_network)
}







#' @title Read spatial networks
#' @name read_spatial_networks
#' @description read spatial networks from list
#' @keywords internal
read_spatial_networks = function(spatial_network,
                                 provenance = NULL) {

  if(is.null(spatial_network)) {
    wrap_msg('No spatial networks are provided')
    return(NULL)
  }

  # if not a list, return directly
  if(!inherits(spatial_network, 'list')) {
    try_val = try(as.list(spatial_network), silent = TRUE)
    if(inherits(try_val, 'try-error')) {
      return(list(spatial_network))
    } else {
      spatial_network = try_val
    }
  }


  # list reading
  obj_list = list()
  spat_unit_list = c()
  name_list = c()


  # read nesting
  if(depth(spatial_network == 1L)) {
    obj_names = names(spatial_network)
    spat_unit_list = rep('cell', length(spatial_network)) # assume

    for(obj_i in seq_along(spatial_network)) {

      network = spatial_network[[obj_i]]
      name = if(is.null(obj_names[[obj_i]])) paste0('sn_', obj_i) else obj_names[[obj_i]]

      obj_list = append(obj_list, network)
      name_list = c(name_list, name)

    }

  } else if(depth(spatial_network == 2L)) {

    spat_unit_names = names(spatial_network)
    for(unit_i in seq_along(spatial_network)) {

      obj_names = names(spatial_network[[unit_i]])
      for(obj_i in seq_along(spatial_network[[unit_i]])) {

        network = spatial_network[[unit_i]][[obj_i]]
        spat_unit = if(is.null(spat_unit_names[[unit_i]])) paste0('unit_', unit_i) else spat_unit_names[[unit_i]]
        name = if(is.null(obj_names[[obj_i]])) paste0('sn_', obj_i) else obj_names[[obj_i]]

        obj_list = append(obj_list, network)
        spat_unit_list = c(spat_unit_list, spat_unit)
        name_list = c(name_list, name)

      }
    }

  } else {
    stop(wrap_txt('Invalid nesting depth'))
  }


  # create spatialNetworkObj return list
  if(length(obj_list > 0)) {

    return_list = lapply(seq_along(obj_list), function(obj_i) {

      if(inherits(obj_list[[obj_i]], 'spatialNetworkObj')) {
        return(obj_list[[obj_i]])
      } else {

        return(
          create_spat_net_obj(
            name = name_list[[obj_i]],
            method = name_list[[obj_i]], # assumed
            spat_unit = spat_unit_list[[obj_i]],
            provenance = if(is.null(provenance)) spat_unit_list[[obj_i]] else provenance, # assumed
            networkDT = obj_list[[obj_i]],
            networkDT_before_filter = NULL,
            cellShapeObj = NULL,
            crossSectionObjects = NULL,
            parameters = NULL,
            outputObj = NULL,
            misc = NULL
          )
        )

      }
    })
    return(return_list)

  } else {
    warning('No objects found in spatial networks input list')
  }



    # for(spat_unit in names(spatial_network)) {
    #   for(name in names(spatial_network[[spat_unit]])) {
    #
    #     # first check if corresponding expression matrix exists
    #     if(!is.null(slot(gobject, 'expression')[[spat_unit]])) {
    #
    #
    #
    #       network = spatial_network[[spat_unit]][[name]]
    #
    #       if(any(c('data.frame', 'data.table') %in% class(network))) {
    #         if(all(c('to', 'from', 'weight', 'sdimx_begin', 'sdimy_begin', 'sdimx_end', 'sdimy_end') %in% colnames(network))) {
    #           spatial_network_Obj = new('spatialNetworkObj',
    #                                     name = name,
    #                                     networkDT = network)
    #           gobject = set_spatialNetwork(gobject = gobject,
    #                                        spat_unit = spat_unit,
    #                                        name = name,
    #                                        spatial_network = spatial_network_Obj)
    #         } else {
    #           warning('\n spatial unit: ', spat_unit, ' with network name: ', name, ' does not have all necessary column names, see details
    #                   and will not be added to the Giotto object \n')
    #         }
    #       } else {
    #         warning('\n spatial unit: ', spat_unit, ' with network name: ', name, ' is not a data.frame or data.table
    #                 and will not be added to the Giotto object \n')
    #       }
    #
    #     }
    #   }
    # }


}


# See function spatShift in generics.R
#' @name shift_spatial_network
#' @title Shift spatial network
#' @description Shift spatial network coordinates
#' @param spatnet spatial network data.table
#' @param dx distance to shift on x axis
#' @param dy distance to shift on y axis
#' @param dz distance to shift on z axis
#' @param copy_obj copy/duplicate object (default = TRUE)
#' @keywords internal
shift_spatial_network = function(spatnet, dx = 0, dy = 0, dz = 0, copy_obj = TRUE) {
  sdimx_begin = sdimx_end = sdimy_begin = sdimy_end = sdimz_begin = sdimz_end = NULL

  # if 3D info present
  is3D = FALSE
  if(all(c('sdimz_begin', 'sdimz_end') %in% colnames(spatnet))) is3D = TRUE

  if(copy_obj) spatnet = data.table::copy(spatnet)

  spatnet[, `:=`(sdimx_begin = sdimx_begin + dx,
                 sdimx_end = sdimx_end + dx,
                 sdimy_begin = sdimy_begin + dy,
                 sdimy_end = sdimy_end + dy)]
  if(is3D) {
    spatnet[, `:=`(sdimz_begin = sdimz_begin + dz,
                   sdimz_end = sdimz_end + dz)]
  }
  return(spatnet)
}


#### Giotto spatial enrichment ####


#' @title Read spatial enrichment
#' @name read_spatial_enrichment
#' @description read spatial enrichment results from list
#' @keywords internal
read_spatial_enrichment = function(spatial_enrichment,
                                   provenance = NULL) {

  if(is.null(spatial_enrichment)) {
    message('No spatial enrichment results are provided')
    return(NULL)
  }

  # return directly if not list
  if(!inherit(spatial_enrichment, 'list')) {
    try_val = try(as.list(spatial_enrichment), silent = TRUE)
    if(inherits(try_val, 'try-error')) {
      return(list(spatial_enrichment))
    } else {
      spatial_enrichment = try_val
    }
  }

  # list reading
  obj_list = list()
  spat_unit_list = c()
  feat_type_list = c()
  name_list = c()



  if(depth(spatial_enrichment) == 1L) {        # ------------------------ 1 #
    obj_names = names(spatial_enrichment)

    for(obj_i in seq_along(spatial_enrichment)) {

      enr = spatial_enrichment[[obj_i]]
      name = if(is.null(obj_names[[obj_i]])) paste0('enr_', obj_i) else obj_names[[obj_i]]

      obj_list = append(obj_list, enr)
      name_list = c(name_list, name)

    }
    feat_type_list = rep('rna', length(obj_list)) # assume
    spat_unit_list = rep('cell', length(obj_list)) # assume

  } else if(depth(spatial_enrichment) == 2L) { # ------------------------ 2 #

    feat_type_names = names(spatial_enrichment)
    for(feat_i in seq_along(spatial_enrichment)) {

      obj_names = names(spatial_enrichment[[feat_i]])
      for(obj_i in seq_along(spatial_enrichment[[feat_i]])) {

        enr = spatial_enrichment[[feat_i]][[obj_i]]
        feat_type = if(is.null(feat_type_names[[feat_i]])) paste0('feat_', feat_i) else feat_type_names[[feat_i]]
        name = if(is.null(obj_names[[obj_i]])) paste0('enr_', obj_i) else obj_names[[obj_i]]

        obj_list = append(obj_list, enr)
        feat_type_list = c(feat_type_list, feat_type)
        name_list = c(name_list, name)

      }
    }
    spat_unit_list = rep('cell', length(obj_list)) # assume

  } else if(depth(spatial_enrichment) == 3L) { # ------------------------ 3 #

    spat_unit_names = names(spatial_enrichment)
    for(unit_i in seq_along(spatial_enrichment)) {

      feat_type_names = names(spatial_enrichment[[unit_i]])
      for(feat_i in seq_along(spatial_enrichment[[unit_i]])) {

        obj_names = names(spatial_enrichment[[unit_i]][[feat_i]])
        for(obj_i in seq_along(spatial_enrichment[[unit_i]][[feat_i]])) {

          enr = spatial_enrichment[[unit_i]][[feat_i]][[obj_i]]
          spat_unit = if(is.null(spat_unit_names[[unit_i]])) paste0('unit_', unit_i) else spat_unit_names[[unit_i]]
          feat_type = if(is.null(feat_type_names[[feat_i]])) paste0('feat_', feat_i) else feat_type_names[[feat_i]]
          name = if(is.null(obj_names[[obj_i]])) paste0('enr_', obj_i) else obj_names[[obj_i]]

          obj_list = append(obj_list, enr)
          spat_unit_list = c(spat_unit_list, spat_unit)
          feat_type_list = c(feat_type_list, feat_type)
          name_list = c(name_list, name)

        }
      }
    }

  } else {
    stop(wrap_txt('Unexpected list depth', errWidth = TRUE))
  }



  # create spatEnrObj return list
  if(length(obj_list > 0)) {

    return_list = lapply(seq_along(obj_list), function(obj_i) {

      if(inherits(obj_list[[obj_i]], 'spatEnrObj')) {
        return(obj_list[[obj_i]])
      } else {

        return(
          create_spat_enr_obj(
            name = name_list[[obj_i]],
            method = name_list[[obj_i]], # assumed
            enrichDT = obj_list[[obj_i]],
            spat_unit = spat_unit_list[[obj_i]],
            feat_type = feat_type_list[[obj_i]],
            provenance = if(is.null(provenance)) spat_unit_list[[obj_i]] else provenance, # assumed
            misc = NULL
          )
        )
      }
    })
    return(return_list)

  } else {
    warning('No objects found in spatial enrichment input list')
  }

  # for(spat_unit in names(spatial_enrichment)) {
  #   for(name in names(spatial_enrichment[[spat_unit]])) {
  #
  #     # first check if corresponding expression matrix exists
  #     if(!is.null(gobject@expression[[spat_unit]])) {
  #
  #
  #       # TODO: use fread if it's an existing path
  #
  #       spat_enrich = spatial_enrichment[[spat_unit]][[name]]
  #
  #       if(nrow(spat_enrich) != ncol(gobject@expression[[spat_unit]][[1]][[1]])) {
  #         stop('\n spatial enrichment for spatial unit: ', spat_unit, ' and name: ', name, ' does not match the corresponding expression data \n')
  #       } else {
  #         gobject@spatial_enrichment[[spat_unit]][[name]] = spat_enrich
  #       }
  #
  #     }
  #
  #   }
  #
  # }

}




#' @name check_spatial_enrichment
#' @description check the spatial enrichment information within the gobject
#' @keywords internal
#' @noRd
check_spatial_enrichment = function(gobject) {

  avail_se = list_spatial_enrichments(gobject = gobject)
  used_su = unique(avail_se$spat_unit)

  if(!is.null(used_su)) {

    for(su_i in used_su) {
      IDs = spatIDs(gobject, spat_unit = su_i)

      su_se = avail_se[spat_unit == su_i,]
      lapply(seq(nrow(avail_se)), function(obj_i) {
        se_obj = get_spatial_enrichment(gobject = gobject,
                                        spat_unit = su_i,
                                        feat_type = su_se$feat_type[[obj_i]],
                                        enrichm_name = su_se$name[[obj_i]],
                                        output = 'spatEnrObj',
                                        copy_obj = FALSE,
                                        set_defaults = FALSE)
        if(!setequal(spatIDs(se_obj), IDs)) {
          warning(wrap_txt('spat_unit:', su_i,
                           'feat_type:', su_se$feat_type[[obj_i]],
                           'name:', su_se$name[[obj_i]], '\n',
                           'Spatial enrichment IDs are not all found in gobject IDs'))
        }
        # print(paste(su_i, su_se_name[[obj_i]])) # debug
      })
    }
  }

}






#### Giotto dimension reduction ####


#' @name evaluate_dimension_reduction
#' @description evaluate dimension reduction input into dimObj matrix
#' @keywords internal
#' @noRd
evaluate_dimension_reduction = function(dimension_reduction) {

  # object level
  if(inherits(dimension_reduction, 'dimObj')) {
    dimension_reduction[] = evaluate_dimension_reduction(dimension_reduction[])
    return(dimension_reduction)
  }

  # coordinates slot matrix
  dimension_reduction = try(as.matrix(dimension_reduction), silent = TRUE)
  if(inherits(dimension_reduction, 'try-error')) {
    stop(wrap_txt('Dimension reduction coordinate input must be coercible to matrix',
                  errWidth = TRUE))
  }
  return(dimension_reduction)

}





#' @title Read dimensional reduction
#' @name read_dimension_reduction
#' @description read dimension reduction results from list
#' @keywords internal
#' @noRd
read_dimension_reduction = function(dimension_reduction,
                                    reduction = c('cells', 'feats'),
                                    provenance = NULL) {

  reduction = match.arg(reduction, choices = c('cells', 'feats'))

  if(is.null(dimension_reduction)) {
    message('No dimension reduction results are provided')
    return(NULL)
  }

  # if not list, return directly
  if(!inherits(dimension_reduction, 'list')) {
    try_val = try(as.list(dimension_reduction), silent = TRUE)
    if(inherits(try_val, 'try-error')) {
      return(list(dimension_reduction))
    } else {
      dimension_reduction = try_val
    }
  }


  # list reading
  obj_list = list()
  spat_unit_list = c()
  feat_type_list = c()
  method_list = c()
  name_list = c()


  # read nesting
  if(depth(dimension_reduction) == 1L) {        # ------------------------ 1 #
    obj_names = names(dimension_reduction)


    for(obj_i in seq_along(dimension_reduction)) {

      dr = dimension_reduction[[obj_i]]
      name = if(is.null(obj_names[[obj_i]])) paste0('dimRed_', obj_i) else obj_names[[obj_i]]
      method = name # assume

      obj_list = append(obj_list, dr)
      name_list = c(name_list, name)
      method_list = c(method_list, method)

    }
    feat_type_list = rep('rna', length(obj_list)) # assume
    spat_unit_list = rep('cell', length(obj_list)) # assume

  } else if(depth(dimension_reduction) == 2L) { # ------------------------ 2 #

    feat_type_names = names(dimension_reduction)
    for(feat_i in seq_along(dimension_reduction)) {

      obj_names = names(dimension_reduction[[feat_i]])
      for(obj_i in seq_along(dimension_reduction[[feat_i]])) {

        dr = dimension_reduction[[feat_i]][[obj_i]]
        feat_type = if(is.null(feat_type_names[[feat_i]])) paste0('feat_', feat_i) else feat_type_names[[feat_i]]
        name = if(is.null(obj_names[[obj_i]])) paste0('dimRed_', obj_i) else obj_names[[obj_i]]
        method = name # assume

        obj_list = append(obj_list, dr)
        feat_type_list = c(feat_type_list, feat_type)
        name_list = c(name_list, name)
        method_list = c(method_list, method)

      }
    }
    spat_unit_list = rep('cell', length(obj_list)) # assume

  } else if(depth(dimension_reduction) == 3L) { # ------------------------ 3 #

    spat_unit_names = names(dimension_reduction)
    for(unit_i in seq_along(dimension_reduction)) {

      feat_type_names = names(dimension_reduction[[unit_i]])
      for(feat_i in seq_along(dimension_reduction[[unit_i]])) {

        obj_names = names(dimension_reduction[[unit_i]][[feat_i]])
        for(obj_i in seq_along(dimension_reduction[[unit_i]][[feat_i]])) {

          dr = dimension_reduction[[unit_i]][[feat_i]][[obj_i]]
          feat_type = if(is.null(feat_type_names[[feat_i]])) paste0('feat_', feat_i) else feat_type_names[[feat_i]]
          spat_unit = if(is.null(spat_unit_names[[unit_i]])) paste0('unit_', unit_i) else spat_unit_names[[unit_i]]
          name = if(is.null(obj_names[[obj_i]])) paste0('dimRed_', obj_i) else obj_names[[obj_i]]
          method = name # assume

          obj_list = append(obj_list, dr)
          feat_type_list = c(feat_type_list, feat_type)
          spat_unit_list = c(spat_unit_list, spat_unit)
          name_list = c(name_list, name)
          method_list = c(method_list, method)

        }
      }
    }

  } else if(depth(dimension_reduction) == 4L) { # ------------------------ 4 #

    spat_unit_names = names(dimension_reduction)
    for(unit_i in seq_along(dimension_reduction)) {

      feat_type_names = names(dimension_reduction[[unit_i]])
      for(feat_i in seq_along(dimension_reduction[[unit_i]])) {

        method_names = names(dimension_reduction[[unit_i]][[feat_i]])
        for(method_i in seq_along(dimension_reduction[[unit_i]][[feat_i]])) {

          obj_names = names(dimension_reduction[[unit_i]][[feat_i]][[method_i]])
          for(obj_i in seq_along(dimension_reduction[[unit_i]][[feat_i]][[method_i]])) {

            dr = dimension_reduction[[unit_i]][[feat_i]][[method_i]][[obj_i]]
            feat_type = if(is.null(feat_type_names[[feat_i]])) paste0('feat_', feat_i) else feat_type_names[[feat_i]]
            spat_unit = if(is.null(spat_unit_names[[unit_i]])) paste0('unit_', unit_i) else spat_unit_names[[unit_i]]
            name = if(is.null(obj_names[[obj_i]])) paste0('dimRed_', obj_i) else obj_names[[obj_i]]
            method = if(is.null(method_names[[method_i]])) paste0('method_', method_i) else method_names[[method_i]]

            obj_list = append(obj_list, dr)
            feat_type_list = c(feat_type_list, feat_type)
            spat_unit_list = c(spat_unit_list, spat_unit)
            name_list = c(name_list, name)
            method_list = c(method_list, method)
          }
        }
      }
    }

  } else {
    stop(wrap_txt('Unexpected list depth', errWidth = TRUE))
  }


  if(length(obj_list > 0L)) {

    return_list = lapply(seq_along(obj_list), function(obj_i) {

      if(inherits(obj_list[[obj_i]], 'dimObj')) {
        return(obj_list[[obj_i]])
      } else {

        return(
          create_dim_obj(name = name_list[[obj_i]],
                         reduction = reduction,
                         coordinates = obj_list[[obj_i]],
                         reduction_method = method_list[[obj_i]],
                         spat_unit = spat_unit_list[[obj_i]],
                         feat_type = feat_type_list[[obj_i]],
                         provenance = if(is.null(provenance)) spat_unit_list[[obj_i]] else provenance, # assumed
                         misc = NULL)
        )
      }
    })
    return(return_list)

  } else {
    warning(wrap_txt('No objects found in dimension reduction input list'))
  }

}







#' @name check_dimension_reduction
#' @keywords internal
#' @noRd
check_dimension_reduction = function(gobject) {

  # check that all spatIDs of coordinates setequals with gobject cell_ID for the particular spat_unit
  avail_dr = list_dim_reductions(gobject = gobject)
  used_su = unique(avail_dr$spat_unit)

  if(!is.null(used_su)) {

    for(su_i in used_su) {
      IDs = spatIDs(gobject, spat_unit = su_i)

      su_dr = avail_dr[spat_unit == su_i,]
      lapply(seq(nrow(avail_dr)), function(obj_i) {
        dr_obj = get_dimReduction(gobject = gobject,
                                  spat_unit = su_i,
                                  feat_type = su_dr$feat_type[[obj_i]],
                                  reduction = su_dr$data_type[[obj_i]],
                                  reduction_method = su_dr$dim_type[[obj_i]],
                                  name = su_dr$name[[obj_i]],
                                  output = 'dimObj',
                                  set_defaults = FALSE)

        # if matrix has no IDs, regenerate from gobject IDs
        if(is.null(spatIDs(dr_obj))) {
          if(nrow(dr[]) == length(IDs)) {
            # if nrow of matrix and number of gobject cell_IDs for the spat unit
            # match, then try guessing then set data back to replace
            warning(wrap_txt('data_type:', su_dr$data_type[[obj_i]],
                             'spat_unit:', su_i,
                             'feat_type:', su_dr$feat_type[[obj_i]],
                             'dim_type:', su_dr$dim_type[[obj_i]],
                             'name:', su_dr$name[[obj_i]], '\n',
                             'Dimension reduction has no cell_IDs.
                             Guessing based on existing expression cell_IDs'))
            rownames(dr_obj[]) = IDs
            ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
            gobject = set_dimReduction(gobject = gobject,
                                       dimObject = dr_obj,
                                       set_defaults = FALSE,
                                       initialize = TRUE,
                                       verbose = FALSE)
            ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
          } else {
            # if number of values do NOT match, throw error
            stop(wrap_txt('data_type:', su_dr$data_type[[obj_i]],
                          'spat_unit:', su_i,
                          'feat_type:', su_dr$feat_type[[obj_i]],
                          'dim_type:', su_dr$dim_type[[obj_i]],
                          'name:', su_dr$name[[obj_i]], '\n',
                          'Dimension reduction has no cell_IDs.
                          Number of rows also does not match expression columns'))
          }
        }

        # if does not have all the IDS seen in the gobject
        if(!all(spatIDs(dr_obj) %in% IDs)) {
          warning(wrap_txt('data_type:', su_dr$data_type[[obj_i]],
                           'spat_unit:', su_i,
                           'feat_type:', su_dr$feat_type[[obj_i]],
                           'dim_type:', su_dr$dim_type[[obj_i]],
                           'name:', su_dr$name[[obj_i]], '\n',
                           'Dimension reduction coord names are not all found in gobject IDs'))
        }
          # print(paste(su_i, su_dr$name[[obj_i]])) # debug
      })
    }
  }
  return(gobject)
}






#### Giotto nearest network ####




#' @title Evaluate nearest networks
#' @name evaluate_nearest_networks
#' @description Evaluate nearest networks input
#' @keywords internal
#' @noRd
evaluate_nearest_networks = function(nn_network) {

  if(inherits(nn_network, 'nnNetObj')) {
    nn_network[] = evaluate_nearest_networks(nn_network = nn_network[])
    return(nn_network)
  } else if(inherits(nn_network, 'igraph')) {
    return(nn_network)
  } else if(inherits(nn_network, 'data.frame')) {
    if(all(c('from', 'to', 'weight', 'distance', 'shared', 'rank') %in% colnames(nn_network))) {
      nn_network = igraph::graph_from_data_frame(nn_network)
    } else {
      stop(wrap_txt('Unable to coerce data.frame type object to nnNetObj igraph
                    Needed columns: from, to, weight, distance, shared, rank'))
    }
  }

}






#' @title Read nearest neighbor networks
#' @name read_nearest_networks
#' @description read nearest network results from list
#' @keywords internal
#' @noRd
read_nearest_networks = function(nn_network,
                                 provenance = NULL) {

  if(is.null(nn_network)) {
    message('No nearest network results are provided')
    return(NULL)
  }

  # return directly if not list
  if(!inherits(nn_network, 'list')) {
    try_val = try(as.list(nn_network), silent = TRUE)
    if(inherits(try_val, 'try-error')) {
      return(list(nn_network))
    } else {
      nn_network = try_val
    }
  }


  # list reading
  obj_list = list()
  spat_unit_list = c()
  feat_type_list = c()
  method_list = c()
  name_list = c()

  # read nesting
  if(depth(nn_network) == 1L) {       # ------------------------ 1 #

    obj_names = names(nn_network)
    for(obj_i in seq_along(nn_network)) {

      nn = nn_network[[obj_i]]
      name = if(is.null(obj_names[[obj_i]])) paste0('nn_', obj_i) else obj_names[[obj_i]]
      method = name # assume

      obj_list = append(obj_list, nn)
      name_list = c(name_list, name)
      method_list = c(method_list, method)
    }
    spat_unit_list = rep('cell', length(obj_list)) # assume
    feat_type_list = rep('rna', length(obj_list)) # assume

  } else if(depth(nn_network) == 2L) { # ------------------------ 2 #

    feat_type_names = names(nn_network)
    for(feat_i in seq_along(nn_network)) {

      obj_names = names(nn_network[[feat_i]])
      for(obj_i in seq_along(nn_network[[feat_i]])) {

        nn = nn_network[[feat_i]][[obj_i]]
        feat_type = if(is.null(feat_type_names[[feat_i]])) paste0('feat_', feat_i) else feat_type_names[[feat_i]]
        name = if(is.null(obj_names[[obj_i]])) paste0('nn_', obj_i) else obj_names[[obj_i]]
        method = name # assume

        obj_list = append(obj_list, nn)
        name_list = c(name_list, name)
        feat_type_list = c(feat_type_list, feat_type)
        method_list = c(method_list, method)
      }
    }
    spat_unit_list = rep('cell', length(obj_list)) # assume

  } else if(depth(nn_network) == 3L) { # ------------------------ 3 #

    spat_unit_names = names(nn_network)
    for(unit_i in seq_along(nn_network)) {

      feat_type_names = names(nn_network[[unit_i]])
      for(feat_i in seq_along(nn_network[[unit_i]])) {

        obj_names = names(nn_network[[unit_i]][[feat_i]])
        for(obj_i in seq_along(nn_network[[unit_i]][[feat_i]])) {

          nn = nn_network[[unit_i]][[feat_i]][[obj_i]]
          spat_unit = if(is.null(spat_unit_names[[unit_i]])) paste0('unit_', unit_i) else spat_unit_names[[unit_i]]
          feat_type = if(is.null(feat_type_names[[feat_i]])) paste0('feat_', feat_i) else feat_type_names[[feat_i]]
          name = if(is.null(obj_names[[obj_i]])) paste0('nn_', obj_i) else obj_names[[obj_i]]
          method = name # assume

          obj_list = append(obj_list, nn)
          name_list = c(name_list, name)
          feat_type_list = c(feat_type_list, feat_type)
          spat_unit_list = c(spat_unit_list, spat_unit)
          method_list = c(method_list, method)
        }
      }
    }

  } else if(depth(nn_network) == 4L) { # ------------------------ 4 #

    spat_unit_names = names(nn_network)
    for(unit_i in seq_along(nn_network)) {

      feat_type_names = names(nn_network[[unit_i]])
      for(feat_i in seq_along(nn_network[[unit_i]])) {

        method_name = names(nn_network[[unit_i]][[feat_i]])
        for(method_i in seq_along(nn_network[[unit_i]][[feat_i]])) {

          obj_name = names(nn_network[[unit_i]][[feat_i]][[method_i]])
          for(obj_i in seq_along(nn_network[[unit_i]][[feat_i]][[method_i]])) {

            nn = nn_network[[unit_i]][[feat_i]][[method_i]][[obj_i]]
            spat_unit = if(is.null(spat_unit_names[[unit_i]])) paste0('unit_', unit_i) else spat_unit_names[[unit_i]]
            feat_type = if(is.null(feat_type_names[[feat_i]])) paste0('feat_', feat_i) else feat_type_names[[feat_i]]
            name = if(is.null(obj_names[[obj_i]])) paste0('nn_', obj_i) else obj_names[[obj_i]]
            method = if(is.null(method_names[[method_i]])) paste0('method_', method_i) else method_names[[method_i]]

            obj_list = append(obj_list, nn)
            spat_unit_list = c(spat_unit_list, spat_unit)
            feat_type_list = c(feat_type_list, feat_type)
            method_list = c(method_list, method)
            name_list = c(name_list, name)
          }
        }
      }
    }

  } else {
    stop(wrap_txt('Unexpected list depth', errWidth = TRUE))
  }


  if(length(obj_list > 0L)) {

    return_list = lapply(seq_along(obj_list), function(obj_i) {

      if(inherits(obj_list[[obj_i]], 'nnNetObj')) {
        return(obj_list[[obj_i]])
      } else {

        return(
          create_nn_net_obj(name = name_list[[obj_i]],
                            nn_type = method_list[[obj_i]],
                            igraph = obj_list[[obj_i]],
                            spat_unit = spat_unit_list[[obj_i]],
                            feat_type = feat_type_list[[obj_i]],
                            provenance = provenance,
                            misc = NULL)
        )
      }
    })
    return(return_list)

  } else {
    warning(wrap_txt('No objects found in nearest neighbor network input list'))
  }

}





#' @keywords internal
#' @noRd
check_nearest_networks = function(gobject) {

  avail_nn = list_nearest_networks(gobject = gobject)
  used_su = unique(avail_nn$spat_unit)

  if(!is.null(used_su)) {

    for(su_i in used_su) {
      IDs = spatIDs(gobject, spat_unit = su_i)

      su_nn = avail_nn[spat_unit == su_i,]
      lapply(seq(nrow(avail_nn)), function(obj_i) {
        nn_obj = get_NearestNetwork(gobject = gobject,
                                    spat_unit = su_i,
                                    feat_type = su_nn$feat_type[[obj_i]],
                                    nn_network_to_use = su_nn$nn_type[[obj_i]],
                                    network_name = su_nn$name[[obj_i]],
                                    output = 'nnNetObj',
                                    set_defaults = FALSE)
        if(!all(spatIDs(nn_obj) %in% IDs)) {
          warning(wrap_txt('spat_unit:', su_i,
                           'feat_type:', su_nn$feat_type[[obj_i]],
                           'nn_type:', su_nn$nn_type[[obj_i]],
                           'name:', su_nn$name[[obj_i]], '\n',
                           'Nearest network vertex names are not all found in gobject IDs'))
        }
        # print(paste(su_i, su_nn$name[[obj_i]])) # debug
      })
    }
  }

}








#### Giotto spatial info ####


#' @title Evaluate spatial info
#' @name evaluate_spatial_info
#' @description Evaluate spatial information input into a SpatVector for
#' giottoPolygon creation
#' @param spatial_info spatial information to evaluate
#' @param skip_eval_dfr (default FALSE) skip evaluation of data.frame like input
#' @param copy_dt (default TRUE) if segmdfr is provided as dt, this determines
#' whether a copy is made
#' @param cores how many cores to use
#' @param verbose be verbose
#' @return SpatVector
#' @keywords internal
evaluate_spatial_info = function(spatial_info,
                                 skip_eval_dfr = FALSE,
                                 copy_dt = TRUE,
                                 cores = determine_cores(),
                                 verbose = TRUE) {

  # data.table vars
  geom = NULL

  ## 1. load or read spatial information data ##
  ## 1.1 read from file
  if(inherits(spatial_info, 'character')) {
    if(!file.exists(path.expand(spatial_info))) stop('path to spatial information does not exist')
    spatial_info = data.table::fread(input = spatial_info, nThread = cores)

  ## 1.2 data.frame-like input
  } else if (inherits(spatial_info, 'data.table')) {
    if(isTRUE(copy_dt)) spatial_info = data.table::copy(spatial_info)
  } else if(inherits(spatial_info, 'data.frame')) {
    spatial_info = data.table::setDT(spatial_info)

  ## 1.3 SpatVector input
  } else if(inherits(spatial_info, 'SpatVector')) {
    if(!'poly_ID' %in% names(spatial_info)) {
      stop(wrap_txt('SpatVector input requires the poly_ID attribute'))
    }
    return(spatial_info)

  ## 1.4 Other inputs
  } else {
    spatial_info = try(data.table::as.data.table(spatial_info), silent = TRUE)
    if(inherits(spatial_info, 'try-error')) {
      stop(wrap_txt('If spatial information is provided then it needs to be a',
                    'file path or a data.frame-like object',
                    errWidth = TRUE))
    }
  }
  # TODO add compatibility for HDF5s from specific outputs?


  # 2. data.frame info evaluation
  if(!isTRUE(skip_eval_dfr)) {
    spatial_info = evaluate_gpoly_dfr(input_dt = spatial_info,
                                      verbose = verbose)
  }


  # 3. add other necessary cols for the input data.table
  nr_of_cells_vec = seq_along(unique(spatial_info$poly_ID))
  names(nr_of_cells_vec) = unique(spatial_info$poly_ID)
  new_vec = nr_of_cells_vec[as.character(spatial_info$poly_ID)]
  spatial_info[, geom := new_vec]

  spatial_info[, c('part', 'hole') := list(1, 0)]
  spatial_info = spatial_info[, c('geom', 'part', 'x', 'y', 'hole', 'poly_ID'), with = FALSE]


  # 4. create spatvector
  spatial_info = dt_to_spatVector_polygon(spatial_info,
                                          include_values = TRUE)


  return(spatial_info)


  # ## 2. check and name columns ##
  # nr_cols = ncol(spatial_info)
  #
  # if(nr_cols < 4) stop('Spatial information needs to have at least 4 columns: \n',
  #                      'x, y, (z) information columns \n',
  #                      'cell ID and polygon point column \n')
  #
  # column_classes = lapply(spatial_info, FUN = class)
  #
  # # 3D or 2D data
  # if(all(column_classes[1:3] == 'numeric')) {
  #   colnames(spatial_info)[1:5] = c('sdimx', 'sdimy', 'sdimz', 'cell_ID', 'point')
  # } else if(all(column_classes[1:2] == 'numeric')){
  #   colnames(spatial_info)[1:4] = c('sdimx', 'sdimy', 'cell_ID', 'point')
  # } else {
  #   stop('First 3 or 2 columns need to be numeric for 3D and 2D data respectively')
  # }



  # ## 3. check cell ID ##
  # spatial_locs_cell_IDs = spatial_locs[['cell_ID']]
  #
  # spatial_info_cell_IDs = spatial_info[['cell_ID']]
  #
  # if(all(spatial_info_cell_IDs %in% spatial_locs_cell_IDs)) {
  #   return(spatial_info)
  # } else {
  #   stop('cell IDs in spatial information are missing in the spatial locations slot')
  # }

}






#' @name check_spatial_info
#' @keywords internal
#' @noRd
check_spatial_info = function(gobject) {

  avail_sinfo = list_spatial_info(gobject)
  avail_slocs = list_spatial_locations(gobject)

  common_su = intersect(avail_sinfo$spat_info, avail_slocs$spat_unit)

  # If there are any shared spatial units, match IDs
  if(length(common_su) != 0) {

    for(su_i in common_su) {

      # get spat_info
      sinfo = get_polygon_info(gobject = gobject,
                               polygon_name = su_i,
                               return_giottoPolygon = TRUE)

      # get spatlocs
      su_sloc = avail_slocs[spat_unit == su_i]
      lapply(seq(nrow(su_sloc)), function(obj_i) {
        spatlocs = get_spatial_locations(gobject = gobject,
                                         spat_unit = su_i,
                                         spat_loc_name = su_sloc$name[[obj_i]],
                                         output = 'spatLocsObj',
                                         set_defaults = FALSE,
                                         copy_obj = FALSE)
        if(!all(spatIDs(spatlocs) %in% spatIDs(sinfo))) {
          warning(wrap_txt('spat_unit:', su_i,
                           'spatloc name:', su_sloc[[obj_i]], '\n',
                           'cell IDs in spatial locations are missing from spatial polygon info'))
        }
        # print(paste(su_i, su_sloc$name[[obj_i]])) # debug
      })
    }
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
                              cores = determine_cores(),
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
                               cores = determine_cores(),
                               verbose = TRUE) {

  # create minimum giotto
  gobject = new('giotto',
                expression_feat = expression_feat,
                offset_file = offset_file,
                instructions = instructions,
                OS_platform = .Platform[['OS.type']])

  # gobject = giotto(expression = list(),
  #                  expression_feat = expression_feat,
  #                  spatial_locs = spatial_locs,
  #                  spatial_info = NULL,
  #                  cell_metadata = cell_metadata,
  #                  feat_metadata = feat_metadata,
  #                  feat_info = feat_info,
  #                  cell_ID = NULL,
  #                  feat_ID = NULL,
  #                  spatial_network = NULL,
  #                  spatial_grid = NULL,
  #                  spatial_enrichment = NULL,
  #                  dimension_reduction = NULL,
  #                  nn_network = NULL,
  #                  images = NULL,
  #                  largeImages = NULL,
  #                  parameters = NULL,
  #                  offset_file = offset_file,
  #                  instructions = instructions,
  #                  OS_platform = .Platform[['OS.type']],
  #                  join_info = NULL)


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
    wrap_msg("Consider to install these (optional) packages to run all possible",
             "Giotto commands for spatial analyses: ", extra_not_installed_packages)
    wrap_msg("Giotto does not automatically install all these packages as they",
             "are not absolutely required and this reduces the number of dependencies")
  }


  ## if cores is not set, then set number of cores automatically, but with limit of 10
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
    for(expr_i in seq_along(expression_data)) {
      gobject = set_expression_values(gobject = gobject,
                                      values = expression_data[[expr_i]],
                                      set_defaults = FALSE)
    }

    print('test 1')


    # Set up gobject cell_ID and feat_ID slots based on expression matrices
    gobject = init_cell_and_feat_IDs(gobject)

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

    # parse spatial_loc param input for any spat_unit/name info
    spatial_location_data = read_spatial_location_data(spat_loc_list = spatial_locs,
                                                       cores = cores,
                                                       verbose = verbose)

    # set spatial location data
    for(spatloc_i in seq_along(spatial_location_data)) {
      gobject = set_spatial_locations(gobject, spatlocs = spatial_location_data[[spatloc_i]])
    }

    # slot(gobject, 'spatial_locs') = spatial_location_data

    # 1. ensure spatial locations and expression matrices have the same cell IDs
    # 2. give cell IDs if not provided
    check_spatial_location_data(gobject) # modifies by reference

  } else {

    if(verbose == TRUE) warning(wrap_txt('No spatial locations have been provided, dummy locations will be created'))

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

      dummySpatLocObj = create_spat_locs_obj(name = 'raw',
                                             coordinates = spatial_locs,
                                             spat_unit = spat_unit)

      ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
      set_spatial_locations(gobject, spatlocs = dummySpatLocObj)
      ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

    }
  }

  if(verbose) message('finished spatial location data')


  ## spatial info ##
  ## ------------ ##
  ## place to store segmentation info in polygon format style

  if(!is.null(spatial_info)) {
    gobject = addGiottoPolygons(gobject = gobject,
                                gpolygons = spatial_info)
  }



  ## cell metadata ##
  ## ------------- ##
  if(is.null(cell_metadata)) {
    # initialize metadata
    gobject = init_cell_metadata(gobject)
  } else {
    gobject = set_cell_metadata(gobject = gobject,
                                metadata = cell_metadata)
  }

  if(verbose) message('finished cell metadata')

  ## feat metadata ##
  ## ------------- ##
  if(is.null(feat_metadata)) {
    # initialize metadata
    gobject = init_feat_metadata(gobject)
  } else {
    gobject = set_feature_metadata(gobject = gobject,
                                   metadata = feat_metadata)
  }

  if(verbose) message('finished feature metadata')

  ## feature info ##
  ## ------------ ##
  ## place to store individual feature info
  if(!is.null(feat_info)) {
    gobject = addGiottoPoints(gobject = gobject,
                              gpoints = feat_info)
  }




  ### OPTIONAL:
  ## spatial network
  spatial_network_list = read_spatial_networks(spatial_network)
  for(sn_i in seq_along(spatial_network_list)) {
    gobject = set_spatialNetwork(gobject = gobject,
                                 spatial_network = spatial_network_list[[sn_i]],
                                 set_defaults = FALSE)
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
  spatial_enrichment = read_spatial_enrichment(spatial_enrichment = spatial_enrichment)
  for(se_i in seq_along(spatial_enrichment)) {
    gobject = set_spatial_enrichment(gobject = gobject,
                                     spatenrichment = spatial_enrichment[[se_i]],
                                     set_defaults = FALSE)
  }




  ## dimension reduction
  dimension_reduction = read_dimension_reduction(dimension_reduction = dimension_reduction)
  for(dr_i in seq_along(dimension_reduction)) {
    gobject = set_dimReduction(gobject = gobject,
                               dimObject = dimension_reduction[[dr_i]],
                               set_defaults = FALSE)
  }


  # NN network
  nn_network = read_nearest_networks(nn_network = nn_network)
  for(nn_i in seq_along(nn_network)) {
    gobject = set_NearestNetwork(gobject = gobject,
                                 nn_network = nn_network[[nn_i]],
                                 set_defaults = FALSE)
  }




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
    visium_dir = path.expand(visium_dir)
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
#' @param largeImages_list_params image params when loading largeImages as list
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

  # data.table vars
  poly_ID = cell_ID = feat_ID = x = y = NULL

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

      # gobject = set_spatial_locations(gobject = gobject,
      #                                 spat_unit = polygon_info,
      #                                 spat_loc_name = 'raw',
      #                                 spatlocs = centroidsDT_loc,
      #                                 verbose = FALSE)

      locsObj = create_spat_locs_obj(name = 'raw',
                                     coordinates = centroidsDT_loc,
                                     spat_unit = polygon_info,
                                     provenance = polygon_info,
                                     misc = NULL)

      ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
      gobject = set_spatial_locations(gobject, spatlocs = locsObj)
      ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

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
      # initialize metadata
      gobject = init_cell_metadata(gobject)
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
      # initialize metadata
      gobject = init_feat_metadata(gobject)
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
  ## spatial network - external input
  if(!is.null(spatial_network)) {
    if(is.null(spatial_network_name) | length(spatial_network) != length(spatial_network_name)) {
      stop('\n each spatial network must be given a unique name \n')
    } else {

      for(network_i in 1:length(spatial_network)) {

        networkname = spatial_network_name[[network_i]]
        network     = spatial_network[[network_i]]

        if(any(c('data.frame', 'data.table') %in% class(network))) {
          if(all(c('to', 'from', 'weight', 'sdimx_begin', 'sdimy_begin', 'sdimx_end', 'sdimy_end') %in% colnames(network))) {

            # create spatialNetworkObj from data.table
            network = data.table::setDT(network)
            # most info will be missing
            warning('spatial_network ', network_i, ' provided as data.table/frame object. Provenance and spat_unit will be assumed: "', names(slot(gobject, 'spatial_info'))[[1]], '"\n')
            spatial_network_Obj = create_spat_net_obj(name = networkname,
                                                      networkDT = network,
                                                      spat_unit = names(slot(gobject, 'spatial_info'))[[1]],
                                                      provenance = names(slot(gobject, 'spatial_info'))[[1]]) # assumed

            ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
            gobject = set_spatialNetwork(gobject, spatial_network = spatial_network_Obj)
            ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

          } else {
            stop('\n network ', networkname, ' does not have all necessary column names, see details\n')
          }
        } else if(inherits(network, 'spatialNetworkObj')) {

          ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
          gobject = set_spatialNetwork(gobject, spatial_network = network)
          ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

        } else {
          stop('\n network ', networkname, ' is not a, spatialNetworkObj, data.frame, or data.table\n')
        }
      }
    }
  }


  ## spatial grid - external input
  if(!is.null(spatial_grid)) {
    if(is.null(spatial_grid_name) | length(spatial_grid) != length(spatial_grid_name)) {
      stop('\n each spatial grid must be given a unique name \n')
    } else {

      for(grid_i in 1:length(spatial_grid)) {

        gridname = spatial_grid_name[[grid_i]]
        grid     = spatial_grid[[grid_i]]

        if(inherits(grid, c('data.table', 'data.frame'))) {
          if(all(c('x_start', 'y_start', 'x_end', 'y_end', 'gr_name') %in% colnames(grid))) {
            if(!inherits(grid, 'data.table')) grid = data.table::setDT(grid)
            # Assume first spat_info and first expression_feat as spat_unit/provenance and feat_type respectively
            warning('spatial_grid ', grid_i, ' provided as data.table/frame object. Provenance and spat_unit will be assumed: "', names(slot(gobject, 'spatial_info'))[[1]], '"\n')
            # most info will be missing
            grid = create_spat_grid_obj(name = gridname,
                                        gridDT = grid,
                                        spat_unit = names(slot(gobject, 'spatial_info'))[[1]],
                                        provenance = names(slot(gobject, 'spatial_info'))[[1]],
                                        feat_type = expression_feat[[1]])

            ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
            gobject = set_spatialGrid(gobject = gobject, spatial_grid = grid)
            ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

          } else {
            stop('\n grid ', gridname, ' does not have all necessary column names, see details \n')
          }
        } else if(inherits(grid, 'spatialGridObj')) {

          ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
          gobject = set_spatialGrid(gobject, spatial_grid = grid)
          ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

        } else {
          stop('\n grid ', gridname, ' is not a spatialGridObj, data.frame, or data.table \n')
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

        # glargeImage = createGiottoLargeImage(raster_object = im,
        #                                      negative_y = FALSE,
        #                                      name = im_name)

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




# ********************************************************************** #
# Spatial-method specific object creation convenience functions moved to
# convenience_create.R
# ********************************************************************** #



#### logging of giotto functions ####


# Determine the name of the function n levels above the current evaluation frame,
# where n is toplevel - 1
#' @keywords internal
#' @noRd
get_prev_fname = function(toplevel = 3L, verbose = FALSE) {
  as.character(sys.call(-toplevel)[[1]])
}





#' @title Log args used
#' @name get_args
#' @keywords internal
#' @noRd
get_args <- function(toplevel = 2L, verbose = FALSE) {

  nframes = sys.nframe()

  if(isTRUE(verbose)) {
    cat('\n number of frames: ')
    print(nframes)
    cat('\n')
  }


  cl = sys.call(-toplevel)

  if(isTRUE(verbose)) {
    cat('\n system call: ')
    print(cl)
    cat('\n')
  }


  # function name
  fname = as.character(cl[[1]])

  if(length(fname) > 1) {
    fname = fname[[3]]
  }

  if(isTRUE(verbose)) {
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
#' @noRd
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



#### provenance ####

#' @title Check provenance info matches across list of S4 subobjects
#' @name prov_match
#' @param ... list of s4 subobjects to match provenance
#' @param verbose print discovered provenance info
#' @keywords internal
#' @return NULL
prov_match = function(..., verbose = FALSE) {

  if(length(match.call()) == 2) {
    if(isS4(...)) stop('more than one object must be provided')
    s4_list = as.list(...)
  } else if(length(match.call()) > 2) {
    s4_list = list(...)
  }

  if(!all(unlist(lapply(s4_list, inherits, 'provData'))))
    stop('Not all items to test contain provenance information\n')

  prov_list = lapply(s4_list, prov)

  if(isTRUE(verbose)) print(unlist(prov_list))

  length(unique(prov_list)) == 1

}




#### joining giotto object ####

#' @title join_expression_matrices
#' @name join_expression_matrices
#' @keywords internal
#' @noRd
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
#' @noRd
join_spatlocs = function(dt_list) {

  final_list = do.call('rbind', dt_list) # breaks DT reference
  return(final_list)
}

#' @title join_cell_meta
#' @name join_cell_meta
#' @keywords internal
#' @noRd
join_cell_meta = function(dt_list) {

  final_list = do.call('rbind', dt_list)
  return(final_list)

}

#' @title join_feat_meta
#' @name join_feat_meta
#' @keywords internal
#' @noRd
join_feat_meta = function(dt_list) {

  feat_ID = NULL

  comb_meta = do.call('rbind', c(dt_list, fill = TRUE))
  comb_meta = unique(comb_meta)

  # find feat_IDs with multiple metadata versions
  dup_feats = unique(comb_meta[duplicated(feat_ID), feat_ID])

  # pick first entry of duplicates
  comb_meta = comb_meta[!duplicated(feat_ID)]

  if(length(dup_feats) > 0) {
    warning(wrap_txt('feature metadata: multiple versions of metadata for:\n',
                     dup_feats,
                     '\n First entry will be selected for joined object.'))
  }

  return(comb_meta)

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
# @param dry_run (experimental) Works best for join_method 'shift' or 'no_change'.
#' @param verbose be verbose
#' Preview where each gobject will be in space with bounding polygons
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
                             # dry_run = FALSE,
                             verbose = TRUE) {

  # define for data.table
  sdimz = cell_ID = sdimx = sdimy = name = NULL

  n_gobjects = length(gobject_list)

  ## check general input params
  if(n_gobjects == 0) stop(wrap_txt('A list of Giotto objects to be joined must be provided.'))
  if(n_gobjects == 1) stop(wrap_txt('Only one gobject provided in gobject_list.'))
  if(!is.vector(gobject_names) | !is.character(gobject_names)) stop(wrap_txt('gobject_names need to be a vector with unique names for the giotto objects'))
  if(n_gobjects != length(gobject_names)) stop(wrap_txt('each giotto object in the list needs to have a unique (short) name'))
  if(is.null(join_method)) wrap_msg('No join_method given. Defaulting to "shift"')


  ## determine join method
  join_method = match.arg(arg = join_method, choices = c('shift', 'z_stack', 'no_change'))
  if(isTRUE(verbose)) message('Join method:', join_method)


  # **** For shift workflow ****
  if(join_method == 'shift') {
    # Make sure enough x_shift and y_shift values are given to cover all gobjects
    if(!is.null(x_shift)) if(length(x_shift) != n_gobjects) stop(wrap_txt('A numeric vector with an x_shift value for each gobject in gobject_list must be given.\n'))
    if(!is.null(y_shift)) if(length(y_shift) != n_gobjects) stop(wrap_txt('A numeric vector with a y_shift value for each gobject in gobject_list must be given.\n'))

    # Set defaults if no shift params are given
    if(is.null(x_shift) & is.null(y_shift) & is.null(x_padding) & is.null(y_padding)) {
      wrap_msg('No xy shift or specific padding values given. Using defaults: x_padding = 1000')
      x_padding = 1000
    }
    # Assign default padding values if NULL
    if(is.null(x_padding)) x_padding = 0
    if(is.null(y_padding)) y_padding = 0
  }



  # **** For no_change workflow ****
  if(join_method == 'no_change') {
    join_method = 'shift'
    x_shift = rep(0, n_gobjects)
    y_shift = rep(0, n_gobjects)
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


  # TODO # **** dry run ****
  # if(isTRUE(dry_run)) {
  #   if(isTRUE(verbose)) wrap_msg('dry_run = TRUE:
  #                                Spatial preview of join operation.')
  #   # Detect sources of bounds info
  #   avail_bound_info = list(
  #     avail_img = lapply(gobject_list, list_images),
  #     avail_spat_info = lapply(gobject_list, list_spatial_info),
  #     avail_feat_info = lapply(gobject_list, list_feature_info),
  #     avail_spatlocs = lapply(gobject_list, list_spatial_locations)
  #   )
  #   avail_bound = lapply(avail_bound_info, function(avail) {
  #     isTRUE(!is.null(unlist(avail)))
  #   })
  #   if(is.null(unlist(avail_bound))) stop(wrap_txt('dry_run error: No shared sources of bounds info
  #                                              Previewing from heterogenous sources not yet implemented'))
  #
  #   bound_to_use = avail_bound_info[names(avail_bound_info)[1]]
  #
  #   # get bound info
  #   if(names(bound_to_use) == 'avail_img') {
  #
  #   }
  #   if(names(bound_to_use) == 'avail_spat_info') {
  #
  #   }
  #   if(names(bound_to_use) == 'avail_feat_info') {
  #
  #   }
  #   if(names(bound_to_use) == 'avail_spatlocs') {
  #
  #   }
  #
  # }


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
  if(verbose == TRUE) wrap_msg('start updating objects')

  all_feat_ID_list = list()
  all_cell_ID_list = list()
  all_image_list = list()
  all_largeImage_list = list()

  xshift_list = list()
  yshift_list = list()

  all_spatinfo_list = list()


  if(verbose) wrap_msg('A) Update giotto Objects \n')

  for(gobj_i in 1:n_gobjects) {

    gobj = gobject_list[[gobj_i]]
    gname = gobject_names[[gobj_i]]


    ## 0. update cell ID and feat ID
    if(verbose) wrap_msg('0. Update cell and feature IDs \n')

    for(spat_unit in names(gobj@cell_ID)) {
      old_cell_ID = get_cell_id(gobject = gobj, spat_unit = spat_unit)
      new_cell_ID = paste0(gname, '-', old_cell_ID)
      all_cell_ID_list[[spat_unit]][[gobj_i]] = new_cell_ID
      gobj = set_cell_id(gobject = gobj,
                         spat_unit = spat_unit,
                         cell_IDs = new_cell_ID,
                         set_defaults = FALSE)
    }



    for(feat_type in names(gobj@feat_ID)) {
      all_feat_ID_list[[feat_type]][[gobj_i]] = get_feat_id(gobject = gobj, feat_type = feat_type)
    }



    ## 1. update expression and all feature IDs
    if(verbose) wrap_msg('1. Update expression IDs \n')

    # provide unique cell ID name
    for(spat_unit in names(gobj@expression)) {

      for(feat_type in names(gobj@expression[[spat_unit]])) {

        for(matr in names(gobj@expression[[spat_unit]][[feat_type]])) {
          colnames(gobj@expression[[spat_unit]][[feat_type]][[matr]][]) = gobj@cell_ID[[spat_unit]]
        }

        #all_feat_ID_list[[feat_type]][[gobj_i]] = gobj@feat_ID[[feat_type]]
      }

    }




    ## 2. update images
    # change individual names
    if(verbose) wrap_msg('2. Update images \n')

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

          # save extent info
          gobj@largeImages[[imname]]@extent = terra::ext(gobj@largeImages[[imname]]@raster_object)[]
        }

        all_largeImage_list[[imname]] = gobj@largeImages[[imname]]

      }


    }




    ## 3. update spatial location
    if(verbose) wrap_msg('3. Update spatial locations \n')

    # add padding to x-axis
    # update cell ID

    # If no images were present
    if(length(xshift_list) == 0) xshift_list = ((seq_along(gobject_list) - 1) * x_padding)

    available_locs = list_spatial_locations(gobj)

    for(spat_unit_i in available_locs[['spat_unit']]) {

      for(locs_i in available_locs[spat_unit == spat_unit_i, name]) {
        spat_obj = get_spatial_locations(gobj,
                                         spat_unit = spat_unit_i,
                                         spat_loc_name = locs_i,
                                         output = 'spatLocsObj',
                                         copy_obj = TRUE,
                                         set_defaults = FALSE)
        myspatlocs = slot(spat_obj, 'coordinates')

        if(join_method == 'z_stack') {
          myspatlocs[, sdimz := z_vals[gobj_i]]
          myspatlocs[, cell_ID := get_cell_id(gobj, spat_unit = spat_unit_i) ]
          myspatlocs = myspatlocs[,.(sdimx, sdimy, sdimz, cell_ID)]

        } else if(join_method == 'shift') {


          # shift for x-axis
          if(is.null(x_shift)) {
            add_to_x = xshift_list[[gobj_i]]
          } else {
            x_shift_i = x_shift[[gobj_i]]
            add_to_x = x_shift_i + (x_padding * (gobj_i - 1))
          }

          if(verbose) cat('Spatial locations: for ',locs_i, ' add_to_x = ', add_to_x, '\n')


          myspatlocs[, sdimx := sdimx + add_to_x]
          myspatlocs[, cell_ID := get_cell_id(gobj, spat_unit = spat_unit_i) ]
        }

        # shift for y-axis
        if(!is.null(y_shift)) {
          y_shift_i = y_shift[[gobj_i]]
          add_to_y = y_shift_i + (y_padding * (gobj_i - 1))

          if(verbose) cat('Spatial locations: for ',locs_i, ' add_to_y = ', add_to_y, '\n')

          myspatlocs[, sdimy := sdimy + add_to_y]
        }

        slot(spat_obj, 'coordinates') = myspatlocs

        ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
        gobj = set_spatial_locations(gobject = gobj, spatlocs = spat_obj)
        ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

      }
    }





    ## 4. cell metadata
    # * (feat metadata happens during joined object creation)
    # rbind metadata
    # create capture area specific names
    if(verbose) wrap_msg('4. rbind cell metadata')

    for(spat_unit in names(gobj@cell_metadata)) {
      for(feat_type in names(gobj@cell_metadata[[spat_unit]])) {
        gobj@cell_metadata[[spat_unit]][[feat_type]]@metaDT[['cell_ID']] = gobj@cell_ID[[spat_unit]]
        gobj@cell_metadata[[spat_unit]][[feat_type]]@metaDT[['list_ID']] = gname
      }
    }



    ## 5. prepare spatial information
    if(verbose) wrap_msg('5. prepare spatial information')

    spatinfo_vector = vector()
    for(spat_info in names(gobj@spatial_info)) {

      # spatVector
      if(verbose) wrap_msg('-- 5.1. spatVector')
      poly_ids = gobj@spatial_info[[spat_info]]@spatVector$poly_ID
      gobj@spatial_info[[spat_info]]@spatVector$poly_ID = paste0(gname,'-',poly_ids)

      # spatVectorCentroids
      if(verbose) wrap_msg('-- 5.2. spatVectorCentroids')
      if(!is.null(gobj@spatial_info[[spat_info]]@spatVectorCentroids)) {
        poly_ids = gobj@spatial_info[[spat_info]]@spatVectorCentroids$poly_ID
        gobj@spatial_info[[spat_info]]@spatVectorCentroids$poly_ID = paste0(gname,'-',poly_ids)
      }


      # overlaps??
      if(verbose) wrap_msg('-- 5.3. overlaps')
      # TODO


      if(join_method == 'shift') {

        if(verbose) wrap_msg('-- 5.4. shift')

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

        gobj@spatial_info[[spat_info]]@spatVector = terra::shift(x = gobj@spatial_info[[spat_info]]@spatVector,
                                                                 dx = add_to_x,
                                                                 dy = add_to_y)
        if(!is.null(gobj@spatial_info[[spat_info]]@spatVectorCentroids)) {
          gobj@spatial_info[[spat_info]]@spatVectorCentroids = terra::shift(x = gobj@spatial_info[[spat_info]]@spatVectorCentroids, dx = add_to_x, dy = add_to_y)
        }


      }

      spatinfo_vector = c(spatinfo_vector, spat_info)
      all_spatinfo_list[[gobj_i]] = spatinfo_vector
    }


    ## 6. prepare feature information
    if(verbose) wrap_msg('6. prepare feature information \n')

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
  if(verbose) wrap_msg('B) Prepare to create new Giotto object \n')

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





  ## 3. merge updated data  ##
  ## ------------------------ ##
  if(verbose) wrap_msg('C) Merge updated data \n')

  first_obj = updated_object_list[[1]]

  if(verbose) wrap_msg('1. cell and feature IDs \n')
  ## cell IDs
  for(spat_unit in names(all_cell_ID_list)) {
    combined_cell_ID = unlist(all_cell_ID_list[[spat_unit]])
    comb_gobject@cell_ID[[spat_unit]] = combined_cell_ID
  }

  ## feat IDs
  for(feat_type in names(all_feat_ID_list)) {
    combined_feat_ID = unique(unlist(all_feat_ID_list[[feat_type]]))
    comb_gobject@feat_ID[[feat_type]] = combined_feat_ID
  }



  ## expression and feat IDs
  ## if no expression matrices are provided, then just combine all feature IDs
  if(verbose) wrap_msg('2. expression data \n')

  avail_expr = list_expression(gobject = first_obj)

  if(is.null(avail_expr)) {

    ## feat IDS
    for(feat in first_features) {
      combined_feat_ID = unique(unlist(all_feat_ID_list[[feat]]))
      ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
      comb_gobject = set_feat_id(gobject = comb_gobject,
                                 feat_type = feat,
                                 feat_IDs = combined_feat_ID,
                                 set_defaults = FALSE)
      ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
    }

    # Moved de novo feature metadata generation as a catch to the end of the fxn
    # Done through init_feat_meta()

    # S4_feat_metadata = create_feat_meta_obj(spat_unit = spat_unit,
    #                                         feat_type = feat_type,
    #                                         metaDT = data.table::data.table(feat_ID = combined_feat_ID))

    # comb_gobject = set_feature_metadata(gobject = comb_gobject,
    #                                     S4_feat_metadata,
    #                                     set_defaults = FALSE)


  } else {

    for(exprObj_i in seq(nrow(avail_expr))) {

      expr_list = lapply(updated_object_list, function(gobj) {
        get_expression_values(gobject = gobj,
                              spat_unit = avail_expr$spat_unit[[exprObj_i]],
                              feat_type = avail_expr$feat_type[[exprObj_i]],
                              values = avail_expr$name[[exprObj_i]],
                              output = 'exprObj',
                              set_defaults = FALSE)
      })

      if(!prov_match(expr_list)) warning(wrap_txt('expression: provenance mismatch'))

      combmat = join_expression_matrices(matrix_list = lapply(expr_list, function(expr) expr[]))
      expr_list[[1]][] = combmat[['matrix']]

      ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
      comb_gobject = set_expression_values(gobject = comb_gobject,
                                           values = expr_list[[1]],
                                           set_defaults = FALSE)

      comb_gobject = set_feat_id(gobject = comb_gobject,
                                 feat_type = avail_expr$feat_type[[exprObj_i]],
                                 feat_IDs = combmat[['sort_all_feats']],
                                 set_defaults = FALSE)
      ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

      # Moved de novo feat metadata generation to end of fxn as a catch

    }

  }




  ## spatial locations
  if(verbose) wrap_msg('3. spatial locations \n')

  available_locs = list_spatial_locations(first_obj)

  for(slObj_i in seq(nrow(available_locs))) {

    sl_list = lapply(updated_object_list, function(gobj) {
      get_spatial_locations(gobject = gobj,
                            spat_unit = available_locs$spat_unit[[slObj_i]],
                            spat_loc_name = available_locs$name[[slObj_i]],
                            output = 'spatLocsObj',
                            copy_obj = FALSE)
    })

    if(!prov_match(sl_list)) warning(wrap_txt('spatial locations: provenance mismatch'))

    combspatlocs = join_spatlocs(dt_list = lapply(sl_list, function(sl) sl[]))
    sl_list[[1]][] = combspatlocs

    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
    comb_gobject = set_spatial_locations(comb_gobject,
                                         spatlocs = sl_list[[1]],
                                         set_defaults = FALSE)
    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  }





  ## cell metadata
  if(isTRUE(verbose)) wrap_msg('4. cell metadata \n')

  for(spat_unit in names(first_obj@cell_metadata)) {

     for(feat_type in names(first_obj@cell_metadata[[spat_unit]])) {

      savelist = list()
      for(gobj_i in seq_along(updated_object_list)) {
        cellmeta = updated_object_list[[gobj_i]]@cell_metadata[[spat_unit]][[feat_type]][]
        savelist[[gobj_i]] = cellmeta
      }
      combcellmeta = join_cell_meta(dt_list = savelist)

      S4_cell_meta = get_cell_metadata(gobject = first_obj,
                                       spat_unit = spat_unit,
                                       feat_type = feat_type,
                                       copy_obj = TRUE,
                                       set_defaults = FALSE,
                                       output = 'cellMetaObj')
      S4_cell_meta[] = combcellmeta

      ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
      comb_gobject = set_cell_metadata(gobject = comb_gobject,
                                       metadata = S4_cell_meta,
                                       set_defaults = FALSE)
      ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

    }
  }


  ## feat metadata
  # check if any already exists in first_obj
  # otherwise, skip and generate feat_metadata de novo at end
  avail_featmeta = list_feat_metadata(gobject = first_obj)
  if(!is.null(avail_featmeta)) {
    if(isTRUE(verbose)) message('   feature metadata \n')
    for(fmObj_i in seq(nrow(avail_featmeta))) {

      fm_list = lapply(updated_object_list, function(gobj) {
        get_feature_metadata(gobject = gobj,
                             spat_unit = avail_featmeta$spat_unit[[fmObj_i]],
                             feat_type = avail_featmeta$feat_type[[fmObj_i]],
                             output = 'featMetaObj',
                             copy_obj = TRUE)
      })

      if(!prov_match(fm_list)) warning(wrap_txt('feature metadata: provenance mismatch'))

      comb_fm = join_feat_meta(dt_list = lapply(fm_list, function(fm) fm[]))
      fm_list[[1]][] = comb_fm

      ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
      comb_gobject = set_feature_metadata(gobject = comb_gobject,
                                          metadata = fm_list[[1]],
                                          set_defaults = FALSE)
      ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

    }

  }





  ## spatial info
  if(isTRUE(verbose)) wrap_msg('5. spatial polygon information \n')

    available_spat_info = unique(unlist(all_spatinfo_list))

  if(isTRUE(verbose)) {
    wrap_msg ('available_spat_info: \n')
    print(available_spat_info)
  }

  for(spat_info in available_spat_info) {

    savelist_vector = list()
    savelist_centroids = list()
    for(gobj_i in seq_along(updated_object_list)) {

      #print(updated_object_list[[gobj_i]]@spatial_info[[spat_info]])

      spat_information_vector = updated_object_list[[gobj_i]]@spatial_info[[spat_info]]@spatVector
      savelist_vector[[gobj_i]] = spat_information_vector

      spat_information_centroids = updated_object_list[[gobj_i]]@spatial_info[[spat_info]]@spatVectorCentroids
      savelist_centroids[[gobj_i]] = spat_information_centroids

      # TODO: add overlaps

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
  if(verbose) wrap_msg('6. spatial feature/points information \n')



  for(feat in comb_gobject@expression_feat) {

    savelist_vector = list()

    for(gobj_i in seq_along(updated_object_list)) {

      if(is.null(updated_object_list[[gobj_i]]@feat_info)) {
        spat_point_vector = NULL
      } else {
        spat_point_vector = updated_object_list[[gobj_i]]@feat_info[[feat]]@spatVector
      }

      savelist_vector[[gobj_i]] = spat_point_vector

      # TODO: add network

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


  ## If no feature_metadata exists, then generate now
  if(is.null(list_cell_metadata(comb_gobject))) {
    comb_gobject = init_feat_metadata()
  }





  ## images
  if(verbose) wrap_msg('7. images \n')

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



