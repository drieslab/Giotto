

# Functions to read data into Giotto, typically as a nested list



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
                          transpose = FALSE,
                          feat_type = 'rna',
                          expression_matrix_class = c('dgCMatrix', 'HDF5Matrix', 'rhdf5'),
                          h5_file = NULL) {

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
  
  if(expression_matrix_class[1] == 'HDF5Matrix') {
    require(HDF5Array)
    spM = methods::as(spM, 'HDF5Matrix')
  }
  
  if(expression_matrix_class[1] == 'rhdf5') {
    if(is.null(h5_file)) {
      stop(wrap_txt('h5_file can not be NULL, please provide a file name',
                    errWidth = TRUE))
    }
    
    rhdf5::h5createGroup(h5_file, paste0('expression/',feat_type))
    
    spM = as.matrix(DT[,-1])
    colnames(spM) = colnames(DT[,-1])
    rownames(spM) = DT[[1]]
    
    HDF5Array::writeHDF5Array(spM, 
                              h5_file,
                              name = paste0('expression/',feat_type,'/raw'),
                              with.dimnames=TRUE)
    
    spM = paste0('expression/',feat_type,'/raw')
  }
  
  return(spM)
}











#' @title Read expression data
#' @name readExprData
#' @description Read a nested list of expression data inputs in order to
#' generate a list of giotto-native exprObj
#' @param expr_list (nested) list of expression input data
#' @param sparse (boolean, default = TRUE) read matrix data in a sparse manner
#' @param cores number of cores to use
#' @param default_feat_type default feature type to use
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
#' @export
readExprData = function(data_list,
                        sparse = TRUE,
                        cores = determine_cores(),
                        default_feat_type = NULL,
                        verbose = TRUE,
                        provenance = NULL,
                        expression_matrix_class = c('dgCMatrix', 'HDF5Matrix', 'rhdf5'),
                        h5_file = NULL) {

  read_expression_data(
    expr_list = data_list,
    sparse = sparse,
    cores = cores,
    default_feat_type = default_feat_type,
    verbose = verbose,
    provenance = provenance,
    expression_matrix_class = expression_matrix_class,
    h5_file = h5_file
  )

}


#' @keywords internal
#' @noRd
read_expression_data = function(expr_list = NULL,
                                sparse = TRUE,
                                cores = determine_cores(),
                                default_spat_unit = NULL,
                                default_feat_type = NULL,
                                verbose = TRUE,
                                provenance = NULL,
                                expression_matrix_class = c('dgCMatrix', 'HDF5Matrix', 'rhdf5'),
                                h5_file = NULL) {

  # import box characters
  ch = box_chars()

  # Check if input exists
  if(is.null(expr_list)) return(NULL)

  # if not list make list with default name of raw
  if(!inherits(expr_list, 'list')) {
    expr_list = list('raw' = expr_list) # single matrix or path (expected)
  }


  # Set default feature type if missing
  if(is.null(default_spat_unit)) default_spat_unit = 'cell'
  if(is.null(default_feat_type)) default_feat_type = 'rna'


  # 1. get depth of list

  list_depth = depth(expr_list)



  # too much information
  if(list_depth > 3L) {
    stop('Depth of expression list is more than 3, only 3 levels are possible:
       0)', ch$s, '.
       1)', ch$s, ch$b, 'spatial unit (e.g. cell)
       2)', ch$s, ch$s, ch$b, 'feature (e.g. RNA)
       3)', ch$s, ch$s, ch$s, ch$b, 'data type (e.g. raw)\n')
  }



  # list reading
  obj_list = list()
  spat_unit_list = c()
  feat_type_list = c()
  name_list = c()


  # read nesting
  if(list_depth == 1L) {
    if(isTRUE(verbose)) message('list depth of 1')

    obj_names = names(expr_list)
    if(is.null(obj_names) & isTRUE(verbose))
      wrap_msg('No list names for objects. Setting defaults.')

    for(obj_i in seq_along(expr_list)) {

      ex = expr_list[[obj_i]]
      name = if(is_empty_char(obj_names[[obj_i]])) paste0('data_', obj_i) else obj_names[[obj_i]]

      obj_list[[length(obj_list) + 1L]] = ex
      name_list = c(name_list, name)
    }
    spat_unit_list = rep(default_spat_unit, length(obj_list)) # assume
    feat_type_list = rep(default_feat_type, length(obj_list)) # assume

  } else if(list_depth == 2L) {
    if(isTRUE(verbose)) message('list depth of 2')

    feat_type_names = names(expr_list)
    if(is.null(feat_type_names) & isTRUE(verbose))
      wrap_msg('No list names for feat_type. Setting defaults.')

    for(feat_i in seq_along(expr_list)) {

      obj_names = names(expr_list[[feat_i]])
      if(is.null(obj_names) & isTRUE(verbose))
        wrap_msg('No list names for objects. Setting defaults.')

      for(obj_i in seq_along(expr_list[[feat_i]])) {

        ex = expr_list[[feat_i]][[obj_i]]
        name = if(is_empty_char(obj_names[[obj_i]])) paste0('data_', obj_i) else obj_names[[obj_i]]
        feat_type = if(is_empty_char(feat_type_names[[feat_i]])) paste0('feat_', feat_i) else feat_type_names[[feat_i]]

        obj_list[[length(obj_list) + 1L]] = ex
        name_list = c(name_list, name)
        feat_type_list = c(feat_type_list, feat_type)
      }
    }
    spat_unit_list = rep(default_spat_unit, length(obj_list)) # assume

  } else if(list_depth == 3L) {
    if(isTRUE(verbose)) message('list depth of 3')

    spat_unit_names = names(expr_list)
    if(is.null(spat_unit_names) & isTRUE(verbose))
      wrap_msg('No list names for spat_unit. Setting defaults.')

    for(unit_i in seq_along(expr_list)) {

      feat_type_names = names(expr_list[[unit_i]])
      if(is.null(feat_type_names) & isTRUE(verbose))
        wrap_msg('No list names for feat_type. Setting defaults.')

      for(feat_i in seq_along(expr_list[[unit_i]])) {

        obj_names = names(expr_list[[unit_i]][[feat_i]])
        if(is.null(obj_names) & isTRUE(verbose))
          wrap_msg('No list names for objects. Setting defaults.')

        for(obj_i in seq_along(expr_list[[unit_i]][[feat_i]])) {

          ex = expr_list[[unit_i]][[feat_i]][[obj_i]]
          name = if(is_empty_char(obj_names[[obj_i]])) paste0('data_', obj_i) else obj_names[[obj_i]]
          feat_type = if(is_empty_char(feat_type_names[[feat_i]])) paste0('feat_', feat_i) else feat_type_names[[feat_i]]
          spat_unit = if(is_empty_char(spat_unit_names[[unit_i]])) paste0('unit_', unit_i) else spat_unit_names[[unit_i]]

          obj_list[[length(obj_list) + 1L]] = ex
          name_list = c(name_list, name)
          feat_type_list = c(feat_type_list, feat_type)
          spat_unit_list = c(spat_unit_list, spat_unit)
        }
      }
    }
  } else {
    stop(wrap_txt('Unexpected list depth', errWidth = TRUE))
  }


  if(length(obj_list) > 0L) {
    
    if(expression_matrix_class[1] == 'rhdf5') {
      if(file.exists(h5_file)) {
        wrap_txt("h5_file already exists, contents will be overwritten")
        file.remove(h5_file)}
      
      rhdf5::h5createFile(h5_file)
      rhdf5::h5createGroup(h5_file,"expression")
    }

    return_list = lapply(seq_along(obj_list), function(obj_i) {

      if(inherits(obj_list[[obj_i]], 'exprObj')) {
        warning(wrap_txt('List item [', obj_i,']: Not possible to read exprObj.
                         Returning without modifications', sep = ''))
        return(obj_list[[obj_i]])
      } else {

        # Get data from collection lists
        name = name_list[[obj_i]]
        exprMat = obj_list[[obj_i]]
        spat_unit = spat_unit_list[[obj_i]]
        feat_type = feat_type_list[[obj_i]]

        if(isTRUE(verbose)) {
          wrap_msg('\nList item [', obj_i, ']:',
                   '\nspat_unit: ', spat_unit,
                   '\nfeat_type: ', feat_type,
                   '\nname: ', name,
                   sep = ''
          )
        }

        return(
          createExprObj(
            name = name,
            expression_data = exprMat,
            spat_unit = spat_unit,
            feat_type = feat_type,
            provenance = if(is_empty_char(provenance)) spat_unit else provenance, # assumed
            misc = NULL,
            expression_matrix_class = expression_matrix_class,
            h5_file = h5_file
          )
        )
      }
    })
    return(return_list)

  } else {
    warning('No objects found in expression data input list')
  }


}









#### Giotto cell metadata ####




#' @title Read cell metadata
#' @name read_cell_metadata
#' @description read cell metadata from list
#' @param data_list nested list of cell metadata information
#' @param provenance provenance information (optional)
#' @param verbose be verbose
#' @export
readCellMetadata = function(data_list,
                            default_spat_unit = NULL,
                            default_feat_type = NULL,
                            provenance = NULL,
                            verbose = TRUE) {
  read_cell_metadata(metadata = data_list,
                     default_spat_unit = default_spat_unit,
                     default_feat_type = default_feat_type,
                     provenance = provenance,
                     verbose = verbose)
}




#' @title Read cell metadata
#' @name read_cell_metadata
#' @description read cell metadata from list
#' @param metadata nested list of cell metadata information
#' @param provenance provenance information (optional)
#' @param verbose be verbose
#' @keywords internal
read_cell_metadata = function(metadata,
                              default_spat_unit = NULL,
                              default_feat_type = NULL,
                              provenance = NULL,
                              verbose = TRUE) {

  # data.table vars
  cell_ID = NULL

  # check if input exists
  if(is.null(metadata)) return(NULL)

  # set defaults if missing
  if(is.null(default_spat_unit)) default_spat_unit = 'cell'
  if(is.null(default_feat_type)) default_feat_type = 'rna'

  # if not list make list
  if(!inherits(metadata, 'list')) {
    meta_list = list()
    meta_list[[default_spat_unit]][[default_feat_type]] = metadata
    metadata = meta_list
  }


  # 1. get detph of list
  list_depth = depth(metadata)


  # list reading
  obj_list = list()
  spat_unit_list = c()
  feat_type_list = c()


  # read nesting
  if(list_depth == 1L) {

    if(isTRUE(verbose)) message('list depth of 1')

    feat_type_names = names(metadata)
    if(is.null(feat_type_names) & isTRUE(verbose))
      wrap_msg('No list names for feat_type. Setting defaults.')

    for(feat_i in seq_along(metadata)) {

      meta = metadata[[feat_i]]
      feat_type = if(is_empty_char(feat_type_names[[feat_i]])) paste0('feat_', feat_i) else feat_type_names[[feat_i]]

      obj_list[[length(obj_list) + 1L]] = meta
      feat_type_list = c(feat_type_list, feat_type)
    }
    spat_unit_list = rep(default_spat_unit, length(obj_list)) # assume

  } else if(list_depth == 2L) {
    if(isTRUE(verbose)) message('list depth of 2')

    spat_unit_names = names(metadata)
    if(is.null(spat_unit_names) & isTRUE(verbose))
      wrap_msg('No list names for spat_unit. Setting defaults.')

    for(unit_i in seq_along(metadata)) {

      feat_type_names = names(metadata[[unit_i]])
      if(is.null(feat_type_names) & isTRUE(verbose))
        wrap_msg('No list names for feat_type. Setting defaults.')

      for(feat_i in seq_along(metadata[[unit_i]])) {

        meta = metadata[[unit_i]][[feat_i]]
        spat_unit = if(is_empty_char(spat_unit_names[[unit_i]])) paste0('unit_', unit_i) else spat_unit_names[[unit_i]]
        feat_type = if(is_empty_char(feat_type_names[[feat_i]])) paste0('feat_', feat_i) else feat_type_names[[feat_i]]

        obj_list[[length(obj_list) + 1L]] = meta
        spat_unit_list = c(spat_unit_list, spat_unit)
        feat_type_list = c(feat_type_list, feat_type)
      }
    }
  } else {
    stop(wrap_txt('Unexpected list depth', errWidth = TRUE))
  }

  if(length(obj_list) > 0L) {

    return_list = lapply(seq_along(obj_list), function(obj_i) {

      if(inherits(obj_list[[obj_i]], 'cellMetaObj')) {
        warning(wrap_txt('List item [', obj_i, ']: Not possible to read cellMetaObj.
                         Returning without modifications', sep = ''))
        return(obj_list[[obj_i]])
      } else {

        # Get data from collection lists
        meta_data = obj_list[[obj_i]]
        spat_unit = spat_unit_list[[obj_i]]
        feat_type = feat_type_list[[obj_i]]

        if(isTRUE(verbose)) {
          wrap_msg('\nList item [', obj_i, ']:',
                   '\nspat_unit: ', spat_unit,
                   '\nfeat_type: ', feat_type,
                   sep = ''
          )
        }

        return(
          createCellMetaObj(
            metadata = meta_data,
            spat_unit = spat_unit,
            feat_type = feat_type,
            provenance = if(is_empty_char(provenance)) spat_unit else provenance, # assumed
            verbose = TRUE,
            col_desc = NA_character_ # unknown
          )
        )
      }
    })
    return(return_list)

  } else {
    warning('No objects found in cell metadata input list')
  }

}








#### Giotto feat metadata ####




#' @title Read feature metadata
#' @name read_feature_metadata
#' @description read feature metadata from listt
#' @param data_list nested list of feature metadata information
#' @param provenance provenance information (optional)
#' @param verbose be verbose
#' @export
readFeatMetadata = function(data_list,
                            default_spat_unit = NULL,
                            default_feat_type = NULL,
                            provenance = NULL,
                            verbose = TRUE) {
  read_feature_metadata(metadata = data_list,
                        default_spat_unit = NULL,
                        default_feat_type = NULL,
                        provenance = provenance,
                        verbose = verbose)
}




#' @keywords internal
#' @noRd
read_feature_metadata = function(metadata,
                                 default_spat_unit = NULL,
                                 default_feat_type = NULL,
                                 provenance = NULL,
                                 verbose = TRUE) {

  # data.table vars
  cell_ID = NULL

  # check if input exists
  if(is.null(metadata)) return(NULL)

  # set defaults if missing
  if(is.null(default_spat_unit)) default_spat_unit = 'cell'
  if(is.null(default_feat_type)) default_feat_type = 'rna'

  # if not list make list
  if(!inherits(metadata, 'list')) {
    meta_list = list()
    meta_list[[default_spat_unit]][[default_feat_type]] = metadata
    metadata = meta_list
  }


  # 1. get detph of list
  list_depth = depth(metadata)


  # list reading
  obj_list = list()
  spat_unit_list = c()
  feat_type_list = c()


  # read nesting
  if(list_depth == 1L) {

    if(isTRUE(verbose)) message('list depth of 1')

    feat_type_names = names(metadata)
    if(is.null(feat_type_names) & isTRUE(verbose))
      wrap_msg('No list names for feat_type. Setting defaults.')

    for(feat_i in seq_along(metadata)) {

      meta = metadata[[feat_i]]
      feat_type = if(is_empty_char(feat_type_names[[feat_i]])) paste0('feat_', feat_i) else feat_type_names[[feat_i]]

      obj_list[[length(obj_list) + 1L]] = meta
      feat_type_list = c(feat_type_list, feat_type)
    }
    spat_unit_list = rep(default_spat_unit, length(obj_list)) # assume

  } else if(list_depth == 2L) {
    if(isTRUE(verbose)) message('list depth of 2')

    spat_unit_names = names(metadata)
    if(is.null(spat_unit_names) & isTRUE(verbose))
      wrap_msg('No list names for spat_unit. Setting defaults.')

    for(unit_i in seq_along(metadata)) {

      feat_type_names = names(metadata[[unit_i]])
      if(is.null(feat_type_names) & isTRUE(verbose))
        wrap_msg('No list names for feat_type. Setting defaults.')

      for(feat_i in seq_along(metadata[[unit_i]])) {

        meta = metadata[[unit_i]][[feat_i]]
        spat_unit = if(is_empty_char(spat_unit_names[[unit_i]])) paste0('unit_', unit_i) else spat_unit_names[[unit_i]]
        feat_type = if(is_empty_char(feat_type_names[[feat_i]])) paste0('feat_', feat_i) else feat_type_names[[feat_i]]

        obj_list[[length(obj_list) + 1L]] = meta
        spat_unit_list = c(spat_unit_list, spat_unit)
        feat_type_list = c(feat_type_list, feat_type)
      }
    }
  } else {
    stop(wrap_txt('Unexpected list depth', errWidth = TRUE))
  }

  if(length(obj_list) > 0L) {

    return_list = lapply(seq_along(obj_list), function(obj_i) {

      if(inherits(obj_list[[obj_i]], 'featMetaObj')) {
        warning(wrap_txt('List item [', obj_i, ']: Not possible to read featMetaObj.
                         Returning without modifications', sep = ''))
        return(obj_list[[obj_i]])
      } else {

        # Get data from collection lists
        meta_data = obj_list[[obj_i]]
        spat_unit = spat_unit_list[[obj_i]]
        feat_type = feat_type_list[[obj_i]]

        if(isTRUE(verbose)) {
          wrap_msg('\nList item [', obj_i, ']:',
                   '\nspat_unit: ', spat_unit,
                   '\nfeat_type: ', feat_type,
                   sep = ''
          )
        }

        return(
          createFeatMetaObj(
            metadata = meta_data,
            spat_unit = spat_unit,
            feat_type = feat_type,
            provenance = if(is_empty_char(provenance)) spat_unit else provenance, # assumed
            verbose = TRUE,
            col_desc = NA_character_ # unknown
          )
        )
      }
    })
    return(return_list)

  } else {
    warning('No objects found in feature metadata input list')
  }

}









#### Giotto locations ####









#' @title Read spatial location data
#' @name readSpatLocsData
#' @description read spatial locations/coordinates from nested list and generate
#' list of Giotto spatLocsObj
#' @inheritParams read_data_params
#' @param data_list (nested) list of spatial locations input data
#' @param cores how many cores to use
#' @return list of spatLocsObj
#' @export
readSpatLocsData = function(data_list,
                            default_spat_unit = NULL,
                            provenance = NULL,
                            cores = determine_cores(),
                            verbose = TRUE) {

  spatLocsObj_list = read_spatial_location_data(
    spat_loc_list = data_list,
    default_spat_unit = default_spat_unit,
    cores = cores,
    provenance = provenance,
    verbose = verbose
  )

  return(spatLocsObj_list)
}



#' @noRd
read_spatial_location_data = function(spat_loc_list,
                                      default_spat_unit = NULL,
                                      provenance = NULL,
                                      cores = determine_cores(),
                                      verbose = TRUE) {

  # data.table vars
  cell_ID = NULL

  if(is.null(spat_loc_list)) return(NULL)

  #  if not list make list with default name of raw
  if(!inherits(spat_loc_list, 'list')) {
    spat_loc_list = list(raw = spat_loc_list) # single matrix or path (expected)
  }


  if(is.null(default_spat_unit)) default_spat_unit = 'cell'


  # 1. get depth of list

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
    if(is.null(obj_names) & isTRUE(verbose))
      wrap_msg('No list names for objects. Setting defaults.')

    for(obj_i in seq_along(spat_loc_list)) {

      spatlocs = spat_loc_list[[obj_i]]
      name = if(is_empty_char(obj_names[[obj_i]])) paste0('coord_', obj_i) else obj_names[[obj_i]]

      obj_list[[length(obj_list) + 1L]] = spatlocs
      name_list = c(name_list, name)
    }
    spat_unit_list = rep(default_spat_unit, length(obj_list)) # add default region = 'cell'

    # for list with 2 depth, expect name info and spat_unit info
  } else if(list_depth == 2L) {
    if(isTRUE(verbose)) message('list depth of 2')

    spat_unit_names = names(spat_loc_list)
    if(is.null(spat_unit_names) & isTRUE(verbose))
      wrap_msg('No list names for spat_unit. Setting defaults.')

    for(unit_i in seq_along(spat_loc_list)) {

      obj_names = names(spat_loc_list[[unit_i]])
      if(is.null(obj_names) & isTRUE(verbose))
        wrap_msg('No list names for object. Setting defaults.')

      for(obj_i in seq_along(spat_loc_list[[unit_i]])) {

        spatlocs = spat_loc_list[[unit_i]][[obj_i]]
        spat_unit = if(is_empty_char(spat_unit_names[[unit_i]])) paste0('unit_', unit_i) else spat_unit_names[[unit_i]]
        name = if(is_empty_char(obj_names[[obj_i]])) paste0('coord_', obj_i) else obj_names[[obj_i]]

        obj_list[[length(obj_list) + 1L]] = spatlocs
        spat_unit_list = c(spat_unit_list, spat_unit)
        name_list = c(name_list, name)

      }
    }

  } else {
    stop('unexpected list_depth error')
  }


  # create spatLocObj return lst
  if(length(obj_list) > 0) {

    return_list = lapply(seq_along(obj_list), function(obj_i) {

      if(inherits(obj_list[[obj_i]], 'spatLocsObj')) {
        warning(wrap_txt('\nList item [', obj_i,']: Not possible to read spatLocsObj.
                         Returning without modifications', sep = ''))
        return(obj_list[[obj_i]])
      } else {

        # Get data from collection lists
        name = name_list[[obj_i]]
        coordinates = obj_list[[obj_i]]
        spat_unit = spat_unit_list[[obj_i]]

        if(isTRUE(verbose)) {
          wrap_msg('\nList item [', obj_i, ']:',
                   '\nspat_unit: ', spat_unit,
                   '\nname: ', name,
                   sep = ''
          )
        }

        return(
          createSpatLocsObj(
            name = name,
            coordinates = coordinates,
            spat_unit = spat_unit,
            provenance = if(is_empty_char(provenance)) spat_unit else provenance, # assumed
            misc = NULL,
            verbose = verbose
          )
        )
      }
    })
    return(return_list)
  } else {
    warning(wrap_txt('No objects found in spatial locations input list'))
  }

}






#### Giotto spatial network ####











#' @title Read spatial networks
#' @name readSpatNetData
#' @inheritParams read_data_params
#' @description read spatial networks data from list
#' @param data_list (nested) list of spatial network input data
#' @export
readSpatNetData = function(data_list,
                           default_spat_unit = NULL,
                           provenance = NULL,
                           verbose = TRUE) {
  read_spatial_networks(spatial_network = data_list,
                        default_spat_unit = default_spat_unit,
                        provenance = provenance,
                        verbose = verbose)
}




#' @keywords internal
#' @noRd
read_spatial_networks = function(spatial_network,
                                 default_spat_unit = NULL,
                                 provenance = NULL,
                                 verbose = TRUE) {

  if(is.null(spatial_network)) {
    wrap_msg('No spatial networks are provided')
    return(NULL)
  }

  # if not a list, make list
  if(!inherits(spatial_network, 'list')) {
    return(list(spatial_network))
  }

  # set defaults
  if(is.null(default_spat_unit)) default_spat_unit = 'cell'

  # list reading
  obj_list = list()
  spat_unit_list = c()
  name_list = c()


  # read nesting
  if(depth(spatial_network) == 1L) {
    if(isTRUE(verbose)) wrap_msg('list depth of 1')

    obj_names = names(spatial_network)
    if(is.null(obj_names) & isTRUE(verbose))
      wrap_msg('No list names for objects. Setting defaults.')

    for(obj_i in seq_along(spatial_network)) {

      network = spatial_network[[obj_i]]
      name = if(is_empty_char(obj_names[[obj_i]])) paste0('sn_', obj_i) else obj_names[[obj_i]]

      obj_list[[length(obj_list) + 1L]] = network
      name_list = c(name_list, name)
    }
    spat_unit_list = rep(default_spat_unit, length(obj_list)) # assumed

  } else if(depth(spatial_network) == 2L) {
    if(isTRUE(verbose)) wrap_msg('list depth of 2')

    spat_unit_names = names(spatial_network)
    if(is.null(spat_unit_names) & isTRUE(verbose))
      wrap_msg('No list names for spat_unit. Setting defaults.')

    for(unit_i in seq_along(spatial_network)) {

      obj_names = names(spatial_network[[unit_i]])
      if(is.null(obj_names) & isTRUE(verbose))
        wrap_msg('No list names for objects. Setting defaults.')
      for(obj_i in seq_along(spatial_network[[unit_i]])) {

        network = spatial_network[[unit_i]][[obj_i]]
        spat_unit = if(is_empty_char(spat_unit_names[[unit_i]])) paste0('unit_', unit_i) else spat_unit_names[[unit_i]]
        name = if(is_empty_char(obj_names[[obj_i]])) paste0('sn_', obj_i) else obj_names[[obj_i]]

        obj_list[[length(obj_list) + 1L]] = network
        spat_unit_list = c(spat_unit_list, spat_unit)
        name_list = c(name_list, name)

      }
    }

  } else {
    stop(wrap_txt('Invalid nesting depth'))
  }


  # create spatialNetworkObj return list
  if(length(obj_list) > 0) {

    return_list = lapply(seq_along(obj_list), function(obj_i) {

      if(inherits(obj_list[[obj_i]], 'spatialNetworkObj')) {
        warning(wrap_txt('\nList item [', obj_i,']: Not possible to read spatialNetworkObj.
                         Returning without modifications', sep = ''))
        return(obj_list[[obj_i]])
      } else {

        # Get data from collection lists
        name = name_list[[obj_i]]
        method = name_list[[obj_i]] # assumed
        spat_unit = spat_unit_list[[obj_i]]
        networkDT = obj_list[[obj_i]]

        if(isTRUE(verbose)) {
          wrap_msg('\nList item [', obj_i, ']:',
                   '\nspat_unit: ', spat_unit,
                   # no method nesting to report
                   '\nname: ', name,
                   sep = ''
          )
        }

        return(
          createSpatNetObj(
            name = name,
            method = method,
            spat_unit = spat_unit,
            provenance = if(is_empty_char(provenance)) spat_unit else provenance, # assumed
            network = networkDT,
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

}















#### Giotto spatial enrichment ####







#' @title Read spatial enrichment
#' @name readSpatEnrichData
#' @description read spatial enrichment results from list
#' @inheritParams read_data_params
#' @param data_list (nested) list of spatial enrichment input data
#' @export
readSpatEnrichData = function(data_list,
                              default_spat_unit = NULL,
                              default_feat_type = NULL,
                              provenance = NULL,
                              verbose = TRUE) {
  read_spatial_enrichment(spatial_enrichment = data_list,
                          default_spat_unit = default_spat_unit,
                          default_feat_type = default_feat_type,
                          provenance = provenance,
                          verbose = verbose)
}




#' @keywords internal
#' @noRd
read_spatial_enrichment = function(spatial_enrichment,
                                   default_spat_unit = NULL,
                                   default_feat_type = NULL,
                                   provenance = NULL,
                                   verbose = TRUE) {

  if(is.null(spatial_enrichment)) {
    message('No spatial enrichment results are provided')
    return(NULL)
  }

  # if not list make list
  if(!inherits(spatial_enrichment, 'list')) {
    return(list(spatial_enrichment))
  }

  # set default vals
  if(is.null(default_spat_unit)) default_spat_unit = 'cell'
  if(is.null(default_feat_type)) default_feat_type = 'rna'

  # list reading
  obj_list = list()
  spat_unit_list = c()
  feat_type_list = c()
  name_list = c()
  method_list = c()



  if(depth(spatial_enrichment) == 1L) {        # ------------------------ 1 #
    if(isTRUE(verbose)) wrap_msg('list depth of 1')

    obj_names = names(spatial_enrichment)
    if(is.null(obj_names) & isTRUE(verbose))
      wrap_msg('No list names for object. Setting defaults.')

    for(obj_i in seq_along(spatial_enrichment)) {

      enr = spatial_enrichment[[obj_i]]
      name = if(is_empty_char(obj_names[[obj_i]])) paste0('enr_', obj_i) else obj_names[[obj_i]]
      method = name # assume

      obj_list[[length(obj_list) + 1L]] = enr
      name_list = c(name_list, name)
      method_list = c(method_list, method)

    }
    feat_type_list = rep(default_feat_type, length(obj_list)) # assume
    spat_unit_list = rep(default_spat_unit, length(obj_list)) # assume

  } else if(depth(spatial_enrichment) == 2L) { # ------------------------ 2 #
    if(isTRUE(verbose)) wrap_msg('list depth of 2')

    feat_type_names = names(spatial_enrichment)
    if(is.null(feat_type_names) & isTRUE(verbose))
      wrap_msg('No list names for feat_type. Setting defaults.')
    for(feat_i in seq_along(spatial_enrichment)) {

      obj_names = names(spatial_enrichment[[feat_i]])
      if(is.null(obj_names) & isTRUE(verbose))
        wrap_msg('No list names for object. Setting defaults.')

      for(obj_i in seq_along(spatial_enrichment[[feat_i]])) {

        enr = spatial_enrichment[[feat_i]][[obj_i]]
        feat_type = if(is_empty_char(feat_type_names[[feat_i]])) paste0('feat_', feat_i) else feat_type_names[[feat_i]]
        name = if(is_empty_char(obj_names[[obj_i]])) paste0('enr_', obj_i) else obj_names[[obj_i]]
        method = name # assume

        obj_list[[length(obj_list) + 1L]] = enr
        feat_type_list = c(feat_type_list, feat_type)
        name_list = c(name_list, name)
        method_list = c(method_list, method)

      }
    }
    spat_unit_list = rep(default_spat_unit, length(obj_list)) # assume

  } else if(depth(spatial_enrichment) == 3L) { # ------------------------ 3 #
    if(isTRUE(verbose)) wrap_msg('list depth of 3')

    spat_unit_names = names(spatial_enrichment)
    if(is.null(spat_unit_names) & isTRUE(verbose))
      wrap_msg('No list names for spat_unit. Setting defaults.')

    for(unit_i in seq_along(spatial_enrichment)) {

      feat_type_names = names(spatial_enrichment[[unit_i]])
      if(is.null(feat_type_names) & isTRUE(verbose))
        wrap_msg('No list names for feat_type. Setting defaults.')

      for(feat_i in seq_along(spatial_enrichment[[unit_i]])) {

        obj_names = names(spatial_enrichment[[unit_i]][[feat_i]])
        if(is.null(obj_names) & isTRUE(verbose))
          wrap_msg('No list names for object. Setting defaults.')

        for(obj_i in seq_along(spatial_enrichment[[unit_i]][[feat_i]])) {

          enr = spatial_enrichment[[unit_i]][[feat_i]][[obj_i]]
          spat_unit = if(is_empty_char(spat_unit_names[[unit_i]])) paste0('unit_', unit_i) else spat_unit_names[[unit_i]]
          feat_type = if(is_empty_char(feat_type_names[[feat_i]])) paste0('feat_', feat_i) else feat_type_names[[feat_i]]
          name = if(is_empty_char(obj_names[[obj_i]])) paste0('enr_', obj_i) else obj_names[[obj_i]]
          method = name # assume

          obj_list[[length(obj_list) + 1L]] = enr
          spat_unit_list = c(spat_unit_list, spat_unit)
          feat_type_list = c(feat_type_list, feat_type)
          name_list = c(name_list, name)
          method_list = c(method_list, method)

        }
      }
    }

  } else if(depth(spatial_enrichment) == 4L) { # ------------------------ 4 #
    if(isTRUE(verbose)) wrap_msg('list depth of 4')

    spat_unit_names = names(spatial_enrichment)
    if(is.null(spat_unit_names) & isTRUE(verbose))
      wrap_msg('No list names for spat_unit. Setting defaults.')

    for(unit_i in seq_along(spatial_enrichment)) {

      feat_type_names = names(spatial_enrichment[[unit_i]])
      if(is.null(feat_type_names) & isTRUE(verbose))
        wrap_msg('No list names for feat_type. Setting defaults.')

      for(feat_i in seq_along(spatial_enrichment[[unit_i]])) {

        method_names = names(spatial_enrichment[[unit_i]][[feat_i]])
        if(is.null(method_names) & isTRUE(verbose))
          wrap_msg('No list names for method. Setting defaults.')
        for(method_i in seq_along(spatial_enrichment[[unit_i]][[feat_i]])) {

          obj_names = names(spatial_enrichment[[unit_i]][[feat_i]][[method_i]])
          if(is.null(obj_names) & isTRUE(verbose))
            wrap_msg('No list names for object. Setting defaults.')

          for(obj_i in seq_along(spatial_enrichment[[unit_i]][[feat_i]][[method_i]])) {

            enr = spatial_enrichment[[unit_i]][[feat_i]][[method_i]][[obj_i]]
            feat_type = if(is_empty_char(feat_type_names[[feat_i]])) paste0('feat_', feat_i) else feat_type_names[[feat_i]]
            spat_unit = if(is_empty_char(spat_unit_names[[unit_i]])) paste0('unit_', unit_i) else spat_unit_names[[unit_i]]
            name = if(is_empty_char(obj_names[[obj_i]])) paste0('enr_', obj_i) else obj_names[[obj_i]]
            method = if(is_empty_char(method_names[[method_i]])) paste0('method_', method_i) else method_names[[method_i]]

            obj_list[[length(obj_list) + 1L]] = enr
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



  # create spatEnrObj return list
  if(length(obj_list) > 0) {

    return_list = lapply(seq_along(obj_list), function(obj_i) {

      if(inherits(obj_list[[obj_i]], 'spatEnrObj')) {
        warning(wrap_txt('\nList item [', obj_i,']: Not possible to read spatEnrObj.
                         Returning without modifications', sep = ''))
        return(obj_list[[obj_i]])
      } else {

        # Get data from collection lists
        name = name_list[[obj_i]]
        method = method_list[[obj_i]]
        enrichDT = obj_list[[obj_i]]
        spat_unit = spat_unit_list[[obj_i]]
        feat_type = feat_type_list[[obj_i]]

        if(isTRUE(verbose)) {
          wrap_msg('\nList item [', obj_i, ']:',
                   '\nspat_unit: ', spat_unit,
                   '\nfeat_type: ', feat_type,
                   '\nmethod: ', method,
                   '\nname: ', name,
                   sep = ''
          )
        }

        return(
          createSpatEnrObj(
            name = name,
            method = method,
            enrichment_data = enrichDT,
            spat_unit = spat_unit,
            feat_type = feat_type,
            provenance = if(is_empty_char(provenance)) spat_unit else provenance, # assumed
            misc = NULL,
            verbose = verbose
          )
        )
      }
    })
    return(return_list)

  } else {
    warning('No objects found in spatial enrichment input list')
  }

}







#### Giotto dimension reduction ####





#' @title Read dimensional reduction data
#' @name readDimReducData
#' @description read dimension reduction results from list
#' @inheritParams read_data_params
#' @param reduction whether dim reduction was performed on 'cels' or 'feats'
#' @param data_list (nested) list of dimension reduction input data
#' @export
readDimReducData = function(data_list,
                            default_spat_unit = NULL,
                            default_feat_type = NULL,
                            reduction = c('cells', 'feats'),
                            provenance = NULL,
                            verbose = TRUE) {
  reduction = match.arg(reduction, choices = c('cells', 'feats'))

  read_dimension_reduction(dimension_reduction = data_list,
                           default_spat_unit = default_spat_unit,
                           default_feat_type = default_feat_type,
                           reduction = reduction,
                           provenance = provenance,
                           verbose = verbose)
}



#' @keywords internal
#' @noRd
read_dimension_reduction = function(dimension_reduction,
                                    default_spat_unit = NULL,
                                    default_feat_type = NULL,
                                    reduction = c('cells', 'feats'),
                                    provenance = NULL,
                                    verbose = TRUE) {

  reduction = match.arg(reduction, choices = c('cells', 'feats'))

  if(is.null(dimension_reduction)) {
    message('No dimension reduction results are provided')
    return(NULL)
  }

  # if not list, make list
  if(!inherits(dimension_reduction, 'list')) {
    return(list(dimension_reduction))
  }

  # set defaults
  if(is.null(default_spat_unit)) default_spat_unit = 'cell'
  if(is.null(default_feat_type)) default_feat_type = 'rna'


  # list reading
  obj_list = list()
  spat_unit_list = c()
  feat_type_list = c()
  method_list = c()
  name_list = c()


  # read nesting
  if(depth(dimension_reduction) == 1L) {        # ------------------------ 1 #
    if(isTRUE(verbose)) wrap_msg('list depth of 1')

    obj_names = names(dimension_reduction)
    if(is.null(obj_names) & isTRUE(verbose))
      wrap_msg('No list names for object. Setting defaults.')

    for(obj_i in seq_along(dimension_reduction)) {

      dr = dimension_reduction[[obj_i]]
      name = if(is_empty_char(obj_names[[obj_i]])) paste0('dimRed_', obj_i) else obj_names[[obj_i]]
      method = name # assume

      obj_list[[length(obj_list) + 1]] = dr
      name_list = c(name_list, name)
      method_list = c(method_list, method)

    }
    feat_type_list = rep(default_feat_type, length(obj_list)) # assume
    spat_unit_list = rep(default_spat_unit, length(obj_list)) # assume

  } else if(depth(dimension_reduction) == 2L) { # ------------------------ 2 #
    if(isTRUE(verbose)) wrap_msg('list depth of 2')

    feat_type_names = names(dimension_reduction)
    if(is.null(feat_type_names) & isTRUE(verbose))
      wrap_msg('No list names for feat_type. Setting defaults.')

    for(feat_i in seq_along(dimension_reduction)) {

      obj_names = names(dimension_reduction[[feat_i]])
      if(is.null(obj_names) & isTRUE(verbose))
        wrap_msg('No list names for object. Setting defaults.')

      for(obj_i in seq_along(dimension_reduction[[feat_i]])) {

        dr = dimension_reduction[[feat_i]][[obj_i]]
        feat_type = if(is_empty_char(feat_type_names[[feat_i]])) paste0('feat_', feat_i) else feat_type_names[[feat_i]]
        name = if(is_empty_char(obj_names[[obj_i]])) paste0('dimRed_', obj_i) else obj_names[[obj_i]]
        method = name # assume

        obj_list[[length(obj_list) + 1]] = dr
        feat_type_list = c(feat_type_list, feat_type)
        name_list = c(name_list, name)
        method_list = c(method_list, method)

      }
    }
    spat_unit_list = rep(default_spat_unit, length(obj_list)) # assume

  } else if(depth(dimension_reduction) == 3L) { # ------------------------ 3 #
    if(isTRUE(verbose)) wrap_msg('list depth of 3')

    spat_unit_names = names(dimension_reduction)
    if(is.null(spat_unit_names) & isTRUE(verbose))
      wrap_msg('No list names for spat_unit. Setting defaults.')

    for(unit_i in seq_along(dimension_reduction)) {

      feat_type_names = names(dimension_reduction[[unit_i]])
      if(is.null(feat_type_names) & isTRUE(verbose))
        wrap_msg('No list names for feat_type. Setting defaults.')

      for(feat_i in seq_along(dimension_reduction[[unit_i]])) {

        obj_names = names(dimension_reduction[[unit_i]][[feat_i]])
        if(is.null(obj_names) & isTRUE(verbose))
          wrap_msg('No list names for object. Setting defaults.')

        for(obj_i in seq_along(dimension_reduction[[unit_i]][[feat_i]])) {

          dr = dimension_reduction[[unit_i]][[feat_i]][[obj_i]]
          feat_type = if(is_empty_char(feat_type_names[[feat_i]])) paste0('feat_', feat_i) else feat_type_names[[feat_i]]
          spat_unit = if(is_empty_char(spat_unit_names[[unit_i]])) paste0('unit_', unit_i) else spat_unit_names[[unit_i]]
          name = if(is_empty_char(obj_names[[obj_i]])) paste0('dimRed_', obj_i) else obj_names[[obj_i]]
          method = name # assume

          obj_list[[length(obj_list) + 1]] = dr
          feat_type_list = c(feat_type_list, feat_type)
          spat_unit_list = c(spat_unit_list, spat_unit)
          name_list = c(name_list, name)
          method_list = c(method_list, method)

        }
      }
    }

  } else if(depth(dimension_reduction) == 4L) { # ------------------------ 4 #
    if(isTRUE(verbose)) wrap_msg('list depth of 4')

    spat_unit_names = names(dimension_reduction)
    if(is.null(spat_unit_names) & isTRUE(verbose))
      wrap_msg('No list names for spat_unit. Setting defaults.')

    for(unit_i in seq_along(dimension_reduction)) {

      feat_type_names = names(dimension_reduction[[unit_i]])
      if(is.null(feat_type_names) & isTRUE(verbose))
        wrap_msg('No list names for feat_type. Setting defaults.')

      for(feat_i in seq_along(dimension_reduction[[unit_i]])) {

        method_names = names(dimension_reduction[[unit_i]][[feat_i]])
        if(is.null(method_names) & isTRUE(verbose))
          wrap_msg('No list names for method. Setting defaults.')

        for(method_i in seq_along(dimension_reduction[[unit_i]][[feat_i]])) {

          obj_names = names(dimension_reduction[[unit_i]][[feat_i]][[method_i]])
          if(is.null(obj_names) & isTRUE(verbose))
            wrap_msg('No list names for object. Setting defaults.')

          for(obj_i in seq_along(dimension_reduction[[unit_i]][[feat_i]][[method_i]])) {

            dr = dimension_reduction[[unit_i]][[feat_i]][[method_i]][[obj_i]]
            feat_type = if(is_empty_char(feat_type_names[[feat_i]])) paste0('feat_', feat_i) else feat_type_names[[feat_i]]
            spat_unit = if(is_empty_char(spat_unit_names[[unit_i]])) paste0('unit_', unit_i) else spat_unit_names[[unit_i]]
            name = if(is_empty_char(obj_names[[obj_i]])) paste0('dimRed_', obj_i) else obj_names[[obj_i]]
            method = if(is_empty_char(method_names[[method_i]])) paste0('method_', method_i) else method_names[[method_i]]

            obj_list[[length(obj_list) + 1]] = dr
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


  if(length(obj_list) > 0L) {

    return_list = lapply(seq_along(obj_list), function(obj_i) {

      if(inherits(obj_list[[obj_i]], 'dimObj')) {
        warning(wrap_txt('\nList item [', obj_i,']: Not possible to read dimObj.
                         Returning without modifications', sep = ''))
        return(obj_list[[obj_i]])
      } else {

        name = name_list[[obj_i]]
        coordinates = obj_list[[obj_i]]
        method = method_list[[obj_i]]
        spat_unit = spat_unit_list[[obj_i]]
        feat_type = feat_type_list[[obj_i]]

        if(isTRUE(verbose)) {
          wrap_msg('\nList item [', obj_i, ']:',
                   '\nspat_unit: ', spat_unit,
                   '\nfeat_type: ', feat_type,
                   '\nmethod: ', method,
                   '\nname: ', name,
                   sep = ''
          )
        }

        return(
          createDimObj(
            name = name,
            reduction = reduction,
            coordinates = coordinates,
            method = method,
            spat_unit = spat_unit,
            feat_type = feat_type,
            provenance = if(is_empty_char(provenance)) spat_unit else provenance, # assumed
            misc = NULL
          )
        )
      }
    })
    return(return_list)

  } else {
    warning(wrap_txt('No objects found in dimension reduction input list'))
  }

}












#### Giotto nearest network ####





#' @title Read nearest neighbor network data
#' @name readNearestNetData
#' @inheritParams read_data_params
#' @description read nearest network results from list
#' @export
readNearestNetData = function(data_list,
                              default_spat_unit = NULL,
                              default_feat_type = NULL,
                              provenance = NULL,
                              verbose = TRUE) {
  read_nearest_networks(nn_network = data_list,
                        default_spat_unit = default_spat_unit,
                        default_feat_type = default_feat_type,
                        provenance = provenance,
                        verbose = verbose)
}



#' @keywords internal
#' @noRd
read_nearest_networks = function(nn_network,
                                 default_spat_unit = NULL,
                                 default_feat_type = NULL,
                                 provenance = NULL,
                                 verbose = TRUE) {

  if(is.null(nn_network)) {
    message('No nearest network results are provided')
    return(NULL)
  }

  # if not list make list
  if(!inherits(nn_network, 'list')) {
    return(list(nn_network))
  }

  # set defaults
  if(is.null(default_spat_unit)) default_spat_unit = 'cell'
  if(is.null(default_feat_type)) default_feat_type = 'rna'


  # list reading
  obj_list = list()
  spat_unit_list = c()
  feat_type_list = c()
  method_list = c()
  name_list = c()

  # read nesting
  if(depth(nn_network, sig = 'igraph') == 1L) {       # ------------------------ 1 #
    if(isTRUE(verbose)) wrap_msg('list depth of 1')

    obj_names = names(nn_network)
    if(is.null(obj_names) & isTRUE(verbose))
      wrap_msg('No list names for object. Setting defaults.')

    for(obj_i in seq_along(nn_network)) {

      nn = nn_network[[obj_i]]
      name = if(is_empty_char(obj_names[[obj_i]])) paste0('nn_', obj_i) else obj_names[[obj_i]]
      method = name # assume

      obj_list[[length(obj_list) + 1L]] = nn
      name_list = c(name_list, name)
      method_list = c(method_list, method)
    }
    spat_unit_list = rep(default_spat_unit, length(obj_list)) # assume
    feat_type_list = rep(default_feat_type, length(obj_list)) # assume

  } else if(depth(nn_network, sig = 'igraph') == 2L) { # ------------------------ 2 #
    if(isTRUE(verbose)) wrap_msg('list depth of 2')

    feat_type_names = names(nn_network)
    if(is.null(feat_type_names) & isTRUE(verbose))
      wrap_msg('No list names for feat_type. Setting defaults.')

    for(feat_i in seq_along(nn_network)) {

      obj_names = names(nn_network[[feat_i]])
      if(is.null(obj_names) & isTRUE(verbose))
        wrap_msg('No list names for object. Setting defaults.')

      for(obj_i in seq_along(nn_network[[feat_i]])) {

        nn = nn_network[[feat_i]][[obj_i]]
        feat_type = if(is_empty_char(feat_type_names[[feat_i]])) paste0('feat_', feat_i) else feat_type_names[[feat_i]]
        name = if(is_empty_char(obj_names[[obj_i]])) paste0('nn_', obj_i) else obj_names[[obj_i]]
        method = name # assume

        obj_list[[length(obj_list) + 1L]] = nn
        name_list = c(name_list, name)
        feat_type_list = c(feat_type_list, feat_type)
        method_list = c(method_list, method)
      }
    }
    spat_unit_list = rep(default_spat_unit, length(obj_list)) # assume

  } else if(depth(nn_network, sig = 'igraph') == 3L) { # ------------------------ 3 #
    if(isTRUE(verbose)) wrap_msg('list depth of 3')

    spat_unit_names = names(nn_network)
    if(is.null(spat_unit_names) & isTRUE(verbose))
      wrap_msg('No list names for spat_unit. Setting defaults.')

    for(unit_i in seq_along(nn_network)) {

      feat_type_names = names(nn_network[[unit_i]])
      if(is.null(feat_type_names) & isTRUE(verbose))
        wrap_msg('No list names for feat_type. Setting defaults.')

      for(feat_i in seq_along(nn_network[[unit_i]])) {

        obj_names = names(nn_network[[unit_i]][[feat_i]])
        if(is.null(obj_names) & isTRUE(verbose))
          wrap_msg('No list names for object. Setting defaults.')

        for(obj_i in seq_along(nn_network[[unit_i]][[feat_i]])) {

          nn = nn_network[[unit_i]][[feat_i]][[obj_i]]
          spat_unit = if(is_empty_char(spat_unit_names[[unit_i]])) paste0('unit_', unit_i) else spat_unit_names[[unit_i]]
          feat_type = if(is_empty_char(feat_type_names[[feat_i]])) paste0('feat_', feat_i) else feat_type_names[[feat_i]]
          name = if(is_empty_char(obj_names[[obj_i]])) paste0('nn_', obj_i) else obj_names[[obj_i]]
          method = name # assume

          obj_list[[length(obj_list) + 1L]] = nn
          name_list = c(name_list, name)
          feat_type_list = c(feat_type_list, feat_type)
          spat_unit_list = c(spat_unit_list, spat_unit)
          method_list = c(method_list, method)
        }
      }
    }

  } else if(depth(nn_network, sig = 'igraph') == 4L) { # ------------------------ 4 #
    if(isTRUE(verbose)) wrap_msg('list depth of 4')

    spat_unit_names = names(nn_network)
    if(is.null(spat_unit_names) & isTRUE(verbose))
      wrap_msg('No list names for spat_unit. Setting defaults.')

    for(unit_i in seq_along(nn_network)) {

      feat_type_names = names(nn_network[[unit_i]])
      if(is.null(feat_type_names) & isTRUE(verbose))
        wrap_msg('No list names for feat_type. Setting defaults.')

      for(feat_i in seq_along(nn_network[[unit_i]])) {

        method_names = names(nn_network[[unit_i]][[feat_i]])
        if(is.null(method_names) & isTRUE(verbose))
          wrap_msg('No list names for method. Setting defaults.')

        for(method_i in seq_along(nn_network[[unit_i]][[feat_i]])) {

          obj_names = names(nn_network[[unit_i]][[feat_i]][[method_i]])
          if(is.null(obj_names) & isTRUE(verbose))
            wrap_msg('No list names for object. Setting defaults.')

          for(obj_i in seq_along(nn_network[[unit_i]][[feat_i]][[method_i]])) {

            nn = nn_network[[unit_i]][[feat_i]][[method_i]][[obj_i]]
            spat_unit = if(is_empty_char(spat_unit_names[[unit_i]])) paste0('unit_', unit_i) else spat_unit_names[[unit_i]]
            feat_type = if(is_empty_char(feat_type_names[[feat_i]])) paste0('feat_', feat_i) else feat_type_names[[feat_i]]
            name = if(is_empty_char(obj_names[[obj_i]])) paste0('nn_', obj_i) else obj_names[[obj_i]]
            method = if(is_empty_char(method_names[[method_i]])) paste0('method_', method_i) else method_names[[method_i]]

            obj_list[[length(obj_list) + 1L]] = nn
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


  if(length(obj_list) > 0L) {

    return_list = lapply(seq_along(obj_list), function(obj_i) {
      if(inherits(obj_list[[obj_i]], 'nnNetObj')) {
        warning(wrap_txt('\nList item [', obj_i,']: Not possible to read nnNetObj.
                         Returning without modifications', sep = ''))
        return(obj_list[[obj_i]])
      } else {

        # Get data from collection lists
        name = name_list[[obj_i]]
        nn_type = method_list[[obj_i]]
        igraph = obj_list[[obj_i]]
        spat_unit = spat_unit_list[[obj_i]]
        feat_type = feat_type_list[[obj_i]]

        if(isTRUE(verbose)) {
          wrap_msg('\nList item [', obj_i, ']:',
                   '\nspat_unit: ', spat_unit,
                   '\nfeat_type: ', feat_type,
                   '\nmethod: ', nn_type,
                   '\nname: ', name,
                   sep = ''
          )
        }

        return(
          createNearestNetObj(
            name = name,
            nn_type = nn_type,
            network = igraph,
            spat_unit = spat_unit,
            feat_type = feat_type,
            provenance = if(is_empty_char(provenance)) spat_unit else provenance, # assume
            misc = NULL
          )
        )
      }
    })
    return(return_list)

  } else {
    warning(wrap_txt('No objects found in nearest neighbor network input list'))
  }

}












#### Giotto spatial info ####



#' @title Read list of polygons information
#' @name readPolygonData
#' @description Function extract list of polygons when given raw input as either
#' mask or tabular data. Calls the respective createGiottoPolygons functions. \cr
#' If a \code{giottoPolygon} object is passed then no edits will be made other
#' than updating the \code{name} slot if the list is named.
#' @param input what type of input is being used. When set to 'guess', uses
#' 'mask' if \code{polygonlist} is of type character and 'table' when
#' \code{polygonlist} is dataframe-like
#' @param default_name default name to assign if \code{polygonlist} is not a
#' list. If \code{polygonlist} is an unnamed list then \code{default_name} will
#' be used as part of the template for generating indexed default names.
#' @param polygon_mask_list_params parameters for when using 'mask' workflow
#' @param polygon_dfr_list_params parameters for when using 'table' workflow
#' @param calc_centroids whether centroids should be calculated during polygon
#' creation
#' @param verbose be verbose
#' @export
readPolygonData = function(data_list,
                           default_name = 'cell',
                           input = 'guess',
                           polygon_mask_list_params = NULL,
                           polygon_dfr_list_params = NULL,
                           calc_centroids = FALSE,
                           verbose = TRUE) {

  if(is.null(data_list)) {
    message('No polygon data/spatial info is provided')
    return(NULL)
  }



  # Setup gpolygon creation settings: mask
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
                                    fix_multipart = TRUE)
  }

  # gpoly creation settings: dataframe
  if(is.null(polygon_dfr_list_params)) {
    polygon_dfr_list_params = list()
  }

  # set centroids param
  polygon_mask_list_params$calc_centroids =
    polygon_dfr_list_params$calc_centroids =
    calc_centroids



  # pass to internal
  extract_polygon_list(
    polygonlist = data_list,
    input = input,
    default_name = default_name,
    polygon_mask_list_params = polygon_mask_list_params,
    polygon_dfr_list_params = polygon_dfr_list_params,
    verbose = verbose
  )


}






#' @title Extract list of polygons
#' @name extract_polygon_list
#' @description Function extract list of polygons when given raw input as either
#' mask or tabular data. Calls the respective createGiottoPolygons functions. \cr
#' If a \code{giottoPolygon} object is passed then no edits will be made other
#' than updating the \code{name} slot if the list is named.
#' @param input what type of input is being used. When set to 'guess', uses
#' 'mask' if \code{polygonlist} is of type character and 'table' when
#' \code{polygonlist} is dataframe-like
#' @param default_name default name to assign if \code{polygonlist} is not a
#' list. If \code{polygonlist} is an unnamed list then \code{default_name} will
#' be used as part of the template for generating indexed default names.
#' @param polygon_mask_list_params parameters for when using 'mask' workflow
#' @param polygon_dfr_list_params parameters for when using 'table' workflow
#' @param verbose be verbose
#' @keywords internal
#' @noRd
extract_polygon_list = function(polygonlist,
                                input = 'guess',
                                default_name = 'cell',
                                polygon_mask_list_params,
                                polygon_dfr_list_params,
                                verbose = TRUE) {

  named_list = FALSE

  # if polygonlist is not a named list
  # try to make list and give default names
  if(!inherits(polygonlist, 'list')) {

    try_val = try(as.list(polygonlist), silent = TRUE)
    if(inherits(try_val, 'try-error')) {
      if(isTRUE(verbose)) (wrap_msg('polygonlist is not a list'))
      polygonlist = list(polygonlist)
    } else polygonlist = try_val
    if(length(polygonlist) == 1L) {
      names(polygonlist) = default_name
    } else {
      polygonlist_l = length(polygonlist)
      names(polygonlist) = c(default_name, paste0('info', 1:(polygonlist_l-1)))
    }
  } else if(is.null(names(polygonlist))) {
    # if it is list
    # test if it has names
    if(isTRUE(verbose)) wrap_msg('polygonlist is a list without names')

    if(length(polygonlist) == 1L) {
      names(polygonlist) = default_name
    } else {
      polygonlist_l = length(polygonlist)
      names(polygonlist) = c(default_name, paste0('info', 1:(polygonlist_l-1)))

    }
  } else {

    if(isTRUE(verbose)) wrap_msg('polygonlist is a list with names')
    named_list = TRUE

  }

  # make sure cell is one of the names
  #all_names = names(polygonlist)
  #if(!any('cell' %in% all_names)) {
  #  stop(" Information about 'cell' needs to be provided")
  #}


  final_list = list()

  for(poly_i in seq_along(polygonlist)) {

    name_polyinfo = names(polygonlist)[[poly_i]]
    polyinfo = polygonlist[[poly_i]]

    if(isTRUE(verbose)) wrap_msg('  [', name_polyinfo, '] Process polygon info...')

    if((is.character(polyinfo) & input == 'guess') | input == 'mask') {

      parameters = c(list(name = name_polyinfo,
                          maskfile = polyinfo),
                     polygon_mask_list_params)
      if(isTRUE(verbose)) print(parameters)

      poly_results = do.call(what = 'createGiottoPolygonsFromMask', args = parameters)

    } else if((inherits(polyinfo, 'data.frame') & input == 'guess') | input == 'table') {

      parameters = c(list(name = name_polyinfo,
                          segmdfr = polyinfo),
                     polygon_dfr_list_params)
      if(isTRUE(verbose)) print(parameters)

      poly_results = do.call(what = 'createGiottoPolygonsFromDfr', args = parameters)

    } else if(inherits(polyinfo, 'giottoPolygon')) {

      # Override name slot ONLY if giottoPolygon provided as named list
      if(isTRUE(named_list)) slot(polyinfo, 'name') = name_polyinfo
      else name_polyinfo = polyinfo@name
      poly_results = polyinfo

    } else {

      stop(wrap_txt('Polygon can only be extracted from a mask file or from a correctly formatted data.frame'))

    }

    final_list[[name_polyinfo]] = poly_results

  }


  return(final_list)

}







#' @title Add giotto polygons to giotto object
#' @name addGiottoPolygons
#' @description Adds Giotto polygon to an existing Giotto object
#' @param gobject giotto object
#' @param gpolygons list of giotto polygon objects,
#' see \code{\link{createGiottoPolygonsFromMask}} and \code{\link{createGiottoPolygonsFromDfr}}
#' @return giotto object
#' @concept polygon
#' @export
addGiottoPolygons = function(gobject,
                             gpolygons) {

  # check input
  guard_against_notgiotto(gobject)

  if(!inherits(gpolygons, 'list')) {
    stop('gpolygons needs to be a list of one or more giottoPolygon objects')
  }


  # add each giottoPoint object to the giotto object
  for(gp_i in 1:length(gpolygons)) {

    gp = gpolygons[[gp_i]]

    # check if giottoPoint object
    if(!inherits(gp, 'giottoPolygon')) {
      stop('gpolygons needs to be a list of one or more giottoPolygon objects', '\n',
           'number ', gp_i, ' is not a giottoPolygon object \n')
    }


    gobject@spatial_info[[gp@name]] = gp

  }

  return(gobject)

}











#### Giotto spatial feature info ####




#' @title Read feature information
#' @name readFeatData
#' @description Function to read lists of feature information data and output
#' a list of generated giottoPoints objects
#' @inheritParams read_data_params
#' @export
readFeatData = function(data_list,
                        verbose = TRUE) {

  if(is.null(data_list)) {
    message('No feature info is provided')
    return(NULL)
  }


  extract_points_list(pointslist = data_list,
                      verbose = verbose)
}





#' @title Extract list of giotto points objects
#' @name extract_points_list
#' @description to extract list of giotto points
#' @param pointslist list of inputs from which to create giotto points objects
#' @param verbose be verbose
#' @keywords internal
#' @noRd
extract_points_list = function(pointslist,
                               verbose = TRUE) {

  named_list = FALSE

  # if pointslist is not a named list
  # try to make list and give default names
  if(!inherits(pointslist, 'list')) {

    if(isTRUE(verbose)) wrap_msg('pointslist is not a list')

    try_val = try(as.list(pointslist), silent = TRUE)
    if(inherits(try_val, 'try-error')) {
      pointslist = list(pointslist)
    } else pointslist = try_val
    if(length(pointslist) == 1L) {
      names(pointslist) = 'rna'
    } else {
      pointslist_l = length(pointslist)
      names(pointslist) = c('rna', paste0('feat', seq(pointslist_l-1L)))
    }
  } else if(is.null(names(pointslist))) {
    # if it is list
    # test if it has names
    if(isTRUE(verbose)) wrap_msg('pointslist is a list without names')
    if(length(pointslist) == 1L) {
      names(pointslist) = 'rna'
    } else {
      pointslist_l = length(pointslist)
      names(pointslist) = c('rna', paste0('feat', seq(pointslist_l-1L)))

    }
  } else {

    if(isTRUE(verbose)) wrap_msg('pointslist is a named list')
    named_list = TRUE

  }

  # make sure rna is one of the names
  #all_names = names(pointslist)
  #if(!any('rna' %in% all_names)) {
  #  stop(" Information about 'rna' needs to be provided")
  #}


  final_list = list()

  for(point_i in seq_along(pointslist)) {

    name_pointinfo = names(pointslist)[[point_i]]
    pointinfo = pointslist[[point_i]]

    if(isTRUE(verbose)) wrap_msg('  [', name_pointinfo, '] Process point info...')

    if(inherits(pointinfo, 'giottoPoints')) {

      if(isTRUE(named_list)) slot(pointinfo, 'feat_type') = name_pointinfo
      else name_pointinfo = pointinfo@feat_type
      point_results = pointinfo

    } else if(inherits(pointinfo, 'data.frame')) {

      point_results = createGiottoPoints(x = pointinfo,
                                         feat_type = name_pointinfo)

    } else if(inherits(pointinfo, 'character')) {

      stop(wrap_txt('Giotto points can not yet be created directly from a file path'))

    } else {

      stop(wrap_txt('Giotto points can only be created from a correctly formatted data.frame-like object'))

    }

    final_list[[name_pointinfo]] = point_results

  }

  return(final_list)

}






# add Giotto points object to existing Giotto object
# cell IDs needs to match

#' @title Add giotto points object to giotto object
#' @name addGiottoPoints
#' @description Adds Giotto points to an existing Giotto object
#' @param gobject giotto object
#' @param gpoints list of giotto point objects, see \code{\link{createGiottoPoints}}
#' @return giotto object
#' @concept polygon
#' @export
addGiottoPoints = function(gobject,
                           gpoints) {

  # check input
  if(!inherits(gobject, 'giotto')) {
    stop('gobject needs to be a giotto object')
  }

  if(!inherits(gpoints, 'list')) {
    stop('gpoints needs to be a list of one or more giottoPoints objects')
  }

  # available features types
  feat_types = gobject@expression_feat


  # add each giottoPoint object to the giotto object
  for(gp_i in 1:length(gpoints)) {

    gp = gpoints[[gp_i]]

    # check if giottoPoint object
    if(!inherits(gp, 'giottoPoints')) {
      stop('gpoints needs to be a list of one or more giottoPoints objects', '\n',
           'number ', gp_i, ' is not a giottoPoints object \n')
    }

    # check if feature type is available
    if(!gp@feat_type %in% feat_types) {
      stop(gp@feat_type, ' is not a feature type in the giotto object \n')
    }

    # check if features match
    gobject_feats = gobject@feat_ID[[gp@feat_type]]
    gpoints_feats = unique(gp@spatVector[['feat_ID']][[1]])

    extra_feats = gpoints_feats[!gpoints_feats %in% gobject_feats]
    if(length(extra_feats) > 0) {
      warning(length(extra_feats), ' too many features, these features are not in the original giotto object: \n',
              paste(extra_feats, ' '), ' \n you may want to remove them')
    }

    missing_feats = gobject_feats[!gobject_feats %in% gpoints_feats]
    if(length(missing_feats) > 0) {
      warning(length(missing_feats), ' missing features, these features are not found in the giotto points object: \n',
              paste(missing_feats, ' '), ' \n you may want to add them')
    }

    gobject@feat_info[[gp@feat_type]] = gp

  }

  return(gobject)

}


#' Add sub cellular 3D coordinates to Giotto object
#'
#' @param gobject  A Giotto object.
#' @param coords A \link{data.frame} or `spatVector` with at least xyz coordinates and feature ids.
#' @param feat_type a character. The feat_type must previously exist in the Giotto object. Default = "rna".
#'
#' @return A Giotto object with a `spatVector` object in the feat_info slot
#' @export
addGiottoPoints3D <- function (gobject, coords, feat_type = "rna")
{
  # verify gobject class
  if (!inherits(gobject, "giotto")) {
    stop("gobject needs to be a giotto object")
  }
  # available features types
  feat_types = gobject@expression_feat
  if(!feat_type %in% feat_types) {
    stop(feat_type, ' is not a feature type in the giotto object \n')
  }

  if (inherits(coords, "data.frame")) {
    spatvec = terra::vect(as.matrix(coords[,1:2]), type = "points", atts = coords)
    names(spatvec)[4] = 'feat_ID'

    g_points = create_giotto_points_object(feat_type = feat_type,
                                           spatVector = spatvec)
  }
  else if (inherits(coords, "spatVector")) {
    g_points = create_giotto_points_object(feat_type = feat_type,
                                           spatVector = coords)
  }
  else {
    stop("Class ", class(coords), " is not supported")
  }

  gobject@feat_info[[g_points@feat_type]] = g_points

  return(gobject)
}









