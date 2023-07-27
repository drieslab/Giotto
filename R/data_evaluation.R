
# EVALUATION FUNCTIONS
# --------------------------------------------------------------------------- #
# Internal functions to evaluate raw inputs into formats that are directly
# compatible with Giotto's functionality
#
# Often called from the class constructor functions
# --------------------------------------------------------------------------- #












# Expression ####

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
#' @noRd
evaluate_expr_matrix = function(inputmatrix,
                                sparse = TRUE,
                                cores = determine_cores(),
                                feat_type = 'rna',
                                expression_matrix_class = c('dgCMatrix', 'HDF5Matrix', 'rhdf5'),
                                h5_file = NULL) {


  if(inherits(inputmatrix, 'character')) {
    inputmatrix = path.expand(inputmatrix)
    mymatrix = readExprMatrix(inputmatrix, cores = cores,
                              expression_matrix_class = expression_matrix_class,
                              feat_type = feat_type,
                              h5_file = h5_file)
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

    mymatrix = methods::as(as.matrix(inputmatrix), 'sparseMatrix')

  } else if(inherits(inputmatrix, 'exprObj')) {

    inputmatrix[] = evaluate_expr_matrix(inputmatrix[], sparse = sparse, cores = cores, expression_matrix_class = expression_matrix_class)
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





# Metadata ####

#' @param metadata metadata input
#' @param cores cores to use if reading in the information
#' @keywords internal
#' @noRd
evaluate_cell_metadata = function(metadata,
                                  cores = determine_cores(),
                                  verbose = TRUE) {

  # data.table vars cell_ID

  # Get data as data.table
  if(!any(class(metadata) %in% c('data.table', 'data.frame', 'matrix', 'character'))) {
    stop(wrap_txt('metadata needs to be a data.table or data.frame-like object',
                  'or a path to one of these',
                  errWidth = TRUE))
  }
  if(inherits(metadata, 'character')) {
    metadata = path.expand(metadata)
    if(!file.exists(metadata)) stop(wrap_txt('path to metadata does not exist',
                                             errWidth = TRUE))
    metadata = data.table::fread(input = metadata, nThread = cores)
  } else {
    metadata = tryCatch(
      data.table::setDT(metadata),
      error = function(e) data.table::as.data.table(metadata)
    )
  }


  # assign cell_ID col
  if('cell_ID' %in% colnames(metadata)) {
    data.table::setcolorder(metadata, neworder = c('cell_ID')) # set as first
    metadata[, cell_ID := as.character(cell_ID)] # ensure character

    # ensure unique entries
    if(any(metadata[, duplicated(cell_ID)])) {
      stop(wrap_txt('Cell metadata: duplicates found in cell_ID column.',
                    errWidth = TRUE))
    }

  } else {

    warning(wrap_txt('Cell metadata input: no col named cell_ID.
                     Setting temporary NA values'))
    # set temporary NA values
    metadata[, cell_ID := NA_character_]
    # re-order so that cell_ID is the first column
    data.table::setcolorder(metadata, neworder = "cell_ID")

  }

  return(metadata)
}









#' @keywords internal
#' @noRd
evaluate_feat_metadata = function(metadata,
                                  cores = determine_cores(),
                                  verbose = TRUE) {

  # data.table vars
  feat_ID = NULL

  # Get data as data.table
  if(!any(class(metadata) %in% c('data.table', 'data.frame', 'matrix', 'character'))) {
    stop(wrap_txt('metadata needs to be a data.table or data.frame-like object',
                  'or a path to one of these',
                  errWidth = TRUE))
  }
  if(inherits(metadata, 'character')) {
    metadata = path.expand(metadata)
    if(!file.exists(metadata)) stop(wrap_txt('path to metadata does not exist',
                                             errWidth = TRUE))
    metadata = data.table::fread(input = metadata, nThread = cores)
  } else {
    metadata = tryCatch(
      data.table::setDT(metadata),
      error = function(e) data.table::as.data.table(metadata)
    )
  }


  # assign feat_ID col
  if('feat_ID' %in% colnames(metadata)) {
    data.table::setcolorder(metadata, neworder = c('feat_ID')) # set as first
    metadata[, feat_ID := as.character(feat_ID)] # ensure character

    # ensure unique entries
    if(any(metadata[, duplicated(feat_ID)])) {
      stop(wrap_txt('Feature metadata: duplicates found in feat_ID column.',
                    errWidth = TRUE))
    }

  } else {

    warning(wrap_txt('Feature metadata input: no col named feat_ID.
                     Setting temporary NA values'))
    # set temporary NA values
    metadata[, feat_ID := NA_character_]

  }

  return(metadata)
}









# Spatial ####

#' @title evaluate_spatial_locations_OLD
#' @name evaluate_spatial_locations_OLD
#' @description Evaluate spatial location input
#' @param spatial_locs spatial locations to evaluate
#' @param cores how many cores to use
#' @param dummy_n number of rows to create dummy spaial locations
#' @param expr_matrix expression matrix to compare the cell IDs with
#' @return data.table
#' @keywords internal
#' @noRd
evaluate_spatial_locations_OLD = function(spatial_locs,
                                          cores = 1,
                                          dummy_n = 2,
                                          expr_matrix = NULL,
                                          verbose = TRUE) {

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

      if(isTRUE(verbose)) {
        wrap_msg('There are non numeric or integer columns for the spatial location input at column position(s): ', non_numeric_indices,
                 '\n The first non-numeric column will be considered as a cell ID to test for consistency with the expression matrix',
                 '\n Other non numeric columns will be removed\n\n')
      }


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
#' @noRd
evaluate_spatial_locations = function(spatial_locs,
                                      cores = determine_cores(),
                                      verbose = TRUE) {

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

    if(isTRUE(verbose)) wrap_msg(
      'There are non numeric or integer columns for the spatial location input at column position(s): ', non_numeric_indices,
      '\n The first non-numeric column will be considered as a cell ID to test for consistency with the expression matrix',
      '\n Other non numeric columns will be removed'
    )


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










#' @title Evaluate spatial network
#' @name evaluate_spatial_network
#' @description function to evaluate a spatial network
#' @keywords internal
#' @noRd
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










#' @name evaluate_spatial_enrichment
#' @description evaluate spatial enrichment input into data.table format
#' compatible with spatEnrObj
#' @keywords internal
#' @noRd
evaluate_spatial_enrichment = function(spatial_enrichment,
                                       provenance = NULL,
                                       cores = determine_cores(),
                                       verbose = TRUE) {

  # data.table vars
  cell_ID = NULL

  if(!any(class(unlist(spatial_enrichment)) %in%
          c('data.table', 'data.frame', 'matrix', 'character'))) {
    stop(wrap_txt('spatial enrichment needss to be a data.table or data.frame-like',
                  'object or a path to one of these', errWidth = TRUE))
  }
  if(inherits(spatial_enrichment, 'character')) {
    if(!file.exists(path.expand(spatial_enrichment))) {
      stop(wrap_txt('path to spatial enrichment info does not exist'))
    }

    spatial_enrichment = data.table::fread(input = spatial_enrichment,
                                           nThread = cores)
  } else {
    spatial_enrichment = tryCatch(
      data.table::setDT(spatial_enrichment),
      error = function(e) data.table::as.data.table(spatial_enrichment)
    )
  }

  # check which columns are numeric (contain enrichment info)
  column_classes = lapply(spatial_enrichment, FUN = class)
  non_numeric_classes = column_classes[!column_classes %in% c('numeric', 'integer')]


  potential_cell_IDs = NULL

  # find non-numeric cols (possible cell_ID col)
  if(length(non_numeric_classes) > 0L) {

    non_numeric_indices = which(!column_classes %in% c('numeric', 'integer'))

    if(isTRUE(verbose)) {
      wrap_msg('There are non numeric or integer columns for the spatial enrichment',
               'input at column position(s):', non_numeric_indices,
               '\nThe first non-numeric column will be considered as a cell ID to',
               'test for consistency with the expression matrix.
               Other non-numeric columns will be removed.')
    }

    potential_cell_IDs = spatial_enrichment[[names(non_numeric_classes)[[1L]]]]

    # subset to only numeric cols for testing
    spatial_enrichment = spatial_enrichment[, -non_numeric_indices, with = FALSE]

  }

  # check number of columns: too few
  if(ncol(spatial_enrichment) < 1L) {
    stop(wrap_txt('There has to be at least 2 columns (1 for cell IDs, and',
                  'at least one other for enrichment data', errWidth = TRUE))
  }

  # Assign first non-numeric as cell_ID
  if(!is.null(potential_cell_IDs)) {
    spatial_enrichment[, cell_ID := potential_cell_IDs]
  }

  return(spatial_enrichment)

}







# Embeddings ####




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








#' @title Evaluate nearest networks
#' @name evaluate_nearest_networks
#' @description Evaluate nearest networks input into igraph for input into
#' nnNetObj
#' @keywords internal
#' @details Minimal input is a data.frame-like input containing 'from', 'to',
#' and 'distance' information
#' @noRd
evaluate_nearest_networks = function(nn_network) {

  # data.table vars
  weight = distance = NULL

  if(inherits(nn_network, 'nnNetObj')) {
    nn_network[] = evaluate_nearest_networks(nn_network = nn_network[])
    return(nn_network)

  } else if(inherits(nn_network, 'igraph')) {
    v_attr = igraph::list.vertex.attributes(nn_network)
    e_attr = igraph::list.edge.attributes(nn_network)

    if(!'name' %in% v_attr) stop(wrap_txt(
      'nearest network igraph input MUST have vertex attribute "name".
      Discovered vertex attributes:', v_attr, errWidth = TRUE
    ))

    if(!'distance' %in% e_attr) stop(wrap_txt(
      'nearest network igraph input MUST have edge attribute "distance".
      Discovered edge attributes:', e_attr, errWidth = TRUE
    ))

    if(!'weight' %in% e_attr) {
      igDT = data.table::setDT(igraph::as_data_frame(nn_network))
      igDT[, weight := 1/(1 + distance)]
      data.table::setcolorder(igDT, neworder = c('from', 'to', 'weight', 'distance'))
      nn_network = igraph::graph_from_data_frame(igDT)
    }

    return(nn_network)

  } else if(inherits(nn_network, 'data.frame')) {
    if(!inherits(nn_network, 'data.table')) nn_network = data.table::setDT(nn_network)

    # if minimal input not given, throw error
    if(!all(c('from', 'to', 'distance') %in% colnames(nn_network))) {
      stop(wrap_txt('Unable to coerce data.frame type object to nnNetObj igraph
                    Needed columns: from, to, distance', errWidth = TRUE))
    }

    # generate weights
    nn_network[, weight := 1/(1 + distance)]
    data.table::setcolorder(nn_network, neworder = c('from', 'to', 'weight', 'distance'))

    nn_network = igraph::graph_from_data_frame(nn_network)
    return(nn_network)

  }

}







# Spatial/Polygon info ####

#' @describeIn createGiottoPolygonsFromDfr Examines provided data.frame type object
#' for columns that should correspond to x/y vertices and the polygon ID. Returns
#' a data.table with those key columns renamed to 'x', 'y', and 'poly_ID' if necessary.
#' @keywords internal
#' @noRd
evaluate_gpoly_dfr = function(input_dt,
                              verbose = TRUE) {

  x = y = poly_ID = NULL

  # data.frame like object needs to have 2 coordinate columns and
  # at least one other column as the poly_ID
  if(ncol(input_dt) < 3) stop('At minimum, columns for xy coordinates and poly ID are needed.\n')

  col_classes = sapply(input_dt, class)


  # 1. detect poly_ID
  ## find poly_ID as either first character col or named column
  ## if neither exist, pick the 3rd column
  if('poly_ID' %in% colnames(input_dt)) {
    poly_ID_col = which(colnames(input_dt) == 'poly_ID')
  } else {
    poly_ID_col = which(col_classes == 'character')
    if(length(poly_ID_col) < 1) poly_ID_col = 3 # case if no char found: default to 3rd
    else poly_ID_col = poly_ID_col[[1]] # case if char is found
  }


  # 2. detect x and y
  ## find first two numeric cols as x and y respectively or named column
  ## if neither exist, pick the 1st and 2nd cols respectively for x and y
  if(all(c('x','y') %in% colnames(input_dt))) {
    x_col = which(colnames(input_dt) == 'x')
    y_col = which(colnames(input_dt) == 'y')
  } else {
    x_col = which(col_classes == 'numeric')
    if(length(x_col) < 2) x_col = 1 # case if no/too few num found: default to 1st
    else x_col = x_col[[1]] # case if num found
    y_col = which(col_classes == 'numeric')
    if(length(y_col) < 2) y_col = 2 # case if no/too few num found: default to 2nd
    else y_col = y_col[[2]] # case if num found
  }


  # 3. print selections and ensure correct data type
  if(isTRUE(verbose)) message(paste0('  Selecting col "',colnames(input_dt[, poly_ID_col, with = FALSE]),'" as poly_ID column'))
  colnames(input_dt)[poly_ID_col] = 'poly_ID'
  if(!input_dt[, inherits(poly_ID, 'character')]) input_dt[, poly_ID := as.character(poly_ID)]

  if(isTRUE(verbose)) message(paste0('  Selecting cols "',colnames(input_dt[, x_col, with = FALSE]),'" and "', colnames(input_dt[, y_col, with = FALSE]),'" as x and y respectively'))
  colnames(input_dt)[x_col] = 'x'
  colnames(input_dt)[y_col] = 'y'
  if(!input_dt[, inherits(x, 'numeric')]) input_dt[, x := as.numeric(x)]
  if(!input_dt[, inherits(y, 'numeric')]) input_dt[, y := as.numeric(y)]

  input_dt
}







#' @keywords internal
#' @param input_sv SpatVector to evaluate
#' @param verbose be verbose
#' @return list of SpatVector and unique_IDs
#' @noRd
evaluate_gpoly_spatvector = function(input_sv,
                                     verbose = TRUE) {

  # determine sv type
  sv_type = terra::geomtype(input_sv)

  if(sv_type != 'polygons') {
    stop('SpatVector is of type "', sv_type, '" instead of "polygons"')
  }


  col_classes = sapply(
    sample(x = input_sv, size = 1L),
    FUN = class
  )


  # 1. detect poly_ID
  ## find poly_ID as either first character col or named column
  ## if neither exist, pick the 1st column
  sv_names = names(input_sv)
  if('poly_ID' %in% sv_names) {
    poly_ID_col = which(sv_names == 'poly_ID')
  } else {
    poly_ID_col = which(col_classes == 'character')
    if(length(poly_ID_col) < 1L) poly_ID_col = 1L # case if no char found: default to 1st
    else poly_ID_col = poly_ID_col[[1L]] # case if char is found
  }


  # 2. print selections and ensure correct data type
  if(isTRUE(verbose)) wrap_msg('  Selecting attribute "', names(input_sv[[poly_ID_col]]), '" as poly_ID', sep = '')
  sv_names[[poly_ID_col]] = 'poly_ID'
  terra::set.names(input_sv, sv_names)

  # strip crs info
  terra::set.crs(input_sv, NULL)

  unique_IDs = NULL
  if(col_classes[[poly_ID_col]] != 'character') {
    poly_ID_vals = terra::as.list(input_sv)$poly_ID
    poly_ID_vals = as.character(poly_ID_vals)

    input_sv$poly_ID = poly_ID_vals
    unique_IDs = unique(poly_ID_vals)
  }

  return_list = list(spatvector = input_sv,
                     unique_IDs = unique_IDs)

  return_list
}








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
#' @return list of SpatVector and unique polygon IDs that it contains
#' @keywords internal
#' @noRd
evaluate_spatial_info = function(spatial_info,
                                 skip_eval_dfr = FALSE,
                                 copy_dt = TRUE,
                                 cores = determine_cores(),
                                 verbose = TRUE) {

  # data.table vars
  geom = poly_ID = NULL

  ## 1. load or read spatial information data ##
  ## 1.1 read from file
  if(inherits(spatial_info, 'character')) {
    spatial_info = path.expand(spatial_info)
    if(!file.exists(spatial_info)) stop('path to spatial information does not exist')

    if(file_extension(spatial_info) %in% c('shp', 'geojson', 'wkt')) {
      spatial_info = terra::vect(spatial_info)
      spatial_info = evaluate_gpoly_spatvector(spatial_info)
      return(spatial_info)
    } else {
      spatial_info = data.table::fread(input = spatial_info, nThread = cores)
    }


    ## 1.2 data.frame-like input
  } else if (inherits(spatial_info, 'data.table')) {
    if(isTRUE(copy_dt)) spatial_info = data.table::copy(spatial_info)
  } else if(inherits(spatial_info, 'data.frame')) {
    spatial_info = data.table::setDT(spatial_info)

    ## 1.3 SpatVector input
  } else if(inherits(spatial_info, 'SpatVector')) {
    spatial_info = evaluate_gpoly_spatvector(spatial_info)
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

  # get unique IDs
  unique_IDs = spatial_info[, unique(poly_ID)]

  # 4. create spatvector
  spatial_info = dt_to_spatVector_polygon(spatial_info,
                                          include_values = TRUE)

  return_list = list(spatvector = spatial_info,
                     unique_IDs = unique_IDs)

  return(return_list)


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









# Feature info ####

#' @title Evaluate feature info
#' @name evaluate_feat_info
#' @description Evaluate spatial feature information input
#' @param spatial_feat_info spatial feature information to evaluate
#' @param cores how many cores to use
#' @param feat_ID feature IDs to check with
#' @return data.table
#' @keywords internal
#' @noRd
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




