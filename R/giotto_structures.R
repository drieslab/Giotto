

## ** cell shape polygons ####



#' @title Convert polygon to raster
#' @name polygon_to_raster
#' @description function to convert terra SpatVector Polygon shape into a terra SpatRaster
#' @keywords internal
polygon_to_raster = function(polygon, field = NULL) {

  pol_xmax = terra::xmax(polygon)
  pol_ymax = terra::ymax(polygon)
  r = terra::rast(polygon, ncols = pol_xmax, nrows = pol_ymax)

  if(is.null(field)) {
    field = names(polygon)[1]
  }

  poly_rast = terra::rasterize(x = polygon, r, field = field)

  return(poly_rast)

}



## extension of spatVector object
## name should match the cellular structure

#' @title Create a giotto polygon object
#' @name create_giotto_polygon_object
#' @keywords internal
create_giotto_polygon_object = function(name = 'cell',
                                        spatVector = NULL,
                                        spatVectorCentroids = NULL,
                                        overlaps = NULL) {


  # create minimum giotto
  g_polygon = giottoPolygon(name = name,
                            spatVector = NULL,
                            spatVectorCentroids = NULL,
                            overlaps = NULL)

  ## 1. check spatVector object
  if(!methods::is(spatVector, 'SpatVector')) {
    stop("spatVector needs to be a SpatVector object from the terra package")
  }

  g_polygon@spatVector = spatVector


  ## 2. centroids need to be of similar length as polygons
  if(!is.null(spatVectorCentroids)) {
    if(!methods::is(spatVectorCentroids, 'SpatVector')) {
      stop("spatVectorCentroids needs to be a spatVector object from the terra package")
    }

    l_centroids = nrow(terra::values(spatVectorCentroids))
    l_polygons = nrow(terra::values(spatVector))

    if(l_centroids == l_polygons) {
      g_polygon@spatVectorCentroids = spatVectorCentroids
    } else {
      stop('number of centroids does not equal number of polygons')
    }

  }

  ## 3. overlaps info
  g_polygon@overlaps = overlaps


  # provide name
  g_polygon@name = name

  # giotto polygon object
  return(g_polygon)
}





#' @title Identify background range polygons
#' @name identify_background_range_polygons
#' @description function to remove background polygon based on largest range
#' @keywords internal
identify_background_range_polygons = function(spatVector) {

  # define for data.table
  x = y = geom = V1 = NULL

  # identify polygon with the largest average range for x and y
  gDT = data.table::as.data.table(terra::geom(spatVector))

  range_geom_x = gDT[, max(x)-min(x), by = geom]
  range_geom_y = gDT[, max(y)-min(y), by = geom]
  range_geom = rbind(range_geom_x, range_geom_y)
  range_geom = range_geom[, mean(V1), by = geom]
  data.table::setorder(range_geom, -V1)

  # get original mask id for identified 'background' polygon
  backgr_polygon_id = range_geom[1, ][['geom']]
  values = terra::values(spatVector)
  poly_id = values[backgr_polygon_id, 1]

  return(poly_id)

}



#' @title Create segmentation polygons
#' @name create_segm_polygons
#' @description creates giotto polygons from segmentation mask data
#' @return giotto polygon
#' @keywords internal
create_segm_polygons = function(maskfile,
                                name = 'cell',
                                poly_IDs = NULL,
                                flip_vertical = TRUE,
                                shift_vertical_step = TRUE,
                                flip_horizontal = TRUE,
                                shift_horizontal_step = TRUE,
                                remove_background_polygon = FALSE) {


  if(!file.exists(maskfile)) {
    stop('path : ', maskfile, ' does not exist \n')
  }

  terra_rast = create_terra_spatRaster(maskfile)
  rast_dimensions = dim(terra_rast)

  terra_polygon = terra::as.polygons(x = terra_rast, value = TRUE)
  names(terra_polygon) = 'mask'


  ## flip axes ##
  if(flip_vertical == TRUE) {
    terra_polygon = terra::flip(terra_polygon, direction = 'vertical')
  }

  if(flip_horizontal == TRUE) {
    terra_polygon = terra::flip(terra_polygon, direction = 'horizontal')
  }

  ## shift values ##
  if(shift_vertical_step == TRUE) {
    shift_vertical_step = rast_dimensions[2]
  } else if(is.numeric(shift_vertical_step)) {
    shift_vertical_step = shift_vertical_step
  } else {
    shift_vertical_step = 0
  }
  if(shift_horizontal_step == TRUE) {
    shift_horizontal_step = rast_dimensions[1]
  } else if(is.numeric(shift_horizontal_step)) {
    shift_horizontal_step = shift_horizontal_step
  } else {
    shift_horizontal_step = 0
  }

  terra_polygon = terra::shift(terra_polygon,
                               dx = shift_horizontal_step,
                               dy = shift_vertical_step)

  # remove background polygon
  if(remove_background_polygon == TRUE) {

    mask_id = identify_background_range_polygons(terra_polygon)
    terra_polygon = terra::subset(x = terra_polygon, terra_polygon[['mask']] != mask_id)

  }

  # provide own cell_ID name
  if(!is.null(poly_IDs)) {
    if(length(poly_IDs) != nrow(terra::values(terra_polygon))) {
      stop('length poly_IDs does not equal number of found polygons \n')
    }
    terra_polygon$poly_ID = poly_IDs
  } else {
    terra_polygon$poly_ID = paste0(name, '_', 1:nrow(terra::values(terra_polygon)))
  }


  g_polygon = create_giotto_polygon_object(name = name,
                                           spatVector = terra_polygon,
                                           spatVectorCentroids = NULL)
  return(g_polygon)

}


#' @title Calculate polygon centroids
#' @name calculate_centroids_polygons
#' @description calculates centroids from selected polygons
#' @keywords internal
calculate_centroids_polygons = function(gpolygon,
                                        name = 'centroids',
                                        append_gpolygon = TRUE) {

  terra_polygon_centroids = terra::centroids(gpolygon@spatVector)

  if(append_gpolygon == TRUE) {

    gpolygon = create_giotto_polygon_object(name = gpolygon@name,
                                            spatVector = gpolygon@spatVector,
                                            spatVectorCentroids = terra_polygon_centroids)

  } else {

    gpolygon = create_giotto_polygon_object(name = name,
                                            spatVector = terra_polygon_centroids,
                                            spatVectorCentroids = NULL)

  }

  return(gpolygon)

}



#' @title Split multi-part polygons
#' @name fix_multipart_geoms
#' @description function to split geoms (polygons) that have multiple parts
#' @keywords internal
fix_multipart_geoms = function(spatVector) {

  # data.table variables
  x = y = geom = part = NULL

  spatVecDT = spatVector_to_dt(spatVector)
  uniq_multi = unique(spatVecDT[part == 2]$geom)

  # geoms to keep
  tokeepDT = spatVecDT[!geom %in% uniq_multi]
  tokeepDT = tokeepDT[,.(x, y, geom)]

  # rename
  total_geoms = length(unique(tokeepDT$geom))

  uniq_geom_vec = 1:total_geoms
  names(uniq_geom_vec) = unique(tokeepDT$geom)
  tokeepDT[, geom := uniq_geom_vec[[as.character(geom)]], by = 1:nrow(tokeepDT)]

  new_list = list()
  add_i = 1
  for(multi in uniq_multi) {

    tosplit = spatVecDT[geom == multi]

    intern_list = list()
    for(part_i in unique(tosplit$part)) {

      tempsplit = tosplit[part == part_i]
      tempsplit = tempsplit[,.(x,y,geom)]
      tempsplit[, geom := (total_geoms+add_i)]

      add_i = add_i + 1

      intern_list[[part_i]] = tempsplit

    }

    final_intern = do.call('rbind', intern_list)

    new_list[[multi]] = final_intern
  }

  final_new = do.call('rbind', new_list)

  finalDT = rbind(tokeepDT[,.(x, y, geom)], final_new)

  #return(finalDT)

  test = createGiottoPolygonsFromDfr(segmdfr = finalDT)

  return(test@spatVector)

}




## segmMaskToPolygon

#' @title Create giotto polygons from mask file
#' @name createGiottoPolygonsFromMask
#' @description Creates Giotto polygon object from a mask file (e.g. segmentation results)
#' @param maskfile path to mask file
#' @param mask_method how the mask file defines individual segmentation annotations
#' @param name name for polygons
#' @param remove_background_polygon try to remove background polygon (default: FALSE)
#' @param background_algo algorithm to remove background polygon
#' @param fill_holes fill holes within created polygons
#' @param poly_IDs unique names for each polygon in the mask file
#' @param flip_vertical flip mask figure in a vertical manner
#' @param shift_vertical_step shift vertical (boolean or numerical)
#' @param flip_horizontal flip mask figure in a horizontal manner
#' @param shift_horizontal_step shift horizontal (boolean or numerical)
#' @param calc_centroids calculate centroids for polygons
#' @param fix_multipart try to split polygons with multiple parts (default: TRUE)
#' @param remove_unvalid_polygons remove unvalid polygons (default: TRUE)
#' @return a giotto polygon object
#' @concept mask polygon
#' @export
createGiottoPolygonsFromMask = function(maskfile,
                                        mask_method = c('guess', 'single', 'multiple'),
                                        name = 'cell',
                                        remove_background_polygon = FALSE,
                                        background_algo = c('range'),
                                        fill_holes = TRUE,
                                        poly_IDs = NULL,
                                        flip_vertical = TRUE,
                                        shift_vertical_step = TRUE,
                                        flip_horizontal = TRUE,
                                        shift_horizontal_step = TRUE,
                                        calc_centroids = FALSE,
                                        fix_multipart = TRUE,
                                        remove_unvalid_polygons = TRUE) {

  # define for .()
  x = NULL
  y = NULL
  geom = NULL
  part = NULL

  # select background algo
  background_algo = match.arg(background_algo, choices = 'range')

  # check if mask file exists
  if(!file.exists(maskfile)) {
    stop('path : ', maskfile, ' does not exist \n')
  }

  # mask method
  # single: single mask value for all segmented cells
  # multiple: multiple mask values and thus a unique value for each segmented cell
  mask_method = match.arg(mask_method, choices = c('guess', 'single', 'multiple'))

  # create polygons from mask
  terra_rast = create_terra_spatRaster(maskfile)
  rast_dimensions = dim(terra_rast)
  terra_polygon = terra::as.polygons(x = terra_rast, value = TRUE)

  # fill holes
  if(fill_holes == TRUE) {
    terra_polygon = terra::fillHoles(terra_polygon)
  }

  # remove unvalid polygons
  if(remove_unvalid_polygons == TRUE) {
    valid_index = terra::is.valid(terra_polygon)
    terra_polygon = terra_polygon[valid_index]
  }


  spatVecDT = spatVector_to_dt(terra_polygon)

  # guess mask method
  if(mask_method == 'guess') {
    uniq_geoms = length(unique(spatVecDT$geom))
    uniq_parts = length(unique(spatVecDT$part))
    mask_method = ifelse(uniq_geoms > uniq_parts, 'multiple', 'single')
  }

  if(mask_method == 'multiple') {
    if(is.null(poly_IDs)) {
      spatVecDT[, geom := paste0(name,geom)]
    }
    g_polygon = createGiottoPolygonsFromDfr(segmdfr = spatVecDT[,.(x, y, geom)])
    terra_polygon = g_polygon@spatVector
  } else if(mask_method == 'single') {
    if(is.null(poly_IDs)) {
      spatVecDT[, part := paste0(name,part)]
    }
    g_polygon = createGiottoPolygonsFromDfr(segmdfr = spatVecDT[,.(x, y, part)])
    terra_polygon = g_polygon@spatVector
  }

  #names(terra_polygon) = 'mask'

  ## flip axes ##
  if(flip_vertical == TRUE) {
    terra_polygon = terra::flip(terra_polygon, direction = 'vertical')
  }

  if(flip_horizontal == TRUE) {
    terra_polygon = terra::flip(terra_polygon, direction = 'horizontal')
  }

  ## shift values ##
  if(shift_vertical_step == TRUE) {
    shift_vertical_step = rast_dimensions[1] # nrows of raster
  } else if(is.numeric(shift_vertical_step)) {
    shift_vertical_step = shift_vertical_step
  } else {
    shift_vertical_step = 0
  }
  if(shift_horizontal_step == TRUE) {
    shift_horizontal_step = rast_dimensions[2] # ncols of raster
  } else if(is.numeric(shift_horizontal_step)) {
    shift_horizontal_step = shift_horizontal_step
  } else {
    shift_horizontal_step = 0
  }

  print(shift_horizontal_step)
  print(shift_vertical_step)

  terra_polygon = terra::shift(terra_polygon,
                               dx = shift_horizontal_step,
                               dy = shift_vertical_step)


  # remove background polygon
  if(remove_background_polygon == TRUE) {

    if(background_algo == 'range') {
      backgr_poly_id = identify_background_range_polygons(terra_polygon)
      print(backgr_poly_id)
    }

    terra_polygon = terra::subset(x = terra_polygon, terra_polygon[['poly_ID']] != backgr_poly_id)

  }


  # provide own cell_ID name
  if(!is.null(poly_IDs)) {

    if(remove_unvalid_polygons == TRUE) {
      poly_IDs = poly_IDs[valid_index]
    }

    if(length(poly_IDs) != nrow(terra::values(terra_polygon))) {
      stop('length cell_IDs does not equal number of found polyongs \n')
    }
    terra_polygon$poly_ID = as.character(poly_IDs)
  } else {
    terra_polygon$poly_ID = paste0(name, '_', 1:nrow(terra::values(terra_polygon)))
  }


  g_polygon = create_giotto_polygon_object(name = name,
                                           spatVector = terra_polygon,
                                           spatVectorCentroids = NULL)

  # add centroids
  if(calc_centroids == TRUE) {
    g_polygon = calculate_centroids_polygons(gpolygon = g_polygon,
                                             name = 'centroids',
                                             append_gpolygon = TRUE)
  }


  return(g_polygon)

}




## segmDfrToPolygon

#' @title Create giotto polygons from dataframe
#' @name createGiottoPolygonsFromDfr
#' @description Creates Giotto polygon object from a structured dataframe-like object.
#' Three of the columns should correspond to x/y vertices and the polygon ID.
#' Additional columns are set as attributes
#' @param segmdfr data.frame-like object with polygon coordinate information (x, y, poly_ID)
#' with x and y being vertex information for the polygon referenced by poly_ID. See details
#' for how columns are selected for coordinate and ID information.
#' @param name name for the \code{giottoPolygon} object
#' @param calc_centroids (default FALSE) calculate centroids for polygons
#' @param verbose be verbose
#' @return giotto polygon object
#' @details When determining which column within the tabular data is intended to
#' provide polygon information, Giotto first checks the column names for 'x', 'y',
#' and 'poly_ID'. If any of these are discovered, they are directly selected. If
#' this is not discovered then Giotto checks the data type of the columns and selects
#' the first `'character'` type column to be 'poly_ID' and the first two `'numeric'`
#' columns as 'x' and 'y' respectively. If this is also unsuccessful then poly_ID
#' defaults to the 3rd column. 'x' and 'y' then default to the 1st and 2nd columns.
#' @concept polygon
#' @export
createGiottoPolygonsFromDfr = function(segmdfr,
                                       name = 'cell',
                                       calc_centroids = FALSE,
                                       verbose = TRUE) {

  # define for data.table
  geom = NULL

  input_dt = data.table::as.data.table(segmdfr)

  # data.frame like object needs to have 2 coordinate columns and
  # at least one other column as the feat_ID
  if(ncol(input_dt) < 3) stop('At minimum, columns for xy coordinates and poly ID are needed.\n')
  col_classes = sapply(input_dt, class)
  ## find poly_ID as either first character col or named column
  ## if neither exist, pick the 3rd column
  if('poly_ID' %in%  colnames(input_dt)) {
    poly_ID_col = which(colnames(input_dt) == 'poly_ID')
  } else {
    poly_ID_col = which(col_classes == 'character')
    if(length(poly_ID_col) < 1) poly_ID_col = 3 # case if no char found: default to 3rd
    else poly_ID_col = poly_ID_col[[1]] # case if char is found
  }

  if(isTRUE(verbose)) message(paste0('  Selecting col "',colnames(input_dt[, poly_ID_col, with = FALSE]),'" as poly_ID column'))
  colnames(input_dt)[poly_ID_col] = 'poly_ID'
  if(!inherits(input_dt$poly_ID, 'character')) {
    input_dt$poly_ID = as.character(input_dt$poly_ID) # ensure char
  }

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

  if(isTRUE(verbose)) message(paste0('  Selecting cols "',colnames(input_dt[, x_col, with = FALSE]),'" and "', colnames(input_dt[, y_col, with = FALSE]),'" as x and y respectively'))
  colnames(input_dt)[x_col] = 'x'
  colnames(input_dt)[y_col] = 'y'
  if(!inherits(input_dt$x, 'numeric')) {
    input_dt$x = as.numeric(input_dt$x) # ensure numeric
  }
  if(!inherits(input_dt$y, 'numeric')) {
    input_dt$y = as.numeric(input_dt$y) # ensure numeric
  }

  #pl = ggplot()
  #pl = pl + geom_polygon(data = input_dt[100000:200000], aes(x = x, y = y, group = poly_ID))
  #print(pl)

  # add other colnames for the input data.table
  nr_of_cells_vec = 1:length(unique(input_dt$poly_ID))
  names(nr_of_cells_vec) = unique(input_dt$poly_ID)
  new_vec = nr_of_cells_vec[as.character(input_dt$poly_ID)]
  input_dt[, geom := new_vec]

  input_dt[, c('part', 'hole') := list(1, 0)]
  input_dt = input_dt[, c('geom', 'part', 'x', 'y', 'hole', 'poly_ID'), with =F]


  #pl = ggplot()
  ##pl = pl + geom_polygon(data = input_dt[100000:200000], aes(x = x, y = y, group = geom))
  #print(pl)



  # create spatvector
  spatvector = dt_to_spatVector_polygon(input_dt,
                                        include_values = TRUE)

  hopla = spatVector_to_dt(spatvector)

  #pl = ggplot()
  #pl = pl + geom_polygon(data = hopla[100000:200000], aes(x = x, y = y, group = geom))
  #print(pl)



  g_polygon = create_giotto_polygon_object(name = name,
                                           spatVector = spatvector,
                                           spatVectorCentroids = NULL)

  # add centroids
  if(calc_centroids == TRUE) {
    g_polygon = calculate_centroids_polygons(gpolygon = g_polygon,
                                             name = 'centroids',
                                             append_gpolygon = TRUE)
  }

  return(g_polygon)

}



#' @title Extract list of polygons
#' @name extract_polygon_list
#' @description  to extract list of polygons
#' @keywords internal
extract_polygon_list = function(polygonlist,
                                polygon_mask_list_params,
                                polygon_dfr_list_params) {

  # if polygonlist is not a names list
  # try to make list and give default names
  if(!is.list(polygonlist)) {

    wrap_msg('polygonlist is not a list')

    polygonlist = as.list(polygonlist)
    if(length(polygonlist) == 1) {
      names(polygonlist) = 'cell'
    } else {
      polygonlist_l = length(polygonlist)
      names(polygonlist) = c('cell', paste0('info', 1:(polygonlist_l-1)))
    }
  }

  # if it is list
  # test if it has names
  if(is.null(names(polygonlist))) {

    wrap_msg('polygonlist is a list without names')

    if(length(polygonlist) == 1) {
      names(polygonlist) = 'cell'
    } else {
      polygonlist_l = length(polygonlist)
      names(polygonlist) = c('cell', paste0('info', 1:(polygonlist_l-1)))

    }
  }

  # make sure cell is one of the names
  #all_names = names(polygonlist)
  #if(!any('cell' %in% all_names)) {
  #  stop(" Information about 'cell' needs to be provided")
  #}


  final_list = list()

  for(poly_i in 1:length(polygonlist)) {

    name_polyinfo = names(polygonlist)[[poly_i]]
    polyinfo = polygonlist[[poly_i]]

    wrap_msg('Process polygon information for: ', name_polyinfo)

    if(is.character(polyinfo)) {
      poly_results = do.call('createGiottoPolygonsFromMask', c(name = name_polyinfo,
                                                               maskfile = polyinfo,
                                                               polygon_mask_list_params))

    } else if(inherits(polyinfo, 'data.frame')) {

      poly_results = do.call('createGiottoPolygonsFromDfr', c(name = name_polyinfo,
                                                              segmdfr = polyinfo,
                                                              polygon_dfr_list_params))

    } else if(inherits(polyinfo, 'giottoPolygon')) {

      poly_results = polyinfo
      name_polyinfo = polyinfo@name

    } else {

      stop('Polygon can only be extraxted from a mask file or from a correctly formatted data.frame')

    }

    final_list[[name_polyinfo]] = poly_results

  }

  print(final_list)
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
  if(!inherits(gobject, 'giotto')) {
    stop('gobject needs to be a giotto object')
  }

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




# helper functions

# convert spatVector to data.table

#' @title Convert spatVector to data.table
#' @name spatVector_to_dt
#' @description  convert spatVector to data.table
#' @keywords internal
spatVector_to_dt = function(spatvector,
                            include_values = TRUE) {

  # define for :=
  geom = NULL

  DT_geom = data.table::as.data.table(terra::geom(spatvector))

  if(include_values == TRUE) {
    DT_values = data.table::as.data.table(terra::values(spatvector))
    DT_values[, geom := 1:nrow(DT_values)]
    DT_full = data.table::merge.data.table(DT_geom, DT_values, by = 'geom')
    return(DT_full)
  } else {
    return(DT_geom)
  }
}


#' @title Convert data.table to polygon spatVector
#' @name dt_to_spatVector_polygon
#' @description convert data.table to spatVector for polygons
#' @keywords internal
dt_to_spatVector_polygon = function(dt,
                                    include_values = TRUE,
                                    specific_values = NULL) {


  all_colnames = colnames(dt)
  geom_values = c('geom', 'part', 'x', 'y', 'hole')
  other_values = all_colnames[!all_colnames %in% geom_values]

  if(include_values == TRUE) {

    if(!is.null(specific_values)) {
      other_values = other_values[other_values %in% specific_values]
    }


    spatVec = terra::vect(x = as.matrix(dt[,geom_values, with = FALSE]),
                          type = 'polygons', atts = unique(dt[,other_values, with = FALSE]))

  } else {

    spatVec = terra::vect(x = as.matrix(dt[,geom_values, with = FALSE]),
                          type = 'polygons', atts = NULL)

  }

  return(spatVec)

}


#' @title Convert spline to polygon
#' @name spline_poly
#' @description spline polynomial to smooth polygon
#' @param xy xy
#' @param vertices vertices
#' @param k k
#' @param ... additional params to pass
#' @keywords internal
spline_poly <- function(xy, vertices = 20, k = 3, ...) {
  # Assert: xy is an n by 2 matrix with n >= k.

  # Wrap k vertices around each end.
  n <- dim(xy)[1]
  if (k >= 1) {
    data <- rbind(xy[(n-k+1):n,], xy, xy[1:k, ])
  } else {
    data <- xy
  }

  # Spline the x and y coordinates.
  data.spline <- stats::spline(1:(n+2*k), data[,1], n=vertices, ...)
  x <- data.spline$x
  x1 <- data.spline$y
  x2 <- stats::spline(1:(n+2*k), data[,2], n=vertices, ...)$y

  # Retain only the middle part.
  cbind(x1, x2)[k < x & x <= n+k, ]
}




#' @title smoothGiottoPolygons
#' @name smoothGiottoPolygons
#' @description Smooths Giotto polygon object
#' @param gpolygon giotto polygon object
#' @param vertices number of vertices
#' @param k k
#' @param set_neg_to_zero set negative values to zero (default: TRUE)
#' @param ... additional params to pass to \code{spline}
#' @return Smoothed Giotto polygon object with reduced vertices
#' @concept polygon
#' @seealso \code{\link[stats]{spline}}
#' @export
smoothGiottoPolygons = function(gpolygon,
                                vertices = 20,
                                k = 3,
                                set_neg_to_zero = TRUE,
                                ...) {

  # define for .()
  x = NULL
  y = NULL

  # define for data.table [] subsetting
  geom = NULL

  polygDT = spatVector_to_dt(gpolygon@spatVector)

  # store other values
  all_colnames = colnames(polygDT)
  geom_values = c('geom', 'part', 'x', 'y', 'hole')
  other_values = all_colnames[!all_colnames %in% geom_values]
  other_values_uniq_dt = unique(polygDT[,c('geom', 'part', 'hole', other_values), with =F])

  # apply smoothing to each polygon
  comb = lapply(1:length(unique(polygDT$geom)), FUN = function(z) {

    polygMat = as.matrix(polygDT[geom == z,.(x, y)])

    # adjust k to maximum value
    max_k = nrow(polygMat)
    if(k >= max_k) {
      cat('k will be set to ', max_k)
      k = max_k
    }

    polygDT_smooth = data.table::as.data.table(spline_poly(polygMat, vertices = vertices, k = k, ...))
    polygDT_smooth[, geom := z]

  })
  comb_res = do.call('rbind', comb)


  # add other columns back
  comb_res = data.table::merge.data.table(comb_res, other_values_uniq_dt, by = 'geom')
  comb_res = comb_res[, c('geom', 'part', 'x1', 'x2', 'hole', other_values), with = F]
  colnames(comb_res)[3:4] = c('x', 'y')

  if(set_neg_to_zero == TRUE) {
    comb_res[, x := ifelse(x < 0, 0, x)]
    comb_res[, y := ifelse(y < 0, 0, y)]
  }

  new_spatvec = dt_to_spatVector_polygon(comb_res)

  for(ID in new_spatvec$poly_ID) {
    bool = terra::is.valid(new_spatvec[new_spatvec$poly_ID == ID])
    if(bool == FALSE) {
      print(ID)
      #plot(new_spatvec[new_spatvec$poly_ID == ID])
      #orig_spatvector = gpolygon@spatVector
      #new_spatvec[new_spatvec$poly_ID == ID] = orig_spatvector[orig_spatvector$poly_ID == ID]
    }
  }

  new_gpolygon = create_giotto_polygon_object(name = gpolygon@name,
                                              spatVector = new_spatvec,
                                              spatVectorCentroids = gpolygon@spatVectorCentroids)

  return(new_gpolygon)

}




### * simple polygon generation ####


#' @title Spatial polygons stamp
#' @name polyStamp
#' @description Takes a given stamp polygon and places it at each spatial location
#' provided.
#' @param stamp_dt data.table with x and y vertices for a polygon to be stamped.
#' Column names are expected to be 'x' and 'y' respectively
#' @param spatlocs spatial locations with x and y coordinates where polygons should
#' be stamped. Column names are expected to be 'sdimx' and 'sdimy' respectively.
#' @return returns a data.table of polygon vertices
#' @export
polyStamp = function(stamp_dt,
                     spatlocs) {

  # define polys relative to centroid
  stamp_centroid = c(x = mean(stamp_dt[['x']]),
                     y = mean(stamp_dt[['y']]))
  rel_vertices = data.table::data.table(x = stamp_dt$x - stamp_centroid[['x']],
                                        y = stamp_dt$y - stamp_centroid[['y']])

  # generate poly vertices around given spatlocs
  poly_dt = apply(X = spatlocs,
                  MARGIN = 1,
                  function(r) {
                    return(data.table::data.table(x = rel_vertices[['x']] + as.numeric(r[['sdimx']]),
                                                  y = rel_vertices[['y']] + as.numeric(r[['sdimy']]),
                                                  poly_ID = as.character(r[['cell_ID']])))
                  })

  return(do.call(rbind, poly_dt))

}


#' @title Generate circle polygon vertices
#' @name circleVertices
#' @description Generates vertex coordinates for a circle around (0,0) with the
#' given radius. Modified from \pkg{packcircles}.
#' @param radius radius of circle to be drawn
#' @param npoints number of vertices to generate
#' @seealso polyStamp rectVertices
#' @return a data.table of circle vertices
#' @export
circleVertices = function(radius,
                          npoints = 25) {
  a = seq(0, 2*pi, length.out = npoints + 1)
  x = radius * cos(a)
  y = radius * sin(a)
  m = data.table::data.table(x = x, y = y)
  return(m)
}


#' @title Generate rectangular polygon vertices
#' @name rectVertices
#' @description Generates vertex coordinates for a rectangle with dimensions given
#' through \code{dims} param.
#' @param dims named vector in the style of c(x = \code{numeric}, y = \code{numeric})
#' that defines the width (x) and height (y) of the generated rectangle polygon.
#' @seealso polyStamp circleVertices
#' @return a data.table of rectangle vertices
#' @export
rectVertices = function(dims) {
  if(length(dims) == 1) xdim = ydim = dims
  else xdim = dims[['x']] ; ydim = dims[['y']]

  m = data.table::data.table(x = c(0,0,xdim,xdim),
                             y = c(0,ydim,ydim,0))
  return(m)
}



## ** feature points ####



#' @title Create feature network object
#' @name create_featureNetwork_object
#' @param name name to assign the created feature network object
#' @param network_datatable network data.table object
#' @param network_lookup_id network lookup id
#' @param full fully connected status
#' @keywords internal
create_featureNetwork_object = function(name = 'feat_network',
                                        network_datatable = NULL,
                                        network_lookup_id = NULL,
                                        full = NULL) {


  # create minimum giotto points object
  f_network = featureNetwork(name = name,
                             network_datatable = NULL,
                             network_lookup_id = NULL,
                             full = NULL)

  ## 1. check network data.table object
  if(!methods::is(network_datatable, 'data.table')) {
    stop("network_datatable needs to be a network data.table object")
  }
  f_network@network_datatable = network_datatable

  ## 2. provide network fully connected status
  f_network@full = full

  ## 3. provide feature network name
  f_network@name = name

  ## 4. network lookup id
  f_network@network_lookup_id = network_lookup_id

  # giotoPoints object
  return(f_network)

}




#' @title Create giotto points object
#' @name create_giotto_points_object
#' @param feat_type feature type
#' @param spatVector terra spatVector object containing point data
#' @param networks feature network object
#' @keywords internal
create_giotto_points_object = function(feat_type = 'rna',
                                       spatVector = NULL,
                                       networks = NULL) {


  # create minimum giotto points object
  g_points = giottoPoints(feat_type = feat_type,
                          spatVector = NULL,
                          networks = NULL)

  ## 1. check terra spatVector object
  if(!inherits(spatVector, 'SpatVector')) {
    stop("spatVector needs to be a spatVector object from the terra package")
  }

  g_points@spatVector = spatVector

  ## 2. provide feature id
  g_points@feat_type = feat_type

  ## 3. feature_network object
  g_points@networks = networks

  # giotoPoints object
  return(g_points)

}



#' @title Create terra spatvector object from a data.frame
#' @name create_spatvector_object_from_dfr
#' @description create terra spatvector from a data.frame where cols 1 and 2 must
#' be x and y coordinates respectively. Additional columns are set as attributes
#' to the points where the first additional (col 3) should be the feat_ID.
#' @param x data.frame object
#' @param verbose be verbose
#' @keywords internal
create_spatvector_object_from_dfr = function(x,
                                             verbose = TRUE) {


  x = data.table::as.data.table(x)

  # data.frame like object needs to have 2 coordinate columns and
  # at least one other column as the feat_ID
  if(ncol(x) < 3) stop('At minimum, columns for xy coordinates and feature ID are needed.\n')
  col_classes = sapply(x, class)
  ## find feat_ID as either first character col or named column
  ## if not detected, select 3rd column
  if('feat_ID' %in% colnames(x)) {
    feat_ID_col = which(colnames(x) == 'feat_ID')
  } else {
    feat_ID_col = which(col_classes == 'character')
    if(length(feat_ID_col) < 1) feat_ID_col = 3 # case if no char found: default to 3rd
    else feat_ID_col = feat_ID_col[[1]] # case if char is found
  }

  if(isTRUE(verbose)) message(paste0('  Selecting col "',colnames(x[, feat_ID_col, with = FALSE]),'" as feat_ID column'))
  colnames(x)[feat_ID_col] = 'feat_ID'
  if(!inherits(x$feat_ID, 'character')) {
    x$feat_ID = as.character(x$feat_ID) # ensure char
  }

  ## find first two numeric cols as x and y respectively or named column
  ## if not detected select 1st and 2nd cols for x and y respectively
  if(all(c('x','y') %in% colnames(x))) {
    x_col = which(colnames(x) == 'x')
    y_col = which(colnames(x) == 'y')
  } else {
    x_col = which(col_classes == 'numeric')
    if(length(x_col) < 2) x_col = 1 # case if no/too few num found: default to 1st
    else x_col = x_col[[1]] # case if num found
    y_col = which(col_classes == 'numeric')
    if(length(y_col) < 2) y_col = 2 # case if no/too few num found: default to 2nd
    else y_col = y_col[[2]] # case if num found
  }


  if(isTRUE(verbose)) message(paste0('  Selecting cols "',colnames(x[, x_col, with = FALSE]),'" and "', colnames(x[, y_col, with = FALSE]),'" as x and y respectively'))
  colnames(x)[x_col] = 'x'
  colnames(x)[y_col] = 'y'
  if(!inherits(x$x, 'numeric')) x$x = as.numeric(x$x) # ensure numeric
  if(!inherits(x$y, 'numeric')) x$y = as.numeric(x$y) # ensure numeric


  ## select location and attribute dataframes
  # Use unique() to set column order
  ordered_colnames = unique(c('feat_ID','x','y', colnames(x)))
  x = x[, ordered_colnames, with = FALSE]
  loc_dfr = x[,2:3]
  att_dfr = x[,-c(2:3)]

  spatvec = terra::vect(as.matrix(loc_dfr), type = 'points', atts = att_dfr)

  # will be given and is a unique numerical barcode for each feature
  spatvec[['feat_ID_uniq']] = 1:nrow(spatvec)

  return(spatvec)

}


# create Giotto points from data.frame or spatVector

#' @title Create giotto points object
#' @name createGiottoPoints
#' @description Creates Giotto point object from a structured dataframe-like object
#' @param x spatVector or data.frame-like object with points coordinate information (x, y, feat_ID)
#' @param feat_type feature type
#' @param verbose be verbose
#' @return giottoPoints
#' @concept polygon
#' @export
createGiottoPoints = function(x,
                              feat_type = 'rna',
                              verbose = TRUE) {

  if(inherits(x, 'data.frame')) {

    spatvec = create_spatvector_object_from_dfr(x = x,
                                                verbose = verbose)
    g_points = create_giotto_points_object(feat_type = feat_type,
                                           spatVector = spatvec)

  } else if(inherits(x, 'spatVector')) {

    g_points = create_giotto_points_object(feat_type = feat_type,
                                           spatVector = x)

  } else {

    stop('Class ', class(x), ' is not supported')

  }

  return(g_points)

}


# data.table to spatVector


#' @title Convert point data data.table to spatVector
#' @name dt_to_spatVector_points
#' @description  data.table to spatVector for points
#' @param dt data.table
#' @param include_values include additional values from data.table as attributes paired with created terra spatVector [boolean]
#' @param specific_values specific values to include as attributes if include_values == TRUE
#' @keywords internal
dt_to_spatVector_points = function(dt,
                                   include_values = TRUE,
                                   specific_values = NULL) {


  all_colnames = colnames(dt)
  geom_values = c('geom', 'part', 'x', 'y', 'hole')
  other_values = all_colnames[!all_colnames %in% geom_values]

  if(include_values == TRUE) {

    if(!is.null(specific_values)) {
      other_values = other_values[other_values %in% specific_values]
    }


    spatVec = terra::vect(x = as.matrix(dt[,geom_values, with = F]),
                          type = 'points', atts = dt[,other_values, with = F])

  } else {

    spatVec = terra::vect(x = as.matrix(dt[,geom_values, with = F]),
                          type = 'points', atts = NULL)

  }

  return(spatVec)

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



#' @title Extract list of giotto points objects
#' @name extract_points_list
#' @description to extract list of giotto points
#' @param pointslist list of inputs from which to create giotto points objects
#' @keywords internal
extract_points_list = function(pointslist) {

  # if pointslist is not a named list
  # try to make list and give default names
  if(!is.list(pointslist)) {
    pointslist = as.list(pointslist)
    if(length(pointslist) == 1) {
      names(pointslist) = 'rna'
    } else {
      pointslist_l = length(pointslist)
      names(pointslist) = c('rna', paste0('feat', 1:(pointslist_l-1)))
    }
  }

  # if it is list
  # test if it has names
  if(is.null(names(pointslist))) {

    if(length(pointslist) == 1) {
      names(pointslist) = 'rna'
    } else {
      pointslist_l = length(pointslist)
      names(pointslist) = c('rna', paste0('feat', 1:(pointslist_l-1)))

    }
  }

  # make sure rna is one of the names
  #all_names = names(pointslist)
  #if(!any('rna' %in% all_names)) {
  #  stop(" Information about 'rna' needs to be provided")
  #}


  final_list = list()

  for(point_i in 1:length(pointslist)) {

    name_pointinfo = names(pointslist)[[point_i]]
    pointinfo = pointslist[[point_i]]


    if(inherits(pointinfo, 'giottoPoints')) {

      point_results = pointinfo
      name_pointinfo = pointinfo@feat_type

    } else if(inherits(pointinfo, 'data.frame')) {

      point_results = createGiottoPoints(x = pointinfo,
                                         feat_type = name_pointinfo)

    } else if(inherits(pointinfo, 'character')) {

      stop('Giotto points can not yet be created directly from a file path')

    } else {

      stop('Giotto points can only be created from a correctly formatted data.frame-like object')

    }


    final_list[[name_pointinfo]] = point_results

  }

  return(final_list)


}



#' @title Create kNN spatial feature network using dbscan
#' @name createSpatialFeaturesKNNnetwork_dbscan
#' @description  to create a feature kNN spatial network using dbscan
#' @param gobject giotto object
#' @param feat_type feature type
#' @param name name to assign generated feature network
#' @param k number of neighbors for kNN to find
#' @param maximum_distance network maximum distance allowed
#' @param minimum_k minimum neighbors allowed
#' @param add_feat_ids whether to add feature information [boolean]
#' @param verbose be verbose
#' @param ... additional parameters to pass to \code{\link[dbscan]{kNN}}
#' @keywords internal
createSpatialFeaturesKNNnetwork_dbscan = function(gobject,
                                                  feat_type = NULL,
                                                  name = "knn_feats_network",
                                                  k = 4,
                                                  maximum_distance = NULL,
                                                  minimum_k = 0,
                                                  add_feat_ids = FALSE,
                                                  verbose = TRUE,
                                                  ...) {

  # define for data.table
  from_feat = from = to_feat = to = from_to_feat = NULL

  ## 1. specify feat_type
  if(is.null(feat_type)) {
    gobject@feat_info[[1]]@feat_type
  }

  ## 2. get spatial feature info and convert to matrix
  if(verbose == TRUE) cat('Convert feature spatial info to matrix \n')
  featDT = spatVector_to_dt(gobject@feat_info[[feat_type]]@spatVector)
  spatial_locations_matrix = as.matrix(featDT[, c('x', 'y', NULL), with = F])

  # store lookup table to keep information about unique ID
  # important with multiple joined objects where row id is not always equal to unique gene
  network_id_lookup_table = data.table::data.table(row = 1:nrow(featDT),
                                                   id = featDT$feat_ID_uniq)

  ## 3. create kNN network
  if(verbose == TRUE) cat('Create kNN network with dbscan \n')
  knn_spatial = dbscan::kNN(x = spatial_locations_matrix,
                            k = k,
                            ...)

  knn_sptial.norm = data.table::data.table(from = rep(1:nrow(knn_spatial$id), k),
                                           to = as.vector(knn_spatial$id),
                                           #weight = 1/(1 + as.vector(knn_spatial$dist)),
                                           distance = as.vector(knn_spatial$dist))

  ## 3. keep minimum and filter
  if(verbose == TRUE) cat('Filter output for distance and minimum neighbours \n')
  knn_sptial.norm[, rank := 1:.N, by = 'from']

  if(minimum_k != 0) {
    filter_bool = knn_sptial.norm$rank <= minimum_k
  } else {
    filter_bool = rep(TRUE, nrow(knn_sptial.norm))
  }


  if(!is.null(maximum_distance)) {
    maximum_distance_bool = knn_sptial.norm$distance <= maximum_distance
    filter_bool = filter_bool + maximum_distance_bool
    filter_bool[filter_bool > 0] = 1
    filter_bool = as.logical(filter_bool)
  }


  knn_sptial.norm = knn_sptial.norm[filter_bool]

  ## 3. add feature information and sort
  if(add_feat_ids == TRUE) {

    if(verbose == TRUE) cat('Add feat IDs and sort output \n')

    featDT_vec = featDT$feat_ID; names(featDT_vec) = featDT$feat_ID_uniq

    knn_sptial.norm[, from_feat := featDT_vec[from]]
    knn_sptial.norm[, to_feat := featDT_vec[to]]
    knn_sptial.norm[, from_to_feat := paste0(from_feat,'--',to_feat)]

    knn_sptial.norm = sort_combine_two_DT_columns(DT = knn_sptial.norm,
                                                  column1 = 'from_feat', column2 = 'to_feat',
                                                  myname = 'comb_feat')
  }


  knn_sptial.norm_object = create_featureNetwork_object(name = name,
                                                        network_datatable = knn_sptial.norm,
                                                        network_lookup_id = network_id_lookup_table,
                                                        full = FALSE)

  return(knn_sptial.norm_object)

}





#' @title Create kNN spatial feature network
#' @name createSpatialFeaturesKNNnetwork
#' @description Calculates the centroid locations for the giotto polygons
#' @param gobject giotto object
#' @param method kNN algorithm method
#' @param feat_type feature type to build feature network
#' @param name name of network
#' @param k number of neighbors
#' @param maximum_distance maximum distance bewteen features
#' @param minimum_k minimum number of neighbors to find
#' @param add_feat_ids add feature id names (default = FALSE, increases object size)
#' @param verbose be verbose
#' @param return_gobject return giotto object (default: TRUE)
#' @param toplevel_params toplevel value to pass when updating giotto params
#' @param ... additional parameters to pass to \code{\link[dbscan]{kNN}}
#' @return If \code{return_gobject = TRUE} a giotto object containing the network
#'   will be returned. If \code{return_gobject = FALSE} the network will be returned
#'   as a datatable.
#' @concept feature
#' @export
createSpatialFeaturesKNNnetwork = function(gobject,
                                           method = c('dbscan'),
                                           feat_type = NULL,
                                           name = "knn_feats_network",
                                           k = 4,
                                           maximum_distance = NULL,
                                           minimum_k = 0,
                                           add_feat_ids = FALSE,
                                           verbose = TRUE,
                                           return_gobject = TRUE,
                                           toplevel_params = 2,
                                           ...) {


  # 1. select feat_type
  if(is.null(feat_type)) {
    feat_type = gobject@expression_feat[[1]]
  }

  # 2. select method
  method = match.arg(method, choices = c('dbscan'))


  if(method == 'dbscan') {

    knn_feat_network_obj = createSpatialFeaturesKNNnetwork_dbscan(gobject = gobject,
                                                                  feat_type = feat_type,
                                                                  name = name,
                                                                  k = k,
                                                                  maximum_distance = maximum_distance,
                                                                  minimum_k = minimum_k,
                                                                  add_feat_ids = add_feat_ids,
                                                                  verbose = verbose,
                                                                  ...)
  }




  if(return_gobject == TRUE) {

    network_names = names(gobject@feat_info[[feat_type]]@networks)

    if(name %in% network_names) {
      cat('\n ', name, ' has already been used, will be overwritten \n')

    }

    gobject@feat_info[[feat_type]]@networks[[name]] = knn_feat_network_obj


    ## update parameters used ##
    gobject = update_giotto_params(gobject,
                                   description = '_featNetwork',
                                   return_gobject = TRUE,
                                   toplevel = toplevel_params)
    return(gobject)


  } else {
    return(knn_feat_network_obj@network_datatable)
  }

}






## * ####
## ** giotto structure functions ####


#' @title addSpatialCentroidLocationsLayer
#' @name addSpatialCentroidLocationsLayer
#' @description Calculates the centroid locations for the polygons within one selected layer
#' @param gobject giotto object
#' @param poly_info polygon information
#' @param feat_type feature type
#' @param spat_loc_name name to give to the created spatial locations
#' @param return_gobject return giotto object (default: TRUE)
#' @return If \code{return_gobject = TRUE} the giotto object containing the calculated
#'   polygon centroids will be returned. If \code{return_gobject = FALSE} only the
#'   generated polygon centroids will be returned.
#' @concept centroid
#' @export
addSpatialCentroidLocationsLayer = function(gobject,
                                            poly_info = 'cell',
                                            feat_type = NULL,
                                            spat_loc_name = 'raw',
                                            return_gobject = TRUE) {

  # define for .()
  x = NULL
  y = NULL
  poly_ID = NULL

  # Set feat_type and spat_unit
  poly_info = set_default_spat_unit(gobject = gobject,
                                    spat_unit = poly_info)
  # feat_type = set_default_feat_type(gobject = gobject,
  #                                   spat_unit = poly_info,
  #                                   feat_type = feat_type)
  feat_type = slot(gobject, 'expression_feat')[[1]] # Specifically preferable over set_default function
  #There may be no existing data in expression slot to find feat_type nesting from

  gpoly = get_polygon_info(gobject, polygon_name = poly_info, return_giottoPolygon = TRUE)

  extended_spatvector = calculate_centroids_polygons(gpolygon = gpoly,
                                                     name = 'centroids',
                                                     append_gpolygon = TRUE)

  centroid_spatvector = spatVector_to_dt(extended_spatvector@spatVectorCentroids)

  # this could be 3D
  spatial_locs = centroid_spatvector[, .(x, y, poly_ID)]
  colnames(spatial_locs) = c('sdimx', 'sdimy', 'cell_ID')

  if(return_gobject == TRUE) {

    # spatial location
    spat_locs_names = list_spatial_locations_names(gobject,
                                                   spat_unit = poly_info)
    if(spat_loc_name %in% spat_locs_names) {
      message('spatial locations for polygon information layer ', poly_info,
          ' and name ', spat_loc_name, ' already exists and will be replaced\n')
    }

    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
    gobject = set_spatial_locations(gobject = gobject,
                                    spat_unit = poly_info,
                                    spat_loc_name = spat_loc_name,
                                    spatlocs = spatial_locs)
    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###


    # cell ID
    gpoly_IDs = gpoly@spatVector$poly_ID
    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
    gobject = set_cell_id(gobject,
                          spat_unit = poly_info,
                          cell_IDs = gpoly_IDs)
    # gobject@cell_ID[[poly_info]] = gobject@spatial_info[[poly_info]]@spatVector$poly_ID
    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###


    # cell metadata
    # new spatial locations come with new cell and feature metadata
    for(type in feat_type) {
      cm = create_cell_meta_obj(metaDT = data.table::data.table(cell_ID = gpoly_IDs),
                                col_desc = c('unique IDs for each cell'),
                                spat_unit = poly_info,
                                feat_type = type,
                                provenance = poly_info)

      gfeat_IDs = get_feat_id(gobject, feat_type = type)
      fm = create_feat_meta_obj(metaDT = data.table::data.table(feat_ID = gfeat_IDs),
                                col_desc = c('unique IDs for each feature'),
                                spat_unit = poly_info,
                                feat_type = type,
                                provenance = poly_info)

      ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
      gobject = set_cell_metadata(gobject, metadata = cm, verbose = FALSE)
      gobject = set_feature_metadata(gobject, metadata = fm, verbose = FALSE)
      # gobject@cell_metadata[[poly_info]][[type]] = data.table::data.table(cell_ID = gobject@spatial_info[[poly_info]]@spatVector$poly_ID)
      # gobject@feat_metadata[[poly_info]][[type]] = data.table::data.table(feat_ID = gobject@feat_ID[[type]])
      ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
    }


    # add centroids information
    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
    gobject = set_polygon_info(gobject,
                               polygon_name = poly_info,
                               gpolygon = extended_spatvector,
                               verbose = FALSE)
    # gobject@spatial_info[[poly_info]] = extended_spatvector
    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###


    return(gobject)

  } else {
    return(spatial_locs)
  }

}


#' @title addSpatialCentroidLocations
#' @name addSpatialCentroidLocations
#' @description Calculates the centroid locations for the polygons within one or more selected layers
#' @param gobject giotto object
#' @param poly_info polygon information
#' @param feat_type feature type
#' @param spat_loc_name name to give to the created spatial locations
#' @param return_gobject return giotto object (default: TRUE)
#' @param verbose be verbose
#' @return If \code{return_gobject = TRUE} the giotto object containing the calculated
#'   polygon centroids will be returned. If \code{return_gobject = FALSE} only the
#'   generated polygon centroids will be returned.
#' @concept centroid
#' @export
addSpatialCentroidLocations = function(gobject,
                                       poly_info = 'cell',
                                       feat_type = NULL,
                                       spat_loc_name = 'raw',
                                       return_gobject = TRUE,
                                       verbose = TRUE) {


  potential_polygon_names = list_spatial_info_names(gobject)

  return_list = list()

  for(poly_layer in unique(poly_info)) {

    if(!poly_layer %in% potential_polygon_names) {
      warning('Polygon info layer with name ', poly_layer, ' has not been found and will be skipped')
    } else {

      if(verbose == TRUE) {
        cat('Start centroid calculation for polygon information layer: ', poly_layer, '\n')
      }

      if(return_gobject == TRUE) {
        gobject = addSpatialCentroidLocationsLayer(gobject = gobject,
                                                   poly_info = poly_layer,
                                                   feat_type = feat_type,
                                                   spat_loc_name = spat_loc_name,
                                                   return_gobject = return_gobject)
      } else {

        return_list[[poly_layer]] = addSpatialCentroidLocationsLayer(gobject = gobject,
                                                                     poly_info = poly_layer,
                                                                     feat_type = feat_type,
                                                                     spat_loc_name = spat_loc_name,
                                                                     return_gobject = return_gobject)

      }

    }
  }

  if(return_gobject == TRUE) {
    return(gobject)
  } else {
    return_list
  }

}





## * calculate overlap between cellular structures and features ####


## ** raster way ####


#' @name polygon_to_raster
#' @description  convert polygon to raster
#' @keywords internal
polygon_to_raster = function(polygon, field = NULL) {

  pol_xmax = terra::xmax(polygon)
  pol_xmin = terra::xmin(polygon)
  ncols = abs(pol_xmax-pol_xmin)

  pol_ymax = terra::ymax(polygon)
  pol_ymin = terra::ymin(polygon)
  nrows = abs(pol_ymax-pol_ymin)

  r = terra::rast(polygon, ncols = ncols, nrows = nrows)

  if(is.null(field)) {
    field = names(polygon)[1]
  }

  # ensure that field is numerical
  polygon$poly_i = 1:nrow(unique(polygon[[field]]))
  poly_rast = terra::rasterize(x = polygon, r, field = 'poly_i')

  poly_ID_vector = polygon[[field]][,1]; names(poly_ID_vector) = polygon[['poly_i']][,1]

  return(list('raster' = poly_rast, 'ID_vector' = poly_ID_vector))

}



#' @title calculateOverlapRaster
#' @name calculateOverlapRaster
#' @description calculate overlap between cellular structures (polygons) and features (points)
#' @param gobject giotto object
#' @param name_overlap name for the overlap results (default to feat_info parameter)
#' @param spatial_info polygon information
#' @param poly_ID_names (optional) list of poly_IDs to use
#' @param feat_info feature information
#' @param feat_subset_column feature info column to subset features with
#' @param feat_subset_ids ids within feature info column to use for subsetting
#' @param return_gobject return giotto object (default: TRUE)
#' @param verbose be verbose
#' @return giotto object or spatVector with overlapping information
#' @details Serial overlapping function.
#' @concept overlap
#' @export
calculateOverlapRaster = function(gobject,
                                  name_overlap = NULL,
                                  spatial_info = NULL,
                                  poly_ID_names = NULL,
                                  feat_info = NULL,
                                  feat_subset_column = NULL,
                                  feat_subset_ids = NULL,
                                  return_gobject = TRUE,
                                  verbose = TRUE) {

  # define for :=
  poly_ID = NULL
  poly_i = NULL
  ID = NULL
  x = NULL
  y = NULL
  feat_ID = NULL
  feat_ID_uniq = NULL

  # set defaults if not provided
  if(is.null(feat_info)) {
    feat_info = names(gobject@feat_info)[[1]]
  }

  if(is.null(name_overlap)) {
    name_overlap = feat_info
  }

  if(is.null(spatial_info)) {
    spatial_info = names(gobject@spatial_info)[[1]]
  }


  # spatial vector
  if(verbose) cat('1. convert polygon to raster \n')
  spatvec = gobject@spatial_info[[spatial_info]]@spatVector

  # subset spatvec
  if(!is.null(poly_ID_names)) {
    spatvec = spatvec[spatvec$poly_ID %in% poly_ID_names]
  }

  # spatial vector to raster
  spatrast_res = polygon_to_raster(spatvec, field = 'poly_ID')
  spatrast = spatrast_res[['raster']]
  ID_vector = spatrast_res[['ID_vector']]

  # point vector
  pointvec = gobject@feat_info[[feat_info]]@spatVector

  # subset points if needed
  # e.g. to select transcripts within a z-plane
  if(!is.null(feat_subset_column) & !is.null(feat_subset_ids)) {
    bool_vector = pointvec[[feat_subset_column]][[1]] %in% feat_subset_ids
    pointvec = pointvec[bool_vector]
  }

  ## overlap between raster and point
  if(verbose) cat('2. overlap raster and points \n')
  overlap_test = terra::extract(x = spatrast, y = pointvec)

  # add poly_ID information
  if(verbose) cat('3. add polygon information \n')
  overlap_test_dt = data.table::as.data.table(overlap_test)
  overlap_test_dt[, poly_ID := ID_vector[poly_i]]

  # add point information
  if(verbose) cat('4. add points information \n')
  pointvec_dt = spatVector_to_dt(pointvec)

  pointvec_dt_x = pointvec_dt$x ; names(pointvec_dt_x) = pointvec_dt$geom
  pointvec_dt_y = pointvec_dt$y ; names(pointvec_dt_y) = pointvec_dt$geom
  pointvec_dt_feat_ID = pointvec_dt$feat_ID ; names(pointvec_dt_feat_ID) = pointvec_dt$geom
  pointvec_dt_feat_ID_uniq = pointvec_dt$feat_ID_uniq ; names(pointvec_dt_feat_ID_uniq) = pointvec_dt$geom

  overlap_test_dt[, x := pointvec_dt_x[ID]]
  overlap_test_dt[, y := pointvec_dt_y[ID]]
  overlap_test_dt[, feat_ID := pointvec_dt_feat_ID[ID]]
  overlap_test_dt[, feat_ID_uniq := pointvec_dt_feat_ID_uniq[ID]]

  if(verbose) cat('5. create overlap polygon information \n')
  overlap_test_dt_spatvector = terra::vect(x = as.matrix(overlap_test_dt[, c('x', 'y'), with = F]),
                                           type = "points",
                                           atts = overlap_test_dt[, c('poly_ID', 'feat_ID', 'feat_ID_uniq'), with = F])
  names(overlap_test_dt_spatvector) = c('poly_ID', 'feat_ID', 'feat_ID_uniq')



  if(return_gobject == TRUE) {
    if(is.null(name_overlap)) {
      name_overlap = feat_info
    }

    gobject@spatial_info[[spatial_info]]@overlaps[[name_overlap]] = overlap_test_dt_spatvector
    return(gobject)

  } else {
    return(overlap_test_dt_spatvector)
  }

}



#' @title Overlap points -- single polygon
#' @name overlap_points_single_polygon
#' @description  overlap for a single polygon
#' @keywords internal
overlap_points_single_polygon = function(spatvec,
                                         poly_ID_name,
                                         pointvec_dt) {

  # define for data.table
  x = y = NULL

  ## extract single polygon and get spatextent
  one_polygon_spatvector = spatvec[spatvec$poly_ID == poly_ID_name]
  ext_limits = terra::ext(one_polygon_spatvector)

  ## extract potential features (points) based on limits
  one_polygon_pointvec_dt = pointvec_dt[x > terra::xmin(ext_limits) & x < terra::xmax(ext_limits)][y > terra::ymin(ext_limits) & y < terra::ymax(ext_limits)]
  one_polygon_pointsvec = dt_to_spatVector_points(one_polygon_pointvec_dt)

  ## calculate intersection between single polygon and points
  one_polygon_overlap = terra::intersect(x = one_polygon_spatvector,
                                         y = one_polygon_pointsvec)

  if(nrow(one_polygon_overlap) > 0) {
    return(one_polygon_overlap)
  } else {
    return(NULL)
  }

}





#' @title calculateOverlapPolygonImages
#' @name calculateOverlapPolygonImages
#' @description calculate overlap between cellular structures (polygons) and images (intensities)
#' @param gobject giotto object
#' @param name_overlap name for the overlap results (default to feat_info parameter)
#' @param spatial_info polygon information
#' @param poly_ID_names (optional) list of poly_IDs to use
#' @param image_names names of the images with raw data
#' @param return_gobject return giotto object (default: TRUE)
#' @param verbose be verbose
#' @return giotto object or data.table with overlapping information
#' @concept overlap
#' @export
calculateOverlapPolygonImages = function(gobject,
                                         name_overlap = 'images',
                                         spatial_info = 'cell',
                                         poly_ID_names = NULL,
                                         image_names = NULL,
                                         return_gobject = TRUE,
                                         verbose = TRUE) {

  if(is.null(image_names)) {
    stop('image_names = NULL, you need to provide the names of the images you want to use,
          see showGiottoImageNames() for attached images')
  }

  ## get polygon information
  poly_info = get_polygon_info(gobject = gobject,
                               polygon_name = spatial_info,
                               return_giottoPolygon = T)


  # calculate centroids for poly_info if not present
  if(is.null(poly_info@spatVectorCentroids)) {
    poly_info = calculate_centroids_polygons(gpolygon = poly_info,
                                             name = 'centroids',
                                             append_gpolygon = T)
  }


  potential_large_image_names = list_images_names(gobject, img_type = 'largeImage')

  for(img_name in image_names) {

    if(verbose) cat('Start process for image: ', img_name, '\n')

    if(!img_name %in% potential_large_image_names) {
      warning('image with the name ', img_name, ' was not found and will be skipped \n')
    } else {

      intensity_image = get_giottoLargeImage(gobject = gobject, name = img_name)

      extract_intensity = data.table::as.data.table(terra::extract(x = intensity_image@raster_object,
                                                            y = poly_info@spatVector))

      poly_ID_vector = poly_info@spatVector$poly_ID
      names(poly_ID_vector) = 1:length(poly_ID_vector)
      extract_intensity[, ID := poly_ID_vector[ID]]
      colnames(extract_intensity)[2] = img_name

      poly_info@overlaps[[name_overlap]][[img_name]] = extract_intensity

      return_list = list()

      if(return_gobject) {

        gobject = set_polygon_info(gobject = gobject,
                                   polygon_name = spatial_info,
                                   gpolygon = poly_info)
      } else {

        return_list[[img_name]] = pol_infoy

      }

    }
  }

  if(return_gobject) {
    return(gobject)
  } else {
    return(return_list)
  }
}









## ** polygon way ####


#' @title Overlap points per polgyon
#' @name overlap_points_per_polygon
#' @description Loop to overlap each single polygon
#' @keywords internal
#' @seealso \code{\link{overlap_points_single_polygon}}
overlap_points_per_polygon = function(spatvec,
                                      pointvec,
                                      poly_ID_names,
                                      verbose = TRUE) {

  # spatial polygon
  spatvec = spatvec[terra::is.valid(spatvec)]

  # points polygon
  pointvec_dt = spatVector_to_dt(pointvec)

  # get polygon names
  unique_cell_names = unique(spatvec$poly_ID)
  poly_ID_names = poly_ID_names[poly_ID_names %in% unique_cell_names]


  #final_vect = terra::vect()
  final_list = list(); i = 1
  for(poly_ID_name in poly_ID_names) {

    if(verbose == TRUE) print(poly_ID_name)

    result = overlap_points_single_polygon(spatvec = spatvec,
                                           poly_ID_name = poly_ID_name,
                                           pointvec_dt = pointvec_dt)

    if(!is.null(result)) {
      final_list[[i]] = result
      i = i+1
      #final_vect = rbind(final_vect, result)
    }


  }

  final_vect = do.call('rbind', final_list)

  return(final_vect)

}



#' @title calculateOverlapSerial
#' @name calculateOverlapSerial
#' @description calculate overlap between cellular structures (polygons) and features (points)
#' @param gobject giotto object
#' @param name_overlap name for the overlap results (default to feat_info parameter)
#' @param spatial_info polygon information
#' @param feat_info feature information
#' @param poly_ID_names list of poly_IDs to use
#' @param polygon_group_size number of polygons to process per group
#' @param return_gobject return giotto object (default: TRUE)
#' @param verbose be verbose
#' @return giotto object or spatVector with overlapping information
#' @details Serial overlapping function that works on groups of polygons at a time.
#'   Number of polygons per group is defined by \code{polygon_group_size} param
#' @concept overlap
#' @export
calculateOverlapSerial = function(gobject,
                                  name_overlap = NULL,
                                  spatial_info = 'cell',
                                  feat_info = 'rna',
                                  poly_ID_names = 'all',
                                  polygon_group_size = 500,
                                  return_gobject = TRUE,
                                  verbose = FALSE) {

  # spatial polygon
  spatvec = gobject@spatial_info[[spatial_info]]@spatVector

  # points polygon
  pointvec = gobject@feat_info[[feat_info]]@spatVector

  if(length(poly_ID_names) == 1) {
    if(poly_ID_names == 'all') {
      poly_ID_names = unique(spatvec$poly_ID)
    }
  }


  total_polygons = length(poly_ID_names)
  total_nr_groups = ceiling(total_polygons/polygon_group_size)
  groupnames = cut(1:total_polygons,
                   breaks = total_nr_groups,
                   labels = 1:total_nr_groups)
  names(poly_ID_names) = groupnames


  final_result = list()
  for(i in 1:total_nr_groups) {

    print((total_nr_groups-i))

    selected_poly_ID_names = poly_ID_names[names(poly_ID_names) == i]
    selected_spatvec = spatvec[spatvec$poly_ID %in% selected_poly_ID_names]

    # print(selected_spatvec)

    spatvec_result = overlap_points_per_polygon(spatvec = selected_spatvec,
                                                pointvec = pointvec,
                                                poly_ID_names = selected_poly_ID_names,
                                                verbose = verbose)

    final_result[[i]] = spatvec_result

  }

  final_result = do.call('rbind', final_result)


  if(return_gobject == TRUE) {

    if(is.null(name_overlap)) {
      name_overlap = feat_info
    }

    gobject@spatial_info[[spatial_info]]@overlaps[[name_overlap]] = final_result
    return(gobject)

  } else {
    return(final_result)
  }

}



#' @title Overlap points per polygon -- wrapped
#' @name overlap_points_per_polygon_wrapped
#' @description overlap wrapped polygons
#' @keywords internal
overlap_points_per_polygon_wrapped = function(spatvec_wrapped,
                                              pointvec_wrapped,
                                              poly_ID_names) {

  unwrap_spatvec = terra::vect(spatvec_wrapped)
  unwrap_pointvec = terra::vect(pointvec_wrapped)

  if(length(poly_ID_names) == 1) {
    if(poly_ID_names == 'all') {
      poly_ID_names = unique(unwrap_spatvec$poly_ID)
    }
  }

  intersect_res = overlap_points_per_polygon(spatvec = unwrap_spatvec,
                                             pointvec = unwrap_pointvec,
                                             poly_ID_names = poly_ID_names,
                                             verbose = FALSE)

  return(terra::wrap(intersect_res))

}



#' @title calculateOverlapParallel
#' @name calculateOverlapParallel
#' @description calculate overlap between cellular structures (polygons) and features (points)
#' @param gobject giotto object
#' @param name_overlap name for the overlap results (default to feat_info parameter)
#' @param spatial_info polygon information
#' @param feat_info feature information
#' @param poly_ID_names list of poly_IDs to use
#' @param polygon_group_size number of polygons to process per parallelization group
#' @param return_gobject return giotto object (default: TRUE)
#' @param verbose be verbose
#' @return giotto object or spatVector with overlapping information
#' @details parallel follows the future approach. This means that plan(multisession) does not work,
#' since the underlying terra objects are internal C pointers. plan(multicore) is also not supported for
#' Rstudio users.
#' @concept overlap
#' @export
calculateOverlapParallel = function(gobject,
                                    name_overlap = NULL,
                                    spatial_info = 'cell',
                                    feat_info = 'rna',
                                    poly_ID_names = 'all',
                                    polygon_group_size = 500,
                                    return_gobject = TRUE,
                                    verbose = TRUE) {

  # spatial polygon
  spatvec = gobject@spatial_info[[spatial_info]]@spatVector

  # points polygon
  pointvec = gobject@feat_info[[feat_info]]@spatVector


  if(length(poly_ID_names) == 1) {
    if(poly_ID_names == 'all') {
      poly_ID_names = unique(spatvec$poly_ID)
    }
  }

  total_polygons = length(poly_ID_names)
  total_nr_groups = ceiling(total_polygons/polygon_group_size)
  groupnames = cut(1:total_polygons,
                   breaks = total_nr_groups,
                   labels = 1:total_nr_groups)
  names(poly_ID_names) = groupnames

  # wrap SpatVector for points
  pointvec_wrap = terra::wrap(pointvec)

  # wrap SpatVectors for polygons
  spatvec_wrap_list = list()
  for(i in 1:total_nr_groups) {
    selected_poly_ID_names = poly_ID_names[names(poly_ID_names) == i]
    selected_spatvec = spatvec[spatvec$poly_ID %in% selected_poly_ID_names]
    spatvec_wrap_list[[i]] = terra::wrap(selected_spatvec)
  }


  # first intersect in parallel on wrapped terra objects
  result1 = lapply_flex(X = 1:length(spatvec_wrap_list),

                                 FUN = function(x) {
                                   test = overlap_points_per_polygon_wrapped(spatvec_wrapped = spatvec_wrap_list[[x]],
                                                                             pointvec_wrapped = pointvec_wrap,
                                                                             poly_ID_names = 'all')
                                 })

  # unwrap overlap results
  final_result = lapply(X = 1:length(result1), FUN = function(x) {
    terra::vect(result1[x][[1]])
  })

  # rbind all results together
  final_result = do.call('rbind', final_result)


  if(return_gobject == TRUE) {

    if(is.null(name_overlap)) {
      name_overlap = feat_info
    }

    gobject@spatial_info[[spatial_info]]@overlaps[[name_overlap]] = final_result
    return(gobject)

  } else {
    return(final_result)
  }

}







# * aggregate ####

#' @title overlapToMatrix
#' @name overlapToMatrix
#' @description create a count matrix based on overlap results from \code{\link{calculateOverlapRaster}}, \code{\link{calculateOverlapSerial}}, or \code{\link{calculateOverlapParallel}}
#' @param gobject giotto object
#' @param name name for the overlap count matrix
#' @param poly_info polygon information
#' @param feat_info feature information
#' @param return_gobject return giotto object (default: TRUE)
#' @return giotto object or count matrix
#' @concept overlap
#' @export
overlapToMatrix = function(gobject,
                           name = 'raw',
                           poly_info = 'cell',
                           feat_info = 'rna',
                           return_gobject = TRUE) {

  # define for data.table
  poly_ID = NULL

  overlap_spatvec = get_polygon_info(gobject = gobject,
                                     polygon_name = poly_info,
                                     polygon_overlap = feat_info)

  if(is.null(overlap_spatvec)) {
    cat('overlap between ', poly_info, ' and ', feat_info, ' has not been found \n')
    stop('Run calculateOverlap() first')
  }

  dtoverlap = spatVector_to_dt(overlap_spatvec)
  dtoverlap = dtoverlap[!is.na(poly_ID)] # removes points that have no overlap with any polygons
  #dtoverlap[, poly_ID := ifelse(is.na(poly_ID), 'no_overlap', poly_ID), by = 1:nrow(dtoverlap)]
  aggr_dtoverlap = dtoverlap[, .N, by = c('poly_ID', 'feat_ID')]


  # get all feature and cell information
  all_feats = gobject@feat_ID[[feat_info]]
  missing_feats = all_feats[!all_feats %in% unique(aggr_dtoverlap$feat_ID)]

  all_ids = gobject@cell_ID[[poly_info]]
  missing_ids = all_ids[!all_ids %in% unique(aggr_dtoverlap$poly_ID)]

  # create missing cell values, only if there are missing cell IDs!
  if(!length(missing_ids) == 0) {
    first_feature = aggr_dtoverlap[['feat_ID']][[1]]
    missing_dt = data.table::data.table(poly_ID = missing_ids, feat_ID = first_feature, N = 0)
    aggr_dtoverlap = rbind(aggr_dtoverlap, missing_dt)
  }

  if(!length(missing_feats) == 0) {
    first_cell = aggr_dtoverlap[['poly_ID']][[1]]
    missing_dt = data.table::data.table(poly_ID = first_cell, feat_ID = missing_feats, N = 0)
    aggr_dtoverlap = rbind(aggr_dtoverlap, missing_dt)
  }


  # TODO: creating missing feature values

  # create matrix
  overlapmatrixDT = data.table::dcast(data = aggr_dtoverlap,
                                      formula = feat_ID~poly_ID,
                                      value.var = 'N', fill = 0)
  overlapmatrix = dt_to_matrix(overlapmatrixDT)

  overlapmatrix = overlapmatrix[match(gobject@feat_ID[[feat_info]], rownames(overlapmatrix)),
                                match(gobject@cell_ID[[poly_info]], colnames(overlapmatrix))]

  overlapExprObj = create_expr_obj(name = name,
                                   exprMat = overlapmatrix,
                                   spat_unit = poly_info,
                                   feat_type = feat_info,
                                   provenance = poly_info)

  if(return_gobject == TRUE) {
    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
    gobject = set_expression_values(gobject, values = overlapExprObj)
    ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

    return(gobject)
  } else {
    return(overlapExprObj)
  }

}




#' @title overlapToMatrixMultiPoly
#' @name overlapToMatrixMultiPoly
#' @description create a count matrix based on overlap results from \code{\link{calculateOverlapRaster}}, \code{\link{calculateOverlapSerial}}, or \code{\link{calculateOverlapParallel}}
#' and aggregate information from multiple polygon layers (e.g. z-stacks) together
#' @param gobject giotto object
#' @param name name for the overlap count matrix
#' @param poly_info vector with polygon information
#' @param feat_info feature information
#' @param new_poly_info name for new aggregated polygon information
#' @param return_gobject return giotto object (default: TRUE)
#' @return giotto object or count matrix
#' @concept overlap
#' @export
overlapToMatrixMultiPoly = function(gobject,
                                    name = 'raw',
                                    poly_info = 'cell',
                                    feat_info = 'rna',
                                    new_poly_info = 'multi',
                                    return_gobject = TRUE) {


  # define for data.table
  i = j = x = NULL

  result_list = list()
  cell_ids_list = list()

  for(poly_info_i in 1:length(poly_info)) {

    poly_info_set = poly_info[[poly_info_i]]

    expr_names = list_expression_names(gobject = gobject,
                                       spat_unit = poly_info_set,
                                       feat_type = feat_info)

    # check if matrix already exists, if not try to make it
    if(!name %in% expr_names) {
      gobject = overlapToMatrix(gobject = gobject,
                                poly_info = poly_info_set,
                                feat_info = feat_info,
                                name = name)
    }

    testmat = get_expression_values(gobject = gobject,
                                    spat_unit = poly_info_set,
                                    feat_type = feat_info,
                                    values = name)

    featnames = dimnames(testmat)[[1]]
    names(featnames) = 1:length(featnames)

    colnames = dimnames(testmat)[[2]]
    names(colnames) = 1:length(colnames)

    testmat_DT = data.table::as.data.table(Matrix::summary(testmat))
    testmat_DT[, i := featnames[i]]
    testmat_DT[, j := colnames[j]]

    result_list[[poly_info_i]] = testmat_DT

    # cell ids
    #cell_ids = gobject@cell_ID[[poly_info_set]]

    #cell_ids_list[[poly_info_i]] = cell_ids

  }

  final_DT = data.table::rbindlist(result_list)
  final_DT_aggr = final_DT[, sum(x), by = .(i, j)]


  #combined_cell_IDs = sort(unique(unlist(cell_ids_list)))


  # create matrix
  overlapmatrixDT = data.table::dcast(data = final_DT_aggr,
                                      formula = i~j,
                                      value.var = 'V1', fill = 0)
  overlapmatrix = dt_to_matrix(overlapmatrixDT)

  #print(overlapmatrix[1:4, 1:4])


  #combined_cell_IDs = combined_cell_IDs[combined_cell_IDs %in% colnames(overlapmatrix)]

  #overlapmatrix = overlapmatrix[match(gobject@feat_ID[[feat_info]], rownames(overlapmatrix)),
  #                              match(combined_cell_IDs, colnames(overlapmatrix))]

  #print(overlapmatrix[1:4, 1:4])

  if(return_gobject == TRUE) {
    gobject@expression[[new_poly_info]][[feat_info]][[name]] = overlapmatrix
    gobject@cell_ID[[new_poly_info]] = colnames(overlapmatrix)

    gobject@cell_metadata[[new_poly_info]][[feat_info]] = data.table::data.table(cell_ID = colnames(overlapmatrix))
    gobject@feat_metadata[[new_poly_info]][[feat_info]] = data.table::data.table(feat_ID = rownames(overlapmatrix))

    return(gobject)
  } else {
    return(overlapmatrix)
  }

}





#' @title overlapImagesToMatrix
#' @name overlapImagesToMatrix
#' @description create a count matrix based on overlap results from \code{\link{calculateOverlapPolygonImages}}
#' @param gobject giotto object
#' @param name name for the overlap count matrix
#' @param poly_info polygon information
#' @param feat_info feature information
#' @param image_names names of images you used
#' @param spat_locs_name name for spatial centroids / locations associated with matrix
#' @param return_gobject return giotto object (default: TRUE)
#' @return giotto object or data.table with aggregated information
#' @concept overlap
#' @export
overlapImagesToMatrix = function(gobject,
                                 name = 'raw',
                                 poly_info = 'cell',
                                 feat_info = 'protein',
                                 image_names = NULL,
                                 spat_locs_name = 'raw',
                                 return_gobject = TRUE) {

  ## get polygon information
  polygon_info = get_polygon_info(gobject = gobject,
                                  polygon_name = poly_info,
                                  return_giottoPolygon = T)


  poly_info_image_overlap_names = names(polygon_info@overlaps$images)

  overlap_aggr_list = list()

  for(img_overlap in poly_info_image_overlap_names) {
    overlap_DT = polygon_info@overlaps[['images']][[img_overlap]]
    aggr_overlap_DT = overlap_DT[, mean(get(img_overlap)), by = 'ID']
    aggr_overlap_DT[, feat := img_overlap]
    colnames(aggr_overlap_DT) = c('poly_ID', 'mean_intensity', 'feat_ID')
    overlap_aggr_list[[img_overlap]] = aggr_overlap_DT
  }

  aggr_comb = do.call('rbind', overlap_aggr_list)



  if(return_gobject) {

    cell_IDs = unique(aggr_comb$poly_ID)
    feat_IDs = unique(aggr_comb$feat_ID)

    # create cell and feature metadata
    gobject@cell_metadata[[poly_info]][[feat_info]] = data.table::data.table(cell_ID = cell_IDs)
    gobject@feat_metadata[[poly_info]][[feat_info]] = data.table::data.table(feat_ID = feat_IDs)

    # add feat_ID and cell_ID
    gobject@feat_ID[[feat_info]] = feat_IDs
    gobject@cell_ID[[poly_info]] = cell_IDs

    # add spatial locations
    centroidsDT = gobject@spatial_info[[poly_info]]@spatVectorCentroids
    centroidsDT = spatVector_to_dt(centroidsDT)
    centroidsDT_loc = centroidsDT[, .(poly_ID, x, y)]
    colnames(centroidsDT_loc) = c('cell_ID', 'sdimx', 'sdimy')

    gobject = set_spatial_locations(gobject = gobject,
                                    spat_unit = poly_info,
                                    spat_loc_name = ,
                                    spatlocs = centroidsDT_loc,
                                    verbose = FALSE)

    # create matrix
    overlapmatrixDT = data.table::dcast(data = aggr_comb,
                                        formula = feat_ID~poly_ID,
                                        value.var = 'mean_intensity', fill = 0)
    overlapmatrix = dt_to_matrix(overlapmatrixDT)

    overlapmatrix = overlapmatrix[match(gobject@feat_ID[[feat_info]], rownames(overlapmatrix)),
                                  match(gobject@cell_ID[[poly_info]], colnames(overlapmatrix))]

    gobject = set_expression_values(gobject = gobject,
                                    spat_unit = poly_info,
                                    feat_type = feat_info,
                                    name = name,
                                    values = overlapmatrix)

  } else {
    return(aggr_comb)
  }

}


# * combine metadata ####

#' @title combineCellData
#' @name combineCellData
#' @description combine cell data information
#' @param gobject giotto object
#' @param feat_type feature type
#' @param include_spat_locs include information about spatial locations
#' @param spat_loc_name spatial location name
#' @param include_poly_info include information about polygon
#' @param poly_info polygon information name
#' @return data.table with combined spatial information
#' @concept combine cell metadata
#' @export
combineCellData = function(gobject,
                           feat_type = 'rna',
                           include_spat_locs = TRUE,
                           spat_loc_name = 'raw',
                           include_poly_info = TRUE,
                           poly_info = 'cell') {

  # combine
  # 1. spatial morphology information ( = polygon)
  # 2. cell metadata

  # specify feat_type
  # Set feat_type and spat_unit
  poly_info = set_default_spat_unit(gobject = gobject,
                                    spat_unit = poly_info)
  feat_type = set_default_feat_type(gobject = gobject,
                                    spat_unit = poly_info,
                                    feat_type = feat_type)


  # get spatial locations
  if(include_spat_locs == TRUE) {
    spat_locs_dt = get_spatial_locations(gobject = gobject,
                                         spat_unit = poly_info,
                                         spat_loc_name = spat_loc_name,
                                         output = 'data.table',
                                         copy_obj = TRUE)
  } else {
    spat_locs_dt = NULL
  }


  # get spatial cell information
  if(include_poly_info == TRUE) {
    # get spatial cell information
    spatial_cell_info_spatvec = get_polygon_info(gobject = gobject,
                                                 polygon_name = poly_info,
                                                 return_giottoPolygon = FALSE)
    spatial_cell_info_dt = spatVector_to_dt(spatial_cell_info_spatvec,
                                            include_values = TRUE)
    data.table::setnames(spatial_cell_info_dt, old = 'poly_ID', new = 'cell_ID')
  } else {
    spatial_cell_info_dt = NULL
  }




  # combine prior information if wanted
  if(!is.null(spat_locs_dt) & !is.null(spatial_cell_info_dt)) {
    comb_dt = data.table::merge.data.table(spat_locs_dt,
                                           spatial_cell_info_dt,
                                           by = 'cell_ID')
  } else if(!is.null(spat_locs_dt)) {
    comb_dt = spat_locs_dt
  } else if(!is.null(spatial_cell_info_dt)) {
    comb_dt = spatial_cell_info_dt
  } else {
    comb_dt = NULL
  }


  res_list = list()
  for(feat in unique(feat_type)) {


    # get spatial cell metadata
    cell_meta = get_cell_metadata(gobject = gobject,
                                  spat_unit = poly_info,
                                  feat_type = feat,
                                  output = 'data.table')

    # merge
    if(!is.null(comb_dt)) {
      spatial_cell_info_meta = data.table::merge.data.table(comb_dt, cell_meta, by = 'cell_ID')
    } else {
      spatial_cell_info_meta = cell_meta
    }

    spatial_cell_info_meta[, 'feat' := feat]

    res_list[[feat]] = spatial_cell_info_meta

  }

  return(res_list)



}


#' @title combineFeatureData
#' @name combineFeatureData
#' @description combine feature data information
#' @param gobject giotto object
#' @param feat_type feature type
#' @param spat_unit spatial unit
#' @param sel_feats selected features (default: NULL or no selection)
#' @return data.table with combined spatial feature information
#' @concept combine feature metadata
#' @export
combineFeatureData = function(gobject,
                              feat_type = NULL,
                              spat_unit = NULL,
                              sel_feats = NULL) {

  # define for data.table [] subsetting
  feat_ID = NULL

  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = spat_unit)
  feat_type = set_default_feat_type(gobject = gobject,
                                    spat_unit = spat_unit,
                                    feat_type = feat_type)

  res_list = list()
  for(feat in unique(feat_type)) {
    for(spat in unique(spat_unit)) {

      # feature meta
      # feat_meta = gobject@feat_metadata[[spat_unit]][[feat]]

      feat_meta = get_feature_metadata(gobject = gobject,
                                       spat_unit = spat_unit,
                                       feat_type = feat,
                                       output = 'data.table')

      if(!is.null(sel_feats[[feat_type]])) {
        selected_features = sel_feats[[feat_type]]
        feat_meta = feat_meta[feat_ID %in% selected_features]
      }


      # feature info
      feat_info_spatvec = get_feature_info(gobject = gobject,
                                           feat_type = feat)
      feat_info = spatVector_to_dt(feat_info_spatvec)
      if(!is.null(sel_feats[[feat_type]])) {
        selected_features = sel_feats[[feat_type]]
        feat_info = feat_info[feat_ID %in% selected_features]
      }

      comb_dt = data.table::merge.data.table(x = feat_meta,
                                             y = feat_info,
                                             by = 'feat_ID')

      comb_dt[, 'feat' := feat]
      comb_dt[, 'spat_unit' := spat]

    }

    res_list[[feat]] = comb_dt

  }

  return(res_list)

}


#' @title combineFeatureOverlapData
#' @name combineFeatureOverlapData
#' @description combine feature data information
#' @param gobject giotto object
#' @param feat_type feature type
#' @param sel_feats selected features (default: NULL or no selection)
#' @param poly_info polygon information name
#' @return data.table with combined spatial polygon information
#' @concept combine feature metadata
#' @export
combineFeatureOverlapData = function(gobject,
                                     feat_type = 'rna',
                                     sel_feats = NULL,
                                     poly_info = c('cell')) {

  # define for data.table [] subsetting
  feat_ID = NULL

  poly_info = set_default_spat_unit(gobject = gobject,
                                    spat_unit = poly_info)
  feat_type = set_default_feat_type(gobject = gobject,
                                    spat_unit = poly_info,
                                    feat_type = feat_type)


  res_list = list()
  for(feat in unique(feat_type)) {

    for(spat in unique(poly_info)) {

      # feature meta
      # feat_meta = gobject@feat_metadata[[feat]][[spat]]

      feat_meta = get_feature_metadata(gobject = gobject,
                                       spat_unit = spat,
                                       feat_type = feat,
                                       output = 'data.table')

      print(feat_meta)

      if(!is.null(sel_feats[[feat_type]])) {
        selected_features = sel_feats[[feat_type]]
        feat_meta = feat_meta[feat_ID %in% selected_features]
      }

      # overlap poly and feat info
      poly_list = list()
      for(poly in poly_info) {
        feat_overlap_info_spatvec = get_polygon_info(gobject = gobject,
                                                     polygon_name = poly,
                                                     polygon_overlap = feat)
        feat_overlap_info = spatVector_to_dt(feat_overlap_info_spatvec)

        if(!is.null(sel_feats[[feat_type]])) {
          selected_features = sel_feats[[feat_type]]
          feat_overlap_info = feat_overlap_info[feat_ID %in% selected_features]
        }

        feat_overlap_info[, poly_info := poly]
        poly_list[[poly]] = feat_overlap_info
      }

      poly_list_res = do.call('rbind', poly_list)

      comb_dt = data.table::merge.data.table(x = feat_meta,
                                             y = poly_list_res,
                                             by = 'feat_ID')

    }

    print(comb_dt)

    comb_dt[, 'feat' := feat]
    res_list[[feat]] = comb_dt

  }

  return(res_list)

}




