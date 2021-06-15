


## ** cell shape polygons ####

## Giotto polygon class

#' @title S4 giotto polygon Class
#' @description Giotto class to store and operate on polygon-like data
#' @keywords giotto, polygon, class
#' @slot name name of polygon shapes
#' @slot spatVector terra spatVector to store polygon shapes
#' @slot spatVectorCentroids centroids of polygon shapes
#' @slot overlaps information about overlapping points and polygons
#' @details
#'
#' @export
giottoPolygon <- setClass(
  Class = "giottoPolygon",

  slots = c(
    name = "ANY",
    spatVector = "ANY",
    spatVectorCentroids = "ANY",
    overlaps = "ANY"
  ),

  prototype = list(
    name = NULL,
    spatVector = NULL,
    spatVectorCentroids = NULL,
    overlaps = NULL
  )
)


## extension of spatVector object
## name should match the cellular structure

#' @name create_giotto_polygon_object
#' @keywords internal
create_giotto_polygon_object = function(name = 'cell',
                                        spatVector = NULL,
                                        spatVectorCentroids = NULL) {


  # create minimum giotto
  g_polygon = giottoPolygon(name = name,
                            spatVector = NULL,
                            spatVectorCentroids = NULL)


  ## 1. check magick image object
  if(!methods::is(spatVector, 'SpatVector')) {
    stop("spatVector needs to be a spatVector object from the terra package")
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

  # provide name
  g_polygon@name = name

  # image object
  return(g_polygon)
}




#' @name identify_background_range_polygons
#' @description function to remove background polygon based on largest range
#' @keywords internal
identify_background_range_polygons = function(spatVector) {

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




#' @name create_segm_polygons
#' @description creates giotto polygons from segmentation mask data
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

  terra_rast = terra::rast(maskfile)
  rast_dimensions = dim(terra_rast)

  terra_polygon = terra::as.polygons(x = terra_rast, value = TRUE)
  names(terra_polygon) = 'mask'


  ## flip axes ##
  if(flip_vertical == TRUE) {
    terra_polygon = flip(terra_polygon, direction = 'vertical')
  }

  if(flip_horizontal == TRUE) {
    terra_polygon = flip(terra_polygon, direction = 'horizontal')
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




#' @name fix_multipart_geoms
#' @description function to split geoms (polygons) that have multiple parts
#' @keywords internal
fix_multipart_geoms = function(spatVector) {

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

#' @name createGiottoPolygonsFromMask
#' @description Creates Giotto polygon object from a mask file (e.g. segmentation results)
#' @param maskfile path to mask file
#' @param name name for polygons
#' @param remove_background_polygon try to remove background polygon (default: FALSE)
#' @param background_algo algorithm to remove background polygon
#' @param poly_IDs unique nanes for each polygon in the mask file
#' @param flip_vertical flip mask figure in a vertical manner
#' @param shift_vertical_step shift vertical (boolean or numerical)
#' @param flip_horizontal flip mask figure in a horizontal manner
#' @param shift_horizontal_step shift horizontal (boolean or numerical)
#' @param calc_centroids calculate centroids for polygons
#' @param fix_multipart try to split polygons with multiple parts (default: TRUE)
#' @return
#' @keywords mask polygon
#' @export
createGiottoPolygonsFromMask = function(maskfile,
                                        mask_method = c('guess', 'single', 'multiple'),
                                        name = 'cell',
                                        remove_background_polygon = FALSE,
                                        background_algo = c('range'),
                                        poly_IDs = NULL,
                                        flip_vertical = TRUE,
                                        shift_vertical_step = TRUE,
                                        flip_horizontal = TRUE,
                                        shift_horizontal_step = TRUE,
                                        calc_centroids = FALSE,
                                        fix_multipart = TRUE) {

  # select background algo
  background_algo = match.arg(background_algo, choices = 'range')

  # check if mask file exists
  if(!file.exists(maskfile)) {
    stop('path : ', maskfile, ' does not exist \n')
  }

  # mask method
  # single: single mask value for all segmented cells
  # multiple: multiple maks values and thus a unique value for each segmented cell
  mask_method = match.arg(mask_method, choices = c('guess', 'single', 'multiple'))

  # create polygons from mask
  terra_rast = terra::rast(maskfile)
  rast_dimensions = dim(terra_rast)
  terra_polygon = terra::as.polygons(x = terra_rast, value = TRUE)
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

    if(background_algo == 'range') {
      backgr_poly_id = identify_background_range_polygons(terra_polygon)
    }

    terra_polygon = terra::subset(x = terra_polygon, terra_polygon[['polyID']] != backgr_poly_id)

  }


  # provide own cell_ID name
  if(!is.null(poly_IDs)) {
    if(length(poly_IDs) != nrow(terra::values(terra_polygon))) {
      stop('length cell_IDs does not equal number of found polyongs \n')
    }
    terra_polygon$poly_ID = poly_IDs
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

#' @name createGiottoPolygonsFromDfr
#' @description Creates Giotto polygon object from a structured dataframe-like object
#' @param segmdfr data.frame-like object with polygon coordinate information (x, y, ID)
#' @param name name for polygons
#' @param calc_centroids calculate centroids for polygons
#' @return
#' @keywords polygon
#' @export
createGiottoPolygonsFromDfr = function(segmdfr,
                                       name = 'cell',
                                       calc_centroids = FALSE) {

  # input data.frame-like object
  # columns: x y cell_ID
  segmdt = data.table::as.data.table(segmdfr)
  input_dt = segmdt[,c(1:3), with = F]
  colnames(input_dt) = c('x', 'y', 'poly_ID')


  nr_of_cells_vec = 1:length(unique(input_dt$poly_ID))
  names(nr_of_cells_vec) = unique(input_dt$poly_ID)

  # add other colnaes for the input data.table
  input_dt[, geom := nr_of_cells_vec[[get('poly_ID')]] , by = 1:nrow(input_dt)]
  #input_dt[, mask := nr_of_cells_vec[[get('poly_ID')]], by = 1:nrow(input_dt)]
  input_dt[, c('part', 'hole') := list(1, 0)]
  input_dt = input_dt[, c('geom', 'part', 'x', 'y', 'hole', 'poly_ID'), with =F]

  # create spatvector
  spatvector = dt_to_spatVector_polygon(input_dt,
                                        include_values = TRUE)


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




#' @name fix_multipart_geoms
#' @description  to extract list of polygons
#' @keywords internal
extract_polygon_list = function(polygonlist) {

  # if polygonlist is not a names list
  # try to make list and give default names
  if(!is.list(polygonlist)) {
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

    if(length(polygonlist) == 1) {
      names(polygonlist) = 'cell'
    } else {
      polygonlist_l = length(polygonlist)
      names(polygonlist) = c('cell', paste0('info', 1:(polygonlist_l-1)))

    }
  }

  # make sure cell is one of the names
  all_names = names(polygonlist)
  if(!any('cell' %in% all_names)) {
    stop(" Information about 'cell' needs to be provided")
  }


  final_list = list()

  for(poly_i in 1:length(polygonlist)) {

    name_polyinfo = names(polygonlist)[[poly_i]]
    polyinfo = polygonlist[[poly_i]]

    if(is.character(polyinfo)) {

      poly_results = createGiottoPolygonsFromMask(name = name_polyinfo,
                                                  maskfile = polyinfo)

    } else if(inherits(polyinfo, 'data.frame')) {

      poly_results = createGiottoPolygonsFromDfr(name = 'cell',
                                                 segmdfr = polyinfo)

    } else if(inherits(polyinfo, 'giottoPolygon')) {

      poly_results = polyinfo

    } else {

      stop('Polygon can only be extraxted from a mask file or from a correctly formatted data.frame')

    }

    final_list[[name_polyinfo]] = poly_results

  }

  return(final_list)

}




#' @name addGiottoPolygons
#' @description Adds Giotto polygon to an existing Giotto object
#' @param gobject giotto object
#' @param gpolygons list of giotto polygon objects
#' @return
#' @keywords polygon
#' @export
addGiottoPolygons = function(gobject,
                             gpolygons) {

  cat('Does not exist yet')

}




# helper functions

# convert spatVector to data.table

#' @name spatVector_to_dt
#' @description  convert spatVector to data.table
#' @keywords internal
spatVector_to_dt = function(spatvector,
                            include_values = TRUE) {

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


    spatVec = terra::vect(x = as.matrix(dt[,geom_values, with = F]),
                          type = 'polygons', atts = unique(dt[,other_values, with = F]))

  } else {

    spatVec = terra::vect(x = as.matrix(dt[,geom_values, with = F]),
                          type = 'polygons', atts = NULL)

  }

  return(spatVec)

}


#' @name spline_poly
#' @description spline polynomial to smooth polygon
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
  data.spline <- spline(1:(n+2*k), data[,1], n=vertices, ...)
  x <- data.spline$x
  x1 <- data.spline$y
  x2 <- spline(1:(n+2*k), data[,2], n=vertices, ...)$y

  # Retain only the middle part.
  cbind(x1, x2)[k < x & x <= n+k, ]
}



#' @name smoothGiottoPolygons
#' @description Smooths Giotto polygon object
#' @param gpolygon giotto polygon object
#' @param vertices number of vertices
#' @param k k
#' @param set_neg_to_zero set negative values to zero (default: TRUE)
#' @return
#' @keywords polygon
#' @export
smoothGiottoPolygons = function(gpolygon,
                                vertices = 20,
                                k = 3,
                                set_neg_to_zero = TRUE,
                                ...) {

  polygDT = spatVector_to_dt(gpolygon@spatVector)

  # store other values
  all_colnames = colnames(polygDT)
  geom_values = c('geom', 'part', 'x', 'y', 'hole')
  other_values = all_colnames[!all_colnames %in% geom_values]
  other_values_uniq_dt = unique(polygDT[,c('geom', 'part', 'hole', other_values), with =F])

  # apply smoothing to each polygon
  comb = lapply(1:length(unique(polygDT$geom)), FUN = function(z) {

    polygMat = as.matrix(polygDT[geom == z,.(x, y)])

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
      plot(new_spatvec[new_spatvec$poly_ID == ID])
      #orig_spatvector = gpolygon@spatVector
      #new_spatvec[new_spatvec$poly_ID == ID] = orig_spatvector[orig_spatvector$poly_ID == ID]
    }
  }

  new_gpolygon = create_giotto_polygon_object(name = gpolygon@name,
                                              spatVector = new_spatvec,
                                              spatVectorCentroids = gpolygon@spatVectorCentroids)

  return(new_gpolygon)

}






## ** feature points ####

# giotto class for points


#' @title S4 giotto points Class
#' @description Giotto class to store and operate on points data
#' @keywords giotto, points, class
#' @slot feat_type name of feature type
#' @slot spatVector terra spatVector to store point shapes
#' @details
#'
#' @export
giottoPoints <- setClass(
  Class = "giottoPoints",

  slots = c(
    feat_type = "ANY",
    spatVector = "ANY"
  ),

  prototype = list(
    feat_type = NULL,
    spatVector = NULL
  )
)


#' @name create_giotto_points_object
#' @keywords internal
create_giotto_points_object = function(feat_type = 'rna',
                                       spatVector = NULL) {


  # create minimum giotto points object
  g_points = giottoPoints(feat_type = feat_type,
                          spatVector = NULL)

  ## 1. check terra spatVector object
  if(!methods::is(spatVector, 'SpatVector')) {
    stop("spatVector needs to be a spatVector object from the terra package")
  }

  g_points@spatVector = spatVector

  ## 2. provide feature id
  g_points@feat_type = feat_type

  # image object
  return(g_points)

}



#' @name create_giotto_points_object
#' @description create terra spatvector from a data.frame
#' @keywords internal
create_spatvector_object_from_dfr = function(x) {

  # data.frame like object needs to have 2 coordinate columns and
  # at least one other column as the feat_ID
  loc_dfr = data.table::as.data.table(x[, 1:2])
  att_dfr = data.table::as.data.table(x[, -c(1:2)])

  spatvec = terra::vect(as.matrix(x[,1:2]), type = 'points', atts = att_dfr)
  names(spatvec)[1] = 'feat_ID'

  # will be given and is a unique numerical barcode for each feature
  spatvec[['feat_ID_uniq']] = 1:nrow(spatvec)

  return(spatvec)

}


# create Giotto points from data.frame or spatVector


#' @name createGiottoPoints
#' @description Creates Giotto point object from a structured dataframe-like object
#' @param x data.frame-like object with points coordinate information (x, y, feat ID)
#' @param feat_type feature type
#' @return
#' @keywords polygon
#' @export
createGiottoPoints = function(x,
                              feat_type = 'rna') {

  if(inherits(x, 'data.frame')) {

    spatvec = create_spatvector_object_from_dfr(x = x)
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


#' @name dt_to_spatVector_points
#' @description  data.table to spatVector for points
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


#' @name addGiottoPoints
#' @description Adds Giotto points to an existing Giotto object
#' @param gobject giotto object
#' @param gpoints list of giotto point objects
#' @return
#' @keywords polygon
#' @export
addGiottoPoints = function(gobject,
                           gpoints) {

  cat('Does not exist yet')

}





#' @name extract_points_list
#' @description  to extract list of giotto points
#' @keywords internal
extract_points_list = function(pointslist) {

  # if polygonlist is not a names list
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


    point_results = createGiottoPoints(x = pointinfo,
                                       feat_type = name_pointinfo)

    final_list[[name_pointinfo]] = point_results

  }

  return(final_list)


}






## ** giotto structure functions ####



#' @name addSpatialCentroidLocations
#' @description Calculates the centroid locations for the giotto polygons
#' @param gobject giotto object
#' @param poly_info polygon information
#' @param return_gobject return giotto object (default: TRUE)
#' @return
#' @keywords centroid
#' @export
addSpatialCentroidLocations = function(gobject,
                                       poly_info = 'cell',
                                       return_gobject = TRUE) {


  extended_spatvector = calculate_centroids_polygons(gpolygon = gobject@spatial_info[[poly_info]],
                                                     name = 'centroids',
                                                     append_gpolygon = TRUE)

  centroid_spatvector = spatVector_to_dt(extended_spatvector@spatVectorCentroids)

  spatial_locs = centroid_spatvector[, .(x,y,poly_ID)]
  colnames(spatial_locs) = c('sdimx', 'sdimy', 'cell_ID')

  if(return_gobject == TRUE) {

    # add spatial locations
    gobject@spatial_locs[[poly_info]] = spatial_locs

    # add centroids information
    gobject@spatial_info[[poly_info]] = extended_spatvector

    return(gobject)

  } else {
    return(spatial_locs)
  }

}





## calculate overlap between cellular structures and features

#' @name calculateOverlap
#' @description calculate overlap between cellular structures (polygons) and features (points)
#' @param gobject giotto object
#' @param name_overlap name for the overlap results (default to feat_info parameter)
#' @param poly_info polygon information
#' @param feat_info feature information
#' @param x_step x-direction step to travel over the polygon landscape
#' @param y_step y-direction step to travel over the polygon landscape
#' @param return_gobject return giotto object (default: TRUE)
#' @return giotto object or spatVector with overlapping information
#' @keywords overlap
#' @export
calculateOverlap = function(gobject,
                            poly_info = 'cell',
                            feat_info = 'rna',
                            x_step = 200,
                            y_step = 200,
                            return_gobject = TRUE,
                            name_overlap = NULL) {


  polvec = gobject@spatial_info[[poly_info]]@spatVector
  pointsvec = gobject@feat_info[[feat_info]]@spatVector

  ## compute windows to look for overlap
  ## create data.table with x and y beginnings and ends
  myext = terra::ext(polvec)
  range_x = xmax(myext) - xmin(myext)
  range_y = ymax(myext) - ymin(myext)

  xrep = ceiling(range_x / x_step)
  yrep = ceiling(range_y / y_step)

  start_x = xmin(myext)
  end_x = start_x + (xrep * x_step)

  start_y = ymin(myext)
  end_y = start_y + (yrep * y_step)

  dt_steps = data.table::data.table(xmin = rep(seq(start_x, end_x, x_step), yrep),
                                    xmax = rep(seq(x_step, (end_x+x_step), x_step), yrep),
                                    ymin = rep(seq(start_y, end_y, y_step), each = xrep),
                                    ymax = rep(seq(y_step, (end_y+y_step), y_step), each = xrep))

  print(dt_steps)


  all_cell_IDs = unique(polvec$poly_ID)
  remain_cell_IDs = all_cell_IDs

  resultsteplist = list()

  i = 1
  for(row in 1:nrow(dt_steps)) {

    print(row)

    crop_polvec = terra::crop(x = polvec,
                              y = terra::ext(dt_steps[row]$xmin, dt_steps[row]$xmax,
                                      dt_steps[row]$ymin, dt_steps[row]$ymax))


    ## only continue if crop_polvec overlaps at least one polygon
    if(length(crop_polvec) > 0) {
      select_cell_IDs = unique(crop_polvec$poly_ID)
      select_cell_IDs = select_cell_IDs[select_cell_IDs %in% remain_cell_IDs]


      if(length(select_cell_IDs) > 0) {

        ## create subset based on selected cell IDs
        subpolvec = terra::subset(polvec, polvec$poly_ID %in% select_cell_IDs)
        subpolvecDT = spatVector_to_dt(subpolvec)

        ## subset points based on range of selected cell IDs polygons
        range_x = range(subpolvecDT$x)
        range_y = range(subpolvecDT$y)

        pointsvecDT = spatVector_to_dt(pointsvec)
        bool_filter = pointsvecDT$x > range_x[1] & pointsvecDT$x < range_x[2] & pointsvecDT$y > range_y[1] & pointsvecDT$y < range_y[2]
        subpointsvec = pointsvec[bool_filter]


        if(length(subpointsvec) > 0) {

          cat('subpolvec ', nrow(subpolvec), '\n')
          cat('subpointsvec ', nrow(subpointsvec), '\n')

          subtestsect = terra::intersect(x = subpolvec, y = subpointsvec)
          resultsteplist[[i]] = subtestsect
          i = i + 1
        }

      }

      remain_cell_IDs = remain_cell_IDs[!remain_cell_IDs %in% select_cell_IDs]

    }


  }

  final_result = do.call('c', resultsteplist)

  if(return_gobject == TRUE) {

    if(is.null(name_overlap)) {
      name_overlap = feat_info
    }

    gobject@spatial_info[[poly_info]]@overlaps[[name_overlap]] = final_result
    return(gobject)

  } else {
    return(final_result)
  }



}





#' @name overlapToMatrix
#' @description create a count matrix based on overlap results from \code{\link{calculateOverlap}}
#' @param gobject giotto object
#' @param name name for the overlap count matrix
#' @param poly_info polygon information
#' @param feat_info feature information
#' @param return_gobject return giotto object (default: TRUE)
#' @return giotto object or count matrix
#' @keywords overlap
#' @export
overlapToMatrix = function(gobject,
                           name = 'raw',
                           poly_info = 'cell',
                           feat_info = 'rna',
                           return_gobject = TRUE) {


  overlap_spatvec = select_polygon_info(gobject = gobject,
                                        polygon_name = poly_info,
                                        polygon_overlap = feat_info)

  if(is.null(overlap_spatvec)) {
    cat('overlap between ', poly_info, ' and ', feat_info, ' has not been found \n')
    stop('Run calculateOverlap() first')
  }

  dtoverlap = spatVector_to_dt(overlap_spatvec)
  aggr_dtoverlap = dtoverlap[, .N, by = c('poly_ID', 'feat_ID')]


  # get all feature and cell information
  all_feats = gobject@feat_ID[[feat_info]]
  missing_feats = all_feats[!all_feats %in% unique(aggr_dtoverlap$feat_ID)]

  all_ids = gobject@cell_ID
  missing_ids = all_ids[!all_ids %in% unique(aggr_dtoverlap$poly_ID)]

  # create missing cell values
  first_feature = aggr_dtoverlap[['feat_ID']][[1]]
  missing_dt = data.table::data.table(poly_ID = missing_ids, feat_ID = first_feature, N = 0)
  aggr_dtoverlap = rbind(aggr_dtoverlap, missing_dt)

  # TODO: creating missing feature values

  # create matrix
  overlapmatrixDT = data.table::dcast(data = aggr_dtoverlap,
                                      formula = feat_ID~poly_ID,
                                      value.var = 'N', fill = 0)
  overlapmatrix = Giotto:::dt_to_matrix(overlapmatrixDT)

  overlapmatrix = overlapmatrix[match(gobject@feat_ID[[feat_info]], rownames(overlapmatrix)),
                                match(gobject@cell_ID, colnames(overlapmatrix))]


  if(return_gobject == TRUE) {
    gobject@expression[[feat_info]][[name]] = overlapmatrix
    return(gobject)
  } else {
    return(overlapmatrix)
  }

}






#' @name combineCellData
#' @description combine cell data information
#' @param gobject giotto object
#' @param feat_type feature type
#' @param include_spat_locs include information about spatial locations
#' @param spat_loc_name spatial location name
#' @param include_poly_info include information about polygon
#' @param poly_info polygon information name
#' @return data.table with combined spatial information
#' @keywords combine cell metadata
#' @export
combineCellData = function(gobject,
                           feat_type = 'rna',
                           include_spat_locs = TRUE,
                           spat_loc_name = 'cell',
                           include_poly_info = TRUE,
                           poly_info = 'cell') {

  # combine
  # 1. spatial morphology information ( = polygon)
  # 2. cell metadata

  # specify feat_type
  if(is.null(feat_type)) {
    feat_type = gobject@expression_feat[[1]]
  }

  # get spatial locations
  if(include_spat_locs == TRUE) {
    spat_locs_dt = select_spatial_locations(gobject = gobject,
                                            spat_loc_name = spat_loc_name)
  } else {
    spat_locs_dt = NULL
  }


  # get spatial cell information
  if(include_poly_info == TRUE) {
    # get spatial cell information
    spatial_cell_info_spatvec = select_polygon_info(gobject = gobject,
                                                    polygon_name = poly_info)
    spatial_cell_info_dt = spatVector_to_dt(spatial_cell_info_spatvec,
                                            include_values = TRUE)
    data.table::setnames(spatial_cell_info_dt, old = 'poly_ID', new = 'cell_ID')
  } else {
    spatial_cell_info_dt = NULL
  }


  # combine prior information if wanted
  if(!is.null(spat_locs_dt) & !is.null(spatial_cell_info_dt)) {
    comb_dt = data.table::merge.data.table(spat_locs_dt, spatial_cell_info_dt, by = 'cell_ID')
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
    cell_meta = pDataDT(gobject, feat_type = feat)

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


#' @name combineFeatureData
#' @description combine feature data information
#' @param gobject giotto object
#' @param feat_type feature type
#' @param sel_feats selected features (default: NULL or no selection)
#' @return data.table with combined spatial feature information
#' @keywords combine feature metadata
#' @export
combineFeatureData = function(gobject,
                              feat_type = 'rna',
                              sel_feats = NULL) {


  # specify feat_type
  if(is.null(feat_type)) {
    feat_type = gobject@expression_feat[[1]]
  }


  res_list = list()
  for(feat in unique(feat_type)) {

    # feature meta
    feat_meta = gobject@feat_metadata[[feat]]
    if(!is.null(sel_feats[[feat_type]])) {
      selected_features = sel_feats[[feat_type]]
      feat_meta = feat_meta[feat_ID %in% selected_features]
    }


    # feature info
    feat_info_spatvec = select_feature_info(gobject = gobject,
                                            feat_name = feat)
    feat_info = spatVector_to_dt(feat_info_spatvec)
    if(!is.null(sel_feats[[feat_type]])) {
      selected_features = sel_feats[[feat_type]]
      feat_info = feat_info[feat_ID %in% selected_features]
    }

    comb_dt = data.table::merge.data.table(x = feat_meta,
                                           y = feat_info,
                                           by = 'feat_ID')

    comb_dt[, 'feat' := feat]

    res_list[[feat]] = comb_dt

  }

  return(res_list)

}



#' @name combineFeatureOverlapData
#' @description combine feature data information
#' @param gobject giotto object
#' @param feat_type feature type
#' @param sel_feats selected features (default: NULL or no selection)
#' @param poly_info polygon information name
#' @return data.table with combined spatial polygon information
#' @keywords combine feature metadata
#' @export
combineFeatureOverlapData = function(gobject,
                                     feat_type = 'rna',
                                     sel_feats = NULL,
                                     poly_info = c('cell')) {


  # specify feat_type
  if(is.null(feat_type)) {
    feat_type = gobject@expression_feat[[1]]
  }


  res_list = list()
  for(feat in unique(feat_type)) {

    # feature meta
    feat_meta = gobject@feat_metadata[[feat]]
    if(!is.null(sel_feats[[feat_type]])) {
      selected_features = sel_feats[[feat_type]]
      feat_meta = feat_meta[feat_ID %in% selected_features]
    }

    # overlap poly and feat info
    poly_list = list()
    for(poly in poly_info) {
      feat_overlap_info_spatvec = select_polygon_info(gobject = gobject,
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
    comb_dt[, 'feat' := feat]

    res_list[[feat]] = comb_dt

  }

  return(res_list)

}




#' @name showGiottoSpatialInfo
#' @description show the available giotto spatial polygon information
#' @param gobject giotto object
#' @keywords show
#' @export
showGiottoSpatialInfo = function(gobject) {

  for(info in names(gobject@spatial_info)) {

    cat("For Spatial info: ", info, "\n\n")
    print(gobject@spatial_info[[info]])
    cat("-----------------------------")
    cat("\n \n")
  }

}


#' @name showGiottoFeatInfo
#' @description show the available giotto spatial feature information
#' @param gobject giotto object
#' @keywords show
#' @export
showGiottoFeatInfo = function(gobject) {

  for(info in names(gobject@feat_info)) {

    cat("For Feature info: ", info, "\n\n")
    print(gobject@feat_info[[info]])
    cat("-----------------------------")
    cat("\n \n")
  }

}

