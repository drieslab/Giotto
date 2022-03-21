

# giottoImage class ####

#' @title S4 giottoImage Class
#' @description Framework of giotto object to store and work with spatial expression data
#' @keywords giotto, object
#' @slot name name of Giotto image
#' @slot mg_object magick image object
#' @slot minmax minimum and maximum of associated spatial location coordinates
#' @slot boundaries x and y coordinate adjustments (default to 0)
#' @slot scale_factor image scaling relative to spatial locations
#' @slot resolution spatial location units covered per pixel
#' @slot file_path file path to the image if given
#' @slot OS_platform Operating System to run Giotto analysis on
#' @details
#' [\strong{mg_object}] Core object is any image that can be read by the magick package
#'
#' [\strong{boundaries}] Boundary adjustments can be used to manually or
#' automatically through a script adjust the image with the spatial data.
#'
#'
#' @export
giottoImage <- setClass(
  Class = "giottoImage",

  slots = c(
    name = "ANY",
    mg_object = "ANY",
    minmax = "ANY",
    boundaries = "ANY",
    scale_factor = "ANY",
    resolution = "ANY",
    file_path = "ANY",
    OS_platform = "ANY"
  ),

  prototype = list(
    name = NULL,
    mg_object = NULL,
    minmax = NULL,
    boundaries = NULL,
    scale_factor = NULL,
    resolution = NULL,
    file_path = NULL,
    OS_platform = NULL
  )
)


#' show method for giottoImage class
#' @param object giottoImage object
#' @aliases show,giottoImage-method
#' @docType methods
#' @rdname show-methods

setMethod(
  f = "show",
  signature = "giottoImage",
  definition = function(object) {

    cat("An object of class '",  class(object), "' with name ", object@name, "\n \n")

    cat("Min and max values are: \n",
        "Max on x-axis: ", object@minmax[['xmax_sloc']], "\n",
        "Min on x-axis: ", object@minmax[['xmin_sloc']], "\n",
        "Max on y-axis: ", object@minmax[['ymax_sloc']], "\n",
        "Min on y-axis: ", object@minmax[['ymin_sloc']], "\n",
        "\n")

    cat("Boundary adjustment are: \n",
        "Max adjustment on x-axis: ", object@boundaries[['xmax_adj']], "\n",
        "Min adjustment on x-axis: ", object@boundaries[['xmin_adj']], "\n",
        "Max adjustment on y-axis: ", object@boundaries[['ymax_adj']], "\n",
        "Min adjustment on y-axis: ", object@boundaries[['ymin_adj']], "\n",
        "\n")

    cat("Boundaries are: \n",
        "Image x-axis max boundary: ", object@minmax[['xmax_sloc']] + object@boundaries[['xmax_adj']], "\n",
        "Image x-axis min boundary: ", object@minmax[['xmin_sloc']] - object@boundaries[['xmin_adj']], "\n",
        "Image y-axis max boundary: ", object@minmax[['ymax_sloc']] + object@boundaries[['ymax_adj']], "\n",
        "Image y-axis min boundary: ", object@minmax[['ymin_sloc']] - object@boundaries[['ymin_adj']], "\n",
        "\n")

    cat("Scale factor: \n")
    print(object@scale_factor)

    cat("\n Resolution: \n")
    print(object@resolution)
    
    cat("\n File Path: \n")
    print(object@file_path)

    # print(object@mg_object)

  }
)


# giottoLargeImage class ####

#' @title S4 giottoLargeImage Class
#' @description class to handle images too large to load in normally through magick
#' @keywords giotto, object, image
#' @slot name name of large Giotto image
#' @slot raster_object terra raster object
#' @slot overall_extent terra extent object covering the original extent of image
#' @slot scale_factor image scaling relative to spatial locations
#' @slot resolution spatial location units covered per pixel
#' @slot max_intensity value to set as maximum intensity in color scaling
#' @slot min_intensity minimum value found
#' @slot is_int values are integers
#' @slot file_path file path to the image if given
#' @slot OS_platform Operating System to run Giotto analysis on
#' @export
giottoLargeImage <- setClass(
  Class = "giottoLargeImage",

  slots = c(
    name = "ANY",
    raster_object = "ANY",
    overall_extent = "ANY",
    scale_factor = "ANY",
    resolution = "ANY",
    max_intensity = "ANY",
    min_intensity = "ANY",
    is_int = "ANY",
    file_path = "ANY",
    OS_platform = "ANY"
  ),

  prototype = list(
    name = NULL,
    raster_object = NULL,
    overall_extent = NULL,
    scale_factor = NULL,
    resolution = NULL,
    max_intensity = NULL,
    min_intensity = NULL,
    is_int = NULL,
    file_path = NULL,
    OS_platform = NULL
  )
)


#' show method for giottoLargeImage class
#' @param object giottoLargeImage object
#' @aliases show,giottoLargeImage-method
#' @docType methods
#' @rdname show-methods

setMethod(
  f = "show",
  signature = "giottoLargeImage",
  definition = function(object) {

    cat("An object of class '",  class(object), "' with name ", object@name, "\n \n")

    cat("Image boundaries are: \n")
    print(terra::ext(object@raster_object)[1:4])

    cat("Original image boundaries are: \n")
    print(object@overall_extent[1:4])

    cat("\n Scale factor: \n")
    print(object@scale_factor)

    cat("\n Resolution: \n")
    print(object@resolution)

    cat('\n Estimated maximum intensity is: ', object@max_intensity, ' \n',
        'Estimated minimum intensity is: ', object@min_intensity, ' \n')
    
    if(object@is_int == TRUE) cat('Values are integers')
    if(object@is_int == FALSE) cat('Values are floating point')
    
    cat('\n File path is: ', object@file_path, '\n')

  }
)


# giottoImage creation ####


#' @title createGiottoImage
#' @name createGiottoImage
#' @description Creates a giotto image that can be added to a Giotto object and/or used to add an image to the spatial plotting functions
#' @param gobject giotto object
#' @param spat_unit spatial unit
#' @param spatial_locs spatial locations (alternative if gobject = NULL)
#' @param spat_loc_name name of spatial locations within gobject
#' @param mg_object magick image object
#' @param name name for the image
#' @param image_transformations vector of sequential image transformations
#' @param negative_y Map image to negative y spatial values if TRUE during automatic alignment. Meaning that origin is in upper left instead of lower left.
#' @param do_manual_adj flag to use manual adj values instead of automatic alignment when given a gobject or spatlocs
#' @param xmax_adj adjustment of the maximum x-value to align the image
#' @param xmin_adj adjustment of the minimum x-value to align the image
#' @param ymax_adj adjustment of the maximum y-value to align the image
#' @param ymin_adj adjustment of the minimum y-value to align the image
#' @param scale_factor scaling of image dimensions relative to spatial coordinates
#' @param x_shift shift image in positive x direction
#' @param x_shift shift image in positive y direction
#' @param scale_x independently scale image in x direction
#' @param scale_y independently scale image in y direction
#' @param order perform scaling or adjustments and shifts first
#' @param xmin_set values to override image minmax spatial anchors when doing adjustments
#' @param xmax_set values to override image minmax spatial anchors when doing adjustments
#' @param ymin_set values to override image minmax spatial anchors when doing adjustments
#' @param ymax_set values to override image minmax spatial anchors when doing adjustments
#' @param verbose be verbose
#' @details image_transformations: transformation options from magick library
#' [\strong{flip_x_axis}] flip x-axis (\code{\link[magick]{image_flop}})
#' [\strong{flip_y_axis}] flip y-axis (\code{\link[magick]{image_flip}})
#' Example: image_transformations = c(flip_x_axis, flip_y_axis); first flip x-axis and then y-axis
#' @return a giottoImage object
#' @export
createGiottoImage = function(gobject = NULL,
                             spat_unit = NULL,
                             spatial_locs = NULL,
                             spat_loc_name = NULL,
                             mg_object,
                             name = 'image',
                             image_transformations = NULL,
                             negative_y = TRUE,
                             do_manual_adj = FALSE,
                             xmax_adj = 0,
                             xmin_adj = 0,
                             ymax_adj = 0,
                             ymin_adj = 0,
                             scale_factor = 1,
                             x_shift = NULL,
                             y_shift = NULL,
                             scale_x = NULL,
                             scale_y = NULL,
                             order = c('first_scale','first_adj'),
                             xmin_set = NULL,
                             xmax_set = NULL,
                             ymin_set = NULL,
                             ymax_set = NULL,
                             verbose = TRUE) {

  # Check params
  order = match.arg(order, choices = c('first_scale','first_adj'))
  scale_factor = c(x = scale_factor, y = scale_factor)

  # create minimum giotto
  g_image = giottoImage(name = name,
                        mg_object = NULL,
                        minmax = NULL,
                        boundaries = NULL,
                        scale_factor = NULL,
                        resolution = NULL,
                        file_path = NULL,
                        OS_platform = .Platform[['OS.type']])


  ## 1.a. check magick image object
  if(!methods::is(mg_object, 'magick-image')) {
    if(file.exists(mg_object)) {
      g_image@file_path = mg_object
      mg_object = try(magick::image_read(mg_object))
      if(class(mg_object) == 'try-error') {
        stop(mg_object, ' can not be read by magick::image_read() \n')
      }
    } else {
      stop("mg_object needs to be an image object 'magick-image' from the magick package or \n
           an existing path that can be read by magick::image_read()")
    }
  }

  ## 1.b. check colorspace
  info = magick::image_info(mg_object)
  mg_colorspace = info$colorspace
  if(mg_colorspace == 'Gray') {
    mg_object = magick::image_convert(mg_object, colorspace = 'rgb')
  }

  ## 1.c. perform transformations if found
  if(!is.null(image_transformations)) {
    for(transf in image_transformations) {
      if(transf == 'flip_x_axis') {
        mg_object = magick::image_flop(mg_object)
      } else if(transf == 'flip_y_axis') {
        mg_object = magick::image_flop(mg_object)
      } else {
        cat(transf, ' is not a supported transformation, see details \n')
      }
    }
  }

  g_image@mg_object = mg_object

  ## 2. spatial minmax and adjustments -- manual OR by image dimensions (auto)
  if(verbose == TRUE) {
    if(do_manual_adj == TRUE) cat('do_manual_adj == TRUE \n','Boundaries will be adjusted by given values.\n')
  }
  # If spatlocs or gobject supplied, minmax values will always be generated
  # If do_manual_adj == TRUE, bypass followup automatic boundary value generation
  if(!is.null(gobject)) {

    # Get spatial locations (or automatically take first available)
    spatlocs = get_spatial_locations(gobject = gobject,
                                     spat_unit = spat_unit,
                                     spat_loc_name = spat_loc_name)

    # spatlocs = gobject@spatial_locs[['raw']]

    # Find g_image minmax (spatial) from spatial_locs in gobject
    my_xmin = min(spatlocs$sdimx)
    my_xmax = max(spatlocs$sdimx)
    my_ymin = min(spatlocs$sdimy)
    my_ymax = max(spatlocs$sdimy)

    if(do_manual_adj == FALSE) {
      # find automatic adjustment values
      img_minmax = get_img_minmax(mg_img = mg_object,
                                  negative_y = negative_y)
      adj_values = get_adj_rescale_img(img_minmax = img_minmax,
                                       spatial_locs = spatlocs,
                                       scale_factor = scale_factor)
      # Automatic g_image@boundaries values
      xmax_adj = as.numeric(adj_values[['xmax_adj_orig']])
      xmin_adj = as.numeric(adj_values[['xmin_adj_orig']])
      ymax_adj = as.numeric(adj_values[['ymax_adj_orig']])
      ymin_adj = as.numeric(adj_values[['ymin_adj_orig']])
    }


  } else if(!is.null(spatial_locs)) {
    spatlocs = spatial_locs
    if(!all(c('sdimx','sdimy') %in% colnames(spatlocs))) {
      stop('spatial_locs needs to be data.frame-like object with a sdimx and sdimy column')
    }
    # Find g_image minmax (spatial) from spatial_locs argument
    my_xmin = min(spatlocs$sdimx)
    my_xmax = max(spatlocs$sdimx)
    my_ymin = min(spatlocs$sdimy)
    my_ymax = max(spatlocs$sdimy)

    if(do_manual_adj == FALSE) {
      #find auto adjustment values
      img_minmax = get_img_minmax(mg_img = mg_object,
                                  negative_y = negative_y)
      adj_values = get_adj_rescale_img(img_minmax = img_minmax,
                                       spatial_locs = spatlocs,
                                       scale_factor = scale_factor)
      # Automatic g_image@boundaries values
      xmax_adj = as.numeric(adj_values[['xmax_adj_orig']])
      xmin_adj = as.numeric(adj_values[['xmin_adj_orig']])
      ymax_adj = as.numeric(adj_values[['ymax_adj_orig']])
      ymin_adj = as.numeric(adj_values[['ymin_adj_orig']])
    }


  } else {
    if(verbose == TRUE) {
      warning('gobject or spatial locations are not provided \n',
              'Arbitrary values will be given \n')
    }
    # Default g_image@minmax values if no spatial_locs provided
    my_xmin = 0; my_xmax = 10; my_ymin = 0; my_ymax = 10

  }

  # Set minmax and boundary values for return
  g_image@minmax = c('xmax_sloc' = my_xmax,
                     'xmin_sloc' = my_xmin,
                     'ymax_sloc' = my_ymax,
                     'ymin_sloc' = my_ymin)

  ## if do_manual == TRUE, boundary values are those given or defaulted as arguments
  ## if do_manual == FALSE, boundary values are taken from automatic processes above
  g_image@boundaries = c('xmax_adj' = xmax_adj,
                         'xmin_adj' = xmin_adj,
                         'ymax_adj' = ymax_adj,
                         'ymin_adj' = ymin_adj)

  # scale factor and resolution values for return
  g_image@scale_factor = scale_factor
  g_image@resolution = 1/scale_factor

  # Apply any additional manual adjustments through updateGiottoImage
  if(do_manual_adj == TRUE) {
    if(length(c(x_shift,
                y_shift,
                scale_x,
                scale_y,
                xmin_set,
                xmax_set,
                ymin_set,
                ymax_set)) > 0) {
      g_image = updateGiottoImageMG(giottoImage = g_image,
                                    return_gobject = FALSE,
                                    xmax_adj = xmax_adj,
                                    xmin_adj = xmin_adj,
                                    ymax_adj = ymax_adj,
                                    ymin_adj = ymin_adj,
                                    x_shift = x_shift,
                                    y_shift = y_shift,
                                    scale_factor = scale_factor,
                                    order = order,
                                    xmin_set = xmin_set,
                                    xmax_set = xmax_set,
                                    ymin_set = ymin_set,
                                    ymax_set = ymax_set,
                                    verbose = FALSE)
    }
  }

  # image object
  return(g_image)
}



# giottoLargeImage creation ####


#' @title createGiottoLargeImage
#' @name createGiottoLargeImage
#' @description Creates a large giotto image that can be added to a Giotto subcellular object. Generates deep copy of SpatRaster
#' @param raster_object terra SpatRaster image object
#' @param name name for the image
#' @param negative_y Map image to negative y spatial values if TRUE. Meaning that origin is in upper left instead of lower left.
#' @param extent SpatExtent object to assign spatial extent. Takes priority unless use_rast_ext is TRUE.
#' @param use_rast_extent Use extent from input raster object
#' @param image_transformations vector of sequential image transformations - under construction
#' @param xmax_adj adjustment of the maximum x-value to align the image
#' @param xmin_adj adjustment of the minimum x-value to align the image
#' @param ymax_adj adjustment of the maximum y-value to align the image
#' @param ymin_adj adjustment of the minimum y-value to align the image
#' @param scale_factor scaling of image dimensions relative to spatial coordinates
#' @param verbose be verbose
#' @return a giottoLargeImage object
#' @export
createGiottoLargeImage = function(raster_object,
                                  name = 'image',
                                  negative_y = TRUE,
                                  extent = NULL,
                                  use_rast_ext = FALSE,
                                  image_transformations = NULL,
                                  xmax_bound = NULL,
                                  xmin_bound = NULL,
                                  ymax_bound = NULL,
                                  ymin_bound = NULL,
                                  scale_factor = 1,
                                  verbose = TRUE) {

  # create minimum giotto
  g_imageL = giottoLargeImage(name = name,
                              raster_object = NULL,
                              overall_extent = NULL,
                              scale_factor = NULL,
                              resolution = NULL,
                              file_path = NULL,
                              OS_platform = .Platform[['OS.type']])


  ## 1. check raster object and load as SpatRaster if necessary
  if(!methods::is(raster_object, 'SpatRaster')) {
    if(file.exists(raster_object)) {
      g_imageL@file_path = raster_object
      raster_object = try(suppressWarnings(terra::rast(x = raster_object)))
      if(class(raster_object) == 'try-error') {
        stop(raster_object, ' can not be read by terra::rast() \n')
      }
    } else {
      stop("raster_object needs to be a'SpatRaster' object from the terra package or \n
           an existing path that can be read by terra::rast()")
    }
  }

  # Prevent updates to original raster object input
  raster_object = terra::deepcopy(raster_object)


  ## 2. image bound spatial extent
  if(use_rast_ext == TRUE) {
    extent = terra::ext(raster_object)
    if(verbose == TRUE) cat('use_rast_ext == TRUE, extent from input raster_object will be used.')
  }

  # By extent object (priority)
  if(!is.null(extent)) {
    if(methods::is(extent, 'SpatExtent')) {
      terra::ext(raster_object) = extent
    } else {
      stop('extent argument only accepts terra SpatExtent objects')
    }
  } else { # OR by manual OR by image dimensions (auto)

    # Check if manual adj values were given
    # Assign default values for any that were not manually given
    if(all(is.null(xmax_bound),
           is.null(xmin_bound),
           is.null(ymax_bound),
           is.null(ymin_bound))) {
      im_dim = dim(raster_object)[2:1]

      # Apply scale_factor
      im_dim = im_dim * scale_factor

      # Automatic extent values
      xmax_bound = im_dim[1]
      xmin_bound = 0
      if(negative_y == TRUE) {
        ymax_bound = 0
        ymin_bound = -im_dim[2]
      } else if(negative_y == FALSE) {
        ymax_bound = im_dim[2]
        ymin_bound = 0
      }

    } else {
      # Manual extent values
      if(is.null(xmax_bound) == TRUE) xmax_bound = 1
      if(is.null(xmin_bound) == TRUE) xmin_bound = 0
      if(negative_y == TRUE) {
        if(is.null(ymax_bound) == TRUE) ymax_bound = 0
        if(is.null(ymin_bound) == TRUE) ymin_bound = -1
      } else if(negative_y == FALSE) {
        if(is.null(ymax_bound) == TRUE) ymax_bound = 1
        if(is.null(ymin_bound) == TRUE) ymin_bound = 0
      }

    }
    terra::ext(raster_object) = c(xmin_bound,xmax_bound,ymin_bound,ymax_bound)
  }


  ## 3. Assign raster_object to giottoLargeImage
  g_imageL@raster_object = raster_object

  ## 4. scale factor and resolution values
  g_imageL@resolution = terra::res(g_imageL@raster_object) # (x,y)
  names(g_imageL@resolution) = c('x','y')
  g_imageL@scale_factor = (1/g_imageL@resolution)

  ## 5. Get image characteristics by sampling
  sampleValues = na.omit(terra::spatSample(raster_object,
                                           size = 5000, # Defines the rough maximum of pixels allowed when resampling
                                           method = 'regular',
                                           value = TRUE))
  if(nrow(sampleValues) == 0) {
    if(verbose == TRUE) cat('No values discovered when sampling for image characteristics')
  } else {
    # get intensity range
    intensityRange = range(sampleValues)
    g_imageL@max_intensity = intensityRange[2]
    g_imageL@min_intensity = intensityRange[1]
    # find out if image is int or floating point
    is_int = identical(sampleValues, round(sampleValues))
    if(is_int == TRUE) {
      g_imageL@is_int = TRUE
    } else {
      g_imageL@is_int = FALSE
    }
  }


  ## 6. extent object
  g_imageL@overall_extent = terra::ext(raster_object)

  ## 7. return image object
  return(g_imageL)
}


#' @title createGiottoLargeImageList
#' @name createGiottoLargeImageList
#' @description Creates a list of large giotto images that can be added to a Giotto object. Generates deep copy of SpatRaster
#' @param raster_objects vector of image paths or terra SpatRaster image objects
#' @param names vector of names for the images
#' @param negative_y Map image to negative y spatial values if TRUE. Meaning that origin is in upper left instead of lower left.
#' @param extent SpatExtent object to assign spatial extent. Takes priority unless use_rast_ext is TRUE.
#' @param use_rast_extent Use extent from input raster object
#' @param image_transformations vector of sequential image transformations - under construction
#' @param xmax_adj adjustment of the maximum x-value to align the image
#' @param xmin_adj adjustment of the minimum x-value to align the image
#' @param ymax_adj adjustment of the maximum y-value to align the image
#' @param ymin_adj adjustment of the minimum y-value to align the image
#' @param scale_factor scaling of image dimensions relative to spatial coordinates
#' @param verbose be verbose
#' @details See \code{\link{createGiottoLargeImage}}
#' @return a list with giottoLargeImage objects
#' @export
createGiottoLargeImageList = function(raster_objects,
                                      names = 'image',
                                      negative_y = TRUE,
                                      extent = NULL,
                                      use_rast_ext = FALSE,
                                      image_transformations = NULL,
                                      xmax_bound = NULL,
                                      xmin_bound = NULL,
                                      ymax_bound = NULL,
                                      ymin_bound = NULL,
                                      scale_factor = 1,
                                      verbose = TRUE) {

  l_images = length(raster_objects)
  l_image_names = length(unique(names))

  if(l_image_names != l_image_names) {
    stop('length of raster_objects and unique names must be the same')
  }

  result_list = list()

  for(i in 1:l_images) {

    image_res = createGiottoLargeImage(raster_object = raster_objects[[i]],
                                      name = names[[i]],
                                      negative_y = negative_y,
                                      extent = extent,
                                      use_rast_ext = use_rast_ext,
                                      image_transformations = image_transformations,
                                      xmax_bound = xmax_bound,
                                      xmin_bound = xmin_bound,
                                      ymax_bound = ymax_bound,
                                      ymin_bound = ymin_bound,
                                      scale_factor = scale_factor,
                                      verbose = verbose)

    result_list[[i]] = image_res

  }

  return(result_list)

}

# giottoImage or magick tools ####

#' @title convert_mgImage_to_array_DT
#' @name convert_mgImage_to_array_DT
#' @description converts a magick image object to a data.table
#' @param mg_object magick image or Giotto image object
#' @return data.table with image pixel information
#' @keywords internal
convert_mgImage_to_array_DT = function(mg_object) {

  if(methods::is(mg_object, 'giottoImage')) {
    mg_object = mg_object@mg_object
  }

  # data.table variables
  RGB = c.1 = c.2 = c.3 = NULL

  # convert magick object to an array
  num_res = as.numeric(mg_object[[1]])
  num_res_m = data.table::as.data.table(reshape2::melt(num_res))
  colnames(num_res_m) = c('x', 'y', 'c', 'color')
  array_dt = data.table::dcast.data.table(num_res_m, value.var = 'color', formula = 'x+y~c')
  colnames(array_dt) = c('x', 'y', 'c.1', 'c.2', 'c.3')
  array_dt[, RGB := grDevices::rgb(c.1, c.2, c.3)]

  return(array_dt)
}


#' @title estimateImageBg
#' @name estimateImageBg
#' @description helps to estimate which color is the background color of your plot
#' @param mg_object magick image or Giotto image object
#' @param top_color_range top possible background colors to return
#' @return vector of pixel color frequencies and an associated barplot
#' @export
estimateImageBg = function(mg_object, top_color_range = 1:50) {

  if(methods::is(mg_object, 'giottoImage')) {
    mg_object = mg_object@mg_object
  }

  arrayDT = convert_mgImage_to_array_DT(mg_object = mg_object)
  sort_table = sort(table(arrayDT$RGB), decreasing = T)
  graphics::barplot(sort_table[top_color_range], col=names(sort_table[top_color_range]))

  cat('Most abundant pixel colors: \n')
  print(sort_table[top_color_range])
}


#' @title changeImageBg
#' @name changeImageBg
#' @description Function to change the background color of a magick image plot to another color
#' @param mg_object magick image or giotto image object
#' @param bg_color estimated current background color
#' @param perc_range range around estimated background color to include (percentage)
#' @param new_color new background color
#' @param new_name change name of Giotto image
#' @return magick image or giotto image object with updated background color
#' @export
changeImageBg = function(mg_object,
                         bg_color,
                         perc_range = 10,
                         new_color = '#FFFFFF',
                         new_name = NULL) {

  if(methods::is(mg_object, 'giottoImage')) {
    is_g_image = TRUE
    g_image = mg_object
    mg_object = mg_object@mg_object
  } else {
    is_g_image = FALSE
  }

  if(!methods::is(mg_object, 'magick-image')) {
    stop("mg_object needs to be a giottImage or a 'magick-image' object from the magick package")
  }

  # new background color
  new_rbg_color = grDevices::col2rgb(new_color)/255

  # current background limits
  rbgcolors = grDevices::col2rgb(bg_color)/255
  perc_range_min = rbgcolors - (rbgcolors/100)*perc_range
  perc_range_max = rbgcolors + (rbgcolors/100)*perc_range

  # convert magick image to array data.table
  arrayDT = convert_mgImage_to_array_DT(mg_object = mg_object)

  # create new background
  c1_min = perc_range_min[1,1]
  c2_min = perc_range_min[2,1]
  c3_min = perc_range_min[3,1]
  c1_max = perc_range_max[1,1]
  c2_max = perc_range_max[2,1]
  c3_max = perc_range_max[3,1]

  c1_new = new_rbg_color[1,1]
  c2_new = new_rbg_color[2,1]
  c3_new = new_rbg_color[3,1]

  # find background color pixels

  # data.table variables
  c.1 = c.2 = c.3 = NULL

  c1_ind = arrayDT[['c.1']] > c1_min & arrayDT[['c.1']] < c1_max
  c2_ind = arrayDT[['c.2']] > c2_min & arrayDT[['c.2']] < c2_max
  c3_ind = arrayDT[['c.3']] > c3_min & arrayDT[['c.3']] < c3_max
  c_ind = c1_ind*c2_ind*c3_ind

  # data.table variables
  c1 = c2 = c3 = NULL

  # replace old background with new background
  arrayDT[, 'c1' := ifelse(c_ind == T, c1_new, c.1)]
  arrayDT[, 'c2' := ifelse(c_ind == T, c2_new, c.2)]
  arrayDT[, 'c3' := ifelse(c_ind == T, c3_new, c.3)]


  # data.table variables
  x = y = NULL

  # setorder for x and y coordinates
  data.table::setorder(arrayDT, y, x)

  # convert array_dt to array and then to magick image object
  original_width = magick::image_info(mg_object)[2]
  original_heigth = magick::image_info(mg_object)[3]
  myarray = array(as.vector(as.matrix(arrayDT[,.(c1, c2, c3)])), dim = c(original_width, original_heigth, 3))
  new_mg_object = magick::image_read(myarray)

  # return magick or giotto image object
  if(is_g_image == TRUE) {
    if(!is.null(new_name)) g_image$name = new_name
    g_image@mg_object = new_mg_object
    return(g_image)
  } else {
    return(new_mg_object)
  }
}


#' @name createGiottoImageOLD
#' @description Creates a giotto image that can be added to a Giotto object and/or used to add an image to the spatial plotting functions
#' @param gobject giotto object
#' @param spatial_locs spatial locations (alternative if giobject = NULL)
#' @param mg_object magick image object
#' @param name name for the image
#' @param xmax_adj adjustment of the maximum x-value to align the image
#' @param xmin_adj adjustment of the minimum x-value to align the image
#' @param ymax_adj adjustment of the maximum y-value to align the image
#' @param ymin_adj adjustment of the minimum y-value to align the image
#' @return a giotto image object
#' @export
createGiottoImageOLD = function(gobject = NULL,
                             spatial_locs = NULL,
                             mg_object,
                             name = 'image',
                             xmax_adj = 0,
                             xmin_adj = 0,
                             ymax_adj = 0,
                             ymin_adj = 0) {

  if(!methods::is(mg_object, 'magick-image')) {
    if(file.exists(mg_object)) {
      mg_object = try(magick::image_read(mg_object))
      if(class(mg_object) == 'try-error') {
        stop(mg_object, ' can not be read by magick::image_read() \n')
      }
    } else {
      stop("mg_object needs to be an image object 'magick-image' from the magick package or \n
           an existig path that can be read by magick::image_read()")
    }
  }

  # min and max
  if(!is.null(gobject)) {
    spatlocs = gobject@spatial_locs
  } else if(!is.null(spatial_locs)) {
    spatlocs = spatial_locs
  } else {
    stop('gobject or spatial locations need to be provided')
  }

  my_xmin = min(spatlocs$sdimx)
  my_xmax = max(spatlocs$sdimx)
  my_ymin = min(spatlocs$sdimy)
  my_ymax = max(spatlocs$sdimy)

  # image object
  imageObj = list(name = name,
                  mg_object = mg_object,
                  minmax = c('xmax_sloc' = my_xmax, 'xmin_sloc' = my_xmin,
                             'ymax_sloc' = my_ymax, 'ymin_sloc' = my_ymin),
                  boundaries = c('xmax_adj' = xmax_adj, 'xmin_adj' = xmin_adj,
                                 'ymax_adj' = ymax_adj, 'ymin_adj' = ymin_adj))

  class(imageObj) <- append(class(imageObj), 'imageGiottoObj')
  return(imageObj)
}


#' @title addGiottoImageMG
#' @name addGiottoImageMG
#' @description Adds giotto image objects to your giotto object
#' @param gobject giotto object
#' @param images list of giotto image objects, see \code{\link{createGiottoImage}}
#' @param spat_unit spatial unit
#' @param spat_loc_name provide spatial location slot in Giotto to align images. Defaults to first one
#' @param scale_factor provide scale of image pixel dimensions relative to spatial coordinates.
#' @param negative_y Map image to negative y spatial values if TRUE during automatic alignment. Meaning that origin is in upper left instead of lower left.
#' @return an updated Giotto object with access to the list of images
#' @export
addGiottoImageMG = function(gobject,
                            images,
                            spat_unit = NULL,
                            spat_loc_name = NULL,
                            scale_factor = NULL,
                            negative_y = TRUE) {

  # 0. check params
  if(is.null(gobject)) stop('The giotto object that will be updated needs to be provided')

  if(is.null(images)) stop('The giotto image(s) that will be added needs to be provided')

  if(is.null(spat_loc_name)) {
    if(!is.null(gobject@spatial_locs)) {
      spat_loc_name = list_spatial_locations(gobject = gobject, spat_unit = spat_unit)[1,]
    } else {
      spat_loc_name = NULL
      cat('No spatial locations have been found \n')
    }
  }

  ext_scale_factor = FALSE
  if(!is.null(scale_factor)) {

    if(!is.numeric(scale_factor)) stop ('Given scale_factor(s) must be numeric')

    if((length(scale_factor) == length(images)) || length(scale_factor) == 1) {
      cat('scale_factor(s) external to giottoImage have been given and will be used')
      ext_scale_factor = TRUE
    } else {
      stop('if scale_factor is given, it must be a numeric with either a single value or as many values as there are images are provided')
    }
  }

  # 1. expand scale_factors
  if(ext_scale_factor == TRUE) {
    if(length(scale_factor == 1)) {
      scale_factor = rep(scale_factor, length(images))
    }
  }


  # 2. Add image with for loop
  for(image_i in 1:length(images)) {

    im = images[[image_i]]

    if(methods::is(im, 'giottoImage')) {
      im_name = im@name

      all_im_names = names(gobject@images)

      if(im_name %in% all_im_names) {
        cat('\n ', im_name, ' has already been used, will be overwritten \n')
      }

      # 3. Update boundaries if not already done during createGiottoImage() due to lack of spatlocs and gobject
      if(sum(im@boundaries == c(0,0,0,0)) == 4 && sum(im@minmax == c(10,0,10,0)) == 4) {
        if(!is.null(spat_loc_name)) { # A check for the first available spatloc was already done
          spatlocs = get_spatial_locations(gobject = gobject,
                                           spat_unit = spat_unit,
                                           spat_loc_name = spat_loc_name)

          #Find spatial minmax values
          xmin_sloc = min(spatlocs$sdimx)
          xmax_sloc = max(spatlocs$sdimx)
          ymin_sloc = min(spatlocs$sdimy)
          ymax_sloc = max(spatlocs$sdimy)

          #Find adjustment values
          img_minmax = get_img_minmax(mg_img = im@mg_object,
                                      negative_y = negative_y)
          if(ext_scale_factor == TRUE) {
            adj_values = get_adj_rescale_img(img_minmax = img_minmax,
                                             spatial_locs = spatlocs,
                                             scale_factor = scale_factor[[image_i]])
          } else if (ext_scale_factor == FALSE) {
            adj_values = get_adj_rescale_img(img_minmax = img_minmax,
                                             spatial_locs = spatlocs,
                                             scale_factor = im@scale_factor[[spat_loc_name]])
          }

          #Add minmax values to giottoImage@minmax
          im@minmax = c('xmax_sloc' = xmax_sloc,
                        'xmin_sloc' = xmin_sloc,
                        'ymax_sloc' = ymax_sloc,
                        'ymin_sloc' = ymin_sloc)

          #Add adjustment values to giottoImage@boundaries
          im@boundaries = c('xmax_adj' = as.numeric(adj_values[['xmax_adj_orig']]),
                            'xmin_adj' = as.numeric(adj_values[['xmin_adj_orig']]),
                            'ymax_adj' = as.numeric(adj_values[['ymax_adj_orig']]),
                            'ymin_adj' = as.numeric(adj_values[['ymin_adj_orig']]))

          # Inherit external scaling factors if given
          if(ext_scale_factor == TRUE) {
            im@scale_factor[[spat_loc_name]] = scale_factor[[image_i]]
            im@resolution[[spat_loc_name]] = 1/(scale_factor[[image_i]])
          }
          ## Externally given scale_factors will only be written in/used if boundary adj values are not pre-existing
        }

      }

      # 4. Add giottoImage to gobject
      gobject@images[[im_name]] = im

    } else {
      warning('image [',image_i,'] is not a giotto image object')
    }
  }

  return(gobject)

}



#' @title addGiottoImageToSpatPlot
#' @name addGiottoImageToSpatPlot
#' @description Add a giotto image to a spatial ggplot object post creation
#' @param spatpl a spatial ggplot object
#' @param gimage a giotto image, see \code{\link{createGiottoImage}}
#' @return an updated spatial ggplot object
#' @export
addGiottoImageToSpatPlot = function(spatpl = NULL,
                                    gimage = NULL) {


  if(is.null(spatpl) | is.null(gimage)) {
    stop('A spatial ggplot object and a giotto image need to be given')
  }

  # extract min and max from object
  my_xmax = gimage@minmax[1]
  my_xmin = gimage@minmax[2]
  my_ymax = gimage@minmax[3]
  my_ymin = gimage@minmax[4]

  # convert giotto image object into array
  img_array = as.numeric(gimage@mg_object[[1]])

  # extract adjustments from object
  xmax_b = gimage@boundaries[1]
  xmin_b = gimage@boundaries[2]
  ymax_b = gimage@boundaries[3]
  ymin_b = gimage@boundaries[4]

  newpl = spatpl + annotation_raster(img_array,
                                     xmin = my_xmin-xmin_b, xmax = my_xmax+xmax_b,
                                     ymin = my_ymin-ymin_b, ymax = my_ymax+ymax_b)

  # move image to background
  nr_layers = length(newpl$layers)
  newpl$layers = c(newpl$layers[[nr_layers]], newpl$layers[1:(nr_layers-1)])

  return(newpl)

}




#' @title updateGiottoImageMG
#' @name updateGiottoImageMG
#' @description Updates the boundaries of a giotto image attached to a giotto object
#' @param gobject giotto object
#' @param image_name name of giottoImage
#' @param giottoImage giottoImage object
#' @param xmax_adj adjustment of the maximum x-value to align the image
#' @param xmin_adj adjustment of the minimum x-value to align the image
#' @param ymax_adj adjustment of the maximum y-value to align the image
#' @param ymin_adj adjustment of the minimum y-value to align the image
#' @param x_shift shift entire image in positive x direction
#' @param y_shift shift entire image in positive y direction
#' @param scale_factor set scale_x and scale_y at the same time
#' @param scale_x independently scale x axis image mapping from origin
#' @param scale_y independently scale y axis image mapping from origin
#' @param order perform fine adjustments (adjustments and shifts) or scaling first
#' @param xmin_set set image xmin boundary. Applied before adjustments
#' @param xmax_set set image xmax boundary. Applied before adjustments
#' @param ymin_set set image ymin boundary. Applied before adjustments
#' @param ymax_set set image ymax boundary. Applied before adjustments
#' @param return_gobject return a giotto object
#' @param verbose be verbose
#' @return a giotto object or an updated giotto image if return_gobject = F
#' @export
updateGiottoImageMG = function(gobject = NULL,
                               image_name = NULL,
                               giottoImage = NULL,
                               xmax_adj = 0,
                               xmin_adj = 0,
                               ymax_adj = 0,
                               ymin_adj = 0,
                               x_shift = 0,
                               y_shift = 0,
                               scale_factor = NULL,
                               scale_x = 1,
                               scale_y = 1,
                               order = c('first_adj','first_scale'),
                               xmin_set = NULL,
                               xmax_set = NULL,
                               ymin_set = NULL,
                               ymax_set = NULL,
                               return_gobject = TRUE,
                               verbose = TRUE) {


  # 0. Check params
  # Check input image
  if(is.null(gobject)) {
    if(is.null(giottoImage)) stop('Image to be updated must be given as gobject AND image_name OR giottoImage argument(s) \n')
    if(verbose == TRUE) cat('gobject argument not given \n return_gobject set to FALSE \n')
    return_gobject = FALSE
  }
  if(is.null(giottoImage) && is.null(image_name)) stop('The name of the giotto image that will be updated needs to be provided \n')

  if(!is.null(giottoImage)) {
    if(!methods::is(giottoImage, 'giottoImage')) stop('giottoImage argument only accepts giottoImage objects \n')
    if(verbose == TRUE && !is.null(gobject)) cat('giottoImage argument is given and will take priority \n return_gobject set to FALSE \n')
    return_gobject = FALSE
  }

  # Check scalefactors
  if(!is.null(scale_factor)) scale_x = scale_y = scale_factor

  # Check spatial anchor values
  spatAnchor = c('xmax_sloc' = xmax_set,
                 'xmin_sloc' = xmin_set,
                 'ymax_sloc' = ymax_set,
                 'ymin_sloc' = ymin_set)
  if(length(spatAnchor) < 4 && length(spatAnchor) > 0) stop('If set arguments are being used, all four must be given \n')
  if(!is.null(spatAnchor)) {
    if(xmax_set < xmin_set) stop('xmax_set must be greater than xmin_set \n')
    if(ymax_set < ymin_set) stop('ymax_set must be greater than ymin_set \n')
  }

  # Find order of adjust and scaling
  order = match.arg(order, choices = c('first_adj','first_scale'))


  # 1. get giottoImage if necessary
  if(is.null(giottoImage)) {
    if(!is.null(gobject) && !is.null(image_name)) {
      g_image = getGiottoImage(gobject = gobject,
                               image_name = image_name)
    } else {
      stop('either a giottoImage or both the gobject and name of the giottoImage must be given. \n')
    }
  }


  # 2. Find minmax spatial anchor values
  if(is.null(spatAnchor)) {
    spatAnchor = g_image@minmax
  }

  # Perform scale if first_scale
  if(order == 'first_scale') {
    spatAnchor = spatAnchor * c(scale_x, scale_x, scale_y, scale_y)
  }

  # 3. Prepare adjustment values
  # Apply shifts
  xmin_adj = xmin_adj - x_shift
  xmax_adj = xmax_adj + x_shift
  ymin_adj = ymin_adj - y_shift
  ymax_adj = ymax_adj + y_shift

  # Find final bounds
  xmax_final = spatAnchor[['xmax_sloc']] + xmax_adj
  xmin_final = spatAnchor[['xmin_sloc']] - xmin_adj
  ymax_final = spatAnchor[['ymax_sloc']] + ymax_adj
  ymin_final = spatAnchor[['ymin_sloc']] - ymin_adj

  # Perform scale if first_adj
  if(order == 'first_adj') {
    xmax_final = xmax_final * scale_x
    xmin_final = xmin_final * scale_x
    ymax_final = ymax_final * scale_y
    ymin_final = ymin_final * scale_y
  }

  # Find final adj values
  xmax_adj = xmax_final - g_image@minmax[['xmax_sloc']]
  xmin_adj = g_image@minmax[['xmin_sloc']] - xmin_final
  ymax_adj = ymax_final - g_image@minmax[['ymax_sloc']]
  ymin_adj = g_image@minmax[['ymin_sloc']] - ymin_final


  # 4. Update the boundaries
  g_image@boundaries = c('xmax_adj' = xmax_adj,
                         'xmin_adj' = xmin_adj,
                         'ymax_adj' = ymax_adj,
                         'ymin_adj' = ymin_adj)

  # 5. Update the scalefactors for x and y
  x_range = xmax_final - xmin_final
  y_range = ymax_final - ymin_final
  im_dims = magick::image_info(g_image@mg_object)
  x_scalefactor = im_dims[['width']] / x_range
  y_scalefactor = im_dims[['height']] / y_range

  g_image@scale_factor = c('x' = x_scalefactor, 'y' = y_scalefactor)
  g_image@resolution = (1/g_image@scale_factor)

  if(return_gobject == TRUE) {
    gobject@images[[image_name]] = g_image
    return(gobject)
  } else {
    return(g_image)
  }

}



#' @title getGiottoImage
#' @name getGiottoImage
#' @description get a giotto image from a giotto object
#' @param gobject giotto object
#' @param image_name name of giotto image \code{\link{showGiottoImageNames}}
#' @return a giotto image
#' @keywords internal
get_GiottoImage_MG = function(gobject = NULL,
                              image_name = NULL) {

  if(is.null(gobject)) stop('The giotto object holding the giottoImage needs to be provided \n')
  g_image_names = names(gobject@images)
  if(is.null(g_image_names)) stop('No giottoImages have been found \n')

  if(is.null(image_name)) {
    image_name = g_image_names[1]
  }

  if(!image_name %in% g_image_names) stop(image_name, ' was not found among the image names, see showImageNames()')

  g_image = gobject@images[[image_name]]

  return(g_image)
}


#' @title plot_giottoImage_MG
#' @name plot_giottoImage_MG
#' @description get and plot a giottoImage either directly or from a giotto object
#' @param gobject giotto object
#' @param image_name name of giotto image \code{\link{showGiottoImageNames}}
#' @param giottoImage giottoImage object
#' @return plot
#' @keywords internal
plot_giottoImage_MG = function(gobject = NULL,
                               image_name = NULL,
                               giottoImage = NULL) {

  if(!is.null(giottoImage)) {
    graphics::plot(giottoImage@mg_object)
  } else {
    if(is.null(gobject)) stop('The giotto object that will be updated needs to be provided \n')
    if(is.null(image_name)) stop('The name of the giotto image that will be updated needs to be provided \n')

    g_image_names = names(gobject@images)
    if(!image_name %in% g_image_names) stop(image_name, ' was not found among the image names, see showImageNames()')

    graphics::plot(gobject@images[[image_name]]@mg_object)
  }

}






# giottoLargeImage or terra tools ####


#' @title stitchGiottoLargeImage
#' @name stitchGiottoLargeImage
#' @description stitch multiple giottoLargeImages into single giottoLargeImage. Time consuming & save location recommended.
#' @param largeImage_list list of giottoLargeImage objects
#' @param gobject_list list of gobjects containing giottoLargeImages
#' @param largeImage_nameList list of names of giottoLargeImages within gobjects
#' @param FOV_positions dataframe of FOV positions. Values (if any) are directly added to current image mapping
#' @param FOV_xcol column name for FOV position x values
#' @param FOV_ycol column name for FOV position y values
#' @param FOV_inverty make FOV y position values negative
#' @param method method of stitching images (mosaic: average overlapping area values, merge: values get priority by order given)
#' @param round_positions round image positions. May be necessary to run.
#' @param filename file name to write the stitched image to. Defaults to "save_dir/stitch.tif" if save_dir param is found in first gobject Giotto instructions
#' @param dataType (optional) values for dataType are "INT1U", "INT2U", "INT2S", "INT4U", "INT4S", "FLT4S", "FLT8S". The first three letters indicate whether the dataType is integer (whole numbers) of a real number (decimal numbers), the fourth character indicates the number of bytes used (allowing for large numbers and/or more precision), and the "S" or "U" indicate whether the values are signed (both negative and positive) or unsigned (positive values only).
#' @param fileType (optional) image format (e.g. .tif) If not given, defaults to format given in the filename
#' @param dryRun display placeholder bounding boxes where images will be stitched
#' @param overwrite overwrite if filename already used. Defaults to TRUE
#' @param verbose be verbose
#' @return largeGiottoImage object of stitched image
#' @export
stitchGiottoLargeImage = function(largeImage_list = NULL,
                                  gobject_list = NULL,
                                  largeImage_nameList = NULL,
                                  FOV_positions = NULL,
                                  FOV_xcol = NULL,
                                  FOV_ycol = NULL,
                                  FOV_inverty = FALSE,
                                  method = c('mosaic','merge'),
                                  round_positions = FALSE,
                                  filename = NULL,
                                  dataType = NULL,
                                  fileType = NULL,
                                  dryRun = TRUE,
                                  overwrite = FALSE,
                                  verbose = TRUE) {
  ## 0. Check params
  if(!is.null(gobject_list)) {
    # Set default largeImage_nameList
    if(is.null(largeImage_nameList)) {
      largeImage_nameList = rep("image", length(gobject_list))
    }
  }

  # Select method for stitching
  method = match.arg(method, choices = c('mosaic', 'merge'))

  # Check for filename, set default if not found
  if(is.null(filename)) {
    if(!is.null(gobject_list)) {
      save_dir = readGiottoInstructions(gobject_list[[1]], param = 'save_dir')
    } else {
      save_dir = path.expand('~')
    }
    filename = paste0(save_dir, '/stitch.tif')
  }

  # check filename
  if(file.exists(filename)) {
    if(verbose == TRUE) {
      if(overwrite == TRUE) cat('File at',filename,'exists.\n (overwrite == TRUE) Image will be overwritten')
      if(overwrite == FALSE) cat('File at',filename,'exists.\n (overwrite == FALSE) Image will not be overwritten')
    }
  }

  # Match dataType input if given
  dataTypeChoices = c('INT1U','INT2U','INT2S','INT4U','INT4S','FLT4S','FLT8S')
  if(!is.null(dataType)) dataType = match.arg(dataType, choices = dataTypeChoices)
  # Determine compatible dataType from first  giottoLargeImage


  ## 1. Get list of raster objects
  if(is.null(largeImage_list)) {
    if(!is.null(gobject_list)) {
      # For loop to grab giottoLargeImages
      largeImage_list = list()
      for(gobj_i in 1:length(gobject_list)) {
        largeImage_list[[gobj_i]] = getGiottoImage(gobject = gobject_list[[gobj_i]],
                                                   largeImage_name = largeImage_name[[gobj_i]])
      }
    } else {
      stop('giottoLargeImages must be given either as the giottoLargeImage itself or as a giotto object and the giottoLargeImage name')
    }
  }

  # Determine datatype from first giottoLargeImage
  if(is.null(dataType)) {
    dataType = find_terra_writeRaster_dataType(giottoLargeImage = largeImage_list[[1]])
  }

  # For loop to extract raster_objects
  raster_list = list()
  for(img_i in 1:length(largeImage_list)) {
    raster_list[[img_i]] = largeImage_list[[img_i]]@raster_object
  }

  ## 2. Apply FOV shifts (if given)
  if(!is.null(FOV_positions)) {

    # Check if there is an FOV position for every raster object
    if(nrow(FOV_positions) != length(raster_list)) {
      stop('If FOV_positions are given then there must be one set of values for every image being stitched')
    }

    if(FOV_inverty == TRUE) {
      FOV_positions['FOV_ycol'] = -FOV_positions['FOV_ycol']
    }

    # Shift the image extents as specified by POV_positions
    for(rast_i in 1:length(raster_list)) {
      raster_list[[rast_i]] = terra::shift(x = raster_list[[rast_i]],
                                           dx = FOV_positions[rast_i,FOV_xcol],
                                           dy = FOV_positions[rast_i,FOV_ycol])
    }
  }

  ## 3. Perform stitch
  # # Round final extent values (merge and mosaic may only work with integer extents)
  if(round_positions == TRUE) {
    if(verbose == TRUE) cat('round_positions == TRUE \n Image spatial positions will be rounded to integers. \n')
    for(rast_i in 1:length(raster_list)) {
      terra::ext(raster_list[[rast_i]]) = round(terra::ext(raster_list[[rast_i]]))
    }
  }

  if(dryRun == TRUE) {
    # Collect SpatExtents then convert to polygons
    imgBounds_list = list()
    for(rast_i in 1:length(raster_list)) {
      img_ext = terra::ext(raster_list[[rast_i]])
      img_bound_poly = terra::as.polygons(img_ext)
      img_bound_poly$FOV = rast_i
      # Add to imgBounds list
      imgBounds_list[[rast_i]] = img_bound_poly
    }
    imgBounds = do.call(rbind, imgBounds_list)
    terra::plot(imgBounds, 'FOV',
                type = 'classes',
                legend = TRUE,
                mar = c(3,3,2,2),
                plg = list(x = 'topright'))
    return(NULL)

  } else if (dryRun == FALSE) {
    # Create SpatRasterCollection
    rasterSRC = terra::src(raster_list)

    # stitch raster objects
    if(method == 'merge') {
      stitchImg = terra::merge(x = rasterSRC,
                               filename = filename,
                               overwrite = overwrite,
                               wopt = list(datatype = dataType))
    } else if(method == 'mosaic') {
      stitchImg = terra::mosaic(x = rasterSRC,
                                filename = filename,
                                overwrite = overwrite,
                                wopt = list(datatype = dataType,
                                            filetype = fileType))
    }
    stitch_gLargeImg = createGiottoLargeImage(raster_object = stitchImg,
                                              use_rast_ext = TRUE)
    return(stitch_gLargeImg)
  }

}



#' @title get_GiottoLargeImage
#' @name get_GiottoLargeImage
#' @description get a giottoLargeImage from a giottoObject
#' @param gobject giotto object
#' @param largeImage_name name of giottoLargeImage
#' @return a giottoLargeImage
#' @keywords internal
get_GiottoLargeImage = function(gobject = NULL,
                                largeImage_name = NULL) {

  if(is.null(gobject)) stop('The giotto object holding the giottoLargeImage needs to be provided \n')
  g_image_names = names(gobject@largeImages)
  if(is.null(g_image_names)) stop('No giottoLargeImages have been found \n')

  if(is.null(largeImage_name)) {
    largeImage_name = g_image_names[1]
  }

  if(!largeImage_name %in% g_image_names) stop(largeImage_name, ' was not found among the largeImage names. \n') #TODO See showImageNames()

  g_imageL = gobject@largeImages[[largeImage_name]]

  return(g_imageL)
}


#' @title cropGiottoLargeImage
#' @name cropGiottoLargeImage
#' @description crop a giottoLargeImage based on crop_extent argument or given values
#' @param gobject gobject holding the giottoLargeImage
#' @param largeImage_name name of giottoLargeImage within gobject
#' @param giottoLargeImage a giottoLargeImage
#' @param crop_name arbitrary name for cropped giottoLargeImage
#' @param crop_extent terra extent object used to crop the giottoLargeImage
#' @param xmax_crop crop xmax bound
#' @param xmin_crop crop xmin bound
#' @param ymax_crop crop ymax bound
#' @param ymin_crop crop ymin bound
#' @return a giottoLargeImage
#' @export
cropGiottoLargeImage = function(gobject = NULL,
                                largeImage_name = NULL,
                                giottoLargeImage = NULL,
                                crop_name = 'image',
                                crop_extent = NULL,
                                xmax_crop = NULL,
                                xmin_crop = NULL,
                                ymax_crop = NULL,
                                ymin_crop = NULL) {

  ## 0. Check inputs
  if(!is.null(crop_extent)) {
    if(!methods::is(crop_extent, 'SpatExtent')) stop('crop_extent argument only accepts terra extent objects. \n')
  }
  if(!is.null(giottoLargeImage)) {
    if(!methods::is(giottoLargeImage, 'giottoLargeImage')) stop('giottoLargeImage argument only accepts giottoLargeImage objects. \n')
  }

  ## 1. get giottoLargeImage if necessary
  if(is.null(giottoLargeImage)) {
    if(!is.null(gobject) && !is.null(largeImage_name)) {
      giottoLargeImage = getGiottoImage(gobject = gobject,
                                        largeImage_name = largeImage_name)
    } else {
      stop('either a giottoLargeImage or both the gobject and name of the giottoLargeImage must be given. \n')
    }
  }

  raster_object = giottoLargeImage@raster_object

  ## 2. Find crop extent
  crop_bounds = c(xmin_crop,xmax_crop,ymin_crop,ymax_crop)

  if(!is.null(crop_extent)) {
    raster_object = terra::crop(raster_object,
                                crop_extent,
                                snap = 'near')
  } else if(length(crop_bounds == 4)) {
    crop_extent = terra::ext(crop_bounds)

    raster_object = terra::crop(raster_object,
                                crop_extent,
                                snap = 'near')
  } else if(length(crop_bounds) > 1) {
    stop('All four crop bounds must be given.')
  }

  ## 3. Return a cropped giottoLargeImage
  giottoLargeImage@name = crop_name
  giottoLargeImage@raster_object = raster_object
  # The only things updated are the raster object itself and the name.
  # The overall_extent slot must NOT be updated since it records the original extent

  return(giottoLargeImage)
}


#' @title plot_giottoLargeImage
#' @name plot_giottoLargeImage
#' @description Plot a downsampled version of giottoLargeImage. Cropping can increase plot resolution of region of interest.
#' @param gobject giotto object
#' @param largeImage_name name of giottoLargeImage
#' @param giottoLargeImage giottoLargeImage object
#' @param crop_extent extent object to focus on specific region of image
#' @param xmax_crop xmax crop boundary
#' @param xmin_crop xmin crop boundary
#' @param ymax_crop ymax crop boundary
#' @param ymin_crop ymin crop boundary
#' @param max_intensity value to treat as maximum intensity in color scale
#' @return plot
#' @keywords internal
plot_giottoLargeImage = function(gobject = NULL,
                                 largeImage_name = NULL,
                                 giottoLargeImage = NULL,
                                 crop_extent = NULL,
                                 xmax_crop = NULL,
                                 xmin_crop = NULL,
                                 ymax_crop = NULL,
                                 ymin_crop = NULL,
                                 max_intensity = NULL) {

  # Get giottoLargeImage and check and perform crop if needed
  giottoLargeImage = cropGiottoLargeImage(gobject = gobject,
                                          largeImage_name = largeImage_name,
                                          giottoLargeImage = giottoLargeImage,
                                          crop_extent = crop_extent,
                                          xmax_crop = xmax_crop,
                                          xmin_crop = xmin_crop,
                                          ymax_crop = ymax_crop,
                                          ymin_crop = ymin_crop)

  raster_object = giottoLargeImage@raster_object

  # plot
  if(raster_object@ptr$rgb == FALSE) {
    terra::plotRGB(raster_object,
                   axes = TRUE,
                   r = 1,g = 1,b = 1,
                   stretch ='lin',
                   smooth = TRUE,
                   mar = c(3,5,1.5,1),
                   asp = 1)
  } else if(raster_object@ptr$rgb == TRUE) {

    # Determine likely image bitdepth
    if(is.null(max_intensity)) {
      bitDepth = ceiling(log(x = giottoLargeImage@max_intensity, base = 2))
      # Assign discovered bitdepth as max_intensity
      max_intensity = 2^bitDepth-1
    }

    terra::plotRGB(raster_object,
                   axes = TRUE,
                   r = 1,g = 2,b = 3,
                   scale = max_intensity,
                   smooth = TRUE,
                   mar = c(5,5,1,1),
                   asp = 1)
  }

}


#' @title convertGiottoLargeImageToMG
#' @name convertGiottoLargeImageToMG
#' @description convert a giottoLargeImage by downsampling into a normal magick based giottoImage
#' @param gobject gobject containing giottoLargeImage
#' @param largeImage_name name of giottoLargeImage
#' @param mg_name name to assign converted magick image based giottoImage. Defaults to name of giottoLargeImage
#' @param spat_unit spatial unit
#' @param spat_loc_name gobject spatial location name to map giottoImage to (optional)
#' @param crop_extent extent object to focus on specific region of image
#' @param xmax_crop assign crop boundary
#' @param xmin_crop assign crop boundary
#' @param ymax_crop assign crop boundary
#' @param ymin_crop assign crop boundary
#' @param resample_size maximum number of pixels to use when resampling
#' @param max_intensity value to treat as maximum intensity in color scale
#' @param return_gobject return as giotto object
#' @param verbose be verbose
#' @return a giotto object or an updated giotto image if return_gobject = F
#' @export
convertGiottoLargeImageToMG = function(gobject = NULL,
                                       largeImage_name = NULL,
                                       giottoLargeImage = NULL,
                                       mg_name = NULL,
                                       spat_unit = NULL,
                                       spat_loc_name = NULL,
                                       crop_extent = NULL,
                                       xmax_crop = NULL,
                                       xmin_crop = NULL,
                                       ymax_crop = NULL,
                                       ymin_crop = NULL,
                                       resample_size = 500000,
                                       max_intensity = NULL,
                                       return_gobject = TRUE,
                                       verbose = TRUE) {

  # Check params
  if(is.null(gobject)) {
    if(return_gobject == TRUE) stop('gobject must be given if return_gobject == TRUE')
    if(!is.null(spat_loc_name)) stop('if spatial location name is given then gobject containing it must also be given')
  }

  # Set spat_unit
  spat_unit = set_default_spat_unit(gobject = gobject,
                                    spat_unit = spat_unit)

  # Get giottoLargeImage and check and perform crop if needed
  giottoLargeImage = cropGiottoLargeImage(gobject = gobject,
                                          largeImage_name = largeImage_name,
                                          giottoLargeImage = giottoLargeImage,
                                          crop_extent = crop_extent,
                                          xmax_crop = xmax_crop,
                                          xmin_crop = xmin_crop,
                                          ymax_crop = ymax_crop,
                                          ymin_crop = ymin_crop)

  raster_object = giottoLargeImage@raster_object

  # Resample and then convert to Array
  rastSample = terra::spatSample(raster_object,
                                 size = resample_size, # Defines the rough maximum of pixels allowed when resampling
                                 method = 'regular',
                                 as.raster = TRUE)

  imArray = terra::as.array(rastSample)

  # Set max_intensity
  if(is.null(max_intensity)) {
    max_intensity = max(imArray)
  }

  # Read in array as magick image
  mImg = magick::image_read(imArray/max_intensity)

  # Set boundary adj values
    xmin_adj = xmax_adj = ymin_adj = ymax_adj = 0


  # magick object name
  if(is.null(mg_name)) {
    mg_name = giottoLargeImage@name
  }

  # Create giottoImage
  g_image = createGiottoImage(name = mg_name,
                              mg_object = mImg,
                              do_manual_adj = TRUE,
                              xmax_adj = xmax_adj,
                              xmin_adj = xmin_adj,
                              ymax_adj = ymax_adj,
                              ymin_adj = ymin_adj,
                              verbose = FALSE)

  # Set minimax
  if(is.null(spat_loc_name)) {
    current_ext = terra::ext(raster_object)
    g_image@minmax = c(current_ext$xmax,
                       current_ext$xmin,
                       current_ext$ymax,
                       current_ext$ymin)
  } else if(!is.null(spat_loc_name)) {
    spatial_locs = get_spatial_locations(gobject = gobject,
                                         spat_unit = spat_unit,
                                         spat_loc_name = spat_loc_name)
    x_range = range(spatial_locs$sdimx)
    y_range = range(spatial_locs$sdimy)
    g_image@minmax = c(x_range[2],
                       x_range[1],
                       y_range[2],
                       y_range[1])
  }

  names(g_image@minmax) = c('xmax_sloc','xmin_sloc','ymax_sloc','ymin_sloc')

  # Set scalefactor
  im_dims = magick::image_info(g_image@mg_object)
  x_scalefactor = im_dims[['width']] / dim(raster_object)[2]
  y_scalefactor = im_dims[['height']] / dim(raster_object)[1]


  g_image@scale_factor = c('x' = x_scalefactor, 'y' = y_scalefactor)
  g_image@resolution = 1/g_image@scale_factor

  if(return_gobject == TRUE) {
    if(verbose == TRUE) {
      if(mg_name %in% names(gobject@images)) cat('\n ', mg_name, ' has already been used, will be overwritten \n')
    }
    gobject@images[[mg_name]] = g_image
    return(gobject)
  } else if(return_gobject == FALSE) {
    return(g_image)
  }

}


#' @title find_terra_writeRaster_dataType
#' @name find_terra_writeRaster_dataType
#' @description find likely compatible datatype for given image characteristics. Values given in arguments take priority over those found from giottoLargeImage metadata
#' @param giottoLargeImage giottoLargeImage object to determine max_intensity, min_intensity, is_int settings from
#' @param quick_INTU_maxval Treat as maximum intensity to find compatible unsigned integer settings
#' @param max_intensity value given as image maximum intensity
#' @param min_intensity value given as image minimum intensity
#' @param is_int if image is integer (TRUE) or floating point (FALSE)
#' @param signed if image is signed (TRUE) or unsigned (TRUE)
#' @param bitDepth image bitDepth
#' @param verbose be verbose
#' @keywords internal
#' @return datatype for terra writeRaster function
find_terra_writeRaster_dataType = function(giottoLargeImage = NULL,
                                           quick_INTS_maxval = NULL,
                                           max_intensity = NULL,
                                           min_intensity = NULL,
                                           is_int = NULL,
                                           signed = NULL,
                                           bitDepth = NULL,
                                           verbose = TRUE) {

  # 1. Get any missing metadata from giottoLargeImage object if given
  if(!is.null(giottoLargeImage)) {
    if(is.null(max_intensity)) max_intensity = giottoLargeImage@max_intensity
    if(is.null(min_intensity)) min_intensity = giottoLargeImage@min_intensity
    if(is.null(is_int)) is_int = giottoLargeImage@is_int
  }


  if(is.null(quick_INTS_maxval)) {
    if(length(c(max_intensity, min_intensity)) < 2) stop('Not enough metadata is given')

    # Determine if negative values are present

    if(min_intensity < 0) {
      signed = TRUE
    }
  }

  # Set defaults if data still missing
  if(is.null(is_int)) is_int = TRUE
  if(is.null(signed)) signed = FALSE

  ## 2. Determine likely compatible datatype
  dataTypeVerbose = data.frame(bitDepth = c(8,16,16,32,32,32,64),
                               signed = c(FALSE,FALSE,TRUE,FALSE,TRUE,TRUE,TRUE),
                               integer = c(TRUE,TRUE,TRUE,TRUE,TRUE,FALSE,FALSE),
                               dataTypeChoices = c('INT1U','INT2U','INT2S','INT4U','INT4S','FLT4S','FLT8S'),
                               dataTypeVerbose = c('8bit unsigned integer','16bit unsigned integer','16bit signed integer',
                                                   '32bit unsigned integer','32bit signed integer','32bit signed floating point',
                                                   '64bit signed floating point'))


  ## Find Compatible Bitdepth
  # If quick_INTS_maxval argument is set, will be treated as the highest needed (unsigned preferred) bitdepth datatype
  if(is.null(quick_INTS_maxval)) {

    if(max_intensity > 0) {

      max_intensity = max_intensity + 1 # Accounts for 0 occupying 1 of the available values.

      if(signed == FALSE) {
        bitDepth = ceiling(log(x = max_intensity, base = 2))
      } else if(signed == TRUE) {
        intensityMinMax = c(min_intensity, max_intensity)
        intensityMinMax = abs(intensityMinMax)
        bitDepthMinMax = ceiling(log(x = intensityMinMax, base = 2))
        bitDepth = max(bitDepthMinMax) + 1
      }
    } else {
      stop('There are no positive image intensities. \n Manual datatype assignment needed \n')
    }

  } else if(!is.null(quick_INTS_maxval)) {
    if(verbose == TRUE) cat('Selecting compatible datatype for given maximum value \n')
    bitDepth = ceiling(log(x = quick_INTS_maxval, base = 2))
  }

  if(bitDepth > 32 && bitDepth <= 128) {
    bitDepth = 32
    is_int = FALSE
    signed = TRUE
  } else if(bitDepth > 128) {
    bitDepth = 64
    is_int = FALSE
    signed = TRUE
  }

  dataType = NULL
  # Determine datatype settings
  if(is_int == TRUE) {
    if(signed == FALSE) {
      if(bitDepth <= 8) {
        dataType = 'INT1U'
      } else if(bitDepth <= 16) {
        dataType = 'INT2U'
      } else if(bitDepth <= 32) {
        dataType = 'INT4U'
      }
    } else if(signed == TRUE) {
      if(bitDepth <= 16) {
        dataType = 'INT2S'
      } else if(bitDepth <= 32) {
        dataType = 'INT4S'
      }
    }
  } else if(is_int == FALSE) {
    if(bitDepth <= 32) {
      dataType = 'FLT4S'
    } else if(bitDepth == 64) {
      dataType = 'FLT8S'
      # These are very large numbers. Can't actually tell the difference from FLT4S (less than 2^128) to FLT8S unless you add roughly 10^25 to it.
      # This necessary minimum change would be 10^22, but the log used when determining bitDepth further increases the needed difference.
      # Manual assignment of dataType could be more reliable than automatic assignment for these very large values.
    }
  }
  return(dataType)
}



#' @title writeGiottoLargeImage
#' @name writeGiottoLargeImage
#' @description write original resolution to file. Filetype extension should be included in filename argument.
#' @param giottoLargeImage giottoLargeImage object
#' @param gobject giotto object
#' @param largeImage_name name of giottoLargeImage
#' @param max_intensity (optional) image max intensity value from which dataType can be automatically determined
#' @param filename path to write the image to
#' @param dataType (optional) values for dataType are "INT1U", "INT2U", "INT2S", "INT4U", "INT4S", "FLT4S", "FLT8S". The first three letters indicate whether the dataType is integer (whole numbers) of a real number (decimal numbers), the fourth character indicates the number of bytes used (allowing for large numbers and/or more precision), and the "S" or "U" indicate whether the values are signed (both negative and positive) or unsigned (positive values only).
#' @param overwrite Overwrite if filename is already existing
#' @param verbose be verbose
#' @export
writeGiottoLargeImage = function(giottoLargeImage = NULL,
                                 gobject = NULL,
                                 largeImage_name = NULL,
                                 max_intensity = NULL,
                                 filename = NULL,
                                 dataType = NULL,
                                 overwrite = FALSE,
                                 verbose = TRUE) {

  # 0. Check params
  if(!is.null(giottoLargeImage)) {
    if(!methods::is(giottoLargeImage, 'giottoLargeImage')) stop('giottoLargeImage argument only accepts giottoLargeImage objects. \n')
  }
  if(!is.null(max_intensity)) {
    if(!is.numeric(max_intensity)) stop('max_intensity must be a numeric \n')
  }
  if(!is.null(filename)) {
    if(!is.character(filename)) stop('filename must be given as character \n')
    # check filename
    if(file.exists(filename)) {
      if(verbose == TRUE) {
        if(overwrite == TRUE) cat('File at',filename,'exists.\n (overwrite == TRUE) Image will be overwritten')
        if(overwrite == FALSE) cat('File at',filename,'exists.\n (overwrite == FALSE) Image will not be overwritten')
      }
    }
  }

  if(is.null(filename)) stop('Please enter a filename to save the image as. \n')

  filename = path.expand(filename)

  # Match dataType input if given
  dataTypeChoices = c('INT1U','INT2U','INT2S','INT4U','INT4S','FLT4S','FLT8S')
  if(!is.null(dataType)) dataType = match.arg(dataType, choices = dataTypeChoices)


  ## 1. get giottoLargeImage if necessary
  if(is.null(giottoLargeImage)) {
    if(!is.null(gobject) && !is.null(largeImage_name)) {
      giottoLargeImage = getGiottoImage(gobject = gobject,
                                        largeImage_name = largeImage_name)
    } else {
      stop('either a giottoLargeImage or both the gobject and name of the giottoLargeImage must be given. \n')
    }
  }

  raster_object = giottoLargeImage@raster_object

  ## 2. Get likely compatible dataType
  if(is.null(dataType)) {
    dataType = find_terra_writeRaster_dataType(giottoLargeImage = giottoLargeImage,
                                               quick_INTS_maxval = max_intensity)
  }


  ## 3. Write to disk
  if(verbose == TRUE) cat(paste0('Writing image to disk as ', dataType))
  terra::writeRaster(x = raster_object,
                     filename = filename,
                     datatype = dataType,
                     overwrite = overwrite)
}



#' @title updateGiottoLargeImage
#' @name updateGiottoLargeImage
#' @description Updates the boundaries of a giottoLargeImage attached to a giotto object
#' @param gobject giotto object containing giottoLargeImage
#' @param largeImage_name name of giottoLargeImage
#' @param giottoLargeImage giottoLargeImage object
#' @param xmax_adj adjustment of the maximum x-value to align the image
#' @param xmin_adj adjustment of the minimum x-value to align the image
#' @param ymax_adj adjustment of the maximum y-value to align the image
#' @param ymin_adj adjustment of the minimum y-value to align the image
#' @param x_shift shift entire image in positive x direction
#' @param y_shift shift entire image in positive y direction
#' @param scale_factor set scale_x and scale_y at the same time
#' @param scale_x independently scale x axis image mapping from origin
#' @param scale_y independently scale y axis image mapping from origin
#' @param order perform fine adjustments (adjustments and shifts) or scaling first
#' @param xmin_bound set image xmin boundary. Overrides minmax values as spatial anchor.
#' @param xmax_bound set image xmax boundary. Overrides minmax values as spatial anchor.
#' @param ymin_bound set image ymin boundary. Overrides minmax values as spatial anchor.
#' @param ymax_bound set image ymax boundary. Overrides minmax values as spatial anchor.
#' @param return_gobject return a giotto object if TRUE, a giottoLargeImage if FALSE
#' @param verbose be verbose
#' @return a giotto object or an updated giottoLargeImage if return_gobject = FALSE
#' @export
updateGiottoLargeImage = function(gobject = NULL,
                                  largeImage_name = NULL,
                                  giottoLargeImage = NULL,
                                  xmax_adj = 0,
                                  xmin_adj = 0,
                                  ymax_adj = 0,
                                  ymin_adj = 0,
                                  x_shift = 0,
                                  y_shift = 0,
                                  scale_factor = NULL,
                                  scale_x = 1,
                                  scale_y = 1,
                                  order = c('first_adj', 'first_scale'),
                                  xmin_set = NULL,
                                  xmax_set = NULL,
                                  ymin_set = NULL,
                                  ymax_set = NULL,
                                  return_gobject = TRUE,
                                  verbose = TRUE) {


  # 0. Check params
  # Check input image
  if(is.null(gobject)) {
    if(is.null(giottoLargeImage)) stop('Image to be updated must be given as gobject AND largeImage_name OR giottoLargeImage argument(s) \n')
    if(verbose == TRUE) cat('gobject argument not given \n return_gobject set to FALSE \n')
    return_gobject = FALSE
  }
  if(is.null(giottoLargeImage) && is.null(largeImage_name)) stop('The name of the giottoLargeImage that will be updated needs to be provided \n')

  if(!is.null(giottoLargeImage)) {
    if(!methods::is(giottoLargeImage, 'giottoLargeImage')) stop('giottoLargeImage argument only accpts giottoLargeImage objects \n')
    if(verbose == TRUE && !is.null(gobject)) cat('giottoLargeImage argument is given and will take priority \n return_gobject set to FALSE \n')
    return_gobject = FALSE
  }

  # Check scalefactors
  if(!is.null(scale_factor)) scale_x = scale_y = scale_factor

  # Check spatial anchor values
  spatAnchor = c(xmin_set,
                 xmax_set,
                 ymin_set,
                 ymax_set)
  if(length(spatAnchor) < 4 && length(spatAnchor) > 0) stop('If set arguments are being used, all four must be given \n')
  if(!is.null(spatAnchor)) {
    if(xmax_set < xmin_set) stop('xmax_set must be greater than xmin_set \n')
    if(ymax_set < ymin_set) stop('ymax_set must be greater than ymin_set \n')
  }

  # Find order of adjust and scaling
  order = match.arg(order, choices = c('first_adj','first_scale'))

  # 1. get giottoImage if necessary
  if(is.null(giottoLargeImage)) {
    if(!is.null(gobject) && !is.null(largeImage_name)) {
      g_imageL = getGiottoImage(gobject = gobject,
                               largeImage_name = largeImage_name)
    } else {
      stop('either a giottoLargeImage or both the gobject and name of the giottoLargeImage must be given. \n')
    }
  }


  # 2. Find minmax spatial anchor values if set values not supplied
  if(is.null(spatAnchor)) {
    spatAnchor = terra::ext(x = g_imageL@raster_object)[1:4] #(xmin, xmax, ymin, ymax)
    names(spatAnchor) = NULL
  }

  # Perform scale if first_scale
  if(order == 'first_scale') {
    spatAnchor = spatAnchor * c(scale_x, scale_x, scale_y, scale_y)
  }

  # 3. Prepare adjustment values
  # Apply shifts
  xmin_adj = xmin_adj - x_shift
  xmax_adj = xmax_adj + x_shift
  ymin_adj = ymin_adj - y_shift
  ymax_adj = ymax_adj + y_shift

  # Find final bounds
  xmin_final = spatAnchor[1] - xmin_adj
  xmax_final = spatAnchor[2] + xmax_adj
  ymin_final = spatAnchor[3] - ymin_adj
  ymax_final = spatAnchor[4] + ymax_adj

  # Perform scale if first_adj
  if(order == 'first_adj') {
    xmin_final = xmin_final * scale_x
    xmax_final = xmax_final * scale_x
    ymin_final = ymin_final * scale_y
    ymax_final = ymax_final * scale_y
  }


  # 4. Update the boundaries
  if(return_gobject == FALSE) g_imageL@raster_object = terra::deepcopy(g_imageL@raster_object)
  terra::ext(g_imageL@raster_object) = c(xmin_final,
                                         xmax_final,
                                         ymin_final,
                                         ymax_final)


  #5. Update the scalefactors for x and y
  g_imageL@resolution = terra::res(g_imageL@raster_object) #(x,y)
  names(g_imageL@resolution) = c('x','y')
  g_imageL@scale_factor = (1/g_imageL@resolution)


  if(return_gobject == TRUE) {
    gobject@largeImages[[largeImage_name]] = g_imageL
    return(gobject)
  } else {
    return(g_imageL)
  }

}


#' @title addGiottoLargeImage
#' @name addGiottoLargeImage
#' @description Adds giotto image objects to your giotto object
#' @param gobject giotto object
#' @param largeImages list of giottoLargeImage objects
#' @param spat_loc_name provide spatial location slot in Giotto to align images. (optional)
#' @param scale_factor provide scale of image pixel dimensions relative to spatial coordinates.
#' @param negative_y Map image to negative y spatial values if TRUE during automatic alignment. Meaning that origin is in upper left instead of lower left.
#' @return an updated Giotto object with access to the list of images
#' @export
addGiottoLargeImage = function(gobject = NULL,
                               largeImages = NULL,
                               spat_loc_name = NULL,
                               scale_factor = NULL,
                               negative_y = TRUE,
                               verbose = TRUE) {


  # 0. check params
  if(is.null(gobject)) stop('The giotto object that will be updated needs to be provided')

  if(is.null(largeImages)) stop('The giotto large image(s) that will be added needs to be provided')

  ext_scale_factor = FALSE
  if(!is.null(scale_factor)) {

    if(!is.numeric(scale_factor)) stop ('Given scale_factor(s) must be numeric')

    if((length(scale_factor) == length(largeImages)) || length(scale_factor) == 1) {
      cat('scale_factor(s) external to giottoImage have been given and will be used')
      ext_scale_factor = TRUE
    } else {
      stop('if scale_factor is given, it must be a numeric with either a single value or as many values as there are largeImages are provided')
    }
  }

  # 1. expand scale_factors
  if(ext_scale_factor == TRUE) {
    if(length(scale_factor == 1)) {
      scale_factor = rep(scale_factor, length(largeImages))
    }
  }


  # 2. Add image with for loop
  for(image_i in 1:length(largeImages)) {

    im = largeImages[[image_i]]

    if(methods::is(im, 'giottoLargeImage')) {
      im_name = im@name

      all_im_names = names(gobject@largeImages)

      if(im_name %in% all_im_names) {
        cat('\n ', im_name, ' has already been used, will be overwritten \n')
      }

      # Deep copy the raster_object
      im@raster_object = terra::deepcopy(im@raster_object)

      # # 3. Update boundaries if not already done during createGiottoImage() due to lack of spatlocs and gobject
      # if(sum(im@boundaries == c(0,0,0,0)) == 4 && sum(im@minmax == c(10,0,10,0)) == 4) {
      #   if(!is.null(spat_loc_name)) { # A check for the first available spatloc was already done
      #     spatlocs = get_spatial_locations(gobject = gobject,
      #                                      spat_loc_name = spat_loc_name)
      #
      #     #Find spatial minmax values
      #     xmin_sloc = min(spatlocs$sdimx)
      #     xmax_sloc = max(spatlocs$sdimx)
      #     ymin_sloc = min(spatlocs$sdimy)
      #     ymax_sloc = max(spatlocs$sdimy)
      #
      #     #Find adjustment values
      #     img_minmax = get_img_minmax(mg_img = im@mg_object,
      #                                 negative_y = negative_y)
      #     if(ext_scale_factor == TRUE) {
      #       adj_values = get_adj_rescale_img(img_minmax = img_minmax,
      #                                        spatial_locs = spatlocs,
      #                                        scale_factor = scale_factor[[image_i]])
      #     } else if (ext_scale_factor == FALSE) {
      #       adj_values = get_adj_rescale_img(img_minmax = img_minmax,
      #                                        spatial_locs = spatlocs,
      #                                        scale_factor = im@scale_factor[[spat_loc_name]])
      #     }
      #
      #     #Add minmax values to giottoImage@minmax
      #     im@minmax = c('xmax_sloc' = xmax_sloc,
      #                   'xmin_sloc' = xmin_sloc,
      #                   'ymax_sloc' = ymax_sloc,
      #                   'ymin_sloc' = ymin_sloc)
      #
      #     #Add adjustment values to giottoImage@boundaries
      #     im@boundaries = c('xmax_adj' = as.numeric(adj_values[['xmax_adj_orig']]),
      #                       'xmin_adj' = as.numeric(adj_values[['xmin_adj_orig']]),
      #                       'ymax_adj' = as.numeric(adj_values[['ymax_adj_orig']]),
      #                       'ymin_adj' = as.numeric(adj_values[['ymin_adj_orig']]))
      #
      #     # Inherit external scaling factors if given
      #     if(ext_scale_factor == TRUE) {
      #       im@scale_factor[[spat_loc_name]] = scale_factor[[image_i]]
      #       im@resolution[[spat_loc_name]] = 1/(scale_factor[[image_i]])
      #     }
      #     ## Externally given scale_factors will only be written in/used if boundary adj values are not pre-existing
      #   }
      #
      # }

      # 4. Add giottoImage to gobject
      gobject@largeImages[[im_name]] = im

    } else {
      warning('image [',image_i,'] is not a giotto image object')
    }
  }

  return(gobject)

}



# Image Tools ####


#' @name get_img_minmax
#' @title get_img_minmax
#' @param mg_img magick object
#' @param negative_y Map image to negative y spatial values if TRUE during automatic alignment. Meaning that origin is in upper left instead of lower left.
#' @keywords internal
get_img_minmax = function(mg_img,
                          negative_y = TRUE) {
  #Get magick object dimensions. xmin and ymax assumed to be 0.
  info = magick::image_info(mg_img)
  img_xmax = info$width     #width
  img_xmin = 0              #x origin
  if(negative_y == TRUE) {
    img_ymax = 0              #y origin
    img_ymin = -(info$height) #height
  } else if(negative_y == FALSE) {
    img_ymax = info$height
    img_ymin = 0
  }


  return(list('img_xmax' = img_xmax,
              'img_xmin' = img_xmin,
              'img_ymax' = img_ymax,
              'img_ymin' = img_ymin))
}


#' @name get_adj_rescale_img
#' @title get_adj_rescale_img
#' @keywords internal
get_adj_rescale_img = function(img_minmax,
                               spatial_locs,
                               scale_factor = 1) {

  # Expand scale_factor if needed
  if(length(scale_factor) == 1) {
    scale_factor = c(x = scale_factor, y = scale_factor)
  }

  # Spatial minmax
  my_xmin = min(spatial_locs$sdimx)
  my_xmax = max(spatial_locs$sdimx)
  my_ymin = min(spatial_locs$sdimy)
  my_ymax = max(spatial_locs$sdimy)

  # Find scaled image adjustments based on scaled spatlocs
  xmin_adj_scaled = (my_xmin*scale_factor[['x']]) - (img_minmax$img_xmin)
  xmin_adj_orig = xmin_adj_scaled/scale_factor[['x']]

  xmax_adj_scaled = (img_minmax$img_xmax) - (my_xmax*scale_factor[['x']])
  xmax_adj_orig = xmax_adj_scaled/scale_factor[['x']]

  ymin_adj_scaled = (my_ymin*scale_factor[['y']]) - (img_minmax$img_ymin)
  ymin_adj_orig = ymin_adj_scaled/scale_factor[['y']]

  ymax_adj_scaled = (img_minmax$img_ymax) - (my_ymax*scale_factor[['y']])
  ymax_adj_orig = ymax_adj_scaled/scale_factor[['y']]

  # return scaled adjustments
  return(c('xmin_adj_orig' = xmin_adj_orig,
           'xmax_adj_orig' = xmax_adj_orig,
           'ymin_adj_orig' = ymin_adj_orig,
           'ymax_adj_orig' = ymax_adj_orig))

}



#' @title plotGiottoImage
#' @name plotGiottoImage
#' @description plot a giottoImage or giottoLargeImage
#' @param gobject gobject containing giottoImage or giottoLargeImage
#' @param image_name name of giottoImage
#' @param largeImage_name name of giottoLargeImage
#' @param giottoImage giottoImage object
#' @param giottoLargeImage giottoLargeImage object
#' @param largeImage_crop_params_list list of parameters for focusing on a specified region of a giottoLargeImage for potentially better plotting resolution
#' @param largeImage_max_intensity assign override value to treat as maximum intensity in color scale when plotting giottoLargeImage
#' @return a plot
#' @export
plotGiottoImage = function(gobject = NULL,
                           image_name = NULL,
                           largeImage_name = NULL,
                           giottoImage = NULL,
                           giottoLargeImage = NULL,
                           largeImage_crop_params_list = NULL,
                           largeImage_max_intensity = NULL) {

  # Check params
  if(!is.null(gobject)) {
    if(!is.null(image_name) && !is.null(largeImage_name)) stop('Only one of a giottoImage or a giottoLargeImage can be plotted at the same time. \n')
  }
  if(!is.null(giottoImage) && !is.null(giottoLargeImage)) stop('Only one of a giottoImage or a giottoLargeImage can be plotted at the same time. \n')

  # Select plotting function
  if(!is.null(image_name) || !is.null(giottoImage)) {
    plot_giottoImage_MG(gobject = gobject,
                        image_name = image_name,
                        giottoImage = giottoImage)
  }
  if(!is.null(largeImage_name) || !is.null(giottoLargeImage)) {
    plot_giottoLargeImage(gobject = gobject,
                          largeImage_name = largeImage_name,
                          giottoLargeImage = giottoLargeImage,
                          crop_extent = largeImage_crop_params_list$crop_extent,
                          xmax_crop = largeImage_crop_params_list$xmax_crop,
                          xmin_crop = largeImage_crop_params_list$xmin_crop,
                          ymax_crop = largeImage_crop_params_list$ymax_crop,
                          ymin_crop = largeImage_crop_params_list$ymin_crop,
                          max_intensity = largeImage_max_intensity)
  }
}


#' @title getGiottoImage
#' @name getGiottoImage
#' @description get giottoImage or giottoLargeImage from gobject
#' @param gobject giotto object
#' @param image_name name of a giottoImage \code{\link{showGiottoImageNames}}
#' @param largeImage_name name of a giottoLargeImage
#' @return a giottoImage or giottoLargeImage
#' @export
getGiottoImage = function(gobject = NULL,
                          image_name= NULL,
                          largeImage_name = NULL) {

  # Check params
  if(!is.null(image_name) && !is.null(largeImage_name)) stop('Get only one of a giottoImage or a giottoLargeImage at a time. \n')
  if(is.null(image_name) && is.null(largeImage_name)) image_name = 'image'

  # Select get function
  if(!is.null(image_name)) {
    g_img = get_GiottoImage_MG(gobject = gobject,
                               image_name = image_name)
  }
  if(!is.null(largeImage_name)) {
    g_img = get_GiottoLargeImage(gobject = gobject,
                                 largeImage_name = largeImage_name)
  }
  return(g_img)
}


#' @title addGiottoImage
#' @name addGiottoImage
#' @description Adds lists of giottoImages and giottoLargeImages to gobjects
#' @param gobject gobject to add images objects to
#' @param images list of giottoImages to add
#' @param largeImages list of giottoLargeImages to add
#' @param spat_loc_name provide spatial location slot in Giotto to align giottoImages. Defaults to first one
#' @param scale_factor provide scale of image pixel dimensions relative to spatial coordinates.
#' @param negative_y Map image to negative y spatial values if TRUE during automatic alignment. Meaning that origin is in upper left instead of lower left.
#' @return an updated Giotto object with access to the list of images
#' @export
addGiottoImage = function(gobject = NULL,
                          images = NULL,
                          largeImages = NULL,
                          spat_loc_name = NULL,
                          scale_factor = NULL,
                          negative_y = TRUE) {

  if(!is.null(images) && !is.null(largeImages)) stop('Can only add one type of image to giotto object at a time')
  if(!is.null(images)) {
    addGiottoImageMG(gobject = gobject,
                     images = images,
                     spat_loc_name = spat_loc_name,
                     scale_factor = scale_factor,
                     negative_y = negative_y)
  } else if(!is.null(largeImages)) {
    addGiottoLargeImage(gobject = gobject,
                        largeImages = largeImages,
                        spat_loc_name = spat_loc_name,
                        scale_factor = scale_factor,
                        negative_y = negative_y)
  }
}



#' @title updateGiottoImage
#' @name updateGiottoImage
#' @description Updates the boundaries of a giottoImage or giottoLargeImage attached to a giotto object
#' @param gobject gobject containing giottoImage or giottoLargeImage
#' @param image_name name of giottoImage
#' @param largeImage_name name of giottoLargeImage
#' @param xmax_adj adjustment of the maximum x-value to align the image
#' @param xmin_adj adjustment of the minimum x-value to align the image
#' @param ymax_adj adjustment of the maximum y-value to align the image
#' @param ymin_adj adjustment of the minimum y-value to align the image
#' @param x_shift shift entire image in positive x direction
#' @param y_shift shift entire image in positive y direction
#' @param scale_x scale x axis image mapping from origin
#' @param scale_y scale y axis image mapping from origin
#' @param order perform fine adjustments (adjustments and shifts) or scaling first
#' @param xmin_set set image xmin boundary. Overrides minmax values as spatial anchor.
#' @param xmax_set set image xmax boundary. Overrides minmax values as spatial anchor.
#' @param ymin_set set image ymin boundary. Overrides minmax values as spatial anchor.
#' @param ymax_set set image ymax boundary. Overrides minmax values as spatial anchor.
#' @param return_gobject return a giotto object
#' @return a giotto object or an updated giotto image if return_gobject = F
#' @export
updateGiottoImage = function(gobject = NULL,
                             image_name = NULL,
                             largeImage_name = NULL,
                             xmax_adj = 0,
                             xmin_adj = 0,
                             ymax_adj = 0,
                             ymin_adj = 0,
                             x_shift = 0,
                             y_shift = 0,
                             scale_x = 1,
                             scale_y = 1,
                             order = c('first_adj', 'first_scale'),
                             xmax_set = NULL,
                             xmin_set = NULL,
                             ymax_set = NULL,
                             ymin_set = NULL,
                             return_gobject = TRUE) {


  # 0. Check params
  if(is.null(gobject)) stop('The giotto object that will be updated needs to be provided \n')
  if(is.null(image_name) && is.null(largeImage_name)) stop('The name of the giotto image that will be updated needs to be provided \n')

  order = match.arg(order, choices = c('first_adj','first_scale'))

  if(!is.null(image_name) && !is.null(largeImage_name)) stop('Adjust only giottoImage OR giottoLargeImage at any one time. \n')


  # 2. Select adjustment function
  if(!is.null(image_name)) {
    out = updateGiottoImageMG(gobject = gobject,
                              image_name = image_name,
                              xmax_adj = xmax_adj,
                              xmin_adj = xmin_adj,
                              ymax_adj = ymax_adj,
                              ymin_adj = ymin_adj,
                              x_shift = x_shift,
                              y_shift = y_shift,
                              scale_x = scale_x,
                              scale_y = scale_y,
                              order = order,
                              xmax_set = xmax_set,
                              xmin_set = xmin_set,
                              ymax_set = ymax_set,
                              ymin_set = ymin_set,
                              return_gobject = return_gobject)
  } else if(!is.null(largeImage_name)) {
    out = updateGiottoLargeImage(gobject = gobject,
                                 largeImage_name = largeImage_name,
                                 xmax_adj = xmax_adj,
                                 xmin_adj = xmin_adj,
                                 ymax_adj = ymax_adj,
                                 ymin_adj = ymin_adj,
                                 x_shift = x_shift,
                                 y_shift = y_shift,
                                 scale_x = scale_x,
                                 scale_y = scale_y,
                                 order = order,
                                 xmax_set = xmax_set,
                                 xmin_set = xmin_set,
                                 ymax_set = ymax_set,
                                 ymin_set = ymin_set,
                                 return_gobject = return_gobject)
  }
  return(out)
}



#' @title reconnectImage
#' @name reconnectImage
#' @description reconnect dead image pointers using filepaths
#' @param gobject giotto object
#' @param image_name name of image to be reconnected
#' @param largeImage_name name of largeImage to be reconnected
#' @return a giotto object with updated image pointer
reconnectImage = function(gobject,
                          image_name = NULL,
                          largeImage_name = NULL,
                          path = NULL) {
  
  # Check params
  if(!is.null(image_name) && !is.null(largeImage_name)) {
    stop('Only one image type can be reconnected at once')
  }
  
  # Load new image pointer and update
  if(!is.null(image_name)) {
    mg_object = magick::image_read(path = path)
    gobject@images[[image_name]]@mg_object = mg_object
  }
  if(!is.null(largeImage_name)) {
    # terra::rast() - Needs more work. Need a way to save the extents because those are part of the raster obj.
    stop('largeImage reconnection still not supported')
  }
  return(gobject)
}







