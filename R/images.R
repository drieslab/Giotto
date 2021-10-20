

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
    OS_platform = "ANY"
  ),

  prototype = list(
    name = NULL,
    mg_object = NULL,
    minmax = NULL,
    boundaries = NULL,
    scale_factor = NULL,
    resolution = NULL,
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
    
    cat("Scale factor(s) for spatlocs are: \n")
    for(name in names(object@scale_factor)) {
      cat(' ',name,': ',object@scale_factor[[name]][1], "\n")
    }
    
    cat("Resolution(s) for spatlocs are: \n")
    for(name in names(object@resolution)) {
      cat(' ',name,': ',object@resolution[[name]][1], "\n")
    }

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
    OS_platform = "ANY"
  ),
  
  prototype = list(
    name = NULL,
    raster_object = NULL,
    overall_extent = NULL,
    scale_factor = NULL,
    resolution = NULL,
    max_intensity = NULL,
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
    
    cat('\n Maximum intensity is: ', object@max_intensity, " \n")
    
  }
)


# giottoImage creation ####


#' @title createGiottoImage
#' @name createGiottoImage
#' @description Creates a giotto image that can be added to a Giotto object and/or used to add an image to the spatial plotting functions
#' @param gobject giotto object
#' @param spatial_locs spatial locations (alternative if gobject = NULL)
#' @param spat_loc_name name of spatial locations within gobject
#' @param mg_object magick image object
#' @param name name for the image
#' @param image_transformations vector of sequential image transformations
#' @param do_manual_adj flag to use manual adj values instead of automatic alignment when given a gobject or spatlocs
#' @param xmax_adj adjustment of the maximum x-value to align the image
#' @param xmin_adj adjustment of the minimum x-value to align the image
#' @param ymax_adj adjustment of the maximum y-value to align the image
#' @param ymin_adj adjustment of the minimum y-value to align the image
#' @param scale_factor scaling of image dimensions relative to spatial coordinates
#' @param verbose be verbose
#' @details image_transformations: transformation options from magick library
#' [\strong{flip_x_axis}] flip x-axis (\code{\link[magick]{image_flop}})
#' [\strong{flip_y_axis}] flip y-axis (\code{\link[magick]{image_flip}})
#' Example: image_transformations = c(flip_x_axis, flip_y_axis); first flip x-axis and then y-axis
#' @return a giottoImage object
#' @export
createGiottoImage = function(gobject = NULL,
                             spatial_locs = NULL,
                             spat_loc_name = NULL,
                             mg_object,
                             name = 'image',
                             image_transformations = NULL,
                             do_manual_adj = FALSE,
                             xmax_adj = 0,
                             xmin_adj = 0,
                             ymax_adj = 0,
                             ymin_adj = 0,
                             scale_factor = 1,
                             verbose = TRUE) {


  # create minimum giotto
  g_image = giottoImage(name = name,
                        mg_object = NULL,
                        minmax = NULL,
                        boundaries = NULL,
                        scale_factor = NULL,
                        resolution = NULL,
                        OS_platform = .Platform[['OS.type']])


  ## 1.a. check magick image object
  if(!methods::is(mg_object, 'magick-image')) {
    if(file.exists(mg_object)) {
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
  if(do_manual_adj == TRUE) cat('do_manual_adj == TRUE \n','Boundaries will be adjusted by given values.\n')
  #If spatlocs or gobject supplied, minmax values will always be generated
  #If do_manual_adj == TRUE, bypass followup automatic boundary value generation
  if(!is.null(gobject)) {
    
    if(!is.null(gobject@spatial_locs)) {
      spat_loc_name = names(gobject@spatial_locs)[[1]]
    }

    spatlocs = get_spatial_locations(gobject = gobject,
                                     spat_loc_name = spat_loc_name)

    # spatlocs = gobject@spatial_locs[['raw']]

    my_xmin = min(spatlocs$sdimx)
    my_xmax = max(spatlocs$sdimx)
    my_ymin = min(spatlocs$sdimy)
    my_ymax = max(spatlocs$sdimy)
    
    if(do_manual_adj == FALSE) {
      #find adjustment values
      img_minmax = get_img_minmax(mg_img = mg_object)
      adj_values = get_adj_rescale_img(img_minmax = img_minmax,
                                       spatial_locs = spatlocs,
                                       scale_factor = scale_factor)
      
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

    my_xmin = min(spatlocs$sdimx)
    my_xmax = max(spatlocs$sdimx)
    my_ymin = min(spatlocs$sdimy)
    my_ymax = max(spatlocs$sdimy)
    
    if(do_manual_adj == FALSE) {
      #find adjustment values
      img_minmax = get_img_minmax(mg_img = mg_object)
      adj_values = get_adj_rescale_img(img_minmax = img_minmax,
                                       spatial_locs = spatlocs,
                                       scale_factor = scale_factor)
      
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

    my_xmin = 0; my_xmax = 10; my_ymin = 0; my_ymax = 10

  }
  #TODO Add spatloc specific slots
  # minmax and boundary values for return
  g_image@minmax = c('xmax_sloc' = my_xmax,
                     'xmin_sloc' = my_xmin,
                     'ymax_sloc' = my_ymax,
                     'ymin_sloc' = my_ymin)

  g_image@boundaries = c('xmax_adj' = xmax_adj,
                         'xmin_adj' = xmin_adj,
                         'ymax_adj' = ymax_adj,
                         'ymin_adj' = ymin_adj)
  
  # scale factor and resolution values for return
  
  if(!is.null(spat_loc_name)) {
    g_image@scale_factor[[spat_loc_name]] = scale_factor
    g_image@resolution[[spat_loc_name]] = 1/scale_factor
  } else {
    g_image@scale_factor$raw = scale_factor
    g_image@resolution$raw = 1/scale_factor    
  }


  # image object
  return(g_image)
}



# giottoLargeImage creation ####


#' @title createGiottoLargeImage
#' @name createGiottoLargeImage
#' @description Creates a large giotto image that can be added to a Giotto subcellular object
#' @param raster_object terra SpatRaster image object
#' @param name name for the image
#' @param extent SpatExtent object to assign spatial extent. Takes priority.
#' @param image_transformations vector of sequential image transformations - under construction
#' @param xmax_adj adjustment of the maximum x-value to align the image
#' @param xmin_adj adjustment of the minimum x-value to align the image
#' @param ymax_adj adjustment of the maximum y-value to align the image
#' @param ymin_adj adjustment of the minimum y-value to align the image
#' @param scale_factor scaling of image dimensions relative to spatial coordinates
#' @return a giottoLargeImage object
#' @export
createGiottoLargeImage = function(raster_object,
                                  name = 'image',
                                  extent = NULL,
                                  image_transformations = NULL,
                                  xmax_bound = NULL,
                                  xmin_bound = NULL,
                                  ymax_bound = NULL,
                                  ymin_bound = NULL,
                                  scale_factor = 1) {
  
  
  # create minimum giotto
  g_imageL = giottoLargeImage(name = name,
                              raster_object = NULL,
                              overall_extent = NULL,
                              scale_factor = NULL,
                              resolution = NULL,
                              OS_platform = .Platform[['OS.type']])
  
  
  ## 1. check raster object and load as SpatRaster if necessary
  if(!methods::is(raster_object, 'SpatRaster')) {
    if(file.exists(raster_object)) {
      raster_object = try(terra::rast(raster_object))
      if(class(raster_object) == 'try-error') {
        stop(raster_object, ' can not be read by terra::rast() \n')
      }
    } else {
      stop("raster_object needs to be a'SpatRaster' object from the terra package or \n
           an existing path that can be read by terra::rast()")
    }
  }
  

  ## 2. image bound spatial extent
  
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
      ymax_bound = 0
      ymin_bound = -im_dim[2]
    } else {
      # Manual extent values
      if(is.null(xmax_bound) == TRUE) xmax_bound = 1
      if(is.null(xmin_bound) == TRUE) xmin_bound = 0
      if(is.null(ymax_bound) == TRUE) ymax_bound = 0
      if(is.null(ymin_bound) == TRUE) ymin_bound = -1
    }
    terra::ext(raster_object) = c(xmin_bound,xmax_bound,ymin_bound,ymax_bound)
  }
  
  
  ## 3. Assign raster_object to giottoLargeImage
  g_imageL@raster_object = raster_object
  
  ## 4. scale factor and resolution values
  g_imageL@resolution = terra::res(g_imageL@raster_object) # (x,y)
  names(g_imageL@resolution) = c('x','y')
  g_imageL@scale_factor = (1/g_imageL@resolution)
  
  ## 5. Set reasonable max intensity
  g_imageL@max_intensity = max(terra::spatSample(raster_object,
                                               size = 5000, # Defines the rough maximum of pixels allowed when resampling
                                               method = 'regular',
                                               value = TRUE))
  
  ## 6. extent object
  g_imageL@overall_extent = terra::ext(raster_object)
  
  ## 7. return image object
  return(g_imageL)
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


#' @title addGiottoImage
#' @name addGiottoImage
#' @description Adds giotto image objects to your giotto object
#' @param gobject giotto object
#' @param images list of giotto image objects, see \code{\link{createGiottoImage}}
#' @param spat_loc_name provide spatial location slot in Giotto to align images. Defaults to first one
#' @param scale_factor provide scale of image pixel dimensions relative to spatial coordinates.
#' @return an updated Giotto object with access to the list of images
#' @export
addGiottoImage = function(gobject,
                          images,
                          spat_loc_name = NULL,
                          scale_factor = NULL) {
  
  # 0. check params
  if(is.null(gobject)) stop('The giotto object that will be updated needs to be provided')
  
  if(is.null(images)) stop('The giotto image(s) that will be added needs to be provided')
  
  if(is.null(spat_loc_name)) {
    if(!is.null(gobject@spatial_locs)) {
      spat_loc_name = names(gobject@spatial_locs)[[1]]
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
                                           spat_loc_name = spat_loc_name)
          
          #Find spatial minmax values
          xmin_sloc = min(spatlocs$sdimx)
          xmax_sloc = max(spatlocs$sdimx)
          ymin_sloc = min(spatlocs$sdimy)
          ymax_sloc = max(spatlocs$sdimy)
          
          #Find adjustment values
          img_minmax = get_img_minmax(mg_img = im@mg_object)
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


#' @title showGiottoImageNames
#' @name showGiottoImageNames
#' @description Prints the available giotto images that are attached to the Giotto object
#' @param gobject a giotto object
#' @param verbose verbosity of function
#' @return a vector of giotto image names attached to the giotto object
#' @export
showGiottoImageNames = function(gobject,
                          verbose = TRUE) {

  if(is.null(gobject)) stop('A giotto object needs to be provided \n')
  g_image_names = names(gobject@images)

  if(verbose == TRUE) {
    cat('The following images are available: ',
        g_image_names, '\n')
  }

  return(g_image_names)

}


#' @title updateGiottoImage
#' @name updateGiottoImage
#' @description Updates the boundaries of a giotto image attached to a giotto object
#' @param gobject giotto object
#' @param image_name spatial locations
#' @param xmax_adj adjustment of the maximum x-value to align the image
#' @param xmin_adj adjustment of the minimum x-value to align the image
#' @param ymax_adj adjustment of the maximum y-value to align the image
#' @param ymin_adj adjustment of the minimum y-value to align the image
#' @param return_gobject return a giotto object
#' @return a giotto object or an updated giotto image if return_gobject = F
#' @export
updateGiottoImage = function(gobject,
                             image_name,
                             xmax_adj = 0,
                             xmin_adj = 0,
                             ymax_adj = 0,
                             ymin_adj = 0,
                             return_gobject = TRUE) {

  if(is.null(gobject)) stop('The giotto object that will be updated needs to be provided \n')
  if(is.null(image_name)) stop('The name of the giotto image that will be updated needs to be provided \n')

  g_image_names = names(gobject@images)
  if(!image_name %in% g_image_names) stop(image_name, ' was not found among the image names, see showImageNames()')

  # if image name is found, update the boundaries
  gobject@images[[image_name]]@boundaries = c(xmax_adj, xmin_adj, ymax_adj, ymin_adj)

  if(return_gobject == TRUE) {
    return(gobject)
  } else {
    return(gobject@images[[image_name]])
  }

}



#' @title rescaleGiottoImage
#' @name  rescaleGiottoImage
#' @description Change scale_factor and resolution slot values
#' @param gimage giottoImage to rescale
#' @param gobject giottoObject containing giottoImage to rescale
#' @param image_name name of giottoImage in giottoObject
#' @param spatloc_name name of spatlocs that scaling is relative to (defaults to first spatloc name for gobject or, if unavailable, "raw")
#' @param scale_factor scale factor to convert from spatial distance to pixel distance
#' @param resolution pixel distance per unit of real world distance (1/scale_factor)
#' @export
rescaleGiottoImage = function(gimage = NULL,
                              gobject = NULL,
                              image_name,
                              spatloc_name = NULL,
                              scale_factor = NULL,
                              resolution = NULL) {
  # Check Params
  if(is.null(scale_factor) && is.null(resolution)) stop('Either a scale factor or a resolution must be given \n')
  if(!is.null(scale_factor) && !is.null(resolution)) stop('Either a scale factor or a resolution must be given - not both \n')
  
  if(is.null(gimage) && is.null(gobject)) stop('Either a giottoImage or a giottoObject with image name is needed \n')
  
  if(!is.null(gobject) && is.null(image_name)) stop ('image_name must be provided if gobject is given \n')
  
  getFlag = NULL
  if(!is.null(gimage)) {
    getFlag = FALSE
  } else if(!is.null(gobject) && !is.null(image_name)) {
    getFlag = TRUE
  } else {
    stop('giottoImage must be given as either the gimage OR a gobject AND an image_name \n')
  }
  
  # Get giottoImage
  if(getFlag == TRUE) {
    g_img = getGiottoImage(gobject = gobject, image_name = image_name)
    if(is.null(spatloc_name)) {
      if(!is.null(gobject@spatial_locs)) {
        spatloc_name = names(gobject@spatial_locs)[[1]]
      } else {
        cat('No spatial locations have been found. Name default to "raw" \n')
        spatloc_name = 'raw'
      }
    }
  }
  if(getFlag == FALSE) {
    g_img = gimage
    if(is.null(spatloc_name)) {
      spatloc_name = 'raw'
    }
  } 
  
  # Assign values to giottoImage resolution and scale_factor slots
  if(!is.null(scale_factor)) {
    g_img@scale_factor[[spatloc_name]] = scale_factor
    g_img@resolution[[spatloc_name]] = 1/scale_factor
  }
  if(!is.null(resolution)) {
    g_img@scale_factor[[spatloc_name]] = 1/resolution
    g_img@resolution[[spatloc_name]] = resolution
  }
  
  # Return object
  if(getFlag == TRUE) {
    out = gobject
    out@images[[image_name]] = g_img
  }
  if(getFlag == FALSE) {
    out = g_img
  }
  
  return(out)

}



#' @title getGiottoImage
#' @name getGiottoImage
#' @description get a giotto image from a giotto object
#' @param gobject giotto object
#' @param image_name name of giotto image \code{\link{showGiottoImageNames}}
#' @return a giotto image
#' @export
getGiottoImage = function(gobject,
                          image_name) {

  if(is.null(gobject)) stop('The giotto object that will be updated needs to be provided \n')
  if(is.null(image_name)) stop('The name of the giotto image that will be updated needs to be provided \n')

  g_image_names = names(gobject@images)
  if(!image_name %in% g_image_names) stop(image_name, ' was not found among the image names, see showImageNames()')

  return(gobject@images[[image_name]])
}


#' @title plotGiottoImage
#' @name plotGiottoImage
#' @description get plot a giotto image from a giotto object
#' @param gobject giotto object
#' @param image_name name of giotto image \code{\link{showGiottoImageNames}}
#' @return plot
#' @export
plotGiottoImage = function(gobject,
                           image_name) {

  if(is.null(gobject)) stop('The giotto object that will be updated needs to be provided \n')
  if(is.null(image_name)) stop('The name of the giotto image that will be updated needs to be provided \n')

  g_image_names = names(gobject@images)
  if(!image_name %in% g_image_names) stop(image_name, ' was not found among the image names, see showImageNames()')

  graphics::plot(gobject@images[[image_name]]@mg_object)
}






# giottoLargeImage tools ####


#' @title getGiottoLargeImage
#' @name getGiottoLargeImage
#' @description get a giottoLargeImage from a giottoObject
#' @param gobject giotto object
#' @param largeImage_name name of giottoLargeImage
#' @return a giottoLargeImage
#' @export
getGiottoLargeImage = function(gobject = NULL,
                               largeImage_name = NULL) {

  if(is.null(gobject)) stop('The giotto object holding the giottoLargeImage needs to be provided \n')
  if(is.null(largeImage_name)) stop('The name of the giottoLargeImage needs to be provided \n')
  
  g_image_names = names(gobject@largeImages)
  if(!largeImage_name %in% g_image_names) stop(largeImage_name, ' was not found among the largeImage names. \n') #TODO See showImageNames()
    
  giottoLargeImage = gobject@largeImages$largeImage_name
  
  return(giottoLargeImage)
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
      giottoLargeImage = getGiottoLargeImage(gobject = gobject,
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


#' @title plotGiottoLargeImage
#' @name plotGiottoLargeImage
#' @description plot a downsampled version of giottoLargeImage
#' @param gobject giotto object
#' @param largeImage_name name of giottoLargeImage
#' @param zoom_extent extent object to focus on specific region of image
#' @param xmax_zoom assign zoom boundary
#' @param xmin_zoom assign zoom boundary
#' @param ymax_zoom assign zoom boundary
#' @param ymin_zoom assign zoom boundary
#' @param max_intensity value to treat as maximum intensity in color scale
#' @return plot
#' @export
plotGiottoLargeImage = function(gobject = NULL,
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
                   mar = c(5,5,1,1),
                   asp = 1)
  } else if(raster_object@ptr$rgb == TRUE) {
    
    # Determine likely image bitdepth
    if(is.null(max_intensity)) {
      bitDepth = ceiling(log(x = giottoLargeImage@max_intensity, base = 2))
      # Assign discovered bitdepth as max_intensity
      max_intensity = 2^bitDepth
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
#' @param mg_name name to assign converted magick image based giottoImage
#' @param spat_loc_name gobject spatial location name to map giottoImage to
#' @param crop_extent extent object to focus on specific region of image
#' @param xmax_crop assign crop boundary
#' @param xmin_crop assign crop boundary
#' @param ymax_crop assign crop boundary
#' @param ymin_crop assign crop boundary
#' @param resample_size maximum number of pixels to use when resampling
#' @param max_intensity value to treat as maximum intensity in color scale
convertGiottoLargeImageToMG = function(gobject = NULL,
                                       largeImage_name = NULL,
                                       giottoLargeImage = NULL,
                                       mg_name = NULL,
                                       spat_loc_name = NULL, #TODO
                                       crop_extent = NULL,
                                       xmax_crop = NULL,
                                       xmin_crop = NULL,
                                       ymax_crop = NULL,
                                       ymin_crop = NULL,
                                       resample_size = 500000,
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
  
  # Find boundary adj values
  overall_ext = giottoLargeImage@overall_extent
  current_ext = terra::ext(raster_object)
  
  xmin_adj = current_ext$xmin - overall_ext$xmin
  xmax_adj = overall_ext$xmax - current_ext$xmax
  ymin_adj = current_ext$ymin - overall_ext$ymin
  ymax_adj = overall_ext$ymax - current_ext$ymax
  
  names(xmin_adj) = NULL
  names(xmax_adj) = NULL
  names(ymin_adj) = NULL
  names(ymax_adj) = NULL
  
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
                              verbose = FALSE) #TODO
  
  # Set minimax
  #TODO make this compatible with spatlocs. The reference frame is weird right now
  g_image@minmax = c(current_ext$xmax,
                     current_ext$xmin,
                     current_ext$ymax,
                     current_ext$ymin)
  
  names(g_image@minmax) = c('xmax_sloc','xmin_sloc','ymax_sloc','ymin_sloc')
  
  # Set scalefactor
  scale_factor = ((unlist(magick::image_info(mImg)['width']))/dim(raster_object)[2])
  names(scale_factor) = NULL
  
  g_image = rescaleGiottoImage(gimage = g_image,
                               scale_factor = scale_factor)
  
  # return giottoImage
  return(g_image)
}
  




