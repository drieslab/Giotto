

# giottoImage class ####

#' @title S4 giottoImage Class
#' @description Framework of giotto object to store and work with spatial expression data
#' @keywords giotto, object
#' @slot name name of Giotto image
#' @slot mg_object magick image object
#' @slot minmax minimum and maximum of associated spatial location coordinates
#' @slot boundaries x and y coordinate adjustments (default to 0)
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
    OS_platform = "ANY"
  ),

  prototype = list(
    name = NULL,
    mg_object = NULL,
    minmax = NULL,
    boundaries = NULL,
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

    # print(object@mg_object)

  }
)


# giottoImage creation ####


#' @title createGiottoImage
#' @name createGiottoImage
#' @description Creates a giotto image that can be added to a Giotto object and/or used to add an image to the spatial plotting functions
#' @param gobject giotto object
#' @param spatial_locs spatial locations (alternative if giobject = NULL)
#' @param mg_object magick image object
#' @param name name for the image
#' @param image_transformations vector of sequential image transformations
#' @param do_manual_adj flag to use manual adj values instead of automatic alignment when given a gobject or spatlocs
#' @param xmax_adj adjustment of the maximum x-value to align the image
#' @param xmin_adj adjustment of the minimum x-value to align the image
#' @param ymax_adj adjustment of the maximum y-value to align the image
#' @param ymin_adj adjustment of the minimum y-value to align the image
#' @param scale_factor scaling of image dimensions relative to spatial coordinates
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
                             scale_factor = 1) {


  # create minimum giotto
  g_image = giottoImage(name = name,
                        mg_object = NULL,
                        minmax = NULL,
                        boundaries = NULL,
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
           an existig path that can be read by magick::image_read()")
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
  if(do_manual_adj == TRUE) cat('do_manual == TRUE \n','Boundaries will be adjusted by given values.\n')
  #If spatlocs or gobject supplied, minmax values will always be generated
  #If do_manual_adj == TRUE, bypass followup automatic boundary value generation
  if(!is.null(gobject)) {

    spatlocs = get_spatial_locations(gobject = gobject,
                                     spat_loc_name = spat_loc_name)

    #spatlocs = gobject@spatial_locs[['raw']]

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
    if(!all(colnames(spatlocs) %in% c('sdimx', 'sdimy'))) {
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
    warning('gobject or spatial locations are not provided \n',
            'Arbitrary values will be given \n')

    my_xmin = 0; my_xmax = 10; my_ymin = 0; my_ymax = 10

  }
  
  #minmax and boundary values for return
  g_image@minmax = c('xmax_sloc' = my_xmax,
                     'xmin_sloc' = my_xmin,
                     'ymax_sloc' = my_ymax,
                     'ymin_sloc' = my_ymin)

  g_image@boundaries = c('xmax_adj' = xmax_adj,
                         'xmin_adj' = xmin_adj,
                         'ymax_adj' = ymax_adj,
                         'ymin_adj' = ymin_adj)

  # image object
  return(g_image)
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
#' @param spat_loc_name provide spatial location slot in Giotto to align images. Defaults to 'raw'
#' @param scale_factor provide scale of image pixel dimensions relative to spatial coordinates. Defaults to 1:1
#' @return an updated Giotto object with access to the list of images
#' @export
addGiottoImage = function(gobject,
                          images,
                          spat_loc_name = 'raw',
                          scale_factor = 1) {

  if(is.null(gobject)) stop('The giotto object that will be updated needs to be provided')

  for(image_i in 1:length(images)) {

    im = images[[image_i]]

    if(methods::is(im, 'giottoImage')) {
      im_name = im@name

      all_im_names = names(gobject@images)

      if(im_name %in% all_im_names) {
        cat('\n ', im_name, ' has already been used, will be overwritten \n')
      }
      
      #Update boundaries if not already done during createGiottoImage() due to lack of spatlocs and gobject
      if(sum(im@boundaries == c(0,0,0,0)) == 4 && sum(im@minmax == c(10,0,10,0)) == 4) {
        spatlocs = get_spatial_locations(gobject = gobject,
                                         spat_loc_name = spat_loc_name)
        
        #Find spatial minmax values
        xmin_sloc = min(spatlocs$sdimx)
        xmax_sloc = max(spatlocs$sdimx)
        ymin_sloc = min(spatlocs$sdimy)
        ymax_sloc = max(spatlocs$sdimy)
        
        #Find adjustment values
        img_minmax = get_img_minmax(mg_img = im@mg_object)
        adj_values = get_adj_rescale_img(img_minmax = img_minmax,
                                         spatial_locs = spatlocs,
                                         scale_factor = scale_factor)
        
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
        
      }
      
      #Add giottoImage to gobject
      gobject@images[[im_name]] = im
      
    } else {
      warning('image is not a giotto image object')
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



#' @title getGiottoImage
#' @name getGiottoImage
#' @description get get a giotto image from a giotto object
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

