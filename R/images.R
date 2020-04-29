

#' @title convert_mgImage_to_array_DT
#' @name convert_mgImage_to_array_DT
#' @description converts a magick image object to a data.table
#' @param mg_object magick image object
#' @return data.table with image pixel information
convert_mgImage_to_array_DT = function(mg_object) {
  
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
#' @param mg_object magick image object
#' @param top_color_range top possible background colors to return
#' @return vector of pixel color frequencies and an associated barplot
#' @export
#' @examples
#'     estimateImageBg(mg_object)
estimateImageBg = function(mg_object, top_color_range = 1:50) {
  
  arrayDT = convert_mgImage_to_array_DT(mg_object = mg_object)
  sort_table = sort(table(arrayDT$RGB), decreasing = T)
  barplot(sort_table[top_color_range], col=names(sort_table[top_color_range]))
  
  cat('Most abundant pixel colors: \n')
  print(sort_table[top_color_range])
}


#' @title changeImageBg
#' @name changeImageBg
#' @description Function to change the background color of a magick image plot to another color
#' @param mg_object magick image object
#' @param bg_color estimated current background color
#' @param perc_range range around estimated background color to include (percentage)
#' @param new_color new background color 
#' @return vector of pixel color frequencies and an associated barplot
#' @export
#' @examples
#'     changeImageBg(mg_object)
changeImageBg = function(mg_object, bg_color, perc_range = 10, new_color = '#FFFFFF') {
  
  if(!is(mg_img, 'magick-image')) {
    stop("mg_object needs to be an image object 'magick-image'' from the magick package")
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
  c1_ind = arrayDT[['c.1']] > c1_min & arrayDT[['c.1']] < c1_max
  c2_ind = arrayDT[['c.2']] > c2_min & arrayDT[['c.2']] < c2_max
  c3_ind = arrayDT[['c.3']] > c3_min & arrayDT[['c.3']] < c3_max
  c_ind = c1_ind*c2_ind*c3_ind
  
  # replace old background with new background
  arrayDT[, 'c1' := ifelse(c_ind == T, c1_new, c.1)]
  arrayDT[, 'c2' := ifelse(c_ind == T, c2_new, c.2)]
  arrayDT[, 'c3' := ifelse(c_ind == T, c3_new, c.3)]
  
  # setorder for x and y coordinates
  data.table::setorder(arrayDT, y, x)
  
  # convert array_dt to array and then to magick image object
  original_width = magick::image_info(mg_img)[2]
  original_heigth = magick::image_info(mg_img)[3]
  myarray = array(as.vector(as.matrix(arrayDT[,.(c1, c2, c3)])), dim = c(original_width, original_heigth, 3))
  new_mg_object = magick::image_read(myarray)
  return(new_mg_object)
}


#' @title createGiottoImage
#' @name createGiottoImage
#' @description Creates a giotto image that can be added to a Giotto object and/or used to add an image to the spatial plotting functions
#' @param gobject giotto object
#' @param spatial_locs spatial locations
#' @param mg_object magick image object
#' @param name name for the image
#' @param xmax_adj adjustment of the maximum x-value to align the image
#' @param xmin_adj adjustment of the minimum x-value to align the image
#' @param ymax_adj adjustment of the maximum y-value to align the image
#' @param ymin_adj adjustment of the minimum y-value to align the image
#' @return a giotto image object
#' @export
#' @examples
#'     createGiottoImage(mg_object)
createGiottoImage = function(gobject = NULL,
                             spatial_locs = NULL,
                             mg_object,
                             name = 'image',
                             xmax_adj = 0, xmin_adj = 0, ymax_adj = 0, ymin_adj = 0) {
  if(!is(mg_img, 'magick-image')) {
    stop("mg_object needs to be an image object 'magick-image'' from the magick package")
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
                  minmax = c(my_xmax, my_xmin, my_ymax, my_ymin),
                  boundaries = c(xmax_adj, xmin_adj, ymax_adj, ymin_adj))
  
  class(imageObj) <- append(class(imageObj), 'imageGiottoObj')
  return(imageObj)
}


#' @title addImage
#' @name addImage
#' @description Adds giotto image objects to your giotto object
#' @param gobject giotto object
#' @param images list of giotto image objects, see \code{\link{createGiottoImage}}
#' @return an updated Giotto object with access to the list of images 
#' @export
#' @examples
#'     addImage(mg_object)
addImage = function(gobject,
                    images) {
  
  if(is.null(gobject)) stop('The giotto object that will be updated needs to be provided')
  
  for(image_i in 1:length(images)) {
    
    im = images[[image_i]]
    
    if(is(im, 'imageGiottoObj')) {
      im_name = im$name
      
      all_im_names = names(gobject@images)
      if(im_name %in% all_im_names) {
        cat('\n ', im_name, ' has already been used, will be overwritten \n')
      }
      
      gobject@images[[im_name]] = im
    } else {
      warning('image: ', im, ' is not a giotto image object')
    }
  }
  
  return(gobject)
  
}



#' @title addImageToSpatPlot
#' @name addImageToSpatPlot
#' @description Add a giotto image to a spatial ggplot object post creation
#' @param spatpl a spatial ggplot object
#' @param gimage a giotto image, see \code{\link{createGiottoImage}}
#' @return an updated spatial ggplot object
#' @export
#' @examples
#'     addImageToSpatPlot(mg_object)
addImageToSpatPlot = function(spatpl = NULL,
                              gimage = NULL) {
  
  
  if(is.null(spatpl) | is.null(gimage)) {
    stop('A spatial ggplot object and a giotto image need to be given')
  }
  
  # extract min and max from object
  my_xmax = gimage$minmax[1]
  my_xmin = gimage$minmax[2]
  my_ymax = gimage$minmax[3]
  my_ymin = gimage$minmax[4]
  
  # convert giotto image object into array
  img_array = as.numeric(gimage$mg_object[[1]])
  
  # extract adjustments from object
  xmax_b = gimage$boundaries[1]
  xmin_b = gimage$boundaries[2]
  ymax_b = gimage$boundaries[3]
  ymin_b = gimage$boundaries[4]
  
  newpl = spatpl + annotation_raster(img_array,
                                     xmin = my_xmin-xmin_b, xmax = my_xmax+xmax_b,
                                     ymin = my_ymin-ymin_b, ymax = my_ymax+ymax_b)
  
  # move image to background
  nr_layers = length(newpl$layers)
  newpl$layers = c(newpl$layers[[nr_layers]], newpl$layers[1:(nr_layers-1)])
  
  return(newpl)
  
}




