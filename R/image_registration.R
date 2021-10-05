### Image registration and creation of registered Giotto object ####


#' @name trakem2_rigid_transf_extract
#' @title trakem2_rigid_transf_extract
#' @description Extract rigid registration transformation values from FIJI TrakEM2 xml file. Generated through register_virtual_stack_slices.
#' @param inputstring string read in from TrakeEM2 xml file
#' @keywords internal
trakem2_rigid_transf_extract = function(inputstring) {

  #Catch wrong inputs
  if (grepl('^.*trakem2.*', inputstring, ignore.case = T) != 1) {
    stop('xml files must be in TrakeEM2 format')
  }


  #Regex to find the values from the TrakEM2 .xml
  transfExtractPatA = '.*<iict_transform class=\\"mpicbg.trakem2.transform.RigidModel2D\\" data=\\"(.*?)\\" />.*'
  transfExtractPatB = '.*class=\\"mpicbg.trakem2.transform.TranslationModel2D\\" data=\\"(.*?)\\" />.*'
  #Split relevant text into numerical values
  out <- c(sapply(strsplit(regmatches(x = inputstring,
                                      m = regexec(pattern = transfExtractPatA,
                                                  text = inputstring))[[1]][2],
                           split = ' '),
                  function(x) as.numeric(x)),
           sapply(strsplit(regmatches(x = inputstring,
                                      m = regexec(pattern = transfExtractPatB,
                                                  text = inputstring))[[1]][2],
                           split = ' '),
                  function(x) as.numeric(x)))

  if(sum(is.na(out)) == 2) {
    out = rep(0,5)
  }

  out = c(out,0,0)
  out = data.table::data.table(t(matrix(out)))
  colnames(out) = c('Theta','Xtransform','Ytransform','itx','ity','XFinalTransform','YFinalTransform')

  #itx and ity are additional values in the trakem2 xml files that must be added to Xtransform and Ytransform in order to get the final transformation values.
  #only relevant for sampleset with more than 1 slice away from the reference image
  out$XFinalTransform = out$Xtransform + out$itx
  out$YFinalTransform = out$Ytransform + out$ity

  #Multiply theta by -1 due to differences in R and image plotting coordinates
  out$Theta <- -out$Theta

  return(out)
}


#' @title scale_spatial_locations
#' @name scale_spatial_locations
#' @description Scale the X and Y coordinates given by the scale factor input
#' @param spatlocs spatial locations to scale
#' @param scale_factor scaling factor to apply to coordinates
#' @keywords internal
scale_spatial_locations = function(spatlocs,
                                    scale_factor) {

  spatlocs$sdimx = spatlocs$sdimx*scale_factor
  spatlocs$sdimy = spatlocs$sdimy*scale_factor

  return(spatlocs)
}


#' @title rotate_spatial_locations
#' @name rotate_spatial_locations
#' @description Rotate given X Y coordinates by given radians in counter clockwise manner about the coordinate origin
#' @param spatlocs spatial locations to use
#' @param rotateradians radians by which the x and y values will be rotated in a counter clockwise manner
#' @keywords internal
rotate_spatial_locations = function(spatlocs,
                                     rotateradians) {

  xvals = spatlocs$sdimx
  yvals = spatlocs$sdimy

  spatlocs$sdimx = xvals*cos(rotateradians) + yvals*sin(rotateradians)
  spatlocs$sdimy = -xvals*sin(rotateradians) + yvals*cos(rotateradians)

  return(spatlocs)
}

#' @title xy_translate_spatial_locations
#' @name xy_translate_spatial_locations
#' @description Translate given X Y coordinates by given x and y translation values
#' @param spatlocs spatial locations to use
#' @param xtranslate value to translate coordinates in the positive x direction
#' @param ytranslate value to translate coordinates in the positive y direction
#' @keywords internal
xy_translate_spatial_locations = function(spatlocs,
                                           xtranslate,
                                           ytranslate) {

  spatlocs$sdimx = spatlocs$sdimx + xtranslate
  spatlocs$sdimy = spatlocs$sdimy + ytranslate

  return(spatlocs)
}


#' @title rigid_transform_spatial_locations
#' @name rigid_transform_spatial_locations
#' @description Performs appropriate transforms to align spatial locations with registered images.
#' @param spatlocs input spatial locations
#' @param transform_values transformation values to use
#' @param method which method is used for image registration
#' @keywords internal
#Rotation is performed first, followed by XY transform.
rigid_transform_spatial_locations = function(spatlocs,
                                             transform_values,
                                             method) {
  if(method == 'fiji') {
    spatlocsXY = spatlocs[,c('sdimx','sdimy')]
    #These functions must be performed in positive y values
    spatlocsXY$sdimy = -1 * spatlocsXY$sdimy
  
    spatlocsXY = rotate_spatial_locations(spatlocsXY,
                                          transform_values$Theta)
  
    spatlocsXY = xy_translate_spatial_locations(spatlocsXY,
                                                 transform_values$XFinalTransform,
                                                 transform_values$YFinalTransform)
    
    spatlocs$sdimx = spatlocsXY$sdimx
    spatlocs$sdimy = -1 * spatlocsXY$sdimy
  
    return(spatlocs)
    
  } else if(method == 'rvision') {
    
    spatLocsXY = spatlocs[,c('sdimx','sdimy')]
    spatLocsXY = rotate_spatial_locations(spatLocsXY,acos(transform_values[1,1]))
    spatLocsXY = xy_translate_spatial_locations(spatLocsXY,-transform_values[1,3],-transform_values[2,3])
    
    spatlocs$sdimx = spatLocsXY[,1]
    spatlocs$sdimy = spatLocsXY[,2]
    
    return(spatlocs)
    
  } else {
    stop('Image registration method must be provided. Only "fiji" and "rvision" methods currently supported.')
  }
}

#' @title reg_img_minmax_finder
#' @name reg_img_minmax_finder
#' @description finds new minmax boundaries of registration transformed images
#' @param gobject_list list of gobjects to use
#' @param image_unreg name of original unregistered images
#' @param scale_factor scale factors for registered images relative to spatlocs.
#' @param transform_values transformation values to use
#' @param spatloc_unreg unregistered spatloc name to use
#' @param method method of registration
#' @keywords internal
#Automatically account for changes in image size due to alignment
reg_img_minmax_finder = function(gobject_list,
                                 image_unreg,
                                 scale_factor,
                                 transform_values,
                                 spatloc_unreg,
                                 method) {

  #Find image spatial info from original image if possible
  #Check to make sure that image_unreg finds an existing image in each gobject to be registered
  imgPresent = function(gobject,image) {
    imgPresent = FALSE
    if(image %in% showGiottoImageNames(gobject = gobject, verbose = FALSE)) {
      imgPresent = TRUE
    }
    return(imgPresent)
  }

  if(all(as.logical(lapply(gobject_list, imgPresent, image = image_unreg)))) {

    giottoImage_list = lapply(gobject_list, getGiottoImage, image_name = image_unreg)
    image_corners = lapply(giottoImage_list, get_img_corners)

    # Infer image corners of registered images PRIOR TO REGISTRATION
    # scale unreg_image corners to registered image (use reg_scalefactor/unreg_scalefactor as scale factor)
    image_corners = lapply(1:length(gobject_list),
                           function(x) {
                             scale_spatial_locations(spatlocs = image_corners[[x]],
                                                     scale_factor = (scale_factor[[x]]/giottoImage_list[[x]]@scale_factor[[spatloc_unreg]]))
                           })

    # register corners based on transform values (only possible at reg_image scaling)
    image_corners_reg = lapply(1:length(image_corners),
                               function(x) {
                                 rigid_transform_spatial_locations(spatlocs = image_corners[[x]],
                                                                   transform_values = transform_values[[x]],
                                                                   method = method)
                               })
    
    # Return registered corners to spatloc scaling
    image_corners_reg = lapply(1:length(image_corners_reg),
                               function(x) {
                                 scale_spatial_locations(spatlocs = image_corners_reg[[x]],
                                                         scale_factor = (1/scale_factor[[x]]))
                               })

    # combine list then find new image bound minmax
    image_corners_reg = do.call(rbind, image_corners_reg)
    minmaxRegVals = list('xmax_reg' = max(image_corners_reg$sdimx),
                         'xmin_reg' = min(image_corners_reg$sdimx),
                         'ymax_reg' = max(image_corners_reg$sdimy),
                         'ymin_reg' = min(image_corners_reg$sdimy))

    #return the minmax values - already scaled to spatlocs
    return(minmaxRegVals)
  } else {
    warning('Original images must be supplied for registered images to be aligned.\n')
  }
}


#' @title get_img_corners
#' @name get_img_corners
#' @description finds four corner spatial coords of giottoImages or magick-images
#' @param img_object giottoImage or magick-image to use
#' @keywords internal
get_img_corners = function(img_object) {
  if(methods::is(img_object,'giottoImage')) {
    img_dims = Giotto:::get_img_minmax(img_object@mg_object)
  } else if(methods::is(img_object,'magick-image')) {
    img_dims = Giotto:::get_img_minmax(img_object)
  } else {
    stop('img_object must be either a giottoImage or a magick-image \n')
  }

  upper_left = c(img_dims$img_xmin,img_dims$img_ymax)
  lower_left = c(img_dims$img_xmin,img_dims$img_ymin)
  lower_right = c(img_dims$img_xmax,img_dims$img_ymin)
  upper_right = c(img_dims$img_xmax,img_dims$img_ymax)

  imageCorners = rbind(upper_left,
                       lower_left,
                       lower_right,
                       upper_right)
  colnames(imageCorners) = c('sdimx','sdimy')
  imageCorners = as.data.frame(imageCorners)
  return(imageCorners)
}


### exported functions ####


#' @title registerGiottoObjectList
#' @name registerGiottoObjectList
#' @description Wrapper function for registerGiottoObjectListFiji and registerGiottoObjectListRvision
#' @param gobject_list List of gobjects to register
#' @param method Method used to align gobjects. Current options are either using FIJI register_virtual_stack_slices output or rvision
#' @param image_unreg Gobject image slot to use. Defaults to 'image' (optional)
#' @param image_reg_name Arbitrary image slot name for registered images to occupy. Defaults to replacement of 'image' slot (optional)
#' @param image_list RVISION - under construction
#' @param save_dir RVISION - under construction
#' @param spatloc_unreg Unregistered spatial locations to align. Defaults to 'raw' slot (optional)
#' @param spatloc_reg_name Arbitrary name for registered spatial locations. Defaults to replacement of 'raw' slot (optional)
#' @param fiji_xml_files Filepaths to FIJI registration XML outputs
#' @param fiji_registered_images Registered images output by FIJI register_virtual_stack_slices
#' @param scale_factor Scaling to be applied to spatial coordinates
#' @param allow_rvision_autoscale Whether or not to allow rvision to automatically scale the images when performing image registration
#' @param verbose Be verbose
#' @return List of registered giotto objects where the registered images and spatial locations
#' @export
registerGiottoObjectList = function(gobject_list,
                                    method = c('fiji','rvision'),
                                    image_unreg = 'image',
                                    image_reg_name = 'image',
                                    image_list = NULL, #Rvision
                                    save_dir = NULL, #Rvision
                                    spatloc_unreg = 'raw',
                                    spatloc_reg_name = 'raw',
                                    fiji_xml_files,
                                    fiji_registered_images,
                                    scale_factor = NULL,
                                    allow_rvision_autoscale = TRUE, #Rvision
                                    # auto_comp_reg_border = TRUE,
                                    verbose = TRUE) {

  if(method == 'fiji') {
    gobject_list = registerGiottoObjectListFiji(gobject_list = gobject_list,
                                                image_unreg = image_unreg,
                                                image_reg_name = image_reg_name,
                                                registered_images = fiji_registered_images,
                                                spatloc_unreg = spatloc_unreg,
                                                spatloc_reg_name = spatloc_reg_name,
                                                xml_files = fiji_xml_files,
                                                scale_factor = scale_factor,
                                                # auto_comp_reg_border = auto_comp_reg_border,
                                                verbose = verbose)

  } else if (method == 'rvision') {
    gobject_list = registerGiottoObjectListRvision(gobject_list = gobject_list,
                                                   image_list = image_list,
                                                   save_dir = save_dir,
                                                   spatloc_unreg = spatloc_unreg,
                                                   spatloc_reg_name = spatloc_reg_name,
                                                   verbose = verbose)

  } else {
    stop('Invalid method input\n Only fiji and rvision methods are currently supported.')
  }

  return(gobject_list)
}


#' @title registerGiottoObjectListFiji
#' @name registerGiottoObjectListFiji
#' @description Function to spatially align gobject data based on FIJI image registration.
#' @param gobject_list list of gobjects to register
#' @param image_unreg name of original unregistered images. Defaults to 'image' (optional)
#' @param image_reg_name Arbitrary name for registered images to occupy. Defaults to replacement of 'image' (optional)
#' @param registered_images registered images output by FIJI register_virtual_stack_slices
#' @param spatloc_unreg spatial locations to use. Defaults to 'raw' (optional)
#' @param spatloc_reg_name name for registered spatial locations. Defaults to replacement of 'raw' (optional)
#' @param xml_files atomic vector of filepaths to xml outputs from FIJI register_virtual_stack_slices
#' @param scale_factor vector of scaling factors of images used in registration vs spatlocs
#' @param verbose be verbose
#' @return list of registered giotto objects where the registered images and spatial locations
#' @export
registerGiottoObjectListFiji = function(gobject_list,
                                        image_unreg = 'image',
                                        image_reg_name = 'image',
                                        registered_images = NULL,
                                        spatloc_unreg = 'raw',
                                        spatloc_reg_name = 'raw',
                                        xml_files,
                                        scale_factor = NULL,
                                        verbose = TRUE) {

  ## 0. Check Params ##
  if(length(gobject_list) != length(xml_files)) {
    stop('xml spatial transforms must be supplied for every gobject to be registered.\n')
  }
  
  if(is.null(registered_images) == FALSE) {
    # If there are not the same number of registered images as gobjects, stop
    if(length(registered_images) != length(gobject_list)) {
      stop('A registered image should be supplied for every gobject to align \n')
    }
    if(sum(as.logical(lapply(registered_images, methods::is, class2 = 'giottoImage'))) > 0) {
      stop('Registered images should be supplied as either magick-objects or filepaths \n')
    }
  }
  
  if(!is.null(scale_factor)) {
    if(!is.numeric(scale_factor)) {
      stop('scale_factor only accepts numerics')
    }
    if((length(scale_factor) != length(gobject_list)) && (length(scale_factor) != 1)) {
      stop('If more than one scale_factor is given, there must be one for each gobject to be registered. \n')
    }
  }
  
  
  # scale_factors will always be given externally. Registered images do not have gobjects yet.
  # expand scale_factor if given as a single value
  scale_list = c()
  if(length(scale_factor) == 1) {
    scale_list = rep(scale_factor, length(gobject_list))
  } else {
    scale_list = unlist(scale_factor) # ensure atomic vector
  }
  
  
  ## 1. Get spatial coordinates and put in lists #
  spatloc_list = list()
  for(gobj_i in 1:length(gobject_list)) {
    gobj = gobject_list[[gobj_i]]
    spatloc = get_spatial_locations(gobject = gobj,
                                    spat_loc_name = spatloc_unreg)
    #------ Put all spatial location data together
    spatloc_list[[gobj_i]] = spatloc
  }

  
  
  ## 2. read transform xml files into list ##
  transf_list = list()
  for(file_i in 1:length(xml_files)) {
    t_file = xml_files[[file_i]]
    #------ Put all transform files together
    transf_list[[file_i]] = paste(readLines(t_file, warn = F), collapse = '\n')
  }

  # Select useful info out of the TrakEM2 files
  transformsDF = lapply(transf_list,
                        FUN = trakem2_rigid_transf_extract)
  
  
  
  ## 3. apply transformation on spatial locations ##
  # Scale by registered image's scale_factor
  spatloc_list = lapply(1:length(spatloc_list),
                        function(x) {
                          scale_spatial_locations(spatlocs = spatloc_list[[x]],
                                                  scale_factor = scale_list[x])
                        })

  # Register spatial locations
  spatloc_list = lapply(1:length(spatloc_list),
                        function(x) {
                          rigid_transform_spatial_locations(spatlocs = spatloc_list[[x]],
                                                            transform_values = transformsDF[[x]],
                                                            method = 'fiji')
                        })
  
  # Return scaling to spatloc original
  spatloc_list = lapply(1:length(spatloc_list),
                        function(x) {
                          scale_spatial_locations(spatlocs = spatloc_list[[x]],
                                                  scale_factor = 1/(scale_list[x]))
                        })
  
  ## 4. update Giotto slots and names and return list of Giotto object

  #Find new image boundaries for registered images
  #Must have original pre-registration images in the gobject for this to work
  #TODO (optional if just registering spatlocs)
  reg_img_boundaries = reg_img_minmax_finder(gobject_list = gobject_list,
                                             image_unreg = image_unreg,
                                             scale_factor = scale_list,
                                             transform_values = transformsDF,
                                             spatloc_unreg = spatloc_unreg,
                                             method = 'fiji')

  # Gobject updating for loop
  for(gobj_i in 1:length(gobject_list)) {
    gobj = gobject_list[[gobj_i]]

    
    # Params check for conflicting names
    if(verbose == TRUE) {
      if(image_unreg == image_reg_name) {
        cat('Image original name matches output name. Unregistered images renamed "unregistered". \n')
      }
      if(spatloc_unreg == spatloc_reg_name) {
        cat('Spatloc original name matches output name. Unregistered spatlocs renamed "unregistered". \n')
      }
    }

    
    # Update Spatial
    #Rename original spatial locations to 'unregistered' if conflicting with output
    if(spatloc_unreg == spatloc_reg_name) {
      gobj = set_spatial_locations(gobject = gobj,
                                   spat_loc_name = 'unregistered',
                                   spatlocs = get_spatial_locations(gobject = gobj,
                                                                    spat_loc_name = spatloc_unreg))
    }


    #Assign registered spatial locations from spatloc_list to gobject_list
    gobj = set_spatial_locations(gobject = gobj,
                                 spat_loc_name = spatloc_reg_name,
                                 spatlocs = spatloc_list[[gobj_i]])
    


    # Update images
    #If there is an existing image with the image_reg_name, rename it "unregistered"
    #Move the original image to 'unregistered'
    if(image_unreg == image_reg_name) {
      gobj@images$unregistered = gobj@images[[image_unreg]]
    }

    
    #Create a giotto image if there are registered images supplied
    if(!is.null(registered_images)) {
      g_image = createGiottoImage(gobject = gobj,
                                  spatial_locs = spatloc_reg_name,
                                  mg_object = registered_images[[gobj_i]],
                                  name = image_reg_name,
                                  scale_factor = scale_list[[gobj_i]])

      #Add the registered image to the gobj.
      gobj = addGiottoImage(gobject = gobj, images = list(g_image))

      #Automatic adjustment
      if(exists('reg_img_boundaries')){ #TODO
        im_info = gobj@images[[image_reg_name]]@minmax

        #update the giottoImage boundaries
        boundaries <- as.numeric(c(reg_img_boundaries$xmax_reg - im_info[['xmax_sloc']],
                                   im_info[['xmin_sloc']] - reg_img_boundaries$xmin_reg,
                                   reg_img_boundaries$ymax_reg - im_info[['ymax_sloc']],
                                   im_info[['ymin_sloc']] - reg_img_boundaries$ymin_reg))

        names(boundaries) = c('xmax_adj','xmin_adj','ymax_adj','ymin_adj')

        gobj@images[[image_reg_name]]@boundaries = boundaries

      }
    }
    gobject_list[[gobj_i]] = gobj
  } # gobj update loop end
  return(gobject_list)
}

#TODO check if spatloc is actually provided in createGiottoImage() and ignore auto align if not.

#' @title registerGiottoObjectListRvision
#' @name registerGiottoObjectListRvision
#' @description Function to spatially align gobject data based on Rvision image registration.
#' @param gobject_list list of gobjects to register
#' @param image_list Filepaths to unregistered images
#' @param save_dir (Optional) If given, save registered images to this directory
#' @param spatloc_unreg spatial locations to use
#' @param spatloc_reg_name name for registered spatial locations to. Defaults to replacement of spat_unreg (optional)
#' @param verbose be verbose
#' @return list of registered giotto objects where the registered images and spatial locations
#' @export
#Register giotto objects when given raw images and spatial locations
registerGiottoObjectListRvision = function(gobject_list = gobject_list,
                                           image_list = NULL,
                                           save_dir = NULL,
                                           spatloc_unreg = NULL,
                                           spatloc_reg_name = 'raw',
                                           verbose = TRUE) { #Not used
  
  ## 1. get spatial coordinates and put in list ##
  spatloc_list = list()
  for(gobj_i in 1:length(gobject_list)) {
    gobj = gobject_list[[gobj_i]]
    spatloc = get_spatial_locations(gobject = gobj, ###Tag for editing later
                                    spat_loc_name = spat_loc_name)
    # Put all spatial location data together
    spatloc_list[[gobj_i]] = spatloc
  }
  
  ## 2. Load images into list
  if (length(spatloc_list) != length(image_list)) {
    stop('images must be supplied for every gobject to be registered.')
  }
  
  unreg_images <- c()
  color_images <- c()
  for (path in image_list) {
    unreg_images <- append(unreg_images, Rvision::image(filename = path), after = length(unreg_images))
    color_images <- append(color_images, Rvision::image(filename = path), after = length(color_images))
  }
  
  ## 3. Perform preprocessing
  rows <- c()
  cols <- c()
  for (image_i in 1:length(unreg_images)) {
    # Make images grayscale
    Rvision::changeColorSpace(unreg_images[[image_i]], colorspace = "GRAY", target = "self")
    # Retrieve image dimensions
    dims <- Rvision:::dim.Rcpp_Image(unreg_images[[image_i]])
    rows <- append(rows, dims[[1]], after = length(rows))
    cols <- append(cols, dims[[2]], after = length(cols))
  }
  maxes <- c(max(cols),max(rows))
  squmax <- max(maxes)
  rm(dims, maxes)
  
  enddim <- 500
  for (i in 1:length(unreg_images)) {
    # Add border so all images have same square dimensions
    Rvision::border(unreg_images[[i]], squmax-rows[[i]], 0, squmax-cols[[i]], 0, border_color = "white", target = "self")
    Rvision::border(color_images[[i]], squmax-rows[[i]], 0, squmax-cols[[i]], 0, border_color = "white", target = "self")
    # Apply scaling so all images of reasonable size for processing
    unreg_images[[i]] <- Rvision::resize(unreg_images[[i]], height = enddim, width = enddim, target = "new")
    color_images[[i]] <- Rvision::resize(color_images[[i]], height = enddim, width = enddim, target = "new")
  }
  rm(cols,rows)
  
  ## 4. Compute transformations
  # Choose reference image
  refImage <- unreg_images[[base::floor(length(unreg_images)/2)]]
  
  # Compute ECC transforms
  transfs <- base::vector(mode = "list", length = length(unreg_images))
  for (i in 1:length(unreg_images)) {
    transfs[[i]] <- Rvision::findTransformECC(refImage, unreg_images[[i]], warp_mode = "euclidean", filt_size = 101)
  }
  rm(refImage)
  
  ## 5. Apply transform
  reg_images <- c()
  for (i in 1:length(unreg_images)) {
    # Apply scaling
    spatloc_list[[i]] <- scale_spatial_locations(spatloc_list[[i]], enddim/squmax)
    # Apply transform to spatlocs
    spatloc_list[[i]] <- rigid_transform_spatial_locations(spatloc_list[[i]], transfs[[i]], method = 'rvision')
  }
  rm(squmax, enddim)
  
  ## 6. Update giotto object 
  for(gobj_i in 1:length(gobject_list)) {
    gobj = gobject_list[[gobj_i]]
    #Rename original spatial locations to 'unregistered'
    gobj@spatial_locs$unregistered <- gobj@spatial_locs$spat_loc_values
    
    #Assign registered spatial locations from spatloc_list to gobject_list
    gobj@spatial_locs$spat_loc_values <- spatloc_list[[gobj_i]]
    
    gobject_list[[gobj_i]] <- gobj
  }
  
  ## 7. Save transformed images
  if (!is.null(save_dir)) {
    # Apply transform to image
    transf_images <- c()
    for (i in 1:length(unreg_images)) {
      transf_images <- append(transf_images, Rvision::warpAffine(color_images[[i]], transfs[[i]], target = "new"), length(transf_images))
    }
    # Save images to save directory
    for(image_i in 1:length(transf_images)) {
      name <- paste(save_dir, image_i, ".jpg")
      Rvision::write.Image(transf_images[[image_i]], name)
    } 
  }
  
  return(gobject_list)
}

### FIJI related functions ####

#' @title fiji
#' @description \code{fiji} returns path to preferred Fiji executable
#' @rdname runFijiMacro
#' @export
#' @examples
#' # Path to current Fiji executable
#' \donttest{
#' fiji()
#' }
#'
#' \dontrun{
#' # This function was taken and modified from jimpipeline by jefferislab #
#' # Set path to preferred Fiji executable (this will be remembered)
#' # you can also set options(giotto.fiji="/some/path")
#' fiji("/Applications/Fiji.app/Contents/MacOS/ImageJ-macosx")
#' }
fiji = function(fijiPath = NULL) {
  if(!is.null(fijiPath)) {
    if(!file.exists(fijiPath))
      stop("fiji is not at: ", fijiPath)
  } else {
    # do we have an option set?
    fijiPath=getOption('giotto.fiji')
    if(!is.null(fijiPath)) {
      if(!file.exists(fijiPath))
        stop("fiji is not at: ", fijiPath, " as specified by options('giotto.fiji')!")
    } else {
      # look for it in sensible places
      if(!nzchar(fijiPath <- Sys.which('fiji'))) {
        macapp="/Applications/Fiji.app/Contents/MacOS/ImageJ-macosx"
        if(file.exists(macapp))
          fijiPath = macapp
        else
          stop("Unable to find fiji!",
               "Set options('giotto.fiji') to point to the fiji command line executable!")
      }
    }
  }
  options(giotto.fiji=fijiPath)
  normalizePath(fijiPath)
}


#' @title registerImagesFIJI
#' @name registerImagesFIJI
#' @description Wrapper function for Register Virtual Stack Slices plugin in FIJI
#' @param source_img_dir Folder containing images to be registered
#' @param output_img_dir Folder to save registered images to
#' @param transforms_save_dir (jython implementation only) Folder to save transforms to
#' @param ref_img_name (jython implementation only) File name of reference image for the registration
#' @param init_gauss_blur Point detector option: initial image blurring
#' @param steps_per_scale_octave Point detector option
#' @param min_img_size Point detector option
#' @param max_img_size Point detector option
#' @param feat_desc_size Feature descriptor option
#' @param feat_desc_orient_bins Feature descriptor option
#' @param closest_next_closest_Ratio Feature descriptor option
#' @param max_align_err Geometric consensus filter option
#' @param inlier_ratio Geometric consensus filter option
#' @param headless Whether to have ImageJ/Fiji running headless #TODO
#' @param batch Use batch mode #TODO
#' @param MinMem,MaxMem Memory limits
#' @param IncrementalGC Whether to use incremental garbage collection
#' @param Threads Number of threads
#' @param fijiArgs Arguments for ImageJ/FIJI
#' @param javaArgs Arguments for Java
#' @param ijArgs Arguments for ImageJ
#' @param jython Use jython wrapper script
#' @param fijiPath Path to fiji executable (can be set by
#'   \code{options(giotto.fiji="/some/path")})
#' @param DryRun Whether to return the command to be run rather than actually
#'   executing it.
#' @return list of registered giotto objects where the registered images and spatial locations
#' \dontrun{
#' #This function was adapted from runFijiMacro function in jimpipeline by jefferislab #
#' }
#' @export
registerImagesFIJI = function(source_img_dir,
                              output_img_dir,
                              transforms_save_dir,
                              ref_img_name,
                              #Scale Invariant Interest Point Detector Options
                              init_gauss_blur = 1.6,
                              steps_per_scale_octave = 3,
                              min_img_size = 64,
                              max_img_size = 1024,
                              #Feature Descriptor Options
                              feat_desc_size = 8,
                              feat_desc_orient_bins = 8,
                              closest_next_closest_ratio = 0.92,
                              #Geometric Consensus Filter Options
                              max_align_err = 25,
                              inlier_ratio = 0.05,
                              #FIJI Options
                              headless = FALSE,
                              batch = TRUE,
                              MinMem = MaxMem,
                              MaxMem = 2500,
                              IncrementalGC = TRUE,
                              Threads = NULL,
                              fijiArgs = NULL,
                              javaArgs = NULL,
                              ijArgs = NULL,
                              jython = FALSE,
                              fijiPath = fiji(),
                              DryRun = FALSE) {

  #Check if output directory exists. If not, create the directory
  if(!file.exists(output_img_dir)) {
    dir.create(output_img_dir)
  }
  
  #expand the paths of source and output
  source_img_dir = path.expand(source_img_dir)
  output_img_dir = path.expand(output_img_dir)

  
  if(headless) fijiArgs = c(fijiArgs,"--headless")
  fijiArgs=paste(fijiArgs,collapse=" ")
  
  javaArgs=c(paste("-Xms",MinMem,'m',sep=""),paste("-Xmx",MaxMem,'m',sep=""),javaArgs)
  if(IncrementalGC) javaArgs=c(javaArgs,"-Xincgc")
  javaArgs=paste(javaArgs,collapse=" ")

  threadAdjust=ifelse(is.null(Threads),"",paste("run(\"Memory & Threads...\", \"parallel=",Threads,"\");",sep=""))
  
  if(jython == TRUE) {
    #TODO Add check to see if jython script is installed.
    cat('jython implementation requires Headless_RVSS.py in "/Giotto/inst/fiji/" to be copied to "/Applications/Fiji.app/plugins/Scripts/MyScripts/Headless_RVSS.py" \n')
    
    macroCall=paste(" -eval '",
                    threadAdjust,
                    "run(\"Headless RVSS\", \"source_dir=[",
                    source_img_dir,
                    "] target_dir=[",
                    output_img_dir,
                    "] transf_dir=[",
                    transforms_save_dir,
                    "] reference_name=[",
                    ref_img_name,
                    "] init_gauss_blur=",
                    init_gauss_blur,
                    " steps_per_scale_octave=",
                    steps_per_scale_octave,
                    " min_img_size=",
                    min_img_size,
                    " max_img_size=",
                    max_img_size,
                    " feat_desc_size=",
                    feat_desc_size,
                    " feat_desc_orient_bins=",
                    feat_desc_orient_bins,
                    " closest_next_closest_ratio=",
                    closest_next_closest_ratio,
                    " max_align_err=",
                    max_align_err,
                    " minInlierRatio=",
                    inlier_ratio,
                    " interpolate=TRUE\");' ",
                    sep="")
  } else {
    macroCall=paste(" -eval '",
                    threadAdjust,
                    "run(\"Register Virtual Stack Slices\", \"source=[",
                    source_img_dir,
                    "] output=[",
                    output_img_dir,
                    "] feature=Similarity registration=[Rigid                -- translate + rotate                  ] advanced save initial_gaussian_blur=",
                    init_gauss_blur,
                    " steps_per_scale_octave=",
                    steps_per_scale_octave,
                    " minimum_image_size=",
                    min_img_size,
                    " maximum_image_size=",
                    max_img_size,
                    " feature_descriptor_size=",
                    feat_desc_size,
                    " feature_descriptor_orientation_bins=",
                    feat_desc_orient_bins,
                    " closest/next_closest_ratio=",
                    closest_next_closest_ratio,
                    " maximal_alignment_error=",
                    max_align_err,
                    " inlier_ratio=",
                    inlier_ratio,
                    " feature_extraction_model=Similarity registration_model=[Rigid                -- translate + rotate                  ] interpolate\");' ",sep="")
  }
  
  ijArgs=paste(c(ijArgs,ifelse(batch,"-batch","")),collapse=" ")

  cmd<-paste(fijiPath,javaArgs,fijiArgs,"--",macroCall,ijArgs)
  if(DryRun) return(cmd)
  return(0==system(cmd))
}




#TODO - merge jython function into normal register FIJI
#TODO - add in manual rigid registration when given a transforms table

### Under Construction ####

# resizeImagesFIJI = function(fiji = fiji()) {}

#TODO - install FIJI jython registration and resize scripts
# install_FIJI_scripts = function(fiji = fiji()) {}

#TODO These things require a correct set of boundary values
# - Subset images in Giotto using Magick and followup reassignment as the default 'image'
# - Follow this up with potential registration
# - Need a way to determine the pixel distances between spots to get an idea of which regions of image 'belong' to a spot
# - Would be nice to be able to put together an image mask even in magick and apply it to the image to aid with img_reg and take care of jagged lines after image subsetting
# - A shiny app to subset tissue regions would be nice
# The shiny app should be able to select spots in a 2d plane by default
# If given the ability, it should also select spots of a single plane or within a certain range of z values and plot them as a 2D for selection purposes




