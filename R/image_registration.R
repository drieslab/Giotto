### registering giotto object ####

### helper functions ####

#' @name trakem2_rigid_transf_extract
#' @description Extract rigid registration transformation values from FIJI TrakEM2 xml file. Generated through register_virtual_stack_slices.
#' @param inputstring string read in from TrakeEM2 xml file
#' @keywords internal
trakem2_rigid_transf_extract <- function(inputstring) {

  #Catch wrong inputs
  if (grepl('^.*trakem2.*', inputstring, ignore.case = T) != 1) {
    stop('xml files must be in TrakeEM2 format')
  }


  #Regex to find the values from the TrakEM2 .xml
  transfExtractPatA = "(?<=iict_transform class=\"mpicbg.trakem2.transform.RigidModel2D\" data=\")(.*)(?=\")"
  transfExtractPatB = "(?<=class=\"mpicbg.trakem2.transform.TranslationModel2D\" data=\")(.*)(?=\")"
  #Split relevant text into numerical values
  out <- c(sapply(stringr::str_split(stringr::str_extract(string = inputstring, pattern = transfExtractPatA), pattern = ' '), function(x) as.numeric(x)),
           sapply(stringr::str_split(stringr::str_extract(string = inputstring, pattern = transfExtractPatB), pattern = ' '), function(x) as.numeric(x)))

  if(sum(is.na(out)) == 2) {
    out <- rep(0,5)
  }

  out <- c(out,0,0)
  out <- data.table::data.table(t(matrix(out)))
  colnames(out) <- c('Theta','Xtransform','Ytransform','itx','ity','XFinalTransform','YFinalTransform')

  #itx and ity are additional values in the trakem2 xml files that must be added to Xtransform and Ytransform in order to get the final transformation values.
  #only relevant for sampleset with more than 1 slice away from the reference image
  out$XFinalTransform <- out$Xtransform + out$itx
  out$YFinalTransform <- out$Ytransform + out$ity

  #Multiply theta by -1 due to differences in R and image plotting coordinates
  out$Theta <- -out$Theta

  return(out)
}



#' @name scale_spatial_locations
#' @description Scale the X and Y coordinates given by the scale factor input
#' @param spatlocs spatial locations to scale
#' @param scalefactor scaling factor to apply to coordinates
#' @keywords internal
scale_spatial_locations <- function(spatlocs,
                                    scalefactor) {

  spatlocs$sdimx <- spatlocs$sdimx*scalefactor
  spatlocs$sdimy <- spatlocs$sdimy*scalefactor

  return(spatlocs)
}



#' @name rotate_spatial_locations
#' @description Rotate given X Y coordinates by given radians in counter clockwise manner about the coordinate origin
#' @param xvals x coordinates
#' @param yvals y coordinates
#' @param rotateradians radians by which the x and y values will be rotated in a counter clockwise manner
#' @keywords internal
rotate_spatial_locations <- function(xvals,
                           yvals,
                           rotateradians) {

  xprime <- xvals*cos(rotateradians) + yvals*sin(rotateradians)
  yprime <- -xvals*sin(rotateradians) + yvals*cos(rotateradians)

  return(cbind(xprime,yprime))
}

#' @name xy_translate_spatial_locations
#' @description Translate given X Y coordinates by given x and y translation values
#' @param xvals x coordinates
#' @param yvals y coordinates
#' @param xtranslate value to translate coordinates in the positive x direction
#' @param ytranslate value to translate coordinates in the positive y direction
#' @keywords internal
xy_translate_spatial_locations <- function(xvals,
                                           yvals,
                                           xtranslate,
                                           ytranslate) {

  xprime <- xvals + xtranslate
  yprime <- yvals + ytranslate

  return(cbind(xprime,yprime))
}



#' @name rigid_transform_spatial_locations
#' @description Performs appropriate transforms to align spatial locations with registered images.
#' @param spatlocs input spatial locations
#' @param transform_values transformation values to use
#' @keywords internal
#Rotation is performed first, followed by XY transform
rigid_transform_spatial_locations <- function(spatlocs,
                                              transform_values) {

  spatLocsXY <- cbind(spatlocs$sdimx, spatlocs$sdimy) #TODO Check these
  spatLocsXY <- rotate_spatial_locations(spatLocsXY[,1],spatLocsXY[,2],transform_values$Theta)
  spatLocsXY <- xy_translate_spatial_locations(spatLocsXY[,1],spatLocsXY[,2],transform_values$XFinalTransform,transform_values$YFinalTransform)

  spatlocs$sdimx <- spatLocsXY[,1]
  spatlocs$sdimy <- spatLocsXY[,2]

  return(spatlocs)
}

#' @name auto_comp_reg_border
#' @description adjusts for increase in gobject image extent from transformations performed during registration
#' @param gobject gobject to use
#' @param transform_values transformation values to use
#' @keywords internal
#Automatically account for changes in image size due to alignment
auto_comp_reg_border <- function(gobject,transform_values) {

  #Find image spatial info from original image if possible
  if(image %in% showGiotoImageNames(gobj, verbose = FALSE)) {
    im_info <- list(gobj@images[[image]]@boundaries,
                    gobj@images[[image]]@minmax)


  }
}




### exported functions ####


#' @name registerGiottoObjectList
#' @description Wrapper function for registerGiottoObjectListFiji and registerGiottoObjectListRvision
#' @param gobject_list list of gobjects to register
#' @param spat_loc_values spatial locations to use
#' @param spat_loc_name arbitrary name for registered spatial locations. Defaults to replacement of spat_loc_values (optional)
#' @param fiji_xml_files record of transformations performed in xml format exported during image registration using FIJI register_virtual_stack_slices
#' @param fiji_registered_images registered images output by FIJI register_virtual_stack_slices
#' @param registered_image_name arbitrary name for registered images. Defaults to replacement of unregistered images
#' @param scaling scaling to be applied to spatial coordinates
#' @param method method used to align gobjects. Current options are either using FIJI register_virtual_stack_slices output or rvision
#' @param allow_rvision_autoscale Whether or not to allow rvision to automatically scale the images when performing image registration
#' @param verbose
#' @return list of registered giotto objects where the registered images and spatial locations
#' @export
registerGiottoObjectList <- function(gobject_list,
                                     image = 'image',
                                     registered_image_name = 'image',
                                     spat_loc_values = NULL, #?
                                     spat_loc_name = 'raw',
                                     method = c('fiji','rvision'),
                                     fiji_xml_files,
                                     fiji_registered_images,
                                     scaling = 1,
                                     allow_rvision_autoscale = TRUE,
                                     auto_comp_reg_border = TRUE,
                                     verbose = TRUE) {

  if(method == 'fiji') {
    gobject_list = registerGiottoObjectListFiji(gobject_list = gobject_list,
                                                images = image,
                                                registered_image_name = registered_image_name,
                                                spat_loc_values = spat_loc_values,
                                                spat_loc_name = spat_loc_name,
                                                xml_files = fiji_xml_files,
                                                registered_images = fiji_registered_images,
                                                registered_image_name = registered_image_name,
                                                scaling = scaling,
                                                auto_comp_reg_border = auto_comp_reg_border,
                                                verbose = verbose)

  ##TODO:
  # rvision registration not yet implemented

  # } else if (method == 'rvision') {
  #   gobject_list = registerGiottoObjectListRvision(gobject_list = gobject_list,
  #                                                  images = images,
  #                                                  registered_image_name = registered_image_name,
  #                                                  spat_loc_values = spat_loc_values,
  #                                                  spat_loc_name = spat_loc_name,
  #                                                  scaling = scaling,
  #                                                  allow_rvision_autoscale = allow_rvision_autoscale,
  #                                                  auto_comp_reg_border = auto_comp_reg_border,
  #                                                  verbose = verbose)

  } else {
      stop('Only FIJI registration is currently supported')
  }

  return(gobject_list)
}



#' @name registerGiottoObjectListFiji
#' @description Function to spatially align gobject data based on FIJI image registration.
#' @param gobject_list list of gobjects to register
#' @param images slot of original unregistered images (optional)
#' @param spat_loc_values spatial locations to use
#' @param spat_loc_name slot to save spatial locations to. Defaults to replacement of spat_loc_values (optional)
#' @param xml_files record of transformations performed in xml format exported during image registration using FIJI register_virtual_stack_slices
#' @param registered_images registered images output by FIJI register_virtual_stack_slices
#' @param scaling scaling to be applied to spatial coordinates
#' @param verbose
#' @return list of registered giotto objects where the registered images and spatial locations
#' @export
#Register giotto objects when given xml transform output files from FIJI
registerGiottoObjectListFiji = function(gobject_list,
                                        image = 'image',
                                        registered_image_name = 'image',
                                        spat_loc_values = NULL,
                                        spat_loc_name = 'raw',
                                        xml_files,
                                        registered_images = NULL,
                                        scaling = 1,
                                        auto_comp_reg_border,
                                        verbose = TRUE) {

  ## 1. get spatial coordinates and put in list ##
  spatloc_list = list()
  for(gobj_i in 1:length(gobject_list)) {
    gobj = gobject_list[[gobj_i]]
    spatloc = Giotto:::select_spatial_locations(gobject = gobj, ###Tag for editing later
                                       spat_loc_name = spat_loc_name)
    #------ Put all spatial location data together
    spatloc_list[[gobj_i]] = spatloc
  }

  ## 2. read transform xml files ##
  transf_list = list()
  for(file_i in 1:length(xml_files)) {
    t_file = xml_files[[file_i]]
    #------ Put all your transform files together
    transf_list[[file_i]] = paste(readLines(t_file, warn = F), collapse = '\n')
  }

  ## 3. apply transformation on spatial locations ##
  if(length(spatloc_list) != length(transf_list)) {
    stop('xml spatial transforms must be supplied for every gobject to be registered.')
  }

  #Select useful info out of the TrakEM2 files
  transformsDF <- lapply(transf_list,
                         FUN = trakem2_rigid_transf_extract)

  ##Transformation of spatial coordinates
  #Start with scaling of the spatial coordinates
  spatloc_list <- lapply(spatloc_list,
                         FUN = scale_spatial_locations,
                         scalefactor = scaling)

  spatloc_list <- lapply(1:length(spatloc_list),
                         FUN = function(x) {rigid_transform_spatial_locations(spatlocs = spatloc_list[[x]],
                                                                              transform_values = transformsDF[[x]])
                           })

  ## 4. update Giotto slots and names and return list of Giotto object

  #If there is a list of paths to image files in registered images, place them in 'image' slot
  if(is.null(registered_images) == FALSE) {
    #If there are not the same number of registered images as gobjects, stop
    if (length(registered_images) != length(gobject_list)) {
      stop('A registered image should be supplied for every gobject to align')
    }

  }

  for(gobj_i in 1:length(gobject_list)) {
    gobj = gobject_list[[gobj_i]]
    #Assign original spatial locations to 'unregistered' slot
    gobj@spatial_locs$unregistered <- gobj@spatial_locs$spat_loc_values

    #Assign registered spatial locations from spatloc_list to gobject_list
    gobj@spatial_locs$spat_loc_values <- spatloc_list[[gobj_i]]





    #If there is an image in the image slot and the intended slot of the registered image is 'image'
    #Move the original image to 'unregistered'
    if('image' %in% showGiottoImageNames(gobj, verbose = FALSE) && registered_image_name == 'image') {
      if (verbose == TRUE) {
        cat(paste0('Gobject', gobj_i, ': original image in "image" moved to "unregistered" slot. \n'))
      }

      gobj@images$unregistered <- gobj@images$image

    }




    #Create a giotto image if there are registered images supplied
    if(registered_images != NULL) {
      g_image <- createGiottoImage(registered_images[gobj_i])


      #TODO Generate image adjustments if possible
      #Check if either 'image' (default) or other supplied unregistered image slot supplied through image exists in the gobject.

    }




    gobj@images[[registered_image_name]]



  }

  return(gobject_list)
}






#TEMPLATE

#Only necessary to change scaleSpatialDataFactor from 1 if the pixel xy coordinates in the spatlocs object require a scaling operation to match that of the images used in registration.
registerSpatialLocations <- function(spatLocs, registerTransforms, scale_spatloc_factor = 1) {
  if(length(spatLocs) != length(registerTransforms)) {
    stop('A set of spatial transforms must be supplied for every set of spatial locations to be registered.')
  }

  #Select useful info out of the TrakEM2 files
  transfDF <- lapply(registerTransforms, trakem2_rigid_transf_extract)

  ##Transformation of spatial coordinates
  #Start with scaling of the spatial coordinates
  spatLocs <- lapply(spatLocs, scale_spatial_locations, scalefactor = scale_spatloc_factor)

  spatLocs <- lapply(1:length(spatLocs), FUN = function(x) {rigid_transform_spatial_locations(spatlocs = spatLocs[[x]],transform_values = transfDF[[x]])})

  return(spatLocs)
}

#TODO edit rigid_transform_spatial_locations (gobj uses xdim and ydim)
#TODO check on sdimx and sdimy functions (rigid_transform_spatial_locations and scale_spatial_locations)









