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
#' @param spatlocs spatial locations to use
#' @param rotateradians radians by which the x and y values will be rotated in a counter clockwise manner
#' @keywords internal
rotate_spatial_locations <- function(spatlocs,
                                     rotateradians) {

  xvals = spatlocs$sdimx
  yvals = spatlocs$sdimy

  spatlocs$sdimx <- xvals*cos(rotateradians) + yvals*sin(rotateradians)
  spatlocs$sdimy <- -xvals*sin(rotateradians) + yvals*cos(rotateradians)

  return(spatlocs)
}

#' @name xy_translate_spatial_locations
#' @description Translate given X Y coordinates by given x and y translation values
#' @param spatlocs spatial locations to use
#' @param xtranslate value to translate coordinates in the positive x direction
#' @param ytranslate value to translate coordinates in the positive y direction
#' @keywords internal
xy_translate_spatial_locations <- function(spatlocs,
                                           xtranslate,
                                           ytranslate) {

  spatlocs$sdimx = spatlocs$sdimx + xtranslate
  spatlocs$sdimy = spatlocs$sdimy + ytranslate

  return(spatlocs)
}



#' @name rigid_transform_spatial_locations
#' @description Performs appropriate transforms to align spatial locations with registered images.
#' @param spatlocs input spatial locations
#' @param transform_values transformation values to use
#' @keywords internal
#Rotation is performed first, followed by XY transform.
rigid_transform_spatial_locations <- function(spatlocs,
                                              transform_values) {

  spatlocsXY <- spatlocs[,c('sdimx','sdimy')]
  #These functions must be performed in positive y values
  spatlocsXY$sdimy <- -spatlocsXY$sdimy

  spatlocsXY <- rotate_spatial_locations(spatlocsXY,
                                         transform_values$Theta)

  spatlocsXY <- xy_translate_spatial_locations(spatlocsXY,
                                               transform_values$XFinalTransform,
                                               transform_values$YFinalTransform)

  spatlocs$sdimx <- spatlocsXY$sdimx
  spatlocs$sdimy <- -spatlocsXY$sdimy

  return(spatlocs)
}

#' @name reg_img_minmax_finder
#' @description finds new minmax boundaries of registration transformed images
#' @param gobject_list list of gobjects to use
#' @param image image slot to use as original unregistered images
#' @param scalefactor scalefactor for registered images relative to originals. Defaults to 1
#' @param transform_values transformation values to use
#' @keywords internal
#Automatically account for changes in image size due to alignment
reg_img_minmax_finder <- function(gobject_list,
                                  unreg_image_slot,
                                  scalefactor = 1,
                                  transform_values) {

  #Find image spatial info from original image if possible
  #Check to make sure that the unreg_image_slot finds an existing image in each gobject to be registered
  imgPresent <- function(gobject,image) {
    imgPresent = FALSE
    if(image %in% showGiottoImageNames(gobject = gobject, verbose = FALSE)) {
      imgPresent = TRUE
    }
    return(imgPresent)
  }

  if(all(as.logical(lapply(gobject_list, imgPresent, image = unreg_image_slot)))) {

    giottoImage_list = lapply(gobject_list, getGiottoImage, image_name = unreg_image_slot)
    image_corners = lapply(giottoImage_list, get_img_corners)

    #scale the corner coords
    image_corners = lapply(image_corners, scale_spatial_locations, scalefactor = scalefactor)

    #register corners based on transform values
    image_corners_reg = lapply(1:length(image_corners), FUN = function(x) {
      rigid_transform_spatial_locations(spatlocs = image_corners[[x]],
                                        transform_values = transform_values[[x]])
    })

    image_corners_reg = dplyr::bind_rows(image_corners_reg)
    minmaxRegVals = list('xmax_reg' = max(image_corners_reg$sdimx),
                         'xmin_reg' = min(image_corners_reg$sdimx),
                         'ymax_reg' = max(image_corners_reg$sdimy),
                         'ymin_reg' = min(image_corners_reg$sdimy))

    #return the max values
    return(minmaxRegVals)
  } else {
    warning('original images must be supplied for registered images to be aligned')
  }
}


#' @name get_img_corners
#' @description finds four corner spatial coords of giottoImages or magick-images
#' @param img_object giottoImage or magick-image to use
#' @keywords internal
get_img_corners <- function(img_object) {
  if(methods::is(img_object,'giottoImage')) {
    img_dims = Giotto:::get_img_minmax(img_object@mg_object)
  } else if(methods::is(img_object,'magick-image')) {
    img_dims = Giotto:::get_img_minmax(img_object)
  } else {
    stop('img_object must be either a giottoImage or a magick-image /n')
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


#' @name registerGiottoObjectList
#' @description Wrapper function for registerGiottoObjectListFiji and registerGiottoObjectListRvision
#' @param gobject_list list of gobjects to register
#' @param spat_loc_values spatial locations to use
#' @param spat_loc_name arbitrary name for registered spatial locations. Defaults to replacement of spat_loc_values (optional)
#' @param fiji_xml_files record of transformations performed in xml format exported during image registration using FIJI register_virtual_stack_slices
#' @param fiji_registered_images registered images output by FIJI register_virtual_stack_slices
#' @param registered_image_name arbitrary slot name for registered images. Defaults to the slot of unregistered images.
#' @param scaling scaling to be applied to spatial coordinates
#' @param method method used to align gobjects. Current options are either using FIJI register_virtual_stack_slices output or rvision
#' @param allow_rvision_autoscale Whether or not to allow rvision to automatically scale the images when performing image registration
#' @param verbose be verbose
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
                                                image = image,
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
#' @param image slot of original unregistered images
#' @param spat_loc_values spatial locations to use
#' @param spat_loc_name slot to save spatial locations to. Defaults to replacement of spat_loc_values (optional)
#' @param xml_files record of transformations performed in xml format exported during image registration using FIJI register_virtual_stack_slices
#' @param registered_images registered images output by FIJI register_virtual_stack_slices
#' @param scaling scaling to be applied to spatial coordinates
#' @param verbose be verbose
#' @return list of registered giotto objects where the registered images and spatial locations
#' @export
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
    spatloc = Giotto:::get_spatial_locations(gobject = gobj, ###Tag for editing later
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
  transformsDF = lapply(transf_list,
                        FUN = trakem2_rigid_transf_extract)

  ##Transformation of spatial coordinates
  #Start with scaling of the spatial coordinates
  spatloc_list = lapply(spatloc_list,
                        FUN = scale_spatial_locations,
                        scalefactor = scaling)

  spatloc_list = lapply(1:length(spatloc_list),
                        FUN = function(x) {
                          rigid_transform_spatial_locations(spatlocs = spatloc_list[[x]],
                                                            transform_values = transformsDF[[x]])
                        })

  ## 4. update Giotto slots and names and return list of Giotto object

  #Determine if registered images are present and if there is one for each gobject
  if(is.null(registered_images) == FALSE) {
    #If there are not the same number of registered images as gobjects, stop
    if(length(registered_images) != length(gobject_list)) {
      stop('A registered image should be supplied for every gobject to align /n')
    }
    if(sum(as.logical(lapply(registered_images, methods::is, class2 = 'giottoImage'))) > 0)
      stop('Registered images should be supplied as either magick-objects or filepaths /n')
  }

  #Find new image boundaries for registered images
  #Must have original pre-registration images in the gobject for this to work
  #TODO (optional if just registering spatlocs)
  reg_img_boundaries = reg_img_minmax_finder(gobject_list = gobject_list,
                                             unreg_image_slot = image,
                                             scalefactor = scaling,
                                             transform_values = transformsDF)

  #GOBJECT UPDATING FOR LOOP
  for(gobj_i in 1:length(gobject_list)) {
    gobj = gobject_list[[gobj_i]]


    #SPATIAL UPDATES
    #Assign original spatial locations to 'unregistered' slot
    gobj@spatial_locs$unregistered = gobj@spatial_locs[[spat_loc_values]]

    #Assign registered spatial locations from spatloc_list to gobject_list
    gobj@spatial_locs[[spat_loc_name]] = spatloc_list[[gobj_i]]


    #IMAGE UPDATES
    #If there is an image in the image slot and the intended slot of the registered image is 'image'
    #Move the original image to 'unregistered'
    if('image' %in% showGiottoImageNames(gobj, verbose = FALSE) && registered_image_name == 'image') {
      if (verbose == TRUE) {
        cat(paste0('Gobject', gobj_i, ': original image in "image" slot moved to "unregistered" slot. \n'))
      }

      gobj@images$unregistered = gobj@images[[image]]

    }

    #Create a giotto image if there are registered images supplied
    if(!is.null(registered_images)) {
      g_image = createGiottoImage(gobject = gobj,
                                  spatial_locs = spat_loc_name,
                                  mg_object = registered_images[gobj_i],
                                  name = registered_image_name)

      #Add the registered image to the gobj under registered_image_name slot.
      gobj = addGiottoImage(gobject = gobj, images = list(g_image))

      #Automatic adjustment
      if(exists('reg_img_boundaries')){ #TODO
        im_info = gobj@images[[image]]@minmax

        #update the giottoImage boundaries
        boundaries <- as.numeric(c(reg_img_boundaries$xmax_reg - im_info[['xmax_sloc']],
                                   im_info[['xmin_sloc']] - reg_img_boundaries$xmin_reg,
                                   reg_img_boundaries$ymax_reg - im_info[['ymax_sloc']],
                                   im_info[['ymin_sloc']] - reg_img_boundaries$ymin_reg))

        names(boundaries) = c('xmax_adj','xmin_adj','ymax_adj','ymin_adj')

        gobj@images[[registered_image_name]]@boundaries = boundaries

      }
    }
    gobject_list[[gobj_i]] = gobj
  }
  return(gobject_list)
}

#TODO check if spatloc is actually provided in createGiottoImage() and ignore auto align if not.


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
#' @param source Folder containing images to be registered
#' @param output Folder to save registered images to
#' @param transformSave Folder to save transforms to
#' @param referenceImg File name of reference image for the registration
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
#' @param fijiPath Path to fiji executable (can be set by
#'   \code{options(giotto.fiji="/some/path")})
#' @param DryRun Whether to return the command to be run rather than actually
#'   executing it.
#' @return list of registered giotto objects where the registered images and spatial locations
#' \dontrun{
#' #This function was adapted from runFijiMacro function in jimpipeline by jefferislab #
#' }
#' @export
registerImagesFIJI = function(source,
                              output,
                              #Scale Invariant Interest Point Detector Options
                              init_gauss_blur = 1.6,
                              steps_per_scale_octave = 3,
                              min_img_size = 64,
                              max_img_size = 1024,
                              #Feature Descriptor Options
                              feat_desc_size = 8,
                              feat_desc_orient_bins = 8,
                              closest_next_closest_Ratio = 0.92,
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
                              fijiPath = fiji(),
                              DryRun = FALSE) {

  if(headless) fijiArgs = c(fijiArgs,"--headless")
  fijiArgs=paste(fijiArgs,collapse=" ")

  javaArgs=c(paste("-Xms",MinMem,'m',sep=""),paste("-Xmx",MaxMem,'m',sep=""),javaArgs)
  if(IncrementalGC) javaArgs=c(javaArgs,"-Xincgc")
  javaArgs=paste(javaArgs,collapse=" ")

  threadAdjust=ifelse(is.null(Threads),"",paste("run(\"Memory & Threads...\", \"parallel=",Threads,"\");",sep=""))

  macroCall=paste(" -eval '",
                  threadAdjust,
                  "run(\"Register Virtual Stack Slices\", \"source=[",
                  source,
                  "] output=[",
                  output,
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
                  closest_next_closest_Ratio,
                  " maximal_alignment_error=",
                  max_align_err,
                  " inlier_ratio=",
                  inlier_ratio,
                  " feature_extraction_model=Similarity registration_model=[Rigid                -- translate + rotate                  ] interpolate\");' ",sep="")

  ijArgs=paste(c(ijArgs,ifelse(batch,"-batch","")),collapse=" ")

  cmd<-paste(fijiPath,javaArgs,fijiArgs,"--",macroCall,ijArgs)
  if(DryRun) return(cmd)
  return(0==system(cmd))
}
















