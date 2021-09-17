 
### Image reg functions ####

#' @name trakem2_rigid_transf_extract
#' @title trakem2_rigid_transf_extract
#' @description Extract rigid registration transformation values from FIJI TrakEM2 xml file. Generated through register_virtual_stack_slices.
#' @param inputstring string read in from TrakeEM2 xml file
#' @keywords internal
trakem2_rigid_transf_extract <- function(inputstring) {
  
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
#' @title scale_spatial_locations
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
#' @title rotate_spatial_locations
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
#' @title xy_translate_spatial_locations
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



#' @name register_spatial_locations
#' @title register_spatial_locations
#' @description Performs appropriate transforms to align spatial locations with registered images.
#' @param xvals x value spatial input
#' @param yvals y value spatial input
#' @param scalefactor scale spatial coords to pixel coords
#' @param transformXML transformation files to use
#' @keywords internal
#Rotation is performed first, followed by XY transform.
register_spatial_locations <- function(xvals,
                                     yvals,
                                     scalefactor = 1,
                                     transformXML) {
  
  spatlocsXY = cbind(xvals,yvals)
  spatlocsXY = as.data.frame(spatlocsXY)
  colnames(spatlocsXY) = c('sdimx','sdimy')
  
  #Generate transform values from xml file
  transform_values = trakem2_rigid_transf_extract(inputstring = transformXML)
  
  #These functions must be performed in positive y values
  spatlocsXY$sdimy = -spatlocsXY$sdimy
  spatlocsXY = scale_spatial_locations(spatlocs = spatlocsXY,scalefactor = scalefactor)
  
  spatlocsXY = rotate_spatial_locations(spatlocsXY,
                                        transform_values$Theta)
  
  spatlocsXY = xy_translate_spatial_locations(spatlocsXY,
                                              transform_values$XFinalTransform,
                                              transform_values$YFinalTransform)
  
  spatlocsXY = scale_spatial_locations(spatlocs = spatlocsXY,scalefactor = (1/scalefactor))
  spatlocsXY$sdimy = -spatlocsXY$sdimy
  
  return(spatlocsXY)
  
}


### Data Prep for zStack Gobject ####

#' @name join_expression
#' @title join_expression
#' @description joins expression list together while appending z value to cellIDs to ensure unique colnames
#' @param expression_list list of expression values to merge
#' @param z_vals z values to use z stacking expression values
#' @keywords internal
join_expression = function(expression_list,
                          z_vals) {
  
  ## 1. Make each slice's 2D ID unique in 3D
  idList = lapply(expression_list, colnames)
  for(i in 1:length(idList)) {
    idList[[i]] = paste0(idList[[i]],'x',z_vals[[i]])
  }
  for(i in 1:length(expression_list)) {
    colnames(expression_list[[i]]) = idList[[i]]
  }
  
  ## 2. Find all genenames (rownames) and add them as rows to expr sets that are missing them
  feats = lapply(expression_list, rownames)
  feats = Reduce(f = c, feats)
  feats = unique(feats)
  
  # Find the features that are missing from individual count matrices but present in others of the z stack
  notin = lapply(1:length(expression_list), FUN = function(x) {
    feats[!(feats %in% rownames(expression_list[[x]]))]
  })
  
  # Append extra rows so that all expression matrices have the same dimensions
  for(i in 1:length(expression_list)) {
    expression_list[[i]] = rbind(expression_list[[i]], Matrix::Matrix(data = 0,
                                                                      nrow = length(notin[[i]]),
                                                                      ncol = ncol(expression_list[[i]]),
                                                                      dimnames = list(c(notin[[i]]),colnames(expression_list[[i]]))))
  }
  
  ## 3. Reorder all count matrices
  for(i in 1:length(expression_list)) {
    expression_list[[i]] = expression_list[[i]][order(rownames(expression_list[[i]])),]
  }
  
  ## 4. Merge expression
  expr_merge = Reduce(f = cbind, expression_list)
  
  ## 5. Return the combined count matrix
  return(expr_merge)
}


#' @name stack_spatlocs
#' @title stack_spatlocs
#' @description Adds z values and combines list of spatlocs into single z-stacked spatloc table
#' @param spatlocs_list list of spatial locations
#' @param z_vals z values to use for z stacking spatial locations list
#' @keywords internal
stack_spatlocs = function(spatlocs_list,
                         z_vals) {
  
  # 1. Append z value based on z_vals
  for(i in 1:length(spatlocs_list)) {
    spatlocs_list[[i]]$sdimz = z_vals[[i]]
  }
  
  # 2. rbind together all spatlocs
  spatloc_merge = Reduce(f = rbind, spatlocs_list)
  
  return(spatloc_merge)
}  

#' @name createRegZStackGobject
#' @title createRegZStackGobject
#' @description create registered z stacked 3D gobject. Images not supported
#' @param expression_list list of expression values to use
#' @param spatlocs_list list of spatlocs to use
#' @param xvals x value spatial input
#' @param yvals y value spatial input
#' @param scalefactor scale spatial coords to pixel coords
#' @param transformXML transformation files to use
#' @param z_vals z values to use in z stacking. In order of expression_list and spatloc_list.
#' @param instructions instructions for Giotto processing
#' @export
createRegZStackGobject = function(expression_list,
                                  spatlocs_list,
                                  xvals,
                                  yvals,
                                  scalefactor = 1,
                                  transformXML,
                                  z_vals,
                                  instructions = NULL) {
  
  # 1. Register spatlocs
  regSpatlocs = lapply(1:length(spatlocs_list), FUN = function(x) {
    register_spatial_locations(xvals = spatlocs_list[[x]][xvals],
                               yvals = spatlocs_list[[x]][yvals],
                               scalefactor = scalefactor,
                               transformXML = transformXML[[x]])
  })
  
  # 2. combine spatlocs
  stackRegSpatLocs = stack_spatlocs(spatlocs_list = regSpatlocs,
                                   z_vals = z_vals)
  
  # 3. combine expression
  expr_merge = join_expression(expression_list = expression_list,
                              z_vals = z_vals)
  
  # 4. create Giotto object
  zStackGobject = createGiottoObject(raw_exprs = expr_merge,
                                     spatial_locs = stackRegSpatLocs,
                                     instructions = instructions)
  
  # 5. return Giotto object
  return(zStackGobject)
}







