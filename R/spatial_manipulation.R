
#' @include generics.R
#' @importMethodsFrom terra spin
#' @importMethodsFrom terra flip
#' @importMethodsFrom Matrix t
#' @importMethodsFrom terra t
#' @importMethodsFrom terra ext
NULL

# Spatially manipulate objects ####




## flip ####
#' @title Flip an object
#' @name flip-generic
#' @description Flip an object over a designated x or y value depending on
#' direction param input. Note that this behavior may be different from terra's
#' @param x object
#' @param direction character. Direction to flip. Should be either partial match to 'vertical' or 'horizontal'
#' @param x0 x value to flip horizontally over (ignored for vertical). Pass NULL
#' to flip over the extent
#' @param y0 y value to flip vertically over (ignored for horizontal). Pass NULL
#' to flip over the extent
#' @param ... additional args to pass
NULL

#' @describeIn flip-generic Flip a giottoPolygon object
#' @export
setMethod('flip', signature(x = 'giottoPolygon'),
          function(x, direction = 'vertical', x0 = 0, y0 = 0, ...) {
            flip_gpoly(gpoly = x, direction = direction, x0 = x0, y0 = y0)
          })

#' @describeIn flip-generic Flip a giottoPoints object
#' @export
setMethod('flip', signature(x = 'giottoPoints'),
          function(x, direction = 'vertical', x0 = 0, y0 = 0, ...) {
            flip_gpoints(gpoints = x, direction = direction, x0 = x0, y0 = y0)
          })

#' @describeIn flip-generic Flip a spatLocsObj
#' @export
setMethod('flip', signature(x = 'spatLocsObj'),
          function(x, direction = 'vertical', x0 = 0, y0 = 0, ...) {
            flip_spatlocs(sl = x, direction = direction, x0 = x0, y0 = y0)
          })

#' @describeIn flip-generic Flip a spatialNetworkObj
#' @export
setMethod('flip', signature(x = 'spatialNetworkObj'),
          function(x, direction = 'vertical', x0 = 0, y0 = 0, ...) {
            flip_spatnet(sn = x, direction = direction, x0 = x0, y0 = y0)
          })

# TODO apply as instructions for lazy eval after crop/resampling
#' @describeIn flip-generic Flip a giottoLargeImage
#' @export
setMethod('flip', signature(x = 'giottoLargeImage'),
          function(x, direction = 'vertical', ...) {
            x@raster_object = terra::flip(x@raster_object)
            x
          })

#' @describeIn flip-generic Flip a SpatExtent
#' @export
setMethod('flip', signature(x = 'SpatExtent'),
          function(x, direction = 'vertical', x0 = 0, y0 = 0) {
            flip_extent(e = x, direction = direction, x0 = x0, y0 = y0)
          })



## spin ####

#' @title Spin an object
#' @name spin-generic
#' @description Spin (rotate) an object spatially (limited to xy rotations)
#' @param x object
#' @param angle numeric. Angle of rotation in degrees
#' @param x0 numeric. x-coordinate of the center of rotation. Defaults to center x val if not given.
#' @param y0 numeric. y-coordinate of the center of rotation. Defaults to center y val if not given.
NULL

#' @describeIn spin-generic Spin a giottoPolygon object
#' @export
setMethod('spin', signature(x = 'giottoPolygon'),
          function(x, angle, x0 = NULL, y0 = NULL) {
            if(is.null(x0)) x0 = terra::mean(terra::ext(x@spatVector))[1]
            if(is.null(y0)) y0 = terra::mean(terra::ext(x@spatVector))[2]
            return(do_gpoly(x = x, what = 'terra'::'spin', args = list(angle = angle, x0 = x0, y0 = y0)))
            # x@spatVector = terra::spin(x = x@spatVector,
            #                            angle = angle,
            #                            x0 = x0,
            #                            y0 = y0)
            # if(!is.null(x@spatVectorCentroids)) {
            #   x@spatVectorCentroids = terra::spin(x = x@spatVectorCentroids,
            #                                       angle = angle,
            #                                       x0 = x0,
            #                                       y0 = y0)
            # }
            # if(!is.null(x@overlaps)) {
            #   lapply(x@overlaps, function(overlap) {
            #     if(inherits(overlap, 'SpatVector')) {
            #       terra::spin(x = overlap,
            #                   angle = angle,
            #                   x0 = x0,
            #                   y0 = y0)
            #     } else {
            #       overlap
            #     }
            #   })
            # }
            # return(x)
          })

#' @describeIn spin-generic Spin a giottoPoints object
#' @export
setMethod('spin', signature(x = 'giottoPoints'),
          function(x, angle, x0 = NULL, y0 = NULL) {
            if(is.null(x0)) x0 = terra::mean(terra::ext(x@spatVector))[1]
            if(is.null(y0)) y0 = terra::mean(terra::ext(x@spatVector))[2]
            x@spatVector = terra::spin(x@spatVector,
                                       angle = angle,
                                       x0 = x0,
                                       y0 = y0)
            return(x)
          })

#' @describeIn spin-generic Spin a spatLocsObj
#' @param z0 spatLocsObj specific. Numeric. z-coordinate of the center of rotation.
#' Depending on if z data is present, defaults to either 0 or center z val if not given.
#' @param xy_angle spatLocsObj specific. xy plane rotation in degrees.
#' Overrides angle param
#' @param zy_angle spatLocsObj specific. zy plane rotation
#' @param xz_angle spatLocsObj specific. xz plane rotation
#' @export
setMethod('spin', signature(x = 'spatLocsObj'),
          function(x, angle = NULL, x0 = NULL, y0 = NULL, z0 = NULL,
                   xy_angle = NULL, zy_angle = NULL, xz_angle = NULL) {
            if(!is.null(angle)) xy_angle = angle
            if(is.null(xy_angle)) xy_angle = 0
            if(is.null(zy_angle)) zy_angle = 0
            if(is.null(xz_angle)) xz_angle = 0
            angles = c(xy = xy_angle, zy = zy_angle, xz = xz_angle)
            angles = radians(angles)

            if(is.null(x0)) x0 = mean(c(min(x[]$sdimx), max(x[]$sdimx)))
            if(is.null(y0)) y0 = mean(c(min(x[]$sdimy), max(x[]$sdimy)))
            if('sdimz' %in% colnames(x[])) {
              if(is.null(z0)) z0 = mean(c(min(x[]$sdimz), max(x[]$sdimz)))
              else z0 = 0
            }
            x[] = rotate_spatial_locations(spatlocs = x[],
                                           rotateradians = angles,
                                           rcenter = c(x = x0, y = y0, z = z0))
            return(x)
          })




## spatShift ####


#' @describeIn spatShift Shift the locations of a spatLocsObj
#' @export
setMethod('spatShift', signature('spatLocsObj'), function(x, dx = 0, dy = 0, dz = 0,
                                                          copy_obj = TRUE, ...) {
  x[] = shift_spatial_locations(spatlocs = x[], dx = dx, dy = dy, dz = dz, ...)
  return(x)
})
#' @describeIn spatShift Shift the locations of a spatialNetworkObj
#' @export
setMethod('spatShift', signature('spatialNetworkObj'), function(x, dx = 0, dy = 0, dz = 0,
                                                                copy_obj = TRUE, ...) {
  x@networkDT = shift_spatial_network(spatnet = x@networkDT,
                                      dx = dx, dy = dy, dz = dz, ...)
  if(!is.null(x@networkDT_before_filter)) {
    x@networkDT_before_filter = shift_spatial_network(spatnet = x@networkDT_before_filter,
                                                      dx = dx, dy = dy, dz = dz, ...)
  }
  return(x)
})





## transpose ####

# S4 methods
#' @title Transpose
#' @name transpose-generic
#' @param x object to be transposed
#' @aliases t
NULL
# if(!isGeneric('t')) setOldClass('t', where = as.environment("package:Giotto"))

#' @rdname transpose-generic
#' @export
setMethod('t', signature('spatLocsObj'), function(x) {
  sdimy = sdimx = NULL
  x = data.table::copy(x)
  x@coordinates[, c('sdimx', 'sdimy') := .(sdimy, sdimx)]
  return(x)
})
#' @rdname transpose-generic
#' @export
setMethod('t', signature('spatialNetworkObj'), function(x) {
  sdimx_begin = sdimx_end = sdimy_begin = sdimy_end = NULL
  x = data.table::copy(x)
  x@networkDT[, c('sdimx_begin', 'sdimy_begin', 'sdimx_end', 'sdimy_end') := .(sdimy_begin, sdimx_begin, sdimy_end, sdimx_end)]
  if(!is.null(x@networkDT_before_filter)) {
    x@networkDT_before_filter[, c('sdimx_begin', 'sdimy_begin', 'sdimx_end', 'sdimy_end') := .(sdimy_begin, sdimx_begin, sdimy_end, sdimx_end)]
  }
  return(x)
})

# s3 methods
#' @rdname transpose-generic
#' @method t spatLocsObj
#' @export
t.spatLocsObj = function(x) {
  sdimy = sdimx = NULL
  x = data.table::copy(x)
  x@coordinates[, c('sdimx', 'sdimy') := .(sdimy, sdimx)]
  return(x)
}


#' @rdname transpose-generic
#' @method t spatialNetworkObj
#' @export
t.spatialNetworkObj = function(x) {
  sdimx_begin = sdimx_end = sdimy_begin = sdimy_end = NULL
  x = data.table::copy(x)
  x@networkDT[, c('sdimx_begin', 'sdimy_begin', 'sdimx_end', 'sdimy_end') := .(sdimy_begin, sdimx_begin, sdimy_end, sdimx_end)]
  if(!is.null(x@networkDT_before_filter)) {
    x@networkDT_before_filter[, c('sdimx_begin', 'sdimy_begin', 'sdimx_end', 'sdimy_end') := .(sdimy_begin, sdimx_begin, sdimy_end, sdimx_end)]
  }
  return(x)
}






## ext ####

#' @name ext-generic
#' @title Get a SpatExtent
#' @description Get a SpatExtent of an object. This is the spatial minmax x and y
#' that the object is mapped to.
#' @param x spatial object
#' @param ... additional params to pass
NULL

#' @describeIn ext-generic Get extent of spatLocsObj
#' @export
setMethod('ext', signature('spatLocsObj'), function(x, ...) {
  sdimx = sdimy = NULL # dt vars
  terra::ext(c(range(x[][, sdimx]), range(x[][, sdimy])))
})

#' @describeIn ext-generic Get extent of giottoPolygon
#' @export
setMethod('ext', signature('giottoPolygon'), function(x, ...) {
  terra::ext(x@spatVector, ...)
})

#' @describeIn ext-generic Get extent of giottoPoints
#' @export
setMethod('ext', signature('giottoPoints'), function(x, ...) {
  terra::ext(x@spatVector, ...)
})

#' @describeIn ext-generic Get extent of spatialNetworkObj
#' @export
setMethod('ext', signature('spatialNetworkObj'), function(x, ...) {
  sdimx_begin = sdimx_end = sdimy_begin = sdimy_end = NULL # dt vars
  terra::ext(c(x[][, range(c(sdimx_begin, sdimx_end))], x[][, range(c(sdimy_begin, sdimy_end))]))
})

















# internals ####

#  TODO Following 3 functions need wrapper function for working with spatlocs and spat_info
#  Additionally add way to refer to subsets of spatial locations by list_ID


#' @title Scale spatial locations
#' @name scale_spatial_locations
#' @description Simple scaling of spatial locations by given \code{scale_factor}.
#' Values will be scaled from the coordinate origin or coordinates provided through
#' \code{scenter} param.
#' @param spatlocs spatial locations information to scale
#' @param scale_factor scaling factor to apply to coordinates.
#' @param scenter center from which to scale spatial coordinates. Given as vector
#' of xy(z) coordinates.
#' @details \code{scale_factor} either given as a single value where it will be applied to
#' x, y, and z (if available) dimensions or as a vector of named values for 'x',
#' 'y', (and 'z').
#' @keywords internal
scale_spatial_locations = function(spatlocs,
                                   scale_factor = c(x=1,y=1,z=1),
                                   scenter = c(x=0,y=0,z=0)) {

  hasZ = 'sdimz' %in% names(spatlocs)

  if(length(scale_factor) == 1) scale_factor = c(x = scale_factor, y = scale_factor, z = scale_factor)
  if(!all(names(scenter) %in% c('x','y','z'))) stop('scenter value names not recognized')
  if(!all(names(scale_factor) %in% c('x','y','z'))) stop('scale_factor value names not recognized')

  # Adjust for scaling center
  spatlocs$sdimx = spatlocs$sdimx - scenter[['x']]
  spatlocs$sdimy = spatlocs$sdimy - scenter[['y']]

  # Perform scale
  spatlocs$sdimx = spatlocs$sdimx*scale_factor[['x']]
  spatlocs$sdimy = spatlocs$sdimy*scale_factor[['y']]

  if(isTRUE(hasZ)) {
    # Adjust for scaling z center
    spatlocs$sdimz = spatlocs$sdimz - scenter[['z']]

    # Perform z scale
    spatlocs$sdimz = spatlocs$sdimx*scale_factor[['z']]

    # Revert z scaling center adjustments
    spatlocs$sdimz = spatlocs$sdimz + scenter[['z']]
  }

  # Revert scaling center adjustments
  spatlocs$sdimx = spatlocs$sdimx + scenter[['x']]
  spatlocs$sdimy = spatlocs$sdimy + scenter[['y']]

  return(spatlocs)
}


# internal for deprecation
xy_translate_spatial_locations = function(...) {
  .Deprecated(new = 'shift_spatial_locations')

  shift_spatial_locations(...)
}


#' @title Shift spatial locations
#' @name shift_spatial_locations
#' @description Shift given coordinates by given translation values
#' @param spatlocs spatial locations to use
#' @param dx value to shift coordinates in the positive x direction
#' @param dy value to shift coordinates in the positive y direction
#' @param dz value to shift coordinates in the positive z direction
#' @param xtranslate deprecated. use dx
#' @param ytranslate deprecated. use dy
#' @param ztranslate deprecated. use dz
#' @param copy_obj copy/duplicate object (default = TRUE)
#' @keywords internal
shift_spatial_locations = function(spatlocs,
                                   dx = 0,
                                   dy = 0,
                                   dz = 0,
                                   xtranslate = NULL,
                                   ytranslate = NULL,
                                   ztranslate = NULL,
                                   copy_obj = TRUE) {
  sdimx = sdimy = sdimz = NULL

  if(!is.null(xtranslate)) {
    warning(wrap_txt('xtranslate is deprecated. use dx'))
    dx = xtranslate
  }
  if(!is.null(ytranslate)) {
    warning(wrap_txt('ytranslate is deprecated. use dy'))
    dy = ytranslate
  }
  if(!is.null(ztranslate)) {
    warning(wrap_txt('ztranslate is deprecated. use dz'))
    dz = ztranslate
  }

  spatlocs[, sdimx := sdimx + dx]
  spatlocs[, sdimy := sdimy + dy]
  if('sdimz' %in% names(spatlocs)) spatlocs[, sdimz := sdimz + dz]

  return(spatlocs)
}









# See function spatShift in generics.R
#' @name shift_spatial_network
#' @title Shift spatial network
#' @description Shift spatial network coordinates
#' @param spatnet spatial network data.table
#' @param dx distance to shift on x axis
#' @param dy distance to shift on y axis
#' @param dz distance to shift on z axis
#' @param copy_obj copy/duplicate object (default = TRUE)
#' @keywords internal
shift_spatial_network = function(spatnet, dx = 0, dy = 0, dz = 0, copy_obj = TRUE) {
  sdimx_begin = sdimx_end = sdimy_begin = sdimy_end = sdimz_begin = sdimz_end = NULL

  # if 3D info present
  is3D = FALSE
  if(all(c('sdimz_begin', 'sdimz_end') %in% colnames(spatnet))) is3D = TRUE

  if(copy_obj) spatnet = data.table::copy(spatnet)

  spatnet[, `:=`(sdimx_begin = sdimx_begin + dx,
                 sdimx_end = sdimx_end + dx,
                 sdimy_begin = sdimy_begin + dy,
                 sdimy_end = sdimy_end + dy)]
  if(is3D) {
    spatnet[, `:=`(sdimz_begin = sdimz_begin + dz,
                   sdimz_end = sdimz_end + dz)]
  }
  return(spatnet)
}









#' @title Rotate spatial locations
#' @name rotate_spatial_locations
#' @description Rotate given spatlocs by given radians
#' @param spatlocs spatial locations to use
#' @param rotateradians Named vector of radians for rotation along each of the 3 coordinate
#' axes. If only a single value is provided, it will be treated as xy rotation.
#' @param rcenter center of rotation given as vector xy(z) coordinates (defaults to coordinate center)
#' @details Radians are provided through \code{rotateradians} param as a named vector
#' with values for \code{xy} (yaw), \code{zy} (pitch), \code{xz} (roll)
#' @keywords internal
rotate_spatial_locations = function(spatlocs,
                                    rotateradians = c(xy=0,zy=0,xz=0),
                                    rcenter = c(x=0,y=0,z=0)) {

  if(length(rotateradians) == 1) rotateradians = c(xy=rotateradians,zy=0,xz=0)
  if(!all(names(rotateradians) %in% c('xy','zy','xz'))) stop('rotateradians value names not recognized')
  if(!all(names(rcenter) %in% c('x','y','z'))) stop('rcenter value names not recognized')
  hasZ = 'sdimz' %in% names(spatlocs)

  # xy center of rotation adjustment
  spatlocs$sdimx = spatlocs$sdimx - rcenter[['x']]
  spatlocs$sdimy = spatlocs$sdimy - rcenter[['y']]

  xvals = spatlocs$sdimx
  yvals = spatlocs$sdimy

  # Perform rotation XY
  if(rotateradians[['xy']] != 0) {
    spatlocs$sdimx = xvals*cos(rotateradians[['xy']]) + yvals*sin(rotateradians[['xy']])
    spatlocs$sdimy = -xvals*sin(rotateradians[['xy']]) + yvals*cos(rotateradians[['xy']])
  }

  # if z values are available
  if(isTRUE(hasZ)) {
    # z center of rotation adjustment
    spatlocs$sdimz = spatlocs$sdimz - rcenter[['z']]

    zvals = spatlocs$sdimz

    # Perform rotations
    if(rotateradians[['zy']] != 0) {
      spatlocs$sdimz = zvals*cos(rotateradians[['zy']]) + yvals*sin(rotateradians[['zy']])
      spatlocs$sdimy = -zvals*sin(rotateradians[['zy']]) + yvals*cos(rotateradians[['zy']])
    }

    if(rotateradians[['xz']] != 0) {
      spatlocs$sdimx = xvals*cos(rotateradians[['xz']]) + zvals*sin(rotateradians[['xz']])
      spatlocs$sdimz = -xvals*sin(rotateradians[['xz']]) + zvals*cos(rotateradians[['xz']])
    }

    # Revert z center of rotation adjustment
    spatlocs$sdimz = spatlocs$sdimz + rcenter[['z']]
  }

  # Revert xy center of rotation adjustment
  spatlocs$sdimx = spatlocs$sdimx + rcenter[['x']]
  spatlocs$sdimy = spatlocs$sdimy + rcenter[['y']]

  return(spatlocs)
}








#' @name flip_gpoly
#' @title Flip a giottoPolygon object
#' @description Flip a giottoPolygon over a designated x or y value depending on
#' direction param input. Note that this behavior is different from terra's
#' implementation of flip for SpatVectors where flips happen over the extent
#' @param gpoly giottoPolygon
#' @param direction character. Direction to flip. Should be either partial match to 'vertical' or 'horizontal'
#' @param x0 x value to flip horizontally over (ignored for vertical). Pass NULL
#' to flip over the extent
#' @param y0 y value to flip vertically over (ignored for horizontal). Pass NULL
#' to flip over the extent
#' @keywords internal
#' @noRd
flip_gpoly = function(gpoly,
                      direction = 'vertical',
                      x0 = 0,
                      y0 = 0) {
  checkmate::assertClass(gpoly, 'giottoPolygon')
  checkmate::assert_character(direction)
  if(!is.null(x0)) {
    checkmate::assert_numeric(x0)
  }
  if(!is.null(y0)) {
    checkmate::assert_numeric(y0)
  }

  # 1. perform flip
  e = terra::ext(gpoly@spatVector)
  gpoly = do_gpoly(x = gpoly,
                   what = terra::flip,
                   args = list(direction = direction))

  # 2. perform shift to match line of symmetry
  if(grepl(direction, 'vertical') & !is.null(y0)) {
    y_min = as.numeric(e$ymin)
    dy = y0 - y_min
    gpoly = do_gpoly(x = gpoly,
                     what = terra::shift,
                     args = list(dy = 2 * dy))
  }
  if(grepl(direction, 'horizontal') & !is.null(x0)) {
    x_min = as.numeric(e$xmin)
    dx = x0 - x_min
    gpoly = do_gpoly(x = gpoly,
                     what = terra::shift,
                     args = list(dx = 2 * dx))
  }

  # 3. return
  return(gpoly)
}




#' @name flip_gpoints
#' @title Flip a giottoPoints object
#' @description Flip a giottoPoints over a designated x or y value depending on
#' direction param input. Note that this behavior is different from terra's
#' implementation of flip for SpatVectors where flips happen over the extent
#' @param gpoly giottoPoints
#' @param direction character. Direction to flip. Should be either partial match to 'vertical' or 'horizontal'
#' @param x0 x value to flip horizontally over (ignored for vertical). Pass NULL
#' to flip over the extent
#' @param y0 y value to flip vertically over (ignored for horizontal). Pass NULL
#' to flip over the extent
#' @keywords internal
#' @noRd
flip_gpoints = function(gpoints,
                        direction = 'vertical',
                        x0 = 0,
                        y0 = 0) {
  checkmate::assertClass(gpoints, 'giottoPoints')
  checkmate::assert_character(direction)
  if(!is.null(x0)) {
    checkmate::assert_numeric(x0)
  }
  if(!is.null(y0)) {
    checkmate::assert_numeric(y0)
  }

  # !will need to update for networks information!

  # 1. perform flip
  e = terra::ext(gpoints@spatVector)
  gpoints@spatVector = terra::flip(gpoints@spatVector,
                                   direction = direction)

  # 2. perform shift to match line of symmetry
  if(grepl(direction, 'vertical') & !is.null(y0)) {
    y_min = as.numeric(e$ymin)
    dy = y0 - y_min
    gpoints@spatVector = terra::shift(x = gpoints@spatVector,
                                      dy = 2 * dy)
  }
  if(grepl(direction, 'horizontal') & !is.null(x0)) {
    x_min = as.numeric(e$xmin)
    dx = x0 - x_min
    gpoints@spatVector = terra::shift(x = gpoints@spatVector,
                                      dx = 2 * dx)
  }

  # 3. return
  return(gpoints)

}





#' @name flip_spatlocs
#' @param sl spatLocsObj
#' @param direction character. Direction to flip. Should be either partial match to 'vertical' or 'horizontal'
#' @param x0 x value to flip horizontally over (ignored for vertical). Pass NULL
#' to flip over the extent
#' @param y0 y value to flip vertically over (ignored for horizontal). Pass NULL
#' to flip over the extent
#' @keywords internal
#' @noRd
flip_spatlocs = function(sl,
                         direction = 'vertical',
                         x0 = 0,
                         y0 = 0,
                         copy_obj = TRUE) {
  sdimy = sdimx = NULL

  checkmate::assertClass(sl, 'spatLocsObj')
  checkmate::assert_character(direction)
  if(!is.null(x0)) {
    checkmate::assert_numeric(x0)
  }
  if(!is.null(y0)) {
    checkmate::assert_numeric(y0)
  }

  if(isTRUE(copy_obj)) sl = copy(sl)

  if(grepl(direction, 'vertical')) {
    y_min = sl[][, min(sdimy)]
    if(is.null(y0)) y0 = y_min
    sl[][, sdimy := -sdimy + (2 * y0)]
  }
  if(grepl(direction, 'horizontal')) {
    x_min = sl[][, min(sdimx)]
    if(is.null(x0)) x0 = x_min
    sl[][, sdimx := -sdimx + (2 * x0)]
  }

  return(sl)
}




#' @name flip_spatnet
#' @param sn spatialNetworkObj
#' @param direction character. Direction to flip. Should be either partial match to 'vertical' or 'horizontal'
#' @param x0 x value to flip horizontally over (ignored for vertical). Pass NULL
#' to flip over the extent
#' @param y0 y value to flip vertically over (ignored for horizontal). Pass NULL
#' to flip over the extent
#' @keywords internal
#' @noRd
flip_spatnet = function(sn,
                        direction = 'vertical',
                        x0 = 0,
                        y0 = 0,
                        copy_obj = TRUE) {
  sdimy_begin = sdimy_end = sdimx_begin = sdimx_end = NULL

  checkmate::assertClass(sn, 'spatialNetworkObj')
  checkmate::assert_character(direction)
  if(!is.null(x0)) {
    checkmate::assert_numeric(x0)
  }
  if(!is.null(y0)) {
    checkmate::assert_numeric(y0)
  }

  if(isTRUE(copy_obj)) sn = copy(sn)

  if(grepl(direction, 'vertical')) {
    y_min = sn[][, min(sdimy_begin, sdimy_end)]
    if(is.null(y0)) y0 = y_min
    sn[][, c('sdimy_begin', 'sdimy_end') := .(-sdimy_begin + (2 * y0), -sdimy_end + (2 * y0))]
    if(!is.null(sn@networkDT_before_filter))
      sn@networkDT_before_filter[, c('sdimy_begin', 'sdimy_end') := .(-sdimy_begin + (2 * y0), -sdimy_end + (2 * y0))]
  }
  if(grepl(direction, 'horizontal')) {
    x_min = sn[][, min(sdimx_begin, sdimx_end)]
    if(is.null(x0)) x0 = x_min
    sn[][, c('sdimx_begin', 'sdimx_end') := .(-sdimx_begin + (2 * x0), -sdimx_end + (2 * x0))]
    if(!is.null(sn@networkDT_before_filter))
      sn@networkDT_before_filter[, c('sdimx_begin', 'sdimx_end') := .(-sdimx_begin + (2 * x0), -sdimx_end + (2 * x0))]
  }

  return(sn)
}




#' @name flip_extent
#' @title Flip a SpatExtent
#' @param e extent
#' @param direction character. Direction to flip. Should be either partial match to 'vertical' or 'horizontal'
#' @param x0 x value to flip horizontally over (ignored for vertical). Pass NULL
#' to flip over the extent
#' @param y0 y value to flip vertically over (ignored for horizontal). Pass NULL
#' to flip over the extent
#' @keywords internal
#' @noRd
flip_extent = function(e,
                       direction = 'vertical',
                       x0 = 0,
                       y0 = 0) {

  checkmate::assertClass(e, 'SpatExtent')
  checkmate::assert_character(direction)
  if(!is.null(x0)) {
    checkmate::assert_numeric(x0)
  }
  if(!is.null(y0)) {
    checkmate::assert_numeric(y0)
  }

  y_vals = as.numeric(c(e$ymin, e$ymax))
  x_vals = as.numeric(c(e$xmin, e$xmax))

  if(grepl(direction, 'vertical')) {
    if(is.null(y0)) y0 = y_vals[1] # set bound min as line of sym (terra default)
    y_vals = -y_vals + (2 * y0)
  }
  if(grepl(direction, 'horizontal')) {
    if(is.null(x0)) x0 = x_vals[1] # set bound min as line of sym (terra default)
    x_vals = -x_vals + (2 * x0)
  }

  terra::ext(c(sort(x_vals), sort(y_vals)))
}


