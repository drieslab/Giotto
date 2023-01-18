# Methods and Generics ####

# NOTE: initialize generics are in classes.R #


#' @title Cell and feature names
#' @name spatIDs-generic
#' @description Get the cells and feature names of a giotto object or subobject
#' @aliases spatIDs featIDs
#' @param x an object
#' @param spat_unit (optional) specify which spatial unit
#' @param feat_type (optional) specify which feature type
#' @param ... additional parameters to pass
#' @include classes.R
#' @usage
#' spatIDs(x, spat_unit, ...)
#'
#' featIDs(x, feat_type, ...)
#'
#' ## Default S4 method for signatures:
#' ## 'giotto', 'exprObj', 'spatLocsObj', 'cellMetaObj', 'spatialNetworkObj' 'dimObj'
setGeneric('spatIDs', function(x, spat_unit, ...) standardGeneric('spatIDs'))
setGeneric('featIDs', function(x, feat_type, ...) standardGeneric('featIDs'))


setMethod('spatIDs', signature(x = 'giotto', spat_unit = 'missing'),
          function(x, ...) {
            get_cell_id(gobject = x, ...)
          })
setMethod('spatIDs', signature(x = 'giotto', spat_unit = 'character'),
          function(x, spat_unit, ...) {
            get_cell_id(gobject = x, spat_unit, ...)
          })
setMethod('spatIDs', signature(x = c('exprObj'), spat_unit = 'missing'),
          function(x, ...) {
            colnames(x[])
          })
setMethod('spatIDs', signature(x = c('spatLocsObj'), spat_unit = 'missing'),
          function(x, ...) {
            x[]$cell_ID
          })
setMethod('spatIDs', signature(x = c('cellMetaObj'), spat_unit = 'missing'),
          function(x, ...) {
            x[]$cell_ID
          })
setMethod('spatIDs', signature(x = c('spatialNetworkObj'), spat_unit = 'missing'),
          function(x, ...) {
            unique(c(x[]$from, x[]$to))
          })
setMethod('spatIDs', signature(x = 'dimObj', spat_unit = 'missing'),
          function(x, ...) {
            rownames(x@coordinates)
          })
setMethod('spatIDs', signature(x = 'giottoPolygon', spat_unit = 'missing'),
          function(x, ...) {
            unique(x@spatVector$poly_ID)
          })


setMethod('featIDs', signature(x = 'giotto', feat_type = 'missing'),
          function(x, ...) {
            get_feat_id(gobject = x, ...)
          })
setMethod('featIDs', signature(x = 'giotto', feat_type = 'character'),
          function(x, feat_type, ...) {
            get_feat_id(gobject = x, feat_type, ...)
          })
setMethod('featIDs', signature(x = 'exprObj', feat_type = 'missing'),
          function(x, ...) {
            rownames(x[])
          })
setMethod('featIDs', signature(x = 'giottoPoints', feat_type = 'missing'),
          function(x, ...) {
            unique(x@spatVector$feat_ID)
          })
# no features generic for dimObj







# colnames and rownames generics ####
setMethod('colnames', signature(x = 'exprObj'),
          function(x) {
            colnames(x[])
          })
setMethod('rownames', signature(x = 'exprObj'),
          function(x) {
            rownames(x[])
          })




# spin() S4 generic ####

#' @title Spin an object
#' @name spin-generic
#' @description Spin (rotate) an object spatially (limited to xy rotations)
#' @param x object
#' @param angle numeric. Angle of rotation in degrees
#' @param x0 numeric. x-coordinate of the center of rotation. Defaults to center x val if not given.
#' @param y0 numeric. y-coordinate of the center of rotation. Defaults to center y val if not given.
NULL

#' @describeIn spin-generic Spin a giottoPolygon object
#' @importMethodsFrom terra spin
#' @export
setMethod('spin', signature(x = 'giottoPolygon'),
          function(x, angle, x0 = NULL, y0 = NULL) {
            if(is.null(x0)) x0 = terra::mean(terra::ext(x@spatVector))[1]
            if(is.null(y0)) y0 = terra::mean(terra::ext(x@spatVector))[2]
            x@spatVector = terra::spin(x = x@spatVector,
                                       angle = angle,
                                       x0 = x0,
                                       y0 = y0)
            if(!is.null(x@spatVectorCentroids)) {
              x@spatVectorCentroids = terra::spin(x = x@spatVectorCentroids,
                                                  angle = angle,
                                                  x0 = x0,
                                                  y0 = y0)
            }
            if(!is.null(x@overlaps)) {
              lapply(x@overlaps, function(overlap) {
                if(inherits(overlap, 'SpatVector')) {
                  terra::spin(x = overlap,
                              angle = angle,
                              x0 = x0,
                              y0 = y0)
                } else {
                  overlap
                }
              })
            }
            return(x)
          })

#' @describeIn spin-generic Spin a giottoPoints object
#' @importMethodsFrom terra spin
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
#' @importMethodsFrom terra spin
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




# t() S4 generic ####

setMethod('t', signature(x = 'spatLocsObj'), function(x) {
  sdimy = sdimx = NULL
  x = data.table::copy(x)
  x@coordinates[, c('sdimx', 'sdimy') := .(sdimy, sdimx)]
  return(x)
})

# S3 definition
t.spatLocsObj = function(mymatrix) {
  t(mymatrix)
}






#' @title Dimensions of giotto objects
#' @name dims-generic
#' @description Find the dimensions of an object
#' @include classes.R
#' @param x object to check dimensions of
NULL

# nrow() S4 generic ####


# setMethod('nrow', signature = 'giotto', function(x) {
#   avail_exp = list_expression(x)
#   get_expression_values(x,
#                         spat_unit = avail_exp$spat_unit[1],
#                         feat_type = avail_exp$feat_type[1],
#                         values = avail_exp$name[1]) %>%
#     nrow()
# })

#' @describeIn dims-generic Find rows of giottoPoints object
#' @export
setMethod('nrow', signature('giottoPoints'), function(x) terra::nrow(x@spatVector))

#' @describeIn dims-generic Find rows of giottoPolygon object
#' @export
setMethod('nrow', signature('giottoPolygon'), function(x) terra::nrow(x@spatVector))

#' @describeIn dims-generic Find rows of giotto S4s with data.table based \code{coordinates} slots
#' @export
setMethod('nrow', signature('coordDataDT'), function(x) x@coordinates[,.N])

# TODO
# setMethod('dims', signature('coordDataMT'), function(x) nrow(x@coordinates))

#' @describeIn dims-generic Find rows of giotto S4s with data.table based \code{coordinates} slots
#' @export
setMethod('nrow', signature('exprData'), function(x) nrow(x@exprMat))

#' @describeIn dims-generic Find rows of giotto S4s with data.table based \code{coordinates} slots
#' @export
setMethod('nrow', signature('metaData'), function(x) nrow(x@metaDT))

# ncol() generic ####

#' @describeIn dims-generic Find cols of giotto S4s with Matrix based \code{exprMat} slots
#' @export
setMethod('ncol', signature('exprData'), function(x) ncol(x@exprMat))

#' @describeIn dims-generic Find cols of giotto S4s with data.table based \code{metaDT} slots
#' @export
setMethod('ncol', signature('metaData'), function(x) ncol(x@metaDT))

# dim() generic ####

#' @describeIn dims-generic Find dimensions of giotto S4s with Matrix based \code{exprMat} slots
#' @export
setMethod('dim', signature('exprData'), function(x) dim(x@exprMat))

#' @describeIn dims-generic Find dimensions of giotto S4s with data.table based \code{metaDT} slots
#' @export
setMethod('dim', signature('metaData'), function(x) dim(x@metaDT))



# wrap() generic ####

#' @title Wrap giotto terra pointer information
#' @name wrap-generic
#' @description Extension of wrap methods from terra for Giotto's terra-based S4
#' objects. Allows pointer information to be packaged into memory so that it can
#' be passed over a connection (e.g. nodes on a computer cluster)
#' @param x giottoPolygon or giottoPoints
NULL


#' @describeIn wrap-generic Wrap giottoPolygon
#' @importMethodsFrom terra wrap
#' @export
setMethod('wrap', signature(x = 'giottoPolygon'),
          function(x) {
            pgp = new('packedGiottoPolygon')
            pgp@name = x@name
            pgp@packed_spatVector = terra::wrap(x@spatVector)
            if(!is.null(x@spatVectorCentroids)) {
              pgp@packed_spatVectorCentroids = terra::wrap(x@spatVectorCentroids)
            }
            if(!is.null(x@overlaps)) {
              pgp@packed_overlaps = lapply(x@overlaps, function(sv) {
                if(inherits(sv, 'SpatVector')) {
                  terra::wrap(sv)
                } else {
                  sv
                }
              })
            }
            return(pgp)
          }
)


#' @describeIn wrap-generic Wrap giotto
#' @importMethodsFrom terra wrap
#' @export
setMethod('wrap', signature(x = 'giotto'),
          function(x) {
            pg = new('packedGiotto')
            g_slots = methods::slotNames('giotto')
            g_slots = g_slots[!g_slots %in% c('spatial_info', 'feat_info')]
            for(g_slot in g_slots) {
              slot(pg, g_slot) = slot(x, g_slot)
            }
            pg@packed_spatial_info = lapply(x@spatial_info, wrap)
            pg@packed_feat_info = lapply(x@feat_info, wrap)
            return(pg)
          })


#' @describeIn wrap-generic Wrap giottoPoints
#' @importMethodsFrom terra wrap
#' @export
setMethod('wrap', signature(x = 'giottoPoints'),
          function(x) {
            pgp = new('packedGiottoPoints')
            pgp@feat_type = x@feat_type
            pgp@packed_spatVector = terra::wrap(x@spatVector)
            pgp@networks = x@networks
            return(pgp)
          }
)



#' @describeIn wrap-generic Unwrap giottoPolygon
#' @importMethodsFrom terra unwrap
#' @export
setMethod('unwrap', signature(x = 'packedGiottoPolygon'),
          function(x) {
            gp = new('giottoPolygon')
            gp@name = x@name
            gp@spatVector = terra::unwrap(x@packed_spatVector)
            if(!is.null(x@packed_spatVectorCentroids)) {
              gp@spatVectorCentroids = terra::unwrap(x@packed_spatVectorCentroids)
            }
            if(length(x@packed_overlaps) > 0) {
              gp@overlaps = lapply(x@packed_overlaps, function(sv) {
                if(inherits(sv, 'PackedSpatVector')) {
                  terra::unwrap(sv)
                } else {
                  sv
                }
              })
            }
            return(gp)
          }
)



#' @describeIn wrap-generic Unwrap giottoPolygon
#' @importMethodsFrom terra unwrap
#' @export
setMethod('unwrap', signature(x = 'packedGiottoPoints'),
          function(x) {
            gp = new('giottoPoints')
            gp@feat_type = x@feat_type
            gp@spatVector = terra::unwrap(x@packed_spatVector)
            gp@networks = x@networks
            return(gp)
            }
)


#' @describeIn wrap-generic Unwrap giotto
#' @importMethodsFrom terra unwrap
#' @export
setMethod('unwrap', signature(x = 'packedGiotto'),
          function(x) {
            gobj = new('giotto')
            g_slots = methods::slotNames('giotto')
            g_slots = g_slots[!g_slots %in% c('spatial_info', 'feat_info')]
            for(g_slot in g_slots) {
              slot(gobj, g_slot) = slot(x, g_slot)
            }
            gobj@spatial_info = lapply(x@packed_spatial_info, unwrap)
            gobj@feat_info = lapply(x@packed_feat_info, unwrap)
            return(gobj)
          })





# rbind() generic ####

#' @title Combine objects by rows (Giotto-related)
#' @name rbind-generic
#' @description row bind two objects
#' @include classes.R
#' @param x item 1 to rbind
#' @param y item 2 to rbind
#' @param add_list_ID whether to generate a list_ID column when giottoPolygons
#' to append have different names
#' @param \dots additional params to pass to methods
NULL


#' @describeIn rbind-generic Append giottoPolygon objects
#' @importFrom methods rbind2
#' @export
setMethod('rbind2', signature(x = 'giottoPolygon', y = 'giottoPolygon'),
          function(x, y, add_list_ID = TRUE, ...) {

  # determine homo or hetero
  poly_names = c(slot(x, 'name'), slot(y, 'name'))
  homo = identical(poly_names[[1L]], poly_names[[2L]])
  if(!isTRUE(homo)) {
    new_name = paste0(sort(poly_names), collapse = '-')
    return(rbind2_giotto_polygon_hetero(x = x, y = y, new_name = new_name, add_list_ID = add_list_ID))
  } else {
    return(rbind2_giotto_polygon_homo(x = x, y = y))
  }
})


setGeneric('rbind', signature = '...')

setMethod("rbind", "giottoPolygon", function(..., deparse.level = 1) {
  if(nargs() <= 2L) {
    rbind2(...)
  } else {
    xs = list(...)
    rbind2(xs[[1L]], do.call(Recall, xs[-1L]))
  }
})




# plot() S4 Generic ####

#' @title Preview a Giotto spatial object
#' @name plot-generic
#' @description S4 generic for previewing Giotto's image and subcellular objects.
#' @include classes.R
#' @param x giotto image, giottoPolygon, or giottoPoints object
#' @param y Not used.
#' @param \dots additional parameters to pass
#' @aliases plot
#' @family plot
#' @exportMethod plot
setOldClass('plot')

#' @describeIn plot-generic Plot \emph{magick}-based giottoImage object. ... param passes to \code{\link{plot_giottoImage_MG}}
#' @export
setMethod('plot', signature(x = 'giottoImage', y = 'missing'), function(x,y,...) plot_giottoImage_MG(giottoImage = x,...))

#' @describeIn plot-generic Plot \emph{terra}-based giottoLargeImage object. ... param passes to \code{\link{plot_giottoLargeImage}}
#' @export
setMethod('plot', signature(x = 'giottoLargeImage', y = 'missing'), function(x,y,...) plot_giottoLargeImage(giottoLargeImage = x,...))

#' @describeIn plot-generic Plot \emph{terra}-based giottoPolygon object. ... param passes to \code{\link[terra]{plot}}
#' @param point_size size of points when plotting giottoPolygon object centroids
#' @param type what to plot: either 'poly' (default) or polygon 'centroid'
#' @export
setMethod('plot', signature(x = 'giottoPolygon', y = 'missing'), function(x,
                                                                          y,
                                                                          point_size = 0.1,
                                                                          type = c('poly', 'centroid'),
                                                                          ...) {
  type = match.arg(type, choices = c('poly', 'centroid'))
  if(type == 'poly') {
    terra::plot(x = x@spatVector, ...)
  }
  if(type == 'centroid') {
    if(!is.null(x@spatVectorCentroids)) {
      terra::plot(x = x@spatVectorCentroids, cex = point_size, ...)
    } else {
      cat('no centroids calculated\n')
    }
  }

})

#' @describeIn plot-generic \emph{terra}-based giottoPoint object. ... param passes to \code{\link[terra]{plot}}
#' @param point_size size of points when plotting giottoPoints
#' @param feats specific features to plot within giottoPoints object (defaults to NULL, meaning all available features)
#' @export
setMethod('plot', signature(x = 'giottoPoints', y = 'missing'), function(x,
                                                                         y,
                                                                         point_size = 0.1,
                                                                         feats = NULL,
                                                                         ...) {

  if(is.null(feats)) terra::plot(x = x@spatVector, cex = point_size, ...)
  else if(length(feats) == 1) {
    gp = x@spatVector
    x_feat_subset = gp[gp$feat_ID %in% feats]
    terra::plot(x = x_feat_subset, cex = point_size, ...)
  }
  else {
    gp = x@spatVector
    x_feat_subset = gp[gp$feat_ID %in% feats]
    terra::plot(x = x_feat_subset, cex = point_size, 'feat_ID', ...)
  }
})


#' @describeIn plot-generic Plot a spatLocsObj
#' @export
setMethod('plot', signature(x = 'spatLocsObj', y = 'missing'), function(x, ...) {
  l = list(...)
  if(is.null(l$asp)) l$asp = 1
  if(is.null(l$xlab)) l$xlab = ''
  if(is.null(l$ylab)) l$ylab = ''
  if(is.null(l$cex)) l$cex = 0.5
  if(nrow(x) > 10000) {
    if(is.null(l$pch)) l$pch = '.'
  }

  do.call('plot', append(l, list(x = x[]$sdimx, y = x[]$sdimy)))
})



# copy() S4 generic ####
#' @title Copy an entire object
#' @name copy-generic
#' @description S4 generic for Giotto's S4 subobjects to return with full copies of
#' certain subobjects that usually return referenced information.
#' @param x a Giotto S4 class subobject
#' @include classes.R
#' @seealso \code{\link[data.table]{copy}} \code{\link[terra]{deepcopy}}
#' @aliases copy
#' @exportMethod copy
setGeneric('copy', function(x) {
  standardGeneric('copy')
})

#' @describeIn copy-generic Copy \emph{data.table}-based spatial locations object.
#' @export
setMethod('copy', signature(x = 'coordDataDT'), function(x) {
  x@coordinates = data.table::copy(x@coordinates)
  x
})








# prov() S4 generic ####
#' @title Provenance information
#' @name prov-generic
#' @description access and set provenance slot of S4 subobject
#' @param x a Giotto S4 class subobject
#' @param value value to set as provenance
#' @aliases prov prov<-
#' @include classes.R
#' @export
setGeneric('prov', function(x) standardGeneric('prov'))
setGeneric('prov<-', function(x, value) standardGeneric('prov<-'))


#' @describeIn prov-generic Get provenance information
#' @export
setMethod('prov', signature = 'provData', function(x) x@provenance)


#' @describeIn prov-generic Set provenance information
#' @export
setMethod('prov<-', signature = 'provData', function(x, value) {
  x@provenance = value
  x
})




# spatUnit() S4 generic ####
#' @title Spatial unit information
#' @name spatUnit-generic
#' @description access and set spat_unit slot of S4 subobject
#' @param x a Giotto S4 class subobject with spatial unit
#' @param value value to set as spatial unit
#' @aliases spatUnit spatUnit<-
#' @include classes.R
#' @export
setGeneric('spatUnit', function(x) standardGeneric('spatUnit'))
setGeneric('spatUnit<-', function(x, value) standardGeneric('spatUnit<-'))


#' @describeIn spatUnit-generic Get spatial unit information
#' @export
setMethod('spatUnit', signature = 'spatData', function(x) x@spat_unit)


#' @describeIn spatUnit-generic Set spatial unit information
#' @export
setMethod('spatUnit<-', signature = 'spatData', function(x, value) {
  x@spat_unit = value
  x
})




# featType() S4 generic ####
#' @title Feature type information
#' @name featType-generic
#' @description access and set feat_type slot of S4 subobject
#' @param x a Giotto S4 class subobject with feature type
#' @param value value to set as feature type
#' @aliases featType featType<-
#' @include classes.R
#' @export
setGeneric('featType', function(x) standardGeneric('featType'))
setGeneric('featType<-', function(x, value) standardGeneric('featType<-'))


#' @describeIn featType-generic Get feature type information
#' @export
setMethod('featType', signature = 'featData', function(x) x@feat_type)


#' @describeIn featType-generic Set feature type information
#' @export
setMethod('featType<-', signature = 'featData', function(x, value) {
  x@feat_type = value
  x
})




# objName() generic ####
#' @title Giotto object name information
#' @name objName-generic
#' @description access and set name slot fo S4 subobject
#' @param x a Giotto S4 class subobject with name data
#' @param value value to set as object name
#' @aliases objName objName<-
#' @include classes.R
#' @export
setGeneric('objName', function(x) standardGeneric('objName'))
setGeneric('objName<-', function(x, value) standardGeneric('objName<-'))

#' @describeIn objName-generic Get name information
#' @export
setMethod('objName', signature = 'nameData', function(x) x@name)

#' @describeIn objName-generic Set name information
#' @export
setMethod('objName<-', signature = 'nameData', function(x, value) {
  x@name = value
  x
})




#' @title Extract or replace parts of an object
#' @name extract-methods
#' @docType methods
#' @aliases `[` `[<-` `$` `$<-`
#' @description Operators Giotto S4 internal data.tables to extract
#' or replace parts.
#' @param x Giotto S4 object to extract columns from
#' @param i,j indices specifying elements to extract or replace. Indices are
#' numeric or character vectors or empty
#' @param name A literal character string (possibly backtick quoted).
#' @param value value(s) to set
#' This is normally matched to the colnames of the data.table object within the S4.
#' @section \code{`$`} methods:
#' @section \code{`$<-`} methods:
#' @section \code{`[`} methods:
#' @section \code{`[<-`} methods:
#'
NULL

# $ S4 access generic ####

## * coordDataDT ####

#' @rdname extract-methods
#' @section \code{`$`} methods:
#'   Select by colname from giotto S4 data.table coordinates slot.
#' @export
setMethod('$', signature(x = 'coordDataDT'),
          function(x, name) x@coordinates[[name]])


#' @rdname extract-methods
#' @section \code{`$<-`} methods:
#'   Set values by colname into giotto S4 data.table coordinates slot.
#'   Works via data.table methods
#' @export
setMethod('$<-', signature(x = 'coordDataDT'),
          function(x, name, value) {
            if(x@coordinates[,.N] == 0) x@coordinates = data.table::data.table()
            x@coordinates[, (name) := value]
            x
          })

## * metaData ####

#' @rdname extract-methods
#' @section \code{`$`} methods:
#'   Select by colname from giotto S4 data.table metaDT slot.
#' @export
setMethod('$', signature(x = 'metaData'),
          function(x, name) x@metaDT[[name]])


#' @rdname extract-methods
#' @section \code{`$<-`} methods:
#'   Set values by colname into giotto S4 data.table metaDT slot.
#'   Works via data.table methods
#' @export
setMethod('$<-', signature(x = 'metaData'),
          function(x, name, value) {
            if(x@metaDT[,.N] == 0) x@metaDT = data.table::data.table()
            x@metaDT[, (name) := value]
            x
          })


# [ S4 access generic ####

## * coordDataDT ####



#' @rdname extract-methods
#' @section \code{`[`} methods:
#'   Select rows (i) and cols (j) from giotto S4 coordinates slot
#' @export
setMethod('[', signature(x = 'coordDataDT', i = 'missing', j = 'ANY', drop = 'missing'),
          function(x, i, j) {
            x@coordinates = x@coordinates[, j = j, with = FALSE]
            x
          })

#' @rdname extract-methods
#' @export
setMethod('[', signature(x = 'coordDataDT', i = 'ANY', j = 'missing', drop = 'missing'),
          function(x, i, j) {
            x@coordinates = x@coordinates[i = i,]
            x
          })

#' @rdname extract-methods
#' @export
setMethod('[', signature(x = 'coordDataDT', i = 'ANY', j = 'ANY', drop = 'missing'),
          function(x, i, j) {
            x@coordinates = x@coordinates[i = i, j = j, with = FALSE]
            x
          })

#' @name [
#' @rdname extract-methods
#' @aliases [,coordDataDT,missing,missing,missing-method
#' @section \code{`[`} methods:
#'   Return \code{coordinates} slot data.table from giotto S4
#' @export
setMethod('[', signature(x = 'coordDataDT', i = 'missing', j = 'missing', drop = 'missing'),
          function(x, i, j) {
            x@coordinates
          })

#' @name [
#' @rdname extract-methods
#' @aliases [<-,coordDataDT,missing,missing,ANY-method [<-,coordDataDT,missing,missing-method
#' @docType methods
#' @section \code{`[<-`} methods:
#'   Assign to \code{coordinates} slot in giotto S4
setReplaceMethod('[', signature(x = 'coordDataDT', i = 'missing', j = 'missing', value = 'ANY'),
          function(x, i, j, value) {
            x@coordinates = value
            x
          })

# setMethod("[[")

## * metaData ####

#' @rdname extract-methods
#' @section \code{`[`} methods:
#'   Select rows (i) and cols (j) from giotto S4 metaDT slot
#' @export
setMethod('[', signature(x = 'metaData', i = 'missing', j = 'ANY', drop = 'missing'),
          function(x, i, j) {
            x@metaDT = x@metaDT[, j = j, with = FALSE]
            x
          })

#' @rdname extract-methods
#' @export
setMethod('[', signature(x = 'metaData', i = 'ANY', j = 'missing', drop = 'missing'),
          function(x, i, j) {
            x@metaDT = x@metaDT[i = i,]
            x
          })

#' @rdname extract-methods
#' @export
setMethod('[', signature(x = 'metaData', i = 'ANY', j = 'ANY', drop = 'missing'),
          function(x, i, j) {
            x@metaDT = x@metaDT[i = i, j = j, with = FALSE]
            x
          })

#' @rdname extract-methods
#' @section \code{`[`} methods:
#'   Return \code{metaDT} slot data.table from giotto S4
#' @export
setMethod('[', signature(x = 'metaData', i = 'missing', j = 'missing', drop = 'missing'),
          function(x, i, j) {
            x@metaDT
          })

#' @rdname extract-methods
#' @aliases [<-,metaData,missing,missing,ANY-method [<-,metaData,missing,missing-method
#' @docType methods
#' @section \code{`[<-`} methods:
#'   Assign to \code{metaDT} slot in giotto S4
#' @export
setMethod('[<-', signature(x = 'metaData', i = 'missing', j = 'missing', value = 'ANY'),
          function(x, i, j, value) {
            x@metaDT = value
            x
          })


## * dimObj (temp) ####

#' @rdname extract-methods
#' @export
setMethod('[', signature(x = 'dimObj', i = 'ANY', j = 'ANY', drop = 'missing'),
          function(x, i, j) {
            x@coordinates = x@coordinates[i = i, j = j]
            x
          })

#' @rdname extract-methods
#' @section \code{`[`} methods:
#'    Return \code{coordinates} slot matrix from giotto S4 dimObj
#' @export
setMethod('[', signature(x = 'dimObj', i = 'missing', j = 'missing', drop = 'missing'),
          function(x, i, j) {
            x@coordinates
          })

#' @rdname extract-methods
#' @aliases [<-,dimObj,missing,missing,ANY-method [<-,dimObj,missing,missing-method
#' @docType methods
#' @section \code{`[<-`} methods:
#'   Assign to \code{coordinates} slot in giotto S4 dimObj
#' @export
setMethod('[<-', signature(x = 'dimObj', i = 'missing', j = 'missing', value = 'ANY'),
          function(x, i, j, value) {
            x@coordinates = value
            x
          })

## * exprData ####

#' @rdname extract-methods
#' @section \code{`[`} methods:
#'   Select rows (i) and cols (j) from giotto S4 exprMat slot
#' @export
setMethod('[', signature(x = 'exprData', i = 'missing', j = 'ANY', drop = 'missing'),
          function(x, i, j) {
            x@exprMat = x@exprMat[, j = j]
            x
          })

#' @rdname extract-methods
#' @export
setMethod('[', signature(x = 'exprData', i = 'ANY', j = 'missing', drop = 'missing'),
          function(x, i, j) {
            x@exprMat = x@exprMat[i = i,]
            x
          })

#' @rdname extract-methods
#' @export
setMethod('[', signature(x = 'exprData', i = 'ANY', j = 'ANY', drop = 'missing'),
          function(x, i, j) {
            x@exprMat = x@exprMat[i = i, j = j]
            x
          })

#' @rdname extract-methods
#' @section \code{`[`} methods:
#'   Return \code{exprMat} slot Matrix object from giotto S4
#' @export
setMethod('[', signature(x = 'exprData', i = 'missing', j = 'missing', drop = 'missing'),
          function(x, i, j) {
            x@exprMat
          })

#' @rdname extract-methods
#' @aliases [<-,exprData,missing,missing,ANY-method [<-,exprData,missing,missing-method
#' @docType methods
#' @section \code{`[<-`} methods:
#'   Assign to \code{exprMat} slot in giotto S4
#' @export
setMethod('[<-', signature(x = 'exprData', i = 'missing', j = 'missing', value = 'ANY'),
          function(x, i, j, value) {
            x@exprMat = value
            x
          })

# * spatNetData ####
#' @rdname extract-methods
#' @section \code{`[`} methods:
#'   Return \code{spatNetData} slot network data.table object from giotto S4
#' @export
setMethod('[', signature(x = 'spatNetData', i = 'missing', j = 'missing', drop = 'missing'),
          function(x, i, j) {
            x@networkDT
          })

#' @rdname extract-methods
#' @aliases [<-,spatNetData,missing,missing,ANY-method [<-,spatNetData,missing,missing-method
#' @docType methods
#' @section \code{`[<-`} methods:
#'   Assign to \code{networkDT} slot in giotto S4
#' @export
setMethod('[<-', signature(x = 'spatNetData', i = 'missing', j = 'missing', value = 'ANY'),
          function(x, i, j, value) {
            x@networkDT = value
            x
          })


# * nnData ####
#' @rdname extract-methods
#' @section \code{`[`} methods:
#'   Return \code{nnData} slot igraph object from giotto S4
#' @export
setMethod('[', signature(x = 'nnData', i = 'missing', j = 'missing', drop = 'missing'),
          function(x, i, j) {
            x@igraph
          })

#' @rdname extract-methods
#' @aliases [<-,nnData,missing,missing,ANY-method [<-,nnData,missing,missing-method
#' @docType methods
#' @section \code{`[<-`} methods:
#'   Assign to \code{igraph} slot in giotto S4
#' @export
setMethod('[<-', signature(x = 'nnData', i = 'missing', j = 'missing', value = 'ANY'),
          function(x, i, j, value) {
            x@igraph = value
            x
          })


# * enrData ####
#' @rdname extract-methods
#' @section \code{`[`} methods:
#'   Return \code{enrData} slot enrichment data.table object from giotto S4
#' @export
setMethod('[', signature(x = 'enrData', i = 'missing', j = 'missing', drop = 'missing'),
          function(x, i, j) {
            x@enrichDT
          })

#' @rdname extract-methods
#' @aliases [<-,enrData,missing,missing,ANY-method [<-,enrData,missing,missing-method
#' @docType methods
#' @section \code{`[<-`} methods:
#'   Assign to \code{enrichDT} slot in giotto S4
#' @export
setMethod('[<-', signature(x = 'enrData', i = 'missing', j = 'missing', value = 'ANY'),
          function(x, i, j, value) {
            x@enrichDT = value
            x
          })

# * spatGridData ####
#' @rdname extract-methods
#' @section \code{`[`} methods:
#'   Return \code{spatGridData} slot data.table object from giotto S4
#' @export
setMethod('[', signature(x = 'spatGridData', i = 'missing', j = 'missing', drop = 'missing'),
          function(x, i, j) {
            x@gridDT
          })

#' @rdname extract-methods
#' @aliases [<-,spatGridData,missing,missing,ANY-method [<-,spatGridData,missing,missing-method
#' @docType methods
#' @section \code{`[<-`} methods:
#'   Assign to \code{gridDT} slot in giotto S4
#' @export
setMethod('[<-', signature(x = 'spatGridData', i = 'missing', j = 'missing', value = 'ANY'),
          function(x, i, j, value) {
            x@gridDT = value
            x
          })





