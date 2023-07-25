
#' @include giotto_structures.R
NULL

# Methods and Generics ####

# NOTE: initialize generics are in classes.R #



# Giotto Object Settings ####

## spatIDs and featIDs generic ####
#' @title Spatial and feature IDs
#' @name spatIDs-generic
#' @description Get the cell IDs (termed spatial IDs to better reflect when not at
#' the single-cell level) and feature IDs of a giotto object or subobject
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

#' @rdname spatIDs-generic
#' @export
setMethod('spatIDs', signature(x = 'giotto', spat_unit = 'missing'),
          function(x, ...) {
            as.character(get_cell_id(gobject = x, ...))
          })
#' @rdname spatIDs-generic
#' @export
setMethod('spatIDs', signature(x = 'giotto', spat_unit = 'character'),
          function(x, spat_unit, ...) {
            as.character(get_cell_id(gobject = x, spat_unit, ...))
          })
#' @rdname spatIDs-generic
#' @export
setMethod('spatIDs', signature(x = c('exprObj'), spat_unit = 'missing'),
          function(x, ...) {
            as.character(colnames(x[]))
          })
#' @rdname spatIDs-generic
#' @export
setMethod('spatIDs', signature(x = c('spatLocsObj'), spat_unit = 'missing'),
          function(x, ...) {
            as.character(x[]$cell_ID)
          })
#' @rdname spatIDs-generic
#' @export
setMethod('spatIDs', signature(x = c('cellMetaObj'), spat_unit = 'missing'),
          function(x, ...) {
            as.character(x[]$cell_ID)
          })
#' @rdname spatIDs-generic
#' @export
setMethod('spatIDs', signature(x = c('spatialNetworkObj'), spat_unit = 'missing'),
          function(x, ...) {
            as.character(unique(c(x[]$from, x[]$to)))
          })
#' @rdname spatIDs-generic
#' @export
setMethod('spatIDs', signature(x = 'dimObj', spat_unit = 'missing'),
          function(x, ...) {
            as.character(rownames(x@coordinates))
          })
#' @rdname spatIDs-generic
#' @param use_cache use cached IDs if available (gpoly and gpoints only)
#' @export
setMethod('spatIDs', signature(x = 'giottoPolygon', spat_unit = 'missing'),
          function(x, use_cache = TRUE, ...) {
            # getting as list first is more performant
            if(!all(is.na(x@unique_ID_cache)) & isTRUE(use_cache)) as.character(x@unique_ID_cache)
            else as.character(unique(terra::as.list(x@spatVector)$poly_ID))
          })
#' @rdname spatIDs-generic
#' @export
setMethod('spatIDs', signature(x = 'spatEnrObj', spat_unit = 'missing'),
          function(x, ...) {
            as.character(x@enrichDT$cell_ID)
          })
#' @rdname spatIDs-generic
#' @export
setMethod('spatIDs', signature(x = 'nnNetObj', spat_unit = 'missing'),
          function(x, ...) {
            as.character(unique(names(igraph::V(x@igraph))))
          })





# feat IDs

#' @rdname spatIDs-generic
#' @export
setMethod('featIDs', signature(x = 'giotto', feat_type = 'missing'),
          function(x, ...) {
            as.character(get_feat_id(gobject = x, ...))
          })
#' @rdname spatIDs-generic
#' @export
setMethod('featIDs', signature(x = 'giotto', feat_type = 'character'),
          function(x, feat_type, ...) {
            as.character(get_feat_id(gobject = x, feat_type, ...))
          })
#' @rdname spatIDs-generic
#' @export
setMethod('featIDs', signature(x = 'exprObj', feat_type = 'missing'),
          function(x, ...) {
            as.character(rownames(x[]))
          })
#' @rdname spatIDs-generic
#' @export
setMethod('featIDs', signature(x = 'giottoPoints', feat_type = 'missing'),
          function(x, use_cache = TRUE, ...) {
            # getting as list is more performant than directly using `$`
            if(!all(is.na(x@unique_ID_cache)) & isTRUE(use_cache)) as.character(x@unique_ID_cache)
            else as.character(unique(terra::as.list(x@spatVector)$feat_ID))
          })
#' @rdname spatIDs-generic
#' @export
setMethod('featIDs', signature(x = 'spatEnrObj', feat_type = 'missing'),
          function(x, ...) {
            as.character(colnames(x@enrichDT[, -'cell_ID']))
          })
# no features generic for dimObj








## default spat unit ####
#' @title Active spatial unit
#' @name activeSpatUnit-generic
#' @aliases activeSpatUnit activeSpatUnit<-
#' @description Retrieve or set the active spatial unit. This value will be the
#' default spatial unit that the giotto object uses.
#' @inheritParams data_access_params
setGeneric('activeSpatUnit', function(gobject, ...) standardGeneric('activeSpatUnit'))
setGeneric('activeSpatUnit<-', function(gobject, ..., value) standardGeneric('activeSpatUnit<-'))

#' @rdname activeSpatUnit-generic
#' @export
setMethod('activeSpatUnit', signature(gobject = 'giotto'), function(gobject) {
  su_try = try(instructions(gobject, 'active_spat_unit'), silent = TRUE)
  if(inherits(su_try, 'try-error')) su_try = NULL
  return(su_try)
})


#' @rdname activeSpatUnit-generic
#' @export
setMethod('activeSpatUnit<-', signature(gobject = 'giotto', value = 'character'),
          function(gobject, value) {
            instructions(gobject, 'active_spat_unit') = value
            return(gobject)
          })


## default feature type ####
#' @title Active feature type
#' @name activeFeatType-generic
#' @aliases activeFeatType activeFeatType<-
#' @description Retrieve or set the active feature type. This value will be the
#' default feature type that the giotto object uses.
#' @inheritParams data_access_params
setGeneric('activeFeatType', function(gobject, ...) standardGeneric('activeFeatType'))
setGeneric('activeFeatType<-', function(gobject, ..., value) standardGeneric('activeFeatType<-'))

#' @rdname activeFeatType-generic
#' @export
setMethod('activeFeatType', signature(gobject = 'giotto'), function(gobject) {
  ft_try = try(instructions(gobject, 'active_feat_type'), silent = TRUE)
  if(inherits(ft_try, 'try-error')) ft_try = NULL
  return(ft_try)
})


#' @rdname activeFeatType-generic
#' @export
setMethod('activeFeatType<-', signature(gobject = 'giotto', value = 'character'),
          function(gobject, value) {
            instructions(gobject, 'active_feat_type') = value
            return(gobject)
          })





## instructions ####
#' @title Access giotto instructions
#' @name instructions-generic
#' @aliases instructions instructions<-
#' @description Retrieve or set giotto instructions. Specific instructions can
#' be replaced using the \code{field} param. Additionally, when using
#' instructions<-, \code{initialize()} will be called on the giotto object if
#' initialize param is TRUE
#' @inheritParams data_access_params
#' @param param Specific param in instructions to access or modify
#' @param initialize (boolean, default = TRUE) whether to initialize the giotto object
#' @param value value to set
setGeneric('instructions', function(gobject, param, ...) standardGeneric('instructions'))
setGeneric('instructions<-', function(gobject, param, initialize, ..., value) standardGeneric('instructions<-'))



# Get instructions object
#' @rdname instructions-generic
#' @export
setMethod('instructions', signature(gobject = 'giotto', param = 'missing'),
          function(gobject) {
            return(showGiottoInstructions(gobject))
          })

# Set instructions object
#' @rdname instructions-generic
#' @export
setMethod('instructions<-',
          signature(gobject = 'giotto', param = 'missing', initialize = 'missing', value = 'ANY'),
          function(gobject, initialize, value) {
            gobject = replaceGiottoInstructions(gobject,
                                                instructions = value,
                                                init_gobject = TRUE)
            return(gobject)
          })
#' @rdname instructions-generic
#' @export
setMethod('instructions<-',
          signature(gobject = 'giotto', param = 'missing', initialize = 'logical', value = 'ANY'),
          function(gobject, initialize, value) {
            gobject = replaceGiottoInstructions(gobject,
                                                instructions = value,
                                                init_gobject = initialize)
            return(gobject)
          })

# Get specific field
#' @rdname instructions-generic
#' @export
setMethod('instructions', signature(gobject = 'giotto', param = 'character'),
          function(gobject, param) {
            instrs = showGiottoInstructions(gobject = gobject)
            return(readGiottoInstructions(giotto_instructions = instrs, param = param))
          })

# Set specific field
#' @rdname instructions-generic
#' @export
setMethod('instructions<-',
          signature(gobject = 'giotto', param = 'character', initialize = 'missing', value = 'ANY'),
          function(gobject, param, initialize, value) {
            gobject = changeGiottoInstructions(gobject = gobject,
                                               params = param,
                                               new_values = value,
                                               return_gobject = TRUE,
                                               init_gobject = TRUE)
            return(gobject)
          })
#' @rdname instructions-generic
#' @export
setMethod('instructions<-',
          signature(gobject = 'giotto', param = 'character', initialize = 'logical', value = 'ANY'),
          function(gobject, param, initialize, value) {
            gobject = changeGiottoInstructions(gobject = gobject,
                                               params = param,
                                               new_values = value,
                                               return_gobject = TRUE,
                                               init_gobject = initialize)
            return(gobject)
          })














# centroids() S4 generic ####
#' @title centroids-generic
#' @name centroids-generic
#' @description Access centroids information from polygon objects
#' @param x object
#' @aliases centroids
#' @details For giottoPolygon, if centroids already exist, pulls from
#' \code{spatVectorCentroids} slot. Otherwise, generates from
#' \code{spatVector} slot de novo
#' @importMethodsFrom terra centroids
NULL

#' @rdname centroids-generic
#' @export
setMethod('centroids', signature(x = 'giottoPolygon'),
          function(x) {
            if(!is.null(x@spatVectorCentroids)) {
              return(x@spatVectorCentroids)
            } else {
              return(terra::centroids(x@spatVector))
            }
          })




#' @name spatShift
#' @title Spatially shift an object
#' @param x object
#' @param dx numeric. The shift on the x axis
#' @param dy numeric. The shift on the y axis
#' @param dz numeric. The shift on the z axis
#' @param copy_obj Default = TRUE
#' @param ... additional params to pass to methods
#' @description Shift the spatial locations of an object
#' @export
setGeneric('spatShift', function(x, ...) standardGeneric('spatShift'))






#' @title overlaps-generic
#' @name overlaps-generic
#' @description Access list of overlaps information from object
#' @param x object
#' @aliases overlaps
setGeneric('overlaps', function(x, ...) standardGeneric('overlaps'))

#' @describeIn overlaps-generic Get overlaps information from giottoPolygon
#' @param name (optional) name of overlaps information to retrieve
#' @export
setMethod('overlaps', signature(x = 'giottoPolygon'),
          function(x, name = NULL) {
            if(is.null(name)) {
              # return entire list
              return(x@overlaps)
            } else {
              # return named entry
              return(x@overlaps[[name]])
            }
          })










# Object Characteristics ####

#' @title Dimensions of giotto objects
#' @name dims-generic
#' @description Find the dimensions of an object
#' @include classes.R
#' @importMethodsFrom terra nrow ncol
#' @param x object to check dimensions of
NULL






## nrow() S4 generic ####

if(!isGeneric('nrow')) setOldClass('nrow')
if(!isGeneric('ncol')) setOldClass('ncol')
if(!isGeneric('dim')) setOldClass('dim')

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
setMethod('nrow', signature('spatLocsObj'), function(x) nrow(x@coordinates))


# TODO
# setMethod('dims', signature('coordDataMT'), function(x) nrow(x@coordinates))

#' @describeIn dims-generic Find rows of giotto S4s with data.table based \code{coordinates} slots
#' @export
setMethod('nrow', signature('exprData'), function(x) nrow(x@exprMat))

#' @describeIn dims-generic Find rows of giotto S4s with data.table based \code{coordinates} slots
#' @export
setMethod('nrow', signature('metaData'), function(x) nrow(x@metaDT))

#' @describeIn dims-generic Find rows of spatialNetworkObj
#' @export
setMethod('nrow', signature('spatialNetworkObj'), function(x) nrow(x@networkDT))

#' @rdname dims-generic
#' @export
setMethod('nrow', signature('enrData'), function(x) nrow(x@enrichDT))

## ncol() generic ####

#' @describeIn dims-generic Find cols of giotto S4s with Matrix based \code{exprMat} slots
#' @export
setMethod('ncol', signature('exprData'), function(x) ncol(x@exprMat))

#' @describeIn dims-generic Find cols of giotto S4s with data.table based \code{metaDT} slots
#' @export
setMethod('ncol', signature('metaData'), function(x) ncol(x@metaDT))

setMethod('nrow', signature('enrData'), function(x) nrow(x@enrichDT))

## dim() generic ####

#' @describeIn dims-generic Find dimensions of giotto S4s with Matrix based \code{exprMat} slots
#' @export
setMethod('dim', signature('exprData'), function(x) dim(x@exprMat))

#' @describeIn dims-generic Find dimensions of giotto S4s with data.table based \code{metaDT} slots
#' @export
setMethod('dim', signature('metaData'), function(x) dim(x@metaDT))

#' @rdname dims-generic
#' @export
setMethod('dim', signature('enrData'), function(x) dim(x@enrichDT))








## colnames and rownames generics ####
#' @title Row and column names
#' @name row-plus-colnames-generic
#' @aliases colnames rownames
#' @description Retrieve or set the row or column names of an object
#' @param x object
#' @return A character vector of row or col names
if(!isGeneric('colnames')) setOldClass('colnames')
if(!isGeneric('rownames')) setOldClass('rownames')


#' @rdname row-plus-colnames-generic
#' @export
setMethod('colnames', signature(x = 'exprObj'), function(x) colnames(x[]))

#' @rdname row-plus-colnames-generic
#' @export
setMethod('rownames', signature(x = 'exprObj'), function(x) rownames(x[]))







# terra-based object serialization ####
## wrap() generic ####

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
            pgp@unique_ID_cache = x@unique_ID_cache
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
            pgp@unique_ID_cache = x@unique_ID_cache
            pgp@packed_spatVector = terra::wrap(x@spatVector)
            pgp@networks = x@networks
            return(pgp)
          }
)


## unwrap methods ####
# For compatibility before terra 1.6.41, vect will be used

#' @describeIn wrap-generic Unwrap giottoPolygon
#' @importMethodsFrom terra vect
#' @export
setMethod('vect', signature(x = 'packedGiottoPolygon'),
          function(x) {
            gp = new('giottoPolygon')
            gp@name = x@name
            gp@spatVector = terra::vect(x@packed_spatVector)

            # new cache slot
            if(!is.null(attr(x, 'unique_ID_cache'))) {
              gp@unique_ID_cache = x@unique_ID_cache
            } else gp@unique_ID_cache = spatIDs(gp)

            if(!is.null(x@packed_spatVectorCentroids)) {
              gp@spatVectorCentroids = terra::vect(x@packed_spatVectorCentroids)
            }
            if(length(x@packed_overlaps) > 0) {
              gp@overlaps = lapply(x@packed_overlaps, function(sv) {
                if(inherits(sv, 'PackedSpatVector')) {
                  terra::vect(sv)
                } else {
                  sv
                }
              })
            }
            return(gp)
          }
)



#' @describeIn wrap-generic Unwrap giottoPolygon
#' @importMethodsFrom terra vect
#' @export
setMethod('vect', signature(x = 'packedGiottoPoints'),
          function(x) {
            gp = new('giottoPoints')
            gp@feat_type = x@feat_type
            gp@spatVector = terra::vect(x@packed_spatVector)

            # new cache slot
            if(!is.null(attr(x, 'unique_ID_cache'))) {
              gp@unique_ID_cache = x@unique_ID_cache
            } else gp@unique_ID_cache = featIDs(gp)

            gp@networks = x@networks
            return(gp)
            }
)


#' @describeIn wrap-generic Unwrap giotto
#' @importMethodsFrom terra vect
#' @export
setMethod('vect', signature(x = 'packedGiotto'),
          function(x) {
            gobj = new('giotto')
            g_slots = methods::slotNames('giotto')
            g_slots = g_slots[!g_slots %in% c('spatial_info', 'feat_info')]
            for(g_slot in g_slots) {
              slot(gobj, g_slot) = slot(x, g_slot)
            }
            gobj@spatial_info = lapply(x@packed_spatial_info, vect)
            gobj@feat_info = lapply(x@packed_feat_info, vect)
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


if(!isGeneric('rbind')) setGeneric('rbind', signature = '...')

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
NULL

#' @describeIn plot-generic Plot \emph{magick}-based giottoImage object. ... param passes to \code{\link{plot_giottoImage_MG}}
#' @export
setMethod('plot', signature(x = 'giottoImage', y = 'missing'), function(x,y,...) plot_giottoImage_MG(giottoImage = x,...))

#' @describeIn plot-generic Plot \emph{terra}-based giottoLargeImage object. ... param passes to \code{\link{plot_giottoLargeImage}}
#' @importMethodsFrom terra plot
#' @export
setMethod('plot', signature(x = 'giottoLargeImage', y = 'missing'), function(x,y,...) plot_giottoLargeImage(giottoLargeImage = x,...))

#' @describeIn plot-generic Plot \emph{terra}-based giottoPolygon object. ... param passes to \code{\link[terra]{plot}}
#' @importMethodsFrom terra plot
#' @param point_size size of points when plotting giottoPolygon object centroids
#' @param type what to plot: either 'poly' (default) or polygon 'centroid'
#' @export
setMethod('plot', signature(x = 'giottoPolygon', y = 'missing'),
          function(x, point_size = 0.1, type = c('poly', 'centroid'), ...) {
            plot_giotto_polygon(x = x, point_size = point_size, type = type, ...)
          })

#' @describeIn plot-generic \emph{terra}-based giottoPoint object. ... param passes to \code{\link[terra]{plot}}
#' @param point_size size of points when plotting giottoPoints
#' @param feats specific features to plot within giottoPoints object (defaults to NULL, meaning all available features)
#' @param raster default = TRUE, whether to plot points as rasterized plot with
#' size based on \code{size} param
#' @param raster_size Default is 600. Only used when \code{raster} is TRUE
#' @export
setMethod('plot', signature(x = 'giottoPoints', y = 'missing'),
          function(x, point_size = 0, feats = NULL, raster = TRUE, raster_size = 600, ...) {
            plot_giotto_points(x = x, point_size = point_size, feats = feats,
                               raster = raster, raster_size = raster_size, ...)
          })


#' @describeIn plot-generic Plot a spatLocsObj
#' @export
setMethod('plot', signature(x = 'spatLocsObj', y = 'missing'), function(x, ...) {
  l = list(...)
  if(is.null(l$asp)) l$asp = 1
  if(is.null(l$xlab)) l$xlab = ''
  if(is.null(l$ylab)) l$ylab = ''
  if(is.null(l$cex)) l$cex = 0.5
  if(nrow(x) > 10000L) {
    if(is.null(l$pch)) l$pch = '.'
  }

  do.call('plot', append(l, list(x = x[]$sdimx, y = x[]$sdimy)))
})


#' @describeIn plot-generic Plot a spatialNetworkObj
#' @export
setMethod('plot', signature(x = 'spatialNetworkObj', y = 'missing'), function(x, ...) {
  l = list(...)
  if(is.null(l$asp)) l$asp = 1
  if(is.null(l$xlab)) l$xlab = ''
  if(is.null(l$ylab)) l$ylab = ''
  if(is.null(l$cex)) l$cex = 0.5
  if(is.null(l$col)) {
    line_col = 'red'
  } else {
    line_col = l$col
    l$col = NULL
  }
  if(is.null(l$lwd)) {
    line_width = 1L
  } else {
    line_width = l$lwd
    l$lwd = NULL
  }
  if(is.null(l$lty)) {
    line_type = 1L
  } else {
    line_type = l$lty
    l$lty = NULL
  }
  # find nodes
  nodes = unique(rbind(x[][, c('sdimx_begin', 'sdimy_begin')],
                       x[][, c('sdimx_end', 'sdimy_end')],
                       use.names = FALSE))
  if(nrow(nodes) > 10000L) {
    if(is.null(l$pch)) l$pch = '.'
  }
  do.call('plot', append(l, list(x = nodes$sdimx_begin, y = nodes$sdimy_begin)))
  segments(x0 = x[]$sdimx_begin, y0 = x[]$sdimy_begin,
           x1 = x[]$sdimx_end, y1 = x[]$sdimy_end,
           col = line_col, lty = line_type, lwd = line_width)
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
if(!isGeneric('copy')) setGeneric('copy', function(x) {
  standardGeneric('copy')
})

#' @describeIn copy-generic Copy \emph{data.table}-based spatial locations object.
#' @export
setMethod('copy', signature(x = 'coordDataDT'), function(x) {
  x@coordinates = data.table::copy(x@coordinates)
  x
})






# Giotto subnesting ####

## prov() S4 generic ####
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




## spatUnit() S4 generic ####
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
  value = as.character(value)
  x@spat_unit = value
  x
})




## featType() S4 generic ####
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
  value = as.character(value)
  x@feat_type = value
  x
})




## objName() generic ####
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
  value = as.character(value)
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


# setMethod('[', signature(x = 'giotto', i = 'character', j = 'missing', drop = 'missing'),
#           function(x, i, spat_unit = NULL, feat_type = NULL, name = NULL) {
#
#             # set defaults
#             spat_unit = set_default_spat_unit(gobject = x,
#                                               spat_unit = spat_unit)
#             feat_type = set_default_feat_type(gobject = x,
#                                               spat_unit = spat_unit,
#                                               feat_type = feat_type)
#             if(is.null(name)) name = 1L
#
#             switch(i,
#                    'expression' = slot(x, 'expression')[[spat_unit]][[feat_type]][[name]],
#                    'spatial_locs' = slot(x, 'spatial_locs')[[spat_unit]][[name]],
#                    'spatial_info' = slot(x, 'spatial_info')[[name]],
#                    'feat_info' = slot(x, 'feat_info')[[feat_type]])
#           })

# setMethod('[', signature(x = 'giotto', i = 'character', j = 'numeric', drop = 'missing'),
#           function(x, i, j, spat_unit = NULL, feat_type = NULL, name = NULL) {
#
#             avail_data = switch(i,
#                                 'expression' = list_expression(gobject = x,
#                                                                spat_unit = spat_unit,
#                                                                feat_type = feat_type),
#                                 'spatial_locs' = list_spatial_locations(gobject = x,
#                                                                         spat_unit = spat_unit),
#                                 'spatial_info' = list_spatial_info(gobject = x),
#                                 'feat_info' = list_feature_info(gobject = x))
#
#             switch(i,
#                    'expression' = get_expression_values(gobject = x,
#                                                         spat_unit = avail_data$spat_unit[[j]],
#                                                         feat_type = avail_data$feat_type[[j]],
#                                                         values = avail_data$name[[j]],
#                                                         output = 'exprObj',
#                                                         set_defaults = FALSE))
#
#           })

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

# * giottoPoints ####
#' @rdname extract-methods
#' @section \code{`[`} methods:
#'   Return \code{giottoPoints} spatVector slot
#' @export
setMethod('[', signature(x = 'giottoPoints', i = 'missing', j = 'missing', drop = 'missing'),
          function(x, i, j) {
            x@spatVector
          })

#' @rdname extract-methods
#' @aliases [<-,giottoPoints,missing,missing,ANY-method [<-,giottoPoints,missing,missing-method
#' @docType methods
#' @section \code{`[<-`} methods:
#'   Assign to \code{spatVector} slot in giotto S4
#' @export
setMethod('[<-', signature(x = 'giottoPoints', i = 'missing', j = 'missing', value = 'ANY'),
          function(x, i, j, value) {
            x@spatVector = value
            x
          })

# * giottoPolygon ####
#' @rdname extract-methods
#' @section \code{`[`} methods:
#'   Return \code{giottoPolygon} spatVector slot
#' @export
setMethod('[', signature(x = 'giottoPolygon', i = 'missing', j = 'missing', drop = 'missing'),
          function(x, i, j) {
            x@spatVector
          })

#' @rdname extract-methods
#' @aliases [<-,giottoPolygon,missing,missing,ANY-method [<-,giottoPolygon,missing,missing-method
#' @docType methods
#' @section \code{`[<-`} methods:
#'   Assign to \code{spatVector} slot in giottoPolygon
#' @export
setMethod('[<-', signature(x = 'giottoPolygon', i = 'missing', j = 'missing', value = 'ANY'),
          function(x, i, j, value) {
            x@spatVector = value
            x
          })












