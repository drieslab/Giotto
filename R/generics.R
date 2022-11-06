# Methods and Generics ####





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
#' @export
setMethod('plot', signature(x = 'giottoPolygon', y = 'missing'), function(x,y,...) terra::plot(x = x@spatVector,...))

#' @describeIn plot-generic \emph{terra}-based giottoPoint object. ... param passes to \code{\link[terra]{plot}}
#' @param point_size size of points when plotting giottoPoint
#' @param feats specific features to plot within giottoPoint object (defaults to NULL, meaning all available features)
#' @export
setMethod('plot', signature(x = 'giottoPoints', y = 'missing'), function(x,
                                                                         y,
                                                                         point_size = 0.1,
                                                                         feats = NULL,
                                                                         ...) {

  if(is.null(feats)) terra::plot(x = x@spatVector, cex = point_size,...)
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











#' @title Extract or replace parts of an object
#' @name extract-generic
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
#'   TODO
NULL

# $ S4 access generic ####

## * coordDataDT ####

#' @rdname extract-generic
#' @section \code{`$`} methods:
#'   Select by colname from giotto S4 data.table coordinates slot.
#' @export
setMethod('$', signature(x = 'coordDataDT'),
          function(x, name) x@coordinates[[name]])


#' @rdname extract-generic
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

#' @rdname extract-generic
#' @section \code{`$`} methods:
#'   Select by colname from giotto S4 data.table metaDT slot.
#' @export
setMethod('$', signature(x = 'metaData'),
          function(x, name) x@metaDT[[name]])


#' @rdname extract-generic
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

#' @rdname extract-generic
#' @section \code{`[`} methods:
#'   Select rows (i) and cols (j) from giotto S4 coordinates slot
#' @export
setMethod('[', signature(x = 'coordDataDT', i = 'missing', j = 'ANY', drop = 'missing'),
          function(x, i, j) {
            x@coordinates = x@coordinates[, j = j, with = FALSE]
            x
          })

#' @rdname extract-generic
#' @export
setMethod('[', signature(x = 'coordDataDT', i = 'ANY', j = 'missing', drop = 'missing'),
          function(x, i, j) {
            x@coordinates = x@coordinates[i = i,]
            x
          })

#' @rdname extract-generic
#' @export
setMethod('[', signature(x = 'coordDataDT', i = 'ANY', j = 'ANY', drop = 'missing'),
          function(x, i, j) {
            x@coordinates = x@coordinates[i = i, j = j, with = FALSE]
            x
          })

#' @rdname extract-generic
#' @section \code{`[`} methods:
#'   Return \code{coordinates} slot data.table from giotto S4
#' @export
setMethod('[', signature(x = 'coordDataDT', i = 'missing', j = 'missing', drop = 'missing'),
          function(x, i, j) {
            x@coordinates
          })

# setMethod("[[")

## * metaData ####

#' @rdname extract-generic
#' @section \code{`[`} methods:
#'   Select rows (i) and cols (j) from giotto S4 metaDT slot
#' @export
setMethod('[', signature(x = 'metaData', i = 'missing', j = 'ANY', drop = 'missing'),
          function(x, i, j) {
            x@metaDT = x@metaDT[, j = j, with = FALSE]
            x
          })

#' @rdname extract-generic
#' @export
setMethod('[', signature(x = 'metaData', i = 'ANY', j = 'missing', drop = 'missing'),
          function(x, i, j) {
            x@metaDT = x@metaDT[i = i,]
            x
          })

#' @rdname extract-generic
#' @export
setMethod('[', signature(x = 'metaData', i = 'ANY', j = 'ANY', drop = 'missing'),
          function(x, i, j) {
            x@metaDT = x@metaDT[i = i, j = j, with = FALSE]
            x
          })

#' @rdname extract-generic
#' @section \code{`[`} methods:
#'   Return \code{coordinates} slot data.table from giotto S4
#' @export
setMethod('[', signature(x = 'metaData', i = 'missing', j = 'missing', drop = 'missing'),
          function(x, i, j) {
            x@metaDT
          })


## * dimObj (temp) ####

#' @rdname extract-generic
#' @export
setMethod('[', signature(x = 'dimObj', i = 'ANY', j = 'ANY', drop = 'missing'),
          function(x, i, j) {
            x@coordinates = x@coordinates[i = i, j = j]
            x
          })

## * exprData ####

#' @rdname extract-generic
#' @section \code{`[`} methods:
#'   Select rows (i) and cols (j) from giotto S4 exprMat slot
#' @export
setMethod('[', signature(x = 'exprData', i = 'missing', j = 'ANY', drop = 'missing'),
          function(x, i, j) {
            x@exprMat = x@exprMat[, j = j]
            x
          })

#' @rdname extract-generic
#' @export
setMethod('[', signature(x = 'exprData', i = 'ANY', j = 'missing', drop = 'missing'),
          function(x, i, j) {
            x@exprMat = x@exprMat[i = i,]
            x
          })

#' @rdname extract-generic
#' @export
setMethod('[', signature(x = 'exprData', i = 'ANY', j = 'ANY', drop = 'missing'),
          function(x, i, j) {
            x@exprMat = x@exprMat[i = i, j = j]
            x
          })

#' @rdname extract-generic
#' @section \code{`[`} methods:
#'   Return \code{exprMat} slot Matrix object from giotto S4
#' @export
setMethod('[', signature(x = 'exprData', i = 'missing', j = 'missing', drop = 'missing'),
          function(x, i, j) {
            x@exprMat
          })






