# Methods and Generics ####



# nrow() S4 generic ####

# Define base::nrow() for giottoPoints and giottoPolygon objects

#' @title Dimensions of giotto points and polygon objects
#' @name nrow-generic
#' @description Get the number of rows (nrow)
#' @include classes.R
#' @param x giottoPoints or giottoPolygon object
#' @aliases nrow
#' @exportMethod nrow
setOldClass('nrow')

# setMethod('nrow', signature = 'giotto', function(x) {
#   avail_exp = list_expression(x)
#   get_expression_values(x,
#                         spat_unit = avail_exp$spat_unit[1],
#                         feat_type = avail_exp$feat_type[1],
#                         values = avail_exp$name[1]) %>%
#     nrow()
# })

#' @describeIn nrow-generic Find rows of giottoPoints object
#' @export
setMethod('nrow', signature('giottoPoints'), function(x) terra::nrow(x@spatVector))


#' @describeIn nrow-generic Find rows of giottoPolygon object
#' @export
setMethod('nrow', signature('giottoPolygon'), function(x) terra::nrow(x@spatVector))





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


