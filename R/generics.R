# Methods and Generics ####

# Define base::nrow() dispatch for gobjects
# setOldClass('nrow')
# setMethod('nrow', signature = 'giotto', function(x) {
#   avail_exp = list_expression(x)
#   get_expression_values(x,
#                         feat_type = avail_exp$feat_type[1],
#                         spat_unit = avail_exp$spat_unit[1],
#                         values = avail_exp$name[1]) %>%
#     nrow()
# })

#' "plot" S4 method
#' 
#' @name plot-method
#' @aliases plot
#' @family plot
#' @exportMethod plot
setOldClass('plot')

#' @export
setMethod('plot', signature('giottoImage'), function(x,...) plot_giottoImage_MG(giottoImage = x,...))

#' @export
setMethod('plot', signature('giottoLargeImage'), function(x,...) plot_giottoLargeImage(giottoLargeImage = x,...))

#' @export
setMethod('plot', signature('giottoPolygon'), function(x,...) {
  terra::plot(x = x@spatVector,...)
})

#' @export
setMethod('plot', signature('giottoPoints'), function(x,point_size = 0.1,...) {
  terra::plot(x = x@spatVector,cex = point_size,...)
})


