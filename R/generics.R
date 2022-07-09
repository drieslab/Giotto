# Methods ####

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

# Set plotting generics
setOldClass('plot')
setMethod('plot', signature('giottoImage'), function(x,...) plot_giottoImage_MG(giottoImage = x,...))
setMethod('plot', signature('giottoLargeImage'), function(x,...) plot_giottoLargeImage(giottoLargeImage = x,...))
setMethod('plot', signature('giottoPolygon'), function(x) terra::plot(x@spatVector))
setMethod('plot', signature('giottoPoints'), function(x) terra::plot(x@spatVector))


