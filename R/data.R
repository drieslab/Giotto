
#' mini Giotto object for MERFISH data (Vizgen) at subcellular resolution
#'
#' Mini Giotto object created from the Mouse Brain Vizgen data.
#' @name gobject_mini_vizgen
#'
#' @docType data
#'
#' @usage loadGiottoMini(dataset = 'vizgen')
#'
#' @format An object of class \code{"giotto"};
#' see \code{\link[Giotto]{createGiottoObject}} or
#'  \code{\link[Giotto]{createGiottoObjectSubcellular}}.
#'
#' @keywords datasets
#'
#' @references www.vizgen.com
#'
#' @examples
#'
#' \dontrun{
#'
#' gvizg = loadGiottoMini(dataset = 'vizgen')
#'
#' spatDimPlot2D(gobject = gvizg,
#'              spat_unit = 'z1',
#'              show_image = F,
#'              largeImage_name = 'polyT_z1',
#'              cell_color = 'leiden_clus',
#'              spat_point_size = 2,
#'              dim_point_size = 2)}
NULL
