% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/interactivity.R
\name{addPolygonCells}
\alias{addPolygonCells}
\title{Add corresponding polygon IDs to cell metadata}
\usage{
addPolygonCells(
  gobject,
  polygon_name = "selections",
  spat_unit = "cell",
  spat_loc_name = "raw",
  feat_type = "rna",
  polygons = NULL,
  na.label = "no_polygon"
)
}
\arguments{
\item{gobject}{A Giotto object}

\item{polygon_name}{name of polygon selections}

\item{spat_unit}{spatial unit}

\item{spat_loc_name}{name of spatial locations to use}

\item{feat_type}{feature name where metadata will be added}

\item{polygons}{polygon names to plot (e.g. 'polygon_1'). If NULL, plots all available polygons}

\item{na.label}{polygon label for cells located outside of polygons area. Default = "no_polygon"}
}
\value{
A Giotto object with a modified cell_metadata slot that includes the
polygon name where each cell is located or no_polygon label if the cell is not located
within a polygon area
}
\description{
Add corresponding polygon IDs to cell metadata
}
\examples{

\dontrun{
## Plot interactive polygons
my_polygon_coords <- plotInteractivePolygons(my_spatPlot)

## Add polygon coordinates to Giotto object
my_giotto_polygons <- createGiottoPolygonsFromDfr(my_polygon_coords)
my_giotto_object <- addGiottoPolygons(gobject = my_giotto_object,
                                      gpolygons = list(my_giotto_polygons))

## Add polygon IDs to cell metadata
my_giotto_object <- addPolygonCells(my_giotto_object)
}

}
