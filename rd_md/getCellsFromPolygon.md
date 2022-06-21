# `getCellsFromPolygon`

Get cells located within the polygons area


## Description

Get cells located within the polygons area


## Usage

```r
getCellsFromPolygon(
  gobject,
  polygon_slot = "spatial_info",
  cells_loc_slot = "spatial_locs"
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     A Giotto object
`polygon_slot`     |     Slot name where polygon coordinates are stored in Giotto object
`cells_loc_slot`     |     Slot name where cell coordinates are stored in Giotto object


## Value

A `SpatVector` with cell IDs, x,y coordinates, and polygon name where
 each cell is located in.


## Examples

```r
## Plot interactive polygons
my_spatPlot <- spatPlot2D(gobject = my_giotto_object,
show_image = TRUE,
point_alpha = 0.75,
save_plot = FALSE)
my_polygon_coords <- plotInteractivePolygons(my_spatPlot)

## Add polygon coordinates to Giotto object
my_giotto_polygons <- createGiottoPolygonsFromDfr(my_polygon_coords)
my_giotto_object <- addGiottoPolygons(gobject = my_giotto_object,
gpolygons = list(my_giotto_polygons))

## Get cells located within polygons area
getCellsFromPolygon(my_giotto_object)
```


