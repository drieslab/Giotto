# `addCellsFromPolygon`

Add corresponding polygon IDs to cell metadata


## Description

Add corresponding polygon IDs to cell metadata


## Usage

```r
addCellsFromPolygon(gobject, cellsFromPolygon, feat_type = "rna")
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     A Giotto object
`cellsFromPolygon`     |     A `SpatVector` with cell IDs located inside each polygon
`feat_type`     |     feature name where metadata will be added


## Value

A Giotto object with a modified cell_metadata slot that includes the
 polygon name where each cell is located or NA if the cell is not located
 within a polygon area


## Examples

```r
## Plot interactive polygons
my_polygon_coords <- plotInteractivePolygons(my_spatPlot)

## Add polygon coordinates to Giotto object
my_giotto_polygons <- createGiottoPolygonsFromDfr(my_polygon_coords)
my_giotto_object <- addGiottoPolygons(gobject = my_giotto_object,
gpolygons = list(my_giotto_polygons))

## Get cells located within polygons area
my_polygon_cells <- getCellsFromPolygon(my_giotto_object)

my_giotto_object <- addCellsFromPolygon(my_giotto_object, my_polygon_cells)
```


