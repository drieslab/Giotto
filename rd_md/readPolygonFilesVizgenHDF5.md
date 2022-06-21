# `readPolygonFilesVizgenHDF5`

readPolygonFilesVizgenHDF5


## Description

Read and create polygons for all cells, or for only selected FOVs.


## Usage

```r
readPolygonFilesVizgenHDF5(
  boundaries_path,
  fovs = NULL,
  polygon_feat_types = 0:6,
  custom_polygon_names = NULL,
  flip_x_axis = F,
  flip_y_axis = F,
  smooth_polygons = TRUE,
  smooth_vertices = 60,
  set_neg_to_zero = FALSE,
  verbose = TRUE
)
```


## Arguments

Argument      |Description
------------- |----------------
`boundaries_path`     |     path to the cell_boundaries folder
`fovs`     |     subset of fovs to use
`polygon_feat_types`     |     a vector containing the polygon feature types
`custom_polygon_names`     |     a character vector to provide custom polygon names (optional)
`flip_x_axis`     |     flip x axis of polygon coordinates (multiply by -1)
`flip_y_axis`     |     flip y axis of polygon coordinates (multiply by -1)
`smooth_polygons`     |     smooth polygons (default = TRUE)
`smooth_vertices`     |     number of vertices for smoothing
`set_neg_to_zero`     |     set negative values to zero when smoothing
`verbose`     |     be verbose


## Seealso

[`smoothGiottoPolygons`](#smoothgiottopolygons)


