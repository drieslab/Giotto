# `plotStatDelaunayNetwork`

plotStatDelaunayNetwork


## Description

Plots network statistics for a Delaunay network..


## Usage

```r
plotStatDelaunayNetwork(
  gobject,
  feat_type = NULL,
  spat_unit = NULL,
  method = c("deldir", "delaunayn_geometry", "RTriangle"),
  dimensions = "all",
  maximum_distance = "auto",
  minimum_k = 0,
  options = "Pp",
  Y = TRUE,
  j = TRUE,
  S = 0,
  show_plot = NA,
  return_plot = NA,
  save_plot = NA,
  save_param = list(),
  default_save_name = "plotStatDelaunayNetwork",
  ...
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object
`feat_type`     |     feature type
`spat_unit`     |     spatial unit
`method`     |     package to use to create a Delaunay network
`dimensions`     |     which spatial dimensions to use (maximum 2 dimensions)
`maximum_distance`     |     distance cuttof for Delaunay neighbors to consider
`minimum_k`     |     minimum neigbhours if maximum_distance != NULL
`options`     |     (geometry) String containing extra control options for the underlying Qhull command; see the Qhull documentation (../doc/qhull/html/qdelaun.html) for the available options. (default = 'Pp', do not report precision problems)
`Y`     |     (RTriangle) If TRUE prohibits the insertion of Steiner points on the mesh boundary.
`j`     |     (RTriangle) If TRUE jettisons vertices that are not part of the final triangulation from the output.
`S`     |     (RTriangle) Specifies the maximum number of added Steiner points.
`show_plot`     |     show plots
`return_plot`     |     return ggplot object
`save_plot`     |     directly save the plot [boolean]
`save_param`     |     list of saving parameters, see [`showSaveParameters`](#showsaveparameters)
`default_save_name`     |     default save name for saving, don't change, change save_name in save_param
`list()`     |     Other parameters


## Value

giotto object with updated spatial network slot


