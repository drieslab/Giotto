# `spatPlot3D`

spatPlot3D


## Description

Visualize cells according to spatial coordinates


## Usage

```r
spatPlot3D(
  gobject,
  spat_unit = NULL,
  feat_type = NULL,
  sdimx = "sdimx",
  sdimy = "sdimy",
  sdimz = "sdimz",
  spat_enr_names = NULL,
  point_size = 3,
  cell_color = NULL,
  cell_color_code = NULL,
  select_cell_groups = NULL,
  select_cells = NULL,
  show_other_cells = T,
  other_cell_color = "lightgrey",
  other_point_size = 0.5,
  other_cell_alpha = 0.5,
  show_network = F,
  spatial_network_name = "Delaunay_network",
  network_color = NULL,
  network_alpha = 1,
  show_grid = F,
  spatial_grid_name = "spatial_grid",
  grid_color = NULL,
  grid_alpha = 1,
  title = "",
  show_legend = T,
  axis_scale = c("cube", "real", "custom"),
  custom_ratio = NULL,
  x_ticks = NULL,
  y_ticks = NULL,
  z_ticks = NULL,
  show_plot = NA,
  return_plot = NA,
  save_plot = NA,
  save_param = list(),
  default_save_name = "spat3D"
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object
`spat_unit`     |     spatial unit
`feat_type`     |     feature type
`sdimx`     |     x-axis dimension name (default = 'sdimx')
`sdimy`     |     y-axis dimension name (default = 'sdimy')
`sdimz`     |     z-axis dimension name (default = 'sdimy')
`spat_enr_names`     |     names of spatial enrichment results to include
`point_size`     |     size of point (cell)
`cell_color`     |     color for cells (see details)
`cell_color_code`     |     named vector with colors
`select_cell_groups`     |     select subset of cells/clusters based on cell_color parameter
`select_cells`     |     select subset of cells based on cell IDs
`show_other_cells`     |     display not selected cells
`other_cell_color`     |     color of not selected cells
`other_point_size`     |     size of not selected cells
`other_cell_alpha`     |     alpha of not selected cells
`show_network`     |     show underlying spatial network
`spatial_network_name`     |     name of spatial network to use
`network_color`     |     color of spatial network
`network_alpha`     |     opacity of spatial network
`show_grid`     |     show spatial grid
`spatial_grid_name`     |     name of spatial grid to use
`grid_color`     |     color of spatial grid
`grid_alpha`     |     opacity of spatial grid
`title`     |     title of plot
`show_legend`     |     show legend
`axis_scale`     |     the way to scale the axis
`custom_ratio`     |     customize the scale of the plot
`x_ticks`     |     set the number of ticks on the x-axis
`y_ticks`     |     set the number of ticks on the y-axis
`z_ticks`     |     set the number of ticks on the z-axis
`show_plot`     |     show plot
`return_plot`     |     return ggplot object
`save_plot`     |     directly save the plot [boolean]
`save_param`     |     list of saving parameters, see [`showSaveParameters`](#showsaveparameters)
`default_save_name`     |     default save name for saving, don't change, change save_name in save_param


## Value

ggplot


## Seealso

Other spatial visualizations:
 [`spatPlot2D`](#spatplot2d) ,
 [`spatPlot`](#spatplot)


