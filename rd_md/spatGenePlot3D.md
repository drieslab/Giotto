# `spatGenePlot3D`

spatGenePlot3D


## Description

Visualize cells and gene expression according to spatial coordinates


## Usage

```r
spatGenePlot3D(
  gobject,
  spat_unit = NULL,
  feat_type = NULL,
  spat_loc_name = NULL,
  expression_values = c("normalized", "scaled", "custom"),
  genes,
  show_network = FALSE,
  network_color = NULL,
  spatial_network_name = "Delaunay_network",
  edge_alpha = NULL,
  cluster_column = NULL,
  select_cell_groups = NULL,
  select_cells = NULL,
  show_other_cells = T,
  other_cell_color = "lightgrey",
  other_point_size = 1,
  genes_high_color = NULL,
  genes_mid_color = "white",
  genes_low_color = "blue",
  show_grid = FALSE,
  spatial_grid_name = "spatial_grid",
  point_size = 2,
  show_legend = TRUE,
  axis_scale = c("cube", "real", "custom"),
  custom_ratio = NULL,
  x_ticks = NULL,
  y_ticks = NULL,
  z_ticks = NULL,
  show_plot = NA,
  return_plot = NA,
  save_plot = NA,
  save_param = list(),
  default_save_name = "spatGenePlot3D"
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object
`spat_unit`     |     spatial unit
`feat_type`     |     feature type
`spat_loc_name`     |     name of spatial locations to use
`expression_values`     |     gene expression values to use
`genes`     |     genes to show
`show_network`     |     show underlying spatial network
`network_color`     |     color of spatial network
`spatial_network_name`     |     name of spatial network to use
`edge_alpha`     |     alpha of edges
`cluster_column`     |     cluster column to select groups
`select_cell_groups`     |     select subset of cells/clusters based on cell_color parameter
`select_cells`     |     select subset of cells based on cell IDs
`show_other_cells`     |     display not selected cells
`other_cell_color`     |     color of not selected cells
`other_point_size`     |     size of not selected cells
`genes_high_color`     |     color represents high gene expression
`genes_mid_color`     |     color represents middle gene expression
`genes_low_color`     |     color represents low gene expression
`show_grid`     |     show spatial grid
`spatial_grid_name`     |     name of spatial grid to use
`point_size`     |     size of point (cell)
`show_legend`     |     show legend
`axis_scale`     |     the way to scale the axis
`custom_ratio`     |     customize the scale of the plot
`x_ticks`     |     set the number of ticks on the x-axis
`y_ticks`     |     set the number of ticks on the y-axis
`z_ticks`     |     set the number of ticks on the z-axis
`show_plot`     |     show plots
`return_plot`     |     return ggplot object
`save_plot`     |     directly save the plot [boolean]
`save_param`     |     list of saving parameters, see [`showSaveParameters`](#showsaveparameters)
`default_save_name`     |     default save name for saving, don't change, change save_name in save_param


## Details

Description of parameters.


## Value

ggplot


## Seealso

Other spatial gene expression visualizations:
 [`spatGenePlot2D`](#spatgeneplot2d) ,
 [`spatGenePlot`](#spatgeneplot)


