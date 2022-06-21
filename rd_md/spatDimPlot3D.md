# `spatDimPlot3D`

spatDimPlot3D


## Description

Visualize cells according to spatial AND dimension reduction coordinates in plotly mode


## Usage

```r
spatDimPlot3D(
  gobject,
  spat_unit = NULL,
  feat_type = NULL,
  plot_alignment = c("horizontal", "vertical"),
  dim_reduction_to_use = "umap",
  dim_reduction_name = "umap",
  dim1_to_use = 1,
  dim2_to_use = 2,
  dim3_to_use = 3,
  spat_loc_name = NULL,
  sdimx = "sdimx",
  sdimy = "sdimy",
  sdimz = "sdimz",
  spat_enr_names = NULL,
  show_NN_network = FALSE,
  nn_network_to_use = "sNN",
  network_name = "sNN.pca",
  nn_network_color = "lightgray",
  nn_network_alpha = 0.5,
  show_cluster_center = F,
  show_center_label = T,
  center_point_size = 4,
  label_size = 16,
  select_cell_groups = NULL,
  select_cells = NULL,
  show_other_cells = T,
  other_cell_color = "lightgrey",
  other_point_size = 1.5,
  cell_color = NULL,
  color_as_factor = T,
  cell_color_code = NULL,
  dim_point_size = 3,
  show_spatial_network = F,
  spatial_network_name = "Delaunay_network",
  spatial_network_color = "lightgray",
  spatial_network_alpha = 0.5,
  show_spatial_grid = F,
  spatial_grid_name = "spatial_grid",
  spatial_grid_color = NULL,
  spatial_grid_alpha = 0.5,
  spatial_point_size = 3,
  axis_scale = c("cube", "real", "custom"),
  custom_ratio = NULL,
  x_ticks = NULL,
  y_ticks = NULL,
  z_ticks = NULL,
  legend_text_size = 12,
  show_plot = NA,
  return_plot = NA,
  save_plot = NA,
  save_param = list(),
  default_save_name = "spatDimPlot3D"
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object
`spat_unit`     |     spatial unit
`feat_type`     |     feature type
`plot_alignment`     |     direction to align plot
`dim_reduction_to_use`     |     dimension reduction to use
`dim_reduction_name`     |     dimension reduction name
`dim1_to_use`     |     dimension to use on x-axis
`dim2_to_use`     |     dimension to use on y-axis
`dim3_to_use`     |     dimension to use on z-axis
`spat_loc_name`     |     name for spatial locations
`sdimx`     |     = spatial dimension to use on x-axis
`sdimy`     |     = spatial dimension to use on y-axis
`sdimz`     |     = spatial dimension to use on z-axis
`spat_enr_names`     |     names of spatial enrichment results to include
`show_NN_network`     |     show underlying NN network
`nn_network_to_use`     |     type of NN network to use (kNN vs sNN)
`network_name`     |     name of NN network to use, if show_NN_network = TRUE
`nn_network_color`     |     color of nn network
`nn_network_alpha`     |     column to use for alpha of the edges
`show_cluster_center`     |     show the center of each cluster
`show_center_label`     |     provide a label for each cluster
`center_point_size`     |     size of the center point
`label_size`     |     size of the center label
`select_cell_groups`     |     select subset of cells/clusters based on cell_color parameter
`select_cells`     |     select subset of cells based on cell IDs
`show_other_cells`     |     display not selected cells
`other_cell_color`     |     color of not selected cells
`other_point_size`     |     size of not selected cells
`cell_color`     |     color for cells (see details)
`color_as_factor`     |     convert color column to factor
`cell_color_code`     |     named vector with colors
`dim_point_size`     |     size of points in dim. reduction space
`show_spatial_network`     |     show spatial network
`spatial_network_name`     |     name of spatial network to use
`spatial_network_color`     |     color of spatial network
`spatial_network_alpha`     |     alpha of spatial network
`show_spatial_grid`     |     show spatial grid
`spatial_grid_name`     |     name of spatial grid to use
`spatial_grid_color`     |     color of spatial grid
`spatial_grid_alpha`     |     alpha of spatial grid
`spatial_point_size`     |     size of spatial points
`axis_scale`     |     the way to scale the axis
`custom_ratio`     |     customize the scale of the plot
`x_ticks`     |     set the number of ticks on the x-axis
`y_ticks`     |     set the number of ticks on the y-axis
`z_ticks`     |     set the number of ticks on the z-axis
`legend_text_size`     |     size of legend
`show_plot`     |     show plot
`return_plot`     |     return ggplot object
`save_plot`     |     directly save the plot [boolean]
`save_param`     |     list of saving parameters, see [`showSaveParameters`](#showsaveparameters)
`default_save_name`     |     default save name for saving, don't change, change save_name in save_param


## Details

Description of parameters.


## Value

plotly


## Seealso

Other spatial and dimension reduction visualizations:
 [`spatDimPlot2D`](#spatdimplot2d) ,
 [`spatDimPlot`](#spatdimplot)


