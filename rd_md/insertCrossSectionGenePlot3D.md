# `insertCrossSectionGenePlot3D`

insertCrossSectionGenePlot3D


## Description

Visualize cells and gene expression in a virtual cross section according to spatial coordinates


## Usage

```r
insertCrossSectionGenePlot3D(
  gobject,
  spat_loc_name = "raw",
  crossSection_obj = NULL,
  name = NULL,
  spatial_network_name = "Delaunay_network",
  mesh_grid_color = "#1f77b4",
  mesh_grid_width = 3,
  mesh_grid_style = "dot",
  sdimx = "sdimx",
  sdimy = "sdimy",
  sdimz = "sdimz",
  show_other_cells = F,
  axis_scale = c("cube", "real", "custom"),
  custom_ratio = NULL,
  show_plot = NA,
  return_plot = NA,
  save_plot = NA,
  save_param = list(),
  default_save_name = "spatGenePlot3D_with_cross_section",
  ...
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object
`spat_loc_name`     |     name of spatial locations
`crossSection_obj`     |     cross section object as alternative input. default = NULL.
`name`     |     name of virtual cross section to use
`spatial_network_name`     |     name of spatial network to use
`mesh_grid_color`     |     color for the meshgrid lines
`mesh_grid_width`     |     width for the meshgrid lines
`mesh_grid_style`     |     style for the meshgrid lines
`sdimx`     |     x-axis dimension name (default = 'sdimx')
`sdimy`     |     y-axis dimension name (default = 'sdimy')
`sdimz`     |     z-axis dimension name (default = 'sdimy')
`show_other_cells`     |     display not selected cells
`axis_scale`     |     axis_scale
`custom_ratio`     |     custom_ratio
`show_plot`     |     show plots
`return_plot`     |     return ggplot object
`save_plot`     |     directly save the plot [boolean]
`save_param`     |     list of saving parameters from [`all_plots_save_function`](#allplotssavefunction)
`default_save_name`     |     default save name for saving, don't change, change save_name in save_param
`...`     |     parameters for spatGenePlot3D


## Details

Description of parameters.


## Value

ggplot


