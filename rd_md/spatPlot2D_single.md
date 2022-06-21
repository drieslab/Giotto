# `spatPlot2D_single`

spatPlot2D_single


## Description

Visualize cells according to spatial coordinates


## Usage

```r
spatPlot2D_single(
  gobject,
  feat_type = NULL,
  spat_unit = NULL,
  show_image = F,
  gimage = NULL,
  image_name = NULL,
  largeImage_name = NULL,
  spat_loc_name = NULL,
  sdimx = "sdimx",
  sdimy = "sdimy",
  spat_enr_names = NULL,
  cell_color = NULL,
  color_as_factor = T,
  cell_color_code = NULL,
  cell_color_gradient = c("blue", "white", "red"),
  gradient_midpoint = NULL,
  gradient_limits = NULL,
  select_cell_groups = NULL,
  select_cells = NULL,
  point_shape = c("border", "no_border", "voronoi"),
  point_size = 3,
  point_alpha = 1,
  point_border_col = "black",
  point_border_stroke = 0.1,
  show_cluster_center = F,
  show_center_label = F,
  center_point_size = 4,
  center_point_border_col = "black",
  center_point_border_stroke = 0.1,
  label_size = 4,
  label_fontface = "bold",
  show_network = F,
  spatial_network_name = "Delaunay_network",
  network_color = NULL,
  network_alpha = 1,
  show_grid = F,
  spatial_grid_name = "spatial_grid",
  grid_color = NULL,
  show_other_cells = T,
  other_cell_color = "lightgrey",
  other_point_size = 1,
  other_cells_alpha = 0.1,
  coord_fix_ratio = NULL,
  title = NULL,
  show_legend = T,
  legend_text = 8,
  legend_symbol_size = 1,
  background_color = "white",
  vor_border_color = "white",
  vor_max_radius = 200,
  vor_alpha = 1,
  axis_text = 8,
  axis_title = 8,
  show_plot = NA,
  return_plot = NA,
  save_plot = NA,
  verbose = FALSE,
  save_param = list(),
  default_save_name = "spatPlot2D_single"
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object
`feat_type`     |     feature type
`spat_unit`     |     spatial unit
`show_image`     |     show a tissue background image
`gimage`     |     a giotto image
`image_name`     |     name of a giotto image
`largeImage_name`     |     name of a giottoLargeImage
`spat_loc_name`     |     name of spatial locations
`sdimx`     |     x-axis dimension name (default = 'sdimx')
`sdimy`     |     y-axis dimension name (default = 'sdimy')
`spat_enr_names`     |     names of spatial enrichment results to include
`cell_color`     |     color for cells (see details)
`color_as_factor`     |     convert color column to factor
`cell_color_code`     |     named vector with colors
`cell_color_gradient`     |     vector with 3 colors for numeric data
`gradient_midpoint`     |     midpoint for color gradient
`gradient_limits`     |     vector with lower and upper limits
`select_cell_groups`     |     select subset of cells/clusters based on cell_color parameter
`select_cells`     |     select subset of cells based on cell IDs
`point_shape`     |     shape of points (border, no_border or voronoi)
`point_size`     |     size of point (cell)
`point_alpha`     |     transparancy of point
`point_border_col`     |     color of border around points
`point_border_stroke`     |     stroke size of border around points
`show_cluster_center`     |     plot center of selected clusters
`show_center_label`     |     plot label of selected clusters
`center_point_size`     |     size of center points
`label_size`     |     size of labels
`label_fontface`     |     font of labels
`show_network`     |     show underlying spatial network
`spatial_network_name`     |     name of spatial network to use
`network_color`     |     color of spatial network
`network_alpha`     |     alpha of spatial network
`show_grid`     |     show spatial grid
`spatial_grid_name`     |     name of spatial grid to use
`grid_color`     |     color of spatial grid
`show_other_cells`     |     display not selected cells
`other_cell_color`     |     color of not selected cells
`other_point_size`     |     point size of not selected cells
`other_cells_alpha`     |     alpha of not selected cells
`coord_fix_ratio`     |     fix ratio between x and y-axis
`title`     |     title of plot
`show_legend`     |     show legend
`legend_text`     |     size of legend text
`legend_symbol_size`     |     size of legend symbols
`background_color`     |     color of plot background
`vor_border_color`     |     border colorr for voronoi plot
`vor_max_radius`     |     maximum radius for voronoi 'cells'
`vor_alpha`     |     transparancy of voronoi 'cells'
`axis_text`     |     size of axis text
`axis_title`     |     size of axis title
`show_plot`     |     show plot
`return_plot`     |     return ggplot object
`save_plot`     |     directly save the plot [boolean]
`verbose`     |     be verbose
`save_param`     |     list of saving parameters from [`all_plots_save_function`](#allplotssavefunction)
`default_save_name`     |     default save name for saving, don't change, change save_name in save_param


## Details

Description of parameters.


## Value

ggplot


## Seealso

[`spatPlot3D`](#spatplot3d)


