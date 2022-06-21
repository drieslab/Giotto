# `cellProximitySpatPlot2D`

cellProximitySpatPlot2D


## Description

Visualize 2D cell-cell interactions according to spatial coordinates in ggplot mode


## Usage

```r
cellProximitySpatPlot2D(
  gobject,
  spat_unit = NULL,
  feat_type = NULL,
  spat_loc_name = NULL,
  interaction_name = NULL,
  cluster_column = NULL,
  sdimx = "sdimx",
  sdimy = "sdimy",
  cell_color = NULL,
  cell_color_code = NULL,
  color_as_factor = T,
  show_other_cells = F,
  show_network = F,
  show_other_network = F,
  network_color = NULL,
  spatial_network_name = "Delaunay_network",
  show_grid = F,
  grid_color = NULL,
  spatial_grid_name = "spatial_grid",
  coord_fix_ratio = 1,
  show_legend = T,
  point_size_select = 2,
  point_select_border_col = "black",
  point_select_border_stroke = 0.05,
  point_size_other = 1,
  point_alpha_other = 0.3,
  point_other_border_col = "lightgrey",
  point_other_border_stroke = 0.01,
  show_plot = NA,
  return_plot = NA,
  save_plot = NA,
  save_param = list(),
  default_save_name = "cellProximitySpatPlot2D"
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object
`spat_unit`     |     spatial unit
`feat_type`     |     feature type
`spat_loc_name`     |     spatial locations to use
`interaction_name`     |     cell-cell interaction name
`cluster_column`     |     cluster column with cell clusters
`sdimx`     |     x-axis dimension name (default = 'sdimx')
`sdimy`     |     y-axis dimension name (default = 'sdimy')
`cell_color`     |     color for cells (see details)
`cell_color_code`     |     named vector with colors
`color_as_factor`     |     convert color column to factor
`show_other_cells`     |     decide if show cells not in network
`show_network`     |     show spatial network of selected cells
`show_other_network`     |     show spatial network of not selected cells
`network_color`     |     color of spatial network
`spatial_network_name`     |     name of spatial network to use
`show_grid`     |     show spatial grid
`grid_color`     |     color of spatial grid
`spatial_grid_name`     |     name of spatial grid to use
`coord_fix_ratio`     |     fix ratio between x and y-axis
`show_legend`     |     show legend
`point_size_select`     |     size of selected points
`point_select_border_col`     |     border color of selected points
`point_select_border_stroke`     |     stroke size of selected points
`point_size_other`     |     size of other points
`point_alpha_other`     |     opacity of other points
`point_other_border_col`     |     border color of other points
`point_other_border_stroke`     |     stroke size of other points
`show_plot`     |     show plots
`return_plot`     |     return ggplot object
`save_plot`     |     directly save the plot [boolean]
`save_param`     |     list of saving parameters from [`all_plots_save_function`](#allplotssavefunction)
`default_save_name`     |     default save name for saving, don't change, change save_name in save_param


## Details

Description of parameters.


## Value

ggplot


