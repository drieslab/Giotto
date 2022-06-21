# `plot_spat_point_layer_ggplot`

plot_spat_point_layer_ggplot


## Description

creat ggplot point layer for spatial coordinates


## Usage

```r
plot_spat_point_layer_ggplot(
  ggobject,
  sdimx = NULL,
  sdimy = NULL,
  cell_locations_metadata_selected,
  cell_locations_metadata_other,
  cell_color = NULL,
  color_as_factor = T,
  cell_color_code = NULL,
  cell_color_gradient = c("blue", "white", "red"),
  gradient_midpoint = NULL,
  gradient_limits = NULL,
  select_cell_groups = NULL,
  select_cells = NULL,
  point_size = 2,
  point_alpha = 1,
  point_border_col = "lightgrey",
  point_border_stroke = 0.1,
  show_cluster_center = F,
  show_center_label = T,
  center_point_size = 4,
  center_point_border_col = "black",
  center_point_border_stroke = 0.1,
  label_size = 4,
  label_fontface = "bold",
  show_other_cells = T,
  other_cell_color = "lightgrey",
  other_point_size = 1,
  show_legend = TRUE
)
```


## Arguments

Argument      |Description
------------- |----------------
`sdimx`     |     x-axis dimension name (default = 'sdimx')
`sdimy`     |     y-axis dimension name (default = 'sdimy')
`cell_locations_metadata_selected`     |     annotated location from selected cells
`cell_locations_metadata_other`     |     annotated location from non-selected cells
`cell_color`     |     color for cells (see details)
`color_as_factor`     |     convert color column to factor
`cell_color_code`     |     named vector with colors
`cell_color_gradient`     |     vector with 3 colors for numeric data
`gradient_midpoint`     |     midpoint for color gradient
`gradient_limits`     |     vector with lower and upper limits
`select_cell_groups`     |     select subset of cells/clusters based on cell_color parameter
`select_cells`     |     select subset of cells based on cell IDs
`point_size`     |     size of point (cell)
`point_alpha`     |     transparancy of point
`point_border_col`     |     color of border around points
`point_border_stroke`     |     stroke size of border around points
`show_cluster_center`     |     plot center of selected clusters
`show_center_label`     |     plot label of selected clusters
`center_point_size`     |     size of center points
`label_size`     |     size of labels
`label_fontface`     |     font of labels
`show_other_cells`     |     display not selected cells
`other_cell_color`     |     color for not selected cells
`other_point_size`     |     point size for not selected cells
`show_legend`     |     show legend
`gobject`     |     giotto object


## Details

Description of parameters.


## Value

ggplot


