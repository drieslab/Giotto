# `plot_point_layer_ggplot`

plot_point_layer_ggplot


## Description

Visualize cells in point layer according to dimension reduction coordinates


## Usage

```r
plot_point_layer_ggplot(
  ggobject,
  annotated_DT_selected,
  annotated_DT_other,
  cell_color = NULL,
  color_as_factor = T,
  cell_color_code = NULL,
  cell_color_gradient = c("blue", "white", "red"),
  gradient_midpoint = 0,
  gradient_limits = NULL,
  select_cell_groups = NULL,
  select_cells = NULL,
  point_size = 1,
  point_alpha = 1,
  point_border_col = "black",
  point_border_stroke = 0.1,
  show_cluster_center = F,
  show_center_label = T,
  center_point_size = 4,
  center_point_border_col = "black",
  center_point_border_stroke = 0.1,
  label_size = 4,
  label_fontface = "bold",
  edge_alpha = NULL,
  show_other_cells = T,
  other_cell_color = "lightgrey",
  other_point_size = 0.5,
  show_legend = T
)
```


## Arguments

Argument      |Description
------------- |----------------
`annotated_DT_selected`     |     annotated data.table of selected cells
`annotated_DT_other`     |     annotated data.table of not selected cells
`cell_color`     |     color for cells (see details)
`color_as_factor`     |     convert color column to factor
`cell_color_code`     |     named vector with colors
`cell_color_gradient`     |     vector with 3 colors for numeric data
`gradient_midpoint`     |     midpoint for color gradient
`gradient_limits`     |     vector with lower and upper limits
`select_cell_groups`     |     select subset of cells/clusters based on cell_color parameter
`select_cells`     |     select subset of cells based on cell IDs
`point_size`     |     size of point (cell)
`point_alpha`     |     transparency of point
`point_border_col`     |     color of border around points
`point_border_stroke`     |     stroke size of border around points
`show_cluster_center`     |     plot center of selected clusters
`show_center_label`     |     plot label of selected clusters
`center_point_size`     |     size of center points
`label_size`     |     size of labels
`label_fontface`     |     font of labels
`edge_alpha`     |     column to use for alpha of the edges
`show_other_cells`     |     display not selected cells
`other_cell_color`     |     color of not selected cells
`other_point_size`     |     size of not selected cells
`show_legend`     |     show legend
`gobject`     |     giotto object


## Details

Description of parameters.


## Value

ggplot


