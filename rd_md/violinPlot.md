# `violinPlot`

violinPlot


## Description

Creates violinplot for selected clusters


## Usage

```r
violinPlot(
  gobject,
  spat_unit = NULL,
  feat_type = NULL,
  expression_values = c("normalized", "scaled", "custom"),
  feats = NULL,
  genes = NULL,
  cluster_column,
  cluster_custom_order = NULL,
  color_violin = c("feats", "cluster"),
  cluster_color_code = NULL,
  strip_position = c("top", "right", "left", "bottom"),
  strip_text = 7,
  axis_text_x_size = 10,
  axis_text_y_size = 6,
  show_plot = NA,
  return_plot = NA,
  save_plot = NA,
  save_param = list(),
  default_save_name = "violinPlot"
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object
`spat_unit`     |     spatial unit
`feat_type`     |     feature type
`expression_values`     |     expression values to use
`feats`     |     features to plot
`genes`     |     deprecated, use feats argument
`cluster_column`     |     name of column to use for clusters
`cluster_custom_order`     |     custom order of clusters
`color_violin`     |     color violin according to genes or clusters
`cluster_color_code`     |     color code for clusters
`strip_position`     |     position of gene labels
`strip_text`     |     size of strip text
`axis_text_x_size`     |     size of x-axis text
`axis_text_y_size`     |     size of y-axis text
`show_plot`     |     show plot
`return_plot`     |     return ggplot object
`save_plot`     |     directly save the plot [boolean]
`save_param`     |     list of saving parameters, see [`showSaveParameters`](#showsaveparameters)
`default_save_name`     |     default save name for saving, don't change, change save_name in save_param


## Value

ggplot


