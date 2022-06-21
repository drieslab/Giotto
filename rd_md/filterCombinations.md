# `filterCombinations`

filterCombinations


## Description

Shows how many genes and cells are lost with combinations of thresholds.


## Usage

```r
filterCombinations(
  gobject,
  feat_type = NULL,
  spat_unit = NULL,
  expression_values = c("raw", "normalized", "scaled", "custom"),
  expression_thresholds = c(1, 2),
  feat_det_in_min_cells = c(5, 50),
  gene_det_in_min_cells = NULL,
  min_det_feats_per_cell = c(200, 400),
  min_det_genes_per_cell = NULL,
  scale_x_axis = "identity",
  x_axis_offset = 0,
  scale_y_axis = "identity",
  y_axis_offset = 0,
  show_plot = TRUE,
  return_plot = FALSE,
  save_plot = NA,
  save_param = list(),
  default_save_name = "filterCombinations"
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object
`feat_type`     |     feature type
`spat_unit`     |     spatial unit
`expression_values`     |     expression values to use
`expression_thresholds`     |     all thresholds to consider a gene expressed
`feat_det_in_min_cells`     |     minimum # of cells that need to express a feature
`gene_det_in_min_cells`     |     deprecated, use feat_det_in_min_cells
`min_det_feats_per_cell`     |     minimum # of features that need to be detected in a cell
`min_det_genes_per_cell`     |     deprecated, use min_det_feats_per_cell
`scale_x_axis`     |     ggplot transformation for x-axis (e.g. log2)
`x_axis_offset`     |     x-axis offset to be used together with the scaling transformation
`scale_y_axis`     |     ggplot transformation for y-axis (e.g. log2)
`y_axis_offset`     |     y-axis offset to be used together with the scaling transformation
`show_plot`     |     show plot
`return_plot`     |     return only ggplot object
`save_plot`     |     directly save the plot [boolean]
`save_param`     |     list of saving parameters from [`all_plots_save_function`](#allplotssavefunction)
`default_save_name`     |     default save name for saving, don't change, change save_name in save_param


## Details

Creates a scatterplot that visualizes the number of genes and cells that are
 lost with a specific combination of a gene and cell threshold given an arbitrary cutoff
 to call a gene expressed. This function can be used to make an informed decision at the
 filtering step with filterGiotto.


## Value

list of data.table and ggplot object


