# `showClusterHeatmap`

showClusterHeatmap


## Description

Creates heatmap based on identified clusters


## Usage

```r
showClusterHeatmap(
  gobject,
  spat_unit = NULL,
  feat_type = NULL,
  expression_values = c("normalized", "scaled", "custom"),
  feats = "all",
  genes = NULL,
  cluster_column,
  cor = c("pearson", "spearman"),
  distance = "ward.D",
  show_plot = NA,
  return_plot = NA,
  save_plot = NA,
  save_param = list(),
  default_save_name = "showClusterHeatmap",
  ...
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object
`spat_unit`     |     spatial unit
`feat_type`     |     feature type
`expression_values`     |     expression values to use
`feats`     |     vector of features to use, default to 'all'
`genes`     |     deprecated. Replaced by `feats` param
`cluster_column`     |     name of column to use for clusters
`cor`     |     correlation score to calculate distance
`distance`     |     distance method to use for hierarchical clustering
`show_plot`     |     show plot
`return_plot`     |     return ggplot object
`save_plot`     |     directly save the plot [boolean]
`save_param`     |     list of saving parameters, see [`showSaveParameters`](#showsaveparameters)
`default_save_name`     |     default save name for saving, don't change, change save_name in save_param
`...`     |     additional parameters for the Heatmap function from ComplexHeatmap


## Details

Correlation heatmap of selected clusters.


## Value

ggplot


