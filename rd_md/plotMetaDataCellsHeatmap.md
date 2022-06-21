# `plotMetaDataCellsHeatmap`

plotMetaDataCellsHeatmap


## Description

Creates heatmap for numeric cell metadata within aggregated clusters.


## Usage

```r
plotMetaDataCellsHeatmap(
  gobject,
  spat_unit = NULL,
  feat_type = NULL,
  metadata_cols = NULL,
  spat_enr_names = NULL,
  value_cols = NULL,
  first_meta_col = NULL,
  second_meta_col = NULL,
  show_values = c("zscores", "original", "zscores_rescaled"),
  custom_cluster_order = NULL,
  clus_cor_method = "pearson",
  clus_cluster_method = "complete",
  custom_values_order = NULL,
  values_cor_method = "pearson",
  values_cluster_method = "complete",
  midpoint = 0,
  x_text_size = 8,
  x_text_angle = 45,
  y_text_size = 8,
  strip_text_size = 8,
  show_plot = NA,
  return_plot = NA,
  save_plot = NA,
  save_param = list(),
  default_save_name = "plotMetaDataCellsHeatmap"
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object
`spat_unit`     |     spatial unit
`feat_type`     |     feature type
`metadata_cols`     |     annotation columns found in pDataDT(gobject)
`spat_enr_names`     |     spatial enrichment results to include
`value_cols`     |     value columns to use
`first_meta_col`     |     if more than 1 metadata column, select the x-axis factor
`second_meta_col`     |     if more than 1 metadata column, select the facetting factor
`show_values`     |     which values to show on heatmap
`custom_cluster_order`     |     custom cluster order (default = NULL)
`clus_cor_method`     |     correlation method for clusters
`clus_cluster_method`     |     hierarchical cluster method for the clusters
`custom_values_order`     |     custom values order (default = NULL)
`values_cor_method`     |     correlation method for values
`values_cluster_method`     |     hierarchical cluster method for the values
`midpoint`     |     midpoint of show_values
`x_text_size`     |     size of x-axis text
`x_text_angle`     |     angle of x-axis text
`y_text_size`     |     size of y-axis text
`strip_text_size`     |     size of strip text
`show_plot`     |     show plot
`return_plot`     |     return ggplot object
`save_plot`     |     directly save the plot [boolean]
`save_param`     |     list of saving parameters, see [`showSaveParameters`](#showsaveparameters)
`default_save_name`     |     default save name for saving, don't change, change save_name in save_param


## Details

Creates heatmap for the average values of selected value columns in the different annotation groups.


## Value

ggplot or data.table


## Seealso

[`plotMetaDataHeatmap`](#plotmetadataheatmap) for gene expression instead of numeric cell annotation data.


