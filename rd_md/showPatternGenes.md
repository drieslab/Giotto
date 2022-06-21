# `showPatternGenes`

showPatternGenes


## Description

show genes correlated with spatial patterns


## Usage

```r
showPatternGenes(
  gobject,
  spatPatObj,
  dimension = 1,
  top_pos_genes = 5,
  top_neg_genes = 5,
  point_size = 1,
  return_DT = FALSE,
  show_plot = NA,
  return_plot = NA,
  save_plot = NA,
  save_param = list(),
  default_save_name = "showPatternGenes"
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object
`spatPatObj`     |     Output from detectSpatialPatterns
`dimension`     |     dimension to plot genes for.
`top_pos_genes`     |     Top positively correlated genes.
`top_neg_genes`     |     Top negatively correlated genes.
`point_size`     |     size of points
`return_DT`     |     if TRUE, it will return the data.table used to generate the plots
`show_plot`     |     show plot
`return_plot`     |     return ggplot object
`save_plot`     |     directly save the plot [boolean]
`save_param`     |     list of saving parameters, see [`showSaveParameters`](#showsaveparameters)
`default_save_name`     |     default save name for saving, don't change, change save_name in save_param


## Value

ggplot


