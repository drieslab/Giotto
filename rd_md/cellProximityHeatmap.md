# `cellProximityHeatmap`

cellProximityHeatmap


## Description

Create heatmap from cell-cell proximity scores


## Usage

```r
cellProximityHeatmap(
  gobject,
  CPscore,
  scale = T,
  order_cell_types = T,
  color_breaks = NULL,
  color_names = NULL,
  show_plot = NA,
  return_plot = NA,
  save_plot = NA,
  save_param = list(),
  default_save_name = "cellProximityHeatmap"
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object
`CPscore`     |     CPscore, output from cellProximityEnrichment()
`scale`     |     scale cell-cell proximity interaction scores
`order_cell_types`     |     order cell types based on enrichment correlation
`color_breaks`     |     numerical vector of length 3 to represent min, mean and maximum
`color_names`     |     character color vector of length 3
`show_plot`     |     show plot
`return_plot`     |     return ggplot object
`save_plot`     |     directly save the plot [boolean]
`save_param`     |     list of saving parameters from [`all_plots_save_function`](#allplotssavefunction)
`default_save_name`     |     default save name for saving, don't change, change save_name in save_param


## Details

This function creates a heatmap that shows the  spatial proximity
 enrichment or depletion of cell type pairs.


## Value

ggplot heatmap


