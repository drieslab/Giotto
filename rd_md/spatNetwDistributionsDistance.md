# `spatNetwDistributionsDistance`

Plot spatial distance distribution


## Description

This function return histograms displaying the distance distribution for each spatial k-neighbor


## Usage

```r
spatNetwDistributionsDistance(
  gobject,
  spat_unit = NULL,
  spatial_network_name = "spatial_network",
  hist_bins = 30,
  test_distance_limit = NULL,
  ncol = 1,
  show_plot = NA,
  return_plot = NA,
  save_plot = NA,
  save_param = list(),
  default_save_name = "spatNetwDistributionsDistance"
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     Giotto object
`spat_unit`     |     spatial unit
`spatial_network_name`     |     name of spatial network
`hist_bins`     |     number of binds to use for the histogram
`test_distance_limit`     |     effect of different distance threshold on k-neighbors
`ncol`     |     number of columns to visualize the histograms in
`show_plot`     |     show plot
`return_plot`     |     return ggplot object
`save_plot`     |     directly save the plot [boolean]
`save_param`     |     list of saving parameters from [`all_plots_save_function`](#allplotssavefunction)
`default_save_name`     |     default save name for saving, alternatively change save_name in save_param


## Value

ggplot plot


