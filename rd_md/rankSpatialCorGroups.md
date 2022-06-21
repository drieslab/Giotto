# `rankSpatialCorGroups`

rankSpatialCorGroups


## Description

Rank spatial correlated clusters according to correlation structure


## Usage

```r
rankSpatialCorGroups(
  gobject,
  spatCorObject,
  use_clus_name = NULL,
  show_plot = NA,
  return_plot = FALSE,
  save_plot = NA,
  save_param = list(),
  default_save_name = "rankSpatialCorGroups"
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object
`spatCorObject`     |     spatial correlation object
`use_clus_name`     |     name of clusters to visualize (from clusterSpatialCorGenes())
`show_plot`     |     show plot
`return_plot`     |     return ggplot object
`save_plot`     |     directly save the plot [boolean]
`save_param`     |     list of saving parameters, see [`showSaveParameters`](#showsaveparameters)
`default_save_name`     |     default save name for saving, don't change, change save_name in save_param


## Value

data.table with positive (within group) and negative (outside group) scores


