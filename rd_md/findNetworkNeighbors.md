# `findNetworkNeighbors`

Find network neighbors


## Description

Find the spatial neighbors for a selected group of cells within the selected spatial network.


## Usage

```r
findNetworkNeighbors(
  gobject,
  spat_unit = NULL,
  spatial_network_name,
  source_cell_ids = NULL,
  name = "nb_cells"
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     Giotto object
`spat_unit`     |     spatial unit
`spatial_network_name`     |     name of spatial network
`source_cell_ids`     |     cell ids for which you want to know the spatial neighbors
`name`     |     name of the results


## Value

data.table


