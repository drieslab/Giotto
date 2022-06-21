# `createSpatialKNNnetwork`

createSpatialKNNnetwork


## Description

Create a spatial knn network.


## Usage

```r
createSpatialKNNnetwork(
  gobject,
  method = "dbscan",
  spat_unit = NULL,
  feat_type = NULL,
  spat_loc_name = NULL,
  dimensions = "all",
  name = "knn_network",
  k = 4,
  maximum_distance = NULL,
  minimum_k = 0,
  verbose = F,
  return_gobject = TRUE,
  ...
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object
`method`     |     method to create kNN network
`spat_unit`     |     spatial unit
`feat_type`     |     feature type
`spat_loc_name`     |     name of spatial locations
`dimensions`     |     which spatial dimensions to use (default = all)
`name`     |     name for spatial network (default = 'spatial_network')
`k`     |     number of nearest neighbors based on physical distance
`maximum_distance`     |     distance cuttof for nearest neighbors to consider for kNN network
`minimum_k`     |     minimum nearest neigbhours if maximum_distance != NULL
`verbose`     |     verbose
`return_gobject`     |     boolean: return giotto object (default = TRUE)
`list()`     |     additional arguments to the selected method function


## Value

giotto object with updated spatial network slot
 
 dimensions:  default = 'all' which takes all possible dimensions.
 Alternatively you can provide a character vector that specififies the spatial dimensions to use, e.g. c("sdimx', "sdimy")
 or a numerical vector, e.g. 2:3
 
 maximum_distance:  to create a network based on maximum distance only, you also need to set k to a very high value, e.g. k = 100


