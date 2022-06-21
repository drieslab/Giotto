# `createSpatialFeaturesKNNnetwork`

Create kNN spatial feature network


## Description

Calculates the centroid locations for the giotto polygons


## Usage

```r
createSpatialFeaturesKNNnetwork(
  gobject,
  method = c("dbscan"),
  feat_type = NULL,
  name = "knn_feats_network",
  k = 4,
  maximum_distance = NULL,
  minimum_k = 0,
  add_feat_ids = FALSE,
  verbose = TRUE,
  return_gobject = TRUE,
  toplevel_params = 2,
  ...
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object
`method`     |     kNN algorithm method
`feat_type`     |     feature type to build feature network
`name`     |     name of network
`k`     |     number of neighbors
`maximum_distance`     |     maximum distance bewteen features
`minimum_k`     |     minimum number of neighbors to find
`add_feat_ids`     |     add feature id names (default = FALSE, increases object size)
`verbose`     |     be verbose
`return_gobject`     |     return giotto object (default: TRUE)
`toplevel_params`     |     toplevel value to pass when updating giotto params
`...`     |     additional parameters to pass to [`kNN`](#knn)


## Value

If `return_gobject = TRUE` a giotto object containing the network
 will be returned. If `return_gobject = FALSE` the network will be returned
 as a datatable.


