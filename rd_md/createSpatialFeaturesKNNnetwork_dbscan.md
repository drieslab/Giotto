# `createSpatialFeaturesKNNnetwork_dbscan`

Create kNN spatial feature network using dbscan


## Description

to create a feature kNN spatial network using dbscan


## Usage

```r
createSpatialFeaturesKNNnetwork_dbscan(
  gobject,
  feat_type = NULL,
  name = "knn_feats_network",
  k = 4,
  maximum_distance = NULL,
  minimum_k = 0,
  add_feat_ids = FALSE,
  verbose = TRUE,
  ...
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object
`feat_type`     |     feature type
`name`     |     name to assign generated feature network
`k`     |     number of neighbors for kNN to find
`maximum_distance`     |     network maximum distance allowed
`minimum_k`     |     minimum neighbors allowed
`add_feat_ids`     |     whether to add feature information [boolean]
`verbose`     |     be verbose
`...`     |     additional parameters to pass to [`kNN`](#knn)


