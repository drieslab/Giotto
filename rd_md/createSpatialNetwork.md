# `createSpatialNetwork`

createSpatialNetwork


## Description

Create a spatial network based on cell centroid physical distances.


## Usage

```r
createSpatialNetwork(
  gobject,
  name = NULL,
  spat_unit = NULL,
  feat_type = NULL,
  spat_loc_name = NULL,
  dimensions = "all",
  method = c("Delaunay", "kNN"),
  delaunay_method = c("deldir", "delaunayn_geometry", "RTriangle"),
  maximum_distance_delaunay = "auto",
  options = "Pp",
  Y = TRUE,
  j = TRUE,
  S = 0,
  minimum_k = 0,
  knn_method = "dbscan",
  k = 4,
  maximum_distance_knn = NULL,
  verbose = F,
  return_gobject = TRUE,
  ...
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object
`name`     |     name for spatial network (default = 'spatial_network')
`spat_unit`     |     spatial unit
`feat_type`     |     feature type
`spat_loc_name`     |     name of spatial locations to use
`dimensions`     |     which spatial dimensions to use (default = all)
`method`     |     which method to use to create a spatial network. (default = Delaunay)
`delaunay_method`     |     Delaunay method to use
`maximum_distance_delaunay`     |     distance cuttof for nearest neighbors to consider for Delaunay network
`options`     |     (geometry) String containing extra control options for the underlying Qhull command; see the Qhull documentation (../doc/qhull/html/qdelaun.html) for the available options. (default = 'Pp', do not report precision problems)
`Y`     |     (RTriangle) If TRUE prohibits the insertion of Steiner points on the mesh boundary.
`j`     |     (RTriangle) If TRUE jettisons vertices that are not part of the final triangulation from the output.
`S`     |     (RTriangle) Specifies the maximum number of added Steiner points.
`minimum_k`     |     minimum nearest neigbhours if maximum_distance != NULL
`knn_method`     |     method to create kNN network
`k`     |     number of nearest neighbors based on physical distance
`maximum_distance_knn`     |     distance cuttof for nearest neighbors to consider for kNN network
`verbose`     |     verbose
`return_gobject`     |     boolean: return giotto object (default = TRUE)
`list()`     |     Additional parameters for the selected function


## Details

Creates a spatial network connecting single-cells based on their physical distance to each other.
 For Delaunay method, neighbors will be decided by delaunay triangulation and a maximum distance criteria. For kNN method, number of neighbors can be determined by k, or maximum distance from each cell with or without
 setting a minimum k for each cell.
 
 dimensions:  default = 'all' which takes all possible dimensions.
 Alternatively you can provide a character vector that specififies the spatial dimensions to use, e.g. c("sdimx', "sdimy")
 or a numerical vector, e.g. 2:3


## Value

giotto object with updated spatial network slot


