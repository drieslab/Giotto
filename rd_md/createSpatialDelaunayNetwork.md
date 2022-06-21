# `createSpatialDelaunayNetwork`

Create a spatial Delaunay network


## Description

Create a spatial Delaunay network based on cell centroid physical distances.


## Usage

```r
createSpatialDelaunayNetwork(
  gobject,
  name = "Delaunay_network",
  spat_unit = NULL,
  feat_type = NULL,
  spat_loc_name = NULL,
  method = c("deldir", "delaunayn_geometry", "RTriangle"),
  dimensions = "all",
  maximum_distance = "auto",
  minimum_k = 0,
  options = "Pp",
  Y = TRUE,
  j = TRUE,
  S = 0,
  verbose = T,
  return_gobject = TRUE,
  ...
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object
`name`     |     name for spatial network (default = 'delaunay_network')
`spat_unit`     |     spatial unit
`feat_type`     |     feature type
`spat_loc_name`     |     name of spatial locations
`method`     |     package to use to create a Delaunay network
`dimensions`     |     which spatial dimensions to use. Use "sdimx" (spatial dimension x), "sdimy", "sdimz" respectively to refer to X (or the 1st), Y (or the 2nd) and Z(or the 3rd) dimension, see details. (default = all)
`maximum_distance`     |     distance cuttof for Delaunay neighbors to consider. If "auto", "upper wisker" value of the distance vector between neighbors is used; see the boxplot graphics documentation for more details.(default = "auto")
`minimum_k`     |     minimum number of neigbhours if maximum_distance != NULL
`options`     |     (geometry) String containing extra control options for the underlying Qhull command; see the Qhull documentation (../doc/qhull/html/qdelaun.html) for the available options. (default = 'Pp', do not report precision problems)
`Y`     |     (RTriangle) If TRUE prohibits the insertion of Steiner points on the mesh boundary.
`j`     |     (RTriangle) If TRUE jettisons vertices that are not part of the final triangulation from the output.
`S`     |     (RTriangle) Specifies the maximum number of added Steiner points.
`verbose`     |     verbose
`return_gobject`     |     boolean: return giotto object (default = TRUE)
`list()`     |     Other additional parameters


## Details

Creates a spatial Delaunay network as explained in [`delaunayn`](#delaunayn) (default), [`deldir`](#deldir) , or [`triangulate`](#triangulate) .


## Value

giotto object with updated spatial network slot


