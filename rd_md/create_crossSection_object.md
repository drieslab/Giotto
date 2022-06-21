# `create_crossSection_object`

create_crossSection_object


## Description

create a crossSection object


## Usage

```r
create_crossSection_object(
  name = NULL,
  method = NULL,
  thickness_unit = NULL,
  slice_thickness = NULL,
  cell_distance_estimate_method = NULL,
  extend_ratio = NULL,
  plane_equation = NULL,
  mesh_grid_n = NULL,
  mesh_obj = NULL,
  cell_subset = NULL,
  cell_subset_spatial_locations = NULL,
  cell_subset_projection_locations = NULL,
  cell_subset_projection_PCA = NULL,
  cell_subset_projection_coords = NULL
)
```


## Arguments

Argument      |Description
------------- |----------------
`name`     |     name of cress section object. (default = cross_sectino)
`method`     |     method to define the cross section plane.
`thickness_unit`     |     unit of the virtual section thickness. If "cell", average size of the observed cells is used as length unit. If "natural", the unit of cell location coordinates is used.(default = cell)
`slice_thickness`     |     thickness of slice
`cell_distance_estimate_method`     |     method to estimate average distance between neighobring cells. (default = mean)
`extend_ratio`     |     deciding the span of the cross section meshgrid, as a ratio of extension compared to the borders of the vitural tissue section. (default = 0.2)
`plane_equation`     |     a numerical vector of length 4, in the form of c(A,B,C,D), which defines plane Ax+By+Cz=D.
`mesh_grid_n`     |     numer of meshgrid lines to generate along both directions for the cross section plane.
`mesh_obj`     |     object that stores the cross section meshgrid information.
`cell_subset`     |     cells selected by the cross section
`cell_subset_spatial_locations`     |     locations of cells selected by the cross section
`cell_subset_projection_locations`     |     3D projection coordinates of selected cells onto the cross section plane
`cell_subset_projection_PCA`     |     pca of projection coordinates
`cell_subset_projection_coords`     |     2D PCA coordinates of selected cells in the cross section plane


