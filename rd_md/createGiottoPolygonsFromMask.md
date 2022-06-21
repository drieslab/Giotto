# `createGiottoPolygonsFromMask`

Create giotto polygons from mask file


## Description

Creates Giotto polygon object from a mask file (e.g. segmentation results)


## Usage

```r
createGiottoPolygonsFromMask(
  maskfile,
  mask_method = c("guess", "single", "multiple"),
  name = "cell",
  remove_background_polygon = FALSE,
  background_algo = c("range"),
  fill_holes = TRUE,
  poly_IDs = NULL,
  flip_vertical = TRUE,
  shift_vertical_step = TRUE,
  flip_horizontal = TRUE,
  shift_horizontal_step = TRUE,
  calc_centroids = FALSE,
  fix_multipart = TRUE,
  remove_unvalid_polygons = TRUE
)
```


## Arguments

Argument      |Description
------------- |----------------
`maskfile`     |     path to mask file
`mask_method`     |     how the mask file defines individual segmentation annotations
`name`     |     name for polygons
`remove_background_polygon`     |     try to remove background polygon (default: FALSE)
`background_algo`     |     algorithm to remove background polygon
`fill_holes`     |     fill holes within created polygons
`poly_IDs`     |     unique nanes for each polygon in the mask file
`flip_vertical`     |     flip mask figure in a vertical manner
`shift_vertical_step`     |     shift vertical (boolean or numerical)
`flip_horizontal`     |     flip mask figure in a horizontal manner
`shift_horizontal_step`     |     shift horizontal (boolean or numerical)
`calc_centroids`     |     calculate centroids for polygons
`fix_multipart`     |     try to split polygons with multiple parts (default: TRUE)
`remove_unvalid_polygons`     |     remove unvalid polygons (default: TRUE)


## Value

a giotto polygon object


