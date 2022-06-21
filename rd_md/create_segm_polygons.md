# `create_segm_polygons`

Create segmentation polygons


## Description

creates giotto polygons from segmentation mask data


## Usage

```r
create_segm_polygons(
  maskfile,
  name = "cell",
  poly_IDs = NULL,
  flip_vertical = TRUE,
  shift_vertical_step = TRUE,
  flip_horizontal = TRUE,
  shift_horizontal_step = TRUE,
  remove_background_polygon = FALSE
)
```


## Value

giotto polygon


