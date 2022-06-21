# `createGiottoLargeImageList`

createGiottoLargeImageList


## Description

Creates a list of large giotto images that can be added to a Giotto object. Generates deep copy of SpatRaster


## Usage

```r
createGiottoLargeImageList(
  raster_objects,
  names = "image",
  negative_y = TRUE,
  extent = NULL,
  use_rast_ext = FALSE,
  image_transformations = NULL,
  xmax_bound = NULL,
  xmin_bound = NULL,
  ymax_bound = NULL,
  ymin_bound = NULL,
  scale_factor = 1,
  verbose = TRUE
)
```


## Arguments

Argument      |Description
------------- |----------------
`raster_objects`     |     vector of image paths or terra SpatRaster image objects
`names`     |     vector of names for the images
`negative_y`     |     Map image to negative y spatial values if TRUE. Meaning that origin is in upper left instead of lower left.
`extent`     |     SpatExtent object to assign spatial extent. Takes priority unless use_rast_ext is TRUE.
`use_rast_ext`     |     Use extent from input raster object
`image_transformations`     |     vector of sequential image transformations - under construction
`xmax_bound, xmin_bound, ymax_bound, ymin_bound`     |     assign min and max x and y values for image spatial placement
`scale_factor`     |     scaling of image dimensions relative to spatial coordinates
`verbose`     |     be verbose


## Details

See [`createGiottoLargeImage`](#creategiottolargeimage)


## Value

a list with giottoLargeImage objects


