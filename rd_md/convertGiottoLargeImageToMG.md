# `convertGiottoLargeImageToMG`

convertGiottoLargeImageToMG


## Description

convert a giottoLargeImage by downsampling into a normal magick based giottoImage


## Usage

```r
convertGiottoLargeImageToMG(
  gobject = NULL,
  largeImage_name = NULL,
  giottoLargeImage = NULL,
  mg_name = NULL,
  spat_unit = NULL,
  spat_loc_name = NULL,
  crop_extent = NULL,
  xmax_crop = NULL,
  xmin_crop = NULL,
  ymax_crop = NULL,
  ymin_crop = NULL,
  resample_size = 5e+05,
  max_intensity = NULL,
  return_gobject = TRUE,
  verbose = TRUE
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     gobject containing giottoLargeImage
`largeImage_name`     |     name of giottoLargeImage
`giottoLargeImage`     |     alternative input param using giottoLargeImage object instead of through `gobject` and `largeImage_name` params
`mg_name`     |     name to assign converted magick image based giottoImage. Defaults to name of giottoLargeImage
`spat_unit`     |     spatial unit
`spat_loc_name`     |     gobject spatial location name to map giottoImage to (optional)
`crop_extent`     |     extent object to focus on specific region of image
`xmax_crop`     |     assign crop boundary
`xmin_crop`     |     assign crop boundary
`ymax_crop`     |     assign crop boundary
`ymin_crop`     |     assign crop boundary
`resample_size`     |     maximum number of pixels to use when resampling
`max_intensity`     |     value to treat as maximum intensity in color scale
`return_gobject`     |     return as giotto object
`verbose`     |     be verbose


## Value

a giotto object if `return_gobject = TRUE` or an updated giotto
 image object if `return_gobject = FALSE`


