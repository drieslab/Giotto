# `plot_giottoLargeImage`

plot_giottoLargeImage


## Description

Plot a downsampled version of giottoLargeImage. Cropping can increase plot resolution of region of interest.


## Usage

```r
plot_giottoLargeImage(
  gobject = NULL,
  largeImage_name = NULL,
  giottoLargeImage = NULL,
  crop_extent = NULL,
  xmax_crop = NULL,
  xmin_crop = NULL,
  ymax_crop = NULL,
  ymin_crop = NULL,
  max_intensity = NULL,
  asRGB = FALSE
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object
`largeImage_name`     |     name of giottoLargeImage
`giottoLargeImage`     |     giottoLargeImage object
`crop_extent`     |     (optional) extent object to focus on specific region of image
`xmax_crop, xmin_crop, ymax_crop, ymin_crop`     |     (optional) crop min/max x and y bounds
`max_intensity`     |     (optional) value to treat as maximum intensity in color scale
`asRGB`     |     (optional) [boolean] force RGB plotting if not automatically detected


## Value

plot


