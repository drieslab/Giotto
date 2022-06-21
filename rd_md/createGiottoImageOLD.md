# `createGiottoImageOLD`

createGiottoImageOLD


## Description

Creates a giotto image that can be added to a Giotto object and/or
 used to add an image to the spatial plotting functions. Deprecated. See [`createGiottoImage`](#creategiottoimage)


## Usage

```r
createGiottoImageOLD(
  gobject = NULL,
  spatial_locs = NULL,
  mg_object,
  name = "image",
  xmax_adj = 0,
  xmin_adj = 0,
  ymax_adj = 0,
  ymin_adj = 0
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object
`spatial_locs`     |     spatial locations (alternative if `gobject = NULL` )
`mg_object`     |     magick image object
`name`     |     name for the image
`xmax_adj`     |     adjustment of the maximum x-value to align the image
`xmin_adj`     |     adjustment of the minimum x-value to align the image
`ymax_adj`     |     adjustment of the maximum y-value to align the image
`ymin_adj`     |     adjustment of the minimum y-value to align the image


## Value

a giotto image object


