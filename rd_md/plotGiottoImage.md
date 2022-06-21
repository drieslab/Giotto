# `plotGiottoImage`

Plot a giotto image object


## Description

Display a giotto image in the viewer panel. Image object to plot
 can be specified by providing the giotto object containing the image ( `gobject` ),
 the image object name ( `image_name` ), and the image object type ( `image_type` ).
 Alternatively, image objects can be directly plotted through their respective
 associated params.


## Usage

```r
plotGiottoImage(
  gobject = NULL,
  image_name = NULL,
  image_type = NULL,
  giottoImage = NULL,
  giottoLargeImage = NULL,
  largeImage_crop_params_list = NULL,
  largeImage_max_intensity = NULL,
  ...
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     gobject containing giotto image object
`image_name`     |     name of giotto image object
`image_type`     |     type of giotto image object to plot
`giottoImage`     |     giottoImage object to plot directly
`giottoLargeImage`     |     giottoLargeImage object to plot directly
`largeImage_crop_params_list`     |     (optional) named list of params for focusing on a specified region of a giottoLargeImage.
`largeImage_max_intensity`     |     (optional) assign override value to treat as maximum intensity in color scale when plotting giottoLargeImage
`...`     |     additional params to pass to image object specific plotting functions


## Seealso

Other basic image functions:
 [`addGiottoImage`](#addgiottoimage) ,
 [`reconnectGiottoImage`](#reconnectgiottoimage) ,
 [`updateGiottoImage`](#updategiottoimage)


