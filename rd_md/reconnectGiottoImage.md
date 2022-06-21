# `reconnectGiottoImage`

Reconnect images with dead pointers


## Description

reconnect a gobject's dead image pointers using filepaths to
 the original source image files


## Usage

```r
reconnectGiottoImage(
  gobject,
  auto_reconnect = TRUE,
  image_name = NULL,
  largeImage_name = NULL,
  image_path = NULL,
  largeImage_path = NULL,
  verbose = TRUE
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object
`auto_reconnect`     |     automatically reconnect images if TRUE. manual if FALSE
`image_name`     |     names of images to reconnect
`largeImage_name`     |     name of large images to reconnect
`image_path`     |     named list of paths to images to reconnect to giottoImages
`largeImage_path`     |     named list of paths to images to reconnect to giottoLargeImages
`verbose`     |     be verbose


## Details

Inputs can either be given as both image name ( `image_name` / `largeImage_name` )
 and filepath ( `image_path` / `largeImage_path` ) args or as only a
 a named list through a filepath argument alone.
 If `auto_reconnect = TRUE` then no additional params need to be supplied.
 As long as giotto image objects were directly created using filepaths, those
 filepaths are stored within the image objects and will be referenced during
 reconnection. Issues will only arise if giotto image objects were created directly
 from the underlying image handling package objects ( magick or raster objects ) or
 if image files have been moved since the the giotto image object was generated.
 In such cases, use manual reconnection by setting `auto_reconnect = FALSE` .


## Value

a giotto object with updated image pointer


## Seealso

Other basic image functions:
 [`addGiottoImage`](#addgiottoimage) ,
 [`plotGiottoImage`](#plotgiottoimage) ,
 [`updateGiottoImage`](#updategiottoimage)


