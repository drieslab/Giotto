# `find_terra_writeRaster_dataType`

find_terra_writeRaster_dataType


## Description

find likely compatible datatype for given image characteristics.
 Values given in arguments take priority over those found from giottoLargeImage
 metadata


## Usage

```r
find_terra_writeRaster_dataType(
  giottoLargeImage = NULL,
  quick_INTS_maxval = NULL,
  max_intensity = NULL,
  min_intensity = NULL,
  is_int = NULL,
  signed = NULL,
  bitDepth = NULL,
  verbose = TRUE
)
```


## Arguments

Argument      |Description
------------- |----------------
`giottoLargeImage`     |     giottoLargeImage object to determine max_intensity, min_intensity, is_int settings from
`max_intensity, min_intensity`     |     value given as image maximum/minimum intensity
`is_int`     |     if image is integer (TRUE) or floating point (FALSE)
`signed`     |     if image is signed (TRUE) or unsigned (TRUE)
`bitDepth`     |     image bitDepth
`verbose`     |     be verbose
`quick_INTU_maxval`     |     Treat as maximum intensity to find compatible unsigned integer settings


## Value

datatype for terra writeRaster function


