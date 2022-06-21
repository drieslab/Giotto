# `writeGiottoLargeImage`

writeGiottoLargeImage


## Description

Write full resolution image to file. Filetype extension should be
 included in `filename` argument. Be careful if write time and disk space
 needed if image is very large.


## Usage

```r
writeGiottoLargeImage(
  giottoLargeImage = NULL,
  gobject = NULL,
  largeImage_name = NULL,
  filename = NULL,
  dataType = NULL,
  max_intensity = NULL,
  overwrite = FALSE,
  verbose = TRUE
)
```


## Arguments

Argument      |Description
------------- |----------------
`giottoLargeImage`     |     giottoLargeImage object
`gobject`     |     giotto object
`largeImage_name`     |     name of giottoLargeImage
`filename`     |     file name and path to write the image to
`dataType`     |     (optional) values for `dataType` are "INT1U", "INT2U", "INT2S", "INT4U", "INT4S", "FLT4S", "FLT8S". The first three letters indicate whether the dataType is integer (whole numbers) of a real number (decimal numbers), the fourth character indicates the number of bytes used (allowing for large numbers and/or more precision), and the "S" or "U" indicate whether the values are signed (both negative and positive) or unsigned (positive values only).
`max_intensity`     |     (optional) image max intensity value from which `dataType`  can be automatically determined
`overwrite`     |     Overwrite if `filename` is already existing
`verbose`     |     be verbose


