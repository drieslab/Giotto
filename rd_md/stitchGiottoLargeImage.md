# `stitchGiottoLargeImage`

Stitch multiple giottoLargeImage objects into a single giottoLargeImage object


## Description

Function to stitch together multiple field of view (FOV) images into a
 single final image. Images are loaded into Giotto as `giottoLargeImage` and
 stitched based on a set of FOV positions into a single final `giottoLargeImage` .


## Usage

```r
stitchGiottoLargeImage(
  largeImage_list = NULL,
  gobject_list = NULL,
  largeImage_nameList = NULL,
  FOV_positions = NULL,
  FOV_xcol = NULL,
  FOV_ycol = NULL,
  FOV_inverty = FALSE,
  method = c("mosaic", "merge"),
  round_positions = FALSE,
  filename = NULL,
  dataType = NULL,
  fileType = NULL,
  dryRun = TRUE,
  overwrite = FALSE,
  verbose = TRUE
)
```


## Arguments

Argument      |Description
------------- |----------------
`largeImage_list`     |     list of `giottoLargeImage` objects
`gobject_list`     |     list of `gobjects` containing `giottoLargeImages`
`largeImage_nameList`     |     list of names of `giottoLargeImages` within `gobjects`
`FOV_positions`     |     dataframe of FOV positions. Values (if any) are directly added to current image mapping
`FOV_xcol`     |     column name for FOV position x values
`FOV_ycol`     |     column name for FOV position y values
`FOV_inverty`     |     make FOV y position values negative
`method`     |     method of stitching images ( mosaic : average overlapping area values, merge :values get priority by order of images given)
`round_positions`     |     [boolean] round image positions. May be necessary to run.
`filename`     |     file name to write the stitched image to. Defaults to `"save_dir/stitch.tif"`  if `save_dir` param is found in the first `gobject` 's Giotto instructions
`dataType`     |     (optional) values for `dataType` are "INT1U", "INT2U", "INT2S", "INT4U", "INT4S", "FLT4S", "FLT8S". The first three letters indicate whether the `dataType` is integer (whole numbers) of a real number (decimal numbers), the fourth character indicates the number of bytes used (allowing for large numbers and/or more precision), and the "S" or "U" indicate whether the values are signed (both negative and positive) or unsigned (positive values only).
`fileType`     |     (optional) image format (e.g. .tif) If not given, defaults to format given in the filename
`dryRun`     |     [boolean] plot placeholder bounding rectangles where FOV images will be stitched without actually proceeding with the full image stitching and saving process.
`overwrite`     |     [boolean] overwrite if filename to save image as already exists. Defaults to TRUE
`verbose`     |     [boolean] be verbose


## Details

This function is time consuming. Setting a save location through the
  `filename` parameter is also highly recommended as file size will likely be large.
 This function creates a single stitched image from multiple FOV tiles and saves that
 image to disk as it works. When finished, the pointer to that new image is loaded in
 as a `giottoLargeImage` object. list() 
  list("Note:") Dry runs are on by default and `dryRun` param must be set to FALSE
 to proceed with the final stitching operation.


## Value

`largeGiottoImage` object with pointer to stitched image


