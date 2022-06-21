# `addGiottoLargeImage`

addGiottoLargeImage


## Description

Adds giotto image objects to your giotto object


## Usage

```r
addGiottoLargeImage(
  gobject = NULL,
  largeImages = NULL,
  spat_loc_name = NULL,
  scale_factor = NULL,
  negative_y = TRUE,
  verbose = TRUE
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object
`largeImages`     |     list of giottoLargeImage objects
`spat_loc_name`     |     provide spatial location slot in Giotto to align images. (optional)
`scale_factor`     |     provide scale of image pixel dimensions relative to spatial coordinates.
`negative_y`     |     map image to negative y spatial values if TRUE during automatic alignment. Meaning that origin is in upper left instead of lower left.
`verbose`     |     be verbose


## Value

an updated Giotto object with access to the list of images


