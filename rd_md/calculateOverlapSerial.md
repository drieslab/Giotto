# `calculateOverlapSerial`

calculateOverlapSerial


## Description

calculate overlap between cellular structures (polygons) and features (points)


## Usage

```r
calculateOverlapSerial(
  gobject,
  name_overlap = NULL,
  spatial_info = "cell",
  feat_info = "rna",
  poly_ID_names = "all",
  polygon_group_size = 500,
  return_gobject = TRUE,
  verbose = FALSE
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object
`name_overlap`     |     name for the overlap results (default to feat_info parameter)
`spatial_info`     |     polygon information
`feat_info`     |     feature information
`poly_ID_names`     |     list of poly_IDs to use
`polygon_group_size`     |     number of polygons to process per group
`return_gobject`     |     return giotto object (default: TRUE)
`verbose`     |     be verbose


## Details

Serial overlapping function that works on groups of polygons at a time.
 Number of polygons per group is defined by `polygon_group_size` param


## Value

giotto object or spatVector with overlapping information


