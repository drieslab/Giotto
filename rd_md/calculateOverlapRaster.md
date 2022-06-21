# `calculateOverlapRaster`

calculateOverlapRaster


## Description

calculate overlap between cellular structures (polygons) and features (points)


## Usage

```r
calculateOverlapRaster(
  gobject,
  name_overlap = NULL,
  spatial_info = NULL,
  poly_ID_names = NULL,
  feat_info = NULL,
  feat_subset_column = NULL,
  feat_subset_ids = NULL,
  return_gobject = TRUE,
  verbose = TRUE
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object
`name_overlap`     |     name for the overlap results (default to feat_info parameter)
`spatial_info`     |     polygon information
`poly_ID_names`     |     (optional) list of poly_IDs to use
`feat_info`     |     feature information
`feat_subset_column`     |     feature info column to subset features with
`feat_subset_ids`     |     ids within feature info column to use for subsetting
`return_gobject`     |     return giotto object (default: TRUE)
`verbose`     |     be verbose


## Details

Serial overlapping function.


## Value

giotto object or spatVector with overlapping information


