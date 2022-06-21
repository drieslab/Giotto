# `overlapToMatrix`

overlapToMatrix


## Description

create a count matrix based on overlap results from [`calculateOverlapRaster`](#calculateoverlapraster) , [`calculateOverlapSerial`](#calculateoverlapserial) , or [`calculateOverlapParallel`](#calculateoverlapparallel)


## Usage

```r
overlapToMatrix(
  gobject,
  name = "raw",
  poly_info = "cell",
  feat_info = "rna",
  return_gobject = TRUE
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object
`name`     |     name for the overlap count matrix
`poly_info`     |     polygon information
`feat_info`     |     feature information
`return_gobject`     |     return giotto object (default: TRUE)


## Value

giotto object or count matrix


