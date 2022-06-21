# `spatialAEH`

spatialAEH


## Description

Compute spatial variable genes with spatialDE method


## Usage

```r
spatialAEH(
  gobject = NULL,
  feat_type = NULL,
  spat_unit = NULL,
  spat_loc_name = "raw",
  SpatialDE_results = NULL,
  name_pattern = "AEH_patterns",
  expression_values = c("raw", "normalized", "scaled", "custom"),
  pattern_num = 6,
  l = 1.05,
  python_path = NULL,
  return_gobject = TRUE
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     Giotto object
`feat_type`     |     feature type
`spat_unit`     |     spatial unit
`spat_loc_name`     |     name for spatial locations
`SpatialDE_results`     |     results of [`spatialDE`](#spatialde) function
`name_pattern`     |     name for the computed spatial patterns
`expression_values`     |     gene expression values to use
`pattern_num`     |     number of spatial patterns to look for
`l`     |     lengthscale
`python_path`     |     specify specific path to python if required
`return_gobject`     |     show plot


## Details

This function is a wrapper for the SpatialAEH method implemented in the ...


## Value

An updated giotto object


