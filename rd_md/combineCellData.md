# `combineCellData`

combineCellData


## Description

combine cell data information


## Usage

```r
combineCellData(
  gobject,
  feat_type = "rna",
  include_spat_locs = TRUE,
  spat_loc_name = "raw",
  include_poly_info = TRUE,
  poly_info = "cell"
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object
`feat_type`     |     feature type
`include_spat_locs`     |     include information about spatial locations
`spat_loc_name`     |     spatial location name
`include_poly_info`     |     include information about polygon
`poly_info`     |     polygon information name


## Value

data.table with combined spatial information


