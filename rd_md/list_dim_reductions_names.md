# `list_dim_reductions_names`

list_dim_reductions_names


## Description

return the available dimension reductions object names


## Usage

```r
list_dim_reductions_names(
  gobject,
  data_type = "cells",
  spat_unit,
  feat_type,
  dim_type
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object
`data_type`     |     cells or feats dim reduction
`spat_unit`     |     spatial unit
`feat_type`     |     feature type
`dim_type`     |     dimensional reduction type (method)


## Details

function that can be used to find which names have been used


## Value

names pf dimension reduction object


