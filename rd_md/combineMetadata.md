# `combineMetadata`

combineMetadata


## Description

This function combines the cell metadata with spatial locations and
 enrichment results from [`runSpatialEnrich`](#runspatialenrich)


## Usage

```r
combineMetadata(
  gobject,
  spat_unit = NULL,
  feat_type = NULL,
  spat_loc_name = "raw",
  spat_enr_names = NULL
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     Giotto object
`spat_unit`     |     spatial unit
`feat_type`     |     feature type
`spat_loc_name`     |     name of spatial locations to include
`spat_enr_names`     |     names of spatial enrichment results to include


## Value

Extended cell metadata in data.table format.


