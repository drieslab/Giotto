# `create_dimObject`

create_dimObject


## Description

Creates an object that stores a dimension reduction output


## Usage

```r
create_dimObject(
  name = "test",
  spat_unit = "cell",
  feat_type = "rna",
  reduction_method = NULL,
  coordinates = NULL,
  misc = NULL,
  my_rownames = NULL
)
```


## Arguments

Argument      |Description
------------- |----------------
`name`     |     arbitrary name for object
`spat_unit`     |     spatial unit
`feat_type`     |     feature type
`reduction_method`     |     method used to reduce dimensions
`coordinates`     |     accepts the coordinates after dimension reduction
`misc`     |     any additional information will be added to this slot
`my_rownames`     |     rownames


## Value

number of distinct colors


