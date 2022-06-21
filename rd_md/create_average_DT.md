# `create_average_DT`

create_average_DT


## Description

calculates average gene expression for a cell metadata factor (e.g. cluster)


## Usage

```r
create_average_DT(
  gobject,
  spat_unit = NULL,
  feat_type = NULL,
  meta_data_name,
  expression_values = c("normalized", "scaled", "custom")
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object
`spat_unit`     |     spatial unit
`feat_type`     |     feature type
`meta_data_name`     |     name of metadata column to use
`expression_values`     |     which expression values to use


## Value

data.table with average gene epression values for each factor


