# `create_average_detection_DT`

create_average_detection_DT


## Description

calculates average gene detection for a cell metadata factor (e.g. cluster)


## Usage

```r
create_average_detection_DT(
  gobject,
  feat_type = NULL,
  spat_unit = NULL,
  meta_data_name,
  expression_values = c("normalized", "scaled", "custom"),
  detection_threshold = 0
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object
`feat_type`     |     feature type
`spat_unit`     |     spatial unit
`meta_data_name`     |     name of metadata column to use
`expression_values`     |     which expression values to use
`detection_threshold`     |     detection threshold to consider a gene detected


## Value

data.table with average gene epression values for each factor


