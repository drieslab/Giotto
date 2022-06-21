# `addCellStatistics`

addCellStatistics


## Description

adds cells statistics to the giotto object


## Usage

```r
addCellStatistics(
  gobject,
  feat_type = NULL,
  spat_unit = NULL,
  expression_values = c("normalized", "scaled", "custom"),
  detection_threshold = 0,
  return_gobject = TRUE
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object
`feat_type`     |     feature type
`spat_unit`     |     spatial unit
`expression_values`     |     expression values to use
`detection_threshold`     |     detection threshold to consider a gene detected
`return_gobject`     |     boolean: return giotto object (default = TRUE)


## Details

This function will add the following statistics to cell metadata:
   

*  nr_feats:  Denotes in how many features are detected per cell   

*  perc_feats:  Denotes what percentage of features is detected per cell   

*  total_expr:  Shows the total sum of feature expression per cell


## Value

giotto object if return_gobject = TRUE


