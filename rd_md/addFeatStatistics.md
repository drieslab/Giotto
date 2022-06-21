# `addFeatStatistics`

Add feature statistics


## Description

Adds feature statistics to the giotto object


## Usage

```r
addFeatStatistics(
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

This function will add the following statistics to feature metadata:
   

*  nr_cells:  Denotes in how many cells the feature is detected   

*  per_cells:  Denotes in what percentage of cells the feature is detected   

*  total_expr:  Shows the total sum of feature expression in all cells   

*  mean_expr:  Average feature expression in all cells   

*  mean_expr_det:  Average feature expression in cells with detectable levels of the gene


## Value

giotto object if return_gobject = TRUE


