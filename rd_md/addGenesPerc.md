# `addGenesPerc`

addGenesPerc


## Description

calculates the total percentage of (normalized) counts for a subset of selected genes


## Usage

```r
addGenesPerc(
  gobject,
  spat_unit = NULL,
  feat_type = NULL,
  expression_values = c("normalized", "scaled", "custom"),
  genes = NULL,
  vector_name = "gene_perc",
  return_gobject = TRUE
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object
`spat_unit`     |     spatial unit
`feat_type`     |     feature type
`expression_values`     |     expression values to use
`genes`     |     vector of selected genes
`vector_name`     |     column name as seen in [`pDataDT`](#pdatadt)
`return_gobject`     |     boolean: return giotto object (default = TRUE)


## Value

giotto object if `return_gobject = TRUE` , else a vector with % results


