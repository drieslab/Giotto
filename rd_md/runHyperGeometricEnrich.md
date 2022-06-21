# `runHyperGeometricEnrich`

runHyperGeometricEnrich


## Description

Function to calculate gene signature enrichment scores per spatial position using a hypergeometric test.


## Usage

```r
runHyperGeometricEnrich(
  gobject,
  spat_unit = NULL,
  feat_type = NULL,
  sign_matrix,
  expression_values = c("normalized", "scaled", "custom"),
  reverse_log_scale = TRUE,
  logbase = 2,
  top_percentage = 5,
  output_enrichment = c("original", "zscore"),
  p_value = FALSE,
  name = NULL,
  return_gobject = TRUE
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     Giotto object
`spat_unit`     |     spatial unit
`feat_type`     |     feature type
`sign_matrix`     |     Matrix of signature genes for each cell type / process
`expression_values`     |     expression values to use
`reverse_log_scale`     |     reverse expression values from log scale
`logbase`     |     log base to use if reverse_log_scale = TRUE
`top_percentage`     |     percentage of cells that will be considered to have gene expression with matrix binarization
`output_enrichment`     |     how to return enrichment output
`p_value`     |     calculate p-values (boolean, default = FALSE)
`name`     |     to give to spatial enrichment results, default = rank
`return_gobject`     |     return giotto object


## Details

The enrichment score is calculated based on the p-value from the
 hypergeometric test, -log10(p-value).


## Value

data.table with enrichment results


