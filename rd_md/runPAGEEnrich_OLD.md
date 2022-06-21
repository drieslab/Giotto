# `runPAGEEnrich_OLD`

runPAGEEnrich_OLD


## Description

Function to calculate gene signature enrichment scores per spatial position using PAGE.


## Usage

```r
runPAGEEnrich_OLD(
  gobject,
  sign_matrix,
  expression_values = c("normalized", "scaled", "custom"),
  reverse_log_scale = TRUE,
  logbase = 2,
  output_enrichment = c("original", "zscore"),
  p_value = FALSE,
  n_times = 1000,
  name = NULL,
  return_gobject = TRUE
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     Giotto object
`sign_matrix`     |     Matrix of signature genes for each cell type / process
`expression_values`     |     expression values to use
`reverse_log_scale`     |     reverse expression values from log scale
`logbase`     |     log base to use if reverse_log_scale = TRUE
`output_enrichment`     |     how to return enrichment output
`p_value`     |     calculate p-values (boolean, default = FALSE)
`n_times`     |     number of permutations to calculate for p_value
`name`     |     to give to spatial enrichment results, default = PAGE
`return_gobject`     |     return giotto object


## Details

sign_matrix: a binary matrix with genes as row names and cell-types as column names.
 Alternatively a list of signature genes can be provided to makeSignMatrixPAGE, which will create
 the matrix for you. list() 
 
 The enrichment Z score is calculated by using method (PAGE) from
 Kim SY et al., BMC bioinformatics, 2005 as $Z = ((Sm – mu)*m^(1/2)) / delta$ .
 For each gene in each spot, mu is the fold change values versus the mean expression
 and delta is the standard deviation. Sm is the mean fold change value of a specific marker gene set
 and  m is the size of a given marker gene set.


## Value

data.table with enrichment results


## Seealso

[`makeSignMatrixPAGE`](#makesignmatrixpage)


