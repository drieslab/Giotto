# `hyperGeometricEnrich`

hyperGeometricEnrich


## Description

Function to calculate gene signature enrichment scores per spatial position using a hypergeometric test.


## Usage

```r
hyperGeometricEnrich(...)
```


## Arguments

Argument      |Description
------------- |----------------
`...`     |      Arguments passed on to [`runHyperGeometricEnrich`](#runhypergeometricenrich)   list("\n", "    ", list(list(list("gobject")), list("Giotto object")), "\n", "    ", list(list(list("spat_unit")), list("spatial unit")), "\n", "    ", list(list(list("feat_type")), list("feature type")), "\n", "    ", list(list(list("sign_matrix")), list("Matrix of signature genes for each cell type / process")), "\n", "    ", list(list(list("expression_values")), list("expression values to use")), "\n", "    ", list(list(list("reverse_log_scale")), list("reverse expression values from log scale")), 
    "\n", "    ", list(list(list("logbase")), list("log base to use if reverse_log_scale = TRUE")), "\n", "    ", list(list(list("top_percentage")), list("percentage of cells that will be considered to have gene expression with matrix binarization")), "\n", "    ", list(list(list("output_enrichment")), list("how to return enrichment output")), "\n", "    ", list(list(list("p_value")), list("calculate p-values (boolean, default = FALSE)")), "\n", "    ", list(list(list("name")), list("to give to spatial enrichment results, default = rank")), 
    "\n", "    ", list(list(list("return_gobject")), list("return giotto object")), "\n", "  ")


## Seealso

[`runHyperGeometricEnrich`](#runhypergeometricenrich)


