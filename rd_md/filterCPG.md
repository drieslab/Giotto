# `filterCPG`

filterCPG


## Description

Filter Interaction Changed Gene scores.


## Usage

```r
filterCPG(...)
```


## Arguments

Argument      |Description
------------- |----------------
`...`     |      Arguments passed on to [`filterICF`](#filtericf)   list("\n", "    ", list(list(list("cpgObject")), list("ICF (interaction changed feature) score object")), "\n", "    ", list(list(list("min_cells")), list("minimum number of source cell type")), "\n", "    ", list(list(list("min_cells_expr")), list("minimum expression level for source cell type")), "\n", "    ", list(list(list("min_int_cells")), list("minimum number of interacting neighbor cell type")), "\n", "    ", list(list(list("min_int_cells_expr")), list("minimum expression level for interacting neighbor cell type")), 
    "\n", "    ", list(list(list("min_fdr")), list("minimum adjusted p-value")), "\n", "    ", list(list(list("min_spat_diff")), list("minimum absolute spatial expression difference")), "\n", "    ", list(list(list("min_log2_fc")), list("minimum log2 fold-change")), "\n", "    ", list(list(list("min_zscore")), list("minimum z-score change")), "\n", "    ", list(list(list("zscores_column")), list("calculate z-scores over cell types or features")), "\n", "    ", list(list(list("direction")), list("differential expression directions to keep")), 
    "\n", "  ")


## Seealso

[`filterICF`](#filtericf)


