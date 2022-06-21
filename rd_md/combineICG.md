# `combineICG`

combineICG


## Description

Combine ICG scores in a pairwise manner.


## Usage

```r
combineICG(
  cpgObject,
  selected_ints = NULL,
  selected_genes = NULL,
  specific_genes_1 = NULL,
  specific_genes_2 = NULL,
  min_cells = 5,
  min_int_cells = 3,
  min_fdr = 0.05,
  min_spat_diff = 0,
  min_log2_fc = 0.5,
  do_parallel = TRUE,
  verbose = T
)
```


## Arguments

Argument      |Description
------------- |----------------
`cpgObject`     |     ICG (interaction changed gene) score object
`selected_ints`     |     subset of selected cell-cell interactions (optional)
`selected_genes`     |     subset of selected genes (optional)
`specific_genes_1`     |     specific geneset combo (need to position match specific_genes_2)
`specific_genes_2`     |     specific geneset combo (need to position match specific_genes_1)
`min_cells`     |     minimum number of target cell type
`min_int_cells`     |     minimum number of interacting cell type
`min_fdr`     |     minimum adjusted p-value
`min_spat_diff`     |     minimum absolute spatial expression difference
`min_log2_fc`     |     minimum absolute log2 fold-change
`do_parallel`     |     run calculations in parallel with mclapply
`verbose`     |     verbose


## Value

cpgObject that contains the filtered differential gene scores


