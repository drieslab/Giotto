# `filterICF`

filterICF


## Description

Filter Interaction Changed Feature scores.


## Usage

```r
filterICF(
  cpgObject,
  min_cells = 4,
  min_cells_expr = 1,
  min_int_cells = 4,
  min_int_cells_expr = 1,
  min_fdr = 0.1,
  min_spat_diff = 0.2,
  min_log2_fc = 0.2,
  min_zscore = 2,
  zscores_column = c("cell_type", "feats"),
  direction = c("both", "up", "down")
)
```


## Arguments

Argument      |Description
------------- |----------------
`cpgObject`     |     ICF (interaction changed feature) score object
`min_cells`     |     minimum number of source cell type
`min_cells_expr`     |     minimum expression level for source cell type
`min_int_cells`     |     minimum number of interacting neighbor cell type
`min_int_cells_expr`     |     minimum expression level for interacting neighbor cell type
`min_fdr`     |     minimum adjusted p-value
`min_spat_diff`     |     minimum absolute spatial expression difference
`min_log2_fc`     |     minimum log2 fold-change
`min_zscore`     |     minimum z-score change
`zscores_column`     |     calculate z-scores over cell types or features
`direction`     |     differential expression directions to keep


## Value

cpgObject that contains the filtered differential feature scores


