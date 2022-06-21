# `filterGiotto`

filterGiotto


## Description

filter Giotto object based on expression threshold


## Usage

```r
filterGiotto(
  gobject,
  spat_unit = NULL,
  feat_type = NULL,
  expression_values = c("raw", "normalized", "scaled", "custom"),
  expression_threshold = 1,
  feat_det_in_min_cells = 100,
  gene_det_in_min_cells = NULL,
  min_det_feats_per_cell = 100,
  min_det_genes_per_cell = NULL,
  poly_info = "cell",
  verbose = TRUE
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object
`spat_unit`     |     spatial unit
`feat_type`     |     feature type
`expression_values`     |     expression values to use
`expression_threshold`     |     threshold to consider a gene expressed
`feat_det_in_min_cells`     |     minimum # of cells that need to express a feature
`gene_det_in_min_cells`     |     deprecated, use feat_det_in_min_cells
`min_det_feats_per_cell`     |     minimum # of features that need to be detected in a cell
`min_det_genes_per_cell`     |     deprecated, use min_det_feats_per_cell
`poly_info`     |     polygon information to use
`verbose`     |     verbose


## Details

The function [`filterCombinations`](#filtercombinations) can be used to explore the effect of different parameter values.


## Value

giotto object


