# `findScranMarkers_one_vs_all`

findScranMarkers_one_vs_all


## Description

Identify marker feats for all clusters in a one vs all manner based on scran's implementation of findMarkers.


## Usage

```r
findScranMarkers_one_vs_all(
  gobject,
  spat_unit = NULL,
  feat_type = NULL,
  expression_values = c("normalized", "scaled", "custom"),
  cluster_column,
  subset_clusters = NULL,
  pval = 0.01,
  logFC = 0.5,
  min_feats = 10,
  min_genes = NULL,
  verbose = TRUE,
  ...
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object
`spat_unit`     |     spatial unit
`feat_type`     |     feature type
`expression_values`     |     feat expression values to use
`cluster_column`     |     clusters to use
`subset_clusters`     |     subset of clusters to use
`pval`     |     filter on minimal p-value
`logFC`     |     filter on logFC
`min_feats`     |     minimum feats to keep per cluster, overrides pval and logFC
`min_genes`     |     deprecated, use min_feats
`verbose`     |     be verbose
`...`     |     additional parameters for the findMarkers function in scran


## Value

data.table with marker feats


## Seealso

[`findScranMarkers`](#findscranmarkers)


