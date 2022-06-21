# `findMarkers_one_vs_all`

findMarkers_one_vs_all


## Description

Identify marker feats for all clusters in a one vs all manner.


## Usage

```r
findMarkers_one_vs_all(
  gobject,
  feat_type = NULL,
  spat_unit = NULL,
  expression_values = c("normalized", "scaled", "custom"),
  cluster_column,
  subset_clusters = NULL,
  method = c("scran", "gini", "mast"),
  pval = 0.01,
  logFC = 0.5,
  min_feats = 10,
  min_genes = NULL,
  min_expr_gini_score = 0.5,
  min_det_gini_score = 0.5,
  detection_threshold = 0,
  rank_score = 1,
  adjust_columns = NULL,
  verbose = TRUE,
  ...
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object
`feat_type`     |     feature type
`spat_unit`     |     spatial unit
`expression_values`     |     feat expression values to use
`cluster_column`     |     clusters to use
`subset_clusters`     |     selection of clusters to compare
`method`     |     method to use to detect differentially expressed feats
`pval`     |     scran & mast: filter on minimal p-value
`logFC`     |     scan & mast: filter on logFC
`min_feats`     |     minimum feats to keep per cluster, overrides pval and logFC
`min_genes`     |     deprecated, use min_feats
`min_expr_gini_score`     |     gini: filter on minimum gini coefficient for expression
`min_det_gini_score`     |     gini: filter minimum gini coefficient for detection
`detection_threshold`     |     gini: detection threshold for feat expression
`rank_score`     |     gini: rank scores to include
`adjust_columns`     |     mast: column in pDataDT to adjust for (e.g. detection rate)
`verbose`     |     be verbose
`...`     |     additional parameters for the findMarkers function in scran or zlm function in MAST


## Details

Wrapper for all one vs all functions to detect marker feats for clusters.


## Value

data.table with marker feats


## Seealso

[`findScranMarkers_one_vs_all`](#findscranmarkersonevsall) , [`findGiniMarkers_one_vs_all`](#findginimarkersonevsall) and [`findMastMarkers_one_vs_all`](#findmastmarkersonevsall)


