# `findGiniMarkers_one_vs_all`

findGiniMarkers_one_vs_all


## Description

Identify marker feats for all clusters in a one vs all manner based on gini detection and expression scores.


## Usage

```r
findGiniMarkers_one_vs_all(
  gobject,
  feat_type = NULL,
  spat_unit = NULL,
  expression_values = c("normalized", "scaled", "custom"),
  cluster_column,
  subset_clusters = NULL,
  min_expr_gini_score = 0.5,
  min_det_gini_score = 0.5,
  detection_threshold = 0,
  rank_score = 1,
  min_feats = 4,
  min_genes = NULL,
  verbose = TRUE
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
`min_expr_gini_score`     |     filter on minimum gini coefficient on expression
`min_det_gini_score`     |     filter on minimum gini coefficient on detection
`detection_threshold`     |     detection threshold for feat expression
`rank_score`     |     rank scores for both detection and expression to include
`min_feats`     |     minimum number of top feats to return
`min_genes`     |     deprecated, use min_feats
`verbose`     |     be verbose


## Value

data.table with marker feats


## Seealso

[`findGiniMarkers`](#findginimarkers)


