# `findMarkers`

findMarkers


## Description

Identify marker feats for selected clusters.


## Usage

```r
findMarkers(
  gobject,
  spat_unit = NULL,
  feat_type = NULL,
  expression_values = c("normalized", "scaled", "custom"),
  cluster_column = NULL,
  method = c("scran", "gini", "mast"),
  subset_clusters = NULL,
  group_1 = NULL,
  group_2 = NULL,
  min_expr_gini_score = 0.5,
  min_det_gini_score = 0.5,
  detection_threshold = 0,
  rank_score = 1,
  min_feats = 4,
  min_genes = NULL,
  group_1_name = NULL,
  group_2_name = NULL,
  adjust_columns = NULL,
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
`method`     |     method to use to detect differentially expressed feats
`subset_clusters`     |     selection of clusters to compare
`group_1`     |     group 1 cluster IDs from cluster_column for pairwise comparison
`group_2`     |     group 2 cluster IDs from cluster_column for pairwise comparison
`min_expr_gini_score`     |     gini: filter on minimum gini coefficient for expression
`min_det_gini_score`     |     gini: filter minimum gini coefficient for detection
`detection_threshold`     |     gini: detection threshold for feat expression
`rank_score`     |     gini: rank scores to include
`min_feats`     |     minimum number of top feats to return (for gini)
`min_genes`     |     deprecated, use min_feats
`group_1_name`     |     mast: custom name for group_1 clusters
`group_2_name`     |     mast: custom name for group_2 clusters
`adjust_columns`     |     mast: column in pDataDT to adjust for (e.g. detection rate)
`...`     |     additional parameters for the findMarkers function in scran or zlm function in MAST


## Details

Wrapper for all individual functions to detect marker feats for clusters.


## Value

data.table with marker feats


## Seealso

[`findScranMarkers`](#findscranmarkers) , [`findGiniMarkers`](#findginimarkers) and [`findMastMarkers`](#findmastmarkers)


