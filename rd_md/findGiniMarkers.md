# `findGiniMarkers`

findGiniMarkers


## Description

Identify marker feats for selected clusters based on gini detection and expression scores.


## Usage

```r
findGiniMarkers(
  gobject,
  feat_type = NULL,
  spat_unit = NULL,
  expression_values = c("normalized", "scaled", "custom"),
  cluster_column,
  subset_clusters = NULL,
  group_1 = NULL,
  group_1_name = NULL,
  group_2 = NULL,
  group_2_name = NULL,
  min_expr_gini_score = 0.2,
  min_det_gini_score = 0.2,
  detection_threshold = 0,
  rank_score = 1,
  min_feats = 5,
  min_genes = NULL
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
`group_1`     |     group 1 cluster IDs from cluster_column for pairwise comparison
`group_1_name`     |     custom name for group_1 clusters
`group_2`     |     group 2 cluster IDs from cluster_column for pairwise comparison
`group_2_name`     |     custom name for group_2 clusters
`min_expr_gini_score`     |     filter on minimum gini coefficient for expression
`min_det_gini_score`     |     filter on minimum gini coefficient for detection
`detection_threshold`     |     detection threshold for feat expression
`rank_score`     |     rank scores for both detection and expression to include
`min_feats`     |     minimum number of top feats to return
`min_genes`     |     deprecated, use min_feats


## Details

Detection of marker feats using the [https://en.wikipedia.org/wiki/Gini_coefficient](https://en.wikipedia.org/wiki/Gini_coefficient) gini 
 coefficient is based on the following steps/principles per feat:
   

*  1. calculate average expression per cluster   

*  2. calculate detection fraction per cluster   

*  3. calculate gini-coefficient for av. expression values over all clusters   

*  4. calculate gini-coefficient for detection fractions over all clusters   

*  5. convert gini-scores to rank scores   

*  6. for each feat create combined score = detection rank x expression rank x expr gini-coefficient x detection gini-coefficient   

*  7. for each feat sort on expression and detection rank and combined score  
 
 As a results "top gini" feats are feats that are very selectivily expressed in a specific cluster,
 however not always expressed in all cells of that cluster. In other words highly specific, but
 not necessarily sensitive at the single-cell level.
 
 To perform differential expression between custom selected groups of cells
 you need to specify the cell_ID column to parameter cluster_column 
 and provide the individual cell IDs to the parameters group_1 and group_2 
 
 By default group names will be created by pasting the different id names within each selected group.
 When you have many different ids in a single group
 it is recommend to provide names for both groups to group_1_name and group_2_name


## Value

data.table with marker feats


