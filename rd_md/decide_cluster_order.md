# `decide_cluster_order`

decide_cluster_order


## Description

creates order for clusters


## Usage

```r
decide_cluster_order(
  gobject,
  spat_unit = NULL,
  feat_type = NULL,
  expression_values = c("normalized", "scaled", "custom"),
  feats,
  genes = NULL,
  cluster_column = NULL,
  cluster_order = c("size", "correlation", "custom"),
  cluster_custom_order = NULL,
  cor_method = "pearson",
  hclust_method = "ward.D"
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object
`spat_unit`     |     spatial unit
`feat_type`     |     feature type
`expression_values`     |     expression values to use
`feats`     |     features to use (e.g. genes)
`genes`     |     deprecated, use feats
`cluster_column`     |     name of column to use for clusters
`cluster_order`     |     method to determine cluster order
`cluster_custom_order`     |     custom order for clusters
`cor_method`     |     method for correlation
`hclust_method`     |     method for hierarchical clustering


## Details

Calculates order for clusters.


## Value

custom


