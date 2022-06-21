# `createHeatmap_DT`

createHeatmap_DT


## Description

creates order for clusters


## Usage

```r
createHeatmap_DT(
  gobject,
  spat_unit = NULL,
  feat_type = NULL,
  expression_values = c("normalized", "scaled", "custom"),
  feats,
  genes = NULL,
  cluster_column = NULL,
  cluster_order = c("size", "correlation", "custom"),
  cluster_custom_order = NULL,
  cluster_cor_method = "pearson",
  cluster_hclust_method = "ward.D",
  feat_order = c("correlation", "custom"),
  gene_order = NULL,
  feat_custom_order = NULL,
  gene_custom_order = NULL,
  feat_cor_method = "pearson",
  gene_cor_method = NULL,
  feat_hclust_method = "complete",
  gene_hclust_method = NULL
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object
`spat_unit`     |     spatial unit
`feat_type`     |     feature type
`expression_values`     |     expression values to use
`feats`     |     features to use
`genes`     |     deprecated, use feats
`cluster_column`     |     name of column to use for clusters
`cluster_order`     |     method to determine cluster order
`cluster_custom_order`     |     custom order for clusters
`cluster_cor_method`     |     method for cluster correlation
`cluster_hclust_method`     |     method for hierarchical clustering of clusters
`feat_order`     |     method to determine features order
`gene_order`     |     deprecated, use feat_order in the future
`feat_custom_order`     |     custom order for features
`gene_custom_order`     |     deprecated, use feat_custom_order in the future
`feat_cor_method`     |     method for features correlation
`gene_cor_method`     |     deprecated, use feat_cor_method in the future
`feat_hclust_method`     |     method for hierarchical clustering of features
`gene_hclust_method`     |     deprecated, use feat_hclust_method in the future


## Details

Creates input data.tables for plotHeatmap function.


## Value

list


