# `mergeClusters`

mergeClusters


## Description

Merge selected clusters based on pairwise correlation scores and size of cluster.


## Usage

```r
mergeClusters(
  gobject,
  spat_unit = NULL,
  feat_type = NULL,
  expression_values = c("normalized", "scaled", "custom"),
  cluster_column,
  cor = c("pearson", "spearman"),
  new_cluster_name = "merged_cluster",
  min_cor_score = 0.8,
  max_group_size = 20,
  force_min_group_size = 10,
  max_sim_clusters = 10,
  return_gobject = TRUE,
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
`cluster_column`     |     name of column to use for clusters
`cor`     |     correlation score to calculate distance
`new_cluster_name`     |     new name for merged clusters
`min_cor_score`     |     min correlation score to merge pairwise clusters
`max_group_size`     |     max cluster size that can be merged
`force_min_group_size`     |     size of clusters that will be merged with their most similar neighbor(s)
`max_sim_clusters`     |     maximum number of clusters to potentially merge to reach force_min_group_size
`return_gobject`     |     return giotto object
`verbose`     |     be verbose


## Details

Merge selected clusters based on pairwise correlation scores and size of cluster.
 To avoid large clusters to merge the max_group_size can be lowered. Small clusters can
 be forcibly merged with their most similar pairwise cluster by adjusting the
 force_min_group_size parameter. Clusters smaller than this value will be merged
 independent on the provided min_cor_score value. The force_min_group_size might not always
 be reached if clusters have already been merged before list() 
 A giotto object is returned by default, if FALSE then the merging vector will be returned.


## Value

Giotto object


