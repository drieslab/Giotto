# `average_feat_feat_expression_in_groups`

average_feat_feat_expression_in_groups


## Description

calculate average expression per cluster


## Usage

```r
average_feat_feat_expression_in_groups(
  gobject,
  spat_unit = NULL,
  feat_type = NULL,
  cluster_column = "cell_types",
  feat_set_1,
  feat_set_2
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object to use
`spat_unit`     |     spatial unit
`feat_type`     |     feature type
`cluster_column`     |     cluster column with cell type information
`feat_set_1`     |     first specific feat set from feat pairs
`feat_set_2`     |     second specific feat set from feat pairs


## Value

data.table with average expression scores for each cluster


