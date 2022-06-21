# `showSpatialCorGenes`

showSpatialCorGenes


## Description

Shows and filters spatially correlated genes


## Usage

```r
showSpatialCorGenes(
  spatCorObject,
  use_clus_name = NULL,
  selected_clusters = NULL,
  genes = NULL,
  min_spat_cor = 0.5,
  min_expr_cor = NULL,
  min_cor_diff = NULL,
  min_rank_diff = NULL,
  show_top_genes = NULL
)
```


## Arguments

Argument      |Description
------------- |----------------
`spatCorObject`     |     spatial correlation object
`use_clus_name`     |     cluster information to show
`selected_clusters`     |     subset of clusters to show
`genes`     |     subset of genes to show
`min_spat_cor`     |     filter on minimum spatial correlation
`min_expr_cor`     |     filter on minimum single-cell expression correlation
`min_cor_diff`     |     filter on minimum correlation difference (spatial vs expression)
`min_rank_diff`     |     filter on minimum correlation rank difference (spatial vs expression)
`show_top_genes`     |     show top genes per gene


## Value

data.table with filtered information


