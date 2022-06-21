# `do_spatial_knn_smoothing`

do_spatial_knn_smoothing


## Description

smooth gene expression over a kNN spatial network


## Usage

```r
do_spatial_knn_smoothing(
  expression_matrix,
  spatial_network,
  subset_feats = NULL,
  b = NULL
)
```


## Arguments

Argument      |Description
------------- |----------------
`subset_feats`     |     subset of features to use
`b`     |     smoothing factor beteen 0 and 1 (default: automatic)
`gobject`     |     giotto object
`expression_values`     |     gene expression values to use
`spatial_network_name`     |     name of spatial network to use


## Details

This function will smoothen the gene expression values per cell according to
 its neighbors in the selected spatial network. list() 
 b is a smoothening factor that defaults to 1 - 1/k, where k is the median number of
 k-neighbors in the selected spatial network. Setting b = 0 means no smoothing and b = 1
 means no contribution from its own expression.


## Value

matrix with smoothened gene expression values based on kNN spatial network


