# `detectSpatialCorGenes`

detectSpatialCorGenes


## Description

Detect genes that are spatially correlated


## Usage

```r
detectSpatialCorGenes(
  gobject,
  feat_type = NULL,
  spat_unit = NULL,
  method = c("grid", "network"),
  expression_values = c("normalized", "scaled", "custom"),
  subset_feats = NULL,
  subset_genes = NULL,
  spatial_network_name = "Delaunay_network",
  network_smoothing = NULL,
  spatial_grid_name = "spatial_grid",
  min_cells_per_grid = 4,
  cor_method = c("pearson", "kendall", "spearman")
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object
`feat_type`     |     feature type
`spat_unit`     |     spatial unit
`method`     |     method to use for spatial averaging
`expression_values`     |     gene expression values to use
`subset_feats`     |     subset of feats to use
`subset_genes`     |     deprecated, use `subset_feats`
`spatial_network_name`     |     name of spatial network to use
`network_smoothing`     |     smoothing factor beteen 0 and 1 (default: automatic)
`spatial_grid_name`     |     name of spatial grid to use
`min_cells_per_grid`     |     minimum number of cells to consider a grid
`cor_method`     |     correlation method


## Details

For method = network, it expects a fully connected spatial network. You can make sure to create a
 fully connected network by setting minimal_k > 0 in the [`createSpatialNetwork`](#createspatialnetwork) function.
   

*  list("1. grid-averaging: ") list("average gene expression values within a predefined spatial grid")   

*  list("2. network-averaging: ") list("smoothens the gene expression matrix by averaging the expression within one cell\n", " by using the neighbours within the predefined spatial network. b is a smoothening factor\n", " that defaults to 1 - 1/k, where k is the median number of  k-neighbors in the\n", " selected spatial network. Setting b = 0 means no smoothing and b = 1 means no contribution\n", " from its own expression.")  
 The spatCorObject can be further explored with showSpatialCorGenes()


## Value

returns a spatial correlation object: "spatCorObject"


## Seealso

[`showSpatialCorGenes`](#showspatialcorgenes)


