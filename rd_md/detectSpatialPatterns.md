# `detectSpatialPatterns`

detectSpatialPatterns


## Description

Identify spatial patterns through PCA on average expression in a spatial grid.


## Usage

```r
detectSpatialPatterns(
  gobject,
  expression_values = c("normalized", "scaled", "custom"),
  spatial_grid_name = "spatial_grid",
  min_cells_per_grid = 4,
  scale_unit = F,
  ncp = 100,
  show_plot = T,
  PC_zscore = 1.5
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object
`expression_values`     |     expression values to use
`spatial_grid_name`     |     name of spatial grid to use (default = 'spatial_grid')
`min_cells_per_grid`     |     minimum number of cells in a grid to be considered
`scale_unit`     |     scale features
`ncp`     |     number of principal components to calculate
`show_plot`     |     show plots
`PC_zscore`     |     minimum z-score of variance explained by a PC


## Details

Steps to identify spatial patterns:
   

*  1. average gene expression for cells within a grid, see createSpatialGrid   

*  2. perform PCA on the average grid expression profiles   

*  3. convert variance of principlal components (PCs) to z-scores and select PCs based on a z-score threshold


## Value

spatial pattern object 'spatPatObj'


