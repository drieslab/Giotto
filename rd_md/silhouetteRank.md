# `silhouetteRank`

silhouetteRank


## Description

Previously: calculate_spatial_genes_python. This method computes a silhouette score per gene based on the
 spatial distribution of two partitions of cells (expressed L1, and non-expressed L0).
 Here, rather than L2 Euclidean norm, it uses a rank-transformed, exponentially weighted
 function to represent the local physical distance between two cells.
 New multi aggregator implementation can be found at [`silhouetteRankTest`](#silhouetteranktest)


## Usage

```r
silhouetteRank(
  gobject,
  expression_values = c("normalized", "scaled", "custom"),
  metric = "euclidean",
  subset_genes = NULL,
  rbp_p = 0.95,
  examine_top = 0.3,
  python_path = NULL
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object
`expression_values`     |     expression values to use
`metric`     |     distance metric to use
`subset_genes`     |     only run on this subset of genes
`rbp_p`     |     fractional binarization threshold
`examine_top`     |     top fraction to evaluate with silhouette
`python_path`     |     specify specific path to python if required


## Value

data.table with spatial scores


