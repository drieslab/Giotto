# `findMastMarkers`

findMastMarkers


## Description

Identify marker feats for selected clusters based on the MAST package.


## Usage

```r
findMastMarkers(
  gobject,
  feat_type = NULL,
  spat_unit = NULL,
  expression_values = c("normalized", "scaled", "custom"),
  cluster_column,
  group_1 = NULL,
  group_1_name = NULL,
  group_2 = NULL,
  group_2_name = NULL,
  adjust_columns = NULL,
  verbose = FALSE,
  ...
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
`group_1`     |     group 1 cluster IDs from cluster_column for pairwise comparison
`group_1_name`     |     custom name for group_1 clusters
`group_2`     |     group 2 cluster IDs from cluster_column for pairwise comparison
`group_2_name`     |     custom name for group_2 clusters
`adjust_columns`     |     column in pDataDT to adjust for (e.g. detection rate)
`verbose`     |     be verbose
`...`     |     additional parameters for the zlm function in MAST


## Details

This is a minimal convenience wrapper around the [`zlm`](#zlm) 
 from the MAST package to detect differentially expressed feats. Caution: with large datasets
 MAST might take a long time to run and finish


## Value

data.table with marker feats


