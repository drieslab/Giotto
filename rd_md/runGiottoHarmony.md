# `runGiottoHarmony`

runGiottoHarmony


## Description

run UMAP


## Usage

```r
runGiottoHarmony(
  gobject,
  spat_unit = NULL,
  feat_type = NULL,
  vars_use = "list_ID",
  do_pca = FALSE,
  expression_values = c("normalized", "scaled", "custom"),
  reduction = "cells",
  dim_reduction_to_use = "pca",
  dim_reduction_name = NULL,
  dimensions_to_use = 1:10,
  name = NULL,
  feats_to_use = NULL,
  toplevel_params = 2,
  return_gobject = TRUE,
  verbose = NULL,
  ...
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object
`spat_unit`     |     spatial unit
`feat_type`     |     feature type
`vars_use`     |     If meta_data is dataframe, this defined which variable(s) to remove (character vector).
`do_pca`     |     Whether to perform PCA on input matrix.
`expression_values`     |     expression values to use
`reduction`     |     reduction on cells or features
`dim_reduction_to_use`     |     use another dimension reduction set as input
`dim_reduction_name`     |     name of dimension reduction set to use
`dimensions_to_use`     |     number of dimensions to use as input
`name`     |     arbitrary name for Harmony run
`feats_to_use`     |     if dim_reduction_to_use = NULL, which genes to use
`toplevel_params`     |     parameters to extract
`return_gobject`     |     boolean: return giotto object (default = TRUE)
`verbose`     |     be verbose
`...`     |     additional Harmony parameters


## Details

See [`HarmonyMatrix`](#harmonymatrix) for more information about these and other parameters.
 This is a simple wrapper for the HarmonyMatrix function in the Harmony package.


## Value

giotto object with updated Harmony dimension recuction


