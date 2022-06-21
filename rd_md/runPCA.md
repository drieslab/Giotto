# `runPCA`

runPCA


## Description

runs a Principal Component Analysis


## Usage

```r
runPCA(
  gobject,
  spat_unit = NULL,
  feat_type = NULL,
  expression_values = c("normalized", "scaled", "custom"),
  reduction = c("cells", "feats"),
  name = NULL,
  feats_to_use = "hvf",
  genes_to_use = NULL,
  return_gobject = TRUE,
  center = TRUE,
  scale_unit = TRUE,
  ncp = 100,
  method = c("irlba", "exact", "random", "factominer"),
  method_params = list(NA),
  rev = FALSE,
  set_seed = TRUE,
  seed_number = 1234,
  verbose = TRUE,
  ...
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object
`spat_unit`     |     spatial unit
`feat_type`     |     feature type
`expression_values`     |     expression values to use
`reduction`     |     cells or genes
`name`     |     arbitrary name for PCA run
`feats_to_use`     |     subset of features to use for PCA
`genes_to_use`     |     deprecated use feats_to_use
`return_gobject`     |     boolean: return giotto object (default = TRUE)
`center`     |     center data first (default = TRUE)
`scale_unit`     |     scale features before PCA (default = TRUE)
`ncp`     |     number of principal components to calculate
`method`     |     which implementation to use
`method_params`     |     additional parameters
`rev`     |     do a reverse PCA
`set_seed`     |     use of seed
`seed_number`     |     seed number to use
`verbose`     |     verbosity of the function
`...`     |     additional parameters for PCA (see details)


## Details

See [`runPCA`](#runpca) and [`PCA`](#pca) for more information about other parameters.
   

*  feats_to_use = NULL: will use all features from the selected matrix  

*  feats_to_use = <hvg name>: can be used to select a column name of highly variable features, created by (see [`calculateHVF`](#calculatehvf) )  

*  feats_to_use = c('geneA', 'geneB', ...): will use all manually provided features


## Value

giotto object with updated PCA dimension recuction


