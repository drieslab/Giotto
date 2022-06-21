# `runtSNE`

Run tSNE dimensional reduction


## Description

run tSNE


## Usage

```r
runtSNE(
  gobject,
  spat_unit = NULL,
  feat_type = NULL,
  expression_values = c("normalized", "scaled", "custom"),
  reduction = c("cells", "feats"),
  dim_reduction_to_use = "pca",
  dim_reduction_name = NULL,
  dimensions_to_use = 1:10,
  name = NULL,
  feats_to_use = NULL,
  genes_to_use = NULL,
  return_gobject = TRUE,
  dims = 2,
  perplexity = 30,
  theta = 0.5,
  do_PCA_first = F,
  set_seed = T,
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
`dim_reduction_to_use`     |     use another dimension reduction set as input
`dim_reduction_name`     |     name of dimension reduction set to use
`dimensions_to_use`     |     number of dimensions to use as input
`name`     |     arbitrary name for tSNE run
`feats_to_use`     |     if dim_reduction_to_use = NULL, which genes to use
`genes_to_use`     |     deprecated, use feats_to_use
`return_gobject`     |     boolean: return giotto object (default = TRUE)
`dims`     |     tSNE param: number of dimensions to return
`perplexity`     |     tSNE param: perplexity
`theta`     |     tSNE param: theta
`do_PCA_first`     |     tSNE param: do PCA before tSNE (default = FALSE)
`set_seed`     |     use of seed
`seed_number`     |     seed number to use
`verbose`     |     verbosity of the function
`...`     |     additional tSNE parameters


## Details

See [`Rtsne`](#rtsne) for more information about these and other parameters. list() 
   

*  Input for tSNE dimension reduction can be another dimension reduction (default = 'pca')  

*  To use gene expression as input set dim_reduction_to_use = NULL  

*  If dim_reduction_to_use = NULL, genes_to_use can be used to select a column name of highly variable genes (see [`calculateHVF`](#calculatehvf) ) or simply provide a vector of genes  

*  multiple tSNE results can be stored by changing the name of the analysis


## Value

giotto object with updated tSNE dimension recuction


