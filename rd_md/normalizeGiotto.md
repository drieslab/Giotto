# `normalizeGiotto`

normalizeGiotto


## Description

fast normalize and/or scale expresion values of Giotto object


## Usage

```r
normalizeGiotto(
  gobject,
  spat_unit = NULL,
  feat_type = NULL,
  expression_values = "raw",
  norm_methods = c("standard", "pearson_resid", "osmFISH"),
  library_size_norm = TRUE,
  scalefactor = 6000,
  log_norm = TRUE,
  log_offset = 1,
  logbase = 2,
  scale_feats = TRUE,
  scale_genes = NULL,
  scale_cells = TRUE,
  scale_order = c("first_feats", "first_cells"),
  theta = 100,
  update_slot = "scaled",
  verbose = TRUE
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object
`spat_unit`     |     spatial unit
`feat_type`     |     feature type
`expression_values`     |     expression values to use
`norm_methods`     |     normalization method to use
`library_size_norm`     |     normalize cells by library size
`scalefactor`     |     scale factor to use after library size normalization
`log_norm`     |     transform values to log-scale
`log_offset`     |     offset value to add to expression matrix, default = 1
`logbase`     |     log base to use to log normalize expression values
`scale_feats`     |     z-score genes over all cells
`scale_genes`     |     deprecated, use scale_feats
`scale_cells`     |     z-score cells over all genes
`scale_order`     |     order to scale feats and cells
`theta`     |     theta parameter for the pearson residual normalization step
`update_slot`     |     slot or name to use for the results from osmFISH and pearson residual normalization
`verbose`     |     be verbose


## Details

Currently there are two 'methods' to normalize your raw counts data.
 
 A. The standard method follows the standard protocol which can be adjusted using
 the provided parameters and follows the following order: list() 
   

*  1. Data normalization for total library size and scaling by a custom scale-factor.   

*  2. Log transformation of data.   

*  3. Z-scoring of data by genes and/or cells.  
 B. The normalization method as provided by the osmFISH paper is also implemented: list() 
   

*  list("1. First normalize genes, for each gene divide the counts by the total gene count and\n", "multiply by the total number of genes.")   

*  list("2. Next normalize cells, for each cell divide the normalized gene counts by the total\n", "counts per cell and multiply by the total number of cells.")  
 C. The normalization method as provided by Lause/Kobak et al is also implemented: list() 
   

*  1. First calculate expected values based on Pearson correlations.   

*  2. Next calculate z-scores based on observed and expected values.  
 By default the latter two results will be saved in the Giotto slot for scaled expression,
 this can be changed by changing the update_slot parameters


## Value

giotto object


