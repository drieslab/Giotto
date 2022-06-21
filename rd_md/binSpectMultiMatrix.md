# `binSpectMultiMatrix`

binSpectMultiMatrix


## Description

binSpect for a single spatial network and a provided expression matrix


## Usage

```r
binSpectMultiMatrix(
  expression_matrix,
  spatial_networks,
  bin_method = c("kmeans", "rank"),
  subset_feats = NULL,
  kmeans_algo = c("kmeans", "kmeans_arma", "kmeans_arma_subset"),
  nstart = 3,
  iter_max = 10,
  extreme_nr = 50,
  sample_nr = 50,
  percentage_rank = c(10, 30),
  do_fisher_test = TRUE,
  adjust_method = "fdr",
  calc_hub = FALSE,
  hub_min_int = 3,
  get_av_expr = TRUE,
  get_high_expr = TRUE,
  implementation = c("data.table", "simple", "matrix"),
  group_size = "automatic",
  do_parallel = TRUE,
  cores = NA,
  verbose = T,
  knn_params = NULL,
  set.seed = NULL,
  summarize = c("adj.p.value", "p.value")
)
```


## Arguments

Argument      |Description
------------- |----------------
`expression_matrix`     |     expression matrix
`spatial_networks`     |     list of spatial networks in data.table format
`bin_method`     |     method to binarize gene expression
`subset_feats`     |     only select a subset of features to test
`kmeans_algo`     |     kmeans algorithm to use (kmeans, kmeans_arma, kmeans_arma_subset)
`nstart`     |     kmeans: nstart parameter
`iter_max`     |     kmeans: iter.max parameter
`extreme_nr`     |     number of top and bottom cells (see details)
`sample_nr`     |     total number of cells to sample (see details)
`percentage_rank`     |     vector of percentages of top cells for binarization
`do_fisher_test`     |     perform fisher test
`adjust_method`     |     p-value adjusted method to use (see [`p.adjust`](#p.adjust) )
`calc_hub`     |     calculate the number of hub cells
`hub_min_int`     |     minimum number of cell-cell interactions for a hub cell
`get_av_expr`     |     calculate the average expression per gene of the high expressing cells
`get_high_expr`     |     calculate the number of high expressing cells  per gene
`implementation`     |     enrichment implementation (data.table, simple, matrix)
`group_size`     |     number of genes to process together with data.table implementation (default = automatic)
`do_parallel`     |     run calculations in parallel with mclapply
`cores`     |     number of cores to use if do_parallel = TRUE
`verbose`     |     be verbose
`knn_params`     |     list of parameters to create spatial kNN network
`set.seed`     |     set a seed before kmeans binarization
`summarize`     |     summarize the p-values or adjusted p-values


## Value

data.table with results


