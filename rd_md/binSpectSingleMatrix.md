# `binSpectSingleMatrix`

binSpectSingleMatrix


## Description

binSpect for a single spatial network and a provided expression matrix


## Usage

```r
binSpectSingleMatrix(
  expression_matrix,
  spatial_network = NULL,
  bin_matrix = NULL,
  bin_method = c("kmeans", "rank"),
  subset_feats = NULL,
  kmeans_algo = c("kmeans", "kmeans_arma", "kmeans_arma_subset"),
  nstart = 3,
  iter_max = 10,
  extreme_nr = 50,
  sample_nr = 50,
  percentage_rank = 30,
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
  verbose = TRUE,
  set.seed = NULL
)
```


## Arguments

Argument      |Description
------------- |----------------
`expression_matrix`     |     expression matrix
`spatial_network`     |     spatial network in data.table format
`bin_matrix`     |     a binarized matrix, when provided it will skip the binarization process
`bin_method`     |     method to binarize gene expression
`subset_feats`     |     only select a subset of features to test
`kmeans_algo`     |     kmeans algorithm to use (kmeans, kmeans_arma, kmeans_arma_subset)
`nstart`     |     kmeans: nstart parameter
`iter_max`     |     kmeans: iter.max parameter
`extreme_nr`     |     number of top and bottom cells (see details)
`sample_nr`     |     total number of cells to sample (see details)
`percentage_rank`     |     percentage of top cells for binarization
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
`set.seed`     |     set a seed before kmeans binarization


## Details

We provide two ways to identify spatial genes based on gene expression binarization.
 Both methods are identicial except for how binarization is performed.
   

*  list("1. binarize: ") list("Each gene is binarized (0 or 1) in each cell with ", list("kmeans"), " (k = 2) or based on ", list("rank"), " percentile")   

*  list("2. network: ") list("Alll cells are connected through a spatial network based on the physical coordinates")   

*  list("3. contingency table: ") list("A contingency table is calculated based on all edges of neighboring cells and the binarized expression (0-0, 0-1, 1-0 or 1-1)")   

*  list("4. For each gene an odds-ratio (OR) and fisher.test (optional) is calculated")  
 Three different kmeans algorithmes have been implemented:
   

*  list("1. kmeans: ") list("default, see ", list(list("kmeans")), " ")   

*  list("2. kmeans_arma: ") list("from ClusterR, see ", list(list("KMeans_arma")), " ")   

*  list("3. kmeans_arma_subst: ") list("from ClusterR, see ", list(list("KMeans_arma")), ",\n", "   but random subsetting the vector for each gene to increase speed. Change extreme_nr and sample_nr for control.  ")  
 Other statistics are provided (optional):
   

*  Number of cells with high expression (binary = 1)   

*  Average expression of each gene within high expressing cells    

*  Number of hub cells, these are high expressing cells that have a user defined number of high expressing neighbors  
 By selecting a subset of likely spatial genes (e.g. soft thresholding highly variable genes) can accelerate the speed.
 The simple implementation is usually faster, but lacks the possibility to run in parallel and to calculate hub cells.
 The data.table implementation might be more appropriate for large datasets by setting the group_size (number of genes) parameter to divide the workload.


## Value

data.table with results (see details)


