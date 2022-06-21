# `doLeidenSubCluster`

doLeidenSubCluster


## Description

Further subcluster cells using a NN-network and the Leiden algorithm


## Usage

```r
doLeidenSubCluster(
  gobject,
  feat_type = NULL,
  name = "sub_pleiden_clus",
  cluster_column = NULL,
  selected_clusters = NULL,
  hvf_param = list(reverse_log_scale = T, difference_in_cov = 1, expression_values =
    "normalized"),
  hvg_param = NULL,
  hvf_min_perc_cells = 5,
  hvg_min_perc_cells = NULL,
  hvf_mean_expr_det = 1,
  hvg_mean_expr_det = NULL,
  use_all_feats_as_hvf = FALSE,
  use_all_genes_as_hvg = NULL,
  min_nr_of_hvf = 5,
  min_nr_of_hvg = NULL,
  pca_param = list(expression_values = "normalized", scale_unit = T),
  nn_param = list(dimensions_to_use = 1:20),
  k_neighbors = 10,
  resolution = 0.5,
  n_iterations = 500,
  python_path = NULL,
  nn_network_to_use = "sNN",
  network_name = "sNN.pca",
  return_gobject = TRUE,
  verbose = T
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object
`feat_type`     |     feature type
`name`     |     name for new clustering result
`cluster_column`     |     cluster column to subcluster
`selected_clusters`     |     only do subclustering on these clusters
`hvf_param`     |     parameters for calculateHVf
`hvg_param`     |     deprecatd, use hvf_param
`hvf_min_perc_cells`     |     threshold for detection in min percentage of cells
`hvg_min_perc_cells`     |     deprecated, use hvf_min_perc_cells
`hvf_mean_expr_det`     |     threshold for mean expression level in cells with detection
`hvg_mean_expr_det`     |     deprecated, use hvf_mean_expr_det
`use_all_feats_as_hvf`     |     forces all features to be HVF and to be used as input for PCA
`use_all_genes_as_hvg`     |     deprecated, use use_all_feats_as_hvf
`min_nr_of_hvf`     |     minimum number of HVF, or all features will be used as input for PCA
`min_nr_of_hvg`     |     deprecated, use min_nr_of_hvf
`pca_param`     |     parameters for runPCA
`nn_param`     |     parameters for parameters for createNearestNetwork
`k_neighbors`     |     number of k for createNearestNetwork
`resolution`     |     resolution of Leiden clustering
`n_iterations`     |     number of interations to run the Leiden algorithm.
`python_path`     |     specify specific path to python if required
`nn_network_to_use`     |     type of NN network to use (kNN vs sNN)
`network_name`     |     name of NN network to use
`return_gobject`     |     boolean: return giotto object (default = TRUE)
`verbose`     |     verbose


## Details

This function performs subclustering using the Leiden algorithm on selected clusters.
 The systematic steps are:
   

*  1. subset Giotto object   

*  2. identify highly variable fetures   

*  3. run PCA   

*  4. create nearest neighbouring network   

*  5. do Leiden clustering


## Value

giotto object with new subclusters appended to cell metadata


## Seealso

[`doLeidenCluster`](#doleidencluster)


