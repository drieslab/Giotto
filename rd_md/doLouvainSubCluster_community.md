# `doLouvainSubCluster_community`

doLouvainSubCluster_community


## Description

subcluster cells using a NN-network and the Louvain community detection algorithm


## Usage

```r
doLouvainSubCluster_community(
  gobject,
  name = "sub_louvain_comm_clus",
  cluster_column = NULL,
  selected_clusters = NULL,
  hvg_param = list(reverse_log_scale = T, difference_in_cov = 1, expression_values =
    "normalized"),
  hvg_min_perc_cells = 5,
  hvg_mean_expr_det = 1,
  use_all_genes_as_hvg = FALSE,
  min_nr_of_hvg = 5,
  pca_param = list(expression_values = "normalized", scale_unit = T),
  nn_param = list(dimensions_to_use = 1:20),
  k_neighbors = 10,
  resolution = 0.5,
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
`name`     |     name for new clustering result
`cluster_column`     |     cluster column to subcluster
`selected_clusters`     |     only do subclustering on these clusters
`hvg_param`     |     parameters for calculateHVG
`hvg_min_perc_cells`     |     threshold for detection in min percentage of cells
`hvg_mean_expr_det`     |     threshold for mean expression level in cells with detection
`use_all_genes_as_hvg`     |     forces all genes to be HVG and to be used as input for PCA
`min_nr_of_hvg`     |     minimum number of HVG, or all genes will be used as input for PCA
`pca_param`     |     parameters for runPCA
`nn_param`     |     parameters for parameters for createNearestNetwork
`k_neighbors`     |     number of k for createNearestNetwork
`resolution`     |     resolution
`python_path`     |     specify specific path to python if required
`nn_network_to_use`     |     type of NN network to use (kNN vs sNN)
`network_name`     |     name of NN network to use
`return_gobject`     |     boolean: return giotto object (default = TRUE)
`verbose`     |     verbose


## Details

This function performs subclustering using the Louvain community algorithm on selected clusters.
 The systematic steps are:
   

*  1. subset Giotto object   

*  2. identify highly variable genes   

*  3. run PCA   

*  4. create nearest neighbouring network   

*  5. do Louvain community clustering


## Value

giotto object with new subclusters appended to cell metadata


## Seealso

[`doLouvainCluster_community`](#dolouvainclustercommunity)


