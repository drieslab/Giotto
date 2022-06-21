# `clusterCells`

clusterCells


## Description

cluster cells using a variety of different methods


## Usage

```r
clusterCells(
  gobject,
  cluster_method = c("leiden", "louvain_community", "louvain_multinet", "randomwalk",
    "sNNclust", "kmeans", "hierarchical"),
  name = "cluster_name",
  nn_network_to_use = "sNN",
  network_name = "sNN.pca",
  pyth_leid_resolution = 1,
  pyth_leid_weight_col = "weight",
  pyth_leid_part_type = c("RBConfigurationVertexPartition",
    "ModularityVertexPartition"),
  pyth_leid_init_memb = NULL,
  pyth_leid_iterations = 1000,
  pyth_louv_resolution = 1,
  pyth_louv_weight_col = NULL,
  python_louv_random = F,
  python_path = NULL,
  louvain_gamma = 1,
  louvain_omega = 1,
  walk_steps = 4,
  walk_clusters = 10,
  walk_weights = NA,
  sNNclust_k = 20,
  sNNclust_eps = 4,
  sNNclust_minPts = 16,
  borderPoints = TRUE,
  expression_values = c("normalized", "scaled", "custom"),
  genes_to_use = NULL,
  dim_reduction_to_use = c("cells", "pca", "umap", "tsne"),
  dim_reduction_name = "pca",
  dimensions_to_use = 1:10,
  distance_method = c("original", "pearson", "spearman", "euclidean", "maximum",
    "manhattan", "canberra", "binary", "minkowski"),
  km_centers = 10,
  km_iter_max = 100,
  km_nstart = 1000,
  km_algorithm = "Hartigan-Wong",
  hc_agglomeration_method = c("ward.D2", "ward.D", "single", "complete", "average",
    "mcquitty", "median", "centroid"),
  hc_k = 10,
  hc_h = NULL,
  return_gobject = TRUE,
  set_seed = T,
  seed_number = 1234
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object
`cluster_method`     |     community cluster method to use
`name`     |     name for new clustering result
`nn_network_to_use`     |     type of NN network to use (kNN vs sNN)
`network_name`     |     name of NN network to use
`pyth_leid_resolution`     |     resolution for leiden
`pyth_leid_weight_col`     |     column to use for weights
`pyth_leid_part_type`     |     partition type to use
`pyth_leid_init_memb`     |     initial membership
`pyth_leid_iterations`     |     number of iterations
`pyth_louv_resolution`     |     resolution for louvain
`pyth_louv_weight_col`     |     python louvain param: weight column
`python_louv_random`     |     python louvain param: random
`python_path`     |     specify specific path to python if required
`louvain_gamma`     |     louvain param: gamma or resolution
`louvain_omega`     |     louvain param: omega
`walk_steps`     |     randomwalk: number of steps
`walk_clusters`     |     randomwalk: number of clusters
`walk_weights`     |     randomwalk: weight column
`sNNclust_k`     |     SNNclust: k neighbors to use
`sNNclust_eps`     |     SNNclust: epsilon
`sNNclust_minPts`     |     SNNclust: min points
`borderPoints`     |     SNNclust: border points
`expression_values`     |     expression values to use
`genes_to_use`     |     = NULL,
`dim_reduction_to_use`     |     dimension reduction to use
`dim_reduction_name`     |     name of reduction 'pca',
`dimensions_to_use`     |     dimensions to use
`distance_method`     |     distance method
`km_centers`     |     kmeans centers
`km_iter_max`     |     kmeans iterations
`km_nstart`     |     kmeans random starting points
`km_algorithm`     |     kmeans algorithm
`hc_agglomeration_method`     |     hierarchical clustering method
`hc_k`     |     hierachical number of clusters
`hc_h`     |     hierarchical tree cutoff
`return_gobject`     |     boolean: return giotto object (default = TRUE)
`set_seed`     |     set seed
`seed_number`     |     number for seed


## Details

Wrapper for the different clustering methods.


## Value

giotto object with new clusters appended to cell metadata


## Seealso

[`doLeidenCluster`](#doleidencluster) , [`doLouvainCluster_community`](#dolouvainclustercommunity) , [`doLouvainCluster_multinet`](#dolouvainclustermultinet) ,
 [`doLouvainCluster`](#dolouvaincluster) , [`doRandomWalkCluster`](#dorandomwalkcluster) , [`doSNNCluster`](#dosnncluster) ,
 [`doKmeans`](#dokmeans) , [`doHclust`](#dohclust)


