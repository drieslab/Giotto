# `doLouvainCluster`

doLouvainCluster


## Description

cluster cells using a NN-network and the Louvain algorithm.


## Usage

```r
doLouvainCluster(
  gobject,
  version = c("community", "multinet"),
  name = "louvain_clus",
  nn_network_to_use = "sNN",
  network_name = "sNN.pca",
  python_path = NULL,
  resolution = 1,
  weight_col = NULL,
  gamma = 1,
  omega = 1,
  louv_random = F,
  return_gobject = TRUE,
  set_seed = F,
  seed_number = 1234,
  ...
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object
`version`     |     implemented version of Louvain clustering to use
`name`     |     name for cluster
`nn_network_to_use`     |     type of NN network to use (kNN vs sNN)
`network_name`     |     name of NN network to use
`python_path`     |     [community] specify specific path to python if required
`resolution`     |     [community] resolution
`weight_col`     |     weight column name
`gamma`     |     [multinet] Resolution parameter for modularity in the generalized louvain method.
`omega`     |     [multinet] Inter-layer weight parameter in the generalized louvain method
`louv_random`     |     [community] Will randomize the node evaluation order and the community evaluation order to get different partitions at each call
`return_gobject`     |     boolean: return giotto object (default = TRUE)
`set_seed`     |     set seed
`seed_number`     |     number for seed
`list()`     |     additional parameters


## Details

Louvain clustering using the community or multinet implementation of the louvain clustering algorithm.


## Value

giotto object with new clusters appended to cell metadata


## Seealso

[`doLouvainCluster_community`](#dolouvainclustercommunity) and [`doLouvainCluster_multinet`](#dolouvainclustermultinet)


