# `doLouvainCluster_community`

doLouvainCluster_community


## Description

cluster cells using a NN-network and the Louvain algorithm from the community module in Python


## Usage

```r
doLouvainCluster_community(
  gobject,
  name = "louvain_clus",
  nn_network_to_use = "sNN",
  network_name = "sNN.pca",
  python_path = NULL,
  resolution = 1,
  weight_col = NULL,
  louv_random = F,
  return_gobject = TRUE,
  set_seed = F,
  seed_number = 1234
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object
`name`     |     name for cluster
`nn_network_to_use`     |     type of NN network to use (kNN vs sNN)
`network_name`     |     name of NN network to use
`python_path`     |     specify specific path to python if required
`resolution`     |     resolution
`weight_col`     |     weight column to use for edges
`louv_random`     |     Will randomize the node evaluation order and the community evaluation order to get different partitions at each call
`return_gobject`     |     boolean: return giotto object (default = TRUE)
`set_seed`     |     set seed
`seed_number`     |     number for seed


## Details

This function is a wrapper for the Louvain algorithm implemented in Python,
 which can detect communities in graphs of nodes (cells).
 See the [https://python-louvain.readthedocs.io/en/latest/index.html](https://python-louvain.readthedocs.io/en/latest/index.html) readthedocs 
 page for more information.
 
 Set weight_col = NULL to give equal weight (=1) to each edge.


## Value

giotto object with new clusters appended to cell metadata


