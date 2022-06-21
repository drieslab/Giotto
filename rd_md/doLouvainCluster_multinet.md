# `doLouvainCluster_multinet`

doLouvainCluster_multinet


## Description

cluster cells using a NN-network and the Louvain algorithm from the multinet package in R.


## Usage

```r
doLouvainCluster_multinet(
  gobject,
  name = "louvain_clus",
  nn_network_to_use = "sNN",
  network_name = "sNN.pca",
  gamma = 1,
  omega = 1,
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
`gamma`     |     Resolution parameter for modularity in the generalized louvain method.
`omega`     |     Inter-layer weight parameter in the generalized louvain method.
`return_gobject`     |     boolean: return giotto object (default = TRUE)
`set_seed`     |     set seed
`seed_number`     |     number for seed


## Details

See [`glouvain_ml`](#glouvainml) from the multinet package in R for
 more information.


## Value

giotto object with new clusters appended to cell metadata


