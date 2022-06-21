# `doLeidenCluster`

doLeidenCluster


## Description

cluster cells using a NN-network and the Leiden community detection algorithm


## Usage

```r
doLeidenCluster(
  gobject,
  spat_unit = NULL,
  feat_type = NULL,
  name = "leiden_clus",
  nn_network_to_use = "sNN",
  network_name = "sNN.pca",
  python_path = NULL,
  resolution = 1,
  weight_col = "weight",
  partition_type = c("RBConfigurationVertexPartition", "ModularityVertexPartition"),
  init_membership = NULL,
  n_iterations = 1000,
  return_gobject = TRUE,
  set_seed = T,
  seed_number = 1234
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object
`spat_unit`     |     spatial unit
`feat_type`     |     feature type
`name`     |     name for cluster
`nn_network_to_use`     |     type of NN network to use (kNN vs sNN)
`network_name`     |     name of NN network to use
`python_path`     |     specify specific path to python if required
`resolution`     |     resolution
`weight_col`     |     weight column to use for edges
`partition_type`     |     The type of partition to use for optimisation.
`init_membership`     |     initial membership of cells for the partition
`n_iterations`     |     number of interations to run the Leiden algorithm. If the number of iterations is negative, the Leiden algorithm is run until an iteration in which there was no improvement.
`return_gobject`     |     boolean: return giotto object (default = TRUE)
`set_seed`     |     set seed
`seed_number`     |     number for seed


## Details

This function is a wrapper for the Leiden algorithm implemented in python,
 which can detect communities in graphs of millions of nodes (cells),
 as long as they can fit in memory. See the [https://github.com/vtraag/leidenalg](https://github.com/vtraag/leidenalg) leidenalg 
 github page or the [https://leidenalg.readthedocs.io/en/stable/index.html](https://leidenalg.readthedocs.io/en/stable/index.html) readthedocs 
 page for more information.
 
 Partition types available and information:
   

*  list("RBConfigurationVertexPartition: ") list("Implements Reichardt and Bornholdt’s Potts model\n", "   with a configuration null model. This quality function is well-defined only for positive edge weights.\n", "    This quality function uses a linear resolution parameter.")   

*  list("ModularityVertexPartition: ") list("Implements modularity.\n", "   This quality function is well-defined only for positive edge weights. It does ", list("not"), " use the resolution parameter")  
 
 Set weight_col = NULL to give equal weight (=1) to each edge.


## Value

giotto object with new clusters appended to cell metadata


