# `annotateSpatialNetwork`

annotateSpatialNetwork


## Description

Annotate spatial network with cell metadata information.


## Usage

```r
annotateSpatialNetwork(
  gobject,
  feat_type = NULL,
  spat_unit = NULL,
  spatial_network_name = "Delaunay_network",
  cluster_column,
  create_full_network = FALSE
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object
`feat_type`     |     feature type
`spat_unit`     |     spatial unit
`spatial_network_name`     |     name of spatial network to use
`cluster_column`     |     name of column to use for clusters
`create_full_network`     |     convert from reduced to full network representation


## Value

annotated network in data.table format


