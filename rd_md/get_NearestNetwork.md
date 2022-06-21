# `get_NearestNetwork`

Get nearest network


## Description

Get a NN-network from a Giotto object


## Usage

```r
get_NearestNetwork(
  gobject,
  spat_unit = NULL,
  feat_type = NULL,
  nn_network_to_use = NULL,
  network_name = NULL,
  output = c("igraph", "data.table")
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object
`spat_unit`     |     spatial unit
`feat_type`     |     feature type
`nn_network_to_use`     |     kNN or sNN
`network_name`     |     name of NN network to be used
`output`     |     return a igraph or data.table object


## Value

igraph or data.table object


## Seealso

Other expression space nearest network accessor functions:
 [`set_NearestNetwork`](#setnearestnetwork) 
 
 Other functions to get data from giotto object:
 [`get_dimReduction`](#getdimreduction) ,
 [`get_expression_values`](#getexpressionvalues) ,
 [`get_feature_info`](#getfeatureinfo) ,
 [`get_giottoImage`](#getgiottoimage) ,
 [`get_polygon_info`](#getpolygoninfo) ,
 [`get_spatialGrid`](#getspatialgrid) ,
 [`get_spatialNetwork`](#getspatialnetwork) ,
 [`get_spatial_enrichment`](#getspatialenrichment) ,
 [`get_spatial_locations`](#getspatiallocations)


