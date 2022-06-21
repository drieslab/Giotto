# `set_NearestNetwork`

Set nearest network


## Description

Set a NN-network for a Giotto object


## Usage

```r
set_NearestNetwork(
  gobject,
  spat_unit = NULL,
  feat_type = NULL,
  nn_network_to_use = "sNN",
  network_name = "sNN.pca",
  nn_network
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
`nn_network`     |     nearest network


## Value

giotto object


## Seealso

Other expression space nearest network accessor functions:
 [`get_NearestNetwork`](#getnearestnetwork) 
 
 Other functions to set data in giotto object:
 [`set_dimReduction`](#setdimreduction) ,
 [`set_expression_values`](#setexpressionvalues) ,
 [`set_feature_info`](#setfeatureinfo) ,
 [`set_giottoImage`](#setgiottoimage) ,
 [`set_polygon_info`](#setpolygoninfo) ,
 [`set_spatialGrid`](#setspatialgrid) ,
 [`set_spatialNetwork`](#setspatialnetwork) ,
 [`set_spatial_enrichment`](#setspatialenrichment) ,
 [`set_spatial_locations`](#setspatiallocations)


