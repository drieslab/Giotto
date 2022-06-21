# `set_dimReduction`

Set dimension reduction


## Description

Function to set a dimension reduction slot


## Usage

```r
set_dimReduction(
  gobject,
  spat_unit = NULL,
  feat_type = NULL,
  reduction = c("cells", "genes"),
  reduction_method = c("pca", "umap", "tsne"),
  name = "pca",
  dimObject
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object
`spat_unit`     |     spatial unit
`feat_type`     |     feature type
`reduction`     |     reduction on cells or features
`reduction_method`     |     reduction method (e.g. pca)
`name`     |     name of reduction results
`dimObject`     |     dimension object result to set


## Value

giotto object


## Seealso

Other dimensional reduction data accessor functions:
 [`get_dimReduction`](#getdimreduction) 
 
 Other functions to set data in giotto object:
 [`set_NearestNetwork`](#setnearestnetwork) ,
 [`set_expression_values`](#setexpressionvalues) ,
 [`set_feature_info`](#setfeatureinfo) ,
 [`set_giottoImage`](#setgiottoimage) ,
 [`set_polygon_info`](#setpolygoninfo) ,
 [`set_spatialGrid`](#setspatialgrid) ,
 [`set_spatialNetwork`](#setspatialnetwork) ,
 [`set_spatial_enrichment`](#setspatialenrichment) ,
 [`set_spatial_locations`](#setspatiallocations)


