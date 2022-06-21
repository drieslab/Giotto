# `get_dimReduction`

Get dimension reduction


## Description

Function to get a dimension reduction object


## Usage

```r
get_dimReduction(
  gobject,
  spat_unit = NULL,
  feat_type = NULL,
  reduction = c("cells", "feats"),
  reduction_method = c("pca", "umap", "tsne"),
  name = "pca",
  return_dimObj = FALSE
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
`return_dimObj`     |     return full dimension object result


## Value

dim reduction coordinates (default) or dim reduction object


## Seealso

Other dimensional reduction data accessor functions:
 [`set_dimReduction`](#setdimreduction) 
 
 Other functions to get data from giotto object:
 [`get_NearestNetwork`](#getnearestnetwork) ,
 [`get_expression_values`](#getexpressionvalues) ,
 [`get_feature_info`](#getfeatureinfo) ,
 [`get_giottoImage`](#getgiottoimage) ,
 [`get_polygon_info`](#getpolygoninfo) ,
 [`get_spatialGrid`](#getspatialgrid) ,
 [`get_spatialNetwork`](#getspatialnetwork) ,
 [`get_spatial_enrichment`](#getspatialenrichment) ,
 [`get_spatial_locations`](#getspatiallocations)


