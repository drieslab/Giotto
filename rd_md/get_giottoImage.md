# `get_giottoImage`

Get giotto image object


## Description

Get giotto image object from gobject


## Usage

```r
get_giottoImage(
  gobject = NULL,
  image_type = c("image", "largeImage"),
  name = NULL
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object
`image_type`     |     type of giotto image object
`name`     |     name of a giotto image object [`showGiottoImageNames`](#showgiottoimagenames)


## Value

a giotto image object


## Seealso

Other image data accessor functions:
 [`set_giottoImage`](#setgiottoimage) 
 
 Other functions to get data from giotto object:
 [`get_NearestNetwork`](#getnearestnetwork) ,
 [`get_dimReduction`](#getdimreduction) ,
 [`get_expression_values`](#getexpressionvalues) ,
 [`get_feature_info`](#getfeatureinfo) ,
 [`get_polygon_info`](#getpolygoninfo) ,
 [`get_spatialGrid`](#getspatialgrid) ,
 [`get_spatialNetwork`](#getspatialnetwork) ,
 [`get_spatial_enrichment`](#getspatialenrichment) ,
 [`get_spatial_locations`](#getspatiallocations)


