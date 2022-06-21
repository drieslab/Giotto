# `get_spatial_enrichment`

Get spatial enrichment


## Description

Function to get a spatial enrichment data.table


## Usage

```r
get_spatial_enrichment(
  gobject,
  spat_unit = NULL,
  feat_type = NULL,
  enrichm_name = "DWLS"
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object
`spat_unit`     |     spatial unit
`feat_type`     |     feature type
`enrichm_name`     |     name of spatial enrichment results


## Value

data.table with fractions


## Seealso

Other spatial enrichment data accessor functions:
 [`set_spatial_enrichment`](#setspatialenrichment) 
 
 Other functions to get data from giotto object:
 [`get_NearestNetwork`](#getnearestnetwork) ,
 [`get_dimReduction`](#getdimreduction) ,
 [`get_expression_values`](#getexpressionvalues) ,
 [`get_feature_info`](#getfeatureinfo) ,
 [`get_giottoImage`](#getgiottoimage) ,
 [`get_polygon_info`](#getpolygoninfo) ,
 [`get_spatialGrid`](#getspatialgrid) ,
 [`get_spatialNetwork`](#getspatialnetwork) ,
 [`get_spatial_locations`](#getspatiallocations)


