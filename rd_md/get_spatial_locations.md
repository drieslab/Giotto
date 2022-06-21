# `get_spatial_locations`

Get spatial locations


## Description

Function to get a spatial location data.table


## Usage

```r
get_spatial_locations(gobject, spat_unit = NULL, spat_loc_name = NULL)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object
`spat_unit`     |     spatial unit
`spat_loc_name`     |     name of spatial locations (defaults to first name in spatial_locs slot)


## Value

data.table with coordinates


## Seealso

Other spatial location data accessor functions:
 [`set_spatial_locations`](#setspatiallocations) 
 
 Other functions to get data from giotto object:
 [`get_NearestNetwork`](#getnearestnetwork) ,
 [`get_dimReduction`](#getdimreduction) ,
 [`get_expression_values`](#getexpressionvalues) ,
 [`get_feature_info`](#getfeatureinfo) ,
 [`get_giottoImage`](#getgiottoimage) ,
 [`get_polygon_info`](#getpolygoninfo) ,
 [`get_spatialGrid`](#getspatialgrid) ,
 [`get_spatialNetwork`](#getspatialnetwork) ,
 [`get_spatial_enrichment`](#getspatialenrichment)


