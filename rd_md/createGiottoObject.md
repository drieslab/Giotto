# `createGiottoObject`

Create a giotto object


## Description

Function to create a giotto object


## Usage

```r
createGiottoObject(
  expression,
  raw_exprs = NULL,
  expression_feat = "rna",
  spatial_locs = NULL,
  spatial_info = NULL,
  cell_metadata = NULL,
  feat_metadata = NULL,
  feat_info = NULL,
  spatial_network = NULL,
  spatial_grid = NULL,
  spatial_grid_name = NULL,
  spatial_enrichment = NULL,
  dimension_reduction = NULL,
  nn_network = NULL,
  images = NULL,
  offset_file = NULL,
  instructions = NULL,
  cores = NA,
  verbose = TRUE
)
```


## Arguments

Argument      |Description
------------- |----------------
`expression`     |     expression information
`raw_exprs`     |     deprecated, use expression
`expression_feat`     |     available features (e.g. rna, protein, ...)
`spatial_locs`     |     data.table or data.frame with coordinates for cell centroids
`spatial_info`     |     list of giotto polygon objects with spatial information, see [`createGiottoPolygonsFromMask`](#creategiottopolygonsfrommask) and [`createGiottoPolygonsFromDfr`](#creategiottopolygonsfromdfr)
`cell_metadata`     |     cell annotation metadata
`feat_metadata`     |     feature annotation metadata for each unique feature
`feat_info`     |     list of giotto point objects with feature info, see [`createGiottoPoints`](#creategiottopoints)
`spatial_network`     |     list of spatial network(s)
`spatial_grid`     |     list of spatial grid(s)
`spatial_grid_name`     |     list of spatial grid name(s)
`spatial_enrichment`     |     list of spatial enrichment score(s) for each spatial region
`dimension_reduction`     |     list of dimension reduction(s)
`nn_network`     |     list of nearest neighbor network(s)
`images`     |     list of images
`offset_file`     |     file used to stitch fields together (optional)
`instructions`     |     list of instructions or output result from [`createGiottoInstructions`](#creategiottoinstructions)
`cores`     |     how many cores or threads to use to read data if paths are provided
`verbose`     |     be verbose when building Giotto object


## Details

See [http://giottosuite.com/articles/getting_started_gobject.html](http://giottosuite.com/articles/getting_started_gobject.html) for more details
 
 [ Requirements ] To create a giotto object you need to provide at least a matrix with genes as
 row names and cells as column names. This matrix can be provided as a base matrix, sparse Matrix, data.frame,
 data.table or as a path to any of those.
 To include spatial information about cells (or regions) you need to provide a matrix, data.table or data.frame (or path to them)
 with coordinates for all spatial dimensions. This can be 2D (x and y) or 3D (x, y, x).
 The row order for the cell coordinates should be the same as the column order for the provided expression data.
 
 [ Instructions ] Additionally an instruction file, generated manually or with [`createGiottoInstructions`](#creategiottoinstructions) 
 can be provided to instructions, if not a default instruction file will be created
 for the Giotto object.
 
 [ Multiple fields ] In case a dataset consists of multiple fields, like seqFISH+ for example,
 an offset file can be provided to stitch the different fields together. [`stitchFieldCoordinates`](#stitchfieldcoordinates) 
 can be used to generate such an offset file.
 
 [ Processed data ] Processed count data, such as normalized data, can be provided using
 one of the different expression slots (norm_expr, norm_scaled_expr, custom_expr).
 
 [ Metadata ] Cell and gene metadata can be provided using the cell and gene metadata slots.
 This data can also be added afterwards using the [`addFeatMetadata`](#addfeatmetadata) or [`addCellMetadata`](#addcellmetadata) functions.
 
 [ Other information ] Additional information can be provided through the appropriate slots:
   

*  spatial networks   

*  spatial girds   

*  spatial enrichments   

*  dimensions reduction   

*  nearest neighbours networks   

*  images


## Value

giotto object


