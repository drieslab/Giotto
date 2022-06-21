# `exportGiottoViewer`

exportGiottoViewer


## Description

compute highly variable genes


## Usage

```r
exportGiottoViewer(
  gobject,
  feat_type = NULL,
  spat_loc_name = "raw",
  output_directory = NULL,
  spat_enr_names = NULL,
  factor_annotations = NULL,
  numeric_annotations = NULL,
  dim_reductions,
  dim_reduction_names,
  expression_values = c("scaled", "normalized", "custom"),
  dim_red_rounding = NULL,
  dim_red_rescale = c(-20, 20),
  expression_rounding = 2,
  overwrite_dir = TRUE,
  verbose = T
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object
`feat_type`     |     feature types
`spat_loc_name`     |     name of spatial locations to export
`output_directory`     |     directory where to save the files
`spat_enr_names`     |     spatial enrichment results to include for annotations
`factor_annotations`     |     giotto cell annotations to view as factor
`numeric_annotations`     |     giotto cell annotations to view as numeric
`dim_reductions`     |     high level dimension reductions to view
`dim_reduction_names`     |     specific dimension reduction names
`expression_values`     |     expression values to use in Viewer
`dim_red_rounding`     |     numerical indicating how to round the coordinates
`dim_red_rescale`     |     numericals to rescale the coordinates
`expression_rounding`     |     numerical indicating how to round the expression data
`overwrite_dir`     |     overwrite files in the directory if it already existed
`verbose`     |     be verbose


## Details

Giotto Viewer expects the results from Giotto Analyzer in a specific format,
 which is provided by this function. To include enrichment results from list(list("createSpatialEnrich")) 
 include the provided spatial enrichment name (default PAGE or rank)
 and add the gene signature names (.e.g cell types) to the numeric annotations parameter.


## Value

writes the necessary output to use in Giotto Viewer


## Examples

```r
data(mini_giotto_single_cell)
exportGiottoViewer(mini_giotto_single_cell)
```


