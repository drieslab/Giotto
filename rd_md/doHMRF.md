# `doHMRF`

doHMRF


## Description

Run HMRF


## Usage

```r
doHMRF(
  gobject,
  spat_unit = NULL,
  feat_type = NULL,
  expression_values = c("normalized", "scaled", "custom"),
  spatial_network_name = "Delaunay_network",
  spat_loc_name = "raw",
  spatial_genes = NULL,
  spatial_dimensions = c("sdimx", "sdimy", "sdimz"),
  dim_reduction_to_use = NULL,
  dim_reduction_name = "pca",
  dimensions_to_use = 1:10,
  seed = 100,
  name = "test",
  k = 10,
  betas = c(0, 2, 50),
  tolerance = 1e-10,
  zscore = c("none", "rowcol", "colrow"),
  numinit = 100,
  python_path = NULL,
  output_folder = NULL,
  overwrite_output = TRUE
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object
`spat_unit`     |     spatial unit
`feat_type`     |     feature type
`expression_values`     |     expression values to use
`spatial_network_name`     |     name of spatial network to use for HMRF
`spat_loc_name`     |     name of spatial locations
`spatial_genes`     |     spatial genes to use for HMRF
`spatial_dimensions`     |     select spatial dimensions to use, default is all possible dimensions
`dim_reduction_to_use`     |     use another dimension reduction set as input
`dim_reduction_name`     |     name of dimension reduction set to use
`dimensions_to_use`     |     number of dimensions to use as input
`seed`     |     seed to fix random number generator (for creating initialization of HMRF) (-1 if no fixing)
`name`     |     name of HMRF run
`k`     |     number of HMRF domains
`betas`     |     betas to test for. three numbers: start_beta, beta_increment, num_betas e.g. c(0, 2.0, 50)
`tolerance`     |     tolerance
`zscore`     |     zscore
`numinit`     |     number of initializations
`python_path`     |     python path to use
`output_folder`     |     output folder to save results
`overwrite_output`     |     overwrite output folder


## Details

Description of HMRF parameters ...


## Value

Creates a directory with results that can be viewed with viewHMRFresults


