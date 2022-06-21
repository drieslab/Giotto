# `runPatternSimulation`

runPatternSimulation


## Description

Creates a known spatial pattern for selected genes one-by-one and runs the different spatial gene detection tests


## Usage

```r
runPatternSimulation(
  gobject,
  pattern_name = "pattern",
  pattern_colors = c(`in` = "green", out = "red"),
  pattern_cell_ids = NULL,
  gene_names = NULL,
  spatial_probs = c(0.5, 1),
  reps = 2,
  spatial_network_name = "kNN_network",
  spat_methods = c("binSpect_single", "binSpect_multi", "spatialDE", "spark",
    "silhouetteRank"),
  spat_methods_params = list(NA, NA, NA, NA, NA),
  spat_methods_names = c("binSpect_single", "binSpect_multi", "spatialDE", "spark",
    "silhouetteRank"),
  scalefactor = 6000,
  save_plot = T,
  save_raw = T,
  save_norm = T,
  save_dir = "~",
  max_col = 4,
  height = 7,
  width = 7,
  run_simulations = TRUE,
  ...
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object
`pattern_name`     |     name of spatial pattern
`pattern_colors`     |     2 color vector for the spatial pattern
`pattern_cell_ids`     |     cell ids that make up the spatial pattern
`gene_names`     |     selected genes
`spatial_probs`     |     probabilities to test for a high expressing gene value to be part of the spatial pattern
`reps`     |     number of random simulation repetitions
`spatial_network_name`     |     which spatial network to use for binSpectSingle
`spat_methods`     |     vector of spatial methods to test
`spat_methods_params`     |     list of parameters list for each element in the vector of spatial methods to test
`spat_methods_names`     |     name for each element in the vector of spatial elements to test
`scalefactor`     |     library size scaling factor when re-normalizing dataset
`save_plot`     |     save intermediate random simulation plots or not
`save_raw`     |     save the raw expression matrix of the simulation
`save_norm`     |     save the normalized expression matrix of the simulation
`save_dir`     |     directory to save results to
`max_col`     |     maximum number of columns for final plots
`height`     |     height of final plots
`width`     |     width of final plots
`run_simulations`     |     run simulations (default = TRUE)
`list()`     |     additional parameters for renormalization


## Value

data.table with results


