# `run_spatial_sim_tests_one_rep`

run_spatial_sim_tests_one_rep


## Description

runs all spatial tests for 1 probability and 1 rep


## Usage

```r
run_spatial_sim_tests_one_rep(
  gobject,
  pattern_name = "pattern",
  pattern_cell_ids = NULL,
  gene_name = NULL,
  spatial_prob = 0.95,
  show_pattern = FALSE,
  spatial_network_name = "kNN_network",
  spat_methods = c("binSpect_single", "binSpect_multi", "spatialDE", "spark",
    "silhouetteRank"),
  spat_methods_params = list(NA, NA, NA, NA, NA),
  spat_methods_names = c("binSpect_single", "binSpect_multi", "spatialDE", "spark",
    "silhouetteRank"),
  save_plot = FALSE,
  save_raw = FALSE,
  save_norm = FALSE,
  save_dir = "~",
  save_name = "plot",
  run_simulations = TRUE,
  ...
)
```


