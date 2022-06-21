# `cellProximityBarplot`

cellProximityBarplot


## Description

Create barplot from cell-cell proximity scores


## Usage

```r
cellProximityBarplot(
  gobject,
  CPscore,
  min_orig_ints = 5,
  min_sim_ints = 5,
  p_val = 0.05,
  show_plot = NA,
  return_plot = NA,
  save_plot = NA,
  save_param = list(),
  default_save_name = "cellProximityBarplot"
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object
`CPscore`     |     CPscore, output from cellProximityEnrichment()
`min_orig_ints`     |     filter on minimum original cell-cell interactions
`min_sim_ints`     |     filter on minimum simulated cell-cell interactions
`p_val`     |     p-value
`show_plot`     |     show plot
`return_plot`     |     return ggplot object
`save_plot`     |     directly save the plot [boolean]
`save_param`     |     list of saving parameters from [`all_plots_save_function`](#allplotssavefunction)
`default_save_name`     |     default save name for saving, don't change, change save_name in save_param


## Details

This function creates a barplot that shows the  spatial proximity
 enrichment or depletion of cell type pairs.


## Value

ggplot barplot


