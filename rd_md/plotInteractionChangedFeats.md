# `plotInteractionChangedFeats`

Plot interaction changed features


## Description

Create barplot to visualize interaction changed features


## Usage

```r
plotInteractionChangedFeats(
  gobject,
  cpgObject,
  source_type,
  source_markers,
  ICF_feats,
  cell_color_code = NULL,
  show_plot = NA,
  return_plot = NA,
  save_plot = NA,
  save_param = list(),
  default_save_name = "plotInteractionChangedFeats"
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object
`cpgObject`     |     ICF (interaction changed feature) score object
`source_type`     |     cell type of the source cell
`source_markers`     |     markers for the source cell type
`ICF_feats`     |     named character vector of ICF features
`cell_color_code`     |     cell color code for the interacting cell types
`show_plot`     |     show plots
`return_plot`     |     return plotting object
`save_plot`     |     directly save the plot [boolean]
`save_param`     |     list of saving parameters from [`all_plots_save_function`](#allplotssavefunction)
`default_save_name`     |     default save name for saving, don't change, change save_name in save_param


## Value

plot


