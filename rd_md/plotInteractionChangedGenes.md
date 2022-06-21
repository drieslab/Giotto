# `plotInteractionChangedGenes`

Plot interaction changed genes


## Description

Create barplot to visualize interaction changed genes


## Usage

```r
plotInteractionChangedGenes(
  gobject,
  cpgObject,
  source_type,
  source_markers,
  ICG_genes,
  cell_color_code = NULL,
  show_plot = NA,
  return_plot = NA,
  save_plot = NA,
  save_param = list(),
  default_save_name = "plotInteractionChangedGenes"
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object
`cpgObject`     |     ICG (interaction changed gene) score object
`source_type`     |     cell type of the source cell
`source_markers`     |     markers for the source cell type
`ICG_genes`     |     named character vector of ICG genes
`cell_color_code`     |     cell color code for the interacting cell types
`show_plot`     |     show plots
`return_plot`     |     return plotting object
`save_plot`     |     directly save the plot [boolean]
`save_param`     |     list of saving parameters from [`all_plots_save_function`](#allplotssavefunction)
`default_save_name`     |     default save name for saving, don't change, change save_name in save_param


## Value

plot


