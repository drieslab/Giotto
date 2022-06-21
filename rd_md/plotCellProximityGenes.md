# `plotCellProximityGenes`

plotCellProximityGenes


## Description

Create visualization for cell proximity gene scores


## Usage

```r
plotCellProximityGenes(
  gobject,
  cpgObject,
  method = c("volcano", "cell_barplot", "cell-cell", "cell_sankey", "heatmap",
    "dotplot"),
  min_cells = 4,
  min_cells_expr = 1,
  min_int_cells = 4,
  min_int_cells_expr = 1,
  min_fdr = 0.1,
  min_spat_diff = 0.2,
  min_log2_fc = 0.2,
  min_zscore = 2,
  zscores_column = c("cell_type", "feats"),
  direction = c("both", "up", "down"),
  cell_color_code = NULL,
  show_plot = NA,
  return_plot = NA,
  save_plot = NA,
  save_param = list(),
  default_save_name = "plotCellProximityGenes"
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object
`cpgObject`     |     ICG (interaction changed gene) score object
`method`     |     plotting method to use
`min_cells`     |     minimum number of source cell type
`min_cells_expr`     |     minimum expression level for source cell type
`min_int_cells`     |     minimum number of interacting neighbor cell type
`min_int_cells_expr`     |     minimum expression level for interacting neighbor cell type
`min_fdr`     |     minimum adjusted p-value
`min_spat_diff`     |     minimum absolute spatial expression difference
`min_log2_fc`     |     minimum log2 fold-change
`min_zscore`     |     minimum z-score change
`zscores_column`     |     calculate z-scores over cell types or genes
`direction`     |     differential expression directions to keep
`cell_color_code`     |     vector of colors with cell types as names
`show_plot`     |     show plots
`return_plot`     |     return plotting object
`save_plot`     |     directly save the plot [boolean]
`save_param`     |     list of saving parameters from [`all_plots_save_function`](#allplotssavefunction)
`default_save_name`     |     default save name for saving, don't change, change save_name in save_param


## Value

plot


