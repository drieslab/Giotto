# `plotRankSpatvsExpr`

plotRankSpatvsExpr


## Description

Plots dotplot to compare ligand-receptor rankings from spatial and expression information


## Usage

```r
plotRankSpatvsExpr(
  gobject,
  combCC,
  expr_rnk_column = "LR_expr_rnk",
  spat_rnk_column = "LR_spat_rnk",
  midpoint = 10,
  size_range = c(0.01, 1.5),
  xlims = NULL,
  ylims = NULL,
  selected_ranks = c(1, 10, 20),
  show_plot = NA,
  return_plot = NA,
  save_plot = NA,
  save_param = list(),
  default_save_name = "plotRankSpatvsExpr"
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object
`combCC`     |     combined communinication scores from [`combCCcom`](#combcccom)
`expr_rnk_column`     |     column with expression rank information to use
`spat_rnk_column`     |     column with spatial rank information to use
`midpoint`     |     midpoint of colors
`size_range`     |     size ranges of dotplot
`xlims`     |     x-limits, numerical vector of 2
`ylims`     |     y-limits, numerical vector of 2
`selected_ranks`     |     numerical vector, will be used to print out the percentage of top spatial ranks are recovered
`show_plot`     |     show plots
`return_plot`     |     return plotting object
`save_plot`     |     directly save the plot [boolean]
`save_param`     |     list of saving parameters from [`all_plots_save_function`](#allplotssavefunction)
`default_save_name`     |     default save name for saving, don't change, change save_name in save_param


## Value

ggplot


