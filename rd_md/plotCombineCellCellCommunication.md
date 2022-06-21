# `plotCombineCellCellCommunication`

plotCombineCellCellCommunication


## Description

Create visualization for combined (pairwise) cell proximity gene scores


## Usage

```r
plotCombineCellCellCommunication(
  gobject,
  combCCcom,
  selected_LR = NULL,
  selected_cell_LR = NULL,
  detail_plot = T,
  simple_plot = F,
  simple_plot_facet = c("interaction", "genes"),
  facet_scales = "fixed",
  facet_ncol = length(selected_LR),
  facet_nrow = length(selected_cell_LR),
  colors = c("#9932CC", "#FF8C00"),
  show_plot = NA,
  return_plot = NA,
  save_plot = NA,
  save_param = list(),
  default_save_name = "plotCombineCellCellCommunication"
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object
`combCCcom`     |     combined communcation scores, output from combCCcom()
`selected_LR`     |     selected ligand-receptor pair
`selected_cell_LR`     |     selected cell-cell interaction pair for ligand-receptor pair
`detail_plot`     |     show detailed info in both interacting cell types
`simple_plot`     |     show a simplified plot
`simple_plot_facet`     |     facet on interactions or genes with simple plot
`facet_scales`     |     ggplot facet scales paramter
`facet_ncol`     |     ggplot facet ncol parameter
`facet_nrow`     |     ggplot facet nrow parameter
`colors`     |     vector with two colors to use
`show_plot`     |     show plots
`return_plot`     |     return plotting object
`save_plot`     |     directly save the plot [boolean]
`save_param`     |     list of saving parameters from [`all_plots_save_function`](#allplotssavefunction)
`default_save_name`     |     default save name for saving, don't change, change save_name in save_param


## Value

ggplot


