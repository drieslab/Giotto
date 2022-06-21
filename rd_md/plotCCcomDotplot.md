# `plotCCcomDotplot`

plotCCcomDotplot


## Description

Plots dotplot for ligand-receptor communication scores in cell-cell interactions


## Usage

```r
plotCCcomDotplot(
  gobject,
  comScores,
  selected_LR = NULL,
  selected_cell_LR = NULL,
  show_LR_names = TRUE,
  show_cell_LR_names = TRUE,
  cluster_on = c("PI", "LR_expr", "log2fc"),
  cor_method = c("pearson", "kendall", "spearman"),
  aggl_method = c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty",
    "median", "centroid"),
  show_plot = NA,
  return_plot = NA,
  save_plot = NA,
  save_param = list(),
  default_save_name = "plotCCcomDotplot"
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object
`comScores`     |     communinication scores from [`exprCellCellcom`](#exprcellcellcom) or [`spatCellCellcom`](#spatcellcellcom)
`selected_LR`     |     selected ligand-receptor combinations
`selected_cell_LR`     |     selected cell-cell combinations for ligand-receptor combinations
`show_LR_names`     |     show ligand-receptor names
`show_cell_LR_names`     |     show cell-cell names
`cluster_on`     |     values to use for clustering of cell-cell and ligand-receptor pairs
`cor_method`     |     correlation method used for clustering
`aggl_method`     |     agglomeration method used by hclust
`show_plot`     |     show plots
`return_plot`     |     return plotting object
`save_plot`     |     directly save the plot [boolean]
`save_param`     |     list of saving parameters from [`all_plots_save_function`](#allplotssavefunction)
`default_save_name`     |     default save name for saving, don't change, change save_name in save_param


## Value

ggplot


