# `jackstrawPlot`

jackstrawPlot


## Description

identify significant prinicipal components (PCs)


## Usage

```r
jackstrawPlot(
  gobject,
  spat_unit = NULL,
  feat_type = NULL,
  expression_values = c("normalized", "scaled", "custom"),
  reduction = c("cells", "feats"),
  feats_to_use = NULL,
  genes_to_use = NULL,
  center = FALSE,
  scale_unit = FALSE,
  ncp = 20,
  ylim = c(0, 1),
  iter = 10,
  threshold = 0.01,
  verbose = TRUE,
  show_plot = NA,
  return_plot = NA,
  save_plot = NA,
  save_param = list(),
  default_save_name = "jackstrawPlot"
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object
`spat_unit`     |     spatial unit
`feat_type`     |     feature type
`expression_values`     |     expression values to use
`reduction`     |     cells or genes
`feats_to_use`     |     subset of features to use for PCA
`genes_to_use`     |     deprecated, use feats_to_use
`center`     |     center data before PCA
`scale_unit`     |     scale features before PCA
`ncp`     |     number of principal components to calculate
`ylim`     |     y-axis limits on jackstraw plot
`iter`     |     number of interations for jackstraw
`threshold`     |     p-value threshold to call a PC significant
`verbose`     |     show progress of jackstraw method
`show_plot`     |     show plot
`return_plot`     |     return ggplot object
`save_plot`     |     directly save the plot [boolean]
`save_param`     |     list of saving parameters from all_plots_save_function()
`default_save_name`     |     default save name for saving, don't change, change save_name in save_param


## Details

The Jackstraw method uses the [`permutationPA`](#permutationpa) function. By
 systematically permuting genes it identifies robust, and thus significant, PCs.
  list()


## Value

ggplot object for jackstraw method


