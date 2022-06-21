# `signPCA`

signPCA


## Description

identify significant prinicipal components (PCs)


## Usage

```r
signPCA(
  gobject,
  feat_type = NULL,
  spat_unit = NULL,
  name = NULL,
  method = c("screeplot", "jackstraw"),
  expression_values = c("normalized", "scaled", "custom"),
  reduction = c("cells", "feats"),
  pca_method = c("irlba", "factominer"),
  rev = FALSE,
  feats_to_use = NULL,
  genes_to_use = NULL,
  center = T,
  scale_unit = T,
  ncp = 50,
  scree_ylim = c(0, 10),
  jack_iter = 10,
  jack_threshold = 0.01,
  jack_ylim = c(0, 1),
  verbose = TRUE,
  show_plot = NA,
  return_plot = NA,
  save_plot = NA,
  save_param = list(),
  default_save_name = "signPCA"
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object
`feat_type`     |     feature type
`spat_unit`     |     spatial unit
`name`     |     name of PCA object if available
`method`     |     method to use to identify significant PCs
`expression_values`     |     expression values to use
`reduction`     |     cells or genes
`pca_method`     |     which implementation to use
`rev`     |     do a reverse PCA
`feats_to_use`     |     subset of features to use for PCA
`genes_to_use`     |     deprecated, use feats_to_use
`center`     |     center data before PCA
`scale_unit`     |     scale features before PCA
`ncp`     |     number of principal components to calculate
`scree_ylim`     |     y-axis limits on scree plot
`jack_iter`     |     number of interations for jackstraw
`jack_threshold`     |     p-value threshold to call a PC significant
`jack_ylim`     |     y-axis limits on jackstraw plot
`verbose`     |     verbosity
`show_plot`     |     show plot
`return_plot`     |     return ggplot object
`save_plot`     |     directly save the plot [boolean]
`save_param`     |     list of saving parameters from all_plots_save_function()
`default_save_name`     |     default save name for saving, don't change, change save_name in save_param


## Details

Two different methods can be used to assess the number of relevant or significant
 prinicipal components (PC's). list() 
 1. Screeplot works by plotting the explained variance of each
 individual PC in a barplot allowing you to identify which PC provides a significant
 contribution  (a.k.a. 'elbow method'). list() 
 2. The Jackstraw method uses the [`permutationPA`](#permutationpa) function. By
 systematically permuting genes it identifies robust, and thus significant, PCs.
  list()


## Value

ggplot object for scree method and maxtrix of p-values for jackstraw


