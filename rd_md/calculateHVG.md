# `calculateHVG`

calculateHVG


## Description

compute highly variable genes


## Usage

```r
calculateHVG(
  gobject,
  expression_values = c("normalized", "scaled", "custom"),
  method = c("cov_groups", "cov_loess", "var_p_resid"),
  reverse_log_scale = FALSE,
  logbase = 2,
  expression_threshold = 0,
  nr_expression_groups = 20,
  zscore_threshold = 1.5,
  HVGname = "hvf",
  difference_in_cov = 0.1,
  var_threshold = 1.5,
  var_number = NULL,
  show_plot = NA,
  return_plot = NA,
  save_plot = NA,
  save_param = list(),
  default_save_name = "HVGplot",
  return_gobject = TRUE
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object
`expression_values`     |     expression values to use
`method`     |     method to calculate highly variable genes
`reverse_log_scale`     |     reverse log-scale of expression values (default = FALSE)
`logbase`     |     if reverse_log_scale is TRUE, which log base was used?
`expression_threshold`     |     expression threshold to consider a gene detected
`nr_expression_groups`     |     [cov_groups] number of expression groups for cov_groups
`zscore_threshold`     |     [cov_groups] zscore to select hvg for cov_groups
`HVGname`     |     name for highly variable genes in cell metadata
`difference_in_cov`     |     [cov_loess] minimum difference in coefficient of variance required
`var_threshold`     |     [var_p_resid] variance threshold for features for var_p_resid method
`var_number`     |     [var_p_resid] number of top variance features for var_p_resid method
`show_plot`     |     show plot
`return_plot`     |     return ggplot object
`save_plot`     |     directly save the plot [boolean]
`save_param`     |     list of saving parameters from [`all_plots_save_function`](#allplotssavefunction)
`default_save_name`     |     default save name for saving, don't change, change save_name in save_param
`return_gobject`     |     boolean: return giotto object (default = TRUE)


## Details

Currently we provide 2 ways to calculate highly variable genes:
 
 list("1. high coeff of variance (COV) within groups: ")  list() 
 First genes are binned ( list("nr_expression_groups") ) into average expression groups and
 the COV for each gene is converted into a z-score within each bin. Genes with a z-score
 higher than the threshold ( list("zscore_threshold") ) are considered highly variable. list() 
 
 list("2. high COV based on loess regression prediction: ")  list() 
 A predicted COV is calculated for each gene using loess regression (COV~log(mean expression))
 Genes that show a higher than predicted COV ( list("difference_in_cov") ) are considered highly variable. list()


## Value

giotto object highly variable genes appended to gene metadata (fDataDT)


