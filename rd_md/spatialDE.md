# `spatialDE`

spatialDE


## Description

Compute spatial variable genes with spatialDE method


## Usage

```r
spatialDE(
  gobject = NULL,
  feat_type = NULL,
  spat_unit = NULL,
  spat_loc_name = "raw",
  expression_values = c("raw", "normalized", "scaled", "custom"),
  size = c(4, 2, 1),
  color = c("blue", "green", "red"),
  sig_alpha = 0.5,
  unsig_alpha = 0.5,
  python_path = NULL,
  show_plot = NA,
  return_plot = NA,
  save_plot = NA,
  save_param = list(),
  default_save_name = "SpatialDE"
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     Giotto object
`feat_type`     |     feature type
`spat_unit`     |     spatial unit
`spat_loc_name`     |     name for spatial locations
`expression_values`     |     gene expression values to use
`size`     |     size of plot
`color`     |     low/medium/high color scheme for plot
`sig_alpha`     |     alpha value for significance
`unsig_alpha`     |     alpha value for unsignificance
`python_path`     |     specify specific path to python if required
`show_plot`     |     show plot
`return_plot`     |     return ggplot object
`save_plot`     |     directly save the plot [boolean]
`save_param`     |     list of saving parameters, see [`showSaveParameters`](#showsaveparameters)
`default_save_name`     |     default save name for saving, don't change, change save_name in save_param


## Details

This function is a wrapper for the SpatialDE method implemented in the ...


## Value

a list of data.frames with results and plot (optional)


