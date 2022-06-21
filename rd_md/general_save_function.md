# `general_save_function`

general_save_function


## Description

Function to automatically save plots to directory of interest


## Usage

```r
general_save_function(
  gobject,
  plot_object,
  save_dir = NULL,
  save_folder = NULL,
  save_name = NULL,
  default_save_name = "giotto_plot",
  save_format = c("png", "tiff", "pdf", "svg"),
  show_saved_plot = F,
  base_width = NULL,
  base_height = NULL,
  base_aspect_ratio = NULL,
  units = NULL,
  dpi = NULL,
  plot_count = NULL,
  ...
)
```


## Arguments

Argument      |Description
------------- |----------------
`gobject`     |     giotto object
`plot_object`     |     non-ggplot object to plot
`save_dir`     |     directory to save to
`save_folder`     |     folder in save_dir to save to
`save_name`     |     name of plot
`save_format`     |     format (e.g. png, tiff, pdf, ...)
`show_saved_plot`     |     load & display the saved plot
`base_width`     |     width
`base_height`     |     height
`base_aspect_ratio`     |     aspect ratio
`units`     |     units
`dpi`     |     Plot resolution
`plot_count`     |     count number for plot


